#!/usr/bin/env python3

import json
import sys
import logging
import copy
import argparse
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd # type: ignore
from numpy.typing import ArrayLike
from dataclasses import dataclass, field
from collections import defaultdict
from importlib.metadata import version
from math import ceil, log
from statistics import fmean
from matplotlib.collections import PatchCollection
from typing import Type, Dict, List, Tuple, Optional, Union, IO
from smudgeplot.exceptions import *
from smudgeplot.config import PlotConfig, AnalysisConfig

logger = logging.getLogger(__name__)

# ==============================================================================
# CLASS: Coverages
# ==============================================================================

@dataclass
class Coverages:
    """
    Container for k-mer coverage data with aggregation and filtering methods.
    
    This class manages coverage table data and performs local aggregation to
    identify distinct k-mer coverage peaks (smudges) and calculate error rates.
    
    Attributes:
        cov_tab : pandas.DataFrame
            The coverage table data.
        cov2peak : dict
            Mapping from coverage coordinates (covA, covB) to peak identifiers.
        total_kmers : int or None
            Total number of k-mers in the dataset.
        total_genomic_kmers : int or None
            Total number of genomic (non-error) k-mers.
        total_error_kmers : int or None
            Total number of error k-mers.
        error_fraction : float or None
            Fraction of k-mers classified as errors.
    """
    cov_tab: pd.DataFrame = None
    cov2peak: dict = None
    total_kmers: int = None
    total_genomic_kmers: int = None
    total_error_kmers: int = None
    error_fraction: int = None

    def local_aggregation(self, distance: int, noise_filter: int, mask_errors: bool) -> None:
        """
        Cluster coverage points into peaks using Manhattan distance.
        
        Performs local aggregation of k-mer pairs based on their coverage
        coordinates. Points within the specified Manhattan distance are grouped
        into the same peak, with peak assignment based on highest-frequency
        neighbors.
        
        Args:
            distance : int
                Manhattan distance threshold for clustering coverage points.
            noise_filter : int
                Minimum frequency threshold; points below this are filtered as noise.
            mask_errors : bool
                If True, coverage points near the error line are marked as errors.
        
        Notes:
            Updates the `cov2peak` attribute with peak assignments. Error k-mers
            are assigned peak ID -1.
        """
        cov2freq, cov2peak = defaultdict(int), defaultdict(int)
        min_covB = min(self.cov_tab["covB"])
        
        next_peak = 1
        for row in self.cov_tab.itertuples("CovRow"):
            cov2freq[(row.covA, row.covB)] = (
                    row.freq  # a make a frequency dictionary on the fly, because I don't need any value that was not processed yet
                )
            if row.freq < noise_filter:
                break
            
            if mask_errors and row.covB < min_covB + distance:
                cov2peak[(row.covA, row.covB)] = -1
                continue
            
            peak_id = self._assign_peak(row, cov2freq, cov2peak, distance, next_peak)
            cov2peak[(row.covA, row.covB)] = peak_id

            if peak_id == next_peak:
                next_peak += 1
        
        self.cov2peak = cov2peak
    
    def _assign_peak(self, row, cov2freq: dict, cov2peak: dict, 
                     distance: int, next_peak: int) -> int:
        """
        Assign peak ID based on highest neighbor.
        """
        highest_neighbour_coords = (0, 0)
        highest_neighbour_freq = 0
        # for each kmer pair I will retrieve all neighbours (Manhattan distance)
        covB, covA, freq = row.covB, row.covA, row.freq
        for xA in range(covA - distance, covA + distance + 1):
            # for explored A coverage in neiborhood, we explore all possible B coordinates
            distanceB = distance - abs(covA - xA)
            for xB in range(covB - distanceB, covB + distanceB + 1):
                xB, xA = sorted([xA, xB])  # this is to make sure xB is smaller than xA
                # iterating only though those that were assigned already
                # and recording only the one with highest frequency
                if cov2peak[(xA, xB)] and cov2freq[(xA, xB)] > highest_neighbour_freq:
                    highest_neighbour_coords = (xA, xB)
                    highest_neighbour_freq = cov2freq[(xA, xB)]

        has_neighbour = highest_neighbour_freq > 0            
        if has_neighbour:
            return cov2peak[highest_neighbour_coords]
        else:
            return next_peak
            

    def peak_aggregation(self) -> None:
        """
        Assign peak identifiers to coverage table and sort.
        
        Adds a 'smudge' column to the coverage table containing peak identifiers
        from the cov2peak mapping, then sorts the table by covA and covB.
        """
        self.cov_tab["smudge"] = [
            self.cov2peak[(row.covA, row.covB)] for row in self.cov_tab.itertuples()
        ]
        self.cov_tab.sort_values(["covA", "covB"], ascending=True, inplace=True)

    def write_peaks(self) -> None:
        """
        Write peak-annotated coverage data to stderr.
        
        Outputs tab-separated values with columns: covB, covA, freq, peak.
        Calls peak_aggregation() first to ensure data is properly annotated.
        """
        self.peak_aggregation()
        for row in self.cov_tab.itertuples():
            logger.info(f"{row.covB}\t{row.covA}\t{row.freq}\t{row.peak}")

    def count_kmers(self) -> 'KmerStatistics':
        """
        Calculate k-mer statistics including total, genomic, and error counts.
        
        Computes summary statistics for k-mers in the dataset, distinguishing
        between genomic k-mers (smudge != -1), k-mers in identified smudges
        (smudge > 0), and error k-mers (smudge == -1).
        
        Notes:
            Updates the following attributes:
            - total_kmers : Total k-mer count
            - total_genomic_kmers : Genomic k-mer count (excluding errors)
            - total_genomic_kmers_in_smudges : K-mers in identified peaks
            - total_error_kmers : Error k-mer count
            - error_fraction : Proportion of error k-mers
        """
        self.peak_aggregation()

        if self.cov_tab.empty:
            raise InvalidCoverageDataError("Cannot count k-mers from empty coverage table")

        stats = KmerStatistics(
            total_kmers=self.cov_tab["freq"].sum(),
            genomic_kmers=self.cov_tab.loc[self.cov_tab["smudge"] != -1]["freq"].sum(),
            error_kmers=self.cov_tab.loc[self.cov_tab["smudge"] == -1]["freq"].sum(),
            genomic_kmers_in_smudges=self.cov_tab.loc[self.cov_tab["smudge"] > 0]["freq"].sum()
        )
    
        # Store for backward compatibility
        self.total_kmers = stats.total_kmers
        self.total_genomic_kmers = stats.genomic_kmers
        self.total_genomic_kmers_in_smudges = stats.genomic_kmers_in_smudges
        self.total_error_kmers = stats.error_kmers
        self.error_fraction = stats.error_fraction
    
        return stats

# ==============================================================================
# Class: KmerStatistics
# ==============================================================================

@dataclass
class KmerStatistics:
    """Statistics about k-mer coverage."""
    total_kmers: int
    genomic_kmers: int
    error_kmers: int
    genomic_kmers_in_smudges: int = 0
    
    @property
    def error_fraction(self) -> float:
        """Fraction of k-mers that are errors."""
        return self.error_kmers / self.total_kmers if self.total_kmers > 0 else 0.0
    
    @property
    def error_percentage(self) -> float:
        """Percentage of k-mers that are errors."""
        return self.error_fraction * 100
    
    @property
    def smudge_fraction(self) -> float:
        """Fraction of genomic k-mers in identified smudges."""
        return self.genomic_kmers_in_smudges / self.genomic_kmers if self.genomic_kmers > 0 else 0.0
    
    def __str__(self) -> str:
        return (
            f"K-mer Statistics:\n"
            f"  Total k-mers: {self.total_kmers:,}\n"
            f"  Genomic k-mers: {self.genomic_kmers:,}\n"
            f"  Error k-mers: {self.error_kmers:,} ({self.error_percentage:.2f}%)\n"
            f"  K-mers in smudges: {self.genomic_kmers_in_smudges:,} ({self.smudge_fraction*100:.2f}%)"
        )

# ==============================================================================
# CLASS: Smudges
# ==============================================================================

@dataclass
class Smudges:
    """
    Analyzer for k-mer coverage smudges (peaks) and genome structure inference.
    
    This class identifies and characterizes coverage smudges corresponding to
    different genomic structures (e.g., AAAB, ABB) and estimates haploid coverage.
    
    Parameters
    ----------
    cov_tab : pandas.DataFrame
        Coverage table with peak annotations.
    total_genomic_kmers : int
        Total number of genomic (non-error) k-mers.
    
    Attributes
    ----------
    cov_tab : pandas.DataFrame
        The coverage table data.
    total_genomic_kmers : int
        Total genomic k-mer count used for normalization.
    cov : float or None
        Estimated haploid coverage (1n coverage).
    centrality_df : pandas.DataFrame or None
        DataFrame of tested coverages and their centrality scores.
    final_smudge_container : dict or None
        Container for final smudge assignments.
    smudge_tab : pandas.DataFrame or None
        Summary table of identified smudges with sizes.
    """

    cov_tab: pd.DataFrame = None
    total_genomic_kmers: int = None
    cov: float = None
    centrality_df: pd.DataFrame = None
    final_smudge_container: pd.DataFrame = None
    smudge_tab: pd.DataFrame = None
    local_agg_smudge_container: pd.DataFrame = None

    def get_centrality_df(self, min_c: float, max_c: float, smudge_size_cutoff=AnalysisConfig.min_smudge_size_default) -> None:
        """
        Estimate optimal haploid coverage through iterative grid search.
        
        Performs a three-stage grid search with increasing precision to find
        the haploid coverage that best centers the observed smudges on their
        theoretical positions.
        
        Args:
            min_c : float
                Initial minimum coverage to test.
            max_c : float
                Initial maximum coverage to test.
            smudge_size_cutoff : float, optional
                Minimum relative size for smudges to be considered (default: 0.02).
        
        Notes:
            Updates the `cov` and `centrality_df` attributes. Progress is written
            to stderr showing the best coverage at each precision level.
        """
        grid_params = [(0.05, 0.05, 2), (-1.9, 1.9, 0.2), (-0.19, 0.19, 0.01)]
        results = []

        for i, params in enumerate(grid_params):
            cov_list = np.arange(int(min_c) + params[0], int(max_c) + params[1], params[2])
            best_cov, centralities = self.get_best_coverage(cov_list, smudge_size_cutoff)

            results.append({"covs": cov_list, "centralities": centralities, "best_cov": best_cov})

            min_c, max_c = best_cov, best_cov

            if i > 0:
                logger.info(f"Best coverage to precision of 1/{10**i}: {best_cov:.2f}")

        # just to be sure
        results[-1]["covs"] = np.append(results[-1]["covs"], results[-1]["best_cov"] / 2)
        best_cov, centralities = self.get_best_coverage(
            cov_list=results[-1]["covs"],
            smudge_size_cutoff=smudge_size_cutoff,
            centralities=results[-1]["centralities"],
            last_check=True,
        )

        logger.info(f"Best coverage to precision of 1/{10**i} (just to be sure): {best_cov:.2f}")

        self.cov = best_cov
        self.centrality_df = pd.DataFrame(
            {
                "coverage": np.concatenate([result["covs"] for result in results]),
                "centrality": np.concatenate([result["centralities"] for result in results]),
            }
        )

    def get_best_coverage(self, cov_list: ArrayLike, smudge_size_cutoff=AnalysisConfig.min_smudge_size_default, centralities=None, last_check=False) -> tuple[float, List[float]]:
        """
        Find coverage with minimum centrality from a list of candidates.
        
        Args:
            cov_list : array-like
                List of coverage values to test.
            smudge_size_cutoff : float, optional
                Minimum relative smudge size threshold (default: 0.02).
            centralities : list or None, optional
                Pre-computed centrality values (default: None).
            last_check : bool, optional
                If True, only test the last coverage value (default: False).
        
        Returns:
            best_cov : float
                Coverage value with minimum centrality.
            centralities : list
                List of centrality values for all tested coverages.
        """

        if centralities is None:
            centralities = []

        if last_check:
            to_test = [cov_list[-1]]
        else:
            to_test = cov_list

        for cov in to_test:
            smudge_container = self.get_smudge_container(cov, smudge_size_cutoff)
            centralities.append(get_centrality(smudge_container, cov))
        return cov_list[np.argmin(centralities)], centralities
    
# Too long 
    def get_smudge_container(self, cov: float, smudge_filter: float, method="fishnet") -> dict:
        """
        Extract and label smudges using specified method.
        
        Args:
            cov : float
                Haploid coverage estimate.
            smudge_filter : float
                Minimum relative size for smudges to be retained.
            method : {'fishnet', 'local_aggregation'}, optional
                Method for smudge extraction (default: 'fishnet').
                - 'fishnet': Uses theoretical coverage windows for each structure
                - 'local_aggregation': Uses pre-computed peak assignments
            
        Returns:
            smudge_container : dict of pandas.DataFrame
                Dictionary mapping structure labels (e.g., 'AAAB') to DataFrames
                containing the k-mer pairs belonging to that structure.
        """
        smudge_container = defaultdict(pd.DataFrame)

        if method == "fishnet":
            for Bs in range(1, 9):
                min_cov_Bs, max_cov_Bs = get_cov_limits(Bs, cov)

                cov_tab_isoB = self.cov_tab.loc[
                    (self.cov_tab["smudge"] != -1)  # this is removing the errors
                    & (self.cov_tab["covB"] > min_cov_Bs)  # bot of the fishnet for A
                    & (self.cov_tab["covB"] < max_cov_Bs)  # top of the fishnet for B
                ]

                for As in range(Bs, (AnalysisConfig.max_structure_range - Bs)):
                    min_cov_As, max_cov_As = get_cov_limits(As, cov)

                    cov_tab_iso_smudge = cov_tab_isoB.loc[
                        (cov_tab_isoB["covA"] > min_cov_As) & (cov_tab_isoB["covA"] < max_cov_As)
                    ]
                    if cov_tab_iso_smudge["freq"].sum() / self.total_genomic_kmers > smudge_filter:
                        if not smudge_container["A" * As + "B" * Bs].empty:
                            smudge_container["A" * As + "B" * Bs] = pd.concat(
                                [
                                    smudge_container["A" * As + "B" * Bs],
                                    cov_tab_iso_smudge,
                                ],
                                axis=0,
                                ignore_index=True,
                            )
                        else:
                            smudge_container["A" * As + "B" * Bs] = cov_tab_iso_smudge

        if method == "local_aggregation":
            peak = 1
            while peak <= max(self.cov_tab["smudge"]):
                cov_tab_smudge = self.cov_tab.loc[self.cov_tab["smudge"] == peak]
                covA, covB = get_centre_cov_by_mode(
                    cov_tab_smudge
                ) 
                As, Bs = round(covA / cov), round(covB / cov) # this is approximate labeling (modus of coverage / 1n coverage)

                if (cov_tab_smudge["freq"].sum() / self.total_genomic_kmers) > smudge_filter: # smudge will be added only if passes the filter
                    logger.debug(
                        f'Recording peak ({covA};{covB}) {peak} as {"A" * As}{"B" * Bs} smudge'
                    )
                    if not smudge_container["A" * As + "B" * Bs].empty:
                        smudge_container["A" * As + "B" * Bs] = pd.concat(
                            [smudge_container["A" * As + "B" * Bs], cov_tab_smudge],
                            axis=0,
                            ignore_index=True,
                        )
                    else:
                        smudge_container["A" * As + "B" * Bs] = cov_tab_smudge
                peak += 1
        return smudge_container

    def generate_smudge_table(self, smudge_container: dict) -> None:
        """
        Create summary table of smudge structures and sizes.
        
        Parameters
        ----------
        smudge_container : dict of pandas.DataFrame
            Dictionary mapping structure labels to their coverage data.
        
        Notes
        -----
        Updates the `smudge_tab` attribute with a DataFrame containing:
        - structure: Smudge label (e.g., 'AAAB')
        - size: Total k-mer count in the smudge
        - rel_size: Relative size (fraction of total genomic k-mers)
        """

        annotated_smudges = list(smudge_container.keys())
        smudge_sizes = [smudge_container[smudge]["freq"].sum() for smudge in annotated_smudges]
        smudge_sizes_rel = [round(smudge_size / self.total_genomic_kmers, 4) for smudge_size in smudge_sizes]
        self.smudge_tab = pd.DataFrame(
            {
                "structure": annotated_smudges,
                "size": smudge_sizes,
                "rel_size": smudge_sizes_rel,
            }
        )

    def centrality_plot(self, output: str, fmt="pdf") -> None:
        """
        Generate plot of coverage vs centrality scores.
        
        Parameters
        ----------
        output : str
            Output file prefix.
        fmt : str, optional
            Output format (default: 'pdf'). Can be any matplotlib-supported format.
        
        Notes
        -----
        Creates a plot showing how centrality varies with coverage, saved as
        '{output}_centralities.{fmt}'. Lower centrality indicates better fit.
        """
        fig, axs = plt.subplots(figsize=(8, 8))
        fontsize = 32
        plt.plot(
            self.centrality_df["coverage"],
            self.centrality_df["centrality"],
            "o",
            color="black",
            markersize=4,
        )
        axs.set_xlabel("Coverage")
        axs.set_ylabel("Centrality [(theoretical_center - actual_center) / coverage ]")
        fig.savefig(f"{output}_centralities.{fmt}")

# ==============================================================================
# CLASS: SmudgeplotData
# ==============================================================================
    
@dataclass
class SmudgeplotData:
    """
    Container for smudgeplot visualization data and parameters.
    
    Stores processed coverage and smudge data along with plotting parameters
    and file paths for output generation.
    
    Parameters
    ----------
    cov_tab : pandas.DataFrame
        Coverage table with columns for covA, covB, and freq.
    smudge_tab : pandas.DataFrame
        Smudge summary table with structure labels and sizes.
    cov : float
        Estimated haploid coverage.
    error_fraction : float, optional
        Fraction of k-mers classified as errors (default: 0).
    
    Attributes
    ----------
    cov_tab : pandas.DataFrame
        The coverage data.
    smudge_tab : pandas.DataFrame
        The smudge summary data.
    cov : float
        Haploid coverage estimate.
    error_fraction : float
        Error k-mer fraction.
    lims : dict
        Axis limits for plotting.
    error_string : str or None
        Formatted error percentage string.
    cov_string : str or None
        Formatted coverage string.
    fig_title : str or None
        Figure title with coverage and error info.
    linear_plot_file : str or None
        Output path for linear-scale plot.
    log_plot_file : str or None
        Output path for log-scale plot.
    json_report_file : str or None
        Output path for JSON report.
    """
    cov_tab: pd.DataFrame
    smudge_tab: pd.DataFrame
    cov: float
    error_fraction: float = 0.0
    lims: dict = field(default_factory=dict)
    error_string: Optional[str] = None
    cov_string: Optional[str] = None
    fig_title: Optional[str] = None
    linear_plot_file: Optional[str] = None
    log_plot_file: Optional[str] = None
    json_report_file: Optional[str] = None

    def calc_cov_columns(self) -> None:
        """
        Calculate derived coverage columns for plotting.
        
        Adds two columns to cov_tab:
        - total_pair_cov: Sum of covA and covB
        - minor_variant_rel_cov: Normalized minor coverage (covB / total)
        """
        self.cov_tab["total_pair_cov"] = self.cov_tab[["covA", "covB"]].sum(axis=1)
        self.cov_tab["minor_variant_rel_cov"] = self.cov_tab["covB"] / self.cov_tab["total_pair_cov"]

    def filter_cov_quant(self, cov_filter: float = None, quant_filter: float = None) -> None:
        """
        Filter coverage data by minimum coverage or quantile threshold.
        
        Parameters
        ----------
        cov_filter : float or None, optional
            Minimum coverage threshold; removes k-mers below this (default: None).
        quant_filter : float or None, optional
            Upper quantile threshold (0-100); removes k-mers above this percentile
            of total coverage (default: None).
        
        Notes
        -----
        Filters are applied to cov_tab in place. Both filters can be applied
        simultaneously.
        """
        if cov_filter:
            self.cov_tab = self.cov_tab.loc[
                (self.cov_tab["covA"] >= cov_filter) & (self.cov_tab["covB"] >= cov_filter)
            ]
        if quant_filter:
            # 0 < int(quant_filter) < 100
            upper_quantile = np.percentile(
                a=self.cov_tab["total_pair_cov"],
                q=quant_filter,
                weights=self.cov_tab["freq"],
                method="inverted_cdf",
            )
            self.cov_tab = self.cov_tab.loc[self.cov_tab["total_pair_cov"] < upper_quantile]

    def get_ax_lims(self, upper_ylim: float = None) -> None:
        """
        Calculate appropriate axis limits for plotting.
        
        Parameters
        ----------
        upper_ylim : float or None, optional
            Override for upper y-axis limit (default: None).
        
        Notes
        -----
        Updates the `lims` attribute with xlim (0 to 0.5) and ylim based on
        coverage range. The y-limit is adaptive based on whether coverage
        estimation was successful.
        """

        if self.cov == np.percentile(
            a=self.cov_tab["total_pair_cov"],
            q=AnalysisConfig.quantile_95,
            weights=self.cov_tab["freq"],
            method="inverted_cdf",
        ):
            self.lims["ylim"] = [
                min(self.cov_tab["total_pair_cov"]),
                max(self.cov_tab["total_pair_cov"]),
            ]
        else:
            self.lims["ylim"] = [
                min(self.cov_tab["total_pair_cov"]) - 1,  # or 0?
                min(max(100, AnalysisConfig.min_coverage_floor * self.cov), max(self.cov_tab["total_pair_cov"])),
            ]
        if upper_ylim:
            self.lims["ylim"][1] = upper_ylim
        self.lims["xlim"] = [0, 0.5]

    def def_strings(self, title: str = None, output: str = "smudgeplot", fmt: str = "pdf") -> None: 
        """
        Define output file paths and figure title.
        
        Parameters
        ----------
        title : str or None, optional
            Custom title for the plot (default: None).
        output : str, optional
            Output file prefix (default: 'smudgeplot').
        fmt : str, optional
            Output format (default: 'pdf').
        
        Notes
        -----
        Updates fig_title, linear_plot_file, log_plot_file, and
        json_report_file attributes.
        """

        if title:
            fig_title = str(title)
        else:
            fig_title = "NA"
        self.fig_title = f"{fig_title}\n1n = {self.cov:.0f}\nerr = {self.error_fraction*100:.2f}%"
        self.linear_plot_file = f"{output}_smudgeplot.{fmt}"
        self.log_plot_file = f"{output}_smudgeplot_log10.{fmt}"
        self.json_report_file = f"{output}_smudgeplot_report.json"

# ==============================================================================
# CLASS: CoverageValidator
# ==============================================================================

class CoverageValidator:
    """Validates coverage data integrity."""
    
    @staticmethod
    def validate_coverage_table(cov_tab: pd.DataFrame) -> None:
        """Validate coverage table meets requirements."""
        required_cols = {"covB", "covA", "freq"}

        if not required_cols.issubset(set(cov_tab.columns)):
            raise InvalidCoverageDataError("Coverage table must contain columns covB, covA, freq")

        if cov_tab.empty:
            raise InvalidCoverageDataError("Coverage table cannot be empty")

        try:
            if (cov_tab[list(required_cols)] < 0).any().any():
                raise InvalidCoverageDataError("Coverage values must be non-negative")
        except KeyError:
            raise InvalidCoverageDataError("Coverage table must contain columns covB, covA, freq")


# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

# Too long
def get_centrality(
    smudge_container: Dict[str, pd.DataFrame], 
    cov: float, 
    centre: str = "mode", 
    dist: str = "theoretical_center"
) -> float:
    """
    Calculate centrality score for a given coverage estimate.
    
    Centrality measures how well the observed smudge positions match their
    theoretical positions based on the proposed haploid coverage. Lower
    centrality indicates better fit.
    
    Parameters
    ----------
    smudge_container : dict of pandas.DataFrame
        Dictionary mapping structure labels to coverage data.
    cov : float
        Haploid coverage to test.
    centre : {'mode', 'mean'}, optional
        Method for calculating smudge centers (default: 'mode').
    dist : {'theoretical_center', 'theoretical_edge', 'empirical_edge'}, optional
        Distance metric to use (default: 'theoretical_center').
    
    Returns
    -------
    centrality : float
        Weighted average of distances from theoretical positions. Returns 1
        if no smudges are found.
    
    Notes
    -----
    The centrality is calculated as the frequency-weighted mean of distances
    between observed and expected positions for all smudges.
    """
    centralities: List[float] = []
    freqs: List[float] = []
    for smudge, smudge_tab in smudge_container.items():
        As = smudge.count("A")
        Bs = smudge.count("B")
        kmer_in_the_smudge = smudge_tab["freq"].sum()
        freqs.append(kmer_in_the_smudge)

        if centre == "mode":
            center_A, center_B = get_centre_cov_by_mode(smudge_tab)

        if centre == "mean":
            center_A = sum(smudge_tab["freq"] * smudge_tab["covA"]) / kmer_in_the_smudge
            center_B = sum(smudge_tab["freq"] * smudge_tab["covB"]) / kmer_in_the_smudge

        if dist == "empirical_edge":
            distA = min(
                [
                    abs(smudge_tab["covA"].max() - center_A),
                    abs(center_A - smudge_tab["covA"].min()),
                ]
            )
            distB = min(
                [
                    abs(smudge_tab["covB"].max() - center_B),
                    abs(center_B - smudge_tab["covB"].min()),
                ]
            )

        if dist == "theoretical_edge":
            distA = min(abs(center_A - (cov * (As - AnalysisConfig.max_manhattan_distance))), abs((cov * (As + AnalysisConfig.max_manhattan_distance)) - center_A))
            distB = min(abs(center_B - (cov * (Bs - AnalysisConfig.max_manhattan_distance))), abs((cov * (Bs + AnalysisConfig.max_manhattan_distance)) - center_B))

        if dist == "theoretical_center":
            distA = abs((center_A - (cov * As)) / cov)
            distB = abs((center_B - (cov * Bs)) / cov)

        logger.debug(f"Processing: {As}A{Bs}B; with center: {distA}, {distB} and joint freq {kmer_in_the_smudge}")
        centrality = distA + distB
        centralities.append(centrality)

    if len(centralities) == 0:
        return 1
    return fmean(centralities, weights=freqs)


def generate_plots(
    smudges: Smudges,
    coverages: Coverages,
    cov: float,
    smudge_size_cutoff: float,
    outfile: str,
    title: str,
    fmt: str = None,
    upper_ylim: str = None,
    json_report: bool = False,
    input_params: dict = None,
    palette: str = 'viridis',
    invert_cols: bool = False
) -> None:
    """
    Generate smudgeplot visualizations and optional JSON report.
    
    Creates both linear and log-scale smudgeplots showing k-mer coverage
    distributions and identified genome structures.
    
    Parameters
    ----------
    smudges : Smudges
        Smudges object with coverage analysis results.
    coverages : Coverages
        Coverages object with k-mer coverage data.
    cov : float
        Estimated haploid coverage.
    smudge_size_cutoff : float
        Minimum relative size for smudges to include.
    outfile : str
        Output file prefix.
    title : str
        Plot title.
    fmt : str or None, optional
        Output format (default: None).
    upper_ylim : float or None, optional
        Upper y-axis limit override (default: None).
    json_report : bool, optional
        Whether to generate JSON report (default: False).
    input_params : dict or None, optional
        Input parameters to include in JSON report (default: None).
    palette : str, optional
        Matplotlib colormap name (default: 'viridis').
    invert_cols : bool, optional
        Whether to invert the color palette (default: False).
    """
    smudges.fishnet_smudge_container = smudges.get_smudge_container(cov, smudge_size_cutoff, "fishnet")
    smudges.generate_smudge_table(smudges.fishnet_smudge_container)

    smudgeplot_data = SmudgeplotData(coverages.cov_tab, smudges.smudge_tab, cov, coverages.error_fraction)
    prepare_smudgeplot_data_for_plotting(smudgeplot_data, outfile, title, fmt=fmt, upper_ylim=upper_ylim)

    config = PlotConfig(palette=palette, invert_colours=invert_cols)
    smudgeplot(smudgeplot_data, log=False, config=config)
    smudgeplot(smudgeplot_data, log=True, config=config)

    if json_report:
        write_json_report(smudgeplot_data, input_params)

def write_json_report(smg_data: SmudgeplotData, 
                     input_params: Optional[Dict] = None, 
                     min_size: float = AnalysisConfig.min_report_size_default) -> None:
    """
    Write JSON report with smudgeplot analysis results.
    
    Parameters
    ----------
    smg_data : SmudgeplotData
        Smudgeplot data object with results.
    input_params : dict or None, optional
        Input parameters to include in report (default: None).
    min_size : float, optional
        Minimum size for smudges in top_smudges list (default: 0.03).
    
    Notes
    -----
    Creates JSON file at smg_data.json_report_file with version info,
    parameters, coverage estimate, error fraction, and detected smudges.
    """
    report = {
        "version": version("smudgeplot"),
        "commandline_arguments": sys.argv[1:],
        "input_parameters": input_params,
        "haploid_coverage": float(f"{smg_data.cov:.3f}"),
        "error_fraction": smg_data.error_fraction,
        "top_smudges": [
            {"structure": row.structure, "fraction": row.rel_size}
            for row in smg_data.smudge_tab.itertuples(index=False)
            if row.rel_size > min_size
        ],
        "smudges": [
            {"structure": row.structure, "count": row.size, "fraction": row.rel_size}
            for row in smg_data.smudge_tab.itertuples(index=False)
        ],
    }
    
    try:
        with open(smg_data.json_report_file, "w") as fh:
            json.dump(report, fh, indent=2)
            fh.write("\n")
        logger.info(f"JSON report written to {smg_data.json_report_file}")
    except IOError as e:
        logger.error(f"Failed to write JSON report: {e}")
        raise IOError(f"Failed to write JSON report to {smg_data.json_report_file}: {e}")


def prepare_smudgeplot_data_for_plotting(smudgeplot_data: SmudgeplotData, output: str, title: str, fmt: str = None, upper_ylim: float = None) -> None:
    """
    Prepare SmudgeplotData object for visualization.
    
    Calculates derived columns, applies filters, determines axis limits,
    and sets output file paths.
    
    Parameters
    ----------
    smudgeplot_data : SmudgeplotData
        Data object to prepare.
    output : str
        Output file prefix.
    title : str
        Plot title.
    fmt : str or None, optional
        Output format (default: None).
    upper_ylim : float or None, optional
        Upper y-axis limit (default: None).
    """

    smudgeplot_data.calc_cov_columns()
    smudgeplot_data.filter_cov_quant()
    smudgeplot_data.get_ax_lims(upper_ylim=upper_ylim)
    smudgeplot_data.def_strings(output=output, title=title, fmt=fmt)

def smudgeplot(data: SmudgeplotData, log: bool = False, config = None) -> None:
    """
    Create smudgeplot visualization.
    
    Generates a comprehensive plot showing k-mer pair coverage distribution,
    marginal histograms, identified genome structures, and a legend.
    
    Parameters
    ----------
    data : SmudgeplotData
        Prepared smudgeplot data object.
    log : bool, optional
        Use log10 scale for color mapping (default: False).
    palette : str, optional
        Matplotlib colormap name (default: 'viridis').
    invert_cols : bool, optional
        Invert color palette (default: False).
    
    Notes
    -----
    Creates a multi-panel figure with:
    - Main panel: 2D coverage plot with color-coded frequency
    - Top panel: Histogram of normalized minor coverage
    - Right panel: Histogram of total coverage
    - Legend panel: Color scale and smudge size table
    """

    #setup
    if config is None:
        config = PlotConfig()
    smudgeplot_data = copy.deepcopy(data)
    fig, axes = _setup_smudgeplot_figure(smudgeplot_data.fig_title, config)
    colour_ramp = get_col_ramp(config.palette, delay=16 if log else 0, invert_cols=config.invert_colours)

    #plot components
    _plot_main_smudges(data, axes['main'], colour_ramp, config, log)
    _plot_legend( max(data.cov_tab["freq"]), axes['legend'], colour_ramp, config, log)
    _plot_smudge_sizes(data.smudge_tab, data.cov, data.error_string, axes['size'])
    if data.cov > 0:
        _plot_expected_haplotype_structure(data, axes['main'], config)
    if config.show_histograms:
        _plot_histograms(data, axes, config, log)

    # save
    outfile = data.log_plot_file if log else data.linear_plot_file
    fig.savefig(outfile, dpi=config.dpi)
    plt.close()

def _setup_smudgeplot_figure(fig_title, config) -> Tuple[plt.Figure, Dict[str, plt.Axes]]:
    """Define smudgeplot layout"""
    fig, ((top_ax, legend_ax), (main_ax, size_ax)) = plt.subplots(
        nrows=2, ncols=2, 
        width_ratios=[3, 1], 
        height_ratios=[1, 3], 
        figsize=config.figsize
    )
    size_ax.sharey(main_ax)
    top_ax.sharex(main_ax)
    legend_ax.axis("off")
    size_ax.axis("off")
    top_ax.axis("off")
    plt.subplots_adjust(wspace=0.05, hspace=0.05)
    
    top_ax.set_title(fig_title, 
        fontsize=32, 
        loc="left", 
        y=1.0, 
        pad=-14, 
        weight="bold",
    )

    return fig, {
        'top': top_ax,
        'legend': legend_ax,
        'main': main_ax,
        'size': size_ax
    }

def _plot_histograms(data: SmudgeplotData, axes: dict, config: PlotConfig, log: bool):
    """
    Plot marginal histograms.

    Args:
        data : SmudgeplotData
            Smudgeplot data object.
        axes : Dict
            Dict of matplotlib.axes.Axes with location keys
        config : PlotConfig object
            Data object containing default plotting parameters.
        log : bool, optional
            Whether to use log scale (default: False).

    """
    cov_tab = data.cov_tab

    # Right histogram - total coverage of kmer pair
    plot_hist(data = cov_tab['total_pair_cov'], 
              ax = axes['size'], 
              orientation = 'horizontal',
              bins = max(cov_tab['total_pair_cov']) - min(cov_tab['total_pair_cov']),
              weights = cov_tab['freq'],
              log = log
    )

    # Top histogram - normalised minor kmer coverage
    plot_hist(data = cov_tab["minor_variant_rel_cov"], 
              ax = axes['top'],
              weights= cov_tab['freq'],
              bins = config.hist_bins_minor,
              log = log
    )

def plot_hist(data: ArrayLike, ax: mpl.axes.Axes, weights: ArrayLike, orientation: str ='vertical', bins: int = 50, log: bool = False) -> None:
    """
    Plot histogram on specified axis.
    
    Args:
        data : array-like
            Data values to histogram.
        ax : matplotlib.axes.Axes
            Axis to plot on.
        weights : array-like
            Weights for each data point (typically frequencies).
        orientation : {'vertical', 'horizontal'}, optional
            Histogram orientation (default: 'vertical').
        bins : int, optional
            Number of bins (default: 50).
        log : bool, optional
            Whether to use log scale (default: False).
    """

    ax.hist(data,
        weights = weights,
        bins = bins,
        color = 'firebrick',
        edgecolor='firebrick',
        orientation = orientation
    )

def _plot_main_smudges(data: SmudgeplotData, ax:  mpl.axes.Axes, colour_ramp: list, config: PlotConfig, log: bool =False) -> None:
    """
    Plot k-mer coverage smudges as colored rectangles.
    
    Creates the main 2D coverage plot with color-coded frequencies, where
    each k-mer pair is represented as a rectangle colored by its frequency.
    
    Args:
        data : SmudgeplotData
            Smudgeplot data object.
        ax : matplotlib.axes.Axes
            Axis to plot on.
        colour_ramp : list
            List of hex color codes for frequency gradient.
        config : PlotConfig object
            Data object containing default plotting parameters.
        log : bool, optional
            Use log10 scale for frequencies (default: False).
    """
    cov_tab = data.cov_tab
    mask = cov_tab["covA"] == cov_tab["covB"]
    cov_tab.loc[mask, "freq"] = cov_tab[mask]["freq"] * 2

    if log:
        cov_tab["freq"] = np.log10(cov_tab["freq"])

    cov_tab["col"] = [
        str(colour_ramp[int(i)])
        for i in round((len(colour_ramp) - 1) * cov_tab["freq"] / max(cov_tab["freq"]))
    ]

    ax.plot()
    ax.set_xlim(data.lims["xlim"])
    ax.set_ylim(data.lims["ylim"])
    ax.set_xlabel("Normalized minor kmer coverage: B / (A + B)", fontsize=config.fontsize)
    ax.set_ylabel("Total coverage of the kmer pair: A + B", fontsize=config.fontsize)
    ax.tick_params(axis="both", labelsize=config.legend_fontsize)
    ax.spines[["right", "top"]].set_visible(False)

    min_cov_to_plot = max(data.lims["ylim"][0], min(cov_tab["total_pair_cov"]))

    patches_nested = [
        get_one_coverage(cov_tab, cov, ax) for cov in np.arange(min_cov_to_plot, data.lims["ylim"][1])
    ]
    patches_flat = [x for xs in patches_nested for x in xs]
    ax.add_collection(PatchCollection(patches_flat, match_original=True))


def get_one_coverage(cov_tab: pd.DataFrame, cov:float, ax: mpl.patches.Rectangle) -> list:
    """
    Get patches for all k-mer pairs at a specific coverage.
    
    Args:
        cov_tab : pandas.DataFrame
            Coverage table.
        cov : float
            Total coverage value to extract.
        ax : matplotlib.axes.Axes
            Axis object (used for coordinate system).
    
    Returns:
        patches : list of matplotlib.patches.Rectangle
            List of colored rectangle patches.
    """
    cov_rows = cov_tab[cov_tab["total_pair_cov"] == cov]
    
    if cov_rows.empty:
        return []
    
    width = 1 / (2 * cov)
    lefts = cov_rows["minor_variant_rel_cov"] - width
    rights = np.minimum(
        AnalysisConfig.rel_coverage_tolerance, 
        cov_rows["minor_variant_rel_cov"] + width
    )
    patches = [
        mpl.patches.Rectangle(
            (left, cov - 0.5), right - left, 1,
            linewidth=1, edgecolor=col, facecolor=col
        )
        for left, right, col in zip(
            lefts.values, rights.values, cov_rows["col"].values
        )
    ]
    return patches

def get_one_box(left: float, right: float, cov: float, colour: str, ax: mpl.axes.Axes) -> mpl.patches.Rectangle:
    """
    Create a single colored rectangle patch.
    
    Args:
        left : float
            Left x-coordinate (normalized minor coverage).
        right : float
            Right x-coordinate.
        cov : float
            Total coverage (y-coordinate center).
        colour : str
            Hex color code.
        ax : matplotlib.axes.Axes
            Axis object.
    
    Returns:
        rectangle : matplotlib.patches.Rectangle
            Rectangle patch object.
    """
    width = float(right) - float(left)
    return mpl.patches.Rectangle(
        (float(left), cov - 0.5),
        width,
        1,
        linewidth=1,
        edgecolor=colour,
        facecolor=colour,
    )

def _plot_expected_haplotype_structure(data: SmudgeplotData, ax: mpl.axes.Axes, config: PlotConfig) -> None:
    """
    Places text labels (e.g., 'AAAB', '3A2B') at the positions
    of identified genome structures.
    
    Args:
        data : SmudgeplotData
            Smudgeplot data object.
        ax : matplotlib.axes.Axes
            Axis to annotate.
        config : PlotConfig object
            Data object containing default plotting parameters.
    
    Returns:
        Only labels smudges with relative size > 0.05.
    """
    smudge_tab = data.smudge_tab.copy(deep=True)
    smudge_tab.loc[:, "ploidy"] = smudge_tab["structure"].str.len()
    smudge_tab = smudge_tab.loc[smudge_tab["rel_size"] > 0.05]
    smudge_tab.loc[:, "corrected_minor_variant_cov"] = (
        smudge_tab["structure"].str.count("B") / smudge_tab["ploidy"]
    )
    smudge_tab.loc[:, "label"] = reduce_structure_representation(smudge_tab["structure"])

    for index, row in smudge_tab.iterrows():

        if (smudge_tab["corrected_minor_variant_cov"][index] == 0.5) and config.adjust:
            ha = "right"
        else:
            ha = "center"

        x = row["corrected_minor_variant_cov"]
        y = row["ploidy"] * data.cov
        ax.text(x, y, row["label"], fontsize=config.label_fontsize, va="center_baseline", ha=ha)

def _plot_legend(kmer_max: float, ax: mpl.axes.Axes, colour_ramp: list, config: PlotConfig , log=False) -> None:
    """
    Create colour scale legend.
    
    Args:
        kmer_max : float
            Maximum k-mer frequency.
        ax : matplotlib.axes.Axes
            Axis for legend.
        colour_ramp : list
            List of hex color codes.
        config : PlotConfig object
            Data object containing default plotting parameters.
        log : bool, optional
            Whether scale is logarithmic (default: False).
    """
    if log:    
        ax.set_title("log kmer pairs\n", ha="center", fontsize=28, weight="bold")
        colour_ramp = colour_ramp[16:]
        for i, colour in enumerate(colour_ramp):
            rect = mpl.patches.Rectangle(
                (0, ((i + 2) - 0.01) / 18),
                0.5,
                ((i + 2) + 0.99) / 18,
                linewidth=1,
                edgecolor=colour,
                facecolor=colour,
            )
            ax.add_patch(rect)
    
        for i, n in enumerate([0, 3.5, 4, 4.5, 5, 5.5, 6]):
            ax.text(
                0.66,
                i / 6,
                str(rounding(10 ** (kmer_max * n / 6))),
                fontsize=config.legend_fontsize,
            )
    else:
        ax.set_title("kmer pairs\n", ha="center", fontsize=28, weight="bold")
        for i, colour in enumerate(colour_ramp):
            rect = mpl.patches.Rectangle(
                (0, ((i + 1) - 0.01) / 33),
                0.5,
                ((i + 1) + 0.99) / 33,
                linewidth=1,
                edgecolor=colour,
                facecolor=colour,
            )
            ax.add_patch(rect)

        for i in range(7):
            ax.text(0.66, i / 6, str(rounding(kmer_max * i / 6)), fontsize=config.legend_fontsize)

def _plot_smudge_sizes(smudge_tab: pd.DataFrame, cov: float, error_string: str, ax: mpl.axes.Axes, min_size: float = AnalysisConfig.min_report_size_default) -> None:
    """
    Display table of smudge sizes.
    
    Args:
        smudge_tab : pandas.DataFrame
            Smudge summary table.
        cov : float
            Haploid coverage.
        error_string : str
            Error information string.
        ax : matplotlib.axes.Axes
            Axis to place text.
        min_size : float, optional
            Minimum relative size to display (default: 0.03).
    """

    ax.plot()
    ax.set_title("")

    if cov > 0:
        size_tuples = sorted(
            [
                (smudge, size)
                for smudge, size in zip(
                    reduce_structure_representation(smudge_tab["structure"]).to_list(),
                    round(smudge_tab["rel_size"], 2),
                    strict=True,
                )
            ],
            key=lambda x: x[1],
            reverse=True,
        )
        labels = [f"{size:>3,.2f}   {smudge:<6s}" for smudge, size in size_tuples if size >= min_size]
        label_string = "\n".join(labels)
    else:
        label_string = error_string

    ax.text(0.1, 1, 
            label_string, 
            ha="left", va="top", 
            fontsize=28, 
            transform=ax.transAxes,
        )


def reduce_structure_representation(smudge_labels: pd.Series) -> pd.Series:
    """
    Converts labels like 'AAAAABBB' to '5A3B' for better readability.
    
    Args:
        smudge_labels : pandas.Series
            Series of structure labels.
    
    Returns:
        reduced_labels : pandas.Series
            Series with shortened labels where appropriate.
    
    Notes:
        Only shortens labels longer than 4 characters.
        Some redundancy with smudge2short() ?
    """

    structures_to_adjust = smudge_labels.str.len() > 4
    if not any(structures_to_adjust):
        return smudge_labels
    else:
        As = smudge_labels[structures_to_adjust].str.count("A").map(str)
        Bs = smudge_labels[structures_to_adjust].str.count("B").map(str)
        new_labels = smudge_labels.copy(deep=True)
        new_labels[structures_to_adjust] = As + "A" + Bs + "B"
        return new_labels

def generate_smudge_report(smudges: Smudges, coverages: Coverages, cov: float, args: argparse.Namespace, smudge_size_cutoff: float, print_header: bool) -> None:
    """
    Generate and write smudge analysis report.
    
    Args:
        smudges : Smudges
            Smudges object with analysis results.
        coverages : Coverages
            Coverages object with k-mer statistics.
        cov : float
            Haploid coverage.
        args : argparse.Namespace
            Command-line arguments.
        smudge_size_cutoff : float
            Minimum smudge size threshold.
        print_header : bool
            Whether to print TSV header row.
    """

    smudges.generate_smudge_table(smudges.local_agg_smudge_container)

    logger.info(
        f"Detected smudges / sizes:\n"
        f"  {smudges.smudge_tab['structure'].to_list()}\n"
        f"  {smudges.smudge_tab['size'].to_list()}"
    )

    write_smudge_report(smudges, coverages, cov, args, print_header=print_header)


def write_smudge_report(smudges: Smudges, coverages: Coverages, cov: float, args: argparse.Namespace, print_header: bool) -> None:
    """
    Write tab-separated smudge report file.
    
    Creates a TSV file with k-mer counts and smudge counts for all
    possible structures up to a maximum ploidy.
    
    Args:
        smudges : Smudges
            Smudges object with results.
        coverages : Coverages
            Coverages object with k-mer statistics.
        cov : float
            Haploid coverage.
        args : argparse.Namespace
            Command-line arguments containing input file and output prefix.
        print_header : bool
            Whether to include column headers.
    """

    smudge_dict = create_smudge_dict(AnalysisConfig.min_ploidy, AnalysisConfig.max_ploidy)

    dataset = args.infile.split("/")[-1]
    meta_df = pd.DataFrame.from_dict(
        {
            "dataset": [dataset],
            "total_kmers": [coverages.total_kmers],
            "total_error_kmers": [coverages.total_error_kmers],
        }
    )

    smudges.smudge_tab.loc[:, "label"] = reduce_structure_representation(smudges.smudge_tab["structure"])
    for row in smudges.smudge_tab.itertuples():
        if smudge_dict.get(row.label, "Missing") != "Missing":
            smudge_dict[row.label] = [row.size]
        else:
            logger.info(f"Unexpected smudge label {row.label} excluded from smudge report")

    smudge_df = pd.DataFrame.from_dict(smudge_dict).fillna(0)
    out_df = pd.concat([meta_df, smudge_df], axis=1)

    out_df.to_csv(f"{args.o}.smudge_report.tsv", sep="\t", index=False, header=print_header)
    logger.info(f"Written smudge report to: {args.o}.smudge_report.tsv")


def create_smudge_dict(min_ploidy: int, max_ploidy: int) -> tuple:
    """
    Create dictionary template for all possible smudge structures (e.g., AB, AAB, ABB, AAAB)
    up in the specified ploidy range.
    
    Returns:
        A tuple of:
            smudge_dict : dict
                Dictionary with structure labels as keys.
            sorted_smudges : pandas.Series
                Sorted list of structure labels.
    
    Notes:
        - Structures are created following the pattern where Bs <= As and
        As + Bs <= 2 * max_ploidy.
    """
    smudge_list = []
    for Bs in range(min_ploidy, max_ploidy + 1):
        for As in range(Bs, ((2 * max_ploidy) + 1 - Bs)):
            smudge_list.append("A" * As + "B" * Bs)

    smudge_list.sort()
    sorted_smudges = sorted(smudge_list, key=len)
    smudges_rr = reduce_structure_representation(pd.Series(sorted_smudges))
    smudge_dict = dict.fromkeys(smudges_rr, np.nan)
    return smudge_dict

def local_min(ys: list) -> list:
    """
    Find indices of local minima in a list.
    
    Args:
        ys : list
            List of numeric values.
    
    Returns
        indices : list of int
            Indices where local minima occur.
    
    Notes:
        - A point is considered a local minimum if it's less than or equal to its
        predecessor and strictly less than its successor.
        - Taken from https://stackoverflow.com/a/29614335
    """
    return [
        i
        for i, y in enumerate(ys)
        if ((i == 0) or (ys[i - 1] >= y)) and ((i == len(ys) - 1) or (y < ys[i + 1]))
    ]


def round_up_nice(x: float) -> int:
    """
    Round number for plotting.
    
    Rounds to the nearest multiple of 10^(digits-1) for single-digit numbers,
    or 10^(digits-2) for larger numbers.
    
    Args:
        x : float
            Number to round.
        
    Returns:
        rounded : int
            Rounded value.
    """
    digits = ceil(log(x, 10))
    if digits <= 1:
        multiplier = 10 ** (digits - 1)
    else:
        multiplier = 10 ** (digits - 2)
    return ceil(x / multiplier) * multiplier


def cutoff(kmer_hist: list, boundary: str = 'L') -> None:
    """
    Prints coverage cutoff boundary from k-mer histogram.
    
    Args:
        kmer_hist : list of str
            K-mer histogram.
        boundary : str
            'L' for lower cutoff, 
            'U' for upper cutoff. (USE NOT RECOMMENDED)
    
    Notes:
        If boundary == 'L', finds first local minimum and returns 1.25x that value.
        If boundary == 'U' finds 99.8th percentile of k-mer coverage.
    """
    hist = [int(line.split()[1]) for line in kmer_hist]
    if boundary == "L":
        local_minima = local_min(hist)[0]
        L = max(10, int(round(local_minima * LOCAL_MIN_MULTIPLIER)))
        logger.info(f"{L}")
    elif boundary == 'U':
        logger.info(
            "Warning: We discourage using the original hetmer algorithm.\n"
        )
        # take 99.8 quantile of kmers that are more than one in the read set
        number_of_kmers = np.sum(hist[1:])
        hist_rel_cumsum = [np.sum(hist[1 : i + 1]) / number_of_kmers for i in range(1, len(hist))]
        min(range(len(hist_rel_cumsum)))
        U = round_up_nice(min([i for i, q in enumerate(hist_rel_cumsum) if q > QUANTILE_998]))
        logger.info(f"{U}")


def load_hetmers(file_h: str) -> pd.DataFrame:
    """
    Load heterozygous k-mer coverage file.
    
    Args:
        file_h : str
            File containing tab-separated coverage data.
    
    Returns:
        cov_tab : pandas.DataFrame
            Coverage table sorted by descending frequency:
            ['covB', 'covA', 'freq'].
    """

    try:
        cov_tab = pd.read_csv(file_h, names=["covB", "covA", "freq"], 
                              sep="\t", comment='#')
        CoverageValidator.validate_coverage_table(cov_tab)
        return cov_tab.sort_values("freq", ascending=False)
    except FileNotFoundError:
        logger.error(f"Coverage file not found: {file_h}")
        raise
    except pd.errors.EmptyDataError:
        raise InvalidCoverageDataError("Coverage file is empty")

def get_centre_cov_by_mode(smudge_tab: pd.DataFrame) -> tuple[float, float]:
    """
    Find center coordinates of a single smudge using the mode.
    """
    centre = smudge_tab.loc[smudge_tab["freq"].idxmax()]
    return centre["covA"], centre["covB"]


def get_cov_limits(Xs: int, cov: float) -> tuple[float, float]:
    """
    Calculate coverage limits.
    """
    min_cov = 0 if Xs == 1 else cov * (Xs - 0.5)
    max_cov = cov * (Xs + 0.5)
    return min_cov, max_cov


def get_col_ramp(col_ramp: str = 'viridis', delay: int = 0, invert_cols: bool = False) -> List[str]:
    """
    Generate color ramp for frequency visualization.
    
    Args:
        col_ramp : str, optional
            Matplotlib colormap name (default: 'viridis').
        delay : int, optional
            Number of low colors to extend at the start (default: 0).
        invert_cols : bool, optional
            Reverse the colormap (default: False).
    
    Returns:
        ramp : list of str
            List of hex color codes.
    """
    if invert_cols:
        col_ramp += "_r"
    cmap = plt.get_cmap(col_ramp, 32 - int(delay))
    ramp = [mpl.colors.rgb2hex(cmap(i)) for i in range(cmap.N)]
    ramp = [ramp[0]] * delay + ramp
    return ramp

def rounding(number: float) -> int:
    """
    Rounded value (to nearest 1000, 100, or 10 depending on magnitude).
    """

    if number > 1000:
        return round(number / 1000) * 1000
    elif number > 100:
        return round(number / 100) * 100
    else:
        return round(number / 10) * 10
        
def smudge2short(smudge_long: str) -> str:
    """
    Convert full smudge label to short format.
    Example: 'AAAABBB' -> 4A3B
    """
    return(str(smudge_long.count('A')) + 'A' + str(smudge_long.count('B')) + 'B')