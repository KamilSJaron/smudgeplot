#!/usr/bin/env python3

import sys
import numpy as np
from pandas import read_csv  # type: ignore
from pandas import DataFrame  # type: ignore
from pandas import Series  # type: ignore
from pandas import concat  # type: ignore
from numpy import arange
from numpy import argmin
from numpy import concatenate
from math import log
from math import ceil
from statistics import fmean
from collections import defaultdict
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from itertools import product


class Coverages:
    def __init__(self, cov_tab):
        self.cov_tab = cov_tab
        self.cov2peak = {}
        self.total_kmers = None
        self.total_genomic_kmers = None
        self.total_error_kmers = None
        self.error_fraction = None

    def local_aggregation(self, distance, noise_filter, mask_errors):
        """
        Generate a dictionary that gives us for each combination of coverages a frequency
        """
        cov2freq = defaultdict(int)
        cov2peak = defaultdict(int)

        L = min(self.cov_tab["covB"])  # important only when --mask_errors is on

        ### clustering
        next_peak = 1
        for idx, covB, covA, freq in self.cov_tab.itertuples():
            cov2freq[(covA, covB)] = (
                freq  # a make a frequency dictionary on the fly, because I don't need any value that was not processed yet
            )
            if freq < noise_filter:
                break
            highest_neigbour_coords = (0, 0)
            highest_neigbour_freq = 0
            # for each kmer pair I will retrieve all neibours (Manhattan distance)
            for xA in range(covA - distance, covA + distance + 1):
                # for explored A coverage in neiborhood, we explore all possible B coordinates
                distanceB = distance - abs(covA - xA)
                for xB in range(covB - distanceB, covB + distanceB + 1):
                    xB, xA = sorted([xA, xB])  # this is to make sure xB is smaller than xA
                    # iterating only though those that were assigned already
                    # and recording only the one with highest frequency
                    if cov2peak[(xA, xB)] and cov2freq[(xA, xB)] > highest_neigbour_freq:
                        highest_neigbour_coords = (xA, xB)
                        highest_neigbour_freq = cov2freq[(xA, xB)]
            if highest_neigbour_freq:  # > 0:
                cov2peak[covA, covB] = cov2peak[highest_neigbour_coords]
            else:
                if mask_errors and covB < L + distance:
                    cov2peak[(covA, covB)] = -1  # error line
                else:
                    cov2peak[(covA, covB)] = (
                        next_peak  # if I want to keep info about all locally agregated smudges
                    )
                    next_peak += 1
        self.cov2peak = cov2peak

    def peak_aggregation(self):
        self.cov_tab["smudge"] = [
            self.cov2peak[(covA, covB)] for idx, covB, covA, freq in self.cov_tab.itertuples()
        ]
        self.cov_tab.sort_values(["covA", "covB"], ascending=True, inplace=True)

    def write_peaks(self):
        self.peak_aggregation()
        for idx, covB, covA, freq, peak in self.cov_tab.itertuples():
            sys.stdout.write(f"{covB}\t{covA}\t{freq}\t{peak}\n")
        sys.stdout.flush()

    def count_kmers(self):
        self.peak_aggregation()

        self.total_kmers = self.cov_tab["freq"].sum()
        self.total_genomic_kmers = self.cov_tab.loc[self.cov_tab["smudge"] != -1]["freq"].sum()
        self.total_genomic_kmers_in_smudges = self.cov_tab.loc[self.cov_tab["smudge"] > 0]["freq"].sum()
        self.total_error_kmers = self.cov_tab.loc[self.cov_tab["smudge"] == -1]["freq"].sum()
        self.error_fraction = self.total_error_kmers / self.total_kmers


class Smudges:
    def __init__(self, cov_tab, total_genomic_kmers):
        self.cov_tab = cov_tab
        self.total_genomic_kmers = total_genomic_kmers
        self.cov = None
        self.centrality_df = None
        self.final_smudge_container = None
        self.smudge_tab = None

    def get_centrality_df(self, min_c, max_c, smudge_size_cutoff=0.02):
        grid_params = [(0.05, 0.05, 2), (-1.9, 1.9, 0.2), (-0.19, 0.19, 0.01)]
        results = []

        for i, params in enumerate(grid_params):
            cov_list = arange(int(min_c) + params[0], int(max_c) + params[1], params[2])
            best_cov, centralities = self.get_best_coverage(cov_list, smudge_size_cutoff)

            results.append({"covs": cov_list, "centralities": centralities, "best_cov": best_cov})

            min_c, max_c = best_cov, best_cov

            if i > 0:
                sys.stderr.write(f"Best coverage to precision of 1/{10**i}: {best_cov:.2f}\n")

        # just to be sure
        results[-1]["covs"] = np.append(results[-1]["covs"], results[-1]["best_cov"] / 2)
        best_cov, centralities = self.get_best_coverage(
            cov_list=results[-1]["covs"],
            smudge_size_cutoff=smudge_size_cutoff,
            centralities=results[-1]["centralities"],
            last_check=True,
        )

        sys.stderr.write(f"Best coverage to precision of 1/{10**i} (just to be sure): {best_cov:.2f}\n")

        self.cov = best_cov
        self.centrality_df = DataFrame(
            {
                "coverage": concatenate([result["covs"] for result in results]),
                "centrality": concatenate([result["centralities"] for result in results]),
            }
        )

    def get_best_coverage(self, cov_list, smudge_size_cutoff=0.02, centralities=None, last_check=False):
        if centralities is None:
            centralities = []

        if last_check:
            to_test = [cov_list[-1]]
        else:
            to_test = cov_list

        for cov in to_test:
            smudge_container = self.get_smudge_container(cov, smudge_size_cutoff)
            centralities.append(get_centrality(smudge_container, cov))
        return cov_list[argmin(centralities)], centralities

    def get_smudge_container(self, cov, smudge_filter, method="fishnet"):
        smudge_container = defaultdict(DataFrame)

        if method == "fishnet":
            for Bs in range(1, 9):
                min_cov_Bs, max_cov_Bs = get_cov_limits(Bs, cov)

                cov_tab_isoB = self.cov_tab.loc[
                    (self.cov_tab["smudge"] != -1)  # this is removing the errors
                    & (self.cov_tab["covB"] > min_cov_Bs)  # bot of the fishnet for A
                    & (self.cov_tab["covB"] < max_cov_Bs)  # top of the fishnet for B
                ]

                for As in range(Bs, (17 - Bs)):
                    min_cov_As, max_cov_As = get_cov_limits(As, cov)

                    cov_tab_iso_smudge = cov_tab_isoB.loc[
                        (cov_tab_isoB["covA"] > min_cov_As) & (cov_tab_isoB["covA"] < max_cov_As)
                    ]
                    if cov_tab_iso_smudge["freq"].sum() / self.total_genomic_kmers > smudge_filter:
                        if not smudge_container["A" * As + "B" * Bs].empty:
                            smudge_container["A" * As + "B" * Bs] = concat(
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
                )  ## fmean(cov_tab_smudge["covA"], cov_tab_smudge["freq"])
                As, Bs = round(covA / cov), round(covB / cov)

                # if As==0 or Bs==0: # skip < 1ns
                #    print(f'Peak {peak}, As={As};Bs={Bs}')
                #    peak += 1
                #    continue

                if (cov_tab_smudge["freq"].sum() / self.total_genomic_kmers) > smudge_filter:
                    # sys.stderr.write(
                    #     f'Recording peak ({covA};{covB}) {peak} as {"A" * As}{"B" * Bs} smudge\n'
                    # )
                    if not smudge_container["A" * As + "B" * Bs].empty:
                        smudge_container["A" * As + "B" * Bs] = concat(
                            [smudge_container["A" * As + "B" * Bs], cov_tab_smudge],
                            axis=0,
                            ignore_index=True,
                        )
                    else:
                        smudge_container["A" * As + "B" * Bs] = cov_tab_smudge
                peak += 1

            # smudge_container = {smudge:cov_tab for smudge, cov_tab in smudge_container.items() if (cov_tab["freq"].sum() / self.total_genomic_kmers) > smudge_filter}
        return smudge_container

    def generate_smudge_table(self, smudge_container):
        annotated_smudges = list(smudge_container.keys())
        smudge_sizes = [smudge_container[smudge]["freq"].sum() for smudge in annotated_smudges]
        smudge_sizes_rel = [round(smudge_size / self.total_genomic_kmers, 4) for smudge_size in smudge_sizes]
        self.smudge_tab = DataFrame(
            {
                "structure": annotated_smudges,
                "size": smudge_sizes,
                "rel_size": smudge_sizes_rel,
            }
        )

    def centrality_plot(self, output):
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
        fig.savefig(output + "_centralities.pdf")


class SmudgeplotData:
    def __init__(self, cov_tab, smudge_tab, cov, error_fraction=0):
        self.cov_tab = cov_tab
        self.smudge_tab = smudge_tab
        self.cov = cov
        self.error_fraction = error_fraction
        self.lims = {}
        self.error_string = None
        self.cov_string = None
        self.fig_title = None
        self.linear_plot_file = None
        self.log_plot_file = None

    def calc_cov_columns(self):
        self.cov_tab["total_pair_cov"] = self.cov_tab[["covA", "covB"]].sum(axis=1)
        self.cov_tab["minor_variant_rel_cov"] = self.cov_tab["covB"] / self.cov_tab["total_pair_cov"]

    def filter_cov_quant(self, cov_filter=None, quant_filter=None):
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

    def get_ax_lims(self, upper_ylim=None):
        if self.cov == np.percentile(
            a=self.cov_tab["total_pair_cov"],
            q=95,
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
                min(max(100, 10 * self.cov), max(self.cov_tab["total_pair_cov"])),
            ]
        if upper_ylim:
            self.lims["ylim"][1] = upper_ylim
        self.lims["xlim"] = [0, 0.5]

    def def_strings(self, title=None, output="smudgeplot"):
        # self.error_string = f"err = {round(self.error_fraction*100, 1)} %"
        # self.cov_string = f"1n = {self.cov}"
        if title:
            fig_title = str(title)
        else:
            fig_title = "NA"
        self.fig_title = f"{fig_title}\n\n1n = {self.cov:.0f}\nerr = {self.error_fraction*100:.2f}%"
        self.linear_plot_file = output + "_smudgeplot_py.pdf"
        self.log_plot_file = output + "_smudgeplot_log10_py.pdf"


def get_centrality(smudge_container, cov, centre="mode", dist="theoretical_center"):
    # sys.stderr.write(f"\tTesting coverage: {cov}\n\n")
    centralities = []
    freqs = []
    for smudge, smudge_tab in smudge_container.items():
        As = smudge.count("A")
        Bs = smudge.count("B")
        kmer_in_the_smudge = smudge_tab["freq"].sum()
        freqs.append(kmer_in_the_smudge)

        if centre == "mode":
            center_A, center_B = get_centre_cov_by_mode(smudge_tab)

        if centre == "mean":
            center_A = sum((smudge_tab["freq"] * smudge_tab["covA"])) / kmer_in_the_smudge
            center_B = sum((smudge_tab["freq"] * smudge_tab["covB"])) / kmer_in_the_smudge

        if dist == "empirical_edge":
            distA = min(
                [
                    abs(smudge_tab["covA"].max() - center["covA"]),
                    abs(center["covA"] - smudge_tab["covA"].min()),
                ]
            )
            distB = min(
                [
                    abs(smudge_tab["covB"].max() - center["covB"]),
                    abs(center["covB"] - smudge_tab["covB"].min()),
                ]
            )

        if dist == "theoretical_edge":
            distA = min(abs(center["covA"] - (cov * (As - 0.5))), abs((cov * (As + 0.5)) - center["covA"]))
            distB = min(abs(center["covB"] - (cov * (Bs - 0.5))), abs((cov * (Bs + 0.5)) - center["covB"]))

        if dist == "theoretical_center":
            distA = abs((center_A - (cov * As)) / cov)
            distB = abs((center_B - (cov * Bs)) / cov)

        # sys.stderr.write(f"Processing: {As}A{Bs}B; with center: {distA}, {distB} and joint freq {kmer_in_the_smudge}\n")
        centrality = distA + distB
        centralities.append(centrality)

    if len(centralities) == 0:
        return 1
    return fmean(centralities, weights=freqs)


def generate_plots(smudges, coverages, cov, smudge_size_cutoff, outfile, title):
    smudges.fishnet_smudge_container = smudges.get_smudge_container(cov, smudge_size_cutoff, "fishnet")
    smudges.generate_smudge_table(smudges.fishnet_smudge_container)

    smudgeplot_data = SmudgeplotData(coverages.cov_tab, smudges.smudge_tab, cov, coverages.error_fraction)
    prepare_smudgeplot_data_for_plotting(smudgeplot_data, outfile, title)

    smudgeplot(smudgeplot_data, log=False)
    smudgeplot(smudgeplot_data, log=True)


def prepare_smudgeplot_data_for_plotting(smudgeplot_data, output, title):

    smudgeplot_data.calc_cov_columns()
    smudgeplot_data.filter_cov_quant()
    smudgeplot_data.get_ax_lims()
    smudgeplot_data.def_strings(output=output, title=title)


def smudgeplot(data, log=False):  # I think user arguments need to be passed here
    cov_tab = data.cov_tab.copy(deep=True)
    smudge_tab = data.smudge_tab
    cov = data.cov
    lims = (
        data.lims
    )  # so things like lims can be set using user defined plotting parameters (I imagine palette will be the same although I did not try that one yet)
    fig_title = data.fig_title

    fig, axs = plt.subplots(nrows=2, ncols=2, width_ratios=[3, 1], height_ratios=[1, 3], figsize=(20, 20))
    plt.subplots_adjust(wspace=0.05, hspace=0.05)
    fontsize = 32

    main_ax = axs[1][0]
    legend_ax = axs[0][1]
    size_ax = axs[1][1]
    top_ax = axs[0][0]

    legend_ax.axis("off")
    size_ax.axis("off")
    top_ax.axis("off")

    if log:
        colour_ramp = get_col_ramp(delay=16)
        outfile = data.log_plot_file
    else:
        colour_ramp = get_col_ramp()
        outfile = data.linear_plot_file

    plot_smudges(cov_tab, colour_ramp, main_ax, lims, log=log, fontsize=fontsize)

    if cov > 0:
        plot_expected_haplotype_structure(smudge_tab, cov, main_ax, adjust=True, xmax=0.49)

    plot_legend(ax=legend_ax, kmer_max=max(cov_tab["freq"]), colour_ramp=colour_ramp, log=log)

    plot_hists = False
    if plot_hists:
        ## Not properly aligned with smudgeplot axes
        #right hist - total coverage of kmer pair
        plot_hist(data = cov_tab['total_pair_cov'], ax = size_ax, weights = np.nan(),orientation = 'horizontal')
        #upper hist - normalised minor kmer coverage
        plot_hist(data = cov_tab["minor_variant_rel_cov"], ax = top_ax, weights = np.nan())

    plot_smudge_sizes(smudge_tab, cov, data.error_string, size_ax)

    top_ax.set_title(fig_title, fontsize=42, loc="left", y=1.0, pad=-14, weight="bold")

    fig.savefig(outfile, dpi=200)
    plt.close()

def plot_hist(data, ax, weights, orientation='vertical'):
    ax.hist(data,
        weights = weights,
        bins = 100,
        color = 'firebrick',
        edgecolor='firebrick',
        orientation = orientation
    )

def plot_smudges(cov_tab, colour_ramp, ax, lims, log=False, fontsize=14):
    mask = cov_tab["covA"] == cov_tab["covB"]
    cov_tab.loc[mask, "freq"] = cov_tab[mask]["freq"] * 2

    if log:
        cov_tab["freq"] = np.log10(cov_tab["freq"])

    cov_tab["col"] = [
        str(colour_ramp[int(i)])
        for i in round((len(colour_ramp) - 1) * cov_tab["freq"] / max(cov_tab["freq"]))
    ]

    ax.plot()
    ax.set_xlim(lims["xlim"])
    ax.set_ylim(lims["ylim"])
    ax.set_xlabel("Normalized minor kmer coverage: B / (A + B)", fontsize=fontsize)
    ax.set_ylabel("Total coverage of the kmer pair: A + B", fontsize=fontsize)
    ax.tick_params(axis="both", labelsize=20)
    ax.spines[["right", "top"]].set_visible(False)

    min_cov_to_plot = max(lims["ylim"][0], min(cov_tab["total_pair_cov"]))

    patches_nested = [
        get_one_coverage(cov_tab, cov, ax) for cov in np.arange(min_cov_to_plot, lims["ylim"][1])
    ]
    patches_flat = [x for xs in patches_nested for x in xs]
    ax.add_collection(PatchCollection(patches_flat, match_original=True))


def get_one_coverage(cov_tab, cov, ax):
    cov_row_to_plot = cov_tab[cov_tab["total_pair_cov"] == cov]
    width = 1 / (2 * cov)
    lefts = (cov_row_to_plot["minor_variant_rel_cov"] - width).to_numpy()
    rights = (cov_row_to_plot["minor_variant_rel_cov"].apply(lambda x: min(0.5, x + width))).to_numpy()
    cols = cov_row_to_plot["col"].to_numpy()
    return [get_one_box(left, right, cov, col, ax) for left, right, col in zip(lefts, rights, cols)]


def get_one_box(left, right, cov, colour, ax):
    width = float(right) - float(left)
    return mpl.patches.Rectangle(
        (float(left), cov - 0.5),
        width,
        1,
        linewidth=1,
        edgecolor=colour,
        facecolor=colour,
    )


def plot_expected_haplotype_structure(smudge_tab, cov, ax, adjust=False, xmax=0.49):
    smudge_tab = smudge_tab.copy(deep=True)
    smudge_tab.loc[:, "ploidy"] = smudge_tab["structure"].str.len()
    smudge_tab = smudge_tab.loc[smudge_tab["rel_size"] > 0.05]
    smudge_tab.loc[:, "corrected_minor_variant_cov"] = (
        smudge_tab["structure"].str.count("B") / smudge_tab["ploidy"]
    )
    smudge_tab.loc[:, "label"] = reduce_structure_representation(smudge_tab["structure"])

    for index, row in smudge_tab.iterrows():

        if (smudge_tab["corrected_minor_variant_cov"][index] == 0.5) & adjust:
            ha = "right"
        else:
            ha = "center"

        x = row["corrected_minor_variant_cov"]
        y = row["ploidy"] * cov
        ax.text(x, y, row["label"], fontsize=28, va="center_baseline", ha=ha)


def plot_smudge_sizes(smudge_tab, cov, error_string, ax, min_size=0.03):
    ax.plot()
    ax.set_title("")
    if cov > 0:

        size_tuples = sorted(
            [
                (smudge, size)
                for smudge, size in zip(
                    reduce_structure_representation(smudge_tab["structure"]).to_list(),
                    round(smudge_tab["rel_size"], 2),
                )
            ],
            key=lambda x: x[1],
            reverse=True,
        )

        labels = [f"{size:>3,.2f}   {smudge:<6s}" for smudge, size in size_tuples if size >= min_size]

        label_string = "\n".join(labels)

        ax.text(0, 1, label_string, ha="left", va="top", fontsize=28, transform=ax.transAxes)

    else:
        ax.text(0, 1, error_string, ha="left", va="top", fontsize=28, transform=ax.transAxes)


def reduce_structure_representation(smudge_labels):
    structures_to_adjust = smudge_labels.str.len() > 4
    if not any(structures_to_adjust):
        return smudge_labels
    else:
        As = smudge_labels[structures_to_adjust].str.count("A").map(str)
        Bs = smudge_labels[structures_to_adjust].str.count("B").map(str)
        new_labels = smudge_labels.copy(deep=True)
        new_labels[structures_to_adjust] = As + "A" + Bs + "B"
        return new_labels


def plot_legend(ax, kmer_max, colour_ramp, log=False):
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
                fontsize=20,
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
            ax.text(0.66, i / 6, str(rounding(kmer_max * i / 6)), fontsize=20)


def generate_smudge_report(smudges, coverages, cov, args, smudge_size_cutoff, print_header):
    smudges.local_agg_smudge_container = smudges.get_smudge_container(
        cov, smudge_size_cutoff, "local_aggregation"
    )
    smudges.generate_smudge_table(smudges.local_agg_smudge_container)

    sys.stderr.write(
        f'Detected smudges / sizes : \n\
        \t {smudges.smudge_tab["structure"].to_list()} \n\
        \t {smudges.smudge_tab["size"].to_list()} \n'
    )

    write_smudge_report(smudges, coverages, cov, args, print_header=print_header)


def write_smudge_report(smudges, coverages, cov, args, print_header):
    smudge_dict, sorted_smudges = create_smudge_dict(16)

    dataset = args.infile.split("/")[-1]
    meta_df = DataFrame.from_dict(
        {
            "dataset": [dataset],
            "total_kmers": [coverages.total_kmers],
            "total_error_kmers": [coverages.total_error_kmers],
        }
    )

    smudges.smudge_tab.loc[:, "label"] = reduce_structure_representation(smudges.smudge_tab["structure"])
    for idx, smudge, size, rel_size, label in smudges.smudge_tab.itertuples():
        if smudge_dict.get(label, "Missing") != "Missing":
            smudge_dict[label] = [size]
        else:
            sys.stdout.write(f"Unexpected smudge label {label} excluded from smudge report\n")

    smudge_df = DataFrame.from_dict(smudge_dict).fillna(0)
    out_df = concat([meta_df, smudge_df], axis=1)

    out_df.to_csv(args.o + ".smudge_report.tsv", sep="\t", index=False, header=print_header)
    sys.stderr.write(f"Written smudge report to: {dataset.split('.')[0]}.smudge_report.tsv\n")


def create_smudge_dict(max_ploidy):
    smudge_list = []
    for Bs in range(1, max_ploidy + 1):
        for As in range(Bs, ((2 * max_ploidy) + 1 - Bs)):
            smudge_list.append("A" * As + "B" * Bs)

    smudge_list.sort()
    sorted_smudges = sorted(smudge_list, key=len)
    smudges_rr = reduce_structure_representation(Series(sorted_smudges))
    smudge_dict = dict.fromkeys(smudges_rr, np.nan)
    return smudge_dict, smudges_rr


# taken from https://stackoverflow.com/a/29614335
def local_min(ys):
    return [
        i
        for i, y in enumerate(ys)
        if ((i == 0) or (ys[i - 1] >= y)) and ((i == len(ys) - 1) or (y < ys[i + 1]))
    ]


def round_up_nice(x):
    digits = ceil(log(x, 10))
    if digits <= 1:
        multiplier = 10 ** (digits - 1)
    else:
        multiplier = 10 ** (digits - 2)
    return ceil(x / multiplier) * multiplier


def cutoff(kmer_hist, boundary):
    hist = [int(line.split()[1]) for line in kmer_hist]
    if boundary == "L":
        local_minima = local_min(hist)[0]
        L = max(10, int(round(local_minima * 1.25)))
        sys.stdout.write(str(L))
    else:
        sys.stderr.write(
            "Warning: We discourage using the original hetmer algorithm.\n\tThe updated (recommended) version does not take the argument U\n"
        )
        # take 99.8 quantile of kmers that are more than one in the read set
        number_of_kmers = np.sum(hist[1:])
        hist_rel_cumsum = [np.sum(hist[1 : i + 1]) / number_of_kmers for i in range(1, len(hist))]
        min(range(len(hist_rel_cumsum)))
        U = round_up_nice(min([i for i, q in enumerate(hist_rel_cumsum) if q > 0.998]))
        sys.stdout.write(str(U))
    sys.stdout.flush()


def load_hetmers(file_h):
    cov_tab = read_csv(file_h, names=["covB", "covA", "freq"], sep="\t")
    return cov_tab.sort_values("freq", ascending=False)


def get_centre_cov_by_mode(smudge_tab):
    centre = smudge_tab.loc[smudge_tab["freq"].idxmax()]
    return centre["covA"], centre["covB"]


def get_cov_limits(Xs, cov):
    min_cov = 0 if Xs == 1 else cov * (Xs - 0.5)
    max_cov = cov * (Xs + 0.5)
    return min_cov, max_cov


def get_col_ramp(col_ramp="viridis", delay=0, invert_cols=False):
    if invert_cols:
        col_ramp += "_r"
    cmap = plt.get_cmap(col_ramp, 32 - int(delay))
    ramp = [mpl.colors.rgb2hex(cmap(i)) for i in range(cmap.N)]
    ramp = [ramp[0]] * delay + ramp
    return ramp


def rounding(number):
    if number > 1000:
        return round(number / 1000) * 1000
    elif number > 100:
        return round(number / 100) * 100
    else:
        return round(number / 10) * 10
