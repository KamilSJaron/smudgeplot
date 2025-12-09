#!/usr/bin/env python3

import argparse
import os
import shutil
import sys
import logging
from importlib.metadata import version
from pathlib import Path
import numpy as np
import smudgeplot.smudgeplot as smg
from smudgeplot.config import PlotConfig, AnalysisConfig

def get_binary_path(name: str) -> str:
    """
    Find the path to a bundled binary.

    Searches in order:
    1. Bundled binaries in package (for wheel installs)
    2. System PATH (for development or conda installs)

    Args:
        name: Name of the binary (e.g., 'hetmers', 'extract_kmer_pairs')

    Returns:
        Full path to the binary

    Raises:
        FileNotFoundError: If binary cannot be found
    """
    # First, check for bundled binary in package
    package_dir = Path(__file__).parent
    bundled_binary = package_dir / "bin" / name

    if bundled_binary.exists() and os.access(bundled_binary, os.X_OK):
        return str(bundled_binary)

    # Fall back to system PATH
    system_binary = shutil.which(name)
    if system_binary:
        return system_binary

    raise FileNotFoundError(
        f"Binary '{name}' not found. Please ensure smudgeplot is properly installed. "
        f"Checked locations:\n"
        f"  - Package: {bundled_binary}\n"
        f"  - System PATH: (not found)\n"
        f"\nYou may need to reinstall smudgeplot or install the binaries manually."
    )


def run_binary(name: str, args: str) -> int:
    """
    Run a binary with the given arguments.

    Args:
        name: Name of the binary
        args: Space-separated argument string

    Returns:
        Return code from the binary
    """
    binary_path = get_binary_path(name)
    cmd = f"{binary_path} {args}"
    logger.info(f"Calling: {name} {args}")
    return os.system(cmd)


class Parser:
    def __init__(self):
        self.arguments = None
        argparser = argparse.ArgumentParser(
            # description='Inference of ploidy and heterozygosity structure using whole genome sequencing data',
            usage="""
            smudgeplot <task> [options]

            tasks: cutoff            Calculate meaningful values for lower kmer histogram cutoff.
                   hetmers           Calculate unique kmer pairs from a FastK k-mer database.
                   peak_aggregation  Agregates smudges using local aggregation algorithm; prints assignments to stdout.
                   plot              Generate 2d histogram; infer ploidy and plot a smudgeplot.
                   all               Runs all the steps (with default options)
                   extract           Extract kmer pair sequences from a FastK k-mer database.
            """
        )
        argparser.add_argument(
            "task",
            help="Task to execute; for task specific options execute smudgeplot <task> -h",
        )
        argparser.add_argument(
            "-v",
            "--version",
            action="store_true",
            default=False,
            help="Print the version and exit.",
        )
        # print version is a special case
        if len(sys.argv) > 1:
            if sys.argv[1] in ["-v", "--version"]:
                self.task = "version"
                return
            # the following line either prints help and die; or assign the name of task to variable task
            self.task = argparser.parse_args([sys.argv[1]]).task
        else:
            self.task = ""
        # if the task is known (i.e. defined in this file);
        if hasattr(self, self.task):
            # load arguments of that task
            getattr(self, self.task)()
        else:
            argparser.print_usage()
            if self.task == "":
                logger.info("No task provided")
            else:
                logger.info('"' + self.task + '" is not a valid task name')
            exit(1)

    def cutoff(self):
        """
        Calculate meaningful values for lower kmer histogram cutoff.
        """
        argparser = argparse.ArgumentParser(
            prog="smudgeplot cutoff",
            description="Calculate meaningful values for lower kmer histogram cutoff.",
        )
        argparser.add_argument(
            "infile",
            type=argparse.FileType("r"),
            help='Name of the input kmer histogram file (default "kmer.hist")."',
        )
        argparser.add_argument("boundary", help="Which bounary to compute L (lower) or U (upper).")
        self.arguments = argparser.parse_args(sys.argv[2:])

    def hetmers(self):
        """
        Calculate the frequencies of unique kmer pairs (hetmers) from a FastK database.
        """
        argparser = argparse.ArgumentParser(
            prog="smudgeplot hetmers",
            description="Calculate unique kmer pairs from FastK k-mer database.",
        )
        argparser.add_argument("infile", help="Input FastK database (.ktab) file.")
        argparser.add_argument(
            "-L",
            help="Count threshold below which k-mers are considered erroneous.",
            type=int,
        )
        argparser.add_argument("-t", help="Number of threads (default 4).", type=int, default=4)
        argparser.add_argument(
            "-o",
            help="The pattern used to name the output (kmerpairs).",
            default="kmerpairs",
        )
        argparser.add_argument(
            "-tmp",
            help="Directory where all temporary files will be stored (default /tmp).",
            default=".",
        )
        argparser.add_argument("--verbose", action="store_true", default=False, help="Verbose mode.")
        self.arguments = argparser.parse_args(sys.argv[2:])

    def peak_aggregation(self):
        """
        Aggregate k-mer pairs by local maxima. Locality is definted by user defined manhattan distance (-d) and the individual smudges are reported as incremeting indicies order by heights of peaks (not sizes of smudges). Unassigned k-mer pairs are labelled as 0.
        """
        argparser = argparse.ArgumentParser(
            prog="smudgeplot peak_aggregation",
            description="Aggregates smudges using local aggregation algorithm.")
        argparser.add_argument(
            "infile",
            help="Name of the input smu file with covarages and frequencies.",
        )
        argparser.add_argument(
            "-nf",
            "-noise_filter",
            help="k-mer pairs with frequencies lower than this value will not be aggregated into smudges.",
            type=int,
            default=50,
        )
        argparser.add_argument(
            "-d",
            "-distance",
            help="Manthattan distance of k-mer pairs that are considered neighbouring for the local aggregation purposes.",
            type=int,
            default=5,
        )
        argparser.add_argument(
            "--mask_errors",
            help="All k-mer pairs belonging to smudges with the peak distant less than -d from the error line will be labeled as -1 (errors).",
            action="store_true",
            default=False,
        )
        argparser.add_argument("-title", help="name printed at the top of the smudgeplot (default: infile prefix).", default=None)
        self.arguments = argparser.parse_args(sys.argv[2:])

    def extract(self):
        """
        Extract kmer pair sequences from a FastK k-mer database.
        """
        argparser = argparse.ArgumentParser(
            prog="smudgeplot extract",
            description="Extract kmer pair sequences from a FastK k-mer database.",
        )
        argparser.add_argument("infile", help="Input FastK database (.ktab) file.")
        argparser.add_argument("sma", help="Input annotated k-mer pair file (.sma).")
        argparser.add_argument("-t", help="Number of threads (default 4)", type=int, default=4)
        argparser.add_argument(
            "-o",
            help="The pattern used to name the output (kmerpairs).",
            default="kmerpairs",
        )
        argparser.add_argument(
            "-tmp",
            help="Directory where all temporary files will be stored (default /tmp).",
            default=".",
        )
        argparser.add_argument("--verbose", action="store_true", default=False, help="verbose mode")
        self.arguments = argparser.parse_args(sys.argv[2:])

    def plot(self):
        """
        Given coverage, infer ploidy and plot a smudgeplot.
        """
        argparser = argparse.ArgumentParser(
            prog="smudgeplot plot",
            description="Generate 2d histogram; infer ploidy and plot a smudgeplot."
        )
        argparser.add_argument("infile", help="Mame of the input tsv file with coverages and frequencies.")
        argparser.add_argument(
            "smudgefile",
            help="Name of the input tsv file with sizes of individual smudges.",
        )
        argparser.add_argument("n", help="The expected haploid coverage.", type=float)
        argparser.add_argument(
            "-o",
            help="The pattern used to name the output (smudgeplot).",
            default="smudgeplot",
        )

        argparser = self.add_plotting_arguments(argparser)

        self.arguments = argparser.parse_args(sys.argv[2:])

    def all(self):
        argparser = argparse.ArgumentParser(
            prog="smudgeplot all",
            description="Runs all the steps (with default options).")
        argparser.add_argument(
            "infile",
            help="Name of the input tsv file with covarages and frequencies.",
        )
        argparser.add_argument(
            "-o",
            help="The pattern used to name the output (smudgeplot).",
            default="smudgeplot",
        )
        argparser.add_argument("-cov_min", help="Minimal coverage to explore (default 6)", default=6)
        argparser.add_argument("-cov_max", help="Maximal coverage to explore (default 100)", default=100)
        argparser.add_argument(
            "-cov",
            help="The assumed coverage (no inference of 1n coverage is made).",
            type=float, # this is funny, it seems like the interface rejects floats although here it all looks correct
            default=0.0,
        )
        argparser.add_argument(
            "-d",
            "-distance",
            help="Manthattan distance of k-mer pairs that are considered neighbouring for local aggregation purposes.",
            type=int,
            default=2,
        )
        argparser = self.add_plotting_arguments(argparser)

        self.arguments = argparser.parse_args(sys.argv[2:])

    def add_plotting_arguments(self, argparser):
        argparser.add_argument(
            "-t",
            "--title",
            help="name printed at the top of the smudgeplot (default: infile prefix).",
            default=None,
        )
        argparser.add_argument(
            "-ylim",
            help="The upper limit for the coverage sum (the y axis)",
            type=int,
            default=None,
        )
        argparser.add_argument(
            "-col_ramp",
            help='Palette used for the plot (default "viridis", other sensible options are "magma", "mako" or "grey.colors" - recommended in combination with --invert_cols).',
            default="viridis",
        )
        argparser.add_argument(
            "--invert_cols",
            action="store_true",
            default=False,
            help="Invert the colour palette (default False).",
        )
        argparser.add_argument(
            "--format",
            default="png",
            help="Output format for the plots (default png)",
            choices=["pdf", "png", "svg"],
        )
        argparser.add_argument(
            "--json_report",
            action="store_true",
            default=False,
            help="Generate a JSON format report alongside the plots (default False)",
        )
        return argparser

def main():

    logging.basicConfig(level=logging.INFO,
                    format="%(message)s",
                    handlers=[
                    logging.StreamHandler(sys.stderr)
                    ]
    )
    logger = logging.getLogger(__name__)
    _parser = Parser()

    smdg_v = version("smudgeplot")
    logger.info(f"Running smudgeplot v{smdg_v}")
    if _parser.task == "version":
        exit(0)

    logger.info("Task: " + _parser.task)

    args = _parser.arguments

    if _parser.task == "cutoff":
        smg.cutoff(args.infile, args.boundary)

    if _parser.task == "hetmers":
        # PloidyPlot is expected to be installed in the system
        plot_args = " -o" + str(args.o)
        plot_args += " -e" + str(args.L)
        plot_args += " -T" + str(args.t)
        if args.verbose:
            plot_args += " -v"
        if args.tmp != ".":
            plot_args += " -P" + args.tmp
        plot_args += " " + args.infile

        run_binary("hetmers", plot_args)

    if _parser.task == "extract":
        plot_args = " -o" + str(args.o)
        plot_args += " -T" + str(args.t)
        if args.verbose:
            plot_args += " -v"
        if args.tmp != ".":
            plot_args += " -P" + args.tmp
        plot_args += " " + args.infile
        if args.sma.endswith(".sma"):
            plot_args += " " + args.sma.removesuffix(".sma")
        else:
            plot_args += " " + args.sma

        run_binary("extract_kmer_pairs", plot_args)

    if args.title:
        title=args.title
    else:
        title = ".".join(args.infile.split("/")[-1].split(".")[0:2])

    if _parser.task == "plot":
        smudge_tab = smg.read_csv(args.smudgefile, sep="\t", names=["structure", "size", "rel_size"])
        cov_tab = smg.load_hetmers(args.infile)
        smudgeplot_data = smg.SmudgeplotData(cov_tab, smudge_tab, args.n)
        smg.prepare_smudgeplot_data_for_plotting(smudgeplot_data, args.o, title, upper_ylim=args.ylim, fmt=args.format)
        config = PlotConfig(palette=args.col_ramp, invert_cols=args.invert_cols)
        smg.smudgeplot(smudgeplot_data, config, log=False)
        smg.smudgeplot(smudgeplot_data, config, log=True)

    # test for existence of smudge file
    if not os.path.exists(args.infile):
        logger.info(f"The input file {args.infile} not found. Please provide a valid smudge file.")

    logger.info("\nLoading data")
    coverages = smg.Coverages()
    coverages.cov_tab = smg.load_hetmers(args.infile)
    logger.info("\nMasking errors using local aggregation algorithm")

    if _parser.task == "peak_aggregation":

        coverages.local_aggregation(distance=args.d, noise_filter=args.nf, mask_errors=args.mask_errors)
        coverages.write_peaks()

    if _parser.task == "all":
        coverages.local_aggregation(distance=args.d, noise_filter=AnalysisConfig.task_all_noise_filter, mask_errors=True)
        stats = coverages.count_kmers()
        logger.info(
            stats
        )

        smudges = smg.Smudges(coverages.cov_tab, coverages.total_genomic_kmers)

        if args.cov == 0.0:
            logger.info("\nInferring 1n coverage using grid algorithm")

            smudges.get_centrality_df(args.cov_min, args.cov_max, AnalysisConfig.smudge_size_cutoff)
            np.savetxt(
                args.o + "_centralities.txt",
                np.around(smudges.centrality_df, decimals=6),
                fmt="%.4f",
                delimiter="\t",
            )

            if coverages.error_fraction < AnalysisConfig.error_limit:
                cov = smudges.cov
            else:
                cov = 0

            logger.info("\nCreating centrality plot")
            smudges.centrality_plot(args.o, args.format)
            logger.info(f"\nInferred coverage: {cov:.3f}")

        else:
            cov = args.cov
            logger(f"\nUser defined coverage: {cov:.3f}")

        logger.info("\nCreating smudge report")

        smudges.local_agg_smudge_container = smudges.get_smudge_container(cov, AnalysisConfig.smudge_size_cutoff, "local_aggregation")
        annotated_smudges = list(smudges.local_agg_smudge_container.keys())
        with open(args.o + ".sma", "w") as annotated_smu:
            annotated_smu.write("covB\tcovA\tfreq\tsmudge\n")
            for smudge in annotated_smudges:
                formated_smudge = smg.smudge2short(smudge)
                for row in smudges.local_agg_smudge_container[smudge].itertuples():
                    covB, covA, freq = row.covB, row.covA, row.freq
                    annotated_smu.write(f"{covB}\t{covA}\t{freq}\t{formated_smudge}\n")

        smg.generate_smudge_report(smudges, coverages, cov, args, AnalysisConfig.smudge_size_cutoff, print_header=True)
        logger.info("\nCreating smudgeplots")
        smg.generate_plots(
            smudges,
            coverages,
            cov,
            AnalysisConfig.smudge_size_cutoff,
            args.o,
            title,
            fmt=args.format,
            upper_ylim=args.ylim,
            json_report=args.json_report,
            input_params=vars(args),
            palette=args.col_ramp,
            invert_cols=args.invert_cols
        )

if __name__ == "__main__":
    main()
