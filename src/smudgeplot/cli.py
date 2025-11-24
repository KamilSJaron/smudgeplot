#!/usr/bin/env python3

import argparse
import sys
from importlib.metadata import version
from os import system
import numpy as np
import smudgeplot.smudgeplot as smg


class Parser:
    def __init__(self):
        self.arguments = None
        argparser = argparse.ArgumentParser(
            # description='Inference of ploidy and heterozygosity structure using whole genome sequencing data',
            usage="""
            smudgeplot <task> [options]

            tasks: cutoff            Calculate meaningful values for lower kmer histogram cutoff.
                   hetmers           Calculate unique kmer pairs from a FastK k-mer database.
                   peak_aggregation  Aggregates smudges using local aggregation algorithm.
                   plot              Generate 2d histogram; infer ploidy and plot a smudgeplot.
                   all               Runs all the steps (with default options).
            """
        )
        # removing this for now;
        #        extract   Extract kmer pairs within specified coverage sum and minor covrage ratio ranges
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
            if self.task != "":
                sys.stderr.write('"' + self.task + '" is not a valid task name\n')
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
            help="Name of the input tsv file with covarages and frequencies.",
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
        argparser.add_argument("-cov_min", help="Minimal coverage to explore (default 6).", default=6)
        argparser.add_argument("-cov_max", help="Maximal coverage to explore (default 60).", default=60)
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
            choices=["pdf", "png"],
        )
        argparser.add_argument(
            "--json_report",
            action="store_true",
            default=False,
            help="Generate a JSON format report alongside the plots (default False)",
        )
        return argparser

def fin():
    sys.stderr.write("\nDone!\n")
    exit(0)

def main():
    _parser = Parser()

    smdg_v = version("smudgeplot")
    sys.stderr.write(f"Running smudgeplot v{smdg_v}\n")
    if _parser.task == "version":
        exit(0)

    sys.stderr.write("Task: " + _parser.task + "\n")

    args = _parser.arguments

    if args.title:
        title=args.title
    else:
        title = ".".join(args.infile.split("/")[-1].split(".")[0:2])

    if _parser.task == "cutoff":
        smg.cutoff(args.infile, args.boundary)
        fin()

    if _parser.task == "hetmers":
        # PloidyPlot is expected ot be installed in the system as well as the R library supporting it
        plot_args = " -o" + str(args.o)
        plot_args += " -e" + str(args.L)
        plot_args += " -T" + str(args.t)
        if args.verbose:
            plot_args += " -v"
        if args.tmp != ".":
            plot_args += " -P" + args.tmp
        plot_args += " " + args.infile

        sys.stderr.write("Calling: hetmers (PloidyPlot kmer pair search) " + plot_args + "\n")
        system("hetmers " + plot_args)

        fin()

    if _parser.task == "plot":
        smudge_tab = smg.read_csv(args.smudgefile, sep="\t", names=["structure", "size", "rel_size"])
        cov_tab = smg.load_hetmers(args.infile)
        smudgeplot_data = smg.SmudgeplotData(cov_tab, smudge_tab, args.n)
        smg.prepare_smudgeplot_data_for_plotting(smudgeplot_data, args.o, title, upper_ylim=args.ylim)
        smg.smudgeplot(smudgeplot_data, log=False, palette=args.col_ramp, invert_cols=args.invert_cols)
        smg.smudgeplot(smudgeplot_data, log=True, palette=args.col_ramp, invert_cols=args.invert_cols)

        smg.fin()

    sys.stderr.write("\nLoading data\n")
    coverages = smg.Coverages(smg.load_hetmers(args.infile))
    sys.stderr.write("\nMasking errors using local aggregation algorithm\n")

    if _parser.task == "peak_aggregation":

        coverages.local_aggregation(distance=args.d, noise_filter=args.nf, mask_errors=args.mask_errors)
        coverages.write_peaks()

    if _parser.task == "all":

        coverages.local_aggregation(distance=args.d, noise_filter=1000, mask_errors=True)
        coverages.count_kmers()
        sys.stderr.write(
            f"\t\
            Total kmers: {coverages.total_kmers}\n\t \
            Genomic kmers: {coverages.total_genomic_kmers}\n\t \
            Genomic kmers in smudges: {coverages.total_genomic_kmers_in_smudges}\n\t \
            Sequencing errors: {coverages.total_error_kmers}\n\t \
            Fraction of errors: {round(coverages.total_error_kmers/coverages.total_kmers, 3)}"
        )

        smudge_size_cutoff = (
            0  # 0.01  # this is % of all k-mer pairs smudge needs to have to be considered a valid smudge
        )
        smudges = smg.Smudges(coverages.cov_tab, coverages.total_genomic_kmers)

        if args.cov == 0.0:
            sys.stderr.write("\nInferring 1n coverage using grid algorithm\n")

            smudges.get_centrality_df(args.cov_min, args.cov_max, smudge_size_cutoff)
            np.savetxt(
                args.o + "_centralities.txt",
                np.around(smudges.centrality_df, decimals=6),
                fmt="%.4f",
                delimiter="\t",
            )

            limit = 0.7
            if coverages.error_fraction < limit:
                cov = smudges.cov
            else:
                cov = 0

            sys.stderr.write("\nCreating centrality plot\n")
            smudges.centrality_plot(args.o, args.format)
            sys.stderr.write(f"\nInferred coverage: {cov:.3f}\n")

        else:
            cov = args.cov
            sys.stderr.write(f"\nUser defined coverage: {cov:.3f}\n")

        sys.stderr.write("\nCreating smudge report\n")

        smudges.local_agg_smudge_container = smudges.get_smudge_container(cov, smudge_size_cutoff, "local_aggregation")
        annotated_smudges = list(smudges.local_agg_smudge_container.keys())
        with open(args.o + "_with_annotated_smu.txt", "w") as annotated_smu:
            annotated_smu.write("covB\tcovA\tfreq\tsmudge\n")
            for smudge in annotated_smudges:
                formated_smudge = smg.smudge2short(smudge)
                for idx, covB, covA, freq, smu in smudges.local_agg_smudge_container[smudge].itertuples():
                    annotated_smu.write(f"{covB}\t{covA}\t{freq}\t{formated_smudge}\n")

        smg.generate_smudge_report(smudges, coverages, cov, args, smudge_size_cutoff, print_header=True)
        sys.stderr.write("\nCreating smudgeplots\n")
        smg.generate_plots(
            smudges,
            coverages,
            cov,
            smudge_size_cutoff,
            args.o,
            title,
            fmt=args.format,
            upper_ylim=args.ylim,
            json_report=args.json_report,
            input_params=vars(args),
            palette=args.col_ramp,
            invert_cols=args.invert_cols
        )

    fin()

if __name__ == "__main__":
    main()
