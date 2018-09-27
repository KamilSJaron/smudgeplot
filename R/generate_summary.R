#' @title generate_summary
#'
#' @description
#' \code{generate_summary}
#'
#' @export

generate_summary <- function(.args, .smudge_summary){
    smudge_summary(.args$output, "1n coverage estimates (Coverage of every haplotype; Don't confuse with genome coverage which is (ploidy * 1n coverage).)")
    smudge_summary(.args$output, "* User defined 1n coverage:", .args$n_cov)
    smudge_summary(.args$output, "* Subset 1n coverage estimate:", .smudge_summary$n_subset_est)
    smudge_summary(.args$output, "* Highest peak 1n coverage estimate:", round(.smudge_summary$n_peak_est, 1))
    smudge_summary(.args$output, "1n coverage used in smudgeplot (one of the three above):", round(.smudge_summary$n, 1))

    smudge_summary(.args$output, "* Estimated ploidy:", .smudge_summary$genome_ploidy)

    for_sure_heterozygous <- sapply(.smudge_summary$peak_sizes$structure, function(x){sum(unlist(strsplit(x, split = '')) == 'B') == 1})
    # THIS LINE
    minimal_number_of_heterozygous_loci <- ceiling(sum(.smudge_summary$peak_sizes$abs_size[for_sure_heterozygous]) / .args$kmer_size)

    smudge_summary(.args$output, "* Minimal number of heterozygous loci:", minimal_number_of_heterozygous_loci)
    smudge_summary(.args$output, "Note: This number is NOT an estimate of the total number heterozygous loci, it's merly setting the lower boundary if the inference of heterozygosity peaks is correct.")

    smudge_summary(.args$output, "* Proportion of heterozygosity carried by pairs in different genome copies (table)")
    tab_to_print <- data.frame( genome_copies = 2:max(.smudge_summary$peak_sizes$ploidy),
                                propotion_of_heterozygosity = round(sapply(2:max(.smudge_summary$peak_sizes$ploidy), function(x){ sum(.smudge_summary$peak_sizes$rel_size[.smudge_summary$peak_sizes$ploidy == x]) }), 2))
    smudge_summary(.args$output, paste0(capture.output(tab_to_print), collapse = "\n"))

    carried_by_paralogs <- sum(.smudge_summary$peak_sizes$rel_size[.smudge_summary$peak_sizes$ploidy > .smudge_summary$genome_ploidy])
    smudge_summary(.args$output, "* Proportion of heterozygosity carried by paralogs:", round(carried_by_paralogs, 3))

    smudge_summary(.args$output, "* Summary of all detected peaks (table)")
    tab_to_print <- .smudge_summary$peak_sizes[,c(11,2,3,8,9)]
    tab_to_print[,c(2:5)] <- round(tab_to_print[,c(2:5)], 2)
    colnames(tab_to_print) <- c('peak','kmers [#]', 'kmers [proportion]', 'summit B / (A + B)', 'summit A + B')
    smudge_summary(.args$output, paste0(capture.output(tab_to_print), collapse = "\n"))

    write.table(tab_to_print, paste0(.args$output, '_summary_table.tsv'),
                quote = F, sep = '\t', row.names = F)
}