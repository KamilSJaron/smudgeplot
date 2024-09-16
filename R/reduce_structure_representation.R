#' @title reduce_structure_representation
#'
#' @description
#' \code{reduce_structure_representation} is the function that collapses labels loger than 4 characters into 4 character representations
#'
#' @param .peak_sizes - a table with peak sizes, the truly important column in there is just the "structure" one
#'
#' @export

reduce_structure_representation <- function(long_smudge_labels){
    structures_to_adjust <- sapply(long_smudge_labels, nchar) > 4
    if (any(structures_to_adjust)) {
        decomposed_struct <- strsplit(long_smudge_labels[structures_to_adjust], '')
        As <- sapply(decomposed_struct, function(x){ sum(x == 'A') } )
        Bs <- sapply(decomposed_struct, length) - As
        long_smudge_labels[structures_to_adjust] <- paste0(As, 'A', Bs, 'B')
    }
    long_smudge_labels
}
