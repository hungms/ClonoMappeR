#' Lookup columns to match
#'
#' @param query The query data frame
#' @param reference The reference data frame
#' @param CDRH3 The heavy chain CDR3 sequence
#' @param VH The heavy chain V gene
#' @param JH The heavy chain J gene
#' @param CDRL3 The light chain CDR3 sequence
#' @param VL The light chain V gene
#' @param JL The light chain J gene
#' @return A list of columns to match
#' @import dplyr magrittr stringr data.table
#' @export
lookup_cols <- function(
    query,
    reference,
    CDRH3,
    VH = NA,
    JH = NA,
    CDRL3 = NA,
    VL = NA,
    JL = NA){

    # store gene columns to match
    genes_to_match <- c(VH, JH, VL, JL)
    names(genes_to_match) <- c("VH", "JH", "VL", "JL")
    genes_to_match <- genes_to_match[!is.na(genes_to_match)]

    # store CDR3 columns to match
    CDR3_to_match <- c(CDRH3, CDRL3)
    names(CDR3_to_match) <- c("CDRH3", "CDRL3")
    CDR3_to_match <- CDR3_to_match[!is.na(CDR3_to_match)]

    # check if all columns to match are in query and reference
    cols_to_match <- c(genes_to_match, CDR3_to_match)
    stopifnot(all(cols_to_match %in% colnames(query)))
    stopifnot(all(paste0("ref_", names(cols_to_match)) %in% colnames(reference)))

    return(cols_to_match)
    }