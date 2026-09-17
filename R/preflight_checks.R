#' Query checks
#'
#' @param query The query data frame
#' @param cols_to_match The columns to match
#' @return A data frame with the preflight checks passed
#' @import dplyr magrittr stringr data.table
#' @export
preflight_query <- function(query, cols_to_match){

    # QUERY CHECKS
    #========================================================
    # make sure query is a data frame or data.table
    stopifnot(any(is.data.frame(query), is.data.table(query)))

    # convert query to data.table
    if(is.data.frame(query)){
        query <- as.data.table(query, keep.rownames = TRUE)}

    # change VDJ column names
    for(col in cols_to_match){
        setnames(query, col, names(cols_to_match)[cols_to_match == col])
    }

    # remove unnecessary columns - keep only barcodes and columns to match
    # Check if rn column exists, if not create it from row names
    if(!"rn" %in% names(query)){
        query[, rn := as.character(1:.N)]
    }
    cols_to_keep <- c("rn", names(cols_to_match))
    query <- query[, cols_to_keep, with = FALSE]
    setnames(query, "rn", "barcodes")

    nrow_query <- nrow(query)

    # get CDR3 and gene columns to match
    cdr3_to_match <- cols_to_match[str_detect(names(cols_to_match), "CDR[HL]3")]
    genes_to_match <- cols_to_match[!str_detect(names(cols_to_match), "CDR[HL]3")]

    # FOR CDR3 COLUMNS
    #========================================================
    for(i in names(cdr3_to_match)){

        # remove BCR with less than 5 CDR3 AA
        rows_to_remove <- !str_detect(query[[i]], "^(?=(?:[^A-Z]*[A-Z]){5})[A-Z]+$")
        if(any(rows_to_remove)){
            query <- query[!rows_to_remove]
        }
        
        # remove NAs and empty strings
        query <- query[
            !is.na(get(i)) & 
            get(i) != "" & 
            get(i) != "NA" & 
            get(i) != "None"
        ]
        
        # add length column
        query[, paste0(i, "_length") := nchar(get(i))]
    }
    
    # report
    nrow_query_cdr3 <- nrow(query)
    message(paste0("\nQUERY : ", nrow_query_cdr3, "/", nrow_query, " rows remaining after removing rows with less than 5 CDR3 AA"))
    
    # FOR GENE COLUMNS
    #========================================================
    # remove NAs from genes_to_match
    if(length(genes_to_match) > 0){
        for(i in names(genes_to_match)){
            query <- query[
                str_detect(get(i), "^IG[HKL]") & 
                get(i) != "" & 
                get(i) != "NA" & 
                get(i) != "None" & 
                !is.na(get(i))
            ]
        }
    }

    nrow_query_genes <- nrow(query)
    message(paste0("\nQUERY : ", nrow_query_genes, "/", nrow_query_cdr3, " rows remaining after removing rows with invalid VJ genes"))

    # return query
    return(query)
}




#' Reference checks
#'
#' @param reference The reference data frame
#' @param cols_to_match The columns to match
#' @return A data frame with the preflight check passed
#' @export
preflight_reference <- function(reference, cols_to_match){
    
    # check reference is a data frame
    stopifnot(any(is.data.frame(reference), is.data.table(reference)))

    # convert reference to data.table if it is a data frame
    if(is.data.frame(reference)){
        reference <- as.data.table(reference, keep.rownames = TRUE)
        # Check if rn column exists, if not create it from row names
        if(!"rn" %in% names(reference)){
            reference[, rn := as.character(1:.N)]
        }
    }

    # get number of rows in reference
    nrow_reference <- nrow(reference)

    # get CDR3 and gene columns to match
    cdr3_to_match <- cols_to_match[str_detect(names(cols_to_match), "CDR[HL]3")]
    names(cdr3_to_match) <- paste0("ref_", names(cdr3_to_match))

    genes_to_match <- cols_to_match[!str_detect(names(cols_to_match), "CDR[HL]3")]
    if(length(genes_to_match) > 0){
        names(genes_to_match) <- paste0("ref_", names(genes_to_match))}

    # FOR CDR3 COLUMNS
    #========================================================
    for(i in names(cdr3_to_match)){

        # remove BCR with less than 5 CDR3 AA
        rows_to_remove <- !str_detect(reference[[i]], "^(?=(?:[^A-Z]*[A-Z]){5})[A-Z]+$")
        if(any(rows_to_remove)){
            reference <- reference[!rows_to_remove]
        }
        
        # remove NAs and empty strings
        reference <- reference[
            !is.na(get(i)) & 
            get(i) != "" & 
            get(i) != "NA" & 
            get(i) != "None"
        ]
        
        # add length column
        reference[, paste0(i, "_length") := nchar(get(i))]
    }
    
    # report
    nrow_reference_cdr3 <- nrow(reference)
    message(paste0("\nREFERENCE : ", nrow_reference_cdr3, "/", nrow_reference, " rows remaining after removing rows with less than 5 CDR3 AA"))
    
    # FOR GENE COLUMNS
    #========================================================
    # remove NAs from genes_to_match
    if(length(genes_to_match) > 0){
        for(i in names(genes_to_match)){
            reference <- reference[
                str_detect(get(i), "^IG[HKL]") & 
                get(i) != "" & 
                get(i) != "NA" & 
                get(i) != "None" & 
                !is.na(get(i))
            ]
        }
    }

    nrow_reference_genes <- nrow(reference)
    message(paste0("\nREFERENCE : ", nrow_reference_genes, "/", nrow_reference_cdr3, " rows remaining after removing rows with invalid VJ genes"))

    # return reference
    return(reference)
}

#' Preflight checks
#'
#' @param ncores The number of cores to use
#' @param dist_method The distance method to use
#' @param output_dir The output directory
preflight_checks <- function(ncores, dist_method, output_dir){
    
    # Package version message
    pkg_version <- as.character(packageVersion("ClonoMappeR"))
    message(paste0("Runnning on ClonoMappeR v", pkg_version, "..."))

    # check ncores is valid
    stopifnot(ncores > 0 & is.numeric(ncores))

    # check dist_method is valid
    stopifnot(all(dist_method %in% c("levenshtein", "hamming")))

    # check output_dir is valid
    if(!is.null(output_dir)){
        stopifnot(dir.exists(output_dir))}

}