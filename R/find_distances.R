#' Find hamming distance
#'
#' @param query The query data frame
#' @param reference The reference data frame
#' @param cols_to_match The columns to match
#' @param ncores The number of cores to use
#' @return A data frame with the hamming distance
#' @import future furrr stringdist alakazam dplyr stringr rlang magrittr data.table
#' @export
find_hamming_dist <- function(query, reference, cols_to_match, ncores){

    message(paste0("\nRunning hamming distance matching..."))

    # Setup future backend for parallelization
    future::plan(future::multisession, workers = ncores)

    # get CDR3 and gene columns to match
    cdr3_to_match <- cols_to_match[str_detect(names(cols_to_match), "CDR[HL]3")]
    genes_to_match <- cols_to_match[!str_detect(names(cols_to_match), "CDR[HL]3")]

    # Pre-filter query by gene matching if needed
    if(length(genes_to_match) > 0){
        query[, gene_key := do.call(paste, c(.SD, sep = "|")), .SDcols = names(genes_to_match)]
        reference[, gene_key := do.call(paste, c(.SD, sep = "|")), .SDcols = paste0("ref_", names(genes_to_match))]
        
        # Get unique gene combinations that exist in both query and reference
        valid_gene_keys <- intersect(query$gene_key, reference$gene_key)
        
        if(length(valid_gene_keys) == 0){
            message("No matching gene combinations found between query and reference")
            return(data.table())
        }
        
        # Filter both datasets to only include valid gene combinations
        query <- query[gene_key %in% valid_gene_keys]
        reference <- reference[gene_key %in% valid_gene_keys]
    }

    # Using future + furrr
    output <- furrr::future_map_dfr(1:nrow(reference), function(i) {
        
        # add CDR3 sequence and length columns from reference data frame in query data frame
        tmp <- copy(query)
        for(col in colnames(reference)){
            tmp[[col]] <- reference[i, ][[col]]
        }

        # filter query for matching CDR3 sequence in heavy/light chains
        if(length(genes_to_match) > 0){
            for(col in names(genes_to_match)){
                tmp <- tmp[get(col) == get(paste0("ref_", col))]
            }
        }
        
        # hamming distance requires matching CDR3 length
        for(col in names(cdr3_to_match)){

            # filter query for matching CDR3 length in heavy/light chains
            tmp <- tmp[get(paste0(col, "_length")) == get(paste0("ref_", col, "_length"))]

            # check if there are any matching sequences
            if(nrow(tmp) == 0) {next}

            # calculate hamming distance
            y <- reference[i, ][[paste0("ref_", col)]]
            tmp[, paste0(col, "_dist_method") := "hamming"]

            # query must have valid CDR3 sequence
            cond1 <- !str_detect(y, "^(?=(?:[^A-Z]*[A-Z]){5})[A-Z]+$")
            cond2 <- any(is.na(tmp[[col]]))

            if(any(cond1, cond2)){
                tmp[, paste0(col, "_dist") := NA]
            } else {
                # Use stringdist for vectorized hamming distance
                tmp[, paste0(col, "_dist") := as.numeric(stringdist::stringdist(y, get(col), method = 'hamming') / nchar(y))]
            }
        }

        return(tmp)
    }, .options = furrr::furrr_options(seed = TRUE))

    return(output)
}

#' Find levenshtein distance
#'
#' @param query The query data frame
#' @param reference The reference data frame
#' @param cols_to_match The columns to match
#' @param ncores The number of cores to use
#' @return A data frame with the levenshtein distance
#' @import future furrr stringdist alakazam dplyr stringr rlang magrittr
#' @export
find_levenshtein_dist <- function(query, reference, cols_to_match, ncores){

    message(paste0("\nRunning levenshtein distance matching..."))

    # Setup future backend for parallelization
    future::plan(future::multisession, workers = ncores)

    # get CDR3 and gene columns to match
    cdr3_to_match <- cols_to_match[str_detect(names(cols_to_match), "CDR[HL]3")]
    genes_to_match <- cols_to_match[!str_detect(names(cols_to_match), "CDR[HL]3")]


    # Pre-filter query by gene matching if needed
    if(length(genes_to_match) > 0){
        query[, gene_key := do.call(paste, c(.SD, sep = "|")), .SDcols = names(genes_to_match)]
        reference[, gene_key := do.call(paste, c(.SD, sep = "|")), .SDcols = paste0("ref_", names(genes_to_match))]
        
        # Get unique gene combinations that exist in both query and reference
        valid_gene_keys <- intersect(query$gene_key, reference$gene_key)
        
        if(length(valid_gene_keys) == 0){
            message("No matching gene combinations found between query and reference")
            return(data.table())
        }
        
        # Filter both datasets to only include valid gene combinations
        query <- query[gene_key %in% valid_gene_keys]
        reference <- reference[gene_key %in% valid_gene_keys]
    }

    # Using future + furrr
    output <- furrr::future_map_dfr(1:nrow(reference), function(i) {
        
        # add CDR3 sequence and length columns from reference data frame in query data frame
        tmp <- copy(query)
        for(col in colnames(reference)){
            tmp[[col]] <- reference[i, ][[col]]
        }
            
        # filter query for matching CDR3 sequence in heavy/light chains
        if(length(genes_to_match) > 0){
            for(col in names(genes_to_match)){
                tmp <- tmp[get(col) == get(paste0("ref_", col))]
            }
        }

        # levenshtein distance does not require matching CDR3 length
        for(col in names(cdr3_to_match)){

            # check if there are any matching sequences
            if(nrow(tmp) == 0) {next}

            # calculate levenshtein distance
            y <- reference[i, ][[paste0("ref_", col)]]
            tmp[, paste0(col, "_dist_method") := "levenshtein"]

            # query must have valid CDR3 sequence
            cond1 <- !str_detect(y, "^(?=(?:[^A-Z]*[A-Z]){5})[A-Z]+$")
            cond2 <- any(is.na(tmp[[col]]))
                
            if(any(cond1, cond2)){
                tmp[, paste0(col, "_dist") := NA]
            } else {
                tmp[, paste0(col, "_dist") := as.numeric(lapply(get(col), function(x){
                    stringdist::stringdist(y, x, method = 'lv') / max(nchar(y), nchar(x))
                }))]
            }
        }

        return(tmp)
    }, .options = furrr::furrr_options(seed = TRUE))

    return(output)
}

#' Find minimum distances
#'
#' @param output The output data frame
#' @return A data frame with the minimum 
#' @import dplyr stringr magrittr
#' @export
find_min_distances <- function(output){

    message(paste0("Finding minimum CDR3 distance for each barcode..."))

    # no candidate pairs survived gene matching, so there is nothing to rank
    if(nrow(output) == 0 | !"ref_CDRH3" %in% colnames(output)){
        warning("\nNo public BCR sequences found")
        return(data.table())}

    # remove rows with NA in ref_CDRH3
    output <- output[!is.na(ref_CDRH3)]

    # filter for rows where all dist_method values are the same
    dist_cols <- colnames(output)[str_detect(colnames(output), "CDR[HL]3_dist$")]
    
    # check if there are any matching sequences
    if(nrow(output) > 0){

        # Calculate mean distance for each row
        output[, mean_dist := rowMeans(.SD, na.rm = TRUE), .SDcols = dist_cols]
        
        # Find minimum distance for each barcode
        output <- output[output[, .I[which.min(mean_dist)], by = barcodes]$V1]

        # found matching BCR sequences
        n_methods <- length(unique(output$CDRH3_dist_method))
        message(paste0("Found ", nrow(output)/n_methods, " public BCR sequences in QUERY..."))

        # return output
        return(output)}

    # if no matching sequences, return NULL
    else{
        warning("\nNo public BCR sequences found")
        return(data.table())}
}


