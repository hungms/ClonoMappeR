#' Detect Public BCR Sequences
#' 
#' @param query The query data frame
#' @param reference The reference data frame
#' @param CDRH3 The heavy chain CDR3 sequence
#' @param VH The heavy chain V gene
#' @param JH The heavy chain J gene
#' @param CDRL3 The light chain CDR3 sequence
#' @param VL The light chain V gene
#' @param JL The light chain J gene
#' @param dist_method The distance method to use
#' @param ncores The number of cores to use
#' @param output_dir The output directory
#' @param output_name The output filename, without the .csv extension
#' @return A data frame with the matched CDR3 sequences. Query barcodes with no
#'   reference match are included as rows holding only a barcode, with NA in
#'   every other column.
#' @import dplyr stringr magrittr data.table
#' @export

find_publicBCR <- function(
    query,
    reference,
    CDRH3,
    VH = NA,
    JH = NA,
    CDRL3 = NA,
    VL = NA,
    JL = NA,
    dist_method = c("levenshtein", "hamming"),
    ncores = 1,
    output_dir = NULL,
    output_name = NULL){

    # Determine columns to match
    #========================================================
    if("binding" %in% colnames(reference)){
        if(all(c(TRUE, FALSE) %in% reference$binding)) {
            message("Database contains non-binding BCRs, removing them for simplicity...")
            reference <- reference %>% filter(binding == TRUE)}
    }
    
    colnames(reference) <- paste0("ref_", colnames(reference))
    cols_to_match <- lookup_cols(query, reference, CDRH3, VH, JH, CDRL3, VL, JL)

    # Record every query barcode before filtering, mirroring how preflight_query assigns them
    #========================================================
    if(is.data.table(query)){
        all_barcodes <- if("rn" %in% colnames(query)){
            as.character(query$rn)} else {
            as.character(seq_len(nrow(query)))}}
    else{
        all_barcodes <- rownames(query)}

    # Run preflight checks for query and reference
    #========================================================
    preflight_checks(ncores, dist_method, output_dir)
    query <- preflight_query(query, cols_to_match)
    reference <- preflight_reference(reference, cols_to_match)

    message(paste0("\nMatching ", nrow(query), " BCR sequences in QUERY against ", nrow(reference), " BCR sequences in REFERENCE..."))

    # Run Gene and CDR3 Matching
    #========================================================
    # create list to store outputs  
    output.list <- list()

    # Find hamming distance
    if("hamming" %in% dist_method){
        hamming_output <- find_hamming_dist(query, reference, cols_to_match, ncores)
        output.list[[length(output.list) + 1]] <- find_min_distances(hamming_output)}

    # Find levenshtein distance
    if("levenshtein" %in% dist_method){
        levenshtein_output <- find_levenshtein_dist(query, reference, cols_to_match, ncores)
        output.list[[length(output.list) + 1]] <- find_min_distances(levenshtein_output)}

    # combine outputs
    heavy_name <- paste0(intersect(c("VH", "JH", "CDRH3"), names(cols_to_match)), collapse = "")
    light_name <- intersect(c("VL", "JL", "CDRL3"), names(cols_to_match))
    if(length(light_name) > 0){
        light_name <- paste0("_", paste0(light_name, collapse = ""))}
    else{
        light_name <- ""}
    match_method <- paste0(heavy_name, light_name)

    output <- bind_rows(output.list) %>%
        mutate(
            match_method = match_method)

    # Add back query barcodes with no reference match, leaving every other column as NA
    #========================================================
    unmatched_barcodes <- setdiff(all_barcodes, unique(output$barcodes))
    if(length(unmatched_barcodes) > 0){
        message(paste0("Adding ", length(unmatched_barcodes), " QUERY barcodes without a reference match..."))
        output <- bind_rows(output, data.table(barcodes = unmatched_barcodes))}

    # write output to file
    if(!is.null(output_dir) & !is.null(output_name)){
        filename <- paste0(output_dir, "/", output_name, ".csv")
        write.csv(output, filename, row.names = F)}
    else{
        message("No output directory or output name provided, skipping output to file")}

    # return output
    return(output)
}
