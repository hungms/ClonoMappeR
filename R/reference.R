#' Get reference data frame
#' 
#' @param context The context to get the reference data frame for
#' @param epitope The epitope to get the reference data frame for
#' @param org The organism to get the reference data frame for
#' @param publication The publication to get the reference data frame for
#' @return A data frame with the reference data
#' @import dplyr magrittr stringr data.table
#' @export

get_reference <- function(context, epitope = NULL, org = NULL, publication = NULL){

    # check context is valid
    stopifnot(context %in% c("SarsCoV2", "Tetanus", "Vaccinia", "Measles", "Mumps", "Hpylori", "NP"))

    # read reference data
    reference <- read.csv(system.file("extdata", paste0(context, ".csv"), package = "ClonoMappeR"), header = T, sep = ",")

    # filter by publication
    if(!is.null(publication)){
        existing_publications <- publication %in% unique(reference$publication)
        if(!all(existing_publications)){
            stop("Publication ", paste(publication[!existing_publications], collapse = ", "), " not found in reference data, existing publications:\n", paste(unique(reference$publication), collapse = ", "))
        } else {
            reference <- reference[reference$publication %in% publication, , drop = FALSE]
        }
    }

    # filter by epitope
    if(!is.null(epitope)){
        existing_epitopes <- epitope %in% unique(reference$epitope)
        if(!all(existing_epitopes)){
            stop("Epitope ", paste(epitope[!existing_epitopes], collapse = ", "), " not found in reference data, existing epitopes:\n", paste(unique(reference$epitope), collapse = ", "))
        } else {
            reference <- reference[reference$epitope %in% epitope, , drop = FALSE]
        }
    }

    # filter by org
    if(!is.null(org)){
        existing_orgs <- org %in% unique(reference$org)
        if(!all(existing_orgs)){
            stop("Org ", paste(org[!existing_orgs], collapse = ", "), " not found in reference data, existing orgs:\n", paste(unique(reference$org), collapse = ", "))
        } else {
            reference <- reference[reference$org %in% org, , drop = FALSE]
        }
    }

    # return reference data
    return(reference)
}
