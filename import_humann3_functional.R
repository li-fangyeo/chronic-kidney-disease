#' Import humann3 with relative abundance results o \code{TreeSummarizedExperiment}
#' 
#' The function import human3 results. As, phenotype and endpoints data from FR02 
#' and/or FR07 are loaded and combined. Resulting in an \code{TreeSummarizedExperiment}
#' object containing all.
#' 
#' @param pheno_FR07 a single \code{character} value defining the file
#'   path to the phenotype data corresponding to FR07.
#'   
#' @param pheno_FR02 a single \code{character} value defining the file
#'   path to the phenotype data corresponding to FR02.
#'   
#' @param endpoints a single \code{character} value defining the file
#'   path to the endpoints data corresponding to FR02 and FR07.
#'   
#' @param abundance_table a single \code{character} value defining the file
#'   path to the abundance table.
#'   
#' @param read_counts a single \code{character} value defining the file path to
#'   the samples read counts. 
#'   
#' @return  A \code{TreeSummarizedExperiment} object.

#' @examples
#' \dontrun{
#' # Use the link_paths function to find the match the correct paths linked to
#' # the chosen importer
#' paths <- link_paths("import_humann3_functional", data_paths)
#' # Importing and constructing the TreeSummarizedExperiment object
#' tse <- import_humann3_functional(
#'     pheno_FR07=paths[["pheno_FR07"]],
#'     pheno_FR02=paths[["pheno_FR02"]],
#'     endpoints=paths[["endpoints"]],
#'     abundance_table=paths[["abundance_table"]],
#'     read_counts=paths[["read_counts"]])
#' 
#' tse
#' 
#' }

#' @name import_humann3_functional
#' @export
#' @importFrom stringr str_split
#' @importFrom S4Vectors DataFrame metadata metadata<-
#' @importFrom utils read.table
#' @importFrom SummarizedExperiment rowData colData<- assay assay<-
#' @importFrom SingleCellExperiment altExp<-
#' @importFrom TreeSummarizedExperiment rownames<- rowTree<-
#' @importFrom plyr rbind.fill
#' @importFrom dplyr %>% 
#' @importFrom rlang .data

import_humann3_functional<- function(pheno_FR07=NULL,
                                   pheno_FR02=NULL,
                                   endpoints,
                                   abundance_table,
                                   read_counts) {

    # input checks
    if(!.is_non_empty_string(abundance_table)){
        stop("'abundance_table' must be a single character value, as a path to
        the abundance table data",
             call. = FALSE)
    }
    
    ## Phenotype and endpoints Data
    col_data <- .load_colData(pheno_FR07, pheno_FR02, endpoints, read_counts)
    
    ## Loading Functional Data
    assay <- read.table(abundance_table,
                        sep="\t", header = TRUE, row.names=1,
                        check.names =FALSE, comment.char = "")
	colnames(assay) <- gsub("\\.R1\\.trimmed\\.filtered\\_Abundance", "", colnames(assay))
    common_ids <- intersect(colnames(assay), rownames(col_data))
    tse <- TreeSummarizedExperiment(
        assays=SimpleList(rel_ab=as.matrix(assay[,common_ids])/100),
        colData=DataFrame(col_data[common_ids,]))

    # Adding the TreeSummarizedExperiment object description as metadata
    metadata(tse)$Description <- paste0(
        "The object was constructed: ", Sys.time(),"\n",
        "Penotype data was loaded from: \n",
        "\tFR02 pheno: ", ifelse(is.null(pheno_FR02),"",.clean_path(pheno_FR02)),"\n",
        "\tFR07 pheno: ", ifelse(is.null(pheno_FR07),"",.clean_path(pheno_FR07)),"\n",
        "\tFR endpoints: ", .clean_path(endpoints),"\n",
        "Functional data was loaded from: \n",
        "\t abundance table: ", .clean_path(abundance_table), "\n")
    return(tse)
}
