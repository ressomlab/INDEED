#' GU cirrhosis (CIRR) and GU Hepatocellular carcinoma (HCC) data.
#'
#' A dataset containing the expression levels for each of the 120 subjects (HCC: 60; CIRR: 60)
#' in terms of 39 metabolites.
#'
#' @format A data frame with 39 variables (rows) and 120 subjects (columns).
"Met_GU"


#' Group label.
#'
#' A dataset containing group information (CIRR group: 0 and HCC group: 1).
#'
#' @format A data frame with 1 row and 120 (subjects) columns.
"Met_Group_GU"


#' KEGG ID
#'
#' A dataset containing the KEGG ID for each metabolite.
#'
#' @format A data frame with 39 KEGG ID as rows and 1 column:
"Met_name_GU"



#' P-values obtained by differential expression (DE) analysis.
#'
#' A dataset containing the p-value for each metabolite obtained through DE analysis.
#'
#' @format A data frame with 39 rows and 2 variables:
#' \describe{
#'   \item{KEGG.ID}{KEGG ID}
#'   \item{p.value}{p-value}
#' }
"pvalue_M_GU"


