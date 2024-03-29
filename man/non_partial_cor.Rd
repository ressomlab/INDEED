% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/non_partial_cor.R
\name{non_partial_cor}
\alias{non_partial_cor}
\title{Non-partial correlation analysis}
\usage{
non_partial_cor(
  data = NULL,
  class_label = NULL,
  id = NULL,
  method = "pearson",
  p_val = NULL,
  permutation = 1000,
  permutation_thres = 0.05,
  fdr = TRUE
)
}
\arguments{
\item{data}{This is a p*n dataframe that contains the expression levels for all biomolecules and samples.}

\item{class_label}{This is a 1*n dataframe that contains the class label with 0 for group 1 and 1 for group 2.}

\item{id}{This is a p*1 dataframe that contains the ID for each biomolecule.}

\item{method}{This is a character string indicating which correlation method is to use. The 
options are either "pearson" as the default or "spearman".}

\item{p_val}{This is optional. It is a p*1 dataframe that contains the p-value for each biomolecule from DE analysis.}

\item{permutation}{This is a positive integer representing the desired number of permutations. 
The default is 1000.}

\item{permutation_thres}{This is a integer representing the threshold for the permutation test. 
The default is 0.05 to achieve 95 percent confidence.}

\item{fdr}{This is a boolean value indicating whether to apply multiple testing correction (TRUE) 
or not (FALSE). The default is TRUE. However, if users find the output network is too sparse 
even after relaxing the permutation_thres, it's probably a good idea to turn off the multiple testing correction.}
}
\value{
A list containing an activity score dataframe with "ID", "P_value", "Node_Degree" and 
    "Activity_Score" as columns and a differential network dataframe with the binary and the 
    weight connection values.
}
\description{
A method that integrates differential expression (DE) analysis and differential 
    network (DN) analysis to select biomarker candidates for cancer studies. non_partial_cor is 
    a one step function for users to perform the typical correlation analysis. No pre-processing 
    step is required.
}
\examples{
non_partial_cor(data = Met_GU, class_label = Met_Group_GU, id = Met_name_GU,
                          method = "pearson", p_val = pvalue_M_GU, permutation = 1000, 
                          permutation_thres = 0.05, fdr = TRUE)
}
