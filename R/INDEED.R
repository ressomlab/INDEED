#' INDEED: A package for biomarker candidate prioritization.
#'
#' The INDEED package provides a important functions below:
#' select_rho_partial
#'
#' @section select_rho_partial function:
#' select_rho_partial function preprocess data for partical correlation analysis,
#' the result contains list of preprocessed data and rho values and
#' error plot for user to choose desired rho value

#' @section partial_cor function:
#' partial_cor function performs partical correlation analysis
#' user input preprocessed list from select_rho_partial step
#' and the rho choosing method or rho of their choice
#' and number of permutations (default 1000), p-value is optional
#' the result of score table and differential network will be returned

#' @section non_partial_cor function:
#' non_partial_cor function performs correlation analysis
#' user input data,class label,p-value, sample id,
#' number of permutations, and method(default pearson)
#' p value is optional
#' the result of score table and differential network will be returned
#'
#' @docType package
#' @name INDEED
NULL







