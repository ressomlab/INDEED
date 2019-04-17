#' INDEED: A package for biomarker candidate prioritization.
#'
#' The INDEED package provides important functions below:
#' select_rho_partial, partial_cor, non_partial_cor and network_display.
#'
#' @section select_rho_partial function:
#' select_rho_partial function preprocess data for partical correlation analysis,
#' the result contains list of preprocessed data and rho values and
#' error plot for user to choose desired rho value

#' @section partial_cor function:
#' partial_cor function performs partical correlation analysis
#' based on user input preprocessed list from select_rho_partial step
#' and the rho choosing method or rho of their choice
#' and number of permutations (default 1000), p-value is optional,
#' the result of score table and differential network will be returned

#' @section non_partial_cor function:
#' non_partial_cor function performs correlation analysis
#' based on user input data, class label, p-value, sample id,
#' number of permutations, and method (default pearson)
#' p value is optional,
#' the result of score table and differential network will be returned

#' @section network_display function:
#' A tool to assist in the network visualization of the results from INDEED functions patial_cor()
#' and non_partial_cor().

#'
#' @docType package
#' @name INDEED
NULL







