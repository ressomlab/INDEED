## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval = F----------------------------------------------------------------
#  # install.packages("devtools")
#  devtools::install_github("ressomlab/INDEED")

## ---- eval = F----------------------------------------------------------------
#  # load INDEED
#  library(INDEED)

## ----dataset, eval = F--------------------------------------------------------
#  # Data matrix contains the expression levels of 39 metabolites from 120 subjects
#  # (6 metabolites and 10 subjects are shown)
#  head(Met_GU[, 1:10])
#  # Group label for each subject (40 subjects are shown)
#  Met_Group_GU[1:40]
#  # Metabolite KEGG IDs (10 metabolites are shown)
#  Met_name_GU[1:10]

## ---- eval = F----------------------------------------------------------------
#  result <- non_partial_cor(data = Met_GU, class_label = Met_Group_GU, id = Met_name_GU, method = "pearson",
#                            p_val = pvalue_M_GU, permutation = 1000, permutation_thres = 0.05, fdr = TRUE)

## ---- eval = F----------------------------------------------------------------
#  pre_data <- select_rho_partial(data = Met_GU, class_label = Met_Group_GU, id = Met_name_GU,
#                                 error_curve = TRUE)
#  result <- partial_cor(data_list = pre_data, rho_group1 = 'min', rho_group2 = "min", p_val = pvalue_M_GU,
#                        permutation = 1000, permutation_thres = 0.05, fdr = TRUE)

## ---- eval = F----------------------------------------------------------------
#  result <- non_partial_cor(data = Met_GU, class_label = Met_Group_GU, id = Met_name_GU, method = "pearson",
#                            p_val = pvalue_M_GU, permutation = 1000, permutation_thres = 0.05, fdr = FALSE)
#  network_display(result = result, nodesize= 'Node_Degree', nodecolor= 'Activity_Score',
#                  edgewidth= FALSE, layout= 'nice')

