## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval = F-----------------------------------------------------------
#  # install.packages("devtools")
#  devtools::install_github("ressomlab/INDEED")

## ---- eval = F-----------------------------------------------------------
#  # load INDEED
#  library(INDEED)
#  # Loading required package: glasso

## ---- eval = F-----------------------------------------------------------
#  non_partial_cor(data = Met_GU, class_label = Met_Group_GU, id = Met_name_GU, method = "spearman")

## ---- eval = F-----------------------------------------------------------
#  pre_data <- select_rho_partial(data = Met_GU, class_label = Met_Group_GU, id = Met_name_GU)
#  result <- partial_cor(data_list = pre_data, rho_group1 = 'min', rho_group2 = "min", permutation = 1000, p_val = pvalue_M_GU, permutation_thres = 0.05)

## ---- eval = F-----------------------------------------------------------
#  network_display(result)

## ---- eval = F-----------------------------------------------------------
#  pre_data <- select_rho_partial(data = Met_GU, class_label = Met_Group_GU, id = Met_name_GU)
#  partial_cor(data_list = pre_data, rho_group1 = 'min', rho_group2 = "min", permutation = 500, p_val = pvalue_M_GU, permutation_thres = 0.05)
#  

## ---- eval = F-----------------------------------------------------------
#  result2 <- non_partial_cor(data = Met_GU, class_label = Met_Group_GU, id = Met_name_GU, method = "spearman", permutation_thres = 0.05)

## ---- eval = F-----------------------------------------------------------
#  network_display(results = result2, layout = 'nice', nodesize= 'Node_Degree', nodecolor = 'Activity_Score', edgewidth = 'NO')

