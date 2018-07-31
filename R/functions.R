#' @title Data preprocessing for partial correlaton analysis
#' @description A method that integrates differential expression (DE) analysis
#'   and differential network (DN) analysis to select biomarker candidates for
#'   survival time prediction. pre_partial is the pre-processing step for INDEED
#'   partial differential analysis
#' @param data input matrix of expression from all metabolites from all samples
#' @param class_label a binary array with 0: group 1; 1: group 2.
#' @param id an array of biomolecule ID to label.
#' @examples pre_partial(data=Met_GU,class_label = Met_Group_GU,id=Met_name_GU)
#' @return a list of processed data for next step and rho, error curve for group 1 and 2
#' @import devtools
#' @importFrom glasso glasso
#' @importFrom stats qnorm cor quantile var sd glm
#' @importFrom graphics abline title plot lines
#' @export
pre_partial <- function(data = NULL, class_label = NULL, id = NULL) {
    data_bind <- rbind(data , class_label)
    raw_group_1 <- data_bind[,data_bind[nrow(data_bind),] == 0][1:(nrow(data_bind) - 1),]  # Group 1: p*n1
    raw_group_2 <- data_bind[,data_bind[nrow(data_bind),] == 1][1:(nrow(data_bind) - 1),]  # Group 2: p*n2

    p <- nrow(raw_group_2)
    n_group_2 <- ncol(raw_group_2)
    n_group_1 <- ncol(raw_group_1)

    # Z-transform the data for group-specific normalization
    data_group_1 <- scale(t(raw_group_1)) # Group 1: n1*p
    data_group_2 <- scale(t(raw_group_2)) # Group 2: n2*p
    cov_group_1 <- var(data_group_1)
    cov_group_2 <- var(data_group_2)
    ## Apply grapihcal LASSO given data_group_1 and data_group_2
    #  Group 1 first
    n_fold <- 5 # number of folds
    rho <- exp(seq(log(0.6), log(0.01), length.out = 20))

    #group1
    # draw error curve
    error_group1 <- choose_rho(data_group_1, n_fold, rho)
    par(mfrow=c(1,2))
    plot(rho, error_group1$log.cv, xlab = expression(lambda), ylab = "Error")
	  lines(rho, error_group1$log.cv)
    title(main=paste("Group1 error curve", "using cross validation", sep="\n"))

    # chosse optimal rho
    rho[error_group1$log.cv == min(error_group1$log.cv)] # rho based on minimum rule
    abline(v = rho[error_group1$log.cv == min(error_group1$log.cv)], col = "red", lty = 3)

    # one standard error rule
    abline(h = min(error_group1$log.cv) + error_group1$log.rho[error_group1$log.cv == min(error_group1$log.cv)], col = "blue")

    # rhos that are under blue line
    rho_under_blue <- rho[error_group1$log.cv < min(error_group1$log.cv) + error_group1$log.rho[error_group1$log.cv == min(error_group1$log.cv)]]

    # rhos that are on the right of the red line
    rho_right_red <- rho[rho > rho[error_group1$log.cv == min(error_group1$log.cv)]]

    #rhos that are under blue line and to the right of the red line
    rho_one_std_group1 <- max(intersect(rho_under_blue, rho_right_red ))

    rho_min_rule_group1 <- rho[error_group1$log.cv == min(error_group1$log.cv)]

    #group2
    # draw error curve
    error_group2 <- choose_rho(data_group_2, n_fold, rho)
    plot(rho, error_group2$log.cv, xlab = expression(lambda), ylab = "Error")
	  lines(rho, error_group2$log.cv)
    title(main=paste("Group2 error curve", "using cross validation", sep="\n"))

    # chosse optimal rho
    rho[error_group2$log.cv == min(error_group2$log.cv)] # rho based on minimum rule
    abline(v = rho[error_group2$log.cv == min(error_group2$log.cv)], col = "red", lty = 3)

    # one standard error rule
    abline(h = min(error_group2$log.cv) + error_group2$log.rho[error_group2$log.cv == min(error_group2$log.cv)], col = "blue")

    # rhos that are under blue line
    rho_under_blue <- rho[error_group2$log.cv < min(error_group2$log.cv) + error_group2$log.rho[error_group2$log.cv == min(error_group2$log.cv)]]

    # rhos that are on the right of the red line
    rho_right_red <- rho[rho > rho[error_group2$log.cv == min(error_group2$log.cv)]]

    #rhos that are under blue line and to the right of the red line
    rho_one_std_group2 <- max(intersect(rho_under_blue, rho_right_red )) # export as global

    rho_min_rule_group2 <- rho[error_group2$log.cv == min(error_group2$log.cv)]   # export as global

    rho_df <- data.frame(c(rho_one_std_group1,rho_min_rule_group1), c(rho_one_std_group2,rho_min_rule_group2), row.names=c("one standard error","minimum"))
    colnames(rho_df)<-c('group1',"group2")
    #cov_group_1
    data_list <- list(p=p,rho_table =rho_df,cov_group_1=cov_group_1,cov_group_2=cov_group_2,n_group_1=n_group_1, n_group_2=n_group_2, data_group_1=data_group_1, data_group_2=data_group_2,class_label=class_label,id=id,data=data)
    return(data_list)
    }


#' @title Data preprocessing for partial correlaton analysis
#' @description A method that integrates differential expression (DE) analysis
#'   and differential network (DN) analysis to select biomarker candidates for
#'   survival time prediction. partial_cor is the second step of partial correlation
#'   calculation after the output result from pre_partial function
#' @param data_list, list of pre-processed data from pre_partial function
#' @param rho_group1 rule to choose rho for group 1, "min": minimum rho,
#' "ste" one standard error from minimum, or user can input rho of their choice, default: minimum
#' @param rho_group2 rule to choose rho for group 1, "min": minimum rho,
#' "ste" one standard error from minimum, or user can input rho of their choice, default: minimum
#' @param p_val optional, a dataframe contains p values for each metabolite/molecule
#' @param permutation, a positive integer of desired number of permutations, default 1000
#' @examples preprocess<- pre_partial(data=Met_GU,class_label = Met_Group_GU,id=Met_name_GU)
#'    partial_cor(data_list=preprocess,rho_group1='min',
#'    rho_group2="min",permutation = 1000,p_val=pvalue_M_GU)
#' @return a list containing a score dataframe and a differential network dataframe
#' @import devtools
#' @importFrom glasso glasso
#' @importFrom stats qnorm cor quantile var sd glm
#' @importFrom graphics abline title plot lines par
#' @export
partial_cor <- function(data_list =NULL, rho_group1=NULL,rho_group2=NULL, permutation=1000, p_val = NULL){
  if(missing(data_list)) {stop("please provide data_list from pre_partial function")}
  else{
    #group1
    if (rho_group1 =='min'){ rho_group_1_opt =data_list$rho_table[1,2] }
    else if (rho_group1 =='ste'){ rho_group_1_opt =data_list$rho_table[1,1] }
    else if (is.numeric(rho_group1) & rho_group1>0) {rho_group_1_opt = rho_group1}
    else if (is.numeric(rho_group1) & rho_group1<=0) {stop("please provide data_list from pre_partial function")}
    else {rho_group_1_opt =data_list$rho_table[1,2]} #default is minimum rho if no rule specified and no valid input entered
    #group2
    if (rho_group2 =='min'){ rho_group_2_opt =data_list$rho_table[2,2] }
    else if (rho_group2 =='ste'){ rho_group_2_opt =data_list$rho_table[2,1] }
    else if (is.numeric(rho_group2) & rho_group2>0) {rho_group_2_opt = rho_group2}
    else if (is.numeric(rho_group2) & rho_group2<=0) {stop("please provide data_list from pre_partial function")}
    else {rho_group_2_opt =data_list$rho_table[2,2]} #default is minimum rho if no rule specified and no valid input entered


    pre_group_1 <- glasso(data_list$cov_group_1, rho = rho_group_1_opt)
    thres <- 1e-3
    sum(abs(pre_group_1$wi) > thres)
    pre_group_1$wi[1:10, 1:10]
    # compute partial correlation
    pc_group_1 <- compute_par(pre_group_1$wi)
    # examine the partial correlation matrix
    sum(abs(pc_group_1) > thres)
    pc_group_1[1:10, 1:10]
    pre_group_2 <- glasso(data_list$cov_group_2, rho = rho_group_2_opt)
    # examine the precision matrix
    sum(abs(pre_group_2$wi) > thres)
    pre_group_2$wi[1:10,1:10]
    # compute partial correlation
    pc_group_2 <- compute_par(pre_group_2$wi)
    # examine the partial correlation matrix
    sum(abs(pc_group_2) > thres)
    pc_group_2[1:10,1:10]
    diff <- pc_group_2 - pc_group_1  # from group 1 to group 2
    thres = 1e-3
    sum(abs(diff) > thres)
    diff[1:10, 1:10]

    ## Permutation test using partial correlation
    if(permutation<=0) {stop("please provide a valid number of permutation (positive integer)")}
    else{
      m <- as.numeric(permutation)
      diff_p <- permutation_pc(m, data_list$p, data_list$n_group_1, data_list$n_group_2, data_list$data_group_1, data_list$data_group_2, rho_group_1_opt, rho_group_2_opt)
      p <- data_list$p
    }
    #####
    ##### final calculation
    #####
    thres_left <- 0.025
    thres_right <- 0.975
    significant_thres <- permutation_thres(thres_left, thres_right, p, diff_p)
    rm(thres_left, thres_right)

    # get binary matrix
    significant_thres_p <- significant_thres$positive
    significant_thres_n <- significant_thres$negative
    binary_link <- matrix(0, p, p) # binary connection
    binary_link[diff < significant_thres_n] <- -1
    binary_link[diff > significant_thres_p] <- 1
    weight_link <- matrix(0, p, p) # weight connection
    weight_link[diff < significant_thres_n] <- diff[diff < significant_thres_n]
    weight_link[diff > significant_thres_p] <- diff[diff > significant_thres_p]
    sum(diff < significant_thres_n)
    sum(diff > significant_thres_p)
    binary_link[1:10, 1:10]
    weight_link[1:10, 1:10]
    rowSums(abs(binary_link)) # node degree for differential networks
    rm(diff_p)

    # Convert adjacent matrix into edge list
    i <- rep(seq_len(nrow(binary_link) - 1), times = (nrow(binary_link)-1):1)
    k <- unlist(lapply(2:nrow(binary_link), seq, nrow(binary_link)))
    binary_link_value <- binary_link[upper.tri(binary_link)]
    weight_link_value <- weight_link[upper.tri(weight_link)]
    edge <- cbind("Node1" = i, "Node2" = k, "Binary" = binary_link_value, "Weight" = weight_link_value)
    edge_dn <- edge[which(edge[,3] != 0),]

    if (is.null(p_val) == TRUE) {
        # Calculate p-values using logistic regression if p-values are not provided by users
        pvalue <- pvalue_logit(data_list$data, data_list$class_label, data_list$id)
        p.value <- pvalue$p.value
        row.names(pvalue)<-NULL
    } else {     # If the p-value matrix is provided
        pvalue <- p_val
        p.value <- pvalue$p.value           # Extract p-values from the table provided
        row.names(pvalue)<-NULL

    }

    # trasfer p-value to z-score
    z_score <- abs(qnorm(1 - p.value/2))
    # calculate differntial network score
    dn_score <- compute_dns(binary_link, z_score)


    indeed_df <- cbind(pvalue, rowSums(abs(binary_link)), dn_score )

    colnames(indeed_df) <- c("MetID", "P_value", "Node Degree", "Activity_Score")
    indeed_df$P_value <- lapply(indeed_df$P_value, round, 3)
    indeed_df$Activity_Score <- lapply(indeed_df$Activity_Score, round, 1)
    indeed_df <- as.data.frame(lapply(indeed_df, unlist))
    indeed_df<-indeed_df[order(indeed_df$Activity_Score, decreasing=TRUE), ]

    result_list <-list(activity_score=indeed_df,diff_network=edge_dn)
  }
}


#' @title Non-partial correlaton analysis
#' @description A method that integrates differential expression (DE) analysis
#'   and differential network (DN) analysis to select biomarker candidates for
#'   survival time prediction. non_partial_cor is a one step function for user
#'   to perform analysis, no pre-processing step required
#' @param data, input matrix of expression from all metabolites from all samples
#' @param class_label, a binary array with 0: group 1; 1: group 2.
#' @param id, an array of biomolecule ID to label.
#' @param method a character string indicating which correlation coefficient is
#'    to be computed. One of "pearson" (default) or "spearman".
#' @param p_val optional, a dataframe contains p values for each metabolite/molecule
#' @param permutation, a positive integer of desired number of permutations, default 1000
#' @examples non_partial_cor(data=Met_GU,class_label = Met_Group_GU,id=Met_name_GU,
#'    method="spearman")
#' @return a list of processed data for next step and rho
#' @import devtools
#' @importFrom glasso glasso
#' @importFrom stats qnorm cor quantile var sd glm
#' @importFrom graphics abline title plot lines
#' @export

non_partial_cor <- function(data = NULL, class_label = NULL, id = NULL, method = "pearson",  p_val = NULL,permutation=1000) {
  data_bind <- rbind(data , class_label)
  raw_group_1 <- data_bind[,data_bind[nrow(data_bind),] == 0][1:(nrow(data_bind) - 1),]  # Group 1: p*n1
  raw_group_2 <- data_bind[,data_bind[nrow(data_bind),] == 1][1:(nrow(data_bind) - 1),]  # Group 2: p*n2

  p <- nrow(raw_group_2)
  n_group_2 <- ncol(raw_group_2)
  n_group_1 <- ncol(raw_group_1)

  # Z-transform the data for group-specific normalization
  data_group_1 <- scale(t(raw_group_1)) # Group 1: n1*p
  data_group_2 <- scale(t(raw_group_2)) # Group 2: n2*p
  cov_group_1 <- var(data_group_1)
  cov_group_2 <- var(data_group_2)
  if(missing(method)){method="pearson"}
  else if(method != "spearman"){method ="pearson"}


  cor <- compute_cor(data_group_2, data_group_1, type_of_cor = method)    # default is pearson correlation
  # Get the correlation matrix
  cor_group_2 <- cor$Group2
  cor_group_1 <- cor$Group1

  # examine the correlation matrix
  thres <- 1e-3
  sum(abs(cor_group_2) > thres)
  cor_group_2[1:10, 1:10]
  sum(abs(cor_group_1) > thres)
  cor_group_1[1:10, 1:10]
  rm(thres)

  # Build differential correlation networks
  diff <- cor_group_2 - cor_group_1 # from group 1 to group 2
  thres = 1e-3
  sum(abs(diff) > thres)
  diff[1:10, 1:10]

  # Permutation test
  if(permutation<=0) {stop("please provide a valid number of permutation (positive integer)")}
  else{
    m <- as.numeric(permutation)
    diff_p <- permutation_cor(m, p, n_group_1, n_group_2, data_group_1, data_group_2, type_of_cor = method)
  }


  #####
  thres_left <- 0.025
  thres_right <- 0.975
  significant_thres <- permutation_thres(thres_left, thres_right, p, diff_p)
  rm(thres_left, thres_right)

  # get binary matrix
  significant_thres_p <- significant_thres$positive
  significant_thres_n <- significant_thres$negative
  binary_link <- matrix(0, p, p) # binary connection
  binary_link[diff < significant_thres_n] <- -1
  binary_link[diff > significant_thres_p] <- 1
  weight_link <- matrix(0, p, p) # weight connection
  weight_link[diff < significant_thres_n] <- diff[diff < significant_thres_n]
  weight_link[diff > significant_thres_p] <- diff[diff > significant_thres_p]
  sum(diff < significant_thres_n)
  sum(diff > significant_thres_p)
  binary_link[1:10, 1:10]
  weight_link[1:10, 1:10]
  rowSums(abs(binary_link)) # node degree for differential networks
  rm(diff_p)
  # Convert adjacent matrix into edge list
  i <- rep(seq_len(nrow(binary_link) - 1), times = (nrow(binary_link)-1):1)
  k <- unlist(lapply(2:nrow(binary_link), seq, nrow(binary_link)))
  binary_link_value <- binary_link[upper.tri(binary_link)]
  weight_link_value <- weight_link[upper.tri(weight_link)]
  edge <- cbind("Node1" = i, "Node2" = k, "Binary" = binary_link_value, "Weight" = weight_link_value)
  edge_dn <- edge[which(edge[,3] != 0),]

  if (is.null(p_val) == TRUE) {
    # Calculate p-values using logistic regression if p-values are not provided by users
    pvalue <- pvalue_logit(data, class_label, id)
    p.value <- pvalue$p.value
    row.names(pvalue)<-NULL
  } else {     # If the p-value matrix is provided
    pvalue <- p_val
    p.value <- pvalue$p.value           # Extract p-values from the table provided
    row.names(pvalue)<-NULL
  }

  # trasfer p-value to z-score
  z_score <- abs(qnorm(1 - p.value/2))
  # calculate differntial network score
  dn_score <- compute_dns(binary_link, z_score)
  indeed_df <- cbind(pvalue, rowSums(abs(binary_link)), dn_score )

  colnames(indeed_df) <- c("MetID", "P_value", "Node Degree", "Activity_Score")
  indeed_df$P_value <- lapply(indeed_df$P_value, round, 3)
  indeed_df$Activity_Score <- lapply(indeed_df$Activity_Score, round, 1)
  indeed_df <- as.data.frame(lapply(indeed_df, unlist))
  indeed_df<-indeed_df[order(indeed_df$Activity_Score, decreasing=TRUE), ]
  result_list <-list(activity_score=indeed_df,diff_network=edge_dn)
  return(result_list)
}
