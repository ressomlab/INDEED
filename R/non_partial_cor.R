#' @title Non-partial correlaton analysis
#' @description A method that integrates differential expression (DE) analysis
#'   and differential network (DN) analysis to select biomarker candidates for
#'   survival time prediction. non_partial_cor is a one step function for user
#'   to perform analysis, no pre-processing step required
#' @param data input matrix of expression from all metabolites from all samples
#' @param class_label a binary array with 0: group 1; 1: group 2.
#' @param id an array of biomolecule ID to label.
#' @param method a character string indicating which correlation coefficient is
#'    to be computed. One of "pearson" (default) or "spearman".
#' @param p_val optional, a dataframe contains p values for each metabolite/molecule
#' @param permutation a positive integer of desired number of permutations, default 1000
#' @param permutation_thres threshold for permutation, defalut is 0.025 for each side to make 95percent
#' @examples non_partial_cor(data=Met_GU,class_label = Met_Group_GU,id=Met_name_GU,
#'    method="spearman",permutation_thres=0.05,permutation=1000)
#' @return a list of processed data for next step and rho
#' @import devtools
#' @importFrom glasso glasso
#' @importFrom stats qnorm cor quantile var sd glm
#' @importFrom graphics abline title plot lines
#' @export

non_partial_cor <- function(data = NULL, class_label = NULL, id = NULL, method = "pearson",  p_val = NULL,permutation=1000,permutation_thres=0.025) {
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
    thres_left <- permutation_thres
    thres_right <- 1-permutation_thres
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
    binary_link_value <- binary_link[lower.tri(binary_link)]
    weight_link_value <- weight_link[lower.tri(weight_link)]
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

