#' @title Data preprocessing for partial correlaton analysis
#' @description A method that integrates differential expression (DE) analysis
#'   and differential network (DN) analysis to select biomarker candidates for
#'   survival time prediction. partial_cor is the second step of partial correlation
#'   calculation after the output result from select_rho_partial function
#' @param data_list list of pre-processed data from select_rho_partial function
#' @param rho_group1 rule to choose rho for group 1, "min": minimum rho,
#' "ste" one standard error from minimum, or user can input rho of their choice, default: minimum
#' @param rho_group2 rule to choose rho for group 1, "min": minimum rho,
#' "ste" one standard error from minimum, or user can input rho of their choice, default: minimum
#' @param p_val optional, a dataframe contains p values for each metabolite/molecule
#' @param permutation a positive integer of desired number of permutations, default 1000
#' @param permutation_thres threshold for permutation, defalut is 0.025 for each side to make 95percent
#' @examples preprocess<- select_rho_partial(data=Met_GU,class_label =
#'    Met_Group_GU,id=Met_name_GU,error_curve="YES")
#'    partial_cor(data_list=preprocess,rho_group1='min',
#'    rho_group2="min",permutation = 1000,p_val=pvalue_M_GU,permutation_thres=0.05)
#' @return a list containing a score dataframe and a differential network dataframe
#' @import devtools
#' @importFrom glasso glasso
#' @importFrom stats qnorm cor quantile var sd glm
#' @importFrom graphics abline title plot lines par
#' @export
partial_cor <- function(data_list =NULL, rho_group1=NULL,rho_group2=NULL, permutation=1000, p_val = NULL,permutation_thres=0.025){
    if(missing(data_list)) {stop("please provide data_list from select_rho_partial function")}
    else{
        #group1
        if (rho_group1 =='min'){ rho_group_1_opt =data_list$rho_table[1,2] }
        else if (rho_group1 =='ste'){ rho_group_1_opt =data_list$rho_table[1,1] }
        else if (is.numeric(rho_group1) & rho_group1>0) {rho_group_1_opt = rho_group1}
        else if (is.numeric(rho_group1) & rho_group1<=0) {stop("please provide data_list from select_rho_partial function")}
        else {rho_group_1_opt =data_list$rho_table[1,2]} #default is minimum rho if no rule specified and no valid input entered
        #group2
        if (rho_group2 =='min'){ rho_group_2_opt =data_list$rho_table[2,2] }
        else if (rho_group2 =='ste'){ rho_group_2_opt =data_list$rho_table[2,1] }
        else if (is.numeric(rho_group2) & rho_group2>0) {rho_group_2_opt = rho_group2}
        else if (is.numeric(rho_group2) & rho_group2<=0) {stop("please provide data_list from select_rho_partial function")}
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
  
        ##### final calculation
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

        colnames(indeed_df) <- c("MetID", "P_value", "Node_Degree", "Activity_Score")
        indeed_df$P_value <- lapply(indeed_df$P_value, round, 3)
        indeed_df$Activity_Score <- lapply(indeed_df$Activity_Score, round, 1)
        indeed_df <- as.data.frame(lapply(indeed_df, unlist))
        
        indeed_df <- cbind(rownames(indeed_df) , data.frame(indeed_df, row.names=NULL) ) # Recopy dataframe with index to help with ighraph formating
        colnames(indeed_df)[1] <- "Node"    # rename the previous index column as "Node"
 
        indeed_df<-indeed_df[order(indeed_df$Activity_Score, decreasing=TRUE), ]
        row.names(indeed_df) <- NULL      # remove index repeat 

        result_list <-list(activity_score=indeed_df,diff_network=edge_dn)
        return (result_list)
    }
}
