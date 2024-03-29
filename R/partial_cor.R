#' @title Partial correlation analysis
#' @description A method that integrates differential expression (DE) analysis and differential 
#'     network (DN) analysis to select biomarker candidates for cancer studies. partial_cor is the 
#'     second step of the partial correlation calculation after getting the result from select_rho_partial function.
#' @param data_list This is a list of pre-processed data outputted by the select_rho_partial function.
#' @param rho_group1 This is a character string indicating the rule for choosing rho value for group 1, 
#'     "min": minimum rho, "ste": one standard error from minimum, or user can input rho of their choice. The default 
#'     is minimum.
#' @param rho_group2 This is a character string indicating the rule for choosing rho value for group 2, 
#'     "min": minimum rho, "ste": one standard error from minimum, or user can input rho of their choice, the default 
#'     is minimum.
#' @param p_val This is optional. It is a p*1 dataframe that contains the p-value for each biomolecule from DE analysis.
#' @param permutation This is a positive integer representing the desired number of permutations. 
#'     The default is 1000.
#' @param permutation_thres This is a integer representing the threshold for the permutation test. 
#'     The default is 0.05 to achieve 95 percent confidence.
#' @param fdr This is a boolean value indicating whether to apply multiple testing correction (TRUE) 
#'     or not (FALSE). The default is TRUE. However, if users find the output network is too sparse 
#'     even after relaxing the permutation_thres, it's probably a good idea to turn off the multiple testing correction.
#' @examples
#' # step 1: select_rho_partial
#' pre_data <- select_rho_partial(data = Met_GU, class_label = Met_Group_GU, id = Met_name_GU, error_curve = TRUE)
#' # step 2: partial_cor
#' result <- partial_cor(data_list = pre_data, rho_group1 = 'min', rho_group2 = "min", p_val = pvalue_M_GU, permutation = 1000,
#'                       permutation_thres = 0.05, fdr = TRUE)
#' @return A list containing an activity score dataframe with "ID", "P_value", "Node_Degree" and 
#'     "Activity_Score" as columns and a differential network dataframe with the binary and the 
#'     weight connection values.
#' @import devtools
#' @importFrom glasso glasso
#' @importFrom stats qnorm cor quantile var sd glm
#' @importFrom graphics abline title plot lines par
#' @export

partial_cor <- function(data_list = NULL, rho_group1 = NULL, rho_group2 = NULL, p_val = NULL, 
                        permutation = 1000, permutation_thres = 0.05, fdr = TRUE){
    if(missing(data_list)) {stop("please provide data_list from select_rho_partial function")}
    else{
        # group 1
        if (rho_group1 =='min'){ rho_group_1_opt = data_list$rho_table["minimum",1] }
        else if (rho_group1 =='ste'){ rho_group_1_opt = data_list$rho_table["one standard error",1] }
        else if (is.numeric(rho_group1) & rho_group1>0) {rho_group_1_opt = rho_group1}
        else if (is.numeric(rho_group1) & rho_group1<=0) 
        {stop("please provide data_list from select_rho_partial function")}
        # default is minimum rho if no rule specified and no valid input entered
        else {rho_group_1_opt = data_list$rho_table["minimum",1]} 
        # group 2
        if (rho_group2 =='min'){ rho_group_2_opt = data_list$rho_table["minimum",2] }
        else if (rho_group2 =='ste'){ rho_group_2_opt = data_list$rho_table["one standard error",2] }
        else if (is.numeric(rho_group2) & rho_group2>0) {rho_group_2_opt = rho_group2}
        else if (is.numeric(rho_group2) & rho_group2<=0) 
        {stop("please provide data_list from select_rho_partial function")}
        # default is minimum rho if no rule specified and no valid input entered
        else {rho_group_2_opt = data_list$rho_table["minimum",2]} 

        ## Compute precision matrix for group 1
        pre_group_1 <- glasso(data_list$cov_group_1, rho = rho_group_1_opt)
        # thres <- 1e-3
        # sum(abs(pre_group_1$wi) > thres)
        # pre_group_1$wi[1:10, 1:10]
        
        ## Compute partial correlation for group 1
        pc_group_1 <- compute_par(pre_group_1$wi)
        # # examine the partial correlation matrix
        # sum(abs(pc_group_1) > thres)
        # pc_group_1[1:10, 1:10]
        
        ## Compute precision matrix for group 2
        pre_group_2 <- glasso(data_list$cov_group_2, rho = rho_group_2_opt)
        # # examine the precision matrix
        # sum(abs(pre_group_2$wi) > thres)
        # pre_group_2$wi[1:10,1:10]
        
        ## Compute partial correlation for group 2
        pc_group_2 <- compute_par(pre_group_2$wi)
        # # examine the partial correlation matrix
        # sum(abs(pc_group_2) > thres)
        # pc_group_2[1:10,1:10]
        
        ## Differential network
        diff <- pc_group_2 - pc_group_1  # from group 1 to group 2
        # thres = 1e-3
        # sum(abs(diff) > thres)
        # diff[1:10, 1:10]

        ## Permutation test using partial correlation
        if(permutation <= 0) 
            {stop("please provide a valid number of permutation (positive integer)")}
        else{
            m <- as.numeric(permutation)
            diff_p <- permutation_pc(m, data_list$p, data_list$n_group_1, data_list$n_group_2, 
                                     data_list$data_group_1, data_list$data_group_2, 
                                     rho_group_1_opt, rho_group_2_opt)
            p <- data_list$p
            ## Multiple testing step
            # p-value for edges
            pvalue_edge <- compute_pvalue_edge(p, diff, diff_p, m)
            # fdr to adjust multiple testing
            if(fdr == TRUE){
                pvalue_edge_fdr <- compute_pvalue_edge_fdr(p, pvalue_edge)
            }
            else{
                pvalue_edge_fdr <- pvalue_edge
            }
        }
        rm(m)

        ## Get binary and weight matrix
        binary_link <- matrix(0, p, p) # binary connection
        binary_link[pvalue_edge_fdr < permutation_thres] <- 1
        binary_link[(pvalue_edge_fdr < permutation_thres) & (diff < 0)] <- -1
        weight_link <- compute_edge_weights(pvalue_edge_fdr, binary_link)
        # binary_link[1:10, 1:10]
        # weight_link[1:10, 1:10]
        # rowSums(abs(binary_link)) # node degree for differential networks
        # rm(diff_p)

        ## Convert adjacent matrix into edge list
        i <- rep(seq_len(nrow(binary_link) - 1), times = (nrow(binary_link)-1):1)
        k <- unlist(lapply(2:nrow(binary_link), seq, nrow(binary_link)))
        binary_link_value <- binary_link[lower.tri(binary_link)]
        weight_link_value <- weight_link[lower.tri(weight_link)]
        edge <- cbind("Node1" = i, "Node2" = k, "Binary" = binary_link_value, 
                      "Weight" = weight_link_value)
        edge_dn <- edge[which(edge[,3] != 0),]
        edge_dn <- as.data.frame(edge_dn)

        ## Compute p-values
        if (is.null(p_val) == TRUE) {
            # calculate p-values using logistic regression if p-values are not provided by users
            pvalue <- pvalue_logit(data_list$data, data_list$class_label, data_list$id)
            p.value <- pvalue$p.value
            row.names(pvalue)<-NULL
        } else {     # if the p-value matrix is provided
            pvalue <- p_val
            p.value <- pvalue$p.value           # extract p-values from the table provided
            row.names(pvalue)<-NULL
        }

        ## Transfer p-value to z-score
        z_score <- abs(qnorm(1 - p.value/2))
        
        ## calculate differntial network score
        dn_score <- compute_dns(binary_link, z_score)
        indeed_df <- cbind(pvalue, rowSums(abs(binary_link)), dn_score )
        colnames(indeed_df) <- c("ID", "P_value", "Node_Degree", "Activity_Score")
        indeed_df$P_value <- lapply(indeed_df$P_value, round, 3)
        indeed_df$Activity_Score <- lapply(indeed_df$Activity_Score, round, 1)
        indeed_df <- as.data.frame(lapply(indeed_df, unlist))
        ## Recopy dataframe with index to help with ighraph formating
        indeed_df <- cbind(rownames(indeed_df) , data.frame(indeed_df, row.names=NULL) ) 
        colnames(indeed_df)[1] <- "Node"    # rename the previous index column as "Node"
        indeed_df<-indeed_df[order(indeed_df$Activity_Score, decreasing=TRUE), ]
        row.names(indeed_df) <- NULL      # remove index repeat

        ## Return
        result_list <-list(activity_score=indeed_df, diff_network=edge_dn)
        return (result_list)

    }
}
