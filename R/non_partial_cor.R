#' @title Non-partial correlation analysis
#' @description A method that integrates differential expression (DE) analysis and differential 
#'     network (DN) analysis to select biomarker candidates for cancer studies. non_partial_cor is 
#'     a one step function for users to perform the typical correlation analysis. No pre-processing 
#'     step is required.
#' @param data This is a p*n dataframe that contains the expression levels for all biomolecules and samples.
#' @param class_label This is a 1*n dataframe that contains the class label with 0 for group 1 and 1 for group 2.
#' @param id This is a p*1 dataframe that contains the ID for each biomolecule.
#' @param method This is a character string indicating which correlation method is to use. The 
#'     options are either "pearson" as the default or "spearman".
#' @param p_val This is optional. It is a p*1 dataframe that contains the p-value for each biomolecule from DE analysis.
#' @param permutation This is a positive integer representing the desired number of permutations. 
#'     The default is 1000.
#' @param permutation_thres This is a integer representing the threshold for the permutation test. 
#'     The default is 0.05 to achieve 95 percent confidence.
#' @param fdr This is a boolean value indicating whether to apply multiple testing correction (TRUE) 
#'     or not (FALSE). The default is TRUE. However, if users find the output network is too sparse 
#'     even after relaxing the permutation_thres, it's probably a good idea to turn off the multiple testing correction.
#' @examples non_partial_cor(data = Met_GU, class_label = Met_Group_GU, id = Met_name_GU,
#'                           method = "pearson", p_val = pvalue_M_GU, permutation = 1000, 
#'                           permutation_thres = 0.05, fdr = TRUE)
#' @return A list containing an activity score dataframe with "ID", "P_value", "Node_Degree" and 
#'     "Activity_Score" as columns and a differential network dataframe with the binary and the 
#'     weight connection values.
#' @import devtools
#' @importFrom glasso glasso
#' @importFrom stats qnorm cor quantile var sd glm
#' @importFrom graphics abline title plot lines
#' @export

non_partial_cor <- function(data = NULL, class_label = NULL, id = NULL, method = "pearson",
                            p_val = NULL, permutation = 1000, permutation_thres = 0.05, fdr = TRUE){
    data_bind <- rbind(data, class_label)
    # Group 1: p*n1
    raw_group_1 <- data_bind[,data_bind[nrow(data_bind),] == 0][1:(nrow(data_bind) - 1),]  
    # Group 2: p*n2
    raw_group_2 <- data_bind[,data_bind[nrow(data_bind),] == 1][1:(nrow(data_bind) - 1),]  
    p <- nrow(raw_group_1)
    n_group_1 <- ncol(raw_group_1)
    n_group_2 <- ncol(raw_group_2)

    # Z-transform the data for group-specific normalization
    data_group_1 <- scale(t(raw_group_1)) # Group 1: n1*p
    data_group_2 <- scale(t(raw_group_2)) # Group 2: n2*p
    cov_group_1 <- var(data_group_1)
    cov_group_2 <- var(data_group_2)
    
    # Get the correlation matrix
    if(missing(method)){method = "pearson"}
    else if(method != "spearman"){method = "pearson"}
    # default is pearson correlation
    cor <- compute_cor(data_group_1, data_group_2, type_of_cor = method) 
    cor_group_1 <- cor$Group1
    cor_group_2 <- cor$Group2
    
    # # examine the correlation matrix
    # thres <- 1e-3
    # sum(abs(cor_group_1) > thres)
    # cor_group_1[1:10, 1:10]
    # sum(abs(cor_group_2) > thres)
    # cor_group_2[1:10, 1:10]
    # rm(thres)

    # Build differential correlation networks
    diff <- cor_group_2 - cor_group_1 # from group 1 to group 2
    # thres = 1e-3
    # sum(abs(diff) > thres)
    # diff[1:10, 1:10]

    # Permutation test
    if(permutation <= 0) {stop("please provide a valid number of permutation (positive integer)")}
    else{
        m <- as.numeric(permutation)
        diff_p <- permutation_cor(m, p, n_group_1, n_group_2, data_group_1, data_group_2, 
                                  type_of_cor = method)
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

    # get binary and weight matrix
    binary_link <- matrix(0, p, p) # binary connection
    binary_link[pvalue_edge_fdr < permutation_thres] <- 1
    binary_link[(pvalue_edge_fdr < permutation_thres) & (diff < 0)] <- -1
    weight_link <- compute_edge_weights(pvalue_edge_fdr, binary_link)
    # binary_link[1:10, 1:10]
    # weight_link[1:10, 1:10]
    # rowSums(abs(binary_link)) # node degree for differential networks
    # rm(diff_p)
    
    # Convert adjacent matrix into edge list
    i <- rep(seq_len(nrow(binary_link) - 1), times = (nrow(binary_link)-1):1)
    k <- unlist(lapply(2:nrow(binary_link), seq, nrow(binary_link)))
    binary_link_value <- binary_link[lower.tri(binary_link)]
    weight_link_value <- weight_link[lower.tri(weight_link)]
    edge <- cbind("Node1" = i, "Node2" = k, "Binary" = binary_link_value, 
                  "Weight" = weight_link_value)
    edge_dn <- edge[which(edge[,3] != 0),]
    edge_dn <- as.data.frame(edge_dn)

    # Compute p-values
    if (is.null(p_val) == TRUE) {
        # Calculate p-values using logistic regression if p-values are not provided by users
        pvalue <- pvalue_logit(data, class_label, id)
        p.value <- pvalue$p.value
        row.names(pvalue) <- NULL
    } else {     # If the p-value matrix is provided
        pvalue <- p_val
        p.value <- pvalue$p.value   # Extract p-values from the table provided
        row.names(pvalue) <- NULL
    }

    # trasfer p-value to z-score
    z_score <- abs(qnorm(1 - p.value/2))
    
    # calculate differntial network score
    dn_score <- compute_dns(binary_link, z_score)
    indeed_df <- cbind(pvalue, rowSums(abs(binary_link)), dn_score)
    colnames(indeed_df) <- c("ID", "P_value", "Node_Degree", "Activity_Score")
    indeed_df$P_value <- lapply(indeed_df$P_value, round, 3)
    indeed_df$Activity_Score <- lapply(indeed_df$Activity_Score, round, 1)
    indeed_df <- as.data.frame(lapply(indeed_df, unlist))
    # Recopy dataframe with index to help with ighraph formating
    indeed_df <- cbind(rownames(indeed_df), data.frame(indeed_df, row.names = NULL)) 
    colnames(indeed_df)[1] <- "Node"    # rename the previous index column as "Node"
    indeed_df<-indeed_df[order(indeed_df$Activity_Score, decreasing = TRUE), ]
    row.names(indeed_df) <- NULL      # remove index repeat
    
    # return
    result_list <-list(activity_score = indeed_df, diff_network = edge_dn)
    return(result_list)
}

