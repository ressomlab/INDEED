#' @title Data preprocessing for partial correlation analysis
#' @description A method that integrates differential expression (DE) analysis
#'     and differential network (DN) analysis to select biomarker candidates for
#'     cancer studies. select_rho_partial is the pre-processing step for INDEED
#'     partial differential analysis.
#' @param data This is a p*n dataframe that contains the expression levels for all biomolecules and samples.
#' @param class_label This is a 1*n dataframe that contains the class label with 0 for group 1 and 1 for group 2.
#' @param id This is a p*1 dataframe that contains the ID for each biomolecule.
#' @param error_curve This is a boolean value indicating whether to plot the error curve (TRUE) or not (FALSE). 
#'     The default is TRUE.
#' @examples select_rho_partial(data = Met_GU, class_label = Met_Group_GU, id = Met_name_GU, 
#'     error_curve = TRUE)
#' @return A list of processed data for the next step, and an error curve to select optimal rho value
#'     for graphical lasso.
#' @import devtools
#' @importFrom glasso glasso
#' @importFrom stats qnorm cor quantile var sd glm 
#' @importFrom graphics abline title plot lines 
#' @export

select_rho_partial <- function(data = NULL, class_label = NULL, id = NULL, error_curve = TRUE){
    data_bind <- rbind(data, class_label)
    # Group 1: p*n1
    raw_group_1 <- data_bind[,data_bind[nrow(data_bind),] == 0][1:(nrow(data_bind) - 1),]
    # Group 2: p*n2
    raw_group_2 <- data_bind[,data_bind[nrow(data_bind),] == 1][1:(nrow(data_bind) - 1),]  
    p <- nrow(raw_group_2)
    n_group_1 <- ncol(raw_group_1)
    n_group_2 <- ncol(raw_group_2)
    

    # Z-transform the data for group-specific normalization
    data_group_1 <- scale(t(raw_group_1)) # Group 1: n1*p
    data_group_2 <- scale(t(raw_group_2)) # Group 2: n2*p
    cov_group_1 <- var(data_group_1)
    cov_group_2 <- var(data_group_2)
    
    ## Apply grapihcal LASSO given data_group_1 and data_group_2
    #  Group 1 first
    n_fold <- 5 # number of folds
    rho <- exp(seq(log(1), log(0.01), length.out = 20))

    # group 1
    # draw error curve
    error_group1 <- choose_rho(data_group_1, n_fold, rho)
    if(error_curve == TRUE){
        par(mfrow = c(1, 2))
        plot(rho, error_group1$log.cv, xlab = expression(rho), ylab = "Error")
        lines(rho, error_group1$log.cv)
        title(main = paste("Group 1 error curve", "using cross validation", sep = "\n"))
        # chosse optimal rho
        abline(v = rho[error_group1$log.cv == min(error_group1$log.cv)], col = "red", lty = 3)
        # one standard error rule
        abline(h = min(error_group1$log.cv) + 
                   error_group1$log.rho[error_group1$log.cv == min(error_group1$log.cv)], 
               col = "blue")
    }
    rho[error_group1$log.cv == min(error_group1$log.cv)] # rho based on minimum rule

    # rhos that are under blue line
    rho_under_blue <- rho[error_group1$log.cv < 
                              min(error_group1$log.cv) + 
                              error_group1$log.rho[error_group1$log.cv == 
                                                       min(error_group1$log.cv)]]

    # rhos that are on the right of the red line
    rho_right_red <- rho[rho > rho[error_group1$log.cv == min(error_group1$log.cv)]]

    #rhos that are under blue line and to the right of the red line
    rho_one_std_group1 <- max(intersect(rho_under_blue, rho_right_red ))
    rho_min_rule_group1 <- rho[error_group1$log.cv == min(error_group1$log.cv)]

    # group 2
    # draw error curve
    error_group2 <- choose_rho(data_group_2, n_fold, rho)
    if(error_curve == TRUE){
        plot(rho, error_group2$log.cv, xlab = expression(rho), ylab = "Error")
        lines(rho, error_group2$log.cv)
        title(main = paste("Group 2 error curve", "using cross validation", sep="\n"))
        # choOse optimal rho
        abline(v = rho[error_group2$log.cv == min(error_group2$log.cv)], col = "red", lty = 3)
        # one standard error rule
        abline(h = min(error_group2$log.cv) + 
                   error_group2$log.rho[error_group2$log.cv == 
                                            min(error_group2$log.cv)], col = "blue")
    }
    rho[error_group2$log.cv == min(error_group2$log.cv)] # rho based on minimum rule

    # rhos that are under blue line
    rho_under_blue <- rho[error_group2$log.cv < min(error_group2$log.cv) + 
                              error_group2$log.rho[error_group2$log.cv == 
                                                       min(error_group2$log.cv)]]

    # rhos that are on the right of the red line
    rho_right_red <- rho[rho > rho[error_group2$log.cv == min(error_group2$log.cv)]]

    #rhos that are under blue line and to the right of the red line
    rho_one_std_group2 <- max(intersect(rho_under_blue, rho_right_red )) # export as global
    rho_min_rule_group2 <- rho[error_group2$log.cv == min(error_group2$log.cv)] # export as global

    rho_df <- data.frame(c(rho_one_std_group1, rho_min_rule_group1), c(rho_one_std_group2, 
                                                                       rho_min_rule_group2), 
                         row.names = c("one standard error","minimum"))
    colnames(rho_df)<-c('group1', "group2")
    
    data_list <- list(p = p, rho_table = rho_df, cov_group_1 = cov_group_1, 
                      cov_group_2 = cov_group_2, n_group_1 = n_group_1, n_group_2 = n_group_2, 
                      data_group_1 = data_group_1, data_group_2 = data_group_2, 
                      class_label = class_label, id = id, data = data)
    return(data_list)
}
