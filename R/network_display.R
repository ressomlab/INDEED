#' @title Interactive Network Visualization
#' @description An interactive tool to assist in the visualization of the results from INDEED 
#'     functions non_partial_corr() or patial_corr(). The size and the color of each node can be 
#'     adjusted by users to represent either the Node_Degree, Activity_Score, Z_Score, or P_Value.
#'     The color of the edge is based on the binary value of either 1 corresonding to a positive 
#'     correlation dipicted as green or a negative correlation of -1 dipicted as red. The user also 
#'     has the option of having the width of each edge be proportional to its weight value. The 
#'     layout of the network can also be customized by choosing from the options: 'nice', 'sphere', 
#'     'grid', 'star', and 'circle'. Nodes can be moved and zoomed in on. Each node and edge will 
#'     display extra information when clicked on. Secondary interactions will be highlighted as 
#'     well when a node is clicked on. 
#' @param results This is the result from calling either non_partial_corr() or partial_corr(). 
#' @param nodesize This parameter determines what the size of each node will represent. The options 
#'     are 'Node_Degree', 'Activity_Score','P_Value' and 'Z_Score'. The title of the resulting 
#'     network will identify which parameter was selected to represent the node size. The default 
#'     is P_Value.
#' @param nodecolor This parameter determines what color each node will be based on a yellow to 
#'     blue color gradient. The options are 'Node_Degree', 'Activity_Score', 'P_Value', and '
#'     Z_Score'. A color bar will be created based on which parameter is chosen. The default is
#'     Activity_Score.
#' @param edgewidth This is a 'YES' or 'NO' option as to if the edgewidth should be representative 
#'     of the weight value corresponding to the correlation change between two nodes. The default 
#'     is NO.
#' @param layout User can choose from a a handful of network visualization templates including:
#'     'nice', 'sphere', 'grid', 'star', and 'circle'. The default is nice. 
#' @examples result = non_partial_cor(data = Met_GU, class_label = Met_Group_GU, id = Met_name_GU, 
#'                                    method = "spearman", permutation_thres = 0.05, 
#'                                    permutation = 1000)
#'           network_display(results = result, nodesize = 'P_Value', 
#'           nodecolor = 'Activity_Score', edgewidth = 'NO', layout = 'nice')
#' @return An interactive dipiction of the network resulting from INDEED functions 
#'     non_partial_corr() or patial_corr().
#' @import igraph
#' @import visNetwork
#' @importFrom grDevices topo.colors 
#' @export

network_display <- function(results = NULL, nodesize= 'P_Value', nodecolor= 'Activity_Score', 
                            edgewidth= 'NO', layout= 'nice'){
    
    nodes <- results$activity_score
    links <- results$diff_network
   
    # Adding Z_Score to dataframe
    Z_Score <- abs(qnorm(1 - (nodes$P_value)/2)) # trasfer p-value to z-score
    nodes$zscore = Z_Score

    vis.nodes <- data.frame(id= nodes$Node, name= nodes$ID,font.size = 24, pval= nodes$P_value, 
                            ndegree= nodes$Node_Degree, ascore= nodes$Activity_Score,
                            zscore= nodes$zscore, stringsAsFactors = FALSE)
    vis.links <- data.frame(from=links$Node1, to=links$Node2, binary= links$Binary, 
                            weight= links$Weight)
    
    vis.nodes$shape  <- "dot"  
    vis.nodes$shadow <- TRUE # Nodes will drop shadow
    # Information that will be displayed when hovering over a node 
    vis.nodes$title <- paste0("<p>", paste('ID: ', vis.nodes$name), "<br>","<br>", 
                              paste('Node Degree: ', vis.nodes$ndegree),"<br>", 
                              paste('Activity Score: ', vis.nodes$ascore),"<br>", 
                              paste('P-value: ', vis.nodes$pval),"<br>",
                              paste('Z-score: ', round(vis.nodes$zscore, digits=3)),"</p>")
    
    # Setting up the Node Size
    if (missing(nodesize)){
      vis.nodes$size   <- (rank(-1 * vis.nodes$pval)+1)
        nodesize <- 'p-value significance '
    }
    else if (nodesize == 'Node_Degree'){
        vis.nodes$size   <- ((vis.nodes$ndegree)+1)*5
        nodesize <- 'Node Degree'
    } 
    else if(nodesize == 'Activity_Score'){
        vis.nodes$size   <- ((vis.nodes$ascore)+1)*5
        nodesize <- 'Activity Score'
    }
    else if(nodesize == 'P_Value'){
        vis.nodes$size   <- (rank(-1 * vis.nodes$pval)+1)
        nodesize <- 'p-value significance'
    }
    else if(nodesize == 'Z_Score'){
        vis.nodes$size   <- ((vis.nodes$zscore)+1)*10
        nodesize <- 'Z-Score'
    } else {
        vis.nodes$size   <- (rank(-1 * vis.nodes$pval)+1)
        nodesize <- 'p-value significance'
    }
    
    H <- "Higher"
    L <- "Lower"
    M <- "       "
    
    # Setting Up Node Color
    if (missing(nodecolor)){
        vis.nodes<- vis.nodes[order(vis.nodes$ascore, decreasing=TRUE), ]
        vis.nodes$color.background <- topo.colors(length(vis.nodes$ascore), alpha=1)
        vis.nodes$color.highlight.background <- topo.colors(length(vis.nodes$ascore), alpha=1)
        nodecolor <- 'Activity Score'
    }
    else if (nodecolor == 'Node_Degree'){
        vis.nodes<- vis.nodes[order(vis.nodes$ndegree, decreasing=TRUE), ]
        vis.nodes$color.background <- topo.colors(length(vis.nodes$ndegree), alpha=1)
        vis.nodes$color.highlight.background <- topo.colors(length(vis.nodes$ndegree), alpha=1)
        nodecolor <- 'Node Degree'
    } 
    else if(nodecolor == 'Activity_Score'){
        vis.nodes<- vis.nodes[order(vis.nodes$ascore, decreasing=TRUE), ]
        vis.nodes$color.background <- topo.colors(length(vis.nodes$ascore), alpha=1)
        vis.nodes$color.highlight.background <- topo.colors(length(vis.nodes$ascore), alpha=1)
        nodecolor <- 'Activity Score'
    }
    else if(nodecolor == 'P_Value'){
        vis.nodes<- vis.nodes[order(vis.nodes$pval, decreasing=FALSE), ]
        vis.nodes$color.background <- topo.colors(length(vis.nodes$pval), alpha=1)
        vis.nodes$color.highlight.background <- topo.colors(length(vis.nodes$pval), alpha=1)
        nodecolor <- 'p-value'
        H <- "Lower"
        L <- "Higher"
        M <- "       "
    }
    else if(nodecolor == 'Z_Score'){
      vis.nodes<- vis.nodes[order(vis.nodes$zscore, decreasing=TRUE), ]
      vis.nodes$color.background <- topo.colors(length(vis.nodes$zscore), alpha=1)
      vis.nodes$color.highlight.background <- topo.colors(length(vis.nodes$zscore), alpha=1)
      nodecolor <- 'Z-Score'
    }
    else {
        vis.nodes<- vis.nodes[order(vis.nodes$ascore, decreasing=TRUE), ]
        vis.nodes$color.background <- topo.colors(length(vis.nodes$ascore), alpha=1)
        vis.nodes$color.highlight.background <- topo.colors(length(vis.nodes$ascore), alpha=1)
        nodecolor <- 'Activity Score'
    }
    
    vis.nodes$borderWidth <- 2 # Node border width
    vis.nodes$label  <- vis.nodes$name # Node label
    vis.nodes$color.highlight.border <- "darkred"
    vis.nodes$color.border <- "black"

    # Setting up edge width parameter 
    if (edgewidth != "NO" ){vis.links$width <- abs(vis.links$weight) * 3
    } else {vis.links$width <- 3}
    # Information that will be displayed when hovering over the edge
    vis.links$title <- paste0("<p>", paste('Edge Weight: ', 
                                           round(abs(vis.links$weight), digits= 3)), "</p>")
    vis.links$color[vis.links$binary == 1] <- "green"    # line color  
    vis.links$color[vis.links$binary == -1] <- "red" 
    vis.links$arrowStrikethrough <- FALSE
    vis.links$smooth <- TRUE    # should the edges be curved?
    vis.links$shadow <- TRUE    # edge shadow
    
    # Setting up layout of network 
    if (missing(layout)){l <- "layout_nicely"}
    else if (layout == 'nice'){l <- "layout_nicely"}
    else if (layout == 'sphere'){l <- "layout_on_sphere"}
    else if (layout == 'star'){l <- "layout_as_star"}
    else if (layout == 'grid'){l <- "layout_on_grid"}
    else if (layout == 'circle'){l <- "layout_in_circle"
    } else {l <- "layout_nicely"}

    lnodes <- data.frame(label= c(H, M, L), position = "left", shape= c("circle"), 
                         color= c("blue", "green", "yellow"))
    ledges <- data.frame(color= c("green", "red"), label= c("Positive Change in Correlation", 
                                                            "Negative Change in Correlation"), 
                         font.align= "top", arrows= c("NA", "NA"), width= 4)
    
    net <- visNetwork(vis.nodes, vis.links, width = "100%", height = "800px", main= "INDEED", 
                      submain= paste("Node size represents: ", nodesize)) %>%
        visOptions( highlightNearest= TRUE, nodesIdSelection= TRUE)  %>% 
        visIgraphLayout(layout=l) %>%
        visInteraction( dragView= TRUE, dragNodes= TRUE, zoomView= TRUE, navigationButtons= FALSE, 
                        hideEdgesOnDrag= FALSE, multiselect = TRUE) %>%
        visLegend(addEdges= ledges, addNodes= lnodes, position= "right", useGroups= FALSE, ncol=1, 
                  main= paste("Node color based on ", nodecolor), width= 0.2, zoom = TRUE)
    print(net)
}



