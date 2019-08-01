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
#'     are 'Node_Degree', 'Activity_Score', 'Z_Score', and 'P_Value'. The title of the resulting 
#'     network will identify which parameter was selected to represent the node size. The default 
#'     is Node_Degree.
#' @param nodecolor This parameter determines what color each node will be based on a yellow to 
#'     blue color gradient. The options are 'Node_Degree', 'Activity_Score', and 'P_Value'. A color 
#'     bar will be created based on which parameter is chosen. 
#' @param edgewidth This is a 'YES' or 'NO' option as to if the edgewidth should be representative 
#'     of the weight value corresponding to the correlation between two nodes. 
#' @param layout User can choose from a a handful of network visualization templates including:
#'     'nice', 'sphere', 'grid', 'star', and 'circle'.  
#' @param bingroups Users can choose between 'Activity_Score',  'Node_Degree', and 'P_Value' to 
#'     sort through varying threshold ranges.
#' @examples result = non_partial_cor(data = Met_GU, class_label = Met_Group_GU, id = Met_name_GU, 
#'                                    method = "spearman", permutation_thres = 0.05, 
#'                                    permutation = 1000)
#'           network_display(results = result, layout = 'nice', nodesize = 'Node_Degree', 
#'                           nodecolor = 'Activity_Score', edgewidth = 'NO', 
#'                           bingroups = 'Node_Degree')
#' @return An interactive dipiction of the network resulting from INDEED functions 
#'     non_partial_corr() or patial_corr().
#' @import igraph
#' @import visNetwork
#' @importFrom grDevices topo.colors 
#' @export

network_display <- function(results = NULL, nodesize= 'Node_Degree', nodecolor= 'Activity_Score', 
                            edgewidth= 'NO', layout= 'nice', bingroups= 'Node_Degree'){
    
    nodes <- results$activity_score
    links <- results$diff_network
    
    
    # Setting up groups for thresholding select option 
    if (missing(bingroups)){
        nodes$bins <- cut(nodes$Node_Degree, breaks=c(as.integer(min(nodes$Node_Degree)), 
                                                      as.integer(max(nodes$Node_Degree)/2), 
                                                      as.integer(max(nodes$Node_Degree)/1.5), 
                                                      as.integer(max(nodes$Node_Degree))))
        bingroups <- 'Node_Degree'
    }
    else if (bingroups == 'Node_Degree'){
        nodes$bins <- cut(nodes$Node_Degree, breaks=c(as.integer(min(nodes$Node_Degree)), 
                                                      as.integer(max(nodes$Node_Degree)/2), 
                                                      as.integer(max(nodes$Node_Degree)/1.5), 
                                                      as.integer(max(nodes$Node_Degree))))
    } 
    else if(bingroups == 'Activity_Score'){
        nodes$bins <- cut(nodes$Activity_Score, breaks=c(min(nodes$Activity_Score), 
                                                         max(nodes$Activity_Score)/2, 
                                                         max(nodes$Activity_Score)/1.5, 
                                                         max(nodes$Activity_Score)))
    } 
    else if(bingroups == 'P_Value'){
        nodes$bins <- cut(nodes$P_value, breaks=c(min(nodes$P_value), 
                                                  max(nodes$P_value)/2, 
                                                  max(nodes$P_value)/1.5, 
                                                  max(nodes$P_value)))
    }
    else { 
        nodes$bins <- cut(nodes$Node_Degree, breaks=c(min(nodes$Node_Degree),
                                                      max(nodes$Node_Degree)/2,
                                                      max(nodes$Node_Degree)/1.5,
                                                      max(nodes$Node_Degree)))
        bingroups <- 'Node_Degree'
        
    }
    
    vis.nodes <- data.frame(id= nodes$Node, name= nodes$ID, pval= nodes$P_value, 
                            ndegree= nodes$Node_Degree, ascore= nodes$Activity_Score, 
                            group = nodes$bins, stringsAsFactors = FALSE)
    vis.links <- data.frame(from=links$Node1, to=links$Node2, binary= links$Binary, 
                            weight= links$Weight)
    
    vis.nodes$shape  <- "dot"  
    vis.nodes$shadow <- TRUE # Nodes will drop shadow
    # Information that will be displayed when hovering over a node 
    vis.nodes$title <- paste0("<p>", paste('MetID: ', vis.nodes$name), "<br>","<br>", 
                              paste('Node Degree: ', vis.nodes$ndegree),"<br>", 
                              paste('Activity Score: ', vis.nodes$ascore),"<br>", 
                              paste('P-value: ', vis.nodes$pval), "</p>") 
    
    # Setting up the Node Size
    if (missing(nodesize)){
        vis.nodes$size   <- ((vis.nodes$ndegree)+1)*5
        nodesize <- 'Node_Degree'
    }
    else if (nodesize == 'Node_Degree'){
        vis.nodes$size   <- ((vis.nodes$ndegree)+1)*5
    } 
    else if(nodesize == 'Activity_Score'){
        vis.nodes$size   <- ((vis.nodes$ascore)+1)*5
    }
    else if(nodesize == 'P_Value'){
        vis.nodes$size   <- (rank(-1 * vis.nodes$pval)+1)
        
    }
    else if(nodesize == 'Z_Score'){
        z_score <- abs(qnorm(1 - (vis.nodes$pval)/2)) # trasfer p-value to z-score
        vis.nodes$size   <- ((z_score)+1)*10
        
    } else { 
        vis.nodes$size   <- ((vis.nodes$ndegree)+1)*5
        nodesize <- "Node_Degree"
    }
    
    
    
    # Setting Up Node Color
    if (missing(nodecolor)){
        vis.nodes<- vis.nodes[order(vis.nodes$ascore, decreasing=TRUE), ]
        vis.nodes$color.background <- topo.colors(length(vis.nodes$ascore), alpha=1)
        vis.nodes$color.highlight.background <- topo.colors(length(vis.nodes$ascore), alpha=1)
        nodecolor <- 'Activity_Score'
    }
    else if (nodecolor == 'Node_Degree'){
        vis.nodes<- vis.nodes[order(vis.nodes$ndegree, decreasing=TRUE), ]
        vis.nodes$color.background <- topo.colors(length(vis.nodes$ndegree), alpha=1)
        vis.nodes$color.highlight.background <- topo.colors(length(vis.nodes$ndegree), alpha=1)
        
    } 
    else if(nodecolor == 'Activity_Score'){
        vis.nodes<- vis.nodes[order(vis.nodes$ascore, decreasing=TRUE), ]
        vis.nodes$color.background <- topo.colors(length(vis.nodes$ascore), alpha=1)
        vis.nodes$color.highlight.background <- topo.colors(length(vis.nodes$ascore), alpha=1)
        
    }
    else if(nodecolor == 'P_Value'){
        vis.nodes<- vis.nodes[order(vis.nodes$pval, decreasing=FALSE), ]
        vis.nodes$color.background <- topo.colors(length(vis.nodes$pval), alpha=1)
        vis.nodes$color.highlight.background <- topo.colors(length(vis.nodes$pval), alpha=1)
        
    }
    else {
        vis.nodes<- vis.nodes[order(vis.nodes$ascore, decreasing=TRUE), ]
        vis.nodes$color.background <- topo.colors(length(vis.nodes$ascore), alpha=1)
        vis.nodes$color.highlight.background <- topo.colors(length(vis.nodes$ascore), alpha=1)
        nodecolor <- 'Activity_Score'
    }
    
    
    vis.nodes$borderWidth <- 2 # Node border width
    vis.nodes$label  <- vis.nodes$name  # Node label
    vis.nodes$color.highlight.border <- "darkred"
    vis.nodes$color.border <- "black"
    
    
    wNorm <- scale_range(abs(vis.links$weight))
    
    # Setting up edge width parameter 
    if (edgewidth != "NO" ){vis.links$width <- 15^(wNorm)
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
    
    lnodes <- data.frame(label= c("High", "       ", "Low"), shape= c("circle"), 
                         color= c("blue", "green", "yellow"))
    ledges <- data.frame(color= c("green", "red"), label= c("Positive Correlation", 
                                                            "Negative Correlation"), 
                         font.align= "top", arrows= c("NA", "NA"), width= 4)
    
    
    net <- visNetwork(vis.nodes, vis.links, width = "100%", height = "800px", main= "INDEED 2.0", 
                      submain= paste("Node size represents: ", nodesize, "&", 
                                     "Groups selection represents: ", bingroups)) %>%
        visOptions( highlightNearest= TRUE,selectedBy = "group", nodesIdSelection= TRUE)  %>% 
        visIgraphLayout(layout=l) %>%
        visInteraction( dragView= TRUE, dragNodes= TRUE, zoomView= TRUE, navigationButtons= FALSE, 
                        hideEdgesOnDrag= FALSE, multiselect = TRUE) %>%
        visLegend(addEdges= ledges, addNodes= lnodes, position= "right", useGroups= FALSE, ncol=1, 
                  main= paste("Node color based on ", nodecolor), width= 0.2, stepX = 50, 
                  stepY = 50, zoom = TRUE)
    print(net)
} 



