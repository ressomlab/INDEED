#' @title Interactive Network Visualization
#' @description An interactive tool to assist in the visualization of the results from INDEED functions patial_corr() or non_partial_corr().
#' The size and the color of each node can be adjusted by users to represent either the Node_Degree, Activity_Score, Z_Score, or P_Value.
#' The color of the edge is based on the binary value of either 1 corresonding to a positive correlation dipicted as green or a 
#' negative correlation of -1 dipicted as red. The user also has the option of having the width of each edge be proportional to its weight value. 
#' The layout of the network can also be customized by choosing from the options: 'nice', 'sphere', 'grid', 'star', and circle'. Nodes can be 
#' moved and zoomed in on. Each node and edge will display extra information when clicked on. Secondary interactions will be highlighted 
#' as well when a node is clicked on. 
#' @param results This is the result from calling either partial_corr() or non_partial_corr()
#' @param nodesize This parameter determines what the size of each node will represent. The options are 
#' 'Node_Degree', 'Activity_Score', 'Z_Score', and 'P_Value'. The title of the resulting network will identify which parameter 
#' was selected to represent the node size. The default is Node_Degree
#' @param nodecolor This parameter determines what color each node will be based on a yellow to blue color gradient.  
#' The options are 'Node_Degree', 'Activity_Score', 'Z_Score', and 'P_Value'. A color bar will be created based on which parameter is chosen. 
#' @param edgewidth This is a 'YES' or 'NO' option as to if the edgewidth should be representative of the weight value corresponding 
#' to the correlation between two nodes. 
#' @param layout User can choose from a a handful of network visualization templates including:'nice', 'sphere', 'grid', 'star', and circle'.  
#'
#' @examples result1 = non_partial_cor(data=Met_GU,class_label = Met_Group_GU,
#'                     id=Met_name_GU,method="spearman",permutation_thres = 0.05, permutation = 1000)
#'           network_display(results = result1, layout= 'nice', nodesize= 'Node_Degree', 
#'                           nodecolor= 'Activity_Score', edgewidth= 'NO')
#' @return An interactive dipiction of the network resulting from INDEED functions patial_corr() or non_partial_corr() 
#' @import igraph
#' @import visNetwork
#' @importFrom grDevices topo.colors 
#' @export

network_display <- function(results = NULL, nodesize= 'Node_Degree', nodecolor= 'Activity_Score', edgewidth= 'NO', layout= 'nice'){
  
  nodes <- results$activity_score
  links <- results$diff_network
  vis.nodes <- data.frame(id= nodes$Node, name= nodes$ID, pval= nodes$P_value, ndegree= nodes$Node_Degree, ascore= nodes$Activity_Score, stringsAsFactors = FALSE)
  vis.links <- data.frame(from=links$Node1, to=links$Node2, binary= links$Binary, weight= links$Weight)
  
  vis.nodes$shape  <- "dot"  
  vis.nodes$shadow <- TRUE # Nodes will drop shadow
  # Information that will be displayed when hovering over a node 
  vis.nodes$title <- paste0("<p>", paste('MetID: ', vis.nodes$name), "<br>","<br>", paste('Node Degree: ', vis.nodes$ndegree),"<br>", paste('Activity Score: ', vis.nodes$ascore),"<br>", paste('P-value: ', vis.nodes$pval), "</p>") 
  
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
    pvalNorm <- scale_range(abs(vis.nodes$pval))
    vis.nodes$size   <- (10^(pvalNorm))*5
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
    vis.nodes$color.background <- topo.colors(length(vis.nodes$ascore), alpha=1)
    vis.nodes$color.highlight.background <- topo.colors(length(vis.nodes$ascore), alpha=1)
    nodecolor <- 'Activity_Score'
  }
  else if (nodecolor == 'Node_Degree'){
    vis.nodes$color.background <- topo.colors(length(vis.nodes$ndegree), alpha=1)
    vis.nodes$color.highlight.background <- topo.colors(length(vis.nodes$ndegree), alpha=1)
    
  } 
  else if(nodecolor == 'Activity_Score'){
    vis.nodes$color.background <- topo.colors(length(vis.nodes$ascore), alpha=1)
    vis.nodes$color.highlight.background <- topo.colors(length(vis.nodes$ascore), alpha=1)
    
  }
  else if(nodecolor == 'P_Value'){
    vis.nodes$color.background <- topo.colors(length(vis.nodes$pval), alpha=1)
    vis.nodes$color.highlight.background <- topo.colors(length(vis.nodes$pval), alpha=1)
    
  }
  else if(nodecolor == 'Z_Score'){
    z_score <- abs(qnorm(1 - (vis.nodes$pval)/2)) # trasfer p-value to z-score
    vis.nodes$color.background <- topo.colors(length(z_score), alpha=1)
    vis.nodes$color.highlight.background <- topo.colors(length(z_score), alpha=1)
    
  } else {
    vis.nodes$color.background <- topo.colors(length(vis.nodes$ascore), alpha=1)
    vis.nodes$color.highlight.background <- topo.colors(length(vis.nodes$ascore), alpha=1)
    nodecolor <- 'Activity_Score'
  }
  
  
  
  vis.nodes$color.background <- topo.colors(length(vis.nodes$ascore), alpha=1)
  vis.nodes$color.highlight.background <- topo.colors(length(vis.nodes$ascore), alpha=1)
  
  vis.nodes$borderWidth <- 2 # Node border width
  vis.nodes$label  <- vis.nodes$name  # Node label
  vis.nodes$color.highlight.border <- "darkred"
  vis.nodes$color.border <- "black"
  
  
  wNorm <- scale_range(abs(vis.links$weight))
  
  # Setting up edge width parameter 
  if (edgewidth != "NO" ){vis.links$width <- 15^(wNorm)
  } else {vis.links$width <- 3}
  # Information that will be displayed when hovering over the edge
  vis.links$title <- paste0("<p>", paste('Edge Weight: ', round(abs(vis.links$weight), digits= 3)), "</p>")
  vis.links$color[vis.links$binary == 1] <- "green"    # line color  
  vis.links$color[vis.links$binary == -1] <- "red" 
  vis.links$arrowStrikethrough <- FALSE
  vis.links$smooth <- TRUE    # should the edges be curved?
  vis.links$shadow <- TRUE    # edge shadow
  
  # Setting up layout of network 
  if (missing(layout)){l <- layout_nicely(net)}
  else if (layout == 'nice'){l <- "layout_nicely"}
  else if (layout == 'sphere'){l <- "layout_on_sphere"}
  else if (layout == 'star'){l <- "layout_as_star"}
  else if (layout == 'grid'){l <- "layout_on_grid"}
  else if (layout == 'circle'){l <- "layout_in_circle"
  } else {l <- "layout_nicely"}
  
  lnodes <- data.frame(label= c("High", "Mild", "Low"), shape= c("circle"), color= c("blue", "green", "yellow"))
  ledges <- data.frame(color= c("green", "red"), label= c("Positive Correlation", "Negative Correlation"), font.align= "top", arrows= c("NA", "NA"))
  
  
  
  
  
  
  net <- visNetwork(vis.nodes, vis.links, width = "100%", height = "800px", main= "INDEED 2.0", submain= paste("Node size is reprentative of: ", nodesize)) %>%
    visOptions( highlightNearest= TRUE, nodesIdSelection= TRUE)  %>% # selectedBy = (list(variable= "ndegree", multiple = FALSE) %>%
    visIgraphLayout(layout=l) %>%
    visInteraction( dragView= TRUE, dragNodes= TRUE, zoomView= TRUE, navigationButtons= FALSE, hideEdgesOnDrag= FALSE, multiselect = TRUE) %>%
    visLegend(addEdges= ledges, addNodes= lnodes, position= "right", useGroups= FALSE, ncol=1, main= paste("Node color based on ", nodecolor))
  print(net)
} 





