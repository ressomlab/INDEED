#' @title Network Visualization for INDEED Partial and Non-Partial Correlation

#' @description A tool to assist in the visualization of the results from INDEED functions patial_cor()
#' and non_partial_cor(). The size and the color of each node can be adjusted by users to represent either
#' the Node_Degree, Activity_Score, Z_Score, or P_Value. The color of the edge is based on the binary value of
#' either 1 corresonding to a positive correlation dipicted as green or a negative correlation of -1 dipicted as red. The
#' user also has the option of having the width of each edge be proportional to its weight value. The layout of the network
#' can also be customized by choosing from the options: 'nice', 'sphere', 'grid', 'star', and 'circle'.
#' @param results This is the result from the calling either partial_cor() or non_partial_cor().
#' @param nodesize This parameter determines what the size of each node will represent. The options
#'        are 'Node_Degree', 'Activity_Score', 'Z_Score', and 'P_Value'. The title of the resulting network will
#'        identify which parameter was selected to represent the node size. The default is Node_Degree.
#' @param nodecolor This parameter determines what color each node will be, based on a yellow or red color gradient.
#'          The options are 'Node_Degree', 'Activity_Score', 'Z_Score', and 'P_Value'. A color bar will be created
#'          based on which parameter is chosen.
#' @param edgewidth	 This is a 'YES' or 'NO' option as to if the edgewidth should be representative of the weight value
#'          corresponding to the correlation between two nodes.
#' @param layout User can choose from a a handful of network visualization templates
#'            including:'nice', 'sphere', 'grid', 'star', and 'circle'.
#'
#' @examples
#' result1 = non_partial_cor(data = Met_GU, class_label = Met_Group_GU, id = Met_name_GU,
#'                           method = "spearman", permutation_thres = 0.05, permutation = 1000)
#' network_display(results = result1, layout = 'nice', nodesize = 'Node_Degree',
#'                 nodecolor = 'Activity_Score', edgewidth = 'NO')
#'
#' @return A visual dipiction of the network resulting from INDEED functions partial_cor() or non_partial_cor()
#' @import igraph
#' @importFrom grDevices rainbow col2rgb heat.colors rgb2hsv
#' @importFrom graphics axis legend rect

#' @export


network_display <- function(results = NULL, nodesize = 'Node_Degree', nodecolor = 'Activity_Score',
                            edgewidth = 'NO', layout = 'nice'){
	nodes <- results$activity_score
	links <- results$diff_network
	net <- graph_from_data_frame(d=links, vertices=nodes, directed=T)


	# Setting up the Node Size
	if (missing(nodesize)){V(net)$size <- V(net)$Node_Degree * 2}
	else if (nodesize == 'Node_Degree'){V(net)$size <- (V(net)$Node_Degree + 1) * 2 }
	else if(nodesize == 'Activity_Score'){V(net)$size <- V(net)$Activity_Score * 2}
	else if(nodesize == 'P_Value'){V(net)$size <- V(net)$P_value * 15}
	else if(nodesize == 'Z_Score'){

				z_score <- abs(qnorm(1 - (V(net)$P_value)/2)) # trasfer p-value to z-score
				V(net)$size <- z_score * 10
	} else {
		nodesize = 'Node_Degree'
		V(net)$size <- (V(net)$Node_Degree + 1) * 2

	}

	if (missing(nodecolor)){ncolor <- V(net)$Activity_Score}
	else if (nodecolor == 'Node_Degree'){ncolor <- V(net)$Node_Degree}
	else if(nodecolor == 'Activity_Score'){ncolor <- V(net)$Activity_Score}
	else if(nodecolor == 'P_Value'){ncolor <- V(net)$P_value}
	else if(nodecolor == 'Z_Score'){
			z_score <- abs(qnorm(1 - (V(net)$P_value)/2)) # trasfer p-value to z-score
			ncolor <- z_score
	} else {
			nodecolor = 'Activity_Score'
			ncolor <- V(net)$Activity_Score
	}

	V(net)$color [ncolor >= 0] <- heat.colors(length(ncolor), alpha=1)  # doesnt need to be >=0 make it something useful

	# Label 'MetID' display settings
	V(net)$label <- V(net)$MetID			# Adding MetID as label to each node
	V(net)$label.dist <- 0
	V(net)$label.cex <- 0.8
	V(net)$label.font <- 2
	V(net)$label.degree <- 0
	V(net)$label.color <- "black"


	# Edge display settings
	E(net)$color[E(net)$Binary == 1] <- 'dark green'
	E(net)$color[E(net)$Binary == -1] <- 'dark red'
	E(net)$arrow.size <- 0
	E(net)$arrow.width <- 0
	E(net)$curved <- 0.5

	# scale_range function spreads out Weight values to better distinguish edge widths in network display
	wNorm <- scale_range(abs(E(net)$Weight))
	# Setting up edge width parameter
	if (edgewidth != "NO" ){E(net)$width <- wNorm * 3
	} else {E(net)$width <- 2}

	# Setting up layout of network
	if (missing(layout)){l <- layout_nicely(net)}
	else if (layout == 'nice'){l <- layout_nicely(net)}
	else if (layout == 'sphere'){l <- layout_on_sphere(net)}
	else if (layout == 'star'){l <- layout_as_star(net)}
	else if (layout == 'grid'){l <- layout_on_grid(net)}
	else if (layout == 'circle'){l <- layout_in_circle(net)
	} else {l <- layout_nicely(net)}


	# creating color gradient
	lut <- rev(rainbow(100, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1]))


	# plotting
	par(mfrow=c(1,4))
	layout(matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2), nrow = 6, ncol = 4, byrow = TRUE))



	par(mar = c(1, 1, 1, 1))
	plot(net, layout=l, main= paste0("Node Size Representative of ", nodesize),  frame = FALSE)
	legend('bottomleft', c("Positive", "Negative"), lty = c(1,1), lwd=c(2.5,2.5), col=c("dark green", "dark red"), bty="o",cex=1.5, title= "Correlation")

	par(mar = c(2, 2, 2, 2))
	color.bar(lut, min=0, max=max(ncolor), title= toString(nodecolor), nticks=5)


}


