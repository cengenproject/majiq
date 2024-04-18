
#' Plot an hclust with a color-coding of clusters
#'
#' @param hc hclust object
#' @param h the height for `cutree`
#' @param pal a palette function (that takes `n` as parameter)
#' @param cex the text `cex` factor
#' @param ... other parameters passed to `plot.hclust`
#'
#' @return the `hc` object
#'
#' @examples
plot_hc_clustered <- function(hc, h = .25, pal = ggsci::pal_d3("category20"), cex = 1, ...){
  clust <- cutree(hc, h = h)
  
  tab1 <- table(clust)
  
  nonclust_names <- names(tab1)[tab1 == 1]
  clust_names <- names(tab1)[tab1 > 1]
  
  
  i <- 0; repeat{
    
    palette <- pal(n = length(clust_names) + i)
    too_grey <- which(col2rgb(palette) |> matrixStats::colSds() < 5)
    palette <- palette[-too_grey]
    
    if(length(palette) == length(clust_names)){
      break
    } else{
      i <- i+length(too_grey)
    }
  }
  
  
  
  cols_by_clustname <- c(rep('grey', length(nonclust_names)),
            palette) |>
    setNames(c(nonclust_names, clust_names))
  
  cols_by_node <- cols_by_clustname[as.character(clust)[hc$order]]
  
  cols_by_node2 <- cols_by_node |> replace(cols_by_node == "grey", rgb(0,0,0, alpha = 0))
  
  labs_in_order <- hc$labels[hc$order]
  
  
  
  plot(hc, hang = -1,labels = FALSE, ...)
  
  for(i in seq_along(cols_by_node)){
    rect(i-.5, -.15, i+.5, -.05,col = cols_by_node2[[i]], border = NA)
    mtext(labs_in_order[[i]], side = 1, at = i, col = cols_by_node[[i]],
          cex = cex, las = 2)
  }
  
  return(invisible(hc))
}



