# Find all the marker genes in parallel. This function can reduce the elapsed time to 1/6-1/3.
library(parallel)
library(Seurat)

FindMarkers_parallel <- function(object = NULL, mc.cores = NULL, ...){
  result <- mclapply(as.numeric(levels(object@active.ident)),
                     FUN =  function(x) {FindMarkers(object, ident.1 = x, ident.2 = NULL)},
                     mc.cores = mc.cores)
  RESULT <- result
  roundN <- 1
  while(any(mapply(length, result, SIMPLIFY = TRUE)!=5)){
    if(any(mapply(length, result, SIMPLIFY = TRUE)!=5)){
      recalculate_clusters <- which(mapply(length, result, SIMPLIFY = TRUE)!=5)-1
      print(recalculate_clusters)
      result1 <- mclapply(recalculate_clusters,
                          FUN =  function(x) {FindMarkers(object, ident.1 = x, ident.2 = NULL)},
                          mc.cores = mc.cores)
    }
    print(roundN + 1)
    for(i in 1:length(recalculate_clusters)){
      result[c(recalculate_clusters+1)[i]] <- result1[i]
    }
  }
  all_markers <- do.call(rbind, result)
  all_markers$gene <- unlist(mapply(rownames, result))
  all_markers$cluster <- rep(levels(object@active.ident), times = mapply(dim, result, SIMPLIFY = TRUE)[1,])
  subset_cells.markers <- all_markers
  return(subset_cells.markers)
}
  
 
