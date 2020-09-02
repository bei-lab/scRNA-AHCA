## monocle version 2.99
library(Rtsne)
library(monocle)

fixed_reduceDimension <- function (cds, max_components = 2, reduction_method = c("UMAP", "DDRTree"), auto_param_selection = TRUE, scaling = TRUE, 
    verbose = TRUE, ...){
    extra_arguments <- list(...)
    set.seed(2016)
    if (verbose) 
        message("Retrieving normalized data ...")
    FM <- cds@auxOrderingData$normalize_expr_data
    irlba_pca_res <- cds@normalized_data_projection
    if (is.null(FM)) {
        message("Warning: The cds has not been pre-processed yet. Running preprocessCDS() with default parameters.")
        cds <- preprocessCDS(cds)
        FM <- cds@auxOrderingData$normalize_expr_data
        irlba_pca_res <- cds@normalized_data_projection
    }
    if (is.function(reduction_method)) {
        if (scaling) {
            FM <- as.matrix(Matrix::t(scale(Matrix::t(FM))))
            FM <- FM[!is.na(row.names(FM)), ]
        }
        else FM <- as.matrix(FM)
        reducedDim <- reduction_method(FM, ...)
        colnames(reducedDim) <- colnames(FM)
        monocle:::reducedDimW(cds) <- as.matrix(reducedDim)
        monocle:::reducedDimA(cds) <- as.matrix(reducedDim)
        monocle:::reducedDimS(cds) <- as.matrix(reducedDim)
        monocle:::reducedDimK(cds) <- as.matrix(reducedDim)
        dp <- as.matrix(dist(reducedDim))
        cellPairwiseDistances(cds) <- dp
        gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
        dp_mst <- minimum.spanning.tree(gp)
        minSpanningTree(cds) <- dp_mst
        cds@dim_reduce_type <- "function_passed"
    }
    else {
        reduction_method <- match.arg(reduction_method)
        if (reduction_method == "ICA") {
            .Deprecated(msg = "ICA (used in monocle 1) is not supported anymore, please use RGE (reversed graph embedding) instead!")
        }
        else if (reduction_method == "tSNE") {
            if ("num_dim" %in% names(extra_arguments)) {
                num_dim <- extra_arguments$num_dim
            }
            else {
                num_dim <- 50
            }
            topDim_pca <- irlba_pca_res
            if (verbose) 
                message("Reduce dimension by tSNE ...")
            tsne_res <- Rtsne::Rtsne(as.matrix(topDim_pca), dims = max_components, 
                pca = F, check_duplicates = FALSE, ...)
            tsne_data <- tsne_res$Y[, 1:max_components]
            row.names(tsne_data) <- colnames(tsne_data)
            reducedDimA(cds) <- t(tsne_data)
            cds@auxClusteringData[["tSNE"]]$pca_components_used <- num_dim
            cds@dim_reduce_type <- "tSNE"
            pData(cds)$tsne_1 = reducedDimA(cds)[1, ]
            pData(cds)$tsne_2 = reducedDimA(cds)[2, ]
        }
        else if (reduction_method == c("DDRTree")) {
            message("DDRTree will be eventually deprecated in reduceDimension call and be used in learnGraph function instead. We are calling learnGraph for you now.")
            cds@reducedDimS <- t(cds@normalized_data_projection)
            cds <- partitionCells(cds)
            cds <- learnGraph(cds, rge_method = "DDRTree", do_partition = F, 
                ...)
        }
        else if (reduction_method == c("UMAP")) {
            if (verbose) 
                message("Running Uniform Manifold Approximation and Projection")
            umap_args <- c(list(X = irlba_pca_res, log = F, n_component = as.integer(max_components), 
                verbose = T, return_all = T))
            tmp <- do.call(UMAP, umap_args)
            tmp$embedding_ <- (tmp$embedding_ - min(tmp$embedding_))/max(tmp$embedding_)
            umap_res <- tmp$embedding_
            adj_mat <- Matrix::sparseMatrix(i = tmp$graph_@j, 
                p = tmp$graph_@p, x = -as.numeric(tmp$graph_@x), 
                dims = c(ncol(cds), ncol(cds)), index1 = F, dimnames = list(colnames(cds), 
                  colnames(cds)))
            S <- t(umap_res)
            Y <- S
            W <- t(irlba_pca_res)
            minSpanningTree(cds) <- graph_from_adjacency_matrix(adj_mat, 
                weighted = TRUE)
            A <- S
            colnames(A) <- colnames(FM)
            reducedDimA(cds) <- A
            colnames(S) <- colnames(FM)
            colnames(Y) <- colnames(FM)
            monocle:::reducedDimW(cds) <- W
            monocle:::reducedDimS(cds) <- as.matrix(Y)
            monocle:::reducedDimK(cds) <- S
            cds@dim_reduce_type <- reduction_method
        }
    }
    return(cds)
}
