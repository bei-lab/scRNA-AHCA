#library(Seurat)
#library(scales)
library(rlang)
#library(ggplot2)

Fixed_DoHeatmap <- function (object, features = NULL, cells = NULL, group.by = "ident", 
    group.bar = TRUE, group.colors = NULL, disp.min = -2.5, disp.max = NULL, 
    slot = "scale.data", assay = NULL, label = TRUE, size = 5.5, 
    hjust = 0, angle = 45, raster = TRUE, draw.lines = TRUE, 
    lines.width = NULL, group.bar.height = 0.02, combine = TRUE) 
{
    cells <- cells %||% colnames(x = object)
    if (is.numeric(x = cells)) {
        cells <- colnames(x = object)[cells]
    }
    assay <- assay %||% DefaultAssay(object = object)
    DefaultAssay(object = object) <- assay
    features <- features %||% VariableFeatures(object = object)
    features <- rev(x = unique(x = features))
    disp.max <- disp.max %||% ifelse(test = slot == "scale.data", 
        yes = 2.5, no = 6)
    possible.features <- rownames(x = GetAssayData(object = object, 
        slot = slot))
    if (any(!features %in% possible.features)) {
        bad.features <- features[!features %in% possible.features]
        features <- features[features %in% possible.features]
        if (length(x = features) == 0) {
            stop("No requested features found in the ", slot, 
                " slot for the ", assay, " assay.")
        }
        warning("The following features were omitted as they were not found in the ", 
            slot, " slot for the ", assay, " assay: ", paste(bad.features, 
                collapse = ", "))
    }
    data <- as.data.frame(x = as.matrix(x = t(x = GetAssayData(object = object, 
        slot = slot)[features, cells, drop = FALSE])))
    object <- suppressMessages(expr = StashIdent(object = object, 
        save.name = "ident"))
    group.by <- group.by %||% "ident"
    groups.use <- object[[group.by]][cells, , drop = FALSE]
    plots <- vector(mode = "list", length = ncol(x = groups.use))
    for (i in 1:ncol(x = groups.use)) {
        data.group <- data
        group.use <- groups.use[, i, drop = TRUE]
        if (!is.factor(x = group.use)) {
            group.use <- factor(x = group.use)
        }
        names(x = group.use) <- cells
        if (draw.lines) {
            lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) * 
                0.0025)
            placeholder.cells <- sapply(X = 1:(length(x = levels(x = group.use)) * 
                lines.width), FUN = function(x) {
                return(Seurat:::RandomName(length = 20))
            })
            placeholder.groups <- rep(x = levels(x = group.use), 
                times = lines.width)
            group.levels <- levels(x = group.use)
            names(x = placeholder.groups) <- placeholder.cells
            group.use <- as.vector(x = group.use)
            names(x = group.use) <- cells
            group.use <- factor(x = c(group.use, placeholder.groups), 
                levels = group.levels)
            na.data.group <- matrix(data = NA, nrow = length(x = placeholder.cells), 
                ncol = ncol(x = data.group), dimnames = list(placeholder.cells, 
                  colnames(x = data.group)))
            data.group <- rbind(data.group, na.data.group)
        }
        lgroup <- length(levels(group.use))
        plot <- Seurat:::SingleRasterMap(data = data.group, raster = raster, 
            disp.min = disp.min, disp.max = disp.max, feature.order = features, 
            cell.order = names(x = sort(x = group.use)), group.by = group.use, colors = colorRampPalette(c("navy", "white", "firebrick3"))(100))
        if (group.bar) {
            default.colors <- c(hue_pal()(length(x = levels(x = group.use))))
            cols <- group.colors[1:length(x = levels(x = group.use))] %||% 
                default.colors
            if (any(is.na(x = cols))) {
                cols[is.na(x = cols)] <- default.colors[is.na(x = cols)]
                cols <- Seurat:::Col2Hex(cols)
                col.dups <- sort(x = unique(x = which(x = duplicated(x = substr(x = cols, 
                  start = 1, stop = 7)))))
                through <- length(x = default.colors)
                while (length(x = col.dups) > 0) {
                  pal.max <- length(x = col.dups) + through
                  cols.extra <- hue_pal()(pal.max)[(through + 
                    1):pal.max]
                  cols[col.dups] <- cols.extra
                  col.dups <- sort(x = unique(x = which(x = duplicated(x = substr(x = cols, 
                    start = 1, stop = 7)))))
                }
            }
            group.use2 <- sort(x = group.use)
            if (draw.lines) {
                na.group <- Seurat:::RandomName(length = 20)
                levels(x = group.use2) <- c(levels(x = group.use2), 
                  na.group)
                group.use2[placeholder.cells] <- na.group
                cols <- c(cols, "#FFFFFF")
            }
            pbuild <- ggplot_build(plot = plot)
            names(x = cols) <- levels(x = group.use2)
            y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
            y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) + 
                y.range * 0.015
            y.max <- y.pos + group.bar.height * y.range
            plot <- plot + annotation_raster(raster = t(x = cols[group.use2]), 
                xmin = -Inf, xmax = Inf, ymin = y.pos, ymax = y.max) + 
                coord_cartesian(ylim = c(0, y.max), clip = "off") + 
                scale_color_manual(values = cols)
            if (label) {
                x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
                x.divs <- pbuild$layout$panel_params[[1]]$x.major %||% 
                  pbuild$layout$panel_params[[1]]$x$break_positions()
                x <- data.frame(group = sort(x = group.use), 
                  x = x.divs)
                label.x.pos <- tapply(X = x$x, INDEX = x$group, 
                  FUN = median) * x.max
                label.x.pos <- data.frame(group = names(x = label.x.pos), 
                  label.x.pos)
                plot <- plot + geom_text(stat = "identity", data = label.x.pos, 
                  aes_string(label = "group", x = "label.x.pos"), 
                  y = y.max + y.max * 0.03 * 0.5, angle = angle, 
                  hjust = hjust, size = size)
                plot <- suppressMessages(plot + coord_cartesian(ylim = c(0, 
                  y.max + y.max * 0.002 * max(nchar(x = levels(x = group.use)))* 
                    size), clip = "off"))
            }
        }
        plot <- plot + theme(line = element_blank())
        plots[[i]] <- plot
    }
    if (combine) {
        plots <- patchwork::wrap_plots(plots)
    }
    return(plots)
}
