library(tidyverse)
library(Seurat)
library(edgeR)
library(enrichR)
library(ggpubr)
library(ggrepel)
library(sva)
library(patchwork)
library(data.table)
library(DEqMS)
library(broman)
library(ggrastr)
library(broom)

enrichdbs <- c('GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023')

remove_batch <- function(x, group, batch) {
    # x may be a sparse or dense matrix. It gets coerced into a dense matrix anyways, in order to function with ComBat
    sampleorder <- str_replace(str_sub(str_extract(colnames(x), pattern = '(_)(?!.*_).*'), start = 2L), '-', '_')
    braakvec <- metad[match(sampleorder, metad$case),]$braak_lewy_inferred
    mybatch <- batch
    # group is a variable whose coefficients to preserve for later testing
    cleandata <- ComBat_seq(as.matrix(x), mybatch, group = group, covar_mod = NULL)
    return(cleandata)
}

mdvec <- function(x, voi) {
    sampleorder <- str_replace(str_sub(str_extract(colnames(x), pattern = '(_)(?!.*_).*'), start = 2L), '-', '_')
    res <- metad[match(sampleorder, metad$case),]
    return(res[[voi]])
}

mdvec_prot <- function(x, voi) {
    colnames(x) <- ifelse(grepl('-', colnames(x)), paste0('NBB', str_replace(colnames(x), '-', '_')), colnames(x))
    sampleorder <- colnames(x)
    res <- metad[match(sampleorder, metad$case),]
    return(res[[voi]])
}

get_pseudobulk_matrix <- function(seurat, colgroups = c('cluster', 'group', 'orig.ident'), assay, layer) {
    pseudobulk <- as.matrix(PseudobulkExpression(seurat, group.by = colgroups, assay = assay, return.seurat = FALSE, layer = layer, method = 'aggregate', normalization.method = NULL)[[1]])
    pseudogroup <- str_sub(str_extract(colnames(pseudobulk), pattern = '[^_]*_[^_]*'), end = 2L)
    pb <- split.data.frame(t(pseudobulk), pseudogroup)
    pb <- lapply(pb, t)
    return(pb)
}

make_DGEList <- function(countmat, design, filter = TRUE, cname = 'NA') {
    y <- DGEList(counts=countmat)
    y <- normLibSizes(y, method="TMM")
    if(filter) {
        keep <- filterByExpr(y, design)
        cat(paste0("Cluster ", cname, ": Removing ", table(keep)['FALSE'], ' and keeping ', table(keep)['TRUE'], ' genes\n'))
        y <- y[keep, , keep.lib.sizes = FALSE]
    }
    return(y)
}

plot_MDS <- function(dgelist, is.pseudobulk = TRUE, color.by = 'group', label = FALSE, scale) {
    bulk_MDS <- plotMDS(dgelist, plot = FALSE)
    samplenames <- colnames(bulk_MDS$distance.matrix.squared)
    if (is.pseudobulk) {
        mymds_plot <- data.frame(x = bulk_MDS$x, y = bulk_MDS$y, samplenames = samplenames) %>% 
            separate(samplenames, into = c('cluster', 'group', 'case'), sep = '_') %>% 
            mutate(case = str_replace(case, '-', '_')) %>% 
            left_join(., metad, by = c('case', 'group'))
        pbyes <- paste0('pseudobulk, cluster ', mymds_plot$cluster)
    } else {
        mymds_plot <- data.frame(x = bulk_MDS$x, y = bulk_MDS$y, samplenames = samplenames) %>% 
            separate(samplenames, into = c('group', 'case'), sep = '_') %>% 
            mutate(case = str_replace(case, '-', '_')) %>% 
            left_join(., metad, by = c('case', 'group'))
        pbyes <- 'bulk'
    }
    doPlot <- function(cb) {
        pp <- mymds_plot %>% 
            ggplot(aes(x = x, y = y, color = .data[[cb]])) +
            geom_point(size = 4, alpha = 0.7) +
            theme_FLS() + 
            labs(x = 'MDS 1', y = 'MDS 2', title = paste0('MDS: ', pbyes))
        if(label) {
            pp <- pp + ggrepel::geom_text_repel(aes(label = case))
        }
        if(scale == 'continuous') {
            pp + scale_color_viridis_c()
        } else {
            pp + scale_color_discrete()
        }
        
    }
    if (length(color.by) == 1) {
        doPlot(cb = color.by)
    } else {
        lapply(color.by, doPlot)
    }
}

edgeR_DEpipe <- function(dgelist, design, contrasts = NULL, coef = NULL, alpha.fdr = 0.1, min.logFC = 0.5, cname = 'NA'){

    y <- estimateDisp(dgelist, design)
    fit <- glmQLFit(y, design)
    if(length(contrasts) == 0) {
        if(length(coef) == 0) {
            cat('Error: No coefficient specified\n')
        }
        qlf <- glmQLFTest(fit, coef=coef)
        de_res <- as.data.frame(topTags(qlf, n = Inf)) %>% 
            rownames_to_column('gene') %>% 
            mutate(significant = ifelse(FDR <= alpha.fdr, 'yes', 'no'))
    } else {
        if(ncol(contrasts) == 1) {
            qlf <- glmQLFTest(fit, contrast = contrasts[,1])
            de_res <- as.data.frame(topTags(qlf, n = Inf)) %>% 
                rownames_to_column('gene') %>% 
                mutate(significant = ifelse(FDR <= alpha.fdr, 'yes', 'no'))
        } else {
            qlf_list <- lapply(colnames(contrasts), function(x) {glmQLFTest(fit, contrast = contrasts[,x])})
            names(qlf_list) <- colnames(contrasts)
            de_res <- lapply(qlf_list, function(x) {as.data.frame(topTags(x, n = Inf)) %>% 
                    rownames_to_column('gene') %>% 
                    mutate(significant = ifelse(FDR <= alpha.fdr & abs(logFC) > min.logFC, 'yes', 'no'))}) %>% 
                rbindlist(idcol = 'test')
        }
    }
    cat(paste0('number of significant events for cluster ', cname, ': ', table(de_res$significant)['yes'], '\n'))
    return(de_res %>% arrange(FDR))
}

plot_DE <- function(edger_res, mytitle, n_label = 10, mult = c(0, 0.15)) {
    tops_fdr <- edger_res %>% 
        filter(significant == 'yes') %>% 
        slice_min(., order_by = FDR, n = n_label, with_ties = FALSE) %>% 
        pull(gene)
    tops_logcpm <- edger_res %>% 
        filter(significant == 'yes') %>% 
        slice_max(., order_by = abs(logCPM), n = n_label, with_ties = FALSE) %>% 
        pull(gene)
    edger_res <- edger_res %>% 
        mutate(label_fdr = ifelse(gene %in% tops_fdr, gene, ''),
               label_logcpm= ifelse(gene %in% tops_logcpm, gene, '')) %>% 
        mutate(significant = factor(case_when(significant == 'yes' & logFC > 0 ~ 'up',
                                              significant == 'yes' & logFC < 0 ~ 'down',
                                              TRUE ~ NA)))
    volcano <- edger_res %>% 
        ggplot(aes(y = -log10(PValue), x = logFC, color= significant, label = label_fdr)) +
        geom_point(alpha = 0.65) +
        geom_label_repel(show.legend = FALSE, 
                         max.overlaps = Inf, 
                         box.padding = 0.25, 
                         point.size = NA, 
                         # color = 'black',
                         min.segment.length = 0,
                         point.padding = 0,
                         label.padding = 0.25,
                         force = 10,
                         force_pull = 10) +
        scale_color_discrete(guide = 'none') +
        labs(title = 'Volcano Plot')
    ma <- edger_res %>% 
        ggplot(aes(x = logCPM, y = logFC, color= significant, label = label_logcpm)) +
        geom_point(alpha = 0.65) +
        geom_label_repel(show.legend = FALSE, 
                         max.overlaps = Inf, 
                         box.padding = 0.25, 
                         point.size = NA, 
                         # color = 'black',
                         min.segment.length = 0,
                         point.padding = 0,
                         label.padding = 0.25,
                         force = 10,
                         force_pull = 10) +
        labs(title = 'Mean-Average Plot')
    a <- patchwork::wrap_plots(volcano, ma) 
    a + plot_annotation(title = mytitle, theme = theme(plot.title = element_text(hjust = 0.5)))
}

plot_volcano <- function(edger_res, mytitle, n_label = 10, mult = c(0, 0.15)) {
    p_cutoff <- max(edger_res$PValue[ edger_res$FDR <= 0.1 ], na.rm=TRUE)
    tops_fdr1 <- edger_res %>% 
        ungroup() %>% 
        filter(significant == 'yes') 
    tops_fdr <- tops_fdr1 %>%
        slice_min(., order_by = FDR, n = n_label, with_ties = FALSE) %>%
        pull(gene)
    edger_res <- edger_res %>% 
        mutate(label_fdr = ifelse(gene %in% tops_fdr, gene, '')) %>% 
        mutate(significant = factor(case_when(significant == 'yes' & logFC > 0 ~ 'up',
                                              significant == 'yes' & logFC < 0 ~ 'down',
                                              TRUE ~ NA)))
    volcano <- edger_res %>%
        ggplot(aes(y = -log10(PValue), x = logFC, label = label_fdr)) +
        geom_point(alpha = 0.65, aes(color = significant)) +
        geom_text_repel(data = subset(edger_res, logFC < -0.5), 
                        box.padding = 0.5, 
                        max.overlaps = Inf, 
                        min.segment.length = 0, 
                        segment.size = 0.2,
                        direction = 'y', 
                        nudge_x = -1 + subset(edger_res, logFC < -0.5)$logFC, 
                        hjust = 0.5) +
        geom_text_repel(data = subset(edger_res, logFC > 0.5), 
                        box.padding = 0.5, 
                        max.overlaps = Inf, 
                        min.segment.length = 0, 
                        segment.size = 0.2,
                        direction = 'y', 
                        nudge_x = 1 + subset(edger_res, logFC > 0.5)$logFC, 
                        hjust = 0.5) +
        geom_hline(yintercept = -log10(p_cutoff), linetype = 'dashed', alpha = 0.5) +
        geom_vline(xintercept = 0.5, linetype = 'dashed', alpha = 0.5) +
        geom_vline(xintercept = -0.5, linetype = 'dashed', alpha = 0.5) +
        labs(title = mytitle) +
        scale_y_continuous(expand = expansion(mult = mult)) +
        scale_color_manual(guide = 'none', values = c('#F8766D', '#00BFC4')) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

plot_DE_prot <- function(deqms_res, mytitle, n_label = 10, mult = c(0, 0.15)) {
    p_cutoff <- max(deqms_res$PValue[ deqms_res$sca.adj.pval <= 0.1 ], na.rm=TRUE)
    tops_fdr1 <- deqms_res %>% 
        ungroup() %>% 
        filter(significant == 'yes') 
    tops_fdr <- tops_fdr1 %>%
        slice_min(., order_by = sca.adj.pval, n = n_label, with_ties = FALSE) %>%
        pull(gene)
    deqms_res <- deqms_res %>% 
        mutate(label_fdr = ifelse(gene %in% tops_fdr, gene, '')) %>% 
        mutate(significant = factor(case_when(significant == 'yes' & logFC > 0 ~ 'up',
                                              significant == 'yes' & logFC < 0 ~ 'down',
                                              TRUE ~ NA)))

    volcano <- deqms_res %>% 
        ggplot(aes(y = -log10(P.Value), x = logFC, label = label_fdr)) +
        geom_point(alpha = 0.65, aes(color= significant)) +
        geom_text_repel(data = subset(deqms_res, logFC < -0.5), 
                        box.padding = 0.5, 
                        max.overlaps = Inf, 
                        min.segment.length = 0, 
                        segment.size = 0.2,
                        direction = 'y', 
                        nudge_x = -0.5 + subset(deqms_res, logFC < -0.5)$logFC, 
                        hjust = 0.5) +
        geom_text_repel(data = subset(deqms_res, logFC > 0.5), 
                        box.padding = 0.5, 
                        max.overlaps = Inf, 
                        min.segment.length = 0, 
                        segment.size = 0.2,
                        direction = 'y', 
                        nudge_x = 0.5 + subset(deqms_res, logFC > 0.5)$logFC, 
                        hjust = 0.5) +

        geom_hline(yintercept = -log10(p_cutoff), linetype = 'dashed', alpha = 0.5) +
        geom_vline(xintercept = 0.5, linetype = 'dashed', alpha = 0.5) +
        geom_vline(xintercept = -0.5, linetype = 'dashed', alpha = 0.5) +
        labs(title = mytitle) +
        scale_y_continuous(expand = expansion(mult = mult)) +
        scale_color_manual(guide = 'none', values = c('#F8766D', '#00BFC4')) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

}

make_DE_heatmap <- function(de.df) {
    de_all_heatmap <- de.df %>% 
        group_by(cluster, test, significant) %>% 
        tally() %>% 
        pivot_wider(names_from = significant, values_from = n) %>% 
        dplyr::mutate(sigDE = ifelse(is.na(yes), 0, yes)) %>% 
        dplyr::select(cluster, test, sigDE) %>% 
        pivot_wider(names_from = 'cluster', values_from = 'sigDE') %>% 
        column_to_rownames('test') %>% 
        as.matrix(.)
    pheatmap(de_all_heatmap, cluster_cols = FALSE, cluster_rows = FALSE)
}

DEqMS_pipe <- function(intmat, pepmat, design, contrasts = NULL, coef = NULL, alpha.fdr = 0.1, cname = 'NA'){ 
    colnames(intmat) <- ifelse(grepl('-', colnames(intmat)), paste0('NBB', str_replace(colnames(intmat), '-', '_')), colnames(intmat))
    fit1 <- lmFit(intmat, design) 
    
    if(length(contrasts) == 0) {
        if(length(coef) == 0) {
            cat('Error: No coefficient specified\n')
        }
        fit3 <- eBayes(fit1, robust = TRUE)
        fit3$count <- pepmat[rownames(fit3$coefficients), ]
        fit4 <- spectraCounteBayes(fit3)
        mycoeffs <- colnames(fit4$coefficients)
        de_res <- outputResult(fit4, coef_col = coef) %>% 
            mutate(test = mycoeffs[coef])
    } else {
        fit2 <- contrasts.fit(fit1, contrasts = contrasts)
        fit3 <- eBayes(fit2, robust = TRUE)
        # topTable(fit3, number = Inf)
        fit3$count <- pepmat[rownames(fit3$coefficients), ]
        fit4 <- spectraCounteBayes(fit3)
        mycoeffs <- colnames(fit4$coefficients)
        de_res <- lapply(seq_along(mycoeffs), function(xx) {
            outputResult(fit4, coef_col = xx) %>% 
                mutate(test = mycoeffs[xx])
        }) %>% data.table::rbindlist(.) 
    }
    de_res <- de_res %>% 
        dplyr::select(test, gene, logFC, AveExpr, P.Value, adj.P.Val, min_psms = count, sca.adj.pval) %>% 
        dplyr::mutate(significant = case_when(
            is.nan(sca.adj.pval) & sca.adj.pval < alpha.fdr & abs(logFC) > 0.5 ~ 'yes',
            sca.adj.pval < alpha.fdr & abs(logFC) > 0.5 ~ 'yes',
            TRUE ~ 'no')) %>% 
        arrange(sca.adj.pval)
    cat(paste0('number of significant proteome events: ', table(de_res$significant)['yes'], '\n'))
    return(de_res)
}

### slightly modified after George Hall, UCL; https://github.com/george-hall-ucl/SpatialFeaturePlotBlend ###
SpatialFeaturePlotBlend <- function(object, features, combine = TRUE,
                                    images = Images(object),
                                    feature_1_alt_name = NULL,
                                    feature_2_alt_name = NULL, assay = NULL,
                                    bottom_left = "#000000",
                                    bottom_right = "#FF0000",
                                    top_left = "#00FF00",
                                    top_right = "#FFFF00", ...)  {
    
    # Generate a grid of RGB color values given the requested corner colours.
    gen_color_grid <- function(side_length, bottom_left, bottom_right,
                               top_left, top_right) {
        
        grad_gen <- function(start, end, n = side_length) {
            colfunc <- colorRampPalette(c(start, end))
            return(colfunc(n))
        }
        
        # x_y = "x to y"; "bl" = "bottom left", etc
        bl_tl <- grad_gen(bottom_left, bottom_right)
        br_tr <- grad_gen(top_left, top_right)
        
        l <- lapply(seq_len(length(bl_tl)),
                    function(i) {
                        start <- bl_tl[i]
                        end <- br_tr[i]
                        new_grad <- grad_gen(start, end)
                    })
        
        return(t(matrix(unlist(l), ncol = side_length, nrow = side_length)))
    }
    
    if (length(features) != 2) {
        stop(paste(c("Incorrect number of features. ",
                     "Requires two features, received ",
                     length(features))))
    }
    
    if (!is.null(assay)) {
        DefaultAssay(object) <- assay
    }
    
    blend_plot_theme <- theme(legend.position = "none",
                              plot.title = element_text(hjust = 0.5))
    
    plot_list_outer <- list()
    dat_norm_list <- list()
    
    for (i in images) {
        plot_list <- lapply(features,
                            function(feature) {
                                max_color <- ifelse(feature == features[1],
                                                    bottom_right, top_left)
                                SpatialFeaturePlot(object, feature,
                                                   images = i, ...) +
                                    scale_fill_gradient(low = bottom_left,
                                                        high = max_color) +
                                    ggtitle(feature) +
                                    blend_plot_theme
                            })
        
        cell_barcodes <- Seurat:::CellsByImage(object, images = i,
                                               unlist = TRUE)
        cells_obj_sub <- object[, cell_barcodes]
        images_sub_list <- list(object[[i]])
        names(images_sub_list) <- i
        cells_obj_sub@images <- images_sub_list
        dat <- FetchData(cells_obj_sub, features)
        side_length <- 100
        col_grid <- gen_color_grid(side_length, bottom_left, bottom_right,
                                   top_left, top_right)
        dat_norm <- apply(dat, 2,
                          function(x) {
                              round((side_length - 1) * x / max(x)) + 1
                          })
        colors <- sapply(seq_len(nrow(dat_norm)),
                         function(x) {
                             col_grid[dat_norm[x, 1], dat_norm[x, 2]]
                         })
        
        new_md_column <- paste0(features[1], "_vs_", features[2])
        cells_obj_sub[[new_md_column]] <- colors
        names(colors) <- as.character(colors)
        
        plot_list[[3]] <- SpatialDimPlot(cells_obj_sub, new_md_column,
                                         cols = colors, images = i, ...) +
            ggtitle(paste0(features[1], "_", features[2])) +
            blend_plot_theme
        
        legend_grid <- expand.grid(seq(from = min(dat[, features[1]]),
                                       to = max(dat[, features[1]]),
                                       length.out = side_length),
                                   seq(from = min(dat[, features[2]]),
                                       to = max(dat[, features[2]]),
                                       length.out = side_length))
        colnames(legend_grid) <- features
        legend_colors <- c(col_grid)
        legend_grid$color <- legend_colors
        names(legend_colors) <- legend_colors
        
        legend <- ggplot(legend_grid,
                         aes(x = .data[[features[1]]], y = .data[[features[2]]],
                             color = color)) +
            geom_point(shape = 15, size = 1.9) +
            scale_color_manual(values = legend_colors) +
            coord_cartesian(expand = FALSE) +
            theme(legend.position = "none", aspect.ratio = 1,
                  panel.background = element_blank(),
                  axis.text.x = element_text(angle = 45, hjust = 1)) +
            xlab(ifelse(is.null(feature_1_alt_name),
                        features[1], feature_1_alt_name)) +
            ylab(ifelse(is.null(feature_2_alt_name),
                        features[2], feature_2_alt_name))
        
        plot_list[[4]] <- wrap_plots(ggplot() + theme_void(), legend,
                                     ggplot() + theme_void(), ncol = 1,
                                     heights = c(0.2, 0.6, 0.2))
        
        plot_list_outer[[i]] <- plot_list
        dat_norm_list[[i]] <- dat_norm
    }
    
    if (combine == FALSE) {
        return(list(plot_list_outer, dat_norm_list))
    } else {
        plot_list_outer <- lapply(plot_list_outer,
                                  function(p) {
                                      wrap_plots(p, nrow = 1,
                                                 widths = c(0.28, 0.28,
                                                            0.28, 0.16))
                                  })
        p <- wrap_plots(plot_list_outer, ncol = 1)
        
        return(list(p, dat_norm_list))
    }
}

group_colors <- data.frame(
    group = c('Ctrl', 'iLBD', 'PD', 'AD', 'ADLBP'),
    color = c('#e0ecf4', '#9ebcda', '#8856a7', '#fdae6b', '#e6550d')
    # brocolors('crayons')['Eggplant']
)

reorder_groups <- function(x) {
    temp <- x %>% mutate(group = factor(group, levels = group_colors$group))
    if(NA %in% temp$group) {
        x %>% 
            mutate(group = ifelse(group == 'ADLBD', 'ADLBP', group)) %>% 
            mutate(group = factor(group, levels = group_colors$group))
    } else {
        temp
    }
}

