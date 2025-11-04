library(BASS)
library(Seurat)
library(tidyverse)
library(sctransform)
library(future)
source('../ggplot_theme_FLS.R')
plan("multisession", workers = 24)
plan()
set.seed(42)

##### Get metadata #####
metad <- read_tsv('resources/essential_metad_combined_cohorts.txt') %>% 
    mutate(case = case_when(
        group %in% c('AD', 'ADLBD') ~ paste0('NBB', str_replace(case, '-', '_')),
        TRUE ~ case
    ))
excludesamples <- read_lines('resources/excludesamples.txt')
filterprobes <- read_lines('resources/drop_probes.txt') %>% str_sub(., start = 17L)

##### Load BASS object and extract cluster assignments #####
mybass <- readRDS('resources/mybass_batchcorr.rds')
zlabels <- mybass@results$z # spatial domain labels 
znames <- lapply(mybass@xy, rownames)
zlabels_named <- mapply(function(zlabel, zname) {
    res <- data.frame(barcode = zname, cluster = zlabel)
    return(res)
}, zlabel = zlabels, zname = znames, SIMPLIFY = FALSE) 
names(zlabels_named) <- names(znames) 
zlabels_named <- data.table::rbindlist(zlabels_named, idcol = 'sample') %>% 
    mutate(cluster = paste0('c', cluster)) %>% 
    mutate(joincol = paste0(sample, '_', barcode))

##### Load Spaceranger output and make Seurat object list #####
probefiles <- list.files('data', pattern = 'raw_probe_bc_matrix.h5', recursive = TRUE, full.names = TRUE)
probefiles <- probefiles[grepl('outs', probefiles, ignore.case = FALSE)]
barcodefiles <- list.files('data', pattern = 'barcodes.tsv.gz', recursive = TRUE, full.names = TRUE)
barcodefiles <- barcodefiles[grepl('outs', barcodefiles, ignore.case = FALSE)]
barcodefiles <- barcodefiles[grepl('filtered_feature_bc_matrix', barcodefiles)]
filterprobes <- read_lines('resources/drop_probes.txt') %>% str_sub(., start = 17L)
imagefiles <- str_sub(list.files('data', pattern = 'spatial', recursive = TRUE, full.names = TRUE), end = -24L)
excludesamples <- read_lines('resources/excludesamples.txt')
imagefiles <- imagefiles[grepl('outs', imagefiles, ignore.case = FALSE)]

probefiles <- probefiles[!grepl(paste(excludesamples, collapse = '|'), probefiles)]
imagefiles <- imagefiles[!grepl(paste(excludesamples, collapse = '|'), imagefiles)]
imagefiles <- imagefiles[!grepl(paste(excludesamples, collapse = '|'), imagefiles)]
barcodefiles <- barcodefiles[!grepl(paste(excludesamples, collapse = '|'), barcodefiles)]
## sanity check:
test <- data.frame(image = imagefiles, probes = probefiles, barcodes = barcodefiles)

seurat_list <- mapply(function(probes, barcodes, images) {
    myid <- str_sub(images, start = 6L, end = -14L) # extract sample ID
    print(paste0('----> working on: ', myid))
    bc <- read_tsv(barcodes, col_names = FALSE) %>% pull(X1)
    print(paste0(myid, ': Filtering out ', length(bc), ' spots'))
    mat <- Read10X_h5(probes) 
    print(paste0(myid, ': dims before filtering :', paste0(dim(mat), collapse = ', ')))
    mat <- mat[!rownames(mat) %in% filterprobes,bc]
    mat <- mat[,colSums(mat) != 0]
    print(paste0(myid, ': dims before filtering :', paste0(dim(mat), collapse = ', ')))
    matched_probes <- str_sub(rownames(mat), end = -9L)
    mat_aggr <- list(Matrix.utils::aggregate.Matrix(mat, groupings = matched_probes, FUN = sum))
    assay.names <- 'spatial'
    slice.names <- myid
    image.list <- mapply(Read10X_Image, images, assay = assay.names, slice = slice.names, image.name = "tissue_hires_image.png")
    object.list <- mapply(CreateSeuratObject, mat_aggr, assay = assay.names)
    object.list <- mapply(function(.object, .image, .assay, 
                                   .slice) {
        .image <- .image[Cells(.object)]
        .object[[.slice]] <- .image
        return(.object)
    }, object.list, image.list, assay.names, slice.names)
    
    seurat <- merge(object.list[[1]], y = object.list[-1])
    
    return(seurat)
    
}, probes = probefiles, barcodes = barcodefiles, images = imagefiles, SIMPLIFY = FALSE)
names(seurat_list) <- str_sub(names(seurat_list), start = 6L, end = -29L)
saveRDS(seurat_list, 'resources/seurat_list_raw.rds', compress = FALSE)

##### Concatenate Seurat objects, normalize and add metadata #####
seurat <- merge(seurat_list[[1]], seurat_list[2:28], add.cell.ids = c(names(seurat_list[1]), names(seurat_list[2:28])), merge.data = TRUE, merge.dr = FALSE, project = 'spatial_nigra')
seurat$orig.ident <- str_sub(names(seurat$orig.ident), end = -20L)
seurat$group <- metad[match(seurat$orig.ident, metad$case),]$group
seurat$sex <- metad[match(seurat$orig.ident, metad$case),]$sex
seurat <- seurat %>% 
    NormalizeData(.) %>% 
    ScaleData(.)
table(colnames(seurat) %in% zlabels_named$joincol) # sanity check
seurat$cluster <- zlabels_named[match(colnames(seurat), zlabels_named$joincol),]$cluster

##### Name clusters upon histological examination #####
zlabels_named <- zlabels_named %>% 
    dplyr::mutate(cluster_anno = case_when(
        cluster == 'c1' ~ 'Fibers',
        cluster == 'c2' ~ 'SNtrans',
        cluster == 'c3' ~ 'SNpr',
        cluster == 'c4' ~ 'CerPed',
        cluster == 'c5' ~ 'SNpc'
    ))

##### Run MAGIC #####
reticulate::use_condaenv('~/miniforge3/bin/python3.10')
library(Rmagic)
seurat_magic <- magic(JoinLayers(seurat), n.jobs = 24, verbose = TRUE, seed = 42)
saveRDS(seurat_magic, 'resources/seurat_clustered_annotated_4excludesamples_MAGIC.rds', compress = FALSE)
seurat$cluster_anno <- zlabels_named[match(colnames(seurat), zlabels_named$joincol),]$cluster_anno


##### Get marker genes #####
DefaultAssay(seurat) <- 'spatial'
Idents(seurat) <- seurat$cluster_anno
allmarkers <- FindAllMarkers(seurat, assay = 'spatial', logfc.threshold = 0.5)
write_tsv(allmarkers, 'resources/seurat_spatial_allmarkers_4exludesamples.tsv')
