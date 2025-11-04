library(tidyverse)
library(Seurat)
library(tidyseurat)
library(nicheDE)
library(spacexr)
RhpcBLASctl::omp_set_num_threads(1)
RhpcBLASctl::blas_set_num_threads(1)
mypwd <- c('/home/fstruebi/projects/spatial_nigra_v3/')

##### Get data #####
rctd <- readRDS(paste0(mypwd, 'resources/myRCTD_reps_fullmode.rds'))
deconv_list <- list()
for (i in 1:length(rctd@RCTD.reps)) {
    deconv_list[[i]] <- normalize_weights(rctd@RCTD.reps[[i]]@results$weights)
}
names(deconv_list) <- names(rctd@group_ids) 
avgexp_ref <- readRDS('~/projects/spatial_nigra_v2/resources/macosko_reference_matrix_for_nicheDE.rds')

##### Create Niche-DE object #####
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

seurat_list <- mapply(function(probes, barcodes, images) {
    myid <- str_sub(images, start = 6L, end = -14L) # extract sample ID
    print(paste0('----> working on: ', myid))
    bc <- read_tsv(barcodes, col_names = FALSE) %>% pull(X1)
    print(paste0(myid, ': Filtering out ', length(bc), ' spots'))
    mat <- Read10X_h5(probes) 
    print(paste0(myid, ': dims before filtering :', paste0(dim(mat), collapse = ', ')))
    mat <- mat[!rownames(mat) %in% filterprobes,bc]
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
saveRDS(names(seurat_list), 'resources/nicheDE_list_sigma-1-100-250_listnames.rds')

##### Filter everything for abundant cell types (pre-determined) #####
abndct <- read_tsv(paste0(mypwd, 'resources/top_spacexr_weights.tsv')) %>% 
    dplyr::filter(aggr_weight > 0.005)
cell_types_fil <- abndct$rowname
avgexp_ref_fil <- avgexp_ref[rownames(avgexp_ref) %in% cell_types_fil,]
deconv_list_fil <- list()
for (i in 1:length(deconv_list)) {
    deconv_list_fil[[i]] <- deconv_list[[i]][, colnames(deconv_list[[1]]) %in% cell_types_fil]
}

##### Run Niche-DE #####
stopifnot(names(seurat_list) == names(deconv_list_fil)) # sanity check
nde_list <- list()
for (i in 1:length(seurat_list)) {
    print(paste0('Running Niche-DE on list object ', i, ' of ', length(seurat_list), ' total items'))
    counts <- t(seurat_list[[i]]@assays$spatial$counts)
    counts_fil <- counts[rownames(counts) %in% rownames(deconv_list_fil[[i]]),]
    coords <- GetTissueCoordinates(seurat_list[[i]]) %>% dplyr::select(-cell) %>% as.matrix(.)
    coords_fil <- coords[rownames(coords) %in% rownames(deconv_list_fil[[i]]),]
    nde_list[[i]] <- CreateNicheDEObject(counts_mat = counts_fil, coordinate_mat = coords_fil, library_mat = avgexp_ref_fil, deconv_mat = deconv_list_fil[[i]], sigma = c(1, 100, 250), Int = TRUE)
    nde_list[[i]] <- CalculateEffectiveNiche(nde_list[[i]])
    nde_list[[i]] = niche_DE(nde_list[[i]],num_cores = 56, outfile = "", C = 150, M = 10, gamma = 0.7, print = TRUE, Int = TRUE, batch = TRUE, self_EN = TRUE, G = 10)
}
saveRDS(nde_list, paste0(mypwd, 'resources/nicheDE_list_sigma-1-100-250.rds'))