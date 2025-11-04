library(tidyverse)
library(Seurat)
library(spacexr)
library(Matrix)
library(doParallel)
RhpcBLASctl::omp_set_num_threads(1)
RhpcBLASctl::blas_set_num_threads(1)
setwd('~/projects/spatial_nigra_v3/')

reference <- readRDS('../spatial_nigra_v2/resources/macosko_reference_obj.rds')
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

# Capture output -----------------------------------------------------------------------------------------------------------------
stdout <- vector('character')
con    <- textConnection('stdout', 'wr', local = TRUE)
sink(con)

rctd_list <- mapply(function(probes, barcodes, images) {
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
    
    coords <- GetTissueCoordinates(seurat) %>% dplyr::select(x, y)
    counts <- GetAssayData(seurat, layer = 'counts')
    puck <- SpatialRNA(coords = coords, counts = counts)
    
    return(puck)
    
}, probes = probefiles, barcodes = barcodefiles, images = imagefiles, SIMPLIFY = FALSE)
names(rctd_list) <- str_sub(names(rctd_list), start = 6L, end = -29L)
sink()
close(con)
write_lines(stdout, 'scripts/01_run_RCTD.log1')

## End capture output ------------------------------------------

#### run RCTD on replicactes ####
visium_metad <- read_tsv('resources/essential_metad_combined_cohorts.txt') %>% 
    mutate(case_new = ifelse(group %in% c('AD', 'ADLBD'), paste0('NBB', str_replace(case, '-', '_')), case))
group_ids <- visium_metad[match(names(rctd_list), visium_metad$case_new),] %>%
    mutate(group_ids = as.numeric(as.factor(group))) %>%
    dplyr::pull(group_ids, name = 'group')
myRCTD <- create.RCTD.replicates(spatialRNA.replicates = rctd_list,
                                 reference = reference,
                                 replicate_names = names(rctd_list),
                                 group_ids = group_ids,
                                 max_cores = 30)
saveRDS(myRCTD, 'resources/myRCTD.rds', compress = FALSE)
myRCTD <- run.RCTD.replicates(myRCTD, doublet_mode = 'full')
saveRDS(myRCTD, 'resources/myRCTD_reps_fullmode.rds')