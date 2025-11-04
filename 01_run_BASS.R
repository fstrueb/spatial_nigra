library(BASS)
library(Seurat)
library(tidyverse)

metad <- read_tsv('resources/essential_metad_combined_cohorts.txt') %>% 
    mutate(case = case_when(
        group %in% c('AD', 'ADLBD') ~ paste0('NBB', str_replace(case, '-', '_')),
        TRUE ~ case
    ))
excludesamples <- read_lines('resources/excludesamples.txt')
filterprobes <- read_lines('resources/drop_probes.txt') %>% str_sub(., start = 17L)

spatcoords_files <- list.files('data', recursive = TRUE, full.names = TRUE, pattern = 'tissue_positions.csv')
spatcoords_files <- spatcoords_files[grepl('outs', spatcoords_files, ignore.case = FALSE)]
spatcoords_files <- spatcoords_files[!grepl(paste0(excludesamples, collapse = '|'), spatcoords_files) ]
spatcoords <- lapply(spatcoords_files, function(x) {
    read_csv(x, col_types = cols()) %>%
        dplyr::filter(in_tissue == 1) %>%
        dplyr::select(barcode, array_row, array_col) %>%
        arrange(barcode) %>%
        column_to_rownames('barcode') %>%
        as.matrix(.)
})
names(spatcoords) <- str_sub(spatcoords_files, start = 6L, end = -35L)
spatcoords$PDC069 <- spatcoords$PDC069[rownames(spatcoords$PDC069) != 'CACGGTCTCCTTACGA-1',] # remove single empty spot from list

probefiles <- list.files('data', pattern = 'raw_probe_bc_matrix.h5', recursive = TRUE, full.names = TRUE)
probefiles <- probefiles[grepl('outs', probefiles, ignore.case = FALSE)]
probefiles <- probefiles[!grepl(paste0(excludesamples, collapse = '|'), probefiles) ]
barcodefiles <- list.files('data', pattern = 'barcodes.tsv.gz', recursive = TRUE, full.names = TRUE)
barcodefiles <- barcodefiles[grepl('outs', barcodefiles, ignore.case = FALSE)]
barcodefiles <- barcodefiles[!grepl(paste0(excludesamples, collapse = '|'), barcodefiles) ]
barcodefiles <- barcodefiles[grepl('filtered_feature_bc_matrix', barcodefiles)]

# Capture output -----------------------------------------------------------------------------------------------------------------
stdout <- vector('character')
con    <- textConnection('stdout', 'wr', local = TRUE)
sink(con)

spatmats <- mapply(function(probes, barcodes) {
    bc <- read_tsv(barcodes, col_names = FALSE) %>% pull(X1)
    print(paste0('Filtering out ', length(bc), ' probes'))
    mat <- Read10X_h5(probes)
    mat <- mat[!rownames(mat) %in% filterprobes,bc]
    print(dim(mat))
    mat <- mat[,colSums(mat) != 0]
    print(dim(mat))
    matched_probes <- str_sub(rownames(mat), end = -9L)
    mat_aggr <- Matrix.utils::aggregate.Matrix(mat, groupings = matched_probes, FUN = sum)
    # ## sanity check if required:
    # # mat[grepl('A2M', rownames(mat)),]
    # # mat_aggr['A2M',]
    
    return(mat_aggr)
    
}, probes = probefiles, barcodes = barcodefiles, SIMPLIFY = FALSE)
names(spatmats) <- str_sub(names(spatmats), start = 6L, end = -29L)

all(names(spatmats) == names(spatcoords)) #sanity check

sink()
close(con)
write_lines(stdout, 'scripts/03_run_BASS.log1')

## End capture output ------------------------------------------


# Setup BASS object:
mybass <- createBASSObject(spatmats, spatcoords, C = 20, R = 5, beta_method = 'SW', init_method = 'mclust', nsample = 10000)

# Data pre-processing:
# 1.Library size normalization followed with a log2 transformation
# 2.Select top 3000 spatially expressed genes with SPARK-X
# 3.Dimension reduction with PCA
mybass <- BASS.preprocess(mybass, doLogNormalize = TRUE,
                          geneSelect = "sparkx", nSE = 3000, doPCA = TRUE, doBatchCorrect = TRUE,
                          scaleFeature = FALSE, nPC = 20)

mybass <- BASS.run(mybass)
mybass <- BASS.postprocess(mybass)
saveRDS(mybass, 'resources/mybass_batchcorr.rds', compress = FALSE)
