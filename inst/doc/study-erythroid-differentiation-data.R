## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE
)

## ----setup--------------------------------------------------------------------
library(GenomicRanges)
library(smer)

## ----bim_data-----------------------------------------------------------------
bim_data <- data.frame(
  chromosome = c(1, 1, 1, 2, 2, 2, 3, 3, 3),
  variant_id = c("rs1", "rs2", "rs3", "rs4", "rs5", "rs6", "rs7", "rs8", "rs9"),
  cm_position = c(0, 0, 0, 0, 0, 0, 0, 0, 0),
  bp_position = c(10, 20, 30, 40, 50, 60, 70, 80, 90),
  allele1 = c("A", "A", "A", "G", "C", "C", "T", "T", "A"),
  allele1 = c("G", "G", "G", "A", "T", "T", "A", "A", "G")
)
bim_data$index <- 1:nrow(bim_data)

# DHS intervals
hg19_dhs_regions <- data.frame(
  chromosome = c(1, 2, 3),
  start = c(5, 45, 85),
  stop = c(15, 55, 95)
)

# LD block intervals
hg19_ld_blocks <- data.frame(
  chromosome = c(1, 1, 2, 2, 3, 3, 3),
  start = c(5, 25, 35, 45, 65, 75, 85),
  stop = c(25, 35, 45, 65, 75, 85, 95)
)

## ----dhs_data-----------------------------------------------------------------
# Convert .bim to GRanges object
bim_gr <- GRanges(
  seqnames = paste0("chr", bim_data$chromosome),
  ranges = IRanges(start = bim_data$bp_position, end = bim_data$bp_position),
  variant_id = bim_data$variant_id,
  genome = "hg19"
)

# Convert DHS to GRanges object
dhs_gr <- GRanges(
  seqnames = paste0("chr", hg19_dhs_regions$chromosome),
  ranges = IRanges(start = hg19_dhs_regions$start, end = hg19_dhs_regions$stop),
  genome = "hg19"
)

# Find overlaps of BIM variants and DHS intervals
overlaps <- findOverlaps(bim_gr, dhs_gr, maxgap = 0)

# Extract overlapping variants
dhs_data <- bim_data[queryHits(overlaps), ]
dhs_data <- dhs_data[!duplicated(dhs_data$index), ]

## -----------------------------------------------------------------------------
# Convert to GRanges object
ld_gr <- GRanges(
  seqnames = paste0("chr", hg19_ld_blocks$chromosome),
  ranges = IRanges(start = hg19_ld_blocks$start, end = hg19_ld_blocks$stop),
  genome = "hg19"
)

# Find LD block of bim variants
ld_overlaps <- findOverlaps(query = bim_gr, subject = ld_gr)

## ----write_mask---------------------------------------------------------------
output_file <- tempfile()
gxg_group <- "gxg"
ld_group <- "ld"

gxg_variants <- dhs_data$index - 1 # 0-base index for C++

create_hdf5_file(output_file)

for (j in bim_data$index - 1) { # 0-base index for C++
  # Write DHS mask
  gxg_ds <- sprintf("%s/%d", gxg_group, j)
  write_hdf5_dataset(file_name = output_file,
                     dataset_name = gxg_ds,
                     gxg_variants)

  # Find LD block of focal SNP
  focal_gr <- ld_gr[subjectHits(ld_overlaps[j,])]

  # Find variants in LD block of focal SNP
  focal_ld <- findOverlaps(query = bim_gr, subject = focal_gr)
  ld_data <- bim_data[queryHits(focal_ld),]
  ld_variants <- ld_data$index - 1 # 0-base index for C++

  # Write LD mask
  ld_ds <- sprintf("%s/%d", ld_group, j)
  write_hdf5_dataset(file_name = output_file,
                     dataset_name = ld_ds,
                     ld_variants)
}

dhs_indices <- read_hdf5_dataset(file_name = output_file, dataset_name = gxg_ds)
print(sprintf("DHS indices: %s", paste(dhs_indices, collapse = ", ")))

## ----run_sme, eval = FALSE----------------------------------------------------
# sme_result <- sme(
#   plink_file = "/path/to/plink/data",
#   pheno_file = "/path/to/pheno/data",
#   mask_file = "/path/to/mask/file",
#   gxg_indices = c(1, 2, 3),
#   chunk_size = 250,
#   n_randvecs = 10,
#   n_blocks = 200,
#   n_threads = 6
# )

## ----sessinfo-----------------------------------------------------------------
sessionInfo()

