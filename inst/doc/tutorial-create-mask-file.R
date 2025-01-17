## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(smer)

## -----------------------------------------------------------------------------
hdf5_file <- tempfile()

# Group names
gxg_h5_group <- "gxg"
ld_h5_group <- "ld"

# Data (still in 1-based R indexing)
include_gxg_snps <- 1:10
exclude_ld_snps <- 5:6

# Focal SNP (still in 1-based R indexing)
focal_snp <- 4

# Dataset names
dataset_name_pattern <- "%s/%s"
# 0-based index!
gxg_dataset <- sprintf(dataset_name_pattern, gxg_h5_group, focal_snp - 1)
ld_dataset <- sprintf(dataset_name_pattern, ld_h5_group, focal_snp - 1)

# Create an empty HDF5 file
create_hdf5_file(hdf5_file)

# Write LD data
write_hdf5_dataset(hdf5_file, ld_dataset, exclude_ld_snps - 1) # 0-based index!

# Write GXG data
write_hdf5_dataset(hdf5_file, gxg_dataset, include_gxg_snps - 1)

## -----------------------------------------------------------------------------
ld_read <- read_hdf5_dataset(hdf5_file, ld_dataset)
gxg_read <- read_hdf5_dataset(hdf5_file, gxg_dataset)

print(sprintf("Zero-based indices of SNPs to exclude: %s", str(ld_read)))
print(sprintf("Zero-based indices of SNPs to include: %s", str(gxg_read)))

## ----seesionInfo--------------------------------------------------------------
sessionInfo()

