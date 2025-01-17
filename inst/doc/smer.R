## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE
)

## ----setup--------------------------------------------------------------------
library(smer)
library(dplyr)
library(ggplot2)

## ----run_sme, eval = FALSE----------------------------------------------------
# # File inputs
# plink_file <- "path/to/plink/file"
# pheno_file <- "path/to/pheno/file"
# mask_file <- "path/to/mask/file"
# 
# # Parameter inputs
# chun_ksize <- 10
# n_randvecs <- 10
# n_blocks <- 10
# n_threads <- 5
# 
# # 1-based Indices of SNPs to be analyzed
# n_snps <- 100
# snp_indices <- 1:n_snps
# 
# sme_result <- sme(
#   plink_file,
#   pheno_file,
#   mask_file,
#   snp_indices,
#   chunk_size,
#   n_randvecs,
#   n_blocks,
#   n_threads
# )
# 

## ----assign_data, include = FALSE---------------------------------------------
sme_result <- getting_started

## ----manhattan_plot, fig.alt="Manhattan plot to illustrate the sparse marginal epistasis test"----
sme_result$summary %>%
  ggplot(aes(
  x = index,
  y = -log10(p),
  color = true_gxg_snp
)) +
  geom_point() +
  xlab("Position") +
  labs(color = "Epistatic SNP")

## ----pve_plot, fig.alt="PVE plot to illustrate the sparse marginal epistasis test"----
sme_result$summary %>%
  ggplot(aes(x = true_gxg_snp, y = pve, fill = true_gxg_snp)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.05, color = "grey40", linetype = "dashed") +
  annotate("text", x = 0.8, y =  0.055,
           label = "True per SNP epistatic PVE", color = "black") +
  xlab("Epistatic SNP") +
  ylab("Phenotypic Variance Explained") +
  theme(legend.position = "none")

## ----h2_plot, fig.alt="h2 plot to illustrate the sparse marginal epistasis test"----
sme_result$vc_estimate %>%
  ggplot(aes(x = component, y = vc_estimate, fill = component)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.3, color = "grey40", linetype = "dashed") +
  annotate("text", x = 0.7, y =  0.33,
           label = expression("True " * h^2), color = "black")  +
  xlab("Component") +
  ylab("Variance Component Estimate") +
  theme(legend.position = "none")

## ----seesionInfo--------------------------------------------------------------
sessionInfo()

