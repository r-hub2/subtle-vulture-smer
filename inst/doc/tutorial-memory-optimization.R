## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE
)

## ----setup--------------------------------------------------------------------
library(smer)
library(tidyr)
library(dplyr)
library(knitr)

## ----memory-------------------------------------------------------------------
n_samples <- c(350000)
n_snps <- c(500000)
n_blocks <- c(1, 100, 1000)
n_randvecs <- c(10, 100)
chunk_size <- c(10, 100)

parameters <- crossing(
  n_samples = n_samples,
  n_snps = n_snps,
  n_blocks = n_blocks,
  n_randvecs = n_randvecs,
  chunk_size = chunk_size
)


estimated_memory <- parameters %>%
  mutate(memory_gb = round(
    approximate_memory_requirements(n_samples,
                                    n_snps,
                                    n_blocks,
                                    n_randvecs,
                                    chunk_size),
    2
  ))

kable(estimated_memory)

## ----seesionInfo--------------------------------------------------------------
sessionInfo()

