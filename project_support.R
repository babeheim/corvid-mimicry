
library(cmdstanr)
library(posterior)
library(rethinking)
library(dplyr)
library(readxl)
library(pander)
library(tinyplot)
library(testthat)
library(rotl)
library(ape)
library(MCMCglmm)
library(ggplot2)

panderOptions('table.split.table', Inf)  # This ensures the table is not split

source("R_functions/color_functions.R")
source("R_functions/misc_functions.R")
source("R_functions/estimation_functions.R")

data_path <- "raw_data/2025-02-08/mimicry.xlsx"

m_stan_single <- cmdstan_model("stan/single.stan")
m_stan_double <- cmdstan_model("stan/double.stan")
