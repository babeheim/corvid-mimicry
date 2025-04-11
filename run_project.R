
rm(list = ls())

stan_files <- list.files("stan", full.names = TRUE)
stan_binaries <- stan_files[!grepl("\\.stan$", stan_files)]
if (length(stan_binaries) > 0) file.remove(stan_binaries)

source("project_support.R")

dir_init("./figures")

source("0_init_projects.R")

source("3_explore_data.R")

source("4_fit_stan_models.R")

source("5_fit_glmms.R")

source("6_test_robustness.R")
