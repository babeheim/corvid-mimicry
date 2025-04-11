
rm(list = ls())

source("project_support.R")

cat("test prior\n")
# 20%, 30%, 50%, 70%

# our prior on q:

q_mu = 0.2
q_theta = 4

spp <- read.csv("data/corvid_species.csv")

# first, use every positive evidence of mimicry

dat_double <- list(
  N = nrow(spp),
  R1 = spp$primary_count,
  D1 = na_to_999(spp$mimicry_primary),
  R2 = spp$secondary_count,
  D2 = na_to_999(spp$mimicry_secondary),
  q_mu = q_mu, # vary this as robustness check
  q_theta = q_theta
)

# important check
stopifnot(all.equal(dat_double$D2 == 999, dat_double$R2 == 0))
stopifnot(all.equal(dat_double$D1 == 999, dat_double$R1 == 0))

# we can fit just the secondary lit like so:

fit_stan_double <- m_stan_double$sample(
  data = dat_double, 
  seed = 123, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500
)

fit_stan_double$draws() |> as_draws_rvars() -> post_double_2


q_mu = 0.3
q_theta = 4

spp <- read.csv("data/corvid_species.csv")

# first, use every positive evidence of mimicry

dat_double <- list(
  N = nrow(spp),
  R1 = spp$primary_count,
  D1 = na_to_999(spp$mimicry_primary),
  R2 = spp$secondary_count,
  D2 = na_to_999(spp$mimicry_secondary),
  q_mu = q_mu, # vary this as robustness check
  q_theta = q_theta
)

# important check
stopifnot(all.equal(dat_double$D2 == 999, dat_double$R2 == 0))
stopifnot(all.equal(dat_double$D1 == 999, dat_double$R1 == 0))

# we can fit just the secondary lit like so:

fit_stan_double <- m_stan_double$sample(
  data = dat_double, 
  seed = 123, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500
)

fit_stan_double$draws() |> as_draws_rvars() -> post_double_3


q_mu = 0.5
q_theta = 4

spp <- read.csv("data/corvid_species.csv")

# first, use every positive evidence of mimicry

dat_double <- list(
  N = nrow(spp),
  R1 = spp$primary_count,
  D1 = na_to_999(spp$mimicry_primary),
  R2 = spp$secondary_count,
  D2 = na_to_999(spp$mimicry_secondary),
  q_mu = q_mu, # vary this as robustness check
  q_theta = q_theta
)

# important check
stopifnot(all.equal(dat_double$D2 == 999, dat_double$R2 == 0))
stopifnot(all.equal(dat_double$D1 == 999, dat_double$R1 == 0))

# we can fit just the secondary lit like so:

fit_stan_double <- m_stan_double$sample(
  data = dat_double, 
  seed = 123, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500
)

fit_stan_double$draws() |> as_draws_rvars() -> post_double_5


q_mu = 0.7
q_theta = 4

spp <- read.csv("data/corvid_species.csv")

# first, use every positive evidence of mimicry

dat_double <- list(
  N = nrow(spp),
  R1 = spp$primary_count,
  D1 = na_to_999(spp$mimicry_primary),
  R2 = spp$secondary_count,
  D2 = na_to_999(spp$mimicry_secondary),
  q_mu = q_mu, # vary this as robustness check
  q_theta = q_theta
)

# important check
stopifnot(all.equal(dat_double$D2 == 999, dat_double$R2 == 0))
stopifnot(all.equal(dat_double$D1 == 999, dat_double$R1 == 0))

# we can fit just the secondary lit like so:

fit_stan_double <- m_stan_double$sample(
  data = dat_double, 
  seed = 123, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500
)

fit_stan_double$draws() |> as_draws_rvars() -> post_double_7



png("figures/estimates_prior_check.png", res = 300, units = "in", height = 5, width = 7)

plot(NULL, xlim = c(0, 1), ylim = c(0.5, 4.5), yaxt = "n", ylab = "",
  xlab = "occurrence of mimicry, q", xaxs = "i", yaxs = "i", main = "posterior estimates of q with various prior μ")

abline(v = seq(0, 1, by = 0.1), col = gray(0.5, 0.4))

q_double <- draws_of(post_double_7$q)
q_double_mu <- mean(q_double)
q_double_lb <- HPDI(q_double)[1]
q_double_ub <- HPDI(q_double)[2]

lines(c(q_double_lb, q_double_ub), c(4, 4))
points(q_double_mu, 4, pch = 20)
text(q_double_mu, 4, pos = 3, "μ = 0.7")

q_double <- draws_of(post_double_5$q)
q_double_mu <- mean(q_double)
q_double_lb <- HPDI(q_double)[1]
q_double_ub <- HPDI(q_double)[2]

lines(c(q_double_lb, q_double_ub), c(3, 3))
points(q_double_mu, 3, pch = 20)
text(q_double_mu, 3, pos = 3, "μ = 0.5")

q_double <- draws_of(post_double_3$q)
q_double_mu <- mean(q_double)
q_double_lb <- HPDI(q_double)[1]
q_double_ub <- HPDI(q_double)[2]

lines(c(q_double_lb, q_double_ub), c(2, 2))
points(q_double_mu, 2, pch = 20)
text(q_double_mu, 2, pos = 3, "μ = 0.3")

q_double <- draws_of(post_double_2$q)
q_double_mu <- mean(q_double)
q_double_lb <- HPDI(q_double)[1]
q_double_ub <- HPDI(q_double)[2]

lines(c(q_double_lb, q_double_ub), c(1, 1))
points(q_double_mu, 1, pch = 20)
text(q_double_mu, 1, pos = 3, "μ = 0.2")

dev.off()



cat("test reliability rating\n")

# our prior on q:
q_mu = 0.5
q_theta = 4

spp <- read.csv("data/corvid_species.csv")

# first, use every positive evidence of mimicry

n_full <- sum(spp$mimicry_primary == 1 | spp$mimicry_secondary == 1, na.rm = TRUE)

dat_double <- list(
  N = nrow(spp),
  R1 = spp$primary_count,
  D1 = na_to_999(spp$mimicry_primary),
  R2 = spp$secondary_count,
  D2 = na_to_999(spp$mimicry_secondary),
  q_mu = q_mu, # vary this as robustness check
  q_theta = q_theta
)

# important check
stopifnot(all.equal(dat_double$D2 == 999, dat_double$R2 == 0))
stopifnot(all.equal(dat_double$D1 == 999, dat_double$R1 == 0))

# we can fit just the secondary lit like so:

fit_stan_double <- m_stan_double$sample(
  data = dat_double, 
  seed = 123, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500
)

fit_stan_double$draws() |> as_draws_rvars() -> post_double_full

# now exclude low-reliability ratings

spp$mimicry_primary[which(spp$highest_rating_xc == 1)] <- 0
spp$mimicry_secondary[which(spp$highest_rating_secondary == 1)] <- 0

n_med_high <- sum(spp$mimicry_primary == 1 | spp$mimicry_secondary == 1, na.rm = TRUE)

dat_double <- list(
  N = nrow(spp),
  R1 = spp$primary_count,
  D1 = na_to_999(spp$mimicry_primary),
  R2 = spp$secondary_count,
  D2 = na_to_999(spp$mimicry_secondary),
  q_mu = q_mu, # vary this as robustness check
  q_theta = q_theta
)

# important check
stopifnot(all.equal(dat_double$D2 == 999, dat_double$R2 == 0))
stopifnot(all.equal(dat_double$D1 == 999, dat_double$R1 == 0))

# we can fit just the secondary lit like so:

fit_stan_double <- m_stan_double$sample(
  data = dat_double, 
  seed = 123, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500
)

fit_stan_double$draws() |> as_draws_rvars() -> post_double_med_high

# now exclude medium-reliability ratings

spp$mimicry_primary[which(spp$highest_rating_xc == 2)] <- 0
spp$mimicry_secondary[which(spp$highest_rating_secondary == 2)] <- 0

n_high <- sum(spp$mimicry_primary == 1 | spp$mimicry_secondary == 1, na.rm = TRUE)

dat_double <- list(
  N = nrow(spp),
  R1 = spp$primary_count,
  D1 = na_to_999(spp$mimicry_primary),
  R2 = spp$secondary_count,
  D2 = na_to_999(spp$mimicry_secondary),
  q_mu = q_mu, # vary this as robustness check
  q_theta = q_theta
)

# important check
stopifnot(all.equal(dat_double$D2 == 999, dat_double$R2 == 0))
stopifnot(all.equal(dat_double$D1 == 999, dat_double$R1 == 0))

# we can fit just the secondary lit like so:

fit_stan_double <- m_stan_double$sample(
  data = dat_double, 
  seed = 123, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500
)

fit_stan_double$draws() |> as_draws_rvars() -> post_double_high




png("figures/estimates_reliability_check.png", res = 300, units = "in", height = 4, width = 7)

plot(NULL, xlim = c(0, 1), ylim = c(0.5, 3.5), yaxt = "n", ylab = "",
  xlab = "occurrence of mimicry, q", xaxs = "i", yaxs = "i")

abline(v = seq(0, 1, by = 0.1), col = gray(0.5, 0.4))

q_double <- draws_of(post_double_high$q)
q_double_mu <- mean(q_double)
q_double_lb <- HPDI(q_double)[1]
q_double_ub <- HPDI(q_double)[2]

lines(c(q_double_lb, q_double_ub), c(3, 3))
points(q_double_mu, 3, pch = 20)
text(q_double_mu, 3, pos = 3, paste0("excluding 'low' and 'medium' (n = ", n_high, ")"))

q_double <- draws_of(post_double_med_high$q)
q_double_mu <- mean(q_double)
q_double_lb <- HPDI(q_double)[1]
q_double_ub <- HPDI(q_double)[2]

lines(c(q_double_lb, q_double_ub), c(2, 2))
points(q_double_mu, 2, pch = 20)
text(q_double_mu, 2, pos = 3, paste0("excluding 'low' (n = ", n_med_high, ")"))

q_double <- draws_of(post_double_full$q)
q_double_mu <- mean(q_double)
q_double_lb <- HPDI(q_double)[1]
q_double_ub <- HPDI(q_double)[2]

lines(c(q_double_lb, q_double_ub), c(1, 1))
points(q_double_mu, 1, pch = 20)
text(q_double_mu, 1, pos = 3, paste0("all detections (n = ", n_full, ")"))

dev.off()

