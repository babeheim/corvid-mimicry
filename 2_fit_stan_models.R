
rm(list = ls())

source("project_support.R")

spp <- read.csv("data/corvid_species.csv")

# our prior on q:
q_mu = 0.5
q_theta = 4

# first, let's fit the single-level models

dat_single1 <- list(
  N = nrow(spp),
  R = spp$primary_count,
  D = na_to_999(spp$mimicry_primary),
  q_mu = q_mu,
  q_theta = q_theta
)

fit_stan_single1 <- m_stan_single$sample(
  data = dat_single1, 
  seed = 123, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500
)

fit_stan_single1$draws() |> as_draws_rvars() -> post_single1

post_single1$rho
post_single1$q

samples <- as_draws_rvars(fit_stan_single1$draws())

samples |>
  summarize_draws_rvars(format = FALSE) |>
  bind_rows() |>
  select(parameter = label, everything()) -> out

  tar <- which(out$parameter %in% "lp__")
  out$mean[tar] <- sprintf("%2.2f", as.numeric(out$mean[tar]))
  out$sd[tar] <- sprintf("%2.2f", as.numeric(out$sd[tar]))
  out$lb[tar] <- sprintf("%2.2f", as.numeric(out$lb[tar]))
  out$ub[tar] <- sprintf("%2.2f", as.numeric(out$ub[tar]))
  tar <- which(!(out$parameter %in% "lp__"))
  out$mean[tar] <- sprintf("%2.4f", as.numeric(out$mean[tar]))
  out$sd[tar] <- sprintf("%2.4f", as.numeric(out$sd[tar]))
  out$lb[tar] <- sprintf("%2.4f", as.numeric(out$lb[tar]))
  out$ub[tar] <- sprintf("%2.4f", as.numeric(out$ub[tar]))

write_pandoc_table(out, "figures/single1_estimates.md", style = "rmarkdown")

dat_single2 <- list(
  N = nrow(spp),
  R = spp$secondary_count,
  D = na_to_999(spp$mimicry_secondary),
  q_mu = q_mu,
  q_theta = q_theta
)

fit_stan_single2 <- m_stan_single$sample(
  data = dat_single2, 
  seed = 123, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500
)

fit_stan_single2$draws() |> as_draws_rvars() -> post_single2

post_single2$rho # per-recording rate is extremely low
post_single2$q

samples <- as_draws_rvars(fit_stan_single2$draws())

samples |>
  summarize_draws_rvars(format = FALSE) |>
  bind_rows() |>
  select(parameter = label, everything()) -> out

  tar <- which(out$parameter %in% "lp__")
  out$mean[tar] <- sprintf("%2.2f", as.numeric(out$mean[tar]))
  out$sd[tar] <- sprintf("%2.2f", as.numeric(out$sd[tar]))
  out$lb[tar] <- sprintf("%2.2f", as.numeric(out$lb[tar]))
  out$ub[tar] <- sprintf("%2.2f", as.numeric(out$ub[tar]))
  tar <- which(!(out$parameter %in% "lp__"))
  out$mean[tar] <- sprintf("%2.4f", as.numeric(out$mean[tar]))
  out$sd[tar] <- sprintf("%2.4f", as.numeric(out$sd[tar]))
  out$lb[tar] <- sprintf("%2.4f", as.numeric(out$lb[tar]))
  out$ub[tar] <- sprintf("%2.4f", as.numeric(out$ub[tar]))

write_pandoc_table(out, "figures/single2_estimates.md", style = "rmarkdown")


# now fit the two-level model...

dat_double <- list(
  N = nrow(spp),
  R1 = spp$primary_count,
  D1 = na_to_999(spp$mimicry_primary),
  R2 = spp$secondary_count,
  D2 = na_to_999(spp$mimicry_secondary),
  q_mu = q_mu,
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

fit_stan_double$draws() |> as_draws_rvars() -> post_double

post_double$rho1
post_double$rho2
post_double$q

samples <- as_draws_rvars(fit_stan_double$draws())

samples |>
  summarize_draws_rvars(format = FALSE) |>
  bind_rows() |>
  select(parameter = label, everything()) -> out

  tar <- which(out$parameter %in% "lp__")
  out$mean[tar] <- sprintf("%2.2f", as.numeric(out$mean[tar]))
  out$sd[tar] <- sprintf("%2.2f", as.numeric(out$sd[tar]))
  out$lb[tar] <- sprintf("%2.2f", as.numeric(out$lb[tar]))
  out$ub[tar] <- sprintf("%2.2f", as.numeric(out$ub[tar]))
  tar <- which(!(out$parameter %in% "lp__"))
  out$mean[tar] <- sprintf("%2.4f", as.numeric(out$mean[tar]))
  out$sd[tar] <- sprintf("%2.4f", as.numeric(out$sd[tar]))
  out$lb[tar] <- sprintf("%2.4f", as.numeric(out$lb[tar]))
  out$ub[tar] <- sprintf("%2.4f", as.numeric(out$ub[tar]))

write_pandoc_table(out, "figures/double_estimates.md", style = "rmarkdown")


pr_mimic_single1 <- matrix(NA, ncol = nrow(spp), nrow = 4000)
for (i in 1:ncol(pr_mimic_single1)) {
  if (dat_single1$D[i] == 1) {
    pr_mimic_single1[, i] <- rep(1, 4000)
  } else {
    pr_mimic_single1[, i] <- draws_of(pr_hidden_single(dat_single1$R[i], post_single1))
  }
}

pr_mimic_mu_single1 <- apply(pr_mimic_single1, 2, mean)



pr_mimic_single2 <- matrix(NA, ncol = nrow(spp), nrow = 4000)
for (i in 1:ncol(pr_mimic_single2)) {
  if (dat_single2$D[i] == 1) {
    pr_mimic_single2[, i] <- rep(1, 4000)
  } else {
    pr_mimic_single2[, i] <- draws_of(pr_hidden_single(dat_single2$R[i], post_single2))
  }
}

pr_mimic_mu_single2 <- apply(pr_mimic_single2, 2, mean)


pr_mimic_double <- matrix(NA, ncol = nrow(spp), nrow = 4000)
for (i in 1:ncol(pr_mimic_double)) {
  if (dat_double$D1[i] == 1 | dat_double$D2[i] == 1) {
    pr_mimic_double[, i] <- rep(1, 4000)
  } else {
    pr_mimic_double[, i] <- draws_of(pr_hidden_double(dat_double$R1[i], dat_double$R2[i], post_double))
  }
}

pr_mimic_mu_post <- apply(pr_mimic_double, 2, mean)
pr_mimic_se_post <- apply(pr_mimic_double, 2, sd)



png("figures/estimates_emp.png", res = 300, units = "in", height = 4, width = 7)

plot(NULL, xlim = c(0, 1), ylim = c(0.5, 4.5), yaxt = "n", ylab = "",
  xlab = "occurrence of mimicry, q", xaxs = "i", yaxs = "i")

abline(v = seq(0, 1, by = 0.1), col = gray(0.5, 0.4))

q_mu <- 0.5
q_theta <- 4
q_prior <- rbeta(10000, q_mu * q_theta, (1 - q_mu) * q_theta)

q_prior_mu <- mean(q_prior)
q_prior_lb <- HPDI(q_prior)[1]
q_prior_ub <- HPDI(q_prior)[2]

lines(c(q_prior_lb, q_prior_ub), c(1, 1))
points(q_prior_mu, 1, pch = 20)
text(q_prior_mu, 1, pos = 3, "prior")

q_single2 <- draws_of(post_single2$q)
q_single2_mu <- mean(q_single2)
q_single2_lb <- HPDI(q_single2)[1]
q_single2_ub <- HPDI(q_single2)[2]

lines(c(q_single2_lb, q_single2_ub), c(2, 2))
points(q_single2_mu, 2, pch = 20)
text(q_single2_mu, 2, pos = 3, "secondary sources only (18 mimic spp.)")

q_single1 <- draws_of(post_single1$q)
q_single1_mu <- mean(q_single1)
q_single1_lb <- HPDI(q_single1)[1]
q_single1_ub <- HPDI(q_single1)[2]

lines(c(q_single1_lb, q_single1_ub), c(3, 3))
points(q_single1_mu, 3, pch = 20)
text(q_single1_mu, 3, pos = 3, "primary sources only (22 mimic spp.)")

q_double <- draws_of(post_double$q)
q_double_mu <- mean(q_double)
q_double_lb <- HPDI(q_double)[1]
q_double_ub <- HPDI(q_double)[2]

lines(c(q_double_lb, q_double_ub), c(4, 4))
points(q_double_mu, 4, pch = 20)
text(q_double_mu, 4, pos = 3, "all sources (31 mimic spp.)")

dev.off()



spp$pr_mimic_mu_post <- pr_mimic_mu_post
spp$pr_mimic_se_post <- pr_mimic_se_post

spp$pr_mimic <- paste0(
  sprintf("%0.3f", spp$pr_mimic_mu_post), " (",
  sprintf("%0.3f", spp$pr_mimic_se_post), ")"
)

x <- spp[which(spp$English_name == "White-necked Raven"), c("English_name", "corvid_database_entries", "primary_count", "pr_mimic_mu_post", "pr_mimic")]

x <- spp[spp$pr_mimic_mu_post < 1, c("English_name", "corvid_database_entries", "primary_count", "pr_mimic_mu_post", "pr_mimic")]

out <- arrange(x, pr_mimic_mu_post)[1:10,]
out$pr_mimic_mu_post <- sprintf("%0.3f", out$pr_mimic_mu_post)
out <- select(out, -pr_mimic_mu_post)
write_pandoc_table(out, "figures/least_likely_ten.md", style = "rmarkdown")

out <- arrange(x, desc(pr_mimic_mu_post))[1:10,]
out$pr_mimic_mu_post <- sprintf("%0.3f", out$pr_mimic_mu_post)
out <- select(out, -pr_mimic_mu_post)
write_pandoc_table(out, "figures/most_likely_ten.md", style = "rmarkdown")


cat("use this to generate imputations\n")

spp$mimicry <- as.numeric(spp$mimicry_secondary | spp$mimicry_primary)

set.seed(1) # seed was `1` in submission version

n_imp <- 10

for (i in 1:n_imp) {
  varname <- paste0("mimicry_imputed", i)
  spp[[varname]] <- rbinom(nrow(spp), 1, spp$pr_mimic_mu_post)
}

write.csv(spp, "data/corvid_species_imputation.csv", row.names = FALSE)



cat("create composite evidence score\n")

A <- apply(pr_mimic_double, 2, mean)
A_lb <- apply(pr_mimic_double, 2, HPDI)[1,]
A_ub <- apply(pr_mimic_double, 2, HPDI)[2,]
q <- mean(draws_of(post_double$q))
a <- 1 - mean(draws_of(post_double$rho1)) 
b <- 1 - mean(draws_of(post_double$rho2)) 

tar <- which(A < 1)

log(a) / log(b) # each additional recording is like 0.12 academic papers on the subject!

log(b) / log(a) # 8 recordings per academic paper

# evidence score
evidence_score <- rep(NA, nrow(spp))
# all species expressed in units of k2, aka when k1 = 0
evidence_score[tar] <- log((A[tar] * (1 - q)) / ((1 - A[tar]) * q)) / log(b)

# note: the composite score only really makes sense for un-detected mimics, so the rest aren't plotted

# got it!
png("figures/composite_evidence_score.png", res = 300, units = "in", height = 5, width = 7)

plot(NULL, ylim = c(0, 1), xlim = c(0, 400), pch = 20, ylab = "probability of hidden mimicry", xlab = "composite evidence score (citations)", frame.plot = FALSE)
curve(pr_hidden_single_true(x, q, mean(draws_of(post_double$rho2))), add = TRUE, lty = 2)
for (i in 1:nrow(spp)) {
  if (A[i] < 1) {
    lines(c(evidence_score[i], evidence_score[i]), c(A_lb[i], A_ub[i]), col = gray(0.3, 0.5))
  }
}
points(evidence_score, A, pch = 20)

tar <- which(spp$English_name == "Rook")
text(evidence_score[tar] - 20, A[tar], label = bquote(italic(.(spp$scientific_name_short[tar]))), pos = 3)

tar <- which(spp$English_name == "House Crow")
text(evidence_score[tar] + 40, A[tar] - 0.03, label = bquote(italic(.(spp$scientific_name_short[tar]))), pos = 3)

tar <- which(spp$English_name == "Large-billed Crow")
text(evidence_score[tar] - 55, A[tar] + 0.04, label = bquote(italic(.(spp$scientific_name_short[tar]))), pos = 1)

tar <- which(spp$English_name == "California Scrub-Jay")
text(evidence_score[tar] + 30, A[tar], label = bquote(italic(.(spp$scientific_name_short[tar]))), pos = 3)

tar <- which(spp$English_name == "Transvolcanic Jay")
text(evidence_score[tar] + 40, A[tar], label = bquote(italic(.(spp$scientific_name_short[tar]))), pos = 3)

dev.off()



# got it!
png("figures/composite_evidence_score_log.png", res = 300, units = "in", height = 5, width = 7)

plot(NULL, ylim = c(0, 1), xlim = c(1, 18700), pch = 20, ylab = "probability of hidden mimicry", xlab = "composite evidence score (citations)", frame.plot = FALSE, log = "x")
curve(pr_hidden_single_true(x, q, mean(draws_of(post_double$rho2))), add = TRUE, lty = 2)
for (i in 1:nrow(spp)) {
  if (A[i] < 1) {
    lines(c(evidence_score[i], evidence_score[i]), c(A_lb[i], A_ub[i]), col = gray(0.3, 0.5))
  }
}
points(evidence_score, A, pch = 20)

tar <- which(spp$English_name == "Rook")
text(evidence_score[tar] - 20, A[tar], label = bquote(italic(.(spp$scientific_name_short[tar]))), pos = 3)

tar <- which(spp$English_name == "House Crow")
text(evidence_score[tar] + 40, A[tar] - 0.03, label = bquote(italic(.(spp$scientific_name_short[tar]))), pos = 3)

tar <- which(spp$English_name == "Large-billed Crow")
text(evidence_score[tar] - 55, A[tar] + 0.04, label = bquote(italic(.(spp$scientific_name_short[tar]))), pos = 1)

tar <- which(spp$English_name == "California Scrub-Jay")
text(evidence_score[tar] + 30, A[tar], label = bquote(italic(.(spp$scientific_name_short[tar]))), pos = 3)

tar <- which(spp$English_name == "Transvolcanic Jay")
text(evidence_score[tar] + 40, A[tar], label = bquote(italic(.(spp$scientific_name_short[tar]))), pos = 3)

dev.off()
