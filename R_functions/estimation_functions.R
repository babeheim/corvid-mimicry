
# a draws_of object
summarize_draws_rvars <- function(samples, format = TRUE) {
  out <- list()
  for (i in seq_along(samples)) {
    if (length(grep("_vll", names(samples)[i])) == 0) {
      draws <- draws_of(samples[[i]])
      if (length(dim(samples[[i]])) == 1) {
        if (dim(samples[[i]]) == 1) {
          add <- summarize_samples(draws, format = format)
          add$label <- names(samples)[i]
          out <- c(out, list(add))
        }
        if (1 < dim(samples[[i]]) & dim(samples[[i]]) < 7) {
          add <- vector("list", length(samples[[i]]))
          for (j in seq_len(length(samples[[i]]))) {
            add[[j]] <- summarize_samples(draws[, j], format = format)
            add[[j]]$label <- paste(names(samples)[i], j, sep = "_")
          }
          out <- c(out, add)
        }
      }
    }
  }
  return(out)
}

update_pr_right <- function(n_right, n_left, time, n = 10000, a_prior = 2, b_prior = 4) {
  lambda_right <- rgamma(n, shape = a_prior + n_right, rate = b_prior + time)
  lambda_left <- rgamma(n, shape = a_prior + n_left, rate = b_prior + time)
  lambda_right / (lambda_right + lambda_left)
}

format_estimates <- function(est) {
  calcs <- list()
  for (i in seq_along(est)) {
    if (hasName(est[[i]], "label")) {
      calc_name <- est[[i]]$label
      calcs[[paste0(calc_name, "Mean")]] <- sprintf("%1.2f", est[[i]]$mean)
      calcs[[paste0(calc_name, "SD")]] <- sprintf("%1.2f", est[[i]]$sd)
      calcs[[paste0(calc_name, "LB")]] <- sprintf("%1.2f", est[[i]]$lb)
      calcs[[paste0(calc_name, "UB")]] <- sprintf("%1.2f", est[[i]]$ub)
      calcs[[paste0(calc_name, "PSign")]] <- est[[i]]$psign
    }
  }
  return(calcs)
}


chain_plotter <- function(par_name, s) {

  h <- (0:(length(s) - 1) %/% n_iter) + 1

  plot(NULL, xlim = c(1, n_iter), ylim = range(s), ylab = "", main = par_name)

  n_chains <- max(h)
  chain_cols <- plasma(n_chains)

  for (i in 1:n_chains) {
    points(s[h == i], type = "l", col = col_alpha(chain_cols[i], 0.8))
  }

}


extract_diagnostics <- function(cmdstan_fit) {

  out <- list()

  diag <- cmdstan_fit$cmdstan_diagnose()$stdout
  out$rhat_good <- grepl("Split R-hat values satisfactory all parameters.", diag)
  out$ess_good <- grepl("Effective sample size satisfactory.", diag)
  out$energy_good <- grepl("E-BFMI satisfactory for all transitions.", diag)
  out$divergence_good <- grepl("No divergent transitions found.", diag)

  x <- cmdstan_fit$cmdstan_summary()$stdout
  pattern <- "Sampling took.*\\n"
  m <- gregexpr(pattern, x, perl = TRUE)
  sampling_time <- regmatches(x, m)[[1]]
  sampling_time <- gsub("\\n$", "", sampling_time)

  s <- sampling_time
  s <- gsub("Sampling .*, ", "", s)
  s <- gsub(" total", "", s)
  units <- gsub("\\d.*\\s", "", s)
  s <- as.numeric(gsub("\\s.*$", "", s))
  tar <- which(units  ==  "seconds")
  if (length(tar) > 0) s[tar] <- s[tar] / 60
  tar <- which(units  ==  "hours")
  if (length(tar) > 0) s[tar] <- s[tar] * 60
  out$sampling_time_min <- round(s, 1)

  return(out)
}


summarize_samples_short <- function(samples, digits = 2, label = NA) {
  if (mean(is.finite(samples)) < 0.9) stop("too many Inf")
  samples <- samples[is.finite(samples)]
  fmt <- paste0("%0.", digits, "f")
  out <- list(
    mean = sprintf(fmt, mean(samples)),
    lb = sprintf(fmt, HPDI(samples)[1]),
    ub = sprintf(fmt, HPDI(samples)[2])
  )
  if (!is.na(label)) {
    out <- c(list(label = label), out)
  }
  out <- paste0(out$mean, " (", out$lb, ", ", out$ub, ")")
  return(out)
}

summarize_samples <- function(samples, format = FALSE, label = NA, p = FALSE) {
  if (mean(is.finite(samples)) < 0.9) stop("too many Inf")
  samples <- samples[is.finite(samples)]
  out <- list(
    mean = mean(samples),
    sd = sd(samples),
    lb = HPDI(samples)[1],
    ub = HPDI(samples)[2]
  )
  if (p == TRUE) out <- c(out, psign = psign(samples))
  if (!is.na(label)) {
    out <- c(list(label = label), out)
  }
  if (format == TRUE) {
    out$mean <- sprintf("%1.2f", out$mean)
    out$sd <- sprintf("%1.2f", out$sd)
    out$lb <- sprintf("%1.2f", out$lb)
    out$ub <- sprintf("%1.2f", out$ub)
  }
  return(out)
}

psign <- function(samples, n_digits = 3) {
  if (all(samples == 0)) stop("all samples are zero!")
  if (mean(samples) > 0) output <- round(mean(samples < 0), n_digits)
  if (mean(samples) < 0) output <- round(mean(samples > 0), n_digits)
  float_digits <- paste0("%.", n_digits, "f")
  output <- sprintf(float_digits, output)
  null_entry <- paste0("0.", paste0(rep("0", n_digits), collapse = ""))
  output[output  ==  null_entry] <- paste0("$<0.",
    paste0(rep("0", n_digits - 1), collapse = ""), "1$")
  return(output)
}

