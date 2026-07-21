
############### Section 4.1: A separable five-band DRW model

options(digits = 11)
setwd("~/Desktop/")
dat <- read.table("1540", skip = 0, sep = "", header = FALSE)

bands <- c("u", "g", "r", "i", "z")
k <- length(bands)
nstart_default <- 30L
ncores_default <- max(1L, parallel::detectCores() - 1L)
seed_default <- 2026L

band_data <- lapply(seq_len(k), function(j) {
  cols <- (3L * j - 2L):(3L * j)
  d <- dat[, cols, drop = FALSE]
  names(d) <- c("time", "value", "error")
  d <- d[is.finite(d$time) & is.finite(d$value) & is.finite(d$error), , drop = FALSE]
  d <- d[d$value != -99.990 & d$error > 0, , drop = FALSE]
  d$band <- j
  d
})

median(diff(band_data[[1]][, 1])); median(band_data[[1]][, 2]); median(band_data[[1]][, 3])
median(diff(band_data[[2]][, 1])); median(band_data[[2]][, 2]); median(band_data[[2]][, 3])
median(diff(band_data[[3]][, 1])); median(band_data[[3]][, 2]); median(band_data[[3]][, 3])
median(diff(band_data[[4]][, 1])); median(band_data[[4]][, 2]); median(band_data[[4]][, 3])
median(diff(band_data[[5]][, 1])); median(band_data[[5]][, 2]); median(band_data[[5]][, 3])

obs <- do.call(rbind, band_data)
time_origin <- min(obs$time)
obs$time_centered <- obs$time - time_origin
obs <- obs[order(obs$time_centered, obs$band), , drop = FALSE]
N <- nrow(obs)
unique_times <- sort(unique(obs$time_centered))
n_time <- length(unique_times)
time_groups <- split(seq_len(N), match(obs$time_centered, unique_times))

############### Figure 2

par(mfcol = c(1, 1), font = 1, font.lab = 2, font.axis = 1, cex = 1.1,
    mai = c(0.7, 0.9, 0.5, 0.1), mgp = c(2.5, 0.5, 0), las = 0)
cex.size <- 0.5
plot(band_data[[1]][, 1], band_data[[1]][, 2], col = 1, ylim = c(21.8, 18.7), xlim = c(51000, 54500),
     xlab = "", ylab = "", pch = 1, cex = cex.size)
mtext(side = 1, text = "Observation time (MJD)", line = 2, cex = 1.3)
mtext(side = 2, text = "Magnitude + constant", line = 1.9, cex = 1.5)
title("")
legend("top", rev(c("z-band", "i-band", "r-band", "g-band", "u-band")), col = rev(c(6, 4, 3, 2, 1)),
       pch = c(1 : 5), ncol = 5, bty = "n")

points(band_data[[2]][, 1], band_data[[2]][, 2] + 1.5, col = 2, pch = 2, cex = cex.size)
points(band_data[[3]][, 1], band_data[[3]][, 2] + 1.3, col = 3, pch = 3, cex = cex.size)
points(band_data[[4]][, 1], band_data[[4]][, 2] + 1.1, col = 4, pch = 4, cex = cex.size)
points(band_data[[5]][, 1], band_data[[5]][, 2] + 0.7, col = 6, pch = 5, cex = cex.size)

############### Positive-definite correlation parameterization

z_pairs <- which(lower.tri(matrix(0, k, k)), arr.ind = TRUE)
n_corr <- nrow(z_pairs)

z_to_R <- function(z) {
  if (length(z) != n_corr) stop("Incorrect number of correlation parameters.")
  L <- diag(1, k)
  L[lower.tri(L)] <- z
  C <- tcrossprod(L)
  s <- sqrt(diag(C))
  R <- C / outer(s, s)
  diag(R) <- 1
  dimnames(R) <- list(bands, bands)
  R
}

############### Parameter unpacking

n_internal <- 2L * k + 1L + n_corr
par_names <- c(
  paste0("mu_", bands),
  paste0("log_sigma_", bands),
  "log_tau",
  apply(z_pairs, 1, function(ind) paste0("z_", bands[ind[1]], bands[ind[2]]))
)

unpack_parameters <- function(par) {
  if (length(par) != n_internal) stop("Unexpected parameter-vector length.")
  mu <- par[seq_len(k)]
  sigma <- exp(par[k + seq_len(k)])
  tau <- exp(par[2L * k + 1L])
  z <- par[(2L * k + 2L):length(par)]
  list(
    mu = setNames(as.numeric(mu), bands),
    sigma = setNames(as.numeric(sigma), bands),
    tau = unname(tau),
    z = z,
    R = z_to_R(z)
  )
}

############### Kalman negative log-likelihood

kalman_nll <- function(par, return_details = FALSE, penalty = 1e100, jitter = 1e-10) {
  p <- try(unpack_parameters(par), silent = TRUE)
  if (inherits(p, "try-error")) return(penalty)
  if (any(!is.finite(p$mu)) || any(!is.finite(p$sigma)) || !is.finite(p$tau) || p$tau <= 0) return(penalty)

  D_sigma <- diag(as.numeric(p$sigma), k, k)
  P_inf <- (p$tau / 2) * D_sigma %*% p$R %*% D_sigma
  P_inf <- (P_inf + t(P_inf)) / 2

  m <- rep(0, k)
  P <- P_inf
  previous_time <- unique_times[1]
  nll <- 0

  filtered_means <- if (return_details) matrix(NA_real_, n_time, k, dimnames = list(NULL, bands)) else NULL
  filtered_covariances <- if (return_details) vector("list", n_time) else NULL
  innovations <- if (return_details) vector("list", n_time) else NULL

  for (g in seq_along(unique_times)) {
    current_time <- unique_times[g]

    if (g > 1L) {
      delta <- current_time - previous_time
      if (!is.finite(delta) || delta < 0) return(penalty)
      a <- exp(-delta / p$tau)
      m <- a * m
      Q <- (1 - a^2) * P_inf
      P <- a^2 * P + Q
      P <- (P + t(P)) / 2
    }

    idx <- time_groups[[g]]
    y_g <- obs$value[idx]
    e_g <- obs$error[idx]
    b_g <- obs$band[idx]
    m_obs <- length(idx)

    H_g <- matrix(0, m_obs, k)
    H_g[cbind(seq_len(m_obs), b_g)] <- 1

    residual <- y_g - as.numeric(p$mu[b_g]) - as.numeric(H_g %*% m)
    V_g <- diag(e_g^2, m_obs)
    S_g <- H_g %*% P %*% t(H_g) + V_g
    S_g <- (S_g + t(S_g)) / 2
    diag(S_g) <- diag(S_g) + jitter

    U <- try(chol(S_g), silent = TRUE)
    if (inherits(U, "try-error")) return(penalty)

    z <- forwardsolve(t(U), residual)
    nll <- nll + 0.5 * (m_obs * log(2 * pi) + 2 * sum(log(diag(U))) + sum(z^2))

    PHt <- P %*% t(H_g)
    S_inv_t_PHt <- backsolve(U, forwardsolve(t(U), t(PHt)))
    K_gain <- t(S_inv_t_PHt)
    m <- m + as.numeric(K_gain %*% residual)

    I_k <- diag(k)
    KH <- K_gain %*% H_g
    P <- (I_k - KH) %*% P %*% t(I_k - KH) + K_gain %*% V_g %*% t(K_gain)
    P <- (P + t(P)) / 2

    if (any(!is.finite(m)) || any(!is.finite(P)) || !is.finite(nll)) return(penalty)

    if (return_details) {
      filtered_means[g, ] <- m
      filtered_covariances[[g]] <- P
      innovations[[g]] <- residual
    }

    previous_time <- current_time
  }

  if (!return_details) return(nll)

  list(
    nll = nll,
    parameters = p,
    stationary_covariance = P_inf,
    unique_times = unique_times,
    filtered_state_means = filtered_means,
    filtered_state_covariances = filtered_covariances,
    innovations = innovations
  )
}

############### Starting values and bounds

empirical_mean <- vapply(band_data, function(d) mean(d$value), numeric(1))
empirical_sd <- vapply(band_data, function(d) {
  ans <- sd(d$value)
  if (!is.finite(ans) || ans <= 0) 1 else ans
}, numeric(1))

baseline <- max(obs$time_centered) - min(obs$time_centered)
if (!is.finite(baseline) || baseline <= 0) baseline <- 1000

tau_lower <- 1
tau_upper <- max(1e4, 2 * baseline)

mu_lower <- vapply(band_data, function(d) min(d$value) - 5 * sd(d$value), numeric(1))
mu_upper <- vapply(band_data, function(d) max(d$value) + 5 * sd(d$value), numeric(1))
mu_lower[!is.finite(mu_lower)] <- empirical_mean[!is.finite(mu_lower)] - 10
mu_upper[!is.finite(mu_upper)] <- empirical_mean[!is.finite(mu_upper)] + 10

lower <- c(mu_lower, rep(log(1e-8), k), log(tau_lower), rep(-8, n_corr))
upper <- c(mu_upper, rep(log(10), k), log(tau_upper), rep(8, n_corr))
names(lower) <- names(upper) <- par_names

make_start <- function() {
  tau0 <- exp(runif(1, log(30), log(min(1000, tau_upper))))
  sigma0 <- empirical_sd * sqrt(2 / tau0) * exp(rnorm(k, 0, 0.35))
  sigma0 <- pmax(sigma0, 1e-7)
  z0 <- rnorm(n_corr, 0, 0.35)
  start <- c(
    empirical_mean + rnorm(k, 0, 0.1 * empirical_sd),
    log(sigma0),
    log(tau0),
    z0
  )
  start <- pmin(pmax(start, lower + 1e-8), upper - 1e-8)
  setNames(start, par_names)
}

############### Multi-start MLE

fit_kalman_separable_mdrw <- function(
    nstart = nstart_default,
    ncores = ncores_default,
    maxit = 1500L,
    final_maxit = 5000L,
    final_hessian = TRUE,
    seed = seed_default) {

  nstart <- as.integer(nstart)
  ncores <- as.integer(ncores)
  if (nstart < 1L) stop("nstart must be at least 1.")

  set.seed(seed)
  starts <- replicate(nstart, make_start(), simplify = FALSE)

  cat("============================================================\n")
  cat("Kalman-filter separable multivariate DRW MLE\n")
  cat("Starts:", nstart, "\nCores :", ncores, "\n")
  cat("Observations:", N, "\nUnique epochs:", n_time, "\n")
  cat("============================================================\n")

  fit_one_start <- function(i) {
    run_start <- Sys.time()
    ans <- tryCatch(
      optim(
        par = starts[[i]],
        fn = kalman_nll,
        method = "L-BFGS-B",
        lower = lower,
        upper = upper,
        control = list(maxit = maxit, factr = 1e7, pgtol = 1e-7)
      ),
      error = function(e) list(
        par = starts[[i]], value = Inf, convergence = 999L,
        message = conditionMessage(e)
      )
    )
    elapsed <- as.numeric(difftime(Sys.time(), run_start, units = "mins"))
    cat(sprintf("[%3d/%3d] %.2f min | nll=%s | conv=%d\n",
                i, nstart, elapsed,
                if (is.finite(ans$value)) sprintf("%.6f", ans$value) else "Inf",
                ans$convergence))
    ans$start_id <- i
    ans$elapsed_minutes <- elapsed
    ans
  }

  ids <- seq_len(nstart)
  if (.Platform$OS.type == "windows" || ncores <= 1L) {
    fits <- lapply(ids, fit_one_start)
  } else {
    fits <- parallel::mclapply(
      ids, fit_one_start,
      mc.cores = min(ncores, nstart),
      mc.preschedule = FALSE
    )
  }

  objective_values <- vapply(fits, function(x) {
    if (is.null(x$value) || !is.finite(x$value)) Inf else x$value
  }, numeric(1))

  if (all(!is.finite(objective_values))) stop("All optimization runs failed.")

  best_id <- which.min(objective_values)
  best <- fits[[best_id]]
  cat("\nBest preliminary start:", best_id, "\n")
  cat("Best preliminary nll  :", best$value, "\n")

  final <- optim(
    par = best$par,
    fn = kalman_nll,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    hessian = final_hessian,
    control = list(maxit = final_maxit, factr = 1e5, pgtol = 1e-8)
  )

  details <- kalman_nll(final$par, return_details = TRUE)
  logLik_value <- -final$value
  npar_total <- length(final$par)
  AIC_value <- -2 * logLik_value + 2 * npar_total

  result <- list(
    fit = final,
    means = details$parameters$mu,
    sigma = details$parameters$sigma,
    tau = details$parameters$tau,
    R = details$parameters$R,
    stationary_covariance = details$stationary_covariance,
    logLik = logLik_value,
    AIC = AIC_value,
    npar = npar_total,
    best_start = best_id,
    all_fits = fits,
    objective_values = objective_values,
    time_origin = time_origin,
    filtered_state_means = details$filtered_state_means,
    filtered_state_covariances = details$filtered_state_covariances,
    unique_times = unique_times,
    model = "Kalman separable covariance-based multivariate DRW"
  )

  cat("Convergence:", final$convergence, "\n")
  cat("Message    :", final$message, "\n")
  cat("logLik     :", logLik_value, "\n")
  cat("AIC        :", AIC_value, "\n\n")
  result
}

############### Summaries and diagnostics

print_mle_summary <- function(fit) {
  tab <- data.frame(
    band = bands,
    mean = as.numeric(fit$means),
    sigma_diffusion = as.numeric(fit$sigma),
    common_tau_days = rep(fit$tau, k),
    stationary_SD = as.numeric(fit$sigma * sqrt(fit$tau / 2))
  )
  cat("\nBand-specific estimates\n-----------------------\n")
  print(tab, row.names = FALSE)
  cat("\nEstimated correlation matrix R\n--------------------------------\n")
  print(round(fit$R, 4))
  cat("\nStationary covariance matrix\n----------------------------\n")
  print(round(fit$stationary_covariance, 6))
  cat("\nCommon tau =", fit$tau, "days\n")
  cat("logLik     =", fit$logLik, "\n")
  cat("AIC        =", fit$AIC, "\n")
  cat("Parameters =", fit$npar, "\n")
  invisible(tab)
}

compare_best_starts <- function(fit, n = 10L) {
  values <- sort(fit$objective_values[is.finite(fit$objective_values)])
  values <- head(values, n)
  data.frame(
    rank = seq_along(values),
    negative_log_likelihood = values,
    difference_from_best = values - min(values)
  )
}

hessian_diagnostics <- function(fit) {
  H <- fit$fit$hessian
  if (is.null(H)) return(list(valid = FALSE, message = "No Hessian was stored."))
  H <- (H + t(H)) / 2
  eig <- try(eigen(H, symmetric = TRUE, only.values = TRUE)$values, silent = TRUE)
  if (inherits(eig, "try-error")) return(list(valid = FALSE, message = "Hessian eigendecomposition failed."))
  valid <- all(is.finite(eig)) && min(eig) > 0
  list(
    valid = valid,
    message = if (valid) "The Hessian is positive definite." else
      "The Hessian is not positive definite; Hessian-based SEs are unreliable.",
    minimum_eigenvalue = min(eig),
    maximum_eigenvalue = max(eig),
    nonpositive = sum(eig <= 0),
    eigenvalues = eig
  )
}

extract_parameter_se <- function(fit) {
  info <- hessian_diagnostics(fit)
  if (!info$valid) {
    warning(info$message)
    return(NULL)
  }
  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    stop("Install the numDeriv package to compute delta-method SEs.")
  }

  H <- (fit$fit$hessian + t(fit$fit$hessian)) / 2
  V_internal <- solve(H)

  natural_parameters <- function(par) {
    p <- unpack_parameters(par)
    rho_values <- p$R[lower.tri(p$R)]
    names(rho_values) <- apply(z_pairs, 1, function(ind) {
      paste0("rho_", bands[ind[1]], bands[ind[2]])
    })
    c(
      setNames(p$mu, paste0("mu_", bands)),
      setNames(p$sigma, paste0("sigma_", bands)),
      tau = p$tau,
      rho_values
    )
  }

  estimate <- natural_parameters(fit$fit$par)
  J <- numDeriv::jacobian(natural_parameters, fit$fit$par, method = "Richardson")
  V_natural <- J %*% V_internal %*% t(J)
  V_natural <- (V_natural + t(V_natural)) / 2

  if (any(!is.finite(diag(V_natural))) || any(diag(V_natural) < 0)) {
    warning("The transformed covariance matrix has invalid diagonal entries.")
    return(NULL)
  }

  data.frame(
    parameter = names(estimate),
    estimate = as.numeric(estimate),
    SE = sqrt(diag(V_natural)),
    row.names = NULL
  )
}

############### PSD and coherence

psd_ij <- function(i, j, omega, fit) {
  fit$R[i, j] * fit$sigma[i] * fit$sigma[j] * fit$tau^2 /
    (1 + omega^2 * fit$tau^2)
}

coherence_ij <- function(i, j, omega, fit) {
  if (i == j) return(rep(1, length(omega)))
  rep(fit$R[i, j]^2, length(omega))
}

spectral_matrix <- function(omega, fit) {
  scale <- fit$tau^2 / (1 + omega^2 * fit$tau^2)
  D_sigma <- diag(as.numeric(fit$sigma), k, k)
  S <- scale * D_sigma %*% fit$R %*% D_sigma
  dimnames(S) <- list(bands, bands)
  S
}

check_spectral_validity <- function(
    fit,
    omega_grid = 10^seq(-6, 3, length.out = 1000)) {
  min_eigs <- vapply(omega_grid, function(w) {
    min(eigen(spectral_matrix(w, fit), symmetric = TRUE, only.values = TRUE)$values)
  }, numeric(1))
  data.frame(
    minimum_eigenvalue = min(min_eigs),
    omega_at_minimum = omega_grid[which.min(min_eigs)],
    all_nonnegative = all(min_eigs >= -1e-10)
  )
}

############### Figure 3

plot_marginal_psds <- function(
    fit,
    omega = 10^seq(-5, 1, length.out = 500),
    outfile = NULL) {
  if (!is.null(outfile)) {
    pdf(outfile, width = 7, height = 5)
    on.exit(dev.off(), add = TRUE)
  }
  P <- sapply(seq_len(k), function(j) psd_ij(j, j, omega, fit))
  matplot(
    omega, P, type = "l", log = "xy",
    lty = seq_len(k), lwd = 2,
    xlab = expression(omega), ylab = "Marginal PSD"
  )
  legend("bottomleft", legend = bands, lty = seq_len(k), lwd = 2, bty = "n")
  invisible(data.frame(omega = omega, P))
}

plot_cross_psds <- function(
    fit,
    omega = 10^seq(-5, 1, length.out = 500),
    absolute = TRUE,
    outfile = NULL,
    legend_position = "bottomleft") {

  pairs <- t(combn(seq_along(bands), 2))
  n_pairs <- nrow(pairs)

  cross_psd <- sapply(
    seq_len(n_pairs),
    function(m) {
      i <- pairs[m, 1]
      j <- pairs[m, 2]

      value <- psd_ij(
        i = i,
        j = j,
        omega = omega,
        fit = fit
      )

      if (absolute) abs(value) else value
    }
  )

  pair_names <- apply(
    pairs,
    1,
    function(ind) {
      paste0(bands[ind[1]], "-", bands[ind[2]])
    }
  )

  colnames(cross_psd) <- pair_names

  if (!is.null(outfile)) {
    pdf(outfile, width = 8, height = 6)
    on.exit(dev.off(), add = TRUE)
  }

  if (absolute) {

    matplot(
      omega,
      cross_psd,
      type = "l",
      log = "xy",
      lty = seq_len(n_pairs),
      lwd = 2,
      xlab = expression(omega),
      ylab = expression(abs(S[ij](omega)))
    )

  } else {

    matplot(
      omega,
      cross_psd,
      type = "l",
      log = "x",
      lty = seq_len(n_pairs),
      lwd = 2,
      xlab = expression(omega),
      ylab = expression(S[ij](omega))
    )

    abline(h = 0, lty = 3)
  }

  legend(
    legend_position,
    legend = pair_names,
    lty = seq_len(n_pairs),
    lwd = 2,
    bty = "n",
    ncol = 2
  )

  invisible(
    list(
      omega = omega,
      cross_psd = cross_psd,
      pairs = pairs
    )
  )
}

############### Test fit

fit_test <- fit_kalman_separable_mdrw(
  nstart = 1,
  ncores = min(3L, ncores_default),
  maxit = 400,
  final_maxit = 1000,
  final_hessian = TRUE,
  seed = rpois(1, 1000)
)

print_mle_summary(fit_test)
compare_best_starts(fit_test)
hessian_diagnostics(fit_test)
check_spectral_validity(fit_test)

############### Fit

fit_full <- fit_kalman_separable_mdrw(
  nstart = 30,
  ncores = ncores_default,
  maxit = 1500,
  final_maxit = 5000,
  final_hessian = TRUE,
  seed = rpois(1, 1000)
)

print_mle_summary(fit_full)
compare_best_starts(fit_full)
hessian_diagnostics(fit_full)
check_spectral_validity(fit_full)
extract_covariance_parameter_se(fit_full)

plot_marginal_psds(fit_full)
plot_pairwise_coherence(fit_full)


