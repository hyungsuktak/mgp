############### Flexible covariance-based five-band DRW
############### Multivariate Matern specialization with nu_ij = 1/2
############### Simple valid construction of Apanasovich, Genton & Sun (2012)

options(digits = 11)
setwd("~/Desktop/")
dat <- read.table("1540", skip = 0, sep = "", header = FALSE)

bands <- c("u", "g", "r", "i", "z")
k <- length(bands)
nstart_default <- 5L
ncores_default <- max(1L, parallel::detectCores() - 1L)
seed_default <- rpois(1, 2026L)

band_data <- lapply(seq_len(k), function(j) {
  cols <- (3L * j - 2L):(3L * j)
  d <- dat[, cols, drop = FALSE]
  names(d) <- c("time", "value", "error")
  d <- d[is.finite(d$time) & is.finite(d$value) & is.finite(d$error), , drop = FALSE]
  d <- d[d$value != -99.990 & d$error > 0, , drop = FALSE]
  d$band <- j
  d
})

obs <- do.call(rbind, band_data)
time_origin <- min(obs$time)
obs$time_centered <- obs$time - time_origin
obs <- obs[order(obs$time_centered, obs$band), , drop = FALSE]
N <- nrow(obs)

############### Positive-definite auxiliary correlation V

z_pairs <- which(lower.tri(matrix(0, k, k)), arr.ind = TRUE)
n_corr <- nrow(z_pairs)

z_to_V <- function(z) {
  if (length(z) != n_corr) stop("Incorrect number of correlation parameters.")
  L <- diag(1, k)
  L[lower.tri(L)] <- z
  S <- tcrossprod(L)
  d <- sqrt(diag(S))
  V <- S / outer(d, d)
  diag(V) <- 1
  dimnames(V) <- list(bands, bands)
  V
}

############### Matern-DRW validity mapping
############### alpha_i = 1/tau_i
############### alpha_ij^2 = (alpha_i^2 + alpha_j^2)/2
############### rhoX_ij = v_ij * sqrt(2*tau_i*tau_j/(tau_i^2+tau_j^2))

matern_drw_mapping <- function(tau, V) {
  tau <- as.numeric(tau)
  tau_cross <- matrix(NA_real_, k, k, dimnames = list(bands, bands))
  rhoX <- matrix(NA_real_, k, k, dimnames = list(bands, bands))

  for (i in seq_len(k)) {
    for (j in seq_len(k)) {
      tau_cross[i, j] <- sqrt(2 / (tau[i]^(-2) + tau[j]^(-2)))
      atten <- sqrt(2 * tau[i] * tau[j] / (tau[i]^2 + tau[j]^2))
      rhoX[i, j] <- V[i, j] * atten
    }
  }
  diag(tau_cross) <- tau
  diag(rhoX) <- 1
  list(tau_cross = tau_cross, rhoX = rhoX)
}

############### Parameter unpacking
############### 5 means + 5 stationary SDs + 5 marginal taus + 10 V parameters = 25
############### covariance/dependence parameters alone = 20

n_internal <- 3L * k + n_corr
par_names <- c(
  paste0("mu_", bands),
  paste0("log_sd_", bands),
  paste0("log_tau_", bands),
  apply(z_pairs, 1, function(ind) paste0("z_", bands[ind[1]], bands[ind[2]]))
)

unpack_parameters <- function(par) {
  if (length(par) != n_internal) stop("Unexpected parameter-vector length.")
  mu <- par[seq_len(k)]
  sd_stationary <- exp(par[k + seq_len(k)])
  tau <- exp(par[2L * k + seq_len(k)])
  z <- par[(3L * k + 1L):length(par)]
  V <- z_to_V(z)
  map <- matern_drw_mapping(tau, V)

  list(
    mu = setNames(as.numeric(mu), bands),
    sd = setNames(as.numeric(sd_stationary), bands),
    tau = setNames(as.numeric(tau), bands),
    z = z,
    V = V,
    tau_cross = map$tau_cross,
    rhoX = map$rhoX
  )
}

############### Direct GP covariance and negative log-likelihood

build_covariance <- function(p, add_measurement_error = TRUE, jitter = 1e-10) {
  b <- obs$band
  tt <- obs$time_centered
  dt <- abs(outer(tt, tt, "-"))

  sd_outer <- outer(as.numeric(p$sd[b]), as.numeric(p$sd[b]))
  rho_obs <- p$rhoX[cbind(b[row(dt)], b[col(dt)])]
  tau_obs <- p$tau_cross[cbind(b[row(dt)], b[col(dt)])]

  K <- sd_outer * rho_obs * exp(-dt / tau_obs)
  K <- (K + t(K)) / 2

  if (add_measurement_error) diag(K) <- diag(K) + obs$error^2
  diag(K) <- diag(K) + jitter
  K
}

matern_drw_nll <- function(par, return_details = FALSE, penalty = 1e100, jitter = 1e-10) {
  p <- try(unpack_parameters(par), silent = TRUE)
  if (inherits(p, "try-error")) return(penalty)
  if (any(!is.finite(p$mu)) || any(!is.finite(p$sd)) || any(!is.finite(p$tau)) ||
      any(p$sd <= 0) || any(p$tau <= 0)) return(penalty)

  K <- try(build_covariance(p, TRUE, jitter), silent = TRUE)
  if (inherits(K, "try-error") || any(!is.finite(K))) return(penalty)

  U <- try(chol(K), silent = TRUE)
  if (inherits(U, "try-error")) return(penalty)

  residual <- obs$value - as.numeric(p$mu[obs$band])
  z <- forwardsolve(t(U), residual)
  nll <- 0.5 * (N * log(2 * pi) + 2 * sum(log(diag(U))) + sum(z^2))
  if (!is.finite(nll)) return(penalty)
  if (!return_details) return(nll)

  list(nll = nll, parameters = p, covariance = K, residual = residual)
}

############### Observation-adaptive univariate DRW warm starts
############### Marginal models are exact exponential (Matern 1/2) GPs.
############### Preliminary fits supply starts only; the final joint likelihood
############### is optimized with all marginal and dependence parameters free.

empirical_mean <- vapply(band_data, function(d) mean(d$value), numeric(1))
empirical_sd <- vapply(band_data, function(d) {
  v <- stats::var(d$value) - mean(d$error^2)
  sqrt(max(v, 1e-8))
}, numeric(1))

baseline <- diff(range(obs$time_centered))
if (!is.finite(baseline) || baseline <= 0) stop("Positive time baseline required.")

band_design <- lapply(band_data, function(d) {
  tt <- sort(unique(d$time))
  gaps <- diff(tt)
  gaps <- gaps[is.finite(gaps) & gaps > 0]
  T <- diff(range(tt))
  if (!is.finite(T) || T <= 0 || !length(gaps))
    stop("Each band needs at least two distinct observation times.")
  delta <- as.numeric(stats::quantile(gaps, 0.5, names = FALSE))
  c(cadence = delta, baseline = T)
})

# These are broad numerical bounds, not confidence intervals or prior bounds.
tau_lower <- min(vapply(band_design, function(x) x["cadence"], numeric(1))) / 10
tau_upper <- 10 * baseline

mu_lower <- vapply(band_data, function(d) min(d$value) - 5 * sd(d$value), numeric(1))
mu_upper <- vapply(band_data, function(d) max(d$value) + 5 * sd(d$value), numeric(1))
mu_lower[!is.finite(mu_lower)] <- empirical_mean[!is.finite(mu_lower)] - 10
mu_upper[!is.finite(mu_upper)] <- empirical_mean[!is.finite(mu_upper)] + 10

lower <- c(mu_lower, rep(log(1e-6), k), rep(log(tau_lower), k), rep(-8, n_corr))
upper <- c(mu_upper, rep(log(10), k), rep(log(tau_upper), k), rep(8, n_corr))
names(lower) <- names(upper) <- par_names

clip_start <- function(x) {
  x <- pmin(pmax(x, lower + 1e-8), upper - 1e-8)
  setNames(x, par_names)
}

# Profile the constant mean analytically at each univariate covariance setting.
# This is only a computational warm start; no marginal estimates are fixed.
fit_univariate_warm <- function(j, nstart = 5L, maxit = 500L) {
  d <- band_data[[j]]
  tt <- d$time - min(d$time)
  dt <- abs(outer(tt, tt, "-"))
  y <- d$value
  e2 <- d$error^2
  design <- band_design[[j]]
  lo <- max(tau_lower, design["cadence"])
  hi <- min(tau_upper, max(design["baseline"], lo * 1.01))
  if (hi <= lo) hi <- lo * 1.01
  nll1 <- function(v) {
    s <- exp(v[1]); tau <- exp(v[2])
    if (!is.finite(s) || !is.finite(tau)) return(1e50)
    K <- s^2 * exp(-dt / tau)
    diag(K) <- diag(K) + e2 + 1e-10
    U <- tryCatch(chol(K), error = function(e) NULL)
    if (is.null(U)) return(1e50)
    one <- rep(1, length(y))
    wy <- backsolve(U, forwardsolve(t(U), y))
    w1 <- backsolve(U, forwardsolve(t(U), one))
    mu <- sum(wy) / sum(w1)
    r <- y - mu
    z <- forwardsolve(t(U), r)
    0.5 * (length(y) * log(2*pi) + 2*sum(log(diag(U))) + sum(z^2))
  }
  tau0 <- c(sqrt(lo*hi), if (nstart > 1L) exp(runif(nstart-1L, log(lo), log(hi))))
  runs <- lapply(tau0, function(t0) {
    tryCatch(optim(c(log(empirical_sd[j]), log(t0)), nll1,
                   method = "L-BFGS-B",
                   lower = c(log(1e-6), log(tau_lower)),
                   upper = c(log(10), log(tau_upper)),
                   control = list(maxit = maxit)),
             error = function(e) list(value = Inf, convergence = 999L))
  })
  vals <- vapply(runs, function(x) x$value, numeric(1))
  if (!any(is.finite(vals)))
    return(c(mu = empirical_mean[j], sd = empirical_sd[j], tau = sqrt(lo*hi),
             nll = Inf, convergence = 999))
  best <- runs[[which.min(vals)]]
  s <- exp(best$par[1]); tau <- exp(best$par[2])
  K <- s^2 * exp(-dt / tau)
  diag(K) <- diag(K) + e2 + 1e-10
  U <- chol(K)
  solveK <- function(x) backsolve(U, forwardsolve(t(U), x))
  mu <- sum(solveK(y)) / sum(solveK(rep(1, length(y))))
  c(mu = mu, sd = s, tau = tau, nll = best$value,
    convergence = best$convergence)
}

# Innovation-informed correlations are not used here: asynchronous sampling
# requires a joint state-space calculation to estimate driving correlations.
# The deterministic independence/equicorrelation starts avoid this cost.

# Construct an equicorrelation matrix in the existing L L' parameterization.
V_to_z <- function(V) {
  U <- chol(V)
  L <- t(U)
  L <- sweep(L, 1L, diag(L), "/")
  as.numeric(L[lower.tri(L)])
}

make_joint_start <- function(warm, kind = c("independent", "equicor", "random", "adaptive")) {
  kind <- match.arg(kind)
  if (kind == "adaptive") {
    tau0 <- vapply(band_design, function(x) {
      lo <- max(tau_lower, x["cadence"])
      hi <- min(tau_upper, max(x["baseline"], lo*1.01))
      exp(runif(1L, log(lo), log(hi)))
    }, numeric(1))
    sd0 <- empirical_sd * exp(rnorm(k, 0, 0.25))
    mu0 <- empirical_mean
  } else {
    tau0 <- warm[, "tau"]
    sd0 <- warm[, "sd"]
    mu0 <- warm[, "mu"]
    if (kind == "random") {
      tau0 <- tau0 * exp(rnorm(k, 0, 0.5))
      sd0 <- sd0 * exp(rnorm(k, 0, 0.2))
    }
  }
  z0 <- switch(kind,
    independent = rep(0, n_corr),
    equicor = V_to_z(matrix(0.5, k, k) + diag(0.5, k)),
    random = rnorm(n_corr, 0, 0.4),
    adaptive = rnorm(n_corr, 0, 0.4))
  clip_start(c(mu0, log(sd0), log(tau0), z0))
}

############### Multi-start MLE

fit_covariance_matern_drw <- function(
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
  cat("Computing observation-adaptive univariate DRW warm starts...\n")
  warm <- t(vapply(seq_len(k), function(j) fit_univariate_warm(j), numeric(5)))
  rownames(warm) <- bands
  print(warm)
  kinds <- c("independent", "equicor", "random", "random", "adaptive")
  if (nstart > length(kinds))
    kinds <- c(kinds, rep(c("random", "adaptive"), length.out = nstart-length(kinds)))
  starts <- lapply(kinds[seq_len(nstart)], function(kind) make_joint_start(warm, kind))

  cat("============================================================\n")
  cat("Flexible covariance-based multivariate DRW MLE\n")
  cat("Multivariate Matern, nu_ij = 1/2\n")
  cat("Starts:", nstart, "\nCores :", ncores, "\n")
  cat("Observations:", N, "\n")
  cat("Total fitted parameters:", n_internal, "(including 5 means)\n")
  cat("============================================================\n")

  fit_one_start <- function(i) {
    run_start <- Sys.time()
    ans <- tryCatch(
      optim(
        par = starts[[i]], fn = matern_drw_nll,
        method = "L-BFGS-B", lower = lower, upper = upper,
        control = list(maxit = maxit, factr = 1e7, pgtol = 1e-7)
      ),
      error = function(e) list(par = starts[[i]], value = Inf,
                               convergence = 999L, message = conditionMessage(e))
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
    fits <- parallel::mclapply(ids, fit_one_start,
                               mc.cores = min(ncores, nstart), mc.preschedule = FALSE)
  }

  objective_values <- vapply(fits, function(x)
    if (is.null(x$value) || !is.finite(x$value)) Inf else x$value, numeric(1))
  if (all(!is.finite(objective_values))) stop("All optimization runs failed.")

  best_id <- which.min(objective_values)
  best <- fits[[best_id]]
  cat("\nBest preliminary start:", best_id, "\n")
  cat("Best preliminary nll  :", best$value, "\n")

  final <- optim(
    par = best$par, fn = matern_drw_nll,
    method = "L-BFGS-B", lower = lower, upper = upper,
    hessian = final_hessian,
    control = list(maxit = final_maxit, factr = 1e5, pgtol = 1e-8)
  )

  # Recompute the Hessian at the final parameter vector; negative eigenvalues
  # can indicate nonconvergence or numerical differentiation error.
  details <- matern_drw_nll(final$par, return_details = TRUE)
  p <- details$parameters
  logLik_value <- -final$value
  AIC_value <- -2 * logLik_value + 2 * length(final$par)

  result <- list(
    fit = final,
    means = p$mu,
    sd = p$sd,
    tau = p$tau,
    V = p$V,
    rhoX = p$rhoX,
    tau_cross = p$tau_cross,
    logLik = logLik_value,
    AIC = AIC_value,
    npar = length(final$par),
    best_start = best_id,
    univariate_warm = warm,
    start_types = kinds[seq_len(nstart)],
    all_fits = fits,
    objective_values = objective_values,
    time_origin = time_origin,
    model = "Flexible covariance-based Matern-DRW (nu=1/2)"
  )

  cat("Convergence:", final$convergence, "\n")
  cat("Message    :", final$message, "\n")
  cat("logLik     :", logLik_value, "\n")
  cat("AIC        :", AIC_value, "\n\n")
  result
}

############### Summaries

print_covariance_matern_summary <- function(fit) {
  tab <- data.frame(
    band = bands,
    mean = as.numeric(fit$means),
    stationary_SD = as.numeric(fit$sd),
    tau_days = as.numeric(fit$tau)
  )
  cat("\nBand-specific estimates\n-----------------------\n")
  print(tab, row.names = FALSE)
  cat("\nAuxiliary correlation matrix V\n------------------------------\n")
  print(round(fit$V, 4))
  cat("\nActual zero-lag process correlation R_X\n---------------------------------------\n")
  print(round(fit$rhoX, 4))
  cat("\nDerived cross-timescale matrix tau_ij (days)\n--------------------------------------------\n")
  print(round(fit$tau_cross, 2))
  cat("\nlogLik     =", fit$logLik, "\n")
  cat("AIC        =", fit$AIC, "\n")
  cat("Parameters =", fit$npar, "(including", k, "means)\n")
  invisible(tab)
}

compare_best_starts <- function(fit, n = 10L) {
  values <- sort(fit$objective_values[is.finite(fit$objective_values)])
  values <- head(values, n)
  data.frame(rank = seq_along(values), negative_log_likelihood = values,
             difference_from_best = values - min(values))
}

hessian_diagnostics <- function(fit) {
  H <- fit$fit$hessian
  if (is.null(H)) return(list(valid = FALSE, message = "No Hessian was stored."))
  H <- (H + t(H)) / 2
  eig <- try(eigen(H, symmetric = TRUE, only.values = TRUE)$values, silent = TRUE)
  if (inherits(eig, "try-error")) return(list(valid = FALSE, message = "Hessian eigendecomposition failed."))
  valid <- all(is.finite(eig)) && min(eig) > 0
  list(valid = valid,
       message = if (valid) "The Hessian is positive definite." else
         "The Hessian is not positive definite; Hessian-based SEs are unreliable.",
       minimum_eigenvalue = min(eig), maximum_eigenvalue = max(eig),
       nonpositive = sum(eig <= 0), eigenvalues = eig)
}

############### Covariance / correlation / spectrum

cross_correlation_ij <- function(i, j, lag, fit) {
  fit$rhoX[i, j] * exp(-abs(lag) / fit$tau_cross[i, j])
}

# Fourier convention: S(omega) = integral C(h) exp(-i omega h) dh.
# A common 1/(2*pi) multiplier can be added if desired.
psd_ij <- function(i, j, omega, fit) {
  a <- fit$sd[i] * fit$sd[j] * fit$rhoX[i, j]
  tij <- fit$tau_cross[i, j]
  2 * a * tij / (1 + omega^2 * tij^2)
}

coherence_ij <- function(i, j, omega, fit) {
  if (i == j) return(rep(1, length(omega)))
  Sij <- psd_ij(i, j, omega, fit)
  Sii <- psd_ij(i, i, omega, fit)
  Sjj <- psd_ij(j, j, omega, fit)
  Sij^2 / (Sii * Sjj)
}

phase_ij <- function(i, j, omega, fit) {
  # Cross-spectrum is real under this symmetric covariance construction.
  # Positive S_ij -> phase 0; negative S_ij -> phase pi.
  s <- psd_ij(i, j, omega, fit)
  ifelse(s >= 0, 0, pi)
}

spectral_matrix <- function(omega, fit) {
  S <- matrix(NA_real_, k, k, dimnames = list(bands, bands))
  for (i in seq_len(k)) for (j in seq_len(k)) S[i, j] <- psd_ij(i, j, omega, fit)
  (S + t(S)) / 2
}

check_spectral_validity <- function(fit, omega_grid = 10^seq(-6, 3, length.out = 1000)) {
  min_eigs <- vapply(omega_grid, function(w)
    min(eigen(spectral_matrix(w, fit), symmetric = TRUE, only.values = TRUE)$values), numeric(1))
  data.frame(minimum_eigenvalue = min(min_eigs),
             omega_at_minimum = omega_grid[which.min(min_eigs)],
             all_nonnegative = all(min_eigs >= -1e-10))
}

############### Simple plots for comparison with mDRW

plot_pair_dependence <- function(fit, band1 = "u", band2 = "z",
                                 lag = seq(-2000, 2000, length.out = 1000),
                                 omega = 10^seq(-5, 0, length.out = 500)) {
  i <- match(band1, bands); j <- match(band2, bands)
  if (is.na(i) || is.na(j)) stop("Unknown band name.")
  oldpar <- par(no.readonly = TRUE); on.exit(par(oldpar))
  par(mfrow = c(1, 3))

  cc <- cross_correlation_ij(i, j, lag, fit)
  plot(lag, cc, type = "l", lwd = 2, xlab = "Lag (days)", ylab = "Cross-correlation",
       main = paste0(band1, "-", band2))
  abline(v = 0, lty = 3)

  coh <- coherence_ij(i, j, omega, fit)
  plot(omega, coh, type = "l", log = "x", lwd = 2,
       xlab = expression(omega), ylab = "Squared coherence",
       ylim = c(0, min(1, max(coh, na.rm = TRUE) * 1.05)))

  ph <- phase_ij(i, j, omega, fit)
  plot(omega, ph, type = "l", log = "x", lwd = 2,
       xlab = expression(omega), ylab = "Phase (radians)")
  abline(h = 0, lty = 3)
  invisible(list(lag = lag, cross_correlation = cc, omega = omega,
                 coherence = coh, phase = ph))
}

############### Full fit

fit_full_cov <- fit_covariance_matern_drw(
  nstart = 5,
  ncores = ncores_default,
  maxit = 1500,
  final_maxit = 5000,
  final_hessian = TRUE,
  seed = rpois(1, 2026)
)

print_covariance_matern_summary(fit_full_cov)
compare_best_starts(fit_full_cov)
hessian_diagnostics(fit_full_cov)
check_spectral_validity(fit_full_cov)

# Example dependence figure; change bands as desired.
plot_pair_dependence(fit_full_cov, "u", "z")

