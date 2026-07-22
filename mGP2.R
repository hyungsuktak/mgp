
############### Section 4.2: Fixed-width latent-process DRW models

options(digits = 11)
setwd("~/Desktop/")
dat.h <- read.table("ha.txt", skip = 1, sep = "", header = FALSE)
dat.c <- read.table("cont.txt", skip = 1, sep = "", header = FALSE)

cont <- dat.c[dat.c[, 4] == "sdss_spec_r", ]
ha   <- dat.h

median(diff(cont[, 1])); median(cont[, 2]); median(cont[, 3])
median(diff(ha[, 1])); median(ha[, 2]); median(ha[, 3])

############### Figure 1

par(mfrow = c(2, 1),
    mar = c(2.1, 2.1, 1.2, 0.4),
    oma = c(3, 3, 1, 1), font = 2, font.axis = 2, cex = 1)
cex.size <- 1

plot(cont[, 1], cont[, 2], 
     col = 1, xlim = c(min(ha[, 1], cont[, 1]),  
                       max(ha[, 1], cont[, 1])),
     ylim = c(3, 16),
     xlab = "", ylab = , pch = 1, cex = cex.size, main = "")
legend("topleft", "Continuum light curve", bty = "n")
mtext("Continuum flux", side = 2, line = 2)
arrows(cont[, 1], cont[, 2] - cont[, 3], cont[, 1], cont[, 2] + cont[, 3], 
       angle = 90, code = 3, length = 0.02, col = "gray70")
title("SDSS-RM840", line = 0.5)

plot(dat.h[, 1], dat.h[, 2], col = 1, xlim = c(min(dat.h[, 1], dat.c[, 1]), max(dat.h[, 1], dat.c[, 1])),
     ylim = c(300, 750),
     xlab = "", ylab = "", pch = 1, cex = cex.size, main = "")
legend("topleft", expression(bold(H*alpha~"light curve")), bty = "n")
mtext(expression(bold(H*alpha~"flux")), side = 2, line = 2)
mtext("Observation time (MJD)", side = 1, line = 2.3)
arrows(ha[, 1], ha[, 2] - ha[, 3], ha[, 1], ha[, 2] + ha[, 3], 
       angle = 90, code = 3, length = 0.02, col = "gray70")

############### Fixed transfer-function width (top-hat) and standard deviation (Gaussian).

fixed_width_G  <- 5
fixed_width_TH <- sqrt(12) * fixed_width_G

############### Number of parallel cores.

ncores_default <- max(1L, parallel::detectCores() - 1L)

############### Read and prepare the data

tc <- as.numeric(cont[, 1])
yc <- as.numeric(cont[, 2])
ec <- as.numeric(cont[, 3])

tl <- as.numeric(ha[, 1])
yl <- as.numeric(ha[, 2])
el <- as.numeric(ha[, 3])

y <- c(yc, yl)

nc <- length(yc)
nl <- length(yl)
N  <- nc + nl

############### ============================================================
############### Stationary DRW covariance and PSD implied by
############### dZ(t) = -(1/tau) Z(t) dt + sigma dW(t)
############### sigma: diffusion coefficient / short-timescale amplitude
############### tau:   characteristic timescale
############### Var{Z(t)} = sigma^2 * tau / 2
############### ============================================================

K_drw <- function(t1, t2, sigma, tau) {
  if (!is.finite(sigma) || sigma <= 0) stop("sigma must be positive.")
  if (!is.finite(tau) || tau <= 0) stop("tau must be positive.")
  lag <- abs(outer(t1, t2, "-"))
  (sigma^2 * tau / 2) * exp(-lag / tau)
}

S_drw <- function(omega, sigma, tau) {
  sigma^2 * tau^2 / (1 + omega^2 * tau^2)
}

############### Transfer-function quadrature

quad_rule <- function(type = c("TH", "G"), tau0, width, nq = 31L,
                      gaussian_limit = 5) {
  type <- match.arg(type)
  nq <- as.integer(nq)
  if (nq < 3L) stop("nq must be at least 3.")
  if (!is.finite(width) || width <= 0) stop("width must be positive.")

  if (type == "TH") {
    u <- seq(tau0 - width / 2, tau0 + width / 2, length.out = nq)
    weights <- rep(1, nq)
    weights[c(1L, nq)] <- 0.5
    weights <- weights / sum(weights)
  } else {
    u <- seq(tau0 - gaussian_limit * width,
             tau0 + gaussian_limit * width,
             length.out = nq)
    weights <- dnorm(u, mean = tau0, sd = width)
    weights[c(1L, nq)] <- weights[c(1L, nq)] / 2
    weights <- weights / sum(weights)
  }
  list(u = u, weights = weights)
}

############### Matrix applying quadrature weights at multiple time points

make_W <- function(n_time, weights) {
  nq <- length(weights)
  W <- matrix(0, nrow = n_time, ncol = n_time * nq)
  for (i in seq_len(n_time)) {
    idx <- ((i - 1L) * nq + 1L):(i * nq)
    W[i, idx] <- weights
  }
  W
}

shifted_times <- function(times, delays) {
  as.vector(t(outer(times, delays, "-")))
}


############### Joint covariance of continuum and emission-line observations

build_joint_covariance <- function(par, type = c("TH", "G"), nq = 31L,
                                   width_TH = fixed_width_TH,
                                   width_G = fixed_width_G,
                                   jitter = 1e-8) {
  type <- match.arg(type)

  mu_c  <- unname(par["mu_c"])
  sigma <- exp(unname(par["log_sigma"]))
  tau   <- exp(unname(par["log_tau"]))
  mu_l  <- unname(par["mu_l"])
  alpha <- exp(unname(par["log_alpha"]))
  tau0  <- unname(par["tau0"])
  width <- if (type == "TH") width_TH else width_G

  q <- quad_rule(type, tau0, width, nq)
  u <- q$u; weights <- q$weights
  shifted_line <- shifted_times(tl, u)
  W_line <- make_W(nl, weights)

  Kcc <- K_drw(tc, tc, sigma, tau)
  Kcl <- alpha * K_drw(tc, shifted_line, sigma, tau) %*% t(W_line)
  Kll <- alpha^2 * W_line %*%
    K_drw(shifted_line, shifted_line, sigma, tau) %*% t(W_line)

  Kobs <- rbind(cbind(Kcc, Kcl), cbind(t(Kcl), Kll))
  diag(Kobs) <- diag(Kobs) + c(ec^2, el^2) + jitter

  list(
    K = Kobs,
    mean = c(rep(mu_c, nc), rep(mu_l, nl)),
    sigma = sigma, tau = tau, alpha = alpha, tau0 = tau0,
    width = width, delays = u, weights = weights,
    W_line = W_line, shifted_line = shifted_line
  )
}

############### Negative log-likelihood

neg_loglik_fixed <- function(par, type = c("TH", "G"), nq = 31L,
                             width_TH = fixed_width_TH,
                             width_G = fixed_width_G,
                             penalty = 1e100) {
  type <- match.arg(type)
  obj <- try(build_joint_covariance(par, type, nq, width_TH, width_G),
             silent = TRUE)
  if (inherits(obj, "try-error")) return(penalty)
  R <- try(chol(obj$K), silent = TRUE)
  if (inherits(R, "try-error")) return(penalty)
  residual <- y - obj$mean
  z <- forwardsolve(t(R), residual)
  nll <- 0.5 * sum(z^2) + sum(log(diag(R))) + 0.5 * N * log(2 * pi)
  if (!is.finite(nll)) penalty else nll
}

############### Starting values and bounds

safe_sd <- function(x) {
  ans <- sd(x)
  if (!is.finite(ans) || ans <= 0) 1 else ans
}

yc_sd <- safe_sd(yc); yl_sd <- safe_sd(yl)
baseline <- max(c(tc, tl)) - min(c(tc, tl))
if (!is.finite(baseline) || baseline <= 0) baseline <- 1000

tau_upper <- max(5000, 2 * baseline)
lag_upper <- min(500, max(300, baseline / 2))

lower_fixed <- c(
  mu_c = min(yc) - 3 * yc_sd,
  log_sigma = log(1e-4),
  log_tau = log(1),
  mu_l = min(yl) - 3 * yl_sd,
  log_alpha = log(1e-4),
  tau0 = 0
)

upper_fixed <- c(
  mu_c = max(yc) + 3 * yc_sd,
  log_sigma = log(100),
  log_tau = log(tau_upper),
  mu_l = max(yl) + 3 * yl_sd,
  log_alpha = log(1e4),
  tau0 = lag_upper
)

make_start_fixed <- function() {
  tau_start <- exp(runif(1, log(30), log(min(1000, tau_upper))))
  sigma_center <- yc_sd / sqrt(tau_start / 2)
  sigma_start <- exp(runif(
    1,
    log(max(1e-4, 0.3 * sigma_center)),
    log(max(1e-3, 3 * sigma_center))
  ))
  alpha_center <- yl_sd / max(yc_sd, 1e-6)
  alpha_start <- exp(runif(
    1,
    log(max(1e-4, 0.2 * alpha_center)),
    log(max(1e-3, 5 * alpha_center))
  ))

  start <- c(
    mu_c = mean(yc),
    log_sigma = log(sigma_start),
    log_tau = log(tau_start),
    mu_l = mean(yl),
    log_alpha = log(alpha_start),
    tau0 = runif(1, 1, lag_upper)
  )
  pmin(pmax(start, lower_fixed + 1e-8), upper_fixed - 1e-8)
}

############### Multi-start MLE

fit_one_fixed <- function(type = c("TH", "G"), nstart = 30L,
                          nq = 21L, final_nq = max(31L, nq),
                          ncores = ncores_default,
                          maxit = 1500L, final_maxit = 3000L,
                          width_TH = fixed_width_TH,
                          width_G = fixed_width_G,
                          seed = 2026,
                          final_hessian = TRUE) {
  type <- match.arg(type)
  set.seed(seed)
  starts <- replicate(nstart, make_start_fixed(), simplify = FALSE)
  overall_start <- Sys.time()

  cat("\n============================================================\n")
  cat("Model:", type, "| starts:", nstart, "| nq:", nq,
      "| final nq:", final_nq, "\n")
  cat("Started:", format(overall_start), "\n")
  cat("============================================================\n")

  fit_fun <- function(start, id) {
    run_start <- Sys.time()
    ans <- tryCatch(
      optim(
        par = start, fn = neg_loglik_fixed,
        type = type, nq = nq,
        width_TH = width_TH, width_G = width_G,
        method = "L-BFGS-B", lower = lower_fixed, upper = upper_fixed,
        control = list(maxit = maxit, factr = 1e7, pgtol = 1e-7)
      ),
      error = function(e) list(
        par = start, value = Inf, convergence = 999L,
        message = conditionMessage(e)
      )
    )
    run_end <- Sys.time()
    ans$start_id <- id
    ans$elapsed_minutes <- as.numeric(difftime(run_end, run_start, units = "mins"))
    cat(sprintf("[%3d/%3d] %s -> %s | %.2f min | nll=%s | conv=%d\n",
                id, nstart,
                format(run_start, "%Y-%m-%d %H:%M:%S"),
                format(run_end, "%Y-%m-%d %H:%M:%S"),
                ans$elapsed_minutes,
                if (is.finite(ans$value)) sprintf("%.6f", ans$value) else "Inf",
                ans$convergence))
    ans
  }

  ids <- seq_along(starts)
  if (.Platform$OS.type == "windows" || ncores <= 1L) {
    fits <- lapply(ids, function(i) fit_fun(starts[[i]], i))
  } else {
    fits <- parallel::mclapply(
      ids, function(i) fit_fun(starts[[i]], i),
      mc.cores = ncores, mc.preschedule = FALSE
    )
  }

  vals <- vapply(fits, function(z) {
    if (is.null(z$value) || !is.finite(z$value)) Inf else z$value
  }, numeric(1))
  if (all(!is.finite(vals))) stop("All optimization starts failed.")

  best_id <- which.min(vals)
  best <- fits[[best_id]]
  cat("\nBest initial start:", best_id, "| nll:", best$value, "\n")

  refined <- optim(
    par = best$par, fn = neg_loglik_fixed,
    type = type, nq = final_nq,
    width_TH = width_TH, width_G = width_G,
    method = "L-BFGS-B", lower = lower_fixed, upper = upper_fixed,
    hessian = final_hessian,
    control = list(maxit = final_maxit, factr = 1e5, pgtol = 1e-8)
  )

  p <- refined$par
  width <- if (type == "TH") width_TH else width_G
  estimates <- c(
    mu_c = unname(p["mu_c"]),
    sigma = exp(unname(p["log_sigma"])),
    tau = exp(unname(p["log_tau"])),
    mu_l = unname(p["mu_l"]),
    alpha = exp(unname(p["log_alpha"])),
    tau0 = unname(p["tau0"]),
    width = width
  )

  logLik_value <- -refined$value
  AIC_value <- -2 * logLik_value + 2 * length(p)
  overall_end <- Sys.time()

  cat("Final convergence:", refined$convergence, "\n")
  cat("Final logLik     :", logLik_value, "\n")
  cat("Final AIC        :", AIC_value, "\n")
  cat("Finished         :", format(overall_end), "\n\n")

  list(
    type = type, fit = refined, estimates = estimates,
    logLik = logLik_value, AIC = AIC_value,
    fixed_width = width, best_start = best_id,
    all_fits = fits, initial_objective_values = vals,
    final_nq = final_nq,
    timing = list(start = overall_start, end = overall_end)
  )
}

############### Hessian diagnostics and SEs

hessian_diagnostics_fixed <- function(fit_object) {
  H <- fit_object$fit$hessian
  if (is.null(H)) return(list(valid = FALSE, message = "No Hessian stored."))
  H <- (H + t(H)) / 2
  eig <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  valid <- all(is.finite(eig)) && min(eig) > 0
  list(
    valid = valid,
    message = if (valid) "Hessian is positive definite." else
      "Hessian is not positive definite; SEs are unreliable.",
    eigenvalues = eig,
    minimum_eigenvalue = min(eig),
    maximum_eigenvalue = max(eig),
    nonpositive = sum(eig <= 0),
    condition_number = if (valid) max(eig) / min(eig) else Inf
  )
}

extract_estimates_and_se_fixed <- function(fit_object) {
  info <- hessian_diagnostics_fixed(fit_object)
  p <- fit_object$fit$par
  est <- c(
    mu_c = unname(p["mu_c"]),
    sigma = exp(unname(p["log_sigma"])),
    tau = exp(unname(p["log_tau"])),
    mu_l = unname(p["mu_l"]),
    alpha = exp(unname(p["log_alpha"])),
    tau0 = unname(p["tau0"])
  )
  if (!info$valid) {
    warning(info$message)
    return(data.frame(parameter = names(est), estimate = as.numeric(est), SE = NA_real_))
  }
  H <- (fit_object$fit$hessian + t(fit_object$fit$hessian)) / 2
  se_int <- sqrt(diag(solve(H)))
  se <- c(
    mu_c = unname(se_int["mu_c"]),
    sigma = est["sigma"] * unname(se_int["log_sigma"]),
    tau = est["tau"] * unname(se_int["log_tau"]),
    mu_l = unname(se_int["mu_l"]),
    alpha = est["alpha"] * unname(se_int["log_alpha"]),
    tau0 = unname(se_int["tau0"])
  )
  data.frame(parameter = names(est), estimate = as.numeric(est), SE = as.numeric(se))
}

compare_best_starts <- function(fit_object, n = 10L) {
  vals <- sort(fit_object$initial_objective_values[
    is.finite(fit_object$initial_objective_values)
  ])
  vals <- head(vals, n)
  data.frame(rank = seq_along(vals), negative_log_likelihood = vals,
             difference_from_best = vals - min(vals))
}

############### Conditional prediction

predict_rm_fixed <- function(fit_or_par, type = c("TH", "G"),
                             t_grid, nq = 41L,
                             width_TH = fixed_width_TH,
                             width_G = fixed_width_G,
                             jitter = 1e-8) {
  type <- match.arg(type)
  par <- if (is.list(fit_or_par) && !is.null(fit_or_par$estimates)) {
    fit_or_par$estimates
  } else fit_or_par

  width <- if ("width" %in% names(par)) unname(par["width"]) else
    if (type == "TH") width_TH else width_G

  q <- quad_rule(type, unname(par["tau0"]), width, nq)
  u <- q$u; weights <- q$weights
  shifted_line <- shifted_times(tl, u)
  W_line <- make_W(nl, weights)

  sigma <- unname(par["sigma"])
  tau <- unname(par["tau"])
  alpha <- unname(par["alpha"])

  Kcc <- K_drw(tc, tc, sigma, tau)
  Kcl <- alpha * K_drw(tc, shifted_line, sigma, tau) %*% t(W_line)
  Kll <- alpha^2 * W_line %*%
    K_drw(shifted_line, shifted_line, sigma, tau) %*% t(W_line)
  Kobs <- rbind(cbind(Kcc, Kcl), cbind(t(Kcl), Kll))
  diag(Kobs) <- diag(Kobs) + c(ec^2, el^2) + jitter

  mean_obs <- c(rep(unname(par["mu_c"]), nc), rep(unname(par["mu_l"]), nl))
  R <- chol(Kobs)
  Kinv_res <- backsolve(R, forwardsolve(t(R), y - mean_obs))

  Kstar_c <- cbind(
    K_drw(t_grid, tc, sigma, tau),
    alpha * K_drw(t_grid, shifted_line, sigma, tau) %*% t(W_line)
  )
  continuum_mean <- as.numeric(unname(par["mu_c"]) + Kstar_c %*% Kinv_res)
  Vc <- forwardsolve(t(R), t(Kstar_c))
  continuum_var <- pmax(sigma^2 * tau / 2 - colSums(Vc^2), 0)

  shifted_grid <- shifted_times(t_grid, u)
  W_grid <- make_W(length(t_grid), weights)
  Kstar_l_cont <- alpha * W_grid %*% K_drw(shifted_grid, tc, sigma, tau)
  Kstar_l_line <- alpha^2 * W_grid %*%
    K_drw(shifted_grid, shifted_line, sigma, tau) %*% t(W_line)
  Kstar_l <- cbind(Kstar_l_cont, Kstar_l_line)
  line_mean <- as.numeric(unname(par["mu_l"]) + Kstar_l %*% Kinv_res)
  Vl <- forwardsolve(t(R), t(Kstar_l))
  Kuu <- K_drw(u, u, sigma, tau)
  line_prior_var <- as.numeric(alpha^2 * t(weights) %*% Kuu %*% weights)
  line_var <- pmax(line_prior_var - colSums(Vl^2), 0)

  list(
    time = t_grid,
    continuum_mean = continuum_mean,
    continuum_sd = sqrt(continuum_var),
    line_mean = line_mean,
    line_sd = sqrt(line_var)
  )
}

############### Figure 5: transfer functions

plot_transfer_functions <- function(fit_TH, fit_G, outfile = NULL,
                                    ngrid = 2000L, xlim = NULL) {
  par_TH <- if (is.list(fit_TH)) fit_TH$estimates else fit_TH
  par_G <- if (is.list(fit_G)) fit_G$estimates else fit_G
  tau_TH <- unname(par_TH["tau0"]); w_TH <- unname(par_TH["width"])
  tau_G <- unname(par_G["tau0"]); omega_G <- unname(par_G["width"])

  u_TH <- seq(tau_TH - 0.8 * w_TH, tau_TH + 0.8 * w_TH, length.out = ngrid)
  psi_TH <- ifelse(abs(u_TH - tau_TH) <= w_TH / 2, 1 / w_TH, 0)
  psi_TH <- psi_TH / max(psi_TH)
  u_G <- seq(tau_G - 4 * omega_G, tau_G + 4 * omega_G, length.out = ngrid)
  psi_G <- dnorm(u_G, tau_G, omega_G); psi_G <- psi_G / max(psi_G)
  if (is.null(xlim)) xlim <- range(c(u_TH, u_G))

  if (!is.null(outfile)) { pdf(outfile, 7.2, 7.2); on.exit(dev.off(), add = TRUE) }
  old <- par(no.readonly = TRUE); on.exit(par(old), add = TRUE)
  par(mfrow = c(2, 1), mar = c(4.2, 4.3, 2.2, 1))

  plot(u_TH, psi_TH, type = "l", lwd = 2, xlim = c(115, 165), ylim = c(0, 1.05),
       xlab = "", ylab = "",
       main = "Top-hat response")
  abline(v = tau_TH, lty = 3)
  mtext("Lag (days)", side = 1, line = 2, cex = 1.1)
  mtext(expression(Psi(u) / max(Psi(u))), side = 2, line = 2, cex = 1.1)

  legend("topright", c(sprintf("Lag = %.2f days", tau_TH),
                       sprintf("Width = %.2f days", w_TH)), bty = "n")

  plot(u_G, psi_G, type = "l", lwd = 2, xlim = c(115, 165), ylim = c(0, 1.05),
       xlab = "", ylab = "",
       main = "")
  title("Gaussian response")
  mtext("Lag (days)", side = 1, line = 2, cex = 1.1)
  mtext(expression(Psi(u) / max(Psi(u))), side = 2, line = 2, cex = 1.1)
  abline(v = tau_G, lty = 3)
  legend("topright", c(sprintf("Lag = %.2f days", tau_G),
                       sprintf("SD = %.2f days", omega_G)), bty = "n")
}

############### Figure 6: conditional H-alpha reconstructions

plot_line_reconstructions <- function(fit_TH, fit_G, outfile = NULL,
                                      nq = 41L, ngrid = 800L,
                                      xlim = range(c(tc, tl))) {
  par_TH <- if (is.list(fit_TH)) fit_TH$estimates else fit_TH
  par_G <- if (is.list(fit_G)) fit_G$estimates else fit_G
  t_grid <- seq(xlim[1], xlim[2], length.out = ngrid)
  pred_TH <- predict_rm_fixed(par_TH, "TH", t_grid, nq)
  pred_G <- predict_rm_fixed(par_G, "G", t_grid, nq)

  ylim_line <- range(
    yl - el, yl + el,
    pred_TH$line_mean - 1.96 * pred_TH$line_sd,
    pred_TH$line_mean + 1.96 * pred_TH$line_sd,
    pred_G$line_mean - 1.96 * pred_G$line_sd,
    pred_G$line_mean + 1.96 * pred_G$line_sd,
    finite = TRUE
  )

  if (!is.null(outfile)) { pdf(outfile, 8, 7); on.exit(dev.off(), add = TRUE) }
  old <- par(no.readonly = TRUE); on.exit(par(old), add = TRUE)
  par(mfrow = c(2, 1), mar = c(4, 4.3, 2.2, 1))

  draw_panel <- function(pred, title_text) {
    plot(tl, yl, type = "n", xlim = xlim, ylim = ylim_line,
         xlab = "",
         ylab = "", main = title_text)
    mtext("Observation time (MJD)", side = 1, line = 2, cex = 1.1)
    mtext(expression(H * alpha ~ "flux"), side = 2, line = 2, cex = 1.1)
    polygon(c(pred$time, rev(pred$time)),
            c(pred$line_mean - 1.96 * pred$line_sd,
              rev(pred$line_mean + 1.96 * pred$line_sd)),
            border = NA, col = adjustcolor("gray75", alpha.f = 0.45))
    arrows(tl, yl - el, tl, yl + el, angle = 90, code = 3,
           length = 0.015, col = "gray60")
    points(tl, yl, pch = 16, cex = 0.65)
    lines(pred$time, pred$line_mean, lwd = 2)
  }

  draw_panel(pred_TH, sprintf("Top-hat: lag = %.2f d, width = %.2f d",
                              par_TH["tau0"], par_TH["width"]))
  draw_panel(pred_G, sprintf("Gaussian: lag = %.2f d, SD = %.2f d",
                             par_G["tau0"], par_G["width"]))
  invisible(list(top_hat = pred_TH, gaussian = pred_G))
}

############### Figure 7: frequency responses and filtered PSD ratio

sinc_safe <- function(x) {
  out <- rep(1, length(x)); idx <- abs(x) > 1e-12
  out[idx] <- sin(x[idx]) / x[idx]
  out
}

plot_frequency_domain <- function(fit_TH, fit_G, outfile = NULL,
                                  omega_min = 1e-5, omega_max = 2,
                                  ngrid = 2000L) {
  par_TH <- if (is.list(fit_TH)) fit_TH$estimates else fit_TH
  par_G <- if (is.list(fit_G)) fit_G$estimates else fit_G
  omega <- 10^seq(log10(omega_min), log10(omega_max), length.out = ngrid)

  response_TH <- sinc_safe(omega * unname(par_TH["width"]) / 2)^2
  response_G <- exp(-omega^2 * unname(par_G["width"])^2)
  SZ_TH <- S_drw(omega, unname(par_TH["sigma"]), unname(par_TH["tau"]))
  SZ_G <- S_drw(omega, unname(par_G["sigma"]), unname(par_G["tau"]))
  SL_TH <- unname(par_TH["alpha"])^2 * response_TH * SZ_TH
  SL_G <- unname(par_G["alpha"])^2 * response_G * SZ_G
  tiny <- .Machine$double.xmin
  log10_ratio <- log10(pmax(SL_TH, tiny)) - log10(pmax(SL_G, tiny))

  if (!is.null(outfile)) { pdf(outfile, 7.2, 7.3); on.exit(dev.off(), add = TRUE) }
  old <- par(no.readonly = TRUE); on.exit(par(old), add = TRUE)
  par(mfrow = c(2, 1), mar = c(4.2, 4.6, 1.4, 1))

  plot(omega, response_TH, type = "l", log = "x", lwd = 2,
       ylim = c(0, 1.02), xlab = "",
       ylab = "")
  mtext(expression(omega), side = 1, line = 2, cex = 1.1)
  mtext(expression(abs(hat(Psi)(omega))^2), side = 2, line = 2, cex = 1.1)

  lines(omega, response_G, lwd = 2, lty = 2)
  legend("bottomleft", c("Top-hat", "Gaussian"), lty = c(1, 2),
         lwd = 2, bty = "n")

  plot(omega, log10_ratio, type = "l", log = "x", lwd = 2,
       xlab = "",
       ylab = "")
  abline(h = 0, lty = 2)
  mtext(expression(omega), side = 1, line = 2, cex = 1.1)
  mtext(expression(log[10] * (S[l]^TH * (omega)/S[l]^G * (omega))), side = 2, line = 2, cex = 1.1)

  text(
    x = omega_min * 3,
    y = 4,
    labels = "equal filtered power",
    adj = 0
  )

  invisible(data.frame(
    omega = omega, response_TH = response_TH, response_G = response_G,
    latent_PSD_TH = SZ_TH, latent_PSD_G = SZ_G,
    line_PSD_TH = SL_TH, line_PSD_G = SL_G,
    log10_ratio = log10_ratio
  ))
}

############### Test fit

fit_TH_test <- fit_one_fixed("TH", nstart = 1, nq = 15, final_nq = 21,
                             ncores = min(3L, ncores_default),
                             maxit = 300, final_maxit = 800)

fit_TH_test$estimates

fit_G_test <- fit_one_fixed("G", nstart = 1, nq = 15, final_nq = 21,
                            ncores = min(3L, ncores_default),
                            maxit = 300, final_maxit = 800)

fit_G_test$estimates

############### Fit

fit_TH_fixed <- fit_one_fixed("TH", nstart = 30, nq = 21, final_nq = 41,
                              ncores = ncores_default)

fit_TH_fixed$estimates
extract_estimates_and_se_fixed(fit_TH_fixed)
hessian_diagnostics_fixed(fit_TH_fixed)
compare_best_starts(fit_TH_fixed)

fit_G_fixed <- fit_one_fixed("G", nstart = 30, nq = 21, final_nq = 41,
                             ncores = ncores_default)

fit_G_fixed$estimates
extract_estimates_and_se_fixed(fit_G_fixed)
hessian_diagnostics_fixed(fit_G_fixed)
compare_best_starts(fit_G_fixed)

############### Figures

plot_transfer_functions(fit_TH_fixed, fit_G_fixed)

plot_line_reconstructions(fit_TH_fixed, fit_G_fixed, nq = 51)

frequency_results <- plot_frequency_domain(
  fit_TH_fixed, fit_G_fixed,
  omega_max = 2
)
