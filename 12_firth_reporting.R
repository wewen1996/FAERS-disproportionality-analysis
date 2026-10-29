# Firth logistic regression, grouped binomial likelihood.
# Reference equations: https://search.r-project.org/CRAN/refmans/logistf/html/logistf-package.html
# No tuning parameter. The full information determinant is retained in constrained
# profile fits. Confidence limits invert 2 * (max logL* - profile logL*) = chi2(.95,1).
# This is a local implementation, NOT a call to the unavailable logistf package.
firth_state <- function(beta, X, events, totals) {
  eta <- drop(X %*% beta)
  mu <- plogis(eta)
  w <- totals * mu * (1 - mu)
  information <- crossprod(X, X * w)
  ch <- tryCatch(chol(information), error = function(e) NULL)
  if (is.null(ch)) return(NULL)
  inv <- chol2inv(ch)
  h <- w * rowSums((X %*% inv) * X)
  score <- drop(crossprod(X, events - totals * mu + h * (0.5 - mu)))
  loglik <- sum(events * eta - totals * (pmax(eta, 0) + log1p(exp(-abs(eta))))) + sum(log(diag(ch)))
  list(loglik = loglik, score = score, information = information, covariance = inv)
}

firth_optimize <- function(X, events, totals, initial = NULL, fixed = NULL, maxit = 300L, tolerance = 1e-6) {
  p <- ncol(X)
  beta <- if (is.null(initial)) c(qlogis((sum(events) + .5) / (sum(totals) + 1)), rep(0, p - 1L)) else initial
  free <- seq_len(p)
  if (!is.null(fixed)) { beta[fixed$index] <- fixed$value; free <- setdiff(free, fixed$index) }
  state <- firth_state(beta, X, events, totals)
  if (is.null(state)) stop("Singular Firth information at initialization")
  for (iteration in seq_len(maxit)) {
    if (max(abs(state$score[free])) < tolerance) break
    step <- drop(solve(state$information[free, free, drop = FALSE], state$score[free]))
    if (max(abs(step)) > 3) step <- step * 3 / max(abs(step))
    accepted <- FALSE
    for (half in 0:35) {
      candidate <- beta; candidate[free] <- candidate[free] + step / 2^half
      next_state <- firth_state(candidate, X, events, totals)
      if (!is.null(next_state) && next_state$loglik >= state$loglik - 1e-9) { accepted <- TRUE; break }
    }
    if (!accepted) stop("Firth line search failed")
    beta <- candidate; state <- next_state
  }
  if (max(abs(state$score[free])) > 1e-4) stop("Firth score did not converge")
  list(beta = beta, loglik = state$loglik, covariance = state$covariance,
       max_score = max(abs(state$score[free])), iterations = iteration)
}

firth_profile <- function(X, events, totals, exposure_column = "exposure") {
  qr_x <- qr(X, tol = 1e-9)
  keep <- sort(qr_x$pivot[seq_len(qr_x$rank)])
  removed <- setdiff(colnames(X), colnames(X)[keep])
  X <- X[, keep, drop = FALSE]
  if (!exposure_column %in% colnames(X)) stop("Exposure is not identifiable in the covariate design")
  j <- match(exposure_column, colnames(X))
  fit <- firth_optimize(X, events, totals)
  center <- fit$beta[j]
  cutoff <- qchisq(.95, 1)
  profile_value <- function(value) {
    constrained <- firth_optimize(X, events, totals, initial = fit$beta, fixed = list(index = j, value = value))
    2 * (fit$loglik - constrained$loglik)
  }
  root <- function(direction) {
    distance <- max(.5, sqrt(fit$covariance[j, j]))
    edge <- center + direction * distance
    for (attempt in 1:18) {
      if (profile_value(edge) >= cutoff) break
      distance <- distance * 1.5; edge <- center + direction * distance
    }
    if (profile_value(edge) < cutoff) stop("Could not bracket profile likelihood CI")
    uniroot(function(value) profile_value(value) - cutoff,
      interval = sort(c(center, edge)), tol = 1e-5)$root
  }
  low <- root(-1); high <- root(1)
  statistic <- max(0, profile_value(0))
  list(ROR = exp(center), low = exp(low), high = exp(high), p = pchisq(statistic, 1, lower.tail = FALSE),
    parameters = ncol(X), removed = paste(removed, collapse = ";"), fit = fit,
    endpoint_error = max(abs(c(profile_value(low), profile_value(high)) - cutoff)), X = X)
}

fit_adjusted_ror <- function(data) {
  # Four clinical/reporting covariates, eight slopes at most including exposure.
  # Age is already in YEARS in this cleaned data. A fixed age of 50 and a separate
  # missing indicator preserve reports without pretending missing age is observed.
  d <- data.table::copy(data)
  d[, `:=`(exposure = as.integer(EXPOSURE == "target"),
    age10 = data.table::fifelse(is.na(AGE_YEARS), 0, (AGE_YEARS - 50) / 10),
    age_missing = as.integer(is.na(AGE_YEARS)),
    sex_male = as.integer(SEX_GROUP == "M"), sex_missing = as.integer(SEX_GROUP == "Unknown"),
    calendar_time = calendar_year - 2022,
    reporter_nonhcp = as.integer(REPORTER_GROUP == "non_HCP"),
    reporter_missing = as.integer(REPORTER_GROUP == "Unknown"))]
  columns <- c("exposure", "age10", "age_missing", "sex_male", "sex_missing", "calendar_time", "reporter_nonhcp", "reporter_missing")
  grouped <- d[, .(events = sum(TARGET), totals = .N), by = columns]
  X <- cbind("(Intercept)" = 1, as.matrix(grouped[, ..columns]))
  fit <- firth_profile(X, grouped$events, grouped$totals)
  list(summary = data.table::data.table(adjusted_model = "Firth_logistic", adjusted_model_status = "validated_profile_CI",
    adjusted_n = nrow(d), adjusted_events = sum(d$TARGET), adjusted_ROR = fit$ROR,
    adjusted_ROR_low = fit$low, adjusted_ROR_high = fit$high, adjusted_p = fit$p,
    covariate_columns = fit$parameters - 1L, parameters = fit$parameters,
    events_per_parameter = sum(d$TARGET) / fit$parameters, removed_constant_or_alias_columns = fit$removed,
    max_score = fit$fit$max_score, iterations = fit$fit$iterations, profile_endpoint_error = fit$endpoint_error),
    grouped = grouped, coefficients = data.table::data.table(term = colnames(fit$X), beta = fit$fit$beta))
}
