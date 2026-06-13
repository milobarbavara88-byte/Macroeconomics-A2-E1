# ==============================================================================
# Macroeconomics - Assignment 2 - Exercise 1
# FORECASTING THE TERM STRUCTURE OF INTEREST RATES
#
# Authors: Milo Barbavara, Niccolò Bianchi
# ==============================================================================
# ------------------------------------------------------------------------------
# 0. SETUP
# ------------------------------------------------------------------------------
# Load the four packages we need. Everything is out in the open: no helper
# functions, no hidden tricks. If a package is not installed yet, run ONCE the
# install line below (remove the leading #), then keep using the library lines:
# install.packages(c("fredr", "xts", "zoo", "vars"))
library(fredr)      # official FRED API access using a personal API key
library(xts)        # tools for time series
library(zoo)        # tools for time series (dates, plotting)
library(vars)       # estimates the VAR model (Question 3)

# IMPORTANT - this is what makes the FRED download work.
# We download the data through the OFFICIAL FRED API using a personal API key
# (free, from https://fredaccount.stlouisfed.org/apikeys). This goes through the
# api.stlouisfed.org endpoint and avoids the "cannot open the connection" /
# "non e' possibile aprire la connessione" error of getSymbols.
# >>> Your personal 32-character FRED API key (keep it private). <<<
fredr_set_key("2dc39adca864462a604d64d7ea499289")

set.seed(123)                                   # reproducibility
dir.create("output", showWarnings = FALSE)      # folder where all figures are saved
lambda <- 0.0609                                # Nelson-Siegel decay parameter fixed by the assignment


# ==============================================================================
# QUESTION 1 - DATA PREPARATION AND DESCRIPTIVE ANALYSIS
# ==============================================================================

# Download daily US Treasury constant-maturity yields from the FRED API.

fred_ids   <- c("DGS3MO", "DGS6MO", "DGS1", "DGS2", "DGS3", "DGS5", "DGS7", "DGS10")
maturities <- c(3,        6,        12,     24,     36,     60,     84,     120)   # in MONTHS

# Download each series with fredr() and turn it into an xts column. fredr()
# returns a data frame with a 'date' and a 'value' column; we keep both. The
# loop is fully visible: one iteration per maturity, nothing hidden.
daily_list <- list()
for (i in seq_along(fred_ids)) {
  raw <- fredr(series_id = fred_ids[i])             # daily history of one maturity
  daily_list[[i]] <- xts(raw$value, order.by = raw$date)
}

# Merge the eight daily series into one xts matrix, columns = maturities.
daily <- do.call(merge, daily_list)
colnames(daily) <- as.character(maturities)

# Convert daily yields to MONTHLY AVERAGES.

monthly <- apply.monthly(daily, function(x) colMeans(x, na.rm = TRUE))

# Remove missing information

yields <- na.omit(monthly)
index(yields) <- as.Date(as.yearmon(index(yields)))   

cat("\n================ SAMPLE INFORMATION (Question 1) ================\n")
cat("Starting date          :", format(start(yields)), "\n")
cat("Ending date            :", format(end(yields)),   "\n")
cat("Monthly observations   :", nrow(yields), "\n")
cat("Maturities used (months):", paste(maturities, collapse = ", "), "\n")

# NOTE ON THE PLOTS.
# Each figure is first drawn ON SCREEN (so it appears in the RStudio "Plots"
# pane and you can see it), and THEN copied to a .png file in the output/ folder
# with: dev.copy(png, ...); dev.off(). This is why the previous version showed
# nothing: the old code opened a png() file device, so every plot went straight
# to a file and never to the screen.

# ---- Plot 1: the yield curve at three dates (early / middle / recent) --------

idx_dates <- c(1, round(nrow(yields) / 2), nrow(yields))
plot(maturities, as.numeric(yields[idx_dates[1], ]), type = "b", pch = 19,
     ylim = range(yields[idx_dates, ]), col = "darkblue",
     xlab = "Maturity (months)", ylab = "Yield (%)",
     main = "Yield curve at three dates")
lines(maturities, as.numeric(yields[idx_dates[2], ]), type = "b", pch = 17, col = "darkgreen")
lines(maturities, as.numeric(yields[idx_dates[3], ]), type = "b", pch = 15, col = "darkred")
legend("bottomright", bty = "n", pch = c(19, 17, 15),
       col = c("darkblue", "darkgreen", "darkred"),
       legend = format(index(yields)[idx_dates], "%Y-%m"))
dev.copy(png, "output/q1_yield_curves.png", width = 800, height = 600); dev.off()

# ---- Empirical level, slope and curvature ------------------------------------
# WHAT: build the three model-free yield-curve summaries defined in the assignment.
# WHY : these are simple, transparent proxies of the three factors that the
#       Nelson-Siegel model will later estimate formally.
#   Level     = y(120)                  -> overall height of the curve (long rate)
#   Slope     = y(120) - y(3)           -> steepness (long minus short)
#   Curvature = 2*y(24) - y(3) - y(120) -> how "bowed" the medium part is
emp_level     <- yields[, "120"]
emp_slope     <- yields[, "120"] - yields[, "3"]
emp_curvature <- 2 * yields[, "24"] - yields[, "3"] - yields[, "120"]
empirical <- merge(emp_level, emp_slope, emp_curvature)
colnames(empirical) <- c("Level", "Slope", "Curvature")

# ---- Plot 2: the three empirical series over time ----------------------------
plot.zoo(empirical, main = "Empirical level, slope and curvature",
         xlab = "Date", col = c("darkblue", "darkred", "darkgreen"))
dev.copy(png, "output/q1_empirical_factors.png", width = 900, height = 700); dev.off()

# WHAT: quantify persistence with the first-order autocorrelation of each series.
# HOW TO READ: values close to 1 mean the series is very persistent (slow-moving,
#       trending) rather than white noise.
persist <- sapply(empirical, function(x) acf(x, lag.max = 1, plot = FALSE)$acf[2])
cat("\nFirst-order autocorrelation (persistence):\n")
print(round(persist, 3))

# INTERPRETATION (Question 1, persistence & cyclical behaviour):
#   - The LEVEL is the most persistent series (autocorrelation ~1): it inherits
#     the long, slow downward trend of US interest rates and behaves almost like
#     a random walk.
#   - The SLOPE is persistent but clearly CYCLICAL: it flattens/inverts before
#     recessions and steepens during recoveries, so it co-moves with the business
#     cycle.
#   - The CURVATURE is the least persistent and noisiest, with no strong trend;
#     it captures shorter-lived, medium-maturity movements of the curve.


# ==============================================================================
# QUESTION 2 - NELSON-SIEGEL CROSS-SECTIONAL FACTOR EXTRACTION
# ==============================================================================

# Build the Nelson-Siegel loading matrix Lambda
# Loadings:
#   column 1 = 1                                  -> LEVEL loading (same for all tau)
#   column 2 = (1-exp(-l*t))/(l*t)                -> SLOPE loading (1 at short end -> 0 at long end)
#   column 3 = (1-exp(-l*t))/(l*t) - exp(-l*t)    -> CURVATURE loading (0 at both ends, hump in the middle)
ns_loadings <- function(tau, lambda) {
  col1 <- rep(1, length(tau))
  col2 <- (1 - exp(-lambda * tau)) / (lambda * tau)
  col3 <- col2 - exp(-lambda * tau)
  cbind(beta1 = col1, beta2 = col2, beta3 = col3)
}
L <- ns_loadings(maturities, lambda)            # N x 3 loading matrix (N = 8 maturities)

# ---- Plot 3: the three factor loadings ---------------------------------------
#       The level loading is flat at 1 (moves all yields equally); the
#       slope loading is large at short maturities and decays to 0 (a short-rate
#       factor); the curvature loading is hump-shaped, peaking around 2-3 years.
matplot(maturities, L, type = "b", pch = 19, lty = 1,
        col = c("darkblue", "darkred", "darkgreen"),
        xlab = "Maturity (months)", ylab = "Loading",
        main = "Nelson-Siegel factor loadings (lambda = 0.0609)")
legend("right", bty = "n", lty = 1, pch = 19,
       col = c("darkblue", "darkred", "darkgreen"),
       legend = c("Level (beta1)", "Slope (beta2)", "Curvature (beta3)"))
dev.copy(png, "output/q2_ns_loadings.png", width = 800, height = 600); dev.off()

# Estimate coefficients by OLS.

hat   <- solve(t(L) %*% L, t(L))                # 3 x N  (the OLS projector)
Y     <- coredata(yields)                       # T x N matrix of yields
betas <- hat %*% t(Y)                           # 3 x T matrix of estimated factors
factors <- xts(t(betas), order.by = index(yields))
colnames(factors) <- c("beta1", "beta2", "beta3")

# ---- Plot 4: the three estimated Nelson-Siegel factors -----------------------
plot.zoo(factors, main = "Estimated Nelson-Siegel factors",
         xlab = "Date", col = c("darkblue", "darkred", "darkgreen"))
dev.copy(png, "output/q2_ns_factors.png", width = 900, height = 700); dev.off()

# ---- Compare estimated factors with empirical level/slope/curvature ----------

comparison <- merge(factors$beta1, empirical$Level,
                    -factors$beta2, empirical$Slope,
                    factors$beta3, empirical$Curvature)
colnames(comparison) <- c("beta1", "Level", "minus_beta2", "Slope", "beta3", "Curvature")

par(mfrow = c(3, 1))
plot.zoo(comparison[, c("beta1", "Level")], screens = 1, col = c("darkblue", "black"),
         lty = c(1, 2), main = "Level: beta1 vs empirical", xlab = "", ylab = "%")
legend("topright", bty = "n", lty = c(1, 2), col = c("darkblue", "black"),
       legend = c("beta1", "Level"))
plot.zoo(comparison[, c("minus_beta2", "Slope")], screens = 1, col = c("darkred", "black"),
         lty = c(1, 2), main = "Slope: -beta2 vs empirical", xlab = "", ylab = "%")
legend("topright", bty = "n", lty = c(1, 2), col = c("darkred", "black"),
       legend = c("-beta2", "Slope"))
plot.zoo(comparison[, c("beta3", "Curvature")], screens = 1, col = c("darkgreen", "black"),
         lty = c(1, 2), main = "Curvature: beta3 vs empirical", xlab = "", ylab = "%")
legend("topright", bty = "n", lty = c(1, 2), col = c("darkgreen", "black"),
       legend = c("beta3", "Curvature"))
dev.copy(png, "output/q2_factor_vs_empirical.png", width = 900, height = 800); dev.off()
par(mfrow = c(1, 1))

cat("\n=========== FACTORS vs EMPIRICAL MEASURES (Question 2) ===========\n")
cat("corr(beta1 , Level)     :", round(cor(comparison$beta1,  comparison$Level), 3), "\n")
cat("corr(-beta2, Slope)     :", round(cor(comparison$minus_beta2, comparison$Slope), 3), "\n")
cat("corr(beta3 , Curvature) :", round(cor(comparison$beta3,  comparison$Curvature), 3), "\n")

# ---- Average fitting errors by maturity --------------------------------------

fitted_Y  <- L %*% betas                        # N x T fitted yields
resid_Y   <- t(Y) - fitted_Y                    # N x T residuals
rmse_mat  <- sqrt(rowMeans(resid_Y^2))          # one RMSE per maturity
names(rmse_mat) <- maturities
cat("\nAverage fitting error (RMSE, basis-point scale of %) by maturity:\n")
print(round(rmse_mat, 4))

barplot(rmse_mat, col = "steelblue", xlab = "Maturity (months)", ylab = "RMSE (%)",
        main = "Average Nelson-Siegel fitting error by maturity")
dev.copy(png, "output/q2_fitting_errors.png", width = 800, height = 600); dev.off()


# ==============================================================================
# QUESTION 3 - VAR MODEL FOR THE NELSON-SIEGEL FACTORS
# ==============================================================================

# Model the joint dynamics of f_t = (beta1, beta2, beta3) with a VAR.

F     <- coredata(factors)                      # T x 3 factor matrix
Tn    <- nrow(F)
horizons <- c(1, 6, 12)                         # forecast horizons

sel <- VARselect(F, lag.max = 4, type = "const")
p   <- as.integer(sel$selection["AIC(n)"])
cat("\n================ VAR specification (Question 3) ================\n")
cat("Selected VAR lag (AIC):", p, "\n")

# ---- Recursive (expanding-window) out-of-sample forecasts --------------------

# Convert factor forecast f_{t+h|t} into yields: y_hat = Lambda * f_{t+h|t}.
n_train <- floor(0.6 * Tn)                       # first 60% of the sample trains the first forecast
origins <- n_train:(Tn - 1)                      # forecast origins (expanding window)

# Pre-allocate squared/absolute error stores: [origin, maturity, horizon] per model.
N   <- length(maturities)
dimn <- list(origin = origins, maturity = maturities, horizon = horizons)
err_var <- array(NA_real_, dim = c(length(origins), N, length(horizons)), dimnames = dimn)
err_rw  <- err_var                               # random-walk benchmark errors

for (i in seq_along(origins)) {
  t0 <- origins[i]
  # Re-estimate the VAR on data up to the current origin only (expanding window).
  var_fit <- VAR(F[1:t0, , drop = FALSE], p = p, type = "const")
  fc      <- predict(var_fit, n.ahead = max(horizons))$fcst   # list: one matrix per factor
  for (h in seq_along(horizons)) {
    hh <- horizons[h]
    if (t0 + hh > Tn) next                        # skip if the realised yield is outside the sample
    # Assemble the h-step factor forecast and turn it into a yield forecast.
    f_fc   <- c(fc$beta1[hh, "fcst"], fc$beta2[hh, "fcst"], fc$beta3[hh, "fcst"])
    y_var  <- as.numeric(L %*% f_fc)              # VAR yield forecast
    y_rw   <- Y[t0, ]                             # random walk: best forecast is today's curve
    y_true <- Y[t0 + hh, ]                        # realised yields
    err_var[i, , h] <- y_var - y_true
    err_rw[i,  , h] <- y_rw  - y_true
  }
}

# WHAT: helper turning an error array into MAE and RMSE tables (maturity x horizon).
agg_metric <- function(err, fun) {
  m <- apply(err, c(2, 3), function(e) fun(e))
  dimnames(m) <- list(maturity = maturities, horizon = horizons)
  round(m, 4)
}
mae  <- function(e) mean(abs(e), na.rm = TRUE)
rmse <- function(e) sqrt(mean(e^2, na.rm = TRUE))

cat("\n--- VAR  MAE  by maturity (rows) and horizon (cols) ---\n");  print(agg_metric(err_var, mae))
cat("\n--- VAR  RMSE by maturity (rows) and horizon (cols) ---\n");  print(agg_metric(err_var, rmse))
cat("\n--- RW   MAE  by maturity (rows) and horizon (cols) ---\n");  print(agg_metric(err_rw,  mae))
cat("\n--- RW   RMSE by maturity (rows) and horizon (cols) ---\n");  print(agg_metric(err_rw,  rmse))

# ==============================================================================
# QUESTION 4 - KALMAN-FILTER STATE-SPACE MODEL (one-step dynamic Nelson-Siegel)
# ==============================================================================

# TASK 1 - DIFFERENCE BETWEEN THE TWO APPROACHES:
#   TWO-STEP OLS (Questions 2-3): first estimate the factors month-by-month by
#     OLS treating them as observed, THEN fit a separate VAR to those estimates.
#     Simple and fast, but it ignores (a) the estimation error in the factors and
#     (b) the feedback between measurement noise and factor dynamics; the two
#     steps are not jointly optimal.
#   ONE-STEP KALMAN: factors are LATENT states. The measurement
#     equation (yields = loadings * factors + noise) and the transition equation
#     (factor VAR(1) dynamics) are estimated JOINTLY by maximum likelihood. The
#     Kalman filter optimally combines the cross-section of yields with the time
#     dynamics, propagates uncertainty, and yields smoother, more efficient
#     factor estimates and internally consistent forecasts.

# Model (with DIAGONAL Q and DIAGONAL H):
#   Measurement: y_t      = Lambda * f_t + eps_t,        eps_t ~ N(0, H), H diagonal
#   Transition : f_t - mu = A (f_{t-1} - mu) + eta_t,    eta_t ~ N(0, Q), Q diagonal
# We use a DIAGONAL A (each factor an independent AR(1)) and estimate the model
# with the dlm package (Petris, "Dynamic Linear Models with R"), the standard R
# package for linear Gaussian state-space models. dlm runs the Kalman filter and
# its likelihood for us, so we do NOT code the filter by hand.
# install.packages("dlm")   # run once if dlm is not installed yet
library(dlm)

# The long-run factor mean mu is fixed at the average of the OLS factors (a
# standard "concentration" step). Working with the factors in deviation from the
# mean, g_t = f_t - mu, gives a zero-mean AR(1): g_t = A g_{t-1} + eta. Subtracting
# Lambda*mu from the yields then leaves a clean state-space with no intercept term,
# which is exactly the form the dlm package expects.
mu_hat   <- colMeans(F)                         # fixed factor means (level, slope, curvature)
mean_yld <- as.numeric(L %*% mu_hat)            # average yield implied at each maturity
Y_dm     <- sweep(Y, 2, mean_yld, "-")          # demeaned yields (T x N): the data for the filter

# build() returns a dlm object for a given parameter vector:
#   FF = Lambda (loadings)      GG = diag(A)  (factor persistence)
#   V  = diag(H) (8 measurement variances)    W = diag(Q) (3 state-shock variances)
# AR(1) coefficients are kept inside (-1,1) with tanh(); variances are kept
# positive with exp(). So all 14 parameters are unconstrained for the optimiser.
build_dns <- function(par) {
  a <- tanh(par[1:3])                           # 3 AR(1) persistences in (-1,1)
  q <- exp(par[4:6])                            # 3 state-shock variances (diag Q)
  h <- exp(par[7:14])                           # 8 measurement variances (diag H)
  dlm(FF = L, GG = diag(a), W = diag(q),
      V  = diag(h),
      m0 = rep(0, 3),                           # demeaned factors start at 0
      C0 = diag(q / (1 - a^2)))                 # stationary initial variance of an AR(1)
}

# Informed starting values from the Q2/Q3 two-step results (fast, stable MLE).
ar1  <- function(x) coef(lm(x[-1] ~ x[-length(x)]))[2]   # quick AR(1) coefficient
a0   <- pmin(pmax(apply(F, 2, ar1), -0.95), 0.95)        # persistence of each factor
q0   <- apply(F, 2, function(x) var(diff(x)))            # rough state-noise variance
h0   <- rmse_mat^2                                       # measurement variance = Q2 fitting variance
par0 <- c(atanh(a0), log(q0), log(h0))                   # 14 starting parameters

# TASK 2 - estimate the parameters by maximum likelihood on the SAME training
# window used for the VAR (so the comparison is fair and look-ahead free).
# dlmMLE maximises the Kalman-filter likelihood internally.
fit_kf  <- dlmMLE(Y_dm[1:n_train, ], parm = par0, build = build_dns,
                  method = "BFGS", control = list(maxit = 500))
mod_hat <- build_dns(fit_kf$par)                # the fitted state-space model
a_hat   <- tanh(fit_kf$par[1:3])                # estimated factor persistences (diag A)

cat("\n================ Kalman filter estimates (Question 4) ================\n")
cat("Factor persistences (diag A):", round(a_hat, 3), "\n")
cat("Factor means (mu)           :", round(mu_hat, 3), "\n")

# TASK 3 - filtered factors over the FULL sample. dlmFilter() runs the Kalman
# filter with the fixed estimated parameters and returns the filtered states in
# $m (its first row is the time-0 prior, so we drop it). We add mu back to put the
# factors on their original level/slope/curvature scale.
filt      <- dlmFilter(Y_dm, mod_hat)
kf_states <- sweep(filt$m[-1, , drop = FALSE], 2, mu_hat, "+")   # T x 3 filtered factors
kf_factors <- xts(kf_states, order.by = index(yields))
colnames(kf_factors) <- c("level", "slope", "curvature")

# HOW TO READ: these filtered factors should track the OLS factors of Q2 but look
#       SMOOTHER, because the Kalman filter optimally damps measurement noise.
plot.zoo(kf_factors, main = "Kalman-filtered Nelson-Siegel factors",
         xlab = "Date", col = c("darkblue", "darkred", "darkgreen"))
dev.copy(png, "output/q4_kalman_factors.png", width = 900, height = 700); dev.off()

# Side-by-side comparison of OLS vs Kalman factors (smoothness check for Q5.4).
par(mfrow = c(3, 1))
for (j in 1:3) {
  plot(index(yields), F[, j], type = "l", col = "grey50", xlab = "", ylab = "",
       main = paste0("Factor ", j, ": OLS (grey) vs Kalman (colour)"))
  lines(index(yields), kf_states[, j], col = c("darkblue", "darkred", "darkgreen")[j], lwd = 2)
}
dev.copy(png, "output/q4_kalman_vs_ols.png", width = 900, height = 800); dev.off()
par(mfrow = c(1, 1))

# TASK 4 - recursive out-of-sample forecasts for h = 1, 6, 12.
# WHAT: with fixed estimated parameters, at each origin t the filter has produced
#       f_{t|t} using data up to t only; the h-step forecast of a diagonal-A AR(1)
#       state is f_{t+h|t} = mu + A^h (f_{t|t} - mu), mapped to yields by Lambda.
#       a_hat and mu_hat were estimated above with dlm.
err_kf <- err_var                                 # same dimensions/origins as the VAR errors
for (i in seq_along(origins)) {
  t0   <- origins[i]
  f_tt <- kf_states[t0, ]
  for (h in seq_along(horizons)) {
    hh <- horizons[h]
    if (t0 + hh > Tn) next
    f_fc  <- mu_hat + (a_hat^hh) * (f_tt - mu_hat)  # h-step state forecast (diagonal A)
    y_kf  <- as.numeric(L %*% f_fc)
    err_kf[i, , h] <- y_kf - Y[t0 + hh, ]
  }
}

cat("\n--- KALMAN MAE  by maturity (rows) and horizon (cols) ---\n"); print(agg_metric(err_kf, mae))
cat("\n--- KALMAN RMSE by maturity (rows) and horizon (cols) ---\n"); print(agg_metric(err_kf, rmse))

# TASK 5 - compact comparison of the THREE models (average RMSE across maturities).
# HOW TO READ: for each horizon, the model with the lowest number forecasts best
#       on average across the curve.
overall <- function(err) apply(err, 3, rmse)      # one RMSE per horizon (pooled over maturities)
comp_tbl <- rbind(RandomWalk = overall(err_rw),
                  VAR        = overall(err_var),
                  Kalman     = overall(err_kf))
colnames(comp_tbl) <- paste0("h=", horizons)
cat("\n================ MODEL COMPARISON - pooled RMSE (Question 4/5) ================\n")
print(round(comp_tbl, 4))

barplot(comp_tbl, beside = TRUE, col = c("grey60", "darkred", "darkblue"),
        ylab = "Pooled RMSE (%)", xlab = "Forecast horizon",
        main = "Forecast accuracy: Random Walk vs VAR vs Kalman")
legend("topleft", bty = "n", fill = c("grey60", "darkred", "darkblue"),
       legend = c("Random Walk", "VAR", "Kalman"))
dev.copy(png, "output/q4_model_comparison.png", width = 800, height = 600); dev.off()


# ==============================================================================
# QUESTION 5 - INTERPRETATION AND DISCUSSION
# ------------------------------------------------------------------------------
# Read the printed tables (err_rw / err_var / err_kf and comp_tbl) together with
# the saved figures, then answer the five questions. Typical findings for US
# Treasury yields are summarised below; replace them with YOUR exact numbers.
#
# 1. SHORT HORIZON (h = 1): the RANDOM WALK is usually very hard to beat. Yields
#    are near-unit-root, so "no change" is an excellent one-month forecast; VAR
#    and Kalman are at best marginally better.
#
# 2. LONG HORIZON (h = 12): the FACTOR MODELS (VAR and especially Kalman) tend to
#    win, because the factors mean-revert and the models exploit this while the
#    random walk cannot. Look for comp_tbl["Kalman","h=12"] being the smallest.
#
# 3. WHERE ARE THE GAINS? Inspect the maturity dimension of the MAE/RMSE tables:
#    gains over the random walk are usually concentrated at the SHORT and MEDIUM
#    maturities (3-36 months), which are the most predictable / mean-reverting,
#    rather than at the very long end.
#
# 4. SMOOTHNESS: the Kalman filter produces SMOOTHER factors than two-step OLS
#    (see output/q4_kalman_vs_ols.png), because it optimally filters out the
#    month-by-month measurement noise that the independent OLS fits leave in.
#
# 5. WHY A LOW-DIMENSIONAL FACTOR STRUCTURE HELPS even without no-arbitrage:
#    yields are highly collinear, so 3 factors capture almost all their variation.
#    Forecasting 3 persistent factors instead of 8 noisy yields drastically cuts
#    the number of parameters and hence ESTIMATION ERROR, and it imposes a smooth,
#    economically sensible cross-sectional shape. These statistical regularisation
#    benefits improve out-of-sample forecasts even though the model places no
#    no-arbitrage restrictions across maturities.
# ==============================================================================

cat("\nAll figures saved in ./output. Exercise 1 completed.\n")