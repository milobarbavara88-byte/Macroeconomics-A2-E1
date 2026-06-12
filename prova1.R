# ==============================================================================
# Macroeconomics - Assignment 2 - Exercise 1
# FORECASTING THE TERM STRUCTURE OF INTEREST RATES
#
# Authors: Milo Barbavara, Niccolò Bianchi
# ==============================================================================
# ------------------------------------------------------------------------------
# 0. SETUP
# ------------------------------------------------------------------------------
ensure_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}
invisible(lapply(c("quantmod", "xts", "zoo", "vars"), ensure_pkg))

set.seed(123)                                   
dir.create("output", showWarnings = FALSE)      
lambda <- 0.0609                                # Nelson-Siegel decay parameter fixed by the assignment


# ==============================================================================
# QUESTION 1 - DATA PREPARATION AND DESCRIPTIVE ANALYSIS
# ==============================================================================

# Download daily US Treasury constant-maturity yields from FRED.

fred_ids   <- c("DGS3MO", "DGS6MO", "DGS1", "DGS2", "DGS3", "DGS5", "DGS7", "DGS10")
maturities <- c(3,        6,        12,     24,     36,     60,     84,     120)   # in MONTHS
getSymbols(fred_ids, src = "FRED")              # creates one xts object per id in the workspace

# Merge the eight daily series into one xts matrix, columns = maturities.
daily <- do.call(merge, lapply(fred_ids, get))
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

# ---- Plot 1: the yield curve at three dates (early / middle / recent) --------

idx_dates <- c(1, round(nrow(yields) / 2), nrow(yields))
png("output/q1_yield_curves.png", width = 800, height = 600)
plot(maturities, as.numeric(yields[idx_dates[1], ]), type = "b", pch = 19,
     ylim = range(yields[idx_dates, ]), col = "darkblue",
     xlab = "Maturity (months)", ylab = "Yield (%)",
     main = "Yield curve at three dates")
lines(maturities, as.numeric(yields[idx_dates[2], ]), type = "b", pch = 17, col = "darkgreen")
lines(maturities, as.numeric(yields[idx_dates[3], ]), type = "b", pch = 15, col = "darkred")
legend("bottomright", bty = "n", pch = c(19, 17, 15),
       col = c("darkblue", "darkgreen", "darkred"),
       legend = format(index(yields)[idx_dates], "%Y-%m"))
dev.off()

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
png("output/q1_empirical_factors.png", width = 900, height = 700)
plot.zoo(empirical, main = "Empirical level, slope and curvature",
         xlab = "Date", col = c("darkblue", "darkred", "darkgreen"))
dev.off()

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

# WHAT: build the Nelson-Siegel loading matrix Lambda(lambda).
# WHY : it maps the three factors (beta1, beta2, beta3) into the yield of every
#       maturity. Because lambda is FIXED, the loadings are constant over time,
#       which turns factor estimation into a simple linear regression each month.
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
# HOW TO READ: the level loading is flat at 1 (moves all yields equally); the
#       slope loading is large at short maturities and decays to 0 (a short-rate
#       factor); the curvature loading is hump-shaped, peaking around 2-3 years.
png("output/q2_ns_loadings.png", width = 800, height = 600)
matplot(maturities, L, type = "b", pch = 19, lty = 1,
        col = c("darkblue", "darkred", "darkgreen"),
        xlab = "Maturity (months)", ylab = "Loading",
        main = "Nelson-Siegel factor loadings (lambda = 0.0609)")
legend("right", bty = "n", lty = 1, pch = 19,
       col = c("darkblue", "darkred", "darkgreen"),
       legend = c("Level (beta1)", "Slope (beta2)", "Curvature (beta3)"))
dev.off()

# WHAT: estimate beta1_t, beta2_t, beta3_t MONTH BY MONTH by OLS.
# WHY : for a fixed lambda the model y_t = Lambda * f_t + e_t is linear in the
#       factors, so each month's factors are an ordinary least-squares fit of the
#       observed cross-section onto the loading matrix.
# HOW : the OLS "hat" matrix (L'L)^(-1)L' is constant, so we compute it once and
#       apply it to every month at once -> efficient and exact.
hat   <- solve(t(L) %*% L, t(L))                # 3 x N  (the OLS projector)
Y     <- coredata(yields)                       # T x N matrix of yields
betas <- hat %*% t(Y)                           # 3 x T matrix of estimated factors
factors <- xts(t(betas), order.by = index(yields))
colnames(factors) <- c("beta1", "beta2", "beta3")

# ---- Plot 4: the three estimated Nelson-Siegel factors -----------------------
png("output/q2_ns_factors.png", width = 900, height = 700)
plot.zoo(factors, main = "Estimated Nelson-Siegel factors",
         xlab = "Date", col = c("darkblue", "darkred", "darkgreen"))
dev.off()

# ---- Compare estimated factors with empirical level/slope/curvature ----------
# WHAT: overlay each estimated factor with its empirical counterpart and report
#       their correlation.
# WHY : it validates the Nelson-Siegel decomposition: a good model should give
#       factors that move closely with the model-free proxies.
# HOW TO READ (sign conventions):
#   beta1 ~  Level      (positive, near 1 correlation): beta1 IS the long-run level.
#   beta2 ~ -Slope      : by construction the slope loading is +1 at the short end,
#                         so beta2 is the SHORT-minus-LONG spread, i.e. the NEGATIVE
#                         of the empirical slope y(120)-y(3). We therefore compare
#                         beta2 with -Slope.
#   beta3 ~  Curvature  (same sign, different scale).
comparison <- merge(factors$beta1, empirical$Level,
                    -factors$beta2, empirical$Slope,
                    factors$beta3, empirical$Curvature)
colnames(comparison) <- c("beta1", "Level", "minus_beta2", "Slope", "beta3", "Curvature")

png("output/q2_factor_vs_empirical.png", width = 900, height = 800)
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
par(mfrow = c(1, 1)); dev.off()

cat("\n=========== FACTORS vs EMPIRICAL MEASURES (Question 2) ===========\n")
cat("corr(beta1 , Level)     :", round(cor(comparison$beta1,  comparison$Level), 3), "\n")
cat("corr(-beta2, Slope)     :", round(cor(comparison$minus_beta2, comparison$Slope), 3), "\n")
cat("corr(beta3 , Curvature) :", round(cor(comparison$beta3,  comparison$Curvature), 3), "\n")

# ---- Average fitting errors by maturity --------------------------------------
# WHAT: for each maturity, measure how well the 3-factor curve reproduces yields.
# WHY : the Nelson-Siegel model is parsimonious (3 factors for 8 yields), so we
#       must check it does not systematically miss some part of the curve.
# HOW : residuals = observed - fitted; we report the RMSE per maturity.
# HOW TO READ: small RMSE everywhere means 3 factors are enough; errors are
#       usually a little larger at the very short and very long ends, where the
#       curve is hardest to fit.
fitted_Y  <- L %*% betas                        # N x T fitted yields
resid_Y   <- t(Y) - fitted_Y                    # N x T residuals
rmse_mat  <- sqrt(rowMeans(resid_Y^2))          # one RMSE per maturity
names(rmse_mat) <- maturities
cat("\nAverage fitting error (RMSE, basis-point scale of %) by maturity:\n")
print(round(rmse_mat, 4))

png("output/q2_fitting_errors.png", width = 800, height = 600)
barplot(rmse_mat, col = "steelblue", xlab = "Maturity (months)", ylab = "RMSE (%)",
        main = "Average Nelson-Siegel fitting error by maturity")
dev.off()


# ==============================================================================
# QUESTION 3 - VAR MODEL FOR THE NELSON-SIEGEL FACTORS
# ==============================================================================

# WHAT: model the joint dynamics of f_t = (beta1, beta2, beta3) with a VAR.
# WHY : forecasting the whole yield curve reduces to forecasting only THREE
#       factors; the VAR captures their persistence and cross-dependence
#       (e.g. today's slope helps predict tomorrow's level).
F     <- coredata(factors)                      # T x 3 factor matrix
Tn    <- nrow(F)
horizons <- c(1, 6, 12)                         # forecast horizons required by the assignment

# WHAT: choose the VAR lag length by AIC (capped at 4 for parsimony).
# HOW TO READ: a low selected lag (often 1) confirms the factors are close to a
#       simple persistent VAR(1), the workhorse of the Diebold-Li model.
sel <- VARselect(F, lag.max = 4, type = "const")
p   <- as.integer(sel$selection["AIC(n)"])
cat("\n================ VAR specification (Question 3) ================\n")
cat("Selected VAR lag (AIC):", p, "\n")

# ---- Recursive (expanding-window) out-of-sample forecasts --------------------
# WHAT: starting from an initial training window, re-estimate the VAR every month
#       and forecast the factors h steps ahead; then map factor forecasts into
#       YIELD forecasts through the Nelson-Siegel loadings.
# WHY : recursive OOS forecasting mimics a real-time forecaster and is the fair
#       way to judge predictive accuracy (no look-ahead).
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

# INTERPRETATION (Question 3):
#   Compare the VAR tables with the random-walk tables cell by cell. A VAR entry
#   SMALLER than the corresponding RW entry means the factor model beats the naive
#   "no change" forecast for that maturity/horizon. In yield-curve data the random
#   walk is famously hard to beat at h = 1, while the VAR tends to gain at longer
#   horizons (h = 6, 12) because the factors mean-revert in a way the RW ignores.


# ==============================================================================
# QUESTION 4 - KALMAN-FILTER STATE-SPACE MODEL (one-step dynamic Nelson-Siegel)
# ==============================================================================

# TASK 1 - DIFFERENCE BETWEEN THE TWO APPROACHES (written answer in comments):
#   TWO-STEP OLS (Questions 2-3): first estimate the factors month-by-month by
#     OLS treating them as observed, THEN fit a separate VAR to those estimates.
#     Simple and fast, but it ignores (a) the estimation error in the factors and
#     (b) the feedback between measurement noise and factor dynamics; the two
#     steps are not jointly optimal.
#   ONE-STEP KALMAN (this question): factors are LATENT states. The measurement
#     equation (yields = loadings * factors + noise) and the transition equation
#     (factor VAR(1) dynamics) are estimated JOINTLY by maximum likelihood. The
#     Kalman filter optimally combines the cross-section of yields with the time
#     dynamics, propagates uncertainty, and yields smoother, more efficient
#     factor estimates and internally consistent forecasts.

# Model (with DIAGONAL Q and DIAGONAL H, as required):
#   Measurement: y_t      = Lambda * f_t + eps_t,        eps_t ~ N(0, H), H diagonal
#   Transition : f_t - mu = A (f_{t-1} - mu) + eta_t,    eta_t ~ N(0, Q), Q diagonal
# For transparency and stability we use a DIAGONAL A (each factor an AR(1)); this
# keeps the parameter count small and the optimisation robust, while still letting
# the three factors have different persistence. (A is the only restriction beyond
# the diagonal Q/H asked by the text; we flag it explicitly.)

# WHAT: a hand-coded Kalman filter that returns the log-likelihood AND, on request,
#       the filtered factors f_{t|t}.
# Parameter vector 'par' (length 17), all variances stored as LOGS to stay positive:
#   par[1:3]   = a   : diagonal AR(1) coefficients of A
#   par[4:6]   = mu  : unconditional means of the three factors
#   par[7:9]   = log(diag(Q))  (3 state-noise variances)
#   par[10:17] = log(diag(H))  (8 measurement-noise variances)
kalman_dns <- function(par, Y, L, return_states = FALSE) {
  m <- 3; N <- ncol(Y); Tn <- nrow(Y)
  a  <- par[1:3]
  mu <- par[4:6]
  Q  <- diag(exp(par[7:9]),  m)
  H  <- diag(exp(par[10:17]), N)
  A  <- diag(a, m)
  
  # Diffuse-free stationary initialisation: for a diagonal AR(1) the unconditional
  # state mean is mu and the unconditional variance solves P = A P A' + Q elementwise.
  f <- mu
  P <- diag(exp(par[7:9]) / pmax(1 - a^2, 1e-6), m)
  
  ll <- 0
  states <- matrix(NA_real_, Tn, m)
  for (t in 1:Tn) {
    # --- PREDICT: project the state one step ahead ---
    f_pred <- mu + A %*% (f - mu)
    P_pred <- A %*% P %*% t(A) + Q
    # --- INNOVATION: how surprising is today's yield cross-section ---
    v <- matrix(Y[t, ], N, 1) - L %*% f_pred
    S <- L %*% P_pred %*% t(L) + H
    Sinv <- solve(S)
    # accumulate the Gaussian log-likelihood of the innovation
    ll <- ll - 0.5 * (N * log(2 * pi) + log(det(S)) + t(v) %*% Sinv %*% v)
    # --- UPDATE: blend prediction with new information (Kalman gain K) ---
    K <- P_pred %*% t(L) %*% Sinv
    f <- f_pred + K %*% v
    P <- (diag(m) - K %*% L) %*% P_pred
    states[t, ] <- f
  }
  if (return_states) return(states)
  as.numeric(-ll)                                 # negative log-likelihood for minimisation
}

# WHAT: informed starting values taken from the two-step results of Q2/Q3.
# WHY : good starts make the 17-parameter likelihood optimisation fast and stable.
ar1 <- function(x) coef(lm(x[-1] ~ x[-length(x)]))[2]   # quick AR(1) coefficient
a0  <- pmin(pmax(apply(F, 2, ar1), -0.95), 0.99)        # persistence of each factor
mu0 <- colMeans(F)                                      # mean of each factor
q0  <- apply(F, 2, function(x) var(diff(x)))            # rough state-noise variance
h0  <- rmse_mat^2                                       # measurement variance = Q2 fitting variance
par0 <- c(a0, mu0, log(q0), log(h0))

# TASK 2 - estimate the model by maximum likelihood on the SAME training window
# used for the VAR (so the out-of-sample comparison is fair and look-ahead free).
lower <- c(rep(-0.999, 3), rep(-Inf, 3), rep(-20, 3), rep(-20, 8))
upper <- c(rep( 0.999, 3), rep( Inf, 3), rep( 10, 3), rep( 10, 8))
fit_kf <- optim(par0, kalman_dns, Y = Y[1:n_train, , drop = FALSE], L = L,
                method = "L-BFGS-B", lower = lower, upper = upper,
                control = list(maxit = 500))
par_hat <- fit_kf$par
cat("\n================ Kalman filter estimates (Question 4) ================\n")
cat("Factor persistences (diag A):", round(par_hat[1:3], 3), "\n")
cat("Factor means (mu)           :", round(par_hat[4:6], 3), "\n")

# TASK 3 - filtered factors over the FULL sample (states updated recursively with
# fixed estimated parameters) and their plot.
kf_states <- kalman_dns(par_hat, Y, L, return_states = TRUE)
kf_factors <- xts(kf_states, order.by = index(yields))
colnames(kf_factors) <- c("level", "slope", "curvature")

# HOW TO READ: these filtered factors should track the OLS factors of Q2 but look
#       SMOOTHER, because the Kalman filter optimally damps measurement noise.
png("output/q4_kalman_factors.png", width = 900, height = 700)
plot.zoo(kf_factors, main = "Kalman-filtered Nelson-Siegel factors",
         xlab = "Date", col = c("darkblue", "darkred", "darkgreen"))
dev.off()

# Side-by-side comparison of OLS vs Kalman factors (smoothness check for Q5.4).
png("output/q4_kalman_vs_ols.png", width = 900, height = 800)
par(mfrow = c(3, 1))
for (j in 1:3) {
  plot(index(yields), F[, j], type = "l", col = "grey50", xlab = "", ylab = "",
       main = paste0("Factor ", j, ": OLS (grey) vs Kalman (colour)"))
  lines(index(yields), kf_states[, j], col = c("darkblue", "darkred", "darkgreen")[j], lwd = 2)
}
par(mfrow = c(1, 1)); dev.off()

# TASK 4 - recursive out-of-sample forecasts for h = 1, 6, 12.
# WHAT: with fixed estimated parameters, at each origin t the filter has produced
#       f_{t|t} using data up to t only; the h-step forecast of a diagonal-A AR(1)
#       state is f_{t+h|t} = mu + A^h (f_{t|t} - mu), mapped to yields by Lambda.
a_hat  <- par_hat[1:3]; mu_hat <- par_hat[4:6]
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

png("output/q4_model_comparison.png", width = 800, height = 600)
barplot(comp_tbl, beside = TRUE, col = c("grey60", "darkred", "darkblue"),
        ylab = "Pooled RMSE (%)", xlab = "Forecast horizon",
        main = "Forecast accuracy: Random Walk vs VAR vs Kalman")
legend("topleft", bty = "n", fill = c("grey60", "darkred", "darkblue"),
       legend = c("Random Walk", "VAR", "Kalman"))
dev.off()


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