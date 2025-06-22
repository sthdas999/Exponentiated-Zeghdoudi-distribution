# fit_EZ_carfibres_full.R

# 1. Install and load package
if (!requireNamespace("DataSetsUni", quietly = TRUE)) install.packages("DataSetsUni")
library(DataSetsUni)

# 2. Load breaking stress dataset
x <- DataSetsUni::data_carfibres  # 100 carbon fibre breaking stresses :contentReference[oaicite:3]{index=3}

# 3. Define negative log-likelihood via log-parameters
negloglik_logpar <- function(logp, x) {
  c     <- exp(logp[1])
  theta <- exp(logp[2])
  u <- 1 - (1 + theta*x + (theta^2 * x^2)/(theta + 2)) * exp(-theta * x)
  if (any(u <= 0) || !is.finite(u)) return(Inf)
  ll <- sum(
    log(c) + 3*log(theta) + log(x) + log(1 + x) -
      theta * x - log(theta + 2) + (c - 1)*log(u)
  )
  -ll
}

# 4. Run the optimizer
opt <- optim(
  par     = c(log(1), log(1)),
  fn      = negloglik_logpar,
  x       = x,
  method  = "L-BFGS-B",
  hessian = TRUE
)

# 5. Transform back and compute standard errors
hat_c     <- exp(opt$par[1])
hat_theta <- exp(opt$par[2])
cov_logp  <- solve(opt$hessian)
se_c      <- hat_c * sqrt(cov_logp[1, 1])
se_theta  <- hat_theta * sqrt(cov_logp[2, 2])

# 6. Fit statistics
n       <- length(x)
logLik  <- -opt$value
AIC     <- -2*logLik + 4
BIC     <- -2*logLik + log(n)*2

# 7. Print parameter estimates
cat("===== EZ Fit (carbon fibres) =====\n")
cat(sprintf("c     = %.4f (SE = %.4f)\n", hat_c, se_c))
cat(sprintf("theta = %.4f (SE = %.4f)\n", hat_theta, se_theta))
cat(sprintf("logLik = %.4f | AIC = %.2f | BIC = %.2f\n", logLik, AIC, BIC))
cat("=================================\n")

# 8A. Plot PDF with histogram
hist(x, prob = TRUE, breaks = 30,
     main = " ",
     xlab = "Stress (GPa)",
     xaxt = "n", col = "lightgray", border = "white")

# Center tick labels under each bar
h <- hist(x, plot = FALSE, breaks = 30)
axis(side = 1, at = h$mids, labels = round(h$mids, 2))

# Add EZ PDF curve
curve({
  theta <- hat_theta; c <- hat_c
  c * theta^3 * x * (1 + x) * exp(-theta * x) / (theta + 2) *
    (1 - (1 + theta * x + (theta^2 * x^2)/(theta + 2)) * exp(-theta * x))^(c - 1)
}, from = min(x), to = max(x), add = TRUE,
col = "blue", lwd = 2, n = 200)

legend("topright", legend = "EZ PDF", col = "blue", lwd = 2)

# 8B. Plot CDF (empirical and theoretical)
plot(ecdf(x), verticals = TRUE, do.points = FALSE,
     main = " ",
     xlab = "Stress (GPa)",
     ylab = "Cumulative Probability",
     col = "darkgreen", lwd = 2)

curve({
  theta <- hat_theta; c <- hat_c
  u <- 1 - (1 + theta * x + (theta^2 * x^2)/(theta + 2)) * exp(-theta * x)
  u^c
}, from = min(x), to = max(x), add = TRUE,
col = "green3", lwd = 2, n = 200)

legend("bottomright",
       legend = c("Empirical CDF", "EZ CDF"),
       col = c("darkgreen", "green3"),
       lwd = 2)

