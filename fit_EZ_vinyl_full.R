# fit_EZ_vinyl_full.R

# 1?????? Install/load DataSetsUni for vinyl data
if (!requireNamespace("DataSetsUni", quietly = TRUE)) install.packages("DataSetsUni")
library(DataSetsUni)

# 2?????? Load the vinyl concentration dataset
x <- DataSetsUni::data_vinyl  # vinyl chloride (mg/L) data :contentReference[oaicite:8]{index=8}

# 3?????? Define negative log-likelihood with log-parameters
negloglik_logpar <- function(logp, x) {
  c     <- exp(logp[1])
  theta <- exp(logp[2])
  u     <- 1 - (1 + theta*x + (theta^2 * x^2)/(theta + 2)) * exp(-theta*x)
  if (any(u <= 0) || !is.finite(u)) return(Inf)
  ll <- sum(
    log(c) + 3*log(theta) + log(x) + log(1 + x) -
      theta*x - log(theta + 2) + (c - 1)*log(u)
  )
  -ll
}

# 4?????? Perform MLE optimization
opt <- optim(par = c(log(1), log(1)),
             fn = negloglik_logpar,
             x = x, method = "L-BFGS-B",
             hessian = TRUE)

# 5?????? Back-transform estimates and compute standard errors
hat_c     <- exp(opt$par[1])
hat_theta <- exp(opt$par[2])
cov_logp  <- solve(opt$hessian)
se_c      <- hat_c * sqrt(cov_logp[1, 1])
se_theta  <- hat_theta * sqrt(cov_logp[2, 2])

# 6?????? Compute model-fit statistics
n       <- length(x)
logLik  <- -opt$value
AIC     <- -2*logLik + 4
BIC     <- -2*logLik + log(n)*2

# 7?????? Print results
cat("===== EZ Fit (vinyl chloride data) =====\n")
cat(sprintf("c     = %.4f (SE = %.4f)\n", hat_c, se_c))
cat(sprintf("theta = %.4f (SE = %.4f)\n", hat_theta, se_theta))
cat(sprintf("logLik = %.4f | AIC = %.2f | BIC = %.2f\n", logLik, AIC, BIC))
cat("========================================\n")

# 8?????? Plot histogram, PDF, and CDF overlays (with centered ticks)
h <- hist(x, prob = TRUE, breaks = 30,
          main = "EZ Fit to Vinyl Chloride Data - PDF & CDF",
          xlab = "Vinyl Chloride (mg/L)",
          xaxt = "n")

axis(side = 1, at = h$mids, labels = round(h$mids, 2))

curve({
  theta <- hat_theta; c <- hat_c
  c * theta^3 * x * (1 + x) * exp(-theta*x) / (theta + 2) *
    (1 - (1 + theta*x + (theta^2 * x^2)/(theta + 2)) * exp(-theta*x))^(c - 1)
}, from = min(x), to = max(x),
add = TRUE, col = "blue", lwd = 2, n = 200)

par(new = TRUE)
plot(ecdf(x), verticals = FALSE, do.points = FALSE,
     axes = FALSE, main = " ", xlab = "", ylab = "",
     col = "darkgreen", lwd = 2)

curve({
  theta <- hat_theta; c <- hat_c
  u <- 1 - (1 + theta * x + (theta^2 * x^2)/(theta + 2)) * exp(-theta * x)
  u^c
}, from = min(x), to = max(x), add = TRUE,
col = "darkgreen", lwd = 2, n = 200)

axis(side = 4, at = seq(0, 1, by = 0.2))
mtext("Cumulative Probability", side = 4, line = 3)

legend("topright",
       legend = c("EZ PDF", "Empirical & EZ CDF"),
       col = c("blue", "darkgreen"), lwd = 2)
