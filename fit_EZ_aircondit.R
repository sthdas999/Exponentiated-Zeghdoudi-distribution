# fit_EZ_aircondit_logparam.R

# Load required library
if (!requireNamespace("boot", quietly = TRUE)) install.packages("boot")
library(boot)

# Load the aircondit dataset
data("aircondit", package = "boot")
x <- aircondit$hours

# Negative log-likelihood using log-parametrization
negloglik_logpar <- function(logp, x) {
  c     <- exp(logp[1])
  theta <- exp(logp[2])
  u     <- 1 - (1 + theta * x + (theta^2 * x^2)/(theta + 2)) * exp(-theta * x)
  if (any(u <= 0) || !is.finite(u)) return(Inf)
  ll <- sum(
    log(c) + 3 * log(theta) + log(x) + log(1 + x) -
      theta * x - log(theta + 2) + (c - 1) * log(u)
  )
  -ll
}

# Run optimization using log-parameters
opt <- optim(
  par      = c(log(1), log(1)),  # initial logs
  fn       = negloglik_logpar,
  x        = x,
  method   = "L-BFGS-B",
  hessian  = TRUE
)

# Back-transform estimates
hat_c     <- exp(opt$par[1])
hat_theta <- exp(opt$par[2])

# Compute standard errors via the delta method
cov_logp    <- solve(opt$hessian)
se_log_c    <- sqrt(cov_logp[1, 1])
se_log_theta<- sqrt(cov_logp[2, 2])
se_c        <- hat_c * se_log_c
se_theta    <- hat_theta * se_log_theta

# Model fit statistics
n         <- length(x)
logLik    <- -opt$value
AIC       <- -2 * logLik + 2 * 2
BIC       <- -2 * logLik + log(n) * 2

# Display results
cat("===== EZ Fit to aircondit (log-parametrized) =====\n")
cat(sprintf("c = %.4f (SE = %.4f)\n", hat_c, se_c))
cat(sprintf("theta = %.4f (SE = %.4f)\n", hat_theta, se_theta))
cat(sprintf("Log-likelihood = %.4f\n", logLik))
cat(sprintf("AIC = %.2f    BIC = %.2f\n", AIC, BIC))
cat("===============================================\n")

# After running your MLE fit and computing `hat_c` and `hat_theta`

# Plot the histogram with density scaling
hist(x, prob = TRUE, breaks = 30,
     main = "EZ fit to aircondit",
     xlab = "Time between AC failures")

# Overlay the PDF of the fitted EZ distribution
curve({
  theta <- hat_theta; c <- hat_c
  pdf <- c * theta^3 * x * (1 + x) * exp(-theta * x) / (theta + 2) *
    (1 - (1 + theta * x + (theta^2 * x^2)/(theta + 2)) * exp(-theta * x))^(c - 1)
  pdf
},
from = min(x), to = max(x),
add = TRUE, col = "blue", lwd = 2,
n = 200)

# Overlay the CDF using a secondary y-axis
par(new = TRUE)
plot(ecdf(x), verticals = FALSE, do.points = FALSE,
     axes = FALSE, main = " ", xlab = "", ylab = "", col = "darkgreen", lwd = 2)

# Add the theoretical CDF of EZ
curve({
  theta <- hat_theta; c <- hat_c
  u <- 1 - (1 + theta * x + (theta^2 * x^2)/(theta + 2)) * exp(-theta * x)
  u^c
},
from = min(x), to = max(x),
add = TRUE, col = "darkgreen", lwd = 2,
n = 200)

# Add right-side y-axis for CDF
axis(side = 4, at = seq(0,1, by = 0.2))
mtext("Cumulative probability", side = 4, line = 3)
legend("topright",
       legend = c("EZ PDF", "Empirical & EZ CDF"),
       col = c("blue", "darkgreen"), lwd = 2)

