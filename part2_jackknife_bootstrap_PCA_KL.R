# ====================================================
# PART 2: Jackknife and Bootstrap Methods for PCA
# ====================================================

# Load data
data <- read.table("BodyMeasurements.txt", header = TRUE)

# ====================================================
# A: Jackknife estimation for theta (proportion of variance explained by PC1)
# ====================================================

# Function to calculate theta = λ1 / sum(λi)
calculate_theta <- function(data) {
  Sigma <- cov(data)
  eigenvalues <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
  lambda1 <- eigenvalues[1]
  sum_lambdas <- sum(eigenvalues)
  theta <- lambda1 / sum_lambdas
  return(theta)
}

# Jackknife function
jknife <- function(data, fxn, ci = 0.95) {
  theta <- fxn(data)
  n <- nrow(data)
  partials <- rep(0, n)
  
  for (i in 1:n) {
    partials[i] <- fxn(data[-i, ])
  }
  
  pseudos <- (n * theta) - (n - 1) * partials
  jack.est <- mean(pseudos)
  jack.se <- sd(pseudos) / sqrt(n)
  alpha <- 1 - ci
  CI <- qt(alpha/2, n-1, lower.tail = FALSE) * jack.se
  jack.ci <- c(jack.est - CI, jack.est + CI)
  
  list(est = jack.est, se = jack.se, ci = jack.ci, pseudos = pseudos, partials = partials)
}

data <- read.table("BodyMeasurements.txt", header = TRUE)
results <- jknife(data, calculate_theta)

# Results
cat("Jackknife estimate θ:", results$est, "\n")
cat("Standard error:", results$se, "\n")
cat("95% Confidence Interval:", results$ci, "\n")

# Bias calculation
jack.bias <- (nrow(data) - 1) * (mean(results$partials) - calculate_theta(data))
cat("Bias:", jack.bias, "\n")

# ====================================================
# B: Bootstrap confidence intervals for theta
# ====================================================

# Load data
data <- read.table("BodyMeasurements.txt", header = TRUE)
set.seed(2025)

# Function to calculate theta
calculate_theta <- function(data) {
  Sigma <- cov(data)
  eigenvalues <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
  lambda1 <- eigenvalues[1]
  sum_lambdas <- sum(eigenvalues)
  theta <- lambda1 / sum_lambdas
  return(theta)
}

n <- nrow(data)
B <- 1000
alpha <- 0.05

# Calculate theta from original data
theta_hat <- calculate_theta(data)
cat("Theta estimate:", theta_hat, "\n")

index.star <- matrix(sample(1:n, B*n, replace = TRUE), nrow = B)

# Bootstrap estimates for theta
theta.star <- apply(index.star, 1, function(m) calculate_theta(data[m,]))

# Method 1: Basic bootstrap interval
theta.quantiles <- quantile(theta.star, c(alpha/2, 1-alpha/2), type = 1)
CI1_theta <- unname(2*theta_hat - rev(theta.quantiles))

# Method 2: Percentile bootstrap interval
CI2_theta <- unname(theta.quantiles)

# Method 3: BCa bootstrap interval
theta.jack <- numeric(n)
for (j in 1:n) theta.jack[j] <- calculate_theta(data[-j,])
a_theta <- sum((theta.jack - mean(theta.jack))^3) / 
  (6 * sum((theta.jack - mean(theta.jack))^2)^(3/2))
z0_theta <- qnorm(mean(theta.star < theta_hat))
z <- qnorm(0.975)
alpha1_theta <- pnorm(z0_theta + (z0_theta - z)/(1 - a_theta*(z0_theta - z)))
alpha2_theta <- pnorm(z0_theta + (z0_theta + z)/(1 - a_theta*(z0_theta + z)))
CI3_theta <- quantile(theta.star, c(alpha1_theta, alpha2_theta), type = 1, names = FALSE)

# Method 4: Bootstrap-t for theta
T.star_theta <- numeric(B)
for (b in 1:B) {
  index.star2 <- matrix(sample(1:n, 100*n, replace = TRUE), nrow = 100)
  theta.star2 <- apply(index.star2, 1, function(m) calculate_theta(data[index.star[b,],][m,]))
  T.star_theta[b] <- (theta.star[b] - theta_hat) / sd(theta.star2)
}
CI4_theta <- theta_hat - quantile(T.star_theta, c(1-alpha/2, alpha/2), type = 1, names = FALSE) * sd(theta.star)

# Print results
cat("\n--- 95% Confidence Intervals for theta ---\n")
cat("Basic bootstrap:", CI1_theta, "\n")
cat("Percentile bootstrap:", CI2_theta, "\n")
cat("BCa bootstrap:", CI3_theta, "\n")
cat("Bootstrap-t:", CI4_theta, "\n")

# ====================================================
# C: Hypothesis testing for μ1 = (7/17)μ5 - CORRECTED VERSION
# ====================================================

# Load data
data <- read.table("BodyMeasurements.txt", header = TRUE)

# Data
X1 <- data[,1] # 1st column
X5 <- data[,5] # 5th column
n <- length(X1)
m <- length(X5)
k <- 7/17

X5_scaled <- k * X5

xy.bar <- mean(c(X1, X5_scaled))

# Data transformation under H0
x.tilde <- X1 - mean(X1) + xy.bar
y.tilde <- X5_scaled - mean(X5_scaled) + xy.bar

# Verify transformation
cat("\n=== HYPOTHESIS TESTING: μ1 = (7/17)μ5 ===\n")
cat("Transformation verification:\n")
cat("mean(x.tilde):", mean(x.tilde), "\n")
cat("mean(y.tilde):", mean(y.tilde), "\n")
cat("Equality:", abs(mean(x.tilde) - mean(y.tilde)) < 1e-10, "\n\n")

# OBSERVED STATISTIC - CORRECTED
t_obs <- (mean(X1) - mean(X5_scaled)) / 
  sqrt(var(X1)/n + var(X5_scaled)/m)
# Equivalent: (mean(X1) - k*mean(X5)) / sqrt(var(X1)/n + k^2*var(X5)/m)

cat("Observed t-statistic:", t_obs, "\n")

# Bootstrap simulation under H0
set.seed(2025)
B <- 1000
t.star <- numeric(B)

for(b in 1:B) {
  index.x.star <- sample(1:n, n, replace = TRUE)
  x.star <- x.tilde[index.x.star]
  index.y.star <- sample(1:m, m, replace = TRUE)
  y.star <- y.tilde[index.y.star]
  
  # Bootstrap t-statistic
  t.star[b] <- (mean(x.star) - mean(y.star)) / 
    sqrt(var(x.star)/n + var(y.star)/m)
}

# Calculate p-value
p_value <- 2 * min(mean(t.star <= t_obs), mean(t.star >= t_obs))

# Results
cat("\n=== RESULTS ===\n")
cat("Bootstrap p-value:", p_value, "\n")
cat("Significance level (α): 0.05\n")
cat("Conclusion:", ifelse(p_value < 0.05, "Reject H0", "Fail to reject H0"), "\n\n")

# ====================================================
# D: Bootstrap for P(Z > 55) with 98% CI
# ====================================================

# Load data
data <- read.table("BodyMeasurements.txt", header = TRUE)

# Data - 2nd column
Z <- data[,2]
n <- length(Z)
B <- 1000
alpha <- 0.02  # 98% confidence level

# Function to estimate θ = P(Z > 55)
estimate_theta <- function(z) {
  mean(z > 55)
}

# Simple empirical estimator
theta_hat <- estimate_theta(Z)
cat("\n=== PROBABILITY P(Z > 55) ===\n")
cat("Simple empirical estimator θ = P(Z > 55):", theta_hat, "\n")
cat("Number of observations Z > 55:", sum(Z > 55), "out of", n, "\n")

set.seed(2025)
index.star <- matrix(sample(1:n, B*n, replace = TRUE), nrow = B)

# Bootstrap estimates for theta
theta.star <- apply(index.star, 1, function(m) estimate_theta(Z[m]))

# Method 1: Basic bootstrap interval
theta.quantiles <- quantile(theta.star, c(alpha/2, 1-alpha/2), type = 1)
CI1_theta <- unname(2*theta_hat - rev(theta.quantiles))

# Method 2: Percentile bootstrap interval
CI2_theta <- unname(theta.quantiles)

# Method 3: BCa bootstrap interval
theta.jack <- numeric(n)
for (j in 1:n) theta.jack[j] <- estimate_theta(Z[-j])

a_theta <- sum((theta.jack - mean(theta.jack))^3) / 
  (6 * sum((theta.jack - mean(theta.jack))^2)^(3/2))
z0_theta <- qnorm(mean(theta.star < theta_hat))
z <- qnorm(1 - alpha/2)  # for 98% CI
alpha1_theta <- pnorm(z0_theta + (z0_theta - z)/(1 - a_theta*(z0_theta - z)))
alpha2_theta <- pnorm(z0_theta + (z0_theta + z)/(1 - a_theta*(z0_theta + z)))
CI3_theta <- quantile(theta.star, c(alpha1_theta, alpha2_theta), type = 1, names = FALSE)

# Method 4: Bootstrap-t for theta
T.star_theta <- numeric(B)
for (b in 1:B) {
  index.star2 <- matrix(sample(1:n, 100*n, replace = TRUE), nrow = 100)
  theta.star2 <- apply(index.star2, 1, function(m) estimate_theta(Z[index.star[b,]][m]))
  T.star_theta[b] <- (theta.star[b] - theta_hat) / sd(theta.star2)
}
CI4_theta <- theta_hat - quantile(T.star_theta, c(1-alpha/2, alpha/2), type = 1, names = FALSE) * sd(theta.star)

# Print results
cat("\n--- 98% Confidence Intervals for θ = P(Z > 55) ---\n")
cat("Basic bootstrap:", CI1_theta, "\n")
cat("Percentile bootstrap:", CI2_theta, "\n")
cat("BCa bootstrap:", CI3_theta, "\n")
cat("Bootstrap-t:", CI4_theta, "\n")

# ====================================================
# E: Hypothesis testing with KL divergence bootstrap
# ====================================================

# Load data
data <- read.table("BodyMeasurements.txt", header = TRUE)
set.seed(2025)

# Data - 1st column
X <- data[,1]
n <- length(X)
mu0 <- 52 # H₀: μ₁ = 52

# Step 1: Classic t-test
t_test_classic <- t.test(X, mu = mu0)
p_value_classic <- t_test_classic$p.value
cat("=== CLASSIC t-test ===\n")
cat("t-statistic:", t_test_classic$statistic, "\n")
cat("p-value:", p_value_classic, "\n\n")

# Step 2: Standard Bootstrap - CORRECTED
t_obs <- (mean(X) - mu0) / (sd(X)/sqrt(n))

# Correct transformation under H₀
x.tilde <- X - mean(X) + mu0
set.seed(2025)
B <- 1000
t.star_normal <- numeric(B)
for(b in 1:B) {
  index.star <- sample(1:n, n, replace = TRUE)
  x.star <- x.tilde[index.star]
  t.star_normal[b] <- (mean(x.star) - mu0) / (sd(x.star)/sqrt(n))
}
p_value_bootstrap_normal <- 2 * min(mean(t.star_normal <= t_obs), mean(t.star_normal >= t_obs))
cat("=== STANDARD BOOTSTRAP ===\n")
cat("p-value:", p_value_bootstrap_normal, "\n\n")

# Step 3: Bootstrap with KL divergence
# Function to calculate p_i
p <- function(lambda, x) {
  exp(lambda * x) / sum(exp(lambda * x))
}

# Find optimal lambda
lambda1 <- -0.5
while(sum(p(lambda1, X) * X) > mu0) lambda1 <- 1.5 * lambda1

lambda2 <- 0.5
while(sum(p(lambda2, X) * X) < mu0) lambda2 <- 1.5 * lambda2

lambda0 <- (lambda1 + lambda2) / 2

while(lambda2 - lambda1 + abs(sum(p(lambda0, X) * X) - mu0) > 1e-8) {
  if(sum(p(lambda0, X) * X) < mu0) {
    lambda1 <- lambda0
  } else {
    lambda2 <- lambda0
  }
  lambda0 <- (lambda1 + lambda2) / 2
}

# Optimal probabilities
p0 <- p(lambda0, X)

cat("=== BOOTSTRAP WITH KL ===\n")
cat("Optimal λ:", lambda0, "\n")
cat("Verification - ∑p_iX_i:", sum(p0 * X), "(should be ≈", mu0, ")\n")
cat("∑p_i:", sum(p0), "\n\n")

# Bootstrap with optimal probabilities
index.star_kl <- matrix(sample(1:n, B*n, replace = TRUE, prob = p0), nrow = B)
t.star_kl <- apply(index.star_kl, 1, function(m) (mean(X[m]) - mu0) / (sd(X[m])/sqrt(n)))

p_value_bootstrap_kl <- 2 * min(mean(t.star_kl <= t_obs), mean(t.star_kl >= t_obs))
cat("p-value (Bootstrap KL):", p_value_bootstrap_kl, "\n\n")

# Step 4: Comparison of all methods
cat("=== COMPARISON OF ALL METHODS ===\n")
cat("Classic t-test p-value:", p_value_classic, "\n")
cat("Standard Bootstrap p-value:", p_value_bootstrap_normal, "\n")
cat("KL Bootstrap p-value:", p_value_bootstrap_kl, "\n")