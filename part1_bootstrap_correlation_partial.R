# Load data
data <- read.table("Games.txt", header = TRUE)
head(data)
set.seed(2025)

# Function for simple correlation coefficient X,Y
my.cor_xy <- function(data) {
  cor(data)[1,2]
}

# Function for partial correlation coefficient X,Y | Z
my.partial_cor <- function(data) {
  Sigma <- cov(data) # Covariance matrix
  Omega <- solve(Sigma) # Inverse matrix (precision matrix)
  
  theta <- -Omega[1,2] / sqrt(Omega[1,1] * Omega[2,2])
  return(theta)
}

n <- nrow(data)
B <- 1000
alpha <- 0.05

# Calculate estimators from original data
cor_xy <- my.cor_xy(data)
partial_cor <- my.partial_cor(data)
cat("Correlation coefficient X,Y:", cor_xy, "\n")
cat("Partial correlation coefficient X,Y|Z:", partial_cor, "\n")

index.star <- matrix(sample(1:n, B*n, replace = TRUE), nrow = B)

# Simple correlation coefficient
cor_xy.star <- apply(index.star, 1, function(m) my.cor_xy(data[m,]))

# Method 1: Basic bootstrap interval
cor_xy.quantiles <- quantile(cor_xy.star, c(alpha/2, 1-alpha/2), type = 1)
CI1_cor <- unname(2*cor_xy - rev(cor_xy.quantiles))

# Method 2: Percentile bootstrap interval
CI2_cor <- unname(cor_xy.quantiles)

# Method 3: BCa bootstrap interval (for correlation coefficient)
cor_xy.jack <- numeric(n)
for (j in 1:n) cor_xy.jack[j] <- my.cor_xy(data[-j,])
a_cor <- sum((cor_xy.jack - mean(cor_xy.jack))^3) / (6 * sum((cor_xy.jack - mean(cor_xy.jack))^2)^(3/2))
z0_cor <- qnorm(mean(cor_xy.star < cor_xy))
z <- qnorm(0.975)
alpha1_cor <- pnorm(z0_cor + (z0_cor - z)/(1 - a_cor*(z0_cor - z)))
alpha2_cor <- pnorm(z0_cor + (z0_cor + z)/(1 - a_cor*(z0_cor + z)))
CI3_cor <- quantile(cor_xy.star, c(alpha1_cor, alpha2_cor), type = 1, names = FALSE)

# Method 4: Bootstrap-t for simple correlation
T.star_cor <- numeric(B)
for (b in 1:B) {
  index.star2 <- matrix(sample(1:n, 100*n, replace = TRUE), nrow = 100)
  cor_xy.star2 <- apply(index.star2, 1, function(m) my.cor_xy(data[index.star[b,],][m,]))
  T.star_cor[b] <- (cor_xy.star[b] - cor_xy) / sd(cor_xy.star2)
}
CI4_cor <- cor_xy - quantile(T.star_cor, c(1-alpha/2, alpha/2), type = 1, names = FALSE) * sd(cor_xy.star)

partial_cor.star <- apply(index.star, 1, function(m) my.partial_cor(data[m,]))

# Method 1: Basic bootstrap interval
partial_cor.quantiles <- quantile(partial_cor.star, c(alpha/2, 1-alpha/2), type = 1)
CI1_partial <- unname(2*partial_cor - rev(partial_cor.quantiles))

# Method 2: Percentile bootstrap interval
CI2_partial <- unname(partial_cor.quantiles)

# Method 3: BCa bootstrap interval
partial_cor.jack <- numeric(n)
for (j in 1:n) partial_cor.jack[j] <- my.partial_cor(data[-j,])
a_partial <- sum((partial_cor.jack - mean(partial_cor.jack))^3) / (6 * sum((partial_cor.jack - mean(partial_cor.jack))^2)^(3/2))
z0_partial <- qnorm(mean(partial_cor.star < partial_cor))
alpha1_partial <- pnorm(z0_partial + (z0_partial - z)/(1 - a_partial*(z0_partial - z)))
alpha2_partial <- pnorm(z0_partial + (z0_partial + z)/(1 - a_partial*(z0_partial + z)))
CI3_partial <- quantile(partial_cor.star, c(alpha1_partial, alpha2_partial), type = 1, names = FALSE)

# Method 4: Bootstrap-t
T.star_partial <- numeric(B)
for (b in 1:B) {
  index.star2 <- matrix(sample(1:n, 100*n, replace = TRUE), nrow = 100)
  partial_cor.star2 <- apply(index.star2, 1, function(m) my.partial_cor(data[index.star[b,],][m,]))
  T.star_partial[b] <- (partial_cor.star[b] - partial_cor) / sd(partial_cor.star2)
}
CI4_partial <- partial_cor - quantile(T.star_partial, c(1-alpha/2, alpha/2), type = 1, names = FALSE) * sd(partial_cor.star)

# Print results
cat("\n=== CORRELATION COEFFICIENT X,Y ===\n")
cat("Basic Bootstrap CI:", CI1_cor, "\n")
cat("Percentile Bootstrap CI:", CI2_cor, "\n")
cat("BCa Bootstrap CI:", CI3_cor, "\n")
cat("Bootstrap-t CI:", CI4_cor, "\n")

cat("\n=== PARTIAL CORRELATION COEFFICIENT X,Y|Z ===\n")
cat("Basic Bootstrap CI:", CI1_partial, "\n")
cat("Percentile Bootstrap CI:", CI2_partial, "\n")
cat("BCa Bootstrap CI:", CI3_partial, "\n")
cat("Bootstrap-t CI:", CI4_partial, "\n")