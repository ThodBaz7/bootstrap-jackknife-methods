# Statistical Analysis with Bootstrap and Jackknife Methods

[![R](https://img.shields.io/badge/R-4.0%2B-blue)](https://www.r-project.org/)
[![Bootstrap](https://img.shields.io/badge/Bootstrap-Resampling-green)](https://en.wikipedia.org/wiki/Bootstrapping_(statistics))
[![Jackknife](https://img.shields.io/badge/Jackknife-Estimation-orange)](https://en.wikipedia.org/wiki/Jackknife_resampling)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## 📋 Description
Analysis of statistical datasets using **bootstrap and jackknife resampling methods**. The project includes correlation analysis, partial correlation, PCA variance explanation, and hypothesis testing with various bootstrap techniques (basic, percentile, BCa, bootstrap-t) and KL divergence bootstrap.

This project was implemented as part of the **MSc in Applied Statistics** at the **University of Piraeus**.

## 📊 Datasets
- **Games.txt:** Dataset for correlation and partial correlation analysis (X, Y, Z variables)
- **BodyMeasurements.txt:** Body measurements dataset with 5 variables and 100 observations for PCA and hypothesis testing

## 🗂️ Project Structure
```
statistiki-ergasia/
│
├── part1_bootstrap_correlation_partial.R    # Exercise 1: Bootstrap for correlation coefficients
├── part2_jackknife_bootstrap_PCA_KL.R       # Exercise 2: Jackknife, Bootstrap PCA, KL divergence
├── Εργασία 1.pdf                             # Full assignment report (in Greek)
│
├── data/                                      # Dataset files
│   ├── Games.txt
│   └── BodyMeasurements.txt
│
└── README.md                                   # Project documentation
```

## 🚀 Installation & Execution

### Prerequisites
- **R** (version 3.6 or higher)
- No external R packages required (uses **base R functions only**)

### Execution Steps
```bash
# 1. Clone the repository
git clone https://github.com/ThodBaz7/bootstrap-jackknife-methods.git
cd bootstrap-jackknife-methods

# 2. Run the analysis scripts in R
# Open R or RStudio and run:
source("part1_bootstrap_correlation_partial.R")
source("part2_jackknife_bootstrap_PCA_KL.R")
```

## 📈 Analysis Stages

### 1️⃣ **Exercise 1: Correlation and Partial Correlation**
- Simple correlation coefficient X,Y calculation
- Partial correlation coefficient X,Y|Z (controlling for Z)
- Bootstrap methods for 95% confidence intervals:
  - Basic bootstrap
  - Percentile bootstrap
  - BCa bootstrap (bias-corrected and accelerated)
  - Bootstrap-t
- Comparison of all four methods for both correlation types

### 2️⃣ **Exercise 2: Jackknife and Bootstrap for PCA**

#### **Part A: Jackknife Estimation**
- Calculation of **θ = λ₁ / Σλᵢ** (proportion of variance explained by first principal component)
- Jackknife bias estimation
- Jackknife standard error calculation
- Comparison with original estimate

#### **Part B: Bootstrap Confidence Intervals for θ**
- 95% confidence intervals using all four bootstrap methods
- Method comparison and interpretation

#### **Part C: Hypothesis Testing for μ₁ = (7/17)μ₅**
- Testing equality of means between first and fifth columns
- Bootstrap simulation under H₀
- P-value calculation and interpretation

#### **Part D: Bootstrap for P(Z > 55)**
- 98% confidence intervals for probability estimation
- All four bootstrap methods applied to binary outcome
- Probability estimation for body measurements > 55

#### **Part E: Hypothesis Testing with KL Divergence Bootstrap**
- Classic t-test for H₀: μ₁ = 52
- Standard bootstrap test
- KL divergence bootstrap:
  - Optimal probability weights under H₀
  - Kullback-Leibler divergence minimization
  - Bootstrap sampling with optimal probabilities
- Comparison of all three methods

## 📦 R Code Features
All scripts use **base R only** - no external packages required.

### Key Functions:
- `my.cor_xy()` - Simple correlation calculation
- `my.partial_cor()` - Partial correlation calculation
- `calculate_theta()` - PCA variance proportion
- `jknife()` - Jackknife estimation function
- `p()` - KL divergence probability weights

## 📝 Results Summary

### **Exercise 1 Results:**
- **Simple correlation X,Y:** Nearly perfect positive correlation (all CIs near 1)
- **Partial correlation X,Y|Z:** All CIs contain 0 → spurious correlation explained by Z

### **Exercise 2 Results:**
- **θ estimate:** ~0.875 (first PC explains 87.5% of variance)
- **Jackknife bias:** Negligible (-0.0001)
- **Jackknife standard error:** 0.017
- **KL Bootstrap p-value:** ~0.0488 (consistent with classic t-test)

## 📄 Report
The full assignment report (in Greek) is available at [`Εργασία 1.pdf`](Εργασία%201.pdf)

## 📜 License
MIT License - Free use and modification with attribution.

## 📬 Contact
thbazhba@gmail.com
