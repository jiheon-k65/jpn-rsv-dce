# ===========================================================================
# Title: 20_hcm_costCont_ranCoeff_wtp.r
# Date: 2026-09-02
# ===========================================================================

# ------------------------------------------------------------------------
# SECTION 1: Load library
# ------------------------------------------------------------------------

# Load packages
library(here)
library(dplyr)
library(MASS)

# ------------------------------------------------------------------------
# SECTION 2: Load results from 19.2_hcm_costCont_ranCoeff.r
# ------------------------------------------------------------------------

# Load data from WD
data <- read.csv(here("data", "jpn_dataset_subset_imputed.csv"), header = TRUE)

# Load results from 19.2_hcm_costCont_ranCoeff.r
results <- readRDS(here("output", "19.2_hcm_costCont_ranCoeff_results.rds"))
est_list <- results$est
varcov_list <- results$varcov

# ------------------------------------------------------------------------
# SECTION 3: Define inputs for calculations
# ------------------------------------------------------------------------

M <- 25 # No. of imputations
N <- nrow(data) / M # No. of respondents
S <- 500 # No. of simulations
F <- 100000 # Scaling factor for cost (JPY)

# ------------------------------------------------------------------------
# SECTION 4: Simulate parameter values
# ------------------------------------------------------------------------

# Initialize array to hold simulated coefficients
bCost_array <- array(NA, dim = c(N, S, M))
bHosp_array <- array(NA, dim = c(N, S, M))
bPreg_array <- array(NA, dim = c(N, S, M))
bPret_array <- array(NA, dim = c(N, S, M))

for (m in 1:M) {
    cat("\n===========================================================================\n")
    cat("\nSimulating coefficients for imputation", m, "\n")

    # Extract estimates & variance-covariance matrix
    theta_hat <- est_list[[m]]
    sigma <- varcov_list[[m]]

    # Draw parameter values from multivariate normal distribution
    set.seed(123 + m)
    theta <- mvrnorm(n = S, mu = theta_hat, Sigma = sigma)
    colnames(theta) <- names(theta_hat)

    # Simulate random coefficients
    for (s in 1:S) {
        # Simulate fixed bHosp0.3
        bHosp0.3 <- rep(theta[s, "bHosp0.3"], N)

        # Simulate bCost for each respondent
        mu_bCost <- theta[s, "mu_bCost"]
        sigma_bCost <- theta[s, "sigma_bCost"]
        bCost <- rnorm(N, mean = mu_bCost, sd = abs(sigma_bCost))

        # Simulate bPreg_mv for each respondent
        mu_bPreg_mv <- theta[s, "mu_bPreg_mv"]
        sigma_bPreg_mv <- theta[s, "sigma_bPreg_mv"]
        bPreg_mv <- rnorm(N, mean = mu_bPreg_mv, sd = abs(sigma_bPreg_mv))

        # Simulate bPreterm_mv for each respondent
        mu_bPreterm_mv <- theta[s, "mu_bPreterm_mv"]
        sigma_bPreterm_mv <- theta[s, "sigma_bPreterm_mv"]
        bPreterm_mv <- rnorm(N, mean = mu_bPreterm_mv, sd = abs(sigma_bPreterm_mv))

        # Store simulated coefficients in arrays
        bHosp_array[, s, m] <- bHosp0.3
        bCost_array[, s, m] <- bCost
        bPreg_array[, s, m] <- bPreg_mv
        bPret_array[, s, m] <- bPreterm_mv
    }
}

# Write output to text file
sink(here("output", "20_hcm_costCont_ranCoeff_wtp_output.txt"))
cat("===========================================================================\n")
cat("Title: 20_hcm_costCont_ranCoeff_wtp_output.txt\n")
cat("Date: 2026-09-02\n")
cat("===========================================================================\n")

# Calculate mean coefficients across respondents (Dims: S x M)
bHosp_mean_byN <- apply(bHosp_array, c(2, 3), mean)
bCost_mean_byN <- apply(bCost_array, c(2, 3), mean)
bPreg_mean_byN <- apply(bPreg_array, c(2, 3), mean)
bPret_mean_byN <- apply(bPret_array, c(2, 3), mean)

# Calculate WTP across respondents (Dims: S x M)
wtp_hosp_byN <- (-bHosp_mean_byN) / bCost_mean_byN * F
wtp_preg_byN <- (-bPreg_mean_byN) / bCost_mean_byN * F
wtp_pret_byN <- (-bPret_mean_byN) / bCost_mean_byN * F # Negative attribute

# Calculate mean WTP across simulations & imputations
wtp_hosp_mean <- mean(wtp_hosp_byN)
wtp_preg_mean <- mean(wtp_preg_byN)
wtp_pret_mean <- mean(wtp_pret_byN)

# Calculate 95% uncertainty intervals (UIs)
wtp_preg_low_byN <- quantile(wtp_preg_byN, probs = 0.025)
wtp_preg_upp_byN <- quantile(wtp_preg_byN, probs = 0.975)
wtp_pret_low_byN <- quantile(wtp_pret_byN, probs = 0.025)
wtp_pret_upp_byN <- quantile(wtp_pret_byN, probs = 0.975)
wtp_hosp_low_byN <- quantile(wtp_hosp_byN, probs = 0.025)
wtp_hosp_upp_byN <- quantile(wtp_hosp_byN, probs = 0.975)

# Format results into a data frame
tot_summary <- data.frame(
    Attribute = c("Preg", "Preterm", "Hosp0.3"),
    Mean = c(wtp_preg_mean, wtp_pret_mean, wtp_hosp_mean),
    Lower_UI = c(wtp_preg_low_byN, wtp_pret_low_byN, wtp_hosp_low_byN),
    Upper_UI = c(wtp_preg_upp_byN, wtp_pret_upp_byN, wtp_hosp_upp_byN)
) %>% mutate(
    across(Mean:Upper_UI, ~ round(.x, 0))
)

cat("\nWillingness-to-pay (WTP) for attribute changes (JPY):\n")
print(tot_summary)
cat("\nNote: WTP for Preterm is -ve as it is a negative attribute\n")

# Switch off writing
sink()
