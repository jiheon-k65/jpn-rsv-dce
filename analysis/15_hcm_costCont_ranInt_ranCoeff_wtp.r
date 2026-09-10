# ===========================================================================
# Title: 15_hcm_costCont_ranInt_ranCoeff_wtp.r
# Date: 2026-09-10
# ===========================================================================

# ------------------------------------------------------------------------
# SECTION 1: Load library
# ------------------------------------------------------------------------

# Load packages
library(here)
library(dplyr)
library(MASS)

# ------------------------------------------------------------------------
# SECTION 2: Load data & results
# ------------------------------------------------------------------------

data <- read.csv(here("data", "jpn_dataset_subset_imputed.csv"), header = TRUE)

results <- readRDS(here("output", "14.2_hcm_costCont_ranInt_ranCoeff_results.rds"))
est_list <- results$est
varcov_list <- results$varcov

# ------------------------------------------------------------------------
# SECTION 3: Define inputs for calculations
# ------------------------------------------------------------------------

# Keep only first imputation & unique respondents
df <- data %>% filter(imp == 1) %>% distinct(id, .keep_all = TRUE)
stopifnot(nrow(df) == length(unique(df$id)))

N <- nrow(df) # No. of respondents
M <- max(data$imp) # No. of imputations
S <- 500 # No. of simulations
F <- 100000 # Scaling factor for cost (JPY)

# ------------------------------------------------------------------------
# SECTION 4: Simulate parameter values
# ------------------------------------------------------------------------

# Initialize arrays to hold WTP values
wtp_preg_array <- array(NA, dim = c(N, S, M))
wtp_pret_array <- array(NA, dim = c(N, S, M))
wtp_hosp_array <- array(NA, dim = c(N, S, M))
    
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
        ln_mu_bCost <- theta[s, "ln_mu_bCost"]
        ln_sigma_bCost <- theta[s, "ln_sigma_bCost"]
        bCost <- - rlnorm(N, mean = ln_mu_bCost, sd = abs(ln_sigma_bCost)) # Lognormal

        # Simulate bPreg_mv for each respondent
        mu_bPreg_mv <- theta[s, "mu_bPreg_mv"]
        sigma_bPreg_mv <- theta[s, "sigma_bPreg_mv"]
        bPreg_mv <- rnorm(N, mean = mu_bPreg_mv, sd = abs(sigma_bPreg_mv))

        # Simulate bPreterm_mv for each respondent
        mu_bPreterm_mv <- theta[s, "mu_bPreterm_mv"]
        sigma_bPreterm_mv <- theta[s, "sigma_bPreterm_mv"]
        bPreterm_mv <- rnorm(N, mean = mu_bPreterm_mv, sd = abs(sigma_bPreterm_mv))

        # Calculate WTP for each respondent (Dim: N x S x M)
        wtp_preg <- - bPreg_mv / bCost * F
        wtp_pret <- - bPreterm_mv / bCost * F # Negative attribute
        wtp_hosp <- - bHosp0.3 / bCost * F
        stopifnot(
            length(wtp_preg) == N,
            length(wtp_pret) == N,
            length(wtp_hosp) == N
        )

        # Store WTPs in arrays (Dim: N x S x M)
        wtp_preg_array[, s, m] <- wtp_preg
        wtp_pret_array[, s, m] <- wtp_pret
        wtp_hosp_array[, s, m] <- wtp_hosp
    }
}

# Write output to text file
sink(here("output", "15_hcm_costCont_ranInt_ranCoeff_wtp_output.txt"))
cat("===========================================================================\n")
cat("Title: 15_hcm_costCont_ranInt_ranCoeff_wtp_output.txt\n")
cat("Date: 2026-09-10\n")
cat("===========================================================================\n")

# Calculate medians across respondents (Dim: S x M)
wtp_preg_medByN <- apply(wtp_preg_array, c(2, 3), median)
wtp_pret_medByN <- apply(wtp_pret_array, c(2, 3), median)
wtp_hosp_medByN <- apply(wtp_hosp_array, c(2, 3), median)

# Calculate IQR across respondents (Dim: S x M)
wtp_preg_q1ByN <- apply(wtp_preg_array, c(2, 3), quantile, probs = 0.25)
wtp_preg_q3ByN <- apply(wtp_preg_array, c(2, 3), quantile, probs = 0.75)
wtp_pret_q1ByN <- apply(wtp_pret_array, c(2, 3), quantile, probs = 0.25)
wtp_pret_q3ByN <- apply(wtp_pret_array, c(2, 3), quantile, probs = 0.75)
wtp_hosp_q1ByN <- apply(wtp_hosp_array, c(2, 3), quantile, probs = 0.25)
wtp_hosp_q3ByN <- apply(wtp_hosp_array, c(2, 3), quantile, probs = 0.75)

# Format results into a data frame
tot_summary <- data.frame(
    Attribute = c("Preg", "Preterm", "Hosp0.3"),
    Median = c(mean(wtp_preg_medByN), mean(wtp_pret_medByN), mean(wtp_hosp_medByN)),
    Mean = c(mean(wtp_preg_array), mean(wtp_pret_array), mean(wtp_hosp_array)),
    Q1 = c(mean(wtp_preg_q1ByN), mean(wtp_pret_q1ByN), mean(wtp_hosp_q1ByN)),
    Q3 = c(mean(wtp_preg_q3ByN), mean(wtp_pret_q3ByN), mean(wtp_hosp_q3ByN))
) %>% mutate(
    across(Median:Q3, ~ format(round(.x, 0), big.mark = ","))
)

cat("\nWillingness-to-pay (WTP) for attribute changes in JPY:\n")
print(tot_summary)
cat("\nNote: median & IQR are calculated across respondents, then averaged across simulations & imputations.\n")

# Switch off writing
sink()
cat("\nOutput successfully saved!\n")