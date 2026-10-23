# ===========================================================================
# Title: 15_calc_wtp.r
# Date: 2026-10-23
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
robvarcov_list <- results$robvarcov # Robust variance-covariance matrices

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
wtp_red_array <- array(NA, dim = c(N, S, M))
wtp_pret_array <- array(NA, dim = c(N, S, M))
wtp_hosp_array <- array(NA, dim = c(N, S, M))
    
for (m in 1:M) {
    cat("\n===========================================================================\n")
    cat("\nCalculating WTPs for imputation", m, "\n")

    # Extract estimates & variance-covariance matrix
    theta_hat <- est_list[[m]]
    sigma <- robvarcov_list[[m]]

    # Draw parameter values from multivariate normal distribution
    set.seed(123 + m)
    theta <- mvrnorm(n = S, mu = theta_hat, Sigma = sigma)
    colnames(theta) <- names(theta_hat)

    # Simulate random coefficients
    for (s in 1:S) {
        # Simulate fixed coefficients
        bRedness <- rep(theta[s, "bRedness"], N)

        # Simulate bCost for each respondent
        ln_mu_cost <- theta[s, "ln_mu_cost"]
        ln_sigma_cost <- theta[s, "ln_sigma_cost"]
        bCost <- - rlnorm(N, mean = ln_mu_cost, sd = abs(ln_sigma_cost))

        # Simulate bPreg_mv for each respondent
        mu_preg_mv <- theta[s, "mu_preg_mv"]
        sigma_preg_mv <- theta[s, "sigma_preg_mv"]
        bPreg_mv <- rnorm(N, mean = mu_preg_mv, sd = abs(sigma_preg_mv))

        # Simulate bPreterm_mv for each respondent
        mu_preterm_mv <- theta[s, "mu_preterm_mv"]
        sigma_preterm_mv <- theta[s, "sigma_preterm_mv"]
        bPreterm_mv <- rnorm(N, mean = mu_preterm_mv, sd = abs(sigma_preterm_mv))

        # Simulate bHosp0.3 for each respondent
        mu_hosp0.3 <- theta[s, "mu_hosp0.3"]
        sigma_hosp0.3 <- theta[s, "sigma_hosp0.3"]
        bHosp0.3 <- rnorm(N, mean = mu_hosp0.3, sd = abs(sigma_hosp0.3))

        # Calculate WTP for each respondent (Dim: N x S x M)
        wtp_preg <- - bPreg_mv / bCost * F
        wtp_red <- - bRedness / bCost * F # Negative attribute
        wtp_pret <- - bPreterm_mv / bCost * F # Negative attribute
        wtp_hosp <- - bHosp0.3 / bCost * F
        stopifnot(
            length(wtp_preg) == N,
            length(wtp_red) == N,
            length(wtp_pret) == N,
            length(wtp_hosp) == N
        )

        # Store WTPs in arrays (Dim: N x S x M)
        wtp_preg_array[, s, m] <- wtp_preg
        wtp_red_array[, s, m] <- wtp_red
        wtp_pret_array[, s, m] <- wtp_pret
        wtp_hosp_array[, s, m] <- wtp_hosp
    }
}

# Write output to text file
sink(here("output", "15_calc_wtp_output.txt"))
cat("===========================================================================\n")
cat("Title: 15_calc_wtp_output.txt\n")
cat("Date: 2026-10-23\n")
cat("===========================================================================\n")

# Calculate medians across respondents (Dim: S x M)
wtp_preg_medByN <- apply(wtp_preg_array, c(2, 3), median)
wtp_pret_medByN <- apply(wtp_pret_array, c(2, 3), median)
wtp_red_medByN <- apply(wtp_red_array, c(2, 3), median)
wtp_hosp_medByN <- apply(wtp_hosp_array, c(2, 3), median)

# Calculate 95% UI across medians
wtp_preg_lowByN <- quantile(wtp_preg_medByN, 0.025)
wtp_preg_uppByN <- quantile(wtp_preg_medByN, 0.975)
wtp_red_lowByN <- quantile(wtp_red_medByN, 0.025)
wtp_red_uppByN <- quantile(wtp_red_medByN, 0.975)
wtp_pret_lowByN <- quantile(wtp_pret_medByN, 0.025)
wtp_pret_uppByN <- quantile(wtp_pret_medByN, 0.975)
wtp_hosp_lowByN <- quantile(wtp_hosp_medByN, 0.025)
wtp_hosp_uppByN <- quantile(wtp_hosp_medByN, 0.975)

# Format results into a data frame
tot_summary <- data.frame(
    Attribute = c("Preg", "Redness", "Preterm", "Hosp0.3"),
    Median = c(mean(wtp_preg_medByN), mean(wtp_red_medByN), mean(wtp_pret_medByN), mean(wtp_hosp_medByN)),
    Lower_UI = c(wtp_preg_lowByN, wtp_red_lowByN, wtp_pret_lowByN, wtp_hosp_lowByN),
    Upper_UI = c(wtp_preg_uppByN, wtp_red_uppByN, wtp_pret_uppByN, wtp_hosp_uppByN)
) %>% mutate(
    across(Median:Upper_UI, ~ format(round(.x, 0), big.mark = ","))
)

cat("\nWillingness-to-pay (WTP) for attribute changes in JPY:\n")
print(tot_summary)
cat("\nNote:\n")
cat("1. WTP is calculated as -bAttribute / bCost * 100,000 (JPY)\n")
cat("2. WTP for redness & preterm birth are negative because they are negative attributes\n")
cat("3. Medians are calculated across respondents, then averaged across simulations & imputations\n")
cat("4. 95% UI represents uncertainty in medians across simulations & imputations\n")

# Switch off writing
sink()
cat("\nOutput successfully saved!\n")