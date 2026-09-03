# ===========================================================================
# Title: 23_hcm_costCont_ranCoeff_choiceProb_mabCost.r
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

# Attribute levels for each alternative
risk <- 0
mv_hosp1.2 <- 0
mv_hosp0.9 <- 1
mv_hosp0.3 <- 0
mv_cost <- 0.24 # 24,000 JPY for MV
mv_redness <- 0
mv_preterm <- 0
mv_preg <- 0
mab_hosp1.2 <- 0
mab_hosp0.9 <- 0
mab_hosp0.3 <- 1
mab_cost_base <- 0.45 # 45,000 JPY for mAb
mab_redness <- 0

# Cost grid for MV
cost_grid <- seq(0, 0.7, by = 0.07) # 11 cost values
cost_grid <- sort(c(cost_grid, mab_cost_base)) # Add base case cost

M <- 25 # No. of imputations
N <- nrow(data) / M # No. of respondents
S <- 200 # No. of simulations
K <- length(cost_grid) # No. of cost values for MV

# ------------------------------------------------------------------------
# SECTION 4: Simulate parameter values
# ------------------------------------------------------------------------

# Obtain respondent characteristics
df <- data %>% filter(imp == 1) # Same for all imputations
ageStd <- df$ageStd
sibGrp <- df$sibGrp
eduGrp <- df$eduGrp
incGrp <- df$incGrp

# Initialize array to hold probabilities
prob_mv_array <- array(NA, dim = c(N, S, M, K))
prob_mab_array <- array(NA, dim = c(N, S, M, K))
prob_opt_array <- array(NA, dim = c(N, S, M, K))

for (m in 1:M) {
    cat("\n===========================================================================\n")
    cat("\nCalculating probabilities for imputation", m, "\n")

    # Extract estimates & variance-covariance matrix
    theta_hat <- est_list[[m]]
    sigma <- varcov_list[[m]]

    # Draw parameter values from multivariate normal distribution
    set.seed(123 + m)
    theta <- mvrnorm(n = S, mu = theta_hat, Sigma = sigma)
    colnames(theta) <- names(theta_hat)

    # Simulate random coefficients
    for (s in 1:S) {
        # Isolate fixed coefficients
        bHosp1.2 <- theta[s, "bHosp1.2"]
        bHosp0.9 <- theta[s, "bHosp0.9"]
        bHosp0.3 <- theta[s, "bHosp0.3"]
        bRedness <- theta[s, "bRedness"]
        bHospIntRisk <- theta[s, "bHospIntRisk"]
        bOptIntLV <- theta[s, "bOptIntLV"]
        bOptIntRisk <- theta[s, "bOptIntRisk"]

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

        # Simulate alpha for each respondent
        sigma_pdt <- theta[s, "sigma_pdt"]
        alpha <- rnorm(N, mean = 0, sd = abs(sigma_pdt))

        # Simulate LV for each respondent
        gamma_age <- theta[s, "gamma_age"]
        gamma_sib <- theta[s, "gamma_sib"]
        gamma_edu <- theta[s, "gamma_edu"]
        gamma_inc <- theta[s, "gamma_inc"]
        eta <- rnorm(N, mean = 0, sd = 1) # Random disturbance
        LV <- (
            gamma_age * ageStd +
            gamma_sib * sibGrp +
            gamma_edu * eduGrp +
            gamma_inc * incGrp +
            eta
        )
        stopifnot(length(LV) == N)

        for (k in seq_along(cost_grid)) {
            # Update mab_cost for each iteration
            mab_cost <- cost_grid[k]

            # Calculate utilities for each alternative
            V_mv <- (
                bPreg_mv * mv_preg +
                bRedness * mv_redness +
                bPreterm_mv * mv_preterm +
                (bHosp1.2 + bHospIntRisk * risk) * mv_hosp1.2 +
                (bHosp0.9 + bHospIntRisk * risk) * mv_hosp0.9 +
                (bHosp0.3 + bHospIntRisk * risk) * mv_hosp0.3 +
                bCost * mv_cost +
                alpha
            )
            V_mab <- (
                bRedness * mab_redness +
                (bHosp1.2 + bHospIntRisk * risk) * mab_hosp1.2 +
                (bHosp0.9 + bHospIntRisk * risk) * mab_hosp0.9 +
                (bHosp0.3 + bHospIntRisk * risk) * mab_hosp0.3 +
                bCost * mab_cost +
                alpha
            )
            V_opt <- (
                bOptIntRisk * risk +
                bOptIntLV* LV
            )

            # Calculate choice probabilities for each alternative
            denom <- exp(V_mv) + exp(V_mab) + exp(V_opt)
            prob_mv <- exp(V_mv) / denom
            prob_mab <- exp(V_mab) / denom
            prob_opt <- exp(V_opt) / denom

            # Store probabilities in arrays
            prob_mv_array[, s, m, k] <- prob_mv
            prob_mab_array[, s, m, k] <- prob_mab
            prob_opt_array[, s, m, k] <- prob_opt
        }
    }
}

# Initialize list to hold results for each cost value
summary_list <- vector("list", length = length(cost_grid))

for (k in seq_along(cost_grid)) {
    cat("\n===========================================================================\n")
    cat("\nCalculating mean & 95% UI for cost value", k, "\n")

    # Calculate means across respondents
    prob_mv_mean_byN <- apply(prob_mv_array[, , , k], c(2, 3), mean)
    prob_mab_mean_byN <- apply(prob_mab_array[, , , k], c(2, 3), mean)
    prob_opt_mean_byN <- apply(prob_opt_array[, , , k], c(2, 3), mean)

    summary_list[[k]] <- data.frame(
        Cost = cost_grid[k],
        Alt = c("MV", "mAb", "Opt-out"),
        Mean = c(
            mean(prob_mv_mean_byN),
            mean(prob_mab_mean_byN),
            mean(prob_opt_mean_byN)
        ),
        Lower_UI = c(
            quantile(as.vector(prob_mv_mean_byN), 0.025),
            quantile(as.vector(prob_mab_mean_byN), 0.025),
            quantile(as.vector(prob_opt_mean_byN), 0.025)
        ),
        Upper_UI = c(
            quantile(as.vector(prob_mv_mean_byN), 0.975),
            quantile(as.vector(prob_mab_mean_byN), 0.975),
            quantile(as.vector(prob_opt_mean_byN), 0.975)
        )
    )
}

# Bind results into a single data frame
summary_combined <- bind_rows(summary_list) %>%
    mutate(across(Mean:Upper_UI, ~ round(.x, 5)))

# Save results as CSV
write.csv(
    summary_combined,
    file = here("data", "choiceProb_mabCost.csv"),
    row.names = FALSE
)
cat("\nFile successfully saved as choiceProb_mabCost.csv\n")