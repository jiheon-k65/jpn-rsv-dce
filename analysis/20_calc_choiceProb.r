# ===========================================================================
# Title: 20_calc_choiceProb.r
# Date: 2026-10-28
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

results <- readRDS(here("output", "17.2_hcm_costCont_ranInt_ranCoeff_results.rds"))
est_list <- results$est
robvarcov_list <- results$robvarcov # Robust variance-covariance matrices

# ------------------------------------------------------------------------
# SECTION 3: Define inputs for calculations
# ------------------------------------------------------------------------

# Attribute levels for each alternative
risk <- 0
mv_hosp1.2 <- 0
mv_hosp0.9 <- 1 # 70% efficacy for MV
mv_hosp0.3 <- 0
mv_cost <- 0.24 # 24,000 JPY for MV
mv_redness <- 0
mv_preterm <- 0
mv_preg <- 1 # Pregnant women as recipient for MV
mab_hosp1.2 <- 0
mab_hosp0.9 <- 0
mab_hosp0.3 <- 1 # 90% efficacy for mAb
mab_cost <- 0.45 # 45,000 JPY for mAb
mab_redness <- 0

# Keep only first imputation & unique respondents
df <- data %>% filter(imp == 1) %>% distinct(id, .keep_all = TRUE)
stopifnot(nrow(df) == length(unique(df$id)))

N <- nrow(df) # No. of respondents
M <- max(data$imp) # No. of imputations
S <- 500 # No. of simulations

# ------------------------------------------------------------------------
# SECTION 4: Simulate parameter values
# ------------------------------------------------------------------------

# Obtain respondent characteristics
ageStd <- df$ageStd
sibGrp <- df$sibGrp
eduGrp <- df$eduGrp
incGrp <- df$incGrp

# Initialize matrices to hold probabilities
prob_mv_mat <- matrix(NA, nrow = S, ncol = M)
prob_mab_mat <- matrix(NA, nrow = S, ncol = M)
prob_opt_mat <- matrix(NA, nrow = S, ncol = M)

for (m in 1:M) {
    cat("\n===========================================================================\n")
    cat("\nCalculating choice probabilities for imputation", m, "\n")

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
        bHosp0.9 <- theta[s, "bHosp0.9"]
        bRedness <- theta[s, "bRedness"]
        bHospIntRisk <- theta[s, "bHospIntRisk"]
        bOptIntLV <- theta[s, "bOptIntLV"]

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

        # Simulate bHosp1.2 for each respondent
        mu_hosp1.2 <- theta[s, "mu_hosp1.2"]
        sigma_hosp1.2 <- theta[s, "sigma_hosp1.2"]
        bHosp1.2 <- rnorm(N, mean = mu_hosp1.2, sd = abs(sigma_hosp1.2))

        # Simulate bHosp0.3 for each respondent
        mu_hosp0.3 <- theta[s, "mu_hosp0.3"]
        sigma_hosp0.3 <- theta[s, "sigma_hosp0.3"]
        bHosp0.3 <- rnorm(N, mean = mu_hosp0.3, sd = abs(sigma_hosp0.3))

        # Simulate alpha for each respondent
        sigma_pdt <- theta[s, "sigma_pdt"]
        alpha <- rnorm(N, mean = 0, sd = abs(sigma_pdt))

        # Simulate bOptIntRisk for each respondent
        mu_optIntRisk <- theta[s, "mu_optIntRisk"]
        sigma_optIntRisk <- theta[s, "sigma_optIntRisk"]
        bOptIntRisk <- rnorm(N, mean = mu_optIntRisk, sd = abs(sigma_optIntRisk))

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
        stopifnot(
            length(V_mv) == N,
            length(V_mab) == N,
            length(V_opt) == N
        )

        # Calculate choice probabilities (Dim: N x S x M)
        prob_total <- exp(V_mv) + exp(V_mab) + exp(V_opt)
        prob_mv <- exp(V_mv) / prob_total
        prob_mab <- exp(V_mab) / prob_total
        prob_opt <- exp(V_opt) / prob_total

        # Calculate & store means across respondents (Dim: S x M)
        prob_mv_mat[s, m] <- mean(prob_mv)
        prob_mab_mat[s, m] <- mean(prob_mab)
        prob_opt_mat[s, m] <- mean(prob_opt)
    }
}

# Write output to text file
sink(here("output", "20_calc_choiceProb_output.txt"))
cat("===========================================================================\n")
cat("Title: 20_calc_choiceProb_output.txt\n")
cat("Date: 2026-10-28\n")
cat("===========================================================================\n")

# Calculate means over simulations & imputations
prob_mv_mean <- mean(prob_mv_mat)
prob_mab_mean <- mean(prob_mab_mat)
prob_opt_mean <- mean(prob_opt_mat)

# Calculate 95% uncertainty intervals (UIs)
prob_mv_lowByN <- quantile(as.vector(prob_mv_mat), 0.025)
prob_mv_uppByN <- quantile(as.vector(prob_mv_mat), 0.975)
prob_mab_lowByN <- quantile(as.vector(prob_mab_mat), 0.025)
prob_mab_uppByN <- quantile(as.vector(prob_mab_mat), 0.975)
prob_opt_lowByN <- quantile(as.vector(prob_opt_mat), 0.025)
prob_opt_uppByN <- quantile(as.vector(prob_opt_mat), 0.975)

# Format results into a data frame
tot_summary <- data.frame(
    Alternative = c("Maternal vaccine", "Long-acting mAb", "Opt out"),
    Mean = c(prob_mv_mean, prob_mab_mean, prob_opt_mean),
    Lower_UI = c(prob_mv_lowByN, prob_mab_lowByN, prob_opt_lowByN),
    Upper_UI = c(prob_mv_uppByN, prob_mab_uppByN, prob_opt_uppByN)
) %>% mutate(
    across(Mean:Upper_UI, ~ round(.x, 3))
)
cat("\nAlternative-specific choice probabilities:\n")
print(tot_summary)
cat("\nNote:\n")
cat("1. Means are calculated across respondents, then averaged across simulations & imputations\n")
cat("2. 95% UI represents uncertainty in means across simulations & imputations")

# Switch off writing
sink()
cat("\nOutput successfully saved!\n")