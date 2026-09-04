# ===========================================================================
# Title: 21_hcm_costCont_ranCoeff_choiceProb.r
# Date: 2026-09-04
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

# Attribute levels for each alternative
risk <- 0
mv_hosp1.2 <- 0 # 50% efficacy for MV
mv_hosp0.9 <- 0
mv_hosp0.3 <- 0
mv_cost <- 0.24 # 24,000 JPY for MV
mv_redness <- 0
mv_preterm <- 0
mv_preg <- 0
mab_hosp1.2 <- 0
mab_hosp0.9 <- 1 # 70% efficacy for mAb
mab_hosp0.3 <- 0
mab_cost <- 0.45 # 45,000 JPY for mAb
mab_redness <- 0

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
prob_mv_array <- array(NA, dim = c(N, S, M))
prob_mab_array <- array(NA, dim = c(N, S, M))
prob_opt_array <- array(NA, dim = c(N, S, M))

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
        # Simulate fixed coefficients
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
        eta <- rnorm(N, mean = 0, sd = 1) # Random distur   xsdbance
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

        # Calculate choice probabilities for each alternative
        denom <- exp(V_mv) + exp(V_mab) + exp(V_opt)
        prob_mv <- exp(V_mv) / denom
        prob_mab <- exp(V_mab) / denom
        prob_opt <- exp(V_opt) / denom

        # Store probabilities in arrays
        prob_mv_array[, s, m] <- prob_mv
        prob_mab_array[, s, m] <- prob_mab
        prob_opt_array[, s, m] <- prob_opt
    }
}

# Write output to text file
sink(here("output", "21_hcm_costCont_ranCoeff_choiceProb_output.txt"))
cat("===========================================================================\n")
cat("Title: 21_hcm_costCont_ranCoeff_choiceProb_output.txt\n")
cat("Date: 2026-09-04\n")
cat("===========================================================================\n")

# Calculate means across respondents (Dim: S x M)
# \bar{P}_j^{(m,s)} = E_n[P_{j,n}^{(m,s)}]
prob_mv_mean_byN <- apply(prob_mv_array, c(2, 3), mean)
prob_mab_mean_byN <- apply(prob_mab_array, c(2, 3), mean)
prob_opt_mean_byN <- apply(prob_opt_array, c(2, 3), mean)

# Flatten over simulations & imputations
# \bar{P}_j = E_{m,s}[\bar{P}_j^{(m,s)}]
prob_mv_mean <- mean(prob_mv_mean_byN)
prob_mab_mean <- mean(prob_mab_mean_byN)
prob_opt_mean <- mean(prob_opt_mean_byN)

# Calculate 95% uncertainty intervals (UIs)
# [Q_{0.025}(\bar{P}_j^{(m,s)}, Q_{0.975}(\bar{P}_j^{(m,s)}]
prob_mv_low_byN <- quantile(as.vector(prob_mv_mean_byN), 0.025)
prob_mv_upp_byN <- quantile(as.vector(prob_mv_mean_byN), 0.975)
prob_mab_low_byN <- quantile(as.vector(prob_mab_mean_byN), 0.025)
prob_mab_upp_byN <- quantile(as.vector(prob_mab_mean_byN), 0.975)
prob_opt_low_byN <- quantile(as.vector(prob_opt_mean_byN), 0.025)
prob_opt_upp_byN <- quantile(as.vector(prob_opt_mean_byN), 0.975)

# Format results into a data frame
tot_summary <- data.frame(
    Alternative = c("MV", "mAb", "Opt"),
    Mean = c(prob_mv_mean, prob_mab_mean, prob_opt_mean),
    Lower_UI = c(prob_mv_low_byN, prob_mab_low_byN, prob_opt_low_byN),
    Upper_UI = c(prob_mv_upp_byN, prob_mab_upp_byN, prob_opt_upp_byN)
) %>% mutate(
    across(Mean:Upper_UI, ~ round(.x, 5))
)
cat("\nAlternative-specific choice probabilities:\n")
print(tot_summary)

# Calculate variance across respondents
# E_{m,s}[Var_n[P_{j,n}^{(m,s)}]]
prob_mv_var_byN <- mean(apply(prob_mv_array, c(2, 3), var))
prob_mab_var_byN <- mean(apply(prob_mab_array, c(2, 3), var))
prob_opt_var_byN <- mean(apply(prob_opt_array, c(2, 3), var))

# Calculate variance across simulations
# E_{m,n}[Var_s[P_{j,n}^{(m,s)}]]
prob_mv_var_byS <- mean(apply(prob_mv_array, c(1, 3), var))
prob_mab_var_byS <- mean(apply(prob_mab_array, c(1, 3), var))
prob_opt_var_byS <- mean(apply(prob_opt_array, c(1, 3), var))

# Calculate variance across imputations
# E_{n,s}[Var_m[P_{j,n}^{(m,s)}]]
prob_mv_var_byM <- mean(apply(prob_mv_array, c(1, 2), var))
prob_mab_var_byM <- mean(apply(prob_mab_array, c(1, 2), var))
prob_opt_var_byM <- mean(apply(prob_opt_array, c(1, 2), var))

# Format results into a data frame
var_summary <- data.frame(
    Alternative = c("MV", "mAb", "Opt"),
    Var_byN = c(prob_mv_var_byN, prob_mab_var_byN, prob_opt_var_byN),
    Var_byS = c(prob_mv_var_byS, prob_mab_var_byS, prob_opt_var_byS),
    Var_byM = c(prob_mv_var_byM, prob_mab_var_byM, prob_opt_var_byM)
) %>% mutate(
    across(Var_byN:Var_byM, ~ round(.x, 5))
)
cat("\nAlternative-specific choice probability variances:\n")
print(var_summary)
cat("\nVar_byN: Preference heterogeneity\n")
cat("Var_byS: Parameter uncertainty\n")
cat("Var_byM: Imputation uncertainty\n")

# Switch off writing
sink()
