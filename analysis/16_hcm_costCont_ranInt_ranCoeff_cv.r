# ===========================================================================
# Title: 16_hcm_costCont_ranInt_ranCoeff_cv.r
# Date: 2026-09-23
# ===========================================================================

# ------------------------------------------------------------------------
# SECTION 1: Load library
# ------------------------------------------------------------------------

# Load packages
library(here)
library(dplyr)
library(MASS)
library(ggplot2)

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

# Attribute levels for each alternative
risk <- 0
mv_hosp1.2 <- 0
mv_hosp0.9 <- 1 # 70% efficacy for MV
mv_hosp0.3 <- 0
mv_redness <- 0
mv_preterm <- 0
mv_preg <- 1 # Pregnant women as recipient for MV
mab_hosp1.2 <- 0
mab_hosp0.9 <- 0
mab_hosp0.3 <- 1 # 90% efficacy for mAb
mab_redness <- 0

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

# Obtain respondent characteristics
ageStd <- df$ageStd
sibGrp <- df$sibGrp
eduGrp <- df$eduGrp
incGrp <- df$incGrp

# Initialize arrays to hold CVs
CV_mv_array <- array(NA, dim = c(N, S, M))
CV_mab_array <- array(NA, dim = c(N, S, M))

for (m in 1:M) {
    cat("\n===========================================================================\n")
    cat("\nCalculating CVs for imputation", m, "\n")

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
        bHosp1.2 <- theta[s, "bHosp1.2"]
        bHosp0.9 <- theta[s, "bHosp0.9"]
        bHosp0.3 <- theta[s, "bHosp0.3"]
        bRedness <- theta[s, "bRedness"]
        bHospIntRisk <- theta[s, "bHospIntRisk"]
        bOptIntLV <- theta[s, "bOptIntLV"]
        bOptIntRisk <- theta[s, "bOptIntRisk"]

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
            alpha
        )
        V_mab <- (
            bRedness * mab_redness +
            (bHosp1.2 + bHospIntRisk * risk) * mab_hosp1.2 +
            (bHosp0.9 + bHospIntRisk * risk) * mab_hosp0.9 +
            (bHosp0.3 + bHospIntRisk * risk) * mab_hosp0.3 +
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

        # Calculate CV (Dim: N x S x M)
        CV_mv <- - (V_mv - V_opt) / bCost * F
        CV_mab <- - (V_mab - V_opt) / bCost * F

        # Store CVs in arrays (Dim: N x S x M)
        CV_mv_array[, s, m] <- CV_mv
        CV_mab_array[, s, m] <- CV_mab
    }
}

# Write output to text file
sink(here("output", "16_hcm_costCont_ranInt_ranCoeff_cv_output.txt"))
cat("===========================================================================\n")
cat("Title: 16_hcm_costCont_ranInt_ranCoeff_cv_output.txt\n")
cat("Date: 2026-09-23\n")
cat("===========================================================================\n")

# Calculate median across respondents (Dim: S x M)
CV_mv_medByN <- apply(CV_mv_array, c(2, 3), median)
CV_mab_medByN <- apply(CV_mab_array, c(2, 3), median)

# Calculate 95% UI across medians
CV_mv_lowByN <- quantile(CV_mv_medByN, 0.025)
CV_mv_uppByN <- quantile(CV_mv_medByN, 0.975)
CV_mab_lowByN <- quantile(CV_mab_medByN, 0.025)
CV_mab_uppByN <- quantile(CV_mab_medByN, 0.975)

# Format results into a data frame
tot_summary <- data.frame(
    Alternative = c("MV", "mAb"),
    Median = c(mean(CV_mv_medByN), mean(CV_mab_medByN)),
    Lower_UI = c(CV_mv_lowByN, CV_mab_lowByN),
    Upper_UI = c(CV_mv_uppByN, CV_mab_uppByN)
) %>% mutate(
    across(Median:Upper_UI, ~ format(round(.x, 0), big.mark = ","))
)
cat("\nAlternative-specific compensating variations (CVs) in JPY:\n")
print(tot_summary)
cat("\nNote:\n")
cat("1. CV is calculated as - (V_alternative - V_opt) / bCost * 100,000 (JPY)\n")
cat("2. Medians are calculated across respondents, then averaged across simulations & imputations\n")
cat("3. 95% UI represents uncertainty in medians across simulations & imputations\n")

# Switch off writing
sink()
cat("\nOutput successfully saved!\n")