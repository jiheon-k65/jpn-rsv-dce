# ===========================================================================
# Title: 6_mnl_costContQuad.r
# Date: 2026-07-20
# ===========================================================================

# ------------------------------------------------------------------------
# SECTION 1: Load library & define core setting
# ------------------------------------------------------------------------

# Clear memory
rm(list = ls())

# Load packages
library(here)
library(apollo)
library(ggplot2)

# Initialize code
apollo_initialise()

# Set core controls
apollo_control <- list(
    modelDescr = "MNL; quadratic cost",
    modelName = "6_mnl_costContQuad",
    indivID = "id",
    outputDirectory = here("output")
)

# ------------------------------------------------------------------------
# SECTION 2: Load data & apply transformations
# ------------------------------------------------------------------------

# Load data from WD
data <- read.csv(here("data", "jpn_dataset_subset.csv"), header = TRUE)

# Exclude repeated task
database <- data[data$rptTask != 1, ]

# ------------------------------------------------------------------------
# SECTION 3: Define model parameters
# ------------------------------------------------------------------------

# List all parameters & starting values
apollo_beta <- c(
    ascPdt = 0,
    ascOpt = 0, # Alternative-specific constant for opt-out
    bInf_mab = 0, # Alternative-specific coefficient
    bPreg_mv = 0, # Alternative-specific coefficient
    bNoAe = 0,
    bRedness = 0,
    bPreterm_mv = 0, # Alternative-specific coefficient
    bHosp1.5 = 0,
    bHosp1.2 = 0,
    bHosp0.9 = 0,
    bHosp0.3 = 0,
    bCost = 0, # Linear cost coefficient (per 100,000 JPY)
    bCostSq = 0, # Quadratic cost coefficient (per 100,000 JPY)
    bOptIntRisk = 0,
    bHospIntRisk = 0
)
# List parameters to be fixed at their starting values
apollo_fixed <- c("ascPdt", "bInf_mab", "bNoAe", "bHosp1.5")

# ------------------------------------------------------------------------
# SECTION 4: Group & validate inputs
# ------------------------------------------------------------------------

apollo_inputs <- apollo_validateInputs()

# ------------------------------------------------------------------------
# SECTION 5: Define model & likelihood function
# ------------------------------------------------------------------------

apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
    # Attach inputs & detach after function exit
    apollo_attach(apollo_beta, apollo_inputs)
    on.exit(apollo_detach(apollo_beta, apollo_inputs))

    # Create list of probabilities
    P <- list()

    # List utilities
    V <- list()
    V[["mv"]] <- (
        bPreg_mv * mv_preg +
        bRedness * mv_redness +
        bPreterm_mv * mv_preterm +
        (bHosp1.2 + bHospIntRisk * risk) * mv_hosp1.2 +
        (bHosp0.9 + bHospIntRisk * risk) * mv_hosp0.9 +
        (bHosp0.3 + bHospIntRisk * risk) * mv_hosp0.3 +
        bCost * mv_cost +
        bCostSq * (mv_cost * mv_cost)
    )
    V[["mab"]] <- (
        bRedness * mab_redness +
        (bHosp1.2 + bHospIntRisk * risk) * mab_hosp1.2 +
        (bHosp0.9 + bHospIntRisk * risk) * mab_hosp0.9 +
        (bHosp0.3 + bHospIntRisk * risk) * mab_hosp0.3 +
        bCost * mab_cost +
        bCostSq * (mab_cost * mab_cost)
    )
    V[["opt"]] <- ascOpt + bOptIntRisk * risk

    # Define MNL model settings
    mnl_settings <- list(
        alternatives = c(mv = 1, mab = 2, opt = 3),
        avail = 1, # All alternatives available
        choiceVar = choice_stage2,
        utilities = V
    )

    # Compute probabilities using MNL model
    P[["model"]] <- apollo_mnl(mnl_settings, functionality)

    # Take product across observation for same individual
    P <- apollo_panelProd(P, apollo_inputs, functionality)

    # Prepare & return outputs of function
    P <- apollo_prepareProb(P, apollo_inputs, functionality)
    return(P)
}

# ------------------------------------------------------------------------
# SECTION 6: Estimate model
# ------------------------------------------------------------------------

# Estimate model
fit <- apollo_estimate(
    apollo_beta,
    apollo_fixed,
    apollo_probabilities,
    apollo_inputs
)

# ------------------------------------------------------------------------
# SECTION 7: Output model results
# ------------------------------------------------------------------------

# Show 2-sided p values for coefficients
modelOutput_settings <- list(printPVal = 2)
saveOutput_settings <- list(printPVal = 2)

# Output model to terminal
apollo_modelOutput(fit, modelOutput_settings)

# Save model output to WD
apollo_saveOutput(fit, saveOutput_settings)

# ------------------------------------------------------------------------
# SECTION 8: Plot model coefficients
# ------------------------------------------------------------------------

# Read results
result <- read.csv(here("output", "6_mnl_costContQuad_estimates.csv"))
bCost <- result$Estimate[result$X == "bCost"]
bCostSq <- result$Estimate[result$X == "bCostSq"]

quad_df <- data.frame(x = seq(0, 1, by = 0.001))
quad_df$V <- bCost * quad_df$x + bCostSq * (quad_df$x ^ 2)

pts_df <- data.frame(x = c(0, 0.05, 0.1, 0.3, 1))
pts_df$V <- bCost * pts_df$x + bCostSq * (pts_df$x ^ 2)

# Minimum point
x_min <- -bCost / (2 * bCostSq)
V_min <- bCost * x_min + bCostSq * (x_min ^ 2)

# Plot quadratic utility function
plot <- ggplot(quad_df, aes(x = x, y = V)) +
    geom_line(color = "#9BBB59", linewidth = 0.6) +
    geom_point(
        data = pts_df,
        aes(x = x, y = V),
        color = "#9BBB59",
        size = 1.2
    ) +
    geom_vline(
        xintercept = x_min,
        linetype = "dashed",
        color = "#C0504D",
        linewidth = 0.4,
        alpha = 0.8
    ) +
    geom_point(
        aes(x = x_min, y = V_min),
        color = "#C0504D",
        size = 1.2
    ) +
    annotate(
        "text",
        x = x_min - 0.02,
        y = -0.25,
        label = paste0("Min: (", round(x_min, 2), ", ", round(V_min, 2), ")"),
        hjust = 1,
        vjust = 0.5,
        size = 2
    ) +
    labs(x = "Cost (per 100,000 JPY)", y = "Utility") +
    theme_bw(base_size = 6)

ggsave(
    here("figs", "vCostQuad_linePlot.png"), 
    plot = plot, height = 4, width = 6, units = "cm", dpi = 600
)