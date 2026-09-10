# ===========================================================================
# Title: 20_hcm_costCont_ranInt_ranCoeff_choiceProb_plot.r
# Date: 2026-09-10
# ===========================================================================

# ------------------------------------------------------------------------
# SECTION 1: Load library
# ------------------------------------------------------------------------

# Load packages
library(here)
library(dplyr)
library(ggplot2)

# ------------------------------------------------------------------------
# SECTION 2: Load results from Monte Carlo simulations
# ------------------------------------------------------------------------

df_mv <- read.csv(here("data", "choiceProb_mvCost.csv"))

df_mab <- read.csv(here("data", "choiceProb_mabCost.csv"))

# ------------------------------------------------------------------------
# SECTION 3: Transform data for plotting
# ------------------------------------------------------------------------

# Format cost & probabilities
df_mv_plot <- df_mv %>%
    mutate(
        Cost = Cost * 100000, # Convert to JPY
        across(Mean:Upper_UI, ~ .x * 100) # Convert to %
    )
df_mab_plot <- df_mab %>%
    mutate(
        Cost = Cost * 100000, # Convert to JPY
        across(Mean:Upper_UI, ~ .x * 100) # Convert to %
    )

# Normalize mean probabilities to sum to 100%
norm_meanProb <- function(df) {
    df %>%
        group_by(Cost) %>%
        mutate(
            Mean = {
                raw_mean <- Mean / sum(Mean) * 100
                rounded_mean <- round(raw_mean, 1)
                correction <- round(100 - sum(rounded_mean), 1)
                rounded_mean[which.max(raw_mean)] <-
                    rounded_mean[which.max(raw_mean)] + correction
                rounded_mean
            }
        ) %>%
        ungroup() %>%
        mutate(across(Lower_UI:Upper_UI, ~ round(.x, 1)))
}

df_mv_plot <- norm_meanProb(df_mv_plot)
df_mab_plot <- norm_meanProb(df_mab_plot)

# ------------------------------------------------------------------------
# SECTION 4: Plot choice probabilities by product cost
# ------------------------------------------------------------------------

color_manual <- c(
    "MV" = "#4F81BD",
    "mAb" = "#9BBB59",
    "Opt-out" = "#C0504D"
)

# Function to plot choice probabilities by product cost
plot_choiceProb <- function(df, base_cost, x_label) {
    # Isolate probabilities for annotation
    prob_c0 <- df[df$Cost == 0, ] # 0 JPY
    prob_base <- df[df$Cost == base_cost, ] # Base cost
    # Add comma separator for plotting
    base_cost_label <- format(base_cost, big.mark = ",")
    # Calculate maximum y-axis value for plotting
    y_max <- ceiling(max(df$Upper_UI) / 5) * 5
    # Plot line plot with uncertainty intervals
    plot <- ggplot(df, aes(x = Cost, y = Mean, color = Alt)) +
        geom_line(linewidth = 0.6) +
        geom_ribbon(
            aes(ymin = Lower_UI, ymax = Upper_UI, fill = Alt),
            alpha = 0.2,
            color = NA
        ) +
        geom_vline(
            xintercept = base_cost,
            linetype = "dashed",
            color = "#404040",
            linewidth = 0.4,
            alpha = 0.6
        ) +
        geom_point(
            data = prob_base,
            aes(x = Cost, y = Mean, color = Alt),
            size = 2,
            shape = 18
        ) +
        geom_text(
            data = prob_c0,
            aes(x = Cost, y = Lower_UI, label = sprintf("%.1f", Mean)),
            vjust = 1.25,
            hjust = 0.65,
            size = 2,
            show.legend = FALSE
        ) +
        geom_text(
            data = prob_base,
            aes(x = Cost, y = Lower_UI, label = sprintf("%.1f", Mean)),
            vjust = 1.25,
            hjust = -0.15,
            size = 2,
            show.legend = FALSE
        ) +
        annotate(
            "text",
            x = base_cost * 0.975,
            y = y_max * 0.995,
            label = paste0("Base cost: ", base_cost_label),
            hjust = 1,
            vjust = 0.5,
            size = 2
        ) +
        labs(
            x = x_label,
            y = "Choice Probability"
        ) +
        theme_bw(base_size = 8) +
        theme(
            legend.position = "top",
            legend.margin = margin(t = -1, b = -5),
            legend.key.size = unit(0.4, "cm"),
            legend.title = element_blank(),
            legend.text = element_text(size = 6)
        ) +
        scale_color_manual(values = color_manual) +
        scale_fill_manual(values = color_manual) +
        scale_x_continuous(labels = scales::comma) +
        coord_cartesian(ylim = c(0, y_max))
    # Return plot
    return(plot)
}

# Plot & save figures
choiceProb_mvCost <- plot_choiceProb(
    df_mv_plot,
    base_cost = 24000,
    x_label = "Maternal Vaccine Cost (JPY)"
)
ggsave(
    here("figs", "choiceProb_mvCost_linePlot.png"),
    plot = choiceProb_mvCost,
    width = 10,
    height = 6,
    units = "cm",
    dpi = 600
)
choiceProb_mabCost <- plot_choiceProb(
    df_mab_plot,
    base_cost = 45000,
    x_label = "Monoclonal Antibody Cost (JPY)"
)
ggsave(
    here("figs", "choiceProb_mabCost_linePlot.png"),
    plot = choiceProb_mabCost,
    width = 10,
    height = 6,
    units = "cm",
    dpi = 600
)
cat("\nPlot successfully saved!\n")