# ===========================================================================
# Title: 21.1_hcm_costCont_ranCoeff_choiceProb_byMvCost.r
# Date: 2026-07-23
# ===========================================================================

# ------------------------------------------------------------------------
# SECTION 1: Load library
# ------------------------------------------------------------------------

# Load packages
library(here)
library(ggplot2)

# ------------------------------------------------------------------------
# SECTION 2: Load results from Monte Carlo simulations
# ------------------------------------------------------------------------

df_mvCost <- read.csv(here("data", "choiceProb_mvCost.csv"))

df_mabCost <- read.csv(here("data", "choiceProb_mabCost.csv"))

# ------------------------------------------------------------------------
# SECTION 3: Plot & save results
# ------------------------------------------------------------------------

color_manual <- c(
    "MV" = "#4F81BD",
    "mAb" = "#9BBB59",
    "Opt-out" = "#C0504D"
)

# Function to plot choice probabilities by product cost
plot_choiceProb <- function(df, base_cost, x_label) {
    # Format dataframe for plotting
    df_plot <- df %>%
        mutate(
            Cost = Cost * 100000, # Convert to JPY
            across(Mean:Upper_UI, ~ .x * 100)  # Converts to %
        )
    # Isolate means from base case for plotting
    base_means <- df_plot[df_plot$Cost == base_cost, ]
    # Add comma separator for plotting
    base_cost_label <- format(base_cost, big.mark = ",")
    # Calculate maximum y-axis value for plotting
    y_max <- ceiling(max(df_plot$Upper_UI) / 5) * 5
    # Plot line plot with uncertainty intervals
    plot <- ggplot(df_plot, aes(x = Cost, y = Mean, color = Alt)) +
        geom_line(size = 0.6) +
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
            alpha = 0.8
        ) +
        geom_point(
            data = base_means,
            aes(x = Cost, y = Mean, color = Alt),
            size = 2,
            shape = 18
        ) +
        annotate(
            "text",
            x = base_cost * 0.975,
            y = y_max * 0.995,
            label = paste0("Base case: ", base_cost_label),
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
    df_mvCost,
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
    df_mabCost,
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