# ===========================================================================
# Title: 21_hcm_costCont_ranInt_ranCoeff_choiceProb_comPlot.r
# Date: 2026-09-10
# ===========================================================================

# ------------------------------------------------------------------------
# SECTION 1: Load library
# ------------------------------------------------------------------------

library(here)
library(dplyr)
library(ggplot2)

# ------------------------------------------------------------------------
# SECTION 2: Load and combine choice-probability results
# ------------------------------------------------------------------------

df_mvCost <- read.csv(here("data", "choiceProb_mvCost.csv"))
df_mabCost <- read.csv(here("data", "choiceProb_mabCost.csv"))

df_plot <- bind_rows(
    df_mvCost %>% mutate(
        CostType = "Change in maternal vaccine cost (JPY)",
        base_cost = 24000
    ),
    df_mabCost %>% mutate(
        CostType = "Change in monoclonal antibody cost (JPY)",
        base_cost = 45000
    )
) %>%
    mutate(
        Cost = Cost * 100000,
        across(Mean:Upper_UI, ~ .x * 100),
        CostType = factor(
            CostType,
            levels = c(
                "Change in maternal vaccine cost (JPY)",
                "Change in monoclonal antibody cost (JPY)"
            )
        )
    )

# Normalize mean probabilities to sum to 100% within each cost curve
norm_meanProb <- function(df) {
    df %>%
        group_by(CostType, Cost) %>%
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

df_plot <- norm_meanProb(df_plot)

base_cases <- df_plot %>%
    distinct(CostType, base_cost) %>%
    mutate(
        label = paste0("Base cost: ", format(base_cost, big.mark = ","))
    )

y_max <- ceiling(max(df_plot$Upper_UI, na.rm = TRUE) / 5) * 5

# ------------------------------------------------------------------------
# SECTION 3: Create combined figure
# ------------------------------------------------------------------------

color_manual <- c(
    "MV" = "#4F81BD",
    "mAb" = "#9BBB59",
    "Opt-out" = "#C0504D"
)

prob_c0 <- df_plot %>%
    filter(Cost == 0)

prob_base <- df_plot %>%
    inner_join(
        base_cases[, c("CostType", "base_cost")],
        by = c("CostType", "Cost" = "base_cost")
    )

choiceProb_combined <- ggplot(
    df_plot,
    aes(x = Cost, y = Mean, color = Alt, group = Alt)
) +
    geom_vline(
        data = base_cases,
        aes(xintercept = base_cost),
        linetype = "dashed",
        color = "#404040",
        linewidth = 0.4,
        alpha = 0.6
    ) +
    geom_ribbon(
        aes(ymin = Lower_UI, ymax = Upper_UI, fill = Alt),
        alpha = 0.2,
        color = NA
    ) +
    geom_line(linewidth = 0.6) +
    geom_point(
        data = df_plot %>%
            semi_join(base_cases, by = c("CostType", "Cost" = "base_cost")),
        aes(),
        size = 2,
        shape = 18
    ) +
    geom_text(
        data = prob_c0,
        aes(
            x = Cost,
            y = Lower_UI,
            label = sprintf("%.1f", Mean)
        ),
        vjust = 1.25,
        hjust = 0.65,
        size = 2,
        show.legend = FALSE
    ) +
    geom_text(
        data = prob_base,
        aes(
            x = Cost,
            y = Lower_UI,
            label = sprintf("%.1f", Mean)
        ),
        vjust = 1.25,
        hjust = -0.15,
        size = 2,
        show.legend = FALSE
    ) +
    geom_text(
        data = base_cases,
        aes(
            x = base_cost,
            y = y_max * 0.995,
            label = label
        ),
        hjust = 1.05,
        vjust = 0.5,
        size = 2.5,
        color = "#404040",
        inherit.aes = FALSE
    ) +
    facet_wrap(
        ~ CostType,
        ncol = 1,
        scales = "free_x",
        strip.position = "top"
    ) +
    labs(
        x = NULL,
        y = "Choice Probability"
    ) +
    scale_color_manual(values = color_manual) +
    scale_fill_manual(values = color_manual) +
    scale_x_continuous(labels = scales::comma) +
    coord_cartesian(ylim = c(0, y_max)) +
    theme_bw(base_size = 8) +
    theme(
        legend.position = "top",
        legend.margin = margin(t = -1, b = -5),
        legend.key.size = grid::unit(0.4, "cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        strip.background = element_rect(fill = "grey95", color = "grey60"),
        strip.text = element_text(size = 8),
        panel.spacing.y = grid::unit(0.2, "lines")
    )

# ------------------------------------------------------------------------
# SECTION 4: Save figure
# ------------------------------------------------------------------------

ggsave(
    here("figs", "choiceProb_combCost_linePlot.png"),
    plot = choiceProb_combined,
    width = 10,
    height = 12,
    units = "cm",
    dpi = 600
)
cat("\nPlot successfully saved!\n")
