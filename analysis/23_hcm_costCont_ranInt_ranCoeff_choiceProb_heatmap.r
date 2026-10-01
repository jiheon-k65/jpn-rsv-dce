
# ===========================================================================
# Title: 23_hcm_costCont_ranInt_ranCoeff_choiceProb_heatmap.r
# Date: 2026-10-01
# ===========================================================================

# ------------------------------------------------------------------------
# SECTION 1: Load library
# ------------------------------------------------------------------------

# Load packages
library(here)
library(dplyr)
library(ggplot2)

# ------------------------------------------------------------------------
# SECTION 2: Load data & results
# ------------------------------------------------------------------------

df <- read.csv(here("data", "choiceProb_bothCost.csv"))
df$MV_Cost <- df$MV_Cost * 100000 # Convert MV cost to JPY
df$mAb_Cost <- df$mAb_Cost * 100000 # Convert mAb cost to JPY

df_plot <- df %>% mutate(
    Alt = factor(Alt, levels = c("MV", "mAb", "Opt-out"))
)
# ------------------------------------------------------------------------
# SECTION 6: Plot heatmaps of mean choice probabilities
# ------------------------------------------------------------------------

color_manual <- c(
    "MV" = "#4F81BD",
    "mAb" = "#9BBB59",
    "Opt-out" = "#C0504D"
)

df_plot <- df_plot %>%
    mutate(
        tile_color = mapply(
            function(alt, mean) {
                grDevices::colorRampPalette(
                    c("#FFFFFF", color_manual[[alt]])
                )(101)[round(mean * 100) + 1]
            },
            Alt,
            Mean
        )
    )

combined_heatmap <- ggplot(
    df_plot,
    aes(x = MV_Cost, y = mAb_Cost)
) +
    geom_tile(
        aes(fill = tile_color),
        color = "#FFFFFF",
        linewidth = 0.25
    ) +
    geom_text(
        aes(label = sprintf("%.1f", 100 * Mean)),
        size = 1.75
    ) +
    scale_x_continuous(
        breaks = sort(unique(df_plot$MV_Cost)),
        name = "MV cost (JPY)",
        label = scales::comma
    ) +
    scale_y_continuous(
        breaks = sort(unique(df_plot$mAb_Cost)),
        name = "mAb cost (JPY)",
        label = scales::comma
    ) +
    scale_fill_identity() +
    facet_wrap(~Alt, nrow = 1) +
    coord_fixed(ratio = 1, expand = FALSE) +
    theme_bw(base_size = 6) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(
            colour = "#000000",
            fill = NA,
            linewidth = 0.4
        ),
        strip.background = element_rect(
            colour = "#000000",
            fill = "#F2F2F2",
            linewidth = 0.4
        ),
        strip.text = element_text(size = 8),
        legend.key = element_rect(
            fill = NA,
            colour = NA
        )
    )

ggsave(
    filename = here("figs", "choiceProb_all_heatmaps.png"),
    plot = combined_heatmap,
    height = 6.1,
    width = 10,
    units = "cm",
    dpi = 600
)

cat("\nCombined heatmaps successfully saved as choiceProb_all_heatmaps.png\n")