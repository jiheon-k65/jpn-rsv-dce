# ===========================================================================
# Title: jpn_universal_func.r
# Date: 2026-06-18
# ===========================================================================

# ------------------------------------------------------------------------
# SECTION 1: Functions for data wrangling
# ------------------------------------------------------------------------

# Create dummy columns for categorical factors
create_dummy_cols <- function(x) {
    # Obtain col name
    col_name <- deparse(substitute(x))
    # Obtain no. of levels excluding reference level
    n_dummy_cols <- length(unique(x[x != 0 & !is.na(x)]))
    # Initialize new matrix with 0s
    new_mat <- matrix(0, nrow = length(x), ncol = n_dummy_cols)
    # Insert 1s into the appropriate column
    for (i in 1:length(x)) {
        if (is.na(x[i])) {
            new_mat[i, ] <- NA
        } else if (x[i] != 0) {
            new_mat[i, x[i]] <- 1
        }
    }
    cat("\nDummy columns created for", col_name, "\n")
    return(new_mat)
}

# Add prefix to column names
add_prefix <- function(df, cols, prefix) {
    idx <- colnames(df) %in% cols
    colnames(df)[idx] <- paste0(prefix, colnames(df)[idx])
    cat("\nAdded prefix ", prefix, " to columns: ", paste(cols, collapse = ", "), "\n", sep = "")
    return(df)
}

# Indicate no. of women that passed logical test
indicate_passTest <- function(df) {
    # Obtain df name
    df_name <- deparse(substitute(df))
    # Indicate if answers to repeated questions are the same
    df$passTest <- ifelse(
        df$block <= 10,
        as.integer(df$Q1S1 == df$Q4S1), # Qns 1 & 4 for blocks 1-10
        as.integer(df$Q6S1 == df$Q9S1) # Qns 6 & 9 for blocks 11-20
    )
    num_passes <- sum(df$passTest == 1)
    cat("\nAdded column for logical test results: passTest\n")
    cat("No. of women that passed logical test:", num_passes, "\n")
    return(df)
}

# Indicate no. of women with same responses
indicate_allSame <- function(df) {
    # Obtain df name
    df_name <- deparse(substitute(df))
    # Create vector of column names
    col_names <- paste0("Q", 1:10, "S1")
    # Count unique no. of values across each row
    unique_counts <- apply(df[, col_names], 1, function(x) length(unique(x)))
    # Indicate allSame = 1 if unique_counts = 1
    df$allSame <- as.integer(unique_counts == 1)
    num_same <- sum(df$allSame == 1)
    cat("\nAdded column for same responses: allSame\n")
    cat("No. of women with same responses:", num_same, "\n")
    return(df)
}

# Indicate repeated task
indicate_rptTask <- function(df) {
    # Obtain df name
    df_name <- deparse(substitute(df))
    # Initialize new column as 0s
    df$rptTask <- 0
    # Loop through each "id"
    for (i in unique(df$id)) {
        sub <- df[df$id == i, ]
        if (all(sub$block <= 10)) {
            sub$rptTask <- ifelse(sub$task == 4, 1, 0) # Qn 4 for blocks 1-10
        } else {
            sub$rptTask <- ifelse(sub$task == 9, 1, 0) # Qn 9 for blocks 11-20
        }
        df$rptTask[df$id == i] <- sub$rptTask
    }
    cat("\nAdded column for repeated task: rptTask\n")
    return(df)
}

# Reassign columns as MV (left) & mAb (right)
reassign_columns <- function(df) {
    # Obtain df name
    df_name <- deparse(substitute(df))
    # Identify & sort "alt1" & "alt2" columns
    alt1_cols <- grep("^alt1_", colnames(df), value = TRUE) %>% sort()
    alt2_cols <- grep("^alt2_", colnames(df), value = TRUE) %>% sort()
    # Indicate rows that require swapping
    swap_rows <- (df$alt1_preg == 0)
    # Replace "alt1" columns with "alt2" columns & vice versa for those rows
    df[swap_rows, c(alt1_cols, alt2_cols)] <- df[swap_rows, c(alt2_cols, alt1_cols)]
    # Update "choice_stage1" column (1 = MV, 2 = mAb)
    df[swap_rows, "choice_stage1"] <- ifelse(df[swap_rows, "choice_stage1"] == 1, 2, 1)
    # Update "choice_stage2" column (1 = MV, 2 = mAb, 3 = opt-out)
    df$choice_stage2 <- ifelse(df$choice_stage2 == 2, 3, df$choice_stage1)
    # Drop "preg" & "preterm" columns for "alt2" (mab)
    df$alt2_preg <- NULL
    df$alt2_preterm <- NULL
    # Update column prefixes
    colnames(df) <- sub("^alt1_", "mv_", colnames(df))
    colnames(df) <- sub("^alt2_", "mab_", colnames(df))
    # Print message & return df
    cat("\nReassigned columns as mv (left) & mab (right)\n")
    cat("Dropped mab_preg & mab_preterm columns as not applicable\n")
    return(df)
}

# Indicate product availability in stage 2
indicate_prodAv <- function(df) {
    # Obtain df name
    df_name <- deparse(substitute(df))
    # Create new columns
    df$mv_av <- ifelse(df$choice_stage1 == 1, 1, 0) # MV chosen
    df$mab_av <- ifelse(df$choice_stage1 == 2, 1, 0) # mAb chosen
    # Add "av_opt" for 2nd stage
    df$opt_av <- 1
    # Print message & return df
    cat("\nAdded columns for product availability: mv_av, mab_av, opt_av\n")
    return(df)
}

# Indicate missing as NA
indicate_NA <- function(df, cols, na_val) {
    for (col in cols) {
        df[[col]][df[[col]] == na_val] <- NA
    }
    cat("\nIndicated NA for columns: ", paste(cols, collapse = ", "), "\n", sep = "")
    return(df)
}

# ------------------------------------------------------------------------
# SECTION 2: Functions for data exploration
# ------------------------------------------------------------------------

# Query time at multple cutoffs
query_time <- function(col) {
    counts <- c(
        "<3" = sum(col < 3),
        "3-<4" = sum(col >= 3 & col < 4),
        "4-<5" = sum(col >= 4 & col < 5),
        ">5" = sum(col > 5)
    ) / 10 # No. of rows per respondent
    cat("\nNo. of women by survey completion time (min):\n")
    print(counts)
    invisible(counts)
}

# Query age at multiple cutoffs
query_age <- function(col) {
    counts <- c(
        "<20" = sum(col < 20),
        "20-<30" = sum(col >= 20 & col < 30),
        "30-<40" = sum(col >= 30 & col < 40),
        "40-<50" = sum(col >= 40 & col < 50),
        "50-<60" = sum(col >= 50 & col < 60),
        ">60" = sum(col > 60)
    ) / 10 # No. of rows per respondent
    cat("\nNo. of women by age group:\n")
    print(counts)
    invisible(counts)
}

# Query gestational age at multiple cutoffs
query_gestAge <- function(col) {
    counts <- c(
        "<=4" = sum(col <= 4),
        ">4-42" = sum(col > 4 & col <= 42),
        ">42" = sum(col > 42)
    ) / 10 # No. of rows per respondent
    cat("\nNo. of women by gestational age (week):\n")
    print(counts)
    invisible(counts)
}

# Summarize categorical variables
summ_cat <- function(df, cols) {
    # Keep first observation for each "id"
    df_person <- df %>%
        group_by(id) %>%
        summarise(
            across(all_of(cols), first),
            .groups = "drop"
        ) %>% as.data.frame()
    # Summarize categorical variable(s)
    summary <- lapply(cols, function(col) {
        tbl <- table(df_person[[col]], useNA = "ifany") # Include NA
        prop <- prop.table(tbl) * 100
        data.frame(
            Variable = col,
            Category = names(tbl),
            Count = as.vector(tbl),
            Percentage = sprintf("%.2f", as.vector(prop))
        )
    }) %>% bind_rows() # Combine into one dataframe
    cat("\nSummary of categorical variable(s):\n")
    print(summary, row.names = FALSE)
    invisible(summary)
}

# Summarize continuous variables
summ_cont <- function(df, cols, iqr = NULL) {
    # Keep first observation for each "id"
    df_person <- df %>%
        group_by(id) %>%
        summarise(
            across(all_of(cols), first),
            .groups = "drop"
        ) %>% as.data.frame()
    # Summarize continuous variable(s)
    if (!is.null(iqr)) {
        summary <- lapply(cols, function(col) {
            data.frame(
                Variable = col,
                Count = sum(!is.na(df_person[[col]])),
                Median = sprintf("%.2f", median(df_person[[col]], na.rm = TRUE)),
                Q1 = sprintf("%.2f", quantile(df_person[[col]], 0.25, na.rm = TRUE)),
                Q3 = sprintf("%.2f", quantile(df_person[[col]], 0.75, na.rm = TRUE)),
                Min = sprintf("%.2f", min(df_person[[col]], na.rm = TRUE)),
                Max = sprintf("%.2f", max(df_person[[col]], na.rm = TRUE))
            )
        }) %>% bind_rows() # Combine into one dataframe
    } else {
        summary <- lapply(cols, function(col) {
            data.frame(
                Variable = col,
                Count = sum(!is.na(df_person[[col]])),
                Mean = sprintf("%.2f", mean(df_person[[col]], na.rm = TRUE)),
                SD = sprintf("%.2f", sd(df_person[[col]], na.rm = TRUE)),
                Min = sprintf("%.2f", min(df_person[[col]], na.rm = TRUE)),
                Max = sprintf("%.2f", max(df_person[[col]], na.rm = TRUE))
            )
        }) %>% bind_rows() # Combine into one dataframe
    }
    cat("\nSummary of continuous variable(s):\n")
    print(summary, row.names = FALSE)
    invisible(summary)
}

# Summarize categorical variable by subgroup
summ_cat_bySub <- function(df, col, grp_col) {
    # Keep first observation for each "id"
    df_person <- df %>%
        group_by(id) %>%
        summarise(
            value = first(.data[[col]]),
            group = first(.data[[grp_col]]),
            .groups = "drop"
        ) %>% as.data.frame()
    # Summarize categorical variable by subgroup
    crosstab <- table(
        setNames(
            list(df_person$value, df_person$group),
            c(col, grp_col)
        )
    )
    percentages <- sprintf("%.2f", prop.table(crosstab, margin = 1) * 100)
    summary <- matrix(
        paste0(crosstab, " (", percentages, ")"),
        nrow = nrow(crosstab),
        ncol = ncol(crosstab),
        dimnames = dimnames(crosstab)
    )
    # Calculate p-value using chi-squared / Fisher's exact test
    chi <- suppressWarnings(chisq.test(crosstab))
    use_fisher <- any(chi$expected < 5)
    if (use_fisher) {
        test <- fisher.test(crosstab, workspace = 2e8)
        p <- sprintf("%.3f", test$p.value)
        test_name <- "Fisher's exact test"
    }
    else {
        p <- sprintf("%.3f", chi$p.value)
        test_name <- "chi-squared test"
    }
    names(dimnames(summary)) <- names(dimnames(crosstab))
    title <- glue("Crosstabulation of {col} by {grp_col}")
    cat("\n", title, ":\n", sep = "")
    print(noquote(summary))
    cat("p-value (", test_name, "): ", ifelse(p < 0.001, "<0.001", p), "\n", sep = "")
    invisible(summary)
}

# Summarize continuous variables by subgroup
summ_cont_bySub <- function(df, col, grp_col) {
    # Keep first observation for each "id"
    df_person <- df %>%
        group_by(id) %>%
        summarise(
            value = first(.data[[col]]),
            group = first(.data[[grp_col]]),
            .groups = "drop"
        ) %>% as.data.frame()
    # Summarize continuous variable by subgroup
    summary <- df_person %>%
        group_by(group) %>%
             summarise(
                Count = sum(!is.na(value)),
                Mean = sprintf("%.2f", mean(value, na.rm = TRUE)),
                SD = sprintf("%.2f", sd(value, na.rm = TRUE)),
                Min = sprintf("%.2f", min(value, na.rm = TRUE)),
                Max = sprintf("%.2f", max(value, na.rm = TRUE)),
                .groups = "drop"
            ) %>% as.data.frame()
    # Calculate p-value using t-test
    test <- t.test(value ~ group, data = df_person)
    p <- sprintf("%.3f", test$p.value)
    title <- glue("Summary of {col} by {grp_col}")
    cat("\n", title, ":\n", sep = "")
    print(summary, row.names = FALSE)
    cat("p-value (t-test):", ifelse(p < 0.001, "<0.001", p), "\n")
    invisible(summary)
}

# ------------------------------------------------------------------------
# SECTION 3: Functions for data visualization
# ------------------------------------------------------------------------

# Plot bar chart of imputed values
plot_bar <- function(
        df,
        x, 
        x_title = NULL
    ) {
    # Plot bar plot by subgroup
    plot <- ggplot(df, aes(x = .data[[x]])) +
        geom_bar(fill = "#9BBB59") +
        labs(x = x_title, y = "Count") +
        scale_x_discrete(
            labels = c(
                "0" = "Strongly Disagree",
                "1" = "Disagree",
                "2" = "Agree",
                "3" = "Strongly Agree"
            )) +
        theme_bw(base_size = 6) +
        theme(
            axis.title.x = element_text(margin = margin(t = 4)),
            axis.text.x = element_text(size = 4.5),
        )
    # Save bar plot
    fig_name <- glue("{x}_imputedVal")
    ggsave(
        filename = here("figs", paste0(fig_name, "_barChart.png")),
        plot = plot, height = 6, width = 6, units = "cm", dpi = 600
    )
}

# Plot bar chart by missingness
plot_bar_byMissing <- function(
        df,
        x, 
        fill,
        x_title = NULL
    ) {
    # Plot bar plot by subgroup
    plot <- ggplot(
        df, aes(
            x = .data[[x]],
            fill = factor(
                .data[[fill]],
                levels = c(0, 1),
                labels = c("Available", "Imputed")
            )
        )
    ) +
        geom_bar(position = position_stack(reverse = TRUE)) +
        labs(x = x_title, y = "Count") +
        scale_x_discrete(
            labels = c(
                "0" = "Strongly Disagree",
                "1" = "Disagree",
                "2" = "Agree",
                "3" = "Strongly Agree"
            )) +
        scale_fill_manual(
            values = c(
                "Available" = "#9BBB59",
                "Imputed" = "#C0504D"
            )
        ) +
        theme_bw(base_size = 6) +
        theme(
            axis.title.x = element_text(margin = margin(t = 4)),
            axis.text.x = element_text(size = 4.5),
            legend.position = "top",
            legend.margin = margin(t = -2, b = -2),
            legend.title = element_blank(),
            legend.key.size = unit(0.2, "cm"),
            legend.text = element_text(size = 4)
        )
    # Save bar plot
    fig_name <- glue("{x}_byMissing")
    ggsave(
        filename = here("figs", paste0(fig_name, "_barChart.png")),
        plot = plot, height = 6, width = 6, units = "cm", dpi = 600
    )
}

# Prepare df for line plot
linePrep <- function(df) {
    # Isolate columns of interest
    cols <- c("X", "Estimate", "Std.err.")
    df <- subset(df, select = cols)
    # Drop int term
    df <- df[!grepl("Int", df$X), ]
    # Assign grp to coeff
    df$grp <- character(nrow(df))
    for (i in 1:nrow(df)) {
        if (df$X[i] %in% c("ascPdt", "ascOpt")) {
            df$grp[i] <- "ASC"
        } else if (df$X[i] %in% c("bPreg_mv", "bInf_mab")) {
            df$grp[i] <- "Recipient"
        } else if (grepl("^(bNoAe|bRedness|bPreterm)", df$X[i])) {
            df$grp[i] <- "Adverse Event"
        } else if (grepl("^bHosp", df$X[i])) {
            df$grp[i] <- "Hospitalization"
        } else if (grepl("^bCost", df$X[i])) {
            df$grp[i] <- "Cost"
        }
    }
    # Calculate 95% CI
    df$ci_lower <- df$Estimate - 1.96 * df$Std.err.
    df$ci_upper <- df$Estimate + 1.96 * df$Std.err.
    return(df)
}

# Plot & save line plot
linePlot <- function(df) {
    # Obtain df name
    df_name <- deparse(substitute(df))
    # Fix factor levels
    df$X <- factor(df$X, levels = unique(df$X))
    # Plot line plot
    plot <- ggplot(df, aes(x = X, y = Estimate, group = grp, color = grp)) +
        geom_line(linewidth = 0.5) +
        geom_point(size = 1) +
        geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
        facet_wrap(~grp, scales = "free", nrow = 1) +
        labs(x = "Coefficients", y = "Estimates", color = "Group") +
        theme_bw(base_size = 8) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "top"
        )
    # Save line plot
    cat("\nSaved", paste0(df_name, "_linePlot.png"), "\n")
    ggsave(
        filename = here("figs", paste0(df_name, "_linePlot.png")),
        plot = plot, height = 6, width = 14, units = "cm", dpi = 600
    )
}

# ------------------------------------------------------------------------
# SECTION 4: Functions for coefficient estimation
# ------------------------------------------------------------------------

# Format CCA results
format_cca_results <- function(fit) {
    # Extract AIC & BIC
    aic <- fit$AIC
    bic <- fit$BIC
    # Extract estimates & robust standard errors
    est <- fit$estimate
    rob_se <- fit$robse
    # Calculate t-ratios & 2-sided p-values
    t_rat <- est / rob_se
    p_val <- 2 * pnorm(-abs(t_rat))
    # Create a data frame to store the results
    results <- data.frame(
        Est = sprintf("%.5f", est),
        Rob_SE = sprintf("%.5f", rob_se),
        p_2sided = ifelse(p_val < 0.001, "<0.001", sprintf("%.3f", p_val)),
        row.names = names(est)
    )
    cat("\nEstimates from complete case analysis:\n")
    print(results)
    invisible(results)
    cat("\nAIC:", sprintf("%.0f", aic), "\n")
    cat("\nBIC:", sprintf("%.0f", bic))
}

# Pool coefficients using Rubin's rules
pool_coeff <- function(est_list, varcov_list, M) {
    # Bind estimates from all imputations
    Qhat <- do.call(rbind, est_list) %>% as.data.frame()
    # Pooled estimate
    Qbar <- colMeans(Qhat)
    # Within-imputation variance
    Ubar <- Reduce("+", varcov_list) / M
    stopifnot(
        identical(colnames(Qhat), colnames(Ubar)),
        identical(colnames(Qhat), rownames(Ubar))
    )
    # Between-imputation variance
    B <- cov(Qhat)
    stopifnot(
        identical(colnames(Qhat), colnames(B)),
        identical(colnames(Qhat), rownames(B))
    )
    # Total variance
    T <- Ubar + (1 + 1 / M) * B
    # Standard errors
    SE <- sqrt(diag(T))
    # Relative increase in variance due to nonresponse
    r <- (1 + 1 / M) * diag(B) / diag(Ubar)
    # Degrees of freedom
    nu <- (M - 1) * (1 + 1 / r)^2
    # F statistic
    F <- Qbar^2 / diag(T)
    # Calculate p-values
    p <- pf(F, df1 = 1, df2 = nu, lower.tail = FALSE)
    # Create a data frame to store the results
    results <- data.frame(
        Est = sprintf("%.5f", Qbar),
        Rob_SE = sprintf("%.5f", SE),
        p_2sided = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)),
        row.names = names(Qbar)
    )
    cat("\nPooled estimates using Rubin's rules:\n")
    print(results)
    invisible(results)
}

# Summarize range of AIC & BIC values across imputations
summ_aic_bic <- function(aic_list, bic_list) {
    aic_vals <- unlist(aic_list)
    bic_vals <- unlist(bic_list)
    aic_min <- sprintf("%.0f", min(aic_vals))
    aic_max <- sprintf("%.0f", max(aic_vals))
    bic_min <- sprintf("%.0f", min(bic_vals))
    bic_max <- sprintf("%.0f", max(bic_vals))
    cat("\nAIC values:\n")
    cat("Min:", aic_min, ", Max:", aic_max, "\n")
    cat("\nBIC values:\n")
    cat("Min:", bic_min, ", Max:", bic_max)
}
