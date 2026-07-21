# ===========================================================================
# Title: 3_dataImp.r
# Date: 2026-07-20
# ===========================================================================

# Load packages
library(here)
library(glue)
library(dplyr)
library(ggplot2)
library(mice, warn.conflicts = FALSE)

# Source universal functions
source(here("utils.r"))

# Create output folder if it doesn't exist
dir.create(here("output"), showWarnings = FALSE, recursive = TRUE)

# Write output to text file
sink(here("output", "3_dataImp_output.txt"))
cat("===========================================================================\n")
cat("Title: 3_dataImp_output.txt\n")
cat("Date: 2026-07-20\n")
cat("===========================================================================\n")

# ===========================================================================
# SECTION 1: Load dataset subset
# ===========================================================================

df <- read.csv(here("data", "jpn_dataset_subset.csv"))

cat("\nColumn names of jpn_dataset_subset.csv:\n")
print(colnames(df))

# Further drop unused columns
cols_unused <- c("fluCovVac", "tdapVac", "rsvVac")
df <- df %>% select(-all_of(cols_unused))
cat("\nDropped columns:", cols_unused, "\n")

# ===========================================================================
# SECTION 2: Save dataset for complete case analysis
# ===========================================================================

df_complete <- df[complete.cases(df), ]
num_id_missing <- length(unique(df$id)) - length(unique(df_complete$id))
prop_missing <- num_id_missing / length(unique(df$id)) * 100
cat("\nProportion of unique IDs with missing values:", round(prop_missing, 2), "%\n")

# Calculate mean of standardized indicators
cols_indicator <- c("vacImpt", "vacSafe", "vacEff", "vacBlf")
df_complete$vacConf <- rowMeans(df_complete[cols_indicator])

# Save dataframe for complete case analysis
cat("\nFile successfully saved as jpn_dataset_subset_complete.csv\n")
write.csv(df_complete, file = here("data", "jpn_dataset_subset_complete.csv"), row.names = FALSE)
cat("Column names of jpn_dataset_subset_complete.csv:\n")
print(colnames(df_complete))

# ===========================================================================
# SECTION 3: Investigate missing data patterns
# ===========================================================================

# Isolate columns for imputation
# Removed derived columns to avoid multicollinearity
cols_forImp <- c(
    "id", "age", "numSib", "gestAge", "govTrust",
    "vacImpt", "vacSafe", "vacEff", "vacBlf", "education",
    "income", "time", "propOpt"
)
df_forImp <- df %>% select(all_of(cols_forImp))

# Keep first row of each ID for for imputation
df_forImp <- df_forImp[!duplicated(df_forImp$id), ]
stopifnot(nrow(df_forImp) == length(unique(df$id)))

# Visualize missing data pattern
md.pattern(df_forImp, rotate.names = TRUE)

# Create dummy column for missingness
for (col in cols_indicator) {
    indicator_col <- paste0(col, "_missing")
    df_forImp[[indicator_col]] <- ifelse(is.na(df_forImp[[col]]), 1, 0)
} # Available: 0, Missing: 1
stopifnot (ncol(df_forImp) == length(cols_forImp) + length(cols_indicator))

cols_cat <- c("numSib", "govTrust", "education", "income")
cols_cont <- c("age", "gestAge", "time", "propOpt")

# Summarize vaccine importance missingness
for (col in cols_cat) {
    summ_cat_bySub(df_forImp, "vacImpt_missing", col)
}
for (col in cols_cont) {
    summ_cont_bySub(df_forImp, col, "vacImpt_missing")
}

# Summarize vaccine safety missingness
for (col in cols_cat) {
    summ_cat_bySub(df_forImp, "vacSafe_missing", col)
}
for (col in cols_cont) {
    summ_cont_bySub(df_forImp, col, "vacSafe_missing")
}

# Summarize vaccine effectiveness missingness
for (col in cols_cat) {
    summ_cat_bySub(df_forImp, "vacEff_missing", col)
}
for (col in cols_cont) {
    summ_cont_bySub(df_forImp, col, "vacEff_missing")
}

# Summarize vaccine belief missingness
for (col in cols_cat) {
    summ_cat_bySub(df_forImp, "vacBlf_missing", col)
}
for (col in cols_cont) {
    summ_cont_bySub(df_forImp, col, "vacBlf_missing")
}

# Fit multivariable logistic regression model for vaccine importance missingness
glm_vacImpt <- summary(glm(vacImpt_missing ~
      age +
      gestAge +
      numSib +
      govTrust +
      education +
      income +
      time +
      propOpt,
    family = binomial,
    data = df_forImp))
print(glm_vacImpt)

# Fit multivariable logistic regression model for vaccine safety missingness
glm_vacSafe <- summary(glm(vacSafe_missing ~
      age +
      gestAge +
      numSib +
      govTrust +
      education +
      income +
      time +
      propOpt,
    family = binomial,
    data = df_forImp))
print(glm_vacSafe)

# Fit multivariable logistic regression model for vaccine effectiveness missingness
glm_vacEff <- summary(glm(vacEff_missing ~
      age +
      gestAge +
      numSib +
      govTrust +
      education +
      income +
      time +
      propOpt,
    family = binomial,
    data = df_forImp))
print(glm_vacEff)

# Fit multivariable logistic regression model for vaccine belief missingness
glm_vacBlf <- summary(glm(vacBlf_missing ~
      age +
      gestAge +
      numSib +
      govTrust +
      education +
      income +
      time +
      propOpt,
    family = binomial,
    data = df_forImp))
print(glm_vacBlf)

# ===========================================================================
# SECTION 4: Conduct multiple imputation
# ===========================================================================

# Convert categorical variables to ordered factors
df_forImp[cols_indicator] <- lapply(
  df_forImp[cols_indicator],
  function(x) ordered(x, levels = 0:3)
)

# Omit selected columns from predictor matrix
pred <- make.predictorMatrix(df_forImp)
cols_toOmit <- c("id", paste0(cols_indicator, "_missing"))
pred[, cols_toOmit] <- 0

imp <- mice(
  df_forImp,
  method = "polr", # Ordered data
  predictorMatrix = pred,
  m = 25, # No. of imputations
  maxit = 40, # No. of iterations
  seed = 123,
  printFlag = FALSE
)

# Prepare original dataset for merging
df_toMerge <- df %>% select(-all_of(cols_indicator))

# Initialize list to store imputed datasets
df_imp_list <- vector("list", imp$m)

# Merge imputed values with original dataset & save CSV
for (i in 1:imp$m) {
    # Keep columns of interest for each imputed dataset
    cols_toKeep <- c("id", cols_indicator, paste0(cols_indicator, "_missing"))
    imp_i <- complete(imp, i) %>% select(all_of(cols_toKeep))
    # Merge imputed values with original dataset
    df_imp_i <- left_join(df_toMerge, imp_i, by = "id")
    stopifnot(
        nrow(df_imp_i) == nrow(df_toMerge),
        ncol(df_imp_i) == ncol(df_toMerge) + 2 * length(cols_indicator)
    )
    # Create new column for imputation number
    df_imp_i$imp <- i
    # Convert indicator columns to numeric
    df_imp_i[cols_indicator] <- lapply(
        df_imp_i[cols_indicator], function(x) as.numeric(as.character(x))
    )
    # Calculate mean of indicators
    df_imp_i$vacConf <- rowMeans(df_imp_i[cols_indicator])
    # Save dataset in list
    df_imp_list[[i]] <- df_imp_i
}

# Stack all imputations
df_imputed <- do.call(rbind, df_imp_list)
stopifnot(
    nrow(df_imputed) == nrow(df_toMerge) * imp$m,
    ncol(df_imputed) == ncol(df_toMerge) + 2 * length(cols_indicator) + 2
)

# Save imputed dataset
write.csv(df_imputed, file = here("data", "jpn_dataset_subset_imputed.csv"), row.names = FALSE)
cat("\nFile successfully saved as jpn_dataset_subset_imputed.csv\n")
cat("Column names of jpn_dataset_subset_imputed.csv:\n")
print(colnames(df_imputed))

# Save imputed values
for (col in cols_indicator) {
    imp_val <- imp$imp[[col]]
    title <- glue("{col}_imputedVal.csv")
    write.csv(imp_val, file = here("data", title), row.names = FALSE)
}

# ===========================================================================
# SECTION 5: Visualize imputed values
# ===========================================================================

# Create a new dataframe with one row per unique ID
df_person <- df_imputed[!duplicated(df_imputed$id), ]
df_person_imp1 <- df_person[df_person$imp == 1, ]
stopifnot(nrow(df_person_imp1) == length(unique(df_imputed$id)))

x_title_map <- list(
    vacImpt = "Vaccine Importance",
    vacSafe = "Vaccine Safety",
    vacEff = "Vaccine Effectiveness",
    vacBlf = "Vaccine Compatibility with Beliefs"
)

# Plot & save bar chart of imputed values for each indicator
for (col in cols_indicator) {
    # Flatten imputed values for plotting
    imp_val <- unlist(imp$imp[[col]], use.names = FALSE)
    # Convert to dataframe for plotting
    df_imp_val <- data.frame(imp_val)
    names(df_imp_val) <- col
    # Plot & save bar chart of imputed values
    plot_bar(
        df = df_imp_val,
        x = paste0(col),
        x_title = paste0(x_title_map[[col]], " (Imputed)")
    )
}

# Plot & save bar chart of missingness for each indicator
for (col in cols_indicator) {
    plot_bar_byMissing(
        df = complete(imp, 1), # Use first imputed dataset
        x = col,
        fill = paste0(col, "_missing"),
        x_title = paste0(x_title_map[[col]])
    )
}

# Switch off writing
sink()