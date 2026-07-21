# ===========================================================================
# Title: 2_dataDesc.r
# Date: 2026-07-20
# ===========================================================================

# Load packages
library(here)
library(glue)
library(dplyr)

# Source universal functions
source(here("utils.r"))

# Create output folder if it doesn't exist
dir.create(here("output"), showWarnings = FALSE, recursive = TRUE)

# Write output to text file
sink(here("output", "2_dataDesc_output.txt"))
cat("===========================================================================\n")
cat("Title: 2_dataDesc_output.txt\n")
cat("Date: 2026-07-20\n")
cat("===========================================================================\n")

# ===========================================================================
# SECTION 1: Load dataset
# ===========================================================================

df <- read.csv(here("data", "jpn_dataset.csv"))

cat("\nSummary of variables:\n")
print(summary(df))

# ===========================================================================
# SECTION 2: Drop unused columns
# ===========================================================================

cols_toDrop <- c(
    "chiiki", "marriage", "job", "family", "home",
    "license", "car"
)
df <- df %>% select(-all_of(cols_toDrop))
cat("\nDropped columns:", cols_toDrop, "\n")

# ===========================================================================
# SECTION 3: Explore respondent characteristics
# ===========================================================================

cat_cols <- c(
    "sibling", "numSib", "fluCovVac", "tdapVac", "rsvVac",
    "govTrust", "vacImpt", "vacSafe", "vacEff", "vacBlf", 
    "education", "income", "sibGrp", "eduGrp", "incGrp"
)
summ_cat(df, cat_cols)

query_age(df$age)

query_gestAge(df$gestAge)

query_time(df$time)

# ===========================================================================
# SECTION 4: Investigate respondents that failed test &/or straightlined
# ===========================================================================

summ_cat_bySub(df, "allSame", "passTest")

cont_cols_toIx <- c("age", "gestAge", "time")
cat_cols_toIx <- c("sibGrp", "eduGrp", "incGrp")

for (col in cont_cols_toIx) {
    summ_cont_bySub(df, col, "passTest")
}

for (col in cat_cols_toIx) {
    summ_cat_bySub(df, "passTest", col)
}

for (col in cont_cols_toIx) {
    summ_cont_bySub(df, col, "allSame")
}

for (col in cat_cols_toIx) {
    summ_cat_bySub(df, "allSame", col)
}

# ===========================================================================
# SECTION 5: Drop respondents that failed test &/or straightlined
# ===========================================================================

df_subset <- subset(
    df,
    passTest == 1 & 
    allSame == 0 & 
    gestAge > 4 & gestAge <= 42
)
num_id <- length(unique(df_subset$id))
cat("\nKeep respondents who passed test, did not straightline & have valid gestational age\n")
cat("Number of unique IDs in the subset:", num_id, "\n")

indicator_cols <- c("govTrust", "vacImpt", "vacSafe", "vacEff", "vacBlf")
cat("\nSpearman correlation between attitudinal indicators (complete case):\n")
corr_mat <- cor(
    df_subset[indicator_cols],
    method = "spearman",
    use = "complete.obs"
)
print(round(corr_mat, 3))

# Save dataframe subset
write.csv(df_subset, file = here("data", "jpn_dataset_subset.csv"), row.names = FALSE)
cat("\nFile successfully saved as jpn_dataset_subset.csv\n")
cat("Column names of jpn_dataset_subset.csv:\n")
print(colnames(df_subset))

# Switch off writing
sink()