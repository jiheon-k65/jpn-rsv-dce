# ===========================================================================
# Title: 1_dataClean.r
# Date: 2026-10-18
# ===========================================================================

# Load packages
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)

# Source universal functions
source(here("utils.r"))

# Create output folder if it doesn't exist
dir.create(here("output"), showWarnings = FALSE, recursive = TRUE)

# Write output to text file
sink(here("output", "1_dataClean_output.txt"))
cat("===========================================================================\n")
cat("Title: 1_dataClean_output.txt\n")
cat("Date: 2026-10-18\n")
cat("===========================================================================\n")

# ===========================================================================
# SECTION 1: Create dataframe for DCE design
# ===========================================================================

# Load DCE design as dataframe
raw_design <- read.csv(
    here("data", "DCE_design_relabel.csv"),
    row.names = 1  # Set 1st column as index
)

# Rename columns
colnames(raw_design) <- c("block", "task", "alt", "preg", "hosp", "adv", "cost")

# Separate "adv" into "adv" & "risk"
raw_design <- separate(
    raw_design,
    col = "adv",
    into = c("adv", "risk"),
    sep = " ",
    remove = FALSE
)

# Recode into 0-indexed categories (except hosp & cost)
maps <- list(
    risk = c("L" = 0, "H" = 1),
    preg = c("A1" = 1, "A2" = 0),
    adv = c("C1" = 2, "C2" = 1, "C3" = 0),
    hosp = c("B1" = 1.5, "B2" = 1.2, "B3" = 0.9, "B4" = 0.3),
    cost = c("D1" = 0, "D2" = 0.05, "D3" = 0.1, "D4" = 0.3, "D5" = 1)
) # Cost per 100,000 JPY

for (i in names(maps)) {
    raw_design[[i]] <- as.numeric(maps[[i]][raw_design[[i]]])
}

# Isolate columns to create dummy columns
adv_col <- raw_design$adv
hosp_map <- c(
    "1.5" = 0,
    "1.2" = 1,
    "0.9" = 2,
    "0.3" = 3
)
hosp_col <- hosp_map[as.character(raw_design$hosp)]
cost_map <- c(
    "0" = 0,
    "0.05" = 1,
    "0.1" = 2,
    "0.3" = 3,
    "1" = 4
)
cost_col <- cost_map[as.character(raw_design$cost)]

# Create dummy columns for categorical variables
adv_dummy_cols <- create_dummy_cols(adv_col)
hosp_dummy_cols <- create_dummy_cols(hosp_col)
cost_dummy_cols <- create_dummy_cols(cost_col)

# Rename dummy columns
colnames(adv_dummy_cols) <- c("redness", "preterm")
colnames(hosp_dummy_cols) <- c("hosp1.2", "hosp0.9", "hosp0.3")
colnames(cost_dummy_cols) <- c("cost0.05", "cost0.1", "cost0.3", "cost1")

# Combine dummy columns with original dataframe
raw_design <- cbind(
    raw_design, adv_dummy_cols, hosp_dummy_cols, cost_dummy_cols
)
raw_design$adv <- NULL # Redundant column

# Create dataframe for 1st stage
design_alt1 <- raw_design[raw_design$alt == 1, ]
design_alt2 <- raw_design[raw_design$alt == 2, ]

# Add prefix to attribute columns
cols_noLabel <- c("block", "task", "alt", "risk")
cols_label <- setdiff(colnames(design_alt1), cols_noLabel)

design_alt1 <- add_prefix(design_alt1, cols_label, "alt1_")
design_alt2 <- add_prefix(design_alt2, cols_label, "alt2_")

# Combine alternative 1 & 2 dataframes
design <- cbind(design_alt1, design_alt2[paste0("alt2_", cols_label)])
stopifnot(
    nrow(design) == nrow(design_alt1),
    ncol(design) == ncol(design_alt1) + length(cols_label)
)
design$alt <- NULL # Redundant column
rownames(design) <- NULL # Reset index

# Sort columns alphabetically
design <- design[, order(names(design))]

# Save dataframe for DCE design
write.csv(design, file = here("data", "jpn_design.csv"), row.names = FALSE)
cat("\nFile successfully saved as jpn_design.csv\n")
cat("Column names of jpn_design.csv:\n")
print(colnames(design))

# ===========================================================================
# SECTION 2: Merge DCE data with design
# ===========================================================================

# Load raw data
raw_data <- read.csv(here("data", "rawdata.csv"))

# Isolate DCE columns
dce_colnames <- character()
for (x in 1:10) {
    for (y in 1:2) {
        dce_colnames <- c(dce_colnames, paste0("Q", x, "S", y))
    }
}
raw_dce <- raw_data[, c("RI_ID", "SEG", dce_colnames)]

# Rename selected columns
raw_dce <- raw_dce %>% rename(
    id = RI_ID,
    block = SEG
)

# Recode "block" values to 1-20
raw_dce$block <- ifelse(raw_dce$block > 20, raw_dce$block - 20, raw_dce$block)

# Indicate respondents that passed test & those with same responses
raw_dce <- indicate_passTest(raw_dce) %>% indicate_allSame()

# Separate "QxSy" into "task" & "stage"
raw_dce <- pivot_longer(
    raw_dce,
    cols = starts_with("Q"),
    names_to = "task",
    values_to = "choice"
)
raw_dce$task <- sub("^Q", "", raw_dce$task)
raw_dce <- separate(
    raw_dce,
    col = "task",
    into = c("task", "stage"),
    sep = "S",
    convert = TRUE
)

# Create separate dataframes for 1st & 2nd stages
dce_stage1 <- raw_dce[raw_dce$stage == 1, ]
dce_stage2 <- raw_dce[raw_dce$stage == 2, ]

# Combine 1st & 2nd stage dataframes
cols_toMerge <- c("id", "block", "task", "choice")
dce <- left_join(
    dce_stage1,
    dce_stage2[cols_toMerge],
    by = c("id", "block", "task"),
    suffix = c("", "_stage2")
)
dce$stage <- NULL # Redundant column
colnames(dce)[colnames(dce) == "choice"] <- "choice_stage1"
stopifnot(
    dim(dce)[1] == dim(dce_stage1)[1],
    dim(dce)[2] == dim(dce_stage1)[2]
)

# Merge DCE data with design
dceDesign <- left_join(
    dce,
    design,
    by = c("block", "task")
)
stopifnot(
    dim(dceDesign)[1] == dim(dce)[1],
    dim(dceDesign)[2] == dim(dce)[2] + ncol(design) - 2
)

# Indicate repeated task
dceDesign <- indicate_rptTask(dceDesign)
stopifnot(sum(dceDesign$rptTask) == length(unique(dceDesign$id)))

# Reassign columns as MV (left) & mAb (right)
# Indicate product availability in stage 2
dceDesign <- reassign_columns(dceDesign) %>% indicate_prodAv()

# Proportion of opt-out in stage 2
dceDesign <- dceDesign %>%
    group_by(id) %>%
    mutate(propOpt = mean(choice_stage2 == 3)) %>%
    ungroup()
cat("\nAdded column for proportion of opt-out (propOpt)\n")

# Save dataframe for DCE design & results
write.csv(dceDesign, file = here("data", "jpn_dceDesign.csv"), row.names = FALSE)
cat("\nFile successfully saved as jpn_dceDesign.csv\n")
cat("Column names of jpn_dceDesign.csv:\n")
print(colnames(dceDesign))

# ===========================================================================
# SECTION 3: Explore respondent characteristics
# ===========================================================================

# Remove DCE data & summarize missing data
char_data <- raw_data %>% select(-all_of(dce_colnames))

# Recode missing values
cols_toIndicate_pt1 <- paste0("Q", 16:19)
cols_toIndicate_pt2 <- c("Q12", "Q13")
char_data <- indicate_NA(char_data, c("Q11", cols_toIndicate_pt1), 5)
char_data <- indicate_NA(char_data, cols_toIndicate_pt2, 3)

# Calculate completion time
char_data$time <- as.numeric(
    difftime(
        as.POSIXct(char_data$end, format = "%Y/%m/%d %H:%M:%S"),
        as.POSIXct(char_data$start, format = "%Y/%m/%d %H:%M:%S"),
        units = "mins"
    )
)
cat("\nAdded column for completion time (time)\n")

# Drop redundant columns
cols_toDrop <- c(
    "start", "end", "X1", "SC1", "SC2",
    "SEG", "SEX", "AGE", "CHILD", "NAGE"
)
char_data <- char_data %>% select(-all_of(cols_toDrop))
cat("\nDropped redundant columns:", cols_toDrop, "\n")

# Rename columns
new_colnames <- c(
    "id", "age", "sibling", "numSib", "fluCovVac",
    "tdapVac", "rsvVac", "gestAge", "govTrust", "vacImpt",
    "vacSafe", "vacEff", "vacBlf", "chiiki", "marriage",
    "job", "family", "education", "income", "home",
    "license", "car", "time"
)
colnames(char_data) <- new_colnames

# ===========================================================================
# SECTION 4: Merge DCE design with respondent characteristics
# ===========================================================================

# Merge DCE design with respondent characteristics
dataset <- left_join(
    dceDesign,
    char_data,
    by = "id"
)
stopifnot(
    dim(dataset)[1] == dim(dceDesign)[1],
    dim(dataset)[2] == dim(dceDesign)[2] + ncol(char_data) - 1
)

# ===========================================================================
# SECTION 5: Create new columns for respondent characteristics
# ===========================================================================

# Recode binary categories to 0-indexed values
binary_cols <- c("sibling", "tdapVac", "rsvVac")
for (col in binary_cols) {
    dataset[[col]] <- ifelse(dataset[[col]] == 2, 0, 1)
} # No: 0, Yes: 1

# Recode flu &/or covid vaccine into 0-indexed values
fluCovVac_map <- c(
    "1" = 0, # None
    "2" = 1, # Flu only
    "3" = 2, # Covid-19 only
    "4" = 3 # Both
)
dataset$fluCovVac <- fluCovVac_map[as.character(dataset$fluCovVac)]

# Recode Likert scale to 0-indexed values
scale_cols <- c("govTrust", "vacImpt", "vacSafe", "vacEff", "vacBlf")
scale_map <- c(
    "1" = 3, # Strongly agree
    "2" = 2, # Somewhat agree
    "3" = 1, # Somewhat disagree
    "4" = 0 # Strongly disagree
)
for (col in scale_cols) {
    dataset[[col]] <- scale_map[as.character(dataset[[col]])]
}

# Group education into 2 categories
edu_map <- c(
    "1" = 0, # Junior high school
    "2" = 0, # High school
    "3" = 0, # Vocational school
    "4" = 0, # Junior college / college of technology
    "5" = 1, # Bachelor's degree
    "6" = 1, # Postgraduate degree
    "7" = 0, # Currently studying
    "8" = 0 # Others
)
dataset$eduGrp <- edu_map[as.character(dataset$education)]
cat("\nAdded column for education group (eduGrp)\n")

# Group annual income into 2 categories
inc_map <- c(
    "1" = 0, # <=4 million
    "2" = 0, # >4 to <=6 million
    "3" = 1, # >6 to <=8 million
    "4" = 1, # >8 to <=10 million
    "5" = 1, # >10 to <=12 million
    "6" = 1, # >12 to <=15 million
    "7" = 1 # >15 million
)
dataset$incGrp <- inc_map[as.character(dataset$income)]
cat("\nAdded column for income group (incGrp)\n")

# Save dataframe for DCE design, results & respondent characteristics
dataset <- dataset[, sort(names(dataset))]
write.csv(dataset, file = here("data", "jpn_dataset.csv"), row.names = FALSE)
cat("\nFile successfully saved as jpn_dataset.csv\n")
cat("Column names of jpn_dataset.csv:\n")
print(colnames(dataset))

# Switch off writing
sink()