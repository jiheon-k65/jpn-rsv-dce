# jpn-rsv-dce

Analaysis codes for a discrete choice experiment (DCE) conducted among Japanese pregnant women to understand their preferences for respiratory syncytial virus (RSV) immunization products. This repository contains R scripts to format DCE dataset, conduct descriptive analyses, impute missing values & fit various choice models.

## 1. Folders

- `analysis/` - numbered R scripts for data preparation, multiple imputation & analysis
- `data/` - raw, cleaned & derived datasets to be used by analysis scripts
- `figs/` - exported figures from analysis scripts
- `output/` - text log files from analysis scripts
- `utils.r` - shared helper functions used across analysis scripts
- `renv/`, `renv.lock` - project library management

## 2. Requirements

- R 4.5.x
- `renv`

## 3. Setup

Restore the project library after opening the repository:

```r
renv::restore()
```

The scripts assume the project root is the working directory & use `here()` for script paths

## 4. Workflow

- `1_dataClean.r` - format dataset & re-code respondent characteristics
- `2_dataDesc.r` - describe data & respondents that failed quality checks
- `3_dataImp.r` - describe missingness & conduct multiple imputation
- `4_mnl.r` onwards - fit choice models to data
- `19_mxl_costCont_ranCoeff_choiceProb.r` onwards - calculate choice probabilities

## 5. Choice Models

### Multinomial logit (MNL) models
|Script|Description|Reported|
|----|-----------|:------:|
|`4_mnl.r`|All categorical variables|✅|
|`5_mnl_costCont.r`|Linear cost coefficient||
|`6_mnl_costContQuad.r`|Quadratic cost coefficient||

### Mixed logit (MXL) models
|Script|Description|Reported|
|----|-----------|:------:|
|`7_mxl_costCont.r`|All random coefficients||

### MNL models with interactions
|Script|Description|Reported|
|----|-----------|:------:|
|`8.x_mnl_costCont_ascInt.r`|ASC interactions|✅|
|`9.x_mnl_costCont_pregInt.r`|Pregnant women interactions|✅|
|`10.x_mnl_costCont_pretermInt.r`|Preterm birth interactions|✅|
|`11.x_mnl_costCont_costInt.r`|Cost interactions|✅|
|`12.x_mnl_costCont_fullInt.r`|Combined significant interactions (`8.2`-`11.2`)||
|`13.x_mnl_costCont_parInt.r`|Reduced model from `12.2`||
|`14.x_mnl_costCont_parInt_final.r`|Final reduced model from `13.2`||

### MXL models with interactions
|Script|Description|Reported|
|----|-----------|:------:|
|`15.x_mxl_costCont_parInt.r`|Selected random coefficient with interactions||
|`16.x_mxl_costCont_parInt_final.r`|Final reduced model from `15.2`||

### Hybrid choice models (HCM)
|Script|Description|Reported|
|----|-----------|:------:|
|`17.x_hcm_costCont.r`|All fixed coefficients||
|`18.x_hcm_costCont_ranCoeff.r`|Selected random coefficients|✅|

**Note:** complete-case (`x=1`) & multiple imputation (`x=2`)

## 6. Outputs

- Derived datasets are written back to `data/`
- Figures are written to `figs/`
- Script logs are written to `output/`