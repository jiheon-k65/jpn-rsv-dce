# jpn-rsv-dce

Analysis codes for a discrete choice experiment (DCE) conducted among Japanese pregnant women to understand their preferences for respiratory syncytial virus (RSV) immunization products. This repository contains R scripts to format DCE dataset, conduct descriptive analyses, impute missing values & fit various choice models.

## 1. Folders

- `analysis/` - numbered R scripts for data preparation, multiple imputation & analysis
- `data/` - raw, cleaned & derived datasets to be used by analysis scripts
- `figs/` - exported figures from analysis scripts
- `markdown/` - markdown files for equations & legend
- `output/` - text log files from analysis scripts
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
- `15_hcm_costCont_ranInt_ranCoeff_wtp.r` onwards - model predictions

## 5. Choice Models

### Multinomial logit (MNL) models
|Script|Description|
|----|-----------|
|`4_mnl.r`|All categorical variables|
|`5_mnl_costCont.r`|Linear cost coefficient|
|`6_mnl_costContQuad.r`|Quadratic cost coefficient|

### Mixed logit (MXL) models
|Script|Description|
|----|-----------|
|`7_mxl_costCont.r`|All random coefficients|

### MNL models with interactions
|Script|Description|
|----|-----------|
|`8.x_mnl_costCont_ascInt.r`|ASC interactions|
|`9.x_mnl_costCont_pregInt.r`|Pregnant women interactions|
|`10.x_mnl_costCont_pretermInt.r`|Preterm birth interactions|
|`11.x_mnl_costCont_costInt.r`|Cost interactions|

### Hybrid choice models (HCM)
|Script|Description|
|----|-----------|
|`12.x_hcm_costCont.r`|All fixed coefficients|
|`13.x_hcm_costCont_ranInt.r`|Random intercept|
|`14.x_hcm_costCont_rantInt_ranCoeff.r`|Random intercept & coefficients|

**Note:** complete-case (`x=1`) & multiple imputation (`x=2`)

## 6. Model predictions

- `15_hcm_costCont_ranInt_ranCoeff_wtp.r` - Willingness-to-pay (WTP)
- `16_hcm_costCont_ranInt_ranCoeff_cv.r` - Compensating variation (CV)
- `17_hcm_costCont_ranInt_ranCoeff_choiceProb.r` - Choice probabilities at base cost
- `18_hcm_costCont_ranInt_ranCoeff_choiceProb_mvCost.r` - Choice probabilities across varying MV costs
- `19_hcm_costCont_ranInt_ranCoeff_choiceProb_mabCost.r` - Choice probabilities across varying mAb costs
- `20_hcm_costCont_ranInt_ranCoeff_choiceProb_plot.r` - Plot choice probabilities separately
- `21_hcm_costCont_ranInt_ranCoeff_choiceProb_comPlot.r` - Plot choice probabilities together
- `22_hcm_costCont_ranInt_ranCoeff_choiceProb_bothCosts.r` - Choice probabilities across simultaneously varying MV & mAb costs
- `23_hcm_costCont_ranInt_ranCoeff_choiceProb_heatmap.r` - Plot choice probabilities with simultaneously varying MV & mAb costs as a heatmap

## 7. Outputs

- Derived datasets are written back to `data/`
- Figures are written to `figs/`
- Script logs are written to `output/`