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
- `4_mnl.r` onwards - fit choice models
- `18_calc_wtp.r` onwards - simulate predictions
- `24_plot_choiceProb.r` onwards - visualize predictions

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
|`7_mxl_costCont.r`|All main effects as random coefficients|

### MNL models with interactions
|Script|Description|
|----|-----------|
|`8.x_mnl_costCont_ascInt.r`|ASC interactions|
|`9.x_mnl_costCont_pregInt.r`|Pregnant women interactions|
|`10.x_mnl_costCont_pretermInt.r`|Preterm birth interactions|
|`11.x_mnl_costCont_costInt.r`|Cost interactions|
|`12.x_mnl_costCont_hosp1.2Int.r`|12/1,000 risk (60% efficacy) interactions|
|`13.x_mnl_costCont_hosp0.9Int.r`|9/1,000 risk (70% efficacy) interactions|
|`14.x_mnl_costCont_hosp0.3Int.r`|3/1,000 risk (90% efficacy) interactions|

### Hybrid choice models (HCM)
|Script|Description|
|----|-----------|
|`15.x_hcm_costCont.r`|All fixed coefficients|
|`16.x_hcm_costCont_ranInt.r`|Random intercept|
|`17.x_hcm_costCont_rantInt_ranCoeff.r`|Random intercept & coefficients|

**Note:** complete-case (`x=1`) & multiple imputation (`x=2`)

## 6. Model predictions

- `18_calc_wtp.r` - Willingness-to-pay (WTP)
- `19_calc_cv.r` - Compensating variation (CV)
- `20_calc_choiceProb.r` - Choice probabilities at base cost for both products
- `21_calc_choiceProb_mvCost.r` - Choice probabilities across varying MV costs
- `22_calc_choicProb_mabCost.r` - Choice probabilities across varying mAb costs
- `23_calc_choicProb_bothCost.r` - Choice probabilities across varying MV & mAb costs

## 7. Model visualization

- `24_plot_choiceProb.r` - Plot choice probabilities as individual line plots
- `25_plot_choiceProb_combined.r` - Plot choice probabilities as a combined line plot
- `26_plot_choiceProb_heatmat.r` - Plot choice probabilities as a heatmap

## 8. Outputs

- Derived datasets are written back to `data/`
- Figures are written to `figs/`
- Script logs are written to `output/`