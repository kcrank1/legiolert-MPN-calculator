# Legiolert MPN Calculator

A  R Shiny app for calculating **Legiolert MPN per 100 mL** and approximate 95% confidence limits from Quanti-Tray/Legiolert results. This is intended as an alternative to the IDEXX program for calculating Legiolert 95% Confidence Intervals.

## Inputs

- **Positive small wells:** 0–90
- **Positive large wells:** 0–6
- **Dilution factor:** use `1` for an undiluted sample; `0.1` for a 1:10 dilution, etc.

The app uses the fixed Legiolert tray parameters:

- 90 small wells (0.198 mL)
- 6 large wells (13.7 mL)

Well volumes were back calculated using the known total volume (100 mL) and the existing MPN table, and validated in the lab.

## Statistical method

The implementation follows the MPN equations described by:

Jarvis B, Wilrich C, Wilrich P-T. (2010). Reconsideration of the derivation of Most Probable Numbers, their standard deviations, confidence bounds and rarity values. *Journal of Applied Microbiology*, 109(5), 1660–1667.

DOI: https://doi.org/10.1111/j.1365-2672.2010.04792.x

The authors supplied a handy excel sheet that can be used to do these same calculations (if you input the derived well volumes above), however, on some computers the excel file doesn't work, so this app can be an alternative for Legiolert cases.

The 95% interval is the approximate interval based on ±2 standard deviations of ln(MPN), transformed back to the original scale.

## Shiny app

This app is available at https://qmraswim.shinyapps.io/legiolert-MPN/

## Run locally

Install R and RStudio, then install the packages:

```r
install.packages(c("shiny", "bslib", "DT"))
```

Put `legiolert-MPN-calc.R` in a folder and run:

```r
shiny::runApp()
```

## Notes

The app treats a tray with zero positive wells as `<1 MPN/100 mL` . Currently, no adjustment to that value is made for dilutions. 

Also, a screenshot of the equations used for calculations is included in the www folder, and it is a screenshot from the Jarvis (2010) excel file. 
