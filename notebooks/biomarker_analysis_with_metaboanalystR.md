BiomarkerAnalysis
================
Geovanny Risco
2023-02-13

- <a href="#1-biomarker-analysis-with-metaboanalystr"
  id="toc-1-biomarker-analysis-with-metaboanalystr">1 Biomarker Analysis
  with MetaboAnalystR</a>
  - <a href="#11-installation" id="toc-11-installation">1.1 Installation</a>
  - <a href="#12-importing-data" id="toc-12-importing-data">1.2 Importing
    data</a>
  - <a href="#13-classical-roc-curve-analysis"
    id="toc-13-classical-roc-curve-analysis">1.3 Classical ROC Curve
    Analysis</a>

# 1 Biomarker Analysis with MetaboAnalystR

## 1.1 Installation

``` r
# Install dependencies
metanr_packages <- function(){

  metr_pkgs <- c("impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", "sva", "limma", "KEGGgraph", "siggenes","BiocParallel", "MSnbase", "multtest","RBGL","edgeR","fgsea","devtools","crmn","httr","qs")
  
  list_installed <- installed.packages()
  
  new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
  
  if(length(new_pkgs)!=0){
    
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install(new_pkgs)
    print(c(new_pkgs, " packages added..."))
  }
  
  if((length(new_pkgs)<1)){
    print("No new packages added...")
  }
}
metanr_packages()

# Install MetaboAnalystR
install.packages("devtools")
library(devtools)
devtools::install_github("xia-lab/MetaboAnalystR", build = TRUE, build_vignettes = FALSE)
```

``` r
# Import library
library(MetaboAnalystR)
```

    ## MetaboAnalystR 3.3.0 initialized Successfully !
    ## https://github.com/xia-lab/MetaboAnalystR

## 1.2 Importing data

``` r
# Create objects for storing processed data from biomarker analysis
mSet<-InitDataObjects("conc", "roc", FALSE)
```

    ## Warning: package 'Rserve' was built under R version 4.1.3

    ## Starting Rserve...
    ##  "D:\DOCUME~1\R\WIN-LI~1\4.1\Rserve\libs\x64\Rserve.exe" --no-save

    ## Warning in Cairo::CairoFonts(regular = "Arial:style=Regular", bold =
    ## "Arial:style=Bold", : CairoFonts() has no effect on Windows. Please use
    ## par(family="...") to specify the desired font - see ?par.

    ## [1] "MetaboAnalyst R objects initialized ..."

``` r
# Read in data and fill in the dataSet list
mSet<-Read.TextData(mSet, "https://www.xialab.ca/api/download/metaboanalyst/human_cachexia.csv")
# Sanity check, replace missing values, check if the sample size is too small
mSet<-SanityCheckData(mSet)
```

    ##  [1] "Successfully passed sanity check!"                                                                                
    ##  [2] "Samples are not paired."                                                                                          
    ##  [3] "2 groups were detected in samples."                                                                               
    ##  [4] "Only English letters, numbers, underscore, hyphen and forward slash (/) are allowed."                             
    ##  [5] "<font color=\"orange\">Other special characters or punctuations (if any) will be stripped off.</font>"            
    ##  [6] "All data values are numeric."                                                                                     
    ##  [7] "A total of 0 (0%) missing values were detected."                                                                  
    ##  [8] "<u>By default, missing values will be replaced by 1/5 of min positive values of their corresponding variables</u>"
    ##  [9] "Click the <b>Proceed</b> button if you accept the default practice;"                                              
    ## [10] "Or click the <b>Missing Values</b> button to use other methods."

``` r
mSet<-ReplaceMin(mSet)
mSet<-IsSmallSmplSize(mSet)
```

    ## [1] 0

``` r
mSet<-PreparePrenormData(mSet)
```

``` r
###### *** OPTION 2 FOR NORMALIZATION
# No normalization, and computeS metabolite ratios and includeS the top 20 
mSet<-Normalization(mSet, "NULL", "NULL", "NULL", "C01", ratio=TRUE, ratioNum=20)

# If ratio = TRUE: view the normalized dataset including the top ranked ratios
# The ratios will be towards the end of the matrix 
mSet$dataSet$norm

#If ratio = TRUE: view just the top ranked included ratios
mSet$dataSet$ratio
```

## 1.3 Classical ROC Curve Analysis

``` r
# Set the biomarker analysis mode to perform Classical ROC curve analysis ("univ")
mSet<-SetAnalysisMode(mSet, "univ")
# Prepare data for biomarker analysis
mSet<-PrepareROCData(mSet)

### OPTION 1 Perform univariate ROC curve analysis ###
mSet<-Perform.UnivROC(mSet, feat.nm = "Creatine", version = "1", format="png", dpi=300, isAUC=F, isOpt=T, optMethod="closest.topleft", isPartial=F, measure="sp", cutoff=0.2)
```

    ## Setting levels: control = cachexic, case = control

    ## Setting direction: controls > cases

``` r
# Create box plot showing the concentrations of the selected compound between the groups
mSet<-PlotRocUnivBoxPlot(mSet, "Creatine", version= "1", "png", 72, T, FALSE)
# Perform calculation of feature importance (AUC, p value, fold change)
mSet<-CalculateFeatureRanking(mSet)
```

The results of executing the above code are the following:

- `metaboanalyst_roc_univ.csv` (numeric results for AUC, p-value, FC and
  clusters)
- `Creatine_box_1dpi72.png` (Boxplot for both groups)
- `Creatine_1dpi300.png` (ROC curve)

Here we can see the resulting images:

<div class="figure" style="text-align: center">

<img src="Creatine_box_1dpi72.png" alt="Creatine_box_1dpi72.png" width="200" />

<p class="caption">

Creatine_box_1dpi72.png

</p>

</div>

<div class="figure" style="text-align: center">

<img src="Creatine_1dpi300.png" alt="Creatine_1dpi300.png" width="1800" />

<p class="caption">

Creatine_1dpi300.png

</p>

</div>
