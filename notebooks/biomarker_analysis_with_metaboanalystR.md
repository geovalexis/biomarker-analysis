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
  - <a href="#14-multivariate-roc-curve-explorer"
    id="toc-14-multivariate-roc-curve-explorer">1.4 Multivariate ROC Curve
    Explorer</a>
  - <a href="#15-roc-curve-based-model-creation-and-evaluation"
    id="toc-15-roc-curve-based-model-creation-and-evaluation">1.5 ROC Curve
    Based Model Creation and Evaluation</a>
  - <a href="#16-sweave-report" id="toc-16-sweave-report">1.6 Sweave
    Report</a>

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

## 1.4 Multivariate ROC Curve Explorer

``` r
# Set the biomarker analysis mode to perform Multivariate exploratory ROC curve analysis ("explore")
mSet<-SetAnalysisMode(mSet, "explore")

# Prepare data for biomarker analysis
mSet<-PrepareROCData(mSet)

# Perform multivariate ROC curve analysis, using SVM classification and ranking
mSet<-PerformCV.explore(mSet, cls.method = "svm", rank.method = "svm", lvNum = 2)

### OPTION 1 Comparison plot of ROC curves of all models ###
mSet<-PlotROC(mSet, imgName = "ROC_all_models", format = "png", dpi = 300, mdl.inx= 0, avg.method = "threshold", show.conf = 0, show.holdout = 0, focus="fpr", cutoff=0.5)

# Plot predicted class probabilities for each sample for a selected model, not showing labels of wrongly classified samples
mSet<-PlotProbView(mSet, imgName = "multi_roc_prob", format = "png", dpi = 300, mdl.inx = -1, show = 0, showPred = 0)
```

    ## [1] "1773 duplicates are merged to their average"

``` r
# Plot the predictive accuracy of models with increasing number of features
mSet<-PlotAccuracy(mSet, imgName = "multi_roc_accuracy", format = "png", dpi = 300)

# Plot the most important features of a selected model ranked from most to least important
mSet<-PlotImpVars(mSet, imgName = "multi_roc_impvar", format="png", dpi=300, mdl.inx = -1, measure="freq", feat.num=15)
```

The results of executing the above code are the following:

- `ROC_all_modelsdpi300.png`: ROC curve plot
- `multi_roc_probdpi300.png`: Probabilities plot
- `multi_roc_accuracydpi300.png`: Predictive accuracies plot
- `multi_roc_impvardpi300.png`: ranking of important features plot
- `roc_pred_prob.csv`: numeric data for ROC curve
- `imp_features_cv.csv`: feature importance numeric data

Here we can see the resulting images:

<div class="figure" style="text-align: center">

<img src="ROC_all_modelsdpi300.png" alt="ROC_all_modelsdpi300.png" width="2400" />

<p class="caption">

ROC_all_modelsdpi300.png

</p>

</div>

<div class="figure" style="text-align: center">

<img src="multi_roc_probdpi300.png" alt="multi_roc_probdpi300.png" width="2700" />

<p class="caption">

multi_roc_probdpi300.png

</p>

</div>

<div class="figure" style="text-align: center">

<img src="multi_roc_accuracydpi300.png" alt="multi_roc_accuracydpi300.png" width="2700" />

<p class="caption">

multi_roc_accuracydpi300.png

</p>

</div>

<div class="figure" style="text-align: center">

<img src="multi_roc_impvardpi300.png" alt="multi_roc_impvardpi300.png" width="2400" />

<p class="caption">

multi_roc_impvardpi300.png

</p>

</div>

## 1.5 ROC Curve Based Model Creation and Evaluation

``` r
# Set the biomarker analysis mode to perform ROC Curve Based Model Creation and Evaluation ("test")
mSet<-SetAnalysisMode(mSet, "test")

# Prepare data for biomarker analysis
mSet<-PrepareROCData(mSet)

# Perform calculation of feature importance (AUC, p value, fold change)
mSet<-CalculateFeatureRanking(mSet)

# Manually select a subset of features for ROC analysis to build a classifier 
selected.cmpds <- c("Betaine", "N,N-Dimethylglycine", "Quinolinate", "Glucose")

# Manually select a subset of samples for ROC analysis hold-out data for validation purposes
selected.smpls <- c("PIF_178", "PIF_087", "PIF_090", "PIF_102", "PIF_111", "PIF_112")

# Prepare the custom data for model creation and sample hold-out 
mSet<-SetCustomData(mSet, selected.cmpds, selected.smpls)

# Perform ROC curve analysis, using SVM classification 
mSet<-PerformCV.test(mSet, method = "svm", lvNum = 2)

# Plot the ROC curve for the created model
mSet<-PlotROC(mSet, imgName = "cls_roc_0_", format="png",  dpi=300, mdl.inx = 0, avg.method = "threshold", 0, 0, "fpr", 0.5)

# Plot the predicted class probabilities for each sample using the user-created classifier, not showing labels of wrongly classified samples
mSet<-PlotProbView(mSet, imgName = "cls_prob_0_", format="png",  dpi=300, mdl.inx =-1, show=0, showPred= 0)
```

    ## [1] "3429 duplicates are merged to their average"

    ## Warning in write.table(test.df, file = "roc_pred_prob1.csv", sep = ",", :
    ## appending column names to file

``` r
# Plot the predictive accuracy of the model with increasing number of features
mSet<-PlotTestAccuracy(mSet, imgName = "cls_accu_0_", format="png",  dpi=300)

# Perform permutations tests using the area under the ROC curve as a measure of performance
mSet<-Perform.Permut(mSet, perf.measure = "auroc", perm.num = 500, propTraining = 2/3)

# Plot the results of the permutation tests
mSet<-Plot.Permutation(mSet, imgName = "roc_perm_1_", format="png",  dpi=300)

# View predicted classes of new samples (only applicable if samples with empty class labels were in the uploaded dataset)
mSet <- ROCPredSamplesTable(mSet) # Create table
```

The results files of the above are the following:

- `cls_roc_0_dpi300.png`: ROC curve for the selected model

<img src="cls_roc_0_dpi300.png" width="2400" style="display: block; margin: auto;" />

- `cls_prob_0_dpi300.png`: scatter plot of predicted class probabilities

<img src="cls_roc_0_dpi300.png" width="2400" style="display: block; margin: auto;" />

- `cls_accu_0_dpi300.png`: Box plot of the predictive accuracy

<img src="cls_prob_0_dpi300.png" width="2700" style="display: block; margin: auto;" />

- `roc_perm_1_dpi300.png`: Plot of the permutations tests using the area
  under the ROC curve or the predictive accuracy

<img src="roc_perm_1_dpi300.png" width="2400" style="display: block; margin: auto;" />

To view the example of the results:

``` r
# View the the predictive accuracy results of the model
GetAccuracyInfo(mSet)
```

    ## [1] "The average accuracy based on 100 cross validations is 0.705. The accuracy for hold out data prediction is 0.667(4/6)."

## 1.6 Sweave Report

``` r
PreparePDFReport(mSet, "Human Cachexia")
```

    ## Writing to file Analysis_Report.tex
    ## Processing code chunks with options ...
    ##  1 : keep.source term verbatim (Analysis_Report.Rnw:6)
    ##  2 : keep.source term tex (Analysis_Report.Rnw:104)
    ##  3 : keep.source term tex (Analysis_Report.Rnw:209)
    ##  4 : keep.source term verbatim (Analysis_Report.Rnw:214)
    ## [1] "R version 4.1.1 (2021-08-10)"
    ## 
    ## You can now run (pdf)latex on 'Analysis_Report.tex'

    ## Warning in system(paste(shQuote(texi2dvi), "--version"),
    ## intern = TRUE): running command
    ## '"C:\Users\geova\AppData\Local\Programs\MiKTeX\miktex\bin\x64\texify.exe"
    ## --version' had status 5

    ## [1] 1
