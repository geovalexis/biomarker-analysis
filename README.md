# API for biomarker discovery

## Summary

This work intends to reproduce the biomarker analysis performed by [Metaboanalyst R package](https://github.com/xia-lab/MetaboAnalystR)
using the [Tidymodels](https://www.tidymodels.org/) framework and serving the as an API endpoint using [Plumber](https://www.rplumber.io/).

For a more detailed description of the analysis, please refer to the [original paper](https://doi.org/10.1007%2Fs11306-012-0482-9).

## Description

Metabolites can provide valuable insights into the underlying processes and pathways of the body, and
changes in their level. As such, one common goal of metabolomic studies is biomarker discovery, aiming to identify a metabolite or a set of metabolites capable of classifying conditions or disease, with high sensitivity (true-positive rate) and specificity (true negative rate).

The MetaboAnalystR package provides a set of functions to perform a ROC AUC based analysis
to identify metabolites that can be used as biomarkers for a given dataset. An example of a biomarker discovery
analysis can be found [here](https://rdrr.io/github/simscr/metaboanalyst/f/vignettes/Biomarker_Analysis.Rmd).

In this work, we have reproduced this same analysis with MetaboAnalystR using the Tidymodels framework. The notebooks with the analysis can be found in the `notebooks` folder or also the rendered version an be found:

- [Biomarker Discovery with MetaboAnalystR](https://geovalexis.github.io/biomarker-analysis/notebooks/01-biomarker-analysis-with-metaboanalyst.html)
- [Biomarker Discovery with Tidymodels](https://geovalexis.github.io/biomarker-analysis/notebooks/02-biomarker-analysis-with-tidymodels.html)

## Requirements

The following R packages are required to run the API endpoint:

- tidyverse (>= 1.3)
- tidymodels (>= 1.0)
- plumber (>= 1.2)

## Usage

To perform the analysis, a web server with specific endpoints is provided. It comes with a swagger
documentation that can be used to interact with the API.

If you use RStudio, you can use the `Run API` button to run the server.
Otherwise, you can run the server by using the following command:

```r
plumb(file='R/biomarker-discovery.api.R')$run()
```

The server will be available at `http://localhost:XXXX` and the swagger documentation will be available at `http://localhost:XXXX/__docs__/`.
