library(plumber)
library(tidyverse)
library(tidymodels)

source("feature_importance.R")
source("roc_curve.R")

SEED = 123

#* @apiTitle Biomarker Discovery API
#* @apiDescription Discovery biomarkers based on ROC AUC method

#* Check if the API is up and running
#* @get /is-alive
function() {
  print("I am alive!")
}

#* Ranking of features importance based on ROC AUC
#* @tag "Biomarker Discovery"
#* @param data_url url of the data
#* @param class_feature name of the target feature
#* @param relevant_class relevant class within the target feature
#* @param training_proportion proportion of the training set
#* @serializer unboxedJSON
#* @get /v1/biomarker-discovery/features-importance
function(data_url, class_feature, relevant_class, training_proportion="0.75") {
  # Treat query params
  training_proportion <- as.numeric(training_proportion)
  
  # Read data
  data <- read_csv(data_url)
  
  # Perform analysis
  features_importance <- compute_features_importance(
    data, 
    class_feature, 
    relevant_class, 
    training_proportion, 
    seed=SEED
  )
  return(feature_importances)
}

#* Multivariate ROC curve
#* @tag "Biomarker Discovery"
#* @param data_url url of the data
#* @param class_feature name of the target feature
#* @param relevant_class relevant class within the target feature
#* @param training_proportion proportion of the training set
#* @serializer png
#* @get /v1/biomarker-discovery/roc-curve
function(data_url, class_feature, relevant_class, training_proportion="0.75") {
  # Treat query params
  training_proportion <- as.numeric(training_proportion)
  
  # Read and treat data
  data <- read_csv(data_url)
  data[[class_feature]] <- factor(data[[class_feature]])

  # Perform analysis
  features_importance <- compute_features_importance(
    data, 
    class_feature, 
    relevant_class, 
    training_proportion, 
    seed=SEED
  )  
  n_features_list <- c(3, 5, 10, round(1/5 * total_features), total_features)
  results <- multivariate_roc_curve(
      data, 
      n_features_list, 
      class_feature, 
      relevant_class, 
      training_proportion, 
      features_importance, 
      seed=SEED
  )
  
  # Plot ROC curve for the different models
  plot <- bind_rows(results) %>% 
    ggplot(aes(x = 1 - specificity, y = sensitivity, col = Model)) + 
    geom_path(linewidth = 1.0, alpha = 0.8) +
    geom_abline(lty = 3) + 
    coord_equal() + 
    scale_color_viridis_d(option = "plasma", end = .6)
  print(plot)
  
}
