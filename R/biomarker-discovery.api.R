library(plumber)
library(tidyverse)
library(tidymodels)

source("feature_importance.R")
source("roc_curve.R")

#* @apiTitle Biomarker discovery API
#* @apiDescription

#* Check if the API is up and running
#* @get /is-alive
function() {
  print("I am alive!")
}

#* Multivariate ROC curve
#* @param data_url url for the CSV
#* @serializer png
#* @get /biomarker-discovery/roc-curve
function(data_url) {
  class_feature <- "Muscle loss"
  relevant_class <- "cachexic"
  training_proportion <- 3/4
  
  # Read data
  data <- read_csv(data_url)
  data <- data %>% mutate(class_feature=factor(class_feature))
  
  # Perform analysis
  n_features_list <- c(3, 5, 10, round(1/5 * total_features), total_features)
  results <- multivariate_roc_curve(data, n_features_list, class_feature, relevant_class, training_proportion, seed=123)
  
  # Plot ROC curve for the different models
  plot <- bind_rows(results) %>% 
    ggplot(aes(x = 1 - specificity, y = sensitivity, col = Model)) + 
    geom_path(linewidth = 1.0, alpha = 0.8) +
    geom_abline(lty = 3) + 
    coord_equal() + 
    scale_color_viridis_d(option = "plasma", end = .6)
  print(plot)
  
}




# Programmatically alter your API
#* @plumber
function(pr) {
  pr %>%
    # Overwrite the default serializer to return unboxed JSON
    pr_set_serializer(serializer_unboxed_json())
}
