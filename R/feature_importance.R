library(tidyverse)
library(tidymodels)
library(vip)

compute_features_importance <- function(data, class_feature, relevant_class, training_proportion, seed) {
  set.seed(seed)
  
  # Data splitting
  splits      <- initial_split(data, strata = class_feature, prop = training_proportion)
  train_data <- training(splits)
  test_data  <- testing(splits)
  
  # Set up recipe
  recipe_formula <- formula(paste('`', class_feature, "` ~ ", ".", sep=""))
  rec <- recipe(recipe_formula, data = train_data) %>%
    step_center(all_predictors()) %>%
    step_scale(all_predictors())
  
  #TODO: tune hyperparameters and select best params
  
  rf_spec <-  rand_forest(min_n = 35, trees = 500) %>% 
    # Set up random forest model
    set_engine("ranger", num.threads = parallel::detectCores(), importance = "impurity") %>% 
    set_mode("classification")
  
  # Set up workflow with random forest model
  rf_wf <- workflow() %>%
    add_recipe(rec) %>%
    add_model(rf_spec)
  
  # Train random forest model
  rf_fit <- rf_wf %>% 
    fit(train_data)
  
  feature_importances <- rf_fit %>% 
    extract_fit_parsnip() %>% 
    vi()
  
  return(feature_importances)
}

