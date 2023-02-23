library(tidyverse)
library(tidymodels)
library(vip)

compute_roc_curve <- function(data, class_feature, relevant_class, recipe_formula, training_proportion, seed) {
  # Data splitting
  splits      <- initial_split(data, strata = class_feature, prop = training_proportion)
  train_data <- training(splits)
  test_data  <- testing(splits)
  
  # Set up recipe
  rec <- recipe(recipe_formula, data = train_data) %>%
    step_center(all_predictors()) %>%
    step_scale(all_predictors())
  
  
  #TODO: tune hyperparameters and select best params
  
  rf_spec <-  rand_forest(min_n = 35, trees = 500) %>% 
    # Set up random forest model
    set_engine("ranger", num.threads = parallel::detectCores()) %>% 
    set_mode("classification")
  
  # Set up workflow with random forest model
  rf_wf <- workflow() %>%
    add_recipe(rec) %>%
    add_model(rf_spec)
  
  # Train random forest model
  rf_fit <- rf_wf %>% 
    fit(train_data)
  
  # Predict
  rf_preds <- predict(rf_fit, test_data) %>% 
    bind_cols(predict(rf_fit, test_data, type = "prob")) %>% 
    bind_cols(
      test_data %>% 
        select(class_feature)
    )

  # Calculate ROC curve data
  roc_curve_data <- rf_preds %>% 
    roc_curve(".pred_class", paste(".pred", relevant_class, sep="_"))

  return(roc_curve_data)
}

multivariate_roc_curve <- function(data, n_features_list, class_feature, relevant_class, training_proportion, seed) {
  # Create list to loop over
  total_features <- nrow(feature_importances)
  results <- vector("list", length = length(n_features_list))
  
  # Compute roc curve data for each subset of features
  for (i in 1:length(n_features_list)){
    n_features = n_features_list[[i]]
    selected_features_joined <- feature_importances %>% 
      top_n(n_features) %>% 
      .$Variable %>% 
      paste("`", ., "`", sep="", collapse="+")
    recipe_formula <- paste('`', class_feature, "` ~ ", selected_features_joined, sep="")
    print(paste("Computing ROC curve for the following formula: ", recipe_formula))
    roc_curve_data <- compute_roc_curve(data, class_feature, relevant_class, formula(recipe_formula), training_proportion)
    roc_curve_data <- roc_curve_data %>% 
      mutate(Model=paste("Top ", n_features, "features"))
    results[[i]] <- roc_curve_data
  }
  return(results)
}