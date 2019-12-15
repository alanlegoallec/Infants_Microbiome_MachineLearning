




tune_caret_models <- function(analysis, target, predictors, algo, i, data_train, x_train, y_train, w_train, data_test, x_test, y_test)
{
  if(algo=="glmnet") {
    names_hyper <- c("lambda", "alpha")
    strides_hyper <- list("lambda"=1, "alpha"=1)
    n_extra_leftright <- 5
    n_fine_tuning <- list("lambda"=10, "alpha"=10)
    min_hyper <- list("lambda"=-Inf, "alpha"=-Inf)
    max_hyper <- list("lambda"=Inf, "alpha"=Inf)
    is_int_hyper <- list("lambda"=F, "alpha"=F)
    list_hyper <- list("lambda"=seq(-7,3,by=strides_hyper[["lambda"]]), "alpha"=c(-Inf, seq(-3,3,by = strides_hyper[["alpha"]]), Inf))
    hyper_fun <- list("lambda"=ten_power_x, "alpha"=inv.logit)
    hyper_fun_inv <- list("lambda"=log10, "alpha"=logit)
    transformed_hyper <- transform_hyper(list_hyper, hyper_fun)
  } else if (algo=="gbm") {
    names_hyper <- c("shrinkage", "interaction.depth", "n.minobsinnode")
    #n_trees <- 1000
    n_trees <- 10
    strides_hyper <- list("shrinkage"=1, "interaction.depth"=1, "n.minobsinnode"=4)
    n_extra_leftright <- 3
    n_fine_tuning <- list("shrinkage"=3, "interaction.depth"=2, "n.minobsinnode"=1)
    min_hyper <- list("shrinkage"=-Inf, "interaction.depth"=1, "n.minobsinnode"=1)
    max_hyper <- list("shrinkage"=Inf, "interaction.depth"=Inf, "n.minobsinnode"=floor(ncol(x_train)/3))
    is_int_hyper <- list("shrinkage"=F, "interaction.depth"=T, "n.minobsinnode"=T)
    list_hyper <- list("shrinkage"=seq(-2,0,by=strides_hyper[["shrinkage"]]), "interaction.depth"=c(1,2,3), "n.minobsinnode"=seq(5,10))
    hyper_fun <- list("shrinkage"=ten_power_x, "interaction.depth"=identity, "n.minobsinnode"=identity)
    hyper_fun_inv <- list("shrinkage"=log10, "interaction.depth"=identity, "n.minobsinnode"=identity)
    transformed_hyper <- transform_hyper(list_hyper, hyper_fun)
    transformed_hyper$n.trees <- n_trees
  } else if (algo=="rf") {
    names_hyper <- c("mtry")
    n_partition <- 5
    strides_hyper <- list("mtry"=ceiling(ncol(x_train)/n_partition))
    n_extra_leftright <- 1
    n_fine_tuning <- list("mtry"=4)
    min_hyper <- list("mtry"=1)
    max_hyper <- list("mtry"=ncol(x_train))
    is_int_hyper <- list("mtry"=T)
    list_hyper <- list("mtry"=seq(strides_hyper[["mtry"]],strides_hyper[["mtry"]]*(n_partition-1),by=strides_hyper[["mtry"]]))
    hyper_fun <- list("mtry"=identity)
    hyper_fun_inv <- list("mtry"=identity)
    transformed_hyper <- transform_hyper(list_hyper, hyper_fun)
  } else if (algo=="svmLinear") {
    names_hyper <- c("C")
    strides_hyper <- list("C"=1)
    n_extra_leftright <- 5
    n_fine_tuning <- list("C"=10)
    min_hyper <- list("C"=-Inf)
    max_hyper <- list("C"=Inf)
    is_int_hyper <- list("C"=F)
    list_hyper <- list("C"=seq(-7,3,by=strides_hyper[["C"]]))
    hyper_fun <- list("C"=ten_power_x)
    hyper_fun_inv <- list("C"=log10)
    transformed_hyper <- transform_hyper(list_hyper, hyper_fun)
  } else if (algo=="svmPoly") {
    names_hyper <- c("C", "scale")
    n_poly <- 2
    strides_hyper <- list("C"=1, "scale"=1)
    n_extra_leftright <- 3
    n_fine_tuning <- list("C"=3, "scale"=3)
    min_hyper <- list("C"=-Inf, "scale"=-Inf)
    max_hyper <- list("C"=Inf, "scale"=Inf)
    is_int_hyper <- list("C"=F, "scale"=F)
    list_hyper <- list("C"=seq(-2,0,by=strides_hyper[["C"]]), "scale"=c(-Inf, seq(-2,0,by = strides_hyper[["scale"]]), Inf))
    hyper_fun <- list("C"=ten_power_x, "scale"=inv.logit)
    hyper_fun_inv <- list("C"=log10, "scale"=logit)
    transformed_hyper <- transform_hyper(list_hyper, hyper_fun)
    transformed_hyper$degree <- n_poly
  } else if (algo=="svmRadial") {
    names_hyper <- c("C", "sigma")
    strides_hyper <- list("C"=1, "sigma"=1)
    n_extra_leftright <- 3
    n_fine_tuning <- list("C"=3, "sigma"=3)
    min_hyper <- list("C"=-Inf, "sigma"=-Inf)
    max_hyper <- list("C"=Inf, "sigma"=Inf)
    is_int_hyper <- list("C"=F, "sigma"=F)
    list_hyper <- list("C"=seq(-2,0,by=strides_hyper[["C"]]), "sigma"=c(-Inf, seq(-2,0,by = strides_hyper[["sigma"]]), Inf))
    hyper_fun <- list("C"=ten_power_x, "sigma"=inv.logit)
    hyper_fun_inv <- list("C"=log10, "sigma"=logit)
    transformed_hyper <- transform_hyper(list_hyper, hyper_fun)
  } else if (algo=="knn") {
    names_hyper <- c("k")
    n_partition <- 5
    strides_hyper <- list("k"=4)
    n_extra_leftright <- 3
    n_fine_tuning <- list("k"=3)
    min_hyper <- list("k"=1)
    max_hyper <- list("k"=ncol(x_train))
    is_int_hyper <- list("k"=T)
    list_hyper <- list("k"=c(1, seq(strides_hyper[["k"]],strides_hyper[["k"]]*(n_partition-1),by=strides_hyper[["k"]])))
    hyper_fun <- list("k"=identity)
    hyper_fun_inv <- list("k"=identity)
    transformed_hyper <- transform_hyper(list_hyper, hyper_fun)
  } else if (algo=="nb") {
    names_hyper <- c("adjust", "fL")
    use_kernel = c(TRUE, FALSE)
    strides_hyper <- list("adjust"=0.5, "fL"=0.5)
    n_extra_leftright <- 3
    n_fine_tuning <- list("adjust"=3, "fL"=3)
    min_hyper <- list("adjust"=-Inf, "fL"=-Inf)
    max_hyper <- list("adjust"=Inf, "fL"=Inf)
    is_int_hyper <- list("adjust"=F, "fL"=F)
    list_hyper <- list("adjust"=seq(-2,1,by=strides_hyper[["adjust"]]), "fL"= seq(-2,1,by = strides_hyper[["fL"]]))
    hyper_fun <- list("adjust"=ten_power_x, "fL"=ten_power_x)
    hyper_fun_inv <- list("adjust"=log10, "fL"=log10)
    transformed_hyper <- transform_hyper(list_hyper, hyper_fun)
    transformed_hyper$usekernel <- use_kernel
  }
  hyper_grid <- do.call(expand.grid, args = transformed_hyper)
  hyper_grid_all <- hyper_grid
  print("starting algo")
  names(data_train)[which(names(data_train) == target)] <- "target"
  names(data_test)[which(names(data_test) == target)] <- "target"
  data_train$target <-  make.names(data_train$target)
  data_test$target <-  make.names(data_test$target)
  names_predictors <- names(data_train)[-which(names(data_train) %in% names_variables)]
  if(analysis=="age_at_collection") {
    metric <- "Rsquared"
  } else if (analysis %in% c("country", "HvsFDvsD")) {
    metric <- "Accuracy"
  } else {
    metric <- "ROC" }
  args_trctrl <- list(method = "cv", number = N_CV_Folds, classProbs = TRUE, summaryFunction = twoClassSummary, seeds = cv_seeds)
  trctrl <- do.call(trainControl, args_trctrl)
  hyper_grid <- do.call(expand.grid, args = transformed_hyper)
  args_model <- list("form"= formula(target ~.), "data"=data_train, "weights"=w_train, "method"=algo, "trControl"=trctrl, "metric"=metric, "tuneGrid"=hyper_grid, verbose=FALSE)
  model <- do.call(train, args=args_model)
  results <- model$results
  best_tune <- model$bestTune
  top_performance <- max(results[,metric],na.rm=TRUE)
  print("Model:"); print(model); print("Results:"); print(results); print("Best tune:"); print(best_tune)
  #coarse grid search: keep shifting the grid until the best hyperparameters are not an extremum value
  no_extremum <- FALSE
  while(!no_extremum)
  {
    no_extremum <- TRUE #if no extremum is detected for this "epoch", skip the while() the next time
    #for each hyperparameter, check if the best value is an extremum.
    for (hyper in names_hyper)
    {
      results_filtered <- results
      for (hyper_other in names_hyper[-which(names_hyper==hyper)]) #filter the results: only consider the rows of the results for which the other hyperparameters are best tuned.
      {
        best_hyper_other <- results_filtered[,hyper_other][which(results_filtered[,hyper_other]==best_tune[[hyper_other]])]
        results_filtered <- results_filtered[which(results_filtered[,hyper_other] %in% best_hyper_other),]
      }
      top_2nd_performance <- max(results_filtered[-which(results_filtered[,hyper]==best_tune[[hyper]]),metric],na.rm=TRUE)
      coarse_search_hyper <- best_tune[[hyper]] %in% c(min(results_filtered[,hyper]), max(results_filtered[,hyper])) & top_2nd_performance < top_performance
      while(coarse_search_hyper)
      {
        if(best_tune[[hyper]] == min(results_filtered[,hyper]) & min(hyper_fun_inv[[hyper]](results_filtered[,hyper])) > min_hyper[[hyper]]){
          print(paste("WARNING! best model selected the lowest ", hyper, " value.", sep ="")); print(best_tune[[hyper]])
          list_hyper[[hyper]] <- seq(from=min(hyper_fun_inv[[hyper]](results_filtered[,hyper]))- strides_hyper[[hyper]], by=-strides_hyper[[hyper]], length.out=n_extra_leftright)
          list_hyper[[hyper]][which(list_hyper[[hyper]] < min_hyper[[hyper]])] <- min_hyper[[hyper]]
        } else if (best_tune[[hyper]] == max(results_filtered[,hyper]) & max(hyper_fun_inv[[hyper]](results_filtered[,hyper])) < max_hyper[[hyper]]){
          print(paste("WARNING! best model selected the highest ", hyper, " value.", sep ="")); print(best_tune[[hyper]])
          list_hyper[[hyper]] <- seq(from=max(hyper_fun_inv[[hyper]](results_filtered[,hyper])) + strides_hyper[[hyper]], by=strides_hyper[[hyper]], length.out=n_extra_leftright)
          list_hyper[[hyper]][which(list_hyper[[hyper]] > max_hyper[[hyper]])] <- max_hyper[[hyper]]
        } else {
          coarse_search_hyper <- FALSE
        }
        if(coarse_search_hyper)
        {
          list_hyper[[hyper]] <- unique(list_hyper[[hyper]])
          transformed_hyper <- transform_hyper(list_hyper, hyper_fun)
          if(algo=="gbm"){
            transformed_hyper$n.trees <- n_trees
          } else if (algo=="svmPoly") {
            transformed_hyper$degree <- n_poly
          } else if (algo=="nb"){
            transformed_hyper$usekernel <- use_kernel}
          hyper_grid <- do.call(expand.grid, args = transformed_hyper)
          hyper_grid <- hyper_grid[index_different_rows(hyper_grid, hyper_grid_all),,drop=FALSE]
          hyper_grid_all <- rbind(hyper_grid_all, hyper_grid)
          if(nrow(hyper_grid) < 1){coarse_search_hyper <- FALSE}
        }
        if(coarse_search_hyper)
        {
          args_model <- list("form"= formula(target ~.), "data"=data_train, "weights"=w_train, "method"=algo, "trControl"=trctrl, "metric"=metric, "tuneGrid"=hyper_grid, verbose=FALSE)
          model2 <- do.call(train, args=args_model)
          if(max(model2$results[,metric],na.rm=TRUE) <= top_performance){
            print(paste("No better value was found for ", hyper, ".", sep =""))
            coarse_search_hyper <- FALSE
          } else {
            model <- model2
            best_tune <- model$bestTune
            results <- rbind(results, model$results)
            top_performance <- max(results[,metric],na.rm=TRUE)
            print("Best tune:"); print(best_tune); print(paste("Performance using metric ", metric, " = ", round(top_performance,4), sep=""))
            no_extremum <- FALSE #Keep looping the great no_extremum loop because new parameters ranges were introduced, so it is possible that another hyperparameter is now on an extremum.
            #check if the new best hyperparameter is an extremum, of if looping can stop for this hyperparameter.
            results_filtered <- results
            for (hyper_other in names_hyper[-which(names_hyper==hyper)]) #filter the results: only consider the rows of the results for which the other hyperparameters are best tuned.
            {
              best_hyper_other <- results_filtered[,hyper_other][which(results_filtered[,hyper_other]==best_tune[[hyper_other]])]
              results_filtered <- results_filtered[which(results_filtered[,hyper_other] %in% best_hyper_other),]
            }
            top_2nd_performance <- max(results_filtered[-which(results_filtered[,hyper]==best_tune[[hyper]]),metric],na.rm=TRUE)
            coarse_search_hyper <- best_tune[[hyper]] %in% c(min(results_filtered[,hyper]), max(results_filtered[,hyper])) & top_2nd_performance < top_performance
          } #end of ifelse a better model is found => update best model. otherwise stop coarse search.
        } else {
          print(paste("No better value was found for ", hyper, ".", sep =""))
          coarse_search_hyper <- FALSE
        } #end of if nrows > 0 => update. Otherwise stop coarse search.
      } #end of coarse search loop for a specific parameter.
    } #end of loop over hyperparamters for coarse search loop.
  }
  #zoom grid search: zoom around the best hyperparameters values to fine tune the grid
  for (hyper in names_hyper)
  {
    center_hyper <- hyper_fun_inv[[hyper]](best_tune[[hyper]])
    best_hypers <- hyper_fun_inv[[hyper]](results[,hyper][which(results[,metric]==top_performance)])
    min_seq <- min(best_hypers)
    max_seq <- max(best_hypers)
    #if there is a plateau of best parameters, explore it. Otherwise, look for best second value. In the other direction, explore one stride.
    if (min_seq == max_seq)
    {
      results_filtered_without <- results[-which(results[,hyper]==best_tune[[hyper]]),]
      hyper_2nd_bests <- hyper_fun_inv[[hyper]](results_filtered_without[,hyper][which(results_filtered_without[,metric]==max(results_filtered_without[,metric]))])
      hyper_2nd_best <- hyper_2nd_bests[which.min(abs(hyper_2nd_bests-center_hyper))]
      if(hyper_2nd_best > center_hyper)
      {
        min_seq <- max(min_hyper[[hyper]], center_hyper - strides_hyper[[hyper]])
        max_seq <- hyper_2nd_best
      } else {
        min_seq <- hyper_2nd_best
        max_seq <- min(max_hyper[[hyper]], center_hyper + strides_hyper[[hyper]])
      }
    }
    if(min_seq==-Inf) min_seq <- -100
    if(max_seq==Inf) max_seq <- 100
    list_hyper[[hyper]] <- seq(min_seq, max_seq, by=(max_seq-min_seq)/(n_fine_tuning[[hyper]]+1))
    if(is_int_hyper[[hyper]]){list_hyper[[hyper]] <- unique(round(list_hyper[[hyper]]))}
  }
  transformed_hyper <- transform_hyper(list_hyper, hyper_fun)
  if(algo=="gbm"){
    transformed_hyper$n.trees <- n_trees
  } else if (algo=="svmPoly") {
    transformed_hyper$degree <- n_poly
  } else if (algo=="nb"){
    transformed_hyper$usekernel <- use_kernel}
  hyper_grid <- do.call(expand.grid, args = transformed_hyper)
  hyper_grid <- hyper_grid[index_different_rows(hyper_grid, hyper_grid_all),,drop=FALSE]
  if(nrow(hyper_grid) > 0)
  {
    hyper_grid_all <- rbind(hyper_grid_all, hyper_grid)
    args_model <- list("form"= formula(target ~.), "data"=data_train, "weights"=w_train, "method"=algo, "trControl"=trctrl, "metric"=metric, "tuneGrid"=hyper_grid, verbose=FALSE)
    model2 <- do.call(train, args=args_model)
    if(max(model2$results[,metric],na.rm=TRUE) <= top_performance){
      print("No better value was found during the fine search.")
    } else {
      model <- model2
      best_tune <- model$bestTune
      results <- rbind(results, model$results)
      top_performance <- max(results[,metric],na.rm=TRUE)
    }
  } else {
    print("No better value was found during the fine search.")
  }
  print(results); print("Best tune:"); print(best_tune); print(paste("Performance using metric ", metric, " = ", round(top_performance,4), sep=""))
  return(list("model"=model))
}


