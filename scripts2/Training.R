#Training microbiome

#extract and use info from command line
args = commandArgs(trailingOnly=TRUE)
exclude_other_predictors_from_analysis = FALSE

#default, if no command line received
if(length(args) == 0){args <- c("age_at_collection", "genes+demo", "glmnet", "0", "1")}
if(length(args) != 5)
{
  stop("use five arguments for analysis")
} else {
  analysis <- args[1]
  predictors <- args[2]
  algo <- args[3]
  i <- args[4]
  n_cores <- as.numeric(args[5])
}
if(algo== "same"){algo <- predictors} #when analyzing at the gene level, use the genes that were selected for the algorithm.

path <- ifelse(getwd() == "/Users/Alan", "~/Desktop/Aging/Microbiome/data/", "/n/groups/patel/Alan/Aging/Microbiome/data/")
source(paste(path, "../scripts/Helpers_microbiome.R", sep = ""))
#set target variable
target <- target_of_analysis(analysis)
#define metric
if(analysis=="age_at_collection") {
  metric <- "Rsquared"
} else if (analysis %in% c("country", "HvsFDvsD")) {
  metric <- "Accuracy"
} else {
  metric <- "ROC" }

#print
print("Analysis:"); print(analysis); print("Target:"); print(target); print("Predictors:"); print(predictors)

#load data
if (predictors=="genes"){
  string <- algo
} else if (predictors=="genes+demo"){
  string <- paste(algo, "demo", sep="+")
} else {
  string <- predictors
}
preprocessed_data <- readRDS(paste(path, "preprocessed_data_", analysis, "_", string, "_", i, ".Rda", sep = ""))
data <- preprocessed_data$data_train
names(data)[which(names(data) == target)] <- "target"
x <- preprocessed_data$x_train
y <- preprocessed_data$y_train
ifelse(analysis=="Surv", w <- preprocessed_data$w_train, w <- generate_weights(y))
if(!(analysis %in% c("age_at_collection", "country", "HvsFDvsD")))
{
  data$target <- make.names(data$target)
  y <- make.names(y)
}


#run analysis
if(analysis=="Surv")
{
  training_surv(analysis, target, predictors, algo, i, data, x, y, w)
} else {
  if(target == "age_at_collection" & algo == "nb"){
    print("Naive Bayes is not available for regression.")
  } else if(algo == "glmnet2") {
    model <- tune_glmnet2(analysis, target, predictors, algo, i, data, x, y, w)
  } else if (algo == "gbm2"){
    model <- tune_gbm2(analysis, target, predictors, algo, i, data, x, y, w)
  } else if (algo == "rf2"){
    model <- tune_rf2(analysis, target, predictors, algo, i, data, x, y, w)
  } else {
    model <- tune_caret_models(analysis, target, predictors, algo, i, data, x, y, w)
  }
  if(!(target == "age_at_collection" & algo == "nb"))
  {
    saveRDS(model, paste(path, "model", "_", analysis, "_", predictors, "_", algo, "_", i,".Rda", sep = ""))
  }
}

#print and save the best predictors if possible
best_predictors(model, analysis, predictors, algo, i)

print("done")



