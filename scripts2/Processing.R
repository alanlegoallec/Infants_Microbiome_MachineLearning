#Processing microbiome

#extract and use info from command line
args = commandArgs(trailingOnly=TRUE)
exclude_other_predictors_from_analysis = FALSE

#default, if no command line received
if(length(args) == 0){args <- c("abx_usage", "demographics", "glmnet", "1", "1")}
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

#process
print(analysis)
print(predictors)
print(target)

#define metric
if(analysis=="age_at_collection") {
  metric <- "Rsquared"
} else if (analysis %in% c("country", "HvsFDvsD")) {
  metric <- "Accuracy"
} else {
  metric <- "ROC" }
#load data
if(predictors == "demographics")
{
  names_others <- c(target, c("age_at_collection", "sex", "country.EST", "country.RUS", "country.SWE", "delivery_type", "exclusive_bf", "abx_usage", "time_to_onset"))
  preprocessed_data <- readRDS(paste(path, "preprocessed_data_", analysis, "_", "cags", "_", i, ".Rda", sep = ""))
  preprocessed_data$data_train <- preprocessed_data$data_train[,which(names(preprocessed_data$data_train) %in% names_others)]
  preprocessed_data$data_test <- preprocessed_data$data_test[,which(names(preprocessed_data$data_test) %in% names_others)]
  x <- data.frame(preprocessed_data$x_train)
  index_x <- which(names(x) %in% names_others)
  preprocessed_data$x_train <- preprocessed_data$x_train[,index_x]
  preprocessed_data$x_test <- preprocessed_data$x_test[,index_x]
} else if (predictors == "genes")
{
  preprocessed_data <- readRDS(paste(path, "preprocessed_data_", analysis, "_", algo, "_", i, ".Rda", sep = ""))
} else {
  preprocessed_data <- readRDS(paste(path, "preprocessed_data_", analysis, "_", predictors, "_", i, ".Rda", sep = ""))
}
names(preprocessed_data$data_train)[which(names(preprocessed_data$data_train) == target)] <- "target"
names(preprocessed_data$data_test)[which(names(preprocessed_data$data_test) == target)] <- "target"



data_train <-  preprocessed_data$data_train
x_train <- preprocessed_data$x_train
y_train <- preprocessed_data$y_train
w_train <- preprocessed_data$w_train
data_test <- preprocessed_data$data_test
x_test <- preprocessed_data$x_test
y_test <- preprocessed_data$y_test
#w_test <- preprocessed_data$w_test
w_test <- seq(y_test)*0+1 #TEMPORARY FIX


indexx <- c()
for (cl in names(table(data_train[["target"]])))
{
  indexx <- c(indexx, which(data_train[["target"]]==cl)[1:28])
}
#indexx <- seq(100)
n <- min(30, ncol(x_train))
data_train <-  preprocessed_data$data_train[,1:(n+1)]
x_train <- preprocessed_data$x_train[,1:n]
y_train <- preprocessed_data$y_train
w_train <- preprocessed_data$w_train
data_train <-  preprocessed_data$data_train[indexx,1:(n+1)]
x_train <- preprocessed_data$x_train[indexx,1:n]
y_train <- preprocessed_data$y_train[indexx]
w_train <- preprocessed_data$w_train[indexx]
data_test <- preprocessed_data$data_test[,1:(n+1)]
x_test <- preprocessed_data$x_test[,1:n]
y_test <- preprocessed_data$y_test

w_test <- seq(y_test)*0+1 #TEMPORARY FIX

#run analysis
processing(analysis, target, predictors, algo, i, data_train, x_train, y_train, w_train, data_test, x_test, y_test, w_test)







#run analysis
if(analysis=="Surv")
{
  processing_surv(analysis, target, predictors, algo, i, preprocessed_data$data_train, preprocessed_data$x_train, preprocessed_data$y_train, preprocessed_data$data_test, preprocessed_data$x_test, preprocessed_data$y_test, preprocessed_data$w_test, n_cores)
} else {
  processing(analysis, target, predictors, algo, i, preprocessed_data$data_train, preprocessed_data$x_train, preprocessed_data$y_train, preprocessed_data$w_train, preprocessed_data$data_test, preprocessed_data$x_test, preprocessed_data$y_test, preprocessed_data$w_test)
}
print("done")



