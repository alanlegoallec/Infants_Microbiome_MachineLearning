#Testing microbiome

#extract and use info from command line
args = commandArgs(trailingOnly=TRUE)
exclude_other_predictors_from_analysis = FALSE

#default, if no command line received
if(length(args) == 0){args <- c("country", "cags", "rf2", "10", "test")}
if(length(args) != 5)
{
  stop("use five arguments for analysis")
} else {
  analysis <- args[1]
  predictors <- args[2]
  algo <- args[3]
  i <- args[4]
  set <- args[5]
}
if(algo== "same"){algo <- predictors} #when analyzing at the gene level, use the genes that were selected for the algorithm.

path <- ifelse(getwd() == "/Users/Alan", "~/Desktop/Aging/Microbiome/data/", "/n/groups/patel/Alan/Aging/Microbiome/data/")
source(paste(path, "../scripts/Helpers_microbiome.R", sep = ""))

#set target variable
target <- target_of_analysis(analysis)

#process
print(paste("Analysis is ", analysis, ".", sep=""))
print(paste("Target variable is ", target, ".", sep=""))
print(paste("Predictors is ", predictors, ".", sep=""))
print(paste("Algorithm is ", algo, ".", sep=""))
print(paste("Fold is ", i, ".", sep=""))
print(paste("Set is ", set, ".", sep=""))


#load data
ifelse(predictors=="genes", string <- algo, string <- predictors)
preprocessed_data <- readRDS(paste(path, "preprocessed_data_", analysis, "_", string, "_", i, ".Rda", sep = ""))
data <- preprocessed_data[[paste("data", set, sep = "_")]]
names(data)[which(names(data) == target)] <- "target"
x <- preprocessed_data[[paste("x", set, sep = "_")]]
y <- preprocessed_data[[paste("y", set, sep = "_")]]
ifelse(analysis=="Surv", w <- preprocessed_data[[paste("w", set, sep = "_")]], w <- generate_weights(y))

#run analysis
if(analysis=="Surv")
{
  testing_surv(analysis, target, predictors, algo, i, data, x, y, w, n_cores)
} else {
  if(target == "age_at_collection" & algo == "nb"){
    print("Naive Bayes is not available for regression.")
  } else {
    testing(analysis, target, predictors, algo, i, data, x, y, w, set)
  }
}
print("done")


