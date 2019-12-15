#Processing microbiome

#extract and use info from command line
args = commandArgs(trailingOnly=TRUE)
#default, if no command line received
if(length(args) == 0){args <- c("taxa", "glmnet", "10", "1")}
if(length(args) != 4)
{
  stop("use four arguments for analysis")
} else {
  predictors <- args[1]
  algo <- args[2]
  i <- args[3]
  n_cores <- as.numeric(args[4])
}

path <- ifelse(getwd() == "/Users/Alan", "~/Desktop/Aging/Microbiome/data/", "/n/groups/patel/Alan/Aging/Microbiome/data/")
source(paste(path, "../scripts/Helpers_microbiome.R", sep = ""))

#set target variable
analysis <- "Surv"
target <- target_of_analysis(analysis)

#process
print(algo)
preprocessed_data <- readRDS(paste(path, "preprocessed_data_", analysis, "_", predictors, "_", i, ".Rda", sep = ""))
processing_surv(analysis, target, predictors, algo, i, preprocessed_data$data_train, preprocessed_data$x_train, preprocessed_data$y_train, preprocessed_data$data_test, preprocessed_data$x_test, preprocessed_data$y_test, n_cores)
print("done")


