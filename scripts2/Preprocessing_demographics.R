#Preprocessing: generating the preprocessed dataset for which the predictors are only the demographics variables

#extract and use info from command line
args = commandArgs(trailingOnly=TRUE)
if(length(args) == 0){args <- c("age_at_collection", "2")} #default, if no command line received
if(length(args) != 2)
{
  stop("use two arguments for analysis")
} else {
  analysis <- args[1]
  i <- args[2]
}

#set path
path <- ifelse(getwd() == "/Users/Alan", "~/Desktop/Aging/Microbiome/data/", "/n/groups/patel/Alan/Aging/Microbiome/data/")
source(paste(path, "../scripts/Helpers_microbiome.R", sep = ""))

#set target variable
target <- target_of_analysis(analysis)

#generate preprocessed "demographics" dataset
names_others <- c(target, c("age_at_collection", "sex", "country.EST", "country.RUS", "country.SWE", "delivery_type", "exclusive_bf", "abx_usage", "time_to_onset"))
preprocessed_data <- readRDS(paste(path, "preprocessed_data_", analysis, "_", "cags+demo", "_", i, ".Rda", sep = ""))
preprocessed_data$data_train <- preprocessed_data$data_train[,which(names(preprocessed_data$data_train) %in% names_others)]
preprocessed_data$data_test <- preprocessed_data$data_test[,which(names(preprocessed_data$data_test) %in% names_others)]
x <- data.frame(preprocessed_data$x_train)
index_x <- which(names(x) %in% names_others)
preprocessed_data$x_train <- preprocessed_data$x_train[,index_x]
preprocessed_data$x_test <- preprocessed_data$x_test[,index_x]
saveRDS(preprocessed_data, paste(path, "preprocessed_data_", analysis, "_", "demographics", "_", i, ".Rda", sep = ""))

