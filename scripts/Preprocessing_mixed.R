#Preprocessing: generating the preprocessed dataset for which the predictors are the best predictors between the cags, the pathways and the taxa

#extract and use info from command line
args = commandArgs(trailingOnly=TRUE)
if(length(args) == 0){args <- c("age_at_collection", "1")} #default, if no command line received
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

#generate the dataset
print("start")
preprocessing_mixed(analysis, i)
print("done")
