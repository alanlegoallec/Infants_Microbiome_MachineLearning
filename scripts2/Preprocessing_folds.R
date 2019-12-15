#Preprocessing microbiome

#extract and use info from command line
args = commandArgs(trailingOnly=TRUE)
if(length(args) == 0){args <- c("age_at_collection", "rf2", "2")} #default, if no command line received
if(length(args) != 3)
{
  stop("use three arguments for analysis")
} else {
  analysis <- args[1]
  predictors <- args[2]
  i <- args[3]
}

#set path
path <- ifelse(getwd() == "/Users/Alan", "~/Desktop/Aging/Microbiome/data/", "/n/groups/patel/Alan/Aging/Microbiome/data/")
source(paste(path, "../scripts/Helpers_microbiome.R", sep = ""))

#set target variable
target <- target_of_analysis(analysis)

#preprocess folds
print(target)
print(i)
if (i == "0")
{
  preprocessing_folds_0(analysis, target, predictors)
} else if(analysis == "Surv") {
  preprocessing_folds_surv(analysis, target, predictors, i)
} else {
  preprocessing_folds(analysis, target, predictors, i)
}
print("done")
