#Microbiome postprocessing

#extract and use info from command line
args = commandArgs(trailingOnly=TRUE)
#default, if no command line received
if(length(args) == 0){args <- c("age_at_collection", "demographics", "train")}
if(length(args) != 3)
{
  stop("use three arguments for analysis")
} else {
  analysis <- args[1]
  predictors <- args[2]
  set <- args[3]
}

path <- ifelse(getwd() == "/Users/Alan", "~/Desktop/Aging/Microbiome/data/", "/n/groups/patel/Alan/Aging/Microbiome/data/")
source(paste(path, "../scripts/Helpers_microbiome.R", sep = ""))

print("start")
if (analysis %in% c("country", "HvsFDvsD")){
  post_processing_multinomial(analysis, predictors, set)
} else if (analysis == "age_at_collection"){
  post_processing_regression(analysis, predictors, set)
} else if (analysis == "Surv"){
  post_processing_surv(analysis, predictors, set)
  } else {
  post_processing_binomial(analysis, predictors, set)
}
print("done")


