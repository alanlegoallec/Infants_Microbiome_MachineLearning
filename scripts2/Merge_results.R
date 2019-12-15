#Microbiome merge results

#extract and use info from command line
args = commandArgs(trailingOnly=TRUE)
#default, if no command line received
if(length(args) == 0){args <- c("country", "Cross_Entropy", "train", "Performance")}
print(args)
print(length(args))

if(length(args) != 4)
{
  stop("use 4 arguments for analysis")
} else {
  analysis <- args[1]
  metric <- args[2]
  set <- args[3]
  type <- args[4]
}


path <- ifelse(getwd() == "/Users/Alan", "~/Desktop/Aging/Microbiome/data/", "/n/groups/patel/Alan/Aging/Microbiome/data/")
source(paste(path, "../scripts/Helpers_microbiome.R", sep = ""))

print("starting")
Performance <- merge_results(type, analysis, metric, set)
print(Performance)
print("done")


