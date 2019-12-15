#Summarize results microbiome

path <- ifelse(getwd() == "/Users/Alan", "~/Desktop/Aging/Microbiome/data/", "/n/groups/patel/Alan/Aging/Microbiome/data/")
source(paste(path, "../scripts/Helpers_microbiome.R", sep = ""))

display_train <- FALSE
ifelse(display_train, sets <- c("train", "test"), sets <- c("test"))
display_sd <- FALSE
display_surv <- TRUE

for (analysis in c("age_at_collection"))
{
  print(analysis)
  print("R2S and R2Ssd")
  for (predictors in predictorsS)
  {
    print(predictors)
    for (type in sets)
    {
      print(type)
      print(readRDS(paste(path, "R2S_", type, "_", analysis, "_", predictors, ".Rda", sep = "")))
    }
  }
}

for (analysis in analyses[-c(1,length(analyses))])
{
  print(analysis)
  for (predictors in predictorsS)
  {
    print(predictors)
    for (type in sets)
    {
      print(type)
      print("Performance")
      print(readRDS(paste(path, "Performance_", type, "_", analysis, "_", predictors, ".Rda", sep = "")))
      if(display_sd)
      {
        print("Performance SD")
        print(readRDS(paste(path, "Performance_sd_", type, "_", analysis, "_", predictors, ".Rda", sep = "")))
      }
    }
  }
}

if(display_surv)
{
  for (analysis in c("Surv"))
  {
    print(analysis)
    print("CIS and CISsd")
    for (predictors in predictorsS)
    {
      print(predictors)
      for (type in sets)
      {
        print(type)
        print(readRDS(paste(path, "CIS_", type, "_", analysis, "_", predictors, ".Rda", sep = "")))
      }
    }
  }
}




