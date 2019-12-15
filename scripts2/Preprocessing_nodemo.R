#Preprocessing: generating the preprocessed dataset for which the predictors are only the demographics variables

#set path
path <- ifelse(getwd() == "/Users/Alan", "~/Desktop/Aging/Microbiome/data/", "/n/groups/patel/Alan/Aging/Microbiome/data/")
source(paste(path, "../scripts/Helpers_microbiome.R", sep = ""))


print("start")
for (analysis in analyses[-which(analyses=="Surv")])
{
  target <- target_of_analysis(analysis)
  names_variables <- c("index", "subjectID", "seqID", "age_at_collection", "delivery_type", "sex", "country", "country.EST", "country.RUS", "country.SWE", "exclusive_bf", "abx_usage", "is_mother", "seroconverted_ever", "seroconverted_time", "D_status", "target", "time_to_onset")
  to_remove <- names_variables[which(!(names_variables==target))]
  for (predictors in predictorsS_withdemo)
  {
    if (predictors == "genes+demo")
    {
      for (algo in algos_genes)
      {
        for (i in seq(0,N_CV_Folds))
        {
          preprocessed_data <- readRDS(paste(path, "preprocessed_data_", analysis, "_", algo, "_", i, ".Rda", sep = ""))
          preprocessed_data$data_train <- preprocessed_data$data_train[,-which(names(preprocessed_data$data_train) %in% to_remove)]
          preprocessed_data$data_test <- preprocessed_data$data_train[,-which(names(preprocessed_data$data_test) %in% to_remove)]
          index_x <- which(names(data.frame(preprocessed_data$x_train)) %in% to_remove)
          preprocessed_data$x_train <- preprocessed_data$x_train[,-index_x]
          preprocessed_data$x_test <- preprocessed_data$x_train[,-index_x]
          saveRDS(preprocessed_data, paste(path, "preprocessed_data_", analysis, "_", gsub("\\+.*","", algo), "_", i, ".Rda", sep = ""))
        }
      }
    } else {
      for (i in seq(0,N_CV_Folds))
      {
        preprocessed_data <- readRDS(paste(path, "preprocessed_data_", analysis, "_", predictors, "_", i, ".Rda", sep = ""))
        preprocessed_data$data_train <- preprocessed_data$data_train[,-which(names(preprocessed_data$data_train) %in% to_remove)]
        preprocessed_data$data_test <- preprocessed_data$data_train[,-which(names(preprocessed_data$data_test) %in% to_remove)]
        index_x <- which(names(data.frame(preprocessed_data$x_train)) %in% to_remove)
        preprocessed_data$x_train <- preprocessed_data$x_train[,-index_x]
        preprocessed_data$x_test <- preprocessed_data$x_train[,-index_x]
        saveRDS(preprocessed_data, paste(path, "preprocessed_data_", analysis, "_", gsub("\\+.*","", predictors), "_", i, ".Rda", sep = ""))
      }
    }
  }
}
print("done")
