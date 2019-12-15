#microbiome antibiotics scaring

#set path
path <- ifelse(getwd() == "/Users/Alan", "~/Desktop/Aging/Microbiome/data/", "/n/groups/patel/Alan/Aging/Microbiome/data/")
source(paste(path, "../scripts/Helpers_microbiome.R", sep = ""))
analysis <- "abx_usage"


Sensitivity <- initiate_store(c("all", "more_recent", "less_recent"), algos)
Sensitivity_sd <- Sensitivity
data_full <- readRDS(paste(path, "data_", analysis, ".Rda", sep = ""))
indices <- c()
y_test <- c()
for (i in seq(N_CV_Folds))
{
  indices <- c(indices, readRDS(paste(path, "indices_test_", analysis, "_", i, ".Rda", sep = "")))
  y_test <- c(y_test, readRDS(paste(path, "y_test_", analysis, "_", i, ".Rda", sep = "")))
}
#For each individual, determine the time they took antibiotics.
#I assume they only took antibiotics once, because there is no way for me to estimate if they take it again
#For individuals who were recruited having already taken antibiotics, I take the mean between age 0 and age of first data collection
meta <- read.xlsx(paste(path, "full_metadata_final_cleaned.xlsx", sep = ""))
meta <- meta[!is.na(meta$abx_usage),]
meta <- meta[which(meta$index %in% indices),]
meta$age_abx <- NA
age_abxS <- list()
ids_togetrid <- c()
for (subject_id in unique(meta$subjectID))
{
  meta_s <- meta[which(meta$subjectID == subject_id),]
  if(max(meta_s$abx_usage) == 0)
  {
    ids_togetrid <- c(ids_togetrid, subject_id)
    next
  }
  lower_bound <- ifelse(min(meta_s$abx_usage) == 0, max(meta_s$age_at_collection[which(meta_s$abx_usage == 0)]), 0)
  upper_bound <- min(meta_s$age_at_collection[which(meta_s$abx_usage == 1)])
  age_abx <- (lower_bound+upper_bound)/2
  age_abxS[[subject_id]] <- age_abx
}

for (algo in algos)
{
  pred <- c()
  for (i in seq(N_CV_Folds))
  {
    pred <- c(pred, as.character(readRDS(paste(path, "pred_test_", analysis, "_", algo, "_", i, ".Rda", sep = ""))))
  }
  data <- data.frame(cbind(indices, y_test))
  data$pred <- as.numeric(pred)
  data <- data[which(data$y_test == 1),]
  #all samples
  boot_all <- boot(data = data$pred, statistic = boot_mean, R = 1000)
  Sensitivity["all", algo] <- boot_all$t0
  Sensitivity_sd["all", algo] <- sd(boot_all$t)
  #recent and less recent
  data$age <- data_full$age_at_collection[which(rownames(data_full) %in% data$indices)]
  data$subjectID <- data_full$subjectID[which(rownames(data_full) %in% data$indices)]
  data$age_abx <- NA
  for (subject_id in unique(meta$subjectID))
  {
    data$age_abx[which(data$subjectID == subject_id)] <- age_abxS[[subject_id]]
  }
  data$age_since_abx <- data$age - data$age_abx
  data <- data[order(data$age_since_abx),]
  data <- data[which(data$age_since_abx > 0),] #notes: there is an error for subject_id 397 and 150. They are first recorded as having taken antibiotics. Then later as having never taken antibiotics.
  middle = floor(nrow(data)/2)
  data_more_recent <- data[1:middle,]
  data_less_recent <- data[(middle+1):nrow(data),]
  boot_more_recent <- boot(data = data_more_recent$pred, statistic = boot_mean, R = 1000)
  Sensitivity["more_recent", algo] <- boot_more_recent$t0
  Sensitivity_sd["more_recent", algo] <- sd(boot_more_recent$t)
  boot_less_recent <- boot(data = data_less_recent$pred, statistic = boot_mean, R = 1000)
  Sensitivity["less_recent", algo] <- boot_less_recent$t0
  Sensitivity_sd["less_recent", algo] <- sd(boot_less_recent$t)
}
#print results
print(nrow(data_more_recent))
print(mean(data_more_recent$age_since_abx))
print(sd(data_more_recent$age_since_abx))
print(min(data_more_recent$age_since_abx))
print(max(data_more_recent$age_since_abx))
print(nrow(data_less_recent))
print(mean(data_less_recent$age_since_abx))
print(sd(data_less_recent$age_since_abx))
print(min(data_less_recent$age_since_abx))
print(max(data_less_recent$age_since_abx))
print(Sensitivity)
print(Sensitivity_sd)
saveRDS(Sensitivity, paste(path, "Sensitivity_abx_usage_scaring.Rda", sep = ""))
saveRDS(Sensitivity_sd, paste(path, "Sensitivity_sd_abx_usage_scaring.Rda", sep = ""))


