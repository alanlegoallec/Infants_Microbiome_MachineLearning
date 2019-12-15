Explanations for the pipeline on O2 for the microbiome project

Step 1: preprocessing_parallel.sh to preprocess the data
Step 2: preprocessing_folds_parallel.sh to preprocess each CV fold of the data
Step 3: preprocessing_demographics_and_mixed_parallel.sh to preprocess the datasets for which predictors=demographics, and without all predictors types combined
Step 4: preprocessing_nodemo_parallel.sh to preprocess the datatsets without demographics
Step 5: training_parallel.sh to train the models
Step 6: testing_parallel.sh to generate the predictions on both the training and the testing set
Step 7: postprocessing_parallel.sh to calculate the R2/accuracy/ROC/cross-entropy for each fold and on the entire dataset
Step 8: merge_hyperparameters_parallel.sh to merge the choices of hyperparameters, along the accuracies and sample size
Step 9: merge_results.sh to merge the prediction performances for each target, predictors and algorithm


