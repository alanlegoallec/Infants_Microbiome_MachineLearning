#!/bin/bash
analyses=( "age_at_collection" "abx_usage" "exclusive_bf" "delivery_type" "sex" "country" "HvsFD" "HvsD" "FDvsD" "HvsFDD" "HvsFDvsD" "Surv" )
analyses=( "Surv" )
PREDICTORS=( "taxa" "cags" "pathways" "glmnet" "glmnet2" "gbm" "gbm2" "rf" "rf2" )
#PREDICTORS=( "taxa" )
PREDICTORS_surv=( "taxa" "cags" "pathways" "glmnet2" "gbm2" "rf2" )
#PREDICTORS_surv=( "cags" )
memory=2G
n_cores=1
time=5
for analysis in "${analyses[@]}"
do
echo ANALYSIS
echo $analysis
if [ $analysis = "Surv" ]; then predictors=("${PREDICTORS_surv[@]}"); else predictors=("${PREDICTORS[@]}"); fi
for predictor in "${predictors[@]}"
do
echo PREDICTORS
echo $predictor
job_name="pr-$analysis-$predictor.job"
out_file="./../eo/pr-$analysis-$predictor.out"
err_file="./../eo/pr-$analysis-$predictor.err"
#only run job if the data has not already been saved
if [ ! -f /n/scratch2/al311/Aging/Microbiome/data/indices_test_${analysis}_${predictor}_10.Rda ]; then
   echo "indices_test_${analysis}_${predictor}_10.Rda has not been computed."
   sbatch --error=$err_file --output=$out_file --job-name=$job_name --mem-per-cpu=$memory -c $n_cores -t $time preprocessing.sh $analysis $predictor
else
   echo ok
#   echo "indices_test_${analysis}_${predictor}_10.Rda has already been computed."
fi
done
done

