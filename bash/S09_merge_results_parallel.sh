#!/bin/bash
analyses=( "age_at_collection" "abx_usage" "exclusive_bf" "delivery_type" "sex" "country" "HvsFD" "HvsD" "FDvsD" "HvsFDD" "HvsFDvsD" "Surv" )
analyses=( "Surv" "exclusive_bf" "delivery_type" )
sets=( "train" "test" )
#sets=( "test" )
types=( "Performance" "Performance_sd" )
#types=( "Performance" )
metrics_regression=( "R2" )
metrics_classification=( "ROC" "Cross_Entropy" "Mean_Accuracy" "Accuracy" "Sensitivity" "Specificity" )
metrics_country=( "Cross_Entropy" "Mean_Accuracy" "Accuracy" "SWE" "EST" "FIN" "RUS" )
metrics_HvsFDvsD=( "Cross_Entropy" "Mean_Accuracy" "Accuracy" "H" "FD" "D" )
metrics_survival=( "CI" )
for analysis in "${analyses[@]}"
do
if [ $analysis = "age_at_collection" ]; then
metrics=("${metrics_regression[@]}")
elif [ $analysis = "country" ]; then
metrics=("${metrics_country[@]}")
elif [ $analysis = "HvsFDvsD" ]; then
metrics=("${metrics_HvsFDvsD[@]}")
elif [ $analysis = "Surv" ]; then
metrics=("${metrics_survival[@]}")
else
metrics=("${metrics_classification[@]}")
fi
for metric in "${metrics[@]}"
do
for set in "${sets[@]}"
do
for type in "${types[@]}"
do
job_name="mr-$analysis-$metric-$set-$type.job"
out_file="./../eo/mr-$analysis-$metric-$set-$type.out"
err_file="./../eo/mr-$analysis-$metric-$set-$type.err"
memory=200
n_cores=1
time=1
sbatch --error=$err_file --output=$out_file --job-name=$job_name --mem-per-cpu=$memory -c $n_cores -t $time merge_results.sh $analysis $metric $set $type
done
done
done
done

