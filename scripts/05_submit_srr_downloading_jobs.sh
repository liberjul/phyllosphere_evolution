#!/bin/bash --login

while IFS= read -r line
do
  TARGET=$(echo $line | cut -f1 -d' ')
  REGION=$(echo $line | cut -f2 -d' ')
  LAYOUT=$(echo $line | cut -f3 -d' ')
  sed "s/<target>/$TARGET/" download_sra_files_job.sb > download_sra_files_job_"$TARGET"_"$REGION"_"$LAYOUT".sb
  sed -i "s/<region>/$REGION/" download_sra_files_job_"$TARGET"_"$REGION"_"$LAYOUT".sb
  sed -i "s/<layout>/$LAYOUT/" download_sra_files_job_"$TARGET"_"$REGION"_"$LAYOUT".sb
  sbatch download_sra_files_job_"$TARGET"_"$REGION".sb
done < ../data/metadata/srr_targets_regions.txt
