#!/bin/bash --login
#
# x=0
# for FILE in ../data/SRA_files/Fungi*/*R1_001.cutadapt.fastq.gz
# do
#   echo $x
#   n=`printf %04d $x`
#   PREFIX=$(echo $FILE | cut -f1-4 -d'_')
#   SUFFIX=$(echo ${FILE%cutadapt.fastq.gz}fastq.gz | cut -f5-9 -d'_')
#   # echo $PREFIX"$n"_$SUFFIX
#   mv $FILE $PREFIX"$n"_$SUFFIX
#   if [ -f ${FILE%R1_001.cutadapt.fastq.gz}R2_001.cutadapt.fastq.gz ]
#     then echo "Paired reads found"
#     SUFFIX2=${SUFFIX%R1_001.fastq.gz}R2_001.fastq.gz
#     # echo $PREFIX"$n"_$SUFFIX2
#     mv ${FILE%R1_001.cutadapt.fastq.gz}R2_001.cutadapt.fastq.gz $PREFIX"$n"_$SUFFIX2
#   fi
#   x=$((x+1))
# done
#
# qiime tools import \
#   --type 'FeatureData[Sequence]' \
#   --input-path ../data/reference_data/Basidio_yeast_its1_58s_its2_dedup.fasta \
#   --input-format DNAFASTAFormat \
#   --output-path ../data/reference_data/Basidio_yeast_its1_58s_its2.qza

while IFS= read -r line
do
  TARGET=$(echo $line | cut -f1 -d' ')
  REGION=$(echo $line | cut -f2 -d' ')
  LAYOUT=$(echo $line | cut -f3 -d' ')
  sed "s/<target>/$TARGET/" qiime_import_dada2_denoise_job.sb > qiime_import_dada2_denoise_job_"$TARGET"_"$REGION"_"$LAYOUT".sb
  sed -i "s/<region>/$REGION/" qiime_import_dada2_denoise_job_"$TARGET"_"$REGION"_"$LAYOUT".sb
  sed -i "s/<layout>/$LAYOUT/" qiime_import_dada2_denoise_job_"$TARGET"_"$REGION"_"$LAYOUT".sb
  sbatch qiime_import_dada2_denoise_job_"$TARGET"_"$REGION"_"$LAYOUT".sb
done < ../data/metadata/srr_targets_regions.txt
