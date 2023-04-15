#!/bin/bash --login

module load BBMap/38.63

DATADIR=/hpc/group/bio556l-s23/jal138/phyllosphere_evolution/data/SRA_files
for i in $DATADIR/*_paired/*.fastq.gz
do
  echo $i
  BASE=${i%.fastq.gz}
  reformat.sh in=$i out1="$BASE"_1.fastq.gz out2="$BASE"_2.fastq.gz overwrite=t
done
