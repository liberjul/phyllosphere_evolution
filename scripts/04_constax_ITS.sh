#!/bin/bash --login
########## Define Resources Needed with SBATCH Lines ##########

#SBATCH --time=00:30:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1                 # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=16          # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=32G                   # memory required per node - amount of memory (in bytes)
#SBATCH --job-name=constax_ITS        # you can give your job a name for easier identification (same as -J)
#SBATCH --output=%x-%j.SLURMout

########## Command Lines to Run ##########

conda activate constax

cd ~/Duke_teaching/BIOL557L_S23/

constax -c 0.8 \
-b -t -i ./dna-sequences_ITS.fasta \
-n $SLURM_CPUS_PER_TASK \
-d /mnt/ufs18/rs-022/bonito_lab/CONSTAX_May2020/UNITE_Fungi_tf/sh_general_release_fungi_44343_RepS_10.05.2021.fasta \
-f /mnt/ufs18/rs-022/bonito_lab/CONSTAX_May2020/UNITE_Fungi_tf --mem $SLURM_MEM_PER_NODE -m 5 \
-x ./taxonomy/tax_fun \
-o ./taxonomy/out_fun \
--high_level_db /mnt/ufs18/home-026/liberjul/Bonito_Lab/constax/test_new/SILVA138_NR_blastdb/SILVA138_NRdb.fasta \
--sintax_path /mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 \
--consistent

conda deactivate

scontrol show job $SLURM_JOB_ID     ### write job information to output file
