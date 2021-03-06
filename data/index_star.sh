#!/bin/bash -l
#SBATCH --job-name=star_align_hg38
#SBATCH --output=/cellar/users/hsher/star_index_new-%A_%a.out
#SBATCH --error=/cellar/users/hsher/star_index_new-%A_%a.error
#SBATCH -n 10
#SBATCH --mem-per-cpu=5G
#SBATCH --array=1-20%250





STAR --runThreadN 16 \
--outTmpDir ./_STARtmp-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} \
--runMode genomeGenerate \
--genomeDir /nrnb/users/hsher/star_hg38 \
--genomeFastaFiles /cellar/users/hsher/Data/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--sjdbGTFfile /cellar/users/hsher/Data/reference/Homo_sapiens.GRCh38.98.gtf \
--sjdbOverhang 99
