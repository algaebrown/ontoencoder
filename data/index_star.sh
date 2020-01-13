#!/bin/bash
#SBATCH --job-name=star_align_hg38
#SBATCH --output=/cellar/users/hsher/star_index.out
#SBATCH --error=/cellar/users/hsher/star_index.error
#SBATCH -n 4
#SBATCH --mem-per-cpu=5G
#SBATCH --array=1-4%250
#SBATCH --mail-type=end          
#SBATCH --mail-user=hsher@ucsd.edu





STAR --runThreadN 4 \
--runMode genomeGenerate \
--genomeDir /cellar/users/hsher/Data/star_index \
--genomeFastaFiles /cellar/users/hsher/Data/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--sjdbGTFfile /cellar/users/hsher/Data/reference/Homo_sapiens.GRCh38.98.gtf \
--sjdbOverhang 99