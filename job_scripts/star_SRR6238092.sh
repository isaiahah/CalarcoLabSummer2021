#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --time=12:00:00
#SBATCH --job-name=star_6238092
#SBATCH --output=output_align.out
#SBATCH --error=output_align.err

cd $SLURM_SUBMIT_DIR

module load CCEnv arch/avx512
module load StdEnv/2020
module load star/2.7.8a

STAR --runThreadN 40 --outSAMtype BAM Unsorted --twopassMode Basic --outFileNamePrefix $SCRATCH/maps/SRR6238092_wbcel235.104_all/ --genomeDir $SCRATCH/genomeIndex/wbcel235.104_all/ --readFilesIn $SCRATCH/read_data/SRR6238092/SRR6238092_1.fastq $SCRATCH/read_data/SRR6238092/SRR6238092_2.fastq --outSJfilterReads Unique
