#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --time=0:30:00
#SBATCH --job-name=star_genome_generate
#SBATCH --output=output.out
#SBATCH --error=output.err

cd $SLURM_SUBMIT_DIR

module load CCEnv arch/avx512
module load StdEnv/2020
module load star/2.7.8a



# Generate the genome files for wbcel235.104_nc.gtf
STAR --runMode genomeGenerate --genomeDir $SCRATCH/genomeIndex/wbcel235.104_all --genomeFastaFiles $SCRATCH/wbcel235_genome/I.fa $SCRATCH/wbcel235_genome/II.fa $SCRATCH/wbcel235_genome/III.fa $SCRATCH/wbcel235_genome/IV.fa $SCRATCH/wbcel235_genome/V.fa $SCRATCH/wbcel235_genome/X.fa $SCRATCH/wbcel235_genome/MtDNA.fa --sjdbGTFfile $SCRATCH/annotations/wbcel235.104.gtf --sjdbOverhang 100 --runThreadN 40 --genomeSAindexNbases 12
