#!/bin/bash
#SBATCH -A ausei14399
#SBATCH -p genoa
#SBATCH --job-name=STAR_index          # Job name
#SBATCH --nodes=1                     # Run on a single node
#SBATCH --ntasks=1                    # Single task
#SBATCH --cpus-per-task=8             # Number of CPU cores
#SBATCH --mem=45G                     # Memory (Human hg38 needs ~32-40GB)
#SBATCH --time=04:00:00               # Time limit (usually takes 1-2 hours)
#SBATCH --output=logs/star_index_%j.log    # Standard output and error log



# 1. Activate your environment
source ~/miniconda3/bin/activate star_env

# 2. Set your paths (Update these to your actual paths!)
GENOME_DIR="/projects/prjs1616/gatk_paper_demo/resources/star_hg38_index"
FASTA="/projects/prjs1616/gatk_paper_demo/resources/GRCh38.primary_assembly.genome.fa"
GTF="/projects/prjs1616/gatk_paper_demo/resources/gencode.v44.primary_assembly.annotation.gtf"

# 3. Create the output directory, GENOME_DIR = where the output of genome index
mkdir -p $GENOME_DIR

# 4. Run STAR
# Note: --sjdbOverhang 99 is standard for 100bp reads (ReadLength - 1)
STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir $GENOME_DIR \
     --genomeFastaFiles $FASTA \
     --sjdbGTFfile $GTF \
     --sjdbOverhang 99

echo "STAR index generation complete at $(date)"
