#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=srb67793@uga.edu
#SBATCH --output=OMNIC_QC.%j.out
#SBATCH --error=OMNIC_QC.%j.err
#SBATCH --job-name=OMNIC_QC72024
#SBATCH --partition=batch
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=10GB

HICREADSDIR="/scratch/srb67793/Pelargonium/HiCReads"

#load modules
module load FastQC/0.11.9-Java-11
module load MultiQC/1.14-foss-2022a
#ml fastp

#QC Hi-C Reads
mkdir -p fastqc
mkdir -p multiqc
#fastqc -t 2 $HICREADSDIR/*fastq -o fastqc
multiqc fastqc/ -o multiqc/

#Trim HIC Reads
#fastp \
#    -w 16 \
#    --disable_quality_filtering \
#    --disable_length_filtering \
#    --disable_trim_poly_g \
#    -i $HICREADSDIR/HiCLib1_1.fastq \
#    -o $HICREADSDIR/Hic_Lib1_R1.trimmed.fastq \
#    -I $HICREADSDIR/HiCLib1_2.fastq
#    -O $HICREADSDIR/Hic_Lib1_R2.trimmed.fastq

#fastp \
#    -w 16 \
#    --disable_quality_filtering \
#    --disable_length_filtering \
#    --disable_trim_poly_g \
#    -i $HICREADSDIR/HiCLib2_2.fastq -o $HICREADSDIR/Hic_Lib2_R1.trimmed.fastq -I $HICREADSDIR/HiCLib2_2.fastq -O $HICREADSDIR/Hic_Lib2_R2.trimmed.fastq

#QC PacBio reads
#mkdir -p $OUTDIR/fastqc $OUTDIR/multiqc
#fastqc -t 2 -o $OUTDIR/fastqc $READSDIR/*gz
#multiqc $OUTDIR/fastqc/* -o $OUTDIR/multiqc/
