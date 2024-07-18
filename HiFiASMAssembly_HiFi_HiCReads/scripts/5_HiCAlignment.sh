#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=srb67793@uga.edu
#SBATCH --output=HiCAlign.%j.out
#SBATCH --error=HiCAlign.%j.err
#SBATCH --job-name=HiCAlign_72024
#SBATCH --partition=batch
#SBATCH --time=4-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=180GB

ml SAMtools 
ml BWA

GENOMEDIR="/scratch/srb67793/Pelargonium/Results/HiFiASM_OmniC"
HICDIR="/scratch/srb67793/Pelargonium/HiCReads"

#HiC data steps
#Step 1: The BWA aligner requires that the input sequences are indexed.
#bwa index -a bwtsw $GENOMEDIR/P_citronellum_OmniC.hic.bp.hap1.p_ctg.fa
#bwa index -a bwtsw $GENOMEDIR/P_citronellum_OmniC.hic.bp.hap2.p_ctg.fa

#Step 2: Perform the BWA alignment
#bwa mem -t 24 -5SP \
#    $GENOMEDIR/P_citronellum_OmniC.hic.bp.hap1.p_ctg.fa \
#    $HICDIR/HiCLibs_1.fastq \
#    $HICDIR/HiCLibs_2.fastq > $GENOMEDIR/P_citronellum.asm.hic.hap1.p_ctg_vs_HiC.sam

#bwa mem -t 24 -5SP \
#    $GENOMEDIR/P_citronellum_OmniC.hic.bp.hap2.p_ctg.fa \
#    $HICDIR/HiCLibs_1.fastq \
#    $HICDIR/HiCLibs_2.fastq > $GENOMEDIR/P_citronellum.asm.hic.hap2.p_ctg_vs_HiC.sam

#Step 3: The BWA aligner creates an alignment file in SAM format. Next, convert the SAM alignment to a BAM file
#samtools view -h -b -F 2316 $GENOMEDIR/P_citronellum.asm.hic.hap1.p_ctg_vs_HiC.sam > $GENOMEDIR/P_citronellum.asm.hic.hap1.p_ctg_vs_HiC.bam
#samtools view -h -b -F 2316 $GENOMEDIR/P_citronellum.asm.hic.hap2.p_ctg_vs_HiC.sam > $GENOMEDIR/P_citronellum.asm.hic.hap2.p_ctg_vs_HiC.bam

#Step 4: the BAM file can be sorted which makes lookups easier for other applications.
samtools sort -n  $GENOMEDIR/P_citronellum.asm.hic.hap1.p_ctg_vs_HiC.bam > $GENOMEDIR/P_citronellum.asm.hic.hap1.p_ctg_vs_HiC.sorted.bam
samtools sort -n  $GENOMEDIR/P_citronellum.asm.hic.hap2.p_ctg_vs_HiC.bam > $GENOMEDIR/P_citronellum.asm.hic.hap2.p_ctg_vs_HiC.sorted.bam

#QC interactively using singularity
#hic_qc.py -b G_maculatum_BF73.asm.hic.hap1.p_ctg_vs_HiC.sorted.bam -o G_maculatum_BF73.hap1_vs_HiC
#hic_qc.py -b G_maculatum_BF73.asm.hic.hap2.p_ctg_vs_HiC.sorted.bam -o G_maculatum_BF73.hap2_vs_HiC

