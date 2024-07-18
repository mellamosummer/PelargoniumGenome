#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=srb67793@uga.edu
#SBATCH --output=HiFiASM.%j.out
#SBATCH --error=HiFiASM.%j.err
#SBATCH --job-name=HiFiASM_72024
#SBATCH --partition=batch
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=180GB

#script modified from https://gitlab.com/ficklinlab-public/wa-38-genome/-/tree/main/04-nuclear_assembly?ref_type=heads
#set output directory variable
READSDIR="/scratch/srb67793/Pelargonium/HiFiReads"
OUTDIR="/scratch/srb67793/Pelargonium/Results"
HICREADSDIR="/scratch/srb67793/Pelargonium/HiCReads"

#load modules 
module load hifiasm/0.19.6-GCCcore-11.3.0

#create the phased assembly at contig level using HiFiASM
mkdir -p $OUTDIR/HiFiASM_OmniC
hifiasm -o $OUTDIR/HiFiASM_OmniC/P_citronellum_BF73_OmniC.asm -t 48 $READSDIR --h1 $HICREADSDIR/HiCLibs_1.fastq --h2 $HICREADSDIR/HiCLibs_2.fastq $READSDIR/*gz
