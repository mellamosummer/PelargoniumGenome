#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=srb67793@uga.edu
#SBATCH --output=FCS_Pcitr.%j.out
#SBATCH --error=FCS_Pcitr.%j.err
#SBATCH --job-name=FCS_Pcitr_72024
#SBATCH --partition=highmem_p
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=500GB

#script modified from https://gitlab.com/ficklinlab-public/wa-38-genome/-/tree/main/04-nuclear_assembly?ref_type=heads
#set output directory variable
#HIFIREADS="/scratch/srb67793/Pelargonium/HiFiReads/m64173_210719_082812.ccs.fasta.gz"
OUTDIR="/scratch/srb67793/Pelargonium/Results"

#load modules
#module load hifiasm/0.19.6-GCCcore-11.3.0
#module load BUSCO/5.5.0-foss-2022a
#module load ncbi-FCS/0.5.0-GCCcore-11.3.0

#create the phased assembly at contig level using HiFiASM
#hifiasm -o $OUTDIR -t 24 $HIFIREADS

#make output into fasta file format
#awk '/^S/{print ">"$2;print $3}' Results.bp.p_ctg.gfa > Results.bp.p_ctg.fa
#awk '/^S/{print ">"$2;print $3}' Results.bp.hap1.p_ctg.gfa > Results.bp.hap1.p_ctg.fa
#awk '/^S/{print ">"$2;print $3}' Results.bp.hap2.p_ctg.gfa > Results.bp.hap2.p_ctg.fa

#get assembly statistics using assemblathon

./assemblathon_stats.pl  Results.bp.p_ctg.fa > Results.bp.p_ctg.txt
./assemblathon_stats.pl Results.bp.hap1.p_ctg.fa > Results.bp.hap1.p_ctg.fa.txt
./assemblathon_stats.pl Results.bp.hap2.p_ctg.fa > Results.bp.hap2.p_ctg.fa.txt

#estimate assembly completeness using BUSCO
busco  \
    -f -i $OUTDIR/HiFiASM/HaplotypeFastas/Results.bp.hap1.p_ctg.fa \
    -m genome \
    -o Pelargonium-hap1-busco \
    -l eudicots_odb10 \
    -c 20

  busco  \
    -f -i $OUTDIR/HiFiASM/HaplotypeFastas/Results.bp.hap2.p_ctg.fa \
    -m genome \
    -o Pelargonium-hap2-busco \
    -l eudicots_odb10 \
    -c 20

#FCS contamination report
python $EBROOTNCBIMINFCS/fcs.py screen genome --fasta Results.bp.hap1.p_ctg.fa --out-dir . --gx-db "$GXDB_LOC/gxdb" --tax-id 73188
python $EBROOTNCBIMINFCS/fcs.py screen genome --fasta Results.bp.hap2.p_ctg.fa --out-dir . --gx-db "$GXDB_LOC/gxdb" --tax-id 73188


