#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=srb67793@uga.edu
#SBATCH --output=AssemblyStats.%j.out
#SBATCH --error=AssemblyStats.%j.err
#SBATCH --job-name=AssemblyStats_72024
#SBATCH --partition=batch
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=180GB

#Set paths
OUTDIR="/scratch/srb67793/Pelargonium/Results"

#load modules
module load BUSCO/5.5.0-foss-2022a

#get assembly statistics using assemblathon
#wget https://raw.githubusercontent.com/KorfLab/Assemblathon/master/assemblathon_stats.pl
#wget http://korflab.ucdavis.edu/Unix_and_Perl/FAlite.pm
#chmod 755 assemblathon_stats.pl

#make output into fasta file format -- do interactively
#awk '/^S/{print ">"$2;print $3}' $OUTDIR/HiFiASM_OmniC/P_citronellum_BF73_OmniC.asm.hic.p_ctg.gfa > $OUTDIR/HiFiASM_OmniC/P_citronellum_OmniC.asm.hic.p_ctg.fa
#awk '/^S/{print ">"$2;print $3}' $OUTDIR/HiFiASM_OmniC/P_citronellum_BF73_OmniC.asm.hic.hap1.p_ctg.gfa > $OUTDIR/HiFiASM_OmniC/P_citronellum_OmniC.hic.hap1.p_ctg.fa
#awk '/^S/{print ">"$2;print $3}' $OUTDIR/HiFiASM_OmniC/P_citronellum_BF73_OmniC.asm.hic.hap2.p_ctg.gfa > $OUTDIR/HiFiASM_OmniC/P_citronellum_OmniC.hic.hap2.p_ctg.fa
#export PERL5LIB=/scratch/srb67793/Pelargonium/scripts:$PERL5LIB
#./assemblathon_stats.pl $OUTDIR/HiFiASM_OmniC/P_citronellum_OmniC.asm.hic.p_ctg.fa > $OUTDIR/HiFiASM/P_citronellum_OmniC.asm.hic.p_ctg.stats.txt
#./assemblathon_stats.pl $OUTDIR/HiFiASM_OmniC/P_citronellum_OmniC.asm.hic.hap1.p_ctg.fa > $OUTDIR/HiFiASM/P_citronellum_OmniC.asm.hic.hap1.p_ctg.stats.txt
#./assemblathon_stats.pl $OUTDIR/HiFiASM_OmniC/P_citronellum_OmniC.asm.hic.hap2.p_ctg.fa > $OUTDIR/HiFiASM/P_citronellum_OmniC.asm.hic.hap2.p_ctg.stats.txt

#estimate assembly completeness using BUSCO
busco  \
    -f -i /scratch/srb67793/Pelargonium/Results/HiFiASM_OmniC/P_citronellum_OmniC.hic.bp.hap1.p_ctg.fa \
    -m genome \
    -o P_citronellum_omnic-hap1-busco \
    -l eudicots_odb10 \
    -c 20

  busco  \
    -f -i /scratch/srb67793/Pelargonium/Results/HiFiASM_OmniC/P_citronellum_OmniC.hic.bp.hap1.p_ctg.fa \
    -m genome \
    -o P_citronellum_omnic-hap2-busco \
    -l eudicots_odb10 \
    -c 20

