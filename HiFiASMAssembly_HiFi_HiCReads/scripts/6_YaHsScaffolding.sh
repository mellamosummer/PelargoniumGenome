#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=srb67793@uga.edu
#SBATCH --output=YaHS.%j.out
#SBATCH --error=YaHS.%j.err
#SBATCH --job-name=YaHS_72024
#SBATCH --partition=batch
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10GB

ml SAMtools
ml YaHS

#samtools faidx /scratch/srb67793/Pelargonium/Results/HiFiASM_OmniC/P_citronellum_OmniC.hic.bp.hap1.p_ctg.fa
#samtools faidx /scratch/srb67793/Pelargonium/Results/HiFiASM_OmniC/P_citronellum_OmniC.hic.bp.hap2.p_ctg.fa

yahs -o Results/P_citronellum-hap1 /scratch/srb67793/Pelargonium/Results/HiFiASM_OmniC/P_citronellum_OmniC.hic.bp.hap1.p_ctg.fa  /scratch/srb67793/Pelargonium/Results/HiFiASM_OmniC/P_citronellum.asm.hic.hap1.p_ctg_vs_HiC.sorted.bam
yahs -o Results/P_citronellum-hap2 /scratch/srb67793/Pelargonium/Results/HiFiASM_OmniC/P_citronellum_OmniC.hic.bp.hap2.p_ctg.fa /scratch/srb67793/Pelargonium/Results/HiFiASM_OmniC/P_citronellum.asm.hic.hap2.p_ctg_vs_HiC.sorted.bam 


#wget https://raw.githubusercontent.com/KorfLab/Assemblathon/master/assemblathon_stats.pl
#wget http://korflab.ucdavis.edu/Unix_and_Perl/FAlite.pm
#chmod 755 assemblathon_stats.pl
#./assemblathon_stats.pl WA_38-hap1_scaffolds_final.fa > WA_38-hap1_scaffolds_final.stats.txt
#./assemblathon_stats.pl WA_38-hap2_scaffolds_final.fa > WA_38-hap2_scaffolds_final.stats.txt

