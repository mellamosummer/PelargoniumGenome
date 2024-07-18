#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=srb67793@uga.edu
#SBATCH --output=HiCMaps.%j.out
#SBATCH --error=HiCMaps.%j.err
#SBATCH --job-name=HiCMaps_72024
#SBATCH --partition=batch
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=300GB

#ml SAMtools
#ml YaHS
ml Juicer

#Index the scaffold files
#samtools faidx /scratch/srb67793/Pelargonium/Results/P_citronellum-hap1_scaffolds_final.fa
#samtools faidx /scratch/srb67793/Pelargonium/Results/P_citronellum-hap2_scaffolds_final.fa

#Create Chromosome size files
#cut -f1-2 /scratch/srb67793/Pelargonium/Results/P_citronellum-hap1_scaffolds_final.fa.fai > /scratch/srb67793/Pelargonium/Results/P_citronellum-hap1_scaffolds_final.chrom.sizes
#cut -f1-2 /scratch/srb67793/Pelargonium/Results/P_citronellum-hap2_scaffolds_final.fa.fai > /scratch/srb67793/Pelargonium/Results/P_citronellum-hap2_scaffolds_final.chrom.sizes

#Run the YaHS "juicer pre" command
#juicer pre -a -o /scratch/srb67793/Pelargonium/Results/JBAT/P_citronellum-hap1-JBAT /scratch/srb67793/Pelargonium/Results/JBAT/P_citronellum-hap1.bin /scratch/srb67793/Pelargonium/Results/JBAT/P_citronellum-hap1_scaffolds_final.agp /scratch/srb67793/Pelargonium/Results/HiFiASM_OmniC/P_citronellum_OmniC.hic.bp.hap1.p_ctg.fa.fai > /scratch/srb67793/Pelargonium/Results/JBAT/P_citronellum-hap1-JBAT.log 2>&1
#juicer pre -a -o /scratch/srb67793/Pelargonium/Results/JBAT/P_citronellum-hap2-JBAT /scratch/srb67793/Pelargonium/Results/JBAT/P_citronellum-hap2.bin /scratch/srb67793/Pelargonium/Results/JBAT/P_citronellum-hap2_scaffolds_final.agp /scratch/srb67793/Pelargonium/Results/HiFiASM_OmniC/P_citronellum_OmniC.hic.bp.hap2.p_ctg.fa.fai > /scratch/srb67793/Pelargonium/Results/JBAT/P_citronellum-hap2-JBAT.log 2>&1

#Run Juicer to create the Hi-C files
java -jar $EBROOTJUICER/juicer_tools.jar pre /scratch/srb67793/Pelargonium/Results/JBAT/P_citronellum-hap1-JBAT.txt /scratch/srb67793/Pelargonium/Results/JBAT/P_citronellum-hap1-JBAT.hic <(cat /scratch/srb67793/Pelargonium/Results/JBAT/P_citronellum-hap1-JBAT.log | grep PRE_C_SIZE | awk '{print $2" "$3}')
java -jar $EBROOTJUICER/juicer_tools.jar pre /scratch/srb67793/Pelargonium/Results/JBAT/P_citronellum-hap2-JBAT.txt /scratch/srb67793/Pelargonium/Results/JBAT/P_citronellum-hap2-JBAT.hic <(cat /scratch/srb67793/Pelargonium/Results/JBAT/P_citronellum-hap2-JBAT.log | grep PRE_C_SIZE | awk '{print $2" "$3}')

#Download the following files to a Windows or Mac computer so that they can be loaded into Juicebox 
#Export from Juicebox the updated "assembly" file and "liftover" file and name them (for each haplotype)
