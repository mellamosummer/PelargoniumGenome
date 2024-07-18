#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=srb67793@uga.edu
#SBATCH --output=Mummer.%j.out
#SBATCH --error=Mummer.%j.err
#SBATCH --job-name=Mummer_72024
#SBATCH --partition=batch
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=100GB

ml MUMmer
#nucmer --maxmatch -l 100 -c 500 /scratch/srb67793/Pelargonium/Results/HiFiASM_OmniC/P_citronellum_OmniC.hic.bp.hap1.p_ctg.fa /scratch/srb67793/Pelargonium/Results/HiFiASM_OmniC/P_citronellum_OmniC.hic.bp.hap2.p_ctg.fa
#dnadiff -d out.delta 
#delta-filter -1 /scratch/srb67793/Pelargonium/out.delta > /scratch/srb67793/Pelargonium/out.2delta
#mummerplot /scratch/srb67793/Pelargonium/out.delta -R /scratch/srb67793/Pelargonium/Results/HiFiASM_OmniC/P_citronellum_OmniC.hic.bp.hap1.p_ctg.fa -Q /scratch/srb67793/Pelargonium/Results/HiFiASM_OmniC/P_citronellum_OmniC.hic.bp.hap2.p_ctg.fa --filter --layout
mummerplot --size large -layout --color -f --png /scratch/srb67793/Pelargonium/out.delta > /scratch/srb67793/Pelargonium/out.2delta -p /scratch/srb67793/Pelargonium
