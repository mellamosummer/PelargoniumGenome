#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=srb67793@uga.edu
#SBATCH --output=Contam.%j.out
#SBATCH --error=Contam.%j.err
#SBATCH --job-name=Contam_72024
#SBATCH --partition=batch
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=180GB


#load module
ml ncbi-FCS/0.5.0-GCCcore-11.3.0

#contamination report FCS-GX
OUTDIR="/scratch/srb67793/Pelargonium/Results"
python $EBROOTNCBIMINFCS/fcs.py screen genome --fasta $OUTDIR/HiFiASM_OmniC/P_citronellum_OmniC.hic.bp.hap1.p_ctg.fa --out-dir . --gx-db "$GXDB_LOC/gxdb" --tax-id 73188
python $EBROOTNCBIMINFCS/fcs.py screen genome --fasta $OUTDIR/HiFiASM_OmniC/P_citronellum_OmniC.hic.bp.hap2.p_ctg.fa --out-dir . --gx-db "$GXDB_LOC/gxdb" --tax-id 73188

#cat $GENOMEDIR/G_maculatum_BF73.asm.bp.hap1.p_ctg.fa | python $EBROOTNCBIMINFCS/fcs.py clean genome --action-report $OUTDIR/FCS/G_maculatum_BF73.asm.bp.hap1.p$
#cat $GENOMEDIR/G_maculatum_BF73.asm.bp.hap2.p_ctg.fa | python $EBROOTNCBIMINFCS/fcs.py clean genome --action-report $OUTDIR/FCS/G_maculatum_BF73.asm.bp.hap2.p$



