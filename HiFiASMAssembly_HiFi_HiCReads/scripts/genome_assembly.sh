#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=srb67793@uga.edu
#SBATCH --output=FCSClean.%j.out
#SBATCH --error=FCSClean.%j.err
#SBATCH --job-name=BUSCO_52024
#SBATCH --partition=batch
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=120GB

#script modified from https://gitlab.com/ficklinlab-public/wa-38-genome/-/tree/main/04-nuclear_assembly?ref_type=heads
#set output directory variable
READSDIR="/scratch/srb67793/Geranium/PacBioReads"
OUTDIR="/scratch/srb67793/Geranium/Results"
HICREADSDIR="/scratch/srb67793/Geranium/OmniCReads"

#load modules
module load FastQC/0.11.9-Java-11
#module load MultiQC/1.14-foss-2022a
#module load hifiasm/0.19.6-GCCcore-11.3.0
#module load BUSCO/5.5.0-foss-2022a
#module load ncbi-FCS/0.5.0-GCCcore-11.3.0
#module load MUMmer/4.0.0rc1-GCCcore-11.3.0

#QC Hi-C Reads
fastqc -t 2 $HICREADSDIR/LAWR_OmniC_NA_NA_GAGACGAT_Geranium_maculatum_BF73-Geranium_maculatum_BF73_OmniC_I1420_L6_R1.fastq $HICREADSDIR/LAWR_OmniC_NA_NA_GAGACGAT_Geranium_maculatum_BF73-Geranium_maculatum_BF73_OmniC_I1420_L6_R2.fastq

#fastp \
#    -w 16 \
#    --disable_quality_filtering \
#    --disable_length_filtering \
#    --disable_trim_poly_g \
#    -i hicread1_R1.fastq.gz \
#    -o hicread1_R1.trimmed.fastq.gz \
#    -I hicread2_R2.fastq.gz \
#    -O hicread2_R2.trimmed.fastq.gz

#fastqc -t 2 trimmedhic1 trimmedhic2 trimmedhic3

#QC PacBio reads
#mkdir -p $OUTDIR/fastqc $OUTDIR/multiqc
#fastqc -t 2 -o $OUTDIR/fastqc $READSDIR/*
#multiqc $OUTDIR/fastqc/* -o $OUTDIR/multiqc/

#KMER plot
#ml Jellyfish/2.3.0-GCC-11.3.0
#R01="$READSDIR/m84139_240224_134350_s2.hifi_reads.bc2080.fastq.gz"
#R02="$READSDIR/m84139_240316_071829_s2.hifi_reads.bc2002.fastq.gz"
#R03="$READSDIR/m84139_240224_154318_s3.hifi_reads.bc2080.fastq.gz"
#jellyfish count <(zcat $R01) <(zcat $R02) <(zcat $R03) -t 64 -m 61 -C -o $OUTDIR/jellyfish/G_maculatum_BF73.61.jf --disk -s 4G
#jellyfish histo $OUTDIR/jellyfish/G_maculatum_BF73.61.jf -o $OUTDIR/jellyfish/G_maculatum_BF73.61.hist


#create the phased assembly at contig level using HiFiASM
#mkdir -p $OUTDIR/HiFiASM_DownsampleTest
#hifiasm -o $OUTDIR/HiFiASM_DownsampleTest/G_maculatum_BF73_downsample.asm -t 48 $READSDIR/*half.fa
#    --h1 hiread1 \
#    --h2 hicread2 \

#make output into fasta file format
#awk '/^S/{print ">"$2;print $3}' $OUTDIR/HiFiASM/G_maculatum_BF73.asm.bp.p_ctg.gfa > $OUTDIR/HiFiASM/G_maculatum_BF73.asm.bp.p_ctg.fa
#awk '/^S/{print ">"$2;print $3}' $OUTDIR/HiFiASM/G_maculatum_BF73.asm.bp.hap1.p_ctg.gfa > $OUTDIR/HiFiASM/G_maculatum_BF73.asm.bp.hap1.p_ctg.fa
#awk '/^S/{print ">"$2;print $3}' $OUTDIR/HiFiASM/G_maculatum_BF73.asm.bp.hap2.p_ctg.gfa > $OUTDIR/HiFiASM/G_maculatum_BF73.asm.bp.hap2.p_ctg.fa

#get assembly statistics using assemblathon
#wget https://raw.githubusercontent.com/KorfLab/Assemblathon/master/assemblathon_stats.pl
#wget http://korflab.ucdavis.edu/Unix_and_Perl/FAlite.pm
#chmod 755 assemblathon_stats.pl
#export PERL5LIB=/scratch/srb67793/Geranium/scripts:$PERL5LIB
#./assemblathon_stats.pl $OUTDIR/HiFiASM/G_maculatum_BF73.asm.bp.p_ctg.fa > $OUTDIR/HiFiASM/G_maculatum_BF73.asm.bp.p_ctg.stats.txt
#./assemblathon_stats.pl $OUTDIR/HiFiASM/G_maculatum_BF73.asm.bp.hap1.p_ctg.fa > $OUTDIR/HiFiASM/G_maculatum_BF73.asm.bp.hap1.p_ctg.stats.txt
#./assemblathon_stats.pl $OUTDIR/HiFiASM/G_maculatum_BF73.asm.bp.hap2.p_ctg.fa > $OUTDIR/HiFiASM/G_maculatum_BF73.asm.bp.hap2.p_ctg.stats.txt
#estimate assembly completeness using BUSCO

#busco  \
#    -f -i $OUTDIR/HiFiASM/G_maculatum_BF73.asm.bp.hap1.p_ctg.fa \
#    -m genome \
#    -o G_maculatum_BF73-hapA-busco \
#    -l eudicots_odb10 \
#    -c 20

#  busco  \
#    -f -i $OUTDIR/HiFiASM/G_maculatum_BF73.asm.bp.hap2.p_ctg.fa \
#    -m genome \
#    -o G_maculatum_BF73-hapB-busco \
#    -l eudicots_odb10 \
#    -c 20

#contamination report FCS-GX
#ml ncbi-FCS/0.5.0-GCCcore-11.3.0
GENOMEDIR="/scratch/srb67793/Geranium/Results/HiFiASM"
#python $EBROOTNCBIMINFCS/fcs.py screen genome --fasta ../Results/HiFiASM/G_maculatum_BF73.asm.bp.hap1.p_ctg.fa --out-dir . --gx-db "$GXDB_LOC/gxdb" --tax-id 1663611
#python $EBROOTNCBIMINFCS/fcs.py screen genome --fasta ../Results/HiFiASM/G_maculatum_BF73.asm.bp.hap2.p_ctg.fa --out-dir . --gx-db "$GXDB_LOC/gxdb" --tax-id 1663611

cat $GENOMEDIR/G_maculatum_BF73.asm.bp.hap1.p_ctg.fa | python $EBROOTNCBIMINFCS/fcs.py clean genome --action-report $OUTDIR/FCS/G_maculatum_BF73.asm.bp.hap1.p_ctg.1663611.fcs_gx_report.txt --output $OUTDIR/FCS/G_maculatum_BF73.asm.hap1.clean.fasta --contam-fasta-out $OUTDIR/FCS/G_maculatum_BF73.asm.hap1.contam.fasta
cat $GENOMEDIR/G_maculatum_BF73.asm.bp.hap2.p_ctg.fa | python $EBROOTNCBIMINFCS/fcs.py clean genome --action-report $OUTDIR/FCS/G_maculatum_BF73.asm.bp.hap2.p_ctg.1663611.fcs_gx_report.txt --output $OUTDIR/FCS/G_maculatum_BF73.asm.hap2.clean.fasta --contam-fasta-out $OUTDIR/FCS/G_maculatum_BF73.asm.hap2.contam.fasta

#./assemblathon_stats.pl $OUTDIR/FCS/G_maculatum_BF73.asm.hap1.clean.fasta > $OUTDIR/FCS/G_maculatum_BF73.asm.hap1.clean.stats
#awk '/^S/{print ">"$2;print $3}' $OUTDIR/HiFiASM_DownsampleTest/G_maculatum_BF73_downsample.asm.bp.p_ctg.gfa > $OUTDIR/HiFiASM_DownsampleTest/G_maculatum_BF73_downsample.asm.bp.p_ctg.fa
#awk '/^S/{print ">"$2;print $3}' $OUTDIR/HiFiASM_DownsampleTest/G_maculatum_BF73_downsample.asm.bp.hap1.p_ctg.gfa > $OUTDIR/HiFiASM_DownsampleTest/G_maculatum_BF73_downsample.asm.bp.hap1.p_ctg.fa
#awk '/^S/{print ">"$2;print $3}' $OUTDIR/HiFiASM_DownsampleTest/G_maculatum_BF73_downsample.asm.bp.hap2.p_ctg.gfa > $OUTDIR/HiFiASM_DownsampleTest/G_maculatum_BF73_downsample.asm.bp.hap2.p_ctg.fa

#./assemblathon_stats.pl $OUTDIR/HiFiASM_DownsampleTest/G_maculatum_BF73_downsample.asm.bp.p_ctg.fa > $OUTDIR/HiFiASM_DownsampleTest/G_maculatum_BF73_downsample.asm.bp.p_ctg.stats
#./assemblathon_stats.pl $OUTDIR/HiFiASM_DownsampleTest/G_maculatum_BF73_downsample.asm.bp.hap1.p_ctg.fa > $OUTDIR/HiFiASM_DownsampleTest/G_maculatum_BF73_downsample.asm.bp.hap1.p_ctg.stats
#./assemblathon_stats.pl $OUTDIR/HiFiASM_DownsampleTest/G_maculatum_BF73_downsample.asm.bp.hap2.p_ctg.fa > $OUTDIR/HiFiASM_DownsampleTest/G_maculatum_BF73_downsample.asm.bp.hap2.p_ctg.stats

 # busco  \
 #   -f -i $OUTDIR/HiFiASM_DownsampleTest/G_maculatum_BF73_downsample.asm.bp.hap1.p_ctg.fa \
 #   -m genome \
 #   -o G_maculatum_BF73-hap1-downsample-busco \
 #   -l eudicots_odb10 \
 #   -c 20

#  busco  \
#    -f -i $OUTDIR/HiFiASM_DownsampleTest/G_maculatum_BF73_downsample.asm.bp.hap2.p_ctg.fa \
#    -m genome \
#    -o G_maculatum_BF73-hap2-downsample-busco \
#    -l eudicots_odb10 \
#    -c 20


 # busco  \
 #   -f -i $OUTDIR/FCS/G_maculatum_BF73.asm.hap1.clean.fasta \
 #   -m genome \
 #   -o G_maculatum_BF73-hap1-cleaned-busco \
 #   -l eudicots_odb10 \
 #   -c 20


#HiC data steps
#Step 1: The BWA aligner requires that the input sequences are indexed.
#bwa index -a bwtsw WA_3e.asm.hic.hap1.p_ctg.fa
#bwa index -a bwtsw WA_38.asm.hic.hap2.p_ctg.fa

#Step 2: Perform the BWA alignment
#bwa mem -t 24 -5SP \
#    G_maculatum_BF73.asm.hic.hap1.p_ctg.fa \
#    JNQM_OmniC_NA_NA_CTTGTCGA_Apple_Cosmic_Crisp-Apple_WA38_Cosmic_OmniC_I1161L1_L1_R1.fastq \
#    JNQM_OmniC_NA_NA_CTTGTCGA_Apple_Cosmic_Crisp-Apple_WA38_Cosmic_OmniC_I1161L1_L1_R2.fastq > G_maculatum_BF73.asm.hic.hap1.p_ctg_vs_HiC.sam

#bwa mem -t 24 -5SP \
#    G_maculatum_BF73.asm.hic.hap2.p_ctg.fa \
#    JNQM_OmniC_NA_NA_CTTGTCGA_Apple_Cosmic_Crisp-Apple_WA38_Cosmic_OmniC_I1161L1_L1_R1.fastq \
#    JNQM_OmniC_NA_NA_CTTGTCGA_Apple_Cosmic_Crisp-Apple_WA38_Cosmic_OmniC_I1161L1_L1_R2.fastq > G_maculatum_BF73.asm.hic.hap2.p_ctg_vs_HiC.sam

#Step 3: The BWA aligner creates an alignment file in SAM format. Next, convert the SAM alignment to a BAM file
#samtools view -h -b -F 2316 G_maculatum_BF73.asm.hic.hap1.p_ctg_vs_HiC.sam > G_maculatum_BF73.asm.hic.hap1.p_ctg_vs_HiC.bam
#samtools view -h -b -F 2316 G_maculatum_BF73.asm.hic.hap2.p_ctg_vs_HiC.sam > G_maculatum_BF73.asm.hic.hap2.p_ctg_vs_HiC.bam

#Step 4: the BAM file can be sorted which makes lookups easier for other applications.
#samtools sort -n G_maculatum_BF73.asm.hic.hap1.p_ctg_vs_HiC.bam > G_maculatum_BF73.asm.hic.hap1.p_ctg_vs_HiC.sorted.bam
#samtools sort -n G_maculatum_BF73.asm.hic.hap2.p_ctg_vs_HiC.bam > G_maculatum_BF73.asm.hic.hap2.p_ctg_vs_HiC.sorted.bam

#HiC
#hic_qc.py -b G_maculatum_BF73.asm.hic.hap1.p_ctg_vs_HiC.sorted.bam -o G_maculatum_BF73.hap1_vs_HiC
#hic_qc.py -b G_maculatum_BF73.asm.hic.hap2.p_ctg_vs_HiC.sorted.bam -o G_maculatum_BF73.hap2_vs_HiC

#generate mummer plot
#GENOMEDIR="/scratch/srb67793/Geranium/Results/HiFiASM"
#mkdir $OUTDIR/Mummer
#nucmer $GENOMEDIR/G_maculatum_BF73.asm.bp.hap1.p_ctg.fa $GENOMEDIR/G_maculatum_BF73.asm.bp.hap2.p_ctg.fa -p $OUTDIR/Mummer/nucmer
#delta-filter -1 $OUTDIR/Mummer/nucmer.delta > $OUTDIR/Mummer/nucmer.1delta
#mummerplot --size large -layout --color -f --png $OUTDIR/Mummer/nucmer.1delta -p $OUTDIR/Mummer/nucmer

#index the contigs from the hifiasm step
#samtools faidx G_maculatum_BF73.asm.hic.hap1.p_ctg.fa 
#samtools faidx G_maculatum_BF73.asm.hic.hap2.p_ctg.fa 

# run YaHS for scaffolding
#yahs -o G_maculatum_BF73-hap1 G_maculatum_BF73.asm.hic.hap1.p_ctg.fa G_maculatum_BF73.asm.hic.hap1.p_ctg_vs_HiC.sorted.bam
#yahs -o G_maculatum_BF73-hap2 G_maculatum_BF73.asm.hic.hap2.p_ctg.fa G_maculatum_BF73.asm.hic.hap2.p_ctg_vs_HiC.sorted.bam

#run the assemblathon stats script
#./assemblathon_stats.pl G_maculatum_BF73-hap1_scaffolds_final.fa > G_maculatum_BF73-hap1_scaffolds_final.stats.txt
#./assemblathon_stats.pl G_maculatum_BF73-hap2_scaffolds_final.fa > G_maculatum_BF73-hap2_scaffolds_final.stats.txt

#Index the scaffold files
#samtools faidx G_maculatum_BF73-hap1_scaffolds_final.fa
#samtools faidx G_maculatum_BF73-hap2_scaffolds_final.fa

#Create Chromosome size files
#cut -f1-2 G_maculatum_BF73-hap1_scaffolds_final.fa.fai > G_maculatum_BF73-hap1_scaffolds_final.chrom.sizes
#cut -f1-2 G_maculatum_BF73-hap2_scaffolds_final.fa.fai > G_maculatum_BF73-hap2_scaffolds_final.chrom.sizes

#Run the YaHS "juicer pre" command
#juicer pre \
#    -a \
#    -o G_maculatum_BF73-hap1-JBAT \
#    G_maculatum_BF73-hap1.bin \
#    G_maculatum_BF73-hap1_scaffolds_final.agp \
#    G_maculatum_BF73.asm.hic.hap1.p_ctg.fa.fai > G_maculatum_BF73-hap1-JBAT.log 2>&1

#juicer pre \
#    -a \
#    -o G_maculatum_BF73-hap2-JBAT \
#    G_maculatum_BF73-hap2.bin \
#    G_maculatum_BF73-hap2_scaffolds_final.agp \
#    G_maculatum_BF73.asm.hic.hap2.p_ctg.fa.fai > G_maculatum_BF73-hap2-JBAT.log 2>&1

#Run Juicer to create the Hi-C files
#java -Xmx96G -jar /usr/local/bin/juicer/juicer_tools_1.22.01.jar pre \
#      G_maculatum_BF73-hap1-JBAT.txt G_maculatum_BF73-hap1-JBAT.hic <(cat G_maculatum_BF73-hap1-JBAT.log | grep PRE_C_SIZE | awk '{print $2" "$3}')

#java -Xmx96G -jar /usr/local/bin/juicer/juicer_tools_1.22.01.jar pre \
#      G_maculatum_BF73-hap2-JBAT.txt G_maculatum_BF73-hap2-JBAT.hic <(cat G_maculatum_BF73-hap2-JBAT.log | grep PRE_C_SIZE | awk '{print $2" "$3}')

#Download the following files to a Windows or Mac computer so that they can be loaded into Juicebox 
#Export from Juicebox the updated "assembly" file and "liftover" file and name them (for each haplotype)

#generate an updated FASTA file of our assembl using the YaHS juicer post
#juicer post \
#    -o G_maculatum_BF73-hap1-JBAT \
#   G_maculatum_BF73-hap1-JBAT.review.assembly \
#   G_maculatum_BF73-hap1-JBAT.liftover.agp \
#   G_maculatum_BF73.asm.hic.hap1.p_ctg.fa

#juicer post \
#    -o G_maculatum_BF73-hap2-JBAT \
#   G_maculatum_BF73-hap2-JBAT.review.assembly \
#   G_maculatum_BF73-hap2-JBAT.liftover.agp \
#   G_maculatum_BF73.asm.hic.hap2.p_ctg.fa

#run BUSCO on final assembly
#mkdir busco_input_files
#cd busco_input_files
#mv G_maculatum_BF73-hap1-JBAT.FINAL.fa G_maculatum_BF73-hap2-JBAT.FINAL.fa .  #FINAL files will be after juicebox curation
#cd ..

#busco  \
#    -i busco_input_files \
#    -m genome \
#    -o G_maculatum_BF73-busco \
#    -l eudicots_odb10 \
#    -c 20

#other steps i can try later

#09-chromosome-rename
#Chromosomes were renamed according to that of 'Gala' hapA assembly.


#10-chromosome-reorient
#Chromosomes were reoriented to match 'Gala' hapA assembly.


#11-filtering
#Contaminants (i.e. plastids, bacterial, and virus) were identified and removed. A BUSCO analysis was performed on the filtered assemblies.


#12-kmer-phasing
#Meryl was used to identify parentage of the chromosomes. Chromosomes from the same parents were placed in the same haplome. The 2 haplomes are now designated as hapA (with a Honeycrisp origin) and hapB (with a Enterprise origin)


#13-final_assembly_busco_mummer
#Last sanity check as well as structural variantion and BUSCO analysis.
