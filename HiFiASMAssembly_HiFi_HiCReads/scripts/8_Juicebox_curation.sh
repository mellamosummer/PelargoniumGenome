#generate an updated FASTA file of our assembl using the YaHS juicer post
juicer post \
    -o G_maculatum_BF73-hap1-JBAT \
   G_maculatum_BF73-hap1-JBAT.review.assembly \
   G_maculatum_BF73-hap1-JBAT.liftover.agp \
   G_maculatum_BF73.asm.hic.hap1.p_ctg.fa

juicer post \
    -o G_maculatum_BF73-hap2-JBAT \
   G_maculatum_BF73-hap2-JBAT.review.assembly \
   G_maculatum_BF73-hap2-JBAT.liftover.agp \
   G_maculatum_BF73.asm.hic.hap2.p_ctg.fa
