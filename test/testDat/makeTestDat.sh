# Making the short SNP file:
(for CHR in {1..24}; do cat /projects1/public_data/Datashare_Boston_Jena_June2018.backup/1240K.snp | awk -v chr=$CHR '$2==chr' | head -100; done) > 1240k_eigenstrat_snp_short.snp.txt

# Making a pos file for samtools
cat 1240k_eigenstrat_snp_short.snp.txt | awk '{if($2==23)$2="X"; if($2==24)$2="Y"; print $2, $4}' > 1240k_eigenstrat_snp_short.pos.txt

# Making a short example mpileup on these positions (takes 10 minutes):
BAM_DIR=/projects1/users/schiffels/AncientBritish/bams; samtools mpileup -B -q30 -Q30 -R -f /projects1/Reference_Genomes/Human/hs37d5/hs37d5.fa -l 1240k_eigenstrat_snp_short.pos.txt $BAM_DIR/12880A.bam $BAM_DIR/12881A.bam $BAM_DIR/12883A.bam $BAM_DIR/12885A.bam > AncientBritish.short.pileup.txt

# Running pileupCaller:
pileupCaller --sampleNames 12880A,12881A,12883A,12885A --randomHaploid --singleStrandMode -f 1240k_eigenstrat_snp_short.snp.txt < AncientBritish.short.pileup.txt
