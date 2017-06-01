# using bam files from Jing's mapping; /jfw-lab/Projects/cis_trans_regulation/HyLiTE_newTM1/mapping0
srun --ntasks=272 --mem=350G --nodes=1 --time=48:00:00 --partition=knl-medium --pty /bin/bash

module load picard/2.9.0
module load samtools/1.2 # avoids the error [W::bam_merge_core2] No @HD tag found
module load parallel/20160422
module rm java/1.7.0_55
module load gatk/3.6

picard CreateSequenceDictionary R=TM1.nuclNmt.transcript.fasta O=TM1.nuclNmt.transcript.dict
samtools faidx TM1.nuclNmt.transcript.fasta

ls *.bam | parallel 'picard AddOrReplaceReadGroups I={} O={.}.RG.bam RGID={.} RGPL=illumina RGPU=none RGLB=lib1 RGSM={.}NULL &> {.}.log'
ls *.RG.bam | parallel 'samtools sort -@ 8 {} {.}.sort &>> {.}.log'
ls *.sort.bam | parallel 'picard MarkDuplicates I={} O={.}.dedup.bam REMOVE_DUPLICATES=TRUE ASSUME_SORTED=TRUE METRICS_FILE={.}.dedup.stats &>> {.}.log'
ls *.dedup.bam | parallel 'samtools index {} &>> {.}.log'

ls *.dedup.bam | parallel 'picard ReorderSam I={} O={.}.reorder.bam R=TM1.nuclNmt.transcript.fasta CREATE_INDEX=TRUE'
ls *.reorder.bam | parallel 'samtools index {} &>> {.}.log'

### realign intervals
ls *.reorder.bam | parallel --jobs 10 'gatk -T RealignerTargetCreator -R TM1.nuclNmt.transcript.fasta -I {} -o {.}.realign.intervals -nt 27 2> {.}.RTC.err'
ls *.reorder.bam | parallel --jobs 10 'gatk -T IndelRealigner -R TM1.nuclNmt.transcript.fasta -I {} -targetIntervals {.}.realign.intervals -o {.}.realign.bam -nt 27 2> {.}.indel.err'

### call variants and then joint genotyping
ls *.realign.bam | parallel --jobs 10 'gatk -T HaplotypeCaller -R TM1.nuclNmt.transcript.fasta -I {} -nct 27 -stand_emit_conf 10 -stand_call_conf 30 --emitRefConfidence GVCF -variant_index_type LINEAR -variant_index_parameter 128000 --genotyping_mode DISCOVERY -o {.}.raw_variants.g.vcf 2> {.}.HaplotypeCaller.err'

### join the vcf for each chromosome for all accessions
variants=""
for i in *.g.vcf; do variants+="--variant $i "; done
gatk -T GenotypeGVCFs -R TM1.nuclNmt.transcript.fasta -nt 272 -stand_emit_conf 10 -stand_call_conf 30 $variants --out all.TxMx.vcf 2> joint.genotype.err

### generate hard filter for SNPs and indels #runs very fast
gatk -T SelectVariants -R TM1.nuclNmt.transcript.fasta -V all.TxMx.vcf -selectType SNP -o all.TxMx.raw_snp.vcf 2> all.TxMx.SelectVariants_snp.err
gatk -T SelectVariants -R TM1.nuclNmt.transcript.fasta -V all.TxMx.vcf -selectType INDEL -o all.TxMx.raw_indel.vcf 2> all.TxMx.SelectVariants_indel.err

# apply hard filters to snp 
gatk -T VariantFiltration -R TM1.nuclNmt.transcript.fasta -V all.TxMx.raw_snp.vcf --filterExpression "QD<2.0||FS>60.0||MQ<40.0||SOR>4.0||MQRankSum<-12.5||ReadPosRankSum<-8.0" --filterName all_snp_filter -o all.TxMx.filtered_snp.vcf 2> all.TxMx.VariantFiltration_snp.err


#apply hard filters to indel 

gatk -T VariantFiltration -R TM1.nuclNmt.transcript.fasta -V all.raw_indel.vcf --filterExpression "QD<2.0||FS>200.0||SOR>10.0||ReadPosRankSum<-20.0" --filterName TANG_indel_filter -o all.filtered_indel.vcf 2> all.VariantFiltration_indel.err

 
gatk -T SelectVariants -R TM1.nuclNmt.transcript.fasta --variant all.filtered_snp.vcf --excludeFiltered -o all.PASS_snp.vcf 2> all.PASS_snp.err

gatk -T SelectVariants -R TM1.nuclNmt.transcript.fasta --variant all.filtered_indel.vcf --excludeFiltered -o all.PASS_indel.vcf 2> all.PASS_indel.err

## combine -- idk how to combine the snp/indel passed vcf intelligently yet
#gatk -T CombineVariants -V all.PASS_snp.vcf -V all.PASS_indel.vcf  -o Adansonia.snp.indel.vcf -R TM1.nuclNmt.transcript.fasta --#genotypemergeoption uniquify --filteredrecordsmergetype KEEP_IF_ALL_UNFILTERED

# filter sites
#gatk -T SelectVariants -R TM1.nuclNmt.transcript.fasta --variant Adansonia.minimal.vcf -o Adansonia.filtered.vcf --maxNOCALLfraction 0.25 
