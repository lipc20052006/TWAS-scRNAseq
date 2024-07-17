#01clean data -> sortdup bam
samtools faidx REFERENCE
bwa index REFERENCE
gatk CreateSequenceDictionary -R REFERENCE
bwa mem REFERENCE NAME_1.fq.gz NAME_2.fq.gz -t 40 -M | \
biobambam2/bamsormadup threads=20 tmpfile=./tmp_sort inputformat=sam outputformat=bam >NAME.bwa.sortdup.bam
#02变异检测
gatk4 HaplotypeCaller --num-threads 40 -R REFERENCE -I NAME.bwa.sortdup.bam -O NAME.bwa.sortdup.hc4.g.vcf.gz -ERC GVCF --interval-padding 0 -stand-call-conf 30 --pcr-indel-model CONSERVATIVE
#03分区段合并gvcf
gatk CombineGVCFs -R REFERENCE -V NAME.bwa.sortdup.hc4.g.vcf.gz -O ./combine/combine_CHR:START-END_raw.g.vcf.gz -L CHR:START-END
#04genotyping
gatk GenotypeGVCFs -R REFERENCE  --variant ./combine/combine_CHR:START-END_raw.g.vcf.gz -O ./genotype/combine_CHR:START-END_raw.vcf.gz
#05call SNPs
gatk SelectVariants -R REFERENCE -V ./genotype/combine_CHR:START-END_raw.vcf.gz --select-type SNP -O ./snp/combine_CHR:START-END_raw.snp.vcf.gz
gatk VariantFiltration -R REFERENCE -V ./snp/combine_CHR:START-END_raw.snp.vcf.gz --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name 'SNP_filter' -O ./snp/combine_CHR:START-END_hardfilter.snp.vcf.gz
gatk SelectVariants  -R REFERENCE -V ./snp/combine_CHR:START-END_hardfilter.snp.vcf.gz --exclude-filtered  -O ./snp/combine_CHR:START-END_hardfiltered.snp.vcf.gz
#06call INDELs
gatk SelectVariants -R REFERENCE -V ./genotype/combine_CHR:START-END_raw.vcf.gz --select-type INDEL -O ./indel/combine_CHR:START-END_raw.indel.vcf.gz
gatk VariantFiltration -R REFERENCE -V ./indel/combine_CHR:START-END_raw.indel.vcf.gz --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name 'INDEL_filter' -O ./indel/combine_CHR:START-END_hardfilter.indel.vcf.gz
gatk SelectVariants  -R REFERENCE -V ./indel/combine_CHR:START-END_hardfilter.indel.vcf.gz --exclude-filtered  -O ./indel/combine_CHR:START-END_hardfiltered.indel.vcf.gz
#07合并vcf
gatk MergeVcfs -I combine_CHR:START-END_hardfiltered.indel.vcf.gz -O combine_hardfiltered.XXX.vcf.gz