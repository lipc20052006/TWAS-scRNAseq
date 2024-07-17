
fastqc -t 50 -o ./ ../0.rawdata.re/*fastq.gz &
multiqc ./*_fastqc.zip
ls ./0.rawdata.re/*fq.gz >sample.txt
cat sample.txt |while read line; 
do 
trim_galore --phred33 -q 20 --stringency 3 --gzip --length 40 -o 2.trim/ -j 8 --paired 0.rawdata.re/${line}_r1.fq.gz 0.rawdata.re/${line}_r2.fq.gz 
done

index=./zmV5/bowtie2/zmV5
#直接运行下面脚本
cat sample.txt |while read line
do
bowtie2 -x $index -p 50 -1 2.trim/${line}_r1_val_1.fq.gz -2 2.trim/${line}_r2_val_2.fq.gz | samtools view -b - -q 30 | samtools sort -@ 50 - -o 3.bam/${line}_sorted.q30.bam 
samtools flagstat 3.bam/${line}_sorted.q30.bam  > 3.bam/${line}_sorted.q30.bam.stat
sambamba markdup --overflow-list-size 600000  --tmpdir='./' -r 3.bam/${line}_sorted.q30.bam  3.bam/${line}_rmdup.bam
done

ls *q30.bam.stat|while read id;
do 
echo $id >>q30.stat;
cat $id |grep -E 'total|N/A' |grep -v singletons >>q30.stat;
done
ls *_sorted.bam|while read id;do (samtools flagstat $id  > $id.stat &) ;done
ls *_rmdup.bam|while read id;do (samtools flagstat $id  > $id.stat &) ;done
ls *q30.bam|while read id;do (samtools flagstat $id  > $id.stat &) ;done
ls *rmdup.bam | while read id; do 
bamCoverage -p 10 --binSize 10 --normalizeUsing RPKM -b $id -o ${id%_*}.rpkm.bw &
done
macs2 callpeak -c bZIP-In_rmdup.bam -t bZIP-1_rmdup.bam  -f BAMPE -B -g 2e9 -n bZIP-1 --outdir ../4.peak  2> ../4.peak/bZIP-1.log &
macs2 callpeak -c bZIP-In_rmdup.bam -t bZIP-2_rmdup.bam  -f BAMPE -B -g 2e9 -n bZIP-2 --outdir ../4.peak  2> ../4.peak/bZIP-2.log &
# peak anno
for i in bZIP-1 bZIP-2
do
echo $i
Rscript ../6.chipseek.R "$i"
done
