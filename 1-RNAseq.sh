analysis_dir=$1
config_file=$2
number1=$3
number2=$4
species=$5

cd $analysis_dir 

mkdir -p 0.logs 0.status 

mkdir -p 1.qc
mkdir -p 2.clean_fq
mkdir -p 3.hisat2
mkdir -p 4.featureCounts
mkdir -p 5.starfusion
 
gtf='./maize/ref/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.52.gtf'
hisat_index="./maize/ref/Zea_mays.v5.hisat2_index"
# ls $hisat_index*  
		 
cat $config_file  |while read id;do 
	if((i%${number1}==${number2}))
	then
      arr=(${id})
	  fq1=${arr[1]}
	  fq2=${arr[2]}
      sample=${arr[0]}
	  # raw fastq  : $fq1 $fq2  
	  # clean fastq : ${sample}_1_val_1.fq.gz -2 2.clean_fq/${sample}_2_val_2.fq.gz
	 
	 
##################################################
#####################  fastqc #################### 
##################################################

        biosoft=fastqc
        if [  ! -f   ./0.status/ok.${sample}_${biosoft}.status  ]; ##如果不存在这个文件，就跑下面的代码
 		then
 		 start=$(date +%s.%N)
 		 echo start ${biosoft}  ${sample} `date` >>  ./0.logs/${sample}.log
           ${biosoft} $fq1 $fq2 -o 1.qc 1>./0.logs/logs.${sample}_${biosoft} 2>&1 
            if [ $? -eq 0 ]
			then
				echo  `date` "${biosoft} succeed ${sample}"  >> ./0.logs/${sample}.log
				touch ./0.status/ok.${sample}_${biosoft}.status 
			else
				echo `date` "${biosoft} failed ${sample}"   >> ./0.logs/${sample}.log
			fi # 
      	 dur=$(echo "$(date +%s.%N) - ${start}" | bc)
		 printf "Execution time for ${biosoft} : %.6f seconds" ${dur} >>  ./0.logs/${sample}.log 
		 echo >>  ./0.logs/${sample}.log 
			
		else
      	echo `date` "${biosoft} skip ${sample}"   >> ./0.logs/${sample}.log
        fi # 


##################################################
#####################  trim_galore ################ 
##################################################


	    biosoft=trim_galore
        if [  ! -f   ./0.status/ok.${sample}_${biosoft}.status  ]; 
 		then
 		 start=$(date +%s.%N)
 		 echo start ${biosoft}  ${sample} `date` >>  ./0.logs/${sample}.log
         ${biosoft} -q 20 --length 36 --max_n 3 --stringency 3 --paired -o 2.clean_fq $fq1 $fq2  1>./0.logs/logs.${sample}_${biosoft} 2>&1     
            if [ $? -eq 0 ]
			then
				echo  `date` "${biosoft} succeed ${sample}"  >> ./0.logs/${sample}.log
				touch ./0.status/ok.${sample}_${biosoft}.status 
			else
				echo `date` "${biosoft} failed ${sample}"   >> ./0.logs/${sample}.log
			fi # 
      	dur=$(echo "$(date +%s.%N) - ${start}" | bc)
		 printf "Execution time for ${biosoft} : %.6f seconds" ${dur} >>  ./0.logs/${sample}.log 
			echo >>  ./0.logs/${sample}.log 
			  else
      	echo `date` "${biosoft} skip ${sample}"   >> ./0.logs/${sample}.log
        fi #  


##################################################
#####################  hisat2     ################ 
##################################################

	  biosoft=hisat2
        if [  ! -f   ./0.status/ok.${sample}_${biosoft}.status  ]; ##
 		then
 		 start=$(date +%s.%N)
		 #  
	     fastqc  2.clean_fq/*${sample}_1_val_1.fq.gz   2.clean_fq/*${sample}_2_val_2.fq.gz -o 1.qc 1>./0.logs/logs.${sample}_${biosoft} 2>&1 
 		 echo start ${biosoft}  ${sample} `date` >>  ./0.logs/${sample}.log
		 index=${hisat_index}
         ${biosoft} -p 2  -x  ${index} -1 2.clean_fq/${sample}_1_val_1.fq.gz -2 2.clean_fq/${sample}_2_val_2.fq.gz | samtools sort -@ 2 -o 3.hisat2/${sample}.sort.bam 1>./0.logs/logs.${sample}_${biosoft} 2>&1     
            if [ $? -eq 0 ]
			then
				echo  `date` "${biosoft} succeed ${sample}"  >> ./0.logs/${sample}.log
				touch ./0.status/ok.${sample}_${biosoft}.status 
			else
				echo `date` "${biosoft} failed ${sample}"   >> ./0.logs/${sample}.log
			fi # 
      	dur=$(echo "$(date +%s.%N) - ${start}" | bc)
		 printf "Execution time for ${biosoft} : %.6f seconds" ${dur} >>  ./0.logs/${sample}.log 
			echo >>  ./0.logs/${sample}.log 
			
			  else
      	echo `date` "${biosoft} skip ${sample}"   >> ./0.logs/${sample}.log
        fi # 


##################################################
#####################  featureCounts  ############
##################################################

	  biosoft=featureCounts
        if [  ! -f   ./0.status/ok.${sample}_${biosoft}.status  ]; 
 		then
 		 start=$(date +%s.%N)
 		 echo start ${biosoft}  ${sample} `date` >>  ./0.logs/${sample}.log

      ${biosoft} -T 2 -p -t exon -g gene_id -a $gtf -o 4.featureCounts/${sample}.txt 3.hisat2/${sample}.sort.bam 1>./0.logs/logs.${sample}_${biosoft} 2>&1     
            if [ $? -eq 0 ]
			then
				echo  `date` "${biosoft} succeed ${sample}"  >> ./0.logs/${sample}.log
				touch ./0.status/ok.${sample}_${biosoft}.status 
			else
				echo `date` "${biosoft} failed ${sample}"   >> ./0.logs/${sample}.log
			fi # 
      	dur=$(echo "$(date +%s.%N) - ${start}" | bc)
		 printf "Execution time for ${biosoft} : %.6f seconds" ${dur} >>  ./0.logs/${sample}.log 
			echo >>  ./0.logs/${sample}.log 
			  else
      	echo `date` "${biosoft} skip ${sample}"   >> ./0.logs/${sample}.log
        fi # 

##################################################
#####################  STAR-Fusion  ############
##################################################
notRun(){  
	    biosoft=STAR-Fusion
        if [  ! -f   ./0.status/ok.${sample}_${biosoft}.status  ]; 
 		then
 		 start=$(date +%s.%N)
 		 echo start ${biosoft}  ${sample} `date` >>  ./0.logs/${sample}.log

          ${biosoft} --genome_lib_dir  $db \
             --left_fq  2.clean_fq/${sample}_1_val_1.fq.gz  \
             --right_fq 2.clean_fq/${sample}_2_val_2.fq.gz  \
             --output_dir 5.starfusion/${sample}_outdir
            if [ $? -eq 0 ]
			then
				echo  `date` "${biosoft} succeed ${sample}"  >> ./0.logs/${sample}.log
				touch ./0.status/ok.${sample}_${biosoft}.status 
			else
				echo `date` "${biosoft} failed ${sample}"   >> ./0.logs/${sample}.log
			fi # 
      	dur=$(echo "$(date +%s.%N) - ${start}" | bc)
		 printf "Execution time for ${biosoft} : %.6f seconds" ${dur} >>  ./0.logs/${sample}.log 
			echo >>  ./0.logs/${sample}.log 
			  else
      	echo `date` "${biosoft} skip ${sample}"   >> ./0.logs/${sample}.log
        fi # 
} 

##################################################
##################### TODO   #####################
##################################################

		
    fi #${number1}==${number2}
	i=$((i+1))
done


