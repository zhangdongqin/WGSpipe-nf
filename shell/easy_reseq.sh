#!/bin/bash
raw_reads_dir=""
genome_fasta=""
index=""
status=off                
output=""
thread=""
indel1=""
indel2=""
dbsnp=""

Usage () {                 
    echo "Usage:"
    echo -e "This program is provided by Zdq ,you can communicate any problems with the software through zhangdongqin2@126.com\n "
    echo -e "raw reads must conform to the format of SAMPLE_1.fastq.gz , SAMPLE_2.fastq.gz\n "
    echo -e "easy_rnaseq.py [-h] [-v] [-d </path/to/raw/reads>] [-o </path/to/output>] [-r </path/to/genome fasta>] [-i </path/to/genome index>] [-t <thread number>][-k <known indel annotation>] [-g <1000 indel annotation>] [-s <dbsnp anntation>]"
    echo -e "To run the program, you need to specify the absolute path of fastp, bwa, samtools, gatk4 software in the configure file in the current directory\n"
    echo -e "\t-h help information for porgram"
    echo -e "\t-d directory of the wgs raw reads"
    echo -e "\t-o result output directory of program"
    echo -e "\t-g dbsnp annotation vcf file"
    echo -e "\t-r genome fasta file"
    echo -e "\t-i bwa index name"
    echo -e "\t-t thread numbers of program"
    echo -e "\t-k indel annotation vcf files"
    echo -e "\t-g indel annotation vcf files"
    exit -1
}

while getopts :hvd:o:r:i:t:k:g:s: varname   

do
   case $varname in
    h)
      #echo "$varname"
      Usage
      exit
      ;;
    v)
      #echo "$varname"
      status=on
      #echo "$status"
      exit
      ;;
    d)
      #echo "$varname"
      #echo "$OPTARG"
      raw_reads_dir=$OPTARG                   
      if [ ! -d $raw_reads_dir ];then             
         echo "the source directory $raw_reads_dir not exist!"
         exit
      fi
      ;;
    o)
      #echo "$varname"
      #echo "$OPTARG"
      output=$OPTARG                      
      if [ ! -d  $output ];then               
         echo "the output path $output not exist"
         exit
      fi
      ;;
    r)
      #echo "$varname"
      #echo "$OPTARG"
      genome_fasta=$OPTARG                      
      if [ ! -f  $genome_fasta ];then              
         echo "the genome fasta file $genome_fasta not exist"
         exit
      fi
      ;;
    g)
      #echo "$varname"
      #echo "$OPTARG"
      indel2=$OPTARG                     
      if [ ! -f  $indel2 ];then             
         echo "the indel file $indel2 does not exist"
         exit
      fi
      ;;
    k)
      #echo "$varname"
      #echo "$OPTARG"
      indel1=$OPTARG                     
      if [ ! -f  $indel1 ];then             
         echo "the indel file $indel1 does not exist"
         exit
      fi
      ;;
    s)
      #echo "$varname"
      #echo "$OPTARG"
      dbsnp=$OPTARG                     
      if [ ! -f  $dbsnp ];then             
         echo "the dbsnp file $dbsnp does not exist"
         exit
      fi
      ;;
    i)
      #echo "$varname"
      #echo "$OPTARG"
      index=$OPTARG                      
      if [ ! -f  ${index}.bwt ];then              
         echo "the index file $index not exist"
         exit
      elif [ ! -f  ${index}.pac ];then 
      	 echo "the index file $index not exist"
         exit
      elif [ ! -f  ${index}.ann ];then 
      	 echo "the index file $index not exist"
         exit
      elif [ ! -f  ${index}.amb ];then 
      	 echo "the index file $index not exist"
         exit
      elif [ ! -f  ${index}.sa ];then 
      	 echo "the index file $index not exist"
         exit
      fi
      ;;
    t)
      #echo "$varname"
      #echo "$OPTARG"
      thread=$OPTARG                    
      ;;  
    :)                                              
      echo "$varname"

      echo "the option -$OPTARG require an arguement"        
      exit 1
      ;;
    ?)                                   
      echo "$varname"
      echo "Invaild option: -$OPTARG"          
      Usage
      exit 2
      ;;
    esac
done

source ./configure

echo "start fastp running"
cd $raw_reads_dir
  if [ ! -d ../clean_reads ];then
    mkdir -p ../clean_reads
  fi 
ls *gz|cut -d"_" -f 1 |sort -u |while read id;do
  if [ ! -f  ${id}_1.fastq.gz ];then               
         echo "The ${id}.1.fq.gz file  does not exist "
       exit
    elif [ ! -f  ${id}_2.fastq.gz ];then             
         echo "The ${id}.2.fq.gz file does not exist "
         exit
    fi    
$FASTP -i ${id}_1.fastq.gz -o ../clean_reads/${id}.1.clean.fq.gz -I ${id}_2.fastq.gz -O ../clean_reads/${id}.2.clean.fq.gz
   if [ ! -f  ../clean_reads/${id}.1.clean.fq.gz ];then              
         echo "The fastp  does not executed correctly "
       exit
   elif [ ! -f  ../clean_reads/${id}.2.clean.fq.gz ];then             
         echo "The fastp  does not executed correctly "
         exit
   fi    
done

echo "start bwa mapping"
cd ../clean_reads

ls *gz|cut -d"." -f 1 |sort -u |while read id;do
  if [ ! -f  ${id}.1.clean.fq.gz ];then               
         echo "The ${id}.1.clean.fq.gz file  does not exist "
       exit
  elif [ ! -f  ${id}.2.clean.fq.gz ];then              
         echo "The ${id}.2.clean.fq.gz file does not exist "
         exit
  fi    
$BWA mem -t $thread -M -R "@RG\tID:lane1\tPL:illumina\tLB:library\tSM:$id" $index ${id}.1.clean.fq.gz ${id}.2.clean.fq.gz > $output/${id}.sam
$SAMTOOLS view -b -S $output/${id}.sam > $output/${id}.bam
$SAMTOOLS sort $output/${id}.bam -o $output/${id}.sorted.bam
$SAMTOOLS flagstat $output/${id}.sorted.bam > $output/${id}.sorted.flagstat  
done

cd $output
  if [ ! -d $output/gatk_bam ];then
    mkdir -p $output/gatk_bam
  fi 
  if [ ! -d $output/gatk_vcf ];then
    mkdir -p $output/gatk_vcf
  fi 

echo "start GATK MarkDuplicates"
cd $output
ls *sorted.bam|cut -d"." -f 1 |sort -u |while read id;do
    if [ ! -f  ${id}.sorted.bam ];then               
         echo "The ${id}.sorted.bam file  does not exist "
       exit
    fi
$GATK --java-options "-Xmx10G -Djava.io.tmpdir=./" MarkDuplicates \
-I ${id}.sorted.bam -O $output/gatk_bam/${id}.sorted.MarkDuplicates.bam \
-M $output/gatk_bam/${id}.sorted.bam.metrics > log.mark
$SAMTOOLS index $output/gatk_bam/${id}.sorted.MarkDuplicates.bam
    if [ ! -f  $output/gatk_bam/${id}.sorted.MarkDuplicates.bam ];then               
         echo "The GATK  MarkDuplicates program  for ${id}.sorted.bam does not executed correctly"
       exit
    fi   
done

echo "start GATK BaseRecalibrator"
cd $output/gatk_bam
ls *sorted.MarkDuplicates.bam|cut -d"." -f 1 |sort -u |while read id;do
    if [ ! -f  ${id}.sorted.MarkDuplicates.bam ];then               
         echo "The ${id}.sorted.MarkDuplicates.bam file  does not exist "
       exit
    fi
$GATK --java-options "-Xmx10G -Djava.io.tmpdir=./" BaseRecalibrator \
-R $genome_fasta -I $output/gatk_bam/${id}.sorted.MarkDuplicates.bam \
--known-sites $indel1 \
--known-sites $indel2 \
--known-sites $dbsnp \
-O $output/gatk_bam/${id}.recal.data.table
    if [ ! -f  ${id}.recal.data.table ];then               
         echo "The GATK  BaseRecalibrator program  for ${id}.sorted.MarkDuplicates.bam does not executed correctly "
       exit
    fi
done

echo "start GATK ApplyBQSR"
cd  $output/gatk_bam
ls *sorted.MarkDuplicates.bam|cut -d"." -f 1 |sort -u |while read id;do
    if [ ! -f  ${id}.sorted.MarkDuplicates.bam ];then               
         echo "The ${id}.sorted.MarkDuplicates.bam file  does not exist "
       exit
    fi
$GATK --java-options "-Xmx10G -Djava.io.tmpdir=./" ApplyBQSR \
-R $genome_fasta -I $output/gatk_bam/${id}.sorted.MarkDuplicates.bam \
-bqsr $output/gatk_bam/${id}.recal.data.table \
-O $output/gatk_bam/${id}.sorted.MarkDuplicates.BQSR.bam
    if [ ! -f  ${id}.sorted.MarkDuplicates.BQSR.bam ];then               
         echo "The GATK  ApplyBQSR program  for ${id}.sorted.MarkDuplicates.BQSR.bam does not executed correctly "
       exit
    fi   
done

echo "start GATK HaplotypeCaller"
cd $output/gatk_bam
ls *sorted.MarkDuplicates.BQSR.bam|cut -d"." -f 1 |sort -u |while read id;do
    if [ ! -f ${id}.sorted.MarkDuplicates.BQSR.bam ];then               
         echo "${id}.sorted.MarkDuplicates.BQSR.bam file  does not exist "
       exit
    fi
$GATK --java-options "-Xmx8G -Djava.io.tmpdir=./" HaplotypeCaller -R $genome_fasta \
--emit-ref-confidence GVCF -I $output/gatk_bam/${id}.sorted.MarkDuplicates.BQSR.bam \
-D $dbsnp -O $output/gatk_vcf/${id}.gvcf
    if [ ! -f  $output/gatk_vcf/${id}.gvcf ];then               
         echo "The GATK  HaplotypeCaller program  for ${id}.sorted.MarkDuplicates.BQSR.bam does not executed correctly "
       exit
    fi  
done

echo "start GATK GenotypeGVCFs"
cd $output/gatk_vcf
ls *gvcf|cut -d"." -f 1 |sort -u |while read id;do
    if [ ! -f ${id}.gvcf ];then               
         echo "${id}.gvcf file  does not exist "
       exit
    fi
$GATK --java-options "-Xmx8G -Djava.io.tmpdir=./" GenotypeGVCFs \
-R $genome_fasta -V $output/gatk_vcf/${id}.gvcf \
-O $output/gatk_vcf/${id}.raw.vcf
    if [ ! -f  $output/gatk_vcf/${id}.raw.vcf ];then               
         echo "The GATK  GenotypeGVCFs program  for $output/gatk_vcf/${id}.gvcf does not executed correctly "
       exit
    fi  
done
