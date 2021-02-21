#!/usr/bin/env bash

SCRIPTARGS="$@"
echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> START: " $0 $SCRIPTARGS
echo Running on host `hostname`
echo Strat Time is `date`
echo Directory is `pwd`


###################################################################################

GATK="/home/origene/yinshan.cui/bin/gatk-4.1.9.0/gatk"
###
CONVERT2ANNOVAR="perl /home/origene/yinshan.cui/bin/annovar-master/convert2annovar.pl"
###
TABLE_ANNOVAR="perl /home/origene/yinshan.cui/bin/annovar-master/table_annovar.pl"
###
output="/home/data/wgs_out/other_wgs"
###
gatk_ann_dir="/home/data/genome_resequencing_DB"
###
indel1="${gatk_ann_dir}/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz"
###
indel2="${gatk_ann_dir}/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" 
###
dbsnp="${gatk_ann_dir}/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf"
###
hapmap="${gatk_ann_dir}/resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz"
###
phase="${gatk_ann_dir}/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz"
###
omni="${gatk_ann_dir}/resources_broad_hg38_v0_1000G_omni2.5.hg38.vcf.gz"
###
vqsr_indel="${gatk_ann_dir}/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
###
merge_vcf="/home/data/wgs_out/other_wgs/gatk_vcf/ln229MUFDSW202189279.raw.vcf"
###
genome_fasta="/home/origene/yinshan.cui/database/homo_sapiens/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"
###
humandb="/home/origene/yinshan.cui/bin/annovar-master/humandb"
###

####################################################################################
cd $output

	if [ ! -d $output/annovar ];then
	  mkdir -p $output/annovar
	fi 

	$GATK SelectVariants \
			-select-type SNP \
			-V $merge_vcf \
			-O $output/annovar/merged.snp.raw.vcf


	$GATK SelectVariants \
			-select-type INDEL \
			-V $merge_vcf \
			-O $output/annovar/merged.indel.raw.vcf

#####################################################################################

if [ -f  $output/annovar/merged.snp.raw.vcf ];then               

    $GATK --java-options "-Xmx8G -Djava.io.tmpdir=./" \
			VariantRecalibrator \
			-R $genome_fasta \
			-V $output/annovar/merged.snp.raw.vcf \
			-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
			-resource:omini,known=false,training=true,truth=false,prior=12.0 $omni \
			-resource:1000G,known=false,training=true,truth=false,prior=10.0 $phase \
			-resource:dbsnp,known=true,training=false,truth=false,prior=6.0 $dbsnp \
			--tranches-file $output/annovar/merged.snps.tranches \
			--rscript-file $output/annovar/merged.snps.plots.R \
			-an QD \
			-an MQ \
			-an MQRankSum \
			-an ReadPosRankSum \
			-an FS \
			-an SOR \
			-an DP \
			-mode SNP \
			-O $output/annovar/merged.snps.recal 

    $GATK ApplyVQSR \
			-R ${genome_fasta} \
			-V $output/annovar/merged.snp.raw.vcf \
			--tranches-file $output/annovar/merged.snps.tranches \
			-recal-file $output/annovar/merged.snps.recal \
			-mode SNP \
			-O $output/annovar/merged.VQSR.snps.vcf
fi

if [ -f $output/annovar/merged.indel.raw.vcf ];then

    $GATK --java-options "-Xmx8G -Djava.io.tmpdir=./" \
			VariantRecalibrator \
			-R $genome_fasta \
			-V $output/annovar/merged.indel.raw.vcf \
			-resource:indel,known=false,training=true,truth=true,prior=15.0 $vqsr_indel \
			-resource:dbsnp,known=true,training=false,truth=false,prior=6.0 $dbsnp \
			--tranches-file $output/annovar/merged.indel.tranches \
			--rscript-file $output/annovar/merged.indel.plots.R \
			-an QD \
			-an MQ \
			-an MQRankSum \
			-an ReadPosRankSum \
			-an FS \
			-an SOR \
			-an DP \
			-mode INDEL \
			-O $output/annovar/merged.indel.recal 

    $GATK ApplyVQSR \
			-R ${genome_fasta} \
			-V $output/annovar/merged.indel.raw.vcf \
			--tranches-file $output/annovar/merged.indel.tranches \
			-recal-file $output/annovar/merged.indel.recal  \
			-mode INDEL \
			-O $output/annovar/merged.VQSR.indel.vcf

fi 


			if [ ! -f  $output/annovar/merged.VQSR.snps.vcf ];then               
				echo "$output/annovar/merged.VQSR.snps.vcf does not exist!"
				exit
			fi
			
			if [ ! -f $output/annovar/merged.VQSR.indel.vcf ];then
				echo "$output/annovar/merged.VQSR.indel.vcf does not exist!"
				exit
			fi 


			if [ -f  $output/annovar/merged.VQSR.snps.vcf ];then               

				$CONVERT2ANNOVAR  -format vcf4 $output/annovar/merged.VQSR.snps.vcf > $output/annovar/snp.avinput
			fi
				
			if [ -f $output/annovar/merged.VQSR.indel.vcf ];then

				$CONVERT2ANNOVAR  -format vcf4 $output/annovar/merged.VQSR.indel.vcf > $output/annovar/indel.avinput

			fi 

			cd $output/annovar

			if [ ! -f  $output/annovar/snp.avinput ];then               
				echo "$output/annovar/snp.avinput does not exist!"
				exit
			fi
			
			if [ ! -f $output/annovar/indel.avinput ];then
				echo "$output/annovar/indel.avinput does not exist!"
				exit

			fi 
#########################################################################################


			if [ -f  $output/annovar/snp.avinput ];then               
			  
			  $TABLE_ANNOVAR $output/annovar/snp.avinput $humandb -buildver hg38 \
			    -out snp_annotation \
			    -remove -protocol refGene,cytoBand,esp6500siv2_all \
			    -operation g,r,f \
			    -nastring . -csvout
			fi
				
			if [ -f $output/annovar/indel.avinput ];then
			  
			  $TABLE_ANNOVAR $output/annovar/indel.avinput $humandb -buildver hg38 \
			    -out indel_annotation \
			    -remove -protocol refGene,cytoBand,esp6500siv2_all \
			    -operation g,r,f \
			    -nastring . -csvout
			fi 

