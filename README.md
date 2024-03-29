# WGSpipe-nf
A easy and quick pipeline for WGS/WES data analysis

## Quick Start
1. -STEP01:
	-INSTALL NEXTFLOW FOR PIPELINE: 
	```bash
	curl -s https://get.nextflow.io | bash  OR INSTALL WITH CONDA : conda install nextflow 
	```
	-This pipeline is based one nextflow DSL2 , so you can run nextflow self-update to update NEXTFLOW to new version

2. -STEP02(optional):
	-CONDA ENVIRONMENT DEPLOY(OPTIONAL): REQUIREMENTS( fastqc, fastp, multiqc, bwa, gatk, annovar, samtools )
	```bash
	conda env create -f environment.yaml
	```
	-This code will create a conda environment for wgs/wes analysis named wgs-nf
3. -STEP03:
	-GIT CLONE PIPELINE
	```bash
	git clone https://github.com/doudoufeishangtian/WGSpipe-nf.git && cd WGSpipe-nf && chmod +x main.nf 
	```
4. -STEP04:
	-RUN PIPELINE:
	```bash
	nextflow run WGSpipe-nf --genome genome.fa --reference /path/to/snpindel_annotation --reads 'reads/*_{1,2}.fq.gz' 
	```