# WGSpipe-nf
 A easy and quick pipeline for WGS/WES data analysis
STEP01:
INSTALL NEXTFLOW FOR PIPELINE: 
		curl -s https://get.nextflow.io | bash  OR INSTALL WITH CONDA : conda install nextflow 
This pipeline is based one nextflow DSL2 , so you can run nextflow self-update to update NEXTFLOW to new version

STEP02:
CONDA ENVIRONMENT DEPLOY:
		conda env create -f environment.yaml
This code will create a conda environment for wgs/wes analysis named wgs-nf
step03:
RUN PIPELINE:
		nextflow run WGSpipe-nf --genome genome.fa --reference /path/to/snpindel_annotation --reads 'reads/*_{1,2}.fq.gz'
