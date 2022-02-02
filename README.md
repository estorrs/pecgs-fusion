# pecgs-fusion
cwl wrapper for dinglab fusion pipeline

based on [dinglab fusion pipeline](https://github.com/ding-lab/Fusion_hg38)

#### Running fusion.py

###### Arguments

usage: fusion.py [-h] [--cpu CPU] [--genome-lib-dir GENOME_LIB_DIR]
                 [--genome-db GENOME_DB] [--bwts BWTS]
                 [--integrate-fasta INTEGRATE_FASTA]
                 [--integrate-annotations INTEGRATE_ANNOTATIONS]
                 [--combine-call-script COMBINE_CALL_SCRIPT]
                 [--filter-script FILTER_SCRIPT]
                 sample fq_1 fq_2

positional arguments:

  - sample
    - Sample id
  - fq_1
    - RNA-seq fastq 1
  - fq_2
    - RNA-seq fastq 2
  
  
optional arguments:

  + -h, --help
    - show this help message and exit
  + --cpu CPU
    - Number of CPUs to run on
  + --genome-lib-dir GENOME_LIB_DIR
    - GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir directory location
  + --genome-db GENOME_DB
    - ericscript_db_homosapiens_ensembl84 directory location
  + --bwts BWTS
    - bwts location
  + --integrate-fasta INTEGRATE_FASTA
    - fasta used during integration step
  + --integrate-annotations INTEGRATE_ANNOTATIONS
    - location of annot.ensembl.GRCh38.txt used during integration step
  + --combine-call-script COMBINE_CALL_SCRIPT
    - location of combine_call.pl
  + --filter-script FILTER_SCRIPT
    - location of filter.pl


###### Docker

docker image is available from dockerhub (estorrs/pecgs_fusion:0.0.1)
