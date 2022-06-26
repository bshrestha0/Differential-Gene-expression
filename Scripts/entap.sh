#!/bin/bash
#SBATCH --job-name=entapAHA
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=50G
#SBATCH --qos=general
#SBATCH -o entapAHA%j.out
#SBATCH -e entapAHA%j.err
#SBATCH --partition=general
#SBATCH --mail-user=bikash.shrestha@uconn.edu
#SBATCH --mail-type=END

module load anaconda/2.4.0
module load perl/5.24.0
module load diamond/0.9.36
module load interproscan/5.25-64.0  

prot=/core/labs/Wegrzyn/bikash/transcriptome/DESeq2_prep/Vsearch/AHA/AHA.cluster.pep

/core/labs/Wegrzyn/EnTAP/EnTAP_v0.10.4/EnTAP/EnTAP --runP \
--ini entap_config.ini \
-i $prot \
-d /isg/shared/databases/Diamond/Uniprot/uniprot_sprot.dmnd \
-d /isg/shared/databases/Diamond/RefSeq/plant.protein.faa.205.dmnd \
--threads 16
