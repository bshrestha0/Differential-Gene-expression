# Differential gene expression using DESeq2 for *Passiflora* and *Adenia* species.

First part: de novo transcriptome assembly using RNA-Seq data, homology based approach to select the gene structure from the assembled transcriptome.

Second part: Differential gene expression using DESeq2 using Rstudio

## Part I: The tutorial will be on *Adenia hastifolia* (AHA) and will need to repeated for *Passiflora morifolia* (PMO)

### Running fastQC for quality assessment on raw read and summarize results with MultiQC

`module load fastqc/0.11.7`
`module load MultiQC/1.10.1`

`mkdir FastQC_raw`
`fastqc -t 12 *fq.gz -o ./FastQC_raw`
`multiqc FastQC_raw`

### Running TrimGalore to remove low quality reads and adapters, if present

`module load TrimGalore/0.6.5`
`module load cutadapt/2.7`
`module load fastqc/0.11.7`

`mkdir trimmed`
`trim_galore --phred33 --fastqc --length 50 -o trimmed --cores 10 \
	--paired AHA_1_NG_R1.fq.gz AHA_1_NG_R2.fq.gz --retain_unpaired`

*We will consider the post-trimmed reads as "clean reads!*

### Re-running fastQC for clean reads and summarize using multiqc

`mkdir FastQC_trimmed`

`fastqc -t 12 ./trimmed/*val_1.fq.gz -o FastQC_trimmed`
`fastqc -t 12 ./trimmed/*val_2.fq.gz -o FastQC_trimmed`

`multiqc FastQC_trimmed`

### De novo transcriptome assembly using Trinity (each RNA-Seq library is assembled separately)

`dir=/path/to/trimmed/reads`

`module load trinity/2.8.5`
`module load samtools`

`Trinity --seqType fq \
	--left "$dir/AHA_1_NG_R1_val_1.fq.gz" --right "$dir/AHA_1_NG_R2_val_2.fq.gz" \
	--min_contig_length 300 --CPU 16 --max_memory 100G --output AHA_1_NG.trinity --full_cleanup`

* This is an example for one set of RNA library, you need to run for other RNA libraries.*

The output of the Trinity assembly will be AHA_1_NG.trinity.Trinity.fasta file. For the multiple runs, multiple ".fasta" files will be produced.

A prefix can be added for each assembly so that when the transcriptomes are combined for ORF prediction, the correct transcript is selected.
`sed 's/>/>AHA_1_NG/g' AHA_1_NG.trinity.Trinity.fasta > AHA_1_NG.trinity.fasta` 

### Aligning clean reads against the assembled transcriptome using HiSAT2 to check how good is the assembly

Index the transcriptome before alignment

`module load hisat2/2.2.0`

`SEQ=/path/to/transcriptome/assembly`

`ls ${SEQ}/*.fasta >> list.txt`

`while read line
do
part=${line%%.*}
hisat2-build -f ${line} ${part}.index
done < list.txt`

Now use the index for alignment. 

`module load hisat2/2.2.0

`trimmeddir=/path/to/trimmed/reads`
`indexDir=/path/to/Index`

`while read line
do
hisat2 -q -x ${indexDir}/${line}.index -1 ${trimmeddir}/${line}_R1_val_1.fq.gz -2 ${trimmeddir}/${line}_R2_val_2.fq.gz -S "$i".sam -p 16
done < input.txt`

I am using "input.txt" as text file with sample ID name to loop each line as below:

AHA_1_NG
AHA_2_G
AHA_3_G
AHA_3_NG
AHA_1_NG
AHA_2a_NG
AHA_2_G
AHA_3_G

### Predict open reading frames (ORFs) for the assembled transcriptome using TransDecoder

Combined transcripts for each condition (G or NG) into single fasta file

`cat AHA_1_NG.Trinity.fasta AHA_3__NG.Trinity.fasta >> AHA_NG.fasta`
`cat AHA_2_G.Trinity.fasta AHA_3_G.Trinity.fasta >> AHA_G.fasta`

*Transdecoder with Pfam and Diamond based homology might take few days to complete!*

`file=AHA_NG.fasta`

#generate longest orfs
`TransDecoder.LongOrfs -t ${file}`

#using diamond blast to find the homologs on the longest orfs
`diamond blastp --query ${file}.transdecoder_dir/longest_orfs.pep  \
    --db /isg/shared/databases/Diamond/Uniprot/uniprot_sprot.dmnd  --max-target-seqs 1 \
    --outfmt 6 --evalue 1e-5 --threads 24 --out diamond.blastp`

#Pfam
`hmmscan --cpu 24 --domtblout pfam.domtblout \
        /isg/shared/databases/Pfam/current/Pfam-A.hmm ${file}.transdecoder_dir/longest_orfs.pep`

#retain homologous proteins only
`TransDecoder.Predict -t ${file} --retain_pfam_hits pfam.domtblout --retain_blastp_hits diamond.blastp --cpu 24`

The output of the TransDecoder will be AHA_G.fasta.transdecoder.cds/pep for gland tissues and AHA_NG.fasta.transdecoder.cds/pep for non-gland tissues

We combined transdecoder output for gland and non-gland tissue into a single file and assess completeness running BUSCO
`cat AHA_G.fasta.transdecoder.cds AHA_NG.fasta.transdecoder.cds > AHA.transdecoder.cds`
`cat AHA_G.fasta.transdecoder.pep AHA_NG.fasta.transdecoder.pep > AHA.transdecoder.pep`

### Running BUSCO to assess completeness

`module load busco/5.0.0`
`module load MetaEuk/4.0`
`module load hmmer/3.2.1`
`module load blast`

`busco -c 12 -m protein -i AHA.transdecoder.pep -o busco_before -l eudicots_odb10`

As combined transcriptome from multiple libraries, we will have redundant transcripts, which we will be removed with clustering

### Cluster at 90% sequence identity using VSEARCH to remove redundant transcripts

`module load vsearch/2.4.3`

`vsearch --threads 16 --log LOG90 \
        --cluster_fast AHA.transdecoder.cds \
        --id 0.90 \
        --centroids AHA.cluster.cds \
        --uc clusters_90.uc`
        
### BUSCO completeness post clustering run

`module load busco/5.0.0`
`module load MetaEuk/4.0`
`module load hmmer/3.2.1`
`module load blast`

`busco -c 12 -m transcriptome -i AHA.cluster.cds -o busco_after -l eudicots_odb10`

### Functional annotation of ORFs using Entap

`module load anaconda/2.4.0`
`module load perl/5.24.0`
`module load diamond/0.9.36`
`module load interproscan/5.25-64.0`

`prot=/path/to/AHA.cluster.pep`

`EnTAP --runP --ini entap_config.ini \
-i $prot \
-d /path/to/database/Diamond/Uniprot/uniprot_sprot.dmnd \
-d /path/to/database/Diamond/RefSeq/plant.protein.faa.205.dmnd \
--threads 16`

### Kallisto for counting transcript/gene abundance

Index the transcriptome first

`module load kallisto/0.46.1`

`cluster=/path/to/cluster`	

`kallisto index -i AHA_index $cluster/AHA.cluster.cds`

Count transcript abundance using Kallisto

`module load kallisto/0.46.1`

`trimmed=/path/to/trimmed/reads

`kallisto quant -i ../index/AHA_index \
        -o AHA_1 -b 100 \
        -t 8 --pseudobam \
        $trimmed/AHA_1_NG_R1_val_1.fq.gz $trimmed/AHA_1_NG_R2_val_2.fq.gz`

kallisto quant -i ../index/AHA_index \
        -o AHA_2 -b 100 \
        -t 8 --pseudobam \
        $trimmed/AHA_2_G_R1_val_1.fq.gz $trimmed/AHA_2_G_R2_val_2.fq.gz

kallisto quant -i ../index/AHA_index \
        -o AHA_2a -b 100 \
        -t 8 --pseudobam \
        $trimmed/AHA_2a_NG_R1_val_1.fq.gz $trimmed/AHA_2a_NG_R2_val_2.fq.gz

kallisto quant -i ../index/AHA_index \
        -o AHA_3 -b 100 \
        -t 8 --pseudobam \
        $trimmed/AHA_3_G_R1_val_1.fq.gz $trimmed/AHA_3_G_R2_val_2.fq.gz`

### Check how good clean RNA reads align to the total consensus transcripts

`module load samtools/1.7`

`samtools flagstat /path/to/kallisto_run/AHA_1/pseudoalignments.bam >> 1.txt`
`samtools flagstat /path/to/kallisto_run/AHA_2/pseudoalignments.bam >> 2.txt`
`samtools flagstat /path/to/kallisto_run/AHA_2a/pseudoalignments.bam >> 2a.txt`
`samtools flagstat /path/to/kallisto_run/AHA_3/pseudoalignments.bam >> 3.txt`

DESeq2 requires a data frame (tx2gene) with transcript ID and gene ID. The transcript ID must be the same ID used in the abundance file (i.e. in the Kallisto count process)
To generate the dataframe tx2gene, we will use the trans_map files generated during the Trinity assembly

##Rename the gene/transcript in the gene_map file
`sed 's/TRINITY/AHA_1_NG_TRINITY/g' AHA_1_NG.Trinity.fasta.gene_trans_map > AHA_1_NG_map`

Repeat this for other libraries too.
##Concatenate all the AHA into single file
`cat AHA_*G_map >> AHA_map`

`awk '{ print $2" " $1}' AHA_map | sed '1 i\TXNAME\tGENEID' | sed 's/ /\t/g' > tx2gene_AHA.tsv`


## Part II: DESeq2 analysis using R

Copy abundance.tsv for each RNA library generated from Kallisto and trans_map file (tx2gene_AHA.tsv) to the local working directory

### Load libraries required for the analysis

`library("tximport")`
`library("readr")`
`library("DESeq2")`
`library("dplyr")`

### Creating a dataframe for DESeq2 analysis. 
Check the link for more information: https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html#kallisto

`outputPrefix <- "aha_deseq2"`
`setwd("~/path/to/deseq2/AHA")`
`filenames <- list.files(pattern="*abundance.tsv", full.names=TRUE)`
`tx2gene <- read.table("tx2gene_AHA.tsv", header=TRUE)`
`txi.kallisto <- tximport(filenames, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = TRUE)`
`sampleNames <- c("AHA_1_NG","AHA_2_G","AHA_3_G","AHA_3_NG")`
`colnames(txi.kallisto$counts) <- sampleNames`
`sampleCondition <- c("NG","G","G","NG")`
`sampleTable <- data.frame(Condition = sampleCondition, Sample = sampleNames)`
`ddstxi <- DESeqDataSetFromTximport(txi.kallisto, sampleTable, design = ~ Condition)`
`treatments <-  c("NG","G")`
`colData(ddstxi)$Condition <- factor(colData(ddstxi)$Condition, levels = treatments)`

Running DESeq2
`dds <- DESeq(ddstxi)`
`res <- results(dds)`
`summary(res)`

Export the results to combine with the annotation. Edit the DESeq2 output accordingly to merge with the annotation file.
`write.csv(res, "AHA.deseq2.raw.csv")`

Convert annotation file "AHA.annotation.csv" generated from functional annotation to a format compatible to DESeq2 output so that they can be merge into a single file.
E.g. remove isoform and protein id (_i1_p1) on the TRINITY header in the annotation file. Remove redundant TRINITY header as well.
Also, add column header for TRINITY gene name that match with DESeq2 csv file so that they can be merge together

`aha_edited <- read.csv("AHA.deseq2.raw.csv")`
`aha_anno <- read.csv("AHA.annotation.csv")`
`aha_merge <- merge(aha_edited,aha_anno,by=c('Query_Sequence'), all.x=T)`
`write.csv(aha_merge, "AHA.deseq2.raw.withannotation.csv")`

The summary is for padj value of 0.1 by default so if you want to change it then
`result_05 <- results(dds, alpha=0.05)`

Summarize the number of genes with padjust less than 0.05 then:
`sum(result_05$padj < 0.05, na.rm=TRUE)`

Subsample the genes with padj <0.05 only
`res05 <- res[ which(res$padj < 0.05 ), ]`

Subsample log2foldchange > 1 only
`res05_up <- res05[ which(res05$log2FoldChange > 1 ), ]`

Subsample log2foldchange < -1
`res05_down <- res05[ which(res05$log2FoldChange < -1 ), ]`

`write.csv(res05_up, "AHA.upregulated.csv")`
`write.csv(res05_down, "AHA.downregulated.csv")`

Write csv file for up and down regulated genes. Manual concatenation is ok (AHA.regulated.csv). Also, header "Query_Sequence" can be added so that annotation can be merged. 
`aha_regulated <- read.csv("AHA.regulated.csv")`
`aha.reg.final <- merge(aha_regulated, aha_anno, by=c('Query_Sequence'), all.x=T)`
`write.csv(aha.reg.final, "AHA.UpDown.annotation.csv")`


#You can repeat DESeq2 analysis for PMO similar to AHA

`outputPrefix <- "pmo_deseq2"`
`setwd("~/path/to/deseq2/PMO")`
`filenames <- list.files(pattern="*abundance.tsv", full.names=TRUE)`
`tx2gene <- read.table("tx2gene_PMO.tsv", header=TRUE)`
`txi.kallisto <- tximport(filenames, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = TRUE)`
`sampleNames <- c("PMO_1_NG","PMO_2_G","PMO_2a_NG","PMO_3_G")`
`colnames(txi.kallisto$counts) <- sampleNames`
`sampleCondition <- c("NG","G","NG","G")`
`sampleTable <- data.frame(Condition = sampleCondition, Sample = sampleNames)`
`ddstxi <- DESeqDataSetFromTximport(txi.kallisto, sampleTable, design = ~ Condition)`
`treatments <-  c("NG","G")`
`colData(ddstxi)$Condition <- factor(colData(ddstxi)$Condition, levels = treatments)`

Running DESeq2
`dds <- DESeq(ddstxi)`
`res <- results(dds)`
`summary(res)`
