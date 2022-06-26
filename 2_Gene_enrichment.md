# Gene Enrichment Analysis using GoSeq
Following the protocol- https://github.com/fishercera/TreehopperSeq/blob/master/GoSeq_Walkthrough.md

## File preparation
### First file: Output of EnTAP with GO terms

Staring with "AHA.annotation.csv", this is the output of Entap run that has been modified removing transcripts without any match and isoform ID in the transcript name.

"AHA.annotation.csv" was further modified with removing all the column except transcript ID and three GO terms, Biological, Cellular and Molecular. Save the file as "AHA.annotation.edited.tsv"

Now, the file will be be further edited using the script "GoTermsMap.py" available at: https://github.com/fishercera/TreehopperSeq/blob/master/GoTermsMap.py.
Command usage:
`python GoTermsMap.py AHA.annotation.edited.tsv outfile`

The output file should look like: https://github.com/fishercera/TreehopperSeq/blob/master/GoSeq_Walkthrough.md

As this file contains lines with empty GO information (GO:"), these lines needs to be removed.
`grep -vw 'GO:"' outfile > outfile-edited`

Edit to remove comma and " after the name. 

Add third column with GO term only. First, remove extra information associated with GO term.
`awk -F "-" '{print $1}' outfile-edited  > outfile-go`

Now, outfile-edited and outfile-go can be combined to a single file "AHA_GOterms"", which should appear as:

AHA_1_NG_TRINITY_DN0_c0_g1	GO:0000003-reproduction(L=1)	GO:0000003
AHA_1_NG_TRINITY_DN0_c0_g1	GO:0009987-cellular process(L=1)	GO:0009987
AHA_1_NG_TRINITY_DN0_c0_g1	GO:0032501-multicellular organismal process(L=1)	GO:0032501
AHA_1_NG_TRINITY_DN0_c0_g1	GO:0032502-developmental process(L=1)	GO:0032502
 
*Check whether GO terms do not contain unusal GO term formats. I had a line with unusal format, such as "GO:PFAM", which I deleted manually. It can cause problem later while running goseq function.* 
This represents all the genes from EnTAP output excluding the genes lacking GO terms.

### Second file: Output of Kallisto count with gene length
There were four Kallisto counts (abundance.tsv) per species and I will be using only one. Each of these 
four files have same gene ID and length so it should not matter which one you pick.

The abundance file is modified such that it contains two columns with gene ID and its length, which was saved as "AHA_abundance.tsv".

It should appear as:

target_id	length
AHA_1_NG_TRINITY_DN1443_c0_g1	16419
AHA_1_NG_TRINITY_DN2242_c0_g1	15297
AHA_1_NG_TRINITY_DN2619_c0_g1	14319
AHA_2_G_TRINITY_DN4021_c0_g1	12888

## GO enrichment analysis
A great detail information on GOseq is available in https://bioconductor.org/packages/devel/bioc/vignettes/goseq/inst/doc/goseq.pdf

It requires to prepare a vector with all the genes for which differential gene expression was carried out and the element of vector with 0 and 1 representing differentially expressed or not, respectively.

### Creating vector
Create an empty vector representing all the genes present in the analysis. For this, I am using abundance file but DESeq output file can be used too!

Reading the output of DeSeq2 to get gene names
`resdata <- read.table("AHA.deseq2.raw.tsv", header=TRUE, sep="\t")`

OR

Reading abundance file. However, abundance file does carry isoforms with different lengths so keep the longest one and remove shorter ones.
`
temp <- read.table("AHA_abundance.tsv", header = TRUE, sep = "\t") # reading abundance file
temp %>% distinct %>% group_by(target_id) %>% top_n(1, length) -> resdata
#write.table(resdata, "aha.abundance.longest.tsv", sep = "\t", row.names=FALSE, col.names = TRUE)
`
Then,

`gene.data <- integer(length=length(resdata$target_id))` #Empty vector

`names(gene.data) <- resdata$target_id` #name the items in vector

`UP <- read.table("AHA.up.tsv", header = TRUE, sep = "\t")` #Get the list of upregulated genes (AHA with log2fold >= 1)

`AbdUp <- UP$Query_Sequence` ##get id of the genes

` #use for loop to assign 1 for the upregulated genes in the newly created vector
for (item in AbdUp) {						
  i <- which(names(gene.data)==item)
  gene.data[i] <- 1
}
`
To check whether the assigned (1) upregulated genes are correct
`length(AbdUp) == table(gene.data)[2]`
`table(gene.data)`

Remove if any NA are present.
`genes.list <- na.omit(gene.data) # this will the input data`
`head(genes.list)`

Getting gene length 
`feature.lengths <- read.table("aha.abundance.longest.tsv", header=TRUE, sep="\t")`
`head(feature.lengths)`

Add the length information from abundance file to the newly created vector with DE genes assigned 1
`genes.length.data <- filter(feature.lengths, target_id %in% names(genes.list))`
`head(genes.length.data)`

Make the length data into a vector with names.
`genes.bias.data <- genes.length.data$length`  #for length
`names(genes.bias.data) <- genes.length.data$target_id`  #for gene names


Getting GO information from the EnTAP output. Make sure to use the edited file with only "Biological GO" terminologies. 
`cat.map <- read.table("AHA.GO", header=FALSE, sep="\t")`
`GO_cats <- cat.map[,c(1,3)]`
`head(GO_cats)`

Fit and plot the probability weighting function (pwf) for the length bias
`pwf <- nullp(genes.list, bias.data = genes.bias.data, plot.fit = FALSE)`
Screen output: Warning message: In pcls(G) : initial point very close to some inequality constraints

`plotPWF(pwf = pwf, binsize = 400)` # Bin size can be changed.
Screen output: Warning message: In pcls(G) : initial point very close to some inequality constraints

Perform the Wallenius test (use the goseq function) to get the most enriched go cats
`GO.wall <- goseq(pwf, gene2cat = GO_cats, method = "Wallenius", use_genes_without_cat = FALSE)`

Screen output
Using manually entered categories.
For 12067 genes, we could not find any categories. These genes will be excluded.
To force their use, please run with use_genes_without_cat=TRUE (see documentation).
This was the default behavior for version 1.15.1 and earlier.
Calculating the p-values...
'select()' returned 1:1 mapping between keys and columns

`write.table(GO.wall, "AHA.UP.GOterms.add.tsv", sep="\t")`  #saving the output

Checking top 25 up regulated categories
`sorted.GO <- arrange(GO.wall, GO.wall$over_represented_pvalue)` 
`sorted.GO`
`sorted.GO[,c(1,6)]`