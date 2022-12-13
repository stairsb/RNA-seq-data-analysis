# RNA-seq-data-analysis
This pipeline is for analyzing RNA-seq data

## Table of Contents

## Requirements


## Obtaining the data
For this anaysis we will need a reference genome and the RNA-seq reads. If you already have these downloaded then you can skip this section. I used genomic data on the NCBI SRA database.

Create a file that contains the SRR names corresponding to the genomic data that will be used. The RNA-seq reads are from this paper https://www.pnas.org/doi/10.1073/pnas.2003857117.
```
SRR11651219 RNA-seq of Rhizopus microsporus infected with btl9-13 mutant Mycetohabitans sp. B13 with the plasmid pBtl19-13 biorep 3
SRR11651220 RNA-seq of Rhizopus microsporus infected with btl9-13 mutant Mycetohabitans sp. B13 with the plasmid pBtl19-13 biorep 2
SRR11651221 RNA-seq of Rhizopus microsporus infected with btl9-13 mutant Mycetohabitans sp. B13 with the plasmid pBtl19-13 biorep 1
SRR11651222 RNA-seq of Rhizopus microsporus infected with btl9-13 mutant Mycetohabitans sp. B13 biorep 3
SRR11651223 RNA-seq of Rhizopus microsporus infected with btl9-13 mutant Mycetohabitans sp. B13 biorep 2
SRR11651224 RNA-seq of Rhizopus microsporus infected with btl9-13 mutant Mycetohabitans sp. B13 biorep 1
SRR11651225 RNA-seq of Rhizopus microsporus infected with WT Mycetohabitans sp. B13 biorep 3
SRR11651226 RNA-seq of Rhizopus microsporus infected with WT Mycetohabitans sp. B13 biorep 2
SRR11651227 RNA-seq of Rhizopus microsporus infected with WT Mycetohabitans sp. B13 biorep 1
SRR4063850  Rhizopus microsporus ATCC52814 Standard Draft (Whole genome sequence data)
```
To download reads from the SRA database we will use the program fasterq-dump https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump. 
```
awk '{print $1}' SRR_names.txt |
while read i;
do
        echo $i
        fasterq-dump -O . --progress --threads 8 --split-files $i
done
```
The shorts reads will be downloaded in the directory that you ran the bash script in. Let's move them into a data directory
```
mkdir data
mv *fastq ./data/
```
We also need to download the reference genome and gtf file associated with its annotation found here https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/083/745/GCA_002083745.1_Rhimi_ATCC52814_1/. First we will make a directory to where we will keep our reference genome and supporting files. 
Run this command if you are still in the data directory

`cd ../`
Then create the dirctory you will be working in
```
mkdir reference
cd reference
```
We will use wget to download the files
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/083/745/GCA_002083745.1_Rhimi_ATCC52814_1/GCA_002083745.1_Rhimi_ATCC52814_1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/083/745/GCA_002083745.1_Rhimi_ATCC52814_1/GCA_002083745.1_Rhimi_ATCC52814_1_genomic.gtf.gz
```
There is in issue with gtf format in the downstream analysis. Run this command to fix it.

`awk '$3 != "gene" ' GCA_002083745.1_Rhimi_ATCC52814_1_genomic.gtf > Rmicrosporus52814.gtf`

## Aligning reads
To align the reads will use use the program HISAT2 http://daehwankimlab.github.io/hisat2/. We will use HISAT2 to map out RNA-seq reads to our reference genome and index the reference genome which will speed up our downstream compuational time.

Let's extract the exons from the GTF file and save this file in the hisat2_index directory (This step is needed to get rid of genomic information that will not be useful in this analysis. First move to your base working directory for this project
```
mkdir hisat2_index
cd hisat2_index
hisat2/hisat2_extract_exons.py ../reference/GCA_002083745.1_Rhimi_ATCC52814_1_genomic.gtf >./GCA52814.exon
```
Next, run this code to extract the splice sites from the annotation file and save it in the hisat2_index directory. 

`hisat2/hisat2_extract_splice_sites.py ../reference/GCA_002083745.1_Rhimi_ATCC52814_1_genomic.gtf > ./GCA52814.ss`

To make sure everything has worked so far run this command to view exons. To check the splice file change the input file to the `.ss` file that was created earlier
```
cat GCA52814.exon | head
```
Output:
```
KV921853.1	364	1059	-
KV921853.1	1112	1405	-
KV921853.1	1459	2131	-
KV921853.1	2179	2636	-
KV921853.1	2697	2784	-
KV921853.1	2831	2948	-
KV921853.1	3000	3210	-
KV921853.1	3273	3308	-
KV921853.1	3375	3486	+
KV921853.1	3538	4035	+
```
Here we will create index that cross references the exon and splice site annotations with the reference genome and convert it to a binary file.
```
hisat2/hisat2-build --ss ./GCA52814.ss --exon ./GCA52814.exon ../reference/GCA_002083745.1_Rhimi_ATCC52814_1_genomic.fna ./GCA52814_idx
```
You should have 8 indecies files. Run this command to double check that the previous step worked.
```
ls -1 *idx.?.ht2
```
Output:
```
GCA52814_idx.1.ht2
GCA52814_idx.2.ht2
GCA52814_idx.3.ht2
GCA52814_idx.4.ht2
GCA52814_idx.5.ht2
GCA52814_idx.6.ht2
GCA52814_idx.7.ht2
GCA52814_idx.8.ht2
```
Paste the following into a bash script which will align reads to the index using hisat2, convert the `.sam` file to a `.bam` using `samtools view` http://www.htslib.org/doc/samtools-view.html, sort the `.bam` file using `samtools sort` http://www.htslib.org/doc/samtools-sort.html and use STRINGTIE to pair each aligned read to a particular gene or transcript of interest and generate some coverage data. Before running this script create a directory in your base working directory `mkdir aligned`. Finally move into the `data` directory and run this script.
```
#!/bin/bash

names=`ls -1 *.fastq | cut -c -11 | uniq`
for i in $names
do
        echo "Input: sample $i"
        echo "Input: sample $i.gtf"
echo " Running hisat"
touch results.txt

hisat2/hisat2 --threads 4 --dta -x ../hisa2_index/GCA52814_idx -1 "$i"_1.fastq -2 "$i"_2.fastq -S ../aligned/$i.sam 2>> results.txt

echo "Running samtools"
samtools view -Sb ../aligned/$i.sam > ../aligned/${i}_unsorted.bam

samtools sort -@ 4  ../aligned/${i}_unsorted.bam -o ../aligned/$i.bam 

echo "Running stringtie"
stringtie -p 4 -e -G ../reference/Rmicrosporus52814.gtf -o ../aligned/$i.gtf ../aligned/$i.bam

done
```
If there is issue with the script above replace the relative paths with absolute paths

To double check and make sure that the script work move into the `aligned` directory and run:
```
head -n 3 Nameoffile.gtf
```
Output:
```
KV921853.1	StringTie	transcript	90325	90684	1000	+	.	gene_id "BCV72DRAFT_218780"; transcript_id "gnl|WGS:MCGZ|BCV72DRAFT_mRNA218780"; cov "2.113889"; FPKM "0.335535"; TPM "0.365274";
KV921853.1	StringTie	exon	90325	90684	1000	+	.	gene_id "BCV72DRAFT_218780"; transcript_id "gnl|WGS:MCGZ|BCV72DRAFT_mRNA218780"; exon_number "1"; cov "2.113889";
KV921853.1	StringTie	transcript	78870	81768	1000	+	.	gene_id "BCV72DRAFT_300342"; transcript_id "gnl|WGS:MCGZ|BCV72DRAFT_mRNA300342"; cov "82.708923"; FPKM "13.128287"; TPM "14.291849";
```

Extract the alignment rate for each sample and save this in a new (temporary file) called rate.tmp

`grep 'overall alignment rate' results.txt | cut -c -5 > rate.tmp`

Get sample names into a file

`ls -1 *fastq | cut -c -11 | uniq > names.tmp`

Assign factors to each group

`echo plasmid plasmid plasmid mutant mutant mutant wildtype wildtype wildtype | tr " " "\n" > group.tmp`

Finally incorporate all this information into a single file. You can delete the `.tmp` file when you are done

`paste names.tmp rate.tmp group.tmp > mapping.txt`

Output:
```
SRR11651219     95.29   plasmid
SRR11651220     95.07   plasmid
SRR11651221     95.58   plasmid
SRR11651222     95.26   mutant
SRR11651223     95.53   mutant
SRR11651224     95.14   mutant
SRR11651225     94.87   wildtype
SRR11651226     95.11   wildtype
SRR11651227     94.79   wildtype
```

## Mapping Rate
This code will output a plot that shows us the mapping rate of each of the samples to the reference genome. If you don't have ggplot2 installed it is avaliable on the cran.
```
library(ggplot2)
mapping <- read.delim("mapping.txt",header=F)
View(mapping)
mapping$V1<-factor(mapping$V1, levels = mapping$V1[order(mapping$V3)])

ggplot(mapping, aes(x=V1, y = V2, fill = V3)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text = element_text(angle = 90)) +
  ylab("Hisat2 Mapping Rate %") +
  xlab("Sample") +
  ylim(c(0,100)) +
  scale_color_tron() +
  scale_fill_tron() +
  theme(legend.title = element_blank())
```
![image](https://user-images.githubusercontent.com/111078377/206634110-e921e5cf-a774-4e4a-9b40-a7890b064b7a.png)
## Creating a gene count matrix
Lets create a file containing the sample names for each `.gtf` file
```
ls -1 *.gtf > samplelist.txt
awk '{print "sample"NR, "\t", $0}' < samplelist.txt > Samples.txt
```
Now we will obtain counts among all of our samples found in the  `.gtf` files and create a matrix. The prepDE.py is in the files section above.
```
prepDE.py -i Samples.txt
```
Output:
```
gene_count_matrix.csv
transcript_count_matrix.csv
```
Create a file containing the sample names and their condition named `design.txt`.
```
sample condition
sample1 plasmid
sample2 plasmid
sample3 plasmid
sample4 mutant
sample5 mutant
sample6 mutant
sample7 wildtype
sample8 wildtype
sample9 wildtype
```
## Installing R-packages
If you don't have BiocManager installed then run the following code. 
```
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
```
Here are the packages found on BiocManager that you will needed to have installed. Only run the code below for the packages that you don't already have installed.
```
BiocManager::install("DESeq2")
BiocManager::install("tximport")
BiocManager::install("tximportData")
```
The packages below are ones that are found on cran. Only install ones that you don't already have installed.
```
install.packages("hexbin")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("readr")
install.packages("adegenet")
install.packages("pheatmap")
install.packages("genefilter")
install.packages("RColorBrewer")
install.packages("ggsci")
```
Now we will load in all of the packages that we will be using (Note: you might already have ggplot2 loaded). 
```
library(DESeq2)
library(tximportData)
library(tximport)
library(dplyr)
library(ggplot2)
library(hexbin)
library(readr)
```
## Visualizing the data using DEseq2
First we will read in the data from the files that we created above. The first line of code reads in a gene count matrix and the second contains the sample names we are using to investigate differential gene expression.
```
countData<- as.matrix(read.csv("gene_count_matrix.csv", row.names="gene_id"))
colData <-read.table("design.txt", row.names = "sample", header = TRUE, sep = "")
```
Next we will put data into a DEseq object.
```
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
```
To lower the compuational tax load we will remove genes with 1 or less reads mapped to it since they won't be relivant in this analysis.
```
dds <-dds[ rowSums(counts(dds)) > 1, ]
```
Here we will transform the data using a "regularized-logarithm transformation" function that is part of DESeq2 which is recomended for small datasets.
```
rld <- rlog(dds)
```
As a checkin here is an example of what our transformed data should look like.
```
 head(assay(rld))
 ```
 Output:
 ```
                   sample1   sample2   sample3   sample4   sample5   sample6   sample7   sample8   sample9
BCV72DRAFT_304378 11.575119 11.682836 11.628389 11.730588 11.705235 11.612698 11.756440 11.680450 11.495844
BCV72DRAFT_240065  9.266680  8.238637  8.531259  8.092372  7.913480  7.941938  7.836073  7.786046  7.791481
BCV72DRAFT_236859  3.801870  3.793351  3.775463  3.775711  3.794300  3.775648  3.843152  3.829419  3.865140
BCV72DRAFT_249555 12.955227 12.919595 12.902062 12.926474 12.920059 12.907649 12.972063 12.909273 12.929597
BCV72DRAFT_321890  1.767938  1.756862  1.769829  1.768753  1.763930  1.775355  1.756795  1.763003  1.756849
BCV72DRAFT_329775  3.263824  3.257570  3.243970  3.240432  3.254201  3.249443  3.270572  3.243337  3.263488
```
We can compare the transfored data of 2 samples of our choosing using the following code.
```
df <-as.data.frame(assay(rld)[, 4:5])
View(df)
colnames(df)[1:2] <- c("Sample_4", "Sample_5") 
```
To plot our comparison of the two samples.
```
ggplot(df, aes(x = Sample_4, y = Sample_5)) + 
  geom_hex(bins = 100, colour="orange", fill="black") + 
  coord_fixed()+theme_classic()
```
![image](https://user-images.githubusercontent.com/111078377/207260881-a64664fc-0084-473e-91b2-04297a12ef7b.png)

We could run `plotPCA(rld)` function to plot the data but first lets make some visual changes to the plot using ggplot2.
Prepping the data:
```
library(ggsci)
pca1<-plotPCA(rld, intgroup="condition",returnData=TRUE)
percVar <- format(100 * attr(pca1, "percentVar"),digits=3) 
```
New plot:
```
ggplot(pca1, aes(PC1, PC2, colour=condition)) + 
  geom_point(size=3) +
  coord_fixed() +
  scale_color_tron() +
  scale_fill_tron() +
  ggtitle("Reference Genome DGE Analysis")+
  xlab(paste0("PC1: ",percVar[1],"% variance"))+
  ylab(paste0("PC2: ",percVar[2],"% variance"))+
  theme_classic()
```
![image](https://user-images.githubusercontent.com/111078377/207262649-b34bf570-9732-4cf3-a2ec-e318f680b14d.png)

Create a loadings plot to visualize the genes which had the most variation. We can set a threshold to help capture genes of interest.
```
library(adegenet)
t_rld<-t(assay((rld)))
t_rld[1:3,1:4]
pca<-prcomp(t_rld)
loadings<-as.data.frame(pca$rotation)
plot_load<-loadingplot(loadings$PC1,threshold = 0.1)
```
Using the names function we can check that out output names stored in the object `pca` look like this:
```
names(pca)
"sdev" "rotation" "center" "scale" "x" 
```
Output plot:
![image](https://user-images.githubusercontent.com/111078377/207264632-aab28e8b-0f3f-4e36-a495-483f57359619.png)

Now we can grab the list of gene names for genes that are above out threshold.
```
row.names(dds)[plot_load$var.names]
```
## Differential Gene Expression Analysis
Transforming the data (raw read counts) and fitting it to a statistical model.
```
dds <-DESeq(dds)
```
Extracting key information from the DESeq oject:
```
str(dds)
res <- results(dds)
head(res)
mcols(res, use.names = TRUE)
summary(res)
```
The function plotMA shows the log2 fold changes attributable to a given gene over the mean of normalized counts for all the samples. Points that are colored blue if their adjusted p value is less than 0.1. You can adjust the p-value threshold for graphing purposes, but let's just stick with the default for now. Points that exceed y-lim are plotted as open triangles. You might have to run the following code in your consol to be able to interact with the plot. Stop the interactive session once you have selection all of your points of interest.
```
plotMA( res, ylim = c(-10, 10))
idx <- identify(res$baseMean, res$log2FoldChange)
idx
rownames(res)[idx]
```
![image](https://user-images.githubusercontent.com/111078377/207268308-8f131f5b-0d04-4537-b97b-37a6c2be276c.png)

The `idx` variable cotains the row numbers for genes out interest that we found. The `327` in the function below correlates to one of the genes, we can visualize the expressions levels of this gene for each of the samples we included.
```
sigGene<-row.names(res[327,])
plotCounts(dds, gene = sigGene)
topVarGenes  <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 10 )
```
![image](https://user-images.githubusercontent.com/111078377/207269940-5a098708-8ad0-439a-a543-db3d5c214b3d.png)

Finally, create a heatmap by looking at the amount each gene deviates in a specific sample from the gene’s average across all samples.
```
library(genefilter)
library(pheatmap)
library(RColorBrewer)
mat  <- assay(rld)[ topVarGenes, ]
condition = c("orange", "black","yellow")
names(condition) = c("wiltype", "mutant", "plasmid")
ann_colors = list(condition=condition)
pheatmap(mat, annotation_col = colData, color=colorRampPalette(c( "blue","white", "red"))(20),annotation_colors =ann_colors[2],border_color = 'black',fontsize = 12,scale = "row",angle_col = 45)
```
![image](https://user-images.githubusercontent.com/111078377/207271563-7307bc18-e7f5-4327-912e-05f4334af917.png)

Let's reorder the `res` data set based on adjusted p-values by descending order and save the top 20 most differentially expressed genes to a new file.
```
print(res)
resOrdered <- res[order(res$pvalue),]
resOrderedDF <-as.data.frame(resOrdered)[1:20, ]
write.csv(resOrderedDF, file = "results_reference.csv")
```
To view the 20 gene names:
```
row.names(resOrderedDF)

"BCV72DRAFT_76322"  "BCV72DRAFT_226884" "BCV72DRAFT_231122" "BCV72DRAFT_300595"
"BCV72DRAFT_247525" "BCV72DRAFT_214971" "BCV72DRAFT_54254"  "BCV72DRAFT_61846" 
"BCV72DRAFT_7021"   "BCV72DRAFT_227356" "BCV72DRAFT_293493" "BCV72DRAFT_298427"
"BCV72DRAFT_54911"  "BCV72DRAFT_294654" "BCV72DRAFT_258870" "BCV72DRAFT_314435"
"BCV72DRAFT_226292" "BCV72DRAFT_335919" "BCV72DRAFT_339250" "BCV72DRAFT_244019"
```
Copy the genes into a file named `topgenes.txt`
```
BCV72DRAFT_76322
BCV72DRAFT_226884
BCV72DRAFT_231122
BCV72DRAFT_300595
BCV72DRAFT_247525 
BCV72DRAFT_214971 
BCV72DRAFT_54254  
BCV72DRAFT_61846 
BCV72DRAFT_7021  
BCV72DRAFT_227356 
BCV72DRAFT_293493
BCV72DRAFT_298427
BCV72DRAFT_54911 
BCV72DRAFT_294654
BCV72DRAFT_258870
BCV72DRAFT_314435
BCV72DRAFT_226292
BCV72DRAFT_335919 
BCV72DRAFT_339250
BCV72DRAFT_244019
```
The names for this example are specific to the example data used. These can be used to trace back the genes which where the most differentially expressed and figure out what functional proteins they produce.




