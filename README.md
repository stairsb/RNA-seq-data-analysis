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
This code will output a plot that shows us the mapping rate of each of the samples to the reference genome
```
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








