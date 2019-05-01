
# Data source  
Poly(A)-/Total long RNA-seq data from ENCODE/Cold Spring Harbor Lab were downloaded from Encyclopedia of DNA Elements ([`ENCODE`](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCshlLongRnaSeq/)) Project Consortium 
<br>

# CircRNA identification and analysis  

CircRNA is relative low expressed and the strategies for the circRNA detection tools also vary a lot, resulting differences in predicted output. To reduce the false-positives of circRNA detection, several widely tools are used to get circRNA candidates.  
<br>
  
## `CIRCexplorer` pipeline 
> The raw RNA-seq data were first aligned to human genome GRCh37/hg19-release75 using [STAR](https://github.com/alexdobin/STAR) version 2.5.1b for chimeric detection using:  
<br>

```Bash
STAR --chimSegmentMin 10 --runThreadN 5 --genomeDir <hg19_STAR_index> --readFilesIn <R1.fastq> <R2.fastq>  
```
> The STAR output file `"Chimeric.out.junction"` was then covert and analyzed by [`CIRCexplorer`](https://github.com/YangLab/CIRCexplorer) version 1.1.10 using:  

```Bash
star_parse.py Chimeric.out.junction fusion_junction.txt
  
CIRCexplorer.py -j fusion_junction.txt -g <hg19.fa> -r <ref.txt>  
```
> The output file of circRNA junction was `"CIRCexplorer_circ.txt"`  
<br>

## `CIRI` pipeline
> The raw reads were aligned to human genome GRCh37/hg19-release75 with BWA-MEM ([BWA](https://github.com/lh3/bwa) version 0.7.15) first to generate SAM file using:  

```Bash
bwa mem –t 10 -T 30 <hg19.fa> <R1.fastq> <R2.fastq> > <SAM>
```

> The BWA output SAM file was then analyzed by [`CIRI`](https://sourceforge.net/projects/ciri/) version 2.0.2 using:  

```Bash
perl CIRI2.pl -I <SAM> -O < outfile > -F <hg19.fa> -A < GRCh37.75.gtf > -T 10
```

> Results of detected circRNA junctions was stored in CIRI outfile.  
<br>


## `KNIFE` pipeline
> The raw reads were aligned to human genome GRCh37/hg19-release75 with [BOWTIE2](http://bowtie-bio.sourceforge.net/bowtie2) version 2.2.7 and analyzed by [`KNIFE`](https://github.com/lindaszabo/KNIFE) version 1.4 using:  

```Bash
./completeRun.sh <read_directory> appended <alignment_parent_directory> <dataset_name> 13 sam_large_phred64 circReads 50 0
```

> Results of detected circRNA junctions was stored in <alignment_parent_directory>/<dataset_name>/circReads/combinedReports/circleJuncProbs.txt  
<br>

# Note  
`CIRCerxplorer` and `CIRI` are recommended for their high sensitivity and decent performance in detecting circRNA candidates from RNA-seq data without poly(A)- or RNase R treatment according to paper of *Zeng, X.* et al\[PMID: [28594838](https://www.ncbi.nlm.nih.gov/pubmed/?term=A+comprehensive+overview+and+evaluation+of+circular+RNA+detection+tools)\].  

CircRNAs with at least 2 junction reads find in at least one sample are regarded as high confidence detection.  
  
Junction read counts are normalized to Spliced Reads per Billion Mapped Reads (SRPBM).  

<br><br>

