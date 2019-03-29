# data source:
https://www.ncbi.nlm.nih.gov/sra (download SraAccList.txt with MEF sras, 
http://www.metagenomics.wiki/tools/short-read/ncbi-sra-file-format/prefetch
prefetch it then fastq-dump)
convert files into fastq files with SRA Toolkit
`fastq-dump /home/laskina/ncbi/public/sra/*.sra -O /project/functional-genomics/2019/data/sra/MEF_G3/prefetched`

'/pkg/python-3.6.6-1/lib/python3.6/site-packages/deeptools/plotCoverage.py'

quality control with FastQC

# Bowtie2 Mapping(Chip-Seq, Atac-Seq)
bowtie2 http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-aligner
#mismatches -> filter from SAM file(XM:i:<N>The number of mismatches in the alignment. Only present if SAM record is for an aligned read.)
``` 
    for i in `ls *.fastq | cut -d "." -f 1` ;
    do /package/sequencer/bowtie2/current/bin/bowtie2 -p 50 -x /project/functional-genomics/2019/data/genome/mm9 -U $i.fastq -S ./mapped/$i.sam ;
done ;
 ```
Create bigWig files from the bam files using bamCoverage:
``` 
for i in `ls /project/functional-genomics/2019/group3/MEF/*.bam | cut -d "." -f 1` ;
	do bamCoverage -b $i.bam -bs 25 --normalizeUsing RPKM -ignore chrX --extendReads 200 -o $i.bw ;
done ;
 ```
Quality Control:
plotFingerprint for all TF's (e.g. MEF):
``` 
./plotFingerprint -b /project/functional-genomics/2019/group3/MEF/SRR5077732_sorted.bam /project/functional-genomics/2019/group3/MEF/SRR5077730_sorted.bam /project/functional-genomics/2019/group3/MEF/SRR5077728_sorted.bam /project/functional-genomics/2019/group3/MEF/SRR5077726_sorted.bam /project/functional-genomics/2019/group3/MEF/SRR5077722_sorted.bam /project/functional-genomics/2019/group3/MEF/SRR5077718_sorted.bam /project/functional-genomics/2019/group3/MEF/SRR5077714_sorted.bam /project/functional-genomics/2019/group3/MEF/SRR5077677_sorted.bam /project/functional-genomics/2019/group3/MEF/SRR5077676_sorted.bam /project/functional-genomics/2019/group3/MEF/SRR5077673_sorted.bam -p max/2 -l Cebpb Cebpa Runx1 Fra1 Brg1 Hdac1 p300 cMyc Klf4 WCE_control -plot /project/functional-genomics/2019/group3/MEF/QC/fingerprint_MEF_TF.png 
 ``` 
Create bamSummary/bigWigSummary (e.g. MEF):
``` 
./multiBamSummary bins -b /project/functional-genomics/2019/group3/MEF/*.bam -bs 25 -p max/2 -o /project/functional-genomics/2019/group3/MEF/QC/results_MEF.npz
./multiBigwigSummary bins -b /project/functional-genomics/2019/group3/MEF/bigWig/*.bw -bs 25 -p max/2 -o /project/functional-genomics/2019/group3/MEF/QC/results_bigWig_MEF.npz
```  
Quality Control:
plotCorrelation from multiBamSummary /multiBigWigSummary (e.g. MEF):
``` 
./plotCorrelation -in /project/functional-genomics/2019/group3/MEF/QC/results_MEF.npz -c pearson -p heatmap -l H3K4me3 H3K4me2 H3K4me1 H3K9ac H3K27ac H3K27me3 H3K79me2 H3K36me3 H3K9me3 H3.3 H3 MNase_control WCE_control Klf4 cMyc p300 Hdac1 Brg1 Fra1 Runx1 Cebpa Cebpb MEF_ATAC-seq MEF_NODOX_ATAC-seq -T Correlation_MEF_rm --removeOutliers --colorMap YlOrRd --outFileCorMatrix /project/functional-genomics/2019/group3/MEF/QC/PearsonCorr_Scores_MEF_rm.tab -o /project/functional-genomics/2019/group3/MEF/QC/correlationMatrix_MEF_rm.png
./plotCorrelation -in /project/functional-genomics/2019/group3/MEF/QC/results_bigWig_MEF.npz -c pearson -p heatmap -l H3K4me3 H3K4me2 H3K4me1 H3K9ac H3K27ac H3K27me3 H3K79me2 H3K36me3 H3K9me3 H3.3 H3 MNase_control WCE_control Klf4 cMyc p300 Hdac1 Brg1 Fra1 Runx1 Cebpa Cebpb MEF_ATAC-seq MEF_NODOX_ATAC-seq -T Correlation_MEF_bigWig_rm --removeOutliers --colorMap YlOrRd --outFileCorMatrix /project/functional-genomics/2019/group3/MEF/QC/PearsonCorr_Scores_bigWig_MEF_rm.tab -o /project/functional-genomics/2019/group3/MEF/QC/correlationMatrix_bigWig_MEF_rm.png
``` 

 # peak calling with MACS2
MACS2 https://github.com/taoliu/MACS

MEF:
~~~
MNaseTreatment="SRR5077653 SRR5077645 SRR5077641 SRR5077637 SRR5077633 SRR5077629 SRR5077625"
for i in $MNaseTreatment;
	do macs2 callpeak -t ${i}_filtered.bam -c SRR5077669_filtered.bam -n $i -f BAM -g mm --bw 150 -q 0.005 --outdir peakcalling ;
done ;

WCETreatment="SRR5077732 SRR5077730 SRR5077728 SRR5077726 SRR5077722 SRR5077718 SRR5077714 SRR5077677 SRR5077676 SRR5077665 SRR5077661 SRR5077657 SRR5077649"
for i in $WCETreatment;
	do macs2 callpeak -t ${i}_filtered.bam -c SRR5077673_filtered.bam -n $i -f BAM -g mm --bw 150 -q 0.005 --outdir peakcalling ;
done ;
~~~
Calculating average width:
~~~
touch avgLength.txt
for i in 'ls *.xls';
	do Rscript calcAvgLength.R $i;
done;
~~~
For screenshots of target genes in IGV see the folder ScreenshotsIGV 
# STAR Mapping(Rna-Seq)
star manual http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STAR.posix/doc/STARmanual.pdf
```
STAR --runThreadN 20 --genomeDir /project/functional-genomics/2019/data/genome/STARindex --readFilesIn /project/functional-genomics/2019/data/sra/MEF_G3/prefetched/RNAseq/SRR5077610.fastq --outFileNamePrefix /project/functional-genomics/2019/data/sra/MEF_G3/prefetched/RNAseq/SRR5077610_ --quantMode GeneCounts --outFilterMismatchNmax 3 --outSAMtype BAM SortedByCoordinate --bamRemoveDuplicatesType UniqueIdentical --outWigType wiggle --outWigStrand Unstranded --outWigNorm RPM
```
## convert wiggle files into bigWig
use unique.str1.wig files, as we are only interested in uniquely mapped reads(?). 
`wigToBigWig /home/laskina/srr.wig sizes.genome /home/laskina/srr10.bw`
Problem: UCSC and STAR Chromosome names are not equal -> error by calling, must rename or ignore other chromosomes(remove all of the chromosomes like "NT_16694" from the list and keep 19+XY Chromosomes(Python script))

## count reads per gene 
use htseq: https://htseq.readthedocs.io/en/release_0.11.1/count.html
htseq-count -f bam -r name -o <alignment_files> <gff_file>
* htseq names differ from the one in STAR or annotation, produce count tables with --quantMode in STAR Mapping.

# peak comparison
bedtools https://bedtools.readthedocs.io/en/latest/


creating unified peak widths:
``` bash unified_summits.sh```


creating many different bedfiles:
``` bash peak_comparison.sh```

Results are in peak_comparison.txt

## comparison with promotor regions
downloaded upstream2000.fa.gz from http://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/

creating a sorted bed file from the fasta file:
``` 
python myFastaToBed.py
bedtools sort -i upstream2000.bed > upstream2000_sorted.bed
```

creating bedfiles :
``` bash peaks_in_promotor.sh ```

Results are in peak_comparison.txt

## Visualization

see folder PeakComparisonPlots

# read counts in promotor regions
MEF:  
```
PROMOTOR=/project/functional-genomics/2019/group3/peak_comparison/upstream2000_sorted.bed
IN_FOLDER=/project/functional-genomics/2019/data/sra/MEF_G3/prefetched/mapped
OUT_FOLDER=/project/functional-genomics/2019/group3/read_counts/MEF
FILE_LIST="SRR5077625 SRR5077629 SRR5077633 SRR5077637 SRR5077641 SRR5077645 SRR5077649 SRR5077653 SRR5077657 SRR5077661 SRR5077665 SRR5077669 SRR5077673"
for i in $FILE_LIST;
	do bedtools coverage -a $PROMOTOR -b $IN_FOLDER/${i}_filtered.bam -counts > $OUT_FOLDER/${i}_readcounts.bed;
	echo "$i is done" ;
done ;
```
