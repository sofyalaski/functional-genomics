# data source:
https://www.ncbi.nlm.nih.gov/sra (download SraAccList.txt with MEF sras, 
http://www.metagenomics.wiki/tools/short-read/ncbi-sra-file-format/prefetch
prefetch it then fastq-dump)
convert files into fastq files with SRA Toolkit
`fastq-dump /home/laskina/ncbi/public/sra/*.sra -O /project/functional-genomics/2019/data/sra/MEF_G3/prefetched`

'/pkg/python-3.6.6-1/lib/python3.6/site-packages/deeptools/plotCoverage.py'


# Bowtie2 Mapping(Chip-Seq, Atac-Seq)
bowtie2 http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-aligner
#mismatches -> filter from SAM file(XM:i:<N>The number of mismatches in the alignment. Only present if SAM record is for an aligned read.)
``` 
    for i in `ls *.fastq | cut -d "." -f 1` ;
    do /package/sequencer/bowtie2/current/bin/bowtie2 -p 50 -x /project/functional-genomics/2019/data/genome/mm9 -U $i.fastq -S ./mapped/$i.sam ;
done ;
 ```
filter bam files that habe quality <10 and >2 mismatches, remove duplicates
	`samtools view -Sh -@ 50 -q 10 -F 1024 SRR5077625.sam | grep -e "^[@]" -e "XM:i:[012]" | samtools view -@ 50 -bSo SRR5077625_filtered.bam`
Create bigWig files from the bam files using bamCoverage:
``` 
for i in `ls /project/functional-genomics/2019/group3/MEF/*.bam | cut -d "." -f 1` ;
	do bamCoverage -b $i.bam -bs 25 --normalizeUsing RPKM -ignore chrX --extendReads 200 -o $i.bw ;
done ;
 ```
	
# Peak calling with MACS2
MACS2 https://github.com/taoliu/MACS
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

Calculating the average wodht of peaks:
```
touch avgLength.txt
for i in 'ls *.xls';
	do Rscript calcAvgLength.R $i;
done;
```
# STAR Mapping(Rna-Seq)
star manual http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STAR.posix/doc/STARmanual.pdf
`STAR --runThreadN 20 --genomeDir /project/functional-genomics/2019/data/genome/STARindex --readFilesIn /project/functional-genomics/2019/data/sra/MEF_G3/prefetched/RNAseq/SRR5077610.fastq --outFileNamePrefix /project/functional-genomics/2019/data/sra/MEF_G3/prefetched/RNAseq/SRR5077610_ --outFilterMismatchNmax 3 --outSAMtype BAM SortedByCoordinate --bamRemoveDuplicatesType UniqueIdentical --outWigType wiggle --outWigStrand Unstranded --outWigNorm RPM`

## UCSC doesn't recognize all of chromosome names

remove all of the chromosomes like "NT_16694" from the list and keep 19+XY Chromosomes(Python script filter_wiggle.py)
``` 
for j in `ls *.wig ` ;
	do wigToBigWig rna_bigWig/${j} sizes.genome   rna_bigWig/${j |cut -d "_" -f 1 }.bw ;
done;
```

## convert wiggle files into bigWig
use unique.str1.wig files, as we are only interested in uniquely mapped reads(?). 
`wigToBigWig /home/laskina/srr.wig sizes.genome /home/laskina/srr10.bw`
Problem: UCSC and STAR Chromosome names are not equal -> error by calling, must rename or ignore other chromosomes(done in previous step)

## count reads per gene 
use htseq: https://htseq.readthedocs.io/en/release_0.11.1/count.html
`htseq-count -f bam -r name -o <alignment_files> <gff_file>`



