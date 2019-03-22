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


# STAR Mapping(Rna-Seq)
star manual http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STAR.posix/doc/STARmanual.pdf
`STAR --runThreadN 20 --genomeDir /project/functional-genomics/2019/data/genome/STARindex --readFilesIn /project/functional-genomics/2019/data/sra/MEF_G3/prefetched/RNAseq/SRR5077610.fastq --outFileNamePrefix /project/functional-genomics/2019/data/sra/MEF_G3/prefetched/RNAseq/SRR5077610_ --outFilterMismatchNmax 3 --outSAMtype BAM SortedByCoordinate --bamRemoveDuplicatesType UniqueIdentical --outWigType wiggle --outWigStrand Unstranded --outWigNorm RPM`

## UCSC doesnt recognize all of chromosome names

remove all of the chromosomes like "NT_16694" from the list and keep 19+XY Chromosomes(Python script)

## convert wiggle files into bigWig
use unique.str1.wig files, as we are only interested in uniquely mapped reads(?). 
`wigToBigWig /home/laskina/srr.wig sizes.genome /home/laskina/srr10.bw`
Problem: UCSC and STAR Chromosome names are not equal -> error by calling, must rename or ignore other chromosomes



