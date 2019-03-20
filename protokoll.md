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
`STAR --runThreadN 20 --genomeDir /project/functional-genomics/2019/data/genome/STARindex --readFilesIn /project/functional-genomics/2019/data/sra/MEF_G3/prefetched/RNAseq/SRR5077610.fastq --outFileNamePrefix /project/functional-genomics/2019/data/sra/MEF_G3/prefetched/RNAseq/SRR5077610_ --outFilterMismatchNmax 3 --outSAMtype BAM SortedByCoordinate --bamRemoveDuplicatesType UniqueIdentical --outWigType wiggle --outWigStrand Unstranded --outWigNorm RPM`
