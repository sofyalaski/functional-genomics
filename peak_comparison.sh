MEF_K=/project/functional-genomics/2019/data/sra/MEF_G3/prefetched/mapped/peakcalling/SRR5077676_summits_unified.bed
MEF_M=/project/functional-genomics/2019/data/sra/MEF_G3/prefetched/mapped/peakcalling/SRR5077677_summits_unified.bed
ESC_O=/project/functional-genomics/2019/data/sra/MEF_G3/prefetched/mapped_ESC/peakcalling/SRR5077690_1_summits_unified.bed
ESC_S=/project/functional-genomics/2019/data/sra/MEF_G3/prefetched/mapped_ESC/peakcalling/SRR5077691_1_summits_unified.bed
ESC_K=/project/functional-genomics/2019/data/sra/MEF_G3/prefetched/mapped_ESC/peakcalling/SRR5077692_1_summits_unified.bed
ESC_M=/project/functional-genomics/2019/data/sra/MEF_G3/prefetched/mapped_ESC/peakcalling/SRR5077693_1_summits_unified.bed
folder=/project/functional-genomics/2019/group3/peak_comparison
#KLF4 peak overlaps between MEF and ESC
bedtools intersect -a $MEF_K -b $ESC_K > $folder/KLF4_intersect.bed 
#cMyc peak overlaps between between MEF and ESC
bedtools intersect -a $MEF_M -b $ESC_M > $folder/cMyc_intersect.bed 
#KLF4 unique peaks in MEF
bedtools intersect -a $MEF_K -b $ESC_K -v > $folder/KLF4_unique_MEF.bed 
#KLF4 unique peaks in ESC
bedtools intersect -b $MEF_K -a $ESC_K -v > $folder/KLF4_unique_ESC.bed 
#cMyc unique peaks in MEF
bedtools intersect -a $MEF_M -b $ESC_M -v > $folder/cMyc_unique_MEF.bed 
#cMyc unique peaks in ESC
bedtools intersect -b $MEF_M -a $ESC_M -v > $folder/cMyc_unique_ESC.bed 
#shared peaks between KLF4 and cMyc in MEF
bedtools intersect -a $MEF_K -b $MEF_M > $folder/MEF_K_M_intersect.bed
#shared peaks between Oct4 and SOX2 in ESC
bedtools intersect -a $ESC_O -b $ESC_S > $folder/ESC_O_S_intersect.bed
#shared peaks between Oct4 and KLF4 in ESC
bedtools intersect -a $ESC_O -b $ESC_K > $folder/ESC_O_K_intersect.bed
#shared peaks between Oct4 and cMyc in ESC
bedtools intersect -a $ESC_O -b $ESC_M > $folder/ESC_O_M_intersect.bed
#shared peaks between SOX2 and KLF4 in ESC
bedtools intersect -a $ESC_S -b $ESC_K > $folder/ESC_S_K_intersect.bed
#shared peaks between SOX2 and cMyc in ESC
bedtools intersect -a $ESC_S -b $ESC_M > $folder/ESC_S_M_intersect.bed
#shared peaks between KLF4 and cMyc in ESC
bedtools intersect -a $ESC_K -b $ESC_M > $folder/ESC_K_M_intersect.bed

#shared peaks between Oct4, SOX2 and KLF4 in ESC
bedops --merge $ESC_O $ESC_S > $folder/ESC_O_S_merged.bed
bedtools intersect -a $folder/ESC_O_S_merged.bed -b $ESC_K > $folder/ESC_O_S_K_intersect.bed
#shared peaks between Oct4, SOX2 and cMyc in ESC
bedtools intersect -a $folder/ESC_O_S_merged.bed -b $ESC_M > $folder/ESC_O_S_M_intersect.bed
#shared peaks between SOX2, KLF4 and cMyc in ESC
bedops --merge $ESC_S $ESC_K > $folder/ESC_S_K_merged.bed
bedtools intersect -a $folder/ESC_S_K_merged.bed -b $ESC_M > $folder/ESC_S_K_M_intersect.bed
#shared peaks between Oct4, KLF4 and cMyc in ESC
bedops --merge $ESC_O $ESC_K > $folder/ESC_O_K_merged.bed
bedtools intersect -a $folder/ESC_O_K_merged.bed -b $ESC_M > $folder/ESC_O_K_M_intersect.bed
#shared peaks between all TFs in ESC
bedops --merge $ESC_O $ESC_S $ESC_K >$folder/ESC_O_S_K_merged.bed
bedtools intersect -a $folder/ESC_O_S_K_merged.bed -b $ESC_M > $folder/ESC_all_intersect.bed



