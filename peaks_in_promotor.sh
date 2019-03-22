promotor=/project/functional-genomics/2019/group3/peak_comparison/upstream2000_sorted.bed
MEF_K=/project/functional-genomics/2019/data/sra/MEF_G3/prefetched/mapped/peakcalling/SRR5077676_summits_unified.bed
MEF_M=/project/functional-genomics/2019/data/sra/MEF_G3/prefetched/mapped/peakcalling/SRR5077677_summits_unified.bed
ESC_O=/project/functional-genomics/2019/data/sra/MEF_G3/prefetched/mapped_ESC/peakcalling/SRR5077690_1_summits_unified.bed
ESC_S=/project/functional-genomics/2019/data/sra/MEF_G3/prefetched/mapped_ESC/peakcalling/SRR5077691_1_summits_unified.bed
ESC_K=/project/functional-genomics/2019/data/sra/MEF_G3/prefetched/mapped_ESC/peakcalling/SRR5077692_1_summits_unified.bed
ESC_M=/project/functional-genomics/2019/data/sra/MEF_G3/prefetched/mapped_ESC/peakcalling/SRR5077693_1_summits_unified.bed
folder=/project/functional-genomics/2019/group3/peak_comparison

bedtools intersect -a $promotor -b $MEF_K > $folder/promotor_MEF_K.bed

bedtools intersect -a $promotor -b $MEF_M > $folder/promotor_MEF_M.bed

bedtools intersect -a $promotor -b $ESC_O > $folder/promotor_ESC_O.bed

bedtools intersect -a $promotor -b $ESC_S > $folder/promotor_ESC_S.bed

bedtools intersect -a $promotor -b $ESC_K > $folder/promotor_ESC_K.bed

bedtools intersect -a $promotor -b $ESC_M > $folder/promotor_ESC_M.bed
