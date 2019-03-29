library(VennDiagram)

pdf("venndiagram_KLF4.pdf")
KLF4.venn <- draw.pairwise.venn(area1=24981, area2=20209, cross.area=5569, category=c("ESC","MEF"), fill=c("darkorange3","darkorange3") )
grid.draw(KLF4.venn)
dev.off()

pdf("venndiagram_cMyc.pdf")
cMyc.venn <- draw.pairwise.venn(area1=3482, area2=22039, cross.area=2231, category=c("ESC","MEF"), fill=c("darkolivegreen4","darkolivegreen4") )
grid.draw(cMyc.venn)
dev.off()


pdf("venndiagrams_shared_peaks.pdf")
O_S.venn <- draw.pairwise.venn(area1=29475, area2=29785, cross.area=19042, category=c("Oct4","SOX2"), fill=c("lightsalmon","moccasin") )
grid.draw(O_S.venn)

grid.newpage()
O_K.venn <- draw.pairwise.venn(area1=29475, area2=24981, cross.area=5441, category=c("Oct4","KLF4"), fill=c("lightsalmon","mistyrose4") )
grid.draw(O_K.venn)

grid.newpage()
S_K.venn <- draw.pairwise.venn(area1=29785, area2=24981, cross.area=6625, category=c("SOX2","KLF4"), fill=c("moccasin","mistyrose4") )
grid.draw(S_K.venn) 

grid.newpage()
O_M.venn <- draw.pairwise.venn(area1=29475, area2=3482, cross.area=409, category=c("Oct4","cMyc"), fill=c("lightsalmon","paleturquoise") )
grid.draw(O_M.venn)

grid.newpage()
S_M.venn <- draw.pairwise.venn(area1=29785, area2=3482, cross.area=451, category=c("SOX2","cMyc"), fill=c("moccasin","paleturquoise") )
grid.draw(S_M.venn)

grid.newpage()
K_M.venn <- draw.pairwise.venn(area1=24981, area2=3482, cross.area=681, category=c("KLF4","cMyc"), fill=c("mistyrose4","paleturquoise") )
grid.draw(K_M.venn)

grid.newpage()
O_S_K.venn <- draw.triple.venn(area1=29475, area2=29785, area3=24981, n12=19042, n23=6625, n13=5441, n123=4417, category=c("Oct4","SOX2","KLF4"), fill=c("lightsalmon","moccasin","mistyrose4") ) 
grid.draw(O_S_K.venn)

grid.newpage()
O_S_M.venn <- draw.triple.venn(area1=29475, area2=29785, area3=3482, n12=19042, n23=451, n13=409, n123=259, category=c("Oct4","SOX2","cMyc"), fill=c("lightsalmon","moccasin","paleturquoise") ) 
grid.draw(O_S_M.venn)

grid.newpage()
O_K_M.venn <- draw.triple.venn(area1=29475, area2=24981, area3=3482, n12=5441, n23=681, n13=409, n123=215, category=c("Oct4","KLF4","cMyc"), fill=c("lightsalmon","mistyrose4","paleturquoise") ) 
grid.draw(O_K_M.venn)

grid.newpage()
S_K_M.venn <- draw.triple.venn(area1=29785, area2=24981, area3=3482, n12=6625, n23=681, n13=451, n123=226, category=c("SOX2","KLF4","cMyc"), fill=c("moccasin","mistyrose4","paleturquoise") ) 
grid.draw(S_K_M.venn)

dev.off()
