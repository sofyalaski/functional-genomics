df_mef <- read.table("MEF.txt")
df_esc <- read.table("ESC.txt")
names <- c("peaks", "readfile", "modification","mus","musculus","method", "code")
colnames(df_esc)<- names
colnames(df_mef)<- names

df_mef$HM <- factor(df_mef$code == "HM",levels = c(TRUE, FALSE),labels = c("HM", "TF"))
g1<- ggplot(data=as.data.frame(df_mef$peaks[1:20]), aes(x=seq(1:20), y=df_mef$peaks[1:20], fill = df_mef$HM[1:20])) +
  geom_bar(stat="identity", width = 0.9)+
  scale_y_continuous(name="Number of peaks")+
  scale_x_continuous(name="Modification")+
  ggtitle("Number of Peaks in MEF")+
  geom_vline(xintercept = 11.5, col = "red")+ 
  scale_fill_manual(name = "", values = c("HM" = "lightcyan1", "TF" = "royalblue1")) +
  geom_text(aes(label = df_mef$modification[1:20]),angle = 0, size = 3)


df_esc$HM <- factor(df_esc$code == "HM",levels = c(TRUE, FALSE),labels = c("HM", "TF"))
g2<- ggplot(data=as.data.frame(df_esc$peaks[1:20]), aes(x=seq(1:20), y=df_esc$peaks[1:20], fill = df_esc$HM[1:20])) +
  geom_bar(stat="identity", width = 0.9)+
  scale_y_continuous(name="Number of peaks")+
  scale_x_continuous(name="Modification")+
  ggtitle("Number of Peaks in ESC")+
  geom_vline(xintercept = 11.5, col = "red")+ 
  scale_fill_manual(name = "", values = c("HM" = "lightcyan1", "TF" = "royalblue1")) +
  geom_text(aes(label = df_esc$modification[1:20]),angle = 0, size = 3)


grid.arrange(g1, g2, nrow = 2)

