##################################
library("limma")
library("clusterProfiler")
library("tibble")
library("forcats")
Data<-read.delim(file = "ESCC_TPM.txt",header=TRUE,sep="\t",stringsAsFactors=F)
tt<-bitr(rownames(Data),fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")
df<-merge(tt,Data,by.y="row.names",by.x = colnames(tt)[1])[,-1]
data<-avereps(df[,-1],ID=df[,1])
options(stringsAsFactors = F)

library(pheatmap)
label<-read.csv(file = "ERV_DEG_inter.k=3.consensusClass.csv",header= FALSE,stringsAsFactors=FALSE)
label<-label[order(label[,2]),]
data<-data[,label[order(label[,2]),1]]
data<-t(scale(t(data)))
library(IOBR)
library(tidyr)					   
xcell <- deconvo_tme(eset = data, method = "xcell", arrays = FALSE)
timer <- deconvo_tme(eset = data, method = "timer", group_list = rep("uvm", dim(data)[2]))

library(ggpubr)
library(stringr)
library(forcats)
timer_long <- timer %>% 
  select(ID,everything()) %>% 
  pivot_longer(-ID,names_to = "cell_type",values_to = "fraction") %>% 
  dplyr::mutate(cell_type = gsub("_TIMER","",cell_type),
                cell_type = gsub("_"," ",cell_type))
timer_type<-merge(label,timer,by.x = colnames(label)[1],by.y = colnames(timer)[1],sort=F)
timer_type[,2]<-factor(c(rep("HERV_1",59),rep("HERV_2",63),rep("HERV_3",15)))
colnames(timer_type)[1:2]<-c("ID","Group")
tt<-merge(timer_long,timer_type[,1:2],by.x = "ID",by.y="ID",sort=F)

p3 <- ggplot(tt,aes(fct_reorder(cell_type,fraction),fraction,fill = Group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  scale_fill_manual(values = palette1[c(2:4)])+ 
  theme_bw() + 
  labs(x = NULL, y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=45,hjust = 1),
        axis.text = element_text(color = "black",size = 12))+
  stat_compare_means(aes(group = Group,label = ..p.signif..),
                     method = "kruskal.test",label.y = 0.78)+
					 scale_fill_manual(values=c("MediumAquamarine","SandyBrown","Coral2")) 
##########################################################
xcell_long <- xcell %>% 
  select(ID,everything()) %>% 
  pivot_longer(-ID,names_to = "cell_type",values_to = "fraction") %>% 
  dplyr::mutate(cell_type = gsub("_xCell","",cell_type),
                cell_type = gsub("_"," ",cell_type))
xcell_type<-merge(label,xcell,by.x = colnames(label)[1],by.y = colnames(xcell)[1],sort=F)
xcell_type[,2]<-factor(c(rep("HERV_1",59),rep("HERV_2",63),rep("HERV_3",15)))
colnames(xcell_type)[1:2]<-c("ID","Group")
tt<-merge(xcell_long,xcell_type[,1:2],by.x = "ID",by.y="ID",sort=F)

p3 <- ggplot(tt,aes(fct_reorder(cell_type,fraction),fraction,fill = Group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  scale_fill_manual(values = c("#EECE5A","#27A9A1","#E45826"))+ 
  theme_bw() + 
  labs(x = NULL, y = "Infiltration Level") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=45,hjust = 1),
        axis.text = element_text(color = "black",size = 12))+
  stat_compare_means(aes(group = Group,label = ..p.signif..),
                     method = "kruskal.test",label.y = 2)
#############################################################


tmp<-data.frame(Group = xcell_type[,2],Tcell = xcell_type[,24])
g3<-ggplot(tmp, aes(Group, Tcell,  fill = Group )) + 
  geom_violin(trim = FALSE, position = position_dodge(0.9)) + 
  geom_boxplot(width = 0.2, position = position_dodge(0.9)) + 
  scale_fill_brewer(palette = "Set2") + 
  theme_bw(base_size = 15)+
  stat_compare_means(comparisons = my_comparisons,aes(group = Group,label = ..p.signif..),
                     method = "wilcox.test")+xlab("")+ylab("xcell score(Epithelial cells)")+
					 scale_fill_manual(values=c("MediumAquamarine","SandyBrown","Coral2")) 
					 
tmp<-data.frame(Group = xcell_type[,2],Tcell = xcell_type[,26])
g4<-ggplot(tmp, aes(Group, Tcell,  fill = Group )) + 
  geom_violin(trim = FALSE, position = position_dodge(0.9)) + 
  geom_boxplot(width = 0.2, position = position_dodge(0.9)) + 
  scale_fill_brewer(palette = "Set2") + 
  theme_bw(base_size = 15)+
  stat_compare_means(comparisons = my_comparisons,aes(group = Group,label = ..p.signif..),
                     method = "wilcox.test")+xlab("")+ylab("xcell score(Fibroblasts)")+
					 scale_fill_manual(values=c("MediumAquamarine","SandyBrown","Coral2")) 



tmp<-data.frame(Group = xcell_type[,2],Tcell = xcell_type[,60])
g7<-ggplot(tmp, aes(Group, Tcell,  fill = Group )) + 
  geom_violin(trim = FALSE, position = position_dodge(0.9)) + 
  geom_boxplot(width = 0.2, position = position_dodge(0.9)) + 
  scale_fill_brewer(palette = "Set2") + 
  theme_bw(base_size = 15)+
  stat_compare_means(comparisons = my_comparisons,aes(group = Group,label = ..p.signif..),
                     method = "wilcox.test")+xlab("")+ylab("xcell score(Skeletal muscle)")+
					 scale_fill_manual(values=c("MediumAquamarine","SandyBrown","Coral2")) 

tmp<-data.frame(Group = xcell_type[,2],Tcell = xcell_type[,25])
g8<-ggplot(tmp, aes(Group, Tcell,  fill = Group )) + 
  geom_violin(trim = FALSE, position = position_dodge(0.9)) + 
  geom_boxplot(width = 0.2, position = position_dodge(0.9)) + 
  scale_fill_brewer(palette = "Set2") + 
  theme_bw(base_size = 15)+
  stat_compare_means(comparisons = my_comparisons,aes(group = Group,label = ..p.signif..),
                     method = "wilcox.test")+xlab("")+ylab("xcell score(Erythrocytes)")+
					 scale_fill_manual(values=c("MediumAquamarine","SandyBrown","Coral2")) 

tmp<-data.frame(Group = xcell_type[,2],Tcell = xcell_type[,64])
g9<-ggplot(tmp, aes(Group, Tcell,  fill = Group )) + 
  geom_violin(trim = FALSE, position = position_dodge(0.9)) + 
  geom_boxplot(width = 0.2, position = position_dodge(0.9)) + 
  scale_fill_brewer(palette = "Set2") + 
  theme_bw(base_size = 15)+
  stat_compare_means(comparisons = my_comparisons,aes(group = Group,label = ..p.signif..),
                     method = "wilcox.test")+xlab("")+ylab("xcell score(Th2 cells)")+
					 scale_fill_manual(values=c("MediumAquamarine","SandyBrown","Coral2")) 

tmp<-data.frame(Group = xcell_type[,2],Tcell = xcell_type[,62])
g10<-ggplot(tmp, aes(Group, Tcell,  fill = Group )) + 
  geom_violin(trim = FALSE, position = position_dodge(0.9)) + 
  geom_boxplot(width = 0.2, position = position_dodge(0.9)) + 
  scale_fill_brewer(palette = "Set2") + 
  theme_bw(base_size = 15)+
  stat_compare_means(comparisons = my_comparisons,aes(group = Group,label = ..p.signif..),
                     method = "wilcox.test")+xlab("")+ylab("xcell score(Tgd cells)")+
                     scale_fill_manual(values=c("MediumAquamarine","SandyBrown","Coral2")) 					 
library("cowplot")
ggarrange(g3,g4,g7,g8,g9,g10,nrow=2,ncol = 3)					 

##########################################################
library("limma")
library("clusterProfiler")
library("tibble")
library("forcats")
dd<-read.table(file = "GSVA_KEGG_metabolic.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
res<-data.frame()
for(i in 1:dim(dd)[1]){
   HERV1<-as.numeric(dd[i,1:59])
   HERV2<-as.numeric(dd[i,60:122])
   HERV3<-as.numeric(dd[i,123:137])
   res[i,1]<-rownames(dd)[i]
   res[i,2]<-median(HERV1)
   res[i,3]<-median(HERV2)
   res[i,4]<-median(HERV3)
   res[i,5]<-wilcox.test(HERV1,HERV2)$p.value
   res[i,6]<-wilcox.test(HERV1,HERV3)$p.value
   res[i,7]<-wilcox.test(HERV2,HERV3)$p.value     
}
path<-res[which(res[,5]<0.05 & res[,6]<0.05 & res[,7]<0.05),1]
path<-read.table(file = "path_sig_KEGG.txt",header=FALSE,sep="\t",stringsAsFactors=FALSE)
ff<-dd[path[,1],]
library(pheatmap)
label<-read.csv(file = "ERV_cox_DE.k=3.consensusClass.csv",header= FALSE,stringsAsFactors=FALSE)
label<-label[order(label[,2]),]
dd<-ff[,label[order(label[,2]),1]]
library("pheatmap")
tt_1<-dd[,c(1:59)]
P1<-pheatmap(tt_1)
h_1<-colnames(tt_1)[P1$tree_col[["order"]]]
tt_2<-dd[,c(60:122)]
P2<-pheatmap(tt_2)
h_2<-colnames(tt_2)[P2$tree_col[["order"]]]
tt_3<-dd[,c(123:137)]
P3<-pheatmap(tt_3)
h_3<-colnames(tt_3)[P3$tree_col[["order"]]]
dd<-dd[,c(h_1,h_2,h_3)]

tt_1<-dd[7:13,]
tt_2<-dd[14:22,]
P1<-pheatmap(tt_1)
h_1<-rownames(tt_1)[P1$tree_row[["order"]]]
P2<-pheatmap(tt_2)
h_2<-rownames(tt_2)[P2$tree_row[["order"]]]
dd<-dd[c(rownames(dd)[1:6],h_1,h_2,rownames(dd)[23:28]),]
test<-dd
annotation_col = data.frame(Sample=factor(c(rep("HERV_1",59),rep("HERV_2",63),rep("HERV_3",15))))
annotation_row = data.frame(
                    PathwayClass = factor(path[,2]))
row.names(annotation_col) = colnames(test)
rownames(annotation_row) = rownames(test)
col<-brewer.pal(n = 10, "Paired")[c(2,4,6)]
ann_colors = list(Sample = c(HERV_1="MediumAquamarine",HERV_2="SkyBlue",HERV_3="Coral2"),
                  PathwayClass = brewer.pal(n = 6, "Accent"))
bk <- c(seq(-1,-0.1,by=0.01),seq(0,1,by=0.01))
P<-pheatmap(test, 
         show_colnames =F,
		 annotation_col = annotation_col,
		 annotation_row = annotation_row,
		 ann_colors = ann_colors,
	     cluster_cols = F,
		 cluster_rows = F,
		 breaks = bk,
		 col = c(colorRampPalette(colors = c("#4D9221","white"))(length(bk)*0.50),colorRampPalette(colors = c("white","#C51B7D"))(length(bk)*0.5)),
		 clustering_distance_rows = "euclidean"
		 )

