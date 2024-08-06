TPM<-read.table(file = "ESCC_TPM_ERV.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE)
DD<-cor(TPM)
annotation_row <- data.frame(Group = factor(rep(c("Tumor", "Tumor-adjacent"), each = 137)))
annotation_col <- data.frame(Group = factor(rep(c("Tumor", "Tumor-adjacent"), each = 137)))
row.names(annotation_col) = colnames(DD)
row.names(annotation_row) = colnames(DD)
ann_colors_row <- list(Group = c("Tumor" = "#F67280", "Tumor-adjacent" = "#7FC7E2"))
ann_colors_col <- list(Group = c("Tumor" = "#F67280", "Tumor-adjacent" = "#7FC7E2"))
#bk <- c(seq(-8,-0.1,by=0.01),seq(0,8,by=0.01))
library(pheatmap)
tt_1<-DD[1:137,1:137]
p1<-pheatmap(tt_1)
t_n<-colnames(tt_1)[p1$tree_col[["order"]]]
tt_2<-DD[138:274,138:274]
p2<-pheatmap(tt_2)
n_n<-colnames(tt_2)[p2$tree_col[["order"]]]
DD<-cor(TPM[,c(t_n,n_n)])
P<-pheatmap(DD, 
         show_rownames = F,
         show_colnames =F,
		 annotation_col = annotation_col,
		 annotation_row = annotation_row,
		 annotation_colors = list(Group = c("Tumor" = "#F67280", "Tumor-adjacent" = "#7FC7E2")),
	     cluster_cols = F,
		 cluster_rows = F
         )	
###########################PCA		 
library("FactoMineR")
library("factoextra")
TPM<-read.table(file = "ESCC_TPM_ERV.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE)	 	 
tt<-t(TPM)
#annotation_col = data.frame(Sample=factor(c(rep("HERV_1",92),rep("HERV_2",49)))
tt<-data.frame(tt,group = c(rep("Tumor",137),rep("Tumor-adjacent",137)))
res.pca <- PCA(tt[,-dim(tt)[2]], graph = FALSE)
fviz_pca_ind(res.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = tt$group, # color by groups
             palette = c("#F67280","#7FC7E2"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)+ theme_minimal()
#########################
Data<-read.delim(file = "ESCC_TPM.txt",header=TRUE,sep="\t",stringsAsFactors=F)
DD<-log2(Data+1)
ERV<-DD[grep('^G',rownames(Data)),]
Tumor<-ERV[,1:137]
Normal<-ERV[,138:274]
Tumor_mean<-apply(Tumor,1,mean)
Normal_mean<-apply(Normal,1,mean)
Data<-read.table(file = "ERV_merge_HERV_uniq.bed",header=FALSE,sep="\t",stringsAsFactors=FALSE)[,c(4,8,9)]
Tumor_df<-data.frame(ID = names(Tumor_mean),value = Tumor_mean,type = "Tumor" )
Normal_df<-data.frame(ID = names(Normal_mean),value = Normal_mean, type = "Tumor-adjacent")
df_tumor<-merge(Data,Tumor_df,by.x = colnames(Data)[1],by.y = colnames(Tumor_df)[1])
colnames(df_tumor)<-c("ID","sub_famlily","famlily","expression","type")
df_normal<-merge(Data,Normal_df,by.x = colnames(Data)[1],by.y = colnames(Normal_df)[1])
colnames(df_normal)<-c("ID","sub_famlily","famlily","expression","type")


library(tidyverse) 
library(cowplot) 
library(smplot2)
g1<-ggplot(mapping = aes(x = famlily, y = expression, fill = type)) +
  sm_raincloud(data = df_tumor, position = position_nudge(x = +0.15),
               show.legend = FALSE,
			   boxplot.params = list(outlier.color =NA ),
               point.params = list(size = 2, shape = 21,
                                   color = 'transparent', 
                                   show.legend = TRUE,
                                   position = sdamr::position_jitternudge(nudge.x = 0.06, seed = 10,
                                                                   jitter.width = 0.06))) +
  scale_fill_manual(values = c("#F67280","#7FC7E2")) +
  sm_raincloud(data = df_normal, which_side = 'left',
               show.legend = FALSE,
			   boxplot.params = list(outlier.color =NA),
               position = position_nudge(x = -0.15),
               point.params = list(size = 2, shape = 21,
                                   show.legend = TRUE,
                                   color = 'transparent',
                                   position = sdamr::position_jitternudge(nudge.x = -0.06, seed = 10,
                                                                   jitter.width = 0.06)))
#########################################################################
library(RColorBrewer)
DE<-read.csv(file = "ERV_DESeq2.csv",header=TRUE,stringsAsFactors=F)
DE_up<-DE[which(DE$log2FoldChange>1&DE$padj<0.05),1]
DE_down<-DE[which((DE$log2FoldChange)<(-1)&DE$padj<0.05),1]

Data<-read.table(file = "ERV_merge_HERV_uniq.bed",header=FALSE,sep="\t",stringsAsFactors=FALSE)
ERV_type<-unique(Data[,c(4,9)])
dd<-data.frame(table(ERV_type[which(ERV_type[,1]%in%DE_up),2]))
per <- paste(round(100 * dd[,2] / sum(dd[,2]),2),"%")
pie(dd[,2],labels=per,col=brewer.pal(6,"Set3"))
legend("topright",as.character(dd[,1]),cex=1.2,border=TRUE,fill=brewer.pal(6,"Set3"),bty="n")

dd<-data.frame(table(ERV_type[which(ERV_type[,1]%in%DE_down),2]))
per <- paste(round(100 * dd[,2] / sum(dd[,2]),2),"%")
pie(dd[,2],labels=per,col=brewer.pal(6,"Set3"))
legend("topright",as.character(dd[,1]),cex=1.2,border=TRUE,fill=brewer.pal(6,"Set3"),bty="n")
  
##################################################################
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
df<-read.csv(file = "ERV_DEseq2.csv",header=TRUE,stringsAsFactors=F)
pvalue = 0.05
log2FC = 1
df$group <- case_when(
df$log2FoldChange > log2FC & df$padj < pvalue ~ "Up",
df$log2FoldChange < -log2FC & df$padj < pvalue ~ "Down",
TRUE ~ 'None'
)
df$'-log10(FDR)' <- -log10(df$padj) #新增-log10p列
df$group <- factor(df$group, levels = c("Up","Down","None"))
#ggplot2建立映射：
mycol <- c("#EB4232","#2DB2EB","#d8d8d8")
mytheme <- theme_classic() +
theme(
panel.grid = element_blank(),
axis.title = element_text(size = 15),
axis.text = element_text(size = 14),
legend.text = element_text(size = 14)
)
p <- ggplot(data = df,
aes(x = log2FoldChange, y = -log10(padj), color = group)) + #建立映射
geom_point(size = 2.2)+
ylab(label = "-log(FDR)")+ #绘制散点
scale_colour_manual(name = "", values = alpha(mycol, 0.7)) +
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = 'black',size=12),
        axis.title = element_text(color = 'black',size=15),
		legend.text = element_text(size = 14),
        legend.position = "right")
#####################################################
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
df<-read.csv(file = "DE_sub.csv")
pvalue = 0.05
log2FC = 1
df$group <- case_when(
df$log2FC > log2FC & df$FDR < pvalue ~ "Up",
df$log2FC < -log2FC & df$FDR < pvalue ~ "Down",
TRUE ~ 'None'
)
df$'-log10(FDR)' <- -log10(df$FDR) #新增-log10p列
df$group <- factor(df$group, levels = c("Up","Down","None"))
#ggplot2建立映射：
mycol <- c("#EB4232","#2DB2EB","#d8d8d8")
mytheme <- theme_classic() +
theme(
panel.grid = element_blank(),
axis.title = element_text(size = 15),
axis.text = element_text(size = 14),
legend.text = element_text(size = 14)
)
p <- ggplot(data = df,
aes(x = log2FC, y = -log10(FDR), color = group)) + #建立映射
geom_point(size = 2.2)+
ylab(label = "-log(FDR)")+ #绘制散点
scale_colour_manual(name = "", values = alpha(mycol, 0.7)) +
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = 'black',size=12),
        axis.title = element_text(color = 'black',size=15),
		legend.text = element_text(size = 14),
        legend.position = "right")		

de<-df[which(df$FDR<0.05&abs(df$FD)>1),]
data<-tt[which(tt[,1]%in%de[,1]),]
rownames(data)<-data[,1]
dd<-data[,-1]
dd<-log2(dd+0.01)
DD<-t(scale(t(dd)))
annotation_col <- data.frame(Group = factor(rep(c("Tumor", "Tumor-adjacent"), each = 137)))
row.names(annotation_col) = colnames(dd)
ann_colors_col <- list(Group = c("Tumor" = "#F67280", "Tumor-adjacent" = "#7FC7E2"))
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
P<-pheatmap(DD, 
         show_rownames = T,
         show_colnames =F,
		 annotation_col = annotation_col,
		 annotation_colors = list(Group = c("Tumor" = "#F67280", "Tumor-adjacent" = "#7FC7E2")),
	     cluster_cols = F,
		 cluster_rows = T,
		 border_color = F,
		 breaks=bk,
		 col = c(colorRampPalette(colors = c("#7FC7E2","white"))(length(bk)*0.50),colorRampPalette(colors = c("white","#F67280"))(length(bk)*0.5))
		 #legend_breaks=seq(-14,14,2)
         )	
##########################################################
up<-read.table(file = "cor_bp_up.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
enrichmentNetwork(up[1:100,], colorBy = 'pvalue', colorType = 'pval',pCutoff = -30)
down<-read.table(file = "cor_bp_down.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
enrichmentNetwork(down[1:100,], colorBy = 'pvalue', colorType = 'pval',pCutoff = -10)
		 
