########################################
uniCox_exp<-read.table(file = "uniCox_exp.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
uniCox_group<-read.table(file = "uniCox_group.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
exp<-uniCox_exp[which(uniCox_exp$pvalue<0.05),1]
group<-uniCox_group[which(uniCox_group$pvalue<0.05),1]
inter<-intersect(exp,group)
ERV_DE<-uniCox_exp[which(uniCox_exp[,1]%in%inter),]
ERV_type<-read.table(file = "ERV_type.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
data<-merge(ERV_DE,ERV_type,by.x=colnames(ERV_DE)[1],by.y=colnames(ERV_type)[1])
data_plot<-data.frame(name = data$id,Hazard.Ratio = data$HR,CI_low=data$HR.95L,CI_high = data$HR.95H,Type = data$V9)
col_1<-c("#E45D5A","#277698","#6DB9AC","#E2ABA0")
p1<-ggplot(data_plot,aes(x=Hazard.Ratio,y=name))+
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high,color=Type),width=0.5,size=0.3)+
  geom_point(pch=22,size=2,fill='white',aes(color=Type))+
  geom_point(pch=16,size=1,aes(color=Type))+
  geom_vline(xintercept = 1,color = 'red',size = 0.5,lty='dashed')+
  xlab(bquote(italic('Hazard Ratios')))+ylab(NULL)+
  ggtitle('Univariate Cox Analysis (OS)')+
  theme1+ scale_color_manual(values = col_1)
############################################
ERV_DE<-uniCox_group[which(uniCox_group[,1]%in%inter),]
ERV_type<-read.table(file = "ERV_type.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
data<-merge(ERV_DE,ERV_type,by.x=colnames(ERV_DE)[1],by.y=colnames(ERV_type)[1])
data_plot<-data.frame(name = data$id,Hazard.Ratio = data$HR,CI_low=data$HR.95L,CI_high = data$HR.95H,Type = data$V9)
col_1<-c("#E45D5A","#277698","#6DB9AC","#E2ABA0")
p2<-ggplot(data_plot,aes(x=Hazard.Ratio,y=name))+
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high,color=Type),width=0.5,size=0.3)+
  geom_point(pch=22,size=2,fill='white',aes(color=Type))+
  geom_point(pch=16,size=1,aes(color=Type))+
  geom_vline(xintercept = 1,color = 'red',size = 0.5,lty='dashed')+
  xlab(bquote(italic('Hazard Ratios')))+ylab(NULL)+
  ggtitle('Univariate Cox Analysis (OS)')+
  theme1+scale_color_manual(values = col_1)

library(aplot)
p2%>%insert_right(p1, width = 1)
###############################################
inter<-read.table(file = "gene.txt",header=FALSE,sep="\t",stringsAsFactors=FALSE)
ERV_type<-read.table(file = "ERV_type.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
tt<-merge(inter,ERV_type,by.x=colnames(inter)[1],by.y=colnames(ERV_type)[1])
library("RColorBrewer")
dd<-data.frame(table(tt[,2]))
per <- paste(round(100 * dd[,2] / sum(dd[,2]),2),"%")
pie(dd[,2],labels=per,col=brewer.pal(6,"Set3"))
legend("topright",as.character(dd[,1]),cex=1.2,border=TRUE,fill=brewer.pal(6,"Set3"),bty="n")
uniCox_exp<-read.table(file = "uniCox_exp.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
HR<-uniCox_exp[which(uniCox_exp[,1]%in%inter[,1]),2]
barplot(c(length(which(HR>1)),length(which(HR<1))),las=1,col = brewer.pal(n = 3, name = "RdBu")[1],
names=c("HR>1","HR<1"),ylab = "Number of HERV")
###############################################
library(corrplot) 
TPM<-read.table(file = "ESCC_TPM.txt",header = TRUE, sep="\t",row.names=1,stringsAsFactors = F)
df<-read.table(file = "gene.txt",header=FALSE,sep="\t",stringsAsFactors=F)
dd<-TPM[df[,1],1:137]
data<-cor(t(dd))
library(RColorBrewer)
corrplot(data, method = "square", type = "lower",col = rev(c("#2f3a7d","#4A7A9F","#76B6EA","#A2D5F2","#D4E6F1","#F9E79F","#F7DC6F","#F9CA3F","#F9A825")),
tl.col = "n", tl.cex = 0.8, tl.pos = "n",order = 'AOE',
add = T)
################################################
inter<-read.table(file = "gene.txt",header=FALSE,sep="\t",stringsAsFactors=FALSE)[,1]
Data<-read.table(file = "ERV_cor_down.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
gene<-unique(Data[which(Data[,1]%in%inter),2])
tt<-bitr(gene,fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")
BP<-enrichGO(tt[,2],OrgDb = "org.Hs.eg.db",ont = "BP",pvalueCutoff = 0.05,
qvalueCutoff = 0.05,keyType = "SYMBOL")
write.table(BP,file = "inter_de_cor_bp_down.txt",sep="\t",quote=F,row.names=F)
###################################################
GO_data<-read.table(file = "inter_de_cor_bp_down.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)[1:20,]
GO_data <- GO_data[order(GO_data$qvalue,decreasing = T),]
level <- GO_data[,2] 
GO_data$Description <- factor(GO_data$Description ,level = level) 
p1 <- ggplot(GO_data) +
  aes(x = Description, y = -log10(qvalue), fill = -log10(qvalue)) +
  geom_bar(stat = "identity") + 
  scale_fill_continuous(low="#D4E6F1", high="#76B6EA") +
  coord_flip() + 
  theme(
    axis.title = element_text(size=15,face="plain",color="black"),
    axis.text = element_text(size=12,face="plain",color="black"), 
    axis.text.y = element_text(hjust=1,vjust=0.6), 
    # axis.title.y = element_blank(), 
    legend.title = element_blank(), 
    legend.text = element_text(size = 8, face = "bold"), 
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
    legend.direction = "vertical", 
    # legend.position = c(0.9,0.92), 
    legend.background = element_blank(), 
    panel.background = element_rect(fill = "transparent",colour = "black"), 
    plot.background = element_blank()
  ) + 
  labs(x = "Biological Process");p1
###########################################  
library(tidyr) 
Data<-read.table(file = "qpcr_long.txt",header=TRUE,sep="\t")
ggplot(data = Data, aes(y = value, x = type, fill = type, color = type)) + 
    geom_line(aes(group = id), col = "gray") +
    geom_point(aes(shape = type), size = 2, 
    position = position_dodge(width = 0.2), show.legend = TRUE) +
    facet_wrap(~ gene, nrow = 1, strip.position = "top") +
    labs(y = "log2(2^-△CT*10^5)") +
    scale_color_manual(values=c("#F67280","#7FC7E2")) +
    scale_shape_manual(values = c(16, 17)) +
    theme_classic()+
	stat_compare_means(method = "t.test", paired = TRUE, label = "..p.signif..",label.y = 11,hjust = 0.5)+
	facet_wrap(~gene,nrow=2)+
	theme(axis.text.x = element_text(angle = 45, hjust = 1))
