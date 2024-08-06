############################################
df<-read.table(file = "ERV_cox_DE.txt",header=FALSE,sep="\t",stringsAsFactors=F)
library(ConsensusClusterPlus)
results = ConsensusClusterPlus(as.matrix(df),maxK=6,reps=1000,pItem=0.8,pFeature=1,
                               title="ERV_cox_DE",clusterAlg="km",distance="euclidean",plot="pdf",writeTable = TRUE)


df<-read.table(file = "ERV_cox_DE.txt",header = TRUE, sep="\t",row.names=1,stringsAsFactors = F)
label<-read.csv(file = "ERV_cox_DE.k=3.consensusClass.csv",header= FALSE,stringsAsFactors=FALSE)
label<-label[order(label[,2]),]
dd<-dd[,label[order(label[,2]),1]]

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

test<-t(scale(t(dd)))
annotation_col = data.frame(Sample=factor(c(rep("HERV_1",59),rep("HERV_2",63),rep("HERV_3",15))))
row.names(annotation_col) = colnames(test)
col<-brewer.pal(n = 10, "Paired")[c(2,4,6)]
ann_colors = list(Sample = c(HERV_1=col[1],HERV_2=col[2],HERV_3=col[3]))
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
P<-pheatmap(test, 
         show_rownames = F,
         show_colnames =F,
		 annotation_col = annotation_col,
		 annotation_colors = list(Sample= c("HERV_1" = "#EECE5A", "HERV_2" = "#27A9A1","HERV_3" = "#E45826")),
	     cluster_cols = F,
		 breaks = bk,
		 col = c(colorRampPalette(colors = c("#7FC7E2","white"))(length(bk)*0.50),colorRampPalette(colors = c("white","#F67280"))(length(bk)*0.5)),
		 clustering_distance_rows = "euclidean"
		 )
		 

library("FactoMineR")
library("factoextra") 	 
tt<-t(test)
tt<-data.frame(tt,group = c(rep("HERV_1",59),rep("HERV_2",63),rep("HERV_3",15)))
res.pca <- PCA(tt[,-dim(tt)[2]], graph = FALSE)
fviz_pca_ind(res.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = tt$group, # color by groups
             palette = c("#EECE5A","#27A9A1","#E45826"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)							   
							   

###################################################
TPM<-read.table(file = "TPM.txt",header = TRUE, sep="\t",row.names=1,stringsAsFactors = F)
df<-read.table(file = "gene.txt",header=TRUE,sep="\t",stringsAsFactors=F)
dd<-TPM[which(rownames(TPM)%in%df[,1]),1:23]
label<-read.csv(file = "ERV_DEG_inter.k=3.consensusClass.csv",header= FALSE,stringsAsFactors=FALSE)
label<-label[order(label[,2]),]
dd<-dd[,label[order(label[,2]),1]]

library("pheatmap")
tt_1<-dd[,c(1:11)]
P1<-pheatmap(tt_1)
h_1<-colnames(tt_1)[P1$tree_col[["order"]]]
tt_2<-dd[,c(12:20)]
P2<-pheatmap(tt_2)
h_2<-colnames(tt_2)[P2$tree_col[["order"]]]
tt_3<-dd[,c(21:23)]
P3<-pheatmap(tt_3)
h_3<-colnames(tt_3)[P3$tree_col[["order"]]]
dd<-dd[,c(h_2,h_1,h_3)]

test<-t(scale(t(dd)))
annotation_col = data.frame(Sample=factor(c(rep("HERV_1",9),rep("HERV_2",11),rep("HERV_3",3))))
row.names(annotation_col) = colnames(test)
col<-brewer.pal(n = 10, "Paired")[c(2,4,6)]
ann_colors = list(Sample = c(HERV_1=col[1],HERV_2=col[2],HERV_3=col[3]))
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
P<-pheatmap(test, 
         show_rownames = F,
         show_colnames =F,
		 annotation_col = annotation_col,
		 annotation_colors = list(Sample= c("HERV_1" = "#EECE5A", "HERV_2" = "#27A9A1","HERV_3" = "#E45826")),
	     cluster_cols = F,
		 breaks = bk,
		 col = c(colorRampPalette(colors = c("#7FC7E2","white"))(length(bk)*0.50),colorRampPalette(colors = c("white","#F67280"))(length(bk)*0.5)),
		 clustering_distance_rows = "euclidean"
		 )

library("FactoMineR")
library("factoextra") 	 
tt<-t(test)
#annotation_col = data.frame(Sample=factor(c(rep("HERV_1",92),rep("HERV_2",49)))
tt<-data.frame(tt,group = c(rep("HERV_1",9),rep("HERV_2",11),rep("HERV_3",3)))
res.pca <- PCA(tt[,-dim(tt)[2]], graph = FALSE)
fviz_pca_ind(res.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = tt$group, # color by groups
             palette = c("#EECE5A","#27A9A1","#E45826"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)

####################################################
library(tidyr) 
Data<-read.table(file = "qpcr_long.txt",header=TRUE,sep="\t")
dd<-Data[which(Data[,2]%in%"Tumor"),c(1,3,4)]
data <- pivot_wider(dd, id_cols = 1, names_from = "id", values_from = value)
df<-data.frame(data)
rownames(df)<-df[,1]
df<-df[,-1]
tt<-t(scale(t(df)))
library(ConsensusClusterPlus)
results = ConsensusClusterPlus(as.matrix(tt),maxK=6,reps=1000,pItem=0.8,pFeature=1,
                               title="qpcr",clusterAlg="km",distance="euclidean",plot="pdf",writeTable = TRUE)
							   
label<-read.csv(file = "qpcr.k=3.consensusClass.csv",header= FALSE,stringsAsFactors=FALSE)
label<-label[order(label[,2]),]
dd<-tt[,label[order(label[,2]),1]]

library("pheatmap")
tt_1<-dd[,c(1:15)]
P1<-pheatmap(tt_1)
h_1<-colnames(tt_1)[P1$tree_col[["order"]]]
tt_3<-dd[,c(17:20)]
P3<-pheatmap(tt_3)
h_3<-colnames(tt_3)[P3$tree_col[["order"]]]
dd<-dd[,c(h_3,h_1,"X19")]

test<-dd
annotation_col = data.frame(Sample=factor(c(rep("HERV_1",4),rep("HERV_2",15),rep("HERV_3",1))))
row.names(annotation_col) = colnames(test)
col<-brewer.pal(n = 10, "Paired")[c(2,4,6)]
ann_colors = list(Sample = c(HERV_1=col[1],HERV_2=col[2],HERV_3=col[3]))
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
P<-pheatmap(test, 
         show_rownames = F,
         show_colnames =F,
		 annotation_col = annotation_col,
		 annotation_colors = list(Sample= c("HERV_1" = "#EECE5A", "HERV_2" = "#27A9A1","HERV_3" = "#E45826")),
	     cluster_cols = F,
		 breaks = bk,
		 col = c(colorRampPalette(colors = c("#7FC7E2","white"))(length(bk)*0.50),colorRampPalette(colors = c("white","#F67280"))(length(bk)*0.5)),
		 clustering_distance_rows = "euclidean"
		 )
		 
library("FactoMineR")
library("factoextra") 	 
tt<-t(test)
#annotation_col = data.frame(Sample=factor(c(rep("HERV_1",92),rep("HERV_2",49)))
tt<-data.frame(tt,group = c(rep("HERV_2",15),rep("HERV_3",1),rep("HERV_1",4)))
res.pca <- PCA(tt[,-dim(tt)[2]], graph = FALSE)
fviz_pca_ind(res.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = tt$group, # color by groups
             palette = c("#EECE5A","#27A9A1","#E45826"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)
	
############################################################
library(survival)
library(survminer)
library(RColorBrewer)
library(tibble)
library(ggpp)
Data<-read.table(file = "clin_new.txt",sep="\t",row.names=1,check.names=F,stringsAsFactors = F,header = TRUE)
dat<-data.frame(OS.time=Data[,25],OS=Data[,23],group=Data[,3])
dat<-dat[-which(dat[,2]%in%"loss"),]
dat[,1]<-dat[,1]/30.4
dat[,2]<-as.numeric(dat[,2])
fitd <- survdiff(Surv(OS.time, OS) ~ group,
                 data      = dat,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(OS.time, OS)~ group,
               data      = dat,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)

ps <- pairwise_survdiff(Surv(OS.time, OS)~ group,
                        data            = dat,
                        p.adjust.method = "none") 

mycol <- brewer.pal(n = 10, "Paired")[c(2,4,6,8)]
names(fit$strata) <- gsub("group=", "", names(fit$strata))
p <- ggsurvplot(fit               = fit,
                conf.int          = FALSE, 
                risk.table        = TRUE, 
                risk.table.col    = "strata",
                palette           = mycol, 
                data              = dat,
                #xlim              = c(0,60),
                size              = 1,
                break.time.by     = 12, 
                legend.title      = "",
                xlab              = "Time (months)",
                ylab              = "Overall survival", #Overall survival Disease free survival
                risk.table.y.text = FALSE,
                tables.height     = 0.3) 

p.lab <- paste0("log-rank test P",
                ifelse(p.val < 0.001, " < 0.001", 
                       paste0(" = ",round(p.val, 3))))

p$plot <- p$plot + annotate("text",
                            x = 0, y = 0.55,
                            hjust = 0,
                            fontface = 4,
                            label = p.lab)

addTab <- as.data.frame(as.matrix(ifelse(round(ps$p.value, 3) < 0.001, "<0.001",
                                         round(ps$p.value, 3))))
addTab[is.na(addTab)] <- "-"
df <- tibble(x = 0, y = 0, tb = list(addTab))
p$plot <- p$plot + 
  geom_table(data = df, 
             aes(x = x, y = y, label = tb), 
             table.rownames = TRUE)

pdf(".pdf", width = 4.5, height = 6)

dev.off()	 
############################################
library(survival)
library(survminer)
library(RColorBrewer)
library(tibble)
library(ggpp)
Data<-read.table(file = "clin_new.txt",sep="\t",row.names=1,check.names=F,stringsAsFactors = F,header = TRUE)
#dat<-data.frame(OS.time=Data[,4]/30.4,OS=Data[,2],group=Data[,1])
dat<-data.frame(OS.time=Data[,26],OS=Data[,24],group=Data[,3])
dat<-dat[-which(dat[,2]%in%"loss"),]
dat[,1]<-dat[,1]/30.4
dat[,2]<-as.numeric(dat[,2])
fitd <- survdiff(Surv(OS.time, OS) ~ group,
                 data      = dat,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(OS.time, OS)~ group,
               data      = dat,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)

ps <- pairwise_survdiff(Surv(OS.time, OS)~ group,
                        data            = dat,
                        p.adjust.method = "none") 

mycol <- brewer.pal(n = 10, "Paired")[c(2,4,6,8)]
names(fit$strata) <- gsub("group=", "", names(fit$strata))
p <- ggsurvplot(fit               = fit,
                conf.int          = FALSE,
                risk.table        = TRUE,
                risk.table.col    = "strata",
                palette           = mycol,
                data              = dat,
                #xlim              = c(0,60), 
                size              = 1,
                break.time.by     = 12, 
                legend.title      = "",
                xlab              = "Time (months)",
                ylab              = "Disease free survival", #Overall survival Disease free survival
                risk.table.y.text = FALSE,
                tables.height     = 0.3) 

p.lab <- paste0("log-rank test P",
                ifelse(p.val < 0.001, " < 0.001", c
                       paste0(" = ",round(p.val, 3))))

p$plot <- p$plot + annotate("text",
                            x = 0, y = 0.55, 
                            hjust = 0,
                            fontface = 4,
                            label = p.lab)

addTab <- as.data.frame(as.matrix(ifelse(round(ps$p.value, 3) < 0.001, "<0.001",
                                         round(ps$p.value, 3))))
addTab[is.na(addTab)] <- "-"
df <- tibble(x = 0, y = 0, tb = list(addTab))
p$plot <- p$plot + 
  geom_table(data = df, 
             aes(x = x, y = y, label = tb), 
             table.rownames = TRUE)

pdf(".pdf", width = 4.5, height = 6)

dev.off()	 
