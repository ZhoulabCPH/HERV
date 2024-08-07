#############################
library(ggplot2)
ERV_num<-read.table(file = "ERV_num_bar.txt",header=FALSE,sep="\t",stringsAsFactors=FALSE)
df$Group<-factor(df$Group,levels=c("Tumor", "Tumor-adjacent"),ordered = TRUE)
df$id<-factor(df$id,levels = index)
  ggplot()+
  geom_bar(data = df, stat="identity",position="dodge", aes( y = ERV_num,x =id, fill = Group))+
  labs(x=NULL,y="HERV number")+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = c("#F1A99A","#78A8CE"),name=NULL)+ 
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom"
        )
		
##############################
library(ggpubr)
library(stringr)
library(RColorBrewer)
dd<-data.frame(Number=df[,2],Group = df[,3],id =df[,4])
my_comparisons <- list( c("Tumor", "Tumor-adjacent") )
g3<-ggpaired(dd, x = 'Group', y = "Number",id = "id",palette=c("#F1A99A","#78A8CE"),
          color = "black",fill = 'Group', point.size=1,linetype=2, point.col = 'Group',
         line.color = "gray", line.size = 0.4,
         short.panel.labs = FALSE)+
		 labs(y = "HERV Number")+
  stat_compare_means(aes(label = ..p.signif..),
                     method = "t.test",paired = TRUE, comparisons = my_comparisons)
					 
ERV_dep<-read.table(file = "ERV_dep.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
library(cowplot)
df_1<-merge(df,ERV_dep,by.x=colnames(df)[1],by.y=colnames(ERV_dep)[1])
dd<-data.frame(Depth=df_1[,5],Group = df_1[,3],id = df_1[,4])
my_comparisons <- list( c("Tumor", "Tumor-adjacent") )
g4<-ggpaired(dd, x = 'Group', y = 'Depth',id = "id",palette=c("#F1A99A","#78A8CE"),
          color = "black",fill = 'Group',  point.size=1,linetype=2, point.col = 'Group',
         line.color = "gray", line.size = 0.4,
         short.panel.labs = FALSE)+
		 labs(x = NULL, y = "Depth") +
  stat_compare_means(aes(label = ..p.signif..),
                     method = "t.test",paired = TRUE, comparisons = my_comparisons)
					 
plot_grid(g3,g4,ncol=2)
#####################################
dd<-read.table(file = "circos_data.txt",header=TRUE,sep="\t")
library(dplyr)
dd <- dd %>%
  mutate(famlily_label = case_when(
    famlily == "ERV1" ~ 1,
    famlily == "ERVK" ~ 2,
    famlily == "ERVL" ~ 3,
	famlily == "ERVL-MaLR" ~ 4,
	famlily == "Gypsy" ~ 5,
	famlily == "LTR" ~ 6
  ))

tt<-dd[,c(1:3,5)]
df<-dd[,c(1:3,8)]
col_fun = colorRamp2(c(1:6), c("#F25F5C","#FFE066","#247BA0","#70C1B3","#527C5A","#E9AFA3"))
col = col_fun(df[,4])
library(circlize)
circos.initializeWithIdeogram(species = "hg38",plotType = NULL)
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  circos.rect(xlim[1], 0, xlim[2], 1, col = rand_color(1))
  circos.text(mean(xlim), mean(ylim), chr, cex = 0.7, col = "white",
              facing = "inside", niceFacing = TRUE)
}, track.height = 0.15, bg.border = NA)
#circos.genomicRainfall(dd,col = "#F25F5C",pch = 16,cex = 0.5)
circos.genomicDensity(tt,border = NA,col = "#F25F5C")
circos.genomicTrack(tt, panel.fun = function(region,value,...){
                      circos.genomicPoints(region,value,numeric.column = 1,pch = 16,col = "#FFE066",cex = 0.5)
                    })
#########################################################
dd<-read.table(file = "rsmk_LTE.txt",header=FALSE,sep="\t",stringsAsFactors=FALSE)
per <- paste(round(100 * dd[,1] / sum(dd[,1]),2),"%")
pie(dd[,1],labels=per,col=brewer.pal(6,"Set3"))
legend("topright",as.character(dd[,2]),cex=1.2,border=TRUE,fill=brewer.pal(6,"Set3"),bty="n")

#######################################
Data<-read.delim(file = "ERV_len.txt",header=TRUE,sep="\t")
type<-read.table(file = "ERV_merge_HERV_uniq.bed",header=FALSE,sep="\t",stringsAsFactors=FALSE)[,c(4,8,9)]
df<-merge(Data,type,by.x=colnames(Data)[1],by.y=colnames(type)[1])
colnames(df)<-c("ID","length","sub_famlily","famlily")
ggplot(df, aes(x = length, fill = famlily)) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = c("#F25F5C","#FFE066","#247BA0","#70C1B3","#527C5A","#E9AFA3"))+
   theme_bw()+scale_x_continuous(limits = c(0,5000))+
  theme(
        legend.position = "right",
		panel.grid = element_blank(),
		panel.background = element_blank()
        )
#######################################
Data<-read.delim(file = "ESCC_TPM.txt",header=TRUE,sep="\t",stringsAsFactors=F)
DD<-log2(Data+1)
ERV<-DD[grep('^G',rownames(Data)),]
other<-DD[-grep('^G',rownames(Data)),]
gene_type<-read.table(file = "gene_type.txt",header=FALSE,sep="\t",stringsAsFactors=FALSE)
#tt<-gene_type[which(rownames(ERV)%in%gene_type[,1]),]
loca<-which(gene_type[,2]%in%"protein_coding")
pd<-gene_type[loca,1]
loca<-which(rownames(other)%in%pd)
dat_pd<-other[loca,]
dat_ncRNA<-other[-loca,]
ERV_tumor<-apply(ERV[,1:137],1,median)
ERV_normal<-apply(ERV[,138:274],1,median)
pd_tumor<-apply(dat_pd[,1:137],1,median)
pd_normal<-apply(dat_pd[,138:274],1,median)
ncRNA_tumor<-apply(dat_ncRNA[,1:137],1,median)
ncRNA_normal<-apply(dat_ncRNA[,138:274],1,median)
expression<-c(ERV_tumor,ERV_tumor,pd_tumor,pd_normal,ncRNA_tumor,ncRNA_normal)
label<-c(rep("HERV_tumor",length(ERV_tumor)),rep("HERV_tumor-adjacent",length(ERV_normal)),rep("Pd_tumor",length(pd_tumor)),rep("Pd_tumor-adjacent",length(pd_normal)),rep("ncRNA_tumor",length(ncRNA_tumor)),rep("ncRNA_tumor-adjacent",length(ncRNA_normal)))
dat<-data.frame(Expression = expression,Group = label)
library(IOBR)
library(ggplot2)
library(ggridges)
ggplot(dat, aes(x = Expression, y = Group)) +
  geom_density_ridges(aes(fill = Group)) +
  scale_fill_manual(values = palette1)					 
