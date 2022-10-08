#### input file ####
rm(list=ls())
options(stringsAsFactors = F)

setwd("yuanwang/pan-cancer/ICD")

immuneCell=readRDS("Rdata/SingleCell/immune.rds")
Tcell=readRDS("Rdata/SingleCell/T-ICD.rds")
write.table(as.matrix(immuneCell@assays$RNA@data), 'data/SingleCell/cellphonedb_count.txt', sep='\t', quote=F)

df=immuneCell@meta.data
df$cellType=df$immune
df.T=Tcell@meta.data
for (i in 1:nrow(df)) {
  if (rownames(df)[i] %in% rownames(df.T))
  {df$cellType[i]=df.T[rownames(df)[i],'ICDgroup']}
  else {df$cellType[i]=df$cellType[i]}    
}
meta_data <- cbind(rownames(df), df[,'cellType', drop=F])  
colnames(meta_data)=c('Cell','cell_type')
write.table(meta_data, 'data/SingleCell/cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)

#### plot ####
rm(list = ls())

mypvals <- read.table("out/pvalues.txt",header = T,sep = "\t",stringsAsFactors = F)
mymeans <- read.table("out/means.txt",header = T,sep = "\t",stringsAsFactors = F) 

kp = grepl(pattern = "high.ICD", colnames(mypvals)) | grepl(pattern = "low.ICD", colnames(mypvals))
table(kp)
pos = (1:ncol(mypvals))[kp] 
choose_pvalues <- mypvals[,c(c(1,5,6,8,9),pos  )]
choose_means <- mymeans[,c(c(1,5,6,8,9),pos)]

logi <- apply(choose_pvalues[,5:ncol(choose_pvalues)]<0.05, 1, sum) 
choose_pvalues <- choose_pvalues[logi>=1,]

logi1 <- choose_pvalues$gene_a != ""
logi2 <- choose_pvalues$gene_b != ""
logi <- logi1 & logi2
choose_pvalues <- choose_pvalues[logi,]

choose_means <- choose_means[choose_means$id_cp_interaction %in% choose_pvalues$id_cp_interaction,]

library(tidyverse)
meansdf <- choose_means %>% reshape2::melt()
meansdf <- data.frame(interacting_pair = paste0(meansdf$gene_a,"_",meansdf$gene_b),
                      CC = meansdf$variable,
                      means = meansdf$value)
pvalsdf <- choose_pvalues %>% reshape2::melt()
pvalsdf <- data.frame(interacting_pair = paste0(pvalsdf$gene_a,"_",pvalsdf$gene_b),
                      CC = pvalsdf$variable,
                      pvals = pvalsdf$value)

pvalsdf$joinlab<- paste0(pvalsdf$interacting_pair,"_",pvalsdf$CC)
meansdf$joinlab<- paste0(meansdf$interacting_pair,"_",meansdf$CC)
pldf <- merge(pvalsdf,meansdf,by = "joinlab")

# dotplot
summary((filter(pldf,means >0))$means)
head(pldf)
library(viridis)
pcc =  pldf%>% filter(means >0) %>% 
  ggplot(aes(CC.x,interacting_pair.x) )+ 
  geom_point(aes(color=means,size=-log10(pvals+0.0001)) ) +
  scale_size_continuous(range = c(0,3))+
  scale_color_viridis(option='turbo')+ 
  theme_bw()+ 
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.8))
pdf("yuanwang/pan-cancer/ICD/figure/SingleCell/cellphoneDB-bubble.pdf",width = 8,height = 12)
print(pcc)
dev.off()

library(psych)
library(qgraph)
library(igraph)
library(cols4all)

mynet <- read.delim('out/count_network.txt', check.names = FALSE)
net<- graph_from_data_frame(mynet)
plot(net)
allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB",
            "#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3",
            "#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
            "#40E0D0","#5F9EA0","#FF1493",
            "#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
            "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347",
            "#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
karate_groups <- cluster_optimal(net)
coords <- layout_in_circle(net, order =
                             order(membership(karate_groups))) 

E(net)$width  <- E(net)$count/10 
plot(net, edge.arrow.size=.1, 
     edge.curved=0,
     vertex.color=allcolour,
     vertex.frame.color="#555555",
     vertex.label.color="black",
     layout = coords,
     vertex.label.cex=.7) 
net2 <- net

for (i in 1: length(unique(mynet$SOURCE)) ){
  E(net)[map(unique(mynet$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
  })%>% unlist()]$color <- allcolour[i]
} 

pdf("yuanwang/pan-cancer/ICD/figure/SingleCell/cellphoneDB-interaction.pdf",width = 5,height = 5)
plot(net, edge.arrow.size=.1, 
     edge.curved=0.2, #
     vertex.color=allcolour,
     vertex.frame.color="#555555",
     vertex.label.color="black",
     layout = coords,
     vertex.label.cex=.7) 
dev.off()

length(unique(mynet$SOURCE))
pdf("yuanwang/pan-cancer/ICD/figure/SingleCell/cellphoneDB-single.pdf",width = 12,height = 5)
par(mfrow=c(2,4), mar=c(.3,.3,.3,.3))
for (i in 1: length(unique(mynet$SOURCE)) ){
  net1<-net2
  
  E(net1)$count <- ""
  E(net1)[map(unique(mynet$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
  })%>% unlist()]$count  <- E(net2)[map(unique(mynet$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
  })%>% unlist()]$count 
  
  E(net1)[map(unique(mynet$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
  })%>% unlist()]$color <- allcolour[i]
  
  plot(net1, edge.arrow.size=.1, 
       edge.curved=0.4,
       edge.label = E(net1)$count,
       vertex.color=allcolour,
       vertex.frame.color="#555555",
       vertex.label.color="black",
       layout = coords,
       vertex.label.cex=1
  ) 
}
dev.off()
