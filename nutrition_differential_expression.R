## nutrition specific gene expression data 
exps <- c("-P / +Fe (3h)", "-P / +Fe (6h)", "-P / +Fe (9h)",  "+P / -Fe (3h)", "+P / -Fe (6h)", "+P / -Fe (9h)", "-P / -Fe (3h)", "-P / -Fe (6h)", "-P / -Fe (9h)")
ill_defined <- c("UNIDENTIFIED", "AMBIGUOUS", "RESCUE")

df.geneontology <- read.csv("data/Table annotation_INITIAL.csv", stringsAsFactors = FALSE) #read.table("Table annotation_INITIAL.txt", header = TRUE, sep = "\t", quote = "")

v.map <- df.geneontology$AGI
names(v.map) <- df.geneontology$Transcript.Cluster.ID

v.anova <- read.csv("data/anova_tests.csv", stringsAsFactors = FALSE, header = FALSE)
df.anova <- read.csv("data/anovaResults.csv", stringsAsFactors = FALSE)
agi <- v.map[as.character(df.anova$Cluster)]
df.anova <- df.anova[,-1]

v.anova <- v.anova[,1]
names(df.anova) <- v.anova

m.anova <- as.matrix(df.anova)
rownames(m.anova) <- agi

# subset all datasets with p-value smaller 0.05
m.anova[m.anova > 0.05] <- 10
m.anova[m.anova <= 0.05] <- 1
m.anova[m.anova == 10] <- 0

m.anova.set <- m.anova[,c(49,50,51,87,88,89,109,110,111),drop = FALSE]
m.anova.set <- m.anova.set[,c(1,4,7,2,5,8,3,6,9)]
m.anova.set <- m.anova.set[which(rowSums(m.anova.set) > 0),]  ## subset of diff. expr. genes 

gns.DE <- rownames(m.anova.set)
gns.DE <- gns.DE[!gns.DE %in% ill_defined]
m.anova.set <- m.anova.set[gns.DE, ]


saveRDS(gns.DE, "gns.DE.rds")

### raw gene expression data ### 
df.geneExp <- read.csv("data/geneExpData.csv", header = TRUE, stringsAsFactors = FALSE)
df.geneExp["AGI"] <- v.map[as.character(df.geneExp$Cluster)]
df.geneExp <- df.geneExp[,-1]

v.gnSets <- names(which(table(df.geneExp$AGI) == 1))
df.geneExp <- subset(df.geneExp, df.geneExp$AGI %in% v.gnSets)
v.gnSets <- df.geneExp$AGI

rownames(df.geneExp) <- v.gnSets

df.geneExp <- (df.geneExp[1:13] + df.geneExp[14:26] + df.geneExp[27:39])/3

# compute the standard deviations for the fold changes

df.geneExp <- df.geneExp[,c(1,2,6,10,3,7,11,4,8,12,5,9,13)]

df.geneExp.set <- df.geneExp[,-1]
df.geneExp.set <- df.geneExp.set[,c(1,4,7,10,2,5,8,11,3,6,9,12)]
df.geneExp.set[,2:4] <- df.geneExp.set[,2:4] - df.geneExp.set[,1]
df.geneExp.set[,6:8] <- df.geneExp.set[,6:8] - df.geneExp.set[,5]
df.geneExp.set[,10:12] <- df.geneExp.set[,10:12] - df.geneExp.set[,9] # log fold change

df.geneExp.set <- df.geneExp.set[,c(2:4,6:8,10:12)]

###

v.0PF <- apply(df.geneExp.set[gns.DE,c(1,4,7)],1, mean)
v.P0F <- apply(df.geneExp.set[gns.DE,c(2,5,8)],1, mean)
v.0P0F <- apply(df.geneExp.set[gns.DE,c(3,6,9)],1, mean)



message("plot condition similarity")

library(ggplot2)
library(ggdendro)

hc <- hclust(dist(t(data.frame("0PF" = v.0PF, "P0F" = v.P0F, "0P0F" = v.0P0F))), "ave")
ggdendrogram(hc, rotate = FALSE, size = 4) + labs(title="Dendrogram in ggplot2")

dendr <- dendro_data(hc, type="rectangle")
#your own labels are supplied in geom_text() and label=label
ggplot() +
  geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_text(data=label(dendr), aes(x=x, y=y, label=c("  - Ph / + Fe", "  + Ph / - Fe", "  - Ph / - Fe"), hjust=0), size=4) +
  coord_flip() + scale_y_reverse(expand=c(0.2, 0)) +
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank())





#### gene set enrichment analysis


library("ggplot2")
library("ggdendro")
library("reshape2")

df.de <- df.geneExp.set[gns.DE,c(1,4,7,2,5,8,3,6,9)]
df.de <- na.omit(df.de)
m.de <- as.matrix(df.de, drop = FALSE)
colnames(m.de) <- exps

F_m2 <- m.de

# m <- heatmap((m.de), Colv=F, scale='none')

# extract 
# https://liuyanguu.github.io/post/2018/07/16/how-to-draw-heatmap-with-colorful-dendrogram/

library(curl)       # read file from google drive
library(gplots)     # heatmap.2
library(dendextend) # make and color dendrogram
library(colorspace) # diverge_hcl / rainbow_hcl / heat_hcl color palettes


dev.off()

c_group <- 8 # number of clusters
hc <- hclust(dist(F_m2))
ct <- cutree(hc, k = c_group)

dend1 <- as.dendrogram(hc)
dend1 <- color_branches(dend1, k = c_group, col = rainbow_hcl) # add color to the lines
dend1 <- color_labels(dend1, k = c_group, col = rainbow_hcl)   # add color to the labels

# reorder the dendrogram, must incl. `agglo.FUN = mean`
rMeans <- rowMeans(F_m2, na.rm = T)
dend1 <- reorder(dend1, rowMeans(F_m2, na.rm = T), agglo.FUN = mean)

# get the color of the leaves (labels) for `heatmap.2`
col_labels <- get_leaves_branches_col(dend1)
col_labels <- col_labels[order(order.dendrogram(dend1))]

# if plot the dendrogram alone:
# the size of the labels:
dend1 <- set(dend1, "labels_cex", 0.5)
par(mar = c(1,1,1,14))
plot_horiz.dendrogram(dend1, side = F) # use side = T to horiz mirror if needed

###

## plot the heatmap with the dendrogram above ##
par(cex.main=0.5)                   # adjust font size of titles
heatmap.2(F_m2, main = '',
          # reorderfun=function(d, w) reorder(d, w, agglo.FUN = mean),
          # order by branch mean so the deepest color is at the top
          dendrogram = "row",        # no dendrogram for columns
          Rowv = dend1,              # * use self-made dendrogram
          Colv = "NA",               # make sure the columns follow data's order
          col = diverge_hcl,         # color pattern of the heatmap
          
          trace="none",              # hide trace
          density.info="none",       # hide histogram
          
          margins = c(5,18),         # margin on top(bottom) and left(right) side.
          cexRow=0.4, cexCol = 0.8,      # size of row / column labels
          xlab = "",
          srtCol=90, adjCol = c(0.5,1), # adjust the direction of row label to be horizontal
          # margin for the color key
          # ("bottom.margin", "left.margin", "top.margin", "left.margin" )
          key.par=list(mar=c(10,1,2,1)),
          RowSideColors = col_labels, # to add nice colored strips        
          colRow = col_labels         # add color to label
)

##### gsea of gene clusters ### 

for(j in 1:c_group){
  
   v.gns <- names(which(ct == j))
   go_gsea(v.gns, th = 0.1, ontology = "BP")
   go_gsea(v.gns, th = 0.1, ontology = "MF")
  
}


## all genes gsea 
go_gsea(gns.DE, th = 0.1, ontology = "BP")
go_gsea(gns.DE, th = 0.1, ontology = "MF")


### gene regulatory network  ###

### nutrition differential expression based network ### 

tgs <- intersect(rownames(m.de), gns.DE)
tfs <- intersect(names(v.regulators), gns.DE)

df.grn <- c()

for(j in 1:length(tfs)){
  
  tf <- tfs[j]
  df.grn.j <- data.frame(TF=rep(tf,length(tgs)), TG=tgs)
  df.grn.j["0PF"] <- 0
  df.grn.j["P0F"] <- 0
  df.grn.j["0P0F"] <- 0
  # df.grn.j["mode"] <- 0
  
  i.min <- min(which(m.anova.set[tf, 1:3] == 1))
  
  if(i.min != Inf){
    
    tgs.de <- intersect(tgs, rownames(m.anova.set)[which(rowSums(m.anova.set[, i.min:3]) >= 1)])
    idx = which(df.grn.j$TG %in% tgs.de)
    df.grn.j$`0PF`[idx] <- 1

    # for(i in 1:length(tgs.de)){
    #   idx <- min(which(m.anova.set[tgs.de[i], 1:3] == 1))
    #   df.grn.j$mode[i] <- sign(m.de[tf, i.min]) * sign(m.de[tgs.de[i], idx])
    # }
  
  }
  
  
  i.min <- min(which(m.anova.set[tf, 4:6] == 1))
  
  if(i.min != Inf){
    
    tgs.de <- intersect(tgs, rownames(m.anova.set)[which(rowSums(m.anova.set[, (3 + i.min):6]) >= 1)])
    idx = which(df.grn.j$TG %in% tgs.de)
    df.grn.j$P0F[idx] <- 1

    # for(i in 1:length(tgs.de)){
    #   idx <- min(which(m.anova.set[tgs.de[i], (3 + i.min):6] == 1))
    #   df.grn.j$mode[i] <- sign(m.de[tf, i.min]) * sign(m.de[tgs.de[i], idx])
    # }

  }
  
  
  i.min <- min(which(m.anova.set[tf, 7:9] == 1))
  
  if(i.min != Inf){
    
    tgs.de <- intersect(tgs, rownames(m.anova.set)[which(rowSums(m.anova.set[, (6 + i.min):9]) >= 1)])
    idx = which(df.grn.j$TG %in% tgs.de)
    df.grn.j$`0P0F`[idx] <- 1
    
    # for(i in 1:length(tgs.de)){
    #   idx <- min(which(m.anova.set[tgs.de[i], (6 + i.min):9] == 1))
    #   df.grn.j$mode[i] <- sign(m.de[tf, i.min]) * sign(m.de[tgs.de[i], idx])
    # }

  }
  
  if(length(df.grn.j) > 0){
    # df.grn.j$coDiffExp <- pmax(pmax(df.grn.j$`0PF`, df.grn.j$P0F), df.grn.j$`0P0F`)
    df.grn <- rbind(df.grn, df.grn.j)
  }
}

# write.csv(df.grn, "manuscript/codifferential_expression/df.grn.codifferential.csv", row.names = FALSE)


### create the entire network -> complete 0.5 
df.code <- df.grn
df.code <- df.code[which(rowSums(df.code[,3:5]) > 0),]


# TF        TG      val spec
# 3055 AT5G61430 AT5G54230 1.850967    0
# 2718 AT1G77200 AT5G06750 1.845253    0


### co-diff heatmaps #### 

tfs <-  unique(df.grn$TF)

l.p.tgs <- vector(mode = "list", length = length(tfs))

for(r in 1:length(tfs)){
  
  df <- subset(df.grn, df.grn$TF == tfs[r])
  df <- subset(df, df$TF != df$TG)
  tgs <- unique(df$TG)

  m.exp <- m.de[tgs,]
  
  m.sig <- m.anova.set[tgs,]
  m.sig[m.sig == 1] <- "*"
  m.sig[m.sig == 0] <- ""
  
  v.sig.tf <- m.anova.set[tfs[r],]
  v.sig.tf[v.sig.tf == 1] <- "*"
  v.sig.tf[v.sig.tf == 0] <- ""
  
  if(dim(m.exp)[1] > 1){
    
    dd.col <- as.dendrogram(hclust(dist((m.exp))))
    row.col <- order.dendrogram(dd.col)
  
    m.exp <- m.exp[row.col,]
    m.sig <- m.sig[row.col,]
    
    tgs <- rownames(m.exp)
    
  }
  
  v.exp.tf <- m.de[tfs[r],]
  names(v.exp.tf) <- exps
  
  m.exp <- rbind(m.exp, v.exp.tf)
  m.sig <- rbind(m.sig, v.sig.tf)
  
  m.exp <- as.matrix(m.exp)
  m.sig <- as.matrix(m.sig)
  
  rownames(m.exp) <- rownames(m.sig) <- c(tgs, tfs[r])
  
  max.line <- length(tgs) + 0.5 + 1
  p.tgs <- plot_ratios(cormat=m.exp, m.value=m.sig, v.title = paste("Targets of ", tfs[r], sep =""), v.name = "logFC", v.max = max(m.exp), v.min = min(m.exp), max.line = max.line)
  
  l.p.tgs[[r]] <- p.tgs # plot as 10 x 30
  
}
# 
# p1 <- plot_ratios(cormat=m.test[,1:3], v.title = "", v.name = "", v.max = max(m.test), v.min = min(m.test))
# p2 <- plot_ratios(cormat=m.test[,4:6], v.title = "", v.name = "", v.max = max(m.test), v.min = min(m.test))
# p3 <- plot_ratios(cormat=m.test[,7:9], v.title = "", v.name = "", v.max = max(m.test), v.min = min(m.test))


library(gridExtra)
library(grid)
grid.arrange(l.p.tgs[[1]], l.p.tgs[[2]], l.p.tgs[[3]], l.p.tgs[[4]], l.p.tgs[[5]], l.p.tgs[[6]], l.p.tgs[[7]], l.p.tgs[[8]], l.p.tgs[[9]], l.p.tgs[[10]], ncol = 5,  top = "log foldChange differential expression patterns of putative regulators and targets across all three experimental time series") # plot the three expression b



### tigre 


# install.packages("BiocManager")
BiocManager::install("tigre") 


library(tigre)


times <- c(3, 6, 9)

df.exp <- processRawData(m.de[,1:3], times, experiments=rep(1, each=3), is.logged = TRUE)#, do.normalisation = ifelse(is.logged, TRUE, FALSE))


# Learn the model using only one of the 3 repeats in the data
model <- GPLearn(df.exp,
                 TF=tfs[1], targets = rsk1,
                 useGpdisim=TRUE, quiet=TRUE)

#m.expression.diffExp.set[targets[1], ]

# Display the model parameters
show(model)

GPPlot(model)

m.loglik[i,j] <- - model@model$llscore






