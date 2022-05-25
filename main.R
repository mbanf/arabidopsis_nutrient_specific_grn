## main script - build and visualize the network

folder_data = "data"
folder_tmp = "tmp"
folder_output = "results"

ifelse(!dir.exists(file.path(folder_output)), dir.create(file.path(folder_output)), FALSE)
ifelse(!dir.exists(file.path(folder_tmp)), dir.create(file.path(folder_tmp)), FALSE)


source("utils.R")
install_and_load_libraries()


source("nutrition_differential_expression.R")
source("transcriptional_regulators.R")
source("root_differential_expression.R")

specificity = c("-P+Fe", "+P-Fe", "-P-Fe", 
                "-P+Fe, +P-Fe", "-P+Fe, -P-Fe",
                "+P-Fe, -P-Fe", "-P+Fe, +P-Fe, -P-Fe")

code_specificity = seq(1, length(specificity))
names(code_specificity) <- specificity

cols_specificity = c( "cyan", "tomato", "steelblue", "black", "pink",  "purple", "sienna") 
names(cols_specificity) <- seq(1, length(specificity))

m.grn <- matrix(0, nrow=length(v.regulators), ncol = length(gns.DE))
rownames(m.grn) <- names(v.regulators)
colnames(m.grn) <- gns.DE

m.spec <- matrix(0, nrow=length(v.regulators), ncol = length(gns.DE))
rownames(m.spec) <- names(v.regulators)
colnames(m.spec) <- gns.DE



### creating grn

tgs <- intersect(rownames(m.de), gns.DE)
tfs <- intersect(names(v.regulators), gns.DE)


### gene specificities
gn_spec <- rep(0, length(tgs))
names(gn_spec) <- tgs

for(i in 1:length(tgs)){
  
  gn = tgs[i]
  npf <- min(1,sum(m.anova.set[gn, 1:3]))
  pnf <- min(1,sum(m.anova.set[gn, 4:6]))
  npnf <- min(1,sum(m.anova.set[gn, 7:9])) 
  
  
  if(npf == 1 & pnf == 0 & npnf == 0){
    gn_spec[i] <- 1}
  
  if(npf == 0 & pnf == 1 & npnf == 0){
    gn_spec[i] <- 2}
  
  if(npf == 0 & pnf == 0 & npnf == 1){
    gn_spec[i] <- 3}
  
  if(npf == 1 & pnf == 1 & npnf == 0){
    gn_spec[i] <- 4}
  
  if(npf == 1 & pnf == 0 & npnf == 1){
    gn_spec[i] <- 5}
  
  if(npf == 0 & pnf == 1 & npnf == 1){
    gn_spec[i] <- 6}
  
  if(npf == 1 & pnf == 1 & npnf == 1){
    gn_spec[i] <- 7}
  
}

table(gn_spec)



df.grn <- c()

for(j in 1:length(tfs)){
  
  tf <- tfs[j]
  df.grn.j <- data.frame(TF=rep(tf,length(tgs)), TG=tgs)
  df.grn.j["0PF"] <- 0
  df.grn.j["P0F"] <- 0
  df.grn.j["0P0F"] <- 0
  
  i.min <- min(which(m.anova.set[tf, 1:3] == 1))
  
  if(i.min != Inf){
    tgs.de <- intersect(tgs, rownames(m.anova.set)[which(rowSums(m.anova.set[, i.min:3]) >= 1)])
    idx = which(df.grn.j$TG %in% tgs.de)
    df.grn.j$`0PF`[idx] <- 1
  }
  
  i.min <- min(which(m.anova.set[tf, 4:6] == 1))
  
  if(i.min != Inf){
    tgs.de <- intersect(tgs, rownames(m.anova.set)[which(rowSums(m.anova.set[, (3 + i.min):6]) >= 1)])
    idx = which(df.grn.j$TG %in% tgs.de)
    df.grn.j$P0F[idx] <- 1
  }
  
  i.min <- min(which(m.anova.set[tf, 7:9] == 1))
  if(i.min != Inf){
    tgs.de <- intersect(tgs, rownames(m.anova.set)[which(rowSums(m.anova.set[, (6 + i.min):9]) >= 1)])
    idx = which(df.grn.j$TG %in% tgs.de)
    df.grn.j$`0P0F`[idx] <- 1
  }
  if(length(df.grn.j) > 0){
    df.grn <- rbind(df.grn, df.grn.j)
  }
}

df.code <- df.grn
df.code <- df.code[which(rowSums(df.code[,3:5]) > 0),]
message("# interactions co-diff expressions: ", nrow(df.code))




for(i in 1:nrow(df.code)){
  
  tf = df.code$TF[i]
  tg = df.code$TG[i]
  if(tf != tg){
    
    m.grn[tf, tg] <- 1
    
    if(df.code$`0PF`[i] == 1 & df.code$P0F[i] == 0 & df.code$`0P0F`[i] == 0){
      m.spec[tf, tg] <- 1}
    
    if(df.code$`0PF`[i] == 0 & df.code$P0F[i] == 1 & df.code$`0P0F`[i] == 0){
      m.spec[tf, tg] <- 2}
    
    if(df.code$`0PF`[i] == 0 & df.code$P0F[i] == 0 & df.code$`0P0F`[i] == 1){
      m.spec[tf, tg] <- 3}
    
    if(df.code$`0PF`[i] == 1 & df.code$P0F[i] == 1 & df.code$`0P0F`[i] == 0){
      m.spec[tf, tg] <- 4}
    
    if(df.code$`0PF`[i] == 1 & df.code$P0F[i] == 0 & df.code$`0P0F`[i] == 1){
      m.spec[tf, tg] <- 5}
    
    if(df.code$`0PF`[i] == 0 & df.code$P0F[i] == 1 & df.code$`0P0F`[i] == 1){
      m.spec[tf, tg] <- 6}
    
    if(df.code$`0PF`[i] == 1 & df.code$P0F[i] == 1 & df.code$`0P0F`[i] == 1){
      m.spec[tf, tg] <- 7}
    
  }
}


# create grn from DNA binding data - set threshold at 0.5
# load DNA binding dataset (created with tfbs project https://github.com/mbanf/TFBS)
df.dna <- readRDS(paste(folder_data, "m.motifNet.rds", sep = "/"))
df.dna <- melt(df.dna)
names(df.dna) <- c("TF", "TG", "val")
df.dna <- subset(df.dna, df.dna$val > 0)
nrow(df.dna)
message("# interactions dna binding: ", nrow(df.dna))

for(i in 1:nrow(df.dna)){
  tf <- as.character(df.dna$TF[i])
  tg <- as.character(df.dna$TG[i])
  if(tf != tg){
    val <- df.dna$val[i]
    if(val > 0.5){  # minimum binding confidence
      m.grn[tf, tg] <- m.grn[tf, tg] + val
    }else{
      m.grn[tf, tg] <- 0 # delete even co-diff links! (if the motif is known)
    }
  }
}

for(i in 1:nrow(df.rf.grn)){
  tf <- as.character(df.rf.grn$TF[i])
  tg <- as.character(df.rf.grn$TG[i])
  if(tf != tg){ # evidence for co-diff without binding motif in tf
    m.grn[tf, tg] <- m.grn[tf, tg] + as.numeric(df.rf.grn$val[i])
  }
}

df.grn <- melt(m.grn)
names(df.grn) <- c("TF", "TG", "val")

df.spec <- melt(m.spec)
names(df.spec) <- c("TF", "TG", "spec")

df.grn["spec"] <- df.spec["spec"]

df.grn <- df.grn[order(-df.grn$val),]

hist(df.grn$val) # figure - distribution integrated grn

df.grn <- subset(df.grn, df.grn$spec != 0)
df.grn <- subset(df.grn, df.grn$val > 1.5) # dna binding 

tfs.final = unique(as.character(df.grn$TF))

table(df.grn$spec)
nrow(df.grn)
message("# interactions final: ", nrow(df.grn))
# idx <- sort(unique(df.grn$spec))
# specificity = specificity[idx]
# cols_specificity = cols_specificity[idx]

df.grn.result <- df.grn[,c("TF", "TG", "val")]
df.grn.result["specifity"] <- specificity[df.grn$spec]

message("ensemble based grn with specificity with ", nrow(df.grn.result), " links, containing ", length(unique(df.grn.result$TF)), " TFs and ", length(unique(df.grn.result$TG)), " targets")

write.csv(df.grn.result, paste(folder_output, "df.grn_w_specificity.csv", sep = "/"), row.names = F)

#####  

v.regs <- c("IAA5", "ERF37", "ERF36", "bHLH93","MYB49", "ERF6","NAC100","bHLH23","bHLH39","PLATZ11","MYB76","DAZ3") 
names(v.regs) <- c("AT1G15580", "AT1G77200", "AT3G16280", "AT5G65640", "AT5G54230", "AT4G17490","AT5G61430","AT4G28790","AT3G56980","AT4G17900","AT5G07700","AT4G35700")

rsk1 = "AT2G26290"


library(igraph)

# triangle vertex shape
mytriangle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  
  symbols(x=coords[,1], y=coords[,2], bg=vertex.color,
          stars=cbind(vertex.size, vertex.size, vertex.size),
          add=TRUE, inches=FALSE)
}
# clips as a circle
add_shape("triangle", clip=shapes("circle")$clip,
          plot=mytriangle)


g=graph.edgelist(as.matrix(df.grn[,1:2]))


# Node size - accorind to number of targets or strength!1
hubs = table(df.grn$TF)
hubs <- hubs / max(hubs)
# specifitiy of regulation - distribution per hub node! 
# highlight RSK1 in size

V(g)$size <-  ifelse(names(V(g)) == rsk1, 4, ifelse(names(V(g)) %in% df.grn$TF, hubs[names(V(g))] * 2 + 1, 2))
E(g)$weight <- df.grn$val / max(df.grn$val)
E(g)$color <- cols_specificity[df.grn$spec]
V(g)$color <- cols_specificity[gn_spec[names(V(g))]]
# V(g)$color <- ifelse(names(V(g)) == rsk1, adjustcolor("blue", alpha.f = .7) , ifelse(names(V(g)) %in% tfs.final, adjustcolor("orange", alpha.f = .6), adjustcolor("gray", alpha.f = .6)))

V(g)$shape <- ifelse(names(V(g)) == rsk1, "triangle" , ifelse(names(V(g)) %in% tfs.final, "circle", "square"))


#vertex.label <- ifelse(names(V(g)) %in% tfs.final, paste(names(V(g)), "(", v.regulators[names(V(g))], ")", sep="") , "")  # TODO: this might not work?
#idx <- which(names(V(g)) == rsk1)
#vertex.label[idx] <- "RSK1" #paste(rsk1, "(RSK1)")

vertex.label <- ifelse(names(V(g)) %in% tfs.final, v.regs[names(V(g))] , "")  # TODO: this might not work?
idx <- which(names(V(g)) == rsk1)
vertex.label[idx] <- "RSK1"



# large supplement version 
set.seed(9999) 
plot(g, edge.arrow.size=.1, 
     edge.curved=seq(-0.5, 0.5, length = ecount(g)),
     vertex.label.dist=0.4,
     vertex.label.cex=0.7,
     vertex.label.color = "black",
     vertex.label=vertex.label,edge.width=E(g)$weight)#, main = paste(v.tissues[s], " / ", names(l.grn_treatment[[s]])[i], " / ", v.domains[d], sep =""))

# plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)

# specificity
legend("bottomleft", legend=specificity , col = cols_specificity , bty = "n", pch="-" , pt.cex = 2, cex = 0.8, horiz = FALSE, inset = c(0.0, 0.2))

# vertex type
legend("bottomleft", legend=c("Transcrition factor", "RSK1", "Other gene") , col = "black" , bty = "n", pch=c(16, 17, 15), pt.cex = 1, cex = 0.8, horiz = FALSE, inset = c(0.0, 0.1))
# 13 * 20


# more compacted version
# 
# Other graph layouts: add_layout_, component_wise, layout_as_bipartite, layout_as_star, 
# layout_as_tree, layout_in_circle, layout_nicely, layout_on_grid, layout_on_sphere, layout_randomly, layout_with_dh, layout_with_fr, 
# layout_with_gem, layout_with_graphopt, layout_with_kk, layout_with_lgl, 
# layout_with_mds, layout_with_sugiyama, merge_coords, norm_coords, normalize
# 


# set.seed(9999) 
# coords <- layout_(g, as_star())
# plot(g, 
#      layout = coords,
#      edge.arrow.size=.1, 
#      edge.curved=seq(-0.5, 0.5, length = ecount(g)),
#      vertex.label.dist=0.4,
#      vertex.label.cex=0.7,
#      vertex.label.color = "black",
#      vertex.label=vertex.label,edge.width=E(g)$weight)#, main = paste(v.tissues[s], " / ", names(l.grn_treatment[[s]])[i], " / ", v.domains[d], sep =""))
# 
# legend("bottomleft", legend=specificity , col = cols_specificity , bty = "n", pch="-" , pt.cex = 2, cex = 0.8, horiz = FALSE, inset = c(-0.13, 0.2))
# 
# # vertex type
# legend("bottomleft", legend=c("Transcrition factor", "RSK1", "Other gene") , col = "black" , bty = "n", pch=c(16, 17, 15), pt.cex = 1, cex = 0.8, horiz = FALSE, inset = c(-0.13, 0.1))
# # 13 * 20


#### subplot -P+Fe

# main players!
ERF036 = "AT3G16280"
ERF037 = "AT1G77200"
MYB49 = "AT5G54230"
RSK1 = "AT2G26290"


df <- subset(df.grn, df.grn$spec %in% c(1,4))
g=graph.edgelist(as.matrix(df[,1:2]))


# Node size - accorind to number of targets or strength!1
hubs = table(df$TF)
hubs <- hubs / max(hubs)
# specifitiy of regulation - distribution per hub node! 
# highlight RSK1 in size

V(g)$size <-  ifelse(names(V(g)) == rsk1, 7, ifelse(names(V(g)) %in% df$TF, hubs[names(V(g))] * 4 + 1, 2))
E(g)$weight <- df$val / max(df$val)
E(g)$lty <- 3


idx <- which(df$TG == rsk1)
E(g)$lty[idx] <- 1

# up and down regultion --- boring :) 
# E(g)$color <- ifelse(df.grn$mode == 1, "red", ifelse(df.grn$mode == -1, "green", "gray"))
# E(g)$color <- ifelse(df.grn$mode == 1, "red", ifelse(df.grn$mode == -1, "green", "gray"))

E(g)$color <- "gray" #  adjustcolor("gray", alpha.f = .9)
E(g)$color[idx] <- adjustcolor("brown", alpha.f = .8)


# V(g)$color <- cols_specificity[gn_spec[names(V(g))]]
V(g)$color <- ifelse(names(V(g)) == rsk1, adjustcolor("blue", alpha.f = .7) , ifelse(names(V(g)) %in% tfs.final, adjustcolor("orange", alpha.f = .6), adjustcolor("darkgreen", alpha.f = .6)))
V(g)$shape <- ifelse(names(V(g)) == rsk1, "triangle" , ifelse(names(V(g)) %in% tfs.final, "circle", "square"))

# 
# V(g)$color 

vertex.label <- ifelse(names(V(g)) %in% tfs.final, v.regs[names(V(g))] , "")  # TODO: this might not work?
idx <- which(names(V(g)) == rsk1)
vertex.label[idx] <- "RSK1"

V(g)$label.cex <- 0.7
V(g)[ERF036]$label.cex <- 1.4
V(g)[ERF037]$label.cex <- 1.4
V(g)[MYB49]$label.cex<- 1.4
V(g)[RSK1]$label.cex <- 1.4

set.seed(7)
plot(g, edge.arrow.size=.1,  # 10 x 10 
     edge.curved=seq(-0.5, 0.5, length = ecount(g)),
     vertex.label.dist=0.95,
     vertex.label.color = "black",
     vertex.label=vertex.label,edge.width=E(g)$weight)#, main = paste(v.tissues[s], " / ", names(l.grn_treatment[[s]])[i], " / ", v.domains[d], sep =""))





##### specific expression and p-values .. 

ERF036 = "AT3G16280"
ERF037 = "AT1G77200"
MYB49 = "AT5G54230"
RSK1 = "AT2G26290"

gn = MYB49 #ERF036 # TODO: repeat per gene, make individual or combined box plots!!

### raw gene expression data ### 
df.geneExp <- read.csv(paste(folder_data, "geneExpData.csv", sep = "/"), header = TRUE, stringsAsFactors = FALSE)
df.geneExp["AGI"] <- v.map[as.character(df.geneExp$Cluster)]
df.geneExp <- df.geneExp[,-1]

v.gnSets <- names(which(table(df.geneExp$AGI) == 1))
df.geneExp <- subset(df.geneExp, df.geneExp$AGI %in% v.gnSets)
v.gnSets <- df.geneExp$AGI

rownames(df.geneExp) <- v.gnSets
df.geneExp <- df.geneExp[c(ERF036, ERF037, MYB49, RSK1),]

df3 = df.geneExp[,c(3, 16, 29)] - df.geneExp[,c(2,15,28)]
df6 = df.geneExp[,c(7, 20, 33)] - df.geneExp[,c(6,19,32)]
df9 = df.geneExp[,c(11, 24, 37)] - df.geneExp[,c(10,23,36)]

y = as.numeric(c(df3[gn,], df6[gn,], df9[gn,]))
x = c(rep(1,3), rep(2,3), rep(3,3))

labels = c("3h", "6h", "9h")

boxplot(y~x, ylab = "logFC", xlab = "", xaxt="n", las=2,  frame = FALSE)

text(x=c(1,2,3),  par("usr")[3], labels = labels, srt = 0, pos = 1, xpd = TRUE)

print(m.anova.set[c(ERF036, ERF037, MYB49, RSK1),1:3])

text(x=1,y=1.5,"*",pos=3,cex=1.5) # MYB49
# text(x=2,y=-0.3,"*",pos=3,cex=1.5) # ARSK1 - 

# TODO: add significance to plots, make loop 



### specific gsea 
source("gsea.R")
df <- subset(df.grn, df.grn$spec %in% c(1,4)) # subset specificitys - why 1 and 4?
gn <- as.character(unique(df$TG))

# enrichment
go_gsea(gn, gn.pop = NULL, th = 0.05, mode = "enrichment", bg.mode = "genome", n.genome = 27655) #0PFE vs genome
go_gsea(unique(df.grn$TG), gn.pop = NULL, th = 0.05, mode = "enrichment", bg.mode = "genome", n.genome = 27655) #GRN vs GENOME
go_gsea(gn, gn.pop = unique(df.grn$TG), th = 0.05, mode = "enrichment") #0PFE vs GRN

# depletion
go_gsea(gn, gn.pop = NULL, th = 0.05, mode = "depletion", bg.mode = "genome", n.genome = 27655) #0PFE vs genome
# go_gsea(unique(df.grn$TG), gn.pop = NULL, th = 0.05, mode = "depletion", bg.mode = "genome", n.genome = 27655) #GRN vs GENOME
go_gsea(gn, gn.pop = unique(df.grn$TG), th = 0.05, mode = "depletion") #0PFE vs GRN




specs_grn <- sort(unique(df.grn$spec))
# 
# for(s in 1:length(specs_grn)){
#   
#   specs_grn[s]
#   
#   
# }
# 
# df.grn

# 1 - 1:3
# 2 - 4:6
# 3 - 7:9
# 4 - c(1,2,3,4,5,6)

 
# reduced representation
specs_grn <- c(1,2,3,4)
l.idx <- vector(mode = "list", length = length(specs_grn))
l.idx[[1]] <- c(1,2,3)
l.idx[[2]] <- c(4,5,6)
l.idx[[3]] <- c(7,8,9)
l.idx[[4]] <- c(1,2,3,4,5,6)


# 
# idx = 1:3; s = 1
# 
# idx = 4:6; s = 2
# 
# idx = 7:9; s = 3
# 
# idx = c(1,2,3,4,5,6); s = 4
# 

# idx = 1:9; s = 1 - all representations


tfs <-  unique(df.grn$TF)

l.p.tgs <- vector(mode = "list", length = length(specs_grn))

for(s in 1:length(specs_grn)){
  
  df.grn.s = subset(df.grn, df.grn$spec == specs_grn[s])
  
  tfs <- as.character(unique(df.grn.s$TF))
  tgs <- as.character(unique(df.grn.s$TG))
  # tgs = intersect(tgs, unique(names(which(specs == s))))
  idx <- l.idx[[s]]
  m.exp <- m.de[tgs,idx]
  m.sig <- m.anova.set[tgs,idx]

  ### if true
  m.exp.sig = m.exp * m.sig
  
  ### 
  m.sig[m.sig == 1] <- "*"
  m.sig[m.sig == 0] <- ""
  
  m.sig.tf <- m.anova.set[tfs,idx]
  m.sig.tf[m.sig.tf == 1] <- "*"
  m.sig.tf[m.sig.tf == 0] <- ""
  
  if(dim(m.exp)[1] > 1){
    
    dd.col <- as.dendrogram(hclust(dist((m.exp))))
    row.col <- order.dendrogram(dd.col)
    
    m.exp <- m.exp[row.col,]
    m.sig <- m.sig[row.col,]
    
    tgs <- rownames(m.exp)
    
  }
  
  # plot targets (wo TFs) clustered with dendrogram
  if(FALSE){ # special plot 
    
    dd.col <- as.dendrogram(hclust(dist((m.sig))))
    row.col <- order.dendrogram(dd.col)
    
    m.exp <- m.exp[row.col,]
    m.sig <- m.sig[row.col,]
    
    tgs <- rownames(m.exp)
    
    plot_clean_targets(m.exp, dd.col)
  }
  
  
  m.exp.tf <- m.de[tfs,idx]
  # names(m.exp.tf) <- exps
  
  m.exp <- rbind(m.exp, m.exp.tf)
  m.sig <- rbind(m.sig, m.sig.tf)
  
  m.exp <- as.matrix(m.exp)
  m.sig <- as.matrix(m.sig)
  
  rownames(m.exp) <- rownames(m.sig) <- c(tgs, tfs)
  
  regs = paste(tfs, " (", v.regulators[tfs], ")", sep="")
  
  title = paste(paste("Regulatory module with with", length(tgs), "genes specifically regulated under", names(code_specificity)[s]),
                paste("involved regulators:", paste(regs, collapse=", ")), sep="\n")
  
  p.tgs <- plot_regulons(m.exp, m.sig, 
                          tfs=tfs, tgs=tgs,
                          title = title, v.name = "logFC",  v.max = max(m.exp), v.min = min(m.exp))

  l.p.tgs[[s]] <- p.tgs # plot as 10 x 30
  
}


# 
# p1 <- plot_ratios(cormat=m.test[,1:3], v.title = "", v.name = "", v.max = max(m.test), v.min = min(m.test))
# p2 <- plot_ratios(cormat=m.test[,4:6], v.title = "", v.name = "", v.max = max(m.test), v.min = min(m.test))
# p3 <- plot_ratios(cormat=m.test[,7:9], v.title = "", v.name = "", v.max = max(m.test), v.min = min(m.test))


library(gridExtra)
library(grid)
grid.arrange(l.p.tgs[[1]], l.p.tgs[[2]], l.p.tgs[[3]], l.p.tgs[[4]], ncol = 1,  top = "log foldChange differential expression patterns of putative regulators and targets across all three experimental time series") # plot the three expression b





# ##### gsea of gene clusters ### 
# 
# source("gsea.R")
# 
# for(j in 1:c_group){
#   v.gns <- names(which(ct == j))
#   go_gsea(v.gns, th = 0.1, ontology = "BP")
#   go_gsea(v.gns, th = 0.1, ontology = "MF")
# }
# 
# 
# ## all genes gsea 
# go_gsea(gns.DE, th = 0.1, ontology = "BP")
# go_gsea(gns.DE, th = 0.1, ontology = "MF")
# 
# gns <- names(V(g))
# specs <- gn_spec[gns]
# 
# c_group <- sort(unique(specs))
# 
# for(s in 1:length(code_specificity)){
#   
#   gn = names(which(specs == s))
#   go_gsea(gn, th = 0.05, title = paste("Biological process of", length(gn), "genes specific to", names(code_specificity)[s]), ontology = "BP")
#   go_gsea(gn, th = 0.05, title = paste("Molecular function of", length(gn), "genes specific to", names(code_specificity)[s]),  ontology = "MF")
#   
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ##################
# 
# 
# # devtools::install_github("USCCANA/netplot")
# library(igraph)
# library(netplot)
# 
# 
# vertex.label <- ifelse(names(V(g)) %in% names(v.regulators), paste(names(V(g)), "(", v.regulators[names(V(g))], ")") , "")  # TODO: this might not work?
# idx <- which(names(V(g)) == rsk1)
# vertex.label[idx] <- paste(rsk1, "(RSK1)")
# 
# size = ifelse(names(V(g)) %in% df.grn$TF, 2.5, 2)
# cols = ifelse(names(V(g)) %in% tfs.final, adjustcolor("brown", alpha.f = .6), adjustcolor("gray", alpha.f = .6))
# l <- layout_with_fr(g)
# nplot(g, layout = l, vertex.label = vertex.label,
#       edge.curvature = 0,
#       edge.arrow.size=.01,
#       edge.color = ~ego(mix=0) + alter(mix=1),
#       edge.width=E(g)$weight,
#       vertex.color = cols, vertex.size = size, vertex.label.cex=0.5)
# 
# 
# ans <- nplot(
#   UKfaculty,
#   layout                = l,
#   vertex.color          = viridis::plasma(5)[V(UKfaculty)$Group + 1],
#   vertex.label          = nam,
#   vertex.size.range     = c(.01, .04, 4),
#   vertex.label.col      =  "black",
#   vertex.label.fontface = "bold",
#   bg.col                = "transparent",
#   vertex.label.show     = .5,
#   vertex.label.range    = c(10, 25),
#   edge.width.range      = c(1, 4, 5)
# )
# 
# 
# 
# 
# nplot(x_network, layout = l)
# 
# # Putting two plots in the same page (one using igraph and the other network)
# gridExtra::grid.arrange(
#   nplot(g, layout = l),
#   nplot(g, layout = l), ncol=2, nrow=1
# )
# 
# 
# set.seed(1)
# data("UKfaculty", package = "igraphdata")
# l <- layout_with_fr(UKfaculty)
# 
# plot(UKfaculty, layout = l) # ala igraph
# 
# # Random names
# set.seed(1)
# nam <- sample(babynames::babynames$name, vcount(UKfaculty))
# 
# ans <- nplot(
#   UKfaculty,
#   layout                = l,
#   vertex.color          = viridis::plasma(5)[V(UKfaculty)$Group + 1],
#   vertex.label          = nam,
#   vertex.size.range     = c(.01, .04, 4),
#   vertex.label.col      =  "black",
#   vertex.label.fontface = "bold",
#   bg.col                = "transparent",
#   vertex.label.show     = .5,
#   vertex.label.range    = c(10, 25),
#   edge.width.range      = c(1, 4, 5)
# )
# 
# 
# # Plot it!
# ans
# 
# 
# 
# # w / wo dna binding 
# test <- subset(df.grn, df.grn$TF == "AT3G16280")
# nrow(test[test$val > 6,]) / nrow(test) # binding ratio
# 
# # question: is there a similarily regulated module => similar differentially expressed and dna binding (double significance?)
# # putative real causal "effect" of the regulator => activation, repression => time delay
# 
# # 5 of the 21 genes with dna binding have a similar differential expression pattern! (the also cluster with rsk1 as one of the genes! for this TF)
# # TODO: "correlation between the target genes to identify differential expression similarity rank 
# # TODO: BP analysis of the regulated moduled?
# # TODO: BP analysis of the 21 genes! 
# 
# tgs.db <- as.character(test[test$val > 6,]$TG)
# m.peak = m.anova.set
# colnames(m.peak) <- exps
# colSums(m.peak[tgs.db,]) / length(tgs.db)
# 
# tgs.ndb <- as.character(test[test$val < 6,]$TG)
# colSums(m.peak[tgs.ndb,]) / length(tgs.ndb)
# 
# ### cluster targets together with rsk1
# rsk1 = "AT2G26290"
# 
# library(igraph)
# g=graph.edgelist(as.matrix(df.grn[,1:2]))
# V(g)$size <-  ifelse(names(V(g)) %in% df.grn$TF, 5, 3)
# 
# E(g)$weight <- df.grn$val / max(df.grn$val) * 3
# 
# E(g)$color <- ifelse(df.grn$mode == 1, "red", ifelse(df.grn$mode == -1, "green", "gray"))
# 
# # idx <- which(names(V(g)) == rsk1)
# vertex.label <- ifelse(names(V(g)) %in% names(v.regulators), names(V(g)), "")  # TODO: this might not work?
# # vertex.label[idx] <- paste(rsk1, "(RSK1)")
# 
# 
# ontology = "MF"
# df.GO.annot <- readRDS("data/Athaliana_167.annot.2017.rds")
# df.GO.annot <- subset(df.GO.annot, df.GO.annot$Ontology == ontology)
# n.genes <- length(unique(df.GO.annot$acc.))
# 
# df.GO.annot = go_function(gns.DE, ontology = ontology)
# 
# gns <- unique(as.character(df.GO.annot$acc.))
# na = "no annotation"
# terms <- c(unique(as.character(df.GO.annot$Term)), na)
# term_map = seq(1,length(terms))
# names(term_map) = terms
# 
# values <- list()
# nodes <- names(V(g))
# for(i in 1:length(nodes)){
#   
#   vals = rep(0, length(terms))
#   df = subset(df.GO.annot, df.GO.annot$acc. == nodes[i])
#   if(nrow(df) > 0){
#     idx = as.integer(term_map[unique(as.character(df$Term))])
#     vals[idx] = 1
#   }else{
#     vals[31] = 1
#   }
#   values <- append(values, list(vals))
#   
# }
# 
# # install.packages('colorRamps')
# library(colorRamps)
# cols = list(primary.colors(length(terms),  steps = 3, no.white = TRUE))
# term_color_map = cols[[1]]
# names(term_color_map) = terms
# term_color_map["no annotation"] = "#808080"
# cols[[1]][31] = "#808080"
# 
# # vals <- lapply(1:10, function(x) sample(1:10,3))
# set.seed(1234) 
# plot(g, 
#      vertex.shape="pie", 
#      edge.arrow.size=.2,
#      vertex.pie=values,
#      vertex.pie.color=cols,
#      vertex.label.dist=1,
#      vertex.label.cex=0.6,
#      vertex.label.color = "black",
#      # vertex.size=seq(10,30,length=10), 
#      vertex.label=vertex.label)
# 
# 
# 
# plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
# legend("bottomleft", legend=names(term_color_map)  , col = term_color_map , bty = "n", pch=20 , pt.cex = 2, cex = 1, horiz = FALSE) #, inset = c(0.1, 0.1))
# 
# 
# plot(g, layout = layout, edge.arrow.size=.2, vertex.label=vertex.label,edge.width=E(g)$weight, vertex.label.dist=0.1, vertex.label.color = "black")#, main = paste(v.tissues[s], " / ", names(l.grn_treatment[[s]])[i], " / ", v.domains[d], sep =""))
# 
# 
# 
# 
# cols = unique(col_labels)
# for(j in 1:c_group){
#   v.gns <- names(which(ct == j))
#   idx <- which(names(V(g)) %in% v.gns)
#   V(g)$color[idx] = rep(cols[j], length(idx))
# }
#   
# 
# vertex.label <- ifelse(names(V(g)) %in% names(v.regulators), names(V(g)), "")
# 
# # set RSK1 as node 
# 
# plot(g, edge.arrow.size=.2, vertex.color=V(g)$color, vertex.label=vertex.label, vertex.label.dist=0.1, vertex.label.color = "black", vertex.label.size = 0.5)#, main = paste(v.tissues[s], " / ", names(l.grn_treatment[[s]])[i], " / ", v.domains[d], sep =""))
# 
# 
# # Add a legend
# legend("bottomleft", legend=levels(as.factor(V(network)$carac))  , col = coul , bty = "n", pch=20 , pt.cex = 3, cex = 1.5, text.col=coul , horiz = FALSE, inset = c(0.1, 0.1))
# 
# 
# 
# plot(g, layout = layout, edge.arrow.size=.2, vertex.label=vertex.label,edge.width=E(g)$weight, vertex.label.dist=0.1, vertex.label.color = "black")#, main = paste(v.tissues[s], " / ", names(l.grn_treatment[[s]])[i], " / ", v.domains[d], sep =""))
# 
# 
# 
# 
# new_cols <- c("white", "red", "black")[membership(wc)]
# 
# V(g)$color <- ifelse(names(V(g)) %in% names(v.tf_families), adjustcolor("Orange", alpha.f = .6), adjustcolor("black", alpha.f = .6))
# 
# V(g)$type <- ifelse(names(V(g)) %in% names(v.tf_families), 0, 1)
# 
# E(g)$weight <- df.regulatoryNetwork.selection$rank / max(df.regulatoryNetwork.selection$rank)
# 
# 
# vertex.label <- ifelse(names(V(g)) %in% df.regulatoryNetwork.selection$tf, names(V(g)), "")
# 
# # layout <-layout.fruchterman.reingold(g)
# plot(g, layout = layout, edge.arrow.size=.2, vertex.label=vertex.label,edge.width=E(g)$weight, vertex.label.dist=0.1, vertex.label.color = "black")#, main = paste(v.tissues[s], " / ", names(l.grn_treatment[[s]])[i], " / ", v.domains[d], sep =""))
