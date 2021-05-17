# source("MERIT_DifferentialExpression.R")
# source("MERIT_DNA.R")

m.de <- readRDS("data/m.de_0502.rds")
m.fc <- readRDS("data/m.fc_0502.rds")
v.series_sets <- readRDS("data/v.series_sets_0502.rds")
v.repeats <- readRDS("data/v.repeats_0502.rds")
v.treatment_buildingblocks <- readRDS("data/v.treatments_0502.rds")

# raw expression => filtered 
m.expression <- as.matrix(read.table("data/GSE69995_re-analyzed_data_matrix.txt", row.names = 1, header = TRUE, sep = "\t", quote = ""))


df.root_selection <- read.csv("data/expMetaRoot.txt", sep = ";", header = FALSE, stringsAsFactors = FALSE)
v.series.root <- unique(df.root_selection$V11)
v.series.root <- v.series.root[v.series.root != "characteristics..Original.series.ID"]

v.samples.root <- unique(df.root_selection$V2)
v.samples.root <- v.samples.root[v.samples.root != "Sample.name"]

v.samples.root <- intersect(colnames(m.expression), v.samples.root)
m.expression.root <- m.expression[,v.samples.root] # expression matrix for random forest 

i.set <- which(colnames(m.fc) %in% v.series.root)

v.series.root <- intersect(colnames(m.de), v.series.root)
m.de.root <- m.de[,v.series.root]

m.fc.root <- m.fc[,i.set] # select from log foldchange matrix

v.treatment_buildingblocks.root <- v.treatment_buildingblocks[i.set]
v.repeats.root <- v.repeats[i.set]



# root specific analysis of the REMAINING gene sets!!! 
# m.cor <- cor(t(m.fc.root))
m.fc.root <- m.fc.root[intersect(rownames(m.fc.root), gns.DE),] # 175 genes 
m.expression.root <- m.expression.root[intersect(rownames(m.fc.root), gns.DE),]



tfs <- intersect(rownames(m.fc.root), names(v.regulators))
tgs <- rownames(m.expression.root)



# v.tf_families[tfs.DE.root]
# m.fc.root was previously
source("utils.R")
m.rf <- compute_randomforest_based_GRN(mat.expression=m.fc.root, k="sqrt", nb.trees=1000, set.regulators = tfs, set.genes = tgs, seed=1234, importance.measure = "impurity", n.cpus = 2)
# saveRDS(m.rf, "m.RF_1000_032621.rds") # 032621 - work with the filtered 380 experiments dataset 
# saveRDS(m.rf, "tmp/m.RF.rds")
# m.rf <- readRDS("data/m.RF_500.rds")

# th = quantile(m.rf, 0.9)
p.rf <- ecdf(m.rf)
df.rf.grn <- melt(m.rf)
names(df.rf.grn) <- c("TF", "TG", "val")
df.rf.grn <- subset(df.rf.grn, df.rf.grn$val > 0)
df.rf.grn["val"] <- p.rf(df.rf.grn$val)

df.rsk1 = subset(df.rf.grn, df.rf.grn$TG == rsk1)


# df.rsk1 = subset(df.grn, df.grn$TG == rsk1)



# library(parmigene)
# 
# mi <- knnmi.all(m.fc.root, k=3, noise=1e-12)
# m.clr <- clr(mi) 
# clr.grn <- m.clr[tfs, ]
# library(reshape2)
# df.clr.grn <- melt(clr.grn)
# names(df.clr.grn) <- c("TF", "TG", "val")
# df.clr.grn <- subset(df.clr.grn, df.clr.grn$val > 0)
# 
# th = quantile(m.clr, 0.95)


# wisdowm of the crowds ensemble network / BP network
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3512113/pdf/nihms420148.pdf
# 
# rsk1 = "AT2G26290"
# 
# df.rsk1 = subset(df.clr.grn, df.clr.grn$TG == rsk1)
# 
# 
# mi[tfs, rsk1]
# 
# 
# install.packages("BiocManager")
# BiocManager::install("WGCNA") 
# 
# 
# 
# # run the MFA function from the FactoMineR package
# r.mfa <- FactoMineR::MFA(
#   t(rbind(x1,x2,x3)), # binding the omics types together
#   c(dim(x1)[1], dim(x2)[1], dim(x3)[1]), # specifying the dimensions of each
#   graph=FALSE)
# 
# 
# 
# 
# 
# 
# 
# # correlation network 
# m.cor <- cor(t(m.fc.root))
# m.cor <- abs(m.cor)
# th <- quantile(m.cor, 0.95)
# m.cor[m.cor < th] = 0
# 
# 
# df.cor <- melt(m.cor[tfs, ])
# names(df.cor) <- c("TF", "TG", "val")
# subset(df.cor, df.cor$TG == rsk1)
# 


# 
# TF        TG       val
# 607 AT1G77200 AT2G26290 0.0700513
# 608 AT3G16280 AT2G26290 2.9521025
# 609 AT3G56980 AT2G26290 0.4988810
# 611 AT4G17900 AT2G26290 0.5030344
# 612 AT4G35700 AT2G26290 1.0983832
# 614 AT5G54230 AT2G26290 1.0629806
# 616 AT5G65640 AT2G26290 0.5337378

# https://rdrr.io/bioc/MODA/man/WeightedModulePartitionSpectral.html

# subspace clustering / module detection
# TODO: do a BP analysis on modules

# 
# 
# library(igraph)
# g=graph.edgelist(as.matrix(df.regulatoryNetwork.selection[,1:2]))
# V(g)$size <-  ifelse(names(V(g)) %in% df.regulatoryNetwork.selection$tf, 3, 2)
# V(g)$color <- ifelse(names(V(g)) %in% names(v.tf_families), adjustcolor("Orange", alpha.f = .6), adjustcolor("black", alpha.f = .6))
# V(g)$type <- ifelse(names(V(g)) %in% names(v.tf_families), 0, 1)
# E(g)$color <- ifelse(df.regulatoryNetwork.selection$mode_of_regulation == "activation", "red", "green")
# E(g)$weight <- df.regulatoryNetwork.selection$rank / max(df.regulatoryNetwork.selection$rank)
# 
# 
# vertex.label <- ifelse(names(V(g)) %in% df.regulatoryNetwork.selection$tf, names(V(g)), "")
# 
# # layout <-layout.fruchterman.reingold(g)
# plot(g, layout = layout, edge.arrow.size=.2, vertex.label=vertex.label,edge.width=E(g)$weight, vertex.label.dist=0.1, vertex.label.color = "black")#, main = paste(v.tissues[s], " / ", names(l.grn_treatment[[s]])[i], " / ", v.domains[d], sep =""))
# 
# 
# # TODO: in the dna based set, highlight the RSK1 
# library(igraph)
# 
# # Create a graph adjacency based on correlation distances between genes in  pairwise fashion.
# g <- graph.adjacency(
#   as.matrix(as.dist(m.cor)),
#   mode="undirected",
#   weighted=TRUE,
#   diag=FALSE,
# )
# 
# # TODO: highlight the transcription factors, along with 
# plot(g,  vertex.size=1, 
#      vertex.label.dist=-0.5,
#      vertex.label.color="black",
#      asp=FALSE,
#      vertex.label.cex=0.6)
# 
# # Simplfy the adjacency object
# g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE)
# 
# # Colour negative correlation edges as blue
# E(g)[which(E(g)$weight<0)]$color <- "darkblue"
# 
# # Colour positive correlation edges as red
# E(g)[which(E(g)$weight>0)]$color <- "darkred"
# 
# # Convert edge weights to absolute values
# E(g)$weight <- abs(E(g)$weight)
# 
# # Change arrow size
# # For directed graphs only
# #E(g)$arrow.size <- 1.0
# 
# # Remove edges below absolute Pearson correlation 0.8
# g <- delete_edges(g, E(g)[which(E(g)$weight<0.8)])
# 
# # Remove any vertices remaining that have no edges
# g <- delete.vertices(g, degree(g)==0)
# 
# # Assign names to the graph vertices (optional)
# V(g)$name <- V(g)$name
# 
# # Change shape of graph vertices
# V(g)$shape <- "sphere"
# 
# # Change colour of graph vertices
# V(g)$color <- "skyblue"
# 
# # Change colour of vertex frames
# V(g)$vertex.frame.color <- "white"
# 
# # Scale the size of the vertices to be proportional to the level of expression of each gene represented by each vertex
# # Multiply scaled vales by a factor of 10
# scale01 <- function(x){(x-min(x))/(max(x)-min(x))}
# # vSizes <- (scale01(apply(estrogenMainEffects, 1, mean)) + 1.0) * 10
# 
# # Amplify or decrease the width of the edges
# edgeweights <- E(g)$weight * 2.0
# 
# # Convert the graph adjacency object into a minimum spanning tree based on Prim's algorithm
# mst <- mst(g, algorithm="prim")
# 
# # Plot the tree object
# plot(
#   mst,
#   layout=layout.fruchterman.reingold,
#   edge.curved=TRUE,
#   vertex.size=1, #vSizes,
#   vertex.label.dist=-0.5,
#   vertex.label.color="black",
#   asp=FALSE,
#   vertex.label.cex=0.6,
#   edge.width=edgeweights,
#   edge.arrow.mode=0,
#   main="My first graph")
# 
# 
# 
# 
# 
# 
# 
# 
# 










#     mat.GE.Ath <- as.matrix(mat.GE.Ath)
# library(WGCNA)
# library(genefilter)
# library(parmigene)

# # two gene expression sets - root related
# if(FALSE){
#   
#     mat.GE.Ath <- read.table("gene_expression/GSE69995_re-analyzed_data_matrix.txt", row.names = 1, header = TRUE, sep = "\t", quote = "")
#     mat.GE.Ath <- as.matrix(mat.GE.Ath)
#     map.GE.MicroArray <- read.csv("gene_expression/tpj13175-sup-0003-TableS2.csv", stringsAsFactors = FALSE)
#     
#     v.conditions <- paste(map.GE.MicroArray$Series.Title, map.GE.MicroArray$Sample.Title, sep = ";")
#     colnames(mat.GE.Ath) <- v.conditions
#     
#     m.rootSet <- mat.GE.Ath[,which(map.GE.MicroArray$Tissue == "root")]
#     v.gns <- rownames(m.rootSet)
#   
#     m.pcc <- cor(t(m.rootSet))
#     saveRDS(m.pcc, "m.pcc_big.rds")
#     
#     m.mi <- parmigene::knnmi.all(m.rootSet) # capturing nonlinear expression relationships
#     saveRDS(m.mi, "m.mi_big.rds")
#     
# }else{
