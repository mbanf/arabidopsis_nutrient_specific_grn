
# main, gene regulatory network

library(reshape2)
library(ggplot2)



source("nutrition_differential_expression.R")
source("transcriptional_regulators.R")



##### gsea of gene clusters ### 

for(j in 1:c_group){
  v.gns <- names(which(ct == j))
  go_gsea(v.gns, th = 0.1, ontology = "BP")
  go_gsea(v.gns, th = 0.1, ontology = "MF")
}


## all genes gsea 
go_gsea(gns.DE, th = 0.1, ontology = "BP")
go_gsea(gns.DE, th = 0.1, ontology = "MF")

gns <- names(V(g))
specs <- gn_spec[gns]
c_group <- sort(unique(specs))

for(s in 1:length(code_specificity)){
  
  gn = names(which(specs == s))
  go_gsea(gn, th = 0.05, title = paste("Biological process of", length(gn), "genes specific to", names(code_specificity)[s]), ontology = "BP")
  go_gsea(gn, th = 0.05, title = paste("Molecular function of", length(gn), "genes specific to", names(code_specificity)[s]),  ontology = "MF")
  
}



cols_specificity[gn_spec[]]


### gene regulatory network  ###

### nutrition differential expression based network ### 

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

### TODO: change node types - regulators and targets 





### network

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




message("Runnning gene expression based inference to support putative links of differential expression")

# perform analysis for all putative links with differential gene expression

tfs.DE <- names(v.regulators)
tgs <- rownames(m.fc.root)
tgs <- rownames(m.expression.root)

#### GENIE3 ####
message("running GENIE3")
tfs.DE.root <- (intersect(tfs.DE, rownames(m.fc.root)))
tfs.DE.root <- (intersect(tfs.DE, rownames(m.expression.root)))

# v.tf_families[tfs.DE.root]
# m.fc.root was previously
m.rf <- compute_randomforest_based_GRN(mat.expression=m.fc.root, k="sqrt", nb.trees=1000, set.regulators = tfs.DE.root, set.genes = tgs, seed=1234, importance.measure = "impurity", n.cpus = 4)
# saveRDS(m.rf, "m.RF_1000_032621.rds") # 032621 - work with the filtered 380 experiments dataset 
saveRDS(m.rf, "m.RF_1000.rds")
m.rf <- readRDS("data/m.RF_500.rds")

## test comparison 

df.regulatoryNetwork <- as.data.frame(as.table(m.rf), stringsAsFactors = FALSE) 
names(df.regulatoryNetwork)[1:3] <- c("tf","target","rf")
df.regulatoryNetwork <- subset(df.regulatoryNetwork, df.regulatoryNetwork$rf > 0)
df.regulatoryNetwork <- subset(df.regulatoryNetwork, as.character(df.regulatoryNetwork$tf) != as.character(df.regulatoryNetwork$target))












# setwd("/home/mbanf/Documents/Computational_Biology/Projects/Hatem/")

#setwd("/shared/Labs/Rhee Lab/Everyone/Michael/Hatem-FePNetwork/")
# setwd("Desktop/At_Root_GeneExp/")
rm(list=ls()) # clear workspace 


# CALC
# setwd("/home/mbanf/Documents/Computational_Biology/Projects/GRACE_release/")

# source("http://bioconductor.org/biocLite.R") 
# biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore")) 
# install.packages("WGCNA")


message("loading general pre-curated differential expression and dna binding data")



# selected by hatem (root specific stress)
message("loading root specific annotations and datasets (gene ontology, anova results)")


#### Differential expression data
#df.diffexp <- read.csv("diffexp_highConfidence.csv", stringsAsFactors = FALSE)


#df.diffexp["AGI"] <- v.map[as.character(df.diffexp$Cluster)]

# read hatem's anova datasets



# 
# exps <- c("-P+Fe (3h)", "-P+Fe (6h)", "-P+Fe (9h)",  "+P-Fe (3h)", "+P-Fe (6h)", "+P-Fe (9h)", "-P-Fe (3h)", "-P-Fe (6h)", "-P-Fe (9h)")
# 
# 
# # cor_0PF_P0F <- apply(df.geneExp, 1, function(m) { cor(m[c(1,4,7)] , m[c(2,5,8)])})
# # cor_0PF_0P0F <- apply(df.geneExp, 1, function(m) { cor(m[c(1,4,7)] , m[c(3,6,9)])})
# # cor_P0F_0P0F <- apply(df.geneExp, 1, function(m) { cor(m[c(2,5,8)] , m[c(3,6,9)])})
# # 
# # cor(cor_0PF_P0F, cor_0PF_0P0F)
# # cor(cor_0PF_P0F, cor_P0F_0P0F)
# # 
# wilcox.test(v.0PF, v.P0F)
# wilcox.test(v.0PF, v.0P0F)
# wilcox.test(v.P0F, v.0P0F)


####

# - replace with updated list #

# load additional evidences 

# conservation

####
# df.regulatoryNetwork.selection <- subset(df.regulatoryNetwork, df.regulatoryNetwork$TF %in% rownames(m.anova.set))

# run GENIE3 #

# subset the dna binding set to anova (either regulator or target needs to be differentially expressed)
# df.dna.set <- subset(df.dna, df.dna$TF %in% rownames(m.anova.set) | df.dna$Target %in% rownames(m.anova.set)) # unique
# 
# df.dna.set["DE.TF"] <- (df.dna.set$TF %in% rownames(m.anova.set))
# df.dna.set["DE.TG"] <- (df.dna.set$Target %in% rownames(m.anova.set))
# 
# df.dna.DE <- subset(df.dna.set, df.dna.set$DE.TF == TRUE & df.dna.set$DE.TG == TRUE) # subset to 237 

# subset(df.dna.DE, df.dna.DE$Target == "AT2G26290")


#tfs.DE <- df.dna.DE$TF
#tgs.DE <- df.dna.DE$Target

# message("running eta squared")
# strt<-Sys.time()
# cl<-makeCluster(min(length(tfs.DE),n.cpus))
# registerDoParallel(cl)
# l.regulatoryNetwork <- foreach(r = 1:length(tfs.DE), .packages=c("fume", "reshape2")) %dopar% { 
# 
#   tf <- tfs.DE[r]
#   
#   #df.dna.r <- subset(df.dna.DE, df.dna.DE[,1] == tf)
#   #tgs <- intersect(as.character(df.dna.r[,2]), v.gns)
#   
#   tgs <- v.gns
#   
#   df.regulatoryNetwork <- c()
#   l.treatments <- list()
#   
#   if(tf %in% v.gns){
#     
#     tgs.links <- tgs
#     # tgs.links <- names(which(m.gc[tf, tgs] >= th.min_overlap))
#     
#     if(length(tgs.links) > 0){
#       
#       #df.regulatoryNetwork <- subset(df.dna.r, df.dna.r$Target %in% tgs.links)
#       
#       df.regulatoryNetwork <- data.frame(tf = rep(tf, length(tgs.links)), fam = rep(v.regulators[tf], length(tgs.links)),  target = tgs.links)
#       
#       df.regulatoryNetwork["eta"] <- 0
#       df.regulatoryNetwork["p.val"] <- 1
#       df.regulatoryNetwork["mechanism"] <- 0
#       df.regulatoryNetwork["scc"] <- 0
#       
#       pb <- txtProgressBar(min = 0, max = length(tgs.links), style = 3)  
#       
#       for(j in 1:length(tgs.links)){
#         
#         setTxtProgressBar(pb, j)
#         
#         expr.tf <- m.fc.root[tf, ]
#         expr.tf.inv <- - m.fc.root[tf, ] 
#         expr.tg <- m.fc.root[tgs.links[j], ]
#         
#         # copy the function of eta squared
#         res_orig <- eta_squared_inference(expr.tf, expr.tg, v.treatment_buildingblocks.root, v.repeats.root) 
#         res_inv <- eta_squared_inference(expr.tf.inv, expr.tg, v.treatment_buildingblocks.root, v.repeats.root) 
#         
#         # use eta to select
#         if (res_orig[2] < res_inv[2]) { # if (res_orig[2] < res_inv[2]) {
#           
#           df.regulatoryNetwork$eta[j] <- res_orig[1]
#           df.regulatoryNetwork$p.val[j] <- res_orig[2]
#           df.regulatoryNetwork$mechanism[j] <- "activation"
#           df.regulatoryNetwork$scc[j] <- cor(expr.tf, expr.tg, method = "spearman")
#           
#         } else {
#           
#           df.regulatoryNetwork$eta[j] <- res_inv[1]
#           df.regulatoryNetwork$p.val[j] <- res_inv[2]
#           df.regulatoryNetwork$mechanism[j] <- "repression"
#           df.regulatoryNetwork$scc[j] <- cor(expr.tf, expr.tg, method = "spearman")
#           
#         }
#       }
#       close(pb)
#       
#     }
#   }
#   
#   df.regulatoryNetwork["fdr"] <- p.adjust(df.regulatoryNetwork$p.val,"fdr") 
#   #df.regulatoryNetwork <- subset(df.regulatoryNetwork, df.regulatoryNetwork$p.val <= 0.05)
#   df.regulatoryNetwork
# }
# stopCluster(cl)
# print(Sys.time()-strt)
# 
# df.regulatoryNetwork <- do.call(rbind, lapply(l.regulatoryNetwork, function(m) return(m)))
# #df.regulatoryNetwork["nr.link"] <- seq(1:nrow(df.regulatoryNetwork))
# df.regulatoryNetwork <- df.regulatoryNetwork[order(-df.regulatoryNetwork$eta),]


#### GENIE3 ####
message("running GENIE3")
tfs.DE.root <- (intersect(tfs.DE, rownames(m.fc.root)))
tfs.DE.root <- (intersect(tfs.DE, rownames(m.expression.root)))

# v.tf_families[tfs.DE.root]
# m.fc.root was previously
m.rf <- compute_randomforest_based_GRN(mat.expression=m.fc.root, k="sqrt", nb.trees=1000, set.regulators = tfs.DE.root, set.genes = tgs, seed=1234, importance.measure = "impurity", n.cpus = 4)
# saveRDS(m.rf, "m.RF_1000_032621.rds") # 032621 - work with the filtered 380 experiments dataset 
saveRDS(m.rf, "m.RF_1000.rds")
m.rf <- readRDS("data/m.RF_500.rds")

## test comparison 

df.regulatoryNetwork <- as.data.frame(as.table(m.rf), stringsAsFactors = FALSE) 
names(df.regulatoryNetwork)[1:3] <- c("tf","target","rf")
df.regulatoryNetwork <- subset(df.regulatoryNetwork, df.regulatoryNetwork$rf > 0)
df.regulatoryNetwork <- subset(df.regulatoryNetwork, as.character(df.regulatoryNetwork$tf) != as.character(df.regulatoryNetwork$target))

# df.regulatoryNetwork["eta.norm"] <- df.regulatoryNetwork$eta / max(df.regulatoryNetwork$eta)
#df.regulatoryNetwork["rf.norm"] <- df.regulatoryNetwork$rf / max(df.regulatoryNetwork$rf)

## add meta annotation ##

# diff exp - for the same experiment (1,2) 
tfs.DE <- as.character(unique(df.regulatoryNetwork$tf))

df.regulatoryNetwork["dna_binding.conservation"] <- 0
df.regulatoryNetwork["dna_binding.dapSeq"] <- 0

df.regulatoryNetwork.meta <- c()

for(j in 1:length(tfs.DE)){
  
  df.regulatoryNetwork.j <- subset(df.regulatoryNetwork, df.regulatoryNetwork$tf ==  tfs.DE[j])
  
  # conservation 
  if(tfs.DE[j] %in% df.dna.conservation$TF){
    df.dna.conservation.j <- subset(df.dna.conservation, df.dna.conservation$TF == tfs.DE[j])
    i.set <- which(df.regulatoryNetwork.j$target %in% df.dna.conservation.j$Target)  
    if(length(i.set) > 0)
      df.regulatoryNetwork.j$dna_binding.conservation[i.set] <- 1
  }
  
  # dna binding 
  if(tfs.DE[j] %in% df.dna$TF){
    df.dna.j <- subset(df.dna, df.dna$TF == tfs.DE[j])
    i.set <- which(df.regulatoryNetwork.j$target %in% df.dna.j$Target)  
    if(length(i.set) > 0)
      df.regulatoryNetwork.j$dna_binding.dapSeq[i.set] <- 1
  }
  
  df.regulatoryNetwork.meta <- rbind(df.regulatoryNetwork.meta, df.regulatoryNetwork.j)
  
}

#df.regulatoryNetwork.selection <- subset(df.regulatoryNetwork.meta, rowSums(df.regulatoryNetwork.meta[,9:17]) > 0)
#df.regulatoryNetwork.selection_dna <- subset(df.regulatoryNetwork.selection, df.regulatoryNetwork.selection$dna_binding.dapSeq == 1)


# construct unsupervised integration 

# p.eta <- ecdf(df.regulatoryNetwork$eta)
p.rf <- ecdf(df.regulatoryNetwork$rf)

# df.regulatoryNetwork["p.eta"] <- p.eta(df.regulatoryNetwork$eta)
df.regulatoryNetwork["p.rf"] <- p.rf(df.regulatoryNetwork$rf)

df.regulatoryNetwork["rank"] <- df.regulatoryNetwork$coDiffExp + df.regulatoryNetwork$dna_binding.dapSeq # + df.regulatoryNetwork$p.rf
df.regulatoryNetwork <- df.regulatoryNetwork[order(-df.regulatoryNetwork$rank),]

#df.dna.set <- subset(df.dna, df.dna$TF %in% rownames(m.anova.set) | df.dna$Target %in% rownames(m.anova.set)) # unique
df.regulatoryNetwork.selection <- subset(df.regulatoryNetwork, df.regulatoryNetwork$rank > 1.95) # this way gene expression really assists

# compute correlation matrix for mode of regulation
df.regulatoryNetwork.selection["mode_of_regulation"] <- m.cor[cbind(df.regulatoryNetwork.selection$tf, df.regulatoryNetwork.selection$target)]
df.regulatoryNetwork.selection$mode_of_regulation <- ifelse(df.regulatoryNetwork.selection$mode_of_regulation >= 0, "activation", "repression")


#table(df.regulatoryNetwork.selection$tf)

#table(df.regulatoryNetwork.selection$target)

# v.tf_families[v.regulators]
# 
# table(v.tf_families[unique(df.regulatoryNetwork.selection$tf)])

df.regulatoryNetwork.selection["fam"] <- v.regulators[as.character(df.regulatoryNetwork.selection$tf)]
df.regulatoryNetwork.selection <- df.regulatoryNetwork.selection[,c(1,13,2:12)]
df.regulatoryNetwork.selection <- df.regulatoryNetwork.selection[,-4]

# write.csv(df.regulatoryNetwork.selection, "df.regulatoryNetwork.selection.csv", row.names = FALSE)

message("...finished")


df.regulatoryNetwork.selection <- read.csv("data/df.regulatoryNetwork.selection.csv")
#table(v.tf_families[unique(df.regulatoryNetwork.selection$tf)])

###
# 
# df.dna.DE["eta"] <- NA
# df.dna.DE["p.val"] <- NA
# df.dna.DE["mechanism"] <- NA
# df.dna.DE["scc"] <- NA
# 
# for(i in 1:nrow(df.dna.DE)){
#   
#   df.regulatoryNetwork.i <- unique(subset(df.regulatoryNetwork, df.regulatoryNetwork$TF == df.dna.DE$TF[i] & df.regulatoryNetwork$Target == df.dna.DE$Target[i]))
#   
#   if(nrow(df.regulatoryNetwork.i)){
#     df.dna.DE$eta[i] <- df.regulatoryNetwork.i$eta
#     df.dna.DE$p.val[i] <- df.regulatoryNetwork.i$p.val
#     df.dna.DE$mechanism[i] <- df.regulatoryNetwork.i$mechanism
#     df.dna.DE$scc[i] <- df.regulatoryNetwork.i$scc
#   }
#   
# }
# 
# 
# ####

message("plot detailed regulation heatmaps - for diff exp targets")



# m.test <- as.matrix(df.geneExp.set[tfs.DE,c(1,4,7,2,5,8,3,6,9), drop = FALSE])
# colnames(m.test) <- exps
# m.sig <- m.anova.set[tfs.DE,]
# 
# p.tf <- plot_ratios(cormat=m.test, v.title = "Transcription", v.name = "", v.max = max(m.test), v.min = min(m.test))
# 
# 
# # tfs.DE <- unique(df.dna.DE$TF)
# 
# m.test <- as.matrix(df.geneExp.set[tfs.DE,c(1,4,7,2,5,8,3,6,9)])
# colnames(m.test) <- exps
# m.sig <- m.anova.set[tfs.DE,]
# 
# m.sig[m.sig == 1] <- "*"
# m.sig[m.sig == 0] <- ""
# 
# p.tfs <- plot_ratios(cormat=m.test, m.value=m.sig, v.title = "Transcription Factors", v.name = "logFC", v.max = max(m.test), v.min = min(m.test))


# TODO: target regulation by individual TF plot 

tfs.DE.selection <-  unique(df.regulatoryNetwork.selection$tf)

l.p.tgs <- vector(mode = "list", length = length(tfs.DE.selection))

for(r in 1:length(tfs.DE.selection)){
  
  df.dna.DE.r <- subset(df.regulatoryNetwork.selection, df.regulatoryNetwork.selection$tf == tfs.DE.selection[r])
  df.dna.DE.r <- subset(df.dna.DE.r, df.dna.DE.r$tf != df.dna.DE.r$target)

  tgs.r <- intersect(rownames(m.anova.set), df.dna.DE.r$target)
  df.dna.DE.r <- subset(df.dna.DE.r, df.dna.DE.r$target %in% tgs.r)
  
  m.test <- as.matrix(df.geneExp.set[tgs.r,c(1,4,7,2,5,8,3,6,9)], drop = FALSE)
  colnames(m.test) <- exps
  m.sig <- m.anova.set[tgs.r,]
  m.sig[m.sig == 1] <- "*"
  m.sig[m.sig == 0] <- ""
  
  v.sig.tf <- m.anova.set[tfs.DE.selection[r],]
  v.sig.tf[v.sig.tf == 1] <- "*"
  v.sig.tf[v.sig.tf == 0] <- ""
  
  if(dim(m.test)[1] > 1){
    
    dd.col <- as.dendrogram(hclust(dist((m.test))))
    row.col <- order.dendrogram(dd.col)
  
    require(ggdendro)
    #ddata_y <- dendro_data(t(dd.row)
    
    m.test <- m.test[row.col,]
    m.sig <- m.sig[row.col,]
    
  }
  
  v.exp.tf <- df.geneExp.set[tfs.DE.selection[r],c(1,4,7,2,5,8,3,6,9)]
  names(v.exp.tf) <- exps
  
  m.test <- rbind(m.test, v.exp.tf)
  m.sig <- rbind(m.sig, v.sig.tf)
  
  m.test <- as.matrix(m.test)
  m.sig <- as.matrix(m.sig)
  
  max.line <- length(df.dna.DE.r$target) + 0.5 + 1
  p.tgs <- plot_ratios(cormat=m.test, m.value=m.sig, v.title = paste("Targets of ", tfs.DE.selection[r], sep =""), v.name = "logFC", v.max = max(m.test), v.min = min(m.test), max.line = max.line)
  
  l.p.tgs[[r]] <- p.tgs 
  
  # integrate TF 
  # integrate tgs
  # integrate equal ordered additional evidence
  
  #grid.arrange(l.p.tgs[[1]],l.p.tgs[[2]]
  
}
# 
# p1 <- plot_ratios(cormat=m.test[,1:3], v.title = "", v.name = "", v.max = max(m.test), v.min = min(m.test))
# p2 <- plot_ratios(cormat=m.test[,4:6], v.title = "", v.name = "", v.max = max(m.test), v.min = min(m.test))
# p3 <- plot_ratios(cormat=m.test[,7:9], v.title = "", v.name = "", v.max = max(m.test), v.min = min(m.test))


library(gridExtra)
library(grid)
grid.arrange(l.p.tgs[[1]], l.p.tgs[[2]], l.p.tgs[[3]], l.p.tgs[[4]], l.p.tgs[[5]], l.p.tgs[[6]], l.p.tgs[[7]], l.p.tgs[[8]], l.p.tgs[[9]], l.p.tgs[[10]], ncol = 5,  top = "log foldChange differential expression patterns of putative regulators and targets across all three experimental time series") # plot the three expression b


message("plot the entire gene expression map")


library("ggplot2")
library("ggdendro")
library("reshape2")

df <- df.geneExp.set[gns.DE,c(1,4,7,2,5,8,3,6,9)]
df <- na.omit(df)
m.test <- as.matrix(df, drop = FALSE)
colnames(m.test) <- exps
heatmap((m.test), Colv=F, scale='none')



message("export gene regulatory network graph for cytoscape")

# red would be regulators
# blue would be others


# red and green arrows (activation and repression)
# diamond and circles
# edge rank on link

library(igraph)
g=graph.edgelist(as.matrix(df.regulatoryNetwork.selection[,1:2]))
V(g)$size <-  ifelse(names(V(g)) %in% df.regulatoryNetwork.selection$tf, 3, 2)
V(g)$color <- ifelse(names(V(g)) %in% names(v.tf_families), adjustcolor("Orange", alpha.f = .6), adjustcolor("black", alpha.f = .6))
V(g)$type <- ifelse(names(V(g)) %in% names(v.tf_families), 0, 1)
E(g)$color <- ifelse(df.regulatoryNetwork.selection$mode_of_regulation == "activation", "red", "green")
E(g)$weight <- df.regulatoryNetwork.selection$rank / max(df.regulatoryNetwork.selection$rank)


vertex.label <- ifelse(names(V(g)) %in% df.regulatoryNetwork.selection$tf, names(V(g)), "")

# layout <-layout.fruchterman.reingold(g)
plot(g, layout = layout, edge.arrow.size=.2, vertex.label=vertex.label,edge.width=E(g)$weight, vertex.label.dist=0.1, vertex.label.color = "black")#, main = paste(v.tissues[s], " / ", names(l.grn_treatment[[s]])[i], " / ", v.domains[d], sep =""))

#bg <- bipartite.projection(g)
# l <- layout_as_bipartite(g, types = NULL, hgap = 1, vgap = 1, maxiter = 100)
# 
# # ertex.size=2, vertex.label.dist=1,
# # 
# # 
# # 
# # exportNetworkToCytoscape {WGCNA}	R Documentation
# # Export network to Cytoscape
# # 
# # Description
# # 
# # This function exports a network in edge and node list files in a format suitable for importing to Cytoscape.
# # 
# # Usage
# 
# library(WGCNA)
# exportNetworkToCytoscape(adjMat, edgeFile = NULL, nodeFile = NULL, weighted = TRUE, threshold = 0.5,
#                          nodeNames = NULL, altNodeNames = NULL, nodeAttr = NULL, includeColNames = TRUE)
# 
# 
# 
# message("perform module detection / GO / ")
# 
# df.geneontology.selection <- subset(df.geneontology, df.geneontology$AGI %in% df.regulatoryNetwork.selection$target)
# 
# 
# # GO 
# tb.go <- table(unlist(sapply(df.geneontology$GO.Biological.Process.Term, function(m) strsplit(m, "; "))))
# tb.go.bg <- tb.go[tb.go >= 3]
# 
# tb.go <- table(unlist(sapply(df.geneontology$GO.Cellular.Component.Term, function(m) strsplit(m, "; "))))
# tb.go.bg <- tb.go[tb.go >= 3]
# 
# tb.go <- table(unlist(sapply(df.geneontology$GO.Molecular.Function.Term, function(m) strsplit(m, "; "))))
# tb.go.bg <- tb.go[tb.go >= 3]
# 
# v.go <- subset(df.geneontology, df.geneontology$AGI %in% df.regulatoryNetwork.selection$target)$GO.Biological.Process.Term
# tb.go <- table(unlist(sapply(v.go, function(m) strsplit(m, "; "))))
# tb.go <- tb.go[tb.go >= 3]
# 
# v.go <- subset(df.geneontology, df.geneontology$AGI %in% df.regulatoryNetwork.selection$target)$GO.Cellular.Component.Term
# tb.go <- table(unlist(sapply(v.go, function(m) strsplit(m, "; "))))
# tb.go <- tb.go[tb.go >= 3]
# 
# v.go <- subset(df.geneontology, df.geneontology$AGI %in% df.regulatoryNetwork.selection$target)$GO.Molecular.Function.Term
# tb.go <- table(unlist(sapply(v.go, function(m) strsplit(m, "; "))))
# tb.go <- tb.go[tb.go >= 3]
# 
# 
# 
# 
# 
# df.enrichment <- c()
# 
# for(i in 1:length(tb.go)){
#   
#   hitInSample <- tb.go[i]
#   sampleSize <- length(v.go) #sum(tb.go)
#   
#   hitInPop <- tb.go.bg[names(tb.go)[i]]
#   popSize <- nrow(df.geneontology) # sum(tb.go.bg)
#   
#   failInPop <- popSize - hitInPop
#   
#   fc <- ((hitInSample / sampleSize) / (hitInPop / popSize))
#   
#   if(fc > 1){
#     
#     p.val <- phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE)
#     
#     df.enrichment <- rbind(df.enrichment, data.frame(term = names(tb.go)[i], n.genes = hitInSample,  foldChange = fc, pvalue = p.val ))
#   }
# }
# 
# df.enrichment <- subset(df.enrichment, df.enrichment$pvalue <= 0.05)
# 
# write.csv(df.enrichment, "df.enrichment_CC.csv", row.names = FALSE)
# 
# 
# 
# 
# 
# 
# print("Gene ontology (biological process) - differential expressed genes")
# print(tb.go)
# 
# ## 
# 
# v.go <- subset(df.geneontology, df.geneontology$AGI %in% v.gns.anova[!v.gns.anova %in% v.gns.diff])$GO.Biological.Process.Term
# tb.go <- table(unlist(sapply(v.go, function(m) strsplit(m, "; "))))
# 
# th.min <- 2
# tb.go <- tb.go[tb.go >= th.min]
# 
# df.geneontology.selection$GO.Biological.Process.Term
# 
# tb.BP <- table(df.geneontology.selection$GO.Biological.Process.Term)
# tb.BP <- tb.BP[tb.BP >=3]
# 
# p.treatment <- rep(1, length(v.treatments.l))
# names(p.treatment) <- names(v.treatments.l)
# 
# for(i in 1:length(v.treatments.l)){
#   
#   hitInSample <- v.treatments.l[i]
#   sampleSize <- sum(v.treatments.l)
#   
#   hitInPop <- l.treatments.tissues[[j]][names(v.treatments.l)[i]]
#   popSize <- sum(l.treatments.tissues[[j]])
#   
#   failInPop <- popSize - hitInPop
#   
#   fc <- ((hitInSample / sampleSize) / (hitInPop / popSize))
#   
#   if(fc > 1 & hitInSample >= th.min.samples)
#     p.treatment[i] <- phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE)
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
# 
# 
# 
# 
# ## make plot
# p <- ggheatmap(X)
# 
# ## display plot
# ggheatmap.show(p)
# 
# 
# #subset(df.dna.DE, df.dna.DE$Target == "AT2G26290")
# #subset(df.dna.DE, df.dna.DE$Target == "AT4G01140")
# 
# #source("http://bioconductor.org/biocLite.R") 
# #biocLite("tigre") 
# library(tigre)
# #biocLite("gprege")
# library(gprege)
# 
# ## try http:// if https:// URLs are not supported
# #source("https://bioconductor.org/biocLite.R")
# #biocLite("timecourse")
# #library(timecourse)
# 
# 
# 
# 
# 
# m <- m.anova[intersect(tfs.DE, rownames(m.anova)),c(4,5,6,49,50,51,87,88,89,109,110,111),drop = FALSE]
# idx.regs <- which(rowSums(m) > 0)
# 
# exps <- c("-P+Fe", "+P-Fe", "-P-Fe")
# expSets <- vector(mode = "list", length = 3)
# expSets[[1]] <- c(1,2,3)
# expSets[[2]] <- c(4,5,6)
# expSets[[3]] <- c(7,8,9)
# 
# for(i in 1:3){
#   
#   print(exps[i])
#   
#   idx.gns <- which(rowSums(m.anova.set[,expSets[[i]]]) > 0)
#   
#   print(length(idx.gns)) # number of diff exp genes
#   
#   # highlight 
#   m.anova.set.i <- m.anova.set[,expSets[[i]]] # 
#   
#   df.geneExp.set.i <- df.geneExp.set[,expSets[[i]]]#[!is.na(df.geneExp.set[,expSets[[i]]])]
#   
#   # cluster the profiles - up->down, down-> up, 
#   #  m.tmp <- kmeans(df.geneExp.set[,expSets[[i]]], 4)
#   
#   # expression fold change profiles
#   plot(as.numeric(df.geneExp.set[idx.gns[1],expSets[[i]]]), type = "l", ylim = c(min(df.geneExp.set.i),max(df.geneExp.set.i)))
#   
#   for(j in 2:length(idx.gns)){
#     lines(as.numeric(df.geneExp.set[idx.gns[j],expSets[[i]]]))
#     #df.geneExp.noLog.set[idx.gns[j],expSets[[i]]]  
#   }
#   
#   # filter for specific terms
#   terms <- c("root", "DNA binding", "phosphate", "iron")
#   
#   go.set <- subset(df.geneontology, df.geneontology$AGI %in% names(idx.gns))
#   
#   # general genes with TF annotation
#   tf.set <- subset(go.set, go.set$GO.Molecular.Function.Term %in% c("DNA binding")) # additional binding support ? / coexpression support?
#   
#   # binding support analysis
#   idx.tfs <- which(names(idx.gns) %in% tfs)
#   tf.dna <- names(idx.gns)[idx.tfs]
#   
#   
#   if(length(tf.dna) > 0){
#     
#     for(r in 1:lenth(tf.dna))
#       
#       print(tf.dna[r])
#     
#     plot(as.numeric(df.geneExp.set.i[tf.dna[r],]), type = "l", ylim = c(-3.5,3.5))
#     
#     tgs.r <- names(which(m.dna[tf.dna,] == 1))
#     tgs.r <- intersect(rownames(df.geneExp.set.i), tgs.r)
#     
#     #       
#     #       tf.dna, tgs.r
#     #       
#     m.anova.set[c(tf.dna,tgs.r),]
#     #       
#     #       df.geneExp.set.i[tgs.r,]
#     #       
#     
#     
#     lines(as.numeric(df.geneExp.set.i[tgs.r[1],]), col = "red")
#     lines(as.numeric(df.geneExp.set.i[tgs.r[2],]), col = "blue")
#     
#     #       if(length(tgs.r) > 0){
#     #         for(j in 2:length(idx.gns)){
#     #           lines(as.numeric(df.geneExp.set[idx.gns[j],expSets[[i]]]), col = "")
#     #           #df.geneExp.noLog.set[idx.gns[j],expSets[[i]]]  
#     #         }
#     #       }
#     #       
#     
#     
#   }
#   
#   
#   
#   idx.ifs <- which(tf.set$AGI %in% rownames(m.dna))
#   
#   test <- subset(df.geneontology, grepl(terms, df.geneontology$GO.Biological.Process.Term))
#   
#   df.geneExp.set["AT2G24850",]
#   
#   plot(as.numeric(df.geneExp.set["AT2G24850",expSets[[i]]]), type = "l", ylim = c(-3.5,3.5))
#   
#   
# }
# 
# 
# # install.packages("pgirmess")
# # library(pgirmess)
# # 
# # kruskalmc function in pgirmess
# 
# # test <-  data.frame(exp = c(as.numeric(df.geneExp[1,c(1,14,27)]), as.numeric(df.geneExp[1,c(2,15,28)]), as.numeric(df.geneExp[1,c(3,16,29)])),
# #            set = factor(c(rep("A",3), rep("B",3), rep("C",3))),
# #            number = factor(rep(c(1,2,3),3)))
# 
# # https://www.qbaseplus.com/knowledge/blog/seven-tips-bio-statistical-analysis-gene-expression-data
# # are these values log normalized?
# 
# # two different control condition runs
# 
# # Test 1 - comparison to initial control 
# 
# m.test.global <- c()
# 
# # kruskal better evaluation for correcting 
# m.test.global <- cbind(m.test.global, apply(df.geneExp, 1, function(m) {t.test(as.numeric(m[c(2,15,28)]), as.numeric(m[c(1,14,27)]))$p.value}))
# m.test.global <- cbind(m.test.global, apply(df.geneExp, 1, function(m) {t.test(as.numeric(m[c(6,19,32)]), as.numeric(m[c(1,14,27)]))$p.value}))
# m.test.global <- cbind(m.test.global, apply(df.geneExp, 1, function(m) {t.test(as.numeric(m[c(10,23,36)]), as.numeric(m[c(1,14,27)]))$p.value}))
# 
# for(i in 1:3){
#   
#   m.test.global <- cbind(m.test.global, apply(df.geneExp, 1, function(m) {t.test(as.numeric(m[c(2+i,15+i,28+i)]), as.numeric(m[c(1,14,27)]))$p.value}))
#   m.test.global <- cbind(m.test.global, apply(df.geneExp, 1, function(m) {t.test(as.numeric(m[c(6+i,19+i,32+i)]), as.numeric(m[c(1,14,27)]))$p.value}))
#   m.test.global <- cbind(m.test.global, apply(df.geneExp, 1, function(m) {t.test(as.numeric(m[c(10+i,23+i,36+i)]), as.numeric(m[c(1,14,27)]))$p.value}))
#   
# }
# 
# # Test 2 - comparison to individual control conditions
# m.test.individual <- c()
# for(i in 1:3){
#   
#   m.test.individual <- cbind(m.test.individual, apply(df.geneExp, 1, function(m) {t.test(as.numeric(m[c(2+i,15+i,28+i)]), as.numeric(m[c(2,15,28)]))$p.value}))
#   m.test.individual <- cbind(m.test.individual, apply(df.geneExp, 1, function(m) {t.test(as.numeric(m[c(6+i,19+i,32+i)]), as.numeric(m[c(6,19,32)]))$p.value}))
#   m.test.individual <- cbind(m.test.individual, apply(df.geneExp, 1, function(m) {t.test(as.numeric(m[c(10+i,23+i,36+i)]), as.numeric(m[c(10,23,36)]))$p.value}))
#   
# }
# 
# l.ttests <- vector(mode = "list", length = 2)
# l.ttests[[1]] <- as.matrix(m.test.global)
# l.ttests[[2]] <- as.matrix(m.test.individual)
# names(l.ttests) <- c("global_control", "separate_controls")
# 
# rownames(l.ttests[[1]]) <- v.gnSets
# rownames(l.ttests[[2]]) <- v.gnSets
# 
# # m.anova <- as.matrix(df.anova)
# # rownames(m.anova) <- agi
# 
# # length(which(m.test > 2))
# # length(which(m.test < 0.0001)) # every 1 in 10000
# 
# # notion: no multiple hypothesis correction to retain 
# 
# # pvalues <- m.test[,4]
# # length(which(pvalues < 0.0001))
# # pvalues <- p.adjust(pvalues, "bonferroni")
# # length(which(pvalues < 0.05))
# 
# 
# # hist(df.geneExp[,c(2,15,28)])
# 
# # Test 2 - comparison to individual controls
# 
# 
# # Test 2 - comparison to time complete 
# 
# df.geneExp <- (df.geneExp[1:13] + df.geneExp[14:26] + df.geneExp[27:39])/3
# rownames(df.geneExp) <- v.gnSets
# 
# # compute the standard deviations for the fold changes
# 
# 
# df.geneExp <- df.geneExp[,c(1,2,6,10,3,7,11,4,8,12,5,9,13)]
# 
# df.geneExp.set <- df.geneExp[,-1]
# df.geneExp.set <- df.geneExp.set[,c(1,4,7,10,2,5,8,11,3,6,9,12)]
# df.geneExp.set[,4:6] <- df.geneExp.set[,4:6] - df.geneExp.set[,1]
# df.geneExp.set[,7:9] <- df.geneExp.set[,7:9] - df.geneExp.set[,2]
# df.geneExp.set[,10:12] <- df.geneExp.set[,10:12] - df.geneExp.set[,3]
# 
# df.geneExp.set <- df.geneExp.set[,4:12]
# 
# 
# #### GWAS data 
# l.gwas <- vector(mode = "list", length = 3)
# l.gwas[[1]] <- read.csv("AT-GWAS/fe/Fe.csv", stringsAsFactors = FALSE)
# l.gwas[[2]] <- read.csv("AT-GWAS/Pi/Pi.csv", stringsAsFactors = FALSE)
# l.gwas[[3]] <- read.csv("AT-GWAS/pfe/pfe.csv", stringsAsFactors = FALSE)
# 
# names(l.gwas) <- c("Fe", "Pi", "PFe")
# 
# 
# for(j in 1:length(l.gwas)){
#   
#   l.gwas[[j]]$SNP.relative.position <- gsub("\\..*", "", l.gwas[[j]]$SNP.relative.position)
#   
#   #TF_in_SNP.tmp <- intersect(v.tfs.set[i], l.gwas[[j]]$SNP.relative.position)
#   
# }
# 
# ### DNA binding 
# 
# # # - replace with updated list #
# # df.dna <- read.csv("Yu2016-srep25164-s2.csv", stringsAsFactors = FALSE)
# # df.dna <- df.dna[,c(1,6)]
# # 
# # 
# # tfs <- character()
# # for(i in 1:nrow(df.dna))
# #   tfs <- c(tfs, unlist(strsplit(df.dna$Putative.TF[i], ", ")))
# # 
# # tfs <- unique(tfs)
# # tgs <- unique(df.dna$Target.gene)
# # 
# # m.dna <- matrix(0, nrow = length(tfs), ncol = length(tgs), dimnames = list(tfs, tgs))
# # 
# # for(i in 1:nrow(df.dna)){
# #   
# #   tfs.i <- unlist(strsplit(df.dna$Putative.TF[i], ", "))
# #   tg.i <- df.dna$Target.gene[i]
# # 
# #   m.dna[tfs.i, tg.i] <- 1
# # }
# 
# 
# library(reshape2)
# setNames(melt(m1), c('rows', 'vars', 'values'))
# 
# m.dna <- m.dna[which(rowSums(m.dna) < 2000),]
# 
# tfs <- rownames(m.dna)
# tgs <- colnames(m.dna[,which(colSums(m.dna) > 0)])
# 
# # df.dna_blueprint.m <- subset(df.dna_blueprint, df.dna_blueprint$Target %in% v.enz.m.c)  
# # tb.regs.m.c <- table((unlist(strsplit(df.dna_blueprint.m$TF, ", "))))
# # tb.regs.m.c <- tb.regs.m.c[names(tb.regs.m.c) %in% v.regs.m]
# # tb.regs.m.c <- tb.regs.m.c[tb.regs.m.c >= length(v.enz.m.c)]
# 
# ### 
# v.anova <- read.csv("anova_tests.csv", stringsAsFactors = FALSE, header = FALSE)
# df.anova <- read.csv("anovaResults.csv", stringsAsFactors = FALSE)
# agi <- v.map[as.character(df.anova$Cluster)]
# df.anova <- df.anova[,-1]
# 
# #which(v.map == "AT5G07080")
# 
# #df.anova <- df.anova[,c(49,50,51,87,88,89,109,110,111)]
# #v.anova <- v.anova[,1][c(49,50,51,87,88,89,109,110,111)]
# 
# v.anova <- v.anova[,1]
# names(df.anova) <- v.anova
# 
# m.anova <- as.matrix(df.anova)
# rownames(m.anova) <- agi
# 
# m.anova[m.anova > 0.05] <- 10
# m.anova[m.anova <= 0.05] <- 1
# m.anova[m.anova == 10] <- 0
# 
# # l.anova <- vector(mode = "list", length = 2)
# # l.anova[[1]] <- m.anova
# # l.anova[[2]]
# # 
# # n.DE <- length(which(rowSums(m.anova) > 0))
# 
# # idx <- which(rowSums(m.anova) > 0)
# # m.anova <- m.anova[idx, ]
# 
# 
# # 
# # # pearson correlation set
# # th.pcc.pos <- quantile(m.pcc, 0.95)
# # th.pcc.neg <- quantile(m.pcc, 0.05)
# # 
# # # > th.pcc.neg
# # # 5% 
# # # -0.363049 
# # # > th.pcc.pos
# # # 95% 
# # # 0.4028822 
# # 
# # # generate high confidence coexpression set
# # m.pcc.binary <- m.pcc
# # m.pcc.binary[m.pcc.binary <= th.pcc.neg] <- -1
# # m.pcc.binary[m.pcc.binary > th.pcc.neg] <- 0
# # m.pcc.binary[m.pcc >= th.pcc.pos] <- 1
# # 
# # m.pcc.binary <- m.pcc.binary[which(rowSums(m.pcc.binary) != 0), which(colSums(m.pcc.binary) != 0)]
# # 
# # 
# # 
# 
# 
# 
# 
# ## Finding enriched regulators for these targets
# 
# v.env <- c("AT1G66480", "AT5G07060", "AT5G07080", "AT3G46090", "AT3G03370", "AT3G03380")
# v.env <- c(v.env, "AT5G07080") # early phosphate deficiency (3h)
# 
# 
# # AT5G07080 # early phosphate deficiency (3h)
# # AT1G66480 (6h)
# 
# v.env.onchip <- intersect(rownames(m.pcc.binary),v.env)
# 
# m.pcc[v.env.onchip,v.env.onchip]
# 
# 
# l.ttests[[1]][v.env.onchip,]
# l.ttests[[2]][v.env.onchip,]
# 
# 
# # generate a differential expression set (co-specificity analysis)
# m.anova.jc <- jaccard(m.anova)
# rownames(m.anova.jc)  <- colnames(m.anova.jc) <- rownames(m.anova)
# diag(m.anova.jc) <- 0
# 
# v.gns.overlap <- intersect(rownames(m.anova.jc), rownames(m.pcc.binary))
# 
# m.overlap <- m.pcc.binary[v.gns.overlap,v.gns.overlap] * as.matrix(m.anova.jc[v.gns.overlap,v.gns.overlap])
# 
# 
# 
# tfs.reg <- names(which(rowSums(m.dna[,intersect(colnames(m.dna),v.env.onchip)]) > 0))
# 
# m.dna[tfs.reg,intersect(colnames(m.dna),v.env.onchip)]
# 
# 
# m.pcc.binary
# 
# m.anova[v.env,]
# 
# # diff exp regulator
# 
# colnames(df.geneExp.set) <- colnames(m.anova.set)
# 
# # df.geneExp.noLog <- exp(df.geneExp.set)
# 
# m.anova.set <- m.anova[,c(49,50,51,87,88,89,109,110,111),drop = FALSE]
# m.anova.set <- m.anova.set[,c(1,4,7,2,5,8,3,6,9)]
# 
# m.anova.set <- m.anova.set[which(rowSums(m.anova.set) > 0),]
# v.gns.diffexp <- rownames(m.anova.set)
# 
# df.geneExp.noLog.set <- df.geneExp.noLog
# 
# df.geneExp.set <- df.geneExp.set[v.gns.diffexp,]
# df.geneExp.set <- df.geneExp.set[which(!is.na(rowSums(df.geneExp.set))),]
# 
# 
# table(colSums(m.anova.set))
# sum(m.anova.set)
# 
# df.geneExp.noLog["AT2G40080",]                                     
# 
# #m.anova.set <- m.anova[,c(4,5,6,49,50,51,87,88,89,109,110,111),drop = FALSE]
# 
# m <- m.anova[intersect(tfs, rownames(m.anova)),c(4,5,6,49,50,51,87,88,89,109,110,111),drop = FALSE]
# idx.regs <- which(rowSums(m) > 0)
# 
# exps <- c("-P+Fe", "+P-Fe", "-P-Fe")
# expSets <- vector(mode = "list", length = 3)
# expSets[[1]] <- c(1,2,3)
# expSets[[2]] <- c(4,5,6)
# expSets[[3]] <- c(7,8,9)
# 
# for(i in 1:3){
#   
#   print(exps[i])
#   
#   idx.gns <- which(rowSums(m.anova.set[,expSets[[i]]]) > 0)
#   
#   print(length(idx.gns)) # number of diff exp genes
#   
#   # highlight 
#   m.anova.set.i <- m.anova.set[,expSets[[i]]] # 
#   
#   df.geneExp.set.i <- df.geneExp.set[,expSets[[i]]]#[!is.na(df.geneExp.set[,expSets[[i]]])]
#   
#   # cluster the profiles - up->down, down-> up, 
#   #  m.tmp <- kmeans(df.geneExp.set[,expSets[[i]]], 4)
#   
#   # expression fold change profiles
#   plot(as.numeric(df.geneExp.set[idx.gns[1],expSets[[i]]]), type = "l", ylim = c(min(df.geneExp.set.i),max(df.geneExp.set.i)))
#   
#   for(j in 2:length(idx.gns)){
#     lines(as.numeric(df.geneExp.set[idx.gns[j],expSets[[i]]]))
#     #df.geneExp.noLog.set[idx.gns[j],expSets[[i]]]  
#   }
#   
#   # filter for specific terms
#   terms <- c("root", "DNA binding", "phosphate", "iron")
#   
#   go.set <- subset(df.geneontology, df.geneontology$AGI %in% names(idx.gns))
#   
#   # general genes with TF annotation
#   tf.set <- subset(go.set, go.set$GO.Molecular.Function.Term %in% c("DNA binding")) # additional binding support ? / coexpression support?
#   
#   # binding support analysis
#   idx.tfs <- which(names(idx.gns) %in% tfs)
#   tf.dna <- names(idx.gns)[idx.tfs]
#   
#   
#   if(length(tf.dna) > 0){
#     
#     for(r in 1:lenth(tf.dna))
#       
#       print(tf.dna[r])
#     
#     plot(as.numeric(df.geneExp.set.i[tf.dna[r],]), type = "l", ylim = c(-3.5,3.5))
#     
#     tgs.r <- names(which(m.dna[tf.dna,] == 1))
#     tgs.r <- intersect(rownames(df.geneExp.set.i), tgs.r)
#     
#     #       
#     #       tf.dna, tgs.r
#     #       
#     m.anova.set[c(tf.dna,tgs.r),]
#     #       
#     #       df.geneExp.set.i[tgs.r,]
#     #       
#     
#     
#     lines(as.numeric(df.geneExp.set.i[tgs.r[1],]), col = "red")
#     lines(as.numeric(df.geneExp.set.i[tgs.r[2],]), col = "blue")
#     
#     
#     
#     #       
#     #       
#     #       if(length(tgs.r) > 0){
#     #         for(j in 2:length(idx.gns)){
#     #           lines(as.numeric(df.geneExp.set[idx.gns[j],expSets[[i]]]), col = "")
#     #           #df.geneExp.noLog.set[idx.gns[j],expSets[[i]]]  
#     #         }
#     #         
#     #         
#     #       }
#     #       
#     
#     
#   }
#   
#   
#   
#   idx.ifs <- which(tf.set$AGI %in% rownames(m.dna))
#   
#   test <- subset(df.geneontology, grepl(terms, df.geneontology$GO.Biological.Process.Term))
#   
#   df.geneExp.set["AT2G24850",]
#   
#   plot(as.numeric(df.geneExp.set["AT2G24850",expSets[[i]]]), type = "l", ylim = c(-3.5,3.5))
#   
#   
# }
# 
# #m[idx.regs,]
# 
# for(i in 1:length(idx.regs)){
#   
#   tf <- names(idx.regs)[i]
#   
#   print(m[tf,])
#   
#   tg.i <- names(which(m.dna[tf,] == 1))
#   
#   #m[tf,]
#   
#   #(2/44) / (78 / 28000) =>  16.31702
#   
#   #   > print(colSums(m.i))
#   #   Tuk_-P:+Fe:3hrs-+P:+Fe:3hrs Tuk_+P:-Fe:3hrs-+P:+Fe:3hrs Tuk_-P:-Fe:3hrs-+P:+Fe:3hrs Tuk_-P:+Fe:6hrs-+P:+Fe:6hrs Tuk_+P:-Fe:6hrs-+P:+Fe:6hrs Tuk_-P:-Fe:6hrs-+P:+Fe:6hrs Tuk_-P:+Fe:9hrs-+P:+Fe:9hrs Tuk_+P:-Fe:9hrs-+P:+Fe:9hrs Tuk_-P:-Fe:9hrs-+P:+Fe:9hrs 
#   #   0                           0                           0                           2                           0                           0                           0                           0                           0 
#   #   > colSums(m.anova[,c(49,50,51,87,88,89,109,110,111),drop = FALSE])
#   #   Tuk_-P:+Fe:3hrs-+P:+Fe:3hrs Tuk_+P:-Fe:3hrs-+P:+Fe:3hrs Tuk_-P:-Fe:3hrs-+P:+Fe:3hrs Tuk_-P:+Fe:6hrs-+P:+Fe:6hrs Tuk_+P:-Fe:6hrs-+P:+Fe:6hrs Tuk_-P:-Fe:6hrs-+P:+Fe:6hrs Tuk_-P:+Fe:9hrs-+P:+Fe:9hrs Tuk_+P:-Fe:9hrs-+P:+Fe:9hrs Tuk_-P:-Fe:9hrs-+P:+Fe:9hrs 
#   #   33                          20                          68                          44                          25                          23                          22                          18                          36 
#   #   
#   colSums(m.anova[,c(49,50,51,87,88,89,109,110,111),drop = FALSE])
#   
#   m.i <- m.anova[intersect(tg.i, rownames(m.anova)),c(49,50,51,87,88,89,109,110,111),drop = FALSE]
#   
#   print(colSums(m.i))
#   print(which(rowSums(m.i) > 0))
#   #colSums(m.i)
#   # enrichment
#   
#   if(tf %in% rownames(m.pcc.binary)){
#     
#     idx.pcc <- which(m.pcc.binary[tf, intersect(tg.i, colnames(m.pcc.binary))] != 0)
#     tg.pcc <- m.pcc[tf, intersect(tg.i, colnames(m.pcc.binary))][idx.pcc]
#     print(tg.pcc)  
#   }
#   
#   AT3G11800
#   
#   print("")
#   print("")
#   print("")
#   
# }
# 
# 
# 
# 
# 
# write.csv(m[idx.regs,], "RegulatorAnova1.csv")
# 
# #idx.regs <- which(rowSums(m.anova[intersect(tfs, rownames(m.anova)),c(4,5,6,7),drop = FALSE]) > 0)
# 
# m <- m.anova[intersect(tfs, rownames(m.anova)),c(4,5,6,7),drop = FALSE]
# idx.regs <- which(rowSums(m) > 0)
# 
# for(i in 1:length(idx.regs)){
#   
#   tf <- names(idx.regs)[i]
#   
#   
#   if(m[tf,2] == 1){ 
#     
#     print(i)
#     print(m[tf,])
#     
#     tg.i <- names(which(m.dna[tf,] == 1))
#     
#     #colSums(m.anova[,c(4,5,6,7),drop = FALSE])
#     
#     
#     m.i <- m.anova[intersect(tg.i, rownames(m.anova)),c(4,5,6,7),drop = FALSE]
#     colSums(m.i)
#     
#     tg.i <- names(which(m.i[,2] == 1))
#     
#     print(colSums(m.i))
#     print(which(rowSums(m.i) > 0))
#     #colSums(m.i)
#     # enrichment
#     
#     if(tf %in% rownames(m.pcc.binary)){
#       
#       idx.pcc <- which(m.pcc.binary[tf, intersect(tg.i, colnames(m.pcc.binary))] != 0)
#       tg.pcc <- m.pcc[tf, intersect(tg.i, colnames(m.pcc.binary))][idx.pcc]
#       print(tg.pcc)  
#     }
#   }
#   
#   print("")
#   print("")
#   print("")
#   
# }
# 
# 
# 
# 
# 
# 
# m[idx.regs,]
# 
# write.csv(m[idx.regs,], "RegulatorAnova2.csv")
# 
# table(rowSums(m.anova[intersect(tfs, rownames(m.anova)),c(49,50,51,87,88,89,109,110,111),drop = FALSE]))
# 
# ###
# 
# m.dna
# 
# 
# 
# m.pcc[intersect(rownames(m.pcc), tfs.reg),intersect(colnames(m.dna),v.env.onchip)]
# 
# ###
# tfs
# tgs
# m.dna
# l.ttests
# 
# v.gns.pcc <- rownames(m.pcc.binary)
# 
# # analyze individual sets #4 out of 4
# for(i in 1:4){
#   
#   for(j in 1:3){
#     
#     m.tmp <- l.ttests[[1]][,1:4]
#     tfs.tmp <- intersect(tfs, rownames(m.tmp))
#     #m.tmp[,]
#     
#     m.tmp.reg <- m.tmp[tfs.tmp,]
#     v.tmp.reg <- apply(m.tmp.reg, 1, function(m) length(which(m < 0.05))) # accumulation rank
#     
#     # rank based on differential expression
#     table(v.tmp.reg)
#     
#     idx.tfSet <- (which(v.tmp.reg < 0.05)) 
#     
#     ## a) TFs expressed
#     
#     rank.tf <- m.tmp.reg
#     
#     tfs.tmp[idx.tfSet]
#     
#     tfs.sel <- names(v.tmp.reg[v.tmp.reg > 0])
#     
#     v.tmp.tg <- apply(m.tmp, 1, function(m) length(which(m < 0.05))) # accumulation rank
#     table(v.tmp.tg)
#     
#     tgs.sel <- names(v.tmp.tg[v.tmp.tg > 0])
#     
#     tfs.sel <- intersect(tfs.sel,v.gns.pcc)
#     tgs.sel <- intersect(tgs.sel,v.gns.pcc)
#     
#     table(m.pcc.binary[tfs.sel,tgs.sel])
#     
#     idx.tgSet <- (which(v.tmp.tg < 0.001)) 
#     
#     ## b) targets expressed (at or after TF)
#     
#     
#     
#     ## c) dna binding
#     
#     
#     ## d) pcc binary confirmation
#     
#     
#     # shifted or same column subset matrix #
#     tgs.sel.dna <- intersect(tgs.sel, colnames(m.dna))
#     
#     
#     # ranked based on number of targets or enrichment
#     table(rowSums(m.dna[tfs.sel, tgs.sel.dna]))
#     
#     
#     
#     m.tmp
#     
#     m.dna[tfs.tmp[idx.tfSet],]
#     
#     table(rowSums(m.dna[tfs.tmp[idx.tfSet],]))
#     
#     # up, down regulation 
#     
#     
#     
#     
#     
#     l.ttests[[2]]
#     
#     
#     
#     
#   }
#   
#   
#   
# }
# 
# ## regulatory analysis to establish 
# v.tfs.anova <- intersect(tfs, rownames(m.anova))
# 
# m.overlap.set <- m.overlap
# idx <- which(rowSums(m.overlap.set) != 0)
# m.overlap.set <- m.overlap.set[idx,idx]
# 
# if(FALSE)
#   hist(m.overlap.set[m.overlap.set!=0], breaks = 100)
# 
# v.tfs.set <- intersect(v.tfs, rownames(m.overlap.set))
# 
# 
# m.overlap.regulation <- m.overlap.set[v.tfs.set,]
# rowSums(m.overlap.regulation)
# 
# df.regulation <- data.frame(TF = character(), targets = character(), TF_in_SNP = character(), targets_in_SNP = character())
# 
# for(i in 1:length(v.tfs.set)){
#   
#   print(v.tfs.set[i])
#   targets <- names(which(m.overlap.set[v.tfs.set[i],] != 0 ))
#   
#   #TF_in_SNP <- character()
#   #targets_in_SNP <- character()
#   
#   for(j in 1:length(l.gwas)){
#     
#     
#     TF_in_SNP.tmp <- intersect(v.tfs.set[i], l.gwas[[j]]$SNP.relative.position)
#     targets_in_SNP.tmp <- intersect(targets, l.gwas[[j]]$SNP.relative.position)
#     
#     if(length(TF_in_SNP.tmp) > 0)
#       print("yes")
#     
#     print(targets_in_SNP.tmp)
#     
#     #TF_in_SNP <- paste(TF_in_SNP, paste(TF_in_SNP.tmp, names(l.gwas)[j], collapse = " - "),  collapse = "| ")
#     #targets_in_SNP <- paste(targets_in_SNP, paste(paste(targets_in_SNP.tmp, collapse = "; "), names(l.gwas)[j], collapse = " - "),  collapse = "| ")
#     
#   }
#   
#   print("")
#   #df.geneExp[v.tfs.set[i],]
#   #df.regulation <- rbind(df.regulation, data.frame(TF =  v.tfs.set[i], targets = paste(names(which(m.overlap.set[v.tfs.set[i],] != 0)), collapse = "; ")), TF_in_SNP = TF_in_SNP, targets_in_SNP = targets_in_SNP)
# }
# 
# write.csv(df.regulation, "df.regulation_novel.csv", row.names = FALSE)
# 
# subset(df.geneontology, df.geneontology$AGI %in% v.tfs.set)$GO.Biological.Process.Term
# 
# ###
# 
# # library(igraph)
# # library(parmigene)
# 
# # normalize coexpression around zero mean 
# 
# # test: alternatives for sparse graph generation #
# # m.clr <- parmigene::clr(m.pcc)
# # m.network <- m.pcc # m.clr
# # th <- 0.9 #quantile(m.pcc, 0.99)
# # m.network[m.network <= th] <- 0
# # m.network[m.network > 0] <- 1
# # 
# # g <- graph_from_adjacency_matrix(m.network)
# # wc <- infomap.community(g)
# # 
# # # v.modules <- unique(membership(wc))
# # tb.modules <- table(membership(wc))
# # v.modules <- names(tb.modules[tb.modules > 1])
# 
# 
# p.modules <- rep(1, length(v.modules))
# names(p.modules) <- v.modules
# 
# # p.dna.modules <- rep(1, length(v.modules))
# # names(p.dna.modules) <- v.modules
# 
# # m = 1 - non assigments
# for(m in 2:length(v.modules)){
#   
#   module <- names(which(v.module_membership == v.modules[m]))
#   # module <- names(which(membership(wc) == v.modules[m]))
#   v.gns.anova <- intersect(rownames(m.anova), module)
#   
#   if(length(v.gns.anova) > 0){
#     # disregarding splitting of test combinations
#     hitInSample = length(which(rowSums(m.anova[v.gns.anova,,drop=FALSE]) > 0))
#     
#     hitInPop = n.DE 
#     failInPop = dim(m.anova)[1] - hitInPop 
#     sampleSize = length(v.gns.anova)
#     
#     if(hitInSample > 1){
#       p.modules[m] <- phyper(hitInSample - 1, hitInPop, failInPop, sampleSize, lower.tail= FALSE);  
#     }
#     
#   }
#   
# }
# 
# # significantly enriched (differentially genes) 
# v.sigModules <- names(which(p.modules <= 0.05))
# 
# df.module <- data.frame(id = character(), gns = character(), diffExpGns = character(), GO = character())
# 
# for(m in 1:length(v.sigModules)){
#   
#   print(paste("Module:", m))
#   
#   module <- names(which(v.module_membership == v.sigModules[m]))
#   #module <- names(which(membership(wc) == v.sigModules[m]))
#   v.gns.anova <- intersect(rownames(m.anova), module)
#   
#   # m.anova[v.gns.anova,,drop=FALSE]
#   
#   #print(rowSums(m.anova[v.gns.anova,,drop=FALSE]))
#   
#   v.gns.diff <- names(which(rowSums(m.anova[v.gns.anova,,drop=FALSE]) == 1))
#   m.anova[v.gns.diff,]
#   
#   # correlation scores
#   m.pcc[v.gns.diff,v.gns.diff]
#   
#   ##
#   conds <- c("Ctrl", "+P+Fe(3h)", "+P+Fe(6h)", "+P+Fe(9h)", "-P+Fe(3h)", "-P+Fe(6h)", "-P+Fe(9h)", "+P-Fe(3h)", "+P-Fe(6h)", "+P-Fe(9h)", "-P-Fe(3h)", "-P-Fe(6h)", "-P-Fe(9h)")
#   
#   df.geneExp[v.gns.diff,]
#   
#   # quad plot # 
#   
#   
#   # additional information #
#   v.gwas <- l.gwas[[1]]
#   v.gwas <- unique(v.gwas$SNP.relative.position)
#   
#   v.gwas.m <- intersect(v.gwas, v.gns.anova)
#   
#   if(length(v.gwas.m) > 0){
#     # disregarding splitting of test combinations
#     hitInSample = length(v.gwas.m)
#     hitInPop = length(v.gwas) 
#     failInPop = dim(m.anova)[1] - hitInPop 
#     sampleSize = length(module)
#     
#     if(hitInSample > 1){
#       p.modules[m] <- phyper(hitInSample - 1, hitInPop, failInPop, sampleSize, lower.tail= FALSE);  
#     }
#     
#   }
#   
#   
#   # GO 
#   v.go <- subset(df.geneontology, df.geneontology$AGI %in% v.gns.diff)$GO.Biological.Process.Term
#   tb.go <- table(unlist(sapply(v.go, function(m) strsplit(m, "; "))))
#   
#   th.min <- 2
#   tb.go <- tb.go[tb.go >= th.min]
#   
#   print("Gene ontology (biological process) - differential expressed genes")
#   print(tb.go)
#   
#   ## 
#   
#   v.go <- subset(df.geneontology, df.geneontology$AGI %in% v.gns.anova[!v.gns.anova %in% v.gns.diff])$GO.Biological.Process.Term
#   tb.go <- table(unlist(sapply(v.go, function(m) strsplit(m, "; "))))
#   
#   th.min <- 2
#   tb.go <- tb.go[tb.go >= th.min]
#   
#   #print("Gene ontology (biological process) - remaining genes")
#   #print(tb.go)
#   
#   # dna blueprint 
#   # df.dna.m <- subset(df.dna, df.dna$Target.gene %in% v.gns.diff)
#   # table(df.dna.m$Putative.TF)
#   # 
#   # df.dna.m <- subset(df.dna, df.dna$Target.gene %in% v.gns.anova[!v.gns.anova %in% v.gns.diff])
#   # table(df.dna.m$Putative.TF)
#   # 
#   
#   # fold change of expression (mean gene expression )
#   
#   # nrow(subset(df.dna, df.dna$Putative.TF == "AT2G45660"))
#   # 
#   # v.go <- subset(df.geneontology, df.geneontology$AGI %in% "AT1G13260")$GO.Biological.Process.Term
#   # tb.go <- table(unlist(sapply(v.go, function(m) strsplit(m, "; "))))
#   # 
#   
#   # regulation plot - split into 4 figures per experiment and time series
#   # df.exp.m <- subset(df.exp, df.exp$AGI == "AT1G13260")
#   # plot(as.numeric(df.exp.m[,1:12]), type = "l", ylim = c(0,10))
#   # 
#   # df.exp.m <- subset(df.exp, df.exp$AGI %in% v.gns.diff)
#   # lines(as.numeric(df.exp.m[1,1:12]), col = "red")
#   # lines(as.numeric(df.exp.m[2,1:12]), col = "blue")
#   # 
#   
#   tmp <- data.frame(id = v.sigModules[m], gns = paste(v.gns.anova, collapse = ","), diffExpGns = paste(v.gns.diff, collapse = ","), GO = paste(names(tb.go), collapse = ","))
#   
#   df.module <- rbind(df.module, tmp)
#   
#   
#   
#   print("") 
# }
# 
# write.csv(df.module, paste("moduleAnalysis.csv"))
# 
# #source("https://bioconductor.org/biocLite.R")
# #biocLite("genefilter")
# 
# 
# # expression analysis #\
# df.exp <- read.csv("geneExpData.csv", header = TRUE)
# df.exp["AGI"] <- v.map[as.character(df.exp$ClusterID)]
# df.exp <- df.exp[,-1]
# 
# 
# v.control <- c(1,14,27) # control
# 
# l.timeseries <- vector(mode = "list", length = 3)
# 
# # 3 hours
# l.timeseries[[1]] <- vector(mode = "list", length = 4)
# l.timeseries[[1]][[1]] <- c(2,15,28) # +P+Fe
# l.timeseries[[1]][[2]] <- c(3,16,29) # -P+Fe
# l.timeseries[[1]][[3]] <- c(4,17,30) # +P-Fe
# l.timeseries[[1]][[4]] <- c(5,18,31) # -P-Fe
# 
# # 6 hours
# l.timeseries[[2]] <- vector(mode = "list", length = 4)
# l.timeseries[[2]][[1]] <- c(6,19,32) # +P+Fe
# l.timeseries[[2]][[2]] <- c(7,20,33) # -P+Fe
# l.timeseries[[2]][[3]] <- c(8,21,34) # +P-Fe
# l.timeseries[[2]][[4]] <- c(9,22,35) # -P-Fe
# 
# # 9 hours
# l.timeseries[[3]] <- vector(mode = "list", length = 4)
# l.timeseries[[3]][[1]] <- c(10,23,36) # +P+Fe
# l.timeseries[[3]][[2]] <- c(11,24,37) # -P+Fe
# l.timeseries[[3]][[3]] <- c(12,25,38) # +P-Fe
# l.timeseries[[3]][[4]] <- c(13,26,39) # -P-Fe
# 
# 
# l.DE.conditions <- vector(mode = "list", length = 4)
# l.DE.conditions[[1]] <- matrix(0, nrow = length(df.exp$AGI), ncol = 3, dimnames = list(df.exp$AGI, c(3,6,9)))
# l.DE.conditions[[2]] <- matrix(0, nrow = length(df.exp$AGI), ncol = 3, dimnames = list(df.exp$AGI, c(3,6,9)))
# l.DE.conditions[[3]] <- matrix(0, nrow = length(df.exp$AGI), ncol = 3, dimnames = list(df.exp$AGI, c(3,6,9)))
# l.DE.conditions[[4]] <- matrix(0, nrow = length(df.exp$AGI), ncol = 3, dimnames = list(df.exp$AGI, c(3,6,9)))
# names(l.DE.conditions) <- c("+P+Fe", "-P+Fe", "+P-Fe", "-P-Fe")
# 
# for(i in 1:length(l.DE.conditions)){
#   
#   for(j in 1:3){
#     
#     v.control <- l.timeseries[[j]][[1]]
#     m <- as.matrix(df.exp[,c(v.control,l.timeseries[[j]][[i]])])
#     f <- factor(c(rep("Control", 3), rep("Treatment", 3)))
#     ttests <- rowttests(m, f)
#     
#     # ttests$p.value <- p.adjust(ttests$p.value, "fdr")
#     idx <- which(ttests$p.value < 0.01)
#     l.DE.conditions[[i]][idx, j] <- ttests$dm[idx]
#     l.DE.conditions[[i]][, j] <- ifelse(l.DE.conditions[[i]][, j] > 0, -1, ifelse(l.DE.conditions[[i]][, j] < 0, 1, 0))
#     
#   }
#   
# }
# 
# lapply(l.DE.conditions, table)
# 
# ###
# 
# 
# jaccard <- function(m) {
#   ## common values:
#   A = tcrossprod(m)
#   ## indexes for non-zero common values
#   im = which(A > 0, arr.ind=TRUE)
#   ## counts for each row
#   b = rowSums(m)
#   ## only non-zero values of common
#   Aim = A[im]
#   ## Jacard formula: #common / (#i + #j - #common)
#   J = sparseMatrix(
#     i = im[,1],
#     j = im[,2],
#     x = Aim / (b[im[,1]] + b[im[,2]] - Aim),
#     dims = dim(A)
#   ) 
#   return( J )
# }
# 
# 
# # m.test.global <- c()
# # for(i in 1:12){
# #   
# #   #   
# #   #   Value <- c(1,2,5,3,2,1,1,3,2,1,4,3,6,5,2,6,1,6,5,4,9,6,7,7,5,1,8,9,6,5)
# #   #   Group <- factor(c(rep(1,10),rep(2,10),rep(3,10)))
# #   #   
# #   #   data <- data.frame(categ, resp)
# #   #   
# #   #   kruskal.test(Value ~ Group, data=data)
# #   #   pairwise.wilcox.test(resp, categ, p.adj="holm", exact=F)
# #   #     
# #   m.test.global <- cbind(m.test.global, apply(df.geneExp, 1, function(m) {t.test(as.numeric(m[c(1+i,14+i,27+i)]), as.numeric(m[c(1,14,27)]))$p.value}))
# # 
# #   
# # #   m.test <- cbind(m.test, apply(df.geneExp, 1, function(m) {
# # #                                                           resp<-c(as.numeric(m[c(1+i,14+i,27+i)]), as.numeric(m[c(1,14,27)]))
# # #                                                           categ<-as.factor(rep(c("A","B"),times=1,each=3))
# # #                                                           ifelse(kruskalmc(resp, categ, probs = 0.05, cont="two-tailed")$dif.com$difference, 1, 0)
# # #                                                         }))
# #                                                           
# #   #    m <- df.geneExp[1,]
# #   #     
# #   #     wilcox.test(as.numeric(m[c(1,14,27)]), as.numeric(m[c(1+i,14+i,27+i)]))$p.value))
# #   #   
# #   #   
# #   
# #   #m.test <- cbind(m.test, apply(df.geneExp, 1, function(m) wilcox.test(as.numeric(m[c(1,14,27)]), as.numeric(m[c(1+i,14+i,27+i)]))$p.value))
# #   #m.test <- cbind(m.test, apply(df.geneExp, 1, function(m) t.test((as.numeric(m[c(1+i,14+i,27+i)])), (as.numeric(m[c(1,14,27)])) )$p.value ))
# #   # m.test <- cbind(m.test, apply(df.geneExp, 1, function(m) (mean(as.numeric(m[c(1+i,14+i,27+i)])) - mean(as.numeric(m[c(1,14,27)])) )))
# # }
# 
# 
# # times <- c(0, 3, 6, 9)
# # 
# # df.exp <- processRawData(m.expression.diffExp.set, times, experiments=rep(1:2, each=7), is.logged = TRUE)#, do.normalisation = ifelse(is.logged, TRUE, FALSE))
# # 
# 
# # The probe identifier for TF 'twi'
# # tfs <- (intersect(v.gns.diffexp, rownames(m.dna)))
# # tfs <- intersect(tfs, rownames(m.expression.diffExp.set))
# # 
# # # The probe identifier for the target gene
# # targets <- (intersect(v.gns.diffexp, colnames(m.dna)))
# # targets <- intersect(targets, rownames(m.expression.diffExp.set))
# # 
# # 
# # exps <- c("-P+Fe", "+P-Fe", "-P-Fe")
# # expSets <- vector(mode = "list", length = 3)
# # expSets[[1]] <- c(1,2,3)
# # expSets[[2]] <- c(4,5,6)
# # expSets[[3]] <- c(7,8,9)
# # 
# # for(i in 1:3){
# #   
# # 
# # 
# # m.loglik <- matrix(0, length(tfs), length(targets), dimnames = list(tfs, targets))
# # 
# # for(i in 1:length(tfs)){
# #   
# #   #v.loglik <- numeric(length(targets))
# #   #names(v.loglik) <- targets
# #   
# #   for(j in 1:length(targets)){ # parallelize
# #     
# #     # Learn the model using only one of the 3 repeats in the data
# #     model <- GPLearn(df.exp[,1:7],
# #                      TF=tfs[i], targets = targets[j],
# #                      useGpdisim=TRUE, quiet=TRUE)
# #     
# #     #m.expression.diffExp.set[targets[1], ]
# #     
# #     # Display the model parameters
# #     show(model)
# #     
# #     GPPlot(model)
# #     
# #     m.loglik[i,j] <- - model@model$llscore
# #     
# #   }
# #   
# #   
# # }
# 
# #df.geneExp <- (df.geneExp[1:13] + df.geneExp[14:26] + df.geneExp[27:39])/3
# 
# 
# # library("ggplot2")
# # 
# # N <- 30
# # id <- as.character(1:N) # create ids
# # sexes = c("male", "female")
# # sex <- sample(sexes, size = N/2, replace = TRUE) # create a sample of sex
# # diseases = c("low", "med", "high")
# # disease <- rep(diseases, each = N/3) # disease severity 
# # times = c("Pre", "0", "30", "60")
# # time <- rep(times, times = N) # times measured 
# # t <- 0:3
# # ntimes = length(t)
# # y1 <- c(replicate(N/2, rnorm(ntimes, mean = 10+2*t)), 
# #         replicate(N/2, rnorm(ntimes, mean = 10+4*t)))
# # y2 <- c(replicate(N/2, rnorm(ntimes, mean = 10-2*t)), 
# #         replicate(N/2, rnorm(ntimes, mean = 10-4*t)))
# # y3 <- c(replicate(N/2, rnorm(ntimes, mean = 10+t^2)), 
# #         replicate(N/2, rnorm(ntimes, mean = 10-t^2)))
# # 
# # data <- data.frame(id=rep(id, each=ntimes), sex=rep(sex, each=ntimes), 
# #                    severity=rep(disease, each=ntimes), time=time, 
# #                    Y1=c(y1), Y2=c(y2), Y3=c(y3)) # create data.frame
# # #### factor the variables so in correct order
# # data$sex = factor(data$sex, levels = sexes)
# # data$time = factor(data$time, levels = times)
# # data$severity = factor(data$severity, levels = diseases)
# # head(data)
# 
# 
# # compute the standard deviations for the fold changes
# # 
# # df.geneExp <- df.geneExp[,c(1,2,6,10,3,7,11,4,8,12,5,9,13)]
# # 
# # df.geneExp.set <- df.geneExp[,-1]
# # df.geneExp.set <- df.geneExp.set[,c(1,4,7,10,2,5,8,11,3,6,9,12)]
# # df.geneExp.set[,4:6] <- df.geneExp.set[,4:6] - df.geneExp.set[,1]
# # df.geneExp.set[,7:9] <- df.geneExp.set[,7:9] - df.geneExp.set[,2]
# # df.geneExp.set[,10:12] <- df.geneExp.set[,10:12] - df.geneExp.set[,3]
# # 
# # df.geneExp.set <- df.geneExp.set[,4:12]
# # 
