
#setwd("/shared/Labs/Rhee Lab/Everyone/Michael/Hatem-FePNetwork/")
# setwd("/Users/michaelbanf/Documents/postdoctoral_work/Projects/Hatem")

# setwd("Desktop/At_Root_GeneExp/")
rm(list=ls()) # clear workspace 

# CALC
# setwd("/home/mbanf/Documents/Computational_Biology/Projects/GRACE_release/")


# source("http://bioconductor.org/biocLite.R") 
# biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore")) 
# install.packages("WGCNA")

library(WGCNA)
library(genefilter)
library(parmigene)


# two gene expression sets - root related
if(FALSE){
  
    mat.GE.Ath <- read.table("GSE69995_re-analyzed_data_matrix.txt", row.names = 1, header = TRUE, sep = "\t", quote = "")
    mat.GE.Ath <- as.matrix(mat.GE.Ath)
    map.GE.MicroArray <- read.csv("tpj13175-sup-0003-TableS2.csv", stringsAsFactors = FALSE)
    
    v.conditions <- paste(map.GE.MicroArray$Series.Title, map.GE.MicroArray$Sample.Title, sep = ";")
    colnames(mat.GE.Ath) <- v.conditions
    
    m.rootSet <- mat.GE.Ath[,which(map.GE.MicroArray$Tissue == "root")]
    v.gns <- rownames(m.rootSet)
  
    m.pcc <- cor(t(m.rootSet))
    saveRDS(m.pcc, "m.pcc_big.rds")
    
    m.mi <- parmigene::knnmi.all(m.rootSet) # capturing nonlinear expression relationships
    saveRDS(m.mi, "m.mi_big.rds")
    
}else{
  
  # smaller root set
  df.rootSet <- read.table("Root_AtExp_v3.txt", header = TRUE, sep = "\t")
  v.gns <- toupper(df.rootSet$X)
  df.rootSet <- df.rootSet[,-1]
  
  m.pcc <- cor(t(df.rootSet))
  saveRDS(m.pcc, "m.pcc_small.rds")
  
  m.pcc <- readRDS("m.pcc_small.rds")
  
  m.mi <- parmigene::knnmi.all(df.rootSet)
  saveRDS(m.mi, "m.mi_small.rds")
  
}
  
rownames(m.pcc) <- colnames(m.pcc) <- v.gns

erer

l.association <- vector(mode = "list", length = 4)
l.association[[1]] <- 
l.association[[2]] <- 
l.association[[3]] <- 
l.association[[4]] <- 


# if(FALSE){
#   beta <- 6 # paper to cite?
#   ADJ1 <- abs(m.pcc)^beta
#   dissADJ = 1 - ADJ1
#   hierADJ = hclust(as.dist(dissADJ), method="average" )
#   
#   v.module_membership = cutreeDynamic(hierADJ, method="tree", minClusterSize = 10) 
#   names(v.module_membership) <- hierADJ$labels
#   
#   v.modules  <- unique(v.module_membership)
# }

#### Differential expression data

df.diffexp <- read.csv("diffexp_highConfidence.csv", stringsAsFactors = FALSE)
df.geneontology <- read.csv("Table annotation_INITIAL.csv", stringsAsFactors = FALSE) #read.table("Table annotation_INITIAL.txt", header = TRUE, sep = "\t", quote = "")

v.map <- df.geneontology$AGI
names(v.map) <- df.geneontology$Transcript.Cluster.ID

df.diffexp["AGI"] <- v.map[as.character(df.diffexp$Cluster)]

### raw gene expression data ### 

df.geneExp <- read.csv("geneExpData.csv", header = TRUE, stringsAsFactors = FALSE)
df.geneExp["AGI"] <- v.map[as.character(df.geneExp$Cluster)]
df.geneExp <- df.geneExp[,-1]

v.gnSets <- names(which(table(df.geneExp$AGI) == 1))
df.geneExp <- subset(df.geneExp, df.geneExp$AGI %in% v.gnSets)
v.gnSets <- df.geneExp$AGI



# install.packages("pgirmess")
# library(pgirmess)
# 
# kruskalmc function in pgirmess

# test <-  data.frame(exp = c(as.numeric(df.geneExp[1,c(1,14,27)]), as.numeric(df.geneExp[1,c(2,15,28)]), as.numeric(df.geneExp[1,c(3,16,29)])),
#            set = factor(c(rep("A",3), rep("B",3), rep("C",3))),
#            number = factor(rep(c(1,2,3),3)))

# https://www.qbaseplus.com/knowledge/blog/seven-tips-bio-statistical-analysis-gene-expression-data
# are these values log normalized?

# two different control condition runs

# Test 1 - comparison to initial control 

m.test.global <- c()

# kruskal better evaluation for correcting 
m.test.global <- cbind(m.test.global, apply(df.geneExp, 1, function(m) {t.test(as.numeric(m[c(2,15,28)]), as.numeric(m[c(1,14,27)]))$p.value}))
m.test.global <- cbind(m.test.global, apply(df.geneExp, 1, function(m) {t.test(as.numeric(m[c(6,19,32)]), as.numeric(m[c(1,14,27)]))$p.value}))
m.test.global <- cbind(m.test.global, apply(df.geneExp, 1, function(m) {t.test(as.numeric(m[c(10,23,36)]), as.numeric(m[c(1,14,27)]))$p.value}))

for(i in 1:3){
  
  m.test.global <- cbind(m.test.global, apply(df.geneExp, 1, function(m) {t.test(as.numeric(m[c(2+i,15+i,28+i)]), as.numeric(m[c(1,14,27)]))$p.value}))
  m.test.global <- cbind(m.test.global, apply(df.geneExp, 1, function(m) {t.test(as.numeric(m[c(6+i,19+i,32+i)]), as.numeric(m[c(1,14,27)]))$p.value}))
  m.test.global <- cbind(m.test.global, apply(df.geneExp, 1, function(m) {t.test(as.numeric(m[c(10+i,23+i,36+i)]), as.numeric(m[c(1,14,27)]))$p.value}))
  
}

# Test 2 - comparison to individual control conditions
m.test.individual <- c()
for(i in 1:3){

  m.test.individual <- cbind(m.test.individual, apply(df.geneExp, 1, function(m) {t.test(as.numeric(m[c(2+i,15+i,28+i)]), as.numeric(m[c(2,15,28)]))$p.value}))
  m.test.individual <- cbind(m.test.individual, apply(df.geneExp, 1, function(m) {t.test(as.numeric(m[c(6+i,19+i,32+i)]), as.numeric(m[c(6,19,32)]))$p.value}))
  m.test.individual <- cbind(m.test.individual, apply(df.geneExp, 1, function(m) {t.test(as.numeric(m[c(10+i,23+i,36+i)]), as.numeric(m[c(10,23,36)]))$p.value}))
  
}

l.ttests <- vector(mode = "list", length = 2)
l.ttests[[1]] <- as.matrix(m.test.global)
l.ttests[[2]] <- as.matrix(m.test.individual)
names(l.ttests) <- c("global_control", "separate_controls")

rownames(l.ttests[[1]]) <- v.gnSets
rownames(l.ttests[[2]]) <- v.gnSets

# m.anova <- as.matrix(df.anova)
# rownames(m.anova) <- agi

# length(which(m.test > 2))
# length(which(m.test < 0.0001)) # every 1 in 10000

# notion: no multiple hypothesis correction to retain 

# pvalues <- m.test[,4]
# length(which(pvalues < 0.0001))
# pvalues <- p.adjust(pvalues, "bonferroni")
# length(which(pvalues < 0.05))


# hist(df.geneExp[,c(2,15,28)])

# Test 2 - comparison to individual controls


# Test 2 - comparison to time complete 

df.geneExp <- (df.geneExp[1:13] + df.geneExp[14:26] + df.geneExp[27:39])/3
rownames(df.geneExp) <- v.gnSets

# compute the standard deviations for the fold changes


df.geneExp <- df.geneExp[,c(1,2,6,10,3,7,11,4,8,12,5,9,13)]

df.geneExp.set <- df.geneExp[,-1]
df.geneExp.set <- df.geneExp.set[,c(1,4,7,10,2,5,8,11,3,6,9,12)]
df.geneExp.set[,4:6] <- df.geneExp.set[,4:6] - df.geneExp.set[,1]
df.geneExp.set[,7:9] <- df.geneExp.set[,7:9] - df.geneExp.set[,2]
df.geneExp.set[,10:12] <- df.geneExp.set[,10:12] - df.geneExp.set[,3]

df.geneExp.set <- df.geneExp.set[,4:12]


#### GWAS data 
l.gwas <- vector(mode = "list", length = 3)
l.gwas[[1]] <- read.csv("AT-GWAS/fe/Fe.csv", stringsAsFactors = FALSE)
l.gwas[[2]] <- read.csv("AT-GWAS/Pi/Pi.csv", stringsAsFactors = FALSE)
l.gwas[[3]] <- read.csv("AT-GWAS/pfe/pfe.csv", stringsAsFactors = FALSE)

names(l.gwas) <- c("Fe", "Pi", "PFe")


for(j in 1:length(l.gwas)){
  
  l.gwas[[j]]$SNP.relative.position <- gsub("\\..*", "", l.gwas[[j]]$SNP.relative.position)
  
  #TF_in_SNP.tmp <- intersect(v.tfs.set[i], l.gwas[[j]]$SNP.relative.position)
  
}

### DNA binding 

# - replace with updated list #
df.dna <- read.csv("Yu2016-srep25164-s2.csv", stringsAsFactors = FALSE)
df.dna <- df.dna[,c(1,6)]


tfs <- character()
for(i in 1:nrow(df.dna))
  tfs <- c(tfs, unlist(strsplit(df.dna$Putative.TF[i], ", ")))

tfs <- unique(tfs)
tgs <- unique(df.dna$Target.gene)

m.dna <- matrix(0, nrow = length(tfs), ncol = length(tgs), dimnames = list(tfs, tgs))

for(i in 1:nrow(df.dna)){
  
  tfs.i <- unlist(strsplit(df.dna$Putative.TF[i], ", "))
  tg.i <- df.dna$Target.gene[i]

  m.dna[tfs.i, tg.i] <- 1
}


library(reshape2)
setNames(melt(m1), c('rows', 'vars', 'values'))

m.dna <- m.dna[which(rowSums(m.dna) < 2000),]

tfs <- rownames(m.dna)
tgs <- colnames(m.dna[,which(colSums(m.dna) > 0)])

# df.dna_blueprint.m <- subset(df.dna_blueprint, df.dna_blueprint$Target %in% v.enz.m.c)  
# tb.regs.m.c <- table((unlist(strsplit(df.dna_blueprint.m$TF, ", "))))
# tb.regs.m.c <- tb.regs.m.c[names(tb.regs.m.c) %in% v.regs.m]
# tb.regs.m.c <- tb.regs.m.c[tb.regs.m.c >= length(v.enz.m.c)]

### 
v.anova <- read.csv("anova_tests.csv", stringsAsFactors = FALSE, header = FALSE)
df.anova <- read.csv("anovaResults.csv", stringsAsFactors = FALSE)
agi <- v.map[as.character(df.anova$Cluster)]
df.anova <- df.anova[,-1]

which(v.map == "AT5G07080")

#df.anova <- df.anova[,c(49,50,51,87,88,89,109,110,111)]
#v.anova <- v.anova[,1][c(49,50,51,87,88,89,109,110,111)]

v.anova <- v.anova[,1]
names(df.anova) <- v.anova

m.anova <- as.matrix(df.anova)
rownames(m.anova) <- agi

m.anova[m.anova > 0.05] <- 10
m.anova[m.anova <= 0.05] <- 1
m.anova[m.anova == 10] <- 0

# l.anova <- vector(mode = "list", length = 2)
# l.anova[[1]] <- m.anova
# l.anova[[2]]
# 
# n.DE <- length(which(rowSums(m.anova) > 0))

# idx <- which(rowSums(m.anova) > 0)
# m.anova <- m.anova[idx, ]

print("Generate backbone expression structure")

# pearson correlation set
th.pcc.pos <- quantile(m.pcc, 0.95)
th.pcc.neg <- quantile(m.pcc, 0.05)

# > th.pcc.neg
# 5% 
# -0.363049 
# > th.pcc.pos
# 95% 
# 0.4028822 

# generate high confidence coexpression set
m.pcc.binary <- m.pcc
m.pcc.binary[m.pcc.binary <= th.pcc.neg] <- -1
m.pcc.binary[m.pcc.binary > th.pcc.neg] <- 0
m.pcc.binary[m.pcc >= th.pcc.pos] <- 1

m.pcc.binary <- m.pcc.binary[which(rowSums(m.pcc.binary) != 0), which(colSums(m.pcc.binary) != 0)]







## Finding enriched regulators for these targets

v.env <- c("AT1G66480", "AT5G07060", "AT5G07080", "AT3G46090", "AT3G03370", "AT3G03380")
v.env <- c(v.env, "AT5G07080") # early phosphate deficiency (3h)


# AT5G07080 # early phosphate deficiency (3h)
# AT1G66480 (6h)

v.env.onchip <- intersect(rownames(m.pcc.binary),v.env)

m.pcc[v.env.onchip,v.env.onchip]


l.ttests[[1]][v.env.onchip,]
l.ttests[[2]][v.env.onchip,]


# generate a differential expression set (co-specificity analysis)
m.anova.jc <- jaccard(m.anova)
rownames(m.anova.jc)  <- colnames(m.anova.jc) <- rownames(m.anova)
diag(m.anova.jc) <- 0

v.gns.overlap <- intersect(rownames(m.anova.jc), rownames(m.pcc.binary))

m.overlap <- m.pcc.binary[v.gns.overlap,v.gns.overlap] * as.matrix(m.anova.jc[v.gns.overlap,v.gns.overlap])



tfs.reg <- names(which(rowSums(m.dna[,intersect(colnames(m.dna),v.env.onchip)]) > 0))

m.dna[tfs.reg,intersect(colnames(m.dna),v.env.onchip)]


m.pcc.binary

m.anova[v.env,]

# diff exp regulator

colnames(df.geneExp.set) <- colnames(m.anova.set)

# df.geneExp.noLog <- exp(df.geneExp.set)

m.anova.set <- m.anova[,c(49,50,51,87,88,89,109,110,111),drop = FALSE]
m.anova.set <- m.anova.set[,c(1,4,7,2,5,8,3,6,9)]

m.anova.set <- m.anova.set[which(rowSums(m.anova.set) > 0),]
v.gns.diffexp <- rownames(m.anova.set)

df.geneExp.noLog.set <- df.geneExp.noLog

df.geneExp.set <- df.geneExp.set[v.gns.diffexp,]
df.geneExp.set <- df.geneExp.set[which(!is.na(rowSums(df.geneExp.set))),]


table(colSums(m.anova.set))
sum(m.anova.set)

df.geneExp.noLog["AT2G40080",]                                     
                                     
#m.anova.set <- m.anova[,c(4,5,6,49,50,51,87,88,89,109,110,111),drop = FALSE]

m <- m.anova[intersect(tfs, rownames(m.anova)),c(4,5,6,49,50,51,87,88,89,109,110,111),drop = FALSE]
idx.regs <- which(rowSums(m) > 0)

exps <- c("-P+Fe", "+P-Fe", "-P-Fe")
expSets <- vector(mode = "list", length = 3)
expSets[[1]] <- c(1,2,3)
expSets[[2]] <- c(4,5,6)
expSets[[3]] <- c(7,8,9)

for(i in 1:3){
  
  print(exps[i])
  
  idx.gns <- which(rowSums(m.anova.set[,expSets[[i]]]) > 0)
  
  print(length(idx.gns)) # number of diff exp genes
  
  # highlight 
  m.anova.set.i <- m.anova.set[,expSets[[i]]] # 
  
  df.geneExp.set.i <- df.geneExp.set[,expSets[[i]]]#[!is.na(df.geneExp.set[,expSets[[i]]])]
  
  # cluster the profiles - up->down, down-> up, 
  #  m.tmp <- kmeans(df.geneExp.set[,expSets[[i]]], 4)
  
  # expression fold change profiles
  plot(as.numeric(df.geneExp.set[idx.gns[1],expSets[[i]]]), type = "l", ylim = c(min(df.geneExp.set.i),max(df.geneExp.set.i)))
  
  for(j in 2:length(idx.gns)){
    lines(as.numeric(df.geneExp.set[idx.gns[j],expSets[[i]]]))
    #df.geneExp.noLog.set[idx.gns[j],expSets[[i]]]  
  }
  
  # filter for specific terms
  terms <- c("root", "DNA binding", "phosphate", "iron")
  
  go.set <- subset(df.geneontology, df.geneontology$AGI %in% names(idx.gns))

  # general genes with TF annotation
  tf.set <- subset(go.set, go.set$GO.Molecular.Function.Term %in% c("DNA binding")) # additional binding support ? / coexpression support?
  
  # binding support analysis
  idx.tfs <- which(names(idx.gns) %in% tfs)
  tf.dna <- names(idx.gns)[idx.tfs]
  
  
  if(length(tf.dna) > 0){
  
    for(r in 1:lenth(tf.dna))
      
      print(tf.dna[r])
    
      plot(as.numeric(df.geneExp.set.i[tf.dna[r],]), type = "l", ylim = c(-3.5,3.5))
        
      tgs.r <- names(which(m.dna[tf.dna,] == 1))
      tgs.r <- intersect(rownames(df.geneExp.set.i), tgs.r)

      #       
      #       tf.dna, tgs.r
      #       
      m.anova.set[c(tf.dna,tgs.r),]
      #       
      #       df.geneExp.set.i[tgs.r,]
      #       
      
      
      lines(as.numeric(df.geneExp.set.i[tgs.r[1],]), col = "red")
      lines(as.numeric(df.geneExp.set.i[tgs.r[2],]), col = "blue")
            
      
      
#       
#       
#       if(length(tgs.r) > 0){
#         for(j in 2:length(idx.gns)){
#           lines(as.numeric(df.geneExp.set[idx.gns[j],expSets[[i]]]), col = "")
#           #df.geneExp.noLog.set[idx.gns[j],expSets[[i]]]  
#         }
#         
#         
#       }
#       
      
      
  }
  
  
  
  idx.ifs <- which(tf.set$AGI %in% rownames(m.dna))
  
  test <- subset(df.geneontology, grepl(terms, df.geneontology$GO.Biological.Process.Term))
  
  df.geneExp.set["AT2G24850",]
  
  plot(as.numeric(df.geneExp.set["AT2G24850",expSets[[i]]]), type = "l", ylim = c(-3.5,3.5))
  
  
}

#m[idx.regs,]

for(i in 1:length(idx.regs)){
  
  tf <- names(idx.regs)[i]
  
  print(m[tf,])
  
  tg.i <- names(which(m.dna[tf,] == 1))
  
  #m[tf,]

  #(2/44) / (78 / 28000) =>  16.31702
  
#   > print(colSums(m.i))
#   Tuk_-P:+Fe:3hrs-+P:+Fe:3hrs Tuk_+P:-Fe:3hrs-+P:+Fe:3hrs Tuk_-P:-Fe:3hrs-+P:+Fe:3hrs Tuk_-P:+Fe:6hrs-+P:+Fe:6hrs Tuk_+P:-Fe:6hrs-+P:+Fe:6hrs Tuk_-P:-Fe:6hrs-+P:+Fe:6hrs Tuk_-P:+Fe:9hrs-+P:+Fe:9hrs Tuk_+P:-Fe:9hrs-+P:+Fe:9hrs Tuk_-P:-Fe:9hrs-+P:+Fe:9hrs 
#   0                           0                           0                           2                           0                           0                           0                           0                           0 
#   > colSums(m.anova[,c(49,50,51,87,88,89,109,110,111),drop = FALSE])
#   Tuk_-P:+Fe:3hrs-+P:+Fe:3hrs Tuk_+P:-Fe:3hrs-+P:+Fe:3hrs Tuk_-P:-Fe:3hrs-+P:+Fe:3hrs Tuk_-P:+Fe:6hrs-+P:+Fe:6hrs Tuk_+P:-Fe:6hrs-+P:+Fe:6hrs Tuk_-P:-Fe:6hrs-+P:+Fe:6hrs Tuk_-P:+Fe:9hrs-+P:+Fe:9hrs Tuk_+P:-Fe:9hrs-+P:+Fe:9hrs Tuk_-P:-Fe:9hrs-+P:+Fe:9hrs 
#   33                          20                          68                          44                          25                          23                          22                          18                          36 
#   
  colSums(m.anova[,c(49,50,51,87,88,89,109,110,111),drop = FALSE])
  
  m.i <- m.anova[intersect(tg.i, rownames(m.anova)),c(49,50,51,87,88,89,109,110,111),drop = FALSE]
  
  print(colSums(m.i))
  print(which(rowSums(m.i) > 0))
  #colSums(m.i)
  # enrichment
  
  if(tf %in% rownames(m.pcc.binary)){
  
    idx.pcc <- which(m.pcc.binary[tf, intersect(tg.i, colnames(m.pcc.binary))] != 0)
    tg.pcc <- m.pcc[tf, intersect(tg.i, colnames(m.pcc.binary))][idx.pcc]
    print(tg.pcc)  
  }
  
  AT3G11800
  
  print("")
  print("")
  print("")
  
}





write.csv(m[idx.regs,], "RegulatorAnova1.csv")

#idx.regs <- which(rowSums(m.anova[intersect(tfs, rownames(m.anova)),c(4,5,6,7),drop = FALSE]) > 0)

m <- m.anova[intersect(tfs, rownames(m.anova)),c(4,5,6,7),drop = FALSE]
idx.regs <- which(rowSums(m) > 0)

for(i in 1:length(idx.regs)){
  
  tf <- names(idx.regs)[i]
  
  
  if(m[tf,2] == 1){ 
    
    print(i)
    print(m[tf,])
    
    tg.i <- names(which(m.dna[tf,] == 1))
  
    #colSums(m.anova[,c(4,5,6,7),drop = FALSE])
    
    
    m.i <- m.anova[intersect(tg.i, rownames(m.anova)),c(4,5,6,7),drop = FALSE]
    colSums(m.i)
    
    tg.i <- names(which(m.i[,2] == 1))
    
    print(colSums(m.i))
    print(which(rowSums(m.i) > 0))
    #colSums(m.i)
    # enrichment
    
    if(tf %in% rownames(m.pcc.binary)){
      
      idx.pcc <- which(m.pcc.binary[tf, intersect(tg.i, colnames(m.pcc.binary))] != 0)
      tg.pcc <- m.pcc[tf, intersect(tg.i, colnames(m.pcc.binary))][idx.pcc]
      print(tg.pcc)  
    }
  }
    
  print("")
  print("")
  print("")
  
}






m[idx.regs,]

write.csv(m[idx.regs,], "RegulatorAnova2.csv")

table(rowSums(m.anova[intersect(tfs, rownames(m.anova)),c(49,50,51,87,88,89,109,110,111),drop = FALSE]))

###

m.dna



m.pcc[intersect(rownames(m.pcc), tfs.reg),intersect(colnames(m.dna),v.env.onchip)]

###
tfs
tgs
m.dna
l.ttests

v.gns.pcc <- rownames(m.pcc.binary)

# analyze individual sets #4 out of 4
for(i in 1:4){
  
  for(j in 1:3){
    
    m.tmp <- l.ttests[[1]][,1:4]
    tfs.tmp <- intersect(tfs, rownames(m.tmp))
    #m.tmp[,]
    
    m.tmp.reg <- m.tmp[tfs.tmp,]
    v.tmp.reg <- apply(m.tmp.reg, 1, function(m) length(which(m < 0.05))) # accumulation rank
    
    # rank based on differential expression
    table(v.tmp.reg)
    
    idx.tfSet <- (which(v.tmp.reg < 0.05)) 
    
    ## a) TFs expressed
    
    rank.tf <- m.tmp.reg
  
    tfs.tmp[idx.tfSet]
    
    tfs.sel <- names(v.tmp.reg[v.tmp.reg > 0])
    
    v.tmp.tg <- apply(m.tmp, 1, function(m) length(which(m < 0.05))) # accumulation rank
    table(v.tmp.tg)
  
    tgs.sel <- names(v.tmp.tg[v.tmp.tg > 0])
    
    tfs.sel <- intersect(tfs.sel,v.gns.pcc)
    tgs.sel <- intersect(tgs.sel,v.gns.pcc)
    
    table(m.pcc.binary[tfs.sel,tgs.sel])
    
    idx.tgSet <- (which(v.tmp.tg < 0.001)) 
    
    ## b) targets expressed (at or after TF)
    
    
    
    ## c) dna binding
    
    
    ## d) pcc binary confirmation
    
    
    # shifted or same column subset matrix #
    tgs.sel.dna <- intersect(tgs.sel, colnames(m.dna))
    
    
    # ranked based on number of targets or enrichment
    table(rowSums(m.dna[tfs.sel, tgs.sel.dna]))
    
    
    
    m.tmp
    
    m.dna[tfs.tmp[idx.tfSet],]
    
    table(rowSums(m.dna[tfs.tmp[idx.tfSet],]))
    
    # up, down regulation 
    
    
    
    
    
    l.ttests[[2]]

    
    
    
  }
  
  
  
}

## regulatory analysis to establish 
v.tfs.anova <- intersect(tfs, rownames(m.anova))

m.overlap.set <- m.overlap
idx <- which(rowSums(m.overlap.set) != 0)
m.overlap.set <- m.overlap.set[idx,idx]

if(FALSE)
  hist(m.overlap.set[m.overlap.set!=0], breaks = 100)

v.tfs.set <- intersect(v.tfs, rownames(m.overlap.set))


m.overlap.regulation <- m.overlap.set[v.tfs.set,]
rowSums(m.overlap.regulation)

df.regulation <- data.frame(TF = character(), targets = character(), TF_in_SNP = character(), targets_in_SNP = character())

for(i in 1:length(v.tfs.set)){
  
  print(v.tfs.set[i])
  targets <- names(which(m.overlap.set[v.tfs.set[i],] != 0 ))
  
  #TF_in_SNP <- character()
  #targets_in_SNP <- character()
  
  for(j in 1:length(l.gwas)){
    
    
    TF_in_SNP.tmp <- intersect(v.tfs.set[i], l.gwas[[j]]$SNP.relative.position)
    targets_in_SNP.tmp <- intersect(targets, l.gwas[[j]]$SNP.relative.position)
    
    if(length(TF_in_SNP.tmp) > 0)
      print("yes")
    
    print(targets_in_SNP.tmp)
    
    #TF_in_SNP <- paste(TF_in_SNP, paste(TF_in_SNP.tmp, names(l.gwas)[j], collapse = " - "),  collapse = "| ")
    #targets_in_SNP <- paste(targets_in_SNP, paste(paste(targets_in_SNP.tmp, collapse = "; "), names(l.gwas)[j], collapse = " - "),  collapse = "| ")
    
  }
  
  print("")
  #df.geneExp[v.tfs.set[i],]
  #df.regulation <- rbind(df.regulation, data.frame(TF =  v.tfs.set[i], targets = paste(names(which(m.overlap.set[v.tfs.set[i],] != 0)), collapse = "; ")), TF_in_SNP = TF_in_SNP, targets_in_SNP = targets_in_SNP)
}

write.csv(df.regulation, "df.regulation_novel.csv", row.names = FALSE)

subset(df.geneontology, df.geneontology$AGI %in% v.tfs.set)$GO.Biological.Process.Term

###

# library(igraph)
# library(parmigene)

# normalize coexpression around zero mean 

# test: alternatives for sparse graph generation #
# m.clr <- parmigene::clr(m.pcc)
# m.network <- m.pcc # m.clr
# th <- 0.9 #quantile(m.pcc, 0.99)
# m.network[m.network <= th] <- 0
# m.network[m.network > 0] <- 1
# 
# g <- graph_from_adjacency_matrix(m.network)
# wc <- infomap.community(g)
# 
# # v.modules <- unique(membership(wc))
# tb.modules <- table(membership(wc))
# v.modules <- names(tb.modules[tb.modules > 1])


p.modules <- rep(1, length(v.modules))
names(p.modules) <- v.modules

# p.dna.modules <- rep(1, length(v.modules))
# names(p.dna.modules) <- v.modules

# m = 1 - non assigments
for(m in 2:length(v.modules)){
  
  module <- names(which(v.module_membership == v.modules[m]))
  # module <- names(which(membership(wc) == v.modules[m]))
  v.gns.anova <- intersect(rownames(m.anova), module)
  
  if(length(v.gns.anova) > 0){
    # disregarding splitting of test combinations
    hitInSample = length(which(rowSums(m.anova[v.gns.anova,,drop=FALSE]) > 0))
    
    hitInPop = n.DE 
    failInPop = dim(m.anova)[1] - hitInPop 
    sampleSize = length(v.gns.anova)
    
    if(hitInSample > 1){
      p.modules[m] <- phyper(hitInSample - 1, hitInPop, failInPop, sampleSize, lower.tail= FALSE);  
    }
    
  }
  
}

# significantly enriched (differentially genes) 
v.sigModules <- names(which(p.modules <= 0.05))

df.module <- data.frame(id = character(), gns = character(), diffExpGns = character(), GO = character())

for(m in 1:length(v.sigModules)){
  
  print(paste("Module:", m))
  
  module <- names(which(v.module_membership == v.sigModules[m]))
  #module <- names(which(membership(wc) == v.sigModules[m]))
  v.gns.anova <- intersect(rownames(m.anova), module)
  
  # m.anova[v.gns.anova,,drop=FALSE]
  
  #print(rowSums(m.anova[v.gns.anova,,drop=FALSE]))
 
  v.gns.diff <- names(which(rowSums(m.anova[v.gns.anova,,drop=FALSE]) == 1))
  m.anova[v.gns.diff,]
  
  # correlation scores
  m.pcc[v.gns.diff,v.gns.diff]
  
  ##
  conds <- c("Ctrl", "+P+Fe(3h)", "+P+Fe(6h)", "+P+Fe(9h)", "-P+Fe(3h)", "-P+Fe(6h)", "-P+Fe(9h)", "+P-Fe(3h)", "+P-Fe(6h)", "+P-Fe(9h)", "-P-Fe(3h)", "-P-Fe(6h)", "-P-Fe(9h)")
  
  df.geneExp[v.gns.diff,]
  
  # quad plot # 
  
  
  # additional information #
  v.gwas <- l.gwas[[1]]
  v.gwas <- unique(v.gwas$SNP.relative.position)
  
  v.gwas.m <- intersect(v.gwas, v.gns.anova)
  
  if(length(v.gwas.m) > 0){
    # disregarding splitting of test combinations
    hitInSample = length(v.gwas.m)
    hitInPop = length(v.gwas) 
    failInPop = dim(m.anova)[1] - hitInPop 
    sampleSize = length(module)
    
    if(hitInSample > 1){
      p.modules[m] <- phyper(hitInSample - 1, hitInPop, failInPop, sampleSize, lower.tail= FALSE);  
    }
    
  }
  
  
  # GO 
  v.go <- subset(df.geneontology, df.geneontology$AGI %in% v.gns.diff)$GO.Biological.Process.Term
  tb.go <- table(unlist(sapply(v.go, function(m) strsplit(m, "; "))))
  
  th.min <- 2
  tb.go <- tb.go[tb.go >= th.min]
  
  print("Gene ontology (biological process) - differential expressed genes")
  print(tb.go)
  
  ## 
  
  v.go <- subset(df.geneontology, df.geneontology$AGI %in% v.gns.anova[!v.gns.anova %in% v.gns.diff])$GO.Biological.Process.Term
  tb.go <- table(unlist(sapply(v.go, function(m) strsplit(m, "; "))))
  
  th.min <- 2
  tb.go <- tb.go[tb.go >= th.min]
  
  #print("Gene ontology (biological process) - remaining genes")
  #print(tb.go)

  # dna blueprint 
    # df.dna.m <- subset(df.dna, df.dna$Target.gene %in% v.gns.diff)
    # table(df.dna.m$Putative.TF)
    # 
    # df.dna.m <- subset(df.dna, df.dna$Target.gene %in% v.gns.anova[!v.gns.anova %in% v.gns.diff])
    # table(df.dna.m$Putative.TF)
    # 
  
  # fold change of expression (mean gene expression )

    # nrow(subset(df.dna, df.dna$Putative.TF == "AT2G45660"))
    # 
    # v.go <- subset(df.geneontology, df.geneontology$AGI %in% "AT1G13260")$GO.Biological.Process.Term
    # tb.go <- table(unlist(sapply(v.go, function(m) strsplit(m, "; "))))
    # 
  
  # regulation plot - split into 4 figures per experiment and time series
  # df.exp.m <- subset(df.exp, df.exp$AGI == "AT1G13260")
  # plot(as.numeric(df.exp.m[,1:12]), type = "l", ylim = c(0,10))
  # 
  # df.exp.m <- subset(df.exp, df.exp$AGI %in% v.gns.diff)
  # lines(as.numeric(df.exp.m[1,1:12]), col = "red")
  # lines(as.numeric(df.exp.m[2,1:12]), col = "blue")
  # 
  
  tmp <- data.frame(id = v.sigModules[m], gns = paste(v.gns.anova, collapse = ","), diffExpGns = paste(v.gns.diff, collapse = ","), GO = paste(names(tb.go), collapse = ","))
  
  df.module <- rbind(df.module, tmp)
    
  
  
  print("") 
}

write.csv(df.module, paste("moduleAnalysis.csv"))

#source("https://bioconductor.org/biocLite.R")
#biocLite("genefilter")


# expression analysis #\
df.exp <- read.csv("geneExpData.csv", header = TRUE)
df.exp["AGI"] <- v.map[as.character(df.exp$ClusterID)]
df.exp <- df.exp[,-1]


v.control <- c(1,14,27) # control

l.timeseries <- vector(mode = "list", length = 3)

# 3 hours
l.timeseries[[1]] <- vector(mode = "list", length = 4)
l.timeseries[[1]][[1]] <- c(2,15,28) # +P+Fe
l.timeseries[[1]][[2]] <- c(3,16,29) # -P+Fe
l.timeseries[[1]][[3]] <- c(4,17,30) # +P-Fe
l.timeseries[[1]][[4]] <- c(5,18,31) # -P-Fe

# 6 hours
l.timeseries[[2]] <- vector(mode = "list", length = 4)
l.timeseries[[2]][[1]] <- c(6,19,32) # +P+Fe
l.timeseries[[2]][[2]] <- c(7,20,33) # -P+Fe
l.timeseries[[2]][[3]] <- c(8,21,34) # +P-Fe
l.timeseries[[2]][[4]] <- c(9,22,35) # -P-Fe

# 9 hours
l.timeseries[[3]] <- vector(mode = "list", length = 4)
l.timeseries[[3]][[1]] <- c(10,23,36) # +P+Fe
l.timeseries[[3]][[2]] <- c(11,24,37) # -P+Fe
l.timeseries[[3]][[3]] <- c(12,25,38) # +P-Fe
l.timeseries[[3]][[4]] <- c(13,26,39) # -P-Fe


l.DE.conditions <- vector(mode = "list", length = 4)
l.DE.conditions[[1]] <- matrix(0, nrow = length(df.exp$AGI), ncol = 3, dimnames = list(df.exp$AGI, c(3,6,9)))
l.DE.conditions[[2]] <- matrix(0, nrow = length(df.exp$AGI), ncol = 3, dimnames = list(df.exp$AGI, c(3,6,9)))
l.DE.conditions[[3]] <- matrix(0, nrow = length(df.exp$AGI), ncol = 3, dimnames = list(df.exp$AGI, c(3,6,9)))
l.DE.conditions[[4]] <- matrix(0, nrow = length(df.exp$AGI), ncol = 3, dimnames = list(df.exp$AGI, c(3,6,9)))
names(l.DE.conditions) <- c("+P+Fe", "-P+Fe", "+P-Fe", "-P-Fe")

for(i in 1:length(l.DE.conditions)){
  
  for(j in 1:3){
    
    v.control <- l.timeseries[[j]][[1]]
    m <- as.matrix(df.exp[,c(v.control,l.timeseries[[j]][[i]])])
    f <- factor(c(rep("Control", 3), rep("Treatment", 3)))
    ttests <- rowttests(m, f)
    
    # ttests$p.value <- p.adjust(ttests$p.value, "fdr")
    idx <- which(ttests$p.value < 0.01)
    l.DE.conditions[[i]][idx, j] <- ttests$dm[idx]
    l.DE.conditions[[i]][, j] <- ifelse(l.DE.conditions[[i]][, j] > 0, -1, ifelse(l.DE.conditions[[i]][, j] < 0, 1, 0))
    
  }
  
}

lapply(l.DE.conditions, table)

###


jaccard <- function(m) {
  ## common values:
  A = tcrossprod(m)
  ## indexes for non-zero common values
  im = which(A > 0, arr.ind=TRUE)
  ## counts for each row
  b = rowSums(m)
  ## only non-zero values of common
  Aim = A[im]
  ## Jacard formula: #common / (#i + #j - #common)
  J = sparseMatrix(
    i = im[,1],
    j = im[,2],
    x = Aim / (b[im[,1]] + b[im[,2]] - Aim),
    dims = dim(A)
  ) 
  return( J )
}


# m.test.global <- c()
# for(i in 1:12){
#   
#   #   
#   #   Value <- c(1,2,5,3,2,1,1,3,2,1,4,3,6,5,2,6,1,6,5,4,9,6,7,7,5,1,8,9,6,5)
#   #   Group <- factor(c(rep(1,10),rep(2,10),rep(3,10)))
#   #   
#   #   data <- data.frame(categ, resp)
#   #   
#   #   kruskal.test(Value ~ Group, data=data)
#   #   pairwise.wilcox.test(resp, categ, p.adj="holm", exact=F)
#   #     
#   m.test.global <- cbind(m.test.global, apply(df.geneExp, 1, function(m) {t.test(as.numeric(m[c(1+i,14+i,27+i)]), as.numeric(m[c(1,14,27)]))$p.value}))
# 
#   
# #   m.test <- cbind(m.test, apply(df.geneExp, 1, function(m) {
# #                                                           resp<-c(as.numeric(m[c(1+i,14+i,27+i)]), as.numeric(m[c(1,14,27)]))
# #                                                           categ<-as.factor(rep(c("A","B"),times=1,each=3))
# #                                                           ifelse(kruskalmc(resp, categ, probs = 0.05, cont="two-tailed")$dif.com$difference, 1, 0)
# #                                                         }))
#                                                           
#   #    m <- df.geneExp[1,]
#   #     
#   #     wilcox.test(as.numeric(m[c(1,14,27)]), as.numeric(m[c(1+i,14+i,27+i)]))$p.value))
#   #   
#   #   
#   
#   #m.test <- cbind(m.test, apply(df.geneExp, 1, function(m) wilcox.test(as.numeric(m[c(1,14,27)]), as.numeric(m[c(1+i,14+i,27+i)]))$p.value))
#   #m.test <- cbind(m.test, apply(df.geneExp, 1, function(m) t.test((as.numeric(m[c(1+i,14+i,27+i)])), (as.numeric(m[c(1,14,27)])) )$p.value ))
#   # m.test <- cbind(m.test, apply(df.geneExp, 1, function(m) (mean(as.numeric(m[c(1+i,14+i,27+i)])) - mean(as.numeric(m[c(1,14,27)])) )))
# }


