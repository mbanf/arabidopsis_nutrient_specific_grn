#source("http://www.bioconductor.org/biocLite.R")
#biocLite("tigre")

setwd("Documents/postdoctoral_work/Projects/Hatem/")

library(tigre)
library(reshape2)


# random forest regression 

df.dna_binding <- readRDS("df.dna_binding_dapSeq.rds")

m.dna <- acast(df.dna_binding, tf.at_id~target.at_id)
m.dna[is.na(m.dna)] <- 0
m.dna[m.dna != "0"] <- 1
class(m.dna) <- "numeric"


#### 

#### Differential expression data

df.diffexp <- read.csv("diffexp_highConfidence.csv", stringsAsFactors = FALSE)
df.geneontology <- read.csv("Table annotation_INITIAL.csv", stringsAsFactors = FALSE)

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


v.anova <- read.csv("anova_tests.csv", stringsAsFactors = FALSE, header = FALSE)
df.anova <- read.csv("anovaResults.csv", stringsAsFactors = FALSE)
agi <- v.map[as.character(df.anova$Cluster)]
df.anova <- df.anova[,-1]

v.anova <- v.anova[,1]
names(df.anova) <- v.anova

m.anova <- as.matrix(df.anova)
rownames(m.anova) <- agi

m.anova[m.anova > 0.05] <- 10
m.anova[m.anova <= 0.05] <- 1
m.anova[m.anova == 10] <- 0

m.anova.set <- m.anova[,c(49,50,51,87,88,89,109,110,111),drop = FALSE]
m.anova.set <- m.anova.set[,c(1,4,7,2,5,8,3,6,9)]

m.anova.set <- m.anova.set[which(rowSums(m.anova.set) > 0),]
v.gns.diffexp <- rownames(m.anova.set) 

###



# automatic identification of the best log likelihood - rank, compare with the expression sets we had 
# combine with t tests - prior of tf and target involvement


m.expression <- read.table("GSE69995_re-analyzed_data_matrix.txt", row.names = 1, header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
m.expression <- as.matrix(m.expression)

m.expression.diffExp <- m.expression[intersect(rownames(m.expression), v.gns.diffexp), ]


df.exp_additional <- read.table("geneExpressionAnnotation_hatem.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)


# select the experiment 
df.exp_additional.set <- subset(df.exp_additional, df.exp_additional$series_id ==  "GSE10502")

m.expression.diffExp.set <- m.expression.diffExp[,df.exp_additional.set$sample_id]

#df.exp_additional.set$source.name

times <- c(0, 3, 6, 12, 24, 48, 72)

df.exp <- processRawData(m.expression.diffExp.set, times, experiments=rep(1:2, each=7), is.logged = TRUE)#, do.normalisation = ifelse(is.logged, TRUE, FALSE))


# The probe identifier for TF 'twi'
tfs <- (intersect(v.gns.diffexp, rownames(m.dna)))
tfs <- intersect(tfs, rownames(m.expression.diffExp.set))

# The probe identifier for the target gene
targets <- (intersect(v.gns.diffexp, colnames(m.dna)))
targets <- intersect(targets, rownames(m.expression.diffExp.set))



m.loglik <- matrix(0, length(tfs), length(targets), dimnames = list(tfs, targets))

for(i in 1:length(tfs)){

  #v.loglik <- numeric(length(targets))
  #names(v.loglik) <- targets
  
  for(j in 1:length(targets)){ # parallelize

    # Learn the model using only one of the 3 repeats in the data
    model <- GPLearn(df.exp[,1:7],
                     TF=tfs[i], targets = targets[j],
                     useGpdisim=TRUE, quiet=TRUE)
      
    #m.expression.diffExp.set[targets[1], ]
    
    # Display the model parameters
    show(model)
    
    GPPlot(model)
    
    m.loglik[i,j] <- - model@model$llscore
    
  }
  

}

#
m.dna[tfs[i], targets[j]]
# sum(m.dna[tfs[i],])
