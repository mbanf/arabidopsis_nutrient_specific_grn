## nutrition specific gene expression data 
exps <- c("-P / +Fe (3h)", "-P / +Fe (6h)", "-P / +Fe (9h)",  "+P / -Fe (3h)", "+P / -Fe (6h)", "+P / -Fe (9h)", "-P / -Fe (3h)", "-P / -Fe (6h)", "-P / -Fe (9h)")
ill_defined <- c("UNIDENTIFIED", "AMBIGUOUS", "RESCUE")

df.geneontology <- read.csv(paste(folder_data, "Table annotation_INITIAL.csv", sep = "/"), stringsAsFactors = FALSE) #read.table("Table annotation_INITIAL.txt", header = TRUE, sep = "\t", quote = "")

v.map <- df.geneontology$AGI
names(v.map) <- df.geneontology$Transcript.Cluster.ID

v.anova <- read.csv(paste(folder_data, "anova_tests.csv", sep = "/"), stringsAsFactors = FALSE, header = FALSE)
df.anova <- read.csv(paste(folder_data, "anovaResults.csv", sep = "/"), stringsAsFactors = FALSE)
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


names(exps) <- colnames(m.anova.set)

saveRDS(gns.DE, paste(folder_tmp, "gns.DE.rds", sep = "/")) # needed for a faster TFBS estimation

### raw gene expression data ### 
df.geneExp <- read.csv(paste(folder_data, "geneExpData.csv", sep = "/"), header = TRUE, stringsAsFactors = FALSE)
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


#### gene set enrichment analysis

df.de <- df.geneExp.set[gns.DE,c(1,4,7,2,5,8,3,6,9)]
df.de <- na.omit(df.de)
m.de <- as.matrix(df.de, drop = FALSE)
colnames(m.de) <- exps

F_m2 <- m.de





