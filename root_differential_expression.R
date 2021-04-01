source("MERIT_DifferentialExpression.R")
source("MERIT_DNA.R")

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

# m.cor <- cor(t(m.fc.root))

v.treatment_buildingblocks.root <- v.treatment_buildingblocks[i.set]
v.repeats.root <- v.repeats[i.set]



# root specific analysis of the REMAINING gene sets!!! 



























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
