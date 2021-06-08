# source("MERIT_DifferentialExpression.R")
# source("MERIT_DNA.R")

m.de.rf <- readRDS("data/m.de_0502.rds")
m.fc.rf <- readRDS("data/m.fc_0502.rds")
v.series_sets <- readRDS("data/v.series_sets_0502.rds")
v.repeats <- readRDS("data/v.repeats_0502.rds")
v.treatment_buildingblocks <- readRDS("data/v.treatments_0502.rds")

# raw expression data actually not neede for random forest regression #
# m.expression <- as.matrix(read.table("data/GSE69995_re-analyzed_data_matrix.txt", row.names = 1, header = TRUE, sep = "\t", quote = ""))

df.root_selection <- read.csv("data/expMetaRoot.txt", sep = ";", header = FALSE, stringsAsFactors = FALSE)
v.series.root <- unique(df.root_selection$V11)
v.series.root <- v.series.root[v.series.root != "characteristics..Original.series.ID"]

#v.samples.root <- unique(df.root_selection$V2)
#v.samples.root <- v.samples.root[v.samples.root != "Sample.name"]

#v.samples.root <- intersect(colnames(m.expression), v.samples.root)
# m.expression.root <- m.expression[,v.samples.root] # expression matrix for random forest 

i.set <- which(colnames(m.fc.rf) %in% v.series.root)

v.series.root <- intersect(colnames(m.de.rf), v.series.root)
m.de.root <- m.de.rf[,v.series.root]

m.fc.root <- m.fc.rf[,i.set] # select from log foldchange matrix

v.treatment_buildingblocks.root <- v.treatment_buildingblocks[i.set]
v.repeats.root <- v.repeats[i.set]



# root specific analysis of the REMAINING gene sets!!! 
# m.cor <- cor(t(m.fc.root))
m.fc.root <- m.fc.root[intersect(rownames(m.fc.root), gns.DE),] # 175 genes 
# m.expression.root <- m.expression.root[intersect(rownames(m.fc.root), gns.DE),]



tfs.rf <- intersect(rownames(m.fc.root), names(v.regulators))
tgs.rf <- rownames(m.fc.root)



# v.tf_families[tfs.DE.root]
# m.fc.root was previously
source("utils.R")
m.rf <- compute_randomforest_based_GRN(mat.expression=m.fc.root, k="sqrt", nb.trees=1000, set.regulators = tfs.rf, set.genes = tgs.rf, seed=1234, importance.measure = "impurity", n.cpus = 2)
# saveRDS(m.rf, "m.RF_1000_032621.rds") # 032621 - work with the filtered 380 experiments dataset 
# saveRDS(m.rf, "tmp/m.RF.rds")
# m.rf <- readRDS("data/m.RF_500.rds")

# th = quantile(m.rf, 0.9)
p.rf <- ecdf(m.rf)
df.rf.grn <- melt(m.rf)
names(df.rf.grn) <- c("TF", "TG", "val")
df.rf.grn <- subset(df.rf.grn, df.rf.grn$val > 0)
df.rf.grn["val"] <- p.rf(df.rf.grn$val)
