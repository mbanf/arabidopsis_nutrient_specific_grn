
# using identifications based on iTAK platform to identify regulator families involved 
tf.families = read.table("data/Arabidopsis_thaliana-all.txt", header = FALSE, fill = TRUE, stringsAsFactors = FALSE)[,1:2] 
tf.families$V1 <- gsub("\\..*", "",tf.families$V1)
tf.families <- unique(tf.families)

v.tf_families <- tf.families$V2
names(v.tf_families) <- tf.families$V1

i.set <- which(names(v.tf_families) %in% rownames(m.anova.set)) # transcriptional regulators as part of the differential expression dataset
v.regulators <-  v.tf_families[i.set]#[!grepl("PPC:", v.tf_families[i.set])]

# filter protein kinases etc.
v.regulators <- v.regulators[!grepl("RLK", v.regulators) & !grepl("PPC", v.regulators) & !grepl("AGC", v.regulators) & !grepl("Group", v.regulators)  & !grepl("mTERF", v.regulators) & !grepl("Tify", v.regulators) & !grepl("Aur", v.regulators)  & !grepl("CAMK", v.regulators)  & !grepl("TRAF", v.regulators) & !grepl("MED7", v.regulators)]

# saveRDS(names(v.regulators), "v.regulators.rds")


### dna binding ###
