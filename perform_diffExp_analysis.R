





precompute_and_store_diffExp_and_foldchange_matrices <- function(folder = "/tmp"){
  
  
  filename.annotation = "A:/junkDNA.ai/datasets/experiment_annotation_He_et_al_2015_reannotation.txt"
  filename.geneExpression = "A:/junkDNA.ai/MERIT/datasets/gene_expression/GSE69995_re-analyzed_data_matrix.txt"
  
  
  v.ath.Stress.treatment <- c("GSE5620", "GSE5621", "GSE5622", "GSE5623", "GSE5624", "GSE5625", "GSE5626", "GSE5627", "GSE5628")
  
  
  df.annotation <- read.csv(filename.annotation,  header = TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
  v.colnames_mandatory = c("sample_id", "series_id", "condition_treatment_1", "condition_treatment_2", "condition_tissue","replicate_group_w_control", "sublevel_series_id", "unique_ID")
  if(!all(v.colnames_mandatory %in% names(df.annotation))){
    stop(paste("could not find all mandatory columns in file:", paste(v.colnames_mandatory, collapse = ", ")))
  }
  
  tb.treatments = table(c(df.annotation$condition_treatment_1, df.annotation$condition_treatment_2))
  tb.treatments = tb.treatments[names(tb.treatments) != ""]
  
  # tissue and treatment distributions
  tb.tissues <- table(df.annotation$condition_tissue)
  v.conditionGroups = names(tb.treatments)
  
  # loading gene expressino matrix 
  m.expression <- read.table(filename.geneExpression, row.names = 1, header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
  m.expression <- as.matrix(m.expression)
  df.annotation <- read.csv(filename.annotation,  header = TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
  
  
  tb.conditions = table(df.annotation$condition_treatment_1)
  tb.conditions[names(table(df.annotation$condition_treatment_2))] = tb.conditions[names(table(df.annotation$condition_treatment_2))] +  table(df.annotation$condition_treatment_2)
  print(tb.conditions)
  print(table(df.annotation$condition_tissue))
  
  # write.table(m.expression, "datasets/gene_expression_He_et_al_2015_subset.txt", row.names = F, sep = "\t")
  
  # split the experiments and tissues 
  df.annotation["key.exp"] <- NA
  df.exp.pairs <- unique(df.annotation[,c("tissue","series_id")])
  for(i in 1:nrow(df.exp.pairs)){
    i.set <- which(df.annotation$tissue == df.exp.pairs[i,1] & df.annotation$series_id == df.exp.pairs[i,2])
    df.annotation$key.exp[i.set] <- i
  }
  
  m.expression <- m.expression[,df.annotation$sample_id]
  v.series_ids <- unique(df.annotation$series_id)
  
  
  v.series_ids <- v.series_ids[!v.series_ids %in% v.ath.Stress.treatment]
  v.series_sets <- c()
  
  
  v.exp_unique_ID = c()
  
  ####

  m.de <- c()
  m.fc <- c()
  m.directionality <- c()
  v.treatment_buildingblocks <- c()
  v.repeats <- c()
  
  exp.counter <- 1
  
  message(paste("Compute differential expression based on Welch t-test and prepare treatment gene expression matrix for two-way anova analysis"))
  
  pb <- txtProgressBar(min = 0, max = length(v.series_ids), style = 3)
  for(i in 1:length(v.series_ids)){
    
    setTxtProgressBar(pb, i)
    
    df.annotation.i <- subset(df.annotation, df.annotation$series_id == v.series_ids[i])
    
    v.exp.i <- unique(df.annotation.i$replicate_group_w_control)
    v.exp.i <- v.exp.i[v.exp.i != 0]
    v.set.i <- unique(df.annotation.i$sublevel_series_id)
    
    for(j in 1:length(v.set.i)){
      
      for(k in 1:length(v.exp.i)){
        
        df.annotation.control <- subset(df.annotation.i, df.annotation.i$sublevel_series_id ==  v.set.i[j] & df.annotation.i$replicate_group_w_control == 0)
        df.annotation.j <- subset(df.annotation.i, df.annotation.i$sublevel_series_id ==  v.set.i[j] & df.annotation.i$replicate_group_w_control == v.exp.i[k])
        
        if(nrow(df.annotation.j) > 0){
          
          ## Individual t-test p-values
          X.treatment <- t(m.expression[,df.annotation.j$sample_id])
          X.control <- t(m.expression[,df.annotation.control$sample_id])
          
          ttpv <- numeric(dim(X.treatment)[2])
          names(ttpv) <- colnames(X.treatment)
          
          v.fc <- numeric(dim(X.treatment)[2])
          names(v.fc) <- colnames(X.treatment)
          
          v.directionality <- numeric(dim(X.treatment)[2])
          names(v.directionality) <- colnames(X.treatment)
          
          v.series_sets <- c(v.series_sets, v.series_ids[i])
          v.repeats <- c(v.repeats, seq(1:dim(X.treatment)[1]))
       
         
          v.treatment_buildingblocks <- c(v.treatment_buildingblocks, rep(exp.counter,dim(X.treatment)[1]))
          exp.counter <- exp.counter + 1
          
          for(l in 1:ncol(X.treatment)){ # per gene
            if(sd(X.treatment[,l]) != 0 | sd(X.control[,l]) != 0){
              tt <- t.test(X.treatment[,l], X.control[,l]) # welch t-test per gene 
              ttpv[l] = unlist(tt$p.value)
              v.fc[l] <- mean(X.treatment[,l]) - mean(X.control[,l])
              if(mean(X.treatment[,l]) == mean(X.control[,l])){
                v.directionality[l] = 0
              }else if(mean(X.treatment[,l]) > mean(X.control[,l])){
                v.directionality[l] = 1
              }else{
                v.directionality[l] = -1
              }
            }else{
              ttpv[l] <- 1
            }
          }
          
          
          v.exp_unique_ID = c(v.exp_unique_ID,  unique(df.annotation.j$unique_ID))
          
          
          m.de <- cbind(m.de, ttpv)
          names(v.fc) <- v.series_ids[i]
          m.fc <- cbind(m.fc, v.fc)
          names(v.directionality) <- v.series_ids[i]
          m.directionality <- cbind(m.directionality, v.directionality)
        }  
      }
    }
  }
  close(pb)
  
  ####
  
  df.annotation.control <- subset(df.annotation, df.annotation$series_id == v.ath.Stress.treatment[1])

  pb <- txtProgressBar(min = 0, max = 9, style = 3)
  for(i in 2:9){
    
    df.annotation.i <- subset(df.annotation, df.annotation$series_id == v.ath.Stress.treatment[i])
    
    v.exp.i <- unique(df.annotation.i$replicate_group_w_control)
    v.exp.i <- v.exp.i[v.exp.i != 0]
    v.set.i <- unique(df.annotation.i$sublevel_series_id)
    
    for(j in 1:length(v.set.i)){
      
      for(k in 1:length(v.exp.i)){
      
        df.annotation.control.j <- subset(df.annotation.control, df.annotation.control$sublevel_series_id ==  v.set.i[j] & df.annotation.control$replicate_group_w_control == v.exp.i[k])
        df.annotation.j <- subset(df.annotation.i, df.annotation.i$sublevel_series_id ==  v.set.i[j] & df.annotation.i$replicate_group_w_control == v.exp.i[k])
        
        if(nrow(df.annotation.j) > 0){
          
          ## Individual t-test p-values
          X.treatment <- t(m.expression[,df.annotation.j$sample_id])
          X.control <- t(m.expression[,df.annotation.control.j$sample_id])
          
          ttpv <- numeric(dim(X.treatment)[2])
          names(ttpv) <- colnames(X.treatment)
          
          v.fc <- numeric(dim(X.treatment)[2])
          names(v.fc) <- colnames(X.treatment)
          
          # m.fc.t <- matrix(0, dim(X.treatment)[2], dim(X.treatment)[1])
          
          v.series_sets <- c(v.series_sets, v.ath.Stress.treatment[i])
          
          
          v.exp_unique_ID = c(v.exp_unique_ID,  unique(df.annotation.j$unique_ID))

          
          v.repeats <- c(v.repeats, seq(1:dim(X.treatment)[1]))
          v.treatment_buildingblocks <- c(v.treatment_buildingblocks, rep(exp.counter,dim(X.treatment)[1]))
          exp.counter <- exp.counter + 1
          
          for(l in 1:ncol(X.treatment)){
            
            if(sd(X.treatment[,l]) != 0 | sd(X.control[,l]) != 0){
              
              tt <- t.test(X.treatment[,l],X.control[,l])
              ttpv[l] = unlist(tt$p.value)
              
              v.fc[l] <- mean(X.treatment[,l]) - mean(X.control[,l])
              
              # for(r in 1:dim(X.treatment)[1]){
              #   m.fc.t[l,r] <- X.treatment[r,l] - mean(X.control[,l])  
              # }
              
              if(mean(X.treatment[,l]) == mean(X.control[,l])){
                v.directionality[l] = 0
              }else if(mean(X.treatment[,l]) > mean(X.control[,l])){
                v.directionality[l] = 1
              }else{
                v.directionality[l] = -1
              }
              
            }else{
              ttpv[l] <- 1
            }
          }
          
          m.de <- cbind(m.de, ttpv)
          
          names(v.fc) <- v.ath.Stress.treatment[i]
          # colnames(m.fc.t) <- rep(v.series_ids[i], dim(m.fc.t)[2])
          
          m.fc <- cbind(m.fc, v.fc)
          
          names(v.directionality) <- v.ath.Stress.treatment[i]
          m.directionality <- cbind(m.directionality, v.directionality)
          
        } 
        
      }
    }
  }
  close(pb)
  
  
  ####

  rownames(m.fc) <- rownames(m.directionality) <- rownames(m.de)
  colnames(m.fc) <- colnames(m.directionality) <- colnames(m.de) <- v.exp_unique_ID
  
  
  write.table(m.de, "A:/junkDNA.ai/datasets/m.pvalue_differentialExpression.txt", row.names = F, col.names = F, sep = "\t")
  write.table(m.fc, "A:/junkDNA.ai/datasets/m.foldChange_differentialExpression.txt", row.names = F, col.names = F, sep = "\t")
  
  df.annotation_series = unique(df.annotation[,c("series_id","condition_treatment_1", "condition_treatment_2", "condition_tissue", "unique_ID")])
  write.table(df.annotation_series, "A:/junkDNA.ai/datasets/experiment_annotation.txt", row.names = F, sep = "\t")
  
  genes = rownames(m.de)
  experiment_series_id = colnames(m.de)
  write.table(genes, "A:/junkDNA.ai/datasets/genes.txt", row.names = F, col.names = F, sep = "\t")
  write.table(experiment_series_id, "A:/junkDNA.ai/datasets/experiment_ids.txt", row.names = F, col.names = F, sep = "\t")
  
  
}