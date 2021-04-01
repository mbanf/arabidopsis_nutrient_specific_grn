


df.dna <- c()

# https://www.nature.com/articles/srep25164
if(FALSE){
  
  df.dna.conservation.selection  <- read.csv("data/Yu2016-srep25164-s2.csv", stringsAsFactors = FALSE)
  df.dna.conservation.selection <- df.dna.conservation.selection[,c(1,6)]
  
  df.dna.conservation <- data.frame(TF = character(), Target = character(), stringsAsFactors = FALSE)
  
  pb <- txtProgressBar(min = 0, max = nrow(df.dna.conservation.selection), style = 3)
  
  for(i in 1:nrow(df.dna.conservation.selection)){
    
    setTxtProgressBar(pb, i)
    
    tfs.i <- unlist(strsplit(df.dna.conservation.selection$TF[i], ", "))
    tg.i <- df.dna.conservation.selection$Target[i]
    
    for(j in 1:length(tfs.i)){
      df.dna.conservation <- rbind(df.dna.conservation, data.frame(TF = tfs.i[j], Target = tg.i))  
    }
    
  }
  
  close(pb)
  
  saveRDS(df.dna.conservation, "data/df.dna.conservation.rds")
  
}else{
  df.dna.conservation <- readRDS("data/df.dna.conservation.rds") 
}
df.dna <- df.dna.conservation

# DAP-seq
df <- readRDS("data/df.dna_binding_dapSeq.rds")
colnames(df) <- c("TF", "Target")
df.dna <- rbind(df.dna, df)

# ATREGNET
df <- read.csv("data/AtRegNet.csv", sep = ",", fill = T,  header = T, row.names=NULL,stringsAsFactors=FALSE)
df<- subset(df, df$TargetType == "Confirmed")
df<- subset(df, df$IsTF == "Direct")
df<- df[,c(2,5)]
df <- apply(df,2,toupper)
colnames(df) <- c("TF", "Target")
df.dna <- rbind(df.dna, df)

# ATRM
df <- read.csv("data/Regulations_in_ATRM.csv", sep = "\t", fill = T, header = T, row.names=NULL,stringsAsFactors=FALSE)
df<- df[,c(1,2)]
colnames(df) <- c("TF", "Target")
df.dna <- rbind(df.dna, df)

df.dna <- unique(df.dna)

###

df.dna.de <- subset(df.dna, df.dna$TF %in% names(v.regulators))
v.tfs <- unique(df.dna.de$TF)

df.dna.de <- subset(df.dna.de, df.dna.de$Target %in% gns.DE)
unique(df.dna.de$TF)
length(unique(df.dna.de$Target))

