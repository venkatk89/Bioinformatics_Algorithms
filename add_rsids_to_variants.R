#load necessary libraries
library(data.table) 

#read base data and add columns
base_data <- read.table("QC.GWAS.txt.gz", header = T) 
colnames(base_data) <- c("SNP", "Chr", "Pos", "EA", "NEA", "EAF", "Beta", "SE", "Pvalue", "Neff")
 
#read target data and add columns
dbsnp_hg19 <- fread("/data/dbsnp/hg19/00-All.vcf.gz") 
colnames(dbsnp_hg19) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")

#loop through the base data and replace ids with rsids
i <- 1 
for(i in c(1:nrow(base_data))){
  chr <- base_data$Chr[i]
  pos <- base_data$Pos[i]
  id <- dbsnp_hg19$ID[which(dbsnp_hg19$CHROM == chr & dbsnp_hg19$POS == pos)]
  # if id is not null, substitute it in base data
  if(length(id)){
      base_data$SNP[i] <- id
  }
  # print log for every 100K variants processed 
  if(i %% 100000 == 0){
      print(paste(i,"variants processed"))
  }
}

# write rsid added GWAS data to a txt file
write.table(base_data, 
  "QC.RSID.GWAS.txt", 
  append = FALSE,
  quote = F,
  row.names = F,
  sep="\t")
