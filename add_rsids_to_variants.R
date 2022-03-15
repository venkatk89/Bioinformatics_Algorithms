
#read base data and add columns accordingly
base_data <- read.table("QC.GWAS.txt.gz", header = T) 
colnames(base_data) <- c("SNP", "Chr", "Pos", "EA", "NEA", "EAF", "Beta", "SE", "Pvalue", "Neff")
 
# Preparation of rsid_data
# Use bash to filter only the columns containing chromosome, position, and rsid 
# The following command achieves this for the dbsnp dataset:
## zcat 00-common_all.vcf.gz | cut -f 1,2,3 > chrom_pos_rsid_dbsnp_common 

#read  rsid_data and add columns
dbsnp_hg19 <- read.table("/data/dbsnp/hg19/chrom_pos_rsid_dbsnp_common") 
colnames(dbsnp_hg19) <- c("Chr", "Pos", "ID")

# Merge the two datasets based on chr and pos. This step adds rsids to a new column called ID and removes all SNPs that aren't in dbsnp. 
rsid_added_base_data <- merge(base_data, dbsnp_hg19, by = c("Chr", "Pos"))


# write rsid added GWAS data to a txt file
write.table(rsid_added_base_data, 
  "QC.RSID.GWAS.txt", 
  append = FALSE,
  quote = F,
  row.names = F,
  sep="\t")




# inefficient algorithm, but keeps the ID of non-matching SNPs intact.
#loop through the base data and replace ids with rsids
#i <- 1 
#for(i in c(1:nrow(base_data))){
#  chr <- base_data$Chr[i]
#  pos <- base_data$Pos[i]
#  id <- dbsnp_hg19$ID[which(dbsnp_hg19$CHROM == chr & dbsnp_hg19$POS == pos)]
#  # if id is not null, substitute it in base data
#  if(length(id)){
#      base_data$SNP[i] <- id
#  }
#  # print log for every 100K variants processed 
#  if(i %% 100000 == 0){
#      print(paste(i,"variants processed"))
#  }
#}


