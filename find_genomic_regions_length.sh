# Code to extract the sizes of different functional regions from Gencode GTF file. 
# In this particular example we work with "gencode.v43.basic.annotation.gtf" downloaded from https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.basic.annotation.gtf.gz 



# Divide gencode annotations into coding and non-coding genes

# Get non-coding regions' annotation. (Taken from https://www.gencodegenes.org/pages/biotypes.html)
grep \
-e "Mt_rRNA" \
-e "Mt_tRNA" \
-e "miRNA" \
-e "misc_RNA" \
-e "rRNA" \
-e "scRNA" \
-e "snRNA" \
-e "snoRNA" \
-e "ribozyme" \
-e "sRNA" \
-e "scaRNA" \
-e "lncRNA" \
-e "Mt_tRNA_pseudogene" \
-e "tRNA_pseudogene" \
-e "snoRNA_pseudogene" \
-e "snRNA_pseudogene" \
-e "scRNA_pseudogene" \
-e "rRNA_pseudogene" \
-e "misc_RNA_pseudogene" \
-e "miRNA_pseudogene" \
gencode.v43.basic.annotation.gtf > gencode_basic_non_coding.txt

# Get coding regions' annotation.
grep -v -E "Mt_rRNA|Mt_tRNA|miRNA|misc_RNA|rRNA|scRNA|snRNA|snoRNA|ribozyme|sRNA|scaRNA|lncRNA|Mt_tRNA_pseudogene|tRNA_pseudogene|snoRNA_pseudogene|snRNA_pseudogene|scRNA_pseudogene|rRNA_pseudogene|misc_RNA_pseudogene|miRNA_pseudogene" gencode.v43.basic.annotation.gtf > gencode_basic_coding.txt



# make a bed file for chromomosome sizes: Can be taken from UCSC genome browser. For HG38, content present at the end of this file
awk 'OFS="\t" {print $1, "0", $2}' chromSizes.txt | sort -k1,1 -k2,2n > chromSizes.bed

#sort chromSizes.txt file
cat chromsizes.txt | awk 'OFS="\t" {print $1, $2}' | sort -k1,1 -k2,2n > sorted_chromSizes.txt

# get sum all region lengths in chromSizes.bed
awk -F'\t' '{sum+=($3 - $2 );}END{print sum;}' chromSizes.bed # 3,088,286,401



#sort the annotations in coding and non-coding gtf files by chr, start, end
cat gencode_basic_coding.txt | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k4,4n -k5,5n"}' > gencode_basic_coding_sorted.txt

cat gencode_basic_non_coding.txt | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k4,4n -k5,5n"}' > gencode_basic_non_coding_sorted.txt



# # Work with coding genes first

# extract all regions listed as genic
grep -v "^#" gencode_basic_coding_sorted.txt | awk -F '\t' -v OFS='\t' '$3 == "gene" {print $1, $4-1, $5}'  > coding_gene_all.bed
# merge overlapping genic regions
bedtools merge -i coding_gene_all.bed > coding_gene_sorted.bed

# get intergenic regions as {all_regions - genic_regions}
bedtools complement -i coding_gene_sorted.bed -g ../chrom_sizes_files/sorted_chromSizes.txt > coding_intergenic_sorted.bed

# extract all regions listed as exonic
grep -v "^#" gencode_basic_coding_sorted.txt | awk -F '\t' -v OFS='\t' '$3 == "exon" {print $1, $4-1, $5}'  > coding_exon_all.bed
# merge overlapping exonic regions
bedtools merge -i coding_exon_all.bed > coding_exon_sorted.bed

# get intronic regions as {all_regions - (exons + intergenic)} 
bedtools complement -i <(cat coding_exon_sorted.bed coding_intergenic_sorted.bed | sort -k1,1 -k2,2n) -g ../chrom_sizes_files/sorted_chromSizes.txt > coding_intron_sorted.bed

# get sum all region lengths in exon_sorted.bed
awk -F'\t' '{sum+=($3 - $2 + 1);}END{print sum;}' coding_exon_all.bed # 211,126,192
awk -F'\t' '{sum+=($3 - $2 + 1);}END{print sum;}' coding_exon_sorted.bed # 94,124,660
# get sum all region lengths in genic_sorted.bed
awk -F'\t' '{sum+=($3 - $2 + 1);}END{print sum;}' coding_gene_sorted.bed # 1,351,498,292
# get sum all region lengths in intergenic_sorted.bed
awk -F'\t' '{sum+=($3 - $2 + 1);}END{print sum;}' coding_intergenic_sorted.bed # 1,736,841,784
# get sum all region lengths in intron_sorted.bed
awk -F'\t' '{sum+=($3 - $2 );}END{print sum;}' coding_intron_sorted.bed # 1,257,593,381

# Length of splice sites: No. of Exons*4
wc -l ../coding_genes/coding_exon_sorted.bed # 246,574
# 246574*4 = 986,296

# So lenth of Intron is (Total Intron - Splice sites) ==> (1,257,593,381 - 986,296) = 1,256,607,085




# # Work with non-coding genes second

# extract all regions listed as genic
grep -v "^#" gencode_basic_non_coding_sorted.txt | awk -F '\t' -v OFS='\t' '$3 == "gene" {print $1, $4-1, $5}'  > non_coding_gene_all.bed
# merge overlapping genic regions
bedtools merge -i non_coding_gene_all.bed > non_coding_gene_sorted.bed

# get intergenic regions as {all_regions - genic_regions}
bedtools complement -i non_coding_gene_sorted.bed -g ../chrom_sizes_files/sorted_chromSizes.txt > non_coding_intergenic_sorted.bed

# extract all regions listed as exonic
grep -v "^#" gencode_basic_non_coding_sorted.txt | awk -F '\t' -v OFS='\t' '$3 == "exon" {print $1, $4-1, $5}'  > non_coding_exon_all.bed
# merge overlapping exonic regions
bedtools merge -i non_coding_exon_all.bed > non_coding_exon_sorted.bed


# Length of non-coding splice sites: No. of non-coding Exons*4
wc -l ../non_coding_genes/non_coding_exon_sorted.bed # 73,579
# 73579*4 = 294,316

# get intronic regions as {all_regions - (exons + intergenic)} 
bedtools complement -i <(cat non_coding_exon_sorted.bed non_coding_intergenic_sorted.bed | sort -k1,1 -k2,2n) -g ../chrom_sizes_files/sorted_chromSizes.txt > non_coding_intron_sorted.bed

# get sum all region lengths in exon_sorted.bed
awk -F'\t' '{sum+=($3 - $2 + 1);}END{print sum;}' non_coding_exon_all.bed # 42,487,318
awk -F'\t' '{sum+=($3 - $2 + 1);}END{print sum;}' non_coding_exon_sorted.bed # 34,169,545
# get sum all region lengths in genic_sorted.bed
awk -F'\t' '{sum+=($3 - $2 + 1);}END{print sum;}' non_coding_gene_sorted.bed # 545,287,899
# get sum all region lengths in intergenic_sorted.bed
awk -F'\t' '{sum+=($3 - $2 + 1);}END{print sum;}' non_coding_intergenic_sorted.bed # 2,543,041,645
# get sum all region lengths in intron_sorted.bed
awk -F'\t' '{sum+=($3 - $2 );}END{print sum;}' non_coding_intron_sorted.bed # 511,170,374

# So lenth of nc_Intron is (Total nc_Intron - nc_Splice sites) ==> (511,170,374 - 294,316) = 510,876,058



# Get all upstream regions
grep -v "^#" ../gencode_sorted.gtf | awk -F '\t' -v OFS='\t' '$3 == "transcript" {print $1, $4-1000, $4}'  > upstream_all.bed
grep -v "chrM" upstream_all.bed > upstream_edited.bed
# merge overlapping upstream regions
bedtools merge -i upstream_edited.bed > upstream_sorted.bed
# get sum all region lengths in upstream_sorted.bed
awk -F'\t' '{sum+=($3 - $2 );}END{print sum;}' upstream_sorted.bed # 77,207,458



# Get all downstream regions
grep -v "^#" ../gencode_sorted.gtf | awk -F '\t' -v OFS='\t' '$3 == "transcript" {print $1, $5, $5 + 1000}'  > downstream_all.bed
# sort file
grep -v "chrM" downstream_all.bed | sort -k1,1 -k2,2n > downstream_edited.bed
# merge overlapping downstream regions
bedtools merge -i downstream_edited.bed > downstream_sorted.bed
# get sum all region lengths in downstream_sorted.bed
awk -F'\t' '{sum+=($3 - $2 );}END{print sum;}' downstream_sorted.bed # 77,525,040



# Get total intergenic sites
bedtools complement -i <(cat ../non_coding_genes/non_coding_gene_sorted.bed ../coding_genes/coding_gene_sorted.bed ../upstream/upstream_sorted.bed ../downstream/downstream_sorted.bed | sort -k1,1 -k2,2n) -g ../chrom_sizes_files/sorted_chromSizes.txt > intergenic_sorted.bed
# get sum all region lengths in intergenic_sorted.bed
awk -F'\t' '{sum+=($3 - $2 );}END{print sum;}' intergenic_sorted.bed # 1,256,141,639



# Get total UTR sites
grep -v "^#" ../gencode_sorted.gtf | awk -F '\t' -v OFS='\t' '$3 == "UTR" {print $1, $4-1, $5}'  > UTR_all.bed
# merge overlapping UTR regions
bedtools merge -i UTR_all.bed > UTR_sorted.bed
# get sum all region lengths in UTR_sorted.bed
awk -F'\t' '{sum+=($3 - $2 );}END{print sum;}' UTR_sorted.bed # 48,115,889

# Total length of CDS is (Total Exon - Total UTR) = (34,169,545 + 94,124,660 - 48,115,889) = 80,178,319



# Summary of lengths calculated for GRCH38
#
# Genomic Region  Length (BPs)    
# Total           3,088,286,401
# Exonic          94,124,660
# Intronic        1,256,607,085
# ncRNA_exonic    34,169,545
# ncRNA_intronic  510,876,058
# Splicing        986,296
# ncRNA_Splicing  294,316
# UTRs            48,115,889
# Upstream        77,207,458
# Downstream      77,525,040
# Intergenic      1,256,141,639



#Contents of chromSized.txt for GRCH38
#
# chr1 248956422
# chr2 242193529
# chr3 198295559
# chr4 190214555
# chr5 181538259
# chr6 170805979
# chr7 159345973
# chr8 145138636
# chr9 138394717
# chr10 133797422
# chr11 135086622
# chr12 133275309
# chr13 114364328
# chr14 107043718
# chr15 101991189
# chr16 90338345
# chr17 83257441
# chr18 80373285
# chr19 58617616
# chr20 64444167
# chr21 46709983
# chr22 50818468
# chrX 156040895
# chrY 57227415
# chrM 16569