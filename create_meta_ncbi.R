library(dplyr)



#1. ncbi or lin
data_set = "lin"

#2. with_ovlp, no_olvp, pruned, all 
name = "all_no_ovlp" 

##### dir ####

dir_base = "/Users/sareh/Desktop/surface"

dir_data = file.path(dir_base,"results",data_set)

dir = file.path(dir_data, name)
  
dir_plot = file.path(dir, "plot")

dir_meta_in <- file.path (dir_data, "meta_input") 
dir_meta_out <- file.path (dir_data, "meta_output") 

mle_path <- file.path(dir,"all_mle.csv")

##################### META #####################
dir <- "/Users/sareh/Desktop/surface"
dir_meta <- file.path (dir, "meta_info") 
dir_meta_in <- file.path (dir_meta, "meta_input") 
dir_meta_out <- file.path (dir_meta, "meta_output") 

dir_ncbi <- file.path(dir, "raw_ncbi_files")

##### INPUT #######

#### 1. (all.mle) mle  #####
all.mle <- read.csv(mle_path, sep=" ", header = FALSE)
all.mle$V1<- gsub(".clean", "",all.mle$V1)
colnames(all.mle) <- c("name","codon", "mean.dN", "mean.dS", 
                       "pos.sites", "neg.sites", "mean.dN/mean.dS",
                       "pos.sites/codon","neg.sites/codon")
# clean up all mle
all.mle$virus <- gsub("^([^_]*_[^_]*)_.*$", "\\1",all.mle$name)
all.mle$label <- gsub(".noovlp","", all.mle$name)
all.mle$gene <- gsub("^([^_]*_[^_]*)_", "",all.mle$label)
all.mle$virus <- as.factor(all.mle$virus)
all.mle$name <- gsub(".noovlp","",all.mle$name)
all.mle <- all.mle[,c("name","virus","gene","codon","mean.dN/mean.dS","mean.dN","mean.dS","pos.sites","neg.sites","pos.sites/codon","neg.sites/codon")]
colnames(all.mle)[colnames(all.mle) == 'name'] <- 'file_name'
length(unique(all.mle$virus)) #86
length(unique(all.mle$gene)) #680
#View(all.mle)

names(all.mle)[names(all.mle) == 'codon'] <- 'codon_no_ovlp'


# to get the codon length with ovlp 
mle_withovlp_path = file.path(dir_data,"all_with_ovlp","all_mle.csv")

mle_withovlp <- read.csv(mle_withovlp_path, sep=" ", header = FALSE)

mle_withovlp$V1<- gsub(".clean", "",mle_withovlp$V1)
colnames(mle_withovlp) <- c("name","codon_with_ovlp", "mean.dN", "mean.dS", 
                       "pos.sites", "neg.sites", "mean.dN/mean.dS",
                       "pos.sites/codon","neg.sites/codon")

# clean up all mle
mle_withovlp$virus <- gsub("^([^_]*_[^_]*)_.*$", "\\1",mle_withovlp$name)
mle_withovlp$gene <- gsub("^([^_]*_[^_]*)_", "",mle_withovlp$name)
#View(mle_withovlp)
all_codon <- mle_withovlp[,c("gene","codon_with_ovlp")]
#View(all_codon)

################# from NCBI (process_raw_ncbi.R) 
dir_ncbi <- file.path(dir, "raw_ncbi_files")

# 2. (all_ncbi) ncbi reference genome ####

all_ncbi_path <- file.path(dir_ncbi,"seg_nonseg_ncbi.csv") 
  all_ncbi_accns_path <- file.path(dir_ncbi,"fltr_all_ref_accn.csv")
  all_ncbi <- read.csv(all_ncbi_path, row.names = 1) #[139:98]
  all_ncbi$length.nt <- gsub(" ","",all_ncbi$length.nt)
  all_ncbi$length.nt <- as.numeric(all_ncbi$length.nt)
  colnames(all_ncbi)[colnames(all_ncbi) == 'length.nt'] <- 'genome_length.nt'
  names(all_ncbi)[names(all_ncbi) == 'name'] <- 'name_ncbi'
  #View(all_ncbi) #[1] 139   6
  
 
# 3. (fltr.nbr) ncbi nbr genome ####
fltr.nbr_path <- file.path(dir_ncbi, "fltr_nbr.csv")
  fltr.nbr_accns_path <- file.path(dir_ncbi, "fltr_nbr_accns.csv")
  fltr.nbr <- read.csv(fltr.nbr_path) # colnames don't show up properly
  names(fltr.nbr)[names(fltr.nbr) == 'family'] <- 'family_nbr'
  names(fltr.nbr)[names(fltr.nbr) == 'genus'] <- 'genus_nbr'
  names(fltr.nbr)[names(fltr.nbr) == 'name'] <- 'name_nbr'
  #View(fltr.nbr) #[1] 142973     10
  
################# from UNIPORT 

# 4. (uni) uniport download ####  
# raw downloaded tabular data from UNIPORT
raw_uni_path <- file.path(dir_meta_in,"lin_uniprot.tab")
#raw_uni_path <- file.path(dir_meta,"uniprot_downloaded2022.tsv")
  
uni <- read.csv(raw_uni_path, sep = "\t", header=TRUE)

uni <- uni[,c("Name","Entry","Entry.name","Protein.names",
              "Gene.names","Organism","Length","Organism.ID",
              "Taxonomic.lineage..ALL.","Virus.hosts",
              "Gene.ontology..biological.process.",
              "Gene.ontology..cellular.component.",
              "Gene.ontology..GO.","Gene.ontology.IDs")]

uniport_path <- file.path(dir_meta_in,"uni_clean.csv")
# clean up uni
#uni <- uni[,c("gene","Entry","Entry.Name","Protein.names","Gene.Names","Organism","Length","Organism..ID.","Taxonomic.lineage","Virus.hosts","Gene.Ontology..biological.process.","Gene.Ontology..cellular.component.","Gene.Ontology..GO.","Gene.Ontology.IDs")]
#names(uni)[names(uni) == "From"] <- "gene"
names(uni)[names(uni) == "Name"] <- "gene"
names(uni)[names(uni) == 'Length'] <- 'length_uniport'
#View(uni) # [1] 1060   14

        # sep taxonomy columns 
        library("stringr") 
        taxonomy <- c("virus","realm","kingdom","phylum","class",
                      "order","family","genus","species") # 8 ranks
        tax <- str_split_fixed(uni$Taxonomic.lineage, ",", 13)
        tax<- as.data.frame(tax)
        taxonomy <- c("virus","realm","kingdom","phylum","class",
                      "order","family","genus","species","V1","V2","V3","V4")
        colnames(tax) <- taxonomy
        uni.tax <- cbind(uni,tax)
        #View(uni.tax) ### NEEDS MANUAL CURRATION ###
      
      for_python <- uni.tax[c("gene","virus","realm","kingdom","phylum",
                              "class","order","family","genus","species",
                              "V1","V2","V3","V4")]
      
      write.table(for_python, file.path(dir_meta_in,"for_python_uni_taxonomy.csv"), quote = FALSE, sep = ",",row.names = FALSE)

# 5. (edited_uni_tax) taxonomy ####
# edit it in virus_info_table.py -> python_edited_uni_tax.csv
python_edited_uni_tax_path <- file.path(dir_meta_in, "pyhton_edited_uni_tax.csv")
edited_uni_tax <- read.csv(python_edited_uni_tax_path)
#View(edited_uni_tax)

names(edited_uni_tax)[names(edited_uni_tax) == 'species'] <- 'species_uniport'
names(edited_uni_tax)[names(edited_uni_tax) == 'genus'] <- 'genus_uniport'
names(edited_uni_tax)[names(edited_uni_tax) == 'family'] <- 'family_uniport'
names(edited_uni_tax)[names(edited_uni_tax) == 'order'] <- 'order_uniport'
names(edited_uni_tax)[names(edited_uni_tax) == 'class'] <- 'class_uniport'
names(edited_uni_tax)[names(edited_uni_tax) == 'phylum'] <- 'phylum_uniport'
names(edited_uni_tax)[names(edited_uni_tax) == 'kingdom'] <- 'kingdom_uniport'


#### 6. (tmbed) TMbed predict #### 
tmbed_path <- file.path(dir_meta_in, "tmbed_analysis.csv")
  a <- c("virus","gene","tmbed_len","tmbed_transmem","tmbed_non-transmem", "tmbed_ratio_of_transmem")
  tmbed = read.table(tmbed_path, col.names = a, sep = ",")
  #View(tmbed) 

#### 7. (manual_virus_info) curated virus info #### 
# manually curated database of env/non-env viruses
manual_virus_info_path <- file.path(dir_meta_in,"manual_info.tsv")  
manual_virus_info <- read.csv(manual_virus_info_path, sep = "\t")
names(manual_virus_info)[names(manual_virus_info) == 'name'] <- 'name_manual'
names(manual_virus_info)[names(manual_virus_info) == 'env'] <- 'enveloped'
#View(manual_virus_info)  

#### 8. (manual_gene_non.env) curated SE non-env genes ####
manual_gene_non.env_path <- file.path(dir_meta_in,"manual_non_env_geneSE.csv")
manual_gene_non.env <- read.csv(manual_gene_non.env_path)
#View(manual_gene_non.env)
  # clean it 
  manual_gene_non.env$SE <- gsub(" ", "", manual_gene_non.env$SE)
  manual_gene_non.env$SE <- as.logical(manual_gene_non.env$SE)
  names(manual_gene_non.env)[names(manual_gene_non.env) == 'manual_name'] <- 'manual_gene_name'
 
  table(manual_gene_non.env$SE, useNA = "always")

#### 9. manual_family_info ####  
manual_family_dna_path <- file.path(dir_meta_in,"manual_family_dna.tsv")
manual_family_dna <- read.csv(manual_family_dna_path)

manual_family_dna$DNA <- gsub(" ", "", manual_family_dna$DNA)
table(manual_family_dna$DNA, useNA = "always")
View(manual_family_dna)

# manual_family for each virus name 
manual_family_path <- file.path(dir_meta_in,"family_manual.csv")
manual_family<- read.csv(manual_family_path)

names(manual_family)[names(manual_family) == 'family'] <- 'manual_family'
manual_family$manual_family <- gsub(" ", "", manual_family$manual_family)


# LIN: all virus_info



# Confounders #####

###### 8. aa entropy #####
###### genome 
aa_entropy_genome_path <- file.path(dir_meta_in, "aa_entropy_genome.csv")
aa_entropy_genome <- read.table(aa_entropy_genome_path, sep = " ")
colnames(aa_entropy_genome) <- c("virus","aa_entropy_genome")
#View(aa_entropy_genome)

###### gene
aa_entropy_gene_path <- file.path(dir_meta_in, "aa_entropy_gene.csv")
aa_entropy_gene <- read.table(aa_entropy_gene_path, sep = " ")
  # clean 
  aa_entropy_gene$V1 <- gsub(".fa", "",aa_entropy_gene$V1)
  aa_entropy_gene$V1 <- gsub(".noovlp", "",aa_entropy_gene$V1)
  aa_entropy_gene$virus <- gsub("^([^_]*_[^_]*)_.*$", "\\1",aa_entropy_gene$V1)
  aa_entropy_gene$gene <- gsub("^([^_]*_[^_]*)_", "",aa_entropy_gene$V1)
  colnames(aa_entropy_gene) <- c("file_name","aa_entropy_gene","virus","gene")
  aa_entropy_gene <- aa_entropy_gene[,c("gene","virus","aa_entropy_gene")]

#View(aa_entropy_gene)

###### 9. tree length #####
###### genome #####
tr_genome_path <- file.path(dir_meta_in, "treelength_genome.csv")
tr_genome <- read.csv(tr_genome_path)

tr_genome <- tr_genome[,c("virus_accn","genome_tree_length")]
colnames(tr_genome) <- c("virus","genome_tree_length")
#View(tr_genome)

###### gene #####
tr_gene_path <- file.path(dir_meta_in, "treelength_gene_no_ovlp.csv")
tr_gene <- read.csv(tr_gene_path)
# split
tr_gene$file_name <- gsub(".noovlp", "",tr_gene$file_name)

tr_gene$virus <- gsub("^([^_]*_[^_]*)_.*$", "\\1",tr_gene$file_name)
tr_gene$gene <- gsub("^([^_]*_[^_]*)_", "",tr_gene$file_name)

tr_gene <- tr_gene[,c("gene","nodes","tree_length")]

#View(tr_gene)

# MERGE ######## 

###### virus info ####

virus_info_a <- merge(all_ncbi, manual_virus_info, all.x = TRUE, by.x = "Accession", by.y = "virus" )  
names(virus_info_a)[names(virus_info_a) == 'Accession'] <- 'virus'
#View(virus_info_a)  

# 3. tr_genome (9)
virus_info_b <- merge(virus_info_a, tr_genome, by = "virus", all.x = TRUE)
 
# 4. aa_entropy_genome (8a)
virus_info_c <- merge(virus_info_b, aa_entropy_genome, by = "virus", all.x = TRUE)
#View(virus_info_c)

viruses <- unique(all.mle$virus)
virus_info86 <- virus_info_c[virus_info_c$virus %in% viruses,]

# final version
virus_info <- virus_info86
#View(virus_info) 

write.csv(virus_info, file.path(dir_meta_out,"virus_info.csv"), quote = FALSE)

###### Gene info ####

# 4. build on all.mle 
meta1 <- merge(all.mle, virus_info, by = c('virus'),all.x = TRUE )
names(meta1)[names(meta1) == 'name'] <- 'name_manual'
#View(meta1)

# 5. add GO (UNIPORT) (4)
meta2prep <- merge(meta1, uni, by = "gene", all.x = TRUE)
#View(meta2prep) 

# 6. add edited uni taxonomy data [edited_uni_tax] (5)
meta3prep <- merge(meta2prep, edited_uni_tax, by.x = "gene", by.y = "gene", all.x = TRUE)
#View(meta3prep) 

meta3 <- meta3prep[!duplicated(meta3prep$gene),]   
#View(meta3)

      # GO ---- databases (go description)
      library(AnnotationDbi)
      library(GO.db)
      GO.descrpt <- AnnotationDbi::select(GO.db, keys(GO.db, "GOID"), c("DEFINITION","TERM", "ONTOLOGY"))
      # frequency GO ID
      library("stringr")
      GO.IDs <- meta3$Gene.Ontology.IDs # get all the GOID's of genes studies
      GO.IDs.df <- data.frame(GO.IDs)
      all.id.df <- str_split_fixed(GO.IDs.df$GO.IDs, ";", 70)
      all.id <- as.vector(t(all.id.df))
      all.id <- all.id[all.id != ""] 
      all.id <- gsub(" ","",all.id)
      GO.freq <- as.data.frame(table(as.factor(all.id))) #frequency of GOID's 
      # merge GO Description & GO freq
      GOID.freq <- merge(GO.freq,GO.descrpt,by.x="Var1", by.y="GOID")
      goid_path <- file.path(dir_meta_out, "GOID.freq.csv" )
      write.csv(GOID.freq,goid_path)

# 7. add all_codon
meta3a <- merge(meta3, all_codon, by = "gene")     
#View(meta3a)      
      
# 7. add tmbed (meta3 & tmbed) (6)   
meta4 <- merge(meta3a, tmbed, all.x = TRUE)           
#View(meta4) 

meta4$tmbed_all <- meta4$tmbed_ratio_of_transmem !=0
# meta4 <- meta4 %>% mutate(surface_tmbed_prep = if_else(tmbed_ratio_of_transmem != 0 , TRUE, FALSE))

# 8. add fltr.nbr (3) 
meta.a <- merge(meta4, fltr.nbr, by.x = "virus", by.y = "Representative", all.x = TRUE )
#View(meta.a)

meta.b <- meta.a[!duplicated(meta.a$gene),]
#View(meta.b)

# 9. add treelength 
meta.c <- merge(meta.b, tr_gene, by = "gene", all.x = TRUE)
#View(meta.c)

# 10. add aa entropy 
meta.d <- merge(meta.c, aa_entropy_gene, by = "gene", all.x = TRUE)
#View(meta.d)

# 11. add virus family name (as name)
meta.e <- merge(meta.d, manual_family, by = "name_ncbi", all.x = TRUE)

# 12. add virus family 
# meta.d$family.y + manual_family_dna$Family
meta.f <- merge(meta.e, manual_family_dna, by.x = "manual_family", by.y = "Family", all.x = TRUE )




meta5 <- meta.f

#### split env & non-env 

## Split ####
non_env <- meta5[meta5$enveloped == "FALSE",]
env <- meta5[meta5$enveloped == "TRUE",]
#View(env) #521 #514
#View(non_env) #159

###### surface id ####
    library(data.table)
    
    ### binary surfaces 
    env_surface.id <- c("GO:0055036","GO:0019031")
    # GO:0055036: virion membrane
    # GO:0019031: viral envelope
    non_env_surface.id <- c("GO:0039624","GO:0019028","GO:0039615")
    # GO:0019028: viral capsid 
    # GO:0039615: icosahedral viral capsid

# add surface_go 
    # env
    env$SEgo <- grepl(paste(env_surface.id, collapse = "|"), 
                            env$Gene.Ontology.IDs)
    # non-env
    non_env$SEgo <- grepl(paste(non_env_surface.id, collapse = "|"), 
                                non_env$Gene.Ontology.IDs)

# add tmbed 
     
    # for env viruses if trasmembrane region exsists then its SE 
    env$SEtmbed <- env$tmbed_ratio_of_transmem != 0
    
    non_env$SEtmbed <- NA

    
# env 
    # if no GOID then surface column would be NA and not FALSE 
    env$SEgo[env$Gene.Ontology.IDs == "" ] <- NA
    env$SEgo[env$Gene.Ontology.IDs == NA ] <- NA
    env$SEgo[env$Gene.Ontology.IDs == "NA" ] <- NA
    env$SEgo[is.na(env$Gene.Ontology.IDs) ] <- NA

    table(env$SEgo, env$SEtmbed, useNA = "always")
    
    #View(env)
   
    # half_manual_env_SE.csv 
    # the missing 44 were filled in manually the rest was filled using the code: 
          # if GOID is available then use that! 
          #env$SE <- env$SEgo 
          # missig GOID (NA) [197]
          #table(env$SEgo, useNA = "always") 
          # use TMbed predictions 
          #env$SE[is.na(env$SEgo)] <- env$SEtmbed[is.na(env$SEgo)]
   
      # read it 
      manual_env_SE_path <- file.path(dir_meta_in, "half_manual_env_SE.csv") 
      manual_env_SE <- read.csv(manual_env_SE_path)
      #View(manual_env_SE)
      # merge 
      env <- merge(env, manual_env_SE, by = "gene")
      env$SE <- gsub(" ","", env$SE)
      
    # env is not manually currated 
    env$manual_gene_name <- NA 
    
    table(env$SE, useNA = "always")
    dim(env)
    #View(env)

# non-env
    # if no GOID then surface column would be NA and not FALSE 
    non_env$SEgo[non_env$Gene.Ontology.IDs == "" ] <- NA
    non_env$SEgo[non_env$Gene.Ontology.IDs == NA ] <- NA
    non_env$SEgo[non_env$Gene.Ontology.IDs == "NA" ] <- NA
    non_env$SEgo[is.na(non_env$Gene.Ontology.IDs) ] <- NA

    # add the manual SE (SE) column and gene name (manual_name)   
    non_env <- merge(non_env, manual_gene_non.env, by = "gene", all.x = TRUE)
    dim(non_env)
    #View(non_env)
        
table(env$SE, useNA = "always")
table(non_env$SE, useNA = "always")

# FINAL VERSION
# combine env & non-env 
meta <- rbind(env, non_env) #678
dim(meta) #[1] 678  66
###### add ####
# add a column of sum of pos and neg sites 
meta$pos_neg_sites <- meta$pos.sites + meta$neg.sites

meta$per_pos_neg_sites <- meta$pos_neg_sites / meta$codon_no_ovlp

meta$len_ovlp <- meta$codon_with_ovlp - meta$codon_no_ovlp

###### clean ####
1. YP_009505610.1
2. YP_009505617.1
3. "NP_054708.1"
4. "NP_777383.1"

# remove 2 na (poly genes without mattpeptide entries)
meta <- meta[!is.na(meta$SE),] # 676

# remove two genes with ref alignment off (I don't know how to fix it)
meta <- meta[!(meta$gene == "NP_054708.1" | meta$gene == "NP_777383.1"),]

# remove duplicate rows 
meta <- within(meta, rm("virus.y"))
names(meta)[names(meta) == 'virus.x'] <- 'virus'

# combine family_nbr and family_uniport
meta$family <- ifelse(is.na(meta$family_nbr), 
                      meta$family_uniport, meta$family_nbr)
# combine genus_nbr and genus_uniport 
meta$genus <- ifelse(is.na(meta$genus_nbr), 
                     meta$genus_uniport, meta$genus_nbr)
#a <- meta[,c("family_nbr","family_uniport", "family","manual_family")]

View(meta)

##### SAVE META ####
# ncbi 
meta_path <- "/Users/sareh/Desktop/surface/results/ncbi/meta.csv"

#lin 
# meta_path <- "/Users/sareh/Desktop/surface/results/lin/meta.csv"

#meta_path <- file.path(dir_meta_out,"meta.csv")
write.table(meta, meta_path, quote = FALSE, sep = "\t", row.names = FALSE)














#colnames(meta) <- c("gene", "virus", "file_name", "no_ovlp_codon", 
"meandN.meandS", "mean.dN","mean.dS,pos.sites", 
"neg.sites","pos.sites.codon","neg.sites.codon", 
"virus_name","genome_length.nt", "protein","neighbors",
"segment", "enveloped","genome_tree_length",
"aa_entropy_genome","Entry","Entry.Name","Protein.names",
"Gene.Names","Organism","Length.uniport","Organism.ID",
"Taxonomic.lineage","Virus.hosts","GO.biological.process",
"GO.cellular.component","GO","GOID", "uniport_species",
"uniport.genus","uniport.family","uniport.order",
"uniport.class","uniport.phylum", "uniport.kingdom",
"all_codon","tmbed_len","tmbed_transmem",
"tmbed_non_transmem","tmbed_ratio_of_transmem","tmbed_all",
"Representative","Neighbor","Host","Selected.lineage",
"Taxonomy.name","Segment.name","family","genus",
"name.ncbi", "file_name.y","nodes","gene_tree_length", 
"virus.a","virus.b","aa_entropy_gene","SEgo",
"SEtmbed","SE","manual_name","pos_neg_sites",
"per_pos_neg_sites","len_ovlp")  



#meta <- meta[,c("gene", "virus", "file_name","genome_length.nt", "all_codon",
"no_ovlp_codon","len_ovlp", 
"meandN.meandS", "mean.dN","mean.dS,pos.sites", 
"neg.sites","pos.sites.codon","neg.sites.codon",
"pos_neg_sites", "per_pos_neg_sites", 
"virus_name", "protein","neighbors",
"segment", "enveloped","genome_tree_length",
"aa_entropy_genome","Entry","Entry.Name","Protein.names",
"Gene.Names","Organism","Length.uniport","Organism.ID",
"Taxonomic.lineage","Virus.hosts","GO.biological.process",
"GO.cellular.component","GO","GOID",
"tmbed_len","tmbed_transmem",
"tmbed_non_transmem","tmbed_ratio_of_transmem","tmbed_all",
"Representative","Neighbor","Host","Selected.lineage",
"Taxonomy.name","Segment.name","family","genus",
"name.ncbi","nodes","gene_tree_length","aa_entropy_gene","SEgo",
"SEtmbed","SE","manual_name")]




###### optional 
# replace this with the Neighbor column 
# need to write all neighbours into one cell as a list
all_ncbi <- aggregate(Neighbor~Accession, data = fltr.nbr, FUN = 'list')
meta.nbr <- merge(meta, all_ncbi, by.x = "virus", by.y = "Accession"  ,all.x = TRUE)

saveRDS(dir_meta,file="meta.nbr.Rda")












#### SUMMARY ####
table(meta$enveloped, useNA = "always") # 521 env 
table(meta$SE, useNA = "always") # 138 SE (2 NA)

table(meta$SE, meta$enveloped)
prop.table(table(meta$SE, meta$enveloped))

# table for venn diagram 
env_se <- meta[,c("SE","enveloped")]
pdf(file=file.path(dir,"plot","env_se_stackedbarplot.pdf"))
x <- barplot(table(env_se), 
        col=c(nonSE,SE), 
        border="white", 
        space=0.04, 
        font.axis=2, 
        xlab = "Enveloped virus")

legend("topleft",
       fill=c(SE,nonSE),
       legend=c("SE", "non-SE"))
dev.off()

#https://stackoverflow.com/questions/40919759/stacked-barplot-using-r-base-how-to-add-values-inside-each-stacked-bar



###### virus family #####

# frequency of GENES
meta$family.y <- as.factor(meta$family.y)
gene_family <- table(meta$family.y)
gene_df <- as.data.frame(gene_family)
colnames(gene_df) <- c("family","gene_freq")

# frequency of VIRUSES 
virus_table$virus_family <- as.factor(virus_table$virus_family)
virus_family <- table(virus_table$virus_family)
virus_df <- as.data.frame(virus_family)
colnames(virus_df) <- c("family","virus_freq")

pie_freq <- merge(gene_df,virus_df, by = "family")

# pie_freq_long
colnames(gene_df) <- c("family","freq")
gene_df$name <- "gene"

colnames(virus_df) <- c("family","freq")
virus_df$name <- "virus"

pie_freq_long <- rbind(gene_df,virus_df)
View(pie_freq_long)


######## pie chart #######
# as is 
png(file=file.path(dir,"virus_family_gene_asis.png"),width=800, height=800)
pie(gene_family)
dev.off()

png(file=file.path(dir,"virus_family_virus_asis.png"),width=800, height=800)
pie(virus_family)
dev.off()

# distinctive colors
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
color = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

### (1) frequency of genes from each family 
png(file=file.path(dir,"piechart_gene_family.png"),
    width=800, height=800)
pie(gene_family, cex=0.7, border="grey", main="Virus Families", col=color, cex.main=1.7)
dev.off()

png(file=file.path(dir,"piechart_virus_family.png"),
    width=800, height=800)
pie(virus_family, cex=0.7, border="grey", main="Virus Families", col=color, cex.main=1.7)
dev.off()

###### nested pie chart ####
#https://biostats.w.uib.no/creating-a-multiple-pie-donut-chart/
library(tidyverse)
ggplot(pie_freq, aes(x=virus_freq, y=gene_freq, fill=family))


#https://plotly.com/r/pie-charts/
library(plotly)
library(dplyr)

# gene
fig <- pie_freq %>% plot_ly(labels = ~family, values = ~gene_freq)
fig <- fig %>% add_pie(hole = 0.6)
fig <- fig %>% layout(title = "Donut charts using Plotly",  showlegend = F,
                      xaxis = list(showgrid = FALSE, zeroline = FALSE,
                                   showticklabels = FALSE),
                      yaxis = list(showgrid = FALSE, zeroline = FALSE,
                                   showticklabels = FALSE))
fig

# virus
fig <- pie_freq %>% plot_ly(labels = ~family, values = ~virus_freq)
fig <- fig %>% add_pie(hole = 0.6)
fig <- fig %>% layout(title = "Donut charts using Plotly",  showlegend = F,
                      xaxis = list(showgrid = FALSE, zeroline = FALSE,
                                   showticklabels = FALSE),
                      yaxis = list(showgrid = FALSE, zeroline = FALSE,
                                   showticklabels = FALSE))
fig

# https://biostats.w.uib.no/creating-a-multiple-pie-donut-chart/
pie_freq_long$name <- as.factor(pie_freq_long$name)
str(pie_freq_long)

# create pie chart using plot_ly() function
pdf(file=file.path(dir,"nested_pie_virus.pdf"))
# create pie chart using plot_ly() function

ggplot(pie_freq_long, aes(x = name, y = freq, fill = family)) +
  geom_col() +
  scale_fill_viridis_d() +
  coord_polar("y")

dev.off()


###### VIRUS INFO TABLE #####
# name.x & name.y are different 
meta.virus <- meta[!duplicated(meta$name.y),] # 56 
# should be less than ref_accns of 86 because [length(unique(all.mle$virus))] because some are segmented viruses  

b <- meta.virus[,c('virus','name.y','enveloped',"genome_length.nt" ,"Organism", "Taxonomic.lineage","family.y","genus.y","family.x","class" )]
View(b)

all_virus_info <- b

### table for thesis 
virus_table <- all_virus_info[,c("virus","genome_length.nt",'name.y',"family.y")]
colnames(virus_table) <- c("Accession","genome_length","virus_name","virus_family")
View(virus_table)

##### SAVE virus_table
virus_tbl_path <- file.path(dir,"virus_table.csv")
write.table(virus_table, virus_tbl_path, quote = FALSE, sep = ",", row.names = FALSE)
View(virus_table) 


###### PLOTS ####
pdf(file=file.path(dir,"plot","hist_ovlp_len.pdf"))
hist(meta$len_ovlp, breaks = 30,
     xlab = "length of overlapping genes",
     col = "cadetblue")
dev.off()

pdf(file=file.path(dir,"plot","hist_ovlp_len_log.pdf"))
hist(log(meta$len_ovlp,10), breaks = 20, 
     xlab = "log transformed length of overlapping genes",
     col = "cadetblue")
dev.off()










##### EXAMINE ##### 

#dN/dS scores 
summary(meta$mean.dN.mean.dS)
summary(meta1$mean.dN.mean.dS)
summary(meta2$mean.dN.mean.dS)
summary(meta3$mean.dN.mean.dS)
summary(meta4$mean.dN.mean.dS)
summary(meta5$mean.dN.mean.dS)

#SE
table(env$SEgo,useNA = "always")
table(env$SEtmbed,useNA = "always")
table(env$SEgo, env$SEtmbed, useNA = "always")

table(meta$SEgo, useNA = "always")
table(meta$SEtmbed, useNA = "always")
table(meta$SEtmbed, meta$SEgo, useNA = "always")


# env 
env <- meta[meta$enveloped == "TRUE",]
View(env)

a <- env[,c("virus","gene","SE","SEgo","SEtmbed","Protein.names","Gene.Names")]
View(a)

write.csv(a,"/Users/sareh/Desktop/env", row.names = FALSE, quote = FALSE)

# non_env
non_env <- meta[meta$enveloped == "FALSE",]
View(non_env)

b <- non_env[,c("virus","gene","SE","SEgo","SEtmbed","Protein.names","Gene.Names","manual_name")]
View(b)

write.csv(b,"/Users/sareh/Desktop/non_env",row.names = FALSE, quote = FALSE)


path_old_meta <- "/Users/sareh/Desktop/surface/meta.csv"
old_meta <- read.csv(path_old_meta,sep = "\t")
View(old_meta)

# df
len <- meta3a[,c("codon","genome_length.nt","Length","all_codon")]
View(len)





## GO databse ####
# GO ---- databases (go description)
library(AnnotationDbi)
library(GO.db)
GO.descrpt = select(GO.db, keys(GO.db, "GOID"), c("DEFINITION","TERM", "ONTOLOGY"))
# frequency GO ID
library("stringr")
GO.IDs <- meta$Gene.Ontology.IDs # get all the GOID's of genes studies
GO.IDs.df <- data.frame(GO.IDs)
all.id.df <- str_split_fixed(GO.IDs.df$GO.IDs, ";", 70)
all.id <- as.vector(t(all.id.df))
all.id <- all.id[all.id != ""] 
all.id <- gsub(" ","",all.id)
GO.freq <- as.data.frame(table(as.factor(all.id))) #frequency of GOID's 
# merge GO Description & GO freq
GOID.freq <- merge(GO.freq,GO.descrpt,by.x="Var1", by.y="GOID")
#write.csv(GOID.freq,"GOID.freq.csv")

# top 10 ordered 
GOID.freq.clean <- (GOID.freq[GOID.freq$Freq > 30,])
png(file= file.path(plot_dir,"high_frequency_GOID.pdf"))
barplot(GOID.freq.clean[order(GOID.freq.clean[,2],decreasing=TRUE),][,2],names.arg=GOID.freq.clean[order(GOID.freq.clean[,2],decreasing=TRUE),][,4], col="cornflowerblue",cex.names = 0.5,las=2,main="Frequency of GO IDs")
dev.off()












