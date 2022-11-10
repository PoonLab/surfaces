# Process Raw Downloaded NCBI tables 

library(dplyr)
library(stringi)
library(stringr)
library(plyr)


dir <- "/Users/sareh/Desktop/surface/raw_ncbi_files"

############################
####### OUTPUT file ########
############################
# opening the files processed in this script

# (1) ncbi reference genome (each row is a segment)
all_ncbi_path <- file.path(dir,"seg_nonseg_ncbi.csv")  # the full table 
all_ncbi <- read.csv(all_ncbi_path, row.names = 1) #[139:98] 
View(all_ncbi)

all_ncbi_accns_path <- file.path(dir,"fltr_all_ref_accn.csv") # only accns
fltr_non_seg.ncbi_path <- file.path(dir, "/fltr_non_seg.ncbi.csv") # tbl only nonseg viruses
fltr_seg.ncbi_path <- file.path(dir,"fltr_seg.ncbi.csv") # tbl only segmented 


# (2) ncbi nbr genome
# old:"/Users/sareh/surfaces/r_surfaces/INDEX/2020/virus_info/fltr_partial_nbr.csv"
fltr.nbr_path <- file.path(dir, "fltr_nbr.csv")
fltr.nbr_accns_path <- file.path(dir, "fltr_nbr_accns.csv")
fltr.nbr <- read.csv(fltr.nbr_path) # colnames don't show up properly
View(fltr.nbr)


###########################
####### INPUT file ####### 
###########################
# file paths needed to run this script 

  #raw ncbi viral genome files filtered for human host 
  #https://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=10239&host=human

  # (1)taxid10239.tbl.csv
    taxid.tbl_path <- "/taxid10239.tbl.csv"
       # (1a) seg_ncbi.csv #old = "/Users/sareh/surfaces/r_surfaces/INDEX/virus_info/segment.tbl"
       taxid.tbl_path3 <- "/seg_ncbi.csv" #processed in bash 
       # (1b) nonseg_ncbi.csv #old = "/Users/sareh/surfaces/r_surfaces/INDEX/virus_info/ncbi.csv"
       taxid.tbl_path2 <- "/nonseg_ncbi.csv" #processed in bash 
       
  # (2)taxid10239.acc_lst
    taxid.acc_lst_path <- "/taxid10239.acc_lst"
    
  # (3)taxid10239.nbr #Delete header of this file manually from downloaded file
    taxid.nbr_path <- "/taxid10239.nbr"
    

### Raw files ###
#raw donwloaded files (2020)

## 1. Complete genomes: Viruses (taxid 10239) filtered by host 'human'
taxid.tbl <- read.csv(file.path(dir,taxid.tbl_path), header = FALSE)
#View(taxid.tbl)

  # NON-SEGMENTED Get rid of the lines starting with tab
  nonseg_ncbi <- read.table(file.path(dir,taxid.tbl_path2), fill = TRUE, sep=",")
  View(nonseg_ncbi) #457
  
  # SEGMENTED only lines starting with tab
  seg_ncbi <- read.csv(file.path(dir,taxid.tbl_path3), fill = TRUE, header = F, stringsAsFactors = TRUE)
  View(seg_ncbi) #389

## 2. list of ref accessions (taxid.acc_lst)
taxid.acc_lst <- read.csv(file.path(dir,taxid.acc_lst_path), header = FALSE) 
      View(taxid.acc_lst)

## 3. Genome neighbours data (taxid.nbr)
taxid.nbr <- read.csv(file.path(dir,taxid.nbr_path),sep="\t",stringsAsFactors = T)
      View(taxid.nbr)

      
##### non-seg ##### (nonseg_ncbi)

# cleaning
colnames(nonseg_ncbi) <- c("Accession","RefSeq.type","source.information","segm","length","protein","Neighbors","host","created","updated","empty")
nonseg_ncbi$Neighbors <- as.numeric(nonseg_ncbi$Neighbors)
nonseg_ncbi$protein <- as.numeric(nonseg_ncbi$protein)
nonseg_ncbi$name <- row.names(nonseg_ncbi)

## subsetting 

###remove segmented viruses 
ncbi.tbl1 <- nonseg_ncbi[nonseg_ncbi$Accession!= "-", ]
      dim(ncbi.tbl1) #[1] 368  11
###remove incomplete seq viruses 
ncbi.tbl2 <- ncbi.tbl1[ncbi.tbl1$RefSeq.type == "complete",]
      dim(ncbi.tbl2) #[1] 365  11
###viruses > 100 neighbours 
ncbi.tbl3 <- ncbi.tbl2[ncbi.tbl2$Neighbors >= 100, ]
      dim(ncbi.tbl3) #[1] 167  11
###remove Neighbors == NA 
ncbi.tbl4 <- ncbi.tbl3[!is.na(ncbi.tbl3$Accession),] #[1] 83  11
      dim(ncbi.tbl4)
###viruses proteins 
ncbi.tbl5 <- ncbi.tbl4[ncbi.tbl4$protein >= 1,  ]
      dim(ncbi.tbl5) #[1] 83 11
ncbi.tbl5$segment <- "FALSE" # add a column

###saving the file (non-empty)
fltr_nonseg_ncbi <- ncbi.tbl5[!is.na(ncbi.tbl5$Accession), ]
      dim(fltr_nonseg_ncbi) #[1] 83 11
      View(fltr_nonseg_ncbi)
fltr_non_seg.ncbi_path <- file.path(dir,"fltr_non_seg.ncbi.csv")
write.csv(fltr_nonseg_ncbi,fltr_non_seg.ncbi_path)
fltr_non_seg.ncbi_accns_path <- file.path(dir,"fltr_non_seg.ncbi_accns.csv")
write.csv(fltr_nonseg_ncbi$Accession,fltr_non_seg.ncbi_accns_path)

##### segmented ##### (seg_ncbi)

#Clean/prepare the file
# remove rows with empty neighbour numbers 
seg <- seg_ncbi[!(is.na(seg_ncbi$V5) | seg_ncbi$V5==""), ]
seg$V6 <-  gsub("neighbors:","",seg$V5)
seg$V7 <-  gsub("proteins:","",seg$V4)
seg$V6 <- as.numeric(seg$V6)
seg$V7 <- as.numeric(seg$V7)
colnames(seg) <- c("name","length","Accn","V4","V5","neighbours","protein")
seg$segment <- "TRUE" # add a column 
      View(seg)
      dim(seg) #[1] 313   8

#Viruses with all segments > 100 Neighbours available 
mean.neibour <- seg%>%group_by(name)%>%summarise(Average=mean(neighbours))
neibour.100 <- mean.neibour[mean.neibour$Average>100,] 
#neibour.100$name #15 viruses
fltr_seg.ncbi <- seg[seg$name %in% neibour.100$name,] 
      dim(seg) #[1] 313   7
      dim(fltr_seg.ncbi) #[1] 56  7
      length(unique(fltr_seg.ncbi$Accn)) #[1] 56 (accns for every segment)
      length(unique(fltr_seg.ncbi$name)) #[1] 15 (virus species)
      View(fltr_seg.ncbi)

# save filtered segmented viruses #table [each row a segment]
# OUTPUT : fltr_seg.ncbi_path  <- file.path(dir,"fltr_seg.ncbi.csv")
write.csv(fltr_seg.ncbi, file.path(dir,"fltr_seg.ncbi.csv"))
# OUTPUT : fltr_seg.ncbi_accns_path <-file.path(dir,"fltr_seg.ncbi_accns.csv")
write.csv(fltr_seg.ncbi$Accn,file.path(dir,"fltr_seg.ncbi_accns.csv"), quote=F, row.names = F)


########## ref accns #########
########## combine ########### **combine seg + non-seg viruses info**

      #View(fltr_seg.ncbi) # segmented: 
      #View(fltr_nonseg_ncbi) # non-segmented: 

fltr_seg.ncbi$host <- "human"
colnames(fltr_seg.ncbi)<-c("name","length.nt","Accession","V4","V5","neighbors","protein","segment" )
fltr_seg.ncbi <- fltr_seg.ncbi[,c("name","Accession","length.nt","protein","neighbors","segment")]
colnames(fltr_nonseg_ncbi)<-c("Accession","RefSeq.type","source","remove","length.nt","protein","neighbors","host","created","updated","empty","name","segment")
fltr_nonseg_ncbi <- fltr_nonseg_ncbi[,c("name","Accession","length.nt","protein","neighbors","segment")]
fltr_nonseg_ncbi$length.nt <- gsub("(nt)","",fltr_nonseg_ncbi$length.nt)
fltr_seg.ncbi$length.nt <- gsub("(nt)","",fltr_seg.ncbi$length.nt)
fltr_seg.ncbi$length.nt <- gsub("[()]", "",fltr_seg.ncbi$length.nt)
all_ncbi<- rbind(fltr_seg.ncbi,fltr_nonseg_ncbi) #virus_info

###################
View(all_ncbi) #139 all fltred ref accns 
###################

# writing it out
write.csv(all_ncbi, all_ncbi_path) # table
# all ref accns in a file 
all_ref_accn <- all_ncbi$Accession
write.csv(all_ncbi$Accession,all_ncbi_accns_path, quote=F, row.names = F)

##############################
########## neighbor ########## taxid.nbr <- taxid10239.nbr
##############################

raw.ncbi.nbr <- read.csv(file.path(dir, taxid.nbr_path),sep="\t",stringsAsFactors = T)
      View(raw.ncbi.nbr)
      dim(raw.ncbi.nbr) #[1] 255844      6
      length(unique(raw.ncbi.nbr$Neighbor)) # [1] 249053

### 1. **hmn.nbr**  | Filter human host 
  #"human" in the host column & remove empty host (unknown host)
  hmn.nbr <- raw.ncbi.nbr[grepl("human", raw.ncbi.nbr$Host, fixed = TRUE),]
        dim(raw.ncbi.nbr) #255844
        dim(hmn.nbr) #160,351
        255844-160351 # 95,493 removed 
        length(unique(hmn.nbr$Representative)) #661 (total viruses)
        length(unique(hmn.nbr$Neighbor)) #154519 (there are some repeat nbrs)
        View(hmn.nbr)

### 2. **no.dup.nbr** | Duplicate Neighbour accn
  # Remove duplicate rows of the dataframe using NAME variable
  no.dup.nbr <- distinct(hmn.nbr,Neighbor, .keep_all= TRUE)
        length(unique(no.dup.nbr$Neighbor)) #[1] 154519 checking all are unique
        dim(hmn.nbr) #[1] 160351      6
        dim(no.dup.nbr) #[1] 154,519      6
        160351 - 154519 # 5,832 removed
        View(no.dup.nbr)

### 3. Pick one (more than one ref assigned to a nbr accn)
#example rep ["NC_002640" "NC_001477" "NC_001474" "NC_001475"] all in fltred ref accns and all have the same nbr assigned to them so no point in running it for each ref accn since EXACTLY same nbr genomes are assinged to them!       
#the entries look very similar cant think of a systematic reasoning for which to pick so I just pick the first item
# [less.10] ONLY the one ref accn
no.dup.nbr$Representative <-as.character(no.dup.nbr$Representative)
      dim(no.dup.nbr) #[1] 154,519      6 
### 1 RefSeq ###
  # if length of characters in rep accn < 10 then only one ref-accns
  less.10 <- no.dup.nbr[nchar(no.dup.nbr$Representative)<=10,]
        dim(less.10) #[1] 133,638      6
        length(unique(less.10$Neighbor)) 
        View(less.10)

### >1 RefSeq ###
# if len of characters in ref accns > 10 then there is more than one accn listed as a representative
  more.10 <- no.dup.nbr[nchar(no.dup.nbr$Representative)>10,]
        dim(more.10) # [1] 20,881     6
        length(unique(more.10$Neighbor)) # 
        View(more.10)

more <- unique(more.10$Representative)
length(more) #57        

# representatives where at least one is in the filtered ref accns 
a <- as.list(more)
accn.morethan1 <- c() # accns in fltred ref accn 
picked.morethan1 <- c() # only one of the accns in fltred ref accn
left <- c() # only ones with more than one accn in fltred ref accns
for (f in a) {
  b <- unlist(strsplit(f,","))
  c <- as.character(stri_remove_empty(b, na_empty = FALSE))
  n <- intersect(c, all_ref_accn)
  # print(length(n))
  accn.morethan1 <- c(accn.morethan1,n)
  # just pick first item of each list 
  picked.morethan1 <- c(picked.morethan1, n[1])
  if (length(n)>1){
    left <- c(left,n)
  }
}

length(accn.morethan1) #54
length(unique(accn.morethan1)) #54

length(left) #50
print(left)

length(picked.morethan1) # 52
length(unique(picked.morethan1))

# remove NA
  picked.morethan1.noNA <- picked.morethan1[!is.na(picked.morethan1)]
  length(picked.morethan1.noNA) #[1] 20
  picked.morethan1.noNA

# picked more than one in fltred ncbi ref accns 
  ncbi.morethan1accn <- all_ncbi[all_ncbi$Accession %in% picked.morethan1.noNA,]
  View(ncbi.morethan1accn) # df [1] 24
        

### 4. **nbr.one.rep** | all nbr (no.dup.nbr) with extra column (one rep)
#Naming each nbr accn with the one ref accn we are using
#The df no.dup.nbr add a column with only one representative ref accn (from the "picked one")    

mat = matrix(ncol = 0, nrow = 0)
more10.one.rep=data.frame(mat)

# all nbr df with an added column with only one representative 
for (a in picked.morethan1) {
  b <- more.10[grepl(a, more.10$Representative, fixed = TRUE),]
  b$Accession <- a
  more10.one.rep <- rbind(more10.one.rep, b)
}

### WHY are these duplicate rows created?? ###
dim(more.10) #[1] 20881     6
dim(more10.one.rep) #[1] 688333      6

length(unique(more10.one.rep$Neighbor)) #1] 20142 
# makes sense a little less (20881-20142) since some not in picked.morethan1 becasue non of the accns were in the fltred_ref_accns
more10.one.rep.nodup <- more10.one.rep[!duplicated(more10.one.rep$Neighbor),]

View(more10.one.rep.nodup)

## combine more10 and less10 

less.10$Accession <- less.10$Representative

nbr.one.rep <- rbind(more10.one.rep.nodup, less.10)
dim(nbr.one.rep)

View(nbr.one.rep)

###################

##### 5.**fltr.nbr** | filter for ref accessions #####
### filter for ref accessions (also needs [all_ncbi] from FILES)
      
  # all filtered ref accns in a vector
  all_ref_accn <- all_ncbi$Accession # 139 ref accns
  mat = matrix(ncol = 0, nrow = 0)
  fltr.nbr=data.frame(mat)
  
  for (a in all_ref_accn) {
    b <- nbr.one.rep[grepl(a, nbr.one.rep$Accession, fixed = TRUE),]
    fltr.nbr <- rbind(fltr.nbr, b)
  }
        View(fltr.nbr) 
        dim(fltr.nbr) #[1] 142973      7
        
        length(unique(fltr.nbr$Accession)) #[1] 97
        length(unique(nbr.one.rep$Accession)) #[1] 97
        length(unique(fltr.nbr$Neighbor)) 
        length(unique(nbr.one.rep$Neighbor))

# separate "Selected.lineage"  
##### 23 FAMILIES ####
fltr.nbr[c("family","genus","name")] <- str_split_fixed(fltr.nbr$Selected.lineage,",",3)        

##############
View(fltr.nbr) # all fltred nbr 
##############    
        

# write it out 
write.csv(fltr.nbr,fltr.nbr_path, row.names = F)
# all nbr accns in a file 
write.csv(fltr.nbr$Neighbor,fltr.nbr_accns_path, quote=F, row.names = F)

# all nbr each ref in file
nbr <- fltr.nbr
nbr$Accession <- as.factor(nbr$Accession)
length(unique(fltr.nbr$Accession)) #127
length(unique(fltr.nbr$Accession)) #93

setwd(file.path(dir,"fltr_nbr"))
d_ply(nbr, "Representative", function(x)
write.csv(unlist(x$"Neighbor"), file = paste(x$Accession[1], "csv", sep = "."), row.names = FALSE, quote = F))      
      

#############################
####### Examine accns #######

# [127] nbr_genome (fasta file with multiple entries) 
a <- read.csv("/Users/sareh/Desktop/nbr_genome127.txt")
View(a)

a127 <- a$X
b139 <- unique(all_ncbi$Accession)
c93 <- unique(fltr.nbr$Accession)

length(intersect(a127,b139)) #all 127 in the 139
setdiff(a127,b139) #character(0)

length(intersect(a127,c93)) #all 93 in 127
setdiff(c93, a127) #character(0)



