# **Virus Accessions**
# Getting the virus accessions 
# Automate this from 
# https://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=10239&host=human
# 2022: Viruses - 11553 complete genomes | filtered by host 'human': 527 genomes
#	Download
#		1* table (taxid.tbl) 
#			[/Users/sareh/surfaces/r_surfaces/INDEX/virus_info/raw_ncbi_Human_virus.tbl]
#      ## Complete genomes: Viruses (taxid 10239) filtered by host 'human'
#			"Genome"	"Accession"	"RefSeq type"	"Source information"	"Number of segments"	"Genome length"	"Number of proteins"	"Genome  "Neighbors"	"Host"	"Date completed"	"Date updated"
#		2* list of accessions (taxid.acc_lst) 
#			  ref-seq accn
#		3* Genome neighbours data (taxid.nbr)
#			## Neighbors data for complete genomes: Viruses (taxid 10239)
#			## Columns:	"Representative"	"Neighbor"	"Host"	"Selected lineage"	"Taxonomy name"	"Segment name"

#raw donwloaded files (2020)
taxid.tbl <- "./virus_info/taxid10239.tbl" # 935 total lines | 408 segmented (start with space) | 515 lines not starting with space [ncbi.csv]
taxid.acc_lst <- "./virus_info/taxid10239.acc_lst"
taxid.nbr <- "./virus_info/taxid10239.nbr"

# ncbi_Human_virus.tbl without segmented viruses
# Get rid of the lines starting with tab (expandable rows on the website)
ncbi.tbl <- read.table("/Users/sareh/surfaces/r_surfaces/INDEX/virus_info/ncbi.csv", fill = TRUE, sep=",")
colnames(ncbi.tbl) <- c("Accession","RefSeq.type","source.information","segm","length","protein","Neighbors","host","created","updated","empty")
ncbi.tbl$Neighbors <- as.numeric(ncbi.tbl$Neighbors)
ncbi.tbl$protein <- as.numeric(ncbi.tbl$protein)
dim(ncbi.tbl) #[1] 457  11

## **SUBSETTING**
#remove segmented viruses 
ncbi.tbl1 <- ncbi.tbl[ncbi.tbl$Accession!= "-", ]
dim(ncbi.tbl1) #[1] 368  11
View(ncbi.tbl1)

#remove incomplete seq viruses 
ncbi.tbl2 <- ncbi.tbl1[ncbi.tbl1$RefSeq.type == "complete",]
dim(ncbi.tbl2) #[1] 365  11
View(ncbi.tbl2)

#viruses > 100 neighbours 
ncbi.tbl3 <- ncbi.tbl2[ncbi.tbl2$Neighbors >= 100, ]
dim(ncbi.tbl3) #[1] 167  11
View(ncbi.tbl3) #this included Neighbors == NA 
ncbi.tbl3 <- ncbi.tbl3[!is.na(ncbi.tbl3$Accession),] #[1] 83  11

#viruses > 3 proteins 
ncbi.tbl4 <- ncbi.tbl3[ncbi.tbl3$protein >= 3,  ]
dim(ncbi.tbl4) #[1] 48 11
View(ncbi.tbl4)

#saving the file (non-empty)
ncbi.viruses <- ncbi.tbl4[!is.na(ncbi.tbl4$Accession), ]
dim(ncbi.tbl) #[1] 48 11
View(ncbi.tbl)


## **combine seg + non-seg viruses info**
seg.virus.name$host <- "human"

colnames(seg.virus.name)<-c("x","name","Accession","length.nt","protein","neighbors","source","updated","created","segment" )
colnames(ncbi.viruses.32)<-c("name","Accession","RefSeq.type","source","segment","length.nt","protein","neighbors","host","created","updated","empty")

seg.virus <- seg.virus.name[,c("name","Accession","length.nt","protein","neighbors","updated","created","segment")]
nonseg.virus <- ncbi.viruses.32[,c("name","Accession","length.nt","protein","neighbors","updated","created","segment")]

nonseg.virus$length.nt <- gsub("(nt)","",nonseg.virus$length.nt)

#View(seg.virus) # segmented: seg.virus.name
#View(nonseg.virus) # non-segmented: ncbi.viruses.32

virus_info <- rbind(seg.virus,nonseg.virus)

# filtered nbr table 
nbr.2020 <- read.csv("/Users/sareh/surfaces/r_surfaces/INDEX/virus_info/nbr.2020")

all.fltr.nbr <- read.csv("/Users/sareh/surfaces/r_surfaces/INDEX/virus_info/all.fltr.nbr.csv")







