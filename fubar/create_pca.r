# fubar json output

# INPUT ####
#1. ncbi or lin
data_set = "ncbi"

#2. with_ovlp, no_olvp, pruned, all 
name = "no_ovlp"

## meta ####
meta_path <- file.path(dir_data, "meta.csv")
meta <-  read.csv(meta_path, sep = "\t")
# 3. pick which label to use for SE 
meta$surface <- meta$SE
View(meta)

# CONSISTANT ####
##### colors ####
SE <- "#FFAB91" # red "coral2" 
nonSE  <-  "#80CBC4" # blue "cornflowerblue"

col_env <-  "#79AF97FF" #"Thistle 3"
col_non_env <- "#6A6599FF" # "Navajo White 2"  

col_dna <- "lightpink3"
col_rna <- "cornflowerblue" # 

col_with_ovlp <- "Honeydew 2"  # light green 
col_without_ovlp <- "Aquamarine 4"  # dark green

##### dir ####
dir_base = "/Users/sareh/Desktop/surface"
dir_data = file.path(dir_base,"results",data_set)
dir = file.path(dir_base,"results", data_set, name)
dir_plot = file.path(dir, "plot")

# Getting rcor_wide #### 
## (1) pearson_correlation.csv  ####
# 1. /all/script$ mpirun -np 15 python3 json_pearson_mpi.py 

## (2) pearson_correlation_clean.csv  ####
#rcor_orgnl_path <- file.path(dir,"pearson_correlation.csv")
#rcor_orgnl_clean_path <- file.path(dir, "pearson_correlation_clean.csv")

#rcor_orgnl <- read.csv(rcor_orgnl_path, header = FALSE, sep = " ")
#rcor_orgnl$V1 <- gsub(".noovlp.FUBAR.json","",rcor_orgnl$V1)
#rcor_orgnl$V2 <- gsub(".noovlp.FUBAR.json","",rcor_orgnl$V2)
#rcor_orgnl$virus <- gsub("^([^_]*_[^_]*)_.*$", "\\1",rcor_orgnl$V1)
#rcor_orgnl$gene1 <- gsub("^([^_]*_[^_]*)_", "",rcor_orgnl$V1)
#rcor_orgnl$gene2 <- gsub("^([^_]*_[^_]*)_", "",rcor_orgnl$V2)
#rcor_clean <- rcor_orgnl[, c("gene1","gene2","V3")]
#colnames(rcor_clean) <- c("gene1","gene2", "pearson_correlation")
#write.csv(rcor_clean,rcor_orgnl_clean_path)

## (3) rcor_wide.csv ####
# 3. python3 json_pearson_rcorwide.py


## (4) clean it up  #### 

### [rcor_wide_original] open the file  ####
rcor_wide_path <- file.path(dir,"rcor_wide.csv")
rcor_wide_original<- read.csv(rcor_wide_path, header = TRUE, sep = ",")
  dim(rcor_wide_original) #[1] 678 679
  View(rcor_wide_original)
  
### [rcor_wide1] 1-pearson correlation = 0  ####
rcor_wide1 <- rcor_wide_original
rcor_wide1[rcor_wide1 == "NA"] <- 0.0
rcor_wide1[rcor_wide1 == ""] <- 0.0
rcor_wide1[is.na(rcor_wide1)] <- 0.0
which(is.na(rcor_wide)) # checking 
  View(rcor_wide1)
  dim(rcor_wide1)

### [rcor_wide2] Row names ####
rcor_wide2 <- rcor_wide1
rownames(rcor_wide2) <- rcor_wide2$gene1
rcor_wide2 <- rcor_wide2[,c(2:ncol(rcor_wide2))]  
  View(rcor_wide2)
  dim(rcor_wide2) #[1] 678 678
  
### [rcor_wide3] only the first ####
rcor_wide3 <- rcor_wide2
library(stringr)
rcor_wide3[] <- lapply(rcor_wide3, str_trunc, 17, ellipsis = "")
rcor_wide3[,1:nrow(rcor_wide3)] <- sapply(rcor_wide3[,1:nrow(rcor_wide3)], as.numeric)
  dim(rcor_wide3) #[1] 678 678
  View(rcor_wide3)

  
### [rcor_wide_clean] remove 4 ####
#1. YP_009505610.1 | 2. YP_009505617.1 | 3. NP_054708.1 | 4. NP_777383.1

# remove rows 
rcor_wide_clean <- rcor_wide3[!(rownames(rcor_wide3) == "NP_054708.1" | 
                            rownames(rcor_wide3) == "NP_777383.1" | 
                            rownames(rcor_wide3) == "YP_009505610.1" |
                            rownames(rcor_wide3) == "YP_009505617.1"),]

# remove columns 
rcor_wide_clean <- within(rcor_wide_clean, rm("NP_054708.1","NP_777383.1","YP_009505610.1","YP_009505617.1"))
  dim(rcor_wide_clean) #[1] 674 674
  View(rcor_wide_clean)


  
####### PCA ####
rcor_wide <- rcor_wide_clean

# I don't know why but need to run this again or get ERROR: 
  #Error in svd(x, nu = 0, nv = k) : infinite or missing values in 'x'
rcor_wide[is.na(rcor_wide)] <- 0.0
which(is.na(rcor_wide))
rcor_wide[,1:nrow(rcor_wide)] <- sapply(rcor_wide[,1:nrow(rcor_wide)], as.numeric)
  
rcor.matrix <- as.matrix(rcor_wide)
pca <- prcomp(rcor.matrix) #scale = TRUE
save(pca,file="pca.RData") #Save pca matrix

# add meta ####
pca.df <- as.data.frame(pca$x)
pca.df <- pca.df[,c("PC1","PC2","PC3","PC4")]
pca.df$gene <- rownames(pca.df)
rownames(pca.df) <- seq(1, nrow(pca.df))
meta_pca <- merge(meta, pca.df)

###### variation  ####
pca.var <- pca$sdev^2 
#amount of variation in the original data each pc accounts for 
#percentage 
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
comp1 <- pca.var.per[1]
comp2 <- pca.var.per[2]
comp3 <- pca.var.per[3]
comp4 <- pca.var.per[4]

# plot ####

library("factoextra")
library("FactoMineR")

## screen plot ####
barplot(pca.var.per, main="Scree Plot", xlab = "Principal Component", ylab="Percent Variation", xlim=c(0,20), col = "blue")

pdf(file=(file.path(dir_plot,"screen_plot.pdf")))
pca.data <- PCA(rcor_wide, scale.unit = FALSE, graph = FALSE)
fviz_eig(pca.data, addlabels = TRUE, ylim = c(0, 80),barfill="cadetblue",ncp=5)
dev.off()

## 1. surface ####
meta_path <- file.path(dir_data, "meta.csv")
meta <-  read.csv(meta_path, sep = "\t")
View(meta)

#PCA 1&2
pdf(file=(file.path(dir_plot,"pca_1_2_surfaces.pdf")))
plot(pca$x[,1], pca$x[,2], col= alpha(ifelse(meta$SEgo, SE, nonSE),0.6), main="Principal Component Analysis", xlab= paste("Component 1 (", comp1, "%)", sep = ""), ylab= paste("Component 2 (", comp2, "%)", sep = ""), cex = 2, pch = 20)
legend("topleft", col=c(SE,nonSE), legend = c("Surface Exposed","Non-Surface Exposed"), 
 pch = 20, bty='n', cex=1)
dev.off()

# PCA 3&4
pdf(file=(file.path(dir_plot,"pca_3_4_surfaces.pdf")))
plot(pca$x[,3], pca$x[,4],col= alpha(ifelse(meta$SEgo,SE,nonSE),0.6),main="Principal Component Analysis",xlab=paste("Component 3 (", comp3, "%)", sep = ""),ylab=paste("Component 4 (", comp4, "%)", sep = ""), cex = 2, pch = 20)
legend("bottomright", col=c(SE,nonSE), 
       legend = c("Surface Exposed","Non-Surface Exposed"), pch = 20, bty='n', cex=1)
dev.off()

#### ggplot ####
a <- as.data.frame(pca$x[,1:2])
a$gene <- rownames(a)
c <- meta[,c("gene","surface")]
se <- merge(a,c, by="gene")

pdf(file=(file.path(dir_plot,"pca_ggplot_se.pdf")))
ggplot(se, aes(x = PC1, y = PC2, col = surface)) +
  labs(x=paste("Component 1 (", comp1, "%)", sep = ""),
      y=paste("Component 2 (", comp2, "%)", sep = "")) +
  scale_fill_manual(labels=c("SE", "nonSE")) +
  geom_point() + theme_classic()
dev.off()



## 2. enveloped #####
#coloured by Enveloped/non-enveloped viruse
pdf(file=(file.path(dir_plot,"pca_enveloped1_2.pdf")))
plot(pca$x[,1], pca$x[,2], col= alpha(ifelse(meta$enveloped, col_env,col_non_env),0.8), main="Principal Component Analysis",xlab=paste("Component 1 (", comp1, "%)", sep = ""),ylab=paste("Component 2 (", comp2, "%)", sep = ""), cex = 2, pch = 20)
legend("topleft", col=c(col_env,col_non_env),legend = c("Enveloped","Non-Enveloped"), pch = 20, bty='n', cex=1)
dev.off()

## 3. dna/rna #####
pdf(file=(file.path(dir_plot,"pca_dna.pdf")))
plot(pca$x[,1], pca$x[,2], col= alpha(ifelse(meta$DNA, col_dna, col_rna),0.8), main="Principal Component Analysis",xlab=paste("Component 1 (", comp1, "%)", sep = ""),ylab=paste("Component 2 (", comp2, "%)", sep = ""), cex = 2, pch = 20)
legend("topleft", col=c(col_dna,col_rna),legend = c("DNA","RNA"), pch = 20, bty='n', cex=1)
dev.off()

##4. family ####
# colours 
c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "#6A3D9A", # purple
  "#FF7F00", # orange 
  "gold1", "skyblue2", "#CAB2D6", # lt purple
  "palegreen2", "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise", "green1", "darkorange4", "brown")

meta$manual_family <- as.factor(meta$manual_family)
pdf(file=(file.path(dir_plot,"pca_family.pdf")))
plot(pca$x[,1], pca$x[,2], 
     col= c25[meta$manual_family], 
     main="Principal Component Analysis", 
     xlab=paste("Component 1 (", comp1, "%)", sep = ""), 
     ylab=paste("Component 2 (", comp2, "%)", sep = ""), 
     cex = 1.5, pch = 20)
legend(8,10.5,unique(meta$manual_family), col=c25, pch=20, cex= 0.5,
       title ="Virus Family", bty='n')
dev.off()


#### LOOP ####
meta$manual_family <- as.factor(meta$manual_family)
family = factor(unique(meta$manual_family))

mycol <- rgb(211,211,211, max = 255, alpha = 130, names = "blue50")

for (index in 1:length(family)) {
  c <- rep(mycol,times=25)
  c[index] <- "dodgerblue2"
  plot_name = paste0("pca_family_",family[index], ".pdf")
  
  pdf(file=(file.path(dir_plot, plot_name)))
  plot(pca$x[,1], pca$x[,2], 
       col= c[meta$manual_family], 
       main="Principal Component Analysis", 
       xlab=paste("Component 1 (", comp1, "%)", sep = ""), 
       ylab=paste("Component 2 (", comp2, "%)", sep = ""), 
       cex = 1.5, pch = 20)
  legend(8,10.5,unique(meta$manual_family), col=c, pch=20, cex= 0.5,
         title ="Virus Family", bty='n')
  dev.off()
}

# other colour schemes 
library("viridis")  
c <- viridis(25)
pdf(file=(file.path(dir_plot,"pca_family3.pdf")))
plot(pca$x[,1], pca$x[,2], 
     col= c[meta$manual_family], 
     main="Principal Component Analysis",xlab=paste("Component 1 (", comp1, "%)", sep = ""),ylab=paste("Component 2 (", comp2, "%)", sep = ""), cex = 2, pch = 20)
legend(8,10.5,unique(meta$manual_family), col=c, pch=20, 
       cex= 0.5, title ="Virus Family", bty='n')
dev.off()

##5. treelength ####

col_tr <- (meta$tree_length + abs(min(meta$tree_length)))/max(meta$tree_length + abs(min(meta$tree_length)))

pdf(file=(file.path(dir_plot,"pca_trelength_gene12.pdf")))
plot(pca$x[,1], pca$x[,2], 
     pch = 20, cex = 1.7, col = rgb(1, 0, 0,col_tr),
     main="Principal Component Analysis",
     xlab=paste("Component 1 (", comp1, "%)", sep = ""),
     ylab=paste("Component 2 (", comp2, "%)", sep = ""))
dev.off()

pdf(file=(file.path(dir_plot,"pca_trelength_gene34.pdf")))
plot(pca$x[,3], pca$x[,4], 
     pch = 20, cex = 1.7, col = rgb(1, 0, 0,col_tr),
     main="Principal Component Analysis",
     xlab=paste("Component 1 (", comp1, "%)", sep = ""),
     ylab=paste("Component 2 (", comp2, "%)", sep = ""))
dev.off()

  
# two 
pal = colorRampPalette(c("blue", "red"))

pal <- colorRamp(c("blue", "green", "orange", "red"))
col <- rgb(pal((meta$tree_length - min(meta$tree_length)) / diff(range(meta$tree_length))), max=255)  # 2) interpolate numbers

colors <- c("red", RColorBrewer::brewer.pal(9, "Blues"))

pdf(file=(file.path(dir_plot,"pca_tree_length_gene2.pdf")))
plot(pca$x[,1], pca$x[,2], 
     col=colors, 
     main="Principal Component Analysis",
     xlab=paste("Component 1 (", comp1, "%)", sep = ""),
     ylab=paste("Component 2 (", comp2, "%)", sep = ""), 
     cex = 2, pch = 20)
#legend(8,10.5,unique(meta$manual_family), col=c, pch=20, cex= 0.5, title ="Virus Family", bty='n')
dev.off()

#### ggplot  ####

library(ggplot2)
library(hrbrthemes)

################################### FIRST METHOD 
png(file=file.path(dir_plot,"pca_trelength_gene_ggplot12.png"),width=800, height=700)
ggplot(meta_pca, aes(x=PC1, y=PC2, color=log(meta$tree_length))) + 
  geom_point(size=3, alpha = 0.8) + ggtitle("Principal Component Analysis") +  
  xlab(paste("Component 1 (", comp1, "%)", sep = "")) + 
  ylab(paste("Component 2 (", comp2, "%)", sep = "")) + 
  labs(fill = "Tree length") +
  theme_classic() +
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=20,face="bold"),
        plot.title = element_text(lineheight=2, 
                                  face="bold", size = 20, hjust = 0.5)) +
  scale_colour_gradient(low = "tomato" , high = "dodger blue")
dev.off()


################################### SECOUND METHOD 
a <- as.data.frame(pca$x[,1:2])
a$gene <- rownames(a)
b <- meta[,c("gene","tree_length")]
tr <- merge(a,b, by="gene")

pdf(file=(file.path(dir_plot,"pca_ggplot_treelength.pdf")))
ggplot(tr, aes(x = PC1, y = PC2, col = tree_length)) +
  labs(x=paste("Component 1 (", comp1, "%)", sep = ""),
       y=paste("Component 2 (", comp2, "%)", sep = ""))+
  geom_point() + theme_classic()
dev.off()


## 6. aa_entropy ####

col_aa <- (meta$aa_entropy_gene + abs(min(meta$aa_entropy_gene)))/max(meta$aa_entropy_gene + abs(min(meta$aa_entropy_gene)))

pdf(file=(file.path(dir_plot,"pca_aa_gene12.pdf")))
plot(pca$x[,1], cex = 1.7, pca$x[,2], pch = 20, 
     col = rgb(1, 0, 0, col_aa),
     main="Principal Component Analysis",
     xlab=paste("Component 1 (", comp1, "%)", sep = ""),
     ylab=paste("Component 2 (", comp2, "%)", sep = ""))
dev.off()

pdf(file=(file.path(dir_plot,"pca_aa_gene34.pdf")))
plot(pca$x[,3], cex = 1.7, pca$x[,4], pch = 20, 
     col = rgb(1, 0, 0, col_aa),
     main="Principal Component Analysis",
     xlab=paste("Component 1 (", comp1, "%)", sep = ""),
     ylab=paste("Component 2 (", comp2, "%)", sep = ""))
dev.off()


#### ggplot #####

################################ FIRST METHOD 
png(file=file.path(dir_plot,"pca_aa_gene_ggplot12.png"),width=800, height=700)
ggplot(meta_pca, aes(x=PC1, y=PC2, color=log(meta$aa_entropy_gene))) + 
  geom_point(size=3, alpha = 0.8) + ggtitle("Principal Component Analysis") +  
  xlab(paste("Component 1 (", comp1, "%)", sep = "")) + 
  ylab(paste("Component 2 (", comp2, "%)", sep = "")) + 
  labs(fill = "Tree length") +
  theme_classic() +
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=20,face="bold"),
        plot.title = element_text(lineheight=2, 
                                  face="bold", size = 20, hjust = 0.5)) +
  scale_colour_gradient(low = "tomato" , high = "dodger blue")
dev.off()


################################ SECOUND METHOD 
a <- as.data.frame(pca$x[,1:2])
a$gene <- rownames(a)
b <- meta[,c("gene","aa_entropy_gene")]
aa <- merge(a,b, by="gene")

pdf(file=(file.path(dir_plot,"pca_ggplot_aa_entropy.pdf")))
ggplot(aa, aes(x = PC1, y = PC2, col = aa_entropy_gene)) +
  labs(x=paste("Component 1 (", comp1, "%)", sep = ""),
       y=paste("Component 2 (", comp2, "%)", sep = ""))+
  geom_point() + theme_classic()
dev.off()


## Select colours ####
meta <- meta_pca

selectSE <- c("NP_040154.2", "YP_009513268.1")
#selectnonSE <- c("YP_007947989.1", "YP_009505441.1")

meta$select <- ifelse(meta$gene %in% selectSE, "SE", NA)

meta[291,73] <- "non-SE" # meta[291,1] "YP_007947989.1"
meta[327,73] <- "non-SE" # meta[327,1] "YP_009505441.1"

mycol <- rgb(211,211,211, max = 255, alpha = 50, names = "blue50")
c <- rep(mycol,times=nrow(meta))

c[meta$select == "SE"] <- "red"
c[meta$select == "non-SE"] <- "blue"

pdf(file=(file.path(dir_plot,"pca_select_colour.pdf")))
plot(pca$x[,1], pca$x[,2], 
     col= c, 
     main="Principal Component Analysis", 
     xlab=paste("Component 1 (", comp1, "%)", sep = ""), 
     ylab=paste("Component 2 (", comp2, "%)", sep = ""), 
     cex = 1.5, pch = 20)

legend(8,10.5,unique(meta$select), col=c, pch=20, cex= 0.5,
       title ="Virus Family", bty='n')
dev.off()




## ADD THE COLOURS ON TOP 

mycol <- rgb(211,211,211, max = 255, alpha = 130, names = "blue50")
c <- rep(mycol,times=nrow(meta))

pdf(file=(file.path(dir_plot,"pca_select_colour_added.pdf")))
plot(pca$x[,1], pca$x[,2], 
     col= c, 
     main="Principal Component Analysis", 
     xlab=paste("Component 1 (", comp1, "%)", sep = ""), 
     ylab=paste("Component 2 (", comp2, "%)", sep = ""), 
     cex = 1.5, pch = 20)
legend("topleft", col=c(SE,nonSE), legend = c("Surface Exposed","Non-Surface Exposed"), 
       pch = 20, bty='n', cex=1)

points(-7.866103, 0.6109093 , col = SE, pch = 20, cex = 3)
points( 9.804051,1.297884   , col = SE, pch = 20, cex = 3)
points(-7.868593, 0.519468  , col = nonSE, pch = 20, cex = 3)
points(9.831324, 0.8801241  , col = nonSE, pch = 20, cex = 3)

dev.off()



SE: "NP_040154.2" (-7.866103, 0.6109093)
non-SE: "YP_007947989.1" (-7.868593, 0.519468)

SE: "YP_009513268.1" (9.804051,1.297884)
non-SE: "YP_009505441.1" (9.831324, 0.8801241)


















