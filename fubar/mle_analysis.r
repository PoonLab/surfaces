# INPUT #### 
library(ggplot2)
library(plotly)
library(ggpubr)

#1. ncbi or lin
data_set = "ncbi"

#2. with_ovlp, no_olvp, pruned, all 
name = "no_ovlp"

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

dir_plot = file.path(dir_data, name, "plot")

# info ####

##### Input (meta) ####
meta_path <- file.path(dir_data, "meta.csv")
meta <-  read.csv(meta_path, sep = "\t")

# 3. pick which label to use for SE 
meta$surface <- meta$SE

View(meta)

# PLOTS ##########

## 1. virus info table ####
# name.x & name.y are different 
meta.virus <- meta[!duplicated(meta$name_ncbi),] # 54

all_virus_info <- meta.virus[,c('virus',"name_manual","DNA",'enveloped',
                                "genome_length.nt","manual_family",
                                "manual_order", "neighbors")]

virus_table <- all_virus_info[,c("virus","genome_length.nt",
                                 'name_manual',"manual_family","manual_order","neighbors")]

colnames(virus_table) <- c("Accession","genome length",
                           "virus name","family","order","neighbour genome")
View(virus_table) # for thesis 

# SAVE virus_table
virus_tbl_path <- file.path(dir_plot,"virus_table.csv")
write.table(virus_table, virus_tbl_path, quote = FALSE, sep = ",", row.names = FALSE)

### add nbr (fltr.nbr) ####
fltr.nbr_path <- file.path(dir_base, "raw_ncbi_files", "fltr_nbr.csv")
fltr.nbr_accns_path <- file.path(dir_ncbi, "fltr_nbr_accns.csv")
fltr.nbr <- read.csv(fltr.nbr_path) # colnames don't show up properly
names(fltr.nbr)[names(fltr.nbr) == 'family'] <- 'family_nbr'
names(fltr.nbr)[names(fltr.nbr) == 'genus'] <- 'genus_nbr'
names(fltr.nbr)[names(fltr.nbr) == 'name'] <- 'name_nbr'
View(fltr.nbr) #[1] 142973     10

# NOT COMPLETE 
# ISSUE is fltr nbr only has ref_accn 
# ref_accn is for ever segment


## 1. pie chart ########## 
pie_path <- file.path(dir_plot, "pie_chart.pdf")

###################### organise dataframe 
# frequency of GENES (family)
meta$manual_family <- as.factor(meta$manual_family)
gene_df <- as.data.frame(table(meta$manual_family))
colnames(gene_df) <- c("family","gene_freq")

# frequency of VIRUSES (family)
virus_table$family <- as.factor(virus_table$family)
virus_df <- as.data.frame(table(virus_table$family))
colnames(virus_df) <- c("family","virus_freq")

# frequency of virus ORDER
virus_table$virus_order <- as.factor(virus_table$order)
order_df <- as.data.frame(table(virus_table$virus_order))
colnames(order_df) <- c("order","order_freq")

pie_freq <- merge(gene_df,virus_df, by = "family")

# add more info 
meta.family <- meta[!duplicated(meta$manual_family),]

all_meta.family_freq <- merge(meta.family, pie_freq,
           by.x = "manual_family", by.y = "family")

meta.family_freq <- all_meta.family_freq[,c("manual_family","gene_freq","virus_freq","DNA","enveloped")]

###################### info

virus_family.df <- as.data.frame(virus_family)

# no legend virus 
png(file=file.path(dir_plot,"piechart_virus_family.png"),
    width=800, height=800)
pie(virus_family, cex=1, border="grey", main="Virus Families", col=color, cex.main=1.7)
dev.off()

# with legend virus
png(file=file.path(dir_plot,"piechart_virus_family_legend.png"),
    width=1000, height=800)
pie(virus_family, cex=0.1, border="grey", main="Virus Families", col=color, cex.main=1.7, labels = NA)
legend("topright", inset = .05, title = "Virus Families",legend= virus_family.df$Var1  ,fill = color,cex =1)
dev.off()


### nested pie chart ####
# load library plotly
library(plotly)

nested_pie <- plot_ly(pie_freq) %>%
  add_pie(labels = ~`family`, values = ~`gene_freq`, 
          type = 'pie', hole = 0.7, sort = FALSE, 
          textinfo = "none",
          insidetextorientation = "radial",
          marker = list(line = list(width = 3))) %>%
  add_pie(pie_freq, 
          labels = ~`family`, values = ~`virus_freq`, 
          domain = list(
            x = c(0.15, 0.85),
            y = c(0.15, 0.85)),
          textinfo = "label",
          insidetextorientation = "radial",
          sort = FALSE) %>% 
  layout (font=list(size = 10),
          legend = list(x = 0.75, y = 0.5, 
                        title=list(text='<b> Virus Families </b>'))) 

nested_pie

np_path <- file.path(dir_plot, "nested_pie.svg")
orca(nested_pie, file = np_path)

#### family & order ####

nested_pie <- plot_ly(pie_freq) %>%
  add_pie(labels = ~`family`, values = ~`gene_freq`, 
          type = 'pie', textinfo = "label",
          insidetextorientation = "radial",
          hole = 0.7, sort = T,
          marker = list(line = list(width = 3))) %>%
  add_pie(order_df, labels = ~`order`, values = ~`order_freq`, 
          domain = list(
            x = c(0.15, 0.85),
            y = c(0.15, 0.85)),
          textinfo = "label",
          insidetextorientation = "radial",
          sort = T) %>% 
  layout (font=list(size = 10))

nested_pie

#### dna and family ####


#### env and family ####



# Distribution  ########## 
## 2. neigbours  ########## 
pdf(file=file.path(dir_plot,"dist_neighbours.pdf"))
hist(log(meta$neighbors), breaks = 30, col = "cadetblue",
     main = paste("Distribution of the number of neighbour genomes"),
     xlab = "log transformed number of neighbour genomes")
dev.off()

## 3. genes #####
df <- as.data.frame(table(meta$name_manual))
pdf(file=file.path(dir_plot,"gene_freq.pdf"))
barplot(table(meta$name_manual))
dev.off()

## 2. codon length ########## 
h_without <- hist(log(meta$codon_no_ovlp,10), 
                  main = "Distribution of Codon Length without ovlp", 
                  breaks = 25, xlab = "log transformed codon length")

h_with <- hist(log(meta$codon_with_ovlp,10), 
               main = paste("Distribution of Original Codon Length"), 
               breaks = 25, xlab = "log transformed codon length")

# colours
mycol1 <- rgb(0,0,1,1/4)
  # rgb(224, 238, 224, max = 255, alpha = 50) # with
mycol2 <- rgb(1,0,0,1/4)
  # rgb(69, 139, 116, max = 255, alpha = 50) # without

# plot both hist 
pdf(file=file.path(dir_plot,"with and without ovlp codon distribution3.pdf"))
plot(h_without, col =mycol1, xlab='log transformed codon length',main = "")
plot(h_with, col =  mycol2, add = TRUE, main = "")
#add legend
legend('topright', border = FALSE,bty = "n",
       c('with ovlp','without ovlp'), 
       fill=c(mycol1, mycol2))
dev.off()

### t-test ####
summary(meta$codon_no_ovlp)
summary(meta$codon_with_ovlp)

#Welch Two Sample t-test
t.test(meta$codon_no_ovlp, meta$codon_with_ovlp)

#Two Sample t-test
t.test(meta$codon_no_ovlp, meta$codon_with_ovlp, var.equal = TRUE)

#Paired t-test
t.test(meta$codon_no_ovlp, meta$codon_with_ovlp, paired = TRUE)




# Bar plots ####
## 3. Env/SE bar plot##########  
pdf(file=file.path(dir_plot,"env_se_stackedbarplot.pdf"))
b <- barplot(table(meta$SE,meta$enveloped), col=c(nonSE,SE), border="white", 
             space=0.04, font.axis=2, cex.names=1.5,
             xlab = "virus", names.arg=c("non-enveloped","enveloped"),
             cex.lab=1.5)

legend("topleft", fill=c(SE,nonSE), legend=c("SE", "non-SE"),
       bty = "n", cex=1.5, border = FALSE)

text(x = b , y = c(5,100,110,450), pos = 3,cex = 1.3,
     label = table(meta$enveloped,meta$SE),col="white")

dev.off()

## 4. DNA/SE bar plot########## 
pdf(file=file.path(dir_plot,"dna_se_stackedbarplot.pdf"))
x <- barplot(table(meta$SE, meta$DNA), col=c(nonSE,SE), border="white", 
             space=0.04, font.axis=2, ylim = c(0,500), 
             xlab = "virus", names.arg=c("DNA","RNA"),
             cex.names=1.5, cex.lab=1.5)

legend("topleft", fill=c(SE,nonSE), legend=c("SE", "non-SE"),
       bty = "n", cex=1.5,border = FALSE)

text(x = b , y = c(50,100,175,400), pos = 3,cex = 1.3,
     label = table(meta$DNA,meta$SE),col="white")

dev.off()


## 5. DNA/ENV bar plot########## 
dna_se <- meta[,c("DNA","enveloped")]
pdf(file=file.path(dir_plot,"dna_enveloped_stackedbarplot.pdf"))
x <- barplot(table(dna_se), col=c(col_dna,col_rna), border="white", 
             space=0.04, font.axis=2, ylim = c(0,500), 
             xlab = "virus", names.arg=c("non-enveloped","enveloped"),
             cex.names=1.5, cex.lab=1.5)

legend("topleft", fill=c(col_dna,col_rna), legend=c("RNA","DNA"),
       bty = "n", cex=1.5,border = FALSE)
text(x = b , y = c(1,80,80,350), pos = 3,cex = 1.3,
     label = table(meta$enveloped,meta$DNA), col="white")
dev.off()

## 6. mosaic plot ####

View(data)












# Box plots ####
## 1. DNA dN/dS boxplot##########  

compare_means(mean.dN.mean.dS~DNA, meta, 
              method = "wilcox.test", paired = FALSE,
              group.by = NULL, ref.group = NULL)

p <- ggboxplot(meta, 
               x = "DNA", y = "mean.dN/mean.dS",
               fill = c(col_dna,col_rna))

stat_box_data <- function(y, upper_limit = max(meta$mean.dN.mean.dS) * 0.6) {
  return( 
    data.frame(
      y = -0.000001 * upper_limit,
      label = paste('n =', length(y), '\n'), '\n'))
}

pdf(file=file.path(dir_plot,"dna_rna_env_wilcoxon5.pdf"))
p + stat_compare_means(label.y = 0.95, size = 8) + 
  geom_jitter(shape=20, position=position_jitter(0.3), alpha=.1) + 
  ylab("dN / dS") + xlab("genes") +
  stat_summary(fun.data = stat_box_data, geom = "text", 
               hjust = 0.5, vjust = 1.1) +   
  theme(text = element_text(size=16),
    axis.title.y = element_text(vjust = +1.5, size=16),
    axis.title.x = element_text(vjust = +0.001, size=16)) +
  scale_x_discrete(labels= c("RNA","DNA"),)
dev.off()


## 2. Env dN/dS boxplot ##########  
meta_nona <- meta[!is.na(meta$enveloped),]

compare_means(mean.dN.mean.dS~enveloped, meta_nona, method = "wilcox.test", paired = FALSE,group.by = NULL, ref.group = NULL)

p <- ggboxplot(meta_nona, 
               x = "enveloped", y = "mean.dN/mean.dS",
               fill = c(col_non_env,col_env))

stat_box_data <- function(y, upper_limit = max(meta$mean.dN.mean.dS) * 0.6) {
  return( 
    data.frame(
      y = -0.000001 * upper_limit,
      label = paste('n =', length(y), '\n'), '\n'))
}

pdf(file=file.path(dir_plot,"env&nonenv_wilcoxon5.pdf"))
p + stat_compare_means(label.y = 0.95, size = 8) + 
  geom_jitter(shape=20, position=position_jitter(0.3), alpha=.1) + 
  ylab("dN / dS") + xlab("genes") +
  stat_summary(fun.data = stat_box_data, geom = "text", 
               hjust = 0.5, vjust = 1.1) +   
  theme(text = element_text(size=16),
        axis.title.y = element_text(vjust = +1.5, size=16),
        axis.title.x = element_text(vjust = +0.001, size=16)) +
  scale_x_discrete(labels= c("Non-Enveloped","Enveloped"))
dev.off()

## 3. dN/dS########## 
meta_nona <- meta[!is.na(meta$surface),]

compare_means(mean.dN.mean.dS~surface, meta_nona, method = "wilcox.test", paired = FALSE,group.by = NULL, ref.group = NULL)

p <- ggboxplot(meta_nona, 
               x = "surface", y = "mean.dN.mean.dS", 
               fill= c(nonSE,SE))

stat_box_data <- function(y, upper_limit = max(meta$mean.dN.mean.dS) * 1.15) {
  return( 
    data.frame(
      y = 0.0001 * upper_limit,
      label = paste('n =', length(y), '\n'), '\n'))
}

pdf(file=file.path(dir_plot,"dNdS_wilcoxon5.pdf"))
p + stat_compare_means(label.y = 0.95, size = 8) + 
  geom_jitter(shape=20, position=position_jitter(0.3), alpha=.1) + 
  ylab("dN / dS ") + xlab("genes") + 
  stat_summary(fun.data = stat_box_data, 
               geom = "text", hjust = 0.5, vjust = 0.95) + 
  theme(text = element_text(size=16),
        axis.title.y = element_text(vjust = +1.5, size=16),
        axis.title.x = element_text(vjust = +0.001, size=16)) +
    scale_x_discrete(labels= c("non-SE","SE"))
dev.off()

## 3. dN-dS########## 
meta$dn.minus.ds <- meta$mean.dN - meta$mean.dS

meta_nona <- meta[!is.na(meta$surface),]

compare_means(dn.minus.ds~surface, meta_nona, method = "wilcox.test", paired = FALSE,group.by = NULL, ref.group = NULL)

p <- ggboxplot(meta_nona, 
               x = "surface", y = "dn.minus.ds", 
               fill= c(nonSE,SE))

stat_box_data <- function(y, upper_limit = max(meta$dn.minus.ds) * 1.15) {
  return( 
    data.frame(
      y = 0.0001 * upper_limit,
      label = paste('n =', length(y), '\n'), '\n'))
}

pdf(file=file.path(dir_plot,"dN-dS_wilcoxon5.pdf"))
p + stat_compare_means(label.y = -4, label.x = 1.21,size = 8) + 
  geom_jitter(shape=20, position=position_jitter(0.3), alpha=.1) + 
  ylab("dN - dS ") + xlab("genes") + 
  stat_summary(fun.data = stat_box_data, 
               geom = "text", hjust = 0.95, vjust = 0.7) + 
  theme(text = element_text(size=16),
        axis.title.y = element_text(vjust = +1.5, size=16),
        axis.title.x = element_text(vjust = +0.001, size=16)) +
  scale_x_discrete(labels= c("non-SE","SE"))
dev.off()


## 4. dN ##########  
meta_nona <- meta[!is.na(meta$surface),]

compare_means(mean.dN~surface, meta_nona, method = "wilcox.test", paired = FALSE,group.by = NULL, ref.group = NULL)

p <- ggboxplot(meta_nona, 
               x = "surface", y = "mean.dN", 
               fill= c(nonSE,SE))


stat_box_data <- function(y, upper_limit = max(meta$mean.dN) * 1.15) {
  return( 
    data.frame(
      y = 0.0001 * upper_limit,
      label = paste('n =', length(y), '\n'), '\n'))
}

pdf(file=file.path(dir_plot,"dN_wilcoxon.pdf"))
p + stat_compare_means(label.y = 1.3, label.x = 1.197,size = 8) + 
  geom_jitter(shape=20, position=position_jitter(0.3), alpha=.1) + 
  ylab("mean dN") + xlab("genes") + 
  stat_summary(fun.data = stat_box_data, geom = "text", 
               hjust = 0.5, vjust = 0.8)  + 
  coord_cartesian(ylim = c(0, 2.5)) + 
  theme(text = element_text(size=16),
        axis.title.y = element_text(vjust = +1.5, size=16),
        axis.title.x = element_text(vjust = +0.001, size=16)) +
  scale_x_discrete(labels= c("non-SE","SE"))
dev.off()


## 5. dS ##########
meta_nona <- meta[!is.na(meta$surface),]

compare_means(mean.dS~surface, meta_nona, method = "wilcox.test", paired = FALSE,group.by = NULL, ref.group = NULL)

p <- ggboxplot(meta_nona, 
               x = "surface", y = "mean.dS", 
               fill= c(nonSE,SE))

stat_box_data <- function(y, upper_limit = max(meta$mean.dS) * 1.15 ) {
  return( 
    data.frame(
      y = 0.0001 * upper_limit,
      label = paste('n =', length(y), '\n'), '\n'))
}

pdf(file=file.path(dir_plot,"dS_wilcoxon.pdf"))
p + stat_compare_means(label.y = 6.5, label.x = 1.25,size = 8) + 
  geom_jitter(shape=20, position=position_jitter(0.3), alpha=.1) + 
  ylab("mean dS") + xlab("genes") + 
  stat_summary(fun.data = stat_box_data, geom = "text", 
               hjust = 0.5, vjust = 0.8) + 
  theme(text = element_text(size=16),
        axis.title.y = element_text(vjust = +1.5, size=16),
        axis.title.x = element_text(vjust = +0.001, size=16)) +
  scale_x_discrete(labels= c("non-SE","SE"))
dev.off()

# percentage of sites #####

## 1. % pos ##########  
meta_nona <- meta[!is.na(meta$surface),]

### dN/dS ####
compare_means(pos.sites.codon~surface, meta_nona, method = "wilcox.test", paired = FALSE,group.by = NULL, ref.group = NULL)

p <- ggboxplot(meta_nona, 
               x = "surface", 
               y = "pos.sites.codon", 
               fill= c(nonSE,SE))

stat_box_data <- function(y, upper_limit = max(meta$pos.sites.codon) * 1.15 ) {
  return( 
    data.frame(
      y = 0.0001 * upper_limit,
      label = paste('n =', length(y), '\n'), '\n'))
}

pdf(file=file.path(dir_plot,"per_pos_sites.pdf"))
p + stat_compare_means(label.y = 0.1, label.x = 1.25,size = 8) + 
  geom_jitter(shape=20, position=position_jitter(0.3), alpha=.1) + 
  ylab("Percentage of Positive Sites") + xlab("genes") + 
  stat_summary(fun.data = stat_box_data, geom = "text", 
               hjust = 0.5, vjust = 1.2) + 
  theme(text = element_text(size=16),
        axis.title.y = element_text(vjust = +1.5, size=16),
        axis.title.x = element_text(vjust = +0.001, size=16)) +
  scale_x_discrete(labels= c("non-SE","SE"))
dev.off()


### dna/rna ####

compare_means(pos.sites.codon~DNA, meta_nona, method = "wilcox.test", paired = FALSE,group.by = NULL, ref.group = NULL)

p <- ggboxplot(meta_nona, 
               x = "DNA", 
               y = "pos.sites.codon", 
               fill= c(col_dna,col_rna))

stat_box_data <- function(y, upper_limit = max(meta$pos.sites.codon) * 1.15 ) {
  return( 
    data.frame(
      y = 0.0001 * upper_limit,
      label = paste('n =', length(y), '\n'), '\n'))
}

pdf(file=file.path(dir_plot,"per_pos_sites_dna.pdf"))
p + stat_compare_means(label.y = 0.1, label.x = 1.25,size = 8) + 
  geom_jitter(shape=20, position=position_jitter(0.3), alpha=.1) + 
  ylab("Percentage of Positive Sites") + xlab("genes") + 
  stat_summary(fun.data = stat_box_data, geom = "text", 
               hjust = 0.5, vjust = 1.2) + 
  theme(text = element_text(size=16),
        axis.title.y = element_text(vjust = +1.5, size=16),
        axis.title.x = element_text(vjust = +0.001, size=16)) +
  scale_x_discrete(labels= c("DNA","RNA"))
dev.off()

### env/non-env ####
compare_means(pos.sites.codon~enveloped, meta_nona, method = "wilcox.test", paired = FALSE,group.by = NULL, ref.group = NULL)

p <- ggboxplot(meta_nona, 
               x = "enveloped", 
               y = "pos.sites.codon", 
               fill= c(col_env,col_non_env))

stat_box_data <- function(y, upper_limit = max(meta$pos.sites.codon) * 1.15 ) {
  return( 
    data.frame(
      y = 0.0001 * upper_limit,
      label = paste('n =', length(y), '\n'), '\n'))
}

pdf(file=file.path(dir_plot,"per_pos_sites_enveloped.pdf"))
p + stat_compare_means(label.y = 0.1, label.x = 1.25,size = 8) + 
  geom_jitter(shape=20, position=position_jitter(0.3), alpha=.1) + 
  ylab("Percentage of Positive Sites") + xlab("genes") + 
  stat_summary(fun.data = stat_box_data, geom = "text", 
               hjust = 0.5, vjust = 1.2) + 
  theme(text = element_text(size=16),
        axis.title.y = element_text(vjust = +1.5, size=16),
        axis.title.x = element_text(vjust = +0.001, size=16)) +
  scale_x_discrete(labels= c("enveloped","non-enveloped"))
dev.off()




## 2. % neg ########## 

meta_nona <- meta[!is.na(meta$surface),]

compare_means(neg.sites.codon~surface, meta_nona, method = "wilcox.test", paired = FALSE,group.by = NULL, ref.group = NULL)

p <- ggboxplot(meta_nona, 
               x = "surface", 
               y = "neg.sites.codon", 
               fill= c(nonSE,SE))

stat_box_data <- function(y, upper_limit = max(meta$neg.sites.codon) * 1.15 ) {
  return( 
    data.frame(
      y = 0.0001 * upper_limit,
      label = paste('n =', length(y), '\n'), '\n'))
}

pdf(file=file.path(dir_plot,"per_neg_sites.pdf"))
p + stat_compare_means(label.y = 0.8, label.x = 1.23,size = 7.5) + 
  geom_jitter(shape=20, position=position_jitter(0.3), alpha=.1) + 
  ylab("Percentage of Negative Sites") + 
  xlab("genes") + 
  stat_summary(fun.data = stat_box_data, geom = "text", 
               hjust = 0.5, vjust = 1.2) + 
  theme(text = element_text(size=16),
        axis.title.y = element_text(vjust = +1.5 ,size=16),
        axis.title.x = element_text(vjust = +0.001, size=16)) +
  scale_x_discrete(labels= c("non-SE","SE"))
dev.off()


## 3. % pos & neg ##########
meta_nona <- meta[!is.na(meta$surface),]

compare_means(per_pos_neg_sites~surface, meta_nona, method = "wilcox.test", paired = FALSE,group.by = NULL, ref.group = NULL)

p <- ggboxplot(meta_nona, 
               x = "surface", 
               y = "per_pos_neg_sites", 
               fill= c(nonSE,SE))

stat_box_data <- function(y, upper_limit = max(meta$per_pos_neg_sites) * 1.15 ) {
  return( 
    data.frame(
      y = 0.0001 * upper_limit,
      label = paste('n =', length(y), '\n'), '\n'))
}

pdf(file=file.path(dir_plot,"per_pos_neg_sites.pdf"))
p + stat_compare_means(label.y = 0.8, label.x = 1.2,size = 8) + 
  geom_jitter(shape=20, position=position_jitter(0.3), alpha=.1) + 
  ylab("Percentage of Positive and Negative Sites") + 
  xlab("genes") + 
  stat_summary(fun.data = stat_box_data, geom = "text", 
               hjust = 1, vjust = 1.2) + 
  theme(text = element_text(size=16),
        axis.title.y = element_text(vjust = +1.5 ,size=16),
        axis.title.x = element_text(vjust = +0.001, size=16))+
  scale_x_discrete(labels= c("non-SE","SE"))
dev.off()





## 10. regression ##########  

# Test
result.test <- glm(mean.dN.mean.dS ~ DNA + enveloped + surface, 
                  family="Gamma", data=meta)
summary(result.test)


# DNA 
 # AIC: -730.19
result.dna <- glm(mean.dN.mean.dS ~ DNA, 
                  family="Gamma", data=meta)
summary(result.dna)

# ENV 
  # AIC: -527.63
result.env <- glm(mean.dN.mean.dS ~ enveloped, 
                family="Gamma", data=meta)
summary(result.env)


# DNA & ENV 
  # AIC: -775.67
result.a <- glm(mean.dN.mean.dS ~ DNA + enveloped, 
                family="Gamma", data=meta)
summary(result.a)


# ALL 
  # AIC: -777.88
result.b <- glm(mean.dN.mean.dS ~ DNA + enveloped + surface, 
                family="Gamma", data=meta)
summary(result.b)


# mixed effect model 

library(BSDA)
EDA(residuals(result.b))

library(lme4)

fit <- glmer(mean.dN.mean.dS ~ DNA + enveloped + surface + (1 | name_manual), data=meta, family='Gamma')

summary(fit)



# ALL ####

[] DARKER COLOURS
[] TRANSPARENT COLOURS 

## all ####
pdf(file=file.path(dir_plot, "all.pdf")) 

file=file.path(dir_plot, "all.png")
png(file, width = 1000, height = 600)
ggplot(meta, aes(x=name.x, y=mean.dN.mean.dS, colour=surface)) + 
  geom_point(size=1.2) + #  add alpha = 0.6 for transparency
  xlab(" ") + ylab("dN/dS") +
  scale_color_manual(labels = c("non-SE","SE"), values = c(nonSE,SE)) + 
  facet_grid(~family.y, switch = "x", scales = "free_x", space = "free_x") + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6), 
        legend.position = c(.9, .9), 
        strip.text.x = element_text(size = 3))
dev.off()


### horizontal ####

pdf(file=file.path(dir_plot, "all_horizontal.pdf")) 

file=file.path(dir_plot, "all_h.png")
png(file, width = 900, height = 1000)
ggplot(meta, aes(y=name.x, x=mean.dN.mean.dS, colour=surface)) + 
  geom_point(size=1.8) + #  add alpha = 0.6 for transparency
  xlab("dN/dS") + ylab("") +
  scale_color_manual(labels = c("non-SE","SE"), values = c(nonSE,SE)) + 
  facet_grid(family.y ~., scales = "free", space = "free") +
  theme(strip.text.y = element_text(angle = 0,size = 7), 
        axis.text = element_text(size = 7)) +
  geom_vline(xintercept = 1,col = "navy", lwd = 1)

dev.off()


## env ####
env <- meta[meta$enveloped == TRUE,]

env$surface <- env$SE
pdf(file=file.path(dir_plot, "all_env_genes.pdf")) 
par(mar=c(10, 15, 3, 15)) 

file=file.path(dir_plot, "all_env.png")
png(file, width = 1000, height = 600)
ggplot(env, aes(x=name.x, y=mean.dN.mean.dS, colour=surface)) + 
  geom_point(size=1.2) + #  add alpha = 0.6 for transparency
  xlab("") + ylab("dN/dS") +
  scale_color_manual(labels = c("non-SE","SE"), values = c(nonSE,SE)) + 
  facet_grid(~family.y, switch = "x", scales = "free_x", space = "free_x") + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6), 
        legend.position = c(.9, .9), 
        strip.text.x = element_text(size = 3))
dev.off()

### horizontal ####
pdf(file=file.path(dir_plot, "env_genes_h.pdf")) 
par(mar=c(10, 15, 3, 15)) 

file=file.path(dir_plot, "all_env_h.png")
png(file, width = 1000, height = 600)
ggplot(env, aes(y=name_manual, x=`mean.dN/mean.dS`, colour=surface)) + 
  geom_point(size=1.8, alpha = 0.8) + 
  xlab("dN/dS") + ylab("") +
  scale_color_manual(labels = c("non-SE","SE"), values = c(nonSE,SE)) + 
  facet_grid(manual_family ~., scales = "free", space = "free") +
  theme(strip.text.y = element_text(angle = 0),
        axis.text.y = element_text(size = 10,angle = 10)) +
  geom_vline(xintercept = 1,col = "navy", lwd = 1)
dev.off()


## non-env ####
non_env <- meta[meta$enveloped == FALSE,]
non_env$surface <- non_env$SE

pdf(file=file.path(dir_plot, "all_non_env_genes.pdf")) 
par(mar=c(10, 15, 3, 15)) 

file=file.path(dir_plot, "all_non_env.png")
png(file, width = 1000, height = 600)
ggplot(non_env, aes(x=name.x, y=mean.dN.mean.dS, colour=surface)) + 
  geom_point(size=1.2) + #  add alpha = 0.6 for transparency
  xlab("") + ylab("dN/dS") +
  scale_color_manual(labels = c("non-SE","SE"), values = c(nonSE,SE)) + 
  facet_grid(~family.y, switch = "x", scales = "free_x", space = "free_x") + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6), 
        legend.position = c(.9, .9), 
        strip.text.x = element_text(size = 3))
dev.off()

### horizontal ####
# non_env horizontal 
pdf(file=file.path(dir_plot, "non_env_genes_h.pdf")) 
par(mar=c(10, 15, 3, 15)) 

file=file.path(dir_plot, "all_non_env_h.png")
png(file, width = 1000, height = 600)
ggplot(non_env, aes(y=name_manual, 
                    x=`mean.dN/mean.dS`, colour=surface)) + 
  xlab("dN/dS") + ylab("") +
  scale_color_manual(labels = c("non-SE","SE"), values = c(nonSE,SE)) + 
  facet_grid(manual_family ~., scales = "free", space = "free") +
  theme(strip.text.y = element_text(angle = 0),
        axis.text.y = element_text(size = 14,angle = 15),
        legend.text = element_text(size=10, face="bold"))+
  geom_point() #  geom_jitter(width = 0.01)
dev.off()







## 12. summary ##########  
env <- meta[meta$enveloped == TRUE,]
non_env <- meta[meta$enveloped == FALSE,]

# treelength 
summary(env$tree_length)
summary(non_env$tree_length)

# aa entropy
summary(env$aa_entropy_gene)
summary(non_env$aa_entropy_gene)



## 13. ##########  


## 14. ##########  


## 15. ########## 





# FINAL ####
# multiple plots together 

## 1. ####
pdf(file=file.path(dir_plot,"seq_info.pdf"))
par( mfrow= c(1,2) )

plot(h_without, col =mycol1, xlab='log transformed codon length')
plot(h_with, col =  mycol2, add = TRUE)
#add legend
legend('topright', c('with ovlp','without ovlp'), fill=c(mycol1, mycol2))

x <- barplot(table(env_se), 
             col=c(nonSE,SE), 
             border="white", 
             space=0.04, 
             font.axis=2, 
             xlab = "virus",
             names.arg=c("non-enveloped","enveloped"))

legend("topleft",
       fill=c(SE,nonSE),
       legend=c("SE", "non-SE"),
       bty = "n")

dev.off()



# Pie Chart ####

# virus family 
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

## as is 
png(file=file.path(dir_plot,"virus_family_gene_asis.png"),width=800, height=800)
pie(gene_family)
dev.off()

png(file=file.path(dir_plot,"virus_family_virus_asis.png"),width=800, height=800)
pie(virus_family)
dev.off()

# distinctive colors
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
color = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


### (1) frequency of genes from each family 
png(file=file.path(dir_plot,"piechart_gene_family.png"),
    width=800, height=800)
pie(gene_family, cex=0.7, border="grey", main="Virus Families", col=color, cex.main=1.7)
dev.off()

png(file=file.path(dir_plot,"piechart_virus_family.png"),
    width=800, height=800)
pie(virus_family, cex=1.4, border="grey", main="Virus Families", col=color, cex.main=1.7)
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

## DOSEN'T WORK: CANOT FIND FUNCTION PieDonus
library(ggplot2)
library(webr)
library(dplyr)
require(moonBook)

# run nested pie 
Pie.df <- meta.family_freq %>% group_by(manual_family, enveloped) %>% summarise(n = sum(virus_freq))

PieDonus(Pie.df, aes(manual_family, enveloped, count = n))




##### Codon Length  #####

pdf(file=file.path(dir_plot,"Frequency_of_gene_lengths.pdf"))
hist(all.mle$codon, col = "cadetblue3", breaks=20, main = "Frequency of gene length", xlab = "Number of Codons", horizontal = TRUE)
dev.off()

pdf(file=file.path(dir_plot,"Distribution of gene lengths.pdf"))
boxplot(all.mle$codon, pch = 20, col ="cadetblue2", horizontal = TRUE, main = "Distribution of Codon Length")
dev.off()

pdf(file=file.path(dir_plot,"Distribution of gene lengths hist log.pdf"))
hist(log(all.mle$codon,10), col ="cadetblue2", main = "Distribution of Codon Length without ovlp", breaks = 25, xlab = "log transformed codon length")
dev.off()

pdf(file=file.path(dir_plot,"Distribution of original gene lengths uniport hist log.pdf"))
hist(log(meta$Length,10), col ="cadetblue2", main = "Distribution of Original Codon Length", breaks = 25, xlab = "log transformed codon length")
dev.off()

pdf(file=file.path(dir_plot,"Distribution of genome lengths hist log.pdf"))
hist(log(meta$genome_length.nt,10), col ="cadetblue2", main = "Distribution of Genome Length", breaks = 20, xlab = "log transformed codon length")
dev.off()

pdf(file=file.path(dir_plot,"Distribution of gene lengths hist.pdf"))
hist(all.mle$codon, col ="cadetblue2", main = "Distribution of Codon Length", breaks = 100, xlab = "codon length")
dev.off()

######################
# save both histograms 
######################

h_without <- hist(log(meta$codon,10), plot = FALSE, main = "Distribution of Codon Length without ovlp", breaks = 25, xlab = "log transformed codon length")

hist(log(meta$codon,10), main = "Distribution of Codon Length without ovlp", breaks = 25, xlab = "log transformed codon length")

h_with <- hist(log(meta$Length,10), plot = FALSE, main = "Distribution of Original Codon Length", breaks = 25, xlab = "log transformed codon length")

# colours
mycol1 <- rgb(0, 204, 204, max = 255, alpha = 50, names = "blue50")
mycol2 <- rgb(255, 255, 0, max = 255, alpha = 50, names = "blue50")

# plot both hist 
pdf(file=file.path(dir_plot,"with and without ovlp codon distribution3.pdf"))
plot(h_without, col =col_without_ovlp,alpha = 50, xlab='log transformed codon length')
plot(h_with, col = col_with_ovlp, alpha = 50, add = TRUE) 
#add legend
legend('topright', c('with ovlp','without ovlp'), fill=c(col_with_ovlp, col_without_ovlp))
dev.off()


###################
##### checkes #####

pdf(file=file.path(dir_plot,"length_and_dNdS_noovlp.pdf"))
plot(meta$mean.dN.mean.dS,meta$all_codon, pch = 20, cex = 0.8, col="Steel Blue 3", main="Gene length and dN/dS ratio", xlab="dN/dS" , ylab="Number of Codons")
dev.off()

pdf(file=file.path(dir_plot,"length_and_dNdS_withovlp.pdf"))
plot(meta$mean.dN.mean.dS,meta$all_codon, pch = 20, cex = 0.8, col="Steel Blue 3", main="Gene length and dN/dS ratio", xlab="dN/dS" , ylab="Number of Codons")
dev.off()

pdf(file=file.path(dir_plot,"percentag_of_positive_sites_and_length.pdf"))
plot(log(meta$codon,10), meta$pos.sites.codon, pch=20, cex = 0.8,col = "cadetblue4", main = "Percentage of Positive sites and gene length", xlab = "length", ylab = "Percentage of positive sites" )
dev.off()

pdf(file=file.path(dir_plot,"percentag_of_negative_sites_and_length.pdf"))
plot(log(meta$codon,10), meta$neg.sites.codon, pch=20, col = "cadetblue4", main = "Percentage of Positive sites and gene length", xlab = "length", ylab = "percentage of negative sites" )
dev.off()



#############
#### MLE ####
#############

######################
####### dN/dS  #######

# method.1
pdf(file=file.path(dir_plot,"dN.dS_plot1.pdf"))
labels <- c("SE", "nonSE")
boxplot(TRUE.uni$mean.dN.mean.dS, 
        FALSE.uni$mean.dN.mean.dS, 
        names=labels, las = 1, cex.lab=4, 
        cex.axis=1, col="cadetblue", main = "dN/dS",
        boxwex= 0.5, outline = FALSE )
dev.off()  

# method.2 
boundaries <- boxplot(meta$mean.dN.mean.dS ~ meta$surface , col="#69b3a2")
# Draw the boxplot. Note result is also stored in a object called boundaries
boundaries$stats
# Now you can type boundaries$stats to get the boundaries of the boxes

# Add sample size on top
nbGroup <- nlevels(meta$surface)
text( 
  x=c(1:nbGroup), 
  y=boundaries$stats[nrow(boundaries$stats),] + 0.5, 
  paste("n = ",table(meta$surface),sep="")  
)

pdf(file=file.path(dir_plot,"dN.dS_plot2.pdf"))
boxplot(meta$mean.dN.mean.dS ~ meta$surface,
        las = 1, cex.lab=1, cex.axis=1, 
        col="cadetblue", main = "dN/dS",
        boxwex= 0.5, outline = FALSE,
        xlab = "Surface Exposed", ylab = "dN/dS")
text( 
  x=c(1:2), 
  y=boundaries$stats[nrow(boundaries$stats),] - 0.1, 
  paste("n = ",table(meta$surface),sep="")  
)
dev.off()

####### dN & dS  #######

pdf(file=file.path(dir_plot,"mean_dN.pdf"))
labels <- c("mean.dN of SE", "mean.dN of nonSE")
boxplot(TRUE.uni$mean.dN, FALSE.uni$mean.dN, names=labels, las = 1, cex.lab=2, cex.axis=0.8, col="cyan4", main = "Mean dN Values", boxwex= 0.5)
dev.off()

pdf(file=file.path(dir_plot,"mean_dS.pdf"))
labels <- c("mean.dS of SE", "mean.dS of nonSE")
boxplot(TRUE.uni$mean.dS, FALSE.uni$mean.dS, names=labels, las = 1, cex.lab=2, cex.axis=0.8, col="cyan4", main = "Mean dS Values", boxwex= 0.5)
dev.off()

pdf(file=file.path(dir_plot,"mean_dN.dS.pdf"))
labels <- c("mean.dN of SE", "mean.dN of nonSE", "mean.dS of SE", "mean.dS of nonSE")
boxplot(TRUE.uni$mean.dN, FALSE.uni$mean.dN, TRUE.uni$mean.dS, FALSE.uni$mean.dS, names=labels, las = 1, cex.lab=2, cex.axis=0.8, col="cyan4", main = "Mean dN & dS Values", boxwex= 0.5)
dev.off()

pdf(file=file.path(dir_plot,"mean_dN.dS2.pdf"))
labels <- c("mean.dN of SE", "mean.dN of nonSE", "mean.dS of SE", "mean.dS of nonSE")
boxplot(TRUE.uni$mean.dN, FALSE.uni$mean.dN, TRUE.uni$mean.dS, FALSE.uni$mean.dS, names=labels, las = 1, cex.lab=2, cex.axis=0.8, col="cyan4", main = "Mean dN & dS Values", boxwex= 0.5)
dev.off()

pdf(file=file.path(dir_plot,"percent_pos_neg.pdf"))
labels <- c("% pos of SE", "% pos of nonSE", "% neg of SE", "% neg of nonSE")
boxplot(TRUE.uni$pos.sites.codon, FALSE.uni$pos.sites.codon, TRUE.uni$neg.sites.codon, FALSE.uni$neg.sites.codon, names=labels, las = 1, cex.lab=1.5, cex.axis=0.8, col="cyan4", main = "Percentage of positive & negative codons", boxwex= 0.2, outline = FALSE)
dev.off()

#### pos sites ####
pdf(file=file.path(dir_plot,"percent_pos.pdf"))
labels <- c("% pos of SE", "% pos of nonSE")
boxplot(TRUE.uni$pos.sites.codon, FALSE.uni$pos.sites.codon, names=labels, las = 1, cex.lab=1.5, cex.axis=0.8, col="cyan4", main = "Percentage of positive codons", boxwex= 0.2)
dev.off()

meta_nona <- meta[!is.na(meta$surface),]
compare_means(`pos.sites/codon`~surface, meta_nona, method = "wilcox.test", paired = FALSE,group.by = NULL, ref.group = NULL)

p <- ggboxplot(meta_nona, 
               x = "surface", y = "pos.sites/codon", 
               fill= c(SE,nonSE))


stat_box_data <- function(y, upper_limit = max(meta$`pos.sites/codon`) * 1.15) {
  return( 
    data.frame(
      y = 0.0001 * upper_limit,
      label = paste('n =', length(y), '\n'), '\n'))
}

pdf(file=file.path(dir_plot,"pos.sites.codon_wilcoxon5.pdf"))
p + stat_compare_means(label.y = 0.95, size = 6) + geom_jitter(shape=20, position=position_jitter(0.3), alpha=.1) + ylab("pos.sites.codon") + xlab("Surface Exposed") + stat_summary(fun.data = stat_box_data, geom = "text", hjust = 0.5, vjust = 0.8) +  coord_cartesian(ylim = c(0, 1)) +   theme(
  axis.title.y = element_text(vjust = -0.001),
  axis.title.x = element_text(vjust = +0.001)
) 
dev.off()

##### neg sites #####
pdf(file=file.path(dir_plot,"percent_neg.pdf"))
labels <- c("% neg of SE", "% neg of nonSE")
boxplot(TRUE.uni$neg.sites.codon, FALSE.uni$neg.sites.codon, names=labels, las = 1, cex.lab=1.5, cex.axis=0.8, col="cyan4", main = "Percentage of negative codons", boxwex= 0.2)
dev.off()

meta_nona <- meta[!is.na(meta$surface),]
compare_means(`neg.sites/codon`~surface, meta_nona, method = "wilcox.test", paired = FALSE,group.by = NULL, ref.group = NULL)

p <- ggboxplot(meta_nona, 
               x = "surface", y = "neg.sites/codon", 
               fill= c(SE,nonSE))


stat_box_data <- function(y, upper_limit = max(meta$`neg.sites/codon`) * 1.15) {
  return( 
    data.frame(
      y = 0.0001 * upper_limit,
      label = paste('n =', length(y), '\n'), '\n'))
}

pdf(file=file.path(dir_plot,"per_neg_sites_wilcoxon5.pdf"))
p + stat_compare_means(label.y = 0.95, size = 6) + geom_jitter(shape=20, position=position_jitter(0.3), alpha=.1) + ylab("per_neg_sites") + xlab("Surface Exposed") + stat_summary(fun.data = stat_box_data, geom = "text", hjust = 0.5, vjust = 0.8) +  coord_cartesian(ylim = c(0, 1)) +   theme(
  axis.title.y = element_text(vjust = -0.001),
  axis.title.x = element_text(vjust = +0.001)
) 
dev.off()


#### + & neg ####
pdf(file=file.path(dir_plot,"percent_pos&neg.pdf"))
labels <- c("% pos and neg of SE", "% pos and neg of nonSE")
boxplot(TRUE.uni$per_pos_neg_sites, FALSE.uni$per_pos_neg_sites, names=labels, las = 1, cex.lab=1.5, cex.axis=0.8, col="cyan4", main = "Percentage of positive and negativly selected codons", boxwex= 0.2)
dev.off()

meta_nona <- meta[!is.na(meta$surface),]
compare_means(per_pos_neg_sites~surface, meta_nona, method = "wilcox.test", paired = FALSE,group.by = NULL, ref.group = NULL)

p <- ggboxplot(meta_nona, 
               x = "surface", y = "per_pos_neg_sites", 
               fill= c(SE,nonSE))


stat_box_data <- function(y, upper_limit = max(meta$per_pos_neg_sites) * 1.15) {
  return( 
    data.frame(
      y = 0.0001 * upper_limit,
      label = paste('n =', length(y), '\n'), '\n'))
}

pdf(file=file.path(dir_plot,"per_pos_neg_sites_wilcoxon5.pdf"))
p + stat_compare_means(label.y = 0.95, size = 6) + geom_jitter(shape=20, position=position_jitter(0.3), alpha=.1) + ylab("per_pos_neg_sites") + xlab("Surface Exposed") + stat_summary(fun.data = stat_box_data, geom = "text", hjust = 0.5, vjust = 0.8) +  coord_cartesian(ylim = c(0, 1)) +   theme(
  axis.title.y = element_text(vjust = -0.001),
  axis.title.x = element_text(vjust = +0.001)
) 
dev.off()



## env & non-env ####

boundaries <- boxplot(meta$mean.dN.mean.dS ~ meta$enveloped, col="#69b3a2")
# Draw the boxplot. Note result is also stored in a object called boundaries
boundaries$stats
# Now you can type boundaries$stats to get the boundaries of the boxes

# Add sample size on top
nbGroup <- nlevels(meta$enveloped)
text( 
  x=c(1:nbGroup), 
  y=boundaries$stats[nrow(boundaries$stats),] + 0.5, 
  paste("n = ",table(meta$enveloped),sep="")  
)

pdf(file=file.path(dir_plot,"env&nonenv_plot.pdf"))
boxplot(meta$mean.dN.mean.dS ~ meta$enveloped,
        las = 1, cex.lab=1, cex.axis=1, 
        col="cadetblue", main = "dN/dS",
        boxwex= 0.5, outline = FALSE,
        xlab = "Enveloped virus", ylab = "dN/dS")
text( 
  x=c(1:2), 
  y=boundaries$stats[nrow(boundaries$stats),] - 0.1, 
  paste("n = ",table(meta$enveloped),sep="")  
)
dev.off()


###### significance ####

meta_nona <- meta[!is.na(meta$enveloped),]

compare_means(mean.dN.mean.dS~enveloped, meta_nona, method = "wilcox.test", paired = FALSE,group.by = NULL, ref.group = NULL)

p <- ggboxplot(meta_nona, 
               x = "enveloped", y = "mean.dN/mean.dS",
               fill = c(SE, nonSE))

stat_box_data <- function(y, upper_limit = max(meta$`mean.dN/mean.dS`) * 1.15) {
  return( 
    data.frame(
      y = -0.000001 * upper_limit,
      label = paste('n =', length(y), '\n'), '\n'))
}
pdf(file=file.path(dir_plot,"env&nonenv_wilcoxon5.pdf"))
p + stat_compare_means(label.y = 0.95, size = 6) + geom_jitter(shape=20, position=position_jitter(0.3), alpha=.1) + ylab("dN / dS ") + xlab("Enveloped virus") + stat_summary(fun.data = stat_box_data, geom = "text", hjust = 0.5, vjust = 0.8) +  coord_cartesian(ylim = c(0, 1.4)) +   theme(
  axis.title.y = element_text(vjust = -0.001),
  axis.title.x = element_text(vjust = +0.001)
) 
dev.off()







#########################
###### COMPARISON #######
#########################

# Significance ----
### dN/dS ----
# https://statsandr.com/blog/wilcoxon-test-in-r-how-to-compare-2-groups-under-the-non-normality-assumption/

test <- wilcox.test(TRUE.uni$mean.dN.mean.dS, FALSE.uni$mean.dN.mean.dS)
summary(test)

#independent 2-group Mann-Whitney U Test
#where y is numeric and A is A binary factor
wilcox.test(meta$mean.dN.mean.dS ~ meta$surface)

# where y and x are numeric
wilcox.test(TRUE.uni$mean.dN.mean.dS,FALSE.uni$mean.dN.mean.dS)


## ADD TO PLOT 
library(ggplot2)
library(rtatix)
library(ggpubr)

meta_nona <- meta[!is.na(meta$surface),]

compare_means(mean.dN.mean.dS~surface, meta_nona, method = "wilcox.test", paired = FALSE,group.by = NULL, ref.group = NULL)

p <- ggboxplot(meta_nona, 
               x = "surface", y = "mean.dN.mean.dS", 
               fill= c(SE,nonSE))


stat_box_data <- function(y, upper_limit = max(meta$mean.dN.mean.dS) * 1.15) {
  return( 
    data.frame(
      y = 0.0001 * upper_limit,
      label = paste('n =', length(y), '\n'), '\n'))
}

# , 'mean =', round(mean(y), 4

pdf(file=file.path(dir_plot,"dNdS_wilcoxon5.pdf"))
p + stat_compare_means(label.y = 0.95, size = 6) + 
  geom_jitter(shape=20, position=position_jitter(0.3), alpha=.1) + 
  ylab("dN / dS ") + xlab("Surface Exposed") + 
  stat_summary(fun.data = stat_box_data, 
               geom = "text", hjust = 0.5, vjust = 0.8) +  
  theme(
  axis.title.y = element_text(vjust = -0.001),
  axis.title.x = element_text(vjust = +0.001)
) 
dev.off()

#coord_cartesian(ylim = c(0, 1)) +






######## dN ############
meta$surface <- meta$SEgo

meta_nona <- meta[!is.na(meta$surface),]

compare_means(mean.dN~surface, meta_nona, method = "wilcox.test", paired = FALSE,group.by = NULL, ref.group = NULL)

p <- ggboxplot(meta_nona, 
               x = "surface", y = "mean.dN", 
               fill= c(SE,nonSE))


stat_box_data <- function(y, upper_limit = max(meta$mean.dN) * 1.15) {
  return( 
    data.frame(
      y = 0.0001 * upper_limit,
      label = paste('n =', length(y), '\n'), '\n'))
}

# , 'mean =', round(mean(y), 4

pdf(file=file.path(dir_plot,"dN_wilcoxon.pdf"))
p + stat_compare_means(label.y = 2, label.x = 1.25,size = 6) + geom_jitter(shape=20, position=position_jitter(0.3), alpha=.1) + ylab("mean dN") + xlab("Surface Exposed") + stat_summary(fun.data = stat_box_data, geom = "text", hjust = 0.5, vjust = 0.8)  + coord_cartesian(ylim = c(0, 2.5)) + theme(
  axis.title.y = element_text(vjust = -0.001),
  axis.title.x = element_text(vjust = +0.001)
) 
dev.off()



######## dS ############
meta$surface <- meta$SEgo

meta_nona <- meta[!is.na(meta$surface),]

compare_means(mean.dS~surface, meta_nona, method = "wilcox.test", paired = FALSE,group.by = NULL, ref.group = NULL)

p <- ggboxplot(meta_nona, 
               x = "surface", y = "mean.dS", 
               fill= c(SE,nonSE))


stat_box_data <- function(y, upper_limit = max(meta$mean.dS) - 5 ) {
  return( 
    data.frame(
      y = 0.0001 * upper_limit,
      label = paste('n =', length(y), '\n'), '\n'))
}

# , 'mean =', round(mean(y), 4

pdf(file=file.path(dir_plot,"dS_wilcoxon.pdf"))
p + stat_compare_means(label.y = 2, label.x = 1.25,size = 6) + geom_jitter(shape=20, position=position_jitter(0.3), alpha=.1) + ylab("mean dS") + xlab("Surface Exposed") + stat_summary(fun.data = stat_box_data, geom = "text", hjust = 2, vjust = 4)  + coord_cartesian(ylim = c(1, 6)) + theme(
  axis.title.y = element_text(vjust = +0.001),
  axis.title.x = element_text(vjust = +0.001)
) 
dev.off()






#meta_nona <- meta[!is.na(meta$surface),]

compare_means(mean.dN.mean.dS~DNA, meta_nona, method = "wilcox.test", paired = FALSE,group.by = NULL, ref.group = NULL)

p <- ggboxplot(meta_nona, 
               x = "DNA", y = "mean.dN.mean.dS", 
               fill= c(nonSE,SE))

stat_box_data <- function(y, upper_limit = max(meta$mean.dN.mean.dS) * 1.15) {
  return( 
    data.frame(
      y = 0.0001 * upper_limit,
      label = paste('n =', length(y), '\n'), '\n'))
}

pdf(file=file.path(dir_plot,"TEST.pdf"))
p + stat_compare_means(label.y = 0.95, size = 6) + 
  geom_jitter(shape=20, position=position_jitter(0.3), alpha=.1) + 
  ylab("dN / dS ") + xlab(" ") + 
  stat_summary(fun.data = stat_box_data, 
               geom = "text", hjust = 0.5, vjust = 0.95) +  
  theme(
    axis.title.y = element_text(vjust = -0.001),
    axis.title.x = element_text(vjust = +0.001)) +
  scale_x_discrete(labels= c("non-SE","SE"))
dev.off()









########## EXPLORE ##########






########## ANALYSIS ############

# examin SE 
#Examine missmatch between GOID and tmbed

a <- meta[meta$surface == TRUE & meta$surface_tmbed == FALSE,] #4
#View(a)

b <- meta[meta$surface == FALSE & meta$surface_tmbed == TRUE,]
c <- b[!is.na(b$gene),]
#View(c)

SE_missmatch <- rbind(a,c)
View(SE_missmatch)

# tables 
table(meta$surface,useNA ='always')
table(meta$surface_tmbed,useNA ='always')
table(meta$surface, meta$surface_tmbed,useNA ='always')

fisher.test(table(meta$surface, meta$surface_tmbed))


### dN ----
#independent 2-group Mann-Whitney U Test
#where y is numeric and A is A binary factor
wilcox.test(meta$mean.dN ~ meta$surface)

# where y and x are numeric
wilcox.test(TRUE.uni$mean.dN,FALSE.uni$mean.dN) 

### dS ----
#independent 2-group Mann-Whitney U Test
#where y is numeric and A is A binary factor
wilcox.test(meta$mean.dS ~ meta$surface)

# independent 2-group Mann-Whitney U Test
wilcox.test(TRUE.uni$mean.dS,FALSE.uni$mean.dS) # where y and x are numeric

### neg% ----
#% do chi SQUARE/ fisher
#independent 2-group Mann-Whitney U Test
#where y is numeric and A is A binary factor
wilcox.test(meta$`neg.sites/codon` ~ meta$surface)

# independent 2-group Mann-Whitney U Test
wilcox.test(TRUE.uni$`neg.sites/codon`,FALSE.uni$`neg.sites/codon`) # where y and x are numeric

### chi
#input data is in the form of a table that contains the count value of the variables in the observation
data <- table(meta$`neg.sites/codon`,meta$surface)
chisq.test(data)


### pos% ----
#independent 2-group Mann-Whitney U Test
#where y is numeric and A is A binary factor
wilcox.test(meta$`pos.sites/codon`  ~ meta$surface)

# independent 2-group Mann-Whitney U Test
wilcox.test(TRUE.uni$`pos.sites/codon`, FALSE.uni$`pos.sites/codon`) # where y and x are numeric

### chi
#input data is in the form of a table that contains the count value of the variables in the observation
data <- table(meta$`pos.sites/codon`,meta$surface)
chisq.test(data)


### bionomial Regression ----
#response variable needs to be categorical
#"where your dependent variable is just 0s and 1"

fit <- glm(formula = surface ~ mean.dS + mean.dN + `pos.sites/codon` + `neg.sites/codon` + mean.dN.mean.dS,
           data = meta,
           family = binomial)
summary(fit)


### fitted column ----
#Add logistic fitted values back to dataframe as
#new column 
meta$fitted.values <- fit$fitted.values
head(meta[c("fitted.values","surface","mean.dN","mean.dS","mean.dN/mean.dS","pos.sites/codon","neg.sites/codon")])

### dn & ds ----
fit.mean.1 <- glm(formula = surface ~ mean.dS + mean.dN,
                  data = meta,
                  family = binomial)
summary(fit.mean)

fit.mean.2 <- glm(formula = surface ~ mean.dS - mean.dN,
                  data = meta,
                  family = binomial)
summary(fit.mean)

fit.mean.3 <- glm(formula = surface ~ mean.dN - mean.dS,
                  data = meta,
                  family = binomial)
summary(fit.mean)


# % pos/neg 
fit.percent <- glm(formula = surface ~ `pos.sites/codon` + `neg.sites/codon`,
                   data = meta,
                   family = binomial)
summary(fit.percent)

# mean.dN/mean.dS
fit.ratio <- glm(formula = surface ~ mean.dN.mean.dS,
                 data = meta,
                 family = binomial)
summary(fit.ratio)


### AIC ----

aic.all <- meta[,c("Length","surface","mean.dS","mean.dN","pos.sites","neg.sites","mean.dN/mean.dS","pos.sites/codon","neg.sites/codon","Mass","aa_entropy","gap_per","total_tree_length")]

aic.dn.ds <- meta[,c("surface","mean.dS","mean.dN","pos.sites","neg.sites","mean.dN/mean.dS","pos.sites/codon" ,"neg.sites/codon" )]

#dn.ds aic.dn.ds
multiLinearReg = lm(mean.dN.mean.dS ~ ., data = aic.dn.ds, drop.unused.levels = TRUE)
stepAIC(multiLinearReg, direction = "both")





# NOTES ####

# poly issue
# ALL
all_df <- read.csv("/users/sareh/Desktop/all.txt")
View(all_df)
all <- all_df$Accn
print(all, quote=FALSE)

#SAREH
sareh_half_df <- read.csv("/users/sareh/Desktop/sareh.txt")
sareh_half<- sareh_half_df$Accn
print(sareh_half, quote=FALSE)
#NC_000940 NC_001348 NC_001355 NC_001357 NC_001414 NC_001427 NC_001430 NC_001434 NC_001436 NC_001437 NC_001461 NC_001472 NC_001489 NC_001498 NC_001526 NC_001538 NC_001542 NC_001547 NC_001563 NC_001593 NC_001612 NC_001617 NC_001653 NC_001699 NC_001781 NC_001802 NC_002058 NC_002200 NC_002640 NC_003045 NC_003899 NC_003977 NC_004162 NC_004296 NC_004297 NC_004441 NC_004455 NC_004718 NC_005148 NC_005218 NC_005219 NC_005222 NC_005300 NC_005302 NC_006273

#ART
art_half <- setdiff(all,sareh_half)
print(art_half, quote = FALSE)
#NC_005301 NC_006383 NC_007543 NC_007544 NC_007545 NC_007546 NC_007547 NC_007569 NC_007570 NC_007571 NC_007572 NC_007573 NC_007574 NC_009334 NC_009942 NC_010956 NC_011202 NC_011203 NC_011500 NC_011501 NC_011502 NC_011503 NC_011504 NC_011505 NC_011506 NC_011507 NC_011508 NC_011509 NC_011510 NC_012532 NC_014395 NC_014396 NC_014397 NC_015843 NC_018136 NC_018137 NC_018138 NC_019843 NC_030791 NC_035889 NC_038235 NC_038270 NC_038319 NC_038878 NC_039199 NC_039237 NC_039476 NC_048214

a <- meta[!duplicated(meta$virus),]
View(a)
View(fltr.nbr)

require(dplyr)

b <- fltr.nbr %>%
  group_by(Representative) %>%
  summarise_each(funs(toString))

c <- b[,c("Representative","Neighbor")]

path <- file.path(dir_base,"c.tsv")
write.table(c, path, quote=FALSE, sep = ":")


# TESTING ####

### 1. DNA/RNA ####
dna_meta <- meta[meta$DNA=="TRUE",]
rna_meta <- meta[meta$DNA=="FALSE",]

summary(rna_meta$mean.dN.mean.dS)
summary(dna_meta$mean.dN.mean.dS)

#### codon length 
summary(dna_meta$codon_no_ovlp)
summary(rna_meta$codon_no_ovlp)

#  1 codon_no_ovlp TRUE   FALSE  0.00109 0.0011 0.0011   **       Wilcoxon
compare_means(codon_no_ovlp~DNA, meta_nona, method = "wilcox.test", paired = FALSE,group.by = NULL, ref.group = NULL)


# codon_no_ovlp    0.003685   0.001465   2.515   0.0121 * 
# codon_with_ovlp -0.002826   0.001450  -1.950   0.0516 
result.test <- glm(mean.dN.mean.dS ~ DNA + codon_no_ovlp + codon_with_ovlp, 
                   family="Gamma", data=meta)
summary(result.test)


### 1. env/non-env ####
env_meta <- meta[meta$enveloped=="TRUE",]
non_env_meta <- meta[meta$enveloped=="FALSE",]
dim(env_meta)
dim(non_env_meta)

summary(env_meta$mean.dN.mean.dS)
summary(non_env_meta$mean.dN.mean.dS)

#codon_no_ovlp TRUE   FALSE  0.0000314 0.000031 3.1e-05  ****     Wilcoxon
compare_means(codon_no_ovlp~enveloped, meta_nona, method = "wilcox.test", paired = FALSE,group.by = NULL, ref.group = NULL)

#codon_no_ovlp    0.003649   0.001251   2.918  0.00364 ** 
result.test <- glm(mean.dN.mean.dS ~ enveloped + codon_no_ovlp + codon_with_ovlp, 
                   family="Gamma", data=meta)
summary(result.test)



