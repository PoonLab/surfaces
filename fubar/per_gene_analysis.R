library(ggplot2)

# INFO ####
## plot dir ####
plot_dir = "/users/sareh/Desktop/surface/plot"
mle_dir = file.path(plot_dir,"mle")
pca_dir = file.path(plot_dir,"pca")
virus_dir = file.path(plot_dir,"each_virus")

## meta ####
meta_path <- "/Users/sareh/Desktop/surface/meta_info/meta_output/meta.csv"
meta <-  read.csv(meta_path, sep = "\t")
View(meta)

# pick which label to use for SE 
meta$surface <- meta$SE

## colours ####
SE <- "#FFAB91" # red 
nonSE  <- "#80CBC4" # blue 

# for each virus
length(unique(meta$virus.x)) #84
length(unique(meta$name.x)) #54   
length(unique(meta$name.y)) #55

meta$name.x <- as.factor(meta$name.x)
View(table(meta$name.x))

# sort by virus family 
meta[order(meta$family.y),]


# One Virus ####

## example: HIV ####
virus_accn <- "NC_001802"
name = "Human immunodeficiency virus 1"
# create specific df for virus 
virus_meta <- meta[meta$virus.x == virus_accn,]
View(virus_meta)

# double check SE label 
View(virus_meta[,c("gene", "mean.dN/mean.dS",
                   "manual_name","SE","surface")])
#virus_meta$correcSE <- virus_meta$surface
# correct if needed 
virus_meta$correcSE <- c("FALSE","FALSE","FALSE",
                         "TRUE","TRUE")
virus_meta$manual_name <- c(
  "gag gene protein p17",
  "capsid Gag protein p24",
  "nucleocapsid",
  "Envelope surface glycoprotein gp120",
  "Envelope transmembrane glycoprotein gp41"
)
    
## scatter plots ####
plot_name = paste0(virus_accn , "_plot.pdf")
pdf(file=file.path(virus_dir, plot_name))      

#par(mar=c(8, 4.1, 4.1, 2.1))    
par(mar=c(10, 7, 5, 7)) 
p <- plot(virus_meta$`mean.dN/mean.dS`, 
     col = ifelse(virus_meta$correcSE, SE, nonSE), 
     main= name, cex.main = 1.2,
     xlab= "", ylab="dN/dS", cex.axis= 1.2, 
     cex.lab = 1, cex = 3,
     pch = 16, ylim = c(0.5,3),
     axes=FALSE)

# gene accn 
axis(1, at=1:5, labels=virus_meta$gene, 
     las=1, cex.axis= 0.6, mgp=c(3,0.5,0))

# gene name 
axis(1, at=1:5, labels=FALSE, las=2,
     cex.axis= 0.5,srt = 35, mgp=c(3,2,0))
text(1:nrow(virus_meta), par("usr")[3] - 0.5,
     virus_meta$manual_name, 
     srt = 30, adj = 0.8, 
     xpd = TRUE, cex = 0.7)

axis(2, cex.axis= 0.5)
mtext(text="Genes", side=1, line=7.9)

dev.off() 


## bar plots ####
plot_name = paste0(virus_accn , "_barplot.pdf")
pdf(file=file.path(plot_dir, plot_name))      

par(mar=c(10, 5, 3, 5)) 
p <- barplot(virus_meta$`mean.dN/mean.dS`, 
          col = ifelse(virus_meta$correcSE, SE, nonSE), 
          main= name, cex.main = 1.2,
          xlab= "", ylab="dN/dS", cex.axis= 1.2, 
          cex.lab = 1,
          pch = 16, ylim = c(0,3),
          axes=FALSE)
# gene accn 
axis(1, at=1:5, labels=virus_meta$gene, las=1, 
     cex.axis= 0.6, mgp=c(3,0.5,0), pos = 0.000001)

# gene name 
axis(1, at=1:5, labels=FALSE, las=2,
     cex.axis= 0.5,srt = 35, mgp=c(3,2,0), pos = 0.000001)
text(1:nrow(virus_meta), par("usr")[3] - 0.4,
     virus_meta$manual_name, 
     srt = 30, adj = 0.8, 
     xpd = TRUE, cex = 0.7)

mtext(text="Genes", side=1, line=7.9)

dev.off() 


# LOOP ####
viruses = factor(unique(meta$name.y))
for (index in 1:length(viruses)) {
  name <- viruses[index]
  virus_meta <- meta[meta$name.y == name,]
  
  ### scatter plots ##
  plot_name = paste0(name , "_plot.pdf")
  pdf(file=file.path(plot_dir, "plot", plot_name))      
  #par(mar=c(8, 4.1, 4.1, 2.1))    
  par(mar=c(10, 7, 5, 7)) 
  p <- plot(virus_meta$`mean.dN/mean.dS`, 
            col = ifelse(virus_meta$surface, SE, nonSE), 
            main= name, cex.main = 1.2,
            xlab= "", ylab="dN/dS", cex.axis= 1.2, 
            cex.lab = 1, cex = 1.5,
            pch = 16,
            axes=FALSE)
  # gene accn 
  axis(1, at=1:nrow(virus_meta), labels=virus_meta$gene, 
       las=2, cex.axis= 0.6, mgp=c(3,0.5,0))
  axis(2, cex.axis= 0.5)
  mtext(text="Genes", side=1, line=7.9)
  dev.off()
  
  ### bar plots ##
  plot_name = paste0(name , "_barplot.pdf")
  pdf(file=file.path(plot_dir, "barplot", plot_name))      
  par(mar=c(10, 5, 3, 5)) 
  barplot(virus_meta$`mean.dN/mean.dS`, 
               col = ifelse(virus_meta$surface, SE, nonSE), 
               main= name, cex.main = 1.2,
               xlab= "", ylab="dN/dS", cex.axis= 1.2, 
               cex.lab = 1,
               pch = 16,
               axes=FALSE)
  # gene accn 
  axis(1, at=1:nrow(virus_meta), 
       labels=virus_meta$gene, las=1, 
       cex.axis= 0.6, mgp=c(3,0.5,0), pos = 0.000001)
  mtext(text="Genes", side=1, line=7.9)
  dev.off() 
} 


# ALL ####

[] DARKER COLOURS
[] TRANSPARENT COLOURS 

## all ####
pdf(file=file.path(plot_dir, "all.pdf")) 

file=file.path(plot_dir, "all.png")
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


## horizontal ####

pdf(file=file.path(plot_dir, "all_horizontal.pdf")) 

file=file.path(plot_dir, "all_h.png")
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
pdf(file=file.path(plot_dir, "all_env_genes.pdf")) 
par(mar=c(10, 15, 3, 15)) 

file=file.path(plot_dir, "all_env.png")
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
pdf(file=file.path(plot_dir, "env_genes_h.pdf")) 
par(mar=c(10, 15, 3, 15)) 

file=file.path(plot_dir, "all_env_h.png")
png(file, width = 1000, height = 600)
ggplot(env, aes(y=name.x, x=mean.dN.mean.dS, colour=surface)) + 
  geom_point(size=1.8) + #  add alpha = 0.6 for transparency
  xlab("dN/dS") + ylab("") +
  scale_color_manual(labels = c("non-SE","SE"), values = c(nonSE,SE)) + 
  facet_grid(family.y ~., scales = "free", space = "free") +
  theme(strip.text.y = element_text(angle = 0)) +
  geom_vline(xintercept = 1,col = "navy", lwd = 1)
dev.off()


## non-env ####
non_env <- meta[meta$enveloped == FALSE,]
non_env$surface <- non_env$SE

pdf(file=file.path(plot_dir, "all_non_env_genes.pdf")) 
par(mar=c(10, 15, 3, 15)) 

file=file.path(plot_dir, "all_non_env.png")
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
pdf(file=file.path(plot_dir, "non_env_genes_h.pdf")) 
par(mar=c(10, 15, 3, 15)) 

file=file.path(plot_dir, "all_non_env_h.png")
png(file, width = 1000, height = 600)
ggplot(non_env, aes(y=name.x, x=mean.dN.mean.dS, colour=surface)) + 
  geom_point(size=3) + #  add alpha = 0.6 for transparency
  xlab("Virus") + ylab("dN/dS") +
  scale_color_manual(labels = c("non-SE","SE"), values = c(nonSE,SE)) + 
  facet_grid(family.y ~., scales = "free", space = "free") +
  theme(strip.text.y = element_text(angle = 0))
dev.off()






















