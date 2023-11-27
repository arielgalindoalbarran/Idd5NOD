#install.packages("qtl2", repos="https://rqtl.org/qtl2cran")
#Installing and loading required packages
if(!require(qtl2))
{install.packages("qtl2", repos="http://rqtl.org/qtl2cran")
}
library(qtl2)
##########
#Working directory:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#Write name of phenotype:
phenotype_name="Treg"
#Write name of physical map file:
pmap_name="pmap_genotypecorrected_mm10_cM.csv"
#Write name of genomic map file:
gmap_name= "gmap_genotypecorrected_cM.csv"

##########
#Read yaml file
df <- read_cross2("Treg_bc_mm10_correctedgenotype_cM.yaml")

#Process data
map <- insert_pseudomarkers(df$gmap)
pr <- calc_genoprob(df, map, error_prob=0.002, cores=4)
apr <- genoprob_to_alleleprob(pr)

#Run analysis
out <- scan1(pr, df$pheno, cores=4)


###############################
#Process data for Graphics
###############################
pmap <- read.csv(file=pmap_name, header=TRUE, sep = ",")
LOD <- as.data.frame(out)
#sort 
LOD$id <- rownames(LOD)
LOD <- LOD[order(LOD$id), ]
pmap <- pmap[order(pmap$id), ]
### Add data
LOD$id  <- as.character(pmap$id)
LOD$chr   <- as.character(pmap$chr)
LOD$Position_mm10   <- as.numeric(as.character(pmap$Position_mm10)) 
rownames(LOD) <- NULL

#Re-order
LOD <- LOD[,c(2,1,3,4)]
LOD <- LOD[order(LOD$chr, LOD$Position_mm10), ]
rownames(LOD) <- NULL

#Reorder the list
partial_list <- LOD[LOD[,"chr"]==0,]
for (i in 1:19){
chr_list <- LOD[LOD[,"chr"]==i,]
partial_list <- rbind(partial_list, chr_list )
}
chr_list <- LOD[LOD[,"chr"]=="X",]
partial_list <- rbind(partial_list, chr_list )
chr_list <- LOD[LOD[,"chr"]=="Y",]
partial_list <- rbind(partial_list, chr_list )
LOD <- partial_list
rm(partial_list, chr_list)

#Add sequential column
LOD$sequential <- NA
LOD[1,5] <- LOD[1,4]
LOD[2,5] <- LOD[2,4]
for (i in (seq(from=3,to=nrow(LOD)))) {
  LOD[i,5] <- LOD[i,4]+LOD[(i-1),5]
}
rm(i)
###############################
#      Graphics
###############################

#Graphics by chromosome:

dir.create("graphics")

for (i in 1:19){
  pdf_name <- paste("./graphics/", "chromosome_",i,".pdf", sep="" )
  cairo_pdf(pdf_name, width = 11.69, height = 8.27, bg = "transparent")
  chr <- LOD[LOD[,"chr"]==i,]
  plot(chr[,4], chr[,2], col = "white", xlab = "Position Mpb", ylab = "LOD score", 
       main=bquote("chromosome "  ~  .(i)  ), xlim = c(7.6, (max(chr[,4]))-6) , ylim = c(0.15,max(chr[,2])), xaxt='n', yaxt='n')
  lines(chr[,4], chr[,2], lwd = 2)
  abline(a = NULL, b = NULL, h = 3, v = NULL, col = "lightgray", lty = "dotted", lwd =3)
  axis(1, at = seq(0, 200, by = 20), las=1)
  axis(2, at = seq(0, 4, by = 1), las=1)
  dev.off()
}
rm(i, chr)

#Chromosome X
pdf_name <- paste("./graphics/", "chromosome_X",".pdf", sep="" )
cairo_pdf(pdf_name, width = 11.69, height = 8.27, bg = "transparent")
chr <- LOD[LOD[,"chr"]=="X",]
plot(chr[,4], chr[,2], col = "white", xlab = "Position Mpb", ylab = "LOD score", 
     main="chromosome X", xlim = c(7.6, (max(chr[,4]))-6) , ylim = c(0.15,max(chr[,2])), xaxt='n', yaxt='n')
lines(chr[,4], chr[,2], lwd = 2)
abline(a = NULL, b = NULL, h = 3, v = NULL, col = "lightgray", lty = "dotted", lwd =3)
axis(1, at = seq(0, 200, by = 20), las=1)
axis(2, at = seq(0, 4, by = 1), las=1)
dev.off()
rm(chr)

#Chromosome Y
pdf_name <- paste("./graphics/", "chromosome_Y",".pdf", sep="" )
cairo_pdf(pdf_name, width = 11.69, height = 8.27, bg = "transparent")
chr <- LOD[LOD[,"chr"]=="Y",]
plot(chr[,4], chr[,2], col = "white", xlab = "Position Mpb", ylab = "LOD score", 
     main="chromosome Y", xlim = c(7.6, (max(chr[,4]))-6) , ylim = c(0.15,max(chr[,2])), xaxt='n', yaxt='n')
lines(chr[,4], chr[,2], lwd = 2)
abline(a = NULL, b = NULL, h = 3, v = NULL, col = "lightgray", lty = "dotted", lwd =3)
axis(1, at = seq(0, 200, by = 20), las=1)
axis(2, at = seq(0, 4, by = 1), las=1)
dev.off()
rm(chr)

#Chromosome Mt
pdf_name <- paste("./graphics/", "chromosome_Mt",".pdf", sep="" )
cairo_pdf(pdf_name, width = 11.69, height = 8.27, bg = "transparent")
chr <- LOD[LOD[,"chr"]=="0",]
plot(chr[,4], chr[,2], col = "white", xlab = "Position Mpb", ylab = "LOD score", 
     main="chromosome Mt", xlim = c(7.6, (max(chr[,4]))-6) , ylim = c(0.15,max(chr[,2])), xaxt='n', yaxt='n')
lines(chr[,4], chr[,2], lwd = 2)
abline(a = NULL, b = NULL, h = 3, v = NULL, col = "lightgray", lty = "dotted", lwd =3)
axis(1, at = seq(0, 200, by = 20), las=1)
axis(2, at = seq(0, 4, by = 1), las=1)
dev.off()
rm(chr)


##### ALL Chromosomes

#Add info to chromosome Y:
LOD_seq_plusY <- c(LOD$sequential,(max(LOD$sequential)+1000))
LOD_triple_hi_plusY <- c(LOD[,2], "0.001")
#Graph
pdf_name <- paste("./graphics/", "All_chromosomes",".pdf", sep="" )
cairo_pdf(pdf_name, width = 11.69, height = 8.27, bg = "transparent")
plot(LOD_seq_plusY, LOD_triple_hi_plusY, col = "white", xlab = "Chromosome", ylab = "LOD score", 
     main="LOD score in all chromosomes", axes=FALSE)
lines(LOD_seq_plusY, LOD_triple_hi_plusY, lwd = 3)
abline(a = NULL, b = NULL, h = 3, v = NULL, col = "lightgray", lty = "dotted", lwd =3)
# Make x axis
labs <- c("Mt")
labs <- c(labs,seq(from=1, to=19))
labs <- c(labs, "X")
labs <- c(labs,"Y")
### AXIS names
position_median <- c() 
for (i in 0:19){
  mTEMP <-as.numeric(LOD[LOD[,3]==i,5])
  medianchr <- ((max(mTEMP)-min(mTEMP))/2)+min(mTEMP)
  position_median <- c(position_median, medianchr)
}
#X
mTEMP <-as.numeric(LOD[LOD[,3]=="X",5])
medianchr <- ((max(mTEMP)-min(mTEMP))/2)+min(mTEMP)
position_median <- c(position_median, medianchr)
#Y
mTEMP <-as.numeric((max(LOD$sequential)+500))
position_median <- c(position_median, mTEMP)
#
Axis(side=1, at= position_median, labels=labs, pos=0, las=2, lwd=2, cex.axis=1)
Axis(side=2, at= 0:4, labels=c("", "1", "2", "3", "4"), pos=0, las=1, lwd=2)
dev.off()

#######################
### Add IDD5
#######################
pdf_name <- paste("./graphics/", "chromosome_1_IDD5",".pdf", sep="" )
pdf(pdf_name, width = 11.69, height = 8.27, bg = "transparent")
chr <- LOD[LOD[,"chr"]==1,]
plot(chr[,4], chr[,2], col = "white", xlab = "Position Mpb", ylab = "LOD score", 
     main=bquote("chromosome "  ~  .(1)  ), xlim = c(7.6, (max(chr[,4]))-6) , ylim = c(0.15,3.89), xaxt='n', yaxt='n')
lines(chr[,4], chr[,2], lwd = 2)
abline(a = NULL, b = NULL, h = 3, v = NULL, col = "lightgray", lty = "dotted", lwd =3)
axis(1, at = seq(0, 200, by = 20), las=1)
axis(2, at = seq(0, 4, by = 1), las=1)
#segments(x0=60.695 , y0= 3.15, x1 = 61.118, y1 = 3.15, col = "blue", lty = "solid", lwd = 6)
rect(xleft=60.695, ybottom=3.1, xright=61.118, ytop=3.2, angle = 45, col = "royalblue", lty = "solid", lwd = 5, border="royalblue")
text(60.9, 2.8, "IDD5.1", cex = 0.8)
dev.off()

pdf_name <- paste("./graphics/", "chromosome_1_IDD5_2",".pdf", sep="" )
pdf(pdf_name, width = 11.69, height = 8.27, bg = "transparent")
chr <- LOD[LOD[,"chr"]==1,]
plot(chr[,4], chr[,2], col = "white", xlab = "Position Mpb", ylab = "LOD score", 
     main=bquote("chromosome "  ~  .(1)  ), xlim = c(7.6, (max(chr[,4]))-6) , ylim = c(0.15,3.89), xaxt='n', yaxt='n')
lines(chr[,4], chr[,2], lwd = 2)
abline(a = NULL, b = NULL, h = 3, v = NULL, col = "lightgray", lty = "dotted", lwd =3)
axis(1, at = seq(0, 200, by = 20), las=1)
axis(2, at = seq(0, 4, by = 1), las=1)
rect(xleft=60.695, ybottom=0.1, xright=61.118, ytop=0.3, angle = 45, col = "royalblue", lty = "solid", lwd = 5, border="royalblue")
text(60.9, 0.5, "IDD5.1", cex = 1)
dev.off()


##############################
#Write the table:
##############################
LOD_table <- LOD
#gmap_cM <- read.csv(file="gmap_genotypecorrected_cM.csv", header=TRUE, sep = ",")
#LOD_table$Position_cM <- gmap_cM$Adress
LOD_table <- LOD_table[,c(1,3,4,2)]
colnames(LOD_table) <- c("id", "chr", "Position mm10 (Mbp)", "LOD")

phenotypes <- read.csv(file="Treg_pheno.csv", header=TRUE, sep = ",")



write.table(LOD_table, file = "LODscore_higher3.csv", sep = ",", row.names = FALSE)


###############################
#Phenotype/Genotype 
###############################
top_lod <- LOD_table[LOD[,2]>3,]
rownames(top_lod) <- NULL
top_lod$sequential <- NULL
top_lod <- as.matrix(top_lod)

### Graphics
dir.create("pheno_geno_graphics")
for(i in c(1:nrow(top_lod))) {
  cM_pos= as.numeric(top_lod[i, 3])
  SNP_name= as.vector(top_lod[i,1])
  chrname= as.numeric(top_lod[i,2])
  pdf_name <- paste("./pheno_geno_graphics/", "SNP_",SNP_name,".pdf", sep="" )
  cairo_pdf(pdf_name, width = 11.69, height = 8.27, bg = "transparent")
  g <- maxmarg(pr, map, chr=chrname, pos=cM_pos, return_char=TRUE)  
  par(mar=c(4.1, 4.1, 4.1, 2))
  plot_pxg(g, df$pheno, main= SNP_name, 
           swap_axes=TRUE, xlab=paste(phenotype_name, " %", sep=""))
  dev.off()
   }
rm(i, cM_pos, SNP_name, chrname, pdf_name, g)




