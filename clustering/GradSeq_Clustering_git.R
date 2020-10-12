###################### Part I - Clustering #############################

require(pdc)
require(reshape2)
require(readxl)

require(WGCNA)
options(stringsAsFactors = FALSE);
enableWGCNAThreads()


RNA <- read.csv("WGCNA_input_RNAseq_features.csv", sep=",")
RNA <- transform(RNA, Feature_type = "RNA")
Protein <- read.csv("WGCNA_input_proteomics.csv", sep=",")
Protein <- transform(Protein, Feature_type = "Protein")

#RNA <- RNA[ -which(RNA$Average_Sum_readcount < 8), ]
#Protein <- Protein[ -which(Protein$Average_Sum_intensities < 10), ]

names(RNA) <- names(Protein)
Input <- rbind(Protein, RNA)




###WGCNA: hierarchical clustering of Features w. dynamic branch cutting ###

x <- cbind( Input[2], round(Input[3:(3+(3*16-1))],2 ))
x_av <- cbind( Input[2], round(Input[(3+(3*16)):(3+(3*16+16-1))],2 ))


scale_free_topology_criteria <- function(x){ ### Analysis if (and at which power) the network satisfies a scale free topology: RB2 >= .95 ###
  powers = c(c(1:10), seq(from = 12, to=20, by=2))  # Choose a set of soft-thresholding powers
  sft = pickSoftThreshold(x, powerVector = powers, verbose = 5) # Call the network topology analysis function
  sizeGrWindow(9, 5); par(mfrow = c(1,2)); cex1 = 0.9; #Plots
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence")); text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red");	abline(h=0.90,col="red") # this line corresponds to using an R^2 cut-off of h
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",	main = paste("Mean connectivity")); text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
}
TOM_heatmap <- function(x, type, geneTree, dynamicColors){### Plot: TOM Heatmap ###
  plotTOM = (1-type); diag(plotTOM) = NA; sizeGrWindow(6,6)
  #pdf("TOMplot_total.pdf", width=12, height=12)
  TOMplot(plotTOM, geneTree, dynamicColors, main = "Network heatmap plot, all genes") ### Adjacency heatmap  
  #dev.off()
}



n <- x$ID; x <- as.data.frame(t(x[,-1])); colnames(x) <- n # Format for WGCNA
x <- scale(x, center = TRUE, scale = TRUE); colnames(x) <- n
n <- x_av$ID; x_av <- as.data.frame(t(x_av[,-1])); colnames(x_av) <- n # Format for WGCNA
x_av <- scale(x_av, center = TRUE, scale = TRUE); colnames(x_av) <- n


sampleTree = hclust(dist(x), method = "average")
sizeGrWindow(12,9);par(cex = 0.6);par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
traitColors = numbers2colors(x_av, signed = TRUE, colors=c("white","white", "white", "grey80", "grey50", "black"))


scale_free_topology_criteria(x)
pow = 16 
MinmodSize = 10	


ADJ=adjacency(x,type="signed", power=pow) 
TOM = TOMsimilarity(ADJ)  
type <- TOM  
geneTree = hclust(as.dist(1-type), method="average") # h-Clustering from dissimilarity adjacency matrix ->TOM MEigengenes


dynamicMods = cutreeDynamic(geneTree, distM = (1-type), deepSplit = 1, pamStage = TRUE, pamRespectsDendro = FALSE, minClusterSize = MinmodSize)	# dynamic tree cutting: based on adjacency matrix (ADJ) or TOM
dynamicColors = labels2colors(dynamicMods)

ModuleNumber <- as.data.frame(dynamicMods) 
y <- cbind(dynamicColors, ModuleNumber)
z <- cbind(Input, y)
#write.csv(z, "WGCNA_full_out.csv")


#pdf("dendrogram.pdf", width=10, height=8)
layout(matrix(c(1,2,3), 3, 1, byrow = TRUE), widths=c(1,1,1), heights=c(3,1,5))
par(mar=c(1,5,1,1) + 0.0) 
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",labels = FALSE, cex = 0.5)
plotColorUnderTree(geneTree, dynamicColors, rowLabels= "Dynamic Tree Cut")
plotColorUnderTree(geneTree, t(traitColors), rowLabels=rownames(x_av))
#dev.off()



