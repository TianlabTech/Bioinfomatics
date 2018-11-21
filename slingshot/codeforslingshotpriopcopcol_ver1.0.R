#[may be updated 11.21]refer to http://www.bioconductor.org/packages/devel/bioc/vignettes/slingshot/inst/doc/slingshot.html#constructing-smooth-curves-and-ordering-cells
#refer to https://satijalab.org/seurat/pbmc3k_tutorial.html

#library
library(Seurat)
library(dplyr)
library(slingshot)
library(RColorBrewer)


pbmc.data=read.table("/home/fyh/Desktop/slingshot/run1978_lane12_normalized_P_only OPC_tscan_marker_delete1290.txt",sep = "\t",header = TRUE,row.names = 1)
pbmc<-CreateSeuratObject(raw.data = pbmc.data,project = "10X_PBMC")
#Normalizing the data
pbmc<-NormalizeData(object = pbmc,normalization.method = "LogNormalize",scale.factor = 10000)
#Detection of variable genes across the single cells
pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
                                        x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
#Scaling the data
pbmc <- ScaleData(object = pbmc)
#Perform linear dimensional reduction
pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
#Cluster the cells
pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:4, resolution = 0.2, print.output = 0, save.SNN = TRUE,force.recalc=TRUE)


#slingshot

t1   <- PCAPlot(object = pbmc, dim.1 = 1, dim.2 = 6)    #dim for different PC, also for xlab & ylab
t2   <- cbind(t1$data$x,t1$data$y)    #x,y for coordinates
rownames(t2) <- pbmc@cell.names    #give cellnames as rownames
g1=getLineages(t2,t1$data$ident,start.clus = 1)    #ident for cluster tag

#result for line31
#class: SlingshotDataSet 

# Samples Dimensions
#    1221          2

#lineages: 2 
#Lineage1: 1  0  
#Lineage2: 1  2  

#curves: 0 


crv<-getCurves(g1)
crv

#result for line47
#class: SlingshotDataSet 

# Samples Dimensions
#    1221          2

#lineages: 2 
#Lineage1: 1  0  
#Lineage2: 1  2  

#curves: 2 
#Curve1: Length: 18.115	Samples: 1076.08
#Curve2: Length: 22.374	Samples: 623.81

#plot cells-dot
plot(t2,pch=16,col = brewer.pal(3,"Set1")[t1$data$ident],xlab='PC1',ylab='PC2')
#plot the lineage curve
lines(crv)
#get pseudotime
pseudo_order<-data.frame(slingPseudotime(crv))
#dataframe——col1:cellnames ; col2:pseudotime
crv1_preorder<-data.frame(cbind(row.names(pseudo_order),pseudo_order$curve1))
crv2_preorder<-data.frame(cbind(row.names(pseudo_order),pseudo_order$curve2))
#sort ascending order based on col2
crv1_order<-crv1_preorder[order(crv1_preorder$X2),]
crv2_order<-crv2_preorder[order(crv2_preorder$X2),]
#write dataframe on disk
write.table(crv1_order,file = "/home/fyh/Desktop/slingshot/crv1_order.txt",sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(crv2_order,file = "/home/fyh/Desktop/slingshot/crv2_order.txt",sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE)
