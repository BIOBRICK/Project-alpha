```R
setwd("C://Users//Administrator//Desktop//newscRNA//GSE102066")
options(stringsAsFactors = F)
exprrt<-read.csv(file = "GSE102066_expresion_matrix.csv",row.names = 1,header = T)
exprentropy<-read.csv(file = "signal_gene_entropy.csv",row.names = 1,header = T)
library(Seurat)
library(sscClust)
library(dplyr)
library(Cairo)
pbmc <- CreateSeuratObject(counts = exprrt, project = "gse102066", min.cells = 3, min.features = 200)
pbmc
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 1500)
top10 <- head(VariableFeatures(pbmc), 10)
VariableFeaturePlot(pbmc)
LabelPoints(plot = plot1, points = top10, repel = TRUE)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
print(pbmc[["pca"]], dims = 1:10, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 100, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:5, cells = 100, balanced = TRUE)
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:10)
JackStrawPlot(pbmc, dims = 1:10)
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.86)
head(Idents(pbmc))
pbmc <- RunTSNE(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "tsne",pt.size = 1.5)
cellcluster<-as.data.frame(Idents(pbmc))
table(cellcluster)
cellposition<-as.data.frame(pbmc@reductions[["tsne"]]@cell.embeddings)
write.csv(cellposition,file = "cellposition1.csv")
write.csv(cellcluster,file = "cellcluster1.csv")
###########################################################################################
library(monocle)
monocleexpr<-read.csv(file='GSE102066_expresion_matrix.csv',header = T)
phenodata<-read.csv(file = "pheno.csv",header = T,row.names = 1)
genedata<-read.csv(file='feature.csv',header = T )
ensembl<-as.data.frame(genedata$ensembl)
names(ensembl)<-'ensembl'
monocleexpr<-merge(ensembl,monocleexpr,by = 'ensembl')
row.names(monocleexpr)<-monocleexpr$ensembl

row.names(genedata)<-genedata$ensembl
genedata<-genedata[,-1]

monocleexpr<-monocleexpr[,-1]
monocleexpr<-as.matrix(monocleexpr)
pd <- new("AnnotatedDataFrame", data = phenodata)
fd <- new("AnnotatedDataFrame", data = genedata)
cds <- newCellDataSet(monocleexpr, phenoData = pd, featureData = fd,expressionFamily=negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.03)
print(head(fData(cds)))
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10))

disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.2)
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
plot_ordering_genes(cds)
plot_pc_variance_explained(cds, return_all = F)
cds <- reduceDimension(cds, max_components = 2, num_dim = 10,
                      reduction_method = 'tSNE', verbose = T)
cds <- clusterCells(cds, num_clusters = 7)
plot1<-plot_cell_clusters(cds, 1, 2, color = "Age",cell_size = 1.5)
plot2<-plot_cell_clusters(cds, 1, 2,cell_size = 1.5)
CairoPDF(file = 'expression.pdf',width =12,height = 8)
CombinePlots(plots = list(plot1, plot2))
dev.off()
cds_expressed_genes <-  row.names(subset(fData(cds),
                                   num_cells_expressed >= 10))
clustering_DEG_genes <-
  differentialGeneTest(cds[cds_expressed_genes,],
                       fullModelFormulaStr = '~Cluster',
                       cores = 1)
diff_test_res <- differentialGeneTest(cds[cds_expressed_genes,],
                                      fullModelFormulaStr = "~Age")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.001))
cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds)
cds <- reduceDimension(cds, max_components = 2,
                            method = 'DDRTree')
cds <- orderCells(cds)
CairoPDF(file = 'trajectory_expression.pdf',width = 10,height = 8)
plot_cell_trajectory(cds, color_by = "Age",cell_size = 1.5)
dev.off()
###########################################################################################
cell3d<-read.csv(file = 'cellposition_expression.csv',header = T)
library(plotly)
cell3d$Age <- as.factor(cell3d$Age)
cell3d$idents <- as.factor(cell3d$idents)
plot_ly(cell3d, x = ~tSNE_1, y = ~tSNE_2, z = ~entropy, color = ~Age,size = 45) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'tSNE_1'),
                      yaxis = list(title = 'tSNE_2'),
                      zaxis = list(title = 'Entropy')))

```





```R
library(ggpubr)
library(sscClust)
library(SingleCellExperiment)
library(Cairo)
########################################################################
setwd("C://Users//Administrator//Desktop//cell_analysis//c405")
options(stringsAsFactors = F)

exp.data <-read.table("c405_expression_matrix.csv",header = T,check.names = F,sep = ',')
row.names(exp.data)<-exp.data$geneSymbol
entropy.data <-read.table("c405_signal_gene_entropy.csv",header = T,sep = ',')
row.names(entropy.data)<-entropy.data$geneSymbol
cell.metadata<-read.table("cell_type.csv",header = T,row.names = 1,sep = ',')
genesymbol<-as.data.frame(exp.data$geneSymbol)
names(genesymbol)<-'geneSymbol'

sce <- ssc.run(sce, subsampling=F, k.batch=5,seed = 9975)
sce<-ssc.build(exp.data[,-1],display.name = exp.data$geneSymbol)
colData(sce)<-DataFrame(cell.metadata)

sce<-ssc.run(sce,method.clust = "SNN", SNN.k=10,SNN.method="eigen",seed=9997)
sce<-ssc.run(sce,method.reduction = "pca",method.clust = "SNN", SNN.k=30,SNN.method="eigen",seed=9997)
#sce<-ssc.run(sce,method.clust = 'kmeans',k.batch = 11,seed = 9997)
#sce<-ssc.run(sce,method.reduction = 'pca',method.clust = 'kmeans',k.batch = 11,seed = 9997)
#sce<-ssc.run(sce,method.reduction = 'pca',method.clust = 'kmeans',k.batch = 5,seed = 9997)


ssc.plot.tsne(sce, 
              gene = NULL,
              columns = ("pca.SNN.kauto"),
              reduced.name = "pca.tsne",p.ncol = 2)


p1<-ssc.plot.tsne(sce,columns = 'pca.SNN.kauto',reduced.name = 'pca.tsne')
p2<-ssc.plot.tsne(sce,columns = 'pca.kmeans.k11',splitBy = 'cell_type')
CairoPDF(file = 'c405_expression_cluster.pdf',height = 8,width = 13)
ggarrange(p1,p2,ncol=2,nrow=1)
dev.off()

p3<-ssc.plot.tsne(sce,columns = 'pca.kmeans.k5',splitBy = "time_point")
p4<-ssc.plot.tsne(sce,columns = 'pca.kmeans.k5',splitBy = "tippint_point")
CairoPDF("c405_cell_expression_time_series.pdf",height = 10,width = 12)
ggarrange(p3,p4,ncol=2,nrow=1)
dev.off()

sce_entropy<-ssc.build(entropy.data[,-1],display.name = entropy.data$geneSymbol)
colData(sce_entropy)<-DataFrame(cell.metadata)

sce_entropy <- ssc.run(sce_entropy, subsampling=F, k.batch=5,seed = 9975)

sce_entropy <- ssc.run(sce_entropy,method.reduction = "pca",
                   method.clust = "SNN", SNN.k=10,SNN.method="eigen",seed=9997)

ssc.plot.tsne(sce_entropy, 
              gene = NULL,
              columns = ("pca.SNN.kauto"),
              reduced.name = "pca.tsne",p.ncol = 2)

p<-ssc.plot.tsne(sce_entropy, 
              gene = NULL,
              columns = c("tipping_point"),
              reduced.name = "pca.tsne",p.ncol = 2)

pp<-ssc.plot.tsne(sce_entropy, 
                  gene = NULL,
                  columns = c("cell_type"),
                  reduced.name = "pca.tsne",p.ncol = 2)

ggarrange(p,pp,ncol = 2)



sce_entropy <- ssc.run(sce_entropy,method.reduction = 'pca',method.clust = "SNN", SNN.k=10,SNN.method="eigen",seed=9997)

sce_entropy<-ssc.run(sce_entropy,method.clust = 'kmeans',k.batch = 3 ,seed = 9997)
sce_entropy<-ssc.run(sce_entropy,method.reduction = 'pca',method.clust = 'kmeans',k.batch = 11 ,seed = 9997)
sce_entropy<-ssc.run(sce_entropy,method.reduction = 'pca',method.clust = 'kmeans',k.batch = 5 ,seed = 9997)
sce_entropy<-ssc.run(sce_entropy,method.reduction = 'pca',method.clust = 'kmeans',k.batch = 3 ,seed = 9997)

p5<-ssc.plot.tsne(sce_entropy,columns = 'pca.kmeans.k11')
p6<-ssc.plot.tsne(sce_entropy,columns = 'pca.kmeans.k11',splitBy = 'cell_type')
CairoPDF("c405_entropy_cell_cluster.pdf",height = ,width = 13)
ggarrange(p5,p6,ncol=2,nrow=1)
dev.off()

p7<-ssc.plot.tsne(sce_entropy,columns = 'pca.kmeans.k5',splitBy = 'time_point')
p8<-ssc.plot.tsne(sce_entropy,columns = 'pca.kmeans.k5',splitBy = 'tippint_point')
CairoPDF("c405_entropy_time_series.pdf",height = 10,width = 12)
ggarrange(p7,p8,ncol=2,nrow=1)
dev.off()

signalgene<-read.table("c405_signal_gene.csv",header = T,sep = ',')
signalgene<-as.data.frame(signalgene[,-2])
names(signalgene)<-'geneSymbol'
signalgene_expr<-merge(signalgene,exp.data,by = 'geneSymbol')
row.names(signalgene_expr)<-signalgene_expr$geneSymbol

sce_signal<-ssc.build(signalgene_expr[,-1],display.name = signalgene_expr$geneSymbol)
colData(sce_signal)<-DataFrame(cell.metadata)

sce_signal<-ssc.run(sce_signal,method.clust = 'kmeans',k.batch = 11,seed = 9997)
sce_signal<-ssc.run(sce_signal,method.reduction = 'pca',method.clust = 'kmeans',k.batch = 11,seed = 9997)
sce_signal<-ssc.run(sce_signal,method.reduction = 'pca',method.clust = 'kmeans',k.batch = 5,seed = 9997)


p9<-ssc.plot.tsne(sce_signal,columns = 'pca.kmeans.k11')
p10<-ssc.plot.tsne(sce_signal,columns = 'pca.kmeans.k11',splitBy = 'cell_type')
CairoPDF(file = 'c405_signal_expression_cluster.pdf',height = 8,width = 13)
ggarrange(p9,p10,ncol=2,nrow=1)
dev.off()

p11<-ssc.plot.tsne(sce_signal,columns = 'pca.kmeans.k5',splitBy = "cell_type")
p12<-ssc.plot.tsne(sce_signal,columns = 'pca.kmeans.k5',splitBy = "time_point")
p13<-ssc.plot.tsne(sce_signal,columns = 'pca.kmeans.k5',splitBy = "tippint_point")
CairoPDF("c405_cell_signal_expression_time_series.pdf",height = 6,width = 20)
ggarrange(p11,p12,p13,ncol=3,nrow=1)
dev.off()

rm(list = ls())
```





```R
options(stringsAsFactors = T)
library(sscClust)
library(dplyr)
library(Cairo)
library(reshape2)
library(plotly)
library(ggpubr)

file_names = c('c405','GSE75748','GSE79578')
for(i in 1:length(file_names)){    
    dir_name<-paste0(file_names[i],'_function_analysis')
    dirpath = paste0("C://Users//Administrator//Desktop//cell_analysis2//",dir_name)
    setwd(dirpath)
    expression_name<-paste0(file_names[i],'_expression_matrix.csv')
    entropy_name<-paste0(file_names[i],'_signal_gene_entropy.csv')
    
    entropymatrix<-read.csv(file = entropy_name,header = T)
    cell_metadata<-read.csv(file = "cell_metadata.csv",sep = ",",header =T,row.names = 1)
    
    gene_meta<-read.csv(file = "mart_export.csv",header = T)
    entropymatrix<-entropymatrix %>% distinct(Gene_name,.keep_all = T)
    signaluniquesymbol<-as.data.frame(entropymatrix$Gene_name)
    names(signaluniquesymbol)<-'Gene_name'
    signalgene_metadata<-merge(signaluniquesymbol,gene_meta,by = 'Gene_name')
    signalgene_metadata<-signalgene_metadata %>% distinct(Gene_name,.keep_all = T)
    signalgene_metadata<-signalgene_metadata %>% arrange(desc(Gene_name))
    signaluniquesymbol<-as.data.frame(signalgene_metadata$Gene_name)
    names(signaluniquesymbol)<-'Gene_name'
    row.names(signalgene_metadata)<-signalgene_metadata$Gene_name
    signalgene_metadata<-signalgene_metadata[,-1]
    entropymatrix<-merge(signaluniquesymbol,entropymatrix,by= 'Gene_name')
    entropymatrix<-entropymatrix %>% arrange(desc(Gene_name))
    row.names(entropymatrix)<-entropymatrix$Gene_name
    
    geneexpr<-read.csv(expression_name,sep = ",")
    geneexpr<-geneexpr %>% distinct(Gene_name,.keep_all = T)
    row.names(geneexpr)<-geneexpr$Gene_name
    signalgene_expr<-merge(signaluniquesymbol,geneexpr,by = 'Gene_name')
    row.names(signalgene_expr)<-signalgene_expr$Gene_name
    
    sce_signal <- ssc.build(signalgene_expr[,-1], display.name = signalgene_expr$Gene_name)
    colData(sce_signal) <- DataFrame(cell_metadata)
    
    sce_signal <- ssc.run(sce_signal,method.reduction = "pca",
                          method.clust = "SNN", SNN.k=10,SNN.method="eigen",seed=9997)
    p1<-ssc.plot.tsne(sce_signal, 
                      gene = NULL,
                      columns = ("pca.SNN.kauto"),
                      reduced.name = "pca.tsne",
                      p.ncol = 2)
    print(p1)
    p2<-ssc.plot.tsne(sce_signal, 
                      gene = NULL,
                      columns = c("source"),
                      reduced.name = "pca.tsne",
                      p.ncol = 2)
    print(p2)
    p3<-ssc.plot.tsne(sce_signal, 
                      gene = NULL,
                      columns = c("time"),
                      reduced.name = "pca.tsne",
                      p.ncol = 2)
    print(p3)
    CairoPDF("signal_gene_expression_cluster.pdf",width = 18,height = 8)
    gg1<-ggarrange(p1,p2,p3,ncol = 3)
    print(gg1)
    dev.off()
    
    sce_entropy <- ssc.build(entropymatrix[,-1], display.name = entropymatrix$Gene_name)
    colData(sce_entropy) <- DataFrame(cell_metadata)
    
    sce_entropy <- ssc.run(sce_entropy,method.reduction = "pca",
                           method.clust = "SNN", SNN.k=10,SNN.method="eigen",seed=9997)
    
    p4<-ssc.plot.tsne(sce_entropy, 
                      gene = NULL,
                      columns = ("pca.SNN.kauto"),
                      reduced.name = "pca.tsne",
                      p.ncol = 2)
    
    p5<-ssc.plot.tsne(sce_entropy, 
                      gene = NULL,
                      columns = c("source"),
                      reduced.name = "pca.tsne",
                      p.ncol = 2)

    p6<-ssc.plot.tsne(sce_entropy, 
                      gene = NULL,
                      columns = c("time"),
                      reduced.name = "pca.tsne",
                      p.ncol = 2)
    
    CairoPDF("signal_gene_entropy_cluster.pdf",width = 18,height = 8)
    gg2<-ggarrange(p4,p5,p6,ncol = 3)
    print(gg2)
    dev.off()
    
    entropy_position<-as.data.frame(sce_entropy@int_colData$reducedDims@listData$pca.tsne)
    names(entropy_position)<-c("Dim1","Dim2")
    cell_names<-row.names(cell_metadata)
    row.names(entropy_position)<-cell_names
    entropy_position[,"Entropy"]<-cell_metadata$cell_entropy
    entropy_position[,"Time"]<-cell_metadata$time
    write.csv(entropy_position,"c405_entropy_position.csv")
    entropy_position$Time <- as.factor(entropy_position$Time)
    plot_ly(entropy_position, x = ~Dim1, y = ~Dim2, z = ~Entropy, color = ~Time,size = 45) %>%
        add_markers() %>%
        layout(scene = list(xaxis = list(title = 'Dim 1'),
                            yaxis = list(title = 'Dim 2'),
                            zaxis = list(title = 'Entropy')))
}    
    rm(list = ls())

    
    
    file_names = c('c405','GSE75748','GSE79578')
    dir_name<-paste0(file_names[2],'_function_analysis')
    dirpath = paste0("C://Users//Administrator//Desktop//cell_analysis2//",dir_name)
    setwd(dirpath)
    expression_name<-paste0(file_names[2],'_expression_matrix.csv')
    entropy_name<-paste0(file_names[2],'_signal_gene_entropy.csv')
    
    entropymatrix<-read.csv(file = entropy_name,header = T)
    cell_metadata<-read.csv(file = "cell_metadata.csv",sep = ",",header =T,row.names = 1)
    
    gene_meta<-read.csv(file = "mart_export.csv",header = T)
    entropymatrix<-entropymatrix %>% distinct(Gene_name,.keep_all = T)
    signaluniquesymbol<-as.data.frame(entropymatrix$Gene_name)
    names(signaluniquesymbol)<-'Gene_name'
    signalgene_metadata<-merge(signaluniquesymbol,gene_meta,by = 'Gene_name')
    signalgene_metadata<-signalgene_metadata %>% distinct(Gene_name,.keep_all = T)
    signalgene_metadata<-signalgene_metadata %>% arrange(desc(Gene_name))
    signaluniquesymbol<-as.data.frame(signalgene_metadata$Gene_name)
    names(signaluniquesymbol)<-'Gene_name'
    row.names(signalgene_metadata)<-signalgene_metadata$Gene_name
    signalgene_metadata<-signalgene_metadata[,-1]
    entropymatrix<-merge(signaluniquesymbol,entropymatrix,by= 'Gene_name')
    entropymatrix<-entropymatrix %>% arrange(desc(Gene_name))
    row.names(entropymatrix)<-entropymatrix$Gene_name
    
    geneexpr<-read.csv(expression_name,sep = ",")
    geneexpr<-geneexpr %>% distinct(Gene_name,.keep_all = T)
    row.names(geneexpr)<-geneexpr$Gene_name
    signalgene_expr<-merge(signaluniquesymbol,geneexpr,by = 'Gene_name')
    row.names(signalgene_expr)<-signalgene_expr$Gene_name
    
    sce_signal <- ssc.build(signalgene_expr[,-1], display.name = signalgene_expr$Gene_name)
    colData(sce_signal) <- DataFrame(cell_metadata)
    
    sce_signal <- ssc.run(sce_signal,method.reduction = "pca",
                          method.clust = "SNN", SNN.k=10,SNN.method="eigen",seed=9997)
    p1<-ssc.plot.tsne(sce_signal, 
                      gene = NULL,
                      columns = ("pca.SNN.kauto"),
                      reduced.name = "pca.tsne",
                      p.ncol = 2)
    
    p2<-ssc.plot.tsne(sce_signal, 
                      gene = NULL,
                      columns = c("source"),
                      reduced.name = "pca.tsne",
                      p.ncol = 2)
    
    p3<-ssc.plot.tsne(sce_signal, 
                      gene = NULL,
                      columns = c("time"),
                      reduced.name = "pca.tsne",
                      p.ncol = 2)
    
    CairoPDF("signal_gene_expression_cluster.pdf",width = 18,height = 8)
    ggarrange(p1,p2,p3,ncol = 3)
    dev.off()
    
    sce_entropy <- ssc.build(entropymatrix[,-1], display.name = entropymatrix$Gene_name)
    colData(sce_entropy) <- DataFrame(cell_metadata)
    
    sce_entropy <- ssc.run(sce_entropy,method.reduction = "pca",
                           method.clust = "SNN", SNN.k=10,SNN.method="eigen",seed=9997)
    
    p4<-ssc.plot.tsne(sce_entropy, 
                      gene = NULL,
                      columns = ("pca.SNN.kauto"),
                      reduced.name = "pca.tsne",
                      p.ncol = 2)
    
    p5<-ssc.plot.tsne(sce_entropy, 
                      gene = NULL,
                      columns = c("source"),
                      reduced.name = "pca.tsne",
                      p.ncol = 2)
    
    p6<-ssc.plot.tsne(sce_entropy, 
                      gene = NULL,
                      columns = c("time"),
                      reduced.name = "pca.tsne",
                      p.ncol = 2)
    
    CairoPDF("signal_gene_entropy_cluster.pdf",width = 18,height = 8)
    ggarrange(p4,p5,p6,ncol = 3)
    dev.off()
    
    entropy_position<-as.data.frame(sce_entropy@int_colData$reducedDims@listData$pca.tsne)
    names(entropy_position)<-c("Dim1","Dim2")
    cell_names<-row.names(cell_metadata)
    row.names(entropy_position)<-cell_names
    entropy_position[,"Entropy"]<-cell_metadata$cell_entropy
    entropy_position[,"Time"]<-cell_metadata$time
    write.csv(entropy_position,"entropy_position.csv")
    entropy_position$Time <- as.factor(entropy_position$Time)
    plot_ly(entropy_position, x = ~Dim1, y = ~Dim2, z = ~Entropy, color = ~Time,size = 45) %>%
        add_markers() %>%
        layout(scene = list(xaxis = list(title = 'Dim 1'),
                            yaxis = list(title = 'Dim 2'),
                            zaxis = list(title = 'Entropy')))
    
    scatterplot3d(entropy_position$Dim1,entropy_position$Dim2,entropy_position$Entropy,
                  xlab = 'Dim 1',ylab = 'Dim 2',zlab = 'Entropy',
                  col = )
    
    rm(list = ls())
    
    
    
    
    file_names = c('c405','GSE75748','GSE79578')
    dir_name<-paste0(file_names[3],'_function_analysis')
    dirpath = paste0("C://Users//Administrator//Desktop//cell_analysis2//",dir_name)
    setwd(dirpath)
    expression_name<-paste0(file_names[3],'_expression_matrix.csv')
    entropy_name<-paste0(file_names[3],'_signal_gene_entropy.csv')
    
    entropymatrix<-read.csv(file = entropy_name,header = T)
    cell_metadata<-read.csv(file = "cell_metadata.csv",sep = ",",header =T,row.names = 1)
    
    gene_meta<-read.csv(file = "mart_export.csv",header = T)
    entropymatrix<-entropymatrix %>% distinct(Gene_name,.keep_all = T)
    signaluniquesymbol<-as.data.frame(entropymatrix$Gene_name)
    names(signaluniquesymbol)<-'Gene_name'
    signalgene_metadata<-merge(signaluniquesymbol,gene_meta,by = 'Gene_name')
    signalgene_metadata<-signalgene_metadata %>% distinct(Gene_name,.keep_all = T)
    signalgene_metadata<-signalgene_metadata %>% arrange(desc(Gene_name))
    signaluniquesymbol<-as.data.frame(signalgene_metadata$Gene_name)
    names(signaluniquesymbol)<-'Gene_name'
    row.names(signalgene_metadata)<-signalgene_metadata$Gene_name
    signalgene_metadata<-signalgene_metadata[,-1]
    entropymatrix<-merge(signaluniquesymbol,entropymatrix,by= 'Gene_name')
    entropymatrix<-entropymatrix %>% arrange(desc(Gene_name))
    row.names(entropymatrix)<-entropymatrix$Gene_name
    
    geneexpr<-read.csv(expression_name,sep = ",")
    geneexpr<-geneexpr %>% distinct(Gene_name,.keep_all = T)
    row.names(geneexpr)<-geneexpr$Gene_name
    signalgene_expr<-merge(signaluniquesymbol,geneexpr,by = 'Gene_name')
    row.names(signalgene_expr)<-signalgene_expr$Gene_name
    
    sce_signal <- ssc.build(signalgene_expr[,-1], display.name = signalgene_expr$Gene_name)
    colData(sce_signal) <- DataFrame(cell_metadata)
    
    sce_signal <- ssc.run(sce_signal,method.reduction = "pca",
                          method.clust = "SNN", SNN.k=10,SNN.method="eigen",seed=9997)
    p1<-ssc.plot.tsne(sce_signal, 
                      gene = NULL,
                      columns = ("pca.SNN.kauto"),
                      reduced.name = "pca.tsne",
                      p.ncol = 2)
    
    p2<-ssc.plot.tsne(sce_signal, 
                      gene = NULL,
                      columns = c("source"),
                      reduced.name = "pca.tsne",
                      p.ncol = 2)
    
    p3<-ssc.plot.tsne(sce_signal, 
                      gene = NULL,
                      columns = c("time"),
                      reduced.name = "pca.tsne",
                      p.ncol = 2)
    
    CairoPDF("signal_gene_expression_cluster.pdf",width = 18,height = 8)
    ggarrange(p1,p2,p3,ncol = 3)
    dev.off()
    
    sce_entropy <- ssc.build(entropymatrix[,-1], display.name = entropymatrix$Gene_name)
    colData(sce_entropy) <- DataFrame(cell_metadata)
    
    sce_entropy <- ssc.run(sce_entropy,method.reduction = "pca",
                           method.clust = "SNN", SNN.k=10,SNN.method="eigen",seed=9997)
    
    p4<-ssc.plot.tsne(sce_entropy, 
                      gene = NULL,
                      columns = ("pca.SNN.kauto"),
                      reduced.name = "pca.tsne",
                      p.ncol = 2)
    
    p5<-ssc.plot.tsne(sce_entropy, 
                      gene = NULL,
                      columns = c("source"),
                      reduced.name = "pca.tsne",
                      p.ncol = 2)
    
    p6<-ssc.plot.tsne(sce_entropy, 
                      gene = NULL,
                      columns = c("time"),
                      reduced.name = "pca.tsne",
                      p.ncol = 2)
    
    CairoPDF("signal_gene_entropy_cluster.pdf",width = 18,height = 8)
    ggarrange(p4,p5,p6,ncol = 3)
    dev.off()
    
    entropy_position<-as.data.frame(sce_entropy@int_colData$reducedDims@listData$pca.tsne)
    names(entropy_position)<-c("Dim1","Dim2")
    cell_names<-row.names(cell_metadata)
    row.names(entropy_position)<-cell_names
    entropy_position[,"Entropy"]<-cell_metadata$cell_entropy
    entropy_position[,"Time"]<-cell_metadata$time
    write.csv(entropy_position,"entropy_position.csv")
    entropy_position$Time <- as.factor(entropy_position$Time)
    plot_ly(entropy_position, x = ~Dim1, y = ~Dim2, z = ~Entropy, color = ~Time,size = 45) %>%
        add_markers() %>%
        layout(scene = list(xaxis = list(title = 'Dim 1'),
                            yaxis = list(title = 'Dim 2'),
                            zaxis = list(title = 'Entropy')))
    rm(list = ls())
    

    setwd("C://Users//Administrator//Desktop//cell_analysis2//GSE102066_function_analysis/")
    entropymatrix<-read.csv(file = "GSE102066_signal_gene_entropy.csv",header = T)
    cell_metadata<-read.csv(file = "cell_metadata.csv",sep = ",",header =T,row.names = 1)
    
    gene_meta<-read.csv(file = "mart_export.csv",header = T)
    entropymatrix<-entropymatrix %>% distinct(Gene_name,.keep_all = T)
    signaluniquesymbol<-as.data.frame(entropymatrix$Gene_name)
    names(signaluniquesymbol)<-'Gene_name'
    signalgene_metadata<-merge(signaluniquesymbol,gene_meta,by = 'Gene_name')
    signalgene_metadata<-signalgene_metadata %>% distinct(Gene_name,.keep_all = T)
    signalgene_metadata<-signalgene_metadata %>% arrange(desc(Gene_name))
    signaluniquesymbol<-as.data.frame(signalgene_metadata$Gene_name)
    names(signaluniquesymbol)<-'Gene_name'
    row.names(signalgene_metadata)<-signalgene_metadata$Gene_name
    signalgene_metadata<-signalgene_metadata[,-1]
    entropymatrix<-merge(signaluniquesymbol,entropymatrix,by= 'Gene_name')
    entropymatrix<-entropymatrix %>% arrange(desc(Gene_name))
    row.names(entropymatrix)<-entropymatrix$Gene_name
    
    
    genesymbol<-as.data.frame(gene_meta$Gene_stable_ID)
    names(genesymbol)<-"Gene_stable_ID"
    geneexpr<-read.csv("GSE102066_expresion_matrix.csv",sep = ",")
    geneexpr<-geneexpr %>% distinct(Gene_stable_ID,.keep_all = T)
    row.names(geneexpr)<-geneexpr$Gene_stable_ID
    
    signalgene_expr<-merge(genesymbol,geneexpr,by = 'Gene_stable_ID')
    signalgene_expr<-merge(signalgene_metadata,geneexpr,by = 'Gene_stable_ID')
    signalgene_expr<-signalgene_expr[,c(-1,-3,-4)]
    row.names(signalgene_expr)<-signalgene_expr$gene_short_name
    signalgene_expr<-signalgene_expr[,c(-1)]
    signalgene_expr<-as.data.frame(t(signalgene_expr))
    
    sce_signal <- ssc.build(signalgene_expr[,-1], display.name = signalgene_expr$Gene_stable_ID)
    colData(sce_signal) <- DataFrame(cell_metadata)
    
    sce_signal <- ssc.run(sce_signal,method.reduction = "pca",
                          method.clust = "SNN", SNN.k=10,SNN.method="eigen",seed=9997)
    p1<-ssc.plot.tsne(sce_signal, 
                      gene = NULL,
                      columns = ("pca.SNN.kauto"),
                      reduced.name = "pca.tsne",
                      p.ncol = 2)
    
    p2<-ssc.plot.tsne(sce_signal, 
                      gene = NULL,
                      columns = c("source"),
                      reduced.name = "pca.tsne",
                      p.ncol = 2)
    
    p3<-ssc.plot.tsne(sce_signal, 
                      gene = NULL,
                      columns = c("time"),
                      reduced.name = "pca.tsne",
                      p.ncol = 2)
    
    CairoPDF("signal_gene_expression_cluster.pdf",width = 18,height = 8)
    ggarrange(p1,p2,p3,ncol = 3)
    dev.off()
    
    sce_entropy <- ssc.build(entropymatrix[,-1], display.name = entropymatrix$Gene_name)
    colData(sce_entropy) <- DataFrame(cell_metadata)
    
    sce_entropy <- ssc.run(sce_entropy,method.reduction = "pca",
                           method.clust = "SNN", SNN.k=10,SNN.method="eigen",seed=9997)
    
    p4<-ssc.plot.tsne(sce_entropy, 
                      gene = NULL,
                      columns = ("pca.SNN.kauto"),
                      reduced.name = "pca.tsne",
                      p.ncol = 2)
    
    p5<-ssc.plot.tsne(sce_entropy, 
                      gene = NULL,
                      columns = c("source"),
                      reduced.name = "pca.tsne",
                      p.ncol = 2)
    
    p6<-ssc.plot.tsne(sce_entropy, 
                      gene = NULL,
                      columns = c("time"),
                      reduced.name = "pca.tsne",
                      p.ncol = 2)
    
    CairoPDF("signal_gene_entropy_cluster.pdf",width = 18,height = 8)
    ggarrange(p4,p5,p6,ncol = 3)
    dev.off()
    
    entropy_position<-as.data.frame(sce_entropy@int_colData$reducedDims@listData$pca.tsne)
    names(entropy_position)<-c("Dim1","Dim2")
    cell_names<-row.names(cell_metadata)
    row.names(entropy_position)<-cell_names
    entropy_position[,"Entropy"]<-cell_metadata$cell_entropy
    entropy_position[,"Time"]<-cell_metadata$time
    write.csv(entropy_position,"entropy_position.csv")
    entropy_position$Time <- as.factor(entropy_position$Time)
    plot_ly(entropy_position, x = ~Dim1, y = ~Dim2, z = ~Entropy, color = ~Time,size = 45) %>%
        add_markers() %>%
        layout(scene = list(xaxis = list(title = 'Dim 1'),
                            yaxis = list(title = 'Dim 2'),
                            zaxis = list(title = 'Entropy')))
    
    
    expr_time<-c(rep('day 0',80),rep('day 1',78),rep('day 5',85),
                 rep('day 7',80),rep('day 10',79),rep('day 30',81))
    expr_time<-factor(expr_time,
                      levels = c("day 0","day 1","day 5","day 7","day 10","day 30"))
    expr_time

    signalgene_expr[,'Time']<-expr_time
    signalgene_expr$Time
    
    pdf("every_gene_expression_by_time.pdf")
    for (i in 1:180) {
        p<-ggplot(data = signalgene_expr,aes_string(x = 'Time',
                                                    y = names(signalgene_expr)[i],group = 1))
        p<-p+labs(title = names(signalgene_expr)[i])+
            theme(title = element_text(size = 16,face = 'bold'))+
            theme(plot.title = element_text(hjust = 0.5))
        p<-p+ylab("Expression")
        p<-p + geom_jitter(aes(color = Time), size = 3)+geom_smooth()
        p<-p+theme(legend.position="none",
                   panel.background = element_blank(),
                   panel.border = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   axis.line.y = element_line(colour = "black"))
        print(p)
        }
    dev.off()
    
    entropymatrix<-entropymatrix[,-1]
    entropymatrix<-as.data.frame(t(entropymatrix))
    entropymatrix[,'Time']<-expr_time
    entropymatrix$Time
    pdf("every_gene_entropy_by_time.pdf")
    for (i in 1:182) {
        p1<-ggplot(data = entropymatrix,aes_string(x = 'Time',
                                                    y = names(entropymatrix)[i],group = 1))
        p1<-p1+ylab("Entropy")
        p1<-p1+geom_jitter(aes(color = Time), size = 3)+geom_smooth()
        p1<-p1+theme(legend.position="none",
                   panel.background = element_blank(),
                   panel.border = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   axis.line.y = element_line(colour = "black"))
        p2<-ggplot2
        print(p)
    }
    dev.off()
    rm(list = ls())
    
    #p1<-p1+labs(title = names(entropymatrix)[i])+
    #theme(title = element_text(size = 16,face = 'bold'))+
    #    theme(plot.title = element_text(hjust = 0.5))
```





```R
rm(list = ls())

library(ggplot2)
library(ggpubr)
library(Cairo)
library(dplyr)
##########################################################################################
setwd("C://Users//Administrator//Desktop//cell_analysis2//GSE79578_function_analysis/")
entropymatrix<-read.csv(file = "GSE79578_signal_gene_entropy.csv",header = T)
cell_metadata<-read.csv(file = "cell_metadata.csv",sep = ",",header =T,row.names = 1)

gene_meta<-read.csv(file = "mart_export.csv",header = T)
entropymatrix<-entropymatrix %>% distinct(Gene_name,.keep_all = T)
entropymatrix<-entropymatrix %>% arrange(desc(Gene_name))
signaluniquesymbol<-as.data.frame(entropymatrix$Gene_name)
names(signaluniquesymbol)<-'Gene_name'
signalgene_metadata<-merge(signaluniquesymbol,gene_meta,by = 'Gene_name')
signalgene_metadata<-signalgene_metadata %>% distinct(Gene_name,.keep_all = T)
signalgene_metadata<-signalgene_metadata %>% arrange(desc(Gene_name))
signaluniquesymbol<-as.data.frame(signalgene_metadata$Gene_name)
names(signaluniquesymbol)<-'Gene_name'
row.names(signalgene_metadata)<-signalgene_metadata$Gene_name
signalgene_metadata<-signalgene_metadata[,-1]
entropymatrix<-merge(signaluniquesymbol,entropymatrix,by= 'Gene_name')
entropymatrix<-entropymatrix %>% arrange(desc(Gene_name))
row.names(entropymatrix)<-entropymatrix$Gene_name

geneexpr<-read.csv("GSE79578_expression_matrix.csv",sep = ",")
geneexpr<-geneexpr %>% distinct(Gene_name,.keep_all = T)
geneexpr<-geneexpr %>% arrange(desc(Gene_name))
row.names(geneexpr)<-geneexpr$Gene_name

geneexpr<-merge(signaluniquesymbol,geneexpr,by = 'Gene_name')
geneexpr<-geneexpr %>% arrange(desc(Gene_name))
row.names(geneexpr)<-geneexpr$Gene_name

geneexpr<-geneexpr[,-1]
geneexpr<-as.data.frame(t(geneexpr))

entropymatrix<-entropymatrix[,-1]
entropymatrix<-as.data.frame(t(entropymatrix))

expr_time<-c(rep('2h',82),rep('12h',86),rep('24h',89),
             rep('48h',81))
expr_time<-factor(expr_time,
                  levels = c("2h","12h","24h","48h"))
expr_time

geneexpr[,'Time']<-expr_time
geneexpr$Time

entropymatrix[,'Time']<-expr_time
entropymatrix$Time

rownames(geneexpr)<-NULL
rownames(entropymatrix)<-NULL

pdf("every_gene_expression_by_time.pdf")
for (i in 1:184) {
  p1<-ggplot(data = geneexpr,aes_string(x = 'Time',y = as.character(names(geneexpr)[i]),group = 1))
  
  p1<-p1+ylab("Expression")
  
  p1<-p1+
    theme(title = element_text(size = 16,face = 'bold'))+
    theme(plot.title = element_text(hjust = 0.5))
  
  p1<-p1 + geom_jitter(aes(color = Time), size = 3)+
    stat_summary(fun.y = "mean", size = 1.5, geom = "line")
  
  p1<-p1+theme(legend.position="none",
             panel.background = element_blank(),
             panel.border = element_blank(),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             axis.line.y = element_line(colour = "black"))

  p2<-ggplot(data = entropymatrix,aes_string(x ='Time', y=names(entropymatrix)[i],group = 1))
  
  p2<-p2+labs(title = names(entropymatrix)[i])+
    theme(title = element_text(size = 16,face = 'bold'))+
    theme(plot.title = element_text(hjust = 0.5))
  
  p2<-p2+ylab("Entropy")
  
  p2<-p2 + geom_jitter(aes(color = Time), size = 3)+
    stat_summary(fun.y = "mean", size = 1.5, geom = "line")
  
  p2<-p2+theme(legend.position="none",
             panel.background = element_blank(),
             panel.border = element_blank(),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             axis.line.y = element_line(colour = "black"))
  
  p<-ggarrange(p2,p1,heights=c(1/2, 1/2),ncol = 1, nrow = 2,common.legend = FALSE,align = "v")
  
  print(p)
}

dev.off()

entropymatrix<-entropymatrix[,-1]
entropymatrix<-as.data.frame(t(entropymatrix))
entropymatrix[,'Time']<-expr_time
entropymatrix$Time
pdf("every_gene_entropy_by_time.pdf")
for (i in 1:184) {
  p<-ggplot(data = entropymatrix,aes_string(x = 'Time',
                                            y = names(entropymatrix)[i],group = 1))
  p<-p+labs(title = names(entropymatrix)[i])+
    theme(title = element_text(size = 16,face = 'bold'))+
    theme(plot.title = element_text(hjust = 0.5))
  p<-p+ylab("Entropy")
  p<-p + geom_jitter(aes(color = Time), size = 3)+
    stat_summary(fun.y = "mean", size = 1.5, geom = "line")
  
  p<-p+theme(legend.position="none",
             panel.background = element_blank(),
             panel.border = element_blank(),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             axis.line.y = element_line(colour = "black"))
  print(p)
}
dev.off()
rm(list = ls())

##########################################################################################
setwd("C://Users//Administrator//Desktop//cell_analysis2//GSE102066_function_analysis//")
entropymatrix<-read.csv(file = "GSE102066_signal_gene_entropy.csv",header = T)
cell_metadata<-read.csv(file = "cell_metadata.csv",sep = ",",header =T,row.names = 1)

gene_meta<-read.csv(file = "mart_export.csv",header = T)
entropymatrix<-entropymatrix %>% distinct(Gene_name,.keep_all = T)
signaluniquesymbol<-as.data.frame(entropymatrix$Gene_name)
names(signaluniquesymbol)<-'Gene_name'
signalgene_metadata<-merge(signaluniquesymbol,gene_meta,by = 'Gene_name')
signalgene_metadata<-signalgene_metadata %>% distinct(Gene_name,.keep_all = T)
signalgene_metadata<-signalgene_metadata %>% arrange(desc(Gene_name))
signaluniquesymbol<-as.data.frame(signalgene_metadata$Gene_name)
names(signaluniquesymbol)<-'Gene_name'
row.names(signalgene_metadata)<-signalgene_metadata$Gene_name
signalgene_metadata<-signalgene_metadata[,-1]
entropymatrix<-merge(signaluniquesymbol,entropymatrix,by= 'Gene_name')
entropymatrix<-entropymatrix %>% arrange(desc(Gene_name))
row.names(entropymatrix)<-entropymatrix$Gene_name

geneexpr<-read.csv("GSE102066_expresion_matrix.csv",sep = ",")
geneexpr<-geneexpr %>% distinct(Gene_name,.keep_all = T)
row.names(geneexpr)<-geneexpr$Gene_name

geneexpr<-merge(signaluniquesymbol,geneexpr,by = 'Gene_name')
row.names(geneexpr)<-geneexpr$Gene_name

geneexpr<-geneexpr[,-1]


geneexpr<-as.data.frame(t(geneexpr))

expr_time<-c(rep('2h',82),rep('12h',86),rep('24h',89),
             rep('48h',81))
expr_time<-factor(expr_time,
                  levels = c("2h","12h","24h","48h"))
expr_time

geneexpr[,'Time']<-expr_time
geneexpr$Time

pdf("every_gene_expression_by_time.pdf")
for (i in 1:184) {
  p<-ggplot(data = geneexpr,aes_string(x = 'Time',y = as.character(names(geneexpr)[i]),group = 1))
  
  p<-p+labs(title = names(geneexpr)[i])+theme(title = element_text(size = 16,face = 'bold'))+theme(plot.title = element_text(hjust = 0.5))
  
  p<-p+ylab("Expression")
  
  p<-p + geom_jitter(aes(color = Time), size = 3)+geom_smooth()
  
  p<-p+theme(legend.position="none",
             panel.background = element_blank(),
             panel.border = element_blank(),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             axis.line.y = element_line(colour = "black"))
  
  print(p)
}
dev.off()

entropymatrix<-entropymatrix[,-1]
entropymatrix<-as.data.frame(t(entropymatrix))
entropymatrix[,'Time']<-expr_time
entropymatrix$Time
pdf("every_gene_entropy_by_time.pdf")
for (i in 1:184) {
  p<-ggplot(data = entropymatrix,aes_string(x = 'Time',
                                            y = names(entropymatrix)[i],group = 1))
  p<-p+labs(title = names(entropymatrix)[i])+
    theme(title = element_text(size = 16,face = 'bold'))+
    theme(plot.title = element_text(hjust = 0.5))
  p<-p+ylab("Entropy")
  p<-p + geom_jitter(aes(color = Time), size = 3)+geom_smooth()
  p<-p+theme(legend.position="none",
             panel.background = element_blank(),
             panel.border = element_blank(),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             axis.line.y = element_line(colour = "black"))
  print(p)
}
dev.off()
rm(list = ls())
```





```R
  p2<-ggplot(data = geneexpr_average,aes(x=Time,y=Aaas,group=1))+
    geom_line(linetype = "solid",color="#6495ED",size=1.5)+
    labs(title = names(geneexpr_average)[i],x=NULL,y=NULL)+
    theme(title = element_text(size = 16,face = 'bold'))+
    theme(plot.title = element_text(hjust = 0.5))+
    geom_point()+
    theme(legend.position="none",
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line.y = element_line(colour = "black"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  p<-ggarrange(p2,p1,heights=c(1/5, 4/5),ncol = 1, nrow = 2,common.legend = TRUE,legend="right",align = "v")
```





```R
setwd("C://Users//Administrator//Desktop//paper//osdata")
options(stringsAsFactors = F)
library(survival)
library(parallel)
library(survminer)
library(Cairo)

#Dead:0
#Alive:1

ECSA_entropy<-read.csv("ESCA_signals_gene_entropy_matrix.csv",header = T)
HNSC_entropy<-read.csv("HNSC_signals_gene_entropy_matrix.csv",header = T)
READ_entropy<-read.csv("READ_signals_gene_entropy_matrix.csv",header = T)
UCEC_entropy<-read.csv("UCEC_signals_gene_entropy_matrix.csv",header = T)


ECSA_surv<-Surv(ECSA_entropy$time,ECSA_entropy$vital_status)

log_rank_p<- apply(ECSA_entropy[7:length(names(ECSA_entropy))],2,function(values){
  group<-ifelse(values>=mean(na.omit(values)),"Chaotic","Static")
  kmfit2<- survival::survfit(ECSA_surv~group,data = ECSA_entropy)
  data.survdiff<-survival::survdiff(ECSA_surv~group)
  p.value = 1 - pchisq(data.survdiff$chisq,length(data.survdiff$n)-1)
})

log_rank_p<-log_rank_p[log_rank_p<0.05]

gene_marker<-as.data.frame(sort(log_rank_p))
gene_marker_name<-row.names(gene_marker)
gene_marker_name
write.table(gene_marker,file = 'ECSA_entropy_gene_marker.txt',sep = '\t')

# "PER2"   "DAPK3"  "HCST"   "RAET1E" "ULBP1"  "ULBP3" 
group<-ifelse(ECSA_entropy$HCST>mean(ECSA_entropy$HCST),
              "High","Low")
sfit<-survfit(Surv(time,vital_status)~group,data = ECSA_entropy)
CairoPDF(file = 'ECSA_entropy_PER2.pdf',width = 8,height = 10)
ggsurvplot(sfit,
           legend.title = 'Gene Entropy',
           conf.int = T,
           pval = T,
           title = 'HCST',
           risk.table = T,
           palette = c("#E7B800","#2E9FDF"),
           risk.table.col = 'strata',
           risk.table.height = 0.25,
           xlab = "Time In Days",
           pval.method = TRUE,
           surv.median.line = "hv",
           ncensor.plot = T)
dev.off()


HNSC_surv<-Surv(HNSC$time,HNSC$vital_status)

log_rank_p<- apply(HNSC[7:length(names(HNSC))],2,function(values){
  group<-ifelse(values>=mean(na.omit(values)),"Chaotic","Static")
  kmfit2<- survival::survfit(HNSC_surv~group,data = HNSC)
  data.survdiff<-survival::survdiff(HNSC_surv~group)
  p.value = 1 - pchisq(data.survdiff$chisq,length(data.survdiff$n)-1)
})

log_rank_p<-log_rank_p[log_rank_p<0.05]

gene_marker<-as.data.frame(sort(log_rank_p))
gene_marker_name<-row.names(gene_marker)
gene_marker_name
write.table(gene_marker,file = 'HNSC_gene_marker.txt',sep = '\t')


#"FLNB" "CKS2" "CCNB3" "PLK2" "PKMYT1" "PALB2" "CKS1B"
#"IL20RB" "TNFSF4" "HIST1H2AJ" "RFC2" "ASAP1" 
#"RAD9B" "HMGA2" "KIF23" "XRCC3" "ENO2" "GTSE1"
#"TTI1"

group<-ifelse((HNSC$TTI1)>mean(HNSC$TTI1),"Chaotic","Static")
sfit<-survfit(Surv(time,vital_status)~group,data = HNSC)
CairoPDF(file = 'HNSC_TTI1.pdf',width = 8,height = 10)
ggsurvplot(sfit,conf.int = F,pval = T,title = 'TTI1',risk.table = T,
           palette = c("#E7B800","#2E9FDF"),risk.table.col = 'strata',
           risk.table.height = 0.25,xlab = "Time In Days",pval.method = TRUE,
           surv.median.line = "hv",ncensor.plot.height = TRUE,ncensor.plot = F)
dev.off()


READ_surv<-Surv(READ$time,READ$vital_status)

log_rank_p<- apply(READ[7:length(names(READ))],2,function(values){
  group<-ifelse(values>=mean(na.omit(values)),"Chaotic","Static")
  kmfit2<- survival::survfit(READ_surv~group,data = READ)
  data.survdiff<-survival::survdiff(READ_surv~group)
  p.value = 1 - pchisq(data.survdiff$chisq,length(data.survdiff$n)-1)
})

log_rank_p<-log_rank_p[log_rank_p<0.05]

gene_marker<-as.data.frame(sort(log_rank_p))
gene_marker_name<-row.names(gene_marker)
gene_marker_name
write.table(gene_marker,file = 'READ_gene_marker.txt',sep = '\t')


#"PALB2" "RAD54L" "POLN" "CENPX" "SLC7A11" "ULBP1" "HCST"
#"IDH3G" "RANBP1" "MMP13"
group<-ifelse((READ$MMP13)>mean(READ$MMP13),"Chaotic","Static")
sfit<-survfit(Surv(time,vital_status)~group,data = READ)
CairoPDF(file = 'READ_MMP13.pdf',width = 8,height = 10)
ggsurvplot(sfit,conf.int = F,pval = T,title = 'MMP13',risk.table = T,
           palette = c("#E7B800","#2E9FDF"),risk.table.col = 'strata',
           risk.table.height = 0.25,xlab = "Time In Days",pval.method = TRUE,
           surv.median.line = "hv",ncensor.plot.height = TRUE,ncensor.plot = F)
dev.off()


UCEC_surv<-Surv(UCEC$time,UCEC$vital_status)

log_rank_p<- apply(UCEC[7:length(names(UCEC))],2,function(values){
  group<-ifelse(values>=mean(na.omit(values)),"Chaotic","Static")
  kmfit2<- survival::survfit(UCEC_surv~group,data = UCEC)
  data.survdiff<-survival::survdiff(UCEC_surv~group)
  p.value = 1 - pchisq(data.survdiff$chisq,length(data.survdiff$n)-1)
})

log_rank_p<-log_rank_p[log_rank_p<0.05]

gene_marker<-as.data.frame(sort(log_rank_p))
gene_marker_name<-row.names(gene_marker)
gene_marker_name
write.table(gene_marker,file = 'UCEC_gene_marker.txt',sep = '\t')


#"ENO4" "LTB" "PARVG" "WTIP" "TCL1B" "ELK4" "GDF15" "MUTYH"
#"DSPP" "BIRC7" "TNFSF18" "PPM1B" "IQGAP3" "CTNNBIP1" "JARID2"
#"FSCN1" "S100A7" "RFC2" "GLYCTK" "ERG"

group<-ifelse((UCEC$ERG)>mean(UCEC$ERG),"Chaotic","Static")
sfit<-survfit(Surv(time,vital_status)~group,data = UCEC)
CairoPDF(file = 'UCEC_ERG.pdf',width = 8,height = 10)
ggsurvplot(sfit,conf.int = F,pval = T,title = 'ERG',risk.table = T,
           palette = c("#E7B800","#2E9FDF"),risk.table.col = 'strata',
           risk.table.height = 0.25,xlab = "Time In Days",pval.method = TRUE,
           surv.median.line = "hv",ncensor.plot.height = TRUE,ncensor.plot = F)
dev.off()

```





```R
setwd("C://Users//Administrator//Desktop//paper//osdata")
options(stringsAsFactors = F)
library(ggbiplot)
#先安装PCA可视化包ggbiplot
#library(devtools)
#install_github("vqv/ggbiplot")

library(ggbiplot)

ECSA_death_pca<-read.table(file = "ESCA_entropy_death_PCA.txt",sep = "\t",header = T)

ECSA_death_status<-as.data.frame(ECSA_death_pca$vital_status)
table(ECSA_death_status)
ECSA_death_status<-c(rep('Alive',67),rep("Dead",40))
ECSA_death_status<-factor(ECSA_death_status,levels = c('Alive',"Dead"))

ECSA_death_pca<-ECSA_death_pca[,-1]

ECSA_death.pca <- prcomp(ECSA_death_pca, scale. = TRUE)
ggbiplot(ECSA_death.pca, obs.scale = 1, 
         var.scale = 1,
         groups = ECSA_death_status, 
         ellipse = TRUE, circle = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')


ECSA_transition_pca<-read.table(file = "ESCA_entropy_transition_PCA.txt",sep = "\t",header = T)

ECSA_transition_status<-as.data.frame(ECSA_transition_pca$AJCC)
table(ECSA_transition_status)
ECSA_transition_status<-c(rep('Before',85),rep("After",22))
ECSA_transition_status<-factor(ECSA_transition_status,levels = c('Before',"After"))

ECSA_transition_pca<-ECSA_transition_pca[,-1]
class(ECSA_transition_pca)

ECSA_transition.pca <- prcomp(ECSA_transition_pca, scale. = F)
ggbiplot(ECSA_transition.pca, obs.scale = 1, 
         var.scale = 1,
         groups = ECSA_transition_status, 
         ellipse = TRUE, circle = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')
summary(ECSA_transition.pca)



train<- ECSA_transition_pca ## Choose the train.csv file downloaded from the link above  
library(Rtsne)
## Curating the database for analysis with both t-SNE and PCA
Labels<-train$AJCC
train$AJCC<-as.factor(train$AJCC)
## for plotting
colors = rainbow(length(unique(train$AJCC)))
names(colors) = unique(train$AJCC)

## Executing the algorithm on curated data
tsne<- Rtsne(train[,-1], dims = 2, perplexity=5, verbose=TRUE, max_iter = 30)

## Plotting
plot(tsne$Y,main="tsne")
text(tsne$Y, labels=train$AJCC, col=colors[train$AJCC])
```

```R
#绘制complexheatmap
##################################################################################################
setwd("C://Users//Administrator//Desktop//paper//Figures//Figure5//ESCA")
options(stringsAsFactors = F)
##################################################################################################
library(ComplexHeatmap)
library(viridis)
library(seriation)
library(circlize)
library(ggpubr)
library(Cairo)
library(dplyr)
##################################################################################################

ESCA_entropy<-read.csv(file = "ESCA_signals_gene_entropy_matrix.csv",header = T,check.names = F)
ESCA_expression<-read.csv(file="ESCA_tumor_exprssion_matrix.csv",header = T,check.names = F)
ESCA_annotations<-read.csv(file = "ESCA_annotation.csv",header = T,check.names = F)
ESCA_important_genelist<-read.csv(file = "ESCA_SNE_genelist.csv",header = T,check.names = F)
ESCA_important_gene_expression<-merge(ESCA_important_genelist,ESCA_expression,by = 'Gene')

ESCA_entropy<-ESCA_entropy %>% arrange(desc(ID))
rownames(ESCA_entropy)<-ESCA_entropy$ID
ESCA_entropy<-ESCA_entropy[,-1]


ESCA_important_gene_expression<-ESCA_important_gene_expression %>% arrange(desc(Gene))
rownames(ESCA_important_gene_expression)<-ESCA_important_gene_expression$Gene
ESCA_important_gene_expression<-ESCA_important_gene_expression[,-1]


##################################################################################################
col_fun = colorRamp2(c(-4, -2, 0, 2,4), c("RoyalBlue4", "Orchid2", "white", "Tomato","red"))

ha1<-HeatmapAnnotation(Patient = ESCA_annotations$Patient,
                       AJCC = sample(ESCA_annotations$AJCC),
                       Status = ESCA_annotations$Status,
                       Lasting = anno_barplot(ESCA_annotations$Lasting,axis = T,border = F,baseline = 0),
                       Entropy = anno_points(ESCA_annotations$`Network-Entropy`,border = T,pch = 20 ),
                       annotation_name_side = "right",
                       gap = unit(0.2, "cm"),
                       show_legend = c(FALSE,TRUE,TRUE,TRUE,TRUE)
                       )

ha_stage_left<-HeatmapAnnotation(AJCC = ESCA_annotations$AJCC,
                                 col = list(AJCC = c("Stage I"="ForestGreen","Stage IIA"="LimeGreen","Stage IIB"="Yellow4","Stage IIIA"="IndianRed1","Stage IIIB"="Firebrick","Stage IV"="LightSalmon4")),
                                 annotation_name_side = "left")
ha_stage_right<-HeatmapAnnotation(AJCC = ESCA_annotations$AJCC,
                                  col = list(AJCC = c("Stage I"="ForestGreen","Stage IIA"="LimeGreen","Stage IIB"="Yellow4","Stage IIIA"="IndianRed1","Stage IIIB"="Firebrick","Stage IV"="LightSalmon4")),
                                 annotation_name_side = "right")
##################################################################################################
mat1<-as.matrix(ESCA_entropy)
o1 = seriate(max(mat1) - mat1, method = "BEA_TSP")
mat2<-as.matrix(ESCA_important_gene_expression)
o2 = seriate(max(mat2) - mat2, method = "BEA_TSP")

CairoPDF(file = "ESCA_entropy_vs_expression_all.pdf",width = 16,height = 12)
p1<-Heatmap(scale(mat1),
        col = col_fun,
        show_column_names=F,
        show_row_names = T,
        cluster_rows = F,
        row_names_gp = gpar(fontsize = 2.5),
        cluster_columns = F,
        row_title = "Vital Genes",
        row_names_side = 'left',
        column_title = "ESCA Important Gene Entropy",
        column_title_side = 'top',
        heatmap_legend_param = list(title = 'Entropy'),
        row_order = get_order(o1, 1), 
        column_order = get_order(o, 2),
        bottom_annotation = ha_stage_left)

p2<-Heatmap(scale(ESCA_important_gene_expression),
            col = col_fun,
            show_column_names=F,
            show_row_names = T,
            row_names_gp = gpar(fontsize = 2.5),
            cluster_columns = F,
            cluster_rows = F,
            row_title = "Vital Genes",
            row_names_side = 'left',
            column_title = "ESCA Important Gene Expression",
            heatmap_legend_param = list(title = 'Expression'),
            row_order = get_order(o2, 1),
            column_order = get_order(o2, 2),
            column_title_side = 'top',
            bottom_annotation = ha_stage_right)

ht_list<-p1+p2

draw(ht_list,merge_legends=T,legend_border=T,gap=unit(0.2,"cm"))

dev.off()
```





```R
BiocManager::install("RnBeads")
BiocManager::install("YAPSA")
library(devtools)
install_bitbucket("weischenfeldt/prescient@master")
BiocManager::install("randomSurvivalForest")
BiocManager::install("ComplexHeatmap")
BiocManager::install("genefilter")
BiocManager::install("ggbiplot")

```

