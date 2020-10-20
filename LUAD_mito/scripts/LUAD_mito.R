if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("karyoploteR")

install.packages("bezier")
library(karyoploteR)
library(GenomicFeatures)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ensembldb")
BiocManager::install("EnsDb.Hsapiens.v86")

## test
?plotKaryotype()
custom.genome <- toGRanges(data.frame(chr=c("MT"), start=c(1), end=c(20000)))
kp <- plotKaryotype(genome=custom.genome,zoom=toGRanges("MT:1-20000"))
kpAddBaseNumbers(kp, tick.dist = 1000, add.units = TRUE)
bam = "/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/mito/data/S10T_MT.bam"
kp <- kpPlotBAMCoverage(kp, data=bam,col="gold", border="red")
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.coverage)
?makeTxDbFromGFF
library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86
edb_mt <- addFilter(edb, SeqNameFilter("MT"))
mt_genes = genes(edb_mt)

library(GenomicFeatures)
gtf = "/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/mito/data/MT.gtf"
txdb <- GenomicFeatures::makeTxDbFromGFF(file=gtf, format="auto")
columns(txdb)
all_genes <- genes(txdb)

#genes = exonsBy(txdb, by="gene")
# samtools depth
depth = read_tsv("/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/mito/data/MT_depth.txt",col_names = F)
ggplot(depth,aes(x=X2,y=X3)) + geom_point()


## run
library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86
edb_mt <- addFilter(edb, SeqNameFilter("MT"))
mt_genes = genes(edb_mt)

kp1 <- plotKaryotype(genome=custom.genome,plot.type=2)
kpAddBaseNumbers(kp1, tick.dist = 1000, add.units = TRUE,cex = 0.8)
kpPlotMarkers(kp1, data=mt_genes,labels=mt_genes$symbol ,
              r1=0.5, cex=0.8, adjust.label.position = TRUE,data.panel = 2)
bam = "/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/mito/advanced/advanced_MT.bam"
kp1 <- kpPlotBAMCoverage(kp1, data=bam,col="gold", border="red")
kpAxis(kp1, ymax=kp1$latest.plot$computed.values$max.coverage)
#?kpAddCytobandLabels(kp1)
#kpPlotRegions(kp1, data=mt_genes,avoid.overlapping=TRUE,num.layers=25,layer.margin = 0.01)

df.mt = as.data.frame(mt_genes)
genes.mt = df.mt$symbol

#samtools depth
depth = read_tsv("/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/mito/advanced/bed_res.txt",col_names = F)
depth = depth[depth$X1=="MT",]
write_tsv(depth,"/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/mito/advanced/MT_bedtools_depth.txt")
min(depth$X3)
max(depth$X3)
## early
szl = readRDS("/SGRNJ01/Aftersales/P2018009_Lung/20190624_jiace/rds/new_ident.rds")
lc30 = SubsetData(szl,subset.name = "orig.ident",accept.value = "LC30_2",subset.raw = T)
mito.genes <- grep(pattern = "^MT-", x = rownames(x = lc30@data), value = TRUE)
lc30 = SetAllIdent(lc30,"new_ident")
lc30 = StashIdent(object = lc30, save.name = "cell_type")
lc30@meta.data$type = "early"

## advanced
s76 = SubsetData(mydata,subset.name = "orig.ident",accept.value = "S76T",subset.raw = T)
s76 = SetAllIdent(s76,"big_type")
s76 = SubsetData(s76,ident.remove = c("Epithelial","LUSC"),subset.raw = T)
s76 = StashIdent(object = s76, save.name = "cell_type")
s76@meta.data$type = "advanced"

#combine
combined = MergeSeurat(lc30,s76)
combined = SetAllIdent(combined,"cell_type")
sdp <- SplitDotPlotGG(cca, genes.plot = rev(markers.to.plot), cols.use = c("blue","red"),
                      x.lab.rot = T, plot.legend = T, dot.scale = 8, do.return = T, grouping.var = "type")

#cca
ctrl = lc30
stim = s76
ctrl <- FindVariableGenes(ctrl, do.plot = F)
stim <- FindVariableGenes(stim, do.plot = F)
g.1 <- head(rownames(ctrl@hvg.info), 1000)
g.2 <- head(rownames(stim@hvg.info), 1000)
genes.use <- unique(c(g.1, g.2))
genes.use <- intersect(genes.use, rownames(ctrl@scale.data))
genes.use <- intersect(genes.use, rownames(stim@scale.data))
cca =  RunCCA(lc30,s76, genes.use = genes.use, num.cc = 30)
cca <- AlignSubspace(cca, reduction.type = "cca", grouping.var = "type", 
                                 dims.align = 1:20)
cca = SetAllIdent(cca,"cell_type")




markers.to.plot <- mito.genes
sdp <- SplitDotPlotGG(cca, genes.plot = "CD4", cols.use = c("blue","red"),
 x.lab.rot = T, plot.legend = T,  do.return = T, grouping.var = "type")
FeaturePlot(cca,mito.genes)

x = cca@meta.data
x = x %>% mutate(cell_ident=paste(cell_type,"_",type,sep=""))
cca@meta.data = x
rownames(cca@meta.data)=cca@cell.names
cca = SetAllIdent(cca,"cell_ident")
temp = cca
cca@ident = factor(as.character(cca@ident),levels=rev(levels(cca@ident)))

setwd("/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/mito")
pdf("./plots/mito_gene_exp.pdf",width=14,height=10)
DotPlot(cca,genes.plot = rev(genes.mt),plot.legend = TRUE,x.lab.rot = T) 
dev.off()
#saveRDS(cca,"combined_cca.rds")


##
df.mt[df.mt$symbol=="MT-ND2",]
