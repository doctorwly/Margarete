#!/usr/bin/env Rscript
library(argparser)
library(tidyverse)


argv <- arg_parser('')
argv <- add_argument(argv,"--config", help="config file.If rds not provided,Required")
argv <- add_argument(argv,"--outdir", help="the output dir.Required")
argv <- add_argument(argv,"--is_10X", help="10X or not.if not provided,default=no",default="no")
argv <- add_argument(argv,"--mdir", help="matrix dir")
argv <- add_argument(argv,"--rds",help="exist rds")
argv <- add_argument(argv,"--step",help="steps to do.all steps are 123456789mono",default="12345")
argv <- add_argument(argv,"--gene_number",help="use top n hvg.genes.default=2000",default=2000)
argv <- add_argument(argv,"--type_marker_tsv",help="cell type marker tsv")
argv <- add_argument(argv,"--ident_tsv",help="cluster identity tsv")
argv <- add_argument(argv,"--compare", help="compare group name:G1vsG2,G3vsG4; split by ,")
argv <- add_argument(argv,"--dim", help="use dim 1:dim",default=15)
argv <- add_argument(argv,"--mtmax", help="filter cell percent.mito max",default=0.2)
argv <- add_argument(argv,"--genemin", help="filter nGene  min",default=200)
argv <- add_argument(argv,"--genemax", help="filter nGene  max",default=Inf)
argv <- add_argument(argv,"--umimin", help="filter umi  min",default=0)
argv <- add_argument(argv,"--umimax", help="filter umi  max",default=30000)
argv <- add_argument(argv,"--resolution", help="tSNE resolution",default=0.6)
argv <- add_argument(argv,"--mono_gene",help="monocle order gene list")
argv <- add_argument(argv,"--remove_contamination",help="remove doublets and unknown cluster in step 9",default="N")
argv <- parse_args(argv)

#read args
mdir <- argv$mdir
dim <- argv$dim
resolution <- argv$resolution
outdir <- argv$outdir
type_marker_tsv <- argv$type_marker_tsv
genemin <- argv$genemin
genemax <- argv$genemax
mtmax <- argv$mtmax
rds <- argv$rds
step <- as.character(argv$step)
config <- argv$config
compare <- unlist(strsplit(argv$compare,split=","))
mono_gene <- argv$mono_gene
ident_tsv <- argv$ident_tsv
is_10X <- argv$is_10X
remove_contamination <- argv$remove_contamination
gene_number <- argv$gene_number
origin.cluster <- paste("res.",resolution,sep="")
pWidth <- 1200
pHeight <- 1000

color2 <- c("OrangeRed","SlateBlue3","DarkOrange","GreenYellow","Purple",
"DarkSlateGray","Gold","DarkGreen","DeepPink2","Red4",
"#4682B4","#FFDAB9","#708090","#836FFF","#CDC673",
"#CD9B1D","#FF6EB4","#CDB5CD","#008B8B","#43CD80",
"#6A5ACD","#66CD00","#CDC673","#CDAD00","#CD9B9B","#FF8247","#8B7355","#8B3A62"
                      ,"#68228B","#CDB7B5","#CD853F","#6B8E23","#E6E6FA","#FFDAB9","Green")

color1 <- c("OrangeRed","SlateBlue3","DarkOrange","GreenYellow","Purple","DarkSlateGray","Gold",
  "DeepPink2","Red4","#4682B4","#FFDAB9","#708090","#836FFF","#CDC673","#CD9B1D","#FF6EB4",
  "#CDB5CD","DarkGreen","#008B8B","#43CD80","#483D8B","#66CD00","#CDC673","#CDAD00","#CD9B9B",
  "#FF8247","#8B7355","#8B3A62","#68228B","#CDB7B5","#CD853F","#6B8E23","#696969","#7B68EE",
  "#9F79EE","#B0C4DE","#7A378B","#66CDAA","#EEE8AA","#00FF00","#EEA2AD","#A0522D","#000080",
  "#E9967A","#00CDCD","#8B4500","#DDA0DD","#EE9572","#EEE9E9","#8B1A1A","#8B8378","#EE9A49",
  "#EECFA1","#8B4726","#8B8878","#EEB4B4","#C1CDCD","#8B7500","#0000FF","#EEEED1","#4F94CD",
  "#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B",
  "#7FFFD4","#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B",
  "#1C86EE","#CDCD00","#473C8B","#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9",
  "#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060","#FF6347","#FF7F50","#CD0000",
  "#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32","#FF00FF","#2E8B57","#CD96CD","#48D1CC",
  "#9B30FF","#1E90FF","#CDB5CD","#191970","#E8E8E8","#FFDAB9")

# read 10X gz
read10X_gz <- function (data.dir = NULL, gene.column = 2) 
{
  full.data <- list()
  for (i in seq_along(data.dir)) {
    run <- data.dir[i]
    if (!dir.exists(run)) {
      stop("Directory provided does not exist")
    }
    if (!grepl("\\/$", run)) {
      run <- paste(run, "/", sep = "")
    }
    barcode.loc <- paste0(run, "barcodes.tsv.gz")
    gene.loc <- paste0(run, "features.tsv.gz")
    matrix.loc <- paste0(run, "matrix.mtx.gz")
    if (!file.exists(barcode.loc)) {
      stop("Barcode file missing")
    }
    if (!file.exists(gene.loc)) {
      stop("Gene name file missing")
    }
    if (!file.exists(matrix.loc)) {
      stop("Expression matrix file missing")
    }
    data <- readMM(file = gzfile(matrix.loc))
    cell.names <- readLines(gzfile(barcode.loc))
    gene.names <- readLines(gzfile(gene.loc))
    if (all(grepl(pattern = "\\-1$", x = cell.names))) {
      cell.names <- as.vector(x = as.character(x = sapply(X = cell.names, 
                                                          FUN = ExtractField, field = 1, delim = "-")))
    }
    rownames(x = data) <- make.unique(names = as.character(x = sapply(X = gene.names, 
                                                                      FUN = ExtractField, field = gene.column, delim = "\\t")))
    if (is.null(x = names(x = data.dir))) {
      if (i < 2) {
        colnames(x = data) <- cell.names
      }
      else {
        colnames(x = data) <- paste0(i, "_", cell.names)
      }
    }
    else {
      colnames(x = data) <- paste0(names(x = data.dir)[i], 
                                   "_", cell.names)
    }
    full.data <- append(x = full.data, values = data)
  }
  full.data <- do.call(cbind, full.data)
  return(full.data)
}

#read config file

tryCatch({
  file.config <- read.table(config,sep="\t",stringsAsFactors=FALSE)
  files <- file.config[,1,drop=T]
  samples <- file.config[,2,drop=T]
  library_names <- sapply(files,function(x) unlist(strsplit(x,"\\."))[[1]])
  names(samples) <- library_names
  print(files)
},error=function(e){cat("There is an error(if rds file provided,this error can be ignored):",conditionMessage(e),"\n")})

tryCatch({
  groups <- file.config[,3]
  names(groups) <- samples
},error=function(e){cat("There is an error(if do not run group compare,this error can be ignored):",conditionMessage(e),"\n")})

#check if all_data exist
checkall_data <- function(){
    if (!exists("all_data")){
        print ("Import RDS file...")
        all_data <- readRDS(rds)
        print ("Import done.")
    }
        return (all_data)
}

#create dir
setwd(outdir)
dir.create("rds")
dir.create("png")
dir.create("csv")
dir.create("pdf")


library(Seurat)
library(dplyr)
############################################## 1.merge 
if (grepl(1,step))
{
print("merge...")
setwd(mdir)
library_number <- length(files)

if (is_10X == "no"){
  s1_matrix <- read.table(gzfile(files[1]),sep="\t",header=TRUE,row.names=1)
  all_data <- CreateSeuratObject(raw.data = s1_matrix,project=library_names[[1]])
} else {
  s1_data <- read10X_gz(data.dir=files[1])
  all_data <- CreateSeuratObject(raw.data = s1_data,project=library_names[[1]])
}

if (library_number>1){
  for (index in c(2:library_number)){
      if (is_10X == "no"){
        new_matrix <- read.table(gzfile(files[index]),sep="\t",header=TRUE,row.names=1)
      } else {
        new_matrix <- read10X_gz(data.dir=files[index])
      }

      library_name <-  library_names[[index]]
      new_obj <- CreateSeuratObject(raw.data = new_matrix,project=library_name)
      all_data <- MergeSeurat(object1 = all_data, object2 = new_obj, add.cell.id2 = library_name,do.normalize=FALSE)
  }
}

setwd(outdir)
saveRDS(all_data,"rds/all_data_raw.rds")
}



################################################# 2.pca 
if (grepl(2,step))
{
print("pca...")
all_data <- checkall_data()
setwd(outdir)
# QC and selecting cells for further analysis
mito.genes <- grep(pattern = "^MT-", x = rownames(x = all_data@data), value = TRUE)
percent.mito <- Matrix::colSums(all_data@raw.data[mito.genes,])/Matrix::colSums(all_data@raw.data)
all_data <- AddMetaData(object = all_data, metadata = percent.mito, col.name = "percent.mito")

png("png/nGene_nUMI_mito.png",,width=pWidth,height =pHeight)
VlnPlot(object = all_data, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
dev.off()

# Plot rare subset of cells with an outlier level of high mitochondrial percentage and also low UMI content
png("png/outlier_mito_UMI.png",,width=pWidth,height =pHeight)
par(mfrow = c(1, 2))
GenePlot(object = all_data, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = all_data, gene1 = "nUMI", gene2 = "nGene")
dev.off()

# filter
print (paste("Raw cell: ",dim(all_data@raw.data)[2]))

all_data <- FilterCells(object = all_data, subset.names = c("percent.mito"), low.thresholds = c(-Inf), high.thresholds = c(mtmax))
print (paste("Cell after filter percent.mito: ",dim(all_data@data)[2]))

all_data <- FilterCells(object = all_data, subset.names = c("nGene"), low.thresholds = c(genemin), high.thresholds = c(genemax))
print (paste("Cell after filter nGene: ",dim(all_data@data)[2]))

print (paste("Remain cell Percent: ",dim(all_data@data)[2]/dim(all_data@raw.data)[2]))



# normalize 
all_data <- NormalizeData(object = all_data, normalization.method = "LogNormalize",scale.factor = 10000)

# Detection of variable genes
all_data <- FindVariableGenes(object = all_data, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.1,  y.cutoff = 1)
print ("var.genes length:")
vargene_len <- length(x = all_data@var.genes)
print (vargene_len)


use.gene <- head(rownames(all_data@hvg.info), gene_number)

# Scaling the data and removing unwanted sources of variation,consume too much mem
# You can perform gene scaling on only the HVG, dramaticall_datay improving speed and memory use. Since dimensional reduction is
# run only on HVG, this will not affect downstream results.


all_data <- ScaleData(object = all_data,vars.to.regress = c("nUMI", "percent.mito"),genes.use =use.gene)

# Perform linear dimensional reduction
all_data <- RunPCA(object = all_data, pc.genes = use.gene, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
# A more ad hoc method for determining which PCs to use is to look at a plot of the standard deviations of the principle 
# components and draw your cutoff where there is a clear elbow in the graph. 
png("png/PCElbow.png",width=pWidth,height =pHeight)
print (PCElbowPlot(object = all_data))
dev.off()

#saveRDS(all_data,"rds/all_PCA.rds")
}



################################################# 3.tsne 
if (grepl(3,step))
{
print("tsne...")
all_data <- checkall_data()
use.gene <- head(rownames(all_data@hvg.info), gene_number)
setwd(outdir)
# save.SNN = T saves the SNN so that the clustering algorithm can be rerun
# using the same graph but with a different resolution value (see docs for
# full details
all_data <- FindClusters(object = all_data, reduction.type = "pca", dims.use = 1:dim, resolution = resolution, print.output = 0, save.SNN = TRUE,force.recalc = TRUE)

# Run Non-linear dimensional reduction (tSNE)
all_data <- RunTSNE(object = all_data, dims.use = 1:dim, do.fast = TRUE)

png("png/TSNE_origin.png",,width=pWidth,height =pHeight)
print (TSNEPlot(object = all_data,do.label = TRUE,pt.size=0.5))
dev.off()

pdf("pdf/TSNE_origin.pdf")
print (TSNEPlot(object = all_data,do.label = TRUE,pt.size=0.5))
dev.off()

#saveRDS(all_data, "rds/all_TSNE.rds")
}



################################################# 4.sample_plot 
if (grepl(4,step))
{
print("sample plot...")
all_data <- checkall_data()
setwd(outdir)
# if "samples" exist(config file exist)
if (exists("samples")){
  all_data@meta.data$samples <- samples[all_data@meta.data$orig.ident]
}
saveRDS(all_data, "rds/all_TSNE_samples.rds")
png("png/sample_tsne.png",width=pWidth,height =pHeight)
print (TSNEPlot(all_data,group.by="samples"))
dev.off()

pdf("pdf/sample_tsne.pdf")
print (TSNEPlot(all_data,group.by="samples",do.return=TRUE,pt.size=0.5)+theme(legend.position = "bottom"))
dev.off()



library(reshape2)
freq_table <- prop.table(x=table(all_data@ident,all_data@meta.data[,"samples"]),margin=2)
write.csv(freq_table,"csv/sample_ident_prop.csv")

dat <- melt(freq_table,varnames=c("type","sample"),value.name = "percent" )
dat$type <- as.character(dat$type)
png("png/sample_ident_barplot.png",width=pWidth,height =pHeight)
print (ggplot(dat, aes(x = sample, y = percent,fill = type )) + geom_bar(stat = 'identity') + coord_flip() +  scale_fill_manual(values=color1) ) 
dev.off()

pdf("pdf/sample_ident_barplot.pdf")
print (ggplot(dat, aes(x = sample, y = percent,fill = type )) + geom_bar(stat = 'identity') + coord_flip() +  scale_fill_manual(values=color1))
dev.off()
}


################################################# 5.marker_gene and go

#GO function
do_go <- function(genes){
    library(clusterProfiler)
    library(org.Hs.eg.db)
    converted <- select(org.Hs.eg.db, 
       keys = genes,
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL")

    gene_id <- na.omit(converted$ENTREZID)

    ego <- enrichGO(gene          = gene_id,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
        readable      = TRUE)
    return (ego)
}


if (grepl(5,step))
{
print("marker gene...")
all_data <- checkall_data()
use.gene <- head(rownames(all_data@hvg.info), gene_number)
all_data.markers <- FindAllMarkers(object = all_data, genes.use = use.gene)
#all_data.markers <- rownames_to_column(all_data.markers,"genes")
write.csv(all_data.markers,"csv/all_markers.csv")
top_10 <- all_data.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
write.csv(top_10,"csv/top_10_markers.csv")

for (c_cluster in unique(all_data.markers$cluster)){
    c_markers <- filter(all_data.markers,cluster==c_cluster)
    file_name <- paste0("./csv/",c_cluster,"_markers.tsv")
    write_tsv(c_markers,file_name)
}

## go
print("go...")
library(clusterProfiler)
library(org.Hs.eg.db)
setwd(outdir)
png_dir <- paste(outdir,"/png/go/",sep="")
csv_dir <- paste(outdir,"/csv/",sep="")
dir.create(png_dir)
dir.create(csv_dir)

for (cluster in unique(all_data.markers$cluster)){
  markers <- as.vector(all_data.markers[all_data.markers$cluster==cluster & all_data.markers$avg_logFC >0,]$gene)
  tryCatch({
    ego <- do_go(markers)

    result_file <- paste(csv_dir,cluster,"_go.csv",sep="")
    write.csv(ego@result,result_file)
  
    dot_png <- paste(png_dir,cluster,"_go_dot.png",sep="")
    png(dot_png,width=pWidth,height =pHeight)
    print (dotplot(ego,showCategory=20))
    dev.off()
  
    bar_png <- paste(png_dir,cluster,"_go_bar.png",sep="")
    png(bar_png,width=pWidth,height =pHeight)
    print (barplot(ego,showCategory=20))
    dev.off()
  
    },error=function(e){cat("There is an error:",conditionMessage(e),"\n")})
}

}



################################################# 6.type_marker
if (grepl(6,step))
{

print("type marker...")
library(tidyverse)

all_data <- checkall_data()
setwd(outdir)
png_dir <- paste(outdir,"/png/","cell_type_marker/",sep="")
pdf_dir <- paste0(outdir,"/pdf/","cell_type_marker/")
auto_dir <- paste0(outdir,"/auto_assign/")
auto_png <- paste0(auto_dir,"png/")
dir.create(png_dir)
dir.create(pdf_dir)
dir.create(auto_dir)
dir.create(auto_png)

# read file
marker_file <- read.table(type_marker_tsv,header = TRUE,sep="\t",row.names = 1,stringsAsFactors=FALSE)
cell_name <- rownames(marker_file)
n_cell_name <- length(cell_name)


# reset
origin.cluster <- paste("res.",resolution,sep="")
all_data <- SetAllIdent(object = all_data, id = origin.cluster)
clusters <- sort(unique(all_data@ident))

#auto
setwd(auto_dir)
c = 0 
for (cluster in clusters){
  index = 0
  for (cell in cell_name){
    index = index + 1
    features = unlist(strsplit(marker_file[index,1],","))
    for (feature in features){
      tryCatch({
        dat <- FindMarkers(all_data,genes.use=feature,ident.1=cluster,min.pct = 0,logfc.threshold = -Inf)
        dat$cell_type <- cell
        dat$cluster <- cluster
        dat <- rownames_to_column(dat,var="gene")
        if (c==0){
          all_dat <- dat
          c = c + 1
        } else {
          all_dat <- rbind(all_dat,dat)
          }
        }
        ,error=function(e){print(paste0(feature," not found!")) })
    }
  }
}

all_dat <- mutate(all_dat,pct.diff=pct.1-pct.2)
write_tsv(all_dat,"type_marker_exp.tsv")

# auto plot
library(RColorBrewer)
exp <- all_dat
a <- group_by(exp,cluster,cell_type)
setwd(auto_png)
for (cluster in sort(unique(a$cluster))){
  c = a[a$cluster==cluster,]

  png(paste0(cluster,"_pctdiff.png"),width=1200,height=1000)
  p1 <- ggplot(c,aes(x=interaction(gene,cell_type),y=pct.diff,fill=cell_type)) +geom_bar(stat="identity")+  coord_flip()

  p1 <- p1 + scale_fill_manual(values=color2)
  
  print (p1)
  dev.off()

  png(paste0(cluster,"_logfc.png"),width=1200,height=1000)
  p2 <- ggplot(c,aes(x=interaction(gene,cell_type),avg_logFC,fill=cell_type)) +geom_bar(stat="identity")+    coord_flip()

  p2 <- p2 + scale_fill_manual(values=color2)
  
  print (p2)
  dev.off()
}

# auto assign
setwd(auto_dir)
as <- summarize(a,avg_pct.diff=mean(pct.diff),avg_logfc=mean(avg_logFC))    
as1 <- group_by(ungroup(as),cluster)
as1 <- mutate(as1,pct_rank = rank(avg_pct.diff),
              logfc_rank= rank(avg_logfc),total_rank=pct_rank+logfc_rank)
as2 <- as1 %>% ungroup %>% group_by(cluster) %>% 
  filter(total_rank==max(total_rank)) %>% arrange(as.numeric(cluster))
write_tsv(as2,"auto_cluster_type.tsv")


# feature and vln
setwd(outdir)
index = 0 
for (cell in cell_name){
  index = index + 1
  features = unlist(strsplit(marker_file[index,1],","))
  tryCatch({
    png_path <- paste(png_dir,cell,"_feature.png",sep="")
    png(png_path,width=pWidth,height =pHeight)
    print (FeaturePlot(object = all_data, features.plot = features,cols.use = c("grey", "blue")))
    dev.off()

    pdf_path <- paste(pdf_dir,cell,"_feature.pdf",sep="")
    pdf(pdf_path)
    print (FeaturePlot(object = all_data, features.plot = features,cols.use = c("grey", "blue")))
    dev.off()

    p <- VlnPlot(object = all_data, features.plot = features,do.return=T)
    if (length(features) == 1){
      p <- p + coord_flip()
    }

    png_path <- paste(png_dir,cell,"_vln.png",sep="")
    png(png_path,width=pWidth,height =pHeight)
    print (p)
    dev.off()

    pdf_path <- paste(pdf_dir,cell,"_vln.pdf",sep="")
    pdf(pdf_path)
    print (p)
    dev.off()

  },error=function(e){cat("There is an error:",conditionMessage(e),"\n")})  
}

}




################################################# 7.group compare diff marker
if (grepl(7,step))
{
print("group compare.diff marker...")
all_data <- checkall_data()
setwd(outdir)

tsne_idents <- levels(all_data@ident)
all_data@meta.data$groups <- groups[all_data@meta.data$samples]
all_data@meta.data$ident.group <- paste0(all_data@ident,"_",all_data@meta.data$groups)
print ("saving rds file...")
#saveRDS(all_data,"rds/all_TSNE_groups.rds")
print ("done")

all_data <- SetAllIdent(all_data, id = "ident.group")

dir.create("compare_diff")
setwd("compare_diff")

for (cm in compare){
  comparegroup <- unlist(strsplit(cm,split="vs"))
  dir.create(cm)
  for (l in tsne_idents){
      s1<-paste(l,'_',comparegroup[1],sep='')
      s2<-paste(l,'_',comparegroup[2],sep='')
      tryCatch({
        cm_markers <- FindMarkers(all_data,  ident.1 = s1, ident.2 = s2,print.bar = FALSE,
          min.diff.pct=0.1,genes.use = use.gene)
        output_dir <- paste(outdir,'/compare_diff/',cm,"/",sep="")
        file_name <- paste(s1,"_vs_",s2,".csv",sep="")
        output_file <- paste(output_dir,file_name,sep="")
        write.csv(cm_markers,file=output_file,sep='\t',quote=F,row.names=T)




      },error=function(e){cat("There is an error:",conditionMessage(e),"\n")})
  }
}

}



################################################# 8.group compare conserve marker
if (grepl(8,step))
{
  print("group compare.conserve marker...")
  all_data <- checkall_data()
  setwd(outdir) 
  #Do not rerun if step 7 exist
  if (!grepl(7,step)){
    tsne_idents <- levels(all_data@ident)
    all_data@meta.data$groups <- groups[all_data@meta.data$samples]
    all_data@meta.data$ident.group <- paste0(all_data@ident,"_",all_data@meta.data$groups)
    saveRDS(all_data,"rds/all_TSNE_groups.rds")
  }

  dir.create("compare_conserve")
  setwd("compare_conserve")

  for (cm in compare){
  comparegroup <- unlist(strsplit(cm,split="vs"))
  all_data <- SetAllIdent(all_data, id = "groups")
  new_data <- SubsetData(all_data,ident.use=comparegroup)
  origin.cluster <- paste("res.",resolution,sep="")
  new_data <- SetAllIdent(new_data, id = origin.cluster)

  dir.create(cm)
    for (l in tsne_idents){
      tryCatch({
        s1<-paste(l,'_',comparegroup[1],sep='')
        s2<-paste(l,'_',comparegroup[2],sep='')

        cm_markers <- FindConservedMarkers(new_data, ident.1 = l, grouping.var = "groups" )
        output_dir <- paste(outdir,'/compare_conserve/',cm,"/",sep="")
        file_name <- paste(s1,"_vs_",s2,sep="")
        output_file <- paste(output_dir,file_name,sep="")
        write.table(cm_markers,file=output_file,sep='\t',quote=F,row.names=T)
        },error=function(e){cat("There is an error:",conditionMessage(e),"\n")})
    }
  }

}




################################################# 9.assign cell ident
if (grepl(9,step))
{
print("assign cell ident...")
all_data <- checkall_data()
new_rds_name <- "rds/new_ident.rds"

setwd(outdir)
cell_ident_file <- read.table(ident_tsv,header = TRUE,sep="\t",stringsAsFactors=FALSE)
current_ident <- cell_ident_file[,1]
new_ident <- cell_ident_file[,2]


if (TRUE %in% duplicated(current_ident))
{
  stop ("duplicated cluster names")
}

if (remove_contamination=="Y"){
  bool <- grepl("unknown|doublet",new_ident)
  current_ident <- current_ident[!bool]
  new_ident <- new_ident[!bool]

  all_data <- SubsetData(all_data,ident.use=current_ident)
  new_rds_name <- "rds/new_ident_rm.rds"
}

c_ident <- 0:(length(current_ident)-1)
c_ident <- paste0("C",c_ident)
new_ident <- paste(c_ident,new_ident,sep="_")


origin.cluster <- paste("res.",resolution,sep="")
all_data <- SetAllIdent(object = all_data, id = origin.cluster)

all_data@ident <- plyr::mapvalues(x = all_data@ident, from = current_ident, to = new_ident)
all_data <- StashIdent(object = all_data, save.name = "new_ident")

pdf("pdf/new_ident_tsne.pdf")
print (TSNEPlot(all_data,do.label=TRUE,do.return=TRUE,pt.size=0.5,colors.use=color1)+theme(legend.position = "bottom"))
dev.off()

png("png/new_ident_tsne.png",,width=pWidth,height =pHeight)
print (TSNEPlot(all_data,do.label=TRUE,do.return=TRUE,pt.size=0.5))
dev.off()

#new_prop_ident
library(reshape2)
sample_col = "orig.ident"
cols = colnames(all_data@meta.data)
if ("samples" %in% cols) {
    sample_col="samples"}else if ("sample" %in% cols) {
    sample_col="sample"}
freq_table <- prop.table(x=table(all_data@ident,all_data@meta.data[,sample_col]),margin=2)
dat <- melt(freq_table,varnames=c("type","sample"),value.name = "percent" )
dat$type <- as.character(dat$type)
dat$type <- factor(dat$type,levels = new_ident)

pdf("pdf/sample_ident_barplot_new.pdf")
print (ggplot(dat, aes(x = sample, y = percent,fill = type )) + 
  geom_bar(stat = 'identity') + coord_flip() +  scale_fill_manual(values=color2) ) 
dev.off()

library(pheatmap)
pdf("pdf/sample_ident_heatmap.pdf")
print (pheatmap(freq_table,display_numbers = TRUE,number_format ="%.2f",fontsize_number=9,fontsize=11,main="sample_identity_heatmap") )
dev.off()

print ("save rds")
saveRDS(all_data,new_rds_name)
print ("save done")
}



################################################# 10.monocle
if (grepl("mono",step))
{
print("monocle...")
library(monocle)
all_data <- checkall_data()

cols <- colnames(all_data@meta.data)
if (!is.na(match("clusters",cols))){
  clusters <- "clusters"
} else {
  clusters <- paste0("res.",resolution)
}
print (clusters)

setwd(outdir)
mono <- importCDS(all_data)
mono <- estimateSizeFactors(mono)

# gene list TODO!
if (is.na(mono_gene)){
  hv.genes <- head(rownames(all_data@hvg.info), 200)
  mono_order_genelist <- hv.genes
  } else {
  mono_order_genelist <- mono_gene
  }

# trajectory
mono <- setOrderingFilter(mono, mono_order_genelist)
mono <- reduceDimension(mono, max_components = 2,method = 'DDRTree')
mono <- orderCells(mono)
print ("ordering cells")
saveRDS(mono,"./rds/monocle.rds")

pdf("pdf/monocle.pdf")
print (plot_cell_trajectory(mono, color_by = clusters))
dev.off()
}



################################################# 11.cor
if (grepl("cor",step))
{
print("cor...")
all_data <- checkall_data()

lm22 <- read_csv("/SGRNJ/Database/download/bulk/LM22/LM22.csv")
genes <- lm22$`Gene symbol`

all_data <- SetAllIdent(all_data,id=origin.cluster)
avg_exp <- AverageExpression(all_data)
genes1 <- rownames(avg_exp)
genes2 <- intersect(genes,genes1)

m <- avg_exp[genes2,]
l <- filter(lm22,`Gene symbol` %in% genes2)

clusters <- colnames(m)
cells <- colnames(l)[-1]

result <- tibble(cluster=character(),ident=character(),cor=numeric())

for (cluster in clusters){
  exp <- m[,cluster]
  cors <- c()
  for (cell in cells){
    exp1 <- l[,cell]
    cors[cell] <- cor(exp,exp1,method="spearman")
  }
  max_cor <- max(cors)
  max_ident <- names(cors[cors==max_cor])
  result <- add_row(result,cluster=cluster,ident=max_ident,cor=max_cor)
}

trans <- result$ident
names(trans) <- result$cluster
all_data@meta.data$cor_ident <- trans[all_data@meta.data[,origin.cluster]]
all_data <- SetAllIdent(all_data,id="cor_ident")
TSNEPlot(all_data,,pt.size=0.5)

}


#######################
