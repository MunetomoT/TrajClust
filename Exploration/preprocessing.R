# Load required packages
library(Seurat)
library (tidyverse)
library(dplyr)
library(Matrix)
library(ggplot2)
library(cowplot)
library (naniar)
library(sctransform)
library(scRepertoire)
library(stringr)

# Adjust memory requirements
options (future.globals.maxSize = 4000 * 1024^6)

# Read count matrix
data_dir <- 'filtered_feature_bc_matrix'
list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
data <- Read10X(data.dir = data_dir)

# Remove T cell associated genes
x <- data
rownames (x)
x[grep ("^Trav", rownames (x)),]
y=grep ("^Trav", rownames (x))
x2<-x[-y,]
grep ("^Trbv", rownames (x2))
y=grep ("^Trbv", rownames (x2))
x3<- x2[-y,]
grep ("^Trav", rownames (x3))
grep ("^Trbv", rownames (x3))
data <- x3
rm (x, x2, x3, y)
rm (data_dir)
mito.features <- grep(pattern = "^mt-", x = rownames(x =data), value = TRUE)

# Save into Seurat object
All <- CreateSeuratObject(counts = data)
suppressMessages(library(scRepertoire))

# Load metadata for each batch (abbreivated to show only csv1)
csv1 <- read.csv("csv1.csv",
                 header=TRUE,stringsAsFactors = FALSE)
csv1 <- separate (csv1, barcode,into = c("barcode","sample"))
csv1$sample <- "1"
csv1 <- unite (csv1, barcode,barcode,sample, sep ="-" )

# Process metadata to obtain clone association for each cell
contig_list <- list(csv1, csv2,csv3, csv4,csv5, csv6,csv7, csv8,csv9, csv10,csv11, csv12,csv13, csv14,csv15, csv16,csv17)
combined <- combineTCR(contig_list, 
                       samples = c("M1S", "M1T", "M2S", "M2T", "M3S","M3T","M4S","M4T","M5S","M5T","M6S","M6T","M1C","M4C","M2C","M7S","M7T"), 
                       ID = c("S", "T", "S", "T", "S", "T",
                              "S", "T", "S", "T", "S", "T",
                              "C","C","C",
                              "S","T"),
                       filterMulti = TRUE, cells ="T-AB")
for (i in seq_along(combined)) {
  combined[[i]] <- stripBarcode(combined[[i]], column = 1, connector = "_", num_connects = 3)
} 
All.clean <- combineExpression(combined, All, cloneCall="gene", groupBy = "none", 
                                 cloneTypes = c(None = 0, Single = 1, Small = 5, Medium = 20, Large = 100,Hyperexpanded = 2500))
rm(All)

# Assign sample number from cell barcodes
temp <- list()
for(i in All.clean@meta.data[["barcode"]]){temp <- append(temp, strsplit(i, "-|\\s")[[1]][2])}
All.clean[["Sample"]] <- unlist(temp, use.names=FALSE)
temp2 <- list()
for(i in All.clean@meta.data[["Sample"]]){temp2 <- append(temp2, if (is.null(i)) {0} else if
                                                          (i=="1" || i=="2" || i=="13") {1} else if 
                                                          (i=="3" || i=="4" || i=="15") {2} else if
                                                          (i=="5" || i=="6") {3} else if
                                                          (i=="7" || i=="8" || i=="14") {4} else if
                                                          (i=="9" || i=="10") {5} else if
                                                          (i=="11" || i=="12") {6} else {7}
                                                          )}

# Define percent.mito
percent.mito <- Matrix::colSums (x = GetAssayData(object = All.clean, slot = 'counts')[mito.features,]) / Matrix::colSums(x = GetAssayData(object = All.clean, slot = 'counts'))
All.clean[['percent.mito']] <- percent.mito

# Look at QC plots
VlnPlot(All.clean, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
plot1 <- FeatureScatter(All.clean, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(All.clean, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Clean data
All.clean <- subset(x = All.clean, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mito < 0.1) 

saveRDS (All.clean, file = "realData")
project <- readRDS("realData")

# Remove T cells that have NA
remove <- WhichCells(project, project@meta.data$CTnt == "N.A")
project <- project[,!colnames(project) %in% remove]

# Add clone count metadata to each cell
i <- 1
codes <- as.data.frame(table(project@meta.data$CTnt))
a<-as.list(project@meta.data$CTnt)
for (x in a) {
  a[i] <- as.integer(codes[which(codes$Var1==x),]$Freq)
  i <- i+1
}
a<-unlist(a)
project <- AddMetaData(project, a, col.name="CloneCount")

# Remove clones with incomplete TCR sequences
a <- startsWith(project@meta.data[["CTnt"]], "NA")
subset <- project[, !a]
b <- endsWith(subset@meta.data[["CTnt"]], "NA")
subset <- subset[, !b]
rm (project)

# Trim down
DefaultAssay(subset) <- "RNA"
subset[['integrated']] <- NULL
subset[['SCT']] <- NULL

# Gene markers for defining cell cycle
g2m.genes.m <- c("Hmgb2","Cdk1","Nusap1","Ube2c","Birc5","Tpx2","Top2a","Ndc80","Cks2","Nuf2","Cks1b","Mki67","Tmpo",
                 "Cenpf","Tacc3","Fam64a","Smc4","Ccnb2","Ckap2l","Ckap2","Aurkb","Bub1","Kif11","Anp32e","Tubb4b","Gtse1",
                 "Kif20b","Hjurp","Cdca3","Hn1","Cdc20","Ttk","Cdc25c","Kif2c","Rangap1","Ncapd2","Dlgap5","Cdca2","Cdca8",
                 "Ect2","Kif23","Hmmr","Aurka","Psrc1","Anln","Lbr","Ckap5","Cenpe","Ctcf","Nek2","G2e3","Gas2l3","Cbx5","Cenpa")
s.genes.m <- c("Mcm5","Pcna","Tyms","Fen1","Mcm2","Mcm4","Rrm1","Ung","Gins2","Mcm6","Cdca7","Dtl","Prim1","Uhrf1","Mlf1p","Hells","Rfc2","Rpa2","Nasp",
               "Rad51ap1","Gmnn","Wdr76","Slbp","Ccne2","Ubr7","Pold3","Msh2","Atad2","Rad51","Rrm2","Cdc45","Cdc6","Exo1","Tipin","Dscc1","Blm","Casp8ap2","Usp1",
               "Clspn","Pola1","Chaf1b","Brip1","E2f8")

# Subset to each mouse
M1.CD8 <- subset(subset, Mouse == 1)
M2.CD8 <- subset(subset, Mouse == 2)
M3.CD8 <- subset(subset, Mouse == 3)
M4.CD8 <- subset(subset, Mouse == 4)
M5.CD8 <- subset(subset, Mouse == 5)
M6.CD8 <- subset(subset, Mouse == 6)
M7.CD8 <- subset(subset, Mouse == 7)

# For each mouse, regress cell cycle score and normalise/ logarithmize (abbreviated to show only Mouse 1)

M1.CD8 <- CellCycleScoring(M1.CD8, s.features = s.genes.m, g2m.features = g2m.genes.m, set.ident = TRUE)
M1.CD8$CC.Difference <-M1.CD8$S.Score - M1.CD8$G2M.Score
M1.CD8 <- SCTransform(M1.CD8, variable.features.n=4000, vars.to.regress = c("S.Score", "G2M.Score"), verbose = TRUE)

# Integrate data together
tumour.list <- list (M1.CD8,M2.CD8,M3.CD8,M4.CD8,M5.CD8,M6.CD8, M7.CD8)
rm (M1.CD8,M2.CD8,M3.CD8,M4.CD8,M5.CD8,M6.CD8,M7.CD8)

# Integration of batch (Seurat v3 Integration)
tumour.features <- SelectIntegrationFeatures(object.list = tumour.list, nfeatures = 4000)
tumour.list <- PrepSCTIntegration (object.list = tumour.list, anchor.features = tumour.features, verbose = TRUE)
tumour.anchors <- FindIntegrationAnchors(object.list = tumour.list, normalization.method = "SCT", anchor.features = tumour.features, verbose = TRUE)
rm(tumour.list)
rm(subset)
tumour.integrated <- IntegrateData(anchorset = tumour.anchors, normalization.method = "SCT",
                                   verbose = TRUE)
DefaultAssay(tumour.integrated) <- "integrated"

# Visualise to check results
tumour.integrated <- RunPCA(tumour.integrated, npcs = 50, verbose = TRUE)
ElbowPlot (tumour.integrated, ndims = 50, reduction = "pca")
tumour.integrated <- RunUMAP(tumour.integrated, dims = 1:30, verbose = TRUE)
DimPlot(tumour.integrated, group.by = "DaySite")+DarkTheme()
DimPlot(tumour.integrated, group.by = "Tissue", cols = c("blue2","orange","firebrick1","grey"))+DarkTheme()
tumour.integrated <- FindNeighbors(tumour.integrated, dims = 1:50, verbose = TRUE)
tumour.integrated <- FindClusters(tumour.integrated, verbose = TRUE, resolution=0.3)

# Save files
saveRDS(tumour.integrated, file='RealDataset.rds')
SaveH5Seurat(pro, filename = "RealDataset.h5Seurat", overwrite = TRUE)
Convert("RealDataset.h5Seurat", dest = "RealDataset.h5ad", overwrite=TRUE)
