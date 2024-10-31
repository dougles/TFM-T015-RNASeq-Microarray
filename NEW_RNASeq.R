# instalacion de Bioconductor
install.packages("BiocManager") # gestiona paquetes de bioconductor
BiocManager::install("maftools")
BiocManager::install("TCGAbiolinks", force=TRUE)
BiocManager::install("DESeq2")

library(BiocManager)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2) # F change, Q valor

#Obtener datos base
base_query_TCGA <- GDCquery(project = "TCGA-BRCA",
                            data.category = "Transcriptome Profiling",
                            experimental.strategy = "RNA-Seq",
                            workflow.type = "STAR - Counts",
                            access = "open",
                            sample.type = c("Primary Tumor","Solid Tissue Normal"))
base_query_TCGA <- getResults(base_query_TCGA)


#Preparar infoStudy
#Inicializar variables 
df_final_infoStudy <- data.frame()
maxNumArr <- 0
normal_breast_result<-c()
tumor_breast_result<-c()

#Loop ejecutar desde aqui para la seleccion dinamica 
for (x in 1:20) {

#Preparar infoStudy    
df_pt <- base_query_TCGA[base_query_TCGA$sample_type =="Primary Tumor",]
df_pt_15 <- df_pt[sample(nrow(df_pt), size = 15), ]

df_stn <- base_query_TCGA[base_query_TCGA$sample_type =="Solid Tissue Normal",]
df_stn_15 <- df_stn[sample(nrow(df_stn), size = 15), ]

df_pt_result <- df_pt_15[c("cases","sample_type")]
df_pt_result$name_id <- c("Tumor_breast_1","Tumor_breast_2","Tumor_breast_3","Tumor_breast_4",
                       "Tumor_breast_5","Tumor_breast_6","Tumor_breast_7","Tumor_breast_8",
                       "Tumor_breast_9","Tumor_breast_10","Tumor_breast_11","Tumor_breast_12",
                       "Tumor_breast_13","Tumor_breast_14","Tumor_breast_15")
df_stn_result <- df_stn_15[c("cases","sample_type")]
df_stn_result$name_id <- c("Normal_breast_1","Normal_breast_2","Normal_breast_3","Normal_breast_4","Normal_breast_5",
                         "Normal_breast_6","Normal_breast_7","Normal_breast_8","Normal_breast_9","Normal_breast_10",
                         "Normal_breast_11","Normal_breast_12","Normal_breast_13","Normal_breast_14","Normal_breast_15")
infoStudy <- rbind(df_stn_result, df_pt_result)


# lanzar consulta para: gene expression data; estrategia experimental: RNA-Seq, flujo de trabajo de analisis: STAR-counts, acceso: open
query_TCGA_final <- GDCquery(project = "TCGA-BRCA",
                       data.category = "Transcriptome Profiling",
                       experimental.strategy = "RNA-Seq",
                       workflow.type = "STAR - Counts",
                       access = "open",
                       barcode = as.vector(infoStudy$cases)) # <-lista de obserciones info study
# muestra la consulta por pantalla
getResults(query_TCGA_final)
#Descarga
GDCdownload(query_TCGA_final)
#preprocesamiento de los datos descargados (counts)
SKCM.counts <- GDCprepare(query = query_TCGA_final, summarizedExperiment = TRUE)
#rm(query_TCGA)
counts_data <- assay(SKCM.counts)


# PREPROCESAMIENTO DE LOS DATOS
# obtener el nombre del gen id
gene_names <-SKCM.counts@rowRanges@elementMetadata@listData$gene_name
rownames(counts_data) <- gene_names

# convert counts to DGEList object
dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = infoStudy, design = ~ sample_type)

# prefiltering: removing rows with low gene counts
# keeping rown that have at least 10 reads total
keep <- rowSums(counts(dds)) >=10  ## Eliminando genes de conteso irrelevante
dds <- dds[keep,]

countdata <- counts(dds) # para obtener el log2
countdata_log2 <- log2(countdata +1)

# CALCULOS DE LOS GENES DIFERENCIALMENTE EXPRESADOS
# set the factor level
dds$sample_type <- relevel(dds$sample_type, ref = 'Solid Tissue Normal')
# Run DESeq
dds <- DESeq(dds) # datos para el analisis diferencial
res <- results(dds)
summary(res)

fold.change.Breast <- res$log2FoldChange # para obtener FC: la magnitud del cambio en la expresion de dos condiciones (diferencia logaritmica)
adj.Pval.Breast <- res$padj # pvalor ajustado: mide la significancia estadistica del cambio en la expresion genica
genes.ids.Breast <- rownames(res)

# metodo combinado para los genes activados y reprimidos
activated.genes.breast.1 <- genes.ids.Breast[fold.change.Breast > 3 & adj.Pval.Breast < 0.05]
repressed.genes.breast.1 <- genes.ids.Breast[fold.change.Breast < -3 & adj.Pval.Breast < 0.05]

# plot(fold.change.Breast)
# plot(adj.Pval.Breast)

length(activated.genes.breast.1)
length(repressed.genes.breast.1)

names(fold.change.Breast) <- genes.ids.Breast
log.padj.breast <- -log10(adj.Pval.Breast)
names(log.padj.breast) <- genes.ids.Breast

# plot(fold.change.Breast,log.padj.breast, ylab="-log10(p value)", xlab= "log2 fold change", pch=19, cex=0.5,
#      col = "grey", xlim=c(-15,15))
# points(fold.change.Breast[activated.genes.breast.1],log.padj.breast[activated.genes.breast.1], 
#        pch = 19, cex = 0.5, col = "red")
# points(fold.change.Breast[repressed.genes.breast.1],log.padj.breast[repressed.genes.breast.1], 
#        pch = 19, cex = 0.5, col = "blue")
# 


# OBTENER DOS DATASET NORMAL Y TUMOR CON GENES DIFERENCIALMENTE EXPRESADOS
breast.all.DEG <- genes.ids.Breast[(fold.change.Breast > 3 | fold.change.Breast < -3) & adj.Pval.Breast < 0.05]
length(breast.all.DEG)

colnames(countdata_log2) <-infoStudy$name_id
normal.breast.DEG.table <- countdata_log2[, c("Normal_breast_1", "Normal_breast_8", 
                                              "Normal_breast_2", "Normal_breast_9", 
                                              "Normal_breast_3", "Normal_breast_10", 
                                              "Normal_breast_4", "Normal_breast_11", 
                                              "Normal_breast_5", "Normal_breast_12", 
                                              "Normal_breast_6", "Normal_breast_13", 
                                              "Normal_breast_7", "Normal_breast_14", "Normal_breast_15")]
normal.breast.DEG.table <- normal.breast.DEG.table[rownames(normal.breast.DEG.table) %in% breast.all.DEG,]
normal.breast.DEG.table <- cbind(attr_name = rownames(normal.breast.DEG.table), normal.breast.DEG.table)

tumor.breast.DEG.table <- countdata_log2[, c("Tumor_breast_1", "Tumor_breast_8", 
                                             "Tumor_breast_2", "Tumor_breast_9", 
                                             "Tumor_breast_3", "Tumor_breast_10", 
                                             "Tumor_breast_4", "Tumor_breast_11", 
                                             "Tumor_breast_5", "Tumor_breast_12", 
                                             "Tumor_breast_6", "Tumor_breast_13", 
                                             "Tumor_breast_7", "Tumor_breast_14", "Tumor_breast_15")]

tumor.breast.DEG.table <- tumor.breast.DEG.table[rownames(tumor.breast.DEG.table) %in% breast.all.DEG,]
tumor.breast.DEG.table <- cbind(attr_name = rownames(tumor.breast.DEG.table), tumor.breast.DEG.table)

# write.table(normal.breast.DEG.table, file="../NormalBreast_DEGs_table.csv", sep = "\t", 
#             quote = FALSE, row.names = FALSE, col.names = TRUE)
# write.table(tumor.breast.DEG.table, file="../TumorBreast_DEGs_table.csv", sep = "\t", 
#             quote = FALSE, row.names = FALSE, col.names = TRUE)


  if (length(tumor.breast.DEG.table) > maxNumArr) {
    df_final_infoStudy <- infoStudy
     normal_breast<-normal.breast.DEG.table
     tumor_breast<-tumor.breast.DEG.table
     maxNumArr <- length(tumor.breast.DEG.table)
  }

}


write.table(normal_breast, file="NormalBreast_DEGs_table.csv", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(tumor_breast, file="TumorBreast_DEGs_table.csv", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)

df_final_infoStudy 

write.table(df_final_infoStudy, file="infostudy_DEGs_table.csv", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)

