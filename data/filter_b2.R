# Cell numbers:
# > table(b2$celltype)
# Basal1 Basal2    LP1    LP2    LP3    LP4    ML1    ML2    ML3    ML4    ML5
# 689    583       975    857    788    744   1140    414    386    157     75
#
# > table(b2$broadcelltype)
# Mature luminal     Luminal progenitor   Basal
# 2172               3364                 1272

#Prepare data

# b2_backup <- b2
apply(b2@meta.data, 2, function(x){ !is.numeric(x)})
"Genotype" # RM = reduction mammoplasty
"broadcelltype"
"celltype"
"seurat_clusters"

df <- b2@meta.data

df2 <- data.frame(CellType = b2$broadcelltype,
                  CellStates = b2$celltype,
                  Seurat_clusters = b2$seurat_clusters,
                  Seurat_cellCyclePhase = b2$Phase,
                  Source=ifelse(b2$Genotype=="RM", "Reduction Mammoplasty", "NA_"),
                  Patient=plyr::mapvalues(b2$barcode, c("1","10","17","25"), as.character(1:4)),
                  row.names = rownames(df))

df2$CellType <- plyr::mapvalues(df2$CellType, c("Mature luminal", "Luminal progenitor", "Basal"), c("ML","LP","BC"))
saveRDS(sc_data, "data/sc_Seurat_object.rds")

# Import proteome data
df_mat <- read.delim("data/human_MEC_proteome_combat_matrix.txt", row.names = 1)
df_row <- read.delim("data/human_MEC_proteome_row_info.txt", row.names = 1)
df_col <- read.delim("data/human_MEC_proteome_sample_info.txt", row.names = 1)  %>%
  .[,c("cellType", "patientAgeGroup", "hormoneStatus")]

df_col$cellType <- plyr::mapvalues(df_col$cellType, c("LC","LP","BC"), c("ML","LP","BC"))

hprot_data <- list(mat=df_mat, row=df_row, col=df_col)
saveRDS(hprot_data, "data/human_proteome.rds")
