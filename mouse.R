x=colnames(m.imputed)
library(kazutils)
mprot_data <- list(vals = m.imputed, 
                   rowAnn = data.frame(ProteinIDs=m.ibaq.original$X, 
                                    ProteinNames=m.ibaq.original$pname, 
                                    GeneNames=get_nth_part(rownames(m.imputed), "\\.", 1), 
                                    GeneNamesUnique=rownames(m.imputed)),
                   colAnn = data.frame(cellType=get_nth_part(x,"_",1), 
                                    hormoneTreatment=get_nth_part(x,"_",2), row.names = x))
                   
mprot_data$colAnn$cellType = plyr::mapvalues(mprot_data$colAnn$cellType, c("B", "LP", "LM"), c("BC", "LP", "ML"))

mprot_data$colors <- list(cellType=c("BC"="red", "LM"="darkblue", "LP"="skyblue"), 
                          hormoneTreatment = c("EP"="darkgoldenrod1","E"="darkgoldenrod"))
saveRDS(mprot_data,"data/mouse_proteome.rds")
mprot_data <- readRDS("data/mouse_proteome.rds")
