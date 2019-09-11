DenoisedCounts = DenoisedCounts[-grep("ERCC",rownames(DenoisedCounts)),]
DenoisedCounts = DenoisedCounts[-grep("Hist",rownames(DenoisedCounts)),]
DenoisedCounts = DenoisedCounts[-grep("^Rps",rownames(DenoisedCounts)),]
DenoisedCounts = DenoisedCounts[-grep("^Rpl",rownames(DenoisedCounts)),]
DenoisedCounts = DenoisedCounts[-grep("^Gm",rownames(DenoisedCounts)),]
DenoisedCounts = DenoisedCounts[-grep("^ENSMUSG",rownames(DenoisedCounts)),]
DenoisedCounts = DenoisedCounts[rowSums(DenoisedCounts > 5)>3,]

expr_data = DenoisedCounts[,order(colnames(DenoisedCounts))]
phenotype_clean = sc_sort_info[sc_sort_info$sample_id %in% colnames(expr_data),]
phenotype_clean = phenotype_clean[order(phenotype_clean$sample_id),]
rownames(phenotype_clean) = phenotype_clean$sample_id
phenotype_clean = phenotype_clean[,c("CD16.32", "CD11b", "CD135", "CD127", "CD150", "Sca.1", "cKit")]
phenotype_clean[phenotype_clean<10] = 10



ph_data = phenotype_clean


save(ph_data, expr_data, file="sc_192.RData")