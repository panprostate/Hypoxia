# R script to calculate alternative hypoxia scores for PPCG data
# 19.08.24
# Maria Jakobsdottir <maria.jakobsdottir@manchester.ac.uk>

# Using hypoxia functions where genes are scaled within each sample

#### Libraries ####
library(singscore)
library(qs)
library(SummarizedExperiment)

#### Data ####
rna_dat_unsupervised = qread("/mnt/bmh01-rds/UoOxford_David_W/shared/PPCG_RNA/K10_unsupervised/PPCG_RNAseq_RUVIII_K10.qs")
rna_dat_supervised = qread("/mnt/bmh01-rds/UoOxford_David_W/shared/PPCG_RNA/DEPRECATED_Release_1.3/PPCG_RNAseq_RUVIIIPRPS_Supervised.qs")

#### Gene lists ####
ragnum_gene_list <- structure(list(ENSG = c("ENSG00000105011", "ENSG00000066279", 
                                            "ENSG00000089685", "ENSG00000154473", "ENSG00000138778", "ENSG00000151725", 
                                            "ENSG00000140931", "ENSG00000159147", "ENSG00000143476", "ENSG00000111206", 
                                            "ENSG00000160211", "ENSG00000123485", "ENSG00000073111", "ENSG00000112818", 
                                            "ENSG00000087053", "ENSG00000071539", "ENSG00000122952", "ENSG00000139372", 
                                            "ENSG00000076248", "ENSG00000196419", "ENSG00000148926", "ENSG00000168209", 
                                            "ENSG00000096696", "ENSG00000088340", "ENSG00000135245", "ENSG00000122884", 
                                            "ENSG00000226784", "ENSG00000067225", "ENSG00000177181", "ENSG00000258818", 
                                            "ENSG00000099194", "ENSG00000061656"), 
                                   SYMBOL = c("ASF1B", "ASPM", 
                                              "BIRC5", "BUB3", "CENPE", "CENPU", "CMTM3", "DONSON", "DTL", 
                                              "FOXM1", "G6PD", "HJURP", "MCM2", "MEP1A", "MTMR2", "TRIP13", 
                                              "ZWINT", "TDG", "UNG", "XRCC6", "ADM", "DDIT4", "DSP", "FER1L4", 
                                              "HILPDA", "P4HA1", "PGAM4", "PKM", "RIMKLA", "RNASE4", "SCD", 
                                              "SPAG4")), class = "data.frame", row.names = c(NA, -32L))

west_gene_list <- structure(list(SYMBOL = c("ADAMTS4", "ATF3", "BHLHE40", "BTG2", 
                                            "CSRNP1", "CYR61", "EGR1", "EGR2", "EGR3", "FOSB", "FOSL2", "GEM", 
                                            "JUNB", "KLF10", "KLF6", "LIF", "MCL1", "NR4A3", "PPP1R15A", 
                                            "RHOB", "SELE", "SIK1", "SLC2A14", "SLC2A3", "SOCS3", "THBS1", 
                                            "TIPARP", "ZFP36"), 
                                 ENSG = c("ENSG00000158859", "ENSG00000162772", 
                                          "ENSG00000134107", "ENSG00000159388", "ENSG00000144655", "ENSG00000142871", 
                                          "ENSG00000120738", "ENSG00000122877", "ENSG00000179388", "ENSG00000125740", 
                                          "ENSG00000075426", "ENSG00000164949", "ENSG00000171223", "ENSG00000155090", 
                                          "ENSG00000067082", "ENSG00000128342", "ENSG00000143384", "ENSG00000119508", 
                                          "ENSG00000087074", "ENSG00000143878", "ENSG00000007908", "ENSG00000142178", 
                                          "ENSG00000173262", "ENSG00000059804", "ENSG00000184557", "ENSG00000137801", 
                                          "ENSG00000163659", "ENSG00000128016"), 
                                 Direction = c("Up","Down", "Up", "Down", "Up", "Up", 
                                               "Down", "Up", "Down", "Down", "Up", 
                                               "Down", "Up", "Down", "Up", "Up", 
                                               "Up", "Down", "Down", "Up", "Down", 
                                               "Up", "Up", "Down", "Up", "Down", 
                                               "Down", "Down"), 
                                 Coefficient = c(0.208143369, 
                                                 -0.182281129, 0.145486214, -0.383633755, 0.064308417, 0.40495632, 
                                                 -0.396671892, 0.314855084, -0.298747262, -0.042448357, 0.176084142, 
                                                 -0.138491418, 0.213008876, -0.183908536, 0.208046534, 0.325537409, 
                                                 0.306647257, -0.059548777, -0.006602269, 0.220516453, -0.103118261, 
                                                 0.290011507, 0.196779232, -0.28084907, 0.088841536, -0.162059422, 
                                                 -0.553245162, -0.14218079)), 
                            class = "data.frame", row.names = c(NA,-28L))

buffa_gene_list <- structure(list(SYMBOL = c("VEGFA", "SLC2A1", "PGAM1", "ENO1", "LDHA",
                                             "TPI1", "P4HA1", "MRPS17", "CDKN3", "ADM",
                                             "NDRG1", "TUBB6", "ALDOA", "MIF", "ACOT7",
                                             "MCTS1", "PSRC1", "PSMA7", "ANLN", "TUBA1B",
                                             "MAD2L2", "GPI", "TUBA1C", "MAP7D1", "DDIT4",
                                             "BNIP3", "MRGBP", "HILPDA", "GAPDH", "MRPL13",
                                             "CHCHD2", "YKT6", "PNP", "CORO1C", "SEC61G",
                                             "ANKRD37", "ESRP1", "PFKP", "SHCBP1", "CTSL2",
                                             "KIF20A", "SLC25A32", "UTP11", "SLC16A1", "MRPL15",
                                             "KIF4A", "LRRC42", "PGK1", "HK2", "AK3",
                                             "CA9"),
                                  ENSG = c("ENSG00000112715", "ENSG00000117394", "ENSG00000171314", "ENSG00000074800", "ENSG00000134333",
                                           "ENSG00000111669", "ENSG00000122884", "ENSG00000239789", "ENSG00000100526", "ENSG00000148926",
                                           "ENSG00000104419", "ENSG00000176014", "ENSG00000149925", "ENSG00000240972", "ENSG00000097021",
                                           "ENSG00000232119", "ENSG00000134222", "ENSG00000101182", "ENSG00000011426", "ENSG00000123416",
                                           "ENSG00000116670", "ENSG00000105220", "ENSG00000167553", "ENSG00000116871", "ENSG00000168209",
                                           "ENSG00000176171", "ENSG00000101189", "ENSG00000135245", "ENSG00000111640", "ENSG00000172172",
                                           "ENSG00000106153", "ENSG00000106636", "ENSG00000198805", "ENSG00000110880", "ENSG00000132432",
                                           "ENSG00000186352", "ENSG00000104413", "ENSG00000067057", "ENSG00000171241", "ENSG00000136943",
                                           "ENSG00000112984", "ENSG00000164933", "ENSG00000183520", "ENSG00000155380", "ENSG00000137547",
                                           "ENSG00000090889", "ENSG00000116212", "ENSG00000102144", "ENSG00000159399", "ENSG00000147853",
                                           "ENSG00000107159")),
                             class = "data.frame", 
                             row.names = c(NA,-51L))


#### Functions ####
Ragnum_hypoxia_score_within <- function(expression_matrix, method, annotation, gene_list){
  # helper not in statement
  '%!in%' <- function(x,y)!('%in%'(x,y))
  
  # Stop statements
  if("matrix" %!in% class(expression_matrix)){
    stop("Input data must be an expression matrix")
  }
  
  if(ncol(expression_matrix) < 3){## Think about threshold here
    stop("Not enough samples to calculate signature")
  }
  
  if((annotation %in% c("ENSG","SYMBOL","ID_SYMBOL")) == FALSE){
    stop("Annotation must be one of 'ENSG' or 'SYMBOL' or 'ID_SYMBOL'.")
  }
  
  # Set up output table
  output_data <- as.data.frame(colnames(expression_matrix))
  colnames(output_data) <- c("Sample")
  output_data$Genes_used_Ragnum <- NA
  
  Signature_Genes_Present <- NA
  Signature_Genes_Absent <- NA
  
  ragnum_genes <- gene_list
  
  # Warnings about missing signature genes
  Signature_Genes_Present <- ragnum_genes[ragnum_genes[,annotation] %in% rownames(expression_matrix), annotation]
  Signature_Genes_Absent <- ragnum_genes[ragnum_genes[,annotation] %!in% rownames(expression_matrix), annotation]
  
  if(length(Signature_Genes_Present) < nrow(ragnum_genes)){
    warning(paste("Warning: ",length(Signature_Genes_Absent),"signature genes not present in input matrix; results may be incomplete or unreliable. Missing genes:",Signature_Genes_Absent, sep = " "))
  }
  
  output_data$Genes_used_Ragnum <- toString(Signature_Genes_Present)
  
  # Median centre and scale to unit variance
  expression_matrix_ragnum <- as.matrix(expression_matrix[ Signature_Genes_Present, ])
  expression_matrix_ragnum  <- expression_matrix_ragnum - rowMedians(expression_matrix_ragnum)
  expression_matrix_ragnum_scaled <- (scale((expression_matrix_ragnum), center = F, scale = T))
  
  
  if(method == "Sum"){
    # add column to output table
    output_data$Ragnum_signature_Sum <- NA
    
    # calculate signature
    output_data$Ragnum_signature_Sum <- colSums(expression_matrix_ragnum_scaled)
    
  }else if(method == "Binary"){
    # add column to output table
    output_data$Ragnum_signature_Binary <- NA
    
    for(i in 1:nrow(output_data)){
      ragnum_sig_value_rna <- 0
      for (j in 1:nrow(expression_matrix_ragnum_scaled)) {
        if(expression_matrix_ragnum_scaled[j,i] > 0){
          ragnum_sig_value_rna <- ragnum_sig_value_rna + 1
        }else{
          ragnum_sig_value_rna <- ragnum_sig_value_rna-1
        }
      }
      output_data$Ragnum_signature_Binary[i] <- ragnum_sig_value_rna
    }
    
    
  }else if(method == "Ranked"){
    # add column to output table
    output_data$Ragnum_signature_Ranked <- NA
    
    rna_rank <- rankGenes(expression_matrix)
    
    rna_ragnum_score <- simpleScore(rna_rank, upSet = ragnum_genes[,annotation])
    
    output_data$Ragnum_signature_Ranked <- rna_ragnum_score$TotalScore
    
    
  }else{
    stop("Method must be one of 'Sum', 'Binary', or 'Ranked'")
    
  }
  
  return(output_data)
}

West_hypoxia_score_within <- function(expression_matrix, method, annotation, gene_list){
  # helper not in statement
  '%!in%' <- function(x,y)!('%in%'(x,y))
  
  # Stop statements
  if("matrix" %!in% class(expression_matrix)){
    stop("Input data must be an expression matrix")
  }
  
  if(ncol(expression_matrix) < 3){## Think about threshold here
    stop("Not enough samples to calculate signature")
  }
  
  if((annotation %in% c("ENSG","SYMBOL","ID_SYMBOL")) == FALSE){
    stop("Annotation must be one of 'ENSG' or 'SYMBOL' or 'ID_SYMBOL'.")
  }
  
  # Set up output table
  output_data <- as.data.frame(colnames(expression_matrix))
  colnames(output_data) <- c("Sample")
  output_data$Genes_used_West <- NA
  
  Signature_Genes_Present <- NA
  Signature_Genes_Absent <- NA
  
  west_genes <- gene_list
  
  # Warnings about missing signature genes
  Signature_Genes_Present <- west_genes[west_genes[,annotation] %in% rownames(expression_matrix), annotation]
  Signature_Genes_Absent <- west_genes[west_genes[,annotation] %!in% rownames(expression_matrix), annotation]
  
  if(length(Signature_Genes_Present) < nrow(west_genes)){
    warning(paste("Warning: ",length(Signature_Genes_Absent),"signature genes not present in input matrix; results may be incomplete or unreliable. Missing genes:",Signature_Genes_Absent, sep = " "))
  }
  
  output_data$Genes_used_West <- toString(Signature_Genes_Present)
  
  # Median centre and scale to unit variance
  expression_matrix_west <- as.matrix(expression_matrix[ Signature_Genes_Present, ])
  expression_matrix_west  <- expression_matrix_west - rowMedians(expression_matrix_west)
  expression_matrix_west_scaled <- (scale((expression_matrix_west), center = F, scale = T))
  
  # Up and Down genes
  up_genes = west_genes[grep("Up", west_genes$Direction), annotation]
  down_genes = west_genes[grep("Down", west_genes$Direction), annotation]
  
  if(method == "West"){
    
  }else if(method == "Sum"){
    
    # add column to output table
    output_data$West_signature_Sum <- NA
    
    # calculate signature
    output_data$West_signature_Sum <- colSums(expression_matrix_west_scaled[up_genes[up_genes %in% Signature_Genes_Present],]) - colSums(expression_matrix_west_scaled[down_genes[down_genes %in% Signature_Genes_Present],])
    
  }else if(method == "Binary"){
    
    # add column to output table
    output_data$West_signature_Binary <- NA
    
    for(i in 1:nrow(output_data)){
      west_sig_value_rna <- 0
      for (j in 1:nrow(expression_matrix_west_scaled)) {
        if((expression_matrix_west_scaled[j,i] > 0) & 
           (rownames(expression_matrix_west_scaled)[j] %in% up_genes)){
          west_sig_value_rna <- west_sig_value_rna + 1
        }else if((expression_matrix_west_scaled[j,i] < 0) & 
                 (rownames(expression_matrix_west_scaled)[j] %in% down_genes)){
          west_sig_value_rna <- west_sig_value_rna + 1
        }else{
          west_sig_value_rna <- west_sig_value_rna-1
        }
      }
      output_data$West_signature_Binary[i] <- west_sig_value_rna
    }
    
  }else if(method == "Ranked"){
    
    # add column to output table
    output_data$West_signature_Ranked <- NA
    
    rna_rank <- rankGenes(expression_matrix)
    
    rna_west_score <- simpleScore(rna_rank, 
                                  upSet = west_genes[grep("Up", west_genes$Direction), annotation],
                                  downSet = west_genes[grep("Down", west_genes$Direction), annotation])
    
    output_data$West_signature_Ranked <- rna_west_score$TotalScore
    
    
  }else{
    stop("Method must be one of 'West', 'Sum', 'Binary', or 'Ranked'. ")
  }
  return(output_data)
}

Buffa_hypoxia_score_within <- function(expression_matrix, method, annotation, gene_list){
  # helper not in statement
  '%!in%' <- function(x,y)!('%in%'(x,y))
  
  # Stop statements
  if("matrix" %!in% class(expression_matrix)){
    stop("Input data must be an expression matrix")
  }
  
  if(ncol(expression_matrix) < 3){## Think about threshold here
    stop("Not enough samples to calculate signature")
  }
  
  if((annotation %in% c("ENSG","SYMBOL","ID_SYMBOL")) == FALSE){
    stop("Annotation must be one of 'ENSG' or 'SYMBOL' or 'ID_SYMBOL'.")
  }
  
  # Set up output table
  output_data <- as.data.frame(colnames(expression_matrix))
  colnames(output_data) <- c("Sample")
  output_data$Genes_used_Buffa <- NA
  
  Signature_Genes_Present <- NA
  Signature_Genes_Absent <- NA
  
  buffa_genes <- gene_list
  
  # Warnings about missing signature genes
  Signature_Genes_Present <- buffa_genes[buffa_genes[,annotation] %in% rownames(expression_matrix), annotation]
  Signature_Genes_Absent <- buffa_genes[buffa_genes[,annotation] %!in% rownames(expression_matrix), annotation]
  
  if(length(Signature_Genes_Present) < nrow(buffa_genes)){
    warning(paste("Warning: ",length(Signature_Genes_Absent),"signature genes not present in input matrix; results may be incomplete or unreliable. Missing genes:",Signature_Genes_Absent, sep = " "))
  }
  
  output_data$Genes_used_Buffa <- toString(Signature_Genes_Present)
  
  # Median centre and scale to unit variance
  expression_matrix_buffa <- as.matrix(expression_matrix[ Signature_Genes_Present, ])
  expression_matrix_buffa  <- expression_matrix_buffa - rowMedians(expression_matrix_buffa)
  expression_matrix_buffa_scaled <- (scale((expression_matrix_buffa), center = F, scale = T))
  
  
  if(method == "Sum"){
    # add column to output table
    output_data$Buffa_signature_Sum <- NA
    
    # calculate signature
    output_data$Buffa_signature_Sum <- colSums(expression_matrix_buffa_scaled)
    
  }else if(method == "Binary"){
    # add column to output table
    output_data$Buffa_signature_Binary <- NA
    
    for(i in 1:nrow(output_data)){
      buffa_sig_value_rna <- 0
      for (j in 1:nrow(expression_matrix_buffa_scaled)) {
        if(expression_matrix_buffa_scaled[j,i] > 0){
          buffa_sig_value_rna <- buffa_sig_value_rna + 1
        }else{
          buffa_sig_value_rna <- buffa_sig_value_rna-1
        }
      }
      output_data$Buffa_signature_Binary[i] <- buffa_sig_value_rna
    }
    
    
  }else if(method == "Ranked"){
    # add column to output table
    output_data$Buffa_signature_Ranked <- NA
    
    rna_rank <- rankGenes(expression_matrix)
    
    rna_buffa_score <- simpleScore(rna_rank, upSet = buffa_genes[,annotation])
    
    output_data$Buffa_signature_Ranked <- rna_buffa_score$TotalScore
    
    
  }else{
    stop("Method must be one of 'Sum', 'Binary', or 'Ranked'")
    
  }
  
  return(output_data)
}

#### Running code ####

ragnum_gene_list_t1 <- ragnum_gene_list

ragnum_gene_list_t1$removed.pseudogenes <- NA
ragnum_gene_list_t1$removed.ig <- NA
ragnum_gene_list_t1$removed.noIds <- NA
ragnum_gene_list_t1$removed.NotDetected <- NA

for(i in 1:nrow(ragnum_gene_list_t1)){
  if(ragnum_gene_list_t1$ENSG[i] %in% rna_dat_supervised@elementMetadata$all.ppcg.trans){
    gene_index <- which(rna_dat_supervised@elementMetadata$all.ppcg.trans == ragnum_gene_list_t1$ENSG[i])
    
    ragnum_gene_list_t1$removed.pseudogenes[i] <- rna_dat_supervised@elementMetadata$removed.pseudogene[gene_index]
    ragnum_gene_list_t1$removed.ig[i] <- rna_dat_supervised@elementMetadata$removed.ig.mt.other.genes[gene_index]
    ragnum_gene_list_t1$removed.noIds[i] <- rna_dat_supervised@elementMetadata$removed.no.gene.ids[gene_index]
    ragnum_gene_list_t1$removed.NotDetected[i] <- rna_dat_supervised@elementMetadata$removed.no.expression[gene_index]
  }
}

# Merge gene list with metadata to get correct row names when calculating signature.
ragnum_gene_list_t1 <- merge(ragnum_gene_list_t1, as.data.frame((rna_dat_supervised@elementMetadata)[which(rna_dat_supervised@elementMetadata$all.ppcg.trans %in% ragnum_gene_list_t1$ENSG),]), 
                             by.x = "ENSG", by.y = "all.ppcg.trans",
                             all.x = TRUE)

# Make a column that matches expression matrix gene name
# ragnum_gene_list$ID_SYMBOL <- paste(ragnum_gene_list$entrezgene_id_BioMart, ragnum_gene_list$hgnc_symbol_BioMart, sep = "||")
ragnum_gene_list_t1$ID_SYMBOL <- ragnum_gene_list_t1$hgnc_symbol

# Check status of West genes in PRPS data - 1 genes removed due to low expression (SLC2A14)
west_gene_list_t1 <- west_gene_list

west_gene_list_t1$removed.pseudogenes <- NA
west_gene_list_t1$removed.ig <- NA
west_gene_list_t1$removed.noIds <- NA
west_gene_list_t1$removed.NotDetected <- NA

for(i in 1:nrow(west_gene_list_t1)){
  if(west_gene_list_t1$ENSG[i] %in% rna_dat_supervised@elementMetadata$all.ppcg.trans){
    gene_index <- which(rna_dat_supervised@elementMetadata$all.ppcg.trans == west_gene_list_t1$ENSG[i])
    
    west_gene_list_t1$removed.pseudogenes[i] <- rna_dat_supervised@elementMetadata$removed.pseudogene[gene_index]
    west_gene_list_t1$removed.ig[i] <- rna_dat_supervised@elementMetadata$removed.ig.mt.other.genes[gene_index]
    west_gene_list_t1$removed.noIds[i] <- rna_dat_supervised@elementMetadata$removed.no.gene.ids[gene_index]
    west_gene_list_t1$removed.NotDetected[i] <- rna_dat_supervised@elementMetadata$removed.no.expression[gene_index]
  }
}

# Merge gene list with metadata to get correct row names when calculating signature.
west_gene_list_t1 <- merge(west_gene_list_t1, as.data.frame((rna_dat_supervised@elementMetadata)[which(rna_dat_supervised@elementMetadata$all.ppcg.trans %in% west_gene_list_t1$ENSG),]), 
                           by.x = "ENSG", by.y = "all.ppcg.trans", all.x = T)

# Make a column that matches expression matrix gene name
# west_gene_list$ID_SYMBOL <- paste(west_gene_list$entrezgene_id_BioMart, west_gene_list$hgnc_symbol_BioMart, sep = "||")
west_gene_list_t1$ID_SYMBOL <- west_gene_list_t1$hgnc_symbol

# Check status of Buffa genes in PRPS data - 1 gene removed due to low expression (CA9)
buffa_gene_list_t1 <- buffa_gene_list

buffa_gene_list_t1$removed.pseudogenes <- NA
buffa_gene_list_t1$removed.ig <- NA
buffa_gene_list_t1$removed.noIds <- NA
buffa_gene_list_t1$removed.NotDetected <- NA

for(i in 1:nrow(buffa_gene_list_t1)){
  if(buffa_gene_list_t1$ENSG[i] %in% rna_dat_supervised@elementMetadata$all.ppcg.trans){
    gene_index <- which(rna_dat_supervised@elementMetadata$all.ppcg.trans == buffa_gene_list_t1$ENSG[i])
    
    buffa_gene_list_t1$removed.pseudogenes[i] <- rna_dat_supervised@elementMetadata$removed.pseudogene[gene_index]
    buffa_gene_list_t1$removed.ig[i] <- rna_dat_supervised@elementMetadata$removed.ig.mt.other.genes[gene_index]
    buffa_gene_list_t1$removed.noIds[i] <- rna_dat_supervised@elementMetadata$removed.no.gene.ids[gene_index]
    buffa_gene_list_t1$removed.NotDetected[i] <- rna_dat_supervised@elementMetadata$removed.no.expression[gene_index]
  }
}

# Merge gene list with metadata to get correct row names when calculating signature.
buffa_gene_list_t1 <- merge(buffa_gene_list_t1, as.data.frame((rna_dat_supervised@elementMetadata)[which(rna_dat_supervised@elementMetadata$all.ppcg.trans %in% buffa_gene_list_t1$ENSG),]), 
                            by.x = "ENSG", by.y = "all.ppcg.trans", all.x = T)

# Make a column that matches expression matrix gene name
# buffa_gene_list$ID_SYMBOL <- paste(buffa_gene_list$entrezgene_id_BioMart, buffa_gene_list$hgnc_symbol_BioMart, sep = "||")
buffa_gene_list_t1$ID_SYMBOL <- buffa_gene_list_t1$hgnc_symbol

# Applying hypoxia functions to data
# Ragnum
all_ragnum_sum_within <- Ragnum_hypoxia_score_within(rna_dat_unsupervised, method = "Sum", annotation = "ID_SYMBOL", gene_list = ragnum_gene_list_t1)
all_ragnum_binary_within <- Ragnum_hypoxia_score_within(rna_dat_unsupervised, method = "Binary", annotation = "ID_SYMBOL", gene_list = ragnum_gene_list_t1)
all_ragnum_ranked_within <- Ragnum_hypoxia_score_within(rna_dat_unsupervised, method = "Ranked", annotation = "ID_SYMBOL", gene_list = subset(ragnum_gene_list_t1, !is.na(ID_SYMBOL)))

# West
all_west_sum_within <- West_hypoxia_score_within(rna_dat_unsupervised, method = "Sum", annotation = "ID_SYMBOL", gene_list = west_gene_list_t1)
all_west_binary_within <- West_hypoxia_score_within(rna_dat_unsupervised, method = "Binary", annotation = "ID_SYMBOL", gene_list = west_gene_list_t1)
all_west_ranked_within <- West_hypoxia_score_within(rna_dat_unsupervised, method = "Ranked", annotation = "ID_SYMBOL", gene_list = subset(west_gene_list_t1, !is.na(ID_SYMBOL)))

# Buffa
all_buffa_sum_within <- Buffa_hypoxia_score_within(rna_dat_unsupervised, method = "Sum", annotation = "ID_SYMBOL", gene_list = buffa_gene_list_t1)
all_buffa_binary_within <- Buffa_hypoxia_score_within(rna_dat_unsupervised, method = "Binary", annotation = "ID_SYMBOL", gene_list = buffa_gene_list_t1)
all_buffa_ranked_within <- Buffa_hypoxia_score_within(rna_dat_unsupervised, method = "Ranked", annotation = "ID_SYMBOL", gene_list = subset(buffa_gene_list_t1, !is.na(ID_SYMBOL)))

# Append label to column name
# Ragnum
colnames(all_ragnum_sum_within)[2:3] <- paste0(colnames(all_ragnum_sum_within)[2:3],"_within")
colnames(all_ragnum_binary_within)[2:3] <- paste0(colnames(all_ragnum_binary_within)[2:3],"_within")
colnames(all_ragnum_ranked_within)[2:3] <- paste0(colnames(all_ragnum_ranked_within)[2:3],"_within")

# West
colnames(all_west_sum_within)[2:3] <- paste0(colnames(all_west_sum_within)[2:3],"_within")
colnames(all_west_binary_within)[2:3] <- paste0(colnames(all_west_binary_within)[2:3],"_within")
colnames(all_west_ranked_within)[2:3] <- paste0(colnames(all_west_ranked_within)[2:3],"_within")

# Buffa
colnames(all_buffa_sum_within)[2:3] <- paste0(colnames(all_buffa_sum_within)[2:3],"_within")
colnames(all_buffa_binary_within)[2:3] <- paste0(colnames(all_buffa_binary_within)[2:3],"_within")
colnames(all_buffa_ranked_within)[2:3] <- paste0(colnames(all_buffa_ranked_within)[2:3],"_within")

# Merge results
# Ragnum
hypoxia_signatures <- all_ragnum_sum_within
hypoxia_signatures <- merge(hypoxia_signatures, all_ragnum_binary_within[,c("Sample","Ragnum_signature_Binary_within")], by = "Sample", all.x = TRUE)
hypoxia_signatures <- merge(hypoxia_signatures, all_ragnum_ranked_within[,c("Sample","Ragnum_signature_Ranked_within")], by = "Sample", all.x = TRUE)
# West
hypoxia_signatures <- merge(hypoxia_signatures, all_west_sum_within, by = "Sample", all.x = TRUE)
hypoxia_signatures <- merge(hypoxia_signatures, all_west_binary_within[,c("Sample","West_signature_Binary_within")], by = "Sample", all.x = TRUE)
hypoxia_signatures <- merge(hypoxia_signatures, all_west_ranked_within[,c("Sample","West_signature_Ranked_within")], by = "Sample", all.x = TRUE)
# Buffa
hypoxia_signatures <- merge(hypoxia_signatures, all_buffa_sum_within, by = "Sample", all.x = TRUE)
hypoxia_signatures <- merge(hypoxia_signatures, all_buffa_binary_within[,c("Sample","Buffa_signature_Binary_within")], by = "Sample", all.x = TRUE)
hypoxia_signatures <- merge(hypoxia_signatures, all_buffa_ranked_within[,c("Sample","Buffa_signature_Ranked_within")], by = "Sample", all.x = TRUE)

hypoxia_signatures <- merge(hypoxia_signatures, rna_dat_supervised@colData, by.x = "Sample", by.y = "Sample.id", all.x = TRUE)
dim(hypoxia_signatures)

hypoxia_signatures <- as.data.frame(hypoxia_signatures)

table(hypoxia_signatures$MetastaticSite.LN.BONE.LIVER.BRAIN.LUNG.BLADDER.ADRENAL.other.)
#### Saving results ####

write.table(hypoxia_signatures, paste0("/mnt/bmh01-rds/UoOxford_David_W/b05055gj/PPCG_analyses/manuscript_v2_outputs/hypoxia_scores/alternative_within_sample_scaled_hypoxia_scores_RNAv1.k10_",Sys.Date(),".txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
