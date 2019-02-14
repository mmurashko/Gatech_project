rm(list=ls())
library("tximport")
library("readr")
library("tximportData")
library("DESeq2")
library(ggplot2)
library("lettercase")
library("data.table")
library(grid)
library(gridExtra)
library(biomaRt)
library("mygene")

#args <- c("pathway", "/Users/Matvey/Desktop/GeorgiaTech_project/output_data/6hrs", "samples_6hrs.csv", "_3D_6hrs", "1", "GO:0007219")
#args <- c("set", "/Users/Matvey/Desktop/GeorgiaTech_project/output_data/6hrs", "samples_6hrs.csv", "_3D_6hrs", "0.5", "Users/Matvey/Desktop/GeorgiaTech_project/Set_of_wanted_genes.txt")


args <- commandArgs(TRUE)
#Choose a mode
mode <- as.character(args[1]) #"pathway" or "set"
#Main directory of salmon data
main_directory <- as.character(args[2])
#Names of sample files in one csv file(.csv file must be in the same directory as salmon files)
samples_file <- as.character(args[3])
#Add label to every new file
label <- as.character(args[4])
threshold <- args[5]


if (mode == "pathway"){
  dir.create(file.path(main_directory, paste0("threshold_", gsub('[.]', '_', threshold))))
  directory <- file.path(main_directory, paste0("threshold_", gsub('[.]', '_', threshold)))
} else if (mode == "set"){
  dir.create(file.path(main_directory, "wanted_genes"))
  directory <- file.path(main_directory, "wanted_genes")
}



#directory of reference between genes and transcripts
dir <- system.file("extdata", package="tximportData")


samples <- read.csv2(file.path(main_directory, samples_file), header=TRUE)
rownames(samples) <- samples$name
samples[,]
files <- file.path(main_directory, samples$name, "salmon", "quant.sf")
names(files) <- samples$name


#File that contains reference between genes and transcripts 
tx2gene <- read_csv(file.path(dir, "tx2gene.gencode.v27.csv"))


#Counts and abundance(!) of genes, converts a salmon data for transcripts to a data for genes
txi <- tximport(files, type="salmon", tx2gene=tx2gene)


#Creating a matrix of abundances
abundance <- txi$abundance
abundance <- as.data.frame(abundance)
colnames(abundance) <- c("NT", "CCM1", "CCM2", "CCM3")


#Genes that are not expressed in one of the cases
not_expressed <- abundance[abundance[,1] == 0 & abundance[,2] != 0 | abundance[,1] == 0 & abundance[,3] != 0 |abundance[,1] == 0 & abundance[,4] != 0 |
                          abundance[,1] != 0 & abundance[,2] == 0 | abundance[,1] != 0 & abundance[,3] == 0 |abundance[,1] != 0 & abundance[,4] == 0 , ]


#counting LogFoldChange
abundance <- within(abundance, logFoldChangeCCM1 <- log2(`CCM1` / `NT`)) #log2(`S11_3D_CCM1_2hrs-48185258` / `S11_3D_NT_2hrs-48176275`))
abundance <- within(abundance, logFoldChangeCCM2 <- log2(`CCM2` / `NT`))
abundance <- within(abundance, logFoldChangeCCM3 <- log2(`CCM3` / `NT`))
rownames(abundance) <- gsub("\\..*","",rownames(abundance))
gene_symbols <- queryMany(rownames(abundance), scopes='ensembl.gene', fields=c('symbol'), species='human')
abundance$gene_name <- as.character(lapply(rownames(abundance), function(x){do.call("paste", c(as.list(gene_symbols[gene_symbols$query==x, "symbol"]), sep = ", "))}))
write.csv2(abundance, file = file.path(directory,paste("abundances", label, ".csv", sep="")))


#genes that are always expressed
abundance_no_zeros <- abundance[abundance[,1] != 0 & abundance[,2] != 0 & abundance[,3] != 0 & abundance[,4] != 0, ] #matrix of abundance


#graphing LogFoldChange for all genes
min_CCM1 <- -ceiling(abs(min(abundance_no_zeros$logFoldChangeCCM1)))
max_CCM1 <- ceiling(max(abundance_no_zeros$logFoldChangeCCM1))
min_CCM2 <- -ceiling(abs(min(abundance_no_zeros$logFoldChangeCCM2)))
max_CCM2 <- ceiling(max(abundance_no_zeros$logFoldChangeCCM2))
min_CCM3 <- -ceiling(abs(min(abundance_no_zeros$logFoldChangeCCM3)))
max_CCM3 <- ceiling(max(abundance_no_zeros$logFoldChangeCCM3))

own_graph_CCM1 = ggplot(data = abundance_no_zeros, aes(x = logFoldChangeCCM1)) + geom_histogram(binwidth = 0.1) + scale_x_continuous(breaks = seq(min_CCM1,max_CCM1), limits=c(min_CCM1,max_CCM1))
#own_graph_CCM1
ggsave(file.path(directory, paste("own_graph_CCM1", label, ".png", sep = "")))

own_graph_CCM2 = ggplot(data = abundance_no_zeros, aes(x = logFoldChangeCCM2)) + geom_histogram(binwidth = 0.1) + scale_x_continuous(breaks = seq(min_CCM2,max_CCM2), limits=c(min_CCM2,max_CCM2))
#own_graph_CCM2
ggsave(file.path(directory,paste("own_graph_CCM2", label, ".png", sep = "")))

own_graph_CCM3 = ggplot(data = abundance_no_zeros, aes(x = logFoldChangeCCM3)) + geom_histogram(binwidth = 0.1) + scale_x_continuous(breaks = seq(min_CCM3,max_CCM3), limits=c(min_CCM3,max_CCM3))
#own_graph_CCM3
ggsave(file.path(directory,paste("own_graph_CCM3", label, ".png", sep = "")))


#differentially expressed genes
differential_gene_expression_CCM1 <- abundance_no_zeros[abs(abundance_no_zeros$logFoldChangeCCM1) > threshold, ]
differential_gene_expression_CCM2 <- abundance_no_zeros[abs(abundance_no_zeros$logFoldChangeCCM2) > threshold, ]
differential_gene_expression_CCM3 <- abundance_no_zeros[abs(abundance_no_zeros$logFoldChangeCCM3) > threshold, ]

differential_gene_expression <- abundance_no_zeros[abs(abundance_no_zeros$logFoldChangeCCM1) > threshold|abs(abundance_no_zeros$logFoldChangeCCM2) > threshold|abs(abundance_no_zeros$logFoldChangeCCM3) > threshold, ]

write.csv2(differential_gene_expression, file = file.path(directory,paste("differentialy_expressed_genes", label, "_", gsub('[.]', '_', threshold), ".csv",sep = "")))
write.csv2(differential_gene_expression_CCM1, file = file.path(directory,paste("differentialy_expressed_genes_CCM1", label, "_", gsub('[.]', '_', threshold), ".csv",sep = "")))
write.csv2(differential_gene_expression_CCM2, file = file.path(directory,paste("differentialy_expressed_genes_CCM2", label, "_", gsub('[.]', '_', threshold), ".csv",sep = "")))
write.csv2(differential_gene_expression_CCM3, file = file.path(directory,paste("differentialy_expressed_genes_CCM3", label, "_", gsub('[.]', '_', threshold), ".csv",sep = "")))


#Removing the number of update(.*) from initial data
rownames(not_expressed) <- gsub("\\..*","",rownames(not_expressed))
rownames(abundance_no_zeros) <- gsub("\\..*","",rownames(abundance_no_zeros))


#Genes with no expression in one of the cases saved to .csv file
not_expressed$gene_id <- rownames(not_expressed)
write.csv2(not_expressed, file = file.path(directory,paste("no_expression", label, ".csv", sep="")))


if (mode == "pathway") {
  
}
  #Pathway from Gene Ontology
  GO_id <- args[6]
  
  
  #Matching ensembl id to each gene symbol using gene_names file
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations
  
  gene.data <- getBM(attributes=c('hgnc_symbol', 'go_id'),
                     filters = 'go', values = GO_id, mart = ensembl)
  
  pathway_genes = unique(gene.data$hgnc_symbol)
  
  info_pathway = queryMany(pathway_genes, scopes='symbol', fields=c('all'), species='human')
  pathway_ensembl = unlist(lapply(1:length(info_pathway$ensembl), function(i) info_pathway$ensembl[[i]]$gene))
  
  
  #Genes from the pathway with no expression in one of the cases saved to .csv file
  not_expressed_names = as.character(lapply(pathway_ensembl, function(x){not_expressed[rownames(not_expressed) == x, 5]}))
  genes_with_no_expression <- na.omit(not_expressed[not_expressed_names, ])
  genes_with_no_expression$gene_name <- as.character(lapply(rownames(genes_with_no_expression), function(x){do.call("paste", c(as.list(gene_symbols[gene_symbols$query==x, "symbol"]), sep = ", "))}))
  write.csv2(genes_with_no_expression, file = file.path(directory,paste("no_expression_from_the_pathway", label, ".csv", sep="")))


  #Creating gene_expression matrix and drawing plots
  gene_expression_CCM1 = lapply(pathway_ensembl, function(x){abundance_no_zeros[rownames(abundance_no_zeros) == x, 5]})
  gene_expression_CCM2 = lapply(pathway_ensembl, function(x){abundance_no_zeros[rownames(abundance_no_zeros) == x, 6]})
  gene_expression_CCM3 = lapply(pathway_ensembl, function(x){abundance_no_zeros[rownames(abundance_no_zeros) == x, 7]})

  
  
  gene_expression_matrix <- data.frame(matrix(nrow = length(pathway_ensembl), ncol = 0))
  gene_expression_matrix$expression_CCM1 <- as.numeric(gene_expression_CCM1)
  gene_expression_matrix$expression_CCM2 <- as.numeric(gene_expression_CCM2)
  gene_expression_matrix$expression_CCM3 <- as.numeric(gene_expression_CCM3)
  gene_expression_matrix$gene_id <- as.character(pathway_ensembl)
  gene_expression_matrix_for_graph_CCM1 <- gene_expression_matrix[lapply(gene_expression_CCM1, function(x){length(x)}) == TRUE, ]
  gene_expression_matrix_for_graph_CCM2 <- gene_expression_matrix[lapply(gene_expression_CCM2, function(x){length(x)}) == TRUE, ]
  gene_expression_matrix_for_graph_CCM3 <- gene_expression_matrix[lapply(gene_expression_CCM3, function(x){length(x)}) == TRUE, ]

  graph_pathway_CCM1 = ggplot(data = gene_expression_matrix_for_graph_CCM1, aes(x = expression_CCM1)) + geom_histogram(binwidth = 0.1) + scale_x_continuous(limits=c(min(gene_expression_matrix_for_graph_CCM1$expression_CCM1)*1.1, max(gene_expression_matrix_for_graph_CCM1$expression_CCM1*1.1)))
  ggsave(file.path(directory,paste("pathway_graph_CCM1", label, ".png", sep = "")))

  graph_pathway_CCM2 = ggplot(data = gene_expression_matrix_for_graph_CCM2, aes(x = expression_CCM2)) + geom_histogram(binwidth = 0.1) + scale_x_continuous(limits=c(min(gene_expression_matrix_for_graph_CCM2$expression_CCM2)*1.1, max(gene_expression_matrix_for_graph_CCM2$expression_CCM2*1.1)))
  ggsave(file.path(directory,paste("pathway_graph_CCM2", label, ".png", sep = "")))

  graph_pathway_CCM3 = ggplot(data = gene_expression_matrix_for_graph_CCM3, aes(x = expression_CCM3)) + geom_histogram(binwidth = 0.1) + scale_x_continuous(limits=c(min(gene_expression_matrix_for_graph_CCM3$expression_CCM3)*1.1, max(gene_expression_matrix_for_graph_CCM3$expression_CCM3*1.1)))
  ggsave(file.path(directory,paste("pathway_graph_CCM3", label, ".png", sep = "")))


  #Genes with LogFoldChange > 1
  differential_expression_in_pathway_CCM1 <- gene_expression_matrix_for_graph_CCM1[abs(gene_expression_matrix_for_graph_CCM1$expression_CCM1) > threshold, ]
  write.csv2(differential_expression_in_pathway_CCM1, file = file.path(directory,paste("differential_expression_in_pathway_CCM1", label, "_", gsub('[.]', '_', threshold), ".csv",sep = "")))

  differential_expression_in_pathway_CCM2 <- gene_expression_matrix_for_graph_CCM2[abs(gene_expression_matrix_for_graph_CCM2$expression_CCM2) > threshold, ]
  write.csv2(differential_expression_in_pathway_CCM2, file = file.path(directory,paste("differential_expression_in_pathway_CCM2", label, "_", gsub('[.]', '_', threshold), ".csv",sep = "")))

  differential_expression_in_pathway_CCM3 <- gene_expression_matrix_for_graph_CCM3[abs(gene_expression_matrix_for_graph_CCM3$expression_CCM3) > threshold, ]
  write.csv2(differential_expression_in_pathway_CCM3, file = file.path(directory,paste("differential_expression_in_pathway_CCM3", label,"_", gsub('[.]', '_', threshold), ".csv",sep = "")))


  #Creating a result table with only differentially expressed genes from a pathway
  result_table <- rbind(differential_expression_in_pathway_CCM1, differential_expression_in_pathway_CCM2, differential_expression_in_pathway_CCM3)
  result_table <- setkey(result_table, NULL)
  result_table <- unique(result_table)

  result_table$gene_name <- as.character(lapply(result_table$gene_id, function(x){do.call("paste", c(as.list(gene_symbols[gene_symbols$query==x, "symbol"]), sep = ", "))}))
  colnames(result_table) <- c("CCM1", "CCM2", "CCM3", "gene_id", "gene_name")
  write.csv2(result_table, file = file.path(directory,paste("results_of_differential_expression_in_pathway", label, "_", gsub('[.]', '_', threshold), ".csv",sep = "")))

  
  #Drawing plots
  n <- 1
  for (gene in row.names(result_table)) {
    if (n %% 4 == 1){
      g1 <- ggplot(data = melt(result_table[gene, 1:3]), aes(x = variable, y = value)) + geom_bar(stat="identity") + labs(x = result_table[gene,5], y = "")
      n <- n + 1
    } else if (n %% 4 == 2){
      g2 <- ggplot(data = melt(result_table[gene, 1:3]), aes(x = variable, y = value)) + geom_bar(stat="identity") + labs(x = result_table[gene,5], y = "")
      n <- n + 1
    } else if (n %% 4 == 3){
      g3 <- ggplot(data = melt(result_table[gene, 1:3]), aes(x = variable, y = value)) + geom_bar(stat="identity") + labs(x = result_table[gene,5], y = "")
      n <- n + 1
    } else if (n %% 4 == 0){
      g4 <- ggplot(data = melt(result_table[gene, 1:3]), aes(x = variable, y = value)) + geom_bar(stat="identity") + labs(x = result_table[gene,5], y = "")
      g <- grid.arrange(g1, g2, g3, g4, nrow = 2)
      ggsave(file.path(directory, paste("graphs", n-3 , "-", n, label, "_", gsub('[.]', '_', threshold), ".png",sep="")), g)
      n <- n + 1
    }
  }

  if ((n-1) %% 4 == 1){
    lay <- rbind(c(1,NA),
                 c(NA,NA))
    g <- grid.arrange(g1, layout_matrix = lay)
    ggsave(file.path(directory, paste("graphs", n-1, label, "_", gsub('[.]', '_', threshold),".png",sep="")), g)

  } else if ((n-1) %% 4 == 2){
    lay <- rbind(c(1,2),
                 c(NA,NA))
    g <- grid.arrange(g1, g2, layout_matrix = lay)
    ggsave(file.path(directory, paste("graphs", n-2, "-", n-1, label, "_", gsub('[.]', '_', threshold),".png",sep="")), g)

  } else if ((n-1) %% 4 == 3){
    lay <- rbind(c(1,2),
                 c(3,NA))
    g <- grid.arrange(g1, g2, g3, layout_matrix = lay)
    ggsave(file.path(directory, paste("graphs", n-3 , "-", n-1, label, "_", gsub('[.]', '_', threshold),".png",sep="")), g)

  }

  #Other way to draw plots
  plots <- list()
  x <- 1
  while (x*(x+1) < length(row.names(result_table))) {
    print(x)
    x = x+1
  }
  n <- x*(x+1)
  #n
  lay <-  matrix(1:n,nrow=x, ncol=x+1, byrow = TRUE)
  #lay
  i <- 1
  for (gene in row.names(result_table)) {
    g <- ggplot(data = melt(result_table[gene, 1:3]), aes(x = variable, y = value)) + geom_bar(stat="identity") + labs(x = result_table[gene,5], y = "")  + theme(axis.text.x = element_text(angle = 45, colour="grey20",size=5),axis.text.y = element_text(colour="grey20",size=10))
    plots[[i]]<-g
    i <- i+1
  }


  ml <- marrangeGrob(plots, layout_matrix = lay)
  ml
  ggsave(file.path(directory, paste("graphs", 1, "-", length(row.names(result_table)), label, "_", gsub('[.]', '_', threshold), ".png",sep="")), ml)
  
  
  
  #Enrichment using Fisher's exact test
  print("Fisher's exact test with adjusted set of genes")
  IDs <- abundance$NT != 0 | abundance$CCM1 != 0 | abundance$CCM2 != 0 | abundance$CCM3 != 0
  gene_universe = rownames(abundance[IDs,])
  adjusted_pathway <- pathway_ensembl[pathway_ensembl %in% gene_universe]
  
  #CCM1
  a = length(differential_expression_in_pathway_CCM1[[4]] %in% adjusted_pathway)
  b = length(adjusted_pathway) - a
  c = length(rownames(differential_gene_expression_CCM1)) - a
  d = length(gene_universe) - a - b - c
  m1 = matrix(c(a,b,c,d), byrow = T, 2, 2)
  res1adj = fisher.test(m1, alternative = "greater")
  print(res1adj$p.value)
  
  #CCM2
  a = length(differential_expression_in_pathway_CCM2[[4]] %in% adjusted_pathway)
  b = length(adjusted_pathway) - a
  c = length(rownames(differential_gene_expression_CCM2)) - a
  d = length(gene_universe) - a - b - c
  m2 = matrix(c(a,b,c,d), byrow = T, 2, 2)
  res2adj = fisher.test(m2, alternative = "greater")
  print(res2adj$p.value)
  
  #CCM3
  a = length(differential_expression_in_pathway_CCM3[[4]] %in% adjusted_pathway)
  b = length(adjusted_pathway) - a
  c = length(rownames(differential_gene_expression_CCM3)) - a
  d = length(gene_universe) - a - b - c
  m3 = matrix(c(a,b,c,d), byrow = T, 2, 2)
  res3adj = fisher.test(m3, alternative = "greater")
  print(res3adj$p.value)
  

  enrichment_matrix = matrix(ncol = 0, nrow = 3)
  enrichment_matrix$adjusted_background <- c(res1adj$p.value, res2adj$p.value, res3adj$p.value)
  enrichment <- as.data.frame(enrichment_matrix)
  enrichment
  m1
  m2
  m3
  write.csv2(enrichment, file = file.path(directory,paste("enrichment", label, "_", gsub('[.]', '_', threshold), ".csv",sep = "")))
  
}else if (mode == "set"){
  wanted_genes_file <- args[6]
  wanted_genes <- read.table(wanted_genes_file)
  wanted_genes <- as.character(wanted_genes[,1])

  info_wanted = queryMany(wanted_genes, scopes='symbol', fields=c('all'), species='human')
  wanted_genes_IDs = unlist(lapply(1:length(info_wanted$ensembl), function(i) info_wanted$ensembl[[i]]$gene))
  
  #Drawing plots for wanted genes
  table_of_abundance <- na.exclude(abundance)
  rownames(table_of_abundance) <- gsub("\\..*","",rownames(table_of_abundance))

  wanted_gene_CCM1 = lapply(wanted_genes_IDs, function(x){table_of_abundance[rownames(table_of_abundance) == x, 5]})
  wanted_gene_CCM2 = lapply(wanted_genes_IDs, function(x){table_of_abundance[rownames(table_of_abundance) == x, 6]})
  wanted_gene_CCM3 = lapply(wanted_genes_IDs, function(x){table_of_abundance[rownames(table_of_abundance) == x, 7]})

  wanted_genes_matrix <- data.frame(matrix(nrow = length(wanted_genes_IDs), ncol = 0))
  wanted_genes_matrix$CCM1 <- as.numeric(wanted_gene_CCM1)
  wanted_genes_matrix$CCM2 <- as.numeric(wanted_gene_CCM2)
  wanted_genes_matrix$CCM3 <- as.numeric(wanted_gene_CCM3)
  wanted_genes_matrix$gene_id <- wanted_genes_IDs
  rownames(wanted_genes_matrix) <- wanted_genes_IDs
  wanted_genes_matrix$symbol <- lapply(wanted_genes_matrix$gene_id, function(x){do.call("paste", c(as.list(gene_symbols[gene_symbols$query==x, "symbol"]), sep = ", "))})
  wanted_genes_matrix <- na.omit(wanted_genes_matrix)
  wanted_genes_matrix
  

  wanted_plots <- list()
  x <- 1
  while (x*(x+1) < length(row.names(wanted_genes_matrix))) {
    x = x+1
  }
  n <- x*(x+1)
  #n
  lay <-  matrix(1:n,nrow=x, ncol=x+1, byrow = TRUE)
  #lay
  i <- 1
  for (gene in row.names(wanted_genes_matrix)) {
    single_matrix <- melt(wanted_genes_matrix[gene, 1:3])
    single_matrix$compare <- ifelse(abs(single_matrix$value) < 0.3, "doubtful", "improtant")
    g <- ggplot(data = single_matrix, aes(x = variable, y = value, fill = compare)) + geom_bar(stat="identity") + labs(x = wanted_genes_matrix[gene,5][[1]], y = "")  + theme(axis.text.x = element_text(angle = 45, colour="grey20",size=5),axis.text.y = element_text(colour="grey20",size=10),legend.position="none") + scale_fill_manual(name="Importance", labels = c("doubtful", "improtant"), values = c("doubtful"="grey", "improtant"="green"))
    wanted_plots[[i]]<-g
    i <- i+1
  }


  wanted <- marrangeGrob(wanted_plots, layout_matrix = lay)
  wanted
  ggsave(file.path(directory, paste("wanted_genes", label, ".png",sep="")), wanted)
}
















