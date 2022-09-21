### Generate count matrix for each tissue

counts <- readRDS("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.rds")
counts <- as.data.frame(counts)

gene.info <- rtracklayer::import("gencode.v26.GRCh38.genes.gtf")
protein.coding.genes <- gene.info[gene.info$type == "gene" & gene.info$gene_type == "protein_coding",]

rownames(counts) <- counts$Name
counts <- counts[protein.coding.genes$gene_id,]
counts$length <- sapply(counts$Name, function(x) sum(rtracklayer::width(gene.info[gene.info$gene_id == x & gene.info$type == "exon"])))


sample_attributes <- data.table::fread("GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")

getTisSamp <- function(tissue) {
  sam <- sample_attributes$SAMPID[sample_attributes$SMTSD == tissue & sample_attributes$SMAFRZE == "RNASEQ"]
  counts_tis <- counts[,colnames(counts) %in% sam]
  colnames(counts_tis) <- paste0("GTEX-", sapply(strsplit(colnames(counts_tis), "-"), "[", 2))
  return(counts_tis)
}


process_file <- function(tissue_name, file_name) {
  Tis_DF <- getTisSamp(tissue_name)
  cov <- data.table::fread(paste0("~/moto/hblab/projects/GTEx/v8/GTEx_Analysis_v8_eQTL_covariates/", tis_save, ".v8.covariates.txt"))
  Tis_DF <- Tis_DF[,colnames(Tis_DF) %in% colnames(cov)]
  Tis_DF <- cbind(counts[,c("Name", "Description", "length")], Tis_DF)
  saveRDS(Tis_DF, paste0("reads/", file_name, ".reads.rds"))
}

tissue_ls <- read.csv("tissue_abbreviations.csv")
for (i in 1:nrow(tissue_ls)) {
  process_file(tissue_name = tissue_ls$tissue_name[i], file_name = tissue_ls$file_name[i])
}


### Downsampling

tis_ls <- list.files("reads/", pattern = "reads")
tis_ls <- gsub(".reads.rds", "", tis_ls)

for (tis in tis_ls) {
  counts_DF <- readRDS(paste0("reads/", tis, ".reads.rds"))
  counts_DF <- counts_DF[,-2]
  colnames(counts_DF)[1] <- "gene"
  samples <- colnames(counts_DF)[!colnames(counts_DF) %in% c("gene", "length")]
  
  sample_sums <- apply(counts_DF[,samples], 2, sum)
  size_min <- min(sample_sums)
  
  downSampling <- function(i) { rmultinom(n = 1, size = size_min, prob = counts_DF[,samples[i]]/sample_sums[i]) }
  counts_down <- as.data.frame(sapply(1:length(samples), downSampling))
  colnames(counts_down) <- samples
  counts_down <- cbind(counts_DF[,1:2], counts_down)
  
  write.csv(counts_down, paste0("GTEx_Sampling/counts/", tis, ".counts.down.csv"), row.names = T, quote = F)
}


