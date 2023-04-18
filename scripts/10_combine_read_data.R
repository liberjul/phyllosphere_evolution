library(rbiom)
library(ggplot2)
library(tidyr)
library(tibble)
library(dplyr)
library(stringr)
# library(DESeq2)
library(MASS)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
### Sample metadata
srr_dat_non <- read.delim("../data/metadata/srr_accessions_with_metadata_nonphyllo.txt",
                      stringsAsFactors = F) %>%
  dplyr::rename(SRR_acc=X)
srr_dat_phy <- read.delim("../data/metadata/phyllo/srr_with_host_taxonomy.txt",
                      stringsAsFactors = F)

length(unique(srr_dat_non$BioProject))
length(unique(srr_dat_phy$BioProject))


import_reads_tax<- function(target, region, layout, srr_dat, x = 0, split=NA, np=F, prop=F, rank=NA, phyllo=FALSE){
  if (phyllo){
    if (is.na(split)){
      biom_path <- paste0("../data/otu_tables/phyllo/feature-table_", target, "_", region, "_", layout, ".biom")
      tax_path <- paste0("../data/constax/phyllo/constax_customdb_taxonomy_", target, "_", region, "_", layout, ".txt")
      fname_path <- paste0("../data/metadata/phyllo/fnames_", target, "_", region, "_", layout, ".txt")
    }
    else {
      biom_path <- paste0("../data/otu_tables/phyllo/feature-table_", target, "_", region, "_", layout, "_", split, ".biom")
      tax_path <- paste0("../data/constax/phyllo/constax_customdb_taxonomy_", target, "_", region, "_", layout, "_", split, ".txt")
      fname_path <- paste0("../data/metadata/phyllo/fnames_", target, "_", region, "_", layout, "_", split, ".txt")
    }
  }
  else {
    if (is.na(split)){
      biom_path <- paste0("../data/otu_tables/feature-table_", target, "_", region, "_", layout, ".biom")
      tax_path <- paste0("../data/constax/constax_customdb_taxonomy_", target, "_", region, "_", layout, ".txt")
      fname_path <- paste0("../data/metadata/fnames_", target, "_", region, "_", layout, ".txt")
    }
    else {
      biom_path <- paste0("../data/otu_tables/feature-table_", target, "_", region, "_", layout, "_", split, ".biom")
      tax_path <- paste0("../data/constax/constax_customdb_taxonomy_", target, "_", region, "_", layout, "_", split, ".txt")
      fname_path <- paste0("../data/metadata/fnames_", target, "_", region, "_", layout, "_", split, ".txt")
    }
  }
  print(biom_path)
  print(tax_path)
  print(fname_path)

  otu_dat <- read.biom(biom_path)$counts %>%
    as.matrix(.) %>%
    as.data.frame(.) %>%
    rownames_to_column(var = "OTU_ID") %>%
    mutate(OTU_name = paste("ASV", row_number(), sep="_"))
  if (np){
    otu_dat %>%
      dplyr::select(starts_with("JL_")) %>%
      colSums(.) -> no_passing_reads
    no_passing_reads[no_passing_reads == 0] %>%
      names(.) -> np
    return(np)
  }
  else {
    tax_dat <- read.delim(tax_path, stringsAsFactors = F)
    comb_dat <- left_join(otu_dat, tax_dat, by = "OTU_ID") %>%
      rename_with(~ paste(.x, region, layout, ifelse(phyllo, "phyllo", "nonphyllo"), sep = "_"),
                  starts_with("JL"))
    fnames_dat <- read.delim(fname_path, skip = 1, header = F, stringsAsFactors = F) %>%
      mutate(sample_name = str_extract(V1, "[JL_|JL]\\d*"),
             SRR_acc = str_extract(V1, ".RR\\d*")) %>%
      left_join(srr_dat, by = "SRR_acc")
    filt_dat <- comb_dat %>%
      filter(! High_level_taxonomy %in% c("Mitochondria", "Chloroplast"),
             ! str_detect(High_level_taxonomy,
                          paste(c("Alveolata", "Archaeplastida", "Archaea"),
                                collapse = "|")),
             # Class %in% c("Microbotryomycetes",
             #              "Agaricostilbomycetes",
             #              "Atractiellomycetes"),
             !is.na(Phylum)) %>%
      mutate(across(starts_with("JL"), ~ifelse(prop, .x / sum(.x), .x)))
    taxa_counts <- data.frame(taxon = as.character(c(NA)))
    # print(colnames(taxa_counts))
    for (rank in c("Class", "Order", "Family", "Genus", "Species")){
      sub_taxa_counts <- filt_dat %>%
        group_by(across(matches(rank))) %>%
        summarize(across(starts_with("JL"), ~sum(.x))) %>%
        rename_with(~ "taxon", matches(rank)) %>%
        # mutate(Rank = rank) %>%
        filter(taxon != "")
      # print(sub_taxa_counts)
      taxa_counts <- full_join(taxa_counts, sub_taxa_counts)
    }
    print(read.biom(biom_path)$counts %>%
            as.matrix(.) %>% sum(.))
    return(taxa_counts)

      # pivot_longer(cols = starts_with("JL"),
      #              names_to = "sample_name", values_to = "Reads") %>%
      # left_join(fnames_dat %>%
      #             dplyr::select(sample_name, SRR_acc),
      #           by = "sample_name")

  }
}

# phyllosphere reads
fungi_its1.its2_wide_pe <- import_reads_tax("Fungi", "ITS1-ITS2", "pe", srr_dat_phy, prop=T, phyllo=T)
fungi_its1_wide_pe <- import_reads_tax("Fungi", "ITS1", "pe", srr_dat_phy, prop=T, phyllo=T)
fungi_its2_wide_pe <- import_reads_tax("Fungi", "ITS2", "pe", srr_dat_phy, prop=T, phyllo=T)
fungi_its1.its2_wide_se <- import_reads_tax("Fungi", "ITS1-ITS2", "se", srr_dat_phy, prop=T, phyllo=T)
fungi_its1_wide_se <- import_reads_tax("Fungi", "ITS1", "se", srr_dat_phy, prop=T, phyllo=T)
fungi_its2_wide_se <- import_reads_tax("Fungi", "ITS2", "se", srr_dat_phy, prop=T, phyllo=T)

phy_fungi_wide <- full_join(fungi_its1.its2_wide_pe, fungi_its1_wide_pe) %>%
  full_join(fungi_its2_wide_pe) %>%
  full_join(fungi_its1.its2_wide_se) %>%
  full_join(fungi_its1_wide_se) %>%
  group_by(taxon) %>%
  summarise(across(starts_with("JL"), ~ sum(.x, na.rm = T)))

phy_fungi_wide
samcount_phy <- phy_fungi_wide %>%
  colnames(.) %>%
  str_detect("JL") %>%
  sum(.)
samcount_phy

phy_fungi_wide_norm <- phy_fungi_wide %>%
  mutate(across(starts_with("JL"), ~ .x/samcount_phy))

# non-phyllosphere reads
fungi_its1.its2_wide_paired <- import_reads_tax("Fungi", "ITS1-ITS2", "paired", srr_dat_non, prop=T, phyllo=F)
fungi_its1_wide_paired <- import_reads_tax("Fungi", "ITS1", "paired", srr_dat_non, prop=T, phyllo=F)
fungi_its2_wide_paired <- import_reads_tax("Fungi", "ITS2", "paired", srr_dat_non, prop=T, phyllo=F)
fungi_its2_wide_single <- import_reads_tax("Fungi", "ITS2", "single", srr_dat_non, prop=T, phyllo=F)

non_fungi_wide <- full_join(fungi_its1.its2_wide_paired, fungi_its1_wide_paired) %>%
  full_join(fungi_its2_wide_paired) %>%
  full_join(fungi_its2_wide_single) %>%
  group_by(taxon) %>%
  summarise(across(starts_with("JL"), ~ sum(.x, na.rm = T)))
colnames(non_fungi_wide)
samcount_non <- non_fungi_wide %>%
  colnames(.) %>%
  str_detect("JL") %>%
  sum(.)
samcount_non
non_fungi_wide_norm <- non_fungi_wide %>%
  mutate(across(starts_with("JL"), ~ .x/samcount_non))

cts <- inner_join(phy_fungi_wide_norm, non_fungi_wide_norm, by = "taxon") %>%
  filter(!is.na(taxon)) %>%
  column_to_rownames(var = "taxon") %>%
  as.matrix()
cts
coldata <- data.frame(cnames = colnames(cts)) %>%
  filter(cnames != "taxon") %>%
  mutate(niche = str_extract(cnames, "phyllo|nonphyllo")) %>%
  column_to_rownames(var="cnames")
coldata %>%
  group_by(niche) %>%
  summarize(sample_count = n())


# dds <- DESeqDataSetFromMatrix(countData = cts,
#                               colData = coldata,
#                               design = ~ niche)
# 
# dds <- estimateSizeFactors(dds, type="iterate")
# dds <- estimateDispersionsGeneEst(dds)

# deseq_rowwise <- function(cts, coldata, i){
#   sub_cts <- cts[i,cts[i,] > 0]
#   sub_cts <- t(as.matrix(sub_cts))
#   rownames(sub_cts) <- rownames(cts)[i]
#   
#   
#   sub_coldat <- data.frame(niche = coldata[colnames(sub_cts), "niche"],
#                            row.names = colnames(sub_cts))
#   dds <- DESeqDataSetFromMatrix(countData = sub_cts,
#                                 colData = sub_coldat,
#                                 design = ~ niche)
#   dds <- estimateSizeFactors(dds)
#   dds <- estimateDispersionsGeneEst(dds)
#   dispersions(dds) <- mcols(dds)$dispGeneEst
#   dds <- nbinomWaldTest(dds)
#   res <- results(dds)
#   return(as.data.frame(res))
# }
log2fold_rowwise <- function(cts, coldata, i){
  sub_cts <- cts[i,cts[i,] > 0]
  sub_cts <- t(as.matrix(sub_cts))
  rownames(sub_cts) <- rownames(cts)[i]
  
  sub_coldat <- data.frame(niche = coldata[colnames(sub_cts), "niche"],
                           row.names = colnames(sub_cts))

  phyllo_mean <- mean(sub_cts[1, rownames(sub_coldat)[sub_coldat$niche == "phyllo"]])
  nonphyllo_mean <- mean(sub_cts[1, rownames(sub_coldat)[sub_coldat$niche == "nonphyllo"]])
  nb_test_df <- data.frame(niche = coldata[colnames(sub_cts), "niche"],
                           value = sub_cts[1,],
                           row.names = colnames(sub_cts))
  m <- glm.nb(value ~ niche, data = nb_test_df)
  coefs_df <- summary(m)$coefficients
  
  log2fold = log2(phyllo_mean/nonphyllo_mean)
  res <- data.frame(taxon = c(rownames(cts)[i]),
                    p_mean = c(phyllo_mean),
                    n_mean = c(nonphyllo_mean),
                    log2fold = c(log2fold),
                    glm.nb_est = coefs_df[2,"Estimate"],
                    glm.nb_stderr = coefs_df[2,"Std. Error"],
                    glm.nb_zval = coefs_df[2,"z value"],
                    glm.nb_pval = coefs_df[2,"Pr(>|z|)"]
                    )
  return(res)
}

# deseq_res <- deseq_rowwise(cts, coldata, 1)
# deseq_res
res <- log2fold_rowwise(cts, coldata, 1)
for (i in 2:dim(cts)[1]){
  res <- rbind(res, log2fold_rowwise(cts, coldata, i))
}
res
res %>%
  ggplot(aes(x = log2fold)) +
  geom_histogram()-> g
g
ggsave("../figures/log2fold_distribution.svg", g, height = 5, width = 7)


cts_df <- data.frame(x = as.vector(cts[cts != 0]))
cts_df %>%
  ggplot(aes(x=x))+
  geom_histogram() +
  scale_x_log10() +
  annotate(geom="text", x = 4e-6, y = 4000,
         label=paste0("Mean = ", round(mean(cts_df$x), 5),
                      "\nVariance = ", round(var(cts_df$x), 5))) +
  labs(x = "Normalized read proportion",
       y = "Count") -> g
ggsave("../figures/nonzero_cts_distribution.svg", g, height = 5, width = 7)

res %>%
  ggplot(aes(x = glm.nb_est)) +
  geom_histogram()
res %>%
  ggplot(aes(x = glm.nb_pval)) +
  geom_histogram()
res %>%
  arrange(-1*glm.nb_zval)

res %>%
  arrange(glm.nb_zval)
res %>%
  filter(str_detect(taxon, "mycetes")) %>%
  mutate(taxon = factor(taxon, levels = taxon[order(log2fold)])) %>%
  # arrange(log2fold) %>%
  ggplot(aes(x = taxon, y = log2fold)) +
  geom_col() + 
  theme(axis.text.x.bottom = element_text(angle = -45,
                                          hjust = 0),
        plot.margin = margin(1,1.5,0.5,0.5, "cm")) + 
  labs(x = "Class", y = "Log2 Fold Enrichment in Phyllosphere")

write.csv(res, "../data/stats/glmnb_results_all_fungi.csv")
