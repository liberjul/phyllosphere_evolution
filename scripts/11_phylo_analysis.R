library(phytools)
library(ggtree)
library(treeio)
library(tidytree)
library(ggplot2)
library(tidyr)
library(tibble)
library(dplyr)
library(phyr) # for cor_phylo
library(phylosignal)
library(phylobase)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

tree <- read.beast("../data/trees/20230111_2_mlst_aimania_iqtree_formatted_labeled.nex")

tree@data <- tree@data %>%
  rename(name = `!name`)

tree

enrich_data <- read.csv("../data/stats/glmnb_results_all_fungi.csv")
enrich_data %>%
  filter(taxon %in% tree@data$name) %>%
  mutate(name = taxon)-> filt_enrich

filt_enrich


tree@data <- left_join(tree@data, filt_enrich, by = "name")
tree@data %>%
  filter(name %in% filt_enrich$taxon) %>%
  mutate(clade_label = taxon)-> trait_data_enrich_only
trait_data_enrich_only

g <- ggtree(tree, aes(color = log2fold)) +
  geom_nodepoint(aes(size = 1*!is.na(log2fold))) +
  geom_tiplab(aes(label = str_replace_all(label.x, "_", " ")), size = 1.6,
              fontface = "italic") +
  xlim(0, .8) + 
  scale_size_continuous(limits=c(0,5), breaks = c(0, 1)) +
  theme(legend.position = c(0.11, 0.5)) +
  guides(size = "none") +
  labs(color = "Log2 Fold\nEnrichment\nin Phyllosphere")
g
ggsave("../figures/tree_with_trait_only_data_nodes.svg", g, height = 10, width = 6)





phylotree <- as.phylo(tree)
tree@data %>%
  filter(!is.na(log2fold)) %>%
  rowwise() %>%
  mutate(daughter_tip_ct = length(getDescendants(phylotree, node))) %>%
  arrange(desc(daughter_tip_ct)) -> enrich_data_tip_ct
enrich_data_tip_ct

trait_data <- tree %>% as_tibble() %>%
  rename(label = label.x)
trait_data$log2fold <- NA
trait_data$glm.nb_stderr <- NA
trait_data$inferred <- TRUE

for (i in 1:dim(enrich_data_tip_ct)[1]){
  enrich_data_tip_ct$node[i]
  trait_data$log2fold[trait_data$node == i] <- enrich_data_tip_ct$log2fold[i]
  trait_data$glm.nb_stderr[trait_data$node == i] <- enrich_data_tip_ct$glm.nb_stderr[i]
  trait_data$inferred[trait_data$node == enrich_data_tip_ct$node[i]] <- FALSE
  
  descendents <- offspring(phylotree, enrich_data_tip_ct$node[i])
  trait_data$log2fold[trait_data$node %in% descendents] <- enrich_data_tip_ct$log2fold[i]
  trait_data$glm.nb_stderr[trait_data$node %in% descendents] <- enrich_data_tip_ct$glm.nb_stderr[i]
  trait_data$inferred[trait_data$node %in% descendents] <- TRUE
}
trait_data

tree2 <- tree %>%
  left_join(trait_data, by = "node")
# tree2 <- trait_data %>% as.treedata()
# tree2 %>% as_tibble()
phylotree2 <- as.phylo(tree2)

trait_data
ggtree(tree2, aes(color = log2fold)) +
  geom_nodepoint() +
  geom_tiplab(aes(label = str_replace_all(label.x, "_", " ")), size = 1.6,
              fontface = "italic") +
  xlim(0, .8) + 
  scale_size_continuous(limits=c(0,5), breaks = c(0, 1)) +
  theme(legend.position = c(0.11, 0.5)) +
  guides(size = "none") +
  labs(color = "Log2 Fold\nEnrichment\nin Phyllosphere") -> g
ggsave("../figures/tree_with_trait_all_nodes.svg", g, height = 10, width = 6)

phylotree <- as.phylo(tree)
phylotree$tip.label

trait_data$num_desc <- NA
for (i in 1:dim(trait_data)[1]){
  trait_data$num_desc[i] <- length(offspring(tree, trait_data$node[i]))
}
trait_data
trait_data %>% filter(num_desc == 0) -> tips

trait_data_noinf <- trait_data
trait_data_noinf$log2fold[trait_data$inferred] <- NA
trait_data_noinf$glm.nb_stderr[trait_data$inferred] <- NA

p4d_tree <- as(phylotree, "phylo4d")
tdata(p4d_tree) <- trait_data_noinf %>%
  dplyr::select(log2fold:glm.nb_stderr)
tdata(p4d_tree) %>% view()

# distance-trait correlation

node_dist_mat <- dist.nodes(phylotree)
node_dist_mat
trait_data_enrich_only

trait_diff_mat <- matrix(nrow=dim(trait_data_enrich_only)[1],
                         ncol=dim(trait_data_enrich_only)[1])
colnames(trait_diff_mat) <- trait_data_enrich_only$node
rownames(trait_diff_mat) <- trait_data_enrich_only$node
trait_err_mat <- trait_diff_mat

# trait_data_enrich_only
# filt_enrich

corr_data_long <- data.frame(node.x = c(), node.y = c(),
                             diff = c(), err = c(), phy_dist = c())
for (i in 1:dim(trait_data_enrich_only)[1]){
  for (j in 1:dim(trait_data_enrich_only)[1]){
    if (j > i){
      res <- data.frame(node.x = c(i), node.y = c(j),
                       diff = c(trait_data_enrich_only$log2fold[i] - trait_data_enrich_only$log2fold[j]),
                       err = c(trait_data_enrich_only$glm.nb_stderr[i] + trait_data_enrich_only$glm.nb_stderr[j]),
                       phy_dist = c(node_dist_mat[i,j]))
      corr_data_long <- rbind(corr_data_long, res)
    }
  }  
}
corr_data_long %>%
  mutate(abs_diff = abs(diff)) -> corr_data_long

m0 <- lm(abs_diff ~ phy_dist, data=corr_data_long)
sum_m0 <- summary(m0)
coefs_m0 <- coef(m0)

corr_data_long %>%
  ggplot(aes(x = phy_dist, y = abs_diff)) +
  geom_point() +
  geom_linerange(aes(ymin=abs_diff-err, ymax=abs_diff+err)) +
  geom_abline(slope = coefs_m0[2],
              intercept = coefs_m0[1]) +
  annotate(geom="text", x = 0.05, y = -15,
           label=paste0("Difference = ", round(coefs_m0[1], 2), " + ", round(coefs_m0[2], 4), "*Distance\n",
                        "p = ", round(sum_m0$coefficients[2,4], 4)),
           hjust = 0) +
  labs(x = "Phylogenetic distance (expected subst./site)",
       y = "Abs. Value of Difference in Log2 Fold Enrichment")-> g
g
ggsave("../figures/difference_distance.svg", g, height = 6, width = 6)

m1 <- phyloSignal(p4d_tree)
m1

m2 <- phyloSignalINT(p4d_tree)
tdata(m2)

filt_enrich

m1 <- cor_phylo(variates = ~log2fold,
                species = ~label,
                phy=phylotree,
                meas_errors = glm.nb_stderr,
                data = tips,
)
trait_data
### Genome-scqle data from Li et al. 2021 https://doi.org/10.1016/j.cub.2021.01.074

genome_tree <- read.newick("../data/trees/3_phylorank/phylorank/1644taxa_290genes_bb_1.tre")
genome_tree %>%
  as.treedata() -> gtreedata

gtreedata %>%
  as_tibble() ->gtreedata_tibble
tidytree::MRCA(gtreedata, 10, 11)

genome_taxa_data <- read.csv("../data/trees/Table_S1_Li_etal2021.csv")
genome_taxa_data

uniq_taxa_df <- data.frame(taxon=c(),
                           rank=c(),
                           node=c())
for (rank in c("class", "order", "family.clade")){ #"genus")){
  unique_taxa <- unique(genome_taxa_data[[paste0(rank, ".linked.to.NCBI.record")]])
  unique_taxa <- unique_taxa[unique_taxa != "no_rank"]
  res_df <- data.frame(taxon = unique_taxa,
                       rank = rank,
                       node = NA)
  for (i in 1:length(unique_taxa)){
    print(unique_taxa[i])
    sub_taxa_list <- genome_taxa_data %>%
      filter(!!as.symbol(paste0(rank, ".linked.to.NCBI.record")) == unique_taxa[i])
    print(dim(sub_taxa_list))
    # print(gtreedata_tibble$node[gtreedata_tibble$label %in% sub_taxa_list$new_taxonID_in_Fig.1])
    nodes <- gtreedata_tibble$node[gtreedata_tibble$label %in% sub_taxa_list$new_taxonID_in_Fig.1]
    if (length(nodes) > 1){
      res_df$node[i] <- treeio::MRCA(gtreedata, nodes)
    }
  }
  uniq_taxa_df <- rbind(uniq_taxa_df, res_df)
}
write.csv(uniq_taxa_df, "../data/stats/unique_taxa_node_df.csv")
uniq_taxa_df <- read.csv("../data/stats/unique_taxa_node_df.csv")
enrich_data <- read.csv("../data/stats/glmnb_results_all_fungi.csv")

genome_tree <- read.newick("../data/trees/3_phylorank/phylorank/1644taxa_290genes_bb_1.tre")
genome_tree %>%
  as.treedata() -> gtreedata

gtreedata %>%
  as_tibble() ->gtreedata_tibble
# ggtree(gtreedata)
uniq_taxa_df %>%
  left_join(enrich_data, by = "taxon") -> gtree_taxa_data

gtree_taxa_data

gtreedata_tibble %>%
  left_join(gtree_taxa_data, by = "node") %>% as.treedata() -> gtreedata_w_data

uniq_taxa_df %>%
  filter(rank == "family.clade",
         !is.na(node)) %>%
  .$node -> fam_nodes
ggtree(gtreedata_w_data, aes(color=log2fold)) +
  collapse(fam_nodes) -> g
g
