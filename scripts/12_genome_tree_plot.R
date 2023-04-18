library(phytools)
library(phylobase)
library(ggtree)
library(treeio)
library(tidytree)
library(ggplot2)
library(tidyr)
library(tibble)
library(dplyr)
library(stringr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
uniq_taxa_df <- read.csv("../data/stats/unique_taxa_node_df.csv")
enrich_data <- read.csv("../data/stats/glmnb_results_all_fungi.csv")

genome_tree <- treeio::read.newick("../data/trees/3_phylorank/phylorank/1644taxa_290genes_bb_1.tre")
genome_tree -> gtreedata
gtreedata

gtreedata %>%
  as_tibble() ->gtreedata_tibble
# ggtree(gtreedata, layout = "radial")
uniq_taxa_df %>%
  left_join(enrich_data, by = "taxon") %>%
  dplyr::select(-X.y) %>%
  rename(X = X.x) -> gtree_taxa_data

gtree_taxa_data

ggtree(gtreedata, aes(color = log2fold), layout = "fan")+
  scale_color_distiller(type = "div") +
  scale_alpha_continuous(limits = c(0.5, 1)) -> g
g$data %>%
  left_join(gtree_taxa_data, by = "node") %>%
  mutate(log2fold = case_when(is.na(log2fold) ~ 0,
                              TRUE ~ log2fold))-> taxa_data_coords
g$data <- taxa_data_coords
g$data %>%
  filter(taxon %in% c("Eurotiomycetes",
                      "Dothideomycetes",
                      "Tremellomycetes",
                      "Microbotryomycetes",
                      "Saccharomycetes",
                      "Leotiomycetes",
                      "Agaricomycetes",
                      "Sordariomycetes")) %>%
  dplyr::select(node, taxon) %>%
  rename(taxon_lab = taxon)-> phy_labels
g2 <- g + geom_cladelab(data = phy_labels,
                       mapping = aes(node = node,
                           label = taxon_lab))
  
g2
ggsave("../figures/genome_tree_with_log2fold.png", g2, width = 10, height = 10)
ggsave("../figures/genome_tree_with_log2fold.svg", g2, width = 10, height = 10)

# ggtree(gtreedata) -> g
# g_data <- g$data
g + geom_point(data = taxa_data_coords, aes(x = x, y = y, color = log2fold))
#   geom_nodepoint(data = gtreedata_w_data,
#                  aes(color = log2fold))
# g

gtreedata %>%
  as_tibble() %>%
  left_join(gtree_taxa_data, by = "node") %>%
  dplyr::select(parent, node, label, X, taxon, branch.length,log2fold, glm.nb_stderr) %>%
  as.treedata() -> gtreedata_w_data
gtreedata_w_data

# write.beast(gtreedata_w_data, file="../data/trees/genome_tree_with_log2fold.tre",
#             )

ancestors <- c()
for (i in uniq_taxa_df$node){
  if (! is.na(i)){
    nodes <- treeio::ancestor(gtreedata, i)
    ancestors <- c(ancestors, nodes)
  }
}
uniq_ance <- unique(ancestors)
uniq_ance

gtreedata_w_data %>% as_tibble() %>% filter(parent == node)

gtreedata_w_data %>%
  as_tibble() %>%
  dplyr::filter(!is.na(taxon) | node %in% uniq_ance)%>%
  dplyr::select(-X) %>%
  rownames_to_column(var = "X")-> tree_to_fam
x <- vector(length = length(unique(tree_to_fam$node)))
tree_to_fam_dedup <- data.frame(X = x, parent = NA, node = NA,
                                branch.length = NA, label = NA,
                                taxon = NA, log2fold=NA, glm.nb_stderr=NA)
x <- 1
for (i in unique(tree_to_fam$node)){
  rows <- tree_to_fam[tree_to_fam$node == i,]
  tree_to_fam_dedup[x,] <- rows[1,]
  tree_to_fam_dedup$taxon[x] <- paste(unique(rows$taxon), collapse = "-")
  tree_to_fam_dedup$log2fold[x] <- mean(rows$log2fold)
  tree_to_fam_dedup$glm.nb_stderr[x] <- mean(rows$glm.nb_stderr)
  x <- x + 1
}
tree_to_fam_dedup
tree_to_fam_dedup[duplicated(tree_to_fam_dedup$node),]
tree_to_fam[duplicated(tree_to_fam$node),]

tree_to_fam_dedup %>% as.treedata() -> treedata_to_fam
tree_to_fam

# X_to_keep <- c()
# for (i in unique(tree_to_fam$node)){
#   X_to_keep <- c(X_to_keep, tree_to_fam$X[tree_to_fam$node == i][1])
# }
# tree_to_fam %>%
#   filter(X %in% X_to_keep) %>%
#   as.treedata() %>%
#   as.phylo()  ->treedata_to_fam
# treedata_to_fam %>% as_tibble()
genome_node_dist_mat <- dist.nodes(gtreedata)
genome_node_dist_mat[1645, 1646]
gtreedata$edge.length

trait_correlation_plot <- function(trait_data_enrich, dist_mat, rank_name, plot=T){
  if (dim(trait_data_enrich)[1] < 2){
    print("No rows in data")
  }
  else{
    i <- dim(trait_data_enrich)[1]
    x <- vector(length = ifelse(i%%2 == 1, i*(floor(0.5*i)), i*(0.5*i-1) + 0.5*i))
    corr_data_long <- data.frame(node.x = x, node.y = x,
                                 diff = x, err = x, phy_dist = x,
                                 taxon.x = x, taxon.y = x)
    x <- 1
    for (i in 1:dim(trait_data_enrich)[1]){
      for (j in 1:dim(trait_data_enrich)[1]){
        if (j > i){
          corr_data_long$node.x[x] <- trait_data_enrich$node[i]
          corr_data_long$node.y[x] <- trait_data_enrich$node[j]
          corr_data_long$diff[x] <- trait_data_enrich$log2fold[i] - trait_data_enrich$log2fold[j]
          corr_data_long$phy_dist[x] <- dist_mat[trait_data_enrich$node[i],
                                                 trait_data_enrich$node[j]]
          corr_data_long$err[x] <-trait_data_enrich$glm.nb_stderr[i] + trait_data_enrich$glm.nb_stderr[j]
          corr_data_long$taxon.x[x] <- trait_data_enrich$taxon[i]
          corr_data_long$taxon.y[x] <- trait_data_enrich$taxon[j]
          x <- x + 1
        }
      }  
    }
    corr_data_long %>%
      mutate(abs_diff = abs(diff)) -> corr_data_long
    m0 <- lm(abs_diff ~ phy_dist, data=corr_data_long)
    sum_m0 <- summary(m0)
    coefs_m0 <- coef(m0)
    if (plot){
      corr_data_long %>%
        ggplot(aes(x = phy_dist, y = abs_diff)) +
        geom_point(alpha = 0.1) +
        geom_smooth(method = lm, se = T) +
        annotate(geom="text", x = 0.05, y = max(corr_data_long$abs_diff),
                 vjust=1,
                 label=paste0("Difference = ", round(coefs_m0[1], 2),
                              " + ", round(coefs_m0[2], 4), "*Distance\n",
                              "p = ", formatC(sum_m0$coefficients[2,4], format="e", digits=2),
                              "\nn = ", length(corr_data_long$phy_dist)),
                 hjust = 0) +
        labs(x = paste0("Phylogenetic distance of ", rank_name, "-level clades"),
             y = "Abs. Value of Difference in Log2 Fold Enrichment")-> g
      g
      ggsave(paste0("../figures/genome_difference_distance_", rank_name, "_6x8.svg"), g, height = 6, width = 8)
      ggsave(paste0("../figures/genome_difference_distance_", rank_name, "_6x6.svg"), g, height = 6, width = 6)
      ggsave(paste0("../figures/genome_difference_distance_", rank_name, "_6x4.svg"), g, height = 6, width = 4)
    }
    return(corr_data_long)
  }
}

tree_to_fam %>%
  filter(!is.na(log2fold)) %>% dim(.)
  trait_correlation_plot(genome_node_dist_mat,
                       "all") -> corr_df_all
m0 <- lm(abs_diff ~ phy_dist, data=corr_df_all)
summary(m0)
confint(m0)
tree_to_fam %>%
  filter(!is.na(log2fold),
         str_detect(taxon, "aceae")) %>%
  trait_correlation_plot(genome_node_dist_mat,
                         "family") -> corr_df_family
tree_to_fam %>%
  filter(!is.na(log2fold),
         str_detect(taxon, "ales")) %>%
  trait_correlation_plot(genome_node_dist_mat,
                         "order") -> corr_df_order
tree_to_fam %>%
  filter(!is.na(log2fold),
         str_detect(taxon, "mycetes")) %>%
  trait_correlation_plot(genome_node_dist_mat,
                         "class") -> corr_df_class
# # Type of links in corr_df_all

# Because nodes can represent multiple taxa, this analysis is unlikely to work
link_type <- function(x, y, pattern=F){
  if (pattern){
    a <- str_replace_all(c(x,y), ".*aceae", "family")
    a <- str_replace_all(a, ".*ales", "order")
    a <- str_replace_all(a, ".*mycetes", "class")
    a <- str_replace_all(a, ".*mycotina", "subpylum")
  }
  else{
    a <- str_replace_all(c(x, y), ".clade", "")
  }
  a <- sort(a)
  return(str_to_title(paste(a[1], a[2], sep="-")))
}
uniq_taxa_df %>%
  filter(!is.na(node)) -> node_lookup
unique(node_lookup$node)
node_lookup_new <- data.frame(X = c(), taxon = c(), rank = c(), node = c())
for (i in unique(node_lookup$node)){
  rows <- node_lookup[node_lookup$node == i,]
  node_lookup_new <- rbind(node_lookup_new, rows[1,])
}
duplicated(node_lookup_new$node)
node_lookup_new[duplicated(node_lookup_new$node),] %>%
  arrange(node)

link_type("family.clade", "family", pattern=F)
link_type("Mucorales", "Mucoraceae", pattern=T)
corr_df_all %>%
  left_join(node_lookup_new %>%
              dplyr::select(node, rank)%>%
              rename(node.x = node, rank.x=rank),
            by = "node.x") %>%
  left_join(node_lookup_new %>%
              dplyr::select(node, rank)%>%
              rename(node.y = node, rank.y=rank),
            by = "node.y") %>%
  rowwise()%>%
  mutate(link = link_type(rank.x, rank.y)) -> corr_df_all_links
corr_df_all %>%
  rowwise() %>%
  mutate(link = link_type(taxon.x, taxon.y, pattern=T)) -> corr_df_all_links_pattern

corr_df_all_links
m0 <- lm(abs_diff ~ phy_dist, data=corr_df_all)
sum_m0 <- summary(m0)
coefs_m0 <- coef(m0)
# 
# corr_df_all
# corr_df_all_links
corr_df_all_links %>%
  ggplot(aes(x = phy_dist, y = abs_diff)) +
  geom_point(aes(color = factor(link)), alpha = 0.1) +
  # geom_linerange(aes(ymin=abs_diff-err, ymax=abs_diff+err),
  #                alpha = 0.05) +
  geom_smooth(method = lm, se = T) +
  annotate(geom="text", x = 0.05, y = max(corr_df_all$abs_diff),
           vjust=1,
           label=paste0("Difference = ", round(coefs_m0[1], 2),
                        " + ", round(coefs_m0[2], 4), "*Distance\n",
                        "p = ", formatC(sum_m0$coefficients[2,4], format="e", digits=2),
                        "\nn = ", length(corr_df_all$phy_dist)),
           hjust = 0) +
  labs(x = "Phylogenetic distance of all clades",
       y = "Abs. Value of Difference in Log2 Fold Enrichment",
       color = "Taxa relation") +
  guides(color = guide_legend(override.aes = list(alpha = 1)))-> g
g
ggsave("../figures/genome_difference_distance_all_colored_links_6x8.svg", g, height = 6, width = 8)

corr_df_all_links %>%
  ggplot(aes(x = phy_dist, y = abs_diff)) +
  geom_point(alpha = 0.1) +
  # geom_linerange(aes(ymin=abs_diff-err, ymax=abs_diff+err),
  #                alpha = 0.05) +
  geom_smooth(method = lm, se = T) +
  facet_wrap(~factor(link), nrow =2, scales = "free_x") +
  labs(x = "Phylogenetic distance of all clades",
       y = "Abs. Value of Difference in Log2 Fold Enrichment") ->g
g
ggsave("../figures/genome_difference_distance_all_grid_links_6x8.svg", g, height = 6, width = 8)
# Now using the pattern
corr_df_all_links %>%
  ggplot(aes(x = phy_dist, y = abs_diff)) +
  geom_point(aes(color = factor(link)), alpha = 0.1) +
  # geom_linerange(aes(ymin=abs_diff-err, ymax=abs_diff+err),
  #                alpha = 0.05) +
  geom_smooth(method = lm, se = T) +
  annotate(geom="text", x = 0.05, y = max(corr_df_all$abs_diff),
           vjust=1,
           label=paste0("Difference = ", round(coefs_m0[1], 2),
                        " + ", round(coefs_m0[2], 4), "*Distance\n",
                        "p = ", formatC(sum_m0$coefficients[2,4], format="e", digits=2),
                        "\nn = ", length(corr_df_all$phy_dist)),
           hjust = 0) +
  labs(x = "Phylogenetic distance of all clades",
       y = "Abs. Value of Difference in Log2 Fold Enrichment",
       color = "Taxa relation") +
  guides(color = guide_legend(override.aes = list(alpha = 1)))-> g
g
ggsave("../figures/genome_difference_distance_all_colored_links_by_pattern_6x8.svg", g, height = 6, width = 8)

corr_df_all_links %>%
  ggplot(aes(x = phy_dist, y = abs_diff)) +
  geom_point(alpha = 0.1) +
  # geom_linerange(aes(ymin=abs_diff-err, ymax=abs_diff+err),
  #                alpha = 0.05) +
  geom_smooth(method = lm, se = T) +
  facet_wrap(~factor(link), nrow =2, scales="free_x") +
  labs(x = "Phylogenetic distance of all clades",
       y = "Abs. Value of Difference in Log2 Fold Enrichment") ->g
g
ggsave("../figures/genome_difference_distance_all_grid_links_by_pattern_6x8.svg", g, height = 6, width = 8)

# Why are the results different from the filter by name?
corr_df_family %>%
  rowwise()%>%
  mutate(link = link_type(taxon.x, taxon.y, pattern=T)) -> corr_df_family_links
corr_df_family_links %>%
  filter(link == "Family-Family")
corr_df_family_links %>%
  ggplot(aes(x = phy_dist, y = abs_diff)) +
  geom_point(color = "red", alpha=0.2) +
  geom_point(data = corr_df_all_links_pattern %>% filter(link == "Family-Family"),
             aes(x = phy_dist, y = abs_diff),
             color = "blue", alpha=0.2)


# Permutational test for null model. What if the trait was distributed randomly?
trait_randomization <- function(trait_data_enrich, dist_mat, rank_name, permutations=99, only_actual=F){
  trait_diff_mat <- matrix(nrow=dim(trait_data_enrich)[1],
                           ncol=dim(trait_data_enrich)[1])
  colnames(trait_diff_mat) <- trait_data_enrich$node
  rownames(trait_diff_mat) <- trait_data_enrich$node
  trait_err_mat <- trait_diff_mat
  if (dim(trait_data_enrich)[1] < 2){
    print("No rows in data")
  }
  else{
    if (! only_actual){
      slopes <- vector(length = permutations)
      intercepts <- slopes
      for (p in 1:permutations){
        i <- dim(trait_data_enrich)[1]
        x <- vector(length = ifelse(i%%2 == 1, i*(floor(0.5*i)), i*(0.5*i-1) + 0.5*i))
        corr_data_long <- data.frame(diff = x, phy_dist = x)
        r_log2fold <- sample(trait_data_enrich$log2fold)
        x <- 1
        for (i in 1:dim(trait_data_enrich)[1]){
          for (j in 1:dim(trait_data_enrich)[1]){
            if (j > i){
              corr_data_long$diff[x] <- r_log2fold[i] - r_log2fold[j]
              corr_data_long$phy_dist[x] <- dist_mat[trait_data_enrich$node[i],trait_data_enrich$node[j]]
              x <- x + 1
            }
          }  
        }
        corr_data_long %>%
          mutate(abs_diff = abs(diff)) -> corr_data_long
        m0 <- lm(abs_diff ~ phy_dist, data=corr_data_long)
        coefs_m0 <- coef(m0)
        intercepts[p] <- coefs_m0[1]
        slopes[p] <- coefs_m0[2]
        print(paste(p, coefs_m0))
      }
      i <- dim(trait_data_enrich)[1]
      x <- vector(length = ifelse(i%%2 == 1, i*(floor(0.5*i)), i*(0.5*i-1) + 0.5*i))
      corr_data_long <- data.frame(diff = x, phy_dist = x)
      x <- 1
      for (i in 1:dim(trait_data_enrich)[1]){
        for (j in 1:dim(trait_data_enrich)[1]){
          if (j > i){
            corr_data_long$diff[x] <- trait_data_enrich$log2fold[i] - trait_data_enrich$log2fold[j]
            corr_data_long$phy_dist[x] <- dist_mat[trait_data_enrich$node[i],trait_data_enrich$node[j]]
            x <- x + 1
          }
        }  
      }
      corr_data_long %>%
        mutate(abs_diff = abs(diff)) -> corr_data_long
      m0 <- lm(abs_diff ~ phy_dist, data=corr_data_long)
      coefs_m0 <- coef(m0)
      intercept_act <- coefs_m0[1]
      slope_act <- coefs_m0[2]
      res_df <- data.frame(intercepts = intercepts,
                           slopes = slopes)
      g_int <- ggplot(res_df, aes(x = intercepts)) +
        geom_histogram() + 
        geom_vline(xintercept = intercept_act)
      g_slope <- ggplot(res_df, aes(x = slopes)) +
        geom_histogram() + 
        geom_vline(xintercept = slope_act)
      ggsave(paste0("../figures/diff_distance_permutations_intercept_", rank_name, "_6x6.svg"), g_int, height = 6, width = 6)
      ggsave(paste0("../figures/diff_distance_permutations_slope_", rank_name, "_6x6.svg"), g_slope, height = 6, width = 6)
      return(res_df)
    }
    else {
      i <- dim(trait_data_enrich)[1]
      x <- vector(length = ifelse(i%%2 == 1, i*(floor(0.5*i)), i*(0.5*i-1) + 0.5*i))
      corr_data_long <- data.frame(diff = x, phy_dist = x)
      for (i in 1:dim(trait_data_enrich)[1]){
        for (j in 1:dim(trait_data_enrich)[1]){
          if (j > i){
            corr_data_long$diff[x] <- trait_data_enrich$log2fold[i] - trait_data_enrich$log2fold[j]
            corr_data_long$phy_dist[x] <- dist_mat[trait_data_enrich$node[i],trait_data_enrich$node[j]]
            x <- x + 1
          }
        }  
      }
      corr_data_long %>%
        mutate(abs_diff = abs(diff)) -> corr_data_long
      m0 <- lm(abs_diff ~ phy_dist, data=corr_data_long)
      coefs_m0 <- coef(m0)
      intercept_act <- coefs_m0[1]
      slope_act <- coefs_m0[2]
      return(coefs_m0)
    }
  }
}
tree_to_fam %>%
  filter(!is.na(log2fold),
         str_detect(taxon, "aceae")) %>%
  trait_randomization(genome_node_dist_mat,
                         "family", permutations = 999) -> family_perm_df
write.csv(family_perm_df, "../data/stats/family_permuation_distance_difference.csv")
tree_to_fam %>%
  filter(!is.na(log2fold),
         str_detect(taxon, "aceae")) %>%
  trait_randomization(genome_node_dist_mat,
                      "family", only_actual = T) -> family_coefs
family_perm_df %>%
  ggplot(aes(x=slopes)) +
  geom_histogram() +
  geom_vline(xintercept = family_coefs[2]) -> g1
family_perm_df %>%
  ggplot(aes(x=intercepts)) +
  geom_histogram() +
  geom_vline(xintercept = family_coefs[1]) -> g2
ggsave("../figures/diff_distance_permutations_slope_family_3x4.svg", g1, height = 3, width = 4)
ggsave("../figures/diff_distance_permutations_intercept_family_3x4.svg", g2, height = 3, width = 4)

## Order
tree_to_fam %>%
  filter(!is.na(log2fold),
         str_detect(taxon, "ales")) %>%
  trait_randomization(genome_node_dist_mat,
                      "order", permutations = 999) -> order_perm_df
write.csv(order_perm_df, "../data/stats/order_permuation_distance_difference.csv")

tree_to_fam %>%
  filter(!is.na(log2fold),
         str_detect(taxon, "ales")) %>%
  trait_randomization(genome_node_dist_mat,
                      "order", only_actual = T) -> order_coefs
order_perm_df %>%
  ggplot(aes(x=slopes)) +
  geom_histogram() +
  geom_vline(xintercept = order_coefs[2]) -> g1
order_perm_df %>%
  ggplot(aes(x=intercepts)) +
  geom_histogram() +
  geom_vline(xintercept = order_coefs[1]) -> g2

ggsave("../figures/diff_distance_permutations_slope_order_3x4.svg", g1, height = 3, width = 4)
ggsave("../figures/diff_distance_permutations_intercept_order_3x4.svg", g2, height = 3, width = 4)
# Class
tree_to_fam %>%
  filter(!is.na(log2fold),
         str_detect(taxon, "mycetes")) %>%
  trait_randomization(genome_node_dist_mat,
                      "class", permutations = 999) -> class_perm_df
write.csv(class_perm_df, "../data/stats/class_permuation_distance_difference.csv")
class_perm_df <- read.csv("../data/stats/class_permuation_distance_difference.csv")

tree_to_fam %>%
  filter(!is.na(log2fold),
         str_detect(taxon, "ales")) %>%
  trait_randomization(genome_node_dist_mat,
                      "class", only_actual = T) -> class_coefs
class_perm_df %>%
  ggplot(aes(x=slopes)) +
  geom_histogram() +
  geom_vline(xintercept = class_coefs[2]) -> g1
class_perm_df %>%
  ggplot(aes(x=intercepts)) +
  geom_histogram() +
  geom_vline(xintercept = class_coefs[1]) -> g2

sum(class_perm_df$intercepts < class_coefs[1]) / length(class_perm_df$intercepts)
ggsave("../figures/diff_distance_permutations_slope_class_3x4.svg", g1, height = 3, width = 4)
ggsave("../figures/diff_distance_permutations_intercept_class_3x4.svg", g2, height = 3, width = 4)


# Enrichment bar chart

parent_df <- enrich_data %>%
  dplyr::select(taxon, log2fold, glm.nb_stderr, glm.nb_zval)
uniq_taxa_df %>%
  dplyr::select(node, taxon, rank) %>%
  filter(!is.na(node)) %>%
  rename_with(~paste0("parent.", .x, recycle0=T)) -> parent_nodes
parent_nodes

parent_df %>%
  left_join(uniq_taxa_df, by = "taxon") %>%
  left_join(gtreedata_tibble, by = "node") %>%
  rename(parent.node = parent) %>%
  left_join(parent_nodes,
            by = "parent.node") -> parent_df
enrich_data

uniq_taxa_df
uniq_taxa_df %>%
  rename_with(~paste0("parent.", .x, recycle0=T),
              cols = c(taxon, rank))

gtreedata_tibble$parent[]

enrich_data %>%
  # filter(str_detect(taxon, "ales")) %>%
  filter(str_detect(taxon, "mycetes")) %>%
  arrange(desc(abs(glm.nb_zval))) %>%
  head(20) -> top20_class

top20_class$Phylum <- c("Ascomycota",
                        "Ascomycota",
                        "Basidiomycota",
                        "Basidiomycota",
                        "Ascomycota",
                        "Basidiomycota",
                        "Ascomycota",
                        "Mucuromycota",
                        "Mucuromycota",
                        "Basidiomycota",
                        "Ascomycota",
                        "Chytridiomycota",
                        "Ascomycota",
                        "Chytridiomycota",
                        "Chytridiomycota",
                        "Ascomycota",
                        "Ascomycota",
                        "Ascomycota",
                        "Chytridiomycota",
                        "Basidiomycota")

enrich_data %>%
  filter(str_detect(taxon, "aceae")) %>%
  arrange(desc(abs(glm.nb_zval))) %>%
  head(20)  %>%
  ggplot(aes(x = factor(taxon, levels = taxon[order(log2fold)]), y = log2fold)) +
  geom_col() +
  geom_errorbar(aes(ymin = log2fold - glm.nb_stderr,
                    ymax = log2fold + glm.nb_stderr)) +
  theme(axis.text.x.bottom = element_text(angle = -45,
                                          hjust = 0),
        plot.margin = margin(0,0.5,0,0,unit="in")) +
  labs(x = "Family",
       y = "Log2 fold-enrichment in phyllosphere") -> g
g
ggsave("../figures/barchart_enrich_family_top20.svg", g, height = 5, width = 7)

top20_class %>%
  ggplot(aes(x = factor(taxon, levels = taxon[order(log2fold)]), y = log2fold)) +
  geom_col(aes(fill = Phylum)) +
  geom_errorbar(aes(ymin = log2fold - glm.nb_stderr,
                    ymax = log2fold + glm.nb_stderr)) +
  theme(axis.text.x.bottom = element_text(angle = -45,
                                          hjust = 0)) +
  labs(x = "Class",
       y = "Log2 fold-enrichment in phyllosphere") -> g
g
ggsave("../figures/barchart_enrich_class_top20.svg", g, height = 5, width = 7)
