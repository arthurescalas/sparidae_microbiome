################################################################################
#         _______  _____  _____ ___ ____  _   _ __  __ _____ ____  
#        | ____\ \/ / _ \|  ___|_ _/ ___|| | | |  \/  | ____|  _ \ 
#        |  _|  \  / | | | |_   | |\___ \| |_| | |\/| |  _| | | | |
#        | |___ /  \ |_| |  _|  | | ___) |  _  | |  | | |___| |_| |
#        |_____/_/\_\___/|_|   |___|____/|_| |_|_|  |_|_____|____/ 
#
################################################################################
#
# R script to compute dissimilarity between fish species in terms of:
#   - phylogenetic distance
#   - gut morphology
#   - diet composition
#   - gut µbiome
#
# Arthur Escalas July 2020
# arthur.escalas@gmail.com
################################################################################



################################ LOAD DATA #####################################


#  Load the CORE phyloseq objects ----

ps_data <- readRDS(paste0(dir_data, "phyloseq_object_core.rds"))


#   List the ranks of bacterial taxonomy ----

nms_ranks <- c("Phylum", "Class", "Order", "Family", "Genus", "OTU")

#   list the differents traits ----

nms_traits_morpho <- c("Gut_length", "Relative_gut_length")

# Load the metadata ----

metadata_sp <- read.csv(paste0(dir_data, "table_metadata_species.csv"), row.names = 1)
metadata <- sample_data(ps_data) %>% data.frame()



################# ESTIMATE DISSIMILARITY BETWEEN FISH SPECIES ##################


# ================== COMPUTE FISH PHYLOGENETIC DISSIMILIARITY =================

#     read the phylogenetic tree ----

phylogeny_santini <- read.tree(paste0(dir_data, "phylogeny_sparid_santini_2013.phy"))

#     extract tip nodes, i.e. species  ----

sp_santini <- phylogeny_santini$tip.label
sp_µbiome <- unique(metadata$Species)

#     prune the fish phylogenetic tree ----

phylogeny_pruned <- drop.tip(phy = phylogeny_santini, 
                             tip = sp_santini[! sp_santini %in% sp_µbiome])
saveRDS(phylogeny_pruned, file = paste0(dir_dissimilarity, "phylogeny_pruned.rds"))

#   get the cophenetic distance ----
diss_fish <- cophenetic(phylogeny_pruned) %>% as.dist()

# names of species odered according to tree
nms_fish_sp_ordered_tree <- phylogeny_pruned$tip.label

#     reformat dissimilarity as table  ----

df_diss_fish_phylo <- dendextend::dist_long(diss_fish) %>% 
  rename(distance = "phylo")
df_diss_fish_phylo$sp_comb <- sapply(1:nrow(df_diss_fish_phylo), function(x) { 
  tmp <- sort(c(as.character(df_diss_fish_phylo$rows)[x], 
                as.character(df_diss_fish_phylo$cols)[x]))
  paste0(tmp[1], "__", tmp[2])
})
df_diss_fish_phylo$phylo <- df_diss_fish_phylo$phylo/2
df_diss_fish_phylo <- df_diss_fish_phylo %>% dplyr::select(sp_comb, phylo)



# =================== COMPUTE FISH FUNCTIONAL DISSIMILIARITY ==================


# ----------------------------- DIET DISSIMILARITY -----------------------------

#   Load diet data from fishbase ----

filename <- paste0(dir_data, "table_diet_composition.csv")
data_diet <- read.csv(filename, encoding = "UTF-7", row.names = 1)

#   Compute diet dissimilarity ----

diss_fish_diet  <- vegdist(data_diet, "bray")


# ---------------------- MORHOLOGICAL DISSIMILARITY ----------------------------

diss_fish_morpho <- lapply(nms_traits_morpho, function(trt) {
  dat <- scale(metadata_sp[, trt])
  row.names(dat) <- row.names(metadata_sp)
  vegdist(dat, method = "euclidean")
}) %>% setNames(nms_traits_morpho)


# --------------------- REFORMAT DISSIMILARITY AS TABLE ------------------------

tmp <- dendextend::dist_long(diss_fish_morpho[[1]]) %>% dplyr::select(- distance)

df_diss_diet <- dendextend::dist_long(diss_fish_diet) %>% 
  rename(distance = "Diet_FishBase") %>% 
  dplyr::select(Diet_FishBase)

df_diss_morpho <- lapply(names(diss_fish_morpho), function(x) {
  dendextend::dist_long(diss_fish_morpho[[x]]) %>% rename(distance = x) %>% 
    dplyr::select(all_of(x))
}) %>% bind_cols()

df_diss_fish <- cbind(tmp, df_diss_morpho, df_diss_diet)
df_diss_fish$sp_comb <- sapply(1:nrow(df_diss_fish), function(x) {
  tmp <- sort(c(as.character(df_diss_fish$rows)[x],
                as.character(df_diss_fish$cols)[x]))
  paste0(tmp[1], "__", tmp[2])
})



# ====================== EXPORT DATA ON FISH DISSMILARITY ======================


df_diss_fish <- left_join(df_diss_fish, df_diss_fish_phylo, "sp_comb") %>%
  dplyr::select(rows, cols, sp_comb, phylo, everything()) %>%
  rename(rows = "species_a", cols = "species_b")

write.csv(df_diss_fish, file = paste0(dir_dissimilarity, "table_fish_dissimilarity.csv"),
          row.names = FALSE)



# ================= TEST PHYLOGENETIC CONSERVATISM OF TRAITS ===================
# using the abouheif.moran() test from adephylo package
# The test of Abouheif (1999) is designed to detect phylogenetic autocorrelation 
# in a quantitative trait.
# https://academic.oup.com/bioinformatics/article/26/15/1907/188748

dat <- phylobase::phylo4d(as(phylogeny_pruned, "phylo4"), 
                          metadata_sp[, nms_traits_morpho])
test <- abouheif.moran(dat, method="Abouheif")

# class: krandtest lightkrandtest 
# Monte-Carlo tests
# Call: as.krandtest(sim = matrix(res$result, ncol = nvar, byrow = TRUE), 
#                    obs = res$obs, alter = alter, names = test.names)
# 
# Number of tests:   2 
# 
# Adjustment method for multiple comparisons:   none 
# Permutation number:   999 
#                  Test       Obs  Std.Obs   Alter Pvalue
# 1          Gut_length 0.1889600 2.250782 greater  0.027
# 2 Relative_gut_length 0.3371291 2.818577 greater  0.021

png(paste0(dir_dissimilarity, "plot_phylogenetic_conservatism_of_gut_traits.png"),
    height = 20, width = 20, unit = "cm", res = 200)
table.phylo4d(dat)
graphics.off()


###################### ESTIMATE MICROBIOME DISSIMILARITY #######################


# ============ Make lists of data at different microbial ranks =================

#     for individual samples ----

ls_glom <- lapply(nms_ranks, function(rk) {
  tax_glom(ps_data, taxrank = rk)
}) %>% setNames(nms_ranks)


#     for species ----

ps_sp <- ps_data %>% merge_samples(group = "Species", fun = mean)
ps_sp@otu_table <- ps_sp@otu_table %>% t()

ls_glom_sp <- lapply(nms_ranks, function(rk) {
  tax_glom(ps_sp, taxrank = rk)
}) %>% setNames(nms_ranks)


# ====================== Estimate the dissimilarity ============================

#     for individual samples ----

ls_diss <- list()

for (rk in nms_ranks) {
  ls_diss[[rk]] <- get_hill_numbers_dissimilarities(ls_glom[[rk]], dir_tmp = dir_dissimilarity)
}

ls_diss <- unlist(ls_diss, recursive = FALSE)
names(ls_diss) <- gsub("\\.", "_", names(ls_diss))

saveRDS(ls_diss, paste0(dir_dissimilarity, "list_microbiome_dissimilarity_inter_samples.rds"))

#     for species ----

ls_diss_sp <- list()

for (rk in nms_ranks) {
  ls_diss_sp[[rk]] <- get_hill_numbers_dissimilarities(ls_glom_sp[[rk]], dir_tmp = dir_dissimilarity)
}

ls_diss_sp <- unlist(ls_diss_sp, recursive = FALSE)
names(ls_diss_sp) <- gsub("\\.", "_", names(ls_diss_sp))

saveRDS(ls_diss_sp, paste0(dir_dissimilarity, "list_microbiome_dissimilarity_inter_species.rds"))



# ================== REFORMAT DISSMILARITY AS A TABLE ==========================

#  -----------------------  for individual samples -----------------------------

nm_methods <- names(ls_diss)

#  make a table of dissmilarity between samples 

tmp <- dendextend::dist_long(ls_diss[[1]]) %>% dplyr::select(- distance)

df <- lapply(nm_methods, function(x) {
  dendextend::dist_long(ls_diss[[x]]) %>% rename(distance = x) %>% 
    dplyr::select(x)
}) %>% bind_cols()

df <- cbind(tmp, df)
df$spl_comb <- sapply(1:nrow(df), function(x) { 
  tmp <- sort(c(as.character(df$rows)[x], 
                as.character(df$cols)[x]))
  paste0(tmp[1], "__", tmp[2])
})
df$sp_comb <- sapply(df$spl_comb, function(x) {
  sapply(str_split_fixed(x, "__", 2), function(xx) {
    metadata[metadata$sample_id_fastq == xx, "Species"]
  }) %>% sort() %>% paste0(collapse = "__")
})

df_diss_µ <- df %>% 
  separate(sp_comb, c("species_a", "species_b"), sep = "__", remove = FALSE) %>%
  rename(rows = "sample_a", cols = "sample_b") %>% 
  dplyr::select(sample_a, sample_b, spl_comb, species_a, species_b, sp_comb,
                everything())

write.csv(df_diss_µ, file = paste0(dir_dissimilarity, "table_microbiome_dissimilarity_samples.csv"),
          row.names = FALSE)


#  ------------------------- for species ---------------------------------------

tmp <- dendextend::dist_long(ls_diss_sp[[1]]) %>% dplyr::select(- distance)

df <- lapply(nm_methods, function(x) {
  dendextend::dist_long(ls_diss_sp[[x]]) %>% rename(distance = x) %>% 
    dplyr::select(x)
}) %>% bind_cols()

df <- cbind(tmp, df)
df$sp_comb <- sapply(1:nrow(df), function(x) { 
  tmp <- sort(c(as.character(df$rows)[x], 
                as.character(df$cols)[x]))
  paste0(tmp[1], "__", tmp[2])
})

df_diss_µ_sp <- df %>% 
  rename(rows = "species_a", cols = "species_b") %>% 
  dplyr::select(species_a, species_b, sp_comb, everything())

write.csv(df_diss_µ_sp, file = paste0(dir_dissimilarity, "table_microbiome_dissimilarity_species.csv"),
          row.names = FALSE)




############### ESTIMATE MICROBIOME DISSIMILARITY WITHIN TAXA ##################


# ============================ PREPARE DATA ====================================

#     merge samples at the species level ------

ps <- ps_data %>% merge_samples(group = "Species", fun = mean)
ps@otu_table <- ps@otu_table %>% t()
tab_taxonomy <- ps@tax_table@.Data %>% as.data.frame()

#     make a list of data for various microbial taxa ----
# Levels of taxonomic resolution considered: Phylum, Class and Order
# THese levels were selected because significant relationships with traits were
# observed only at these ranks

ranks <- nms_ranks[c(1:3)]

ls_data <- list()

for (rk in ranks) { # Loop on each taxonomic level
  
  ls_data[[rk]] <- list()
  
  # check which taxa have at least 10 OTUs
  dims <- do.call(c, lapply(split(tab_taxonomy, tab_taxonomy[, rk]), nrow))
  taxa_to_use <- names(which(dims >= 10))
  
  for (taxa in taxa_to_use) {
    
    # Subset the phyloseq object to keep only OTU from the selected taxa
    eval(parse(text = paste0("tmp <- phyloseq::subset_taxa(ps,", rk, " == taxa)")))
    
    # Remove species that do not contain this taxa in their microbiome
    tmp <- subset_samples(tmp, colSums(tmp@otu_table) != 0)
    
    ls_data[[rk]][[taxa]] <- tmp
  }
} 


# ========================= ESTIMATE DISSIMILARITY =============================

ls_diss <- lapply(ls_data, function(X) {
  lapply(X, function(x) {
    get_hill_numbers_dissimilarities(x, dir_tmp = dir_dissimilarity)
  })
})

saveRDS(ls_diss, file = paste0(dir_dissimilarity, 
                               "list_microbiome_dissimilarity_matrices_per_taxa.rds"))


# ================== REFORMAT DISSMILARITY AS A TABLE ==========================

df <- lapply(ls_diss, function(X) {
  lapply(X, function(x) {
    lapply(x, function(xx) dendextend::dist_long(xx)) %>%
      reformat_as_df(new_var_name = "index")
  }) %>% reformat_as_df(new_var_name = "taxa_name")
}) %>%  reformat_as_df(new_var_name = "rank")

df$sp_comb <- sapply(1:nrow(df), function(x) { 
  tmp <- sort(c(as.character(df$rows)[x], 
                as.character(df$cols)[x]))
  paste0(tmp[1], "__", tmp[2])
})

df_diss <- df %>% 
  # separate(sp_comb, c("species_a", "species_b"), sep = "__", remove = FALSE) %>%
  rename(rows = "species_a", cols = "species_b") %>% 
  # Remove NA and "Inf" in dissimilarity values
  filter(distance != "Inf" & ! is.na(distance)) %>% 
  dplyr::select(species_a, species_b, sp_comb, rank, taxa_name, index,
                distance, everything())

write.csv(df_diss, file = paste0(dir_dissimilarity, 
                                 "table_microbiome_dissimilarity_per_taxa.csv"),
          row.names = FALSE)










