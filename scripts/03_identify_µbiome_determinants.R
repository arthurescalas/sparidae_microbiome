################################################################################
#         _______  _____  _____ ___ ____  _   _ __  __ _____ ____  
#        | ____\ \/ / _ \|  ___|_ _/ ___|| | | |  \/  | ____|  _ \ 
#        |  _|  \  / | | | |_   | |\___ \| |_| | |\/| |  _| | | | |
#        | |___ /  \ |_| |  _|  | | ___) |  _  | |  | | |___| |_| |
#        |_____/_/\_\___/|_|   |___|____/|_| |_|_|  |_|_____|____/ 
#
################################################################################
#
# R script to identify gut µbiome determinants
#
# Arthur Escalas July 2020
# arthur.escalas@gmail.com
################################################################################


################################ LOAD DATA #####################################

#     Load the phyloseq objects ----

ps_data <- readRDS(paste0(dir_data, "phyloseq_object_core.rds"))


# Load the metadata ----

metadata_sp <- read.csv(paste0(dir_data, "table_metadata_species.csv"), row.names = 1)
metadata <- sample_data(ps_data) %>% data.frame()

#   Load the dissimilarities ----

# inter fishes
df_diss_fish <- read.csv(paste0(dir_dissimilarity, "table_fish_dissimilarity.csv"))

# inter individuals
df_diss_fish_indiv <- read.csv(paste0(dir_dissimilarity, 
                                      "table_microbiome_dissimilarity_samples.csv"))

ls_diss <- readRDS(paste0(dir_dissimilarity, 
                          "list_microbiome_dissimilarity_inter_samples.rds"))


#   List the ranks of bacterial taxonomy ----

nms_ranks <- c("Phylum", "Class", "Order", "Family", "Genus", "OTU")

#   list the differents traits ----

nms_traits_morpho <- c("Gut_length", "Relative_gut_length")

nms_factors <- c("Species", "Diet_category")

# names of the methodsXrank combinations ----

nms_methods <- names(ls_diss)



# =============== TEST DIFFERENCE BETWEEN GROUPS (Wds test) ===================


theme_set(theme_classic(base_size = 16))

ls_res <- lapply(nms_factors, function(trt) {
  lapply(ls_diss, function(d) {
    tmp <- WdS.test(d, factor(metadata[, trt]), nrep = 999, strata = NULL)
    out <- data.frame(p_value = tmp$p.value, statistic = tmp$statistic)
    return(round(out, 3))
  }) %>% reformat_as_df(new_var_name = "data_type")
}) %>% setNames(nms_factors)


#     format the results as a table ----

df_res <- reformat_as_df(ls_res, new_var_name = "y_var")
df_res$rank = factor(str_split_fixed(df_res$data_type, "_", 2)[,1], level = nms_ranks)
df_res$method = factor(str_split_fixed(df_res$data_type, "_", 2)[,2])


#     Plot the results ----

png(paste0(dir_determinants, "plot_pvalues_wdstest.png"), height = 15, width = 10,
    unit = "cm", res = 200)
ggplot(df_res, aes(x = rank, y = p_value)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color="orangered", size=0.6, alpha=0.6, width = 0.2) +
  geom_abline(aes(intercept = 0.05, slope = 0, color = "red"), show.legend = FALSE) +
  geom_abline(aes(intercept = 0.01, slope = 0, color = "blue"), show.legend = FALSE) +
  facet_grid(y_var ~ .)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.background =element_blank())
graphics.off()

png(paste0(dir_determinants, "plot_pvalues_wdstest_logged.png"), height = 15, width = 10,
    unit = "cm", res = 200)
ggplot(df_res, aes(x = rank, y = log(p_value))) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color="orangered", size=0.6, alpha=0.6, width = 0.2) +
  geom_abline(aes(intercept = log(0.05), slope = 0, color = "red"), show.legend = FALSE) +
  geom_abline(aes(intercept = log(0.01), slope = 0, color = "blue"), show.legend = FALSE) +
  facet_grid(y_var ~ .)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.background = element_blank())
graphics.off()


# ======================== Summarize the results ===============================


out <- df_res %>% group_by(y_var) %>% 
  dplyr::summarize(avg_p = round(mean(p_value, na.rm = TRUE),3), 
                   sd_p = round(sd(p_value, na.rm = TRUE),3), 
                   avg_f = round(mean(statistic, na.rm = TRUE),1),
                   sd_f = round(sd(statistic, na.rm = TRUE),1),
                   num_sign_0.05 = sum(p_value < 0.05, na.rm = TRUE),
                   num_sign_0.01 = sum(p_value < 0.01, na.rm = TRUE)) %>% 
  data.frame()%>% arrange(dplyr::desc(avg_f))

write.csv(out, file = paste0(dir_determinants, "table_results_wdstest_per_y_var.csv"), 
          row.names = FALSE)




###### IDENTIFY TAXA ASSOCIATED WITH SPECIES AND DIET CATEGORIES ###############


# ========================= Run the LEFSE function =============================


# Run the LEFSE analysis ----

ls_res_lefse <- list()
for (trt in nms_factors) {
    ls_res_lefse[[trt]] <- gimme_lefse(ps_data, variable = trt, 
                                             taxrank = "Genus", 
                                             data_transform = "compositional")
}

# result tables ----

tab_taxa <- ps_data@tax_table@.Data %>% data.frame()
tab_taxa <- do.call(rbind, lapply(split(tab_taxa, tab_taxa$Genus), function(x) {
  x[1,]
}))
tab_taxa$feature <- tab_taxa$Genus


# Export the results ----

# Species

tmp <- ls_res_lefse$Species$res_stats_lda %>% 
  filter(num_sign_diff != 0 & kw_pvalues < 0.05) %>% 
  left_join(tab_taxa, by = "feature")

write.csv(tmp, file = paste0(dir_determinants, "table_lefse_Species.csv"), row.names = FALSE)


# Diet category

tmp <- ls_res_lefse$Diet_category$res_stats_lda %>% 
  filter(num_sign_diff != 0 & kw_pvalues < 0.05) %>% 
  left_join(tab_taxa, by = "feature")

write.csv(tmp, file = paste0(dir_determinants, "table_lefse_DietCategory.csv"), row.names = FALSE)




## ASSOCIATION BETWEEN DIET, MORPHOLOGICAL TRAITS, PHYLOGENY AND MICROBIOME ####


# ========================= PREPARE DATA =======================================

#     Microbiome data ----
ps_compo <- ps_data %>% microbiome::transform("compositional")

#   Morphological traits data ----
metamorpho <- metadata[, nms_traits_morpho]

#   Diet data ----
filename <- paste0(dir_data, "table_diet_composition.csv")
data_diet <- read.csv(filename, encoding = "UTF-7", row.names = 1) %>% 
  rownames_to_column("Species")
nm_food_items <- names(data_diet)[-1]

# expand to individuals
metadiet <- left_join(metadata, data_diet, "Species") %>% 
  dplyr::select(sample_id_fastq, nm_food_items) %>% 
  column_to_rownames("sample_id_fastq") %>% as.data.frame()

#   Phylogenetic data ----
# do the eigenvector decomposition of the sparidae phylogeny
# We keep the first 3 eigenvectors as they account for > 10% of the signal each 
#  and for a total of > 60%

phylogeny_pruned <- readRDS(paste0(dir_dissimilarity, "phylogeny_pruned.rds"))
tmp <- PVRdecomp(phylogeny_pruned, type = "newick")
tmp2 <- data.frame(tmp@Eigen$vectors[, 1:3], row.names = phylogeny_pruned$tip.label) %>% 
  setNames(paste0("eigen_phylo_", 1:3)) %>%
rownames_to_column("Species")

# expand to individuals
metaphylo <- left_join(metadata, tmp2, "Species") %>% 
  dplyr::select(sample_id_fastq, names(tmp2)[-1]) %>% 
  column_to_rownames("sample_id_fastq") %>% as.data.frame()

# Merge all the metadata ----
metacap <- apply(cbind(metamorpho, metadiet, metaphylo), 2, as.numeric) %>% 
                   scale() %>% as.data.frame()

nms_diss_index <- names(ls_diss)


# ========================= RUN THE CAP ANALYSIS ===============================

dir_cap <- paste0(dir_determinants, "CAP_analysis/")
dir.create(dir_cap)


# --------------- RUN THE CAP ANALYSIS ON ALL DISSIMILARITY MEASURES -----------

ls_cap <- list()

for (idx in nms_diss_index) {
  
  # identify the rank
  rk <- strsplit(idx, "_")[[1]][1]
  
  if (rk == "OTU") {
    taxa_tab <- ps_compo@otu_table %>% as.matrix()
  } else {
    taxa_tab <- aggregate_taxa(ps_compo, rk) %>% otu_table() %>% as.matrix()
  }
  
  # the dissmilarity matrix associated
  diss <- ls_diss[[idx]]
  
  # CAP
  f <- formula(paste0("diss ~ ", paste(names(metacap), collapse = "+")))
  ls_cap[[idx]] <- capscale(f, data = metacap, sqrt.dist = FALSE, comm = t(taxa_tab))
}


# ---------------------- SIGNIFICANCE OF GLOBAL CAP MODEL ----------------------

#     Test the significance of the model ----

ls_aov_cap <- lapply(ls_cap, function(x) {
    out <- data.frame(anova(x))[1,]
    names(out) <- c("df","sum_square","F_value","p_value")
    out
})


# Transform into a table ----

df_aov_cap_global <- do.call(rbind, ls_aov_cap) %>% 
  data.frame() %>% 
  rownames_to_column("method") %>% 
  separate(method, into = c("rank", "index"), sep = "_", extra = "merge", remove = FALSE) %>%   arrange(index, dplyr::desc(F_value))


# Number and proportion of significant models ?

sum(df_aov_cap_global$p_value < 0.05)
sum(df_aov_cap_global$p_value < 0.05)  / 36

signif_aov <- df_aov_cap_global$method[df_aov_cap_global$p_value < 0.05]

df_aov_cap_global %>% group_by(index) %>% 
  summarize(sum(p_value < 0.05))


# ---------------------- ANALYSE TERMS EFFECTS ---------------------------------

#     Test the significance of terms ----

ls_res_term <- lapply(ls_cap, function(X) {
  out <- data.frame(anova.cca(X, by = "terms", step = 200)) %>% 
    rownames_to_column("y_var") %>% filter(y_var != "Residual")
  names(out) <- c("y_var", "df", "sum_square","F_value","p_value")
  out
}) 

# Summarize the effect of each explaining variable ----

df_res_term <- ls_res_term[signif_aov] %>% reformat_as_df(new_var_name = "method")

df_out <- df_res_term %>% filter(p_value < 0.05) %>% group_by(y_var) %>% 
  dplyr::summarize(num_sign = n(),
                   pct_sign = round(n() / 36 * 100, 1),
                   avg_F = mean(F_value), sd_F = sd(F_value),
  ) %>% arrange(dplyr::desc(pct_sign))

write.csv(df_out, file = paste0(dir_cap, "table_results_CAP_terms_effect.csv"),
          row.names = FALSE)


# --------------------- ANALYSE SIGNIFICANCE OF AXIS ---------------------------


#     Test the significance of axes ----

ls_res_axis <- lapply(ls_cap, function(fit) {
  out <- data.frame(anova.cca(fit, by = "axis", step = 200)) %>% 
    rownames_to_column("CAP_axis") %>% filter(CAP_axis != "Residual")
  names(out) <- c("CAP_axis", "df", "sum_square","F_value","p_value")
  out
}) 


# Variance explained by the axes ----

df_explained_var <- do.call(rbind, lapply(ls_cap, function(x){
  cap_var_props(x)[1:5]})) %>% 
  data.frame() %>% rownames_to_column(var = "method") %>% 
  pivot_longer(cols = 2:6, names_to = "CAP_axis") %>% 
  dplyr::rename("explained_var" = value)

# Summarize the results ----

df_res_axis <- ls_res_axis %>% 
  reformat_as_df(new_var_name = "method") %>% 
  separate(method, into = c("rank", "index"), sep = "_", extra = "merge", remove = FALSE) %>% 
  filter(CAP_axis %in% paste0("CAP", 1:5)) %>% 
  bind_cols(dplyr::select(df_explained_var, explained_var))


# variance explained by the global models that are significant (AOV on global model)

df_res_axis %>% filter(method %in% signif_aov) %>% 
  filter(CAP_axis %in% c("CAP1", "CAP2")) %>% 
  group_by(method) %>% 
  dplyr::summarize(num_sign = sum(p_value < 0.05),
                   tot_expl_var = sum(explained_var)) %>% 
  dplyr::arrange(dplyr::desc(tot_expl_var)) %>% data.frame()



# ========= ANALYSE CAP RESULTS IN DETAILS FOR ONE SPECIFIC METHOD  ============


# Identify the method ----

nm_selected_method  <- "OTU_taxo_q1"
nm_selected_idx <- "taxo_q1"
nms_methods <- paste0(nms_ranks, "_", nm_selected_idx)


# Prepare data for the plots ----

ls_data_for_plot <- list()

for (rk in nms_ranks) {
  
  nm_mthd <- paste0(rk, "_", nm_selected_idx)
  fit   <- ls_cap[[nm_mthd]]
  terms <- ls_res_term[[nm_mthd]]
  ls_data_for_plot[[rk]]$var_expl <- df_explained_var %>% 
    filter(method == nm_mthd & CAP_axis %in% c("CAP1", "CAP2"))
  
  # extract data from the CCA object
  ls_data_for_plot[[rk]]$sple_score <- fit$CCA$wa[, c(1,2)]
  ls_data_for_plot[[rk]]$otu_score  <- fit$CCA$v[, c(1,2)]
  ls_data_for_plot[[rk]]$var_score  <- fit$CCA$biplot[, c(1,2)] * 0.75
  
  # identify significant variables
  ls_data_for_plot[[rk]]$nm_terms <- terms %>% filter(p_value < 0.05) %>% 
    dplyr::select(y_var) %>% unlist()
}



# -------- Find the taxa correlated with the significant variables -------------

ls_correlated_taxa <- list()

for (rk in nms_ranks) {
  
  X <- ls_data_for_plot[[rk]]$otu_score
  
  # identify the thresholds for extreme values
  thds <- apply(X, 2, function(x) {
    quantile(x, probs = seq(0, 1, 0.1))
  })
  # identify the extreme taxa
  ls_extreme_taxa <- lapply(1:2, function(x) {
    names(which(X[,x] < thds[2,x] | X[,x] > thds[10,x]))
  })
  nm_extreme_taxa <- ls_extreme_taxa %>% unlist() %>% unique()
  
  # get the taxa table
  if (rk == "OTU") {
    taxa_tab <- ps_compo@otu_table@.Data %>% as.matrix()
  } else {
    tmp <- aggregate_taxa(ps_compo, rk) 
    taxa_tab <- tmp@otu_table@.Data %>% as.matrix()
  }
  
  # test correlation between taxa and variables
  
  vars_to_test <- ls_data_for_plot[[rk]]$nm_terms
  
  rescor <- lapply(vars_to_test, function(var) {
    df_cor <- do.call(rbind, lapply(nm_extreme_taxa, function(x) {
      dat <- taxa_tab[x,]
      res <- cor.test(dat, metacap[, var])
      data.frame(estimate = round(res$estimate,2), p_value = round(res$p.value,3))
    })) %>% mutate(nm_taxa = nm_extreme_taxa)
    df_cor %>% filter(p_value < 0.05) %>% arrange(estimate)
  }) %>% setNames(vars_to_test)
  
  ls_correlated_taxa[[rk]] <- rescor
} 


# --------- Summarize the microbial taxa associated with dietary items --------

nms_y_vars <- lapply(ls_data_for_plot, function(x) x$nm_term) %>% 
  unlist() %>% unique() %>% sort()

# Signs of the correlations: mostly positive

df_correlated_taxa <- lapply(ls_correlated_taxa, function(X) {
  reformat_as_df(X, new_var_name = "y_var")
}) %>% reformat_as_df(new_var_name = "rank")

df_correlated_taxa %>% group_by(y_var) %>% 
  dplyr::summarize(num_neg = sum(estimate < 0),
                   num_pos = sum(estimate > 0)) %>% 
  arrange(num_pos)


# Number of associate taxa ----

tplte <- matrix(0, length(nms_y_vars), length(nms_ranks))
row.names(tplte) <- nms_y_vars
colnames(tplte) <- nms_ranks

for (rk in nms_ranks) {
  for (item in nms_y_vars) {
    repl <- nrow(ls_correlated_taxa[[rk]][[item]])
    if (! is.null(repl))  tplte[item,rk] <- repl
  }
}

tplte <- tplte[order(rowSums(tplte)),]

write.csv(tplte, file = paste0(dir_cap, 
                               "table_results_CAP_number_taxa_associated_with_vars.csv"),
          row.names = TRUE)


# Taxonomic composition of correlated OTUs -----

ls_taxa <- lapply(ls_correlated_taxa$OTU, function(X) {
  ps_data@tax_table@.Data[X$nm_taxa,] %>% as.data.frame() %>% 
    arrange(Phylum, Class, Order, Family, Genus) %>% dplyr::select(- Domain)
}) 

df_taxa <- reformat_as_df(ls_taxa, new_var_name = "y_var")
for (i in names(df_taxa)) { df_taxa[,i] <- factor(df_taxa[,i])}

ls_table_compo <- lapply(nms_ranks[1:5], function(rk) {
  out <- do.call(cbind, lapply(split(df_taxa, df_taxa$y_var), function(X) { 
    table(X[,rk]) }))
  apply(out, 2, function(x) x / sum(x))
}) %>% setNames(nms_ranks[1:5])



#### FIGURE OF CAP ANALYSIS RESULTS ############################################


# Supplementary figure with all ranks ==========================================

png(paste0(dir_cap, "plot_CAP_analysis_", nm_selected_method, "_all_ranks.png"), height = 15, width = 25,
    unit = "cm", res= 400)
par(mar = c(4,4,1,1), mgp = c(2,0.5,0), las = 1, mfrow = c(2,3), oma = c(1,1,1,1))

for (rk in nms_ranks) {
  
  X <- ls_data_for_plot[[rk]]
  
  tmp <- rbind(X$sple_score, X$otu_score, X$var_score)
  xlim <- round(setrge(tmp[,"CAP1"],10),1)
  ylim <- round(setrge(tmp[,"CAP2"]),1)
  
  # axes labels
  xlab <- paste("CAP 1 (", X$var_expl[1,"explained_var"] * 100, " %)", sep = "")
  ylab <- paste("CAP 2 (", X$var_expl[2,"explained_var"] * 100, " %)", sep = "")
  
  # par(mar = c(3,4.5,1,1), mgp = c(1.5,0.5,0), las = 1)
  plot(NA, xlab = xlab, xlim = xlim, ylim = ylim,
       ylab = ylab, main = rk)
  abline(v = 0, h = 0, lty = "dotted")
  
  # Add data points
  points(X$otu_score, pch = 22, col = "#66666670", bg = "#66666670", cex = 0.5) 
  points(X$sple_score, pch = 21, col = as.character(metadata$Color), 
         bg = as.character(metadata$Color), cex = 0.7) 
  
  # add diet items
  dat_arrows <- X$var_score[X$nm_terms,]
  nm_terms <- gsub("_", " ", X$nm_terms)
  
  if (is.null(dim(dat_arrows))) {
    arrows(x0 = rep(0, length(nm_terms)), y0 = rep(0, length(nm_terms)), 
           x1 = dat_arrows[1], y1 = dat_arrows[2],
           col = "#66666680", length = 0.0, cex = 1.5, lwd = 1)
    text(dat_arrows[1], dat_arrows[2], labels = nm_terms, font = 3, cex = 0.75, 
         adj =c(0.5,0.5))
  } else {
    arrows(x0 = rep(0, length(nm_terms)), y0 = rep(0, length(nm_terms)), 
           x1 = dat_arrows[, 1], y1 = dat_arrows[, 2],
           col = "#66666680", length = 0.0, cex = 1.5, lwd = 1)
    text(dat_arrows, labels = nm_terms, font = 3, cex = 0.75, 
         adj =c(0.5,0.5))
  }
}

dev.off()




# Figure 2 =====================================================================


png(paste0(dir_cap, "figure_CAP_analysis.png"), height = 12, width = 25,
    unit = "cm", res= 400)

layout(matrix(c(1,1,2,2,1,1,2,2,1,1,3,3,1,1,3,3), 4, 4, TRUE))
par(mar = c(3,3,0,2), mgp = c(2.5,0.5,0), las = 1, oma = c(4,2,2,2), 
    font.lab=2, cex.lab =1.2)

# -------------------- CAP analysis plot ------------------

X <- ls_data_for_plot$OTU

tmp <- rbind(X$sple_score, X$otu_score, X$var_score)
xlim <- round(setrge(tmp[,"CAP1"],10),1)
ylim <- round(setrge(tmp[,"CAP2"]),1)

# axes labels
xlab <- paste("CAP 1 (", X$var_expl[1,"explained_var"] * 100, " %)", sep = "")
ylab <- paste("CAP 2 (", X$var_expl[2,"explained_var"] * 100, " %)", sep = "")

par(xpd=NA)
plot(NA, xlab = xlab, xlim = xlim, ylim = ylim,
     ylab = ylab, main = "")
par(xpd=FALSE)
abline(v = 0, h = 0, lty = "dotted")

# Add data points
points(X$otu_score, pch = 22, col = "#66666670", bg = "#66666670", cex = 0.7) 
points(X$sple_score, pch = 21, col = as.character(metadata$Color), 
       bg = as.character(metadata$Color), cex = 1) 

# add diet items
dat_arrows <- X$var_score[X$nm_terms,]
nm_terms <- gsub("_", " ", X$nm_terms)
nm_terms[2] <- paste0("        ", nm_terms[2])

arrows(x0 = rep(0, length(nm_terms)), y0 = rep(0, length(nm_terms)), 
       x1 = dat_arrows[, 1], y1 = dat_arrows[, 2],
       col = "#66666680", length = 0.0, cex = 1.5, lwd = 1)
text(dat_arrows, labels = gsub("_", " ", nm_terms), font = 3, cex = 1.5, 
     adj =c(0.5,0.5))

legend("topleft", legend = gsub("_", " ", names(colors_sparids)), 
       pt.bg = colors_sparids, col = colors_sparids,
       pch = 22, pt.cex = 2, cex = 1, title = "", bty = "n",
       x.intersp = 1, y.intersp = 1.1, text.font = 3)

text(-.80, 0.51, labels = "(A)", font = 2, cex = 1.5, xpd = NA)

# --------------- Number of taxa associated with traits ---------------

tplte <- tplte[order(rowSums(tplte)),]

par(mar = c(3,6,0,1), mgp = c(2,0.5,0), las = 1, xpd = NA)
colvec <- c("#000000", "#333333", "#666666", "#999999", "#CCCCCC", "#FFFFFF")
barplot(t(tplte), horiz = TRUE, las = 1, 
        xlab = "Number of taxa correlated with variable",
        names.arg = c("Annelids","Echinoderms", "Mollusks","Relative gut length","Gut length"), 
        col = colvec,     xlim = c(0,60))
legend("bottomright", legend = colnames(tplte), pt.bg = colvec, pch = 22, pt.cex = 2,
       cex = 1, bty = "n")
text(-10, 6, labels = "(B)", font = 2, cex = 1.5)

# --------------- Taxonomy of associated OTU ---------------

lev <- "Family"
par(mar = c(4,6,1,1), mgp = c(2,0.5,0), las = 1, xpd = NA)
colvec2 <- distinctColorPalette(k = nrow(ls_table_compo[[lev]]), altCol = FALSE, runTsne = FALSE)
dat <- ls_table_compo[[lev]][,c(2,3,1)]
barplot(dat, horiz = TRUE, las = 1, col = colvec2,
        xlab = "Distribution of OTUs correlated with variables into bacterial Families  ",
        names.arg = c("Mollusks", "Relative gut length", "Gut length"))

legend(-0.55,-.9, legend = row.names(dat), pt.bg = colvec2, ncol = 5,
       pch = 22, pt.cex = 2, cex = 1, title = "", bty = "n",
       x.intersp = 1, y.intersp = 1.1)
text(-0.175, 3.5, labels = "(C)", font = 2, cex = 1.5)

dev.off()

