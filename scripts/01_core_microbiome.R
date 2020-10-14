################################################################################
#         _______  _____  _____ ___ ____  _   _ __  __ _____ ____  
#        | ____\ \/ / _ \|  ___|_ _/ ___|| | | |  \/  | ____|  _ \ 
#        |  _|  \  / | | | |_   | |\___ \| |_| | |\/| |  _| | | | |
#        | |___ /  \ |_| |  _|  | | ___) |  _  | |  | | |___| |_| |
#        |_____/_/\_\___/|_|   |___|____/|_| |_|_|  |_|_____|____/ 
#
################################################################################
#
# R script to describe the core microbiome of Sparidae
#
# Arthur Escalas July 2020
# arthur.escalas@gmail.com
################################################################################


# LOAD DATA ====================================================================

# phyloseq object ----

ps_data <- readRDS(paste0(dir_data, "phyloseq_object_core.rds"))


# extract metadata ---

metadata <- sample_data(ps_data) %>% data.frame()


#### DESCRIBE THE CORE MICROBIOME OF SPARIDS #####################################


#### Different numbers taxon for various levels of taxonomic resolution ========

apply(ps_data@tax_table, 2, function(x) length(unique(x[! is.na(x)])))

# Number of taxa in each Phylum ----

tmp <- ps_data@tax_table@.Data %>% data.frame()
tab_taxa_per_phyl <- do.call(rbind, lapply(split(tmp, tmp[, "Phylum"]), function(X) {
  apply(X, 2, function(x) length(unique(x[! is.na(x)])))
}))
tab_taxa_per_phyl


# Most abundant taxa ----

nms_rks <- rank_names(ps_data)[c(2:5)]

ls_top_taxa <- lapply(nms_rks, function(rk) {
  ps_tmp <- aggregate_top_taxa(ps_data, top = 20, level = rk)
  out <- rowSums(ps_tmp@otu_table@.Data) / sum(ps_tmp@otu_table@.Data) *100
  out <- sort(out, decreasing = TRUE)
  out[!names(out) %in% c("Other", "Unknown")][1:10]
}) %>% setNames(nms_rks)

png(paste0(dir_core_microbiome, "barplot_composition_core_top10_global.png"),
    height = 15, width = 20, unit = "cm", res = 400)
par(mfrow = c(2,2), mar = c(4,10,1,1), oma = c(1,1,1,1), mgp = c(2,0.5,0), las = 1)
lapply(nms_rks, function(rk) {
  barplot(rev(ls_top_taxa[[rk]]), horiz = TRUE, main = rk)
})
mtext(text = "Relative abundance in the core microbiome (%)", outer = TRUE, 
      side = 1, line = -1)
dev.off()


#### Composition of the core in each species ===================================

# merge data per species and make data compositional

ps_sp <- merge_samples(ps_data, group = "Species")%>%
  microbiome::transform(transform = "compositional")


# Make the table at the Genus level ----

dat <- ps_sp %>% aggregate_top_taxa(top = 20, level = "Genus")
tab <- dat@otu_table * 100

# keep only the genus representing more than 5% of the abundance in at least one species

mask <- apply(tab[,names_sparids], 1, function(x) sum(x > 5)) != 0
tab <- tab[mask,]

write.csv(tab, row.names = TRUE, 
          file = paste0(dir_core_microbiome, "table_top_Genus.csv"))


# Make the table at the family level ----

dat <- ps_sp %>% aggregate_top_taxa(top = 20, level = "Family")
tab <- dat@otu_table * 100


# keep only the family representing more than 10% of the abundance in at least one species

mask <- apply(tab[,names_sparids], 1, function(x) sum(x > 10)) != 0
tab <- tab[mask,]

write.csv(tab, row.names = TRUE, 
          file = paste0(dir_core_microbiome, "table_top_Family.csv"))
































