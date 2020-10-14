################################################################################
#  _______  _____  _____ ___ ____  _   _ __  __ _____ ____  
# | ____\ \/ / _ \|  ___|_ _/ ___|| | | |  \/  | ____|  _ \ 
# |  _|  \  / | | | |_   | |\___ \| |_| | |\/| |  _| | | | |
# | |___ /  \ |_| |  _|  | | ___) |  _  | |  | | |___| |_| |
# |_____/_/\_\___/|_|   |___|____/|_| |_|_|  |_|_____|____/ 
#
################################################################################
#
# Script to make the fiunal figures
#
# Arthur Escalas July 2019
# arthur.escalas@gmail.com
################################################################################



#### FIGURE OF FISH PHYLOGENY AND TRAITS #######################################

#     Load the CORE phyloseq objects ----

ps_data <- readRDS(paste0(dir_data, "phyloseq_object_core.rds"))


# Load the metadata ----

metadata_sp <- read.csv(paste0(dir_data, "table_metadata_species.csv"), row.names = 1)
metadata <- sample_data(ps_data) %>% data.frame()

# extract metadata

metadata_sp$Diet_category <- gsub("Selective_plankton_feeding",
                                 "Plankton feeder", metadata_sp$Diet_category)
metadata_sp$Diet_category <- gsub("Hunting_macrofauna",
                                 "Macrofauna hunter", metadata_sp$Diet_category)
metadata_sp$Diet_category <- gsub("Grazing_aquatic_plants",
                                 "Plants grazer", metadata_sp$Diet_category)

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

metadata_sp$color_diet <- metadata_sp$Diet_category
metadata_sp$color_diet <- gsub("Plankton feeder", "#0000CD", metadata_sp$color_diet)
metadata_sp$color_diet <- gsub("Macrofauna hunter", "#FF3030", metadata_sp$color_diet)
metadata_sp$color_diet <- gsub("Plants grazer", "#228B22", metadata_sp$color_diet)

n_samples <- metadata %>% group_by(Species) %>% tally() %>% 
  column_to_rownames("Species")

# read the phylogenetic tree ----
tree_fish <- readRDS(paste0(dir_dissimilarity, "phylogeny_pruned.rds"))


# reorder metadata and clean tip labesl

metadata_sp <- metadata_sp[tree_fish$tip.label, ]
n_samples <- n_samples[tree_fish$tip.label, ]
tree_fish$tip.label <- gsub("_", " ", tree_fish$tip.label)



# all in one plot ----

data_food <- read.csv(paste0(dir_data, "table_diet_composition.csv"), row.names = 1)

cols_food_III <- c(Finfish = "firebrick2",
                   Annelids = "deeppink3",
                   Echinoderms = "darkorchid",
                   Benthic_algae = "forestgreen",
                   Cnidarians = "deepskyblue",
                   Planktonic_invertebrates = "deepskyblue3",
                   Sessile_invertebrates = "deepskyblue4", 
                   Mollusks = "gold",
                   Shelled_mollusks = "gold3",
                   Planktonic_crustaceans = "darkorange",
                   Small_benthic_crustaceans = "darkorange2",
                   Large_benthic_crustaceans = "darkorange4")
X <- data_food[row.names(metadata_sp), names(cols_food_III)]


png(paste0(dir_final_plots, "figure_1_all_in_one.png"), unit = "cm",
    height = 15, width = 25, res=600)

layout(matrix(c(1,1,1,2,2,
                1,1,1,2,2,
                1,1,1,2,2), 3,5, TRUE))
# the tree ----
par(mar = c(4,0,1,0), mgp = c(2.5,1,0), xpd = NA, oma = c(2,1,3,1))
plot(tree_fish, lwd = 4, x.lim = c(0,190), cex = 1)
axis(1, at = c(0,20,40,60), labels = c(60,40,20,0), lwd = 1, cex.axis = 1.2)
text(x = 30, y = -0.75, labels = "Millions of years", cex = 1.2, font = 2)
rect(60,0.9, 61.5, 12.1, col = "white", border = "white")

# the traits ----
text(rep(145,12), 1:12, labels = n_samples, font = 1, cex = 1.2)
text(rep(170,12), 1:12, labels = metadata_sp$Diet_category, 
     col = "black", font = 2)
# the traits names
text(c(145, 170), rep(13, 2), font = 2, cex = 1.2, 
      labels = c("N", "Diet\nclassification"), adj = c(0.5,0.5))

# plot the diet ----
par(mar = c(3,0,0,12), mgp = c(2.5,1,0), las = 1, xpd = NA)

barplot(t(X), horiz = TRUE, col = cols_food_III, names.arg = rep("",12), main = "Diet composition",
        xlab = "Proportion in the diet", cex.axis = 1.2, cex.lab = 1.2, xpd = NA,
        font.lab = 2)
legend(1.05,14, legend = gsub("_", " ", names(cols_food_III)), 
       pt.bg = cols_food_III, pt.cex = 3,pch=22,
       col = cols_food_III, bty = "n", cex = 1, y.intersp = 1.5)

dev.off()



