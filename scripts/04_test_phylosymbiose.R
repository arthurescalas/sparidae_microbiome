################################################################################
#         _______  _____  _____ ___ ____  _   _ __  __ _____ ____  
#        | ____\ \/ / _ \|  ___|_ _/ ___|| | | |  \/  | ____|  _ \ 
#        |  _|  \  / | | | |_   | |\___ \| |_| | |\/| |  _| | | | |
#        | |___ /  \ |_| |  _|  | | ___) |  _  | |  | | |___| |_| |
#        |_____/_/\_\___/|_|   |___|____/|_| |_|_|  |_|_____|____/ 
#
################################################################################
#
# R script to test relatinoships between host and microbiome dissimilarity
#
# Arthur Escalas July 2020
# arthur.escalas@gmail.com
################################################################################


################################ LOAD DATA #####################################

#     Load the CORE phyloseq objects ----

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

# microbiome dissimilarity
df_diss_µ <- read.csv(paste0(dir_dissimilarity, "table_microbiome_dissimilarity_samples.csv"))
df_diss_µ_sp <- read.csv(paste0(dir_dissimilarity, "table_microbiome_dissimilarity_species.csv"))


#   List the ranks of bacterial taxonomy ----

nms_ranks <- c("Phylum", "Class", "Order", "Family", "Genus", "OTU")

#   list the differents variables ----

nms_traits_morpho <- c("Gut_length", "Relative_gut_length")
nms_factors <- c("Species", "Diet_category")

# name of explanatory variables
nms_X_vars <- names(df_diss_fish)[c(4:7)]


# names of the methodsXrank combinations ----

nms_methods <- nms_Y_vars <- names(ls_diss)
nms_index <- c("phylo_q0", "phylo_q1", "phylo_q2", "taxo_q0", "taxo_q1", "taxo_q2")


# color vector for species pairs ----

pch_transparency <- 80

colvec <- rep(paste0("#B3B3B3", pch_transparency), nrow(df_diss_µ_sp)) 
colvec[grep("Sarpa_salpa", df_diss_µ_sp$sp_comb)] <- paste0("#228B22", pch_transparency) # green
colvec[str_count(df_diss_µ_sp$sp_comb, "Diplodus") == 2] <- paste0("#00BFFF", pch_transparency) # blue
colvec[str_count(df_diss_µ_sp$sp_comb, "Pagellus") == 2] <- paste0("#FF0777", pch_transparency) # pink
colvec[grep("Boops_boops", df_diss_µ_sp$sp_comb)] <- paste0("#FF5500", pch_transparency) # green
colvec[grep("Boops_boops__Sarpa_salpa", df_diss_µ_sp$sp_comb)] <- paste0("#000000", pch_transparency) # gold



################# TEST THE RELATIONHSIP BETWEEN DISSMILARITIES #################


#### 1. BETWEEN MICROBIOME DISSMILARITY AND FISH PHYLOGENETIC DISTANCE


# --------------------------- Mantel tests ----------------------------------

ls_mant <- list()

for (xvar in nms_X_vars) {
  for (yvar in nms_Y_vars) {
    
    diss <- df_diss_µ_sp[, yvar]
    fit <- mantel(formula = formula(diss ~ df_diss_fish[, xvar]), 
                  data = sys.parent(), nperm = 1000,
                  mrank = FALSE, nboot = 500, pboot = 0.9, cboot = 0.95)
    ls_mant[[xvar]][[yvar]] <- data.frame(t(fit))
  }
}

df_mant <- lapply(ls_mant, function(X) {
  reformat_as_df(X, new_var_name = "data_type")
}) %>% reformat_as_df(new_var_name = 'X_var')
df_mant <- cbind(df_mant, str_split_fixed(df_mant$data_type, "_", 2) %>% 
                   data.frame() %>% setNames(c("rank", "index")))


# --------------------------- explore the results ------------------------------

mask <- df_mant$pval1 < 0.05
df_res_mantel_signif <- df_mant %>% filter(mask) %>% arrange(pval1)

df_mant %>% group_by(X_var) %>% 
  dplyr::summarize(avg_corr = round(mean(mantelr),2),
                   sd_corr  = round(sd(mantelr), 2),
                   avg_pval = round(mean(pval1),2),
                   num_sign = sum(pval1 < 0.05)) %>% data.frame() %>% 
  arrange(dplyr::desc(num_sign), avg_pval)


# --------------------------- plot the relationships ---------------------------

dir_plot <- paste0(dir_phylosymbiose, "plots_relationship_dissimilarity_µbiome_hosts/")
dir.create(dir_plot)

# MANTEL plots ----

nms_for_plot <- c("Host phylogeny", "Gut length", "Relative gut length",
                  "Diet composition") %>% setNames(nms_X_vars)

nms_hn_idx <- c("Phylogenetic dissimilariry, q = 0", 
                "Abundance-weighted phylogenetic dissimilarity, q = 1",
                "Abundance-weighted phylogenetic dissimilarity, q = 2",
                "Taxonomic dissimilarity, q = 0",
                "Abundance-weighted taxonomic dissimilarity, q = 1",
                "Abundance-weighted taxonomic dissimilarity, q = 2"
                ) %>% setNames(nms_index)

for (idx in nms_index) {
  
  png(paste0(dir_plot, "plot_MANTEL_all_X_vars_", idx,".png"),
      height = 25, width = 25, res = 500, unit = "cm")
  
  layout(matrix(1:24, 6, 4, FALSE))
  par(mar = c(2,2,2,1), oma = c(6,7,3,4), las = 1, 
      mgp = c(1.5,0.5,0), xpd = TRUE)
  
  for (xvar in nms_X_vars) {
    for (rk in nms_ranks) {
      
      dat <- df_mant %>% filter(index == idx & X_var == xvar & rank == rk)
      yvar <- paste0(rk, "_", idx)
      
      if (dat$pval1 < 0.05) {
        p_text <- "*"
        p_font <- 2
        p_cex <- 1.3
      } else {
        p_text <- ""
        p_font <- 1
        p_cex <- 1
      }
      
      ylim <- c(0, setrge(df_diss_µ_sp[, yvar], 5)[2])
      xlim <- c(0, setrge(df_diss_fish[, xvar], 5)[2])
      plot(df_diss_µ_sp[, yvar] ~ df_diss_fish[, xvar], ylab = "", xlab = "",
           main = "", cex = 1, pch = 21, bg = colvec, col = colvec, 
           xlim = xlim, ylim = ylim)
      
      if (dat$rank == "Phylum") {
        mtext(text = nms_for_plot[xvar], side = 3, line = 2, font = 2)
      }
      if (dat$X_var == "Diet_FishBase") {
        mtext(text = dat$rank, side = 4, line = 2, las = 0, 
              font = 4, col = "black")
      }
      
      # results of mantel test
      text(x = min(xlim), y = max(ylim) + 0.1*diff(ylim), font = p_font, 
           cex = p_cex,
           labels = paste0("r = ", round(dat$mantelr, 2)), adj = c(0,0))
      text(x = min(xlim) + 0.4*diff(xlim), y = max(ylim) + 0.1*diff(ylim), 
           font = p_font, cex = p_cex,
           labels = paste0("p-val. = ", round(dat$pval1, 3), p_text), adj = c(0,0))
    }
  }
  mtext(text = paste0("Host dissimilarity"), side = 1, line = 3, 
        las = 0, font = 2, outer = TRUE, cex = 1.5)
  mtext(text = "Microbiome dissimilarity", side = 2, line = 4, 
        las = 0, font = 2, outer = TRUE, cex = 1.5)
  mtext(text = nms_hn_idx[idx], side = 2, 
        line = 2, las = 0, font = 2, outer = TRUE, cex = 1)
  
  dev.off()
}



#### 2. BETWEEN MICROBIOME DISSMILARITY PER TAXA AND FISH PHYLOGENETIC DISTANCE 


# Load the dissimilarity data

ls_diss <- readRDS(paste0(dir_dissimilarity, "list_microbiome_dissimilarity_matrices_per_taxa.rds.rds"))
df_diss <- read.csv(paste0(dir_dissimilarity, "table_microbiome_dissimilarity_per_taxa.csv")) %>% 
  filter(index == "phylo_q2")

# Objects for analyses

ranks  <- unique(df_diss$rank) %>% as.vector()
Y_vars <- unique(df_diss$index) %>% as.vector()
X_vars <- nms_X_vars


# ============================= Using Mantel test ==============================

#   Run Mantel test on each combination of rank/method/taxa ----

ls_mant <- list()
ls_data_diss <- list()

for (xvar in X_vars) {   # loop on explanatory variables
  
  for (index in Y_vars) {  # loop onmethods
    
    for (rk in ranks) {   # loop on ranks
      
      df <- df_diss %>% filter(index == index & rank == rk)
      nms_taxa <- unique(df$taxa_name)
      
      for (taxa in nms_taxa) {  # loop on taxa
        
        # Get the dissimilarity matrices for the focus taxa
        
        d_µ <- ls_diss[[rk]][[taxa]][[index]]
        d_µ[d_µ == "Inf"] <- 1
        d_µ[is.na(d_µ)] <- 0
        
        if (length(d_µ) > 25) { # test whether the matrix has more than 5 species
          
          # Recover the phylogenetic distance for the species pairs for which 
          #  we have a microbiome dissmilarity
          
          # Make an empty matrix with the right dimensions
          tmp <- as.matrix(d_µ)
          tplte <- matrix(NA, nrow(tmp), ncol(tmp))
          row.names(tplte) <- colnames(tplte) <- row.names(tmp)
          # fill up the matrix
          for (i in row.names(tmp)) {
            for (j in row.names(tmp)) {
              if (i != j) {
                sp_pair <- paste0(sort(c(i,j)), collapse = "__")
                repl <- df_diss_fish %>% filter(sp_comb == sp_pair) %>% 
                  dplyr::select(all_of(xvar)) %>% unlist()
                tplte[i,j] <- repl
              }
            }
          }
          d_X <- as.dist(tplte)
          
          # Fit the mantel test ----
          
          fit <- mantel(formula = d_µ ~ d_X, 
                        data = sys.parent(), nperm = 1000,
                        mrank = FALSE, nboot = 500, pboot = 0.9, cboot = 0.95)
          
          # Make a table of dissmilarity for later plot ----
          
          tmp <- dendextend::dist_long(d_X) %>% dplyr::select(distance) %>% 
            rename(distance = "fd")
          df_out <- cbind(dendextend::dist_long(d_µ), tmp)
          df_out$sp_comb <- sapply(1:nrow(df_out), function(x) { 
            tmp <- sort(c(as.character(df_out$rows)[x], 
                          as.character(df_out$cols)[x]))
            paste0(tmp[1], "__", tmp[2])
          })

          
          # output the results ----
          ls_mant[[xvar]][[index]][[rk]][[taxa]] <- data.frame(t(fit)) 
          ls_data_diss[[xvar]][[index]][[rk]][[taxa]] <- df_out
          
        }##eo if matrix size
      }
    }
  }
}


#   Reformat results into a table ----

df_res_mantel <- lapply(ls_mant, function(X) {
  lapply(X, function(XX) {
    lapply(XX, function(XXX) { 
      reformat_as_df(XXX, new_var_name = "taxa_name")  
    }) %>%  reformat_as_df(new_var_name = "rank")
  }) %>% reformat_as_df(new_var_name = "index")
}) %>% reformat_as_df(new_var_name = "X_var")


# ---------------------------  Explore the results -----------------------------

# only significant results

df_res_mantel_signif <- df_res_mantel %>% 
  filter(pval1 < 0.05 | pval2 < 0.05) %>% arrange(mantelr)

# Number of significant relationhips per trait/phylo/diet
df_res_mantel_signif %>% group_by(X_var) %>% tally() %>% arrange(n) %>% data.frame()


# Plot the relatinships for each type of data -----

df_plot <- df_res_mantel_signif %>% arrange(rank, taxa_name) %>% data.frame()

nms_for_plot <- c("Host phylogeny", "Gut length", "Relative gut length",
                  "Diet composition") %>% setNames(nms_X_vars)

png(paste0(dir_plot, "plot_MANTEL_per_taxa.png"),
    height = 25, width = 25, res = 600, unit = "cm")

layout(matrix(1:16, 4, 4, TRUE))
par(mar = c(4,4,2,1), oma = c(6,7,3,2), las = 1, 
    mgp = c(2.3,0.7,0), xpd = TRUE, font.lab = 2)

lapply(1:nrow(df_plot), function(X){
  
  dat <- df_plot[X,]
  tmp <- ls_data_diss[[dat$X_var]][[dat$index]][[dat$rank]][[dat$taxa_name]]
  
  pch_transparency <- 80
  colvec <- rep(paste0("#B3B3B3", pch_transparency), nrow(tmp)) 
  colvec[grep("Sarpa_salpa", tmp$sp_comb)] <- paste0("#228B22", pch_transparency) # green
  colvec[str_count(tmp$sp_comb, "Diplodus") == 2] <- paste0("#00BFFF", pch_transparency) # blue
  colvec[str_count(tmp$sp_comb, "Pagellus") == 2] <- paste0("#FF0777", pch_transparency) # pink
  colvec[grep("Boops_boops", tmp$sp_comb)] <- paste0("#FF5500", pch_transparency) # green
  colvec[grep("Boops_boops__Sarpa_salpa", tmp$sp_comb)] <- paste0("#000000", pch_transparency) # gold
  
  ylim <- c(0, setrge(tmp$distance, 5)[2])
  xlim <- c(0, setrge(tmp$fd, 5)[2])
  plot(tmp$distance ~ tmp$fd, col = colvec, bg = colvec,
       cex = 1, main = "", xlim = xlim, ylim = ylim,
       pch = 21, xlab = paste0(nms_for_plot[dat$X_var]), 
       ylab = paste0(dat$rank, ": ", dat$taxa_name))

  p_text <- "*"
  p_font <- 2
  p_cex <- 1
  if (dat$mantelr < 0) {
    pval <- dat$pval2
  } else {
    pval <- dat$pval1
  }
  text(x = min(xlim), y = max(ylim) + 0.1*diff(ylim), font = p_font, 
        cex = p_cex, 
       labels = paste0("r = ", round(dat$mantelr, 2)), adj = c(0,0))
  text(x = min(xlim) + 0.4*diff(xlim), y = max(ylim) + 0.1*diff(ylim), 
       font = p_font, cex = p_cex,
       labels = paste0("p-val. = ", round(pval, 3), p_text), adj = c(0,0))
})

mtext(text = paste0("Host dissimilarity"), side = 1, line = 3, 
      las = 0, font = 2, outer = TRUE, cex = 1.5)
mtext(text = "Microbiome dissimilarity", side = 2, line = 4, 
      las = 0, font = 2, outer = TRUE, cex = 1.5)
mtext(text = "Abundance-weighted phylogenetic dissimilarity, q = 2", side = 2, 
      line = 2, las = 0, font = 2, outer = TRUE, cex = 1)

dev.off()

