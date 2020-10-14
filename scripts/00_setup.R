################################################################################
#         _______  _____  _____ ___ ____  _   _ __  __ _____ ____  
#        | ____\ \/ / _ \|  ___|_ _/ ___|| | | |  \/  | ____|  _ \ 
#        |  _|  \  / | | | |_   | |\___ \| |_| | |\/| |  _| | | | |
#        | |___ /  \ |_| |  _|  | | ___) |  _  | |  | | |___| |_| |
#        |_____/_/\_\___/|_|   |___|____/|_| |_|_|  |_|_____|____/ 
#
################################################################################
#
# R script to set up the project
#
# Arthur Escalas June 2020
# arthur.escalas@gmail.com
################################################################################


# cleaning memory ==============================================================

rm(list = ls()) 

# Define the directories of the project ========================================
# create objects corresponding to the directories and create directories locally

dir_project <- getwd()

dir_scripts  <- paste0(dir_project, "/scripts/")
dir_data     <- paste0(dir_project, "/data/")
dir_analyses <- paste0(dir_project, "/analyses/")
dir_core_microbiome  <- paste0(dir_analyses, "01_core_microbiome/")
dir_dissimilarity    <- paste0(dir_analyses, "02_estimate_dissimilarity/")
dir_determinants     <- paste0(dir_analyses, "03_µbiome_determinants/")
dir_phylosymbiose    <- paste0(dir_analyses, "04_test_phylosymbiose/")
dir_final_plots      <- paste0(dir_analyses, "05_final_plots/")

dir.create(dir_scripts)
dir.create(dir_data)
dir.create(dir_analyses)
dir.create(dir_core_microbiome)
dir.create(dir_dissimilarity)
dir.create(dir_determinants)
dir.create(dir_phylosymbiose)
dir.create(dir_final_plots)

# Soucre functions and packages ------------------------------------------------

# source utility functions
dir_funs <- paste0(dir_scripts, "/functions/")
for (f in list.files(dir_funs, full.names = T)) { source(f) }

# source package list
source(paste0(dir_scripts, "packages.R"))
install_if_not_there(cran_packages, type = "CRAN")
install_if_not_there(bioc_packages, type = "bioconductor")
sapply(c(cran_packages, bioc_packages), require, character.only = TRUE)


# objects for data analysis ----------------------------------------------------

names_sparids <- c("Boops_boops", "Diplodus_annularis", "Diplodus_puntazzo",
                   "Diplodus_sargus", "Diplodus_vulgaris", 
                   "Lithognathus_mormyrus", "Oblada_melanura", 
                   "Pagellus_acarne", "Pagellus_erythrinus", 
                   "Pagrus_pagrus", "Sarpa_salpa", "Sparus_aurata")

colors_sparids <- c("#EE7600","#97FFFF","#00CED1","#00BFFF","#104E8B","#7F7F7F",
                    "#000000","#FF3030","#8B1A1A","#EE1289","#228B22","#FFD700")
names(colors_sparids) <- names_sparids

