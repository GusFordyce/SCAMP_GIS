# ============================================================
# SCAMP GIS Analysis — Project Setup
# ============================================================
# This script documents the one-time setup steps used to
# initialise the SCAMP_GIS project and GitHub repository.
#
# These commands were run in the console (they restart RStudio):
#   usethis::create_project("~/SCAMP_GIS")
#   usethis::use_git()
#   usethis::use_github()
# ============================================================





# --- Packages required for GIS analysis ---
install.packages(c(
  "sf",          # vector spatial data
  "terra",       # raster spatial data (memory-efficient)
  "dplyr",       # data manipulation
  "ggplot2",     # plotting
  "usethis",     # project and git management
  "gitcreds"     # GitHub credentials
))

