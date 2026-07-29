##===================================================================
##=====================================================================PROJECT SETUP
## Creates the standard SCAMP_GIS folder structure.
## Run once after cloning the repo.
##===================================================================

dirs <- c(
  "data/raw",
  "data/processed",
  "scripts",
  "outputs/vectors",
  "outputs/rasters",
  "figures"
)

for (d in dirs) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
    cat("Created:", d, "\n")
  } else {
    cat("Exists: ", d, "\n")
  }
}

cat("\nSCAMP_GIS folder structure ready.\n")
