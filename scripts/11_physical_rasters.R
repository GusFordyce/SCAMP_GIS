# ==============================================================================
# 11_physical_rasters.R
# Write the physical covariates that currently exist only as marine-cube bands
# out as standalone native rasters, so each can be opened on its own.
#
# No recomputation: every layer is lifted from SCAMP_marine_cube_286m.tif
# exactly as it was written there, so the values are identical to the cube.
#
# Also copies the two unmodified third-party base layers (GEBCO bathymetry and
# the COP30 DEM) into outputs/physical/ under clearer names. Both are freely
# redistributable with attribution, which SOURCES_AND_LICENCES carries.
#
# Run from the project root after 08. Takes under a minute.
# ==============================================================================

library(terra)

CUBE <- "outputs/design/SCAMP_marine_cube_286m.tif"
stopifnot("marine cube not found -- run 08 first" = file.exists(CUBE))

cube <- rast(CUBE)
stopifnot("cube is not float -- re-run 08 with the datatype fix" =
            grepl("^FLT", terra::datatype(cube)[1]))

dir.create("outputs/physical", recursive = TRUE, showWarnings = FALSE)

## Band in the cube -> file written. Only layers that do not already exist as a
## standalone native raster are listed; SPM and the salinity layers are already
## written by 04_physical.R at their own resolutions and are not duplicated.
EXPORT <- c(
  depth          = "outputs/physical/phys_bathymetry.tif",
  current_energy = "outputs/physical/phys_current_energy.tif",
  bottomt_mean   = "outputs/physical/phys_bottom_temp_mean.tif",
  bottomt_max    = "outputs/physical/phys_bottom_temp_max.tif",
  bottomt_range  = "outputs/physical/phys_bottom_temp_range.tif",
  sst_mean       = "outputs/physical/phys_surface_temp_mean.tif",
  sst_max        = "outputs/physical/phys_surface_temp_max.tif",
  sst_range      = "outputs/physical/phys_surface_temp_range.tif"
)

missing <- setdiff(names(EXPORT), names(cube))
if (length(missing) > 0) {
  stop("band(s) not in the cube: ", paste(missing, collapse = ", "),
       "\nCube bands are: ", paste(names(cube), collapse = ", "))
}

cat("\n-- Writing physical rasters from the marine cube --\n")
for (nm in names(EXPORT)) {
  r <- cube[[nm]]
  names(r) <- nm
  writeRaster(r, EXPORT[[nm]], filetype = "COG", overwrite = TRUE,
              gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3"), names = nm)
  rng <- as.numeric(global(r, range, na.rm = TRUE))
  cat(sprintf("  %-22s %-44s  %.2f to %.2f\n", nm, basename(EXPORT[[nm]]),
              rng[1], rng[2]))
}

## ---- third-party base layers, copied under clearer names -------------------
## Unmodified clips. Shipped because depth and terrain are the layers people
## need to orient themselves, and both licences permit redistribution with
## attribution.
BASE <- c(
  "outputs/physical/gebco_bathy_solway.tif" = "outputs/physical/phys_bathymetry_gebco_native.tif",
  "outputs/physical/cop30_dem_solway.tif"   = "outputs/physical/phys_elevation_cop30.tif"
)
cat("\n-- Copying base terrain layers --\n")
for (src in names(BASE)) {
  if (!file.exists(src)) {
    cat(sprintf("  MISSING: %s -- skipped\n", src))
    next
  }
  file.copy(src, BASE[[src]], overwrite = TRUE)
  cat(sprintf("  %s -> %s\n", basename(src), basename(BASE[[src]])))
}

cat("\nDone. outputs/physical/ now holds depth, currents, temperature,\n")
cat("turbidity and salinity as standalone rasters.\n")
cat("\nNOT recoverable without re-running 04_physical.R with CMEMS downloads:\n")
cat("  chlorophyll-a, dissolved oxygen, bed shear stress.\n")
cat("These were computed in memory as habitat-model filters and never written.\n")
