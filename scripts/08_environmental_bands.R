# 08_environmental_bands.R
#
# PURPOSE
# Restoration sites are not yet confirmed, so named control sites cannot be
# picked now. This script bands the whole programme area into environmental
# strata in advance so that, once a site is chosen, control matching becomes a
# lookup: a candidate control is a cell in the same band, ranked by covariate
# distance (match_controls()), with balance checked after the fact
# (match_balance()). The same aligned covariate stacks it builds are written
# out as the two partner-handover data cubes (SCAMP_marine_cube_286m.tif,
# SCAMP_terrestrial_cube_100m.tif).
#
# Run top-to-bottom with source(). Two passes:
#   Pass 1 — leave K_MARINE / K_TERR as NA. The script fits k = 2:12, writes
#            the elbow plots, and stops before band assignment.
#   Pass 2 — set K_MARINE / K_TERR below from the elbow plots, re-run.
#
# AUDIT NOTES (see conversation / PROJECT_INVENTORY.md for full detail) —
# three deliberate departures from a literal reading of the brief, all
# editable in the COVARIATE TABLE below:
#   1. scamp_aoi.gpkg carries no marine/land attribute (fid, geom only), so
#      the marine and land programme-area masks are rebuilt here from the same
#      raw boundary files and filters 02_themeD_species.R:47-75 uses, rather
#      than guessed.
#   2. "Dynamic Coast modern MHWS erosion/accretion class" does not exist as
#      named: DC2_MHWS_MODERN.dbf carries only survey/capture metadata, no
#      rate or class field. The actual historical rate (HIST_RATE, m/yr) is on
#      the separate DC2_RCP8_Transects layer. That is what is wired in below,
#      binned at +-DC_STABLE_THRESHOLD_MYR. Confirmed.
#   3. "Sea temperature mean, max and annual range from thetao and bottomT" is
#      read as six bands (surface x3, bottom x3), not three. All six go in the
#      marine cube, but only the three bottomT layers are used for banding:
#      surface and bottom temperature are highly correlated in the shallow
#      Solway, so banding on both would double-weight temperature.
# Coverage percentages in the brief could not be measured without running R;
# outputs/design/coverage_report.csv carries the real figures after pass 1.

suppressMessages({
  library(sf)
  library(terra)
  library(dplyr)
  library(ggplot2)
})

set.seed(1)
terraOptions(memfrac = 0.6)

dir.create("outputs/design", recursive = TRUE, showWarnings = FALSE)
dir.create("data/raw/climate", recursive = TRUE, showWarnings = FALSE)

cat("========== 08_environmental_bands.R ==========\n\n")

##=====================================================================
## 0. CONFIG CONSTANTS
##=====================================================================

K_MARINE <- 7   # chosen from the elbow plot in outputs/design/ (k_selection_marine.png); NA re-runs the k-selection screen only
K_TERR   <- 7   # chosen from the elbow plot in outputs/design/ (k_selection_terrestrial.png); NA re-runs the k-selection screen only

SAMPLE_N           <- 50000   # cells sampled per cube to fit scaling + kmeans
CAT_WEIGHT         <- 0.5     # weight applied to one-hot categorical columns
SEED                <- 1
RARE_LEVEL_MIN_PROP <- 0.02   # categorical levels below this sample share -> "other"
DC_STABLE_THRESHOLD_MYR <- 0.1  # Dynamic Coast HIST_RATE band: +-this = "stable"
MARINE_GAPFILL_MAX_PASSES <- 50
CMEMS_CROP_BUFFER_M <- 20000  # crop CMEMS to marine AOI + this buffer (grid is ~9 km)
REPR_SAMPLE_N      <- 20000  # cells sampled per AOI / per region for representativeness

set.seed(SEED)

##=====================================================================
## 1. COVARIATE TABLE — the single editable source of truth
##=====================================================================
## type: continuous | categorical
## resampling: bilinear | mode  (near-NN under the hood via terra for mode)
## cube: marine | terrestrial
## gap_fill: TRUE only for raw CMEMS-derived continuous marine layers that are
##           known to have nearshore/inlet gaps at source resolution. Finished
##           analytical products (cumulative pressure, pollution exposure) and
##           all categorical / terrestrial layers are never gap-filled.

covariate_table <- tibble::tribble(
  ~name,                     ~path,                                                              ~layer,                           ~field,        ~type,         ~resampling, ~cube,          ~use_for_banding, ~gap_fill, ~units,              ~description,
  "depth",                   "outputs/physical/gebco_bathy_solway.tif",                           NA_character_,                    NA_character_, "continuous",  "bilinear",  "marine",       TRUE,             FALSE,     "m",                 "GEBCO-derived bathymetry, 286 m marine template",
  "current_energy",          "outputs/physical/ke_currents_solway.tif",                           NA_character_,                    NA_character_, "continuous",  "bilinear",  "marine",       TRUE,             TRUE,      "index",             "Kinetic energy of tidal/residual currents",
  "turbidity_spm",           "outputs/physical/spm_mean_solway.tif",                              NA_character_,                    NA_character_, "continuous",  "bilinear",  "marine",       TRUE,             TRUE,      "g/m3",              "Mean suspended particulate matter (ocean colour)",
  "salinity_p90",            "outputs/physical/salinity_p90_mean_solway.tif",                     NA_character_,                    NA_character_, "continuous",  "bilinear",  "marine",       TRUE,             TRUE,      "PSU",               "90th percentile salinity",
  "salinity_p10",            "outputs/physical/sal_p10_solway.tif",                               NA_character_,                    NA_character_, "continuous",  "bilinear",  "marine",       TRUE,             TRUE,      "PSU",               "10th percentile salinity",
  "salinity_days_below20",   "outputs/physical/sal_mean_days_below20_solway.tif",                 NA_character_,                    NA_character_, "continuous",  "bilinear",  "marine",       TRUE,             TRUE,      "days/yr",           "Mean days per year below 20 PSU",
  "energy_class",            "outputs/physical/energy_class_solway.gpkg",                         "energy_class_solway",           "Energy",      "categorical", "mode",      "marine",       TRUE,             FALSE,     "class",             "EUSeaMap energy class ('No energy information' treated as NA)",
  "seabed_substrate",        "outputs/physical/ukseamap_seabed_substrate_solway.gpkg",            "ukseamap_seabed_substrate_solway","sub_class",  "categorical", "mode",      "marine",       TRUE,             FALSE,     "class",             "UKSeaMap seabed substrate ('unknown substrate' treated as NA; rare levels collapsed)",
  "cumulative_pressure",     "outputs/threats/cumulative_pressure_with_pollution.tif",            NA_character_,                    NA_character_, "continuous",  "bilinear",  "marine",       TRUE,             FALSE,     "index",             "Cumulative marine pressure incl. diffuse pollution",
  "pollution_exposure",      "outputs/catchments/marine_pollution_exposure.tif",                  NA_character_,                    NA_character_, "continuous",  "bilinear",  "marine",       TRUE,             FALSE,     "index",             "Sea-path decay exposure from catchment pour points",
  "sst_mean",                "data/processed/physical/thetao",                                    NA_character_,                    NA_character_, "continuous",  "bilinear",  "marine",       FALSE,            TRUE,      "degC",              "Surface temperature, long-term mean (CMEMS thetao, derived). Cube-only: not used for banding (correlated with bottomT)",
  "sst_max",                 "data/processed/physical/thetao",                                    NA_character_,                    NA_character_, "continuous",  "bilinear",  "marine",       FALSE,            TRUE,      "degC",              "Surface temperature, mean of annual maxima (derived). Cube-only: not used for banding (correlated with bottomT)",
  "sst_range",               "data/processed/physical/thetao",                                    NA_character_,                    NA_character_, "continuous",  "bilinear",  "marine",       FALSE,            TRUE,      "degC",              "Surface temperature, annual range: sst_max - winter min (derived). Cube-only: not used for banding (correlated with bottomT)",
  "bottomt_mean",            "data/processed/physical/bottomT",                                   NA_character_,                    NA_character_, "continuous",  "bilinear",  "marine",       TRUE,             TRUE,      "degC",              "Bottom temperature, long-term mean (CMEMS bottomT, derived)",
  "bottomt_max",             "data/processed/physical/bottomT",                                   NA_character_,                    NA_character_, "continuous",  "bilinear",  "marine",       TRUE,             TRUE,      "degC",              "Bottom temperature, mean of annual maxima (derived)",
  "bottomt_range",           "data/processed/physical/bottomT",                                   NA_character_,                    NA_character_, "continuous",  "bilinear",  "marine",       TRUE,             TRUE,      "degC",              "Bottom temperature, annual range: max - min (derived)",

  "elevation",               "outputs/physical/cop30_dem_solway.tif",                             NA_character_,                    NA_character_, "continuous",  "bilinear",  "terrestrial",  TRUE,             FALSE,     "m",                 "COP30 DEM",
  "slope",                   "outputs/catchments/slope_cop30.tif",                                NA_character_,                    NA_character_, "continuous",  "bilinear",  "terrestrial",  TRUE,             FALSE,     "degrees",           "Slope, COP30-derived",
  "northness",               NA_character_,                                                       NA_character_,                    NA_character_, "continuous",  "bilinear",  "terrestrial",  TRUE,             FALSE,     "unitless",          "cos(aspect), derived from COP30 DEM",
  "eastness",                NA_character_,                                                       NA_character_,                    NA_character_, "continuous",  "bilinear",  "terrestrial",  TRUE,             FALSE,     "unitless",          "sin(aspect), derived from COP30 DEM",
  "twi",                     "outputs/catchments/twi_cop30.tif",                                  NA_character_,                    NA_character_, "continuous",  "bilinear",  "terrestrial",  TRUE,             FALSE,     "index",             "Topographic wetness index, COP30-derived (raw, not pre-normalised)",
  "dist_coast",              "data/processed/boundaries/scamp_coastline.gpkg",                    "scamp_coastline",               NA_character_, "continuous",  "bilinear",  "terrestrial",  TRUE,             FALSE,     "m",                 "Distance to coastline",
  "diffuse_pollution_risk",  "outputs/catchments/diffuse_pollution_risk.tif",                     NA_character_,                    NA_character_, "continuous",  "bilinear",  "terrestrial",  TRUE,             FALSE,     "index",             "SCIMAP diffuse pollution risk, 20 m",
  "runoff_risk",             "outputs/catchments/hlcm_source_runoff_risk.tif",                    NA_character_,                    NA_character_, "continuous",  "bilinear",  "terrestrial",  TRUE,             FALSE,     "index",             "HLCM source runoff risk, 20 m",
  "land_cover",              "data/processed/habitats/hlcm_habmos_lca_reclass.tif",               NA_character_,                    NA_character_, "categorical", "mode",      "terrestrial",  TRUE,             FALSE,     "class",             "HLCM EUNIS L2 land cover, reclassified (grassland-intensity split, G4.F); rare levels collapsed",
  "bio1",                    "data/raw/climate",                                                  NA_character_,                    NA_character_, "continuous",  "bilinear",  "terrestrial",  TRUE,             FALSE,     "degC",              "WorldClim v2.1 30 arc-sec annual mean temperature",
  "bio12",                   "data/raw/climate",                                                  NA_character_,                    NA_character_, "continuous",  "bilinear",  "terrestrial",  TRUE,             FALSE,     "mm",                "WorldClim v2.1 30 arc-sec annual precipitation",

  "erosion_accretion_class", "data/raw/habitats/Sand dunes/DYNAMICCOAST2_SCOTLAND_SHP_27700/shp/DC2_RCP8_Transects.shp", NA_character_, "HIST_RATE", "categorical", "mode", "terrestrial", FALSE,            FALSE,     "class",             "Dynamic Coast 2 historical shoreline rate of change (DC2_RCP8_Transects, HIST_RATE), binned at +/-DC_STABLE_THRESHOLD_MYR m/yr and rasterised onto the cells each transect crosses. Cells not crossed by a transect have no value, meaning not assessed. Cube-only; not used for banding."
)

cat(sprintf("Covariate table: %d rows (%d marine, %d terrestrial)\n",
            nrow(covariate_table),
            sum(covariate_table$cube == "marine"),
            sum(covariate_table$cube == "terrestrial")))

##=====================================================================
## 2. AUDIT — confirm each covariate's source exists; skip, don't invent
##=====================================================================

audit_covariate <- function(row) {
  if (is.na(row$path)) return(list(ok = TRUE, note = "derived in-script"))
  if (!file.exists(row$path) && !dir.exists(row$path)) {
    return(list(ok = FALSE, note = sprintf("MISSING: %s", row$path)))
  }
  if (grepl("\\.gpkg$", row$path)) {
    layers <- tryCatch(st_layers(row$path)$name, error = function(e) character(0))
    if (!row$layer %in% layers) {
      return(list(ok = FALSE, note = sprintf("layer '%s' not found in %s", row$layer, row$path)))
    }
    fields <- tryCatch(names(st_read(row$path, layer = row$layer, query =
                         sprintf("SELECT * FROM \"%s\" LIMIT 0", row$layer), quiet = TRUE)),
                        error = function(e) character(0))
    if (!is.na(row$field) && !row$field %in% fields) {
      return(list(ok = FALSE, note = sprintf("field '%s' not found in layer '%s'", row$field, row$layer)))
    }
  }
  if (grepl("\\.shp$", row$path)) {
    shp_layer <- tools::file_path_sans_ext(basename(row$path))
    fields <- tryCatch(
      names(st_read(row$path, query = sprintf("SELECT * FROM \"%s\" LIMIT 0", shp_layer), quiet = TRUE)),
      error = function(e) character(0)
    )
    if (!is.na(row$field) && !row$field %in% fields) {
      return(list(ok = FALSE, note = sprintf("field '%s' not found in %s", row$field, row$path)))
    }
  }
  list(ok = TRUE, note = "found")
}

cat("\n-- Audit: checking covariate sources on disk --\n")
audit_results <- lapply(seq_len(nrow(covariate_table)), function(i) audit_covariate(covariate_table[i, ]))
covariate_table$audit_ok   <- vapply(audit_results, function(x) x$ok, logical(1))
covariate_table$audit_note <- vapply(audit_results, function(x) x$note, character(1))

for (i in seq_len(nrow(covariate_table))) {
  flag <- if (covariate_table$audit_ok[i]) "OK  " else "SKIP"
  cat(sprintf("  [%s] %-24s %s\n", flag, covariate_table$name[i], covariate_table$audit_note[i]))
}

if (any(!covariate_table$audit_ok)) {
  cat("\nWARNING: the above covariates are missing and will be skipped:\n")
  cat(paste0("  - ", covariate_table$name[!covariate_table$audit_ok], "\n"))
}

covariate_table_active <- covariate_table %>% filter(audit_ok)

##=====================================================================
## 3. PROGRAMME AREA — rebuild marine/land masks (scamp_aoi has no split field)
##=====================================================================

cat("\n-- Programme area: rebuilding marine/land masks --\n")
cat("  (scamp_aoi.gpkg carries fid/geom only; replicating 02_themeD_species.R:47-75)\n")

smr <- st_read(
  "data/raw/boundaries/administrative_units_scottish_marine_regions/administrative_units_scottish_marine_regionsPolygon.shp",
  quiet = TRUE
)
las <- st_read(
  "data/raw/boundaries/Local_Authority_Boundaries_-_Scotland/pub_las.shp",
  quiet = TRUE
)

smr_solway_bng <- smr %>% filter(grepl("Solway", objnam, ignore.case = TRUE)) %>% st_transform(27700)
las_dg_bng     <- las %>% filter(grepl("dumfries", local_auth, ignore.case = TRUE)) %>% st_transform(27700)

marine_aoi_sf <- st_sf(geometry = st_union(st_geometry(smr_solway_bng)))
land_aoi_sf   <- st_sf(geometry = st_union(st_intersection(
                   st_geometry(st_buffer(smr_solway_bng, dist = 8047)),
                   st_geometry(las_dg_bng)
                 )))

marine_aoi_v <- vect(marine_aoi_sf)
land_aoi_v   <- vect(land_aoi_sf)

##=====================================================================
## 4. GRIDS
##=====================================================================

cat("\n-- Building target grids --\n")

marine_template <- rast("outputs/physical/gebco_bathy_solway.tif")
if (is.na(crs(marine_template)) || crs(marine_template) == "") crs(marine_template) <- "EPSG:27700"

cop30_dem <- rast("outputs/physical/cop30_dem_solway.tif")
if (is.na(crs(cop30_dem)) || crs(cop30_dem) == "") crs(cop30_dem) <- "EPSG:27700"
terr_template <- rast(ext(cop30_dem), resolution = 100, crs = "EPSG:27700")

cat(sprintf("  Marine template : %d x %d cells, %.0f m res\n",
            ncol(marine_template), nrow(marine_template), res(marine_template)[1]))
cat(sprintf("  Terrestrial grid: %d x %d cells, %.0f m res\n",
            ncol(terr_template), nrow(terr_template), res(terr_template)[1]))

marine_mask_r <- rasterize(marine_aoi_v, marine_template, field = 1, background = NA)
land_mask_r   <- rasterize(land_aoi_v, terr_template, field = 1, background = NA)

##=====================================================================
## 5. HELPER FUNCTIONS
##=====================================================================

## dplyr/rlang do not export `%||%` and base R only gained it in 4.4 --
## define it locally so the script runs on older R.
if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

pct_coverage <- function(r, mask_r) {
  masked <- mask(r, mask_r)
  n_total <- global(!is.na(mask_r), "sum", na.rm = TRUE)[[1]]
  n_data  <- global(!is.na(masked), "sum", na.rm = TRUE)[[1]]
  if (n_total == 0) return(NA_real_)
  100 * n_data / n_total
}

gapfill_focal <- function(r, mask_r, max_passes = MARINE_GAPFILL_MAX_PASSES) {
  r <- mask(r, mask_r)
  n_na_before <- global(is.na(r), "sum", na.rm = TRUE)[[1]]
  if (n_na_before == 0) {
    return(list(filled = r, was_filled = r * 0, passes = 0, n_filled = 0))
  }
  was_na <- is.na(r)
  filled <- focal(r, w = 3, fun = "mean", na.policy = "only", na.rm = TRUE)
  n_remaining <- global(is.na(filled) & !is.na(mask_r), "sum", na.rm = TRUE)[[1]]
  passes <- 1
  while (n_remaining > 0 && passes < max_passes) {
    filled <- focal(filled, w = 3, fun = "mean", na.policy = "only", na.rm = TRUE)
    n_remaining <- global(is.na(filled) & !is.na(mask_r), "sum", na.rm = TRUE)[[1]]
    passes <- passes + 1
  }
  filled <- mask(filled, mask_r)
  was_filled <- was_na & !is.na(filled)
  n_filled <- global(was_filled, "sum", na.rm = TRUE)[[1]]
  list(filled = filled, was_filled = was_filled, passes = passes, n_filled = n_filled)
}

load_continuous <- function(path, template) {
  r <- rast(path)
  if (is.na(crs(r)) || crs(r) == "") crs(r) <- "EPSG:27700"
  project(r, template, method = "bilinear")
}

load_categorical_vector <- function(path, layer, field, template, na_values = character(0)) {
  v <- st_read(path, layer = layer, quiet = TRUE) %>% st_transform(27700)
  vals <- v[[field]]
  vals[vals %in% na_values] <- NA
  v[[field]] <- vals
  vv <- vect(v)
  r <- rasterize(vv, template, field = field)
  r
}

## Category table of a categorical SpatRaster as (code, label), using the
## active category column for the label.
cat_codes <- function(r) {
  ct <- cats(r)[[1]]
  if (is.null(ct) || nrow(ct) == 0) stop(sprintf("'%s' has no category table", names(r)))
  out <- data.frame(code = ct[[1]], label = as.character(ct[[activeCat(r) + 1]]),
                    stringsAsFactors = FALSE)
  out[!is.na(out$label), , drop = FALSE]
}
## Rare levels are collapsed by rebuilding the category table, not by label
## substitution: each rare level's integer code is classify()'d to a single new
## "other" code, then levels() is set explicitly to the kept levels + "other".
## Levels absent from the sample are treated as rare (sample share 0).
collapse_rare_levels <- function(r, template_mask, min_prop = RARE_LEVEL_MIN_PROP) {
  ct <- cats(r)[[1]]
  lab_col <- names(ct)[activeCat(r) + 1]
  cc <- cat_codes(r)
  cat("    cats() before collapse:\n"); print(ct)

  r_num <- r
  levels(r_num) <- NULL
  samp <- spatSample(mask(r_num, template_mask), size = min(SAMPLE_N, ncell(r)),
                      method = "regular", na.rm = TRUE, as.df = TRUE)
  if (nrow(samp) == 0) return(r)
  tab <- prop.table(table(samp[[1]]))
  keep_codes <- as.numeric(names(tab)[tab >= min_prop])
  rare_codes <- setdiff(cc$code, keep_codes)
  if (length(rare_codes) == 0) {
    cat("    No rare levels to collapse\n")
    return(r)
  }
  cat(sprintf("    Collapsing rare levels into 'other': %s\n",
              paste(cc$label[cc$code %in% rare_codes], collapse = ", ")))

  other_code <- max(cc$code) + 1
  out <- classify(r_num, rcl = cbind(rare_codes, other_code))
  kept <- cc[cc$code %in% keep_codes, ]
  new_levels <- data.frame(id = c(kept$code, other_code))
  new_levels[[lab_col]] <- c(kept$label, "other")
  levels(out) <- new_levels
  names(out) <- names(r)
  cat("    cats() after collapse:\n"); print(cats(out)[[1]])
  out
}

## All rasters written to outputs/design/ go through here: COG, DEFLATE, band
## names carried from names(r), and the DEFLATE predictor matched to the data
## type -- PREDICTOR=3 (floating-point) for float layers, PREDICTOR=2
## (horizontal differencing) for integer layers. The type is taken from the
## `datatype` argument when one is given, otherwise from the raster itself.
## A single COG holds one data type, so a mixed stack (the cubes) is written as
## float and takes PREDICTOR=3: hence "any layer float -> 3".
write_cog <- function(r, path, datatype = NULL) {
  dt <- if (!is.null(datatype)) datatype else
    tryCatch(terra::datatype(r), error = function(e) character(0))
  dt <- dt[!is.na(dt) & nzchar(dt)]
  is_float <- if (length(dt) > 0) any(grepl("^FLT", dt)) else
    !all(tryCatch(terra::is.int(r), error = function(e) FALSE))
  predictor <- if (is_float) "PREDICTOR=3" else "PREDICTOR=2"
  cat(sprintf("    COG: %s (DEFLATE, %s)\n", path, predictor))
  args <- list(x = r, filename = path, filetype = "COG", overwrite = TRUE,
               gdal = c("COMPRESS=DEFLATE", predictor), names = names(r))
  if (!is.null(datatype)) args$datatype <- datatype
  do.call(writeRaster, args)
}

## The six SCAMP AOI bays, ordered west to east. This vector fixes both the
## required contents of aois_named.gpkg and the output order of every
## representativeness CSV row and plot legend entry.
AOI_BAYS_WEST_TO_EAST <- c("Loch Ryan", "Luce Bay", "Wigtown Bay",
                           "Fleet Bay", "Kirkcudbright Bay", "Auchencairn Bay")

## AOI bay polygons, EPSG:27700, with an aoi_name column taken from the
## gpkg's own bay_name field. aois_named.gpkg is the ONLY source: there is no
## shapefile fallback and no positional naming, because a bay name assigned by
## row position would silently mislabel every downstream representativeness
## row and station plan without erroring. Anything wrong with the file is a
## hard stop telling the user to rebuild it.
load_aoi_bays <- function(path = "data/processed/boundaries/aois_named.gpkg") {
  rebuild <- function(why) stop(sprintf(paste0(
    "load_aoi_bays: %s\n",
    "  Rebuild %s: one polygon per bay, with a 'bay_name' column naming exactly\n",
    "  these six bays -- %s"),
    why, path, paste(AOI_BAYS_WEST_TO_EAST, collapse = ", ")), call. = FALSE)

  if (!file.exists(path)) rebuild("file not found.")
  aois <- tryCatch(st_read(path, quiet = TRUE), error = function(e)
    rebuild(sprintf("file could not be read (%s).", conditionMessage(e))))

  if (!any(st_geometry_type(aois) %in% c("POLYGON", "MULTIPOLYGON"))) {
    rebuild("layer holds no POLYGON/MULTIPOLYGON features.")
  }
  if (!"bay_name" %in% names(aois)) {
    rebuild(sprintf("no 'bay_name' column (columns present: %s).",
                    paste(names(aois), collapse = ", ")))
  }

  aois$aoi_name <- trimws(as.character(aois$bay_name))
  if (nrow(aois) != length(AOI_BAYS_WEST_TO_EAST) ||
      !setequal(aois$aoi_name, AOI_BAYS_WEST_TO_EAST)) {
    rebuild(sprintf("bay_name does not hold exactly the six expected bays -- found %d row(s): %s.",
                    nrow(aois), paste(sort(aois$aoi_name), collapse = ", ")))
  }

  aois <- st_make_valid(st_transform(aois, 27700))[, "aoi_name"]
  aois <- aois[match(AOI_BAYS_WEST_TO_EAST, aois$aoi_name), ]  # west -> east
  cat(sprintf("  AOI source: %s\n  Bays (west to east): %s\n",
              path, paste(aois$aoi_name, collapse = ", ")))
  aois
}

derive_cmems_metrics <- function(nc_dir, template, mask_r) {
  files <- list.files(nc_dir, pattern = "\\.nc$", full.names = TRUE)
  files <- files[!grepl("\\(1\\)", files)]
  stopifnot("CMEMS NetCDF not found" = length(files) >= 1)
  cat(sprintf("    Reading %d file(s) from %s:\n", length(files), nc_dir))
  cat(paste0("      ", basename(files), "\n"), sep = "")
  stack <- tryCatch(rast(files), error = function(e) {
    stop(sprintf("CMEMS files in %s could not be combined (different grids?): %s",
                 nc_dir, conditionMessage(e)))
  })
  ## Depth dimension: keep ONLY the shallowest level. thetao is 4D
  ## (time x depth); bottomT is 2D and passes through untouched. The level is
  ## identified from terra::depth() where the NetCDF carries it, else from a
  ## "depth=" token in the layer names -- never by layer position. This must
  ## run BEFORE the time de-duplication below: with several depth levels the
  ## same timestamps repeat per level, and !duplicated(tm) would otherwise
  ## keep an arbitrary mixture of depths.
  depth_vals <- tryCatch(terra::depth(stack), error = function(e) numeric(0))
  if (length(depth_vals) != nlyr(stack) || all(is.na(depth_vals))) {
    nm <- names(stack)
    depth_vals <- ifelse(grepl("depth=", nm),
                         suppressWarnings(as.numeric(sub(".*depth=([0-9eE.+-]+).*", "\\1", nm))),
                         NA_real_)
  }
  u_depth <- sort(unique(depth_vals[!is.na(depth_vals)]))
  if (length(u_depth) > 1) {
    shallowest <- u_depth[1]
    cat(sprintf("    Depth dimension: %d levels (%g to %g m); keeping shallowest = %g m\n",
                length(u_depth), u_depth[1], u_depth[length(u_depth)], shallowest))
    stack <- stack[[which(!is.na(depth_vals) & depth_vals == shallowest)]]
  } else if (length(u_depth) == 1) {
    cat(sprintf("    Single depth level (%g m) -- all layers kept\n", u_depth[1]))
  } else {
    cat("    No depth dimension found -- 2D variable, all layers kept\n")
  }

  tm <- time(stack)
  if (all(is.na(tm))) stop(sprintf("CMEMS files in %s carry no time dimension", nc_dir))

  ## Combine along time: drop duplicated / missing time steps, then sort
  keep <- which(!duplicated(tm) & !is.na(tm))
  n_dup <- length(tm) - length(keep)
  stack <- stack[[keep]]; tm <- tm[keep]
  o <- order(tm)
  stack <- stack[[o]]; tm <- tm[o]

  years <- as.numeric(format(tm, "%Y"))
  uyears <- sort(unique(years))
  cat(sprintf("    %d time steps (%d duplicate/undated dropped), years used: %d-%d\n",
              length(tm), n_dup, min(uyears), max(uyears)))

  ## Crop (NO mask) to the marine AOI buffered by CMEMS_CROP_BUFFER_M. The
  ## CMEMS grid is ~9 km, so mask = TRUE at native resolution would discard
  ## most cells covering the Solway's inlets -- their centres fall outside the
  ## AOI even though the cell overlaps it -- leaving holes that gap-filling
  ## then has to invent. Keeping the surrounding cells lets the bilinear
  ## projection interpolate into the inlets properly; the result is masked to
  ## marine_mask_r only after it is on the 286 m template.
  aoi_buffered <- st_buffer(marine_aoi_sf, dist = CMEMS_CROP_BUFFER_M)
  aoi_native <- project(vect(aoi_buffered), crs(stack))
  cropped <- crop(stack, aoi_native)

  overall_mean <- app(cropped, fun = "mean", na.rm = TRUE)

  annual_max <- rast(); annual_min <- rast()
  for (yr in uyears) {
    yr_layers <- cropped[[which(years == yr)]]
    annual_max <- c(annual_max, app(yr_layers, fun = "max", na.rm = TRUE))
    annual_min <- c(annual_min, app(yr_layers, fun = "min", na.rm = TRUE))
  }
  clim_max <- app(annual_max, fun = "mean", na.rm = TRUE)
  clim_min <- app(annual_min, fun = "mean", na.rm = TRUE)
  clim_range <- clim_max - clim_min

  ## Project to the 286 m template first, mask to the programme area second.
  to_template <- function(x) mask(project(x, template, method = "bilinear"), mask_r)
  list(
    mean  = to_template(overall_mean),
    max   = to_template(clim_max),
    range = to_template(clim_range)
  )
}

##=====================================================================
## 6. LOAD + ALIGN MARINE COVARIATES
##=====================================================================

cat("\n-- Loading marine covariates --\n")

marine_rasters <- list()
marine_fill_log <- list()

m_rows <- covariate_table_active %>% filter(cube == "marine", type == "continuous", is.na(layer))
cmems_done <- character(0)

for (i in seq_len(nrow(m_rows))) {
  row <- m_rows[i, ]
  cat(sprintf("  %s\n", row$name))

  if (row$name %in% c("sst_mean", "sst_max", "sst_range")) {
    if (!"sst" %in% cmems_done) {
      sst <- derive_cmems_metrics("data/processed/physical/thetao", marine_template, marine_mask_r)
      marine_rasters[["sst_mean"]]  <- sst$mean
      marine_rasters[["sst_max"]]   <- sst$max
      marine_rasters[["sst_range"]] <- sst$range
      cmems_done <- c(cmems_done, "sst")
    }
    next
  }
  if (row$name %in% c("bottomt_mean", "bottomt_max", "bottomt_range")) {
    if (!"bottomt" %in% cmems_done) {
      bt <- derive_cmems_metrics("data/processed/physical/bottomT", marine_template, marine_mask_r)
      marine_rasters[["bottomt_mean"]]  <- bt$mean
      marine_rasters[["bottomt_max"]]   <- bt$max
      marine_rasters[["bottomt_range"]] <- bt$range
      cmems_done <- c(cmems_done, "bottomt")
    }
    next
  }

  marine_rasters[[row$name]] <- load_continuous(row$path, marine_template)
}

## Categorical marine covariates
cat_rows_m <- covariate_table_active %>% filter(cube == "marine", type == "categorical")
na_lookup <- list(energy_class = "No energy information", seabed_substrate = "unknown substrate")
for (i in seq_len(nrow(cat_rows_m))) {
  row <- cat_rows_m[i, ]
  cat(sprintf("  %s (categorical)\n", row$name))
  r <- load_categorical_vector(row$path, row$layer, row$field, marine_template,
                                na_values = na_lookup[[row$name]] %||% character(0))
  r <- collapse_rare_levels(r, marine_mask_r)
  marine_rasters[[row$name]] <- r
}

## Coverage BEFORE gap-fill
coverage_rows <- list()
for (nm in names(marine_rasters)) {
  coverage_rows[[length(coverage_rows) + 1]] <- data.frame(
    cube = "marine", covariate = nm,
    coverage_before_pct = pct_coverage(marine_rasters[[nm]], marine_mask_r),
    gap_fill_applied = FALSE, passes = NA_integer_, coverage_after_pct = NA_real_
  )
}

## Gap-fill (continuous CMEMS-derived only, per covariate_table$gap_fill)
cat("\n-- Gap-filling marine covariates (CMEMS-derived only) --\n")
fill_masks <- list()
gap_fill_names <- covariate_table_active %>% filter(cube == "marine", gap_fill) %>% pull(name)
for (nm in gap_fill_names) {
  if (!nm %in% names(marine_rasters)) next
  res_fill <- gapfill_focal(marine_rasters[[nm]], marine_mask_r)
  marine_rasters[[nm]] <- res_fill$filled
  fill_masks[[nm]] <- res_fill$was_filled
  n_total <- global(!is.na(marine_mask_r), "sum", na.rm = TRUE)[[1]]
  cat(sprintf("  %-20s filled %d/%0.f cells in %d passes\n", nm, res_fill$n_filled, n_total, res_fill$passes))
  idx <- which(sapply(coverage_rows, function(x) x$covariate == nm))
  coverage_rows[[idx]]$gap_fill_applied <- TRUE
  coverage_rows[[idx]]$passes <- res_fill$passes
  coverage_rows[[idx]]$coverage_after_pct <- pct_coverage(marine_rasters[[nm]], marine_mask_r)
}

## Force each raster's own internal layer name to match its covariate name
## before stacking -- rast(path) keeps the source file's band name otherwise,
## which would break every [["covariate_name"]] lookup downstream.
for (nm in names(marine_rasters)) names(marine_rasters[[nm]]) <- nm
marine_stack <- rast(marine_rasters)
marine_stack <- mask(marine_stack, marine_mask_r)

## n_gapfilled band: per-cell count of covariates that were gap-filled
if (length(fill_masks) > 0) {
  n_gapfilled <- app(rast(fill_masks), fun = "sum", na.rm = TRUE)
} else {
  n_gapfilled <- marine_mask_r * 0
}
names(n_gapfilled) <- "n_gapfilled"

## Expected banding coverage: cells with ALL use_for_banding marine covariates present
banding_names_m <- covariate_table_active %>% filter(cube == "marine", use_for_banding) %>% pull(name)
banding_names_m <- intersect(banding_names_m, names(marine_rasters))
complete_m <- app(marine_stack[[banding_names_m]], fun = function(x) all(!is.na(x)))
pct_banded_marine <- 100 * global(complete_m, "sum", na.rm = TRUE)[[1]] /
                      global(!is.na(marine_mask_r), "sum", na.rm = TRUE)[[1]]
cat(sprintf("\n  Expected marine cube coverage receiving a band: %.1f%%\n", pct_banded_marine))

##=====================================================================
## 7. LOAD + ALIGN TERRESTRIAL COVARIATES
##=====================================================================

cat("\n-- Loading terrestrial covariates --\n")

terr_rasters <- list()
t_rows <- covariate_table_active %>% filter(cube == "terrestrial", type == "continuous",
                                             !name %in% c("northness", "eastness", "dist_coast", "bio1", "bio12"))
for (i in seq_len(nrow(t_rows))) {
  row <- t_rows[i, ]
  cat(sprintf("  %s\n", row$name))
  terr_rasters[[row$name]] <- load_continuous(row$path, terr_template)
}

## Northness / eastness from COP30 DEM aspect
cat("  northness / eastness (derived from COP30 aspect)\n")
aspect <- terrain(cop30_dem, v = "aspect", unit = "radians")
aspect_bng <- project(aspect, terr_template, method = "bilinear")
terr_rasters[["northness"]] <- cos(aspect_bng)
terr_rasters[["eastness"]]  <- sin(aspect_bng)
names(terr_rasters[["northness"]]) <- "northness"
names(terr_rasters[["eastness"]])  <- "eastness"

## Distance to coastline
cat("  dist_coast (derived from scamp_coastline.gpkg)\n")
coastline <- st_read("data/processed/boundaries/scamp_coastline.gpkg", quiet = TRUE) %>% st_transform(27700)
terr_rasters[["dist_coast"]] <- terra::distance(terr_template, vect(coastline))
names(terr_rasters[["dist_coast"]]) <- "dist_coast"

## WorldClim bio1 / bio12
cat("  bio1 / bio12 (WorldClim v2.1, cached under data/raw/climate/)\n")
if (requireNamespace("geodata", quietly = TRUE)) {
  wc <- geodata::worldclim_country("GBR", var = "bio", res = 0.5, path = "data/raw/climate")
  wc_bio1  <- wc[["wc2.1_30s_bio_1"]]
  wc_bio12 <- wc[["wc2.1_30s_bio_12"]]
  terr_rasters[["bio1"]]  <- project(wc_bio1,  terr_template, method = "bilinear")
  terr_rasters[["bio12"]] <- project(wc_bio12, terr_template, method = "bilinear")
  names(terr_rasters[["bio1"]])  <- "bio1"
  names(terr_rasters[["bio12"]]) <- "bio12"
} else {
  warning("geodata package not installed; skipping bio1/bio12")
}

## Categorical: land_cover
cat("  land_cover (categorical)\n")
land_cover_lookup <- read.csv("data/processed/habitats/hlcm_lc_lookup_reclass.csv")
lc_r <- rast("data/processed/habitats/hlcm_habmos_lca_reclass.tif")
if (is.na(crs(lc_r)) || crs(lc_r) == "") crs(lc_r) <- "EPSG:27700"
levels(lc_r) <- data.frame(id = land_cover_lookup$lc_value, class = land_cover_lookup$eunis_code)
lc_r_terr <- project(lc_r, terr_template, method = "mode")
lc_r_terr <- collapse_rare_levels(lc_r_terr, land_mask_r)
terr_rasters[["land_cover"]] <- lc_r_terr

for (nm in names(terr_rasters)) names(terr_rasters[[nm]]) <- nm
terr_stack <- rast(terr_rasters)
terr_stack <- mask(terr_stack, land_mask_r)

coverage_rows_t <- list()
for (nm in names(terr_rasters)) {
  coverage_rows_t[[length(coverage_rows_t) + 1]] <- data.frame(
    cube = "terrestrial", covariate = nm,
    coverage_before_pct = pct_coverage(terr_rasters[[nm]], land_mask_r),
    gap_fill_applied = FALSE, passes = NA_integer_, coverage_after_pct = NA_real_
  )
}

banding_names_t <- covariate_table_active %>% filter(cube == "terrestrial", use_for_banding) %>% pull(name)
banding_names_t <- intersect(banding_names_t, names(terr_rasters))
complete_t <- app(terr_stack[[banding_names_t]], fun = function(x) all(!is.na(x)))
pct_banded_terr <- 100 * global(complete_t, "sum", na.rm = TRUE)[[1]] /
                    global(!is.na(land_mask_r), "sum", na.rm = TRUE)[[1]]
cat(sprintf("\n  Expected terrestrial cube coverage receiving a band: %.1f%%\n", pct_banded_terr))

## Cube-only terrestrial layer: Dynamic Coast erosion/accretion class
dc_row <- covariate_table_active %>% filter(name == "erosion_accretion_class")
dc_class_r <- NULL
if (nrow(dc_row) == 1) {
  cat("\n  erosion_accretion_class (cube-only, not used for banding)\n")
  dc <- st_read(dc_row$path, quiet = TRUE) %>% st_transform(27700)
  dc$band_class <- cut(dc[[dc_row$field]],
                        breaks = c(-Inf, -DC_STABLE_THRESHOLD_MYR, DC_STABLE_THRESHOLD_MYR, Inf),
                        labels = c("eroding", "stable", "accreting"))
  dc_class_r <- rasterize(vect(dc), terr_template, field = "band_class", touches = TRUE)
  dc_class_r <- mask(dc_class_r, land_mask_r)
  names(dc_class_r) <- "erosion_accretion_class"
}

##=====================================================================
## 8. COVERAGE REPORT
##=====================================================================

coverage_report <- bind_rows(coverage_rows, coverage_rows_t) %>%
  left_join(covariate_table %>% select(name, use_for_banding), by = c("covariate" = "name"))

readr_write <- function(df, path) write.csv(df, path, row.names = FALSE)
readr_write(coverage_report, "outputs/design/coverage_report.csv")

low_cov <- coverage_report %>%
  mutate(final_pct = ifelse(gap_fill_applied, coverage_after_pct, coverage_before_pct)) %>%
  filter(use_for_banding, final_pct < 95)
if (nrow(low_cov) > 0) {
  cat("\nWARNING: banding covariates below 95% coverage (will shrink the banded area):\n")
  print(low_cov %>% select(cube, covariate, final_pct))
}

cat(sprintf("\nCoverage report written: outputs/design/coverage_report.csv\n"))
cat(sprintf("Expected banded area -- marine: %.1f%%, terrestrial: %.1f%%\n",
            pct_banded_marine, pct_banded_terr))

##=====================================================================
## 9. SAMPLING, SCALING, ONE-HOT ENCODING
##=====================================================================

build_feature_frame <- function(stack, banding_names, mask_r, sample_n = SAMPLE_N) {
  s <- mask(stack[[banding_names]], mask_r)
  df <- spatSample(s, size = min(sample_n, ncell(s)), method = "regular",
                    na.rm = TRUE, as.df = TRUE)
  df
}

## ONE TRANSFORM RULE, used identically for the training sample (here) and the
## prediction raster (make_z_stack()):
##   - continuous banding columns: z-score with the sample mean and sd
##   - one-hot columns: NOT z-scored; value = CAT_WEIGHT * indicator, i.e.
##     0 or CAT_WEIGHT
## One-hot column names are cov__level, with the level label sanitised to
## [A-Za-z0-9_] so names survive terra::extract()/data.frame round trips.
## onehot_map keeps the raw label for each column so prediction matches by
## label.
onehot_colname <- function(cov, label) paste0(cov, "__", gsub("[^A-Za-z0-9]+", "_", label))

one_hot <- function(df, cat_cols, weight = CAT_WEIGHT) {
  map <- list()
  for (cl in cat_cols) {
    v <- as.character(df[[cl]])
    lv <- sort(unique(na.omit(v)))
    cn <- onehot_colname(cl, lv)
    if (anyDuplicated(cn)) stop(sprintf("one_hot: level labels of '%s' collide after sanitising: %s",
                                        cl, paste(lv[duplicated(cn) | duplicated(cn, fromLast = TRUE)], collapse = ", ")))
    for (j in seq_along(lv)) df[[cn[j]]] <- weight * as.numeric(v == lv[j])
    map[[cl]] <- data.frame(column = cn, covariate = cl, level = lv, stringsAsFactors = FALSE)
    df[[cl]] <- NULL
  }
  list(df = df, onehot_map = bind_rows(map))
}

fit_scaling <- function(df, cont_cols) {
  means <- sapply(df[cont_cols], mean, na.rm = TRUE)
  sds   <- sapply(df[cont_cols], sd, na.rm = TRUE)
  sds[sds == 0 | is.na(sds)] <- 1
  list(means = means, sds = sds)
}

## z-scores the continuous columns only; one-hot columns pass through as-is.
apply_scaling_df <- function(df, scaling) {
  for (cl in names(scaling$means)) {
    df[[cl]] <- (df[[cl]] - scaling$means[[cl]]) / scaling$sds[[cl]]
  }
  df
}

cat_cols_m <- intersect(cat_rows_m$name, banding_names_m)
cat_cols_t <- intersect("land_cover", banding_names_t)

cat("\n-- Fitting scaling + sampling for k-selection --\n")
samp_m <- build_feature_frame(marine_stack, banding_names_m, marine_mask_r)
oh_m <- one_hot(samp_m, cat_cols_m)
onehot_map_m <- oh_m$onehot_map
scaling_m <- fit_scaling(oh_m$df, setdiff(banding_names_m, cat_cols_m))
z_samp_m <- apply_scaling_df(oh_m$df, scaling_m)

samp_t <- build_feature_frame(terr_stack, banding_names_t, land_mask_r)
oh_t <- one_hot(samp_t, cat_cols_t)
onehot_map_t <- oh_t$onehot_map
scaling_t <- fit_scaling(oh_t$df, setdiff(banding_names_t, cat_cols_t))
z_samp_t <- apply_scaling_df(oh_t$df, scaling_t)

##=====================================================================
## 10. K SELECTION (elbow)
##=====================================================================

run_k_screen <- function(z_samp, k_range = 2:12, seed = SEED) {
  set.seed(seed)
  res <- lapply(k_range, function(k) {
    km <- kmeans(z_samp, centers = k, nstart = 10, iter.max = 50)
    data.frame(k = k, tot_withinss = km$tot.withinss)
  })
  bind_rows(res)
}

cat("\n-- K selection: marine --\n")
k_screen_m <- run_k_screen(z_samp_m)
p_m <- ggplot(k_screen_m, aes(k, tot_withinss)) + geom_line() + geom_point() +
  scale_x_continuous(breaks = 2:12) +
  labs(title = "Marine band k-selection", x = "k", y = "Total within-cluster SS") +
  theme_minimal()
ggsave("outputs/design/k_selection_marine.png", p_m, width = 8, height = 5, dpi = 300)

cat("-- K selection: terrestrial --\n")
k_screen_t <- run_k_screen(z_samp_t)
p_t <- ggplot(k_screen_t, aes(k, tot_withinss)) + geom_line() + geom_point() +
  scale_x_continuous(breaks = 2:12) +
  labs(title = "Terrestrial band k-selection", x = "k", y = "Total within-cluster SS") +
  theme_minimal()
ggsave("outputs/design/k_selection_terrestrial.png", p_t, width = 8, height = 5, dpi = 300)

cat("\nElbow plots written to outputs/design/. Set K_MARINE and K_TERR at the\n")
cat("top of this script and re-run for band assignment.\n")

##=====================================================================
## 10b. CONTROL MATCHING + DESIGN FUNCTIONS
##=====================================================================
## Defined here, before the K gate, so pass 2 can call compare_aoi_to_region()
## once the cubes are written. Definitions only -- nothing runs in pass 1.

#' Find candidate control cells for a set of intervention sites
#'
#' @param sites sf POINT, one row per intervention site. Must have site_id,
#'   or one is added (seq_len).
#' @param band_rast SpatRaster, single band ("band"), INT1U.
#' @param z_stack SpatRaster, z-scored covariate stack used to fit the bands
#'   (same layers as make_z_stack() output) -- used for in-band ranking.
#' @param n integer, candidates returned per site.
#' @param min_dist numeric (m), candidates within this distance of ANY site
#'   in `sites` are excluded (buffer against spatial autocorrelation/spillover).
#' @param mask optional sf/SpatVector restricting candidates (e.g. a habitat
#'   extent polygon). NULL = no restriction beyond band_rast's own extent.
#' @param exclude optional sf/SpatVector of areas to exclude (e.g. other
#'   planned intervention footprints, designated-site buffers).
#' @return sf POINT with site_id, band, rank, cov_distance. Sites with zero
#'   eligible candidates in their band produce a warning and no rows.
match_controls <- function(sites, band_rast, z_stack, n = 5, min_dist = 0,
                            mask = NULL, exclude = NULL) {
  stopifnot(inherits(sites, "sf"))
  if (!"site_id" %in% names(sites)) sites$site_id <- seq_len(nrow(sites))

  cand_stack <- c(band_rast, z_stack)
  if (!is.null(mask)) cand_stack <- mask(cand_stack, vect(mask))
  if (!is.null(exclude)) cand_stack <- mask(cand_stack, vect(exclude), inverse = TRUE)

  cand_df <- as.data.frame(cand_stack, xy = TRUE, na.rm = TRUE)
  if (nrow(cand_df) == 0) {
    warning("match_controls: no candidate cells after masking/exclusion")
    return(sites[0, c("site_id"), drop = FALSE])
  }
  cand_pts <- st_as_sf(cand_df, coords = c("x", "y"), crs = crs(band_rast), remove = FALSE)

  if (min_dist > 0) {
    d_to_sites <- st_distance(cand_pts, sites)
    min_d <- apply(d_to_sites, 1, min)
    cand_pts <- cand_pts[as.numeric(min_d) >= min_dist, ]
    cand_df  <- cand_df[as.numeric(min_d) >= min_dist, ]
  }

  cov_names <- setdiff(names(z_stack), "band")
  sites_vals <- terra::extract(z_stack, vect(sites))[, cov_names, drop = FALSE]
  sites_band <- terra::extract(band_rast, vect(sites))[, 2]

  out <- list()
  for (i in seq_len(nrow(sites))) {
    sid <- sites$site_id[i]
    b   <- sites_band[i]
    pool <- cand_df[cand_df$band == b, , drop = FALSE]
    if (nrow(pool) == 0) {
      warning(sprintf("match_controls: no candidates in band %s for site_id %s", b, sid))
      next
    }
    site_v <- as.numeric(sites_vals[i, ])
    pool_mat <- as.matrix(pool[, cov_names, drop = FALSE])
    d <- sqrt(rowSums(sweep(pool_mat, 2, site_v, "-")^2))
    ord <- order(d)[seq_len(min(n, nrow(pool)))]
    res <- st_as_sf(pool[ord, ], coords = c("x", "y"), crs = crs(band_rast))
    res$site_id <- sid
    res$band <- b
    res$rank <- seq_along(ord)
    res$cov_distance <- d[ord]
    out[[length(out) + 1]] <- res[, c("site_id", "band", "rank", "cov_distance")]
  }
  if (length(out) == 0) return(sites[0, c("site_id"), drop = FALSE])
  do.call(rbind, out)
}

#' Check covariate balance between intervention sites and their matched controls
#'
#' Standardised mean difference (SMD) per covariate: (mean_sites - mean_controls)
#' / pooled SD. |SMD| > 0.1 is flagged as imbalanced (Austin 2009 convention).
#'
#' @param sites sf POINT, intervention sites.
#' @param controls sf POINT, output of match_controls().
#' @param z_stack SpatRaster used for matching -- balance is reported in this
#'   same z-scored space, so SMDs are directly comparable to the 0.1 threshold.
#' @param diag_stack optional SpatRaster of covariates NOT used in matching
#'   (e.g. the excluded biodiversity layers), reported separately as
#'   diagnostic-only balance -- imbalance here does not mean matching failed,
#'   since these were deliberately excluded from banding.
#' @return list(matching = data.frame, diagnostic = data.frame or NULL), each
#'   with covariate, mean_sites, mean_controls, smd, imbalanced.
match_balance <- function(sites, controls, z_stack, diag_stack = NULL) {
  smd_table <- function(stack) {
    site_vals <- terra::extract(stack, vect(sites))[, -1, drop = FALSE]
    ctrl_vals <- terra::extract(stack, vect(controls))[, -1, drop = FALSE]
    covs <- names(stack)
    res <- lapply(covs, function(cv) {
      a <- site_vals[[cv]]; b <- ctrl_vals[[cv]]
      pooled_sd <- sqrt((var(a, na.rm = TRUE) + var(b, na.rm = TRUE)) / 2)
      smd <- if (is.na(pooled_sd) || pooled_sd == 0) NA_real_ else
        (mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE)) / pooled_sd
      data.frame(covariate = cv, mean_sites = mean(a, na.rm = TRUE),
                 mean_controls = mean(b, na.rm = TRUE), smd = smd)
    })
    df <- bind_rows(res)
    df$imbalanced <- abs(df$smd) > 0.1
    df
  }

  matching_balance <- smd_table(z_stack[[setdiff(names(z_stack), "band")]])
  if (any(matching_balance$imbalanced, na.rm = TRUE)) {
    cat("WARNING: imbalanced covariates (|SMD| > 0.1) in matching set:\n")
    print(matching_balance[which(matching_balance$imbalanced), c("covariate", "smd")])
  }

  diagnostic_balance <- NULL
  if (!is.null(diag_stack)) {
    diagnostic_balance <- smd_table(diag_stack)
    diagnostic_balance$note <- "diagnostic only -- not used for band assignment"
  }

  list(matching = matching_balance, diagnostic = diagnostic_balance)
}

#' Rank intervention sites within each band by the Layer 4 monitoring score
#'
#' CAVEATS: the monitoring score (build_monitoring_score(), 07_plots.R:3038) is
#' multiplied by the restoration opportunity layer and so is masked to the
#' opportunity footprint -- it can rank candidate intervention sites but
#' CANNOT be used to select controls (controls by design sit outside the
#' footprint or are not being restored). Its diversity and similarity inputs
#' come from species records, so it inherits recorder-effort bias:
#' well-surveyed places score higher regardless of true value.
#'
#' @param sites sf POINT, one row per site. site_id added if absent.
#' @param band_rast SpatRaster, single band ("band").
#' @param monitoring_score_rast SpatRaster or path, e.g.
#'   outputs/habitats/seagrass_monitoring_score.tif.
#' @return sf sites with band, monitoring_score, rank_in_band (1 = highest
#'   score in that band; NA where the score is NA) and n_in_band.
rank_sites_within_band <- function(sites, band_rast, monitoring_score_rast) {
  stopifnot(inherits(sites, "sf"))
  if (is.character(monitoring_score_rast)) monitoring_score_rast <- rast(monitoring_score_rast)
  if (!"site_id" %in% names(sites)) sites$site_id <- seq_len(nrow(sites))

  sites$band <- terra::extract(band_rast, vect(st_transform(sites, crs(band_rast))))[, 2]
  sites$monitoring_score <- terra::extract(
    monitoring_score_rast, vect(st_transform(sites, crs(monitoring_score_rast))))[, 2]

  sites %>%
    group_by(band) %>%
    mutate(rank_in_band = rank(-monitoring_score, ties.method = "min", na.last = "keep"),
           n_in_band = n()) %>%
    ungroup() %>%
    arrange(band, rank_in_band)
}

#' Spatially balanced monitoring-station placement within a site polygon
#'
#' GRTS (generalized random tessellation stratified) sample via spsurvey.
#' Over-sample points are replacements, used in order, if base stations
#' prove inaccessible.
#'
#' @param site_poly sf POLYGON/MULTIPOLYGON, the site (or control) footprint.
#' @param n_base integer, number of base stations.
#' @param n_over integer, number of over-sample (replacement) stations.
#' @return sf POINT, base then over-sample stations (spsurvey's `siteuse`
#'   column distinguishes them).
place_stations <- function(site_poly, n_base, n_over = 0) {
  if (!requireNamespace("spsurvey", quietly = TRUE)) {
    stop("place_stations: package 'spsurvey' is required -- install.packages('spsurvey')")
  }
  stopifnot(inherits(site_poly, "sf"))
  site_poly <- st_make_valid(st_transform(site_poly, 27700))
  set.seed(SEED)
  design <- spsurvey::grts(sframe = site_poly, n_base = n_base,
                           n_over = if (n_over > 0) n_over else NULL)
  out <- design$sites_base
  if (!is.null(design$sites_over)) out <- rbind(out, design$sites_over)
  out
}

#' Representativeness of AOI polygons relative to the whole programme area
#'
#' For each AOI polygon and for the whole programme area: band composition
#' (% of cells per band) and median / IQR of each continuous banding covariate
#' in the cube. The programme area is a regular sample (~sample_n cells inside
#' mask_r); each AOI is the cells inside that polygon, subsampled to sample_n
#' if larger. The full raster is never read into memory.
#'
#' @param aoi_polys sf polygons with aoi_name (see load_aoi_bays()).
#' @param cube_stack SpatRaster cube (layer names = covariate names).
#' @param band_rast SpatRaster band layer on the same grid.
#' @param mask_r programme-area mask on the same grid.
#' @param label "marine" / "terrestrial" -- used in file names.
#' @param sample_n target cells per area (REPR_SAMPLE_N); AOIs larger than
#'   this are randomly subsampled to it.
#' @return invisibly, list(table, samples).
compare_aoi_to_region <- function(aoi_polys, cube_stack, band_rast, mask_r, label,
                                  sample_n = REPR_SAMPLE_N) {
  region_name <- "Programme area"
  cont_covs <- covariate_table %>% filter(type == "continuous", use_for_banding) %>% pull(name)
  cont_covs <- intersect(cont_covs, names(cube_stack))
  stk <- c(mask_r, band_rast, cube_stack[[cont_covs]])
  names(stk) <- c("in_area", "band", cont_covs)

  ## Whole programme area: regular sample sized so ~sample_n land inside the mask
  frac <- global(!is.na(mask_r), "mean", na.rm = TRUE)[[1]]
  size <- min(ncell(stk), ceiling(sample_n / max(frac, 1e-6)))
  reg <- spatSample(stk, size = size, method = "regular", as.df = TRUE, na.rm = FALSE)
  reg <- reg[!is.na(reg$in_area), , drop = FALSE]
  reg$area <- region_name

  ## Each AOI polygon: cells inside the polygon (and inside the programme area)
  aoi_v <- vect(st_transform(aoi_polys, crs(band_rast)))
  aoi_samples <- list()
  for (i in seq_len(nrow(aoi_polys))) {
    nm <- aoi_polys$aoi_name[i]
    vals <- terra::extract(stk, aoi_v[i], ID = FALSE)
    vals <- vals[!is.na(vals$in_area), , drop = FALSE]
    if (nrow(vals) == 0) {
      cat(sprintf("    %s: no %s programme-area cells -- skipped\n", nm, label))
      next
    }
    if (nrow(vals) > sample_n) vals <- vals[sample.int(nrow(vals), sample_n), , drop = FALSE]
    vals$area <- nm
    aoi_samples[[nm]] <- vals
  }
  all_df <- bind_rows(c(list(reg), aoi_samples))
  areas <- c(region_name, names(aoi_samples))

  ## Band composition, with 0% rows for bands absent from an area
  counts <- all_df %>% filter(!is.na(band)) %>% count(area, band, name = "n_cells")
  grid <- expand.grid(area = areas, band = sort(unique(counts$band)), stringsAsFactors = FALSE)
  band_comp <- merge(grid, counts, all.x = TRUE) %>%
    mutate(n_cells = ifelse(is.na(n_cells), 0L, n_cells)) %>%
    group_by(area) %>%
    mutate(pct_cells = 100 * n_cells / sum(n_cells)) %>%
    ungroup() %>%
    mutate(measure = "band_composition", covariate = "band")

  ## Continuous covariate summaries
  long <- bind_rows(lapply(cont_covs, function(cv)
    data.frame(area = all_df$area, covariate = cv, value = all_df[[cv]], stringsAsFactors = FALSE))) %>%
    filter(!is.na(value))
  cov_stats <- long %>%
    group_by(area, covariate) %>%
    summarise(n_cells = n(), median = median(value),
              q25 = unname(quantile(value, 0.25)), q75 = unname(quantile(value, 0.75)),
              .groups = "drop") %>%
    mutate(iqr = q75 - q25, measure = "covariate_summary")

  out <- bind_rows(band_comp, cov_stats) %>%
    mutate(area = factor(area, levels = areas)) %>%
    arrange(measure, area, covariate, band) %>%
    mutate(area = as.character(area)) %>%
    select(measure, area, covariate, band, n_cells, pct_cells, median, q25, q75, iqr)
  csv_path <- sprintf("outputs/design/representativeness_%s.csv", label)
  readr_write(out, csv_path)

  long$area <- factor(long$area, levels = areas)
  p <- ggplot(mapping = aes(x = value)) +
    geom_density(data = filter(long, area == region_name), fill = "grey85", colour = "grey35") +
    geom_density(data = filter(long, area != region_name), aes(colour = area), linewidth = 0.6) +
    facet_wrap(~covariate, scales = "free") +
    labs(title = sprintf("Representativeness of AOI bays -- %s banding covariates", label),
         subtitle = "Grey fill = whole programme area (regular sample); lines = AOI polygons",
         x = NULL, y = "Density", colour = "AOI") +
    theme_minimal()
  png_path <- sprintf("outputs/design/representativeness_%s.png", label)
  ggsave(png_path, p, width = 12, height = 8, dpi = 300)

  cat(sprintf("  %s band composition (%% of cells):\n", label))
  print(as.data.frame(band_comp %>% select(area, band, pct_cells) %>%
                        mutate(pct_cells = round(pct_cells, 1))))
  cat(sprintf("  Written: %s\n  Written: %s\n", csv_path, png_path))
  invisible(list(table = out, samples = long))
}

if (is.na(K_MARINE) || is.na(K_TERR)) {
  cat("\n*** PASS 1 COMPLETE. K_MARINE / K_TERR not set -- stopping before band assignment. ***\n")
} else {

##=====================================================================
## 11. BAND ASSIGNMENT (chunked via terra::predict)
##=====================================================================

nearest_centroid_fun <- function(model, data, ...) {
  ## Column order of `data` follows the SpatRaster's own layer order, which is
  ## not guaranteed to match model$centers' column order (built from the
  ## sampled training data.frame). Realign by name, not position.
  data <- data[, colnames(model$centers), drop = FALSE]
  d <- as.matrix(data)
  out <- rep(NA_integer_, nrow(d))
  ok <- stats::complete.cases(d)
  if (any(ok)) {
    dmat <- d[ok, , drop = FALSE]
    dist2 <- vapply(seq_len(nrow(model$centers)), function(k) {
      rowSums(sweep(dmat, 2, model$centers[k, ], "-")^2)
    }, numeric(nrow(dmat)))
    out[ok] <- max.col(-dist2, ties.method = "first")
  }
  out
}

## Same transform rule as the training sample (see section 9): continuous
## layers z-scored with the sample mean/sd; one-hot layers = CAT_WEIGHT *
## indicator, built by matching the raster's category LABELS to onehot_map
## (not segregate() output order). A training level absent from the raster
## becomes a zero layer (NA where the categorical layer is NA). Output layers
## are ordered as `layers` (the training column order).
make_z_stack <- function(stack, banding_names, cat_cols, scaling, onehot_map, layers, mask_r) {
  s <- mask(stack[[banding_names]], mask_r)
  z_layers <- list()
  for (nm in names(scaling$means)) {
    z_layers[[nm]] <- (s[[nm]] - scaling$means[[nm]]) / scaling$sds[[nm]]
  }
  for (cl in cat_cols) {
    cc <- cat_codes(s[[cl]])
    r_num <- s[[cl]]
    levels(r_num) <- NULL
    cl_map <- onehot_map[onehot_map$covariate == cl, , drop = FALSE]
    for (j in seq_len(nrow(cl_map))) {
      codes <- cc$code[which(cc$label == cl_map$level[j])]
      if (length(codes) == 0) {
        cat(sprintf("    %s: level '%s' not in raster -- zero layer\n", cl, cl_map$level[j]))
        ind <- r_num * 0
      } else {
        ind <- r_num == codes[1]
        for (cd in codes[-1]) ind <- ind | (r_num == cd)
      }
      z_layers[[cl_map$column[j]]] <- CAT_WEIGHT * ind
    }
  }
  missing_cols <- setdiff(layers, names(z_layers))
  if (length(missing_cols) > 0) stop("make_z_stack: no layer built for: ", paste(missing_cols, collapse = ", "))
  z <- rast(z_layers[layers])
  names(z) <- layers
  z
}

check_z_names <- function(z_stack, centers, label) {
  if (!identical(names(z_stack), colnames(centers))) {
    stop(sprintf(paste0("%s z_stack layer names do not match kmeans centre columns.\n",
                        "  z_stack: %s\n  centers: %s"),
                 label, paste(names(z_stack), collapse = ", "),
                 paste(colnames(centers), collapse = ", ")))
  }
  invisible(TRUE)
}

cat("\n-- Assigning bands: marine --\n")
km_m <- kmeans(z_samp_m, centers = K_MARINE, nstart = 25, iter.max = 100)
z_stack_m <- make_z_stack(marine_stack, banding_names_m, cat_cols_m, scaling_m,
                          onehot_map_m, colnames(z_samp_m), marine_mask_r)
check_z_names(z_stack_m, km_m$centers, "Marine")
marine_bands <- predict(z_stack_m, model = list(centers = km_m$centers), fun = nearest_centroid_fun,
                         na.rm = FALSE)
names(marine_bands) <- "band"
marine_bands <- write_cog(marine_bands, "outputs/design/marine_bands.tif", datatype = "INT1U")
names(marine_bands) <- "band"

cat("-- Assigning bands: terrestrial --\n")
km_t <- kmeans(z_samp_t, centers = K_TERR, nstart = 25, iter.max = 100)
z_stack_t <- make_z_stack(terr_stack, banding_names_t, cat_cols_t, scaling_t,
                          onehot_map_t, colnames(z_samp_t), land_mask_r)
check_z_names(z_stack_t, km_t$centers, "Terrestrial")
terr_bands <- predict(z_stack_t, model = list(centers = km_t$centers), fun = nearest_centroid_fun,
                       na.rm = FALSE)
names(terr_bands) <- "band"
terr_bands <- write_cog(terr_bands, "outputs/design/terrestrial_bands.tif", datatype = "INT1U")
names(terr_bands) <- "band"

saveRDS(list(centers = km_m$centers, scaling = scaling_m, layers = colnames(z_samp_m),
             onehot_cols = onehot_map_m$column, onehot_map = onehot_map_m, cat_weight = CAT_WEIGHT),
        "outputs/design/band_model_marine.rds")
saveRDS(list(centers = km_t$centers, scaling = scaling_t, layers = colnames(z_samp_t),
             onehot_cols = onehot_map_t$column, onehot_map = onehot_map_t, cat_weight = CAT_WEIGHT),
        "outputs/design/band_model_terrestrial.rds")

##=====================================================================
## 12. BAND POLYGONS + PROFILES
##=====================================================================

cat("\n-- Writing band polygons and profiles --\n")

bands_to_polygons <- function(band_r, out_path) {
  p <- as.polygons(band_r, dissolve = TRUE)
  st_write(st_as_sf(p), out_path, delete_dsn = TRUE, quiet = TRUE)
}
bands_to_polygons(marine_bands, "outputs/design/marine_bands.gpkg")
bands_to_polygons(terr_bands, "outputs/design/terrestrial_bands.gpkg")

band_profile <- function(band_r, cov_stack, cat_cols) {
  full <- c(band_r, cov_stack)
  df <- as.data.frame(full, na.rm = TRUE)
  cont_cols <- setdiff(names(cov_stack), cat_cols)
  prof <- df %>% group_by(band) %>%
    summarise(across(all_of(cont_cols), ~mean(.x, na.rm = TRUE)), n = n(), .groups = "drop")
  for (cl in cat_cols) {
    mode_tab <- df %>% group_by(band) %>%
      summarise(!!paste0(cl, "_mode") := names(sort(table(.data[[cl]]), decreasing = TRUE))[1],
                .groups = "drop")
    prof <- left_join(prof, mode_tab, by = "band")
  }
  prof
}

profile_m <- band_profile(marine_bands, marine_stack[[banding_names_m]], cat_cols_m)
profile_t <- band_profile(terr_bands, terr_stack[[banding_names_t]], cat_cols_t)
readr_write(profile_m, "outputs/design/marine_band_profiles.csv")
readr_write(profile_t, "outputs/design/terrestrial_band_profiles.csv")

##=====================================================================
## 13. CUBES
##=====================================================================

cat("\n-- Assembling data cubes --\n")

build_cube <- function(cont_stack, cat_layers, extra_layers, band_r, out_path, bands_csv_path) {
  layer_list <- list()
  band_meta <- list()

  for (nm in names(cont_stack)) {
    layer_list[[nm]] <- cont_stack[[nm]]
    band_meta[[length(band_meta) + 1]] <- data.frame(name = nm, type = "continuous")
  }
  for (nm in names(cat_layers)) {
    layer_list[[nm]] <- cat_layers[[nm]]
    band_meta[[length(band_meta) + 1]] <- data.frame(name = nm, type = "categorical")
  }
  for (nm in names(extra_layers)) {
    layer_list[[nm]] <- extra_layers[[nm]]
    band_meta[[length(band_meta) + 1]] <- data.frame(name = nm, type = "cube_only")
  }
  layer_list[["band"]] <- band_r
  band_meta[[length(band_meta) + 1]] <- data.frame(name = "band", type = "band")

  cube <- rast(layer_list)
  names(cube) <- names(layer_list)
  ## A COG carries one data type across all bands. The cube mixes negative
  ## continuous values (depth, aspect) with small integer class codes, and
  ## letting terra choose produced an INT1U file in which every negative
  ## value was lost to nodata. Forced to float.
  write_cog(cube, out_path, datatype = "FLT4S")

  meta <- bind_rows(band_meta)
  meta$band_number <- seq_len(nrow(meta))
  meta <- meta %>%
    left_join(covariate_table %>% select(name, units, resampling, gap_fill, use_for_banding, description),
              by = "name") %>%
    mutate(gap_filled = ifelse(is.na(gap_fill), FALSE, gap_fill),
           gap_fill_method = ifelse(gap_filled, "iterative 3x3 focal mean (na.policy='only')", "none")) %>%
    select(band_number, name, type, units, resampling, gap_filled, gap_fill_method, use_for_banding, description)
  readr_write(meta, bands_csv_path)
  cube
}

## NB: [[ ]], not [ ], for layer subsetting -- [ ] extracts cell values on a SpatRaster.
marine_cat_names <- intersect(cat_rows_m$name, names(marine_stack))
marine_cat_layers <- marine_stack[[marine_cat_names]]
marine_cont_layers <- marine_stack[[setdiff(names(marine_stack), marine_cat_names)]]

marine_cube <- build_cube(
  cont_stack = marine_cont_layers,
  cat_layers = marine_cat_layers,
  extra_layers = list(n_gapfilled = n_gapfilled),
  band_r = marine_bands,
  out_path = "outputs/design/SCAMP_marine_cube_286m.tif",
  bands_csv_path = "outputs/design/SCAMP_marine_cube_286m_bands.csv"
)

terr_cat_names <- intersect(c("land_cover"), names(terr_stack))
terr_cat_layers <- terr_stack[[terr_cat_names]]
terr_cont_layers <- terr_stack[[setdiff(names(terr_stack), terr_cat_names)]]

terr_extra <- list()
if (!is.null(dc_class_r)) terr_extra[["erosion_accretion_class"]] <- dc_class_r

terr_cube <- build_cube(
  cont_stack = terr_cont_layers,
  cat_layers = terr_cat_layers,
  extra_layers = terr_extra,
  band_r = terr_bands,
  out_path = "outputs/design/SCAMP_terrestrial_cube_100m.tif",
  bands_csv_path = "outputs/design/SCAMP_terrestrial_cube_100m_bands.csv"
)

## Categorical lookup CSVs alongside the cubes -- one per categorical cube
## band, taken from cats() of the final (post-collapse) in-memory cube layer so
## the codes match the cube exactly. Extra descriptive columns are joined on
## by label where available.
write_cat_lookup <- function(cube, nm, out_path, extra = NULL, extra_by = NULL) {
  if (!nm %in% names(cube)) return(invisible(NULL))
  lk <- cat_codes(cube[[nm]])
  if (!is.null(extra)) {
    lk <- left_join(lk, extra, by = c("label" = extra_by))
    ## The collapsed rare-level bin has no row in the source lookup; label it
    ## in the first descriptive column so it doesn't ship blank.
    desc_col <- setdiff(names(extra), extra_by)[1]
    if (!is.na(desc_col)) {
      lk[[desc_col]][which(lk$label == "other")] <- sprintf(
        "Classes individually covering less than %g%% of the area, combined",
        RARE_LEVEL_MIN_PROP * 100)
    }
  }
  readr_write(lk, out_path)
  cat(sprintf("  Lookup: %s (%d codes)\n", out_path, nrow(lk)))
}

cat("\n-- Writing categorical lookups --\n")
for (nm in marine_cat_names) {
  write_cat_lookup(marine_cube, nm, sprintf("outputs/design/SCAMP_marine_cube_286m_%s_lookup.csv", nm))
}
write_cat_lookup(terr_cube, "land_cover",
                 "outputs/design/SCAMP_terrestrial_cube_100m_land_cover_lookup.csv",
                 extra = land_cover_lookup %>% select(eunis_code, eunis_desc), extra_by = "eunis_code")
if (!is.null(dc_class_r)) {
  dc_rules <- data.frame(
    class = c("eroding", "stable", "accreting"),
    rule  = c(sprintf("HIST_RATE < -%.2f m/yr", DC_STABLE_THRESHOLD_MYR),
              sprintf("-%.2f <= HIST_RATE <= %.2f m/yr", DC_STABLE_THRESHOLD_MYR, DC_STABLE_THRESHOLD_MYR),
              sprintf("HIST_RATE > %.2f m/yr", DC_STABLE_THRESHOLD_MYR)))
  write_cat_lookup(terr_cube, "erosion_accretion_class",
                   "outputs/design/SCAMP_terrestrial_cube_100m_erosion_accretion_class_lookup.csv",
                   extra = dc_rules, extra_by = "class")
}

cat("\nCubes written:\n")
cat("  outputs/design/SCAMP_marine_cube_286m.tif\n")
cat("  outputs/design/SCAMP_terrestrial_cube_100m.tif\n")

##=====================================================================
## 13b. REPRESENTATIVENESS — AOI bays vs whole programme area
##=====================================================================

cat("\n-- Representativeness: AOI bays vs programme area --\n")
aoi_bays <- load_aoi_bays()
repr_m <- compare_aoi_to_region(aoi_bays, marine_cube, marine_bands, marine_mask_r, "marine")
repr_t <- compare_aoi_to_region(aoi_bays, terr_cube, terr_bands, land_mask_r, "terrestrial")

} # end K_MARINE/K_TERR gate

##=====================================================================
## 14. WORKED EXAMPLE — controls for stand-in intervention sites
##=====================================================================
## Not run by default (needs pass-2 objects marine_bands, z_stack_m). AOI bay
## centroids stand in for intervention sites until real sites are confirmed.

if (FALSE) {
  aoi_bays <- load_aoi_bays()
  example_sites <- st_sf(site_id = aoi_bays$aoi_name,
                         geometry = st_centroid(st_geometry(aoi_bays)))
  site_band <- terra::extract(marine_bands, vect(example_sites))[, 2]
  if (any(is.na(site_band))) {
    cat(sprintf("Dropping centroid(s) outside the banded marine area: %s\n",
                paste(example_sites$site_id[is.na(site_band)], collapse = ", ")))
    example_sites <- example_sites[!is.na(site_band), ]
  }

  example_controls <- match_controls(example_sites, marine_bands, z_stack_m, n = 5, min_dist = 2000)
  example_balance  <- match_balance(example_sites, example_controls, z_stack_m)
  print(example_balance$matching)

  st_write(example_controls, "outputs/design/worked_example_controls.gpkg",
           layer = "controls", delete_dsn = TRUE, quiet = TRUE)
  st_write(example_sites, "outputs/design/worked_example_controls.gpkg",
           layer = "sites", quiet = TRUE)
  cat("Written: outputs/design/worked_example_controls.gpkg (layers: controls, sites)\n")
}

cat("\n========== 08_environmental_bands.R complete ==========\n")

plot_bands <- function(band_path, label, coast_path = "data/processed/boundaries/scamp_coastline.gpkg",
                       aoi_path = "data/processed/boundaries/aois_named.gpkg") {
  b <- rast(band_path)
  cell_km2 <- prod(res(b)) / 1e6

  # coverage table: cells, km2 and % per band
  cov <- freq(b) |>
    transmute(band = value, cells = count,
              km2 = round(count * cell_km2, 1),
              pct = round(100 * count / sum(count), 1))
  print(cov)
  write.csv(cov, sprintf("outputs/design/%s_band_coverage.csv", label), row.names = FALSE)

  b_f <- as.factor(b)
  levels(b_f) <- data.frame(id = cov$band,
                            band = sprintf("Band %d (%.0f%%)", cov$band, cov$pct))
  coast <- st_read(coast_path, quiet = TRUE)
  aois  <- st_read(aoi_path, quiet = TRUE)

  p <- ggplot() +
    geom_spatraster(data = b_f) +
    scale_fill_brewer(palette = "Set2", na.value = NA, na.translate = FALSE, name = NULL) +
    geom_sf(data = coast, colour = "grey30", linewidth = 0.2) +
    geom_sf(data = aois, fill = NA, colour = "black", linewidth = 0.4, linetype = "dashed") +
    geom_sf_text(data = aois, aes(label = bay_name), size = 2.5, nudge_y = 2500) +
    labs(title = sprintf("%s environmental bands (k = %d)", tools::toTitleCase(label), nrow(cov)),
         subtitle = "Dashed boxes: candidate pilot bays. Percent = share of banded area.") +
    theme_minimal() + theme(axis.title = element_blank())

  ggsave(sprintf("outputs/design/%s_bands_map.png", label), p, width = 10, height = 7, dpi = 300)
  p
}

plot_bands("outputs/design/marine_bands.tif", "marine")
plot_bands("outputs/design/terrestrial_bands.tif", "terrestrial")


r <- rast("outputs/design/marine_bands.tif")
r                                    # look at min/max and the layer name
global(!is.na(r), "sum", na.rm = TRUE)   # how many non-NA cells
freq(r)

# compare with what is still in memory
marine_bands
global(!is.na(marine_bands), "sum", na.rm = TRUE)
freq(marine_bands)

global(!is.na(z_stack_m), "sum", na.rm = TRUE)   # non-NA count per layer
scaling_m                     # look for any sd of 0 or NA
freq(terr_bands)              # did terrestrial survive?
r <- marine_stack[["energy_class"]]; levels(r) <- NULL; freq(r == 1)
