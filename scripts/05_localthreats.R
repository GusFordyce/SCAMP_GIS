##=====================================================================
## SCAMP GIS Analysis — Local Threats & Pressures
## 05_localthreats.R
##=====================================================================
Sys.setenv(GFW_TOKEN = "eyJhbGciOiJSUzI1NiIsInR5cCI6IkpXVCIsImtpZCI6ImtpZEtleSJ9.eyJkYXRhIjp7Im5hbWUiOiJTQ0FNUCIsInVzZXJJZCI6NjgxMTEsImFwcGxpY2F0aW9uTmFtZSI6IlNDQU1QIiwiaWQiOjEzMjE0LCJ0eXBlIjoidXNlci1hcHBsaWNhdGlvbiJ9LCJpYXQiOjE3ODY0NTUyNzEsImV4cCI6MjEwMTgxNTI3MSwiYXVkIjoiZ2Z3IiwiaXNzIjoiZ2Z3In0.nrATPYq193xBmlfQIe5VYxxWQ-ugBI04BLOCDYbFBIHP1x0WEaX770cWkS1zfDtmKag1OtTbpR0AwNdm9S40uuwgtDQKUy-w_5sMjWcOBu2yG9iHeQT83dyD1wap_Z4oYfhS9HVpYsW53LJFCaxdKUCPnrxJk4nrbWNqNudEkEcW6GErGsSCtSQF8OBimPpqme9rO0sh6g9Ch9KUwwTxHGg2SYWOT66Pe2EECxeibQ7cNmgbpblkKYnH9YwFnvNWmyUzidR_CJnugHFtX3hOnu2ESjsr8RYi4YxrjhskQfq7yizMgj3UoiTEkxLyWkKZoojHtk-cG-GB9IqyEpvOQMWjF80ijCll2d8fQ453hIN95Jptv85ImxbsIGYnUf_k2iPi6bdTxsMoXDUGf7b_lUHJc7hvsgXr5gYgcw7ZULQEySqhCb91IrlZ5sAMinX5jqPdb_jzeh6r4q-b7L3z6qwF1ZyQBl6YwuudfrSiMbEIQXNm8Rz5RqTDMkDZWdCf")
require(sf)
require(terra)
require(dplyr)
require(ggplot2)
require(ggspatial)
require(patchwork)
require(mapview)
require(whitebox)
##=====================================================================
## 1. SETUP
##=====================================================================

cat("\n===== SETUP =====\n")

Sys.setenv(SHAPE_RESTORE_SHX = "YES")

scamp_aoi_sf <- st_read("data/processed/boundaries/scamp_aoi.gpkg", quiet = TRUE)
cat(sprintf("  SCAMP AOI loaded: CRS %s\n", st_crs(scamp_aoi_sf)$input))

smr_solway_bng <- st_read(
  "data/raw/boundaries/administrative_units_scottish_marine_regions/administrative_units_scottish_marine_regionsPolygon.shp",
  quiet = TRUE
) %>% st_transform(27700)
smr_solway_bng <- smr_solway_bng[grepl("Solway", smr_solway_bng$objnam, ignore.case = TRUE), ]
cat(sprintf("  SMR Solway: %d features\n", nrow(smr_solway_bng)))

las_dg_bng <- st_read(
  "data/raw/boundaries/Local_Authority_Boundaries_-_Scotland/pub_las.shp",
  quiet = TRUE
) %>% st_transform(27700)
las_dg_bng <- las_dg_bng[grepl("dumfries", las_dg_bng$local_auth, ignore.case = TRUE), ]
cat(sprintf("  D&G boundary: %d features\n", nrow(las_dg_bng)))

scamp_coastline <- st_read("data/processed/boundaries/scamp_coastline.gpkg", quiet = TRUE)
cat(sprintf("  Coastline: %d features\n", nrow(scamp_coastline)))

map_theme <- theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text  = element_text(size = 7),
        plot.title = element_text(face = "bold"))

dir.create("outputs/threats", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/threats", showWarnings = FALSE, recursive = TRUE)

aoi_bbox <- st_bbox(scamp_aoi_sf)

## --- Clip designated area boundaries to SCAMP AOI ---

cat("\n===== DESIGNATED AREAS =====\n")

boundary_gpkg <- list(
  ibas            = "data/raw/boundaries/IBAs_UK_-8129902679595446195.gpkg",
  rspb_reserves   = "data/raw/boundaries/RSPB_Public_Reserves_7326421649958477753.gpkg",
  mpas            = "data/raw/boundaries/Marine_Protected_Areas_-6371676791122909031.gpkg",
  ramsar          = "data/raw/boundaries/RAMSAR_Wetlands_of_International_Importance_-1649623871716453965.gpkg",
  sssi            = "data/raw/boundaries/Sites_of_Special_Scientific_Interest_8447690424398606363.gpkg",
  sac             = "data/raw/boundaries/Special_Areas_of_Conservation_5743281901804056965.gpkg",
  nnr             = "data/raw/boundaries/National_Nature_Reserves_1039238937738926456.gpkg",
  spa             = "data/raw/boundaries/Special_Protection_Areas_4781773662199779435.gpkg",
  nature_reserves = "data/raw/boundaries/Nature_Reserves_1930300340402281996.gpkg",
  lnr             = "data/raw/boundaries/Local_Nature_Reserves_7017289517571999860.gpkg",
  harbour_limits  = "data/raw/boundaries/area_management_harbour_limits/area_management_harbour_limitsPolygon.shp"
)

scamp_boundaries <- list()

# Read IBAs
cat("Reading ibas...\n")
raw <- st_read(boundary_gpkg[["ibas"]], quiet = TRUE) %>%
  st_zm(drop = TRUE, what = "ZM") %>% st_transform(27700) %>% st_make_valid()
hits <- lengths(st_intersects(raw, scamp_aoi_sf)) > 0
scamp_boundaries[["ibas"]] <- st_intersection(raw[hits, ], scamp_aoi_sf)
cat(sprintf("  ibas: %d features\n", nrow(scamp_boundaries[["ibas"]])))

# Read RSPB reserves
cat("Reading rspb_reserves...\n")
raw <- st_read(boundary_gpkg[["rspb_reserves"]], quiet = TRUE) %>%
  st_zm(drop = TRUE, what = "ZM") %>% st_transform(27700) %>% st_make_valid()
hits <- lengths(st_intersects(raw, scamp_aoi_sf)) > 0
scamp_boundaries[["rspb_reserves"]] <- st_intersection(raw[hits, ], scamp_aoi_sf)
cat(sprintf("  rspb_reserves: %d features\n", nrow(scamp_boundaries[["rspb_reserves"]])))

# Read MPAs
cat("Reading mpas...\n")
raw <- st_read(boundary_gpkg[["mpas"]], quiet = TRUE) %>%
  st_zm(drop = TRUE, what = "ZM") %>% st_transform(27700) %>% st_make_valid()
hits <- lengths(st_intersects(raw, scamp_aoi_sf)) > 0
scamp_boundaries[["mpas"]] <- st_intersection(raw[hits, ], scamp_aoi_sf)
cat(sprintf("  mpas: %d features\n", nrow(scamp_boundaries[["mpas"]])))

# Read Ramsar
cat("Reading ramsar...\n")
raw <- st_read(boundary_gpkg[["ramsar"]], quiet = TRUE) %>%
  st_zm(drop = TRUE, what = "ZM") %>% st_transform(27700) %>% st_make_valid()
hits <- lengths(st_intersects(raw, scamp_aoi_sf)) > 0
scamp_boundaries[["ramsar"]] <- st_intersection(raw[hits, ], scamp_aoi_sf)
cat(sprintf("  ramsar: %d features\n", nrow(scamp_boundaries[["ramsar"]])))

# Read SSSI
cat("Reading sssi...\n")
raw <- st_read(boundary_gpkg[["sssi"]], quiet = TRUE) %>%
  st_zm(drop = TRUE, what = "ZM") %>% st_transform(27700) %>% st_make_valid()
hits <- lengths(st_intersects(raw, scamp_aoi_sf)) > 0
scamp_boundaries[["sssi"]] <- st_intersection(raw[hits, ], scamp_aoi_sf)
cat(sprintf("  sssi: %d features\n", nrow(scamp_boundaries[["sssi"]])))

# Read SAC
cat("Reading sac...\n")
raw <- st_read(boundary_gpkg[["sac"]], quiet = TRUE) %>%
  st_zm(drop = TRUE, what = "ZM") %>% st_transform(27700) %>% st_make_valid()
hits <- lengths(st_intersects(raw, scamp_aoi_sf)) > 0
scamp_boundaries[["sac"]] <- st_intersection(raw[hits, ], scamp_aoi_sf)
cat(sprintf("  sac: %d features\n", nrow(scamp_boundaries[["sac"]])))

# Read NNR
cat("Reading nnr...\n")
raw <- st_read(boundary_gpkg[["nnr"]], quiet = TRUE) %>%
  st_zm(drop = TRUE, what = "ZM") %>% st_transform(27700) %>% st_make_valid()
hits <- lengths(st_intersects(raw, scamp_aoi_sf)) > 0
scamp_boundaries[["nnr"]] <- st_intersection(raw[hits, ], scamp_aoi_sf)
cat(sprintf("  nnr: %d features\n", nrow(scamp_boundaries[["nnr"]])))

# Read SPA
cat("Reading spa...\n")
raw <- st_read(boundary_gpkg[["spa"]], quiet = TRUE) %>%
  st_zm(drop = TRUE, what = "ZM") %>% st_transform(27700) %>% st_make_valid()
hits <- lengths(st_intersects(raw, scamp_aoi_sf)) > 0
scamp_boundaries[["spa"]] <- st_intersection(raw[hits, ], scamp_aoi_sf)
cat(sprintf("  spa: %d features\n", nrow(scamp_boundaries[["spa"]])))

# Read Nature Reserves
cat("Reading nature_reserves...\n")
raw <- st_read(boundary_gpkg[["nature_reserves"]], quiet = TRUE) %>%
  st_zm(drop = TRUE, what = "ZM") %>% st_transform(27700) %>% st_make_valid()
hits <- lengths(st_intersects(raw, scamp_aoi_sf)) > 0
scamp_boundaries[["nature_reserves"]] <- st_intersection(raw[hits, ], scamp_aoi_sf)
cat(sprintf("  nature_reserves: %d features\n", nrow(scamp_boundaries[["nature_reserves"]])))

# Read LNR
cat("Reading lnr...\n")
raw <- st_read(boundary_gpkg[["lnr"]], quiet = TRUE) %>%
  st_zm(drop = TRUE, what = "ZM") %>% st_transform(27700) %>% st_make_valid()
hits <- lengths(st_intersects(raw, scamp_aoi_sf)) > 0
scamp_boundaries[["lnr"]] <- st_intersection(raw[hits, ], scamp_aoi_sf)
cat(sprintf("  lnr: %d features\n", nrow(scamp_boundaries[["lnr"]])))

# Read Harbour Limits
cat("Reading harbour_limits...\n")
raw <- st_read(boundary_gpkg[["harbour_limits"]], quiet = TRUE) %>%
  st_zm(drop = TRUE, what = "ZM") %>% st_transform(27700) %>% st_make_valid()
hits <- lengths(st_intersects(raw, scamp_aoi_sf)) > 0
scamp_boundaries[["harbour_limits"]] <- st_intersection(raw[hits, ], scamp_aoi_sf)
cat(sprintf("  harbour_limits: %d features\n", nrow(scamp_boundaries[["harbour_limits"]])))

cat(sprintf("\nBoundary clipping complete: %d/%d datasets loaded\n",
            length(scamp_boundaries), length(boundary_gpkg)))


##=====================================================================
## STEP 1 — EXPLORE PUTUHENA GDB
##=====================================================================

cat("\n\n===== STEP 1: PUTUHENA GDB EXPLORATION =====\n")

putuhena_gdb <- "data/raw/Putahena_all_GDB/CartOcean_DS_UK_v1_250801.gdb"
stopifnot("Putuhena GDB not found" = file.exists(putuhena_gdb))

put_layers <- st_layers(putuhena_gdb)
cat(sprintf("  GDB path: %s\n", putuhena_gdb))
cat(sprintf("  Number of layers: %d\n", length(put_layers$name)))
cat("\n  Layer summary:\n")
print(data.frame(
  layer     = put_layers$name,
  features  = put_layers$features,
  geomtype  = put_layers$geomtype
))


##=====================================================================
## STEP 1b — PUTUHENA DATASET COMPARISON & EXTRACTION
##=====================================================================

cat("\n\n===== PUTUHENA: OVERLAP ANALYSIS =====\n")

# Putuhena fields grouped by theme
# A-prefix = activity/pressure, E-prefix = ecological

# OVERLAPPING with other datasets in this script:
#   Fishing effort   → fish12 (<12m shp), GFW API (>15m AIS)
#     Putuhena: A65 demersal_fishing_effort, A66 demersal_fishing_time,
#              A67 total_fishing_effort, A68 total_fishing_time,
#              A95-A100 GFW trawler/non-trawler effort
#   Vessel density   → EMODnet vessel density GeoTIFFs (all/cargo/service/tanker)
#     Putuhena: A01 vessel_density_kde, A04 vessel_perday

# UNIQUE to Putuhena (not covered by other datasets):
#   Underwater noise: A05 MoD, A06 seismic, A07 explosives, A08 multibeam,
#                     A09 sub-bottom profilers, A10 acoustic deterrent, A11 piling
#   Infrastructure:   A61 active cable sites, A62 active cables,
#                     A33 wind turbines (kde), A36 wind turbines (exact)
#   Satellite obs:    A77-A79 infrastructure (oil/wind/unknown),
#                     A80-A83 vessel (AIS fishing/non-fishing, dark fishing/non-fishing)
#   Ecology:          E17-E86+ biodiversity interpolations (Shannon, richness, etc.)

cat("  Overlapping themes: fishing effort, vessel density\n")
cat("  Unique themes: underwater noise, infrastructure, dark vessels, ecology\n")

# --- Read Putuhena for SCAMP AOI — first check years and spatial extent ---

cat("\n  Reading x, y, year columns to check coverage...\n")
put_meta <- st_read(putuhena_gdb,
  query = "SELECT x, y, year FROM ED50_10km2_merged_cube_table_v_010825 LIMIT 100",
  quiet = TRUE)
cat(sprintf("  CRS: %s\n", st_crs(put_meta)$input))
cat(sprintf("  Sample x range: %.1f – %.1f\n", min(put_meta$x), max(put_meta$x)))
cat(sprintf("  Sample y range: %.1f – %.1f\n", min(put_meta$y), max(put_meta$y)))

# Get full year range
put_years <- st_read(putuhena_gdb,
  query = "SELECT DISTINCT year FROM ED50_10km2_merged_cube_table_v_010825 ORDER BY year",
  quiet = TRUE)
cat(sprintf("  Available years: %s\n", paste(put_years$year, collapse = ", ")))
latest_year <- max(put_years$year)
avg_years   <- tail(sort(put_years$year), 3)
cat(sprintf("  Latest year: %d\n", latest_year))
cat(sprintf("  3-year average window: %s\n", paste(avg_years, collapse = ", ")))

# --- Read full Putuhena clipped to SCAMP AOI ---
# Transform AOI bbox to the Putuhena CRS for spatial filtering
put_crs <- st_crs(put_meta)
aoi_put <- st_transform(scamp_aoi_sf, put_crs)
put_bbox <- st_bbox(aoi_put)

# Read all columns for latest 3 years, spatially filtered
cat(sprintf("\n  Reading Putuhena data for years %s within SCAMP bbox...\n",
            paste(avg_years, collapse = ", ")))
put_query <- sprintf(
  "SELECT * FROM ED50_10km2_merged_cube_table_v_010825 WHERE year IN (%s)",
  paste(avg_years, collapse = ", "))
put_raw <- st_read(putuhena_gdb, query = put_query, quiet = TRUE)
cat(sprintf("  Read %d features\n", nrow(put_raw)))

# Spatial crop to SCAMP AOI bbox
put_raw <- st_transform(put_raw, 27700) %>% st_make_valid()
put_crop <- st_crop(put_raw, aoi_bbox)
cat(sprintf("  After bbox crop: %d features\n", nrow(put_crop)))
rm(put_raw)

cat(sprintf("  Years in cropped data: %s\n", paste(sort(unique(put_crop$year)), collapse = ", ")))

# --- COMPARISON: Putuhena fishing vs other fishing datasets ---

cat("\n----- Comparison: Putuhena fishing vs other datasets -----\n")

# Putuhena fishing fields for latest year
put_latest <- put_crop[put_crop$year == latest_year, ]
cat(sprintf("  Putuhena cells in AOI (latest year): %d\n", nrow(put_latest)))

# Total fishing effort from Putuhena (log10 kwh)
put_fish_vals <- put_latest$F_L_A67__Layer_of_total_fishing_effort_log10_kwh
cat(sprintf("  Putuhena total fishing effort (log10 kwh):\n"))
cat(sprintf("    Non-NA: %d / %d\n", sum(!is.na(put_fish_vals)), length(put_fish_vals)))
cat(sprintf("    Range: %.3f – %.3f\n",
            min(put_fish_vals, na.rm = TRUE), max(put_fish_vals, na.rm = TRUE)))
cat(sprintf("    Mean: %.3f | Median: %.3f\n",
            mean(put_fish_vals, na.rm = TRUE), median(put_fish_vals, na.rm = TRUE)))

# Putuhena demersal fishing effort
put_dem_vals <- put_latest$F_L_A65__Layer_of_demersal_fishing_effort_log10_kwh
cat(sprintf("  Putuhena demersal fishing effort (log10 kwh):\n"))
cat(sprintf("    Non-NA: %d / %d\n", sum(!is.na(put_dem_vals)), length(put_dem_vals)))
if (sum(!is.na(put_dem_vals)) > 0) {
  cat(sprintf("    Range: %.3f – %.3f\n",
              min(put_dem_vals, na.rm = TRUE), max(put_dem_vals, na.rm = TRUE)))
}

# Compare spatial extent: Putuhena grid vs fish12 vs GFW
cat(sprintf("\n  Spatial coverage comparison:\n"))
cat(sprintf("    Putuhena (10km² grid): %d cells in AOI\n", nrow(put_latest)))
cat(sprintf("    Fishing <12m (0.05° C-squares): %d cells in AOI\n", nrow(fish12_crop)))
if (exists("gfw_effort")) {
  gfw_cells <- nrow(distinct(gfw_effort, Lon, Lat))
  cat(sprintf("    GFW AIS (0.01° grid): %d cells in AOI\n", gfw_cells))
}

# --- COMPARISON: Putuhena vessel density vs EMODnet ---

cat("\n----- Comparison: Putuhena vessel density vs EMODnet -----\n")

put_vess_vals <- put_latest$F_L_A01__Layer_of_kernel_density_estimation_of_vessel_perday
cat(sprintf("  Putuhena vessel KDE:\n"))
cat(sprintf("    Non-NA: %d / %d\n", sum(!is.na(put_vess_vals)), length(put_vess_vals)))
if (sum(!is.na(put_vess_vals)) > 0) {
  cat(sprintf("    Range: %.4f – %.4f\n",
              min(put_vess_vals, na.rm = TRUE), max(put_vess_vals, na.rm = TRUE)))
}

cat(sprintf("  EMODnet vessel density (all types, mean 2017–2024):\n"))
cat(sprintf("    Resolution: ~%dm (BNG)\n", round(res(vd_all_bng)[1])))
cat(sprintf("    Range: %.2f – %.2f hours/km²/month\n",
            minmax(vd_all_bng)[1], minmax(vd_all_bng)[2]))

# --- UNIQUE PUTUHENA LAYERS: Extract for SCAMP AOI ---

cat("\n\n----- Putuhena unique layers: underwater noise -----\n")

# 3-year mean for noise fields (latest year for static infrastructure)
noise_cols <- c(
  "noise_kde_MoD"           = "F_L_A05__Layer_of_kde_of_noises_of_Ministry_of_Defence",
  "noise_kde_seismic"       = "F_L_A06__Layer_of_kde_of_noises_of_Seismic_Survey",
  "noise_kde_explosives"    = "F_L_A07__Layer_of_kde_of_noises_of_Explosives",
  "noise_kde_multibeam"     = "F_L_A08__Layer_of_kde_of_noises_of_Multibeam_Echosounders",
  "noise_kde_subbottom"     = "F_L_A09__Layer_of_kde_of_noises_of_Sub_bottom_Profilers",
  "noise_kde_ADD"           = "F_L_A10__Layer_of_kde_of_noises_of_Acoustic_Deterrent_Device",
  "noise_kde_piling"        = "F_L_A11__Layer_of_kde_of_noises_of_Piling"
)

cat("  Noise KDE fields (3-year mean across cells):\n")
put_noise_mean <- put_crop %>%
  st_drop_geometry() %>%
  group_by(x, y) %>%
  summarise(across(all_of(unname(noise_cols)), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

# Report which noise types have data in the AOI
noise_summary <- data.frame(
  layer = names(noise_cols),
  field = unname(noise_cols),
  non_na = sapply(unname(noise_cols), function(f) sum(!is.na(put_noise_mean[[f]]))),
  non_zero = sapply(unname(noise_cols), function(f) sum(put_noise_mean[[f]] > 0, na.rm = TRUE)),
  max_val = sapply(unname(noise_cols), function(f) max(put_noise_mean[[f]], na.rm = TRUE))
)
print(noise_summary[, c("layer", "non_na", "non_zero", "max_val")])


cat("\n----- Putuhena unique layers: infrastructure -----\n")

infra_cols <- c(
  "cables_active"     = "F_L_A62__Layer_of_CuOW_active_cable",
  "cable_sites"       = "F_L_A61__Layer_of_CuOW_active_sites",
  "wind_turbines_kde" = "F_L_A33__Layer_of_kernel_density_estimation_of_observed_wind_tur",
  "wind_turbines_n"   = "F_L_A36__Layer_of_exact_amount_observed_wind_turbines",
  "infra_oil"         = "F_L_A78__Layer_of_kde_of_infra_satellite_obs_of_tot_oil",
  "infra_wind"        = "F_L_A77__Layer_of_kde_of_infra_satellite_obs_of_tot_wind",
  "infra_unknown"     = "F_L_A79__Layer_of_kde_of_infra_satellite_obs_of_unknown"
)

# Use latest year for infrastructure
infra_summary <- data.frame(
  layer = names(infra_cols),
  non_na = sapply(unname(infra_cols), function(f) sum(!is.na(put_latest[[f]]))),
  non_zero = sapply(unname(infra_cols), function(f) sum(put_latest[[f]] > 0, na.rm = TRUE)),
  max_val = sapply(unname(infra_cols), function(f) max(put_latest[[f]], na.rm = TRUE))
)
cat("  Infrastructure fields (latest year):\n")
print(infra_summary[, c("layer", "non_na", "non_zero", "max_val")])


cat("\n----- Putuhena unique layers: dark vessels -----\n")

dark_cols <- c(
  "dark_fishing_kde"     = "F_L_A82__Layer_of_kde_of_vessel_satellite_obs_of_dark_fishing",
  "dark_nonfishing_kde"  = "F_L_A83__Layer_of_kde_of_vessel_satellite_obs_of_dark_non_fishin",
  "AIS_fishing_kde"      = "F_L_A80__Layer_of_kde_of_vessel_satellite_obs_of_AIS_fishing",
  "AIS_nonfishing_kde"   = "F_L_A81__Layer_of_kde_of_vessel_satellite_obs_of_AIS_nonfishing"
)

# 3-year mean for vessel activity
dark_mean <- put_crop %>%
  st_drop_geometry() %>%
  group_by(x, y) %>%
  summarise(across(all_of(unname(dark_cols)), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

dark_summary <- data.frame(
  layer = names(dark_cols),
  non_na = sapply(unname(dark_cols), function(f) sum(!is.na(dark_mean[[f]]))),
  non_zero = sapply(unname(dark_cols), function(f) sum(dark_mean[[f]] > 0, na.rm = TRUE)),
  max_val = sapply(unname(dark_cols), function(f) max(dark_mean[[f]], na.rm = TRUE))
)
cat("  Vessel satellite obs (3-year mean):\n")
print(dark_summary[, c("layer", "non_na", "non_zero", "max_val")])


cat("\n===== PUTUHENA COMPARISON COMPLETE =====\n")


##=====================================================================
## STEP 2 — READ INDIVIDUAL THREAT DATASETS
##=====================================================================

cat("\n\n===== STEP 2: INDIVIDUAL THREAT DATASETS =====\n")

# ------------------------------------------------------------------
# 2a. Fishing <12m
# ------------------------------------------------------------------

cat("\n----- 2a. Fishing <12m -----\n")

fish12_path <- "data/raw/threats/local/fishing_under12_data/fishing_under12_dataPolygon.shp"
stopifnot("Fishing <12m shapefile not found" = file.exists(fish12_path))

fish12 <- st_read(fish12_path, quiet = TRUE) %>% st_transform(27700) %>% st_make_valid()
cat(sprintf("  Class: %s\n", paste(class(fish12), collapse = ", ")))
cat(sprintf("  CRS: %s\n", st_crs(fish12)$input))
cat(sprintf("  Features: %d\n", nrow(fish12)))
cat(sprintf("  Columns: %s\n", paste(names(fish12), collapse = ", ")))

fish12_crop <- st_crop(fish12, aoi_bbox)
cat(sprintf("  After bbox crop: %d features\n", nrow(fish12_crop)))

# Investigate data quality — how much is redacted (-999)?
cat(sprintf("\n  Gear categories: %s\n",
            paste(unique(fish12_crop$category), collapse = ", ")))
cat(sprintf("  Records by category:\n"))
print(table(fish12_crop$category))

n_total    <- nrow(fish12_crop)
n_redacted <- sum(fish12_crop$num_vess == -999, na.rm = TRUE)
n_valid    <- n_total - n_redacted
cat(sprintf("\n  Total cells in AOI: %d\n", n_total))
cat(sprintf("  Redacted (<5 vessels): %d (%.0f%%)\n", n_redacted, 100 * n_redacted / n_total))
cat(sprintf("  With valid data: %d (%.0f%%)\n", n_valid, 100 * n_valid / n_total))

# For valid records, summarise values
fish12_valid <- fish12_crop[fish12_crop$num_vess != -999, ]
if (nrow(fish12_valid) > 0) {
  cat(sprintf("\n  Valid records — num_vess range: %d – %d\n",
              min(fish12_valid$num_vess), max(fish12_valid$num_vess)))
  cat(sprintf("  Valid records — avgdays range: %.1f – %.1f\n",
              min(fish12_valid$avgdays), max(fish12_valid$avgdays)))
  cat(sprintf("  Valid records — avgwtkg range: %.1f – %.1f\n",
              min(fish12_valid$avgwtkg), max(fish12_valid$avgwtkg)))
  cat("  Valid records by category:\n")
  print(table(fish12_valid$category))
}

# Even redacted cells confirm fishing presence, so presence/absence by category is useful
# Plot by gear category — presence gives spatial footprint even when values are redacted
p_fish12 <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue", linewidth = 0.3) +
  geom_sf(data = fish12_crop, aes(fill = category), colour = NA, alpha = 0.6) +
  scale_fill_brewer(name = "Gear type", palette = "Set2") +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  labs(title = "Fishing Activity <12m vessels (2018–2022)",
       subtitle = "C-square presence by gear type (0.05° grid)") +
  map_theme + coord_sf(crs = 27700)
ggsave("figures/threats/fishing_under12.png", p_fish12,
       width = 12, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/threats/fishing_under12.png\n")

# Where valid data exists, plot fishing intensity (avg days fished)
if (nrow(fish12_valid) > 0) {
  p_fish12_intensity <- ggplot() +
    geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
    geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue", linewidth = 0.3) +
    geom_sf(data = fish12_valid, aes(fill = avgdays), colour = NA, alpha = 0.7) +
    scale_fill_viridis_c(name = "Avg days\nfished", option = "inferno") +
    geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
    labs(title = "Fishing Intensity <12m vessels (2018–2022)",
         subtitle = "Average days fished per vessel (non-redacted cells only)") +
    map_theme + coord_sf(crs = 27700)
  ggsave("figures/threats/fishing_under12_intensity.png", p_fish12_intensity,
         width = 12, height = 10, units = "in", dpi = 300)
  cat("  Saved: figures/threats/fishing_under12_intensity.png\n")
}


# ------------------------------------------------------------------
# 2b. Fishing effort — Global Fishing Watch (gfwr)
# ------------------------------------------------------------------

cat("\n----- 2b. Fishing effort (Global Fishing Watch) -----\n")
require(gfwr)
?gfw_auth

gfw_token <- gfwr::gfw_auth()
cat(sprintf("  GFW token loaded: %s...\n", substr(gfw_token, 1, 8)))

# Bbox from SCAMP AOI in WGS84 for the API query
aoi_wgs84_bbox <- st_bbox(st_transform(scamp_aoi_sf, 4326))

# Fetch apparent fishing effort (hours) at 0.01° resolution, one year at a time
cat("  Querying GFW API year by year...\n")

gfw_region_sf <- st_transform(scamp_aoi_sf, 4326)

cat("  2020...")
gfw_2020 <- gfw_ais_fishing_hours(
  spatial_resolution = "HIGH", temporal_resolution = "YEARLY",
  group_by = "GEARTYPE", start_date = "2020-01-01", end_date = "2021-01-01",
  region = gfw_region_sf, region_source = "USER_SHAPEFILE",
  key = Sys.getenv("GFW_TOKEN"))
cat(" done\n  2021...")

gfw_2021 <- gfw_ais_fishing_hours(
  spatial_resolution = "HIGH", temporal_resolution = "YEARLY",
  group_by = "GEARTYPE", start_date = "2021-01-01", end_date = "2022-01-01",
  region = gfw_region_sf, region_source = "USER_SHAPEFILE", key = gfw_token)
cat(" done\n  2022...")

gfw_2022 <- gfw_ais_fishing_hours(
  spatial_resolution = "HIGH", temporal_resolution = "YEARLY",
  group_by = "GEARTYPE", start_date = "2022-01-01", end_date = "2023-01-01",
  region = gfw_region_sf, region_source = "USER_SHAPEFILE", key = gfw_token)
cat(" done\n  2023...")

gfw_2023 <- gfw_ais_fishing_hours(
  spatial_resolution = "HIGH", temporal_resolution = "YEARLY",
  group_by = "GEARTYPE", start_date = "2023-01-01", end_date = "2024-01-01",
  region = gfw_region_sf, region_source = "USER_SHAPEFILE", key = gfw_token)
cat(" done\n  2024...")

gfw_2024 <- gfw_ais_fishing_hours(
  spatial_resolution = "HIGH", temporal_resolution = "YEARLY",
  group_by = "GEARTYPE", start_date = "2024-01-01", end_date = "2025-01-01",
  region = gfw_region_sf, region_source = "USER_SHAPEFILE", key = gfw_token)
cat(" done\n")

gfw_effort <- bind_rows(gfw_2020, gfw_2021, gfw_2022, gfw_2023, gfw_2024)
rm(gfw_2020, gfw_2021, gfw_2022, gfw_2023, gfw_2024)

cat(sprintf("  Rows returned: %d\n", nrow(gfw_effort)))
cat(sprintf("  Columns: %s\n", paste(names(gfw_effort), collapse = ", ")))
cat(sprintf("  Gear types: %s\n", paste(unique(gfw_effort$geartype), collapse = ", ")))

# Convert to sf points and project to BNG
gfw_sf <- st_as_sf(gfw_effort, coords = c("Lon", "Lat"), crs = 4326) %>%
  st_transform(27700)
cat(sprintf("  GFW points in AOI: %d\n", nrow(gfw_sf)))

summary(gfw_sf)

# Aggregate total fishing hours per cell across all years and gear types
gfw_total <- gfw_effort %>%
  group_by(Lon, Lat) %>%
  summarise(total_hours = sum(`Apparent Fishing Hours`, na.rm = TRUE), .groups = "drop")
gfw_total_sf <- st_as_sf(gfw_total, coords = c("Lon", "Lat"), crs = 4326) %>%
  st_transform(27700)

p_gfw <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue", linewidth = 0.3) +
  geom_sf(data = gfw_total_sf, aes(colour = total_hours), size = 1.5, alpha = 0.8) +
  scale_colour_viridis_c(name = "Total fishing\nhours (2020–24)", option = "inferno") +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  labs(title = "Apparent Fishing Effort (Global Fishing Watch, 2020–2024)") +
  map_theme + coord_sf(crs = 27700)
ggsave("figures/threats/gfw_fishing_effort.png", p_gfw,
       width = 12, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/threats/gfw_fishing_effort.png\n")

# Summary by gear type
gfw_gear_summary <- gfw_effort %>%
  group_by(geartype) %>%
  summarise(
    total_hours = sum(`Apparent Fishing Hours`, na.rm = TRUE),
    n_cells     = n(),
    .groups     = "drop"
  ) %>%
  arrange(desc(total_hours))
cat("\n  Fishing effort by gear type:\n")
print(gfw_gear_summary)


# ------------------------------------------------------------------
# 2c. Waste discharge points
# ------------------------------------------------------------------

cat("\n----- 2c. Waste discharge points -----\n")

waste_path <- "data/raw/threats/local/EMODnet_HA_WasteDisposal_UWWTD_DischargePoints_20250718/EMODnet_HA_WasteDisposal_UWWTP_DischargePoints_Locations_20250718.shp"
stopifnot("Waste discharge shapefile not found" = file.exists(waste_path))

waste_pts <- st_read(waste_path, quiet = TRUE) %>% st_transform(27700) %>% st_make_valid()
cat(sprintf("  Class: %s\n", paste(class(waste_pts), collapse = ", ")))
cat(sprintf("  CRS: %s\n", st_crs(waste_pts)$input))
cat(sprintf("  Features: %d\n", nrow(waste_pts)))
cat(sprintf("  Columns: %s\n", paste(names(waste_pts), collapse = ", ")))

waste_crop <- st_crop(waste_pts, aoi_bbox)
cat(sprintf("  After bbox crop: %d features\n", nrow(waste_crop)))

p_waste <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue", linewidth = 0.3) +
  geom_sf(data = waste_crop, colour = "brown", size = 2, alpha = 0.7) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  labs(title = "Waste Discharge Points (UWWTD)") +
  map_theme + coord_sf(crs = 27700)
ggsave("figures/threats/waste_discharge.png", p_waste,
       width = 12, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/threats/waste_discharge.png\n")


# ------------------------------------------------------------------
# 2d. Dredge spoil dumping (polygons + points)
# ------------------------------------------------------------------

cat("\n----- 2d. Dredge spoil dumping -----\n")

spoil_pg_path <- "data/raw/threats/local/EMODnet_HA_DredgeSpoilDumping_20251105/EMODnet_HA_DredgeSpoilDumping_pg_Locations_20251105.shp"
spoil_pt_path <- "data/raw/threats/local/EMODnet_HA_DredgeSpoilDumping_20251105/EMODnet_HA_DredgeSpoilDumping_pt_Locations_20251105.shp"
stopifnot("Dredge spoil polygon shapefile not found" = file.exists(spoil_pg_path))
stopifnot("Dredge spoil point shapefile not found" = file.exists(spoil_pt_path))

spoil_pg <- st_read(spoil_pg_path, quiet = TRUE) %>% st_transform(27700) %>% st_make_valid()
cat(sprintf("  Polygons — Class: %s | CRS: %s | Features: %d\n",
            paste(class(spoil_pg), collapse = ", "), st_crs(spoil_pg)$input, nrow(spoil_pg)))
cat(sprintf("  Polygon columns: %s\n", paste(names(spoil_pg), collapse = ", ")))

spoil_pt <- st_read(spoil_pt_path, quiet = TRUE) %>% st_transform(27700) %>% st_make_valid()
cat(sprintf("  Points — Class: %s | CRS: %s | Features: %d\n",
            paste(class(spoil_pt), collapse = ", "), st_crs(spoil_pt)$input, nrow(spoil_pt)))
cat(sprintf("  Point columns: %s\n", paste(names(spoil_pt), collapse = ", ")))

spoil_pg_crop <- st_crop(spoil_pg, aoi_bbox)
spoil_pt_crop <- st_crop(spoil_pt, aoi_bbox)
cat(sprintf("  After bbox crop: %d polygons, %d points\n", nrow(spoil_pg_crop), nrow(spoil_pt_crop)))

p_spoil <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue", linewidth = 0.3) +
  {if (nrow(spoil_pg_crop) > 0) geom_sf(data = spoil_pg_crop, fill = "darkred", colour = NA, alpha = 0.5)} +
  {if (nrow(spoil_pt_crop) > 0) geom_sf(data = spoil_pt_crop, colour = "darkred", size = 2, alpha = 0.7)} +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  labs(title = "Dredge Spoil Dumping Sites") +
  map_theme + coord_sf(crs = 27700)
ggsave("figures/threats/dredge_spoil.png", p_spoil,
       width = 12, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/threats/dredge_spoil.png\n")


# ------------------------------------------------------------------
# 2e. Dredging
# ------------------------------------------------------------------

cat("\n----- 2e. Dredging -----\n")

dredge_path <- "data/raw/threats/local/EMODnet_HA_Dredging_20251103/EMODnet_HA_Dredging_Locations_20251103.shp"
stopifnot("Dredging shapefile not found" = file.exists(dredge_path))

dredge <- st_read(dredge_path, quiet = TRUE) %>% st_transform(27700) %>% st_make_valid()
cat(sprintf("  Class: %s\n", paste(class(dredge), collapse = ", ")))
cat(sprintf("  CRS: %s\n", st_crs(dredge)$input))
cat(sprintf("  Features: %d\n", nrow(dredge)))
cat(sprintf("  Geometry type: %s\n", unique(st_geometry_type(dredge))))
cat(sprintf("  Columns: %s\n", paste(names(dredge), collapse = ", ")))

dredge_crop <- st_crop(dredge, aoi_bbox)
cat(sprintf("  After bbox crop: %d features\n", nrow(dredge_crop)))

p_dredge <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue", linewidth = 0.3) +
  geom_sf(data = dredge_crop, fill = "sienna", colour = "sienna4", alpha = 0.5) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  labs(title = "Dredging Locations") +
  map_theme + coord_sf(crs = 27700)
ggsave("figures/threats/dredging.png", p_dredge,
       width = 12, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/threats/dredging.png\n")


# ------------------------------------------------------------------
# 2f. Vessel density — all types (monthly GeoTIFFs → mean)
# ------------------------------------------------------------------

cat("\n----- 2f. Vessel density — all types -----\n")

vd_all_dir <- "data/raw/threats/local/EMODnet_HA_Vessel_Density_all_2017-24_monthly average/EMODnet_HA_Vessel_Density_all"
vd_all_tifs <- list.files(vd_all_dir, pattern = "\\.tif$", full.names = TRUE)
cat(sprintf("  GeoTIFFs found: %d\n", length(vd_all_tifs)))

vd_all_stack <- rast(vd_all_tifs)
cat(sprintf("  Stack: %d layers | %d x %d | res: %.4f\n",
            nlyr(vd_all_stack), ncol(vd_all_stack), nrow(vd_all_stack), res(vd_all_stack)[1]))
cat(sprintf("  CRS: %s\n", crs(vd_all_stack, describe = TRUE)$code))

aoi_vd <- st_transform(scamp_aoi_sf, crs(vd_all_stack))
vd_all_crop <- crop(vd_all_stack, vect(aoi_vd))
vd_all_mean <- mean(vd_all_crop, na.rm = TRUE)
cat(sprintf("  Mean density range: %.2f to %.2f\n",
            minmax(vd_all_mean)[1], minmax(vd_all_mean)[2]))

vd_all_bng <- project(vd_all_mean, "EPSG:27700", method = "bilinear")

p_vd_all <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  tidyterra::geom_spatraster(data = vd_all_bng) +
  scale_fill_viridis_c(name = "Hours/km²/month", na.value = "transparent", option = "inferno") +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  labs(title = "Vessel Density — All Types (2017–2024 mean)") +
  map_theme + coord_sf(crs = 27700)
ggsave("figures/threats/vessel_density_all.png", p_vd_all,
       width = 12, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/threats/vessel_density_all.png\n")


# ------------------------------------------------------------------
# 2g. Vessel density — cargo
# ------------------------------------------------------------------

cat("\n----- 2g. Vessel density — cargo -----\n")

vd_cargo_dir <- "data/raw/threats/local/EMODnet_HA_Vessel_Density_cargo_hours per sq km per month/EMODnet_HA_Vessel_Density_09"
vd_cargo_tifs <- list.files(vd_cargo_dir, pattern = "\\.tif$", full.names = TRUE)
cat(sprintf("  GeoTIFFs found: %d\n", length(vd_cargo_tifs)))

vd_cargo_stack <- rast(vd_cargo_tifs)
vd_cargo_crop <- crop(vd_cargo_stack, vect(aoi_vd))
vd_cargo_mean <- mean(vd_cargo_crop, na.rm = TRUE)
cat(sprintf("  Mean density range: %.2f to %.2f\n",
            minmax(vd_cargo_mean)[1], minmax(vd_cargo_mean)[2]))

vd_cargo_bng <- project(vd_cargo_mean, "EPSG:27700", method = "bilinear")

p_vd_cargo <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  tidyterra::geom_spatraster(data = vd_cargo_bng) +
  scale_fill_viridis_c(name = "Hours/km²/month", na.value = "transparent", option = "inferno") +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  labs(title = "Vessel Density — Cargo (2017–2024 mean)") +
  map_theme + coord_sf(crs = 27700)
ggsave("figures/threats/vessel_density_cargo.png", p_vd_cargo,
       width = 12, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/threats/vessel_density_cargo.png\n")


# ------------------------------------------------------------------
# 2h. Vessel density — service
# ------------------------------------------------------------------

cat("\n----- 2h. Vessel density — service -----\n")

vd_svc_dir <- "data/raw/threats/local/EMODnet_HA_Vessel_Density_service_hours per sq km per month/EMODnet_HA_Vessel_Density_02"
vd_svc_tifs <- list.files(vd_svc_dir, pattern = "\\.tif$", full.names = TRUE)
cat(sprintf("  GeoTIFFs found: %d\n", length(vd_svc_tifs)))

vd_svc_stack <- rast(vd_svc_tifs)
vd_svc_crop <- crop(vd_svc_stack, vect(aoi_vd))
vd_svc_mean <- mean(vd_svc_crop, na.rm = TRUE)
cat(sprintf("  Mean density range: %.2f to %.2f\n",
            minmax(vd_svc_mean)[1], minmax(vd_svc_mean)[2]))

vd_svc_bng <- project(vd_svc_mean, "EPSG:27700", method = "bilinear")

p_vd_svc <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  tidyterra::geom_spatraster(data = vd_svc_bng) +
  scale_fill_viridis_c(name = "Hours/km²/month", na.value = "transparent", option = "inferno") +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  labs(title = "Vessel Density — Service (2017–2024 mean)") +
  map_theme + coord_sf(crs = 27700)
ggsave("figures/threats/vessel_density_service.png", p_vd_svc,
       width = 12, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/threats/vessel_density_service.png\n")


# ------------------------------------------------------------------
# 2i. Vessel density — tanker
# ------------------------------------------------------------------

cat("\n----- 2i. Vessel density — tanker -----\n")

vd_tank_dir <- "data/raw/threats/local/EMODnet_HA_Vessel_Density_tanker_hours per sq km per month/EMODnet_HA_Vessel_Density_10"
vd_tank_tifs <- list.files(vd_tank_dir, pattern = "\\.tif$", full.names = TRUE)
cat(sprintf("  GeoTIFFs found: %d\n", length(vd_tank_tifs)))

vd_tank_stack <- rast(vd_tank_tifs)
vd_tank_crop <- crop(vd_tank_stack, vect(aoi_vd))
vd_tank_mean <- mean(vd_tank_crop, na.rm = TRUE)
cat(sprintf("  Mean density range: %.2f to %.2f\n",
            minmax(vd_tank_mean)[1], minmax(vd_tank_mean)[2]))

vd_tank_bng <- project(vd_tank_mean, "EPSG:27700", method = "bilinear")

p_vd_tank <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  tidyterra::geom_spatraster(data = vd_tank_bng) +
  scale_fill_viridis_c(name = "Hours/km²/month", na.value = "transparent", option = "inferno") +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  labs(title = "Vessel Density — Tanker (2017–2024 mean)") +
  map_theme + coord_sf(crs = 27700)
ggsave("figures/threats/vessel_density_tanker.png", p_vd_tank,
       width = 12, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/threats/vessel_density_tanker.png\n")


# ------------------------------------------------------------------
# 2j. UK Ports
# ------------------------------------------------------------------

cat("\n----- 2j. UK Ports -----\n")

ports_path <- "data/raw/threats/local/UK_Ports_coastal/UK_Ports_coastal.shp"
stopifnot("UK Ports shapefile not found" = file.exists(ports_path))

ports <- st_read(ports_path, quiet = TRUE) %>% st_transform(27700) %>% st_make_valid()
cat(sprintf("  Class: %s\n", paste(class(ports), collapse = ", ")))
cat(sprintf("  CRS: %s\n", st_crs(ports)$input))
cat(sprintf("  Features: %d\n", nrow(ports)))
cat(sprintf("  Columns: %s\n", paste(names(ports), collapse = ", ")))

ports_crop <- st_crop(ports, aoi_bbox)
cat(sprintf("  After bbox crop: %d features\n", nrow(ports_crop)))

p_ports <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue", linewidth = 0.3) +
  geom_sf(data = ports_crop, colour = "navy", size = 3, shape = 17, alpha = 0.8) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  labs(title = "UK Coastal Ports") +
  map_theme + coord_sf(crs = 27700)
ggsave("figures/threats/uk_ports.png", p_ports,
       width = 12, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/threats/uk_ports.png\n")


# ------------------------------------------------------------------
# 2k. Non-native macroalgae models
# ------------------------------------------------------------------

cat("\n----- 2k. Non-native macroalgae models -----\n")

macro_path <- "data/raw/threats/local/non-native_macroalgae_kuenm_models_review.nc"
stopifnot("Non-native macroalgae NetCDF not found" = file.exists(macro_path))

macro_nc <- rast(macro_path)
cat(sprintf("  Class: %s\n", paste(class(macro_nc), collapse = ", ")))
cat(sprintf("  Layers: %d | %d x %d\n", nlyr(macro_nc), ncol(macro_nc), nrow(macro_nc)))
cat(sprintf("  CRS: %s\n", crs(macro_nc, describe = TRUE)$code))
cat(sprintf("  Resolution: %.4f x %.4f\n", res(macro_nc)[1], res(macro_nc)[2]))
cat(sprintf("  Layer names: %s\n", paste(names(macro_nc), collapse = ", ")))

aoi_macro <- st_transform(scamp_aoi_sf, crs(macro_nc))
macro_crop <- crop(macro_nc, vect(aoi_macro))
cat(sprintf("  After bbox crop: %d x %d\n", ncol(macro_crop), nrow(macro_crop)))

macro_bng <- project(macro_crop, "EPSG:27700", method = "bilinear")

p_macro <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  tidyterra::geom_spatraster(data = macro_bng[[1]]) +
  scale_fill_viridis_c(name = "Suitability", na.value = "transparent") +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  labs(title = sprintf("Non-native Macroalgae — %s", names(macro_bng)[1])) +
  map_theme + coord_sf(crs = 27700)
ggsave("figures/threats/macroalgae_nonnative.png", p_macro,
       width = 12, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/threats/macroalgae_nonnative.png\n")


# ------------------------------------------------------------------
# 2l. WW Wrecks
# ------------------------------------------------------------------

cat("\n----- 2l. WW Wrecks -----\n")

wrecks_gdb <- "data/raw/threats/local/EMODnet_HA_Heritage_WW_Wrecks_20241226/EMODnet_HA_Heritage_WW_Wrecks_20241226.gdb"
stopifnot("WW Wrecks GDB not found" = file.exists(wrecks_gdb))

wrecks_layers <- st_layers(wrecks_gdb)
cat(sprintf("  GDB layers: %s\n", paste(wrecks_layers$name, collapse = ", ")))
cat(sprintf("  Features: %s\n", paste(wrecks_layers$features, collapse = ", ")))
cat(sprintf("  Geom types: %s\n", paste(wrecks_layers$geomtype, collapse = ", ")))

wrecks <- st_read(wrecks_gdb, layer = wrecks_layers$name[1], quiet = TRUE) %>%
  st_transform(27700) %>% st_make_valid()
cat(sprintf("  Class: %s\n", paste(class(wrecks), collapse = ", ")))
cat(sprintf("  CRS: %s\n", st_crs(wrecks)$input))
cat(sprintf("  Features: %d\n", nrow(wrecks)))
cat(sprintf("  Columns: %s\n", paste(names(wrecks), collapse = ", ")))

wrecks_crop <- st_crop(wrecks, aoi_bbox)
cat(sprintf("  After bbox crop: %d features\n", nrow(wrecks_crop)))

p_wrecks <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue", linewidth = 0.3) +
  geom_sf(data = wrecks_crop, colour = "darkgrey", size = 1.5, shape = 4, alpha = 0.7) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  labs(title = "WW Wrecks") +
  map_theme + coord_sf(crs = 27700)
ggsave("figures/threats/ww_wrecks.png", p_wrecks,
       width = 12, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/threats/ww_wrecks.png\n")


# ------------------------------------------------------------------
# 2m. Dumped munitions (polygons + points)
# ------------------------------------------------------------------

cat("\n----- 2m. Dumped munitions -----\n")

mun_pg_path <- "data/raw/threats/local/EMODnet_HA_WasteDisposal_DumpedMunitions_20260609/EMODnet_HA_WasteDisposal_DumpedMunitions_pg_20260609.shp"
mun_pt_path <- "data/raw/threats/local/EMODnet_HA_WasteDisposal_DumpedMunitions_20260609/EMODnet_HA_WasteDisposal_DumpedMunitions_pt_20260609.shp"
stopifnot("Munitions polygon shapefile not found" = file.exists(mun_pg_path))
stopifnot("Munitions point shapefile not found" = file.exists(mun_pt_path))

mun_pg <- st_read(mun_pg_path, quiet = TRUE) %>% st_transform(27700) %>% st_make_valid()
cat(sprintf("  Polygons — Class: %s | CRS: %s | Features: %d\n",
            paste(class(mun_pg), collapse = ", "), st_crs(mun_pg)$input, nrow(mun_pg)))
cat(sprintf("  Polygon columns: %s\n", paste(names(mun_pg), collapse = ", ")))

mun_pt <- st_read(mun_pt_path, quiet = TRUE) %>% st_transform(27700) %>% st_make_valid()
cat(sprintf("  Points — Class: %s | CRS: %s | Features: %d\n",
            paste(class(mun_pt), collapse = ", "), st_crs(mun_pt)$input, nrow(mun_pt)))
cat(sprintf("  Point columns: %s\n", paste(names(mun_pt), collapse = ", ")))

mun_pg_crop <- st_crop(mun_pg, aoi_bbox)
mun_pt_crop <- st_crop(mun_pt, aoi_bbox)
cat(sprintf("  After bbox crop: %d polygons, %d points\n", nrow(mun_pg_crop), nrow(mun_pt_crop)))

p_munitions <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue", linewidth = 0.3) +
  {if (nrow(mun_pg_crop) > 0) geom_sf(data = mun_pg_crop, fill = "black", colour = NA, alpha = 0.4)} +
  {if (nrow(mun_pt_crop) > 0) geom_sf(data = mun_pt_crop, colour = "black", size = 2, shape = 8, alpha = 0.7)} +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  labs(title = "Dumped Munitions") +
  map_theme + coord_sf(crs = 27700)
ggsave("figures/threats/dumped_munitions.png", p_munitions,
       width = 12, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/threats/dumped_munitions.png\n")


##=====================================================================
## CATCHMENTS AND RUNOFF
##=====================================================================

cat("\n\n===== CATCHMENTS AND RUNOFF =====\n")

# ------------------------------------------------------------------
# OpenRivers
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# CEFAS catchments
# ------------------------------------------------------------------

cat("\n----- CEFAS catchments -----\n")

cefas_path <- "data/raw/physical/CEFAS catchments/Recordset_12840/MDRExternalStorage.Recordset_12840.shp"
stopifnot("CEFAS catchments shapefile not found" = file.exists(cefas_path))

cefas <- st_read(cefas_path, quiet = TRUE) %>% st_set_crs(27700) %>% st_make_valid()
cat(sprintf("  Class: %s\n", paste(class(cefas), collapse = ", ")))
cat(sprintf("  CRS: %s\n", st_crs(cefas)$input))
cat(sprintf("  Features: %d\n", nrow(cefas)))
cat(sprintf("  Columns: %s\n", paste(names(cefas), collapse = ", ")))

cefas_hits <- lengths(st_intersects(cefas, scamp_aoi_sf)) > 0
cefas_crop <- cefas[cefas_hits, ]
cat(sprintf("  Catchments draining into SCAMP AOI: %d features\n", nrow(cefas_crop)))
cat(sprintf("  Columns: %s\n", paste(names(cefas_crop), collapse = ", ")))

p_cefas <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = cefas_crop, fill = "lightblue", colour = "steelblue", alpha = 0.4, linewidth = 0.3) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  labs(title = "CEFAS Catchments") +
  map_theme + coord_sf(crs = 27700)
ggsave("figures/threats/cefas_catchments.png", p_cefas,
       width = 12, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/threats/cefas_catchments.png\n")

names(cefas_crop)
head(st_drop_geometry(cefas_crop), 5)


# ------------------------------------------------------------------
# OpenRivers (masked to catchments)
# ------------------------------------------------------------------

cat("\n----- OpenRivers -----\n")

rivers_path <- "data/raw/physical/OpenRivers/Data/oprvrs_gb.gpkg"
stopifnot("OpenRivers GPKG not found" = file.exists(rivers_path))

rivers_layers <- st_layers(rivers_path)
cat(sprintf("  Layers: %s\n", paste(rivers_layers$name, collapse = ", ")))
cat(sprintf("  Features: %s\n", paste(rivers_layers$features, collapse = ", ")))

rivers <- st_read(rivers_path, layer = "watercourse_link", quiet = TRUE) %>%
  st_transform(27700) %>% st_make_valid()
cat(sprintf("  Class: %s\n", paste(class(rivers), collapse = ", ")))
cat(sprintf("  CRS: %s\n", st_crs(rivers)$input))
cat(sprintf("  Features: %d\n", nrow(rivers)))
cat(sprintf("  Columns: %s\n", paste(names(rivers), collapse = ", ")))

catch_u <- st_union(cefas_crop)
rivers_crop <- st_intersection(rivers, catch_u) %>%
  st_make_valid() %>%
  st_collection_extract("LINESTRING")
cat(sprintf("  Rivers clipped to catchments: %d features\n", nrow(rivers_crop)))

p_rivers <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = rivers_crop, colour = "steelblue", linewidth = 0.3, alpha = 0.6) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  labs(title = "OpenRivers — Watercourses") +
  map_theme + coord_sf(crs = 27700)
ggsave("figures/threats/openrivers.png", p_rivers,
       width = 12, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/threats/openrivers.png\n")


##=====================================================================
## STEP 3 — PRESSURE SURFACE CONSTRUCTION
##=====================================================================

cat("\n\n===== STEP 3: PRESSURE SURFACE =====\n")

# --- 3a. Create 1km² reference raster aligned to BNG ---

cat("\n----- 3a. Reference raster (1km² BNG) -----\n")

smr_bbox <- st_bbox(smr_solway_bng)
ref_rast <- rast(
  xmin = floor(smr_bbox["xmin"]),
  xmax = ceiling(smr_bbox["xmax"]),
  ymin = floor(smr_bbox["ymin"]),
  ymax = ceiling(smr_bbox["ymax"]),
  res  = 1000,
  crs  = "EPSG:27700"
)
cat(sprintf("  Reference raster: %d x %d cells at 1km\n", ncol(ref_rast), nrow(ref_rast)))

# Marine mask from SMR boundary
smr_mask <- rasterize(vect(smr_solway_bng), ref_rast, field = 1)
cat(sprintf("  Marine mask cells: %d\n", sum(!is.na(values(smr_mask)))))


# --- 3b. Resample raster pressures to 1km, mask to marine area ---

cat("\n----- 3b. Resample & mask raster pressures -----\n")

normalise_01 <- function(r) {
  v <- values(r)
  vmin <- min(v, na.rm = TRUE)
  vmax <- max(v, na.rm = TRUE)
  if (vmax == vmin) return(r * 0)
  (r - vmin) / (vmax - vmin)
}

# Store raw (original units) and normalised separately
raw_rasters  <- list()
norm_rasters <- list()

# Metadata: layer name, units, source
layer_meta <- list(
  vessel_all     = list(units = "hours/km²/month", source = "EMODnet"),
  vessel_cargo   = list(units = "hours/km²/month", source = "EMODnet"),
  vessel_service = list(units = "hours/km²/month", source = "EMODnet"),
  vessel_tanker  = list(units = "hours/km²/month", source = "EMODnet"),
  macroalgae_inns = list(units = "suitability", source = "kuenm model"),
  fishing_gfw    = list(units = "fishing hours (2020–24)", source = "GFW AIS"),
  fishing_under12 = list(units = "gear-type count", source = "Scottish Gov <12m"),
  noise_MoD       = list(units = "KDE", source = "Putuhena/CartOcean"),
  noise_seismic   = list(units = "KDE", source = "Putuhena/CartOcean"),
  noise_explosives = list(units = "KDE", source = "Putuhena/CartOcean"),
  noise_multibeam = list(units = "KDE", source = "Putuhena/CartOcean"),
  noise_subbottom = list(units = "KDE", source = "Putuhena/CartOcean"),
  noise_ADD       = list(units = "KDE", source = "Putuhena/CartOcean"),
  noise_piling    = list(units = "KDE", source = "Putuhena/CartOcean"),
  infra_cables   = list(units = "count", source = "Putuhena/CartOcean"),
  infra_wind     = list(units = "KDE", source = "Putuhena/CartOcean"),
  infra_oil      = list(units = "KDE", source = "Putuhena/CartOcean"),
  dark_fishing   = list(units = "KDE", source = "Putuhena/CartOcean"),
  dark_nonfishing = list(units = "KDE", source = "Putuhena/CartOcean")
)

# Vessel density — all types
cat("  Resampling vessel density (all)...\n")
vd_all_1km <- resample(vd_all_bng, ref_rast, method = "bilinear")
vd_all_1km <- mask(vd_all_1km, smr_mask)
raw_rasters[["vessel_all"]]  <- vd_all_1km
norm_rasters[["vessel_all"]] <- normalise_01(vd_all_1km)

# Vessel density — cargo
cat("  Resampling vessel density (cargo)...\n")
vd_cargo_1km <- resample(vd_cargo_bng, ref_rast, method = "bilinear")
vd_cargo_1km <- mask(vd_cargo_1km, smr_mask)
raw_rasters[["vessel_cargo"]]  <- vd_cargo_1km
norm_rasters[["vessel_cargo"]] <- normalise_01(vd_cargo_1km)

# Vessel density — service
cat("  Resampling vessel density (service)...\n")
vd_svc_1km <- resample(vd_svc_bng, ref_rast, method = "bilinear")
vd_svc_1km <- mask(vd_svc_1km, smr_mask)
raw_rasters[["vessel_service"]]  <- vd_svc_1km
norm_rasters[["vessel_service"]] <- normalise_01(vd_svc_1km)

# Vessel density — tanker
cat("  Resampling vessel density (tanker)...\n")
vd_tank_1km <- resample(vd_tank_bng, ref_rast, method = "bilinear")
vd_tank_1km <- mask(vd_tank_1km, smr_mask)
raw_rasters[["vessel_tanker"]]  <- vd_tank_1km
norm_rasters[["vessel_tanker"]] <- normalise_01(vd_tank_1km)

# Non-native macroalgae suitability
cat("  Resampling non-native macroalgae...\n")
macro_1km <- resample(macro_bng[[1]], ref_rast, method = "bilinear")
macro_1km <- mask(macro_1km, smr_mask)
raw_rasters[["macroalgae_inns"]]  <- macro_1km
norm_rasters[["macroalgae_inns"]] <- normalise_01(macro_1km)

# GFW fishing effort — rasterise point data
cat("  Rasterising GFW fishing effort...\n")
gfw_total_sf <- st_as_sf(
  gfw_effort %>%
    group_by(Lon, Lat) %>%
    summarise(total_hours = sum(`Apparent Fishing Hours`, na.rm = TRUE), .groups = "drop"),
  coords = c("Lon", "Lat"), crs = 4326
) %>% st_transform(27700)
gfw_rast <- rasterize(vect(gfw_total_sf), ref_rast, field = "total_hours", fun = "sum")
gfw_rast <- mask(gfw_rast, smr_mask)
raw_rasters[["fishing_gfw"]]  <- gfw_rast
norm_rasters[["fishing_gfw"]] <- normalise_01(gfw_rast)

# Fishing <12m — rasterise as gear-type count
cat("  Rasterising fishing <12m...\n")
fish12_rast <- rasterize(vect(fish12_crop), ref_rast, fun = "count")
fish12_rast <- mask(fish12_rast, smr_mask)
raw_rasters[["fishing_under12"]]  <- fish12_rast
norm_rasters[["fishing_under12"]] <- normalise_01(fish12_rast)

# --- Putuhena layers: rasterise polygons to 1km grid ---

cat("  Rasterising Putuhena noise layers...\n")

put_noise_sf <- put_crop %>%
  group_by(x, y) %>%
  summarise(
    noise_MoD        = mean(F_L_A05__Layer_of_kde_of_noises_of_Ministry_of_Defence, na.rm = TRUE),
    noise_seismic    = mean(F_L_A06__Layer_of_kde_of_noises_of_Seismic_Survey, na.rm = TRUE),
    noise_explosives = mean(F_L_A07__Layer_of_kde_of_noises_of_Explosives, na.rm = TRUE),
    noise_multibeam  = mean(F_L_A08__Layer_of_kde_of_noises_of_Multibeam_Echosounders, na.rm = TRUE),
    noise_subbottom  = mean(F_L_A09__Layer_of_kde_of_noises_of_Sub_bottom_Profilers, na.rm = TRUE),
    noise_ADD        = mean(F_L_A10__Layer_of_kde_of_noises_of_Acoustic_Deterrent_Device, na.rm = TRUE),
    noise_piling     = mean(F_L_A11__Layer_of_kde_of_noises_of_Piling, na.rm = TRUE),
    .groups = "drop"
  )

put_noise_names <- c("noise_MoD", "noise_seismic", "noise_explosives",
                     "noise_multibeam", "noise_subbottom", "noise_ADD", "noise_piling")

# noise_MoD
r <- rasterize(vect(put_noise_sf), ref_rast, field = "noise_MoD", fun = "mean")
r <- mask(r, smr_mask)
if (sum(!is.na(values(r)) & values(r) > 0) > 0) {
  raw_rasters[["noise_MoD"]] <- r; norm_rasters[["noise_MoD"]] <- normalise_01(r)
}

# noise_seismic
r <- rasterize(vect(put_noise_sf), ref_rast, field = "noise_seismic", fun = "mean")
r <- mask(r, smr_mask)
if (sum(!is.na(values(r)) & values(r) > 0) > 0) {
  raw_rasters[["noise_seismic"]] <- r; norm_rasters[["noise_seismic"]] <- normalise_01(r)
}

# noise_explosives
r <- rasterize(vect(put_noise_sf), ref_rast, field = "noise_explosives", fun = "mean")
r <- mask(r, smr_mask)
if (sum(!is.na(values(r)) & values(r) > 0) > 0) {
  raw_rasters[["noise_explosives"]] <- r; norm_rasters[["noise_explosives"]] <- normalise_01(r)
}

# noise_multibeam
r <- rasterize(vect(put_noise_sf), ref_rast, field = "noise_multibeam", fun = "mean")
r <- mask(r, smr_mask)
if (sum(!is.na(values(r)) & values(r) > 0) > 0) {
  raw_rasters[["noise_multibeam"]] <- r; norm_rasters[["noise_multibeam"]] <- normalise_01(r)
}

# noise_subbottom
r <- rasterize(vect(put_noise_sf), ref_rast, field = "noise_subbottom", fun = "mean")
r <- mask(r, smr_mask)
if (sum(!is.na(values(r)) & values(r) > 0) > 0) {
  raw_rasters[["noise_subbottom"]] <- r; norm_rasters[["noise_subbottom"]] <- normalise_01(r)
}

# noise_ADD
r <- rasterize(vect(put_noise_sf), ref_rast, field = "noise_ADD", fun = "mean")
r <- mask(r, smr_mask)
if (sum(!is.na(values(r)) & values(r) > 0) > 0) {
  raw_rasters[["noise_ADD"]] <- r; norm_rasters[["noise_ADD"]] <- normalise_01(r)
}

# noise_piling
r <- rasterize(vect(put_noise_sf), ref_rast, field = "noise_piling", fun = "mean")
r <- mask(r, smr_mask)
if (sum(!is.na(values(r)) & values(r) > 0) > 0) {
  raw_rasters[["noise_piling"]] <- r; norm_rasters[["noise_piling"]] <- normalise_01(r)
}

# Putuhena infrastructure (latest year)
cat("  Rasterising Putuhena infrastructure layers...\n")

r <- rasterize(vect(put_latest), ref_rast, field = "F_L_A62__Layer_of_CuOW_active_cable", fun = "mean")
r <- mask(r, smr_mask)
if (sum(!is.na(values(r)) & values(r) > 0) > 0) {
  raw_rasters[["infra_cables"]] <- r; norm_rasters[["infra_cables"]] <- normalise_01(r)
}

r <- rasterize(vect(put_latest), ref_rast,
               field = "F_L_A33__Layer_of_kernel_density_estimation_of_observed_wind_tur", fun = "mean")
r <- mask(r, smr_mask)
if (sum(!is.na(values(r)) & values(r) > 0) > 0) {
  raw_rasters[["infra_wind"]] <- r; norm_rasters[["infra_wind"]] <- normalise_01(r)
}

r <- rasterize(vect(put_latest), ref_rast,
               field = "F_L_A78__Layer_of_kde_of_infra_satellite_obs_of_tot_oil", fun = "mean")
r <- mask(r, smr_mask)
if (sum(!is.na(values(r)) & values(r) > 0) > 0) {
  raw_rasters[["infra_oil"]] <- r; norm_rasters[["infra_oil"]] <- normalise_01(r)
}

# Putuhena dark vessels (3-year mean)
cat("  Rasterising Putuhena dark vessel layers...\n")

put_dark_sf <- put_crop %>%
  group_by(x, y) %>%
  summarise(
    dark_fishing    = mean(F_L_A82__Layer_of_kde_of_vessel_satellite_obs_of_dark_fishing, na.rm = TRUE),
    dark_nonfishing = mean(F_L_A83__Layer_of_kde_of_vessel_satellite_obs_of_dark_non_fishin, na.rm = TRUE),
    .groups = "drop"
  )

r <- rasterize(vect(put_dark_sf), ref_rast, field = "dark_fishing", fun = "mean")
r <- mask(r, smr_mask)
if (sum(!is.na(values(r)) & values(r) > 0) > 0) {
  raw_rasters[["dark_fishing"]] <- r; norm_rasters[["dark_fishing"]] <- normalise_01(r)
}

r <- rasterize(vect(put_dark_sf), ref_rast, field = "dark_nonfishing", fun = "mean")
r <- mask(r, smr_mask)
if (sum(!is.na(values(r)) & values(r) > 0) > 0) {
  raw_rasters[["dark_nonfishing"]] <- r; norm_rasters[["dark_nonfishing"]] <- normalise_01(r)
}

cat(sprintf("\n  Raster pressure layers created: %d\n", length(norm_rasters)))
cat(sprintf("  Layers: %s\n", paste(names(norm_rasters), collapse = ", ")))


# --- 3c. Plot each layer: raw values (original units) + normalised ---

cat("\n----- 3c. Individual pressure layer plots -----\n")

dir.create("figures/threats/layers_raw", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/threats/layers_normalised", showWarnings = FALSE, recursive = TRUE)

# Plot raw layer
plot_raw <- function(lyr_name) {
  meta <- layer_meta[[lyr_name]]
  p <- ggplot() +
    geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
    tidyterra::geom_spatraster(data = raw_rasters[[lyr_name]]) +
    scale_fill_viridis_c(name = meta$units, na.value = "transparent", option = "inferno") +
    geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
    labs(title = sprintf("%s (%s)", lyr_name, meta$source),
         subtitle = sprintf("Original units: %s", meta$units)) +
    map_theme + coord_sf(crs = 27700)
  ggsave(sprintf("figures/threats/layers_raw/%s.png", lyr_name), p,
         width = 12, height = 10, units = "in", dpi = 300)
}

# Plot normalised layer
plot_norm <- function(lyr_name) {
  meta <- layer_meta[[lyr_name]]
  p <- ggplot() +
    geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
    tidyterra::geom_spatraster(data = norm_rasters[[lyr_name]]) +
    scale_fill_viridis_c(name = "Pressure\n(0–1)", na.value = "transparent", option = "magma",
                         limits = c(0, 1)) +
    geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
    labs(title = sprintf("%s — normalised (%s)", lyr_name, meta$source)) +
    map_theme + coord_sf(crs = 27700)
  ggsave(sprintf("figures/threats/layers_normalised/%s.png", lyr_name), p,
         width = 12, height = 10, units = "in", dpi = 300)
}

# Generate all plots
cat(sprintf("  Plotting %d layers (raw + normalised)...\n", length(norm_rasters)))

plot_raw("vessel_all"); plot_norm("vessel_all")
plot_raw("vessel_cargo"); plot_norm("vessel_cargo")
plot_raw("vessel_service"); plot_norm("vessel_service")
plot_raw("vessel_tanker"); plot_norm("vessel_tanker")
plot_raw("macroalgae_inns"); plot_norm("macroalgae_inns")
plot_raw("fishing_gfw"); plot_norm("fishing_gfw")
plot_raw("fishing_under12"); plot_norm("fishing_under12")

# Putuhena layers — only plot if they exist
if ("noise_MoD" %in% names(norm_rasters))        { plot_raw("noise_MoD"); plot_norm("noise_MoD") }
if ("noise_seismic" %in% names(norm_rasters))     { plot_raw("noise_seismic"); plot_norm("noise_seismic") }
if ("noise_explosives" %in% names(norm_rasters))  { plot_raw("noise_explosives"); plot_norm("noise_explosives") }
if ("noise_multibeam" %in% names(norm_rasters))   { plot_raw("noise_multibeam"); plot_norm("noise_multibeam") }
if ("noise_subbottom" %in% names(norm_rasters))   { plot_raw("noise_subbottom"); plot_norm("noise_subbottom") }
if ("noise_ADD" %in% names(norm_rasters))          { plot_raw("noise_ADD"); plot_norm("noise_ADD") }
if ("noise_piling" %in% names(norm_rasters))       { plot_raw("noise_piling"); plot_norm("noise_piling") }
if ("infra_cables" %in% names(norm_rasters))       { plot_raw("infra_cables"); plot_norm("infra_cables") }
if ("infra_wind" %in% names(norm_rasters))         { plot_raw("infra_wind"); plot_norm("infra_wind") }
if ("infra_oil" %in% names(norm_rasters))          { plot_raw("infra_oil"); plot_norm("infra_oil") }
if ("dark_fishing" %in% names(norm_rasters))       { plot_raw("dark_fishing"); plot_norm("dark_fishing") }
if ("dark_nonfishing" %in% names(norm_rasters))    { plot_raw("dark_nonfishing"); plot_norm("dark_nonfishing") }

cat("  All layer plots saved to figures/threats/layers_raw/ and layers_normalised/\n")


# --- 3d. Combine vector threat features, mask to SCAMP AOI ---

cat("\n----- 3d. Unified vector threat dataset -----\n")

# Waste discharge points
vec_waste <- waste_crop %>%
  select(geometry) %>%
  mutate(threat_type = "waste_discharge")

# Dredge spoil — polygons
vec_spoil_pg <- spoil_pg_crop %>%
  select(geometry) %>%
  mutate(threat_type = "dredge_spoil")

# Dredge spoil — points
vec_spoil_pt <- spoil_pt_crop %>%
  select(geometry) %>%
  mutate(threat_type = "dredge_spoil")

# Dredging
vec_dredge <- dredge_crop %>%
  select(geometry) %>%
  mutate(threat_type = "dredging")

# WW Wrecks
wrecks_crop <- st_set_geometry(wrecks_crop, "SHAPE")
vec_wrecks <- wrecks_crop %>%
  select(geometry = SHAPE) %>%
  mutate(threat_type = "ww_wrecks")

# Dumped munitions — polygons
vec_mun_pg <- mun_pg_crop %>%
  select(geometry) %>%
  mutate(threat_type = "dumped_munitions")

# Dumped munitions — points
vec_mun_pt <- mun_pt_crop %>%
  select(geometry) %>%
  mutate(threat_type = "dumped_munitions")

# UK Ports
vec_ports <- ports_crop %>%
  select(geometry) %>%
  mutate(threat_type = "ports")

# Combine and mask to SCAMP AOI
vec_threats_raw <- bind_rows(
  vec_waste, vec_spoil_pg, vec_spoil_pt, vec_dredge,
  vec_wrecks, vec_mun_pg, vec_mun_pt, vec_ports
)

# Mask: keep only features within or intersecting the SCAMP AOI
hits <- lengths(st_intersects(vec_threats_raw, scamp_aoi_sf)) > 0
vec_threats <- vec_threats_raw[hits, ]

# Clip polygon/line features to AOI boundary
vec_threats <- st_intersection(vec_threats, scamp_aoi_sf) %>%
  select(threat_type)

cat(sprintf("  Unified vector threats (masked to AOI): %d features\n", nrow(vec_threats)))
cat("  By type:\n")
print(table(vec_threats$threat_type))

st_write(vec_threats, "outputs/threats/vector_threats_unified.gpkg",
         delete_dsn = TRUE, quiet = TRUE)
cat("  Saved: outputs/threats/vector_threats_unified.gpkg\n")


# --- 3e. Cumulative pressure map ---

cat("\n----- 3e. Cumulative pressure map -----\n")

# Stack all normalised rasters and sum
pressure_stack <- rast(norm_rasters)
cumulative_pressure <- sum(pressure_stack, na.rm = TRUE)
cumulative_pressure <- mask(cumulative_pressure, smr_mask)

cat(sprintf("  Cumulative pressure range: %.2f – %.2f\n",
            minmax(cumulative_pressure)[1], minmax(cumulative_pressure)[2]))
cat(sprintf("  Number of pressure layers summed: %d\n", nlyr(pressure_stack)))

# Split vector threats into points and polygons for plotting
vec_pts <- vec_threats[st_geometry_type(vec_threats) %in% c("POINT", "MULTIPOINT"), ]
vec_polys <- vec_threats[st_geometry_type(vec_threats) %in%
                           c("POLYGON", "MULTIPOLYGON", "GEOMETRYCOLLECTION"), ]

p_cumulative <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  tidyterra::geom_spatraster(data = cumulative_pressure) +
  scale_fill_viridis_c(name = "Cumulative\npressure", na.value = "transparent", option = "magma") +
  {if (nrow(vec_polys) > 0) geom_sf(data = vec_polys, aes(colour = threat_type),
                                     fill = NA, linewidth = 0.5)} +
  {if (nrow(vec_pts) > 0) geom_sf(data = vec_pts, aes(colour = threat_type),
                                   size = 1.5, alpha = 0.8)} +
  scale_colour_brewer(name = "Vector threats", palette = "Set1") +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  labs(title = "Cumulative Marine Pressure Surface",
       subtitle = sprintf("%d raster layers (normalised 0–1) + vector threats",
                          nlyr(pressure_stack))) +
  map_theme + coord_sf(crs = 27700)
ggsave("figures/threats/cumulative_pressure.png", p_cumulative,
       width = 14, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/threats/cumulative_pressure.png\n")

# Save outputs
writeRaster(cumulative_pressure, "outputs/threats/cumulative_pressure_1km.tif",
            overwrite = TRUE)
cat("  Saved: outputs/threats/cumulative_pressure_1km.tif\n")

# Save raw layers as individual TIFs
dir.create("outputs/threats/layers_raw", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/threats/layers_normalised", showWarnings = FALSE, recursive = TRUE)

raw_names <- names(raw_rasters)
writeRaster(raw_rasters[["vessel_all"]], "outputs/threats/layers_raw/vessel_all.tif", overwrite = TRUE)
writeRaster(raw_rasters[["vessel_cargo"]], "outputs/threats/layers_raw/vessel_cargo.tif", overwrite = TRUE)
writeRaster(raw_rasters[["vessel_service"]], "outputs/threats/layers_raw/vessel_service.tif", overwrite = TRUE)
writeRaster(raw_rasters[["vessel_tanker"]], "outputs/threats/layers_raw/vessel_tanker.tif", overwrite = TRUE)
writeRaster(raw_rasters[["macroalgae_inns"]], "outputs/threats/layers_raw/macroalgae_inns.tif", overwrite = TRUE)
writeRaster(raw_rasters[["fishing_gfw"]], "outputs/threats/layers_raw/fishing_gfw.tif", overwrite = TRUE)
writeRaster(raw_rasters[["fishing_under12"]], "outputs/threats/layers_raw/fishing_under12.tif", overwrite = TRUE)

# Putuhena raw layers — only save if they exist
if ("noise_MoD" %in% raw_names) writeRaster(raw_rasters[["noise_MoD"]], "outputs/threats/layers_raw/noise_MoD.tif", overwrite = TRUE)
if ("noise_seismic" %in% raw_names) writeRaster(raw_rasters[["noise_seismic"]], "outputs/threats/layers_raw/noise_seismic.tif", overwrite = TRUE)
if ("noise_explosives" %in% raw_names) writeRaster(raw_rasters[["noise_explosives"]], "outputs/threats/layers_raw/noise_explosives.tif", overwrite = TRUE)
if ("noise_multibeam" %in% raw_names) writeRaster(raw_rasters[["noise_multibeam"]], "outputs/threats/layers_raw/noise_multibeam.tif", overwrite = TRUE)
if ("noise_subbottom" %in% raw_names) writeRaster(raw_rasters[["noise_subbottom"]], "outputs/threats/layers_raw/noise_subbottom.tif", overwrite = TRUE)
if ("noise_ADD" %in% raw_names) writeRaster(raw_rasters[["noise_ADD"]], "outputs/threats/layers_raw/noise_ADD.tif", overwrite = TRUE)
if ("noise_piling" %in% raw_names) writeRaster(raw_rasters[["noise_piling"]], "outputs/threats/layers_raw/noise_piling.tif", overwrite = TRUE)
if ("infra_cables" %in% raw_names) writeRaster(raw_rasters[["infra_cables"]], "outputs/threats/layers_raw/infra_cables.tif", overwrite = TRUE)
if ("infra_wind" %in% raw_names) writeRaster(raw_rasters[["infra_wind"]], "outputs/threats/layers_raw/infra_wind.tif", overwrite = TRUE)
if ("infra_oil" %in% raw_names) writeRaster(raw_rasters[["infra_oil"]], "outputs/threats/layers_raw/infra_oil.tif", overwrite = TRUE)
if ("dark_fishing" %in% raw_names) writeRaster(raw_rasters[["dark_fishing"]], "outputs/threats/layers_raw/dark_fishing.tif", overwrite = TRUE)
if ("dark_nonfishing" %in% raw_names) writeRaster(raw_rasters[["dark_nonfishing"]], "outputs/threats/layers_raw/dark_nonfishing.tif", overwrite = TRUE)

# Save normalised stack
writeRaster(pressure_stack, "outputs/threats/pressure_layers_normalised_1km.tif",
            overwrite = TRUE)
cat("  Saved all raw and normalised TIFs to outputs/threats/\n")


##=====================================================================
## STEP 4 — CATCHMENT-BASED LAND USE PRESSURE
##=====================================================================

cat("\n\n===== STEP 4: CATCHMENT LAND USE PRESSURE =====\n")


# --- 4a. Load and explore HLCM 2022 land cover raster ---

cat("\n----- 4a. Load HLCM 2022 (EUNIS Level 2) -----\n")

hlcm_path <- "data/raw/habitats/General/HLCM-2022-L2_SCOTLAND_TIFF_27700/HLCM_2022_EUNIS_LEVEL2.tif"
stopifnot("HLCM 2022 raster not found" = file.exists(hlcm_path))

hlcm <- rast(hlcm_path)
cat(sprintf("  Full raster: %d x %d cells\n", ncol(hlcm), nrow(hlcm)))
cat(sprintf("  Resolution: %.1f x %.1f m\n", res(hlcm)[1], res(hlcm)[2]))
cat(sprintf("  CRS: %s\n", crs(hlcm, describe = TRUE)$code))

hlcm_cats <- cats(hlcm)
cat(sprintf("  Has category table: %s\n", !is.null(hlcm_cats[[1]])))
if (!is.null(hlcm_cats[[1]])) {
  cat("  Category table:\n")
  print(hlcm_cats[[1]])
}

# Build the land cover lookup directly from the raster's own category table
cat_tbl <- hlcm_cats[[1]]
lc_lookup <- data.frame(
  lc_value   = cat_tbl$VALUE,
  eunis_code = cat_tbl$EUNIS_CODE,
  eunis_desc = cat_tbl$EUNIS_DESC,
  stringsAsFactors = FALSE
)
cat("  Land cover lookup (from raster):\n")
print(lc_lookup)

# Crop to extent of catchments (they extend inland beyond the AOI)
cefas_ext <- ext(vect(cefas_crop))
hlcm_crop <- crop(hlcm, cefas_ext)
cat(sprintf("  Cropped to catchment extent: %d x %d cells\n", ncol(hlcm_crop), nrow(hlcm_crop)))

hlcm_vals <- unique(values(hlcm_crop))
hlcm_vals <- sort(hlcm_vals[!is.na(hlcm_vals)])
cat(sprintf("  Unique land cover values in study area: %s\n", paste(hlcm_vals, collapse = ", ")))
cat(sprintf("  Number of unique classes: %d\n", length(hlcm_vals)))


# --- 4b. NWSS cross-tabulation: how does HLCM classify NWSS woodlands? ---

cat("\n----- 4b. NWSS vs HLCM woodland cross-tabulation -----\n")

habmos_path <- "data/raw/habitats/General/Habitat_Map_of_Scotland_2181603415421412055.gpkg"
stopifnot("HabMoS GPKG not found" = file.exists(habmos_path))

# Explore layers and fields
habmos_layers <- st_layers(habmos_path)
cat(sprintf("  HabMoS layers: %s\n", paste(habmos_layers$name, collapse = ", ")))
cat(sprintf("  Feature counts: %s\n", paste(habmos_layers$features, collapse = ", ")))

habmos_crop <- st_read(habmos_path, layer = habmos_layers$name[1],
                       wkt_filter = st_as_text(st_as_sfc(st_bbox(cefas_crop))),
                       quiet = TRUE) %>%
  st_make_valid()
cat(sprintf("  HabMoS features in study area: %d\n", nrow(habmos_crop)))
cat(sprintf("  HabMoS CRS: %s\n", st_crs(habmos_crop)$input))
cat(sprintf("  HabMoS columns: %s\n", paste(names(habmos_crop), collapse = ", ")))

# Print all column names and first rows to identify classification fields
cat("\n  First 10 rows (no geometry):\n")
print(head(st_drop_geometry(habmos_crop), 10))

# Look for classification columns
class_cols <- grep("type|category|interp|class|habitat|eunis|nvc|annex|phase",
                   names(habmos_crop), ignore.case = TRUE, value = TRUE)
cat(sprintf("\n  Candidate classification columns: %s\n",
            paste(class_cols, collapse = ", ")))

# Print unique values for each candidate column
if (length(class_cols) > 0) {
  cat("\n  Unique values per candidate column:\n")
  cat(sprintf("  %s: %s\n", class_cols[1],
              paste(head(unique(habmos_crop[[class_cols[1]]]), 30), collapse = ", ")))
  if (length(class_cols) > 1) {
    cat(sprintf("  %s: %s\n", class_cols[2],
                paste(head(unique(habmos_crop[[class_cols[2]]]), 30), collapse = ", ")))
  }
  if (length(class_cols) > 2) {
    cat(sprintf("  %s: %s\n", class_cols[3],
                paste(head(unique(habmos_crop[[class_cols[3]]]), 30), collapse = ", ")))
  }
  if (length(class_cols) > 3) {
    cat(sprintf("  %s: %s\n", class_cols[4],
                paste(head(unique(habmos_crop[[class_cols[4]]]), 30), collapse = ", ")))
  }
}

# >>> AFTER RUNNING: set habmos_class_col to the column that distinguishes
#     native woodland from plantation (e.g. "HABITAT", "EUNIS", "INTERP_HAB") <<<
habmos_class_col <- "HABITAT_CODE"

# Filter HabMoS to woodland polygons only (HABITAT_CODE starting with "G")
habmos_wood <- habmos_crop[grepl("^G", habmos_crop[[habmos_class_col]]), ]
cat(sprintf("  HabMoS woodland polygons in study area: %d\n", nrow(habmos_wood)))
cat(sprintf("  Woodland codes present: %s\n",
            paste(sort(unique(habmos_wood[[habmos_class_col]])), collapse = ", ")))

# Extract HLCM categories under woodland polygons
habmos_hlcm <- terra::extract(hlcm_crop, vect(habmos_wood))
hlcm_col <- names(habmos_hlcm)[2]
habmos_hlcm$habmos_code <- habmos_wood[[habmos_class_col]][habmos_hlcm$ID]

habmos_hlcm <- habmos_hlcm[!is.na(habmos_hlcm[[hlcm_col]]), ]
cat(sprintf("  Pixels extracted under woodland polygons: %d\n", nrow(habmos_hlcm)))

# Cross-tabulation: HabMoS woodland code vs HLCM EUNIS class
# terra::extract on a categorical raster returns the EUNIS codes directly
cat("\n  ===== CROSS-TABULATION: HabMoS woodland vs HLCM EUNIS =====\n")
xtab <- table(
  HabMoS = habmos_hlcm$habmos_code,
  HLCM   = habmos_hlcm[[hlcm_col]]
)
print(xtab)

# Row percentages — how each HabMoS woodland type maps to HLCM
cat("\n  Row proportions (% of each HabMoS type mapped to each HLCM class):\n")
xtab_prop <- round(100 * prop.table(xtab, margin = 1), 1)
print(xtab_prop)

# Focused view: just the HLCM woodland columns
wood_hlcm <- c("G1", "G3.4", "G3.F", "G4", "G5")
wood_cols <- colnames(xtab) %in% wood_hlcm
if (sum(wood_cols) > 0) {
  cat("\n  HLCM woodland classes only:\n")
  print(xtab[, wood_cols, drop = FALSE])
  cat("\n  Row proportions (woodland classes only):\n")
  wood_prop <- round(100 * prop.table(xtab[, wood_cols, drop = FALSE], margin = 1), 1)
  print(wood_prop)
}

# SUMMARY OF HABMOS vs HLCM WOODLAND CROSS-TABULATION
#
# G3.4 (Scots pine) is entirely absent from HLCM in this study area.
#
# HabMoS G3.F (conifer plantation) aligns well with HLCM G3.F (~85% of
# woodland pixels agree). However, large areas within HabMoS plantation
# boundaries are classified by HLCM as non-woodland: bare land, windthrow,
# grassland, heathland — likely felled/restocking areas captured by the
# 2022 satellite pass but not yet updated in HabMoS.
#
# Native broadleaved types (G1.xx) align well with HLCM G1 (79–87%).
#
# HabMoS G3.F2 (exotic conifer plantations) is misclassified by HLCM:
# 64% mapped to G1 (broadleaved), only 21% to G3.F.
#
# HabMoS G4.F (mixed forestry plantations) is also misclassified:
# 70% mapped to G1, only 9% to G3.F.
#
# HLCM 2022 is a uniform wall-to-wall remote sensing product. HabMoS is a
# composite of surveys from multiple years (2017–2023). HLCM better reflects
# current land cover state; HabMoS better reflects intended land management.
#
# DECISION: proceed with HLCM 2022 for catchment pressure analysis.
# G3.F vs G1 separation is adequate for the main plantation/native distinction.
# Exotic conifers (G3.F2) and mixed plantations (G4.F) are under-represented
# in HLCM's G3.F class — land-based pressure from these areas may be slightly
# underestimated as a result.


# --- 4b-ii. HabMoS vs HLCM grassland cross-tabulation ---

cat("\n----- 4b-ii. HabMoS vs HLCM grassland cross-tabulation -----\n")

# Filter HabMoS to grassland polygons (HABITAT_CODE starting with "E")
habmos_grass <- habmos_crop[grepl("^E", habmos_crop[[habmos_class_col]]), ]
cat(sprintf("  HabMoS grassland polygons in study area: %d\n", nrow(habmos_grass)))
cat(sprintf("  Grassland codes present: %s\n",
            paste(sort(unique(habmos_grass[[habmos_class_col]])), collapse = ", ")))

# Extract HLCM categories under grassland polygons
habmos_grass_hlcm <- terra::extract(hlcm_crop, vect(habmos_grass))
grass_hlcm_col <- names(habmos_grass_hlcm)[2]
habmos_grass_hlcm$habmos_code <- habmos_grass[[habmos_class_col]][habmos_grass_hlcm$ID]

habmos_grass_hlcm <- habmos_grass_hlcm[!is.na(habmos_grass_hlcm[[grass_hlcm_col]]), ]
cat(sprintf("  Pixels extracted under grassland polygons: %d\n", nrow(habmos_grass_hlcm)))

# Also add HabMoS HABITAT_NAME for readability
habmos_grass_hlcm$habmos_name <- habmos_grass[["HABITAT_NAME"]][habmos_grass_hlcm$ID]

# Cross-tabulation: HabMoS grassland code vs HLCM EUNIS class
cat("\n  ===== CROSS-TABULATION: HabMoS grassland vs HLCM EUNIS =====\n")
xtab_grass <- table(
  HabMoS = habmos_grass_hlcm$habmos_code,
  HLCM   = habmos_grass_hlcm[[grass_hlcm_col]]
)
print(xtab_grass)

cat("\n  Row proportions (% of each HabMoS grassland type mapped to each HLCM class):\n")
xtab_grass_prop <- round(100 * prop.table(xtab_grass, margin = 1), 1)
print(xtab_grass_prop)

# Focused view: just the HLCM grassland columns
grass_hlcm_codes <- c("E1", "E2", "E3", "E4", "E5")
grass_cols <- colnames(xtab_grass) %in% grass_hlcm_codes
if (sum(grass_cols) > 0) {
  cat("\n  HLCM grassland classes only:\n")
  print(xtab_grass[, grass_cols, drop = FALSE])
  cat("\n  Row proportions (grassland classes only):\n")
  grass_prop <- round(100 * prop.table(xtab_grass[, grass_cols, drop = FALSE], margin = 1), 1)
  print(grass_prop)
}

# Reverse view: column proportions — of all pixels HLCM calls E2 (or E3),
# what HabMoS sub-codes are present? Shows dominant HabMoS composition within each HLCM class.
cat("\n  Column proportions (% of each HLCM class composed of each HabMoS type):\n")
xtab_grass_colprop <- round(100 * prop.table(xtab_grass, margin = 2), 1)
print(xtab_grass_colprop)

if (sum(grass_cols) > 0) {
  cat("\n  Column proportions (grassland HLCM classes only):\n")
  grass_colprop <- round(100 * prop.table(xtab_grass[, grass_cols, drop = FALSE], margin = 2), 1)
  print(grass_colprop)
}

cat("\n  HabMoS grassland code descriptions (for reference):\n")
grass_names <- unique(habmos_grass_hlcm[, c("habmos_code", "habmos_name")])
grass_names <- grass_names[order(grass_names$habmos_code), ]
print(grass_names, row.names = FALSE)

# SUMMARY OF HABMOS vs HLCM GRASSLAND CROSS-TABULATION
#
# HLCM E2 (mesic grassland) is a reasonable proxy for improved grassland:
# 86.8% of E2.6 (improved/fertilised) and 92.3% of E2.112 (ryegrass) map to E2.
# But E2 also contains some semi-natural types (e.g. E2.111 Cynosurus-Centaurea).
#
# E2.11 (unbroken pastures) leaks: only 54% to E2, 32% to E3.
# So some agricultural pressure hides in HLCM E3.
#
# HLCM E1 and E3 are predominantly semi-natural (acid/wet grasslands).
# E1 and E3 bleed into each other for some dry grassland types (E1.72, E1.73).
#
# DECISION: HabMoS coverage too sparse for wall-to-wall grassland reclassification.
# Use James Hutton Land Capability for Agriculture (LCA) instead — full coverage,
# land capability classes 1-7 indicate agricultural intensity.
# HabMoS woodland correction (G3.F2, G4.F) retained — good coverage for plantations.


# --- 4b-iii. Reclassify HLCM using LCA (grasslands) and HabMoS (woodlands) ---

cat("\n----- 4b-iii. HLCM reclassification -----\n")

hlcm_int <- as.numeric(hlcm_crop)
NAflag(hlcm_int) <- -32768
cat(sprintf("  HLCM integer raster: %d x %d cells\n", ncol(hlcm_int), nrow(hlcm_int)))

e1_val  <- lc_lookup$lc_value[lc_lookup$eunis_code == "E1"]
e2_val  <- lc_lookup$lc_value[lc_lookup$eunis_code == "E2"]
e3_val  <- lc_lookup$lc_value[lc_lookup$eunis_code == "E3"]
g1_val  <- lc_lookup$lc_value[lc_lookup$eunis_code == "G1"]
g3f_val <- lc_lookup$lc_value[lc_lookup$eunis_code == "G3.F"]

# --- Load James Hutton Land Capability for Agriculture (LCA) ---

cat("\n  Loading James Hutton LCA 1:50k...\n")
lca_path <- "data/raw/physical/Hutton_LCA50K_OpenData/LCA_50K.shp"
lca <- st_read(lca_path, quiet = TRUE)
cat(sprintf("  LCA features: %d\n", nrow(lca)))
cat(sprintf("  LCA CRS: %s (EPSG:%s)\n", st_crs(lca)$Name, st_crs(lca)$epsg))
cat("  LCA column names:\n")
print(names(lca))
cat("  Unique LCCODE values:\n")
print(sort(unique(lca$LCCODE)))

lca_crop <- st_crop(lca, st_bbox(cefas_crop)) %>% st_make_valid()
lca_crop <- lca_crop[lca_crop$Water == "No", ]
cat(sprintf("  LCA features after crop (excl. water): %d\n", nrow(lca_crop)))

# --- Grassland reclassification using LCA ---
# Extract HLCM cell values under each LCA class group.
# Only grassland cells (E1, E2, E3) are reclassified.
# LCCODE is numeric: integer part is the LCA class (e.g. 3.1 = class 3, division 1).

grass_vals <- c(e1_val, e2_val, e3_val)

lca_class <- floor(lca_crop$LCCODE)

# LCA 1-3: Likely improved pasture (value 29)
lca_13 <- lca_crop[which(lca_class >= 1 & lca_class <= 3), ]
cat(sprintf("\n  LCA class 1-3 polygons: %d\n", nrow(lca_13)))
if (nrow(lca_13) > 0) {
  ex_13 <- terra::extract(hlcm_int, vect(lca_13), cells = TRUE)
  hits_13 <- ex_13$cell[ex_13[[2]] %in% grass_vals]
  hlcm_int[hits_13] <- 29
  cat(sprintf("  Grassland pixels reclassified (LCA 1-3): %d\n", length(hits_13)))
}

# LCA 4: Likely improved pasture, marginal (value 30)
lca_4 <- lca_crop[which(lca_class == 4), ]
cat(sprintf("  LCA class 4 polygons: %d\n", nrow(lca_4)))
if (nrow(lca_4) > 0) {
  ex_4 <- terra::extract(hlcm_int, vect(lca_4), cells = TRUE)
  hits_4 <- ex_4$cell[ex_4[[2]] %in% grass_vals]
  hlcm_int[hits_4] <- 30
  cat(sprintf("  Grassland pixels reclassified (LCA 4): %d\n", length(hits_4)))
}

# LCA 5: Possible improved or rough grazing (value 31)
lca_5 <- lca_crop[which(lca_class == 5), ]
cat(sprintf("  LCA class 5 polygons: %d\n", nrow(lca_5)))
if (nrow(lca_5) > 0) {
  ex_5 <- terra::extract(hlcm_int, vect(lca_5), cells = TRUE)
  hits_5 <- ex_5$cell[ex_5[[2]] %in% grass_vals]
  hlcm_int[hits_5] <- 31
  cat(sprintf("  Grassland pixels reclassified (LCA 5): %d\n", length(hits_5)))
}

# LCA 6-7: Likely semi-natural/rough grazing (value 32)
lca_67 <- lca_crop[which(lca_class >= 6 & lca_class <= 7), ]
cat(sprintf("  LCA class 6-7 polygons: %d\n", nrow(lca_67)))
if (nrow(lca_67) > 0) {
  ex_67 <- terra::extract(hlcm_int, vect(lca_67), cells = TRUE)
  hits_67 <- ex_67$cell[ex_67[[2]] %in% grass_vals]
  hlcm_int[hits_67] <- 32
  cat(sprintf("  Grassland pixels reclassified (LCA 6-7): %d\n", length(hits_67)))
}

# --- Woodland reclassification (HabMoS) ---
# Only reclassifying HLCM G1 (broadleaved) pixels where HabMoS identifies plantations.

g4f_desc <- na.omit(habmos_crop[["HABITAT_NAME"]][habmos_crop[[habmos_class_col]] == "G4.F"])[1]

habmos_g3f2 <- habmos_crop[which(habmos_crop[[habmos_class_col]] == "G3.F2"), ]
cat(sprintf("\n  HabMoS G3.F2 polygons: %d\n", nrow(habmos_g3f2)))
if (nrow(habmos_g3f2) > 0) {
  ex_g3f2 <- terra::extract(hlcm_int, vect(habmos_g3f2), cells = TRUE)
  hits_g3f2 <- ex_g3f2$cell[ex_g3f2[[2]] == g1_val]
  hlcm_int[hits_g3f2] <- g3f_val
  cat(sprintf("  G3.F2 pixels reclassified to G3.F: %d\n", length(hits_g3f2)))
}

habmos_g4f <- habmos_crop[which(habmos_crop[[habmos_class_col]] == "G4.F"), ]
cat(sprintf("  HabMoS G4.F polygons: %d\n", nrow(habmos_g4f)))
if (nrow(habmos_g4f) > 0) {
  ex_g4f <- terra::extract(hlcm_int, vect(habmos_g4f), cells = TRUE)
  hits_g4f <- ex_g4f$cell[ex_g4f[[2]] == g1_val]
  hlcm_int[hits_g4f] <- 33
  cat(sprintf("  G4.F pixels reclassified: %d\n", length(hits_g4f)))
}

# Build updated lookup table with new classes
lc_lookup_reclass <- rbind(
  lc_lookup,
  data.frame(
    lc_value   = c(29, 30, 31, 32, 33),
    eunis_code = c("E.lca13", "E.lca4", "E.lca5", "E.lca67", "G4.F"),
    eunis_desc = c("Likely improved pasture (LCA 1-3)",
                   "Likely improved pasture, marginal (LCA 4)",
                   "Possible improved or rough grazing (LCA 5)",
                   "Likely semi-natural/rough grazing (LCA 6-7)",
                   g4f_desc),
    stringsAsFactors = FALSE
  )
)
cat("\n  Reclassified classes added to lookup:\n")
print(lc_lookup_reclass[lc_lookup_reclass$lc_value >= 29, ])

hlcm_reclass <- hlcm_int
cat(sprintf("\n  Reclassified HLCM ready: %d unique values\n",
            length(unique(na.omit(values(hlcm_reclass))))))

# --- Grassland intensity frequency table (EUNIS L2 x LCA class) ---

cat("\n  Grassland intensity cross-tabulation (EUNIS x LCA, area in hectares):\n")

# Build frequency table from the original (pre-reclassification) HLCM
# crossed with LCA classes, for grassland cells only.
grass_int_fresh <- as.numeric(hlcm_crop)
pixel_ha <- prod(res(hlcm_crop)) / 10000

ex_all_13 <- terra::extract(grass_int_fresh, vect(lca_13), cells = TRUE)
ex_all_13$lca_group <- "LCA 1-3"
ex_all_4  <- terra::extract(grass_int_fresh, vect(lca_4), cells = TRUE)
ex_all_4$lca_group  <- "LCA 4"
ex_all_5  <- terra::extract(grass_int_fresh, vect(lca_5), cells = TRUE)
ex_all_5$lca_group  <- "LCA 5"
ex_all_67 <- terra::extract(grass_int_fresh, vect(lca_67), cells = TRUE)
ex_all_67$lca_group <- "LCA 6-7"

ex_all <- rbind(ex_all_13, ex_all_4, ex_all_5, ex_all_67)
names(ex_all)[2] <- "lc_value"
ex_all <- ex_all[ex_all$lc_value %in% c(e1_val, e2_val, e3_val), ]
ex_all <- left_join(ex_all, lc_lookup[, c("lc_value", "eunis_code", "eunis_desc")],
                    by = "lc_value")

grass_lca_freq <- ex_all %>%
  group_by(eunis_code, eunis_desc, lca_group) %>%
  summarise(n_pixels = n(), .groups = "drop") %>%
  mutate(area_ha = round(n_pixels * pixel_ha, 1)) %>%
  arrange(eunis_code, lca_group)
print(grass_lca_freq, n = 50)

# Save grassland intensity raster
dir.create("outputs/catchments", showWarnings = FALSE, recursive = TRUE)
grass_intensity <- hlcm_reclass
grass_mask <- values(grass_intensity) %in% c(e1_val, e2_val, e3_val, 29, 30, 31, 32)
values(grass_intensity)[!grass_mask] <- NA
writeRaster(grass_intensity, "outputs/catchments/grassland_intensity.tif",
            overwrite = TRUE)
cat("  Saved: outputs/catchments/grassland_intensity.tif\n")

st_write(st_transform(cefas_crop[1, ], 4326),
         "outputs/catchments/cefas_catchments_crop.shp",
         delete_dsn = TRUE, quiet = TRUE)
shp_files <- list.files("outputs/catchments",
                        pattern = "cefas_catchments_crop\\.(shp|shx|dbf|prj|cpg)$",
                        full.names = TRUE)
zip("outputs/catchments/cefas_catchments_crop.zip", shp_files, flags = "-j")
cat("  Saved: outputs/catchments/cefas_catchments_crop.zip\n")


# --- 4c. Extract land cover per CEFAS catchment ---

cat("\n----- 4c. Extract land cover per catchment -----\n")

hlcm_extract <- terra::extract(hlcm_reclass, vect(cefas_crop))
cat(sprintf("  Total pixels extracted: %d\n", nrow(hlcm_extract)))

lc_colname <- names(hlcm_extract)[2]
cat(sprintf("  Land cover column name: %s\n", lc_colname))

hlcm_extract$S_Catchmt <- cefas_crop$S_Catchmt[hlcm_extract$ID]
hlcm_extract$S_CatNAME <- cefas_crop$S_CatNAME[hlcm_extract$ID]

hlcm_extract <- hlcm_extract[!is.na(hlcm_extract[[lc_colname]]), ]
cat(sprintf("  Non-NA pixels: %d\n", nrow(hlcm_extract)))

# Integer raster returns numeric values — join to lookup for EUNIS codes/descriptions
names(hlcm_extract)[names(hlcm_extract) == lc_colname] <- "lc_value"

lc_counts <- hlcm_extract %>%
  group_by(S_Catchmt, S_CatNAME, lc_value) %>%
  summarise(n_pixels = n(), .groups = "drop") %>%
  left_join(lc_lookup_reclass, by = "lc_value")

catch_totals <- lc_counts %>%
  group_by(S_Catchmt, S_CatNAME) %>%
  summarise(total_pixels = sum(n_pixels), .groups = "drop")

lc_props <- lc_counts %>%
  left_join(catch_totals, by = c("S_Catchmt", "S_CatNAME")) %>%
  mutate(proportion = n_pixels / total_pixels)

cat(sprintf("  Catchments with land cover data: %d\n", nrow(catch_totals)))
cat("  Pixel counts per catchment:\n")
print(catch_totals)

lc_freq <- lc_counts %>%
  group_by(eunis_code, eunis_desc) %>%
  summarise(total_pixels = sum(n_pixels), .groups = "drop") %>%
  arrange(desc(total_pixels)) %>%
  mutate(pct = round(100 * total_pixels / sum(total_pixels), 1))
cat("\n  Land cover class frequency (all catchments):\n")
print(lc_freq)


# --- 4d. SCIMAP diffuse pollution risk ---

cat("\n----- 4d. SCIMAP diffuse pollution risk -----\n")

# ---- STEP 1: SCIMAP-derived source risk weightings ----

## Crosswalk: Space Intelligence EUNIS classes to SCIMAP risk weightings
## Reference: SCIMAP parameter values (Reaney et al., scimap.org.uk)
## Adjustments noted where EUNIS classes do not directly map to CEH LCM classes
##
## I1      Arable land and market gardens          → 1.00  (SCIMAP: Arable and Horticulture — direct match)
## E.lca13 Likely improved pasture (LCA 1-3)       → 0.30  (SCIMAP: Improved Grassland — direct match)
## E.lca4  Likely improved pasture marginal         → 0.30  (SCIMAP: Improved Grassland — still fertilised/reseeded)
## E.lca5  Possible improved or rough grazing       → 0.20  (No SCIMAP equivalent — interpolated between Improved 0.30 and Rough 0.15)
## E.lca67 Likely semi-natural/rough grazing        → 0.15  (SCIMAP: Rough grassland — direct match)
## E1      Dry grasslands (unclassified by LCA)     → 0.15  (SCIMAP: Acid grassland — unimproved)
## E2      Mesic grasslands (residual after LCA)    → 0.15  (SCIMAP: Neutral Grassland — residual unclassified)
## E3      Seasonally wet and wet grasslands        → 0.15  (SCIMAP: Fen, Marsh and Swamp)
## E4      Alpine and subalpine grasslands          → 0.05  (SCIMAP: Montane Habitats — direct match)
## E5      Woodland fringes and clearings           → 0.10  (Below SCIMAP Neutral Grassland 0.15 — semi-wooded, lower exposure)
## G1      Broadleaved deciduous woodland           → 0.05  (SCIMAP: Broadleaved woodland — direct match)
## G3.F    Highly artificial coniferous plantations → 0.10  (Above SCIMAP Coniferous 0.05 — intensive forestry practices
##                                                           including fertilisation during planting, litter acidification,
##                                                           soil disturbance during harvesting)
## G4      Mixed deciduous and coniferous woodland  → 0.05  (SCIMAP: Broadleaved woodland — treat as semi-natural)
## G4.F    Mixed forestry plantations               → 0.05  (SCIMAP: Coniferous Woodland — not intensively managed like G3.F)
## G5      Lines of trees, felled, early-stage      → 0.10  (Recently felled woodland has exposed bare soil and brash,
##                                                           temporary elevated erosion/runoff risk)
## F4      Temperate shrub heathland                → 0.05  (SCIMAP: Heather — direct match)
## F3      Temperate scrub                          → 0.05  (SCIMAP: nearest equivalent is Heather — low pressure)
## F2      Arctic/alpine scrub                      → 0.05  (SCIMAP: Montane Habitats)
## F9      Riverine scrub                           → 0.05  (Low pressure natural feature)
## D1      Raised and blanket bogs                  → 0.05  (SCIMAP: Bog — direct match)
## D2      Valley mires and poor fens               → 0.05  (Natural feature — nutrient dynamics are baseline not additive)
## J       Built-up                                 → 0.01  (SCIMAP: Urban/Suburban — direct match)
## A2.5    Coastal saltmarshes                      → 0     (SCIMAP: Saltmarsh — not a pollution source)
## B1      Coastal dunes and sandy shores           → 0     (SCIMAP: Supra-littoral Sediment)
## B2      Coastal shingle                          → 0     (SCIMAP: Supra-littoral Sediment)
## B3      Rock cliffs and shores                   → 0     (SCIMAP: Supra-littoral Rock)
## C       Water                                    → 0     (SCIMAP: Freshwater)
## O       Bare land                                → 0     (Context-dependent — includes natural bare ground)
## OW      Windthrow                                → 0.10  (Exposed soil post-windthrow — known erosion risk)
## H2      Screes                                   → 0     (SCIMAP: Inland Rock)
## H3      Inland cliffs and outcrops               → 0     (SCIMAP: Inland Rock)
## NA      Unclassified                             → NA    (Exclude from analysis)

scimap_lookup <- data.frame(
  lc_value = lc_lookup_reclass$lc_value,
  eunis_code = lc_lookup_reclass$eunis_code,
  stringsAsFactors = FALSE
)

risk_map <- c(
  "I1" = 1.00, "E.lca13" = 0.30, "E.lca4" = 0.30, "E.lca5" = 0.20,
  "E.lca67" = 0.15, "E1" = 0.15, "E2" = 0.15, "E3" = 0.15, "E4" = 0.05,
  "E5" = 0.10, "G1" = 0.05, "G3.F" = 0.10, "G3.4" = 0.05, "G4" = 0.05,
  "G4.F" = 0.05, "G5" = 0.10, "F4" = 0.05, "F3" = 0.05, "F2" = 0.05,
  "F9" = 0.05, "D1" = 0.05, "D2" = 0.05, "D4" = 0.05, "J" = 0.01,
  "A2.5" = 0, "B1" = 0, "B2" = 0, "B3" = 0, "C" = 0, "O" = 0,
  "OW" = 0.10, "H2" = 0, "H3" = 0
)

scimap_lookup$risk_weight <- risk_map[scimap_lookup$eunis_code]

cat("  SCIMAP source risk weightings:\n")
print(scimap_lookup[order(-scimap_lookup$risk_weight), ], row.names = FALSE)

# Build classify matrix: from lc_value → risk_weight
rcl_risk <- as.matrix(scimap_lookup[!is.na(scimap_lookup$risk_weight),
                                     c("lc_value", "risk_weight")])
source_risk <- terra::classify(hlcm_reclass, rcl_risk, others = NA)
cat(sprintf("  Source risk raster: %d x %d cells\n", ncol(source_risk), nrow(source_risk)))
cat(sprintf("  Value range: %.2f – %.2f\n",
            min(values(source_risk), na.rm = TRUE),
            max(values(source_risk), na.rm = TRUE)))

plot(source_risk)
# ---- STEP 2: TWI — hydrological connectivity pathway ----
# TWI captures hydrological connectivity based on upslope contributing area and
# local slope. It does not distinguish nutrient pathways (N, P, sediment).
# This is a simplification but provides complete spatial coverage where
# Hutton soil layers have gaps.

cat("\n  --- Step 2: Computing TWI from GEBCO DEM ---\n")

gebco_path <- "data/raw/physical/GEBCO_2026/GEBCO_2026.nc"
gebco_raw <- rast(gebco_path)
cat(sprintf("  GEBCO raw: %d x %d, CRS: %s\n",
            ncol(gebco_raw), nrow(gebco_raw), crs(gebco_raw, describe = TRUE)$name))

cefas_wgs84 <- st_transform(cefas_crop, crs(gebco_raw))
gebco_crop <- crop(gebco_raw, ext(vect(cefas_wgs84)))
cat(sprintf("  GEBCO cropped: %d x %d\n", ncol(gebco_crop), nrow(gebco_crop)))

gebco_bng <- project(gebco_crop, crs(source_risk), method = "bilinear")
cat(sprintf("  GEBCO reprojected to BNG: %d x %d\n", ncol(gebco_bng), nrow(gebco_bng)))

cefas_union <- st_union(cefas_crop)
gebco_bng[gebco_bng < 0] <- 0
gebco_masked <- mask(gebco_bng, vect(cefas_union))
plot(gebco_masked)


dem_path <- "outputs/catchments/dem_catchments.tif"
writeRaster(gebco_masked, dem_path, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
cat(sprintf("  DEM saved: %s\n", dem_path))

dem_vals <- values(gebco_masked)[, 1]
cat(sprintf("  DEM range: %.1f – %.1f m\n",
            min(dem_vals, na.rm = TRUE), max(dem_vals, na.rm = TRUE)))
cat(sprintf("  DEM cells below 0m (sea): %d / %d (%.1f%%)\n",
            sum(dem_vals < 0, na.rm = TRUE), sum(!is.na(dem_vals)),
            100 * sum(dem_vals < 0, na.rm = TRUE) / sum(!is.na(dem_vals))))
cat(sprintf("  DEM quantiles (m): 5%%=%.1f, 25%%=%.1f, 50%%=%.1f, 75%%=%.1f, 95%%=%.1f\n",
            quantile(dem_vals, 0.05, na.rm = TRUE),
            quantile(dem_vals, 0.25, na.rm = TRUE),
            quantile(dem_vals, 0.50, na.rm = TRUE),
            quantile(dem_vals, 0.75, na.rm = TRUE),
            quantile(dem_vals, 0.95, na.rm = TRUE)))

filled_path <- "outputs/catchments/dem_filled.tif"
wbt_fill_depressions(dem = dem_path, output = filled_path)
cat(sprintf("  Filled DEM: %s\n", filled_path))

sca_path <- "outputs/catchments/sca.tif"
wbt_d_inf_flow_accumulation(input = filled_path, output = sca_path, out_type = "Specific Contributing Area")
cat(sprintf("  Specific contributing area: %s\n", sca_path))

slope_path <- "outputs/catchments/slope.tif"
wbt_slope(dem = filled_path, output = slope_path)
cat(sprintf("  Slope: %s\n", slope_path))

twi_path <- "outputs/catchments/twi.tif"
wbt_wetness_index(sca = sca_path, slope = slope_path, output = twi_path)
cat(sprintf("  TWI computed: %s\n", twi_path))

twi_raw <- rast(twi_path)
twi_vals <- values(twi_raw)[, 1]
cat(sprintf("  TWI raw range: %.3f – %.3f\n",
            min(twi_vals, na.rm = TRUE), max(twi_vals, na.rm = TRUE)))
cat(sprintf("  TWI non-NA cells: %d / %d\n",
            sum(!is.na(twi_vals)), length(twi_vals)))
cat(sprintf("  TWI quantiles: 5%%=%.2f, 25%%=%.2f, 50%%=%.2f, 75%%=%.2f, 95%%=%.2f, 99%%=%.2f\n",
            quantile(twi_vals, 0.05, na.rm = TRUE),
            quantile(twi_vals, 0.25, na.rm = TRUE),
            quantile(twi_vals, 0.50, na.rm = TRUE),
            quantile(twi_vals, 0.75, na.rm = TRUE),
            quantile(twi_vals, 0.95, na.rm = TRUE),
            quantile(twi_vals, 0.99, na.rm = TRUE)))
twi_min <- min(twi_vals, na.rm = TRUE)
twi_max <- max(twi_vals, na.rm = TRUE)
twi_norm <- (twi_raw - twi_min) / (twi_max - twi_min)

twi_norm_path <- "outputs/catchments/twi_normalised.tif"
writeRaster(twi_norm, twi_norm_path, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
cat(sprintf("  TWI normalised (0-1) saved: %s\n", twi_norm_path))
cat(sprintf("  TWI normalised range: %.3f – %.3f\n",
            min(values(twi_norm), na.rm = TRUE),
            max(values(twi_norm), na.rm = TRUE)))

twi_resampled <- clamp(resample(twi_norm, source_risk, method = "bilinear"), lower = 0, upper = 1)
twi_res_vals <- values(twi_resampled)[, 1]
src_vals <- values(source_risk)[, 1]
cat(sprintf("  TWI resampled to HLCM grid: %d x %d\n",
            ncol(twi_resampled), nrow(twi_resampled)))
cat(sprintf("  TWI resampled: non-NA = %d, source_risk non-NA = %d\n",
            sum(!is.na(twi_res_vals)), sum(!is.na(src_vals))))
cat(sprintf("  TWI coverage of source_risk cells: %.1f%%\n",
            100 * sum(!is.na(twi_res_vals) & !is.na(src_vals)) / sum(!is.na(src_vals))))


# ---- STEP 3: Diffuse pollution risk ----

cat("\n  --- Step 3: Diffuse pollution risk = source × TWI ---\n")

diffuse_pollution_risk <- source_risk * twi_resampled

cat(sprintf("  Diffuse pollution risk: non-NA = %d, range = %.4f – %.4f\n",
            sum(!is.na(values(diffuse_pollution_risk))),
            min(values(diffuse_pollution_risk), na.rm = TRUE),
            max(values(diffuse_pollution_risk), na.rm = TRUE)))

writeRaster(diffuse_pollution_risk, "outputs/catchments/diffuse_pollution_risk.tif",
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))
cat("  Saved: outputs/catchments/diffuse_pollution_risk.tif\n")


# ---- STEP 3b (optional): Hutton validation ----
# If Hutton soil layers exist, compute Spearman correlation with TWI in overlap areas.

cat("\n  --- Step 3b: Hutton layer validation (optional) ---\n")

hutton_validated <- FALSE

erosion_path <- "data/raw/physical/Hutton_SoilErosionRisk_and_LUT/SoilErosionRisk/SoilErosionRisk_v1.tif"
leach_path <- "data/raw/physical/Hutton_SoilLeachingPotential/SoilLeachingPotential/SoilLeachingPotential_v1.shp"
runoff_path <- "data/raw/physical/Hutton_SoilRunoffRisk/SoilRunoffRisk/SoilRunoffRisk_v1.shp"

if (file.exists(erosion_path) & file.exists(leach_path) & file.exists(runoff_path)) {

  cat("  Hutton layers found — computing Spearman correlations with TWI\n")

  erosion_raw <- rast(erosion_path)
  lut_path <- "data/raw/physical/Hutton_SoilErosionRisk_and_LUT/SoilErosionRisk/Erosion_risk_class_LUT.txt"
  erosion_lut <- read.delim(lut_path)
  erosion_risk_map <- c(
    "Li" = 0.1, "L2" = 0.2, "L3" = 0.3,
    "Lii" = 0.2, "Liii" = 0.3,
    "M1" = 0.4, "M2" = 0.5, "M3" = 0.6,
    "Mi" = 0.4, "Mii" = 0.5, "Miii" = 0.6, "Miv" = 0.7,
    "H" = 0.9, "H1" = 0.7, "H2" = 0.8, "H3" = 0.9
  )
  erosion_lut$risk_01 <- erosion_risk_map[erosion_lut$RiskClass]
  rcl_erosion <- as.matrix(erosion_lut[, c("DATAVALUE", "risk_01")])
  erosion_crop <- crop(erosion_raw, ext(hlcm_crop))
  erosion_risk <- terra::classify(erosion_crop, rcl_erosion)
  erosion_risk <- resample(erosion_risk, source_risk, method = "near")

  leach_sf <- st_read(leach_path, quiet = TRUE)
  leach_crop <- st_crop(leach_sf, st_bbox(cefas_crop)) %>% st_make_valid()
  leach_risk_map <- c("L" = 0.2, "I1" = 0.4, "I2" = 0.5,
                      "H3" = 0.7, "H2" = 0.8, "H1" = 0.9)
  leach_risk <- source_risk * 0 + NA
  ex_leach_L  <- terra::extract(leach_risk, vect(leach_crop[which(leach_crop$SLP == "L"), ]), cells = TRUE)
  leach_risk[ex_leach_L$cell] <- 0.2
  ex_leach_I1 <- terra::extract(leach_risk, vect(leach_crop[which(leach_crop$SLP == "I1"), ]), cells = TRUE)
  leach_risk[ex_leach_I1$cell] <- 0.4
  ex_leach_I2 <- terra::extract(leach_risk, vect(leach_crop[which(leach_crop$SLP == "I2"), ]), cells = TRUE)
  leach_risk[ex_leach_I2$cell] <- 0.5
  ex_leach_H3 <- terra::extract(leach_risk, vect(leach_crop[which(leach_crop$SLP == "H3"), ]), cells = TRUE)
  leach_risk[ex_leach_H3$cell] <- 0.7
  ex_leach_H2 <- terra::extract(leach_risk, vect(leach_crop[which(leach_crop$SLP == "H2"), ]), cells = TRUE)
  leach_risk[ex_leach_H2$cell] <- 0.8
  ex_leach_H1 <- terra::extract(leach_risk, vect(leach_crop[which(leach_crop$SLP == "H1"), ]), cells = TRUE)
  leach_risk[ex_leach_H1$cell] <- 0.9

  runoff_sf <- st_read(runoff_path, quiet = TRUE)
  runoff_crop <- st_crop(runoff_sf, st_bbox(cefas_crop)) %>% st_make_valid()
  runoff_risk <- source_risk * 0 + NA
  ex_run_low <- terra::extract(runoff_risk, vect(runoff_crop[which(runoff_crop$RUN_RISK == "Low"), ]), cells = TRUE)
  runoff_risk[ex_run_low$cell] <- 0.1
  ex_run_med <- terra::extract(runoff_risk, vect(runoff_crop[which(runoff_crop$RUN_RISK == "Medium"), ]), cells = TRUE)
  runoff_risk[ex_run_med$cell] <- 0.5
  ex_run_high <- terra::extract(runoff_risk, vect(runoff_crop[which(runoff_crop$RUN_RISK == "High"), ]), cells = TRUE)
  runoff_risk[ex_run_high$cell] <- 0.9

  twi_v <- values(twi_resampled)[, 1]
  eros_v <- values(erosion_risk)[, 1]
  leach_v <- values(leach_risk)[, 1]
  run_v <- values(runoff_risk)[, 1]

  ok_eros <- !is.na(twi_v) & !is.na(eros_v)
  ok_leach <- !is.na(twi_v) & !is.na(leach_v)
  ok_run <- !is.na(twi_v) & !is.na(run_v)

  cor_eros <- cor(twi_v[ok_eros], eros_v[ok_eros], method = "spearman")
  cor_leach <- cor(twi_v[ok_leach], leach_v[ok_leach], method = "spearman")
  cor_run <- cor(twi_v[ok_run], run_v[ok_run], method = "spearman")

  cat(sprintf("  Spearman correlation TWI vs Erosion:  %.4f  (n = %d overlap cells)\n",
              cor_eros, sum(ok_eros)))
  cat(sprintf("  Spearman correlation TWI vs Leaching: %.4f  (n = %d overlap cells)\n",
              cor_leach, sum(ok_leach)))
  cat(sprintf("  Spearman correlation TWI vs Runoff:   %.4f  (n = %d overlap cells)\n",
              cor_run, sum(ok_run)))

  hutton_validated <- TRUE
  hutton_cors <- c(erosion = cor_eros, leaching = cor_leach, runoff = cor_run)

} else {
  cat("  Hutton layers not found — skipping validation\n")
}


# ---- STEP 4: Weighted flow accumulation (risk routing) ----
# Accumulate diffuse_pollution_risk downstream through the drainage network.
# Each cell contributes its risk value (not 1) to all downstream cells.
# Done at GEBCO DEM resolution (~90K cells) then resampled.

cat("\n  --- Step 4: Weighted flow accumulation ---\n")

flow_dir_path <- "outputs/catchments/flow_dir.tif"
wbt_d8_pointer(dem = filled_path, output = flow_dir_path)
cat(sprintf("  D8 flow direction: %s\n", flow_dir_path))

flow_dir <- rast(flow_dir_path)

risk_at_dem <- resample(diffuse_pollution_risk, flow_dir, method = "average")
risk_dem_vals <- values(risk_at_dem)[, 1]
cat(sprintf("  Risk aggregated to DEM grid: %d x %d, non-NA = %d\n",
            ncol(risk_at_dem), nrow(risk_at_dem),
            sum(!is.na(risk_dem_vals))))
cat(sprintf("  Aggregated risk range: %.4f – %.4f (original: %.4f – %.4f)\n",
            min(risk_dem_vals, na.rm = TRUE), max(risk_dem_vals, na.rm = TRUE),
            min(values(diffuse_pollution_risk), na.rm = TRUE),
            max(values(diffuse_pollution_risk), na.rm = TRUE)))

dem_filled <- rast(filled_path)
fdir_vals <- values(flow_dir)[, 1]
risk_vals <- values(risk_at_dem)[, 1]
elev_vals <- values(dem_filled)[, 1]

cat("  D8 pointer unique values:\n")
print(sort(unique(na.omit(fdir_vals))))

cat(sprintf("  Risk at DEM grid before cleanup: min=%.6f, max=%.6f, NAs=%d, negatives=%d\n",
            min(risk_vals, na.rm = TRUE), max(risk_vals, na.rm = TRUE),
            sum(is.na(risk_vals)), sum(risk_vals < 0, na.rm = TRUE)))

risk_vals[is.na(risk_vals)] <- 0
risk_vals[risk_vals < 0] <- 0
accum_vals <- risk_vals

nr <- nrow(flow_dir)
nc <- ncol(flow_dir)

# D8 pointer encoding (WhiteboxTools): 1=E, 2=SE, 4=S, 8=SW, 16=W, 32=NW, 64=N, 128=NE
d8_dr <- c(`1`=0L, `2`=1L, `4`=1L, `8`=1L, `16`=0L, `32`=-1L, `64`=-1L, `128`=-1L)
d8_dc <- c(`1`=1L, `2`=1L, `4`=0L, `8`=-1L, `16`=-1L, `32`=-1L, `64`=0L, `128`=1L)

valid_cells <- which(!is.na(fdir_vals) & !is.na(elev_vals))
cell_order <- valid_cells[order(-elev_vals[valid_cells])]

cat(sprintf("  Routing risk through %d cells (highest to lowest)...\n", length(cell_order)))
t_start <- Sys.time()

for (ci in cell_order) {
  dir_val <- fdir_vals[ci]
  dir_key <- as.character(as.integer(dir_val))
  if (is.na(dir_val) || !(dir_key %in% names(d8_dr))) next

  r <- ((ci - 1) %/% nc) + 1L
  cc <- ((ci - 1) %% nc) + 1L
  r2 <- r + d8_dr[[dir_key]]
  c2 <- cc + d8_dc[[dir_key]]

  if (r2 >= 1L & r2 <= nr & c2 >= 1L & c2 <= nc) {
    target <- (r2 - 1L) * nc + c2
    accum_vals[target] <- accum_vals[target] + accum_vals[ci]
  }
}

t_elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
cat(sprintf("  Flow accumulation complete: %.1f seconds\n", t_elapsed))

risk_accum <- flow_dir
values(risk_accum) <- accum_vals

accum_valid <- accum_vals[valid_cells]
cat(sprintf("  Accumulated risk range: %.2f – %.2f\n",
            min(accum_valid, na.rm = TRUE), max(accum_valid, na.rm = TRUE)))
cat(sprintf("  Accumulated risk quantiles: 50%%=%.1f, 75%%=%.1f, 90%%=%.1f, 95%%=%.1f, 99%%=%.1f\n",
            quantile(accum_valid, 0.50, na.rm = TRUE),
            quantile(accum_valid, 0.75, na.rm = TRUE),
            quantile(accum_valid, 0.90, na.rm = TRUE),
            quantile(accum_valid, 0.95, na.rm = TRUE),
            quantile(accum_valid, 0.99, na.rm = TRUE)))
n_neg <- sum(accum_valid < 0)
cat(sprintf("  Cells with accum > 0: %d / %d (%.1f%%)\n",
            sum(accum_valid > 0), length(accum_valid),
            100 * sum(accum_valid > 0) / length(accum_valid)))
if (n_neg > 0) {
  cat(sprintf("  WARNING: %d cells have negative accumulated risk (min = %.2f)\n", n_neg, min(accum_valid)))
  cat("  This should not happen — check input risk raster for anomalies.\n")
}


risk_accum_path <- "outputs/catchments/risk_accumulated.tif"
writeRaster(risk_accum, risk_accum_path, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
cat(sprintf("  Saved: %s\n", risk_accum_path))


# ---- STEP 4b: River accumulated risk (SCIMAP raster extraction) ----
# Samples the SCIMAP risk_accum raster at both endpoints of each clipped
# river segment.  The downstream end has higher accumulated risk, so we
# take the max of both endpoints as the segment's value.

cat("\n  --- Step 4b: River accumulated risk (SCIMAP raster extraction) ---\n")

rivers_in_catch <- rivers_crop
n_seg <- nrow(rivers_in_catch)
cat(sprintf("  River segments: %d\n", n_seg))
rivers_in_catch$seg_length_m <- as.numeric(st_length(rivers_in_catch))

seg_first_xy <- t(sapply(seq_len(n_seg), function(i) {
  cc <- st_coordinates(rivers_in_catch$geometry[i])
  cc[1, 1:2]
}))
seg_last_xy <- t(sapply(seq_len(n_seg), function(i) {
  cc <- st_coordinates(rivers_in_catch$geometry[i])
  cc[nrow(cc), 1:2]
}))

pts_first <- vect(seg_first_xy, crs = "EPSG:27700")
pts_last  <- vect(seg_last_xy,  crs = "EPSG:27700")
risk_at_first <- terra::extract(risk_accum, pts_first)[, 2]
risk_at_last  <- terra::extract(risk_accum, pts_last)[, 2]
risk_at_first[is.na(risk_at_first)] <- 0
risk_at_last[is.na(risk_at_last)]   <- 0

rivers_in_catch$accum_risk_raw <- pmax(risk_at_first, risk_at_last)
riv_min <- min(rivers_in_catch$accum_risk_raw, na.rm = TRUE)
riv_max <- max(rivers_in_catch$accum_risk_raw, na.rm = TRUE)
rivers_in_catch$accum_risk <- if (riv_max > riv_min) {
  (rivers_in_catch$accum_risk_raw - riv_min) / (riv_max - riv_min)
} else {
  rep(0, n_seg)
}
cat(sprintf("  Accumulated risk (raw): %.4f – %.2f\n", riv_min, riv_max))
riv_accum <- rivers_in_catch$accum_risk
cat(sprintf("  River accum risk quantiles: 25%%=%.3f, 50%%=%.3f, 75%%=%.3f, 95%%=%.3f, max=%.3f\n",
            quantile(riv_accum, 0.25), quantile(riv_accum, 0.50),
            quantile(riv_accum, 0.75), quantile(riv_accum, 0.95),
            max(riv_accum)))

st_write(rivers_in_catch, "outputs/catchments/rivers_with_accumulated_risk.gpkg",
         delete_dsn = TRUE, quiet = TRUE)
cat("  Saved: outputs/catchments/rivers_with_accumulated_risk.gpkg\n")


# ---- STEP 5: Coastal pour points (geometry-based) ----
# For each named river: most seaward endpoint = pour point.
# Coast reference = SMR boundary (Mean High Water Springs), which wraps
# around estuary mouths so tidal rivers like the Urr/Nith are captured.

cat("\n----- Step 5: Coastal pour points (geometry-based) -----\n")

coast_thresh <- 100  # metres

# All segment endpoints
ep_first <- t(sapply(seq_len(n_seg), function(i) {
  cc <- st_coordinates(rivers_in_catch$geometry[i])
  cc[1, 1:2]
}))
ep_last <- t(sapply(seq_len(n_seg), function(i) {
  cc <- st_coordinates(rivers_in_catch$geometry[i])
  cc[nrow(cc), 1:2]
}))

all_eps <- rbind(
  data.frame(x = ep_first[, 1], y = ep_first[, 2], seg_idx = seq_len(n_seg)),
  data.frame(x = ep_last[, 1],  y = ep_last[, 2],  seg_idx = seq_len(n_seg))
)
all_eps$watercourse_name <- rivers_in_catch$watercourse_name[all_eps$seg_idx]
all_eps$accum_risk_raw   <- rivers_in_catch$accum_risk_raw[all_eps$seg_idx]
all_eps_sf <- st_as_sf(all_eps, coords = c("x", "y"), crs = 27700)

# Distance to SMR polygon (0 if inside the marine area, positive on land)
all_eps_sf$coast_dist <- as.numeric(st_distance(all_eps_sf, smr_solway_bng))
cat(sprintf("  Total endpoints: %d, within %dm of coast: %d\n",
            nrow(all_eps_sf), coast_thresh,
            sum(all_eps_sf$coast_dist <= coast_thresh)))

# Named rivers: most seaward endpoint per watercourse name
named_eps <- all_eps_sf[which(!is.na(all_eps_sf$watercourse_name) &
                              all_eps_sf$coast_dist <= coast_thresh), ]
pour_named <- named_eps %>%
  group_by(watercourse_name) %>%
  slice_min(coast_dist, n = 1, with_ties = FALSE) %>%
  ungroup()

# Unnamed segments: endpoints within threshold of coast
unnamed_eps <- all_eps_sf[which(is.na(all_eps_sf$watercourse_name) &
                                all_eps_sf$coast_dist <= coast_thresh), ]

pour_points <- rbind(pour_named, unnamed_eps)
in_aoi <- st_intersects(pour_points, scamp_aoi_sf, sparse = FALSE)[, 1]
pour_points <- pour_points[which(in_aoi), ]

# Join to catchments
pour_points <- st_join(pour_points, cefas_crop[, c("S_Catchmt", "S_CatNAME")],
                       join = st_nearest_feature)
pour_points <- pour_points[which(!is.na(pour_points$S_Catchmt)), ]

pp_min <- min(pour_points$accum_risk_raw, na.rm = TRUE)
pp_max <- max(pour_points$accum_risk_raw, na.rm = TRUE)
pour_points$accum_risk <- if (pp_max > pp_min) {
  (pour_points$accum_risk_raw - pp_min) / (pp_max - pp_min)
} else {
  rep(0, nrow(pour_points))
}

cat(sprintf("  Pour points: %d (%d named, %d unnamed)\n",
            nrow(pour_points), sum(!is.na(pour_points$watercourse_name)),
            sum(is.na(pour_points$watercourse_name))))
cat(sprintf("  Accum risk (raw): %.4f – %.2f\n", pp_min, pp_max))
cat("  Named pour points (top 15 by risk):\n")
pp_show <- st_drop_geometry(pour_points[order(-pour_points$accum_risk_raw), ])
print(head(pp_show[which(!is.na(pp_show$watercourse_name)),
                   c("watercourse_name", "S_CatNAME", "accum_risk")], 15))

# --- Mapview: check pour point locations before proceeding ---
smr_coast <- st_cast(st_boundary(smr_solway_bng), "MULTILINESTRING")
m_pp <- mapview(cefas_crop, col.regions = "grey90", alpha.regions = 0.3,
                layer.name = "Catchments") +
  mapview(smr_coast, color = "red", lwd = 2, layer.name = "SMR boundary (MHWS)") +
  mapview(scamp_coastline, color = "orange", lwd = 1.5,
          layer.name = "Current coastline (D&G admin)") +
  mapview(rivers_in_catch, color = "steelblue", lwd = 0.5,
          label = "watercourse_name", layer.name = "Rivers") +
  mapview(pour_points, zcol = "accum_risk", cex = 5,
          label = "watercourse_name", layer.name = "Pour points")
print(m_pp)

st_write(pour_points, "outputs/catchments/coastal_pour_points_raw.gpkg",
         delete_dsn = TRUE, quiet = TRUE)
cat("  Saved: outputs/catchments/coastal_pour_points_raw.gpkg\n")

# --- Spatial deduplication: cluster points within 500m, keep highest risk ---
cat("\n  --- Pour point deduplication (500m clusters) ---\n")
cat(sprintf("  Before dedup: %d pour points\n", nrow(pour_points)))

pp_dist <- st_distance(pour_points)
cluster_thresh <- 500  # metres
pp_clust <- cutree(hclust(as.dist(pp_dist), method = "single"), h = cluster_thresh)
pour_points$cluster <- pp_clust

pour_points <- pour_points %>%
  group_by(cluster) %>%
  slice_max(accum_risk_raw, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(-cluster)

cat(sprintf("  After dedup: %d pour points\n", nrow(pour_points)))

# --- Snap pour points to nearest river endpoint (geometry only, keep risk) ---
all_river_eps <- rbind(
  data.frame(x = ep_first[, 1], y = ep_first[, 2]),
  data.frame(x = ep_last[, 1],  y = ep_last[, 2])
)
all_river_eps_sf <- st_as_sf(all_river_eps, coords = c("x", "y"), crs = 27700)

nearest_idx <- st_nearest_feature(pour_points, all_river_eps_sf)
snap_dist <- as.numeric(st_distance(pour_points, all_river_eps_sf[nearest_idx, ],
                                     by_element = TRUE))
st_geometry(pour_points) <- st_geometry(all_river_eps_sf[nearest_idx, ])
cat(sprintf("  Snapped to nearest river endpoint (mean %.0fm, max %.0fm)\n",
            mean(snap_dist), max(snap_dist)))

m_pp2 <- mapview(cefas_crop, col.regions = "grey90", alpha.regions = 0.3,
                 layer.name = "Catchments") +
  mapview(smr_coast, color = "red", lwd = 2, layer.name = "SMR boundary (MHWS)") +
  mapview(rivers_in_catch, color = "steelblue", lwd = 0.5,
          label = "watercourse_name", layer.name = "Rivers") +
  mapview(pour_points, zcol = "accum_risk", cex = 5,
          label = "watercourse_name", layer.name = "Pour points (snapped)")
print(m_pp2)

st_write(pour_points, "outputs/catchments/coastal_pour_points.gpkg",
         delete_dsn = TRUE, quiet = TRUE)
cat("  Saved: outputs/catchments/coastal_pour_points.gpkg\n")


# ---- STEP 6: Per-catchment risk summaries ----

cat("\n  --- Step 6: Per-catchment risk summaries ---\n")

cell_area <- cellSize(diffuse_pollution_risk, unit = "ha")

ex_risk <- terra::extract(diffuse_pollution_risk, vect(cefas_crop), cells = TRUE)
ex_area <- terra::extract(cell_area, vect(cefas_crop), cells = TRUE)
names(ex_risk)[2] <- "risk"
names(ex_area)[2] <- "area_ha"
ex_risk$area_ha <- ex_area$area_ha
ex_risk$S_Catchmt <- cefas_crop$S_Catchmt[ex_risk$ID]

catch_risk <- ex_risk %>%
  group_by(S_Catchmt) %>%
  summarise(
    mean_diffuse_risk = round(weighted.mean(risk, area_ha, na.rm = TRUE), 4),
    total_ha = round(sum(area_ha, na.rm = TRUE), 1),
    .groups = "drop"
  )

ex_source <- terra::extract(source_risk, vect(cefas_crop), cells = TRUE)
ex_src_area <- terra::extract(cell_area, vect(cefas_crop), cells = TRUE)
names(ex_source)[2] <- "src_risk"
names(ex_src_area)[2] <- "area_ha"
ex_source$area_ha <- ex_src_area$area_ha
ex_source$S_Catchmt <- cefas_crop$S_Catchmt[ex_source$ID]

catch_source <- ex_source %>%
  group_by(S_Catchmt) %>%
  summarise(
    mean_source_risk = round(weighted.mean(src_risk, area_ha, na.rm = TRUE), 4),
    .groups = "drop"
  )

pour_max <- pour_points %>%
  st_drop_geometry() %>%
  group_by(S_Catchmt) %>%
  summarise(max_pour_accum = max(accum_risk, na.rm = TRUE), .groups = "drop")
pour_max$max_pour_accum[is.infinite(pour_max$max_pour_accum)] <- NA

catch_risk <- catch_risk %>%
  left_join(catch_source, by = "S_Catchmt") %>%
  left_join(pour_max, by = "S_Catchmt") %>%
  left_join(st_drop_geometry(cefas_crop[, c("S_Catchmt", "S_CatNAME")]),
            by = "S_Catchmt") %>%
  arrange(desc(mean_diffuse_risk))

cat("  Catchment risk summary:\n")
print(catch_risk, n = 50)

write.csv(catch_risk, "outputs/catchments/catchment_risk_summary.csv", row.names = FALSE)
cat("  Saved: outputs/catchments/catchment_risk_summary.csv\n")

cefas_risk <- merge(cefas_crop, catch_risk, by = c("S_Catchmt", "S_CatNAME"), all.x = TRUE)
st_write(cefas_risk, "outputs/catchments/catchments_with_risk.gpkg",
         delete_dsn = TRUE, quiet = TRUE)
cat("  Saved: outputs/catchments/catchments_with_risk.gpkg\n")


# ---- Catchment composition (for downstream maps) ----

catch_composition <- lc_props %>%
  group_by(S_Catchmt, S_CatNAME) %>%
  summarise(
    n_classes = n_distinct(eunis_code),
    total_pixels = sum(n_pixels),
    dominant_class = eunis_code[which.max(proportion)],
    dominant_desc  = eunis_desc[which.max(proportion)],
    dominant_pct   = round(100 * max(proportion), 1),
    .groups = "drop"
  ) %>%
  arrange(S_CatNAME)

print(catch_composition)


# ---- STEP 7: Diffuse pollution risk maps ----

cat("\n  --- Step 7: Diffuse pollution risk maps ---\n")

dir.create("figures/catchments", showWarnings = FALSE, recursive = TRUE)

# Panel 1: Source risk per catchment
p1 <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = cefas_risk, aes(fill = mean_source_risk),
          colour = "grey30", linewidth = 0.3) +
  scale_fill_gradient(name = "Mean\nsource\nrisk",
                      low = "#fee0d2", high = "#67000d", na.value = "grey80") +
  geom_sf(data = rivers_in_catch, colour = "steelblue", linewidth = 0.2, alpha = 0.4) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.6) +
  labs(title = "Source Risk (SCIMAP land use weighting)") +
  map_theme + coord_sf(crs = 27700)

# Panel 2: Diffuse pollution risk per catchment
p2 <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = cefas_risk, aes(fill = mean_diffuse_risk),
          colour = "grey30", linewidth = 0.3) +
  scale_fill_gradient(name = "Mean\ndiffuse\nrisk",
                      low = "#fee0d2", high = "#67000d", na.value = "grey80") +
  geom_sf(data = rivers_in_catch, colour = "steelblue", linewidth = 0.2, alpha = 0.4) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.6) +
  labs(title = "Diffuse Pollution Risk (source × TWI)") +
  map_theme + coord_sf(crs = 27700)

p_rivlen_hist <- ggplot(rivers_in_catch, aes(x = seg_length_m)) +
  geom_histogram(bins = 50, fill = "steelblue", colour = "white") +
  scale_x_log10() +
  labs(title = "River Segment Length Distribution",
       x = "Segment length (m, log scale)", y = "Count") +
  theme_minimal()
ggsave("figures/threats/river_segment_length_hist.png", p_rivlen_hist,
       width = 8, height = 5, units = "in", dpi = 300)
cat("  Saved: figures/threats/river_segment_length_hist.png\n")
cat(sprintf("  Segment length quantiles (m): 5%%=%.0f, 25%%=%.0f, 50%%=%.0f, 75%%=%.0f, 95%%=%.0f\n",
            quantile(rivers_in_catch$seg_length_m, 0.05),
            quantile(rivers_in_catch$seg_length_m, 0.25),
            quantile(rivers_in_catch$seg_length_m, 0.50),
            quantile(rivers_in_catch$seg_length_m, 0.75),
            quantile(rivers_in_catch$seg_length_m, 0.95)))

# Panel 3: River network coloured by accumulated risk, top 10 labelled
rivers_ranked <- rivers_in_catch[order(-rivers_in_catch$accum_risk, na.last = TRUE), ]
rivers_top10 <- head(rivers_ranked[!is.na(rivers_ranked$accum_risk), ], 10)
top10_labels <- st_centroid(rivers_top10)

river_name_col <- "watercourse_name"

rivers_plot <- rivers_in_catch
rivers_plot$accum_risk[rivers_plot$accum_risk == 0] <- NA

p3 <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = cefas_risk, fill = "white", colour = "grey70", linewidth = 0.2) +
  geom_sf(data = rivers_plot[is.na(rivers_plot$accum_risk), ],
          colour = "steelblue", linewidth = 0.3, alpha = 0.4) +
  geom_sf(data = rivers_plot[which(!is.na(rivers_plot$accum_risk)), ],
          aes(colour = accum_risk), linewidth = 0.6) +
  scale_colour_viridis_c(name = "Accumulated\nrisk", option = "inferno",
                         direction = -1, trans = "sqrt") +
  geom_sf_label(data = top10_labels, aes(label = .data[[river_name_col]]),
                size = 2.2, fill = alpha("white", 0.7), linewidth = 0.2) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.6) +
  labs(title = "River Network — Accumulated Risk (top 10 labelled)") +
  map_theme + coord_sf(crs = 27700)

p_pour_hist <- ggplot(pour_points[!is.na(pour_points$accum_risk), ],
                      aes(x = accum_risk)) +
  geom_histogram(bins = 30, fill = "steelblue", colour = "white") +
  labs(title = "Pour Point Accumulated Risk Distribution",
       x = "Accumulated risk", y = "Count") +
  theme_minimal()
ggsave("figures/threats/pour_point_histogram.png", p_pour_hist,
       width = 8, height = 5, units = "in", dpi = 300)
cat("  Saved: figures/threats/pour_point_histogram.png\n")

# Panel 4: Catchments coloured by mean risk, pour points sized by accumulated risk
p4 <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = cefas_risk, aes(fill = mean_diffuse_risk),
          colour = "grey30", linewidth = 0.3) +
  scale_fill_gradient(name = "Mean\ndiffuse\nrisk",
                      low = "#fee0d2", high = "#67000d", na.value = "grey80") +
  geom_sf(data = pour_points, aes(size = accum_risk, colour = accum_risk),
          alpha = 0.7, shape = 16) +
  scale_size_continuous(name = "Pour point\naccumulated\nrisk", range = c(2, 8), trans = "sqrt") +
  scale_colour_viridis_c(name = "Pour point\naccumulated\nrisk", option = "inferno",
                         direction = -1, trans = "sqrt", guide = "none") +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.6) +
  labs(title = "Catchment Risk with Pour Points") +
  map_theme + coord_sf(crs = 27700)

# Panel 5: Pixel-wise accumulated risk (SCIMAP flow accumulation raster)
risk_accum_crop <- mask(crop(risk_accum, ext(vect(cefas_crop))), vect(cefas_crop))
agg_factor <- max(1, round(ncol(risk_accum_crop) / 500))
risk_accum_plot <- terra::aggregate(risk_accum_crop, fact = agg_factor, fun = "max", na.rm = TRUE)
cat(sprintf("  Raster aggregated %dx for plotting: %d x %d -> %d x %d\n",
            agg_factor, nrow(risk_accum_crop), ncol(risk_accum_crop),
            nrow(risk_accum_plot), ncol(risk_accum_plot)))
p5_accum <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  tidyterra::geom_spatraster(data = risk_accum_plot) +
  scale_fill_viridis_c(name = "Accumulated\nrisk", option = "inferno",
                       direction = -1, na.value = "transparent", trans = "sqrt") +
  geom_sf(data = rivers_in_catch, colour = "steelblue", linewidth = 0.2, alpha = 0.3) +
  geom_sf(data = pour_points, colour = "red", size = 1.5) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.6) +
  labs(title = "SCIMAP Accumulated Risk (pixel-wise)",
       subtitle = "DEM weighted flow accumulation — sqrt scale") +
  map_theme + coord_sf(crs = 27700)
ggsave("figures/threats/risk_accum_raster.png", p5_accum,
       width = 12, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/threats/risk_accum_raster.png\n")

m <- mapview(cefas_risk, zcol = "mean_diffuse_risk", layer.name = "Catchment diffuse risk") +
  mapview(rivers_in_catch, zcol = "accum_risk",
          label = "watercourse_name", layer.name = "River accumulated risk") +
  mapview(rivers_top10, zcol = "accum_risk", color = "red", lwd = 3,
          label = "watercourse_name", layer.name = "Top 10 rivers") +
  mapview(pour_points, zcol = "accum_risk", cex = "accum_risk",
          label = "watercourse_name", layer.name = "Pour points (network outlets)")
print(m)

# Panel 5 (optional): TWI vs Hutton validation
if (hutton_validated) {
  cor_subtitle <- sprintf(
    "Spearman r: erosion=%.3f, leaching=%.3f, runoff=%.3f",
    hutton_cors["erosion"], hutton_cors["leaching"], hutton_cors["runoff"]
  )

  twi_eros_df <- data.frame(
    twi = twi_v[ok_eros],
    hutton = eros_v[ok_eros]
  )
  set.seed(42)
  samp_idx <- sample(nrow(twi_eros_df), min(50000, nrow(twi_eros_df)))

  p5 <- ggplot(twi_eros_df[samp_idx, ], aes(x = twi, y = hutton)) +
    geom_hex(bins = 60) +
    scale_fill_viridis_c(name = "Count", option = "inferno", trans = "log10") +
    geom_smooth(method = "lm", colour = "red", linewidth = 0.8, se = FALSE) +
    labs(x = "TWI (normalised)", y = "Hutton erosion risk (0-1)",
         title = "TWI vs Hutton Erosion Risk",
         subtitle = cor_subtitle) +
    theme_minimal()

  p_combined <- (p1 + p2) / (p3 + p4) / (p5 + plot_spacer()) +
    plot_annotation(
      title = "Diffuse Pollution Risk per Catchment",
      subtitle = "SCIMAP source risk × TWI hydrological connectivity — CEFAS/IHU catchments",
      theme = theme(plot.title = element_text(size = 16, face = "bold"),
                    plot.subtitle = element_text(size = 12))
    )
  ggsave("figures/catchments/diffuse_pollution_risk_comparison.png", p_combined,
         width = 16, height = 18, units = "in", dpi = 300)

} else {
  p_combined <- (p1 + p2) / (p3 + p4) +
    plot_annotation(
      title = "Diffuse Pollution Risk per Catchment",
      subtitle = "SCIMAP source risk × TWI hydrological connectivity — CEFAS/IHU catchments",
      theme = theme(plot.title = element_text(size = 16, face = "bold"),
                    plot.subtitle = element_text(size = 12))
    )
  ggsave("figures/catchments/diffuse_pollution_risk_comparison.png", p_combined,
         width = 16, height = 12, units = "in", dpi = 300)
}

cat("  Saved: figures/catchments/diffuse_pollution_risk_comparison.png\n")


# --- 4f. Catchment maps ---

cat("\n----- 4f. Catchment land cover maps -----\n")

cefas_comp <- merge(cefas_crop, catch_composition,
                     by = c("S_Catchmt", "S_CatNAME"), all.x = TRUE)

cefas_extent <- st_bbox(cefas_crop)

p_catch_dominant <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = cefas_comp, aes(fill = dominant_desc),
          colour = "grey30", linewidth = 0.3) +
  scale_fill_brewer(name = "Dominant\nland cover", palette = "Set3") +
  geom_sf(data = rivers_in_catch, colour = "steelblue", linewidth = 0.2, alpha = 0.4) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.6) +
  geom_sf(data = smr_solway_bng, fill = NA, colour = "navy",
          linewidth = 0.4, linetype = "dashed") +
  labs(title = "Dominant Land Cover per Catchment",
       subtitle = "HLCM 2022 EUNIS Level 2 — CEFAS/IHU catchments draining into SCAMP AOI") +
  map_theme + coord_sf(crs = 27700)
ggsave("figures/threats/catchment_dominant_landcover.png", p_catch_dominant,
       width = 14, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/threats/catchment_dominant_landcover.png\n")

p_catch_pour <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = cefas_risk, aes(fill = mean_diffuse_risk),
          colour = "grey30", linewidth = 0.3) +
  scale_fill_gradient(name = "Mean\ndiffuse\nrisk",
                      low = "#fee0d2", high = "#67000d", na.value = "grey80") +
  geom_sf(data = rivers_in_catch, colour = "steelblue", linewidth = 0.3, alpha = 0.5) +
  geom_sf(data = pour_points, aes(size = accum_risk),
          colour = "red", alpha = 0.7, shape = 17) +
  scale_size_continuous(name = "Pour point\naccumulated risk", range = c(1.5, 5)) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.6) +
  labs(title = "Catchment Pour Points — River Discharge to Coast",
       subtitle = "Triangles sized by accumulated diffuse pollution risk at pour point") +
  map_theme + coord_sf(crs = 27700)
ggsave("figures/threats/catchment_pour_points.png", p_catch_pour,
       width = 14, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/threats/catchment_pour_points.png\n")

plantation_pct <- lc_props %>%
  filter(eunis_code == "G3.F") %>%
  select(S_Catchmt, plantation_pct = proportion) %>%
  mutate(plantation_pct = round(plantation_pct * 100, 1))

cefas_plantation <- merge(cefas_crop, plantation_pct, by = "S_Catchmt", all.x = TRUE)
cefas_plantation$plantation_pct[is.na(cefas_plantation$plantation_pct)] <- 0

p_catch_plantation <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = cefas_plantation, aes(fill = plantation_pct),
          colour = "grey30", linewidth = 0.3) +
  scale_fill_viridis_c(name = "Conifer\nplantation\n(% area)",
                       na.value = "grey80", option = "inferno") +
  geom_sf(data = rivers_in_catch, colour = "steelblue", linewidth = 0.2, alpha = 0.4) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.6) +
  labs(title = "Conifer Plantation Cover per Catchment (G3.F)",
       subtitle = "HLCM 2022 — Highly artificial coniferous plantations") +
  map_theme + coord_sf(crs = 27700)
ggsave("figures/threats/catchment_plantation_pct.png", p_catch_plantation,
       width = 14, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/threats/catchment_plantation_pct.png\n")


# --- 4g. Save outputs ---

cat("\n----- 4g. Save catchment outputs -----\n")

st_write(cefas_comp, "outputs/threats/catchment_landcover.gpkg",
         delete_dsn = TRUE, quiet = TRUE)
cat("  Saved: outputs/threats/catchment_landcover.gpkg\n")

write.csv(as.data.frame(lc_props %>% select(-geometry)),
          "outputs/threats/catchment_landcover_proportions.csv",
          row.names = FALSE)
cat("  Saved: outputs/threats/catchment_landcover_proportions.csv\n")

write.csv(as.data.frame(catch_composition),
          "outputs/threats/catchment_composition_summary.csv",
          row.names = FALSE)
cat("  Saved: outputs/threats/catchment_composition_summary.csv\n")

st_write(pour_points, "outputs/threats/pour_points.gpkg",
         delete_dsn = TRUE, quiet = TRUE)
cat("  Saved: outputs/threats/pour_points.gpkg\n")

dim(risk_accum_crop)
# ---- STEP 8: Marine pollution exposure (sea-path decay from pour points) ----
# BFS computes shortest sea-path distance from each pour point (routes around land)
# Exponential decay applied: exposure = risk × exp(-sea_dist / lambda)
# Concept: Halpern et al. (2008, 2015) cost-path routing
# Lambda: 1852 m (1 nautical mile — WFD coastal water body definition)

cat("\n----- Step 8: Marine pollution exposure (sea-path distance decay) -----\n")

gebco_marine <- rast(gebco_path)
aoi_wgs84 <- st_transform(scamp_aoi_sf, crs(gebco_marine))
gebco_aoi <- crop(gebco_marine, ext(vect(aoi_wgs84)))
gebco_aoi_bng <- project(gebco_aoi, "EPSG:27700", method = "bilinear")
marine_raster <- mask(gebco_aoi_bng, vect(scamp_aoi_sf))
marine_raster[marine_raster > 0] <- NA

n_marine <- sum(!is.na(values(marine_raster)))
cell_res_m <- res(marine_raster)[1]
cat(sprintf("  Marine raster: %d x %d, res %.0fm, marine cells: %d\n",
            ncol(marine_raster), nrow(marine_raster), cell_res_m, n_marine))

lambda <- 5000  # 1 nautical mile (WFD coastal water body limit)
max_dist <- lambda * 7  # exp(-7) ≈ 0.001, negligible beyond this (~13 km)

pp_valid <- pour_points[which(!is.na(pour_points$accum_risk)), ]
cat(sprintf("  Lambda: %d m (1 nm, WFD), max search: %.0f m\n", lambda, max_dist))
cat(sprintf("  Computing sea-path decay from %d pour points...\n", nrow(pp_valid)))

sea_mask <- !is.na(values(marine_raster)[, 1])
n_cells <- ncell(marine_raster)
nc <- ncol(marine_raster)
exposure_vals <- rep(0, n_cells)

adj_offsets <- c(-nc - 1L, -nc, -nc + 1L, -1L, 1L, nc - 1L, nc, nc + 1L)
adj_dists <- c(cell_res_m * sqrt(2), cell_res_m, cell_res_m * sqrt(2),
               cell_res_m, cell_res_m,
               cell_res_m * sqrt(2), cell_res_m, cell_res_m * sqrt(2))

max_queue <- min(n_marine * 3L, 500000L)
queue_cell <- integer(max_queue)
queue_dist <- numeric(max_queue)

t_start <- Sys.time()

for (i in seq_len(nrow(pp_valid))) {
  risk_i <- pp_valid$accum_risk[i]
  if (is.na(risk_i) || risk_i <= 0) next

  pt_coords <- st_coordinates(pp_valid[i, ])
  start_cell <- cellFromXY(marine_raster, matrix(pt_coords[1, 1:2], ncol = 2))

  if (is.na(start_cell) || !sea_mask[start_cell]) {
    near_cells <- which(sea_mask)
    if (length(near_cells) == 0) next
    xy_sea <- xyFromCell(marine_raster, near_cells)
    dists_sq <- (xy_sea[, 1] - pt_coords[1, 1])^2 + (xy_sea[, 2] - pt_coords[1, 2])^2
    start_cell <- near_cells[which.min(dists_sq)]
  }

  cell_dist <- rep(Inf, n_cells)
  cell_dist[start_cell] <- 0
  head_idx <- 1L
  tail_idx <- 1L
  queue_cell[1] <- start_cell
  queue_dist[1] <- 0

  while (head_idx <= tail_idx) {
    cc <- queue_cell[head_idx]
    cd <- queue_dist[head_idx]
    head_idx <- head_idx + 1L

    if (cd > cell_dist[cc]) next
    if (cd > max_dist) next

    cur_col <- ((cc - 1L) %% nc) + 1L

    for (j in 1:8) {
      nb <- cc + adj_offsets[j]
      if (nb < 1L || nb > n_cells) next
      if (!sea_mask[nb]) next
      nb_col <- ((nb - 1L) %% nc) + 1L
      if (abs(nb_col - cur_col) > 1L) next

      new_dist <- cd + adj_dists[j]
      if (new_dist < cell_dist[nb]) {
        cell_dist[nb] <- new_dist
        tail_idx <- tail_idx + 1L
        if (tail_idx > length(queue_cell)) {
          queue_cell <- c(queue_cell, integer(max_queue))
          queue_dist <- c(queue_dist, numeric(max_queue))
        }
        queue_cell[tail_idx] <- nb
        queue_dist[tail_idx] <- new_dist
      }
    }
  }

  reached <- which(is.finite(cell_dist) & sea_mask)
  exposure_vals[reached] <- exposure_vals[reached] +
    risk_i * exp(-cell_dist[reached] / lambda)

  if (i %% 10 == 0 || i == nrow(pp_valid)) {
    cat(sprintf("\r  Pour point %d / %d (%d cells reached)", i, nrow(pp_valid), length(reached)))
  }
}
cat("\n")

t_elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
cat(sprintf("  Sea-path decay computed: %.1f seconds\n", t_elapsed))

exposure_rast <- rast(marine_raster)
v <- rep(NA_real_, n_cells)
v[sea_mask] <- exposure_vals[sea_mask]
values(exposure_rast) <- v

exp_vals <- exposure_vals[sea_mask]
exp_min <- min(exp_vals, na.rm = TRUE)
exp_max <- max(exp_vals, na.rm = TRUE)
exposure_norm <- (exposure_rast - exp_min) / (exp_max - exp_min)

cat(sprintf("  Exposure range (raw): %.6f – %.6f\n", exp_min, exp_max))
cat(sprintf("  Exposure quantiles: 50%%=%.4f, 75%%=%.4f, 90%%=%.4f, 95%%=%.4f\n",
            quantile(exp_vals, 0.50, na.rm = TRUE),
            quantile(exp_vals, 0.75, na.rm = TRUE),
            quantile(exp_vals, 0.90, na.rm = TRUE),
            quantile(exp_vals, 0.95, na.rm = TRUE)))

writeRaster(exposure_norm, "outputs/catchments/marine_pollution_exposure.tif",
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))
cat("  Saved: outputs/catchments/marine_pollution_exposure.tif\n")

p_exposure <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = cefas_risk, fill = "grey95", colour = "grey70", linewidth = 0.2) +
  tidyterra::geom_spatraster(data = exposure_norm) +
  scale_fill_viridis_c(name = "Pollution\nexposure", option = "inferno",
                       direction = -1, na.value = "transparent") +
  geom_sf(data = pour_points, colour = "red", size = 1.5, alpha = 0.5) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.6) +
  labs(title = "Marine Diffuse Pollution Exposure (sea-path distance decay)",
       subtitle = sprintf("Lambda = %d m (1 nm WFD), %d pour points", lambda, nrow(pp_valid))) +
  map_theme + coord_sf(crs = 27700)
ggsave("figures/threats/marine_pollution_exposure.png", p_exposure,
       width = 12, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/threats/marine_pollution_exposure.png\n")


urr_idx <- which(grepl("Urr", rivers_in_catch$watercourse_name, ignore.case = TRUE))
cat(sprintf("\n  Urr Water: %d segments found, accum_risk range: %.4f – %.4f\n",
            length(urr_idx),
            min(rivers_in_catch$accum_risk[urr_idx], na.rm = TRUE),
            max(rivers_in_catch$accum_risk[urr_idx], na.rm = TRUE)))
print(st_drop_geometry(rivers_in_catch[urr_idx, c("watercourse_name", "accum_risk")]))

urr_rivers <- rivers_in_catch[urr_idx, ]
urr_bbox <- st_bbox(urr_rivers)
p_urr <- ggplot() +
  geom_sf(data = cefas_risk, fill = "grey95", colour = "grey70", linewidth = 0.2) +
  geom_sf(data = rivers_in_catch, colour = "steelblue", linewidth = 0.3, alpha = 0.4) +
  geom_sf(data = urr_rivers, aes(colour = accum_risk), linewidth = 1.2) +
  scale_colour_viridis_c(name = "Accumulated\nrisk", option = "inferno", direction = -1) +
  coord_sf(xlim = c(urr_bbox["xmin"], urr_bbox["xmax"]),
           ylim = c(urr_bbox["ymin"], urr_bbox["ymax"]),
           crs = 27700) +
  labs(title = "Urr Water — accumulated diffuse pollution risk") +
  map_theme
ggsave("figures/threats/urr_water_risk.png", p_urr,
       width = 10, height = 8, units = "in", dpi = 300)
cat("  Saved: figures/threats/urr_water_risk.png\n")



# ---- STEP 9: Integrate marine pollution into cumulative threat exposure ----

cat("\n----- Step 9: Integrate marine pollution into cumulative threats -----\n")

# Resample marine pollution exposure to the 1km reference grid
pollution_1km <- resample(exposure_norm, ref_rast, method = "bilinear")
pollution_1km <- mask(pollution_1km, smr_mask)
cat(sprintf("  Marine pollution resampled to 1km: %d x %d\n",
            nrow(pollution_1km), ncol(pollution_1km)))

# Add to normalised layers and rebuild cumulative pressure (excluding noise)
norm_rasters[["diffuse_pollution"]] <- pollution_1km
raw_rasters[["diffuse_pollution"]]  <- resample(exposure_rast, ref_rast, method = "bilinear") %>%
  mask(smr_mask)

noise_layers <- grep("^noise_", names(norm_rasters), value = TRUE)
stack_layers <- setdiff(names(norm_rasters), noise_layers)
cat(sprintf("  Excluding %d noise layers: %s\n", length(noise_layers),
            paste(noise_layers, collapse = ", ")))

pressure_stack_full <- rast(norm_rasters[stack_layers])
cumulative_pressure_full <- sum(pressure_stack_full, na.rm = TRUE)
cumulative_pressure_full <- mask(cumulative_pressure_full, smr_mask)

# Fill coastal edge gaps: iteratively extend nearest values until no NA cells remain
smr_valid <- !is.na(values(smr_mask))
n_na <- sum(is.na(values(cumulative_pressure_full)) & smr_valid)
pass <- 0L
while (n_na > 0) {
  pass <- pass + 1L
  cumulative_pressure_full <- focal(cumulative_pressure_full, w = 3, fun = "mean",
                                     na.policy = "only", na.rm = TRUE)
  cumulative_pressure_full <- mask(cumulative_pressure_full, smr_mask)
  n_na_new <- sum(is.na(values(cumulative_pressure_full)) & smr_valid)
  if (n_na_new == n_na) break
  n_na <- n_na_new
}
cat(sprintf("  Coastal gap fill: %d passes, %d NA cells remaining\n", pass, n_na))

cat(sprintf("  Layers in updated stack: %d (%s)\n",
            nlyr(pressure_stack_full), paste(names(norm_rasters), collapse = ", ")))
cat(sprintf("  Cumulative pressure range (with pollution): %.2f – %.2f\n",
            minmax(cumulative_pressure_full)[1], minmax(cumulative_pressure_full)[2]))

p_cumulative_full <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  tidyterra::geom_spatraster(data = cumulative_pressure_full) +
  scale_fill_viridis_c(name = "Cumulative\npressure", na.value = "transparent", option = "magma") +
  {if (nrow(vec_polys) > 0) geom_sf(data = vec_polys, aes(colour = threat_type),
                                     fill = NA, linewidth = 0.5)} +
  {if (nrow(vec_pts) > 0) geom_sf(data = vec_pts, aes(colour = threat_type),
                                   size = 1.5, alpha = 0.8)} +
  scale_colour_brewer(name = "Vector threats", palette = "Set1") +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  labs(title = "Cumulative Marine Pressure Surface (including diffuse pollution)",
       subtitle = sprintf("%d raster layers (normalised 0–1) + vector threats",
                          nlyr(pressure_stack_full))) +
  map_theme + coord_sf(crs = 27700)
ggsave("figures/threats/cumulative_pressure_with_pollution.png", p_cumulative_full,
       width = 14, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/threats/cumulative_pressure_with_pollution.png\n")

writeRaster(cumulative_pressure_full, "outputs/threats/cumulative_pressure_with_pollution_1km.tif",
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))
cat("  Saved: outputs/threats/cumulative_pressure_with_pollution_1km.tif\n")


# ---- STEP 10: Individual threat figures ----

cat("\n----- Step 10: Individual threat figures -----\n")

dir.create("figures/threats/individual", showWarnings = FALSE, recursive = TRUE)

# Raster threat layers
for (lyr_name in names(norm_rasters)) {
  p <- ggplot() +
    geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
    tidyterra::geom_spatraster(data = norm_rasters[[lyr_name]]) +
    scale_fill_viridis_c(name = "Normalised\npressure (0–1)",
                         na.value = "transparent", option = "inferno") +
    geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
    labs(title = gsub("_", " ", tools::toTitleCase(lyr_name))) +
    map_theme + coord_sf(crs = 27700)
  ggsave(sprintf("figures/threats/individual/%s.png", lyr_name), p,
         width = 12, height = 10, units = "in", dpi = 300)
  cat(sprintf("  Saved: figures/threats/individual/%s.png\n", lyr_name))
}

# Cumulative pressure
ggsave("figures/threats/individual/cumulative_pressure.png", p_cumulative_full,
       width = 12, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/threats/individual/cumulative_pressure.png\n")

# River accumulated risk
p_rivers_risk <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = cefas_risk, fill = "grey95", colour = "grey70", linewidth = 0.2) +
  geom_sf(data = rivers_in_catch[which(!is.na(rivers_in_catch$accum_risk) &
                                        rivers_in_catch$accum_risk > 0), ],
          aes(colour = accum_risk), linewidth = 0.6) +
  scale_colour_viridis_c(name = "Accumulated\nrisk", option = "inferno",
                         direction = -1, trans = "sqrt") +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.6) +
  labs(title = "River Accumulated Diffuse Pollution Risk (SCIMAP)") +
  map_theme + coord_sf(crs = 27700)
ggsave("figures/threats/individual/river_accumulated_risk.png", p_rivers_risk,
       width = 12, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/threats/individual/river_accumulated_risk.png\n")

# Pour points
p_pour <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = cefas_risk, aes(fill = mean_diffuse_risk),
          colour = "grey30", linewidth = 0.3) +
  scale_fill_gradient(name = "Catchment\nmean risk",
                      low = "#fee0d2", high = "#67000d", na.value = "grey80") +
  geom_sf(data = pour_points, aes(size = accum_risk, colour = accum_risk),
          alpha = 0.7, shape = 16) +
  scale_size_continuous(name = "Pour point\nrisk", range = c(2, 8), trans = "sqrt") +
  scale_colour_viridis_c(name = "Pour point\nrisk", option = "inferno",
                         direction = -1, trans = "sqrt", guide = "none") +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.6) +
  labs(title = "Coastal Pour Points — Accumulated Diffuse Pollution Risk") +
  map_theme + coord_sf(crs = 27700)
ggsave("figures/threats/individual/pour_points.png", p_pour,
       width = 12, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/threats/individual/pour_points.png\n")

cat(sprintf("\n  Total individual figures saved: %d\n",
            length(names(norm_rasters)) + 3))


cat("\n\n===== 05_localthreats.R COMPLETE =====\n")
cat("  TODO: INNS pressure layer (non-native species beyond macroalgae)\n")
