##=====================================================================
## SCAMP GIS — Physical, Oceanographic & Climate Data
## 04_physical.R
##=====================================================================
## Loads CMEMS reanalysis layers, crops to SCAMP AOI, and produces
## a Solway-specific seagrass potential layer using EA methodology.
##=====================================================================
install.packages('whitebox')
require(sf)
require(terra)
require(dplyr)
require(tidyr)
require(readr)
require(ggplot2)
require(ggspatial)
require(ncdf4)
require(jsonlite)

##=====================================================================
## 0. SETUP — Boundaries, AOI, directories, CMEMS helper
##=====================================================================

cat("\n=== 0. SETUP ===\n")

## ---- Load boundaries ----

aoi_path <- "data/processed/boundaries/scamp_aoi.gpkg"
if (!file.exists(aoi_path)) stop("SCAMP AOI not found: ", aoi_path)
scamp_aoi_sf <- st_read(aoi_path, quiet = TRUE)
if (st_crs(scamp_aoi_sf)$epsg != 27700) scamp_aoi_sf <- st_transform(scamp_aoi_sf, 27700)

scamp_wgs <- st_transform(scamp_aoi_sf, 4326)
bb <- st_bbox(scamp_wgs)
cat(sprintf("  AOI bounding box (WGS84): %.4f to %.4f lon, %.4f to %.4f lat\n",
            bb["xmin"], bb["xmax"], bb["ymin"], bb["ymax"]))

smr_path <- "data/raw/boundaries/administrative_units_scottish_marine_regions/administrative_units_scottish_marine_regionsPolygon.shp"
if (file.exists(smr_path)) {
  smr_solway_bng <- st_read(smr_path, quiet = TRUE) %>%
    filter(grepl("Solway", objnam, ignore.case = TRUE)) %>%
    st_transform(27700)
  cat(sprintf("  SMR Solway loaded: %d features\n", nrow(smr_solway_bng)))
} else {
  warning("SMR shapefile not found: ", smr_path)
  smr_solway_bng <- st_sf(geometry = st_sfc(crs = 27700))
}

las_path <- "data/raw/boundaries/Local_Authority_Boundaries_-_Scotland/pub_las.shp"
if (file.exists(las_path)) {
  las_dg_bng <- st_read(las_path, quiet = TRUE) %>%
    filter(grepl("dumfries", local_auth, ignore.case = TRUE)) %>%
    st_transform(27700)
  cat(sprintf("  D&G boundary loaded: %d features\n", nrow(las_dg_bng)))
} else {
  warning("Local authority shapefile not found: ", las_path)
  las_dg_bng <- st_sf(geometry = st_sfc(crs = 27700))
}

coast_path <- "data/processed/boundaries/scamp_coastline.gpkg"
if (file.exists(coast_path)) {
  scamp_coastline <- st_read(coast_path, quiet = TRUE)
  if (st_crs(scamp_coastline)$epsg != 27700) scamp_coastline <- st_transform(scamp_coastline, 27700)
  cat(sprintf("  Coastline loaded: %d features\n", nrow(scamp_coastline)))
} else {
  warning("Coastline not found: ", coast_path)
  scamp_coastline <- st_sf(geometry = st_sfc(crs = 27700))
}

## ---- Output directories ----

dir.create("data/processed/physical",  showWarnings = FALSE, recursive = TRUE)
dir.create("data/threats_global/cmems", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/physical",          showWarnings = FALSE, recursive = TRUE)
dir.create("figures/physical",          showWarnings = FALSE, recursive = TRUE)

## ---- Map theme (reused across plots) ----

map_theme <- theme_minimal(base_size = 10) +
  theme(
    legend.position  = "right",
    plot.title       = element_text(size = 12, face = "bold"),
    plot.subtitle    = element_text(size = 9, colour = "grey30"),
    panel.grid.major = element_line(colour = "grey85", linewidth = 0.2),
    axis.text        = element_text(size = 7)
  )

## ---- Copernicus Marine CLI ----

cmems_exe <- Sys.which("copernicusmarine")
if (nchar(cmems_exe) == 0) {
  cmems_exe <- "copernicusmarine"
  warning("copernicusmarine not on PATH — commands may fail")
}

download_cmems <- function(dataset_id, variable, start_date, end_date,
                           min_depth = NULL, max_depth = NULL,
                           output_dir = "data/threats_global/cmems/") {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  cmd <- paste(
    shQuote(cmems_exe),
    "subset",
    paste0("--dataset-id ", dataset_id),
    paste0("--variable ", variable),
    paste0("--start-datetime ", start_date, "T00:00:00"),
    paste0("--end-datetime ", end_date, "T00:00:00"),
    paste0("--minimum-longitude ", bb["xmin"]),
    paste0("--maximum-longitude ", bb["xmax"]),
    paste0("--minimum-latitude ", bb["ymin"]),
    paste0("--maximum-latitude ", bb["ymax"]),
    if (!is.null(min_depth)) paste0("--minimum-depth ", min_depth) else "",
    if (!is.null(max_depth)) paste0("--maximum-depth ", max_depth) else "",
    paste0("-o ", shQuote(output_dir)),
    "--force-download",
    sep = " "
  )
  cat("  Downloading:", variable, "from", dataset_id, "\n")
  cat("  Command:", cmd, "\n")
  result <- system(cmd)
  if (result != 0) warning("Download failed for ", variable, " (exit code ", result, ")")
  invisible(result)
}

cat("\n  Setup complete.\n")

##=====================================================================
## CMEMS Dataset IDs — VERIFIED
##=====================================================================
## Browse/update at: https://data.marine.copernicus.eu/products
##
## Physical oceanography (Global reanalysis, ~9km, monthly):
##   cmems_mod_glo_phy_my_0.083deg_P1M-m
##   Variables: bottomT (no depth), thetao (0.5m), so (0.5m),
##              uo (0.5m), vo (0.5m)
##   Time: 2010-01-01 to 2023-12-31
##
## Biogeochemistry (NW Shelf reanalysis, 7km, monthly, one ID per var):
##   cmems_mod_nws_bgc-chl_my_7km-3D_P1M-m   — chl  (0.5m)
##   cmems_mod_nws_bgc-o2_my_7km-3D_P1M-m    — o2   (0.5m)
##   Time: 2003-01-01 to 2023-12-31
##
## Ocean colour (Atlantic, MY L3, 1km, daily):
##   cmems_obs-oc_atl_bgc-transp_my_l3-multi-1km_P1D
##   Variables: KD490, SPM (no depth)
##   Time: 2020-01-01 to 2023-12-31
##   NOTE: daily 1km — download in yearly chunks to avoid timeouts
##=====================================================================

phy_id  <- "cmems_mod_glo_phy_my_0.083deg_P1M-m"
bgc_chl <- "cmems_mod_nws_bgc-chl_my_7km-3D_P1M-m"
bgc_o2  <- "cmems_mod_nws_bgc-o2_my_7km-3D_P1M-m"
oc_id   <- "cmems_obs-oc_atl_bgc-transp_my_l3-multi-1km_P1D"

all_ids <- c(phy_id, bgc_chl, bgc_o2, oc_id)

## Verify dataset IDs exist before downloading
check_dataset <- function(dataset_id) {
  cmd <- paste(shQuote(cmems_exe), "describe --contains", dataset_id)
  result <- system(cmd, intern = TRUE, ignore.stderr = TRUE)
  found <- any(grepl(dataset_id, result, fixed = TRUE))
  if (!found) warning("Dataset ID not found: ", dataset_id)
  found
}

cat("\n=== Verifying CMEMS dataset IDs ===\n")
id_status <- sapply(all_ids, function(id) {
  ok <- check_dataset(id)
  cat(sprintf("  %s: %s\n", id, ifelse(ok, "OK", "NOT FOUND")))
  ok
})


##=====================================================================
## 1. DOWNLOAD — Global Physics (monthly)
##=====================================================================

cat("\n=== 1. Global Physics ===\n")

if (id_status[phy_id]) {

  ## Bottom temperature — covariate stack for site matching (20-year record)
  download_cmems(phy_id, "bottomT", "2003-01-01", "2023-12-31",
                 output_dir = "data/processed/physical/bottomT/")

  ## SST — environmental risk mapping (20-year record, flag months >=25C)
  download_cmems(phy_id, "thetao", "2003-01-01", "2023-12-31",
                 min_depth = 0.5, max_depth = 0.5,
                 output_dir = "data/processed/physical/thetao/")

  ## Salinity — seagrass potential model (6-year mean, >10 PSU threshold)
  download_cmems(phy_id, "so", "2018-01-01", "2023-12-31",
                 min_depth = 0.5, max_depth = 0.5,
                 output_dir = "data/processed/physical/so/")

  ## Eastward current — covariate stack for site matching
  download_cmems(phy_id, "uo", "2018-01-01", "2023-12-31",
                 min_depth = 0.5, max_depth = 0.5,
                 output_dir = "data/processed/physical/uo/")

  ## Northward current — covariate stack for site matching
  download_cmems(phy_id, "vo", "2018-01-01", "2023-12-31",
                 min_depth = 0.5, max_depth = 0.5,
                 output_dir = "data/processed/physical/vo/")

} else {
  cat("  SKIPPED — dataset ID not found\n")
}


##=====================================================================
## 1b. DOWNLOAD — IBI Daily Salinity (2013–2023)
##=====================================================================

cat("\n=== 1b. IBI Daily Salinity ===\n")

ibi_sal_id <- "cmems_mod_ibi_phy-sal_my_0.027deg_P1D-m"

ibi_sal_ok <- check_dataset(ibi_sal_id)
cat(sprintf("  %s: %s\n", ibi_sal_id, ifelse(ibi_sal_ok, "OK", "NOT FOUND")))

if (ibi_sal_ok) {
  ## Download year-by-year to avoid timeouts (daily data is large)
  for (yr in 2013:2023) {
    cat(sprintf("  Year %d...\n", yr))
    download_cmems(ibi_sal_id, "so",
                   paste0(yr, "-01-01"), paste0(yr, "-12-31"),
                   min_depth = 0.5, max_depth = 0.5,
                   output_dir = "data/processed/physical/so_daily_ibi/")
  }
} else {
  cat("  SKIPPED — dataset ID not found\n")
}


##=====================================================================
## 2. DOWNLOAD — NW Shelf Biogeochemistry (monthly, 2003–2023)
##=====================================================================

cat("\n=== 2. NW Shelf Biogeochemistry ===\n")

## Chlorophyll-a — environmental risk mapping (20-year record, flag 95th percentile exceedances)
if (id_status[bgc_chl]) {
  download_cmems(bgc_chl, "chl", "2003-01-01", "2023-12-31",
                 min_depth = 0.5, max_depth = 0.5,
                 output_dir = "data/processed/physical/chl/")
} else {
  cat("  SKIPPED chl — dataset ID not found\n")
}

## Dissolved oxygen — environmental risk mapping (20-year record, flag 5th percentile lows)
if (id_status[bgc_o2]) {
  download_cmems(bgc_o2, "o2", "2003-01-01", "2023-12-31",
                 min_depth = 0.5, max_depth = 0.5,
                 output_dir = "data/processed/physical/o2/")
} else {
  cat("  SKIPPED o2 — dataset ID not found\n")
}


##=====================================================================
## 3. DOWNLOAD — Atlantic Ocean Colour (daily 1km, 2020–2023)
##=====================================================================
## Daily 1km data — downloaded in yearly chunks to avoid timeouts

cat("\n=== 3. Atlantic Ocean Colour ===\n")

if (id_status[oc_id]) {

  ## KD490 — seagrass potential model (6-year mean, turbidity threshold)
  download_cmems(oc_id, "KD490", "2020-01-01", "2020-12-31",
                 output_dir = "data/processed/physical/KD490/")
  download_cmems(oc_id, "KD490", "2021-01-01", "2021-12-31",
                 output_dir = "data/processed/physical/KD490/")
  download_cmems(oc_id, "KD490", "2022-01-01", "2022-12-31",
                 output_dir = "data/processed/physical/KD490/")
  download_cmems(oc_id, "KD490", "2023-01-01", "2023-12-31",
                 output_dir = "data/processed/physical/KD490/")

  ## SPM — seagrass potential model (6-year mean, turbidity threshold)
  download_cmems(oc_id, "SPM", "2020-01-01", "2020-12-31",
                 output_dir = "data/processed/physical/SPM/")
  download_cmems(oc_id, "SPM", "2021-01-01", "2021-12-31",
                 output_dir = "data/processed/physical/SPM/")
  download_cmems(oc_id, "SPM", "2022-01-01", "2022-12-31",
                 output_dir = "data/processed/physical/SPM/")
  download_cmems(oc_id, "SPM", "2023-01-01", "2023-12-31",
                 output_dir = "data/processed/physical/SPM/")

} else {
  cat("  SKIPPED — dataset ID not found\n")
}

cat("\n=== CMEMS downloads complete ===\n")


##=====================================================================
## 3b. BED SHEAR STRESS — Kinetic Energy (currents)
##=====================================================================

cat("\n=== 3b. Bed Shear Stress (KE currents) ===\n")

ke_path <- "data/raw/physical/ke_currents_atlantic_withconf/ke_currents_atlantic.tif"
stopifnot("KE currents raster not found" = file.exists(ke_path))

ke_raw <- rast(ke_path)
cat(sprintf("  Raw raster: %d x %d | CRS: %s\n",
            ncol(ke_raw), nrow(ke_raw), crs(ke_raw, describe = TRUE)$code))
cat(sprintf("  Extent: %.4f to %.4f (x), %.4f to %.4f (y)\n",
            ext(ke_raw)[1], ext(ke_raw)[2], ext(ke_raw)[3], ext(ke_raw)[4]))

aoi_ke <- vect(st_transform(scamp_aoi_sf, crs(ke_raw)))
ke_crop <- crop(ke_raw, aoi_ke, mask = TRUE)

cat(sprintf("  Cropped raster: %d x %d\n", ncol(ke_crop), nrow(ke_crop)))
cat(sprintf("  Value range: %.4f to %.4f\n",
            global(ke_crop, "min", na.rm = TRUE)[[1]],
            global(ke_crop, "max", na.rm = TRUE)[[1]]))

ke_bng <- project(ke_crop, "EPSG:27700")

dir.create("outputs/physical", showWarnings = FALSE, recursive = TRUE)
writeRaster(ke_bng, "outputs/physical/ke_currents_solway.tif", overwrite = TRUE)
cat("  Saved: outputs/physical/ke_currents_solway.tif\n")

p_ke <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
          linewidth = 0.3) +
  tidyterra::geom_spatraster(data = ke_bng) +
  scale_fill_viridis_c(name = expression(N/m^2), na.value = NA, option = "inferno") +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  annotation_north_arrow(location = "tl", style = north_arrow_minimal(),
                         height = unit(1, "cm"), width = unit(1, "cm")) +
  labs(title = "Bed Shear Stress — Kinetic Energy (currents)",
       subtitle = sprintf("Range: %.2f – %.2f N/m²",
                          global(ke_bng, "min", na.rm = TRUE)[[1]],
                          global(ke_bng, "max", na.rm = TRUE)[[1]])) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "right",
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 9, colour = "grey30")) +
  coord_sf(crs = 27700)

dir.create("figures/physical", showWarnings = FALSE, recursive = TRUE)
ggsave("figures/physical/ke_currents_solway.png", p_ke,
       width = 14, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/physical/ke_currents_solway.png\n")

cat("\n=== Bed Shear Stress complete ===\n")


##=====================================================================
## 4. BATHYMETRY — GEBCO 2026
##=====================================================================

cat("\n=== 4. GEBCO 2026 Bathymetry ===\n")

bathy_path <- "data/raw/physical/GEBCO_2026/GEBCO_2026.nc"
bathy_mask <- NULL

stopifnot("GEBCO 2026 not found" = file.exists(bathy_path))

## terra misreads GEBCO NetCDF axis ordering — extract with ncdf4 directly
nc <- nc_open(bathy_path)
lon_all <- ncvar_get(nc, "lon")
lat_all <- ncvar_get(nc, "lat")

aoi_wgs <- st_bbox(st_transform(scamp_aoi_sf, 4326))
ix <- which(lon_all >= aoi_wgs["xmin"] & lon_all <= aoi_wgs["xmax"])
iy <- which(lat_all >= aoi_wgs["ymin"] & lat_all <= aoi_wgs["ymax"])

elev_mat <- ncvar_get(nc, "elevation",
                      start = c(min(ix), min(iy)),
                      count = c(length(ix), length(iy)))
nc_close(nc)

lon_sub <- lon_all[ix]
lat_sub <- lat_all[iy]
res_deg <- abs(lon_all[2] - lon_all[1])

bathy_crop <- rast(t(elev_mat[, ncol(elev_mat):1]),
                   extent = ext(min(lon_sub) - res_deg/2, max(lon_sub) + res_deg/2,
                                min(lat_sub) - res_deg/2, max(lat_sub) + res_deg/2),
                   crs = "EPSG:4326")

cat(sprintf("  Resolution: %.4f deg\n", res_deg))
cat(sprintf("  Value range: %d to %d m\n",
            global(bathy_crop, "min", na.rm = TRUE)[[1]],
            global(bathy_crop, "max", na.rm = TRUE)[[1]]))

aoi_v <- vect(st_transform(scamp_aoi_sf, 4326))
bathy_crop <- mask(bathy_crop, aoi_v)
bathy_mask <- project(bathy_crop, "EPSG:27700")
names(bathy_mask) <- "depth"

writeRaster(bathy_mask, "outputs/physical/gebco_bathy_solway.tif",
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))
cat("  Saved: outputs/physical/gebco_bathy_solway.tif\n")

bathy_df <- as.data.frame(bathy_mask, xy = TRUE)

p_bathy <- ggplot() +
  geom_tile(data = bathy_df, aes(x = x, y = y, fill = depth)) +
  scale_fill_viridis_c(name = "Depth (m)", option = "mako", direction = -1,
                       na.value = "transparent") +
  geom_sf(data = las_dg_bng, fill = NA, colour = "grey50", linewidth = 0.3) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  annotation_north_arrow(location = "tl", style = north_arrow_minimal(),
                         height = unit(1, "cm"), width = unit(1, "cm")) +
  labs(title = "GEBCO 2026 Bathymetry (SCAMP AOI)",
       subtitle = sprintf("Resolution: %.0fm", res(bathy_mask)[1])) +
  map_theme + coord_sf(crs = 27700)

ggsave("figures/physical/gebco_bathy_solway.png", p_bathy,
       width = 14, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/physical/gebco_bathy_solway.png\n")

gc()


##=====================================================================
## 5. ENERGY CLASS — EUSeaMap 2023
##=====================================================================

cat("\n=== 5. EUSeaMap 2023 Energy ===\n")

energy_gdb <- "data/raw/physical/EUSeaMap_2023_Energy/EUSeaMap_2023_Energy.gdb"

stopifnot("EUSeaMap 2023 gdb not found" = file.exists(energy_gdb))

cat("  Layers in gdb:\n")
print(st_layers(energy_gdb))

energy_raw <- st_read(energy_gdb, quiet = TRUE)
cat(sprintf("  CRS: %s\n", st_crs(energy_raw)$input))
cat(sprintf("  Columns: %s\n", paste(names(energy_raw), collapse = ", ")))
cat(sprintf("  Rows: %d\n", nrow(energy_raw)))

energy_bng <- st_transform(energy_raw, 27700) %>% st_make_valid()
energy_clip <- st_intersection(energy_bng, scamp_aoi_sf)
cat(sprintf("  Clipped rows: %d\n", nrow(energy_clip)))

cat("\n  Unique values per column:\n")
for (col in names(energy_clip)[!names(energy_clip) %in% c("geometry", "Shape")]) {
  vals <- unique(energy_clip[[col]])
  if (length(vals) <= 30) {
    cat(sprintf("    %s: %s\n", col, paste(vals, collapse = " | ")))
  } else {
    cat(sprintf("    %s: %d unique values (too many to list)\n", col, length(vals)))
  }
}

st_write(energy_clip, "outputs/physical/energy_class_solway.gpkg",
         delete_dsn = TRUE, quiet = TRUE)
cat("  Saved: outputs/physical/energy_class_solway.gpkg\n")

energy_col <- names(energy_clip)[grepl("energy|class", names(energy_clip),
                                        ignore.case = TRUE)][1]
cat(sprintf("  Energy column detected: %s\n", energy_col))

p_energy <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
          linewidth = 0.3) +
  geom_sf(data = energy_clip, aes(fill = .data[[energy_col]]),
          colour = NA, alpha = 0.7) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  scale_fill_brewer(palette = "Set2", name = energy_col) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  annotation_north_arrow(location = "tl", style = north_arrow_minimal(),
                         height = unit(1, "cm"), width = unit(1, "cm")) +
  labs(title = "EUSeaMap 2023 Energy Classification (SCAMP AOI)",
       subtitle = sprintf("%d features | Column: %s", nrow(energy_clip), energy_col)) +
  map_theme + coord_sf(crs = 27700)

ggsave("figures/physical/energy_class_solway.png", p_energy,
       width = 14, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/physical/energy_class_solway.png\n")

gc()


##=====================================================================
## 6. UKSeaMap v2026
##=====================================================================

cat("\n=== 6. UKSeaMap v2026 ===\n")

ukseamap_path <- "data/raw/physical/UKSeaMap-Predictive-Map-v2026-1/UKSeaMap_Predictive_Map_v2026_1/UKSeaMap_Predictive_Map_v2026_1.gpkg"

stopifnot("UKSeaMap gpkg not found" = file.exists(ukseamap_path))

ukseamap_lyrs <- st_layers(ukseamap_path)
cat("  Layers:\n")
print(ukseamap_lyrs)

## Transform AOI to UKSeaMap CRS once (WGS 84), then clip each layer
aoi_wgs84 <- st_transform(scamp_aoi_sf, 4326) %>% st_make_valid()

ukseamap_layers <- list()
for (lyr in ukseamap_lyrs$name) {
  cat(sprintf("\n  Reading: %s\n", lyr))
  raw <- st_read(ukseamap_path, layer = lyr, quiet = TRUE)
  clipped <- st_intersection(raw, aoi_wgs84) %>% st_make_valid()
  clipped <- st_transform(clipped, 27700)
  cat(sprintf("    Clipped: %d features | Columns: %s\n",
              nrow(clipped), paste(names(clipped), collapse = ", ")))

  short_name <- gsub("UKSeaMap_|_v2026_1|Habitat_Descriptor_", "", lyr)
  ukseamap_layers[[short_name]] <- clipped

  out_file <- sprintf("outputs/physical/ukseamap_%s_solway.gpkg", tolower(short_name))
  st_write(clipped, out_file, delete_dsn = TRUE, quiet = TRUE)
  cat(sprintf("    Saved: %s\n", out_file))
}

## Plot each layer
for (nm in names(ukseamap_layers)) {
  dat <- ukseamap_layers[[nm]]
  geom_cols <- c("geometry", "Shape", "geom")
  class_cols <- names(dat)[!names(dat) %in% c(geom_cols, "updatedate", "changetype",
                                                "jnccid", "fid")]
  fill_col <- class_cols[1]
  if (is.na(fill_col) || is.null(fill_col)) next

  n_vals <- length(unique(dat[[fill_col]]))
  cat(sprintf("\n  Plotting %s — column: %s (%d classes)\n", nm, fill_col, n_vals))

  p <- ggplot() +
    geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
    geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
            linewidth = 0.3) +
    geom_sf(data = dat, aes(fill = .data[[fill_col]]),
            colour = NA, alpha = 0.7) +
    geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
    annotation_scale(location = "bl", width_hint = 0.2) +
    annotation_north_arrow(location = "tl", style = north_arrow_minimal(),
                           height = unit(1, "cm"), width = unit(1, "cm")) +
    labs(title = sprintf("UKSeaMap %s (SCAMP AOI)", nm),
         subtitle = sprintf("%d features | %s", nrow(dat), fill_col)) +
    map_theme + coord_sf(crs = 27700)

  if (n_vals <= 12) {
    p <- p + scale_fill_brewer(palette = "Set3", name = fill_col)
  } else {
    p <- p + scale_fill_viridis_d(name = fill_col, option = "turbo")
  }

  fig_file <- sprintf("figures/physical/ukseamap_%s_solway.png", tolower(nm))
  ggsave(fig_file, p, width = 14, height = 10, units = "in", dpi = 300)
  cat(sprintf("    Saved: %s\n", fig_file))
}

cat("\n  UKSeaMap processing complete\n")

gc()


##=====================================================================
## 7. SEA LEVEL RISE — UKCP18 (RCP4.5)
##=====================================================================

cat("\n=== 7. UKCP18 Sea Level Rise ===\n")

slr_path <- "data/raw/physical/seaLevelAnom_marine-sim_rcp45_ann_2007-2100.nc"

stopifnot("UKCP18 sea level file not found" = file.exists(slr_path))

slr_nc <- nc_open(slr_path)
cat("  NetCDF structure:\n")
cat(sprintf("    Variables: %s\n", paste(names(slr_nc$var), collapse = ", ")))
cat(sprintf("    Dimensions: %s\n",
            paste(sapply(slr_nc$dim, function(d) sprintf("%s [%d]", d$name, d$len)),
                  collapse = ", ")))
for (v in names(slr_nc$var)) {
  cat(sprintf("    %s: units=%s, dims=%s\n", v,
              slr_nc$var[[v]]$units,
              paste(sapply(slr_nc$var[[v]]$dim, `[[`, "name"), collapse = ", ")))
}
nc_close(slr_nc)

slr_raw <- rast(slr_path)
cat(sprintf("  Layers: %d\n", nlyr(slr_raw)))
cat(sprintf("  CRS: %s\n", crs(slr_raw, describe = TRUE)$code))
cat(sprintf("  Resolution: %.4f x %.4f\n", res(slr_raw)[1], res(slr_raw)[2]))

slr_crs <- crs(slr_raw)
if (is.na(slr_crs) || slr_crs == "") {
  crs(slr_raw) <- "EPSG:4326"
  cat("  CRS missing — assigned EPSG:4326\n")
}

aoi_slr <- vect(st_transform(scamp_aoi_sf, crs(slr_raw)))
slr_crop <- crop(slr_raw, aoi_slr)
slr_crop <- mask(slr_crop, aoi_slr)
slr_bng <- project(slr_crop, "EPSG:27700")

writeRaster(slr_bng, "outputs/physical/ukcp18_sea_level_solway.tif",
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))
cat("  Saved: outputs/physical/ukcp18_sea_level_solway.tif\n")

slr_last <- slr_bng[[nlyr(slr_bng)]]
slr_df <- as.data.frame(slr_last, xy = TRUE)
names(slr_df)[3] <- "slr"

p_slr <- ggplot() +
  geom_tile(data = slr_df, aes(x = x, y = y, fill = slr)) +
  scale_fill_viridis_c(name = "Sea level\nanomaly (m)", option = "plasma",
                       na.value = "transparent") +
  geom_sf(data = las_dg_bng, fill = NA, colour = "grey50", linewidth = 0.3) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  annotation_north_arrow(location = "tl", style = north_arrow_minimal(),
                         height = unit(1, "cm"), width = unit(1, "cm")) +
  labs(title = "UKCP18 Sea Level Anomaly — RCP4.5 (SCAMP AOI)",
       subtitle = sprintf("Final year layer (%d of %d) | Resolution: %.2f deg",
                          nlyr(slr_raw), nlyr(slr_raw), res(slr_raw)[1])) +
  map_theme + coord_sf(crs = 27700)

ggsave("figures/physical/ukcp18_sea_level_solway.png", p_slr,
       width = 14, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/physical/ukcp18_sea_level_solway.png\n")

gc()

cat("\n=== Manual dataset processing complete ===\n")


##=====================================================================
## 8. SEAGRASS POTENTIAL MODEL (EA v3 methodology)
##=====================================================================

cat("\n=== 8. Seagrass Potential Model ===\n")

dir.create("outputs/habitats", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/habitats", showWarnings = FALSE, recursive = TRUE)

## bathy_mask already exists from section 4 (GEBCO, EPSG:27700)
bathy_r <- bathy_mask
cat(sprintf("  Bathymetry grid: %d x %d | res: %.0fm\n",
            ncol(bathy_r), nrow(bathy_r), res(bathy_r)[1]))

## ---- 8a. Compute SPM mean from CMEMS daily NetCDFs ----

cat("\n-- 8a. SPM mean (2020-2023) --\n")

spm_files <- sort(list.files("data/processed/physical/SPM", pattern = "\\.nc$",
                              full.names = TRUE))
cat(sprintf("  SPM files found: %d\n", length(spm_files)))
stopifnot("No SPM NetCDF files found" = length(spm_files) > 0)

spm_stack <- rast(spm_files)
cat(sprintf("  SPM stack: %d layers\n", nlyr(spm_stack)))

spm_mean_wgs <- app(spm_stack, fun = "mean", na.rm = TRUE)
spm_mean_bng <- project(spm_mean_wgs, bathy_r, method = "bilinear")
names(spm_mean_bng) <- "spm_mean"

writeRaster(spm_mean_bng, "outputs/physical/spm_mean_solway.tif",
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))
cat(sprintf("  SPM mean range: %.2f to %.2f g/m3\n",
            global(spm_mean_bng, "min", na.rm = TRUE)[[1]],
            global(spm_mean_bng, "max", na.rm = TRUE)[[1]]))
cat("  Saved: outputs/physical/spm_mean_solway.tif\n")

spm_df <- as.data.frame(spm_mean_bng, xy = TRUE)

p_spm <- ggplot() +
  geom_tile(data = spm_df, aes(x = x, y = y, fill = spm_mean)) +
  scale_fill_viridis_c(name = "SPM (g/m³)", option = "inferno",
                       na.value = "transparent") +
  geom_sf(data = las_dg_bng, fill = NA, colour = "grey50", linewidth = 0.3) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  annotation_north_arrow(location = "tl", style = north_arrow_minimal(),
                         height = unit(1, "cm"), width = unit(1, "cm")) +
  labs(title = "Mean Suspended Particulate Matter (SCAMP AOI)",
       subtitle = sprintf("CMEMS 2020-2023 daily mean | Range: %.1f–%.1f g/m³",
                          global(spm_mean_bng, "min", na.rm = TRUE)[[1]],
                          global(spm_mean_bng, "max", na.rm = TRUE)[[1]])) +
  map_theme + coord_sf(crs = 27700)

ggsave("figures/physical/spm_mean_solway.png", p_spm,
       width = 14, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/physical/spm_mean_solway.png\n")

## ---- 8b. Compute salinity: mean of annual 90th percentiles ----

cat("\n-- 8b. Salinity — mean of annual P90 (2018-2023) --\n")

so_file <- list.files("data/processed/physical/so", pattern = "\\.nc$",
                       full.names = TRUE)[1]
stopifnot("No salinity NetCDF found" = !is.na(so_file) && file.exists(so_file))

so_stack <- rast(so_file)
cat(sprintf("  Salinity stack: %d monthly layers\n", nlyr(so_stack)))

## Extract year from layer names or timestamps
so_time <- time(so_stack)
so_years <- as.numeric(format(so_time, "%Y"))
unique_years <- sort(unique(so_years))
cat(sprintf("  Years: %s\n", paste(unique_years, collapse = ", ")))

## Compute 90th percentile per year
annual_p90 <- rast()
for (yr in unique_years) {
  yr_layers <- so_stack[[which(so_years == yr)]]
  p90 <- app(yr_layers, fun = function(x) quantile(x, 0.9, na.rm = TRUE))
  names(p90) <- paste0("p90_", yr)
  annual_p90 <- c(annual_p90, p90)
  cat(sprintf("    %d: %d months -> P90 range %.1f–%.1f\n", yr, nlyr(yr_layers),
              global(p90, "min", na.rm = TRUE)[[1]],
              global(p90, "max", na.rm = TRUE)[[1]]))
}

## Mean of annual P90s
sal_mean_p90 <- app(annual_p90, fun = "mean", na.rm = TRUE)
sal_mean_bng <- project(sal_mean_p90, bathy_r, method = "bilinear")
names(sal_mean_bng) <- "sal_p90_mean"

writeRaster(sal_mean_bng, "outputs/physical/salinity_p90_mean_solway.tif",
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))
cat(sprintf("  Salinity P90 mean range: %.1f to %.1f PSU\n",
            global(sal_mean_bng, "min", na.rm = TRUE)[[1]],
            global(sal_mean_bng, "max", na.rm = TRUE)[[1]]))
cat("  Saved: outputs/physical/salinity_p90_mean_solway.tif\n")

sal_df <- as.data.frame(sal_mean_bng, xy = TRUE)

p_sal <- ggplot() +
  geom_tile(data = sal_df, aes(x = x, y = y, fill = sal_p90_mean)) +
  scale_fill_viridis_c(name = "Salinity (PSU)", option = "cividis",
                       na.value = "transparent") +
  geom_sf(data = las_dg_bng, fill = NA, colour = "grey50", linewidth = 0.3) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  annotation_north_arrow(location = "tl", style = north_arrow_minimal(),
                         height = unit(1, "cm"), width = unit(1, "cm")) +
  labs(title = "Salinity — Mean of Annual 90th Percentiles (SCAMP AOI)",
       subtitle = sprintf("CMEMS 2018-2023 monthly | Range: %.1f–%.1f PSU",
                          global(sal_mean_bng, "min", na.rm = TRUE)[[1]],
                          global(sal_mean_bng, "max", na.rm = TRUE)[[1]])) +
  map_theme + coord_sf(crs = 27700)

ggsave("figures/physical/salinity_p90_mean_solway.png", p_sal,
       width = 14, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/physical/salinity_p90_mean_solway.png\n")

## ---- STEP 1: ENERGY FILTER ----

cat("\n-- Step 1: Energy filter --\n")

## energy_clip or ukseamap_layers$Energy_Regime may already be in memory from section 5/6
## Read from the EUSeaMap energy gpkg produced in section 5
energy_sf <- st_read("outputs/physical/energy_class_solway.gpkg", quiet = TRUE)
if (st_crs(energy_sf)$epsg != 27700) energy_sf <- st_transform(energy_sf, 27700)

cat(sprintf("  Energy features: %d\n", nrow(energy_sf)))
energy_cols <- names(energy_sf)[grepl("energy|class", names(energy_sf), ignore.case = TRUE)]
cat(sprintf("  Energy columns: %s\n", paste(energy_cols, collapse = ", ")))
for (ec in energy_cols) {
  cat(sprintf("    %s: %s\n", ec, paste(sort(unique(energy_sf[[ec]])), collapse = " | ")))
}

## If wave and current are separate, both must be "Low"
wave_col <- energy_cols[grepl("wave", energy_cols, ignore.case = TRUE)]
curr_col <- energy_cols[grepl("current", energy_cols, ignore.case = TRUE)]

if (length(wave_col) > 0 && length(curr_col) > 0) {
  cat(sprintf("  Separate wave (%s) and current (%s) columns\n", wave_col[1], curr_col[1]))
  energy_low <- energy_sf %>%
    filter(grepl("Low", .data[[wave_col[1]]], ignore.case = TRUE) &
           grepl("Low", .data[[curr_col[1]]], ignore.case = TRUE))
} else {
  energy_low <- energy_sf %>%
    filter(grepl("Low", .data[[energy_cols[1]]], ignore.case = TRUE))
}

cat(sprintf("  Low-energy polygons: %d of %d\n", nrow(energy_low), nrow(energy_sf)))

energy_mask <- rasterize(vect(energy_low), bathy_r, field = 1, background = 0)
names(energy_mask) <- "energy"

cat(sprintf("  Energy mask: %d pixels pass (%.1f%%)\n",
            global(energy_mask, "sum", na.rm = TRUE)[[1]],
            global(energy_mask, "mean", na.rm = TRUE)[[1]] * 100))

## ---- STEP 2: DEPTH FILTER (turbidity-dependent thresholds) ----

cat("\n-- Step 2: Depth filter (turbidity-dependent) --\n")

cat(sprintf("  SPM mean range: %.1f to %.1f g/m3\n",
            global(spm_mean_bng, "min", na.rm = TRUE)[[1]],
            global(spm_mean_bng, "max", na.rm = TRUE)[[1]]))

## Low turbidity (SPM < 10): depth between -10m and +5m
## Higher turbidity (SPM >= 10): depth between -5m and +5m
## Where SPM is NA (no satellite coverage), assume low turbidity
low_turb   <- ifel(is.na(spm_mean_bng), 1, spm_mean_bng < 10)
depth_low  <- (bathy_r >= -10) & (bathy_r <= 5)
depth_high <- (bathy_r >= -5)  & (bathy_r <= 5)

depth_mask <- ifel(low_turb == 1, depth_low, depth_high)
names(depth_mask) <- "depth"

n_spm_na <- global(is.na(spm_mean_bng), "sum", na.rm = TRUE)[[1]]
cat(sprintf("  SPM NA pixels (assumed low turbidity): %d\n", n_spm_na))
cat(sprintf("  Low turbidity pixels (SPM < 10 or NA): %d | Higher: %d\n",
            global(low_turb, "sum", na.rm = TRUE)[[1]],
            global(low_turb == 0, "sum", na.rm = TRUE)[[1]]))
cat(sprintf("  Depth mask: %d pixels pass (%.1f%%)\n",
            global(depth_mask, "sum", na.rm = TRUE)[[1]],
            global(depth_mask, "mean", na.rm = TRUE)[[1]] * 100))

## ---- STEP 3: SALINITY FILTER ----

cat("\n-- Step 3: Salinity filter --\n")

cat(sprintf("  Salinity P90 mean range: %.1f to %.1f PSU\n",
            global(sal_mean_bng, "min", na.rm = TRUE)[[1]],
            global(sal_mean_bng, "max", na.rm = TRUE)[[1]]))

## Where salinity is NA (no satellite coverage in nearshore), assume passes
sal_mask <- ifel(is.na(sal_mean_bng), 1, ifel(sal_mean_bng >= 10, 1, 0))
names(sal_mask) <- "salinity"

n_sal_na <- global(is.na(sal_mean_bng), "sum", na.rm = TRUE)[[1]]
cat(sprintf("  Salinity NA pixels (assumed pass): %d\n", n_sal_na))
cat(sprintf("  Salinity mask: %d pixels pass (%.1f%%)\n",
            global(sal_mask, "sum", na.rm = TRUE)[[1]],
            global(sal_mask, "mean", na.rm = TRUE)[[1]] * 100))

## ---- STEP 4: COMBINE ----

cat("\n-- Step 4: Combine masks --\n")

sg_potential <- energy_mask * depth_mask * sal_mask
names(sg_potential) <- "potential"
cat(sprintf("  Combined: %d pixels pass all three filters (%.1f%%)\n",
            global(sg_potential, "sum", na.rm = TRUE)[[1]],
            global(sg_potential, "mean", na.rm = TRUE)[[1]] * 100))

## ---- STEP 5: EXCLUSIONS ----

cat("\n-- Step 5: Exclusions --\n")

## Exclude saltmarsh areas
saltmarsh_path <- "data/processed/species/saltmarsh_harmonised.gpkg"
stopifnot("Saltmarsh harmonised not found" = file.exists(saltmarsh_path))
sm_harm <- st_read(saltmarsh_path, quiet = TRUE)
if (st_crs(sm_harm)$epsg != 27700) sm_harm <- st_transform(sm_harm, 27700)
sm_rast <- rasterize(vect(sm_harm), bathy_r, field = 1, background = 0)
sg_potential <- ifel(sm_rast == 1, 0, sg_potential)
cat(sprintf("  After saltmarsh exclusion: %d pixels\n",
            global(sg_potential, "sum", na.rm = TRUE)[[1]]))

writeRaster(sg_potential, "outputs/habitats/seagrass_potential_raster.tif",
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))
cat("  Saved: outputs/habitats/seagrass_potential_raster.tif\n")

## ---- STEP 6: VECTORISE AND SUMMARISE ----

cat("\n-- Step 6: Vectorise and summarise --\n")

sg_poly <- as.polygons(sg_potential) %>% st_as_sf()
sg_poly <- sg_poly[sg_poly$potential == 1, ]
sg_poly <- st_transform(sg_poly, 27700) %>% st_make_valid()
sg_poly <- st_cast(sg_poly, "POLYGON")
sg_poly$area_ha <- as.numeric(st_area(sg_poly)) / 10000

cat(sprintf("  Potential seagrass polygons: %d\n", nrow(sg_poly)))
cat(sprintf("  Total potential area: %.1f ha\n", sum(sg_poly$area_ha)))
cat(sprintf("  Area range: %.2f to %.2f ha\n",
            min(sg_poly$area_ha), max(sg_poly$area_ha)))

st_write(sg_poly, "outputs/habitats/seagrass_potential_solway.gpkg",
         delete_dsn = TRUE, quiet = TRUE)
cat("  Saved: outputs/habitats/seagrass_potential_solway.gpkg\n")

## ---- Plot ----

cat("\n-- Seagrass potential map --\n")

p_sg_potential <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
          linewidth = 0.3) +
  geom_sf(data = sg_poly, fill = "seagreen", colour = NA, alpha = 0.6) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  annotation_north_arrow(location = "tl", style = north_arrow_minimal(),
                         height = unit(1, "cm"), width = unit(1, "cm")) +
  labs(title = "Seagrass Potential — EA v3 Methodology (SCAMP AOI)",
       subtitle = sprintf("%d polygons | %.1f ha total | Excludes saltmarsh",
                          nrow(sg_poly), sum(sg_poly$area_ha))) +
  map_theme + coord_sf(crs = 27700)

ggsave("figures/habitats/seagrass_potential_solway.png", p_sg_potential,
       width = 14, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/habitats/seagrass_potential_solway.png\n")

## Diagnostic summary
cat("\n-- Diagnostic summary --\n")
cat(sprintf("  Step 1 (energy):   %d pixels\n", global(energy_mask, "sum", na.rm = TRUE)[[1]]))
cat(sprintf("  Step 2 (depth):    %d pixels\n", global(depth_mask, "sum", na.rm = TRUE)[[1]]))
cat(sprintf("  Step 3 (salinity): %d pixels\n", global(sal_mask, "sum", na.rm = TRUE)[[1]]))
cat(sprintf("  Step 4 (combined): %d pixels\n",
            global(energy_mask * depth_mask * sal_mask, "sum", na.rm = TRUE)[[1]]))
cat(sprintf("  Step 5 (final):    %d pixels\n", global(sg_potential, "sum", na.rm = TRUE)[[1]]))

cat("\n=== Seagrass Potential Model complete ===\n")


##=====================================================================
## 9. NATIVE OYSTER POTENTIAL MODEL
##=====================================================================
##
## Sources: Hughes et al. (2023), Preston et al. (2020),
##          Environment Agency, Marine Scotland.
##=====================================================================

cat("\n\n=== 9. NATIVE OYSTER POTENTIAL MODEL ===\n")

## Reference grid: bathymetry raster (EPSG:27700)
oy_bathy <- rast("outputs/physical/gebco_bathy_solway.tif")
cat(sprintf("  Reference grid: %d x %d | res: %.0fm | CRS: EPSG:%s\n",
            ncol(oy_bathy), nrow(oy_bathy), res(oy_bathy)[1],
            crs(oy_bathy, describe = TRUE)$code))


## ---- STEP 1: DEPTH FILTER ----

cat("\n-- Step 1: Depth filter (0m to -80m) --\n")

oy_depth_mask <- ifel(oy_bathy >= -80 & oy_bathy <= 0, 1, 0)
names(oy_depth_mask) <- "depth"

cat(sprintf("  Depth range in AOI: %.0f to %.0f m\n",
            global(oy_bathy, "min", na.rm = TRUE)[[1]],
            global(oy_bathy, "max", na.rm = TRUE)[[1]]))
cat(sprintf("  Pixels passing (0 to -80m): %d (%.1f%%)\n",
            global(oy_depth_mask, "sum", na.rm = TRUE)[[1]],
            global(oy_depth_mask, "mean", na.rm = TRUE)[[1]] * 100))

p_oy_depth <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
          linewidth = 0.3) +
  tidyterra::geom_spatraster(data = as.factor(oy_depth_mask)) +
  scale_fill_manual(values = c("0" = "transparent", "1" = "#2166ac"),
                    labels = c("Excluded", "0 to -80m"), name = "Depth",
                    na.translate = FALSE) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  labs(title = "Oyster Potential — Step 1: Depth Filter",
       subtitle = sprintf("Retained: %d pixels (0 to -80m)",
                          global(oy_depth_mask, "sum", na.rm = TRUE)[[1]])) +
  map_theme + coord_sf(crs = 27700)

ggsave("figures/habitats/oyster_potential_step1_depth.png", p_oy_depth,
       width = 14, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/habitats/oyster_potential_step1_depth.png\n")


## ---- STEP 2: SUBSTRATE FILTER (BGS Predictive Seabed Sediments) ----
## UKSeaMap EUNIS A5.4 does not exist within the SCAMP AOI despite known
## present and historical oyster sites — use BGS instead.

cat("\n-- Step 2: Substrate filter (BGS Predictive Seabed Sediments) --\n")

bgs_path <- "data/raw/physical/Predictive_SBS_UK_V1_GeoPackage/BGS_Predictive_Seabed_Sediments_UK_v1.gpkg"
stopifnot("BGS sediment geopackage not found" = file.exists(bgs_path))

bgs_layers <- st_layers(bgs_path)
cat("  Available layers:\n")
for (i in seq_along(bgs_layers$name)) {
  cat(sprintf("    %d: %s (%d features, %s)\n",
              i, bgs_layers$name[i], bgs_layers$features[i],
              bgs_layers$geomtype[i]))
}

bgs_sf <- st_read(bgs_path, layer = bgs_layers$name[1], quiet = TRUE)
cat(sprintf("\n  Loaded layer: %s | %d features | CRS: %s\n",
            bgs_layers$name[1], nrow(bgs_sf),
            st_crs(bgs_sf)$input))
cat(sprintf("  Columns: %s\n", paste(names(bgs_sf), collapse = ", ")))

## Print unique values of every non-geometry column
geom_col <- attr(bgs_sf, "sf_column")
class_cols <- setdiff(names(bgs_sf), geom_col)
for (cc in class_cols) {
  vals <- sort(unique(bgs_sf[[cc]]))
  cat(sprintf("\n  %s (%d unique):\n", cc, length(vals)))
  for (v in vals) cat(sprintf("    - %s\n", v))
}

## Crop to SCAMP AOI — transform AOI to sediment CRS to avoid large reprojection
aoi_bgs <- st_transform(scamp_aoi_sf, st_crs(bgs_sf))
bgs_clip <- st_intersection(bgs_sf, aoi_bgs)
cat(sprintf("\n  Features within AOI: %d of %d\n", nrow(bgs_clip), nrow(bgs_sf)))

## Print unique sediment classes within AOI
cat("\n  === SEDIMENT CLASSES WITHIN SCAMP AOI ===\n")
for (cc in class_cols) {
  vals <- sort(unique(bgs_clip[[cc]]))
  cat(sprintf("\n  %s (%d unique):\n", cc, length(vals)))
  for (v in vals) {
    n <- sum(bgs_clip[[cc]] == v, na.rm = TRUE)
    cat(sprintf("    - %s  [%d polygons]\n", v, n))
  }
}

cat("\n  >>> STOP: Review the sediment classes above.\n")
cat("  >>> Confirm which classes to assign to each substrate tier before proceeding.\n")

## Quick diagnostic plot of sediment classes in AOI
bgs_clip_bng <- st_transform(bgs_clip, 27700)
sed_col <- class_cols[1]

p_bgs <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
          linewidth = 0.3) +
  geom_sf(data = bgs_clip_bng, aes(fill = .data[[sed_col]]),
          colour = NA, alpha = 0.7) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  labs(title = "BGS Predictive Seabed Sediments — SCAMP AOI",
       subtitle = sprintf("Column: %s | %d classes | %d polygons",
                          sed_col, length(unique(bgs_clip[[sed_col]])),
                          nrow(bgs_clip)),
       fill = sed_col) +
  map_theme + coord_sf(crs = 27700)

ggsave("figures/habitats/oyster_step2_bgs_sediments.png", p_bgs,
       width = 14, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/habitats/oyster_step2_bgs_sediments.png\n")

## Build substrate mask from Folk_Class tiers
excluded_classes <- c("Sand", "sandy Mud")

bgs_clip_bng$substrate_pass <- ifelse(bgs_clip_bng$Folk_Class %in% excluded_classes, 0, 1)

cat(sprintf("\n  Optimal (pass = 1): %d polygons\n",
            sum(bgs_clip_bng$substrate_pass == 1)))
cat(sprintf("  Excluded (pass = 0): %d polygons (%s)\n",
            sum(bgs_clip_bng$substrate_pass == 0),
            paste(excluded_classes, collapse = ", ")))

## Rasterise to reference grid
oy_substrate_mask <- rasterize(vect(bgs_clip_bng), oy_bathy,
                               field = "substrate_pass", background = 0)
names(oy_substrate_mask) <- "substrate"

cat(sprintf("  Substrate mask: %d pixels pass (%.1f%%)\n",
            global(oy_substrate_mask, "sum", na.rm = TRUE)[[1]],
            global(oy_substrate_mask, "mean", na.rm = TRUE)[[1]] * 100))

## Plot mask
p_oy_sub <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
          linewidth = 0.3) +
  tidyterra::geom_spatraster(data = as.factor(oy_substrate_mask)) +
  scale_fill_manual(values = c("0" = "#d73027", "1" = "#1a9850"),
                    labels = c("Excluded (Sand / sandy Mud)",
                               "Optimal (gravel/mud/rock mixes)"),
                    name = "Substrate", na.translate = FALSE) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  labs(title = "Oyster Potential — Step 2: Substrate Filter (BGS Folk Class)",
       subtitle = sprintf("Pass: %d pixels (%.1f%%) | Excluded: Sand, sandy Mud",
                          global(oy_substrate_mask, "sum", na.rm = TRUE)[[1]],
                          global(oy_substrate_mask, "mean", na.rm = TRUE)[[1]] * 100)) +
  map_theme + coord_sf(crs = 27700)

ggsave("figures/habitats/oyster_potential_step2_substrate.png", p_oy_sub,
       width = 14, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/habitats/oyster_potential_step2_substrate.png\n")


## ---- STEP 3: SALINITY FILTER (min of monthly means, IBI daily) ----

cat("\n-- Step 3: Salinity filter (min-of-monthly-means >= 20 PSU) --\n")

## Uses mean_min_monthly_bng computed in step 3b below
## If pre-computed raster exists, load it; otherwise step 3b must run first
sal_min_monthly_path <- "outputs/physical/sal_min_monthly_mean_solway.tif"
if (file.exists(sal_min_monthly_path)) {
  mean_min_monthly_bng <- rast(sal_min_monthly_path)
  mean_min_monthly_bng <- resample(mean_min_monthly_bng, oy_bathy, method = "bilinear")
  cat("  Loaded pre-computed min-of-monthly-means salinity\n")
} else {
  cat("  WARNING: sal_min_monthly_mean_solway.tif not found\n")
  cat("  Run step 3b (daily salinity metrics) first, or run the download\n")
  cat("  Setting salinity mask to all-pass (NA)\n")
  mean_min_monthly_bng <- oy_bathy * NA
}
names(mean_min_monthly_bng) <- "sal_min_monthly"

n_sal_na_before <- global(is.na(mean_min_monthly_bng), "sum", na.rm = TRUE)[[1]]
cat(sprintf("  NA pixels before gap-fill: %d\n", n_sal_na_before))

## Gap-fill nearshore NAs with nearest-neighbour interpolation
if (n_sal_na_before > 0) {
  sal_filled <- focal(mean_min_monthly_bng, w = 3, fun = "mean", na.policy = "only",
                      na.rm = TRUE)
  n_remaining <- global(is.na(sal_filled), "sum", na.rm = TRUE)[[1]]
  passes <- 1
  while (n_remaining > 0 && passes < 50) {
    sal_filled <- focal(sal_filled, w = 3, fun = "mean", na.policy = "only",
                        na.rm = TRUE)
    n_remaining <- global(is.na(sal_filled), "sum", na.rm = TRUE)[[1]]
    passes <- passes + 1
  }
  mean_min_monthly_bng <- mask(sal_filled, mean_min_monthly_bng, inverse = TRUE,
                               updatevalue = NA)
  mean_min_monthly_bng <- cover(mean_min_monthly_bng, sal_filled)
  ## Re-mask to bathymetry extent
  mean_min_monthly_bng <- mask(mean_min_monthly_bng, oy_bathy)
  n_sal_na_after <- global(is.na(mean_min_monthly_bng), "sum", na.rm = TRUE)[[1]]
  cat(sprintf("  Gap-filled NAs with nearest-neighbour (%d passes)\n", passes))
  cat(sprintf("  NA pixels after gap-fill: %d\n", n_sal_na_after))
} else {
  cat("  No NAs to gap-fill\n")
}

cat(sprintf("  Min-of-monthly-means range: %.1f to %.1f PSU\n",
            global(mean_min_monthly_bng, "min", na.rm = TRUE)[[1]],
            global(mean_min_monthly_bng, "max", na.rm = TRUE)[[1]]))

## Hard filter: min-of-monthly-means >= 20 PSU
oy_sal_mask <- ifel(is.na(mean_min_monthly_bng), 0,
                    ifel(mean_min_monthly_bng >= 20, 1, 0))
names(oy_sal_mask) <- "salinity"

## Optimum band: 25-35 PSU
oy_sal_optimum <- ifel(is.na(mean_min_monthly_bng), 0,
                       ifel(mean_min_monthly_bng >= 25 & mean_min_monthly_bng <= 35, 1, 0))
names(oy_sal_optimum) <- "sal_optimum"

n_sal_na <- global(is.na(mean_min_monthly_bng), "sum", na.rm = TRUE)[[1]]
cat(sprintf("  Salinity NA pixels (assumed pass): %d\n", n_sal_na))
cat(sprintf("  Hard filter (min-monthly >= 20 PSU): %d pixels pass (%.1f%%)\n",
            global(oy_sal_mask, "sum", na.rm = TRUE)[[1]],
            global(oy_sal_mask, "mean", na.rm = TRUE)[[1]] * 100))
cat(sprintf("  Optimum band (25-35 PSU): %d pixels\n",
            global(oy_sal_optimum, "sum", na.rm = TRUE)[[1]]))

p_oy_sal <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
          linewidth = 0.3) +
  tidyterra::geom_spatraster(data = mean_min_monthly_bng) +
  scale_fill_viridis_c(name = "PSU", na.value = NA, option = "cividis") +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  labs(title = "Oyster Potential — Step 3: Salinity (Min of Monthly Means)",
       subtitle = sprintf("Hard threshold: >= 20 PSU | Optimum: 25-35 PSU | Pass: %d pixels",
                          global(oy_sal_mask, "sum", na.rm = TRUE)[[1]])) +
  map_theme + coord_sf(crs = 27700)

ggsave("figures/habitats/oyster_potential_step3_salinity.png", p_oy_sal,
       width = 14, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/habitats/oyster_potential_step3_salinity.png\n")


## ---- STEP 3b: DAILY SALINITY METRICS (IBI 0.027°, 2013–2023) ----

cat("\n-- Step 3b: Daily salinity metrics (IBI model) --\n")

sal_daily_dir <- "data/processed/physical/so_daily_ibi"
sal_daily_files <- list.files(sal_daily_dir, pattern = "\\.nc$", full.names = TRUE)

if (length(sal_daily_files) > 0) {

  sal_daily_stack <- rast(sal_daily_files)
  cat(sprintf("  Daily salinity layers: %d\n", nlyr(sal_daily_stack)))

  ## Crop to AOI in source CRS
  aoi_sal_d <- vect(st_transform(scamp_aoi_sf, crs(sal_daily_stack)))
  sal_daily_crop <- crop(sal_daily_stack, aoi_sal_d, mask = TRUE)

  sal_d_time <- time(sal_daily_crop)
  sal_d_years <- as.numeric(format(sal_d_time, "%Y"))
  sal_d_months <- as.numeric(format(sal_d_time, "%m"))
  unique_sal_years <- sort(unique(sal_d_years))
  cat(sprintf("  Years: %s (%d years)\n",
              paste(range(unique_sal_years), collapse = "–"), length(unique_sal_years)))


  ## --- 3b-i: Mean days per year with PSU < 20 ---

  cat("\n  3b-i: Mean days below 20 PSU per year...\n")

  annual_days_below20 <- rast()
  for (yr in unique_sal_years) {
    yr_idx <- which(sal_d_years == yr)
    yr_layers <- sal_daily_crop[[yr_idx]]
    days_below <- app(yr_layers, fun = function(x) sum(x < 20, na.rm = TRUE))
    names(days_below) <- paste0("days_below20_", yr)
    annual_days_below20 <- c(annual_days_below20, days_below)
  }
  mean_days_below20 <- app(annual_days_below20, fun = "mean", na.rm = TRUE)
  names(mean_days_below20) <- "mean_days_below20"

  mean_days_below20_bng <- project(mean_days_below20, oy_bathy, method = "bilinear")

  cat(sprintf("  Mean days < 20 PSU range: %.1f to %.1f days/yr\n",
              global(mean_days_below20_bng, "min", na.rm = TRUE)[[1]],
              global(mean_days_below20_bng, "max", na.rm = TRUE)[[1]]))

  p_sal_days <- ggplot() +
    geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
    geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
            linewidth = 0.3) +
    tidyterra::geom_spatraster(data = mean_days_below20_bng) +
    scale_fill_viridis_c(name = "days/yr", na.value = NA, option = "inferno",
                         direction = -1) +
    geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
    labs(title = "Oyster Potential — Step 3b-i: Mean Days Below 20 PSU per Year",
         subtitle = sprintf("IBI daily salinity %d–%d | Range: %.1f–%.1f days/yr",
                            min(unique_sal_years), max(unique_sal_years),
                            global(mean_days_below20_bng, "min", na.rm = TRUE)[[1]],
                            global(mean_days_below20_bng, "max", na.rm = TRUE)[[1]])) +
    map_theme + coord_sf(crs = 27700)

  ggsave("figures/habitats/oyster_sal_days_below20.png", p_sal_days,
         width = 14, height = 10, units = "in", dpi = 300)
  cat("  Saved: figures/habitats/oyster_sal_days_below20.png\n")


  ## --- 3b-ii: Minimum of monthly means per year ---

  cat("\n  3b-ii: Minimum of monthly mean salinity per year...\n")

  annual_min_monthly <- rast()
  for (yr in unique_sal_years) {
    yr_months <- sort(unique(sal_d_months[sal_d_years == yr]))
    monthly_means <- rast()
    for (mo in yr_months) {
      mo_idx <- which(sal_d_years == yr & sal_d_months == mo)
      mo_mean <- app(sal_daily_crop[[mo_idx]], fun = "mean", na.rm = TRUE)
      names(mo_mean) <- paste0("mean_", yr, "_", sprintf("%02d", mo))
      monthly_means <- c(monthly_means, mo_mean)
    }
    yr_min <- app(monthly_means, fun = "min", na.rm = TRUE)
    names(yr_min) <- paste0("min_monthly_", yr)
    annual_min_monthly <- c(annual_min_monthly, yr_min)
  }
  mean_min_monthly <- app(annual_min_monthly, fun = "mean", na.rm = TRUE)
  names(mean_min_monthly) <- "mean_min_monthly_sal"

  mean_min_monthly_bng <- project(mean_min_monthly, oy_bathy, method = "bilinear")

  cat(sprintf("  Mean of annual min-monthly salinity: %.1f to %.1f PSU\n",
              global(mean_min_monthly_bng, "min", na.rm = TRUE)[[1]],
              global(mean_min_monthly_bng, "max", na.rm = TRUE)[[1]]))

  p_sal_min_monthly <- ggplot() +
    geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
    geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
            linewidth = 0.3) +
    tidyterra::geom_spatraster(data = mean_min_monthly_bng) +
    scale_fill_viridis_c(name = "PSU", na.value = NA, option = "cividis") +
    geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
    labs(title = "Oyster Potential — Step 3b-ii: Min of Monthly Mean Salinity",
         subtitle = sprintf("Mean of annual minimums %d–%d | Range: %.1f–%.1f PSU",
                            min(unique_sal_years), max(unique_sal_years),
                            global(mean_min_monthly_bng, "min", na.rm = TRUE)[[1]],
                            global(mean_min_monthly_bng, "max", na.rm = TRUE)[[1]])) +
    map_theme + coord_sf(crs = 27700)

  ggsave("figures/habitats/oyster_sal_min_monthly_mean.png", p_sal_min_monthly,
         width = 14, height = 10, units = "in", dpi = 300)
  cat("  Saved: figures/habitats/oyster_sal_min_monthly_mean.png\n")


  ## --- 3b-iii: P10 (10th percentile) salinity ---

  cat("\n  3b-iii: P10 (10th percentile) salinity...\n")

  sal_p10 <- app(sal_daily_crop, fun = function(x) {
    quantile(x, probs = 0.10, na.rm = TRUE)
  })
  names(sal_p10) <- "sal_p10"

  sal_p10_bng <- project(sal_p10, oy_bathy, method = "bilinear")

  cat(sprintf("  P10 salinity range: %.1f to %.1f PSU\n",
              global(sal_p10_bng, "min", na.rm = TRUE)[[1]],
              global(sal_p10_bng, "max", na.rm = TRUE)[[1]]))

  p_sal_p10 <- ggplot() +
    geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
    geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
            linewidth = 0.3) +
    tidyterra::geom_spatraster(data = sal_p10_bng) +
    scale_fill_viridis_c(name = "PSU", na.value = NA, option = "cividis") +
    geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
    labs(title = "Oyster Potential — Step 3b-iii: P10 Salinity",
         subtitle = sprintf("10th percentile of all daily values %d–%d | Range: %.1f–%.1f PSU",
                            min(unique_sal_years), max(unique_sal_years),
                            global(sal_p10_bng, "min", na.rm = TRUE)[[1]],
                            global(sal_p10_bng, "max", na.rm = TRUE)[[1]])) +
    map_theme + coord_sf(crs = 27700)

  ggsave("figures/habitats/oyster_sal_p10.png", p_sal_p10,
         width = 14, height = 10, units = "in", dpi = 300)
  cat("  Saved: figures/habitats/oyster_sal_p10.png\n")

  ## Save processed rasters
  writeRaster(mean_days_below20_bng,
              "outputs/physical/sal_mean_days_below20_solway.tif",
              overwrite = TRUE, gdal = c("COMPRESS=LZW"))
  writeRaster(mean_min_monthly_bng,
              "outputs/physical/sal_min_monthly_mean_solway.tif",
              overwrite = TRUE, gdal = c("COMPRESS=LZW"))
  writeRaster(sal_p10_bng,
              "outputs/physical/sal_p10_solway.tif",
              overwrite = TRUE, gdal = c("COMPRESS=LZW"))
  cat("  Saved processed salinity rasters to outputs/physical/\n")

} else {
  cat("  WARNING: No daily salinity files found in ", sal_daily_dir, "\n")
  cat("  Run download section 1b first\n")
}


## ---- STEP 4: CURRENT VELOCITY FILTER ----

cat("\n-- Step 4: Current velocity filter --\n")

uo_file <- list.files("data/processed/physical/uo", pattern = "\\.nc$",
                       full.names = TRUE)[1]
vo_file <- list.files("data/processed/physical/vo", pattern = "\\.nc$",
                       full.names = TRUE)[1]
stopifnot("uo NetCDF not found" = !is.na(uo_file) && file.exists(uo_file))
stopifnot("vo NetCDF not found" = !is.na(vo_file) && file.exists(vo_file))

uo_stack <- rast(uo_file)
vo_stack <- rast(vo_file)
cat(sprintf("  uo layers: %d | vo layers: %d\n", nlyr(uo_stack), nlyr(vo_stack)))

uo_mean <- app(uo_stack, fun = "mean", na.rm = TRUE)
vo_mean <- app(vo_stack, fun = "mean", na.rm = TRUE)

## Crop to AOI in source CRS before reprojecting
aoi_uv <- vect(st_transform(scamp_aoi_sf, crs(uo_mean)))
uo_crop <- crop(uo_mean, aoi_uv, mask = TRUE)
vo_crop <- crop(vo_mean, aoi_uv, mask = TRUE)

## Current speed
speed_wgs <- sqrt(uo_crop^2 + vo_crop^2)
names(speed_wgs) <- "current_speed"

speed_r <- project(speed_wgs, oy_bathy, method = "bilinear")
names(speed_r) <- "current_speed"

cat(sprintf("  Current speed range: %.3f to %.3f m/s\n",
            global(speed_r, "min", na.rm = TRUE)[[1]],
            global(speed_r, "max", na.rm = TRUE)[[1]]))

## Hard filter: speed <= 0.8 m/s passes (NA assumed pass)
oy_current_mask <- ifel(is.na(speed_r), 1, ifel(speed_r <= 0.8, 1, 0))
names(oy_current_mask) <- "current"

## Quality classification: optimum / marginal / excluded
## 0.25-0.45 = optimum (1), 0.45-0.8 = marginal (0.5), >0.8 = excluded (0)
oy_current_quality <- ifel(is.na(speed_r), NA,
  ifel(speed_r >= 0.25 & speed_r <= 0.45, 1,
    ifel(speed_r > 0.45 & speed_r <= 0.8, 0.5, 0)))
names(oy_current_quality) <- "current_quality"

cat(sprintf("  Hard filter (<= 0.8 m/s): %d pixels pass (%.1f%%)\n",
            global(oy_current_mask, "sum", na.rm = TRUE)[[1]],
            global(oy_current_mask, "mean", na.rm = TRUE)[[1]] * 100))

p_oy_current <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
          linewidth = 0.3) +
  tidyterra::geom_spatraster(data = speed_r) +
  scale_fill_viridis_c(name = "m/s", na.value = NA, option = "plasma") +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  labs(title = "Oyster Potential — Step 4: Current Speed",
       subtitle = sprintf("Range: %.3f–%.3f m/s | Hard threshold: <= 0.8 m/s | Pass: %d pixels",
                          global(speed_r, "min", na.rm = TRUE)[[1]],
                          global(speed_r, "max", na.rm = TRUE)[[1]],
                          global(oy_current_mask, "sum", na.rm = TRUE)[[1]])) +
  map_theme + coord_sf(crs = 27700)

ggsave("figures/habitats/oyster_potential_step4_current.png", p_oy_current,
       width = 14, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/habitats/oyster_potential_step4_current.png\n")


## ---- STEP 5: SUSPENDED SEDIMENT FILTER ----

cat("\n-- Step 5: SPM filter (< 60 mg/l) --\n")

oy_spm <- rast("outputs/physical/spm_mean_solway.tif")
oy_spm_r <- resample(oy_spm, oy_bathy, method = "bilinear")
names(oy_spm_r) <- "spm"

cat(sprintf("  SPM range: %.1f to %.1f mg/l\n",
            global(oy_spm_r, "min", na.rm = TRUE)[[1]],
            global(oy_spm_r, "max", na.rm = TRUE)[[1]]))

## NA assumed pass (nearshore satellite gaps)
oy_spm_mask <- ifel(is.na(oy_spm_r), 1, ifel(oy_spm_r < 60, 1, 0))
names(oy_spm_mask) <- "spm"

n_spm_na <- global(is.na(oy_spm_r), "sum", na.rm = TRUE)[[1]]
cat(sprintf("  SPM NA pixels (assumed pass): %d\n", n_spm_na))
cat(sprintf("  SPM filter (< 60 mg/l): %d pixels pass (%.1f%%)\n",
            global(oy_spm_mask, "sum", na.rm = TRUE)[[1]],
            global(oy_spm_mask, "mean", na.rm = TRUE)[[1]] * 100))

p_oy_spm <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
          linewidth = 0.3) +
  tidyterra::geom_spatraster(data = oy_spm_r) +
  scale_fill_viridis_c(name = "mg/l", na.value = NA, option = "inferno") +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  labs(title = "Oyster Potential — Step 5: Suspended Sediment",
       subtitle = sprintf("Range: %.1f–%.1f mg/l | Threshold: < 60 mg/l | Pass: %d pixels",
                          global(oy_spm_r, "min", na.rm = TRUE)[[1]],
                          global(oy_spm_r, "max", na.rm = TRUE)[[1]],
                          global(oy_spm_mask, "sum", na.rm = TRUE)[[1]])) +
  map_theme + coord_sf(crs = 27700)

ggsave("figures/habitats/oyster_potential_step5_spm.png", p_oy_spm,
       width = 14, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/habitats/oyster_potential_step5_spm.png\n")


## ---- STEP 6: TEMPERATURE FILTERS ----

cat("\n-- Step 6: Temperature filters (bottom temp 2003-2023) --\n")

bt_file <- list.files("data/processed/physical/bottomT", pattern = "\\.nc$",
                       full.names = TRUE)
bt_file <- bt_file[!grepl("\\(1\\)", bt_file)][1]
stopifnot("Bottom temp NetCDF not found" = !is.na(bt_file) && file.exists(bt_file))

bt_stack <- rast(bt_file)
cat(sprintf("  Bottom temp layers: %d\n", nlyr(bt_stack)))

bt_time <- time(bt_stack)
bt_years <- as.numeric(format(bt_time, "%Y"))
bt_months <- as.numeric(format(bt_time, "%m"))
unique_bt_years <- sort(unique(bt_years))
cat(sprintf("  Years: %s\n", paste(range(unique_bt_years), collapse = "–")))

## Crop to AOI in source CRS
aoi_bt <- vect(st_transform(scamp_aoi_sf, crs(bt_stack)))
bt_crop <- crop(bt_stack, aoi_bt, mask = TRUE)

## 6a: Summer max — mean of annual maximum monthly SST
cat("\n  6a: Summer max (mean of annual max monthly SST)...\n")
annual_max <- rast()
for (yr in unique_bt_years) {
  yr_layers <- bt_crop[[which(bt_years == yr)]]
  yr_max <- app(yr_layers, fun = "max", na.rm = TRUE)
  names(yr_max) <- paste0("max_", yr)
  annual_max <- c(annual_max, yr_max)
}
summer_max <- app(annual_max, fun = "mean", na.rm = TRUE)
names(summer_max) <- "summer_max"

summer_max_bng <- project(summer_max, oy_bathy, method = "bilinear")

cat(sprintf("  Summer max range: %.1f to %.1f °C\n",
            global(summer_max_bng, "min", na.rm = TRUE)[[1]],
            global(summer_max_bng, "max", na.rm = TRUE)[[1]]))

## Summer max filter: >= 15 AND < 28
oy_summer_mask <- ifel(is.na(summer_max_bng), 1,
                       ifel(summer_max_bng >= 15 & summer_max_bng < 28, 1, 0))
names(oy_summer_mask) <- "summer_temp"

cat(sprintf("  Summer max filter (15-28°C): %d pixels pass\n",
            global(oy_summer_mask, "sum", na.rm = TRUE)[[1]]))

## 6b: Winter min — mean of annual minimum monthly SST
cat("\n  6b: Winter min (mean of annual min monthly SST)...\n")
annual_min <- rast()
for (yr in unique_bt_years) {
  yr_layers <- bt_crop[[which(bt_years == yr)]]
  yr_min <- app(yr_layers, fun = "min", na.rm = TRUE)
  names(yr_min) <- paste0("min_", yr)
  annual_min <- c(annual_min, yr_min)
}
winter_min <- app(annual_min, fun = "mean", na.rm = TRUE)
names(winter_min) <- "winter_min"

winter_min_bng <- project(winter_min, oy_bathy, method = "bilinear")

cat(sprintf("  Winter min range: %.1f to %.1f °C\n",
            global(winter_min_bng, "min", na.rm = TRUE)[[1]],
            global(winter_min_bng, "max", na.rm = TRUE)[[1]]))

## Winter min filter: > 5°C
oy_winter_mask <- ifel(is.na(winter_min_bng), 1,
                       ifel(winter_min_bng > 5, 1, 0))
names(oy_winter_mask) <- "winter_temp"

cat(sprintf("  Winter min filter (> 5°C): %d pixels pass\n",
            global(oy_winter_mask, "sum", na.rm = TRUE)[[1]]))

## Combined temperature mask
oy_temp_mask <- oy_summer_mask * oy_winter_mask
names(oy_temp_mask) <- "temperature"

cat(sprintf("  Combined temp filter: %d pixels pass (%.1f%%)\n",
            global(oy_temp_mask, "sum", na.rm = TRUE)[[1]],
            global(oy_temp_mask, "mean", na.rm = TRUE)[[1]] * 100))

p_oy_temp <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
          linewidth = 0.3) +
  tidyterra::geom_spatraster(data = summer_max_bng) +
  scale_fill_viridis_c(name = "°C", na.value = NA, option = "magma") +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  labs(title = "Oyster Potential — Step 6: Temperature",
       subtitle = sprintf("Summer max: %.1f–%.1f°C | Winter min: %.1f–%.1f°C | Temp pass: %d px",
                          global(summer_max_bng, "min", na.rm = TRUE)[[1]],
                          global(summer_max_bng, "max", na.rm = TRUE)[[1]],
                          global(winter_min_bng, "min", na.rm = TRUE)[[1]],
                          global(winter_min_bng, "max", na.rm = TRUE)[[1]],
                          global(oy_temp_mask, "sum", na.rm = TRUE)[[1]])) +
  map_theme + coord_sf(crs = 27700)

ggsave("figures/habitats/oyster_potential_step6_temperature.png", p_oy_temp,
       width = 14, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/habitats/oyster_potential_step6_temperature.png\n")


## ---- STEP 7: DISSOLVED OXYGEN FILTER ----

cat("\n-- Step 7: Dissolved oxygen filter (> 3.5 mg/l) --\n")

o2_file <- list.files("data/processed/physical/o2", pattern = "\\.nc$",
                       full.names = TRUE)
if (length(o2_file) > 0 && file.exists(o2_file[1])) {
  o2_stack <- rast(o2_file[1])
  cat(sprintf("  DO layers: %d\n", nlyr(o2_stack)))

  o2_mean_wgs <- app(o2_stack, fun = "mean", na.rm = TRUE)
  aoi_o2 <- vect(st_transform(scamp_aoi_sf, crs(o2_mean_wgs)))
  o2_crop <- crop(o2_mean_wgs, aoi_o2, mask = TRUE)
  o2_r <- project(o2_crop, oy_bathy, method = "bilinear")
  names(o2_r) <- "do_mean"

  cat(sprintf("  DO range: %.2f to %.2f mg/l\n",
              global(o2_r, "min", na.rm = TRUE)[[1]],
              global(o2_r, "max", na.rm = TRUE)[[1]]))

  oy_do_mask <- ifel(is.na(o2_r), 1, ifel(o2_r > 3.5, 1, 0))
  names(oy_do_mask) <- "dissolved_oxygen"

  cat(sprintf("  DO filter (> 3.5 mg/l): %d pixels pass (%.1f%%)\n",
              global(oy_do_mask, "sum", na.rm = TRUE)[[1]],
              global(oy_do_mask, "mean", na.rm = TRUE)[[1]] * 100))

  p_oy_do <- ggplot() +
    geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
    geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
            linewidth = 0.3) +
    tidyterra::geom_spatraster(data = o2_r) +
    scale_fill_viridis_c(name = "mg/l", na.value = NA, option = "viridis") +
    geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
    labs(title = "Oyster Potential — Step 7: Dissolved Oxygen",
         subtitle = sprintf("Range: %.2f–%.2f mg/l | Threshold: > 3.5 mg/l | Pass: %d px",
                            global(o2_r, "min", na.rm = TRUE)[[1]],
                            global(o2_r, "max", na.rm = TRUE)[[1]],
                            global(oy_do_mask, "sum", na.rm = TRUE)[[1]])) +
    map_theme + coord_sf(crs = 27700)

  ggsave("figures/habitats/oyster_potential_step7_do.png", p_oy_do,
         width = 14, height = 10, units = "in", dpi = 300)
  cat("  Saved: figures/habitats/oyster_potential_step7_do.png\n")

} else {
  cat("  WARNING: Dissolved oxygen data not available — setting DO mask to all 1s\n")
  oy_do_mask <- oy_bathy * 0 + 1
  names(oy_do_mask) <- "dissolved_oxygen"
}


## ---- STEP 8: CHLOROPHYLL / FOOD FILTER ----

cat("\n-- Step 8: Chlorophyll filter (growth > 0.5, reproductive > 1.68 µg/l) --\n")

chl_file <- list.files("data/processed/physical/chl", pattern = "\\.nc$",
                        full.names = TRUE)
if (length(chl_file) > 0 && file.exists(chl_file[1])) {
  chl_stack <- rast(chl_file[1])
  cat(sprintf("  Chl layers: %d\n", nlyr(chl_stack)))

  chl_mean_wgs <- app(chl_stack, fun = "mean", na.rm = TRUE)
  aoi_chl <- vect(st_transform(scamp_aoi_sf, crs(chl_mean_wgs)))
  chl_crop <- crop(chl_mean_wgs, aoi_chl, mask = TRUE)
  chl_r <- project(chl_crop, oy_bathy, method = "bilinear")
  names(chl_r) <- "chl_mean"

  cat(sprintf("  Chl range: %.3f to %.3f µg/l\n",
              global(chl_r, "min", na.rm = TRUE)[[1]],
              global(chl_r, "max", na.rm = TRUE)[[1]]))

  ## Growth threshold (for combined mask)
  oy_chl_mask <- ifel(is.na(chl_r), 1, ifel(chl_r > 0.5, 1, 0))
  names(oy_chl_mask) <- "chlorophyll"

  ## Reproductive threshold (quality layer)
  oy_chl_repro <- ifel(is.na(chl_r), NA, ifel(chl_r > 1.68, 1, 0))
  names(oy_chl_repro) <- "chl_reproductive"

  cat(sprintf("  Growth filter (> 0.5 µg/l): %d pixels pass (%.1f%%)\n",
              global(oy_chl_mask, "sum", na.rm = TRUE)[[1]],
              global(oy_chl_mask, "mean", na.rm = TRUE)[[1]] * 100))
  cat(sprintf("  Reproductive threshold (> 1.68 µg/l): %d pixels\n",
              global(oy_chl_repro, "sum", na.rm = TRUE)[[1]]))

  p_oy_chl <- ggplot() +
    geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
    geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
            linewidth = 0.3) +
    tidyterra::geom_spatraster(data = chl_r) +
    scale_fill_viridis_c(name = "µg/l", na.value = NA, option = "mako") +
    geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
    labs(title = "Oyster Potential — Step 8: Chlorophyll",
         subtitle = sprintf("Range: %.3f–%.3f µg/l | Growth (>0.5): %d px | Repro (>1.68): %d px",
                            global(chl_r, "min", na.rm = TRUE)[[1]],
                            global(chl_r, "max", na.rm = TRUE)[[1]],
                            global(oy_chl_mask, "sum", na.rm = TRUE)[[1]],
                            global(oy_chl_repro, "sum", na.rm = TRUE)[[1]])) +
    map_theme + coord_sf(crs = 27700)

  ggsave("figures/habitats/oyster_potential_step8_chl.png", p_oy_chl,
         width = 14, height = 10, units = "in", dpi = 300)
  cat("  Saved: figures/habitats/oyster_potential_step8_chl.png\n")

} else {
  cat("  WARNING: Chlorophyll data not available — setting chl mask to all 1s\n")
  oy_chl_mask <- oy_bathy * 0 + 1
  names(oy_chl_mask) <- "chlorophyll"
  oy_chl_repro <- oy_bathy * NA
  names(oy_chl_repro) <- "chl_reproductive"
}


## ---- STEP 8b: BED SHEAR STRESS FILTER ----

cat("\n-- Step 8b: Bed shear stress filter (KE currents) --\n")

oy_ke <- resample(ke_bng, oy_bathy, method = "bilinear")
names(oy_ke) <- "ke_shear"

cat(sprintf("  KE shear range: %.4f to %.4f\n",
            global(oy_ke, "min", na.rm = TRUE)[[1]],
            global(oy_ke, "max", na.rm = TRUE)[[1]]))

## Investigate two thresholds
ke_mask_130 <- ifel(is.na(oy_ke), 1, ifel(oy_ke < 130, 1, 0))
ke_mask_10  <- ifel(is.na(oy_ke), 1, ifel(oy_ke < 10, 1, 0))

n_pass_130 <- global(ke_mask_130, "sum", na.rm = TRUE)[[1]]
n_pass_10  <- global(ke_mask_10, "sum", na.rm = TRUE)[[1]]
n_total    <- global(!is.na(oy_ke) | is.na(oy_ke), "sum", na.rm = TRUE)[[1]]

cat(sprintf("  Threshold < 130: %d pixels pass (%.1f%%)\n", n_pass_130,
            n_pass_130 / n_total * 100))
cat(sprintf("  Threshold < 10:  %d pixels pass (%.1f%%)\n", n_pass_10,
            n_pass_10 / n_total * 100))

## Side-by-side diagnostic: both thresholds
p_ke_130 <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
          linewidth = 0.3) +
  tidyterra::geom_spatraster(data = as.factor(ke_mask_130)) +
  scale_fill_manual(values = c("0" = "transparent", "1" = "#5e3c99"),
                    labels = c("Excluded", "< 130"), name = "KE shear",
                    na.translate = FALSE) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  labs(title = "Bed Shear Stress — threshold < 130",
       subtitle = sprintf("%d pixels pass", n_pass_130)) +
  map_theme + coord_sf(crs = 27700)

p_ke_10 <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
          linewidth = 0.3) +
  tidyterra::geom_spatraster(data = as.factor(ke_mask_10)) +
  scale_fill_manual(values = c("0" = "transparent", "1" = "#e66101"),
                    labels = c("Excluded", "< 10"), name = "KE shear",
                    na.translate = FALSE) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  labs(title = "Bed Shear Stress — threshold < 10",
       subtitle = sprintf("%d pixels pass", n_pass_10)) +
  map_theme + coord_sf(crs = 27700)

p_ke_compare <- p_ke_130 + p_ke_10 +
  patchwork::plot_annotation(
    title = "Oyster Potential — Step 8b: Bed Shear Stress Comparison",
    subtitle = sprintf("KE range: %.2f–%.2f | <130: %d px | <10: %d px",
                       global(oy_ke, "min", na.rm = TRUE)[[1]],
                       global(oy_ke, "max", na.rm = TRUE)[[1]],
                       n_pass_130, n_pass_10))

ggsave("figures/habitats/oyster_potential_step8b_shear_compare.png", p_ke_compare,
       width = 20, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/habitats/oyster_potential_step8b_shear_compare.png\n")

## Also save the raw KE values for reference
p_ke_raw <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
          linewidth = 0.3) +
  tidyterra::geom_spatraster(data = oy_ke) +
  scale_fill_viridis_c(name = "KE", na.value = NA, option = "turbo") +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  labs(title = "Bed Shear Stress — KE currents (raw values)",
       subtitle = sprintf("Range: %.2f–%.2f",
                          global(oy_ke, "min", na.rm = TRUE)[[1]],
                          global(oy_ke, "max", na.rm = TRUE)[[1]])) +
  map_theme + coord_sf(crs = 27700)

ggsave("figures/habitats/oyster_potential_step8b_shear_raw.png", p_ke_raw,
       width = 14, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/habitats/oyster_potential_step8b_shear_raw.png\n")

## Use < 10 as default mask — review comparison plot and adjust if needed
oy_shear_mask <- ke_mask_130
names(oy_shear_mask) <- "bed_shear"
cat(sprintf("  Active threshold: < 130 | %d pixels pass\n", n_pass_10))


## ---- STEP 9: COMBINE ALL MASKS ----

cat("\n-- Step 9: Combine all masks --\n")

oy_potential <- oy_depth_mask * oy_substrate_mask * oy_sal_mask *
                oy_current_mask * oy_spm_mask * oy_temp_mask *
                oy_do_mask * oy_chl_mask * oy_shear_mask
names(oy_potential) <- "potential"

cat(sprintf("  Step 1 (depth):      %d pixels\n", global(oy_depth_mask, "sum", na.rm = TRUE)[[1]]))
cat(sprintf("  Step 2 (substrate):  %d pixels\n", global(oy_substrate_mask, "sum", na.rm = TRUE)[[1]]))
cat(sprintf("  Step 3 (salinity):   %d pixels\n", global(oy_sal_mask, "sum", na.rm = TRUE)[[1]]))
cat(sprintf("  Step 4 (current):    %d pixels\n", global(oy_current_mask, "sum", na.rm = TRUE)[[1]]))
cat(sprintf("  Step 5 (SPM):        %d pixels\n", global(oy_spm_mask, "sum", na.rm = TRUE)[[1]]))
cat(sprintf("  Step 6 (temp):       %d pixels\n", global(oy_temp_mask, "sum", na.rm = TRUE)[[1]]))
cat(sprintf("  Step 7 (DO):         %d pixels\n", global(oy_do_mask, "sum", na.rm = TRUE)[[1]]))
cat(sprintf("  Step 8 (chl):        %d pixels\n", global(oy_chl_mask, "sum", na.rm = TRUE)[[1]]))
cat(sprintf("  Step 8b (shear):     %d pixels\n", global(oy_shear_mask, "sum", na.rm = TRUE)[[1]]))
cat(sprintf("  COMBINED:            %d pixels pass all filters (%.1f%%)\n",
            global(oy_potential, "sum", na.rm = TRUE)[[1]],
            global(oy_potential, "mean", na.rm = TRUE)[[1]] * 100))

p_oy_combined <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
          linewidth = 0.3) +
  tidyterra::geom_spatraster(data = as.factor(oy_potential)) +
  scale_fill_manual(values = c("0" = "transparent", "1" = "#d6604d"),
                    labels = c("Excluded", "All filters pass"), name = "Suitability",
                    na.translate = FALSE) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  labs(title = "Oyster Potential — Step 9: Combined Mask",
       subtitle = sprintf("All 8 filters applied | %d pixels pass (%.1f%%)",
                          global(oy_potential, "sum", na.rm = TRUE)[[1]],
                          global(oy_potential, "mean", na.rm = TRUE)[[1]] * 100)) +
  map_theme + coord_sf(crs = 27700)

ggsave("figures/habitats/oyster_potential_step9_combined.png", p_oy_combined,
       width = 14, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/habitats/oyster_potential_step9_combined.png\n")

writeRaster(oy_potential, "outputs/habitats/oyster_potential_raster.tif",
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))
cat("  Saved: outputs/habitats/oyster_potential_raster.tif\n")


## ---- STEP 10: VECTORISE AND SUMMARISE ----

cat("\n-- Step 10: Vectorise and summarise --\n")

oy_poly <- as.polygons(oy_potential) %>% st_as_sf()
oy_poly <- oy_poly[oy_poly$potential == 1, ]
oy_poly <- st_transform(oy_poly, 27700) %>% st_make_valid()
oy_poly <- st_cast(oy_poly, "POLYGON")
oy_poly$area_ha <- as.numeric(st_area(oy_poly)) / 10000

cat(sprintf("  Potential oyster polygons: %d\n", nrow(oy_poly)))
cat(sprintf("  Total potential area: %.1f ha\n", sum(oy_poly$area_ha)))

st_write(oy_poly, "outputs/habitats/oyster_potential_solway.gpkg",
         delete_dsn = TRUE, quiet = TRUE)
cat("  Saved: outputs/habitats/oyster_potential_solway.gpkg\n")


## ---- STEP 11: VALIDATION ----

cat("\n-- Step 11: Validation against known records --\n")

## Load known and historical oyster records
oy_known_path <- "data/processed/species/oyster_known_extent.gpkg"
oy_hist_path  <- "data/processed/species/oyster_historical_extent.gpkg"
oy_pts_path   <- "data/processed/species/oyster_observation_points.gpkg"

has_known <- file.exists(oy_known_path)
has_hist  <- file.exists(oy_hist_path)
has_pts   <- file.exists(oy_pts_path)

if (has_known) {
  oy_known_sf <- st_read(oy_known_path, quiet = TRUE) %>% st_make_valid()
  cat(sprintf("  Known reef extent: %d features\n", nrow(oy_known_sf)))
} else {
  oy_known_sf <- st_sf(geometry = st_sfc(crs = 27700))
  cat("  Known reef extent: not found\n")
}

if (has_hist) {
  oy_hist_sf <- st_read(oy_hist_path, quiet = TRUE) %>% st_make_valid()
  cat(sprintf("  Historical reef extent: %d features\n", nrow(oy_hist_sf)))
} else {
  oy_hist_sf <- st_sf(geometry = st_sfc(crs = 27700))
  cat("  Historical reef extent: not found\n")
}

if (has_pts) {
  oy_pts_sf <- st_read(oy_pts_path, quiet = TRUE) %>% st_make_valid()
  cat(sprintf("  Observation points: %d\n", nrow(oy_pts_sf)))
} else {
  oy_pts_sf <- st_sf(geometry = st_sfc(crs = 27700))
  cat("  Observation points: not found\n")
}

## Validation plot: modelled potential + known/historical records
p_oy_validate <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
          linewidth = 0.3) +
  geom_sf(data = scamp_coastline, colour = "grey30", linewidth = 0.3) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  geom_sf(data = oy_poly, fill = "#d6604d", colour = NA, alpha = 0.5) +
  {if (nrow(oy_known_sf) > 0)
    geom_sf(data = oy_known_sf, colour = "#e7298a", fill = "#e7298a",
            alpha = 0.6, linewidth = 0.5)} +
  {if (nrow(oy_hist_sf) > 0)
    geom_sf(data = oy_hist_sf, colour = "#a6761d", fill = "#a6761d",
            alpha = 0.4, linewidth = 0.5)} +
  {if (nrow(oy_pts_sf) > 0)
    geom_sf(data = oy_pts_sf, colour = "black", shape = 4,
            size = 2, stroke = 1, alpha = 0.8)} +
  annotation_scale(location = "bl", width_hint = 0.2) +
  annotation_north_arrow(location = "tl", style = north_arrow_minimal(),
                         height = unit(1, "cm"), width = unit(1, "cm")) +
  labs(title = "Native Oyster Potential — Validation",
       subtitle = sprintf(
         "Modelled: %.1f ha | Known reefs: %d | Historical: %d | Obs points: %d",
         sum(oy_poly$area_ha), nrow(oy_known_sf), nrow(oy_hist_sf), nrow(oy_pts_sf))) +
  map_theme + coord_sf(crs = 27700)

ggsave("figures/habitats/oyster_potential_validation.png", p_oy_validate,
       width = 14, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/habitats/oyster_potential_validation.png\n")


## ---- STEP 12: FUTURE RISK SCREENING (placeholder) ----

cat("\n-- Step 12: Future risk screening --\n")
cat("  [placeholder — awaiting Bio-ORACLE present/2050 layers]\n")
cat("  Required: SST, salinity, and pH for present and 2050 scenarios\n")
cat("  Thresholds: SST > 28°C, salinity < 20 PSU, pH < 6.9\n")


## Diagnostic summary
cat("\n-- Native Oyster Potential: diagnostic summary --\n")
cat(sprintf("  Step 1  (depth 0 to -80m):       %d pixels\n", global(oy_depth_mask, "sum", na.rm = TRUE)[[1]]))
cat(sprintf("  Step 2  (substrate BGS Folk):      %d pixels\n", global(oy_substrate_mask, "sum", na.rm = TRUE)[[1]]))
cat(sprintf("  Step 3  (sal min-monthly >= 20):   %d pixels\n", global(oy_sal_mask, "sum", na.rm = TRUE)[[1]]))
cat(sprintf("  Step 4  (current <= 0.8 m/s):     %d pixels\n", global(oy_current_mask, "sum", na.rm = TRUE)[[1]]))
cat(sprintf("  Step 5  (SPM < 60 mg/l):          %d pixels\n", global(oy_spm_mask, "sum", na.rm = TRUE)[[1]]))
cat(sprintf("  Step 6  (temp 15-28/> 5°C):       %d pixels\n", global(oy_temp_mask, "sum", na.rm = TRUE)[[1]]))
cat(sprintf("  Step 7  (DO > 3.5 mg/l):          %d pixels\n", global(oy_do_mask, "sum", na.rm = TRUE)[[1]]))
cat(sprintf("  Step 8  (chl > 0.5 µg/l):         %d pixels\n", global(oy_chl_mask, "sum", na.rm = TRUE)[[1]]))
cat(sprintf("  Step 8b (bed shear < 10):         %d pixels\n", global(oy_shear_mask, "sum", na.rm = TRUE)[[1]]))
cat(sprintf("  Step 9  (combined):               %d pixels\n", global(oy_potential, "sum", na.rm = TRUE)[[1]]))
cat(sprintf("  Step 10 (vectorised):             %d polygons | %.1f ha\n",
            nrow(oy_poly), sum(oy_poly$area_ha)))

cat("\n=== Native Oyster Potential Model complete ===\n")
