# 07_plots.R — Report-quality figures for SCAMP GIS
# Organised by habitat section to match the report structure.
# Saves static figures to figures/report draft/ and an interactive leaflet HTML.

# ============================================================
# 0. SETUP
# ============================================================
library(sf)
library(terra)
library(ggplot2)
library(tidyterra)
library(patchwork)
library(ggspatial)
library(dplyr)
library(rosm)
library(leaflet)
library(htmlwidgets)

cat("\n========== 07_plots.R ==========\n")
cat("Loading base layers...\n")

# --- Boundaries ---
scamp_aoi <- st_read("data/processed/boundaries/scamp_aoi.gpkg", quiet = TRUE) %>%
  st_transform(27700)

smr_raw <- st_read(
  "data/raw/boundaries/administrative_units_scottish_marine_regions/administrative_units_scottish_marine_regionsPolygon.shp",
  quiet = TRUE
) %>% st_transform(27700)
smr <- smr_raw[grepl("Solway", smr_raw$objnam, ignore.case = TRUE), ]

dg_raw <- st_read(
  "data/raw/boundaries/Local_Authority_Boundaries_-_Scotland/pub_las.shp",
  quiet = TRUE
) %>% st_transform(27700)
dg <- dg_raw[grepl("dumfries", dg_raw$local_auth, ignore.case = TRUE), ]

coastline <- st_read("data/processed/boundaries/scamp_coastline.gpkg", quiet = TRUE) %>%
  st_transform(27700)

harbour <- st_read("data/processed/boundaries/harbour_limits_scamp.gpkg", quiet = TRUE) %>%
  st_transform(27700)

# --- Bathymetry for ocean background (used in study area variant) ---
bathy <- rast("outputs/physical/gebco_bathy_solway.tif")
bathy[bathy > 0] <- NA

# --- Designated sites ---
desig_files <- list(
  SAC    = "data/processed/boundaries/sac_scamp.gpkg",
  SPA    = "data/processed/boundaries/spa_scamp.gpkg",
  SSSI   = "data/processed/boundaries/sssi_scamp.gpkg",
  Ramsar = "data/processed/boundaries/ramsar_scamp.gpkg",
  MPA    = "data/processed/boundaries/mpas_scamp.gpkg",
  NNR    = "data/processed/boundaries/nnr_scamp.gpkg",
  LNR    = "data/processed/boundaries/lnr_scamp.gpkg",
  IBA    = "data/processed/boundaries/ibas_scamp.gpkg",
  RSPB   = "data/processed/boundaries/rspb_reserves_scamp.gpkg"
)
desig_layers <- list()
for (nm in names(desig_files)) {
  fp <- desig_files[[nm]]
  if (file.exists(fp)) {
    d <- st_read(fp, quiet = TRUE) %>% st_transform(27700)
    d$designation <- nm
    desig_layers[[nm]] <- d
  }
}
cat(sprintf("  Loaded %d designation types\n", length(desig_layers)))

# --- Output directory ---
out_dir <- "figures/report draft"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# THEME & HELPER FUNCTIONS
# ============================================================

theme_scamp <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.major = element_line(colour = "grey85", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      axis.text      = element_text(size = 7),
      axis.title     = element_blank(),
      plot.title     = element_text(face = "bold", size = 11),
      plot.subtitle  = element_text(size = 9, colour = "grey30"),
      legend.position  = "bottom",
      legend.key.size  = unit(0.4, "cm"),
      legend.text      = element_text(size = 7),
      legend.title     = element_text(size = 8, face = "bold")
    )
}

# Base map: OSM tiles for land context + place names; boundaries overlaid
base_map <- function(title = NULL, subtitle = NULL, tile_type = "osm", zoom = 9) {
  ggplot() +
    annotation_map_tile(type = tile_type, zoom = zoom, progress = "none",
                        quiet = TRUE) +
    geom_sf(data = dg, fill = NA, colour = "grey40", linewidth = 0.4) +
    geom_sf(data = smr, fill = NA, colour = "steelblue", linewidth = 0.5,
            linetype = "dotted") +
    geom_sf(data = scamp_aoi, fill = NA, colour = "black", linewidth = 0.6,
            linetype = "dashed") +
    annotation_scale(location = "bl", width_hint = 0.2,
                     style = "bar", text_cex = 0.7, line_width = 0.5,
                     bar_cols = c("black", "white")) +
    annotation_north_arrow(location = "tr", which_north = "true",
                           height = unit(1, "cm"), width = unit(0.8, "cm"),
                           style = north_arrow_orienteering()) +
    coord_sf(crs = 27700) +
    labs(title = title, subtitle = subtitle) +
    theme_scamp()
}

save_plot <- function(p, filename, width = 10, height = 8, dpi = 300) {
  path <- file.path(out_dir, filename)
  ggsave(path, p, width = width, height = height, dpi = dpi, bg = "white")
  cat(sprintf("  Saved: %s\n", path))
}

# Load cluster boxes saved by 03_habitats.R
load_clusters <- function(path) {
  if (!file.exists(path)) {
    warning(sprintf("Cluster file not found: %s — run 03_habitats.R first", path))
    return(st_sf(name = character(), geometry = st_sfc(crs = 27700)))
  }
  st_read(path, quiet = TRUE)
}


# ============================================================
# 1. OVERVIEW
# ============================================================
cat("\n===== 1. OVERVIEW =====\n")

# 1a. Study area — OSM version (place names, roads, land context)
p_study_osm <- base_map("SCAMP Study Area",
                         "Solway Firth, Dumfries & Galloway")
save_plot(p_study_osm, "01a_study_area_osm.png")

# 1b. Study area — bathymetry version (ocean depth detail)
p_study_bathy <- ggplot() +
  geom_spatraster(data = bathy) +
  scale_fill_gradientn(
    colours = c("#08306b", "#2171b5", "#6baed6", "#c6dbef", "#deebf7"),
    na.value = NA, name = "Depth (m)"
  ) +
  annotation_map_tile(type = "osm", zoom = 9, progress = "none", quiet = TRUE,
                      alpha = 0.7) +
  geom_sf(data = dg, fill = NA, colour = "grey40", linewidth = 0.4) +
  geom_sf(data = smr, fill = NA, colour = "steelblue", linewidth = 0.5,
          linetype = "dotted") +
  geom_sf(data = scamp_aoi, fill = NA, colour = "black", linewidth = 0.6,
          linetype = "dashed") +
  annotation_scale(location = "bl", width_hint = 0.2,
                   style = "bar", text_cex = 0.7, bar_cols = c("black", "white")) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         height = unit(1, "cm"), width = unit(0.8, "cm"),
                         style = north_arrow_orienteering()) +
  coord_sf(crs = 27700) +
  labs(title = "SCAMP Study Area",
       subtitle = "Solway Firth, Dumfries & Galloway — with bathymetry") +
  theme_scamp()
save_plot(p_study_bathy, "01b_study_area_bathymetry.png")

# 1c. Designated sites
desig_colours <- c(
  SAC = "#1b9e77", SPA = "#d95f02", SSSI = "#7570b3", Ramsar = "#e7298a",
  MPA = "#66a61e", NNR = "#e6ab02", LNR = "#a6761d", IBA = "#666666"
)
desig_all <- bind_rows(desig_layers) %>% st_make_valid()

p_desig <- base_map("Designated Sites", "SCAMP AOI") +
  geom_sf(data = desig_all, aes(fill = designation), colour = NA, alpha = 0.35) +
  geom_sf(data = harbour, fill = NA, colour = "red", linewidth = 0.5,
          linetype = "longdash") +
  scale_fill_manual(values = desig_colours, name = "Designation") +
  guides(fill = guide_legend(nrow = 1, order = 2))

save_plot(p_desig, "02_designated_sites.png")


# ============================================================
# 2. SEAGRASS
# ============================================================
cat("\n===== 2. SEAGRASS =====\n")

sg_known <- st_read("data/processed/seagrass_polygons_harmonised.gpkg", quiet = TRUE) %>%
  st_transform(27700)
sg_potential <- st_read("outputs/habitats/seagrass_potential_solway.gpkg", quiet = TRUE) %>%
  st_transform(27700)

ind_all <- st_read("outputs/species/indicator_records.gpkg", quiet = TRUE) %>%
  st_transform(st_crs(smr))
sg_indicators <- ind_all[which(ind_all$habitat == "Seagrass"), ]
sg_indicators <- sg_indicators[lengths(st_intersects(sg_indicators, smr)) > 0, ]
cat(sprintf("  Seagrass indicators (within SMR): %d\n", nrow(sg_indicators)))

sg_boxes <- load_clusters("data/processed/habitats/seagrass_cluster_boxes.gpkg")
cat(sprintf("  Seagrass clusters: %d\n", nrow(sg_boxes)))

# 2a. Combined: known + potential + indicators + harbour limits + cluster boxes
p_sg_combined <- base_map("Seagrass",
                           "Known extent, modelled potential & Zostera indicators") +
  geom_sf(data = sg_potential, aes(fill = "Modelled potential"), colour = NA,
          alpha = 0.4) +
  geom_sf(data = sg_known, aes(fill = "Known extent"), colour = "#1b7837",
          alpha = 0.6) +
  geom_sf(data = harbour, fill = NA, colour = "red", linewidth = 0.4,
          linetype = "longdash", show.legend = FALSE) +
  {if (nrow(sg_indicators) > 0)
    geom_sf(data = sg_indicators, aes(colour = "Zostera indicators"),
            size = 1.2, alpha = 0.7, shape = 17)} +
  geom_sf(data = sg_boxes, fill = NA, colour = "red", linewidth = 0.6) +
  geom_sf_text(data = sg_boxes, aes(label = name), colour = "red",
               size = 2.5, fontface = "bold") +
  scale_fill_manual(values = c("Known extent" = "#1b7837",
                                "Modelled potential" = "#a6d854"),
                    name = NULL) +
  scale_colour_manual(values = c("Zostera indicators" = "#4d004b"), name = NULL) +
  guides(fill = guide_legend(order = 1), colour = guide_legend(order = 2))
save_plot(p_sg_combined, "03_seagrass_combined.png")

# 2b. Individual: known extent
p_sg_known <- base_map("Seagrass — Known Extent") +
  geom_sf(data = sg_known, fill = "#1b7837", colour = "#1b7837", alpha = 0.7)
save_plot(p_sg_known, "04_seagrass_known.png")

# 2c. Individual: modelled potential
p_sg_pot <- base_map("Seagrass — Modelled Potential") +
  geom_sf(data = sg_potential, fill = "#a6d854", colour = NA, alpha = 0.5) +
  geom_sf(data = harbour, fill = NA, colour = "red", linewidth = 0.4,
          linetype = "longdash")
save_plot(p_sg_pot, "05_seagrass_potential.png")

# 2d. Individual: Zostera indicators
p_sg_ind <- base_map("Seagrass — Zostera Indicator Species",
                      "Filtered to SMR boundary") +
  {if (nrow(sg_indicators) > 0)
    geom_sf(data = sg_indicators, aes(colour = latin_name), size = 1.5, alpha = 0.7,
            shape = 17)} +
  scale_colour_brewer(palette = "Dark2", name = "Species") +
  guides(colour = guide_legend(nrow = 2))
save_plot(p_sg_ind, "06_seagrass_indicators.png")

# 2e. Seagrass cluster sub-plots (zoomed to each cluster box)
for (i in seq_len(nrow(sg_boxes))) {
  bb <- st_bbox(sg_boxes[i, ])
  p_cl <- base_map(sprintf("Seagrass — %s", sg_boxes$name[i]), zoom = 13) +
    geom_sf(data = sg_potential, aes(fill = "Modelled potential"), colour = NA,
            alpha = 0.4) +
    geom_sf(data = sg_known, aes(fill = "Known extent"), colour = "#1b7837",
            alpha = 0.6) +
    {if (nrow(sg_indicators) > 0)
      geom_sf(data = sg_indicators, aes(colour = "Zostera indicators"),
              size = 2.5, alpha = 0.7, shape = 17)} +
    scale_fill_manual(values = c("Known extent" = "#1b7837",
                                  "Modelled potential" = "#a6d854"),
                      name = NULL) +
    scale_colour_manual(values = c("Zostera indicators" = "#4d004b"), name = NULL) +
    coord_sf(xlim = c(bb["xmin"], bb["xmax"]),
             ylim = c(bb["ymin"], bb["ymax"]), crs = 27700)
  save_plot(p_cl, sprintf("07_seagrass_cluster_%02d.png", i), width = 8, height = 8)
}

# 2f. Seagrass suitability patchwork (annex) — individual physical parameters
cat("  Building seagrass suitability annex panels...\n")

suit_rasters <- list(
  "Bathymetry (m)"       = "outputs/physical/gebco_bathy_solway.tif",
  "SPM (mg/l)"           = "outputs/physical/spm_mean_solway.tif",
  "Salinity P90 (PSU)"   = "outputs/physical/salinity_p90_mean_solway.tif",
  "Current energy"       = "outputs/physical/ke_currents_solway.tif",
  "Suitability (combined)" = "outputs/habitats/seagrass_potential_raster.tif"
)

sg_panels <- list()
for (nm in names(suit_rasters)) {
  fp <- suit_rasters[[nm]]
  if (!file.exists(fp)) next
  r <- rast(fp)
  p <- ggplot() +
    annotation_map_tile(type = "osm", zoom = 9, progress = "none", quiet = TRUE) +
    geom_spatraster(data = r) +
    geom_sf(data = dg, fill = NA, colour = "grey40", linewidth = 0.2) +
    geom_sf(data = scamp_aoi, fill = NA, colour = "black", linewidth = 0.3,
            linetype = "dashed") +
    scale_fill_viridis_c(na.value = NA, name = nm) +
    coord_sf(crs = 27700) +
    labs(title = nm) +
    theme_scamp(base_size = 8) +
    theme(legend.key.width = unit(0.6, "cm"))
  sg_panels <- c(sg_panels, list(p))
}

if (length(sg_panels) > 0) {
  ncols <- min(length(sg_panels), 3)
  p_sg_annex <- wrap_plots(sg_panels, ncol = ncols) +
    plot_annotation(
      title = "Seagrass Physical Suitability — Input Parameters",
      theme = theme(plot.title = element_text(face = "bold", size = 14))
    )
  save_plot(p_sg_annex, "08_seagrass_suitability_annex.png",
            width = 7 * ncols, height = 7 * ceiling(length(sg_panels) / ncols))
}


# ============================================================
# 3. SALTMARSH
# ============================================================
cat("\n===== 3. SALTMARSH =====\n")

sm_known <- st_read("data/processed/species/saltmarsh_harmonised.gpkg", quiet = TRUE) %>%
  st_transform(27700)
sm_ind <- st_read("data/processed/species/saltmarsh_indicator_points.gpkg", quiet = TRUE) %>%
  st_transform(27700)

sm_boxes <- load_clusters("data/processed/habitats/saltmarsh_cluster_boxes.gpkg")
cat(sprintf("  Saltmarsh clusters: %d\n", nrow(sm_boxes)))

# 3a. Combined: known extent + indicator species + cluster boxes
p_sm_combined <- base_map("Saltmarsh", "Known extent & indicator species") +
  geom_sf(data = sm_known, aes(fill = "Known extent"), colour = "#238443",
          alpha = 0.6) +
  {if (nrow(sm_ind) > 0)
    geom_sf(data = sm_ind, aes(colour = "Indicator species"),
            size = 1.2, alpha = 0.6, shape = 16)} +
  geom_sf(data = sm_boxes, fill = NA, colour = "red", linewidth = 0.6) +
  geom_sf_text(data = sm_boxes, aes(label = name), colour = "red",
               size = 2.5, fontface = "bold") +
  scale_fill_manual(values = c("Known extent" = "#78c679"), name = NULL) +
  scale_colour_manual(values = c("Indicator species" = "#004529"), name = NULL)
save_plot(p_sm_combined, "09_saltmarsh_combined.png")

# 3b. Individual: known extent
p_sm_known <- base_map("Saltmarsh — Known Extent") +
  geom_sf(data = sm_known, fill = "#78c679", colour = "#238443", alpha = 0.7)
save_plot(p_sm_known, "10_saltmarsh_known.png")

# 3c. Individual: indicator species
p_sm_ind_plot <- base_map("Saltmarsh — Indicator Species") +
  {if (nrow(sm_ind) > 0)
    geom_sf(data = sm_ind, aes(colour = latin_name), size = 1.5, alpha = 0.7,
            shape = 16)} +
  scale_colour_brewer(palette = "Set2", name = "Species") +
  guides(colour = guide_legend(nrow = 2))
save_plot(p_sm_ind_plot, "11_saltmarsh_indicators.png")

# 3d. Cluster sub-plots (zoomed versions of combined)
for (i in seq_len(nrow(sm_boxes))) {
  bb <- st_bbox(sm_boxes[i, ])
  p_cl <- base_map(sprintf("Saltmarsh — %s", sm_boxes$name[i]), zoom = 14) +
    geom_sf(data = sm_known, aes(fill = "Known extent"), colour = "#238443",
            alpha = 0.6) +
    {if (nrow(sm_ind) > 0)
      geom_sf(data = sm_ind, aes(colour = "Indicator species"),
              size = 2.5, alpha = 0.6, shape = 16)} +
    scale_fill_manual(values = c("Known extent" = "#78c679"), name = NULL) +
    scale_colour_manual(values = c("Indicator species" = "#004529"), name = NULL) +
    coord_sf(xlim = c(bb["xmin"], bb["xmax"]),
             ylim = c(bb["ymin"], bb["ymax"]), crs = 27700)
  save_plot(p_cl, sprintf("12_saltmarsh_cluster_%02d.png", i), width = 8, height = 8)
}


# ============================================================
# 4. NATIVE OYSTERS
# ============================================================
cat("\n===== 4. NATIVE OYSTERS =====\n")

oy_known <- st_read("data/processed/species/oyster_known_extent.gpkg", quiet = TRUE) %>%
  st_transform(27700)

oy_hist_path <- "data/processed/species/oyster_historical_locations.gpkg"
oy_hist <- if (file.exists(oy_hist_path)) {
  st_read(oy_hist_path, quiet = TRUE) %>% st_transform(27700)
} else {
  cat("  WARNING: oyster_historical_extent.gpkg not found — skipping\n")
  st_sf(geometry = st_sfc(crs = 27700))
}

oy_potential <- st_read("outputs/habitats/oyster_potential_solway.gpkg", quiet = TRUE) %>%
  st_transform(27700)

oy_indicators <- ind_all[which(ind_all$habitat == "Oyster"), ]
oy_indicators <- oy_indicators[lengths(st_intersects(oy_indicators, smr)) > 0, ]
cat(sprintf("  Oyster indicators (within SMR): %d\n", nrow(oy_indicators)))

oy_boxes <- load_clusters("data/processed/habitats/oyster_cluster_boxes.gpkg")
cat(sprintf("  Oyster clusters: %d\n", nrow(oy_boxes)))

# 4a. Combined
p_oy_combined <- base_map("Native Oysters",
                           "Known, historical, modelled potential & indicator species") +
  geom_sf(data = oy_potential, aes(fill = "Modelled potential"), colour = NA,
          alpha = 0.4) +
  {if (nrow(oy_hist) > 0)
    geom_sf(data = oy_hist, aes(colour = "Historical"), shape = 17,
            size = 2, alpha = 0.6)} +
  geom_sf(data = oy_known, aes(fill = "Known extent"), colour = "#54278f",
          alpha = 0.7) +
  {if (nrow(oy_indicators) > 0)
    geom_sf(data = oy_indicators, aes(colour = "Indicator species"),
            size = 1.2, alpha = 0.7, shape = 18)} +
  {if (nrow(oy_boxes) > 0)
    geom_sf(data = oy_boxes, fill = NA, colour = "red", linewidth = 0.6)} +
  {if (nrow(oy_boxes) > 0)
    geom_sf_text(data = oy_boxes, aes(label = name), colour = "red",
                 size = 2.5, fontface = "bold")} +
  scale_fill_manual(values = c("Known extent" = "#756bb1",
                                "Historical" = "#fc9272",
                                "Modelled potential" = "#bcbddc"),
                    name = NULL) +
  scale_colour_manual(values = c("Indicator species" = "#e6550d"), name = NULL)
save_plot(p_oy_combined, "13_oyster_combined.png")

# 4b–e. Individuals
p_oy_known <- base_map("Native Oysters — Known Extent") +
  geom_sf(data = oy_known, fill = "#756bb1", colour = "#54278f", alpha = 0.7)
save_plot(p_oy_known, "14_oyster_known.png")

if (nrow(oy_hist) > 0) {
  p_oy_hist <- base_map("Native Oysters — Historical Extent") +
    geom_sf(data = oy_hist, colour = "#de2d26", shape = 17, size = 2, alpha = 0.6)
  save_plot(p_oy_hist, "15_oyster_historical.png")
}

p_oy_pot <- base_map("Native Oysters — Modelled Potential") +
  geom_sf(data = oy_potential, fill = "#bcbddc", colour = NA, alpha = 0.5)
save_plot(p_oy_pot, "16_oyster_potential.png")

p_oy_ind <- base_map("Native Oysters — Indicator Species",
                      "Filtered to SMR boundary") +
  {if (nrow(oy_indicators) > 0)
    geom_sf(data = oy_indicators, aes(colour = latin_name), size = 1.5, alpha = 0.7,
            shape = 18)} +
  scale_colour_brewer(palette = "Set1", name = "Species") +
  guides(colour = guide_legend(nrow = 2, order = 2))
save_plot(p_oy_ind, "17_oyster_indicators.png")

# 4f. Cluster sub-plots (zoomed versions of combined)
for (i in seq_len(nrow(oy_boxes))) {
  bb <- st_bbox(oy_boxes[i, ])
  p_cl <- base_map(sprintf("Native Oysters — %s", oy_boxes$name[i]), zoom = 14) +
    geom_sf(data = oy_potential, aes(fill = "Modelled potential"), colour = NA,
            alpha = 0.4) +
    {if (nrow(oy_hist) > 0)
      geom_sf(data = oy_hist, aes(colour = "Historical"), shape = 17,
              size = 2, alpha = 0.6)} +
    geom_sf(data = oy_known, aes(fill = "Known extent"), colour = "#54278f",
            alpha = 0.7) +
    {if (nrow(oy_indicators) > 0)
      geom_sf(data = oy_indicators, aes(colour = "Indicator species"),
              size = 2.5, alpha = 0.7, shape = 18)} +
    scale_fill_manual(values = c("Known extent" = "#756bb1",
                                  "Historical" = "#fc9272",
                                  "Modelled potential" = "#bcbddc"),
                      name = NULL) +
    scale_colour_manual(values = c("Indicator species" = "#e6550d"), name = NULL) +
    coord_sf(xlim = c(bb["xmin"], bb["xmax"]),
             ylim = c(bb["ymin"], bb["ymax"]), crs = 27700)
  save_plot(p_cl, sprintf("18_oyster_cluster_%02d.png", i), width = 8, height = 8)
}

# 4g. Oyster suitability patchwork (annex) — individual physical parameters
cat("  Building oyster suitability annex panels...\n")

oy_suit_rasters <- list(
  "Bathymetry (m)"              = "outputs/physical/gebco_bathy_solway.tif",
  "SPM (mg/l)"                  = "outputs/physical/spm_mean_solway.tif",
  "Salinity min monthly (PSU)"  = "outputs/physical/sal_min_monthly_mean_solway.tif",
  "Current energy"              = "outputs/physical/ke_currents_solway.tif",
  "Suitability (combined)"      = "outputs/habitats/oyster_potential_raster.tif"
)

oy_panels <- list()
for (nm in names(oy_suit_rasters)) {
  fp <- oy_suit_rasters[[nm]]
  if (!file.exists(fp)) next
  r <- rast(fp)
  p <- ggplot() +
    annotation_map_tile(type = "osm", zoom = 9, progress = "none", quiet = TRUE) +
    geom_spatraster(data = r) +
    geom_sf(data = dg, fill = NA, colour = "grey40", linewidth = 0.2) +
    geom_sf(data = scamp_aoi, fill = NA, colour = "black", linewidth = 0.3,
            linetype = "dashed") +
    scale_fill_viridis_c(na.value = NA, name = nm) +
    coord_sf(crs = 27700) +
    labs(title = nm) +
    theme_scamp(base_size = 8) +
    theme(legend.key.width = unit(0.6, "cm"))
  oy_panels <- c(oy_panels, list(p))
}

if (length(oy_panels) > 0) {
  ncols <- min(length(oy_panels), 3)
  p_oy_annex <- wrap_plots(oy_panels, ncol = ncols) +
    plot_annotation(
      title = "Native Oyster Physical Suitability — Input Parameters",
      theme = theme(plot.title = element_text(face = "bold", size = 14))
    )
  save_plot(p_oy_annex, "19_oyster_suitability_annex.png",
            width = 7 * ncols, height = 7 * ceiling(length(oy_panels) / ncols))
}


# ============================================================
# 5. SAND DUNES
# ============================================================
cat("\n===== 5. SAND DUNES =====\n")

dune_known <- st_read("outputs/habitats/sand_dunes_habmos_solway.gpkg", quiet = TRUE) %>%
  st_transform(27700)
dune_restore <- st_read("outputs/habitats/sand_dune_restoration_potential.gpkg",
                        quiet = TRUE) %>%
  st_transform(27700)

dune_indicators <- ind_all[which(ind_all$habitat == "Sand dunes"), ]
cat(sprintf("  Sand dune indicators: %d\n", nrow(dune_indicators)))

# NOTE: Shoreline change transect data is computed in 03_habitats.R but not saved
# to disk. To include it here, save transects_sf from 03_habitats.R first.

dune_boxes <- load_clusters("data/processed/habitats/dune_cluster_boxes.gpkg")
cat(sprintf("  Sand dune clusters: %d\n", nrow(dune_boxes)))

# 5a. Combined: known + indicators + restoration potential + cluster boxes
p_dune_combined <- base_map("Sand Dunes",
                             "Known extent, indicator species & restoration potential") +
  geom_sf(data = dune_restore, aes(fill = factor(tier)), colour = NA, alpha = 0.4) +
  geom_sf(data = dune_known, fill = "#d9b38c", colour = "#8c6d46", alpha = 0.7) +
  {if (nrow(dune_indicators) > 0)
    geom_sf(data = dune_indicators, colour = "#7a0177", size = 1.2, alpha = 0.7,
            shape = 15)} +
  geom_sf(data = dune_boxes, fill = NA, colour = "red", linewidth = 0.6) +
  geom_sf_text(data = dune_boxes, aes(label = name), colour = "red",
               size = 2.5, fontface = "bold") +
  scale_fill_manual(values = c("1" = "#fdae61", "2" = "#fee090"),
                    name = "Restoration tier",
                    labels = c("1" = "Tier 1 (near dunes)", "2" = "Tier 2"))
save_plot(p_dune_combined, "20_sand_dunes_combined.png")

# 5b–d. Individuals
p_dune_known <- base_map("Sand Dunes — Known Extent (HabMoS)") +
  geom_sf(data = dune_known, fill = "#d9b38c", colour = "#8c6d46", alpha = 0.7)
save_plot(p_dune_known, "21_sand_dunes_known.png")

p_dune_ind <- base_map("Sand Dunes — Indicator Species") +
  {if (nrow(dune_indicators) > 0)
    geom_sf(data = dune_indicators, aes(colour = latin_name), size = 1.5, alpha = 0.7,
            shape = 15)} +
  scale_colour_brewer(palette = "Set2", name = "Species") +
  guides(colour = guide_legend(nrow = 2, order = 2))
save_plot(p_dune_ind, "22_sand_dunes_indicators.png")

p_dune_restore <- base_map("Sand Dunes — Restoration Potential") +
  geom_sf(data = dune_restore, aes(fill = factor(tier)), colour = NA, alpha = 0.5) +
  scale_fill_manual(values = c("1" = "#fdae61", "2" = "#fee090"),
                    name = "Tier",
                    labels = c("1" = "Tier 1", "2" = "Tier 2"))
save_plot(p_dune_restore, "23_sand_dunes_restoration.png")

# 5e. Cluster sub-plots
for (i in seq_len(nrow(dune_boxes))) {
  bb <- st_bbox(dune_boxes[i, ])
  p_cl <- base_map(sprintf("Sand Dunes — %s", dune_boxes$name[i]), zoom = 14) +
    geom_sf(data = dune_restore, aes(fill = factor(tier)), colour = NA, alpha = 0.4) +
    geom_sf(data = dune_known, fill = "#d9b38c", colour = "#8c6d46", alpha = 0.7) +
    {if (nrow(dune_indicators) > 0)
      geom_sf(data = dune_indicators, colour = "#7a0177", size = 2.5, alpha = 0.7,
              shape = 15)} +
    scale_fill_manual(values = c("1" = "#fdae61", "2" = "#fee090"),
                      name = "Tier", labels = c("1" = "Tier 1", "2" = "Tier 2")) +
    coord_sf(xlim = c(bb["xmin"], bb["xmax"]),
             ylim = c(bb["ymin"], bb["ymax"]), crs = 27700)
  save_plot(p_cl, sprintf("24_sand_dunes_cluster_%02d.png", i), width = 8, height = 8)
}


# ============================================================
# 5b. COASTAL WOODLANDS
# ============================================================
cat("\n===== 5b. COASTAL WOODLANDS =====\n")

hlcm_path <- "data/raw/habitats/General/HLCM-2022-L2_SCOTLAND_TIFF_27700/HLCM_2022_EUNIS_LEVEL2.tif"
habmos_path <- "data/raw/habitats/General/Habitat_Map_of_Scotland_2181603415421412055.gpkg"

if (!file.exists(hlcm_path)) {
  cat("  WARNING: HLCM 2022 raster not found — skipping woodland plots\n")
} else {

  # --- Load HLCM and build lookup from raster category table ---
  hlcm_raw <- rast(hlcm_path)
  cat_tbl <- cats(hlcm_raw)[[1]]
  lc_lookup <- data.frame(
    lc_value   = cat_tbl$VALUE,
    eunis_code = cat_tbl$EUNIS_CODE,
    eunis_desc = cat_tbl$EUNIS_DESC,
    stringsAsFactors = FALSE
  )

  aoi_buf <- st_buffer(scamp_aoi, dist = 10000)
  hlcm_crop <- crop(hlcm_raw, ext(vect(aoi_buf)))
  hlcm_int <- as.numeric(hlcm_crop)
  NAflag(hlcm_int) <- -32768
  cat(sprintf("  HLCM cropped: %d x %d cells\n", ncol(hlcm_int), nrow(hlcm_int)))

  # --- Replicate HabMoS woodland reclassification from 05_localthreats.R ---
  g1_val  <- lc_lookup$lc_value[lc_lookup$eunis_code == "G1"]
  g3f_val <- lc_lookup$lc_value[lc_lookup$eunis_code == "G3.F"]

  habmos_wood <- NULL
  if (file.exists(habmos_path)) {
    habmos_all <- st_read(habmos_path, layer = st_layers(habmos_path)$name[1],
                          wkt_filter = st_as_text(st_as_sfc(st_bbox(aoi_buf))),
                          quiet = TRUE) %>% st_make_valid()

    habmos_g3f2 <- habmos_all[which(habmos_all$HABITAT_CODE == "G3.F2"), ]
    if (nrow(habmos_g3f2) > 0) {
      ex <- terra::extract(hlcm_int, vect(habmos_g3f2), cells = TRUE)
      hits <- ex$cell[ex[[2]] == g1_val]
      if (length(hits) > 0) hlcm_int[hits] <- g3f_val
      cat(sprintf("  G3.F2 → G3.F reclassified: %d pixels\n", length(hits)))
    }

    habmos_g4f <- habmos_all[which(habmos_all$HABITAT_CODE == "G4.F"), ]
    if (nrow(habmos_g4f) > 0) {
      ex <- terra::extract(hlcm_int, vect(habmos_g4f), cells = TRUE)
      hits <- ex$cell[ex[[2]] == g1_val]
      if (length(hits) > 0) hlcm_int[hits] <- 33
      cat(sprintf("  G4.F reclassified: %d pixels\n", length(hits)))
    }

    habmos_wood <- habmos_all[grepl("^G", habmos_all$HABITAT_CODE), ]
    cat(sprintf("  HabMoS woodland polygons: %d\n", nrow(habmos_wood)))
  } else {
    cat("  WARNING: HabMoS not found — using raw HLCM without reclassification\n")
  }

  # --- Define woodland classes and colours ---
  wood_classes <- data.frame(
    eunis_code = c("G1", "G3.F", "G4", "G5", "E5", "G4.F"),
    label = c("Broadleaved woodland", "Conifer plantation",
              "Mixed woodland", "Felled / early-stage",
              "Woodland fringes", "Mixed forestry plantation"),
    colour = c("#1b7837", "#7570b3", "#66c2a5", "#a6761d", "#b3de69", "#fc8d62"),
    stringsAsFactors = FALSE
  )

  wood_classes$lc_value <- sapply(wood_classes$eunis_code, function(code) {
    if (code == "G4.F") return(33L)
    v <- lc_lookup$lc_value[lc_lookup$eunis_code == code]
    if (length(v) == 0) return(NA_integer_)
    v
  })
  wood_classes <- wood_classes[!is.na(wood_classes$lc_value), ]

  rcl <- cbind(from = wood_classes$lc_value, to = seq_len(nrow(wood_classes)))
  hlcm_wood <- terra::classify(hlcm_int, rcl, others = NA)
  levels(hlcm_wood) <- data.frame(value = seq_len(nrow(wood_classes)),
                                   label = wood_classes$label)
  wood_pal <- setNames(wood_classes$colour, wood_classes$label)

  cat(sprintf("  Woodland classes present: %s\n",
              paste(wood_classes$label, collapse = ", ")))

  # --- Woodland indicator records ---
  wood_indicators <- st_sf(geometry = st_sfc(crs = 27700))
  if (file.exists("outputs/species/indicator_records.gpkg")) {
    all_ind <- st_read("outputs/species/indicator_records.gpkg", quiet = TRUE) %>%
      st_transform(27700)
    wood_indicators <- all_ind[all_ind$habitat == "Woodland", ]
    wood_indicators <- wood_indicators[lengths(st_intersects(wood_indicators, aoi_buf)) > 0, ]
    cat(sprintf("  Woodland indicator records: %d\n", nrow(wood_indicators)))
  }

  # --- 5b-1. HLCM reclassified woodland map ---
  cat("  Plotting HLCM woodland...\n")
  p_wood_hlcm <- base_map("Coastal Woodlands",
                           "HLCM 2022 reclassified with HabMoS corrections") +
    geom_spatraster(data = hlcm_wood, alpha = 0.85) +
    scale_fill_manual(values = wood_pal, name = "Woodland type",
                      na.translate = FALSE) +
    {if (nrow(wood_indicators) > 0)
      geom_sf(data = wood_indicators, aes(shape = "Woodland indicators"),
              colour = "#d62728", size = 1.5, alpha = 0.7)} +
    {if (nrow(wood_indicators) > 0)
      scale_shape_manual(values = c("Woodland indicators" = 17), name = NULL)}
  save_plot(p_wood_hlcm, "24a_woodland_hlcm.png")

  # --- 5b-2. HabMoS woodland map (comparison) ---
  if (!is.null(habmos_wood) && nrow(habmos_wood) > 0) {
    cat("  Plotting HabMoS woodland...\n")

    habmos_wood$wood_type <- case_when(
      grepl("^G1", habmos_wood$HABITAT_CODE) ~ "Native broadleaved",
      habmos_wood$HABITAT_CODE == "G3.F"     ~ "Conifer plantation",
      habmos_wood$HABITAT_CODE == "G3.F2"    ~ "Exotic conifer plantation",
      grepl("^G3", habmos_wood$HABITAT_CODE) ~ "Other coniferous",
      habmos_wood$HABITAT_CODE == "G4.F"     ~ "Mixed forestry plantation",
      grepl("^G4", habmos_wood$HABITAT_CODE) ~ "Mixed woodland",
      habmos_wood$HABITAT_CODE == "G5"       ~ "Felled / early-stage",
      TRUE                                   ~ habmos_wood$HABITAT_CODE
    )

    habmos_pal <- c(
      "Native broadleaved"        = "#1b7837",
      "Conifer plantation"        = "#7570b3",
      "Exotic conifer plantation" = "#e7298a",
      "Other coniferous"          = "#a6cee3",
      "Mixed woodland"            = "#66c2a5",
      "Mixed forestry plantation" = "#fc8d62",
      "Felled / early-stage"      = "#a6761d"
    )
    habmos_pal <- habmos_pal[names(habmos_pal) %in% unique(habmos_wood$wood_type)]

    p_wood_habmos <- base_map("Coastal Woodlands",
                               "Habitat Map of Scotland (HabMoS)") +
      geom_sf(data = habmos_wood, aes(fill = wood_type), colour = NA, alpha = 0.7) +
      scale_fill_manual(values = habmos_pal, name = "Woodland type")
    save_plot(p_wood_habmos, "24b_woodland_habmos.png")
  }
}


# ============================================================
# 6. SPECIES
# ============================================================
cat("\n===== 6. SPECIES =====\n")

# --- Helper: downscale 5 km grid → 500 m raster, masked to coastline ---
# Nearest-neighbour preserves original values; mask clips to boundary.
downscale_grid <- function(grid_sf, field, mask_sf, exclude_sf = NULL, fact = 10) {
  v <- vect(grid_sf)
  template <- rast(ext(v), resolution = 5000, crs = crs(v))
  r <- rasterize(v, template, field = field)
  r_fine <- disagg(r, fact = fact, method = "near")
  r_out <- mask(r_fine, vect(mask_sf))
  if (!is.null(exclude_sf)) r_out <- mask(r_out, vect(exclude_sf), inverse = TRUE)
  r_out
}

# --- Load species grids ---
rich_marine <- st_read("outputs/species/richness_grid_marine.gpkg", quiet = TRUE) %>%
  st_transform(27700)
rich_terr <- st_read("outputs/species/richness_grid_terrestrial.gpkg", quiet = TRUE) %>%
  st_transform(27700)
cat(sprintf("  Marine grid: %d cells\n", nrow(rich_marine)))
cat(sprintf("  Terrestrial grid: %d cells\n", nrow(rich_terr)))

# --- 6a–f. Paired marine/terrestrial metric maps (downscaled + coastline-masked) ---
metric_specs <- list(
  list(col = "species_richness", title = "Species Richness",
       opt = "mako", dir = -1, lab = "Richness", num = "25"),
  list(col = "record_count", title = "Recording Effort",
       opt = "inferno", dir = -1, lab = "Records", num = "26"),
  list(col = "median_year", title = "Median Recording Year",
       opt = NULL, dir = 1, lab = "Median year", num = "27",
       diverging = TRUE, mid = 2010),
  list(col = "most_recent_year", title = "Most Recent Record",
       opt = "plasma", dir = 1, lab = "Year", num = "28"),
  list(col = "beta_diversity", title = "Beta Diversity",
       opt = "turbo", dir = -1, lab = "Beta\ndiversity", num = "29"),
  list(col = "community_similarity", title = "Community Similarity",
       opt = "cividis", dir = -1, lab = "Similarity", num = "30")
)

for (spec in metric_specs) {
  for (realm in c("marine", "terrestrial")) {
    grid_sf <- if (realm == "marine") rich_marine else rich_terr
    if (!spec$col %in% names(grid_sf)) {
      cat(sprintf("  SKIP: %s not in %s grid\n", spec$col, realm))
      next
    }

    mask_sf <- if (realm == "marine") smr else dg
    excl_sf <- if (realm == "marine") NULL else smr
    r <- downscale_grid(grid_sf, spec$col, mask_sf, exclude_sf = excl_sf)

    rlabel <- if (realm == "marine") "Marine" else "Terrestrial"
    p <- base_map(sprintf("%s %s", rlabel, spec$title),
                  "5 km grid, nearest-neighbour to 500 m, coastline-masked") +
      geom_spatraster(data = r)

    if (!is.null(spec$diverging) && spec$diverging) {
      p <- p + scale_fill_gradient2(
        low = "#2166ac", mid = "white", high = "#b2182b",
        midpoint = spec$mid, name = spec$lab, na.value = NA)
    } else {
      p <- p + scale_fill_viridis_c(option = spec$opt, direction = spec$dir,
                                     name = spec$lab, na.value = NA)
    }

    suffix <- if (realm == "marine") "a" else "b"
    slug <- gsub(" ", "_", tolower(spec$title))
    save_plot(p, sprintf("%s%s_%s_%s.png", spec$num, suffix, slug, realm))
  }
}

# --- Load INNS records ---
inns <- st_read("outputs/species/inns_records.gpkg", quiet = TRUE) %>%
  st_transform(27700)
inns <- inns[lengths(st_intersects(inns, scamp_aoi)) > 0, ]
cat(sprintf("  INNS records in AOI: %d\n", nrow(inns)))

# --- Build 5 km grid for INNS analysis ---
grid_5km <- st_make_grid(scamp_aoi, cellsize = 5000, square = TRUE) %>%
  st_sf(cell_id = seq_along(.), geometry = .)
grid_5km <- grid_5km[lengths(st_intersects(grid_5km, scamp_aoi)) > 0, ]

inns_gridded <- st_join(inns, grid_5km, join = st_within) %>%
  st_drop_geometry() %>%
  filter(!is.na(cell_id))
cat(sprintf("  INNS records gridded: %d\n", nrow(inns_gridded)))

# --- 6g. INNS detections per grid cell ---
if (nrow(inns_gridded) > 0) {
  det_counts <- inns_gridded %>%
    group_by(cell_id) %>%
    summarise(n_detections = n(), .groups = "drop")
  det_grid <- grid_5km %>%
    left_join(det_counts, by = "cell_id") %>%
    filter(!is.na(n_detections))
  r_det <- downscale_grid(det_grid, "n_detections", scamp_aoi)

  p_det <- base_map("INNS Detections per Grid Cell",
                     "Total INNS records per 5 km cell, downscaled to 500 m") +
    geom_spatraster(data = r_det) +
    scale_fill_viridis_c(option = "rocket", direction = -1,
                         name = "Detections", na.value = NA)
  save_plot(p_det, "31_inns_detections_grid.png")
}

# --- 6h. Year of first INNS detection per grid cell ---
if (nrow(inns_gridded) > 0 && "year" %in% names(inns_gridded)) {
  first_det <- inns_gridded %>%
    filter(!is.na(year)) %>%
    group_by(cell_id) %>%
    summarise(first_year = min(year), .groups = "drop")
  first_grid <- grid_5km %>%
    left_join(first_det, by = "cell_id") %>%
    filter(!is.na(first_year))
  r_first <- downscale_grid(first_grid, "first_year", scamp_aoi)

  p_first <- base_map("Year of First INNS Detection",
                       "Earliest INNS record per 5 km cell, downscaled to 500 m") +
    geom_spatraster(data = r_first) +
    scale_fill_viridis_c(option = "plasma", direction = 1,
                         name = "First\ndetection\nyear", na.value = NA)
  save_plot(p_first, "32_inns_first_detection.png")
}

# --- 6i. Cumulative INNS spread — small multiples ---
if (nrow(inns_gridded) > 0 && "year" %in% names(inns_gridded)) {
  time_bands <- list(
    "Pre-2000"          = c(-Inf, 1999),
    "2000–2009"    = c(2000, 2009),
    "2010–2019"    = c(2010, 2019),
    "2020–present" = c(2020, Inf)
  )

  first_by_cell <- inns_gridded %>%
    filter(!is.na(year)) %>%
    group_by(cell_id) %>%
    summarise(first_year = min(year), .groups = "drop")

  cum_cells  <- integer(length(time_bands))
  cum_spp    <- integer(length(time_bands))
  band_panels <- list()

  for (b in seq_along(time_bands)) {
    band_range <- time_bands[[b]]
    cumulative_ids <- first_by_cell %>%
      filter(first_year <= band_range[2]) %>% pull(cell_id)
    new_ids <- first_by_cell %>%
      filter(first_year >= band_range[1], first_year <= band_range[2]) %>%
      pull(cell_id)

    cum_cells[b] <- length(cumulative_ids)
    cum_spp[b] <- inns_gridded %>%
      filter(!is.na(year), year <= band_range[2]) %>%
      summarise(n = n_distinct(latin_name)) %>% pull(n)

    band_grid <- grid_5km %>%
      filter(cell_id %in% cumulative_ids) %>%
      mutate(status = if_else(cell_id %in% new_ids,
                              "Newly occupied", "Previously occupied"))

    p_b <- ggplot() +
      annotation_map_tile(type = "osm", zoom = 9, progress = "none", quiet = TRUE) +
      geom_sf(data = scamp_aoi, fill = NA, colour = "black", linewidth = 0.4,
              linetype = "dashed")
    prev_cells <- band_grid %>% filter(status == "Previously occupied")
    new_cells  <- band_grid %>% filter(status == "Newly occupied")
    if (nrow(prev_cells) > 0)
      p_b <- p_b + geom_sf(data = prev_cells, fill = "grey60", colour = NA)
    if (nrow(new_cells) > 0)
      p_b <- p_b + geom_sf(data = new_cells, fill = "#d7191c", colour = NA)
    p_b <- p_b +
      coord_sf(crs = 27700) +
      labs(title = names(time_bands)[b],
           subtitle = sprintf("%d cells, %d species (cumulative)",
                              length(cumulative_ids), cum_spp[b])) +
      theme_scamp(base_size = 8) +
      theme(legend.position = "none")
    band_panels[[b]] <- p_b
  }

  p_inset <- ggplot(data.frame(band = seq_along(time_bands),
                                label = names(time_bands),
                                cells = cum_cells),
                    aes(x = band, y = cells)) +
    geom_line(colour = "#d7191c", linewidth = 0.8) +
    geom_point(colour = "#d7191c", size = 2) +
    scale_x_continuous(breaks = seq_along(time_bands), labels = names(time_bands)) +
    labs(x = NULL, y = "Cells") +
    theme_minimal(base_size = 7) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "white", colour = "grey60",
                                         linewidth = 0.3))

  band_panels[[length(time_bands)]] <- band_panels[[length(time_bands)]] +
    inset_element(p_inset, left = 0.55, bottom = 0.02, right = 0.98, top = 0.38)

  p_spread <- wrap_plots(band_panels, ncol = 2) +
    plot_annotation(
      title = "Cumulative INNS Spread across SCAMP AOI",
      subtitle = "Grey = previously occupied, red = newly occupied in this time band",
      theme = theme(plot.title = element_text(size = 14, face = "bold"),
                    plot.subtitle = element_text(size = 10, colour = "grey30")))
  save_plot(p_spread, "33_inns_cumulative_spread.png", width = 14, height = 12)
}

# --- 6j. Habitat-builder per-species maps ---
cat("  Habitat-builder maps...\n")
ind_val_path <- "outputs/species/indicator_species_validated.csv"
ind_rec_path <- "outputs/species/indicator_records.gpkg"

if (file.exists(ind_val_path) && file.exists(ind_rec_path)) {
  ind_validated <- read.csv(ind_val_path, stringsAsFactors = FALSE)
  all_ind <- st_read(ind_rec_path, quiet = TRUE) %>% st_transform(27700)

  if ("functional_group" %in% names(ind_validated)) {
    builders <- ind_validated[ind_validated$functional_group == "habitat_builder", ]
    n_built <- 0L

    for (bi in seq_len(nrow(builders))) {
      sp_latin  <- builders$latin_name[bi]
      sp_common <- builders$common_name[bi]
      sp_data   <- all_ind[all_ind$latin_name == sp_latin, ]
      if (nrow(sp_data) == 0) next

      sp_label <- if (!is.na(sp_common) && sp_common != "") sp_common else sp_latin

      p_sp <- base_map(sp_label, sp_latin)
      if ("year" %in% names(sp_data)) {
        p_sp <- p_sp +
          geom_sf(data = sp_data, aes(colour = year), size = 2, alpha = 0.8) +
          scale_colour_viridis_c(option = "plasma", direction = 1, name = "Year")
      } else {
        p_sp <- p_sp +
          geom_sf(data = sp_data, colour = "#d62728", size = 2, alpha = 0.8)
      }

      fname_clean <- gsub("[^a-zA-Z0-9]", "_", tolower(sp_label))
      save_plot(p_sp, sprintf("34_%s.png", fname_clean))
      n_built <- n_built + 1L
    }
    cat(sprintf("  %d habitat-builder species mapped\n", n_built))
  } else {
    cat("  WARNING: functional_group column not found in validated CSV\n")
  }
} else {
  cat("  WARNING: indicator data not found — skipping habitat-builder maps\n")
}

# --- 6k. PMF records ---
pmf_pts_path  <- "data/processed/species/pmf_points_scamp.gpkg"
pmf_poly_path <- "data/processed/species/pmf_polygons_scamp.gpkg"
pmf_pts   <- NULL
pmf_polys <- NULL

p_pmf <- base_map("Priority Marine Features", "GeMS records")
if (file.exists(pmf_poly_path)) {
  pmf_polys <- st_read(pmf_poly_path, quiet = TRUE) %>% st_transform(27700)
  if (nrow(pmf_polys) > 0)
    p_pmf <- p_pmf +
      geom_sf(data = pmf_polys, fill = "#df65b0", colour = "#ce1256", alpha = 0.3)
} else {
  cat("  WARNING: pmf_polygons_scamp.gpkg not found\n")
}
if (file.exists(pmf_pts_path)) {
  pmf_pts <- st_read(pmf_pts_path, quiet = TRUE) %>% st_transform(27700)
  if (nrow(pmf_pts) > 0)
    p_pmf <- p_pmf +
      geom_sf(data = pmf_pts, colour = "#ce1256", size = 0.8, alpha = 0.5)
} else {
  cat("  WARNING: pmf_points_scamp.gpkg not found\n")
}
save_plot(p_pmf, "35_pmf_records.png")

# --- 6l. INNS interactive leaflet (per-species, toggleable) ---
if (nrow(inns) > 0 && "latin_name" %in% names(inns)) {
  cat("  Building INNS interactive leaflet...\n")
  inns_wgs <- st_transform(inns, 4326)
  species_list <- sort(unique(inns_wgs$latin_name))
  sp_colours <- setNames(
    hcl.colors(length(species_list), palette = "Dynamic"),
    species_list)

  m_inns <- leaflet() %>%
    addProviderTiles(providers$CartoDB.Positron, group = "CartoDB") %>%
    addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
    addProviderTiles(providers$OpenStreetMap, group = "OSM") %>%
    addPolygons(data = st_transform(scamp_aoi, 4326), fill = FALSE,
                color = "black", weight = 2, dashArray = "5,5", group = "AOI")

  sp_group_labels <- character(length(species_list))
  for (si in seq_along(species_list)) {
    sp <- species_list[si]
    sp_data <- inns_wgs[inns_wgs$latin_name == sp, ]

    sp_label <- sp
    if ("common_name" %in% names(sp_data)) {
      cn <- unique(na.omit(sp_data$common_name))
      if (length(cn) > 0 && cn[1] != "")
        sp_label <- sprintf("%s (%s)", cn[1], sp)
    }
    sp_group_labels[si] <- sp_label

    popup_text <- sprintf("<b>%s</b>", sp_label)
    if ("year" %in% names(sp_data))
      popup_text <- paste0(popup_text, "<br>Year: ", sp_data$year)
    if ("realm" %in% names(sp_data))
      popup_text <- paste0(popup_text, "<br>Realm: ", sp_data$realm)

    m_inns <- m_inns %>%
      addCircleMarkers(data = sp_data, radius = 3,
                       color = sp_colours[sp], fillOpacity = 0.7, stroke = FALSE,
                       group = sp_label, popup = popup_text)
  }

  all_groups <- c("AOI", sp_group_labels)
  m_inns <- m_inns %>%
    addLayersControl(
      baseGroups = c("CartoDB", "Satellite", "OSM"),
      overlayGroups = all_groups,
      options = layersControlOptions(collapsed = TRUE))

  for (grp in sp_group_labels) m_inns <- m_inns %>% hideGroup(grp)

  inns_leaflet_path <- file.path(out_dir, "inns_interactive.html")
  saveWidget(m_inns, inns_leaflet_path, selfcontained = TRUE)
  cat(sprintf("  Saved: %s\n", inns_leaflet_path))
}


# ============================================================
# 7. THREATS
# ============================================================
cat("\n===== 7. THREATS =====\n")

# 7a. Fishing pressure
fish_gfw_path <- "outputs/threats/layers_raw/fishing_gfw.tif"
fish_u12_path <- "outputs/threats/layers_raw/fishing_under12.tif"

fish_panels <- list()

if (file.exists(fish_gfw_path)) {
  fish_gfw <- rast(fish_gfw_path)
  p_gfw <- ggplot() +
    annotation_map_tile(type = "osm", zoom = 9, progress = "none", quiet = TRUE) +
    geom_spatraster(data = fish_gfw) +
    geom_sf(data = dg, fill = NA, colour = "grey40", linewidth = 0.3) +
    geom_sf(data = smr, fill = NA, colour = "steelblue", linewidth = 0.4,
            linetype = "dotted") +
    geom_sf(data = scamp_aoi, fill = NA, colour = "black", linewidth = 0.5,
            linetype = "dashed") +
    scale_fill_viridis_c(option = "inferno", na.value = NA,
                         name = "Fishing hours") +
    annotation_scale(location = "bl", style = "bar", bar_cols = c("black", "white"),
                     text_cex = 0.6) +
    coord_sf(crs = 27700) +
    labs(title = "GFW AIS fishing (>15m, 2020-24)") +
    theme_scamp()
  fish_panels <- c(fish_panels, list(p_gfw))
}

if (file.exists(fish_u12_path)) {
  fish_u12 <- rast(fish_u12_path)
  u12_max <- max(values(fish_u12), na.rm = TRUE)
  u12_breaks <- seq(1, ceiling(u12_max))

  p_u12 <- ggplot() +
    annotation_map_tile(type = "osm", zoom = 9, progress = "none", quiet = TRUE) +
    geom_spatraster(data = fish_u12) +
    geom_sf(data = dg, fill = NA, colour = "grey40", linewidth = 0.3) +
    geom_sf(data = smr, fill = NA, colour = "steelblue", linewidth = 0.4,
            linetype = "dotted") +
    geom_sf(data = scamp_aoi, fill = NA, colour = "black", linewidth = 0.5,
            linetype = "dashed") +
    scale_fill_stepsn(colours = c("#fee5d9", "#fcae91", "#fb6a4a", "#de2d26",
                                   "#a50f15", "#67000d"),
                      breaks = u12_breaks,
                      na.value = NA, name = "Gear types present") +
    annotation_scale(location = "bl", style = "bar", bar_cols = c("black", "white"),
                     text_cex = 0.6) +
    coord_sf(crs = 27700) +
    labs(title = "Under-12m fishing effort",
         subtitle = "Count of gear categories per cell") +
    theme_scamp()
  fish_panels <- c(fish_panels, list(p_u12))
}

if (length(fish_panels) > 0) {
  p_fishing <- wrap_plots(fish_panels, ncol = length(fish_panels)) +
    plot_annotation(title = "Fishing Pressure",
                    theme = theme(plot.title = element_text(face = "bold", size = 14)))
  save_plot(p_fishing, "36_fishing_pressure.png",
            width = 7 * length(fish_panels), height = 8)
}

# 7b. Diffuse pollution — combined (river risk + marine exposure)
rivers_path   <- "outputs/catchments/rivers_with_accumulated_risk.gpkg"
exposure_path <- "outputs/catchments/marine_pollution_exposure.tif"
pp_path       <- "outputs/catchments/coastal_pour_points.gpkg"

poll_panels <- list()

if (file.exists(rivers_path)) {
  rivers_risk <- st_read(rivers_path, quiet = TRUE) %>% st_transform(27700)
  risk_col <- if ("accum_risk_raw" %in% names(rivers_risk)) "accum_risk_raw" else
    if ("accum_risk" %in% names(rivers_risk)) "accum_risk" else
    names(rivers_risk)[grep("risk", names(rivers_risk), ignore.case = TRUE)[1]]

  p_rivers <- base_map("River Accumulated Risk (SCIMAP)") +
    geom_sf(data = rivers_risk, aes(colour = .data[[risk_col]]), linewidth = 0.6) +
    scale_colour_viridis_c(option = "plasma", name = "Accumulated\nrisk")

  if (file.exists(pp_path)) {
    pp <- st_read(pp_path, quiet = TRUE) %>% st_transform(27700)
    p_rivers <- p_rivers +
      geom_sf(data = pp, colour = "red", size = 2.5, shape = 17)
  }
  poll_panels <- c(poll_panels, list(p_rivers))
  save_plot(p_rivers, "37_river_accumulated_risk.png")
}

if (file.exists(exposure_path)) {
  marine_exp <- rast(exposure_path)
  p_exposure <- ggplot() +
    annotation_map_tile(type = "osm", zoom = 9, progress = "none", quiet = TRUE) +
    geom_spatraster(data = marine_exp) +
    geom_sf(data = dg, fill = NA, colour = "grey40", linewidth = 0.3) +
    geom_sf(data = smr, fill = NA, colour = "steelblue", linewidth = 0.4,
            linetype = "dotted") +
    geom_sf(data = scamp_aoi, fill = NA, colour = "black", linewidth = 0.5,
            linetype = "dashed") +
    scale_fill_viridis_c(option = "magma", na.value = NA,
                         name = "Exposure") +
    annotation_scale(location = "bl", style = "bar", bar_cols = c("black", "white"),
                     text_cex = 0.6) +
    coord_sf(crs = 27700) +
    labs(title = "Marine Pollution Exposure (BFS decay)") +
    theme_scamp()
  poll_panels <- c(poll_panels, list(p_exposure))
  save_plot(p_exposure, "38_marine_pollution_exposure.png")
}

if (length(poll_panels) > 0) {
  p_pollution <- wrap_plots(poll_panels, ncol = length(poll_panels)) +
    plot_annotation(title = "Diffuse Pollution",
                    theme = theme(plot.title = element_text(face = "bold", size = 14)))
  save_plot(p_pollution, "39_diffuse_pollution_combined.png",
            width = 7 * length(poll_panels), height = 8)
}

# 7c. Cumulative pressure + vector threats
cum_path <- "outputs/threats/cumulative_pressure_with_pollution_1km.tif"
vec_path <- "outputs/threats/vector_threats_unified.gpkg"

threat_shape_map <- c(
  waste_discharge  = 16,
  dredge_spoil     = 17,
  dredging         = 15,
  ww_wrecks        = 8,
  dumped_munitions = 13,
  ports            = 18
)

vec_threats <- NULL
vec_pts     <- NULL
vec_polys   <- NULL

if (file.exists(vec_path)) {
  vec_threats <- st_read(vec_path, quiet = TRUE) %>% st_transform(27700)
  vec_pts <- vec_threats[st_geometry_type(vec_threats) %in%
                           c("POINT", "MULTIPOINT"), ]
  vec_polys <- vec_threats[st_geometry_type(vec_threats) %in%
                             c("POLYGON", "MULTIPOLYGON", "LINESTRING",
                               "MULTILINESTRING"), ]
}

if (file.exists(cum_path)) {
  cum_press <- rast(cum_path)

  # Combined: cumulative + vector overlaid
  p_cum <- ggplot() +
    annotation_map_tile(type = "osm", zoom = 9, progress = "none", quiet = TRUE) +
    geom_spatraster(data = cum_press) +
    geom_sf(data = dg, fill = NA, colour = "grey40", linewidth = 0.3) +
    geom_sf(data = smr, fill = NA, colour = "steelblue", linewidth = 0.4,
            linetype = "dotted") +
    geom_sf(data = scamp_aoi, fill = NA, colour = "black", linewidth = 0.5,
            linetype = "dashed") +
    scale_fill_viridis_c(option = "turbo", na.value = NA,
                         name = "Cumulative\npressure") +
    annotation_scale(location = "bl", style = "bar", bar_cols = c("black", "white"),
                     text_cex = 0.6) +
    annotation_north_arrow(location = "tr", which_north = "true",
                           height = unit(1, "cm"), width = unit(0.8, "cm"),
                           style = north_arrow_orienteering()) +
    coord_sf(crs = 27700) +
    labs(title = "Cumulative Pressure & Vector Threats") +
    theme_scamp()

  if (!is.null(vec_polys) && nrow(vec_polys) > 0) {
    p_cum <- p_cum +
      geom_sf(data = vec_polys, fill = NA, colour = "white", linewidth = 0.4)
  }
  if (!is.null(vec_pts) && nrow(vec_pts) > 0) {
    p_cum <- p_cum +
      geom_sf(data = vec_pts, aes(shape = threat_type), colour = "white",
              size = 2, alpha = 0.9) +
      scale_shape_manual(values = threat_shape_map, name = "Vector threat") +
      guides(shape = guide_legend(override.aes = list(size = 3)))
  }
  save_plot(p_cum, "40_cumulative_pressure_combined.png")

  # Cumulative only
  p_cum_only <- ggplot() +
    annotation_map_tile(type = "osm", zoom = 9, progress = "none", quiet = TRUE) +
    geom_spatraster(data = cum_press) +
    geom_sf(data = dg, fill = NA, colour = "grey40", linewidth = 0.3) +
    geom_sf(data = smr, fill = NA, colour = "steelblue", linewidth = 0.4,
            linetype = "dotted") +
    geom_sf(data = scamp_aoi, fill = NA, colour = "black", linewidth = 0.5,
            linetype = "dashed") +
    scale_fill_viridis_c(option = "turbo", na.value = NA,
                         name = "Cumulative\npressure") +
    annotation_scale(location = "bl", style = "bar", bar_cols = c("black", "white"),
                     text_cex = 0.6) +
    coord_sf(crs = 27700) +
    labs(title = "Cumulative Pressure Surface") +
    theme_scamp()
  save_plot(p_cum_only, "41_cumulative_pressure.png")
}

# Vector threats standalone
if (!is.null(vec_threats) && nrow(vec_threats) > 0) {
  p_vec <- base_map("Vector Threats")
  if (nrow(vec_polys) > 0) {
    p_vec <- p_vec +
      geom_sf(data = vec_polys, aes(fill = threat_type), colour = NA, alpha = 0.4) +
      scale_fill_brewer(palette = "Set2", name = "Threat (area)")
  }
  if (nrow(vec_pts) > 0) {
    p_vec <- p_vec +
      geom_sf(data = vec_pts, aes(shape = threat_type), size = 2.5, alpha = 0.8) +
      scale_shape_manual(values = threat_shape_map, name = "Threat (point)") +
      guides(shape = guide_legend(override.aes = list(size = 3)))
  }
  save_plot(p_vec, "42_vector_threats.png")
}


# ============================================================
# 8. INTERACTIVE LEAFLET
# ============================================================
cat("\n===== 8. INTERACTIVE LEAFLET =====\n")

to_wgs <- function(x) st_transform(x, 4326)

m <- leaflet() %>%
  addProviderTiles(providers$CartoDB.Positron, group = "CartoDB") %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
  addProviderTiles(providers$OpenStreetMap, group = "OSM")

# --- Overview ---
m <- m %>%
  addPolygons(data = to_wgs(scamp_aoi), fill = FALSE, color = "black", weight = 2,
              dashArray = "5,5", group = "AOI") %>%
  addPolygons(data = to_wgs(dg), fill = FALSE, color = "grey50", weight = 1,
              group = "D&G boundary") %>%
  addPolygons(data = to_wgs(smr), fill = FALSE, color = "steelblue", weight = 1,
              dashArray = "3,3", group = "SMR Solway") %>%
  addPolygons(data = to_wgs(harbour), fill = FALSE, color = "red", weight = 1,
              dashArray = "4,4", group = "Harbour limits")

# Designated sites
desig_colours_lf <- c(
  SAC = "#1b9e77", SPA = "#d95f02", SSSI = "#7570b3", Ramsar = "#e7298a",
  MPA = "#66a61e", NNR = "#e6ab02", LNR = "#a6761d", IBA = "#666666"
)
for (nm in names(desig_layers)) {
  d_wgs <- to_wgs(desig_layers[[nm]])
  m <- m %>%
    addPolygons(data = d_wgs, fillColor = desig_colours_lf[nm], fillOpacity = 0.25,
                color = desig_colours_lf[nm], weight = 1,
                group = paste0("Designated: ", nm),
                popup = ~paste0("<b>", nm, "</b>"))
}

# --- Seagrass ---
m <- m %>%
  addPolygons(data = to_wgs(sg_known), fillColor = "#1b7837", fillOpacity = 0.5,
              color = "#1b7837", weight = 1, group = "Seagrass: Known") %>%
  addPolygons(data = to_wgs(sg_potential), fillColor = "#a6d854", fillOpacity = 0.3,
              color = "#a6d854", weight = 1, group = "Seagrass: Potential")
if (nrow(sg_indicators) > 0) {
  m <- m %>%
    addCircleMarkers(data = to_wgs(sg_indicators), radius = 3,
                     color = "#4d004b", fillOpacity = 0.6, stroke = FALSE,
                     group = "Seagrass: Indicators",
                     popup = ~paste0("<b>", latin_name, "</b>"))
}

# --- Saltmarsh ---
m <- m %>%
  addPolygons(data = to_wgs(sm_known), fillColor = "#78c679", fillOpacity = 0.5,
              color = "#238443", weight = 1, group = "Saltmarsh: Known")
if (nrow(sm_ind) > 0) {
  m <- m %>%
    addCircleMarkers(data = to_wgs(sm_ind), radius = 3,
                     color = "#004529", fillOpacity = 0.6, stroke = FALSE,
                     group = "Saltmarsh: Indicators",
                     popup = ~paste0("<b>", latin_name, "</b>"))
}

# --- Oysters ---
m <- m %>%
  addPolygons(data = to_wgs(oy_known), fillColor = "#756bb1", fillOpacity = 0.5,
              color = "#54278f", weight = 1, group = "Oyster: Known") %>%
  addPolygons(data = to_wgs(oy_potential), fillColor = "#bcbddc", fillOpacity = 0.3,
              color = "#bcbddc", weight = 1, group = "Oyster: Potential")
if (nrow(oy_hist) > 0) {
  m <- m %>%
    addCircleMarkers(data = to_wgs(oy_hist), radius = 4,
                     color = "#de2d26", fillColor = "#fc9272", fillOpacity = 0.6,
                     stroke = TRUE, weight = 1, group = "Oyster: Historical")
}
if (nrow(oy_indicators) > 0) {
  m <- m %>%
    addCircleMarkers(data = to_wgs(oy_indicators), radius = 3,
                     color = "#e6550d", fillOpacity = 0.6, stroke = FALSE,
                     group = "Oyster: Indicators",
                     popup = ~paste0("<b>", latin_name, "</b>"))
}

# --- Sand Dunes ---
m <- m %>%
  addPolygons(data = to_wgs(dune_known), fillColor = "#d9b38c", fillOpacity = 0.5,
              color = "#8c6d46", weight = 1, group = "Sand Dunes: Known") %>%
  addPolygons(data = to_wgs(dune_restore), fillColor = "#fdae61", fillOpacity = 0.3,
              color = "#fdae61", weight = 1, group = "Sand Dunes: Restoration",
              popup = ~paste0("Tier ", tier))
if (nrow(dune_indicators) > 0) {
  m <- m %>%
    addCircleMarkers(data = to_wgs(dune_indicators), radius = 3,
                     color = "#7a0177", fillOpacity = 0.6, stroke = FALSE,
                     group = "Sand Dunes: Indicators",
                     popup = ~paste0("<b>", latin_name, "</b>"))
}

# --- Species ---
pal_marine <- colorNumeric("mako", rich_marine$species_richness, reverse = TRUE)
m <- m %>%
  addPolygons(data = to_wgs(rich_marine),
              fillColor = ~pal_marine(species_richness),
              fillOpacity = 0.5, color = NA, weight = 0,
              group = "Species: Marine richness",
              popup = ~paste0("Richness: ", species_richness))

pal_terr <- colorNumeric("viridis", rich_terr$species_richness, reverse = TRUE)
m <- m %>%
  addPolygons(data = to_wgs(rich_terr),
              fillColor = ~pal_terr(species_richness),
              fillOpacity = 0.5, color = NA, weight = 0,
              group = "Species: Terrestrial richness",
              popup = ~paste0("Richness: ", species_richness))

if (nrow(inns) > 0) {
  realm_cols <- c(Marine = "#2166ac", Terrestrial = "#b2182b", Freshwater = "#1b7837")
  pal_inns <- colorFactor(realm_cols, inns$realm)
  m <- m %>%
    addCircleMarkers(data = to_wgs(inns), radius = 2,
                     color = ~pal_inns(realm),
                     fillOpacity = 0.5, stroke = FALSE,
                     group = "Species: INNS",
                     popup = ~paste0("<b>", latin_name, "</b><br>", realm))
}

if (!is.null(pmf_pts) && nrow(pmf_pts) > 0) {
  m <- m %>%
    addCircleMarkers(data = to_wgs(pmf_pts), radius = 2,
                     color = "#ce1256", fillOpacity = 0.5, stroke = FALSE,
                     group = "Species: PMF",
                     popup = ~paste0("<b>", common_name, "</b>"))
}
if (!is.null(pmf_polys) && nrow(pmf_polys) > 0) {
  m <- m %>%
    addPolygons(data = to_wgs(pmf_polys), fillColor = "#df65b0", fillOpacity = 0.25,
                color = "#ce1256", weight = 1, group = "Species: PMF")
}

# --- INNS detection grid ---
if (exists("det_grid") && nrow(det_grid) > 0) {
  pal_det <- colorNumeric("YlOrRd", det_grid$n_detections)
  m <- m %>%
    addPolygons(data = to_wgs(det_grid),
                fillColor = ~pal_det(n_detections),
                fillOpacity = 0.5, color = NA, weight = 0,
                group = "Species: INNS detections",
                popup = ~paste0("Detections: ", n_detections))
}

# --- INNS first detection year ---
if (exists("first_grid") && nrow(first_grid) > 0) {
  pal_first <- colorNumeric("plasma", first_grid$first_year)
  m <- m %>%
    addPolygons(data = to_wgs(first_grid),
                fillColor = ~pal_first(first_year),
                fillOpacity = 0.5, color = NA, weight = 0,
                group = "Species: INNS first detection",
                popup = ~paste0("First detection: ", first_year))
}

# --- Habitat-builder indicator records ---
if (exists("all_ind") && exists("builders") && nrow(all_ind) > 0) {
  builder_recs <- all_ind[all_ind$latin_name %in% builders$latin_name, ]
  if (nrow(builder_recs) > 0) {
    builder_recs$popup_text <- paste0("<b>", builder_recs$latin_name, "</b>")
    if ("common_name" %in% names(builder_recs))
      builder_recs$popup_text <- paste0(builder_recs$popup_text,
                                         "<br>", builder_recs$common_name)
    if ("year" %in% names(builder_recs))
      builder_recs$popup_text <- paste0(builder_recs$popup_text,
                                         "<br>Year: ", builder_recs$year)
    m <- m %>%
      addCircleMarkers(data = to_wgs(builder_recs), radius = 3,
                       color = "#d62728", fillOpacity = 0.6, stroke = FALSE,
                       group = "Species: Habitat builders",
                       popup = ~popup_text)
  }
}

# --- Woodland (HabMoS) ---
if (exists("habmos_wood") && !is.null(habmos_wood) && nrow(habmos_wood) > 0) {
  habmos_wood_wgs <- to_wgs(habmos_wood)
  habmos_wood_wgs$popup_text <- habmos_wood_wgs$HABITAT_CODE
  if ("wood_type" %in% names(habmos_wood_wgs))
    habmos_wood_wgs$popup_text <- habmos_wood_wgs$wood_type
  habmos_wood_pal <- c(
    "Native broadleaved" = "#1b7837", "Conifer plantation" = "#7570b3",
    "Exotic conifer plantation" = "#e7298a", "Other coniferous" = "#a6cee3",
    "Mixed woodland" = "#66c2a5", "Mixed forestry plantation" = "#fc8d62",
    "Felled / early-stage" = "#a6761d")
  wood_col_var <- if ("wood_type" %in% names(habmos_wood_wgs)) "wood_type" else "HABITAT_CODE"
  pal_wood <- colorFactor(
    palette = unname(habmos_wood_pal),
    domain = unique(habmos_wood_wgs[[wood_col_var]]))
  habmos_wood_wgs$fill_colour <- pal_wood(habmos_wood_wgs[[wood_col_var]])
  m <- m %>%
    addPolygons(data = habmos_wood_wgs,
                fillColor = ~fill_colour,
                fillOpacity = 0.4, color = NA, weight = 0,
                group = "Woodland: HabMoS",
                popup = ~popup_text)
}

# --- Threats ---

# River risk (pre-compute popup to avoid .data issues in leaflet formulas)
if (exists("rivers_risk") && nrow(rivers_risk) > 0) {
  rivers_wgs <- to_wgs(rivers_risk)
  rivers_wgs$risk_value <- rivers_wgs[[risk_col]]
  rivers_wgs$popup_text <- paste0("Risk: ", round(rivers_wgs$risk_value, 2))
  pal_risk <- colorNumeric("plasma", rivers_wgs$risk_value)
  m <- m %>%
    addPolylines(data = rivers_wgs, color = ~pal_risk(risk_value),
                 weight = 2, opacity = 0.7, group = "Threats: River risk",
                 popup = ~popup_text)
}

# Pour points (pre-compute popup)
if (!exists("pp_path")) pp_path <- "outputs/catchments/coastal_pour_points.gpkg"
if (file.exists(pp_path)) {
  pp_wgs <- st_read(pp_path, quiet = TRUE) %>% st_transform(4326)
  pp_wgs$popup_text <- "Pour point"
  if ("watercourse_name" %in% names(pp_wgs))
    pp_wgs$popup_text <- pp_wgs$watercourse_name
  pp_risk_col <- if ("accum_risk_raw" %in% names(pp_wgs)) "accum_risk_raw" else
    if ("accum_risk" %in% names(pp_wgs)) "accum_risk" else NULL
  if (!is.null(pp_risk_col))
    pp_wgs$popup_text <- paste0(pp_wgs$popup_text,
                                 "<br>Risk: ", round(pp_wgs[[pp_risk_col]], 2))
  m <- m %>%
    addCircleMarkers(data = pp_wgs, radius = 4,
                     color = "red", fillColor = "red", fillOpacity = 0.7,
                     group = "Threats: Pour points", popup = ~popup_text)
}

# Vector threats
if (!is.null(vec_pts) && nrow(vec_pts) > 0) {
  m <- m %>%
    addCircleMarkers(data = to_wgs(vec_pts), radius = 3,
                     color = "orange", fillOpacity = 0.6, stroke = FALSE,
                     group = "Threats: Vector (points)",
                     popup = ~paste0("<b>", threat_type, "</b>"))
}
if (!is.null(vec_polys) && nrow(vec_polys) > 0) {
  m <- m %>%
    addPolygons(data = to_wgs(vec_polys), fillColor = "orange", fillOpacity = 0.2,
                color = "orange", weight = 1, group = "Threats: Vector (areas)",
                popup = ~paste0("<b>", threat_type, "</b>"))
}

# Raster threats (aggregate for performance, project to WGS84)
raster_threat_specs <- list(
  list(var = "fish_gfw",   name = "Threats: Fishing GFW",          pal = "inferno"),
  list(var = "fish_u12",   name = "Threats: Fishing under-12m",    pal = "YlOrRd"),
  list(var = "marine_exp", name = "Threats: Marine exposure",      pal = "magma"),
  list(var = "cum_press",  name = "Threats: Cumulative pressure",  pal = "turbo")
)
for (rt in raster_threat_specs) {
  if (exists(rt$var, inherits = FALSE)) {
    r_obj <- get(rt$var)
    r_agg <- terra::aggregate(r_obj, fact = 4, fun = "mean", na.rm = TRUE)
    r_wgs <- terra::project(r_agg, "EPSG:4326")
    vals <- terra::values(r_wgs)[!is.na(terra::values(r_wgs))]
    if (length(vals) > 0) {
      pal_r <- colorNumeric(rt$pal, range(vals), na.color = "transparent")
      m <- m %>%
        addRasterImage(raster::raster(r_wgs), colors = pal_r, opacity = 0.7,
                       group = rt$name)
      cat(sprintf("  Leaflet raster added: %s\n", rt$name))
    }
  }
}

# --- Legends ---

# Habitats — fixed-colour layers
m <- m %>%
  addLegend(position = "bottomright", group = "Seagrass: Known",
            colors = c("#1b7837", "#a6d854", "#4d004b"),
            labels = c("Known extent", "Modelled potential", "Indicator records"),
            title = "Seagrass", opacity = 0.8) %>%
  addLegend(position = "bottomright", group = "Saltmarsh: Known",
            colors = c("#78c679", "#004529"),
            labels = c("Known extent", "Indicator records"),
            title = "Saltmarsh", opacity = 0.8) %>%
  addLegend(position = "bottomright", group = "Oyster: Known",
            colors = c("#756bb1", "#bcbddc", "#de2d26", "#e6550d"),
            labels = c("Known extent", "Modelled potential", "Historical", "Indicator records"),
            title = "Oyster", opacity = 0.8) %>%
  addLegend(position = "bottomright", group = "Sand Dunes: Known",
            colors = c("#d9b38c", "#fdae61", "#7a0177"),
            labels = c("Known extent", "Restoration potential", "Indicator records"),
            title = "Sand Dunes", opacity = 0.8)

# Species — continuous palettes
m <- m %>%
  addLegend(position = "bottomright", group = "Species: Marine richness",
            pal = pal_marine, values = rich_marine$species_richness,
            title = "Marine richness", opacity = 0.8) %>%
  addLegend(position = "bottomright", group = "Species: Terrestrial richness",
            pal = pal_terr, values = rich_terr$species_richness,
            title = "Terrestrial richness", opacity = 0.8)

if (nrow(inns) > 0) {
  m <- m %>%
    addLegend(position = "bottomright", group = "Species: INNS",
              pal = pal_inns, values = inns$realm,
              title = "INNS by realm", opacity = 0.8)
}

m <- m %>%
  addLegend(position = "bottomright", group = "Species: PMF",
            colors = "#ce1256", labels = "PMF record",
            title = "PMF", opacity = 0.8) %>%
  addLegend(position = "bottomright", group = "Species: Habitat builders",
            colors = "#d62728", labels = "Builder record",
            title = "Habitat builders", opacity = 0.8)

if (exists("det_grid") && nrow(det_grid) > 0) {
  m <- m %>%
    addLegend(position = "bottomright", group = "Species: INNS detections",
              pal = pal_det, values = det_grid$n_detections,
              title = "INNS detections", opacity = 0.8)
}

if (exists("first_grid") && nrow(first_grid) > 0) {
  m <- m %>%
    addLegend(position = "bottomright", group = "Species: INNS first detection",
              pal = pal_first, values = first_grid$first_year,
              title = "First detection year", opacity = 0.8)
}

# Woodland
if (exists("habmos_wood") && !is.null(habmos_wood) && nrow(habmos_wood) > 0) {
  wood_types <- unique(habmos_wood_wgs[[wood_col_var]])
  m <- m %>%
    addLegend(position = "bottomright", group = "Woodland: HabMoS",
              pal = pal_wood, values = wood_types,
              title = "Woodland type", opacity = 0.8)
}

# Threats
if (exists("rivers_risk") && nrow(rivers_risk) > 0) {
  m <- m %>%
    addLegend(position = "bottomleft", group = "Threats: River risk",
              pal = pal_risk, values = rivers_wgs$risk_value,
              title = "River risk", opacity = 0.8)
}

m <- m %>%
  addLegend(position = "bottomleft", group = "Threats: Pour points",
            colors = "red", labels = "Pour point",
            title = "Pour points", opacity = 0.8) %>%
  addLegend(position = "bottomleft", group = "Threats: Vector (points)",
            colors = "orange", labels = "Vector threat",
            title = "Vector threats", opacity = 0.8)

# Raster threat legends (stored during the addRasterImage loop above)
# Re-create palettes for legends
for (rt in raster_threat_specs) {
  if (exists(rt$var, inherits = FALSE)) {
    r_obj <- get(rt$var)
    r_agg <- terra::aggregate(r_obj, fact = 4, fun = "mean", na.rm = TRUE)
    r_wgs <- terra::project(r_agg, "EPSG:4326")
    vals <- terra::values(r_wgs)[!is.na(terra::values(r_wgs))]
    if (length(vals) > 0) {
      pal_r <- colorNumeric(rt$pal, range(vals), na.color = "transparent")
      m <- m %>%
        addLegend(position = "bottomleft", group = rt$name,
                  pal = pal_r, values = range(vals),
                  title = gsub("Threats: ", "", rt$name), opacity = 0.8)
    }
  }
}

# --- Layer control ---
overlay_groups <- c(
  "AOI", "D&G boundary", "SMR Solway", "Harbour limits",
  paste0("Designated: ", names(desig_layers)),
  "Seagrass: Known", "Seagrass: Potential", "Seagrass: Indicators",
  "Saltmarsh: Known", "Saltmarsh: Indicators",
  "Oyster: Known", "Oyster: Historical", "Oyster: Potential", "Oyster: Indicators",
  "Sand Dunes: Known", "Sand Dunes: Restoration", "Sand Dunes: Indicators",
  "Woodland: HabMoS",
  "Species: Marine richness", "Species: Terrestrial richness",
  "Species: INNS", "Species: INNS detections", "Species: INNS first detection",
  "Species: Habitat builders", "Species: PMF",
  "Threats: River risk", "Threats: Pour points",
  "Threats: Vector (points)", "Threats: Vector (areas)",
  "Threats: Fishing GFW", "Threats: Fishing under-12m",
  "Threats: Marine exposure", "Threats: Cumulative pressure"
)

m <- m %>%
  addLayersControl(
    baseGroups = c("CartoDB", "Satellite", "OSM"),
    overlayGroups = overlay_groups,
    options = layersControlOptions(collapsed = FALSE)
  )

hide_by_default <- overlay_groups[!overlay_groups %in%
  c("AOI", "D&G boundary", "SMR Solway")]
for (grp in hide_by_default) {
  m <- m %>% hideGroup(grp)
}

leaflet_path <- file.path(out_dir, "scamp_interactive_map.html")
saveWidget(m, leaflet_path, selfcontained = FALSE)
cat(sprintf("  Saved: %s\n", leaflet_path))

cat("\n========== 07_plots.R COMPLETE ==========\n")
cat(sprintf("All figures saved to: %s/\n", out_dir))


# ============================================================
# PRIORITY MAPS — PLAN (to be implemented below)
# ============================================================
#
# Goal: produce aggregate spatial layers that help prioritise where to
# protect, where to restore, and where to monitor for biodiversity uplift.
#
# ---------------------
# LAYER 0: PRESSURE COUNT (management complexity)
# ---------------------
# For each 1km grid cell in the SMR, count how many of the ~15-20 normalised
# pressure layers (from 05_localthreats.R) the cell falls in the top 50%.
# This gives a "breadth of challenge" metric — a cell scoring 8/15 faces
# pressures across many fronts regardless of individual intensity. Distinct
# from cumulative pressure (which sums intensity). A high pressure count
# means more different things to manage.
#
# ---------------------
# LAYER 1: REFERENCE SITES
# ---------------------
# Named reference sites per habitat — these define "good condition":
#   - Sand dunes:  RSPB Mersehead & Caerlaverock
#   - Saltmarsh:   RSPB Mersehead & Crook of Baldoon
#   - Oyster:      Loch Ryan (Heriot-Watt may provide reference data)
#   - Seagrass:    SFP monitoring data if available; otherwise top 10%
#                  richness + bottom 10% cumulative pressure where known
#                  seagrass is present
#
# For each reference site, extract the community composition profile
# (species list, richness, community similarity baseline). These profiles
# become the benchmark that other cells are compared against in layers 2-4.
# Where no named reference exists, use top/bottom 10% filter as fallback.
# Load RSPB boundaries from data/processed/boundaries/rspb_reserves_scamp.gpkg.
#
# ---------------------
# LAYER 2: PROTECTION PRIORITY
# ---------------------
# Cells with high existing biodiversity value under threat, not already
# formally protected. Score from:
#   - Community similarity to reference sites (whole-seascape Sorensen,
#     already computed in 02_themeD_species.R) — high similarity = high
#     current ecological value
#   - Richness relative to reference site richness
#   - Pressure count (Layer 0) — how many different threats need managing
#   - Protection gap — outside existing SAC/SPA/SSSI polygons
#
# Interpretation of similarity + richness combinations:
#   - High similarity + high richness = genuinely healthy → PROTECT
#   - High similarity + low richness  = may be misleading → investigate
#   - Low similarity  + high richness = unique community  → PROTECT (different)
#   - Low similarity  + low richness  = degraded          → candidate for
#                                                           restoration instead
#
# ---------------------
# LAYER 3a: ACTIVE RESTORATION PRIORITY
# ---------------------
# Where to deploy cultch, plant seagrass, restore dune morphology — spending
# capital that you don't want wasted.
#   - Restoration opportunity: potential extent minus known extent
#     (seagrass: requires intertidal bathy fix; oyster: historical locations
#     that pass depth/substrate/salinity filters from 04_physical.R)
#   - LOW current pressure — maximises return on investment
#   - Proximity to source populations (coastalNet connectivity)
#   - Physical suitability (passes all environmental filters)
#
# ---------------------
# LAYER 3b: NATURAL RECOVERY SUPPORT PRIORITY
# ---------------------
# Where reducing pressures could accelerate recovery that's already plausible.
#   - Same restoration opportunity layer as 3a
#   - MEDIUM pressure band — too high and intervention fails; too low and
#     there are no pressures to reduce, so you can't actively accelerate
#     recovery through management
#   - Connectivity to healthy source populations (coastalNet dispersal)
#
# ---------------------
# LAYER 4: MONITORING PRIORITY
# ---------------------
# Where to detect recovery — prioritise where expected uplift is greatest,
# avoiding ceiling effects from monitoring already-healthy ecosystems.
#   - Hyperlocal beta diversity (queen-contiguity) is HIGH but richness
#     is LOW — the cell's neighbours are diverse but the cell itself is
#     depauperate. Strong recovery signal expected if pressures reduce.
#   - Low community similarity to reference sites — more room for
#     measurable convergence toward the reference profile
#   - Connectivity to healthy source populations (coastalNet) — cells
#     downstream of thriving habitat most likely to show colonisation
#   - Pressure reduction potential — links to Layer 3b areas
#
# ---------------------
# LAYER 5: COMBINED ACTION MAP
# ---------------------
# Classify each cell by dominant score:
#   protect / restore (active) / restore (support) / monitor / data gap
# Add as interactive leaflet layer with popup showing all component scores.
# Produce static ggplot version for report, facetable by habitat.
#
# ---------------------
# DEPENDENCIES
# ---------------------
# - coastalNet package (https://github.com/jorgeassis/coastalNet):
#   load and explore in 04_physical.R, save connectivity matrix for SMR
#   coastal cells. Feeds into restoration (which sites reachable from
#   surviving populations) and monitoring (which sites downstream of
#   sources, most likely to show colonisation).
# - Intertidal bathymetry fix in 04_physical.R: GEBCO (~460m) cannot
#   resolve narrow intertidal estuaries like the Rough Firth where known
#   seagrass exists. Need higher-res DTM or adjusted upper depth limit
#   with SMR masking to prevent flagging terrestrial areas.
# - RSPB reserve boundaries: data/processed/boundaries/rspb_reserves_scamp.gpkg
# - Community similarity matrices: outputs/species/community_similarity_matrix_*.rds
# - Normalised pressure rasters: outputs/threats/normalised_layers/
# ============================================================
