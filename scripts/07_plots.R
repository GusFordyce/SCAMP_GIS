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
  IBA    = "data/processed/boundaries/ibas_scamp.gpkg"
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

# Cluster detection: 2km buffer, union, decompose — returns sf of bounding boxes
# Replicates the approach from 03_habitats.R
build_cluster_boxes <- function(sf_data, pad_m = 2000, merge_dist = 1000) {
  if (nrow(sf_data) == 0) return(st_sf(name = character(), geometry = st_sfc(crs = 27700)))
  centroids <- suppressWarnings(st_centroid(sf_data))
  buf <- st_buffer(centroids, dist = 2000) %>% st_union() %>% st_cast("POLYGON")
  ids <- st_intersects(centroids, buf)
  sf_data$cluster <- sapply(ids, `[`, 1)

  boxes <- do.call(rbind, lapply(sort(unique(sf_data$cluster)), function(cl) {
    bb <- st_bbox(st_buffer(sf_data[which(sf_data$cluster == cl), ], dist = pad_m))
    st_sf(
      name = sprintf("Cluster %d", cl),
      geometry = st_as_sfc(bb, crs = 27700)
    )
  }))

  # Merge overlapping boxes
  merge_buf <- st_buffer(boxes, dist = merge_dist) %>% st_union() %>%
    st_cast("POLYGON")
  grp <- sapply(st_intersects(boxes, merge_buf), `[`, 1)
  boxes <- do.call(rbind, lapply(sort(unique(grp)), function(g) {
    bb <- st_bbox(boxes[grp == g, ])
    st_sf(name = sprintf("Cluster %d", g),
          geometry = st_as_sfc(bb, crs = 27700))
  }))

  # Order west to east by centroid
  boxes <- boxes[order(st_coordinates(suppressWarnings(st_centroid(boxes)))[, 1]), ]
  boxes$name <- sprintf("Cluster %d", seq_len(nrow(boxes)))
  boxes
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
  st_transform(27700)
sg_indicators <- ind_all[which(ind_all$habitat == "Seagrass"), ]
sg_indicators <- sg_indicators[lengths(st_intersects(sg_indicators, smr)) > 0, ]
cat(sprintf("  Seagrass indicators (within SMR): %d\n", nrow(sg_indicators)))

# Build cluster boxes from all seagrass data (same approach as 03_habitats.R)
sg_all_for_clusters <- bind_rows(sg_known, sg_potential) %>% st_make_valid()
sg_boxes <- build_cluster_boxes(sg_all_for_clusters)
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
  p_cl <- base_map(sprintf("Seagrass — %s", sg_boxes$name[i]), zoom = 12) +
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

# Build saltmarsh cluster boxes
sm_boxes <- build_cluster_boxes(sm_known)
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
  p_cl <- base_map(sprintf("Saltmarsh — %s", sm_boxes$name[i]), zoom = 12) +
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

oy_hist_path <- "data/processed/species/oyster_historical_extent.gpkg"
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

# Build oyster cluster boxes from known + historical
oy_for_clusters <- bind_rows(
  if (nrow(oy_known) > 0) oy_known else NULL,
  if (nrow(oy_hist) > 0) oy_hist else NULL
) %>% st_make_valid()
oy_boxes <- if (nrow(oy_for_clusters) > 0) {
  build_cluster_boxes(oy_for_clusters)
} else {
  st_sf(name = character(), geometry = st_sfc(crs = 27700))
}
cat(sprintf("  Oyster clusters: %d\n", nrow(oy_boxes)))

# 4a. Combined
p_oy_combined <- base_map("Native Oysters",
                           "Known, historical, modelled potential & indicator species") +
  geom_sf(data = oy_potential, aes(fill = "Modelled potential"), colour = NA,
          alpha = 0.4) +
  {if (nrow(oy_hist) > 0)
    geom_sf(data = oy_hist, aes(fill = "Historical"), colour = "#de2d26",
            alpha = 0.4)} +
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
    geom_sf(data = oy_hist, fill = "#fc9272", colour = "#de2d26", alpha = 0.5)
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
  p_cl <- base_map(sprintf("Native Oysters — %s", oy_boxes$name[i]), zoom = 12) +
    geom_sf(data = oy_potential, aes(fill = "Modelled potential"), colour = NA,
            alpha = 0.4) +
    {if (nrow(oy_hist) > 0)
      geom_sf(data = oy_hist, aes(fill = "Historical"), colour = "#de2d26",
              alpha = 0.4)} +
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

# Build dune cluster boxes
dune_boxes <- build_cluster_boxes(dune_known)
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
  p_cl <- base_map(sprintf("Sand Dunes — %s", dune_boxes$name[i]), zoom = 12) +
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
# 6. SPECIES
# ============================================================
cat("\n===== 6. SPECIES =====\n")

# 6a. Marine species richness
rich_marine <- st_read("outputs/species/richness_grid_marine.gpkg", quiet = TRUE) %>%
  st_transform(27700)

p_rich_marine <- base_map("Marine Species Richness", "5 km grid") +
  geom_sf(data = rich_marine, aes(fill = species_richness), colour = NA, alpha = 0.7) +
  scale_fill_viridis_c(option = "mako", direction = -1, name = "Richness",
                       na.value = NA)
save_plot(p_rich_marine, "25_marine_richness.png")

# 6b. Terrestrial species richness
rich_terr <- st_read("outputs/species/richness_grid_terrestrial.gpkg", quiet = TRUE) %>%
  st_transform(27700)

p_rich_terr <- base_map("Terrestrial Species Richness", "5 km grid") +
  geom_sf(data = rich_terr, aes(fill = species_richness), colour = NA, alpha = 0.7) +
  scale_fill_viridis_c(option = "viridis", direction = -1, name = "Richness",
                       na.value = NA)
save_plot(p_rich_terr, "26_terrestrial_richness.png")

# 6c. INNS distribution — cropped to SCAMP AOI
inns <- st_read("outputs/species/inns_records.gpkg", quiet = TRUE) %>%
  st_transform(27700)
inns <- inns[lengths(st_intersects(inns, scamp_aoi)) > 0, ]
cat(sprintf("  INNS records in AOI: %d\n", nrow(inns)))

p_inns <- base_map("Invasive Non-Native Species", "Records within SCAMP AOI") +
  geom_sf(data = inns, aes(colour = realm), size = 0.8, alpha = 0.5) +
  scale_colour_manual(values = c(Marine = "#2166ac", Terrestrial = "#b2182b",
                                  Freshwater = "#1b7837"),
                      name = "Realm")
save_plot(p_inns, "27_inns_distribution.png")

# 6d. PMF records
pmf_pts_path  <- "data/processed/species/pmf_points_scamp.gpkg"
pmf_poly_path <- "data/processed/species/pmf_polygons_scamp.gpkg"

pmf_pts   <- NULL
pmf_polys <- NULL

p_pmf <- base_map("Priority Marine Features", "GeMS records")

if (file.exists(pmf_poly_path)) {
  pmf_polys <- st_read(pmf_poly_path, quiet = TRUE) %>% st_transform(27700)
  cat(sprintf("  PMF polygons loaded: %d features\n", nrow(pmf_polys)))
  if (nrow(pmf_polys) > 0) {
    p_pmf <- p_pmf +
      geom_sf(data = pmf_polys, fill = "#df65b0", colour = "#ce1256", alpha = 0.3)
  }
} else {
  cat("  WARNING: pmf_polygons_scamp.gpkg not found — run 02_themeD_species.R first\n")
}

if (file.exists(pmf_pts_path)) {
  pmf_pts <- st_read(pmf_pts_path, quiet = TRUE) %>% st_transform(27700)
  cat(sprintf("  PMF points loaded: %d features\n", nrow(pmf_pts)))
  if (nrow(pmf_pts) > 0) {
    p_pmf <- p_pmf +
      geom_sf(data = pmf_pts, colour = "#ce1256", size = 0.8, alpha = 0.5)
  }
} else {
  cat("  WARNING: pmf_points_scamp.gpkg not found — run 02_themeD_species.R first\n")
}

save_plot(p_pmf, "28_pmf_records.png")


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
  save_plot(p_fishing, "29_fishing_pressure.png",
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
  save_plot(p_rivers, "31_river_accumulated_risk.png")
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
  save_plot(p_exposure, "32_marine_pollution_exposure.png")
}

if (length(poll_panels) > 0) {
  p_pollution <- wrap_plots(poll_panels, ncol = length(poll_panels)) +
    plot_annotation(title = "Diffuse Pollution",
                    theme = theme(plot.title = element_text(face = "bold", size = 14)))
  save_plot(p_pollution, "30_diffuse_pollution_combined.png",
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
  save_plot(p_cum, "33_cumulative_pressure_combined.png")

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
  save_plot(p_cum_only, "34_cumulative_pressure.png")
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
  save_plot(p_vec, "35_vector_threats.png")
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
    addPolygons(data = to_wgs(oy_hist), fillColor = "#fc9272", fillOpacity = 0.3,
                color = "#de2d26", weight = 1, group = "Oyster: Historical")
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

# --- Threats ---
if (exists("rivers_risk")) {
  pal_risk <- colorNumeric("plasma", rivers_risk[[risk_col]])
  m <- m %>%
    addPolylines(data = to_wgs(rivers_risk),
                 color = ~pal_risk(rivers_risk[[risk_col]]),
                 weight = 2, opacity = 0.7, group = "Threats: River risk",
                 popup = ~paste0("Risk: ", round(.data[[risk_col]], 2)))
}

if (file.exists(pp_path)) {
  pp_wgs <- st_read(pp_path, quiet = TRUE) %>% st_transform(4326)
  pp_risk_col <- if ("accum_risk_raw" %in% names(pp_wgs)) "accum_risk_raw" else
    if ("accum_risk" %in% names(pp_wgs)) "accum_risk" else NULL
  m <- m %>%
    addCircleMarkers(data = pp_wgs, radius = 4,
                     color = "red", fillColor = "red", fillOpacity = 0.7,
                     group = "Threats: Pour points",
                     popup = ~paste0(
                       if ("watercourse_name" %in% names(pp_wgs)) watercourse_name
                       else "Pour point",
                       if (!is.null(pp_risk_col))
                         paste0("<br>Risk: ", round(.data[[pp_risk_col]], 2))
                       else ""))
}

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

# --- Layer control ---
overlay_groups <- c(
  "AOI", "D&G boundary", "SMR Solway", "Harbour limits",
  paste0("Designated: ", names(desig_layers)),
  "Seagrass: Known", "Seagrass: Potential", "Seagrass: Indicators",
  "Saltmarsh: Known", "Saltmarsh: Indicators",
  "Oyster: Known", "Oyster: Historical", "Oyster: Potential", "Oyster: Indicators",
  "Sand Dunes: Known", "Sand Dunes: Restoration", "Sand Dunes: Indicators",
  "Species: Marine richness", "Species: Terrestrial richness", "Species: INNS",
  "Species: PMF",
  "Threats: River risk", "Threats: Pour points",
  "Threats: Vector (points)", "Threats: Vector (areas)"
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
saveWidget(m, leaflet_path, selfcontained = TRUE)
cat(sprintf("  Saved: %s\n", leaflet_path))

cat("\n========== 07_plots.R COMPLETE ==========\n")
cat(sprintf("All figures saved to: %s/\n", out_dir))
