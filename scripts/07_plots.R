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
library(vegan)
library(gstat)
library(tidyr)

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

# Fast base map: vector boundaries only, no OSM tile fetch. Use for
# diagnostic/priority-map plots where annotation_map_tile()'s network
# round-trip (fresh tiles per call) is not worth the land-context detail.
base_map_fast <- function(title = NULL, subtitle = NULL) {
  ggplot() +
    geom_sf(data = dg, fill = "grey92", colour = "grey40", linewidth = 0.3) +
    geom_sf(data = smr, fill = NA, colour = "steelblue", linewidth = 0.4,
            linetype = "dotted") +
    geom_sf(data = scamp_aoi, fill = NA, colour = "black", linewidth = 0.5,
            linetype = "dashed") +
    coord_sf(crs = 27700) +
    labs(title = title, subtitle = subtitle) +
    theme_minimal()
}

# Woodland HABITAT_CODE -> readable label. Shared by the static plot and
# the leaflet map so they can't drift apart. Bounded fallback ("Other
# woodland"), not raw text — a fixed palette (scale_fill_manual /
# colorFactor) needs a fixed, known set of category values, so falling
# back to free-text HABITAT_NAME would break it on any unmapped code.
# Comprehensive grepl coverage (^G5 not =="G5") so G5.1/G5.8 etc. don't
# leak through as raw codes the way they did before.
classify_wood_type <- function(habitat_code) {
  case_when(
    grepl("^G1", habitat_code) ~ "Native broadleaved",
    habitat_code == "G3.F"     ~ "Conifer plantation",
    habitat_code == "G3.F2"    ~ "Exotic conifer plantation",
    grepl("^G3", habitat_code) ~ "Other coniferous",
    habitat_code == "G4.F"     ~ "Mixed forestry plantation",
    grepl("^G4", habitat_code) ~ "Mixed woodland",
    grepl("^G5", habitat_code) ~ "Felled / early-stage",
    TRUE                       ~ "Other woodland"
  )
}

wood_type_pal <- c(
  "Native broadleaved"        = "#1b7837",
  "Conifer plantation"        = "#7570b3",
  "Exotic conifer plantation" = "#e7298a",
  "Other coniferous"          = "#a6cee3",
  "Mixed woodland"            = "#66c2a5",
  "Mixed forestry plantation" = "#fc8d62",
  "Felled / early-stage"      = "#a6761d",
  "Other woodland"            = "#999999"
)

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

    dir.create("data/processed/habitats", recursive = TRUE, showWarnings = FALSE)
    st_write(habmos_wood, "data/processed/habitats/woodland_habmos_solway.gpkg",
             delete_layer = TRUE, quiet = TRUE)
    cat("  Saved: data/processed/habitats/woodland_habmos_solway.gpkg\n")

    ## Save the full HabMoS-corrected HLCM raster + its lc_value lookup —
    ## this is a general land-cover product (all EUNIS classes: grassland,
    ## heathland etc., not just woodland), reusable beyond just extracting
    ## woodland from it
    writeRaster(hlcm_int, "data/processed/habitats/hlcm_habmos_corrected.tif",
                overwrite = TRUE, gdal = c("COMPRESS=LZW"))
    write.csv(lc_lookup, "data/processed/habitats/hlcm_lc_lookup.csv", row.names = FALSE)
    cat("  Saved: data/processed/habitats/hlcm_habmos_corrected.tif\n")
    cat("  Saved: data/processed/habitats/hlcm_lc_lookup.csv\n")
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

    habmos_wood$wood_type <- classify_wood_type(habmos_wood$HABITAT_CODE)

    habmos_pal <- wood_type_pal[names(wood_type_pal) %in% unique(habmos_wood$wood_type)]

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
  ## Computed unconditionally — don't rely on wood_type having already
  ## been set earlier in the session (that fragile check is exactly what
  ## let raw HABITAT_CODE values like "G5.1"/"G5.8" leak into the legend)
  habmos_wood_wgs$wood_type <- classify_wood_type(habmos_wood_wgs$HABITAT_CODE)
  habmos_wood_wgs$popup_text <- habmos_wood_wgs$wood_type

  pal_wood <- colorFactor(
    palette = unname(wood_type_pal),
    levels = names(wood_type_pal))
  habmos_wood_wgs$fill_colour <- pal_wood(habmos_wood_wgs$wood_type)
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
  wood_types <- unique(habmos_wood_wgs$wood_type)
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

cat("\n========== Static + interactive plots complete ==========\n")
cat(sprintf("All figures saved to: %s/\n", out_dir))


# ============================================================
# PRIORITY MAPS
# ============================================================

cat("\n\n========== PRIORITY MAPS ==========\n")

dir.create(file.path(out_dir, "priority"), showWarnings = FALSE, recursive = TRUE)

# ============================================================
# LAYER 0: PRESSURE COUNT (management complexity)
# ============================================================

cat("\n-- Layer 0: Pressure count --\n")

norm_stack_path <- "outputs/threats/pressure_layers_normalised_1km.tif"
stopifnot("Normalised pressure stack not found" = file.exists(norm_stack_path))

norm_stack <- rast(norm_stack_path)
cat(sprintf("  Pressure layers: %d (%s)\n",
            nlyr(norm_stack), paste(names(norm_stack), collapse = ", ")))

pressure_count <- app(norm_stack, fun = function(x) {
  sum(x > 0.5, na.rm = TRUE)
})
names(pressure_count) <- "pressure_count"
pressure_count <- mask(pressure_count, vect(smr))

cat(sprintf("  Pressure count range: %d to %d\n",
            global(pressure_count, "min", na.rm = TRUE)[[1]],
            global(pressure_count, "max", na.rm = TRUE)[[1]]))

pressure_count_df <- as.data.frame(pressure_count, xy = TRUE)

p_pressure_count <- base_map("Layer 0 — Pressure Count",
  sprintf("Number of pressure layers in top 50%% per cell (of %d layers)",
          nlyr(norm_stack))) +
  geom_tile(data = pressure_count_df,
            aes(x = x, y = y, fill = pressure_count)) +
  scale_fill_viridis_c(name = "Pressure\ncount", option = "inferno",
                       direction = -1, na.value = NA)

save_plot(p_pressure_count, "priority/L0_pressure_count.png")

writeRaster(pressure_count, "outputs/threats/pressure_count_1km.tif",
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))
cat("  Saved: outputs/threats/pressure_count_1km.tif\n")


# ============================================================
# LAYER 1: REFERENCE SITES
# ============================================================

cat("\n-- Layer 1: Reference sites --\n")

## Load RSPB reserves
rspb <- desig_layers[["RSPB"]]
if (is.null(rspb)) {
  rspb_path <- "data/processed/boundaries/rspb_reserves_scamp.gpkg"
  if (file.exists(rspb_path)) {
    rspb <- st_read(rspb_path, quiet = TRUE) %>% st_transform(27700)
  }
}
rspb <- st_make_valid(rspb)
rspb_name_col <- names(rspb)[grepl("name|Name|NAME|SITE", names(rspb))][1]
cat(sprintf("  RSPB reserves: %d features\n", nrow(rspb)))
if (!is.na(rspb_name_col)) {
  cat(sprintf("  Sites: %s\n",
              paste(unique(rspb[[rspb_name_col]]), collapse = ", ")))
}

## Load habitat layers
sg_known_path <- "data/processed/seagrass_polygons_harmonised.gpkg"
sm_path <- "data/processed/species/saltmarsh_harmonised.gpkg"
oy_known_path <- "data/processed/species/oyster_known_extent.gpkg"
dune_path <- "outputs/habitats/sand_dunes_habmos_solway.gpkg"

sg_ref <- if (file.exists(sg_known_path))
  st_read(sg_known_path, quiet = TRUE) %>% st_transform(27700)
sm_ref <- if (file.exists(sm_path))
  st_read(sm_path, quiet = TRUE) %>% st_transform(27700)
oy_ref <- if (file.exists(oy_known_path))
  st_read(oy_known_path, quiet = TRUE) %>% st_transform(27700)
dune_ref <- if (file.exists(dune_path))
  st_read(dune_path, quiet = TRUE) %>% st_transform(27700)

## Woodland: saved by Section 5b (HabMoS + HLCM 2022 reclassification)
wood_path <- "data/processed/habitats/woodland_habmos_solway.gpkg"
wood_ref <- NULL
if (file.exists(wood_path)) {
  habmos_wood_saved <- st_read(wood_path, quiet = TRUE)
  if (st_crs(habmos_wood_saved)$epsg != 27700) {
    habmos_wood_saved <- st_transform(habmos_wood_saved, 27700)
  }
  wood_ref <- habmos_wood_saved[grepl("^G1", habmos_wood_saved$HABITAT_CODE), ]
  cat(sprintf("  Native woodland (G1) polygons: %d\n", nrow(wood_ref)))
} else {
  cat("  WARNING: woodland_habmos_solway.gpkg not found — run Section 5b first\n")
}

## Reference site container, filled in one habitat at a time below
ref_sites <- list()

## ---- Helper: clip a habitat layer to the RSPB reserve(s) that contain it ----
## The reference site is the habitat itself, not the whole reserve (which
## includes farmland, grassland, open water etc. with no direct bearing on
## that habitat)
clip_habitat_to_rspb <- function(hab_data, rspb, rspb_name_col) {
  hits <- lengths(st_intersects(rspb, hab_data)) > 0
  if (!any(hits)) return(NULL)
  rspb_with_hab <- rspb[hits, ]
  name_keep <- if (!is.na(rspb_name_col)) rspb_name_col else character(0)
  hab_clipped <- st_intersection(hab_data, rspb_with_hab[, name_keep, drop = FALSE])
  hab_clipped <- st_collection_extract(hab_clipped, "POLYGON")
  hab_clipped[st_is_valid(hab_clipped), ]
}

## ---- RSPB x Saltmarsh ----
ref_sites[["Saltmarsh"]] <- clip_habitat_to_rspb(sm_ref, rspb, rspb_name_col)
cat(sprintf("  Saltmarsh reference: %d polygons clipped to RSPB reserves\n",
            nrow(ref_sites[["Saltmarsh"]])))

## ---- RSPB x Sand dunes ----
ref_sites[["Sand dunes"]] <- clip_habitat_to_rspb(dune_ref, rspb, rspb_name_col)
cat(sprintf("  Sand dune reference: %d polygons clipped to RSPB reserves\n",
            nrow(ref_sites[["Sand dunes"]])))

## ---- RSPB x Woodland ----
ref_sites[["Woodland"]] <- clip_habitat_to_rspb(wood_ref, rspb, rspb_name_col)
cat(sprintf("  Woodland reference: %d polygons clipped to RSPB reserves\n",
            nrow(ref_sites[["Woodland"]])))

## ---- Oyster reference: Loch Ryan known extent ----
## The only site oysters occur at in the SCAMP area, so the entire known
## extent defines "good condition"
ref_sites[["Oyster"]] <- oy_ref
cat(sprintf("  Oyster reference: all known extent (%d polygons)\n", nrow(oy_ref)))


## ============================================================
## DATA-DRIVEN CANDIDATE SITES — Saltmarsh & Sand dunes
## ============================================================
## Terrestrial pressure data is comparatively sparse (only SCIMAP diffuse
## pollution risk — no equivalent to the ~15-20 layer marine pressure
## stack), so these use richness alone as an exploratory check for
## high-value sites the RSPB-derived references might miss. These are
## CANDIDATES to flag for review, not formal reference sites.

## ---- Saltmarsh candidates ----
sm_grid <- rich_terr[lengths(st_intersects(rich_terr, sm_ref)) > 0, ]
sm_rich_q75 <- quantile(sm_grid$species_richness, 0.75, na.rm = TRUE)
sm_candidates <- sm_grid[sm_grid$species_richness >= sm_rich_q75, ]
cat(sprintf("  Saltmarsh candidate cells: %d of %d (richness >= %.0f)\n",
            nrow(sm_candidates), nrow(sm_grid), sm_rich_q75))

## ---- Sand dune candidates ----
dune_grid <- rich_terr[lengths(st_intersects(rich_terr, dune_ref)) > 0, ]
dune_rich_q75 <- quantile(dune_grid$species_richness, 0.75, na.rm = TRUE)
dune_candidates <- dune_grid[dune_grid$species_richness >= dune_rich_q75, ]
cat(sprintf("  Dune candidate cells: %d of %d (richness >= %.0f)\n",
            nrow(dune_candidates), nrow(dune_grid), dune_rich_q75))


## ============================================================
## SEAGRASS REFERENCE — rarefied + interpolated richness surface
## ============================================================
## Raw richness (species count per grid cell) is biased by recording effort
## — a heavily-visited cell accumulates more species regardless of true
## diversity. The existing marine richness grid also only contains cells
## with at least one record, leaving gaps in under-sampled coastal areas
## (e.g. narrow firths) where seagrass is known to occur.
##
## Fix: rarefy richness to a common sample size per cell (Hurlbert's
## rarefaction via vegan), then interpolate (IDW) across the full SMR
## extent at finer resolution to produce a continuous, effort-corrected
## surface with no gaps. Run each step below in turn and inspect before
## moving to the next.

## ---- Step 1: load records + build fine grid across full SMR ----

records_path <- "data/processed/species/all_species_records.gpkg"
stopifnot("Species records not found" = file.exists(records_path))
stopifnot("Cumulative pressure raster not found" = file.exists(cum_path))

cum_r <- rast(cum_path)
all_records <- st_read(records_path, quiet = TRUE)
if (st_crs(all_records)$epsg != 27700) all_records <- st_transform(all_records, 27700)

marine_records <- all_records[lengths(st_intersects(all_records, smr)) > 0, ]
cat(sprintf("  Marine records within SMR: %d\n", nrow(marine_records)))

## Fine grid (500m) covering the FULL SMR extent — downscaled from the
## 5km grid used for the main richness maps, and not limited to cells
## that happen to contain a record
rich_fine_grid <- st_make_grid(smr, cellsize = 500, square = TRUE) %>%
  st_sf(cell_id = seq_along(.), geometry = .)
rich_fine_grid <- rich_fine_grid[lengths(st_intersects(rich_fine_grid, smr)) > 0, ]
cat(sprintf("  Fine grid cells (500m) across SMR: %d\n", nrow(rich_fine_grid)))

## ---- Step 2: species x cell abundance matrix ----

rec_gridded <- st_join(marine_records, rich_fine_grid, join = st_within) %>%
  st_drop_geometry() %>%
  filter(!is.na(cell_id), !is.na(latin_name), latin_name != "")

sp_counts <- rec_gridded %>% count(cell_id, latin_name, name = "n")
sp_wide <- sp_counts %>%
  pivot_wider(id_cols = cell_id, names_from = latin_name,
             values_from = n, values_fill = 0)
cell_ids <- sp_wide$cell_id
sp_mat <- as.matrix(sp_wide[, -1, drop = FALSE])
record_totals <- rowSums(sp_mat)
raw_richness <- rowSums(sp_mat > 0)

cat(sprintf("  Cells with >=1 record: %d | record range: %d to %d\n",
            length(cell_ids), min(record_totals), max(record_totals)))

## ---- Diagnostic: does raw richness just track recording effort? ----
## If cells with more records simply show more species regardless of true
## diversity, raw richness is an effort artefact, not a biodiversity signal
effort_check <- data.frame(record_count = record_totals, raw_richness = raw_richness)
effort_cor <- cor(effort_check$record_count, effort_check$raw_richness,
                  method = "spearman")
cat(sprintf("  Spearman correlation (raw richness vs record count): %.2f\n", effort_cor))

p_effort_bias <- ggplot(effort_check, aes(x = record_count, y = raw_richness)) +
  geom_point(alpha = 0.4, colour = "steelblue") +
  geom_smooth(method = "loess", colour = "firebrick", se = FALSE) +
  labs(title = "Does raw richness just track recording effort?",
       subtitle = sprintf("Spearman r = %.2f — %s", effort_cor,
                          if (effort_cor > 0.5) "strong effort bias" else "weak effort bias"),
       x = "Records per cell", y = "Raw species richness") +
  theme_minimal()

save_plot(p_effort_bias, "priority/L1e0_raw_richness_vs_effort.png", width = 7, height = 6)

## ---- Spatial: raw richness and record count across the fine grid ----
effort_pts <- rich_fine_grid[match(cell_ids, rich_fine_grid$cell_id), ]
effort_pts <- suppressWarnings(st_centroid(effort_pts))
effort_pts$raw_richness <- raw_richness
effort_pts$record_count <- record_totals

p_raw_richness_map <- base_map("Raw species richness (pre-rarefaction)",
  "Species count per 500m cell — before correcting for recording effort") +
  geom_sf(data = effort_pts, aes(colour = raw_richness), size = 1.2) +
  scale_colour_viridis_c(name = "Raw\nrichness", option = "viridis")

save_plot(p_raw_richness_map, "priority/L1e0a_raw_richness_map.png")

p_record_count_map <- base_map("Recording effort",
  "Records per 500m cell — compare against raw richness map above") +
  geom_sf(data = effort_pts, aes(colour = record_count), size = 1.2) +
  scale_colour_viridis_c(name = "Records", option = "inferno", trans = "log1p")

save_plot(p_record_count_map, "priority/L1e0b_record_count_map.png")

## ---- Step 3: rarefy to a common sample size ----
## IMPORTANT: the target must NOT be the dataset minimum, nor a low
## percentile — record counts are extremely right-skewed (a mass of
## cells with only 1-2 incidental records, a handful with thousands), so
## both collapse to a near-useless target (sample=2 caps richness at 2
## for almost every cell). Instead use a FIXED, ecologically sensible
## minimum sample size. Cells below it are excluded from rarefaction
## entirely (NA) rather than distorted — the IDW step (Step 4) then
## fills those gaps from nearby reliable cells, which is exactly what
## it's for.

cat(sprintf("  Record count distribution: min=%d p10=%.0f p25=%.0f median=%.0f max=%d\n",
            min(record_totals), quantile(record_totals, 0.10),
            quantile(record_totals, 0.25), median(record_totals),
            max(record_totals)))

## Retained-cell count at a range of candidate fixed targets, to inspect
## before committing to one
for (candidate in c(5, 10, 15, 20, 25, 50)) {
  n_at <- sum(record_totals >= candidate)
  cat(sprintf("    target=%3d records -> %5d of %d cells retained (%.0f%%)\n",
              candidate, n_at, length(record_totals),
              100 * n_at / length(record_totals)))
}

## ---- Rarefaction curve: find the knee point, don't guess ----
## The first attempt below is WRONG and must not be used for picking a
## target: at each sample size it restricted to a different, shrinking
## set of "eligible" cells (those with >= n records). As n rises, that
## filter keeps only richer, more heavily-recorded cells — so the mean
## climbs almost linearly with n not because richness fails to saturate,
## but because the eligible-cell population itself is shifting toward
## higher-diversity sites. This confounds sampling bias with the very
## thing we're trying to remove it.
##
## Correct approach: fix ONE constant set of well-sampled cells (those
## with enough records to be evaluated at the largest candidate size),
## then rarefy that SAME set at increasing sample sizes. This isolates
## the true within-cell accumulation curve — sampling bias is removed
## by holding the population fixed, not by letting it drift.

curve_sizes <- c(2, 5, 10, 15, 20, 30, 40, 50, 75, 100, 150, 200)
curve_min_records <- max(curve_sizes)
common_cells <- which(record_totals >= curve_min_records)
cat(sprintf("  Fixed cell set for curve: %d cells with >= %d records (same cells at every sample size)\n",
            length(common_cells), curve_min_records))

curve_long <- do.call(rbind, lapply(curve_sizes, function(n) {
  r <- vegan::rarefy(sp_mat[common_cells, , drop = FALSE], sample = n)
  data.frame(cell = common_cells, sample_size = n, richness = as.numeric(r))
}))

curve_summary <- curve_long %>%
  group_by(sample_size) %>%
  summarise(mean_richness = mean(richness), .groups = "drop")

cat("\n  Rarefaction curve on fixed cell set (sample size | mean richness):\n")
print(curve_summary, row.names = FALSE)

p_rarecurve <- ggplot(curve_long, aes(x = sample_size, y = richness, group = cell)) +
  geom_line(alpha = 0.1, colour = "grey40") +
  geom_line(data = curve_summary, aes(x = sample_size, y = mean_richness, group = 1),
            colour = "steelblue", linewidth = 1.3) +
  labs(title = "Rarefaction curve — pick the target from where this bends",
       subtitle = sprintf("Fixed set of %d cells (>= %d records) | thin grey = individual cells, blue = mean",
                          length(common_cells), curve_min_records),
       x = "Sample size (records)", y = "Rarefied richness") +
  theme_minimal()

save_plot(p_rarecurve, "priority/L1e_rarefaction_curve.png", width = 7, height = 6)

## ---- Interpretation of the fixed-cell-set curve ----
## 1. NO PLATEAU: mean richness is still climbing steeply from n=150
##    (~95) to n=200 (~120) — no knee, no flattening within the tested
##    range. Even the best-recorded cells (>= 200 records) are nowhere
##    near having their species pool exhausted. Whatever target we pick,
##    rarefied richness remains a lower bound, not "true" richness — the
##    curve doesn't hand us a clean saturation point the way a textbook
##    example would.
## 2. THE FAN IS THE SIGNAL, NOT THE BEND: near the origin, individual
##    cell lines (thin grey) are bundled tight together — at low n every
##    cell looks similarly (un)diverse because too few records have been
##    drawn to tell rich sites from poor ones. The lines only start
##    visibly separating from roughly n=50-100 onward. That spread is
##    rarefaction actually discriminating a genuinely diverse cell from
##    a genuinely poor one, rather than just reporting effort.
## 3. DECISION: target=50 is a pragmatic compromise, not a mathematically
##    "correct" answer forced by the data. From the retained-cell table
##    above: target=50 keeps 489 cells (roughly where the fan starts
##    opening up) vs target=100's 268 cells or target=200's 143 cells.
##    Fewer anchor points at higher targets would leave the IDW
##    interpolation (Step 4) with bigger gaps to fill in under-sampled
##    coastal areas — exactly the problem this whole exercise is meant
##    to fix. Cells below the target are excluded from rarefaction
##    entirely (not distorted) and get filled by interpolation instead.

rare_target <- 50
keep <- record_totals >= rare_target
cat(sprintf("  Rarefaction target: %d records/cell (fixed) | %d of %d cells retained (%.0f%%)\n",
            rare_target, sum(keep), length(record_totals),
            100 * sum(keep) / length(record_totals)))

rarefied_richness <- vegan::rarefy(sp_mat[keep, , drop = FALSE], sample = rare_target)

rich_pts <- rich_fine_grid[match(cell_ids[keep], rich_fine_grid$cell_id), ]
rich_pts <- suppressWarnings(st_centroid(rich_pts))
rich_pts$rarefied_richness <- as.numeric(rarefied_richness)
rich_pts$raw_richness <- raw_richness[keep]
rich_pts$record_count <- record_totals[keep]

cat(sprintf("  Cells with rarefied richness (input to interpolation): %d\n",
            nrow(rich_pts)))

## ---- Diagnostic: raw vs rarefied richness ----
## Cells that shift most from the 1:1 line are where effort correction
## matters — typically high-record cells whose raw richness was inflated
rare_cor <- cor(rich_pts$raw_richness, rich_pts$rarefied_richness, method = "spearman")
cat(sprintf("  Spearman correlation (raw vs rarefied richness): %.2f\n", rare_cor))

p_rarefaction_effect <- ggplot(st_drop_geometry(rich_pts),
                               aes(x = raw_richness, y = rarefied_richness,
                                   colour = record_count)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(alpha = 0.6) +
  scale_colour_viridis_c(name = "Records\nper cell", option = "inferno") +
  labs(title = "Effect of rarefaction on species richness",
       subtitle = sprintf("Rarefied to %d records/cell | dashed line = no change",
                          rare_target),
       x = "Raw richness (species count)", y = "Rarefied richness") +
  theme_minimal()

save_plot(p_rarefaction_effect, "priority/L1e1_rarefaction_effect.png",
          width = 7, height = 6)

## ---- Spatial: rarefied richness itself ----
p_rarefied_richness_map <- base_map("Rarefied species richness",
  sprintf("Effort-corrected richness (rarefied to %d records/cell) — compare against raw richness map",
          rare_target)) +
  geom_sf(data = rich_pts, aes(colour = rarefied_richness), size = 1.2) +
  scale_colour_viridis_c(name = "Rarefied\nrichness", option = "viridis")

save_plot(p_rarefied_richness_map, "priority/L1e1a_rarefied_richness_map.png")

## ---- Spatial: where does rarefaction change richness the most? ----
rich_pts$rarefaction_shift <- rich_pts$raw_richness - rich_pts$rarefied_richness

p_rarefaction_map <- base_map("Where rarefaction changes the richness picture most",
  "Raw richness minus rarefied richness — high values = effort-inflated cells") +
  geom_sf(data = rich_pts, aes(colour = rarefaction_shift), size = 1.2) +
  scale_colour_viridis_c(name = "Raw −\nrarefied", option = "magma")

save_plot(p_rarefaction_map, "priority/L1e1b_rarefaction_shift_map.png")

## ---- Step 4: IDW interpolation across the full SMR extent ----
## Fills gaps between sparsely-sampled cells so coastal areas with real
## but under-recorded seagrass habitat aren't excluded purely for lack
## of survey effort

## terra::interpolate() builds its prediction grid internally as a plain
## data.frame of x/y columns (not sp/sf), then calls predict(model,
## newdata) — so the model must be built with a `~x+y` locations formula
## on a plain data.frame, NOT an sp/sf object. Passing sp/sf locations
## makes model-fitting work, but predict() then fails on the plain-df
## newdata terra hands it internally, since gstat can't call coordinates()
## on something that was never told to expect a data.frame.
rich_coords <- st_coordinates(rich_pts)
rich_df <- st_drop_geometry(rich_pts)
rich_df$x <- rich_coords[, 1]
rich_df$y <- rich_coords[, 2]

rich_template <- rast(vect(smr), resolution = 300)
idw_model <- gstat::gstat(formula = rarefied_richness ~ 1,
                          locations = ~x + y, data = rich_df, nmax = 12,
                          set = list(idp = 2))
rich_surface <- terra::interpolate(rich_template, idw_model, debug.level = 0)
rich_surface <- rich_surface[["var1.pred"]]
rich_surface <- mask(rich_surface, vect(smr))
names(rich_surface) <- "rarefied_richness"

cat(sprintf("  Interpolated richness surface range: %.2f to %.2f\n",
            global(rich_surface, "min", na.rm = TRUE)[[1]],
            global(rich_surface, "max", na.rm = TRUE)[[1]]))

writeRaster(rich_surface, "outputs/species/rarefied_richness_surface_marine.tif",
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))
cat("  Saved: outputs/species/rarefied_richness_surface_marine.tif\n")

## ---- Step 5: sample richness + pressure within known seagrass footprint ----

sg_fine <- rich_fine_grid[lengths(st_intersects(rich_fine_grid, sg_ref)) > 0, ]
sg_fine_c <- suppressWarnings(st_centroid(sg_fine))
sg_fine$rarefied_richness <- terra::extract(rich_surface, vect(sg_fine_c))[[2]]
sg_fine$cum_pressure <- terra::extract(cum_r, vect(sg_fine_c))[[2]]
sg_fine <- sg_fine[!is.na(sg_fine$rarefied_richness) & !is.na(sg_fine$cum_pressure), ]

cat(sprintf("  Fine cells overlapping known seagrass: %d\n", nrow(sg_fine)))

## ---- Step 6: top 75% richness + bottom 25% pressure -> reference cells ----

rich_q75 <- quantile(sg_fine$rarefied_richness, 0.75, na.rm = TRUE)
pres_q25 <- quantile(sg_fine$cum_pressure, 0.25, na.rm = TRUE)

sg_ref_cells <- sg_fine[sg_fine$rarefied_richness >= rich_q75 &
                        sg_fine$cum_pressure <= pres_q25, ]

if (nrow(sg_ref_cells) == 0) {
  sg_ref_cells <- sg_fine[sg_fine$rarefied_richness >= rich_q75, ]
  cat("  (relaxed: top 25% rarefied richness only — no cells passed both filters)\n")
}

ref_sites[["Seagrass"]] <- sg_ref_cells
cat(sprintf("  Seagrass reference cells: %d (rarefied richness >= %.2f, pressure <= %.2f)\n",
            nrow(sg_ref_cells), rich_q75, pres_q25))

cat(sprintf("\n  Reference site habitats defined: %s\n",
            paste(names(ref_sites), collapse = ", ")))

## ---- L1a: Saltmarsh reference sites ----

p_ref_sm <- base_map("Layer 1a — Saltmarsh Reference Sites",
  sprintf("RSPB saltmarsh within reserve boundary (n=%d) | Candidate cells: top 25%% richness (n=%d, no pressure filter)",
          if (!is.null(ref_sites[["Saltmarsh"]])) nrow(ref_sites[["Saltmarsh"]]) else 0,
          if (!is.null(sm_candidates)) nrow(sm_candidates) else 0)) +
  geom_sf(data = sm_ref, aes(colour = "Known saltmarsh"),
          fill = "grey70", alpha = 0.3, linewidth = 0.3) +
  {if (!is.null(sm_candidates))
    geom_sf(data = sm_candidates, aes(colour = "Data-driven candidate"),
            fill = NA, linewidth = 0.6, linetype = "dashed")} +
  geom_sf(data = ref_sites[["Saltmarsh"]], aes(colour = "RSPB reference"),
          fill = "#66a61e", alpha = 0.5, linewidth = 0.3) +
  scale_colour_manual(name = NULL,
    values = c("Known saltmarsh" = "grey40", "RSPB reference" = "#2d6a1e",
               "Data-driven candidate" = "#e7298a"),
    guide = guide_legend(override.aes = list(
      fill = c("grey70", "#66a61e", NA),
      linetype = c("solid", "solid", "dashed"),
      alpha = c(0.3, 0.5, 1))))

save_plot(p_ref_sm, "priority/L1a_reference_saltmarsh.png")


## ---- L1b: Sand dune reference sites ----

p_ref_dunes <- base_map("Layer 1b — Sand Dune Reference Sites",
  sprintf("RSPB dunes within reserve boundary (n=%d) | Candidate cells: top 25%% richness (n=%d, no pressure filter)",
          if (!is.null(ref_sites[["Sand dunes"]])) nrow(ref_sites[["Sand dunes"]]) else 0,
          if (!is.null(dune_candidates)) nrow(dune_candidates) else 0)) +
  geom_sf(data = dune_ref, aes(colour = "Known dunes"),
          fill = "grey70", alpha = 0.3, linewidth = 0.3) +
  {if (!is.null(dune_candidates))
    geom_sf(data = dune_candidates, aes(colour = "Data-driven candidate"),
            fill = NA, linewidth = 0.6, linetype = "dashed")} +
  geom_sf(data = ref_sites[["Sand dunes"]], aes(colour = "RSPB reference"),
          fill = "#e6ab02", alpha = 0.5, linewidth = 0.3) +
  scale_colour_manual(name = NULL,
    values = c("Known dunes" = "grey40", "RSPB reference" = "#b8860b",
               "Data-driven candidate" = "#e7298a"),
    guide = guide_legend(override.aes = list(
      fill = c("grey70", "#e6ab02", NA),
      linetype = c("solid", "solid", "dashed"),
      alpha = c(0.3, 0.5, 1))))

save_plot(p_ref_dunes, "priority/L1b_reference_sand_dunes.png")


## ---- L1c: Woodland reference sites ----

if (!is.null(ref_sites[["Woodland"]]) && !is.null(wood_ref)) {
  p_ref_wood <- base_map("Layer 1c — Woodland Reference Sites",
    sprintf("RSPB reserves containing native woodland (n=%d)",
            nrow(ref_sites[["Woodland"]]))) +
    geom_sf(data = wood_ref, aes(colour = "Native woodland (G1)"),
            fill = "grey70", alpha = 0.3, linewidth = 0.1) +
    geom_sf(data = ref_sites[["Woodland"]], aes(colour = "RSPB reference"),
            fill = "#1b7837", alpha = 0.4, linewidth = 0.5) +
    scale_colour_manual(name = NULL,
      values = c("Native woodland (G1)" = "grey40", "RSPB reference" = "#0a5020"),
      guide = guide_legend(override.aes = list(
        fill = c("grey70", "#1b7837"), alpha = c(0.3, 0.4))))

  save_plot(p_ref_wood, "priority/L1c_reference_woodland.png")
}


## ---- L1d: Oyster reference sites ----

p_ref_oy <- base_map("Layer 1d — Oyster Reference Sites",
  sprintf("All known extent (%d polygons) — no RSPB reserve",
          if (!is.null(oy_ref)) nrow(oy_ref) else 0)) +
  geom_sf(data = ref_sites[["Oyster"]], aes(colour = "Known extent (reference)"),
          fill = "#a6761d", alpha = 0.5, linewidth = 0.5) +
  scale_colour_manual(name = NULL,
    values = c("Known extent (reference)" = "black"),
    guide = guide_legend(override.aes = list(
      fill = "#a6761d", alpha = 0.5)))

save_plot(p_ref_oy, "priority/L1d_reference_oyster.png")


## ---- L1e: Seagrass reference sites ----

p_ref_sg <- base_map("Layer 1e — Seagrass Reference Sites",
  sprintf("Top 25%% rarefied richness (500m, IDW-interpolated) + bottom 25%% pressure (n=%d)",
          if (!is.null(ref_sites[["Seagrass"]])) nrow(ref_sites[["Seagrass"]]) else 0)) +
  geom_sf(data = sg_ref, aes(colour = "Known seagrass"),
          fill = "grey70", alpha = 0.3, linewidth = 0.3) +
  geom_sf(data = ref_sites[["Seagrass"]], aes(colour = "Reference cells"),
          fill = "#1b9e77", alpha = 0.6, linewidth = 0.5) +
  scale_colour_manual(name = NULL,
    values = c("Known seagrass" = "grey40", "Reference cells" = "black"),
    guide = guide_legend(override.aes = list(
      fill = c("grey70", "#1b9e77"), alpha = c(0.3, 0.6))))

save_plot(p_ref_sg, "priority/L1e_reference_seagrass.png")


## ---- L1e2: Rarefied richness surface (diagnostic) ----

p_rich_surface <- ggplot() +
  geom_sf(data = dg, fill = "grey92", colour = "grey40", linewidth = 0.3) +
  geom_sf(data = smr, fill = NA, colour = "steelblue", linewidth = 0.4, linetype = "dotted") +
  geom_tile(data = rich_surf_df, aes(x = x, y = y, fill = rarefied_richness), alpha = 0.9) +
  scale_fill_viridis_c(name = "Rarefied\nrichness", option = "viridis", na.value = NA) +
  geom_sf(data = sg_ref, fill = NA, colour = "red", linewidth = 0.4) +
  geom_sf(data = scamp_aoi, fill = NA, colour = "black", linewidth = 0.5, linetype = "dashed") +
  coord_sf(crs = 27700) +
  labs(title = "Layer 1e — Rarefied Richness Surface",
       subtitle = "IDW-interpolated (500m grid, 300m raster) — input to seagrass reference selection") +
  theme_minimal()

save_plot(p_rich_surface, "priority/L1e2_rarefied_richness_surface.png")


## ---- Marine beta diversity surface (interpolated from native 5km) ----
## TODO (deferred — revisit later): this interpolates the RAW, un-rarefied
## beta_diversity column as-is. It is NOT yet effort-corrected the way
## richness is. The proper fix (Beck, Holloway & Schwanghart 2013,
## Methods in Ecology & Evolution — "Undersampling and the measurement of
## beta diversity") is to rarefy BEFORE comparing composition: for each
## 5km cell and its queen-contiguous neighbours, randomly subsample both
## down to a common record depth (reuse rare_target = 50 for consistency
## with the richness rarefaction), compute Sorensen dissimilarity from the
## subsampled sets, repeat ~100x and average — there's no closed-form
## shortcut for this the way vegan::rarefy() has for richness, so genuine
## Monte Carlo resampling is required. That also needs the marine
## species-by-cell presence/absence matrix rebuilt at the native 5km grid
## (not saved to disk from 02_themeD_species.R — only the final
## beta_diversity column was kept). Skipped for now — using the raw
## column directly below. Revisit before finalising Layer 3b.

cat("\n-- Marine beta diversity surface (interpolated, NOT yet effort-corrected) --\n")

beta_pts <- rich_marine %>% filter(!is.na(beta_diversity))
beta_pts <- suppressWarnings(st_centroid(beta_pts))

## terra::interpolate() passes a plain data.frame of grid x/y as newdata
## internally — gstat's predict() needs the model built with a ~x+y
## locations formula on a plain data.frame to match, not an sf object
## (same fix as the richness IDW model)
beta_coords <- st_coordinates(beta_pts)
beta_df <- st_drop_geometry(beta_pts)
beta_df$x <- beta_coords[, 1]
beta_df$y <- beta_coords[, 2]

beta_idw_model <- gstat::gstat(formula = beta_diversity ~ 1,
                               locations = ~x + y, data = beta_df, nmax = 12,
                               set = list(idp = 2))
beta_surface <- terra::interpolate(rich_template, beta_idw_model, debug.level = 0)
beta_surface <- beta_surface[["var1.pred"]]
beta_surface <- mask(beta_surface, vect(smr))
names(beta_surface) <- "beta_diversity"

cat(sprintf("  Interpolated beta diversity surface range: %.3f to %.3f\n",
            global(beta_surface, "min", na.rm = TRUE)[[1]],
            global(beta_surface, "max", na.rm = TRUE)[[1]]))

writeRaster(beta_surface, "outputs/species/beta_diversity_surface_marine.tif",
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))
cat("  Saved: outputs/species/beta_diversity_surface_marine.tif\n")

beta_surf_df <- as.data.frame(beta_surface, xy = TRUE)

p_beta_surface <- ggplot() +
  geom_sf(data = dg, fill = "grey92", colour = "grey40", linewidth = 0.3) +
  geom_sf(data = smr, fill = NA, colour = "steelblue", linewidth = 0.4, linetype = "dotted") +
  geom_tile(data = beta_surf_df, aes(x = x, y = y, fill = beta_diversity), alpha = 0.9) +
  scale_fill_viridis_c(name = "Beta\ndiversity", option = "magma", na.value = NA) +
  geom_sf(data = scamp_aoi, fill = NA, colour = "black", linewidth = 0.5, linetype = "dashed") +
  coord_sf(crs = 27700) +
  labs(title = "Marine Beta Diversity Surface (interpolated from native 5km)",
       subtitle = "NOT yet effort-corrected — see TODO note in script") +
  theme_minimal()

save_plot(p_beta_surface, "priority/L3b_beta_diversity_surface.png")


# ============================================================
# TERRESTRIAL (DUNE): RAREFIED RICHNESS + BETA DIVERSITY SURFACES
# ============================================================
# Mirrors the marine rarefaction/interpolation pipeline exactly, applied
# to terrestrial records within the SCAMP AOI instead of marine records
# within smr. Uses scamp_aoi (the actual project boundary), NOT las_dg_bng
# (the wider Dumfries & Galloway council boundary used elsewhere for map
# land-context fill) — las_dg_bng extends well beyond the SCAMP AOI, which
# was letting the interpolation bleed north past the project's own border.

cat("\n-- Terrestrial: rarefied + interpolated richness surface --\n")

terr_records <- all_records[lengths(st_intersects(all_records, scamp_aoi)) > 0, ]
cat(sprintf("  Terrestrial records within SCAMP AOI: %d\n", nrow(terr_records)))

terr_fine_grid <- st_make_grid(scamp_aoi, cellsize = 500, square = TRUE) %>%
  st_sf(cell_id = seq_along(.), geometry = .)
terr_fine_grid <- terr_fine_grid[lengths(st_intersects(terr_fine_grid, scamp_aoi)) > 0, ]
cat(sprintf("  Fine grid cells (500m) across SCAMP AOI: %d\n", nrow(terr_fine_grid)))

terr_rec_gridded <- st_join(terr_records, terr_fine_grid, join = st_within) %>%
  st_drop_geometry() %>%
  filter(!is.na(cell_id), !is.na(latin_name), latin_name != "")

terr_sp_counts <- terr_rec_gridded %>% count(cell_id, latin_name, name = "n")
terr_sp_wide <- terr_sp_counts %>%
  pivot_wider(id_cols = cell_id, names_from = latin_name, values_from = n, values_fill = 0)
terr_cell_ids <- terr_sp_wide$cell_id
terr_sp_mat <- as.matrix(terr_sp_wide[, -1, drop = FALSE])
terr_record_totals <- rowSums(terr_sp_mat)

cat(sprintf("  Cells with >=1 record: %d | record range: %d to %d\n",
            length(terr_cell_ids), min(terr_record_totals), max(terr_record_totals)))

## ---- Rarefaction curve: terrestrial-specific, don't reuse marine's 50 ----
## Terrestrial recording effort (e.g. bird survey coverage) can have a very
## different distribution shape than marine OBIS/GBIF records — same
## diagnostic as the marine curve, but must be re-derived from this data,
## not assumed. See marine section for full reasoning on why this needs a
## FIXED cell set (not a shifting one) and why the curve, not a round
## number, should drive the choice.

for (candidate in c(5, 10, 15, 20, 25, 50)) {
  n_at <- sum(terr_record_totals >= candidate)
  cat(sprintf("    target=%3d records -> %5d of %d cells retained (%.0f%%)\n",
              candidate, n_at, length(terr_record_totals),
              100 * n_at / length(terr_record_totals)))
}

terr_curve_sizes <- c(2, 5, 10, 15, 20, 30, 40, 50, 75, 100, 150, 200)
terr_curve_min_records <- max(terr_curve_sizes)
terr_common_cells <- which(terr_record_totals >= terr_curve_min_records)
cat(sprintf("  Fixed cell set for curve: %d cells with >= %d records (same cells at every sample size)\n",
            length(terr_common_cells), terr_curve_min_records))

terr_curve_long <- do.call(rbind, lapply(terr_curve_sizes, function(n) {
  r <- vegan::rarefy(terr_sp_mat[terr_common_cells, , drop = FALSE], sample = n)
  data.frame(cell = terr_common_cells, sample_size = n, richness = as.numeric(r))
}))

terr_curve_summary <- terr_curve_long %>%
  group_by(sample_size) %>%
  summarise(mean_richness = mean(richness), .groups = "drop")

cat("\n  Terrestrial rarefaction curve on fixed cell set (sample size | mean richness):\n")
print(terr_curve_summary, row.names = FALSE)

p_terr_rarecurve <- ggplot(terr_curve_long, aes(x = sample_size, y = richness, group = cell)) +
  geom_line(alpha = 0.1, colour = "grey40") +
  geom_line(data = terr_curve_summary, aes(x = sample_size, y = mean_richness, group = 1),
            colour = "steelblue", linewidth = 1.3) +
  labs(title = "Terrestrial rarefaction curve — pick the target from where this bends",
       subtitle = sprintf("Fixed set of %d cells (>= %d records) | thin grey = individual cells, blue = mean",
                          length(terr_common_cells), terr_curve_min_records),
       x = "Sample size (records)", y = "Rarefied richness") +
  theme_minimal()

save_plot(p_terr_rarecurve, "priority/L3_terrestrial_rarefaction_curve.png", width = 7, height = 6)

## Chosen from the curve + retention table above: 100 records/cell
terr_rare_target <- 100
terr_keep <- terr_record_totals >= terr_rare_target
cat(sprintf("  Rarefaction target: %d records/cell | %d of %d cells retained (%.0f%%)\n",
            terr_rare_target, sum(terr_keep), length(terr_record_totals),
            100 * sum(terr_keep) / length(terr_record_totals)))

terr_rarefied_richness <- vegan::rarefy(terr_sp_mat[terr_keep, , drop = FALSE], sample = terr_rare_target)

terr_rich_pts <- terr_fine_grid[match(terr_cell_ids[terr_keep], terr_fine_grid$cell_id), ]
terr_rich_pts <- suppressWarnings(st_centroid(terr_rich_pts))
terr_rich_pts$rarefied_richness <- as.numeric(terr_rarefied_richness)
terr_rich_pts$record_count <- terr_record_totals[terr_keep]

cat(sprintf("  Cells with rarefied richness (input to interpolation): %d\n", nrow(terr_rich_pts)))

terr_rich_coords <- st_coordinates(terr_rich_pts)
terr_rich_df <- st_drop_geometry(terr_rich_pts)
terr_rich_df$x <- terr_rich_coords[, 1]
terr_rich_df$y <- terr_rich_coords[, 2]

terr_rich_template <- rast(vect(scamp_aoi), resolution = 300)
terr_idw_model <- gstat::gstat(formula = rarefied_richness ~ 1,
                               locations = ~x + y, data = terr_rich_df, nmax = 12,
                               set = list(idp = 2))
terr_rich_surface <- terra::interpolate(terr_rich_template, terr_idw_model, debug.level = 0)
terr_rich_surface <- terr_rich_surface[["var1.pred"]]
terr_rich_surface <- mask(terr_rich_surface, vect(scamp_aoi))
names(terr_rich_surface) <- "rarefied_richness"

cat(sprintf("  Interpolated terrestrial richness surface range: %.2f to %.2f\n",
            global(terr_rich_surface, "min", na.rm = TRUE)[[1]],
            global(terr_rich_surface, "max", na.rm = TRUE)[[1]]))

writeRaster(terr_rich_surface, "outputs/species/rarefied_richness_surface_terrestrial.tif",
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))
cat("  Saved: outputs/species/rarefied_richness_surface_terrestrial.tif\n")

p_terr_rich_surface <- base_map_fast("Terrestrial Rarefied Richness Surface",
  "IDW-interpolated (500m grid, 300m raster) — dune restoration priority input") +
  geom_spatraster(data = terr_rich_surface, aes(fill = rarefied_richness), alpha = 0.9) +
  scale_fill_viridis_c(name = "Rarefied\nrichness", option = "viridis", na.value = NA) +
  geom_sf(data = dune_ref, fill = NA, colour = "red", linewidth = 0.4)

save_plot(p_terr_rich_surface, "priority/L3_terrestrial_rarefied_richness_surface.png")


## ---- Terrestrial beta diversity surface (interpolated from native 5km) ----
## Same caveat as marine: un-rarefied, effort-bias not yet corrected —
## see the TODO note in the marine beta diversity section above

cat("\n-- Terrestrial beta diversity surface (interpolated, NOT yet effort-corrected) --\n")

terr_beta_pts <- rich_terr %>% filter(!is.na(beta_diversity))
terr_beta_pts <- suppressWarnings(st_centroid(terr_beta_pts))

terr_beta_coords <- st_coordinates(terr_beta_pts)
terr_beta_df <- st_drop_geometry(terr_beta_pts)
terr_beta_df$x <- terr_beta_coords[, 1]
terr_beta_df$y <- terr_beta_coords[, 2]

terr_beta_idw_model <- gstat::gstat(formula = beta_diversity ~ 1,
                                    locations = ~x + y, data = terr_beta_df, nmax = 12,
                                    set = list(idp = 2))
terr_beta_surface <- terra::interpolate(terr_rich_template, terr_beta_idw_model, debug.level = 0)
terr_beta_surface <- terr_beta_surface[["var1.pred"]]
terr_beta_surface <- mask(terr_beta_surface, vect(scamp_aoi))
names(terr_beta_surface) <- "beta_diversity"

cat(sprintf("  Interpolated terrestrial beta diversity surface range: %.3f to %.3f\n",
            global(terr_beta_surface, "min", na.rm = TRUE)[[1]],
            global(terr_beta_surface, "max", na.rm = TRUE)[[1]]))

writeRaster(terr_beta_surface, "outputs/species/beta_diversity_surface_terrestrial.tif",
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))
cat("  Saved: outputs/species/beta_diversity_surface_terrestrial.tif\n")

p_terr_beta_surface <- base_map_fast("Terrestrial Beta Diversity Surface (interpolated from native 5km)",
  "NOT yet effort-corrected — see TODO note in marine section") +
  geom_spatraster(data = terr_beta_surface, aes(fill = beta_diversity), alpha = 0.9) +
  scale_fill_viridis_c(name = "Beta\ndiversity", option = "magma", na.value = NA)

save_plot(p_terr_beta_surface, "priority/L3_terrestrial_beta_diversity_surface.png")

cat("\n-- Terrestrial richness + beta diversity surfaces complete --\n")


## ---- L1f: Combined overview ----

ref_combined <- do.call(rbind, lapply(names(ref_sites), function(hab) {
  st_sf(habitat = hab, geometry = st_geometry(ref_sites[[hab]]))
}))

p_ref_all <- base_map("Layer 1 — All Reference Sites",
  "RSPB reserves matched to habitats + data-driven references") +
  geom_sf(data = ref_combined, aes(fill = habitat),
          colour = "black", alpha = 0.6, linewidth = 0.3) +
  scale_fill_manual(name = "Habitat",
    values = c("Saltmarsh" = "#66a61e", "Sand dunes" = "#e6ab02",
               "Woodland" = "#1b7837", "Oyster" = "#a6761d",
               "Seagrass" = "#1b9e77"))

save_plot(p_ref_all, "priority/L1f_reference_sites_all.png")

cat("\n-- Layers 0 & 1 complete --\n")


# ============================================================
# RESTORATION OPPORTUNITY + CONNECTIVITY (feeds Layer 3a/3b)
# ============================================================

cat("\n-- Restoration opportunity + coastalNet connectivity --\n")

## ---- Seagrass restoration opportunity: potential minus known extent ----
## st_difference is safe here regardless of whether the potential layer
## already excludes known sites (04_physical.R force-includes known
## seagrass into the potential layer) — differencing out known extent
## either strips it back out, or is a no-op if already excluded.

sg_potential_path <- "outputs/habitats/seagrass_potential_solway.gpkg"
stopifnot("Seagrass potential layer not found" = file.exists(sg_potential_path))
sg_potential_poly <- st_read(sg_potential_path, quiet = TRUE) %>% st_transform(27700)

sg_opportunity <- st_difference(st_union(sg_potential_poly), st_union(sg_ref))
sg_opportunity <- st_sf(habitat = "Seagrass", geometry = st_sfc(sg_opportunity, crs = 27700))
sg_opportunity <- st_cast(sg_opportunity, "POLYGON")
cat(sprintf("  Seagrass restoration opportunity polygons: %d\n", nrow(sg_opportunity)))

## ---- Oyster restoration opportunity: potential minus known extent ----

oy_potential_path <- "outputs/habitats/oyster_potential_solway.gpkg"
stopifnot("Oyster potential layer not found" = file.exists(oy_potential_path))
oy_potential_poly <- st_read(oy_potential_path, quiet = TRUE) %>% st_transform(27700)

oy_opportunity <- st_difference(st_union(oy_potential_poly), st_union(oy_ref))
oy_opportunity <- st_sf(habitat = "Oyster", geometry = st_sfc(oy_opportunity, crs = 27700))
oy_opportunity <- st_cast(oy_opportunity, "POLYGON")
cat(sprintf("  Oyster restoration opportunity polygons: %d\n", nrow(oy_opportunity)))

## ---- Join coastalNet connectivity score onto opportunity cells ----
## Hexagon-level scores (max connectivity from any existing source, per
## 04_physical.R Section 10j) — every cell within the same hexagon gets
## the same score, since that's the native resolution of the connectivity
## data (8.45km hexagon edge)

sg_conn_hex_path <- "outputs/physical/seagrass_connectivity_hexagons.gpkg"
oy_conn_hex_path <- "outputs/physical/oyster_connectivity_hexagons.gpkg"
stopifnot("Seagrass connectivity hexagons not found — run 04_physical.R Section 10 first" =
            file.exists(sg_conn_hex_path))
stopifnot("Oyster connectivity hexagons not found — run 04_physical.R Section 10 first" =
            file.exists(oy_conn_hex_path))

sg_conn_hex <- st_read(sg_conn_hex_path, quiet = TRUE) %>% st_transform(27700)
oy_conn_hex <- st_read(oy_conn_hex_path, quiet = TRUE) %>% st_transform(27700)

sg_opportunity <- st_join(sg_opportunity, sg_conn_hex[, "connectivity_score"])
oy_opportunity <- st_join(oy_opportunity, oy_conn_hex[, "connectivity_score"])

cat(sprintf("  Seagrass opportunity cells with connectivity score: %d of %d\n",
            sum(!is.na(sg_opportunity$connectivity_score)), nrow(sg_opportunity)))
cat(sprintf("  Oyster opportunity cells with connectivity score: %d of %d\n",
            sum(!is.na(oy_opportunity$connectivity_score)), nrow(oy_opportunity)))

st_write(sg_opportunity, "outputs/habitats/seagrass_restoration_opportunity.gpkg",
         delete_layer = TRUE, quiet = TRUE)
st_write(oy_opportunity, "outputs/habitats/oyster_restoration_opportunity.gpkg",
         delete_layer = TRUE, quiet = TRUE)
cat("  Saved: outputs/habitats/seagrass_restoration_opportunity.gpkg\n")
cat("  Saved: outputs/habitats/oyster_restoration_opportunity.gpkg\n")

## ---- Plot: restoration opportunity coloured by connectivity score ----

p_sg_opportunity <- base_map_fast("Seagrass Restoration Opportunity x Connectivity",
  "Potential minus known extent, coloured by coastalNet connectivity score") +
  geom_sf(data = sg_ref, fill = "grey60", colour = NA, alpha = 0.5) +
  geom_sf(data = sg_opportunity, aes(fill = connectivity_score),
          colour = NA, alpha = 0.85) +
  scale_fill_viridis_c(name = "Connectivity\nscore", option = "magma",
                       na.value = "grey85")

save_plot(p_sg_opportunity, "priority/L3_seagrass_opportunity_connectivity.png")

p_oy_opportunity <- base_map_fast("Oyster Restoration Opportunity x Connectivity",
  "Potential minus known extent, coloured by coastalNet connectivity score") +
  geom_sf(data = oy_ref, fill = "grey60", colour = NA, alpha = 0.5) +
  geom_sf(data = oy_opportunity, aes(fill = connectivity_score),
          colour = NA, alpha = 0.85) +
  scale_fill_viridis_c(name = "Connectivity\nscore", option = "magma",
                       na.value = "grey85")

save_plot(p_oy_opportunity, "priority/L3_oyster_opportunity_connectivity.png")

cat("\n-- Restoration opportunity + connectivity complete --\n")


# ============================================================
# LAYER 2: PROTECTION PRIORITY
# ============================================================
# High-value cells (per existing whole-seascape richness + community
# similarity, already computed in 02_themeD_species.R) that aren't
# already covered by a formal designation.
#
# Quadrant interpretation (median split):
#   High similarity + high richness = typical, healthy community -> PROTECT
#   High similarity + low richness  = may be misleading -> INVESTIGATE
#   Low similarity  + high richness = unique community -> PROTECT (different)
#   Low similarity  + low richness  = degraded -> RESTORATION CANDIDATE
#
# Pressure count (Layer 0) is marine-only (no equivalent terrestrial
# pressure stack exists) — included as an informational sort column for
# marine cells, NOT used for terrestrial classification.

cat("\n-- Layer 2: Protection priority --\n")

score_protection_priority <- function(grid_sf, pressure_r = NULL, desig_sf) {
  grid_sf <- grid_sf %>% filter(!is.na(species_richness), !is.na(community_similarity))

  sim_med <- median(grid_sf$community_similarity, na.rm = TRUE)
  rich_med <- median(grid_sf$species_richness, na.rm = TRUE)

  grid_sf <- grid_sf %>%
    mutate(
      high_sim = community_similarity >= sim_med,
      high_rich = species_richness >= rich_med,
      value_category = case_when(
        high_sim & high_rich   ~ "Protect (healthy)",
        high_sim & !high_rich  ~ "Investigate",
        !high_sim & high_rich  ~ "Protect (unique)",
        TRUE                   ~ "Restoration candidate"
      )
    )

  if (!is.null(pressure_r)) {
    grid_sf$pressure_count <- terra::extract(pressure_r, vect(grid_sf), fun = "mean", na.rm = TRUE)[[2]]
  } else {
    grid_sf$pressure_count <- NA_real_
  }

  ## Protection gap: not already covered by any existing designation
  grid_sf$protected <- lengths(st_intersects(grid_sf, desig_sf)) > 0
  grid_sf$protection_gap <- !grid_sf$protected

  grid_sf %>%
    mutate(protection_priority = value_category %in% c("Protect (healthy)", "Protect (unique)") &
                                  protection_gap)
}

rich_marine_l2 <- score_protection_priority(rich_marine, pressure_count, desig_all)
cat(sprintf("  Marine — protection priority cells: %d of %d\n",
            sum(rich_marine_l2$protection_priority), nrow(rich_marine_l2)))

rich_terr_l2 <- score_protection_priority(rich_terr, NULL, desig_all)
cat(sprintf("  Terrestrial — protection priority cells: %d of %d (no pressure data available)\n",
            sum(rich_terr_l2$protection_priority), nrow(rich_terr_l2)))

st_write(rich_marine_l2, "outputs/species/protection_priority_marine.gpkg",
         delete_layer = TRUE, quiet = TRUE)
st_write(rich_terr_l2, "outputs/species/protection_priority_terrestrial.gpkg",
         delete_layer = TRUE, quiet = TRUE)
cat("  Saved: outputs/species/protection_priority_marine.gpkg\n")
cat("  Saved: outputs/species/protection_priority_terrestrial.gpkg\n")

value_colours <- c("Protect (healthy)" = "#1b7837", "Protect (unique)" = "#7570b3",
                   "Investigate" = "#e6ab02", "Restoration candidate" = "#d95f02")

## ---- L2a: Marine value category ----

p_l2a <- base_map_fast("Layer 2a — Marine Ecological Value",
  "Community similarity x richness quadrant (median split)") +
  geom_sf(data = rich_marine_l2, aes(fill = value_category), colour = NA, alpha = 0.75) +
  scale_fill_manual(name = "Category", values = value_colours)

save_plot(p_l2a, "priority/L2a_marine_value_category.png")


## ---- L2b: Marine protection priority ----

p_l2b <- base_map_fast("Layer 2b — Marine Protection Priority",
  "High-value cells not covered by SAC/SPA/SSSI/Ramsar/MPA/NNR/LNR/IBA/RSPB") +
  geom_sf(data = rich_marine_l2, aes(fill = protection_priority), colour = NA, alpha = 0.75) +
  scale_fill_manual(name = "Protection\npriority",
    values = c("TRUE" = "#d7191c", "FALSE" = "grey85")) +
  geom_sf(data = desig_all, fill = NA, colour = "black", linewidth = 0.3, linetype = "dotted")

save_plot(p_l2b, "priority/L2b_marine_protection_priority.png")


## ---- L2c: Terrestrial value category ----

p_l2c <- base_map_fast("Layer 2c — Terrestrial Ecological Value",
  "Community similarity x richness quadrant (median split)") +
  geom_sf(data = rich_terr_l2, aes(fill = value_category), colour = NA, alpha = 0.75) +
  scale_fill_manual(name = "Category", values = value_colours)

save_plot(p_l2c, "priority/L2c_terrestrial_value_category.png")


## ---- L2d: Terrestrial protection priority ----

p_l2d <- base_map_fast("Layer 2d — Terrestrial Protection Priority",
  "High-value cells not covered by SAC/SPA/SSSI/Ramsar/MPA/NNR/LNR/IBA/RSPB") +
  geom_sf(data = rich_terr_l2, aes(fill = protection_priority), colour = NA, alpha = 0.75) +
  scale_fill_manual(name = "Protection\npriority",
    values = c("TRUE" = "#d7191c", "FALSE" = "grey85")) +
  geom_sf(data = desig_all, fill = NA, colour = "black", linewidth = 0.3, linetype = "dotted")

save_plot(p_l2d, "priority/L2d_terrestrial_protection_priority.png")

cat("\n-- Layer 2 complete --\n")


# ============================================================
# LAYER 3a: ACTIVE RESTORATION PRIORITY (marine — seagrass & oyster)
# ============================================================
# Score = suitable (0/1) + inverse cumulative pressure (0-1) +
#         connectivity (0-1), summed to a 0-3 scale.
# Higher = physically suitable AND low pressure AND well-connected to an
# existing source population — the ideal active-restoration target.

cat("\n-- Layer 3a: Active restoration priority (marine) --\n")

## ---- Cumulative pressure: normalise 0-1 and invert ----

cum_r <- rast(cum_path)
cum_min <- global(cum_r, "min", na.rm = TRUE)[[1]]
cum_max <- global(cum_r, "max", na.rm = TRUE)[[1]]
pressure_inv <- 1 - ((cum_r - cum_min) / (cum_max - cum_min))
names(pressure_inv) <- "pressure_inv"

cat(sprintf("  Pressure inverse range: %.2f to %.2f (1 = lowest pressure)\n",
            global(pressure_inv, "min", na.rm = TRUE)[[1]],
            global(pressure_inv, "max", na.rm = TRUE)[[1]]))

## ---- Helper: build the composite score for one habitat ----

build_restoration_score <- function(suitable_path, known_sf, conn_hex, habitat_label) {
  suitable_r <- rast(suitable_path)
  names(suitable_r) <- "suitable"

  ## Restoration-opportunity mask: suitable AND not already known extent
  ## (suitable_r alone force-includes known sites — see 04_physical.R —
  ## so the final score must be masked to exclude already-occupied cells;
  ## this isn't a restoration candidate, it already has the habitat)
  known_r <- rasterize(vect(known_sf), suitable_r, field = 1, background = 0)
  opportunity_r <- suitable_r * (1 - known_r)

  pressure_r <- resample(pressure_inv, suitable_r, method = "bilinear")

  conn_r <- rasterize(vect(conn_hex), suitable_r, field = "connectivity_score")
  conn_r[is.na(conn_r)] <- 0
  names(conn_r) <- "connectivity"

  score <- (suitable_r + pressure_r + conn_r) * opportunity_r
  score <- mask(score, vect(smr))
  names(score) <- "restoration_score"

  cat(sprintf("  %s restoration score range: %.2f to %.2f\n",
              habitat_label,
              global(score, "min", na.rm = TRUE)[[1]],
              global(score, "max", na.rm = TRUE)[[1]]))

  list(suitable = suitable_r, opportunity = opportunity_r, pressure = pressure_r,
       connectivity = conn_r, score = score)
}

sg_l3 <- build_restoration_score("outputs/habitats/seagrass_potential_raster.tif",
                                 sg_ref, sg_conn_hex, "Seagrass")
oy_l3 <- build_restoration_score("outputs/habitats/oyster_potential_raster.tif",
                                 oy_ref, oy_conn_hex, "Oyster")

writeRaster(sg_l3$score, "outputs/habitats/seagrass_restoration_score.tif",
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))
writeRaster(oy_l3$score, "outputs/habitats/oyster_restoration_score.tif",
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))
cat("  Saved: outputs/habitats/seagrass_restoration_score.tif\n")
cat("  Saved: outputs/habitats/oyster_restoration_score.tif\n")

## ---- L3a-1: Seagrass restoration score ----
## geom_spatraster() (tidyterra) plots the SpatRaster directly instead of
## going through as.data.frame()+geom_tile() — avoids the tile-width
## auto-detection artifact (visible seam lines) that shows up once NA
## masking creates irregular gaps in the data.frame's x/y spacing

p_sg_l3 <- base_map_fast("Layer 3a — Seagrass Active Restoration Priority",
  "Suitable (0/1) + inverse pressure (0-1) + connectivity (0-1), range 0-3") +
  geom_spatraster(data = sg_l3$score, aes(fill = restoration_score), alpha = 0.85) +
  scale_fill_viridis_c(name = "Restoration\nscore", option = "plasma", na.value = NA) +
  geom_sf(data = sg_ref, fill = NA, colour = "black", linewidth = 0.3)

save_plot(p_sg_l3, "priority/L3a_seagrass_restoration_score.png")


## ---- L3a-2: Oyster restoration score ----

p_oy_l3 <- base_map_fast("Layer 3a — Oyster Active Restoration Priority",
  "Suitable (0/1) + inverse pressure (0-1) + connectivity (0-1), range 0-3") +
  geom_spatraster(data = oy_l3$score, aes(fill = restoration_score), alpha = 0.85) +
  scale_fill_viridis_c(name = "Restoration\nscore", option = "plasma", na.value = NA) +
  geom_sf(data = oy_ref, fill = NA, colour = "black", linewidth = 0.3)

save_plot(p_oy_l3, "priority/L3a_oyster_restoration_score.png")

cat("\n-- Layer 3a complete --\n")


# ============================================================
# LAYER 3b: NATURAL RECOVERY SUPPORT PRIORITY (marine — seagrass & oyster)
# ============================================================
# Score = pressure_sum_score (0-1, peaks at medium pressure) +
#         pressure_count_score (0-1, fewer concurrent pressures = better) +
#         connectivity (0-1) +
#         diversity_score (0-1, distinct + locally depauperate = high)
# Summed to a 0-4 scale, restricted to the restoration-opportunity
# footprint (suitable AND not already occupied by known habitat).
#
# diversity_score = (1 - richness_norm) * (0.5 + 0.5 * beta_norm)
#   High richness always suppresses the score regardless of beta (a
#   thriving cell isn't a recovery candidate). Low richness + high beta
#   (distinct from richer neighbours) = top priority (1). Low richness +
#   low beta (similar to equally-poor neighbours) = still a candidate,
#   not excluded (0.5).
#
# NOTE: beta_diversity here is the un-rarefied interpolated surface built
# earlier (see TODO in that section) — not yet effort-corrected.

cat("\n-- Layer 3b: Natural recovery support priority (marine) --\n")

## ---- Normalise the four component layers to 0-1 ----

cum_norm <- (cum_r - cum_min) / (cum_max - cum_min)

count_min <- global(pressure_count, "min", na.rm = TRUE)[[1]]
count_max <- global(pressure_count, "max", na.rm = TRUE)[[1]]
count_norm <- (pressure_count - count_min) / (count_max - count_min)

rich_min <- global(rich_surface, "min", na.rm = TRUE)[[1]]
rich_max <- global(rich_surface, "max", na.rm = TRUE)[[1]]
rich_norm <- (rich_surface - rich_min) / (rich_max - rich_min)

beta_min <- global(beta_surface, "min", na.rm = TRUE)[[1]]
beta_max <- global(beta_surface, "max", na.rm = TRUE)[[1]]
beta_norm <- (beta_surface - beta_min) / (beta_max - beta_min)

cat(sprintf("  Normalised ranges — pressure sum: %.2f-%.2f | count: %.2f-%.2f | richness: %.2f-%.2f | beta: %.2f-%.2f\n",
            global(cum_norm, "min", na.rm = TRUE)[[1]], global(cum_norm, "max", na.rm = TRUE)[[1]],
            global(count_norm, "min", na.rm = TRUE)[[1]], global(count_norm, "max", na.rm = TRUE)[[1]],
            global(rich_norm, "min", na.rm = TRUE)[[1]], global(rich_norm, "max", na.rm = TRUE)[[1]],
            global(beta_norm, "min", na.rm = TRUE)[[1]], global(beta_norm, "max", na.rm = TRUE)[[1]]))

## ---- Helper: build the composite recovery-support score for one habitat ----

build_recovery_score <- function(suitable_r, known_sf, conn_hex, habitat_label) {
  ## Restoration-opportunity mask: suitable AND not already known extent
  ## (suitable_r alone force-includes known sites — see 04_physical.R —
  ## so it must be masked out here to match "opportunity" = unoccupied)
  known_r <- rasterize(vect(known_sf), suitable_r, field = 1, background = 0)
  opportunity_r <- suitable_r * (1 - known_r)

  pressure_sum_r <- 1 - 2 * abs(resample(cum_norm, suitable_r, method = "bilinear") - 0.5)
  names(pressure_sum_r) <- "pressure_sum_score"

  pressure_count_r <- 1 - resample(count_norm, suitable_r, method = "bilinear")
  names(pressure_count_r) <- "pressure_count_score"

  conn_r <- rasterize(vect(conn_hex), suitable_r, field = "connectivity_score")
  conn_r[is.na(conn_r)] <- 0
  names(conn_r) <- "connectivity"

  rich_r <- resample(rich_norm, suitable_r, method = "bilinear")
  beta_r <- resample(beta_norm, suitable_r, method = "bilinear")
  diversity_r <- (1 - rich_r) * (0.5 + 0.5 * beta_r)
  names(diversity_r) <- "diversity_score"

  score <- (pressure_sum_r + pressure_count_r + conn_r + diversity_r) * opportunity_r
  score <- mask(score, vect(smr))
  names(score) <- "recovery_score"

  cat(sprintf("  %s recovery score range: %.2f to %.2f\n",
              habitat_label,
              global(score, "min", na.rm = TRUE)[[1]],
              global(score, "max", na.rm = TRUE)[[1]]))

  list(opportunity = opportunity_r, pressure_sum = pressure_sum_r,
       pressure_count = pressure_count_r, connectivity = conn_r,
       diversity = diversity_r, score = score)
}

sg_l3b <- build_recovery_score(sg_l3$suitable, sg_ref, sg_conn_hex, "Seagrass")
oy_l3b <- build_recovery_score(oy_l3$suitable, oy_ref, oy_conn_hex, "Oyster")

writeRaster(sg_l3b$score, "outputs/habitats/seagrass_recovery_score.tif",
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))
writeRaster(oy_l3b$score, "outputs/habitats/oyster_recovery_score.tif",
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))
cat("  Saved: outputs/habitats/seagrass_recovery_score.tif\n")
cat("  Saved: outputs/habitats/oyster_recovery_score.tif\n")

## ---- L3b-1: Seagrass recovery score ----

p_sg_l3b <- base_map_fast("Layer 3b — Seagrass Natural Recovery Support Priority",
  "Pressure sum (peaked) + pressure count + connectivity + diversity, range 0-4") +
  geom_spatraster(data = sg_l3b$score, aes(fill = recovery_score), alpha = 0.85) +
  scale_fill_viridis_c(name = "Recovery\nscore", option = "cividis", na.value = NA) +
  geom_sf(data = sg_ref, fill = NA, colour = "black", linewidth = 0.3)

save_plot(p_sg_l3b, "priority/L3b_seagrass_recovery_score.png")


## ---- L3b-2: Oyster recovery score ----

p_oy_l3b <- base_map_fast("Layer 3b — Oyster Natural Recovery Support Priority",
  "Pressure sum (peaked) + pressure count + connectivity + diversity, range 0-4") +
  geom_spatraster(data = oy_l3b$score, aes(fill = recovery_score), alpha = 0.85) +
  scale_fill_viridis_c(name = "Recovery\nscore", option = "cividis", na.value = NA) +
  geom_sf(data = oy_ref, fill = NA, colour = "black", linewidth = 0.3)

save_plot(p_oy_l3b, "priority/L3b_oyster_recovery_score.png")

cat("\n-- Layer 3b complete --\n")


# ============================================================
# LAYER 4: MONITORING PRIORITY (marine — seagrass & oyster)
# ============================================================
# Score = diversity_score (0-1, same formula as 3b: distinct + locally
#         depauperate = high) + similarity_score (0-1, LOW community
#         similarity = more room to detect convergence toward a typical
#         profile) + connectivity (0-1) + recovery_score_norm (0-1,
#         Layer 3b's score rescaled — "links to Layer 3b areas": monitor
#         where pressure reduction could plausibly show uplift).
# Summed to a 0-4 scale, restricted to the same restoration-opportunity
# footprint as 3a/3b — avoids "ceiling effects" from monitoring habitat
# that's already established/healthy, per the original plan.

cat("\n-- Layer 4: Monitoring priority (marine) --\n")

## ---- Marine community similarity surface (interpolated from native 5km) ----
## Same caveat as beta_diversity above: un-rarefied, effort-bias not yet
## corrected — see the TODO note there for what a proper fix would need

sim_pts <- rich_marine %>% filter(!is.na(community_similarity))
sim_pts <- suppressWarnings(st_centroid(sim_pts))

sim_coords <- st_coordinates(sim_pts)
sim_df <- st_drop_geometry(sim_pts)
sim_df$x <- sim_coords[, 1]
sim_df$y <- sim_coords[, 2]

sim_idw_model <- gstat::gstat(formula = community_similarity ~ 1,
                              locations = ~x + y, data = sim_df, nmax = 12,
                              set = list(idp = 2))
sim_surface <- terra::interpolate(rich_template, sim_idw_model, debug.level = 0)
sim_surface <- sim_surface[["var1.pred"]]
sim_surface <- mask(sim_surface, vect(smr))
names(sim_surface) <- "community_similarity"

sim_min <- global(sim_surface, "min", na.rm = TRUE)[[1]]
sim_max <- global(sim_surface, "max", na.rm = TRUE)[[1]]
sim_norm <- (sim_surface - sim_min) / (sim_max - sim_min)

cat(sprintf("  Community similarity surface range: %.3f to %.3f\n", sim_min, sim_max))

writeRaster(sim_surface, "outputs/species/community_similarity_surface_marine.tif",
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))
cat("  Saved: outputs/species/community_similarity_surface_marine.tif\n")

## ---- Helper: build the monitoring-priority score for one habitat ----

build_monitoring_score <- function(l3b_result, habitat_label) {
  similarity_r <- 1 - resample(sim_norm, l3b_result$opportunity, method = "bilinear")
  names(similarity_r) <- "similarity_score"

  recovery_norm_r <- l3b_result$score / 4
  names(recovery_norm_r) <- "recovery_score_norm"

  score <- (l3b_result$diversity + similarity_r + l3b_result$connectivity + recovery_norm_r) *
    l3b_result$opportunity
  score <- mask(score, vect(smr))
  names(score) <- "monitoring_score"

  cat(sprintf("  %s monitoring score range: %.2f to %.2f\n",
              habitat_label,
              global(score, "min", na.rm = TRUE)[[1]],
              global(score, "max", na.rm = TRUE)[[1]]))

  list(similarity = similarity_r, recovery_norm = recovery_norm_r, score = score)
}

sg_l4 <- build_monitoring_score(sg_l3b, "Seagrass")
oy_l4 <- build_monitoring_score(oy_l3b, "Oyster")

writeRaster(sg_l4$score, "outputs/habitats/seagrass_monitoring_score.tif",
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))
writeRaster(oy_l4$score, "outputs/habitats/oyster_monitoring_score.tif",
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))
cat("  Saved: outputs/habitats/seagrass_monitoring_score.tif\n")
cat("  Saved: outputs/habitats/oyster_monitoring_score.tif\n")

## ---- L4-1: Seagrass monitoring score ----

p_sg_l4 <- base_map_fast("Layer 4 — Seagrass Monitoring Priority",
  "Diversity + inverse similarity + connectivity + recovery potential, range 0-4") +
  geom_spatraster(data = sg_l4$score, aes(fill = monitoring_score), alpha = 0.85) +
  scale_fill_viridis_c(name = "Monitoring\nscore", option = "mako", na.value = NA) +
  geom_sf(data = sg_ref, fill = NA, colour = "black", linewidth = 0.3)

save_plot(p_sg_l4, "priority/L4_seagrass_monitoring_score.png")


## ---- L4-2: Oyster monitoring score ----

p_oy_l4 <- base_map_fast("Layer 4 — Oyster Monitoring Priority",
  "Diversity + inverse similarity + connectivity + recovery potential, range 0-4") +
  geom_spatraster(data = oy_l4$score, aes(fill = monitoring_score), alpha = 0.85) +
  scale_fill_viridis_c(name = "Monitoring\nscore", option = "mako", na.value = NA) +
  geom_sf(data = oy_ref, fill = NA, colour = "black", linewidth = 0.3)

save_plot(p_oy_l4, "priority/L4_oyster_monitoring_score.png")

cat("\n-- Layer 4 complete --\n")


# ============================================================
# LAYER 5: COMBINED ACTION MAP (marine — seagrass & oyster)
# ============================================================
# Classifies each cell by whichever of Restore (active) / Restore
# (support) / Monitor scores highest, after rescaling each to a
# comparable 0-1 intensity (dividing by its own max possible sum: 3a/3,
# 3b/4, 4/4).
#
# NOTE: "Protect" is not a class here yet. Layer 2 (Protection Priority)
# operates on a different spatial domain — the whole marine+terrestrial
# 5km grid, not the seagrass/oyster-specific fine suitability rasters
# these other layers use — so folding it in cleanly needs the "how does
# protection actually fit into this" discussion, flagged separately.
# Revisit once that's settled.

cat("\n-- Layer 5: Combined action map (marine) --\n")

build_action_map <- function(l3, l3b, l4, habitat_label) {
  restore_active_norm <- l3$score / 3
  restore_support_norm <- l3b$score / 4
  monitor_norm <- l4$score / 4

  action_stack <- c(restore_active_norm, restore_support_norm, monitor_norm)
  names(action_stack) <- c("Restore (active)", "Restore (support)", "Monitor")

  dominant_idx <- terra::which.max(action_stack)

  action_labels <- data.frame(id = 1:3, action = names(action_stack))
  levels(dominant_idx) <- action_labels
  names(dominant_idx) <- "action"

  ## No signal at all (all three exactly 0) -> mask out rather than an
  ## arbitrary first-place tie
  no_signal <- sum(action_stack) == 0
  dominant_idx <- mask(dominant_idx, no_signal, maskvalue = 1)

  cat(sprintf("  %s action map: %d cells classified\n",
              habitat_label, sum(!is.na(values(dominant_idx)))))

  list(stack = action_stack, action = dominant_idx)
}

sg_l5 <- build_action_map(sg_l3, sg_l3b, sg_l4, "Seagrass")
oy_l5 <- build_action_map(oy_l3, oy_l3b, oy_l4, "Oyster")

writeRaster(sg_l5$action, "outputs/habitats/seagrass_action_map.tif",
            overwrite = TRUE, datatype = "INT1U")
writeRaster(oy_l5$action, "outputs/habitats/oyster_action_map.tif",
            overwrite = TRUE, datatype = "INT1U")
cat("  Saved: outputs/habitats/seagrass_action_map.tif\n")
cat("  Saved: outputs/habitats/oyster_action_map.tif\n")

action_colours <- c("Restore (active)" = "#d95f02",
                    "Restore (support)" = "#e6ab02",
                    "Monitor" = "#1b9e77")

## ---- L5-1: Seagrass combined action map ----

p_sg_l5 <- base_map_fast("Layer 5 — Seagrass Combined Action Map",
  "Dominant action per cell — whichever of restore(active)/restore(support)/monitor scores highest") +
  geom_spatraster(data = sg_l5$action, aes(fill = action), alpha = 0.85) +
  scale_fill_manual(name = "Action", values = action_colours, na.value = NA) +
  geom_sf(data = sg_ref, fill = NA, colour = "black", linewidth = 0.3)

save_plot(p_sg_l5, "priority/L5_seagrass_action_map.png")


## ---- L5-2: Oyster combined action map ----

p_oy_l5 <- base_map_fast("Layer 5 — Oyster Combined Action Map",
  "Dominant action per cell — whichever of restore(active)/restore(support)/monitor scores highest") +
  geom_spatraster(data = oy_l5$action, aes(fill = action), alpha = 0.85) +
  scale_fill_manual(name = "Action", values = action_colours, na.value = NA) +
  geom_sf(data = oy_ref, fill = NA, colour = "black", linewidth = 0.3)

save_plot(p_oy_l5, "priority/L5_oyster_action_map.png")

cat("\n-- Layer 5 complete --\n")


# ============================================================
# SAND DUNES — LAYER 3 (merged, no a/b split) + LAYER 4 + LAYER 5
# ============================================================
# Marine split Layer 3 into 3a (active restoration, LOW pressure) and 3b
# (recovery support, MEDIUM pressure) because pressure level was the
# thing distinguishing them. Dunes have no pressure data at all yet, so
# there's nothing to split on — one merged "restoration_score" layer for
# now. Once dune pressure data exists, split into proper 3a/3b then.
#
# coastalNet connectivity doesn't apply either — dune vegetation
# disperses by wind/rhizome, not ocean currents. Two dune-appropriate
# proxies replace it:
#   dist_score    = 1 - normalised distance to nearest existing dune
#                   (closer to an existing system = better seed/sediment
#                   source — same role connectivity played for marine)
#   builder_score = normalised count of habitat-builder species records
#                   (marram grass etc.) within or within 2km of each
#                   cell — direct evidence dune-building process is
#                   already active nearby
#
# Built as separate inspectable steps — plot and check each one before
# combining, not one giant block.

cat("\n-- Sand dunes: Layer 3 components --\n")

if (!exists("dune_suitable_r")) {
  dune_suitable_r <- rast("outputs/habitats/sand_dune_potential_raster.tif")
  cat("  Loaded dune_suitable_r from disk\n")
}
if (!exists("ind_validated")) {
  ind_validated <- read.csv("outputs/species/indicator_species_validated.csv", stringsAsFactors = FALSE)
  cat("  Loaded ind_validated from disk\n")
}
if (!exists("dune_indicators")) {
  ind_all_tmp <- st_read("outputs/species/indicator_records.gpkg", quiet = TRUE) %>% st_transform(27700)
  dune_indicators <- ind_all_tmp[ind_all_tmp$habitat == "Sand dunes", ]
  cat("  Loaded dune_indicators from disk\n")
}

## ---- Step 1: distance to nearest existing dune ----

dune_dist_r <- terra::distance(terr_rich_template, vect(dune_ref))
dist_min <- global(dune_dist_r, "min", na.rm = TRUE)[[1]]
dist_max <- global(dune_dist_r, "max", na.rm = TRUE)[[1]]
dune_dist_score <- 1 - (dune_dist_r - dist_min) / (dist_max - dist_min)
names(dune_dist_score) <- "dist_score"

cat(sprintf("  Distance-to-existing-dune score range: %.2f to %.2f\n",
            global(dune_dist_score, "min", na.rm = TRUE)[[1]],
            global(dune_dist_score, "max", na.rm = TRUE)[[1]]))

p_dune_dist <- base_map_fast("Sand Dune: Distance-to-Existing Score",
  "1 - normalised distance to nearest known dune system") +
  geom_spatraster(data = dune_dist_score, aes(fill = dist_score), alpha = 0.85) +
  scale_fill_viridis_c(name = "Distance\nscore", option = "viridis", na.value = NA) +
  geom_sf(data = dune_ref, fill = NA, colour = "red", linewidth = 0.4)

save_plot(p_dune_dist, "priority/L3_dune_dist_score.png")


## ---- Step 2: habitat-builder species proximity (within 2km) ----

builder_species <- ind_validated[ind_validated$functional_group == "habitat_builder", ]
dune_builders <- dune_indicators[dune_indicators$latin_name %in% builder_species$latin_name, ]
cat(sprintf("  Dune habitat-builder records: %d\n", nrow(dune_builders)))

dune_builders$one <- 1
builder_count_r <- rasterize(vect(dune_builders), terr_rich_template, field = "one",
                             fun = "sum", background = 0)

builder_focal_matrix <- focalMat(terr_rich_template, d = 2000, type = "circle")
builder_density_r <- focal(builder_count_r, w = builder_focal_matrix, fun = "sum", na.rm = TRUE)

bd_min <- global(builder_density_r, "min", na.rm = TRUE)[[1]]
bd_max <- global(builder_density_r, "max", na.rm = TRUE)[[1]]
dune_builder_score <- (builder_density_r - bd_min) / (bd_max - bd_min)
names(dune_builder_score) <- "builder_score"

cat(sprintf("  Habitat-builder proximity score range: %.2f to %.2f\n",
            global(dune_builder_score, "min", na.rm = TRUE)[[1]],
            global(dune_builder_score, "max", na.rm = TRUE)[[1]]))

p_dune_builder <- base_map_fast("Sand Dune: Habitat-Builder Proximity Score",
  "Normalised count of habitat-builder species records within 2km") +
  geom_spatraster(data = dune_builder_score, aes(fill = builder_score), alpha = 0.85) +
  scale_fill_viridis_c(name = "Builder\nscore", option = "inferno", na.value = NA) +
  geom_sf(data = dune_ref, fill = NA, colour = "black", linewidth = 0.3) +
  geom_sf(data = dune_builders, colour = "cyan", size = 1, alpha = 0.7)

save_plot(p_dune_builder, "priority/L3_dune_builder_score.png")


## ---- Step 3: combine into merged Layer 3 (restoration_score) ----
## diversity_score reuses the richness x beta combination already built
## and inspected in the terrestrial richness/beta section above

terr_rich_norm_min <- global(terr_rich_surface, "min", na.rm = TRUE)[[1]]
terr_rich_norm_max <- global(terr_rich_surface, "max", na.rm = TRUE)[[1]]
terr_rich_norm <- (terr_rich_surface - terr_rich_norm_min) / (terr_rich_norm_max - terr_rich_norm_min)

terr_beta_norm_min <- global(terr_beta_surface, "min", na.rm = TRUE)[[1]]
terr_beta_norm_max <- global(terr_beta_surface, "max", na.rm = TRUE)[[1]]
terr_beta_norm <- (terr_beta_surface - terr_beta_norm_min) / (terr_beta_norm_max - terr_beta_norm_min)

dune_diversity_r <- (1 - terr_rich_norm) * (0.5 + 0.5 * terr_beta_norm)
names(dune_diversity_r) <- "diversity_score"

dist_score_r      <- resample(dune_dist_score, dune_suitable_r, method = "bilinear")
builder_score_r   <- resample(dune_builder_score, dune_suitable_r, method = "bilinear")
diversity_score_r <- resample(dune_diversity_r, dune_suitable_r, method = "bilinear")

dune_l3_score <- (dune_suitable_r + dist_score_r + builder_score_r + diversity_score_r) *
  (dune_suitable_r > 0)
dune_l3_score <- mask(dune_l3_score, vect(scamp_aoi))
names(dune_l3_score) <- "restoration_score"

cat(sprintf("  Dune restoration_score range: %.2f to %.2f (0-4: suitable + dist + builder + diversity)\n",
            global(dune_l3_score, "min", na.rm = TRUE)[[1]],
            global(dune_l3_score, "max", na.rm = TRUE)[[1]]))

writeRaster(dune_l3_score, "outputs/habitats/sand_dune_restoration_score.tif",
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))
cat("  Saved: outputs/habitats/sand_dune_restoration_score.tif\n")

p_dune_l3 <- base_map_fast("Layer 3 — Sand Dune Restoration Priority",
  "Suitable (tier-weighted) + distance + habitat-builder proximity + diversity, range 0-4") +
  geom_spatraster(data = dune_l3_score, aes(fill = restoration_score), alpha = 0.85) +
  scale_fill_viridis_c(name = "Restoration\nscore", option = "plasma", na.value = NA) +
  geom_sf(data = dune_ref, fill = NA, colour = "black", linewidth = 0.3)

save_plot(p_dune_l3, "priority/L3_dune_restoration_score.png")

cat("\n-- Sand dunes Layer 3 complete --\n")


# ============================================================
# CLUSTER SUB-PLOTS: aggregate priority layers, zoomed per named site
# ============================================================
# Same pattern as the original habitat cluster sub-plots earlier in this
# script (e.g. "2e. Seagrass cluster sub-plots" ~line 268) — bbox each
# named cluster, zoom the viewport via coord_sf(xlim/ylim), one PNG per
# cluster. Uses base_map_fast() rather than base_map() (no OSM tile
# fetch) since this produces many plots in a row.

cat("\n-- Cluster sub-plots: aggregate priority layers --\n")

plot_cluster_subplots <- function(boxes, raster, known_sf, value_col, layer_label,
                                  output_prefix, fill_option = "viridis", buffer_m = 500) {
  for (i in seq_len(nrow(boxes))) {
    bb <- st_bbox(st_buffer(boxes[i, ], buffer_m))
    cluster_name <- boxes$name[i]
    fname_clean <- gsub("[^a-zA-Z0-9]", "_", tolower(cluster_name))

    ## Crop to the cluster extent FIRST, then plot the small cropped piece
    ## — handing geom_spatraster() the whole-SMR raster and relying on
    ## coord_sf(xlim/ylim) to zoom in afterward triggered a bogus
    ## "16777216 Tb" allocation error on one cluster (float precision
    ## edge case in the resample/crop-window size mismatch)
    raster_crop <- crop(raster, ext(bb["xmin"], bb["xmax"], bb["ymin"], bb["ymax"]))

    p_cl <- base_map_fast(sprintf("%s — %s", layer_label, cluster_name)) +
      geom_spatraster(data = raster_crop, aes(fill = .data[[value_col]]), alpha = 0.85) +
      scale_fill_viridis_c(name = value_col, option = fill_option, na.value = NA) +
      geom_sf(data = known_sf, fill = NA, colour = "black", linewidth = 0.4) +
      coord_sf(xlim = c(bb["xmin"], bb["xmax"]), ylim = c(bb["ymin"], bb["ymax"]), crs = 27700)

    save_plot(p_cl, sprintf("%s_%02d_%s.png", output_prefix, i, fname_clean),
              width = 8, height = 8)
  }
  cat(sprintf("  %d cluster sub-plots saved: %s\n", nrow(boxes), layer_label))
}

plot_cluster_action_subplots <- function(boxes, raster, known_sf, layer_label, output_prefix,
                                         buffer_m = 500) {
  for (i in seq_len(nrow(boxes))) {
    bb <- st_bbox(st_buffer(boxes[i, ], buffer_m))
    cluster_name <- boxes$name[i]
    fname_clean <- gsub("[^a-zA-Z0-9]", "_", tolower(cluster_name))

    raster_crop <- crop(raster, ext(bb["xmin"], bb["xmax"], bb["ymin"], bb["ymax"]))

    p_cl <- base_map_fast(sprintf("%s — %s", layer_label, cluster_name)) +
      geom_spatraster(data = raster_crop, aes(fill = action), alpha = 0.85) +
      scale_fill_manual(name = "Action", values = action_colours, na.value = NA) +
      geom_sf(data = known_sf, fill = NA, colour = "black", linewidth = 0.4) +
      coord_sf(xlim = c(bb["xmin"], bb["xmax"]), ylim = c(bb["ymin"], bb["ymax"]), crs = 27700)

    save_plot(p_cl, sprintf("%s_%02d_%s.png", output_prefix, i, fname_clean),
              width = 8, height = 8)
  }
  cat(sprintf("  %d cluster action sub-plots saved: %s\n", nrow(boxes), layer_label))
}

## ---- Seagrass ----

plot_cluster_subplots(sg_boxes, sg_l3$score, sg_ref, "restoration_score",
                      "Seagrass Restoration Priority", "priority/L3a_seagrass_cluster", "plasma")
plot_cluster_subplots(sg_boxes, sg_l3b$score, sg_ref, "recovery_score",
                      "Seagrass Recovery Support Priority", "priority/L3b_seagrass_cluster", "cividis")
plot_cluster_subplots(sg_boxes, sg_l4$score, sg_ref, "monitoring_score",
                      "Seagrass Monitoring Priority", "priority/L4_seagrass_cluster", "mako")
plot_cluster_action_subplots(sg_boxes, sg_l5$action, sg_ref,
                             "Seagrass Combined Action Map", "priority/L5_seagrass_cluster")

## ---- Oyster ----

plot_cluster_subplots(oy_boxes, oy_l3$score, oy_ref, "restoration_score",
                      "Oyster Restoration Priority", "priority/L3a_oyster_cluster", "plasma")
plot_cluster_subplots(oy_boxes, oy_l3b$score, oy_ref, "recovery_score",
                      "Oyster Recovery Support Priority", "priority/L3b_oyster_cluster", "cividis")
plot_cluster_subplots(oy_boxes, oy_l4$score, oy_ref, "monitoring_score",
                      "Oyster Monitoring Priority", "priority/L4_oyster_cluster", "mako")
plot_cluster_action_subplots(oy_boxes, oy_l5$action, oy_ref,
                             "Oyster Combined Action Map", "priority/L5_oyster_cluster")

## ---- Sand dune (only Layer 3 exists so far — 4/5 pending review) ----

plot_cluster_subplots(dune_boxes, dune_l3_score, dune_ref, "restoration_score",
                      "Sand Dune Restoration Priority", "priority/L3_dune_cluster", "plasma")

cat("\n-- Cluster sub-plots complete --\n")


# ============================================================
# PRIORITY MAPS — PLAN (remaining layers)
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


# ============================================================
# HYDROBASINS SUB-CATCHMENT DELINEATION
# ============================================================
# HydroBASINS as an alternative/complementary catchment delineation to
# the existing CEFAS catchments, for prioritisation. 12 levels available
# (1 = coarse, 12 = very granular) — levels 7, 9, 10 and 11 were tried
# and compared against level 12; level 12 gave the most useful
# sub-catchment granularity and was adopted as the working layer.

cat("\n-- HydroBASINS sub-catchment delineation (level 12) --\n")

library(mapview)

hybas_path <- "data/raw/physical/hybas_lake_eu_lev01-12_v1c/hybas_lake_eu_lev12_v1c.shp"
stopifnot("HydroBASINS level 12 shapefile not found" = file.exists(hybas_path))

hybas_lv12 <- st_read(hybas_path, quiet = TRUE) %>% st_transform(27700)
cat(sprintf("  HydroBASINS level 12 loaded: %d features\n", nrow(hybas_lv12)))

## Bbox crop first (cheap) to cut down feature count before the more
## expensive st_intersection topology operation
las_dg_bbox <- st_as_sfc(st_bbox(las_dg_bng))
hybas_lv12_crop <- st_crop(hybas_lv12, las_dg_bbox) %>% st_make_valid()
cat(sprintf("  After bbox crop: %d features\n", nrow(hybas_lv12_crop)))

hybas_lv12_dg <- st_intersection(hybas_lv12_crop, las_dg_bng) %>% st_make_valid()
cat(sprintf("  After st_intersection with D&G boundary: %d features\n", nrow(hybas_lv12_dg)))

st_write(hybas_lv12_dg, "outputs/catchments/hybas_lev12_dg.gpkg",
         delete_layer = TRUE, quiet = TRUE)
cat("  Saved: outputs/catchments/hybas_lev12_dg.gpkg\n")

## ---- Load existing layers for comparison ----

current_catchments <- st_read("outputs/catchments/cefas_catchments_crop.gpkg", quiet = TRUE) %>%
  st_transform(27700)
cat(sprintf("  Current (CEFAS) catchments: %d features\n", nrow(current_catchments)))

os_rivers <- st_read("outputs/catchments/rivers_with_accumulated_risk.gpkg", quiet = TRUE) %>%
  st_transform(27700)
os_rivers_aoi <- st_intersection(os_rivers, scamp_aoi) %>% st_make_valid()
cat(sprintf("  OS rivers masked to SCAMP AOI: %d of %d features\n",
            nrow(os_rivers_aoi), nrow(os_rivers)))

## ---- Simple mapview comparison ----

mapview(hybas_lv12_dg, color = "purple", col.regions = "purple", alpha.regions = 0.1,
        lwd = 2, layer.name = "HydroBASINS Lvl 12") +
  mapview(current_catchments, color = "steelblue", col.regions = "lightblue",
          alpha.regions = 0.2, lwd = 1, layer.name = "Current (CEFAS) Catchments") +
  mapview(las_dg_bng, color = "black", col.regions = NA, lwd = 2, alpha.regions = 0,
          layer.name = "D&G Council Boundary") +
  mapview(scamp_aoi, color = "red", col.regions = NA, lwd = 2, alpha.regions = 0,
          layer.name = "SCAMP AOI") +
  mapview(os_rivers_aoi, color = "darkblue", lwd = 1, layer.name = "OS Rivers (SCAMP AOI)")

cat("\n-- HydroBASINS level 12 delineation complete --\n")


# ============================================================
# KNOWN SITES + RESTORATION POTENTIAL — ALL FIVE HABITATS
# ============================================================
# For sharing with D&G Council / Technical Advisory Group: a single
# high-quality static map instead of the leaflet (too many layers,
# several not loading reliably). Eight classes: known extent for all
# five habitats, plus potential/restoration-opportunity extent for the
# three habitats that have a suitability model (seagrass, oyster, sand
# dune — saltmarsh and woodland don't have one yet).
#
# "Potential" here means restoration OPPORTUNITY (potential minus known
# extent), not raw suitability — so known and potential polygons for the
# same habitat never overlap on the map. Seagrass/oyster opportunity
# layers were already built earlier in this script; sand dune's
# restoration-potential layer already excludes known dune extent at
# source (built via st_difference in 03_habitats.R).

cat("\n-- Known sites + restoration potential map (all habitats) --\n")

## ---- Load/verify each layer, with disk fallback ----

if (!exists("sg_known")) {
  sg_known <- st_read("data/processed/seagrass_polygons_harmonised.gpkg", quiet = TRUE) %>%
    st_transform(27700)
}
if (!exists("sm_known")) {
  sm_known <- st_read("data/processed/species/saltmarsh_harmonised.gpkg", quiet = TRUE) %>%
    st_transform(27700)
}
if (!exists("oy_known_sf")) {
  oy_known_sf <- st_read("data/processed/species/oyster_known_extent.gpkg", quiet = TRUE) %>%
    st_transform(27700) %>% st_make_valid()
}
if (!exists("dune_known")) {
  dune_known <- st_read("outputs/habitats/sand_dunes_habmos_solway.gpkg", quiet = TRUE) %>%
    st_transform(27700)
}
## Woodland: derived from the full HabMoS-corrected HLCM raster, filtered
## to woodland (G-prefixed) EUNIS classes, rather than the separately
## HabMoS-derived polygon set — per instruction to use HLCM as the base
## classification (it covers all land-use types, not just woodland,
## should this map need other classes later) and filter out non-woodland
if (!exists("hlcm_int")) {
  hlcm_int <- rast("data/processed/habitats/hlcm_habmos_corrected.tif")
  cat("  Loaded hlcm_int (HabMoS-corrected HLCM) from disk\n")
}
if (!exists("lc_lookup")) {
  lc_lookup <- read.csv("data/processed/habitats/hlcm_lc_lookup.csv", stringsAsFactors = FALSE)
  cat("  Loaded lc_lookup from disk\n")
}

woodland_lc_values <- lc_lookup$lc_value[grepl("^G", lc_lookup$eunis_code)]
hlcm_wood_mask <- ifel(hlcm_int %in% woodland_lc_values, 1, NA)
habmos_wood <- as.polygons(hlcm_wood_mask) %>% st_as_sf() %>% st_transform(27700) %>%
  st_make_valid()
cat(sprintf("  Woodland extent from HLCM (filtered to %d woodland classes): %d polygons\n",
            length(woodland_lc_values), nrow(habmos_wood)))

if (!exists("sg_opportunity")) {
  sg_potential_poly <- st_read("outputs/habitats/seagrass_potential_solway.gpkg", quiet = TRUE) %>%
    st_transform(27700)
  sg_opportunity <- st_difference(st_union(sg_potential_poly), st_union(sg_known))
  sg_opportunity <- st_sf(geometry = st_sfc(sg_opportunity, crs = 27700)) %>% st_cast("POLYGON")
}
if (!exists("oy_opportunity")) {
  oy_potential_poly <- st_read("outputs/habitats/oyster_potential_solway.gpkg", quiet = TRUE) %>%
    st_transform(27700)
  oy_opportunity <- st_difference(st_union(oy_potential_poly), st_union(oy_known_sf))
  oy_opportunity <- st_sf(geometry = st_sfc(oy_opportunity, crs = 27700)) %>% st_cast("POLYGON")
}
if (!exists("dune_restore")) {
  dune_restore <- st_read("outputs/habitats/sand_dune_restoration_potential.gpkg", quiet = TRUE) %>%
    st_transform(27700)
}

## ---- Mask every layer to scamp_aoi explicitly, then combine ----
## Don't trust each source script's own clipping extent — e.g. habmos_wood
## is clipped to a 10km BUFFER around scamp_aoi (Section 5b, for HLCM
## cropping headroom), not scamp_aoi itself, so left unmasked it would
## show woodland up to 10km outside the actual project boundary. Applying
## one uniform mask here is safer than auditing each source separately.

mask_to_aoi <- function(x, label) {
  x <- st_make_valid(x)
  out <- suppressWarnings(st_intersection(x, scamp_aoi))
  cat(sprintf("  %s: %d features before mask -> %d after\n", label, nrow(x), nrow(out)))
  out
}

sg_known_aoi      <- mask_to_aoi(sg_known, "Seagrass known")
sm_known_aoi      <- mask_to_aoi(sm_known, "Saltmarsh known")
oy_known_aoi      <- mask_to_aoi(oy_known_sf, "Oyster known")
dune_known_aoi    <- mask_to_aoi(dune_known, "Sand dune known")
habmos_wood_aoi   <- mask_to_aoi(habmos_wood, "Woodland known")
sg_opportunity_aoi  <- mask_to_aoi(sg_opportunity, "Seagrass potential")
oy_opportunity_aoi  <- mask_to_aoi(oy_opportunity, "Oyster potential")
dune_restore_aoi    <- mask_to_aoi(dune_restore, "Sand dune potential")

## ---- Combine into one sf with a single "class" column ----

known_sites_all <- rbind(
  st_sf(class = "Seagrass (known)",   geometry = st_geometry(sg_known_aoi)),
  st_sf(class = "Saltmarsh (known)",  geometry = st_geometry(sm_known_aoi)),
  st_sf(class = "Oyster (known)",     geometry = st_geometry(oy_known_aoi)),
  st_sf(class = "Sand dune (known)",  geometry = st_geometry(dune_known_aoi)),
  st_sf(class = "Woodland (known)",   geometry = st_geometry(habmos_wood_aoi)),
  st_sf(class = "Seagrass (potential)",  geometry = st_geometry(sg_opportunity_aoi)),
  st_sf(class = "Oyster (potential)",    geometry = st_geometry(oy_opportunity_aoi)),
  st_sf(class = "Sand dune (potential)", geometry = st_geometry(dune_restore_aoi))
)
known_sites_all$class <- factor(known_sites_all$class, levels = c(
  "Seagrass (known)", "Seagrass (potential)",
  "Saltmarsh (known)",
  "Oyster (known)", "Oyster (potential)",
  "Sand dune (known)", "Sand dune (potential)",
  "Woodland (known)"
))

cat(sprintf("  Combined known/potential layer: %d features across %d classes\n",
            nrow(known_sites_all), nlevels(known_sites_all$class)))
print(table(known_sites_all$class))

## ---- Distinct 8-colour palette (Brewer Dark2 — designed for exactly
## this: maximum distinguishability between adjacent categorical fills) ----

known_sites_pal <- c(
  "Seagrass (known)"      = "#1b9e77",
  "Seagrass (potential)"  = "#d95f02",
  "Saltmarsh (known)"     = "#7570b3",
  "Oyster (known)"        = "#e7298a",
  "Oyster (potential)"    = "#66a61e",
  "Sand dune (known)"     = "#e6ab02",
  "Sand dune (potential)" = "#a6761d",
  "Woodland (known)"      = "#666666"
)

## ---- Historical oyster locations (red triangles) ----

if (!exists("oy_hist")) {
  oy_hist_path <- "data/processed/species/oyster_historical_locations.gpkg"
  oy_hist <- if (file.exists(oy_hist_path)) {
    st_read(oy_hist_path, quiet = TRUE) %>% st_transform(27700)
  } else {
    cat("  WARNING: oyster_historical_locations.gpkg not found — skipping\n")
    st_sf(geometry = st_sfc(crs = 27700))
  }
}
oy_hist_aoi <- if (nrow(oy_hist) > 0) mask_to_aoi(oy_hist, "Oyster historical") else oy_hist

## ---- Plot on CartoDB (Positron) basemap ----

p_known_sites <- base_map("Solway Firth - SCAMP Priority Habitats & Restoration Potential",
                          "Seagrass, saltmarsh, native oyster, sand dune, woodland",
                          tile_type = "cartolight") +
  geom_sf(data = known_sites_all, aes(fill = class), colour = NA, alpha = 0.85) +
  scale_fill_manual(values = known_sites_pal, name = NULL,
                    guide = guide_legend(override.aes = list(alpha = 1))) +
  {if (nrow(oy_hist_aoi) > 0)
    geom_sf(data = oy_hist_aoi, aes(colour = "Historical oyster sites"),
            shape = 17, size = 2.5, alpha = 0.8)} +
  {if (nrow(oy_hist_aoi) > 0)
    scale_colour_manual(values = c("Historical oyster sites" = "#de2d26"), name = NULL,
                        guide = guide_legend(override.aes = list(shape = 17, size = 3, alpha = 1)))} +
  theme(legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15))

save_plot(p_known_sites, "priority/known_sites_all_habitats.png", width = 14, height = 11)

st_write(known_sites_all, "outputs/habitats/known_sites_all_habitats.gpkg",
         delete_layer = TRUE, quiet = TRUE)
cat("  Saved: outputs/habitats/known_sites_all_habitats.gpkg\n")

cat("\n-- Known sites + restoration potential map complete --\n")


# ============================================================
# LAND USE INTENSITY (SCIMAP RUNOFF RISK) PER HYDROBASINS LEVEL 12
# ============================================================
# Reuses the exact SCIMAP-derived land-cover -> runoff source risk
# weighting already established in 05_localthreats.R (Step 1), applied to
# our own HabMoS-corrected HLCM raster — but NOW also replicates that
# script's grassland intensity split first: James Hutton Land Capability
# for Agriculture (LCA) data divides plain E1/E2/E3 grassland into four
# risk tiers based on agricultural capability class (LCA 1-3 = likely
# improved pasture -> highest runoff risk, down to LCA 6-7 = likely
# semi-natural/rough grazing -> lowest). Without this, all grassland gets
# one generic weight regardless of how intensively it's actually managed
# — a real precision loss for a runoff-risk metric specifically.
#
# Also fixes a related gap: our own G4.F reclassification (Section 5b)
# writes lc_value 33 into the raster but never added a matching row to
# lc_lookup, so those pixels would silently classify to NA below.
#
# Area-weighted mean per HydroBASINS level 12 catchment via
# terra::extract(weights=TRUE) — the built-in mechanism for fraction-of-
# cell-covered weighting where catchment boundaries don't align to raster
# cell edges.

cat("\n-- Land use intensity (SCIMAP runoff risk) per HydroBASINS lvl 12 --\n")

if (!exists("hlcm_int")) {
  hlcm_int <- rast("data/processed/habitats/hlcm_habmos_corrected.tif")
  cat("  Loaded hlcm_int from disk\n")
}
if (!exists("lc_lookup")) {
  lc_lookup <- read.csv("data/processed/habitats/hlcm_lc_lookup.csv", stringsAsFactors = FALSE)
  cat("  Loaded lc_lookup from disk\n")
}
if (!exists("hybas_lv12_dg")) {
  hybas_lv12_dg <- st_read("outputs/catchments/hybas_lev12_dg.gpkg", quiet = TRUE) %>%
    st_transform(27700)
  cat("  Loaded hybas_lv12_dg from disk\n")
}

## ---- Grassland intensity split via James Hutton LCA (replicates
## 05_localthreats.R Step 4b-iii) ----

cat("\n  Loading James Hutton LCA 1:50k...\n")
lca_path <- "data/raw/physical/Hutton_LCA50K_OpenData/LCA_50K.shp"
stopifnot("Hutton LCA shapefile not found" = file.exists(lca_path))

lca <- st_read(lca_path, quiet = TRUE)
if (st_crs(lca)$epsg != 27700) lca <- st_transform(lca, 27700)

lca_crop <- st_crop(lca, st_bbox(hlcm_int)) %>% st_make_valid()
lca_crop <- lca_crop[lca_crop$Water == "No", ]
cat(sprintf("  LCA features after crop (excl. water): %d\n", nrow(lca_crop)))

e1_val <- lc_lookup$lc_value[lc_lookup$eunis_code == "E1"]
e2_val <- lc_lookup$lc_value[lc_lookup$eunis_code == "E2"]
e3_val <- lc_lookup$lc_value[lc_lookup$eunis_code == "E3"]
grass_vals <- c(e1_val, e2_val, e3_val)

lca_class <- floor(lca_crop$LCCODE)
hlcm_reclass <- hlcm_int

## LCA 1-3 -> 29, LCA 4 -> 30, LCA 5 -> 31, LCA 6-7 -> 32
lca_groups <- list(
  list(range = 1:3, new_val = 29, label = "LCA 1-3"),
  list(range = 4,   new_val = 30, label = "LCA 4"),
  list(range = 5,   new_val = 31, label = "LCA 5"),
  list(range = 6:7, new_val = 32, label = "LCA 6-7")
)

for (g in lca_groups) {
  lca_sub <- lca_crop[lca_class %in% g$range, ]
  if (nrow(lca_sub) == 0) next
  ex <- terra::extract(hlcm_reclass, vect(lca_sub), cells = TRUE)
  hits <- ex$cell[ex[[2]] %in% grass_vals]
  if (length(hits) > 0) hlcm_reclass[hits] <- g$new_val
  cat(sprintf("  Grassland pixels reclassified (%s): %d\n", g$label, length(hits)))
}

## Extend lookup with the new grassland-intensity codes + the G4.F=33
## row our own Section 5b never added
lc_lookup_reclass <- rbind(
  lc_lookup,
  data.frame(
    lc_value = c(29, 30, 31, 32, 33),
    eunis_code = c("E.lca13", "E.lca4", "E.lca5", "E.lca67", "G4.F"),
    eunis_desc = c("Likely improved pasture (LCA 1-3)",
                  "Likely improved pasture, marginal (LCA 4)",
                  "Possible improved or rough grazing (LCA 5)",
                  "Likely semi-natural/rough grazing (LCA 6-7)",
                  "Mixed forestry plantation (HabMoS G4.F)"),
    stringsAsFactors = FALSE
  )
)

writeRaster(hlcm_reclass, "data/processed/habitats/hlcm_habmos_lca_reclass.tif",
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))
write.csv(lc_lookup_reclass, "data/processed/habitats/hlcm_lc_lookup_reclass.csv",
          row.names = FALSE)
cat("  Saved: data/processed/habitats/hlcm_habmos_lca_reclass.tif\n")
cat("  Saved: data/processed/habitats/hlcm_lc_lookup_reclass.csv\n")

## ---- SCIMAP source risk weighting (reused from 05_localthreats.R Step 1) ----

risk_map <- c(
  "I1" = 1.00, "E.lca13" = 0.30, "E.lca4" = 0.30, "E.lca5" = 0.20,
  "E.lca67" = 0.15, "E1" = 0.15, "E2" = 0.15, "E3" = 0.15, "E4" = 0.05,
  "E5" = 0.10, "G1" = 0.05, "G3.F" = 0.10, "G3.4" = 0.05, "G4" = 0.05,
  "G4.F" = 0.05, "G5" = 0.10, "F4" = 0.05, "F3" = 0.05, "F2" = 0.05,
  "F9" = 0.05, "D1" = 0.05, "D2" = 0.05, "D4" = 0.05, "J" = 0.01,
  "A2.5" = 0, "B1" = 0, "B2" = 0, "B3" = 0, "C" = 0, "O" = 0,
  "OW" = 0.10, "H2" = 0, "H3" = 0
)

lc_lookup_reclass$risk_weight <- risk_map[lc_lookup_reclass$eunis_code]
cat(sprintf("  EUNIS classes matched to a risk weight: %d of %d\n",
            sum(!is.na(lc_lookup_reclass$risk_weight)), nrow(lc_lookup_reclass)))
cat("  Unmatched EUNIS codes (will be NA in source_risk):\n")
print(lc_lookup_reclass$eunis_code[is.na(lc_lookup_reclass$risk_weight)])

## ---- Apply weighting to the grassland-intensity-reclassified raster ----

rcl_risk <- as.matrix(lc_lookup_reclass[!is.na(lc_lookup_reclass$risk_weight),
                                        c("lc_value", "risk_weight")])
source_risk <- terra::classify(hlcm_reclass, rcl_risk, others = NA)
names(source_risk) <- "risk_weight"

cat(sprintf("  Source risk raster: %d x %d cells | range %.2f to %.2f\n",
            ncol(source_risk), nrow(source_risk),
            global(source_risk, "min", na.rm = TRUE)[[1]],
            global(source_risk, "max", na.rm = TRUE)[[1]]))

writeRaster(source_risk, "outputs/catchments/hlcm_source_runoff_risk.tif",
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))
cat("  Saved: outputs/catchments/hlcm_source_runoff_risk.tif\n")

## ---- Area-weighted mean per HydroBASINS level 12 catchment ----

hybas_lv12_dg <- st_make_valid(hybas_lv12_dg)
catchment_risk_extract <- terra::extract(source_risk, vect(hybas_lv12_dg),
                                         fun = mean, na.rm = TRUE, weights = TRUE, ID = TRUE)

hybas_lv12_dg$mean_runoff_risk <- catchment_risk_extract$risk_weight

cat(sprintf("  Catchments with a runoff risk value: %d of %d\n",
            sum(!is.na(hybas_lv12_dg$mean_runoff_risk)), nrow(hybas_lv12_dg)))
cat(sprintf("  Mean runoff risk range across catchments: %.3f to %.3f\n",
            min(hybas_lv12_dg$mean_runoff_risk, na.rm = TRUE),
            max(hybas_lv12_dg$mean_runoff_risk, na.rm = TRUE)))

st_write(hybas_lv12_dg, "outputs/catchments/hybas_lev12_runoff_risk.gpkg",
         delete_layer = TRUE, quiet = TRUE)
cat("  Saved: outputs/catchments/hybas_lev12_runoff_risk.gpkg\n")

## ---- Plot in the same style as the known-sites habitat map above ----

p_catchment_risk <- base_map("Diffuse runoff risk (Land Cover) by sub-catchment",
                             "area-weighted mean of land-cover-based runoff source risk",
                             tile_type = "cartolight") +
  geom_sf(data = hybas_lv12_dg, aes(fill = mean_runoff_risk), colour = "grey40", linewidth = 0.1) +
  scale_fill_viridis_c(name = "Mean\nrunoff risk", option = "plasma", na.value = "grey85") +
  ## Redraw the AOI boundary on top — base_map() draws it early, so the
  ## catchment fill layer above would otherwise cover it
  geom_sf(data = scamp_aoi, fill = NA, colour = "black", linewidth = 0.6,
          linetype = "dashed") +
  theme(legend.text = element_text(size = 13),
        legend.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 13))

save_plot(p_catchment_risk, "priority/hybas12_runoff_risk.png", width = 14, height = 11)

cat("\n-- Land use intensity per catchment complete --\n")


# ============================================================
# AOI EXPLORATION MAPS
# ============================================================
# Multi-layer spatial planning map for user-defined AOI bboxes.
# Run sequentially — inspect each plot in the R graphics device
# before deciding whether to save.

cat("\n\n========== AOI EXPLORATION MAPS ==========\n")

# ---- Load AOI bboxes ----
aois_main <- st_read("data/processed/boundaries/aois_possibole/input.shp", quiet = TRUE) %>%
  st_transform(27700)
cat(sprintf("  Main AOIs loaded: %d features\n", nrow(aois_main)))

# Confirmed bay names, in upload order: 1 = Loch Ryan, 2 = Fleet Bay,
# 3 = Kirkcudbright Bay, 4 = Auchencairn Bay.
aois_main$bay_name <- c("Loch Ryan", "Fleet Bay", "Kirkcudbright Bay", "Auchencairn Bay")[seq_len(nrow(aois_main))]

# ---- Load Luce Bay / Wigtown Bay bboxes ----
# Two new boxes, order not stated by the user. Luce Bay sits west of
# Wigtown Bay on the Solway coast, so sort by centroid easting (west to
# east) to assign names rather than assuming upload order. VERIFY against
# the printed table below before trusting downstream maps.
aois_lw <- st_read("data/processed/boundaries/aois_luce_wigtown/input.shp", quiet = TRUE) %>%
  st_transform(27700)
cat(sprintf("  Luce/Wigtown AOIs loaded: %d features\n", nrow(aois_lw)))

lw_easting <- st_coordinates(st_centroid(aois_lw))[, 1]
aois_lw <- aois_lw[order(lw_easting), ]
aois_lw$bay_name <- c("Luce Bay", "Wigtown Bay")[seq_len(nrow(aois_lw))]

cat("  Luce/Wigtown assignment (verify centroid positions look right):\n")
print(st_coordinates(st_centroid(aois_lw)) %>%
        cbind(bay_name = aois_lw$bay_name))

# ---- Combine all bays ----
common_cols <- intersect(names(aois_main), names(aois_lw))
aois <- rbind(aois_main[, common_cols], aois_lw[, common_cols])
aois <- aois[order(match(aois$bay_name,
                         c("Loch Ryan", "Luce Bay", "Wigtown Bay", "Fleet Bay",
                           "Kirkcudbright Bay", "Auchencairn Bay"))), ]
rownames(aois) <- NULL

cat(sprintf("  Combined AOIs: %d features\n", nrow(aois)))
print(aois[, "bay_name"])

st_write(aois, "data/processed/boundaries/aois_named.gpkg", delete_layer = TRUE, quiet = TRUE)
cat("  Saved named AOIs: data/processed/boundaries/aois_named.gpkg\n")

# ---- Overview map: cumulative pressure + catchment risk + rivers + AOI labels ----

if (!exists("rivers_in_catch"))
  rivers_in_catch <- st_read("outputs/catchments/rivers_with_accumulated_risk_hydrosheds.gpkg", quiet = TRUE)
if (!exists("hybas_risk"))
  hybas_risk <- st_read("outputs/catchments/catchments_with_risk.gpkg", quiet = TRUE)

cum_fine_path <- "outputs/threats/cumulative_pressure_with_pollution.tif"
if (!file.exists(cum_fine_path))
  cum_fine_path <- "outputs/threats/cumulative_pressure_with_pollution_1km.tif"
cum_pressure <- rast(cum_fine_path)

riv_plot_all <- rivers_in_catch[!is.na(rivers_in_catch$accum_risk_raw) &
                                  rivers_in_catch$accum_risk_raw > 0, ]
riv_plot_all <- st_intersection(riv_plot_all, scamp_aoi)

p_overview <- base_map("Solway Overview — Cumulative Pressure, Catchment Runoff Risk & AOIs",
                       tile_type = "cartolight", zoom = 10) +
  tidyterra::geom_spatraster(data = cum_pressure, alpha = 0.6) +
  scale_fill_viridis_c(name = "Cumulative\nmarine pressure",
                       na.value = "transparent", option = "magma") +
  ggnewscale::new_scale_fill() +
  geom_sf(data = hybas_risk, aes(fill = mean_diffuse_risk),
          colour = "grey50", linewidth = 0.15, alpha = 0.4) +
  scale_fill_gradient(name = "Catchment\nmean runoff risk",
                      low = "#f5deb3", high = "#a50026", na.value = "transparent",
                      labels = function(x) ifelse(x == min(x), "Low", ifelse(x == max(x), "High", ""))) +
  ggnewscale::new_scale_colour() +
  geom_sf(data = riv_plot_all, aes(colour = accum_risk_raw), linewidth = 0.8) +
  scale_colour_gradient(name = "River\naccum. risk",
                        low = "steelblue", high = "orangered", trans = "log1p",
                        breaks = range(riv_plot_all$accum_risk_raw, na.rm = TRUE),
                        labels = c("Low", "High")) +
  geom_sf(data = aois, fill = NA, colour = "white", linewidth = 1.2) +
  geom_sf(data = aois, fill = NA, colour = "black", linewidth = 0.6) +
  theme(legend.position = "right", legend.box = "vertical",
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14),
        axis.text = element_text(size = 12))
print(p_overview)

# ---- Load data layers (disk fallback for anything not in session) ----

if (!exists("sg_known"))
  sg_known <- st_read("data/processed/seagrass_polygons_harmonised.gpkg", quiet = TRUE) %>%
    st_transform(27700)
if (!exists("sm_known"))
  sm_known <- st_read("data/processed/species/saltmarsh_harmonised.gpkg", quiet = TRUE) %>%
    st_transform(27700)
if (!exists("oy_known_sf"))
  oy_known_sf <- st_read("data/processed/species/oyster_known_extent.gpkg", quiet = TRUE) %>%
    st_transform(27700) %>% st_make_valid()
if (!exists("dune_known"))
  dune_known <- st_read("outputs/habitats/sand_dunes_habmos_solway.gpkg", quiet = TRUE) %>%
    st_transform(27700)

if (!exists("sg_opportunity")) {
  sg_pot <- st_read("outputs/habitats/seagrass_potential_solway.gpkg", quiet = TRUE) %>%
    st_transform(27700)
  sg_opportunity <- st_difference(st_union(sg_pot), st_union(sg_known))
  sg_opportunity <- st_sf(geometry = st_sfc(sg_opportunity, crs = 27700)) %>% st_cast("POLYGON")
}
if (!exists("oy_opportunity")) {
  oy_pot <- st_read("outputs/habitats/oyster_potential_solway.gpkg", quiet = TRUE) %>%
    st_transform(27700)
  oy_opportunity <- st_difference(st_union(oy_pot), st_union(oy_known_sf))
  oy_opportunity <- st_sf(geometry = st_sfc(oy_opportunity, crs = 27700)) %>% st_cast("POLYGON")
}
if (!exists("dune_restore"))
  dune_restore <- st_read("outputs/habitats/sand_dune_restoration_potential.gpkg", quiet = TRUE) %>%
    st_transform(27700)

if (!exists("rivers_in_catch"))
  rivers_in_catch <- st_read("outputs/catchments/rivers_with_accumulated_risk_hydrosheds.gpkg", quiet = TRUE)
if (!exists("hybas_risk"))
  hybas_risk <- st_read("outputs/catchments/catchments_with_risk.gpkg", quiet = TRUE)
if (!exists("pour_points"))
  pour_points <- st_read("outputs/catchments/coastal_pour_points.gpkg", quiet = TRUE)

cum_fine_path <- "outputs/threats/cumulative_pressure_with_pollution.tif"
if (!file.exists(cum_fine_path))
  cum_fine_path <- "outputs/threats/cumulative_pressure_with_pollution_1km.tif"
cum_pressure <- rast(cum_fine_path)
cat(sprintf("  Cumulative pressure: %s (%d x %d)\n",
            cum_fine_path, ncol(cum_pressure), nrow(cum_pressure)))

if (!exists("vec_threats_aoi"))
  vec_threats_aoi <- st_read("outputs/threats/vector_threats_unified.gpkg", quiet = TRUE)

cat("  All layers loaded.\n")

# ---- Helper: clip features to an AOI extent ----
clip_to_aoi <- function(sf_obj, aoi_buf) {
  hits <- lengths(st_intersects(sf_obj, aoi_buf)) > 0
  if (sum(hits) == 0) return(NULL)
  sf_obj[hits, ]
}

# ---- AOI 3: Kirkcudbright ----

aoi_i <- aois[3, ]
aoi_bbox <- st_bbox(aoi_i)
aoi_buf <- st_buffer(aoi_i, 500)
cat(sprintf("\n---- AOI 3: Kirkcudbright ----\n"))

# Woodland from HLCM (crop to AOI, filter to G-prefixed EUNIS classes)
if (!exists("hlcm_int"))
  hlcm_int <- rast("data/processed/habitats/hlcm_habmos_corrected.tif")
if (!exists("lc_lookup"))
  lc_lookup <- read.csv("data/processed/habitats/hlcm_lc_lookup.csv", stringsAsFactors = FALSE)

woodland_lc_values <- lc_lookup$lc_value[grepl("^G", lc_lookup$eunis_code)]
hlcm_aoi <- crop(hlcm_int, ext(vect(aoi_buf)))
hlcm_wood_aoi <- ifel(hlcm_aoi %in% woodland_lc_values, 1, NA)
wood_polys <- as.polygons(hlcm_wood_aoi) %>% st_as_sf() %>% st_transform(27700) %>% st_make_valid()

# Clip all layers to AOI
sg_k  <- clip_to_aoi(sg_known, aoi_buf)
sm_k  <- clip_to_aoi(sm_known, aoi_buf)
oy_k  <- clip_to_aoi(oy_known_sf, aoi_buf)
dn_k  <- clip_to_aoi(dune_known, aoi_buf)
sg_p  <- clip_to_aoi(sg_opportunity, aoi_buf)
oy_p  <- clip_to_aoi(oy_opportunity, aoi_buf)
dn_p  <- clip_to_aoi(dune_restore, aoi_buf)
riv_aoi   <- clip_to_aoi(rivers_in_catch, aoi_buf)
hybas_aoi <- clip_to_aoi(hybas_risk, aoi_buf)
vt_aoi_i  <- clip_to_aoi(vec_threats_aoi, aoi_buf)

cum_aoi <- tryCatch(
  crop(cum_pressure, ext(vect(st_transform(aoi_buf, crs(cum_pressure))))) %>%
    mask(vect(st_transform(aoi_buf, crs(cum_pressure)))),
  error = function(e) NULL)
if (!is.null(cum_aoi)) cum_aoi <- project(cum_aoi, "EPSG:27700", method = "bilinear")

riv_plot <- if (!is.null(riv_aoi) && "accum_risk_raw" %in% names(riv_aoi))
  riv_aoi[!is.na(riv_aoi$accum_risk_raw) & riv_aoi$accum_risk_raw > 0, ] else NULL

vt_pts <- vt_pgs <- NULL
if (!is.null(vt_aoi_i)) {
  vt_pts <- vt_aoi_i[st_geometry_type(vt_aoi_i) %in% c("POINT", "MULTIPOINT"), ]
  vt_pgs <- vt_aoi_i[st_geometry_type(vt_aoi_i) %in% c("POLYGON", "MULTIPOLYGON", "GEOMETRYCOLLECTION"), ]
}

# ---- Biodiversity surfaces (load from disk) ----
if (!exists("rich_surface"))
  rich_surface <- rast("outputs/species/rarefied_richness_surface_marine.tif")
if (!exists("beta_surface"))
  beta_surface <- rast("outputs/species/beta_diversity_surface_marine.tif")

# Combine known habitats for panel B
hab_known_list <- list()
if (!is.null(sg_k)) hab_known_list[["Seagrass"]]  <- st_sf(habitat = "Seagrass", geometry = st_geometry(sg_k))
if (!is.null(sm_k)) hab_known_list[["Saltmarsh"]]  <- st_sf(habitat = "Saltmarsh", geometry = st_geometry(sm_k))
if (!is.null(oy_k)) hab_known_list[["Oyster"]]     <- st_sf(habitat = "Oyster", geometry = st_geometry(oy_k))
if (!is.null(dn_k)) hab_known_list[["Sand dune"]]  <- st_sf(habitat = "Sand dune", geometry = st_geometry(dn_k))
if (nrow(wood_polys) > 0) hab_known_list[["Woodland"]] <- st_sf(habitat = "Woodland", geometry = st_geometry(wood_polys))
hab_known_sf <- if (length(hab_known_list) > 0) bind_rows(hab_known_list) else NULL

hab_pot_list <- list()
if (!is.null(sg_p)) hab_pot_list[["Seagrass"]] <- st_sf(habitat = "Seagrass", geometry = st_geometry(sg_p))
if (!is.null(oy_p)) hab_pot_list[["Oyster"]]   <- st_sf(habitat = "Oyster", geometry = st_geometry(oy_p))
if (!is.null(dn_p)) hab_pot_list[["Sand dune"]] <- st_sf(habitat = "Sand dune", geometry = st_geometry(dn_p))
hab_pot_sf <- if (length(hab_pot_list) > 0) bind_rows(hab_pot_list) else NULL

hab_fill_pal <- c("Seagrass" = "#2ca02c", "Saltmarsh" = "#8c564b", "Oyster" = "#9467bd",
                  "Sand dune" = "#e6ab02", "Woodland" = "#0b3d0b")

# Historical oyster locations
if (!exists("oy_hist")) {
  oy_hist_path <- "data/processed/species/oyster_historical_locations.gpkg"
  oy_hist <- if (file.exists(oy_hist_path))
    st_read(oy_hist_path, quiet = TRUE) %>% st_transform(27700) else st_sf(geometry = st_sfc())
}
oy_hist_aoi <- clip_to_aoi(oy_hist, aoi_buf)

# Shared coord_sf for all panels
aoi_coord <- coord_sf(xlim = c(aoi_bbox["xmin"], aoi_bbox["xmax"]),
                      ylim = c(aoi_bbox["ymin"], aoi_bbox["ymax"]), crs = 27700)
aoi_base <- annotation_map_tile(type = "cartolight", zoom = 12, progress = "none", quiet = TRUE)
panel_theme <- theme_scamp() + theme(legend.position = "right",
                                      legend.key.size = unit(0.6, "cm"),
                                      legend.justification = "top",
                                      legend.box.just = "top",
                                      plot.title = element_text(size = 15, face = "bold"),
                                      plot.margin = margin(2, 2, 2, 2))

# ---- Panel A: Threats ----
p_threats <- ggplot() + aoi_base

if (!is.null(cum_aoi))
  p_threats <- p_threats + tidyterra::geom_spatraster(data = cum_aoi, alpha = 0.6) +
    scale_fill_viridis_c(name = "Cumulative\npressure", na.value = "transparent", option = "magma")

if (!is.null(hybas_aoi) && "mean_diffuse_risk" %in% names(hybas_aoi))
  p_threats <- p_threats + ggnewscale::new_scale_fill() +
    geom_sf(data = hybas_aoi, aes(fill = mean_diffuse_risk),
            colour = "grey50", linewidth = 0.15, alpha = 0.35) +
    scale_fill_gradient(name = "Catchment\nmean risk", low = "#f5deb3", high = "#a50026",
                        na.value = "transparent",
                        guide = guide_colourbar(title.position = "top"),
                        labels = function(x) ifelse(x == min(x), "Low", ifelse(x == max(x), "High", "")))

if (!is.null(riv_plot) && nrow(riv_plot) > 0)
  p_threats <- p_threats + ggnewscale::new_scale_colour() +
    geom_sf(data = riv_plot, aes(colour = accum_risk_raw), linewidth = 1.2) +
    scale_colour_gradient(name = "River\naccum. risk", low = "steelblue", high = "orangered",
                          trans = "log1p",
                          breaks = range(riv_plot$accum_risk_raw, na.rm = TRUE),
                          labels = c("Low", "High"))

if (!is.null(vt_aoi_i) && nrow(vt_aoi_i) > 0) {
  vt_aoi_i$threat_label <- gsub("_", " ", tools::toTitleCase(vt_aoi_i$threat_type))
  vt_pts_l <- vt_aoi_i[st_geometry_type(vt_aoi_i) %in% c("POINT", "MULTIPOINT"), ]
  vt_pgs_l <- vt_aoi_i[st_geometry_type(vt_aoi_i) %in% c("POLYGON", "MULTIPOLYGON", "GEOMETRYCOLLECTION"), ]
  if (nrow(vt_pgs_l) > 0)
    p_threats <- p_threats + geom_sf(data = vt_pgs_l, aes(shape = threat_label),
                                     fill = NA, colour = "black", linewidth = 0.8)
  if (nrow(vt_pts_l) > 0)
    p_threats <- p_threats + geom_sf(data = vt_pts_l, aes(shape = threat_label),
                                     colour = "black", size = 3, stroke = 1.2) +
      scale_shape_manual(name = "Vector threats",
                         values = c("Dredge spoil" = 0, "Dredging" = 1,
                                    "Dumped munitions" = 2, "Ports" = 5,
                                    "Waste discharge" = 6, "Ww wrecks" = 4),
                         guide = guide_legend(override.aes = list(size = 3.5),
                                              keywidth = unit(0.6, "cm"),
                                              keyheight = unit(0.5, "cm")))
}

p_threats <- p_threats +
  geom_sf(data = smr, fill = NA, colour = "black", linewidth = 0.5) +
  geom_sf(data = aoi_i, fill = NA, colour = "black", linewidth = 0.6, linetype = "dashed") +
  aoi_coord + labs(title = "A. Threats & Pressures") + panel_theme

# ---- Panel B: Habitats ----
p_habs <- ggplot() + aoi_base

if (!is.null(hab_known_sf))
  p_habs <- p_habs + geom_sf(data = hab_known_sf, aes(fill = habitat),
                               colour = "grey20", linewidth = 0.4, alpha = 0.6) +
    scale_fill_manual(name = "Known habitat", values = hab_fill_pal)

if (!is.null(hab_pot_sf))
  p_habs <- p_habs + ggnewscale::new_scale_fill() +
    geom_sf(data = hab_pot_sf, aes(fill = habitat),
            colour = "black", linewidth = 0.6, alpha = 0.2, linetype = "longdash") +
    scale_fill_manual(name = "Potential habitat", values = hab_fill_pal)

if (!is.null(riv_plot) && nrow(riv_plot) > 0)
  p_habs <- p_habs + geom_sf(data = riv_plot, colour = "steelblue", linewidth = 0.6)

if (!is.null(oy_hist_aoi) && nrow(oy_hist_aoi) > 0)
  p_habs <- p_habs + geom_sf(data = oy_hist_aoi, colour = "#de2d26",
                               shape = 17, size = 2.5, alpha = 0.7) +
    annotate("text", x = mean(st_coordinates(oy_hist_aoi)[, 1]),
             y = min(st_coordinates(oy_hist_aoi)[, 2]) - 800,
             label = "Historical\noyster beds", colour = "#de2d26", size = 2.5,
             fontface = "italic")

p_habs <- p_habs +
  geom_sf(data = smr, fill = NA, colour = "black", linewidth = 0.5) +
  geom_sf(data = aoi_i, fill = NA, colour = "black", linewidth = 0.6, linetype = "dashed") +
  aoi_coord + labs(title = "B. Known & Potential Habitats") + panel_theme

# ---- Panel C: Species Richness ----
rich_aoi <- crop(rich_surface, ext(vect(st_transform(aoi_buf, crs(rich_surface)))))
rich_aoi_bng <- project(rich_aoi, "EPSG:27700", method = "bilinear")

p_rich <- ggplot() + aoi_base +
  tidyterra::geom_spatraster(data = rich_aoi_bng, alpha = 0.75) +
  scale_fill_viridis_c(name = "Rarefied\nrichness", na.value = "transparent", option = "mako",
                       direction = -1)
if (!is.null(hybas_aoi))
  p_rich <- p_rich + geom_sf(data = hybas_aoi, fill = NA, colour = "grey60",
                              linewidth = 0.15, linetype = "dotted")
if (!is.null(riv_plot) && nrow(riv_plot) > 0)
  p_rich <- p_rich + geom_sf(data = riv_plot, colour = "steelblue", linewidth = 0.4, alpha = 0.5)
p_rich <- p_rich +
  geom_sf(data = smr, fill = NA, colour = "black", linewidth = 0.5) +
  geom_sf(data = aoi_i, fill = NA, colour = "black", linewidth = 0.6, linetype = "dashed") +
  aoi_coord + labs(title = "C. Species Richness") + panel_theme

# ---- Panel D: Beta Diversity ----
beta_aoi <- crop(beta_surface, ext(vect(st_transform(aoi_buf, crs(beta_surface)))))
beta_aoi_bng <- project(beta_aoi, "EPSG:27700", method = "bilinear")

p_beta <- ggplot() + aoi_base +
  tidyterra::geom_spatraster(data = beta_aoi_bng, alpha = 0.75) +
  scale_fill_viridis_c(name = "Beta\ndiversity", na.value = "transparent", option = "inferno")
if (!is.null(hybas_aoi))
  p_beta <- p_beta + geom_sf(data = hybas_aoi, fill = NA, colour = "grey60",
                               linewidth = 0.15, linetype = "dotted")
if (!is.null(riv_plot) && nrow(riv_plot) > 0)
  p_beta <- p_beta + geom_sf(data = riv_plot, colour = "steelblue", linewidth = 0.4, alpha = 0.5)
p_beta <- p_beta +
  geom_sf(data = smr, fill = NA, colour = "black", linewidth = 0.5) +
  geom_sf(data = aoi_i, fill = NA, colour = "black", linewidth = 0.6, linetype = "dashed") +
  aoi_coord + labs(title = "D. Beta Diversity") + panel_theme

# ---- Combine four panels ----
p3_combo <- (p_threats + p_habs) / (p_rich + p_beta) +
  patchwork::plot_layout(widths = c(1, 1)) +
  patchwork::plot_annotation(
    title = "Kirkcudbright Bay — Spatial Planning Overview",
    theme = theme(plot.title = element_text(face = "bold", size = 15))
  ) &
  theme(plot.margin = margin(4, 4, 4, 4))
print(p3_combo)

# ---- Export individual panels for PowerPoint assembly ----
ppt_dir <- "figures/ppt_panels"
if (!dir.exists(ppt_dir)) dir.create(ppt_dir, recursive = TRUE)

ggsave(file.path(ppt_dir, "overview_pressure_catchments.png"),
       p_overview, width = 14, height = 11, dpi = 300, bg = "white")
ggsave(file.path(ppt_dir, "aoi3_a_threats.png"),
       p_threats, width = 8, height = 7, dpi = 300, bg = "white")
ggsave(file.path(ppt_dir, "aoi3_b_habitats.png"),
       p_habs, width = 8, height = 7, dpi = 300, bg = "white")
ggsave(file.path(ppt_dir, "aoi3_c_richness.png"),
       p_rich, width = 8, height = 7, dpi = 300, bg = "white")
ggsave(file.path(ppt_dir, "aoi3_d_beta_diversity.png"),
       p_beta, width = 8, height = 7, dpi = 300, bg = "white")
ggsave(file.path(ppt_dir, "aoi3_combined.png"),
       p3_combo, width = 16, height = 14, dpi = 300, bg = "white")
cat(sprintf("  Exported %d panels to %s/\n", 6, ppt_dir))

# ============================================================
# BATCH EXPORT: FOUR PANELS PER AOI
# ============================================================

hab_fill_pal <- c("Seagrass" = "#2ca02c", "Saltmarsh" = "#8c564b", "Oyster" = "#9467bd",
                  "Sand dune" = "#e6ab02", "Woodland" = "#0b3d0b")

for (idx in seq_len(nrow(aois))) {

  aoi_label <- aois$bay_name[idx]
  cat(sprintf("\n==== %s (AOI %d of %d) ====\n", aoi_label, idx, nrow(aois)))

  aoi_dir <- file.path("figures/ppt_panels", aoi_label)
  if (!dir.exists(aoi_dir)) dir.create(aoi_dir, recursive = TRUE)

  aoi_i  <- aois[idx, ]
  aoi_bbox <- st_bbox(aoi_i)
  aoi_buf  <- st_buffer(aoi_i, 500)

  # ---- Clip all layers ----
  sg_k  <- clip_to_aoi(sg_known, aoi_buf)
  sm_k  <- clip_to_aoi(sm_known, aoi_buf)
  oy_k  <- clip_to_aoi(oy_known_sf, aoi_buf)
  dn_k  <- clip_to_aoi(dune_known, aoi_buf)
  sg_p  <- clip_to_aoi(sg_opportunity, aoi_buf)
  oy_p  <- clip_to_aoi(oy_opportunity, aoi_buf)
  dn_p  <- clip_to_aoi(dune_restore, aoi_buf)
  riv_aoi   <- clip_to_aoi(rivers_in_catch, aoi_buf)
  hybas_aoi <- clip_to_aoi(hybas_risk, aoi_buf)
  vt_aoi_i  <- clip_to_aoi(vec_threats_aoi, aoi_buf)
  oy_hist_aoi <- clip_to_aoi(oy_hist, aoi_buf)

  cum_aoi <- tryCatch(
    crop(cum_pressure, ext(vect(st_transform(aoi_buf, crs(cum_pressure))))) %>%
      mask(vect(st_transform(aoi_buf, crs(cum_pressure)))),
    error = function(e) NULL)
  if (!is.null(cum_aoi)) cum_aoi <- project(cum_aoi, "EPSG:27700", method = "bilinear")

  riv_plot <- if (!is.null(riv_aoi) && "accum_risk_raw" %in% names(riv_aoi))
    riv_aoi[!is.na(riv_aoi$accum_risk_raw) & riv_aoi$accum_risk_raw > 0, ] else NULL

  # Woodland from HLCM
  hlcm_aoi_i <- tryCatch(crop(hlcm_int, ext(vect(aoi_buf))), error = function(e) NULL)
  wood_polys <- if (!is.null(hlcm_aoi_i)) {
    hlcm_w <- ifel(hlcm_aoi_i %in% woodland_lc_values, 1, NA)
    tryCatch(as.polygons(hlcm_w) %>% st_as_sf() %>% st_transform(27700) %>% st_make_valid(),
             error = function(e) st_sf(geometry = st_sfc(crs = 27700)))
  } else st_sf(geometry = st_sfc(crs = 27700))

  # Combine habitat sf objects
  hab_known_list <- list()
  if (!is.null(sg_k)) hab_known_list[["Seagrass"]]  <- st_sf(habitat = "Seagrass", geometry = st_geometry(sg_k))
  if (!is.null(sm_k)) hab_known_list[["Saltmarsh"]]  <- st_sf(habitat = "Saltmarsh", geometry = st_geometry(sm_k))
  if (!is.null(oy_k)) hab_known_list[["Oyster"]]     <- st_sf(habitat = "Oyster", geometry = st_geometry(oy_k))
  if (!is.null(dn_k)) hab_known_list[["Sand dune"]]  <- st_sf(habitat = "Sand dune", geometry = st_geometry(dn_k))
  if (nrow(wood_polys) > 0) hab_known_list[["Woodland"]] <- st_sf(habitat = "Woodland", geometry = st_geometry(wood_polys))
  hab_known_sf <- if (length(hab_known_list) > 0) bind_rows(hab_known_list) else NULL

  hab_pot_list <- list()
  if (!is.null(sg_p)) hab_pot_list[["Seagrass"]] <- st_sf(habitat = "Seagrass", geometry = st_geometry(sg_p))
  if (!is.null(oy_p)) hab_pot_list[["Oyster"]]   <- st_sf(habitat = "Oyster", geometry = st_geometry(oy_p))
  if (!is.null(dn_p)) hab_pot_list[["Sand dune"]] <- st_sf(habitat = "Sand dune", geometry = st_geometry(dn_p))
  hab_pot_sf <- if (length(hab_pot_list) > 0) bind_rows(hab_pot_list) else NULL

  # Shared panel settings
  aoi_coord <- coord_sf(xlim = c(aoi_bbox["xmin"], aoi_bbox["xmax"]),
                        ylim = c(aoi_bbox["ymin"], aoi_bbox["ymax"]), crs = 27700)
  aoi_base <- geom_sf(data = dg, fill = "grey95", colour = "grey60", linewidth = 0.2)
  panel_theme <- theme_scamp() + theme(legend.position = "right",
                                        legend.key.size = unit(0.3, "cm"),
                                        legend.justification = "top",
                                        legend.box.just = "top",
                                        plot.title = element_text(size = 10, face = "bold"),
                                        plot.margin = margin(2, 2, 2, 2))

  # ---- Panel A: Threats ----
  p_threats <- ggplot() + aoi_base

  if (!is.null(cum_aoi))
    p_threats <- p_threats + tidyterra::geom_spatraster(data = cum_aoi, alpha = 0.6) +
      scale_fill_viridis_c(name = "Cumulative\npressure", na.value = "transparent", option = "magma")

  if (!is.null(hybas_aoi) && "mean_diffuse_risk" %in% names(hybas_aoi))
    p_threats <- p_threats + ggnewscale::new_scale_fill() +
      geom_sf(data = hybas_aoi, aes(fill = mean_diffuse_risk),
              colour = "grey50", linewidth = 0.15, alpha = 0.35) +
      scale_fill_gradient(name = "Catchment\nmean risk", low = "#f5deb3", high = "#a50026",
                          na.value = "transparent",
                          guide = guide_colourbar(title.position = "top"),
                          labels = function(x) ifelse(x == min(x), "Low", ifelse(x == max(x), "High", "")))

  if (!is.null(riv_plot) && nrow(riv_plot) > 0)
    p_threats <- p_threats + ggnewscale::new_scale_colour() +
      geom_sf(data = riv_plot, aes(colour = accum_risk_raw), linewidth = 1.2) +
      scale_colour_gradient(name = "River\naccum. risk", low = "steelblue", high = "orangered",
                            trans = "log1p",
                            breaks = range(riv_plot$accum_risk_raw, na.rm = TRUE),
                            labels = c("Low", "High"))

  if (!is.null(vt_aoi_i) && nrow(vt_aoi_i) > 0) {
    vt_aoi_i$threat_label <- gsub("_", " ", tools::toTitleCase(vt_aoi_i$threat_type))
    vt_pts_l <- vt_aoi_i[st_geometry_type(vt_aoi_i) %in% c("POINT", "MULTIPOINT"), ]
    vt_pgs_l <- vt_aoi_i[st_geometry_type(vt_aoi_i) %in% c("POLYGON", "MULTIPOLYGON", "GEOMETRYCOLLECTION"), ]
    if (nrow(vt_pgs_l) > 0)
      p_threats <- p_threats + geom_sf(data = vt_pgs_l, aes(shape = threat_label),
                                       fill = NA, colour = "black", linewidth = 0.8)
    if (nrow(vt_pts_l) > 0)
      p_threats <- p_threats + geom_sf(data = vt_pts_l, aes(shape = threat_label),
                                       colour = "black", size = 3, stroke = 1.2) +
        scale_shape_manual(name = "Vector threats",
                           values = c("Dredge spoil" = 0, "Dredging" = 1,
                                      "Dumped munitions" = 2, "Ports" = 5,
                                      "Waste discharge" = 6, "Ww wrecks" = 4),
                           guide = guide_legend(override.aes = list(size = 3.5),
                                                keywidth = unit(0.6, "cm"),
                                                keyheight = unit(0.5, "cm")))
  }

  p_threats <- p_threats +
    geom_sf(data = smr, fill = NA, colour = "black", linewidth = 0.5) +
    geom_sf(data = aoi_i, fill = NA, colour = "black", linewidth = 0.6, linetype = "dashed") +
    aoi_coord + labs(title = paste0(aoi_label, " — A. Threats & Pressures")) + panel_theme

  # ---- Panel B: Habitats ----
  p_habs <- ggplot() + aoi_base

  if (!is.null(hab_known_sf))
    p_habs <- p_habs + geom_sf(data = hab_known_sf, aes(fill = habitat),
                                 colour = "grey20", linewidth = 0.4, alpha = 0.6) +
      scale_fill_manual(name = "Known habitat", values = hab_fill_pal)

  if (!is.null(hab_pot_sf))
    p_habs <- p_habs + ggnewscale::new_scale_fill() +
      geom_sf(data = hab_pot_sf, aes(fill = habitat),
              colour = "black", linewidth = 0.6, alpha = 0.2, linetype = "longdash") +
      scale_fill_manual(name = "Potential habitat", values = hab_fill_pal)

  if (!is.null(riv_plot) && nrow(riv_plot) > 0)
    p_habs <- p_habs + geom_sf(data = riv_plot, colour = "steelblue", linewidth = 0.6)

  if (!is.null(oy_hist_aoi) && nrow(oy_hist_aoi) > 0)
    p_habs <- p_habs + geom_sf(data = oy_hist_aoi, colour = "#de2d26",
                                 shape = 17, size = 2.5, alpha = 0.7) +
      annotate("text", x = mean(st_coordinates(oy_hist_aoi)[, 1]),
               y = min(st_coordinates(oy_hist_aoi)[, 2]) - 800,
               label = "Historical\noyster beds", colour = "#de2d26", size = 2.5,
               fontface = "italic")

  p_habs <- p_habs +
    geom_sf(data = smr, fill = NA, colour = "black", linewidth = 0.5) +
    geom_sf(data = aoi_i, fill = NA, colour = "black", linewidth = 0.6, linetype = "dashed") +
    aoi_coord + labs(title = paste0(aoi_label, " — B. Known & Potential Habitats")) + panel_theme

  # ---- Panel C: Species Richness ----
  rich_aoi <- tryCatch(
    crop(rich_surface, ext(vect(st_transform(aoi_buf, crs(rich_surface))))) %>%
      project("EPSG:27700", method = "bilinear"),
    error = function(e) NULL)

  p_rich <- ggplot() + aoi_base
  if (!is.null(rich_aoi))
    p_rich <- p_rich + tidyterra::geom_spatraster(data = rich_aoi, alpha = 0.75) +
      scale_fill_viridis_c(name = "Rarefied\nrichness", na.value = "transparent", option = "mako",
                           direction = -1)
  if (!is.null(hybas_aoi))
    p_rich <- p_rich + geom_sf(data = hybas_aoi, fill = NA, colour = "grey60",
                                linewidth = 0.15, linetype = "dotted")
  if (!is.null(riv_plot) && nrow(riv_plot) > 0)
    p_rich <- p_rich + geom_sf(data = riv_plot, colour = "steelblue", linewidth = 0.4, alpha = 0.5)
  p_rich <- p_rich +
    geom_sf(data = smr, fill = NA, colour = "black", linewidth = 0.5) +
    geom_sf(data = aoi_i, fill = NA, colour = "black", linewidth = 0.6, linetype = "dashed") +
    aoi_coord + labs(title = paste0(aoi_label, " — C. Species Richness")) + panel_theme

  # ---- Panel D: Beta Diversity ----
  beta_aoi <- tryCatch(
    crop(beta_surface, ext(vect(st_transform(aoi_buf, crs(beta_surface))))) %>%
      project("EPSG:27700", method = "bilinear"),
    error = function(e) NULL)

  p_beta <- ggplot() + aoi_base
  if (!is.null(beta_aoi))
    p_beta <- p_beta + tidyterra::geom_spatraster(data = beta_aoi, alpha = 0.75) +
      scale_fill_viridis_c(name = "Beta\ndiversity", na.value = "transparent", option = "inferno")
  if (!is.null(hybas_aoi))
    p_beta <- p_beta + geom_sf(data = hybas_aoi, fill = NA, colour = "grey60",
                                 linewidth = 0.15, linetype = "dotted")
  if (!is.null(riv_plot) && nrow(riv_plot) > 0)
    p_beta <- p_beta + geom_sf(data = riv_plot, colour = "steelblue", linewidth = 0.4, alpha = 0.5)
  p_beta <- p_beta +
    geom_sf(data = smr, fill = NA, colour = "black", linewidth = 0.5) +
    geom_sf(data = aoi_i, fill = NA, colour = "black", linewidth = 0.6, linetype = "dashed") +
    aoi_coord + labs(title = paste0(aoi_label, " — D. Beta Diversity")) + panel_theme

  # ---- Save individual panels ----
  ggsave(file.path(aoi_dir, "a_threats.png"), p_threats, width = 8, height = 7, dpi = 300, bg = "white")
  ggsave(file.path(aoi_dir, "b_habitats.png"), p_habs, width = 8, height = 7, dpi = 300, bg = "white")
  ggsave(file.path(aoi_dir, "c_richness.png"), p_rich, width = 8, height = 7, dpi = 300, bg = "white")
  ggsave(file.path(aoi_dir, "d_beta_diversity.png"), p_beta, width = 8, height = 7, dpi = 300, bg = "white")
  cat(sprintf("  Saved 4 panels to %s/\n", aoi_dir))

}
cat("\n==== Batch export complete ====\n")

# ============================================================
# OVERVIEW COMPANION MAPS: RICHNESS & BETA DIVERSITY
# ============================================================

cat("\n==== Overview companion maps ====\n")

if (!exists("rich_surface"))
  rich_surface <- rast("outputs/species/rarefied_richness_surface_marine.tif")
if (!exists("beta_surface"))
  beta_surface <- rast("outputs/species/beta_diversity_surface_marine.tif")

rich_bng <- project(rich_surface, "EPSG:27700", method = "bilinear")
beta_bng <- project(beta_surface, "EPSG:27700", method = "bilinear")

p_richness_overview <- base_map("Solway Overview — Rarefied Species Richness",
                                subtitle = "IDW-interpolated marine survey data",
                                tile_type = "cartolight", zoom = 10) +
  tidyterra::geom_spatraster(data = rich_bng, alpha = 0.75) +
  scale_fill_viridis_c(name = "Rarefied\nrichness", na.value = "transparent",
                       option = "mako", direction = -1) +
  geom_sf(data = aois, fill = NA, colour = "white", linewidth = 1.2) +
  geom_sf(data = aois, fill = NA, colour = "black", linewidth = 0.6) +
  theme(legend.position = "right", legend.box = "vertical",
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14),
        axis.text = element_text(size = 12))
print(p_richness_overview)

ggsave(file.path(ppt_dir, "overview_richness.png"), p_richness_overview,
       width = 14, height = 11, dpi = 300, bg = "white")
cat("  Saved: figures/ppt_panels/overview_richness.png\n")

p_beta_overview <- base_map("Solway Overview — Beta Diversity",
                            subtitle = "IDW-interpolated marine survey data",
                            tile_type = "cartolight", zoom = 10) +
  tidyterra::geom_spatraster(data = beta_bng, alpha = 0.75) +
  scale_fill_viridis_c(name = "Beta\ndiversity", na.value = "transparent",
                       option = "inferno") +
  geom_sf(data = aois, fill = NA, colour = "white", linewidth = 1.2) +
  geom_sf(data = aois, fill = NA, colour = "black", linewidth = 0.6) +
  theme(legend.position = "right", legend.box = "vertical",
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14),
        axis.text = element_text(size = 12))
print(p_beta_overview)

ggsave(file.path(ppt_dir, "overview_beta_diversity.png"), p_beta_overview,
       width = 14, height = 11, dpi = 300, bg = "white")
cat("  Saved: figures/ppt_panels/overview_beta_diversity.png\n")

# ============================================================
# HEXAGON-LEVEL COASTAL CONNECTIVITY
# ============================================================
# Plots the actual coastalNet hexagons and their pairwise dispersal
# connections directly — no bay-level aggregation. AOI polygons are
# used only to identify which hexagons sit within which named bay
# (for a label/outline), not to collapse the connectivity data.
# The previous bay-centroid flow-line version produced misleading
# results (e.g. Loch Ryan showing strong "connectivity" to Wigtown
# but not Luce) because collapsing many hex-to-hex pairs into a
# single mean/max per bay pair discards the real spatial pattern.
# ============================================================

cat("\n==== Hexagon-level coastal connectivity ====\n")

if (!requireNamespace("coastalNet", quietly = TRUE)) {
  remotes::install_github("jorgeassis/coastalNet")
}
library(coastalNet)

# ---- Load saved connectivity data (pairwise matrix + pairs table) ----
sg_conn <- readRDS("outputs/physical/coastalnet_seagrass_connectivity.rds")
oy_conn <- readRDS("outputs/physical/coastalnet_oyster_connectivity.rds")

# ---- Load hexagon grid (native coastalNet CRS, do not transform yet —
# mapConnectivity() builds lines from these geometries directly) ----
sg_hex <- st_read("outputs/physical/seagrass_connectivity_hexagons.gpkg", quiet = TRUE)
oy_hex <- st_read("outputs/physical/oyster_connectivity_hexagons.gpkg", quiet = TRUE)

# ---- Identify which hexagons fall in which named bay (context only) ----
aois_hexcrs <- st_transform(aois[, "bay_name"], st_crs(sg_hex))
sg_hex <- st_join(sg_hex, aois_hexcrs, join = st_intersects, largest = TRUE)
oy_hex <- st_join(oy_hex, aois_hexcrs, join = st_intersects, largest = TRUE)

cat(sprintf("  Seagrass hexagons within a named bay: %d of %d\n",
            sum(!is.na(sg_hex$bay_name)), nrow(sg_hex)))
cat(sprintf("  Oyster hexagons within a named bay: %d of %d\n",
            sum(!is.na(oy_hex$bay_name)), nrow(oy_hex)))
print(table(sg_hex$bay_name, useNA = "ifany"))

# ---- Build hexagon-to-hexagon connection lines (coastalNet's own logic) ----
mappedSeagrass <- mapConnectivity(obj = sg_hex, featureName = "ID",
                                  connectivityPairs = sg_conn$connectivityPairs)
mappedOyster <- mapConnectivity(obj = oy_hex, featureName = "ID",
                                connectivityPairs = oy_conn$connectivityPairs)

# ---- Reproject everything to BNG for plotting on the same base map ----
sg_hex_bng   <- st_transform(sg_hex, 27700)
oy_hex_bng   <- st_transform(oy_hex, 27700)
sg_lines_bng <- st_transform(mappedSeagrass$lineConnections, 27700)
oy_lines_bng <- st_transform(mappedOyster$lineConnections, 27700)

bay_labels <- st_centroid(aois[, "bay_name"])

hex_status_pal <- c("Existing" = "#1b7837", "Potential" = "#ff7f00", "Neither" = NA)

# ---- Curve straight hex-to-hex lines into arcs ----
# Many pairs share nearly the same path between hexagon centroids, so
# straight lines stack on top of each other and become indistinguishable.
# Bowing each line into a quadratic-bezier arc (perpendicular offset
# proportional to line length) visually separates overlapping paths.
curve_sf_lines <- function(lines_sf, curvature = 0.3, n_pts = 30) {
  arcs <- lapply(seq_len(nrow(lines_sf)), function(i) {
    coords <- st_coordinates(lines_sf[i, ])[, 1:2, drop = FALSE]
    p1 <- coords[1, ]
    p2 <- coords[nrow(coords), ]
    mid <- (p1 + p2) / 2
    d <- p2 - p1
    len <- sqrt(sum(d^2))
    if (len == 0) return(st_linestring(rbind(p1, p2)))
    perp <- c(-d[2], d[1]) / len
    ctrl <- mid + perp * curvature * len
    t <- seq(0, 1, length.out = n_pts)
    x <- (1 - t)^2 * p1[1] + 2 * (1 - t) * t * ctrl[1] + t^2 * p2[1]
    y <- (1 - t)^2 * p1[2] + 2 * (1 - t) * t * ctrl[2] + t^2 * p2[2]
    st_linestring(cbind(x, y))
  })
  st_sf(st_drop_geometry(lines_sf), geometry = st_sfc(arcs, crs = st_crs(lines_sf)))
}

sg_lines_bng <- curve_sf_lines(sg_lines_bng)
oy_lines_bng <- curve_sf_lines(oy_lines_bng)

# Draw weakest connections first, strongest last (on top) so strong
# connections aren't buried under a mass of near-zero clutter.
sg_lines_bng <- sg_lines_bng[order(sg_lines_bng$Value), ]
oy_lines_bng <- oy_lines_bng[order(oy_lines_bng$Value), ]

# ---- Plot: Seagrass hexagon connectivity ----
p_sg_hex <- base_map("Hexagon Connectivity — Seagrass",
                     subtitle = "coastalNet forward dispersal probability (180-day window)",
                     tile_type = "cartolight", zoom = 10) +
  geom_sf(data = aois, fill = NA, colour = "steelblue", linewidth = 0.4, linetype = "dashed") +
  geom_sf(data = sg_hex_bng, aes(fill = status), colour = "grey40", linewidth = 0.15, alpha = 0.35) +
  scale_fill_manual(name = "Hexagon status", values = hex_status_pal, na.value = NA) +
  ggnewscale::new_scale_colour() +
  geom_sf(data = sg_lines_bng, aes(colour = Value, linewidth = Value), alpha = 0.85,
          arrow = grid::arrow(length = unit(0.12, "cm"), type = "closed", angle = 20)) +
  scale_colour_viridis_c(name = "Connectivity\nprobability", option = "magma", trans = "log1p") +
  scale_linewidth_continuous(range = c(0.4, 2.5), guide = "none") +
  theme(legend.position = "right", legend.box = "vertical",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        axis.text = element_text(size = 12))
print(p_sg_hex)

ggsave(file.path(ppt_dir, "hexagon_connectivity_seagrass.png"), p_sg_hex,
       width = 14, height = 11, dpi = 300, bg = "white")
cat("  Saved: figures/ppt_panels/hexagon_connectivity_seagrass.png\n")

# ---- Plot: Oyster hexagon connectivity ----
p_oy_hex <- base_map("Hexagon Connectivity — Native Oyster",
                     subtitle = "coastalNet forward dispersal probability (180-day window)",
                     tile_type = "cartolight", zoom = 10) +
  geom_sf(data = aois, fill = NA, colour = "steelblue", linewidth = 0.4, linetype = "dashed") +
  geom_sf(data = oy_hex_bng, aes(fill = status), colour = "grey40", linewidth = 0.15, alpha = 0.35) +
  scale_fill_manual(name = "Hexagon status", values = hex_status_pal, na.value = NA) +
  ggnewscale::new_scale_colour() +
  geom_sf(data = oy_lines_bng, aes(colour = Value, linewidth = Value), alpha = 0.85,
          arrow = grid::arrow(length = unit(0.12, "cm"), type = "closed", angle = 20)) +
  scale_colour_viridis_c(name = "Connectivity\nprobability", option = "magma", trans = "log1p") +
  scale_linewidth_continuous(range = c(0.4, 2.5), guide = "none") +
  theme(legend.position = "right", legend.box = "vertical",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        axis.text = element_text(size = 12))
print(p_oy_hex)

ggsave(file.path(ppt_dir, "hexagon_connectivity_oyster.png"), p_oy_hex,
       width = 14, height = 11, dpi = 300, bg = "white")
cat("  Saved: figures/ppt_panels/hexagon_connectivity_oyster.png\n")

cat("\n==== All overview and connectivity maps complete ====\n")


