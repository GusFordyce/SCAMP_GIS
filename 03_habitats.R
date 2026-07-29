##=====================================================================
## SCAMP GIS — Theme D: Habitats
## 03_habitats.R
##=====================================================================
## One habitat at a time. This section: SEAGRASS.
## Independent script — loads its own packages and boundaries.
##=====================================================================
install.packages(c('maptiles', 'tidyterra'))
require(sf)
require(dplyr)
require(ggplot2)
require(ggspatial)
require(cowplot)
require(maptiles)
require(tidyterra)

##=====================================================================
## 0. SETUP — Boundaries and AOI
##=====================================================================

smr <- st_read(
  "data/raw/boundaries/administrative_units_scottish_marine_regions/administrative_units_scottish_marine_regionsPolygon.shp",
  quiet = TRUE
)
las <- st_read(
  "data/raw/boundaries/Local_Authority_Boundaries_-_Scotland/pub_las.shp",
  quiet = TRUE
)

smr_solway <- smr %>% filter(grepl("Solway", objnam, ignore.case = TRUE))
las_dg     <- las %>% filter(grepl("dumfries", local_auth, ignore.case = TRUE))

smr_solway_bng <- st_transform(smr_solway, 27700)
las_dg_bng     <- st_transform(las_dg, 27700)

smr_buf       <- st_buffer(smr_solway_bng, dist = 8047)
coastal_strip <- st_intersection(st_geometry(smr_buf), st_geometry(las_dg_bng))
scamp_aoi     <- st_union(c(st_geometry(smr_solway_bng), coastal_strip))
scamp_aoi_sf  <- st_sf(geometry = scamp_aoi, crs = 27700)

scamp_coastline <- st_intersection(
  st_cast(st_boundary(las_dg_bng), "MULTILINESTRING"),
  scamp_aoi_sf
)

dir.create("outputs/habitats",  showWarnings = FALSE, recursive = TRUE)
dir.create("figures/habitats",  showWarnings = FALSE, recursive = TRUE)

map_theme <- theme_minimal(base_size = 10) +
  theme(
    legend.position  = "right",
    plot.title       = element_text(size = 12, face = "bold"),
    plot.subtitle    = element_text(size = 9, colour = "grey30"),
    panel.grid.major = element_line(colour = "grey85", linewidth = 0.2),
    axis.text        = element_text(size = 7)
  )


##=====================================================================
## 1. SEAGRASS
##=====================================================================

cat("\n=== 1. SEAGRASS ===\n")

## ---- 1a. GeMS seagrass — habitat polygons & points ----

gems_path <- "data/raw/GEMS-PMF_SCOTLAND_GPKG_4326/gems.gpkg"

gems_hab_py <- st_read(gems_path, layer = "gems_habitat_polygon_dataset",
                        quiet = TRUE) %>%
  filter(PMF == "Seagrass beds") %>%
  st_transform(27700) %>%
  st_make_valid() %>%
  st_intersection(scamp_aoi_sf) %>%
  mutate(source = "GeMS habitat polygon")

gems_hab_pt <- st_read(gems_path, layer = "gems_habitat_point_dataset",
                        quiet = TRUE) %>%
  filter(PMF == "Seagrass beds") %>%
  st_transform(27700) %>%
  st_make_valid() %>%
  st_intersection(scamp_aoi_sf) %>%
  mutate(source = "GeMS habitat point")

cat(sprintf("  GeMS seagrass polygons in AOI: %d\n", nrow(gems_hab_py)))
cat(sprintf("  GeMS seagrass points in AOI:   %d\n", nrow(gems_hab_pt)))

## Combine into one GeMS layer: cast points to geometry collection won't work,
## so keep polygons and points as separate sf objects tagged by source
gems_sg_py <- gems_hab_py %>% select(source, Shape)
st_geometry(gems_sg_py) <- "geometry"
gems_sg_pt <- gems_hab_pt %>% select(source, Shape)
st_geometry(gems_sg_pt) <- "geometry"

## Cluster polygons spatially: buffer by 2 km, union, decompose into clusters
gems_centroids <- st_centroid(gems_sg_py)
gems_buf <- st_buffer(gems_centroids, dist = 2000) %>% st_union() %>%
  st_cast("POLYGON")
cluster_ids <- st_intersects(gems_centroids, gems_buf)
gems_sg_py$cluster <- sapply(cluster_ids, `[`, 1)

cluster_stats <- gems_sg_py %>%
  mutate(area_m2 = as.numeric(st_area(geometry))) %>%
  st_drop_geometry() %>%
  group_by(cluster) %>%
  summarise(n = n(), total_ha = sum(area_m2) / 10000, .groups = "drop")

cluster_centres <- gems_sg_py %>%
  group_by(cluster) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  st_centroid() %>%
  bind_cols(st_coordinates(.) %>% as.data.frame()) %>%
  left_join(cluster_stats, by = "cluster") %>%
  mutate(label = sprintf("n=%d\n%.2f ha", n, total_ha))

## Offset label positions outward from cluster centres for readability
aoi_centre <- st_coordinates(st_centroid(scamp_aoi_sf))
cluster_centres <- cluster_centres %>%
  mutate(
    dx = X - aoi_centre[1, "X"],
    dy = Y - aoi_centre[1, "Y"],
    len = sqrt(dx^2 + dy^2),
    nudge_x = (dx / len) * 15000,
    nudge_y = (dy / len) * 15000,
    label_x = X + nudge_x,
    label_y = Y + nudge_y
  )

p_gems_sg <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
          linewidth = 0.3) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  {if (nrow(gems_sg_py) > 0)
    geom_sf(data = gems_sg_py, fill = "#1b9e77", colour = "#1b9e77", alpha = 0.5)} +
  geom_segment(data = cluster_centres,
               aes(x = label_x, y = label_y, xend = X, yend = Y),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
               colour = "grey20", linewidth = 0.4) +
  geom_label(data = cluster_centres,
             aes(x = label_x, y = label_y, label = label),
             size = 2.8, fill = "white", colour = "grey20",
             label.padding = unit(0.2, "lines"), fontface = "bold") +
  labs(title = "GeMS Seagrass Beds (SCAMP AOI)",
       subtitle = sprintf("Polygons: %d across %d clusters | Points: %d",
                          nrow(gems_sg_py), nrow(cluster_stats), nrow(gems_sg_pt))) +
  map_theme + coord_sf(crs = 27700)

ggsave("figures/habitats/seagrass_gems.png", p_gems_sg,
       width = 14, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/habitats/seagrass_gems.png\n")

## Diagnostics: geometry types and area range
cat("\n  -- GeMS polygon diagnostics --\n")
cat("  Geometry types:\n")
print(table(st_geometry_type(gems_sg_py)))
gems_areas <- as.numeric(st_area(gems_sg_py))
cat(sprintf("  Area range: %.2f — %.2f m²  (median %.2f m²)\n",
    min(gems_areas), max(gems_areas), median(gems_areas)))
cat(sprintf("  Polygons with zero area: %d\n", sum(gems_areas == 0)))
cat(sprintf("  Polygons > 1 ha: %d\n", sum(gems_areas > 10000)))

p_gems_areas <- ggplot(data.frame(area_m2 = gems_areas), aes(x = area_m2)) +
  geom_histogram(bins = 50, fill = "#1b9e77", colour = "white") +
  scale_x_log10(labels = scales::comma) +
  labs(title = "GeMS Seagrass polygon areas",
       x = "Area (m², log scale)", y = "Count") +
  map_theme

ggsave("figures/habitats/seagrass_gems_area_hist.png", p_gems_areas,
       width = 8, height = 5, units = "in", dpi = 200)
cat("  Saved: figures/habitats/seagrass_gems_area_hist.png\n")

## Zoomed inset: Loch Ryan (approx BNG 306000, 560000 centre)
p_gems_lochRyan <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
          linewidth = 0.3) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  {if (nrow(gems_sg_py) > 0)
    geom_sf(data = gems_sg_py, fill = "#1b9e77", colour = "#1b9e77", alpha = 0.6)} +
  #{if (nrow(gems_sg_pt) > 0)
    #geom_sf(data = gems_sg_pt, colour = "#d95f02", size = 2.5, alpha = 0.8)} +
  coord_sf(crs = 27700,
           xlim = c(200000, 215000), ylim = c(555000, 575000)) +
  labs(title = "GeMS Seagrass — Loch Ryan inset",
       subtitle = sprintf("Green = polygons (%d) | Orange = points (%d)",
                          nrow(gems_sg_py), nrow(gems_sg_pt))) +
  map_theme

ggsave("figures/habitats/seagrass_gems_lochRyan.png", p_gems_lochRyan,
       width = 10, height = 10, units = "in", dpi = 200)
cat("  Saved: figures/habitats/seagrass_gems_lochRyan.png\n")


## ---- 1b. WCMC seagrass polygons ----

wcmc_sg_py <- st_read(
  "data/raw/habitats/WCMC013-014_SeagrassPtPy2021_v7_1/014_001_WCMC013-014_SeagrassPtPy2021_v7_1/01_Data/WCMC013_014_Seagrasses_Py_v7_1.shp",
  quiet = TRUE
) %>%
  st_transform(27700) %>%
  st_make_valid() %>%
  st_intersection(scamp_aoi_sf) %>%
  mutate(source = "WCMC polygon") %>%
  select(source, geometry)

wcmc_sg_pt <- st_read(
  "data/raw/habitats/WCMC013-014_SeagrassPtPy2021_v7_1/014_001_WCMC013-014_SeagrassPtPy2021_v7_1/01_Data/WCMC013_014_Seagrasses_Pt_v7_1.shp",
  quiet = TRUE
) %>%
  st_transform(27700) %>%
  st_make_valid() %>%
  st_intersection(scamp_aoi_sf) %>%
  mutate(source = "WCMC point") %>%
  select(source, geometry)

cat(sprintf("  WCMC seagrass polygons in AOI: %d\n", nrow(wcmc_sg_py)))
cat(sprintf("  WCMC seagrass points in AOI:   %d\n", nrow(wcmc_sg_pt)))

cat("\n  -- WCMC polygon diagnostics --\n")
cat("  Geometry types:\n")
print(table(st_geometry_type(wcmc_sg_py)))
wcmc_areas <- as.numeric(st_area(wcmc_sg_py))
cat(sprintf("  Area range: %.2f — %.2f m²  (median %.2f m²)\n",
    min(wcmc_areas), max(wcmc_areas), median(wcmc_areas)))

p_wcmc_py_only <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
          linewidth = 0.3) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  geom_sf(data = wcmc_sg_py, fill = "#7570b3", colour = "#7570b3", alpha = 0.6) +
  labs(title = "WCMC Seagrass polygons only (SCAMP AOI)",
       subtitle = sprintf("Polygons: %d", nrow(wcmc_sg_py))) +
  map_theme + coord_sf(crs = 27700)

ggsave("figures/habitats/seagrass_wcmc_only.png", p_wcmc_py_only,
       width = 12, height = 8, units = "in", dpi = 200)
cat("  Saved: figures/habitats/seagrass_wcmc_only.png\n")

## WCMC cluster statistics (same method as GeMS)
wcmc_centroids <- st_centroid(wcmc_sg_py)
wcmc_buf <- st_buffer(wcmc_centroids, dist = 2000) %>% st_union() %>%
  st_cast("POLYGON")
wcmc_cluster_ids <- st_intersects(wcmc_centroids, wcmc_buf)
wcmc_sg_py$cluster <- sapply(wcmc_cluster_ids, `[`, 1)

wcmc_cluster_stats <- wcmc_sg_py %>%
  mutate(area_m2 = as.numeric(st_area(geometry))) %>%
  st_drop_geometry() %>%
  group_by(cluster) %>%
  summarise(n = n(), total_ha = sum(area_m2) / 10000, .groups = "drop")

wcmc_cluster_centres <- wcmc_sg_py %>%
  group_by(cluster) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  st_centroid() %>%
  bind_cols(st_coordinates(.) %>% as.data.frame()) %>%
  left_join(wcmc_cluster_stats, by = "cluster") %>%
  mutate(label = sprintf("n=%d\n%.2f ha", n, total_ha))

cat("\n  -- WCMC cluster statistics --\n")
print(wcmc_cluster_stats)
cat(sprintf("  WCMC clusters: %d\n", nrow(wcmc_cluster_stats)))

## Zoom panels: Loch Ryan + each WCMC cluster location
zoom_panels <- data.frame(
  name = c("Loch Ryan", "Kirkcudbright Bay",
           "Auchencairn Bay & Rough Island"),
  xmin = c(200000, 248000, 278000),
  xmax = c(215000, 262000, 295000),
  ymin = c(555000, 546000, 548000),
  ymax = c(575000, 560000, 565000),
  stringsAsFactors = FALSE
)

cat(sprintf("  Producing %d zoom panels with insets...\n", nrow(zoom_panels)))

## Build bounding boxes for the inset map
zoom_boxes <- do.call(rbind, lapply(seq_len(nrow(zoom_panels)), function(i) {
  zp <- zoom_panels[i, ]
  st_sf(
    name = zp$name,
    geometry = st_sfc(st_polygon(list(matrix(c(
      zp$xmin, zp$ymin,
      zp$xmax, zp$ymin,
      zp$xmax, zp$ymax,
      zp$xmin, zp$ymax,
      zp$xmin, zp$ymin
    ), ncol = 2, byrow = TRUE))), crs = 27700)
  )
}))

## Simple check: plot all bboxes on the AOI
p_bboxes <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
          linewidth = 0.3) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  geom_sf(data = zoom_boxes, fill = NA, colour = "red", linewidth = 0.8) +
  geom_sf_label(data = zoom_boxes, aes(label = name), colour = "red",
                size = 3, fontface = "bold") +
  labs(title = "Seagrass cluster bounding boxes",
       subtitle = sprintf("%d clusters defined", nrow(zoom_boxes))) +
  map_theme + coord_sf(crs = 27700)

ggsave("figures/habitats/seagrass_cluster_bboxes.png", p_bboxes,
       width = 14, height = 10, units = "in", dpi = 200)
cat("  Saved: figures/habitats/seagrass_cluster_bboxes.png\n")

## Inset map (reused in each panel)
p_inset_base <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.15) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
          linewidth = 0.15) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.2) +
  geom_sf(data = zoom_boxes, fill = NA, colour = "red", linewidth = 0.6) +
  coord_sf(crs = 27700) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.4)
  )

zoom_final <- list()
for (i in seq_len(nrow(zoom_panels))) {
  zp <- zoom_panels[i, ]

  ## Highlight current box in the inset
  current_box <- zoom_boxes[i, ]
  p_inset_i <- p_inset_base +
    geom_sf(data = current_box, fill = "red", colour = "red",
            alpha = 0.25, linewidth = 0.8)

  ## Fetch satellite tiles for this panel's extent
  tile_bbox <- st_sf(geometry = st_sfc(st_polygon(list(matrix(c(
    zp$xmin, zp$ymin,
    zp$xmax, zp$ymin,
    zp$xmax, zp$ymax,
    zp$xmin, zp$ymax,
    zp$xmin, zp$ymin
  ), ncol = 2, byrow = TRUE))), crs = 27700))
  sat_tiles <- get_tiles(tile_bbox, provider = "Esri.WorldImagery",
                         zoom = 14, crop = TRUE)

  ## Main zoom panel with satellite basemap
  p_main <- ggplot() +
    geom_spatraster_rgb(data = sat_tiles) +
    geom_sf(data = scamp_aoi_sf, fill = NA, colour = "white", linewidth = 0.5) +
    geom_sf(data = gems_sg_py, fill = "#1b9e77", colour = "white", alpha = 0.6,
            linewidth = 0.3) +
    geom_sf(data = wcmc_sg_py, fill = "#7570b3", colour = "white", alpha = 0.6,
            linewidth = 0.3) +
    geom_sf(data = gems_sg_pt, colour = "#d95f02", size = 2.5, alpha = 0.9, shape = 16) +
    geom_sf(data = wcmc_sg_pt, colour = "#a6761d", size = 3, alpha = 0.9, shape = 17) +
    coord_sf(crs = 27700,
             xlim = c(zp$xmin, zp$xmax), ylim = c(zp$ymin, zp$ymax)) +
    labs(title = zp$name,
         subtitle = "Green = GeMS | Purple = WCMC | Orange/Brown = points") +
    map_theme + theme(legend.position = "none")

  ## Compose with cowplot: main panel + small inset in top-right corner
  p_combined <- ggdraw(p_main) +
    draw_plot(p_inset_i, x = 0.68, y = 0.65, width = 0.3, height = 0.3)

  zoom_final[[i]] <- p_combined
}

p_zoom_grid <- plot_grid(plotlist = zoom_final, ncol = 3)

ggsave("figures/habitats/seagrass_wcmc_gems_zoom_panels.png", p_zoom_grid,
       width = 18, height = 7,
       units = "in", dpi = 300)
cat("  Saved: figures/habitats/seagrass_wcmc_gems_zoom_panels.png\n")


p_wcmc_gems <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
          linewidth = 0.3) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  ## Polygons — each source as own layer, mapping fill for legend
  geom_sf(data = gems_sg_py, aes(fill = "GeMS polygon"), colour = NA, alpha = 0.5) +
  geom_sf(data = wcmc_sg_py, aes(fill = "WCMC polygon"), colour = NA, alpha = 0.5) +
  ## GeMS cluster arrows
  geom_segment(data = cluster_centres,
               aes(x = label_x, y = label_y, xend = X, yend = Y),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
               colour = "grey20", linewidth = 0.4) +
  geom_label(data = cluster_centres,
             aes(x = label_x, y = label_y, label = label),
             size = 2.5, fill = "white", colour = "#1b9e77",
             label.padding = unit(0.2, "lines"), fontface = "bold") +
  ## Points — each source as own layer
  geom_sf(data = gems_sg_pt, aes(colour = "GeMS point"), shape = 16,
          size = 2, alpha = 0.8) +
  geom_sf(data = wcmc_sg_pt, aes(colour = "WCMC point"), shape = 17,
          size = 3, alpha = 0.8) +
  scale_fill_manual(
    values = c("GeMS polygon" = "#1b9e77", "WCMC polygon" = "#7570b3"),
    name = "Polygons") +
  scale_colour_manual(
    values = c("GeMS point" = "#d95f02", "WCMC point" = "#a6761d"),
    name = "Points") +
  labs(title = "Seagrass — GeMS + WCMC (SCAMP AOI)",
       subtitle = sprintf(
         "GeMS poly: %d | GeMS pts: %d | WCMC poly: %d | WCMC pts: %d",
         nrow(gems_sg_py), nrow(gems_sg_pt), nrow(wcmc_sg_py), nrow(wcmc_sg_pt))) +
  map_theme + coord_sf(crs = 27700)

ggsave("figures/habitats/seagrass_gems_wcmc.png", p_wcmc_gems,
       width = 14, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/habitats/seagrass_gems_wcmc.png\n")


## ---- 1c. Indicator species occurrence records (Zostera) ----

ind_path <- "outputs/species/indicator_records.gpkg"
ind_zostera <- st_sf(geometry = st_sfc(crs = 27700))

if (file.exists(ind_path)) {
  ind_zostera <- st_read(ind_path, quiet = TRUE) %>%
    filter(grepl("Zostera", latin_name, ignore.case = TRUE)) %>%
    st_transform(27700) %>%
    st_make_valid() %>%
    st_intersection(scamp_aoi_sf) %>%
    mutate(source = "Indicator records") %>%
    select(source, geom)
  cat(sprintf("  Indicator Zostera records in AOI: %d\n", nrow(ind_zostera)))
} else {
  cat("  WARNING: indicator_records.gpkg not found — skipping\n")
}

p_indicator_pts <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
          linewidth = 0.3) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  geom_sf(data = zoom_boxes, fill = NA, colour = "red", linewidth = 0.6,
          linetype = "dashed") +
  geom_sf(data = gems_sg_pt, aes(colour = "GeMS point"), shape = 16,
          size = 2, alpha = 0.7) +
  geom_sf(data = wcmc_sg_pt, aes(colour = "WCMC point"), shape = 17,
          size = 2.5, alpha = 0.7) +
  {if (nrow(ind_zostera) > 0)
    geom_sf(data = ind_zostera, aes(colour = "Indicator record"), shape = 4,
            size = 2.5, alpha = 0.8, stroke = 1)} +
  scale_colour_manual(
    values = c("GeMS point" = "#d95f02", "WCMC point" = "#a6761d",
               "Indicator record" = "#e6ab02"),
    name = "Source") +
  labs(title = "Seagrass points — GeMS vs WCMC vs Indicator records",
       subtitle = sprintf("GeMS: %d | WCMC: %d | Indicator: %d",
                          nrow(gems_sg_pt), nrow(wcmc_sg_pt), nrow(ind_zostera))) +
  map_theme + coord_sf(crs = 27700)

ggsave("figures/habitats/seagrass_indicator_vs_gems_wcmc_pts.png", p_indicator_pts,
       width = 14, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/habitats/seagrass_indicator_vs_gems_wcmc_pts.png\n")

## Cluster indicator records (same method as GeMS/WCMC: 2km buffer, union, decompose)
ind_buf <- st_buffer(ind_zostera, dist = 2000) %>% st_union() %>%
  st_cast("POLYGON")
ind_cluster_ids <- st_intersects(ind_zostera, ind_buf)
ind_zostera$cluster <- sapply(ind_cluster_ids, `[`, 1)

ind_cluster_stats <- ind_zostera %>%
  st_drop_geometry() %>%
  group_by(cluster) %>%
  summarise(n = n(), .groups = "drop")

cat("\n  -- Indicator cluster statistics --\n")
print(ind_cluster_stats)
cat(sprintf("  Indicator clusters: %d\n", nrow(ind_cluster_stats)))

## Build bboxes per indicator cluster
ind_cluster_boxes <- do.call(rbind, lapply(sort(unique(ind_zostera$cluster)), function(cl) {
  bb <- st_bbox(ind_zostera[ind_zostera$cluster == cl, ])
  st_sf(
    cluster = cl,
    n = sum(ind_zostera$cluster == cl),
    geometry = st_sfc(st_polygon(list(matrix(c(
      bb["xmin"], bb["ymin"],
      bb["xmax"], bb["ymin"],
      bb["xmax"], bb["ymax"],
      bb["xmin"], bb["ymax"],
      bb["xmin"], bb["ymin"]
    ), ncol = 2, byrow = TRUE))), crs = 27700)
  )
}))

## Filter to NEW indicator clusters only (no overlap with existing zoom_boxes)
ind_new <- ind_cluster_boxes[lengths(st_intersects(ind_cluster_boxes, zoom_boxes)) == 0, ]

## Sort new clusters west-to-east by centroid easting
ind_new <- ind_new[order(st_coordinates(st_centroid(ind_new))[, 1]), ]
ind_new$rank <- seq_len(nrow(ind_new))

ind_new_coords <- st_coordinates(st_centroid(ind_new))
cat(sprintf("  New indicator clusters (not overlapping WCMC boxes): %d\n", nrow(ind_new)))
cat("  New clusters west to east:\n")
for (j in seq_len(nrow(ind_new))) {
  cat(sprintf("    rank %d: E=%.0f  N=%.0f  (n=%d)\n",
              j, ind_new_coords[j, 1], ind_new_coords[j, 2], ind_new$n[j]))
}

## Helper: build bbox polygon from xmin/ymin/xmax/ymax
make_box <- function(xmin, ymin, xmax, ymax, name) {
  st_sf(name = name,
        geometry = st_sfc(st_polygon(list(matrix(c(
          xmin, ymin, xmax, ymin, xmax, ymax, xmin, ymax, xmin, ymin
        ), ncol = 2, byrow = TRUE))), crs = 27700))
}

## Rank 1 = Wigtown Bay: expand ~2000m all directions
bb1 <- st_bbox(ind_new[ind_new$rank == 1, ])
box_wigtown <- make_box(bb1["xmin"] - 2000, bb1["ymin"] - 2000,
                        bb1["xmax"] + 2000, bb1["ymax"] + 2000, "Wigtown Bay")

## Rank 2 = Fleet Bay: expand ~2000m all directions
bb2 <- st_bbox(ind_new[ind_new$rank == 2, ])
box_fleet <- make_box(bb2["xmin"] - 2000, bb2["ymin"] - 2000,
                      bb2["xmax"] + 2000, bb2["ymax"] + 2000, "Fleet Bay")

## Rank 3 = merge into Auchencairn: union of existing bbox + cluster 3, expand south 2000m
bb_auch_old <- st_bbox(zoom_boxes[zoom_boxes$name == "Auchencairn Bay & Rough Island", ])
bb3 <- st_bbox(ind_new[ind_new$rank == 3, ])
box_auch <- make_box(min(bb_auch_old["xmin"], bb3["xmin"]),
                     min(bb_auch_old["ymin"], bb3["ymin"]) - 2000,
                     max(bb_auch_old["xmax"], bb3["xmax"]),
                     max(bb_auch_old["ymax"], bb3["ymax"]),
                     "Auchencairn Bay & Rough Island")

## Kirkcudbright Bay: expand southward ~1000m, westward ~1000m
bb_kirk <- st_bbox(zoom_boxes[zoom_boxes$name == "Kirkcudbright Bay", ])
box_kirk <- make_box(bb_kirk["xmin"] - 1000, bb_kirk["ymin"] - 1000,
                     bb_kirk["xmax"], bb_kirk["ymax"], "Kirkcudbright Bay")

## Ranks 4+ = Nith Estuary: combined bbox of all remaining new clusters
bb_nith <- st_bbox(ind_new[ind_new$rank >= 4, ])
box_nith <- make_box(bb_nith["xmin"], bb_nith["ymin"],
                     bb_nith["xmax"], bb_nith["ymax"], "Nith Estuary")

## Loch Ryan: keep as-is
box_loch <- zoom_boxes[zoom_boxes$name == "Loch Ryan", ]

all_boxes <- rbind(box_loch, box_wigtown, box_fleet, box_kirk, box_auch, box_nith)
cat(sprintf("  Final cluster boxes: %d\n", nrow(all_boxes)))

p_all_bboxes <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
          linewidth = 0.3) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  geom_sf(data = gems_sg_pt, aes(colour = "GeMS point"), shape = 16,
          size = 2, alpha = 0.7) +
  geom_sf(data = wcmc_sg_pt, aes(colour = "WCMC point"), shape = 17,
          size = 2.5, alpha = 0.7) +
  {if (nrow(ind_zostera) > 0)
    geom_sf(data = ind_zostera, aes(colour = "Indicator record"), shape = 4,
            size = 2.5, alpha = 0.8, stroke = 1)} +
  scale_colour_manual(
    values = c("GeMS point" = "#d95f02", "WCMC point" = "#a6761d",
               "Indicator record" = "#e6ab02"),
    name = "Source") +
  geom_sf(data = all_boxes, fill = NA, colour = "red", linewidth = 0.8) +
  geom_sf_text(data = all_boxes, aes(label = name), colour = "red",
               size = 3, fontface = "bold") +
  labs(title = "Seagrass cluster bounding boxes — all sources",
       subtitle = sprintf("%d clusters defined", nrow(all_boxes))) +
  map_theme + coord_sf(crs = 27700)

ggsave("figures/habitats/seagrass_all_cluster_bboxes.png", p_all_bboxes,
       width = 14, height = 10, units = "in", dpi = 200)
cat("  Saved: figures/habitats/seagrass_all_cluster_bboxes.png\n")

## Update zoom_boxes for downstream use
zoom_boxes <- all_boxes

## Individual cluster zoom plots
for (i in seq_len(nrow(all_boxes))) {
  box_i <- all_boxes[i, ]
  bb_i  <- st_bbox(box_i)

  gems_py_i  <- gems_sg_py[lengths(st_intersects(gems_sg_py, box_i)) > 0, ]
  gems_pt_i  <- gems_sg_pt[lengths(st_intersects(gems_sg_pt, box_i)) > 0, ]
  wcmc_py_i  <- wcmc_sg_py[lengths(st_intersects(wcmc_sg_py, box_i)) > 0, ]
  wcmc_pt_i  <- wcmc_sg_pt[lengths(st_intersects(wcmc_sg_pt, box_i)) > 0, ]
  ind_i      <- ind_zostera[lengths(st_intersects(ind_zostera, box_i)) > 0, ]

  p_i <- ggplot() +
    geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
    geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
            linewidth = 0.3) +
    geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
    {if (nrow(gems_py_i) > 0)
      geom_sf(data = gems_py_i, aes(fill = "GeMS polygon"), colour = NA, alpha = 0.5)} +
    {if (nrow(wcmc_py_i) > 0)
      geom_sf(data = wcmc_py_i, aes(fill = "WCMC polygon"), colour = NA, alpha = 0.5)} +
    {if (nrow(gems_pt_i) > 0)
      geom_sf(data = gems_pt_i, aes(colour = "GeMS point"), shape = 16,
              size = 3, alpha = 0.8)} +
    {if (nrow(wcmc_pt_i) > 0)
      geom_sf(data = wcmc_pt_i, aes(colour = "WCMC point"), shape = 17,
              size = 3, alpha = 0.8)} +
    {if (nrow(ind_i) > 0)
      geom_sf(data = ind_i, aes(colour = "Indicator record"), shape = 4,
              size = 3, alpha = 0.8, stroke = 1)} +
    scale_fill_manual(
      values = c("GeMS polygon" = "#1b9e77", "WCMC polygon" = "#7570b3"),
      name = "Polygons", drop = FALSE) +
    scale_colour_manual(
      values = c("GeMS point" = "#d95f02", "WCMC point" = "#a6761d",
                 "Indicator record" = "#e6ab02"),
      name = "Points", drop = FALSE) +
    coord_sf(xlim = c(bb_i["xmin"], bb_i["xmax"]),
             ylim = c(bb_i["ymin"], bb_i["ymax"]), crs = 27700) +
    labs(title = box_i$name,
         subtitle = sprintf("GeMS poly: %d | GeMS pt: %d | WCMC poly: %d | WCMC pt: %d | Indicator: %d",
                            nrow(gems_py_i), nrow(gems_pt_i), nrow(wcmc_py_i),
                            nrow(wcmc_pt_i), nrow(ind_i))) +
    map_theme

  fname <- paste0("figures/habitats/seagrass_cluster_",
                  gsub("[^a-z0-9]+", "_", tolower(box_i$name)), ".png")
  ggsave(fname, p_i, width = 12, height = 10, units = "in", dpi = 300)
  cat(sprintf("  Saved: %s\n", fname))
}

## Combine GeMS + WCMC + Indicator into unified polygon and point layers
sg_polygons <- rbind(
  gems_sg_py %>% select(source, geometry),
  wcmc_sg_py %>% select(source, geometry)
)

ind_pts <- ind_zostera %>% select(source, geom)
st_geometry(ind_pts) <- "geometry"

sg_points <- rbind(
  gems_sg_pt %>% select(source, geometry),
  wcmc_sg_pt %>% select(source, geometry),
  ind_pts
)

## Assign each feature to a cluster (st_within for points, st_intersects for polygons)
pt_hits <- st_within(sg_points, zoom_boxes)
sg_points$cluster <- zoom_boxes$name[sapply(pt_hits, function(x) if (length(x) == 0) NA_integer_ else x[1])]

py_hits <- st_intersects(sg_polygons, zoom_boxes)
sg_polygons$cluster <- zoom_boxes$name[sapply(py_hits, function(x) if (length(x) == 0) NA_integer_ else x[1])]

cat(sprintf("\n  Combined seagrass polygons: %d (GeMS: %d, WCMC: %d)\n",
            nrow(sg_polygons), nrow(gems_sg_py), nrow(wcmc_sg_py)))
cat(sprintf("  Combined seagrass points:   %d (GeMS: %d, WCMC: %d, Indicator: %d)\n",
            nrow(sg_points), nrow(gems_sg_pt), nrow(wcmc_sg_pt), nrow(ind_zostera)))
cat(sprintf("  Points unassigned to cluster: %d\n", sum(is.na(sg_points$cluster))))
cat(sprintf("  Polygons unassigned to cluster: %d\n", sum(is.na(sg_polygons$cluster))))

cat("\n  -- Points per cluster and source --\n")
print(sg_points %>% st_drop_geometry() %>% count(cluster, source))
cat("\n  -- Polygons per cluster and source --\n")
print(sg_polygons %>% st_drop_geometry() %>% count(cluster, source))


## ---- 1d. SFP Zostera polygons ----

sfp_dir  <- "data/raw/habitats/SFP Zostera polygons"
sfp_shps <- list.files(sfp_dir, pattern = "\\.shp$",
                        recursive = TRUE, full.names = TRUE)
sfp_shps <- sfp_shps[!grepl("__MACOSX", sfp_shps)]

cat(sprintf("  SFP shapefiles found: %d\n", length(sfp_shps)))

sfp_list <- list()
for (shp in sfp_shps) {
  z <- tryCatch(
    st_read(shp, quiet = TRUE) %>%
      st_transform(27700) %>%
      st_make_valid(),
    error = function(e) {
      cat(sprintf("    SKIP: %s — %s\n", basename(shp), e$message))
      NULL
    }
  )
  if (!is.null(z) && nrow(z) > 0) {
    z$site   <- basename(dirname(shp))
    z$file   <- tools::file_path_sans_ext(basename(shp))
    z$species <- case_when(
      grepl("marina", z$file, ignore.case = TRUE) &
        grepl("noltei|noltii", z$file, ignore.case = TRUE) ~ "Z. marina & Z. noltei",
      grepl("marina", z$file, ignore.case = TRUE) ~ "Z. marina",
      grepl("noltei|noltii", z$file, ignore.case = TRUE) ~ "Z. noltei",
      TRUE ~ "Zostera sp."
    )
    sfp_list[[shp]] <- z
  }
}

sfp_all <- bind_rows(sfp_list) %>%
  st_intersection(scamp_aoi_sf) %>%
  mutate(source = "SFP survey") %>%
  select(source, site, species, geometry)

cat(sprintf("  SFP Zostera polygons in AOI: %d across %d sites\n",
    nrow(sfp_all), length(unique(sfp_all$site))))

## SFP polygon centroids with arrows
sfp_centroids <- sfp_all %>%
  mutate(area_m2 = as.numeric(st_area(geometry))) %>%
  group_by(site) %>%
  summarise(
    geometry = st_union(geometry),
    n = n(),
    total_ha = sum(area_m2) / 10000,
    .groups = "drop"
  ) %>%
  st_centroid() %>%
  bind_cols(st_coordinates(.) %>% as.data.frame())

aoi_ctr <- st_coordinates(st_centroid(scamp_aoi_sf))
sfp_centroids <- sfp_centroids %>%
  mutate(
    dx = X - aoi_ctr[1, "X"],
    dy = Y - aoi_ctr[1, "Y"],
    len = sqrt(dx^2 + dy^2),
    nudge_x = (dx / len) * 15000,
    nudge_y = (dy / len) * 15000,
    label_x = X + nudge_x,
    label_y = Y + nudge_y,
    label = sprintf("%s\nn=%d (%.2f ha)", site, n, total_ha)
  )

p_sfp <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
          linewidth = 0.3) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  geom_sf(data = zoom_boxes, fill = NA, colour = "red", linewidth = 0.6,
          linetype = "dashed") +
  geom_sf(data = sfp_all, aes(fill = species), colour = NA, alpha = 0.6) +
  scale_fill_brewer(palette = "Set2", name = "Species") +
  geom_segment(data = sfp_centroids,
               aes(x = label_x, y = label_y, xend = X, yend = Y),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
               colour = "grey20", linewidth = 0.4) +
  geom_label(data = sfp_centroids,
             aes(x = label_x, y = label_y, label = label),
             size = 2.2, fill = alpha("white", 0.7), colour = "grey20",
             label.padding = unit(0.2, "lines"), fontface = "bold") +
  labs(title = "SFP Zostera survey polygons (SCAMP AOI)",
       subtitle = sprintf("%d polygons across %d sites",
                          nrow(sfp_all), length(unique(sfp_all$site)))) +
  map_theme + coord_sf(crs = 27700)

ggsave("figures/habitats/seagrass_sfp_aoi.png", p_sfp,
       width = 14, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/habitats/seagrass_sfp_aoi.png\n")

## Per-cluster comparison: GeMS + WCMC + SFP polygons
for (i in seq_len(nrow(zoom_boxes))) {
  box_i <- zoom_boxes[i, ]
  bb_i  <- st_bbox(box_i)

  gems_py_i <- gems_sg_py[lengths(st_intersects(gems_sg_py, box_i)) > 0, ]
  wcmc_py_i <- wcmc_sg_py[lengths(st_intersects(wcmc_sg_py, box_i)) > 0, ]
  sfp_py_i  <- sfp_all[lengths(st_intersects(sfp_all, box_i)) > 0, ]

  p_i <- ggplot() +
    geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
    geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
            linewidth = 0.3) +
    geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
    {if (nrow(gems_py_i) > 0)
      geom_sf(data = gems_py_i, aes(fill = "GeMS"), colour = NA, alpha = 0.5)} +
    {if (nrow(wcmc_py_i) > 0)
      geom_sf(data = wcmc_py_i, aes(fill = "WCMC"), colour = NA, alpha = 0.5)} +
    {if (nrow(sfp_py_i) > 0)
      geom_sf(data = sfp_py_i, aes(fill = "SFP survey"), colour = NA, alpha = 0.5)} +
    scale_fill_manual(
      values = c("GeMS" = "#1b9e77", "WCMC" = "#7570b3", "SFP survey" = "#d95f02"),
      name = "Source") +
    coord_sf(xlim = c(bb_i["xmin"], bb_i["xmax"]),
             ylim = c(bb_i["ymin"], bb_i["ymax"]), crs = 27700) +
    labs(title = sprintf("%s — polygon comparison", box_i$name),
         subtitle = sprintf("GeMS: %d | WCMC: %d | SFP: %d",
                            nrow(gems_py_i), nrow(wcmc_py_i), nrow(sfp_py_i))) +
    map_theme

  fname <- paste0("figures/habitats/seagrass_polygons_",
                  gsub("[^a-z0-9]+", "_", tolower(box_i$name)), ".png")
  ggsave(fname, p_i, width = 12, height = 10, units = "in", dpi = 300)
  cat(sprintf("  Saved: %s\n", fname))
}

## Per-cluster comparison: SFP polygons + all points
for (i in seq_len(nrow(zoom_boxes))) {
  box_i <- zoom_boxes[i, ]
  bb_i  <- st_bbox(box_i)

  sfp_py_i  <- sfp_all[lengths(st_intersects(sfp_all, box_i)) > 0, ]
  gems_pt_i <- gems_sg_pt[lengths(st_intersects(gems_sg_pt, box_i)) > 0, ]
  wcmc_pt_i <- wcmc_sg_pt[lengths(st_intersects(wcmc_sg_pt, box_i)) > 0, ]
  ind_i     <- ind_zostera[lengths(st_intersects(ind_zostera, box_i)) > 0, ]

  p_i <- ggplot() +
    geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
    geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
            linewidth = 0.3) +
    geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
    {if (nrow(sfp_py_i) > 0)
      geom_sf(data = sfp_py_i, aes(fill = "SFP survey"), colour = NA, alpha = 0.5)} +
    scale_fill_manual(values = c("SFP survey" = "#d95f02"), name = "Polygons") +
    {if (nrow(gems_pt_i) > 0)
      geom_sf(data = gems_pt_i, aes(colour = "GeMS point"), shape = 16,
              size = 3, alpha = 0.8)} +
    {if (nrow(wcmc_pt_i) > 0)
      geom_sf(data = wcmc_pt_i, aes(colour = "WCMC point"), shape = 17,
              size = 3, alpha = 0.8)} +
    {if (nrow(ind_i) > 0)
      geom_sf(data = ind_i, aes(colour = "Indicator record"), shape = 4,
              size = 3, alpha = 0.8, stroke = 1)} +
    scale_colour_manual(
      values = c("GeMS point" = "#1b9e77", "WCMC point" = "#7570b3",
                 "Indicator record" = "#e6ab02"),
      name = "Points") +
    coord_sf(xlim = c(bb_i["xmin"], bb_i["xmax"]),
             ylim = c(bb_i["ymin"], bb_i["ymax"]), crs = 27700) +
    labs(title = sprintf("%s — SFP polygons + points", box_i$name),
         subtitle = sprintf("SFP poly: %d | GeMS pt: %d | WCMC pt: %d | Indicator: %d",
                            nrow(sfp_py_i), nrow(gems_pt_i), nrow(wcmc_pt_i),
                            nrow(ind_i))) +
    map_theme

  fname <- paste0("figures/habitats/seagrass_sfp_pts_",
                  gsub("[^a-z0-9]+", "_", tolower(box_i$name)), ".png")
  ggsave(fname, p_i, width = 12, height = 10, units = "in", dpi = 300)
  cat(sprintf("  Saved: %s\n", fname))
}


## ---- 1e. EA Seagrass Potential ----

ea_sg <- st_read(
  "data/raw/habitats/Seagrass_Potential.gpkg/Seagrass_Potential.gpkg",
  quiet = TRUE
) %>%
  st_transform(27700) %>%
  st_make_valid() %>%
  st_intersection(scamp_aoi_sf) %>%
  mutate(source = "EA potential") %>%
  select(source, shape)

cat(sprintf("  EA seagrass potential sites in AOI: %d\n", nrow(ea_sg)))


## ---- 1f. INTERMEDIATE MAP — all seagrass sources ----

cat("\n  Building intermediate seagrass map...\n")

## Collect all polygon sources
all_sg_py <- bind_rows(
  gems_sg_py,
  wcmc_sg_py,
  sfp_all %>% select(source, geometry)
)

## Collect all point sources
all_sg_pt <- bind_rows(
  gems_sg_pt,
  wcmc_sg_pt,
  ind_zostera
)

source_colours <- c(
  "GeMS habitat polygon" = "#1b9e77",
  "WCMC polygon"         = "#7570b3",
  "SFP survey"           = "#e7298a",
  "EA potential"         = "#66a61e",
  "GeMS habitat point"   = "#d95f02",
  "WCMC point"           = "#a6761d",
  "Indicator records"    = "#e6ab02"
)

p_seagrass_sources <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50",
          linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
          linewidth = 0.3) +
  geom_sf(data = scamp_coastline, colour = "grey30", linewidth = 0.3) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  geom_sf(data = ea_sg, aes(fill = source), colour = NA, alpha = 0.3) +
  geom_sf(data = all_sg_py, aes(fill = source), colour = NA, alpha = 0.6) +
  geom_sf(data = all_sg_pt, aes(colour = source), size = 1.2, alpha = 0.7) +
  scale_fill_manual(values = source_colours, name = "Source (polygons)") +
  scale_colour_manual(values = source_colours, name = "Source (points)") +
  annotation_scale(location = "bl", width_hint = 0.2) +
  annotation_north_arrow(location = "tl", style = north_arrow_minimal(),
                         height = unit(1, "cm"), width = unit(1, "cm")) +
  labs(title = "Seagrass — All data sources (SCAMP AOI)",
       subtitle = sprintf("Polygons: %d | Points: %d | EA potential: %d",
                           nrow(all_sg_py), nrow(all_sg_pt), nrow(ea_sg))) +
  map_theme +
  coord_sf(crs = 27700)

ggsave("figures/habitats/seagrass_all_sources.png", p_seagrass_sources,
       width = 14, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/habitats/seagrass_all_sources.png\n")


## ---- 1g. HARMONISED SEAGRASS LAYER ----

cat("\n  Harmonising seagrass extents...\n")

## SFP is the highest-confidence polygon source — it takes priority

## 1) Remove points that overlap, touch or are within 100m of SFP polygons
sfp_union <- st_union(sfp_all) %>% st_make_valid()
sfp_buf100 <- st_buffer(sfp_union, dist = 100)
pts_near_sfp <- lengths(st_intersects(sg_points, sfp_buf100)) > 0
cat(sprintf("  Points within 100m of SFP polygons (removed): %d\n", sum(pts_near_sfp)))
sg_points_clean <- sg_points[!pts_near_sfp, ]

## 2) Remove GeMS/WCMC polygons that overlap SFP polygons (SFP takes priority)
other_polys <- sg_polygons
other_overlap_sfp <- lengths(st_intersects(other_polys, sfp_union)) > 0
cat(sprintf("  GeMS/WCMC polygons overlapping SFP (removed): %d\n", sum(other_overlap_sfp)))
other_polys_clean <- other_polys[!other_overlap_sfp, ]

## 3) Flag potentially unmapped sites: GeMS/WCMC polygons that are near points
##    but do NOT overlap SFP
if (nrow(other_polys_clean) > 0 && nrow(sg_points_clean) > 0) {
  other_buf <- st_buffer(other_polys_clean, dist = 100)
  pts_near_other <- lengths(st_intersects(sg_points_clean, other_buf)) > 0
  unmapped_flag_pts <- sg_points_clean[pts_near_other, ] %>%
    mutate(flag = "potential_unmapped_site")
  cat(sprintf("  Points near non-SFP polygons (potential unmapped sites): %d\n",
              nrow(unmapped_flag_pts)))
} else {
  unmapped_flag_pts <- sg_points_clean[0, ] %>% mutate(flag = character(0))
  cat("  Points near non-SFP polygons (potential unmapped sites): 0\n")
}

## 4) Build harmonised polygon layer: SFP + surviving GeMS/WCMC
harm_polys <- rbind(
  sfp_all %>% select(source, geometry),
  other_polys_clean %>% select(source, geometry)
)

cat(sprintf("\n  -- Harmonised polygon layer --\n"))
cat(sprintf("  SFP polygons: %d\n", nrow(sfp_all)))
cat(sprintf("  Surviving GeMS/WCMC polygons: %d\n", nrow(other_polys_clean)))
cat(sprintf("  Total harmonised polygons: %d\n", nrow(harm_polys)))

cat(sprintf("\n  -- Harmonised point layer (DRAFT) --\n"))
cat(sprintf("  Original points: %d\n", nrow(sg_points)))
cat(sprintf("  After removing SFP overlaps: %d\n", nrow(sg_points_clean)))
cat(sprintf("  Flagged as potential unmapped: %d\n", nrow(unmapped_flag_pts)))

## 5) Mark flagged points in the clean layer
sg_points_clean$flag <- NA_character_
if (nrow(unmapped_flag_pts) > 0) {
  flag_idx <- lengths(st_intersects(sg_points_clean,
                                     st_buffer(other_polys_clean, dist = 100))) > 0
  sg_points_clean$flag[flag_idx] <- "potential_unmapped_site"
}

## 6) Save harmonised outputs
st_write(harm_polys,
         "data/processed/seagrass_polygons_harmonised.gpkg",
         delete_dsn = TRUE, quiet = TRUE)
cat("  Saved: data/processed/seagrass_polygons_harmonised.gpkg\n")

st_write(sg_points_clean,
         "data/processed/seagrass_points_DRAFT.gpkg",
         delete_dsn = TRUE, quiet = TRUE)
cat("  Saved: data/processed/seagrass_points_DRAFT.gpkg\n")

if (nrow(ea_sg) > 0) {
  st_write(ea_sg,
           "data/processed/seagrass_ea_potential.gpkg",
           delete_dsn = TRUE, quiet = TRUE)
  cat("  Saved: data/processed/seagrass_ea_potential.gpkg\n")
}


## ---- 1h. HARMONISED MAPS — per cluster ----

for (i in seq_len(nrow(zoom_boxes))) {
  box_i <- zoom_boxes[i, ]
  bb_i  <- st_bbox(box_i)

  harm_py_i <- harm_polys[lengths(st_intersects(harm_polys, box_i)) > 0, ]
  pts_i     <- sg_points_clean[lengths(st_intersects(sg_points_clean, box_i)) > 0, ]

  sfp_i    <- harm_py_i[harm_py_i$source == "SFP survey", ]
  unmapped_i <- harm_py_i[harm_py_i$source != "SFP survey", ]

  p_i <- ggplot() +
    geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
    geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
            linewidth = 0.3) +
    geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
    {if (nrow(sfp_i) > 0)
      geom_sf(data = sfp_i, aes(fill = "Known extent"),
              colour = NA, alpha = 0.6)} +
    {if (nrow(unmapped_i) > 0)
      geom_sf(data = unmapped_i, aes(fill = "Potentially unmapped"),
              colour = NA, alpha = 0.6)} +
    scale_fill_manual(
      values = c("Known extent" = "#1b9e77", "Potentially unmapped" = "#e41a1c"),
      name = "Polygons") +
    {if (nrow(pts_i) > 0)
      geom_sf(data = pts_i, aes(colour = source, shape = source),
              size = 2.5, alpha = 0.8)} +
    scale_colour_manual(
      values = c("GeMS habitat point" = "#d95f02", "WCMC point" = "#a6761d",
                 "Indicator records" = "#e6ab02"),
      name = "Points (draft)") +
    scale_shape_manual(
      values = c("GeMS habitat point" = 16, "WCMC point" = 17,
                 "Indicator records" = 4),
      name = "Points (draft)") +
    coord_sf(xlim = c(bb_i["xmin"], bb_i["xmax"]),
             ylim = c(bb_i["ymin"], bb_i["ymax"]), crs = 27700) +
    labs(title = sprintf("%s — harmonised seagrass", box_i$name),
         subtitle = sprintf("Known: %d | Potentially unmapped: %d | Points (draft): %d",
                            nrow(sfp_i), nrow(unmapped_i), nrow(pts_i))) +
    map_theme

  fname <- paste0("figures/habitats/seagrass_harmonised_",
                  gsub("[^a-z0-9]+", "_", tolower(box_i$name)), ".png")
  ggsave(fname, p_i, width = 12, height = 10, units = "in", dpi = 300)
  cat(sprintf("  Saved: %s\n", fname))
}


##=====================================================================
## INVESTIGATION — Diagnostics for seagrass layers
##=====================================================================

cat("\n\n=== INVESTIGATION: Seagrass layer diagnostics ===\n")

cat("\n-- GeMS habitat polygons (seagrass, AOI) --\n")
cat(sprintf("  Rows: %d | CRS: %s\n", nrow(gems_hab_py), st_crs(gems_hab_py)$input))
if (nrow(gems_hab_py) > 0) {
  cat("  Columns:", paste(names(gems_hab_py), collapse = ", "), "\n")
  cat("  Area range (m²):",
      sprintf("%.1f — %.1f\n",
              min(as.numeric(st_area(gems_hab_py))),
              max(as.numeric(st_area(gems_hab_py)))))
}

cat("\n-- GeMS habitat points (seagrass, AOI) --\n")
cat(sprintf("  Rows: %d\n", nrow(gems_hab_pt)))

cat("\n-- WCMC seagrass polygons (AOI) --\n")
cat(sprintf("  Rows: %d\n", nrow(wcmc_sg_py)))

cat("\n-- WCMC seagrass points (AOI) --\n")
cat(sprintf("  Rows: %d\n", nrow(wcmc_sg_pt)))

cat("\n-- Indicator Zostera records (AOI) --\n")
cat(sprintf("  Rows: %d\n", nrow(ind_zostera)))

cat("\n-- SFP Zostera polygons (AOI) --\n")
cat(sprintf("  Rows: %d across %d sites\n",
    nrow(sfp_all), length(unique(sfp_all$site))))
if (nrow(sfp_all) > 0) {
  cat("  By site:\n")
  print(table(sfp_all$site, sfp_all$species))
  cat(sprintf("  Total area: %.1f m²\n",
      sum(as.numeric(st_area(sfp_all)))))
}

cat("\n-- EA seagrass potential (AOI) --\n")
cat(sprintf("  Rows: %d\n", nrow(ea_sg)))

cat("\n-- Harmonised outputs --\n")
if (length(known_extent) > 0 && !st_is_empty(known_extent)) {
  cat(sprintf("  Known extent total area: %.1f m² (%.2f ha)\n",
      as.numeric(st_area(known_extent)),
      as.numeric(st_area(known_extent)) / 10000))
}
cat(sprintf("  Expansion candidate points: %d\n", nrow(expansion_pts)))
cat(sprintf("  EA potential sites: %d\n", nrow(ea_sg)))

cat("\n=== Seagrass section complete ===\n")
