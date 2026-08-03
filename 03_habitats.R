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
  "data/raw/habitats/Seagrass/WCMC013-014_SeagrassPtPy2021_v7_1/014_001_WCMC013-014_SeagrassPtPy2021_v7_1/01_Data/WCMC013_014_Seagrasses_Py_v7_1.shp",
  quiet = TRUE
) %>%
  st_transform(27700) %>%
  st_make_valid() %>%
  st_intersection(scamp_aoi_sf) %>%
  mutate(source = "WCMC polygon") %>%
  select(source, geometry)

wcmc_sg_pt <- st_read(
  "data/raw/habitats/Seagrass/WCMC013-014_SeagrassPtPy2021_v7_1/014_001_WCMC013-014_SeagrassPtPy2021_v7_1/01_Data/WCMC013_014_Seagrasses_Pt_v7_1.shp",
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

## Build zoom boxes from GeMS + WCMC polygon clusters (data-driven, not hardcoded)
## Combine all seagrass polygons with their cluster IDs
gems_cluster_polys <- gems_sg_py %>%
  mutate(src = "GeMS") %>%
  group_by(cluster) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  mutate(src = "GeMS")

wcmc_cluster_polys <- wcmc_sg_py %>%
  mutate(src = "WCMC") %>%
  group_by(cluster) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  mutate(src = "WCMC")

## Build bbox per cluster with 2km buffer
all_cluster_polys <- rbind(gems_cluster_polys, wcmc_cluster_polys)

zoom_boxes <- do.call(rbind, lapply(seq_len(nrow(all_cluster_polys)), function(i) {
  bb <- st_bbox(st_buffer(all_cluster_polys[i, ], dist = 2000))
  st_sf(
    name = sprintf("Cluster_%s_%d", all_cluster_polys$src[i],
                   all_cluster_polys$cluster[i]),
    geometry = st_sfc(st_polygon(list(matrix(c(
      bb["xmin"], bb["ymin"],
      bb["xmax"], bb["ymin"],
      bb["xmax"], bb["ymax"],
      bb["xmin"], bb["ymax"],
      bb["xmin"], bb["ymin"]
    ), ncol = 2, byrow = TRUE))), crs = 27700)
  )
}))

## Merge overlapping boxes (clusters that are near each other)
zoom_buf <- st_buffer(zoom_boxes, dist = 1000) %>% st_union() %>%
  st_cast("POLYGON")
zoom_group <- sapply(st_intersects(zoom_boxes, zoom_buf), `[`, 1)

zoom_boxes <- do.call(rbind, lapply(sort(unique(zoom_group)), function(g) {
  members <- zoom_boxes[zoom_group == g, ]
  bb <- st_bbox(members)
  st_sf(
    name = sprintf("Cluster %d", g),
    geometry = st_sfc(st_polygon(list(matrix(c(
      bb["xmin"], bb["ymin"],
      bb["xmax"], bb["ymin"],
      bb["xmax"], bb["ymax"],
      bb["xmin"], bb["ymax"],
      bb["xmin"], bb["ymin"]
    ), ncol = 2, byrow = TRUE))), crs = 27700)
  )
}))

cat(sprintf("  Data-driven zoom boxes: %d\n", nrow(zoom_boxes)))

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

## Expand each new indicator cluster box by 2km and add to zoom_boxes
for (j in seq_len(nrow(ind_new))) {
  bb_j <- st_bbox(st_buffer(ind_new[j, ], dist = 2000))
  new_box <- make_box(bb_j["xmin"], bb_j["ymin"], bb_j["xmax"], bb_j["ymax"],
                      sprintf("Indicator_%d", j))
  zoom_boxes <- rbind(zoom_boxes, new_box)
}

## Re-merge overlapping boxes after adding indicator clusters
zb_buf <- st_buffer(zoom_boxes, dist = 1000) %>% st_union() %>%
  st_cast("POLYGON")
zb_group <- sapply(st_intersects(zoom_boxes, zb_buf), `[`, 1)

all_boxes <- do.call(rbind, lapply(sort(unique(zb_group)), function(g) {
  members <- zoom_boxes[zb_group == g, ]
  bb <- st_bbox(members)
  st_sf(
    name = sprintf("Cluster %d", g),
    geometry = st_sfc(st_polygon(list(matrix(c(
      bb["xmin"], bb["ymin"],
      bb["xmax"], bb["ymin"],
      bb["xmax"], bb["ymax"],
      bb["xmin"], bb["ymax"],
      bb["xmin"], bb["ymin"]
    ), ncol = 2, byrow = TRUE))), crs = 27700)
  )
}))

## Assign geographic names by centroid easting (west to east)
all_boxes <- all_boxes[order(st_coordinates(st_centroid(all_boxes))[, 1]), ]
place_names <- c("Loch Ryan", "Wigtown Bay", "Fleet Bay","Kirkcudbright Bay",
                 "Auchencairn Bay & Rough Island", "Nith Estuary")
if (nrow(all_boxes) <= length(place_names)) {
  all_boxes$name <- place_names[seq_len(nrow(all_boxes))]
} else {
  all_boxes$name <- c(place_names,
    sprintf("Cluster %d", seq(length(place_names) + 1, nrow(all_boxes))))
}

cat(sprintf("  Final cluster boxes: %d\n", nrow(all_boxes)))
for (k in seq_len(nrow(all_boxes))) {
  cat(sprintf("    %s\n", all_boxes$name[k]))
}

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

sfp_dir  <- "data/raw/habitats/Seagrass/SFP Zostera polygons"
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
  "data/raw/habitats/Seagrass/Seagrass_Potential.gpkg/Seagrass_Potential.gpkg",
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


## ---- 1i. SOLWAY SEAGRASS POTENTIAL (from physical model) ----

cat("\n  Reading Solway seagrass potential model output...\n")

solway_sg_path <- "outputs/habitats/seagrass_potential_solway.gpkg"
stopifnot("Solway seagrass potential not found — run 04_physical.R section 8 first" =
            file.exists(solway_sg_path))

solway_sg <- st_read(solway_sg_path, quiet = TRUE) %>% st_make_valid()
cat(sprintf("  Solway potential polygons: %d | Total area: %.1f ha\n",
            nrow(solway_sg),
            sum(as.numeric(st_area(solway_sg))) / 10000))


## ---- 1j. SEASCAPE MAP — all seagrass layers ----

cat("\n  Building seascape map (all seagrass layers)...\n")

p_seascape <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
          linewidth = 0.3) +
  geom_sf(data = scamp_coastline, colour = "grey30", linewidth = 0.3) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  geom_sf(data = solway_sg, aes(fill = "Solway potential (modelled)"),
          colour = NA, alpha = 0.3) +
  geom_sf(data = ea_sg, aes(fill = "EA potential (national)"),
          colour = NA, alpha = 0.3) +
  geom_sf(data = harm_polys[harm_polys$source == "SFP survey", ],
          aes(fill = "Known extent (SFP)"), colour = NA, alpha = 0.6) +
  geom_sf(data = harm_polys[harm_polys$source != "SFP survey", ],
          aes(fill = "Potentially unmapped"), colour = NA, alpha = 0.6) +
  geom_sf(data = sg_points_clean, aes(colour = source, shape = source),
          size = 1.2, alpha = 0.7) +
  scale_fill_manual(
    values = c("Known extent (SFP)"           = "#1b9e77",
               "Potentially unmapped"          = "#e41a1c",
               "EA potential (national)"       = "#66a61e",
               "Solway potential (modelled)"   = "#7570b3"),
    name = "Polygons") +
  scale_colour_manual(
    values = c("GeMS habitat point" = "#d95f02", "WCMC point" = "#a6761d",
               "Indicator records" = "#e6ab02"),
    name = "Points (draft)") +
  scale_shape_manual(
    values = c("GeMS habitat point" = 16, "WCMC point" = 17,
               "Indicator records" = 4),
    name = "Points (draft)") +
  annotation_scale(location = "bl", width_hint = 0.2) +
  annotation_north_arrow(location = "tl", style = north_arrow_minimal(),
                         height = unit(1, "cm"), width = unit(1, "cm")) +
  labs(title = "Seagrass seascape — SCAMP AOI",
       subtitle = sprintf(
         "Known: %d | Unmapped: %d | EA potential: %d | Solway model: %d | Points: %d",
         nrow(harm_polys[harm_polys$source == "SFP survey", ]),
         nrow(harm_polys[harm_polys$source != "SFP survey", ]),
         nrow(ea_sg), nrow(solway_sg), nrow(sg_points_clean))) +
  map_theme +
  coord_sf(crs = 27700)

ggsave("figures/habitats/seagrass_seascape.png", p_seascape,
       width = 14, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/habitats/seagrass_seascape.png\n")


## ---- 1k. CLUSTER MAPS — seascape per cluster ----

cat("\n  Building per-cluster seascape maps...\n")

for (i in seq_len(nrow(zoom_boxes))) {
  box_i <- zoom_boxes[i, ]
  bb_i  <- st_bbox(box_i)

  solway_i  <- solway_sg[lengths(st_intersects(solway_sg, box_i)) > 0, ]
  ea_i      <- ea_sg[lengths(st_intersects(ea_sg, box_i)) > 0, ]
  harm_i    <- harm_polys[lengths(st_intersects(harm_polys, box_i)) > 0, ]
  pts_i     <- sg_points_clean[lengths(st_intersects(sg_points_clean, box_i)) > 0, ]

  sfp_i     <- harm_i[harm_i$source == "SFP survey", ]
  unmapped_i <- harm_i[harm_i$source != "SFP survey", ]

  p_ci <- ggplot() +
    geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
    geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
            linewidth = 0.3) +
    geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
    {if (nrow(solway_i) > 0)
      geom_sf(data = solway_i, aes(fill = "Solway potential (modelled)"),
              colour = NA, alpha = 0.3)} +
    {if (nrow(ea_i) > 0)
      geom_sf(data = ea_i, aes(fill = "EA potential (national)"),
              colour = NA, alpha = 0.3)} +
    {if (nrow(sfp_i) > 0)
      geom_sf(data = sfp_i, aes(fill = "Known extent (SFP)"),
              colour = NA, alpha = 0.6)} +
    {if (nrow(unmapped_i) > 0)
      geom_sf(data = unmapped_i, aes(fill = "Potentially unmapped"),
              colour = NA, alpha = 0.6)} +
    scale_fill_manual(
      values = c("Known extent (SFP)"           = "#1b9e77",
                 "Potentially unmapped"          = "#e41a1c",
                 "EA potential (national)"       = "#66a61e",
                 "Solway potential (modelled)"   = "#7570b3"),
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
    labs(title = sprintf("%s — seagrass seascape", box_i$name),
         subtitle = sprintf(
           "Known: %d | Unmapped: %d | EA: %d | Solway model: %d | Points: %d",
           nrow(sfp_i), nrow(unmapped_i), nrow(ea_i),
           nrow(solway_i), nrow(pts_i))) +
    map_theme

  fname <- paste0("figures/habitats/seagrass_seascape_",
                  gsub("[^a-z0-9]+", "_", tolower(box_i$name)), ".png")
  ggsave(fname, p_ci, width = 12, height = 10, units = "in", dpi = 300)
  cat(sprintf("  Saved: %s\n", fname))
}

cat("\n=== Seagrass seascape maps complete ===\n")


##=====================================================================
## 2. SALTMARSH
##=====================================================================

cat("\n\n=== 2. SALTMARSH ===\n")

dir.create("data/processed/species", showWarnings = FALSE, recursive = TRUE)

## ---- 2a. Read and clip Saltmarsh Survey of Scotland ----

cat("\n-- 2a. Saltmarsh Survey of Scotland --\n")

sss_path <- "data/raw/habitats/Saltmarsh/Saltmarsh_Survey_of_Scotland_460138621685825751.gpkg"
stopifnot("Saltmarsh Survey of Scotland not found" = file.exists(sss_path))

cat("  Layers:\n")
print(st_layers(sss_path))

sss_raw <- st_read(sss_path, quiet = TRUE)
cat(sprintf("  CRS: %s | Rows: %d\n", st_crs(sss_raw)$input, nrow(sss_raw)))
cat(sprintf("  Columns: %s\n", paste(names(sss_raw), collapse = ", ")))

if (st_crs(sss_raw)$epsg != 27700) sss_raw <- st_transform(sss_raw, 27700)
sss_raw <- st_make_valid(sss_raw)
sss_clip <- st_intersection(sss_raw, scamp_aoi_sf)
st_geometry(sss_clip) <- "geometry"
sss_clip$source <- "Saltmarsh Survey of Scotland"
cat(sprintf("  Clipped to AOI: %d features\n", nrow(sss_clip)))

## ---- 2b. Read and clip WCMC Saltmarsh ----

cat("\n-- 2b. WCMC Saltmarsh --\n")

wcmc_sm_path <- "data/raw/habitats/Saltmarsh/WCMC027_Saltmarsh_v6_1/WCMC027_Saltmarsh_v6_1/01_Data/WCMC027_Saltmarshes_Py_v6_1.shp"
stopifnot("WCMC Saltmarsh not found" = file.exists(wcmc_sm_path))

wcmc_sm_raw <- st_read(wcmc_sm_path, quiet = TRUE)
wcmc_sm_crs <- st_crs(wcmc_sm_raw)
cat(sprintf("  CRS: %s | Rows: %d\n", wcmc_sm_crs$input, nrow(wcmc_sm_raw)))
cat(sprintf("  Columns: %s\n", paste(names(wcmc_sm_raw), collapse = ", ")))

aoi_wcmc_sm <- st_transform(scamp_aoi_sf, wcmc_sm_crs) %>% st_make_valid()
wcmc_sm_clip <- st_intersection(st_make_valid(wcmc_sm_raw), aoi_wcmc_sm)
wcmc_sm_clip <- st_transform(wcmc_sm_clip, 27700)
wcmc_sm_clip$source <- "WCMC"
cat(sprintf("  Clipped to AOI: %d features\n", nrow(wcmc_sm_clip)))

## ---- 2c. Cluster analysis — SSS ----

cat("\n-- 2c. Spatial clustering — SSS --\n")

sss_centroids <- st_centroid(sss_clip)
sss_buf <- st_buffer(sss_centroids, dist = 2000) %>% st_union() %>%
  st_cast("POLYGON")
sss_cluster_ids <- st_intersects(sss_centroids, sss_buf)
sss_clip$cluster <- sapply(sss_cluster_ids, `[`, 1)

sss_cluster_stats <- sss_clip %>%
  st_drop_geometry() %>%
  group_by(cluster) %>%
  summarise(n = n(), .groups = "drop")

cat(sprintf("  SSS clusters: %d\n", nrow(sss_cluster_stats)))
print(sss_cluster_stats)

## Build cluster bboxes for SSS
sss_cluster_boxes <- do.call(rbind, lapply(sort(unique(sss_clip$cluster)), function(cl) {
  bb <- st_bbox(sss_clip[sss_clip$cluster == cl, ])
  st_sf(
    cluster = cl,
    n = sum(sss_clip$cluster == cl),
    name = paste0("SSS_", cl),
    geometry = st_as_sfc(bb, crs = 27700)
  )
}))

## ---- 2d. Cluster analysis — WCMC ----

cat("\n-- 2d. Spatial clustering — WCMC --\n")

wcmc_sm_centroids <- st_centroid(wcmc_sm_clip)
wcmc_sm_buf <- st_buffer(wcmc_sm_centroids, dist = 2000) %>% st_union() %>%
  st_cast("POLYGON")
wcmc_sm_cluster_ids <- st_intersects(wcmc_sm_centroids, wcmc_sm_buf)
wcmc_sm_clip$cluster <- sapply(wcmc_sm_cluster_ids, `[`, 1)

wcmc_sm_cluster_stats <- wcmc_sm_clip %>%
  st_drop_geometry() %>%
  group_by(cluster) %>%
  summarise(n = n(), .groups = "drop")

cat(sprintf("  WCMC clusters: %d\n", nrow(wcmc_sm_cluster_stats)))
print(wcmc_sm_cluster_stats)

## Build cluster bboxes for WCMC
wcmc_sm_cluster_boxes <- do.call(rbind, lapply(sort(unique(wcmc_sm_clip$cluster)), function(cl) {
  bb <- st_bbox(wcmc_sm_clip[wcmc_sm_clip$cluster == cl, ])
  st_sf(
    cluster = cl,
    n = sum(wcmc_sm_clip$cluster == cl),
    name = paste0("WCMC_", cl),
    geometry = st_as_sfc(bb, crs = 27700)
  )
}))

## Combined cluster boxes (union of both sources)
sm_all_polys <- rbind(
  sss_clip %>% select(source, geometry),
  wcmc_sm_clip %>% select(source, geometry)
)


sm_centroids_all <- st_centroid(sm_all_polys)
sm_buf_all <- st_buffer(sm_centroids_all, dist = 2000) %>% st_union() %>%
  st_cast("POLYGON")
sm_cluster_ids_all <- st_intersects(sm_centroids_all, sm_buf_all)
sm_all_polys$cluster <- sapply(sm_cluster_ids_all, `[`, 1)

sm_cluster_boxes <- do.call(rbind, lapply(sort(unique(sm_all_polys$cluster)), function(cl) {
  bb <- st_bbox(sm_all_polys[sm_all_polys$cluster == cl, ])
  st_sf(
    cluster = cl,
    n = sum(sm_all_polys$cluster == cl),
    name = paste0("Saltmarsh_", cl),
    geometry = st_as_sfc(bb, crs = 27700)
  )
}))

cat(sprintf("\n  Combined saltmarsh clusters: %d\n", nrow(sm_cluster_boxes)))

## ---- 2e. Per-cluster plots — SSS vs WCMC comparison ----

cat("\n-- 2e. Per-cluster comparison plots --\n")

for (i in seq_len(nrow(sm_cluster_boxes))) {
  box_i <- sm_cluster_boxes[i, ]

  sss_i  <- sss_clip[lengths(st_intersects(sss_clip, box_i)) > 0, ]
  wcmc_i <- wcmc_sm_clip[lengths(st_intersects(wcmc_sm_clip, box_i)) > 0, ]

  p_i <- ggplot() +
    geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
    geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
            linewidth = 0.3) +
    {if (nrow(wcmc_i) > 0)
      geom_sf(data = wcmc_i, aes(fill = "WCMC"), colour = NA, alpha = 0.5)} +
    {if (nrow(sss_i) > 0)
      geom_sf(data = sss_i, aes(fill = "SSS"), colour = NA, alpha = 0.5)} +
    geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
    scale_fill_manual(name = "Source",
                      values = c("SSS" = "darkorange", "WCMC" = "steelblue")) +
    annotation_scale(location = "bl", width_hint = 0.2) +
    annotation_north_arrow(location = "tl", style = north_arrow_minimal(),
                           height = unit(1, "cm"), width = unit(1, "cm")) +
    labs(title = sprintf("Saltmarsh Cluster %d (%s)", box_i$cluster, box_i$name),
         subtitle = sprintf("SSS: %d | WCMC: %d", nrow(sss_i), nrow(wcmc_i))) +
    map_theme +
    coord_sf(crs = 27700,
             xlim = c(st_bbox(box_i)["xmin"], st_bbox(box_i)["xmax"]),
             ylim = c(st_bbox(box_i)["ymin"], st_bbox(box_i)["ymax"]))

  fname <- sprintf("figures/habitats/saltmarsh_cluster_%02d.png", box_i$cluster)
  ggsave(fname, p_i, width = 12, height = 10, units = "in", dpi = 300)
  cat(sprintf("  Saved: %s\n", fname))
}

## ---- 2f. Overview plot — both sources ----

cat("\n-- 2f. Overview plot --\n")

p_sm_overview <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
          linewidth = 0.3) +
  {if (nrow(wcmc_sm_clip) > 0)
    geom_sf(data = wcmc_sm_clip, aes(fill = "WCMC"), colour = NA, alpha = 0.5)} +
  {if (nrow(sss_clip) > 0)
    geom_sf(data = sss_clip, aes(fill = "SSS"), colour = NA, alpha = 0.5)} +
  geom_sf(data = sm_cluster_boxes, fill = NA, colour = "red", linewidth = 0.6) +
  geom_sf_text(data = sm_cluster_boxes, aes(label = name), colour = "red", size = 2.5) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  scale_fill_manual(name = "Source",
                    values = c("SSS" = "darkorange", "WCMC" = "steelblue")) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  annotation_north_arrow(location = "tl", style = north_arrow_minimal(),
                         height = unit(1, "cm"), width = unit(1, "cm")) +
  labs(title = "Saltmarsh — SSS & WCMC (SCAMP AOI)",
       subtitle = sprintf("SSS: %d | WCMC: %d | Clusters: %d",
                          nrow(sss_clip), nrow(wcmc_sm_clip), nrow(sm_cluster_boxes))) +
  map_theme + coord_sf(crs = 27700)

ggsave("figures/habitats/saltmarsh_overview.png", p_sm_overview,
       width = 14, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/habitats/saltmarsh_overview.png\n")

## ---- 2g. Union into harmonised saltmarsh polygon layer ----

cat("\n-- 2g. Harmonised saltmarsh layer --\n")

sm_combined <- rbind(
  sss_clip %>% mutate(source = "SSS") %>% select(source, geometry),
  wcmc_sm_clip %>% mutate(source = "WCMC") %>% select(source, geometry)
)

## Where polygons overlap, merge them; retain non-overlapping polygons as-is
sm_unioned <- st_union(sm_combined) %>% st_make_valid()
sm_harmonised_sf <- st_sf(
  habitat = "Saltmarsh",
  geometry = st_cast(sm_unioned, "POLYGON")
) %>% st_make_valid()

cat(sprintf("  Input features: %d (SSS: %d, WCMC: %d)\n",
            nrow(sm_combined), sum(sm_combined$source == "SSS"),
            sum(sm_combined$source == "WCMC")))
cat(sprintf("  Harmonised polygons: %d (overlapping merged, non-overlapping retained)\n",
            nrow(sm_harmonised_sf)))

st_write(sm_harmonised_sf, "data/processed/species/saltmarsh_harmonised.gpkg",
         delete_dsn = TRUE, quiet = TRUE)
cat("  Saved: data/processed/species/saltmarsh_harmonised.gpkg\n")
cat(sprintf("  Harmonised area: %.1f ha\n",
            sum(as.numeric(st_area(sm_harmonised_sf))) / 10000))

## Per-cluster harmonised maps — showing merged result
for (i in seq_len(nrow(sm_cluster_boxes))) {
  box_i <- sm_cluster_boxes[i, ]
  harm_i <- sm_harmonised_sf[lengths(st_intersects(sm_harmonised_sf, box_i)) > 0, ]

  p_hi <- ggplot() +
    geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
    geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
            linewidth = 0.3) +
    {if (nrow(harm_i) > 0)
      geom_sf(data = harm_i, fill = "darkgreen", colour = "darkgreen",
              alpha = 0.5, linewidth = 0.3)} +
    geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
    annotation_scale(location = "bl", width_hint = 0.2) +
    annotation_north_arrow(location = "tl", style = north_arrow_minimal(),
                           height = unit(1, "cm"), width = unit(1, "cm")) +
    labs(title = sprintf("Saltmarsh Cluster %d — Harmonised (%s)", box_i$cluster, box_i$name),
         subtitle = sprintf("Merged polygons: %d | Area: %.1f ha",
                            nrow(harm_i),
                            sum(as.numeric(st_area(harm_i))) / 10000)) +
    map_theme +
    coord_sf(crs = 27700,
             xlim = c(st_bbox(box_i)["xmin"], st_bbox(box_i)["xmax"]),
             ylim = c(st_bbox(box_i)["ymin"], st_bbox(box_i)["ymax"]))

  fname <- sprintf("figures/habitats/saltmarsh_cluster_%02d_harmonised.png", box_i$cluster)
  ggsave(fname, p_hi, width = 12, height = 10, units = "in", dpi = 300)
  cat(sprintf("  Saved: %s\n", fname))
}

## ---- 2h. Saltmarsh indicator species ----

cat("\n-- 2h. Saltmarsh indicator species --\n")

ind_path <- "outputs/species/indicator_records.gpkg"
stopifnot("Indicator records not found" = file.exists(ind_path))

ind_all <- st_read(ind_path, quiet = TRUE)
cat(sprintf("  Total indicator records: %d\n", nrow(ind_all)))
cat(sprintf("  Columns: %s\n", paste(names(ind_all), collapse = ", ")))

ind_sm <- ind_all %>%
  filter(habitat == "Saltmarsh") %>%
  st_transform(27700) %>%
  st_make_valid() %>%
  st_intersection(scamp_aoi_sf)

cat(sprintf("  Saltmarsh indicator records in AOI: %d\n", nrow(ind_sm)))
if (nrow(ind_sm) > 0) {
  cat("  Species breakdown:\n")
  print(ind_sm %>% st_drop_geometry() %>% count(latin_name, sort = TRUE))
}

## ---- 2i. Plot saltmarsh indicators with harmonised layer ----

cat("\n-- 2i. Saltmarsh indicator map --\n")

p_sm_ind <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
          linewidth = 0.3) +
  geom_sf(data = sm_harmonised_sf, fill = "darkgreen", colour = NA, alpha = 0.4) +
  {if (nrow(ind_sm) > 0)
    geom_sf(data = ind_sm, aes(colour = latin_name), size = 1.2, alpha = 0.7)} +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  annotation_north_arrow(location = "tl", style = north_arrow_minimal(),
                         height = unit(1, "cm"), width = unit(1, "cm")) +
  labs(title = "Saltmarsh — Harmonised Extent & Indicator Species (SCAMP AOI)",
       subtitle = sprintf("Harmonised polygons: %.1f ha | Indicator records: %d",
                          sum(as.numeric(st_area(sm_harmonised_sf))) / 10000,
                          nrow(ind_sm)),
       colour = "Indicator species") +
  map_theme + coord_sf(crs = 27700) +
  guides(colour = guide_legend(override.aes = list(size = 3)))

ggsave("figures/habitats/saltmarsh_indicators.png", p_sm_ind,
       width = 14, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/habitats/saltmarsh_indicators.png\n")

## Per-cluster indicator plots
for (i in seq_len(nrow(sm_cluster_boxes))) {
  box_i <- sm_cluster_boxes[i, ]
  harm_i <- sm_harmonised_sf[lengths(st_intersects(sm_harmonised_sf, box_i)) > 0, ]
  ind_i <- ind_sm[lengths(st_intersects(ind_sm, box_i)) > 0, ]

  p_ci <- ggplot() +
    geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
    geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
            linewidth = 0.3) +
    {if (nrow(harm_i) > 0)
      geom_sf(data = harm_i, fill = "darkgreen", colour = NA, alpha = 0.4)} +
    {if (nrow(ind_i) > 0)
      geom_sf(data = ind_i, aes(colour = latin_name), size = 1.5, alpha = 0.7)} +
    geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
    annotation_scale(location = "bl", width_hint = 0.2) +
    annotation_north_arrow(location = "tl", style = north_arrow_minimal(),
                           height = unit(1, "cm"), width = unit(1, "cm")) +
    labs(title = sprintf("Saltmarsh Cluster %d — Indicators", box_i$cluster),
         subtitle = sprintf("Polygons: %.1f ha | Indicator pts: %d",
                            sum(as.numeric(st_area(harm_i))) / 10000, nrow(ind_i)),
         colour = "Species") +
    map_theme +
    coord_sf(crs = 27700,
             xlim = c(st_bbox(box_i)["xmin"], st_bbox(box_i)["xmax"]),
             ylim = c(st_bbox(box_i)["ymin"], st_bbox(box_i)["ymax"])) +
    guides(colour = guide_legend(override.aes = list(size = 3)))

  fname <- sprintf("figures/habitats/saltmarsh_cluster_%02d_indicators.png", box_i$cluster)
  ggsave(fname, p_ci, width = 12, height = 10, units = "in", dpi = 300)
  cat(sprintf("  Saved: %s\n", fname))
}

## Save indicator points for reference
if (nrow(ind_sm) > 0) {
  st_write(ind_sm, "data/processed/species/saltmarsh_indicator_points.gpkg",
           delete_dsn = TRUE, quiet = TRUE)
  cat("  Saved: data/processed/species/saltmarsh_indicator_points.gpkg\n")
}

cat("\n=== Saltmarsh section complete ===\n")


##=====================================================================
## 3. OYSTERS
##=====================================================================

cat("\n\n=== 3. OYSTERS ===\n")

gems_path <- "data/raw/GEMS-PMF_SCOTLAND_GPKG_4326/gems.gpkg"

## ---- 3a. GeMS native oyster — polygons & points ----

cat("\n-- 3a. GeMS native oyster beds --\n")

gems_oy_py <- st_read(gems_path, layer = "gems_habitat_polygon_dataset",
                       quiet = TRUE) %>%
  filter(PMF == "Native oysters") %>%
  st_transform(27700) %>%
  st_make_valid() %>%
  st_intersection(scamp_aoi_sf) %>%
  mutate(source = "GeMS habitat polygon")

gems_oy_pt <- st_read(gems_path, layer = "gems_habitat_point_dataset",
                       quiet = TRUE) %>%
  filter(PMF == "Native oysters") %>%
  st_transform(27700) %>%
  st_make_valid() %>%
  st_intersection(scamp_aoi_sf) %>%
  mutate(source = "GeMS habitat point")

cat(sprintf("  GeMS oyster polygons in AOI: %d\n", nrow(gems_oy_py)))
cat(sprintf("  GeMS oyster points in AOI:   %d\n", nrow(gems_oy_pt)))

## Plot GeMS oyster data
p_gems_oy <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
          linewidth = 0.3) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  {if (nrow(gems_oy_py) > 0)
    geom_sf(data = gems_oy_py, fill = "#e7298a", colour = "#e7298a", alpha = 0.5)} +
  {if (nrow(gems_oy_pt) > 0)
    geom_sf(data = gems_oy_pt, colour = "#d95f02", shape = 16, size = 2, alpha = 0.7)} +
  labs(title = "GeMS Native Oyster Beds (SCAMP AOI)",
       subtitle = sprintf("Polygons: %d | Points: %d",
                          nrow(gems_oy_py), nrow(gems_oy_pt))) +
  map_theme + coord_sf(crs = 27700)

ggsave("figures/habitats/oyster_gems.png", p_gems_oy,
       width = 14, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/habitats/oyster_gems.png\n")

## Cluster GeMS polygons (2km buffer on centroids)
if (nrow(gems_oy_py) > 1) {
  gems_oy_centroids <- st_centroid(gems_oy_py)
  gems_oy_buf <- st_buffer(gems_oy_centroids, dist = 2000) %>% st_union() %>%
    st_cast("POLYGON")
  gems_oy_cids <- st_intersects(gems_oy_centroids, gems_oy_buf)
  gems_oy_py$cluster <- sapply(gems_oy_cids, `[`, 1)
} else if (nrow(gems_oy_py) == 1) {
  gems_oy_py$cluster <- 1L
}

## Build cluster bboxes from GeMS polygons
if (nrow(gems_oy_py) > 0) {
  oy_boxes <- do.call(rbind, lapply(sort(unique(gems_oy_py$cluster)), function(cl) {
    bb <- st_bbox(st_buffer(gems_oy_py[gems_oy_py$cluster == cl, ], dist = 2000))
    st_sf(
      name = sprintf("GeMS_%d", cl),
      geometry = st_sfc(st_polygon(list(matrix(c(
        bb["xmin"], bb["ymin"], bb["xmax"], bb["ymin"],
        bb["xmax"], bb["ymax"], bb["xmin"], bb["ymax"],
        bb["xmin"], bb["ymin"]
      ), ncol = 2, byrow = TRUE))), crs = 27700)
    )
  }))
  cat(sprintf("  GeMS oyster clusters: %d\n", nrow(oy_boxes)))
} else {
  oy_boxes <- st_sf(name = character(0),
                    geometry = st_sfc(crs = 27700))
  cat("  No GeMS oyster polygons — starting with empty cluster set\n")
}

## Plot cluster bboxes
p_oy_boxes_gems <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
          linewidth = 0.3) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  {if (nrow(gems_oy_py) > 0)
    geom_sf(data = gems_oy_py, fill = "#e7298a", colour = NA, alpha = 0.5)} +
  {if (nrow(gems_oy_pt) > 0)
    geom_sf(data = gems_oy_pt, colour = "#d95f02", shape = 16, size = 2, alpha = 0.7)} +
  {if (nrow(oy_boxes) > 0)
    geom_sf(data = oy_boxes, fill = NA, colour = "red", linewidth = 0.8)} +
  {if (nrow(oy_boxes) > 0)
    geom_sf_text(data = oy_boxes, aes(label = name), colour = "red",
                 size = 3, fontface = "bold")} +
  labs(title = "Oyster cluster bounding boxes — GeMS",
       subtitle = sprintf("%d clusters defined", nrow(oy_boxes))) +
  map_theme + coord_sf(crs = 27700)

ggsave("figures/habitats/oyster_gems_clusters.png", p_oy_boxes_gems,
       width = 14, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/habitats/oyster_gems_clusters.png\n")


## ---- 3b. Historical oyster locations ----

cat("\n-- 3b. Historical oyster locations --\n")

hist_oy_path <- "data/raw/habitats/Oyster/historical_oyster_distribution/01_Oyster_Habitat_Locations.shp"
stopifnot("Historical oyster shapefile not found" = file.exists(hist_oy_path))

hist_oy <- st_read(hist_oy_path, quiet = TRUE) %>%
  st_transform(27700) %>%
  st_make_valid() %>%
  st_intersection(scamp_aoi_sf) %>%
  mutate(source = "Historical")

cat(sprintf("  Historical oyster locations in AOI: %d\n", nrow(hist_oy)))
if (nrow(hist_oy) > 0) {
  cat(sprintf("  Geometry types: %s\n",
              paste(unique(st_geometry_type(hist_oy)), collapse = ", ")))
}

## Separate polygons and points from historical data
hist_oy_geom_types <- st_geometry_type(hist_oy)
hist_oy_py <- hist_oy[hist_oy_geom_types %in% c("POLYGON", "MULTIPOLYGON"), ]
hist_oy_pt <- hist_oy[hist_oy_geom_types %in% c("POINT", "MULTIPOINT"), ]
cat(sprintf("  Historical polygons: %d | points: %d\n",
            nrow(hist_oy_py), nrow(hist_oy_pt)))

## Plot historical overlaid on GeMS
p_hist_oy <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
          linewidth = 0.3) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  {if (nrow(gems_oy_py) > 0)
    geom_sf(data = gems_oy_py, aes(fill = "GeMS polygon"), colour = NA, alpha = 0.5)} +
  {if (nrow(hist_oy_py) > 0)
    geom_sf(data = hist_oy_py, aes(fill = "Historical polygon"), colour = NA, alpha = 0.4)} +
  {if (nrow(hist_oy_pt) > 0)
    geom_sf(data = hist_oy_pt, aes(colour = "Historical point"), shape = 17,
            size = 2, alpha = 0.7)} +
  {if (nrow(gems_oy_pt) > 0)
    geom_sf(data = gems_oy_pt, aes(colour = "GeMS point"), shape = 16,
            size = 2, alpha = 0.7)} +
  scale_fill_manual(
    values = c("GeMS polygon" = "#e7298a", "Historical polygon" = "#a6761d"),
    name = "Polygons") +
  scale_colour_manual(
    values = c("GeMS point" = "#d95f02", "Historical point" = "#7570b3"),
    name = "Points") +
  labs(title = "Native Oyster — GeMS + Historical (SCAMP AOI)",
       subtitle = sprintf("GeMS py: %d | GeMS pt: %d | Hist py: %d | Hist pt: %d",
                          nrow(gems_oy_py), nrow(gems_oy_pt),
                          nrow(hist_oy_py), nrow(hist_oy_pt))) +
  map_theme + coord_sf(crs = 27700)

ggsave("figures/habitats/oyster_gems_historical.png", p_hist_oy,
       width = 14, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/habitats/oyster_gems_historical.png\n")

## Cluster historical locations
hist_oy_all_pts <- bind_rows(
  if (nrow(hist_oy_pt) > 0) hist_oy_pt %>% select(source, geometry),
  if (nrow(hist_oy_py) > 0) st_centroid(hist_oy_py) %>% select(source, geometry)
)

if (nrow(hist_oy_all_pts) > 1) {
  hist_buf <- st_buffer(hist_oy_all_pts, dist = 2000) %>% st_union() %>%
    st_cast("POLYGON")
  hist_cids <- st_intersects(hist_oy_all_pts, hist_buf)
  hist_oy_all_pts$cluster <- sapply(hist_cids, `[`, 1)

  hist_oy_boxes <- do.call(rbind, lapply(sort(unique(hist_oy_all_pts$cluster)), function(cl) {
    bb <- st_bbox(st_buffer(hist_oy_all_pts[hist_oy_all_pts$cluster == cl, ], dist = 2000))
    st_sf(
      name = sprintf("Hist_%d", cl),
      geometry = st_sfc(st_polygon(list(matrix(c(
        bb["xmin"], bb["ymin"], bb["xmax"], bb["ymin"],
        bb["xmax"], bb["ymax"], bb["xmin"], bb["ymax"],
        bb["xmin"], bb["ymin"]
      ), ncol = 2, byrow = TRUE))), crs = 27700)
    )
  }))
} else if (nrow(hist_oy_all_pts) == 1) {
  hist_oy_all_pts$cluster <- 1L
  bb <- st_bbox(st_buffer(hist_oy_all_pts, dist = 2000))
  hist_oy_boxes <- st_sf(
    name = "Hist_1",
    geometry = st_sfc(st_polygon(list(matrix(c(
      bb["xmin"], bb["ymin"], bb["xmax"], bb["ymin"],
      bb["xmax"], bb["ymax"], bb["xmin"], bb["ymax"],
      bb["xmin"], bb["ymin"]
    ), ncol = 2, byrow = TRUE))), crs = 27700)
  )
} else {
  hist_oy_boxes <- st_sf(name = character(0), geometry = st_sfc(crs = 27700))
}

cat(sprintf("  Historical oyster clusters: %d\n", nrow(hist_oy_boxes)))

## Merge historical boxes with GeMS boxes — combine overlapping, rbind otherwise
if (nrow(hist_oy_boxes) > 0 && nrow(oy_boxes) > 0) {
  combined_boxes <- rbind(oy_boxes, hist_oy_boxes)
  cb_buf <- st_buffer(combined_boxes, dist = 1000) %>% st_union() %>%
    st_cast("POLYGON")
  cb_group <- sapply(st_intersects(combined_boxes, cb_buf), `[`, 1)

  oy_boxes <- do.call(rbind, lapply(sort(unique(cb_group)), function(g) {
    bb <- st_bbox(combined_boxes[cb_group == g, ])
    st_sf(
      name = sprintf("Cluster %d", g),
      geometry = st_sfc(st_polygon(list(matrix(c(
        bb["xmin"], bb["ymin"], bb["xmax"], bb["ymin"],
        bb["xmax"], bb["ymax"], bb["xmin"], bb["ymax"],
        bb["xmin"], bb["ymin"]
      ), ncol = 2, byrow = TRUE))), crs = 27700)
    )
  }))
} else if (nrow(hist_oy_boxes) > 0) {
  oy_boxes <- hist_oy_boxes
}

cat(sprintf("  Merged oyster clusters (GeMS + Historical): %d\n", nrow(oy_boxes)))

## Plot merged clusters
p_oy_boxes_hist <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
          linewidth = 0.3) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  {if (nrow(gems_oy_py) > 0)
    geom_sf(data = gems_oy_py, aes(fill = "GeMS polygon"), colour = NA, alpha = 0.5)} +
  {if (nrow(hist_oy_py) > 0)
    geom_sf(data = hist_oy_py, aes(fill = "Historical polygon"), colour = NA, alpha = 0.4)} +
  scale_fill_manual(
    values = c("GeMS polygon" = "#e7298a", "Historical polygon" = "#a6761d"),
    name = "Polygons") +
  {if (nrow(oy_boxes) > 0)
    geom_sf(data = oy_boxes, fill = NA, colour = "red", linewidth = 0.8)} +
  {if (nrow(oy_boxes) > 0)
    geom_sf_text(data = oy_boxes, aes(label = name), colour = "red",
                 size = 3, fontface = "bold")} +
  labs(title = "Oyster cluster bounding boxes — GeMS + Historical",
       subtitle = sprintf("%d clusters", nrow(oy_boxes))) +
  map_theme + coord_sf(crs = 27700)

ggsave("figures/habitats/oyster_gems_hist_clusters.png", p_oy_boxes_hist,
       width = 14, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/habitats/oyster_gems_hist_clusters.png\n")

## Per-cluster maps (GeMS + Historical)
for (i in seq_len(nrow(oy_boxes))) {
  box_i <- oy_boxes[i, ]
  bb_i  <- st_bbox(box_i)

  gpy_i  <- gems_oy_py[lengths(st_intersects(gems_oy_py, box_i)) > 0, ]
  gpt_i  <- gems_oy_pt[lengths(st_intersects(gems_oy_pt, box_i)) > 0, ]
  hpy_i  <- hist_oy_py[lengths(st_intersects(hist_oy_py, box_i)) > 0, ]
  hpt_i  <- hist_oy_pt[lengths(st_intersects(hist_oy_pt, box_i)) > 0, ]

  p_i <- ggplot() +
    geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
    geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
            linewidth = 0.3) +
    geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
    {if (nrow(gpy_i) > 0) geom_sf(data = gpy_i, aes(fill = "GeMS polygon"),
                                   colour = NA, alpha = 0.5)} +
    {if (nrow(hpy_i) > 0) geom_sf(data = hpy_i, aes(fill = "Historical polygon"),
                                   colour = NA, alpha = 0.4)} +
    scale_fill_manual(
      values = c("GeMS polygon" = "#e7298a", "Historical polygon" = "#a6761d"),
      name = "Polygons") +
    {if (nrow(gpt_i) > 0) geom_sf(data = gpt_i, aes(colour = "GeMS point"),
                                   shape = 16, size = 2.5, alpha = 0.7)} +
    {if (nrow(hpt_i) > 0) geom_sf(data = hpt_i, aes(colour = "Historical point"),
                                   shape = 17, size = 2.5, alpha = 0.7)} +
    scale_colour_manual(
      values = c("GeMS point" = "#d95f02", "Historical point" = "#7570b3"),
      name = "Points") +
    coord_sf(xlim = c(bb_i["xmin"], bb_i["xmax"]),
             ylim = c(bb_i["ymin"], bb_i["ymax"]), crs = 27700) +
    labs(title = sprintf("%s — oyster (GeMS + Historical)", box_i$name),
         subtitle = sprintf("GeMS py: %d pt: %d | Hist py: %d pt: %d",
                            nrow(gpy_i), nrow(gpt_i), nrow(hpy_i), nrow(hpt_i))) +
    map_theme

  fname <- paste0("figures/habitats/oyster_cluster_gems_hist_",
                  gsub("[^a-z0-9]+", "_", tolower(box_i$name)), ".png")
  ggsave(fname, p_i, width = 12, height = 10, units = "in", dpi = 300)
  cat(sprintf("  Saved: %s\n", fname))
}


## ---- 3c. Indicator species records (oyster) ----

cat("\n-- 3c. Indicator species records (oyster) --\n")

ind_path <- "outputs/species/indicator_records.gpkg"
stopifnot("Indicator records not found" = file.exists(ind_path))

ind_all <- st_read(ind_path, quiet = TRUE)
ind_oy <- ind_all %>%
  filter(habitat == "Oyster") %>%
  filter(grepl('ostrea', latin_name, ignore.case=T))%>%
  st_transform(27700) %>%
  st_make_valid() %>%
  st_intersection(scamp_aoi_sf) %>%
  mutate(source = "Indicator records")

cat(sprintf("  Oyster indicator records in AOI: %d\n", nrow(ind_oy)))
if (nrow(ind_oy) > 0) {
  cat("  Species breakdown:\n")
  print(ind_oy %>% st_drop_geometry() %>% count(latin_name, sort = TRUE))
}

## Plot indicators overlaid on previous layers
p_oy_ind <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
          linewidth = 0.3) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  {if (nrow(gems_oy_py) > 0)
    geom_sf(data = gems_oy_py, aes(fill = "GeMS polygon"), colour = NA, alpha = 0.5)} +
  {if (nrow(hist_oy_py) > 0)
    geom_sf(data = hist_oy_py, aes(fill = "Historical polygon"), colour = NA, alpha = 0.4)} +
  scale_fill_manual(
    values = c("GeMS polygon" = "#e7298a", "Historical polygon" = "#a6761d"),
    name = "Polygons") +
  {if (nrow(gems_oy_pt) > 0)
    geom_sf(data = gems_oy_pt, aes(colour = "GeMS point"), shape = 16,
            size = 2, alpha = 0.7)} +
  {if (nrow(hist_oy_pt) > 0)
    geom_sf(data = hist_oy_pt, aes(colour = "Historical point"), shape = 17,
            size = 2, alpha = 0.7)} +
  {if (nrow(ind_oy) > 0)
    geom_sf(data = ind_oy, aes(colour = "Indicator record"), shape = 4,
            size = 2, alpha = 0.8, stroke = 1)} +
  scale_colour_manual(
    values = c("GeMS point" = "#d95f02", "Historical point" = "#7570b3",
               "Indicator record" = "#e6ab02"),
    name = "Points") +
  labs(title = "Native Oyster — All sources (SCAMP AOI)",
       subtitle = sprintf("GeMS py: %d pt: %d | Hist py: %d pt: %d | Indicator: %d",
                          nrow(gems_oy_py), nrow(gems_oy_pt),
                          nrow(hist_oy_py), nrow(hist_oy_pt), nrow(ind_oy))) +
  map_theme + coord_sf(crs = 27700)

ggsave("figures/habitats/oyster_all_sources.png", p_oy_ind,
       width = 14, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/habitats/oyster_all_sources.png\n")

## Cluster indicator records
if (nrow(ind_oy) > 1) {
  ind_oy_buf <- st_buffer(ind_oy, dist = 2000) %>% st_union() %>%
    st_cast("POLYGON")
  ind_oy_cids <- st_intersects(ind_oy, ind_oy_buf)
  ind_oy$cluster <- sapply(ind_oy_cids, `[`, 1)

  ind_oy_boxes <- do.call(rbind, lapply(sort(unique(ind_oy$cluster)), function(cl) {
    bb <- st_bbox(st_buffer(ind_oy[ind_oy$cluster == cl, ], dist = 2000))
    st_sf(
      name = sprintf("Ind_%d", cl),
      geometry = st_sfc(st_polygon(list(matrix(c(
        bb["xmin"], bb["ymin"], bb["xmax"], bb["ymin"],
        bb["xmax"], bb["ymax"], bb["xmin"], bb["ymax"],
        bb["xmin"], bb["ymin"]
      ), ncol = 2, byrow = TRUE))), crs = 27700)
    )
  }))
} else if (nrow(ind_oy) == 1) {
  ind_oy$cluster <- 1L
  bb <- st_bbox(st_buffer(ind_oy, dist = 2000))
  ind_oy_boxes <- st_sf(
    name = "Ind_1",
    geometry = st_sfc(st_polygon(list(matrix(c(
      bb["xmin"], bb["ymin"], bb["xmax"], bb["ymin"],
      bb["xmax"], bb["ymax"], bb["xmin"], bb["ymax"],
      bb["xmin"], bb["ymin"]
    ), ncol = 2, byrow = TRUE))), crs = 27700)
  )
} else {
  ind_oy_boxes <- st_sf(name = character(0), geometry = st_sfc(crs = 27700))
}

cat(sprintf("  Indicator oyster clusters: %d\n", nrow(ind_oy_boxes)))

## Merge indicator boxes with existing oyster boxes
if (nrow(ind_oy_boxes) > 0 && nrow(oy_boxes) > 0) {
  combined_boxes <- rbind(oy_boxes, ind_oy_boxes)
  cb_buf <- st_buffer(combined_boxes, dist = 1000) %>% st_union() %>%
    st_cast("POLYGON")
  cb_group <- sapply(st_intersects(combined_boxes, cb_buf), `[`, 1)

  oy_boxes <- do.call(rbind, lapply(sort(unique(cb_group)), function(g) {
    bb <- st_bbox(combined_boxes[cb_group == g, ])
    st_sf(
      name = sprintf("Cluster %d", g),
      geometry = st_sfc(st_polygon(list(matrix(c(
        bb["xmin"], bb["ymin"], bb["xmax"], bb["ymin"],
        bb["xmax"], bb["ymax"], bb["xmin"], bb["ymax"],
        bb["xmin"], bb["ymin"]
      ), ncol = 2, byrow = TRUE))), crs = 27700)
    )
  }))
} else if (nrow(ind_oy_boxes) > 0) {
  oy_boxes <- ind_oy_boxes
}

cat(sprintf("  Merged oyster clusters (all sources): %d\n", nrow(oy_boxes)))

## Plot merged clusters with all data
p_oy_boxes_all <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
          linewidth = 0.3) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  {if (nrow(gems_oy_py) > 0)
    geom_sf(data = gems_oy_py, aes(fill = "GeMS polygon"), colour = NA, alpha = 0.5)} +
  {if (nrow(hist_oy_py) > 0)
    geom_sf(data = hist_oy_py, aes(fill = "Historical polygon"), colour = NA, alpha = 0.4)} +
  scale_fill_manual(
    values = c("GeMS polygon" = "#e7298a", "Historical polygon" = "#a6761d"),
    name = "Polygons") +
  {if (nrow(gems_oy_pt) > 0)
    geom_sf(data = gems_oy_pt, aes(colour = "GeMS point"), shape = 16,
            size = 1.5, alpha = 0.7)} +
  {if (nrow(hist_oy_pt) > 0)
    geom_sf(data = hist_oy_pt, aes(colour = "Historical point"), shape = 17,
            size = 1.5, alpha = 0.7)} +
  {if (nrow(ind_oy) > 0)
    geom_sf(data = ind_oy, aes(colour = "Indicator record"), shape = 4,
            size = 2, alpha = 0.8, stroke = 1)} +
  scale_colour_manual(
    values = c("GeMS point" = "#d95f02", "Historical point" = "#7570b3",
               "Indicator record" = "#e6ab02"),
    name = "Points") +
  {if (nrow(oy_boxes) > 0)
    geom_sf(data = oy_boxes, fill = NA, colour = "red", linewidth = 0.8)} +
  {if (nrow(oy_boxes) > 0)
    geom_sf_text(data = oy_boxes, aes(label = name), colour = "red",
                 size = 3, fontface = "bold")} +
  labs(title = "Oyster cluster bounding boxes — all sources",
       subtitle = sprintf("%d clusters", nrow(oy_boxes))) +
  map_theme + coord_sf(crs = 27700)

ggsave("figures/habitats/oyster_all_clusters.png", p_oy_boxes_all,
       width = 14, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/habitats/oyster_all_clusters.png\n")

## Per-cluster maps (all sources)
for (i in seq_len(nrow(oy_boxes))) {
  box_i <- oy_boxes[i, ]
  bb_i  <- st_bbox(box_i)

  gpy_i  <- gems_oy_py[lengths(st_intersects(gems_oy_py, box_i)) > 0, ]
  gpt_i  <- gems_oy_pt[lengths(st_intersects(gems_oy_pt, box_i)) > 0, ]
  hpy_i  <- hist_oy_py[lengths(st_intersects(hist_oy_py, box_i)) > 0, ]
  hpt_i  <- hist_oy_pt[lengths(st_intersects(hist_oy_pt, box_i)) > 0, ]
  ind_i  <- ind_oy[lengths(st_intersects(ind_oy, box_i)) > 0, ]

  p_i <- ggplot() +
    geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
    geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
            linewidth = 0.3) +
    geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
    {if (nrow(gpy_i) > 0) geom_sf(data = gpy_i, aes(fill = "GeMS polygon"),
                                   colour = NA, alpha = 0.5)} +
    {if (nrow(hpy_i) > 0) geom_sf(data = hpy_i, aes(fill = "Historical polygon"),
                                   colour = NA, alpha = 0.4)} +
    scale_fill_manual(
      values = c("GeMS polygon" = "#e7298a", "Historical polygon" = "#a6761d"),
      name = "Polygons") +
    {if (nrow(gpt_i) > 0) geom_sf(data = gpt_i, aes(colour = "GeMS point"),
                                   shape = 16, size = 2.5, alpha = 0.7)} +
    {if (nrow(hpt_i) > 0) geom_sf(data = hpt_i, aes(colour = "Historical point"),
                                   shape = 17, size = 2.5, alpha = 0.7)} +
    {if (nrow(ind_i) > 0) geom_sf(data = ind_i, aes(colour = "Indicator record"),
                                   shape = 4, size = 2.5, alpha = 0.8, stroke = 1)} +
    scale_colour_manual(
      values = c("GeMS point" = "#d95f02", "Historical point" = "#7570b3",
                 "Indicator record" = "#e6ab02"),
      name = "Points") +
    coord_sf(xlim = c(bb_i["xmin"], bb_i["xmax"]),
             ylim = c(bb_i["ymin"], bb_i["ymax"]), crs = 27700) +
    labs(title = sprintf("%s — oyster (all sources)", box_i$name),
         subtitle = sprintf("GeMS py:%d pt:%d | Hist py:%d pt:%d | Ind:%d",
                            nrow(gpy_i), nrow(gpt_i), nrow(hpy_i),
                            nrow(hpt_i), nrow(ind_i))) +
    map_theme

  fname <- paste0("figures/habitats/oyster_cluster_all_",
                  gsub("[^a-z0-9]+", "_", tolower(box_i$name)), ".png")
  ggsave(fname, p_i, width = 12, height = 10, units = "in", dpi = 300)
  cat(sprintf("  Saved: %s\n", fname))
}


## ---- 3d. Harmonised oyster layer ----

cat("\n-- 3d. Harmonised oyster layer --\n")

## Known reef extent: GeMS polygons (highest confidence)
oy_known <- gems_oy_py %>% select(source, Shape)
cat(sprintf("  Known reef extent polygons (GeMS): %d\n", nrow(oy_known)))

## Historical reef polygons: flag as historical
if (nrow(hist_oy_py) > 0) {
  oy_historical <- hist_oy_py %>% select(source, geometry)
  cat(sprintf("  Historical reef polygons: %d\n", nrow(oy_historical)))
} else {
  oy_historical <- st_sf(source = character(0), geometry = st_sfc(crs = 27700))
  cat("  No historical reef polygons\n")
}

## All oyster points: GeMS + historical + indicator
oy_points <- bind_rows(
  if (nrow(gems_oy_pt) > 0) gems_oy_pt %>% select(source, Shape),
  if (nrow(hist_oy_pt) > 0) hist_oy_pt %>% select(source, geometry),
  if (nrow(ind_oy) > 0) ind_oy %>% select(source, geom)
)
cat(sprintf("  All oyster observation points: %d\n", nrow(oy_points)))

## Save harmonised outputs
if (nrow(oy_known) > 0) {
  st_write(oy_known, "data/processed/species/oyster_known_extent.gpkg",
           delete_dsn = TRUE, quiet = TRUE)
  cat("  Saved: data/processed/species/oyster_known_extent.gpkg\n")
}

if (nrow(oy_historical) > 0) {
  st_write(oy_historical, "data/processed/species/oyster_historical_extent.gpkg",
           delete_dsn = TRUE, quiet = TRUE)
  cat("  Saved: data/processed/species/oyster_historical_extent.gpkg\n")
}

if (nrow(oy_points) > 0) {
  st_write(oy_points, "data/processed/species/oyster_observation_points.gpkg",
           delete_dsn = TRUE, quiet = TRUE)
  cat("  Saved: data/processed/species/oyster_observation_points.gpkg\n")
}


## ---- 3e. Solway Native Oyster Potential (modelled) ----

cat("\n-- 3e. Solway Native Oyster Potential (modelled) --\n")

## Load modelled oyster potential from 04_physical.R
solway_oy_pot_path <- "outputs/habitats/oyster_potential_solway.gpkg"
stopifnot("Solway oyster potential not found — run 04_physical.R section 9 first" =
            file.exists(solway_oy_pot_path))

solway_oy_pot <- st_read(solway_oy_pot_path, quiet = TRUE) %>% st_make_valid()
if (st_crs(solway_oy_pot)$epsg != 27700) {
  solway_oy_pot <- st_transform(solway_oy_pot, 27700)
}
cat(sprintf("  Solway oyster potential polygons: %d\n", nrow(solway_oy_pot)))
cat(sprintf("  Total potential area: %.1f ha\n",
            sum(as.numeric(st_area(solway_oy_pot))) / 10000))

## --- Whole-AOI seascape map: all oyster data + modelled potential ---

p_oy_full <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
          linewidth = 0.3) +
  geom_sf(data = scamp_coastline, colour = "grey30", linewidth = 0.3) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  geom_sf(data = solway_oy_pot, aes(fill = "Modelled potential"),
          colour = NA, alpha = 0.3) +
  {if (nrow(oy_known) > 0)
    geom_sf(data = oy_known, aes(fill = "Known extent (GeMS)"),
            colour = NA, alpha = 0.6)} +
  {if (nrow(oy_historical) > 0)
    geom_sf(data = oy_historical, aes(fill = "Historical extent"),
            colour = NA, alpha = 0.4)} +
  scale_fill_manual(
    values = c("Modelled potential" = "#d6604d",
               "Known extent (GeMS)"  = "#e7298a",
               "Historical extent"    = "#a6761d"),
    name = "Polygons") +
  {if (nrow(oy_points) > 0)
    geom_sf(data = oy_points, aes(colour = source, shape = source),
            size = 1.5, alpha = 0.7)} +
  scale_colour_manual(
    values = c("GeMS habitat point" = "#d95f02", "Historical" = "#7570b3",
               "Indicator records" = "#e6ab02"),
    name = "Points") +
  scale_shape_manual(
    values = c("GeMS habitat point" = 16, "Historical" = 17,
               "Indicator records" = 4),
    name = "Points") +
  geom_sf(data = oy_boxes, fill = NA, colour = "red", linewidth = 0.6,
          linetype = "dashed") +
  geom_sf_text(data = oy_boxes, aes(label = name), colour = "red",
               size = 3, fontface = "bold") +
  annotation_scale(location = "bl", width_hint = 0.2) +
  annotation_north_arrow(location = "tl", style = north_arrow_minimal(),
                         height = unit(1, "cm"), width = unit(1, "cm")) +
  labs(title = "Native Oyster — Full Seascape (SCAMP AOI)",
       subtitle = sprintf(
         "Modelled: %d poly | Known: %d | Historical: %d | Points: %d",
         nrow(solway_oy_pot),
         nrow(oy_known), nrow(oy_historical), nrow(oy_points))) +
  map_theme + coord_sf(crs = 27700)

ggsave("figures/habitats/oyster_full_seascape.png", p_oy_full,
       width = 14, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/habitats/oyster_full_seascape.png\n")

## --- Per-cluster maps ---

for (i in seq_len(nrow(oy_boxes))) {
  box_i <- oy_boxes[i, ]
  bb_i  <- st_bbox(box_i)

  sol_i <- solway_oy_pot[lengths(st_intersects(solway_oy_pot, box_i)) > 0, ]
  kn_i  <- oy_known[lengths(st_intersects(oy_known, box_i)) > 0, ]
  hi_i  <- oy_historical[lengths(st_intersects(oy_historical, box_i)) > 0, ]
  pt_i  <- oy_points[lengths(st_intersects(oy_points, box_i)) > 0, ]

  p_i <- ggplot() +
    geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
    geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
            linewidth = 0.3) +
    geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
    {if (nrow(sol_i) > 0)
      geom_sf(data = sol_i, aes(fill = "Modelled potential"),
              colour = NA, alpha = 0.3)} +
    {if (nrow(kn_i) > 0)
      geom_sf(data = kn_i, aes(fill = "Known extent (GeMS)"),
              colour = NA, alpha = 0.6)} +
    {if (nrow(hi_i) > 0)
      geom_sf(data = hi_i, aes(fill = "Historical extent"),
              colour = NA, alpha = 0.4)} +
    scale_fill_manual(
      values = c("Modelled potential" = "#d6604d",
                 "Known extent (GeMS)"  = "#e7298a",
                 "Historical extent"    = "#a6761d"),
      name = "Polygons") +
    {if (nrow(pt_i) > 0)
      geom_sf(data = pt_i, aes(colour = source, shape = source),
              size = 2.5, alpha = 0.8)} +
    scale_colour_manual(
      values = c("GeMS habitat point" = "#d95f02", "Historical" = "#7570b3",
                 "Indicator records" = "#e6ab02"),
      name = "Points") +
    scale_shape_manual(
      values = c("GeMS habitat point" = 16, "Historical" = 17,
                 "Indicator records" = 4),
      name = "Points") +
    coord_sf(xlim = c(bb_i["xmin"], bb_i["xmax"]),
             ylim = c(bb_i["ymin"], bb_i["ymax"]), crs = 27700) +
    labs(title = sprintf("%s — oyster full seascape", box_i$name),
         subtitle = sprintf(
           "Modelled: %d | Known: %d | Hist: %d | Points: %d",
           nrow(sol_i), nrow(kn_i), nrow(hi_i), nrow(pt_i))) +
    map_theme

  fname <- paste0("figures/habitats/oyster_full_seascape_",
                  gsub("[^a-z0-9]+", "_", tolower(box_i$name)), ".png")
  ggsave(fname, p_i, width = 12, height = 10, units = "in", dpi = 300)
  cat(sprintf("  Saved: %s\n", fname))
}


## ---- 3f. Broadscale Habitat Classification ----

cat("\n-- 3f. Broadscale Habitat Classification --\n")

## TODO: Define subtidal mixed sediments (EUNIS A5.4) classification
## User to define broadscale habitat type source and filtering criteria

cat("  [placeholder — awaiting habitat classification input]\n")


cat("\n=== Oyster section complete ===\n")


##=====================================================================
## 4. SAND DUNES
##=====================================================================

cat("\n\n##=== 4. SAND DUNES ===##\n")


## ---- 4a. HabMoS Sand Dune Habitat ----

cat("\n-- 4a. HabMoS Sand Dune Habitat --\n")

habmos_path <- "data/raw/habitats/General/Habitat_Map_of_Scotland_2181603415421412055.gpkg"
stopifnot("HabMoS geopackage not found" = file.exists(habmos_path))

habmos_layers <- st_layers(habmos_path)
cat("  Available layers:\n")
for (i in seq_along(habmos_layers$name)) {
  cat(sprintf("    %d: %s (%d features, %s)\n",
              i, habmos_layers$name[i], habmos_layers$features[i],
              habmos_layers$geomtype[i]))
}

habmos_sf <- st_read(habmos_path, layer = habmos_layers$name[1], quiet = TRUE)
cat(sprintf("\n  Loaded: %s | %d features | CRS: %s\n",
            habmos_layers$name[1], nrow(habmos_sf), st_crs(habmos_sf)$input))
cat(sprintf("  Columns: %s\n", paste(names(habmos_sf), collapse = ", ")))

## Print unique values of classification columns to identify sand dune classes
geom_col_hm <- attr(habmos_sf, "sf_column")
class_cols_hm <- setdiff(names(habmos_sf), geom_col_hm)
for (cc in class_cols_hm) {
  vals <- sort(unique(habmos_sf[[cc]]))
  if (length(vals) <= 50) {
    cat(sprintf("\n  %s (%d unique):\n", cc, length(vals)))
    for (v in vals) cat(sprintf("    - %s\n", v))
  } else {
    cat(sprintf("\n  %s (%d unique — too many to list, showing first 30):\n", cc, length(vals)))
    for (v in vals[1:30]) cat(sprintf("    - %s\n", v))
  }
}

## Clip to SCAMP AOI (transform AOI to source CRS)
aoi_habmos <- st_transform(scamp_aoi_sf, st_crs(habmos_sf))
habmos_clip <- st_intersection(habmos_sf, aoi_habmos)
cat(sprintf("\n  Features within AOI: %d of %d\n", nrow(habmos_clip), nrow(habmos_sf)))

## Filter to sand dune habitats
## Search for dune-related classes across all columns
dune_pattern <- "[Dd]une|[Mm]achair|[Ss]trandline|[Ss]and.?[Dd]une"
dune_hits <- rep(FALSE, nrow(habmos_clip))
for (cc in class_cols_hm) {
  if (is.character(habmos_clip[[cc]]) || is.factor(habmos_clip[[cc]])) {
    dune_hits <- dune_hits | grepl(dune_pattern, habmos_clip[[cc]], ignore.case = TRUE)
  }
}
habmos_dunes <- habmos_clip[dune_hits, ]
print(unique(habmos_dunes$HABITAT_NAME))
habmos_dunes_bng <- st_transform(habmos_dunes, 27700) %>% st_make_valid()

cat(sprintf("  Sand dune features in AOI: %d\n", nrow(habmos_dunes_bng)))

## Print breakdown of dune classes
for (cc in class_cols_hm) {
  if (is.character(habmos_dunes_bng[[cc]]) || is.factor(habmos_dunes_bng[[cc]])) {
    dune_vals <- sort(unique(habmos_dunes_bng[[cc]]))
    dune_vals <- dune_vals[grepl(dune_pattern, dune_vals, ignore.case = TRUE)]
    if (length(dune_vals) > 0) {
      cat(sprintf("\n  Dune classes in %s:\n", cc))
      for (v in dune_vals) {
        n <- sum(habmos_dunes_bng[[cc]] == v, na.rm = TRUE)
        cat(sprintf("    - %s  [%d features]\n", v, n))
      }
    }
  }
}

habmos_dunes_bng$area_ha <- as.numeric(st_area(habmos_dunes_bng)) / 10000
cat(sprintf("  Total sand dune area: %.1f ha\n", sum(habmos_dunes_bng$area_ha)))

## Save clipped dune polygons
dir.create("outputs/habitats", showWarnings = FALSE, recursive = TRUE)
st_write(habmos_dunes_bng, "outputs/habitats/sand_dunes_habmos_solway.gpkg",
         delete_dsn = TRUE, quiet = TRUE)
cat("  Saved: outputs/habitats/sand_dunes_habmos_solway.gpkg\n")

## Simple plot of sand dune polygons
p_dunes_simple <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
          linewidth = 0.3) +
  geom_sf(data = scamp_coastline, colour = "grey30", linewidth = 0.3) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  geom_sf(data = habmos_dunes_bng, fill = "#f4a582", colour = "#d6604d",
          alpha = 0.6, linewidth = 0.3) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  annotation_north_arrow(location = "tl", style = north_arrow_minimal(),
                         height = unit(1, "cm"), width = unit(1, "cm")) +
  labs(title = "Sand Dune Habitat — HabMoS (SCAMP AOI)",
       subtitle = sprintf("%d polygons | %.1f ha",
                          nrow(habmos_dunes_bng),
                          sum(habmos_dunes_bng$area_ha))) +
  map_theme + coord_sf(crs = 27700)

ggsave("figures/habitats/sand_dunes_habmos.png", p_dunes_simple,
       width = 14, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/habitats/sand_dunes_habmos.png\n")

## Cluster sand dune polygons into hotspot boxes
dune_centroids <- st_centroid(habmos_dunes_bng)
dune_buf <- st_buffer(dune_centroids, dist = 2000) %>% st_union() %>%
  st_cast("POLYGON")
dune_cids <- st_intersects(dune_centroids, dune_buf)
habmos_dunes_bng$cluster <- sapply(dune_cids, `[`, 1)

dune_boxes <- do.call(rbind, lapply(sort(unique(habmos_dunes_bng$cluster)), function(cl) {
  bb <- st_bbox(st_buffer(habmos_dunes_bng[habmos_dunes_bng$cluster == cl, ], dist = 1000))
  st_sf(
    name = sprintf("Dune_%d", cl),
    geometry = st_sfc(st_polygon(list(matrix(c(
      bb["xmin"], bb["ymin"],
      bb["xmax"], bb["ymin"],
      bb["xmax"], bb["ymax"],
      bb["xmin"], bb["ymax"],
      bb["xmin"], bb["ymin"]
    ), ncol = 2, byrow = TRUE))), crs = 27700)
  )
}))

## Merge overlapping boxes
db_buf <- st_buffer(dune_boxes, dist = 1000) %>% st_union() %>%
  st_cast("POLYGON")
db_group <- sapply(st_intersects(dune_boxes, db_buf), `[`, 1)

dune_boxes <- do.call(rbind, lapply(sort(unique(db_group)), function(g) {
  members <- dune_boxes[db_group == g, ]
  merged_bb <- st_bbox(members)
  st_sf(
    name = sprintf("Dune_%d", g),
    geometry = st_sfc(st_polygon(list(matrix(c(
      merged_bb["xmin"], merged_bb["ymin"],
      merged_bb["xmax"], merged_bb["ymin"],
      merged_bb["xmax"], merged_bb["ymax"],
      merged_bb["xmin"], merged_bb["ymax"],
      merged_bb["xmin"], merged_bb["ymin"]
    ), ncol = 2, byrow = TRUE))), crs = 27700)
  )
}))

## Name boxes west-to-east
dune_boxes <- dune_boxes[order(st_coordinates(st_centroid(dune_boxes))[, 1]), ]
dune_boxes$name <- sprintf("Dune cluster %d", seq_len(nrow(dune_boxes)))

cat(sprintf("  Sand dune clusters: %d\n", nrow(dune_boxes)))
for (j in seq_len(nrow(dune_boxes))) {
  n_in <- sum(lengths(st_intersects(habmos_dunes_bng, dune_boxes[j, ])) > 0)
  cat(sprintf("    %s: %d polygons\n", dune_boxes$name[j], n_in))
}

## Plot dunes with cluster boxes
p_dunes_clusters <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
          linewidth = 0.3) +
  geom_sf(data = scamp_coastline, colour = "grey30", linewidth = 0.3) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  geom_sf(data = habmos_dunes_bng, fill = "#f4a582", colour = "#d6604d",
          alpha = 0.6, linewidth = 0.3) +
  geom_sf(data = dune_boxes, fill = NA, colour = "red", linewidth = 0.8,
          linetype = "dashed") +
  geom_sf_text(data = dune_boxes, aes(label = name), colour = "red",
               size = 3, fontface = "bold") +
  annotation_scale(location = "bl", width_hint = 0.2) +
  annotation_north_arrow(location = "tl", style = north_arrow_minimal(),
                         height = unit(1, "cm"), width = unit(1, "cm")) +
  labs(title = "Sand Dune Clusters — HabMoS (SCAMP AOI)",
       subtitle = sprintf("%d polygons | %d clusters | %.1f ha",
                          nrow(habmos_dunes_bng), nrow(dune_boxes),
                          sum(habmos_dunes_bng$area_ha))) +
  map_theme + coord_sf(crs = 27700)

ggsave("figures/habitats/sand_dunes_clusters.png", p_dunes_clusters,
       width = 14, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/habitats/sand_dunes_clusters.png\n")


## ---- 4b. Sand Dune Indicator Species ----

cat("\n-- 4b. Sand Dune Indicator Species --\n")

ind_path <- "outputs/species/indicator_records.gpkg"
stopifnot("Indicator records not found" = file.exists(ind_path))

ind_all_dune <- st_read(ind_path, quiet = TRUE)

## Print available habitat categories to identify dune-related ones
cat("  Habitat categories in indicator records:\n")
hab_counts <- ind_all_dune %>%
  st_drop_geometry() %>%
  count(habitat, sort = TRUE)
print(hab_counts)

## Filter to sand dune indicators
dune_ind <- ind_all_dune %>%
  filter(grepl("dune|machair| strandline", habitat, ignore.case = TRUE)) %>%
  st_transform(27700) %>%
  st_make_valid() %>%
  st_intersection(scamp_aoi_sf) %>%
  filter(grepl("Leymus|Arvicola|Ammophila", latin_name, ignore.case=TRUE))

cat(sprintf("  Sand dune indicator records in AOI: %d\n", nrow(dune_ind)))

if (nrow(dune_ind) > 0) {
  cat("  Species breakdown:\n")
  print(dune_ind %>% st_drop_geometry() %>% count(latin_name, sort = TRUE))
}


## ---- 4c. Sand Dune Maps ----

cat("\n-- 4c. Sand Dune Maps --\n")

## AOI-wide map: HabMoS dunes + indicator species
p_dunes <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
          linewidth = 0.3) +
  geom_sf(data = scamp_coastline, colour = "grey30", linewidth = 0.3) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  {if (nrow(habmos_dunes_bng) > 0)
    geom_sf(data = habmos_dunes_bng, fill = "#f4a582", colour = NA, alpha = 0.5)} +
  {if (nrow(dune_ind) > 0)
    geom_sf(data = dune_ind, aes(colour = latin_name), shape = 16,
            size = 2, alpha = 0.8)} +
  scale_colour_brewer(palette = "Dark2", name = "Indicator species") +
  geom_sf(data = dune_boxes, fill = NA, colour = "red", linewidth = 0.6,
          linetype = "dashed") +
  geom_sf_text(data = dune_boxes, aes(label = name), colour = "red",
               size = 3, fontface = "bold") +
  annotation_scale(location = "bl", width_hint = 0.2) +
  annotation_north_arrow(location = "tl", style = north_arrow_minimal(),
                         height = unit(1, "cm"), width = unit(1, "cm")) +
  labs(title = "Sand Dunes — HabMoS Habitat + Indicator Species (SCAMP AOI)",
       subtitle = sprintf("HabMoS dune polygons: %d (%.1f ha) | Indicator records: %d",
                          nrow(habmos_dunes_bng),
                          sum(habmos_dunes_bng$area_ha, na.rm = TRUE),
                          nrow(dune_ind))) +
  map_theme + coord_sf(crs = 27700)

ggsave("figures/habitats/sand_dunes_seascape.png", p_dunes,
       width = 14, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/habitats/sand_dunes_seascape.png\n")

## Per-cluster maps
for (i in seq_len(nrow(dune_boxes))) {
  box_i <- dune_boxes[i, ]
  bb_i  <- st_bbox(box_i)

  dune_poly_i <- habmos_dunes_bng[lengths(st_intersects(habmos_dunes_bng, box_i)) > 0, ]
  dune_ind_i  <- dune_ind[lengths(st_intersects(dune_ind, box_i)) > 0, ]

  if (nrow(dune_poly_i) == 0 && nrow(dune_ind_i) == 0) next

  p_i <- ggplot() +
    geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
    geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
            linewidth = 0.3) +
    geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
    {if (nrow(dune_poly_i) > 0)
      geom_sf(data = dune_poly_i, fill = "#f4a582", colour = "#d6604d",
              alpha = 0.5, linewidth = 0.3)} +
    {if (nrow(dune_ind_i) > 0)
      geom_sf(data = dune_ind_i, aes(colour = latin_name), shape = 16,
              size = 3, alpha = 0.8)} +
    scale_colour_brewer(palette = "Dark2", name = "Indicator species") +
    coord_sf(xlim = c(bb_i["xmin"], bb_i["xmax"]),
             ylim = c(bb_i["ymin"], bb_i["ymax"]), crs = 27700) +
    labs(title = sprintf("%s — Sand Dunes", box_i$name),
         subtitle = sprintf("HabMoS polygons: %d | Indicator records: %d",
                            nrow(dune_poly_i), nrow(dune_ind_i))) +
    map_theme

  fname <- paste0("figures/habitats/sand_dunes_cluster_",
                  gsub("[^a-z0-9]+", "_", tolower(box_i$name)), ".png")
  ggsave(fname, p_i, width = 12, height = 10, units = "in", dpi = 300)
  cat(sprintf("  Saved: %s\n", fname))
}

## ---- 4d. Sand Dunes Potential — Solway ----

cat("\n-- 4d. Sand Dunes Potential — Solway --\n")

## 1. Intertidal Substrate Foreshore
foreshore_path <- "data/raw/habitats/Sand dunes/Intertidal+Substrate+Foreshore+(England+and+Scotland)/Intertidal Substrate Foreshore (England and Scotland)/DEFR00000009.shp"
stopifnot("Intertidal Substrate Foreshore not found" = file.exists(foreshore_path))

foreshore_sf <- st_read(foreshore_path, quiet = TRUE) %>% st_make_valid()
cat(sprintf("  Foreshore loaded: %d features | CRS: %s\n",
            nrow(foreshore_sf), st_crs(foreshore_sf)$input))
cat(sprintf("  Columns: %s\n", paste(names(foreshore_sf), collapse = ", ")))

## CRS missing from .prj — coordinates are BNG, assign EPSG:27700
if (is.na(st_crs(foreshore_sf))) {
  cat("  CRS is NA — assigning EPSG:27700 (coordinates consistent with BNG)\n")
  st_crs(foreshore_sf) <- 27700
}

## Clip to AOI
foreshore_clip <- st_intersection(foreshore_sf, scamp_aoi_sf) %>% st_make_valid()
cat(sprintf("  Foreshore in AOI: %d features\n", nrow(foreshore_clip)))

## Quick check plot: clipped foreshore vs AOI
p_fs_check <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  geom_sf(data = foreshore_clip, fill = "#fee08b", colour = NA, alpha = 0.4) +
  labs(title = "Foreshore CRS check — clipped to SCAMP AOI") +
  map_theme + coord_sf(crs = 27700)

ggsave("figures/habitats/foreshore_crs_check.png", p_fs_check,
       width = 14, height = 10, units = "in", dpi = 300)
cat("  Saved: figures/habitats/foreshore_crs_check.png (check alignment)\n")

## Print unique substrate classes
geom_col_fs <- attr(foreshore_clip, "sf_column")
class_cols_fs <- setdiff(names(foreshore_clip), geom_col_fs)
for (cc in class_cols_fs) {
  if (is.character(foreshore_clip[[cc]]) || is.factor(foreshore_clip[[cc]])) {
    vals <- sort(unique(foreshore_clip[[cc]]))
    if (length(vals) <= 30) {
      cat(sprintf("\n  %s (%d unique):\n", cc, length(vals)))
      for (v in vals) {
        n <- sum(foreshore_clip[[cc]] == v, na.rm = TRUE)
        cat(sprintf("    - %s  [%d]\n", v, n))
      }
    }
  }
}

## 2. Dynamic Coast — MHWS 1970
mhws1970_path <- "data/raw/habitats/Sand dunes/DYNAMICCOAST2_SCOTLAND_SHP_27700/shp/DC2_MHWS_1970.shp"
stopifnot("DC2_MHWS_1970 not found" = file.exists(mhws1970_path))

mhws_1970 <- st_read(mhws1970_path, quiet = TRUE) %>% st_make_valid()
cat(sprintf("\n  MHWS 1970 loaded: %d features | CRS: %s\n",
            nrow(mhws_1970), st_crs(mhws_1970)$input))
if (st_crs(mhws_1970)$epsg != 27700) mhws_1970 <- st_transform(mhws_1970, 27700)

aoi_dc <- st_transform(scamp_aoi_sf, st_crs(mhws_1970))
mhws_1970_clip <- st_intersection(mhws_1970, aoi_dc) %>% st_make_valid()
cat(sprintf("  MHWS 1970 in AOI: %d features\n", nrow(mhws_1970_clip)))

## 3. Dynamic Coast — MHWS Modern
mhws_modern_path <- "data/raw/habitats/Sand dunes/DYNAMICCOAST2_SCOTLAND_SHP_27700/shp/DC2_MHWS_MODERN.shp"
stopifnot("DC2_MHWS_MODERN not found" = file.exists(mhws_modern_path))

mhws_modern <- st_read(mhws_modern_path, quiet = TRUE) %>% st_make_valid()
cat(sprintf("\n  MHWS Modern loaded: %d features | CRS: %s\n",
            nrow(mhws_modern), st_crs(mhws_modern)$input))
if (st_crs(mhws_modern)$epsg != 27700) mhws_modern <- st_transform(mhws_modern, 27700)

mhws_modern_clip <- st_intersection(mhws_modern, st_transform(scamp_aoi_sf, st_crs(mhws_modern))) %>%
  st_make_valid()
cat(sprintf("  MHWS Modern in AOI: %d features\n", nrow(mhws_modern_clip)))

## 2b. RCP8.5 Future Erosion 2050
ero2050_path <- "data/raw/habitats/Sand dunes/DYNAMICCOAST2_SCOTLAND_SHP_27700/shp/DC2_RCP8_FUTURE_ERO_2050_PUB.shp"
stopifnot("DC2_RCP8_FUTURE_ERO_2050 not found" = file.exists(ero2050_path))
ero_2050 <- st_read(ero2050_path, quiet = TRUE) %>% st_make_valid()
cat(sprintf("  Erosion 2050 loaded: %d features | CRS: %s\n",
            nrow(ero_2050), st_crs(ero_2050)$input))
if (st_crs(ero_2050)$epsg != 27700) ero_2050 <- st_transform(ero_2050, 27700)
ero_2050_clip <- st_intersection(ero_2050, scamp_aoi_sf) %>% st_make_valid()
cat(sprintf("  Erosion 2050 in AOI: %d features\n", nrow(ero_2050_clip)))

## 4. COP30 DEM
cop30_path <- "data/raw/physical/rasters_COP30/output_hh.tif"
stopifnot("COP30 DEM not found" = file.exists(cop30_path))

cop30_raw <- rast(cop30_path)
cat(sprintf("\n  COP30 DEM: %d x %d | res: %.2f | CRS: %s\n",
            ncol(cop30_raw), nrow(cop30_raw), res(cop30_raw)[1],
            crs(cop30_raw, describe = TRUE)$code))

aoi_cop <- vect(st_transform(scamp_aoi_sf, crs(cop30_raw)))
cop30_crop <- crop(cop30_raw, aoi_cop, mask = TRUE)
cop30_bng <- project(cop30_crop, "EPSG:27700")
names(cop30_bng) <- "elevation"

cat(sprintf("  COP30 cropped: %d x %d | Elevation range: %.1f to %.1f m\n",
            ncol(cop30_bng), nrow(cop30_bng),
            global(cop30_bng, "min", na.rm = TRUE)[[1]],
            global(cop30_bng, "max", na.rm = TRUE)[[1]]))

writeRaster(cop30_bng, "outputs/physical/cop30_dem_solway.tif",
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))
cat("  Saved: outputs/physical/cop30_dem_solway.tif\n")

## Erosion test zones (WGS84 → BNG)
erosion_zones_wgs <- list(
  list(name = "Caerlaverock",  xmin = -3.57,  xmax = -3.42, ymin = 54.94, ymax = 55.00),
  list(name = "Auchencairn",   xmin = -4.456, xmax = -4.34, ymin = 54.822, ymax = 54.90),
  list(name = "Luce Bay",      xmin = -5.00,  xmax = -4.77, ymin = 54.78,  ymax = 54.88)
)

erosion_boxes <- do.call(rbind, lapply(erosion_zones_wgs, function(z) {
  st_sf(
    name = z$name,
    geometry = st_sfc(st_polygon(list(matrix(c(
      z$xmin, z$ymin,
      z$xmax, z$ymin,
      z$xmax, z$ymax,
      z$xmin, z$ymax,
      z$xmin, z$ymin
    ), ncol = 2, byrow = TRUE))), crs = 4326)
  )
}))
erosion_boxes <- st_transform(erosion_boxes, 27700)

cat(sprintf("  Erosion test zones: %d\n", nrow(erosion_boxes)))
for (i in seq_len(nrow(erosion_boxes))) {
  bb_i <- st_bbox(erosion_boxes[i, ])
  cat(sprintf("    %s: %.0f–%.0f E, %.0f–%.0f N\n",
              erosion_boxes$name[i], bb_i["xmin"], bb_i["xmax"],
              bb_i["ymin"], bb_i["ymax"]))
}

ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue",
          linewidth = 0.3) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
  geom_sf(data = erosion_boxes, fill = NA, colour = "red", linewidth = 0.8) +
  geom_sf_text(data = erosion_boxes, aes(label = name), colour = "red",
               size = 3, fontface = "bold") +
  labs(title = "Erosion Test Zones") +
  map_theme + coord_sf(crs = 27700)

# --- Simplify modern MHWS: compare tolerances ---
mhws_smooth_10 <- st_simplify(mhws_modern_clip, dTolerance = 10)
mhws_smooth_20 <- st_simplify(mhws_modern_clip, dTolerance = 20)
mhws_smooth_50 <- st_simplify(mhws_modern_clip, dTolerance = 100)

bb_caer <- st_bbox(erosion_boxes[erosion_boxes$name == "Caerlaverock", ])
ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue", linewidth = 0.3) +
  geom_sf(data = mhws_modern_clip, colour = "black", linewidth = 0.5) +
  geom_sf(data = mhws_smooth_10, colour = "blue", linewidth = 0.4) +
  geom_sf(data = mhws_smooth_20, colour = "orange", linewidth = 0.4) +
  geom_sf(data = mhws_smooth_50, colour = "red", linewidth = 0.4) +
  coord_sf(xlim = c(bb_caer["xmin"], bb_caer["xmax"]),
           ylim = c(bb_caer["ymin"], bb_caer["ymax"]), crs = 27700) +
  labs(title = "MHWS Modern Simplify — Caerlaverock",
       subtitle = "Black = original, Blue = 10m, Orange = 20m, Red = 50m") +
  map_theme

bb_auch <- st_bbox(erosion_boxes[erosion_boxes$name == "Auchencairn", ])
ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue", linewidth = 0.3) +
  geom_sf(data = mhws_modern_clip, colour = "black", linewidth = 0.5) +
  geom_sf(data = mhws_smooth_10, colour = "blue", linewidth = 0.4) +
  geom_sf(data = mhws_smooth_20, colour = "orange", linewidth = 0.4) +
  geom_sf(data = mhws_smooth_50, colour = "red", linewidth = 0.4) +
  coord_sf(xlim = c(bb_auch["xmin"], bb_auch["xmax"]),
           ylim = c(bb_auch["ymin"], bb_auch["ymax"]), crs = 27700) +
  labs(title = "MHWS Modern Simplify — Auchencairn",
       subtitle = "Black = original, Blue = 10m, Orange = 20m, Red = 50m") +
  map_theme

bb_luce <- st_bbox(erosion_boxes[erosion_boxes$name == "Luce Bay", ])
ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue", linewidth = 0.3) +
  geom_sf(data = mhws_modern_clip, colour = "black", linewidth = 0.5) +
  geom_sf(data = mhws_smooth_10, colour = "blue", linewidth = 0.4) +
  geom_sf(data = mhws_smooth_20, colour = "orange", linewidth = 0.4) +
  geom_sf(data = mhws_smooth_50, colour = "red", linewidth = 0.4) +
  coord_sf(xlim = c(bb_luce["xmin"], bb_luce["xmax"]),
           ylim = c(bb_luce["ymin"], bb_luce["ymax"]), crs = 27700) +
  labs(title = "MHWS Modern Simplify — Luce Bay",
       subtitle = "Black = original, Blue = 10m, Orange = 20m, Red = 50m") +
  map_theme

# --- Perpendicular transects from smoothed MHWS (50m tolerance) ---
# TEST: Luce Bay only — clip MHWS to Luce Bay bbox before processing
luce_bbox_geom <- st_as_sfc(st_bbox(erosion_boxes[erosion_boxes$name == "Luce Bay", ]))
st_crs(luce_bbox_geom) <- 27700
mhws_test <- st_intersection(mhws_modern_clip, luce_bbox_geom) %>% st_make_valid()
cat(sprintf("Luce Bay test: %d MHWS features\n", nrow(mhws_test)))

mhws_smooth <- st_simplify(mhws_test, dTolerance = 100)
mhws_smooth_merged <- st_line_merge(st_union(st_cast(st_geometry(mhws_smooth), "LINESTRING")))
coords_all <- st_coordinates(mhws_smooth_merged)

transect_length <- 1000
bearing_offset <- 50
sample_spacing <- 50

transects_list <- lapply(unique(coords_all[, "L1"]), function(lid) {
  cc <- coords_all[coords_all[, "L1"] == lid, 1:2, drop = FALSE]
  seg_dx <- diff(cc[, 1])
  seg_dy <- diff(cc[, 2])
  seg_len <- sqrt(seg_dx^2 + seg_dy^2)
  cum_len <- c(0, cumsum(seg_len))
  total_len <- cum_len[length(cum_len)]
  if (total_len < sample_spacing) return(NULL)

  sample_dists <- seq(0, total_len, by = sample_spacing)
  seg_idx <- findInterval(sample_dists, cum_len, rightmost.closed = TRUE)
  seg_idx <- pmin(seg_idx, length(seg_len))

  frac <- (sample_dists - cum_len[seg_idx]) / seg_len[seg_idx]
  frac[is.nan(frac)] <- 0
  px <- cc[seg_idx, 1] + frac * seg_dx[seg_idx]
  py <- cc[seg_idx, 2] + frac * seg_dy[seg_idx]

  d_back <- pmax(0, sample_dists - bearing_offset)
  d_fwd  <- pmin(total_len, sample_dists + bearing_offset)
  seg_b <- pmin(findInterval(d_back, cum_len, rightmost.closed = TRUE), length(seg_len))
  seg_f <- pmin(findInterval(d_fwd, cum_len, rightmost.closed = TRUE), length(seg_len))
  frac_b <- (d_back - cum_len[seg_b]) / seg_len[seg_b]
  frac_f <- (d_fwd  - cum_len[seg_f]) / seg_len[seg_f]
  frac_b[is.nan(frac_b)] <- 0
  frac_f[is.nan(frac_f)] <- 0
  bx <- cc[seg_b, 1] + frac_b * seg_dx[seg_b]
  by_ <- cc[seg_b, 2] + frac_b * seg_dy[seg_b]
  fx <- cc[seg_f, 1] + frac_f * seg_dx[seg_f]
  fy <- cc[seg_f, 2] + frac_f * seg_dy[seg_f]

  bearing <- atan2(fx - bx, fy - by_)
  perp <- bearing + pi / 2

  x1 <- px + transect_length * sin(perp)
  y1 <- py + transect_length * cos(perp)
  x2 <- px - transect_length * sin(perp)
  y2 <- py - transect_length * cos(perp)

  cbind(x1, y1, x2, y2, px, py)
})

trans_mat <- do.call(rbind, Filter(Negate(is.null), transects_list))

transects_sf <- st_sf(
  id = seq_len(nrow(trans_mat)),
  geometry = st_sfc(lapply(seq_len(nrow(trans_mat)), function(i) {
    st_linestring(matrix(trans_mat[i, 1:4], ncol = 2, byrow = TRUE))
  }), crs = 27700)
)
sample_pts_sf <- st_sf(
  id = seq_len(nrow(trans_mat)),
  geometry = st_sfc(lapply(seq_len(nrow(trans_mat)), function(i) {
    st_point(trans_mat[i, 5:6])
  }), crs = 27700)
)
cat(sprintf("Built %d transects at %dm spacing, %dm either side\n",
            nrow(transects_sf), sample_spacing, transect_length))

# --- Luce Bay transect test plot ---
ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue", linewidth = 0.3) +
  geom_sf(data = mhws_test, colour = "grey40", linewidth = 0.3) +
  geom_sf(data = mhws_smooth, colour = "black", linewidth = 0.6) +
  geom_sf(data = transects_sf, colour = "red", linewidth = 0.3, alpha = 0.6) +
  geom_sf(data = sample_pts_sf, colour = "blue", size = 0.5) +
  coord_sf(xlim = c(bb_luce["xmin"], bb_luce["xmax"]),
           ylim = c(bb_luce["ymin"], bb_luce["ymax"]), crs = 27700) +
  labs(title = "Perpendicular Transects — Luce Bay (test)",
       subtitle = "Grey = original, Black = smoothed (50m), Red = transects at 50m spacing") +
  map_theme

# --- Filter MHWS: remove saltmarsh zones and non-sandy foreshores ---

# Sandy foreshore mask: filter FORE_DESC for sand-containing classes
cat("\nFORE_DESC classes in AOI:\n")
print(sort(unique(foreshore_clip$FORE_DESC)))

sandy_foreshore <- foreshore_clip[grepl("SAND", foreshore_clip$FORE_DESC,
                                        ignore.case = TRUE), ]
cat(sprintf("  Sandy foreshore polygons: %d of %d\n",
            nrow(sandy_foreshore), nrow(foreshore_clip)))

# Strict sand-only foreshore: FORE_DESC == "SAND" exactly
sandy_strict <- foreshore_clip[foreshore_clip$FORE_DESC == "SAND", ]
cat(sprintf("  Strict sand-only foreshore polygons: %d of %d\n",
            nrow(sandy_strict), nrow(foreshore_clip)))

# Saltmarsh mask: load harmonised saltmarsh layer
sm_harmonised_sf <- st_read("data/processed/species/saltmarsh_harmonised.gpkg", quiet = TRUE)
sm_union <- st_union(sm_harmonised_sf) %>% st_make_valid()
cat(sprintf("  Saltmarsh polygons loaded: %d\n", nrow(sm_harmonised_sf)))

# Rule 1: Buffer sandy foreshore by 200m, keep MHWS segments that intersect
sandy_union <- st_union(sandy_foreshore) %>% st_buffer(200) %>% st_make_valid()
cat(sprintf("  Sandy foreshore buffer (200m) created\n"))

# Rule 2: Saltmarsh union for removal (200m buffer)
sm_union <- st_union(sm_harmonised_sf) %>% st_buffer(200) %>% st_make_valid()
cat(sprintf("  Saltmarsh buffer (200m) created\n"))

# Load 1890s MHWS
mhws1890_path <- "data/raw/habitats/Sand dunes/DYNAMICCOAST2_SCOTLAND_SHP_27700/shp/DC2_MHWS_1890.shp"
stopifnot("DC2_MHWS_1890 not found" = file.exists(mhws1890_path))
mhws_1890 <- st_read(mhws1890_path, quiet = TRUE) %>% st_make_valid()
cat(sprintf("  MHWS 1890 loaded: %d features | CRS: %s\n",
            nrow(mhws_1890), st_crs(mhws_1890)$input))
if (st_crs(mhws_1890)$epsg != 27700) mhws_1890 <- st_transform(mhws_1890, 27700)
mhws_1890_clip <- st_intersection(mhws_1890, scamp_aoi_sf) %>% st_make_valid()
cat(sprintf("  MHWS 1890 in AOI: %d features\n", nrow(mhws_1890_clip)))

# Filter: keep segments within 500m of sandy foreshore, then remove saltmarsh
filter_mhws <- function(mhws, sandy_buf, saltmarsh, label) {
  hits <- lengths(st_intersects(mhws, sandy_buf)) > 0
  within_sand <- mhws[hits, ]
  cat(sprintf("  %s: %d -> %d after sandy foreshore filter (500m buffer)\n",
              label, nrow(mhws), nrow(within_sand)))
  minus_sm <- st_difference(st_union(within_sand), saltmarsh) %>%
    st_make_valid() %>%
    st_cast("MULTILINESTRING") %>%
    st_cast("LINESTRING")
  result <- st_sf(geometry = minus_sm, crs = 27700)
  cat(sprintf("  %s: %d features after removing saltmarsh\n", label, nrow(result)))
  result
}

mhws_1890_filt <- filter_mhws(mhws_1890_clip, sandy_union, sm_union, "MHWS 1890")
mhws_1970_filt <- filter_mhws(mhws_1970_clip, sandy_union, sm_union, "MHWS 1970")
mhws_modern_filt <- filter_mhws(mhws_modern_clip, sandy_union, sm_union, "MHWS Modern")

# --- Vectorised transects on filtered modern MHWS ---
mhws_smooth <- st_simplify(mhws_modern_filt, dTolerance = 150)
mhws_smooth_merged <- st_line_merge(st_union(st_cast(st_geometry(mhws_smooth), "LINESTRING")))
coords_all <- st_coordinates(mhws_smooth_merged)

transect_length <- 200
bearing_offset <- 50
sample_spacing <- 50

transects_list <- lapply(unique(coords_all[, "L1"]), function(lid) {
  cc <- coords_all[coords_all[, "L1"] == lid, 1:2, drop = FALSE]
  seg_dx <- diff(cc[, 1])
  seg_dy <- diff(cc[, 2])
  seg_len <- sqrt(seg_dx^2 + seg_dy^2)
  cum_len <- c(0, cumsum(seg_len))
  total_len <- cum_len[length(cum_len)]
  if (total_len < sample_spacing) return(NULL)

  sample_dists <- seq(0, total_len, by = sample_spacing)
  seg_idx <- findInterval(sample_dists, cum_len, rightmost.closed = TRUE)
  seg_idx <- pmin(seg_idx, length(seg_len))

  frac <- (sample_dists - cum_len[seg_idx]) / seg_len[seg_idx]
  frac[is.nan(frac)] <- 0
  px <- cc[seg_idx, 1] + frac * seg_dx[seg_idx]
  py <- cc[seg_idx, 2] + frac * seg_dy[seg_idx]

  d_back <- pmax(0, sample_dists - bearing_offset)
  d_fwd  <- pmin(total_len, sample_dists + bearing_offset)
  seg_b <- pmin(findInterval(d_back, cum_len, rightmost.closed = TRUE), length(seg_len))
  seg_f <- pmin(findInterval(d_fwd, cum_len, rightmost.closed = TRUE), length(seg_len))
  frac_b <- (d_back - cum_len[seg_b]) / seg_len[seg_b]
  frac_f <- (d_fwd  - cum_len[seg_f]) / seg_len[seg_f]
  frac_b[is.nan(frac_b)] <- 0
  frac_f[is.nan(frac_f)] <- 0
  bx <- cc[seg_b, 1] + frac_b * seg_dx[seg_b]
  by_ <- cc[seg_b, 2] + frac_b * seg_dy[seg_b]
  fx <- cc[seg_f, 1] + frac_f * seg_dx[seg_f]
  fy <- cc[seg_f, 2] + frac_f * seg_dy[seg_f]

  bearing <- atan2(fx - bx, fy - by_)
  perp <- bearing + pi / 2

  x1 <- px + transect_length * sin(perp)
  y1 <- py + transect_length * cos(perp)
  x2 <- px - transect_length * sin(perp)
  y2 <- py - transect_length * cos(perp)

  cbind(x1, y1, x2, y2, px, py)
})

trans_mat <- do.call(rbind, Filter(Negate(is.null), transects_list))

transects_sf <- st_sf(
  id = seq_len(nrow(trans_mat)),
  geometry = st_sfc(lapply(seq_len(nrow(trans_mat)), function(i) {
    st_linestring(matrix(trans_mat[i, 1:4], ncol = 2, byrow = TRUE))
  }), crs = 27700)
)
sample_pts_sf <- st_sf(
  id = seq_len(nrow(trans_mat)),
  geometry = st_sfc(lapply(seq_len(nrow(trans_mat)), function(i) {
    st_point(trans_mat[i, 5:6])
  }), crs = 27700)
)
cat(sprintf("Built %d transects at %dm spacing, %dm either side\n",
            nrow(transects_sf), sample_spacing, transect_length))

# Remove transects that don't intersect sandy foreshore
hits_fs <- lengths(st_intersects(transects_sf, sandy_foreshore)) > 0
cat(sprintf("Transects intersecting sandy foreshore: %d of %d\n", sum(hits_fs), nrow(transects_sf)))
transects_sf <- transects_sf[hits_fs, ]
sample_pts_sf <- sample_pts_sf[hits_fs, ]

# --- Interactive mapview: transects + filtered MHWS ---
library(mapview)
mapview(transects_sf, color = "red", layer.name = "Transects") +
  mapview(mhws_modern_filt, color = "black", lwd = 2, layer.name = "MHWS Modern (filtered)") +
  mapview(mhws_1970_filt, color = "orange", layer.name = "MHWS 1970 (filtered)") +
  mapview(mhws_1890_filt, color = "brown", layer.name = "MHWS 1890 (filtered)") +
  mapview(sm_harmonised_sf, col.regions = "green", layer.name = "Saltmarsh (excluded)") +
  mapview(sandy_foreshore, col.regions = "yellow", layer.name = "Sandy foreshore") +
  mapview(ero_2050_clip, col.regions = "red", alpha.regions = 0.4, layer.name = "Erosion 2050") +
  mapview(habmos_dunes_bng, col.regions = "purple", layer.name = "Sand Dune Habitat (HabMoS)")

# ============================================================
# TEST: Overlap correction methods comparison
# ============================================================

# --- Step 1: Identify crossing transects ---
cross_mat <- st_intersects(transects_sf, transects_sf)
transects_sf$n_crosses <- lengths(cross_mat) - 1  # subtract self
cat(sprintf("Transects with overlaps: %d of %d\n",
            sum(transects_sf$n_crosses > 0), nrow(transects_sf)))

# --- Step 2: Find a test area with both overlapping and non-overlapping transects ---
# Pick the single most-overlapping transect and build a window around it
worst_idx <- which.max(transects_sf$n_crosses)
worst_xy <- st_coordinates(sample_pts_sf[worst_idx, ])
cat(sprintf("Highest overlap transect: id %d with %d crossings at %.0f, %.0f\n",
            worst_idx, transects_sf$n_crosses[worst_idx], worst_xy[1], worst_xy[2]))
test_radius <- 2000
test_bbox <- st_bbox(c(
  xmin = worst_xy[1] - test_radius,
  xmax = worst_xy[1] + test_radius,
  ymin = worst_xy[2] - test_radius,
  ymax = worst_xy[2] + test_radius
), crs = 27700)
test_area <- st_as_sfc(test_bbox)

# Select transects within test area
in_test <- lengths(st_intersects(transects_sf, test_area)) > 0
test_transects <- transects_sf[in_test, ]
test_pts <- sample_pts_sf[in_test, ]
cat(sprintf("Test area transects: %d (overlapping: %d, clean: %d)\n",
            nrow(test_transects),
            sum(test_transects$n_crosses > 0),
            sum(test_transects$n_crosses == 0)))

# --- Step 3: Determine seaward direction using COP30 DEM ---
# For each transect endpoint pair, the seaward end has lower elevation
# Extract DEM elevation at both ends of each transect

get_seaward_direction <- function(transects, dem) {
  coords_list <- lapply(seq_len(nrow(transects)), function(i) {
    cc <- st_coordinates(transects[i, ])
    data.frame(x1 = cc[1, 1], y1 = cc[1, 2], x2 = cc[2, 1], y2 = cc[2, 2])
  })
  coords_df <- do.call(rbind, coords_list)

  end1 <- vect(cbind(coords_df$x1, coords_df$y1), crs = "EPSG:27700")
  end2 <- vect(cbind(coords_df$x2, coords_df$y2), crs = "EPSG:27700")
  elev1 <- terra::extract(dem, end1)[, 2]
  elev2 <- terra::extract(dem, end2)[, 2]

  # Seaward direction vector: from high end to low end
  sea_dx <- ifelse(elev1 > elev2, coords_df$x2 - coords_df$x1, coords_df$x1 - coords_df$x2)
  sea_dy <- ifelse(elev1 > elev2, coords_df$y2 - coords_df$y1, coords_df$y1 - coords_df$y2)
  len <- sqrt(sea_dx^2 + sea_dy^2)
  data.frame(sea_dx = sea_dx / len, sea_dy = sea_dy / len,
             elev1 = elev1, elev2 = elev2)
}

# --- Step 3b: Calculate signed shoreline change per transect ---
calc_shoreline_change <- function(transects, mhws_mod, mhws_70, pts, dem) {
  sea_dir <- get_seaward_direction(transects, dem)

  results <- data.frame(
    id = seq_len(nrow(transects)),
    change_m = NA_real_
  )
  for (i in seq_len(nrow(transects))) {
    tr <- transects[i, ]
    pt_xy <- st_coordinates(pts[i, ])

    int_mod <- st_intersection(tr, mhws_mod)
    int_70 <- st_intersection(tr, mhws_70)
    if (nrow(int_mod) == 0 || nrow(int_70) == 0) next

    mod_pts <- st_coordinates(st_cast(st_geometry(int_mod), "POINT"))
    d_mod <- sqrt((mod_pts[, 1] - pt_xy[1])^2 + (mod_pts[, 2] - pt_xy[2])^2)
    mod_xy <- mod_pts[which.min(d_mod), 1:2]

    pts_70 <- st_coordinates(st_cast(st_geometry(int_70), "POINT"))
    d_70 <- sqrt((pts_70[, 1] - pt_xy[1])^2 + (pts_70[, 2] - pt_xy[2])^2)
    xy_70 <- pts_70[which.min(d_70), 1:2]

    # Vector from 1970 to modern intersection
    shift_dx <- mod_xy[1] - xy_70[1]
    shift_dy <- mod_xy[2] - xy_70[2]
    shift_dist <- sqrt(shift_dx^2 + shift_dy^2)

    # Dot product with seaward direction: positive = accretion, negative = erosion
    dot <- shift_dx * sea_dir$sea_dx[i] + shift_dy * sea_dir$sea_dy[i]
    results$change_m[i] <- sign(dot) * shift_dist
  }
  results
}

cat("Calculating signed shoreline change (test area)...\n")
test_change <- calc_shoreline_change(test_transects, mhws_modern_filt, mhws_1970_filt,
                                     test_pts, cop30_bng)
test_transects$change_m <- test_change$change_m
test_transects$has_change <- !is.na(test_change$change_m)
cat(sprintf("Transects with valid change measurement: %d of %d\n",
            sum(test_transects$has_change), nrow(test_transects)))

# --- Step 4: Assign overlap groups ---
test_valid <- test_transects[test_transects$has_change, ]
cross_valid <- st_intersects(test_valid, test_valid)

visited <- rep(FALSE, nrow(test_valid))
group_id <- rep(NA_integer_, nrow(test_valid))
g <- 0L
for (i in seq_len(nrow(test_valid))) {
  if (visited[i]) next
  g <- g + 1L
  stack <- i
  while (length(stack) > 0) {
    cur <- stack[1]
    stack <- stack[-1]
    if (visited[cur]) next
    visited[cur] <- TRUE
    group_id[cur] <- g
    neighbours <- cross_valid[[cur]]
    stack <- c(stack, neighbours[!visited[neighbours]])
  }
}
test_valid$group <- group_id
test_valid$is_overlap <- lengths(cross_valid) > 1

cat(sprintf("Overlap groups: %d | Groups with overlaps: %d\n",
            max(group_id),
            length(unique(group_id[test_valid$is_overlap]))))

# --- Step 5: Three methods (test area) ---
method_a_mean <- mean(test_valid$change_m, na.rm = TRUE)
method_a_sd <- sd(test_valid$change_m, na.rm = TRUE)
method_a_n <- nrow(test_valid)

clean_only <- test_valid[!test_valid$is_overlap, ]
method_b_mean <- mean(clean_only$change_m, na.rm = TRUE)
method_b_sd <- sd(clean_only$change_m, na.rm = TRUE)
method_b_n <- nrow(clean_only)

group_means <- tapply(test_valid$change_m, test_valid$group, mean, na.rm = TRUE)
method_c_mean <- mean(group_means, na.rm = TRUE)
method_c_sd <- sd(group_means, na.rm = TRUE)
method_c_n <- length(group_means)

cat("\n===== TEST AREA RESULTS (signed, DEM-corrected) =====\n")
cat(sprintf("Total valid transects: %d | Overlapping: %d | Clean: %d\n",
            nrow(test_valid), sum(test_valid$is_overlap), sum(!test_valid$is_overlap)))
cat(sprintf("  A (no correction):    mean = %+.2f m | SD = %.2f | n = %d\n",
            method_a_mean, method_a_sd, method_a_n))
cat(sprintf("  B (remove overlaps):  mean = %+.2f m | SD = %.2f | n = %d\n",
            method_b_mean, method_b_sd, method_b_n))
cat(sprintf("  C (group & average):  mean = %+.2f m | SD = %.2f | n = %d\n",
            method_c_mean, method_c_sd, method_c_n))
cat(sprintf("  Bias A vs B: %+.2f m | Bias C vs B: %+.2f m\n",
            method_a_mean - method_b_mean, method_c_mean - method_b_mean))
cat("=====================================================\n")

# --- Step 6: Test area plots ---
test_valid$method_label <- ifelse(test_valid$is_overlap, "Overlapping", "Non-overlapping")

ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey92", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = smr_solway_bng, fill = "aliceblue", colour = "steelblue", linewidth = 0.3) +
  geom_sf(data = mhws_1970_filt, colour = "orange", linewidth = 0.4) +
  geom_sf(data = mhws_modern_filt, colour = "black", linewidth = 0.5) +
  geom_sf(data = test_valid[!test_valid$is_overlap, ],
          aes(colour = change_m), linewidth = 0.6) +
  geom_sf(data = test_valid[test_valid$is_overlap, ],
          aes(colour = change_m), linewidth = 0.6, linetype = "dashed") +
  scale_colour_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0,
                         name = "Change (m)") +
  coord_sf(xlim = c(st_bbox(test_area)["xmin"], st_bbox(test_area)["xmax"]),
           ylim = c(st_bbox(test_area)["ymin"], st_bbox(test_area)["ymax"]),
           crs = 27700) +
  labs(title = "TEST: Signed Shoreline Change — Modern vs 1970",
       subtitle = "Solid = non-overlapping, Dashed = overlapping | Red = erosion, Blue = accretion") +
  map_theme

# Barplot: overlapping vs non-overlapping
overlap_summary <- data.frame(
  Type = c("Non-overlapping", "Overlapping"),
  Mean = c(mean(test_valid$change_m[!test_valid$is_overlap], na.rm = TRUE),
           mean(test_valid$change_m[test_valid$is_overlap], na.rm = TRUE)),
  SD = c(sd(test_valid$change_m[!test_valid$is_overlap], na.rm = TRUE),
         sd(test_valid$change_m[test_valid$is_overlap], na.rm = TRUE)),
  N = c(sum(!test_valid$is_overlap), sum(test_valid$is_overlap))
)

ggplot(overlap_summary, aes(x = Type, y = Mean, fill = Type)) +
  geom_col(width = 0.5) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.15) +
  geom_text(aes(label = sprintf("n = %d", N)), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("Non-overlapping" = "steelblue", "Overlapping" = "tomato")) +
  labs(title = "TEST: Overlapping vs Non-overlapping Transects",
       subtitle = "Signed change (+ accretion, - erosion) | MHWS 1970 to Modern",
       y = "Change (m)", x = NULL) +
  theme_minimal() +
  theme(legend.position = "none")

# Method comparison barplot
method_df <- data.frame(
  Method = c("A: No correction", "B: Remove overlaps", "C: Group & average"),
  Mean_change = c(method_a_mean, method_b_mean, method_c_mean),
  SD = c(method_a_sd, method_b_sd, method_c_sd),
  N = c(method_a_n, method_b_n, method_c_n)
)

ggplot(method_df, aes(x = Method, y = Mean_change, fill = Method)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = Mean_change - SD, ymax = Mean_change + SD), width = 0.2) +
  geom_hline(yintercept = method_b_mean, linetype = "dashed", colour = "grey40") +
  geom_text(aes(label = sprintf("n = %d", N)), vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a")) +
  labs(title = "TEST: Method Comparison — Mean Signed Shoreline Change",
       subtitle = "+ accretion, - erosion | Dashed = reference (non-overlapping only)",
       y = "Mean shoreline change (m)", x = NULL) +
  theme_minimal() +
  theme(legend.position = "none")

# ============================================================
# FULL AOI: Shoreline change analysis
# ============================================================

cat("\n===== FULL AOI SHORELINE CHANGE =====\n")
cat("Calculating signed shoreline change (full AOI)...\n")
full_change <- calc_shoreline_change(transects_sf, mhws_modern_filt, mhws_1970_filt,
                                     sample_pts_sf, cop30_bng)
transects_sf$change_m <- full_change$change_m
transects_sf$has_change <- !is.na(full_change$change_m)
cat(sprintf("Valid change measurements: %d of %d transects\n",
            sum(transects_sf$has_change), nrow(transects_sf)))

# Overlap groups (full AOI)
full_valid <- transects_sf[transects_sf$has_change, ]
full_pts_valid <- sample_pts_sf[transects_sf$has_change, ]
cross_full <- st_intersects(full_valid, full_valid)
full_valid$is_overlap <- lengths(cross_full) > 1

visited <- rep(FALSE, nrow(full_valid))
group_id <- rep(NA_integer_, nrow(full_valid))
g <- 0L
for (i in seq_len(nrow(full_valid))) {
  if (visited[i]) next
  g <- g + 1L
  stack <- i
  while (length(stack) > 0) {
    cur <- stack[1]
    stack <- stack[-1]
    if (visited[cur]) next
    visited[cur] <- TRUE
    group_id[cur] <- g
    neighbours <- cross_full[[cur]]
    stack <- c(stack, neighbours[!visited[neighbours]])
  }
}
full_valid$group <- group_id

cat(sprintf("Full AOI: %d valid transects | %d overlapping | %d groups\n",
            nrow(full_valid), sum(full_valid$is_overlap), max(group_id)))

# Method A: no correction
full_a <- full_valid
full_a$method <- "A: No correction"

# Method B: remove overlapping transects
full_b <- full_valid[!full_valid$is_overlap, ]
full_b$method <- "B: Overlaps removed"

# Method C: group & average — assign group mean to each transect
group_means_full <- tapply(full_valid$change_m, full_valid$group, mean, na.rm = TRUE)
full_c <- full_valid
full_c$change_m <- group_means_full[as.character(full_c$group)]
full_c$method <- "C: Group corrected"

# --- Build thick shoreline ribbons coloured by change ---
# Buffer sample points into short segments along the MHWS, coloured by change value
build_ribbon <- function(valid_sf, pts_sf, buffer_dist = 100) {
  pts_with_change <- st_sf(
    change_m = valid_sf$change_m,
    geometry = pts_sf$geometry[match(seq_len(nrow(valid_sf)),
                                     seq_len(nrow(pts_sf)))],
    crs = 27700
  )
  pts_with_change$geometry <- st_buffer(pts_with_change$geometry, buffer_dist)
  pts_with_change
}

ribbon_a <- build_ribbon(full_a, full_pts_valid)
ribbon_b <- build_ribbon(full_b, full_pts_valid[!full_valid$is_overlap, ])
ribbon_c <- build_ribbon(full_c, full_pts_valid)

# Interactive mapview: all three methods on one map
# Symmetric breaks centred on zero so grey = 0
change_max <- max(abs(full_valid$change_m), na.rm = TRUE)
change_breaks <- seq(-change_max, change_max, length.out = 21)
change_pal <- colorRampPalette(c("red", "grey90", "blue"))

mapview(ribbon_a, zcol = "change_m", layer.name = "A: No correction",
        col.regions = change_pal, at = change_breaks, alpha.regions = 0.6) +
  mapview(ribbon_b, zcol = "change_m", layer.name = "B: Overlaps removed",
          col.regions = change_pal, at = change_breaks, alpha.regions = 0.6) +
  mapview(ribbon_c, zcol = "change_m", layer.name = "C: Group corrected",
          col.regions = change_pal, at = change_breaks, alpha.regions = 0.6) +
  mapview(mhws_modern_filt, color = "black", lwd = 1, layer.name = "MHWS Modern") +
  mapview(mhws_1970_filt, color = "orange", lwd = 1, layer.name = "MHWS 1970")

# Full AOI summary stats
full_a_mean <- mean(full_a$change_m, na.rm = TRUE)
full_b_mean <- mean(full_b$change_m, na.rm = TRUE)
full_c_groups <- tapply(full_c$change_m, full_valid$group, mean, na.rm = TRUE)
full_c_mean <- mean(full_c_groups, na.rm = TRUE)

cat(sprintf("\n===== FULL AOI RESULTS =====\n"))
cat(sprintf("  A (no correction):   mean = %+.2f m | SD = %.2f | n = %d\n",
            full_a_mean, sd(full_a$change_m, na.rm = TRUE), nrow(full_a)))
cat(sprintf("  B (remove overlaps): mean = %+.2f m | SD = %.2f | n = %d\n",
            full_b_mean, sd(full_b$change_m, na.rm = TRUE), nrow(full_b)))
cat(sprintf("  C (group & average): mean = %+.2f m | SD = %.2f | n = %d groups\n",
            full_c_mean, sd(full_c_groups, na.rm = TRUE), length(full_c_groups)))
cat(sprintf("  Bias A vs B: %+.2f m | Bias C vs B: %+.2f m\n",
            full_a_mean - full_b_mean, full_c_mean - full_b_mean))
cat("============================\n")

# Full AOI barplots
full_overlap_summary <- data.frame(
  Type = c("Non-overlapping", "Overlapping"),
  Mean = c(mean(full_valid$change_m[!full_valid$is_overlap], na.rm = TRUE),
           mean(full_valid$change_m[full_valid$is_overlap], na.rm = TRUE)),
  SD = c(sd(full_valid$change_m[!full_valid$is_overlap], na.rm = TRUE),
         sd(full_valid$change_m[full_valid$is_overlap], na.rm = TRUE)),
  N = c(sum(!full_valid$is_overlap), sum(full_valid$is_overlap))
)

ggplot(full_overlap_summary, aes(x = Type, y = Mean, fill = Type)) +
  geom_col(width = 0.5) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.15) +
  geom_text(aes(label = sprintf("n = %d", N)), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("Non-overlapping" = "steelblue", "Overlapping" = "tomato")) +
  labs(title = "Full AOI: Overlapping vs Non-overlapping Transects",
       subtitle = "Signed change (+ accretion, - erosion) | MHWS 1970 to Modern",
       y = "Change (m)", x = NULL) +
  theme_minimal() +
  theme(legend.position = "none")

full_method_df <- data.frame(
  Method = c("A: No correction", "B: Remove overlaps", "C: Group & average"),
  Mean_change = c(full_a_mean, full_b_mean, full_c_mean),
  SD = c(sd(full_a$change_m, na.rm = TRUE),
         sd(full_b$change_m, na.rm = TRUE),
         sd(full_c_groups, na.rm = TRUE)),
  N = c(nrow(full_a), nrow(full_b), length(full_c_groups))
)

ggplot(full_method_df, aes(x = Method, y = Mean_change, fill = Method)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = Mean_change - SD, ymax = Mean_change + SD), width = 0.2) +
  geom_hline(yintercept = full_b_mean, linetype = "dashed", colour = "grey40") +
  geom_text(aes(label = sprintf("n = %d", N)), vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a")) +
  labs(title = "Full AOI: Method Comparison — Mean Signed Shoreline Change",
       subtitle = "+ accretion, - erosion | Dashed = reference (non-overlapping only)",
       y = "Mean shoreline change (m)", x = NULL) +
  theme_minimal() +
  theme(legend.position = "none")

# ============================================================
# Foreshore slope gradient (dissipative vs reflective)
# ============================================================

cat("\n===== FORESHORE SLOPE GRADIENT =====\n")

# For each sandy foreshore polygon:
# 1. Extract DEM pixels within the polygon
# 2. Determine seaward direction (lowest elevation edge)
# 3. Calculate slope as elevation range / horizontal distance (landward to seaward)

sandy_foreshore$slope_deg <- NA_real_
sandy_foreshore$elev_range_m <- NA_real_
sandy_foreshore$horiz_dist_m <- NA_real_
sandy_foreshore$mean_elev_m <- NA_real_

for (i in seq_len(nrow(sandy_foreshore))) {
  poly_i <- sandy_foreshore[i, ]

  dem_crop <- terra::crop(cop30_bng, vect(poly_i), mask = TRUE)
  vals <- values(dem_crop, na.rm = TRUE)
  if (length(vals) < 2) next

  cells_xy <- xyFromCell(dem_crop, which(!is.na(values(dem_crop))))
  cells_elev <- vals[!is.na(vals)]

  if (nrow(cells_xy) < 2) next

  elev_max <- max(cells_elev)
  elev_min <- min(cells_elev)
  elev_range <- elev_max - elev_min

  if (elev_range < 0.1) {
    sandy_foreshore$slope_deg[i] <- 0
    sandy_foreshore$elev_range_m[i] <- elev_range
    sandy_foreshore$horiz_dist_m[i] <- NA
    sandy_foreshore$mean_elev_m[i] <- mean(cells_elev)
    next
  }

  # Landward point: highest elevation pixel
  land_idx <- which.max(cells_elev)
  land_xy <- cells_xy[land_idx, ]

  # Seaward point: lowest elevation pixel
  sea_idx <- which.min(cells_elev)
  sea_xy <- cells_xy[sea_idx, ]

  horiz_dist <- sqrt((land_xy[1] - sea_xy[1])^2 + (land_xy[2] - sea_xy[2])^2)

  if (horiz_dist < 1) {
    sandy_foreshore$slope_deg[i] <- NA
    next
  }

  slope_rad <- atan2(elev_range, horiz_dist)
  sandy_foreshore$slope_deg[i] <- slope_rad * 180 / pi
  sandy_foreshore$elev_range_m[i] <- elev_range
  sandy_foreshore$horiz_dist_m[i] <- horiz_dist
  sandy_foreshore$mean_elev_m[i] <- mean(cells_elev)
}

valid_slopes <- sum(!is.na(sandy_foreshore$slope_deg))
cat(sprintf("Foreshore slopes calculated: %d of %d polygons\n",
            valid_slopes, nrow(sandy_foreshore)))
cat(sprintf("  Slope range: %.2f° to %.2f°\n",
            min(sandy_foreshore$slope_deg, na.rm = TRUE),
            max(sandy_foreshore$slope_deg, na.rm = TRUE)))
cat(sprintf("  Median slope: %.2f° | Mean slope: %.2f°\n",
            median(sandy_foreshore$slope_deg, na.rm = TRUE),
            mean(sandy_foreshore$slope_deg, na.rm = TRUE)))

# Classify: dissipative (< ~2°), intermediate (2-6°), reflective (> ~6°)
sandy_foreshore$beach_type <- cut(sandy_foreshore$slope_deg,
                                   breaks = c(-Inf, 2, 6, Inf),
                                   labels = c("Dissipative", "Intermediate", "Reflective"))
cat(sprintf("  Dissipative (< 2°): %d | Intermediate (2-6°): %d | Reflective (> 6°): %d\n",
            sum(sandy_foreshore$beach_type == "Dissipative", na.rm = TRUE),
            sum(sandy_foreshore$beach_type == "Intermediate", na.rm = TRUE),
            sum(sandy_foreshore$beach_type == "Reflective", na.rm = TRUE)))

# Mapview: slope gradient
mapview(sandy_foreshore[!is.na(sandy_foreshore$slope_deg), ],
        zcol = "slope_deg", layer.name = "Foreshore slope (°)",
        col.regions = colorRampPalette(c("blue", "yellow", "red")),
        alpha.regions = 0.7) +
  mapview(sandy_foreshore[!is.na(sandy_foreshore$beach_type), ],
          zcol = "beach_type", layer.name = "Beach type",
          col.regions = c("Dissipative" = "blue", "Intermediate" = "yellow",
                          "Reflective" = "red"),
          alpha.regions = 0.5) +
  mapview(mhws_modern_filt, color = "black", lwd = 1, layer.name = "MHWS Modern")

# --- Rule 3: Filter MHWS to dissipative beaches only ---
dissipative_fs <- sandy_foreshore[sandy_foreshore$beach_type == "Dissipative" &
                                    !is.na(sandy_foreshore$beach_type), ]
cat(sprintf("\nDissipative foreshore polygons: %d of %d sandy foreshore\n",
            nrow(dissipative_fs), nrow(sandy_foreshore)))

dissipative_union <- st_union(dissipative_fs) %>% st_buffer(200) %>% st_make_valid()

filter_mhws_dissipative <- function(mhws, diss_buf, saltmarsh, label) {
  hits <- lengths(st_intersects(mhws, diss_buf)) > 0
  within_diss <- mhws[hits, ]
  cat(sprintf("  %s: %d -> %d after dissipative foreshore filter\n",
              label, nrow(mhws), nrow(within_diss)))
  minus_sm <- st_difference(st_union(within_diss), saltmarsh) %>%
    st_make_valid() %>%
    st_cast("MULTILINESTRING") %>%
    st_cast("LINESTRING")
  result <- st_sf(geometry = minus_sm, crs = 27700)
  cat(sprintf("  %s: %d features after removing saltmarsh\n", label, nrow(result)))
  result
}

mhws_1890_diss <- filter_mhws_dissipative(mhws_1890_clip, dissipative_union, sm_union, "MHWS 1890")
mhws_1970_diss <- filter_mhws_dissipative(mhws_1970_clip, dissipative_union, sm_union, "MHWS 1970")
mhws_modern_diss <- filter_mhws_dissipative(mhws_modern_clip, dissipative_union, sm_union, "MHWS Modern")

cat(sprintf("\nFiltered to dissipative beaches:\n"))
cat(sprintf("  MHWS 1890: %d features\n", nrow(mhws_1890_diss)))
cat(sprintf("  MHWS 1970: %d features\n", nrow(mhws_1970_diss)))
cat(sprintf("  MHWS Modern: %d features\n", nrow(mhws_modern_diss)))

# Mapview: dissipative-filtered MHWS
mapview(dissipative_fs, col.regions = "lightblue", layer.name = "Dissipative foreshore") +
  mapview(mhws_1890_diss, color = "brown", lwd = 1, layer.name = "MHWS 1890 (dissipative)") +
  mapview(mhws_1970_diss, color = "orange", lwd = 1, layer.name = "MHWS 1970 (dissipative)") +
  mapview(mhws_modern_diss, color = "black", lwd = 2, layer.name = "MHWS Modern (dissipative)") +
  mapview(sm_harmonised_sf, col.regions = "green", alpha.regions = 0.3, layer.name = "Saltmarsh (excluded)")

# ============================================================
# Calculating Accommodation Space
# ============================================================

cat("\n===== CALCULATING ACCOMMODATION SPACE =====\n")

# --- 1. Coastal Defences — seawalls only (Dynamic Coast) ---
cat("\n-- Coastal Defences (seawalls) --\n")

defence_line_path <- "data/raw/habitats/Sand dunes/DYNAMICCOAST2_SCOTLAND_SHP_27700/shp/DC2_Defences_Line.shp"
stopifnot("DC2_Defences_Line not found" = file.exists(defence_line_path))
defence_line <- st_read(defence_line_path, quiet = TRUE) %>% st_make_valid()
cat(sprintf("  Defences (line) loaded: %d features | CRS: %s\n",
            nrow(defence_line), st_crs(defence_line)$input))
cat(sprintf("  Columns: %s\n", paste(names(defence_line), collapse = ", ")))
if (st_crs(defence_line)$epsg != 27700) defence_line <- st_transform(defence_line, 27700)
defence_line_clip <- st_intersection(defence_line, scamp_aoi_sf) %>% st_make_valid()
cat(sprintf("  Defences (line) in AOI: %d features\n", nrow(defence_line_clip)))

# Print all categorical columns to inspect defence types
geom_col_dl <- attr(defence_line_clip, "sf_column")
for (cc in setdiff(names(defence_line_clip), geom_col_dl)) {
  if (is.character(defence_line_clip[[cc]]) || is.factor(defence_line_clip[[cc]])) {
    vals <- sort(unique(defence_line_clip[[cc]]))
    if (length(vals) <= 30) {
      cat(sprintf("\n  %s (%d unique):\n", cc, length(vals)))
      for (v in vals) {
        n <- sum(defence_line_clip[[cc]] == v, na.rm = TRUE)
        cat(sprintf("    - %s  [%d]\n", v, n))
      }
    }
  }
}

# Filter for seawalls only
seawall_col <- names(defence_line_clip)[sapply(names(defence_line_clip), function(cc) {
  if (cc == geom_col_dl) return(FALSE)
  any(grepl("seawall|sea wall", defence_line_clip[[cc]], ignore.case = TRUE))
})]
cat(sprintf("\n  Column(s) containing 'seawall': %s\n", paste(seawall_col, collapse = ", ")))

if (length(seawall_col) > 0) {
  seawalls_clip <- defence_line_clip[grepl("seawall|sea wall", defence_line_clip[[seawall_col[1]]],
                                            ignore.case = TRUE), ]
} else {
  cat("  WARNING: No 'seawall' found in any column — keeping all defences\n")
  seawalls_clip <- defence_line_clip
}
cat(sprintf("  Seawalls in AOI: %d features\n", nrow(seawalls_clip)))

# --- 2. Road Networks (OS OpenRoads — RoadLink, tiles NX/NW/NY) ---
cat("\n-- Road Networks (OS OpenRoads) --\n")

road_tiles <- c("NX", "NW", "NY")
road_list <- lapply(road_tiles, function(tile) {
  fpath <- sprintf("data/raw/physical/OpenRoads/data/%s_RoadLink.shp", tile)
  stopifnot(file.exists(fpath))
  sf_obj <- st_read(fpath, quiet = TRUE) %>% st_make_valid()
  cat(sprintf("  %s_RoadLink: %d features | CRS: %s\n",
              tile, nrow(sf_obj), st_crs(sf_obj)$input))
  if (st_crs(sf_obj)$epsg != 27700) sf_obj <- st_transform(sf_obj, 27700)
  sf_obj
})
roads_all <- do.call(rbind, road_list)
cat(sprintf("  Combined road features: %d\n", nrow(roads_all)))
cat(sprintf("  Columns: %s\n", paste(names(roads_all), collapse = ", ")))

roads_clip <- st_intersection(roads_all, scamp_aoi_sf) %>% st_make_valid()
cat(sprintf("  Roads in AOI: %d features\n", nrow(roads_clip)))

geom_col_rd <- attr(roads_clip, "sf_column")
for (cc in setdiff(names(roads_clip), geom_col_rd)) {
  if (is.character(roads_clip[[cc]]) || is.factor(roads_clip[[cc]])) {
    vals <- sort(unique(roads_clip[[cc]]))
    if (length(vals) <= 30) {
      cat(sprintf("\n  %s (%d unique):\n", cc, length(vals)))
      for (v in vals) {
        n <- sum(roads_clip[[cc]] == v, na.rm = TRUE)
        cat(sprintf("    - %s  [%d]\n", v, n))
      }
    }
  }
}

# --- 3. Buildings (OS OpenMap Local — tiles NX/NW/NY) ---
cat("\n-- Buildings (OS OpenMap Local) --\n")

building_paths <- c(
  "data/raw/physical/opmplc_essh_nx/OS OpenMap Local (ESRI Shape File) NX/data/NX_Building.shp",
  "data/raw/physical/opmplc_essh_nw/OS OpenMap Local (ESRI Shape File) NW/data/NW_Building.shp",
  "data/raw/physical/opmplc_essh_ny/OS OpenMap Local (ESRI Shape File) NY/data/NY_Building.shp"
)

building_list <- lapply(seq_along(building_paths), function(j) {
  fpath <- building_paths[j]
  tile <- c("NX", "NW", "NY")[j]
  stopifnot(file.exists(fpath))
  sf_obj <- st_read(fpath, quiet = TRUE) %>% st_make_valid()
  cat(sprintf("  %s_Building: %d features | CRS: %s\n",
              tile, nrow(sf_obj), st_crs(sf_obj)$input))
  if (st_crs(sf_obj)$epsg != 27700) sf_obj <- st_transform(sf_obj, 27700)
  sf_obj
})
buildings_all <- do.call(rbind, building_list)
cat(sprintf("  Combined building features: %d\n", nrow(buildings_all)))
cat(sprintf("  Columns: %s\n", paste(names(buildings_all), collapse = ", ")))

buildings_clip <- st_intersection(buildings_all, scamp_aoi_sf) %>% st_make_valid()
cat(sprintf("  Buildings in AOI: %d features\n", nrow(buildings_clip)))

geom_col_bl <- attr(buildings_clip, "sf_column")
for (cc in setdiff(names(buildings_clip), geom_col_bl)) {
  if (is.character(buildings_clip[[cc]]) || is.factor(buildings_clip[[cc]])) {
    vals <- sort(unique(buildings_clip[[cc]]))
    if (length(vals) <= 30) {
      cat(sprintf("\n  %s (%d unique):\n", cc, length(vals)))
      for (v in vals) {
        n <- sum(buildings_clip[[cc]] == v, na.rm = TRUE)
        cat(sprintf("    - %s  [%d]\n", v, n))
      }
    }
  }
}

cat("\n===== Accommodation space data loaded — inspect summaries above =====\n")

# --- Land cover: inspect HabMoS full dataset structure ---
cat("\n-- HabMoS data structure (full AOI clip) --\n")
cat(sprintf("  Features: %d | CRS: %s\n", nrow(habmos_clip), st_crs(habmos_clip)$input))
cat(sprintf("  Columns: %s\n", paste(names(habmos_clip), collapse = ", ")))

print(unique(habmos_clip$HABITAT_NAME))

# ============================================================
# Accommodation Space: distance from filtered MHWS to nearest constraint
# ============================================================

cat("\n===== ACCOMMODATION SPACE CALCULATION =====\n")

# Sample points along the dissipative-filtered modern MHWS at 50m spacing
mhws_diss_merged <- st_line_merge(st_union(st_cast(st_geometry(mhws_modern_diss), "LINESTRING")))
diss_parts <- st_cast(mhws_diss_merged, "LINESTRING")
diss_coords <- st_coordinates(mhws_diss_merged)

accom_pts_list <- lapply(unique(diss_coords[, "L1"]), function(lid) {
  cc <- diss_coords[diss_coords[, "L1"] == lid, 1:2, drop = FALSE]
  seg_len <- sqrt(diff(cc[, 1])^2 + diff(cc[, 2])^2)
  cum_len <- c(0, cumsum(seg_len))
  total_len <- cum_len[length(cum_len)]
  if (total_len < 50) return(NULL)
  sample_dists <- seq(0, total_len, by = 50)
  seg_idx <- pmin(findInterval(sample_dists, cum_len, rightmost.closed = TRUE), length(seg_len))
  frac <- (sample_dists - cum_len[seg_idx]) / seg_len[seg_idx]
  frac[is.nan(frac)] <- 0
  px <- cc[seg_idx, 1] + frac * diff(cc[, 1])[seg_idx]
  py <- cc[seg_idx, 2] + frac * diff(cc[, 2])[seg_idx]
  cbind(px, py)
})
accom_mat <- do.call(rbind, Filter(Negate(is.null), accom_pts_list))
accom_pts <- st_sf(
  id = seq_len(nrow(accom_mat)),
  geometry = st_sfc(lapply(seq_len(nrow(accom_mat)), function(i) {
    st_point(accom_mat[i, ])
  }), crs = 27700)
)
cat(sprintf("Sampled %d points along dissipative-filtered MHWS\n", nrow(accom_pts)))

# Distance to nearest seawall
dist_seawall <- rep(NA_real_, nrow(accom_pts))
if (nrow(seawalls_clip) > 0) {
  dist_seawall <- as.numeric(st_distance(accom_pts, st_union(seawalls_clip)))
}
cat(sprintf("  Seawall distances: %d valid\n", sum(!is.na(dist_seawall))))

# Distance to nearest road
dist_road <- as.numeric(st_distance(accom_pts, st_union(roads_clip)))
cat(sprintf("  Road distances: %d valid\n", sum(!is.na(dist_road))))

# Distance to nearest building
dist_building <- as.numeric(st_distance(accom_pts, st_union(buildings_clip)))
cat(sprintf("  Building distances: %d valid\n", sum(!is.na(dist_building))))

# Accommodation space = distance to nearest constraint
accom_pts$dist_seawall <- dist_seawall
accom_pts$dist_road <- dist_road
accom_pts$dist_building <- dist_building
accom_pts$accom_space_m <- pmin(dist_seawall, dist_road, dist_building, na.rm = TRUE)
accom_pts$nearest_constraint <- c("Seawall", "Road", "Building")[
  apply(cbind(dist_seawall, dist_road, dist_building), 1, which.min)
]

cat(sprintf("\nAccommodation space summary:\n"))
cat(sprintf("  Mean: %.0f m | Median: %.0f m | Range: %.0f–%.0f m\n",
            mean(accom_pts$accom_space_m, na.rm = TRUE),
            median(accom_pts$accom_space_m, na.rm = TRUE),
            min(accom_pts$accom_space_m, na.rm = TRUE),
            max(accom_pts$accom_space_m, na.rm = TRUE)))
cat(sprintf("  Nearest constraint breakdown:\n"))
for (ctype in c("Seawall", "Road", "Building")) {
  n <- sum(accom_pts$nearest_constraint == ctype, na.rm = TRUE)
  med <- median(accom_pts$accom_space_m[accom_pts$nearest_constraint == ctype], na.rm = TRUE)
  cat(sprintf("    %s: %d points (median %.0f m)\n", ctype, n, med))
}

# Mapview: accommodation space
accom_pal <- colorRampPalette(c("red", "yellow", "green"))
mapview(accom_pts, zcol = "accom_space_m", layer.name = "Accommodation space (m)",
        col.regions = accom_pal, cex = 3, alpha.regions = 0.8) +
  mapview(accom_pts, zcol = "nearest_constraint", layer.name = "Nearest constraint",
          col.regions = c("Seawall" = "black", "Road" = "grey50", "Building" = "brown"),
          cex = 3, alpha.regions = 0.7) +
  mapview(seawalls_clip, color = "black", lwd = 2, layer.name = "Seawalls") +
  mapview(roads_clip, color = "grey50", lwd = 0.5, layer.name = "Roads") +
  mapview(mhws_modern_diss, color = "blue", lwd = 1, layer.name = "MHWS Modern (dissipative)")

# Filter: minimum 100m accommodation space
accom_viable <- accom_pts[accom_pts$accom_space_m >= 100, ]
cat(sprintf("\nViable accommodation (>= 100m): %d of %d points (%.0f%%)\n",
            nrow(accom_viable), nrow(accom_pts),
            100 * nrow(accom_viable) / nrow(accom_pts)))

# ============================================================
# Sea Level Rise — RCP4.5 projection to 2050
# ============================================================

cat("\n===== SEA LEVEL RISE (RCP4.5) =====\n")

slr_path <- "data/raw/physical/seaLevelAnom_marine-sim_rcp45_ann_2007-2100.nc"
stopifnot("Sea level NetCDF not found" = file.exists(slr_path))

slr_nc <- rast(slr_path)
cat(sprintf("  SLR NetCDF loaded: %d layers | %d x %d | CRS: %s\n",
            nlyr(slr_nc), ncol(slr_nc), nrow(slr_nc), crs(slr_nc, describe = TRUE)$code))
cat(sprintf("  Layer names (first 10): %s\n",
            paste(head(names(slr_nc), 10), collapse = ", ")))
cat(sprintf("  Time info: %s\n", paste(head(time(slr_nc), 10), collapse = ", ")))

# 846 layers = 9 percentiles x 94 years (2007-2100)
# Identify all percentile bands
slr_names <- names(slr_nc)
percentiles <- unique(gsub("_\\d+$", "", slr_names))
cat(sprintf("  Percentile bands: %s\n", paste(percentiles, collapse = ", ")))
n_percentiles <- length(percentiles)
n_years <- nlyr(slr_nc) / n_percentiles
cat(sprintf("  %d percentiles x %.0f years\n", n_percentiles, n_years))

# Diagnose NetCDF structure: print ALL unique prefixes and sample layer names
cat("  First 20 layer names:\n")
cat(paste("   ", head(slr_names, 20), collapse = "\n"), "\n")

# Find the median (percentile=50) band — search for "50" in any format
p50_idx <- grep("50", slr_names)
cat(sprintf("  Layers matching '50': %d indices\n", length(p50_idx)))
if (length(p50_idx) > 0) {
  cat(sprintf("  Sample names: %s\n", paste(head(slr_names[p50_idx], 5), collapse = ", ")))
}

# Strategy: with 9 percentiles x 94 years, layers are likely ordered as
# [p1_y1, p1_y2, ..., p1_y94, p2_y1, ..., p9_y94]
# OR [p1_y1, p2_y1, ..., p9_y1, p1_y2, ..., p9_y94]
# Print values from several layers to find non-zero ones
cat("  Checking values in every 94th layer (one per percentile band):\n")
for (i in 1:min(9, nlyr(slr_nc))) {
  lyr_idx <- (i - 1) * 94 + 44  # year 44 = 2050 for each percentile block
  if (lyr_idx <= nlyr(slr_nc)) {
    vals <- values(slr_nc[[lyr_idx]])
    cat(sprintf("    Layer %d (%s): min=%.4f max=%.4f\n",
                lyr_idx, names(slr_nc)[lyr_idx],
                min(vals, na.rm = TRUE), max(vals, na.rm = TRUE)))
  }
}
cat("  Checking layers 44, 45, 46, 47, 48 directly:\n")
for (i in 44:48) {
  vals <- values(slr_nc[[i]])
  cat(sprintf("    Layer %d (%s): min=%.4f max=%.4f\n",
              i, names(slr_nc)[i],
              min(vals, na.rm = TRUE), max(vals, na.rm = TRUE)))
}

# For now, pick the layer with the largest max value near index 44*5 (median band, 2050)
# We'll refine once we see the diagnostic output
best_idx <- 1
best_max <- 0
for (i in seq_len(nlyr(slr_nc))) {
  mx <- minmax(slr_nc[[i]])[2]
  if (!is.na(mx) && mx > best_max) {
    best_max <- mx
    best_idx <- i
  }
}
cat(sprintf("  Layer with highest max value: %d (%s) max=%.4f\n",
            best_idx, names(slr_nc)[best_idx], best_max))

# Select median percentile for 2050 — use p50 match if valid, fallback to brute-force
p50_prefix <- unique(gsub("_\\d+$", "", slr_names))
p50_match <- p50_prefix[grep("50", p50_prefix)]
if (length(p50_match) > 0) {
  p50_layers <- grep(paste0("^", fixed = FALSE,
                            gsub("([\\(\\)\\[\\]])", "\\\\\\1", p50_match[1])), slr_names)
  year_idx <- 2050 - 2007 + 1
  if (length(p50_layers) >= year_idx) {
    slr_2050_idx <- p50_layers[year_idx]
  } else {
    cat("  WARNING: year index out of range, using layer with highest max\n")
    slr_2050_idx <- best_idx
  }
} else {
  cat("  WARNING: No percentile=50 found, using layer with highest max\n")
  slr_2050_idx <- best_idx
}
slr_2050 <- slr_nc[[slr_2050_idx]]
cat(sprintf("  FINAL selected layer %d: %s\n", slr_2050_idx, names(slr_2050)))
cat(sprintf("  Value range: %.4f to %.4f\n",
            minmax(slr_2050)[1], minmax(slr_2050)[2]))

# Crop to AOI vicinity and reproject to BNG
aoi_slr <- vect(st_transform(scamp_aoi_sf, crs(slr_2050)))
slr_2050_crop <- crop(slr_2050, ext(aoi_slr) + 1)
slr_2050_bng <- project(slr_2050_crop, "EPSG:27700", method = "bilinear")
cat(sprintf("  SLR 2050 reprojected to BNG: %d x %d | res: %.0f m\n",
            ncol(slr_2050_bng), nrow(slr_2050_bng), res(slr_2050_bng)[1]))

# Disaggregate to ~1km and mask to near-coast low-elevation areas
slr_fine <- disagg(slr_2050_bng, fact = 10, method = "bilinear")
cat(sprintf("  SLR disaggregated: %d x %d | res: %.0f m\n",
            ncol(slr_fine), nrow(slr_fine), res(slr_fine)[1]))

# Resample DEM to match SLR grid, mask to low elevation (< 20m) to keep coastal zone only
dem_match <- resample(cop30_bng, slr_fine, method = "bilinear")
slr_coastal <- mask(slr_fine, ifel(dem_match < 20, 1, NA))
cat(sprintf("  SLR masked to coastal zone (DEM < 20m): %d non-NA pixels\n",
            sum(!is.na(values(slr_coastal)))))

# Convert to polygons for mapview
slr_2050_sf <- as.polygons(slr_coastal) %>% st_as_sf()
names(slr_2050_sf)[1] <- "slr_m"
slr_2050_sf <- slr_2050_sf[!is.na(slr_2050_sf$slr_m), ]

# Mapview: SLR 2050 + accommodation space + existing dunes
mapview(slr_2050_sf, zcol = "slr_m", layer.name = "Sea Level Rise 2050 (m)",
        col.regions = colorRampPalette(c("lightyellow", "orange", "darkred")),
        alpha.regions = 0.5) +
  mapview(accom_viable, zcol = "accom_space_m", layer.name = "Accommodation >= 100m",
          col.regions = colorRampPalette(c("yellow", "green")),
          cex = 3, alpha.regions = 0.8) +
  mapview(habmos_dunes_bng, col.regions = "purple", alpha.regions = 0.5,
          layer.name = "Existing Sand Dunes (HabMoS)") +
  mapview(mhws_modern_diss, color = "black", lwd = 1, layer.name = "MHWS Modern (dissipative)")

# ----- Alternative SLR dataset: Local Authority projections -----
slr_la_path <- "data/raw/physical/sea_level_rise_projections_local_authority_v1_-4803453840407794973.gpkg"
stopifnot("SLR LA GPKG not found" = file.exists(slr_la_path))

cat("\n===== SLR LOCAL AUTHORITY GPKG =====\n")
cat(sprintf("  Layers: %s\n", paste(st_layers(slr_la_path)$name, collapse = ", ")))

slr_la <- st_read(slr_la_path, quiet = TRUE)
cat(sprintf("  Features: %d | Fields: %d\n", nrow(slr_la), ncol(slr_la)))
cat(sprintf("  CRS: %s\n", st_crs(slr_la)$input))
cat("  Column names:\n")
cat(paste("   ", names(slr_la), collapse = "\n"), "\n")
cat("  Summary:\n")
print(summary(slr_la))
cat("  First 5 rows:\n")
print(head(slr_la, 5))

# Unique RCP/percentile combos
cat("  Unique RCPandPercentile values:\n")
print(sort(unique(slr_la$RCPandPercentile)))

# Filter to RCP4.5 median (50th percentile)
slr_la_rcp45 <- slr_la[slr_la$RCPandPercentile == "rcp45_50" & !is.na(slr_la$RCPandPercentile), ]
cat(sprintf("  RCP4.5 median features: %d\n", nrow(slr_la_rcp45)))

# Reproject to BNG and clip to SCAMP AOI
slr_la_bng <- st_transform(slr_la_rcp45, 27700)
slr_la_clip <- st_intersection(slr_la_bng, scamp_aoi_sf)
cat(sprintf("  Clipped to SCAMP AOI: %d features\n", nrow(slr_la_clip)))
cat(sprintf("  Local authorities: %s\n", paste(unique(slr_la_clip$NAME), collapse = ", ")))
cat(sprintf("  SLR 2050s range: %d–%d cm\n",
            min(slr_la_clip$seaLevelAnom_2050s, na.rm = TRUE),
            max(slr_la_clip$seaLevelAnom_2050s, na.rm = TRUE)))

mapview(slr_la_clip, zcol = "seaLevelAnom_2050s", layer.name = "SLR 2050s (cm, RCP4.5 median)",
        col.regions = colorRampPalette(c("lightyellow", "orange", "darkred")),
        alpha.regions = 0.5) +
  mapview(accom_viable, zcol = "accom_space_m", layer.name = "Accommodation >= 100m",
          col.regions = colorRampPalette(c("yellow", "green")),
          cex = 3, alpha.regions = 0.8) +
  mapview(habmos_dunes_bng, col.regions = "purple", alpha.regions = 0.5,
          layer.name = "Existing Sand Dunes (HabMoS)") +
  mapview(mhws_modern_diss, color = "black", lwd = 1, layer.name = "MHWS Modern (dissipative)")

# ============================================================
# 2050 Accommodation Space: reduce by SLR inland shift (18cm)
# ============================================================

cat("\n===== 2050 SLR ACCOMMODATION ADJUSTMENT =====\n")
slr_rise_m <- 0.18

# Calculate local slope at each accommodation point from DEM
# Use terrain() for slope in radians, then inland_shift = slr / tan(slope)
slope_rad <- terrain(cop30_bng, v = "slope", unit = "radians")
pt_slope <- terra::extract(slope_rad, vect(accom_pts))[, 2]

# Inland shift: how far 0.18m rise moves the waterline at each point's slope
# Cap minimum slope at 0.5 degrees to avoid extreme shifts on flat terrain
min_slope_rad <- 0.5 * pi / 180
pt_slope_capped <- pmax(pt_slope, min_slope_rad, na.rm = TRUE)
inland_shift_m <- slr_rise_m / tan(pt_slope_capped)
inland_shift_m[is.na(inland_shift_m)] <- 0

accom_pts$inland_shift_m <- inland_shift_m
accom_pts$accom_2050_m <- pmax(accom_pts$accom_space_m - inland_shift_m, 0)

cat(sprintf("  SLR rise: %.2f m\n", slr_rise_m))
cat(sprintf("  Inland shift: mean=%.1f m, median=%.1f m, max=%.1f m\n",
            mean(inland_shift_m, na.rm = TRUE),
            median(inland_shift_m, na.rm = TRUE),
            max(inland_shift_m, na.rm = TRUE)))
cat(sprintf("\n  Current accommodation: mean=%.0f m, median=%.0f m\n",
            mean(accom_pts$accom_space_m, na.rm = TRUE),
            median(accom_pts$accom_space_m, na.rm = TRUE)))
cat(sprintf("  2050 accommodation:   mean=%.0f m, median=%.0f m\n",
            mean(accom_pts$accom_2050_m, na.rm = TRUE),
            median(accom_pts$accom_2050_m, na.rm = TRUE)))

accom_viable_current <- sum(accom_pts$accom_space_m >= 100, na.rm = TRUE)
accom_viable_2050 <- sum(accom_pts$accom_2050_m >= 100, na.rm = TRUE)
cat(sprintf("\n  Viable sites (>= 100m): current %d → 2050 %d (lost %d, %.1f%%)\n",
            accom_viable_current, accom_viable_2050,
            accom_viable_current - accom_viable_2050,
            100 * (accom_viable_current - accom_viable_2050) / accom_viable_current))

# Mapview: current vs 2050 accommodation
mapview(accom_pts, zcol = "accom_space_m", layer.name = "Accommodation Current (m)",
        col.regions = colorRampPalette(c("red", "yellow", "green")),
        cex = 3, alpha.regions = 0.8) +
  mapview(accom_pts, zcol = "accom_2050_m", layer.name = "Accommodation 2050 (m)",
          col.regions = colorRampPalette(c("red", "yellow", "green")),
          cex = 3, alpha.regions = 0.6) +
  mapview(accom_pts, zcol = "inland_shift_m", layer.name = "Inland Shift (m)",
          col.regions = colorRampPalette(c("lightyellow", "orange", "darkred")),
          cex = 3, alpha.regions = 0.6) +
  mapview(mhws_modern_diss, color = "blue", lwd = 1.5, layer.name = "MHWS Modern") +
  mapview(habmos_dunes_bng, col.regions = "purple", alpha.regions = 0.4,
          layer.name = "Existing Sand Dunes")

# ============================================================
# Ecological filter: proximity to dune systems & indicator species
# ============================================================

cat("\n===== ECOLOGICAL FILTER =====\n")

# Start with 2050-viable accommodation points (>= 100m)
accom_2050_viable <- accom_pts[accom_pts$accom_2050_m >= 100, ]
cat(sprintf("  Starting points (2050 accommodation >= 100m): %d\n", nrow(accom_2050_viable)))

# Rule 1: proximity to existing dune systems (source populations for colonisation)
# Distance from each viable point to nearest existing dune polygon
dunes_union <- st_union(habmos_dunes_bng)
accom_2050_viable$dist_to_dune_m <- as.numeric(st_distance(accom_2050_viable, dunes_union))
cat(sprintf("  Distance to nearest dune system: mean=%.0f m, median=%.0f m, max=%.0f m\n",
            mean(accom_2050_viable$dist_to_dune_m),
            median(accom_2050_viable$dist_to_dune_m),
            max(accom_2050_viable$dist_to_dune_m)))

# Rule 2: observations of dune indicator species outside known dune systems
# Remove indicator records in the sea (DEM elevation > 0 = on land)
ind_elev <- terra::extract(cop30_bng, vect(dune_ind))[, 2]
dune_ind_land <- dune_ind[!is.na(ind_elev) & ind_elev > 0, ]
cat(sprintf("  Indicator records on land (DEM > 0m): %d / %d\n", nrow(dune_ind_land), nrow(dune_ind)))
# Remove records that fall within existing dune polygons
dune_ind_outside <- dune_ind_land[lengths(st_intersects(dune_ind_land, habmos_dunes_bng)) == 0, ]
cat(sprintf("  Dune indicator records total: %d | outside known dunes: %d\n",
            nrow(dune_ind), nrow(dune_ind_outside)))
if (nrow(dune_ind_outside) > 0) {
  cat("  Species outside known dunes:\n")
  print(dune_ind_outside %>% st_drop_geometry() %>% count(latin_name, sort = TRUE))
}

# Combine: point is ecologically viable if within 2km of existing dune system
# OR within 1km of an indicator species observation outside known dunes
near_dune_system <- accom_2050_viable$dist_to_dune_m <= 2000

near_indicator <- rep(FALSE, nrow(accom_2050_viable))
if (nrow(dune_ind_outside) > 0) {
  ind_union <- st_union(dune_ind_outside)
  dist_to_ind <- as.numeric(st_distance(accom_2050_viable, ind_union))
  near_indicator <- dist_to_ind <= 1000
  accom_2050_viable$dist_to_indicator_m <- dist_to_ind
}

accom_2050_viable$eco_viable <- near_dune_system | near_indicator
accom_2050_viable$eco_reason <- ifelse(
  near_dune_system & near_indicator, "Both",
  ifelse(near_dune_system, "Near dune system",
         ifelse(near_indicator, "Near indicator spp", "Too remote")))

cat(sprintf("\n  Ecologically viable: %d / %d (%.0f%%)\n",
            sum(accom_2050_viable$eco_viable),
            nrow(accom_2050_viable),
            100 * sum(accom_2050_viable$eco_viable) / nrow(accom_2050_viable)))
cat("  Breakdown:\n")
for (r in c("Near dune system", "Near indicator spp", "Both", "Too remote")) {
  n <- sum(accom_2050_viable$eco_reason == r)
  if (n > 0) cat(sprintf("    %s: %d\n", r, n))
}

# Assign tiers: Tier 1 = near existing dunes (with or without indicators), Tier 2 = everything else
accom_2050_viable$tier <- ifelse(
  accom_2050_viable$eco_reason %in% c("Both", "Near dune system"), 1, 2)
cat(sprintf("\n  Tier 1 (near existing dunes): %d\n", sum(accom_2050_viable$tier == 1)))
cat(sprintf("  Tier 2 (other viable): %d\n", sum(accom_2050_viable$tier == 2)))

# Mapview: tiered viability
mapview(accom_2050_viable, zcol = "tier", layer.name = "Suitability Tier",
        col.regions = c("1" = "darkgreen", "2" = "orange"),
        cex = 3, alpha.regions = 0.8) +
  mapview(habmos_dunes_bng, col.regions = "purple", alpha.regions = 0.4,
          layer.name = "Existing Dune Systems") +
  mapview(mhws_modern_diss, color = "blue", lwd = 1, layer.name = "MHWS Modern")

# ============================================================
# Build restoration site polygons (MHWS to nearest constraint)
# ============================================================

cat("\n===== BUILDING RESTORATION POLYGONS =====\n")

# For each viable point, build a rectangle from the MHWS inland to the accommodation distance
# Direction: perpendicular to MHWS, landward (determined by DEM elevation)
# Width: 50m (matching sample spacing), centred on each point

# Find landward bearing: vectorised MHWS bearing + DEM check
coords_viable <- st_coordinates(accom_2050_viable)
half_width <- 25
mhws_all_coords <- st_coordinates(st_union(mhws_modern_diss))
n_mhws <- nrow(mhws_all_coords)
n_pts <- nrow(accom_2050_viable)

# Nearest MHWS vertex for each viable point (vectorised distance matrix approach in chunks)
nearest_idx <- integer(n_pts)
chunk_size <- 500
for (ch in seq(1, n_pts, by = chunk_size)) {
  ch_end <- min(ch + chunk_size - 1, n_pts)
  dx_mat <- outer(coords_viable[ch:ch_end, 1], mhws_all_coords[, 1], `-`)
  dy_mat <- outer(coords_viable[ch:ch_end, 2], mhws_all_coords[, 2], `-`)
  nearest_idx[ch:ch_end] <- max.col(-sqrt(dx_mat^2 + dy_mat^2))
}

# Local shore bearing from ±5 vertices around nearest point
idx_before <- pmax(1, nearest_idx - 5)
idx_after <- pmin(n_mhws, nearest_idx + 5)
dx <- mhws_all_coords[idx_after, 1] - mhws_all_coords[idx_before, 1]
dy <- mhws_all_coords[idx_after, 2] - mhws_all_coords[idx_before, 2]
shore_bearing <- atan2(dy, dx)

# Two perpendicular candidates
perp1 <- shore_bearing + pi / 2
perp2 <- shore_bearing - pi / 2

# Sample DEM 100m in each direction (single batch extract)
probe1 <- cbind(coords_viable[, 1] + 100 * cos(perp1), coords_viable[, 2] + 100 * sin(perp1))
probe2 <- cbind(coords_viable[, 1] + 100 * cos(perp2), coords_viable[, 2] + 100 * sin(perp2))
elev1 <- terra::extract(cop30_bng, vect(probe1, crs = "EPSG:27700"))[, 2]
elev2 <- terra::extract(cop30_bng, vect(probe2, crs = "EPSG:27700"))[, 2]
elev1[is.na(elev1)] <- -999
elev2[is.na(elev2)] <- -999

# Pick the perpendicular with higher elevation (landward)
landward_bearing <- ifelse(elev1 >= elev2, perp1, perp2)
cat(sprintf("  Computed landward bearings for %d points\n", n_pts))

# Build rectangles: 50m wide (25m either side), extending landward by accom_2050_m
depths <- accom_2050_viable$accom_2050_m
cat(sprintf("  Depths: %d values | NA: %d | <= 0: %d | range: %.1f–%.1f\n",
            length(depths), sum(is.na(depths)), sum(depths <= 0, na.rm = TRUE),
            min(depths, na.rm = TRUE), max(depths, na.rm = TRUE)))
cat(sprintf("  Bearings: %d values | NA: %d\n",
            length(landward_bearing), sum(is.na(landward_bearing))))

# Vectorised rectangle construction
cx <- coords_viable[, 1]
cy <- coords_viable[, 2]
perp_angle <- landward_bearing + pi / 2

# Four corners per rectangle (seaward-left, seaward-right, landward-right, landward-left, close)
hw_cos_p <- half_width * cos(perp_angle)
hw_sin_p <- half_width * sin(perp_angle)
d_cos_b <- depths * cos(landward_bearing)
d_sin_b <- depths * sin(landward_bearing)

rect_polys <- vector("list", n_pts)
for (i in seq_len(n_pts)) {
  corners <- matrix(c(
    cx[i] + hw_cos_p[i],              cy[i] + hw_sin_p[i],
    cx[i] - hw_cos_p[i],              cy[i] - hw_sin_p[i],
    cx[i] - hw_cos_p[i] + d_cos_b[i], cy[i] - hw_sin_p[i] + d_sin_b[i],
    cx[i] + hw_cos_p[i] + d_cos_b[i], cy[i] + hw_sin_p[i] + d_sin_b[i],
    cx[i] + hw_cos_p[i],              cy[i] + hw_sin_p[i]
  ), ncol = 2, byrow = TRUE)
  rect_polys[[i]] <- st_polygon(list(corners))
}
cat(sprintf("  Built %d rectangle polygons\n", n_pts))
cat(sprintf("  accom_2050_viable columns: %s\n", paste(names(accom_2050_viable), collapse = ", ")))
cat(sprintf("  accom_2050_viable rows: %d\n", nrow(accom_2050_viable)))
cat(sprintf("  rect_polys length: %d\n", length(rect_polys)))

accom_2050_viable$tier <- ifelse(
  accom_2050_viable$eco_reason %in% c("Both", "Near dune system"), 1, 2)

restoration_polys <- st_sf(
  tier = accom_2050_viable$tier,
  accom_2050_m = accom_2050_viable$accom_2050_m,
  nearest_constraint = accom_2050_viable$nearest_constraint,
  geometry = st_sfc(rect_polys, crs = 27700)
)
cat(sprintf("  Built %d restoration rectangles\n", nrow(restoration_polys)))

# --- Method A: Rectangle approach (MHWS to constraint) ---

# Dissolve overlapping rectangles by tier
restoration_t1 <- st_union(restoration_polys[restoration_polys$tier == 1, ]) %>%
  st_make_valid() %>% st_cast("POLYGON")
restoration_t2 <- st_union(restoration_polys[restoration_polys$tier == 2, ]) %>%
  st_make_valid() %>% st_cast("POLYGON")

rect_sites <- st_sf(
  tier = c(rep(1, length(restoration_t1)), rep(2, length(restoration_t2))),
  geometry = c(restoration_t1, restoration_t2),
  crs = 27700
)

# Clip to land (remove sea) using DEM > 0 as land mask
dem_land <- ifel(cop30_bng > 0, 1, NA)
land_poly <- as.polygons(dem_land) %>% st_as_sf() %>% st_union() %>% st_make_valid()
rect_sites <- st_intersection(rect_sites, land_poly) %>% st_make_valid()
cat(sprintf("  Rectangles after sea clip: %d\n", nrow(rect_sites)))

# Remove areas covered by existing dunes
dunes_union <- st_union(habmos_dunes_bng) %>% st_make_valid()
rect_sites <- st_difference(rect_sites, dunes_union) %>% st_make_valid()
rect_sites <- rect_sites[!st_is_empty(rect_sites), ]
cat(sprintf("  Rectangles after dune removal: %d\n", nrow(rect_sites)))

rect_sites$area_ha <- as.numeric(st_area(rect_sites)) / 10000
cat(sprintf("  Method A (rectangles):\n"))
cat(sprintf("    Tier 1: %.1f ha | Tier 2: %.1f ha\n",
            sum(rect_sites$area_ha[rect_sites$tier == 1]),
            sum(rect_sites$area_ha[rect_sites$tier == 2])))

# --- Method B: Buffer approach (circles around points) ---

# Buffer each viable point by half its accommodation distance, dissolve by tier
accom_2050_viable$buffer_radius <- accom_2050_viable$accom_2050_m / 2
buff_t1 <- accom_2050_viable[accom_2050_viable$tier == 1, ]
buff_t2 <- accom_2050_viable[accom_2050_viable$tier == 2, ]

buff_t1_poly <- st_buffer(buff_t1, buff_t1$buffer_radius) %>% st_union() %>%
  st_make_valid() %>% st_cast("POLYGON")
buff_t2_poly <- st_buffer(buff_t2, buff_t2$buffer_radius) %>% st_union() %>%
  st_make_valid() %>% st_cast("POLYGON")

buff_sites <- st_sf(
  tier = c(rep(1, length(buff_t1_poly)), rep(2, length(buff_t2_poly))),
  geometry = c(buff_t1_poly, buff_t2_poly),
  crs = 27700
)

# Clip to land
buff_sites <- st_intersection(buff_sites, land_poly) %>% st_make_valid()
cat(sprintf("  Buffers after sea clip: %d\n", nrow(buff_sites)))

# Remove existing dunes
buff_sites <- st_difference(buff_sites, dunes_union) %>% st_make_valid()
buff_sites <- buff_sites[!st_is_empty(buff_sites), ]
cat(sprintf("  Buffers after dune removal: %d\n", nrow(buff_sites)))

buff_sites$area_ha <- as.numeric(st_area(buff_sites)) / 10000
cat(sprintf("  Method B (buffers):\n"))
cat(sprintf("    Tier 1: %.1f ha | Tier 2: %.1f ha\n",
            sum(buff_sites$area_ha[buff_sites$tier == 1]),
            sum(buff_sites$area_ha[buff_sites$tier == 2])))

# --- Mapview comparison ---
mapview(rect_sites[rect_sites$tier == 1, ], col.regions = "darkgreen",
        alpha.regions = 0.5, layer.name = "A: Rect Tier 1") +
  mapview(rect_sites[rect_sites$tier == 2, ], col.regions = "olivedrab",
          alpha.regions = 0.5, layer.name = "A: Rect Tier 2") +
  mapview(buff_sites[buff_sites$tier == 1, ], col.regions = "dodgerblue",
          alpha.regions = 0.4, layer.name = "B: Buffer Tier 1") +
  mapview(buff_sites[buff_sites$tier == 2, ], col.regions = "lightskyblue",
          alpha.regions = 0.4, layer.name = "B: Buffer Tier 2") +
  mapview(habmos_dunes_bng, col.regions = "purple", alpha.regions = 0.4,
          layer.name = "Existing Dune Systems") +
  mapview(mhws_modern_diss, color = "black", lwd = 1, layer.name = "MHWS Modern")

# ============================================================
# Save outputs
# ============================================================

cat("\n===== SAVING SAND DUNE OUTPUTS =====\n")

# Intermediate outputs
dir.create("data/processed/sand_dunes", showWarnings = FALSE, recursive = TRUE)

st_write(mhws_modern_diss, "data/processed/sand_dunes/mhws_modern_dissipative.gpkg",
         delete_dsn = TRUE, quiet = TRUE)
cat("  Saved: data/processed/sand_dunes/mhws_modern_dissipative.gpkg\n")

st_write(accom_pts, "data/processed/sand_dunes/accommodation_space_points.gpkg",
         delete_dsn = TRUE, quiet = TRUE)
cat("  Saved: data/processed/sand_dunes/accommodation_space_points.gpkg\n")

st_write(accom_2050_viable, "data/processed/sand_dunes/viable_points_2050.gpkg",
         delete_dsn = TRUE, quiet = TRUE)
cat("  Saved: data/processed/sand_dunes/viable_points_2050.gpkg\n")

st_write(sandy_foreshore, "data/processed/sand_dunes/sandy_foreshore.gpkg",
         delete_dsn = TRUE, quiet = TRUE)
cat("  Saved: data/processed/sand_dunes/sandy_foreshore.gpkg\n")

# Final output: buffered restoration sites (Method B)
st_write(buff_sites, "outputs/habitats/sand_dune_restoration_potential.gpkg",
         delete_dsn = TRUE, quiet = TRUE)
cat("  Saved: outputs/habitats/sand_dune_restoration_potential.gpkg\n")

cat(sprintf("  Total restoration area: %.1f ha (Tier 1: %.1f ha, Tier 2: %.1f ha)\n",
            sum(buff_sites$area_ha),
            sum(buff_sites$area_ha[buff_sites$tier == 1]),
            sum(buff_sites$area_ha[buff_sites$tier == 2])))

# ============================================================
# Combined habitat GPKG: all known + potential layers
# ============================================================

cat("\n===== COMBINED HABITAT GPKG =====\n")

combined_path <- "outputs/habitats/scamp_all_habitats.gpkg"

# Known habitats
sg_known <- st_read("data/processed/seagrass_polygons_harmonised.gpkg", quiet = TRUE)
if (st_crs(sg_known)$epsg != 27700) sg_known <- st_transform(sg_known, 27700)
sg_known$habitat <- "Seagrass"
sg_known$type <- "Known"
sg_known <- st_sf(habitat = sg_known$habitat, type = sg_known$type,
                  geometry = st_geometry(sg_known), crs = 27700)

sm_known <- st_read("data/processed/species/saltmarsh_harmonised.gpkg", quiet = TRUE)
if (st_crs(sm_known)$epsg != 27700) sm_known <- st_transform(sm_known, 27700)
sm_known$habitat <- "Saltmarsh"
sm_known$type <- "Known"
sm_known <- st_sf(habitat = sm_known$habitat, type = sm_known$type,
                  geometry = st_geometry(sm_known), crs = 27700)

oy_known <- st_read("data/processed/species/oyster_known_extent.gpkg", quiet = TRUE)
if (st_crs(oy_known)$epsg != 27700) oy_known <- st_transform(oy_known, 27700)
oy_known$habitat <- "Native Oyster"
oy_known$type <- "Known"
oy_known <- st_sf(habitat = oy_known$habitat, type = oy_known$type,
                  geometry = st_geometry(oy_known), crs = 27700)

dune_known <- habmos_dunes_bng
dune_known$habitat <- "Sand Dune"
dune_known$type <- "Known"
dune_known <- st_sf(habitat = dune_known$habitat, type = dune_known$type,
                    geometry = st_geometry(dune_known), crs = 27700)

# Potential restoration sites
sg_pot <- st_read("outputs/habitats/seagrass_potential_solway.gpkg", quiet = TRUE)
if (st_crs(sg_pot)$epsg != 27700) sg_pot <- st_transform(sg_pot, 27700)
sg_pot$habitat <- "Seagrass"
sg_pot$type <- "Potential"
sg_pot <- st_sf(habitat = sg_pot$habitat, type = sg_pot$type,
                geometry = st_geometry(sg_pot), crs = 27700)

oy_pot <- st_read("outputs/habitats/oyster_potential_solway.gpkg", quiet = TRUE)
if (st_crs(oy_pot)$epsg != 27700) oy_pot <- st_transform(oy_pot, 27700)
oy_pot$habitat <- "Native Oyster"
oy_pot$type <- "Potential"
oy_pot <- st_sf(habitat = oy_pot$habitat, type = oy_pot$type,
                geometry = st_geometry(oy_pot), crs = 27700)

dune_pot <- buff_sites
dune_pot$habitat <- "Sand Dune"
dune_pot$type <- ifelse(dune_pot$tier == 1, "Potential Tier 1", "Potential Tier 2")
dune_pot <- st_sf(habitat = dune_pot$habitat, type = dune_pot$type,
                  geometry = st_geometry(dune_pot), crs = 27700)

# Cast all to MULTIPOLYGON for consistency
all_habitats <- rbind(
  st_cast(sg_known, "MULTIPOLYGON"),
  st_cast(sm_known, "MULTIPOLYGON"),
  st_cast(oy_known, "MULTIPOLYGON"),
  st_cast(dune_known, "MULTIPOLYGON"),
  st_cast(sg_pot, "MULTIPOLYGON"),
  st_cast(oy_pot, "MULTIPOLYGON"),
  st_cast(dune_pot, "MULTIPOLYGON")
)
all_habitats$area_ha <- as.numeric(st_area(all_habitats)) / 10000

st_write(all_habitats, combined_path, delete_dsn = TRUE, quiet = TRUE)
cat(sprintf("  Saved: %s\n", combined_path))

cat("\n  Summary by habitat and type:\n")
hab_summary <- aggregate(area_ha ~ habitat + type, data = st_drop_geometry(all_habitats), FUN = sum)
hab_summary <- hab_summary[order(hab_summary$habitat, hab_summary$type), ]
for (i in seq_len(nrow(hab_summary))) {
  cat(sprintf("    %-15s %-20s %.1f ha\n",
              hab_summary$habitat[i], hab_summary$type[i], hab_summary$area_ha[i]))
}
cat(sprintf("  Total features: %d\n", nrow(all_habitats)))

