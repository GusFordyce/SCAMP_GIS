# rebuild_clusters.R — Quick rebuild of cluster boxes from intermediate outputs
# Run this instead of the full 03_habitats.R when you just need the cluster gpkg files.
# Reads saved harmonised/processed data and applies the same clustering algorithm.

library(sf)
library(dplyr)

cat("========== Rebuilding cluster boxes ==========\n\n")

dir.create("data/processed/habitats", recursive = TRUE, showWarnings = FALSE)

aoi <- st_read("data/processed/boundaries/scamp_aoi.gpkg", quiet = TRUE)

make_box <- function(bb, name) {
  st_sf(name = name,
        geometry = st_sfc(st_polygon(list(matrix(c(
          bb["xmin"], bb["ymin"], bb["xmax"], bb["ymin"],
          bb["xmax"], bb["ymax"], bb["xmin"], bb["ymax"],
          bb["xmin"], bb["ymin"]
        ), ncol = 2, byrow = TRUE))), crs = 27700))
}

cluster_features <- function(sf_data, buf_dist = 2000) {
  centroids <- suppressWarnings(st_centroid(sf_data))
  buf <- st_buffer(centroids, dist = buf_dist) %>% st_union() %>% st_cast("POLYGON")
  ids <- st_intersects(centroids, buf)
  sapply(ids, `[`, 1)
}

merge_boxes <- function(boxes, merge_dist = 1000) {
  if (nrow(boxes) <= 1) return(boxes)
  buf <- st_buffer(boxes, dist = merge_dist) %>% st_union() %>% st_cast("POLYGON")
  grp <- sapply(st_intersects(boxes, buf), `[`, 1)
  do.call(rbind, lapply(sort(unique(grp)), function(g) {
    bb <- st_bbox(boxes[grp == g, ])
    make_box(bb, sprintf("Cluster %d", g))
  }))
}


# ============================================================
# 1. SEAGRASS
# ============================================================
cat("--- Seagrass ---\n")
sg_polys <- st_read("data/processed/seagrass_polygons_harmonised.gpkg", quiet = TRUE)

sg_polys$cluster <- cluster_features(sg_polys, buf_dist = 2000)

sg_boxes <- do.call(rbind, lapply(sort(unique(sg_polys$cluster)), function(cl) {
  bb <- st_bbox(st_buffer(sg_polys[sg_polys$cluster == cl, ], dist = 2000))
  make_box(bb, sprintf("Cluster %d", cl))
}))
sg_boxes <- merge_boxes(sg_boxes, merge_dist = 1000)

# Add indicator clusters (Zostera records not overlapping existing boxes)
ind_path <- "outputs/species/indicator_records.gpkg"
if (file.exists(ind_path)) {
  ind_zostera <- st_read(ind_path, quiet = TRUE) %>%
    filter(grepl("Zostera", latin_name, ignore.case = TRUE)) %>%
    st_transform(27700) %>% st_make_valid() %>%
    st_intersection(aoi)
  ind_zostera <- st_geometry(ind_zostera) %>% st_sf()

  if (nrow(ind_zostera) > 0) {
    ind_zostera$cluster <- cluster_features(ind_zostera, buf_dist = 2000)
    ind_boxes <- do.call(rbind, lapply(sort(unique(ind_zostera$cluster)), function(cl) {
      bb <- st_bbox(ind_zostera[ind_zostera$cluster == cl, ])
      make_box(bb, sprintf("Ind_%d", cl))
    }))
    # Keep only indicator clusters that don't overlap existing boxes
    ind_new <- ind_boxes[lengths(st_intersects(ind_boxes, sg_boxes)) == 0, ]
    if (nrow(ind_new) > 0) {
      for (j in seq_len(nrow(ind_new))) {
        bb_j <- st_bbox(st_buffer(ind_new[j, ], dist = 2000))
        sg_boxes <- rbind(sg_boxes, make_box(bb_j, sprintf("Indicator_%d", j)))
      }
      sg_boxes <- merge_boxes(sg_boxes, merge_dist = 1000)
    }
  }
}

# Geographic names west-to-east
sg_boxes <- sg_boxes[order(st_coordinates(suppressWarnings(st_centroid(sg_boxes)))[, 1]), ]
sg_place_names <- c("Loch Ryan", "Wigtown Bay", "Fleet Bay", "Kirkcudbright Bay",
                    "Auchencairn Bay & Rough Island", "Nith Estuary")
if (nrow(sg_boxes) <= length(sg_place_names)) {
  sg_boxes$name <- sg_place_names[seq_len(nrow(sg_boxes))]
} else {
  sg_boxes$name <- c(sg_place_names,
    sprintf("Cluster %d", seq(length(sg_place_names) + 1, nrow(sg_boxes))))
}

st_write(sg_boxes, "data/processed/habitats/seagrass_cluster_boxes.gpkg",
         delete_layer = TRUE, quiet = TRUE)
cat(sprintf("  %d clusters -> data/processed/habitats/seagrass_cluster_boxes.gpkg\n",
            nrow(sg_boxes)))
for (k in seq_len(nrow(sg_boxes))) cat(sprintf("    %s\n", sg_boxes$name[k]))


# ============================================================
# 2. SALTMARSH
# ============================================================
cat("\n--- Saltmarsh ---\n")
sm_polys <- st_read("data/processed/species/saltmarsh_harmonised.gpkg", quiet = TRUE)

sm_polys$cluster <- cluster_features(sm_polys, buf_dist = 2000)

sm_cluster_boxes <- do.call(rbind, lapply(sort(unique(sm_polys$cluster)), function(cl) {
  bb <- st_bbox(sm_polys[sm_polys$cluster == cl, ])
  st_sf(cluster = cl,
        n = sum(sm_polys$cluster == cl),
        name = paste0("Saltmarsh_", cl),
        geometry = st_as_sfc(bb, crs = 27700))
}))

st_write(sm_cluster_boxes, "data/processed/habitats/saltmarsh_cluster_boxes.gpkg",
         delete_layer = TRUE, quiet = TRUE)
cat(sprintf("  %d clusters -> data/processed/habitats/saltmarsh_cluster_boxes.gpkg\n",
            nrow(sm_cluster_boxes)))


# ============================================================
# 3. OYSTER
# ============================================================
cat("\n--- Oyster ---\n")

# Stage 1: known extent polygons
oy_known <- st_read("data/processed/species/oyster_known_extent.gpkg", quiet = TRUE)
if (nrow(oy_known) > 1) {
  oy_known$cluster <- cluster_features(oy_known, buf_dist = 2000)
} else if (nrow(oy_known) == 1) {
  oy_known$cluster <- 1L
}

if (nrow(oy_known) > 0) {
  oy_boxes <- do.call(rbind, lapply(sort(unique(oy_known$cluster)), function(cl) {
    bb <- st_bbox(st_buffer(oy_known[oy_known$cluster == cl, ], dist = 2000))
    make_box(bb, sprintf("Known_%d", cl))
  }))
} else {
  oy_boxes <- st_sf(name = character(0), geometry = st_sfc(crs = 27700))
}
cat(sprintf("  Known extent clusters: %d\n", nrow(oy_boxes)))

# Stage 2: merge with historical locations
hist_path <- "data/processed/species/oyster_historical_locations.gpkg"
if (file.exists(hist_path)) {
  hist_oy <- st_read(hist_path, quiet = TRUE)
  hist_geom_types <- st_geometry_type(hist_oy)
  hist_pts <- hist_oy[hist_geom_types %in% c("POINT", "MULTIPOINT"), ]
  hist_py  <- hist_oy[hist_geom_types %in% c("POLYGON", "MULTIPOLYGON"), ]

  hist_all_pts <- bind_rows(
    if (nrow(hist_pts) > 0) st_sf(geometry = st_geometry(hist_pts)),
    if (nrow(hist_py) > 0) st_sf(geometry = st_geometry(suppressWarnings(st_centroid(hist_py))))
  )

  if (nrow(hist_all_pts) > 1) {
    hist_all_pts$cluster <- cluster_features(hist_all_pts, buf_dist = 2000)
    hist_boxes <- do.call(rbind, lapply(sort(unique(hist_all_pts$cluster)), function(cl) {
      bb <- st_bbox(st_buffer(hist_all_pts[hist_all_pts$cluster == cl, ], dist = 2000))
      make_box(bb, sprintf("Hist_%d", cl))
    }))
  } else if (nrow(hist_all_pts) == 1) {
    bb <- st_bbox(st_buffer(hist_all_pts, dist = 2000))
    hist_boxes <- make_box(bb, "Hist_1")
  } else {
    hist_boxes <- st_sf(name = character(0), geometry = st_sfc(crs = 27700))
  }

  if (nrow(hist_boxes) > 0 && nrow(oy_boxes) > 0) {
    combined <- rbind(oy_boxes, hist_boxes)
    oy_boxes <- merge_boxes(combined, merge_dist = 1000)
  } else if (nrow(hist_boxes) > 0) {
    oy_boxes <- hist_boxes
  }
  cat(sprintf("  After historical merge: %d clusters\n", nrow(oy_boxes)))
}

# Stage 3: merge with indicator records (Ostrea)
if (file.exists(ind_path)) {
  ind_oy <- st_read(ind_path, quiet = TRUE) %>%
    filter(habitat == "Oyster") %>%
    filter(grepl("ostrea", latin_name, ignore.case = TRUE)) %>%
    st_transform(27700) %>% st_make_valid() %>%
    st_intersection(aoi)

  if (nrow(ind_oy) > 1) {
    ind_oy$cluster <- cluster_features(ind_oy, buf_dist = 2000)
    ind_oy_boxes <- do.call(rbind, lapply(sort(unique(ind_oy$cluster)), function(cl) {
      bb <- st_bbox(st_buffer(ind_oy[ind_oy$cluster == cl, ], dist = 2000))
      make_box(bb, sprintf("Ind_%d", cl))
    }))
  } else if (nrow(ind_oy) == 1) {
    bb <- st_bbox(st_buffer(ind_oy, dist = 2000))
    ind_oy_boxes <- make_box(bb, "Ind_1")
  } else {
    ind_oy_boxes <- st_sf(name = character(0), geometry = st_sfc(crs = 27700))
  }

  if (nrow(ind_oy_boxes) > 0 && nrow(oy_boxes) > 0) {
    combined <- rbind(oy_boxes, ind_oy_boxes)
    oy_boxes <- merge_boxes(combined, merge_dist = 1000)
  } else if (nrow(ind_oy_boxes) > 0) {
    oy_boxes <- ind_oy_boxes
  }
  cat(sprintf("  After indicator merge: %d clusters\n", nrow(oy_boxes)))
}

st_write(oy_boxes, "data/processed/habitats/oyster_cluster_boxes.gpkg",
         delete_layer = TRUE, quiet = TRUE)
cat(sprintf("  %d clusters -> data/processed/habitats/oyster_cluster_boxes.gpkg\n",
            nrow(oy_boxes)))


# ============================================================
# 4. SAND DUNES
# ============================================================
cat("\n--- Sand Dunes ---\n")
dunes <- st_read("outputs/habitats/sand_dunes_habmos_solway.gpkg", quiet = TRUE)

dunes$cluster <- cluster_features(dunes, buf_dist = 2000)

dune_boxes <- do.call(rbind, lapply(sort(unique(dunes$cluster)), function(cl) {
  bb <- st_bbox(st_buffer(dunes[dunes$cluster == cl, ], dist = 1000))
  make_box(bb, sprintf("Dune_%d", cl))
}))
dune_boxes <- merge_boxes(dune_boxes, merge_dist = 1000)

dune_boxes <- dune_boxes[order(st_coordinates(suppressWarnings(st_centroid(dune_boxes)))[, 1]), ]
dune_boxes$name <- sprintf("Dune cluster %d", seq_len(nrow(dune_boxes)))

st_write(dune_boxes, "data/processed/habitats/dune_cluster_boxes.gpkg",
         delete_layer = TRUE, quiet = TRUE)
cat(sprintf("  %d clusters -> data/processed/habitats/dune_cluster_boxes.gpkg\n",
            nrow(dune_boxes)))
for (j in seq_len(nrow(dune_boxes))) {
  n_in <- sum(lengths(st_intersects(dunes, dune_boxes[j, ])) > 0)
  cat(sprintf("    %s: %d polygons\n", dune_boxes$name[j], n_in))
}

# Also save historical oyster extent (polygon subset) if not yet saved
hist_ext_path <- "data/processed/species/oyster_historical_extent.gpkg"
if (!file.exists(hist_ext_path) && file.exists(hist_path)) {
  hist_oy <- st_read(hist_path, quiet = TRUE)
  hist_py <- hist_oy[st_geometry_type(hist_oy) %in% c("POLYGON", "MULTIPOLYGON"), ]
  if (nrow(hist_py) > 0) {
    oy_historical <- st_sf(source = "Historical", geometry = st_geometry(hist_py))
    st_write(oy_historical, hist_ext_path, delete_layer = TRUE, quiet = TRUE)
    cat(sprintf("\n  Also saved: %s (%d polygons)\n", hist_ext_path, nrow(oy_historical)))
  }
}

cat("\n========== Cluster boxes rebuilt ==========\n")
