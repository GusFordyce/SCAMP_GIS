##=====================================================================
## 00_rebuild_missing.R
## Regenerates outputs/catchments/rivers_with_accumulated_risk_hydrosheds.gpkg
## in a fresh session, without running the whole of 05_localthreats.R.
##
## Code taken from 05_localthreats.R:
##   L1021-1038  OpenRivers clipped to catchments  -> rivers_crop
##   L2495-2618  Step 4 (revised) HydroSHEDS routing -> risk_accum
##               (skipped if risk_accumulated_hydrosheds.tif already exists)
##   L2633-2689  Step 4b river accumulated risk     -> the gpkg
##
## Land cover: nothing to rebuild. Use hlcm_habmos_lca_reclass.tif and
## hlcm_lc_lookup_reclass.csv (07_plots.R L3840) wherever the corrected HLCM
## raster is needed. Woodland classes are identical and the lookup adds G4.F.
##
## Run from the project root:  source("scripts/00_rebuild_missing.R")
##=====================================================================

suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(dplyr)
})
Sys.setenv(SHAPE_RESTORE_SHX = "YES")
terraOptions(memfrac = 0.6, progress = 0)
t0 <- Sys.time()

out_path        <- "outputs/catchments/rivers_with_accumulated_risk_hydrosheds.gpkg"
risk_accum_path <- "outputs/catchments/risk_accumulated_hydrosheds.tif"

## ---- Catchments (CEFAS crop, saved by 05_localthreats.R L1910) ----------
cefas_candidates <- c("outputs/catchments/cefas_catchments_crop.shp",
                      "outputs/catchments/cefas_catchments_crop.gpkg")
cefas_file <- cefas_candidates[file.exists(cefas_candidates)][1]
if (is.na(cefas_file)) stop("cefas_catchments_crop not found in outputs/catchments/")
cefas_crop <- st_read(cefas_file, quiet = TRUE) %>% st_transform(27700) %>% st_make_valid()
catch_u    <- st_union(cefas_crop)
cat(sprintf("Catchments: %d features from %s\n", nrow(cefas_crop), cefas_file))

## ---- Rivers clipped to catchments (L1021-1038) -------------------------
rivers <- st_read("data/raw/physical/OpenRivers/Data/oprvrs_gb.gpkg",
                  layer = "watercourse_link", quiet = TRUE,
                  wkt_filter = st_as_text(st_as_sfc(st_bbox(catch_u)))) %>%
  st_transform(27700) %>% st_make_valid()
rivers_crop <- st_intersection(rivers, catch_u) %>%
  st_make_valid() %>%
  st_collection_extract("LINESTRING")
rm(rivers); invisible(gc())
cat(sprintf("Rivers clipped to catchments: %d segments\n", nrow(rivers_crop)))

## ---- Accumulated risk raster -------------------------------------------
if (file.exists(risk_accum_path)) {

  risk_accum <- rast(risk_accum_path)
  cat("Loaded existing", risk_accum_path, "\n")

} else {

  cat("risk_accumulated_hydrosheds.tif not found: running Step 4 (revised)\n")

  diffuse_pollution_risk <- rast("outputs/catchments/diffuse_pollution_risk.tif")

  las_dg_bng <- st_read("data/raw/boundaries/Local_Authority_Boundaries_-_Scotland/pub_las.shp",
                        quiet = TRUE) %>% st_transform(27700)
  las_dg_bng <- las_dg_bng[grepl("dumfries", las_dg_bng$local_auth, ignore.case = TRUE), ]
  cefas_land <- st_intersection(st_make_valid(catch_u), st_union(las_dg_bng))

  dir_raw <- rast("data/raw/physical/eu_dir_3s/eu_dir_3s.tif")
  cefas_land_buf_wgs <- st_transform(st_buffer(cefas_land, 2000), crs(dir_raw))
  dir_crop <- crop(dir_raw, ext(vect(cefas_land_buf_wgs)))

  risk_on_dir <- terra::project(diffuse_pollution_risk, dir_crop, method = "average")
  risk_vals <- values(risk_on_dir)[, 1]
  fdir_vals <- values(dir_crop)[, 1]
  nr <- nrow(dir_crop); nc <- ncol(dir_crop); ncell_dir <- ncell(dir_crop)

  risk_vals[is.na(risk_vals)] <- 0
  risk_vals[risk_vals < 0] <- 0
  accum_vals <- risk_vals

  d8_dr <- c(`1`=0L, `2`=1L, `4`=1L, `8`=1L, `16`=0L, `32`=-1L, `64`=-1L, `128`=-1L)
  d8_dc <- c(`1`=1L, `2`=1L, `4`=0L, `8`=-1L, `16`=-1L, `32`=-1L, `64`=0L, `128`=1L)

  valid_cells <- which(!is.na(fdir_vals))
  r_all <- ((valid_cells - 1L) %/% nc) + 1L
  c_all <- ((valid_cells - 1L) %% nc) + 1L
  dir_key <- as.character(as.integer(fdir_vals[valid_cells]))
  has_dir <- dir_key %in% names(d8_dr)

  target_idx <- rep(NA_integer_, ncell_dir)
  r2 <- r_all; c2 <- c_all
  r2[has_dir] <- r_all[has_dir] + d8_dr[dir_key[has_dir]]
  c2[has_dir] <- c_all[has_dir] + d8_dc[dir_key[has_dir]]
  in_bounds <- has_dir & r2 >= 1L & r2 <= nr & c2 >= 1L & c2 <= nc
  tgt <- rep(NA_integer_, length(valid_cells))
  tgt[in_bounds] <- (r2[in_bounds] - 1L) * nc + c2[in_bounds]
  target_idx[valid_cells] <- tgt

  in_degree <- tabulate(target_idx[!is.na(target_idx)], nbins = ncell_dir)

  # Kahn's algorithm: sum every upstream contributor before passing downstream
  queue <- integer(length(valid_cells))
  zero_indeg <- valid_cells[in_degree[valid_cells] == 0L]
  head_i <- 1L; tail_i <- length(zero_indeg)
  queue[seq_len(tail_i)] <- zero_indeg

  cat(sprintf("Routing risk through %d cells...\n", length(valid_cells)))
  while (head_i <= tail_i) {
    ci <- queue[head_i]; head_i <- head_i + 1L
    tgt_ci <- target_idx[ci]
    if (!is.na(tgt_ci)) {
      accum_vals[tgt_ci] <- accum_vals[tgt_ci] + accum_vals[ci]
      in_degree[tgt_ci] <- in_degree[tgt_ci] - 1L
      if (in_degree[tgt_ci] == 0L) {
        tail_i <- tail_i + 1L
        queue[tail_i] <- tgt_ci
      }
    }
  }
  n_unresolved <- sum(in_degree[valid_cells] > 0L)
  if (n_unresolved > 0) warning(n_unresolved, " cells never resolved (possible DIR cycle)")

  risk_accum_dir <- dir_crop
  values(risk_accum_dir) <- accum_vals
  risk_accum <- terra::project(risk_accum_dir, "EPSG:27700", method = "bilinear")
  writeRaster(risk_accum, risk_accum_path, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
  cat("Saved", risk_accum_path, "\n")
  rm(risk_vals, fdir_vals, accum_vals, target_idx, in_degree, queue); invisible(gc())
}

## ---- Step 4b: river accumulated risk (L2633-2689) ----------------------
rivers_in_catch <- rivers_crop
n_seg <- nrow(rivers_in_catch)
rivers_in_catch$seg_length_m <- as.numeric(st_length(rivers_in_catch))

# first and last vertex of each segment, vectorised
cc <- st_coordinates(rivers_in_catch)
seg_id <- cc[, ncol(cc)]
first_i <- !duplicated(seg_id)
last_i  <- !duplicated(seg_id, fromLast = TRUE)
pts_first <- vect(cc[first_i, 1:2, drop = FALSE], crs = "EPSG:27700")
pts_last  <- vect(cc[last_i,  1:2, drop = FALSE], crs = "EPSG:27700")

risk_at_first <- terra::extract(risk_accum, pts_first)[, 2]
risk_at_last  <- terra::extract(risk_accum, pts_last)[, 2]
risk_at_first[is.na(risk_at_first)] <- 0
risk_at_last[is.na(risk_at_last)]   <- 0

rivers_in_catch$accum_risk_raw <- pmax(risk_at_first, risk_at_last)
rivers_in_catch$accum_risk_log <- log1p(rivers_in_catch$accum_risk_raw)
lmin <- min(rivers_in_catch$accum_risk_log, na.rm = TRUE)
lmax <- max(rivers_in_catch$accum_risk_log, na.rm = TRUE)
rivers_in_catch$accum_risk <- if (lmax > lmin) {
  (rivers_in_catch$accum_risk_log - lmin) / (lmax - lmin)
} else rep(0, n_seg)

st_write(rivers_in_catch, out_path, delete_dsn = TRUE, quiet = TRUE)

## ---- Check ------------------------------------------------------------
chk <- st_read(out_path, quiet = TRUE)
cat(sprintf("\nSaved %s\n  %d features | EPSG:%s | accum_risk range %.3f-%.3f\n",
            out_path, nrow(chk), st_crs(chk)$epsg,
            min(chk$accum_risk), max(chk$accum_risk)))
cat(sprintf("Done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
