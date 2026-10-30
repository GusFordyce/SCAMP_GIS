##=====================================================================
## 09_handover_package.R — partner handover package build
##=====================================================================
## Packages the derived SCAMP outputs for partners who do not use R:
## vectors in one GeoPackage, rasters as COGs at native resolution, and the
## two aligned cubes from outputs/design/. No code is shipped; METHODS.md is
## the technical record.
##
## The curated manifest decides what ships. This script does not decide —
## it audits (stage 1) and builds (stage 2). Anything it cannot verify it
## reports rather than guesses, because a partner cannot check our working.
##
## Run top-to-bottom with source(). Two passes:
##   STAGE 1 — audit only. Writes outputs/handover/gap_report.md and
##             nothing else. Fix the manifest, re-run until the report is
##             clean.
##   STAGE 2 — build outputs/handover/SCAMP_GIS_Handover_v1/.
##
## Hand-written inputs, never generated here:
##   outputs/handover/README.md
##   outputs/handover/METHODS.md
##   scripts/handover_metadata.R
##   outputs/handover/manifest_curated.csv  (or the draft .xlsx)
##=====================================================================

suppressMessages({
  library(sf)
  library(terra)
})

##=====================================================================
## 0. CONFIG
##=====================================================================

STAGE <- 2     # 1 = audit, 2 = build
ZIP   <- FALSE   # stage 2 only: also write SCAMP_GIS_Handover_v1.zip

## Fixed because README.md documents this exact layout.
PKG_NAME <- "SCAMP_GIS_Handover_v1"

HANDOVER_DIR <- "outputs/handover"
PKG_DIR      <- file.path(HANDOVER_DIR, PKG_NAME)
DESIGN_DIR   <- "outputs/design"

PATH_README    <- file.path(HANDOVER_DIR, "README.md")
PATH_METHODS   <- file.path(HANDOVER_DIR, "METHODS.md")
PATH_METADATA  <- "scripts/handover_metadata.R"
PATH_MANIFEST  <- file.path(HANDOVER_DIR, "manifest_curated.csv")
PATH_MANIFEST_XLSX <- file.path(HANDOVER_DIR, "SCAMP_manifest_draft.xlsx")
PATH_AOI       <- "data/processed/boundaries/scamp_aoi.gpkg"
PATH_INVENTORY <- "PROJECT_INVENTORY.md"
PATH_LC_LOOKUP <- "data/processed/habitats/hlcm_lc_lookup_reclass.csv"

TARGET_EPSG   <- 27700
SIZE_WARN_MB  <- 500

## Roots searched when a manifest path is missing and for the orphan check.
SCAN_ROOTS <- c("outputs", "data/processed")

## Spatial file types the orphan check (stage 1c) considers.
SPATIAL_EXT <- c("tif", "tiff", "gpkg", "shp")

terraOptions(memfrac = 0.5, progress = 0)

if (!STAGE %in% c(1, 2)) {
  stop("\n\nSTAGE must be 1 (audit) or 2 (build); it is currently ", STAGE, ".\n",
       call. = FALSE)
}

cat(sprintf("========== 09_handover_package.R (STAGE %s) ==========\n\n", STAGE))

##=====================================================================
## 1. SMALL UTILITIES
##=====================================================================

## Warnings are collected rather than emitted, so the final report is a
## single place to look instead of a scrollback hunt.
WARNINGS <- character(0)
warn_note <- function(...) {
  msg <- paste0(...)
  WARNINGS <<- c(WARNINGS, msg)
  cat("  ! ", msg, "\n", sep = "")
  invisible(msg)
}

stop_clear <- function(...) stop(paste0("\n\n", paste0(...), "\n"), call. = FALSE)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a

blank <- function(x) is.na(x) | !nzchar(trimws(as.character(x)))

size_mb <- function(paths) {
  s <- file.size(paths)
  s[is.na(s)] <- 0
  round(sum(s) / 1024^2, 3)
}

## Forward slashes and no leading "./", so a manifest path and a list.files()
## path for the same file compare equal.
fwd <- function(p) sub("^\\./", "", gsub("\\\\", "/", p))

ext_of <- function(p) tolower(sub(".*\\.", "", basename(p)))

## Em-dashes do not survive every partner's spreadsheet import, and a dash
## separating a label from its detail reads as a colon anyway. The first
## " \u2014 " in a string becomes ": " and any later one ", ". En-dashes in
## number ranges are left alone. Applied to generated text, not to code.
de_emdash <- function(x) {
  x <- sub(" \u2014 ", ": ", x, fixed = TRUE)
  gsub(" \u2014 ", ", ", x, fixed = TRUE)
}

base_noext <- function(p) sub("\\.[^.]*$", "", basename(p))

## Semicolon-separated manifest cells -> trimmed character vector.
split_semi <- function(x) {
  if (length(x) == 0 || blank(x)) return(character(0))
  out <- trimws(unlist(strsplit(as.character(x), ";", fixed = TRUE)))
  out[nzchar(out)]
}

## Nested handover_meta lookup that tolerates fields the analyst has not
## added yet, so a config edit cannot crash the build.
meta_get <- function(..., default = NA_character_) {
  keys <- c(...)
  x <- handover_meta
  for (k in keys) {
    if (!is.list(x) || is.null(x[[k]])) return(default)
    x <- x[[k]]
  }
  if (length(x) == 0) default else x
}

## Any TODO left in the config is a reporting matter, not a build stopper —
## a package with a placeholder contact is still better than no package,
## provided nobody can miss the placeholder.
find_todos <- function(x, path = "handover_meta") {
  out <- character(0)
  if (is.list(x)) {
    nms <- names(x) %||% rep("", length(x))
    for (i in seq_along(x)) {
      lbl <- if (nzchar(nms[i])) paste0(path, "$", nms[i]) else sprintf("%s[[%d]]", path, i)
      out <- c(out, find_todos(x[[i]], lbl))
    }
  } else if (is.character(x) && any(grepl("TODO", x, fixed = TRUE))) {
    out <- c(out, sprintf("%s: %s", path, paste(x[grepl("TODO", x, fixed = TRUE)], collapse = " | ")))
  }
  out
}

##=====================================================================
## 2. REQUIRED INPUTS
##=====================================================================

missing_inputs <- character(0)
if (!file.exists(PATH_README))   missing_inputs <- c(missing_inputs, PATH_README)
if (!file.exists(PATH_METHODS))  missing_inputs <- c(missing_inputs, PATH_METHODS)
if (!file.exists(PATH_METADATA)) missing_inputs <- c(missing_inputs, PATH_METADATA)

manifest_src <- NULL
if (file.exists(PATH_MANIFEST)) {
  manifest_src <- PATH_MANIFEST
} else if (file.exists(PATH_MANIFEST_XLSX)) {
  manifest_src <- PATH_MANIFEST_XLSX
} else {
  missing_inputs <- c(missing_inputs, paste0(PATH_MANIFEST, " (or ", PATH_MANIFEST_XLSX, ")"))
}

if (length(missing_inputs) > 0) {
  stop_clear(
    "09_handover_package.R cannot run. These hand-written inputs are missing:\n",
    paste0("  - ", missing_inputs, collapse = "\n"),
    "\n\nThis script never generates them. README.md and METHODS.md are the\n",
    "technical record and are copied into the package verbatim; the manifest\n",
    "decides what ships; handover_metadata.R carries everything a person has\n",
    "to decide (title, abstract, owner, contact, licence, keywords)."
  )
}

source(PATH_METADATA)
if (!exists("handover_meta") || !is.list(handover_meta)) {
  stop_clear(PATH_METADATA, " was sourced but did not define a list called handover_meta.")
}
cat(sprintf("Metadata config: %s\n", PATH_METADATA))

TODOS <- find_todos(handover_meta)

##=====================================================================
## 3. MANIFEST
##=====================================================================

MANIFEST_COLS <- c("include", "layer_name", "path", "theme", "geom_type",
                   "resolution_m", "crs_epsg", "units", "description",
                   "derived_from", "drop_fields", "exclude_reason", "review_note")

read_manifest <- function(src) {
  if (grepl("\\.csv$", src, ignore.case = TRUE)) {
    m <- utils::read.csv(src, stringsAsFactors = FALSE, colClasses = "character",
                         check.names = FALSE, na.strings = c("NA", ""))
  } else {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop_clear("The manifest is an .xlsx (", src, ") but the readxl package is not\n",
                 "installed. Either install.packages(\"readxl\") or export sheet\n",
                 "\"manifest\" to ", PATH_MANIFEST, ".")
    }
    sheets <- readxl::excel_sheets(src)
    if (!"manifest" %in% sheets) {
      stop_clear(src, " has no sheet called \"manifest\". Sheets found: ",
                 paste(sheets, collapse = ", "), ".")
    }
    m <- as.data.frame(readxl::read_excel(src, sheet = "manifest", col_types = "text"),
                       stringsAsFactors = FALSE)
  }
  names(m) <- trimws(names(m))
  missing_cols <- setdiff(MANIFEST_COLS, names(m))
  if (length(missing_cols) > 0) {
    stop_clear("The manifest (", src, ") is missing required column(s):\n  ",
               paste(missing_cols, collapse = ", "),
               "\n\nExpected columns:\n  ", paste(MANIFEST_COLS, collapse = ", "))
  }
  for (cl in MANIFEST_COLS) {
    m[[cl]] <- trimws(as.character(m[[cl]]))
    m[[cl]][is.na(m[[cl]])] <- ""
  }
  m[, MANIFEST_COLS, drop = FALSE]
}

as_include <- function(x) {
  v <- toupper(trimws(as.character(x)))
  v %in% c("TRUE", "T", "YES", "Y", "1")
}

manifest <- read_manifest(manifest_src)
manifest$include_flag <- as_include(manifest$include)
manifest$row_id <- seq_len(nrow(manifest))

cat(sprintf("Manifest: %s — %d rows, %d marked include\n",
            manifest_src, nrow(manifest), sum(manifest$include_flag)))

dup_names <- manifest$layer_name[manifest$include_flag & nzchar(manifest$layer_name)]
dup_names <- dup_names[duplicated(dup_names)]
if (length(dup_names) > 0) {
  warn_note("Duplicate layer_name among include = TRUE rows: ",
            paste(unique(dup_names), collapse = ", "))
}

##=====================================================================
## 4. PATH EXPANSION
##=====================================================================
## Wildcard rows let the manifest name a set (layers_raw/*.tif) without
## listing 14 files. Each expanded file needs its own published name: if
## layer_name itself carries a *, the matched text is substituted into it,
## otherwise the matched text is appended. Stage 1b prints the result so
## the naming is reviewed, not assumed.

wildcard_match <- function(pattern, path) {
  rx <- gsub("\\*", "(.*)", gsub("([.|()\\^{}+$?\\[\\]])", "\\\\\\1", fwd(pattern)))
  m <- regmatches(fwd(path), regexec(paste0("^", rx, "$"), fwd(path)))[[1]]
  if (length(m) < 2) "" else paste(m[-1], collapse = "_")
}

expand_row <- function(row) {
  p <- fwd(row$path)
  if (!grepl("*", p, fixed = TRUE)) {
    return(data.frame(row_id = row$row_id, layer_name = row$layer_name,
                      file = p, wildcard = FALSE, stringsAsFactors = FALSE))
  }
  hits <- fwd(Sys.glob(p))
  hits <- hits[!dir.exists(hits)]
  if (length(hits) == 0) {
    return(data.frame(row_id = row$row_id, layer_name = row$layer_name,
                      file = character(0), wildcard = TRUE, stringsAsFactors = FALSE))
  }
  hits <- sort(hits)
  nms <- vapply(hits, function(h) {
    tok <- wildcard_match(p, h)
    tok <- gsub("[^A-Za-z0-9]+", "_", tok)
    tok <- gsub("^_+|_+$", "", tok)
    if (!nzchar(tok)) tok <- base_noext(h)
    if (grepl("*", row$layer_name, fixed = TRUE)) {
      sub("*", tok, row$layer_name, fixed = TRUE)
    } else if (length(hits) == 1) {
      row$layer_name
    } else {
      paste0(row$layer_name, "_", tok)
    }
  }, character(1), USE.NAMES = FALSE)
  data.frame(row_id = row$row_id, layer_name = nms, file = hits,
             wildcard = TRUE, stringsAsFactors = FALSE)
}

expanded <- do.call(rbind, lapply(split(manifest, manifest$row_id),
                                  function(r) expand_row(r)))
expanded <- expanded[order(expanded$row_id), , drop = FALSE]

## Carry the manifest attributes onto each expanded file.
shipped <- merge(expanded,
                 manifest[, c("row_id", "include", "include_flag", "theme", "geom_type",
                              "resolution_m", "crs_epsg", "units", "description",
                              "derived_from", "drop_fields", "exclude_reason",
                              "review_note", "path")],
                 by = "row_id", all.x = TRUE)
shipped <- shipped[order(shipped$row_id, shipped$file), , drop = FALSE]
names(shipped)[names(shipped) == "path"] <- "manifest_path"

## Type is taken from the extension, not the geom_type column: a manifest
## typo must not route a raster into the GeoPackage.
shipped$type <- ifelse(ext_of(shipped$file) %in% c("tif", "tiff"), "raster",
                ifelse(ext_of(shipped$file) %in% c("gpkg", "shp", "geojson"), "vector",
                ifelse(ext_of(shipped$file) %in% c("csv"), "table", "other")))

inc <- shipped[shipped$include_flag, , drop = FALSE]

##=====================================================================
## 5. PROBING LAYERS
##=====================================================================
## Everything here opens headers only. terra::rast() is lazy and st_layers()
## reads the GeoPackage catalogue, so nothing loads cell values or geometry.

probe_raster <- function(path) {
  r <- try(terra::rast(path), silent = TRUE)
  if (inherits(r, "try-error")) {
    return(list(ok = FALSE, msg = trimws(conditionMessage(attr(r, "condition")))))
  }
  cr <- terra::crs(r, describe = TRUE)
  epsg <- suppressWarnings(as.integer(cr$code[1]))
  list(ok = TRUE, nlyr = terra::nlyr(r), res = terra::res(r)[1],
       epsg = epsg, lonlat = terra::is.lonlat(r),
       datatype = terra::datatype(r), names = names(r),
       nodata = terra::NAflag(r), ncol = terra::ncol(r), nrow = terra::nrow(r),
       categorical = terra::is.factor(r))
}

probe_vector <- function(path) {
  lyr <- try(sf::st_layers(path), silent = TRUE)
  if (inherits(lyr, "try-error")) {
    return(list(ok = FALSE, msg = trimws(conditionMessage(attr(lyr, "condition")))))
  }
  epsg <- vapply(lyr$crs, function(crs_i) {
    e <- try(sf::st_crs(crs_i)$epsg, silent = TRUE)
    if (inherits(e, "try-error") || is.null(e)) NA_integer_ else as.integer(e)
  }, integer(1))
  list(ok = TRUE, layers = lyr$name,
       features = as.integer(lyr$features),
       geomtype = vapply(lyr$geomtype, function(g) paste(unlist(g), collapse = "/"), character(1)),
       epsg = epsg)
}

##=====================================================================
## 6. FIELD DICTIONARIES FROM PROJECT_INVENTORY.md
##=====================================================================
## Parsed rather than transcribed so the inventory stays the single source
## of truth. Section 7 tables are "| Field | Meaning | Units |" under
## headings that name the layers they describe.

parse_field_dictionaries <- function(path) {
  empty <- data.frame(layer_key = character(0), field = character(0),
                      regex = logical(0), description = character(0),
                      units = character(0), stringsAsFactors = FALSE)
  if (!file.exists(path)) return(empty)
  ln <- readLines(path, warn = FALSE)
  start <- grep("^##\\s*7\\.\\s", ln)
  if (length(start) == 0) return(empty)
  end <- grep("^##\\s*8\\.\\s", ln)
  end <- if (length(end) > 0) min(end[end > start[1]]) - 1L else length(ln)
  seg <- ln[start[1]:end]

  out <- empty
  keys <- character(0)
  for (l in seg) {
    if (grepl("^###\\s", l)) {
      keys <- unlist(regmatches(l, gregexpr("`[^`]+`", l)))
      keys <- gsub("`", "", keys)
      keys <- keys[grepl("\\.(gpkg|tif|shp)$", keys)]
      next
    }
    if (length(keys) == 0 || !grepl("^\\|", l)) next
    cells <- trimws(strsplit(sub("\\|\\s*$", "", sub("^\\s*\\|", "", l)), "|", fixed = TRUE)[[1]])
    if (length(cells) < 3) next
    if (grepl("^[-: ]+$", cells[1])) next                      # separator row
    if (tolower(cells[1]) %in% c("field", "layer")) next        # header row

    desc  <- cells[2]
    units <- cells[3]
    if (units %in% c("—", "-", "--", "")) units <- ""

    ## "`richness_<Phylum>` (34 marine / 35 terrestrial)" describes a family
    ## of columns, so it becomes one regex entry instead of 69 literals.
    if (grepl("<", cells[1], fixed = TRUE)) {
      stem <- sub("<.*$", "", gsub("`", "", sub("\\s*\\(.*$", "", cells[1])))
      if (nzchar(stem)) {
        out <- rbind(out, data.frame(layer_key = keys, field = paste0("^", stem),
                                     regex = TRUE, description = desc, units = units,
                                     stringsAsFactors = FALSE))
      }
      next
    }

    flds <- gsub("`", "", unlist(regmatches(cells[1], gregexpr("`[^`]+`", cells[1]))))
    if (length(flds) == 0) next
    for (k in keys) {
      out <- rbind(out, data.frame(layer_key = k, field = flds, regex = FALSE,
                                   description = desc, units = units,
                                   stringsAsFactors = FALSE))
    }
  }
  out[!duplicated(out[, c("layer_key", "field")]), , drop = FALSE]
}

FIELD_DICT <- parse_field_dictionaries(PATH_INVENTORY)
if (nrow(FIELD_DICT) == 0) {
  warn_note("No field dictionaries parsed from ", PATH_INVENTORY,
            " — data_dictionary.csv descriptions will be mostly blank.")
}

## The inventory states these layers share another layer's dictionary.
DICT_ALIASES <- c(
  "richness_grid_marine.gpkg"      = "protection_priority_marine.gpkg",
  "richness_grid_terrestrial.gpkg" = "protection_priority_terrestrial.gpkg",
  "hybas_lev12_runoff_risk.gpkg"   = "catchments_with_risk.gpkg",
  "catchment_landcover.gpkg"       = "catchments_with_risk.gpkg",
  "coastal_pour_points_raw.gpkg"   = "coastal_pour_points.gpkg",
  "pour_points.gpkg"               = "rivers_with_accumulated_risk.gpkg"
)

## Fields the inventory lists in section 7.6 without a description column,
## plus GDAL housekeeping. Transcribed from the inventory text, not invented;
## anything not here stays blank and is reported.
## layer scopes an entry: NA applies it to any layer, otherwise it is a regular
## expression matched against the layer key, so a short field name shared with
## another layer cannot pick up a description written for this one.
FIELD_DICT_EXTRA <- data.frame(
  field = c("fid", "geom", "geometry", "area_ha", "habitat", "type", "class",
            "potential", "tier", "connectivity_score", "threat_type",
            "cell_id", "n_indicator", "accum_risk", "accum_risk_raw",
            "seg_length_m", "watercourse_name", "band", "aoi_name",
            ## pass-through attributes from coastalNet, OS Open Rivers and
            ## CEFAS, plus SCAMP fields not yet in the inventory
            "ID",
            "IDSEQ",
            "HEX",
            "status",
            "S_CatNAME",
            "n_classes",
            "total_pixels",
            "dominant_class",
            "dominant_desc",
            "dominant_pct",
            "id",
            "flow_direction",
            "length",
            "fictitious",
            "form",
            "watercourse_name_alternative",
            "start_node",
            "end_node",
            "accum_risk_log",
            "bay_name",
            "local_auth",
            "code",
            "hectares",
            ## fields added by the habitat, indicator and designated-site
            ## rebuilds; all three carry a "label" column, so each entry is
            ## scoped to the layer it belongs to
            "status",
            "label",
            "functional_group",
            "label",
            "designation",
            "site_name",
            "label"),
  description = c(
    "GeoPackage feature identifier",
    "Feature geometry", "Feature geometry",
    "Feature area",
    "Habitat name",
    "Known, Potential, or for sand dune Potential Tier 1 or Potential Tier 2",
    "Habitat by known/potential label",
    "Modelled suitability flag (1 = suitable)",
    "1 = within 2 km of existing dunes; 2 = further away, near dune indicator-species records only",
    "Highest coastalNet probability of larvae arriving from a hexagon with known habitat. Only Potential hexagons are scored (Existing for saltmarsh, which has no potential model); polygons take the score of the hexagon they lie in. Blank means not scored, not zero connectivity",
    "Threat category",
    "Grid cell identifier",
    "Distinct indicator species recorded since 2021 in the 5 km cell and its eight neighbours",
    "Rescaled accumulated diffuse pollution risk at this location",
    "Accumulated diffuse pollution risk at this location, unscaled",
    "Length of the river segment",
    "River name from OS OpenRivers",
    "Environmental band assigned by k-means clustering",
    "Candidate pilot area name",
    "coastalNet hexagon identifier",
    "coastalNet sequential hexagon index",
    "coastalNet H3 hexagon code",
    "Habitat status of the hexagon: Existing, Potential or Neither",
    "Shellfish water catchment name (CEFAS)",
    "Number of distinct land cover classes in the catchment",
    "Land cover cells counted in the catchment",
    "Most extensive land cover class code in the catchment; blank for sub-catchments wholly outside the study area, where land cover was not mapped",
    "Most extensive land cover class, description; blank for sub-catchments wholly outside the study area",
    "Share of the catchment under the most extensive class",
    "Feature identifier from the source dataset",
    "Direction of flow along the segment, from OS Open Rivers",
    "Segment length from OS Open Rivers",
    "OS Open Rivers flag: segment interpolated through a lake or estuary",
    "Watercourse form from OS Open Rivers, e.g. inland river, tidal",
    "Alternative watercourse name from OS Open Rivers",
    "Upstream node identifier from OS Open Rivers",
    "Downstream node identifier from OS Open Rivers",
    "Accumulated diffuse pollution risk, log-transformed before rescaling",
    "Candidate pilot bay name",
    "Local authority name",
    "Local authority code",
    "Area of the local authority polygon",
    "Known extent or modelled potential",
    "Habitat and type combined, e.g. Seagrass - Known; classify on this field",
    "Habitat builders (species that form the habitat) or associated community",
    "Habitat and functional group combined; classify on this field",
    "Designation type: SAC, SPA, SSSI, Ramsar, MPA, NNR, LNR or RSPB",
    "Site name as given by the designating body",
    "Designation type and site name combined"),
  units = c("", "", "", "ha", "", "", "", "1 = suitable", "", "0-1", "",
            "", "count", "0-1", "index", "m", "", "", "",
            "",        # ID
            "",        # IDSEQ
            "",        # HEX
            "",        # status
            "",        # S_CatNAME
            "count",   # n_classes
            "count",   # total_pixels
            "",        # dominant_class
            "",        # dominant_desc
            "%",       # dominant_pct
            "",        # id
            "",        # flow_direction
            "m",       # length
            "",        # fictitious
            "",        # form
            "",        # watercourse_name_alternative
            "",        # start_node
            "",        # end_node
            "",        # accum_risk_log
            "",        # bay_name
            "",        # local_auth
            "",        # code
            "ha",      # hectares
            "",        # status            (all_habitats)
            "",        # label             (all_habitats)
            "",        # functional_group  (indicator_density_by_habitat)
            "",        # label             (indicator_density_by_habitat)
            "",        # designation       (designated_sites)
            "",        # site_name         (designated_sites)
            ""),       # label             (designated_sites)
  ## ID, IDSEQ, HEX and status are coastalNet pass-throughs that appear on all
  ## three connectivity hexagon layers, and status, label, designation and
  ## site_name are generic enough to belong to more than one layer, so all of
  ## them are scoped by pattern rather than by layer name. NA applies anywhere.
  layer = c(rep(NA_character_, 19),
            rep("connectivity_hexagons\\.gpkg$", 4),   # ID, IDSEQ, HEX, status
            rep(NA_character_, 19),
            rep("scamp_all_habitats\\.gpkg$", 2),              # status, label
            rep("indicator_density_by_habitat\\.gpkg$", 2),    # functional_group, label
            rep("designated_sites\\.gpkg$", 3)),              # designation, site_name, label
  stringsAsFactors = FALSE
)
stopifnot(nrow(FIELD_DICT_EXTRA) == 49)

## A regex FIELD_DICT entry covers a family of columns, so the inventory
## description lists the examples instead of describing the column that
## matched. The template names the column, with the part of the field name
## after the regex stem substituted in.
FIELD_DICT_REGEX_DESC <- c("^richness_" = "Species richness for the phylum %s")

lookup_field <- function(layer_key, field) {
  key <- layer_key
  if (key %in% names(DICT_ALIASES)) key <- DICT_ALIASES[[key]]
  d <- FIELD_DICT[FIELD_DICT$layer_key == key, , drop = FALSE]
  hit <- d[!d$regex & tolower(d$field) == tolower(field), , drop = FALSE]
  if (nrow(hit) == 0) {
    rx <- d[d$regex, , drop = FALSE]
    if (nrow(rx) > 0) {
      ok <- vapply(rx$field, function(p) grepl(p, field), logical(1))
      hit <- rx[ok, , drop = FALSE]
      if (nrow(hit) > 0) {
        ## Stems are literal text behind a "^", so the suffix is a substring,
        ## not a sub() that would have to escape the field name.
        stem <- sub("^\\^", "", hit$field[1])
        sfx <- substring(field, nchar(stem) + 1L)
        tpl <- if (hit$field[1] %in% names(FIELD_DICT_REGEX_DESC))
                 FIELD_DICT_REGEX_DESC[[hit$field[1]]] else NA_character_
        hit$description[1] <- if (!is.na(tpl)) sprintf(tpl, sfx)
                              else paste0(sub(" \u2014 .*$", "", hit$description[1]),
                                          " \u2014 ", sfx)
      }
    }
  }
  if (nrow(hit) == 0) {
    ## Drop entries scoped to another layer before matching on field name, so
    ## the rivers layer's lowercase id cannot take the hexagons' ID entry.
    lay <- FIELD_DICT_EXTRA$layer
    in_scope <- is.na(lay) | vapply(seq_along(lay),
                                    function(i) !is.na(lay[i]) && grepl(lay[i], key),
                                    logical(1))
    ex <- FIELD_DICT_EXTRA[in_scope, , drop = FALSE]
    hit <- ex[tolower(ex$field) == tolower(field), , drop = FALSE]
  }
  if (nrow(hit) == 0) return(c(description = "", units = ""))
  c(description = as.character(hit$description[1]), units = as.character(hit$units[1]))
}

##=====================================================================
## 7. SOURCE CATALOGUE
##=====================================================================
## Providers, URLs and attribution for the inputs named in derived_from.
## Licence is stated per dataset where it has been confirmed; anything not
## confirmed falls back to LICENCE_DEFAULT, which tells the recipient to
## check the provider rather than naming a licence we cannot stand behind.
## A licence or README file found on disk overrides the default, and an
## unmatched derived_from entry is reported rather than silently defaulted.
## url = NA means "not confirmed" and lands in the gap report.

src_row <- function(key, aliases, name, provider, url = NA_character_,
                    attribution = NA_character_, licence = NA_character_) {
  data.frame(key = key, aliases = paste(c(key, aliases), collapse = "|"),
             name = name, provider = provider, url = url,
             attribution = attribution, licence = licence,
             stringsAsFactors = FALSE)
}

SOURCE_CATALOGUE <- do.call(rbind, list(
  src_row("gbif", c("global biodiversity information facility"),
          "GBIF occurrence records", "GBIF Secretariat", "https://www.gbif.org"),
  src_row("obis", c("ocean biodiversity information system"),
          "OBIS occurrence records", "OBIS, IOC-UNESCO", "https://obis.org"),
  src_row("nbn", c("nbn atlas", "national biodiversity network"),
          "NBN Atlas occurrence records", "National Biodiversity Network Trust",
          "https://nbnatlas.org"),
  src_row("gems", c("gems", "pmf", "priority marine features", "gems-pmf"),
          "GeMS / Scottish Priority Marine Features", "NatureScot",
          "https://opendata.nature.scot",
          licence = "Open Government Licence v3.0 (NatureScot)"),
  src_row("jncc_inns", c("jncc inns list", "b6 inns"),
          "JNCC B6 invasive non-native species list", "JNCC", "https://jncc.gov.uk",
          licence = "Open Government Licence v3.0 (JNCC)"),
  src_row("hlcm", c("hlcm", "habitat land cover map", "space intelligence"),
          "Habitat and Land Cover Map of Scotland 2022 (EUNIS level 2)",
          "Space Intelligence / NatureScot", "https://opendata.nature.scot",
          licence = "Open Government Licence v3.0 (NatureScot)"),
  src_row("habmos", c("habmos", "habitat map of scotland"),
          "Habitat Map of Scotland (HabMoS)", "NatureScot",
          "https://opendata.nature.scot",
          licence = "Open Government Licence v3.0 (NatureScot)"),
  src_row("sfp_zostera", c("solway firth partnership", "sfp zostera"),
          "Solway Firth Partnership Zostera surveys", "Solway Firth Partnership",
          "https://www.solwayfirthpartnership.co.uk"),
  src_row("wcmc_seagrass", c("wcmc seagrass", "wcmc013", "unep-wcmc seagrass"),
          "Global distribution of seagrasses v7.1", "UNEP-WCMC",
          "https://data.unep-wcmc.org"),
  src_row("wcmc_saltmarsh", c("wcmc saltmarsh", "wcmc027"),
          "Global distribution of saltmarsh v6.1", "UNEP-WCMC",
          "https://data.unep-wcmc.org"),
  src_row("ea_seagrass_potential", c("environment agency seagrass", "seagrass_potential"),
          "Seagrass Potential", "Environment Agency", "https://www.data.gov.uk"),
  src_row("saltmarsh_survey", c("saltmarsh survey of scotland"),
          "Saltmarsh Survey of Scotland", "NatureScot", "https://opendata.nature.scot",
          licence = "Open Government Licence v3.0 (NatureScot)"),
  src_row("historical_oyster", c("historical oyster", "oyster historical"),
          "Historical native oyster habitat locations", "Published historical records"),
  src_row("dynamic_coast", c("dynamic coast", "dc2"),
          "Dynamic Coast 2 (shoreline change)", "Dynamic Coast / Scottish Government",
          "https://www.dynamiccoast.com"),
  src_row("intertidal_foreshore", c("intertidal substrate foreshore"),
          "Intertidal Substrate Foreshore", "Defra / Scottish Government",
          licence = "Open Government Licence v3.0"),
  src_row("cmems", c("copernicus marine", "cmems", "thetao", "bottomt", "ocean colour"),
          "Copernicus Marine Service products (physics, biogeochemistry, ocean colour)",
          "E.U. Copernicus Marine Service", "https://marine.copernicus.eu",
          attribution = "E.U. Copernicus Marine Service Information",
          licence = "Copernicus Marine Service licence (free use with attribution)"),
  src_row("gebco", c("gebco", "gebco 2026"),
          "GEBCO 2026 global bathymetric grid", "GEBCO Compilation Group",
          "https://www.gebco.net",
          attribution = "GEBCO Compilation Group (2026) GEBCO 2026 Grid",
          licence = "Public domain (GEBCO terms of use)"),
  src_row("cop30", c("cop30", "copernicus dem", "copernicus glo-30"),
          "Copernicus DEM GLO-30", "European Space Agency / Copernicus",
          "https://spacedata.copernicus.eu/collections/copernicus-digital-elevation-model",
          licence = "Copernicus DEM licence (free use with attribution)"),
  src_row("euseamap", c("euseamap", "euseamap 2023", "emodnet seabed habitats"),
          "EUSeaMap 2023 (energy class)", "EMODnet Seabed Habitats",
          "https://emodnet.ec.europa.eu/en/seabed-habitats",
          licence = "CC BY 4.0 (EMODnet)"),
  src_row("ukseamap", c("ukseamap", "ukseamap v2026"),
          "UKSeaMap predictive habitat map v2026", "JNCC", "https://jncc.gov.uk",
          licence = "Open Government Licence v3.0 (JNCC)"),
  src_row("bgs_sbs", c("bgs", "predictive seabed sediments", "british geological survey"),
          "Predictive Seabed Sediments UK v1", "British Geological Survey",
          "https://www.bgs.ac.uk"),
  src_row("ke_currents", c("kinetic energy currents", "ke_currents"),
          "Kinetic energy of currents, Atlantic", "Published dataset"),
  src_row("ukcp18", c("ukcp18", "sea level projection"),
          "UKCP18 marine sea level projections (RCP4.5)", "Met Office",
          "https://www.metoffice.gov.uk/research/approach/collaboration/ukcp",
          licence = "Open Government Licence v3.0 (Met Office)"),
  src_row("hydrosheds", c("hydrobasins", "hydrosheds", "hybas"),
          "HydroBASINS level 12 / HydroSHEDS flow direction", "WWF / HydroSHEDS",
          "https://www.hydrosheds.org"),
  src_row("cefas_catchments", c("cefas catchments", "cefas"),
          "CEFAS shellfish water catchments", "CEFAS", "https://www.cefas.co.uk"),
  src_row("os_openrivers", c("os openrivers", "openrivers"),
          "OS Open Rivers", "Ordnance Survey", "https://osdatahub.os.uk/downloads/open",
          attribution = "Contains OS data (c) Crown copyright and database right 2026",
          licence = "Open Government Licence v3.0 (Ordnance Survey)"),
  src_row("os_openmap", c("os openmap local", "openmap local", "opmplc", "os openroads", "openroads"),
          "OS OpenMap Local / OS Open Roads", "Ordnance Survey",
          "https://osdatahub.os.uk/downloads/open",
          attribution = "Contains OS data (c) Crown copyright and database right 2026",
          licence = "Open Government Licence v3.0 (Ordnance Survey)"),
  src_row("hutton", c("james hutton", "hutton", "soil erosion risk",
                      "soil leaching", "soil runoff", "lca 50k"),
          "James Hutton Institute soil and land capability datasets",
          "James Hutton Institute",
          licence = "James Hutton Institute Open Data Licence"),
  src_row("coastalnet", c("coastalnet", "assis"),
          "coastalNet oceanographic connectivity database", "Assis et al. (2025)",
          "https://github.com/jorgeassis/coastalNet",
          licence = "CC BY 4.0 (Assis et al. 2025, Scientific Data)"),
  src_row("cartocean", c("cartocean", "putuhena", "putahena"),
          "CartOcean UK v1 pressure datasets", "Putuhena et al."),
  src_row("emodnet_ha", c("emodnet", "emodnet human activities", "vessel density",
                          "dredge spoil", "dumped munitions", "uwwtd", "ww wrecks",
                          "dredging"),
          "EMODnet Human Activities layers (vessel density, dredging, disposal, wrecks)",
          "EMODnet Human Activities", "https://emodnet.ec.europa.eu/en/human-activities",
          licence = "CC BY 4.0 (EMODnet)"),
  src_row("gfw", c("global fishing watch"),
          "Global Fishing Watch apparent fishing effort", "Global Fishing Watch",
          "https://globalfishingwatch.org"),
  src_row("fishing_u12", c("fishing under 12", "fishing_under12", "inshore fishing"),
          "Scottish inshore fishing activity, vessels under 12 m", "Marine Directorate",
          "https://marine.gov.scot",
          licence = "Open Government Licence v3.0"),
  src_row("macroalgae_models", c("non-native macroalgae", "kuenm"),
          "Non-native macroalgae distribution models", "Published models"),
  src_row("smr", c("scottish marine regions", "marine region"),
          "Scottish Marine Regions", "Marine Directorate", "https://marine.gov.scot",
          licence = "Open Government Licence v3.0"),
  src_row("la_boundaries", c("local authority boundaries", "pub_las"),
          "Scottish local authority boundaries", "Scottish Government",
          "https://spatialdata.gov.scot",
          licence = "Open Government Licence v3.0"),
  src_row("harbour_limits", c("harbour limits"),
          "Harbour authority limits", "Marine Directorate", "https://marine.gov.scot",
          licence = "Open Government Licence v3.0"),
  src_row("designations", c("sac", "spa", "sssi", "ramsar", "mpa", "nnr", "lnr",
                            "designated sites", "nature reserves"),
          "Statutory designated sites (SAC, SPA, SSSI, Ramsar, MPA, NNR, LNR)",
          "NatureScot / JNCC", "https://opendata.nature.scot",
          licence = "Open Government Licence v3.0 (NatureScot)"),
  src_row("rspb", c("rspb reserves"),
          "RSPB public reserves", "RSPB", "https://www.rspb.org.uk"),
  src_row("ibas", c("iba", "important bird areas", "birdlife"),
          "Important Bird and Biodiversity Areas", "BirdLife International",
          "http://datazone.birdlife.org"),
  src_row("worldclim", c("worldclim", "bio1", "bio12"),
          "WorldClim v2.1 bioclimatic variables", "WorldClim",
          "https://www.worldclim.org/data/worldclim21.html",
          attribution = "Fick & Hijmans (2017) WorldClim 2")
))

## Standard attribution statements for open-licence sources whose catalogue
## entry above has none. Only fills blanks; a stated attribution is kept.
ATTRIB_FILL <- c(
  gems             = "Contains NatureScot information licensed under the Open Government Licence v3.0",
  habmos           = "Contains NatureScot information licensed under the Open Government Licence v3.0",
  hlcm             = "Contains NatureScot information licensed under the Open Government Licence v3.0",
  saltmarsh_survey = "Contains NatureScot information licensed under the Open Government Licence v3.0",
  designations     = "Contains NatureScot information licensed under the Open Government Licence v3.0",
  jncc_inns        = "Contains JNCC data (c) copyright and database right 2026, licensed under the Open Government Licence v3.0",
  ukseamap         = "Contains JNCC data (c) copyright and database right 2026, licensed under the Open Government Licence v3.0",
  smr              = "Contains public sector information licensed under the Open Government Licence v3.0",
  harbour_limits   = "Contains public sector information licensed under the Open Government Licence v3.0",
  fishing_u12      = "Contains public sector information licensed under the Open Government Licence v3.0",
  intertidal_foreshore = "Contains public sector information licensed under the Open Government Licence v3.0",
  la_boundaries    = "Contains public sector information licensed under the Open Government Licence v3.0",
  ukcp18           = "(c) Crown copyright 2026, Met Office, licensed under the Open Government Licence v3.0",
  euseamap         = "EMODnet Seabed Habitats, EUSeaMap 2023, licensed under CC BY 4.0",
  emodnet_ha       = "EMODnet Human Activities, licensed under CC BY 4.0",
  coastalnet       = "Assis et al. (2025) Scientific Data 12, 737, licensed under CC BY 4.0",
  cop30            = "(c) DLR e.V. 2010-2014 and (c) Airbus Defence and Space GmbH 2014-2018 provided under COPERNICUS by the European Union and ESA; all rights reserved"
)
for (k in names(ATTRIB_FILL)) {
  i <- SOURCE_CATALOGUE$key == k &
       (is.na(SOURCE_CATALOGUE$attribution) | !nzchar(SOURCE_CATALOGUE$attribution))
  SOURCE_CATALOGUE$attribution[i] <- ATTRIB_FILL[[k]]
}

## Licences, attributions and links confirmed by the originator on 2 October
## 2026. These override the catalogue entries above. NA leaves a value as is.
SOURCE_CONFIRMED <- list(
  sfp_zostera = list(attribution = "Solway Firth Partnership data"),
  wcmc_seagrass = list(
    licence = "CC BY 4.0",
    attribution = "UNEP-WCMC (2026) Ocean+ Habitats. Available at: https://habitats.oceanplus.org",
    url = "https://habitats.oceanplus.org"),
  historical_oyster = list(
    provider = "Thurstan et al. (2024), Nature Sustainability",
    licence = "CC BY 4.0",
    attribution = "Thurstan, R.H., McCormick, H., Preston, J. et al. (2024) Records reveal the vast historical extent of European oyster reef ecosystems. Nature Sustainability 7, 1719-1729",
    url = "https://doi.org/10.1038/s41893-024-01441-4"),
  bgs_sbs = list(
    licence = "Open Government Licence v3.0",
    attribution = "Contains British Geological Survey materials (c) UKRI 2026"),
  ke_currents = list(
    provider = "EMODnet Seabed Habitats",
    licence = "CC BY 4.0 (EMODnet)",
    attribution = "EMODnet Seabed Habitats, licensed under CC BY 4.0",
    url = "https://emodnet.ec.europa.eu/en/seabed-habitats"),
  dynamic_coast = list(
    licence = "Provider terms: Dynamic Coast analyses must not be used for property-level assessments",
    attribution = "Rennie, A.F. et al. (2021) Dynamic Coast: The National Overview (2021). CREW, Edinburgh",
    url = "https://www.crew.ac.uk/sites/www.crew.ac.uk/files/publication/CREW_DC2_SYNOPSIS_FINAL.pdf"),
  intertidal_foreshore = list(
    url = "https://www.data.gov.uk/dataset/6efcebae-874e-4691-bf46-53057bdebda1/intertidal-substrate-foreshore-england-and-scotland"),
  hutton = list(
    attribution = "Soil Survey of Scotland Staff (1984-1987) Land Capability for Agriculture maps of Scotland at a scale of 1:50 000. Macaulay Institute for Soil Research, Aberdeen",
    url = "https://doi.org/10.5281/zenodo.6322760"),
  gfw = list(
    licence = "CC BY-NC 4.0",
    attribution = "Global Fishing Watch apparent fishing effort, licensed under CC BY-NC 4.0"),
  cartocean = list(
    licence = "CC BY 4.0",
    attribution = "Putuhena, H., Williams, T.J., Sturt, F. et al. (2025) Integrated geospatial datasets to inform marine spatial planning and impact assessment in waters surrounding the United Kingdom. Scientific Data 12, 1845",
    url = "https://doi.org/10.1038/s41597-025-05950-5"),
  hydrosheds = list(
    licence = "CC BY 4.0",
    attribution = "Lehner, B. & Grill, G. (2013) Global river hydrography and network routing: baseline data and new approaches to study the world's large river systems. Hydrological Processes 27, 2171-2186"),
  worldclim = list(licence = "CC BY 4.0"),
  rspb = list(
    licence = "CC BY 4.0",
    attribution = "Contains public sector information licensed under the CC BY 4.0. Copyright (c) Royal Society for the Protection of Birds (RSPB)"),
  ibas = list(
    licence = "BirdLife International terms: non-commercial use only; no redistribution. Not included in this package")
)
for (k in names(SOURCE_CONFIRMED)) {
  i <- SOURCE_CATALOGUE$key == k
  for (fld in names(SOURCE_CONFIRMED[[k]])) {
    SOURCE_CATALOGUE[[fld]][i] <- SOURCE_CONFIRMED[[k]][[fld]]
  }
}

LICENCE_DEFAULT <- "Not confirmed: see the provider's terms of use"

## Licence and README files on disk override the default. Found by scan so
## a newly downloaded dataset's terms are picked up without a code edit.
find_licence_files <- function(roots = "data/raw") {
  pat <- "licen|readme|terms|copyright"
  out <- character(0)
  for (rt in roots) {
    if (!dir.exists(rt)) next
    f <- list.files(rt, recursive = TRUE, full.names = TRUE, all.files = FALSE)
    out <- c(out, f[grepl(pat, basename(f), ignore.case = TRUE)])
  }
  fwd(out)
}
LICENCE_FILES <- find_licence_files()

## Match a licence file to a catalogue entry by the directory it sits in.
licence_file_for <- function(key) {
  if (length(LICENCE_FILES) == 0) return(NA_character_)
  al <- strsplit(SOURCE_CATALOGUE$aliases[SOURCE_CATALOGUE$key == key], "|", fixed = TRUE)[[1]]
  al <- gsub("[^a-z0-9]+", "", tolower(al))
  al <- al[nchar(al) >= 4]
  if (length(al) == 0) return(NA_character_)
  norm <- gsub("[^a-z0-9]+", "", tolower(LICENCE_FILES))
  hit <- LICENCE_FILES[vapply(norm, function(n) any(vapply(al, function(a) grepl(a, n, fixed = TRUE), logical(1))), logical(1))]
  if (length(hit) == 0) NA_character_ else paste(hit, collapse = "; ")
}

norm_key <- function(x) gsub("[^a-z0-9]+", " ", tolower(trimws(x)))

## Whole-word matching on space-padded strings, so a short alias like "sac"
## cannot match a longer word that happens to contain it.
match_source <- function(term) {
  t <- paste0(" ", norm_key(term), " ")
  if (!nzchar(trimws(t))) return(NA_character_)
  for (i in seq_len(nrow(SOURCE_CATALOGUE))) {
    al <- norm_key(strsplit(SOURCE_CATALOGUE$aliases[i], "|", fixed = TRUE)[[1]])
    al <- paste0(" ", al[nzchar(al)], " ")
    if (any(vapply(al, function(a) grepl(a, t, fixed = TRUE) || grepl(t, a, fixed = TRUE),
                   logical(1)))) {
      return(SOURCE_CATALOGUE$key[i])
    }
  }
  NA_character_
}

##=====================================================================
## 8. OCCURRENCE PROVENANCE
##=====================================================================
## Traced by reading 02_themeD_species.R: query_nbn() returns early for the
## bulk query (no species list), so the richness / beta / similarity /
## protection-priority chain is GBIF + OBIS only. The named-species queries
## for INNS and indicators do call NBN. The trace is verified against the
## scripts at runtime, and against the data where a source field exists,
## because an unreported record source is a licensing problem as well as a
## provenance one.

OCCURRENCE_TRACE <- data.frame(
  pattern = c("richness_grid", "rarefied_richness", "beta_diversity",
              "community_similarity", "protection_priority",
              "indicator_density", "indicator_records", "inns_records",
              "sand_dune_restoration_score"),
  sources = c("GBIF;OBIS", "GBIF;OBIS", "GBIF;OBIS", "GBIF;OBIS", "GBIF;OBIS",
              "GBIF;OBIS;NBN", "GBIF;OBIS;NBN", "GBIF;OBIS;NBN",
              "GBIF;OBIS;NBN"),
  evidence = c(
    "02_themeD_species.R: bulk GBIF + OBIS download; query_nbn() skips bulk queries (no species list)",
    "07_plots.R: IDW interpolation of the 5 km richness grid (GBIF + OBIS)",
    "07_plots.R: IDW interpolation of grid beta diversity (GBIF + OBIS)",
    "07_plots.R: IDW interpolation of grid community similarity (GBIF + OBIS)",
    "07_plots.R: derived from the 5 km richness grid plus designations (GBIF + OBIS)",
    "02_themeD_species.R: indicator_records = GBIF + OBIS + query_nbn(indicator_names)",
    "02_themeD_species.R: bind_rows(ind_gbif, ind_obis, ind_nbn)",
    "02_themeD_species.R: bind_rows(inns_gbif, inns_obis, inns_nbn)",
    "07_plots.R: habitat-builder density component uses indicator_records (GBIF + OBIS + NBN)"),
  stringsAsFactors = FALSE
)

## Confirm the trace still describes the code, rather than trusting a
## comment written against an earlier version of it.
verify_trace <- function() {
  notes <- character(0)
  sp <- "02_themeD_species.R"
  if (!file.exists(sp)) {
    return("02_themeD_species.R not found — the occurrence trace could not be re-verified.")
  }
  ln <- readLines(sp, warn = FALSE)
  if (!any(grepl("query_nbn", ln, fixed = TRUE))) {
    notes <- c(notes, "query_nbn() no longer appears in 02_themeD_species.R — re-trace NBN use.")
  }
  if (!any(grepl("Skipping NBN for bulk query", ln, fixed = TRUE))) {
    notes <- c(notes, paste("The bulk-query NBN skip could not be confirmed in",
                            "02_themeD_species.R. If bulk NBN downloads are now enabled,",
                            "the richness / beta / similarity / protection-priority layers",
                            "include NBN records and the trace below is wrong."))
  }
  for (v in c("inns_nbn", "ind_nbn")) {
    if (!any(grepl(v, ln, fixed = TRUE))) {
      notes <- c(notes, sprintf("%s no longer appears in 02_themeD_species.R — re-trace NBN use.", v))
    }
  }
  notes
}

## Distinct values of a source field, read with OGR SQL so no geometry or
## full attribute table is loaded.
distinct_sources <- function(path, layer) {
  q <- sprintf('SELECT DISTINCT "source" FROM "%s"', layer)
  v <- try(suppressWarnings(sf::st_read(path, query = q, quiet = TRUE)), silent = TRUE)
  if (inherits(v, "try-error") || is.null(v) || nrow(v) == 0) return(character(0))
  if (inherits(v, "sf")) v <- sf::st_drop_geometry(v)
  sort(unique(stats::na.omit(as.character(v[[1]]))))
}

traced_sources_for <- function(layer_name, file) {
  hit <- OCCURRENCE_TRACE[vapply(OCCURRENCE_TRACE$pattern,
                                 function(p) grepl(p, layer_name, fixed = TRUE) ||
                                             grepl(p, basename(file), fixed = TRUE),
                                 logical(1)), , drop = FALSE]
  if (nrow(hit) == 0) return(NULL)
  hit[1, , drop = FALSE]
}

##=====================================================================
## 9. STAGE 1 — AUDIT
##=====================================================================

if (STAGE == 1) {

  dir.create(HANDOVER_DIR, recursive = TRUE, showWarnings = FALSE)
  rpt <- character(0)
  add <- function(...) rpt <<- c(rpt, paste0(...))

  add("# SCAMP handover — gap report")
  add("")
  add(sprintf("Generated %s by `scripts/09_handover_package.R` (stage 1).", format(Sys.time(), "%Y-%m-%d %H:%M")))
  add(sprintf("Manifest: `%s` — %d rows, %d marked include.",
              fwd(manifest_src), nrow(manifest), sum(manifest$include_flag)))
  add("")
  add("Stage 1 writes nothing into the package. Work through the sections below,")
  add("amend the manifest, and re-run until (a), (c), (e) and (f) are clear.")
  add("")

  ##---- (a) missing paths ------------------------------------------------
  add("## (a) Included rows whose path is missing")
  add("")
  all_candidates <- unlist(lapply(SCAN_ROOTS, function(rt) {
    if (!dir.exists(rt)) character(0) else fwd(list.files(rt, recursive = TRUE, full.names = TRUE))
  }))
  cand_base <- basename(all_candidates)

  nearest <- function(target, n = 3) {
    if (length(all_candidates) == 0) return(character(0))
    d <- as.numeric(utils::adist(tolower(basename(target)), tolower(cand_base),
                                 costs = list(ins = 1, del = 1, sub = 1)))
    all_candidates[order(d)][seq_len(min(n, length(all_candidates)))]
  }

  inc_rows <- manifest[manifest$include_flag, , drop = FALSE]
  missing_any <- FALSE
  for (i in seq_len(nrow(inc_rows))) {
    r <- inc_rows[i, ]
    files <- expanded$file[expanded$row_id == r$row_id]
    if (length(files) > 0 && all(file.exists(files))) next
    missing_any <- TRUE
    add(sprintf("- **%s** — `%s`%s", r$layer_name, fwd(r$path),
                if (grepl("*", r$path, fixed = TRUE)) " (wildcard matched nothing)" else " not found"))
    for (s in nearest(r$path)) add(sprintf("    - closest on disk: `%s`", s))
  }
  if (!missing_any) add("None. Every included path resolves to a file on disk.")
  add("")

  ##---- (b) wildcard expansion ------------------------------------------
  add("## (b) Wildcard rows and what they expand to")
  add("")
  wc <- manifest[grepl("*", manifest$path, fixed = TRUE), , drop = FALSE]
  if (nrow(wc) == 0) {
    add("No wildcard paths in the manifest.")
  } else {
    for (i in seq_len(nrow(wc))) {
      r <- wc[i, ]
      e <- expanded[expanded$row_id == r$row_id, , drop = FALSE]
      add(sprintf("**%s** (`%s`, include = %s) — %d file(s)",
                  r$layer_name, fwd(r$path), r$include, nrow(e)))
      if (nrow(e) == 0) {
        add("  - *nothing matched*")
      } else {
        for (j in seq_len(nrow(e))) {
          add(sprintf("  - `%s` -> published as **%s** (%s MB)",
                      e$file[j], e$layer_name[j], format(size_mb(e$file[j]), nsmall = 2)))
        }
      }
      add("")
    }
    add("Published names for wildcard rows are generated from the text the `*`")
    add("matched. Put a `*` in `layer_name` to control where that text lands.")
  }
  add("")

  ##---- (c) uncovered spatial files -------------------------------------
  add("## (c) Spatial files under outputs/ that no manifest row covers")
  add("")
  outputs_files <- character(0)
  for (rt in c("outputs", DESIGN_DIR)) {
    if (!dir.exists(rt)) next
    f <- fwd(list.files(rt, recursive = TRUE, full.names = TRUE))
    outputs_files <- c(outputs_files, f[ext_of(f) %in% SPATIAL_EXT])
  }
  outputs_files <- unique(outputs_files)
  outputs_files <- outputs_files[!grepl(sprintf("^%s/", PKG_DIR), outputs_files)]
  covered <- unique(fwd(expanded$file))
  orphans <- setdiff(outputs_files, covered)
  if (!dir.exists(DESIGN_DIR)) {
    add(sprintf("> `%s` does not exist yet, so no cube or environmental band files", DESIGN_DIR))
    add("> could be checked. Run `scripts/08_environmental_bands.R` first.")
    add("")
  }
  if (length(orphans) == 0) {
    add("None. Every spatial file under `outputs/` appears in the manifest.")
  } else {
    add(sprintf("%d file(s) are not referenced by any manifest row, included or", length(orphans)))
    add("excluded. Each needs a row — if only to record why it is not shipped.")
    add("")
    for (o in sort(orphans)) add(sprintf("- `%s` (%s MB)", o, format(size_mb(o), nsmall = 2)))
  }
  add("")

  ##---- (d) actual vs declared ------------------------------------------
  add("## (d) Included layers: what is actually on disk")
  add("")
  add("| layer_name | type | CRS (EPSG) | resolution / features | size MB | flags |")
  add("|---|---|---|---|---|---|")
  probe_cache <- list()
  for (i in seq_len(nrow(inc))) {
    row <- inc[i, ]
    if (!file.exists(row$file)) next
    flags <- character(0)
    crs_actual <- NA_integer_
    detail <- ""
    if (row$type == "raster") {
      p <- probe_raster(row$file)
      probe_cache[[row$file]] <- p
      if (!p$ok) {
        flags <- c(flags, paste0("UNREADABLE: ", p$msg))
      } else {
        crs_actual <- p$epsg
        detail <- sprintf("%s m, %d band(s)", format(signif(p$res, 6)), p$nlyr)
        dec_res <- suppressWarnings(as.numeric(row$resolution_m))
        if (!is.na(dec_res) && !isTRUE(all.equal(dec_res, p$res, tolerance = 0.02))) {
          flags <- c(flags, sprintf("resolution_m says %s, raster is %s",
                                    format(dec_res), format(signif(p$res, 6))))
        }
        if (isTRUE(p$lonlat)) flags <- c(flags, "geographic CRS — will not be reprojected (rasters ship at native CRS)")
      }
    } else if (row$type == "vector") {
      p <- probe_vector(row$file)
      probe_cache[[row$file]] <- p
      if (!p$ok) {
        flags <- c(flags, paste0("UNREADABLE: ", p$msg))
      } else {
        crs_actual <- p$epsg[1]
        detail <- paste(sprintf("%s: %d feat (%s)", p$layers, p$features, p$geomtype),
                        collapse = "; ")
        if (length(p$layers) > 1) {
          flags <- c(flags, sprintf("%d layers — ships as %s_<layer>", length(p$layers), row$layer_name))
        }
        if (any(is.na(p$epsg))) flags <- c(flags, "undefined CRS — cannot reproject")
        if (any(!is.na(p$epsg) & p$epsg != TARGET_EPSG)) {
          flags <- c(flags, sprintf("EPSG:%s — will be reprojected to %d",
                                    paste(unique(p$epsg[!is.na(p$epsg)]), collapse = "/"), TARGET_EPSG))
        }
      }
    } else {
      detail <- "table"
    }
    dec_crs <- suppressWarnings(as.integer(row$crs_epsg))
    if (!is.na(dec_crs) && !is.na(crs_actual) && dec_crs != crs_actual) {
      flags <- c(flags, sprintf("crs_epsg says %d, file is %d", dec_crs, crs_actual))
    }
    add(sprintf("| %s | %s | %s | %s | %s | %s |",
                row$layer_name, row$type,
                ifelse(is.na(crs_actual), "?", crs_actual),
                detail, format(size_mb(row$file), nsmall = 2),
                ifelse(length(flags) == 0, "", paste(flags, collapse = "; "))))
  }
  add("")

  ##---- (e) blank description / units -----------------------------------
  add("## (e) Included rows with an empty description or units")
  add("")
  need <- inc_rows[blank(inc_rows$description) | blank(inc_rows$units), , drop = FALSE]
  if (nrow(need) == 0) {
    add("None.")
  } else {
    add("Stage 2 stops on an empty `description`. Empty `units` is only a warning,")
    add("but a unitless index should still say so.")
    add("")
    for (i in seq_len(nrow(need))) {
      what <- c(if (blank(need$description[i])) "description",
                if (blank(need$units[i])) "units")
      add(sprintf("- **%s** (`%s`) — missing %s",
                  need$layer_name[i], fwd(need$path[i]), paste(what, collapse = " and ")))
    }
  }
  add("")

  ##---- (f) occurrence provenance ---------------------------------------
  add("## (f) Occurrence provenance of the biodiversity layers")
  add("")
  declared <- names(meta_get("occurrence_sources", default = list()))
  if (length(declared) == 0) declared <- character(0)
  add(sprintf("`handover_meta$occurrence_sources` declares: **%s**",
              if (length(declared)) paste(declared, collapse = ", ") else "nothing"))
  add("")
  trace_notes <- verify_trace()
  for (n in trace_notes) add(sprintf("> %s", n))
  if (length(trace_notes) > 0) add("")
  add("| layer | sources (traced from code) | sources (from the `source` field) | evidence |")
  add("|---|---|---|---|")
  observed_all <- character(0)
  traced_all <- character(0)
  for (i in seq_len(nrow(inc))) {
    row <- inc[i, ]
    tr <- traced_sources_for(row$layer_name, row$file)
    if (is.null(tr)) next
    traced <- split_semi(tr$sources)
    traced_all <- c(traced_all, traced)
    obs <- character(0)
    if (row$type == "vector" && file.exists(row$file)) {
      pv <- probe_cache[[row$file]]
      if (is.null(pv)) pv <- probe_vector(row$file)
      if (isTRUE(pv$ok)) for (L in pv$layers) obs <- c(obs, distinct_sources(row$file, L))
    }
    obs <- sort(unique(obs))
    observed_all <- c(observed_all, obs)
    add(sprintf("| %s | %s | %s | %s |", row$layer_name,
                paste(traced, collapse = ", "),
                if (length(obs)) paste(obs, collapse = ", ") else "no `source` field",
                tr$evidence))
  }
  add("")
  found <- sort(unique(toupper(c(traced_all, observed_all))))
  extra <- setdiff(found, toupper(declared))
  if (length(extra) > 0) {
    add(sprintf("**Mismatch: %s appears in the shipped layers but is not declared in `handover_meta$occurrence_sources`.**",
                paste(extra, collapse = ", ")))
    add("")
    if ("NBN" %in% extra) {
      add("NBN Atlas records reach the indicator and INNS layers through the")
      add("named-species queries. Either add NBN to `occurrence_sources` with its")
      add("licence and attribution, or drop those layers from the manifest. Do not")
      add("edit the lineage text to remove the mention.")
      add("")
    }
  } else if (length(found) > 0) {
    add("No mismatch: every source traced into a shipped layer is declared.")
    add("")
  }
  unused <- setdiff(toupper(declared), found)
  if (length(unused) > 0) {
    add(sprintf("Declared but not traced into any shipped layer: %s.",
                paste(unused, collapse = ", ")))
    add("")
  }

  ##---- (g) size --------------------------------------------------------
  add("## (g) Estimated package size")
  add("")
  inc_files <- unique(inc$file[file.exists(inc$file)])
  cube_files <- if (dir.exists(DESIGN_DIR)) {
    c(Sys.glob(file.path(DESIGN_DIR, "SCAMP_*_cube_*.tif")),
      Sys.glob(file.path(DESIGN_DIR, "SCAMP_*_cube_*_*.csv")))
  } else character(0)
  est <- size_mb(c(inc_files, setdiff(fwd(cube_files), inc_files)))
  add(sprintf("- Included manifest files: %d, %s MB", length(inc_files),
              format(size_mb(inc_files), nsmall = 1)))
  add(sprintf("- Cubes and sidecars from `%s`: %d file(s), %s MB", DESIGN_DIR,
              length(cube_files), format(size_mb(cube_files), nsmall = 1)))
  add(sprintf("- **Estimated total: %s MB** (COG compression and dropped fields will", format(est, nsmall = 1)))
  add("  move this; GeoPackage overhead is small.")
  add("")
  if (est > SIZE_WARN_MB) {
    add(sprintf("> **Over %d MB.** Check whether the largest rasters are all needed at", SIZE_WARN_MB))
    add("> native resolution, or whether the package should be split by theme.")
    add("")
  }
  big <- inc_files[order(-file.size(inc_files))]
  big <- big[seq_len(min(10, length(big)))]
  if (length(big) > 0) {
    add("Largest included files:")
    add("")
    for (b in big) add(sprintf("- `%s` — %s MB", b, format(size_mb(b), nsmall = 1)))
    add("")
  }

  ##---- figures referenced by the documentation -------------------------
  ## README.md and METHODS.md are shipped verbatim, so a figure they link to
  ## but that is not on disk becomes a broken image in the partner's copy.
  add("## Figures referenced by README.md and METHODS.md")
  add("")
  fig_on_disk <- basename(fwd(list.files(HANDOVER_DIR, pattern = "\\.png$",
                                         full.names = TRUE, recursive = FALSE,
                                         ignore.case = TRUE)))
  fig_refs <- character(0)
  for (doc in c(PATH_README, PATH_METHODS)) {
    txt <- readLines(doc, warn = FALSE)
    hits <- unlist(regmatches(txt, gregexpr("figures/[^\\s\\)\\]\"'`]+", txt, perl = TRUE)))
    if (length(hits) > 0) {
      fig_refs <- c(fig_refs, sprintf("%s|%s", basename(doc), hits))
    }
  }
  if (length(fig_refs) == 0) {
    add("Neither document references a `figures/` path. Nothing to check.")
  } else {
    missing_figs <- character(0)
    for (fr in unique(fig_refs)) {
      parts <- strsplit(fr, "|", fixed = TRUE)[[1]]
      want <- basename(parts[2])
      if (!tolower(want) %in% tolower(fig_on_disk)) {
        missing_figs <- c(missing_figs, sprintf("- `%s` referenced by **%s** — not found in `%s/`",
                                                parts[2], parts[1], HANDOVER_DIR))
      }
    }
    if (length(missing_figs) == 0) {
      add(sprintf("All %d referenced figure(s) are present in `%s/` and will be copied to `figures/`.",
                  length(unique(fig_refs)), HANDOVER_DIR))
    } else {
      add(sprintf("%d referenced figure(s) are missing. Stage 2 copies .png files from", length(missing_figs)))
      add(sprintf("`%s/` into the package's `figures/` folder; anything not there ships", HANDOVER_DIR))
      add("as a broken image link.")
      add("")
      for (mf in missing_figs) add(mf)
    }
    if (length(fig_on_disk) > 0) {
      add("")
      add(sprintf("Figures on disk ready to ship: %s", paste(sprintf("`%s`", fig_on_disk), collapse = ", ")))
    }
  }
  add("")

  ##---- config TODOs ----------------------------------------------------
  add("## Metadata config")
  add("")
  if (length(TODOS) == 0) {
    add(sprintf("No TODO placeholders left in `%s`.", PATH_METADATA))
  } else {
    add(sprintf("%d TODO placeholder(s) in `%s`:", length(TODOS), PATH_METADATA))
    add("")
    for (t in TODOS) add(sprintf("- %s", t))
  }
  add("")

  ##---- anything else worth knowing -------------------------------------
  if (length(WARNINGS) > 0) {
    add("## Other warnings raised while auditing")
    add("")
    for (w in WARNINGS) add(sprintf("- %s", w))
    add("")
  }

  out_rpt <- file.path(HANDOVER_DIR, "gap_report.md")
  writeLines(rpt, out_rpt)
  cat(sprintf("\nStage 1 complete. Report written: %s\n", out_rpt))
  cat("Nothing was written into the package. Review the report, amend the\n")
  cat("manifest, then set STAGE <- 2.\n")
}

##=====================================================================
## 10. STAGE 2 — BUILD
##=====================================================================

if (STAGE == 2) {

  if (!requireNamespace("geometa", quietly = TRUE)) {
    stop_clear("The geometa package is required to write the GEMINI 2.3 metadata\n",
               "record and is not installed. install.packages(\"geometa\")")
  }

  ##---- 10.1 pre-flight checks -----------------------------------------
  cat("\n-- Pre-flight --\n")

  missing_paths <- inc$file[!file.exists(inc$file)]
  zero_rows <- manifest[manifest$include_flag, , drop = FALSE]
  zero_rows <- zero_rows[!zero_rows$row_id %in% expanded$row_id[file.exists(expanded$file)], , drop = FALSE]
  if (length(missing_paths) > 0 || nrow(zero_rows) > 0) {
    stop_clear("Stage 2 stopped: included manifest paths are missing.\n",
               if (nrow(zero_rows) > 0) paste0(
                 "Rows resolving to no file:\n",
                 paste0("  - ", zero_rows$layer_name, " -> ", fwd(zero_rows$path), collapse = "\n"), "\n") else "",
               if (length(missing_paths) > 0) paste0(
                 "Files not found:\n",
                 paste0("  - ", missing_paths, collapse = "\n"), "\n") else "",
               "\nRun stage 1 for the closest matches on disk.")
  }

  no_desc <- inc[blank(inc$description), , drop = FALSE]
  if (nrow(no_desc) > 0) {
    stop_clear("Stage 2 stopped: included layers with an empty description.\n",
               "A partner cannot use a layer nobody has described.\n\n",
               paste0("  - ", no_desc$layer_name, " (", fwd(no_desc$file), ")", collapse = "\n"))
  }

  no_units <- inc[blank(inc$units) & inc$type != "vector", , drop = FALSE]
  if (nrow(no_units) > 0) {
    warn_note(nrow(no_units), " raster/table layer(s) have no units: ",
              paste(no_units$layer_name, collapse = ", "))
  }
  for (t in TODOS) warn_note("Unresolved config placeholder — ", t)

  cube_tifs <- fwd(Sys.glob(file.path(DESIGN_DIR, "SCAMP_*_cube_*.tif")))
  if (length(cube_tifs) == 0) {
    stop_clear("Stage 2 stopped: no cubes found at\n  ",
               file.path(DESIGN_DIR, "SCAMP_*_cube_*.tif"), "\n\n",
               "README.md documents rasters/cubes/ as part of the package. Run\n",
               "scripts/08_environmental_bands.R (both passes) first.")
  }

  ##---- 10.2 folders ---------------------------------------------------
  if (dir.exists(PKG_DIR)) unlink(PKG_DIR, recursive = TRUE)
  for (d in c("", "metadata", "figures", "vectors", "rasters/native", "rasters/cubes", "tables")) {
    dir.create(file.path(PKG_DIR, d), recursive = TRUE, showWarnings = FALSE)
  }
  cat(sprintf("Building into %s\n", PKG_DIR))

  ## One row per shipped file, accumulated as the build proceeds, so
  ## layer_register.csv describes what was actually written.
  REGISTER <- data.frame(layer_name = character(0), file = character(0),
                         theme = character(0), type = character(0),
                         resolution_m = character(0), crs_epsg = character(0),
                         units = character(0), description = character(0),
                         derived_from = character(0), stringsAsFactors = FALSE)
  reg_add <- function(layer_name, file, theme, type, resolution_m, crs_epsg,
                      units, description, derived_from) {
    REGISTER <<- rbind(REGISTER, data.frame(
      layer_name = layer_name, file = fwd(file), theme = theme, type = type,
      resolution_m = as.character(resolution_m), crs_epsg = as.character(crs_epsg),
      units = units, description = description, derived_from = derived_from,
      stringsAsFactors = FALSE))
  }

  ##---- 10.3 vectors ---------------------------------------------------
  cat("\n-- Vectors --\n")
  GPKG <- file.path(PKG_DIR, "vectors", "SCAMP_vectors.gpkg")
  vec_rows <- inc[inc$type == "vector", , drop = FALSE]
  written_layers <- character(0)

  for (i in seq_len(nrow(vec_rows))) {
    row <- vec_rows[i, ]
    lyrs <- sf::st_layers(row$file)
    multi <- length(lyrs$name) > 1
    for (L in lyrs$name) {
      out_name <- if (multi) paste0(row$layer_name, "_", L) else row$layer_name
      x <- sf::st_read(row$file, layer = L, quiet = TRUE)

      ## Check, do not assume: several project layers are EPSG:4326.
      cur <- sf::st_crs(x)
      if (is.na(cur)) {
        warn_note("Skipped ", out_name, ": source has no CRS defined (", row$file, ")")
        next
      }
      if (is.na(cur$epsg) || cur$epsg != TARGET_EPSG) {
        x <- sf::st_transform(x, TARGET_EPSG)
        cat(sprintf("  %s: reprojected from %s to EPSG:%d\n", out_name,
                    ifelse(is.na(cur$epsg), "a non-EPSG CRS", paste0("EPSG:", cur$epsg)),
                    TARGET_EPSG))
      }

      ## Field names in the manifest may not match the source's capitalisation
      ## (inspireid vs InspireID), so match case-insensitively.
      drop <- split_semi(row$drop_fields)
      if (length(drop) > 0) {
        gcol <- attr(x, "sf_column")
        drop_lc <- tolower(drop)
        keep <- names(x)[!(tolower(names(x)) %in% drop_lc) | names(x) == gcol]
        absent <- drop[!(drop_lc %in% tolower(names(x)))]
        if (length(absent) > 0) {
          warn_note(out_name, ": drop_fields names field(s) not present — ",
                    paste(absent, collapse = ", "))
        }
        x <- x[, union(keep, gcol), drop = FALSE]
      }

      ## BirdLife's terms do not allow Important Bird Area boundaries to be
      ## passed on, so they are removed from the designated sites layer.
      if (out_name == "bnd_designated_sites" && "designation" %in% names(x)) {
        n0 <- nrow(x)
        x <- x[is.na(x$designation) | x$designation != "IBA", ]
        if (nrow(x) < n0) cat(sprintf("  %s: %d Important Bird Area feature(s) removed (BirdLife terms)\n",
                                      out_name, n0 - nrow(x)))
      }

      sf::st_write(x, GPKG, layer = out_name, append = FALSE, quiet = TRUE)
      written_layers <- c(written_layers, out_name)
      cat(sprintf("  %s (%d features, %d fields)\n", out_name, nrow(x),
                  ncol(x) - 1L))

      reg_add(out_name, "vectors/SCAMP_vectors.gpkg", row$theme, "vector",
              "", TARGET_EPSG, row$units,
              if (multi) sprintf("%s (source layer: %s)", row$description, L) else row$description,
              row$derived_from)
      rm(x)
    }
  }
  if (length(written_layers) == 0) warn_note("No vector layers were written.")

  ##---- 10.4 native rasters --------------------------------------------
  cat("\n-- Native rasters --\n")

  ## COG options: DEFLATE with a horizontal predictor, 3 for floating point
  ## and 2 for integers (predictor 3 on integer data makes files larger).
  pick_datatype <- function(dt) {
    order <- c("INT1U", "INT2U", "INT2S", "INT4U", "INT4S", "FLT4S", "FLT8S")
    dt <- dt[dt %in% order]
    if (length(dt) == 0) return(NA_character_)
    order[max(match(dt, order))]
  }
  default_nodata <- function(dt) {
    if (is.na(dt)) return(NA)
    switch(dt, INT1U = 255, INT2U = 65535, INT4U = 4294967295,
           INT2S = -32768, INT4S = -2147483648, -9999)
  }

  ras_rows <- inc[inc$type == "raster", , drop = FALSE]
  for (i in seq_len(nrow(ras_rows))) {
    row <- ras_rows[i, ]
    r <- terra::rast(row$file)                     # lazy: header only
    dt <- pick_datatype(terra::datatype(r))
    is_float <- !is.na(dt) && grepl("^FLT", dt)
    nd <- suppressWarnings(terra::NAflag(r))
    nd <- if (length(nd) == 0 || all(is.na(nd)) || all(is.nan(nd))) default_nodata(dt) else nd[1]

    theme_dir <- gsub("[^A-Za-z0-9_-]+", "_", row$theme)
    if (!nzchar(theme_dir)) theme_dir <- "other"
    out_dir <- file.path(PKG_DIR, "rasters", "native", theme_dir)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    out_file <- file.path(out_dir, paste0(row$layer_name, ".tif"))

    ## Band names are carried on the object; setting them does not read cells.
    ## A single-band raster takes the layer name, so GIS software does not show
    ## a leftover processing name (e.g. "elevation", "focal_mean") as its band.
    nm <- names(r)
    if (any(blank(nm))) nm[blank(nm)] <- paste0("band_", which(blank(nm)))
    if (length(nm) == 1) nm <- row$layer_name
    names(r) <- nm

    ## Overviews are what GIS software draws when zoomed out. GDAL's default
    ## (cubic) invents in-between values, so class codes stop being classes and
    ## 0-1 indices go negative. Classes use nearest, continuous layers average.
    ## Full-resolution values are not affected.
    is_class <- !is_float ||
      grepl("action_map|potential_raster|bands_raster|landcover|pressure_count",
            row$layer_name)
    gdal_opts <- c("COMPRESS=DEFLATE",
                   paste0("PREDICTOR=", if (is_float) 3 else 2),
                   "OVERVIEWS=AUTO",
                   paste0("RESAMPLING=", if (is_class) "NEAREST" else "AVERAGE"),
                   "BIGTIFF=IF_SAFER", "NUM_THREADS=ALL_CPUS")
    args <- list(x = r, filename = out_file, filetype = "COG",
                 gdal = gdal_opts, overwrite = TRUE, NAflag = nd)
    if (!is.na(dt)) args$datatype <- dt
    do.call(terra::writeRaster, args)              # streams by block

    res_m <- terra::res(r)[1]
    cr <- terra::crs(r, describe = TRUE)
    epsg <- suppressWarnings(as.integer(cr$code[1]))
    cat(sprintf("  %s/%s.tif — %s m, %d band(s), %s, nodata %s, %s MB\n",
                theme_dir, row$layer_name, format(signif(res_m, 6)),
                terra::nlyr(r), dt, format(nd), format(size_mb(out_file), nsmall = 2)))

    reg_add(row$layer_name, file.path("rasters/native", theme_dir, paste0(row$layer_name, ".tif")),
            row$theme, "raster", signif(res_m, 6), ifelse(is.na(epsg), "", epsg),
            row$units, row$description, row$derived_from)

    ## Action maps are categorical with three classes; partners need the
    ## labels without opening a RAT-aware tool.
    if (grepl("action_map", row$layer_name) || grepl("action_map", basename(row$file))) {
      lk <- data.frame(value = 1:3,
                       label = c("Restore (active)", "Restore (support)", "Monitor"),
                       stringsAsFactors = FALSE)
      lk_file <- file.path(out_dir, paste0(row$layer_name, "_lookup.csv"))
      utils::write.csv(lk, lk_file, row.names = FALSE)
      reg_add(paste0(row$layer_name, "_lookup"),
              file.path("rasters/native", theme_dir, basename(lk_file)),
              row$theme, "table", "", "", "",
              sprintf("Class labels for %s.", row$layer_name), row$derived_from)
      cat(sprintf("    + %s\n", basename(lk_file)))
    }

    ## The refined land cover raster is coded; its crosswalk travels with it.
    ## Matched by exact stem: hlcm_source_runoff_risk.tif is a continuous risk
    ## surface, not a coded land cover raster, and must not pick up the lookup.
    if (grepl("hlcm_habmos_lca_reclass", basename(row$file), ignore.case = TRUE) &&
        file.exists(PATH_LC_LOOKUP)) {
      ## Named <layer>_lookup.csv like the other class lookups, and given a row
      ## for value 0, which the raster uses for cells outside the mapped area.
      lk <- utils::read.csv(PATH_LC_LOOKUP, stringsAsFactors = FALSE,
                            check.names = FALSE)
      if (!any(lk[[1]] == 0, na.rm = TRUE)) {
        lk0 <- lk[1, , drop = FALSE]
        lk0[1, ] <- ""
        lk0[[1]] <- 0
        lk0[[ncol(lk)]] <- "Outside the mapped area (no data)"
        lk <- rbind(lk0, lk)
      }
      lk_file <- file.path(out_dir, paste0(row$layer_name, "_lookup.csv"))
      utils::write.csv(lk, lk_file, row.names = FALSE)
      reg_add(paste0(row$layer_name, "_lookup"),
              file.path("rasters/native", theme_dir, basename(lk_file)),
              row$theme, "table", "", "", "",
              sprintf("Class codes and names for %s. Value 0 is outside the mapped area.",
                      row$layer_name),
              row$derived_from)
      cat(sprintf("    + %s\n", basename(lk_file)))
    }
    rm(r)
  }

  ##---- 10.5 cubes -----------------------------------------------------
  cat("\n-- Cubes --\n")
  cube_csvs <- fwd(Sys.glob(file.path(DESIGN_DIR, "SCAMP_*_cube_*_*.csv")))
  for (f in c(cube_tifs, cube_csvs)) {
    dest <- file.path(PKG_DIR, "rasters", "cubes", basename(f))
    if (grepl("\\.tif$", f)) {
      ## Re-written as a clean COG rather than copied byte for byte: values,
      ## band names and the existing (nearest-neighbour) overviews carry over,
      ## and anything else a tool has appended to the source file does not.
      if (file.exists(dest)) file.remove(dest)
      sf::gdal_utils("translate", f, dest, quiet = TRUE,
                     options = c("-of", "COG", "-co", "COMPRESS=DEFLATE",
                                 "-co", "PREDICTOR=3", "-co", "RESAMPLING=NEAREST",
                                 "-co", "BIGTIFF=IF_SAFER"))
    } else {
      file.copy(f, dest, overwrite = TRUE)
    }
    cat(sprintf("  %s (%s MB)\n", basename(f), format(size_mb(dest), nsmall = 2)))
    if (grepl("\\.tif$", f)) {
      cr <- terra::rast(dest)
      reg_add(base_noext(f), file.path("rasters/cubes", basename(f)),
              "design", "raster cube", signif(terra::res(cr)[1], 6),
              suppressWarnings(as.integer(terra::crs(cr, describe = TRUE)$code[1])),
              "see _bands.csv",
              sprintf("Aligned multi-band covariate cube, %d bands. Band names, units and gap-fill status in %s_bands.csv.",
                      terra::nlyr(cr), base_noext(f)),
              "see SOURCES_AND_LICENCES.md")
      rm(cr)
    } else {
      reg_add(base_noext(f), file.path("rasters/cubes", basename(f)),
              "design", "table", "", "", "",
              if (grepl("_bands\\.csv$", f)) "Band table for the matching cube: band number, name, units, resampling, gap-fill status and description."
              else "Code to class-name lookup for a categorical cube band.",
              "")
    }
  }

  ##---- 10.6 tables ----------------------------------------------------
  ## Band profiles and representativeness summaries, but only where the
  ## manifest asked for them.
  cat("\n-- Tables --\n")
  tab_rows <- inc[inc$type == "table", , drop = FALSE]
  if (nrow(tab_rows) == 0) {
    cat("  (no table rows in the manifest)\n")
  }
  for (i in seq_len(nrow(tab_rows))) {
    row <- tab_rows[i, ]
    dest <- file.path(PKG_DIR, "tables", paste0(row$layer_name, ".csv"))
    file.copy(row$file, dest, overwrite = TRUE)
    cat(sprintf("  %s\n", basename(dest)))
    reg_add(row$layer_name, file.path("tables", basename(dest)), row$theme, "table",
            row$resolution_m, row$crs_epsg, row$units, row$description, row$derived_from)
  }

  ##---- 10.7 data dictionary -------------------------------------------
  cat("\n-- Data dictionary --\n")
  DICT <- data.frame(layer = character(0), field = character(0), type = character(0),
                     description = character(0), units = character(0),
                     stringsAsFactors = FALSE)
  blank_dict <- character(0)

  ## layer_name -> source basename, so the inventory dictionary can be found.
  src_for_layer <- stats::setNames(basename(vec_rows$file), vec_rows$layer_name)

  for (L in written_layers) {
    hdr <- sf::st_read(GPKG, query = sprintf('SELECT * FROM "%s" LIMIT 0', L), quiet = TRUE)
    flds <- setdiff(names(hdr), attr(hdr, "sf_column"))
    if (L %in% names(src_for_layer)) {
      key <- src_for_layer[[L]]
    } else {
      ## Multi-layer source: layer_name had the source layer name appended, so
      ## the longest manifest layer_name that prefixes L is the right source.
      hit <- names(src_for_layer)[vapply(names(src_for_layer),
                                         function(n) startsWith(L, n), logical(1))]
      key <- if (length(hit)) src_for_layer[[hit[which.max(nchar(hit))]]] else ""
    }
    for (f in flds) {
      d <- lookup_field(key, f)
      DICT <- rbind(DICT, data.frame(
        layer = L, field = f, type = class(hdr[[f]])[1],
        description = d[["description"]], units = d[["units"]],
        stringsAsFactors = FALSE))
      if (!nzchar(d[["description"]])) blank_dict <- c(blank_dict, paste0(L, ".", f))
    }
  }
  ## Em-dashes out on the way to disk, not when the text is built, so there is
  ## one place to look for what the recipient actually reads.
  DICT[] <- lapply(DICT, de_emdash)
  utils::write.csv(DICT, file.path(PKG_DIR, "data_dictionary.csv"), row.names = FALSE, na = "")
  cat(sprintf("  %d field(s) across %d layer(s); %d without a description\n",
              nrow(DICT), length(written_layers), length(blank_dict)))

  ##---- 10.9 sources and licences --------------------------------------
  cat("\n-- SOURCES_AND_LICENCES.md --\n")
  src_terms <- unique(unlist(lapply(inc$derived_from, split_semi)))
  src_terms <- src_terms[nzchar(src_terms)]
  matched <- vapply(src_terms, match_source, character(1))
  unmatched <- src_terms[is.na(matched)]
  if (length(unmatched) > 0) {
    warn_note("derived_from entries not in the source catalogue (listed with the ",
              "default licence and no URL): ", paste(unmatched, collapse = "; "))
  }

  sl <- character(0)
  sl_add <- function(...) sl <<- c(sl, paste0(...))
  sl_add("# Sources and licences")
  sl_add("")
  sl_add(sprintf("SCAMP GIS Handover v%s. Generated %s.",
                 meta_get("version", default = "1"), format(Sys.Date(), "%d %B %Y")))
  sl_add("")
  sl_add(sprintf("Every layer in this package is derived. The package itself is licensed as %s.",
                 meta_get("licence", default = LICENCE_DEFAULT)))
  sl_add("The input datasets below keep their own licences and attribution")
  sl_add("requirements, which carry through to anything made from these layers.")
  sl_add("")
  sl_add("## 1. Input datasets feeding the shipped layers")
  sl_add("")
  sl_add("| Dataset | Provider | Licence | Attribution | URL |")
  sl_add("|---|---|---|---|---|")

  emit_source <- function(key, fallback_name) {
    if (is.na(key)) {
      sl_add(sprintf("| %s | not recorded | %s | not recorded | not recorded |",
                     fallback_name, LICENCE_DEFAULT))
      return(invisible())
    }
    s <- SOURCE_CATALOGUE[SOURCE_CATALOGUE$key == key, , drop = FALSE][1, ]
    lf <- licence_file_for(key)
    ## The recipient has no copy of our data/raw tree, so name the licence
    ## document rather than pointing at a path on this machine.
    lic <- if (!is.na(s$licence)) s$licence
           else if (!is.na(lf)) sprintf("See the provider's licence document (%s)",
                                        paste(unique(basename(strsplit(lf, "; ", fixed = TRUE)[[1]])),
                                              collapse = ", "))
           else LICENCE_DEFAULT
    sl_add(sprintf("| %s | %s | %s | %s | %s |", s$name, s$provider, lic,
                   ifelse(is.na(s$attribution), "not recorded", s$attribution),
                   ifelse(is.na(s$url), "not recorded", s$url)))
  }

  ## Occurrence databases come from the config, which carries their DOIs.
  occ <- meta_get("occurrence_sources", default = list())
  for (nm in names(occ)) {
    o <- occ[[nm]]
    sl_add(sprintf("| %s | %s | %s | %s | %s |",
                   o$name %||% nm, o$provider %||% "not recorded",
                   o$licence %||% LICENCE_DEFAULT, o$attribution %||% "not recorded",
                   o$url %||% "not recorded"))
  }
  occ_keys <- vapply(names(occ), match_source, character(1))
  for (k in setdiff(unique(matched[!is.na(matched)]), occ_keys)) {
    emit_source(k, k)
  }
  for (u in unmatched) emit_source(NA_character_, u)
  sl_add("")

  sl_add("## 2. Source data not included")
  sl_add("")
  sl_add("These third-party datasets are used by the analysis but are not")
  sl_add("redistributed here, either because they are unmodified clips of someone")
  sl_add("else's data or because they are raw records. Download them from the")
  sl_add("provider if you need them.")
  sl_add("")
  sl_add("| Layer | Why not included | Source | URL |")
  sl_add("|---|---|---|---|")

  ## This section is about other people's data, so it lists only the rows whose
  ## exclude_reason marks them as unmodified third-party clips or raw records.
  ## Superseded intermediates and our own diagnostics are also excluded from
  ## the package, but they are not someone else's data to point a partner at,
  ## and listing them here would read as though they were.
  exc <- manifest[!manifest$include_flag, , drop = FALSE]
  exc <- exc[grepl("^\\s*(Unmodified|Raw)\\b", exc$exclude_reason, ignore.case = TRUE), ,
             drop = FALSE]
  if (nrow(exc) == 0) {
    sl_add("| *none* | | | |")
    warn_note("SOURCES_AND_LICENCES.md section 2 is empty: no excluded row has an ",
              "exclude_reason starting \"Unmodified\" or \"Raw\". Check the manifest ",
              "wording if third-party layers were meant to be listed there.")
  }
  for (i in seq_len(nrow(exc))) {
    k <- NA_character_
    for (term in c(split_semi(exc$derived_from[i]), exc$layer_name[i], base_noext(exc$path[i]))) {
      k <- match_source(term)
      if (!is.na(k)) break
    }
    s <- if (is.na(k)) NULL else SOURCE_CATALOGUE[SOURCE_CATALOGUE$key == k, , drop = FALSE][1, ]
    sl_add(sprintf("| %s | %s | %s | %s |",
                   ifelse(nzchar(exc$layer_name[i]), exc$layer_name[i], basename(exc$path[i])),
                   ifelse(blank(exc$exclude_reason[i]), "not recorded", exc$exclude_reason[i]),
                   ifelse(is.null(s), "not recorded", s$provider),
                   ifelse(is.null(s) || is.na(s$url), "not recorded", s$url)))
  }
  sl_add("")
  ## No list of licence files: they are named in the section 1 table, and the
  ## paths they sit at here mean nothing to the recipient.
  writeLines(de_emdash(sl), file.path(PKG_DIR, "SOURCES_AND_LICENCES.md"))
  cat(sprintf("  %d input dataset(s); %d third-party layer(s) listed as not included\n",
              length(unique(c(names(occ), matched[!is.na(matched)], unmatched))), nrow(exc)))

  ##---- 10.9b figures --------------------------------------------------
  ## README.md and METHODS.md link to these; without them the recipient's
  ## copy has broken images. Listed non-recursively so the package folder
  ## being built underneath HANDOVER_DIR is never scanned into itself.
  cat("\n-- Figures --\n")
  fig_files <- fwd(list.files(HANDOVER_DIR, pattern = "\\.png$", full.names = TRUE,
                              recursive = FALSE, ignore.case = TRUE))
  if (length(fig_files) == 0) {
    cat("  (no figures found in outputs/handover/)\n")
    warn_note("No figures copied. README.md and METHODS.md reference ",
              "figures/environmental_bands_explainer.png, so those links will be ",
              "broken in the delivered package if a figure was expected.")
  }
  for (f in fig_files) {
    dest <- file.path(PKG_DIR, "figures", basename(f))
    file.copy(f, dest, overwrite = TRUE)
    cat(sprintf("  %s (%s MB)\n", basename(f), format(size_mb(f), nsmall = 2)))
    reg_add(base_noext(f), file.path("figures", basename(f)), "documentation",
            "figure", NA_character_, NA_character_, "n/a",
            "Explanatory figure referenced in README.md and METHODS.md",
            NA_character_)
  }

  ##---- 10.9c layer register -------------------------------------------
  ## Written last: every section that ships a file adds a row to REGISTER,
  ## the figures above included, so an earlier write would miss them.
  REGISTER[] <- lapply(REGISTER, de_emdash)
  utils::write.csv(REGISTER, file.path(PKG_DIR, "layer_register.csv"),
                   row.names = FALSE, na = "")
  cat(sprintf("\nlayer_register.csv: %d row(s)\n", nrow(REGISTER)))

  ##---- 10.10 documentation -------------------------------------------
  file.copy(PATH_README,  file.path(PKG_DIR, "README.md"),  overwrite = TRUE)
  file.copy(PATH_METHODS, file.path(PKG_DIR, "METHODS.md"), overwrite = TRUE)
  cat("\nREADME.md and METHODS.md copied verbatim.\n")

  ##---- 10.11 GEMINI 2.3 metadata -------------------------------------
  cat("\n-- GEMINI 2.3 metadata --\n")

  ## Derived, never hand-entered: bounding box, resolution range, CRS,
  ## formats, dates and the lineage summary.
  aoi_bbox <- NULL
  if (file.exists(PATH_AOI)) {
    aoi <- sf::st_read(PATH_AOI, quiet = TRUE)
    aoi_bbox <- sf::st_bbox(sf::st_transform(sf::st_geometry(aoi), 4326))
    cat(sprintf("  WGS84 bbox from %s: %.4f, %.4f, %.4f, %.4f\n", PATH_AOI,
                aoi_bbox["xmin"], aoi_bbox["ymin"], aoi_bbox["xmax"], aoi_bbox["ymax"]))
    rm(aoi)
  } else {
    warn_note("Study area not found at ", PATH_AOI,
              " — the GEMINI bounding box falls back to the shipped rasters.")
    ext_all <- NULL
    for (f in list.files(file.path(PKG_DIR, "rasters"), "\\.tif$", recursive = TRUE, full.names = TRUE)) {
      e <- sf::st_bbox(terra::rast(f))
      ext_all <- if (is.null(ext_all)) e else
        c(xmin = min(ext_all["xmin"], e["xmin"]), ymin = min(ext_all["ymin"], e["ymin"]),
          xmax = max(ext_all["xmax"], e["xmax"]), ymax = max(ext_all["ymax"], e["ymax"]))
    }
    if (!is.null(ext_all)) {
      g <- sf::st_as_sfc(sf::st_bbox(c(xmin = ext_all[["xmin"]], ymin = ext_all[["ymin"]],
                                       xmax = ext_all[["xmax"]], ymax = ext_all[["ymax"]]),
                                     crs = TARGET_EPSG))
      aoi_bbox <- sf::st_bbox(sf::st_transform(g, 4326))
    }
  }

  res_vals <- suppressWarnings(as.numeric(REGISTER$resolution_m[nzchar(REGISTER$resolution_m)]))
  res_vals <- res_vals[!is.na(res_vals) & res_vals > 0]
  res_min <- if (length(res_vals)) min(res_vals) else NA_real_
  res_max <- if (length(res_vals)) max(res_vals) else NA_real_

  ## Readable list of occurrence sources: "GBIF, OBIS and the NBN Atlas".
  occ_nice <- c(GBIF = "GBIF", OBIS = "OBIS", NBN = "the NBN Atlas")[names(occ)]
  occ_nice[is.na(occ_nice)] <- names(occ)[is.na(occ_nice)]
  occ_txt <- if (length(occ_nice) > 1)
    paste(paste(occ_nice[-length(occ_nice)], collapse = ", "), "and", occ_nice[length(occ_nice)])
  else paste(occ_nice, collapse = "")

  lineage <- paste0(
    "All layers are derived products. Open third-party datasets (listed in ",
    "SOURCES_AND_LICENCES.md) were clipped to the study area, harmonised to ",
    "EPSG:27700 and combined through habitat suitability modelling, pressure ",
    "normalisation and accumulation, SCIMAP-style diffuse pollution routing, ",
    "occurrence-record quality control and interpolation, larval connectivity ",
    "from published dispersal probabilities (coastalNet), and k-means ",
    "environmental banding. Raster layers are supplied at the resolution at ",
    "which they were produced, except the physical depth, temperature and ",
    "current rasters, which are taken from the marine cube. The two cubes are ",
    "resampled and, for some marine covariates, gap-filled onto a common grid. ",
    "Occurrence records derive from ",
    occ_txt,
    ". Full method descriptions, parameters and references are in METHODS.md, ",
    "distributed with this package. No analysis code is distributed."
  )

  ## geometa's class and method names have shifted between releases. Rather
  ## than pin this script to one version, try the known spellings and report
  ## what did not take, so a partial record is still written and the gaps are
  ## visible instead of silent.
  m_call <- function(obj, candidates, ...) {
    for (m in candidates) {
      if (!is.null(obj[[m]]) && is.function(obj[[m]])) {
        res <- try(obj[[m]](...), silent = TRUE)
        if (!inherits(res, "try-error")) return(invisible(TRUE))
      }
    }
    warn_note("GEMINI: none of these geometa methods worked: ",
              paste(candidates, collapse = ", "))
    invisible(FALSE)
  }

  geometa_ns <- asNamespace("geometa")
  new_iso <- function(candidates, ...) {
    for (cl in candidates) {
      if (!exists(cl, envir = geometa_ns, inherits = FALSE)) next
      gen <- get(cl, envir = geometa_ns)
      o <- try(gen$new(...), silent = TRUE)
      if (!inherits(o, "try-error")) return(o)
    }
    warn_note("GEMINI: no usable geometa class among ", paste(candidates, collapse = ", "))
    NULL
  }

  ## Element-level guard: a GEMINI element the installed geometa cannot set
  ## becomes a named warning, not a lost record.
  safe <- function(label, expr) {
    res <- try(expr, silent = TRUE)
    if (inherits(res, "try-error")) {
      warn_note("GEMINI element not set (", label, "): ",
                trimws(conditionMessage(attr(res, "condition"))))
      return(invisible(FALSE))
    }
    invisible(TRUE)
  }

  rp_from <- function(cfg, default_role = "pointOfContact") {
    rp <- geometa::ISOResponsibleParty$new()
    if (!is.null(cfg$individual) && !is.na(cfg$individual)) rp$setIndividualName(cfg$individual)
    rp$setOrganisationName(cfg$organisation %||% "not recorded")
    if (!is.null(cfg$position) && !is.na(cfg$position)) rp$setPositionName(cfg$position)
    rp$setRole(cfg$role %||% default_role)
    ci <- geometa::ISOContact$new()
    addr <- geometa::ISOAddress$new()
    addr$setEmail(cfg$email %||% "not recorded")
    ci$setAddress(addr)
    rp$setContactInfo(ci)
    rp
  }

  gemini_path <- file.path(PKG_DIR, "metadata", paste0(PKG_NAME, "_gemini.xml"))
  validation_msg <- "not run"

  gemini_ok <- try({
    md <- geometa::ISOMetadata$new()

    ## --- metadata about the metadata ---
    safe("fileIdentifier", md$setFileIdentifier(meta_get("identifier", "code", default = PKG_NAME)))
    safe("characterSet",   md$setCharacterSet("utf8"))
    safe("metadata language", md$setLanguage(meta_get("language", default = "eng")))
    safe("metadata date",  md$setDateStamp(Sys.Date()))
    safe("standard name",  md$setMetadataStandardName("ISO 19115 (UK GEMINI)"))
    safe("standard version", md$setMetadataStandardVersion("2.3"))
    m_call(md, c("addHierarchyLevel", "setHierarchyLevel"), "series")   # resource type
    safe("metadata point of contact",
         md$addContact(rp_from(meta_get("metadata_contact", default = list()), "pointOfContact")))

    ## --- spatial reference system ---
    safe("reference system", {
      rs <- geometa::ISOReferenceSystem$new()
      rs$setReferenceSystemIdentifier(
        geometa::ISOReferenceIdentifier$new(code = as.character(TARGET_EPSG), codeSpace = "EPSG"))
      md$addReferenceSystemInfo(rs)
    })

    ## --- identification ---
    ident <- geometa::ISODataIdentification$new()

    ct <- geometa::ISOCitation$new()
    safe("title", ct$setTitle(meta_get("title", default = PKG_NAME)))
    m_call(ct, c("addAlternateTitle", "setAlternateTitle"),
           meta_get("alternative_title", default = PKG_NAME))

    ## Dataset reference date: the publication date if the analyst has set one,
    ## otherwise the build date as a creation date.
    pub <- meta_get("publication_date", default = NA_character_)
    safe("dataset reference date", {
      d1 <- geometa::ISODate$new()
      d1$setDate(if (!is.na(pub)) as.Date(pub) else Sys.Date())
      d1$setDateType(if (!is.na(pub)) "publication" else "creation")
      ct$addDate(d1)
      if (!is.na(pub)) {
        d2 <- geometa::ISODate$new(); d2$setDate(Sys.Date()); d2$setDateType("creation")
        ct$addDate(d2)
      }
    })

    ## Unique resource identifier.
    safe("unique resource identifier", {
      cs <- meta_get("identifier", "codespace", default = NA_character_)
      uid <- if (is.na(cs)) {
        geometa::ISOMetaIdentifier$new(code = meta_get("identifier", "code", default = PKG_NAME))
      } else {
        geometa::ISOMetaIdentifier$new(code = meta_get("identifier", "code", default = PKG_NAME),
                                       codeSpace = cs)
      }
      ct$addIdentifier(uid)
    })
    safe("citation", ident$setCitation(ct))

    safe("abstract", ident$setAbstract(meta_get("abstract", default = "not recorded")))
    safe("topic category", {
      for (tc in meta_get("topic_categories", default = "environment")) ident$addTopicCategory(tc)
    })
    safe("dataset language", ident$addLanguage(meta_get("language", default = "eng")))
    safe("dataset character set", ident$addCharacterSet("utf8"))

    ## The package holds both vector layers and rasters.
    for (srt in c("vector", "grid")) m_call(ident, "addSpatialRepresentationType", srt)

    safe("keywords", {
      kw_free <- geometa::ISOKeywords$new()
      for (k in meta_get("keywords", default = character(0))) kw_free$addKeyword(k)
      kw_free$setKeywordType("theme")
      ident$addKeywords(kw_free)
    })

    ## GEMINI 2.3 requires at least one INSPIRE theme, cited to the GEMET
    ## thesaurus rather than supplied as a free-text keyword.
    safe("INSPIRE theme keywords", {
      kw_insp <- geometa::ISOKeywords$new()
      for (k in meta_get("inspire_themes", default = "Habitats and biotopes")) kw_insp$addKeyword(k)
      kw_insp$setKeywordType("theme")
      th <- geometa::ISOCitation$new()
      th$setTitle("GEMET - INSPIRE themes, version 1.0")
      thd <- geometa::ISODate$new()
      thd$setDate(as.Date("2008-06-01")); thd$setDateType("publication")
      th$addDate(thd)
      kw_insp$setThesaurusName(th)
      ident$addKeywords(kw_insp)
    })

    safe("frequency of update", {
      mi <- geometa::ISOMaintenanceInformation$new()
      mi$setMaintenanceFrequency(meta_get("frequency", default = "notPlanned"))
      ident$addResourceMaintenance(mi)
    })

    safe("constraints", {
      lc <- geometa::ISOLegalConstraints$new()
      lc$addUseLimitation(meta_get("use_constraints", default = "not recorded"))
      lc$addAccessConstraint("otherRestrictions")
      lc$addUseConstraint("otherRestrictions")
      m_call(lc, "addOtherConstraint", meta_get("access_limitations", default = "not recorded"))
      m_call(lc, "addOtherConstraint",
             trimws(sprintf("%s (%s). %s", meta_get("licence", default = LICENCE_DEFAULT),
                            meta_get("licence_url", default = ""),
                            meta_get("attribution", default = ""))))
      ident$addResourceConstraints(lc)
    })

    safe("responsible organisation",
         ident$addPointOfContact(rp_from(meta_get("responsible_org", default = list()), "owner")))

    ## Spatial resolution as a range: the finest and coarsest shipped raster.
    if (!is.na(res_min)) {
      safe("spatial resolution", {
        for (v in unique(c(res_min, res_max))) {
          rr <- geometa::ISOResolution$new()
          rr$setDistance(geometa::ISODistance$new(value = v, uom = "m", useUomURI = FALSE))
          ident$addSpatialResolution(rr)
        }
      })
    } else {
      warn_note("GEMINI: no raster resolutions available, spatial resolution omitted.")
    }

    safe("extent", {
      exg <- geometa::ISOExtent$new()
      if (!is.null(aoi_bbox)) {
        bb <- geometa::ISOGeographicBoundingBox$new(
          minx = as.numeric(aoi_bbox["xmin"]), miny = as.numeric(aoi_bbox["ymin"]),
          maxx = as.numeric(aoi_bbox["xmax"]), maxy = as.numeric(aoi_bbox["ymax"]))
        exg$addGeographicElement(bb)
      }
      ident$addExtent(exg)
    })
    md$addIdentificationInfo(ident)

    ## --- distribution ---
    safe("distribution format", {
      dist <- geometa::ISODistribution$new()
      for (fmt in list(c("GeoPackage", "1.3"), c("Cloud Optimised GeoTIFF", "1.0"),
                       c("CSV", "RFC 4180"))) {
        ft <- geometa::ISOFormat$new()
        ft$setName(fmt[1]); ft$setVersion(fmt[2])
        dist$addFormat(ft)
      }
      rl <- meta_get("resource_locator", default = NA_character_)
      if (!is.na(rl)) {
        dto <- geometa::ISODigitalTransferOptions$new()
        or <- geometa::ISOOnlineResource$new()
        or$setLinkage(rl)
        or$setName(PKG_NAME)
        or$setProtocol("WWW:LINK-1.0-http--link")
        dto$addOnlineResource(or)
        m_call(dist, c("setDigitalTransferOptions", "addDigitalTransferOptions"), dto)
      }
      md$setDistributionInfo(dist)
    })

    ## --- lineage ---
    safe("lineage", {
      dq <- geometa::ISODataQuality$new()
      sc <- new_iso(c("ISOScope", "ISODataQualityScope"))
      if (!is.null(sc)) {
        sc$setLevel("series")
        m_call(dq, "setScope", sc)
      }
      lin <- geometa::ISOLineage$new()
      lin$setStatement(lineage)
      dq$setLineage(lin)
      md$addDataQualityInfo(dq)
    })

    md$save(gemini_path)
    TRUE
  }, silent = TRUE)

  if (inherits(gemini_ok, "try-error")) {
    warn_note("GEMINI record could not be built: ",
              trimws(conditionMessage(attr(gemini_ok, "condition"))),
              " — check the geometa version and handover_metadata.R.")
  } else {
    cat(sprintf("  %s written\n", gemini_path))
    v <- try(utils::capture.output(md$validate(strict = FALSE)), silent = TRUE)
    if (inherits(v, "try-error")) {
      validation_msg <- paste("ISO 19139 validation could not run:",
                              trimws(conditionMessage(attr(v, "condition"))))
    } else {
      validation_msg <- if (length(v) == 0) "ISO 19139 validation passed with no messages." else
        paste(v, collapse = " | ")
    }
    cat(sprintf("  ISO validation: %s\n", validation_msg))
  }

  ##---- 10.12 checksums -----------------------------------------------
  cat("\n-- Checksums --\n")
  pkg_files <- list.files(PKG_DIR, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  pkg_files <- pkg_files[!grepl("checksums\\.txt$", pkg_files)]
  md5 <- tools::md5sum(pkg_files)
  chk <- sprintf("%s  %s", md5, fwd(sub(paste0("^", PKG_DIR, "/?"), "", fwd(names(md5)))))
  writeLines(chk, file.path(PKG_DIR, "checksums.txt"))
  cat(sprintf("  %d file(s) hashed\n", length(chk)))

  ##---- 10.13 final report --------------------------------------------
  pkg_files <- list.files(PKG_DIR, recursive = TRUE, full.names = TRUE)
  sizes <- file.size(pkg_files)
  rel <- fwd(sub(paste0("^", PKG_DIR, "/?"), "", fwd(pkg_files)))
  ord <- order(-sizes)

  cat("\n==================== HANDOVER PACKAGE ====================\n")
  cat(sprintf("%s\n", PKG_DIR))
  cat(sprintf("%d files, %s MB total\n\n", length(pkg_files),
              format(round(sum(sizes) / 1024^2, 1), nsmall = 1)))
  for (i in ord) {
    cat(sprintf("  %9.2f MB  %s\n", sizes[i] / 1024^2, rel[i]))
  }
  cat(sprintf("\nVector layers: %d in vectors/SCAMP_vectors.gpkg\n", length(written_layers)))
  cat(sprintf("Native rasters: %d\n", nrow(ras_rows)))
  cat(sprintf("Cubes: %d (+%d sidecar CSVs)\n", length(cube_tifs), length(cube_csvs)))
  cat(sprintf("Register rows: %d\n", nrow(REGISTER)))

  cat(sprintf("\nData dictionary: %d field(s), %d without a description\n",
              nrow(DICT), length(blank_dict)))
  if (length(blank_dict) > 0) {
    cat("  Undescribed fields (add them to PROJECT_INVENTORY.md section 7):\n")
    for (b in blank_dict) cat(sprintf("    - %s\n", b))
  }

  cat(sprintf("\nGEMINI: %s\n", basename(gemini_path)))
  cat(sprintf("  ISO 19139 validation: %s\n", validation_msg))
  cat("  GEMINI schematron validation has NOT been run. ISO 19139 validation\n")
  cat("  checks the XML against the ISO schema only; it does not check the UK\n")
  cat("  GEMINI 2.3 profile rules. Run the record through the schematron\n")
  cat("  separately before publishing it to a UK discovery service.\n")

  if (length(WARNINGS) == 0) {
    cat("\nNo warnings.\n")
  } else {
    cat(sprintf("\n%d WARNING(S):\n", length(WARNINGS)))
    for (w in WARNINGS) cat(sprintf("  - %s\n", w))
  }

  ##---- 10.14 optional zip --------------------------------------------
  if (ZIP) {
    cat("\n-- Zipping --\n")
    zip_path <- normalizePath(file.path(HANDOVER_DIR, paste0(PKG_NAME, ".zip")),
                              winslash = "/", mustWork = FALSE)
    ## Zip from the handover directory so the archive holds one top-level
    ## folder rather than the whole project path.
    wd <- getwd()
    setwd(HANDOVER_DIR)
    z <- try(utils::zip(zipfile = basename(zip_path),
                        files = file.path(PKG_NAME, list.files(PKG_NAME, recursive = TRUE)),
                        flags = "-r9Xq"), silent = TRUE)
    setwd(wd)
    if (inherits(z, "try-error") || !file.exists(zip_path)) {
      warn_note("Zip failed (is a zip executable on PATH?). The package folder is complete.")
    } else {
      cat(sprintf("  %s (%s MB)\n", zip_path, format(size_mb(zip_path), nsmall = 1)))
    }
  }

  cat("\nStage 2 complete.\n")
}
