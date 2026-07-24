##=====================================================================
## SCAMP GIS Analysis — Theme D: Species
## 02_species.R
##=====================================================================
## Produces four outputs:
##   1. Species richness per grid cell (all taxa, effort-adjusted)
##   2. INNS occurrence records
##   3. Indicator species occurrence records (per habitat)
##   4. Priority Marine Features (PMF) occurrence records
##
## Data sources: GBIF, OBIS, NBN Atlas
## Grid: 5 km (appropriate for mixed-precision records)
##
## Outputs (all in outputs/species/):
##   - richness_grid.gpkg        (richness + effort per cell)
##   - inns_records.gpkg         (INNS occurrence points)
##   - indicator_records.gpkg    (indicator species points)
##   - pmf_records.gpkg          (PMF species points)
##   - species_summary.csv       (presence/absence per species)
##=====================================================================

require(rgbif)
require(robis)
require(galah)
require(sf)
require(terra)
require(dplyr)
require(tidyr)
require(readr)
require(stringr)

##=====================================================================
## 1. SETUP
##=====================================================================

## --- Defining seascape ---
smr <- st_read(
  "data/raw/boundaries/administrative_units_scottish_marine_regions/administrative_units_scottish_marine_regionsPolygon.shp",
  quiet = TRUE
)

las <- st_read(
  "data/raw/boundaries/Local_Authority_Boundaries_-_Scotland/pub_las.shp",
  quiet = TRUE
)

cat("=== Scottish Marine Regions ===\n")
str(smr)

cat("\n=== Local Authority Boundaries ===\n")
str(las)

smr_solway <- smr %>% filter(grepl("Solway", objnam, ignore.case = TRUE))
las_dg     <- las %>% filter(grepl("dumfries", local_auth, ignore.case = TRUE))

plot(st_geometry(smr_solway))
plot(las_dg)

scamp_aoi <- st_union(
  c(st_geometry(st_transform(smr_solway, 27700)), st_geometry(las_dg))
)

plot(st_geometry(scamp_aoi))

Sys.setenv(SHAPE_RESTORE_SHX = "YES")

biosphere <- st_read("data/raw/boundaries/Biosphere_SCAMP_Zone.shp", quiet = TRUE)
mcnabb    <- st_read("data/raw/boundaries/DG_McNabb_5milezone.shp",  quiet = TRUE)

cat("=== Biosphere SCAMP Zone ===\n")
str(biosphere)

cat("\n=== McNabb 5 Mile Zone ===\n")
str(mcnabb)

plot(st_geometry(biosphere))
plot(st_geometry(mcnabb))

biosphere_res <- st_read(
  "data/raw/boundaries/Biosphere_Reserves_-2410630408604070948.gpkg",
  quiet = TRUE
) %>%
  st_transform(27700)

gsab_scamp <- biosphere_res %>%
  filter(lengths(st_intersects(., scamp_aoi)) > 0) %>%
  st_intersection(scamp_aoi)

cat("=== GSAB intersecting SCAMP AOI ===\n")
str(gsab_scamp)
plot(st_geometry(gsab_scamp))

## Try UTM 30N first — the northings strongly suggest this
mcnabb_test <- st_set_crs(mcnabb, 27700)

## Transform to WGS84 and check the coordinates make sense
## Should be roughly -5 to -3 longitude, 54.5 to 55.2 latitude
st_bbox(st_transform(mcnabb_test, 4326))




aoi_new <- st_intersection(
  scamp_aoi,
  st_union(st_geometry(st_transform(mcnabb, 27700)))
)

plot(st_geometry(aoi_new))

gsab_scamp_new <- biosphere_res %>%
  filter(lengths(st_intersects(., aoi_new)) > 0) %>%
  st_intersection(aoi_new)

cat("=== GSAB clipped to McNabb AOI ===\n")
str(gsab_scamp_new)
plot(st_geometry(gsab_scamp_new))



galah_config(atlas = "United Kingdom", email = Sys.getenv("GALAH_EMAIL"))

## Create output directory
dir.create("outputs/species", showWarnings = FALSE, recursive = TRUE)

## SCAMP bounding box — Solway Firth region
scamp_bbox <- c(xmin = -5.1, xmax = -3.0, ymin = 54.6, ymax = 55.2)

scamp_wkt <- sprintf(
  "POLYGON((%f %f,%f %f,%f %f,%f %f,%f %f))",
  scamp_bbox["xmin"], scamp_bbox["ymin"],
  scamp_bbox["xmax"], scamp_bbox["ymin"],
  scamp_bbox["xmax"], scamp_bbox["ymax"],
  scamp_bbox["xmin"], scamp_bbox["ymax"],
  scamp_bbox["xmin"], scamp_bbox["ymin"]
)

## 5 km grid in OSGB36 for richness calculation
scamp_extent <- ext(-5.1, -3.0, 54.6, 55.2)
grid_template <- rast(scamp_extent, res = 5000 / 111320, crs = "EPSG:4326")

##=====================================================================
## 2. SPECIES LISTS
##=====================================================================

## --- INNS target list (Annex A.6) ---
inns_list <- tribble(
  ~latin_name,                ~common_name,
  "Spartina anglica",         "Common cord-grass",
  "Crassostrea gigas",        "Pacific oyster",
  "Crepidula fornicata",      "Slipper limpet",
  "Dikerogammarus villosus",  "Killer shrimp",
  "Dreissena polymorpha",     "Zebra mussel",
  "Corella eumyota",          "Orange-tipped sea squirt",
  "Didemnum vexillum",        "Carpet sea squirt",
  "Styela clava",             "Stalked sea squirt",
  "Sargassum muticum",        "Japanese wireweed",
  "Ocinebrellus inornatus",   "Japanese oyster drill",
  "Urosalpinx cinerea",       "American oyster drill"
) %>% mutate(category = "INNS")

## --- Indicator species (Annexes A.1–A.5) ---
indicator_list <- tribble(
  ~latin_name,                ~common_name,              ~habitat,
  # Seagrass
  "Zostera marina",           "Common eelgrass",         "Seagrass",
  "Zostera noltii",           "Dwarf eelgrass",          "Seagrass",
  "Dicentrarchus labrax",     "European bass",           "Seagrass",
  "Clupea harengus",          "Herring",                 "Seagrass",
  # Oyster
  "Ostrea edulis",            "Native oyster",           "Oyster",
  "Mytilus edulis",           "Blue mussel",             "Oyster",
  "Lanice conchilega",        "Sand mason worm",         "Oyster",
  "Arenicola marina",         "Lugworm",                 "Oyster",
  "Cerastoderma edule",       "Common cockle",           "Oyster",
  "Hediste diversicolor",     "Ragworm",                 "Oyster",
  "Corophium volutator",      "Mud shrimp",              "Oyster",
  "Scrobicularia plana",      "Peppery furrow shell",    "Oyster",
  "Peringia ulvae",           "Mud snail",               "Oyster",
  "Somateria mollissima",     "Eider duck",              "Oyster",
  "Haematopus ostralegus",    "Oystercatcher",           "Oyster",
  "Asterias rubens",          "Common starfish",         "Oyster",
  "Carcinus maenas",          "Shore crab",              "Oyster",
  "Ocenebra erinacea",        "European sting winkle",   "Oyster",
  # Saltmarsh
  "Salicornia europaea",      "Glasswort",               "Saltmarsh",
  "Puccinellia maritima",     "Saltmarsh grass",         "Saltmarsh",
  "Atriplex portulacoides",   "Sea purslane",            "Saltmarsh",
  "Armeria maritima",         "Thrift",                  "Saltmarsh",
  "Limonium vulgare",         "Sea lavender",            "Saltmarsh",
  "Plantago maritima",        "Sea plantain",            "Saltmarsh",
  "Phragmites australis",     "Common reed",             "Saltmarsh",
  "Tringa totanus",           "Redshank",                "Saltmarsh",
  "Aristotelia brizella",     "Saltmarsh micromoth",     "Saltmarsh",
  "Gortyna borelii lunata",   "Fisher's estuarine moth", "Saltmarsh",
  # Sand dunes
  "Epidalea calamita",        "Natterjack toad",         "Sand dunes",
  "Lacerta agilis",           "Sand lizard",             "Sand dunes",
  "Cicindela hybrida",        "Northern dune tiger beetle","Sand dunes",
  "Petalophyllum ralfsii",    "Petalwort",               "Sand dunes",
  "Sternula albifrons",       "Little tern",             "Sand dunes",
  "Arvicola amphibius",       "Water vole",              "Sand dunes",
  "Meles meles",              "Badger",                  "Sand dunes",
  "Sciurus vulgaris",         "Red squirrel",            "Sand dunes",
  # Coastal woodland
  "Lutra lutra",              "Otter",                   "Woodland"
) %>% mutate(category = "Indicator")

## --- PMF species relevant to the Solway ---
## Key subset; full PMF ranges come from the GeMS download (Theme C)
pmf_list <- tribble(
  ~latin_name,                  ~common_name,
  "Ostrea edulis",              "Native oyster",
  "Zostera marina",             "Common eelgrass",
  "Arctica islandica",          "Ocean quahog",
  "Modiolus modiolus",          "Horse mussel",
  "Nucella lapillus",           "Dog whelk",
  "Cetorhinus maximus",         "Basking shark",
  "Squalus acanthias",          "Spurdog",
  "Raja clavata",               "Thornback ray",
  "Raja montagui",              "Spotted ray",
  "Leucoraja naevus",           "Cuckoo ray",
  "Gadus morhua",               "Atlantic cod",
  "Hippocampus guttulatus",     "Long-snouted seahorse",
  "Phocoena phocoena",          "Harbour porpoise",
  "Tursiops truncatus",         "Bottlenose dolphin",
  "Halichoerus grypus",         "Grey seal",
  "Phoca vitulina",             "Harbour seal",
  "Anguilla anguilla",          "European eel",
  "Salmo salar",                "Atlantic salmon",
  "Petromyzon marinus",         "Sea lamprey",
  "Lampetra fluviatilis",       "River lamprey",
  "Clupea harengus",            "Herring",
  "Dicentrarchus labrax",       "European bass"
) %>% mutate(category = "PMF")

cat(sprintf("Species lists compiled: %d INNS, %d indicators, %d PMF\n",
    nrow(inns_list), nrow(indicator_list), nrow(pmf_list)))

##=====================================================================
## 3. QUERY FUNCTIONS
##=====================================================================

## Standardised query function for GBIF
query_gbif <- function(species_names = NULL, wkt, limit = 500) {
  cat("Querying GBIF...\n")

  if (is.null(species_names)) {
    ## Bulk download — all taxa in bounding box
    cat("  Bulk query (all taxa)...\n")
    recs <- rgbif::occ_data(
      geometry = wkt,
      hasCoordinate = TRUE,
      limit = 100000
    )
    if (!is.null(recs$data) && nrow(recs$data) > 0) {
      out <- recs$data %>%
        transmute(
          latin_name = species,
          longitude = decimalLongitude,
          latitude = decimalLatitude,
          year = year,
          source = "GBIF"
        ) %>%
        filter(!is.na(latin_name))
      cat(sprintf("  %d records retrieved\n", nrow(out)))
      return(out)
    }
    return(tibble())
  }

  ## Per-species query
  results <- list()
  for (i in seq_along(species_names)) {
    sp <- species_names[i]
    cat(sprintf("  [%d/%d] %s", i, length(species_names), sp))
    tryCatch({
      recs <- rgbif::occ_data(
        scientificName = sp,
        geometry = wkt,
        hasCoordinate = TRUE,
        limit = limit
      )
      if (!is.null(recs$data) && nrow(recs$data) > 0) {
        results[[sp]] <- recs$data %>%
          transmute(
            latin_name = species,
            longitude = decimalLongitude,
            latitude = decimalLatitude,
            year = year,
            source = "GBIF"
          )
        cat(sprintf(" — %d records\n", nrow(results[[sp]])))
      } else {
        cat(" — 0 records\n")
      }
    }, error = function(e) cat(sprintf(" — error: %s\n", e$message)))
  }
  bind_rows(results)
}

## Standardised query function for OBIS
query_obis <- function(species_names = NULL, wkt, limit = 500) {
  cat("Querying OBIS...\n")

  if (is.null(species_names)) {
    ## Bulk download — all taxa in bounding box
    cat("  Bulk query (all taxa)...\n")
    recs <- robis::occurrence(geometry = wkt)
    if (!is.null(recs) && nrow(recs) > 0) {
      out <- recs %>%
        transmute(
          latin_name = species,
          longitude = decimalLongitude,
          latitude = decimalLatitude,
          year = date_year,
          source = "OBIS"
        ) %>%
        filter(!is.na(latin_name))
      cat(sprintf("  %d records retrieved\n", nrow(out)))
      return(out)
    }
    return(tibble())
  }

  ## Per-species query
  results <- list()
  for (i in seq_along(species_names)) {
    sp <- species_names[i]
    cat(sprintf("  [%d/%d] %s", i, length(species_names), sp))
    tryCatch({
      recs <- robis::occurrence(scientificname = sp, geometry = wkt)
      if (!is.null(recs) && nrow(recs) > 0) {
        recs <- head(recs, limit)
        results[[sp]] <- recs %>%
          transmute(
            latin_name = species,
            longitude = decimalLongitude,
            latitude = decimalLatitude,
            year = date_year,
            source = "OBIS"
          )
        cat(sprintf(" — %d records\n", nrow(results[[sp]])))
      } else {
        cat(" — 0 records\n")
      }
    }, error = function(e) cat(sprintf(" — error: %s\n", e$message)))
  }
  bind_rows(results)
}

## Standardised query function for NBN Atlas
query_nbn <- function(species_names = NULL, bbox) {
  cat("Querying NBN Atlas...\n")

  if (is.null(species_names)) {
    ## NBN bulk download not practical without species filter
    ## Use GBIF + OBIS for richness; NBN for targeted queries only
    cat("  Skipping NBN for bulk query (use GBIF + OBIS)\n")
    return(tibble())
  }

  results <- list()
  for (i in seq_along(species_names)) {
    sp <- species_names[i]
    cat(sprintf("  [%d/%d] %s", i, length(species_names), sp))
    tryCatch({
      recs <- galah_call() %>%
        identify(sp) %>%
        filter(
          decimalLongitude >= bbox["xmin"],
          decimalLongitude <= bbox["xmax"],
          decimalLatitude  >= bbox["ymin"],
          decimalLatitude  <= bbox["ymax"]
        ) %>%
        atlas_occurrences()
      if (!is.null(recs) && nrow(recs) > 0) {
        results[[sp]] <- recs %>%
          transmute(
            latin_name = scientificName,
            longitude = decimalLongitude,
            latitude = decimalLatitude,
            year = as.integer(substr(eventDate, 1, 4)),
            source = "NBN"
          )
        cat(sprintf(" — %d records\n", nrow(results[[sp]])))
      } else {
        cat(" — 0 records\n")
      }
    }, error = function(e) cat(sprintf(" — error: %s\n", e$message)))
  }
  bind_rows(results)
}

##=====================================================================
## 4. DEDUPLICATION
##=====================================================================
## NBN feeds into GBIF so records can appear in multiple databases.
## Deduplicate on species + coordinates (rounded to 3dp ~100m) + year

deduplicate_records <- function(df) {
  n_before <- nrow(df)
  df <- df %>%
    mutate(
      lon_r = round(longitude, 3),
      lat_r = round(latitude, 3)
    ) %>%
    distinct(latin_name, lon_r, lat_r, year, .keep_all = TRUE) %>%
    select(-lon_r, -lat_r)
  cat(sprintf("Deduplication: %d → %d records (%d duplicates removed)\n",
      n_before, nrow(df), n_before - nrow(df)))
  df
}

##=====================================================================
## 5. MAP 1 — SPECIES RICHNESS PER GRID CELL
##=====================================================================
## Bulk download all taxa from GBIF + OBIS for the Solway bbox,
## then count unique species per 5 km grid cell with effort layer

cat("\n========== MAP 1: SPECIES RICHNESS ==========\n\n")

richness_gbif <- query_gbif(species_names = NULL, wkt = scamp_wkt)
richness_obis <- query_obis(species_names = NULL, wkt = scamp_wkt)

richness_all <- bind_rows(richness_gbif, richness_obis) %>%
  filter(!is.na(latitude), !is.na(longitude), !is.na(latin_name)) %>%
  deduplicate_records()

## Convert to sf
richness_sf <- richness_all %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

## Create 5 km grid in OSGB36 then back to WGS84
grid_osgb <- st_make_grid(
  st_transform(st_as_sfc(st_bbox(richness_sf)), 27700),
  cellsize = 5000, square = TRUE
) %>%
  st_sf(cell_id = seq_along(.), geometry = .) %>%
  st_transform(4326)

## Spatial join: assign each record to a grid cell
records_gridded <- st_join(richness_sf, grid_osgb, join = st_within)

## Calculate richness and effort per cell
richness_grid <- records_gridded %>%
  st_drop_geometry() %>%
  group_by(cell_id) %>%
  summarise(
    species_richness = n_distinct(latin_name),
    record_count = n(),
    .groups = "drop"
  ) %>%
  filter(!is.na(cell_id))

## Join back to grid polygons
richness_grid_sf <- grid_osgb %>%
  left_join(richness_grid, by = "cell_id") %>%
  mutate(
    species_richness = replace_na(species_richness, 0L),
    record_count = replace_na(record_count, 0L)
  ) %>%
  filter(record_count > 0)  # drop empty cells

st_write(richness_grid_sf, "outputs/species/richness_grid.gpkg",
         delete_dsn = TRUE)
cat(sprintf("Saved richness grid: %d cells, max richness = %d\n",
    nrow(richness_grid_sf), max(richness_grid_sf$species_richness)))

rm(richness_gbif, richness_obis, richness_all, richness_sf,
   records_gridded, richness_grid)
gc()

##=====================================================================
## 6. MAP 2 — INNS OCCURRENCE RECORDS
##=====================================================================

cat("\n========== MAP 2: INNS OCCURRENCES ==========\n\n")

inns_names <- unique(inns_list$latin_name)

inns_gbif <- query_gbif(inns_names, wkt = scamp_wkt)
inns_obis <- query_obis(inns_names, wkt = scamp_wkt)
inns_nbn  <- query_nbn(inns_names, bbox = scamp_bbox)

inns_records <- bind_rows(inns_gbif, inns_obis, inns_nbn) %>%
  filter(!is.na(latitude), !is.na(longitude)) %>%
  deduplicate_records() %>%
  left_join(inns_list %>% select(latin_name, common_name), by = "latin_name")

if (nrow(inns_records) > 0) {
  inns_sf <- inns_records %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
  st_write(inns_sf, "outputs/species/inns_records.gpkg", delete_dsn = TRUE)
}

cat(sprintf("INNS: %d records across %d species\n",
    nrow(inns_records), n_distinct(inns_records$latin_name)))

## INNS summary
inns_summary <- inns_list %>%
  left_join(
    inns_records %>% count(latin_name, name = "n_records"),
    by = "latin_name"
  ) %>%
  mutate(n_records = replace_na(n_records, 0L)) %>%
  arrange(desc(n_records))

cat("\nINNS presence in Solway:\n")
print(inns_summary, n = Inf)

rm(inns_gbif, inns_obis, inns_nbn); gc()

##=====================================================================
## 7. MAP 3 — INDICATOR SPECIES OCCURRENCE RECORDS
##=====================================================================

cat("\n========== MAP 3: INDICATOR SPECIES ==========\n\n")

indicator_names <- unique(indicator_list$latin_name)

ind_gbif <- query_gbif(indicator_names, wkt = scamp_wkt)
ind_obis <- query_obis(indicator_names, wkt = scamp_wkt)
ind_nbn  <- query_nbn(indicator_names, bbox = scamp_bbox)

indicator_records <- bind_rows(ind_gbif, ind_obis, ind_nbn) %>%
  filter(!is.na(latitude), !is.na(longitude)) %>%
  deduplicate_records() %>%
  left_join(
    indicator_list %>% select(latin_name, common_name, habitat),
    by = "latin_name"
  )

if (nrow(indicator_records) > 0) {
  indicator_sf <- indicator_records %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
  st_write(indicator_sf, "outputs/species/indicator_records.gpkg",
           delete_dsn = TRUE)
}

cat(sprintf("Indicators: %d records across %d species\n",
    nrow(indicator_records), n_distinct(indicator_records$latin_name)))

## Indicator summary — flag species with no records
indicator_summary <- indicator_list %>%
  left_join(
    indicator_records %>% count(latin_name, name = "n_records"),
    by = "latin_name"
  ) %>%
  mutate(
    n_records = replace_na(n_records, 0L),
    flag = if_else(n_records == 0,
                   "NO RECORDS — review inclusion", "")
  ) %>%
  arrange(habitat, desc(n_records))

cat("\nFlagged indicator species (zero records in Solway):\n")
indicator_summary %>%
  filter(n_records == 0) %>%
  select(common_name, latin_name, habitat, flag) %>%
  print(n = Inf)

rm(ind_gbif, ind_obis, ind_nbn); gc()

##=====================================================================
## 8. MAP 4 — PMF SPECIES OCCURRENCE RECORDS
##=====================================================================

cat("\n========== MAP 4: PMF SPECIES ==========\n\n")

pmf_names <- unique(pmf_list$latin_name)

pmf_gbif <- query_gbif(pmf_names, wkt = scamp_wkt)
pmf_obis <- query_obis(pmf_names, wkt = scamp_wkt)
pmf_nbn  <- query_nbn(pmf_names, bbox = scamp_bbox)

pmf_records <- bind_rows(pmf_gbif, pmf_obis, pmf_nbn) %>%
  filter(!is.na(latitude), !is.na(longitude)) %>%
  deduplicate_records() %>%
  left_join(pmf_list %>% select(latin_name, common_name), by = "latin_name")

if (nrow(pmf_records) > 0) {
  pmf_sf <- pmf_records %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
  st_write(pmf_sf, "outputs/species/pmf_records.gpkg", delete_dsn = TRUE)
}

cat(sprintf("PMF: %d records across %d species\n",
    nrow(pmf_records), n_distinct(pmf_records$latin_name)))

## PMF summary
pmf_summary <- pmf_list %>%
  left_join(
    pmf_records %>% count(latin_name, name = "n_records"),
    by = "latin_name"
  ) %>%
  mutate(n_records = replace_na(n_records, 0L)) %>%
  arrange(desc(n_records))

cat("\nPMF species in Solway:\n")
print(pmf_summary, n = Inf)

rm(pmf_gbif, pmf_obis, pmf_nbn); gc()

##=====================================================================
## 9. COMBINED SUMMARY TABLE
##=====================================================================

cat("\n========== COMBINED SUMMARY ==========\n\n")

combined_summary <- bind_rows(
  inns_summary %>% mutate(habitat = "Multiple") %>%
    select(latin_name, common_name, habitat, category, n_records),
  indicator_summary %>%
    select(latin_name, common_name, habitat, category, n_records),
  pmf_summary %>% mutate(habitat = "Marine") %>%
    select(latin_name, common_name, habitat, category, n_records)
) %>%
  ## Some species appear in multiple lists (e.g., oyster is indicator + PMF)
  ## Keep the row with highest record count per species per category
  group_by(latin_name, category) %>%
  slice_max(n_records, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(category, habitat, desc(n_records))

write_csv(combined_summary, "outputs/species/species_summary.csv")
cat("Saved: outputs/species/species_summary.csv\n")

## Final report
cat(sprintf(
  "\n===== COMPLETE =====
  Richness grid:  outputs/species/richness_grid.gpkg
  INNS records:   outputs/species/inns_records.gpkg (%d records)
  Indicators:     outputs/species/indicator_records.gpkg (%d records)
  PMF records:    outputs/species/pmf_records.gpkg (%d records)
  Summary table:  outputs/species/species_summary.csv

  Species with NO records in Solway: %d (review before finalising monitoring plan)\n",
  nrow(inns_records), nrow(indicator_records), nrow(pmf_records),
  sum(combined_summary$n_records == 0)
))
