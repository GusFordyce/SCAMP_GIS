##=====================================================================
## SCAMP GIS Analysis — Theme D: Species
## 02_species.R
##=====================================================================
## Produces four outputs:
##   1. Species richness per grid cell — marine and terrestrial
##      (richness, effort, temporal recency [median + most recent],
##      sampling completeness, beta diversity, community similarity,
##      plus per-phylum taxonomic breakdown)
##   2. INNS occurrence records
##   3. Indicator species occurrence records (per habitat)
##   4. Priority Marine Features (PMF) occurrence records
##
## Data sources: GBIF, OBIS, NBN Atlas
## Grid: 5 km (appropriate for mixed-precision records)
##
## Outputs (all in outputs/species/):
##   - richness_grid_marine.gpkg       (marine metrics per cell)
##   - richness_grid_terrestrial.gpkg  (terrestrial metrics per cell)
##   - inns_records.gpkg               (INNS occurrence points)
##   - indicator_records.gpkg          (indicator species points)
##   - pmf_records.gpkg                (PMF species points)
##   - species_summary.csv             (presence/absence per species)
##=====================================================================
install.packages('worrms')

require(rgbif)
require(robis)
require(galah)
require(sf)
require(terra)
require(dplyr)
require(tidyr)
require(readr)
require(stringr)
require(vegan)
require(cowplot)
if (!requireNamespace("worrms", quietly = TRUE)) install.packages("worrms")
require(worrms)

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

smr_solway_bng <- st_transform(smr_solway, 27700)
las_dg_bng     <- st_transform(las_dg, 27700)

## Buffer marine region by 5 miles and intersect with D&G to get coastal land strip
## the provided land strip had no CRS associated with it
smr_buf       <- st_buffer(smr_solway_bng, dist = 8047)
coastal_strip <- st_intersection(st_geometry(smr_buf), st_geometry(las_dg_bng))

scamp_aoi    <- st_union(c(st_geometry(smr_solway_bng), coastal_strip))
scamp_aoi_sf <- st_sf(geometry = scamp_aoi)

dir.create("data/processed/boundaries", recursive = TRUE, showWarnings = FALSE) #save to processed data
st_write(st_sf(geometry = scamp_aoi), "data/processed/boundaries/scamp_aoi.gpkg",
         delete_dsn = TRUE)

plot(st_geometry(las_dg_bng), col = "grey80", border = "grey40")
plot(scamp_aoi, col = adjustcolor("red", alpha.f = 0.4), border = "red", add = TRUE)

Sys.setenv(SHAPE_RESTORE_SHX = "YES")

biosphere <- st_read("data/raw/boundaries/Biosphere_SCAMP_Zone.shp", quiet = TRUE)
mcnabb    <- st_read("data/raw/boundaries/DG_McNabb_5milezone.shp",  quiet = TRUE)

cat("=== Biosphere SCAMP Zone ===\n")
str(biosphere) #linestring only

cat("\n=== McNabb 5 Mile Zone ===\n")
str(mcnabb) #no crs

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
plot(st_geometry(gsab_scamp)) #biosphere polygon created

## --- Clip remaining boundary datasets to SCAMP AOI ---
boundary_gpkg <- list(
  ibas            = "data/raw/boundaries/IBAs_UK_-8129902679595446195.gpkg",
  rspb_reserves   = "data/raw/boundaries/RSPB_Public_Reserves_7326421649958477753.gpkg",
  mpas            = "data/raw/boundaries/Marine_Protected_Areas_-6371676791122909031.gpkg",
  ramsar          = "data/raw/boundaries/RAMSAR_Wetlands_of_International_Importance_-1649623871716453965.gpkg",
  sssi            = "data/raw/boundaries/Sites_of_Special_Scientific_Interest_8447690424398606363.gpkg",
  sac             = "data/raw/boundaries/Special_Areas_of_Conservation_5743281901804056965.gpkg",
  nnr             = "data/raw/boundaries/National_Nature_Reserves_1039238937738926456.gpkg",
  spa             = "data/raw/boundaries/Special_Protection_Areas_4781773662199779435.gpkg",
  nature_reserves = "data/raw/boundaries/Nature_Reserves_1930300340402281996.gpkg",
  lnr             = "data/raw/boundaries/Local_Nature_Reserves_7017289517571999860.gpkg",
  harbour_limits  = "data/raw/boundaries/area_management_harbour_limits/area_management_harbour_limitsPolygon.shp"
)

scamp_boundaries <- list()

for (nm in names(boundary_gpkg)) {
  cat(sprintf("Reading %s...\n", nm))
  tryCatch({
    raw <- st_read(boundary_gpkg[[nm]], quiet = TRUE) %>%
      st_zm(drop = TRUE, what = "ZM") %>%
      st_transform(27700) %>%
      st_make_valid()

    hits <- lengths(st_intersects(raw, scamp_aoi)) > 0
    if (sum(hits) == 0) {
      cat(sprintf("  %s: no features intersect SCAMP AOI — skipped\n", nm))
      next
    }

    clipped <- st_intersection(raw[hits, ], scamp_aoi)
    scamp_boundaries[[nm]] <- clipped

    out_path <- sprintf("data/processed/boundaries/%s_scamp.gpkg", nm)
    st_write(clipped, out_path, delete_dsn = TRUE, quiet = TRUE)
    cat(sprintf("  %s: %d features clipped → %s\n", nm, nrow(clipped), out_path))
  }, error = function(e) {
    cat(sprintf("  %s: ERROR — %s\n", nm, e$message))
  })
}

cat(sprintf("\nBoundary clipping complete: %d/%d datasets intersect SCAMP AOI\n",
    length(scamp_boundaries), length(boundary_gpkg)))

## --- Publication map: SCAMP AOI & designated areas ---
require(ggplot2)
install.packages('ggspatial')
require(ggspatial)

designation_labels <- c(
  ibas            = "Important Bird Areas",
  rspb_reserves   = "RSPB Reserves",
  mpas            = "Marine Protected Areas",
  ramsar          = "Ramsar Wetlands",
  sssi            = "Sites of Special Scientific Interest",
  sac             = "Special Areas of Conservation",
  nnr             = "National Nature Reserves",
  spa             = "Special Protection Areas",
  nature_reserves = "Nature Reserves",
  lnr             = "Local Nature Reserves",
  harbour_limits  = "Harbour Limits"
)

designation_colours <- c(
  "Important Bird Areas"              = "#E69F00",
  "RSPB Reserves"                     = "#56B4E9",
  "Marine Protected Areas"            = "#0072B2",
  "Ramsar Wetlands"                   = "#009E73",
  "Sites of Special Scientific Interest" = "#F0E442",
  "Special Areas of Conservation"     = "#CC79A7",
  "National Nature Reserves"          = "#D55E00",
  "Special Protection Areas"          = "#332288",
  "Nature Reserves"                   = "#88CCEE",
  "Local Nature Reserves"             = "#44AA99",
  "Harbour Limits"                    = "#999999",
  "Galloway & Southern Ayrshire Biosphere" = "#882255"
)

designations_sf <- bind_rows(
  lapply(names(scamp_boundaries), function(nm) {
    st_sf(
      designation = designation_labels[[nm]],
      geometry = st_geometry(scamp_boundaries[[nm]])
    )
  })
)

if (nrow(gsab_scamp) > 0) {
  designations_sf <- bind_rows(
    designations_sf,
    st_sf(
      designation = "Galloway & Southern Ayrshire Biosphere",
      geometry = st_geometry(gsab_scamp)
    )
  )
}

desig_fill    <- designations_sf %>% filter(designation != "Galloway & Southern Ayrshire Biosphere")
desig_outline <- designations_sf %>% filter(designation == "Galloway & Southern Ayrshire Biosphere")

p_designations <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey90", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = desig_fill, aes(fill = designation),
          colour = NA, alpha = 0.45) +
  geom_sf(data = desig_outline, fill = NA, colour = "#882255",
          linewidth = 0.6, linetype = "dashed") +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black",
          linewidth = 0.6, linetype = "solid") +
  scale_fill_manual(values = designation_colours, name = "Designation") +
  annotation_scale(location = "bl", width_hint = 0.25,
                   style = "ticks", line_width = 0.5) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         height = unit(1, "cm"), width = unit(1, "cm"),
                         style = north_arrow_minimal) +
  labs(title = "SCAMP Project Area — Designated Sites",
       subtitle = "Solway Coast and Marine Project, Dumfries & Galloway") +
  theme_minimal(base_size = 10) +
  theme(
    legend.position  = "right",
    legend.key.size  = unit(0.4, "cm"),
    legend.text      = element_text(size = 7),
    legend.title     = element_text(size = 8, face = "bold"),
    plot.title       = element_text(size = 12, face = "bold"),
    plot.subtitle    = element_text(size = 9, colour = "grey30"),
    panel.grid.major = element_line(colour = "grey85", linewidth = 0.2),
    axis.text        = element_text(size = 7)
  ) +
  coord_sf(crs = 27700)

dir.create("figures", showWarnings = FALSE)
ggsave("figures/scamp_designated_areas.pdf", p_designations,
       width = 16, height = 10, units = "in")
ggsave("figures/scamp_designated_areas.png", p_designations,
       width = 16, height = 10, units = "in", dpi = 300)

cat("Saved: figures/scamp_designated_areas.pdf & .png\n")
print(p_designations)



## --- Setting up APIs for species---
galah_config(atlas = "United Kingdom", email = Sys.getenv("gus.fordyce@noc.ac.uk"))

## GBIF credentials — reads GBIF_USER, GBIF_PWD, GBIF_EMAIL from .Renviron.
## To set these, run usethis::edit_r_environ() and add:
##   GBIF_USER=your_username
##   GBIF_PWD=your_password
##   GBIF_EMAIL=your_email
## Then restart R for the changes to take effect.
usethis::edit_r_environ()

options(
  gbif_user  = Sys.getenv("GBIF_USER"),
  gbif_pwd   = Sys.getenv("GBIF_PWD"),
  gbif_email = Sys.getenv("GBIF_EMAIL")
)

## Create output directory
dir.create("data/processed/species", showWarnings = FALSE, recursive = TRUE)

## 5 km grid in OSGB36 derived from SCAMP AOI
scamp_extent <- ext(st_bbox(scamp_aoi)[c("xmin", "xmax", "ymin", "ymax")])
grid_template <- rast(scamp_extent, res = 5000, crs = "EPSG:27700")

scamp_bbox <- st_bbox(st_transform(scamp_aoi, 4326))
scamp_wkt <- st_as_text(st_as_sfc(scamp_bbox)) #for API queries

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

cat(sprintf("Species lists compiled: %d INNS, %d indicators\n",
    nrow(inns_list), nrow(indicator_list)))

##=====================================================================
## 3. QA/QC CONFIGURATION
##=====================================================================
## Species occurrence records from online databases often contain errors
## that can distort spatial analysis — for example, records plotted at
## country centroids instead of the actual observation location, or
## coordinates that have been swapped or negated during data entry.
##
## We apply quality control in three stages:
##
##   Stage 1 — At query time, we ask each database to exclude records it
##     has already flagged as having geographic problems. This catches the
##     most obvious errors before the data even reaches us.
##
##   Stage 2 — After downloading, we apply our own filters: removing
##     records that are too old (pre-1970), too imprecise (location
##     uncertainty > 10 km), or from inappropriate sources (e.g. fossils,
##     zoo specimens rather than wild observations).
##
##   Stage 3 — We run the CoordinateCleaner package, which checks for
##     common data artefacts that slip through database-level QC. For
##     example, records sitting exactly on a capital city or a museum
##     are likely geocoding errors, not real observations.
##
## Every record we exclude is logged with the reason, so the decisions
## are auditable and we can revisit thresholds if needed.
##
## References:
##   GBIF data quality: https://data-blog.gbif.org/post/gbif-filtering-guide/
##   OBIS QC pipeline:  https://manual.obis.org/dataquality.html
##   CoordinateCleaner: Zizka et al. 2019, Methods Ecol Evol 10:744-751

install.packages("CoordinateCleaner", quiet = TRUE)
require(CoordinateCleaner)

## -- GBIF issue flags --
## GBIF attaches machine-readable issue codes to each record. These nine
## flags all indicate the coordinates cannot be trusted — the location
## might be (0,0), outside the stated country, or mangled during
## reprojection from the original coordinate system.
gbif_bad_flags <- c(
  "ZERO_COORDINATE",
  "COORDINATE_OUT_OF_RANGE",
  "COORDINATE_INVALID",
  "COUNTRY_COORDINATE_MISMATCH",
  "COORDINATE_REPROJECTION_FAILED",
  "COORDINATE_REPROJECTION_SUSPICIOUS",
  "PRESUMED_NEGATED_LATITUDE",
  "PRESUMED_NEGATED_LONGITUDE",
  "PRESUMED_SWAPPED_COORDINATE"
)
gbif_bad_pattern <- paste(gbif_bad_flags, collapse = "|")

## -- OBIS quality flags --
## OBIS runs its own QC pipeline and attaches flags in a similar way.
## We exclude records with missing or impossible coordinates, and those
## flagged as being on land (since OBIS is a marine database, an on-land
## flag means the coordinates are almost certainly wrong).
obis_bad_flags <- c("NO_COORD", "ZERO_COORD", "LON_OUT_OF_RANGE",
                     "LAT_OUT_OF_RANGE", "ON_LAND")
obis_bad_pattern <- paste(obis_bad_flags, collapse = "|")

## -- Basis of record --
## GBIF records come from different sources. Fossils tell us about
## historical distribution, not current presence. Living specimens
## (zoos, aquaria, botanic gardens) reflect captive locations, not
## wild populations. We keep everything else: human observations,
## machine observations, preserved museum specimens, etc.
bor_exclude <- c("FOSSIL_SPECIMEN", "LIVING_SPECIMEN")

## -- Bogus coordinate uncertainty values --
## Some data publishers use placeholder values (e.g. 301m, 3036m) for
## coordinate uncertainty instead of leaving the field blank. These are
## not real measurements, so we treat them the same as genuinely
## imprecise records and exclude them.
bogus_uncertainty <- c(301, 3036, 999, 9999)

## -- QA exclusion log --
## Every record we remove gets added here with a reason code. At the
## end of the script, write_qa_log() saves this to a CSV so we have
## a full audit trail of what was excluded and why.
qa_log <- list()

## Helper function: takes a dataset and a TRUE/FALSE mask indicating
## which records to keep. Records that fail the mask are logged with
## a reason code and removed. Returns the filtered dataset.
qa_filter <- function(df, keep_mask, reason, source_db) {
  n_excluded <- sum(!keep_mask)
  if (n_excluded > 0) {
    excluded <- df[!keep_mask, ] %>%
      mutate(qa_reason = reason, qa_source = source_db)
    qa_log[[length(qa_log) + 1]] <<- excluded
    cat(sprintf("    QA [%s]: %d records excluded (%s)\n",
        source_db, n_excluded, reason))
  }
  df[keep_mask, ]
}

##=====================================================================
## 3a. QUERY FUNCTIONS
##=====================================================================

## --- GBIF (Global Biodiversity Information Facility) ---
##
## GBIF is the largest aggregator of species occurrence records worldwide.
## It collects data from museums, citizen science platforms (e.g. iRecord),
## and national databases. Many NBN Atlas records also flow into GBIF, so
## we deduplicate later to avoid double-counting.
##
## Two modes:
##   - Bulk (species_names = NULL): downloads ALL species in the study area
##     using the async occ_download() endpoint, which has no record cap.
##     GBIF applies our QA filters server-side before building the download.
##     Used for the all-taxa richness map.
##
##   - Per-species: queries one species at a time using occ_data(), which
##     returns up to 'limit' records per species. QA filters are applied
##     after download since occ_data() has fewer server-side options.
##     Used for INNS and indicator species maps.
##
query_gbif <- function(species_names = NULL, wkt, limit = 500) {
  cat("Querying GBIF...\n")

  if (is.null(species_names)) {
    ## -- Bulk download --
    ## We ask GBIF to build a download file containing all georeferenced
    ## records within the study area WKT polygon. The predicates tell GBIF
    ## to pre-filter: only records with coordinates, no known geospatial
    ## issues, confirmed present (not absence records), and excluding
    ## fossils and captive specimens.
    cat("  Bulk download (occ_download, server-side QA)...\n")
    dl_key <- rgbif::occ_download(
      pred("hasCoordinate", TRUE),
      pred("hasGeospatialIssue", FALSE),
      pred("occurrenceStatus", "PRESENT"),
      pred_not(pred_in("basisOfRecord", bor_exclude)),
      pred_within(wkt),
      format = "SIMPLE_CSV"
    )

    ## GBIF queues the download and prepares it asynchronously. This
    ## typically takes 2-10 minutes depending on how many records match.
    cat(sprintf("  Download key: %s — waiting for GBIF to prepare...\n", dl_key))
    rgbif::occ_download_wait(dl_key, curlopts = list())

    ## Retrieve the prepared zip file and import it into R as a data frame.
    dl_path <- rgbif::occ_download_get(dl_key, overwrite = TRUE)
    recs    <- rgbif::occ_download_import(dl_path)

    if (is.null(recs) || nrow(recs) == 0) return(tibble())
    cat(sprintf("  %d raw records downloaded\n", nrow(recs)))

    ## Stage 2: even after server-side filtering, we apply additional
    ## checks that GBIF doesn't offer as predicates.

    ## Remove records where the stated location precision is worse than
    ## 10 km, or where the uncertainty value is a known placeholder.
    ## Records with no uncertainty value are kept — many legitimate
    ## records simply lack this metadata.
    recs <- qa_filter(recs,
      is.na(recs$coordinateUncertaintyInMeters) |
        (recs$coordinateUncertaintyInMeters < 10000 &
         !recs$coordinateUncertaintyInMeters %in% bogus_uncertainty),
      "coord_uncertainty_>10km_or_bogus", "GBIF_bulk")

    ## Remove records from before 1970 or with no date. Older records
    ## are less relevant to current species distributions and more
    ## likely to have imprecise or estimated coordinates.
    recs <- qa_filter(recs,
      !is.na(recs$year) & recs$year >= 1970,
      "year_pre1970_or_missing", "GBIF_bulk")

    ## Keep the columns we need — core schema shared across databases
    ## plus extra fields for marine/terrestrial splitting and metrics.
    out <- recs %>%
      transmute(
        latin_name = species,
        longitude = decimalLongitude,
        latitude = decimalLatitude,
        year = year,
        coord_uncertainty = coordinateUncertaintyInMeters,
        basis_of_record = basisOfRecord,
        kingdom = kingdom,
        phylum = phylum,
        class = class,
        source = "GBIF"
      ) %>%
      filter(!is.na(latin_name))
    cat(sprintf("  %d records after QA\n", nrow(out)))
    return(out)
  }

  ## -- Per-species queries --
  ## Loop through each species name and query GBIF individually.
  ## This uses the simpler occ_data() endpoint which returns results
  ## immediately but has a per-query record limit.
  results <- list()
  for (i in seq_along(species_names)) {
    sp <- species_names[i]
    cat(sprintf("  [%d/%d] %s", i, length(species_names), sp))
    tryCatch({
      ## Request records for this species within the study area polygon.
      ## hasCoordinate = TRUE excludes records with no location.
      recs <- rgbif::occ_data(
        scientificName = sp,
        geometry = wkt,
        hasCoordinate = TRUE,
        limit = limit
      )
      if (!is.null(recs$data) && nrow(recs$data) > 0) {
        d <- recs$data
        n_raw <- nrow(d)

        ## Check GBIF's issue flags — remove any record that has been
        ## flagged with a coordinate problem (e.g. country mismatch,
        ## zero coordinates, reprojection failure).
        if ("issues" %in% names(d)) {
          d <- qa_filter(d,
            !grepl(gbif_bad_pattern, d$issues),
            "gbif_issue_flag", "GBIF")
        }

        ## Remove fossils and zoo/aquarium specimens.
        if ("basisOfRecord" %in% names(d)) {
          d <- qa_filter(d,
            !d$basisOfRecord %in% bor_exclude,
            "basis_of_record_excluded", "GBIF")
        }

        ## Remove records with location precision worse than 10 km
        ## or with known placeholder uncertainty values.
        if ("coordinateUncertaintyInMeters" %in% names(d)) {
          d <- qa_filter(d,
            is.na(d$coordinateUncertaintyInMeters) |
              (d$coordinateUncertaintyInMeters < 10000 &
               !d$coordinateUncertaintyInMeters %in% bogus_uncertainty),
            "coord_uncertainty_>10km_or_bogus", "GBIF")
        }

        ## Remove records from before 1970 or with no date.
        d <- qa_filter(d,
          !is.na(d$year) & d$year >= 1970,
          "year_pre1970_or_missing", "GBIF")

        ## Standardise to common schema with extra fields.
        results[[sp]] <- d %>%
          transmute(
            latin_name = species,
            longitude = decimalLongitude,
            latitude = decimalLatitude,
            year = year,
            coord_uncertainty = if ("coordinateUncertaintyInMeters" %in% names(d))
              coordinateUncertaintyInMeters else NA_real_,
            basis_of_record = if ("basisOfRecord" %in% names(d))
              basisOfRecord else NA_character_,
            kingdom = if ("kingdom" %in% names(d)) kingdom else NA_character_,
            phylum = if ("phylum" %in% names(d)) phylum else NA_character_,
            class = if ("class" %in% names(d)) class else NA_character_,
            source = "GBIF"
          )
        cat(sprintf(" — %d/%d records after QA\n", nrow(results[[sp]]), n_raw))
      } else {
        cat(" — 0 records\n")
      }
    }, error = function(e) cat(sprintf(" — error: %s\n", e$message)))
  }
  bind_rows(results)
}

## --- OBIS (Ocean Biodiversity Information System) ---
##
## OBIS is the marine equivalent of GBIF — it aggregates species occurrence
## records from ocean surveys, monitoring programmes, and research cruises.
## OBIS runs its own QC pipeline and attaches quality flags to each record.
##
## Two modes (same as GBIF): bulk for all taxa, per-species for target lists.
##
## QA is handled by the internal .obis_qa() helper, which:
##   1. Checks OBIS quality flags — removes records flagged as having
##      missing/impossible coordinates or being located on land (which
##      for a marine database means the location is wrong).
##   2. Removes records from before 1970 or with no date.
##   3. Standardises columns to the common schema.
##
query_obis <- function(species_names = NULL, wkt, limit = 500) {
  cat("Querying OBIS...\n")

  ## Internal helper that applies QA to any set of OBIS records.
  ## Used by both the bulk and per-species paths below.
  .obis_qa <- function(recs, label) {
    if (is.null(recs) || nrow(recs) == 0) return(tibble())
    n_raw <- nrow(recs)

    ## Check the OBIS quality flags column. Remove records flagged with
    ## missing coordinates, impossible coordinate values, or on-land
    ## locations (which indicate a geocoding error for marine records).
    if ("flags" %in% names(recs)) {
      recs <- qa_filter(recs,
        is.na(recs$flags) | !grepl(obis_bad_pattern, recs$flags),
        "obis_spatial_flag", label)
    }

    ## Remove records from before 1970 or with no date.
    if ("date_year" %in% names(recs)) {
      recs <- qa_filter(recs,
        !is.na(recs$date_year) & recs$date_year >= 1970,
        "year_pre1970_or_missing", label)
    }

    ## Standardise to the common schema shared with GBIF and NBN.
    recs %>%
      transmute(
        latin_name = species,
        longitude = decimalLongitude,
        latitude = decimalLatitude,
        year = date_year,
        coord_uncertainty = if ("coordinateUncertaintyInMeters" %in% names(recs))
          coordinateUncertaintyInMeters else NA_real_,
        basis_of_record = if ("basisOfRecord" %in% names(recs))
          basisOfRecord else NA_character_,
        kingdom = if ("kingdom" %in% names(recs)) kingdom else NA_character_,
        phylum = if ("phylum" %in% names(recs)) phylum else NA_character_,
        class = if ("class" %in% names(recs)) class else NA_character_,
        source = "OBIS"
      ) %>%
      filter(!is.na(latin_name))
  }

  if (is.null(species_names)) {
    ## Bulk query: download all species records within the study area.
    cat("  Bulk query (all taxa)...\n")
    recs <- robis::occurrence(geometry = wkt)
    out <- .obis_qa(recs, "OBIS_bulk")
    cat(sprintf("  %d records after QA\n", nrow(out)))
    return(out)
  }

  ## Per-species queries: loop through each species individually.
  results <- list()
  for (i in seq_along(species_names)) {
    sp <- species_names[i]
    cat(sprintf("  [%d/%d] %s", i, length(species_names), sp))
    tryCatch({
      ## Request all OBIS records for this species within the study area.
      recs <- robis::occurrence(scientificname = sp, geometry = wkt)
      if (!is.null(recs) && nrow(recs) > 0) {
        ## Cap at 'limit' records per species to keep query times manageable.
        recs <- head(recs, limit)
        n_raw <- nrow(recs)
        ## Apply the same QA checks as the bulk path.
        results[[sp]] <- .obis_qa(recs, "OBIS")
        cat(sprintf(" — %d/%d records after QA\n", nrow(results[[sp]]), n_raw))
      } else {
        cat(" — 0 records\n")
      }
    }, error = function(e) cat(sprintf(" — error: %s\n", e$message)))
  }
  bind_rows(results)
}

## --- NBN Atlas (National Biodiversity Network, via galah package) ---
##
## The NBN Atlas holds UK-specific occurrence records, including many
## records from local recording groups and county recorders that may
## not appear in GBIF. However, a large proportion of NBN records also
## feed into GBIF, so we deduplicate across databases later.
##
## NBN Atlas does not support practical bulk downloads without a species
## filter, so this function only works in per-species mode. It is used
## for the INNS and indicator species queries, not the all-taxa richness
## map (which relies on GBIF + OBIS).
##
## QA approach:
##   - Server-side: the galah package lets us exclude records that NBN
##     has flagged with coordinate problems (invalid, zero, or country
##     mismatch) before they are downloaded.
##   - Post-download: we extract the year from the event date and remove
##     records from before 1970.
##
query_nbn <- function(species_names = NULL, bbox) {
  cat("Querying NBN Atlas...\n")

  if (is.null(species_names)) {
    ## Bulk download is not practical through the NBN Atlas API without
    ## a species filter — the query would time out or return incomplete
    ## results. The all-taxa richness map uses GBIF + OBIS instead.
    cat("  Skipping NBN for bulk query (use GBIF + OBIS)\n")
    return(tibble())
  }

  ## Per-species queries: loop through each species and query NBN.
  results <- list()
  for (i in seq_along(species_names)) {
    sp <- species_names[i]
    cat(sprintf("  [%d/%d] %s", i, length(species_names), sp))
    tryCatch({
      ## Build the query: search for this species within the bounding box,
      ## and tell NBN to exclude records it has flagged with coordinate
      ## problems before sending them to us.
      recs <- galah_call() %>%
        identify(sp) %>%
        filter(
          decimalLongitude >= bbox["xmin"],
          decimalLongitude <= bbox["xmax"],
          decimalLatitude  >= bbox["ymin"],
          decimalLatitude  <= bbox["ymax"]
        ) %>%
        galah_filter(
          assertions != "COORDINATE_INVALID",
          assertions != "ZERO_COORDINATE",
          assertions != "COUNTRY_COORDINATE_MISMATCH"
        ) %>%
        atlas_occurrences()

      if (!is.null(recs) && nrow(recs) > 0) {
        n_raw <- nrow(recs)

        ## Extract the observation year from the event date string.
        recs <- recs %>%
          mutate(yr = as.integer(substr(eventDate, 1, 4)))

        ## Remove records from before 1970 or with no date.
        recs <- qa_filter(recs,
          !is.na(recs$yr) & recs$yr >= 1970,
          "year_pre1970_or_missing", "NBN")

        ## Standardise to the common schema with extra fields.
        results[[sp]] <- recs %>%
          transmute(
            latin_name = scientificName,
            longitude = decimalLongitude,
            latitude = decimalLatitude,
            year = yr,
            coord_uncertainty = if ("coordinateUncertaintyInMeters" %in% names(recs))
              coordinateUncertaintyInMeters else NA_real_,
            basis_of_record = if ("basisOfRecord" %in% names(recs))
              basisOfRecord else NA_character_,
            kingdom = if ("kingdom" %in% names(recs)) kingdom else NA_character_,
            phylum = if ("phylum" %in% names(recs)) phylum else NA_character_,
            class = if ("class" %in% names(recs)) class else NA_character_,
            source = "NBN"
          )
        cat(sprintf(" — %d/%d records after QA\n", nrow(results[[sp]]), n_raw))
      } else {
        cat(" — 0 records\n")
      }
    }, error = function(e) cat(sprintf(" — error: %s\n", e$message)))
  }
  bind_rows(results)
}

##=====================================================================
## 3b. SPATIAL QA (CoordinateCleaner)
##=====================================================================
## After downloading and filtering on database flags, we run a final
## spatial check using the CoordinateCleaner package. This catches a
## common class of errors that database-level QC often misses:
## records whose coordinates are technically valid but are actually
## geocoding artefacts rather than real observation locations.
##
## The tests we run:
##   - capitals:     Is the record suspiciously close to a capital city?
##                   Many records get assigned to a country or province
##                   capital when the true location is unknown.
##   - centroids:    Same idea, but for country/province geographic
##                   centres — another common default for unknown locations.
##   - equal:        Does latitude equal longitude? This almost never
##                   happens in nature and usually indicates a data entry
##                   error or placeholder.
##   - gbif:         Is the record at GBIF headquarters in Copenhagen?
##                   A surprising number of records end up here due to
##                   default coordinate values in data processing.
##   - institutions: Is the record at a known museum, herbarium, or
##                   university? These are likely specimen storage
##                   locations, not where the organism was observed.
##   - zeros:        Is the record near (0, 0) in the Gulf of Guinea?
##                   This is the most common coordinate default worldwide.
##
## Tests we deliberately skip:
##   - seas:   Flags points in the ocean. We skip this because our study
##             area spans both land and sea — marine records are expected.
##   - outliers: Flags geographic outliers per species. We skip this
##             because many of our target species have small sample sizes
##             in the Solway, and the test needs at least 7 records per
##             species to work reliably.

cc_clean <- function(df) {
  if (nrow(df) == 0) return(df)
  cat("  Running CoordinateCleaner...\n")

  ## CoordinateCleaner expects columns named 'decimalLongitude',
  ## 'decimalLatitude', and 'species' — we create these from our
  ## standardised column names.
  cc_input <- df %>%
    mutate(decimalLongitude = longitude, decimalLatitude = latitude,
           species = latin_name)

  ## Run all six tests. The result is the original data frame with
  ## extra TRUE/FALSE columns for each test, plus a .summary column
  ## that is TRUE if the record passed ALL tests.
  cc_result <- clean_coordinates(
    cc_input,
    lon = "decimalLongitude",
    lat = "decimalLatitude",
    species = "species",
    tests = c("capitals", "centroids", "equal", "gbif",
              "institutions", "zeros"),
    value = "spatialvalid",
    verbose = FALSE
  )

  ## Records where .summary is FALSE failed at least one test.
  flagged <- !cc_result$.summary
  n_flagged <- sum(flagged)

  if (n_flagged > 0) {
    ## Log which specific test each record failed, so we can see in the
    ## QA log whether most exclusions are e.g. centroid hits vs institution
    ## hits — this helps judge whether thresholds need adjusting.
    test_cols <- c(".cap", ".cen", ".equ", ".gbf", ".inst", ".zer")
    test_names <- c("near_capital", "near_centroid", "equal_lat_lon",
                     "near_gbif_hq", "near_institution", "near_zero")
    for (j in seq_along(test_cols)) {
      if (test_cols[j] %in% names(cc_result)) {
        test_fail <- !cc_result[[test_cols[j]]]
        if (sum(test_fail) > 0) {
          qa_log[[length(qa_log) + 1]] <<- df[test_fail, ] %>%
            mutate(qa_reason = paste0("cc_", test_names[j]),
                   qa_source = "CoordinateCleaner")
        }
      }
    }
    cat(sprintf("    CoordinateCleaner: %d records flagged\n", n_flagged))
  }

  ## Return only the records that passed all tests.
  df[!flagged, ]
}

##=====================================================================
## 4. DEDUPLICATION & QA LOG
##=====================================================================
## The same observation can appear in multiple databases — for example,
## a record submitted to the NBN Atlas often also flows into GBIF. If
## we don't remove these duplicates, the same sighting would be counted
## multiple times in our species richness and abundance analyses.
##
## We identify duplicates by rounding coordinates to 3 decimal places
## (~100 m precision) and matching on species name + rounded location +
## year. This is slightly conservative — it may merge genuinely separate
## records that happen to be very close — but that is preferable to
## inflating counts with database-level duplicates.

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

## Writes the accumulated QA exclusion log to CSV. This is called once
## at the end of the script after all maps have been produced. The CSV
## contains every record that was excluded during processing, along with
## the reason code and which database it came from. This makes the QA
## decisions fully auditable — if a species appears absent and we want
## to check whether records existed but were filtered out, this is
## where to look.
write_qa_log <- function() {
  if (length(qa_log) == 0) {
    cat("QA log: no records excluded\n")
    return(invisible(NULL))
  }
  qa_all <- bind_rows(qa_log) %>%
    select(any_of(c("latin_name", "longitude", "latitude", "year",
                     "source", "qa_reason", "qa_source")))
  write_csv(qa_all, "data/processed/species/qa_excluded.csv")
  cat(sprintf("QA log: %d excluded records written to data/processed/species/qa_excluded.csv\n",
      nrow(qa_all)))

  ## Print a summary table showing how many records were excluded
  ## by each QA step, so we can quickly see which filters are having
  ## the biggest impact.
  qa_summary <- qa_all %>% count(qa_source, qa_reason, name = "n_excluded")
  cat("\nQA exclusion summary:\n")
  print(qa_summary, n = Inf)
}

##=====================================================================
## 5. MAP 1 — SPECIES RICHNESS PER GRID CELL (MARINE / TERRESTRIAL)
##=====================================================================
## Bulk download all taxa from GBIF + OBIS, split into marine and
## terrestrial, calculate 7 per-cell metrics, produce paired plots.

cat("\n========== MAP 1: SPECIES RICHNESS ==========\n\n")

richness_gbif <- query_gbif(species_names = NULL, wkt = scamp_wkt)
richness_obis <- query_obis(species_names = NULL, wkt = scamp_wkt)

richness_all <- bind_rows(richness_gbif, richness_obis) %>%
  filter(!is.na(latitude), !is.na(longitude), !is.na(latin_name)) %>%
  cc_clean() %>%
  deduplicate_records()

##=====================================================================
## 5a. TAXONOMIC NAME CLEANING
##=====================================================================
## Resolve species names to accepted taxonomy before calculating any
## richness or diversity metrics. Marine names go through WoRMS,
## terrestrial names through the GBIF backbone.

cat("\n--- Taxonomic name cleaning ---\n")
n_before_taxclean <- nrow(richness_all)
unique_names <- unique(richness_all$latin_name)
cat(sprintf("Unique species names to resolve: %d\n", length(unique_names)))

## We need a preliminary spatial split to decide which names go to
## WoRMS vs GBIF backbone. Use a quick point-in-polygon test.
tmp_pts <- st_as_sf(richness_all,
                    coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(27700)
tmp_in_sea  <- lengths(st_intersects(tmp_pts, smr_solway_bng)) > 0
tmp_in_land <- lengths(st_intersects(tmp_pts, las_dg_bng)) > 0
richness_all$.tmp_realm <- case_when(
  tmp_in_sea  ~ "marine",
  tmp_in_land ~ "terrestrial",
  TRUE        ~ "terrestrial"
)
rm(tmp_pts, tmp_in_sea, tmp_in_land)

marine_names <- richness_all %>%
  filter(.tmp_realm == "marine") %>%
  pull(latin_name) %>%
  unique()
terrestrial_names <- richness_all %>%
  filter(.tmp_realm == "terrestrial") %>%
  pull(latin_name) %>%
  unique()

## --- WoRMS resolution for marine names ---
cat(sprintf("Resolving %d marine names via WoRMS...\n", length(marine_names)))
worms_lookup <- tibble(original = character(), accepted = character(),
                       aphia_id = integer(), worms_status = character())

batch_size <- 50
for (i in seq(1, length(marine_names), by = batch_size)) {
  batch <- marine_names[i:min(i + batch_size - 1, length(marine_names))]
  cat(sprintf("  WoRMS batch %d–%d of %d\n",
      i, min(i + batch_size - 1, length(marine_names)), length(marine_names)))
  tryCatch({
    res <- worrms::wm_records_names(batch, marine_only = FALSE)
    for (j in seq_along(batch)) {
      r <- res[[j]]
      if (is.null(r) || nrow(r) == 0 || all(is.na(r$valid_AphiaID))) {
        worms_lookup <- bind_rows(worms_lookup, tibble(
          original = batch[j], accepted = batch[j],
          aphia_id = NA_integer_, worms_status = "unmatched"
        ))
      } else {
        best <- r %>%
          filter(status %in% c("accepted", "valid") |
                 !is.na(valid_AphiaID)) %>%
          slice(1)
        if (nrow(best) == 0) best <- slice(r, 1)
        worms_lookup <- bind_rows(worms_lookup, tibble(
          original = batch[j],
          accepted = if (!is.na(best$valid_name)) best$valid_name else batch[j],
          aphia_id = best$valid_AphiaID,
          worms_status = best$status
        ))
      }
    }
    Sys.sleep(0.5)
  }, error = function(e) {
    cat(sprintf("  WoRMS batch error: %s\n", e$message))
    for (nm in batch) {
      worms_lookup <<- bind_rows(worms_lookup, tibble(
        original = nm, accepted = nm,
        aphia_id = NA_integer_, worms_status = "error"
      ))
    }
  })
}

## --- GBIF backbone resolution for terrestrial names ---
cat(sprintf("Resolving %d terrestrial names via GBIF backbone...\n",
    length(terrestrial_names)))
gbif_lookup <- tibble(original = character(), accepted = character(),
                      gbif_key = integer(), gbif_status = character(),
                      gbif_rank = character())

for (i in seq_along(terrestrial_names)) {
  nm <- terrestrial_names[i]
  if (i %% 200 == 0) cat(sprintf("  GBIF backbone: %d / %d\n",
                                   i, length(terrestrial_names)))
  tryCatch({
    bb <- rgbif::name_backbone(name = nm, rank = "species",
                                kingdom = NULL, strict = FALSE)
    if (!is.null(bb) && "species" %in% names(bb) && !is.na(bb$species)) {
      gbif_lookup <- bind_rows(gbif_lookup, tibble(
        original = nm,
        accepted = bb$species,
        gbif_key = bb$usageKey,
        gbif_status = bb$status,
        gbif_rank = bb$rank
      ))
    } else if (!is.null(bb) && "canonicalName" %in% names(bb) &&
               !is.na(bb$canonicalName)) {
      gbif_lookup <- bind_rows(gbif_lookup, tibble(
        original = nm,
        accepted = bb$canonicalName,
        gbif_key = bb$usageKey,
        gbif_status = bb$status,
        gbif_rank = if ("rank" %in% names(bb)) bb$rank else NA_character_
      ))
    } else {
      gbif_lookup <- bind_rows(gbif_lookup, tibble(
        original = nm, accepted = nm,
        gbif_key = NA_integer_, gbif_status = "unmatched",
        gbif_rank = NA_character_
      ))
    }
  }, error = function(e) {
    gbif_lookup <<- bind_rows(gbif_lookup, tibble(
      original = nm, accepted = nm,
      gbif_key = NA_integer_, gbif_status = "error",
      gbif_rank = NA_character_
    ))
  })
}

## --- Build combined lookup and apply to records ---
combined_lookup <- bind_rows(
  worms_lookup %>% transmute(original, accepted,
                              source_db = "WoRMS", status = worms_status),
  gbif_lookup %>% transmute(original, accepted,
                             source_db = "GBIF_backbone", status = gbif_status)
) %>%
  distinct(original, .keep_all = TRUE)

## Flag synonyms, unmatched, and genus-only records
combined_lookup <- combined_lookup %>%
  mutate(
    is_synonym  = (original != accepted & status != "unmatched" & status != "error"),
    is_unmatched = (status %in% c("unmatched", "error"))
  )

## Save the lookup table for reference
write_csv(combined_lookup, "data/processed/species/taxonomic_lookup.csv")

## Apply accepted names to records
richness_all <- richness_all %>%
  left_join(combined_lookup %>% select(original, accepted),
            by = c("latin_name" = "original")) %>%
  mutate(
    original_name = latin_name,
    latin_name = if_else(!is.na(accepted), accepted, latin_name)
  ) %>%
  select(-accepted)

## Drop genus-only records (single-word names with no species epithet)
genus_only <- !grepl(" ", richness_all$latin_name)
n_genus_dropped <- sum(genus_only)
richness_all <- richness_all %>% filter(!genus_only)

## Re-deduplicate after name changes (synonyms may now match)
n_after_rename <- nrow(richness_all)
richness_all <- deduplicate_records(richness_all)

## Summary
n_synonyms  <- sum(combined_lookup$is_synonym, na.rm = TRUE)
n_unmatched <- sum(combined_lookup$is_unmatched, na.rm = TRUE)

cat(sprintf(
  "\nTaxonomic cleaning summary:
  Records before:       %d
  Genus-only dropped:   %d
  Synonyms resolved:    %d unique names
  Unmatched names:      %d (kept as-is, flagged)
  Records after:        %d
  Lookup saved:         data/processed/species/taxonomic_lookup.csv\n",
  n_before_taxclean, n_genus_dropped, n_synonyms,
  n_unmatched, nrow(richness_all)
))

## Clean up temporary column
richness_all <- richness_all %>% select(-.tmp_realm)

## --- Marine / terrestrial classification (spatial split) ---
## Classify each record by whether it falls on land (D&G boundary)
## or sea (Solway Marine Region). Records on the exact coastline
## (inside both polygons) are duplicated into both realms.
richness_pts <- st_as_sf(richness_all,
                         coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(27700)

in_land <- lengths(st_intersects(richness_pts, las_dg_bng)) > 0
in_sea  <- lengths(st_intersects(richness_pts, smr_solway_bng)) > 0

richness_all$realm <- case_when(
  in_land & in_sea  ~ "coastal",
  in_land           ~ "terrestrial",
  in_sea            ~ "marine",
  TRUE              ~ NA_character_
)

n_coastal <- sum(richness_all$realm == "coastal", na.rm = TRUE)
cat(sprintf("Spatial split: %d marine, %d terrestrial, %d coastal (included in both)\n",
    sum(richness_all$realm == "marine", na.rm = TRUE),
    sum(richness_all$realm == "terrestrial", na.rm = TRUE),
    n_coastal))

## Coastal records contribute to both realms — duplicate them
coastal_as_marine <- richness_all %>%
  filter(realm == "coastal") %>%
  mutate(realm = "marine")
coastal_as_terr <- richness_all %>%
  filter(realm == "coastal") %>%
  mutate(realm = "terrestrial")

richness_all <- richness_all %>%
  filter(realm %in% c("marine", "terrestrial")) %>%
  bind_rows(coastal_as_marine, coastal_as_terr)

rm(richness_pts, in_land, in_sea, coastal_as_marine, coastal_as_terr)

cat(sprintf("After coastal duplication: %d marine, %d terrestrial\n",
    sum(richness_all$realm == "marine"),
    sum(richness_all$realm == "terrestrial")))

## --- Coastline extraction ---
## Intersect D&G boundary with SCAMP AOI to get the coastline.
scamp_coastline <- st_intersection(
  st_cast(st_boundary(las_dg_bng), "MULTILINESTRING"),
  st_transform(scamp_aoi_sf, 27700)
)
st_write(scamp_coastline, "data/processed/boundaries/scamp_coastline.gpkg",
         delete_dsn = TRUE, quiet = TRUE)

## --- Build grid (once, shared by both realms) ---
richness_sf <- richness_all %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

grid_osgb <- st_make_grid(
  st_transform(st_as_sfc(st_bbox(richness_sf)), 27700),
  cellsize = 5000, square = TRUE
) %>%
  st_sf(cell_id = seq_along(.), geometry = .) %>%
  st_transform(4326)

records_gridded <- st_join(richness_sf, grid_osgb, join = st_within) %>%
  st_drop_geometry()

## --- Save cleaned species data for sharing ---
write_csv(richness_all, "data/processed/species/all_species_records.csv")
st_write(richness_sf, "data/processed/species/all_species_records.gpkg",
         delete_dsn = TRUE, quiet = TRUE)
cat(sprintf("Saved %d cleaned records to data/processed/species/all_species_records.csv & .gpkg\n",
    nrow(richness_all)))

## --- Queen-contiguity neighbour lookup for beta diversity ---
grid_bng <- st_transform(grid_osgb, 27700)
grid_touches <- st_touches(grid_bng)

## --- Helper: compute metrics for one realm ---
## Metrics: richness, effort, temporal (median year + most recent),
## sampling completeness, beta diversity, community similarity.
compute_grid_metrics <- function(df, realm_label) {
  df_realm <- df %>%
    filter(realm == realm_label, !is.na(cell_id),
           !is.na(latin_name), latin_name != "")

  if (nrow(df_realm) == 0) {
    cat(sprintf("  No %s records — skipping\n", realm_label))
    return(NULL)
  }

  ## Species richness
  richness <- df_realm %>%
    group_by(cell_id) %>%
    summarise(species_richness = n_distinct(latin_name), .groups = "drop")

  ## Recording effort (total, GBIF/OBIS split)
  effort <- df_realm %>%
    group_by(cell_id) %>%
    summarise(
      record_count = n(),
      gbif_records = sum(source == "GBIF"),
      obis_records = sum(source == "OBIS"),
      .groups = "drop"
    )

  ## Temporal recency (median year + most recent record)
  temporal <- df_realm %>%
    group_by(cell_id) %>%
    summarise(
      median_year     = median(year, na.rm = TRUE),
      most_recent_year = max(year, na.rm = TRUE),
      .groups = "drop"
    )

  ## Species-by-cell matrix (shared by completeness, beta, similarity)
  sp_cell <- df_realm %>%
    count(cell_id, latin_name, name = "abundance")
  sp_wide <- sp_cell %>%
    pivot_wider(id_cols = cell_id, names_from = latin_name,
                values_from = abundance, values_fill = 0L)
  cell_ids_wide <- sp_wide$cell_id
  sp_mat <- as.matrix(sp_wide[, -1])

  ## Sampling completeness (rarefied richness at 10th percentile)
  row_totals <- rowSums(sp_mat)
  rarefy_n <- max(1L, as.integer(quantile(row_totals, 0.10)))
  rarefied <- rep(NA_real_, nrow(sp_mat))
  can_rarefy <- row_totals >= rarefy_n
  if (any(can_rarefy)) {
    rarefied[can_rarefy] <- vegan::rarefy(sp_mat[can_rarefy, , drop = FALSE],
                                           sample = rarefy_n)
  }
  completeness <- tibble(
    cell_id = cell_ids_wide,
    rarefied_richness = rarefied,
    rarefy_n = rarefy_n
  )

  ## Beta diversity (mean Sorensen to queen-contiguous neighbours)
  pa_mat <- (sp_mat > 0) * 1L
  row_lookup <- setNames(seq_along(cell_ids_wide), cell_ids_wide)
  beta_vals <- rep(NA_real_, nrow(pa_mat))
  for (i in seq_len(nrow(pa_mat))) {
    cid <- cell_ids_wide[i]
    nb_cids <- grid_touches[[cid]]
    nb_rows <- na.omit(row_lookup[as.character(nb_cids)])
    if (length(nb_rows) == 0) next
    pair_mat <- pa_mat[c(i, nb_rows), , drop = FALSE]
    if (nrow(pair_mat) < 2) next
    d <- vegan::vegdist(pair_mat, method = "bray", binary = TRUE)
    beta_vals[i] <- mean(as.matrix(d)[1, -1])
  }
  beta_div <- tibble(cell_id = cell_ids_wide, beta_diversity = beta_vals)

  ## Community composition similarity (mean Jaccard to all other cells)
  if (nrow(pa_mat) >= 2) {
    jaccard_dist <- vegan::vegdist(pa_mat, method = "bray", binary = TRUE)
    jaccard_sim  <- 1 - as.matrix(jaccard_dist)
    diag(jaccard_sim) <- NA
    mean_similarity <- rowMeans(jaccard_sim, na.rm = TRUE)

    saveRDS(jaccard_sim,
            sprintf("outputs/species/community_similarity_matrix_%s.rds",
                    realm_label))
    cat(sprintf("  Saved %s Jaccard similarity matrix (%d x %d)\n",
        realm_label, nrow(jaccard_sim), ncol(jaccard_sim)))
  } else {
    mean_similarity <- NA_real_
  }
  community_sim <- tibble(
    cell_id = cell_ids_wide,
    community_similarity = mean_similarity
  )

  ## Richness per taxonomic group (top phyla as separate columns)
  taxon_richness <- df_realm %>%
    filter(!is.na(phylum), phylum != "") %>%
    group_by(cell_id, phylum) %>%
    summarise(spp = n_distinct(latin_name), .groups = "drop") %>%
    pivot_wider(id_cols = cell_id, names_from = phylum,
                values_from = spp, values_fill = 0L,
                names_prefix = "richness_")

  ## Assemble all metrics
  metrics <- richness %>%
    left_join(effort, by = "cell_id") %>%
    left_join(temporal, by = "cell_id") %>%
    left_join(completeness, by = "cell_id") %>%
    left_join(beta_div, by = "cell_id") %>%
    left_join(community_sim, by = "cell_id") %>%
    left_join(taxon_richness, by = "cell_id")

  ## Join to grid polygons
  grid_sf <- grid_osgb %>%
    left_join(metrics, by = "cell_id") %>%
    filter(!is.na(species_richness))

  grid_sf
}

## --- Compute metrics for both realms ---
grid_marine      <- compute_grid_metrics(records_gridded, "marine")
grid_terrestrial <- compute_grid_metrics(records_gridded, "terrestrial")

## --- Save outputs ---
if (!is.null(grid_marine)) {
  st_write(grid_marine, "outputs/species/richness_grid_marine.gpkg",
           delete_dsn = TRUE)
  cat(sprintf("Marine grid: %d cells, max richness = %d\n",
      nrow(grid_marine), max(grid_marine$species_richness)))
}
if (!is.null(grid_terrestrial)) {
  st_write(grid_terrestrial, "outputs/species/richness_grid_terrestrial.gpkg",
           delete_dsn = TRUE)
  cat(sprintf("Terrestrial grid: %d cells, max richness = %d\n",
      nrow(grid_terrestrial), max(grid_terrestrial$species_richness)))
}

## --- Shared theme and map furniture for all plots ---
scamp_aoi_bng  <- st_transform(scamp_aoi_sf, 27700)
smr_solway_bng <- st_transform(smr_solway, 27700)

map_theme <- theme_minimal(base_size = 10) +
  theme(
    legend.position  = "right",
    plot.title       = element_text(size = 12, face = "bold"),
    plot.subtitle    = element_text(size = 9, colour = "grey30"),
    panel.grid.major = element_line(colour = "grey85", linewidth = 0.2),
    axis.text        = element_text(size = 7)
  )

## Marine panels clip to SMR boundary, terrestrial to D&G boundary
make_panel <- function(grid_sf, fill_var, fill_label, title, subtitle,
                       realm, scale_type = "viridis",
                       viridis_opt = "viridis", direction = -1,
                       midpoint = NULL, low = NULL, high = NULL) {
  grid_bng_plot <- st_transform(grid_sf, 27700)

  clip_boundary <- if (realm == "marine") smr_solway_bng else las_dg_bng
  bg_fill  <- if (realm == "marine") "lightblue" else "grey90"
  bg_color <- if (realm == "marine") "steelblue" else "grey50"

  p <- ggplot() +
    geom_sf(data = clip_boundary, fill = bg_fill, colour = bg_color,
            linewidth = 0.3) +
    geom_sf(data = grid_bng_plot, aes(fill = .data[[fill_var]]),
            colour = NA) +
    geom_sf(data = scamp_coastline, colour = "grey30", linewidth = 0.3) +
    geom_sf(data = scamp_aoi_bng, fill = NA, colour = "black",
            linewidth = 0.5)

  if (scale_type == "diverging" && !is.null(midpoint)) {
    p <- p + scale_fill_gradient2(
      name = fill_label, midpoint = midpoint,
      low = if (is.null(low)) "#2166ac" else low,
      mid = "white",
      high = if (is.null(high)) "#b2182b" else high,
      na.value = "grey80"
    )
  } else {
    p <- p + scale_fill_viridis_c(option = viridis_opt, name = fill_label,
                                   direction = direction, na.value = "grey80")
  }

  p +
    annotation_scale(location = "bl", width_hint = 0.25,
                     style = "ticks", line_width = 0.5) +
    annotation_north_arrow(location = "tr", which_north = "true",
                           height = unit(1, "cm"), width = unit(1, "cm"),
                           style = north_arrow_minimal) +
    labs(title = title, subtitle = subtitle) +
    map_theme +
    coord_sf(crs = 27700)
}

## Helper to save paired marine/terrestrial plots
save_paired <- function(p_marine, p_terrestrial, filename, w = 16, h = 10) {
  p_combined <- plot_grid(p_marine, p_terrestrial, ncol = 2,
                          labels = c("Marine", "Terrestrial"),
                          label_size = 11)
  ggsave(paste0("figures/", filename, ".pdf"), p_combined,
         width = w, height = h, units = "in")
  ggsave(paste0("figures/", filename, ".png"), p_combined,
         width = w, height = h, units = "in", dpi = 300)
  cat(sprintf("Saved: figures/%s.pdf & .png\n", filename))
  print(p_combined)
}

## --- 1. Species richness ---
if (!is.null(grid_marine) && !is.null(grid_terrestrial)) {
  p_rich_m <- make_panel(grid_marine, "species_richness",
                         "Species\nrichness", "Marine", "GBIF + OBIS",
                         realm = "marine")
  p_rich_t <- make_panel(grid_terrestrial, "species_richness",
                         "Species\nrichness", "Terrestrial", "GBIF + OBIS",
                         realm = "terrestrial")
  save_paired(p_rich_m, p_rich_t, "map1a_richness")
}

## --- 2. Recording effort ---
if (!is.null(grid_marine) && !is.null(grid_terrestrial)) {
  p_eff_m <- make_panel(grid_marine, "record_count",
                        "Records", "Marine — Recording Effort",
                        "Total records per cell",
                        realm = "marine", viridis_opt = "inferno")
  p_eff_t <- make_panel(grid_terrestrial, "record_count",
                        "Records", "Terrestrial — Recording Effort",
                        "Total records per cell",
                        realm = "terrestrial", viridis_opt = "inferno")
  save_paired(p_eff_m, p_eff_t, "map1b_effort")
}

## --- 3a. Temporal recency — median recording year ---
## Diverging scale centred on 2010: blue = older, red = more recent
if (!is.null(grid_marine) && !is.null(grid_terrestrial)) {
  p_med_m <- make_panel(grid_marine, "median_year",
                        "Median\nyear", "Marine — Median Recording Year",
                        "Diverging around 2010",
                        realm = "marine", scale_type = "diverging",
                        midpoint = 2010)
  p_med_t <- make_panel(grid_terrestrial, "median_year",
                        "Median\nyear", "Terrestrial — Median Recording Year",
                        "Diverging around 2010",
                        realm = "terrestrial", scale_type = "diverging",
                        midpoint = 2010)
  save_paired(p_med_m, p_med_t, "map1c_median_year")
}

## --- 3b. Temporal recency — most recent record ---
## Sequential scale: darker = more recently recorded
if (!is.null(grid_marine) && !is.null(grid_terrestrial)) {
  p_rec_m <- make_panel(grid_marine, "most_recent_year",
                        "Most recent\nyear", "Marine — Most Recent Record",
                        "Year of latest observation per cell",
                        realm = "marine", viridis_opt = "plasma")
  p_rec_t <- make_panel(grid_terrestrial, "most_recent_year",
                        "Most recent\nyear", "Terrestrial — Most Recent Record",
                        "Year of latest observation per cell",
                        realm = "terrestrial", viridis_opt = "plasma")
  save_paired(p_rec_m, p_rec_t, "map1d_most_recent")
}

## --- 4. Sampling completeness (terrestrial only) ---
if (!is.null(grid_terrestrial)) {
  p_rare_t <- make_panel(grid_terrestrial, "rarefied_richness",
                         "Rarefied\nrichness", "Terrestrial — Sampling Completeness",
                         sprintf("Rarefied at n = %d",
                                 grid_terrestrial$rarefy_n[1]),
                         realm = "terrestrial", viridis_opt = "mako")
  ggsave("figures/map1e_sampling_completeness.pdf", p_rare_t,
         width = 10, height = 10, units = "in")
  ggsave("figures/map1e_sampling_completeness.png", p_rare_t,
         width = 10, height = 10, units = "in", dpi = 300)
  cat("Saved: figures/map1e_sampling_completeness.pdf & .png\n")
  print(p_rare_t)
}

## --- 5. Beta diversity (mean Jaccard to neighbours) ---
if (!is.null(grid_marine) && !is.null(grid_terrestrial)) {
  p_beta_m <- make_panel(grid_marine, "beta_diversity",
                         "Beta\ndiversity", "Marine — Beta Diversity",
                         "Mean Sorensen dissimilarity to neighbours",
                         realm = "marine", viridis_opt = "turbo")
  p_beta_t <- make_panel(grid_terrestrial, "beta_diversity",
                         "Beta\ndiversity", "Terrestrial — Beta Diversity",
                         "Mean Sorensen dissimilarity to neighbours",
                         realm = "terrestrial", viridis_opt = "turbo")
  save_paired(p_beta_m, p_beta_t, "map1f_beta_diversity")
}

## --- 6. Community composition similarity ---
if (!is.null(grid_marine) && !is.null(grid_terrestrial)) {
  p_sim_m <- make_panel(grid_marine, "community_similarity",
                        "Mean\nsimilarity", "Marine — Community Similarity",
                        "Mean Sorensen similarity to all other cells",
                        realm = "marine", viridis_opt = "cividis")
  p_sim_t <- make_panel(grid_terrestrial, "community_similarity",
                        "Mean\nsimilarity", "Terrestrial — Community Similarity",
                        "Mean Sorensen similarity to all other cells",
                        realm = "terrestrial", viridis_opt = "cividis")
  save_paired(p_sim_m, p_sim_t, "map1g_community_similarity")
}

## --- 7. Taxonomic group breakdown ---
## Identify the top phyla in each realm and produce per-phylum richness maps
make_taxon_panels <- function(grid_sf, df_records, realm_label) {
  realm_records <- df_records %>%
    filter(realm == realm_label, !is.na(phylum), phylum != "")
  top_phyla <- realm_records %>%
    count(phylum, sort = TRUE) %>%
    head(6) %>%
    pull(phylum)

  richness_cols <- paste0("richness_", top_phyla)
  present_cols <- richness_cols[richness_cols %in% names(grid_sf)]

  if (length(present_cols) == 0) return(invisible(NULL))

  panels <- lapply(present_cols, function(col) {
    phylum_name <- sub("^richness_", "", col)
    make_panel(grid_sf, col, "Species",
               phylum_name, paste(realm_label, "richness"),
               realm = realm_label, viridis_opt = "viridis")
  })

  n_panels <- length(panels)
  ncols <- min(3L, n_panels)
  p_grid <- plot_grid(plotlist = panels, ncol = ncols)
  fname <- sprintf("map1h_taxon_breakdown_%s", realm_label)
  h <- if (n_panels <= 3) 8 else 16
  ggsave(paste0("figures/", fname, ".pdf"), p_grid,
         width = 16, height = h, units = "in")
  ggsave(paste0("figures/", fname, ".png"), p_grid,
         width = 16, height = h, units = "in", dpi = 300)
  cat(sprintf("Saved: figures/%s.pdf & .png (%d phyla)\n", fname, n_panels))
  print(p_grid)
}

if (!is.null(grid_marine)) {
  make_taxon_panels(grid_marine, records_gridded, "marine")
}
if (!is.null(grid_terrestrial)) {
  make_taxon_panels(grid_terrestrial, records_gridded, "terrestrial")
}

rm(richness_gbif, richness_obis, richness_all, richness_sf,
   records_gridded, grid_marine, grid_terrestrial)
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
  cc_clean() %>%
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

## Plot INNS occurrences — points coloured by species on SCAMP AOI background
if (nrow(inns_records) > 0) {
  p_inns <- ggplot() +
    geom_sf(data = las_dg_bng, fill = "grey90", colour = "grey50", linewidth = 0.3) +
    geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
    geom_sf(data = inns_sf, aes(colour = common_name), size = 1.5, alpha = 0.7) +
    scale_colour_brewer(palette = "Set3", name = "INNS Species") +
    annotation_scale(location = "bl", width_hint = 0.25,
                     style = "ticks", line_width = 0.5) +
    annotation_north_arrow(location = "tr", which_north = "true",
                           height = unit(1, "cm"), width = unit(1, "cm"),
                           style = north_arrow_minimal) +
    labs(title = "INNS Occurrence Records",
         subtitle = "Invasive Non-Native Species in SCAMP AOI") +
    theme_minimal(base_size = 10) +
    theme(
      legend.position  = "right",
      legend.key.size  = unit(0.4, "cm"),
      legend.text      = element_text(size = 7),
      plot.title       = element_text(size = 12, face = "bold"),
      plot.subtitle    = element_text(size = 9, colour = "grey30"),
      panel.grid.major = element_line(colour = "grey85", linewidth = 0.2)
    ) +
    coord_sf(crs = 27700)

  ggsave("figures/map2_inns_occurrences.pdf", p_inns,
         width = 16, height = 10, units = "in")
  ggsave("figures/map2_inns_occurrences.png", p_inns,
         width = 16, height = 10, units = "in", dpi = 300)
  cat("Saved: figures/map2_inns_occurrences.pdf & .png\n")
  print(p_inns)
}

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
  cc_clean() %>%
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

## Plot indicator species — faceted by habitat type so each panel shows
## occurrences for one habitat (seagrass, oyster, saltmarsh, etc.)
if (nrow(indicator_records) > 0) {
  p_indicators <- ggplot() +
    geom_sf(data = las_dg_bng, fill = "grey90", colour = "grey50", linewidth = 0.3) +
    geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5) +
    geom_sf(data = indicator_sf, aes(colour = common_name), size = 1.2, alpha = 0.7) +
    facet_wrap(~habitat, ncol = 3) +
    annotation_scale(location = "bl", width_hint = 0.2,
                     style = "ticks", line_width = 0.5) +
    labs(title = "Indicator Species Occurrence Records",
         subtitle = "By habitat type — SCAMP AOI",
         colour = "Species") +
    theme_minimal(base_size = 9) +
    theme(
      legend.position  = "bottom",
      legend.key.size  = unit(0.3, "cm"),
      legend.text      = element_text(size = 6),
      plot.title       = element_text(size = 12, face = "bold"),
      plot.subtitle    = element_text(size = 9, colour = "grey30"),
      panel.grid.major = element_line(colour = "grey85", linewidth = 0.2),
      strip.text       = element_text(size = 9, face = "bold")
    ) +
    guides(colour = guide_legend(ncol = 4, override.aes = list(size = 2.5))) +
    coord_sf(crs = 27700)

  ggsave("figures/map3_indicator_species.pdf", p_indicators,
         width = 16, height = 12, units = "in")
  ggsave("figures/map3_indicator_species.png", p_indicators,
         width = 16, height = 12, units = "in", dpi = 300)
  cat("Saved: figures/map3_indicator_species.pdf & .png\n")
  print(p_indicators)
}

rm(ind_gbif, ind_obis, ind_nbn); gc()

##=====================================================================
## 8. MAP 4 — PMF SPECIES (from GeMS data)
##=====================================================================

cat("\n========== MAP 4: PMF SPECIES (GeMS) ==========\n\n")

gems_points <- read_csv(
  "data/raw/species/GeMS___Scottish_Priority_Marine_Features_-species points.csv",
  show_col_types = FALSE
)
gems_polys <- st_read(
  "data/raw/species/GeMS___Scottish_Priority_Marine_Features_-species_polygons.gpkg",
  quiet = TRUE
)

cat("GeMS species points:", nrow(gems_points), "\n")
cat("GeMS species polygons:", nrow(gems_polys), "\n")
str(gems_points, max.level = 1)
str(gems_polys, max.level = 1)

scamp_aoi_wgs <- st_transform(st_sf(geometry = scamp_aoi), 4326)

## Points — filter to those with coords, convert to sf, clip to SCAMP AOI
gems_pts_sf <- gems_points %>%
  filter(!is.na(LATITUDE), !is.na(LONGITUDE)) %>%
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = 4326) %>%
  st_intersection(scamp_aoi_wgs)

## Polygons — transform and clip to SCAMP AOI
gems_polys_clip <- gems_polys %>%
  st_transform(4326) %>%
  st_make_valid() %>%
  {
    hits <- lengths(st_intersects(., scamp_aoi_wgs)) > 0
    st_intersection(.[hits, ], scamp_aoi_wgs)
  }

pmf_records <- bind_rows(
  gems_pts_sf %>% st_drop_geometry() %>%
    transmute(latin_name = SPECIES_NAME, common_name = COMMON_NAME, source = "GeMS_points"),
  gems_polys_clip %>% st_drop_geometry() %>%
    transmute(latin_name = SPECIES_NAME, common_name = COMMON_NAME, source = "GeMS_polygons")
)

cat(sprintf("PMF in SCAMP AOI: %d records (%d points, %d polygon features) across %d species\n",
    nrow(pmf_records), nrow(gems_pts_sf), nrow(gems_polys_clip),
    n_distinct(pmf_records$latin_name)))

## Save clipped outputs
st_write(gems_pts_sf, "data/processed/species/pmf_points_scamp.gpkg",
         delete_dsn = TRUE, quiet = TRUE)
if (nrow(gems_polys_clip) > 0) {
  st_write(gems_polys_clip, "data/processed/species/pmf_polygons_scamp.gpkg",
           delete_dsn = TRUE, quiet = TRUE)
}

## PMF summary
pmf_summary <- pmf_records %>%
  count(latin_name, common_name, name = "n_records") %>%
  mutate(category = "PMF") %>%
  arrange(desc(n_records))

cat("\nPMF species in Solway (GeMS):\n")
print(pmf_summary, n = Inf)

## Plot PMF species — points and polygons from the GeMS dataset.
## Polygons shown as semi-transparent fills, points overlaid on top.
p_pmf <- ggplot() +
  geom_sf(data = las_dg_bng, fill = "grey90", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = scamp_aoi_sf, fill = NA, colour = "black", linewidth = 0.5)

if (exists("gems_polys_clip") && nrow(gems_polys_clip) > 0) {
  gems_polys_bng <- st_transform(gems_polys_clip, 27700)
  p_pmf <- p_pmf +
    geom_sf(data = gems_polys_bng, fill = adjustcolor("#0072B2", alpha.f = 0.3),
            colour = "#0072B2", linewidth = 0.3)
}

if (exists("gems_pts_sf") && nrow(gems_pts_sf) > 0) {
  gems_pts_bng <- st_transform(gems_pts_sf, 27700)
  p_pmf <- p_pmf +
    geom_sf(data = gems_pts_bng, aes(colour = SPECIES_NAME),
            size = 1.2, alpha = 0.7)
}

p_pmf <- p_pmf +
  annotation_scale(location = "bl", width_hint = 0.25,
                   style = "ticks", line_width = 0.5) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         height = unit(1, "cm"), width = unit(1, "cm"),
                         style = north_arrow_minimal) +
  labs(title = "Priority Marine Features (PMF) — GeMS Records",
       subtitle = "Species points and habitat polygons within SCAMP AOI",
       colour = "PMF Species") +
  theme_minimal(base_size = 10) +
  theme(
    legend.position  = "right",
    legend.key.size  = unit(0.3, "cm"),
    legend.text      = element_text(size = 6),
    plot.title       = element_text(size = 12, face = "bold"),
    plot.subtitle    = element_text(size = 9, colour = "grey30"),
    panel.grid.major = element_line(colour = "grey85", linewidth = 0.2)
  ) +
  guides(colour = guide_legend(ncol = 2, override.aes = list(size = 2.5))) +
  coord_sf(crs = 27700)

ggsave("figures/map4_pmf_gems.pdf", p_pmf,
       width = 16, height = 10, units = "in")
ggsave("figures/map4_pmf_gems.png", p_pmf,
       width = 16, height = 10, units = "in", dpi = 300)
cat("Saved: figures/map4_pmf_gems.pdf & .png\n")
print(p_pmf)

rm(gems_points, gems_polys, gems_pts_sf, gems_polys_clip); gc()

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
    select(any_of(c("latin_name", "common_name", "habitat", "category", "n_records")))
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
  Marine grid:    outputs/species/richness_grid_marine.gpkg
  Terrestrial:    outputs/species/richness_grid_terrestrial.gpkg
  Coastline:      data/processed/boundaries/scamp_coastline.gpkg
  INNS records:   outputs/species/inns_records.gpkg (%d records)
  Indicators:     outputs/species/indicator_records.gpkg (%d records)
  PMF points:     data/processed/species/pmf_points_scamp.gpkg
  PMF polygons:   data/processed/species/pmf_polygons_scamp.gpkg
  PMF total:      %d records
  Summary table:  outputs/species/species_summary.csv

  Species with NO records in Solway: %d (review before finalising monitoring plan)\n",
  nrow(inns_records), nrow(indicator_records), nrow(pmf_records),
  sum(combined_summary$n_records == 0)
))

## Write QA exclusion log
write_qa_log()
