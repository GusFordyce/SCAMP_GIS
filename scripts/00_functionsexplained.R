##########################################################################
##  00_functionsexplained.R
##  SCAMP GIS — Annotated reference for all custom functions
##
##  This script documents every function defined in 02_themeD_species.R.
##  Each function is reproduced with line-by-line comments explaining
##  what the code does, why it does it, and how each step connects to
##  the ones before and after it.
##
##  This file is for reference only — it is not intended to be sourced.
##########################################################################


## =====================================================================
##  FUNCTIONS IN SPECIES SCRIPT (02_themeD_species.R)
## =====================================================================


## ---------------------------------------------------------------------
##  1. qa_filter
## ---------------------------------------------------------------------
##  SUMMARY: Applies a logical mask to a data frame, keeping rows that
##  pass (TRUE) and logging rows that fail (FALSE) into a global QA
##  exclusion list. This is the single gateway through which all QA
##  exclusions flow — every record removed from the pipeline passes
##  through this function, so the audit trail is complete.
##
##  Parameters:
##    df         — data frame of occurrence records
##    keep_mask  — logical vector (same length as nrow(df)); TRUE = keep
##    reason     — character string labelling why records were removed
##    source_db  — which database the records came from (GBIF, OBIS, NBN)
##
##  Returns: the subset of df where keep_mask is TRUE
## ---------------------------------------------------------------------

qa_filter <- function(df, keep_mask, reason, source_db) {
  ## Count how many records will be excluded by this filter step.
  ## This is used both for the console progress message and to decide
  ## whether we need to log anything at all.
  n_excluded <- sum(!keep_mask)

  if (n_excluded > 0) {
    ## Extract the rows that FAILED the filter — these are the records
    ## we are about to remove. We attach two metadata columns so the QA
    ## log records not just what was excluded but why and from where.
    excluded <- df[!keep_mask, ] %>%
      mutate(qa_reason = reason, qa_source = source_db)

    ## Append this batch of excluded records to the global qa_log list.
    ## The <<- operator writes to the qa_log variable in the parent
    ## (global) environment, not a local copy — this is how exclusions
    ## accumulate across multiple calls throughout the script.
    qa_log[[length(qa_log) + 1]] <<- excluded

    ## Print a progress line so the user can see which filter just ran,
    ## how many records it removed, and from which database. This is
    ## helpful during long processing runs.
    cat(sprintf("    QA [%s]: %d records excluded (%s)\n",
        source_db, n_excluded, reason))
  }

  ## Return only the rows that passed the filter. Downstream code
  ## receives a cleaned data frame and never sees the excluded records.
  df[keep_mask, ]
}


## ---------------------------------------------------------------------
##  2. query_gbif
## ---------------------------------------------------------------------
##  SUMMARY: Queries the Global Biodiversity Information Facility for
##  species occurrence records within the SCAMP study area. Operates
##  in two modes: (a) bulk download of ALL species via the asynchronous
##  occ_download() endpoint (no record cap — used for the all-taxa
##  richness map), or (b) per-species queries via occ_data() (used for
##  targeted INNS and indicator species maps). Both paths apply multi-
##  stage QA: server-side predicates where possible, then post-download
##  filters for coordinate uncertainty, date range, issue flags, and
##  basis-of-record. Returns a standardised tibble with the 10-column
##  common schema shared across all three databases.
##
##  Parameters:
##    species_names — NULL for bulk mode; character vector for per-species
##    wkt           — WKT polygon string defining the study area (WGS84)
##    limit         — max records per species in per-species mode (default 500)
##
##  Returns: tibble with columns latin_name, longitude, latitude, year,
##           coord_uncertainty, basis_of_record, kingdom, phylum, class, source
## ---------------------------------------------------------------------

query_gbif <- function(species_names = NULL, wkt, limit = 500) {
  ## Print a header so the console log shows when the GBIF query phase starts.
  cat("Querying GBIF...\n")

  if (is.null(species_names)) {
    ## ---- BULK DOWNLOAD PATH ----
    ## When no species list is provided, we want every georeferenced
    ## record inside the study area. The occ_download() endpoint is the
    ## only way to do this without hitting record caps — it asks GBIF
    ## to build a complete extract server-side.
    cat("  Bulk download (occ_download, server-side QA)...\n")

    ## Compose a set of predicates that GBIF evaluates before building
    ## the download file. This filters out problem records at the source:
    ##   - hasCoordinate TRUE: exclude records with no lat/lon
    ##   - hasGeospatialIssue FALSE: exclude records GBIF has flagged
    ##   - occurrenceStatus PRESENT: exclude absence-only records
    ##   - basisOfRecord not in bor_exclude: exclude fossils & captive
    ##   - pred_within(wkt): restrict to the study area polygon
    dl_key <- rgbif::occ_download(
      pred("hasCoordinate", TRUE),
      pred("hasGeospatialIssue", FALSE),
      pred("occurrenceStatus", "PRESENT"),
      pred_not(pred_in("basisOfRecord", bor_exclude)),
      pred_within(wkt),
      format = "SIMPLE_CSV"
    )

    ## occ_download() returns immediately with a download key — the file
    ## is prepared asynchronously on GBIF's servers. We block here until
    ## it is ready, typically 2-10 minutes.
    cat(sprintf("  Download key: %s — waiting for GBIF to prepare...\n", dl_key))
    rgbif::occ_download_wait(dl_key, curlopts = list())

    ## Once ready, download the zip file and import it as a data frame.
    ## overwrite = TRUE replaces any previous download file with the
    ## same key, avoiding stale cache issues.
    dl_path <- rgbif::occ_download_get(dl_key, overwrite = TRUE)
    recs    <- rgbif::occ_download_import(dl_path)

    ## Guard against an empty result — possible if no records match the
    ## predicates (unlikely but defensive). Return an empty tibble so
    ## downstream code doesn't need special NULL handling.
    if (is.null(recs) || nrow(recs) == 0) return(tibble())
    cat(sprintf("  %d raw records downloaded\n", nrow(recs)))

    ## ---- POST-DOWNLOAD QA (STAGE 2) ----
    ## GBIF's predicates don't cover everything we want. Two more filters:

    ## Filter 1: coordinate uncertainty.
    ## Remove records whose stated precision is worse than 10 km or whose
    ## uncertainty is a known bogus placeholder (301, 3036, 999, 9999).
    ## Records with NA uncertainty are KEPT — many legitimate records
    ## simply lack this metadata.
    recs <- qa_filter(recs,
      is.na(recs$coordinateUncertaintyInMeters) |
        (recs$coordinateUncertaintyInMeters < 10000 &
         !recs$coordinateUncertaintyInMeters %in% bogus_uncertainty),
      "coord_uncertainty_>10km_or_bogus", "GBIF_bulk")

    ## Filter 2: temporal range.
    ## Remove records from before 1970 or with no year. Older records
    ## are less relevant to current species distributions and more likely
    ## to have imprecise coordinates.
    recs <- qa_filter(recs,
      !is.na(recs$year) & recs$year >= 1970,
      "year_pre1970_or_missing", "GBIF_bulk")

    ## Standardise to the 10-column common schema. transmute() selects
    ## and renames in one step — columns not listed are dropped. We
    ## keep kingdom/phylum/class for the taxonomic breakdown maps.
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
      ## Drop records where the species name is missing — these would
      ## be uncountable in richness calculations.
      filter(!is.na(latin_name))
    cat(sprintf("  %d records after QA\n", nrow(out)))
    return(out)
  }

  ## ---- PER-SPECIES PATH ----
  ## When a species list is provided, query each one individually
  ## via occ_data(). This endpoint returns results immediately but
  ## has a per-query record limit (default 500).
  results <- list()
  for (i in seq_along(species_names)) {
    sp <- species_names[i]
    ## Print progress so the user can monitor long loops.
    cat(sprintf("  [%d/%d] %s", i, length(species_names), sp))
    tryCatch({
      ## Query GBIF for this species within the study area polygon.
      ## hasCoordinate = TRUE is applied server-side as a first filter.
      recs <- rgbif::occ_data(
        scientificName = sp,
        geometry = wkt,
        hasCoordinate = TRUE,
        limit = limit
      )
      if (!is.null(recs$data) && nrow(recs$data) > 0) {
        d <- recs$data
        n_raw <- nrow(d)

        ## QA Step 1: Check GBIF's issue flags column.
        ## Each record can carry comma-separated issue codes. If any
        ## code matches our bad-issue regex (gbif_bad_pattern), the
        ## record is excluded — it has a known coordinate problem.
        if ("issues" %in% names(d)) {
          d <- qa_filter(d,
            !grepl(gbif_bad_pattern, d$issues),
            "gbif_issue_flag", "GBIF")
        }

        ## QA Step 2: basis of record.
        ## Remove fossils and zoo/aquarium specimens, same as the
        ## bulk path's server-side predicate.
        if ("basisOfRecord" %in% names(d)) {
          d <- qa_filter(d,
            !d$basisOfRecord %in% bor_exclude,
            "basis_of_record_excluded", "GBIF")
        }

        ## QA Step 3: coordinate uncertainty.
        ## Same logic as the bulk path — exclude >10 km or bogus values.
        if ("coordinateUncertaintyInMeters" %in% names(d)) {
          d <- qa_filter(d,
            is.na(d$coordinateUncertaintyInMeters) |
              (d$coordinateUncertaintyInMeters < 10000 &
               !d$coordinateUncertaintyInMeters %in% bogus_uncertainty),
            "coord_uncertainty_>10km_or_bogus", "GBIF")
        }

        ## QA Step 4: temporal filter.
        ## Remove pre-1970 or date-less records.
        d <- qa_filter(d,
          !is.na(d$year) & d$year >= 1970,
          "year_pre1970_or_missing", "GBIF")

        ## Standardise to common schema. Conditional column extraction
        ## handles cases where optional fields are absent from the GBIF
        ## response — different species can return different column sets.
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
    ## Wrap each species query in tryCatch so one failure doesn't abort
    ## the entire loop — we just log the error and move on.
    }, error = function(e) cat(sprintf(" — error: %s\n", e$message)))
  }
  ## Combine all per-species results into a single data frame.
  ## bind_rows handles differing column sets gracefully, filling
  ## missing columns with NA.
  bind_rows(results)
}


## ---------------------------------------------------------------------
##  3. query_obis  (with internal helper .obis_qa)
## ---------------------------------------------------------------------
##  SUMMARY: Queries the Ocean Biodiversity Information System for
##  marine species occurrence records within the SCAMP study area.
##  Like query_gbif, it has bulk and per-species modes. QA is handled
##  by an internal helper function (.obis_qa) that checks OBIS quality
##  flags, applies the 1970 date filter, and standardises columns to
##  the common schema. Returns the same 10-column tibble as query_gbif.
##
##  Parameters:
##    species_names — NULL for bulk mode; character vector for per-species
##    wkt           — WKT polygon string defining the study area (WGS84)
##    limit         — max records per species in per-species mode (default 500)
##
##  Returns: tibble with columns latin_name, longitude, latitude, year,
##           coord_uncertainty, basis_of_record, kingdom, phylum, class, source
## ---------------------------------------------------------------------

query_obis <- function(species_names = NULL, wkt, limit = 500) {
  ## Console header for the OBIS query phase.
  cat("Querying OBIS...\n")

  ## -- Internal helper: .obis_qa --
  ## Applies QA filters and column standardisation to a batch of OBIS
  ## records. Extracted as a helper because the same logic is used for
  ## both the bulk and per-species paths — avoids duplicating 30+ lines.
  .obis_qa <- function(recs, label) {
    ## If the query returned nothing, return an empty tibble immediately.
    ## This prevents errors from trying to filter or transmute a NULL.
    if (is.null(recs) || nrow(recs) == 0) return(tibble())
    ## Capture the raw count before any filtering, for the progress log.
    n_raw <- nrow(recs)

    ## QA Step 1: OBIS quality flags.
    ## OBIS stores QC flags in a 'flags' column as comma-separated codes.
    ## We check for flags indicating missing coordinates, impossible
    ## coordinate values, or on-land locations. For a marine database,
    ## "on land" means the geolocation is wrong.
    if ("flags" %in% names(recs)) {
      recs <- qa_filter(recs,
        is.na(recs$flags) | !grepl(obis_bad_pattern, recs$flags),
        "obis_spatial_flag", label)
    }

    ## QA Step 2: temporal filter.
    ## Same cutoff as GBIF — remove records before 1970 or without a
    ## year. OBIS uses 'date_year' instead of 'year' as the column name.
    if ("date_year" %in% names(recs)) {
      recs <- qa_filter(recs,
        !is.na(recs$date_year) & recs$date_year >= 1970,
        "year_pre1970_or_missing", label)
    }

    ## Standardise to the common 10-column schema.
    ## OBIS calls the year column 'date_year', so we map it to 'year'.
    ## Optional columns are checked with conditional if/else to handle
    ## cases where OBIS didn't return that field.
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
      ## Remove records with missing species names, same as GBIF path.
      filter(!is.na(latin_name))
  }

  if (is.null(species_names)) {
    ## ---- BULK PATH ----
    ## Download all marine species records within the study area. The
    ## robis::occurrence() function handles pagination internally —
    ## unlike GBIF's occ_data(), it fetches all matching records.
    cat("  Bulk query (all taxa)...\n")
    recs <- robis::occurrence(geometry = wkt)
    ## Apply the shared QA helper to the full result set.
    out <- .obis_qa(recs, "OBIS_bulk")
    cat(sprintf("  %d records after QA\n", nrow(out)))
    return(out)
  }

  ## ---- PER-SPECIES PATH ----
  ## Query each species individually. This is used for INNS and
  ## indicator species maps where we have a defined species list.
  results <- list()
  for (i in seq_along(species_names)) {
    sp <- species_names[i]
    cat(sprintf("  [%d/%d] %s", i, length(species_names), sp))
    tryCatch({
      ## Query OBIS for this species within the study area.
      recs <- robis::occurrence(scientificname = sp, geometry = wkt)
      if (!is.null(recs) && nrow(recs) > 0) {
        ## Cap at 'limit' records to keep query times manageable.
        ## robis fetches all records, so we truncate after download.
        recs <- head(recs, limit)
        n_raw <- nrow(recs)
        ## Apply the same QA as the bulk path.
        results[[sp]] <- .obis_qa(recs, "OBIS")
        cat(sprintf(" — %d/%d records after QA\n", nrow(results[[sp]]), n_raw))
      } else {
        cat(" — 0 records\n")
      }
    ## One species failing shouldn't kill the loop.
    }, error = function(e) cat(sprintf(" — error: %s\n", e$message)))
  }
  ## Combine all per-species results into one data frame.
  bind_rows(results)
}


## ---------------------------------------------------------------------
##  4. query_nbn  (with internal helper .nbn_fetch)
## ---------------------------------------------------------------------
##  SUMMARY: Queries the NBN Atlas for UK-specific species occurrence
##  records using the direct REST API at records-ws.nbnatlas.org. NBN
##  holds records from local recording groups and county recorders that
##  may not be in GBIF. Only operates in per-species mode (returns empty
##  for bulk). Paginates through results 1000 at a time, applies spatial
##  and assertion filters server-side via Solr query syntax, retries on
##  rate-limiting (HTTP 412/429), and standardises output to the common
##  10-column schema. Includes a 2-second delay between species to
##  respect NBN's rate limits.
##
##  Parameters:
##    species_names — NULL returns empty; character vector for per-species
##    bbox          — named numeric vector with xmin, ymin, xmax, ymax (WGS84)
##
##  Returns: tibble with columns latin_name, longitude, latitude, year,
##           coord_uncertainty, basis_of_record, kingdom, phylum, class, source
## ---------------------------------------------------------------------

query_nbn <- function(species_names = NULL, bbox) {
  ## Console header for the NBN query phase.
  cat("Querying NBN Atlas (REST API)...\n")

  ## Ensure jsonlite is available — it's needed to parse the JSON
  ## responses from the NBN REST API. Install it if missing.
  if (!requireNamespace("jsonlite", quietly = TRUE)) install.packages("jsonlite")

  if (is.null(species_names)) {
    ## NBN Atlas doesn't support practical bulk downloads without a
    ## species filter, so we skip and return empty. The all-taxa richness
    ## map relies on GBIF + OBIS instead.
    cat("  Skipping NBN for bulk query (use GBIF + OBIS)\n")
    return(tibble())
  }

  ## The base URL for NBN's occurrence search endpoint.
  nbn_base <- "https://records-ws.nbnatlas.org/occurrences/search"
  ## Number of records to fetch per page. NBN paginates large result sets;
  ## we loop until we've retrieved all pages.
  page_size <- 1000

  ## Build a Solr spatial filter string from the bounding box.
  ## NBN Atlas uses Solr's range syntax: field:[min TO max]. We
  ## constrain both latitude and longitude to the study area.
  spatial_fq <- sprintf(
    "decimalLatitude:[%f TO %f] AND decimalLongitude:[%f TO %f]",
    bbox["ymin"], bbox["ymax"], bbox["xmin"], bbox["xmax"]
  )

  ## Build a Solr assertion filter to exclude records with known
  ## coordinate problems. The minus sign means "NOT" — we exclude
  ## records that have any of these three assertion flags.
  qa_fq <- paste0(
    "-assertions:COORDINATE_INVALID",
    " AND -assertions:ZERO_COORDINATE",
    " AND -assertions:COUNTRY_COORDINATE_MISMATCH"
  )

  ## -- Internal helper: .nbn_fetch --
  ## Fetches all pages of results for a single species from NBN.
  ## Uses a repeat loop that increments the start index until a
  ## page returns fewer records than page_size (meaning we've
  ## reached the end) or the response is empty/NULL.
  .nbn_fetch <- function(sp) {
    ## Accumulator list — each page's data frame is appended here,
    ## then combined at the end. More efficient than growing a single
    ## data frame row by row.
    all_recs <- list()
    ## Start index for pagination — NBN uses 0-based indexing.
    start <- 0
    repeat {
      ## Construct the full URL with query parameters:
      ##   q=       — search by LSID (taxon key) OR scientific name
      ##   fq=      — spatial filter (bbox) and assertion filter
      ##   pageSize — records per page
      ##   startIndex — offset for pagination
      ##   sort/dir — sort by year descending (newest first)
      ## All dynamic values are URL-encoded to handle spaces and special chars.
      url <- sprintf(
        "%s?q=lsid:%s OR scientificName:\"%s\"&fq=%s&fq=%s&pageSize=%d&startIndex=%d&sort=year&dir=desc",
        nbn_base,
        utils::URLencode(sp, reserved = TRUE),
        utils::URLencode(sp, reserved = TRUE),
        utils::URLencode(spatial_fq, reserved = TRUE),
        utils::URLencode(qa_fq, reserved = TRUE),
        page_size, start
      )

      ## Fetch and parse the JSON response. tryCatch returns NULL on
      ## any error (network timeout, malformed JSON, etc.) so the loop
      ## can exit gracefully rather than crashing.
      resp <- tryCatch(
        jsonlite::fromJSON(url, flatten = TRUE),
        error = function(e) NULL
      )

      ## Exit conditions: no response, no occurrences key in the JSON,
      ## or the occurrences array is empty — we've exhausted all records.
      if (is.null(resp) || is.null(resp$occurrences) ||
          length(resp$occurrences) == 0) break
      recs <- resp$occurrences
      ## Double-check it's actually a non-empty data frame.
      if (!is.data.frame(recs) || nrow(recs) == 0) break

      ## Append this page to the accumulator.
      all_recs[[length(all_recs) + 1]] <- recs

      ## If this page had fewer records than page_size, we've reached
      ## the last page — no need to request more.
      if (nrow(recs) < page_size) break

      ## Move the start index forward by one page.
      start <- start + page_size
    }
    ## If no pages returned any records, return NULL (not an empty tibble)
    ## so the caller can distinguish "no data" from "empty after filtering".
    if (length(all_recs) == 0) return(NULL)
    ## Combine all pages into a single data frame.
    bind_rows(all_recs)
  }

  ## ---- Main per-species loop ----
  results <- list()
  for (i in seq_along(species_names)) {
    sp <- species_names[i]
    ## Progress indicator: species number and name.
    cat(sprintf("  [%d/%d] %s", i, length(species_names), sp))

    recs <- NULL
    success <- FALSE
    ## Retry loop: try up to 2 times. If the first attempt gets rate-
    ## limited (HTTP 412 or 429), wait 10 seconds and try once more.
    for (attempt in 1:2) {
      tryCatch({
        ## Call the pagination helper to fetch all pages for this species.
        recs <- .nbn_fetch(sp)
        success <- TRUE
      }, error = function(e) {
        ## If the error message mentions rate limiting and this is the
        ## first attempt, wait and retry. Otherwise log and move on.
        if (grepl("412|429|Too Many", e$message) && attempt == 1) {
          cat(" — rate limited, retrying in 10s...")
          Sys.sleep(10)
        } else {
          cat(sprintf(" — error: %s\n", e$message))
        }
      })
      ## If the fetch succeeded, no need for another attempt.
      if (success) break
    }

    if (success && !is.null(recs) && nrow(recs) > 0) {
      ## We have records — now standardise them to the common schema.
      n_raw <- nrow(recs)
      ## Get the column names once for all the defensive checks below.
      cn <- names(recs)

      ## ---- Extract coordinates ----
      ## NBN can return coordinates under different column names
      ## depending on the record source. We check for both variants.
      lon_col <- if ("decimalLongitude" %in% cn) "decimalLongitude" else
                 if ("longitude" %in% cn) "longitude" else NULL
      lat_col <- if ("decimalLatitude" %in% cn) "decimalLatitude" else
                 if ("latitude" %in% cn) "latitude" else NULL

      ## Force coordinates to numeric — NBN sometimes returns them as
      ## character strings, which would break spatial operations.
      if (!is.null(lon_col)) {
        recs[[lon_col]] <- as.numeric(recs[[lon_col]])
      }
      if (!is.null(lat_col)) {
        recs[[lat_col]] <- as.numeric(recs[[lat_col]])
      }

      ## ---- Extract year ----
      ## Prefer the explicit 'year' column; fall back to parsing the
      ## first 4 characters of 'eventDate' (e.g. "2019-06-15" → 2019);
      ## if neither exists, set to NA so temporal filters don't crash.
      if ("year" %in% cn) {
        recs$yr <- as.integer(recs$year)
      } else if ("eventDate" %in% cn) {
        recs$yr <- as.integer(substr(as.character(recs$eventDate), 1, 4))
      } else {
        recs$yr <- NA_integer_
      }

      ## ---- Temporal filter ----
      ## Remove records from before 1970 or with no date, using the
      ## shared qa_filter() function so they are logged for auditing.
      recs <- qa_filter(recs,
        !is.na(recs$yr) & recs$yr >= 1970,
        "year_pre1970_or_missing", "NBN")

      ## ---- Extract species name ----
      ## NBN can use 'scientificName', 'species', or 'taxonConceptID'.
      ## We check in order of preference and fall back to the queried
      ## name if none are present.
      name_col <- if ("scientificName" %in% cn) "scientificName" else
                  if ("species" %in% cn) "species" else
                  if ("taxonConceptID" %in% cn) "taxonConceptID" else NULL

      ## ---- Build standardised output ----
      ## Assemble the 10-column common schema tibble. Every field
      ## uses conditional extraction to handle missing columns
      ## gracefully — if a column doesn't exist, NA is used instead.
      results[[sp]] <- tibble(
        latin_name       = if (!is.null(name_col)) recs[[name_col]] else sp,
        longitude        = if (!is.null(lon_col)) recs[[lon_col]] else NA_real_,
        latitude         = if (!is.null(lat_col)) recs[[lat_col]] else NA_real_,
        year             = recs$yr,
        coord_uncertainty = if ("coordinateUncertaintyInMeters" %in% cn)
          as.numeric(recs$coordinateUncertaintyInMeters) else NA_real_,
        basis_of_record  = if ("basisOfRecord" %in% cn)
          recs$basisOfRecord else NA_character_,
        kingdom          = if ("kingdom" %in% cn) recs$kingdom else NA_character_,
        phylum           = if ("phylum" %in% cn) recs$phylum else NA_character_,
        ## NBN sometimes spells it 'classs' (three s's) — check both.
        class            = if ("classs" %in% cn) recs$classs else
                           if ("class" %in% cn) recs$class else NA_character_,
        source           = "NBN"
      )
      cat(sprintf(" — %d/%d records after QA\n", nrow(results[[sp]]), n_raw))
    } else if (success) {
      cat(" — 0 records\n")
    }

    ## Rate-limiting pause: wait 2 seconds between species to avoid
    ## triggering NBN Atlas's rate limits (HTTP 412/429).
    Sys.sleep(2)
  }
  ## Combine all species results into a single data frame.
  bind_rows(results)
}


## ---------------------------------------------------------------------
##  5. cc_clean
## ---------------------------------------------------------------------
##  SUMMARY: Runs the CoordinateCleaner package's spatial QA tests on
##  a set of occurrence records. CoordinateCleaner catches a specific
##  class of errors that database-level QC often misses: records whose
##  coordinates are technically valid but are actually geocoding artefacts
##  (e.g. coordinates placed at a country capital, country centroid,
##  GBIF headquarters, or the origin at 0,0). Six tests are run; sea
##  and outlier tests are deliberately skipped because the study area
##  spans land and sea, and many species have small local sample sizes.
##  Failed records are logged per-test to the QA log.
##
##  Parameters:
##    df — data frame with longitude, latitude, and latin_name columns
##
##  Returns: the subset of df that passed all six CoordinateCleaner tests
## ---------------------------------------------------------------------

cc_clean <- function(df) {
  ## Short-circuit if the input is empty — nothing to clean.
  if (nrow(df) == 0) return(df)
  cat("  Running CoordinateCleaner...\n")

  ## CoordinateCleaner expects specific column names: 'decimalLongitude',
  ## 'decimalLatitude', and 'species'. Our common schema uses 'longitude',
  ## 'latitude', and 'latin_name', so we create the expected aliases.
  ## We keep the original columns too (mutate, not rename) because the
  ## output needs to stay in our schema.
  cc_input <- df %>%
    mutate(decimalLongitude = longitude, decimalLatitude = latitude,
           species = latin_name)

  ## Run the six coordinate tests. value = "spatialvalid" means the
  ## result includes TRUE/FALSE columns for each individual test plus
  ## a .summary column that is TRUE only if the record passed ALL tests.
  ## verbose = FALSE suppresses per-test console output (we print our
  ## own summary below).
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

  ## Identify which records failed at least one test — these will be
  ## removed from the output.
  flagged <- !cc_result$.summary
  n_flagged <- sum(flagged)

  if (n_flagged > 0) {
    ## Log the failures broken down by individual test, not just the
    ## aggregate. This lets us inspect the QA CSV and see whether most
    ## exclusions are e.g. centroid hits (common, often harmless) vs
    ## GBIF-HQ hits (rare, clearly bogus), which helps judge whether
    ## thresholds need adjusting.
    test_cols <- c(".cap", ".cen", ".equ", ".gbf", ".inst", ".zer")
    test_names <- c("near_capital", "near_centroid", "equal_lat_lon",
                     "near_gbif_hq", "near_institution", "near_zero")
    ## Loop through each test column and log the records that failed it.
    for (j in seq_along(test_cols)) {
      ## Check the column exists — some versions of CoordinateCleaner
      ## may name them slightly differently.
      if (test_cols[j] %in% names(cc_result)) {
        ## FALSE means the record failed this particular test.
        test_fail <- !cc_result[[test_cols[j]]]
        if (sum(test_fail) > 0) {
          ## Append to the global QA log with a specific reason code
          ## like "cc_near_capital" or "cc_near_institution".
          qa_log[[length(qa_log) + 1]] <<- df[test_fail, ] %>%
            mutate(qa_reason = paste0("cc_", test_names[j]),
                   qa_source = "CoordinateCleaner")
        }
      }
    }
    cat(sprintf("    CoordinateCleaner: %d records flagged\n", n_flagged))
  }

  ## Return only the records that passed ALL tests.
  df[!flagged, ]
}


## ---------------------------------------------------------------------
##  6. deduplicate_records
## ---------------------------------------------------------------------
##  SUMMARY: Removes duplicate occurrence records that appear across
##  multiple databases. The same observation often exists in both GBIF
##  and NBN Atlas (because NBN feeds into GBIF). Duplicates are
##  identified by matching on species name + coordinates rounded to
##  3 decimal places (~100 m) + year. This is slightly conservative
##  (may merge genuinely separate records that are very close) but
##  preferable to inflating counts with cross-database duplicates.
##
##  Parameters:
##    df — data frame with latin_name, longitude, latitude, year columns
##
##  Returns: deduplicated data frame (rounded helper columns removed)
## ---------------------------------------------------------------------

deduplicate_records <- function(df) {
  ## Capture the count before deduplication for the progress message.
  n_before <- nrow(df)

  df <- df %>%
    ## Round coordinates to 3 decimal places. At UK latitudes, 0.001
    ## degrees is approximately 100 metres — so records within ~100 m
    ## of each other are treated as the same location.
    mutate(
      lon_r = round(longitude, 3),
      lat_r = round(latitude, 3)
    ) %>%
    ## Keep the first occurrence of each species+location+year
    ## combination. .keep_all = TRUE retains all columns from the
    ## first match (not just the grouping columns).
    distinct(latin_name, lon_r, lat_r, year, .keep_all = TRUE) %>%
    ## Remove the temporary rounded-coordinate columns — they were
    ## only needed for the deduplication matching and shouldn't
    ## persist in the output.
    select(-lon_r, -lat_r)

  ## Report how many duplicates were found and removed.
  cat(sprintf("Deduplication: %d → %d records (%d duplicates removed)\n",
      n_before, nrow(df), n_before - nrow(df)))
  df
}


## ---------------------------------------------------------------------
##  7. write_qa_log
## ---------------------------------------------------------------------
##  SUMMARY: Writes the accumulated QA exclusion log to a CSV file for
##  auditing. Called once at the end of the script after all maps have
##  been produced. The CSV contains every record that was excluded during
##  processing, along with the reason code and source database. Also
##  prints a summary table showing how many records were excluded by
##  each QA step, so the user can quickly see which filters have the
##  biggest impact.
##
##  Parameters: none (reads from the global qa_log list)
##
##  Returns: invisible NULL
## ---------------------------------------------------------------------

write_qa_log <- function() {
  ## If no records were excluded during the entire run (unlikely but
  ## possible with a tiny test dataset), just say so and exit.
  if (length(qa_log) == 0) {
    cat("QA log: no records excluded\n")
    return(invisible(NULL))
  }

  ## Combine all exclusion batches into a single data frame.
  ## Each element of qa_log is a data frame from a different qa_filter()
  ## or cc_clean() call. bind_rows handles varying column sets.
  qa_all <- bind_rows(qa_log) %>%
    ## Select only the columns we want in the output CSV. any_of()
    ## means "include these if they exist, silently skip if not" —
    ## this handles edge cases where some batches lack certain columns.
    select(any_of(c("latin_name", "longitude", "latitude", "year",
                     "source", "qa_reason", "qa_source")))

  ## Write the exclusion log to CSV. This file is the definitive audit
  ## trail — if someone questions why a species appears absent, this
  ## is where to check if records existed but were filtered out.
  write_csv(qa_all, "data/processed/species/qa_excluded.csv")
  cat(sprintf("QA log: %d excluded records written to data/processed/species/qa_excluded.csv\n",
      nrow(qa_all)))

  ## Print a cross-tabulation of exclusion counts by source and reason.
  ## This gives a quick overview of where records are being lost:
  ## e.g. "GBIF_bulk / coord_uncertainty: 500" tells you that 500
  ## GBIF records had uncertainty > 10 km.
  qa_summary <- qa_all %>% count(qa_source, qa_reason, name = "n_excluded")
  cat("\nQA exclusion summary:\n")
  ## n = Inf ensures the full table is printed even if it exceeds
  ## tibble's default row limit.
  print(qa_summary, n = Inf)
}


## ---------------------------------------------------------------------
##  8. compute_grid_metrics
## ---------------------------------------------------------------------
##  SUMMARY: Computes seven per-cell metrics for a single realm (marine
##  or terrestrial) and joins them to the spatial grid. This is the
##  analytical core of the species richness map. The metrics are:
##    1. Species richness (count of distinct species per cell)
##    2. Recording effort (total records, GBIF/OBIS split)
##    3. Temporal recency — median recording year per cell
##    4. Temporal recency — year of most recent record per cell
##    5. Beta diversity (mean Sorensen dissimilarity to queen-contiguous
##       neighbour cells)
##    6. Community composition similarity (mean Sorensen similarity to
##       ALL other cells, with full matrix saved as RDS)
##    7. Per-phylum richness breakdown (one column per phylum)
##
##  Parameters:
##    df          — data frame of gridded records with realm, cell_id,
##                  latin_name, year, source, phylum columns
##    realm_label — "marine" or "terrestrial"
##
##  Returns: sf object (grid polygons with metric columns joined),
##           or NULL if no records for this realm
## ---------------------------------------------------------------------

compute_grid_metrics <- function(df, realm_label) {
  ## Filter the full dataset to just this realm. Also exclude rows
  ## with no cell assignment (outside the grid), no species name, or
  ## empty species name strings — any of these would produce NAs or
  ## errors in the summarise steps below.
  df_realm <- df %>%
    filter(realm == realm_label, !is.na(cell_id),
           !is.na(latin_name), latin_name != "")

  ## If the realm has no records (e.g. no marine records found), return
  ## NULL. The calling code checks for NULL before proceeding to plots.
  if (nrow(df_realm) == 0) {
    cat(sprintf("  No %s records — skipping\n", realm_label))
    return(NULL)
  }

  ## ---- Metric 1: Species richness ----
  ## Count the number of distinct species in each grid cell.
  ## n_distinct() counts unique values, so duplicate records of the
  ## same species in the same cell are counted only once.
  richness <- df_realm %>%
    group_by(cell_id) %>%
    summarise(species_richness = n_distinct(latin_name), .groups = "drop")

  ## ---- Metric 2: Recording effort ----
  ## Count total records and the breakdown by database source. This
  ## helps identify cells where richness might be inflated or deflated
  ## by uneven sampling effort.
  effort <- df_realm %>%
    group_by(cell_id) %>%
    summarise(
      record_count = n(),
      gbif_records = sum(source == "GBIF"),
      obis_records = sum(source == "OBIS"),
      .groups = "drop"
    )

  ## ---- Metrics 3 & 4: Temporal recency ----
  ## Two complementary measures of how recently each cell was surveyed:
  ##   - median_year: the "typical" recording date (robust to outliers)
  ##   - most_recent_year: when the cell was last visited
  ## Together these distinguish cells with steady recording history
  ## from those with only old records or a single recent visit.
  temporal <- df_realm %>%
    group_by(cell_id) %>%
    summarise(
      median_year      = median(year, na.rm = TRUE),
      most_recent_year = max(year, na.rm = TRUE),
      .groups = "drop"
    )

  ## ---- Build species-by-cell matrix (shared by beta & similarity) ----
  ## First, count records per species per cell. Then pivot to a wide
  ## matrix where rows are cells and columns are species. This matrix
  ## is used for both beta diversity and community similarity.
  sp_cell <- df_realm %>%
    count(cell_id, latin_name, name = "abundance")
  sp_wide <- sp_cell %>%
    pivot_wider(id_cols = cell_id, names_from = latin_name,
                values_from = abundance, values_fill = 0L)
  ## Store the cell IDs separately — they'll be our row labels for
  ## the numeric matrix operations below.
  cell_ids_wide <- sp_wide$cell_id
  ## Convert to a plain numeric matrix (dropping the cell_id column).
  ## vegan functions need a matrix, not a tibble.
  sp_mat <- as.matrix(sp_wide[, -1])
  ## Convert abundance to presence/absence (1/0). Sorensen dissimilarity
  ## operates on binary data — we don't want abundance to influence it.
  pa_mat <- (sp_mat > 0) * 1L

  ## ---- Metric 5: Beta diversity ----
  ## For each cell, compute the mean Sorensen dissimilarity to its
  ## queen-contiguous neighbours (8 surrounding cells in a regular grid).
  ## High beta diversity means the cell's species composition differs
  ## strongly from its neighbours — it's an ecological boundary or hotspot.

  ## Create a lookup from cell_id to row index in the matrix, so we
  ## can quickly find which matrix rows correspond to a cell's neighbours.
  row_lookup <- setNames(seq_along(cell_ids_wide), cell_ids_wide)
  beta_vals <- rep(NA_real_, nrow(pa_mat))

  for (i in seq_len(nrow(pa_mat))) {
    ## Get this cell's ID and find its neighbour cell IDs from the
    ## pre-computed grid_touches adjacency list (built earlier in the
    ## script using st_touches on the grid polygons).
    cid <- cell_ids_wide[i]
    nb_cids <- grid_touches[[cid]]
    ## Look up which rows in the matrix correspond to those neighbours.
    ## na.omit handles neighbours that have no species records (and
    ## therefore aren't in the matrix).
    nb_rows <- na.omit(row_lookup[as.character(nb_cids)])
    ## Need at least one neighbour with data to compute dissimilarity.
    if (length(nb_rows) == 0) next

    ## Build a submatrix: focal cell as row 1, neighbours as remaining
    ## rows. vegdist computes pairwise distances between all rows.
    pair_mat <- pa_mat[c(i, nb_rows), , drop = FALSE]
    if (nrow(pair_mat) < 2) next

    ## Compute Sorensen dissimilarity (method = "bray" with binary = TRUE
    ## is equivalent to Sorensen). The distance matrix is symmetric, so
    ## row 1 vs all others gives us focal-to-each-neighbour distances.
    d <- vegan::vegdist(pair_mat, method = "bray", binary = TRUE)
    ## Take the mean distance from the focal cell to all its neighbours.
    beta_vals[i] <- mean(as.matrix(d)[1, -1])
  }
  beta_div <- tibble(cell_id = cell_ids_wide, beta_diversity = beta_vals)

  ## ---- Metric 6: Community composition similarity ----
  ## Compute the mean Sorensen SIMILARITY (1 - dissimilarity) between
  ## each cell and ALL other cells, not just neighbours. This gives a
  ## landscape-level view: cells with high similarity share species
  ## composition with the broader study area, while low-similarity
  ## cells are compositionally distinct.
  if (nrow(pa_mat) >= 2) {
    ## Compute the full pairwise Sorensen distance matrix.
    sorensen_dist <- vegan::vegdist(pa_mat, method = "bray", binary = TRUE)
    ## Convert distance to similarity (1 - distance).
    sorensen_sim  <- 1 - as.matrix(sorensen_dist)
    ## Set the diagonal to NA — a cell's similarity to itself is
    ## always 1, which would inflate the mean.
    diag(sorensen_sim) <- NA
    ## Mean similarity of each cell to all other cells.
    mean_similarity <- rowMeans(sorensen_sim, na.rm = TRUE)

    ## Save the full similarity matrix as an RDS file. This is a large
    ## object (n_cells x n_cells) that downstream analyses might use
    ## for clustering, ordination, or other community analyses.
    saveRDS(sorensen_sim,
            sprintf("outputs/species/community_similarity_matrix_%s.rds",
                    realm_label))
    cat(sprintf("  Saved %s Sorensen similarity matrix (%d x %d)\n",
        realm_label, nrow(sorensen_sim), ncol(sorensen_sim)))
  } else {
    ## With fewer than 2 cells, similarity is undefined.
    mean_similarity <- NA_real_
  }
  community_sim <- tibble(
    cell_id = cell_ids_wide,
    community_similarity = mean_similarity
  )

  ## ---- Metric 7: Per-phylum richness ----
  ## Count distinct species per phylum per cell, then pivot so each
  ## phylum becomes its own column (e.g. richness_Chordata,
  ## richness_Arthropoda). This feeds the taxonomic breakdown maps.
  taxon_richness <- df_realm %>%
    filter(!is.na(phylum), phylum != "") %>%
    group_by(cell_id, phylum) %>%
    summarise(spp = n_distinct(latin_name), .groups = "drop") %>%
    pivot_wider(id_cols = cell_id, names_from = phylum,
                values_from = spp, values_fill = 0L,
                names_prefix = "richness_")

  ## ---- Assemble all metrics ----
  ## Join the seven metric tables into a single data frame using
  ## cell_id as the key. left_join ensures cells present in the
  ## richness table (which has all cells with any records) get NA
  ## for metrics where they have no data.
  metrics <- richness %>%
    left_join(effort, by = "cell_id") %>%
    left_join(temporal, by = "cell_id") %>%
    left_join(beta_div, by = "cell_id") %>%
    left_join(community_sim, by = "cell_id") %>%
    left_join(taxon_richness, by = "cell_id")

  ## ---- Join to grid polygons ----
  ## Attach the metrics to the spatial grid (grid_osgb) so we can
  ## map them. Filter out cells with no species richness — these are
  ## grid squares that fall within the study area but had no records.
  grid_sf <- grid_osgb %>%
    left_join(metrics, by = "cell_id") %>%
    filter(!is.na(species_richness))

  grid_sf
}


## ---------------------------------------------------------------------
##  9. make_panel
## ---------------------------------------------------------------------
##  SUMMARY: Creates a single ggplot2 map panel for one realm (marine
##  or terrestrial). Each panel shows a gridded metric (species richness,
##  beta diversity, etc.) overlaid on the study area with coastline,
##  boundary polygon, scale bar, and north arrow. Marine panels clip to
##  the Solway Marine Region (SMR) boundary with a blue sea background;
##  terrestrial panels clip to the Dumfries & Galloway (D&G) local
##  authority boundary with a grey land background. Supports two scale
##  types: viridis (sequential) and diverging (gradient2 with a midpoint).
##
##  Parameters:
##    grid_sf     — sf object with the metric to plot
##    fill_var    — character: column name in grid_sf to map to fill colour
##    fill_label  — character: label for the colour scale legend
##    title       — character: plot title
##    subtitle    — character: plot subtitle
##    realm       — "marine" or "terrestrial" (determines clipping boundary)
##    scale_type  — "viridis" (default) or "diverging"
##    viridis_opt — viridis palette name (default "viridis")
##    direction   — viridis direction: -1 (default, reversed) or 1
##    midpoint    — numeric: centre of the diverging scale (required if diverging)
##    low, high   — optional colours for diverging scale endpoints
##
##  Returns: a ggplot object (not yet saved or printed)
## ---------------------------------------------------------------------

make_panel <- function(grid_sf, fill_var, fill_label, title, subtitle,
                       realm, scale_type = "viridis",
                       viridis_opt = "viridis", direction = -1,
                       midpoint = NULL, low = NULL, high = NULL) {
  ## Transform the grid data to British National Grid (EPSG:27700) for
  ## plotting. All map elements need to be in the same CRS — we use
  ## BNG because it produces correctly proportioned maps of the UK.
  grid_bng_plot <- st_transform(grid_sf, 27700)

  ## Choose the clipping boundary and background colours based on realm.
  ## Marine maps show the sea region with a blue background; terrestrial
  ## maps show the local authority area with a grey background. These
  ## objects (smr_solway_bng, las_dg_bng) are defined earlier in the
  ## main script.
  clip_boundary <- if (realm == "marine") smr_solway_bng else las_dg_bng
  bg_fill  <- if (realm == "marine") "lightblue" else "grey90"
  bg_color <- if (realm == "marine") "steelblue" else "grey50"

  ## Build the base map with four layers:
  p <- ggplot() +
    ## Layer 1: Background boundary polygon — provides geographic context
    ## and fills the non-grid area with the appropriate colour.
    geom_sf(data = clip_boundary, fill = bg_fill, colour = bg_color,
            linewidth = 0.3) +
    ## Layer 2: Grid cells coloured by the metric. colour = NA removes
    ## cell borders so the colour fill reads clearly.
    geom_sf(data = grid_bng_plot, aes(fill = .data[[fill_var]]),
            colour = NA) +
    ## Layer 3: Coastline for geographic reference. Drawn on top of the
    ## grid so it's visible.
    geom_sf(data = scamp_coastline, colour = "grey30", linewidth = 0.3) +
    ## Layer 4: SCAMP AOI boundary — the overall study area outline.
    ## fill = NA makes it transparent; the black border marks the boundary.
    geom_sf(data = scamp_aoi_bng, fill = NA, colour = "black",
            linewidth = 0.5)

  ## Choose the colour scale based on scale_type.
  if (scale_type == "diverging" && !is.null(midpoint)) {
    ## Diverging scale: used for metrics where a midpoint has meaning
    ## (e.g. median year centred on 2010 — blue for older, red for newer).
    ## scale_fill_gradient2 maps values below/above the midpoint to
    ## different colour ramps.
    p <- p + scale_fill_gradient2(
      name = fill_label, midpoint = midpoint,
      low = if (is.null(low)) "#2166ac" else low,
      mid = "white",
      high = if (is.null(high)) "#b2182b" else high,
      na.value = "grey80"
    )
  } else {
    ## Sequential scale: used for most metrics (richness, effort, etc.).
    ## Viridis palettes are perceptually uniform and colourblind-safe.
    ## direction = -1 reverses the palette (high values are darker).
    p <- p + scale_fill_viridis_c(option = viridis_opt, name = fill_label,
                                   direction = direction, na.value = "grey80")
  }

  ## Add map furniture and apply the shared theme.
  p +
    ## Scale bar in the bottom-left corner for geographic reference.
    annotation_scale(location = "bl", width_hint = 0.25,
                     style = "ticks", line_width = 0.5) +
    ## North arrow in the top-right corner.
    annotation_north_arrow(location = "tr", which_north = "true",
                           height = unit(1, "cm"), width = unit(1, "cm"),
                           style = north_arrow_minimal) +
    ## Title and subtitle for the panel.
    labs(title = title, subtitle = subtitle) +
    ## Apply the shared minimal theme defined earlier in the script.
    map_theme +
    ## Force the coordinate system to BNG — ensures correct aspect ratio.
    coord_sf(crs = 27700)
}


## ---------------------------------------------------------------------
##  10. save_paired
## ---------------------------------------------------------------------
##  SUMMARY: Places a marine and terrestrial map panel side by side
##  using cowplot::plot_grid, then saves the combined figure as both
##  PDF (for publication) and PNG (for presentations/reports) at 300 dpi.
##  Also prints the combined plot to the RStudio viewer so the user
##  can inspect it immediately.
##
##  Parameters:
##    p_marine       — ggplot object for the marine panel
##    p_terrestrial  — ggplot object for the terrestrial panel
##    filename       — base name (without extension) for the output files
##    w              — figure width in inches (default 16)
##    h              — figure height in inches (default 10)
##
##  Returns: invisible (called for side effects: file saves + print)
## ---------------------------------------------------------------------

save_paired <- function(p_marine, p_terrestrial, filename, w = 16, h = 10) {
  ## Arrange the two panels side by side in a 1-row, 2-column layout.
  ## The labels appear above each panel to identify the realm.
  p_combined <- plot_grid(p_marine, p_terrestrial, ncol = 2,
                          labels = c("Marine", "Terrestrial"),
                          label_size = 11)

  ## Save as PDF — vector format, ideal for publication and printing.
  ## The figures/ directory is created by the project setup script.
  ggsave(paste0("figures/", filename, ".pdf"), p_combined,
         width = w, height = h, units = "in")

  ## Save as PNG — raster format at 300 dpi, suitable for reports,
  ## presentations, and web display.
  ggsave(paste0("figures/", filename, ".png"), p_combined,
         width = w, height = h, units = "in", dpi = 300)

  ## Confirm the save in the console log.
  cat(sprintf("Saved: figures/%s.pdf & .png\n", filename))

  ## Print to the RStudio plot viewer so the user can see the result
  ## immediately without opening the saved file.
  print(p_combined)
}


## ---------------------------------------------------------------------
##  11. make_taxon_panels
## ---------------------------------------------------------------------
##  SUMMARY: Creates a multi-panel figure showing species richness
##  broken down by taxonomic group (phylum). Identifies the top 6 phyla
##  by record count in the given realm, creates a map panel for each
##  using make_panel(), arranges them in a grid layout (3 columns),
##  and saves the result as PDF and PNG. This reveals which taxonomic
##  groups drive the overall richness pattern — e.g. whether a cell's
##  high richness comes mainly from plants, insects, birds, etc.
##
##  Parameters:
##    grid_sf      — sf object with richness_<phylum> columns from
##                   compute_grid_metrics()
##    df_records   — the full gridded records data frame (used to
##                   determine which phyla are most common)
##    realm_label  — "marine" or "terrestrial"
##
##  Returns: invisible NULL (called for side effects: file saves + print)
## ---------------------------------------------------------------------

make_taxon_panels <- function(grid_sf, df_records, realm_label) {
  ## Filter the records to this realm and exclude records with no phylum
  ## information — we can't break down by taxon if the taxon is unknown.
  realm_records <- df_records %>%
    filter(realm == realm_label, !is.na(phylum), phylum != "")

  ## Identify the 6 most common phyla by total record count. Using
  ## record count (not species count) as the ranking criterion means
  ## phyla with many records get mapped first — these are the groups
  ## most likely to show meaningful spatial patterns.
  top_phyla <- realm_records %>%
    count(phylum, sort = TRUE) %>%
    head(6) %>%
    pull(phylum)

  ## Map phylum names to the column names in the grid_sf data frame.
  ## compute_grid_metrics() prefixed them with "richness_" when
  ## creating the per-phylum columns.
  richness_cols <- paste0("richness_", top_phyla)
  ## Only keep columns that actually exist in the grid — a phylum might
  ## have records but zero grid cells with data if all records fell
  ## outside the grid.
  present_cols <- richness_cols[richness_cols %in% names(grid_sf)]

  ## If no phyla have grid data, there's nothing to plot.
  if (length(present_cols) == 0) return(invisible(NULL))

  ## Create one map panel per phylum using make_panel(). Each panel
  ## shows the spatial distribution of richness for that phylum alone.
  panels <- lapply(present_cols, function(col) {
    ## Strip the "richness_" prefix to get the phylum name for the title.
    phylum_name <- sub("^richness_", "", col)
    ## Reuse make_panel() with the phylum-specific richness column.
    make_panel(grid_sf, col, "Species",
               phylum_name, paste(realm_label, "richness"),
               realm = realm_label, viridis_opt = "viridis")
  })

  ## Arrange all panels in a grid with up to 3 columns. If there are
  ## more than 3 panels, they wrap to additional rows.
  n_panels <- length(panels)
  ncols <- min(3L, n_panels)
  p_grid <- plot_grid(plotlist = panels, ncol = ncols)

  ## Scale the figure height based on how many rows of panels there are.
  ## 3 or fewer panels = single row = 8 inches; more = two rows = 16 inches.
  fname <- sprintf("map1g_taxon_breakdown_%s", realm_label)
  h <- if (n_panels <= 3) 8 else 16

  ## Save as PDF and PNG, same as save_paired().
  ggsave(paste0("figures/", fname, ".pdf"), p_grid,
         width = 16, height = h, units = "in")
  ggsave(paste0("figures/", fname, ".png"), p_grid,
         width = 16, height = h, units = "in", dpi = 300)

  ## Console confirmation and preview.
  cat(sprintf("Saved: figures/%s.pdf & .png (%d phyla)\n", fname, n_panels))
  print(p_grid)
}
