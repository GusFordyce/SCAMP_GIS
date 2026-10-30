# ==============================================================================
# handover_metadata.R
# Hand-entered metadata for the SCAMP GIS handover package.
# Sourced by scripts/09_handover_package.R. Edit values here, not in the script.
# Derived automatically by the script (do not enter here): WGS84 bounding box,
# spatial resolution range, CRS, metadata date, lineage.
# Keys match what 09_handover_package.R reads via meta_get().
# ==============================================================================

handover_meta <- list(

  title             = "SCAMP Biodiversity Monitoring Strategy: GIS Data Handover",
  alternative_title = "SCAMP GIS Handover v1",
  version           = "1.0",
  language          = "eng",

  identifier = list(
    code      = "DGC-SCAMP-GIS-HANDOVER-v1",
    codespace = "Dumfries and Galloway Council"
  ),

  abstract = paste(
    "Derived spatial datasets produced for the Solway Coast and Marine Project",
    "(SCAMP) Biodiversity Monitoring Strategy, covering the Solway Scottish Marine",
    "Region and a five-mile coastal land strip within Dumfries and Galloway.",
    "The package includes known habitat extent for seagrass, native oyster,",
    "saltmarsh and sand dune; modelled restoration potential for seagrass, native",
    "oyster and sand dune; coastal woodland as classes of a land cover map; the",
    "physical marine environment; marine pressures and catchment diffuse pollution",
    "risk; biodiversity surfaces derived from public occurrence records; larval",
    "connectivity; restoration, recovery and monitoring priority scores; and",
    "environmental bands for selecting matched control sites.",
    "Vector data are supplied in a single GeoPackage and rasters as Cloud-Optimised",
    "GeoTIFFs, with aligned marine (286 m) and terrestrial (100 m) multi-band",
    "datasets. Intended for monitoring design and strategic planning; biodiversity",
    "layers are not suitable for baselining."
  ),

  # ---- dates and maintenance -------------------------------------------------
  publication_date = "2026-10-02",
  frequency        = "notPlanned",

  # ---- responsible parties ---------------------------------------------------
  # Owner of the resource (GEMINI "responsible organisation").
  responsible_org = list(
    organisation = "Dumfries and Galloway Council",
    email        = "nick.chisholm@dumgal.gov.uk",
    role         = "owner"
  ),
  # Originator. Not currently written by the script; recorded for completeness
  # and stated in README / METHODS.
  originator = list(
    individual   = "Alexander Fordyce",
    organisation = "Independent consultant",
    role         = "originator"
  ),
  metadata_contact = list(
    organisation = "Dumfries and Galloway Council",
    email        = "nick.chisholm@dumgal.gov.uk",
    role         = "pointOfContact"
  ),

  # ---- licence and constraints -----------------------------------------------
  licence     = "Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)",
  licence_url = "https://creativecommons.org/licenses/by-nc/4.0/",
  attribution = paste(
    "Attribute as: Dumfries and Galloway Council (2026) SCAMP GIS Data Handover",
    "v1.0, prepared by Alexander Fordyce."
  ),
  use_constraints = paste(
    "Third-party source data retain their original licences and attribution",
    "requirements, listed in SOURCES_AND_LICENCES.md. Biodiversity layers derive",
    "from opportunistic occurrence records, reflect recorder effort as well as",
    "ecology, and are not suitable for baselining or condition assessment.",
    "Suitability models are screening tools; field verification is required",
    "before site selection."
  ),
  access_limitations = paste(
    "No limitations on public access. Onward distribution at the discretion of",
    "Dumfries and Galloway Council."
  ),

  # ---- classification --------------------------------------------------------
  topic_categories = c("biota", "environment", "oceans", "inlandWaters", "elevation"),

  inspire_themes = c(
    "Habitats and biotopes",
    "Species distribution",
    "Sea regions",
    "Land cover",
    "Elevation"
  ),

  keywords = c(
    "Solway Firth", "Dumfries and Galloway", "SCAMP", "biodiversity monitoring",
    "habitat restoration", "seagrass", "native oyster", "Ostrea edulis",
    "saltmarsh", "sand dunes", "coastal woodland", "diffuse pollution",
    "cumulative pressure", "connectivity", "environmental stratification",
    "control site selection", "Scotland"
  ),

  # ---- occurrence data sources -----------------------------------------------
  # Named list: the names are what the stage 1 provenance check compares
  # against (GBIF, OBIS, NBN). NBN feeds only the indicator and INNS layers and
  # the sand dune restoration score (via habitat-builder records).
  occurrence_sources = list(
    GBIF = list(
      name        = "GBIF occurrence records",
      provider    = "GBIF Secretariat",
      licence     = "CC0 1.0 / CC BY 4.0 / CC BY-NC 4.0 (per record)",
      attribution = "GBIF.org (2026) GBIF occurrence data, accessed via the GBIF API",
      url         = "https://www.gbif.org"
    ),
    OBIS = list(
      name        = "OBIS occurrence records",
      provider    = "OBIS, IOC-UNESCO",
      licence     = "CC0 1.0 / CC BY 4.0 / CC BY-NC 4.0 (per dataset)",
      attribution = "OBIS (2026) Ocean Biodiversity Information System, IOC-UNESCO",
      url         = "https://obis.org"
    ),
    NBN = list(
      name        = "NBN Atlas occurrence records (indicator and INNS species only)",
      provider    = "National Biodiversity Network Trust",
      licence     = "OGL / CC0 / CC BY / CC BY-NC (per dataset)",
      attribution = "NBN Atlas occurrence download, accessed 2026",
      url         = "https://nbnatlas.org"
    )
  ),

  # ---- distribution ----------------------------------------------------------
  resource_locator = NA_character_   # add the transfer link on Friday if wanted
)
