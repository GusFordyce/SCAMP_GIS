##===================================================================
##=====================================================================BOUNDARIES
##===================================================================

require(sf)

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
