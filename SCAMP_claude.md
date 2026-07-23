# CLAUDE.md — SCAMP GIS Analysis

## Project overview

This is the GIS analysis and biodiversity monitoring design for the Solway Coast and Marine Project (SCAMP). The project spans five habitat types (seagrass, native oyster reef, saltmarsh, sand dunes, coastal woodland) across 210 miles of the Solway coastline in southwest Scotland. The analysis harmonises ~30 spatial datasets to produce habitat maps, threat overlays, species inventories, and a spatially balanced monitoring network with BACI design.

Local workspace: `C:/Users/gusford/SCAMP_GIS`

## Language and environment

- R (latest stable), run in RStudio on Windows
- Git version control, GitHub remote

## Core packages — always prefer these

| Purpose | Package | Notes |
|---------|---------|-------|
| Vector spatial | `sf` | Never `sp`. Use `st_read()`, `st_write()`, `st_crop()`, `st_join()` |
| Raster spatial | `terra` | Never `raster`. Memory-efficient. Use `rast()`, `crop()`, `mask()`, `writeRaster()` |
| Data wrangling | `dplyr`, `tidyr` | Pipe with `%>%` (magrittr), not `\|>` |
| String ops | `stringr` | |
| Plotting | `ggplot2` + `tidyterra` | `geom_spatraster()`, `geom_sf()` |
| Map furniture | `ggspatial` | Scale bars, north arrows |
| Basemaps | `rnaturalearth` | Coastlines, boundaries |
| NetCDF | `ncdf4` | For Copernicus and Bio-ORACLE layers |
| Tables | `data.table` | Only when performance requires it |
| Reading files | `readr`, `readxl` | |

## Coding style

### Naming
- **snake_case** for all variables and functions: `fishery_bathy_shallow`, `process_esm()`
- Descriptive names that read as what the object *is*: `ocean_mask`, `coast_aoi`, `habs_aligned`
- Prefix related objects: `fishery_bathy_full`, `fishery_bathy_shallow`, `fishery_bathy_deep`

### Assignment and booleans
- Use `<-` for assignment, never `=`
- `TRUE`/`FALSE` in formal code; `T`/`F` acceptable in quick interactive calls

### Pipes
- Use `%>%` (magrittr pipe), not the native `|>`
- Chain dplyr verbs and simple transformations
- Break long chains with one verb per line, indented

### Script structure
- Section headers using `#===` style:
```r
##===================================================================
##=====================================================================SECTION NAME
```
- Inline comments explaining *why*, not just *what*
- Keep commented-out install lines at top of script for reproducibility
- Use `require()` for loading packages at script top

### Functions
- Write processing functions for any operation applied to multiple files
- Include `cat()` progress messages for batch operations
- Return the result, not just side effects

### File I/O
- Use **relative paths** within the project, never `file.choose()` or `choose.dir()`
- Organise project as: `data/`, `scripts/`, `outputs/`, `figures/`
- Raster output: GeoTIFF with compression:
```r
writeRaster(x, "outputs/file.tif", overwrite = TRUE,
            gdal = c("COMPRESS=LZW", "PREDICTOR=2", "TILED=YES"),
            datatype = "FLT4S")
```
- Vector output: GeoPackage (`.gpkg`) preferred over shapefile

### Plotting
- `ggplot2` with `theme_minimal()` as base
- Save with `ggsave()`, both PDF (vector) and PNG (raster backup)
- Use `geom_sf()` for vector layers, `geom_spatraster()` for rasters
- Always include scale bar and north arrow on maps via `ggspatial`
- Publication-quality: `width = 16, height = 10` for landscape maps

## Memory and performance

- Use `terra` throughout — it is disk-backed and memory-efficient
- For large rasters, process in chunks or use `focal()` with `na.policy`
- Compress output rasters with LZW
- Use `data.table` for large tabular operations when dplyr is too slow
- Clean up large intermediates with `rm()` and `gc()` when memory is tight
- Avoid loading entire NetCDF files when only a subset of layers is needed

## Coordinate reference systems

- **Marine/oceanographic data**: WGS84 `EPSG:4326` (Copernicus, Bio-ORACLE, EMODnet)
- **Scottish terrestrial/local analysis**: OSGB36 `EPSG:27700` (OS grid)
- Always reproject to a common CRS before spatial operations
- Document CRS choices in comments

## SCAMP-specific context

### Five priority habitats
1. **Seagrass** — *Zostera marina* (subtidal), *Z. noltii* (intertidal)
2. **Native oyster** — *Ostrea edulis*, reference site at Loch Ryan
3. **Saltmarsh** — NVC communities, managed realignment sites
4. **Sand dunes** — natterjack toad (*Epidalea calamita*) as key indicator
5. **Coastal woodland** — NWSS baseline, to be confirmed with D&G Woodlands

### Key data sources
- NMPi (GeMS habitat records, PMF species, designations)
- NBN Atlas Scotland (species occurrence records)
- Copernicus Marine (SST, salinity, currents — reanalysis and forecast)
- NatureScot (HabMoS, Dynamic Coast, NWSS, designated sites)
- Bio-ORACLE v3.0 (future climate projections)
- EMODnet (bathymetry, biodiversity occurrence)
- SEPA (water quality — WFD classifications)

### Monitoring design
- BACI (Before-After Control-Impact) framework
- Covariate matching for control sites (MatchIt package in R)
- Spatially balanced sampling (GRTS via spsurvey package)
- 10-year time horizon with Tier 1 (intensive) and Tier 2 (annual) campaigns

## Token and output rules

- **Be concise.** Produce working code, not tutorials.
- **No preamble.** Skip "Sure, here's..." or "Let me help you with..."
- **Comment the intent**, not the syntax. Assume R proficiency.
- **One task per response** unless asked for more.
- **Show the code**, don't describe what you would write.
- When modifying existing scripts, show only the changed section with enough context to locate it.
- If a task is ambiguous, state your assumption and proceed rather than asking.
