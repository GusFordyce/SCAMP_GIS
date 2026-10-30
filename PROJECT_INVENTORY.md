# SCAMP GIS — Data Inventory

Solway Coast & Marine Project, Dumfries & Galloway.
Inventory of input datasets, their uses, and the derived layers intended for a
data cube handover to partners.

**Study area (`scamp_aoi`)** — Solway Scottish Marine Region + an 8,047 m (5 mile)
buffer intersected with the Dumfries & Galloway local authority boundary, giving
marine region plus coastal land strip. Defined in `02_themeD_species.R:74-82`,
saved to `data/processed/boundaries/scamp_aoi.gpkg`.

**Working CRS** — EPSG:27700 (British National Grid) throughout, with exceptions
noted in the tables below (species grids, connectivity hexagons and the AOI boxes
are EPSG:4326 and need reprojecting before they join a BNG cube).

---

## 1. Scripts and how they link

### 1.1 What each script does

**`01_project_setup.R`** (25 lines) — Creates the folder structure. Run once after
cloning. No data dependencies.

**`03_boundaries.R`** (21 lines) — Reads the marine region and local authority
shapefiles and prints their structure. Exploratory only; nothing downstream depends
on it. The boundary work that matters happens at the top of `02_themeD_species.R`.

**`02_themeD_species.R`** (2,538 lines) — *Species and the study area.*
Despite the name, this is the script that **defines the project AOI** and so must
run before everything else. It buffers the Solway marine region by 5 miles,
intersects with Dumfries & Galloway, and writes `scamp_aoi.gpkg`. It then runs the
full species pipeline: downloads occurrence records from GBIF, OBIS and NBN;
applies QA/QC (quality flags, date filters, `CoordinateCleaner` spatial checks,
taxonomic name resolution via WoRMS and the GBIF backbone); deduplicates; and
builds a 5 km grid carrying species richness, beta diversity (mean Sørensen
dissimilarity to queen-contiguous neighbours), community similarity, sampling
effort and recency. It also maps INNS, indicator species and Priority Marine
Features.

**`03_habitats.R`** (4,365 lines) — *Known habitat extent.*
Harmonises the five habitats from multiple overlapping sources: seagrass (Solway
Firth Partnership surveys, WCMC, GeMS), saltmarsh (Saltmarsh Survey of Scotland,
WCMC), native oyster (GeMS, historical records), sand dunes (HabMoS, Dynamic
Coast) and woodland (HLCM). For each it resolves source conflicts, produces a
single harmonised layer, and clusters features into named bay-level bounding boxes
used for the zoomed sub-plots. Section 4 also builds the COP30 DEM for the Solway
and writes it to `outputs/physical/`.

**`04_physical.R`** (2,356 lines) — *Physical environment and suitability models.*
Downloads and summarises CMEMS oceanographic data (temperature, salinity,
currents, chlorophyll, oxygen, suspended matter), clips GEBCO bathymetry,
EUSeaMap energy class, UKSeaMap biotopes and UKCP18 sea level projections to the
AOI. Sections 8 and 9 build the **seagrass potential model** (following Environment
Agency v3 methodology) and the **native oyster potential model** by stacking
depth, substrate, salinity, current, turbidity, temperature, oxygen and
chlorophyll envelopes. Section 10 runs **coastalNet** larval dispersal
connectivity between habitat hexagons.

**`scripts/05_localthreats.R`** (3,681 lines) — *Pressures and diffuse pollution.*
Two halves. The marine half reads the Putuhena CartOcean database and EMODnet
layers to build 14 normalised 1 km pressure surfaces (shipping, fishing, noise,
infrastructure, dark vessel activity, INNS) plus a unified vector threats layer.
The catchment half implements **SCIMAP-style diffuse pollution risk**: EUNIS land
cover is crosswalked to source risk weights, multiplied by TWI hydrological
connectivity, then accumulated downstream along HydroSHEDS flow direction to give
river accumulated risk, coastal pour points and per-catchment summaries at
HydroBASINS level 12. Step 8 converts pour point loads into a marine pollution
exposure surface by sea-path decay. Step 9 combines everything into the cumulative
pressure surface.

**`scripts/07_plots.R`** (4,684 lines) — *Analysis and all final maps.*
Reads from every prior script. Note it is **not purely a plotting script** — it
also generates analysis outputs that other sections and scripts consume: the
rarefied richness and beta diversity interpolated surfaces, protection priority
grids, restoration / recovery / monitoring scores, combined action maps, and an
alternative HydroBASINS sub-catchment delineation. On top of that it produces the
priority layers L0–L5, the per-AOI four-panel spatial planning maps, the Solway
overview maps and the hexagon connectivity maps.

**Supporting files**

- `00_setup.R` — package install list only, not part of the pipeline.
- `00_functionsexplained.R` (1,244 lines) — line-by-line annotated reference for the functions defined in `02_themeD_species.R`. Documentation; explicitly not meant to be sourced.
- `scripts/rebuild_clusters.R` (254 lines) — utility to regenerate the habitat cluster boxes from saved intermediates without re-running all of `03_habitats.R`.

### 1.2 How they link

```
                    02_themeD_species.R
                    (defines the AOI)
                            |
              scamp_aoi.gpkg + species grids
                            |
        +-------------------+---------------------+
        |                   |                     |
        v                   v                     v
  03_habitats.R  <---->  04_physical.R     05_localthreats.R
  harmonised            potential models    pressure layers
  habitat extent   -->  + connectivity      + catchment risk
        |                   |                     |
        +-------------------+---------------------+
                            |
                            v
                      07_plots.R
              (analysis surfaces + all figures)
```

The important linkages:

1. **Everything depends on `scamp_aoi.gpkg`**, written by `02_themeD_species.R:81`.
   Every other script clips to it. It must exist before anything else runs.

2. **`03_habitats.R` and `04_physical.R` are mutually dependent.** `04_physical.R`
   reads the harmonised known-extent layers (`seagrass_polygons_harmonised.gpkg`,
   `oyster_known_extent.gpkg`, `saltmarsh_harmonised.gpkg`) that `03_habitats.R`
   produces, and uses them to train and validate the potential models. But
   `03_habitats.R` sections 1e and 3e then read `seagrass_potential_solway.gpkg`
   and `oyster_potential_solway.gpkg` back from `04_physical.R` to draw the
   combined known-plus-potential figures. The practical run order is:
   `03_habitats.R` (harmonisation sections) → `04_physical.R` (full) →
   `03_habitats.R` (plotting sections). This is the single biggest obstacle to a
   clean re-run and is the first item on the improvement list.

3. **`05_localthreats.R` is largely independent.** It only needs `scamp_aoi.gpkg`
   and `scamp_coastline.gpkg` from the boundary step, plus raw data. It can run in
   parallel with the habitat and physical scripts.

4. **`07_plots.R` consumes everything and also produces.** It reads 30+ derived
   layers from `outputs/` and `data/processed/`, and writes new analysis outputs
   back into `outputs/species/`, `outputs/habitats/` and `outputs/catchments/`.
   Some of those it then reads again later in the same script. It is effectively
   "script 6 and 7 merged".

5. **Disk fallbacks make sections re-runnable.** Most scripts guard loads with
   `if (!exists("x")) x <- st_read(...)`, so an individual section can be run in a
   fresh session provided the upstream intermediate files are on disk. This is
   convenient interactively but hides the true dependency graph — the same object
   may come from memory or from disk depending on session state.

### 1.3 Minimum path to reproduce the headline outputs

To regenerate the cumulative pressure map, the AOI four-panel maps and the
overview maps from a clean checkout:

1. `02_themeD_species.R` — sections 1 (boundaries/AOI) and 5 (richness grids)
2. `03_habitats.R` — habitat harmonisation sections
3. `04_physical.R` — sections 4–10 (bathymetry through connectivity)
4. `scripts/05_localthreats.R` — all steps, in particular Step 4 (revised) which
   writes `rivers_with_accumulated_risk_hydrosheds.gpkg`, currently missing
5. `scripts/07_plots.R` — from the biodiversity surface section onward

---

## 2. Input data

### 2.1 Boundaries

| Dataset | Path (`data/raw/boundaries/`) | Format | Source | Used for |
|---|---|---|---|---|
| Scottish Marine Regions | `administrative_units_scottish_marine_regions/...Polygon.shp` | Shapefile | Marine Scotland | Defines Solway marine region; AOI basis; marine mask |
| Local Authority Boundaries | `Local_Authority_Boundaries_-_Scotland/pub_las.shp` | Shapefile | Scottish Government | D&G land extent; coastal strip |
| Harbour limits | `area_management_harbour_limits/...Polygon.shp` | Shapefile | Marine Scotland | Port exclusion in seagrass suitability |
| SAC / SPA / SSSI / Ramsar / MPA / NNR / LNR | individual `.gpkg` files | GeoPackage | NatureScot / JNCC | Designated site context; protection priority |
| RSPB reserves | `RSPB_Public_Reserves_*.gpkg` | GeoPackage | RSPB | Reference site identification |
| IBAs | `IBAs_UK_*.gpkg` | GeoPackage | BirdLife | Designated site context |
| Biosphere zone | `Biosphere_SCAMP_Zone.shp` | Shapefile | Galloway & Southern Ayrshire Biosphere | Context |
| 5-mile zone | `DG_McNabb_5milezone.shp` | Shapefile | Project-supplied | Context (no CRS defined in source) |

### 2.2 Species

| Dataset | Path | Format | Source | Used for |
|---|---|---|---|---|
| GBIF / OBIS / NBN occurrences | downloaded via API in `02_themeD_species.R` | API → CSV/GPKG | GBIF, OBIS, NBN Atlas | Richness, beta diversity, INNS, indicators |
| PMF species points | `data/raw/species/GeMS___Scottish_Priority_Marine_Features_-species points.csv` | CSV | NatureScot GeMS | Priority Marine Feature mapping |
| PMF species polygons | `data/raw/species/...species_polygons.gpkg` | GeoPackage | NatureScot GeMS | Priority Marine Feature mapping |
| GeMS full database | `data/raw/GEMS-PMF_SCOTLAND_GPKG_4326/gems.gpkg` | GeoPackage | NatureScot GeMS | Seagrass/oyster known extent |
| JNCC INNS list | `data/raw/species/jncc_b6_inns_list.csv` | CSV | JNCC | INNS species validation |

### 2.3 Habitats

| Dataset | Path (`data/raw/habitats/`) | Format | Source | Used for |
|---|---|---|---|---|
| HLCM 2022 EUNIS L2 | `General/HLCM-2022-L2_SCOTLAND_TIFF_27700/HLCM_2022_EUNIS_LEVEL2.tif` | GeoTIFF, 20 m | Space Intelligence | Land cover; SCIMAP runoff weights; woodland |
| Habitat Map of Scotland (HabMoS) | `General/Habitat_Map_of_Scotland_*.gpkg` | GeoPackage | NatureScot | Sand dunes, saltmarsh cross-check |
| SFP Zostera polygons | `Seagrass/SFP Zostera polygons/` | Shapefiles | Solway Firth Partnership | Seagrass known extent (primary) |
| WCMC seagrass v7.1 | `Seagrass/WCMC013-014_SeagrassPtPy2021_v7_1/` | Shapefile | UNEP-WCMC | Seagrass known extent (supplementary) |
| Seagrass Potential | `Seagrass/Seagrass_Potential.gpkg` | GeoPackage | Environment Agency | Seagrass suitability comparison |
| Saltmarsh Survey of Scotland | `Saltmarsh/Saltmarsh_Survey_of_Scotland_*.gpkg` | GeoPackage | NatureScot | Saltmarsh known extent |
| WCMC saltmarsh v6.1 | `Saltmarsh/WCMC027_Saltmarsh_v6_1/` | Shapefile | UNEP-WCMC | Saltmarsh supplementary |
| Historical oyster locations | `Oyster/historical_oyster_distribution/01_Oyster_Habitat_Locations.shp` | Shapefile | Historical records | Oyster restoration targeting |
| Dynamic Coast 2 (MHWS 1890/1970/Modern, defences, RCP8.5 2050) | `Sand dunes/DYNAMICCOAST2_SCOTLAND_SHP_27700/shp/` | Shapefiles | Dynamic Coast / Scottish Government | Sand dune erosion and accommodation space |
| Intertidal Substrate Foreshore | `Sand dunes/Intertidal+Substrate+Foreshore.../DEFR00000009.shp` | Shapefile | Defra / Scottish Government | Sandy foreshore identification |

### 2.4 Physical and oceanographic

| Dataset | Path (`data/raw/physical/` unless noted) | Format | Source | Used for |
|---|---|---|---|---|
| CMEMS Global Physics (temp, salinity, currents, bottom T) | `data/processed/physical/{thetao,so,uo,vo,bottomT}/` | NetCDF, ~0.083° monthly | Copernicus Marine | Habitat suitability envelopes |
| CMEMS IBI daily salinity 2013–2023 | `data/processed/physical/so_daily_ibi/` | NetCDF, 0.028° daily | Copernicus Marine | Low-salinity exposure for oyster |
| CMEMS NW Shelf biogeochemistry (chl, O₂) | `data/processed/physical/{chl,o2}/` | NetCDF, ~7 km monthly | Copernicus Marine | Habitat suitability |
| CMEMS Ocean Colour (SPM, KD490) | `data/processed/physical/{SPM,KD490}/` | NetCDF, 1 km daily | Copernicus Marine | Turbidity / light availability |
| GEBCO 2026 bathymetry | `GEBCO_2026/GEBCO_2026.nc` | NetCDF, global 15 arc-sec | GEBCO | Depth; reference grid for fine-resolution outputs |
| Copernicus COP30 DEM | `rasters_COP30/output_hh.tif` | GeoTIFF, ~19 m | Copernicus | Terrain, TWI, flow routing |
| EUSeaMap 2023 Energy | `EUSeaMap_2023_Energy/EUSeaMap_2023_Energy.gdb` | File GDB | EMODnet | Wave/current energy class |
| UKSeaMap v2026 | `UKSeaMap-Predictive-Map-v2026-1/...gpkg` | GeoPackage | JNCC | Biozones, substrate, salinity/energy regime |
| BGS Predictive Seabed Sediments v1 | `Predictive_SBS_UK_V1_GeoPackage/...gpkg` | GeoPackage | British Geological Survey | Substrate suitability |
| Kinetic energy currents | `ke_currents_atlantic_withconf/ke_currents_atlantic.tif` | GeoTIFF | Derived/published | Bed shear stress proxy |
| UKCP18 sea level RCP4.5 | `seaLevelAnom_marine-sim_rcp45_ann_2007-2100.nc` | NetCDF | Met Office UKCP18 | Sea level rise projection |
| HydroBASINS level 12 | `hybas_lake_eu_lev01-12_v1c/hybas_lake_eu_lev12_v1c.shp` | Shapefile | HydroSHEDS | Catchment delineation (primary) |
| CEFAS catchments | `CEFAS catchments/...shp` | Shapefile | CEFAS | Catchment delineation (superseded) |
| OS OpenRivers | `OpenRivers/Data/oprvrs_gb.gpkg` | GeoPackage | Ordnance Survey | River network, accumulated risk |
| OS OpenRoads / OpenMap Local | `OpenRoads/`, `opmplc_essh_{nw,nx,ny}/` | Shapefiles | Ordnance Survey | Land use context |
| Hutton soil erosion risk | `Hutton_SoilErosionRisk_and_LUT/` | GeoTIFF + LUT | James Hutton Institute | Diffuse pollution validation |
| Hutton soil leaching / runoff | `Hutton_SoilLeachingPotential/`, `Hutton_SoilRunoffRisk/` | Shapefiles | James Hutton Institute | Diffuse pollution validation |
| Hutton LCA 50K | `Hutton_LCA50K_OpenData/LCA_50K.shp` | Shapefile | James Hutton Institute | Grassland intensity split |
| EU DIR 3s | `eu_dir_3s.tif` | GeoTIFF | EU | Flow direction |
| coastalNet connectivity DB | `oceanographicConnectivity.RData` | RData | Assis et al. 2025 | Larval dispersal connectivity |

### 2.5 Threats and pressures

| Dataset | Path (`data/raw/threats/local/` unless noted) | Format | Source | Used for |
|---|---|---|---|---|
| Putuhena CartOcean UK v1 | `data/raw/Putahena_all_GDB/CartOcean_DS_UK_v1_250801.gdb` | File GDB | Putuhena et al. | Noise, cables, oil, dark vessel activity |
| EMODnet vessel density (all/cargo/service/tanker) | `EMODnet_HA_Vessel_Density_*/` | GeoTIFF monthly | EMODnet | Shipping pressure |
| EMODnet dredge spoil dumping | `EMODnet_HA_DredgeSpoilDumping_20251105/` | Shapefile | EMODnet | Vector threat |
| EMODnet dredging | `EMODnet_HA_Dredging_20251103/` | Shapefile | EMODnet | Vector threat |
| EMODnet dumped munitions | `EMODnet_HA_WasteDisposal_DumpedMunitions_20260609/` | Shapefile | EMODnet | Vector threat |
| EMODnet UWWTD discharge points | `EMODnet_HA_WasteDisposal_UWWTD_DischargePoints_20250718/` | Shapefile | EMODnet | Vector threat |
| EMODnet WWI/WWII wrecks | `EMODnet_HA_Heritage_WW_Wrecks_20241226/` | File GDB | EMODnet | Vector threat |
| UK ports | `UK_Ports_coastal/UK_Ports_coastal.shp` | Shapefile | — | Vector threat; seagrass exclusion |
| Fishing under 12 m | `fishing_under12_data/...Polygon.shp` | Shapefile | Marine Scotland | Fishing pressure |
| Global Fishing Watch | API in `05_localthreats.R` | API | GFW | Fishing effort |
| Non-native macroalgae models | `non-native_macroalgae_kuenm_models_review.nc` | NetCDF | Published models | INNS pressure |

---

## 3. Derived layers — data cube candidates

All EPSG:27700 unless stated. These exist on disk now and need no re-run.

### 3.1 Species (`outputs/species/`)

| Layer | Type | Structure | Created by |
|---|---|---|---|
| `rarefied_richness_surface_marine.tif` | Raster | 532×227, 300 m | `07_plots.R` |
| `rarefied_richness_surface_terrestrial.tif` | Raster | 548×254, 300 m | `07_plots.R` |
| `beta_diversity_surface_marine.tif` | Raster | 532×227, 300 m | `07_plots.R` |
| `beta_diversity_surface_terrestrial.tif` | Raster | 548×254, 300 m | `07_plots.R` |
| `community_similarity_surface_marine.tif` | Raster | 532×227, 300 m | `07_plots.R` |
| `community_similarity_surface_terrestrial.tif` | Raster | 498×302, 300 m | `07_plots.R` |
| `richness_grid_marine.gpkg` | Polygon | 216 features, 41 fields, **EPSG:4326**, 5 km grid | `02_themeD_species.R` |
| `richness_grid_terrestrial.gpkg` | Polygon | 219 features, 41 fields, **EPSG:4326**, 5 km grid | `02_themeD_species.R` |
| `protection_priority_marine.gpkg` | Polygon | 216 features, 48 fields | `07_plots.R` |
| `protection_priority_terrestrial.gpkg` | Polygon | 219 features, 48 fields | `07_plots.R` |
| `indicator_records.gpkg` | Point | 13,447 features, **EPSG:4326** | `02_themeD_species.R` |
| `inns_records.gpkg` | Point | 18,565 features, **EPSG:4326** | `02_themeD_species.R` |
| `indicator_density_grid.gpkg` | Polygon | 2 layers: habitat_builders (331), associated_community (453), **EPSG:4326** | `02_themeD_species.R` |

Richness is Hurlbert-rarefied to a common sample size (`vegan::rarefy`), then IDW-interpolated.
Beta diversity is mean Sørensen dissimilarity to queen-contiguous neighbours, IDW-interpolated
(note: raw, not effort-corrected — flagged as deferred in `07_plots.R:2330`).

### 3.2 Habitats (`outputs/habitats/`)

| Layer | Type | Structure | Created by |
|---|---|---|---|
| `scamp_all_habitats.gpkg` | MultiPolygon | 3,786 features | `03_habitats.R` |
| `known_sites_all_habitats.gpkg` | Mixed | 3,890 features | `07_plots.R` |
| `seagrass_potential_solway.gpkg` | Polygon | 70 features | `04_physical.R` |
| `oyster_potential_solway.gpkg` | Polygon | 101 features | `04_physical.R` |
| `sand_dunes_habmos_solway.gpkg` | Mixed | 3,409 features, 16 fields | `03_habitats.R` |
| `sand_dune_restoration_potential.gpkg` | Mixed | 100 features | `03_habitats.R` |
| `seagrass_restoration_opportunity.gpkg` | Polygon | 113 features | `07_plots.R` |
| `oyster_restoration_opportunity.gpkg` | Polygon | 102 features | `07_plots.R` |
| `seagrass_{potential,restoration,recovery,monitoring}_*.tif` | Raster | 1677×842, 100 m | `04_physical.R` / `07_plots.R` |
| `oyster_{potential,restoration,recovery,monitoring}_*.tif` | Raster | 586×294, 286 m | `04_physical.R` / `07_plots.R` |
| `sand_dune_{potential,restoration,recovery,monitoring}_*.tif` | Raster | 8612×4358, 19 m | `03_habitats.R` / `07_plots.R` |
| `seagrass_action_map.tif` / `oyster_action_map.tif` | Raster | as above | `07_plots.R` |

### 3.3 Physical (`outputs/physical/`)

| Layer | Type | Structure | Created by |
|---|---|---|---|
| `gebco_bathy_solway.tif` | Raster | 586×294, 286 m | `04_physical.R` |
| `cop30_dem_solway.tif` | Raster | 8612×4358, 19 m | `04_physical.R` |
| `ke_currents_solway.tif` | Raster | 815×408, 206 m | `04_physical.R` |
| `spm_mean_solway.tif` | Raster | 1677×842, 100 m | `04_physical.R` |
| `salinity_p90_mean_solway.tif` | Raster | 1677×842, 100 m | `04_physical.R` |
| `sal_p10_solway.tif`, `sal_min_monthly_mean_solway.tif`, `sal_mean_days_below20_solway.tif` | Raster | 586×294, 286 m | `04_physical.R` |
| `ukcp18_sea_level_solway.tif` | Raster | 15×7, ~11 km, **846 bands** (9 percentiles × 94 years) | `04_physical.R` |
| `ukseamap_{habitats,biozones,seabed_substrate,salinity_regime,energy_regime}_solway.gpkg` | Polygon | 866–10,906 features | `04_physical.R` |
| `energy_class_solway.gpkg` | Polygon | 730 features | `04_physical.R` |
| `{seagrass,oyster,saltmarsh}_connectivity_hexagons.gpkg` | Polygon | 19 features each, 5 fields, **EPSG:4326**, 8.45 km hexagons | `04_physical.R` |
| `coastalnet_{seagrass,oyster,saltmarsh}_connectivity.rds` | R object | Pairwise probability matrices | `04_physical.R` |

Connectivity is coastalNet forward direct dispersal probability, 180-day window,
stepping-stone disabled (see note at `04_physical.R:2091`).

### 3.4 Threats (`outputs/threats/`)

| Layer | Type | Structure | Created by |
|---|---|---|---|
| `cumulative_pressure_with_pollution.tif` | Raster | 587×294, 286 m — **current version** | `05_localthreats.R` Step 9 |
| `cumulative_pressure_with_pollution_1km.tif` | Raster | 160×68, 1 km — superseded | `05_localthreats.R` |
| `cumulative_pressure_1km.tif` | Raster | 160×68, 1 km — no pollution | `05_localthreats.R` |
| `pressure_layers_normalised_1km.tif` | Raster | 160×68, 1 km, **14 bands** | `05_localthreats.R` |
| `pressure_count_1km.tif` | Raster | 160×68, 1 km | `05_localthreats.R` |
| `layers_raw/*.tif` (14 files) | Raster | 160×68, 1 km each | `05_localthreats.R` |
| `vector_threats_unified.gpkg` | Mixed | 239 features | `05_localthreats.R` |
| `catchment_landcover.gpkg` | Mixed | 71 features, 25 fields | `05_localthreats.R` |

The 14 normalised pressure bands are: vessel_all, vessel_cargo, vessel_service,
vessel_tanker, macroalgae_inns, fishing_gfw, fishing_under12, noise_MoD,
noise_multibeam, noise_subbottom, infra_cables, infra_oil, dark_fishing,
dark_nonfishing. Diffuse pollution is added as a 15th input at the fine
resolution in Step 9.

### 3.5 Catchments and runoff (`outputs/catchments/`)

| Layer | Type | Structure | Created by |
|---|---|---|---|
| `catchments_with_risk.gpkg` | Mixed | 71 features, 26 fields | `05_localthreats.R` |
| `hybas_lev12_runoff_risk.gpkg` | Mixed | 71 features, 19 fields | `05_localthreats.R` |
| `rivers_with_accumulated_risk.gpkg` | MultiLineString | 7,650 features, 12 fields | `05_localthreats.R` |
| `coastal_pour_points.gpkg` | Point | 158 features, 21 fields | `05_localthreats.R` |
| `diffuse_pollution_risk.tif` | Raster | 9227×3193, 20 m | `05_localthreats.R` |
| `hlcm_source_runoff_risk.tif` | Raster | 9227×3193, 20 m | `05_localthreats.R` |
| `marine_pollution_exposure.tif` | Raster | 587×294, 286 m | `05_localthreats.R` Step 8 |
| `twi.tif`, `slope.tif`, `sca.tif`, `flow_dir.tif`, `risk_accumulated.tif` | Raster | 580×336, 290 m (GEBCO-derived) | `05_localthreats.R` |
| `*_cop30.tif` variants | Raster | 8291×4385, 19 m (COP30-derived) | `05_localthreats.R` |
| `{erosion,leaching,runoff,sediment,nutrient_n,nutrient_p}_*.tif` | Raster | 8343×4563, 20 m | `05_localthreats.R` |

Diffuse pollution risk follows SCIMAP methodology: EUNIS land cover → source risk
weight × TWI hydrological connectivity, accumulated downstream via HydroSHEDS flow
direction. The land-cover-to-risk crosswalk is documented inline at
`05_localthreats.R:2012-2052`.

---

## 4. Notes for the data cube

**Grid alignment.** Outputs sit on several different grids:

| Grid | Resolution | Extent origin | Layers |
|---|---|---|---|
| Pressure grid | 1 km | 173067, 507998 | threat layers, pressure stack |
| GEBCO-derived fine | 286 m | 169792/169860, 504634 | cumulative pressure, oyster, bathymetry, salinity |
| Species IDW | 300 m | 173068, 507998 | richness, beta diversity, similarity |
| Seagrass/SPM | 100 m | 169860, 504634 | seagrass scores, SPM, salinity p90 |
| HLCM/land | 20 m | 163060/195920, 530220 | land cover, diffuse pollution |
| COP30 terrain | 19 m | 173211/194397, 505223 | DEM, dune scores, terrain derivatives |

A cube will need a single target grid. The 286 m GEBCO-derived grid is the most
common choice for marine layers and is what the current cumulative pressure uses.

**Reprojection needed.** These are EPSG:4326 and must be reprojected to 27700
before joining a BNG cube: `richness_grid_{marine,terrestrial}.gpkg`,
`indicator_records.gpkg`, `inns_records.gpkg`, `indicator_density_grid.gpkg`,
`*_connectivity_hexagons.gpkg`, `all_species_records.gpkg`, and the AOI shapefiles
under `data/processed/boundaries/aois_*/`.

**Licence review before redistribution.** Check terms for: UNEP-WCMC seagrass and
saltmarsh, Environment Agency Seagrass_Potential, NatureScot GeMS/PMF, Putuhena
CartOcean, James Hutton Institute soil layers, Ordnance Survey products, and
GEBCO. Several are research/non-commercial or require attribution; some may not
permit onward sharing of derived products.

**Missing or broken on disk.**

| Item | Status |
|---|---|
| `data/processed/habitats/hlcm_habmos_corrected.tif` | Referenced by `05_localthreats.R` and `07_plots.R` (woodland in AOI panels) but absent. Only `hlcm_habmos_lca_reclass.tif` exists. |
| `data/processed/habitats/hlcm_lc_lookup.csv` | Referenced but absent. Only `hlcm_lc_lookup_reclass.csv` exists. |
| `outputs/catchments/rivers_with_accumulated_risk_hydrosheds.gpkg` | **Blocking.** Written by `05_localthreats.R:2687`, read by `07_plots.R` for the overview and AOI threat panels. Absent. The older `rivers_with_accumulated_risk.gpkg` on disk is from the superseded COP30 D8 method, not a substitute. |
| `outputs/catchments/hybas_lev12_dg.gpkg` | Written and read within `07_plots.R` (lines 3538 / 3780). Absent, but self-generated on a full run. |
| `data/processed/boundaries/aois_named.gpkg` | Present (64 KB) but fails to open — likely locked by an open R session or written empty. Regenerate from `07_plots.R`. |
| `data/processed/species/pmf_points.gpkg`, `pmf_polygons.gpkg` (and `_scamp` variants) | Referenced by `02_themeD_species.R`, absent. |
| `data/processed/species/taxonomic_lookup.csv`, `qa_excluded.csv` | QA outputs referenced, absent. |
| `data/processed/species/oyster_historical_extent.gpkg` | Written conditionally by `rebuild_clusters.R`, absent. |

**Superseded layers.** `cumulative_pressure_with_pollution_1km.tif` and
`cumulative_pressure_1km.tif` are earlier versions of the current fine-resolution
`cumulative_pressure_with_pollution.tif`. CEFAS catchment outputs are superseded by
the HydroBASINS level 12 versions. The `*_cop30.tif` terrain set is an alternative
to the GEBCO-derived set — both retained, COP30 is the finer.

---

## 5. Figures

Final report figures are written to `figures/report draft/` (numbered 01–42 for the
main sequence) and `figures/report draft/priority/` (107 files, priority layers
L0–L5). Partner-facing panels are in `figures/ppt_panels/`, including per-AOI
four-panel sets (threats, habitats, richness, beta diversity) in subfolders named
by bay, plus the Solway overview maps and hexagon connectivity maps.

AOI bays, in order as stored: Loch Ryan, Fleet Bay, Kirkcudbright Bay, Auchencairn
Bay (`aois_possibole`), plus Luce Bay and Wigtown Bay (`aois_luce_wigtown`).

---

## 6. Figure manifest

All figures come from `scripts/07_plots.R` unless stated. `save_plot()` is defined
at `07_plots.R:122` and writes to `out_dir <- "figures/report draft"` at 300 dpi.

### 6.1 `figures/report draft/` — main numbered series

| File(s) | Shows | Code | Layers drawn |
|---|---|---|---|
| `01_study_area.png`, `01a_study_area_osm.png` | SCAMP AOI on OSM tiles | L190–192 | `scamp_aoi`, `smr`, `dg` |
| `01b_study_area_bathymetry.png` | AOI over GEBCO bathymetry | L195–217 | `gebco_bathy_solway.tif`, boundaries |
| `02_designated_sites.png` | SAC/SPA/SSSI/Ramsar/MPA/NNR/LNR/RSPB/IBA overlay | L226–233 | all `*_scamp.gpkg` designation layers |
| `03`–`06` seagrass combined / known / potential / indicators | Seagrass extent, modelled potential, Zostera indicator records | L256–297 | `seagrass_polygons_harmonised`, `seagrass_potential_solway`, `indicator_records` | 
| `07_seagrass_cluster_NN.png` | Per-cluster zoom, one per seagrass cluster box, titled by bay name | L300–316 (loop) | `seagrass_cluster_boxes`, known + potential seagrass |
| `08_seagrass_suitability_annex.png` | Multi-panel suitability criteria annex | L330–356 | seagrass model inputs (depth, substrate, salinity, energy) |
| `09`–`11` saltmarsh combined / known / indicators | Saltmarsh extent and indicator species | L375–400 | `saltmarsh_harmonised`, `saltmarsh_indicator_points` |
| `12_saltmarsh_cluster_NN.png` (14) | Per-cluster saltmarsh zooms | L405–415 (loop) | `saltmarsh_cluster_boxes` |
| `13`–`17` oyster combined / known / historical / potential / indicators | Native oyster: current extent, historical beds, modelled potential, Ostrea records | L446–492 | `oyster_known_extent`, `oyster_historical_locations`, `oyster_potential_solway` |
| `18_oyster_cluster_NN.png` | Per-cluster oyster zooms | L497–515 (loop) | `oyster_cluster_boxes` |
| `19_oyster_suitability_annex.png` | Oyster suitability criteria annex | L530–555 | oyster model inputs |
| `20`–`23` sand dunes combined / known / indicators / restoration | Dune extent from HabMoS, indicators, restoration potential | L581–614 | `sand_dunes_habmos_solway`, `sand_dune_restoration_potential` |
| `24_sand_dunes_cluster_NN.png` | Per-cluster dune zooms | L619–629 (loop) | `dune_cluster_boxes` |
| `24a_woodland_hlcm.png` | Coastal woodland from HLCM EUNIS G-classes (broadleaved, conifer plantation, mixed, felled) | L747–757 | `hlcm_habmos_corrected.tif` reclassified |
| `24b_woodland_habmos.png` | Same, from HabMoS polygons | L767–771 | `woodland_habmos_solway.gpkg` |
| `25a`/`25b` … `30a`/`30b` | **Paired marine (`a`) / terrestrial (`b`) metric series.** One pair per metric, all from the same loop: 5 km grid downscaled nearest-neighbour to 500 m and coastline-masked. | L801–847 (`metric_specs` loop) | `richness_grid_marine`, `richness_grid_terrestrial` |
| — `25` | Species Richness (mako) | | `species_richness` |
| — `26` | Recording Effort (inferno) | | `record_count` |
| — `27` | Median Recording Year (diverging, midpoint 2010) | | `median_year` |
| — `28` | Most Recent Record (plasma) | | `most_recent_year` |
| — `29` | Beta Diversity (turbo) | | `beta_diversity` |
| — `30` | Community Similarity (cividis) | | `community_similarity` |
| `31_inns_detections_grid.png` | INNS detections per grid cell | L876–881 | `inns_records.gpkg` |
| `32_inns_first_detection.png` | Year of first INNS detection per cell | L895–900 | `inns_records.gpkg` |
| `33_inns_cumulative_spread.png` | Cumulative INNS spread through time bands, with inset | L939–982 | `inns_records.gpkg` |
| `34_<species>.png` (11) | One map per indicator species — common eelgrass, dwarf eelgrass, glasswort, common reed, saltmarsh grass, sea lavender, sea purslane, thrift, marram grass, lyme grass, native oyster | L1006–1017 (loop) | `indicator_records.gpkg` filtered by species |
| `35_pmf_records.png` | Priority Marine Feature records from GeMS | L1034–1051 | `pmf_points`, `pmf_polygons` |
| `36_fishing_pressure.png` | GFW effort + under-12 m fishing, two panels | L1122–1170 | `fishing_gfw.tif`, `fishing_under12.tif` |
| `37_river_accumulated_risk.png` | SCIMAP river accumulated diffuse risk | L1187–1197 | `rivers_with_accumulated_risk_hydrosheds.gpkg` |
| `38_marine_pollution_exposure.png` | Sea-path decay exposure from pour points | L1202–1218 | `marine_pollution_exposure.tif` |
| `39_diffuse_pollution_combined.png` | Catchment risk + rivers + pour points + marine exposure | L1225 | catchment and river risk layers |
| `40_cumulative_pressure_combined.png` | Cumulative pressure with diffuse pollution included | L1259–1289 | `cumulative_pressure_with_pollution.tif` |
| `41_cumulative_pressure.png` | Cumulative marine pressure alone | L1292–1307 | `cumulative_pressure_with_pollution.tif` |
| `42_vector_threats.png` | Point/polygon threats by type | L1312–1324 | `vector_threats_unified.gpkg` |

**Numbering collisions.** `31` and `32` are each used twice — once for INNS
(`31_inns_detections_grid`, `32_inns_first_detection`, written by the INNS block)
and once for threats (`31_river_accumulated_risk`, `32_marine_pollution_exposure`,
apparently from an earlier numbering). The threats versions duplicate `37` and `38`.

**Absent from disk** despite being in the code: `03`–`06` (seagrass combined,
known, potential, indicators), `19_oyster_suitability_annex`, `35_pmf_records`.
Seagrass clusters `01`–`06` and oyster/dune clusters beyond `01` are also missing.
`25_marine_richness.png` is a legacy file predating the `25a`/`25b` convention.

### 6.2 `figures/report draft/priority/` — prioritisation layers

The framework runs L0 → L5. Marine layers (seagrass, oyster) use the fine
suitability rasters; terrestrial (dune) uses its own merged Layer 3.

| File(s) | Shows | Code |
|---|---|---|
| `L0_pressure_count.png` | Count of distinct pressures per cell — management complexity | L1773–1781 |
| `L1a`–`L1e` reference sites | Least-impacted examples per habitat: saltmarsh, sand dunes, woodland, oyster, seagrass | L2214–2308 |
| `L1f_reference_sites_all.png`, `L1_reference_sites.png` | All reference sites combined | L2563–2572 |
| `L1e0_raw_richness_vs_effort.png` | Scatter: raw richness against record count, showing effort bias | L1964–1973 |
| `L1e0a_raw_richness_map.png`, `L1e0b_record_count_map.png` | Pre-rarefaction richness and effort maps | L1981–1993 |
| `L1e_rarefaction_curve.png` | Rarefaction curves for a fixed cell set at increasing sample sizes | L2054–2064 |
| `L1e1_rarefaction_effect.png` | Raw vs rarefied richness per cell | L2114–2126 |
| `L1e1a_rarefied_richness_map.png` | Rarefied richness, 5 km cells | L2130–2136 |
| `L1e1b_rarefaction_shift_map.png` | Where rarefaction changes the picture most (raw − rarefied) | L2141–2146 |
| `L1e2_rarefied_richness_surface.png` | IDW-interpolated rarefied richness surface | L2313–2325 |
| `L2a`–`L2d` | Marine and terrestrial ecological value category, then protection priority (high value × outside designation) | L2737–2776 |
| `L3a_{seagrass,oyster}_restoration_score.png` | Active restoration priority: suitable + low pressure + connected | L2859–2870 |
| `L3b_{seagrass,oyster}_recovery_score.png` | Natural recovery support: mid-pressure peak + diversity + connectivity | L2975–2986 |
| `L3b_beta_diversity_surface.png` | Marine beta diversity IDW surface | L2377–2388 |
| `L3_{seagrass,oyster}_opportunity_connectivity.png` | Restoration opportunity cells coloured by coastalNet score | L2642–2660 |
| `L3_terrestrial_*` (3) | Terrestrial rarefaction curve, richness surface, beta surface | L2458–2552 |
| `L3_dune_dist_score.png`, `L3_dune_builder_score.png`, `L3_dune_restoration_score.png`, `L3a_dune_restoration_score.png`, `L3b_dune_recovery_score.png` | Dune Layer 3 components and combined | L3229–3302 |
| `L4_{seagrass,oyster}_monitoring_score.png`, `L4_dune_monitoring_score.png` | Monitoring priority | L3076–3087 |
| `L5_{seagrass,oyster}_action_map.png` | Dominant action class per cell | L3156–3167 |
| `L3a_`/`L3b_`/`L4_`/`L5_` `seagrass_cluster_NN_<bay>.png` (6 each) | Per-bay zooms: Loch Ryan, Wigtown Bay, Fleet Bay, Kirkcudbright Bay, Auchencairn Bay & Rough Island, Nith Estuary | L3332–3359 (loop) |
| `L3a_`/`L3b_`/`L4_`/`L5_` `oyster_cluster_NN.png` (7 each) | Per-cluster oyster zooms | same loop |
| `L3_dune_cluster_NN.png` (13) | Per-cluster dune zooms | same loop |
| `known_sites_all_habitats.png` | All five habitats, known + potential, for TAG circulation | L3720–3737 |
| `hybas12_runoff_risk.png` | Diffuse runoff risk by HydroBASINS level 12 sub-catchment | L3901–3915 |

### 6.3 `figures/ppt_panels/` — partner-facing panels

Written with direct `ggsave()` rather than `save_plot()`.

| File(s) | Shows | Code |
|---|---|---|
| `overview_pressure_catchments.png` | Solway overview: cumulative pressure + catchment runoff risk + rivers (Low/High) + AOI boxes, CartoDB Positron base, scale bar, north arrow | L3985–4020 |
| `overview_richness.png` | Same style, IDW rarefied richness surface (mako) | L4514–4528 |
| `overview_beta_diversity.png` | Same style, beta diversity surface (inferno) | L4532–4546 |
| `aoi3_a_threats.png` | AOI 3 Kirkcudbright — cumulative pressure, catchment risk, river risk, vector threat symbols, SMR boundary | L4136–4180 |
| `aoi3_b_habitats.png` | Known habitats filled, potential dashed, historical oyster markers, rivers | L4182–4210 |
| `aoi3_c_richness.png` | Rarefied richness with catchment outlines and rivers for context | L4212–4228 |
| `aoi3_d_beta_diversity.png` | Beta diversity, same context layers | L4230–4246 |
| `aoi3_combined.png` | The four panels via patchwork | L4248–4258 |
| `<Bay>/a_threats.png`, `b_habitats.png`, `c_richness.png`, `d_beta_diversity.png` | Same four panels regenerated per AOI into a folder named by bay | L4299–4490 (batch loop) |
| `hexagon_connectivity_seagrass.png` | coastalNet hexagons by status (Existing/Potential/Neither) with curved directional arrows between hexagon pairs, coloured and weighted by dispersal probability | L4639–4660 |
| `hexagon_connectivity_oyster.png` | Same for native oyster | L4661–4682 |

Existing subfolders are `AOI_1`, `AOI_2`, `AOI_4` and `Kirkcudbright` — these
predate the `bay_name` attribute being written onto the AOI polygons, so a re-run
will produce `Loch Ryan`, `Luce Bay`, `Wigtown Bay`, `Fleet Bay`,
`Kirkcudbright Bay`, `Auchencairn Bay` instead.

### 6.4 Other figure locations

| Folder | Contents | Script |
|---|---|---|
| `figures/` (root) | `map1a`–`map1h`, `map2a`–`map2b`, `map3` — paired marine/terrestrial species panels as PNG + PDF | `02_themeD_species.R:1467` (`save_paired`) |
| `figures/inns_species/{marine,terrestrial,freshwater}/` | One map per INNS species (109 files) | `02_themeD_species.R:1799` |
| `figures/indicator_species/habitat_builders/` | Per-species habitat-builder maps (11) | `02_themeD_species.R` |
| `figures/habitats/` | Source comparison and model-step diagnostics (106) | `03_habitats.R` |
| `figures/physical/` | CMEMS summaries, connectivity (12) | `04_physical.R` |
| `figures/threats/`, `.../individual/`, `.../layers_raw/`, `.../layers_normalised/` | Per-pressure diagnostics (68) | `05_localthreats.R` |
| `figures/catchments/` | Diffuse pollution comparison | `05_localthreats.R` |

---

## 7. Field dictionaries

Field names read with `LIMIT 0` queries — no geometry or rows loaded.

### 7.1 `protection_priority_marine.gpkg` / `protection_priority_terrestrial.gpkg`

216 / 219 polygons, 48 fields, EPSG:27700, 5 km cells.

| Field | Meaning | Units |
|---|---|---|
| `cell_id` | Grid cell identifier | — |
| `species_richness` | Distinct species recorded in cell | count |
| `record_count` | Total occurrence records | count |
| `gbif_records`, `obis_records` | Records by source database | count |
| `median_year`, `most_recent_year` | Recording recency | year |
| `beta_diversity` | Mean Sørensen dissimilarity to queen-contiguous neighbours | 0–1 |
| `community_similarity` | Mean Sørensen similarity to all other cells | 0–1 |
| `richness_<Phylum>` (34 marine / 35 terrestrial) | Per-phylum richness — Annelida, Arthropoda, Chordata, Echinodermata, Mollusca, Bryozoa, Cnidaria, Rhodophyta, Porifera, Tracheophyta, Ascomycota, Basidiomycota etc. | count |
| `high_sim`, `high_rich` | Above-threshold flags for similarity and richness | logical |
| `value_category` | Ecological value class derived from the two flags | category |
| `pressure_count` | Number of distinct marine pressures in the cell | count |
| `protected` | Cell intersects a designated site | logical |
| `protection_gap` | High value but not protected | logical |
| `protection_priority` | Final priority flag | logical |

`richness_grid_marine.gpkg` / `richness_grid_terrestrial.gpkg` (EPSG:4326, 41
fields) carry the same fields up to the per-phylum block, without the
`high_sim` → `protection_priority` decision columns.

### 7.2 `catchments_with_risk.gpkg`

71 polygons, 26 fields. HydroBASINS level 12 with risk attributes appended.

| Field | Meaning | Units |
|---|---|---|
| `HYBAS_ID`, `NEXT_DOWN`, `NEXT_SINK`, `MAIN_BAS` | HydroBASINS topology identifiers | — |
| `DIST_SINK`, `DIST_MAIN` | Distance to sink / main stem | km |
| `SUB_AREA`, `UP_AREA` | Sub-basin and upstream contributing area | km² |
| `PFAF_ID` | Pfafstetter code | — |
| `SIDE`, `LAKE`, `ENDO`, `COAST`, `ORDER`, `SORT` | HydroBASINS attributes (bank side, lake, endorheic, coastal, order) | — |
| `local_auth`, `code` | Local authority name and code | — |
| `hectares`, `total_ha` | Area of the sub-catchment covered by the land cover and risk rasters; 0 where the sub-catchment lies wholly outside the study area | ha |
| `mean_runoff_risk` | Mean SCIMAP land-cover source risk weight | 0–1 |
| `mean_diffuse_risk` | Mean source risk × TWI connectivity | index |
| `mean_source_risk` | Mean source risk before connectivity weighting | 0–1 |
| `max_pour_accum` | Highest accumulated risk at any pour point | index |
| `total_pour_accum_raw` | Summed raw accumulated risk across pour points | index |
| `n_pour_points` | Pour points draining this catchment | count |

`hybas_lev12_runoff_risk.gpkg` (19 fields) is the same layer before the pour-point
and diffuse-risk joins — `mean_runoff_risk` is its only risk field.

### 7.3 `coastal_pour_points.gpkg`

158 points, 21 fields. River outlets where accumulated catchment risk reaches the coast.

| Field | Meaning | Units |
|---|---|---|
| `watercourse_name` | River name from OS OpenRivers | — |
| `accum_risk_raw` | Accumulated diffuse risk at outlet, unscaled | index |
| `accum_risk` | Rescaled accumulated risk | 0–1 |
| `area_sq_km` | Contributing area | km² |
| `HYBAS_ID`, `S_CatNAME` | Catchment identifiers | — |
| `inspireid`, `feature`, `jrsdtn`, `nation`, `objnam`, `nobjnm`, `entitycode`, `instancecd`, `inform`, `scamax`, `txtdsc`, `sordat`, `sorind`, `copyright`, `metadata` | Pass-through attributes from the source coastline / watercourse feature dataset | — |

`coastal_pour_points_raw.gpkg` (204 points) is the pre-filter version.

### 7.4 `vector_threats_unified.gpkg`

239 features, 1 field. Deliberately minimal — source attributes are dropped on merge.

| Field | Meaning | Values |
|---|---|---|
| `threat_type` | Threat category | `dredge_spoil`, `dredging`, `dumped_munitions`, `ports`, `waste_discharge`, `ww_wrecks` |

### 7.5 `catchment_landcover.gpkg`

71 polygons, 25 fields. HydroBASINS attributes (as 7.2) plus composition:

| Field | Meaning | Units |
|---|---|---|
| `n_classes` | Distinct EUNIS land-cover classes present | count |
| `total_pixels` | HLCM pixels in catchment | count (20 m pixels) |
| `dominant_class` | Most extensive EUNIS code | — |
| `dominant_desc` | Its description | — |
| `dominant_pct` | Share of catchment under the dominant class | % |

### 7.6 Remaining `outputs/` layers

| Layer | Fields |
|---|---|
| `rivers_with_accumulated_risk*.gpkg` | `id`, `flow_direction`, `length` (m), `fictitious`, `form`, `watercourse_name`, `watercourse_name_alternative`, `start_node`, `end_node`, `seg_length_m`, `accum_risk_raw`, `accum_risk` |
| `outputs/threats/pour_points.gpkg` | As rivers, plus `S_Catchmt`, `S_CatNAME`; no `seg_length_m` |
| `{seagrass,oyster}_potential_solway.gpkg` | `potential` (1 = suitable), `area_ha` |
| `{seagrass,oyster}_restoration_opportunity.gpkg` | `habitat`, `connectivity_score` (coastalNet max probability, 0–1) |
| `sand_dune_restoration_potential.gpkg` | `tier` (restoration tier), `area_ha` |
| `sand_dunes_habmos_solway.gpkg` | HabMoS schema: `HABITAT_CODE`, `HABITAT_NAME`, `EUNIS_COMPONENT`, `COMPONENT_NAME`, `HAB_PROP`, `HAB_AREA`, `POLYGON_ID`, `MOSAIC`, `HAB_NUMBER`, `HABMOS_SOURCE`, `SURVEY_ID`, `SURVEY_TYPE`, `SURVEY_DATE`, `DATA_OWNER`, `SYMBOL`, `area_ha` |
| `scamp_all_habitats.gpkg` | `habitat`, `type` (known/potential), `area_ha` |
| `known_sites_all_habitats.gpkg` | `class` (habitat × known/potential label) |
| `{seagrass,oyster,saltmarsh}_connectivity_hexagons.gpkg` | `ID`, `IDSEQ`, `HEX` (coastalNet hexagon ids), `status` (Existing/Potential/Neither), `connectivity_score` (0–1) |
| `ukseamap_habitats_solway.gpkg` | `sal_class`, `en_class`, `bz_class`, `sub_class`, `fulldetail`, `mhc_l2`–`mhc_l4`, `mhc_code`, `mhc_name`, `eunis07_l2`–`l4`, `eunis_code`, `eunis_name`, `in_offshore`, `country`, `jnccid`, `updatedate`, `changetype` |
| `ukseamap_{biozones,energy_regime,salinity_regime,seabed_substrate}_solway.gpkg` | One class field (`bz_class` / `en_class` / `sal_class` / `sub_class`) plus `jnccid`, `updatedate`, `changetype` |
| `energy_class_solway.gpkg` | `Energy`, `Shape_Length`, `Shape_Area` |
| `indicator_records.gpkg` | `latin_name`, `year`, `coord_uncertainty` (m), `basis_of_record`, `kingdom`, `phylum`, `class`, `source`, `common_name`, `habitat` |
| `inns_records.gpkg` | As indicator records, with `realm` and `source_list` instead of `habitat` |
| `indicator_density_grid.gpkg` | Two layers (`habitat_builders`, `associated_community`): `cell_id`, `n_indicator` |
| `cefas_catchments_crop.gpkg` | `S_Catchmt`, `S_CatNAME`, `COUNT_S_ID`, `Shape_Leng`, `Shape_Area` (all stored as character — superseded layer) |

---

## 8. Sizes

Measured file sizes, MB.

| Group | Layers | Size |
|---|---|---|
| Species (3.1) | 13 | 8.7 MB |
| Habitats (3.2) | 21 | 35.6 MB |
| Physical (3.3) | 17 | 61.6 MB |
| Threats (3.4) | 21 | 1.7 MB |
| Catchments (3.5), excluding COP30 duplicates | 16 | 126.5 MB |
| **Curated cube total** | **88** | **≈ 234 MB** |

Largest single contributors: `diffuse_pollution_risk.tif` 55.6 MB,
`cop30_dem_solway.tif` 32.6 MB, `known_sites_all_habitats.gpkg` 11.7 MB,
`sand_dune_*_score.tif` 2.8–3.0 MB each, `ukseamap_habitats_solway.gpkg` 9.4 MB.

**The whole `outputs/` tree is 745.6 MB.** The difference is almost entirely the
COP30 terrain duplicates (`*_cop30.tif`: 499.5 MB across 8 files) plus the
intermediate pollution component rasters (`erosion`, `leaching`, `runoff`,
`sediment`, `nutrient_n`, `nutrient_p`, `grassland_intensity`: 55.7 MB). Neither
set is needed in a partner handover if the final `diffuse_pollution_risk.tif` and
`risk_accumulated.tif` are included.

A single GeoPackage holding only the vector layers from section 3 would be
roughly 60 MB before compression; rasters cannot go in a GeoPackage efficiently
and are better shipped as Cloud-Optimised GeoTIFFs alongside it. Expect
**≈ 234 MB uncompressed, ≈ 120–150 MB** as a zipped handover once the COP30 and
intermediate sets are excluded.
