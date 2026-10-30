# SCAMP GIS — Reproducibility Improvement Tasks

Ordered by priority. P1 items block partners from re-running the pipeline at all;
P2 items make it reliable; P3 items make it pleasant.

Companion document: [PROJECT_INVENTORY.md](PROJECT_INVENTORY.md).

---

## P1 — Blocking

### 1. Regenerate the missing river risk layer

`outputs/catchments/rivers_with_accumulated_risk_hydrosheds.gpkg` is read by
`07_plots.R` for the Solway overview map and every AOI threats panel, but is not on
disk. The file present (`rivers_with_accumulated_risk.gpkg`) comes from the
superseded COP30 D8 method and is not interchangeable.

- Re-run `05_localthreats.R` Step 4 (revised, HydroSHEDS DIR) to produce it
- Confirm `07_plots.R` reads the intended file and not the older one

### 2. Resolve the missing HLCM land cover files

`data/processed/habitats/hlcm_habmos_corrected.tif` and `hlcm_lc_lookup.csv` are
read by both `05_localthreats.R` and `07_plots.R` (woodland classification in the
AOI habitat panels) but are absent. Only `hlcm_habmos_lca_reclass.tif` and
`hlcm_lc_lookup_reclass.csv` exist.

- Determine whether `_reclass` supersedes `_corrected`, or whether `_corrected` is a genuinely missing build step
- Either regenerate `_corrected`, or update the ~6 call sites to point at `_reclass`
- Whichever is authoritative, note it in the inventory

### 3. Break the `03_habitats.R` ↔ `04_physical.R` circular dependency

`04_physical.R` needs harmonised habitats from `03_habitats.R`; `03_habitats.R`
sections 1e and 3e need potential models from `04_physical.R`. Neither runs cleanly
end-to-end from a clean checkout.

- Split `03_habitats.R` into `03a_habitats_harmonise.R` (produces harmonised extent, no dependency on 04) and `03b_habitats_plots.R` (combined known + potential figures, runs after 04)
- Resulting order becomes strictly linear: 02 → 03a → 04 → 03b → 05 → 07

### 4. Restore or remove the other missing intermediates

| File | Action |
|---|---|
| `data/processed/species/pmf_points.gpkg`, `pmf_polygons.gpkg` (+ `_scamp` variants) | Re-run `02_themeD_species.R` section 8, or drop the references |
| `data/processed/species/taxonomic_lookup.csv`, `qa_excluded.csv` | QA audit trail — regenerate, these are the provenance record for excluded records |
| `data/processed/species/oyster_historical_extent.gpkg` | Conditional write in `rebuild_clusters.R`; confirm still needed |
| `data/processed/boundaries/aois_named.gpkg` | 64 KB but will not open. Delete and regenerate from `07_plots.R` |

---

## P2 — Reliability

### 5. Split the analysis out of `07_plots.R`

At 4,684 lines, `07_plots.R` both computes and plots. It generates the rarefied
richness surfaces, beta diversity surfaces, protection priority grids, restoration
/ recovery / monitoring scores, action maps and an alternative catchment
delineation — then reads several of them back. A partner wanting the richness
surface has to run a plotting script to get it.

- Extract analysis into `06_priority_layers.R` (Layers 0–5, biodiversity surfaces, priority scoring)
- Leave `07_plots.R` as presentation only: read finished layers, draw maps
- This also makes the outputs in section 3 of the inventory traceable to a script whose name matches what it does

### 6. Pin the environment

There is no lockfile, and `00_setup.R` is a bare `install.packages()` call listing
6 packages while the pipeline actually loads 20+ (`sf`, `terra`, `tidyterra`,
`patchwork`, `ggspatial`, `ggnewscale`, `vegan`, `gstat`, `rosm`, `leaflet`,
`coastalNet`, `CoordinateCleaner`, `rgbif`, `robis`, …).

- Add `renv` and commit `renv.lock`
- Record R version (currently 4.6.0) and GDAL/PROJ versions
- `coastalNet` installs from GitHub — pin the commit SHA, not just the repo

### 7. Add a single entry point

- Add `run_all.R` that sources scripts in the correct order with clear stage markers
- Guard each stage so a failure stops the run rather than silently continuing with stale objects
- Alternatively adopt `targets`, which would make the dependency graph explicit and skip work that is already current

### 8. Centralise paths and parameters

Paths are hard-coded as string literals throughout, some repeated a dozen times.
Analysis constants (5 km species grid, 1 km pressure grid, 2 km cluster buffer,
8,047 m AOI buffer, 180-day dispersal window, rarefaction target) are scattered
inline.

- Create `R/config.R` holding every input path and every tunable constant
- Source it at the top of each script
- Partners can then re-point the project at their own data directory by editing one file

### 9. Make the disk-fallback pattern explicit

The `if (!exists("x")) x <- st_read(...)` idiom means an object may come from memory
or disk depending on session state, so results can differ between a fresh run and an
interactive one.

- Prefer always loading from disk at section start, or
- Keep the fallback but log which branch was taken, so the run log records provenance

### 10. Decide on one target grid

Outputs currently sit on six grids (1 km, 300 m, 286 m, 100 m, 20 m, 19 m) with
three different extent origins. For the data cube this has to be resolved.

- Choose the target grid (the 286 m GEBCO-derived grid is the natural candidate for marine layers)
- Add a `harmonise_to_grid()` helper and a final step that emits the aligned stack
- Reproject the EPSG:4326 layers listed in the inventory before they enter the cube

### 11. Remove or clearly quarantine superseded code

`05_localthreats.R` contains sections explicitly labelled DEFUNCT (Step 4 custom D8,
Step 5 nearest-endpoint pour points, Step 6 CEFAS summaries) alongside their REVISED
replacements. `07_plots.R` has a "HUTTON SOIL RISK — DIAGNOSTIC" block and a
"PRIORITY MAPS — PLAN" comment block.

- Delete superseded blocks (git history preserves them), or move to `scripts/archive/`
- Same for the `_cop30` duplicate terrain outputs and the CEFAS catchment outputs if the HydroBASINS versions are authoritative

---

## P3 — Quality of life

### 12. Write a README

There is no README. A partner cloning this has no entry point.

- What the project is, who it is for
- Folder layout, run order, expected runtime
- Where raw data must be placed (it is not in the repo — see item 14)
- Link to `PROJECT_INVENTORY.md` and this file

### 13. Document licences properly

No licence or attribution notes appear anywhere in the scripts. Several sources are
restricted or attribution-bound: UNEP-WCMC, Environment Agency, NatureScot GeMS,
Putuhena CartOcean, James Hutton Institute, Ordnance Survey, GEBCO, Copernicus.

- Add `DATA_LICENCES.md` recording source, licence, attribution string and whether derived products may be redistributed
- Flag any layer that must be excluded from the partner data cube

### 14. Fix repository hygiene

`git status` currently shows `data/`, `outputs/`, `figures/`, `rosm.cache/` and two
`.zip` downloads as untracked, plus a stray `ChatGPT Image ....png` in the root.

- Add a `.gitignore` covering `data/`, `outputs/`, `figures/`, `rosm.cache/`, `*.zip`, `.Rproj.user/`
- Remove the stray image and the two `0006*.zip` GBIF download archives from the working tree
- Decide how partners obtain raw data (download script, shared drive, DOI) since it cannot live in git

### 15. Add validation checks

- After each script, assert the expected outputs exist and are non-empty
- Check CRS is 27700 where expected before writing
- Assert feature counts and raster dimensions against known values, so a silent upstream change is caught

### 16. Tidy naming

- `02_themeD_species.R` both defines the AOI and does species work; the AOI step belongs in a `01_boundaries.R` that actually runs (the current `03_boundaries.R` is exploratory and misnumbered)
- Numbering skips 06
- `Putahena_all_GDB` is spelled `Putuhena` in the literature
- `aois_possibole` appears to be a typo for "possible"
- Scripts live in two places (root and `scripts/`) — consolidate into `scripts/`

### 17. Record method provenance in the outputs

- Write a sidecar `.json` or `.md` next to each derived layer recording the script, date, input files and key parameters
- For the data cube specifically, partners will want to know that richness is Hurlbert-rarefied and IDW-interpolated, that beta diversity is queen-contiguity Sørensen and **not** effort-corrected, and that connectivity is 180-day forward direct dispersal with stepping-stone disabled

---

## Suggested sequence

1. **Make it run** — items 1, 2, 4 (restore missing files), then 3 (break the circular dependency)
2. **Make it repeatable** — items 6, 7, 8 (environment, entry point, config)
3. **Make it shareable** — items 10, 13, 14 (grid, licences, repo hygiene), then 12 (README)
4. **Make it maintainable** — items 5, 11, 16 (split, prune, rename) as capacity allows
