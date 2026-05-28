# Bridger Data Request — Template Schema

Four CSV files, loosely relational. If a GeoPackage (`.gpkg`) is easier on your end we accept that instead; same columns, with geometry in `latitude_wgs84` / `longitude_wgs84` replaced by a point geometry column in WGS84 (EPSG:4326). Dates/times are ISO-8601 in UTC.

| File | Grain | Primary key | Foreign keys | One row = |
|------|-------|-------------|--------------|-----------|
| `01_surveys.csv` | one per survey (flight pass / campaign segment) | `survey_id` | — | One contiguous flight covering a set of wells |
| `02_site_passes.csv` | one per (survey × well flown over) | `site_pass_id` | `survey_id`, `site_id` | A well that was flown, detection or not |
| `03_detections.csv` | one per plume detection | `plume_id` | `site_pass_id`, `survey_id`, `site_id` | One emission plume observed |
| `04_well_metadata.csv` | one per well | `site_id` | — | Stable well attributes |

**Why all four.** Detections alone bias us toward the high end. Site passes tell us which wells were *looked at and not detected* — essential for an unbiased supplement to our 400-well ground campaign. Survey metadata gives us the POD at the moment of each pass.

---

## File 1 — `01_surveys.csv`

| column | type | unit | required | description |
|--------|------|------|----------|-------------|
| `survey_id` | string | — | yes | Unique id for the flight/survey |
| `survey_date_utc` | date | UTC | yes | Date of the flight (YYYY-MM-DD) |
| `basin` | string | — | yes | Marcellus / San Juan / Permian |
| `region` | string | — | yes | State or sub-basin (PA / NM / etc.) |
| `aircraft_id` | string | — | no | Tail number or platform label |
| `flight_start_utc` | datetime | UTC | yes | ISO-8601 |
| `flight_end_utc` | datetime | UTC | yes | ISO-8601 |
| `coverage_wkt` | WKT polygon | WGS84 | yes | Polygon of the area surveyed; used to resolve which un-listed wells are inside the footprint |
| `gcn_ppm_m` | float | ppm-m | yes | Gas concentration noise (Thorpe 2024 Combined model input) |
| `mean_wind_speed_ms` | float | m/s | yes | Flight-averaged, 10-m reference height |
| `mean_wind_direction_deg` | float | deg | no | 0–360, clockwise from north |
| `mean_altitude_ft_agl` | float | ft AGL | no | Flight-average height above ground |
| `atmospheric_stability_class` | string | — | no | Pasquill class (A–F) if measured/modeled |
| `mdl_kgph_at_3ms` | float | kg/h | yes | Headline PoD90 at standard 3 m/s wind for this survey |
| `pod_model_version` | string | — | yes | e.g. `GML2.0` |
| `notes` | string | — | no | Free text |

## File 2 — `02_site_passes.csv`

| column | type | unit | required | description |
|--------|------|------|----------|-------------|
| `site_pass_id` | string | — | yes | Unique id for this (survey × site) event |
| `survey_id` | string | — | yes | FK to `01_surveys.survey_id` |
| `site_id` | string | — | yes | Your internal well id |
| `api_number` | string | — | yes* | 14-digit API-number (*or best available identifier; we'll match to state databases) |
| `latitude_wgs84` | float | deg | yes | Well latitude |
| `longitude_wgs84` | float | deg | yes | Well longitude |
| `pass_timestamp_utc` | datetime | UTC | yes | Moment the aircraft passed over |
| `wind_speed_ms_at_pass` | float | m/s | yes | Local wind at pass, not flight average |
| `wind_direction_deg_at_pass` | float | deg | no | — |
| `altitude_ft_agl_at_pass` | float | ft | no | — |
| `plume_detected` | bool | — | yes | TRUE/FALSE for this pass |
| `pod_at_pass` | float | 0–1 | yes* | PoD at the pass conditions (*if Bridger computes this per-pass) |
| `mdl_kgph_at_pass` | float | kg/h | no | Pass-specific PoD90 |
| `n_passes_in_survey` | int | count | no | If the site was crossed multiple times in one survey |
| `notes` | string | — | no | — |

Wells in the survey footprint that are not listed here will be treated as *not flown*. Please include every pass, including those with no detection.

## File 3 — `03_detections.csv`

| column | type | unit | required | description |
|--------|------|------|----------|-------------|
| `plume_id` | string | — | yes | Unique plume id |
| `site_pass_id` | string | — | yes | FK to `02_site_passes.site_pass_id` |
| `survey_id` | string | — | yes | FK to `01_surveys.survey_id` |
| `detection_timestamp_utc` | datetime | UTC | yes | — |
| `latitude_wgs84` | float | deg | yes | Plume centroid latitude |
| `longitude_wgs84` | float | deg | yes | Plume centroid longitude |
| `horizontal_position_uncertainty_m` | float | m | no | 1-σ position uncertainty |
| `associated_site_id` | string | — | yes | Nearest/attributed well (your best-guess attribution) |
| `associated_api_number` | string | — | yes* | 14-digit API; same caveat as above |
| `emission_rate_kgph` | float | kg/h | yes | Best-estimate emission rate |
| `emission_rate_lo_kgph` | float | kg/h | yes | Lower end of reported uncertainty |
| `emission_rate_hi_kgph` | float | kg/h | yes | Upper end of reported uncertainty |
| `emission_rate_uncertainty_basis` | string | — | yes | e.g. `95% CI via bootstrap`, `±1σ from inversion`, etc. |
| `wind_speed_ms` | float | m/s | yes | Wind at the plume, used in the flux retrieval |
| `wind_direction_deg` | float | deg | no | — |
| `gcn_ppm_m` | float | ppm-m | yes | Used in the PoD model |
| `altitude_ft_agl` | float | ft | no | — |
| `plume_length_m` | float | m | no | Downwind extent |
| `plume_image_filename` | string | — | no | Reference to an attached raster/plot, if provided |
| `classification` | string | — | no | e.g. production / tank / compressor / pipeline — if attributed to a component type |
| `notes` | string | — | no | — |

Multiple detections on the same site_pass_id are fine (e.g. two plumes at one well).

## File 4 — `04_well_metadata.csv`

If Bridger already has this information, great. If not, we can join externally from state records — but including any Bridger-known values reduces ambiguity.

| column | type | unit | required | description |
|--------|------|------|----------|-------------|
| `site_id` | string | — | yes | Matches `02_site_passes.site_id` |
| `api_number` | string | — | yes* | — |
| `basin` | string | — | yes | Marcellus / San Juan / Permian |
| `region` | string | — | yes | — |
| `operator` | string | — | no | — |
| `latitude_wgs84` | float | deg | yes | — |
| `longitude_wgs84` | float | deg | yes | — |
| `spud_date` | date | — | no | — |
| `first_production_date` | date | — | no | — |
| `well_type` | string | — | no | oil / gas / CBM |
| `lift_type` | string | — | no | plunger / rod_pump / gas_lift / natural / ESP |
| `gas_production_mcfpd_12mo` | float | Mcf/d | no | 12-month average daily gas production |
| `oil_production_bopd_12mo` | float | bbl/d | no | 12-month average daily oil production |
| `water_production_bblpd_12mo` | float | bbl/d | no | — |
| `gor_scf_per_bbl` | float | scf/bbl | no | Gas-to-oil ratio |
| `marginal_status` | bool / yes-no | — | no | Marginal-well classification at survey time |
| `notes` | string | — | no | — |

---

## Delivery format

Acceptable, in order of preference:

1. `.gpkg` (one file per table, or one file with multiple layers).
2. `.csv` with headers (this template).
3. `.parquet` if volume is very large.
4. ZIP bundle containing any of the above plus this `README.md`.

CRS: WGS84 (EPSG:4326). All timestamps UTC. Units as specified per column.

## Minimal vs. full delivery

If a full delivery is not feasible, the minimum useful subset is:

- `03_detections.csv` (all required columns)
- `01_surveys.csv` (all required columns, one row per survey)
- `02_site_passes.csv` with at least `site_id`, `survey_id`, `plume_detected`, `latitude_wgs84`, `longitude_wgs84`

Without `02_site_passes.csv` we cannot compute unbiased well-level emission rates — only an upper-tail catalog.
