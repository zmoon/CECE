# MEGAN Biogenic Emissions

## Overview

CECE provides two C++ MEGAN biogenic emission schemes:

- **`megan`** — Single-species isoprene scheme, ported from HEMCO. Supports
  two emission methods:
  - `"megan21"` (default) — configurable MEGAN2.1 isoprene calculation. The
    existing `"native"` spelling remains an equivalent compatibility alias.
  - `"hemco_3_12_1"` — source-pinned HEMCO 3.12.1 isoprene calculation with
    cold-started, evolving per-cell temperature, PAR, and LAI histories by
    default. Set `history_mode: fixed` to use the configured scalar histories.
- **`megan3`** — CECE's 19-class biogenic emission and chemical-speciation scheme. Its C++ runtime uses bulk activity factors; the available multilayer canopy helpers are not yet wired into that emission calculation. This is not a claim of complete or independently validated upstream MEGAN3 science.

Both schemes execute CECE code; `hemco_3_12_1` does not launch HEMCO, and
`megan3` does not launch a standalone MEGAN3 model. The HEMCO method is a
separate calculation path, not another name for `megan21`.

Use the [standalone example](../examples/cece_config_megan3.yaml) and the
[method-selection instructions](#standalone-driver-setup) below. Both schemes
are supported through their C++ implementations.

---

## MEGAN (Single-Species Isoprene)

The original scheme computes isoprene emissions using activity factors for LAI, temperature (light-dependent and light-independent pathways), PAR (via PCEEA), leaf age, soil moisture, and CO₂ inhibition. Ported from HEMCO's `hcox_megan_mod.F90`.

### Registration Names

- Native C++: `"megan"`

### Emission methods

| `megan_method` value | Description |
|---|---|
| `"megan21"` (default) | Configurable MEGAN2.1 isoprene |
| `"native"` | Compatibility alias for `"megan21"`; unchanged calculation |
| `"hemco_3_12_1"` | Source-pinned HEMCO 3.12.1 isoprene calculation; evolving histories by default, optional fixed scalar histories |

### Configuration (MEGAN2.1 mode)

```yaml
physics_schemes:
  - name: megan
    language: cpp
    options:
      megan_method: megan21   # optional; "native" remains accepted
      beta: 0.13
      ldf: 1.0
      aef: 1.0e-9
      co2_concentration: 400.0
      day_of_year: 171       # match the standalone example's 20 June 2021 date
```

### Configuration (HEMCO 3.12.1 source-conformance mode)

```yaml
physics_schemes:
  - name: megan
    language: cpp
    options:
      megan_method: hemco_3_12_1
      aef: 1.0e-9                 # or: hemco_aef_source: gridded, with an aef_isoprene import field
      hemco_normalize_lai: false  # default true; needs the pft_* import fields
      hemco_co2_inhibition: true
      hemco_co2_ppm: 390.0
      hemco_par_direct_history_wm2: 30.0
      hemco_par_diffuse_history_wm2: 48.0
      hemco_temperature_history_k: 288.15
      history_mode: evolving   # default; the CECE clock supplies dt and day of year
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `megan_method` | string | `"megan21"` | `"megan21"`, its `"native"` alias, or `"hemco_3_12_1"` |
| `aef` | double | required unless `hemco_aef_source: gridded` | Scalar effective AEF [kg m⁻² s⁻¹] (`aef_isop` alias) |
| `hemco_aef_source` | string | `scalar` if `aef` is set, else required | `scalar` uses `aef`; `gridded` reads the `aef_isoprene` import field [kg m⁻² s⁻¹] (HEMCO `MEGAN_AEF_ISOP`) and rejects a scalar `aef` |
| `hemco_isoprene_scaling` | double | `1.0` | HEMCO `Isoprene scaling`, a nonnegative multiplier on the isoprene flux (applies to scalar and gridded AEF) |
| `hemco_normalize_lai` | bool | `true` | HEMCO `Normalize LAI`: divide current and previous LAI by the sum of the 15 vegetated PFT fractions and cap at 6 m² m⁻². Requires the `pft_*` import fields; set `false` when supplying already-normalized LAI |
| `species_name` | string | `"isoprene"` | This source-conformance path accepts only `"isoprene"` or HEMCO name `"ISOP"` |
| `hemco_co2_inhibition` | bool | required | Explicitly enable or disable the HEMCO CO₂ inhibition option |
| `hemco_co2_ppm` | double | required if enabled | CO₂ concentration [ppm] in the finite range 150-1250 (`hemco_3_12_1` only) |
| `history_mode` | string | `"evolving"` | `"evolving"` initializes HEMCO-like per-cell histories and updates after emissions; `"fixed"` uses configured values |
| `hemco_par_direct_history_wm2` | double | `30.0` | Initial direct PAR history [W m⁻²] in evolving mode, fixed direct history otherwise |
| `hemco_par_diffuse_history_wm2` | double | `48.0` | Initial diffuse PAR history [W m⁻²] in evolving mode, fixed diffuse history otherwise |
| `hemco_temperature_history_k` | double | `288.15` | Initial/fixed historical temperature [K] |
| `hemco_solar_angle_mode` | string | `clock` | `clock` computes HEMCO's `SOLAR_ANGLE` PAR-response sine from the clock's UTC day of year, local hour and grid latitude; `solar_cosine` then only gates day/night (`> 0`). `supplied` uses `solar_cosine` directly as the sine. |
| `hemco_day_of_year` | integer | required in fixed mode | Fixed-history day of year. Ignored (with a warning) in evolving mode, where the CECE clock supplies it |

`co2_concentration` remains a compatibility alias for `hemco_co2_ppm`, and
`aef_isop` is accepted when `aef` is absent. Source mode requires a finite,
nonnegative effective AEF. The obsolete `hemco_par_avg_umol` and
`hemco_t_avg_15_k` options are rejected; use the separate direct/diffuse PAR
histories and `hemco_temperature_history_k` instead.

### Import Fields

| Field Name | Units | Description |
| --- | --- | --- |
| `temperature` | K | Surface air temperature |
| `leaf_area_index` | m²/m² | Effective current LAI |
| `leaf_area_index_prev` | m²/m² | Required exact effective previous-day `PMISOLAI` only in fixed-history mode; evolving mode initializes it from current LAI |
| `par_direct` | W/m² | Direct PAR |
| `par_diffuse` | W/m² | Diffuse PAR |
| `solar_cosine` | — | Finite value in [-1, 1]. With `hemco_solar_angle_mode: clock` it is HEMCO's `SUNCOS` day/night gate (nonpositive is night); with `supplied` it is the effective sine of solar elevation |
| `soil_moisture_root` | fraction | Root-zone soil moisture (optional) |
| `aef_isoprene` | kg m⁻² s⁻¹ | Gridded isoprene AEF; required with `hemco_aef_source: gridded` (finite, nonnegative) |
| `pft_ndlf_evgn_tmpt_tree`, `pft_ndlf_evgn_borl_tree`, `pft_ndlf_decd_borl_tree`, `pft_bdlf_evgn_trop_tree`, `pft_bdlf_evgn_tmpt_tree`, `pft_bdlf_decd_trop_tree`, `pft_bdlf_decd_tmpt_tree`, `pft_bdlf_decd_borl_tree`, `pft_bdlf_evgn_shrb`, `pft_bdlf_decd_tmpt_shrb`, `pft_bdlf_decd_borl_shrb`, `pft_c3_arct_grss`, `pft_c3_narc_grss`, `pft_c4_grss`, `pft_crop` | fraction | CLM4 vegetated PFT fractions (HEMCO `CLM4_PFT_*`, excluding bare ground); all 15 are required when `hemco_normalize_lai: true` |

Field names can be remapped with the scheme's `input_mapping`. As in HEMCO, a
cell is zeroed when the unnormalized `LAI * AEF` is not positive, and a cell
whose PFT sum is zero keeps its raw LAI. The 12-hour LAI history stores raw
LAI; normalization is applied to the current and history values each step.

In `clock` mode the scheme needs the clock's day of year and the grid latitude
and longitude (supplied by the driver). Local hour uses the `local_time`
service offsets when it is enabled, and otherwise HEMCO's 15-degree longitude
bins. The scheme throws if the coordinates or day of year are missing.

All required source-conformance fields and the output must have the same
horizontal extents and one level. `leaf_area_index_prev` is required only in
`history_mode: fixed`; evolving mode initializes previous-day LAI from current
LAI and advances it as a 12-hour running average.
The scheme fails before launching its kernel if those shapes differ.

For `hemco_3_12_1`, a scalar AEF must already include any upstream PFT/AEF
processing; a gridded AEF is used as given and is then multiplied by
`hemco_isoprene_scaling`. With `hemco_normalize_lai: false`, supply both LAI
fields already normalized. With a scalar AEF, the configured cold-start/fixed histories are spatially
uniform. The default evolving mode stores per-cell state in double precision and updates
the 5-day temperature/direct-PAR/diffuse-PAR means and 12-hour temperature/LAI
means after each emission calculation. The 12-hour temperature mean mirrors
HEMCO's `T_PREVDAY`, which HEMCO 3.12.1 passes to its light-dependent
temperature function but documents as unused; no emission term reads it here
either. A cold start uses 288.15 K, 30 W/m² direct PAR, 48 W/m² diffuse PAR, and
current LAI; restart I/O is not yet implemented. PAR remains configured in W/m²
and is converted internally to µmol/m²/s using the HEMCO factor 4.766.

Evolving mode requires the scheduler clock: the running-mean weight is
`dt / tau`, where `dt` is the scheme's refresh interval (`refresh_interval_seconds`,
or the base `timestep_seconds`), matching HEMCO's emission time step, and the day
of year comes from the clock. Calling the scheme without that context fails with
an error; use `history_mode: fixed` for context-free calls. If `dt` exceeds the
tau of a running mean (for example a scheme interval above 12 hours for the LAI
history) the weight is clamped to 1, so that history equals the current sample,
and a warning is logged once.

### Effective-history settings for the configurable C++ schemes

The following scalar options apply to `megan_method: megan21` (including its
`native` alias) and `name: megan3`, both with `language: cpp`:

| Option | Default | Meaning |
|---|---|---|
| `temperature_history_k` | 297 | Finite positive effective historical temperature, K |
| `par_history_wm2` | 78 | Finite nonnegative effective total historical PAR, W/m² |
| `days_between_lai` | 30 | Finite positive interval for leaf-age calculation, days |
| `day_of_year` | 180 | Integer in 1-366 |
| `leaf_age_uses_temperature_history` | false | Use history temperature for leaf age instead of current temperature |

For `megan21`, these settings remain supplied scalar histories. For `megan3`
they are the `history_mode: fixed` settings. In the default evolving mode,
`temperature_history_k` and `par_history_wm2` only seed the cold start (the PAR
value seeds the diffuse history, direct starts at zero; without them the HEMCO
values 288.15 K, 30 and 48 W/m² are used), while `day_of_year`,
`days_between_lai`, and `leaf_age_uses_temperature_history` are ignored with a
warning: the clock supplies the day, leaf age uses a one-day interval and the
history temperature. `megan3` evolving mode requires the scheduler clock like
the HEMCO method above.
Changing the method name from `native` to `megan21` does not alter defaults or
calculations.

### Export Fields

| Field Name | Units | Description |
| --- | --- | --- |
| `isoprene_emissions` | kg/m²/s | Isoprene emission flux |

---

## MEGAN3 (Multi-Species, Multi-Class)

The C++ MEGAN3-class scheme computes 19 class totals using bulk activity factors and converts them to mechanism-specific output species via a configurable speciation engine. Its separate canopy helpers do not currently affect those runtime totals.

### Registration Names

- Native C++: `"megan3"`

### 19 Emission Classes

| Class | Description |
| --- | --- |
| ISOP | Isoprene |
| MBO | 2-methyl-3-buten-2-ol |
| MT_PINE | Monoterpenes (α-pinene type) |
| MT_ACYC | Monoterpenes (acyclic, e.g., myrcene) |
| MT_CAMP | Monoterpenes (camphene type) |
| MT_SABI | Monoterpenes (sabinene type) |
| MT_AROM | Monoterpenes (aromatic, e.g., p-cymene) |
| NO | Nitric oxide (from soil, via export state) |
| SQT_HR | Sesquiterpenes (high reactivity) |
| SQT_LR | Sesquiterpenes (low reactivity) |
| MEOH | Methanol |
| ACTO | Acetone |
| ETOH | Ethanol |
| ACID | Organic acids |
| LVOC | Low-volatility organic compounds |
| OXPROD | Oxidation products |
| STRESS | Stress-induced emissions |
| OTHER | Other VOCs |
| CO | Carbon monoxide |

### Gamma Factors

For each vegetation emission class (soil NO is handled separately), the scheme uses:

- **γ_T_LI** — Light-independent temperature response (exponential β formulation)
- **γ_T_LD** — Light-dependent temperature response (Guenther et al. 2012)
- **γ_PAR** — PAR response via PCEEA algorithm
- **γ_LAI** — Leaf area index correction
- **γ_age** — Leaf age (new/growing/mature/old fractions)
- **γ_SM** — Currently neutral (1), including when `soil_moisture_root` is supplied
- **γ_CO₂** — CO₂ inhibition (Possell or Wilkinson)
- **γ_stress** — Wind and temperature stress (optional)

The `enable_aq_stress` setting is parsed but is not applied by
`Megan3Scheme::Run`. These limitations describe the current bulk C++ calculation;
the separate canopy helpers below are also not part of that runtime calculation.

Combined via LDF partitioning:

```
emission[class] = NORM_FAC × AEF × γ_LAI × γ_age × γ_SM × γ_CO₂ × [(1-LDF)×γ_T_LI + LDF×γ_PAR×γ_T_LD] × γ_stress
```

### Canopy Model (MEGCANOPY)

A separate, tested 5-layer Gaussian quadrature canopy helper provides the following calculations. These are not yet integrated into the C++ `Megan3Scheme::Run` emission calculation:
- Beer-Lambert PAR extinction through the canopy
- Sunlit/shaded leaf fractions at each layer
- Leaf temperature via energy balance
- Canopy-integrated emission activity factor

All light-dependent factors are zero when solar cosine ≤ 0 (nighttime).

### Speciation

After computing 19 class totals, the speciation engine converts them to mechanism-specific output species using scale factors from a YAML MAP file:

```
output[species] = (Σ class_total[c] × scale_factor[c→s]) × MW[s]
```

The vegetation-class totals and their AEF inputs are amount fluxes in kmol
class m⁻² s⁻¹. The loader multiplies each SPC molecular weight in kg mol⁻¹ by
1000; the resulting numerical value is also the molecular weight in kg
kmol⁻¹ used by the calculation. Mechanism-species outputs are therefore mass
fluxes in kg m⁻² s⁻¹. Divide a mass-basis AEF by its class molecular weight
before supplying it to MEGAN3.

See [Speciation Configuration](#speciation-configuration) below for the YAML format.

### Soil NO Handling

The NO emission class reads `soil_nox_emissions` from the export state
(produced by BDSNP or another soil NO scheme). If the field is missing, a
warning is logged and NO is set to zero. Canonical `soil_no_method: bdsnp`
requires the complete effective-input contract documented in
[BDSNP Soil NO Emissions](soil_nox.md), including both 24-layer fields. The
compact example below selects `yl95` so it does not imply those effective
inputs are generated by MEGAN3.

`soil_nox_emissions` has a mass-flux contract of kg NO m⁻² s⁻¹. MEGAN3
converts that field to an NO amount flux before the generic speciation engine
applies target-species molecular weights. This preserves the incoming mass for
the one-to-one `NO` mapping and permits other configured mappings to retain the
normal speciation behavior. Soil NO is independent of the vegetation LAI gate:
a positive soil flux remains positive when `leaf_area_index` is zero, while
LAI-dependent biogenic VOC classes remain zero there.

### Configuration

```yaml
physics_schemes:
  - name: bdsnp
    options:
      soil_no_method: yl95
  - name: megan3
    options:
      mechanism_file: data/speciation/spc_cb6.yaml
      speciation_file: data/speciation/map_cb6.yaml
      speciation_dataset: MEGAN
      co2_concentration: 415.0
      co2_method: possell
      enable_wind_stress: false
      enable_temp_stress: false
      emission_classes:
        ISOP:
          ldf: 0.9996
          ct1: 95.0
          cleo: 2.0
          beta: 0.13
          anew: 0.05
          agro: 0.6
          amat: 1.0
          aold: 0.9
          default_aef: 1.0e-9  # kmol ISOP m-2 s-1
        MT_PINE:
          ldf: 0.10
          ct1: 80.0
          cleo: 1.83
          beta: 0.10
          anew: 2.0
          agro: 1.8
          amat: 1.0
          aold: 1.05
          default_aef: 3.0e-10  # kmol MT_PINE m-2 s-1
        # ... remaining classes
      input_mapping:
        temperature: T2M
        leaf_area_index: LAI
        par_direct: PARDR
        par_diffuse: PARDF
        solar_cosine: COSZS
        soil_moisture_root: GWETROOT
        wind_speed: U10M
```

### Import Fields

| Field Name | Units | Description |
| --- | --- | --- |
| `temperature` | K | Surface air temperature |
| `leaf_area_index` | m²/m² | Current month LAI |
| `leaf_area_index_prev` | m²/m² | Previous month LAI (optional) |
| `par_direct` | W/m² | Direct PAR |
| `par_diffuse` | W/m² | Diffuse PAR |
| `solar_cosine` | — | Cosine of solar zenith angle |
| `soil_moisture_root` | fraction | Root-zone soil moisture (optional) |
| `wind_speed` | m/s | Wind speed (optional, for stress) |
| `AEF_<CLASS>` | kmol class/m²/s | Per-class gridded amount-basis AEF (optional) |

### Export Fields

Dynamic — one field per mechanism species with `MEGAN_` prefix:

| Field Name | Units | Description |
| --- | --- | --- |
| `MEGAN_ISOP` | kg/m²/s | Isoprene (mechanism species) |
| `MEGAN_TERP` | kg/m²/s | Terpenes (mechanism species) |
| `MEGAN_<NAME>` | kg/m²/s | Other mechanism species |

The current C++ MEGAN3 speciation engine writes these `MEGAN_`-prefixed names
directly. Select them in `output.fields`; an `output_mapping` entry does not
rename the speciation engine's output fields.

---

## Speciation Configuration

MEGAN3 uses two YAML files for chemical mechanism speciation:

### SPC File (Mechanism Species)

Defines the target mechanism species and their molecular weights. Uses the MICM/OpenAtmos format:

```yaml
name: CB6_AE7
species:
  - name: ISOP
    molecular weight [kg mol-1]: 0.06812
  - name: TERP
    molecular weight [kg mol-1]: 0.13623
  - name: PAR
    molecular weight [kg mol-1]: 0.01443
  # ... up to 36 species for CB6
```

### MAP File (Speciation Mappings)

Defines how the 19 MEGAN emission classes map to mechanism species with per-class scale factors. Uses a dataset-oriented format that supports multiple emission sources:

```yaml
mechanism: CB6_AE7
datasets:
  MEGAN:
    ISOP:
      ISOP: 1.0
    TERP:
      MT_PINE: 0.5
      MT_ACYC: 0.3
      MT_CAMP: 0.1
      MT_SABI: 0.05
      MT_AROM: 0.05
    MEOH:
      MEOH: 1.0
    SESQ:
      SQT_HR: 0.7
      SQT_LR: 0.3
    "NO":
      "NO": 1.0
    CO:
      CO: 1.0
```

Each entry under a mechanism species name lists the contributing emission classes and their fractional scale factors. The speciation engine computes:

```
output[TERP] = (class_total[MT_PINE] × 0.5 + class_total[MT_ACYC] × 0.3 + ...) × MW[TERP]
```

### Supported Mechanisms

CECE ships with speciation files for:

| Mechanism | SPC File | MAP File |
| --- | --- | --- |
| CB6_AE7 | `data/speciation/spc_cb6.yaml` | `data/speciation/map_cb6.yaml` |
| RACM2 | `data/speciation/spc_racm2.yaml` | `data/speciation/map_racm2.yaml` |
| SAPRC07 | `data/speciation/spc_saprc07.yaml` | `data/speciation/map_saprc07.yaml` |
| CRACMM2 | `data/speciation/spc_cracmm.yaml` | `data/speciation/map_cracmm.yaml` |

To switch mechanisms at runtime, change `mechanism_file` and `speciation_file` in the YAML config — no recompilation needed.

### Adding a New Mechanism

1. Create an SPC file with species names and molecular weights (kg/mol)
2. Create a MAP file with a `MEGAN` dataset section mapping emission classes to your mechanism species
3. Point `mechanism_file` and `speciation_file` to your new files

---

## Scheme Ordering

BDSNP must run before MEGAN3 so that `soil_nox_emissions` is available. This
ordering applies to both canonical `bdsnp` and the compact `yl95` fallback:

```yaml
physics_schemes:
  - name: bdsnp        # Runs first, writes soil_nox_emissions
  - name: megan3       # Runs second, reads soil_nox_emissions for NO class
```

## Standalone-driver setup

The single [MEGAN example](../examples/cece_config_megan3.yaml) provides the
driver, data streams, MEGAN3 class settings, and output configuration. Replace
every `/path/to` value with your files, and run from a directory where the SPC
and MAP paths resolve. Its input data and gridspec are not bundled. This is a
configuration template, not an executable HEMCO reference experiment.

The example uses exact-grid passthrough and a one-level global GEOS 4°×5° grid.
All streams and the gridspec must have identical ordered coordinates:
longitude centers `-180, -175, ..., 175` and latitude centers
`-89, -86, -82, ..., 86, 89`. A uniform `lat_min/lat_max/ny` grid does not
reproduce the polar half cells. Use passthrough only for matching meshes;
prepare compatible input fields or configure supported regridding otherwise.

To run MEGAN2.1 instead, replace the example's entire `physics_schemes` list
with the [MEGAN2.1 configuration](#configuration-megan21-mode) above. To run
the source-pinned method, use the
[HEMCO configuration](#configuration-hemco-3121-source-conformance-mode).
Both alternatives use the same canonical stream field names. For either
single-species method, replace `output.fields` with:

<!-- cece-validate: context output -->
```yaml
fields:
  - name: isoprene_emissions
    attributes:
      units: "kg m-2 s-1"
      long_name: "CECE isoprene emission flux"
```

Changing the method does not by itself make meteorological fields HEMCO
effective inputs. In particular, for an executed-reference comparison, replace
the example's `COSZS` data with the matched MEGAN PAR-response sine combined
with HEMCO's outer day/night gate, as described below. A meteorological zenith
cosine with the same variable name is not necessarily equivalent.

Use separate output directories when comparing configurations. Schedule the
two single-species methods in separate runs: both register as `megan`, and
the current dispatcher selects the first configured instance with a requested
name. `megan3` and `megan` can coexist because their registration names differ.
Do not use `megan21` as a scheme name; it is a `megan_method` value.

The example's class coefficients and AEFs are illustrative, not a global
vegetation inventory. For a matched isoprene AEF comparison, the single-species
`aef: 1.0e-9` kg isoprene m⁻² s⁻¹ corresponds to MEGAN3 ISOP
`default_aef: 1.4679976512037582e-11` kmol ISOP m⁻² s⁻¹ with the configured
68.12 kg kmol⁻¹ molecular weight. Do not use the same numerical AEF for both
mass and amount bases. Other MEGAN3 classes require their own appropriate AEFs.

The standalone executable accepts the configuration path as its argument:

```console
./build/cece_standalone_driver examples/cece_config_megan3.yaml
```

Adjust the executable path to your build. Soil NO is zero in the standalone
MEGAN3 example unless a preceding soil-NO scheme supplies the export field.
For coupling, add the necessary soil-NO input streams and schedule BDSNP
before MEGAN3 as described above; MEGAN3 does not synthesize BDSNP inputs.

## HEMCO source settings and independent-reference validation

The source-pinned method follows
[HEMCO 3.12.1 `hcox_megan_mod.F90`](https://github.com/geoschem/HEMCO/blob/07da3c29fd85abc3824cb6288578b0b68c2395a3/src/Extensions/hcox_megan_mod.F90).
Its source SHA-256 is
`a298e4003210c7dba86c53cdd37f85a868dcb3a89b3de56ab175257e04614f31`.
The device-callable implementation is in
`include/cece/physics/cece_megan.hpp`; the independent Python
transcription generates the scalar vectors in `tests/data/hemco_megan/`.

| Quantity | Source-pinned treatment |
|---|---|
| ISOP parameters | `LDF=1`, `CT1=95`, `CEO=2`, leaf-age weights `0.05/0.60/1.00/0.90` |
| Instantaneous PAR | Separate direct and diffuse W/m² fields, each converted by `4.766` |
| PTOA | `3000 + 99 cos(2π(DOY-10)/365)` µmol m⁻² s⁻¹ |
| PAR history | Separate `PARDR_DAVG` and `PARDF_DAVG`, W/m² |
| Temperature history | `T_DAVG`, K, used for light-dependent temperature response and leaf age |
| Previous LAI | Effective previous-day `PMISOLAI` after upstream storage/preprocessing |
| LAI interval | One day |
| CO₂ inhibition | Explicit switch; Possell and Hewitt (2011) when enabled |
| Normalization | Full `CALC_NORM_FAC` expression (`0.9899364002107353`) |

These equations do not imply identical behavior to the configurable MEGAN2.1
or MEGAN3 paths: normalization, LDF settings, leaf-age settings, and
low-solar-elevation treatment can differ. Compare those paths both at their
defaults and with matched effective histories. Differences alone establish
neither a bug nor a scientific improvement.

To compare with an independently executed HEMCO reference, match what reaches
the emission calculation rather than merely matching NetCDF variable names:

1. Preserve the pinned source, build settings, executable hash, configuration,
   input hashes, and output diagnostics for the reference run.
2. Account for HEMCO's `REAL(sp)` input storage, promoting those effective
   values to CECE's double precision. Match effective AEF and both LAI fields,
   including PFT normalization and isoprene scaling where enabled. For a
   no-restart case with LAI normalization disabled, initial previous-day LAI
   equals current LAI; an arbitrary input `LAI_PREV` is not interchangeable.
3. Match the actual history state used for emissions. CECE stores evolving
   histories in double precision; HEMCO stores them in `REAL(sp)`. For smooth
   quantities the difference is negligible (a test run of the 5-day temperature
   mean over 60 hourly steps differed by about 2e-4 K), but MEGAN's leaf-age
   function branches on exact equality of current and previous LAI. With a
   constant LAI that is not exactly representable in single precision, a double
   history stays exactly equal (leaf-age factor 0.95) while a single-precision
   history never does (factor of about 1.0), a roughly 5% emission difference in
   such cells. Check for this when comparing against an executed HEMCO run.
   At cold start, previous LAI is current LAI. Temperature history is 288.15 K,
   direct PAR is 30 W/m², and diffuse PAR is 48 W/m². A restart written after
   emissions is not automatically the state used for the preceding calculation.
   Choose the correct date/DOY and CO₂ settings;
   for example, 20 June 2021 is DOY 171, and 390 ppm is a case choice.
4. Match both solar calculations: HEMCO's `ExtState%SUNCOS` gates day/night,
   while MEGAN `SOLAR_ANGLE` supplies the PAR-response sine using local time
   from `HcoClock_GetLocal`. Without a TIMEZONES file, that clock uses
   15-degree longitude bins. CECE's default `hemco_solar_angle_mode: clock`
   reproduces this: supply HEMCO's `SUNCOS` as `solar_cosine` and CECE derives
   the PAR sine itself. For reference comparisons that start from an
   externally computed PAR sine, use `supplied`. Preserve both upstream fields
   and local-time settings for each reference case; `HCO_SUNCOS` alone is not
   the PAR-response sine.
5. Compare ordered coordinates, masks, mass units, and single-/multi-rank
   output before interpreting flux differences. HEMCO mean diagnostics use
   `REAL(sp)` storage and accumulation. For a single 3600-second interval,
   report raw CECE-double differences separately from the explicit HEMCO
   diagnostic cast/accumulate/divide projection. Do not call raw double values
   bitwise identical to single-precision diagnostics.

Source-conformance tests do not launch HEMCO. An independently executed,
controlled reference comparison is a separate validation level; neither
establishes real-meteorology, restart parity, or complete upstream MEGAN3
parity. Future work includes multi-compound `CALC_AEF` generation from PFTs, a
timezone-grid local time for the standard solar angle, multilayer canopy integration, and multi-timestep/coupled
validation. Validation reports, plots,
and run-specific inputs/outputs belong with the external experiment evidence.

### Regression tests

With tests enabled in the configured build:

```console
cmake --build build --target test_hemco_megan_runtime test_hemco_megan_global test_megan3
ctest --test-dir build -R 'HEMCO3121|HemcoMeganGlobal|Megan3' --output-on-failure
python tests/test_megan_global_parity.py
```

The historical Python filename is retained, but its synthetic global tests
establish source invariants, not independent HEMCO execution. The scalar
oracle contains 16 cases, including mixed-precision previous-LAI coverage.
To regenerate it deliberately:

```console
python scripts/generate_hemco_megan_oracle.py \
  > tests/data/hemco_megan/hemco_3_12_1_megan_reference.csv
```

Review the regenerated diff and run the regression tests before accepting it.
