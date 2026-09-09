# Material data inventory

Snapshot: **2026-09-09 UTC**, refreshed after the metal/water curves, silicon tensor, dry-air model, 316 reference-state and copper/aluminum resistivity additions, current shared working tree on main, including uncommitted material additions. This is an inventory of populated source records, not a claim that every bundle compiles, every property is measured, or every material can run a complete simulation.

Focused consumer evidence now includes actual source compilation, storage, material resolution and steady conduction for NIST 304/316/6061 and IAPWS liquid water. The water check covers 20–25 °C and 60–65 °C at exactly 0.1 MPa, retains liquid-phase receipts and refuses unsupported states. These manufactured temperature/flux and energy-balance checks exercise the declared interpolants; they do not establish experimental accuracy, flowing-water behavior, transient heat or phase changes. The executable cases are in [matdb_pack_cli.rs](../xtask/tests/matdb_pack_cli.rs).

The sourced catalog lives in **[data/matdb/seed-v1](../data/matdb/seed-v1/README.md)** as static text files. It currently occupies **0.794138 MB** (0.757349 MiB), including its README and license notice. The TSV data and manifests alone occupy **0.672951 MB**. There are **169 source bundles**, including **151 bulk-material bundles**.

A bundle identifies a source and condition, not a unique chemical material. Copper, lead, stainless steel, and other materials have multiple bundles for different sources, phases, grades, processing, or measurements. Conversely, one bundle can contain several specimens or conditions. The directory therefore does **not** establish an exact deduplicated count of distinct materials, nor 150 complete material cards. The tables below use the exact, auditable source-bundle count.

## Where the values live and how code reads them

| Location | Role |
| --- | --- |
| [data/matdb/seed-v1/](../data/matdb/seed-v1/README.md) | Human-readable material source records, units, conditions, observations, and licenses |
| `data/matdb/seed-v1/<bundle>/manifest.tsv` | Pack ID, source-file path/profile, citation, redistribution decision, and optional property/axis mappings |
| `properties.tsv` / `interface.tsv` | Scalar values and curves, source observations, validity conditions, and uncertainty; interfaces also identify the ordered surfaces and system |
| `species.tsv` | Gas identity, molar mass, phase/EOS association, reference pressure, and reference convention |
| `contact.tsv` | Authored contact-law parameters with pair, geometry, identification and validity context |
| [xtask/src/matdb_pack.rs](../xtask/src/matdb_pack.rs) | Offline compiler: validates supported source profiles, normalizes units, and creates deterministic binary packs |
| [crates/fs-matdb/src/pack.rs](../crates/fs-matdb/src/pack.rs) | Normalized material-property pack format and loading |
| [crates/fs-matdb-store/src/lib.rs](../crates/fs-matdb-store/src/lib.rs) | FrankenSQLite storage, transactional ingestion, indexing, and discovery |
| [crates/fs-cli/src/discover.rs](../crates/fs-cli/src/discover.rs) | Reads supplied binary packs and ingests them into an in-memory store for the discovery command |
| [crates/fs-material/src/state_point.rs](../crates/fs-material/src/state_point.rs) | Resolves selected material-card properties for physical consumers |
| [crates/fs-matdb/src/contact_pack.rs](../crates/fs-matdb/src/contact_pack.rs) | Separate contact-card loader; reads `contact.tsv` directly |

The regular material path is source TSV → offline compiler → normalized pack → store/direct typed query → material resolution → consumer. Binary packs are generated artifacts; the seed README says they are deliberately not committed. There is no single prebuilt, comprehensive material database file in this directory. The runtime does not automatically treat every source directory as a complete, loaded specimen.

The silicon reference now follows that path through all 36 stiffness components
and density into the real oriented tetrahedral operator. Its focused numerical
test checks forces, energy and mass at two orientations against an independent
cubic-crystal formula. This is a sourced engineering reference at 25 °C with
explicit cross-source/unknown-condition limits, not a qualified wafer model.

The dry-air source now resolves through the real compiler and reopened store
into the existing gas model and cylinder acoustic-loss calculation. Four
temperature/pressure states match independent source-equation/property and
loss references; a synthetic reference-viscosity change alters the result.
The application domain is explicitly 273.15–313.15 K, 80–110 kPa, dry,
USSA-1976 composition. These are five model parameters, not five new measured
thermophysical properties, and the result does not qualify humidity transport.

The complete 2024-T3 and 316 reference profiles now pass the same actual plate
consumer through compilation, persistent storage/reopen, discovery and typed
resolution. All six selected properties and their receipts reach plate mass,
stiffness and thermal damping; 128-step pressure traces are deterministic and
change when thermal damping is removed. The test checks mass and Zener damping
against independent equations, and modal frequency against the existing 20%
coarse-mesh continuum band (observed difference about 2.3%). The profiles use
different explicit temperatures, 300 K and 293.15 K. The 316 data are a
cross-source engineering approximation with pressure and product form unknown;
this evidence does not establish a controlled same-temperature substitution,
measured acoustic accuracy, or a finite heating trajectory. The eight-case
`scripts/e2e/material_sources.sh --run` remote run passed on 2026-09-09 UTC.

The copper and EC-H19 aluminum resistivity curves also pass compilation,
persistent storage/reopen, discovery and the source-resolved uniform conductor
adapter into the existing circuit DAE. Eight source knots and two interior
interpolation states check resistance, voltage, Joule energy, supplied energy,
geometry scaling and deterministic replay. At 20 °C, a 1 m conductor with
1 mm² cross-section has resistance 0.017241 ohm for copper and 0.028264 ohm for
aluminum. At 2 A for 1 s, the circuit dissipates 0.068964 J and 0.113056 J,
respectively. The nine-case remote runner passed on 2026-09-09 UTC (receipt
`frankensim-material-sources.UKHQHK`, 9 passed, zero failed/ignored).
The handbook approximations and respective 10–30 °C / 0–30 °C coverage remain
explicit. This proves the fixed-state DC path, not temperature feedback,
ampacity, contact resistance, AC behavior or experimental qualification.

There are also numerical presets in code. For example, [ThermoelasticZener handbook constructors in visco.rs](../crates/fs-material/src/visco.rs) embed aluminum/structural-steel constants, and [rc_section in fs-solid](../crates/fs-solid/src/fiber.rs) constructs an example steel law with fixed parameters. These are outside this source-data census. The dry-air adapter reuses the existing gas equations and explicit conductivity-model choice in [gas.rs](../crates/fs-material/src/gas.rs); adding its parameter pack does not replace every existing ambient preset automatically. This inventory is not an exhaustive audit of constants, examples or test fixtures throughout the repository.

## Size and populated records

Sizes are logical file-content bytes, not filesystem allocation, Git object size, remote build cache size, or compressed size. MB means 1,000,000 bytes; MiB means 1,048,576 bytes.

| Measured item | Bytes | MB | MiB |
| --- | ---: | ---: | ---: |
| All files in data/matdb/seed-v1 | 794,138 | 0.794138 | 0.757349 |
| All TSV files, including manifests | 672,951 | 0.672951 | 0.641776 |
| Manifests | 158,140 | 0.158140 | 0.150814 |
| Source TSV plus axis-convention TSV | 514,811 | 0.514811 | 0.490962 |
| README and license notice | 121,187 | 0.121187 | 0.115573 |

There are **341 files**: 339 TSV files and 2 Markdown files. The TSV set contains 169 manifests, 169 referenced source files and one [instrument-axis-convention.tsv](../data/matdb/seed-v1/instrument-axis-convention.tsv). No PDF, compiled pack, or SQLite file is present under this directory in this snapshot. Downloaded reference PDFs in session scratch storage are excluded.

| Bundle category | Bundles | Distinct populated property names, summed per bundle | Scalar claims | Curve claims | Curve knots |
| --- | ---: | ---: | ---: | ---: | ---: |
| Bulk material/condition | 151 | 828 | 864 | 59 | 614 |
| Ordered interface system | 8 | 23 | 29 | 0 | 0 |
| Gas species metadata | 7 | 0 | 0 | 0 | 0 |
| Authored contact-law card | 3 | 0 | 0 | 0 | 0 |

The bulk and interface sources contain **952 property claims**: 893 scalars and 59 curves with 614 retained knots. That is 1,507 scalar values or curve ordinate entries, before counting temperatures, validity endpoints, uncertainty fields, species metadata, or contact parameters. Five curves replaced ten isolated 6061-T6/OFHC endpoint claims, liquid water added five curves with 75 knots, and two stainless instantaneous-expansion curves add 22 knots derived from the published relative-length fits. The silicon tensor adds 36 explicitly addressed entries derived from only three cubic constants, plus a separately sourced density. Dry air and dilute water vapor each add five parameters for the existing gas model, not separate measured thermodynamic/transport output properties. The 316 engineering reference adds six scalar inputs at 20 °C, four evaluated from existing NIST equations and two supplier facts; it is not six new measurements. Copper and EC-H19 aluminum add two electrical volume-resistivity curves with eight knots, retaining the handbook reference approximations. These entries include measurements, source fits, handbook values and authored/model inputs; they are not all independent experimental measurements.

Across the 151 bulk bundles, distinct populated property names range from **1 to 18**, with **median 4** and **mean 5.48**. The sum is 828 populated bundle/property-name pairs. Across bulk and interface bundles together there are 304 distinct property identifier spellings after each manifest's explicit mapping. This last count does not merge legacy spelling aliases or prove semantic equivalence. The silicon tensor's 36 coordinate-bearing `stiffness` keys count as one property name here, alongside density; they are not a scalar isotropic stiffness or 36 independent measurements.

A **metric** below means a distinct populated property name in a bundle. It does not count a manifest declaration without a value, an observation, a validity axis, or an uncertainty row as another material property. Several condition-specific scalar claims for one name still count as one metric; a curve counts as one metric regardless of its knot count. Units shown are the source-record units, before compiler normalization. Byte size includes every file within that bundle's directory. Species fields and contact-law parameters are listed separately because they are not scalar/curve property claims.

## Bulk material bundles

| Source bundle | Metrics | Scalar claims | Curve claims | Knots | KiB | Populated property names [source units] |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| [water-vapor-sutherland-ambient](../data/matdb/seed-v1/water-vapor-sutherland-ambient/manifest.tsv) | 5 | 5 | 0 | 0 | 4.31 | `molar_mass` [kg/mol]; `heat_capacity_ratio` [1]; `sutherland_reference_viscosity` [Pa*s]; `sutherland_reference_temperature` [K absolute]; `sutherland_temperature` [K interval]; dilute-vapor model parameters, constant Cp and extrapolated steam transport, not pure-phase validity |
| [copper-annealed-iacs-nbs-hb100](../data/matdb/seed-v1/copper-annealed-iacs-nbs-hb100/manifest.tsv) | 1 | 0 | 1 | 3 | 2.28 | `electrical_resistivity` [Ohm*m]; standard annealed 100% IACS copper, source-derived linear volume-resistivity reference at 10–30 °C |
| [aluminum-ec-h19-nbs-hb109](../data/matdb/seed-v1/aluminum-ec-h19-nbs-hb109/manifest.tsv) | 1 | 0 | 1 | 5 | 2.28 | `electrical_resistivity` [Ohm*m]; EC-H19 wire, five source table entries with linear interpolation at 0–30 °C |
| [stainless-316-20c-engineering-reference](../data/matdb/seed-v1/stainless-316-20c-engineering-reference/manifest.tsv) | 6 | 6 | 0 | 0 | 6.13 | `density` [kg/m3]; `young_modulus` [Pa]; `poisson_ratio` [1]; `specific_heat_capacity` [J/kg/K]; `thermal_conductivity` [W/m/K]; `linear_thermal_expansion_coefficient` [K-1]; explicit cross-source engineering approximation at exactly 20 °C, pressure unknown |
| [air-dry-ussa1976](../data/matdb/seed-v1/air-dry-ussa1976/manifest.tsv) | 5 | 5 | 0 | 0 | 4.80 | `molar_mass` [kg/mol]; `heat_capacity_ratio` [1]; `sutherland_reference_viscosity` [Pa*s]; `sutherland_reference_temperature` [K absolute]; `sutherland_temperature` [K interval]; bounded dry-air model parameters, with explicit conductivity-model choice |
| [silicon-cubic-25c-nasa-rp1057](../data/matdb/seed-v1/silicon-cubic-25c-nasa-rp1057/manifest.tsv) | 2 | 37 | 0 | 0 | 14.76 | `density` [kg/m3]; `stiffness` [MPa], 36 explicit engineering-Voigt components at 25 °C derived from three cubic constants; cross-source engineering reference, not a qualified wafer |
| [aisi-1045-cold-drawn](../data/matdb/seed-v1/aisi-1045-cold-drawn/manifest.tsv) | 3 | 3 | 0 | 0 | 2.19 | `tensile_elongation_50mm` [%]; `ultimate_tensile_strength` [MPa]; `yield_strength` [MPa] |
| [aisi-4140-rc33](../data/matdb/seed-v1/aisi-4140-rc33/manifest.tsv) | 7 | 14 | 0 | 0 | 5.50 | `charpy_v_notch_impact_energy` [J]; `double_shear_ultimate_strength` [GPa]; `double_shear_yield_strength` [GPa]; `tensile_elongation_2in` [%]; `tensile_reduction_of_area` [%]; `ultimate_tensile_strength` [GPa]; `yield_strength_0p2_offset` [GPa] |
| [aisi-52100-cvm-hot-hardness](../data/matdb/seed-v1/aisi-52100-cvm-hot-hardness/manifest.tsv) | 8 | 15 | 0 | 0 | 8.59 | `carbon_mass_fraction` [%]; `chromium_mass_fraction` [%]; `manganese_mass_fraction` [%]; `phosphorus_mass_fraction` [%]; `retained_austenite_volume_fraction` [%]; `rockwell_c_scale_reading` [1]; `silicon_mass_fraction` [%]; `sulfur_mass_fraction` [%] |
| [aisi-9310-cvm-carburized](../data/matdb/seed-v1/aisi-9310-cvm-carburized/manifest.tsv) | 12 | 13 | 0 | 0 | 3.84 | `carbon_mass_fraction` [%]; `carburized_case_depth` [mm]; `case_rockwell_c_scale_reading` [1]; `chromium_mass_fraction` [%]; `copper_mass_fraction` [%]; `core_rockwell_c_scale_reading` [1]; `manganese_mass_fraction` [%]; `molybdenum_mass_fraction` [%]; `nickel_mass_fraction` [%]; `phosphorus_mass_fraction` [%]; `silicon_mass_fraction` [%]; `sulfur_mass_fraction` [%] |
| [aluminum-1100-nist-cryogenic](../data/matdb/seed-v1/aluminum-1100-nist-cryogenic/manifest.tsv) | 1 | 2 | 0 | 0 | 2.19 | `thermal-conductivity` [W/m/K] |
| [aluminum-2024-t3-nasa-tn-d6448](../data/matdb/seed-v1/aluminum-2024-t3-nasa-tn-d6448/manifest.tsv) | 7 | 7 | 0 | 0 | 2.77 | `density` [kg/m3]; `linear_thermal_expansion_coefficient` [K-1]; `loss_factor_thermoelastic_peak` [1]; `poisson_ratio` [1]; `specific_heat_capacity` [J/kg/K]; `thermal_conductivity` [W/m/K]; `young_modulus` [MPa] |
| [aluminum-2024-t3-sheet-mil-hdbk-5j](../data/matdb/seed-v1/aluminum-2024-t3-sheet-mil-hdbk-5j/manifest.tsv) | 13 | 13 | 0 | 0 | 4.03 | `compression_modulus` [GPa]; `compressive_yield_l_a_basis` [MPa]; `compressive_yield_lt_a_basis` [MPa]; `density` [kg/m3]; `elongation_lt_minimum` [%]; `poisson_ratio` [1]; `shear_modulus` [GPa]; `shear_ultimate_a_basis` [MPa]; `tensile_ultimate_l_a_basis` [MPa]; `tensile_ultimate_lt_a_basis` [MPa]; `tensile_yield_l_a_basis` [MPa]; `tensile_yield_lt_a_basis` [MPa]; `young_modulus` [GPa] |
| [aluminum-2024-t4-damping-nasa-tn-d2893](../data/matdb/seed-v1/aluminum-2024-t4-damping-nasa-tn-d2893/manifest.tsv) | 7 | 7 | 0 | 0 | 2.72 | `loss_factor_1500hz` [1]; `loss_factor_150hz` [1]; `loss_factor_15hz` [1]; `loss_factor_300hz` [1]; `loss_factor_30hz` [1]; `loss_factor_700hz` [1]; `loss_factor_70hz` [1] |
| [aluminum-6061-t6-cryogenic](../data/matdb/seed-v1/aluminum-6061-t6-cryogenic/manifest.tsv) | 3 | 0 | 3 | 72 | 5.62 | `specific-heat-capacity` [J/kg/K]; `thermal-conductivity` [W/m/K]; `young-modulus` [GPa] |
| [aluminum-7075-t6-nasa-cr-71699](../data/matdb/seed-v1/aluminum-7075-t6-nasa-cr-71699/manifest.tsv) | 1 | 0 | 1 | 6 | 2.03 | `thermal-conductivity` [W/m/K] |
| [aluminum-fusion-nasa-tp3287](../data/matdb/seed-v1/aluminum-fusion-nasa-tp3287/manifest.tsv) | 2 | 2 | 0 | 0 | 1.81 | `latent-heat-fusion` [J/kg]; `melting-point` [K] |
| [aluminum-liquid-nasa-tp3287](../data/matdb/seed-v1/aluminum-liquid-nasa-tp3287/manifest.tsv) | 2 | 0 | 2 | 10 | 2.77 | `specific-enthalpy-reference-29815k` [J/kg]; `specific-heat-capacity` [J/kg/K] |
| [aluminum-pure-nasa-cr-71699](../data/matdb/seed-v1/aluminum-pure-nasa-cr-71699/manifest.tsv) | 1 | 0 | 1 | 13 | 2.64 | `thermal-conductivity` [W/m/K] |
| [aluminum-pure-nbs-c447](../data/matdb/seed-v1/aluminum-pure-nbs-c447/manifest.tsv) | 2 | 2 | 0 | 0 | 1.44 | `density` [kg/m3]; `melting_point` [K] |
| [aluminum-solid-nasa-tp3287](../data/matdb/seed-v1/aluminum-solid-nasa-tp3287/manifest.tsv) | 2 | 0 | 2 | 20 | 2.90 | `specific-enthalpy-reference-29815k` [J/kg]; `specific-heat-capacity` [J/kg/K] |
| [ash-white-fpl-gtr282](../data/matdb/seed-v1/ash-white-fpl-gtr282/manifest.tsv) | 14 | 14 | 0 | 0 | 3.82 | `density` [kg/m3]; `er_over_el` [1]; `et_over_el` [1]; `glr_over_el` [1]; `glt_over_el` [1]; `modulus_of_elasticity_bending` [MPa]; `nu_lr` [1]; `nu_lt` [1]; `nu_rl` [1]; `nu_rt` [1]; `nu_tl` [1]; `nu_tr` [1]; `specific_gravity` [1]; `young_modulus_longitudinal` [MPa] |
| [baldcypress-fpl-gtr282](../data/matdb/seed-v1/baldcypress-fpl-gtr282/manifest.tsv) | 13 | 13 | 0 | 0 | 3.73 | `density` [kg/m3]; `er_over_el` [1]; `et_over_el` [1]; `glr_over_el` [1]; `glt_over_el` [1]; `grt_over_el` [1]; `modulus_of_elasticity_bending` [MPa]; `nu_lr` [1]; `nu_lt` [1]; `nu_rt` [1]; `nu_tr` [1]; `specific_gravity` [1]; `young_modulus_longitudinal` [MPa] |
| [balsa-fpl-gtr282](../data/matdb/seed-v1/balsa-fpl-gtr282/manifest.tsv) | 14 | 14 | 0 | 0 | 3.78 | `er_over_el` [1]; `et_over_el` [1]; `glr_over_el` [1]; `glt_over_el` [1]; `grt_over_el` [1]; `modulus_of_elasticity_bending` [MPa]; `nu_lr` [1]; `nu_lt` [1]; `nu_rl` [1]; `nu_rt` [1]; `nu_tl` [1]; `nu_tr` [1]; `specific_gravity` [1]; `young_modulus_longitudinal` [MPa] |
| [basswood-american-fpl-gtr282](../data/matdb/seed-v1/basswood-american-fpl-gtr282/manifest.tsv) | 14 | 14 | 0 | 0 | 3.83 | `density` [kg/m3]; `er_over_el` [1]; `et_over_el` [1]; `glr_over_el` [1]; `glt_over_el` [1]; `modulus_of_elasticity_bending` [MPa]; `nu_lr` [1]; `nu_lt` [1]; `nu_rl` [1]; `nu_rt` [1]; `nu_tl` [1]; `nu_tr` [1]; `specific_gravity` [1]; `young_modulus_longitudinal` [MPa] |
| [birch-yellow-fpl-gtr282](../data/matdb/seed-v1/birch-yellow-fpl-gtr282/manifest.tsv) | 15 | 15 | 0 | 0 | 3.97 | `density` [kg/m3]; `er_over_el` [1]; `et_over_el` [1]; `glr_over_el` [1]; `glt_over_el` [1]; `grt_over_el` [1]; `modulus_of_elasticity_bending` [MPa]; `nu_lr` [1]; `nu_lt` [1]; `nu_rl` [1]; `nu_rt` [1]; `nu_tl` [1]; `nu_tr` [1]; `specific_gravity` [1]; `young_modulus_longitudinal` [MPa] |
| [brass-70-30-nbs-c447-hardness](../data/matdb/seed-v1/brass-70-30-nbs-c447-hardness/manifest.tsv) | 4 | 4 | 0 | 0 | 2.67 | `hardness_brinell_annealed` [1]; `hardness_brinell_cold_rolled_11pct` [1]; `hardness_rockwell_b_cold_worked_37pct` [1]; `hardness_rockwell_b_soft_sheet` [1] |
| [brass-c26000-nist-cryogenic](../data/matdb/seed-v1/brass-c26000-nist-cryogenic/manifest.tsv) | 1 | 1 | 0 | 0 | 1.88 | `thermal-conductivity` [W/m/K] |
| [brass-cartridge-c26000-mil-hdbk-698a](../data/matdb/seed-v1/brass-cartridge-c26000-mil-hdbk-698a/manifest.tsv) | 4 | 4 | 0 | 0 | 2.17 | `density` [kg/m3]; `poisson_ratio_derived` [1]; `shear_modulus` [MPa]; `young_modulus` [MPa] |
| [brass-muntz-c28000-mil-hdbk-698a](../data/matdb/seed-v1/brass-muntz-c28000-mil-hdbk-698a/manifest.tsv) | 4 | 4 | 0 | 0 | 1.83 | `density_printed_range_low` [kg/m3]; `melting_point_liquidus` [K]; `melting_point_solidus` [K]; `thermal_conductivity` [W/m/K] |
| [brass-yellow-half-hard-damping-nasa-tn-d1467](../data/matdb/seed-v1/brass-yellow-half-hard-damping-nasa-tn-d1467/manifest.tsv) | 2 | 2 | 0 | 0 | 1.91 | `specimen_loss_factor_10ksi_500hz` [1]; `specimen_loss_factor_10ksi_50hz` [1] |
| [brick-clay-cfast-sp1041](../data/matdb/seed-v1/brick-clay-cfast-sp1041/manifest.tsv) | 5 | 5 | 0 | 0 | 5.58 | `density` [kg/m3]; `emissivity` [1]; `nominal-thickness` [m]; `specific-heat-capacity` [J/kg/K]; `thermal-conductivity` [W/m/K] |
| [brick-common-cfast-sp1041](../data/matdb/seed-v1/brick-common-cfast-sp1041/manifest.tsv) | 5 | 5 | 0 | 0 | 5.60 | `density` [kg/m3]; `emissivity` [1]; `nominal-thickness` [m]; `specific-heat-capacity` [J/kg/K]; `thermal-conductivity` [W/m/K] |
| [bronze-gunmetal-damping-rsic508](../data/matdb/seed-v1/bronze-gunmetal-damping-rsic508/manifest.tsv) | 1 | 1 | 0 | 0 | 2.11 | `specific_damping_capacity_5ksi_shear` [%] |
| [calcium-silicate-board-cfast-sp1041](../data/matdb/seed-v1/calcium-silicate-board-cfast-sp1041/manifest.tsv) | 5 | 5 | 0 | 0 | 5.71 | `density` [kg/m3]; `emissivity` [1]; `nominal-thickness` [m]; `specific-heat-capacity` [J/kg/K]; `thermal-conductivity` [W/m/K] |
| [cane-arundo-damping-scielo-mr20170795](../data/matdb/seed-v1/cane-arundo-damping-scielo-mr20170795/manifest.tsv) | 4 | 4 | 0 | 0 | 2.54 | `logarithmic_decrement_printed_range_high` [1]; `logarithmic_decrement_printed_range_low` [1]; `storage_modulus_printed_range_high` [MPa]; `storage_modulus_printed_range_low` [MPa] |
| [cane-arundo-reedblank-mdpi-ma18122759](../data/matdb/seed-v1/cane-arundo-reedblank-mdpi-ma18122759/manifest.tsv) | 4 | 4 | 0 | 0 | 2.36 | `density_printed_range_high` [kg/m3]; `density_printed_range_low` [kg/m3]; `young_modulus_longitudinal_printed_range_high` [MPa]; `young_modulus_longitudinal_printed_range_low` [MPa] |
| [cane-arundo-reedheel-mdpi-ma13204566](../data/matdb/seed-v1/cane-arundo-reedheel-mdpi-ma13204566/manifest.tsv) | 4 | 4 | 0 | 0 | 2.40 | `young_modulus_longitudinal` [MPa]; `young_modulus_longitudinal_water_soaked` [MPa]; `young_modulus_transverse` [MPa]; `young_modulus_transverse_water_soaked` [MPa] |
| [carbon-steel-cast-c033-nbs-c447](../data/matdb/seed-v1/carbon-steel-cast-c033-nbs-c447/manifest.tsv) | 2 | 2 | 0 | 0 | 1.43 | `hardness_brinell_annealed_1700f` [1]; `hardness_brinell_as_cast` [1] |
| [cedar-western-red-fpl-gtr282](../data/matdb/seed-v1/cedar-western-red-fpl-gtr282/manifest.tsv) | 13 | 13 | 0 | 0 | 3.73 | `density` [kg/m3]; `er_over_el` [1]; `et_over_el` [1]; `glr_over_el` [1]; `glt_over_el` [1]; `grt_over_el` [1]; `modulus_of_elasticity_bending` [MPa]; `nu_lr` [1]; `nu_lt` [1]; `nu_rt` [1]; `nu_tr` [1]; `specific_gravity` [1]; `young_modulus_longitudinal` [MPa] |
| [cellulose-insulation-cfast-sp1041](../data/matdb/seed-v1/cellulose-insulation-cfast-sp1041/manifest.tsv) | 5 | 5 | 0 | 0 | 5.73 | `density` [kg/m3]; `emissivity` [1]; `nominal-thickness` [m]; `specific-heat-capacity` [J/kg/K]; `thermal-conductivity` [W/m/K] |
| [cement-mortar-cfast-sp1041](../data/matdb/seed-v1/cement-mortar-cfast-sp1041/manifest.tsv) | 5 | 5 | 0 | 0 | 5.58 | `density` [kg/m3]; `emissivity` [1]; `nominal-thickness` [m]; `specific-heat-capacity` [J/kg/K]; `thermal-conductivity` [W/m/K] |
| [cherry-black-fpl-gtr282](../data/matdb/seed-v1/cherry-black-fpl-gtr282/manifest.tsv) | 14 | 14 | 0 | 0 | 3.82 | `density` [kg/m3]; `er_over_el` [1]; `et_over_el` [1]; `glr_over_el` [1]; `glt_over_el` [1]; `modulus_of_elasticity_bending` [MPa]; `nu_lr` [1]; `nu_lt` [1]; `nu_rl` [1]; `nu_rt` [1]; `nu_tl` [1]; `nu_tr` [1]; `specific_gravity` [1]; `young_modulus_longitudinal` [MPa] |
| [concrete-lightweight-cfast-sp1041](../data/matdb/seed-v1/concrete-lightweight-cfast-sp1041/manifest.tsv) | 5 | 5 | 0 | 0 | 5.65 | `density` [kg/m3]; `emissivity` [1]; `nominal-thickness` [m]; `specific-heat-capacity` [J/kg/K]; `thermal-conductivity` [W/m/K] |
| [concrete-normalweight-cfast-sp1041](../data/matdb/seed-v1/concrete-normalweight-cfast-sp1041/manifest.tsv) | 5 | 5 | 0 | 0 | 5.69 | `density` [kg/m3]; `emissivity` [1]; `nominal-thickness` [m]; `specific-heat-capacity` [J/kg/K]; `thermal-conductivity` [W/m/K] |
| [concrete-nsc-mix-iv-nistir6475](../data/matdb/seed-v1/concrete-nsc-mix-iv-nistir6475/manifest.tsv) | 4 | 4 | 0 | 0 | 5.89 | `compressive-strength` [MPa]; `density` [kg/m3]; `dynamic-youngs-modulus` [Pa]; `static-compressive-youngs-modulus` [MPa] |
| [copper-c11000-mil-hdbk-698a](../data/matdb/seed-v1/copper-c11000-mil-hdbk-698a/manifest.tsv) | 5 | 5 | 0 | 0 | 2.13 | `density_printed_range_high` [kg/m3]; `density_printed_range_low` [kg/m3]; `melting_point_liquidus` [K]; `melting_point_solidus` [K]; `thermal_conductivity` [W/m/K] |
| [copper-fusion-nasa-tp3287](../data/matdb/seed-v1/copper-fusion-nasa-tp3287/manifest.tsv) | 2 | 2 | 0 | 0 | 1.80 | `latent-heat-fusion` [J/kg]; `melting-point` [K] |
| [copper-liquid-nasa-tp3287](../data/matdb/seed-v1/copper-liquid-nasa-tp3287/manifest.tsv) | 2 | 0 | 2 | 6 | 2.62 | `specific-enthalpy-reference-29815k` [J/kg]; `specific-heat-capacity` [J/kg/K] |
| [copper-pure-nasa-cr-71699](../data/matdb/seed-v1/copper-pure-nasa-cr-71699/manifest.tsv) | 1 | 0 | 1 | 15 | 2.58 | `thermal-conductivity` [W/m/K] |
| [copper-pure-nbs-c447](../data/matdb/seed-v1/copper-pure-nbs-c447/manifest.tsv) | 5 | 5 | 0 | 0 | 2.63 | `density` [kg/m3]; `hardness_brinell_annealed` [1]; `hardness_brinell_cold_drawn_56pct` [1]; `hardness_rockwell_b_cold_drawn_29pct` [1]; `melting_point` [K] |
| [copper-solid-nasa-tp3287](../data/matdb/seed-v1/copper-solid-nasa-tp3287/manifest.tsv) | 2 | 0 | 2 | 28 | 3.12 | `specific-enthalpy-reference-29815k` [J/kg]; `specific-heat-capacity` [J/kg/K] |
| [cottonwood-eastern-fpl-gtr282](../data/matdb/seed-v1/cottonwood-eastern-fpl-gtr282/manifest.tsv) | 14 | 14 | 0 | 0 | 3.83 | `density` [kg/m3]; `er_over_el` [1]; `et_over_el` [1]; `glr_over_el` [1]; `glt_over_el` [1]; `modulus_of_elasticity_bending` [MPa]; `nu_lr` [1]; `nu_lt` [1]; `nu_rl` [1]; `nu_rt` [1]; `nu_tl` [1]; `nu_tr` [1]; `specific_gravity` [1]; `young_modulus_longitudinal` [MPa] |
| [douglas-fir-coast-fpl-gtr282](../data/matdb/seed-v1/douglas-fir-coast-fpl-gtr282/manifest.tsv) | 15 | 15 | 0 | 0 | 3.98 | `density` [kg/m3]; `er_over_el` [1]; `et_over_el` [1]; `glr_over_el` [1]; `glt_over_el` [1]; `grt_over_el` [1]; `modulus_of_elasticity_bending` [MPa]; `nu_lr` [1]; `nu_lt` [1]; `nu_rl` [1]; `nu_rt` [1]; `nu_tl` [1]; `nu_tr` [1]; `specific_gravity` [1]; `young_modulus_longitudinal` [MPa] |
| [face-g-cdtrf-g-2023-v1](../data/matdb/seed-v1/face-g-cdtrf-g-2023-v1/manifest.tsv) | 6 | 7 | 0 | 0 | 6.69 | `cyclohexane_component_volume_fraction` [%]; `diisobutylene_component_volume_fraction` [%]; `isooctane_component_volume_fraction` [%]; `n_heptane_component_volume_fraction` [%]; `reported_calculated_research_octane_number` [1]; `toluene_component_volume_fraction` [%] |
| [glass-borosilicate-duran-ntrs-19860021558](../data/matdb/seed-v1/glass-borosilicate-duran-ntrs-19860021558/manifest.tsv) | 8 | 8 | 0 | 0 | 5.30 | `density` [kg/m3]; `glass-transition-temperature` [K]; `mean-linear-thermal-expansion-coefficient` [K^-1]; `poisson-ratio` [1]; `specific-heat-capacity` [J/kg/K]; `thermal-conductivity` [W/m/K]; `thermal-shock-resistance-temperature-difference` [K]; `young-modulus` [GPa] |
| [glass-fiber-insulation-cfast](../data/matdb/seed-v1/glass-fiber-insulation-cfast/manifest.tsv) | 5 | 5 | 0 | 0 | 5.37 | `density` [kg/m3]; `emissivity` [1]; `nominal-thickness` [m]; `specific-heat-capacity` [J/kg/K]; `thermal-conductivity` [W/m/K] |
| [glass-soda-lime-srm-1826b-nist](../data/matdb/seed-v1/glass-soda-lime-srm-1826b-nist/manifest.tsv) | 2 | 2 | 0 | 0 | 2.33 | `density` [kg/m3]; `linear-thermal-expansion-coefficient` [K^-1] |
| [gray-cast-iron-s2-s](../data/matdb/seed-v1/gray-cast-iron-s2-s/manifest.tsv) | 15 | 15 | 0 | 0 | 4.99 | `carbon_equivalent_ce` [%]; `carbon_mass_fraction` [%]; `copper_mass_fraction` [%]; `eutectic_colony_areal_density` [cm-2]; `graphite_area_fraction` [%]; `manganese_mass_fraction` [%]; `maximum_graphite_flake_length` [um]; `molybdenum_mass_fraction` [%]; `phosphorus_mass_fraction` [%]; `primary_dendrite_area_fraction` [%]; `silicon_mass_fraction` [%]; `sulfur_mass_fraction` [%]; `thermal_conductivity` [W/m/K]; `tin_mass_fraction` [%]; `ultimate_tensile_strength` [MPa] |
| [gypsum-board-5-8-cfast](../data/matdb/seed-v1/gypsum-board-5-8-cfast/manifest.tsv) | 5 | 5 | 0 | 0 | 5.47 | `density` [kg/m3]; `emissivity` [1]; `nominal-thickness` [m]; `specific-heat-capacity` [J/kg/K]; `thermal-conductivity` [W/m/K] |
| [gypsum-board-type-x-5-8-cfast](../data/matdb/seed-v1/gypsum-board-type-x-5-8-cfast/manifest.tsv) | 5 | 5 | 0 | 0 | 5.37 | `density` [kg/m3]; `emissivity` [1]; `nominal-thickness` [m]; `specific-heat-capacity` [J/kg/K]; `thermal-conductivity` [W/m/K] |
| [hemlock-western-fpl-gtr282](../data/matdb/seed-v1/hemlock-western-fpl-gtr282/manifest.tsv) | 13 | 13 | 0 | 0 | 3.73 | `density` [kg/m3]; `er_over_el` [1]; `et_over_el` [1]; `glr_over_el` [1]; `glt_over_el` [1]; `grt_over_el` [1]; `modulus_of_elasticity_bending` [MPa]; `nu_lr` [1]; `nu_lt` [1]; `nu_rt` [1]; `nu_tr` [1]; `specific_gravity` [1]; `young_modulus_longitudinal` [MPa] |
| [inconel-x750-nasa-cr-71699](../data/matdb/seed-v1/inconel-x750-nasa-cr-71699/manifest.tsv) | 1 | 0 | 1 | 15 | 2.10 | `thermal-conductivity` [W/m/K] |
| [iron-pure-nasa-cr-71699](../data/matdb/seed-v1/iron-pure-nasa-cr-71699/manifest.tsv) | 1 | 0 | 1 | 12 | 2.16 | `thermal-conductivity` [W/m/K] |
| [jinshan-n42-pristine-temperature](../data/matdb/seed-v1/jinshan-n42-pristine-temperature/manifest.tsv) | 5 | 8 | 0 | 0 | 7.80 | `coercivity_temperature_coefficient` [K-1]; `intrinsic_coercive_field_strength` [kA/m]; `maximum_magnetic_energy_product` [kJ/m3]; `remanence_temperature_coefficient` [K-1]; `remanent_flux_density` [T] |
| [kim-baek-2026-y30-afcp-demagnetization](../data/matdb/seed-v1/kim-baek-2026-y30-afcp-demagnetization/manifest.tsv) | 1 | 2 | 0 | 0 | 5.13 | `application_model_maximum_demagnetization_fraction` [%] |
| [larch-western-fpl-gtr282](../data/matdb/seed-v1/larch-western-fpl-gtr282/manifest.tsv) | 13 | 13 | 0 | 0 | 3.73 | `density` [kg/m3]; `er_over_el` [1]; `et_over_el` [1]; `glr_over_el` [1]; `glt_over_el` [1]; `grt_over_el` [1]; `modulus_of_elasticity_bending` [MPa]; `nu_lr` [1]; `nu_lt` [1]; `nu_rt` [1]; `nu_tr` [1]; `specific_gravity` [1]; `young_modulus_longitudinal` [MPa] |
| [lead-cast-expansion-nbs-rp500](../data/matdb/seed-v1/lead-cast-expansion-nbs-rp500/manifest.tsv) | 2 | 5 | 0 | 0 | 6.04 | `density` [kg/m3]; `mean-linear-expansion-coefficient-from-20c` [K^-1] |
| [lead-fusion-nasa-tp3287](../data/matdb/seed-v1/lead-fusion-nasa-tp3287/manifest.tsv) | 2 | 2 | 0 | 0 | 1.88 | `latent-heat-fusion` [J/kg]; `melting-point` [K] |
| [lead-liquid-nasa-tp3287](../data/matdb/seed-v1/lead-liquid-nasa-tp3287/manifest.tsv) | 2 | 0 | 2 | 16 | 2.90 | `specific-enthalpy-reference-29815k` [J/kg]; `specific-heat-capacity` [J/kg/K] |
| [lead-liquid-transport-nasa-cr144016](../data/matdb/seed-v1/lead-liquid-transport-nasa-cr144016/manifest.tsv) | 4 | 0 | 4 | 18 | 5.89 | `density` [kg/m3]; `dynamic-viscosity` [Pa*s]; `surface-tension` [N/m]; `thermal-conductivity` [W/m/K] |
| [lead-pure-nbs-c447](../data/matdb/seed-v1/lead-pure-nbs-c447/manifest.tsv) | 4 | 4 | 0 | 0 | 2.23 | `density` [kg/m3]; `hardness_brinell_die_cast` [1]; `hardness_vickers_extruded` [1]; `melting_point` [K] |
| [lead-solid-conductivity-nbs-rp668](../data/matdb/seed-v1/lead-solid-conductivity-nbs-rp668/manifest.tsv) | 1 | 0 | 1 | 4 | 2.93 | `thermal-conductivity` [W/m/K] |
| [lead-solid-nasa-tp3287](../data/matdb/seed-v1/lead-solid-nasa-tp3287/manifest.tsv) | 2 | 0 | 2 | 20 | 2.91 | `specific-enthalpy-reference-29815k` [J/kg]; `specific-heat-capacity` [J/kg/K] |
| [lip-perioral-invivo-mdpi-ma17153654](../data/matdb/seed-v1/lip-perioral-invivo-mdpi-ma17153654/manifest.tsv) | 3 | 3 | 0 | 0 | 1.76 | `young_modulus_effective_cheek_left` [kPa]; `young_modulus_effective_cheek_right` [kPa]; `young_modulus_effective_upper_lip` [kPa] |
| [mahogany-african-fpl-gtr282](../data/matdb/seed-v1/mahogany-african-fpl-gtr282/manifest.tsv) | 14 | 14 | 0 | 0 | 3.80 | `er_over_el` [1]; `et_over_el` [1]; `glr_over_el` [1]; `glt_over_el` [1]; `grt_over_el` [1]; `modulus_of_elasticity_bending` [MPa]; `nu_lr` [1]; `nu_lt` [1]; `nu_rl` [1]; `nu_rt` [1]; `nu_tl` [1]; `nu_tr` [1]; `specific_gravity` [1]; `young_modulus_longitudinal` [MPa] |
| [mahogany-honduras-fpl-gtr282](../data/matdb/seed-v1/mahogany-honduras-fpl-gtr282/manifest.tsv) | 14 | 14 | 0 | 0 | 3.81 | `er_over_el` [1]; `et_over_el` [1]; `glr_over_el` [1]; `glt_over_el` [1]; `grt_over_el` [1]; `modulus_of_elasticity_bending` [MPa]; `nu_lr` [1]; `nu_lt` [1]; `nu_rl` [1]; `nu_rt` [1]; `nu_tl` [1]; `nu_tr` [1]; `specific_gravity` [1]; `young_modulus_longitudinal` [MPa] |
| [maple-red-fpl-gtr282](../data/matdb/seed-v1/maple-red-fpl-gtr282/manifest.tsv) | 14 | 14 | 0 | 0 | 3.81 | `density` [kg/m3]; `er_over_el` [1]; `et_over_el` [1]; `glr_over_el` [1]; `glt_over_el` [1]; `modulus_of_elasticity_bending` [MPa]; `nu_lr` [1]; `nu_lt` [1]; `nu_rl` [1]; `nu_rt` [1]; `nu_tl` [1]; `nu_tr` [1]; `specific_gravity` [1]; `young_modulus_longitudinal` [MPa] |
| [maple-sugar-fpl-gtr282](../data/matdb/seed-v1/maple-sugar-fpl-gtr282/manifest.tsv) | 14 | 14 | 0 | 0 | 3.82 | `density` [kg/m3]; `er_over_el` [1]; `et_over_el` [1]; `glr_over_el` [1]; `glt_over_el` [1]; `modulus_of_elasticity_bending` [MPa]; `nu_lr` [1]; `nu_lt` [1]; `nu_rl` [1]; `nu_rt` [1]; `nu_tl` [1]; `nu_tr` [1]; `specific_gravity` [1]; `young_modulus_longitudinal` [MPa] |
| [maple-sycamore-danihelova2022](../data/matdb/seed-v1/maple-sycamore-danihelova2022/manifest.tsv) | 4 | 4 | 0 | 0 | 2.38 | `density` [kg/m3]; `dynamic_young_modulus_longitudinal` [MPa]; `log_decrement_longitudinal` [1]; `loss_factor_longitudinal` [1] |
| [music-wire-nbs-c447](../data/matdb/seed-v1/music-wire-nbs-c447/manifest.tsv) | 3 | 3 | 0 | 0 | 2.24 | `density` [kg/m3]; `tensile_strength_diameter_0p01in` [MPa]; `tensile_strength_diameter_0p20in` [MPa] |
| [music-wire-spring-moduli-fuchs-1968](../data/matdb/seed-v1/music-wire-spring-moduli-fuchs-1968/manifest.tsv) | 2 | 2 | 0 | 0 | 1.72 | `shear_modulus` [MPa]; `young_modulus` [MPa] |
| [n0602-001-nitrile-jp8-compatibility](../data/matdb/seed-v1/n0602-001-nitrile-jp8-compatibility/manifest.tsv) | 8 | 10 | 0 | 0 | 7.26 | `absorbed_fuel_volume_fraction` [%]; `jp8_alkane_fuel_polymer_partition_coefficient` [1]; `jp8_aromatic_fuel_polymer_partition_coefficient` [1]; `jp8_aromatic_to_alkane_partition_ratio` [1]; `jp8_volume_swell_aromatic_fraction_r_squared` [1]; `jp8_volume_swell_per_aromatic_volume_fraction` [1]; `jp8_volume_swell_zero_aromatic_intercept` [%]; `tga_semivolatile_mass_fraction` [%] |
| [naca-tn-2680-isooctane-flame-speed](../data/matdb/seed-v1/naca-tn-2680-isooctane-flame-speed/manifest.tsv) | 2 | 16 | 0 | 0 | 19.29 | `maximum_laminar_flame_speed` [cm/s]; `minimum_reported_fuel_mole_fraction_purity` [%] |
| [napc-pe-5-l-1274-gear-oil](../data/matdb/seed-v1/napc-pe-5-l-1274-gear-oil/manifest.tsv) | 3 | 3 | 0 | 0 | 1.71 | `flash_point_temperature` [K]; `reported_specific_gravity` [1]; `total_acid_number_as_koh_mass_per_oil_mass` [mg/g] |
| [napc-pe-5-l-1307-1553-gear-oil](../data/matdb/seed-v1/napc-pe-5-l-1307-1553-gear-oil/manifest.tsv) | 4 | 6 | 0 | 0 | 3.21 | `flash_point_temperature` [K]; `pour_point_temperature` [K]; `reported_specific_gravity` [1]; `total_acid_number_as_koh_mass_per_oil_mass` [mg/g] |
| [nasa-cr-115153-water-ethylene-glycol](../data/matdb/seed-v1/nasa-cr-115153-water-ethylene-glycol/manifest.tsv) | 9 | 11 | 0 | 0 | 6.19 | `density` [kg/m3]; `sodium_benzoate_mass_fraction_lower_bound` [%]; `sodium_benzoate_mass_fraction_upper_bound` [%]; `sodium_nitrite_mass_fraction_lower_bound` [%]; `sodium_nitrite_mass_fraction_upper_bound` [%]; `specific_heat_capacity` [J/kg/K]; `thermal_conductivity` [W/m/K]; `water_mass_fraction_lower_bound` [%]; `water_mass_fraction_upper_bound` [%] |
| [nasa-cr-195445-omc-ps200-rotary-coating](../data/matdb/seed-v1/nasa-cr-195445-omc-ps200-rotary-coating/manifest.tsv) | 4 | 7 | 0 | 0 | 9.14 | `ps200_baf2_caf2_eutectic_feedstock_mass_fraction` [%]; `ps200_bonded_chromium_carbide_feedstock_mass_fraction` [%]; `ps200_silver_feedstock_mass_fraction` [%]; `surface_roughness_rms` [um] |
| [nasa-cr-4538-tempel-24n208-m19](../data/matdb/seed-v1/nasa-cr-4538-tempel-24n208-m19/manifest.tsv) | 3 | 3 | 0 | 0 | 5.48 | `lamination_thickness` [m]; `nominal_silicon_mass_fraction` [%]; `specific_hysteresis_loss_rating` [W/kg] |
| [nasa-tn-d-8184-m19-material-deck](../data/matdb/seed-v1/nasa-tn-d-8184-m19-material-deck/manifest.tsv) | 6 | 5 | 1 | 14 | 8.94 | `core_loss_frequency_power_law_exponent` [1]; `core_loss_reference_flux_density` [T]; `core_loss_reference_frequency` [Hz]; `lamination_thickness` [m]; `magnetic_flux_density` [T]; `specific_core_loss` [W/kg] |
| [nasa-uam-cooltherm-ep2000-180c-cure](../data/matdb/seed-v1/nasa-uam-cooltherm-ep2000-180c-cure/manifest.tsv) | 2 | 2 | 0 | 0 | 2.89 | `highest_completed_post_cure_temperature` [degC]; `manufacturer_recommended_final_cure_temperature` [degC] |
| [nasa-uam-mw16c-polyimide-magnet-wire](../data/matdb/seed-v1/nasa-uam-mw16c-polyimide-magnet-wire/manifest.tsv) | 2 | 2 | 0 | 0 | 3.12 | `thermal_endurance_reference_duration` [h]; `thermal_endurance_reference_temperature` [degC] |
| [nasa-uam-nomex-410-slot-liner](../data/matdb/seed-v1/nasa-uam-nomex-410-slot-liner/manifest.tsv) | 1 | 1 | 0 | 0 | 1.92 | `selected_slot_liner_thickness` [mm] |
| [ngyc-n42-sintered-nickel-coated](../data/matdb/seed-v1/ngyc-n42-sintered-nickel-coated/manifest.tsv) | 3 | 4 | 0 | 0 | 3.24 | `coercive_field_strength` [kA/m]; `maximum_magnetic_energy_product` [kJ/m3]; `remanent_flux_density` [mT] |
| [nickel-pure-nasa-cr-71699](../data/matdb/seed-v1/nickel-pure-nasa-cr-71699/manifest.tsv) | 1 | 0 | 1 | 13 | 2.03 | `thermal-conductivity` [W/m/K] |
| [nist-srm-1720-northern-continental-air](../data/matdb/seed-v1/nist-srm-1720-northern-continental-air/manifest.tsv) | 4 | 4 | 0 | 0 | 6.27 | `information_argon_amount_fraction` [%]; `information_carbon_monoxide_amount_fraction_lower_bound` [%]; `information_carbon_monoxide_amount_fraction_upper_bound` [%]; `information_oxygen_amount_fraction` [%] |
| [nist-srm-2728-auto-emission-reference-gas](../data/matdb/seed-v1/nist-srm-2728-auto-emission-reference-gas/manifest.tsv) | 4 | 4 | 0 | 0 | 6.28 | `information_total_other_hydrocarbons_propane_equivalent_amount_fraction` [%]; `nominal_carbon_dioxide_amount_fraction` [%]; `nominal_carbon_monoxide_amount_fraction` [%]; `nominal_propane_amount_fraction` [%] |
| [ofhc-copper-rrr100](../data/matdb/seed-v1/ofhc-copper-rrr100/manifest.tsv) | 2 | 0 | 2 | 48 | 4.29 | `specific-heat-capacity` [J/kg/K]; `thermal-conductivity` [W/m/K] |
| [osb-aspen-pu-1992-mill7-fpl-gtr282](../data/matdb/seed-v1/osb-aspen-pu-1992-mill7-fpl-gtr282/manifest.tsv) | 6 | 6 | 0 | 0 | 4.22 | `bending-moe-parallel` [GPa]; `bending-moe-perpendicular` [GPa]; `bending-mor-parallel` [MPa]; `bending-mor-perpendicular` [MPa]; `internal-bond-strength` [MPa]; `specific-gravity` [1] |
| [osb-southern-pine-biblis-1989-mill1-fpl-gtr282](../data/matdb/seed-v1/osb-southern-pine-biblis-1989-mill1-fpl-gtr282/manifest.tsv) | 6 | 6 | 0 | 0 | 4.32 | `bending-moe-parallel` [GPa]; `bending-moe-perpendicular` [GPa]; `bending-mor-parallel` [MPa]; `bending-mor-perpendicular` [MPa]; `internal-bond-strength` [MPa]; `specific-gravity` [1] |
| [peek-nasa-thermic-plate](../data/matdb/seed-v1/peek-nasa-thermic-plate/manifest.tsv) | 3 | 9 | 0 | 0 | 4.83 | `density` [kg/m3]; `specific_heat_capacity` [J/kg/K]; `thermal_conductivity` [W/m/K] |
| [pennzane-shf-x-2000-bearing-oil](../data/matdb/seed-v1/pennzane-shf-x-2000-bearing-oil/manifest.tsv) | 5 | 7 | 0 | 0 | 2.55 | `density` [g/mL]; `flash_point_temperature` [degC]; `kinematic_viscosity` [mm^2/s]; `pour_point_temperature` [degC]; `viscosity_index_scale_reading` [1] |
| [phosphor-bronze-5a-mil-hdbk-698a](../data/matdb/seed-v1/phosphor-bronze-5a-mil-hdbk-698a/manifest.tsv) | 3 | 3 | 0 | 0 | 1.98 | `density` [kg/m3]; `shear_modulus` [MPa]; `young_modulus` [MPa] |
| [phosphor-bronze-c51000-mil-hdbk-698a](../data/matdb/seed-v1/phosphor-bronze-c51000-mil-hdbk-698a/manifest.tsv) | 4 | 4 | 0 | 0 | 1.98 | `density_printed_range_low` [kg/m3]; `melting_point_liquidus` [K]; `melting_point_solidus` [K]; `thermal_conductivity` [W/m/K] |
| [phosphor-bronze-c51000-nist-mono177](../data/matdb/seed-v1/phosphor-bronze-c51000-nist-mono177/manifest.tsv) | 2 | 2 | 0 | 0 | 1.92 | `poisson_ratio` [1]; `young_modulus` [MPa] |
| [plywood-douglas-fir-fpl-gtr282](../data/matdb/seed-v1/plywood-douglas-fir-fpl-gtr282/manifest.tsv) | 6 | 6 | 0 | 0 | 4.52 | `bending-moe` [GPa]; `bending-mor` [MPa]; `fiber-stress-proportional-limit` [MPa]; `glue-line-shear-strength` [MPa]; `rail-shear-strength` [MPa]; `specific-gravity` [1] |
| [plywood-southern-pine-fpl-gtr282](../data/matdb/seed-v1/plywood-southern-pine-fpl-gtr282/manifest.tsv) | 6 | 6 | 0 | 0 | 4.54 | `bending-moe` [GPa]; `bending-mor` [MPa]; `fiber-stress-proportional-limit` [MPa]; `glue-line-shear-strength` [MPa]; `rail-shear-strength` [MPa]; `specific-gravity` [1] |
| [ptfe-teflon-nist-cryogenic](../data/matdb/seed-v1/ptfe-teflon-nist-cryogenic/manifest.tsv) | 2 | 4 | 0 | 0 | 2.53 | `specific_heat_capacity` [J/kg/K]; `thermal_conductivity` [W/m/K] |
| [pvc-cryogenic-nist](../data/matdb/seed-v1/pvc-cryogenic-nist/manifest.tsv) | 2 | 3 | 0 | 0 | 3.85 | `specific-heat-capacity` [J/kg/K]; `thermal-conductivity` [W/m/K] |
| [redwood-old-growth-fpl-gtr282](../data/matdb/seed-v1/redwood-old-growth-fpl-gtr282/manifest.tsv) | 13 | 13 | 0 | 0 | 3.74 | `density` [kg/m3]; `er_over_el` [1]; `et_over_el` [1]; `glr_over_el` [1]; `glt_over_el` [1]; `grt_over_el` [1]; `modulus_of_elasticity_bending` [MPa]; `nu_lr` [1]; `nu_lt` [1]; `nu_rt` [1]; `nu_tr` [1]; `specific_gravity` [1]; `young_modulus_longitudinal` [MPa] |
| [rheolube-2000-pennzane-grease](../data/matdb/seed-v1/rheolube-2000-pennzane-grease/manifest.tsv) | 3 | 3 | 0 | 0 | 1.89 | `density` [g/cm3]; `nlgi_consistency_grade` [1]; `oil_separation_mass_fraction` [%] |
| [rosewood-brazilian-fpl-gtr282](../data/matdb/seed-v1/rosewood-brazilian-fpl-gtr282/manifest.tsv) | 3 | 3 | 0 | 0 | 2.35 | `modulus_of_elasticity_bending` [MPa]; `specific_gravity` [1]; `young_modulus_longitudinal` [MPa] |
| [rosewood-indian-fpl-gtr282](../data/matdb/seed-v1/rosewood-indian-fpl-gtr282/manifest.tsv) | 3 | 3 | 0 | 0 | 2.35 | `modulus_of_elasticity_bending` [MPa]; `specific_gravity` [1]; `young_modulus_longitudinal` [MPa] |
| [sapele-fpl-gtr190](../data/matdb/seed-v1/sapele-fpl-gtr190/manifest.tsv) | 7 | 7 | 0 | 0 | 2.86 | `bending_modulus_of_elasticity` [GPa]; `compression_parallel_max` [MPa]; `modulus_of_rupture` [MPa]; `shear_parallel_max` [MPa]; `side_hardness` [N]; `specific_gravity_green_basis` [1]; `work_to_maximum_load` [kJ/m3] |
| [sjolund-2020-y30-catalog-model-inputs](../data/matdb/seed-v1/sjolund-2020-y30-catalog-model-inputs/manifest.tsv) | 5 | 5 | 0 | 0 | 7.63 | `coercive_field_strength` [kA/m]; `intrinsic_coercive_field_strength` [kA/m]; `maximum_magnetic_energy_product` [kJ/m3]; `model_relative_permeability` [1]; `remanent_flux_density` [mT] |
| [spruce-engelmann-fpl-gtr282](../data/matdb/seed-v1/spruce-engelmann-fpl-gtr282/manifest.tsv) | 15 | 15 | 0 | 0 | 3.97 | `density` [kg/m3]; `er_over_el` [1]; `et_over_el` [1]; `glr_over_el` [1]; `glt_over_el` [1]; `grt_over_el` [1]; `modulus_of_elasticity_bending` [MPa]; `nu_lr` [1]; `nu_lt` [1]; `nu_rl` [1]; `nu_rt` [1]; `nu_tl` [1]; `nu_tr` [1]; `specific_gravity` [1]; `young_modulus_longitudinal` [MPa] |
| [spruce-norway-danihelova2022](../data/matdb/seed-v1/spruce-norway-danihelova2022/manifest.tsv) | 4 | 4 | 0 | 0 | 2.34 | `density` [kg/m3]; `dynamic_young_modulus_longitudinal` [MPa]; `log_decrement_longitudinal` [1]; `loss_factor_longitudinal` [1] |
| [spruce-sitka-fpl-gtr282](../data/matdb/seed-v1/spruce-sitka-fpl-gtr282/manifest.tsv) | 15 | 15 | 0 | 0 | 3.97 | `density` [kg/m3]; `er_over_el` [1]; `et_over_el` [1]; `glr_over_el` [1]; `glt_over_el` [1]; `grt_over_el` [1]; `modulus_of_elasticity_bending` [MPa]; `nu_lr` [1]; `nu_lt` [1]; `nu_rl` [1]; `nu_rt` [1]; `nu_tl` [1]; `nu_tr` [1]; `specific_gravity` [1]; `young_modulus_longitudinal` [MPa] |
| [spruce-sitka-loss-qiu2026](../data/matdb/seed-v1/spruce-sitka-loss-qiu2026/manifest.tsv) | 2 | 2 | 0 | 0 | 2.45 | `loss_factor_longitudinal` [1]; `loss_factor_longitudinal_mode5` [1] |
| [stainless-17-4ph-h1025-bar-mil-hdbk-5j](../data/matdb/seed-v1/stainless-17-4ph-h1025-bar-mil-hdbk-5j/manifest.tsv) | 9 | 9 | 0 | 0 | 3.16 | `compression_modulus` [GPa]; `compressive_yield_l_s_basis` [MPa]; `elongation_l_s_basis` [%]; `poisson_ratio` [1]; `shear_modulus` [GPa]; `shear_ultimate_s_basis` [MPa]; `tensile_ultimate_l_s_basis` [MPa]; `tensile_yield_l_s_basis` [MPa]; `young_modulus` [GPa] |
| [stainless-301-annealed-mil-hdbk-5j](../data/matdb/seed-v1/stainless-301-annealed-mil-hdbk-5j/manifest.tsv) | 18 | 17 | 1 | 4 | 5.99 | `compression_modulus_l` [GPa]; `compressive_yield_l_s_basis` [MPa]; `compressive_yield_lt_s_basis` [MPa]; `density` [kg/m3]; `elongation_lt_s_basis` [%]; `poisson_ratio` [1]; `shear_modulus` [GPa]; `shear_ultimate_s_basis` [MPa]; `tensile_ultimate_l_s_basis` [MPa]; `tensile_ultimate_lt_s_basis` [MPa]; `tensile_yield_fraction` [1]; `tensile_yield_fraction_1033k` [1]; `tensile_yield_fraction_478k` [1]; `tensile_yield_fraction_589k` [1]; `tensile_yield_fraction_811k` [1]; `tensile_yield_l_s_basis` [MPa]; `tensile_yield_lt_s_basis` [MPa]; `young_modulus_l` [GPa] |
| [stainless-301-full-hard-mil-hdbk-5j](../data/matdb/seed-v1/stainless-301-full-hard-mil-hdbk-5j/manifest.tsv) | 13 | 13 | 0 | 0 | 4.24 | `compressive_yield_l_b_basis` [MPa]; `compressive_yield_lt_b_basis` [MPa]; `density` [kg/m3]; `elongation_lt_minimum` [%]; `poisson_ratio` [1]; `shear_modulus` [GPa]; `shear_ultimate_b_basis` [MPa]; `tensile_ultimate_l_b_basis` [MPa]; `tensile_ultimate_lt_b_basis` [MPa]; `tensile_yield_l_b_basis` [MPa]; `tensile_yield_lt_b_basis` [MPa]; `young_modulus_l` [GPa]; `young_modulus_lt` [GPa] |
| [stainless-304-nist-cryogenic](../data/matdb/seed-v1/stainless-304-nist-cryogenic/manifest.tsv) | 5 | 0 | 5 | 54 | 8.35 | `linear-expansion-relative-to-293k` [1]; `linear-thermal-expansion-coefficient` [K^-1, model-derived]; `specific-heat-capacity` [J/kg/K]; `thermal-conductivity` [W/m/K]; `young-modulus` [GPa] |
| [stainless-304a-nasa-cr-71699](../data/matdb/seed-v1/stainless-304a-nasa-cr-71699/manifest.tsv) | 1 | 0 | 1 | 7 | 2.06 | `thermal-conductivity` [W/m/K] |
| [stainless-316-nist-cryogenic](../data/matdb/seed-v1/stainless-316-nist-cryogenic/manifest.tsv) | 5 | 0 | 5 | 55 | 8.53 | `linear-expansion-relative-to-293k` [1]; `linear-thermal-expansion-coefficient` [K^-1, model-derived]; `specific-heat-capacity` [J/kg/K]; `thermal-conductivity` [W/m/K]; `young-modulus` [GPa] |
| [stainless-347-nasa-cr-71699](../data/matdb/seed-v1/stainless-347-nasa-cr-71699/manifest.tsv) | 1 | 0 | 1 | 13 | 2.07 | `thermal-conductivity` [W/m/K] |
| [steel-4130-sheet-normalized-mil-hdbk-5j](../data/matdb/seed-v1/steel-4130-sheet-normalized-mil-hdbk-5j/manifest.tsv) | 10 | 10 | 0 | 0 | 3.58 | `compression_modulus` [GPa]; `compressive_yield_s_basis` [MPa]; `density` [kg/m3]; `elongation_t_minimum` [%]; `poisson_ratio` [1]; `shear_modulus` [GPa]; `shear_ultimate_s_basis` [MPa]; `tensile_ultimate_s_basis` [MPa]; `tensile_yield_s_basis` [MPa]; `young_modulus` [GPa] |
| [sweetgum-fpl-gtr282](../data/matdb/seed-v1/sweetgum-fpl-gtr282/manifest.tsv) | 15 | 15 | 0 | 0 | 3.97 | `density` [kg/m3]; `er_over_el` [1]; `et_over_el` [1]; `glr_over_el` [1]; `glt_over_el` [1]; `grt_over_el` [1]; `modulus_of_elasticity_bending` [MPa]; `nu_lr` [1]; `nu_lt` [1]; `nu_rl` [1]; `nu_rt` [1]; `nu_tl` [1]; `nu_tr` [1]; `specific_gravity` [1]; `young_modulus_longitudinal` [MPa] |
| [teak-fpl-gtr190](../data/matdb/seed-v1/teak-fpl-gtr190/manifest.tsv) | 7 | 7 | 0 | 0 | 2.81 | `bending_modulus_of_elasticity` [GPa]; `compression_parallel_max` [MPa]; `modulus_of_rupture` [MPa]; `shear_parallel_max` [MPa]; `side_hardness` [N]; `specific_gravity_green_basis` [1]; `work_to_maximum_load` [kJ/m3] |
| [tin-fusion-nasa-tp3287](../data/matdb/seed-v1/tin-fusion-nasa-tp3287/manifest.tsv) | 2 | 2 | 0 | 0 | 1.80 | `latent-heat-fusion` [J/kg]; `melting-point` [K] |
| [tin-liquid-nasa-tp3287](../data/matdb/seed-v1/tin-liquid-nasa-tp3287/manifest.tsv) | 2 | 0 | 2 | 12 | 2.85 | `specific-enthalpy-reference-29815k` [J/kg]; `specific-heat-capacity` [J/kg/K] |
| [tin-pure-nbs-c447](../data/matdb/seed-v1/tin-pure-nbs-c447/manifest.tsv) | 2 | 2 | 0 | 0 | 1.42 | `density` [kg/m3]; `melting_point` [K] |
| [tin-solid-nasa-tp3287](../data/matdb/seed-v1/tin-solid-nasa-tp3287/manifest.tsv) | 2 | 0 | 2 | 16 | 2.81 | `specific-enthalpy-reference-29815k` [J/kg]; `specific-heat-capacity` [J/kg/K] |
| [tissue-density-stated-plos-pcbi1004907](../data/matdb/seed-v1/tissue-density-stated-plos-pcbi1004907/manifest.tsv) | 1 | 1 | 0 | 0 | 1.31 | `density_stated_soft_tissue` [kg/m3] |
| [titanium-a110at-nasa-cr-71699](../data/matdb/seed-v1/titanium-a110at-nasa-cr-71699/manifest.tsv) | 1 | 0 | 1 | 7 | 2.19 | `thermal-conductivity` [W/m/K] |
| [titanium-pure-nasa-cr-71699](../data/matdb/seed-v1/titanium-pure-nasa-cr-71699/manifest.tsv) | 1 | 0 | 1 | 9 | 2.05 | `thermal-conductivity` [W/m/K] |
| [torrent-2018-m19-steinmetz-inputs](../data/matdb/seed-v1/torrent-2018-m19-steinmetz-inputs/manifest.tsv) | 10 | 10 | 0 | 0 | 12.17 | `equation_4_flux_density_exponent_n` [1]; `equation_4_frequency_exponent_a` [1]; `equation_4_reported_k_h_numeric` [1]; `equation_4_reported_output_scale_numeric` [1]; `equation_5_flux_density_exponent_z` [1]; `equation_5_frequency_exponent_x` [1]; `equation_5_reported_k_f_numeric` [1]; `equation_5_reported_output_scale_numeric` [1]; `equation_5_sheet_thickness_e` [m]; `equation_5_thickness_exponent_y` [1] |
| [twomass-if72-standard-plos-pone0187486](../data/matdb/seed-v1/twomass-if72-standard-plos-pone0187486/manifest.tsv) | 6 | 6 | 0 | 0 | 2.27 | `model_fold_length` [cm]; `model_mass_lower` [g]; `model_mass_upper` [g]; `model_stiffness_coupling` [N/m]; `model_stiffness_lower` [N/m]; `model_stiffness_upper` [N/m] |
| [urethane-rigid-foam-insulation-cfast](../data/matdb/seed-v1/urethane-rigid-foam-insulation-cfast/manifest.tsv) | 5 | 5 | 0 | 0 | 5.44 | `density` [kg/m3]; `emissivity` [1]; `nominal-thickness` [m]; `specific-heat-capacity` [J/kg/K]; `thermal-conductivity` [W/m/K] |
| [vocalfold-felid-rheometry-plos-pone0027029](../data/matdb/seed-v1/vocalfold-felid-rheometry-plos-pone0027029/manifest.tsv) | 6 | 6 | 0 | 0 | 2.43 | `loss_tangent_lion` [1]; `loss_tangent_tiger` [1]; `shear_loss_modulus_printed_range_high` [Pa]; `shear_loss_modulus_printed_range_low` [Pa]; `shear_storage_modulus_printed_range_high` [Pa]; `shear_storage_modulus_printed_range_low` [Pa] |
| [vocalfold-human-indentation-jvoice-pmc12180296](../data/matdb/seed-v1/vocalfold-human-indentation-jvoice-pmc12180296/manifest.tsv) | 3 | 3 | 0 | 0 | 1.91 | `young_modulus_effective_inferior` [kPa]; `young_modulus_effective_medial` [kPa]; `young_modulus_effective_superior` [kPa] |
| [vocalfold-porcine-aspiration-mdpi-s21092923](../data/matdb/seed-v1/vocalfold-porcine-aspiration-mdpi-s21092923/manifest.tsv) | 4 | 4 | 0 | 0 | 2.18 | `young_modulus_dynamic_day_printed_range_high` [kPa]; `young_modulus_dynamic_day_printed_range_low` [kPa]; `young_modulus_dynamic_printed_range_high` [kPa]; `young_modulus_dynamic_printed_range_low` [kPa] |
| [walnut-black-fpl-gtr282](../data/matdb/seed-v1/walnut-black-fpl-gtr282/manifest.tsv) | 15 | 15 | 0 | 0 | 3.97 | `density` [kg/m3]; `er_over_el` [1]; `et_over_el` [1]; `glr_over_el` [1]; `glt_over_el` [1]; `grt_over_el` [1]; `modulus_of_elasticity_bending` [MPa]; `nu_lr` [1]; `nu_lt` [1]; `nu_rl` [1]; `nu_rt` [1]; `nu_tl` [1]; `nu_tr` [1]; `specific_gravity` [1]; `young_modulus_longitudinal` [MPa] |
| [water-liquid-iapws-sr6-08](../data/matdb/seed-v1/water-liquid-iapws-sr6-08/manifest.tsv) | 5 | 0 | 5 | 75 | 8.43 | `density` [kg/m3]; `dynamic-viscosity` [Pa*s]; `specific-enthalpy` [J/kg]; `specific-heat-capacity` [J/kg/K]; `thermal-conductivity` [W/m/K] |
| [waterborne-preservative-treated-lumber-fpl-gtr282](../data/matdb/seed-v1/waterborne-preservative-treated-lumber-fpl-gtr282/manifest.tsv) | 4 | 4 | 0 | 0 | 4.97 | `post-treatment-redrying-standard-temperature-limit` [K]; `post-treatment-redrying-temperature-threshold` [K]; `waterborne-retention-high-reference` [kg/m3]; `waterborne-retention-low-threshold` [kg/m3] |
| [wo2018-125520-formulation-8-5w30](../data/matdb/seed-v1/wo2018-125520-formulation-8-5w30/manifest.tsv) | 11 | 12 | 0 | 0 | 9.02 | `cold_cranking_simulator_dynamic_viscosity` [mPa*s]; `high_temperature_high_shear_dynamic_viscosity` [mPa*s]; `infineum_p6003_component_mass_fraction` [%]; `kinematic_viscosity` [mm^2/s]; `mini_rotary_viscometer_dynamic_viscosity` [mPa*s]; `noack_mass_loss_fraction` [%]; `pour_point_temperature` [degC]; `spectrasyn_4_component_mass_fraction` [%]; `spectrasyn_elite_150_component_mass_fraction` [%]; `synesstic_5_component_mass_fraction` [%]; `viscosity_index_scale_reading` [1] |
| [yellow-poplar-fpl-gtr282](../data/matdb/seed-v1/yellow-poplar-fpl-gtr282/manifest.tsv) | 15 | 15 | 0 | 0 | 3.98 | `density` [kg/m3]; `er_over_el` [1]; `et_over_el` [1]; `glr_over_el` [1]; `glt_over_el` [1]; `grt_over_el` [1]; `modulus_of_elasticity_bending` [MPa]; `nu_lr` [1]; `nu_lt` [1]; `nu_rl` [1]; `nu_rt` [1]; `nu_tl` [1]; `nu_tr` [1]; `specific_gravity` [1]; `young_modulus_longitudinal` [MPa] |
| [zinc-fusion-nasa-tp3287](../data/matdb/seed-v1/zinc-fusion-nasa-tp3287/manifest.tsv) | 2 | 2 | 0 | 0 | 1.80 | `latent-heat-fusion` [J/kg]; `melting-point` [K] |
| [zinc-liquid-nasa-tp3287](../data/matdb/seed-v1/zinc-liquid-nasa-tp3287/manifest.tsv) | 2 | 0 | 2 | 10 | 2.74 | `specific-enthalpy-reference-29815k` [J/kg]; `specific-heat-capacity` [J/kg/K] |
| [zinc-pure-nbs-c447](../data/matdb/seed-v1/zinc-pure-nbs-c447/manifest.tsv) | 2 | 2 | 0 | 0 | 1.42 | `density` [kg/m3]; `melting_point` [K] |
| [zinc-solid-nasa-tp3287](../data/matdb/seed-v1/zinc-solid-nasa-tp3287/manifest.tsv) | 2 | 0 | 2 | 14 | 2.70 | `specific-enthalpy-reference-29815k` [J/kg]; `specific-heat-capacity` [J/kg/K] |

## Interface bundles

These properties belong to an ordered pair of surfaces and a declared environment/history. They are not bulk properties of either constituent in isolation.

| Source bundle | Metrics | Scalar claims | Curve claims | Knots | KiB | Populated property names [source units] |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| [deaconescu-2020-ptfe-cf10-steel-bore-interface](../data/matdb/seed-v1/deaconescu-2020-ptfe-cf10-steel-bore-interface/manifest.tsv) | 1 | 7 | 0 | 0 | 13.38 | `fluid-film-thickness` [m] |
| [mahle-wo2019072721-ptfe-nickel-sic-coated-bore-interface](../data/matdb/seed-v1/mahle-wo2019072721-ptfe-nickel-sic-coated-bore-interface/manifest.tsv) | 11 | 11 | 0 | 0 | 15.35 | `patent_claim_honed_liner_arithmetic_mean_roughness_approximate` [m]; `patent_claim_liner_layer_thickness_lower_exclusive_bound` [m]; `patent_claim_liner_layer_thickness_preferred` [m]; `patent_claim_liner_layer_thickness_upper_exclusive_bound` [m]; `patent_claim_liner_silicon_carbide_mass_fraction_approximate_lower_bound` [%]; `patent_claim_liner_silicon_carbide_mass_fraction_approximate_upper_bound` [%]; `patent_claim_silicon_carbide_particle_diameter_lower_exclusive_bound` [m]; `patent_claim_silicon_carbide_particle_diameter_preferred` [m]; `patent_claim_silicon_carbide_particle_diameter_upper_exclusive_bound` [m]; `patent_embodiment_liner_silicon_carbide_mass_fraction_lower_bound` [%]; `patent_embodiment_liner_silicon_carbide_mass_fraction_upper_bound` [%] |
| [nasa-52100-dry-air-interface](../data/matdb/seed-v1/nasa-52100-dry-air-interface/manifest.tsv) | 1 | 1 | 0 | 0 | 2.30 | `kinetic-friction-coefficient` [1] |
| [nasa-52100-gxl320a-vacuum-interface](../data/matdb/seed-v1/nasa-52100-gxl320a-vacuum-interface/manifest.tsv) | 3 | 3 | 0 | 0 | 3.96 | `kinetic-friction-coefficient` [1]; `kinetic-friction-coefficient-observed-maximum` [1]; `kinetic-friction-coefficient-observed-minimum` [1] |
| [nasa-tn-d-2223-4340-high-lead-bronze-journal](../data/matdb/seed-v1/nasa-tn-d-2223-4340-high-lead-bronze-journal/manifest.tsv) | 1 | 1 | 0 | 0 | 3.05 | `maximum-demonstrated-unit-bearing-load` [Pa] |
| [yilmaz-2026-a2017-seiken-llc-ra005-wetting](../data/matdb/seed-v1/yilmaz-2026-a2017-seiken-llc-ra005-wetting/manifest.tsv) | 1 | 1 | 0 | 0 | 3.24 | `static-contact-angle` [1] |
| [yilmaz-2026-a2017-seiken-llc-ra3-wetting](../data/matdb/seed-v1/yilmaz-2026-a2017-seiken-llc-ra3-wetting/manifest.tsv) | 1 | 1 | 0 | 0 | 3.24 | `static-contact-angle` [1] |
| [zhang-2021-carbon-ptfe-cr-piston-rod](../data/matdb/seed-v1/zhang-2021-carbon-ptfe-cr-piston-rod/manifest.tsv) | 4 | 4 | 0 | 0 | 5.86 | `single-seal-instroke-friction-force-at-source-observed-maximum` [N]; `single-seal-instroke-friction-force-at-source-observed-minimum` [N]; `single-seal-outstroke-friction-force-at-source-observed-maximum` [N]; `single-seal-outstroke-friction-force-at-source-observed-minimum` [N] |

## Species records

Each row has two populated numerical metadata fields, molar mass and reference pressure. Identity, phase, EOS and reference convention are additional categorical metadata. These records do not supply heat capacity or a temperature-dependent thermodynamic model.

| Bundle | Species | Molar mass | Reference pressure | Phase / EOS | KiB |
| --- | --- | --- | --- | --- | ---: |
| [argon](../data/matdb/seed-v1/argon/manifest.tsv) | `Ar` | 39.94800 g/mol | 100 kPa | gas / ideal-gas | 0.60 |
| [carbon-dioxide](../data/matdb/seed-v1/carbon-dioxide/manifest.tsv) | `CO2` | 44.00950 g/mol | 100 kPa | gas / ideal-gas | 0.61 |
| [carbon-monoxide](../data/matdb/seed-v1/carbon-monoxide/manifest.tsv) | `CO` | 28.01010 g/mol | 100 kPa | gas / ideal-gas | 0.60 |
| [methane](../data/matdb/seed-v1/methane/manifest.tsv) | `CH4` | 16.04246 g/mol | 100 kPa | gas / ideal-gas | 0.61 |
| [nitrogen](../data/matdb/seed-v1/nitrogen/manifest.tsv) | `N2` | 28.01340 g/mol | 100 kPa | gas / ideal-gas | 0.60 |
| [oxygen](../data/matdb/seed-v1/oxygen/manifest.tsv) | `O2` | 31.99880 g/mol | 100 kPa | gas / ideal-gas | 0.60 |
| [water-vapor](../data/matdb/seed-v1/water-vapor/manifest.tsv) | `H2O` | 18.01528 g/mol | 100 kPa | gas / ideal-gas | 0.61 |

## Contact-law cards

These three cards have three authored law parameters each. They are Estimate-class inputs, not measured identification data; read their pair, geometry, validity and graze advisory before use.

| Bundle | Law | K [N/m^alpha] | alpha | chi [s/m] | KiB |
| --- | --- | ---: | ---: | ---: | ---: |
| [contact-cane-reed-on-mouthpiece-lay](../data/matdb/seed-v1/contact-cane-reed-on-mouthpiece-lay/manifest.tsv) | `penalty-power` | 1.0e7 | 2.0 | 1.0e-4 | 1.26 |
| [contact-jawari-bone-bridge](../data/matdb/seed-v1/contact-jawari-bone-bridge/manifest.tsv) | `penalty-power` | 1.0e8 | 2.0 | 0.0 | 1.27 |
| [contact-string-on-fret-nickel](../data/matdb/seed-v1/contact-string-on-fret-nickel/manifest.tsv) | `penalty-power` | 5.0e7 | 1.5 | 0.0 | 1.23 |

## What the counts do not establish

A count of populated properties is not a completeness score against a simulation's requirements. Many records are exact-temperature points or narrow source conditions. A source may supply density while lacking conductivity, heat capacity, mechanics, loss, optics or a constitutive model. A temperature curve for one property does not extend another property's domain. Material names do not make separately sourced specimens compatible.

The metric counts also include composition fractions, test ratings, manufacturing limits and other source-specific quantities stored as properties. Eighteen populated names do not necessarily mean eighteen constitutive inputs useful to a particular solver; the listed names show what is actually present.

For example, the two stainless cryogenic packs now contain four property names each, with continuous engineering approximations over declared intervals. They still lack density, an instantaneous expansion-coefficient law, process-specific specimen binding, and full coupled-model inputs. The two solid-lead packs added in this snapshot describe different specimens: RP500 sample 1144 has density and four interval-mean expansion measurements, while RP668 L.S. has a conductivity estimate with an explicitly assumed absolute calibration. They cannot silently be assembled into one condition-matched lead specimen.

The current seed contains no NASA-9 or kinetics source profiles, even though the offline compiler has code for those formats. The seven species records must not be mistaken for seven complete gas models. Likewise, the three contact cards follow their separate direct-loader path; the common bulk-pack compiler command is not their loader.

To inspect a specific bundle, start with its linked manifest, then its referenced TSV and the [seed README](../data/matdb/seed-v1/README.md). To compile a bulk bundle, the command shape is:

```bash
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/stainless-304-nist-cryogenic/manifest.tsv \
  --out /path/to/304.fsmatpk
```

Use the repository's required remote execution policy for Cargo work. For discovery, pass generated packs to the CLI with a declared request such as [stainless-thermomechanical.json](../examples/material-discovery/stainless-thermomechanical.json). Successful discovery establishes the requested data/model coverage, not physical validation of a full run.

This table was counted directly from every current manifest and its referenced source, using populated scalar/curve records and explicit property-name mappings. It was not generated by executing every compiler or solver. It is a dated working-tree snapshot: later source edits require recounting before quoting its totals as current.
