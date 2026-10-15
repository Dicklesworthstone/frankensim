# FrankenTextiles materials in FrankenSim

## What is shared, and what is deliberately not duplicated

This integration was based on the actual FrankenTextiles implementation at
`d9779504623f7bd8b7fd01cc70c5f6270bfac6e6`, not merely its design aspirations.
Relevant source owners are `ftx-pattern/src/materials.rs`, `material_bending.rs`,
`material_response.rs`, `material_thermal.rs`, `material_catalog.rs`, and
`ftx-drape/src/material.rs`; the governing plan is
`COMPREHENSIVE_PLAN_FOR_THE_DESIGN_OF_FRANKENTEXTILES.md` in that repository.

The reusable distinction is **characterized sheet versus constituent fiber**.
A fabric's areal mass, membrane stiffness, directional bending, friction and
thermal response are independently observed quantities. Fiber composition is
not a license to derive all of them. A fiber modulus in Pa is not an effective
membrane modulus in N/m; a coupon rigidity in N m2 is not a sheet bending
coefficient in N m. We import source data into the existing immutable material
system and lower it into existing physical operators. No second material
database, DKT solver, eigensolver, caloric chart, or heat integrator is introduced.

| FrankenTextiles capability | FrankenSim integration / boundary |
|---|---|
| Canonical catalog records with methods, sources, specimen and conditioning | Offline adapter to the existing `matdb-pack` compiler; exact input bytes and qualifiers retained |
| Effective warp/weft/shear membrane constants | `fs_material::sheet::QuadraticMembrane`: Green-strain energy, work-conjugate stress and consistent tangent |
| Independently characterized directional bending | `ResolvedSheet::bending` and `fs_couple::sheet_plate::SheetPlate` feed D directly to the original DKT and modal owners |
| Per-piece grain | Explicit per-reference-triangle rotation, preserving engineering shear/twist and energy |
| Temperature-dependent conductivity and heat capacity | Separate source selection; `SheetCaloric` uses the existing Cp-to-enthalpy constructor; conductivity is re-evaluated at actual thermal state |
| Nonlinear membrane tables, locking, hysteresis and relaxation | Retained in the raw input, not silently flattened. Mechanical import refuses until a matching law is selected/implemented |
| Nonzero preferred curvature | Not a flat-plate substitution; explicit refusal for the bending import |
| Skin/self friction | Pair-specific keys preserved; not assigned to arbitrary contact partners |
| Air permeability / vapor resistance | Test-property units preserved; not interpreted as Darcy permeability, diffusion or a moisture model |
| Complete cloth/drape/contact solver | Not vendored: its dependency graph and separate evolving histories would duplicate existing owners |

The review also found that FTX's response-curve implementation explicitly
credits a prior FrankenSim response-curve owner. Copying the entire module back
would add a competing implementation, not useful integration. Its richer signed
curve/domain semantics remain useful for a later nonlinear sheet consumer.

This work advances the effective-sheet part of `frankensim-7sga6.5.7`; it does
not close that entire material-system bead or claim complete cotton, polyester,
felt and leather constitutive families.

## Import a canonical record, with explicit source authority

`scripts/import_franken_textiles.py` uses only Python's standard library during
offline source preparation. Python and FrankenTextiles are **not simulation
runtime dependencies**. Input may be `ftx-textile-catalog/1` or one canonical
Material record. Legacy/snapshot formats must first be resolved by their FTX
owner. Material selection is by exact ID, with optional source SHA-256 pinning.

The following synthetic fixture exercises the complete pipeline without
pretending to be a measured commercial fabric:

```sh
python3 scripts/import_franken_textiles.py \
  data/matdb/textile-examples/synthetic-sheet.json \
  synthetic-independent-sheet /tmp/textile-sheet \
  --citation 'Synthetic FrankenTextiles-format integration specimen' \
  --license 'Synthetic numerical fixture; repository license applies' \
  --redistribution-terms 'Authored synthetic observations may be redistributed' \
  --allow-estimates --interpolate-temperature

cargo run -p xtask -- matdb-pack \
  --manifest /tmp/textile-sheet/manifest.tsv \
  --out /tmp/textile-sheet/material.pack

cargo run -p fs-couple --example textile_sheet -- \
  /tmp/textile-sheet/material.pack 293.15 333.15 5 0.3
```

The output directory must not exist. The adapter validates before creating
files and writes the manifest last. It retains `original.json`, `selection.json`
and the native source TSV. Exact original JSON is additionally bound into the
source envelope, so methods, conditioning, SD, composition and unused fields
remain tied to the compiled pack's source identity. Hash equality establishes
byte identity, not empirical truth or permission to redistribute.

The source compiler remains the only pack compiler. Every import requires
caller-supplied citation, license and redistribution terms. An imported
`estimate` requires `--allow-estimates`. Published/estimated labels are not
promoted to measured observation backing. Means, bounds and sample SD remain
available in source metadata, but SD is never converted into a confidence
interval or standard error without a sample-size/statistical model.

### Real partial-source records stay partial

For example, FTX's
`crates/ftx-pattern/src/material_catalog/thermal-measurements.json` contains
`fsri.cotton-batting.dried.2026-10-14`: through-thickness conductivity at two
reported temperatures, without a matching dried heat-capacity curve or
qualified sheet mass/geometry. The adapter accepts its conductivity without
manufacturing the missing capabilities:

```sh
python3 scripts/import_franken_textiles.py \
  ../franken_textiles/crates/ftx-pattern/src/material_catalog/thermal-measurements.json \
  fsri.cotton-batting.dried.2026-10-14 /tmp/dried-batting \
  --properties thermal_conductivity \
  --citation 'FTX dated transcription of FSRI Cotton_Batting dried HFM table' \
  --license 'Attributed factual observations; no page or image reproduction' \
  --redistribution-terms 'Only attributed factual observations are retained'

cargo run -p xtask -- matdb-pack \
  --manifest /tmp/dried-batting/manifest.tsv --out /tmp/dried-batting/material.pack
cargo run -p fs-couple --example textile_sheet -- \
  /tmp/dried-batting/material.pack --conductivity-at 288.15
```

This reads the exact reported 15 C knot. It does not authorize interpolation,
transient storage, garments, or isotropic conductivity. Rights and factual
provenance must be checked for the exact selected input in an actual deployment.

## Source and numerical semantics

Frozen scalar properties require the explicit `ftx_characterization=0` query
coordinate. This names the retained record's fixed construction/conditioning;
it is not a physical temperature or humidity range. Temperature tables use
Kelvin, retain their actual bounded support, and default to **exact knots only**.
`--interpolate-temperature` explicitly admits a linear estimate and additionally
requires `--allow-estimates`; missing data and extrapolation remain errors.
A selected thermal response and supplied companion areal mass/thickness must
name the same specimen and conditioning. Polynomial thermal fits are not sampled
silently. Nonlinear/history-dependent mechanical records cannot be converted
into equilibrium scalar models by ignoring their richer payload.

`ResolvedSheet` retains the entire original `ResolvedMaterialStatePoint`, not
just mechanical coefficients. Consumers independently request complete membrane,
bending or thermal capabilities. Missing D12 or D66 is not set to a convenient
Poisson/shear default. Bending matrices must be positive definite. Material
rotation uses engineering curvature `[k11,k22,2*k12]`; membrane strain uses
`[E11,E22,2*E12]` and work-conjugate `[S11,S22,S12]` resultants.

The plate adapter uses `density = areal_mass / thickness` only to represent
reference inertia in the existing DKT data layout. The conserved mass is
`areal_mass * reference_area`; it is not a new fiber-density law. The original
uniform-through-thickness rotary inertia approximation remains explicit.
Uniform pretension is a scenario load, not inferred from a material label.

`SheetCaloric` pins the source Cp claim and retains all conditioning coordinates
while changing only temperature. It includes every interior Cp source knot,
then uses `EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity` and the existing
single-phase temperature inversion. Constant or admitted Kelvin-linear Cp is
required. Its interpolation allowance is separate from source uncertainty and
time error. Conductivity can be used independently without Cp.

## What the executable actually computes

The example constructs the actual DKT eigenproblem of a small clamped reference
sheet and reports its native modes, reference mass and membrane strain energy.
Its coarse 2x2 mesh is an integration demonstration, not a converged textile
frequency prediction. Bending/mass/grain changes alter the original operators;
no fabricated equivalent Young's modulus or gain-normalized sound is used.

The same material bundle supplies a two-volume through-thickness slab with
reference area 0.01 m2, each volume containing half the sheet mass and thickness.
The outer faces contact explicit hot/cold reservoirs with film coefficient
10 W/(m2 K); edges are adiabatic. The interface uses the two half-cell Fourier
resistances with each side's actual source-qualified k(T). The original thermal
kernel advances enthalpy using source Cp and its existing conservation and CFL
checks, with a declared maximum step of 0.01 seconds. The example checks aggregate
energy before each CSV row. Thermal chart support and time/spatial refinement
remain separate obligations. Earlier stdout rows may exist if a later step fails.

This is **not** a large-deformation cloth simulation, stress-dependent porosity,
moisture transport, textile compression, full thermomechanical work balance,
thermal modal-history remapping, empirical timbre model or arbitrary fiber-mixture
homogenizer. Such capabilities require their own qualified source laws and common
history coupling, rather than guessed conversions. All dynamics here reuse
FrankenSim's existing solvers.
