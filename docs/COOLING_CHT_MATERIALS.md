# Material data in `cooling-cht`

`frankensim --json cooling-cht scene.json` accepts the same gas declarations and
canonical material-card packs as the shared physics libraries. Changing the
declaration or pack changes the coefficients used by both the existing momentum
and heat solvers. Geometry, boundary conditions and solver code stay the same.

Every source declaration selects its mode explicitly. With `"mode": "frozen"`,
the full source is evaluated at its declared reference state and its coefficients
are held fixed during the run. Gases, liquids and solid thermal sources can instead use
`"mode": "temperature-dependent"` to reevaluate their original source laws at
local solved temperatures in the steady coupled flow/heat solve described below.
An evolving fluid can also drive the native unsteady flow/heat march with fixed
solid storage, as described in the transient section.
Neither mode introduces compressibility, models melting, or establishes phase
stability. Source errors, missing properties and out-of-domain reference states
refuse before solving; the evolving mode also checks every trial and endpoint.

## Atmospheres

Replace the scene's `fluid` field with, for example:

```json
{
  "mode": "frozen",
  "gas": "mixture 300 100000 mole 3 N2 0.77 O2 0.21 Ar 0.02"
}
```

The common gas reader owns temperature in kelvin, pressure in pascals,
composition normalization and the named-species model's validity limits.
Both `mole` and `mass` composition bases are supported. The legacy moist-air
record is also supported, for example `"300 100000 0.5"` for T, p and relative
humidity. These convenience models carry their existing engineering-model
limitations; they are not newly measured source data.

Density, heat capacity, conductivity and viscosity come from that same gas
state. Pressure therefore changes mass flow and kinematic viscosity. Replacing
nitrogen/oxygen with helium changes conductivity and caloric properties as
well. With frozen-property buoyancy, the gas's original reference temperature and `1/T`
expansion coefficient are used; conflicting scene overrides refuse.

Use `"mode": "temperature-dependent"` on the same gas declaration to evaluate
the actual local density, viscosity, conductivity and heat capacity throughout
the steady or unsteady buoyancy/heat solve. Gas temperature is already part of the shared
grammar, so **do not add `temperature_axis`** to a gas declaration. For example:

```json
{
  "mode": "temperature-dependent",
  "gas": "mixture 300 100000 mole 2 N2 0.8 O2 0.2"
}
```

Pressure and species fractions remain fixed. A heated humid-air record retains
its water mole fraction, so its relative humidity changes with temperature;
the model does not add water to keep relative humidity fixed. Cooling into
supersaturation refuses. Named mixtures retain their existing temperature and
pressure limits. Local density supplies `(rho(T) - rho_reference) * gravity`
to the original momentum solver, without a second expansion coefficient.

For additional species, supply the existing canonical source packs:

```json
{
  "mode": "frozen",
  "gas": "thermochemical 300 100000 mole 1 canonical-species-id 1",
  "phase": "frozen-gas",
  "sources": [
    {
      "species": "species.fsmatpk",
      "nasa9": "nasa9.fsmatpk",
      "transport": "transport.fsmcdpk"
    }
  ]
}
```

The species identity must match the original pack association. All paths are
relative to the scene. The original source loader checks species, NASA-9 and
transport associations; no species name chooses transport data. Alternatively,
`"phase": "liquid-water-limit"` requires an explicit `"water_species"` and uses
the original liquid-water admission policy. Source files are limited to
16 MiB each and 64 MiB in total per scene.

Temperature-dependent source-pack gases retain their original NASA-9 and
transport models and reevaluate them at every trial temperature. This requires
`"phase": "frozen-gas"`; a one-time `liquid-water-limit` check cannot admit
future cell states. Every original source domain remains enforced, and the
same density/heat-capacity departure limits below apply.

## Liquids

Use a canonical material-card pack with the four liquid transport properties:
density, specific heat, thermal conductivity, and the explicitly selected
dynamic or kinematic viscosity. The existing liquid resolver enforces its
property quantity conventions, including static thermal conductivity.

```json
{
  "mode": "frozen",
  "pack": "coolant.fsmcdpk",
  "viscosity": "dynamic",
  "state": [
    {"name": "T", "value": "300 K", "kind": "legacy"},
    {"name": "pressure", "value": "100000 Pa", "kind": "pressure"}
  ]
}
```

The state grammar is shared with `electrothermal`; use the coordinate names
and quantity kinds required by the source pack. Kinematic viscosity is converted
with the same sourced density. A **frozen** liquid with buoyancy requires explicit
`expansion_per_k` and `reference_temperature_k` at scene level. A liquid never
inherits the legacy gas approximation `1/T`.

## Solid thermal materials

An entry in the scene's `materials` array can use a source instead of authored
coefficients:

```json
{
  "name": "heat-spreader",
  "source": {
    "mode": "frozen",
    "pack": "spreader.fsmcdpk",
    "state": [
      {"name": "T", "value": "300 K", "kind": "legacy"}
    ]
  }
}
```

The source must contain positive isotropic thermal conductivity, density and
specific heat. The steady solve uses its conductivity; transient storage uses
the product of its density and specific heat. This path does not infer
anisotropic axes or a temperature table. An explicit `source` cannot be mixed
with `conductivity_w_m_k`, `conductivity_table`, or
`volumetric_heat_capacity_j_m3_k` overrides. Surface emissivity remains an
independent explicitly authored surface property in frozen mode.

## Coupled temperature-dependent fluid and solids

The existing steady solver closes gas or liquid and solid material feedback on
one temperature field. Select the evolving mode on the fluid and **every** solid
source. Liquids and solids name the absolute-temperature coordinate already
present in their state:

```json
{
  "fluid": {
    "mode": "temperature-dependent",
    "temperature_axis": "T",
    "pack": "coolant.fsmcdpk",
    "viscosity": "dynamic",
    "state": [
      {"name": "T", "value": "300 K", "kind": "absolute-temperature"},
      {"name": "pressure", "value": "100000 Pa", "kind": "pressure"}
    ]
  },
  "materials": [{
    "name": "heat-spreader",
    "source": {
      "mode": "temperature-dependent",
      "temperature_axis": "T",
      "pack": "spreader.fsmcdpk",
      "state": [
        {"name": "T", "value": "300 K", "kind": "absolute-temperature"}
      ]
    }
  }],
  "gravity_m_s2": [0, 0, -9.81]
}
```

These fields belong to an ordinary `frankensim.cooling-cht.v1` scene with its
existing geometry, thermal boundaries, contacts and solver controls. Coordinate
names and kinds must match each pack; explicit legacy `T` and dimensional
absolute kelvin remain supported. Each material's reference temperature supplies
its cells' initial numerical guess. Pressure, composition and other source
coordinates remain fixed throughout the solve.

The retained fluid law supplies density, heat capacity, viscosity and thermal
conductivity together. Its actual density departure from the reference drives
buoyancy; no independent expansion coefficient is used. An authored
`expansion_per_k` therefore refuses, as does a scene `reference_temperature_k`
that differs from the fluid's original reference. Omitting gravity gives the
same temperature feedback under the declared forced-flow boundaries.

Inertia and continuity retain the fluid's reference density under the existing
Boussinesq approximation. Heat advection uses reference heat capacity unless
the source-enthalpy option below is selected. Every trial and
endpoint must keep each of those two source values within **5%** of its reference.
Viscosity and conductivity use their actual local values. Solid conductivity is
reevaluated at each solid cell; the complete original source bundle, including
its other required properties, must remain supported. The source claims are
pinned: leaving their domains refuses instead of switching tables, clamping
endpoints or silently using reference constants.

The existing conservative heat operator and variable-viscosity momentum operator
solve the coupled state. Acceptance checks the unrelaxed temperature defect,
fluid and solid transport defects, and buoyancy-force defect. With
`r = min(solver.tolerance, 1e-6)`, the temperature tolerance is `r * T_reference`,
the relative transport tolerance is `r`, and the force tolerance is
`r * rho_reference * max(1, |gravity_x|, |gravity_y|, |gravity_z|)` in N/m³.
The next temperature trial is relaxed by one half; this does not relax the
acceptance checks. `solver.max_iterations` bounds the coupled passes and each
momentum solve, and the existing wall-time cancellation budget remains active.

This steady mode supports laminar gas or liquid flow and isotropic sourced solids.
Mixed frozen/authored/evolving material assignments, radiation, LVEL turbulence
and compact components refuse explicitly in the steady source-feedback solve.
The distinct unsteady path below also supports source-integrated solid heat
storage at fixed reference mass. Phase change, deformation and arbitrary
anisotropic source conductivity remain outside these CLI paths.

Results identify `solver: "fv-simplec-material"`. `material_feedback` reports
the accepted defects and tolerances, the solid-to-fluid heat transfer, and
temperature/conductivity ranges for actual liquid and solid endpoint samples.
Each endpoint `states_identity` hashes the original cell index and complete
resolved material-state identity in deterministic cell order; it summarizes all
cells rather than substituting a body-average query for the solved field.
For gases, the fluid endpoint reports temperature, density, heat capacity,
dynamic viscosity and conductivity ranges; `samples_identity` hashes the actual
numerical samples and retained thermochemical source receipts in cell order.
Source or convergence failures publish no partial result.

JSON results include `material_inputs`, with the exact fluid coefficients and
reference gas declaration, original thermochemical pack identities when used,
and material-card/state identities for source-backed liquids and solids. These
identify reference inputs; evolving endpoint evidence is reported separately in
`material_feedback`. Neither promotes estimated numerical results to
experimental validation.

## Source-law fluid enthalpy

An evolving fluid can use its original caloric law for transported heat by
adding a top-level `fluid_enthalpy` object. This also selects source-law fluid
storage in a native unsteady scene. For a liquid material pack, declare the
supported temperature span, the caloric interpolation accuracy and a local
heat-flux convergence tolerance explicitly:

```json
{
  "fluid_enthalpy": {
    "minimum_temperature_k": 290,
    "maximum_temperature_k": 315,
    "maximum_interpolation_error_j_kg": 0.001,
    "maximum_flux_defect_w": 1e-7
  }
}
```

The span must contain the fluid's original reference temperature, every inlet
and potential backflow temperature, and every thermal trial. Bounds must fit
the original source support. Liquid `Cp(T)` is integrated into the existing
bounded caloric chart from the same material pack, fixed pressure/composition
coordinates and pinned source claims. The declared interpolation tolerance
controls its specific-enthalpy approximation in J/kg; an accuracy request that
exceeds the chart's finite knot budget refuses. Density and transport still
come from the original liquid law at their native material time level.

For a gas, supply the same temperature and flux fields but **omit
`maximum_interpolation_error_j_kg`**. The original gas model supplies its
analytic enthalpy increment: the engineering gas models retain their constant
mixture Cp, while thermochemical gases use the original NASA-9 caloric source.
No formation-energy subtraction, species substitution or source reconciliation
is introduced. Original temperature domains and humid-air saturation checks
remain active.

The native heat operator transports donor `h(T)` with conservative first-order
upwind fluxes on the projected flow faces. Nonlinear acceptance checks the
source-law heat defect on every face and in every fluid cell; global
cancellation cannot conceal a local error. This option retains the existing
fixed-density Boussinesq model and its **5% source density and Cp departure
limits**, and does not admit latent phase change. The original source reference
temperature supplies a common numerical enthalpy datum.

`material_feedback.fluid_enthalpy` identifies the selected advection law,
reference temperature, declared span and tolerances, liquid chart identity
(or `null` for the original analytic gas law), and final face/cell flux defects.
Existing material input and endpoint receipts remain available.

### Exact original isobaric Cp

To integrate the original liquid Cp segments directly, add
`basis: "isobaric-cp"` to the same `fluid_enthalpy` object:

```json
{
  "fluid_enthalpy": {
    "basis": "isobaric-cp",
    "pressure_axis": "pressure",
    "minimum_temperature_k": 290,
    "maximum_temperature_k": 315,
    "maximum_flux_defect_w": 1e-7
  }
}
```

This explicit option uses the existing `IsobaricHeatCapacityLaw` and pins the
Cp claim selected by the complete original fluid state. It retains the source's
full fixed pressure/composition context and every Cp breakpoint across the
declared interval. The native caloric owner integrates those segments directly,
without an interpolated enthalpy chart or a second material fit. Omit
`maximum_interpolation_error_j_kg` with this basis: no chart approximation
tolerance applies. Omitting `basis` preserves the existing liquid chart path
and its required interpolation tolerance.

A liquid must name its original pressure coordinate with `pressure_axis`.
The source state must carry positive typed absolute temperature and pressure,
and the selected Cp claim must explicitly support that pressure over the
complete interval. An unsupported extra coordinate, legacy untyped pressure,
stress quantity or incompatible semantic kind refuses. The selected specific
heat must represent single-phase isobaric Cp at the held conditions; this
declaration does not convert Cv or an apparent phase-change capacity into Cp.

A gas may also select `basis: "isobaric-cp"`, while omitting both
`pressure_axis` and `maximum_interpolation_error_j_kg`. Its existing record
owns pressure and composition. Preparation calls the original gas owner's
whole-interval caloric admission, including capacity bounds and crossed
thermochemical region seams. An unreconciled NASA-9 enthalpy jump refuses;
this adapter does not adjust or reconcile source coefficients.

Exact liquid results report `basis: "isobaric-cp"`, held `pressure_pa`
and `caloric_identity`, which retains the prepared original source receipts.
Their `chart_identity` and interpolation tolerance are `null`.
Original material input and endpoint identities remain unchanged. Exact gas
results report their basis and held pressure alongside the existing source
evidence. The same 5% density/Cp bounds, upwind operator and nonlinear heat-flux
acceptance apply.

### Transient source enthalpy

For an unsteady scene, add explicit nonlinear and total-work controls inside
the same object:

```json
{
  "fluid_enthalpy": {
    "minimum_temperature_k": 290,
    "maximum_temperature_k": 315,
    "maximum_interpolation_error_j_kg": 0.001,
    "maximum_flux_defect_w": 1e-7,
    "transient_controls": {
      "maximum_passes_per_step": 40,
      "maximum_energy_solves": 4000,
      "temperature_tolerance_k": 1e-7,
      "maximum_storage_defect_w": 1e-6,
      "maximum_energy_closure_j": 1e-6,
      "temperature_relaxation": 1
    }
  }
}
```

These controls are required for transient source enthalpy and refuse on a
steady scene. Transient `fluid_enthalpy` supports authored or frozen solid
storage, or every solid evolving through its original source storage interval.
The latter shares one native thermal endpoint with fluid calorimetry and retains
the additional `transient.solid_storage` conductivity/storage controls.
Counts must be positive integers, with at most 1,000,000 passes
per step and 10,000,000 energy solves for the whole march. All tolerances must
be finite and nonnegative, and relaxation must be in `(0, 1]`. Choose them for
the scene's temperature, heat and energy scales; the examples are explicit
settings, not automatic error estimates.

Every accepted thermal step stores `rho_reference * V * [h(T_new)-h(T_old)]`
using the original source increment. Backward-Euler storage linearizations
all retain the same accepted old state and the same native momentum step.
Acceptance checks the unrelaxed temperature defect, each cell's storage-rate
defect, local advective defects and actual energy closure independently.
`maximum_energy_solves` counts every linear heat solve across all steps,
including inner endpoint-radiation trials. Exhaustion or nonconvergence
publishes no partial march.

The transient caloric report includes actual `energy_solves`, the largest
accepted per-step work and defects, and all declared limits. Each emitted
transient history record includes `stored_energy_change_j`, `source_energy_j`,
`boundary_outflow_j` and `closure_j`, with positive boundary energy directed
outward. These are the native physical heat records, including thermostat and
radiation heat. The existing history sampling limit still applies.

When evolving solid storage is also selected, each accepted native heat step
contains both the fluid's original enthalpy increment and each solid's
fixed-reference-mass Cp integral. Solid conductivity is evaluated at that same
thermal endpoint. Native interface/contact fluxes exchange equal and opposite
heat; they do not create a second fluid/solid heat budget. The ordinary
`stored_energy_change_j`, `source_energy_j`, `boundary_outflow_j` and
`closure_j` records retain the combined physical ledger.
`material_feedback.solid_storage` continues to report original source
identities, solid endpoint defects and solid stored-energy increments.
The combined report identifies `solid_coupling: "joint-implicit-endpoint"`.
The fluid caloric report keeps its explicit per-step and whole-march work caps;
a source, acceptance or budget failure publishes no partial result.


## Unsteady flow with temperature-dependent fluid properties

An evolving gas or liquid can now use the existing native transient momentum
and heat-storage calculation. Add this block to the scene:

```json
{
  "transient": {
    "flow": "unsteady",
    "time_step_s": 0.025,
    "steps": 40,
    "scheme": "bdf2",
    "initial_temperature_k": 300,
    "inner_iterations": 200,
    "inner_tolerance": 1e-8
  }
}
```

Select `mode: "temperature-dependent"` on the fluid. Solids can retain authored
or explicitly `mode: "frozen"` thermal coefficients and storage; frozen packs
use their original density-times-specific-heat storage. Alternatively, put every
solid source in temperature-dependent mode and declare its supported storage
temperature range as described below. That path integrates the original heat
capacity and solves solid conductivity at each accepted heat endpoint. Authored
`conductivity_table` entries still refuse in these material-fluid paths; they
cannot silently become constant or replace a sourced solid response.

At each accepted temperature, the original fluid source resolves local density,
viscosity, heat capacity and conductivity together. The next native momentum
step uses `(rho(T) - rho_reference) * gravity` and local viscosity; the heat
step uses the same projected volume fluxes, local conductivity and reference
`rho * Cp` storage/advection unless `fluid_enthalpy` selects the source-law
storage and transport above. Both source density and Cp must remain within
5% of their original reference at the initial state and every endpoint.
Pressure, composition, source choices and source domains retain the same rules
as in the steady path.

This retains **first-order fluid-material/thermal coupling**, even when the
momentum scheme is BDF2. The viscosity update preserves native velocity,
pressure and BDF history. This coupling order is reported separately from the
requested momentum scheme; time refinement is still required to assess error.
The original inlet and power schedules, material contacts, internal fans,
hydraulic resistances and optional surface-radiation owner remain available.
The existing transient face-fan and compact-component restrictions still apply.

The fluid law requires `transient.flow: "unsteady"` and marched energy.
Frozen-flow transient energy, `energy: "steady-on-mean-flow"`, LVEL and SGS
cannot represent this evolving material path and refuse. An authored
`expansion_per_k` also refuses, because actual source density owns buoyancy.
For an enclosed scene without another initial temperature, the declared fluid
reference temperature supplies the initial condition.

Results identify `solver: "fv-simplec-unsteady-material"` and retain the native
flow-step records, thermal storage closure, probes, time averages and VTK
output. `material_feedback` reports `coupling: "explicit-first-order"`,
`applied_time_s` and `endpoint_time_s`. `fluid_applied` describes the source
samples used during the final step; `fluid_endpoint` describes the checked
sources at its accepted final temperature. They are deliberately distinct
epochs. Source, departure-bound, convergence or cancellation failure returns
no partial march.

These additions retain the existing laminar, fixed-geometry Boussinesq
approximation. They do not supply compressible gas dynamics, evolving mass,
latent heat or experimental validation.

### Temperature-dependent solid heat storage

An evolving solid source supplies the same density, heat capacity and
conductivity bundle used by the steady path, plus an explicit single-phase
temperature interval:

```json
{
  "name": "solid",
  "source": {
    "mode": "temperature-dependent",
    "temperature_axis": "T",
    "storage_temperature_range_k": [290, 350],
    "pack": "solid.fsmcdpk",
    "state": [
      {"name": "T", "value": "300 K", "kind": "absolute-temperature"},
      {"name": "pressure", "value": "100 kPa", "kind": "pressure"}
    ]
  }
}
```

The selected source must support the complete declared interval. Preparation
retains the original heat-capacity claim and its breakpoints. Each voxel stores
`rho_reference * volume * integral(Cp(T) dT)`; density changes cannot manufacture
or remove solid mass in a fixed voxel. The source integrator evaluates every
crossed segment, so the stored-energy change does not substitute endpoint Cp
or an arbitrary temperature average. The reference density comes from the same
original solid snapshot as the other thermal properties.

The native backward-Euler heat equation iterates solid conductivity and the
source-integrated storage response together, retaining the original fluid
fluxes, thermal contacts and radiation calculation. Optional
`transient.solid_storage` sets its acceptance controls independently of the
momentum solve. The defaults are:

```json
{
  "maximum_iterations": 100,
  "temperature_tolerance_k": 1e-9,
  "conductivity_relative_tolerance": 1e-9,
  "storage_tolerance_j": 1e-8,
  "temperature_relaxation": 1.0
}
```

Acceptance checks the unrelaxed temperature and conductivity defects, the true
per-cell caloric-storage defect, and the corrected whole-step heat balance.
`material_feedback.solid_storage` reports the controls, original storage
identities, reference densities, declared intervals, and accepted step records.
`solid_endpoint` retains the actual endpoint material samples;
`solid_coupling: "implicit-endpoint"` distinguishes this solve from the explicit
fluid coefficient update. The overall time discretization remains first order.

All solids must use this evolving source path when any one does; mixed
evolving/frozen solid storage is refused. A missing caloric interval, unsupported
source state, failed nonlinear acceptance or cancellation publishes no partial
march. The model supplies single-phase sensible storage in fixed reference
geometry, without latent heat, phase inventory changes, deformation or inferred
mechanical work.
