# Material data in `cooling-cht`

`frankensim --json cooling-cht scene.json` accepts the same gas declarations and
canonical material-card packs as the shared physics libraries. Changing the
declaration or pack changes the coefficients used by both the existing momentum
and heat solvers. Geometry, boundary conditions and solver code stay the same.

Every source declaration selects its mode explicitly. With `"mode": "frozen"`,
the full source is evaluated at its declared reference state and its coefficients
are held fixed during the run. Liquids and solid thermal sources can instead use
`"mode": "temperature-dependent"` to reevaluate their original source laws at
local solved temperatures in the steady coupled flow/heat solve described below.
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
well. With buoyancy, the gas's original reference temperature and `1/T`
expansion coefficient are used; conflicting scene overrides refuse.

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

## Coupled temperature-dependent liquid and solids

The existing steady solver can close liquid and solid material feedback on one
temperature field. Select the evolving mode on the liquid and **every** solid
source, naming the absolute-temperature coordinate already present in its state:

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

The retained liquid law supplies density, heat capacity, viscosity and thermal
conductivity together. Its actual density departure from the reference drives
buoyancy; no independent expansion coefficient is used. An authored
`expansion_per_k` therefore refuses, as does a scene `reference_temperature_k`
that differs from the liquid's original reference. Omitting gravity gives the
same temperature feedback under the declared forced-flow boundaries.

Inertia, continuity and heat advection retain the liquid's reference density and
heat capacity under the existing Boussinesq approximation. Every trial and
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

This mode supports steady laminar liquid flow and isotropic sourced solids.
Gas evolution, mixed frozen/authored/evolving material assignments, transient
storage, radiation, LVEL turbulence and compact components refuse explicitly.
Existing frozen scenes retain their original supported features. Solid storage,
phase change, deformation and arbitrary anisotropic conductivity are outside
this evolving CLI path.

Results identify `solver: "fv-simplec-material"`. `material_feedback` reports
the accepted defects and tolerances, the solid-to-fluid heat transfer, and
temperature/conductivity ranges for actual liquid and solid endpoint samples.
Each endpoint `states_identity` hashes the original cell index and complete
resolved material-state identity in deterministic cell order; it summarizes all
cells rather than substituting a body-average query for the solved field.
Source or convergence failures publish no partial result.

JSON results include `material_inputs`, with the exact fluid coefficients and
reference gas declaration, original thermochemical pack identities when used,
and material-card/state identities for source-backed liquids and solids. These
identify reference inputs; evolving endpoint evidence is reported separately in
`material_feedback`. Neither promotes estimated numerical results to
experimental validation.
