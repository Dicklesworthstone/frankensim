# Material data in `cooling-cht`

`frankensim --json cooling-cht scene.json` accepts the same gas declarations and
canonical material-card packs as the shared physics libraries. Changing the
declaration or pack changes the coefficients used by both the existing momentum
and heat solvers. Geometry, boundary conditions and solver code stay the same.

Every source declaration currently requires `"mode": "frozen"`: the full source
is evaluated at the explicitly declared reference state, then its coefficients
are held fixed during this CHT run. This is a constant-property approximation;
it does **not** refresh the material at local solved temperatures, introduce
compressibility, model melting, or establish that a different phase is stable.
Use the evolving material/enthalpy APIs for those calculations. Source errors,
missing properties and out-of-domain reference states refuse before solving.

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
with the same sourced density. A liquid with buoyancy requires explicit
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
independent explicitly authored surface property.

JSON results include `material_inputs`, with the exact fluid coefficients and
reference gas declaration, original thermochemical pack identities when used,
and material-card/state identities for source-backed liquids and solids. These
identify reference inputs; they do not promote the run's estimated numerical
results to experimental validation.
