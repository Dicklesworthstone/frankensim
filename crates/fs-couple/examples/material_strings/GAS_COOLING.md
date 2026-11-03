# Source-gas convection on the native string thermal timeline

A v5 contact-thermal profile may explicitly select `boundary gas-reservoir`.
Existing v1-v4 adiabatic profiles retain their interpretation. The instrument's
original `gas_record` supplies reservoir temperature, pressure, composition
and the existing temperature-dependent transport law. This is a prescribed
infinite reservoir, not a finite room-air energy or flow simulation.

After the last `cell` row, and before `dry-contacts`, add:

```text
convection 2
convection-cell 10 cylinder-crossflow 0.000001 0.0003 1
convection-cell 20 laminar-flat-plate 0.000001 0.002 1
```

The columns are cell ID, exterior geometry/regime, exposed area [m2],
characteristic length [m], and a geometry-specific driving input. Forced
geometries use prescribed relative free-stream speed [m/s]; natural geometries
use explicit gravity magnitude [m/s2] in quiescent gas. The aiding mixed plate
requires three driving fields: signed upward speed [m/s], gravity [m/s2], and
the caller-admitted maximum relative density departure. For example:

```text
convection 1
convection-cell 10 natural-sphere 0.000001 0.001 9.81
```

```text
convection 1
convection-cell 10 aiding-vertical-plate 0.000001 0.002 0.2 9.81 0.1
```

| Geometry | Characteristic length | Driving fields |
| --- | --- | --- |
| `cylinder-crossflow` | Cylinder diameter | Positive speed |
| `laminar-flat-plate` | Streamwise extent | Positive speed |
| `turbulent-flat-plate` | Streamwise extent | Positive speed |
| `sphere` | Sphere diameter | Positive speed |
| `natural-vertical-plate` | Plate height | Nonnegative gravity |
| `natural-horizontal-cylinder` | Cylinder diameter | Nonnegative gravity |
| `natural-horizontal-plate-up` | Full plate area/perimeter | Nonnegative gravity |
| `natural-horizontal-plate-down` | Full plate area/perimeter | Nonnegative gravity |
| `natural-sphere` | Sphere diameter | Nonnegative gravity |
| `aiding-vertical-plate` | Vertical plate height | Signed upward speed, nonnegative gravity, density-departure bound in [0,1) |

Both sphere geometries admit exposed area only up to pi D². Horizontal face
orientation is relative to gravity; contact-masked exposed area must not replace
the FULL plate area in its area/perimeter characteristic length. The turbulent
plate uses the existing average relation including its laminar leading-edge
correction, with its original Re/Pr domain. It does not infer a transition from
surface roughness or free-stream turbulence. Zero speed never changes a
forced-only card into a natural one; the geometry declares which model is
intended. The mixed plate explicitly retains the pure natural limit as its
speed reaches zero and the pure forced limit when gravity is zero.
Every other v5 record is the existing v4 record, in the same order: dry contacts,
films, film surfaces, series heat, material heat budget, fiber/string heat, links.
`boundary adiabatic` remains available in v5 and has no `convection` block.
These rows are a syntax example, not measured instrument geometry.

## One physical state and one heat step

Each original nonlinear thermal trial evaluates forced plate/cylinder, natural and mixed-convection properties at
the arithmetic wall/bulk film temperature. The original `IsobaricGasLaw` preserves composition
and pressure, and resolves density, viscosity, conductivity and Cp together.
For humid gas, water mole fraction stays fixed rather than inventing water to
hold relative humidity constant. All original source/phase errors propagate.

The existing `fs-convection` card computes the Nusselt number and h = Nu k/L:
Churchill-Bernstein crossflow, the laminar/turbulent average flat-plate laws,
the existing oriented natural-convection laws, Manning-Qureshi aiding mixed
convection, or Whitaker's forced-sphere law.
The forced sphere uses the original bulk gas for
Re, Pr and k, and queries the same gas at wall temperature only for its
bulk/surface viscosity correction. `GasConvectionSample::property_gas` retains
the coefficient state; `surface_gas` retains the forced sphere's second query.

Natural convection also obtains thermal expansion from that original gas's
thermal transport owner. The signed full wall/bulk temperature difference
selects stable versus unstable horizontal convection. The existing geometry
selector owns every Rayleigh/Prandtl bound and its empirical branch switches.
At exact equilibrium, wall heat and entropy production are zero. For a natural
plate/cylinder this is represented by absent `nusselt` and `coefficient`, without
extrapolating a Ra=0 card; it makes no claim about the heat derivative there.
Nonzero driving outside a card's range still refuses. The natural sphere retains
its sourced Nu=2 conductive limit, including at zero gravity; no parallel
stagnant-gas conductance should be added to it.

The aiding mixed plate uses the original film gas's expansion, Reynolds,
Rayleigh and Prandtl numbers with the full wall/bulk contrast. A positive speed
must aid upward buoyancy at a hotter wall; a colder wall needs downward speed.
The shared `fs-material::fluid_transport` admission checks the actual signed
expansion and requires `|beta (T_wall - T_bulk)|` not to exceed the authored
bound. There is no default bound. The complete authored temperature interval
must stay on the admitted side of the reservoir for a fixed nonzero stream;
each subsequent thermal trial is checked again. Opposing and crossflow cases
need their own model and cannot use this vertical aiding row.

Manning-Qureshi supplies an average laminar isothermal-plate relation with
strict Re < 500000 and Ra < 1000000000 bounds. The shared card's Pr interval
[0.003,100] is an implementation admission policy spanning its source comparison
points, not a continuously validated experimental interval. Its theory assumes
constant properties; mean-film sampling is this adapter's declared approximation,
and a small density departure alone does not establish weak variation of every
transport property. The native route retains the same card in its pure limits,
without switching to unrelated flat-plate or natural-convection formulas.
The original relation and evidence remain in `fs-convection`:
Manning and Qureshi (1994), Eqs. (4), (8), (10), (11),
https://doi.org/10.1115/1.2910940 .

At exact resting mixed equilibrium, the shared sample has zero heat and entropy
with absent Nu and h. With a nonzero stream at equilibrium, its finite forced
coefficient remains available. If the wall differs from the reservoir while
both drivers vanish, this card refuses; it supplies no invented diffusion floor.
World-frame library callers use `AidingVerticalPlateFrame` to reject inclined
gravity or oblique flow before constructing the same signed settings.

The sphere card requires Re in [3.5, 76000], Pr in [0.71, 380] and
mu_bulk/mu_surface in [1, 3.2]. A hotter sphere in a colder ordinary gas can
violate the viscosity range, and the USSA dry-air fit at 300 K has Pr below
0.71. Both cases refuse without clamping; a source-backed nitrogen atmosphere
and a colder sphere provide an admitted example. The entire thermal-source
temperature interval must remain admissible during profile preparation.
These limits come from Whitaker (1972), Eq. (33), Table 1,
https://doi.org/10.1002/aic.690180219 .

No second correlation formula, fitted cooling coefficient or material-name
switch is introduced. `GasConvectionBoundary` explicitly selects a
`ForcedGasBoundary`, `NaturalGasBoundary` or `AidingVerticalPlateGasBoundary`;
their shared sample retains the
actual gas evaluation, correlation evidence, signed heat and instantaneous
entropy production. These are correlation-rung models, not resolved CFD or
experimental qualification of an instrument.

The full temperature drop supplies Q_wall = h A (T_bulk - T_wall). It enters the
existing `StepCoefficients.sources_w` once as external heat; -Q_wall is heat to
the prescribed reservoir. Original conduction, dry/film heat, caloric storage,
material-temperature updates, motion and PCM keep their existing owners.
The film-surface callback starts from this same source map, so adding retained
film transfer does not discard the environmental term. Rejected trials change
neither the accepted thermal state nor any external reservoir state.
The native caloric owner records external heat; this adapter does not add a
separate time-integrated reservoir entropy ledger. The shared sample reports
instantaneous boundary entropy production only.

## Physical boundaries

Exposure is explicit, fixed, and unique per cell. It cannot exceed the cell's
original surface. The author must exclude contact-covered area and choose
geometry/regime consistent with the real object. There is no automatic
moving-contact exposure mask or replacement of actual bow speed by a mean wind.
Correlation domain failures refuse instead of selecting a different regime.

The native profile currently supports legacy moist-air and named-mixture gas
records. A thermochemical snapshot alone cannot reconstruct its source model,
so the frontend refuses that case. The shared boundary API accepts an original
`IsobaricGasLaw::try_from_model` law when a caller retains that model. It never
turns a thermochemical snapshot into a constant-property gas.

Radiation, opposing or turbulent mixed convection, evaporation, chemistry, rarefaction, finite-gas
energy and temperature-dependent elastic storage are not supplied by this path.
The mechanical approximation and source-temperature splitting remain those of
the original performer. Time, spatial and material-model validation remain
necessary; this extension does not establish violin/cello sound fidelity.

## Focused native checks

```sh
cargo test --locked -p fs-couple --lib air_path::input::convection::tests
cargo test --locked -p fs-couple --example material_strings gas_convection
cargo check --locked -p fs-couple --bin material_strings --example material_strings
```

The shared tests compare the original source evaluation against independent
correlation arithmetic, pressure scaling, mixture substitution, humid-gas
composition, natural face orientation, mixed pure limits and density/alignment
refusal, equilibrium and passive heat/entropy signs. Native tests cover
forced/natural/mixed grammar, complete source-interval admission,
original finite-cell energy balance, temperature-dependent pressure and
unchanged stream/chunk semantics. The accompanying patch bundle records which
checks were actually executed; authored tests are not passing-test claims.
