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
characteristic length [m], and prescribed relative free-stream speed [m/s].
Cylinder length means diameter; flat-plate length is the streamwise extent.
Every other v5 record is the existing v4 record, in the same order: dry contacts,
films, film surfaces, series heat, material heat budget, fiber/string heat, links.
`boundary adiabatic` remains available in v5 and has no `convection` block.
These rows are a syntax example, not measured instrument geometry.

## One physical state and one heat step

Each original nonlinear thermal trial evaluates the gas at the arithmetic
wall/bulk film temperature. The original `IsobaricGasLaw` preserves composition
and pressure, and resolves density, viscosity, conductivity and Cp together.
For humid gas, water mole fraction stays fixed rather than inventing water to
hold relative humidity constant. All original source/phase errors propagate.

The existing `fs-convection` card computes the Nusselt number and h = Nu k/L:
Churchill-Bernstein crossflow or the laminar average isothermal-flat-plate law.
No second correlation formula, fitted cooling coefficient or material-name
switch is introduced. The public `ForcedGasBoundary` result retains the
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

Radiation, natural convection, evaporation, chemistry, rarefaction, finite-gas
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
composition and passive heat/entropy signs. Native tests cover grammar,
original finite-cell energy balance, temperature-dependent pressure and
unchanged stream/chunk semantics. The accompanying patch bundle records which
checks were actually executed; authored tests are not passing-test claims.
