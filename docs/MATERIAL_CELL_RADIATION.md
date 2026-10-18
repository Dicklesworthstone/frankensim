# Material radiation on finite thermal and circuit histories

`fs_conduction::transient::deposition::radiation::MaterialRadiation` connects
an existing gray-diffuse enclosure to the shared `fs-thermal` cell history.
Each surface names its existing receiving cell, original caloric curve,
material card and complete physical condition point. It does not introduce
another temperature store, radiation equation or circuit solver.

## One material state for every participating property

`SurfaceEmissivity::from_thermal_state` resolves the canonical
`hemispherical-total-emissivity` property through the common material-state
resolver. Temperature and the explicitly declared liquid-mass fraction come
from the original enthalpy chart. Pressure, lot and other supplied conditions,
quantity descriptors, property selections and source receipts remain intact.
A phase-dependent emissivity may change on a constant-temperature latent
plateau; it does not need an invented Kelvin increment.

`ThermalMaterialPoint` now lives in `fs_material::thermal::source`. The original
`fs_couple::cavity_phs::piston::thermal::source` path re-exports those exact types,
so electrical heaters, conductive links and radiation share the same condition
update and static-law admission. The original temperature-only emissivity
constructors remain available as convenience paths, not substitutes for the
full condition point.

Every radiation trial replaces the complete emissivity list on the original
`GrayDiffuseEnclosure`, preserving its mesh traces, areas, names and view-factor
matrix. The existing radiosity solver evaluates the resulting operator. An old
prepared factorization remains tied to its old emissivities and is never
silently reused for changed material coefficients.

## Finite-cell use

Create one `RadiatingCell` binding per enclosure surface. Surface indices follow
matrix order; receiving `CellId`s need not. Bind the same borrowed caloric object
used to initialize the thermal timeline, not a separate curve that happens to
agree at one temperature. Multiple distinct traces may feed one uniform cell;
their summed area must fit its exposed area after thermal contacts and reservoir
coverage. These checks do not prove the geometric trace-to-cell correspondence
or physical non-overlap; those remain explicit modeling obligations.

```rust,ignore
use fs_conduction::transient::deposition::radiation::MaterialRadiation;

let radiation = MaterialRadiation {
    enclosure: &existing_enclosure,
    cells: &radiating_cells,
    power_tolerance_w: 1e-8,
    radiosity_tolerance_w_m2: 1e-8,
};
let (next_thermal, report) = radiation.step_candidate(
    cx, &thermal, target_s, controls, &cell_energy_guess_j,
)?;
// Radiation is already included in next_thermal. Do not apply it again.
```

The standalone interval uses the existing backward-Euler cell-energy residual,
Newton--Krylov driver and fallible finite-difference adapter. Each trial starts
from the same accepted enthalpy history. Conductive links, authored heat,
convection and reservoirs remain with the original thermal owner. Radiation
and convection may act on the same exposed surface because they are different
heat-transfer mechanisms. Explicit events still delimit physical intervals:
process start-time events first, split at interior events, and process target
events after accepting the returned left-limit state.

Fresh final material/radiosity evaluation, original radiosity residual and
closed-enclosure power checks, local cell energy gates, and total storage versus
nonradiative external heat all precede acceptance. Every source evaluation,
including numerical probes and the final check, consumes the declared budget.
Cancellation, unsupported material, exhausted work or failed balance returns
no thermal endpoint. Residual allowances are not discretization-error bounds.

## Compose with the existing finite circuit network

The same source result can be added to the existing circuit/enthalpy residual
without introducing another coupled time stepper:

```rust,ignore
// In circuit_enthalpy::step_candidate's existing material callback:
let mut coefficients = material_network.resolve(trial_thermal, cancelled)?;
let emission = radiation.resolve(cx, trial_thermal)?;
emission.add_to(&mut coefficients.sources_w)?;
// Return these coefficients and retain emission in the concrete source payload.
```

The radiation source adds signed watts to the correct cells without replacing
independent sources. The circuit driver alone deposits Joule heat and advances
capacitor charge/inductor flux. All owners evaluate against the same trial
enthalpy. The cross-crate regression changes only an emissivity value and checks
that both temperatures and finite electrical-storage history change while total
energy remains balanced. No radiation watts should also be added as external
electrical work.

## Physical scope and focused checks

This is a closed opaque gray-diffuse enclosure with fixed geometry/view factors
and uniform receiving thermal cells. A source must supply static total emissivity
in `(0,1]`; missing data, dynamic-only laws and unsupported conditions refuse.
It is not a spectral optical model, participating-gas radiation, a moving mesh,
a phase-front solver, or an automatic visible-color/emissivity conversion.
The enclosure's original view-factor evidence remains authoritative; synthetic
infinite-plate fixtures do not qualify finite-cube geometric view factors.

Thermal rates and material coefficients use endpoint evaluation, so the coupled
method is generally first order. Check time resolution and the numerical
Jacobian scale independently; energy conservation and algebraic convergence do
not establish uniqueness, passivity or measured-material accuracy.

```sh
cargo test -p fs-conduction --test radiation_material_state --test deposition_radiation -- --test-threads=1
cargo test -p fs-couple --test circuit_material_radiation -- --test-threads=1
cargo test -p fs-couple --lib cavity_phs::piston::thermal:: -- --test-threads=1
```

Native Rust execution was unavailable during authoring. Independent Python
radiosity, implicit cell and circuit/thermal equations supply numerical reference
checks, not native test passes or measured-material qualification.
