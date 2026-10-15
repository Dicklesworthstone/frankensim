# Surface material and external heat exchange in recirculating cooling

The existing executable accepts `--radiation SURFACE.txt` before `--solid`.
The solid, coolant, pump and surface finish are separate data inputs to their
existing owners, not a new combined material database or a substitute solver.
The previous commands and CSV columns remain unchanged without this option.

```sh
cargo run -p fs-airflow --example recirculating_cooling -- \
  --radiation crates/fs-airflow/examples/recirculating_cooling/surface_low_emissivity.txt \
  --solid crates/fs-airflow/examples/recirculating_cooling/solid_low_capacity.txt \
  crates/fs-airflow/examples/pressure_cooling/constant.txt \
  crates/fs-airflow/examples/pressure_cooling/pump.txt \
  0.1 1 2 20 0 299 302 20 80
```

Replace only `surface_low_emissivity.txt` with `surface_high_emissivity.txt` to
change the declared emissivity from 0.1 to 0.9. The two descriptions are otherwise
identical. They are synthetic fixtures, not measured surface finishes. Source
strings and licenses are retained but do not authenticate a measurement.

The same 0.06 x 0.10 x 0.01 m solid retains two top-face coolant patches. The
new external patch occupies its bottom face (0.006 m2), with independently
specified convection to air and radiation to a black isothermal reservoir.
Other faces remain adiabatic. The example adds no obstruction or coating
thickness; a finish changes emissivity only, not bulk heat capacity or k(T).

## Surface description

Required `key=value` fields are `material`, `source`, `license`, `emissivity`,
`minimum_surface_k`, `maximum_surface_k`, `convection_w_m2_k`,
`air_temperature_k`, and `radiation_temperature_k`. `#` begins a comment.
The importer rejects duplicates, unknown fields, missing data, invalid values
and descriptions above 128 KiB. It creates an ordinary material card containing
one constant scalar hemispherical-total-emissivity claim and resolves that
claim through `SurfaceEmissivity`. Emissivity must lie in (0,1].

The temperature span must include the initial 300 K surface and every actual
surface state admitted by the existing radiation owner. Environmental reservoir
temperatures need not belong to that surface-property span. The explicit air
coefficient is positive; zero does not get replaced with a fabricated small
value. The present ambient-patch owner requires an existing positive uniform
Robin convection row. Pure-radiation/vacuum boundaries are not added here.

## One time interval and separate heat destinations

`advance_stored_solid_with_pump_and_radiation` calls the existing implicit
`EnthalpyBackwardEuler::advance_with_ambient_radiation` at every solid response.
Every radiation, pump-pressure, coolant and reservoir trial uses the same old
solid enthalpy and physical duration. The original radiative solver checks the
freshly reassembled physical nodal balance and nonlinear energy residual.

Only original convection on coolant regions enters the returning coolant.
Radiation never heats that coolant directly. Original convection on the named
external patch belongs to its separately declared air reservoir. The joint
storage balance is therefore

    delta E_solid + delta H_tank
      = dt * (body_power - Neumann_out - external_convection_out - radiation_out)
        + explicit_tank_heat

Radiative and convective losses may have opposite signs. Hot radiative
surroundings can supply heat even when body power is zero. The implementation
uses the accepted nonlinear radiative heat, not a stale secant estimate, for
this balance. Both histories publish together only after all existing gates
and a final cancellation check. Invalid or unowned extra Robin regions refuse.

The CSV appends `radiation_out_j`, `external_convection_out_j`,
`physical_solid_residual_j`, and `enthalpy_solves`. Positive external heat is
outward, negative heats the solid. Existing `solid_heat_w` remains coolant-only
heat; `solid_external_heat_j` now subtracts both external boundary losses.
Complete applied boundaries, separated original traces, emissivity receipts,
physical residuals and cumulative Newton/Krylov work remain on library results.

The original per-step pump/branch work limits remain active. A new 65536 cap
bounds actual enthalpy solves across all inlet, pressure, wall and radiation
trials, with at most 48 per radiation response. Numerical tolerances never
widen material support. Hydraulic heat stays a bounded omission (0.01 W and
0.01*dt J), not an automatically deposited source.

## Physical and verification limits

Radiation uses the existing area-mean gray-patch law with view factor one and
constant source-resolved emissivity. It is not pointwise T(x)^4 quadrature,
temperature-dependent emissivity, enclosure exchange, occlusion, participating
media or fluid radiation. Solid geometry/reference mass stay fixed. Reservoir
pressure is fixed with volume free to change, and pipes have no storage or
residence delay. No complete mechanical-energy certificate is claimed.

```sh
cargo test -p fs-airflow --test stored_cooling --test stored_radiation
cargo test -p fs-airflow --example recirculating_cooling
```

Five new library and three example regressions exercise the actual FEM path,
material-only substitutions, signed external heat, original command admission,
heat partition, and failure/cancellation rollback. Native execution was not
available locally because Cargo/rustc are absent. The existing cooling workflow
now requests these tests; its result must be inspected separately. Independent
NumPy spatial-equation comparisons are not execution of the Rust/FEM/IQN code.
