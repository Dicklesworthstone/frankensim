# Coupled bulk materials and electrical interfaces

`MaterialSlab::step_magnetothermal_midpoint_with_interfaces` connects ordered
interface response to the existing local-material/current/heat/flow iteration.
Both 2D and 3D use the same prepared-state midpoint driver. Every nonlinear trial
starts from the same accepted populations. No separate current, fluid, heat,
contact-resistance or material-database implementation is introduced.

The public types are in
`material::electrothermal::magnetic::transport::midpoint::interfaces`.
`InterfaceMidpointOptions` carries the existing material/coordinate/source controls,
independent resistance and heat-partition checks, and explicit wired/isolated
circuit selection. The isolated mode requires an empty electrode list and retains
the existing per-component gauges and current-continuity check.

## One evaluation epoch, ordered receiving cells

`ElectricHeatInput.contacts` supplies the fixed face sites, their order, and
callback templates. Its resistance and partition values are NOT applied as a
fallback: only `InterfaceResponse` values enter the current calculation.
The callback receives `(index, template, [first_bulk_sample, second_bulk_sample])`.
Those are the exact sampled objects used for electrical resistivity, viscosity
and conductivity at this epoch, including their Kelvin temperatures and concrete
source payloads. It does not independently re-query a bulk material.

The temperatures are the adjacent **cell-centered** receiving states. A physical
surface trace, film temperature, pressure/roughness/history law or other interface
query convention must be specified by the caller, not inferred by the numerical
adapter. The example explicitly uses their arithmetic mean as its film-temperature
approximation and requires both individual temperatures to be within table bounds.

`InterfaceResponse<I>` retains the original concrete interface evaluation as I,
its specific electrical resistance in ohm m2, and the explicit fraction of
interface-generated heat deposited in the first cell. Missing source data must
return an error. It is not interpreted as zero resistance or insulation.
The caller's original error type is retained with contact index and temperatures.
Neither bulk nor interface payloads need Clone.

An existing `ResolvedElectricalInterface` can be retained directly as I: resolve
it with `fs_material::electrical_interface::resolve_electrical_interface_state`,
verify its ordered assignments using `require_surface_materials`, and place
`resolved.resistance_ohm_m2()` and the original resolved object in the response.
Keep source selection, pressure, composition, finish and history coordinates
explicit in that original query. Heat partition remains an independent input;
electrical resistance does not determine thermal resistance or heat partition.
No new source pack or source-backed measurement is claimed by the synthetic demo.

## Acceptance and energy ownership

All existing bulk transport, temperature/velocity, local Joule-energy, magnetic
impulse/work and current/power criteria still apply. Resistance and heat-partition
convergence are checked separately at the actual old/new midpoint, even when
current is zero. Contact discrepancies cannot cancel across different faces.
A zero versus positive contact resistance has relative discrepancy one; two
zeros match. Tolerances govern numerical consistency, not material uncertainty.

The existing current kernel owns I-squared-R generation and bulk/contact heat
partition exactly once. Its cell powers and reciprocal magnetic forces enter
the existing thermal/flow kernels. Do not add a second contact-heating source.
Applied interfaces are returned separately from endpoint support evaluations.
Endpoint interface errors, loss of conduction on a contact side, failed balances,
nonconvergence and observed cancellation leave the accepted slab unchanged.
Callback side effects/history are caller-owned and must not advance during trials.

This is a fixed electrical-face model, not contact mechanics: it does not add a
hydrodynamic wall, adhesion, wear, thermal contact resistance or interface storage.
The underlying fixed-density/Cp Boussinesq and imposed-field approximations remain.
There is no melting, induction, Hall, deforming contact geometry, unconditional
stability, second-order split accuracy or exact kinetic-plus-thermal energy claim.
Reported mechanical work remains evaluated Lorentz power times dt, not a measured
kinetic-energy change. Source objects alone do not confer experimental validation.

## Existing executable

Append an interface file to the existing electrical midpoint command:

```sh
cargo run -p fs-lbm --example material_transport_feedback -- \
  crates/fs-lbm/examples/material_transport_feedback/electrical_constant.txt \
  100 --electrical 10 0.5 --midpoint --interfaces \
  crates/fs-lbm/examples/material_transport_feedback/interface_constant.txt

cargo run -p fs-lbm --example material_transport_feedback -- \
  crates/fs-lbm/examples/material_transport_feedback/electrical_constant.txt \
  100 --electrical 10 0.5 --midpoint --interfaces \
  crates/fs-lbm/examples/material_transport_feedback/interface_temperature_dependent.txt

cargo test -p fs-lbm --test material_interface_midpoint
cargo test -p fs-lbm --example material_transport_feedback
```

Interface rows contain `film_kelvin resistance_ohm_m2 heat_fraction_to_first`.
The bounded example-only reader accepts 2..256 knots and at most 64 KiB, with
piecewise-linear interpolation and no extrapolation. The two synthetic files
change ONLY contact resistance. Bulk data, two interior interface sites, electrode
geometry, external force and voltage schedule remain identical. The voltage is
set to zero halfway through, shorting the ideal electrodes without resetting
motion/temperature; induced current, braking and heating may continue.

The new CSV reports actual current, Joule heat, Lorentz work, velocity, temperature,
applied interface resistance and convergence. All prior invocation forms and
explicit/midpoint CSV layouts remain unchanged when the new option is absent.
Eight core and four example tests are authored. Native compilation/tests were
not available during implementation (Cargo/rustc absent). Independent NumPy/SciPy
circuit/lattice references cover 24 coupled cases and four data-file trajectories;
those checks are not execution of repository Rust or experimental validation.
