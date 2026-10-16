# Coupled material transport at the magnetic midpoint

`MaterialSlab::step_magnetothermal_midpoint_with_material_transport` and
`step_magnetothermal_midpoint_isolated_with_material_transport` connect local
viscosity, thermal conductivity and electrical resistivity to the SAME within-step
temperature/velocity iteration in both 2D and 3D. They reuse the existing magnetic
midpoint driver, material sampler, current/contact/force owner, physical source
conversion and thermal/flow kernels. The original frozen-transport midpoint APIs
use that same driver without extra material queries; explicit APIs stay explicit.

For each candidate midpoint, the callback returns the existing
`ElectricalTransport<F>` at each fluid cell's Kelvin temperature, where F is any
`FluidTransport` implementation. Its transport snapshot is extracted once. The
same prepared object supplies resistivity to the current solve and k/mu to the
actual local collision rates. Concrete source-backed objects retain their own
receipts; synthetic coefficients do not acquire measurement authority. Electrode
cells still transport heat and momentum; their electrical resistivity is unused.
Fixed-temperature walls are separate reservoirs and are not sampled. Numerical
potential gauges in isolated conductors keep all their physical contributions.

Every trial restarts from the SAME accepted populations. Neither iterating nor
re-evaluating material data advances physical time or repeatedly deposits heat.
The accepted result retains the source-producing current/force/heat fields, the
applied material objects and their evaluation coordinates. Endpoint material
objects are checked before publication and retained separately; they never replace
the properties that actually drove the step. Missing data, invalid storage changes,
failed electrical/thermal balances, nonconvergence or observed cancellation leave
the accepted slab intact. Captured callback side effects remain caller-owned.

## Independent acceptance criteria

`material::electrothermal::magnetic::transport::midpoint::MaterialMidpointOptions`
contains the original `MagneticMidpointOptions` and a local transport tolerance.
Temperature and velocity coordinate mismatches must meet their K and m/s limits.
A fresh electrical solve at the actual old/new midpoint checks the sums of
absolute local Joule-energy, magnetic-impulse and magnetic-work differences.
Opposing local errors cannot cancel in these checks.

In addition, every fluid cell's positive thermal conductivity and viscosity must
satisfy `abs(applied - checked) / max(applied, checked)` within
`transport_relative_tolerance` (finite, in [0,1]). Zero requests exact agreement;
one admits any finite positive coefficient difference and is intentionally weak.
The returned `TransportResidual` reports separate maxima for k and mu. This gate
matters even with ZERO current: electrical residuals alone say nothing about
whether a nonlinear thermal or momentum coefficient has converged.

The material callback must evaluate the same stated law on repeated queries.
Midpoint iteration can refuse discontinuous laws or exhaust its budget; that is
not proof of physical runaway. The iteration budget remains 1..128 with at most
two current solves per trial. Each trial and current recheck samples all fluid
cells; acceptance adds one endpoint material evaluation, not another current solve.

## Existing executable

The existing four-column electrical tables, prescribed-heating command and
explicit electrical command remain unchanged. Add a final `--midpoint` flag:

```sh
cargo run -p fs-lbm --example material_transport_feedback -- \
  crates/fs-lbm/examples/material_transport_feedback/electrical_constant.txt \
  100 --electrical 10 0.5 --midpoint
cargo run -p fs-lbm --example material_transport_feedback -- \
  crates/fs-lbm/examples/material_transport_feedback/electrical_temperature_dependent.txt \
  100 --electrical 10 0.5 --midpoint
cargo test -p fs-lbm --test material_transport_midpoint
cargo test -p fs-lbm --example material_transport_feedback
```

Only k(T) and mu(T) differ between these synthetic tables. Geometry, density, Cp,
resistivity(T), drive and imposed field are identical. Midpoint CSV adds actual
nonlinear trial counts and the largest coefficient-consistency residual. It labels
mechanical energy as midpoint work; the original explicit CSV and lagged-work
label remain unchanged. The voltage becomes zero halfway through, SHORTING the
ideal electrodes, not opening the circuit. Motion, heat and short-circuit braking
history are retained. Midpoint feedback can change even the first accepted current;
there is no requirement that it equal the current evaluated at the initial state.

## Physical boundaries and verification

Density and Cp remain equal to the retained uniform storage snapshot. Changing
them requires conservative mass/enthalpy evolution, not in-place coefficient
assignment. Geometry, expansion/buoyancy parameters, contacts, electrical drives,
imposed B and external forces are fixed within a call. Only out-of-plane Bz is
admitted in 2D, with explicit extrusion depth. Electrical exterior faces insulate
even where the fluid has periodic lateral boundaries. The isolated method requires
no electrodes and an explicit local current-continuity tolerance in amperes.

The work account is source Lorentz power at the evaluation velocity times dt,
NOT the measured kinetic-energy change. This is not magnetic induction, Hall
physics, electrochemistry, latent heat, exact kinetic-plus-thermal conservation,
a second-order claim for the split lattice integrator, unconditional stability,
or experimental calibration. Midpoint/endpoint support checks do not establish
continuous source support between those points. Existing low-Mach and transport
accuracy limits apply. Cancellation is polled around whole kernels and between
material queries; callbacks must return. Cell limits are not total memory limits.

Seven core regressions and three additional example regressions are authored.
Native Cargo commands could not start in the implementation environment because
Cargo/rustc are absent: Rust compilation, tests and executable runs remain
unverified. Independent NumPy/SciPy equations passed 24 warmed-state comparisons
against simultaneous nonlinear roots and four complete table-driven trajectories
(including voltage reversal and zero-field checks). Those references are not
execution of repository Rust or validation against measured physical data.
