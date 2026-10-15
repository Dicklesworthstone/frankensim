# One material state for current, heat diffusion and fluid momentum

`MaterialSlab::step_magnetothermal_with_material_transport` and its
`step_magnetothermal_isolated_with_material_transport` counterpart connect the
existing moving-current and local thermal-flow owners in 2D and 3D. Zero imposed
B gives Ohmic heating with local transport. Existing frozen-transport and midpoint
entry points remain unchanged; they do not silently enable this new behavior.

The callback receives `(cell, kelvin, retained_transport)` and returns
`material::electrothermal::magnetic::transport::ElectricalTransport<F>`, where F
is any existing `FluidTransport` implementation. Pair the original concrete
material evaluation with `Some(resistivity_ohm_m)` or explicit insulation (`None`).
Missing or unsupported data must be an error, never guessed insulation. There
is no material-name dispatch, automatic mixing law or new material database.

Every fluid cell is sampled once at each evaluation epoch. Unlike the old
resistivity-only callback, this includes ideal electrode cells because they still
transport heat and momentum. Their resistivity is unused; None is sufficient.
Fixed-temperature wall reservoirs are not queried. The isolated path retains all
material contributions at numerical gauge cells. Concrete F objects, their exact
transport snapshots and applied versus endpoint epochs remain in the result.

The electrical pass reads resistance from those same applied samples, without
re-querying the material object. The original current owner supplies bulk/contact
Joule heat and reciprocal Lorentz force. The shared physical source conversion
and local-transport kernel use the corresponding k/mu. Existing contact partitions,
external forces, buoyancy, independently measured wall exchanges, low-Mach checks,
current/power checks and atomic failure semantics remain active. No second current,
flow or heat discretization is introduced. Endpoint material support and lattice
admission are checked before the coupled step replaces the original slab.

## Model boundaries

Sources and coefficients are evaluated at the accepted INITIAL state and held
for one step. This is not the within-step midpoint algorithm. Refine time for
strong feedback; endpoint support alone does not prove support throughout the
intervening trajectory. Density and Cp must match the retained uniform storage
snapshot bit-for-bit. Changing them requires conservative mass/enthalpy evolution,
not another coefficient assignment. Expansion/buoyancy and geometry stay fixed.

B is imposed, stationary within the step. This does not solve magnetic induction,
Hall effects, electrochemistry, phase changes or field-source dynamics. In 2D only
out-of-plane Bz is supported, with explicit extrusion depth. Electrical exterior
faces insulate even where fluid boundaries are periodic. Wired conductors require
voltage references; isolated stepping requires an empty electrode list.

The returned mechanical energy is initial Lorentz power times dt, not the measured
change of kinetic energy. Electrical-to-heat/work accounting is not exact total
energy conservation of the split time integrator. External mechanical forcing and
viscous/numerical dissipation are not silently converted into electrical heat.
Cancellation is polled between material queries and around complete thermal steps;
callbacks must return. Original, tentative and thermal-owner scratch coexist.
The cell limit bounds grid size, not total bytes. No speed or stability claim.

## Existing executable, extended table

`material_transport_feedback` keeps its original prescribed-heating command and
three-column temperature table. Electrical mode requires a fourth, positive
resistivity column at EVERY temperature knot. All three local properties use the
same interpolation bracket. Mixed row formats, unsupported temperatures and
missing electrical data refuse. The existing 64 KiB/256-knot bounds remain.
These example tables are synthetic, not production source packs or calibrated
properties of a named substance.

```sh
cargo run -p fs-lbm --example material_transport_feedback -- \
  crates/fs-lbm/examples/material_transport_feedback/electrical_constant.txt \
  100 --electrical 10 0.5
cargo run -p fs-lbm --example material_transport_feedback -- \
  crates/fs-lbm/examples/material_transport_feedback/electrical_temperature_dependent.txt \
  100 --electrical 10 0.5
cargo test -p fs-lbm --test material_electromagnetic_transport
cargo test -p fs-lbm --example material_transport_feedback
```

The two tables have identical density, Cp and resistivity(T); ONLY k(T) and mu(T)
differ. They therefore begin with the same electrical current. Subsequent local
transport changes temperature and velocity, which feed back into later current
calculations. CSV reports the actual current, Joule power, magnetic work, fluid
velocity, temperature and applied k/mu ranges, plus coupling balance.

The voltage becomes zero halfway through, with B and the external force retained.
That SHORTS the ideal electrodes; it does not open the circuit. Moving fluid may
continue generating current, braking and Joule heat after voltage switch-off.
No state is reset. The original prescribed-heating mode remains separate and
retains its original source schedule and CSV layout.

Six core integration tests plus three new example tests are authored. Native
compilation/tests/example execution were not available during implementation
(Cargo/rustc absent). Independent circuit/lattice equations passed 24 cases with
local heat and current checks. Such reference equations are not execution of
repository Rust or experimental validation.
