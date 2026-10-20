# Material-driven gas networks with moving inertial walls

`fs_couple::cavity_phs::piston::exchange::graph::step_candidate` composes the
existing unsplit nozzle/material network with existing `GasPiston` dynamics.
`piston::material::graph::step_candidate` accepts `MaterialGasPiston` instead,
retaining each original `fs_mbd_material::ResolvedBody`, regional density claims,
source conditions, canonical inertia and accepted port history. No numeric mass
override or new material resolver is introduced.

Supply disjoint maps of moving chambers and finite fixed-volume chambers,
physical `SnapshotNozzleLink`s, one shared duration, complete endpoint gas controls,
and exactly one `MovingGasDrive` per moving chamber. Every drive uses the same
physical interval. Return is a `MovingGasCandidate` containing all moving owners,
all fixed-volume owners and the whole-network account; inputs are immutable.

```rust,ignore
use fs_couple::cavity_phs::piston::material::graph::step_candidate;

let proposal = step_candidate(
    &material_walls, &finite_reservoirs, &links, duration_s,
    &endpoint_controls, &wall_drives, network_controls, || false,
)?;
let (next_walls, next_reservoirs, account) = proposal.into_parts();
// Publish both maps together. Internal enthalpy is already in the gas states.
let pressure_pa = next_walls[&chamber_id].gas().state().state().pressure;
let sound_speed_m_s = next_walls[&chamber_id].gas().state().state().sound_speed;
```

## One physical interval, original constitutive owners

1. Every nozzle sees the same accepted pressures, compositions and material laws.
   The existing snapshot owner resolves all mass/species/enthalpy fluxes together,
   retaining choking and original source-backed throat states. A newly arriving
   species cannot traverse another outlet within this explicit flow stage.
2. Each wall advances its reconstructed gas with the original heat/work/momentum
   solver, retaining the same mass inventory and species during that substep.
   Wall displacement changes volume, pressure, temperature and the common
   acoustic/transport medium used by the next interval's nozzle calculation.
3. Final gas storage and wall kinetic energy are independently compared with
   exterior-pressure work, applied-force work and externally supplied heat:

   `sum(delta U_gas) + sum(delta K_wall) + sum(p_ext delta V)
      - sum(F_external delta x) - sum(Q_external) = 0`.

The balance includes every finite chamber, not just each pair in isolation.
Transported enthalpy is internal, including upstream flow work and thermalized
jet energy; it must not be deposited a second time as prescribed heat. The
control-volume convention is the usual open-system first law, summarized in
[Cantera's control-volume equations](https://cantera.org/stable/reference/reactors/controlreactor.html).
No Cantera code or runtime dependency is used.

Original flow withdrawal, source, equilibrium-crossing and energy gates remain.
Original wall temperature/volume brackets, impulse tolerance and compressibility
clock admission remain. A failure in the last wall or final cancellation returns
no endpoints, including any previously solved trial walls. Link reporting order
and input order do not change physical updates. Closed paths retain ordinary
closed-piston results and the original two-port audit.

An actively exchanging chamber's existing wall/heat-only port schema cannot
account for mass transport. Its real wall and heat observations are retained,
but its two-port audit is explicitly unavailable, including balanced throughflow
with zero net mass. The returned flow accounts and whole-system balance cover
that boundary. No stream-port or whole-system passivity certificate is invented.

## Runnable material-only substitutions

From the repository root:

```sh
cargo run --release -p fs-couple --example material_moving_gas_network -- \
  crates/fs-material/examples/material_gas_spring/mixture-a.tsv \
  crates/fs-material/examples/material_gas_spring/mixture-b.tsv \
  1000 40 0.00005 > light-wall.csv

# Change only the first wall's density input (kg/m3), not its geometry or solver.
cargo run --release -p fs-couple --example material_moving_gas_network -- \
  crates/fs-material/examples/material_gas_spring/mixture-a.tsv \
  crates/fs-material/examples/material_gas_spring/mixture-b.tsv \
  8000 40 0.00005 > dense-wall.csv

cargo test -p fs-couple --lib cavity_phs::piston::exchange::graph::
cargo test -p fs-couple --lib cavity_phs::piston::material::graph::
cargo test -p fs-couple --example material_moving_gas_network
```

The example reuses the existing synthetic NASA-9/transport loader and files.
The density input is lowered into an explicit synthetic material card, resolved
at 300 K through the existing regional inertia compiler. Both walls are 0.1 x
0.1 x 0.01 m cuboids with 0.01 m2 pressure areas; the second wall retains density
8000 kg/m3. Initial gas volumes are 0.001 m3 at 400 K / 100000 Pa before explicit
150 J supply and 30 J chamber preparation. Initial wall velocity is 0.02 m/s.
A common finite supply feeds both moving chambers; both discharge into a finite
sink. Each chamber receives 100 W for the first half of the requested intervals,
then zero. CSV retains motion, gas inventories, thermal/transport/acoustic state,
source-derived wall mass and step/history energy residuals. Swapping the gas files
changes only their declared mixture data. These are synthetic sources, not
measured properties of named gases or solids.

## Limits and verification

This is a first-order transfer-then-motion split. Both stages represent the SAME
elapsed interval, not two elapsed intervals. It is not a simultaneous implicit
network solve: a wall-induced pressure reversal drives the next interval, and
external events should split the physical timestep. Refine time independently
of algebraic tolerances. Uniform ideal-gas source/support limits still apply.

Walls are independent constant-area rigid translations, not two faces of one
shared piston, crankshaft or deforming/melting body. Solid source state is frozen;
gas heating does not silently heat or re-resolve wall density. This composition
has no reactions, resolved jet momentum, spatial pipe acoustics, finite thermal
wall stores or electrical network coupling. Existing dedicated thermal/reactive
owners remain separate. Sound speed is an exposed physical property, not a live
acoustic render. Energy acceptance is not an entropy-production or stability proof.

Twelve new Rust regressions cover shared flows and actual motion, original-owner
parity, symmetry, source shutoff, material density changes, source retention,
first-order refinement, cancellation and failure atomicity. Native compilation
and tests were not run successfully during authoring: Cargo/rustc are absent.
Independent NumPy/SciPy equations ran 522 coupled steps; velocity errors versus
DOP853 decreased 0.00046696 -> 0.00023373 -> 0.00011693 m/s (difference ratio
1.99674), with maximum history energy discrepancy 4.06e-10 J. Those references
are not Rust execution, empirical material validation or DSR certification.
