# Two material-defined textile walls and one shared gas cavity

`textile_cavity` joins the existing source-bound membrane equilibrium and
prestressed modal basis to the existing gas, Helmholtz and reciprocal modal
network owners. Both sheets and the gas remain independent data inputs. The
second sheet is not a prescribed oscillator or read-only output: its displaced
volume changes the pressure loading the first sheet and the common neck.

```sh
cargo run -p fs-couple --example textile_cavity -- \
  /path/to/sheet-a.pack /path/to/sheet-b.pack \
  crates/fs-couple/examples/textile_cavity/air.txt
```

The same material pack may occupy both positions. Compile FrankenTextiles
records with the existing importer and `matdb-pack` path described in
`FRANKENTEXTILES_INTEGRATION.md`. Each pack needs characterized reference areal
mass, thickness and all four effective membrane properties; bending, thermal
conductivity and heat capacity are NOT inferred or required. Supplied integration
fixtures are synthetic, not empirical textile calibration.

Change only `air.txt` to `co2.txt` to replace part of the nitrogen with carbon
dioxide while keeping oxygen fraction, total pressure, temperature, geometry,
sheet data and drive fixed. The original `fs-material` gas-record parser owns
both legacy `T P RH` and explicit mass/mole-fraction mixture inputs; this example
adds no species database, normalization or constitutive mixing rule. Supported
named species and domain/condensation refusals are exactly that owner's.

## Input, with explicit SI units

The bounded file contains these records, in this order, with no trailing fields:

```text
textile-cavity-v1
gas mixture 293.15 101325 mole 2 N2 0.79 O2 0.21
cavity 0.001 0.005 0.02 2000
surface 293.15 1.02 0.3 0.1
drive 48000 2400 240 -0.001
limits 0.0001 0.01
```

`cavity` specifies reference volume [m3], circular neck radius and length [m],
and the explicit acoustic resistance [Pa s/m3]. Original end correction,
inertia and compliance are evaluated once by the existing Helmholtz model.

`surface` specifies material query temperature [K], imposed perimeter stretch
ratio, material grain [rad], and the caller's maximum Green-strain component.
The fixed reference scene has two opposing 0.1 m square, eight-triangle walls,
with three free Cartesian centre coordinates each. Reference volume is the
volume at the prestretched equilibrium, not before stretching. Both sides of
each wall have the same equilibrium static pressure; there is no imposed static
pressure difference requiring an additional preload solve.

`drive` specifies mechanical samples per second, total frames, pulse frames,
and signed centre-node force on the first wall [N]. The pulse is held for exactly
the stated number of steps. The second wall receives no external force. Its
motion and the neck flow arise from reciprocal pressure feedback.

`limits` specifies maximum Cartesian nodal displacement increment [m] and maximum
fractional cavity-volume change. Coordinate ceilings derived from modal row
norms enforce displacement in the original component transaction. A conservative
total-energy ceiling enforces volume change via compression energy; it can
refuse high mechanical energy even when instantaneous volume changes cancel.
The input caps these limits at 1 mm and 1%, respectively. These numerical
operating bounds are not empirical accuracy certificates.

Only a completed run emits CSV: both wall increments, actual cavity-pressure
increment, outward neck flow, fractional volume change, stored incremental
energy, external work, disclosed loss and energy residual. Stored energy excludes
the already equilibrated membrane preload; it includes structural perturbation,
neck kinetic and gas compression energy. Compression is counted once. Native
Rust execution was unavailable in the authoring environment; no completed CSV,
listening test or empirical material validation is claimed.

## Library reuse

`fs_couple::cavity_modal::shared::{CavitySurface, SharedCavity}` accepts 1..64
original modal surfaces (at most 4096 total coordinates including the neck),
with caller-owned clock and state/energy/work budgets. All original source
receipts and geometric eigenvectors remain with the upstream material/structure
owner. The executable retains both original packs, cards, patches and vibration
objects throughout simulation. No new runtime dependency is introduced.

The original Helmholtz lowering supplies each signed shape; the original
whole-column modal solver uses one collective compression connection. Thus
opposite volumes cancel and every inter-surface cross term remains. This is NOT
one spring per surface. `step` accepts surface generalized forces and an optional
nonzero exterior neck-pressure value, whose work belongs to the external source
account. `observe` reads accepted endpoint pressure/flow/storage; it is not a
second time step and its diagnostic errors do not roll back a prior step.

The model remains fixed-basis, small-signal, uniform-pressure, single-neck gas
loading. It does not model nonlinear textile hysteresis, wrinkles, moving mesh,
heat/moisture exchange, internal standing waves or broadband radiation. Structural
damping and direct pressure transfers are explicitly zero in this example; neck
loss is explicit. The CSV pressure is INTERNAL cavity pressure, not microphone
sound. Original acoustic receivers may consume the neck flow separately.

Focused native checks:

```sh
cargo test -p fs-couple --lib cavity_modal::
cargo test -p fs-solid --test membrane
cargo test -p fs-couple --example textile_cavity
```
