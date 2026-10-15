# Material reeds between two acoustic loads

`two_sided_reed` derives two independent blade modes with the existing DKT plate
and modal owners, connects them to upstream and downstream acoustic networks,
and observes the downstream terminal volume flow through the shared retarded
Rayleigh receiver. The WAV uses a fixed 0.005 Pa full-scale, never peak normalization.

From the repository root, supply NEW output paths:

```sh
cargo run --release -p fs-couple --example two_sided_reed -- short-tract.wav 0.12 900
cargo run --release -p fs-couple --example two_sided_reed -- long-tract.wav 0.24 900
cargo run --release -p fs-couple --example two_sided_reed -- heavy-blade.wav 0.12 1800
```

The first two runs change only upstream length in metres. The third changes
only the second blade's density in kg/m3. Mass and resonance follow the plate
reduction; there is no frequency input or compensating pressure adjustment.
The canonical `two_sided_reed.gesture` data supplies 5 Pa for 960 samples, then
releases for the remainder of a 4800-sample, 48 kHz window. Both blades, both
ducts and all existing load memories keep evolving during that tail. This is a
coupling experiment, not a calibrated sustained oboe/clarinet preset. The plate
sections, 0.7 upstream reflectance and effective shared contact are explicitly
authored numerical data. The downstream terminal uses the original geometry/gas-
derived compact baffled radiation load. Distributed viscothermal loss is not
selected in this example.

## Library composition

`bernoulli_aperture::two_sided::combined_inlet_state` derives the summed inlet
impedance and checks matching coherent inlet gas states. Construct either
`DynamicAperture` or `PairedAperture` with that load, then pass it into
`TwoSidedAperture::new`. Original material-bound plates, spatial closure and
independent material histories stay in those existing valve objects. Each side
is an `ApertureDuctSpec` with its original graph, complete regional gas mapping,
explicit viscothermal selections and existing memory budget. Far sections may
have different frozen gas states; the two valve inlets must match.

Each inlet points AWAY from the valve. For through-flow Q:

```text
b_up = a_up - Z_up Q
b_down = a_down + Z_down Q
P_up = P_source + a_up + b_up
P_down = a_down + b_down
DeltaP = P_source + 2(a_up-a_down) - (Z_up+Z_down)Q
```

The original nonlinear valve solves this combined load. Returned pressure
observations are mapped back to the actual two faces. The complete work account
includes the explicit series pressure source, independent mechanical forces,
valve/contact/material losses and BOTH acoustic networks' waves and reactive
loads. Failed trials publish none of the three participants. `advance_block`
retains exact accepted prefixes for cancellation and budget exhaustion.

`two_sided::observer::TwoSidedPressure` binds one selected passive outlet to a
`CircularOutletReceiver` using its actual radius and gas. It implements the
existing `PressureRenderer`, so the common PCM/decimation/ensemble adapters can
consume it. Physical propagation latency stays in the signal. The receiver is
one-way prescribed-flow radiation; it does not invent or double an acoustic load.
The ordinary constructors retain their finite-difference/delay approximation
and declared frequency band. A failed callback refuses continuation.

`two_sided::observer::schedule::ScheduledTwoSidedPressure` attaches an existing
`GestureSchedule` and `AperturePerformanceConfig` to that unadvanced pressure
source. It reuses the SAME pressure compiler and candidate-control evaluator as
scalar and paired aperture performance, not a second sampler or envelope law.
Controls apply at their exact mechanical sample even inside a callback or an
output-decimation group. Independent mechanical blade forces are preserved and
may be changed between callbacks; pressure release resets no physical state.
Nonzero pre-held pressure, late attachment, incompatible clocks and incomplete
control budgets refuse rather than silently replacing inputs. Inspect the
original schedule, pending/applied controls and retained source through read-only
accessors. The example uses this path instead of host-block pressure conditionals.

The series-load approximation is described by Fritz and Wolfe, JASA 118 (2005),
*How do clarinet players adjust the resonances of their vocal tracts for different
playing effects*, https://arxiv.org/abs/physics/0505195 . It is not an anatomical
lung model, compressible mean-flow simulation, live gas-mixing law, time-varying
tract geometry or empirical instrument-fidelity certificate. Component model,
resolution, contact and material applicability limits remain unchanged.

## Focused native checks

```sh
cargo test -p fs-vfit --test waveguide_pair
cargo test -p fs-couple --test two_sided_aperture
cargo test -p fs-couple --test two_sided_pressure
cargo check -p fs-couple --example two_sided_reed
```

These 14 authored tests cover network publication, original single/paired
mechanics, two-sided work/flow, return causality, material mass changes,
refusals, resume, exterior pressure and canonical phrase timing across callback
sizes. Native compilation/tests and the example WAV have NOT been executed in
the implementation environment: Cargo/rustc are absent. The existing GitHub job
has begun but dependency bootstrap failed; no native pass is inferred from a
pending downstream test step. Independent Python/SciPy references passed 360
simultaneous-equation cases and 2304 complete two-network steps; these are
equation checks, not Rust execution, DKT validation or measured timbre.
