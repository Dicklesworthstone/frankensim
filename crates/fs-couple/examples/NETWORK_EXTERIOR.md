# Coupled openings: exterior feedback and microphone from one physical solve

`driving_point::network::radiation::RadiatingNetwork` connects the original
`fs-duct` graph and the original `fs-bem` closed-body radiation solver.
This is a **harmonic, frequency-domain** construction, not a time-domain
reed/jet oscillator, PCM renderer or measured instrument model.

For every declared opening, a unit outward volume flow becomes a uniform normal
velocity over its physical cap panels. One complete BEM batch supplies every
surface pressure. Area-weighted pressure averages form the full acoustic load
matrix, including all cross-opening entries. The same BEM solutions produce
pressure and Cartesian pressure gradients at the declared exterior receivers.

The original network Schur solve determines the actual opening flows under that
full load. Those flows coherently combine the original exterior fields; internal
port pressure is not used as microphone pressure. All phasors use exp(-i omega t).
Network currents point inward, opening flows outward, hence `p = -Z U_network`.
The same area weights appear in flow normalization and mechanical/acoustic work.
`exterior_surface_power_w` is the raw mean surface-work diagnostic, not a measured
far-field-power integral or time-window energy certificate.

## Construction and use

1. Declare a `NetworkDrivingPoint::new_with_coupled_ports` with one active inlet,
   independent scalar loads where explicitly needed, and the exact ordered labels
   for all coupled openings. These sets must be disjoint and cover every port.
2. Supply a closed outward `SpherePanels` mesh with retained triangles, including
   actual opening caps. `OpeningPatch` assigns disjoint original panel indices to
   each coupled label in that order. Unassigned panels are rigid.
3. Construct `RadiatingNetwork` with the original coherent `GasState`, exterior
   receiver positions, minimum clearance, source `ContentHash`, work limits and
   cancellation gate. Geometry/source/medium are retained immutably.
4. Call `solve(omega, loss_model, injected_volume_flow, gate)`. Its result retains
   the computed boundary matrix, original unit-flow network response, actual
   inlet phasor, actual outward flows and actual receiver pressure/gradients.

The gas supplies duct viscothermal propagation and exterior density/sound speed
without a separate air default. Rebuild with a new admitted gas state to compare
pressure, temperature or composition changes. Receiver relocation does not alter
the boundary matrix. Do not add an independent termination/end correction for the
same exterior cap or a second 1/r, delay or volume-flow source after observation.
Other independent loads do not contribute an invented exterior sound source.

```bash
cargo run --release -p fs-couple --example network_exterior > exterior.csv
```

This example has a synthetic Y graph and two triangular caps on a rigid faceted
sphere. Branch circles have the same area as the caps; this is a scalar port
coupling fixture, **not** proof of three-dimensional junction or real instrument
geometry equivalence. CSV values are harmonic pressure and flow, not audio samples.
There is no default gain adjustment or fitting to a recording.

## Limits and refinement

The original Burton-Miller solver, mesh wavelength check and finite receiver
geometry admission are reused. A coarse or active complete impedance matrix
refuses in the duct owner's full Hermitian passivity test. Entries are not
symmetrized, diagonalized, clipped or repaired to pass it. The original dense
solver has bounded but synchronous factorization; cancellation is checked between
operators, not promised during LU. Field and logical-work caps cover the declared
complete opening basis and receivers; they are not numerical error bounds.

This construction assumes rigid walls, one uniform opening coordinate per cap,
one uniform gas and fixed geometry. The caller owns interior/exterior reference
plane correspondence. It does not infer evanescent junction modes, elastic wall
motion, moving keys, external gas gradients, or a passive causal time realization.
Mesh/port-basis refinement, centroid receiver quadrature error and experimental
validation remain separate. Do not use a single sample as a broadband load.

Native Rust tests were added but not executed in the authoring environment:

```bash
cargo test -p fs-couple --lib driving_point::network -- --test-threads=1
cargo check -p fs-couple --example network_exterior
```
