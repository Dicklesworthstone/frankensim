# Fit spatial thermal measurements, not invented vertex targets

`ScanPowerObjective<O = ScanObservation>`, `ScanPowerStudy<F, O = ScanObservation>`
and `ConstrainedScanStudy<F, O = ScanObservation>` now accept the original sealed
`ScanObservationSpec`. Existing vertex-objective callers keep their default
parameter. For spatial data, supply `Vec<ScanSpatialObservation<'m>>` instead.
Type inference selects the spatial observation layout from that vector.

`ScanSpatialObservation`, `ScanReadout`, and `ScanObservable` live under
`fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective`
(the first two are in its `spatial` module). A spatial row can represent an
explicit vertex, a P1 point in a specified tetrahedron, a volume mean over
specified tetrahedra, or a normalized nonnegative sensor-weight row. A vertex
row in a spatial observation vector also permits mixed point/footprint data.

Build each readout against the **same mesh object** the factory returns. Merely
matching vertex counts is insufficient. The immutable readout borrows that
mesh; every candidate retains its geometry and weights. Powers, source paths
and source widths can still be optimized. Readout locations are not variables.

The actual nodal temperature, specific enthalpy, or liquid mass fraction is
sampled first. The row then forms **one prediction and one fitting residual**.
Its transpose distributes that residual back to the nodes before the local
material-chart derivative. Assigning the same measured mean to separate vertex
targets is a different optimization problem and is not performed here. Neither
temperature nor fraction is inferred from averaged enthalpy.

Both L-BFGS and SQP use this same engine. Optional checkpointing shares replay
sections across the objective and all attached point/peak constraints. Existing
power/energy/geometry limits, exact off-state coordinates, source normalization,
material/phase/energy acceptance, cancellation, and evaluation accounting stay
in force. State and peak constraints keep their existing **nodal** meaning;
choosing an averaged fitting sensor does not replace a hot-spot limit by a mean.

A volume average of liquid **mass** fraction is not generally a bulk mass
fraction. Use declared reference-mass weights for the latter. Temperature data
on a perfectly isothermal latent plateau cannot identify powers on their own;
additional informative observations are required. No identifiability, global
optimum, sensor dynamics, or experimental-calibration guarantee is inferred.

## CSV-driven executable consumer

From the repository root:

```bash
cargo run --manifest-path examples/scan_power_fit/Cargo.toml \
  --example spatial_measurement_fit -- \
  examples/scan_power_fit/examples/spatial_measurements.csv 0 0
```

This executable deliberately fixes a **synthetic** unit-cube, unit-density
specimen with eight vertices and six tetrahedra, isotropic conductivity 2 W/(m K),
adiabatic boundaries, initial specific enthalpy 150 J/kg and temperature 310 K.
Its declared latent plateau runs from 100 to 200 J/kg. A normalized absorbed-power
source with half widths [1.4, 1.2, 1.1] m sits at [0.25, 0.5, 0.5] m for the first
second and [0.75, 0.5, 0.5] m for the second. It fits two nonnegative power values;
the shutdown at 2 seconds is fixed at exactly zero. The library API has no such
fixture restriction: applications supply their own physical factory and mesh.

Input is a bounded, unquoted CSV dialect: at most 128 KiB and 256 nonempty
measurement rows. Blank lines and `#` comment lines after the header are ignored.
The header is exactly:

```text
time_s,observable,target,scale,readout,selection
```

`observable` explicitly selects `temperature_k`, `specific_enthalpy_j_kg`, or
`liquid_mass_fraction`; target and positive scale use those units. Fractions use
[0,1], not percent. `time_s` lies in (0,2]. The fixed grid is the sorted union of
measurement times and the 1/2-second source switches, without tolerance merging.
Measurement order and duplicate observations are retained. Event-left source
sampling preserves the heat in an interval ending at a switch.

Selections (all IDs zero-based) are:

- `vertex,3`: one original vertex.
- `tet,0:1:2:3:4`: tetrahedron 0, homogeneous barycentric weights [1,2,3,4].
- `volume,0;1;2;3;4;5`: volume mean over the listed increasing tetrahedron IDs.
- `weighted,0:1;3:2;7:1`: normalized sensor weights on distinct increasing vertices.

No readout is snapped to a nearest node. Invalid selections, foreign mesh
readouts, nonfinite input, unsupported quantities, and exhausted input allowances
refuse. Underlying source, material and nonlinear failures propagate; invalid
physics is never replaced by a fitting penalty.

The supplied CSV asks for mean enthalpies 152/157 J/kg and final mean liquid
fraction 0.57. On this reference specimen its analytical schedule is 2/5 W.
These means do not require a spatially uniform final enthalpy field. The fitting
example reads external target values; it does not generate targets inside the
optimizer. Replaying the exported schedule consumes a separately reserved final
evaluation. Output includes the actual unsnapped power schedule plus each
prediction and standardized residual. A converged noisy/inconsistent-data fit
can retain nonzero mismatch, which is reported rather than called exact recovery.
Failed convergence or feasibility checks produce an unsuccessful process exit.

## Focused native checks

```bash
cargo test --manifest-path examples/scan_power_fit/Cargo.toml --test spatial_design
cargo test --manifest-path examples/scan_power_fit/Cargo.toml --example spatial_measurement_fit
```

Authored regressions cover direct generic-evaluator equality, spatial L-BFGS
recovery, complete source/geometry coordinate gradients, full/replayed equality,
SQP point/peak limits and zero-power derivatives, vertex reduction, foreign
meshes/cancellation, and CSV admission/order/fit behavior. Native compilation and
execution were unavailable in the implementation session (Cargo/rustc and the
DSR/RCH runners were absent). Independent Python reference results, when cited,
do not establish a native Rust pass.
