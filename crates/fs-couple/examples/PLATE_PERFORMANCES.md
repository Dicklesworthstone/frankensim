# Geometry-driven plate audio

```sh
cargo run -p fs-couple --bin music_render -- \
  plate crates/fs-couple/examples/plate-mesh.performance \
  /tmp/plate-mesh.wav --block 37
```

This path takes an actual flat triangular mesh, a section for every triangle,
explicit supports, a nodal force footprint and a physical force history. It
assembles the existing `fs-plate` DKT pencil, harvests modes with `fs-modal`,
and hosts the resulting `CompactBody` objects in the existing sample scheduler.
The caller does not assign pitches, modal masses, mode shapes or radiation areas.
The checked-in example has two different material/thickness regions and timed
release/re-excitation of one distributed actuator. All material numbers are
authored demonstration inputs, not measurements or a named material data card.

`plate-mesh.performance` is editable UTF-8. Records have the following fixed
order. Each repeated row occurs exactly as often as its preceding count says;
indices are zero-based. ASCII whitespace separates fields. Extra fields, missing
or trailing records, and unknown versions refuse. There are no ignored comments.

```text
frankensim-plate-performance-v2
audio RATE_HZ SAMPLE_COUNT FULL_SCALE_PA
ambient GAS_RECORD
observer DISTANCE_M
limits MAX_MODES NYQUIST_FRACTION MAX_ABS_FORCE_N MAX_ABS_PRESSURE_PA MAX_EVENTS
mechanics SUPPORT PRETENSION_N_M DAMPING_RATIO LOW_HZ HIGH_HZ RETAINED_MODES
sections COUNT
section THICKNESS_M DENSITY_KG_M3 E1_PA E2_PA NU12 G12_PA ANGLE_RAD fixture LABEL
nodes COUNT
node X_M Y_M
triangles COUNT
triangle NODE_A NODE_B NODE_C SECTION_INDEX
supports COUNT
support NODE_INDEX
footprint COUNT
weight NODE_INDEX UNIT_TOTAL_FORCE_WEIGHT
initial_force_n FORCE_N
events COUNT
force SAMPLE_INDEX FORCE_N
```

`GAS_RECORD` is the shared gas declaration read by
`fs_material::gas_species::input::read_gas_record`, as in the bowed, lip and
voice files: `T_K P_PA RH` (moist air), `mixture T_K P_PA mole|mass COUNT
SPECIES FRACTION...`, or `thermochemical T_K P_PA mole|mass COUNT ID
FRACTION...` with `--gas-component` source packs on the command line. Its
density scales the pressure and its sound speed sets the time of flight. A
record that does not parse or resolve refuses the file; there is no air
default. Version 1 files, which authored a bare fluid density, refuse.

`SUPPORT` is `simply-supported` (zero transverse displacement) or `clamped`
(zero displacement and slopes), applied to exactly the listed support nodes.
Unlisted edges remain unconstrained by this declaration. Triangles must be
counterclockwise, conforming and nonoverlapping; this parser is not a general
triangle-soup topology certifier. Shared nodes impose perfect bonding. Each
section is homogeneous through its centered thickness. Its material axis 1 is
rotated counterclockwise from mesh x by `ANGLE_RAD`. No out-of-plane constants,
laminate offsets, interfacial slip, geometric nonlinearity or thermal evolution
are inferred. Use equal E1/E2 and G12 = E/(2(1+NU12)) for isotropic elasticity.

Every section names where its material comes from (bead pexfp.32). Bare numbers
must end in an explicit `fixture LABEL`, or the file refuses; the example's two
regions are labelled `fixture authored-illustrative-plate-not-measured`. A
pack-backed section instead reads `section material isotropic THICKNESS_M PACK
STATE [AXIS=VALUE ...]` or `section material orthotropic12|orthotropic21
THICKNESS_M ANGLE_RAD PACK STATE [AXIS=VALUE ...]`. PACK is a compiled
material-card pack relative to this file, STATE the state its card was compiled at
(its `MaterialStateId` process), and the axes are the query point in the card's
normalized units. The section's density and elasticity resolve through
`fs_material::state_point` (`young_modulus`/`poisson_ratio`, or `young_modulus_1`,
`young_modulus_2`, `shear_modulus_12`, `poisson_ratio_12`). Another state or a
missing property refuses by name. The sidecar's `section_sources` lists each
section's pack, hash, state and point, or its fixture label.

`LOW_HZ`/`HIGH_HZ` define a search window, not imposed frequencies. The existing
solver returns the positive eigenmodes in that window; at most `RETAINED_MODES`
are kept. The sidecar reports requested and actual retained counts separately.
The window need not contain the fundamental or all significant modes. Mesh and
modal convergence remain the caller's responsibility. The existing reduction
interprets the authored damping ratio as its Rayleigh targets at omega0 and
4 omega0; it does not discover material loss from the elastic constants.

Force weights are signed and must sum to one; they are never silently
normalized. A weight at a supported displacement DOF is reacted by the support,
not reassigned to a free node. `FORCE_N` is the signed total force on this fixed
footprint. It starts on an initially resting plate; no settled preload is
implied. Events apply before their named sample, retain source order at equal
times, and hold until the next event. Zero force releases the actuator while
vibration continues. Events at or beyond SAMPLE_COUNT refuse.

The command requires 48 kHz and accepts only `--block` as an override. A source
is capped at 4 MiB, 512 nodes, 2048 triangles, 64 sections, 64 retained modes,
16384 events and 600 seconds. File-declared mode/event ceilings may be smaller.
These are bounded workload sizes, not a calibrated time/memory or realtime
certificate. Modal reduction is offline and presently not cancellable. Audio
output staging is block-sized, uses the existing PCM16 encoder, counts clips,
and never peak-normalizes. Existing output/sidecar files are never overwritten.

Pressure uses the existing signed-area compact baffled monopole observation,
`rho A a(t - r/c) / (2 pi r)`, with `rho` and `c` from the ambient record. The
retarded acceleration is evaluated exactly inside the recorded exact-ZOH step,
not interpolated; before the first wavefront arrives the output is zero. There
is no radiation reaction load, atmospheric absorption, resolved directivity, or
whole-fluid energy balance in this input path. Large/noncompact radiators need a more complete acoustic model.
The source hash binds exact bytes; relocation does not change it, but whitespace
changes do. The sidecar is numerical/software provenance, not a certificate of
experimental validity, continuum accuracy, or source material truth.

Focused checks (native execution required):

```sh
cargo test -p fs-couple --test plate_render
cargo test -p fs-couple --lib render::plate::file
cargo test -p fs-couple --bin music_render
cargo test -p fs-couple --test music_render_plate
```
