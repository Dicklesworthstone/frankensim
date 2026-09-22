# Crowned, reinforced piano soundboards

`piano_board_import import-crowned` preserves the selected OBJ midsurface's
actual heights and prepares a **3-D orthotropic CST/DKT shell**, reinforced
by bonded eccentric ribs and bridges. It uses `fs-plate::shell`, the existing
sparse/modal owners, and the existing physical piano runtime. No new oscillator
frequencies, output EQ, samples, or mesh processing in the audio loop.

## Import and play supplied crown geometry

The existing OBJ material sidecar remains the input. Select the actual
soundboard midsurface, not the outside of a solid cabinet. Triangulate curved
faces before importing: the shared OBJ reader intentionally refuses nonplanar
n-gons rather than guessing their interior surface.

```sh
cargo run --release -p fs-couple --example piano_board_import -- \
  inspect piano.obj
cargo run --release -p fs-couple --example piano_board_import -- \
  import-crowned piano.obj materials.fspi crowned.fss
cargo run --release -p fs-couple --example grand_piano -- \
  --preset steinway-d --board-geometry crowned.fss --render piano.wav
```

The main `grand_piano` CLI accepts both flat FSB and crowned FSS through
`--board-geometry`. With `--preset steinway-d`, the supplied file overrides
**only the board**: source strings, per-key felt and shanks remain enabled.
It skips generation of the preset board and never falls back to that board
when supplied geometry fails. Without a preset, existing custom-scale/hammer
behavior is unchanged. Modal CSV `--board` still excludes a preset.

All existing main-CLI controls compose with the supplied shell: `--midi`,
`--midi-half-pedal`, `--performance`, `--concert-pitch`, `--raw-tensions`,
`--scale`, `--hammers`, `--dampers`, `--microphone`, `--board-band-hz`, `--modes`,
`--sample-rate` and `--substeps`, under their existing admission rules. This
allows measured/custom hammer and damper cards to be used with actual crown
geometry rather than being restricted to the convenience renderer defaults.
For example, with the named geometry and performance files already supplied:

```sh
cargo run --release -p fs-couple --example grand_piano -- \
  --preset steinway-d --board-geometry crowned.fss \
  --midi performance.mid --midi-half-pedal --concert-pitch 442 \
  --duration 30 --modes 128 --render performance.wav
```

The main preset normally adjusts physical string tension to A4=440 Hz (or the
supplied concert pitch); `--raw-tensions` retains source tensions. The existing
`piano_board_import render-steinway crowned.fss piano.wav 6 [performance.mid]`
convenience path also accepts crowns, but preserves raw source tensions and
its fixed budgets. Both retain all 88 source courses, including silent strings;
each required bridge station must exist. `--note` alone does not remove other
strings from the resonator. Mesh-generation and geometry-export controls refuse
when the main CLI receives a supplied board, avoiding export of the wrong board.

The sidecar must explicitly declare `support,clamped` and `pretension,0`.
The board coordinate system is the sidecar's orthonormal `frame`; `units` is
metres per source unit. Each source node's height along frame `u cross v` is
retained. The source `flatness` field still controls numerical admission of
the temporary reference-plane chart; it is **not** permission to flatten the
structural crown. Heights must remain within +/-50 mm, and facet normals must
satisfy `normal.z >= 0.95` in the board frame. Steep shells, overhangs and solid
meshes refuse instead of silently becoming a shallow board.

A starting editable asset is available from the source-derived Model D:

```sh
cargo run --release -p fs-couple --example grand_piano -- \
  --preset steinway-d --dump-geometry model-d.fsb
cargo run --release -p fs-couple --example piano_board_import -- \
  export model-d.fsb model-d.obj model-d.fspi
```

Supply the actual crown heights in the OBJ while retaining node identities and
physical section assignments, then import with `import-crowned`. Merely raising
a flat panel by a constant changes its datum, not its curvature. Do not describe
an authored crown function as a measured Steinway surface. No external mesh or
specimen measurement is bundled by this feature.

## Structural coupling

Every facet retains its own supplied thickness, density, longitudinal and
transverse moduli, Poisson ratio and shear modulus. The sidecar grain angle is
in the board XY chart. Its direction is projected into each facet and rotated
into that element's tangent frame before stiffness assembly. It is not reset
to each triangle's first edge.

Beam reinforcement uses both transverse bending planes, axial stiffness and
Saint-Venant torsion. The legacy physical beam card supplies E, G, A, Iy, J,
eccentricity and density. This import path explicitly declares a rectangular
section reconstruction: `height = sqrt(12 Iy / A)`, `width = A / height`, and
`Iz = A width^2 / 12`. It preserves the supplied Iy and J. Arbitrary section
shapes are not inferred. The generic `fs-plate::shell::stiffened` API admits
independently supplied Iy and Iz for other section shapes.

The beam's centroidal endpoints are displaced from the shell midsurface along
area-weighted nodal normals. The **same** rigid-offset map
`displacement_centroid = displacement_surface + rotation cross offset`
transforms stiffness and inertia. This preserves reciprocal force/motion work,
the sign of eccentricity, axial/bending coupling, and physical beam mass.
Simply adding EAe^2 to a bending coefficient would miss these shell couplings.

The shell's six nodal coordinates include all three translations and rotations.
Clamped support elimination and modal analysis are performed by the existing
owners. The entire admitted frequency slice is retained; over-budget slices
refuse rather than losing modes silently. Source piano string forces consume
the resulting mass-normalized vertical bridge shapes through the existing
reciprocal moving-boundary coupling.

## Native format and explicit bridge bearing points

The native header is `frankensim-crowned-board-si-v1`. Rows follow the ordinary
FSB format, except for 3-D nodes and the required declaration:

```text
frankensim-crowned-board-si-v1
beam-section,rectangular-from-area-inertia
source,mixed,Your geometry and material attribution
node,0,0.0,0.0,0.0
# ... all nodes, triangle sections, fixed nodes, beams and bridge stations ...
support,clamped
pretension,0
damping,0.01
```

Node identifiers remain contiguous and zero-based. The ordinary triangle and
bridge rows refer to the same node/triangle topology. Optional
`bridge_arm,key,dx,dy,dz` gives an SI offset from the interpolated midsurface
station to the force-bearing point. Its vertical displacement is
`u_z + (rotation cross arm)_z`; the resulting single work-conjugate projection
is used for both force and motion. Absent arms retain the original midsurface
port, not an invented top-of-bridge point. A vertical force applied along a
purely vertical arm has no torque, correctly. Horizontal arms admit moments.

## Limits that affect the sound

This is a **linearization about a supplied initial geometry**. It does not
solve string downbearing, static deformation, residual crown stress, glue slip,
rim flexibility or prestress-dependent tangent stiffness. A nonzero pretension
request refuses. Initial crown and a preloaded equilibrium are different models;
see Mamou-Mani, Frelat and Besnainou, *Numerical simulation of piano soundboard
under downbearing*, JASA 123 (2008), DOI 10.1121/1.2836787.

The pressure observer remains a **projected flat-baffle approximation**. At each
actual shell quadrature point, signed normal displacement times true area is
preserved, including in-plane displacement contributions. The equivalent source
is relocated to z=0 for the existing Rayleigh observer. The structural mesh
remains fully 3-D. This acoustic source relocation is not an exact exterior
solution; source heights, finite-rim diffraction, both radiating sides, lid,
room scattering and air backreaction require a separate acoustic model.

The convenience `render-steinway` currently retains board modes through 400 Hz,
at most 24 partials per string, four mechanical substeps and 48 kHz output. The
main CLI exposes these controls, subject to the shared 128-board-mode and
512-partial-per-string ceilings and the output band. A larger budget is not
proof of spatial/time convergence. Both use the existing 2 Pa PCM full scale
and report clipping without normalization. No full-audible-band, calibrated-SPL,
real-time or measured-Steinway fidelity claim is made.

## Focused checks

```sh
cargo test --release -p fs-plate --lib shell::stiffened -- --test-threads=1
cargo test --release -p fs-couple --example piano_board_import -- --test-threads=1
cargo test --release -p fs-couple --example grand_piano -- --test-threads=1
```

The tests cover rigid-motion invariance, analytical eccentric beam energy,
physical mass, source-frame/unit covariance, actual crown-induced operator and
mode changes, signed normal volume conservation, preserved Model D ribs and
88 stations, unsupported-data refusals, and deterministic physical piano WAVs
that change when the supplied crown changes. Main-CLI regressions also compose
tuning, source shanks, custom felt and spatial dampers through persistent audio
blocks. Native execution status belongs to the test run, not to this description
of the authored tests.
