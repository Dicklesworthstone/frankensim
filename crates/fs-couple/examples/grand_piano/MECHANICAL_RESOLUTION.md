# Resolve the soundboard between ribs as well as along them

`piano_geometry` exports the existing source-derived Model D soundboard at an
explicit mechanical resolution. It calls the same `steinway_d` construction as
the renderer, with its existing outline, seventeen ribs, bridges, cut-off bar,
material laws and either existing per-key bridge map. There is no second
soundboard model, external mesher, assigned frequency or output equalization.

The original `--mesh-divisions` control subdivides transverse sweep rows but
does not uniformly resolve the spaces between them. The new optional
`--rib-normal-step-m` bounds the distance **perpendicular to the ribs** between
successive rows in metres. All original outline, rib and bridge feature
stations are retained exactly. The existing conforming triangulation attaches
beams to actual panel edges; it does not decorate a mesh with unattached ribs.

## Export and play

```sh
cargo run --release -p fs-couple --example piano_geometry -- \
  --divisions 32 --rib-normal-step-m 0.02 --bridge-map estimated \
  --panel-rotary-mass consistent-cubic --stiffener-mass consistent-eccentric \
  --output resolved-board.fsb

cargo run --release -p fs-couple --example grand_piano -- \
  --preset steinway-d --board-geometry resolved-board.fsb \
  --edge-cubic-board-mass --equilibrate-board-mass \
  --note 69 --duration 0.2 --render resolved-note.wav
```

`--bridge-map rt0425` instead selects the already-existing source table projected
to the reconstructed structural bridges, including its four extrapolated end
keys. Keep that choice unchanged across mesh comparisons; it is independent of
resolution. `--format obj --output resolved-board.obj` exports the same builder's
three-dimensional visualization. That OBJ is not a complete piano enclosure or
a certified watertight assembly of the overlapping rib and panel solids.

The optional `--panel-rotary-mass lumped|consistent-cubic` and
`--stiffener-mass lumped|consistent-hermite|consistent-eccentric` write the
existing physical inertia selectors into FSB. They do not alter section
thickness, density, elastic constants, supports or string tensions. The shared
board reader selects the corresponding mass before eigenanalysis, reduction,
bridge coupling and acoustic motion. Cubic panel inertia requires the renderer
flag `--edge-cubic-board-mass`, as in the example above. These are explicit
physical approximation choices, not calibrated voicings or output filters.
Omitting them preserves the original file bytes and reader defaults. OBJ
cannot record mechanical operators, so combining these selections with OBJ
refuses rather than silently ignoring them. Existing selector rows also refuse
instead of being duplicated or overwritten.

The exporter defaults to eight divisions, original row spacing, the historical
estimated key map and FSB output. Without the spacing option, divisions remain
in 4..32; explicitly resolved meshes admit 4..64. Every output must be fresh.
The complete geometry is constructed before output is opened. Existing files
are never overwritten; a write or disk failure may leave a partial new file.

The production exporter needs only Rust's standard library, so geometry export
can also run without materializing the simulation workspace dependencies:

```sh
rustc --edition=2024 -O crates/fs-couple/examples/piano_geometry.rs -o /tmp/piano_geometry
/tmp/piano_geometry --divisions 16 --rib-normal-step-m 0.0125 --output fine.fsb
```

## Physical scope and admission

This is refinement of the same **source-defined fields**, not subdivision of
an arbitrary imported FSB with a frozen piecewise-constant material assignment.
Panel thickness, rib taper and bridge sections are evaluated at the new element
stations by the existing field definitions. Their discretized mass and stiffness
can therefore change; no mass renormalization, material refit or new measurement
claim is made. The panel outline and physical bridge locations are unchanged.

Row counts are bounded before allocating the expanded station list or mesh,
using the maximum possible vertices per row. Requests beyond the existing
20,000-node/40,000-triangle geometry capacity refuse rather than coarsening the
requested resolution or dropping structural features. The original constitutive
constants, source labels, supports, damping and legacy constructors are retained.

Normal row spacing is **not** a bound on every triangle diameter or aspect ratio.
Acoustic quadrature, selected board frequency band, eigensolver budgets, retained
mode limits, string physics and audio rate are separate controls. A finer mesh
may reveal more modes and still be refused by the existing spectrum limits.
Do not raise a cap or declare a converged high-frequency response merely because
a geometry export succeeded. Bridge mobility, mesh/cutoff convergence and paired
listening comparisons remain required by `frankensim-5c0sj`.

## Focused checks

```sh
rustc --edition=2024 --test crates/fs-couple/examples/grand_piano/steinway_resolution.rs -o /tmp/piano-resolution-tests
/tmp/piano-resolution-tests
rustc --edition=2024 --test crates/fs-couple/examples/piano_geometry.rs -o /tmp/piano-geometry-tests
/tmp/piano-geometry-tests
cargo test --release -p fs-couple --example grand_piano resolution
```

The first target tests physical spacing, feature retention and capacity
admission. The second tests exporter options, physical inertia rows and exclusive
byte-preserving writes; full builder tests stay with `grand_piano` and its actual geometry
reader. These commands have not been executed in the development runtime,
which has no Rust toolchain. Independent Python geometry references covered
seven resolutions and both bridge maps, but do not constitute native Rust,
mechanical eigenanalysis, renderer or listening validation.

Publication continuation: the core resolution and its tests were recovered from
the saved patch without source changes. The exporter additionally connects the
already-existing panel and beam inertia choices, with two focused regressions.
Native execution remains unverified in this runtime; no workflow, dependency,
lockfile or default voicing changes are required by this implementation.
