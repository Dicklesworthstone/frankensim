# fs-couple: legacy modules

This folder holds the fs-couple modules that no live code uses. Most are hand-coupled pairings that
came before the composition kernel (fs-ports, fs-system). They are kept as regression references: a
kernel graph is checked against the module it replaces (owner ruling, 2026-11-11). They are not
deleted, and new work does not build on them.

Each module keeps its module path. `lib.rs` loads it with `#[path = "legacy/..."]`, so
`fs_couple::rigid_damping_heat` is still `fs_couple::rigid_damping_heat`.

- A one-file module is `legacy/X.rs`.
- A module with submodules is `legacy/X/mod.rs`, with its files beside it.
- Three modules whose root names its children by `#[path = "X/..."]` are `legacy/X.rs` beside
  `legacy/X/`: bulk_viscous_heating, plate_impact and porous_momentum.

"No live code" was checked on 2026-11-13 at 7e31441ff2. With every module here switched off:

- fs-couple's library, its binaries and its live examples still compile;
- so do the crates and apps that depend on fs-couple: fs-airflow, fs-ascent, fs-cli,
  fs-euler-disc-e2e, fs-flutter-e2e, fs-opdsl, xtask, fs-wasm, engine-drum, engine-piano,
  render-gpu-piano and material_convective_network. Three targets of these failed, each with an
  error of its own that it also gives with every module switched on: engine-drum's `bench`
  example, engine-piano's `acoustic_medium` test, and the material_convective_network binary.

Five modules that no live code uses stay in `src/` all the same, because the instrument-claim
registry cites them as evidence or the source manifest pins them: bakeoff, pickup, speaker,
cinematic_stem and traveling_wave_line.

A module that live code starts to use moves back to `src/`.

The tests of these modules are in `tests/legacy/`, and their examples in `examples/legacy/`.
