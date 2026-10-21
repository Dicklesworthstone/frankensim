# Closed guitar body modes from a shell assembly

**Superseded by the shared module** `crates/fs-couple/examples/shell_body_modes.rs` (table format
`frankensim.shell-assembly.v1`, with `drive` records). `gen_shell_assembly.py` now writes that format by default
(`format=guitar` writes the first format this example reads); on the same 24 mm Esteso geometry the shared example and
this one print identical frequencies, net volumes, Q columns and effective masses (largest relative difference 0 over
every compared column). Use the shared module; this file stays until its removal is approved. The generator also takes
`wall_layers_min=N` (at least N element layers through every deep brace and bar; 5 is converged to 0.3 percent, 1 locks).
`audit_boundaries.py` (with the owner-side `shell_body_modes/audit_welds.py`) checks a table for unwelded seams and cracked
components.

`guitar_body_modes` reads a `frankensim.guitar-shell-assembly.v1` table, builds one
`fs_couple::shear_shell::assembly::MaterialShellAssembly` from its components (top, brace
walls, back, ladder bars, rib strips), ties them with the rigid work-conjugate mortar, runs
the one assembled eigensolve (`reduce`), and puts the enclosed air on the same modes: one
collective compression spring on the shell modes' cavity volume rows plus the sound-hole
plug (the library's `HelmholtzModalCavity`). It prints the in-vacuo shell modes, the
air-coupled modes with the hole open, the same with the hole sealed, and the effective mass
of every coupled mode at three bridge string positions. No frequency, loss or gain is
edited anywhere.

```bash
cargo run --release -p fs-couple --example guitar_body_modes -- \
  data/instruments/geometry/guitar-esteso-1964-zimmer2024/shell-assembly-h12.tsv PACK_DIR \
  --modes 10 --fmax 600 --shapes top-back-shapes.tsv
```

`PACK_DIR/<name>.fsmatpk` must exist for each pack name in the table (`spruce`, `rosewood`,
`cedar`); build them with `xtask matdb-pack --manifest ABSOLUTE_PATH/manifest.tsv --out
ABSOLUTE_PATH/<name>.fsmatpk` from `data/matdb/seed-v1/guitar-*-shell/`.

Things learned by running it on the Esteso table:

- `--modes N` must not exceed the number of assembled modes inside `[--fmin, --fmax]` (the
  library certifies the whole window, then keeps the first N); a closed box has few modes below
  500 Hz, so ask for about 10 up to 600 Hz rather than 36 up to 800 Hz.
- The rigid tie requires the slave nodes to sit EXACTLY on the master skin: slave position =
  master node - side x t/2 x nodal director. A vertical offset under a tilted director (a gap of
  64 micrometres, well inside the declared tolerance) made the assembled mass matrix lose
  positive definiteness (`FS-SPARSE-DIRECT-PIVOT-BREAKDOWN` at shift NaN). `gen_shell_assembly.py`
  places bar and rib nodes this way; `check_assembly.py` must report gaps of about 1e-9 m.
- The cost is set by the number of modes in the window and by the node count: the 1647-node
  (24 mm) table takes about 1 minute single-threaded, the 5337-node (12 mm) table many minutes.

What it does NOT model (named in the table's `meta` rows and in `unknowns.tsv`): bridge,
neck and heel, tail and heel blocks, linings, rosette, fingerboard, strings, the central back
strip and patches, glue compliance, finish. Loss is not in the packs, so decay and Q are not
predicted. The wood constants are handbook values of other specimens; rosewood is a declared
surrogate.

Tools next to this file (Python, standard library plus numpy and Shewchuk's `triangle`; run
with `uv run --python 3.12 --with triangle --with numpy python3 -I ...`):

- `gen_shell_assembly.py SPECIMEN_DIR OUT.tsv H_MM [key=value ...]` writes the table from the
  retained plan tables (options: `back_arch=flat`, `brace_h_scale`, `back_bar_h_scale`,
  `top_t_scale`, `top_pack`, `back_pack`, `brace_pack`, `backbar_pack`, `rib_pack`).
- `check_assembly.py TABLE` prints facet-versus-director alignment, every tie's gap and the
  enclosed volume.
- `gen_wood_packs.py REPO_ROOT` writes the three shell-ready packs from the committed
  FPL-GTR-282 seed packs and the Bremaud, Gril and Thibaut (2011) ratios.
- `classify_modes.py SHAPES.tsv` counts antinodes and the top/back phase of each coupled mode.
