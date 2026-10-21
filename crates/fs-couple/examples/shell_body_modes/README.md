# Closed body modes from a shell assembly (il Cannone violin)

`shell_body_modes` reads a `frankensim.shell-assembly.v1` table, builds one
`fs_couple::shear_shell::assembly::MaterialShellAssembly` from its components, ties them with the rigid
work-conjugate mortar, runs the one assembled eigensolve (`reduce`), and puts the enclosed air on the same
modes: one collective compression spring on the shell modes' cavity volume rows plus the neck plug
(the library's `HelmholtzModalCavity`). It prints the in-vacuo shell modes, the air-coupled modes with the neck
open and sealed, and the effective mass of every coupled mode at the table's `drive` points. No frequency, loss
or gain is edited anywhere. It is the guitar example's path (`guitar_body_modes`) with the table's plates in
several uniform-thickness zones (a component's group is its name up to the first `-`) and the drive points
named by the table instead of fixed in the code.

```bash
cargo run --release -p fs-couple --example shell_body_modes -- \
  data/instruments/geometry/violin-cannone-guarneri-1743-esrf2026/shell-assembly-h8.tsv PACK_DIR \
  --modes 24 --fmin 20 --fmax 700 --shapes shapes.tsv
```

`PACK_DIR/<name>.fsmatpk` must exist for the pack names in the table (`spruce`, `maple`):

```bash
cargo run -p xtask -- matdb-pack --manifest ABS/data/matdb/seed-v1/guitar-spruce-sitka-fpl-gtr282-shell/manifest.tsv --out ABS/PACK_DIR/spruce.fsmatpk
cargo run -p xtask -- matdb-pack --manifest ABS/data/matdb/seed-v1/violin-maple-sycamore-sonderegger2013-shell/manifest.tsv --out ABS/PACK_DIR/maple.fsmatpk
```

Tools next to this file (Python; run with `uv run --python 3.12 --with triangle --with numpy python3 -I ...`):

- `gen_violin_shell_assembly.py SPECIMEN_DIR OUT.tsv H_MM [key=value ...]` writes the table from the specimen
  tables (options: `scale_xz`, `thickness_shift_mm` as declared sensitivities; `bar`, `post`, `bridge` = 0 or 1;
  `zone_step_top_mm`, `zone_step_back_mm`, `min_patch`, `edge`). Every choice that is not a table number is a
  `meta` row of the output.
- `check_shell_assembly.py TABLE name=kg/m3 ...` prints per component the facet-against-director alignment,
  triangle quality and mass, every tie's gap (about 1e-9 m when the generator is right) and the cavity volume.
- `gen_violin_maple_pack.py REPO_ROOT OUT_ROOT` writes the shell-ready maple pack from the committed
  Sonderegger 2013 and Laemmlein 2020 packs.

Modelled: top plate with its two f-holes, back plate, one rib strip, bass bar (wall shell on the top), soundpost
(tube shell between the inner skins), bridge (shell on the outer skin at its two feet), enclosed air with the two
f-holes as one neck of their conductance. Not modelled: neck and upper block, tail block, tailpiece, strings,
chinrest, linings and corner blocks, purfling, varnish, glue compliance (named in the table's `omitted` row and
in `body-mode-preregistration.tsv`).
