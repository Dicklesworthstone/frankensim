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

`--shapes` writes the displacement of the top and back nodes along their directors; `--node-shapes PATH` writes the three
global displacement components of every node of every part (m/sqrt(kg), mass-normalized) for the first 36 air-coupled modes.

`PACK_DIR/<name>.fsmatpk` must exist for the pack names in the table (`spruce`, `maple`):

```bash
cargo run -p xtask -- matdb-pack --manifest ABS/data/matdb/seed-v1/guitar-spruce-sitka-fpl-gtr282-shell/manifest.tsv --out ABS/PACK_DIR/spruce.fsmatpk
cargo run -p xtask -- matdb-pack --manifest ABS/data/matdb/seed-v1/violin-maple-sycamore-sonderegger2013-shell/manifest.tsv --out ABS/PACK_DIR/maple.fsmatpk
```

Tools next to this file (Python; run with `uv run --python 3.12 --with triangle --with numpy python3 -I ...`):

- `gen_violin_shell_assembly.py SPECIMEN_DIR OUT.tsv H_MM [key=value ...]` writes the table from the specimen
  tables (options: `scale_xz`, `thickness_shift_mm` as declared sensitivities; `bar`, `post`, `bridge` = 0 or 1;
  `zone_step_top_mm`, `zone_step_back_mm`, `min_patch`, `edge`; `fanfix=1` repairs the vertices where a component meets itself,
  which some scaled meshes need and the base table does not; `cavity_volume=exact` records the air volume enclosed by the
  inner skins, 2.0610 L for the base table, instead of the default facet sum, 1.9653 L, which is origin dependent for open surfaces;
  `uniform_t_mm` gives every plate zone one thickness, a weld diagnostic). With every option at its default the output is
  byte-identical to the committed `shell-assembly-h8.tsv`. Every choice that is not a table number is a `meta` row of the output.
- `check_shell_assembly.py TABLE name=kg/m3 ...` prints per component the facet-against-director alignment,
  triangle quality and mass, every tie's gap (about 1e-9 m when the generator is right) and the cavity volume, and
  refuses (as the owner's code does) a component that is not edge-connected, a vertex with several fans, or neighbouring
  facets whose normals differ by more than 0.9 in cosine.
- `classify_violin_modes.py TABLE OUTPUT.txt SHAPES.tsv` applies the identification rules of `body-mode-preregistration.tsv`
  to a run and prints each named mode against its band.
- `make_results.py OUT.tsv UNUSED RUN=OUTPUT.txt ...` and `make_outcome.py OUT.tsv RUN=CLASSIFY_OUTPUT.txt ...` write
  `body-mode-results.tsv` (the air-coupled mode list of each run) and `body-mode-outcome.tsv` (the registered rules per run);
  a run whose name starts with `cracked` is labelled superseded.
- `gen_violin_maple_pack.py REPO_ROOT OUT_ROOT [shear=fpl]` writes the shell-ready maple pack from the committed
  Sonderegger 2013 and Laemmlein 2020 packs; `shear=fpl` writes the declared shear sensitivity pack (Wood Handbook
  sugar maple G_LR and G_LT ratios times E_L) under the name `...-fplshear-shell`.

Modelled: top plate with its two f-holes, back plate, one rib strip, bass bar (wall shell on the top), soundpost
(tube shell between the inner skins), bridge (shell on the outer skin at its two feet), enclosed air with the two
f-holes as one neck of their conductance. Not modelled: neck and upper block, tail block, tailpiece, strings,
chinrest, linings and corner blocks, purfling, varnish, glue compliance (named in the table's `omitted` row and
in `body-mode-preregistration.tsv`).

The owner's code refuses, before any frequency: neighbouring facets or a facet and its director more than about 26 degrees apart (the post is a 16-gon),
a component that is not one edge-connected body with one fan at every vertex, a node that is a slave in one tie and a master in another, a node twice in one tie
(rings are cut into two chains) and a table of 72 components (63 and 61 were accepted; the exact bound is not stated here, and the scaled variants are merged to fit).

## Independent check of the shell and tie machinery (no specimen data)

`gen_check_strip.py REPO_ROOT KIND OUT.tsv NX [K=..]` writes free spruce strips (500 x 25 x 3 mm, grain along the length, the guitar spruce shell pack) as tables for this example and prints the
closed-form free-free Euler-Bernoulli bending frequencies f_n = (beta_n L)^2 / (2 pi L^2) sqrt(E_L I / rho A), beta_n L = 4.73004, 7.85320, 10.99561, 14.13717 (roots of cos x cosh x = 1), and a Timoshenko
estimate. KIND `plate` is one component, `split` the same strip as two components tied end to end at the midsurface, `tbeam` the strip with a 12 x 3 mm wall standing on its upper skin along the whole length,
tied on the skin (side +1, full rotations): the junction that every body table uses for ribs, bar, post and bridge. Run 2026-10-22 with the binary of 8a0fb7ec1 (drilling penalty 1e-3), trj, free in vacuo:

- plate, closed form 63.520, 175.095, 343.256, 567.420 Hz (Timoshenko estimate 63.48, 174.76, 341.97, 563.93): shell 20 / 40 / 80 / 160 elements along the length: 63.68, 63.53, 63.49, 63.48 (mode 1);
  176.80, 175.18, 174.82, 174.73 (mode 2); modes 4 and 5 of the 80 and 160 element runs (the third listed mode, 309 Hz, lies near the first-torsion estimate of a thin strip and is not compared): 342.09 and 564.08 (80), 341.72 and 563.01 (160). At 160 elements the shell is within 0.1 percent of the
  Timoshenko estimate on all four bending modes and 0.06, 0.2, 0.4 and 0.8 percent below the Euler-Bernoulli value on modes 1, 2, 4 and 5 (the shear and rotary-inertia effect it leaves out).
- split: the same strip as two components with one midsurface tie gives the frequencies of `plate` at 20 and 40 elements to the printed digit (63.68, 176.80, 306.58 and 63.53, 175.18, 308.38).
- tbeam: closed form (parallel axis, E_L for flange and web) 299.96 Hz for the first vertical bending mode; shell 304.52 (20), 303.36 (40), 298.29 (80, 4 layers in the wall), 297.18 (160, 4 layers),
  295.53 (80, 8 layers), 294.24 (160, 8 layers). The wall's shear flexibility is not in the Euler-Bernoulli value (a Timoshenko estimate with the wall as the shear area is 286 to 295 Hz), so the skin-tied
  T-junction sits 1.9 percent below the closed form and still falling slowly: not a 40 percent stiffness defect.
