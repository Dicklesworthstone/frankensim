# Cello body: the 'Castelbarco' proportions on the il Cannone tables

A cello body for the `material_strings` shell route. It is a composite, not a
scan. The il Cannone violin specimen tables
(`../violin-cannone-guarneri-1743-esrf2026/`) are mapped to cello dimensions:

- The back length and the three bout widths come from the Library of Congress
  record of the cello by Antonio Stradivari, Cremona, 1697, "Castelbarco"
  (<https://www.loc.gov/item/2022560102/>; the repo table
  `../loc-stradivari-family-bouts/dimensions.tsv`).
- The f-holes are placed from image 0021 of the same record's images
  (<https://www.loc.gov/resource/music.musihas-200154818/?sp=21>). It is a
  measured drawing, "Cello the Castelbarco Antonius Stradivarius 1697", by John
  Montgomery, December 2005, in caliper millimetres.
- The thicknesses, arching heights, ribs, bass bar, soundpost and bridge size
  are the 4/4 column of Hans Johannsson's cello table.
- The bridge has the typical cello bridge mass of 19 g (Zhang, Woodhouse and
  Stoppani 2016).
- The measured-outline body below takes the plan outlines of the back and of
  the top from the same record's photographs (back, front and strung front, images
  0001, 0002 and 0019 of resource music.musihas-200154818), scaled by the printed
  bouts, and places the f-holes by their traced shape (image 0021) at the notch
  height the photographs show. The default body keeps the stretched il Cannone
  outline.

`PROVENANCE.tsv` gives the source of every mapped value and names each
ASSUMPTION. `specimen.tsv` gives the identity and the source ids of
`data/instruments/ledger.tsv`.

## Files

- **Specimen tables**, in the il Cannone schema:
  - `top-plate-*`, `back-plate-*`, `ribs-sections.tsv`, `bass-bar-sections.tsv`,
    `soundpost-sections.tsv`, `bridge-*`, `string-contacts.tsv`,
    `part-volumes.tsv`, `dimensions.tsv` and `PROVENANCE.tsv`;
  - written by `cello_specimen.py`, then `cello_bridge_19g.py`.
- **`shell-assembly-h12.5.tsv`**: the 12.5 mm shell table (39 components, 92
  ties), written by `gen_violin_shell_assembly.py ... instrument=cello`.
- **`shell-assembly-h12.5-bar-r24.tsv`**: the table the live cello pages are
  rendered on. It comes from the same scripts with two changes:
  `bar_slide_mm=12.27` and `rib_layers=24` in place of 6. With
  `bar_slide_mm=12.27`, `cello_specimen.py` slides the bass bar 12.27 mm toward
  the neck as one piece, so its ends sit 83.5 and 83.4 mm from the top's ends
  instead of the mapping's 71.2 and 95.7. Alan Goldblatt's cello specification
  sheet prints 80 to 85 mm, and a maker's forum post (Cliff Green) gives 80.
  Both are makers' figures, not measurements of this cello.
- **`outline-measured-loc-photographs.tsv`**: the plan outlines of the back and of the top measured on the
  photographs (z above the bottom of the back, left and right x at 1 mm rows; the front to z 730 where the
  fingerboard hides the plate edge, the back to z 748 where the neck heel hides it). Derived numbers only; no
  image is kept. Scale: the mean of the three printed bouts (346 / 234 / 440 mm), to which the measured widths agree within 0.2 percent on the back photograph and 0.9 percent on the front.
- **`shell-assembly-h12.5-outline-r24.tsv`**: the table of the measured-outline body, from the same scripts with
  `outline=outline-measured-loc-photographs.tsv fholes=fholes-trace-image0021.tsv fhstop=photo` and
  `rib_layers=24` (38 components, 17733 nodes, 30711 triangles, 87 ties). The outline, the ring, the arching, the
  thickness zones, the bar, the post and the bridge all go through one map of five landmarks (bottom edge, lower
  ledge, waist, upper ledge, top edge) between the il Cannone back and the measured back. What a printed or
  photographed distance fixes is placed, not mapped: the traced f-holes (the notch line 19.0 mm above the front's
  lower ledge, which is 405.0 mm below the top edge: the printed stop of 419 / 420 mm disagrees with the
  photographs by about 14 mm and is not used), the bridge plane on the notch line, and the bass bar, slid until its
  ends are equidistant from the top plate's ends (83.4 mm, inside Goldblatt's 80 to 85 mm).
- **`PROVENANCE-outline-r24.tsv`**: the provenance of that body's specimen, with the landmark map and its check
  rows. The specimen tables themselves are not kept: check N24 makes them again.
- **`top-plate-fhole-walls-outline-r24.tsv`**: the f-holes' bores of that specimen, DECLARED and not measured
  (no CT of this cello exists). Each hole's bore is the mean normal thickness of the specimen's
  `top-plate-maps.tsv` along its traced ring (4.39926 and 4.42334 mm), and the pair's keeps the holes' parallel
  straight-bore conductance (4.41104 mm). Written by `cello_fhole_walls.py`, which also writes the raster the
  solve reads: the specimen's rows moved -0.15 mm in z onto the solver's grid, a translation that changes no
  conductance.
- **`top-plate-fhole-conductance-outline-r24.tsv`**: the output of `crates/fs-bem/examples/fhole_conductance.rs`
  on that raster and those walls: both holes together 100.968172 mm, the Richardson limit of subdivisions 1 to 8
  (760128 cells and 71 layers at subdivision 8). It is computed and not fitted. The plate of the solve is flat,
  rigid and infinite.
- **`shell-assembly-h12.5-outline-r24-bem.tsv`**: the measured-outline table with that conductance as its neck
  (`neck=bem`, the specimen made with `fhole_walls=` and `fhole_conductance=`). It differs from
  `shell-assembly-h12.5-outline-r24.tsv`, whose neck is the equivalent-ellipse estimate of 107.13 mm, only in its
  `neck` meta row and its `cavity` row (neck radius 55.37 mm in place of 58.60 mm, length 4.411 mm). The body
  profiles `assemble-shell` writes from the two tables are the same; only the cavity record differs.
- **`cello-{C,G,D,A}.assembly`**: the `assemble-shell` options, one bridge per
  string. Their comments name the table by its earlier file name,
  `cello-h12.5-s8-v5.tsv`.
- **`cello-{C,G,D,A}-bar-r24.assembly`**: the `assemble-shell` options the live
  pages' four profiles were made with, on `shell-assembly-h12.5-bar-r24.tsv`.
  Each is `cello-X.assembly` with one row changed, `projection-work
  17592186044416` (2^44) in place of 2^41, which `assemble-shell` refuses for
  this table. The comment above that row still gives the 2^41 reasoning: the
  files are kept byte for byte as run.
- **`assemble-shell-D.txt`**: what `assemble-shell` prints for the D assembly.
- **`cello-shell-D3.material-strings`**: the bowed D3 played on this body. Its
  body and cavity records are the values in `assemble-shell-D.txt`.

All the scripts are in `crates/fs-couple/examples/shell_body_modes/`.

## Reproduce the bytes

```bash
uv run --python 3.12 --with triangle --with numpy python3 -I \
  crates/fs-couple/examples/shell_body_modes/check_cello_shell.py "$(mktemp -d)/check"
```

Check B (this table) runs only on Linux x86-64, and check A (il Cannone's h6
table) only on macOS arm64: with the same numpy and triangle versions, the
triangulation differs by architecture.

Check B24 reproduces the table of the live pages, also only on Linux x86-64
(md5 b59a6549f59166e59b1ac8b2b0bb3970). Its three steps:

```bash
s=crates/fs-couple/examples/shell_body_modes
w=$(mktemp -d)
py="uv run --python 3.12 --with numpy==2.5.4 --with triangle==20250106 python3 -I"
$py $s/cello_specimen.py data/instruments/geometry/violin-cannone-guarneri-1743-esrf2026 "$w/raw" bar_slide_mm=12.27
$py $s/cello_bridge_19g.py "$w/raw" "$w/specimen"
$py $s/gen_violin_shell_assembly.py "$w/specimen" "$w/table.tsv" 12.5 neck=ellipse foot_band_mm=8 \
  fanfix=move crossing_snap_mm=2.5 bridge_dp_mm=0.5 rib_layers=24 instrument=cello
```

Check N24 reproduces `PROVENANCE-outline-r24.tsv`, `top-plate-fhole-walls-outline-r24.tsv` and the solver's raster on
any architecture, and the tables `shell-assembly-h12.5-outline-r24.tsv` (md5 acec5e667a6035f18b6f3f544bd662b5) and
`shell-assembly-h12.5-outline-r24-bem.tsv` (md5 babb46607d79881925866e273983ae0e) on Linux x86-64. On the body alone
(`shell_body_modes`, 256 modes, 15 to 5000 Hz, the registered spruce and maple packs) this body's lowest
open-sound-hole mode is at 82.87 Hz where the live pages' body has 82.12 Hz, its in-vacuo modes 1 to 100 lie
1.3 percent lower (median), and its highest purely bridge mode below 1 kHz is at 785 Hz (779 Hz): a body with
no strings, so these are checks against measured cellos, not matches.

## Render the D3

```bash
demo=$(mktemp -d)
d=data/instruments/geometry/cello-castelbarco-1697-shell
p=crates/fs-couple/examples/material_strings/packs
cp $d/cello-shell-D3.material-strings $p/violin-spruce-norway-viala2020-shell.fsmatpk \
  $p/violin-maple-sycamore-sonderegger2013-shell.fsmatpk "$demo/"
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/cello-ks510-mansour-2016-percival-2013/manifest.tsv \
  --out "$demo/cello-ks510.fsmatpk"
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/bow-hair-horse-pitteroff-woodhouse-acustica1998/manifest.tsv \
  --out "$demo/bow-hair-horse-pitteroff-woodhouse-acustica1998.fsmatpk"
cargo run --release -p fs-couple --example material_strings -- assemble-shell \
  $d/shell-assembly-h12.5.tsv $d/cello-D.assembly "$demo/cello-D.body"
cargo run --release -p fs-couple --example material_strings -- \
  "$demo/cello-shell-D3.material-strings" "$demo/d3.wav" --body-profile "$demo/cello-D.body"
```
