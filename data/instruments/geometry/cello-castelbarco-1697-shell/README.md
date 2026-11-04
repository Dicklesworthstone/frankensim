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
- **`cello-{C,G,D,A}.assembly`**: the `assemble-shell` options, one bridge per
  string. Their comments name the table by its earlier file name,
  `cello-h12.5-s8-v5.tsv`.
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
