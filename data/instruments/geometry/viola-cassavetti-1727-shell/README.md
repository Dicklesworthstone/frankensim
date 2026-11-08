# Viola body: the 'Cassavetti' proportions on the il Cannone tables

A viola body for the `material_strings` shell route. It is a composite, not a
scan. The il Cannone violin specimen tables
(`../violin-cannone-guarneri-1743-esrf2026/`) are mapped to viola dimensions:

- The back length and the three bout widths come from the Library of Congress
  record of the viola by Antonio Stradivari, Cremona, 1727, "Cassavetti"
  (<https://www.loc.gov/item/2022560103/>; the repo table
  `../loc-stradivari-family-bouts/dimensions.tsv`).
- The f-holes are placed from image 0015 of the same record's images
  (<https://www.loc.gov/resource/music.musihas-200154830/?sp=15>), John
  Montgomery's measured drawing of this viola, in millimetres.
- The arching heights, the back's interior thickness, the bass bar section, the soundpost
  and the bridge size are the Medium column of Hans Johannsson's viola table. That column's
  body length, 405-415 mm, holds this back.
- The top plate's interior thickness (47 printed belly points, 2.0 to 3.1 mm, mean 2.42 mm, as a
  Delaunay map with linear interpolation, `gibson-belly-points-normalised.tsv`), the rib height
  (36.3125 mm, the mean of eight printed heights), the rib wall (1.0 mm) and the back's edge thickness
  (4.75 mm) are the printed measurements of Stradivari's 'Gibson' viola, c.1734 (Roger Hargrave,
  The Strad, September 1986, pp. 113-117; ledger rows `hargrave-strad-1986-09-gibson-viola` and
  `hargrave-strad-1986-09-gibson-viola-belly-thickness`). They replace Johannsson's flank 2.6 and
  island 3.2 mm, rib height 38 mm, rib wall 1.2 mm and back edge 4 mm. The top plate's edge stays
  Johannsson's 4.0 mm: the Gibson's edge thickness is printed for the back only. One instrument, a
  same-period, same-size Stradivari viola, not the Cassavetti.

`PROVENANCE.tsv` gives the source of every mapped value and names each
ASSUMPTION. `specimen.tsv` gives the identity and the source ids of
`data/instruments/ledger.tsv`.

## What this body does not meet

- **f-hole checks.** The holes keep il Cannone's scaled shape, so these printed
  values are checks the placed holes miss (rows `f_hole_*` of `PROVENANCE.tsv`):

  | quantity | placed | printed |
  |---|---|---|
  | lower eye to the outline, bass / treble | 12.18 / 9.36 mm | 17.7 / 15.7 mm |
  | notch gap | 82.65 mm | 85.7 mm |
  | upper eye diameters | 5.56 / 5.90 mm | 8.0 / 8.0 mm |
  | lower eye diameters | 7.53 / 7.21 mm | 9.7 / 9.25 mm |

- **No printed bridge mass or bar length.** No source here prints a viola bridge
  mass. The bridge volume follows from Johannsson's bridge size on the il Cannone
  silhouette, which gives 2.81 g at the maple pack's density: an estimate, not a
  measurement. Johannsson's bass bar length cell is blank for Medium, so the length
  is il Cannone's through the plan scale. Both are ASSUMPTION rows.

## Files

- **Specimen tables**, in the il Cannone schema: `top-plate-*`, `back-plate-*`,
  `ribs-sections.tsv`, `bass-bar-sections.tsv`, `soundpost-sections.tsv`,
  `bridge-*`, `string-contacts.tsv`, `part-volumes.tsv`, `dimensions.tsv` and
  `PROVENANCE.tsv`, written by `cello_specimen.py ... instrument=viola`.
- **`shell-assembly-h6.5.tsv`**: the 6.5 mm shell table (37 components, 12087
  nodes, 18971 triangles, 75 ties; the generator's 0.2 mm thickness zones of the Gibson map), md5
  `c0736e8087d3`, written by
  `gen_violin_shell_assembly.py ... 6.5 neck=ellipse fanfix=move instrument=viola`. The same recipe without
  `gibson=` gives the earlier table (25 components, 11796 nodes, 19002 triangles, 48 ties, md5
  `0aa927b1bf74`, Johannsson's values).
  This is the finest step under the v4 body profile caps: the rib has 3056 of the
  3072 nodes a component may have.
- **`viola-{C,G,D,A}.assembly`**: the `assemble-shell` options, one bridge per
  string.
- **`assemble-shell-{C,G,D,A}.txt`**: what `assemble-shell` prints for each
  assembly.
- **`viola-shell-{C3,G3,D4,A4}.material-strings`**: the four bowed open strings
  played on this body. Their body and cavity records are the values in
  `assemble-shell-C.txt`.

The scripts are in `crates/fs-couple/examples/shell_body_modes/`. The string pack
is `data/matdb/seed-v1/viola-percival-2013-guettler-2012/`.

## The cards

Each card uses the records of the cello shell cards, with these viola choices:

- **Clock:** 48 kHz to 144 kHz, the violin pages' clock
  (`audio 144000 432000 4 128`). At 48 kHz and 64 bending modes, D4 refused on the
  coupled stiffness bound and A4 on the bow modes' time-domain admission; at 48 and
  32 modes both refused on the coupled stiffness bound. That bound does not depend
  on the time step. At 64 modes it is 8.65e4 rad/s for C3, 1.10e5 for G3, 1.64e5 for
  D4 and 2.17e5 for A4, against guards of 1.36e5 at 48 kHz and 4.07e5 at 144 kHz.
- **Modes:** 64 bending modes on every string, and 279 body modes, the whole 15 Hz-5 kHz window of this
  table with the audio-band maple (the render's own count at setup). The earlier body held 261 modes with
  the registered maple and 254 with the audio-band maple, and its cards kept the lowest 209.
- **String:** speaking length 0.375 m (Guettler) plus the 1/6 afterlength, 0.0625 m
  (Fouilhé and Houssay, SMAC 2013, p. 60), so the span is 0.4375 m and the bridge
  is at station 6/7.
- **Bow:** at beta 1/11, with the bow forces of the live viola cards; friction
  `rosin-wg2025-eq13`.

The clock did not brighten the low strings: over the sustain, the mean of H4-H10
stayed at -36.3 dB on C3 and went from -31.2 to -31.4 dB on G3, relative to the
strongest harmonic.

## The string pack

The pack is a stated approximation, not one set of measured strings:

- Mass per unit length and tension are Guettler's (mid-range mass; tension = mass x
  the printed speed squared, as on the live cards).
- Diameters are from a different string set, Percival's viola-I (2013, Table 3.2),
  of unknown make, counted from photographs and believed accurate to 10 percent.
- **Bending stiffness is derived. No source prints a viola EI.**
  EI = B T L^2 / pi^2, with Percival's B, L and f0 (Tables 3.2 and 3.3) and
  T = Guettler's mass x (2 L f0)^2. The same derivation overstates Mansour's printed
  KS510 cello EI by 1.36x (A), 3.80x (D) and 5.54x (C), so these values are likely
  high by a similar factor.

## Status

These rows and the render table below are the earlier body's (Johannsson's thicknesses, the registered maple, 209 modes), graded on 2026-11-04 against the live viola pages of that day; the sample cards in this directory are the current body's and are not graded here, the site pages carry that body's grades.

Graded on 2026-11-04 against the live viola pages, from the renders below:

| note | inside, of 14 | live page before |
|---|---|---|
| A4 | 8 | 4 |
| D4 | 3 | 2 |
| G3 | 7 | 7 |
| C3 | refused | 5 (reduced body, rosin law) |

C3 is refused because its strongest line is at 395 Hz, not the named note. C3 and
G3 both collapse onto a line near 393 Hz; D4 and A4 do not. That fault is open. The
C3 card is kept so that the refusal can be reproduced.

## Reproduce the bytes

The table is architecture-dependent: these bytes come from Linux x86-64 with numpy
2.5.4 and triangle 20250106.

```bash
g=data/instruments/geometry
x=crates/fs-couple/examples/shell_body_modes
out=$(mktemp -d)
uv run --python 3.12 --with triangle==20250106 --with numpy==2.5.4 python3 -I \
  $x/cello_specimen.py $g/violin-cannone-guarneri-1743-esrf2026 "$out/spec" instrument=viola \
  gibson=$g/viola-cassavetti-1727-shell/gibson-belly-points-normalised.tsv
uv run --python 3.12 --with triangle==20250106 --with numpy==2.5.4 python3 -I \
  $x/gen_violin_shell_assembly.py "$out/spec" "$out/viola-h6.5.tsv" 6.5 neck=ellipse fanfix=move instrument=viola
```

Leave out `gibson=` to rebuild the earlier body (md5 `0aa927b1bf74`). The default f-hole seed is the one this table was made with, and the trap it can
hit does not fire at h6.5. `hole_seed=median-run` is only for building at h6: there
the default stops with "chain bar of top is not inside one component" and the
option clears it. At h6.5 the option moves the mesh slightly (md5 `76d2230b15a8`),
which would be a different body.

## Render a note

```bash
demo=$(mktemp -d)
d=data/instruments/geometry/viola-cassavetti-1727-shell
p=crates/fs-couple/examples/material_strings/packs
cp $d/viola-shell-*.material-strings $p/violin-spruce-norway-viala2020-shell.fsmatpk \
  $p/violin-maple-sycamore-sonderegger2013-shell-laemmlein.fsmatpk "$demo/"
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/viola-percival-2013-guettler-2012/manifest.tsv \
  --out "$demo/viola-percival-2013-guettler-2012.fsmatpk"
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/bow-hair-horse-pitteroff-woodhouse-acustica1998/manifest.tsv \
  --out "$demo/bow-hair-horse-pitteroff-woodhouse-acustica1998.fsmatpk"
cargo run --release -p fs-couple --example material_strings -- assemble-shell \
  $d/shell-assembly-h6.5.tsv $d/viola-A.assembly "$demo/viola-A.body"
cargo run --release -p fs-couple --example material_strings -- \
  "$demo/viola-shell-A4.material-strings" "$demo/a4.wav" --body-profile "$demo/viola-A.body"
```

Use the assembly of the same string: C for C3, G for G3, D for D4, A for A4. The
string pack's hash is `4ec14974a0c3208f`.

The earlier body's graded renders came from the binary with sha256 prefix `78740edfad4e3054`, in
about one hour each at 1.5 GB, four at once on a loaded host:

| note | profile md5 | WAV md5 |
|---|---|---|
| C3 | `eef6a8835399` | `7e699b9bdf35` |
| G3 | `3c0b4624bf36` | `32063d919150` |
| D4 | `68b07b8259a0` | `25273eac3636` |
| A4 | `01ea12e4fb4b` | `d9a731c7d638` |
