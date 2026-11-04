# double-bass-johannsson-shell

The body of the posted double-bass notes E1, A1, D2 and G2: a shell model of a four-string double bass built from a maker's published specification table and mapped onto the il Cannone specimen tables. **It is not a scanned bass.** The outline, arching shape, bar profile and bridge silhouette are the il Cannone's (ledger row `esrf-il-cannone-ct-2026`), scaled to the table's body length (1130 mm) and three bout widths (530 / 380 / 650 mm). Every mapped value, its source and every assumption is in `PROVENANCE.tsv` (specimen) and `PROVENANCE-generator.tsv` (the generator's choices); the sources are in the ledger rows `johannsson-bass-specification-shell-model` and `brown-2004-double-bass-dissertation`.

## Files

| File | What it is |
|---|---|
| `specimen.tsv` | identity, rights and reference state (`frankensim.specimen.v1`) |
| `PROVENANCE.tsv` | every mapped quantity with its source or labelled assumption; written by `bass_specimen.py` |
| `PROVENANCE-generator.tsv` | the generator's options and choices (mesh size, rib layers, foot tolerance, bridge mass, packs, what is not represented) |
| `specimen-files.md5` | md5 of the 14 specimen tables `bass_specimen.py` writes |
| `bass-h19-s2d.tsv` | the shell-assembly table (`frankensim.shell-assembly.v1`): 48 components, 11358 nodes, 17777 triangles, 114 ties, cavity 116.0666 L |
| `bass-h19-{E,A,D,G}-lb.assembly` | the `assemble-shell` options of each note (window, projection work, the string's bridge point) |
| `profiles.md5` | md5 of the four body profiles `assemble-shell` writes from them (2.6 MB each, not stored) |
| `bass-shell-lb-{E1,A1,D2,G2}.material-strings` | the earlier performance cards: the hyperbolic friction law, an afterlength of a quarter of the speaking length, the bow at 1/11 of the open-string length. **The posted notes are not rendered from these** |
| `bass-shell-lb-eq13-al6-E1.material-strings`, `bass-shell-lb-eq13-al6la-{A1,D2,G2}.material-strings` | the cards the posted E1, A1, D2 and G2 pages were rendered from (2026-11-04): the measured-rosin friction row (`rosin-wg2025-eq13`), an afterlength of a sixth of the speaking length (E1) or a sixth x 116/115 (A1, D2, G2), each note's own bow force, 256 body modes, bytes as rendered |

The tools are in `crates/fs-couple/examples/shell_body_modes/`: `bass_specimen.py` (specimen) and `check_bass_shell.py` (the checks below).

## What the posted notes carry, and what they do not

- **Wood.** Top, bar and post are the Norway spruce shell pack, back, ribs and bridge the sycamore maple shell pack, in their `-declared-lowband` variants (compiled files in `crates/fs-couple/examples/material_strings/packs/`, sources in `data/matdb/seed-v1/`, ledger row `bass-declared-lowband-shell-packs`). The assembled body has modes at 94.8, 113.9, 123.0 and 126.0 Hz in vacuo, below the parents' bands (100 Hz spruce, 150 Hz maple). **Body wood properties below 150 Hz (maple) and 100 Hz (spruce) are a declared extrapolation of the measured packs, not measurements.** Every value is the parent's; only the band's lower edge moved to 90 Hz.
- **Strings.** Four-quarter Spirocore rows at a 1.100 m speaking length (the table's string length is about 1070 mm; the 30 mm is not modelled), with no intrinsic string loss (no loss is printed for these strings). The string pack `bass-spirocore.fsmatpk` named by the cards is `crates/fs-couple/examples/material_strings/packs/bass-spirocore.fsmatpk`, compiled from `data/matdb/seed-v1/bass-strings-thomastik-spirocore-orchestra-catalogue-2019` (ledger row `thomastik-spirocore-orchestra-bass-strings`).
- **Named positions and values that are not published:** the bass bar is moved 18.1 mm toward the centre to clear the bass-side f-hole (the one unpublished position); the bridge is a uniform shell of 152.6 g (no bass bridge mass is published); the interior plate thickness is the mean of three printed stations; the rib taper 220-160 is replaced by its midpoint 190; the afterlength of the earlier `bass-shell-lb-*` cards is the card template's 1/4 of the speaking length; the live cards use a sixth (E1) or a sixth x 116/115 (A1, D2, G2), the 116/115 being the cello's printed ratio, borrowed because no double-bass value is printed. There is no tailpiece.
- **Not represented:** see `not_represented` in `specimen.tsv`.

## Checks

```bash
python3 -I crates/fs-couple/examples/shell_body_modes/check_bass_shell.py [--bin material_strings] [--generator]
```

- **C** (numpy only, any machine): `bass_specimen.py` regenerates the 14 specimen tables from the committed il Cannone tables, byte for byte (`specimen-files.md5`), and the committed `PROVENANCE.tsv`.
- **T** (standard library): table counts and md5, the cards' common body record and cavity row against the table, the declared-lowband packs named by the cards and assemblies.
- **E** (`--bin`, `material_strings` from the R5C pin or later): `material_strings assemble-shell bass-h19-s2d.tsv bass-h19-E-lb.assembly OUT.profile` writes the profile of `profiles.md5`, and prints the body length, width and cavity numbers the card carries.
- **B** (`--generator`, Linux x86-64 only): the table from the regenerated specimen. The triangulation differs by architecture with the same numpy and triangle versions, so the table reproduces only on the architecture that made it; the check says so and skips elsewhere.

## Reproducing a note

```bash
material_strings bass-shell-lb-eq13-al6-E1.material-strings OUT.wav --body-profile bass-h19-E-lb.profile
```

with the profile from check E, the declared-lowband packs and `bass-spirocore.fsmatpk` beside the card, and the PCM scale of the card (8 Pa). The posted renders (E1, A1, D2 and G2 of 2026-11-04) ran on the R5C pin plus the friction code landed as 2d08240c9 and 1635e3998, with the body limited to its lowest 256 modes by that binary. The earlier renders ran on `material_strings` built from the R5C pin plus the band-refusal naming landed with this directory; that change alters only what a refusing run prints (checked on this body with every retained mode in band: the `shell_body_modes` output and shapes file are byte-identical with and without it; no full render was repeated on the unpatched binary).
