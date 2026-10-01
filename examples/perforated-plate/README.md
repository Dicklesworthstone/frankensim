# Corpus body: a perforated plate (mesher complexity gate)

A 126 × 126 × 2 mm aluminium plate with a 20 × 20 grid of 3 mm square holes,
rotated 35° about z and 21° about x, so every one of its 8,652 facets is
oblique. It dissipates 10 W and is cooled on every face by h = 10 W/m²/K.
It is the real-scale body that the scaling work on 2026-09-30 was measured
on (bead q61wp.46). Real CAD exports carry 10⁴ facets and more, and before
that day nothing in the corpus went past 256.

```bash
python3 examples/perforated-plate/generate_perforated_stl.py 20 /tmp/perforated.stl --rotate
cargo run -p fs-cli --bin frankensim -- --json import \
  examples/perforated-plate/perforated-plate-rotated.fsim /tmp/perforated.stl \
  "${WORK}/perf.db" --unit m --max-hole-edges 0
```

The generator is deterministic, so the project pins the STL's source hash
and the freshness lane generates the ~3 MB file instead of tracking it.

## What the probe found

The same family at N = 1..28, axis-aligned and rotated, run through the real
binary, exposed five costs the small corpus had hidden:

| Fix | Commit | Effect |
|-----|--------|--------|
| Mesh once per solve: uncertainty re-solves re-meshed identical geometry | b3e882d52 | 6–10× |
| Incremental edge index in segment recovery | 9b869b336 | 1.4–2× |
| `orient3d` decides axis-aligned coplanarity at the filter | 4237dc1c6 | 2–3× on flat bodies |
| Incremental face index in facet recovery | 766f9f002 | 2.4× at 16k facets |
| ASCII STL parsed at f64 instead of f32 | 2f49027e4 | rotated bodies mesh at all |

Conduction stage, axis-aligned: 5,644 facets 199 s → 5.3 s, and 16,588
facets in 23 s. Rotated: every size refused on f32-rounding slivers before;
now 5,644 facets take 10.4 s and this body 23.8 s. Every performance change
was checked to leave stage receipts byte-identical. The lane's gate is 120 s
for this solve, about 5× the measurement, because the host is shared.

The probe also turned a misleading refusal into an answer. Small plates run
hotter than the aa6061 card's 398.15 K validity ceiling, and they used to
refuse with "Armijo backtracking failed". They now refuse with
`cli-solve-conduction-material-span`, which names the span and the
temperature the solver was driven to.
