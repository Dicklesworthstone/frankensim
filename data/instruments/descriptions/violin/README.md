# Material-driven violin renders

The default `render.sh` now runs the existing **four-string material bank** on
the il Cannone shell assembly and its enclosed air. It uses the installable
`material_strings` executable, which is also the research example of that name.
The earlier `arco-open-*.performance` cards remain as source records for bow
speed, force, attack/release and duration. Their single-string, zero-bending-
stiffness mechanics and one-way body oscillator are not used by this recipe.

```bash
cargo build --release -p fs-couple --bin material_strings
data/instruments/descriptions/violin/render.sh \
  target/release/material_strings /path/to/new-violin-run
```

Each `G3`, `D4`, `A4` and `E5` directory contains its complete prepared
performance, body profile, original binary material packs, native construction
and rendering logs, and `render.wav`. Set `NOTES="G3"` to render one note.
Existing note directories refuse; completed or partial work is never overwritten.
The native solver remains responsible for admitting all source data and physics.

The recipe preserves the committed 144 kHz clock, 30 body modes and 8 mm shell
mesh rung in `cannone-h8.assembly`. Bowing one string retains the other three
strings in the reciprocal mechanical assembly. Bow position is explicitly one
tenth of that string's speaking length from its bridge. Geometry, tension and
PCM scale do not change when another material pack is selected.

## Material and atmosphere experiments

The common launcher `apps/instrument-lab/tools/render_material_strings.py`
accepts any existing `.material-strings` performance and optional body profile.
It is usable for supplied cello, bass or other string assemblies without adding
instrument-specific physics or inventing their missing geometry.

Pass an exact shared gas record to change the atmosphere. For example, this
keeps the gas temperature and pressure fixed while adding helium:

```bash
NOTES=G3 data/instruments/descriptions/violin/render.sh \
  target/release/material_strings /path/to/new-helium-run \
  --gas-record 'mixture 293.15 101325 mole 3 N2 0.7 O2 0.2 He 0.1'
```

The one material-gas state reaches the existing string air loading/loss,
cavity and pressure-observation owners. Solid temperature remains the separate
`state` record. Thermochemical gas records can supply the native
`--gas-component SPECIES NASA9 TRANSPORT`, `--gas-water` and `--gas-phase`
options; no new gas-property calculator is introduced here.

Use `--pack REFERENCED_PATH=SOURCE_PACK` to substitute an explicitly compiled
material pack, including a body wood, core/section or interface. For example:

```bash
NOTES=G3 data/instruments/descriptions/violin/render.sh \
  target/release/material_strings /path/to/new-section-run \
  --pack violin-string-dominant-g-pickering1985-section.fsmatpk=/path/to/other-section.fsmatpk
```

The native source resolver requires every property that this selected
construction consumes. A replacement lacking its density, stiffness, loss,
torsion or directional properties refuses; no gold/silver/wood constants are
guessed. A wound string's effective section pack is not interchangeable with a
bulk alloy unless its construction assumptions actually apply. Original source
identities and queries are retained in the native log. The string is **not
retuned** after substitution, and the sound is not peak-normalized.

## Listen and inspect in the existing Instrument Lab

With the lab's usual NumPy, Pillow and Matplotlib dependencies installed, add
`--analysis` and the actual revision of the compiled renderer:

```bash
NOTES=G3 data/instruments/descriptions/violin/render.sh \
  target/release/material_strings /path/to/new-listening-run \
  --analysis /path/to/analysis --git-sha ACTUAL_RENDERER_COMMIT \
  --variant cannone-air
python3 apps/instrument-lab/tools/build_instruments.py --analysis /path/to/analysis
```

This calls the original `instrgrade/decompose.py listen` path, producing the
standard audio, spectrogram, scalogram and `pair.json` consumed by the app's
instrument cards. The distinct slot `G3-cannone-air` preserves earlier grades.
New configurations are explicitly ungraded: no recording is used to tune them
or to imply validation. The listening export retains the declared pressure
scale instead of independently normalizing every variant. Compare variants at
the **same** `--full-scale-pa`; increasing that common value can provide PCM
headroom without altering the mechanical calculation. Analysis still applies
its documented high-pass and onset crop; the original WAV stays in the run.

## Current physical limits

This is a reachable source-based assembly, not a claim of recording-equivalent
violin sound. The mesh and modal truncation need convergence evidence; the
existing compact volume/cavity observer is not a converged exterior wave solve.
The committed data combine il Cannone geometry with independently sourced
spruce/maple properties and Pickering's string sections, rather than one
co-characterized instrument. The default bow remains a point contact with the
base file's authored Stribeck law; neither thermal rosin nor hair noise is added.
The full source details and additional physical options remain in
`crates/fs-couple/examples/material_strings/README.md` and `BODY_PROFILES.md`.
