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

## Apply a measured varnish to the existing violin

`cannone-h8-varnish-memory.assembly` adds an explicitly authored 20 µm film
to the positive, outward faces of all exterior top, back and rib regions.
It uses the measured modern oil/copal specimen
`varnish-oil-copal-lammlein2021-s1` and advances its passive relaxation
histories through the existing coupled renderer. This is a mechanical
surrogate on the il Cannone model, not a recovered Guarneri or Stradivari
recipe. The uncoated assembly remains the comparison.

Compile the ordinary source pack, then select the coated assembly:

```bash
varnish_dir="$(mktemp -d)"
varnish=varnish-oil-copal-lammlein2021-s1
cargo run -p xtask -- matdb-pack \
  --manifest "data/matdb/seed-v1/$varnish/manifest.tsv" \
  --out "$varnish_dir/$varnish.fsmatpk"
NOTES=G3 \
ASSEMBLY_OPTIONS="$PWD/crates/fs-couple/examples/material_strings/cannone-h8-varnish-memory.assembly" \
data/instruments/descriptions/violin/render.sh \
  target/release/material_strings /path/to/new-varnished-violin-run \
  --pack "$varnish.fsmatpk=$varnish_dir/$varnish.fsmatpk"
```

Run the ordinary recipe into another new directory for the bare comparison.
Keep bow input, strings, mesh, mode count, atmosphere and pressure scale equal.
For another film family/specimen, edit a copy of the assembly options and
compile its corresponding pack; merely changing a file label does not change
its retained source identity. Multiple coating rows can represent separate
ground and outer films. Historical oil penetration is not added as a second
external thickness.

The film changes the physical laminate stiffness, mass, mass centroid and
material damping before body modes, bridge reactions, cavity coupling and
pressure are formed. A one-sided film retains membrane–bending coupling; it
is not a scalar stiffness multiplier or audio equalizer. The source conditioning
is 293.15 K and RH 0.35. The declared Poisson ratio 0.30 and small-strain
envelope 0.001 are modeling assumptions, not source measurements.

The hereditary recipe uses equilibrium film stiffness and adds the complete
projected positive Maxwell branches, including cross-mode terms. Film mass is
included once; wood/glue, interface and film losses retain separate owners.
Its explicit `prony-continuation` policy accepts the fitted passive model
outside the finite measured band during a transient, while retaining the
original source band and separate sampled storage/loss error certificates.
It does not claim measured behavior outside that band.

The current 16-candidate O1 fit has 13 positive relaxation times. Combining
equal times across regions before factorization requires at most 390
histories for the 30-mode body. The existing 4 GiB prepared resource
declaration admits that model analytically; the independent 512-history cap
does not guarantee that every larger selection fits the byte budget.

For the previous frozen-reference approximation, select
`cannone-h8-varnished.assembly` instead. That option freezes storage at 1 kHz
and checks each retained mode against the authored 25% film-storage drift
allowance. The allowance is not a bound on frequency, sound or total
prediction error. Its source-dependent modal loss is passive but does not
advance hereditary film states.

### Inspect the varnish's structural response

Add `--prepare-only` to the launcher command above to stage a new run, then
request the structural response using those staged files. No WAV positional
argument is needed for this analysis:

```bash
target/release/material_strings \
  /path/to/new-varnished-violin-run/G3/performance.material-strings \
  --body-profile /path/to/new-varnished-violin-run/G3/body.profile \
  --body-response-only /path/to/new-varnished-violin-run/G3/body-response.tsv 101 2000 128 0
```

The logarithmic TSV evaluates each film's original law at the common drive
frequency. It reports complex bridge mobilities, signed exterior/cavity
volume velocities, supplied power, separate wood/glue, interface and film
dissipation, and power/solve residuals. Bare and hereditary images are
supported; frozen-film images refuse this exact-source analysis. All samples
must pass source support and solve checks before a new file is created.

This command projects the same source body without constructing the string
bank, cavity dynamics, bow performer or prepared time propagator. The staged
performance still supplies body dimensions, material state, mode count and
input limits. It is the unloaded retained body's response, normalized to a
1 N peak bridge force, with no string or air-pressure feedback and no
per-element finite-amplitude strain certificate.

### Inspect the loaded bridge-to-pressure response

To retain string loading, cavity feedback, material histories and selected
exterior radiation reaction, analyze the complete connected bank without
preparing playback or creating a WAV:

```bash
target/release/material_strings \
  /path/to/new-varnished-violin-run/G3/performance.material-strings \
  --body-profile /path/to/new-varnished-violin-run/G3/body.profile \
  --loaded-body-response-only /path/to/new-varnished-violin-run/G3/loaded-response.tsv 101 2000 128 0
```

The logarithmic report includes complex receiver pressure in Pa/N, every
bridge mobility, shell volume velocities, coating heat and stored energy,
other component/material/connection/radiation loss, and solve/power residuals.
It uses the original coupled harmonic solver and the same native pressure
observer as the time calculation. All retained strings remain connected;
bow/contact traction is released for this linear transfer. This is not a
harmonic approximation to nonlinear bowed playing.

The loaded response evaluates the actual passive Prony film realization
used for playback; the structural report evaluates the original fitted film
law. Both enforce source support at every frequency. Coating heat comes from
its original history range, with no fabricated allocation among films whose
equal-time factors were combined. The query does not change accepted states
or histories.

This command uses the same physical bank construction as playback, with an
explicit harmonic-only policy in the original coupled owner. It creates no
coupled time-transition or heat-integral matrices and performs no held-step
bow-compliance solve. Shell reduction and original component initialization
still run; connections, materials and acoustic loads retain their existing
physical admission checks. Admitted shell BEM and skin loads remain available;
declaring the analysis resource budget does not select the prepared playback
backend's narrower acoustic compatibility. Rebuilding a source preserves the
analysis policy, and attempts to advance this bank's clock refuse.

The existing explicit `prepared-linear` record remains required as the mode,
byte and work allowance; its multiply-work allowance covers the entire sweep.
It does not request a time exponential in this mode. To export the same loaded
response before an ordinary WAV render, supply the new-WAV positional path and
use `--loaded-body-response` instead. That option still prepares playback.

Unsupported nonlinear/contact attachments, residual-flexibility extensions
without full direct force columns, and additional frontend microphone filters
refuse. Pressure retains the selected observer's existing physical scope;
the compact default is not a converged exterior-wave solution.

Compare bare and coated runs with the same inputs and pressure scale,
including pressure phase, bridge mobility and dissipated power at common
supported frequencies. No independent audio normalization or sound-pressure
gain is introduced. Both response commands validate all samples before a
create-new report write.

Evolving coating strain checks, moisture/cure/aging effects, penetration and
frequency-dependent Poisson ratio remain outside these coating options.

The [source-data guide](../../../matdb/seed-v1/VIOLIN_VARNISHES.md) distinguishes
historical observations from modern DMA fits. The executable
[section study](../../../../docs/VIOLIN_VARNISH_PHYSICS.md) quantifies directional
frequency and decay estimates using the same material and laminate owners.
No coated full-violin audio or mesh/mode/time convergence result is supplied
by the addition of this configuration.

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
