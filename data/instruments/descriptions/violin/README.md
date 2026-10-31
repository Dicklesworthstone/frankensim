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

The `gas` state reaches the existing string air loading/loss and exterior
pressure-observation owners. The cavity inherits this state when the performance
has no `cavity_gas` declaration. An existing independent cavity atmosphere stays
unchanged when only `--gas-record` is supplied.

Use `--cavity-gas-record` to replace or add the interior medium of a description
whose native cavity supports separate media. The default violin's compact
`physical-neck` declaration uses one medium for both ends and therefore keeps
the same gas inside and outside. A differing medium requires a supplied spatial
cavity or an exterior-owned mouth, rather than a change to that default alone.

For example, the existing manufactured spatial-cavity fixture supports helium
inside and air outside. First compile its `wire.fsmatpk` and `closed-shell.fsmatpk`
sources as described in the material-string example README, then run:

```bash
python3 apps/instrument-lab/tools/render_material_strings.py \
  --renderer target/release/material_strings \
  --performance crates/fs-couple/examples/material_strings/regional-cavity.material-strings \
  --body-profile crates/fs-couple/examples/material_strings/closed-shell-cavity.body \
  --pack-root /path/to/compiled-example-packs \
  --out /path/to/new-cavity-gas-run \
  --gas-record 'mixture 293.15 101325 mole 2 N2 0.8 O2 0.2' \
  --cavity-gas-record 'mixture 293.15 101325 mole 1 He 1'
```

These are independently declared acoustic media at the same static pressure;
the native model does not add mixing or thermal equilibration through the
openings. A differing static pressure or a cavity declaration without a modeled
cavity refuses natively. Molecular-law data remain unchanged and must admit
the selected medium. Solid temperature remains the separate `state` record.
Both atmosphere records are retained in a prepared Instrument Lab experiment.
Thermochemical gas records can supply the native
`--gas-component SPECIES NASA9 TRANSPORT`, `--gas-water` and `--gas-phase`
options for the species in either region; the same native source registry
resolves both media.

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

The strings' tension loss eta_F stays a stated unknown, entered as 0 (`unknowns.tsv`,
`u-string-loss-spectrum`). The one violin row of
[`tension-loss.tsv`](../../../../crates/fs-couple/examples/material_strings/tension-loss.tsv),
Galluzzo, Woodhouse and Mansour 2017, measured a plain steel monofilament E. The
Dominant E here is wound steel (Pickering 1985, the table's note 1), and eta_F is
the loss of a construction, so that row is not applied.

## Resolve the finite bow contact width

`cannone-ribbon-10mm.material-strings` selects a **rigid, uniformly loaded
ribbon** on the same four-string, shell and cavity assembly. Its only change
from the original Cannone performance is `bow_ribbon 0.010 5`. The 10 mm
physical width comes from the violin-bow ribbon in Table I of
[Pitteroff and Woodhouse (1998), Part I](https://euphonics.org/wp-content/uploads/2022/03/Pitteroff_1.pdf),
also retained in `bow-hair-horse-pitteroff-woodhouse-acustica1998`. It is a
published bow-family dimension, not a measurement of il Cannone's or Iowa's
actual bow. The contact is perpendicular to the string; this is the width
along the string, without an inferred tilt correction.

Select this data through the existing violin renderer:

```bash
NOTES=G3 \
MATERIAL_STRINGS_PERFORMANCE="$PWD/data/instruments/descriptions/violin/cannone-ribbon-10mm.material-strings" \
data/instruments/descriptions/violin/render.sh \
  target/release/material_strings /path/to/new-ribbon-violin-run
```

The original bow card still supplies the speed ramp, **total** normal force
and release. Five midpoint cells divide that total force equally; each has
its own transverse/torsional contact map, slip velocity and traction in the
simultaneous string/body solve. The width is neither an output filter nor a
gain. Body, strings, atmosphere, receiver and pressure scale are unchanged.

The cell count is a numerical resolution, not a fitted physical parameter.
For a spatial refinement check, make separate copies with
`bow_ribbon 0.010 3`, `bow_ribbon 0.010 5` and `bow_ribbon 0.010 9`, select
each through `MATERIAL_STRINGS_PERFORMANCE`, and keep every other input and
the common pressure scale fixed. Compare unnormalized pressure histories;
neither three rungs nor a changing spectrum alone establishes convergence.

This resolves finite contact breadth only. It retains the original authored
Stribeck law, rigid bow support and uniform load distribution. It does not
claim measured hair compliance, a rosin thermal law, individual hair
engagement or improved acoustic fidelity. The source horsehair measurements
do not supply the full frequency/temperature-dependent constitutive data
needed to replace those assumptions. The point-contact recipe remains the
default and is available for the controlled comparison.

## Add source-selected hair compliance

`cannone-hair-compliance.material-strings` adds the published approximate axial
hair response from Pitteroff and Woodhouse (1998), Part I, pp. 546–549, to the
same four-string Cannone assembly. The source selects 50 engaged hairs of
nominal diameter 0.19 mm and free length 0.65 m, driven at midspan. Three
spatial cells retain their own slip and series elastic state in the same
coupled solver and common gesture timeline. This is a declared bow-family
model, not a measurement of Cannone's bow or the reference recording's bow.

```bash
hair_dir="$(mktemp -d)"
hair=bow-hair-horse-pitteroff-woodhouse-acustica1998
cargo run -p xtask -- matdb-pack \
  --manifest "data/matdb/seed-v1/$hair/manifest.tsv" \
  --out "$hair_dir/$hair.fsmatpk"
NOTES=G3 \
MATERIAL_STRINGS_PERFORMANCE="$PWD/data/instruments/descriptions/violin/cannone-hair-compliance.material-strings" \
data/instruments/descriptions/violin/render.sh \
  target/release/material_strings /path/to/new-hair-compliance-run \
  --pack "$hair.fsmatpk=$hair_dir/$hair.fsmatpk"
```

The selected source point is static and of unstated humidity; the file declares
its quasistatic continuation explicitly. Young modulus and geometry determine
about 61 kN/m total axial stiffness. The 10 kg/s total dashpot is the paper's
published reduced-model approximation, normalized to an equivalent material
viscosity with that same geometry. The paper explicitly extrapolates its
low-frequency damping measurements; this does not establish the real hair's
broadband response. Source records preserve that distinction.

Changing these source data changes contact force, the original string/body
motion and pressure directly. Neither compliance nor loss is applied as an
audio filter. Normal force is still prescribed uniformly, and the existing
Stribeck contact law is retained as an authored assumption. This option does
not resolve bow-stick flexure, distributed hair-wave motion, rosin temperature
or a measured rosin shear law. It is an additional physical rung for controlled
comparison, not evidence that the violin now matches a recording.

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


### Controlled bridge pulse and ringdown

`cannone-dominant-bridge-pulse.material-strings` supplies a prescribed
1 N bridge-0 load for the first 64 mechanical samples, then releases it.
At the existing 144 kHz clock this is a 0.444 ms rectangular pulse, with
impulse 0.000444 N·s. All bows have zero normal force and velocity, while all
four strings remain reciprocally connected to the body and enclosed air.
The pulse is an authored controlled input, not a measured hammer waveform or
an ideal impulse.

Use the common launcher directly so the released-bow schedule is retained.
With the ordinary source packs available in the recipe's existing pack
directory, render the bare assembly into a new directory:

```bash
materials=crates/fs-couple/examples/material_strings
geometry=data/instruments/geometry/violin-cannone-guarneri-1743-esrf2026
python3 apps/instrument-lab/tools/render_material_strings.py \
  --renderer target/release/material_strings \
  --performance "$materials/cannone-dominant-bridge-pulse.material-strings" \
  --pack-root "$materials/packs" \
  --shell-table "$geometry/shell-assembly-h8.tsv" \
  --assembly-options "$materials/cannone-h8.assembly" \
  --full-scale-pa 1 \
  --out /path/to/new-bare-bridge-pulse
```

For the coated comparison, use another new output directory and select
`cannone-h8-varnish-memory.assembly`. Supply the compiled varnish pack with
the same `--pack` binding shown above. Optionally add
`--coating-strain-guard 536870912 100000000000000`; these explicit example
allowances must admit the complete run. Keep the pulse, strings, geometry,
retained modes, atmosphere, receiver and common PCM scale equal. If the native
writer reports clipping, raise that common scale for both renders; do not
independently normalize their audio.

The original coupled solver advances the forced attack and the subsequent
free decay without resetting displacement, velocity or material histories
at release. Varnish-induced changes to mass, section stiffness and relaxation
can therefore change bridge motion, decay, body/air exchange and pressure in
one physical calculation. The native log retains the actual physical ports,
initial newtons, initialization and authored event count.

The finite pulse has a frequency-dependent spectrum, so its raw pressure FFT
is not itself a normalized bridge transfer. Use the loaded harmonic report
for supported Pa/N comparisons, and compare the actual transient waveforms
and ringdown at equal input and pressure scale. The rectangular input also
contains frequencies outside a finite measured varnish band; the hereditary
recipe's explicit continuation assumption still applies. This fixture adds
a reproducible experiment input, not a completed audio or convergence result.

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

The logarithmic v2 report includes native and selected receiver pressure in Pa/N, every
bridge mobility, shell volume velocities, coating heat and stored energy,
other component/material/connection/radiation loss, and solve/power residuals.
It uses the original coupled harmonic solver and retains the native pressure
observer from the time calculation. All retained strings remain connected;
bow/contact traction is released for this linear transfer. This is not a
harmonic approximation to nonlinear bowed playing.

Replace the final `0` with `all` to obtain every bridge-drive column in one
sweep. Rows are ordered by frequency and then declared bridge; `drive_bridge`
identifies the 1 N peak input. A single dynamic-stiffness factorization serves
all columns at each frequency. The library also accepts coherent complex
forces across bridges and computes their interference in power and coating
heat; those quadratic quantities cannot be reconstructed by adding isolated
column losses.

When the performance explicitly selects `shell_observer`, this report uses
its original reciprocal closed-shell BEM observation owner. With an explicit
`one_way_observer` and profile acoustic skin, it uses the original ruled
Burton–Miller producer, all retained body modes and the actual aperture-flow
shares. Both evaluate the producer at the forcing frequency, preserve raw
propagation phase, and share one acoustic sample across bridge columns. They
do not build a FIR, add playback latency or taper, or add another reaction or
loss debit. The report records the actual receiver, formulation, source
identities, pass band and loaded versus one-way mode scope, and includes the
three complex pressure gradients. The native pressure remains a separate
column. The default violin recipe still uses its compact pressure observer;
an exterior microphone requires these explicit performance/profile records.

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
without full direct force columns, and plate microphone selections refuse.
Receiver boundary/quadrature work is part of the same complete sweep allowance;
its frequency-sample limit and pass band also remain binding. Pressure retains
the selected observer's existing physical scope; neither the compact default
nor partial/one-way skin feedback is certified as a converged exterior-wave
solution.

Compare bare and coated runs with the same inputs and pressure scale,
including pressure phase, bridge mobility and dissipated power at common
supported frequencies. No independent audio normalization or sound-pressure
gain is introduced. Both response commands validate all samples before a
create-new report write.


### Drive several bridges coherently

The loaded harmonic command accepts physical force phasors at multiple
original bridge ports. Select `coherent` and repeat
`--body-response-force BRIDGE REAL_N IMAG_N` for each actuator:

```bash
target/release/material_strings \
  /path/to/new-varnished-violin-run/G3/performance.material-strings \
  --body-profile /path/to/new-varnished-violin-run/G3/body.profile \
  --loaded-body-response-only /path/to/new-varnished-violin-run/G3/coherent-response.tsv 101 2000 32 coherent \
  --body-response-force 0 0.01 0 \
  --body-response-force 1 -0.01 0 \
  --body-response-strain 536870912 1000000000000
```

This example applies opposite 0.01 N peak forces at declared ports 0 and 1.
Their original positions and directions determine the resulting bridge
motion; the command does not infer a torque or move a load to another
footprint. Real and imaginary components use `exp(-i omega t)`. A
nonzero imaginary component selects relative phase. The force pattern stays
constant across the frequency sweep.

The complete force vector enters one original coupled solve at each frequency.
Pressure, coating dissipation, stored energy and strain are evaluated from
that combined response. In particular, adding separate unit-column loss
values would omit the interference terms. All original strings, cavity
loading, material histories and the selected receiver retain their existing
owners. Duplicate actuator records add in authored order; all unmentioned
ports have zero force. The source validates actual port existence.

Coherent reports use **loaded v4** and record the summed physical force
vector in metadata. They report velocity in m/s, pressure in Pa, power in W,
stored energy in J and receiver gradients in Pa/m. The strain column
`coating_strain_sufficient_common_amplitude_scale` is a dimensionless
multiplier of the entire authored force pattern, preserving its relative
phases. It is not a force in newtons. The individual-film scale columns have
the same meaning. Existing numeric-bridge and `all` selections keep their
unit-force v2/v3 formats.

The common launcher can perform this analysis directly while staging the same
geometry and source packs. With the ordinary packs already available:

```bash
materials=crates/fs-couple/examples/material_strings
geometry=data/instruments/geometry/violin-cannone-guarneri-1743-esrf2026
python3 apps/instrument-lab/tools/render_material_strings.py \
  --renderer target/release/material_strings \
  --performance "$materials/cannone-dominant-g3.material-strings" \
  --pack-root "$materials/packs" \
  --shell-table "$geometry/shell-assembly-h8.tsv" \
  --assembly-options "$materials/cannone-h8.assembly" \
  --loaded-response-only 101 2000 32 coherent \
  --body-response-force 0 0.01 0 \
  --body-response-force 1 -0.01 0 \
  --body-response-strain 536870912 1000000000000 \
  --out /path/to/new-bare-coherent-response
```

The new directory contains the staged performance, body profile, material
packs, `loaded-response.tsv` and `loaded-response.log`. For the coated
comparison, select `cannone-h8-varnish-memory.assembly`, supply its compiled
varnish pack with the earlier `--pack` binding, and use another new output
directory. Keep the force vector, frequency grid, geometry, retained modes,
atmosphere and receiver equal. The allowances above are explicit examples,
not runtime estimates.

This selects a steady, released harmonic experiment. The performance's bow
schedule and transient force events are not replayed. Playback decimation,
contact-thermal and stop profiles, the transient strain guard, and audio/lab
exports are incompatible with this launcher mode. Harmonic strain bounds use
`--body-response-strain`. Source frequency support and the existing retained
linear model limitations remain binding; no full-violin convergence or
historical-recipe fidelity follows from this interface.

### Check coating applicability while rendering sound

The native renderer and common violin launcher accept
`--coating-strain-guard MAX_BYTES MAX_WORK`. For an already staged coated run:

```bash
target/release/material_strings \
  /path/to/new-varnished-violin-run/G3/performance.material-strings \
  /path/to/new-varnished-violin-run/G3/guarded.wav \
  --body-profile /path/to/new-varnished-violin-run/G3/body.profile \
  --coating-strain-guard 536870912 100000000000000
```

These are explicit example allowances, not a measured runtime or a guarantee
that an arbitrary geometry fits. The renderer admits the complete initial
state plus every mechanical endpoint before creating output; no film, mode or
endpoint is dropped to fit. Add the same option to `render.sh` to collect the
per-film peak bounds, limit ratios and sample locations in each `render.log`.
The guard also runs before native pressure decimation and with contact heat.
It is currently a native playback option; the launcher refuses combining it
with a lab-experiment export whose packet cannot retain the setting.

An unresolved authored small-strain envelope stops the render, retains its
diagnostic and prevents continuation. Already accepted mechanics and any
incomplete WAV prefix remain explicit. The source envelope is a modeling
assumption, not a measured damage threshold, and a conservative bound can
refuse without proving actual overstrain. The check observes the represented
shell displacement, including retained static offsets; it does not establish
intersample, omitted-mode, continuum or missing-prestress accuracy. Passing
does not establish a perceptual improvement or historical violin fidelity.

### Inspect coating strain and the supported drive scale

The loaded and unloaded response commands can now report local film-strain
bounds. For a complete bridge matrix with this observation:

```bash
target/release/material_strings \
  /path/to/new-varnished-violin-run/G3/performance.material-strings \
  --body-profile /path/to/new-varnished-violin-run/G3/body.profile \
  --loaded-body-response-only /path/to/new-varnished-violin-run/G3/strain-response.tsv 101 2000 32 all \
  --body-response-strain 536870912 1000000000000
```

The extra option declares observer bytes and total preparation/sweep work.
It prepares the original MITC observation once per coated triangle and reuses
it at every frequency. Actual coupled body coordinates feed the observer;
no bridge-velocity reconstruction or different shell basis is introduced.
The full triangle and signed film thickness participate, with coherent phase,
plane-stress thickness strain and correctly converted tensor shear.

The strain-enabled loaded v3 report retains pressure, mobility and losses,
and adds each film's Frobenius bound, tensor-component bounds, ratio to its
authored envelope, sufficient force scale and ruling bound coefficients.
`coating_strain_sufficient_peak_force_N` gives the common sufficient scale of
that row's 1 N drive under these film bounds. A bound ratio above one does not
prove material failure: the spatial bound is conservative and the source limit
is an authored small-strain assumption. The ruling Bernstein coefficient is
not an observed maximum location. `inf` means no finite cap from zero film
strain or absence of films, not unlimited instrument strength.

For a chosen common force factor, linear pressure and displacement scale with
that factor, while cycle-mean power and heat scale with its square. The same
source bands and all other model applicability limits still matter. This
allows pressure comparisons to carry an explicit coating-amplitude check;
it does not establish the instrument's overall allowed sound level.

These are bounds on the retained linear harmonic increment, subject to the
original floating-point/condensation roundoff. Unmodeled preload/cure strain,
transients, missing modes, continuum error and 3D interlaminar recovery remain
outside their scope. Evolving strain enforcement during nonlinear playback,
moisture/cure/aging effects, penetration and frequency-dependent Poisson ratio
remain outside these coating options.

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
