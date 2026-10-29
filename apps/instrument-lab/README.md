# Struck — FrankenSim instrument lab

A public demo of FrankenSim's physically simulated drum kit and concert grand.
Live at https://frankensim.org/music and https://frankensim-struck.vercel.app.

Two tracks, always labelled:

- **Offline reference** — clips rendered natively by the `percussion` and `grand_piano`
  examples in `crates/fs-couple/examples`, with every clip's command, render time and energy
  closure parsed from its own log. One fixed playback gain per comparison; nothing is
  normalized, equalized or reverberated.
- **Live** — WebAssembly engines that compute sound in the browser: `engine-kit` (the
  `fs-vkplate` drum-kit solver itself on a coarser grid, one Worker per piece) and
  `engine-piano` (a lean kernel that steps the grand piano's own equations for the sounding
  keys only, validated against the reference `Instrument`). `engine-drum` is the earlier modal
  drum (FE heads, eigen-solver, prepared linear impact engine), no longer on the page.

## Layout

| Path | What |
|---|---|
| `index.html`, `src/` | the page (Vite, TypeScript, three.js; no framework) |
| `src/live/` | live section UIs and typed clients for the engine workers |
| `public/live/` | engine wasm + module workers + the AudioWorklet sink (static, unbundled; a build serves them under `live/<hash>/`, see Rebuilding) |
| `engine-kit/`, `engine-piano/`, `engine-drum/` | the live engines (standalone cargo workspaces) |
| `render-gpu/` | offline kit renders on an NVIDIA GPU when one is present (one persistent CUDA kernel stepping the exported fs-vkplate instrument; batches, chunks, both GPUs), on the CPU otherwise; see its README for the measured agreement and speed |
| `public/live/mix-worklet.js`, `lane.js` | `fs-live-mix`: sums frame-stamped lanes from several engine Workers in lock-step with the AudioContext clock |
| `tools/strandbar.py` | builds the Strandbar studio loop and the live engines' score and levels (`data/strandbar.json`) |
| `tools/build_data.py` | turns a render campaign into `public/{audio,img,data,viz}` |
| `tools/diagnose.py` | traces each real-vs-ours gap to a link of the signal chain and scores ablation renders: `public/data/diagnosis.json` |
| `tools/verisim.py` | grades renders against real recordings: `public/data/verisim.json` plus `vs-*` spectrograms and level-matched clips |
| `tools/percussion-viz-export.patch` | inert-unless-`FSIM_VIZ_DIR` export of mode shapes and modal coordinates used by the solver-state replays |

`public/{audio,img,data,viz}` and `public/live/*/engine*` are generated and not committed.

## Rebuilding

The workspace needs its sibling repos at the `constellation.lock` pins. If the shared sibling
checkouts have drifted, build from an export with pinned siblings beside it
(`git archive HEAD` for frankensim, `git archive <pin>` for each sibling).

1. Build the examples (release) with the viz patch applied in the export, then render the
   campaign; the exact command of every clip is in its "Run evidence" drawer on the page.
2. `npm run data -- <campaign-dir> public`
3. `engine-kit/build-wasm.sh` (it also writes `engine_levels.json`), `engine-piano/build-wasm.sh`
   (run it in the export: it links fs-couple by path) and, for the older drum section,
   `engine-drum/build-wasm.sh`; they install into `public/live/`.
4. `npm run build`, then `vercel deploy --prod --yes` (prebuilt `dist/`; see `vercel.json`).

The live files are never bundled, and each Worker loads its glue, wasm and data relative to its own
URL. So a build (`vite.config.ts`) serves the whole of `public/live` under one directory named by a
hash of its contents, `live/<hash>/`, and the page's `LIVE_BASE_URL` points there: a page always gets
Workers, glue and wasm from one build, and `vercel.json` and `cloudflare/_headers` cache `live/`
immutably. The build fails if what it copied differs from what it hashed. A page left open across a
deploy asks for its old `live/<hash>/`, which is gone, and says to reload. The dev server serves
`public/live` at `live/`.

The page's data, audio, images, replays and MIDI rolls keep their paths and a one-hour cache, but a build
hands the page each file's content hash and `util.asset()` asks for `path?v=<hash>`: a file a deploy
changed gets a new URL, so no visitor pairs new code with an hour-old `strandbar.json` or recording.

## Native wind and string renderer

The instrument-family section can edit and render wind, resolved-flow and string
performances with their original material packs using the **native `music_render wind`,
`music_render resolved` and `material_strings` solvers**. Build the local executables and the app, then serve that app through the optional loopback
companion. From the repository root:

```sh
cargo build --release -p fs-couple --bin music_render --bin material_strings
npm --prefix apps/instrument-lab ci
npm --prefix apps/instrument-lab run build
node apps/instrument-lab/tools/serve-native-wind.mjs \
  --engine "$(pwd)/target/release/music_render" \
  --strings-engine "$(pwd)/target/release/material_strings"
```

Open `http://127.0.0.1:8787`. `FSIM_MUSIC_RENDER` and `FSIM_MATERIAL_STRINGS` can
supply the respective absolute executable paths; use the actual Cargo output
paths if your target directory differs. `--dist` selects an existing app build
and `--port` changes the loopback port. `--timeout-ms` sets the execution allowance
for each render, from 1 to 600000 milliseconds (default 60000). For a longer
resolved-flow calculation, use an explicit allowance such as `--timeout-ms 300000`;
the panel displays it and Cancel still terminates the native calculation.
Either executable can be configured on its
own. The panel reports availability for the selected performer and refuses if
that executable is absent. Starting without executables leaves the viewer
available. Static hosted builds do not acquire a native rendering service.

The editor starts from the repository's `reed-two-medium-tone-hole.performance`:
its internal mixture, exterior mixture, tone hole and axial gas change reach the
original source parser and physics. It is synthetic integration data, not a
measured clarinet or oboe. The starter selector also opens the canonical
`flue-nolle.performance`, whose `jet_transfer` record chooses the same native
Nolle/Price source. Its mechanical and output clocks remain part of the native
description. Both starters are numerical examples, not measured instruments.

The same panel accepts the original reed, flue, plate-valve and material-modal
valve source formats, including `frankensim-material-modal-dg-valve-v1`, and
`frankensim-material-strings-v1`, `-v2` and `-v3` with the string performer selected.
Use **Attach source files…** to supply the original binary material packs or
artifacts referenced by the performance. Each displayed path must match the
source exactly; edit it to retain a nested path such as `materials/reed.fsmcdpk`.
A different wall or reed pack therefore reaches the native material resolver
without translating its properties into browser parameters. The shared
`dg_wall` declarations use this route too.

The panel retains selected bytes and shows their SHA-256 identities. Each WAV
records the submitted file identities in its playback details; changing a
source path, file, performer, profile or performance marks the existing WAV as an older result.
Keep the declared pressure-to-PCM scale fixed when comparing material variants.
For a DG wind, resolved-flow or material-string performance whose source clock is above 48 kHz,
select **Use native decimation**.
This passes the native `--decimate` option: physical integration retains its original
timestep and the existing native downsampler produces the 48 kHz WAV. The adapter never
rewrites the source clock or chooses a lower-resolution physics model.
The Nolle starter already declares separate mechanical and output clocks,
which its native observer handles directly.
Playback never peak-normalizes or retunes the native result. Original model
limits and missing-source refusals still apply; attaching a pack does not
qualify an instrument as measured or realistic.

### Resolved gas fields coupled to a duct

Select **Resolved flow — physical gas fields and spatial duct** and open a
`frankensim-resolved-dg-performance-v1` performance. The local server invokes
`music_render resolved` with the original source bytes. This lane evolves the
declared compressible gas fields and their coupling to the spatial duct; it does
not add a jet-transfer source to the same flow. The native source describes the
grid, stationary solids, initial physical density/velocity/temperature, incoming
characteristic history, duct/material walls and pressure observation.

Change those records in the source editor, or replace the original wall and gas
packs through the shared attachment controls. Thermochemical compositions use
the same explicit species/NASA-9/transport associations as other native lanes.
The native gas owner resolves the changing field states from those sources;
the browser supplies no frozen replacement coefficients. `observation` selects
inlet diagnostic pressure or the native baffled receiver, and `observer_gas`
retains the declared propagation model. Keep the receiver and pressure scale
fixed when comparing two material or atmosphere experiments.

The **Resolved flow — physical blowing pressure** starter opens
`resolved-dg-pressure.performance`. Its `boundary pressure` records hold a
pressure offset in pascals relative to the original source-face gas. The
original characteristic reconstruction uses the returning wave to determine
normal velocity; changing gas composition, calorics or static pressure therefore
changes actual source flow at the same prescribed blowing pressure. The same
resolved field, DG feedback and exterior pressure observer remain active.
This is an ideal pressure boundary within the existing local-linear domain,
not finite airway storage or a measured flute. Explicit `boundary characteristic`
and legacy `boundary COUNT` remain incoming-wave controls in m/s. Optional
`dg_molecular` uses the same source-bound DG loss reader as the wind performers.

The **Resolved flow — upstream acoustic feedback** starter opens
`resolved-dg-upstream.performance`. `upstream_dg` attaches a second original
finite material duct to the source zone. Its pressure and flow are solved with
the downstream bore on the same physical interval, retaining both ducts' wave,
wall and terminal histories. `boundary upstream-pressure` declares the ideal
series pressure rise between the upstream duct and the resolved zone; it does
not prescribe the resulting face pressure. Changing upstream geometry or its
material response therefore changes the actual source feedback. The source work
has a separate signed account. The observer remains on the downstream duct.
This is acoustic coupling without mean mass/enthalpy advection into either duct;
the short closed-tube fixture does not qualify an anatomical airway or flute.

The canonical `resolved-dg.performance` starter is a small exploratory physical
field fixture. It does not describe a complete flute or reed instrument.
Native admission retains the stated grid spacing, substeps and physical clock;
the Lab never coarsens them to fit the execution allowance. Complete experiments
save and restore this performer, the unchanged source files and gas associations.
Temperature and composition changes can alter the original material-derived
physical timestep while the declared observation clock and source duration stay
fixed. The native interval projection records its clock and bandwidth limits in
the provenance; refine both clocks for the intended physical frequency range.
The primary resolved gas offers dry-air and source-backed thermochemical models.
For a model change, author the thermochemical `energy_solve` row in the source
editor, or remove it when returning to dry air; the controls preserve those
numerical choices instead of supplying implicit tolerances. Common-pressure
edits include both the resolved primary gas and its exterior region.
The existing native pressure receipt enables the same calibrated A/B comparison
as the other performers, including refusal of comparisons with mismatched scales
or clipped outputs.

### Open a complete material experiment

**Open experiment…** restores a performance together with its original binary
material sources, nested source paths, body/contact/thermal profiles, gas source
associations and native-decimation choice. **Save experiment** writes that
complete editable configuration as `*.native.json`. This is the existing
`frankensim-native-wind-input-v1` request packet used by the local renderer;
the browser and service share its reader. No material property or physical
record is translated into a browser-specific model.

Use **Edit text…** beside an attached profile to change its original body
geometry, finger/stop paths or contact-thermal controls. **Apply source edit**
keeps the file's path and profile association and marks previous audio as stale;
render again to hear the change. **Cancel edit**, or applying unchanged text,
keeps every original byte, including line endings. Apply or cancel before
saving or rendering so the experiment always contains the displayed sources.
Binary bulk, interface and model packs use **Replace…** at their existing path.

The material-string launcher can prepare this packet from the same source
assembly used for native research renders. For the committed four-string
il Cannone assembly, from the repository root:

```sh
NOTES=G3 data/instruments/descriptions/violin/render.sh \
  target/release/material_strings /path/to/new-violin-input \
  --prepare-only --decimate
```

Open `/path/to/new-violin-input/G3/experiment.native.json` in the panel. The
native shell converter still prepares the body from its source geometry;
`--prepare-only` skips audio rendering. The packet carries the resulting body
profile and all original material packs, so the browser needs no manual
reattachment of those files. Render it, keep the result as A, replace a material
pack at its existing path or edit the gas record, and render B. Save each variant
to reopen the same complete configuration later. The declared pressure scale,
tensions, mechanical clock and source bytes remain part of the experiment.

For other supplied string assemblies, use
`tools/render_material_strings.py --performance ... --out ... --prepare-only`
with the appropriate `--pack-root`, `--body-profile`, `--stop-profile` and
`--contact-thermal` inputs. A prepared body profile needs no executable during
this staging step; generating a shell profile still requires `--renderer`.
`--stop-profile` includes the original finger, fret and material bow-contact
data, including moving contacts. Material-derived normal contacts retain their
receiver bulk `.fsmatpk`, ordered interface `.fsintpk` and normal-law
`.fsmodelpk` sources in the prepared experiment. For a source using another
filename, bind that exact referenced path with
`--pack contacts/normal.model-data=/path/to/original-model-pack`; its bytes
still go to the original native decoder. No contact coefficient is inferred
by the launcher or editor. Source-backed gases retain every original
`--gas-component SPECIES NASA9 TRANSPORT` triple and require explicit
`--gas-water none|ID` and `--gas-phase frozen|liquid-water` choices for export.
`--lab-experiment` saves the packet alongside a normal native audio render.

Choose `--decimate` explicitly when the performance's mechanical clock is above
48 kHz. This reaches the existing native downsampler and leaves the mechanical
integration unchanged. Experiment files contain input data; native coefficient
admission, convergence limits and specimen qualifications still apply when
rendering. Opening or saving works on the hosted viewer too; rendering requires
the local companion described above.

### Flute and reed atmosphere propagation

**Regional atmosphere** edits the native performance's temperature (kelvin),
absolute pressure (pascals), relative humidity or declared composition directly.
Choose a region, change its values, then **Apply atmosphere** and render again.
The same controls work with material-string exterior and cavity gases, flute/reed
bore and exterior gases, DG chimneys, prescribed DG gas-field stations and
explicit cylinder/cone/hole/duct-section gas overrides. Every
other source record, including geometry, gestures, mechanical clock, observer
law and pressure-to-PCM scale, is preserved. Apply or discard pending atmosphere
edits before rendering, saving or switching performances.

Add or remove species with their native IDs and explicit mole or mass fractions.
Fractions must sum to one; changing the fraction basis changes the declaration,
without silently converting its numbers. Choosing **Attached thermochemical
species** retains the exact species IDs and the source associations selected
below. Attach every required species, NASA-9 and transport pack; the native
material owner resolves the physical state and enforces its source domains.
Choosing a different gas model is an explicit edit to the original source.
The browser supplies no density, sound speed, loss or pitch coefficients.

For a common static-pressure comparison, select **Apply this absolute pressure
to every editable region**. Only explicit pressure fields are changed; legacy
segment declarations continue to inherit ambient pressure and `exterior same`
continues to use its original alias. Spatial field stations inherit that common
pressure; their coordinates, temperatures and fractions are unchanged by the
pressure action. Pack references remain untouched. Boundary conditions requiring a
common pressure are still enforced by the native solver.

For a `dg_bore` with an explicit `gas_field`, select a `gas_station` to edit its
temperature and mole fractions independently. Station coordinates and the
header's species columns are retained exactly. Change the common pressure in
the enclosing `gas`/`ambient` record. Changes to species IDs, their column order,
or station geometry belong in the full performance editor, with matching edits
to all affected rows and original source associations.

The **Wind bore — temperature and composition gradients** starter opens
`flue-dg-field-molecular.performance`. Build and attach the six N2/CO2 source
packs using the [canonical axial-gas recipe](../../crates/fs-couple/examples/reed-axial-gas-field.md#build-the-supplied-source-packs-and-render),
enable attached thermochemical sources, and declare their species/NASA-9/transport
triples, frozen phase and no water species. The field reaches the native bore's
integrated pressure storage, local molecular reservoirs and wall diffusion;
the exterior atmosphere stays independent. The fixture's relaxation rates are
manufactured, and its jet is interim. It demonstrates source wiring, without
claiming measured breath-field or complete instrument accuracy.

The **Flute — resolved duct and material atmosphere** starter opens the original
`flue-dg-baffled.performance`. Its spatial duct, outlet and rectangular mouth
use the native solver, with separate bore and exterior gases. Native decimation
is selected for its 96 kHz mechanical clock.

The flute starter and `paired-dg-chimney.performance` declare
`observer_gas stokes-kirchhoff MAX_PATHS MAX_FILTER_TAPS MAX_KERNEL_TERMS`
after `observation`. This applies the existing material-dependent propagation
law to each physical source-to-microphone path, using the complete `dg_exterior`
state. Temperature, pressure and mixture changes therefore reach density,
sound speed, viscous and thermal absorption, and the associated dispersion.
The flute mouth and outlet retain their relative propagation delay. The three
integer limits bound preparation work; they are numerical budgets rather than
acoustic parameters.

For an explicit comparison, change only `stokes-kirchhoff` to `lossless`, keeping
the rest of the source and PCM scale fixed. The observer selection changes
exterior propagation; it does not add molecular-relaxation dynamics inside the
bore or qualify the interim jet model as a measured flute. The browser sends the
original records unchanged and adds no equalization or normalization.

The **Flue tube — molecular gas relaxation (exploration)** starter opens
`flue-dg-molecular.performance`. Its `dg_molecular iso9613-terrestrial-air`
declaration supplies the bore relaxation model; its `observer_gas molecular`
record explicitly selects the same physical law for exterior propagation.
Atmosphere edits retain those declarations. The ISO terrestrial-air law has a
bounded composition domain, so enriched oxygen or additional gases may require
an explicitly supplied compatible molecular law before the native solver can
admit the change. The browser neither substitutes a different law nor invents
relaxation coefficients. This is a generic flue-tube integration fixture with an
interim jet closure, not a measured flute or an instrument-fidelity claim.

The **Flute — Nolle jet, viscous windway and molecular gas** starter opens
`flue-dg-nolle-molecular.performance`. The characteristic-network and spatial
duct performers now share the original windway and jet selection path, so the
Nolle/Price transfer can drive the molecular bore without substituting another
source model. The windway uses the inlet gas's density and viscosity; the bore
and microphone paths retain their separately declared media. This 96 kHz
fixture enables native decimation. The Nolle source remains a small-signal
transfer model with its existing channel/profile domain; it does not establish
resolved turbulent-flow or measured-flute accuracy.

The **Conical flue — material gas dynamics (exploration)** starter opens
`flue-dg-conical-molecular.performance`. Its `dg_bore` and `dg_segment cone`
records retain the taper while the shared material law drives molecular and
thermal pressure dynamics at the original wall quadrature points. Edit the
cone dimensions in the source or use **Regional atmosphere** for the uniform
internal mixture and explicit exterior. These edits preserve the geometry,
selected relaxation law, excitation, mechanical clock and pressure scale.
Native decimation is selected for the source's 96 kHz mechanics.

This is a generic tapered-bore integration fixture. The spatial wall model uses
the local cylindrical approximation. Its uniform medium can be replaced by an
explicit source-backed axial field, as demonstrated by the field starter above;
live molecular-history transfer remains unsupported. The same gas law still
enforces its source domain. The example's interim
jet does not provide measured oboe, saxophone or flute qualification.

### Material reeds and compliant walls

The **Material reeds** starters open the original paired-blade descriptions.
`paired-walls.performance` uses `blade.matpack` and `wall.matpack`; each blade's
mesh, thickness regions and source material drive motion and swept flow while
the material walls return pressure to the same upstream/downstream solve.
Its 768 kHz mechanics already have a declared 48 kHz pressure observer, so the
starter leaves additional decimation off.

`paired-dg-chimney.performance` uses `blade.matpack` with the existing spatial
duct, side chimney and separately declared bore/chimney/exterior gases. That
source runs at 96 kHz, so selecting it visibly enables native decimation while
preserving the source bytes. The [material-modal source instructions](../../crates/fs-couple/examples/material_modal_valve/README.md)
give the pack-building commands and physical declarations. These are
manufactured integration specimens, not measured clarinet or oboe models.

The **Material reeds — full directional relaxation** starter opens
`paired-tensor-dg.performance`. Build `tensor-blade.matpack` from the supplied
`tensor-blade.tsv` and manifest using the same material-modal instructions,
then attach it. `loss tensor-maxwell` resolves the source strain bound and all
21 entries of each symmetric Mandel branch, including explicit zero claims.
The original material card also supplies equilibrium elasticity. Changes to
out-of-plane stiffness and relaxation times reach the complete constrained
material history, reed motion, returning bore pressure and microphone. The
example retains 12 area points and 72 history coordinates; refinement must
respect the exact owner's explicit capacity. Native decimation retains 96 kHz
mechanics. These manufactured source coefficients are not wet-cane measurements.

The **Material reed — three-dimensional material memory** starter opens
`solid-hereditary-dg.performance`. Build and attach `solid-valve.matpack` using
the material-modal instructions. `solid_memory` retains the original tensor
branches and relaxation times in the existing exact aperture solver, alongside
the complete selected equilibrium modes, pressure work, contact and bore state.
Relaxed and unrelaxed initial history are explicit data. Material and source-time
changes therefore alter the physical trajectory and microphone output through
the original constitutive memory. Point, projection, source-strain, temporal and
exact-state limits remain active; refinement that exceeds them is refused.
This tiny manufactured Tet4 specimen does not establish measured reed accuracy.

The **Material reed — original solid surface contact** starter opens
`solid-surface-dg.performance` with the same `solid-valve.matpack`.
`solid_closure spatial` selects original exterior faces, physical boundary
edges and independently declared nodal clearances. Every retained solid mode
changes the local signed gap, even when its average opening is zero. The same
native contact potential and positive P1 slit integral determine partial
closure, reaction and flow; the original pressure-work map and hereditary
material histories stay coupled to the bore and microphone. Change those
source data to compare local contact without changing renderer code. The
current geometry is fixed-reference Tet4 with parallel face normals; it does
not resolve moving-channel gas flow or establish measured reed accuracy.

The **Material reed — three-dimensional solid** starter opens
`solid-dg.performance`. Build `solid-valve.matpack` from `solid-valve.tsv`
and its manifest in the same material-modal example directory, then attach it.
The source declares the tetrahedral mesh, constrained degrees of freedom,
pressure-loaded exterior faces, material regions and tensor property names.
The existing solid owner derives all selected modes and the full coupled
viscosity; zero-opening modes remain in that same aperture/duct solve.
Material, geometry and support changes therefore reach the native microphone.
This Tet4 example is a small manufactured soft-solid specimen. Its material
response is frozen at the declared reference frequency, with explicit storage
and loss dispersion checks at retained modal frequencies. It does not evolve
hereditary solid histories or qualify measured cane, clarinet or oboe accuracy.
Native decimation retains the original 96 kHz mechanics.

### Material-string performances

The **Strings — force-loaded contact and nonlinear material** starter opens
`nonlinear-stopped-preload.material-strings`. Attach `wire.fsmatpk`,
`sleeve.fsmatpk`, `body.fsmatpk` and `stopped-pluck.stops`, then select that
stop file as the **Finger, fret and contact profile**. The same static solver
now includes both nonlinear string storage and the original compliant normal
potential, with reciprocal motion of any coupled contact receiver. Final
source, force, penetration and energy limits apply before initialization is
accepted. The supplied support is manufactured; source-bound `material-stop`
profiles use the same path. This specimen is not a measured instrument.

The **Strings — force-loaded nonlinear pizzicato** starter opens
`nonlinear-preload.material-strings`. Attach `wire.fsmatpk`, `sleeve.fsmatpk`
and `body.fsmatpk` using the [native string source instructions](../../crates/fs-couple/examples/material_strings/README.md).
Its `forces preload` declaration solves the full coupled static force balance
with geometric stretching, then releases the force at sample zero. The same
path accepts longitudinal-wave storage. Changing a material changes both the
initial deformation and subsequent motion through the original constitutive
law; no displacement shape, settling damping or retuning is inserted by the
browser. The source retains its bonded section, two transverse planes, torsion,
body/cavity and gas-loading states. This is a manufactured integration specimen.

The **Strings — full tensor body memory** starter opens the canonical
`tensor-body.material-strings`. Attach `wire.fsmatpk`, `tensor-blade.fsmatpk`
and `tensor-body.body`, then select `tensor-body.body` as the **Body
geometry/material profile** under **String specimen profiles**. The
[complete tensor-memory recipe](../../crates/fs-couple/examples/material_strings/BODY_PROFILES.md#complete-three-dimensional-tensor-memory)
uses `wire-manifest.tsv` for the wire and the material-reed fixture's
`tensor-blade-manifest.tsv` for the body; normalize the latter with
`--material-state synthetic solid manufactured 0` and keep its output filename
`tensor-blade.fsmatpk` so it resolves the performance's original source path.
The body retains all six tensor strain histories through the shared
plane-stress elimination, including out-of-plane effects on bending and
reciprocal string reactions. Change `branch_33` or `relaxation_time` in the
source and rebuild its attached pack to compare the resulting pressure.
This 48 kHz fixture needs no decimation. Its coarse, fixed-state cantilever
and material are manufactured data, not measured wood or a complete
violin, cello or bass.

The **Strings — material plate and exterior microphone** starter exposes the
canonical `rayleigh-observed.material-strings` performance. Attach its original
`wire.fsmatpk`, `directional-body.fsmatpk` and `orthotropic.body`, selecting the
latter as the body profile. The native receiver observes signed plate motion at
the specified microphone position and propagates each path through the declared
gas, while the same plate retains its reciprocal acoustic load. Atmosphere
controls edit that original gas record. This is a synthetic flat-plate specimen
in an infinite baffle, not a measured violin-family body.

The **Strings — independent cavity and exterior gases** starter opens
`regional-cavity.material-strings`. Attach `wire.fsmatpk`, `closed-shell.fsmatpk`
and `closed-shell-cavity.body`, selecting the last file as the body profile.
The performance's `cavity_gas` record supplies the existing interior acoustic
field independently of `gas`, which supplies the exterior and string drag.
Use **Regional atmosphere** to change either composition or temperature; use
the common-pressure option when changing pressure because this quiescent model
requires equal interior and exterior reference pressures. Omitting `cavity_gas`
retains the original inherited-gas behavior. Declaring it without the supporting
interior field refuses instead of silently ignoring the material change.
This manufactured closed-shell fixture is not a measured violin-family body.

Select **Strings — materials, bow/contact and shared body** in **Physical performer**,
or choose a canonical string starter. `pluck-and-bow.material-strings` connects
pluck and bow gestures to shared strings and a body; `relaxing-layered-strings.material-strings`
retains each bonded layer's original material relaxation spectrum. Attach the
original packs named in the selected source. The
[native string instructions](../../crates/fs-couple/examples/material_strings/README.md)
show how to build those packs from the supplied manifests. Both starters are
synthetic integration fixtures, not measured violin, cello or bass models.

String, body, region, fiber and contact-layer material paths accept full
`NormalizedMaterialCardPack` artifacts (`FSMCDPK`) as well as legacy
`NormalizedPack` claim artifacts (`FSMATPK`). These are the same full material
cards used by material reeds and duct walls. The native decoder reads the wire
format, independent of filename suffix. Use **Replace…** on an attached file to
keep its declared path when comparing a different material.

Full cards retain their original chemistry, phase, process, revision, source
identity and associated model records. The native string resolver uses the
declared temperature in kelvin and angular frequency in radians per second,
including compatible typed material axes. Missing claims, incompatible axis
quantities or unavailable constitutive laws still refuse; retaining a model
record does not enable an arbitrary new law. No browser conversion to copied
material constants is involved.

The **Material shear layer** starter uses the original
`shear-layer-bow.material-strings`. Attach `wire.fsmatpk`, `body.fsmatpk` and
`manufactured-shear-layer.fsmatpk`. The same native instructions explain how
geometric stretching or longitudinal waves retain the layer's elastic/plastic
history and compliant-drive storage in the coupled solve. A contact-thermal
profile additionally requires genuine caloric and transport claims for its
declared material; this isothermal example does not supply those missing data.

The optional **String specimen profiles** selectors associate exact attached paths with
the native body, stop and contact-thermal profiles. The companion forwards them
as `--body-profile`, `--stop-profile` and `--contact-thermal`. Attach every material
or geometry artifact referenced by a profile too; the same frozen byte map serves
the profile, material, contact, body and radiation loaders. These selectors add no
browser coefficients or alternate synthesis model. Strings use the same gas-source
controls below, with their atmosphere record kept in the performance.

Explicit string decimation runs the existing pressure observer through the shared
native downsampler and writes 48 kHz PCM, retaining the original mechanical clock,
gesture times, pressure scale and filter delay. The source duration must contain
a whole number of output frames. Contact-thermal profiles keep the complete
mechanical-rate feedback while this observer produces 48 kHz PCM. The native
CLI can also retain a full-mechanical-rate `--bow-trace` before decimation.
Bow traces and observer-store caches are not exposed by the web adapter.

### Source-backed gas mixtures

Enable **Use attached thermochemical gas sources** to use the same species
registry as the native `--gas-component` workflow. Attach each original species
pack, NASA-9 model pack, and **pure transport material-card pack**, then select
their three paths for each component. The native loader preserves the full
`NormalizedMaterialCardPack` identity and provenance; legacy `NormalizedPack`
transport claims also remain supported. Choose the phase policy and enter the canonical
water species ID, or explicitly enter `none` when water is absent.

Temperature, pressure, composition and canonical species IDs stay in the
performance's existing `thermochemical` records. The browser does not decode
material coefficients, assign gas identities from filenames, or rewrite those
records. The same attached species may be used by the inlet, regional and
exterior gas declarations; native source identity, thermodynamic, transport and
domain checks remain authoritative. The liquid-water option retains the
existing saturation admission check; it does not add condensation dynamics.

The companion forwards only the declared component triples as
`--gas-component SPECIES NASA9 TRANSPORT`, plus explicit `--gas-water`
and `--gas-phase` options. Every path must also name an attached file. The native
gas loader reads those components from the same frozen attachment byte map as
wall, reed, string and body materials. Unlisted files cannot fall through to the host
filesystem. Changing an association or policy marks the existing WAV as stale,
just as replacing the source pack bytes does.

Gas records stay inline in the performance; the panel does not pass a separate
`--gas-record` override.

The companion passes the exact UTF-8 bytes to the selected executable as
`music_render wind INPUT.performance OUT.wav --self-contained` or
`material_strings INPUT.performance OUT.wav --self-contained`, with a private working directory
and fixed arguments (plus explicitly selected decimation, gas associations and string profiles).
Attached files are written byte-for-byte beneath that
private directory at their declared relative paths. Each is explicitly allowed
with `--attached-source RELATIVE_PATH`; the native loader snapshots only these
files before parsing, and both material and radiation callbacks read only that
byte map. Unlisted paths stay refused, including when a referenced file exists
elsewhere on the host. Rebuild the selected executable for these options; older binaries
refuse rather than weaken the file policy.

The companion does not parse physical parameters, retune notes, change
PCM gain, synthesize fallback audio or run a shell. Playback becomes available
only after successful native execution and a complete 48 kHz mono PCM16 WAV.
Source or solver refusals remain visible; cancellation terminates the native
child and removes that job's temporary files.

The local adapter admits one render at a time, at most 1 MiB of source, 32 attached
files (8 MiB each, 32 MiB combined), 16 MiB of WAV, 256 KiB of process diagnostics
and 60 seconds per request. A versioned JSON/base64 envelope carries attachments;
its complete encoded size is bounded before decoding, and the original text
request remains supported. These are local
execution budgets, not physical acceptance criteria. Host/Origin checks and a
per-session request token restrict execution to the same local app. Only files
inside the selected built app are served; there is no command or filesystem API.

```sh
node --test apps/instrument-lab/tests/native-wind-server.test.mjs
node --test apps/instrument-lab/tools/test-native-wind-client.mjs
```

The server suite uses an explicitly synthetic subprocess to check byte/argument
preservation, real HTTP responses, refusal handling, cancellation and limits.
It verifies the adapter rather than native Rust physics or instrument realism.
Its executable-fixture checks require a temporary directory that permits local
execution. No rendered audio or built app assets are committed by this workflow.

## Grading against real recordings

The "Real vs ours" section compares renders with University of Iowa Musical Instrument Samples
(<https://theremin.music.uiowa.edu/MIS.html>, free to use without restriction): a Steinway Model B
(`Piano.{pp,mf,ff}.A3`, `Piano.mf.C6`) and an anechoic splash plus 13/17/18-inch crashes
(`splash.stick.normal.*`, `*crash.stick.normal.mf`). Convert each to 48 kHz float WAV. The campaign's
splash renders come from the kit's solver, `perc/vk-splash.wav` and `perc/vk-splash-soft.wav`
(`kit splash OUT.wav 2 2.6` and `kit splash OUT.wav 0.6 2.6`, see "The drum kit"); a pair's label
names what its render simulated. Then:

    uv run --with numpy --with pillow python tools/verisim.py <campaign-dir> <wav48-dir> public

Both signals of a pair get the same 40–70 Hz zero-phase high-pass and onset alignment. The script
asserts that both notes start within 12 ms of their windows. Every gap is reported beside a ruler
built only from real-vs-real distances: the Iowa piano across dynamics, and the other real cymbals
against the splash. The page's findings list is computed from `verisim.json`, so it cannot drift
from the numbers.

## Diagnosing the gaps

`tools/diagnose.py` (same arguments as `verisim.py`, whose measurement code it imports) asks
WHY each gap exists. It uses the A3 bridge trace (`--bridge-trace-csv`) to separate what the
board does from what radiation to the microphone does. It also scores ablation renders, which
switch one mechanism on at a time, whenever they are present in the campaign:

    # the same A3 mf heard where the Iowa microphones sat (stereo; see "Listening point" below)
    grand_piano ...a3-mf flags... --microphone 0.275,0.556,0.25 --microphone-right 0.867,0.147,0.233 --render a3-mf-close.wav
    # ...and by stereo pairs 0.1 m either side of the bass-side point, along x and along y
    grand_piano ...a3-mf flags... --microphone 0.175,0.556,0.25 --microphone-right 0.375,0.556,0.25 --render a3-mf-px.wav
    grand_piano ...a3-mf flags... --microphone 0.275,0.456,0.25 --microphone-right 0.275,0.656,0.25 --render a3-mf-py.wav
    # the former estimated bridge stations on the same mesh, with the same air and receiver
    piano_geometry --divisions 16 --rib-normal-step-m 0.03 --bridge-map estimated --output board-est16.fsb
    grand_piano ...a3-mf flags... --board-geometry board-est16.fsb --render a3-mf-estmap.wav
    # bare-board modes and per-key bridge displacements of the graded board (A3's driving-point mobility)
    grand_piano --preset steinway-d --board-geometry board-rt16.fsb --board-band-hz 1600 --dump-board board-modes.csv
    # longitudinal string waves on A3 only; every other key `linear,<key>` (88 rows)
    #   waves,57,166505,0.2,10000,8,0.01,0.01   EA = E*A, E = 200 GPa, A from the scale's own EI
    grand_piano ...a3-{mf,ff} flags... --string-stretching axial-a3.fspx --render a3-{mf,ff}-waves.wav

A card counts as confirmed only when its experiment moves the measured gap. A mechanism that has
moved into the graded renders is reported as fixed or partly fixed: the A3 renders
(`a3-{pp,mf,ff,fff}`) carry `--string-stretching axial-a3.fspx`, whose waves drive the bridge axially,
and when the campaign also has the same strokes on linear strings (`a3-{mf,ff,fff}-linear`) the card
compares the two. Every piano render the page grades (A3, C6), and every ablation above except in
the one thing it changes, couples each string to its bridge at the point RT-0425 publishes for it,
carries two-sided air loading and listens at the default point, 1 m above the middle of the board:

    piano_geometry --divisions 16 --rib-normal-step-m 0.03 --bridge-map rt0425 --output board-rt16.fsb  # 11,772 area samples
    grand_piano --preset steinway-d --board-geometry board-rt16.fsb --board-band-hz 2400 \
      --board-reduction pod,128,64,16,8000 --rayleigh-loading 2,90,2620,213 \
      --rt0425-hammer-stiffness --rt0425-hammer-dissipation --rt0425-string-damping ...

The hammers follow RT-0425's published law: its stiffness K_H acts on each unison string (the preset's
former felt conversion divided it across the unison, which its own comment calls too soft) with the
published dissipation R_H, and the strings take RT-0425's published losses. The same mf/ff/8 m/s strokes
with the former felt law, `a3-{mf,ff,fff}-legacyfelt` (the graded command without its two hammer flags),
let the ff-brightness card show what the hammer law alone changes; the card says so only when the
commands differ in nothing else.

The former estimated stations (`--bridge-map estimated`) were spaced evenly along each bridge. With
the published string lengths they put most agraffes over the soundboard itself (A3's 0.77 m string,
running straight back from its contact, would start 0.38 m inside the board's front edge) and
reversed the bass order. The reduced board keeps a coordinate at 2.60 kHz, so the air band runs to
2.62 kHz. Of the meshes tried, 12 to 15 divisions left the air cells under six per wavelength there;
16 divisions with a 0.03 m step is the coarsest that passes. The polarization ablation adds bridge
directions that reach 2.62 kHz, so it runs its air to 2.64 kHz (215 samples).

Listening point: the Iowa recordings used two Neumann KM 84 cardioids, "8 in above center bass
strings" (left channel) and "8 in above center treble strings" (right). In the board frame those are
the centroids of the bass and treble strings' speaking lengths, 0.203 m above the strings: the
published contacts, the agraffe line the treble strings imply, and the bass agraffes continuing that
line. The points are estimates, and ours are omnidirectional. The page does not grade there. Real
microphones at those points sit above the strings and the cast-iron plate, which covers much of the
board, usually under a lid; the model's receivers look at a bare board in an infinite baffle. Heard
there, the A3 partial envelopes (0.05-1.05 s) lie 11.3 dB (bass side) and 13.3 dB (treble side) from
the recording on the matching channel, against 10.7 dB at the default point. The real piano's own two
microphones differ by 8.4 dB. `a3-mf-close` keeps those positions as the receiver experiment.

The graded board is the 128-coordinate POD reduction (`pod,128,64,16,8000`; crates/fs-couple/examples/
grand_piano/BOARD_REDUCTION.md), the same as the live piano's largest board and the Strandbar studio stem.
Against its own 166-mode slice it is closer than the former pivoted reduction (`128,64,200,400,800,1600,2300`)
at A3's partials above 1.46 kHz: driving-point |Y| rms 0.59 against 1.83 dB, and microphone pressure
0.51 against 1.04 dB. The pivoted board missed partial 8 by +4.8 dB and partial 9 by -7.3 dB. Those
errors happened to suit the recording: on the pivoted board A3 mf's envelope distance was 7.7 dB, on
the POD board it is 10.7. A3's 2-8 kHz brightness moves by 0.4-1.3 dB across pp-fff, and the early
low resonance and the attack rise do not move (frankensim-a8lld).

The page's demonstration clips (Satie, the silent strings, the octaves, A0) are rendered with the same
board, laws and air. `tools/piano_performances.py OUT_DIR` writes their performances (and the A3/C6
strokes); its docstring lists each clip's receiver and duration. The board drawing comes from
`grand_piano --preset steinway-d --dump-obj board.obj --dump-geometry board.fsb --dump-scale scale-et.csv
--dump-board board-modes-400.csv`.

When the campaign also has the same A3 mf stroke without air on that board (`a3-mf-noair`), the
low-resonance and ringing cards compare the two. Each render first samples and fits the board's
complete 128-mode air reaction (15 minutes or more on the 16-division board; see the example's
RAYLEIGH_LOADING.md).

## The drum kit

Every kit piece comes from one solver, `crates/fs-vkplate` (its own workspace, std only; see its
README for the physics). Render each piece, then grade them against the Salamander Drumkit
(Alexander Holm, CC BY-SA 3.0; `archive.org/details/SalamanderDrumkit`, extract `OH/`):

    cd crates/fs-vkplate && cargo build --release --example kit
    for p in kick snare snare-off tom-12 tom-14; do kit $p $KIT/$p.wav; done        # 1-5 min each
    kit hihat-closed $KIT/hihat-closed.wav 0 0.6; kit hihat-open $KIT/hihat-open.wav 0 2.5
    kit crash $KIT/crash.wav 0 3; kit ride $KIT/ride.wav 0 3                       # ~20 min each
    uv run --with numpy --with pillow python tools/kitgrade.py $KIT $SALAMANDER/OH public

The five drums render from the sourced material cards by default: PET head film, hickory stick,
the film loss a stated unknown, and the kick's preset felt beater. See the crate README's usage.
Add `--preset-heads` for the recording-fitted preset heads, the drums' earlier default.
The card-built drums carry the head-radiation load, which costs time. On ts1, 0.3 s of tom-12
took 20 s on cards against 2.4 s preset, and 0.3 s of kick took 90 s.

`kitgrade.py` writes `data/kit.json`, level-matched excerpts, spectrograms and the two grooves
(the same two bars on the real kit and on FrankenSim's). Its ruler is real-vs-real: another stroke
of the same real drum at the same dynamic, and the same drum one dynamic away. Octave-band balance
is floored at -50 dB, because the real overheads carry room noise and cymbal bleed where a render
is silent. `diagnose.py` takes `$KIT $SALAMANDER/OH` as optional fourth and fifth arguments to run
the cymbal-shimmer experiment on the kit's crash.

## Larger live soundboards

The live piano embeds a 32-coordinate board; `public/live/piano/boards/` holds 48/64/96/128-coordinate
boards (not committed, ~17 MB) that the browser fetches only when its timing probe says it can afford
them (`boards.js`). The score lanes (Strandbar, the player) choose board and mode cap together on a
sustained probe (`chooseQuality`); the keyboard section chooses the largest board whose 12-key probe runs
1.4x faster than real time (`chooseBoard`) and leaves the mode budget to its load governor. For the lanes,
a Chromium desktop (8+ logical cores, `userAgentData.mobile` false, no data saving) starts fetching board 96,
their usual pick there, while the embedded board is timed, so the fetch is off the start-up path. Each board the probe times is built once: the Worker keeps a copy
(`PianoEngine.copy()`) from before the timing and plays the chosen board's copy, which renders the same
samples as a new engine (building board 96 takes ~0.9 s, a copy under a millisecond); an engine that cannot
copy (a wasm from an older build) builds the chosen board again. On a device too slow for even the
embedded board at the 700-mode start floor, the piano falls back to its 16 exact low modes, started at
the mode cap its own timing predicts (not below 700) rather than the full 2,400; the load governor takes
it lower from there. The boards are the offline
binary's POD reductions of the same FE slice of the 16-division board the graded renders use
(`board-rt16.fsb`); generate them natively with the engine's `gen_board` example:

    gen_board public/live/piano/boards/model-d-board-2400-r64.bin 2400 pod,64,32,16,8000
    # r48: pod,48,24,16,8000  r96: pod,96,48,...  r128: pod,128,64,... (the studio render's reduction)
    # a fourth argument picks the mesh: 16,0.03 (default) or preset (the unrefined 8 divisions)

Each keeps N/2 exact low modes, plus the leading proper orthogonal decomposition of every key's static
and damped bridge responses at 16 frequencies, log-spaced from the first unprotected mode to 8 kHz.
`examples/board_admittance.rs` measures a reduction against its complete slice
(crates/fs-couple/examples/grand_piano/BOARD_REDUCTION.md has the table). The former boards
(`N,N/2,200,400,800,1600,2300`, largest-residual pivots) sampled nothing above 2.3 kHz. String partials
still drive the slice there, through its mass line, so the pressure those boards radiated above 2 kHz
swung with size. The live A3 mf read -41.7/-37.6/-45.4/-45.0/-35.2 dB at 2-8 kHz for 32...128
coordinates (real -33.2); on the POD boards it reads -34.7/-36.8/-35.4/-36.8/-36.1.

The 16-division board carries 11,772 surface points to the 8-division board's 3,918, so the files are
2.8x larger (the wasm with its embedded 32-coordinate board 1.8 MB gzipped, r128 5.8 MB); building an
engine takes about as long.

Measured on the Strandbar piano part against the studio stem (RMS octave-balance error, 63 Hz-8 kHz,
same 16-division board and 128-coordinate POD reduction, contacts, hammer law and string losses):
32 coordinates 1.52 dB, 48 0.79, 64 1.17, 96 0.50, 128 0.20. These include the silent strings' bank
(`Config::silent_bank_hz`, 600 Hz in the browser: every partial of every string up to the band moves,
driven by the bridge and drained by its damper). Holding the silent strings as rigid chords it was
2.51, 1.92, 2.06, 1.80, 1.74, with 125-250 Hz 2.5-2.8 dB high on every board; on the pivoted boards
5.99, 3.54, 3.03, 2.77, 1.76. On an A3 mf strike the bank takes the early decay from -7.1 to -14.9 dB/s
(offline -14.8) and the 2-8 kHz share from 2.4 dB below the offline render to 0.6 dB above it. C6, all
of whose partials lie above the band, keeps 1.8 dB less energy at 2-8 kHz and a 1.4 dB/s slower early
decay than with every silent partial to 2.4 kHz moving; a 1200 Hz band (1,473 modes) brings it within
1.0 dB and 0 dB/s, but costs up to 25 % more per sample (V8, Strandbar; 14 % on 128 coordinates),
and it moves the 128-coordinate error only from 0.20 to 0.15 dB. In V8 the 600 Hz bank costs
1.31-1.36x per sample on the sustained probe on every board (its gather and kick per course are
board-length); `boards.js`'s cost model includes it.

The browser plays RT-0425's published laws (`wasm.rs`): its string losses, its hammer stiffness K_H on
each unison string and its rate loss R_H in place of the estimated crush and Prony terms
(`Config::rt0425_*`, the reference's `--rt0425-*` flags). With nothing but one key on the bridge,
`examples/validate.rs 800 256 1 --rt0425-string-damping --rt0425-hammer-stiffness --rt0425-hammer-dissipation`
agrees with the reference engine to 4e-12 in bridge velocity (8 keys, 3 speeds, 1 s); run it with and
without the flags after changing the reference's felt, contact or string-loss code, whose fast paths
the live kernel re-implements.

### Live piano strings from canonical material packs

In the live piano's start panel, expand **Choose string materials**, enable the custom construction,
select the original `NormalizedMaterialCardPack` files and the concentric geometry CSV, then enter the
exact temperature axis, temperature in kelvin, solid-viscosity property, cyclic-frequency axis and
query frequency in hertz. The displayed file order defines `material_index`, starting at zero. All
inputs are required for this path and are fixed when Start constructs the engine. The keyboard uses
the embedded 88-key scale, so its geometry must cover every MIDI key from 21 through 108.

The construction file begins with these four exact-form lines: the format, the tension policy, the
sourced tensile-limit property every core is checked against, and the column names. Each subsequent row
is one layer of one key. `core` rows give the **outer** radius of a perfectly bonded core layer, from the
centre outward with strictly increasing radii. `wrap` rows follow a key's core rows and give a helical
wrap's wire diameter and pitch in metres, innermost first; `close` selects the tightest pitch at which
neighbouring turns do not interpenetrate. The last field is the pack index. For a separately imported
scale containing only MIDI keys 33 and 69, a wrapped core for key 33 and a plain two-layer core for 69:

```csv
frankensim-piano-materials-v2
tension,target-pitch
strength,ultimate_tensile_strength
key,layer,outer_radius_m,wire_diameter_m,pitch_m,material_index
69,core,0.0004,,,0
69,core,0.0006,,,1
33,core,0.0008,,,0
33,wrap,,0.0012,close,2
```

The core carries the whole tension, EA, EI and bending loss; a wrap adds only its wire mass (no wrap
stiffness, slip or wrap loss is modelled) and sets the finished radius for gas drag. `tension` is
`target-pitch` (re-solve each tension so the new section keeps the installed fundamental),
`installed-tension` (keep the tension; the pitch follows the material), or `yield-fraction,<x>` with
`0 < x <= 1` (tension is that fraction of the core's capacity under the declared limit). At a fixed pitch
the stress `4 L^2 f^2 rho` does not depend on the wire diameter, so a dense, weak metal cannot reach a
concert scale: if any string exceeds its sourced limit, construction refuses with every overstressed key,
its required stress, its limit and the highest pitch its section reaches at that length (`fMaxHz`).
Appending `,non-physical-override` to the `strength` line lets such strings play and labels them so; real
strings would break. Drawing work-hardens wire, so a limit claim holds at the wire diameters it states
(a typed `wire_diameter` axis): each core's limit is the claim evaluated at the core's own outer
diameter, from the core's packs or from the one pack named by `,pack=<index>`. A core outside the stated
diameters refuses as missing data unless the line declares `,extrapolate-diameter`, which takes the
nearest stated diameter's value; a claim that states no diameter refuses unless the line declares
`,diameter-free`, and such a limit cannot set `yield-fraction` tension. Each core's receipt prints the
resolved limit, the claim rows it came from and, when the limit is not the claim's value at that
diameter, what it is instead. Every core (or the named strength pack) must supply the declared limit
claim, or construction refuses as missing data; the older `frankensim-piano-materials-v1` format names no
limit and refuses for that reason.

Every course of the selected scale needs exactly one complete section; missing or extra keys refuse.
Names never imply a radius, an alloy or a source coefficient. Up to 128 nonempty packs are admitted,
with 16 MiB combined bytes and 256 KiB of UTF-8 geometry. Axes and the viscosity property must be
explicit, trimmed names of at most 128 UTF-8 bytes; the two axes must differ, and both numeric query
values must be finite and positive. Core packs must supply coherent, dimensioned density, Young's
modulus, nonnegative solid Kelvin–Voigt viscosity and the declared tensile limit at that state, with a
source band covering the simulated string partials; wrap packs need density. Missing, ambiguous or
out-of-domain source data refuse in Rust.

For programmatic playback, `LivePianoOptions` and `PianoLane` share `PianoEngineOptions`. This example
accepts user-selected files and an existing physical scale CSV. It assumes the supplied packs declare
typed axes named `T` and `frequency`, and a viscosity property named `solid_viscosity`; use the exact
names in your own cards. Call it with an audio context resumed from a user gesture:

```ts
import { LivePianoClient, type PianoEngineOptions } from './src/live/pianoClient';

async function startMaterialPiano(
  context: AudioContext, packFiles: File[], geometryFile: File, scaleCsv: string,
) {
  const options: PianoEngineOptions = {
    boardCoordinates: 'auto', substeps: 2, scaleCsv,
    stretchTuning: { concertAHz: 440 },
    gas: 'mixture 293.15 101325 mole 2 N2 0.8 O2 0.2',
    stringMaterials: {
      packs: await Promise.all(packFiles.map((file) => file.arrayBuffer())),
      geometry: await geometryFile.text(),
      temperatureAxis: 'T', temperatureK: 293.15,
      viscosityProperty: 'solid_viscosity', frequencyAxis: 'frequency', frequencyHz: 1000,
    },
  };
  return LivePianoClient.create(context, options);
  // The same options can be passed as new PianoLane(mix, { options }).
}
```

This models **perfectly bonded solid/concentric cores, optionally under mass-only helical wraps, at a
fixed material state**. Original source integrals determine line density, axial rigidity, flexural
rigidity and bending viscosity in the live strings, silent-string bank, bridge and hammer coupling. The
source path replaces the legacy string loss law; it adds no unexplained material damping. Length, unison
and duplex geometry remain the chosen scale's data; tension follows the declared policy, and optional
stretch tuning is applied to the material-derived partials and then checked against the tensile limit
again. Wrap stiffness, interlayer slip and temperature evolution need additional physical owners.
The same atmosphere controls the Rayleigh receiver's pressure and propagation delay and the retained
source-bound string modes' cylinder drag, using their explicit finished radii. This frozen modal drag
does not include fluid added mass, viscous memory or loading of the unresolved residual chord.
Effective scales without material geometry do not acquire an inferred radius or gas drag.
`piano.info.stringMaterials` reports the actual source state identities and section coefficients
installed in playback.

The same source packet reaches timing probes, copied engines, larger soundboards and 1×/2×/4× mechanical
clocks (up to 192 kHz). With no `stringMaterials`, the existing constructor routes remain available.
A requested material construction requires WASM rebuilt with `withStringMaterials`; an older build
or a material-owner refusal is reported to the caller. Omitting `gas` on the material path selects the
explicit standard declaration `293.15 101325 0`.

Arbitrary source-backed gases use the same canonical species/NASA-9/transport triples as the live kit,
including with `stringMaterials`. Supply both fields in `PianoEngineOptions`:

```ts
const atmosphere = {
  gas: 'thermochemical 300 100000 mole 1 custom-alpha 1',
  gasSources: {
    components: [{ speciesPack, nasa9Pack, transportPack }], // original Uint8Array bytes
    waterSpecies: null,
    phase: 'frozen' as const,
  },
};
```

Species IDs and gas properties come from those sources. The gas is frozen at preparation; source
models survive copies and `piano.info.gasSourcePacks` reports their identities. Rebuild piano WASM
for `withGasSources` and `withStringMaterialsAndGasSources`. Source limits and the alternative
`liquid-water` policy (with an explicit water ID) are the same as the kit below. Invalid sources and
older WASM refuse instead of substituting air.

## The live kit

`engine-kit/build-wasm.sh` compiles `fs_vkplate::live` (wasm-bindgen) into `public/live/kit/`. Each
piece gets its own Worker (`public/live/kit/worker.js`), which times a stroke at grid scales 2.5-8
(the snare's 6-8: its 5x and 4x sound no closer to the 1x grid than 6x; the hats' also 4.5x: on 5x
they ring a tonal line near 558 Hz, which 4.5x weakens or drops; the kick's 6-8, the rack tom's 5-6
and the floor tom's 4-6, the grids past which they sound no closer to the 1x grid, so their CPU
goes to the pieces that do) and keeps the finest that
runs 1.7x faster than real time on the device. It times from the coarse
end and stops at the first grid that is too slow, so the costly fine grids are never built. It builds
each grid once and keeps a copy from before the timing stroke, and the chosen grid's copy plays: a
snare's or hi-hat's build is almost all its settle (wire preload, pedal; 0.05-0.9 s), a copy takes
under a millisecond and plays the same samples. An engine that cannot copy (a wasm from an older
build) gets the chosen grid built again, which plays the same samples. Each
grid is timed over the first 112 ms of a stroke, the snare's over 450 ms: its wires touch the head
from ~85 ms on and its steps then cost ~1.6x as much, so the short test promised ~1.5x the speed it
kept under rapid tapping, and the snare dropped audio. The player (`src/player/`) times its snare
the same way, at 6x or else 8x, against 1.5 (its lanes render a second ahead). A silent piece is not stepped. Strokes land where the head or cymbal is touched; the Worker sends
the surface velocity at ~45 Hz and, after each stroke, the first 8 ms every 0.25 ms, which the
page replays 60x slower.
Levels: one gain per piece, set on the 5x grid: its 50 ms loudness matched to the Salamander
overheads, capped so a standard stroke peaks at 0.55 (the live cymbals' attacks are too sharp to
match their loudness). A stroke's loudness moves with the grid, by piece and not monotonically (the
crash is 5 dB louder at 4x than at 5x, the closed hi-hat 5 dB quieter; the drums move by under
2 dB, the snare by up to 2.9 dB at 2x), so `build-wasm.sh` also writes `engine_levels.json`: each
piece's standard stroke per grid 2-10x, its loudness (largest 400 ms RMS, so a cymbal's wash counts)
and peak, measured with that wasm by `tools/kit_levels.mjs` and stamped with its SHA-256. The
Worker scales the page's gain so the stroke on the grid it picked is as loud as on 5x
(`public/live/kit/levels.js`): a cut always, a boost only while the stroke still peaks at or under
0.55. So the closed hi-hat stays up to 3 dB quiet at 3.5-4x, and the hats and cymbals at 6x and
coarser, whose strokes already peak over 0.55, are not raised (up to 4.3 dB quiet). Beside any
other wasm the table is ignored and the gain is used as given. The same applies wherever a kit Worker plays
(Strandbar's gains are also measured on 5x).

### Custom gas sources in the live kit

Open **Atmosphere → Import gas material sources** to supply any gas supported by
your canonical source data. Each component needs its original
`NormalizedSpeciesPack`, NASA-9 `NormalizedModelPack`, and pure-species transport
`NormalizedMaterialCardPack`. The three files explicitly associate transport
with that species. The species identity and NASA-9 pack ID must agree; names do
not select substitute coefficients. Transport claims must resolve unambiguously
under `SingleClaimOnly` at the declared state.

Enter the shared record with kelvin, pascals and dimensionless fractions:

```text
thermochemical 300 100000 mole 2 custom-alpha 0.7 custom-beta 0.3
```

The example IDs are placeholders for IDs in your packs, not bundled gases.
Declare the water species ID or `none`, and choose frozen composition or the
existing liquid-water saturation check. Mass fractions are also supported.
The original source domains determine admitted temperatures and pressures;
the five-species convenience model's range is not applied to imported sources.
At most 128 source triples, 16 MiB per file, 64 MiB total and an 8 KiB record
are accepted. Applying changes reconstructs every piece from rest.

Programmatic lanes use the same optional `GasSourceOptions`:

```ts
import { KitLane, type GasSourceOptions } from './src/live/mix';

const gasSources: GasSourceOptions = {
  phase: 'frozen',
  waterSpecies: null,
  components: [{ speciesPack, nasa9Pack, transportPack }], // Uint8Array bytes
};
const lane = new KitLane(mix, 'tom-12', {
  gain: 0.03,
  gasRecord: 'thermochemical 300 100000 mole 1 custom-alpha 1',
  gasSources,
  // Optional `materials` also binds the surfaces, striker and bending loss.
});
```

The original `fs-material` thermochemical evaluator supplies one coherent gas
snapshot to cavity stiffness, added mass, radiation loading, ports and microphone
propagation. `Live::with_resolved_atmosphere` retains the complete `GasRecord`,
including source models, transport claims and canonical pack receipts; clones
and timing retries retain it. `info.gas.sourceSpecies` reports the number of
retained species source triples. The shared `gas_species::input::packs` adapter
is file-independent and available to other native consumers too.

The kit remains a fixed-composition, stationary-medium model with lumped
adiabatic cavities and its existing loading/radiation approximations. This adds
neither chemical reactions nor evolving gas temperature. Imported solid sources
and gas sources compose, including sourced bending viscosity. Source-dependent
playback keeps the declared pressure gain; the standard-material grid level
table is not applied. Rebuild through `engine-kit/build-wasm.sh` to expose the
new source constructors; an older bundle or refused source cannot fall back to
ordinary air or discard a selected solid-material model.

## Strandbar

The studio loop is three native renders that start 0.15 s before the first downbeat:
`grand_piano --performance strandbar-piano.csv` (Steinway D preset on `--board-geometry board-rt16.fsb`,
`--board-band-hz 2400 --board-reduction pod,128,64,16,8000` as the live piano's largest board, one substep,
`--rt0425-hammer-stiffness --rt0425-hammer-dissipation --rt0425-string-damping`, as the live piano plays) and
fs-vkplate `perform` for the kick and the closed hi-hat. Long fs-vkplate
parts can be rendered as bar-sized simulations, each with a pre-roll of real strokes, and stitched
(the stitched kick matched a continuous render to 0.9 % RMS). The hi-hat's banded Airy solve is
memory-bound: run one or two such processes per L3 cache (`taskset`, `-C target-cpu=native`).

    uv run --with numpy --with soundfile python tools/strandbar.py ARRANGEMENT.json STEMS_DIR RECORDING.wav public

It balances the three parts against the recording's low/middle/high band energies, folds the
tail onto the loop's start, writes `audio/strandbar-studio.flac` and `data/strandbar.json` (the
same score for the live piano and kit lanes), and gives each live part the RMS of its studio part
by running the live engines' WebAssembly in Node (`tools/strandbar_live_rms.mjs`, board 32 and the 5x
grid). One master level serves both versions: the studio loop peaks at -1 dBFS at most and the live
parts' sum at -2 dBFS, so the boards and grids a device picks stay under full scale.

On the page the live lanes render two seconds ahead in 512-frame messages. Each drum keeps the finest
grid whose timing reaches 1.4x real time (`src/strandbar.ts`): the kick from 2-8x, the hi-hat from
2.5-4.5x (4.5x drops the 558 Hz line 5x rings), 6x and then 8x; 5x is skipped, where a closed stroke
keeps a 250 Hz component ringing ~4 s even with the rims' friction. Strandbar does not restart
a lane on a coarser grid while it plays, so the lists end at 8x for slow devices.

## Player piano

`src/player/` plays any MIDI file on the live piano and kit: a paper roll over a tracker bar, the
88 keys, the soundboard's modal field in 3D and each drum piece's surface velocity in 3D.
`src/midi/` reads the file: `midi-salvage.ts` is copied from the owner's jazz chord editor
(`jazz_chord_progression_editor_html`), `smf.ts` reads the result under that editor's M0 laws plus
the controllers a player needs, `chords.ts` names chords with its reverse-T1 templates, and
`arrange.ts` maps the file (pitched channels to the piano, GM channel 10 to kit pieces, velocity to
hammer or stick speed, CC64 to the sustain pedal, CC66 to sostenuto and CC67 to una corda, the
grand's own pedals) and counts every fold, merge, thinning and skipped drum. Seeking or resuming
restores the dampers sostenuto had caught with silent key presses before re-striking held keys.
`orchestrate.ts` turns a General MIDI arrangement into one pianist's part first: the file's own
channel volume and expression (CC7, CC11) shape every strike, parts that double or echo another
part are left out, chiptune ornaments are thinned, parts written at one flat velocity get dynamics
by role (before their channel volume and expression), and the page lists what happened to every
part. The shelf's second pack (`public/midi/game-pack*`) holds video-game files of that
kind, supplied by the owner.

The shelf's rolls are that editor's chord-chart library played through its own bass and comping
generators, with drum parts written per groove (from the frankensim root; the output in
`public/midi/` is committed):

    bun apps/instrument-lab/tools/standards/build-standards.ts ../jazz_chord_progression_editor_html apps/instrument-lab/public
    uv run --with mido python apps/instrument-lab/tools/standards/verify_standards.py apps/instrument-lab/public

The lanes take a whole score at once (`{type:'score', events}` with absolute frames, sorted once in
the Worker and fed to the queue only as far as the look-ahead), and a stop or a jump flushes what
they rendered ahead (`mix.flush({epoch})` + `lane.resync(epoch)`; audio, snapshots and fields carry
the epoch they were computed in). The piano lane sends `engine.snapshot()` per 1/30 s of audio
when asked (`snapshotHz`), and a kit lane its fields per `fieldFrames` of audio.
