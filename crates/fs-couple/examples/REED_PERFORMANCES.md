# Authored reed and duct performances

`fs_couple::render::schedule::reed::ReedPerformance::from_bytes` loads a finite
physical pressure phrase, not a synthesizer preset. It implements the shared
`PressureRenderer` interface and can be placed directly inside a
`DecimatedRenderer` or a heterogeneous `PressureEnsemble`. Retain that complete
object when resuming; extracting/reconstructing a bare voice would lose the
control, filter and finite-horizon contract.

## Source format

`reed-duct.performance` is a complete example. Header records are in this order:

```
frankensim-reed-performance-v1
audio 48000 4801 20000
ambient 293.15 101325 0
reed 0.0004 0.013 6000 0 0 0.35
listener_m 1
termination unflanged
segments 1
cylinder 0.0022 0.5
compile_limits 100000 1024
schedule
```

Follow the `schedule` line immediately with the exact canonical bytes of one
`GestureSchedule` blowing-pressure track. The example deliberately uses a
700 Hz control clock against 48 kHz mechanics, interrupted pressure ramps,
release, a second attack and final release. Existing pressure sampling is held
between control ticks; tick `k` applies at `ceil(k * mechanics_rate / control_rate)`.
Neither callback boundaries nor release reset the vibrating reed/bore state.
The track owns the complete mouth-pressure history. There is no additional
fixture attack envelope layered over it.

`audio` contains solver rate in Hz, source samples, and positive PCM full scale
in Pa. Rates are positive and at most 192 kHz; duration is positive and bounded
by both 600 seconds and 28,800,000 source samples. A separate observation clock
requires explicit integer decimation; it never changes these mechanics.
`ambient` is temperature in K, pressure in Pa and relative humidity in [0,1].
The shared moist-air owner's further domain limits apply.

`reed` contains rest opening [m], width [m], closing pressure [Pa], mass [kg],
stiffness [N/m] and damping ratio. Opening, width and closing pressure must be
positive. Mass, stiffness and damping must be nonnegative. Nonzero mass invokes
the existing implicit massive-reed/lay-contact dynamics; zero mass retains the
quasistatic aperture. Zero stiffness retains the existing derived-stiffness
convention. These are authored primitive parameters, **not** an identified reed
species, laminate or temperature-dependent solid material card.

An optional `slit_length LENGTH_M` record immediately after `reed` selects the
shared material-gas channel law. For example:

```
reed 0.0004 0.013 6000 0.00001 100 0.35
slit_length 0.001
```

The length is the physical streamwise channel length in metres, supplied from
geometry; the example value is illustrative, not a measured clarinet reed.
The same `ViscousSlit` owner used by material plate valves then solves
`dp = 12 mu L U/(w h^3) + rho U|U|/(2 w^2 h^2)` at every trial opening and
pressure. It includes the moving reed's swept flow and lay contact in the same
sample. Both massless and massive reeds support this selection. The exact
regional or axial-field **inlet** gas supplies `mu` and `rho`, so composition,
temperature and static pressure changes propagate into the excitation as well
as the bore. The ambient default, first cell midpoint and exterior receiver gas
are never substituted for that inlet state. The CLI sidecar records the length
and actual inlet properties. Omitting the record preserves the inviscid path.

This is quasisteady, fully developed parallel-plate resistance plus separated
jet loss. It regularizes the small-gap flow, including nonzero flow below the
legacy Bernoulli pressure cutoff, but it does not supply unsteady channel
inertia, entrance development, turbulence, or measured-instrument validation.
The optional channel uses the shared strict root solver; fast-Newton mode
explicitly records a fallback rather than using the inviscid Jacobian.

`reed-viscous-slit.performance` is a complete runnable demonstration with a
massive reed, the 1 mm illustrative channel and the same pressure phrase as
`reed-duct.performance`:

```bash
cargo run -p fs-couple --bin music_render -- wind \
  crates/fs-couple/examples/reed-viscous-slit.performance reed-viscous-slit.wav
```

It retains the diagnostic `listener_m` observation and authored mechanics; it
is not a reconstructed or validated clarinet. Change `ambient`, add regional
segment gases, or edit `slit_length` to change the physical model through data.
The WAV and its provenance JSON must name new output paths.

`listener_m` is the distance used by the existing compact-jet observation.
`termination` is one of `closed`, `ideal-open`, `unflanged`, `flanged`.
Each of the declared 1..256 segments is one of:

```
cylinder RADIUS_M LENGTH_M
cone INLET_RADIUS_M OUTLET_RADIUS_M LENGTH_M
hole HOLE_RADIUS_M CHIMNEY_HEIGHT_M BORE_RADIUS_M OPEN_FRACTION
```

Radii and lengths are finite positive SI values. A hole must be narrower than
its stated main bore; opening is in [0,1], never silently clamped. Fractions
between zero and one use the existing open/closed admittance interpolation.
The first and last records must be axial segments. The same existing Bessel/TMM
reflectance owner handles the entire ordered duct. Its FIR construction,
passivity and low-`ka` radiation admission can still refuse a source. No cutoff,
loss model or fit tolerance is relaxed to accept supplied geometry.

`compile_limits` contains the allowed source visits (at most 16,777,216) and
stored controls (at most 262,144). Input is capped at 1 MiB, with at most 16,384
pressure events and 262,144 sampled control ticks. Counts are bounded before
count-driven decoding/compilation. Headers reject missing/extra/trailing fields;
noncanonical schedules and unsupported tracks refuse. Every command must start
by the last sampled control tick. A ramp may extend beyond the finite window;
its unobserved continuation is not fabricated.

## Observation and scope

`listener_m` selects the current `ReedBoreVoice` **bore-pressure plus compact-jet
proxy** with an empty plate bank. The physical bell/side-hole observation below
is a separate supported selection. Neither observation supplies a complete
instrument material reconstruction or a moving fingering mechanism.
The existing constructor's 5 Pa traveling-wave seed is
retained, not silently replaced with a new initialization law. Passive bore
filtering and structural/contact dynamics stay in their original owners.

### Independent room gas with open tone holes

Replace the `listener_m` record with the following records to use the existing
baffled microphone and an independently authored room atmosphere:

```
baffled_outlet 0 0 0.5 4 16 1000
exterior mixture 303.15 100000 mole 3 N2 0.7 O2 0.2 CO2 0.1
observer_law stokes-kirchhoff
termination flanged
```

The receiver coordinates are metres; the last three numbers select radial
rings, angular points and the observation band in Hz. The `exterior` record
uses the same gas-material resolver as the bore and must have the bore's
static pressure. Start from `reed-baffled-mixture.performance`, whose common
pressure is 100000 Pa, and edit only the exterior data to compare atmospheres.
`reed-two-medium-tone-hole.performance` supplies a complete example with a
partly open chimney, two internal gas regions and a third exterior state:

```bash
cargo run -p fs-couple --bin music_render -- wind \
  crates/fs-couple/examples/reed-two-medium-tone-hole.performance two-medium.wav
```

Open and partially vented `hole` records now support this independent exterior.
Each chimney keeps its own bore-region gas for propagation, internal junction
corrections, wall coupling and viscothermal loss. At its opening the shared
flanged radiation owner matches pressure and volume flow to the room gas,
including the internal higher-mode impedance. The same resulting hole load
changes reed feedback and the radiating volume flow heard by the microphone.
Closed holes remain sealed. Bell and hole propagation to the receiver use the
same declared room gas. No instrument-specific retuning or code change is
required for a material substitution.

This is linear acoustic matching across prescribed stationary gases, not a
simulation of breath mixing or airflow. Existing restrictions on mutual
coupling across different internal hole media remain. The piecewise regional
renderer supplies bell and hole pressure; the continuous gas-field observer
currently supplies bell pressure only. The existing flanged-aperture and
one-dimensional bore approximations remain, without a measured-instrument
accuracy claim. Native regressions check open and vented holes through reed
motion, physical microphone pressure, block-independent WAVs and sidecars.

Superposition with another source requires a physically compatible observation
and common time origin supplied by the caller; matching scalar units and sample
clocks alone cannot prove that. Keep this limitation when using ensembles.

The finite source refuses a request extending past its remaining physical
window before applying any pressure commands. Callback-shape admission also
precedes controls. Cancellation at the common WAV render boundary preserves
all source/filter history for exact continuation with the same objects. A late
physical failure retains the renderer's poison-on-failure semantics.

`tests/reed_performance.rs` compares the file adapter against direct physical
stepping with independent per-sample control application. It covers interrupted
ramps, release, callback partitions, massive mechanics, static hole/cone geometry,
a changed gas state, finite-window refusal and decimated cancellation/resume.
These are implementation regressions, not measured-instrument validation or a
real-time throughput claim.

## WAV and ensemble commands

Run the supplied phrase on its declared 48 kHz clock:

```bash
cargo run --release -p fs-couple --bin music_render -- \
  wind crates/fs-couple/examples/reed-duct.performance /tmp/reed-phrase.wav \
  --block 37
```

For a supplied higher-rate performance, add `--decimate`. The high-rate source
must already satisfy the existing full-band radiation/realization domain. For
example, doubling the mechanics rate of the supplied 2.2 mm-radius unflanged
bore without changing the physical source violates the low-`ka` load limit.
The regression's 96 kHz source explicitly supplies a 1.1 mm-radius bore instead;
no geometry is changed automatically. Noninteger rate ratios and incomplete
output intervals refuse. The declared window includes filter startup delay;
no truncated tail is normalized, padded or flushed after that window.

The same source can enter the existing file-driven ensemble:

```bash
cargo run --release -p fs-couple --bin music_render -- \
  ensemble /tmp/reed-and-plate.wav --full-scale-pa 20000 \
  --reed crates/fs-couple/examples/reed-duct.performance \
  --plate crates/fs-couple/examples/plate-mesh.performance --block 37
```

This command is a scalar-pressure composition example, not a claim that the
reed bore and plate exterior are one calibrated microphone. Source observation
compatibility remains the caller's responsibility. `--reed`, `--bow`, `--modal`
and `--plate` may be combined in explicit summation order. Every part must have
the same physical duration; the ensemble's scale is applied only after mixing.
The source's own PCM scale is recorded but is never used as a hidden gain.
The reed observation limitation is included in the part's provenance.

`wind` preserves the source scale; `ensemble` requires a separate explicit scale.
Both use the existing incremental PCM16 writer, clip counts and content hash.
Input/model/clock refusal precedes output creation, and existing WAV/sidecar
files are never overwritten. The old `music_render reed` fixture is unchanged.
