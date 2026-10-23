# Finite-body piano acoustics from supplied geometry

`piano_exterior response` connects the physical piano's existing modal motion
to the existing Helmholtz boundary-element solver. Unlike the ordinary Rayleigh
observer, this calculation has no infinite baffle: supply both sides and edges
of the soundboard and any rigid cabinet/lid surfaces that should scatter sound.
A lid is part of the same boundary solve, not a second sound source or an EQ.
`piano_exterior render` then connects those transfers to the actual played piano.

```sh
cargo run --release -p fs-couple --example piano_exterior -- \
  response settled.fss strings.csv acoustic-body.obj acoustic.fspe response.csv
cargo run --release -p fs-couple --example piano_exterior -- \
  render settled.fss strings.csv acoustic-body.obj acoustic.fspe piano.wav 6 score.mid
```

`steinway-d` may replace the scale CSV to select all 88 raw source Model D
courses. An explicit CSV preserves its supplied tensions; use the SAME scale
that produced a settled board. Flat and crowned/preloaded structural files are
supported. The existing source felt, shank and reciprocal string bank prepare
the loaded basis, including silent strings. No source course is discarded.

## Acoustic specification

Every singleton row below is mandatory; numbers illustrate the schema, not a
measured Steinway or a promise this bandwidth will pass on an arbitrary mesh.
Object/group labels must match the supplied OBJ exactly.

```text
frankensim-piano-exterior-si-v1
source,estimated,Explicit acoustic idealization; replace with actual attribution
obj-scale-m,1
obj-origin,0,0,0
max-skin-offset-m,0.015
medium,1.204,343
band-hz,40,800,41
board-band-hz,400
fit-order,8
min-panels-per-wavelength,6
full-scale-pa,2
receiver-m,0.6,1.0,3.0
moving,soundboard_skin
rigid,cabinet
rigid,lid
```

The exterior medium is exactly one of two rows. `ambient,<gas record>` is the
shared gas declaration, resolved through the gas owner. `medium,RHO,C` is a bare
medium with no thermochemistry, and every report marks it "bare medium, not a
gas record". The grand_piano microphone, Rayleigh loading and imported-board
renders use the committed reference record `examples/moist-air-20c-65rh.gas`
(20 C, 101325 Pa, 65 % RH) unless `--gas-record FILE` overrides it. The Model D
soundboard reads the same record for its wood's equilibrium moisture content.
`examples/dry-air-20c.gas` remains available as an explicit record.

`response` and `admittance` accept **1 to 64 `receiver-m` rows**, in the order
used by the output indices. Coincident receivers remain distinct rows. Optional
`receiver-pattern` indices follow the same order; see `DIRECTIONAL_RECEIVERS.md`.
All receivers share the source boundary solve. Near-field evaluation also reuses
each frequency's prepared observation rows across modal inputs. An array does not repeat piano preparation
or change the source impedance. `render` and `render-loaded` accept only one or
two receivers and reject a larger array before structural or acoustic preparation.

Positions and medium values are coherent SI in the structural board frame. Transform OBJ
coordinates by `p_m = obj-scale-m * (p_obj - obj-origin)`; origin is in source
OBJ units. Axes must already agree with the board axes. This does not infer
orientation, centimetres/metres, acoustic impedances or mechanical material
cards from visual MTL properties. No referenced material file is opened.

Every OBJ region must match exactly one `moving`/`rigid` label; unused labels,
unmapped regions and ambiguous assignments refuse. A moving skin gets its
normal velocity from the SAME structural eigenvectors and mass-loading basis
as the string-bridge mechanics. Facet-normal projection must lie in a real
structural triangle and within the declared skin offset, at every acoustic
corner and three quadrature points. Through-thickness rotation contributes
`theta cross arm`. The crowned path retains all Cartesian motion and the
actual equilibrium geometry, not the flattened Rayleigh field. The acoustic
OBJ itself must describe the intended equilibrium skin: the offset allowance
is not permission to invent thickness or silently deform a supplied asset.

## Geometry and numerical admission

Supply outward-oriented closed components. A finite board needs top, bottom
and edge surfaces; a floating rigid lid must also have thickness and closure.
Exact coincident OBJ seam vertices are identified, but no tolerance welding,
hole filling, face flipping or invented cabinet is performed. Duplicate faces,
inconsistent/open edge uses and nonpositive component volume refuse. These
checks are NOT a general self-intersection, nonmanifold-vertex, component-overlap
or fluid-accessibility certificate; those remain input responsibilities. Do not
supply intersecting, nested or inaccessible closed pieces as independent bodies.

The current dense BEM budget is 2048 acoustic triangles and 128 input modes.
Structural motion projection uses a source-triangle bounding tree and keeps a
3,000,000-work cap on tree visits, triangle bounds and exact facet checks;
the old acoustic-panel × structural-triangle estimate no longer rejects a
spatially sparse, separately resolved pair. For acoustic meshes of at least
1024 panels, independent frequencies use at most four scoped workers. The
768 MiB budget covers their estimated dense matrices plus one shared static
geometry table. Sweeps with at least two plain-CBIE frequencies prepare the
exact frequency-independent triangle integrals once (16 bytes per panel pair).
At 2048 panels this admits three workers and a 64 MiB shared table. A single
frequency or a Burton–Miller-only grid needs no table. Formulation selection,
dynamic triangle integration and the modal factorization still run at each
frequency; there is no interpolation between frequencies. The same preparation
is used by harmonic admittance and loaded-playback fitting, which also reuse
their modal velocity fields and receiver admission. Results return in the
original frequency order. Each requested frequency must
meet the declared minimum panels/wavelength (at least six); negative radiation
power beyond roundoff and nonfinite diagnostics refuse. The reported condition
number lower bound is not a condition-number certificate. Refine the input
mesh or choose a narrower explicitly declared band when resolution is refused;
no modes or frequencies are silently dropped.

The default centroid receivers must conservatively lie outside the body's
enclosing sphere. Explicit `receiver-evaluation,near-field` uses actual
point-to-triangle clearance and triangle-integrated observation instead; see
`NEAR_FIELD_RECEIVERS.md`. Both require at least two 48 kHz propagation samples
and at most 0.5 seconds of guaranteed travel time. Receivers below the soundboard
are legal. Near-field arrays retain the existing **2,000,000 kernel-evaluation
budget across all receivers and panels at each frequency**; the 64-receiver
count is a ceiling, not a promise every mesh/array will fit that work budget.
The BEM evaluates finite-distance pressure, including spreading and propagation;
no extra 1/r or gain is applied.

The CSV contains receiver-major/input transfer values at each frequency:
pressure [Pa] per unit generalized modal acceleration [m sqrt(kg)/s^2], with
`exp(-i omega t)` phasors. The boundary input is `v_n = i shape / omega` for
unit acceleration. All receiver responses share one boundary factorization
and modal batch at a frequency.

## Render through the supplied geometry

The optional MIDI filename selects the existing channel-1 importer, default
velocity 127 -> 4.5 m/s, and the existing sustain controls. Without MIDI the
example strikes key 69 at 2 m/s. MIDI supplies gestures, never material data,
string frequencies or missing courses. Duration must be 0.05..60 seconds.

The render wrapper defaults to frequency-based string retention through
21.6 kHz (at most 512 partials per segment), four mechanical substeps per 48 kHz
sample, source wool/Prony felt and the existing shank image. `--modes 24` retains
the smaller preview budget; any count-limited segments are reported.
The complete board mode set through `board-band-hz` is required. This is a
bounded offline rendering image, not a full-band Steinway fidelity claim.

Before any mechanical step, the existing fs-vfit owner fits one stable proper
transfer per loaded-board mode and receiver. Only even-index frequency values
are supplied as fit samples; odd values check the realized digital transfer.
The adapter removes only the whole-sample part of the guaranteed flight phase,
leaving the fractional phase in the fitted and held-out transfer. It conjugates
the negative-time BEM convention and Tustin-warps the fit frequencies. It restores physical Pa
scaling before rendering. A sampled row's peak magnitude normalizes only its
error metric; it never normalizes the output audio.

Held-out peak-normalized error must be <=0.15 at every odd frequency, and RMS
must be <=0.05 for each mode/receiver. Failed, unstable, improper or nonfinite
fits refuse instead of relaxing these bounds. These are sampled transfer-fit
checks including the residual fractional propagation, NOT whole-band or spatial
convergence certificates. The exact integer delay has unit magnitude. The
decimator still has its own numerical approximation and declared latency:
52 output frames at 4x mechanics (44 anti-alias plus 8 derivative), or 48 at 1x.

Every accepted mechanical substep's loaded-board end-velocity enters the shared
`Decimator::new_differentiating` observer. Its bandlimited derivative reconstructs
SI acceleration before both receivers consume the same output. This removes
the former finite difference's high-frequency attenuation and implicit half-
mechanics-sample phase shift. Each receiver has its own filter histories and an
integer-sample delay line. No second piano, per-channel gain,
control clock, future-sample padding or latency compensation is added. Pressure
is encoded through the existing PCM16 owner with explicit `full-scale-pa` and
counted clipping; there is no peak normalization.

The continuous-acceleration regression drives the actual pressure observer and
receiver runtime at 440, 4800, 15000 and 21600 Hz, at mechanics ratios 1, 3, 4
and 16. Both independently delayed receiver channels must remain below 0.0005
relative RMS error. This checks the numerical observation path; it does not
establish structural, BEM, whole-render or listening accuracy.

The flight-phase treatment avoids the former unchecked two-tap interpolation loss (about 5.1 dB at
15 kHz for a half-sample flight delay). Fractional phase now consumes some of
the requested fit order; insufficient fits refuse under the same error gates.
The runtime regression checks analytic delayed pressure at 4.8, 15 and 20 kHz
for four fractional offsets, including actual state-space and delay playback.

The report prints the acoustic band, resolution, fit error, propagation and
filter delays, peak pressure, clipping, and the unchanged mechanical energy
accounting. An error creates no candidate WAV before successful preparation and
rendering. An OS write failure can still leave a partial newly created file;
outputs are never overwritten and no filesystem transaction is claimed.

## Physical boundaries and assets

This is one-way linear exterior acoustics on fixed geometry about the supplied
structural equilibrium. Radiation pressure does not feed back into mechanics.
The selected rigid cabinet/lid does not flex or absorb; no room, air absorption,
string-direct radiation, nonlinear air or above-band accuracy is claimed. A
hammer attack contains energy outside a narrow fitted band; a stable fit there
does not validate the complete rendered attack. Mesh resolution and the
structural/modal truncation need independent convergence and measurement
comparisons before a realism claim.

A public asset lead is seavenois's CC0 Steinway D274 on BlendSwap:
https://blendswap.com/blend/7279 . The author describes it as simplified; its
render geometry is not automatically an admissible structural/acoustic mesh.
No external model bytes or unverified factory measurements are included by
this feature. A licensed surface still needs explicit acoustic part mapping,
physical structural cards, coherent scale and mesh preparation. Attribution
and measured/published/estimated labels do not themselves verify the data.
