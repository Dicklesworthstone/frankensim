# One exterior gas for instrument reaction and pressure

The one-way `DuctPressure` / `TwoSidedPressure::with_radiation` selection leaves
mechanical terminal loads alone. That remains useful for prescribed-flow studies.
It is not sufficient when an exterior-material change should also change the
load seen by the instrument. `bind_baffled_environment` makes that connection
explicit before either the source or its receiver advances.

The function uses the actual selected terminal radii, the exterior `GasState`,
the declared mechanical clock and the observer's physical frequency band to
construct the existing `BaffledRadiationLoad`. All ports are admitted by the
original `DuctPressure` geometry/gas/clock checks. Only selected memoryless
endpoints may be replaced; an existing reactive terminal is not silently
retuned. Every candidate passes before any graph changes. Unselected loads,
regional internal gases, geometry and propagation impedances stay unchanged.

## Material reeds

`TwoSidedPressure::from_baffled_environment` accepts an existing unadvanced
single- or paired-blade valve, its two original duct sources, optional original
wall-material assignments, original limits, optional fractional policies, a
`DuctRadiationSpec`, an output callback bound and a cancellation gate. It returns
`(TwoSidedPressure, Vec<BaffledOutletBinding>)`. The performer directly uses the
shared canonical pressure scheduler and PCM interfaces. Boundary records expose
the exact installed load, previous memoryless condition, selected port and
original internal gas.

Network 0 is upstream and 1 downstream. For example, changing only
`radiation.exterior` before this construction changes both the radiating load
and the microphone's propagation state. It does not overwrite either bore gas,
retune blade mass/stiffness, reset material history or normalize the sound.
Changing only that field on the *one-way* constructor changes observation alone.
The distinction is deliberate, not two conflicting definitions of atmosphere.

The constructor delegates material walls and selected fractional/viscothermal
propagation to `TwoSidedAperture::with_material_walls`. A selected fractional
band must cover the declared exterior band. The mechanical clock must match the
original valve exactly; output may be an integer factor slower through the
existing differentiating receiver. Source gestures remain on mechanical samples.

## Jet-driven sources

`FlueVoice::from_baffled_environment` takes the existing jet card, receptivity and
seed parameters, geometry, `FlueNetworkSpec`, optional material walls, optional
fractional policy, radiation declaration, stream identity and gate. Its source
clock comes from `radiation.reconstruction.mechanics_rate_hz()`, not a second
rate parameter. It returns `(FlueVoice, BaffledEnvironment)`.

To construct the matching receiver, borrow the returned voice's actual graph and
regional gas mapping as a `DuctView` and pass `environment.radiation().clone()`
to `DuctPressure::new`. For every accepted source sample, pass all original
terminal volume flows to `push_frame`. `None` is an incomplete multirate packet,
not silence. The shared observer then supplies physical pressure to the existing
PCM path. `FlueNetworkSpec::with_baffled_environment` exposes the same cold binding
for compositions that already own their source factory; those callers must use
the declared mechanical clock when constructing the voice.

## Physical scope and verification

This uses the existing **independent compact baffled-piston loads**. Each load
retains its ka <= 0.5 domain and fixed 5% analogue/discrete response checks. The
observer does not enlarge that reaction-model band. Multiple openings interfere
in the observed pressure, but their mutual radiation reaction is not modeled.
The geometry remains stationary, coplanar and infinitely baffled. Pressure-
matched differing gas states are a small-signal interface, not mean-flow gas
mixing, thermal equilibration or live atmospheric evolution.

`cargo test -p fs-couple --test baffled_environment` covers boundary coefficients,
local/exterior identity, transactional refusals, actual reed/jet feedback,
multirate owner equivalence, callback partitions and fractional selection.
The nine tests were authored but not executed in the implementation container,
which has no Rust toolchain. No measured instrument-fidelity claim is made.
Independent midpoint references executed 3,072 loaded reed/network steps and
360 frequency cases, with maximum combined energy residual 1.14e-21 J. Those
references are separate equations, not execution of this Rust implementation.

## Plate-derived comparison example

The existing DKT-reduced two-blade example now accepts an optional exterior
absolute temperature. All original three-argument commands still use 293.15 K.
Both internal ducts stay at 293.15 K; the optional fourth argument changes the
external radiation load AND pressure propagation, not the blade geometry or
source-pressure schedule. Existing destinations refuse; audio is not normalized.

```sh
cargo run --release -p fs-couple --example two_sided_reed -- exterior293.wav 0.12 900 293.15
cargo run --release -p fs-couple --example two_sided_reed -- exterior333.wav 0.12 900 333.15
```

These are coupling experiments with authored plate/contact data, not calibrated
oboe presets. The updated example's WAV output has not been executed or auditioned
in the implementation container. Its source still reports derived blade mechanics,
represented duct lengths, clipping and final stored energy.
