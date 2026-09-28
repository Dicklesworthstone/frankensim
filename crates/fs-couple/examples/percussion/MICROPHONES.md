# Close and directional percussion microphones

`--microphone-spec INPUT.frm` selects one or two fixed physical microphones for
`splash-mic`, `drum-mic`, `drum-stretch-mic`, `drum-modal-mic`, `snare-mic`,
`snare-off-mic`, and `hihat-mic`. Each location is tested against the **actual
retained source triangles**, not an enclosing sphere. A receiver can therefore
sit just above a drumhead, below the resonant head, beside a cymbal, or between
the separated reference hi-hat shells when it is geometrically exterior and
satisfies the declared clearance. Nothing relocates a microphone to make it fit.

The existing `fs-bem::near_field` owner integrates the boundary solution's Green
representation adaptively at each point. Directional receivers additionally use
its analytic particle velocity. The same source-mode boundary solves supply all
microphones and optional radiation loading; changing an observer does not
change head/stick/cymbal motion, contact, gas, or the acoustic force matrix.

## Complete input

UTF-8, at most 8 KiB; blank lines and `#` comments are allowed:

```text
frankensim-microphones-v1
near_field,0.001,0.000001,14,2000000
receiver,0.04,0,0.105
receiver,0.04,0,-0.105
pattern,0,0.5,0,0,-1
pattern,1,0.5,0,0,1
```

The included `close-drum.frm` contains these **declared placements**, not measured
microphone response or instrument calibration. They are 22.45 mm beyond the top
and bottom reference planes of the nominal 6.5-inch-deep drum. Other geometries,
including cymbals or hi-hats, require their own explicit positions and axes.

`near_field` occurs exactly once and gives minimum triangle clearance [m],
relative quadrature tolerance, maximum subdivision depth, and maximum kernel
work **per receiver per frequency**. The existing owner admits tolerance
1e-12..1e-3, depth 0..16 and 80..20,000,000 evaluations. These are quadrature
controls, not bounds on the source model or boundary-solve error. All selected
receivers must pass geometry admission before any boundary source solve or
physical time advance. Exhausted quadrature refuses instead of using centroid
values or weakening the request.

One or two `receiver,x,y,z` rows specify metre coordinates in the instrument's
reference frame, in WAV channel order. Coincident points are allowed for
comparing different patterns. Each optional `pattern` gives a zero-based
receiver index, pressure fraction alpha, and an explicit unit front axis.
Duplicate patterns, out-of-range indices, malformed data and extra records
refuse. Omission selects an ideal omnidirectional pressure microphone.

The front axis points **from the microphone toward its front source**. It is
never inferred from the body centre or silently normalized. Alpha=1 is omni,
alpha=0.5 cardioid, alpha=0 figure-eight; intermediate values are admitted.
For co-located pressure p and particle velocity v, the normalized observation is

```text
output = alpha*p - (1-alpha)*rho*c*(front_axis dot v).
```

The directional output is **Pa-equivalent at unit on-axis plane-wave
sensitivity**, not raw scalar pressure and not microphone voltage. Orientation
and reactive near-field velocity enter the observation itself, not a guessed
source direction or a distance-based bass EQ. There is no diaphragm, housing,
electronics, frequency-dependent polar response, noise or microphone-body
backreaction model. Pressure/velocity are linear exterior fields of the fixed
reference acoustic scene, even when the mechanical source is nonlinear.

## Use with the existing physical performance

```sh
cargo run --release -p fs-couple --example percussion -- \
  drum-mic 4800 20 --analytic-newton \
  --microphone-spec crates/fs-couple/examples/percussion/close-drum.frm
```

This is an invocation, **not a claimed completed render or fidelity result**.
The existing acoustic bandwidth, panel-resolution, work and independent
holdout-fit checks still apply. Use `--radiation-spec` to request appropriate
source-preserving acoustic refinement and a supported band; close microphone
admission does not make an underresolved source mesh adequate.

The file supplies every receiver. Combining it with positional microphone XYZ
or `--microphone-right` refuses rather than picking one silently. CSV mechanics
and far-field `-wav` commands refuse the new option. Without it, the original
finite-point/far-field receiver paths and their admission remain unchanged.

Existing first/second sticks, shaft flexure, force programs, head stretching,
snares, material memory, compliant mutes, and hi-hat gap-gas selection still
construct their same mechanical model. `--radiation-feedback` keeps its existing
restrictions and source addresses. Explicit prescribed vent-flow source columns
remain source columns; a microphone never becomes an extra mechanical mode.

## Timing, amplitude and limits

The nearest triangle distance provides a lower bound on flight time. When it
is at least two output samples, that interval is peeled from the frequency
response and realized by the existing propagation line. At shorter distances,
**no explicit delay is inserted**: the entire physical propagation phase stays
in the fitted transfer. This is not zero flight time, an artificial two-sample
latency, or a moved microphone. Stable finite-band fitting is an approximation,
not an exact wavefront or all-frequency phase certificate.

The original 48 kHz output, sixteen mechanical substeps, step-average
acceleration convention and decimator remain. All receiver filters consume the
same accepted mechanics block. PCM full scale remains the supplied physical Pa
or Pa-equivalent value; there is no channel normalization or second 1/r gain.
The shared pressure fitter canonicalizes polarity as well as amplitude before
identification and restores the signed scale on the resulting transfer.

Exterior admission requires retained closed, consistently outward components.
It checks exact-coordinate edge closure, component-wise exterior winding, and
face/edge/vertex clearance; this is not a global mesh-intersection certificate.
The reference surface is stationary. A microphone admitted in a hi-hat gap is
not certified to stay exterior as the real pair closes, and scattering geometry
is not recomputed during motion. Moving microphones and calibrated capsule
models remain separate capabilities.

Focused checks in the existing percussion example test target:

```sh
cargo test --release -p fs-couple --example percussion acoustics::receivers -- --test-threads=1
```

Tests cover physical close-head geometry, direct Green-field/pattern identities,
short-flight handling, complete input and conflict refusals, and real BEM-to-PCM
stereo with and without source loading against an independently advanced
mechanical trajectory. Native execution and measured comparisons must be
reported separately; passing input or geometry checks is not instrument realism.
