# Prepared percussion execution

## Nonlinear cymbal and stretching-drum execution

`--prepared-nonlinear` selects reusable Gonzalez execution **after the same
physical construction**. It works with `splash`, `drum`, `drum-stretch`, and their
`-wav` and `-mic` forms. Omitting the flag retains the original numerical path.
The flag changes neither geometry nor material parameters, modal bandwidth,
strike position/velocity, felt conditioning, cavity volume or sample clocks.

```sh
cargo run --release -p fs-couple --example percussion -- \
  splash 4096 --prepared-nonlinear > splash-prepared.csv
cargo run --release -p fs-couple --example percussion -- \
  drum-stretch 4096 --prepared-nonlinear --strike-speed-m-s 4 > drum-stretch-prepared.csv
(set -C; cargo run --release -p fs-couple --example percussion -- \
  drum-stretch-mic 48000 20 0.08 0.05 0.35 --prepared-nonlinear \
  --strike-speed-m-s 4 --strike-position-m 0.06 0.01 > drum-stretch-prepared.wav)
```

The existing BEM receiver preparation, decimation, propagation and PCM encoder
are unchanged. Pressure remains an observed physical response, not a fabricated
strike waveform or normalized modal sum. This flag also works with far-field WAV
export. Acoustic baking is still offline and limited to its declared band.

Library callers can use `ImpactSystem::prepare()` and feed the result directly
to `ImpactPressureRenderer`. It retains a current strike and its material history;
`into_reference()` moves back without resetting either. Physical acceptance is
shared with the reference: rejected solves never advance motion, felt state or
the clock. Numerical cancellation is additionally polled inside Newton/Jacobian
work. The read-only reference observations remain available on the prepared host.

Preparation removes solver/host scratch allocation, not the dense Newton cost.
Storage callbacks must also avoid allocation. This is not a native deadline or
allocation certificate, and does not expand the existing geometric/material
validity range. `drum-modal` and `snare` explicitly refuse this flag rather than
silently replacing their distinct linear-body/contact realization.

## Prepared linear-head image

`drum-modal` compiles the existing two-head drum into the retained exact-ZOH
modal/contact owners. This removes the all-coordinate finite-difference Newton
solve for this already-linear body model. `drum` remains the Gonzalez reference;
`splash` stays nonlinear, with its original shell and felt history.

```sh
cargo run --release -p fs-couple --example percussion -- drum-modal 4096 > drum-modal.csv
(set -C; cargo run --release -p fs-couple --example percussion -- \
  drum-modal-mic 48000 20 0.08 0.05 0.35 > drum-modal-mic.wav)
(set -C; cargo run --release -p fs-couple --example percussion -- \
  drum-modal-wav 48000 20 > drum-modal-far.wav)
```

The microphone command takes audio frames, PCM full-scale in pascals, and an
optional complete x/y/z position. The WAV command retains the existing far-field
receiver. There is no peak normalization. The source geometry, film eigensolve,
contact tip, damping, gas volume, BEM boundary, withheld radiation fit, observation
decimator, propagation and PCM encoder are shared with the reference commands.
The output clocks remain 500 kHz for mechanics CSV and 768 kHz mechanics / 48 kHz
pressure for audio. Acoustic preparation is still an offline BEM/vector-fit task.
The 40--1640 Hz fitted acoustic band is unchanged, not expanded by a faster solver.

## What changes numerically

Each isolated linear body uses `ModalAcousticTimeModel`'s exact held-force
transition. `CoupledModalSystem` factors the bilateral coupling matrix once.
`ContactModalSystem` or `MultiContactModalSystem` condenses that whole linear
network and solves only the scalar/joint contact reactions. The same contact
potential, gaps, weights and body basis are retained. A moving free striker is
still a mass with initial kinetic energy, not an authored force pulse.

For the sealed volume, the caller supplies a coordinate reference area A:

```
x = (sum modal_area_i * q_i) / A
k = bulk_modulus * A^2 / cavity_volume
H_air = k*x^2/2 = bulk_modulus*(sum modal_area_i*q_i)^2/(2*cavity_volume)
```

Thus A cancels from physical storage and pressure. It is not an inferred piston
or a radiation correction. The example declares the clear-span disk area for A.

This finite-step port coupling is not the exact exponential of the complete
coupled system and is not bit-identical to Gonzalez. Both require time-refinement
checks against the same physical equations. The new image explicitly refuses
nonlinear shells, felt memory, free-coordinate drag and single ports spanning
more than two bodies; it does not silently erase them or fall back to a cheaper
model. Many independent two-body ports and distributed joint contacts are legal
within the original owners' work limits. A broader shell/air/hardware model needs
its appropriate image, not more relaxation of these admissions.

## Public composition

`render::plate::impact::linear::LinearImpactSystem` accepts the original
`ImpactBody`, `Obstacle` and `VolumeSpring` data plus explicit owner budgets.
It retains initial vibration, reports complete energy and contact residuals,
and supports transactional cancellation and lifetime-budget extension.

The existing `ImpactPressureRenderer<'a, M = ImpactSystem>` now accepts both
images through `ImpactSource`. Its source-specific diagnostic type prevents a
contact residual in newtons from being mislabeled as a Gonzalez state residual.
Both use the unchanged directional radiation, `DecimatedRenderer` and
`render_pressure_pcm16` APIs. No additional renderer or WAV format is introduced.

## Evidence and remaining work

Focused native regressions compare the two mechanical images under refinement,
check volume energy and area-scale invariance, simultaneous contacts, exact
silence without coupling, pause/refusal/budget replay, and the full prepared
pressure/decimation/PCM route. An example regression constructs the actual
existing drumhead meshes and compares onset states and source maps.

These Rust tests and the full examples were **not executed** in the authoring
environment: no Cargo/Rust toolchain was available. Independent Python/SciPy
calculations are separate numerical evidence, not native execution or timing.
At 96/192/384 kHz a synthetic coupled impact's scaled error against a continuous
DOP853 reference decreases from 3.159e-5 to 7.897e-6 to 1.974e-6. The maximum
prepared-image energy residual in that comparison is 1.339e-17 J.

No hard-real-time qualification follows from removing dense Newton. Existing
contact evaluation still allocates, the modal/contact hot path still has work to
optimize, and there is no native callback benchmark. These remain linear head
mechanics with small-signal sealed air, not a complete snare, a measured material
identification, full-band cymbal nonlinearity or radiation-force back-coupling.
