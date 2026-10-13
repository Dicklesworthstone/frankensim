# Exterior pressure feedback on a drum vent

`--coupled-vent-radiation` connects the existing head/mouth BEM boundary to the
same joint mechanical time owner that advances the heads, wire contacts,
interior acoustic modes and neck slug. It implies `--radiation-feedback`.
This is a stationary, linear exterior impedance coupled to the retained
nonlinear mechanics, not a resolved exterior fluid simulation or a calibrated
instrument.

The accompanying `--cavity-neck-internal` declaration is mandatory:

```text
--cavity-neck-internal radius_m length_internal_m resistance_internal_Pa_s_m3 azimuth_rad z_m
```

The length includes the duct and any declared **interior** end correction; the
resistance includes only retained internal losses. Both **exclude exterior
radiation**, which the BEM-derived load supplies. No correction is guessed from
the radius, subtracted from an effective length, or added a second time. The
numbers below are explicit illustrative inputs, not measurements of a specimen.
The axial coordinate is measured from the top head, as for `--cavity-neck`.

From the workspace root, this selects both stretching heads and a bounded
40..500 Hz acoustic demonstration:

```sh
(set -C; cargo run --release -p fs-couple --example percussion -- \
  drum-stretch-mic 192 20 --cavity-modes --coupled-vent-radiation \
  --cavity-neck-internal 0.003 0.008 1000 0.35 0.08 \
  --radiation-spec crates/fs-couple/examples/percussion/coupled-vent.radiation \
  --analytic-newton --microphone-right -0.12,0.05,0.4 > coupled-vent.wav)
```

The run still has to pass its native modal, contact, BEM resolution, multiport
passivity and independent fit-audit gates. This command and the new integration regressions have not been
executed in the editing environment, which lacks Cargo/Rust. A declared band is
not a convergence result; failed fitting does not switch to prescribed flow or
invent a replacement impedance.

## One physical source map and one playing clock

Each head's original modal normal-velocity row is retained. The additional mouth
row is normalized by its actual polygon area so its integrated outward flow is
`Q = volume_weight * p_neck`, with the same neck coordinate used by the cavity.
The original closed computational boundary is preserved; a prescribed Neumann
mouth closes it acoustically but is not a rigid cap. The geometry adds no neck
inertance or damping.

For every source, BEM pressure is integrated against every source row to obtain
the **full** impedance. Head-mouth cross-reactions are retained, including their
signs. The full Hermitian power form must pass the existing check; positive
self-resistances alone are not sufficient. The passive realization then attaches
through `CavityCoupling::attach_exterior_radiation`, not a parallel solver or a
one-step-delayed pressure force.

The resisting mouth pressure obeys
`F_neck = -volume_weight * p_exterior`, hence
`F_neck * p_neck = -p_exterior * Q`. Its acoustic storage and dissipation are
already part of the joint energy ledger. They must not be added again as an
external force or a second loss. Internal cavity pressure and the mouth's load
pressure remain distinct from pressure observed at a microphone.

The same BEM scene supplies the mono or stereo pressure filters. Both receiver
channels observe one loaded trajectory. Preparing filters neither strikes the
drum nor advances its clock; failed/cancelled mechanical trials retain the
existing rollback contract. Analytic Newton, optional impact subdivision and
player-force programs remain execution controls on that same state.

## Preserved options and refusal boundaries

`--cavity-neck` keeps its original effective-impedance, zero-exterior-pressure
semantics. `--prescribed-vent-radiation` with that declaration keeps the existing
one-way observation model. It cannot be combined with feedback. A coupled
request cannot reuse the old effective declaration, omit a mouth source, expose
an interior standing-wave coordinate as an exterior source, or use a source
column twice.

Coupled vent audio is available on the nonlinear-capable `drum`, `drum-stretch`,
`snare` and `snare-off` microphone/WAV families. It does not silently convert the
linear-only `drum-modal` image. Snare wire count, contacts, second-stick addresses
and any selected head-stretching law remain explicit. The compact-neck chart is
checked against the requested acoustic upper frequency, not only the cavity
basis or the default observer window.

The current mouth must lie inside one fixed planar sidewall strip. Moving
apertures and elastic barrels with necks remain refused. The boundary, modal,
source-count, fit-order, pole, time-step and work limits are unchanged. Turbulent
jets, finite-amplitude exterior flow, measured material/neck identification,
full-band convergence and real-time qualification remain outside this increment.

## Focused native checks

```sh
cargo test --release -p fs-couple --example percussion coupled_selection -- --test-threads=1
cargo test --release -p fs-couple --example percussion feedback:: -- --test-threads=1
```

The new checks cover explicit declaration/admission, omitted and duplicate mouth
columns, accidental interior sources, out-of-chart frequency bands, the actual
two-stick stretching-head/snare/cavity/neck composition, work-conjugate exterior
pressure, total energy and refusal/retry. Two checks additionally execute actual
BEM projection and loaded stereo PCM preparation. The declared passive-load
fixture isolates mechanical coupling; it is not evidence that the separate BEM
fit succeeds. Native execution of all these checks is still required.
