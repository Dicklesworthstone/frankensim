# Played geometric string extension

`piano_exterior render` and `render-loaded` accept
`--string-stretching strings.fsps`. This selects the existing same-tick
Kirchhoff--Carrier string/hammer/bridge mechanics. It is not stretch tuning,
MIDI pitch bend, an output effect, or an independent nonlinear oscillator.

```sh
cargo run --release -p fs-couple --example piano_exterior -- \
  render-loaded settled.fss strings.csv body.obj acoustics.fspe piano.wav 6 \
  --string-stretching strings.fsps --modes 128 --substeps 8 \
  --hammer-footprints faces.fshp --dampers pads.fspd --midi score.mid
```

The optional controls compose with the same complete supplied hammer materials,
point/finite hammer faces, jack-force CSV, MIDI mapping, sustain, sostenuto,
una corda, spatial dampers and one-way or passive radiation feedback. Existing
inputs without this option retain the original linear-string image.

## Complete per-key SI input

```text
frankensim-piano-string-stretching-v1
# Illustrative two-key study, NOT identified Steinway axial rigidities.
stretch,69,150000,0.2
linear,72
```

The file must cover **every key in the supplied scale exactly once**, including
unplayed sympathetic courses. A full preset needs 88 rows. `linear,key` explicitly
retains the original linear potential for that key. A `stretch` row supplies
key, positive finite axial rigidity **EA in newtons**, and a finite continuous
slope limit in `(0,0.3]`. The example is valid only with a two-key 69/72 scale.
Comments after `#` and blank lines are permitted. The file is bounded to 64 KiB.

EA is neither installed tension T, bending rigidity EI, nor spring rate EA/L.
In particular, an effective winding mass and bending rigidity do not determine
the axial rigidity of a wrapped string. The program does not infer or invent
that missing material input. One course's supplied EA applies to all its unison
members and existing duplex segments; the geometry and individual detuning
remain those of the original scale.

Missing files, duplicate/extra/missing keys, nonfinite values and invalid slope
limits refuse before structural eigenanalysis or acoustic BEM preparation.
An explicitly all-linear file is legal and retains the old numerical path.
`response` and `admittance` reject the option: their present linear harmonic
experiments cannot silently stand in for an amplitude-dependent response.

## What actually changes

The existing string owner supplies the extension channel. The piano adapter
includes both the fixed-interface string displacement and the moving bridge's
chord slope. Its work-conjugate forces act on the **same string and board
coordinates** as the existing contact solution. The nonlinear extension force
is resolved against all enabled hammer sites, their original felt/Prony
histories and the shared board before accepting one mechanical tick.

The unloaded small-amplitude frequencies, modal mass normalization, unison
layout, material loss and acoustic boundary projection are unchanged. During a
strike, geometric extension stores additional energy and changes physical
string tension. The complete energy balance includes that storage. No stored
energy is turned into an extra dissipative loss or an output pitch envelope.

Hammer/jack clocks, damper flow and acoustic splitting advance once per tick,
not once per nonlinear trial. The engine retains its original convergence,
material, energy and sample-rollback checks. A slope or convergence refusal
does not clip the string or publish a weaker linear substitute. Raising
substeps is an explicit caller choice, never automatic clock repair.

Both exterior receivers observe the same nonlinear mechanical trajectory.
Radiation-loaded playback retains the fitted acoustic storage and reaction;
selecting extension neither changes the fit tolerances nor bypasses a failed
passive fit. The report identifies whether nonlinear channels are active and
the selected material path.

## Boundaries

This is averaged, moderate-slope axial extension of the existing single
transverse-polarization string image. It does not resolve longitudinal waves,
a second transverse polarization, tension-dependent winding slip, full 3-D
hammer/string contact, or measured piano string materials. The acoustic fit
remains bounded to its admitted sampled band even when a nonlinear attack
contains higher-frequency energy. Retained modes, substeps and a supplied EA
are not a realism, convergence or real-time certificate.

The exterior regression target exercises strict complete-scale admission,
bitwise all-linear jack/pedal behavior, and actual one-way and radiation-loaded
stereo against direct calls to the same mechanical owner. It requires nonzero
extension, changed physical motion, felt/damper loss and the original combined
energy bound; native execution status is reported separately from these tests'
existence.
