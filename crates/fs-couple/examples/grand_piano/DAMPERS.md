# Finite-footprint piano damping

`--dampers estimated` or `--dampers pads.fspd` selects spatial viscous pads
on the existing speaking strings and loaded soundboard. Omitting the option
keeps the original point-drag image. Both are approximations, not a full
falling-damper action or an identified wool-felt contact model.

```bash
cargo run --release -p fs-couple --example grand_piano -- \
  --preset steinway-d --dampers estimated --note 69 \
  --duration 6 --render release.wav

cargo run --release -p fs-couple --example grand_piano -- \
  --preset steinway-d --dampers pads.fspd --midi performance.mid \
  --duration 30 --render performance.wav
```

## What changes physically

The original point at 35% of string length is blind to partials with a node
there, including the 20th pinned-string partial. A finite pad instead integrates
local velocity **squared** across its span. Opposite-moving portions do not
cancel their dissipation, as they would when damping only the average velocity.
The full modal cross terms and moving-bridge displacement contribution remain.
Each speaking member of a unison is damped independently; duplex segments have
no direct pad force but continue to interact through the shared bridge.

Key hold, sostenuto capture and sustain use the existing engine controls,
including MIDI and SI CSV. A held or latched key lifts its pads. Sustain travel
retains the existing `(1-travel)^2` drag scaling; full sustain removes pad drag.
The spatial model replaces the point model, rather than stacking on top of it.
Explicit pad/free rows own the upper break and individual drag coefficients.
There is no output envelope, resonance reset, tail cutoff or pressure gain.
The original contact solver, loss accounting and failed-sample rollback remain
in use. Microphones receive the changed physical board velocities.

## Per-key input in SI units

The header is `frankensim-piano-dampers-v1`. Every key in the supplied string
scale must have exactly one row, even when the performance never strikes it:

```text
frankensim-piano-dampers-v1
# Example ONLY for a two-key scale, with key 69 longer than 0.14 m.
pad,69,0.12,0.14,0.4
free,90
```

`pad,key,start_m,end_m,drag_ns_m` specifies endpoints measured from the fixed
end of the speaking string. Require `0 < start < end < string length` and
finite positive drag up to 1e6 N s/m. Drag is the TOTAL coefficient per speaking
string across the pad, not a coefficient per metre and not divided by unison
count. The local dissipation is `c/(end-start) * integral(v(x)^2 dx)`.
`free,key` explicitly declares no damper on that course. A full preset therefore
needs 88 rows; the two-row example is not a complete full-keyboard specification.
Blank lines and `#` comments are accepted. Duplicate, missing, unknown, invalid
or out-of-span rows refuse without choosing defaults. Reads are bounded to 1 MiB.
The option requires `--render` and shares the existing direct input/output path
collision checks. A file literally named `estimated` can be selected as
`./estimated`; the bare word is reserved for the following estimate.

`estimated` uses a centre at 35% of speaking length, width equal to the smaller
of 40 mm and 8% of length, and 0.4 N s/m per string. Keys above MIDI 88 are free.
These are explicit estimates, NOT Steinway dimensions or measured felt losses.
Supplied geometry, tensions, hammer cards and acoustic observer are unchanged.

## Numerical scope

Positive midpoint quadrature retains at least eight stations per footprint,
with additional stations tied to the highest retained spatial partial. The
limits are 128 stations per string and four million projection terms; they
refuse rather than dropping modes. This guard is not a spatial error bound.
The entire shared-board port sequence is composed forward then backward using
exact dissipative rank-one velocity flows. This gives a passive, second-order
split approximation, not the exact full coupled damping exponential. Stored
energy does not change when installing a viscous specification. Preparation is
cold; stepping adds no allocation in this damping kernel.

The render reports actual pad/station counts and damper energy loss (already
included in total losses). No real-time, full-band, measured-decay, falling-pad,
preload, stick-slip, collision or hysteretic damper-contact claim is made.
Native regressions cover nodal damping, moving-bridge geometry, energy, pedal
and duplex behavior, contact-failure retry and MIDI/CSV render preparation.
Native Rust test execution and audible comparison are still required.
