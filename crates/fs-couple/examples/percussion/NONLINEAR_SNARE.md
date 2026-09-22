# Stretching heads, complete snare bank, one nonlinear solve

Add `--head-stretching` to any `snare` or `snare-off` command, including its
`-wav` and `-mic` forms. Both heads then use the existing geometry-derived,
statically relaxed von Karman membrane potential. All 20 strands, 160 wire
coordinates and 240 distributed wire/head contacts remain in the same physical
state. The wires' existing nonzero Hunt–Crossley loss is retained. Omitting the
flag keeps the original linear-head prepared modal/contact image unchanged.

```sh
# Nonlinear-head mechanical onset with all wires and spatial enclosed air.
(set -C; cargo run --release -p fs-couple --example percussion -- \
  snare 4096 --head-stretching --analytic-newton --impact-substeps 8 511 \
  --cavity-modes --strike-speed-m-s 2 --strike-position-m 0.06 0.01 \
  > stretching-snare.csv)

# The same mechanical selection with two existing finite-point receivers.
(set -C; cargo run --release -p fs-couple --example percussion -- \
  snare-mic 4800 20 0.08 0.05 0.35 --head-stretching --analytic-newton \
  --impact-substeps 8 511 --cavity-modes --cavity-drag-per-s 25 \
  --microphone-right -0.12,0.05,0.4 > stretching-snare.wav)
```

These are example invocations, not a claimed completed render or calibrated
sound. The explicit internal-solve budget can still be exhausted; no force,
energy, membrane-slope or fit tolerance is relaxed to manufacture an output.
`--prepared-nonlinear` retains scratch for the original finite-difference
Newton image. `--analytic-newton` additionally selects exact storage and
contact-loss derivatives. Without either option, the allocating nonlinear
reference uses the same contact-loss equation. These numerical choices do not
replace the nonlinear heads with a linear model. On ordinary snare commands
without `--head-stretching`, nonlinear preparation continues to refuse rather
than silently changing their physical model.

## Contact loss inside the actual mechanical equation

The existing `fs-dcontact` contact energy still supplies conservative forces.
Its existing dissipative increment is evaluated at trial midpoint penetration
and the full Gonzalez discrete effort. A point contributes positive resisting
momentum flow `b * weight * K * penetration^alpha * max(chi * rate, -1)`, where
`rate` is the same effort projected through the contact row `b`. Each point's
conjugate power is nonnegative. The unloading cap and all supplied contact
coefficients are unchanged; no new felt, wire or restitution law is introduced.

`fs-phs::StepWorkspace` includes this flow in its original implicit equation,
not a previous-sample force, endpoint correction or a second integration stage.
Finite-difference probes include the full flow. Analytic Newton differentiates
both penetration and the corrected discrete effort, including the cap branch.
The resulting contact dissipation enters the same accepted energy ledger as
modal loss, fixed mufflers, gas drag and any Kelvin/felt loss. It is not reported
as external player work. Invalid or active dissipative trials refuse; the
original complete-window Hamiltonian/work check remains independent.

Stick, head and wire addresses are unchanged. A second stick follows the heads,
all wires follow it, and distributed air inertia is appended last. Supplied drum
geometry/materials, independent stick-force files, fixed mufflers, sealed gas
drag and the existing bounded substep recovery compose with this selection.
All player ports advance once per accepted mechanical tick. A refused tick
leaves physical state, history, controls and output clocks unchanged.

The nonlinear host's explicit total-coordinate ceiling is now 256 so the full
wire bank fits; it never truncates strands to fit the old 64-coordinate limit.
This is a bounded work envelope, not a wider physical frequency band. Head
reduction, contact sampling, acoustic source count and all energy/validity
limits remain unchanged. Dense Newton/LU can be expensive at this size; no
measured real-time, deadline or full-polyphony claim is made.

## Observations and limits

Stretching CSV adds the two head-slope diagnostics and total stretching energy.
That energy is already in total storage. The existing one-way BEM observes head
motion, not an invented direct wire or gas audio channel. Mono/stereo receiver
positions, pressure scaling, sample counts and propagation keep their meanings.
Changing `snare` to `snare-off` selects the original larger wire clearance; it
does not delete the wires or replace their equations with a mute envelope.

A neck remains available for mechanics CSV with `--cavity-modes`; vented WAV
and microphone output still refuse because aperture radiation is missing. The
head-mode and acoustic-fit windows remain unchanged. Nonlinear harmonics may
exceed them. No modal/boundary convergence, measured specimen agreement,
room/lid scattering, radiation feedback or full-band sound adequacy is implied.

## Focused regression coverage

Four time-owner regressions exercise implicit static-resistance equivalence,
analytic state/effort derivatives, nonlinear work balance and exact retry.
Four contact/impact regressions compare the existing loss law, passivity and
piecewise tangents, rebound and reciprocal motion, reference/analytic execution,
refusal/clock preservation and the bounded wire-sized state. Four example
regressions cover CLI selection, complete 20-strand construction, actual wire
exchange with stretching heads, and two driven sticks with cavity/substeps.
The small runtime fixture explicitly declares installed wire interference to
exercise contact immediately; the production 20-micrometre gap is unchanged.
Construction coverage alone does not establish a completed BEM/audio render.

```sh
cargo test --release -p fs-phs --lib prepared::dissipation::tests
cargo test --release -p fs-couple --lib render::plate::impact::contact_loss_tests
cargo test --release -p fs-couple --example percussion nonlinear_snare::tests -- --test-threads=1
```

Native execution is required to establish these authored regressions pass.
The editing environment has no Cargo/rustc; independent arithmetic or source
inspection is not a substitute for those native results.
