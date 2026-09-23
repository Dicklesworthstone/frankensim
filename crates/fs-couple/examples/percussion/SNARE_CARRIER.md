# Force-driven snare engagement with a moving wire carrier

`--snare-carrier carrier.fsc` adds one finite-mass translating rail carrying
BOTH ends of the complete snare bank. An authored force can withdraw the wires
from the resonant head and drive them back while their existing vibration,
head motion, contact loss and cavity state continue. It does not rewrite gaps,
teleport endpoints, clear modes or substitute an output mute envelope.

The option works on all six `snare` / `snare-off` mechanics, WAV and microphone
commands. It independently composes with `--snare-spec`, linear or stretching
wires, `--head-stretching`, two sticks and their separate force files, fixed
mufflers, enclosed air/drag, admitted neck mechanics and existing stereo or
explicit prescribed-vent observers. Compliant mutes on snares remain refused.
Omitting the carrier keeps the original fixed-end bank and numerical image.

```sh
# Illustrative moving carrier; no hardware calibration is claimed.
(set -C; cargo run --release -p fs-couple --example percussion -- \
  snare 10000 --snare-carrier crates/fs-couple/examples/percussion/estimated-snare-carrier.fsc \
  --analytic-newton --impact-substeps 8 511 --strike-speed-m-s 0 \
  > carrier.csv)

# The SAME ongoing mechanics through two existing spatial pressure receivers.
(set -C; cargo run --release -p fs-couple --example percussion -- \
  snare-mic 4800 20 0.08 0.05 0.35 \
  --snare-carrier crates/fs-couple/examples/percussion/estimated-snare-carrier.fsc \
  --head-stretching --analytic-newton --impact-substeps 8 511 --cavity-modes \
  --microphone-right -0.12,0.05,0.4 > moving-snare.wav)
```

These are supported entry points, not claims of completed native renders or
real-time performance. Dense coupled Newton can be expensive for all 20 wires.
Every existing energy, contact, slope, clock, acoustic-fit and work budget still
applies. The explicitly selected substep budget can be exhausted; nothing is
retuned or clamped to guarantee a recording.

## Complete SI input

A UTF-8 file of at most 64 KiB begins with `frankensim-snare-carrier-v1`.
Exactly one `support` and one `initial` record are required. All force records
are required in increasing time order. Blank lines and `#` comments are allowed.
Unknown, duplicate, missing, nonfinite or invalid records refuse without defaults.

```text
frankensim-snare-carrier-v1
support,0.02,50,0.05,0.005,0.2
initial,-0.00002,0
force,0,0
force,0.002,0.05
force,0.008,0.05
force,0.012,-0.05
force,0.018,-0.05
force,0.02,0
```

`support,mass_kg,stiffness_N_m,damping_Ns_m,maximum_travel_m,maximum_wire_slope`
supplies the rail/actuator effective mass ONLY; the code adds the actual wire
mass once. Mass is positive; ground spring and drag are nonnegative. The spring
rests at zero carrier position. The symmetric travel limit is positive and the
relative-wire slope bound is in `(0,0.3]`. These are refusal limits, not stops,
clamps, contact springs or controlled trajectories. A supplied nonlinear wire
card's slope bound also applies; the stricter of the two is used.

`initial,position_m,velocity_m_s` initializes the carrier AND every wire in the
same rigid translation. Relative wire modes start at rest. Positive position
is downward, away from the resonant head. Thus an unchanged gap g has initially
available separation `g + position_m`. `snare-off` still selects its supplied
larger baseline gap; the carrier then moves relative to that installation.
Negative initial separation declares contact energy, not a solved equilibrium.
There is no inferred gravity load, latch, spring preload or hand position.

`force,time_s,force_N` supplies positive downward withdrawal force or negative
upward engagement force on the rail. The existing force-program implementation
integrates piecewise-linear forces over each mechanical interval. At least two
knots, nonnegative strictly increasing times and ZERO first/last force are
required; the force is zero outside the program. Its whole duration must fit
inside the requested render. No silent truncation or instantaneous impulses.
A force reversal is not a prescribed impact or disengagement timestamp.

## Reciprocal mechanics rather than moving contact gaps

For each mass-normalized relative sine coordinate q_n, the physical wire field
is `w(s) = X + sum(phi_n(s) q_n)`. Its kinetic energy contains cross terms
between rail velocity and relative motion. The implementation retains them
exactly with `z_n = q_n + c_n X` and `z_C = sqrt(M_r) X`, where
`c_n = integral(mu phi_n ds)` and
`M_r = M_rail + sum(mu L) - sum(c_n^2)`. The unresolved sine tail still carries
its rigid-translation mass; it is not discarded or added as another resonance.
The complete kinetic energy is then `sum(zdot^2)/2` in the original time owner.

The existing prestressed-beam/Kirchhoff-Carrier potentials act on the recovered
relative coordinates. Their transpose supplies the carrier reaction. Contact
rows use the same transformation to observe absolute wire motion. Relative
modal damping likewise becomes a reciprocal positive resistance. A force F on
the rail enters only the exact `F/sqrt(M_r)` port and contributes its conjugate
work to the existing ledger. No guessed inertia, duplicate tension spring,
additional contact law or second time integrator is introduced.

The rail adds exactly ONE coordinate to the full bank, after all relative wire
coordinates and before air inertia. Head/stick coordinates remain unchanged.
The bank is one coupled body for mass normalization, not a homogenized wire:
every strand retains all its own relative modes, contact rows and stretching
energy. The fixed-end per-strand body indices do not apply inside this combined
body; use `support_observation` and `support_force_port` instead. Wire and rail
coordinates have zero direct cavity-volume and acoustic-source weights; they
excite exterior sound through head motion.

Carrier selection requires the existing nonlinear-capable coupled solver even
when its wire energy is quadratic. `--prepared-nonlinear`, `--analytic-newton`
and `--impact-substeps` select its numerical realization. A linear-only modal
image is never allowed to silently drop the off-diagonal coupling. Every
accepted-state travel/slope failure or cancelled solve preserves state, felt
history, all staged force programs and the mechanical clock.

CSV reports `snare_carrier_position_m`, `snare_carrier_velocity_m_s` and the
existing wire slope/tension/stretch-energy columns. Stretching energy is ALREADY
included in total storage. `player_work_j` includes ALL supplied player ports,
not just the rail. WAV still uses the original pressure observer and scaling;
no direct rail sound or output normalization is appended.

## Scope and focused checks

This is a one-axis rail with both endpoints translating together, not a fully
resolved throw-off lever. It does not model tilting endplates, changing endpoint
spacing/tension, latches, wire sliding, torsion or coil contact. Ground spring
and damping are supplied effective properties, not measured hardware. Radiation
remains one-way and band-limited; neither acoustic feedback nor specimen realism
is established by adding a moving carrier.

Six library tests cover spatial kinetic energy, force/work scaling, conservative
derivatives, transformed damping, contact and force-driven withdrawal, exact
retry and physical limits. Four example tests cover input/command admission,
all 20 carried strands with either wire/head law, three player ports on one
cavity/contact/substep timeline, and unchanged unselected behavior.

```sh
cargo test --release -p fs-couple --lib render::plate::impact::supported::tests
cargo test --release -p fs-couple --example percussion snare::carrier::tests -- --test-threads=1
```

The new tests require native execution. Source review and independent arithmetic
checks are not Rust test results or evidence of a completed audio render.
