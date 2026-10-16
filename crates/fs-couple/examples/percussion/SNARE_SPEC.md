# Supplied snare banks and geometric wire stretching

`--snare-spec wires.fsn` supplies the physical wire-bank inputs to every
`snare` / `snare-off` command, including `-wav` and `-mic`. The same head meshes,
reciprocal distributed contacts, enclosed air and acoustic observers remain in
use. It is not a preset sound, a sampled buzz, or an output equalizer.

Without this option, the original estimated 20-strand bank remains unchanged.
`estimated-snare.fsn` reproduces those same estimated values and explicitly
selects linear wires. It is not a manufacturer specification or calibration.

```sh
# Supplied reference-equivalent bank, original linear-head modal mechanics.
(set -C; cargo run --release -p fs-couple --example percussion -- \
  snare 4096 --snare-spec crates/fs-couple/examples/percussion/estimated-snare.fsn \
  --strike-speed-m-s 2 > supplied-snare.csv)
```

## Complete SI format

The UTF-8 file is at most 64 KiB, starts with `frankensim-snare-spec-v1`, and
has exactly one of each following record. Empty lines and `#` comments are
allowed. Record order after the header is arbitrary. Unknown, duplicate,
missing, nonfinite, malformed or physically invalid inputs refuse; no value
is completed from the stock bank.

```text
frankensim-snare-spec-v1
bank,20,8,12,0.30,0.04
coil,0.00015,0.00055,0.00085,7800
mechanics,0.7,0.000001,4
contact,0.00002,0.003,500000000,1.5,0.05
stretching,off
```

`bank,strands,modes_per_strand,contact_cells,length_m,width_m` retains the
existing centered, parallel fixed-end layout. The width separates the first
and last strand; a one-strand bank is centered. Existing bounds remain: 1–24
strands, 1–16 modes per strand, and at least one more contact cell than modes,
up to 20 cells. The combined host also counts every head, striker and air mode
against its existing capacity. It refuses excess work rather than dropping wires.

`coil,metal_wire_radius_m,coil_centerline_radius_m,pitch_m,density_kg_m3`
identifies the mass per axial metre through the existing helical arc-length
calculation. Zero coil radius represents a straight wire. Metal radius, pitch
and density are positive. **This geometry does not identify installed tension,
transverse flexural rigidity, axial rigidity, damping or contact stiffness.**

`mechanics,tension_per_strand_n,bending_per_strand_n_m2,damping_per_s` supplies
those independent effective properties. Tension is positive; bending and modal
viscous rate are nonnegative. The original prestressed-beam frequencies include
both tension and bending. No MIDI pitch or observer gain enters this calculation.

`contact,engaged_gap_m,disengaged_gap_m,stiffness_per_length,exponent,chi_s_m`
uses the existing line-distributed power law and Hunt–Crossley loss. Stiffness
has units N/m^(exponent+1); physical contact weights are metres, not normalized
fractions. The exponent is at least one, stiffness and chi are nonnegative.
The disengaged gap must be larger than the engaged gap. `snare-off` uses that
SUPPLIED gap, not a hardcoded stock offset. Negative engaged gap explicitly
means initial interference/contact energy, not a solved preload equilibrium.

`stretching,off` selects the original linear filament. Alternatively,
`stretching,axial_rigidity_n,maximum_slope` adds geometric tension modulation
using the existing `fs-nlmodal` Kirchhoff–Carrier stress channel. Effective
axial rigidity E*A is supplied in newtons, not N/m or initial tension. It is
nonnegative; zero is useful for controlled nonlinear-image comparisons. The
slope validity bound is in (0,0.3]. For example, `stretching,100,0.2` declares
an illustrative 100 N axial rigidity; it is NOT measured coil data or inferred
from a bulk steel modulus. Supply identified effective data for such a claim.

## Actual nonlinear wire mechanics

For the retained physical sine amplitudes Q_n and k_n=n*pi/L, geometric added
strain is `sum(k_n^2 * Q_n^2)/4`. Tension becomes `T0 + E*A*added_strain` and
stretching storage is `E*A*L*added_strain^2/2`. This one positive quartic channel
couples a strand's partials through its physical extension. It does not add
another copy of initial tension, modify bending, or retune the audio output.
Each strand keeps its own state; it is not a single homogenized noise source.

The existing impact owner advances wires, drumheads, contact loss, sticks and
air together. Preparation, analytic tangents and internal substeps operate on
that SAME storage. No new time integrator or material family is introduced.
The continuous-span bound `sum(abs(Q_n)*k_n)` is checked initially and before
accepting every step. It is conservative when modal slopes cancel, but cannot
miss a peak between contact stations. Rejection does not clip displacement,
clear vibration/history, advance a player-force program, or consume an output
tick. Internal Newton trials may extrapolate; this is a model-validity bound,
not a temporal error certificate or between-step overshoot guarantee.

Wire and head nonlinearity are independently selected. A supplied stretching
wire bank uses nonlinear mechanics even without `--head-stretching`; that
flag additionally selects both existing nonlinear head potentials. Ordinary
linear-wire/linear-head snares retain the original prepared modal/contact image.
`--prepared-nonlinear`, `--analytic-newton` and `--impact-substeps` can be used
when either physical nonlinearity is selected. A linear-only solver never
silently drops a supplied quartic channel.

```sh
# Here wires.fsn explicitly includes a stretching,E*A,slope row.
(set -C; cargo run --release -p fs-couple --example percussion -- \
  snare 4096 --snare-spec wires.fsn --analytic-newton --impact-substeps 8 511 \
  --head-stretching --cavity-modes --strike-speed-m-s 2 > stretching-wires.csv)

# The same supplied bank in the existing sealed stereo microphone path.
(set -C; cargo run --release -p fs-couple --example percussion -- \
  snare-mic 4800 20 0.08 0.05 0.35 --snare-spec wires.fsn \
  --analytic-newton --impact-substeps 8 511 --head-stretching --cavity-modes \
  --microphone-right -0.12,0.05,0.4 > stretching-wires.wav)
```

These are entry points, not claims of completed native renders. Two stick
launches and independent force files, supplied drum specifications, fixed
mufflers, cavity losses and admitted necks keep their existing behavior.
The explicit prescribed-vent observer remains one-way; this wire model does
not supply radiation backreaction. Compliant mutes on snares remain unsupported
rather than being silently omitted.

CSV adds `snare_max_slope_bound`, `snare_max_tension_n` and
`snare_stretching_energy_j` only when wire stretching is selected. The last is
already included in total storage; do not add it again. Head and wire body
indices remain unchanged. Wires drive exterior sound through the head; no
direct wire radiation channel or synthetic buzz is appended.

## Scope and checks

The wire law remains a fixed-end, uniform, one-polarization, moderate-slope
model with quasi-static averaged axial tension. It does not resolve individual
coils, torsion, longitudinal inertia, inter-wire friction, plasticity, end-plate
compliance or a full throw-off linkage. `--snare-carrier` independently adds a
force-driven finite-mass translating support carrying both wire endpoints; see
[SNARE_CARRIER.md](SNARE_CARRIER.md). Input values apply uniformly across the bank. Modal/contact resolution and acoustic bandwidth still need convergence
and specimen validation; nonlinear harmonics can exceed the retained band.
No calibrated realism, full-band adequacy or real-time deadline is claimed.

Five core regressions exercise continuous extension/energy, the original bending
and zero-rigidity limits, analytic derivatives, finite-amplitude period refinement,
contact backreaction/work, and exact retry/slope refusal. Five example regressions
exercise full input admission, reference-equivalent motion, the complete bank
with either head law, and supplied nonlinear wire/head/cavity/player composition.

```sh
cargo test --release -p fs-couple --lib render::plate::impact::string::tests
cargo test --release -p fs-couple --example percussion snare::spec::tests -- --test-threads=1
```

These native tests require execution. Independent period/tangent arithmetic
checks and source inspection alone do not establish that the Rust tests pass.
