# Experimental ball-mill selection and daughter distributions

`quartz-assiros-petrakis-2017-austin` supplies **three experimental calibrations
and 24 claims**: seven scalar Austin coefficients and one exact size-class curve
per calibration. `comminution::austin::resolve_austin_grinding` connects them to
the existing `SelectionBreakageModel`, so the same population-balance solver
now calculates time-dependent particle-size distributions from a real test source.

## Source and model scope

Petrakis, Stamboliadis and Komnitsas (2017),
[*Identification of Optimal Mill Operating Parameters during Grinding of Quartz
with the Use of Population Balance Modeling*](https://doi.org/10.14356/kona.2017007),
KONA 34, 213–223, provides the Table 2 fits and equations 11–15. The article is
CC-BY-4.0. Assiros quartz was tested at ball fillings of 20%, 10% and 6.7%, using
345 g dry charges, a 204 × 166 mm laboratory mill, 25.4 mm steel balls and 66 rpm.
Powder filling was 4%.

These are the **original Table 2 fits**. The paper reports a mismatch between
their modeled and measured product PSDs. Table 3's improved fit does not report
a complete amplitude, so it is not substituted or reconstructed. Beta was
fixed to 5.80; parameter covariance and test temperature/pressure are unstated.
No numerical uncertainty is manufactured.

The source selection fit spans parent sizes 0.212–3.35 mm. Daughter parameters
were fitted with the 0.850–0.600 mm feed fraction. This profile deliberately
uses only the six source sieve sizes **0.850, 0.600, 0.425, 0.300, 0.212 and
0.150 mm**, with the final node an unresolved-fines sink. Its 0–240 s
calculation window is a declared fitted-reference scope within the test-time
range, not a claim that every time and product size was measured or validated.

## Equations and the last class

The selection frequency is

`S_i = a_T * (x_i / 1 mm)^alpha / (1 + (x_i / mu)^lambda)`.

The compiler converts the printed amplitude from inverse minutes to inverse
seconds and mu from millimetres to metres. The source class coordinate is the
upper-size argument of this equation. A geometric-mean diameter is not silently
substituted, and `a_T` is not itself `S(1 mm)` when the rolloff correction differs
from one.

For descending source class indices, the cumulative daughter law is

`B_ij = phi*(x_(i-1)/x_j)^gamma + (1-phi)*(x_(i-1)/x_j)^beta`, for `i > j`.

Thus `B_(j+1),j = 1`. Each resolved daughter fraction is the difference
`B_ij - B_(i+1),j`. Using `x_i/x_j` in place of the published index would lose
part of each selected parent's mass. The final class receives the remaining
cumulative fraction directly; its selection rate is zero. The sink node is a
representative label, not a resolved distribution or mean size of its contents.
No coefficients are clipped or rescaled to force closure.

Class identity and discretization travel in the material source. A caller cannot
substitute an arbitrary grid, extend the calculation to smaller active particles,
or invent an unresolved-fines selection rate. Those changes require another
explicit source model.

## Run the source-backed PSD calculation

```bash
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/quartz-assiros-petrakis-2017-austin/manifest.tsv \
  --out /tmp/quartz-assiros-petrakis-2017-austin.fsmatpk
cargo run -p fs-material --example austin_grinding_reference -- \
  --material /tmp/quartz-assiros-petrakis-2017-austin.fsmatpk
```

The example starts from a normalized 0.850 mm source class and prints modeled
mass fractions after 30, 60, 120 and 240 seconds for each original ball filling.
It reports mass closure and the existing solver's nominal Poisson-tail bound.
That numerical bound does not include experimental uncertainty or model error.

## Grind to a specified source-bin fraction

`SelectionBreakageModel::time_to_target` now solves the contact time needed to
reach a cumulative product mass fraction from an actual `ParticleStream` feed.
The target names the coarsest included source class; every class with an equal
or smaller declared diameter is included. This avoids inventing an interpolated
sieve cut or an internal distribution for the unresolved-fines sink. The input
flow scale supplies the PSD weights, not a claim about batch mass or throughput.

The original daughter matrix admits only transfers toward finer classes, and
the inverse checks that no active target class can lose mass outside the target
set. Under the existing constant-rate population balance, cumulative target mass
is therefore nondecreasing. Every trial calls the original `prepare` and
`transform`; the returned product preserves their solids, element and water
accounts, source conditions and evidence.

The search starts at zero and is bounded by `maximum_contact_time`, a total
`maximum_forward_evaluations` budget and the original per-trial uniformization
budget. A nonzero result retains actual below-target and meeting-target points
whose time separation is at most `time_tolerance`. Each comparison includes twice
the original Poisson-tail bound plus an explicit `fraction_roundoff_allowance`.
The latter is a caller assumption about source/arithmetic rounding, not a
computed interval certificate. These allowances do not cover the source's
experimental uncertainty or its reported model mismatch. An ambiguous threshold,
unreachable target, source refusal, exhausted budget or cancellation returns no
accepted time. The result does not claim an exact earliest crossing.

The example also solves for 50% of the initial top-class mass to enter
`size-001` and finer source bins, using a one-microsecond time-bracket tolerance.
This is the original model's calculated threshold, not a measured grinding time.
The integration test independently checks it against top-class survival
`exp(-S_0*t)`, whose 50% crossing occurs at `ln(2)/S_0`.

Optionally append `--net-power-w DECLARED_WATTS` to the example command to quote
batch energy. `GrindTargetSolution::energy_at_constant_net_power` computes
`E = P_net*t` at both retained time endpoints and `E/m_dry` at the meeting time.
The power is an explicit external input: this source pack does not provide a
power-draw measurement or law. The example uses the admitted 0.345 kg dry charge,
and the API refuses a supplied mass that conflicts with an admitted
`dry-charge-mass`. It does not infer batch mass from PSD flow units. The quote
excludes idle time, auxiliaries and drive losses and does not recalibrate the
grinding rates when a different external power is supplied.

The resolver admits the entire time interval and every exact class through the
ordinary typed material queries. All eight claims must retain one observation
and identical support. Other filling values, mill dimensions, ball sizes, speed,
charge mass, wet grinding and plant qualification refuse. Original observations,
claims, caveats, pack hash and receipts survive in
`SelectionBreakageModel::austin_calibration_source()` and its prepared passes.

`xtask/tests/mining_austin_sources.rs` checks the printed Table 2 selection rates,
unit conversion, independent analytical first and second population solutions,
all-column mass closure, growing sink inventory, all three target-crossing
brackets and condition refusals. Unit tests additionally check an exact
three-class solution, mixed-feed target inversion, secondary daughter breakage,
element transport, declared-power energy, unresolved threshold and work-budget
refusals, source mixing, incomplete time support and invalid grid/coefficient
refusal.

This supplies a bounded laboratory PSD model. It does not establish plant
throughput, power demand, scale-up, grinding-media consumption, liberation or
recovery, and it does not derive a full PSD from a lone Bond index or t10 value.
