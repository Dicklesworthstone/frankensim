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

The resolver admits the entire time interval and every exact class through the
ordinary typed material queries. All eight claims must retain one observation
and identical support. Other filling values, mill dimensions, ball sizes, speed,
charge mass, wet grinding and plant qualification refuse. Original observations,
claims, caveats, pack hash and receipts survive in
`SelectionBreakageModel::austin_calibration_source()` and its prepared passes.

`xtask/tests/mining_austin_sources.rs` checks the printed Table 2 selection rates,
unit conversion, independent analytical first and second population solutions,
all-column mass closure, growing sink inventory and condition refusals. Unit
tests additionally check an exact three-class solution, element transport,
source mixing, incomplete time support and invalid grid/coefficient refusal.

This supplies a bounded laboratory PSD model. It does not establish plant
throughput, power demand, scale-up, grinding-media consumption, liberation or
recovery, and it does not derive a full PSD from a lone Bond index or t10 value.
