# Continuous mill, classifier, recycle and purge

`fs_material::comminution::population_balance::circuit::ContinuousGrindingCircuit`
connects the continuous perfectly mixed population-balance mill to an original
classifier card and an explicit underflow recycle/purge split. It retains
size-resolved dry solids and every assayed element, with a separate carrier
water account. The result contains both the existing `CircuitSolution` stream
accounts and the original `ContinuousMillSolution`, including source queries
and class/water hold-up.

## Physical inputs and closure

The caller supplies fresh feed (including makeup water), common residence
time, original selection frequencies and strictly finer daughter fractions,
classifier conditions/card, return fraction and circuit closure tolerances.
For fresh class flow f_i, discharge p_i, returned fraction s, classifier
underflow partition r_i, selection frequency S_i and residence tau, the steady
equation is:

```text
(1 - s*r_i + tau*S_i) * p_i = f_i + tau * sum_j(b_ij*S_j*p_j)
```

Only coarser parents j contribute, so a deterministic triangular solve closes
every dry-mass and elemental channel. The implementation evaluates the escape
coefficient as `purged + returned*overflow_i`, retaining independently
represented small outlet fractions. Water obeys its own classifier/recycle
balance. A fed class or water inventory with no outlet refuses; unforced
trapped circulation is set to zero under the existing circuit convention.

The predicted discharge determines returned underflow and actual mill feed.
The original continuous mill is then evaluated on that circulating feed, and
the original classifier on its discharge. Actual throughput, residence and
solids fraction must fit their source support. Caller labels cannot override
these solved coordinates. The original circuit closure owner checks both
local reconstruction and the combined external overflow/purge accounts.
Excessive circulation or failed balances return no accepted circuit result.

## Runnable example

```bash
cargo run -p fs-material --example continuous_grinding_circuit
```

The example explicitly declares synthetic selection and classifier data;
it is not a measured ore calibration. Fresh dry flow is 10 kg/s, with 10 s
prescribed residence. The coarse class is 200 micrometres, the fine class
100 micrometres; all selected coarse mass transfers to fines at 0.1 /s.

| Case | Mill dry feed, kg/s | Recycle/fresh ratio | Overflow, kg/s | Purge, kg/s | Dry hold-up, kg |
| --- | ---: | ---: | ---: | ---: | ---: |
| All underflow returned | 33.333333333 | 2.333333333 | 10.000000000 | 0 | 333.333333333 |
| 75% returned, 25% purged | 21.888790821 | 1.188879082 | 6.037069726 | 3.962930274 | 218.887908208 |

Overflow alone does not equal fresh feed when purge also leaves. Copper and
water close separately. The example prints the reported closure residuals
alongside these streams, so the purge cannot disappear from the accounting.

This API controls residence and reports the resulting hold-up. It does not
solve residence from fixed hold-up, infer mill power or equipment capacity,
calibrate selection/daughter data from a Bond work index, predict liberation,
or attach a generic ore name to an unsupported classifier card. Those inputs
remain explicit and source-qualified.

Seven focused regressions cover independent two-/three-class answers, mixed
assays, purge and no-return identity, source and circulation limits, fed and
unforced traps, water trapping, permutation invariance and residence response.
