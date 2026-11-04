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

`at_residence_time` controls residence and reports the resulting hold-up.
`at_dry_holdup` instead takes a positive target dry mass and solves residence
inside an explicit time bracket. Every trial passes through the same original
continuous circuit, including its source and balance checks. The caller
declares mass tolerances and a finite iteration budget. An unbracketed target,
exhausted budget, unresolved numerical interval or unsupported operating state
refuses; the solver does not expand the bracket or relax the source domain.

A fully returned breakable class can impose a positive minimum inventory even
as residence approaches zero. Positive target mass alone does not guarantee
a feasible steady circuit; the evaluated bracket must actually contain it.

Holding dry inventory fixed makes throughput affect residence. For the full
return example above, the exact two-class relation is
`H/f = tau*(0.4 + 0.1*tau)/(0.4*(0.05 + 0.1*tau))`.
With all fresh solids in the coarse class, H = 1000/3 kg and fresh feed
10 kg/s give 10 s residence. Doubling fresh feed
to 20 kg/s reduces the solved residence to about 3.5941 s. The returned circuit
still reports its actual circulating flow, product size distribution, assays
and water; fresh throughput is not substituted for circulating throughput.

Neither continuous-mill control calibrates selection/daughter data from a
Bond work index, predicts liberation, or attaches a generic ore name to an
unsupported classifier card. These inputs remain explicit and source-qualified.

## Bond energy and power at the solved operating state

`population_balance::circuit::power::ContinuousBondCircuit` takes the original
continuous circuit and an independently declared `BondGrindingCard`. It offers
the same prescribed-residence and fixed-dry-hold-up controls. Each call solves
the original circuit internally before evaluating fresh-feed F80 and the P80
of **overflow plus purge**. The original work-index card must admit that actual
size reduction. Estimated net pinion power is specific energy times **fresh
dry throughput**; circulating ore and carrier water are not additional fresh
tonnage. The original circuit, sources, streams and hold-up diagnostics remain
available from the returned point.

The example uses an explicitly authored 15 kWh/metric-tonne work index and
prints the resulting power alongside the fixed-inventory response. In its
two-bin illustration, the discrete P80 remains 100 micrometres even though
fine-product fraction changes. Estimated power therefore changes from
158.162338 kW at 10 kg/s fresh flow to 316.324676 kW at 20 kg/s.

The three-bin analytical regression exercises the case that an unchanged P80
would miss. Its synthetic rates, classifier and 15 kWh/metric-tonne index give:

| Dry inventory, kg | Fresh dry flow, kg/s | Solved residence, s | External P80, micrometres | Estimated pinion power, kW |
| ---: | ---: | ---: | ---: | ---: |
| 752.760129631 | 10 | 40 | 100 | 270 |
| 752.760129631 | 20 | 15.775944376 | 200 | 223.675324 |

The discrete P80 jump lowers the Bond estimate in this particular fixture.
It does not imply a plant power-draw prediction. A smooth monotone inverse
from power to capacity is not available from these discrete quantiles, and
the adapter does not manufacture one. Selection frequencies, daughter
distributions, work-index test applicability, equipment power draw, drive
efficiency and auxiliary loads require their own physical inputs. A finer
overflow cannot conceal an unsupported or unreduced combined purge/product.

The Bond relation and fresh-feed accounting follow the existing
[GMG 2021 grinding-circuit guideline](https://gmggroup.org/wp-content/uploads/2024/07/GUIDELINE__Determining-the-Bond-Efficiency-of-Industrial-Grinding-Circuits_2021.pdf),
sections 4–5; the population-balance fixture supplies no measured calibration.

Seven focused regressions cover independent two-/three-class answers, mixed
assays, purge and no-return identity, source and circulation limits, fed and
unforced traps, water trapping, permutation invariance and residence response.
Seven additional fixed-hold-up regressions check independent residence and
product-size answers, carrier-water exclusion from dry inventory, purge and
no-recycle limits, source admission, infeasible targets, iteration exhaustion,
overflowed tolerance and unresolved floating-point intervals.
Six Bond-adapter tests add changing-P80 power, purge, fresh-tonnage accounting,
original work-index/source restrictions and preserved physical-owner refusals.

## Startup, feed changes and finite mill inventory

`ContinuousPopulationBalanceMill::advance_inventory` now advances one local
mill's `transient::MillInventory` through a constant-feed interval. Starting
empty, changing feed, and washing out an existing charge retain material held
inside the mill. The caller carries the returned inventory into the next
interval; every size class keeps its explicit elemental masses and water has
its own extensive account. For dry or elemental class mass `H_i`:

```text
dH_i/dt = f_i - (S_i + 1/tau)*H_i + sum_j(b_ij*S_j*H_j)
```

The zero-accumulation limit is exactly the existing steady mill equation.
Water follows `dH_water/dt = f_water - H_water/tau`. Constant source rates,
daughters, feed and common residence are declared for each interval. Positive
backward-Euler substeps solve from coarser to finer classes and accumulate
discharge as `dt*H_new/tau`. The result returns final inventory, instantaneous
outlet, integrated feed/discharge and actual dry, water and elemental closure
residuals. Account values are never adjusted to force closure.

The caller supplies maximum step, finite substep budget and mass-balance
tolerances. Time integration is first order; balance admission is not a bound
on trajectory error. Cancellation or an unsupported, unrepresentable or
unbalanced operation returns no partial response and leaves the inputs intact.
Actual feed flow, residence and defined feed solids fraction enter the original
source query. Empty-feed washout requires source support at zero dry flow; if
the source requires feed solids fraction, the undefined empty-feed ratio causes
refusal. Inventory composition is never substituted for feed composition.

```bash
cargo run -p fs-material --example transient_grinding
```

This synthetic example starts empty, feeds 10 kg/s dry ore for 30 seconds,
doubles feed for 30 seconds, then washes out for 30 seconds at a declared
10-second common residence. It prints retained and discharged solids, water,
copper and balance residuals. The G1 tests compare startup/washout and a driven
two-size system to independent analytical solutions, demonstrate first-order
refinement and retain the original steady solution as a fixed point. G0/G3/G4
checks cover separate elemental/water balances, feed changes, replay, source
restrictions, cancellation and bounded work.

This is one local perfectly mixed inventory. Dynamic classifier recycle,
size-dependent discharge, liberation, varying kinetic coefficients, equipment
power draw and plant calibration remain separate physical inputs/models. The
source-backed Bond test rows do not supply any of those missing calibrations.
