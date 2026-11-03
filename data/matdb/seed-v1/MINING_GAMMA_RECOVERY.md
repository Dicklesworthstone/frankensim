# Gamma-distributed recovery in batch tests and stirred-tank trains

`fs_material::recovery::RecoveryKinetics` and
`fs_material::recovery::continuous::ContinuousRecoveryCard` accept separately
sourced gamma-rate declarations. They describe a heterogeneous recoverable
population whose first-order rates follow a gamma distribution. A mineral name,
reagent identity or mean rate does not supply its shape, scale or calibration.

## Equations and physical meaning

Let `a > 0` be dimensionless gamma shape, `theta` the rate scale in s⁻¹ and
`Rinf` the recoverable fraction of the entering constituent. Mean rate is
`a * theta`; variance is `a * theta²`. Batch recovery at contact time `t` is

```text
R(t) = Rinf * [1 - (1 + theta*t)^(-a)].
```

This is the gamma Laplace transform tabulated by Vinnett and Waters,
[Minerals 10(10), 913 (2020), Table 1](https://doi.org/10.3390/min10100913).
Shape one is an exponential distribution of rates, rather than one fixed rate.

For `N` equal, full, perfectly mixed tanks with total active volume `V` and
actual carrier flow `Q`, mean residence is `tau = V/Q`. Each rate class keeps
its rate through every tank: its surviving fraction is
`(1 + k*tau/N)^(-N)`. The surviving mixture is the expectation of that
expression over the original gamma distribution. Equivalently,

```text
b = theta*tau/N
S = integral[0,infinity] y^(N-1) * exp(-y) * (1+b*y)^(-a) / (N-1)! dy
R = Rinf * (1-S)
unrecovered = (1-Rinf) + Rinf*S.
```

The Erlang residence distribution and first-order tank balance follow
[Fogler, Elements of Chemical Reaction Engineering, §17.4](https://websites.umich.edu/~elements/5e/17chap/Ch17_Web-Additional%20Material.pdf),
equations 17-4/17-5 and 5-15. Combining those equations with the gamma Laplace
transform gives the integral above. Evaluating batch recovery at `tau`, or
resetting the rate distribution after each tank, gives a different model.

## Source-card declarations

Both owners consume the existing `ConstitutiveModelCard`; no pack-schema
extension or fitted default is introduced.

| Item | Required declaration |
| --- | --- |
| Batch law | `gamma-rate-recovery`, version 1 |
| Stationary train law | `cstr-gamma-rate-recovery`, version 1 |
| `ultimate-recovery` | Dimensionless fraction in `[0,1]` |
| `gamma-shape` | Positive, finite, dimensionless `a` |
| `gamma-rate-scale` | Finite inverse-second `theta`; zero explicitly denotes the zero-rate limit |
| State | Schema version zero and `ZeroInternalState` |
| Sources | Original nonempty source hashes and valid provenance |
| Batch support | Finite typed `elapsed-time` bounds in seconds |
| Train support | Finite typed `mean-residence-time` bounds in seconds, and dimensionless integer `tank-count` bounds within `1..=64` |

Every additional validity coordinate, such as pH, grind size, temperature or
reagent dosage, must have finite typed bounds and be supplied at evaluation.
The owner supplies actual contact time, or derives actual mean residence and
tank count; conflicting duplicate coordinates refuse. A stationary train card
cannot use the batch `elapsed-time` axis. There is no automatic conversion of
finite batch-test support into stationary local-law authority.

## Runnable authored comparison

```bash
cargo run -p fs-material --example continuous_recovery
cargo run -p fs-material --example continuous_recovery -- 0.39071314
```

The example compares single, uniform and gamma populations at the same mean
rate `1/(6 h)` and `Rinf=0.9`, using one and eight tanks. Gamma shape is `0.5`
and rate scale is `2/(6 h)`. The default volume and flow give 24 h mean
residence. Gamma recovery is approximately **50.5636%** with one tank and
**58.8573%** with eight. Doubling flow halves residence and gives approximately
**40.8923%** and **48.5375%**. These are authored mathematical scenarios, not
measured ore or plant performance. The output retains original card hashes,
constituent flows, closure residuals and gamma integration diagnostics.

## Numerical diagnostics and limits

The gamma train evaluator integrates in log residence, where both response
integrands are smooth even when the rate density is singular at zero (`a<1`).
Log-concavity supplies tangent bounds on the omitted infinite tails. Adaptive
Simpson integration has at most 65,536 quadrature evaluations, depth 40, and
separately bounded mode and tail searches. Recovered and surviving probabilities
are integrated independently and checked for closure; the larger probability
then uses the exact complement of the smaller to retain tiny positive tails.

`ContinuousTreatment::gamma_integration` reports the numerical absolute-error
estimate, nominal omitted-tail bound, independent probability-closure residual
and quadrature work. These concern recoverable-population fractions before
`Rinf` scaling; they are not certified intervals, fitted-parameter uncertainty
or model uncertainty. Analytic single/uniform responses and explicit zero-rate
or zero-recoverable responses have no gamma integration report.

Missing sources, wrong units, unsupported conditions, invalid geometry,
stopped through-flow, exhausted work and unresolved or unrepresentable positive
responses refuse. No partial treatment is returned. This steady dilute-
constituent model supplies no reagent limitation, froth loss, entrainment,
adsorption, transient hold-up, solid/liquid slip, volume change, liberation or
plant scale-up calibration. Those effects require their own declared models.
