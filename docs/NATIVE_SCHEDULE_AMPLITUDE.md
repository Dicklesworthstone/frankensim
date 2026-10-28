# Native transient workload amplitude

The `power-schedule-scale` uncertainty target changes the amplitude of one
named volume region's original `power-schedules` declaration. Its coordinate
is a finite, nonnegative dimensionless multiplier, not watts or duty:

```lisp
(uniform :name "pulse" :target power-schedule-scale :entity "metal"
  :low 0.5 :high 1.5)
```

Every sample scales all delivered-watt entries of that original schedule.
Switch times, off intervals, source attribution, static power and duty,
other regions, initial state and solver controls remain unchanged. Since a
schedule already supplies delivered power, duty is not applied again.
For the example's 20 W, 0.25 s startup pulse, the support gives 10..30 W
and 2.5..7.5 J without changing the pulse duration or subsequent off period.
Each candidate starts from the immutable base project, never the previous
sample. Zero disables that region's scheduled heating without disabling its
storage, conduction or boundary exchange. An overflowing or underflowing
positive pulse refuses rather than being clipped or silently erased.

The target requires exactly one existing schedule on a volumetric conduction
region. A missing transient, missing/ambiguous schedule or simultaneous static
`power` coordinate for that same overridden region refuses. An invalid base
cannot be repaired by sampling a zero amplitude. Other model/domain refusals
remain those of the ordinary native physical solve.

```bash
cargo run -p fs-cli --bin frankensim -- --json study examples/native-uncertainty/pulse-amplitude.fsim pulse-amplitude.db
```

The ordinary Monte Carlo, replicated QMC/LHS, independent-input Sobol and
Bernoulli-compliance paths consume the same physical coordinate. Final and
sampled-peak observations retain their original meanings. Coordinate-secant and nominal-adjoint
mean controls use the same amplitude coordinate. Retained source-independent
resume reconstructs the same scaled schedule.

The example's probability laws, capacity and fan data are illustrative.
Schedule scaling does not assert that real workload statistics are uniform,
that temperature is monotone in amplitude, or that a sampled peak bounds
continuous-time temperature. No material/source authority is added.

## Size the pulse against the full transient

The same coordinate is an explicit native design control:

```lisp
:control (power-schedule-scale :entity "metal" :minimum 0 :maximum 3 :tolerance 0.001)
:temperature-tolerance 0.01K
:assumption nondecreasing-temperature-with-schedule-scale
```

The native sizing driver checks the minimum-amplitude baseline, then the
upper bound, and retains a measured feasible endpoint while bisecting the
bracket. It maximizes amplitude, not static authored watts or pulse duration.
The original schedule must contain nonzero heating. Even zero amplitude
cannot cure an initial state already above the temperature limit.

Every candidate runs the ordinary complete native trajectory. Its decision
uses the maximum over the original initial state and accepted fine-grid
states; the cooled final temperature is reported separately. Optional finite
`operating-cases` evaluate that entire trajectory under every authored case.
Only a completed group can establish feasibility. Cases can prescribe another
region's schedule amplitude, including zero; they cannot overwrite the search
coordinate. These finite conditions are not a probability law or an interval
robustness guarantee.

```bash
cargo run -p fs-cli --bin frankensim -- --json study examples/native-uncertainty/pulse-design.fsim pulse-design.db --budget 2
```

The ordinary `study --resume`, `report` and `package` routes retain the original
schedule and all completed physical observations. Budgets count actual native
case evaluations. Reports use `power_schedule_scale`, `schedule_region`,
`schedule_scale_tolerance` and `bracket_width_scale`; they never label this
multiplier as watts or fan speed. The selected exported project contains the
actually evaluated scaled watt history at its governing operating case.
Completion requires both coordinate and temperature tolerances, or a passing
upper bound. Measured nonmonotonicity refuses; global monotonicity remains the
explicit caller assumption.

## One-probe trajectory adjoint calibration

Version 3 admits `:mean-control (nominal-adjoint :max-solves 1)` for this
coordinate with both capacity and total-enthalpy storage, and final or selected
sampled-peak observations. It reuses the existing complete accepted-trajectory
adjoint, including every admitted material, contact and boundary feedback.
The amplitude derivative contracts each interval's absolute delivered-watt
partial with that interval's **original** scheduled watts:

`dJ/da = sum_i (dJ/dP_i) * P_original,i`.

Neither the timestep nor static duty is multiplied again. No division by the
sampled amplitude occurs, so zero amplitude does not erase the underlying
workload or turn a nonzero physical response into zero. The retained original
and sampled source, region, event times, scaled watts and complete native fine
grid must match; future intervals must have zero influence on an earlier
selected peak. Off periods contribute zero without borrowing another region's
power. The contraction makes no extra physical or transpose solve.

The existing one calibration solve supplies all supported controls before
probability samples start. Its frozen coefficients and original workload are
recovered from retained evidence, without the source files. Raw physical
samples and their numerical pass events remain unchanged. A selected maximum
branch is a local discrete derivative; ties, a switching active peak, numerical
error and physical-model uncertainty retain their existing limitations.

```bash
cargo run -p fs-cli --bin frankensim -- --json study examples/native-uncertainty/pulse-amplitude-adjoint.fsim pulse-adjoint.db
```
