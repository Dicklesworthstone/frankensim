# Derivative-free fitting of supplied thermal measurements

`spatial_measurement_fit` can use the existing feasibility-ranked differential
 evolution instead of SQP. This is useful when an initial state lies exactly on a
phase-chart corner or when cold-start liquid fractions are flat across small
power changes. The original physical solver and all declared requirements stay
in force. The search requests values only, with no tangent tape or transpose
solve. It does not replace physical failures with large objective penalties.

From the repository root:

```bash
cargo run --manifest-path examples/scan_power_fit/Cargo.toml \
  --example spatial_measurement_fit -- \
  examples/scan_power_fit/examples/specimens/multiple_measurements.csv \
  --specimen examples/scan_power_fit/examples/specimens/rectangular.scan \
  --heaters examples/scan_power_fit/examples/specimens/multiple.heaters \
  --population 40 2400 --seed 17 --objective-target 1e-8 --substeps 4
```

`--population MAX_W ATTEMPTS` declares a finite physical search box
`0 <= P <= MAX_W` for every **selected** power. Original reference-watt scales
are used to convert that box to optimizer coordinates. Fixed powers (including
complete fixed heaters) are unchanged and are not implicitly limited by this
box. Initial selected powers outside the box refuse; they are not clipped.
ATTEMPTS must be 1..100000. The original DE policy uses F=0.8, CR=0.9 and
max(24, 4*selected-power-count) members. The seed defaults to 17.

`--objective-target LOSS` defaults to 1e-8. It bounds the actual sum of declared
quadratic or Huber losses, not each raw residual and not a physical uncertainty.
`--seed` and `--objective-target` without `--population` refuse. Existing
commands retain their SQP behavior. Instantaneous/spatial observations,
exposure means, first-order sensor response and custom specimens/heaters all
use the same input and observation owners as before.

Add `--requirements FILE` for native power, absorbed-energy, sampled state/peak,
slew and speed limits. Constraints are evaluated during selection, not after
an unconstrained search. The existing bounded requirement-input limits remain
in force, including their conservative row/work admission; population execution
allocates no gradient arrays or dense KKT matrix. Fixed paths remain fixed, so
power search cannot repair overspeed motion. Limits are not inferred from data.

The physical work allowance is ATTEMPTS plus **one reserved final replay**.
Failed physical evaluations propagate. After a successful search advance the
best-ranked sources are replayed against every original requirement, and the
objective and feasibility score must reproduce bit-for-bit. Export separately
checks every fixed power, path, clock, interpolation and shutdown. Search attempts
and physical evaluations are reported separately. The replay is a fresh solve,
not a cached fitness result.

The process exits successfully only if the replay is feasible within the original
scaled tolerance 1e-8 **and** meets the declared objective target. On budget
exhaustion it still prints the replayed best-so-far, then exits unsuccessfully
unless that target is met. A feasible but mismatched result is not target recovery.
No stationarity, infeasibility proof, unique recovery or global-optimum claim is
attached to this population search. Unsampled peaks and experimental validity
remain outside the admitted model.

A uniform one-kilogram fixture starting at h=80 J/kg, with solidus h=100 J/kg
and latent heat 100 J/kg, needs 30 W then 20 W to reach fractions 0.1 and 0.3
at 1 and 2 seconds. At the solidus itself it needs 10 W then 20 W. With a
20 J total-energy cap on the latter fixture, final fraction cannot exceed 0.2:
matching the second target is impossible, even though a feasible fit exists.
The focused Rust regressions exercise these cases through the actual model;
native execution of the new tests has not been verified in this session.

## Parallel candidate evaluations

Add `--workers N` (1..128) to opt into the existing TilePool path. The default
remains serial, and the flag without `--population` refuses. Each worker runs
its own complete original thermal trajectory from its own execution context;
the shared immutable specimen, observations, source program and requirements
are not reconstructed differently per worker. There is no global solver mutex.

The original physical study owns **one** evaluation budget. Workers draw from
that shared allowance, not ATTEMPTS each. The executor drains launched work
before returning, and failed physical evaluations propagate without a best-fit
success report. The one final replay remains serial and reserved. No new
checkpoint CLI or retry policy is introduced by this mode.

Batches contain at most 24 candidates. Candidate generation and publication
use the original logical ordering and frozen parent population, not completion
order. Target checks occur between batches; a parallel run can finish extra
already-launched evaluations after the target first becomes attainable and can
therefore return a different best result than an early-stopped serial run.
Fixed-budget equivalence is covered by focused native regressions, which remain
unexecuted here. No native parity or speedup result is claimed.

The 64 MiB executor lease covers dispatch bookkeeping, **not** the aggregate
memory used by simultaneous physical solves. Size the worker count for those
solver workspaces. Reported worker count is requested parallelism, not a claim
that all requested workers are busy in every wave.

An explicit cold-start example, beginning with both powers off:

```bash
cargo run --manifest-path examples/scan_power_fit/Cargo.toml \
  --example spatial_measurement_fit -- \
  examples/scan_power_fit/examples/specimens/cold_measurements.csv 0 0 \
  --specimen examples/scan_power_fit/examples/specimens/cold.scan \
  --population 40 2400 --seed 17 --objective-target 1e-8 --workers 4
```

The expected uniform-specimen schedule is 30 W then 20 W. These input files
make the phase-transition workflow callable without editing Rust source; they
are manufactured data, not an experimental calibration dataset.
