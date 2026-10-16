# Event-resolved heat transfer on prescribed separation histories

`FractureCooling::advance_prescribed_to` resolves the thermal consequences of
complete cohesive failure and directed opening/reclosure thresholds inside
one explicitly linear, fixed-local-frame jump segment. The ordinary
`advance_to` / `synchronize` endpoint coupling remains unchanged.

Supply one `prescribed::Endpoint` per bound interface, an absolute target
in seconds, and `prescribed::Config`. The retained jumps are the start.
`advance_prescribed_history` accepts strictly ordered `Frame` values, each
borrowing a complete endpoint array. Holds repeat the same jumps at a later
time; each actual reversal needs a frame. It publishes the whole history
atomically, not a partially successful prefix. Clone and resume at supplied
frame boundaries with the same numerical grid for deterministic replay.

## What changes physically

The actual cohesive law locates complete failure. Normal gap thresholds use
the same interpolated jumps as the material updates. Each numerical crossing
bracket supplies a cut in the existing thermal march. A failed open interface
restores ambient exposure; a failed closed interface uses the explicitly
supplied reclosure conductance. Reclosure never repairs cohesive strength.
Failure can change conductance without opening a surface, and an initially
intact interface can fail, open and reclose within one nonproportional segment.
Comparing only initial/final contact flags would miss that open interval.

An affine normal gap crosses each directed threshold at most once. The
effective opening radius is convex along a straight segment; starting inside
the failure surface, it has at most one outward crossing. These restricted
kinematics justify the search. They do not establish event completeness for
an arbitrary reconstructed equilibrium trajectory or an unreported reversal.

## Numerical and model boundaries

Bisection evaluates the physical owner's law without aging the accepted
history. Returned brackets are numerical locations, not interval-verified
physical certificates. Cuts use their upper endpoints; another simultaneous
cut can detect a crossing earlier inside its bracket. The declared event-time
tolerance controls this numerical delay separately from the thermal timestep.
A converged event location does not bound thermal discretization or model error.

Whole-history event-query, crossing-storage and thermal-substep allowances are
shared across all motion segments and event cuts. The existing thermal ceiling
can only be tightened by `advance_observed_limited`. Counts cover predicate
queries and thermal steps; they are not a total CPU-time or byte certificate.
Interface sweeps and cloned states are bounded by the existing interface/cell
budgets and retained crossing count. There is no full-field history log.

Unrelated thermal schedule events keep their ordering and exact times. They
execute before mechanical updates at a shared endpoint. Conflicting scheduled
changes to mechanically owned contacts still refuse. The whole thermal energy
balance is recomputed from final storage, incoming material and external heat;
fracture energy remains separately reported rather than silently injected as heat.

This is prescribed motion, not a structural equilibrium, stability or dynamic
snap-back solver. Contact geometry/areas/frames and cohesive material laws are
fixed. No thermal feedback into that prescribed trajectory, simultaneous
healing/fracture, evolving contact discovery or calibrated print prediction
is claimed. Inserting motion frames can change the thermal timestep grid.

## Example and validation status

```sh
cargo run -p fs-conduction --example prescribed_debonding
```

The synthetic two-cell example opens, recloses and reopens a bond, checks
permanent failure and thermal/fracture accounting, compares four thermal
step sizes against the independent two-body exponential solution, and checks
restart at an existing motion-frame boundary.

Eleven focused Rust regressions accompany the implementation. Neither they
nor the example have been run by Rust in the authoring environment.
DSR (including its documented
absolute fallback), RCH and Cargo were unavailable. Executed Python/SciPy
reference equations are not native Rust validation or material calibration.
