# Minimum material under explicit eroded-load compliance caps

```sh
cargo run --release -p fs-topopt --features cutfem-marquee \
  --example worst_load_sdf3 -- --minimum-material 1000 1500 200 250000
```

The first two numbers are authored physical compliance caps for the example's
independent Y and Z load laws, in its dimensionless units. They are illustrative
requirements, not calibrated ratings. The last two numbers are the accepted
update allowance and one cumulative Krylov-iteration budget. This mode never
computes or loosens a cap from the initial design. The existing worst-load
command without `--minimum-material` retains its original objective and behavior.

## Actual objective and constraints

`envelope::limits::ComplianceLimitedDesign3` minimizes **dilated material-volume
fraction**, using the existing projected augmented-Lagrangian solver. Every
independent load is solved on the **eroded material**. For caps `b_i`, its smooth
inequality is the existing unnormalized log-sum-exp envelope of `C_i/b_i`, minus
one. In exact arithmetic the envelope bounds the maximum from above; the extra
`temperature*ln(case_count)` is conservatism, not a violation allowance. All caps,
load ordering, temperature, geometry and material parameters remain fixed.

`ComplianceCaps3` distinguishes requirements from references that merely scale
an objective. Feasible-incumbent retention compares each physical `C_i <= b_i`
directly, without division, and also checks the smooth constraint at the AL
numerical tolerance. A good average, a small normalized violation, or a rounded
ratio cannot waive a physical cap. Strong smoothing can conservatively refuse
an otherwise physically admissible design; it is not silently reduced.

The best incumbent has the smallest dilated volume among accepted cap-feasible
points, not the lowest envelope. Accepted inner steps can be infeasible. Terminal
KKT residuals and stopping status belong to the final AL iterate, which may not
be the incumbent. A convergence status alone is not a design-feasibility claim.
Both constraint and objective gradients are checked at initialization using
independently solved numerical probes of the actual design map.

## Feasible starts and final re-solves

`limits::start::controlled_compliance_feasible_start` first tests the incoming
raw design. When needed, it tests raw one (without overriding protected voids),
then bisects the segment from the incoming point toward that endpoint. Every
probe uses actual eroded physics and preserves the caller's nominal stiffness.
The default is at most 16 bisections plus two endpoint solves, with a 0.001
reserve below the smooth boundary. The reserve does not alter the physical caps.
An already feasible start requires only its initial evaluation.

A rejected full-density endpoint is NOT a global infeasibility proof. The report
separates that refusal from numerical/cancellation stops. Late interruption
retains a previously solved cap-feasible endpoint when one exists; that seed is
not mislabeled as an optimized solution. The example emits it before returning
an initialization failure. All restoration, gradient-probe, rejected-trial and
final-recheck work remains on the same physical budget.

`recheck_best` freshly solves the retained raw design and checks every cap again.
It cannot mutate the optimizer trajectory or cached incumbent. Its work is
reported in the physical solve counters rather than optimization callback
counts. This is a fresh numerical solve, not an independent solver or proof.
The example exports the solved feasible prefix before attempting that recheck,
then reports the recheck outcome and returns nonzero on a late failure.

## Output and scope

Output includes restored-start diagnostics, accepted material/constraint history,
every load's compliance and cap, all three density fields and both independent
eroded displacement fields. The model is fixed-grid linear elasticity with
homogeneous supports and density-independent dead loads. It does not establish
continuum load guarantees, buckling or stress safety, manufactured length scales,
global optimality, or cross-process restart. `run` supports in-memory segments
without resetting AL state or budgets.

```sh
cargo test --release -p fs-topopt --features cutfem-marquee \
  --lib sdf3::robust::envelope::limits
cargo check --release -p fs-topopt --features cutfem-marquee --example worst_load_sdf3
```

Native Rust execution was unavailable in the editing environment. Independent
dense body-fitted 3-D P1 FEM/SciPy checks exercise the formulation, not these
Rust implementations or production CutFEM.
