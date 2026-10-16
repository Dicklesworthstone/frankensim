# Stress-constrained material design under prescribed motion

```sh
cargo run -p fs-topopt --features cutfem-marquee --release \
  --example motion_stress_sdf3 -- 40 250000 0.1
```

The fixed implicit slab is stretched and compressed by two independent
reference actuator laws, with independently solved body forces. Arguments are
accepted updates, cumulative Krylov allowance and the stress cap. The selected
functional is the existing unweighted sampled-peak norm of both relaxed and
physical stresses. It bounds retained numerical samples, not continuum stress.
No density-changing interpolation, new stress integrator or optimizer is added.

`evaluate_stress_with_prescribed_motion` accepts the existing force `LoadCase`
family and one `Option<&PrescribedMotion3>` per case. `None` keeps that case's
homogeneous support. Motions apply only to retained embedded support patches;
strong box clamps remain homogeneous. All laws are pure and fixed in the
reference configuration. Surface and material model changes require a new study.

The equation is K(s) u = f + b_g(s). The stress derivative includes the direct
stress-measure term plus z^T (db_g/ds - dK/ds u), followed by the original
SIMP/projection/filter pullback. Freezing the imposed-motion load gives a wrong
derivative even though the forward field still solves its current equation.
Preparation is shared across all primal and adjoint solves at each density.
Successful probes and returned errors restore incoming material. Original
all-None calls retain the legacy arithmetic and counted work.

`StressDesignStudy3::new_with_prescribed_motion` runs the same projected AL,
stress measure, bounds and feasible archive. Split `run` calls retain the same
search state. `restore_with_prescribed_motion` rebuilds saved accepted and
incumbent endpoints under the original motion laws before restoring the existing
checkpoint. Consumers must bind those laws, geometry and budgets externally;
matching endpoint values is not proof of arbitrary callback equivalence.

`case_compliances` remains external force work f^T u. It excludes b_g^T u:
the latter is not a physical actuator-work law. The example reports external
work, all per-case maxima, density fields and the least-volume feasible accepted
endpoint before propagating late failures. Its optimizer stop/KKT residual
belongs to the latest AL iterate, not necessarily the selected feasible point.

Pure displacement-controlled stress can decrease when material is removed.
A stress cap alone does not enforce actuator force, stiffness, connectivity or
unique material recovery. This extension does not combine the response-limit
and stress constraints into one optimizer, nor add follower loads, moving-boundary
shape derivatives or stress-specific adaptive certification. The existing
ersatz-turnover admission remains unchanged, including for mixed loading.

Nine native regressions cover affine stress, the full mixed-load gradient,
zero-weight peak inclusion, legacy equivalence, rollback, input admission,
minimum-volume descent and restoration. They have not run in the editing
environment: Cargo is absent. An executed independent 1-D Nitsche rod oracle
checks the derivative equation, not the Rust/3-D implementation.
