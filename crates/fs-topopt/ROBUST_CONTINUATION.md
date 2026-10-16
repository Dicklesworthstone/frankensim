# Robust 3-D design continuation and eroded error assessment

Run the curved implicit cantilever with three material stages, independent
loads, and a final refined-grid assessment:

```sh
cargo run --release -p fs-topopt --features cutfem-marquee \
  --example elastic_sdf3 -- --robust-continuation 200 250000
```

The positional arguments are accepted updates **per stage** and one cumulative
Krylov-iteration allowance. The schedule is `(penal,beta)=(1,1),(2,2),(3,8)`.
The example fixes its compliance scale from a charged physical baseline and
prints the scale, threshold offset, filter radius, raw floor, material cap,
optimizer policy, geometry limits and spent work. Stage completion is not a
convergence or optimality claim.

## Stage boundaries carry raw controls, not old physics

`robust::continuation::controlled_sdf3_robust_continuation` restores the
**new dilated** volume cap before checking eroded-compliance and dilated-volume
sensitivities with independently solved numerical probes. Each admitted model
gets a fresh projected augmented-Lagrangian state. The next stage receives only
the lowest-compliance feasible accepted raw design from that model; its old
multipliers, gradients and displacement fields are not transferred.

`RobustContinuationReport3::last` is the feasible incumbent of the **latest
admitted model**, with matching `params`. It is not a minimum across different
models. `selected_iteration` identifies its row in the stage history, while
`outcome` describes the final optimizer iterate, which may be different.

Failed transition preparation preserves the previous admitted result. A late
interruption after an accepted optimizer step retains the aligned stage prefix.
Gradient-check rejection, restoration failure, optimizer resource limits and
schedule completion have distinct stops. All solves and probes share the
caller's physical budget; optimizer evaluation limits apply separately to each
new stage. `controlled_sdf3_robust_feasible_start` also exposes the dilated-cap
restoration without an elasticity solve.

## Assess the eroded material, not the borrowed nominal state

Continuation and robust evaluation leave the borrowed study's incoming nominal
parameters and stiffness scales unchanged. Do not pass its nominal operator
and an eroded displacement field to the ordinary goal estimator.

For an adaptive background, call
`estimate_robust_reference_compliance_enrichment` with the retained raw design,
its explicit nominal `params`, robustness, reference load laws, and a compatible
enriched operator. It freshly solves eroded coarse fields, temporarily installs
their eroded stiffness scales, and calls the existing reference-load DWR path.
The fine operator inherits **physical eroded stiffness** rather than refiltered
fine raw controls. Success and returned errors restore the incoming study.
The reference load laws, geometry, base material and homogeneous supports must
represent the same physical problem on both grids.

The example prints stage histories, selected three-field densities, independent
eroded displacement fields, the coarse/fine compliance difference, and bounded
Dörfler marks. It exports the feasible prefix before returning a late failure.
The refined solve is a diagnostic probe; it does **not** reoptimize a fine-grid
design. Neither its correction nor its marking mass is a continuum bound, a
minimum-feature-size guarantee, or an estimate of dilated-volume error.

Focused native checks (not executed in the connector editing environment):

```sh
cargo test --release -p fs-topopt --features cutfem-marquee --test sdf3
cargo test --release -p fs-topopt --features cutfem-marquee --lib sdf3::robust::goal::tests
cargo check --release -p fs-topopt --features cutfem-marquee --example elastic_sdf3
```
