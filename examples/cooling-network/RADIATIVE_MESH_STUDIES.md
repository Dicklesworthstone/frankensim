# Mesh studies with radiative cooling

```bash
frankensim --json cooling-network \
  examples/cooling-network/adaptive-radiative-contact-hotspot.json
```

Steady `cooling-network` requests can now combine `radiation` with
`mesh_convergence`. Both the existing uniform ladder and the explicit
`strategy=goal-recovery` local strategy run the COMPLETE radiating producer.
This supersedes the mesh-study exclusion in the earlier radiation documents.
The example combines a 20 W localized source, two k(T) materials, matching
contact resistance and a fan-driven split/merge air network. Its physical inputs
are illustrative, not measured hardware. The 0.5 K mesh-comparison tolerance is
an explicit example policy, not an accuracy guarantee.

## Preserve the model while changing the mesh

Radiation remains the declared mean-temperature patch Robin closure described
in `RADIATIVE_COOLING.md`. Refining triangles does NOT split a named patch into
new independent radiators. Each patch retains its name, complete surface,
emissivity, surroundings temperature and provenance. Its mean is evaluated on
the refined trace. Pointwise integration of T(x)^4, enclosure radiosity and
occlusion are not introduced or approximated by this feature.

The original P1 volumetric source is prolonged without renormalizing its support.
Each child inherits its parent material law. Cooling faces retain their region
ownership, and matching contact traces remain separate unknowns with unchanged
resistance. The ordinary parser re-admits the complete refined request. Each
mesh enforces the existing nonlinear radiation, convective exchange, contact
and whole-system heat checks; radiative heat never enters an air branch.

## Goal marking includes the radiation feedback

With `strategy=goal-recovery`, the marking field comes from the existing TOTAL
radiating adjoint, including the radiation coefficient's dependence on the
patch mean, mixed-air transport and nonlinear material Jacobian. It is not the
adjoint of the convection-only problem or a frozen radiation coefficient.
Recovery scores choose cells; they are not Kelvin error estimates or DWR bounds.

`objective.gradient` remains false for a mesh study. Internal marking still
executes its required adjoint, but does not publish that work as requested
primal sensitivities: ordinary gradient fields and `radiation.adjoint` remain
null. The actual work remains visible in `total_adjoint_sweeps` and radiation's
reconstruction count. Uniform-only studies do not run a marker adjoint.

Two or more consecutive objective comparisons must pass. Local refinement must
also complete the existing global uniform-refinement confirmation. A failed
confirmation resets the streak; insufficient mesh, refinement, radiation,
derivative or wall-time budget refuses instead of publishing partial success.
Observed agreement is not a continuum error bound, maximum-temperature
certificate, experimentally validated model, or guarantee of monotone refinement.

## Results and independent replay

The final `radiation` report describes the final solved mesh. Per-mesh
`solid_solves` and the study's `total_solid_solves` include every completed inner
radiation FEM solve and any adjoint reconstruction solve, not just outer air
iterations. Adjoint sweeps remain a separate count, not every Krylov iteration.
The original wall deadline and numerical limits apply throughout the ladder.

`mesh_convergence.resolved_request` retains the final mesh, source, materials,
contact faces and unchanged radiation declaration. Submit that object to the
ordinary `cooling-network` command to reproduce the published field. Its
`gradient=false` replay does no marker work; work counts can therefore differ
without changing the physical temperature field. No recursive mesh study is
embedded in the resolved request. Transient or nested design mesh studies remain
excluded; this adds neither space-time adaptation nor mesh derivatives.

## Focused checks

```bash
cargo test -p fs-cli --test cooling_radiation_mesh
cargo test -p fs-cli --test cooling_adaptive_mesh
cargo test -p fs-cli --test cooling_mesh_convergence
cargo test -p fs-cli --test cooling_radiation
```

Six new actual-command Rust regressions cover both strategies, exact refined
field replay, preserved patch/material/source declarations, signed heating by
hot surroundings, patch-order determinism, vanishing emissivity, marker-free
uniform execution and budget/global-confirmation refusals. These tests were
written but NOT executed in the authoring environment, which lacks Rust/RCH.

Independent sparse NumPy/SciPy FEM and direct-transpose calculations, with
analytic air elimination rather than the Rust nested iteration, exercised the
same illustrative model. At the explicit 0.5 K tolerance, the uniform ladder
used 12, 96 and 768 cells, with peaks 329.419419, 329.133530 and 329.300586 K.
The local ladder used 12, 18, 32 and 256 cells, with peaks 329.419419, 329.485771,
329.609442 and 329.586494 K. Its last comparison was a complete uniform probe,
changing the peak by 0.022949 K. These are independent references, NOT Rust
executions, runtime speedups, or a proof that either reported peak is accurate.

On that final local mesh, omitting radiation gives 332.601378 K. A separate
350 K surroundings control loses -3.271406 W radiatively and sends 23.271406 W
to air, closing the 20 W source without clipping. Sixteen signed-source local
transfer checks preserved integrated power and both 0.01 m2 patch areas. A
vanishing-emissivity control differed from no radiation by 4.13e-10 K. No
universal cell reduction or physical validation is inferred from these checks.
