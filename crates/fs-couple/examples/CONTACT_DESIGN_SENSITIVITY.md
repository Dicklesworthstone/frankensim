# Contact-loaded inverse-design sensitivities

```bash
cargo run -p fs-couple --example contact_design_sensitivity
```

The example settles a physically supported 40-gram mass against a resonator under
2 N, then evaluates a target receiver displacement and its design derivatives.
It returns a dimensionless objective, receiver displacement in metres, gradients
with respect to physical actuator newtons, contact gap/stiffness and support
stiffness, and the recomputed adjoint residual. Its numbers describe an authored
reduced model, not an experimentally identified contact pair or an optimized
instrument. This is a derivative-producing consumer, not another optimizer.

## Native use

First call the existing `initialize_static_equilibrium` or
`initialize_contact_equilibrium` on `CoupledModalSystem`. Then construct
`coupled::equilibrium::sensitivity::EquilibriumLinearization` with that immutable
network, its actual external forces, the complete contact set, and explicit
preparation/query/margin budgets. A stale load, omitted reacting contact, moving
state, unresolved free support, or near-switch contact refuses. The source
network is immutably borrowed, so a prepared tangent cannot silently survive a
state or parameter mutation through these APIs.

The static equations are differentiated, not nonlinear iteration history or
simulation timesteps. Write their residual as

```text
R(q,p) = K0 q + sum_i b_i r_i(b_i^T q - gap_i) - f
```

with bilateral rest loads included in `K0`'s existing primal response. The
linearization applies and solves

```text
K = dR/dq = K0 + sum_i tangent_i b_i b_i^T.
K^T lambda = dJ/dq
 dJ/dp = explicit_dJ/dp - (dR/dp)^T lambda
```

`solve` is the tangent/adjoint solve and `parameter_pullback` is the **positive
residual pullback**, not the total gradient. This is the same IFT convention as
`fs-adjoint::ift`; no Krylov iteration or nonlinear sweep is differentiated.
The factorization stays in the existing `fs-la` owner. Its base inverse is the
existing elastic/supported-free static response, plus a contact-sized Woodbury
correction. No full modal matrix or finite-difference primal solves are used.
The complete coordinate residual is checked after each solve, not only the
small contact-space system.

Partials cover generalized held forces, ordinary modal angular frequencies,
bilateral stiffness/rest extension, contact stiffness/gap/quadrature weight,
and signed attachment columns. Static damping partials are zero; these goals
cannot identify damping. The exponent and uncertainty of a contact law are not
included in this parameter family. Explicit free coordinates report no angular-
frequency partial; differentiating a physical mass requires the corresponding
coordinate and attachment-normalization chain rule.

## Physical observations and force units

`displacement_objective` evaluates
`sum 0.5*weight*((observed_m-target_m)/scale_m)^2`. Every observation and additive
actuator uses an explicit `ModalAttachment`. Normalization scales and weights
are declared; they are not automatically estimated noise levels. One adjoint
serves all observations and physical actuator derivatives. The result includes
explicit target, scale, weight and observation-shape partials alongside the
mechanical residual pullback.

A shape used both as an observation and as a contact/connection parameter needs
**both** contributions. The observation term is direct; the mechanical term is
minus the residual pullback. For a signed connection column `B=left-right`, left
shape parameters take the column partial, while right shape parameters take its
negative. The regression tests deliberately perturb a shared right-side shape
to catch either sign errors or dropping the observation term.

Physical actuator gradients are per newton, not per mass-normalized force.
For a 0.04 kg unit translation the shape is 5 per sqrt(kg): projection into the
mechanical residual and projection of the adjoint back to physical force must
both use it. No scalar `1` is silently substituted.

## Admission and scope

The contact law owner computes the stationary force and analytic differential.
A caller-provided distance from contact switching is enforced. This is not a
certificate that a finite optimization step remains in the same activity region;
re-solve and re-admit at every candidate. Primal force tolerances and derivative
residuals are retained separately. These local floating-point results are not
interval-certified derivatives or material-validation evidence.

Preparation screens `n*(k+p+1)^2+(k+p+1)^3` terms, and ordinary queries screen
`n*(k+p+1)+(k+p+1)^2`, where n is retained coordinates, k bilateral connections,
and p contacts. Physical objective visits additionally count targets/actuators.
These bound admitted structural work, not exact flop counts or wall-clock time.
Cancellation and all failures return no partial gradient family and never alter
the underlying equilibrium. Repeated query count is owned by the caller.

The scope is a fixed, mass-normalized static model: not transient audio gradients,
frictional equilibrium, contact-only grounding of a free body, mode extraction,
CAD/material uncertainty, or a real-time guarantee. Dynamics and existing audio
rendering paths are unchanged.

```bash
cargo test -p fs-dcontact --lib opening::static_response
cargo test -p fs-couple --test equilibrium_sensitivity
cargo test -p fs-couple --test modal_contact_preload --test supported_free_preload
```

Native compilation and these Rust tests have not been executed in the authoring
environment, which lacks Cargo/rustc. Independent numerical references test the
derivation but do not execute Rust or the deterministic math implementation.
