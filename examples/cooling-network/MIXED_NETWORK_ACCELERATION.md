# Accelerated mixed-network cooling

`cooling-network` now uses bounded vector IQN-ILS for the complete solid/air
reference-temperature interface, including upstream heating, bypasses and
split/merge mixing. This closes the gap between the earlier independent-branch
accelerator and the actual mixed-network producer used by this command.

```bash
frankensim --json cooling-network examples/cooling-network/stiff-mixed-contact.json
```

The example retains a split solid, finite-resistance matching contact, hot
supply, cold bypass and downstream mixing. It deliberately sets strong film
coupling and only twelve primal/adjoint sweeps. Its coefficients and contact
are illustrative solver inputs, not measured hardware or validated data.
The independent slab/NTU equations give a first-wall mean of approximately
328.4882733363 K; the actual-command regression checks this reference.

## One physical map, bounded vector history

A map evaluation performs one shared solid solve and the existing graph air
march. One history spans all interface regions. The inverse least-squares
update retains individual residual components instead of collapsing opposing
components into a signed branch average. The existing `fs-couple` producer
uses twice-orthogonalized, rank-filtered QR and at most eight retained secants.
No dense network Jacobian or normal-equation solve is introduced.

The command uses rank tolerance 1e-10. Its existing `tolerances.relaxation`
sets startup and rank-zero fallback. If an extrapolation proposes nonpositive
absolute reference temperatures, the history is discarded and the declared
relaxation is used instead. Coordinates are not clipped. Non-finite arithmetic
or an actual solid/material/transport refusal remains an error.

A proposed update is NEVER itself evidence of convergence. The next fresh
solid and air evaluations must pass the existing unrelaxed temperature and
independent branch-local watt checks. Returned references are the exact ones
used by the final solid solve, not the next air proposal. The separate
whole-domain source/Robin/contact accounting is unchanged.

## Tangents, adjoints and transient endpoints

The explicit accelerated library methods also solve the forward and transpose
implicit interface equations. They reuse the real FEM and transport derivative
producers, including K'(T), direct heat-objective terms and source-load
pullbacks. They do not differentiate the primal iteration history. Each
individual tangent/adjoint solve starts its own fresh secant history and must
satisfy the original equation residual before returning a derivative.

The cooling command uses the accelerated adjoint. Its existing contact
bilinear-form projection therefore receives the same total coupled gradient.
Effective-h design searches inherit this path; fan trials still recompute
hydraulics and correlation-derived coefficients as before.

Transient endpoints inherit the same accelerated map. Every map evaluation
still holds the previous accepted physical temperature field fixed. Every new
timestep, adaptive trial, workload or fan candidate starts NEW history. Neither
a rejected coupling proposal nor a discarded adaptive trial advances physical
history. This addition does not implement transient derivatives.

Results include `coupling_solver` with method, scope, history/rank policy,
fallback factor and the adjoint method (null when no adjoint was requested).
The request's coupling, derivative, linear and wall-time budgets are not raised.
Iteration counts and final floating-point rounding may differ from stationary
runs. Existing UQ executable-bound checkpoints continue to reject a different
binary; old samples are never silently relabeled as new-executable results.

## Library compatibility and checks

`solve_coupled_transport` / `_from` and `CoupledLinearization::apply` /
`pullback` retain their stationary behavior. The explicit new alternatives are
`solve_coupled_transport_iqn` / `_iqn_from` and `apply_iqn` / `pullback_iqn`,
which accept `IqnIlsConfig`. A reference-only accelerated warm start drops
secant history: it is continuation, not bitwise replay of an interrupted solve.

```bash
cargo test -p fs-airflow --test mixed_network_iqn
cargo test -p fs-airflow --test mixed_network_iqn_adjoint
cargo test -p fs-cli --test cooling_mixed_iqn
```

The ten new tests cover stiff split/merge coupling, real FEM/air derivatives,
full transpose identities, contact gradients, nonlinear transient regression,
positive-temperature fallback, cancellation and exhausted budgets. The Rust
tests were not executed in the authoring environment, which has no Rust
compiler. Independent Python thermal-algebra checks are not Rust benchmarks:
they needed five map evaluations for the powered split/merge fixture and four
for each slab primal, tangent and adjoint; stationary controls exhausted twelve.
No general convergence, continuous-time peak, continuum-error or physical
validation claim follows from these numerical solver changes.
