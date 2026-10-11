# Adaptive radiation on the native `.fsim` product path

A steady native project with declared ambient surface radiation can select
`adaptive` solver fidelity with prescribed convection or independent,
fixed-flow airflow-convection branches. This is the same import → solve →
QoI → retained-report path used by the ordinary project, not a separate
cooling-network laboratory command.

The enriched comparison keeps each patch's immutable emissivity card, query,
validity domain and reservoir. It evaluates the radiative secant at each
field's own area-mean wall temperature. With airflow, both the solved field
and the prolonged coarse field are independently marched through the original
convective air paths. Only convection exchanges heat with air; radiative heat
remains assigned to its black reservoir, including when a hotter reservoir
heats the solid.

The existing conduction owner supplies the complete transposed tangent,
including nonlinear conductivity, matching contact, radiative coefficient
and weighted-reference derivatives, and cross-wall air-reference feedback.
There is no frozen-Robin fallback and no replacement primal hidden inside the
goal comparison. The physical reference residual and dual residual must both
meet the configured linear tolerance. A field adequate for a visual report
can therefore refuse a tighter adaptive comparison; tighten the physical
radiation convergence controls rather than removing the residual check.

## Reading the retained evidence

Each `adaptive.history` row retains the signed dual-weighted change, finite
nonlinear remainder, and maximum-location remainder. Their signed sum is the
observed maximum-temperature change; their absolute sum is the marking and
stopping estimate. Contributions are assigned to parent tetrahedra using the
uniform enrichment lineage before the existing marked-refinement step.

The row's `algebraic_balance.comparison_model` is
`area-mean-radiation-full-transpose`; `air_reference_feedback` records whether
air paths participate. The legacy `coupled_adjoint` field is not reused: that
field specifically describes IQN-ILS, whereas this complete transpose uses
the conduction owner's bounded FGMRES path.

These are Estimated, same-model discrete comparisons. They are not certified
continuum error bounds or physical validation. The nonlinear model's rigorous
algebraic balance remains unsupported, rather than becoming a zero bound.
Unmeasured model-form uncertainty remains unmeasured. At most 64 selected
uniform Robin traces are admitted, and four nodal feedback factors per trace
plus the reference matrix must fit the declared feedback-entry budget. This
is not an allocator/RSS certificate. Cancellation and execution deadlines
remain active.

Combined natural-convection/radiation adaptive goals remain explicitly
refused; use base or ladder fidelity for that combination. This change does
not add enclosure radiosity, occlusion, temperature-dependent emissivity,
pointwise fourth-power quadrature, recirculation, or flow redistribution.

## Native regression commands

```sh
cargo test -p fs-cli --test solve radiation_product
cargo test -p fs-conduction mean_feedback_goal::radiation
```

The product regressions cover the independently balanced slab, a hotter
reservoir, convection-only air heat, repeated retained outputs, nonlinear
material/contact coupling, and missing-emissivity refusal under both base
and adaptive fidelity. The core tests separately exercise the total tangent
and own-state residual comparison.
