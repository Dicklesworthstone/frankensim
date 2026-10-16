# Adaptive radiation on the native `.fsim` product path

A steady native project with declared ambient surface radiation can select
`adaptive` solver fidelity with prescribed convection, Churchill-Chu
vertical-plate natural convection, or independent fixed-flow airflow-convection
branches. This is the same import → solve → QoI → retained-report path used by
the ordinary project, not a separate cooling-network laboratory command.

The enriched comparison keeps each patch's immutable emissivity card, query,
validity domain and reservoir. It evaluates the radiative secant at each
field's own area-mean wall temperature. With airflow, both the solved field
and the prolonged coarse field are independently marched through the original
convective air paths. Only convection exchanges heat with air; radiative heat
remains assigned to its black reservoir, including when a hotter reservoir
heats the solid.

With natural convection, every declared natural law is evaluated at BOTH
fields' own wall means under the same envelope pressure and coefficient
multiplier as the primal. The reference tangent includes the signed heated-
or cooled-wall derivative of the actual Churchill-Chu card, including the
film-temperature dependence of air density and expansion. Ordinary prescribed
convection rows and independent radiative reservoirs keep their own values.
Natural laws on nonradiating surfaces still participate in the same full
transpose. The approximate field can cross from the heated to the cooled
branch; its own residual and the finite nonlinear remainder are retained.
Exact equilibrium and out-of-domain Rayleigh states still refuse.

The existing conduction owner supplies the complete transposed tangent,
including nonlinear conductivity, matching contact, natural and radiative
coefficient derivatives, weighted-reference derivatives, and cross-wall
air-reference feedback when present. There is no frozen-Robin fallback and no
replacement primal hidden inside the goal comparison. The physical reference
residual and dual residual must both meet the configured linear tolerance. A
field adequate for a visual report can therefore refuse a tighter adaptive
comparison; the acceptance gates are not relaxed to admit that field.

## Reading the retained evidence

Each `adaptive.history` row retains the signed dual-weighted change, finite
nonlinear remainder, and maximum-location remainder. Their signed sum is the
observed maximum-temperature change; their absolute sum is the stopping
estimate. Absolute nodal contributions are assigned to parent tetrahedra using
the uniform enrichment lineage for the existing marked-refinement step.

The row's `algebraic_balance.comparison_model` is
`area-mean-radiation-full-transpose`; `air_reference_feedback` records whether
air paths participate, and `natural_coefficient_feedback` in the history row
records the presence of natural-convection laws. The legacy `coupled_adjoint`
field is not reused: that field specifically describes IQN-ILS, whereas this
complete transpose uses the conduction owner's bounded FGMRES path.

These are Estimated, same-model discrete comparisons. They are not certified
continuum error bounds or physical validation. The nonlinear model's rigorous
algebraic balance remains unsupported, rather than becoming a zero bound.
Unmeasured model-form uncertainty remains unmeasured. At most 64 selected
uniform Robin traces are admitted, and four nodal feedback factors per trace
plus the reference matrix must fit the declared feedback-entry budget. This
is not an allocator/RSS certificate. Cancellation and execution deadlines
remain active.

Simultaneous natural and forced-airflow convection remains refused; no mixed
convection model is inferred. This change does not add enclosure radiosity,
occlusion, temperature-dependent emissivity, pointwise fourth-power quadrature,
recirculation, or flow redistribution. Previously accepted nonradiating and
prescribed/airflow-radiating paths retain their numerical and receipt semantics;
the newly admitted natural/radiating adaptive combination had no completed
checkpoint under the old driver.

## Native regression commands

```sh
cargo test -p fs-cli --lib natural_radiation
cargo test -p fs-cli --test solve radiation_product
cargo test -p fs-conduction mean_feedback_goal::radiation
```

The law regressions cover both wall branches, independent comparison-state
coefficients, pressure and coefficient scaling, permuted trace ordering, and
atomic refusal for missing/duplicate laws, changed ambients and invalid states.
The product regressions cover the independently balanced slab, a hotter
reservoir, convection-only air heat, repeated retained outputs, nonlinear
material/contact coupling, and missing-emissivity refusal under both base
and adaptive fidelity. The core tests separately exercise the total tangent
and own-state residual comparison.
