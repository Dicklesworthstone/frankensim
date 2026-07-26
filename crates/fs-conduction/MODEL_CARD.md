# Model card: steady heat conduction (`fs-conduction`)

Fields mirror `fs_evidence::ModelCard` so this document and a runtime card
cannot drift apart. Beads `frankensim-extreal-program-f85xj.5.1`, `.5.3`, and
`.5.4`.

| field | value |
| --- | --- |
| `name` | `steady-conduction-p1-feec` |
| `version` | `0.0.1` (the crate version; `fs_conduction::VERSION`) |
| `ambition` | `[S]` Solid — established mathematics, default path, no feature flag |
| `discrepancy_rel` | **not claimed**. See "Discrepancy" below: this crate reports residuals and observed orders, not a model-form band |
| `calibration` | none. No parameter of this model is fitted to data; the only empirical inputs are `fs-matdb` conductivity claims, and each carries its own receipt and uncertainty |

## What is modelled

Fourier conduction in a rigid solid at steady state:

```text
  ∇·(k(T, x) ∇T) + f = 0
  q = −k(T, x) ∇T
```

with `k` a symmetric positive-definite second-order tensor that may depend on
temperature, and `f` a volumetric source in W/m³.

## Assumptions

Every one of these is a modelling choice the solve inherits, not a numerical
detail:

1. **Steady state.** `∂(ρ c_p T)/∂t = 0`. Any transient is out of model. No
   heat capacity or density appears anywhere in the crate — their absence is
   structural, not an omission to be patched at the call site.
2. **Rigid, non-deforming solid.** No thermal expansion feedback, no
   deformation-dependent conductivity, no advective transport term `ρ c_p u·∇T`.
   The domain is fixed.
3. **Fourier's law with a symmetric positive-definite tensor.** No hyperbolic
   (Cattaneo) conduction, no size effects, no ballistic transport. Symmetry and
   positive definiteness are CHECKED at construction (Sylvester's criterion on
   the leading principal minors), not assumed.
4. **Local thermal equilibrium**, a single temperature field. No two-temperature
   or porous-medium formulations.
5. **`k` depends on temperature only** — never on position independently, on
   gradient, on history, or on damage state. A spatially varying material is
   expressed by solving per-region or by the per-element multiplier that
   `assemble_operator_scaled` provides.
6. **Base boundary data is exact and known.** `T_D`, `q_n`, `h`, and `T_ref`
   are inputs. The optional radiation path computes a card-backed linearized
   coefficient or a diffuse-gray surface exchange from an admitted view-factor
   matrix; it does not generate that matrix or compute a convective
   correlation.
7. **The Robin reference temperature is a declared property of the row.**
   `fs-scenario`'s Robin row carries only `h`, so `T_ref` is named at the
   lowering call. A correlation that supplies `h` must state the ambient it was
   fitted against; this crate will not infer one.
8. **Contact is never implicit.** A coincident duplicated P1 trace must be
   bound to one positive card-backed area-specific resistance; otherwise the
   contact solve refuses. Perfect contact is not inferred.
9. **Gray-diffuse radiation is surface-only.** Each named trace has one
   hemispherical-total emissivity, surfaces are opaque and diffuse, and the
   enclosure matrix is fixed. Participating media, spectral/specular effects,
   and geometry-derived visibility are outside this model.

## Validity domain

| axis | bound | set by |
| --- | --- | --- |
| `T` | the sampled conductivity span `[low, high]` K | `ConductivityTable::from_claims`, from the declared grid; the intersection over components is `ConductivityModel::temperature_span()` |
| radiation `T` | the emissivity card's validity domain; for a linearized row, additionally `|T_s - T_mean| <= max_departure` | `SurfaceEmissivity::from_card` and `LinearizedSurfaceRadiation` |
| view factors | each row closes within its admitted absolute tolerance and `A_i F_ij = A_j F_ji` within its admitted relative tolerance | `ViewFactorMatrix::admit` |
| everything else | unconstrained | no other axis is constrained by this crate |

Outside the temperature span the model REFUSES
(`ConductionError::OutsideTemperatureSpan`) rather than extrapolating — during
assembly, during flux recovery, and inside the line search, where an
inadmissible trial iterate is a rejected step rather than an error.

For an `fs-matdb`-backed model the span is additionally bounded by the source
claims' own `ValidityDomain`: a grid point outside it never produces a knot,
because the query refuses first. The retained receipts record exactly which
claims were considered, which was selected under which policy, and how each
value was evaluated.

## Known failure modes

1. **Pure Neumann / no Robin.** With neither a Dirichlet nor a Robin row the
   steady operator is singular up to an additive constant. A Krylov method
   returns something that looks like a temperature field; this crate refuses
   (`SingularPureNeumann`).
2. **A `k(T)` strong enough to break Newton.** The Armijo line search backtracks
   and, if it cannot find a step, refuses with `LineSearchFailed` naming the
   iteration and the smallest step tried. It does not return the last iterate as
   if it had converged.
3. **An iterate leaving the sampled span.** Common when a source is large
   relative to the sampled range. The guard fires as a rejected step; if the
   converged solution genuinely lies outside the span, the solve refuses. Widen
   the grid — do not widen the claim.
4. **Slivers and degenerate elements.** Refused up front
   (`DegenerateElement`) against a relative volume floor, before `fs-feec`'s
   own assertion can fire.
5. **Faceted curved boundaries.** A curved surface meshed with flat faces sits a
   chord sagitta inside the true geometry. On the annular fixture this is the
   DOMINANT error in the computed conductance (≈0.35% at the finest grid tested)
   and it is a geometry error, not a discretization error — refining the
   circumferential direction is what fixes it.
6. **Coarse boundary-layer resolution under a large Biot number.** The P₁ space
   cannot resolve a thermal boundary layer thinner than an element. Nothing in
   the report detects this; the residual will be small and the answer wrong.
7. **A linearized radiation row used too far from its mean.** The constructor
   requires a finite departure budget and evaluation outside it refuses. Inside
   the domain, the report still exposes the measured difference from the full
   `T⁴` law; the linearization is not silently relabelled exact.
8. **An invalid enclosure matrix or exhausted outer fixed point.** Non-finite,
   out-of-range, non-closing, nonreciprocal, or area-mismatched view factors
   refuse before solving. Coupled conduction refuses if its declared
   surface-temperature iteration budget is exhausted.

## Discrepancy

**No model-form discrepancy band is claimed** (`discrepancy_rel` is not set).
That is deliberate. Quantifying it requires validation against an external
corpus with a recorded validity domain — `fs-vvreg` territory, and the honest
reason nothing in this workspace is L4.

What IS reported, per solve, and exactly what each means:

| quantity | what it establishes |
| --- | --- |
| `ConductionReport::final_residual` / `residual_threshold` | the nonlinear residual met the DECLARED stop rule. Not an error bound |
| `LinearSolveEvidence::true_relative_residual` | `‖b − Ax‖₂/‖b‖₂`, RECOMPUTED by this crate, not the Krylov recurrence estimate |
| `EnergyBalance::closure_w` | consistency between the assembled operator and an independent post-integration of the same boundary data — identically `−Σ_free r_i`. It does not check the data itself |
| `material_provenance` / `material_receipts` | whether every conductivity number carried an `fs-matdb` receipt, and how many travel with the solve |
| `LinearizedRadiationPoint::discrepancy_w_m2` | pointwise difference between the declared linearization and the full surface-to-ambient `T⁴` expression; not a global model-form band |
| `RadiosityReport::linear_residual_max_w_m2` / `relative_energy_closure()` | algebraic radiosity residual and closed-enclosure heat balance for the admitted matrix; neither validates the view-factor generator or surface assumptions |
| `CoupledRadiationReport::surface_update_history_k` / `final_threshold_k` | whether the partitioned outer fixed point met its declared temperature rule; not a nonlinear error bound |

The G1 orders (`tests/mms.rs`) are OBSERVED convergence rates on the fixture
ladders, gated by `fs_mms::OrderGate`. They are not proven orders and not error
bounds. Read `CONTRACT.md`'s "Evidence scope" before quoting one — in
particular, P₁ on these meshes reproduces cubic solutions at the nodes, which
makes any ladder run on a low-degree manufactured solution a measurement of the
INTERPOLATION error rather than of the scheme.

## Verification status

| tier | status | where |
| --- | --- | --- |
| G0 algebraic laws | green | `tests/conformance.rs` plus `tests/contact.rs` and `tests/radiation.rs` (card receipts, view-factor row/reciprocity admission, unequal-area three-surface enclosure balance, typed refusals) |
| G1 manufactured solutions | green, 5 ladders | `tests/mms.rs` |
| G2 canonical benchmarks | partial | `tests/analytic.rs` and `tests/radiation.rs`: slab, slab+source, Dirichlet–Robin slab, cylindrical/spherical shells, straight fin, and two-surface radiosity with a caller-declared limiting matrix — closed forms, not a community benchmark suite; no geometry kernel computes the parallel-plate view factor |
| G3 metamorphic | not run | no metamorphic battery exists for this crate |
| G4 cancellation | green | `tests/conformance.rs` cancellation drills plus radiation radiosity/coupling refusal checks |
| G5 determinism | partial | same-ISA conduction replay and snapshot-resume are bitwise; calling the pure radiosity solve twice is not separate determinism evidence, and there is no registered golden or cross-ISA audit |

## Maturity

Against `docs/MATURITY_LEVELS.md`, the capability `thermal.conduction-solve`
now meets the **L2 (numerically verified)** bar as far as this crate's own
evidence goes: named Gauntlet tiers are green, the evidence refs resolve to
currently-passing tests, and the checks are against independent closed forms and
convergence orders rather than a golden hash alone.

Promoting the registry entry is a governed event and is NOT done by this bead:
`capability-maturity.json` is outside its file boundary. The evidence a
promotion would cite is listed above.

L3 (integrated workflow) additionally wants an e2e lane crossing a crate
boundary under real types; the conjugate-coupling and QoI beads are where that
lands. L4 wants an external validation corpus, which does not exist yet.
