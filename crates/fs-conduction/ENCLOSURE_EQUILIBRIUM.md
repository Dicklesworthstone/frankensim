# Enclosure equilibrium and reusable radiosity

`solve_with_gray_diffuse_enclosure` uses one prepared factorization of the
immutable emissivity/view-factor matrix. Independent temperature right-hand
sides retain the original pivot and FMA order. `GrayDiffuseEnclosure::solve`
remains available for a single evaluation; `prepare_radiosity` exposes reuse.
Neither operation changes the gray-diffuse area-mean physical model.

Every initial and radiation-loaded solid solve must pass its original residual
threshold. A small-step termination alone is not enclosure equilibrium.
After the surface-temperature update criterion passes, the driver recomputes
radiosity at the actual candidate field, rebuilds those Neumann loads, and
reassembles the original solid/contact equations at that same temperature.
The norm of the free residual must also pass the original inner threshold.
Heterogeneous and temperature-dependent material assignments, source, contact
and prescribed boundary data remain unchanged. There is no extra solid solve.

A failed final-state check continues the existing bounded iteration. It does
not increase the iteration cap, rescale the residual allowance or return a
partial equilibrium. Cancellation and physical input failures propagate.

The returned radiation report retains both:

* `radiosity`: the loads actually applied to the last solid solve;
* `equilibrium.nonlinear_radiosity`: radiation recomputed at the returned field.

`equilibrium` also contains the actual nonlinear free residual, its unchanged
threshold, and the largest applied-versus-recomputed patch heat discrepancy.
`equilibrium_checks` counts completed final-state residual assemblies, including
failed acceptance checks. These are not additional solid solves or Krylov
iterations. Patch heat discrepancy is a diagnostic, not an independent patch
heat tolerance. Prescribed rows are reactions, not free equilibrium equations.

This is a discrete-equation acceptance check, not a temperature-error bound,
continuum certificate, independently validated view-factor model, or experiment.
The library report gained `equilibrium` and `equilibrium_checks`; downstream
code constructing report literals must supply them. Solver input types did not
change. This does not add a native project grammar or another cooling command.

```sh
cargo test -p fs-conduction --test enclosure_solver --test radiation
```

Four prepared-factor regressions and three final-state equilibrium regressions
join the two existing solid-acceptance tests. They exercise the public solver,
full-field replay, nonlinear assignments and a temperature-only candidate with
closed total radiation but large actual equilibrium error. The one-shot-factor
comparison now uses the same stronger physical acceptance on both paths.

Native compilation and these tests have not been executed in the authoring
environment: Cargo/rustc are absent. Independent Python checks ran 186 exact
factor-reuse comparisons and 12 two-slab equilibrium comparisons. On one
illustrative linear case the first temperature-only candidate has a 624.524 W
free residual and 13.3224 K interface error despite zero summed radiation;
the residual-gated reference agrees with a simultaneous solve within 3.01e-8 K.
These references are not execution of the Rust implementation or benchmarks.
