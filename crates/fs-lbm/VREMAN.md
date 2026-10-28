# Direction-aware Vreman SGS transport

`SgsModel::Vreman { cv }` selects the Vreman (2004) algebraic eddy-viscosity
closure through the existing explicit SGS owners. Nothing enables it by default.
`cv` is the paper's LINEAR coefficient, not a Smagorinsky coefficient and not
its square root. `cv = 0` retains the exact ordinary molecular stepping path.

## Physics and numerical evaluation

Primary source: A. W. Vreman, *An eddy-viscosity subgrid-scale model for turbulent
shear flow: Algebraic theory and applications*, Physics of Fluids 16 (2004),
3670-3681, doi:10.1063/1.1785131, equations (5)-(8).
Author-hosted paper: https://www.vremanresearch.nl/Vreman-PF2004-subgridmodel.pdf

Using this repository's convention `G[i][j] = du_i/dx_j` (the transpose of the
paper's alpha), the model is

```text
beta = G diag(Delta_x^2, Delta_y^2, Delta_z^2) G^T
B_beta = beta_11 beta_22 - beta_12^2
       + beta_11 beta_33 - beta_13^2
       + beta_22 beta_33 - beta_23^2
nu_t = cv sqrt(B_beta / (G:G))
```

The implementation instead evaluates `B_beta` as the sum of the squares of
nine weighted 2-by-2 gradient minors (Cauchy-Binet). This avoids subtracting
nearly equal Gram products, so near-rank-one gradients retain small positive
transport without a negative-invariant clamp. Explicit fused product residuals
improve the minors. Power-of-two gradient normalization and exponent-tracked
length products avoid fourth-power intermediate overflow and avoid losing
extreme-aspect-ratio terms merely by normalizing all widths by their maximum.
Unrepresentable nonzero transport refuses instead of silently becoming laminar.

The local kernel has fixed-size stack storage, no heap allocation, no test
filter, no spatial averaging and no wall-distance lookup. This is an algorithmic
property, not a measured speedup or a complete runtime benchmark.

## Existing solvers, one transport owner

- `sgs3::eddy_viscosity` passes equal widths for scalar-filter lattice calls;
  `sgs3::vreman_eddy_viscosity` accepts three widths explicitly.
- Generic D3Q19/D3Q27 `Grid3` and D3Q19 `BoundaryGrid3` use the existing SGS
  preparation/streaming path. Shear changes; independently owned bulk and
  higher-order relaxation remain as before.
- `NsfGrid3::with_sgs` uses the same modeled shear in momentum and total-energy
  stress work, and the existing `rho cp nu_t / Pr_t` heat flux. No additional
  dissipation-heating source is introduced. Existing low-Mach, open-boundary
  and transport-envelope restrictions remain.
- `simple_unsteady_with_sgs` and `march_conjugate_with_sgs` use the actual
  `filter_width_factor * width[axis]` of each graded cell. They retain the
  native MAC gradient reconstruction, full symmetric momentum stress, and
  thermal conductivity with `Pr_t = 0.9`. No scalar volume-equivalent width
  replaces these directional widths. Existing Smagorinsky/WALE widths and
  arithmetic are unchanged. Velocity scratch is reused, and Vreman does not
  allocate a redundant per-cell three-width cache.

The existing eddy/molecular ratio gates still refuse excess transport. LVEL
and equilibrium thermal wall laws are not silently overlaid on this LES path.

## Complete cooling application

```sh
cargo run --release -p fs-cli --bin frankensim -- \
  --json cooling-cht examples/cooling-cht/vreman-graded-obstacle.json
```

This is a 256-cell graded 24 x 12 x 8 mm duct containing a 0.5 W aluminium
obstacle, ramping to a 2 m/s inlet over 40 coupled steps of 0.5 ms. The authored
coefficient is 0.07. This short startup illustrates the complete application;
it is not a developed turbulent spectrum or a converged cooling prediction.

The JSON selector is `solver.sgs.model: "vreman"`, with the same five required
controls as the other SGS models. Results add `filter_width_definition:
"directional-cell-widths"` and minimum/maximum widths by axis under `flow.sgs`.
The scalar min/max fields bound all three physical widths, not their geometric
mean. Existing model output layouts are retained. VTK exports the actual final
velocity, temperature, eddy viscosity and turbulent conductivity.

## Checks and limits

```sh
cargo test -p fs-lbm --lib sgs3
cargo test -p fs-lbm --lib vreman_tests
cargo test -p fs-lbm --test vreman_sgs
cargo test -p fs-cli --test cooling_cht sgs::
```

Six kernel groups cover the 13 laminar structural classes, directional-axis
ownership, isotropic rotations, quadratic filter refinement, nearly rank-one
transport and extreme dimensional scaling/refusal. Three FV groups add native
gradients on a graded mesh, cache reuse, caps, exact-zero evolution and live
momentum/heat coupling. Two NSF groups cover the zero control and combined
shear/bulk/heat ownership, closed energy and worker-count identity. Two lattice
groups cover both stencils, unchanged collision rates, changed populations and
mass. CLI coverage adds the graded example and expands existing dispatch and
zero-control tests to Vreman.

Authoring-session evidence: 1,597 finite nonzero random cases from an independent
Python transcription agreed with an 800-digit Gram-tensor oracle to maximum
relative difference 4.45e-16; all 13 laminar classes and the named near-rank-one
and extreme-scaling controls passed. Rust compilation, these integration tests
and the complete example have NOT been executed in this environment because
rustc/cargo are unavailable. The mathematical checks are not native execution.

Vreman's algebraic zero-dissipation properties do not supply transition physics,
roughness, non-equilibrium/separated wall treatment, a subgrid kinetic-energy
transport equation or a validated turbulent heat flux. A coefficient choice
and a conservative run do not establish turbulent accuracy. Developed-flow
benchmarks, mesh/time refinement and empirical heat-transfer validation remain
necessary. Acceptance-render use continues to refuse.
