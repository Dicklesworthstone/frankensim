# CONTRACT: fs-lbm

Lattice Boltzmann core (D2Q9 BGK) with the lattice-scaling assistant.

## Purpose and layer

Layer L3 (FLUX). Depends only on `fs-evidence` (the `Color` for the
Evidence-typed scaling plan). Pure, deterministic (fixed cell order).

## Public types and semantics

- D2Q9 constants: `Q` (9), `CS2` (`1/3`), `MACH_LIMIT` (`0.3`).
- `equilibrium(rho, ux, uy) -> [f64; 9]` — the D2Q9 equilibrium distribution
  (recovers `Σf = ρ`, `Σeᵢfᵢ = ρu`).
- `Lbm::channel(nx, ny, tau, gx)` — a body-force-driven channel (periodic x,
  bounce-back y-walls). `step`/`run` (collide + Guo forcing + stream +
  bounce-back); `density`, `velocity` (Guo-corrected), `total_mass`,
  `viscosity` (`ν = (τ−½)/3`), `x_velocity_profile`.
- `plan_scaling(reynolds, char_length_lu, u_lattice) -> ScalingPlan { tau,
  viscosity, u_lattice, mach, tau_margin, stable }` — the lattice-scaling
  assistant. `ScalingPlan::color()` (verified when comfortably stable, else
  estimated). Panics on non-positive Reynolds / length.
- `poiseuille_analytic(gx, viscosity, ny, y)` — the analytic reference profile.

## Invariants

- The equilibrium recovers its density + momentum moments exactly.
- MASS is conserved by a step (collision, forcing, streaming, bounce-back all
  conserve mass).
- Steady Poiseuille channel flow matches the analytic parabola to a few percent
  (halfway bounce-back resolves the quadratic profile).
- `plan_scaling` derives `τ = 3ν + ½`, flags `stable` iff `τ > ½` AND
  `Mach < MACH_LIMIT`.

## Error model

Total functions; the only panic is a nonsensical scaling request (non-positive
Reynolds / length).

## Determinism class

Fully deterministic: fixed cell iteration order; no RNG.

## Cancellation behavior

None here (a step is synchronous); polling at tile boundaries under `Cx` is the
production kernel's concern.

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/lbm.rs` (7 cases): equilibrium moments; mass conservation; Poiseuille
flow matches the analytic parabola (symmetric, centered); the scaling assistant
derives τ + flags stability + colors the plan; it rejects a high-Mach plan and
nonsense inputs; determinism.

## No-claim boundaries

- v0 is D2Q9 BGK on a DENSE grid with a body force + bounce-back walls. The
  full core — D3Q19/D3Q27, sparse FrankenVDB tiles, CUMULANT / central-moment
  collision (BGK's high-Re replacement), interpolated Bouzidi curved boundaries
  sampled from SDF charts, momentum-exchange drag/lift, and the bandwidth
  roofline / fs-tilelang kernels — is staged.
- The scaling assistant covers the `τ`/`ν`/`Mach` core; consuming fs-regime's
  dimensionless groups and emitting a full `dx`/`dt` unit conversion with
  Evidence provenance is the fuller deliverable.
