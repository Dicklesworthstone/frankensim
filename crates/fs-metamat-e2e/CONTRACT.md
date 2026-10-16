# CONTRACT: fs-metamat-e2e

MetamatCert — a certified stiffness-density frontier for a porous metamaterial.
Layer L4 (ASCENT).

## Purpose and layer

Composes `fs-lattice` (homogenization + Voigt bound), `fs-sos` (PSD certificate),
`fs-evidence` (Verified). Deps point downward.

## Public types and semantics

- `CellPoint { r, density, c11, specific_stiffness, stable, admissible }`.
- `run_campaign(n, &radii) -> MetamatReport` — homogenizes each porosity,
  certifies stability (PSD) and admissibility (≤ Voigt), and reports the frontier.
- `default_radii()` — the porosity sweep.

## Invariants

- STABILITY: every effective tensor is certified positive-definite (`is_psd`).
- ADMISSIBILITY: every `C₁₁` respects the Voigt upper bound `ρ·C₁₁ˢᵒˡⁱᵈ` — a
  violation would mean the homogenizer is wrong (certifying the certifier).
- `C₁₁` and density both fall monotonically with porosity; the Voigt bound proves
  the solid cell maximizes specific stiffness (`C₁₁/ρ ≤ C₁₁ˢᵒˡⁱᵈ`).
- The frontier color is `Verified` iff all-stable and all-admissible.
- Deterministic (fixed FE homogenization; no RNG).

## Error model

Panics only on an empty radius list.

## Determinism class

Fully deterministic (G5).

## Cancellation behavior

None (a synchronous batch).

## Unsafe boundary

None; `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/metamat.rs` (3): the frontier is stable + admissible + solid-optimal;
porosity reduces stiffness; determinism.

## No-claim boundaries

A fixture-scale square-symmetry holed-plate homogenizer (single material, void by
density scaling); the frontier is a 1-parameter porosity sweep, not a full
inverse-design; the Voigt bound is the elementary mixture bound (Hashin-Shtrikman
is tighter).
