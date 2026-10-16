# CONTRACT: fs-symmetry

Symmetry harvesting (plan addendum, Proposal 13): harvest cyclic symmetry as
both correctness and speed via isotypic (DFT) block-diagonalization.

## Purpose and layer

Layer L1 (spectral/linear-algebra numerics). No dependencies — pure Rust,
including a minimal internal complex type and a naive DFT.

## Public types and semantics

- `CyclicGroup { n }` — `character(irrep, element)` is the `Cₙ` character table
  (roots of unity).
- `circulant_matvec(first_row, x)` — the circulant `C[i][j] = first_row[(j−i)
  mod n]` times `x`.
- `solve_circulant(first_row, rhs)` — solve `C x = rhs` by isotypic
  block-diagonalization `x = IDFT(DFT(rhs) / eigenvalues)` (`O(n²)` vs the
  `O(n³)` dense solve — the `k`-fold ≈ `k`× win); errors `Singular` on a zero
  eigenvalue.
- `cyclic_residual(v, k_fold) -> SymmetryResidual` — the certified asymmetry
  residual `||v − rotate(v, len/k)||` (absolute + relative + `is_exact`).
- `symmetrize(v, k_fold) -> (sym, asym)` — the isotypic projection onto the
  symmetric subspace and its remainder; `sym` is exactly `k`-fold symmetric.
- `symmetrized_solve(first_row, rhs, k_fold) -> PerturbationBound` — solve the
  symmetric part of `rhs` and return a CERTIFIED bound `asymmetry_residual /
  λ_min` that provably contains `||x_full − symmetric_solution||`.
- `SymmetryError` — `EmptyInput` / `LengthMismatch` / `ZeroFold` /
  `NotDivisible` / `Singular`.

## Invariants

- `solve_circulant` returns the exact full-system solution (verified by
  `circulant_matvec(first_row, x) == rhs`) — the block-diagonal solve is
  bit-consistent with the full solve up to floating tolerance.
- `symmetrize`'s symmetric part is exactly `k`-fold symmetric (zero residual).
- The perturbation `correction_bound` CONTAINS the true correction the
  asymmetric remainder induces — approximate symmetry never yields a silent
  wrong answer, only a symmetrized solve plus a certified correction bound.

## Error model

Structured `SymmetryError` values; no panics.

## Determinism class

Fully deterministic: transforms, solves, residuals, and bounds are pure
functions of the inputs.

## Cancellation behavior

None (synchronous pure functions).

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/symmetry.rs` (Proposal 13, 8 cases): cyclic character table; the
block-diagonal solve inverts a NON-symmetric circulant and scales to a 12-node
ring; singular + shape rejection; certified asymmetry residual + error paths;
isotypic symmetrization (symmetric part is exactly symmetric); the certified
perturbation bound contains the true correction; determinism.

## No-claim boundaries

- v1 covers the CYCLIC group `Cₙ` and circulant operators (the algebraic
  signature of cyclic symmetry). Dihedral `Dₙ` character tables, and detection
  via graph automorphism (FrankenNetworkx) + geometric hashing, are follow-ons.
- The DFT is a naive `O(n²)` transform; a production build uses fs-la's FFT and
  the geometric-algebra (motor/screw) group representation from fs-ga.
- The perturbation bound treats the asymmetry as a right-hand-side remainder;
  an operator that is only APPROXIMATELY circulant (asymmetry in the operator,
  not just the RHS) is a further generalization.
- Symmetry detection here is by residual on a supplied field; harvesting a
  DECLARED group from the type system (Proposal 13 interfaces) is the caller's.
