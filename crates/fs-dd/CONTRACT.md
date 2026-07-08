# fs-dd — CONTRACT

Domain decomposition (plan §8.9, bead tfz.11): BDDC substructuring with
corners-primal coarse spaces, sheaf-derived edge enrichment, and
CCD-aligned partitioning metrics.

## Purpose and layer
Layer L3. The extreme-core-count solver path; single-machine p-MG covers
the committed scale (the plan's staging note).

## Public types and semantics
- `Decomposition`: the structured 2-D model problem — `s × s` square
  subdomains of `m × m` cells, variable per-cell coefficients
  (`uniform`, `checkerboard` fixtures), 5-point stencil with averaged
  edge coefficients, Dirichlet rim. `apply_global` is the whole-system
  oracle.
- `Bddc::new(decomp, with_edge_modes)`: substructures, factors interior
  blocks, local free-boundary Schur blocks, and the coarse Galerkin
  matrix ONCE; `with_edge_modes` adds one average constraint per open
  interface (the sheaf-derived enrichment) and deflates the local
  corrections against them so the two spaces do not fight.
- `schur_apply` (never assembled), `precondition` (weighted local +
  coarse corrections), `solve_cg` (preconditioned CG with a Lanczos
  condition estimate), `ccd_locality` (the topological locality metric
  for island-aligned partitions), `gamma_len`/`coarse_dim`.

## Invariants
- The preconditioner is SPD and symmetric on probes (dd-001).
- Condition numbers track the BDDC signature: `κ/(1+log(H/h))²` stays
  in a bounded band across H/h ∈ {4, 8, 16} (dd-002).
- Checkerboard 1e6 coefficient jumps stay within 3× the uniform
  iteration count (subdomain-aligned jumps — the BDDC-friendly case,
  noted; dd-003).
- The MEASURED coarse-space comparison (dd-004): on the adversarial
  jump fixture the sheaf-edge enrichment strictly improves the
  condition estimate (≈20 vs ≈33 recorded) at comparable iterations;
  on the uniform fixture the trade (comparable κ, slightly more
  iterations for a 33-vs-9 coarse dimension) is ledgered honestly.
- Deflation lesson (recorded): NAIVE edge enrichment DEGRADES BDDC —
  local corrections must be projected away from the coarse edge space.

## Error model
Construction panics on non-SPD local blocks (a mesh/coefficient bug,
not a runtime condition). All solves are total.

## Determinism class
Fully deterministic: fixed assembly order, dense factorizations,
fixed-seed test probes. Subdomain iteration order is row-major.

## Cancellation behavior
Each subdomain factorization/solve is an independent unit (the scoped
resumability contract); the fs-exec tile integration lands when this
crate leaves its feature gate.

## Unsafe boundary
No `unsafe` anywhere in this crate.

## Feature flags
`bddc` ([M], default OFF) gates everything; `sheaf-coarse` (implies
`bddc`, pulls fs-geom) gates the Bet-11 cross-checks.

## Conformance tests
tests/bddc.rs — dd-001 G0 preconditioner properties; dd-002 log²(H/h)
scaling; dd-003 jump robustness; dd-004 the measured sheaf-vs-corners
comparison; dd-005 the sheaf cross-check (the 2×2 subdomain adjacency
is a 4-cycle whose 1-D harmonic space is exactly the mode the corner
constraint pins — the sheaf explains WHY corners are primal) +
CCD-locality metrics.

## No-claim boundaries
- The model problem is the structured 2-D 5-point Laplacian; unstructured
  meshes, elasticity, and FETI-DP variants are follow-up scope (the
  bead names FETI-DP as infrastructure-sharing, not shipped here).
- The sheaf-harmonic coarse space is the EDGE-AVERAGE enrichment framed
  and cross-checked by Bet 11's machinery; adaptive eigenmode
  enrichment (adaptive-BDDC class) is the [M] growth path.
- `ccd_locality` is a topological metric, NOT a wall-clock claim — no
  performance numbers without benchmarks (AGENTS perf rules).
- Condition estimates come from CG's Lanczos tridiagonal — estimates,
  not certified enclosures.
