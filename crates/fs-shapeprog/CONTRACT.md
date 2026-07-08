# CONTRACT: fs-shapeprog

Generative geometry program synthesis: a typed constructive-geometry DSL with a
certified rewrite engine — the discrete-invention medium.

## Purpose and layer

Layer L2 (geometry). No dependencies — self-contained AST, SDF evaluator,
rewrite engine, and s-expression parser.

## Public types and semantics

- `Geom` — a constructive-geometry program: `Empty`, `Primitive{shape, size}`
  (`Sphere`/`Cube`), `Union`/`Intersect`/`Difference`, `Offset{child, radius}`,
  `Translate{child, t}`. Builders `sphere`/`cube`/`union`/`offset`/`translate`.
- `sdf(p)` — the signed distance (union = min, intersect = max, difference =
  `max(a, −b)`, offset = `child − radius`, empty = `+∞`).
- `to_sexpr` / `parse` — a round-tripping s-expression syntax.
- `canonical` / `canonical_hash` — commutative operands sorted; equivalent
  programs share a content hash (archive/ledger dedup).
- `simplify(&Geom, tiny_offset_tol) -> Simplified { program, rewrites,
  max_error }` — rewrite to a fixpoint under geometric identities; each
  `Rewrite` carries a `Certificate` (`Exact` or `Approximate{bound}`).
- `max_sdf_discrepancy(a, b, samples)` — the rewrite-safety check.
- `linear_repeat` / `stochastic_repeat` — shape-grammar productions (seeded,
  reproducible).
- `ParseError`.

## Invariants

- SAFETY: every rewrite preserves the SDF within its certificate — exact
  identities (offset composition, union/difference/translate identities,
  transform distribution) leave the SDF unchanged; `drop-tiny-offset` changes it
  by at most `|radius|`. Verified by `max_sdf_discrepancy` over a sample grid.
- Round-trip: `parse(g.to_sexpr()) == g`.
- Canonicalization: commutative-equivalent programs share `canonical_hash`.
- Grammar derivations are reproducible from their seed.

## Error model

Structured `ParseError`; no panics.

## Determinism class

Fully deterministic: SDF, rewrites, canonical form, and grammar derivations are
pure functions of the program + seed.

## Cancellation behavior

None (synchronous pure functions).

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/shapeprog.rs` (8 cases): DSL round-trip; SDF semantics; exact rewrites
preserve the SDF exactly (sampled); identity + distribution rewrites; a
certified-approximate rewrite stays within its bound; canonicalization
deduplicates commutative programs; grammar derivations reproducible from seeds;
malformed programs rejected.

## No-claim boundaries

- The rewrite engine is a rewrite-to-fixpoint simplifier over a fixed identity
  set; a full geometric E-GRAPH with equality saturation, program
  mutation/crossover for evolutionary search, and INVERSE FITTING (mesh/SDF →
  program via segmentation + parameter fitting) are the fuller deliverable.
- Programs lower to `fs-rep-frep` Region DAGs with parameter-Jacobian hooks
  (program-level adjoints) — that lowering + `Dual<N>` adjoint agreement is a
  downstream integration; here programs are SDF-valued directly.
- Parameters are plain `f64`; the `Qty`-dimensioned typed parameters and the
  F-rep/GA operator families are staged.
