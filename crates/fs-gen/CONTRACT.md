# CONTRACT: fs-gen

> Status: PARTIAL and moonshot-gated. The crate is empty on the default
> feature set. All public generative APIs are behind `proposal-gen`.

## Purpose and layer

Proposal-only generative models for ASCENT candidate generation. Layer:
**L4 ASCENT**. This crate is [M]-tagged and may propose seeds,
mutations, graph candidates, or archive-aware acquisition batches; it
must not certify physics, feasibility, safety, or final rankings.

## Public types and semantics

- `Proposal<T>` — quarantines a generated payload. The payload is
  private; the only non-logging exit is `promote(validate)`, which
  returns the payload only after caller-supplied validation succeeds.
- `ModelCard` — records model identity, corpus hash, and determinism
  class. It is provenance for a proposal, not scientific evidence.
- `ShapePrior` — fits a KDE-style prior to finite, fixed-dimension
  corpus vectors and proposes bandwidth-jittered design vectors.
- `MutationKernel` — fits a dominant covariance direction and proposes
  corpus-shaped mutation directions.
- `GraphGenerator` — fits degree-biased attachment weights and proposes
  simple undirected edge lists.
- `acquire` — ranks prior proposals by `density * min-distance-to-archive`.

## Invariants

1. Proposals do not expose payloads to certified paths without an
   explicit validator.
2. Corpus-fit generators reject empty, zero-dimensional, non-finite, or
   ragged vector corpora rather than silently truncating dimensions.
3. Degenerate covariance falls back to a deterministic unit direction
   rather than returning a zero "principal" axis.
4. Graph proposals are self-loop-free, duplicate-free, fill every
   requested feasible edge count, and cannot request more edges than a
   simple undirected graph can hold.
5. Archive-aware acquisition rejects archive rows with dimensions or
   non-finite values that do not match the fitted prior.
6. Sampling is deterministic per seed.
7. Corpus hashes include row counts and row lengths, so differently
   shaped corpora do not alias through byte concatenation.

## Error model

Structured panics for modeling errors: invalid corpora, invalid density
queries, zero-node graph generators, and impossible graph edge counts.
`Proposal::promote` returns `Rejected` for validator failure.

## Determinism class

Bit-deterministic on one ISA for fixed corpus, seed, and feature set.
Cross-ISA numerical equivalence is not claimed for the KDE/covariance
floating-point path.

## Cancellation behavior

None. Generators are synchronous and small in v0.

## Unsafe boundary

None. The workspace `unsafe_code = "deny"` lint applies.

## Feature flags

- `proposal-gen` — enables the [M]-tagged proposal APIs. The default
  feature set intentionally exposes no generator API.

## Conformance tests

`tests/gen.rs` (feature `proposal-gen`): proposal promotion/rejection,
shape-prior determinism and finite fixed-dimensional outputs,
mutation-kernel determinism including degenerate covariance,
graph-generator simple-edge invariants, feasible-count fill, and
impossible-edge guard, malformed corpus guards, archive-aware
acquisition determinism, and malformed archive rejection.

## No-claim boundaries

- No certified validity, feasibility, safety, or objective-quality claim.
- No learned diffusion, flow, transformer, neural, or LLM generator is
  implemented.
- No ledger database integration yet; callers provide the corpus bytes and
  validation machinery.
- No statistical calibration claim for proposal density; this is a
  candidate generator, not a surrogate model.
