# CONTRACT: fs-ladder

The fidelity-graph and legacy-ladder registry: the shared substrate that the
speculation proposer, query planner, Goodhart guard, discrepancy probes, and
tolerance allocation walk.

## Purpose and layer

Layer L3 (FLUX-adjacent). Pure abstraction with one lower utility dependency:
the safe-Rust, dependency-free `fs-blake3` identity owner. It does not depend
on L6 `fs-plan` or on concrete physics solvers. Those producers are represented
by typed artifact references and resolved through a consumer-supplied trait.

The native v2 model is a graph. Nodes identify model implementations and their
governing model cards. Directed edges are contextual evidence that a target is
more informative than a source for a declared QoI/regime. Cost is orthogonal:
the target may be cheaper, and expense never creates authority. The original
ordered ladder remains supported exactly as a path-graph specialization.

## Public types and semantics

- `ModelId`, `ModelCardRef`, `CostModelRef`, `DiscrepancyModelRef`,
  `TransferRef`, `QueryEvidenceRef`, `EdgeId`, and `FidelityGraphId` are
  non-confusable wrappers over exact 32-byte roots.
- `QoiId` and `RegimeAxis` are bounded visible-ASCII semantic names.
- `ClosedInterval`, `QoiSelector`, `ContextClause`, and
  `ContextPredicateSet` encode deterministic QoI/regime predicates. Axes
  inside a clause are conjunctive; clauses are disjunctive. Declaration order
  and exact duplicate clauses are nonsemantic. An empty set is explicit
  unknown, never universal.
- `ValidityDomain` says where an edge comparison has evidence.
  `Informativeness` says where the target is evidenced as more informative.
- `FidelityNode` binds an exact model identity to an exact model-card identity.
- `FidelityEdge` binds source/target, an `fs-plan` cost-model reference, an
  `fs-evidence` discrepancy-model reference, a transfer reference, validity,
  and informativeness. Native self-loops refuse. Legacy embeddings retain the
  old advisory scalar and state discrepancy evidence as `UnknownLegacy`
  instead of inventing an artifact.
- `FidelityGraph` validates card-node closure, owns canonical sorted nodes and
  edges, emits/decodes bounded canonical transport, and exposes a
  domain-separated BLAKE3 identity. Decode refuses unknown schema versions,
  malformed/trailing bytes, noncanonical encodings, and body/identity
  mismatches.
- `QueryContext` binds QoI, exact finite regime coordinates, problem size, a
  seconds budget, and the maximum admitted relative discrepancy.
  `EdgeEvidenceResolver` is the adapter boundary where a higher layer resolves
  the exact referenced cost/discrepancy artifacts. The resolver supplies a
  numeric assessed relative discrepancy (or explicit unknown); the graph
  derives `Adequate`/`Inadequate` against the query tolerance. Resolutions with
  mismatched references are treated as unresolved.
- `best_model_for(start, context, resolver)` selects a reachable contextual
  graph maximum inside budget. A unique maximum is named directly. Multiple
  incomparable maxima remain named as incomparable; cost then model identity
  is an operational replay tie-break only.
- `cheapest_adequate(start, context, resolver)` selects the cheapest reachable
  source only when every applicable pairwise comparison is resolved and
  within the exact query tolerance. Unknown, contradictory, and inadequate
  assessments do not widen into adequacy; multiple cost estimates compose by
  conservative maximum.
- `ModelRecommendation` carries the selected model, predicted cost, exact
  start/QoI/regime/problem-size/budget/tolerance context, edge path, graph
  identity, matched specificity, every considered/resolved edge, unresolved
  edge identities, incomparable maxima, and selection basis.
- `GraphTransfers` generalizes prolongation/restriction to exact graph edges.
- `LadderDescriptor` and `EmbeddedLadderGraph` retain every legacy rung field
  and move each runtime transfer to the corresponding path edge.
- `Rung { index, name, relative_cost, note }` remains the legacy fidelity
  declaration; index 0 is the coarsest rung and the scalar cost is advisory.
- `Transfer` (trait) — `prolongate(coarse) -> fine` (rung k→k+1) and
  `restrict(fine) -> coarse` (rung k+1→k). Consumers supply the operator
  (fs-feec transfer, a correlation model, …).
- `Ladder` — builder `Ladder::new(kernel, base_name, base_cost, base_note)`
  then `.then(transfer, name, cost, note)` per finer rung (keeps
  `transfers.len() == rungs.len()-1` by construction). Query: `rung(k)`,
  `top()`, `bottom()`, `adjacent_rungs(k) -> AdjacentRungs { coarser, finer }`
  (empty in the correct direction at the ends), `prolongate(from, coarse)`,
  `restrict(from, fine)`.
- `LadderRegistry` — `register(ladder)`, `ladder(kernel)`, `kernels()`;
  `cht()` seeds the Proposal-7 conjugate-heat-transfer ladder
  (`correlation-Nu` → `RANS` → `LES`); `cht_graph()` migrates the same
  declaration and transfers into the graph model.
- `Refine1d` — a concrete 1D coarsen/refine-by-2 demonstrator: prolongation is
  linear interpolation (`n → 2n-1`), restriction is injection (`2n-1 → n`), so
  `restrict ∘ prolongate = identity` and `prolongate ∘ restrict` is an
  idempotent projection.
- `LadderError` — `NoSuchRung` / `AtTop` (prolongate at the finest rung) /
  `AtBottom` (restrict at the coarsest rung) / `NoKernel`, each with a teaching
  `Display`.

## Invariants

- Graph node/card identities are exact and unique; every edge endpoint
  resolves before insertion.
- Native edges carry separate cost, discrepancy, transfer, validity, and
  informativeness axes. None implies another.
- Empty validity or informativeness matches no query. Missing, mismatched, or
  nonfinite resolved evidence fails closed.
- Query traversal is bounded by 4,096 nodes and 32,768 edges. Context-specific
  reverse edges may coexist for different QoIs/regimes; the finite
  best-candidate update reaches a deterministic fixed point.
- Canonical graph identity binds schema version, graph/legacy identity,
  model/card identities and labels, exact legacy metadata, every edge
  reference, QoI selectors, interval endpoints by IEEE bits, and canonical
  clause order.
- Rung indices are `0..len`, strictly increasing in a legacy embedding;
  `transfers[k]` bridges rung `k` and `k+1`.
- `adjacent_rungs` returns `coarser = None` at the bottom and `finer = None`
  at the top; out-of-range indices are structured errors, never panics.
- For a consistent transfer, `restrict(prolongate(x)) == x` (bit-exact for
  `Refine1d`).

## Error model

`GraphError` refuses invalid names/numbers/intervals, duplicate axes/nodes/
edges/transfers, missing endpoints/edges/transfers, self-loops, bounded-limit
violations, schema drift, malformed/noncanonical transport, and identity
mismatch. `QueryRefusal` names unknown starts, absent applicable evidence, and
absent adequate evidence. Legacy `LadderError` retains the existing teaching
refusals for out-of-range rungs, boundary transfers, and unknown kernels.

## Determinism class

Fully deterministic for identical graph, context, and resolver outputs. Maps
and clauses canonicalize through `BTreeMap` and canonical byte ordering;
semantically equivalent signed zero normalizes to positive zero. Graph/edge
identities are domain-separated BLAKE3 roots. Query path selection uses total
float order, then exact edge/model identity; every non-epistemic tie-break is
recorded. Rung resolution and `Refine1d` transfers remain pure and
bit-identical on replay (G5).

## Cancellation behavior

None — operations are bounded, synchronous pure functions with explicit
node/edge/clause/axis/transport limits (no `Cx`). The crate does not run
physics solves or fit cost/discrepancy models.

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/ladder.rs` (Proposal 3, 12 cases): total rung ordering; out-of-range
rung structured error; adjacency empty in the correct direction at the ends;
the G0 approximation property (`restrict∘prolongate = identity`,
`prolongate∘restrict` idempotent); prolongate-at-top / restrict-at-bottom
structured errors; registry register/resolve + unknown-kernel error; the CHT
ladder (`correlation-Nu`→`RANS`→`LES`); `Refine1d` degenerate lengths; G5
determinism.

`tests/fidelity_graph.rs` (f85xj.10.1, 16 cases):

- G0 construction refusal for self-loops, missing card nodes, and duplicates;
- G3 identity invariance under node/axis/clause order and duplicate clauses;
- G3 signed-zero canonicalization in interval and query contexts;
- identity movement for QoI, regime-domain, cost, and discrepancy mutations;
- canonical encode/decode, schema/trailing-byte refusal, and edge-root binding;
- a four-node cooling graph selecting different models for mean and maximum
  temperature, including a calibrated cheaper target that outranks RANS only
  in its evidenced context;
- cheapest-adequate pairwise discrepancy selection;
- conjunctive adequacy across every applicable pairwise comparison, with
  contradictory or unresolved evidence refusing;
- explicit unknown informativeness and out-of-domain refusal;
- exact resolver-reference binding;
- refusal to traverse from an over-budget source even when its target is cheap;
- deterministic, visibly operational tie-breaking for incomparable maxima;
- evidence-specificity monotonicity and edge-removal path-depth monotonicity;
- lossless ladder descriptor and transfer embedding; and
- migration of the CHT instance while the original registry keeps working.

## No-claim boundaries

- Typed roots assert exact binding, not truth, provenance, or promotion. In
  particular, current `fs-evidence::DiscrepancyModel` does not yet expose a
  durable native artifact identity; the empirical edge-population successor
  must produce retained canonical discrepancy artifacts before native
  production edges can claim that binding.
- `EdgeEvidenceResolver` output is query evidence supplied by a higher layer.
  This crate verifies exact reference equality and finite costs; it does not
  verify the scientific correctness, freshness, signature, or authority of
  those artifacts.
- A graph edge is contextual comparison evidence, not a theorem that the
  target is globally superior. Multiple maxima remain scientifically
  incomparable even when an operational tie-break selects one.
- The legacy total-order predicate preserves old behavior but carries
  `UnknownLegacy` discrepancy evidence. Its scalar `relative_cost` remains an
  advisory hint, not an `fs-plan` fit and not epistemic rank.
- The graph does not run solves, populate empirical edges, establish cooling
  validation, or integrate `certify_or_escalate`; those remain f85xj.10.2,
  f85xj.10.3, and f85xj.10.4.
- `Refine1d` is a 1D demonstrator, not a production multigrid transfer. Real
  per-kernel operators are injected.
