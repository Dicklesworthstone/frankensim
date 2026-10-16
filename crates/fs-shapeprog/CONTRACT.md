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
- `to_sexpr` / `parse` — a round-tripping s-expression syntax. Parsed numeric
  atoms must be finite; signed zero remains valid.
- `canonical` / `canonical_hash` — commutative operands sorted; equivalent
  programs share a content hash (archive/ledger dedup).
- `simplify(&Geom, tiny_offset_tol) -> Simplified { program, rewrites,
  max_error, first_lossy_rewrite, pass_bounds, refusal }` — transactional
  rewrite to a true fixpoint. A successful result composes a root-level uniform
  SDF error bound; a refusal returns the original program, empty rewrite/bound
  traces, and zero identity error. `Rewrite` records its zero-based pass,
  deterministic root-to-node `RewritePathStep` path, local `Certificate`, and
  accumulated subtree bound.
- `Certificate::Exact` — bit-equivalent interpreter results for admitted finite
  evaluations under the documented signed-zero/non-finite policy.
- `Certificate::Approximate{bound}` — a finite nonnegative local envelope. The
  current lossy rule drops an admitted offset with local bound `2*|radius|`:
  correctly rounded subtraction may cross a rounding-cell boundary by more than
  `|radius|`, while nearest rounding and the representable unshifted value give
  the universal factor-two envelope.
- `SimplifyRefusal` / `BoundOperation` — structured fail-closed diagnostics for
  non-finite program parameters, invalid tolerance, unrepresentable outward
  bound propagation, or deterministic pass-budget exhaustion. Invalid tolerance
  diagnostics retain exact IEEE bits so signed values and NaN payloads replay
  without NaN-equality ambiguity.
- `max_sdf_discrepancy(a, b, samples)` — the rewrite-safety check, returned as
  an outward binary64 upper bound. Empty or non-finite evidence and
  unrepresentable arithmetic return `+∞` as a refusal sentinel; matching `+∞`
  is agreement only for structurally empty SDFs. For a non-structural
  evaluation, every visited branch value and translated coordinate must be
  finite: a finite selected Boolean root cannot launder an overflowing inactive
  branch into admissible evidence.
- `linear_repeat` / `stochastic_repeat` — shape-grammar productions (seeded,
  reproducible).
- `ParseError`.

## Invariants

- SAFETY: for every successful simplification and every admitted finite
  evaluation, `|sdf(original)-sdf(program)| <= max_error`.
- COMPOSITION: sequential lossy effects on one evaluation path are added with
  outward rounding. Independent `Union`, `Intersect`, and `Difference` branch
  effects use `max` under their executable finite-value 1-Lipschitz hypotheses;
  the discrepancy checker enforces those hypotheses recursively rather than
  checking only the selected root. `Translate` carries a child bound with
  factor one. A retained rounded `Offset` carries a nonzero child bound through
  the real affine factor one and adds the range-free `2*|radius|`
  nearest-rounding envelope; a zero child bound remains exactly zero because
  both rounded evaluations are identical. Pass-level root bounds are
  themselves added outward; `max_error` is not a global maximum over local
  rewrites.
- FLOATING EVALUATION ORDER: consecutive offsets are preserved. The
  real-arithmetic reassociation
  `offset(offset(a,r1),r2) -> offset(a,r1+r2)` is not `Exact` for binary64
  because it replaces two rounded subtractions with one subtraction by a
  rounded parameter sum. Exact rules are limited to interpreter-bit-equivalent
  empty identities and translation distribution over admitted finite values.
- FAIL CLOSED: a negative/non-finite tolerance, non-finite program parameter,
  non-finite/unrepresentable required bound, or failure to reach a syntactic
  fixed point rolls back the entire transaction. No partial simplified program
  receives a finite certificate.
- TRACEABILITY: the first lossy rewrite index, pass bounds, deterministic paths,
  and accumulated subtree bounds are sufficient to locate and replay the bound
  composition. Successful simplification is deterministic, monotone in the
  nonnegative tolerance, and idempotent at its returned fixed point.
- Round-trip: `parse(g.to_sexpr()) == g` for finite-parameter programs.
- Canonicalization: commutative-equivalent programs share `canonical_hash`.
- Grammar derivations are reproducible from their seed.

## Error model

Structured `ParseError` and `SimplifyRefusal` cover modeled syntax, admission,
arithmetic, and fixed-point failures; there are no intentional public-input
panics. Refusal is transactional: `program` equals the original input and
`max_error=0` is the exact distance of that rollback identity, not authority for
an attempted partial rewrite. `max_sdf_discrepancy` retains `+∞` as its
evidence-refusal sentinel.

## Determinism class

Fully deterministic: SDF, rewrite order/path/pass traces, outward-bound
composition, refusal selection, canonical form, and grammar derivations are pure
functions of the program, tolerance, and seed.

## Cancellation behavior

None (synchronous pure functions).

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/shapeprog.rs` provides focused G0/G3 coverage for DSL round-trip and SDF
semantics; exact identity/distribution rewrites; the two mandatory historical
false-certificate witnesses (sequential loss hidden by global max and
catastrophic consecutive-offset reassociation); the rounded-subtraction
threshold requiring `2*|radius|`; sequential versus branch versus unary scaling
and rounded-affine contexts; branch permutation/selection; tolerance neighbours
and monotonicity; signed zero, subnormal, sign-reversed, overflow, and
invalid-input cases; the retained-offset rounding-cell witness that falsifies
global factor-two child propagation; transactional pass-limit refusal;
deterministic trace replay; fixed-point idempotence; a node-count assertion for
the single-pass evidence walk; direct discrepancy verification;
canonicalization; seeded grammar replay; parser rejection; absorbed and exact
one-ulp outward-discrepancy subtraction; and discrepancy-evidence refusal,
including masked non-finite Boolean branches.

## No-claim boundaries

- The rewrite engine is a rewrite-to-fixpoint simplifier over a fixed identity
  set; a full geometric E-GRAPH with equality saturation, program
  mutation/crossover for evolutionary search, and INVERSE FITTING (mesh/SDF →
  program via segmentation + parameter fitting) are the fuller deliverable.
- `max_sdf_discrepancy` is a fail-closed finite-sample falsifier, not a continuum
  proof. The compositional certificate is only about the implemented SDF
  algebra under admitted finite binary64 evaluations. Admission requires every
  non-structural intermediate branch value and translated coordinate to remain
  finite; selected-root finiteness alone is deliberately insufficient. The
  certificate does not claim chart validity, topology preservation outside
  those algebraic rules, exact real arithmetic, cross-radix behavior, or
  meaningful results for NaN/infinite SDF evaluations. Approximate bounds are
  conservative envelopes, not tight-error estimates.
- Exact means bit-equivalent only inside the stated finite-evaluation policy.
  Signed zero is valid input and lossy zero-radius drops may be numerically
  zero-bounded without being relabeled as an exact reassociation. Non-finite
  program parameters are refused before rewrite admission.
- `tiny_offset_tol` is a strict local radius-admission threshold, not a global
  error budget. The returned `max_error` can exceed that threshold because a
  local context-free envelope is `2*|radius|`, sequential effects add, and
  retained rounded parents contribute their own envelopes. Consumers must use
  the returned certificate or refusal and must not infer
  `max_error <= tiny_offset_tol`.
- A retained rounded `Offset` is not claimed globally Lipschitz on the binary64
  lattice: tiny adjacent child values can straddle a much wider rounding cell
  after a common shift. Its range-free `E + 2*|radius|` bound is therefore
  intentionally conservative. Tighter propagation requires a certified child
  range/ulp analysis. General affine/scaling constructors must register their
  own outward absolute-Lipschitz transformer and rounding envelope before they
  may carry a child approximation; the engine does not infer such hypotheses
  heuristically.
- AST traversal and s-expression production are recursive and do not yet carry
  an explicit node/depth or memory budget. Adversarially deep caller-built trees
  therefore remain outside the resource-safety claim even though rewrite-pass
  count is bounded and arithmetic refusal is transactional.
- Programs lower to `fs-rep-frep` Region DAGs with parameter-Jacobian hooks
  (program-level adjoints) — that lowering + `Dual<N>` adjoint agreement is a
  downstream integration; here programs are SDF-valued directly.
- Parameters are plain `f64`; the `Qty`-dimensioned typed parameters and the
  F-rep/GA operator families are staged.
