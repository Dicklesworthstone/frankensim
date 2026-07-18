# CONTRACT: fs-conform

The restriction-map plugin conformance SDK (plan addendum, Proposal 7):
certify chart-to-chart converters into tiers by the sheaf axioms. Owns risk R6.

## Purpose and layer

Layer L2 (restriction maps / Rep Router edges). Pure Rust; runtime arithmetic
uses L0 `fs-math`, and `fs-propcheck` is test-only.

## Public types and semantics

- `Converter` (object-safe trait) — `id`, `source_dim`, `target_dim`,
  `apply` (source→target), `adjoint` (declared transpose), `declared_error`.
- `ManufacturedCase { input, exact_output }`; `Composition { after, direct,
  probes }` (a functoriality witness); `ConformanceSuite { adjoint_pairs,
  manufactured, composition, identity, tolerance }` (`identity` is an optional
  witness — probes on which a converter claiming to be the identity must act as
  one; `None` for non-identity converters).
- Axiom checks: `check_functoriality` (`after∘self == direct`), `check_identity`
  (`id(x) == x`), `check_adjoint` (`⟨Ax,y⟩ == ⟨x,Aᵀy⟩`),
  `check_tolerance_honesty` (exact norm ≤ declared + suite tolerance) →
  `(honest, outward_measured)`. Each check requires a nonempty witness collection
  and rejects malformed dimensions, non-finite evidence, and arithmetic outside
  the robust evidence rung. Dot products use exact-classified f64 products and
  exact-checked double-double accumulation. Norm booleans compare exact DD
  coordinate-difference squares in fixed full-range positive/negative
  superaccumulators, without taking a square root. The separately reported
  measurement uses a power-of-two-scaled DD sum and outward f64 square-root
  bound. Detectable loss beyond that reporting rung refuses instead of rounding
  a counterexample away.
- `certify(&converter, &suite) -> ConformanceReport` — awards a `Tier`
  (`Rejected` / `Bronze` / `Silver` / `Gold`); reaches a tier above `Rejected`
  ONLY with nonempty adjoint and manufactured evidence and by passing every
  supplied axiom, with the level set by the exact admitted bound `declared +
  suite tolerance` after that bound honestly contains the exact manufactured
  error. Charging tolerance to the tier prevents a loose suite policy from
  laundering an understated declaration into Gold.
  Optional composition and identity witnesses may be absent, but a supplied
  witness must contain at least one probe. `ConformanceReport::certified()`.

## Invariants

- A converter that fails ANY axiom (functoriality witness, adjoint, tolerance
  honesty) is `Rejected` — never certified.
- Tolerances and declared errors are finite and non-negative. Every witness and
  converter result has exactly its declared dimension and contains only finite
  values. Invalid evidence fails closed; malformed manufactured evidence reports
  positive infinity rather than a misleading finite error. Admitted finite
  evidence inside the reporting rung retains a finite measured error even when
  ordinary f64 accumulation would erase a counterexample; unrepresentable DD
  arithmetic refuses. The public f64 measurement is rounded outward; exact
  perfect-square norms avoid an unnecessary successor step. Boolean norm
  decisions instead compare exact squared values and retain an f64-absorbed
  `declared + tolerance` DD sum, so a coarse public f64 projection cannot erase
  a valid sub-ULP admitted bound.
- A converter whose exact manufactured error exceeds its declared error plus the
  explicit suite tolerance FAILS tolerance honesty (a dishonest converter
  cannot buy a tier by understating beyond that admitted tolerance).
- Tier is monotone in the exact admitted bound (`declared + suite tolerance`)
  among conformant converters.
- The same `certify` path applies to first-party and third-party converters
  (R6: identical severity).

## Error model

The SDK expresses malformed, non-finite, or numerically inconclusive evidence as
`false` or a `Tier::Rejected` report rather than manufacturing a certificate.
The current `Converter` callback surface is arbitrary caller-supplied Rust: it
can panic, allocate, block, mutate interior state, or otherwise fail outside
this return-value model. This crate does not yet contain those callback faults.

## Determinism class

SDK control flow and evidence arithmetic are deterministic for the same finite
callback transcript. The object-safe `Converter` trait does not enforce purity,
state isolation, or deterministic callback results, so a current tier is not by
itself G5 evidence about an untrusted implementation.

## Cancellation behavior

None in the current synchronous API. There is no bounded cancellation point
inside a converter callback and no way to preempt a callback that fails to
return.

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

Two private arithmetic unit tests directly cover cross-word `u128` placement,
a 69-limb carry chain, capacity-overflow refusal, signed DD tails, and
maximum-finite exact squares. `tests/conform.rs` adds sixteen fixed Proposal 7
regression tests, two generated G0 laws, and one generated G3 conversion-path
relation (512 cases each, deterministic seeds and integrated shrinking). The
fixed tests cover adjoint
consistency, tolerance honesty, correct and incorrect functoriality witnesses,
identity recognition and false-identity rejection, exact admitted-bound tiering
that charges suite tolerance, and uniform first-party/third-party R6 severity.
They also reject empty or dimensionally malformed witnesses, non-finite
policy/evidence, overflowing arithmetic, and robust scale-disparate arithmetic
that would otherwise erase a nonzero failure, while accepting the same evidence
once its conservative bound includes the retained residual. The generated G0
laws exercise exact
functor composition and identity action between the fixed pins. The G3 SDK
harness declares `restriction-map-direct-vs-composed` with
`conversion_path_independence` and applies it to the test-local `Mtx`
implementation of the public `Converter::apply` trait (seed
`0xC0F0_4A48_0003`, 512 cases, exact component tolerance): a manufactured
direct matrix converter versus the staged `f(g(x))` route over nonzero,
exactly representable small-integer cases. The joint shrink surface retains the
matrix/probe input and route transform. This is conformance coverage of the SDK
trait harness, not a production geometry-conversion adopter.

## No-claim boundaries

- Converters are modeled as finite-dimensional linear operators for the
  conformance harness; the SDK surface (the `Converter` trait) is what a real
  chart-to-chart trace/conversion operator implements.
- The generated G0 laws use small, exactly representable 2x2 integer matrices
  and probes. They do not certify arbitrary nonlinear converters or general
  floating-point associativity.
- Norm booleans use a fixed full-binary64 exact superaccumulator. Dot products
  currently require their exact running sum to remain representable as DD and
  refuse otherwise; the public measurement path likewise refuses when its
  exact two-component squared sum cannot be retained. A successor certifier
  must generalize exact/binned dot accumulation, expose typed rung/refusal
  evidence and admitted work budgets, and freeze the arithmetic receipt before
  claiming theorem-strength arbitrary-f64 certification across the full API.
- The G3 path-independence SDK harness covers the same test-local `Mtx`
  finite-dimensional linear fixture family. Exact component equality is
  justified by its bounded integer arithmetic; no production chart or
  geometry converter is invoked, and the result does not extend that tolerance
  to nonlinear, approximate, or floating-roundoff-sensitive chart conversion
  paths. The declaration is not evidence of a passing batch run.
- The suite is SUPPLIED here (probes, manufactured cases, composition witness);
  AUTO-GENERATING it from a chart pair's sheaf axioms is the generator's job (a
  downstream producer feeding this harness).
- The current callback ABI does not isolate panics, allocation failure,
  nontermination, interior mutation, nondeterminism, or time-of-check versus
  time-of-use implementation changes. Certification is therefore evidence about
  the observed callback transcript, not yet a fault-contained or replay-bound
  capability for arbitrary third-party code. A successor must use an admitted,
  fallible execution protocol; bind implementation identity, budgets, and the
  canonical transcript into the receipt; and prove G4/G5 containment before
  removing this no-claim boundary.
- Stamping the certified tier onto every ledger entry the converter touches is
  fs-ledger's integration; this crate produces the tier.
- Adjoint consistency is checked against the converter's OWN declared adjoint;
  cross-checking it against Proposal 1's ledger adjoint is the ledger's wiring.
