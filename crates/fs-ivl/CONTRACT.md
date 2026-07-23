# CONTRACT: fs-ivl

## Purpose and layer
Certified arithmetic (plan §6.4): outward-rounded interval arithmetic and
affine (noise-symbol) forms. This crate defines what "certified" means
system-wide: every operation's postcondition is ENCLOSURE — the result
contains the true image of the inputs. Layer: **L1**. Depends on
fs-math (strict functions, declared budgets, nudging primitives) and the
lower-layer `fs-evidence` useful-bound vocabulary; the
in-house `eft`/`dd` reference machinery lives THERE at L0 as a single
implementation shared with fs-la — recorded relocation, beads
6ys.8/6ys.12). The independent high-precision development oracle is the
isolated, non-production `tools/oracle` MPFR harness; it is not part of this
crate's dependency graph. A second isolated, non-production
`tools/cert-kernel` workspace supplies the diverse implementation: its basic
operations use exact residuals and dyadic comparisons while its `exp` and
`ln` use bounded positive series. It does not reuse fs-ivl's outward-nudge
mechanism and is likewise absent from the production dependency graph.

## Public types and semantics
- `INTERVAL_SEMANTICS_VERSION = 1` — package-version-independent identity for
  the outward-operation/refusal policy. Certificate consumers additionally
  bind fs-math's strict-core semantic version and retained golden fingerprint.
- `Interval` — closed [lo, hi], ±∞ endpoints allowed, NaN rejected at
  construction, INVARIANT lo ≤ hi. Std operator traits (+, −, ×, ÷, unary −)
  plus `sqrt/exp/ln/tanh/sin/cos`, `abs`, set ops (`hull`, `intersect`,
  `encloses`, `contains`), `width` (rounded up — never understates),
  `midpoint` (representative only, NOT certified central).
- Directed rounding: computed endpoints are nudged outward via
  fs-math `next_up`/`next_down` — basic ops 1 ULP (rigorous: IEEE
  correctly-rounded ops err ≤ ½ ULP), elementary functions by fs-math's
  DECLARED constants directly (exp/ln/sin/cos 3, tanh 5, sqrt 1 after
  the exact hardware sqrt), with no mirrored budget literals in the interval
  implementation. NO global rounding-mode state anywhere (thread-safe,
  SIMD-mixing-safe, grep-lintable). When a basic operation on two finite
  endpoints overflows, the rounded infinity is not treated as an exact
  singleton: positive overflow becomes `[f64::MAX,+∞]`, negative overflow
  becomes `[−∞,−f64::MAX]`, and endpoint hulls compose these one-sided
  enclosures. Actual infinity operands retain extended-real semantics.
- `Interval::WHOLE` = [−∞, +∞]: the "no useful enclosure" answer. Division
  by a zero-containing interval returns WHOLE (documented; not a panic —
  certified callers branch on it). 0·∞ ambiguities in mul/div also return
  WHOLE (conservative, never wrong).
- sin/cos: non-monotone handling by conservative critical-point counting
  against an interval enclosure of π; "possibly contains a critical point"
  counts as containment. Width ≥ 2π (or beyond `det::TRIG_DOMAIN`, or
  unbounded) → [−1, 1] (honest, always correct).
- `TaylorModel1` — univariate Taylor models: f64 polynomial in (x−c)
  plus a RIGOROUS interval remainder; the containment law extended to
  FUNCTIONS (f(x) ∈ P(x−c) + rem for all x in the domain, tested).
  Constants admit orders 0..=169; the identity variable admits 1..=169.
  The maximum is numerical and operational, not cosmetic: elementary
  remainders use `1/(order+1)!`; `1/170!` remains normal in binary64 while
  `1/171!` is subnormal. The cap also limits live centered-domain powers to
  339 intervals, multiplication to 28,900 coefficient pairs, and one `exp` or
  `sin` composition to about 4.9 million pairs. Requests above the cap,
  unbounded domains, and non-finite scalar inputs are rejected before
  allocation or arithmetic.
  A degree-k coefficient's rounding interval is multiplied by the
  outward-rounded centered-domain power `(domain−c)^k` before absorption;
  flat absorption is not generally sound because it loses both degree and sign
  dependence. Sequential domain powers keep each coefficient-absorption pass
  linear in the emitted vector length (no linear-time claim is made for the
  polynomial multiplication itself).
  Elementary compositions
  (exp, sin) carry Lagrange remainders with declared-budget slack.
  Remainders shrink SUPERLINEARLY under subdivision (>20× per halving
  at order 5, tested) and beat plain interval excess by ≥1e6 on
  dependency-problem fixtures (x − x², tested). NOTE (measured, kept as
  documentation): expressions whose subterms are monotone and
  single-occurrence give interval arithmetic near-zero excess — TMs pay
  off exactly where the dependency problem lives.
  `bound_with_usefulness` preserves the full-domain enclosure and applies one
  caller-declared contextual width threshold. Non-finite width becomes
  `NoUsefulBound(LipschitzBlowup)`; a finite over-threshold enclosure retains
  the caller's closed cause and E09 reformulation.
- `newton::{newton_roots, krawczyk_step, RootBox, lipschitz_bound}` —
  certified roots: `Certified` ONLY under strict-interior contraction
  (existence + uniqueness); double roots come back `Possible`, never
  falsely certified (tested); empty Newton intersections certify
  ABSENCE. `lipschitz_bound` = outward-rounded derivative-enclosure
  magnitude (∞ when unbounded — never understated); the primitive the
  fs-geom certified-Lipschitz chart contract consumes.
  `RootSearchReport` also retains the evaluated-box width trajectory.
  `bound_with_usefulness` returns `None` for a complete no-candidate search,
  classifies a complete root hull contextually, and makes an incomplete search
  an absorbing `NoUsefulBound(BudgetExhausted)` even if the current hull is
  narrow.
- `AffineCtx`/`Affine` — affine forms `c₀ + Σ cᵢ·εᵢ + [−err, +err]`.
  Symbol identity IS correlation (same-context symbols cancel; fresh
  symbols don't). All rounding and the first-order mul residue
  (rad·rad) are absorbed into `err` with upward rounding. Reference
  operator impls (`&x + &y`), `scale`, `to_interval` (rigorous collapse),
  `radius`. Symbol ids are deterministic (context counter — replay-stable).

- `expansion` — Shewchuk floating-point expansion arithmetic (exact sums,
  scaled products, expansion×expansion products, exact-diff two-component
  constructors, sign extraction). For valid finite inputs whose result remains
  representable as a finite expansion, every op is error-free: the output's
  exact sum equals the exact real result (residual-law tested). Arithmetic
  panics rather than emit a non-finite component when a result leaves that
  domain.
  `two_diff`/`diff_expansion` reject non-finite operands and finite operands
  whose difference overflows rather than emit a NaN-bearing invalid expansion
  that downstream sign checks could mistake for exact zero. `expansion_sign`
  independently rejects every non-finite component, so invalid low-level
  arithmetic fails closed at the certificate decision boundary.
- `predicates` — adaptive-precision EXACT `orient2d`/`orient3d`/`incircle`/
  `insphere` (+ `*_with_stage` telemetry variants): fast float evaluation
  guarded by Shewchuk's proven stage-A error bounds (arrangement-specific
  constants with ε = 2⁻⁵³), escalating to exact expansion evaluation of
  the standard difference determinants (orient2d runs the full faithful
  B/C/D adaptive ladder). `orient2d_sos` and `orient3d_sos` implement the
  Edelsbrunner–Mücke tie-breaking ladder (total, antisymmetric, deterministic):
  2D term-by-term; 3D as the exact sign of the leading κ-coefficient of a
  moment-curve perturbation `pᵢ + κ·(sᵢ, sᵢ², sᵢ³)` (consistent — moment-curve
  points are always in general position), antisymmetry supplied by sort-by-index
  parity. Public predicate boundaries reject non-finite coordinates and panic
  before returning a sign if their filter or exact tail leaves the finite
  certified arithmetic domain.

## Invariants
1. **Containment law (G0)**: for every op and every point tuple inside the
   inputs, the true result lies inside the output. Tested: 100k arithmetic
   checks vs the dd oracle, 20k rewrite pairs (x·(y+z) vs x·y+x·z both
   contain truth), 50k intervals × 6 elementary functions vs an independent
   platform-libm point oracle, 20k point checks over 4000 random expression
   DAGs, Rump's polynomial containing the exact −54767/66192 (wide-and-
   right beats narrow-and-wrong — asserted width > 1 documents the honest
   outcome), and sign-complete finite-overflow cases for add/sub/mul/div that
   distinguish one-sided finite-real overflow from actual extended infinity.
2. Affine linear terms cancel EXACTLY on shared symbols: x − x collapses to
   width < 1e−13 (tightness ratio vs plain intervals > 1e10, tested);
   first-order identity-zero DAGs stay O(radius²) wide, independent of
   center magnitude.
3. Trig enclosures capture extrema: intervals straddling π/2 + kπ (sin) or
   kπ (cos) include ±1 exactly.
4. Enclosures never claim outside function ranges (sin/cos ⊆ [−1,1],
   tanh ⊆ [−1,1], exp ≥ 0, sqrt ≥ 0).
5. **Independent finite audit (G0)**: the isolated `tools/oracle` lane checks
   add/sub/mul/div, negation, absolute value, hull, intersection, certified
   width, and `exp`/`ln`/`sqrt`/`sin`/`cos`/`tanh` enclosures against exact-f64
   inputs evaluated by 256-bit MPFR. Its deterministic 4,096-generated-sample
   run per family passed as part of a 327,049-observation report with zero
   failures. The lane is development-only and absent from the production
   workspace dependency graph.
6. **Diverse-implementation cross-check (G0)**: `tools/cert-kernel`
   deterministically compares independently enclosed add/sub/mul/div/sqrt/
   exp/ln results with fs-ivl, retains per-family non-overlap counts and width
   ratio quantiles, and must detect a seeded one-ULP enclosure shrink. A green
   finite corpus establishes agreement only for those observations; it cannot
   promote either implementation to certificate authority or discharge the
   exhaustive corpus obligation. Any non-overlap fails the lane and is
   adjudicated with exact rational arithmetic where applicable, otherwise the
   256-bit MPFR oracle. Neither side is silently widened to manufacture
   agreement.

## Error model
Domain violations that admit NO enclosure panic with structured messages:
NaN/inverted constructor endpoints, sqrt of entirely-negative interval,
ln with hi ≤ 0. Partial-domain overlaps degrade gracefully (sqrt clips at
0; ln returns lower bound −∞). Division by zero-containing intervals
returns `WHOLE`, never panics. Indeterminate extended-real endpoint operations
such as `+∞ + −∞` and `0 × ∞` also return `WHOLE`. Representable operations
that merely overflow binary64 remain sound through one-sided infinite
enclosures rather than false singleton infinities.

Exact coordinate-difference expansion constructors panic on non-finite
operands and on finite differences outside the representable expansion
domain. Low-level expansion sum/product operations retain their documented
valid-finite-expansion input precondition and panic if an intermediate result
leaves the finite representable domain. Sign extraction panics if any component
is non-finite rather than treating NaN as exact zero.

Geometric predicates panic on non-finite coordinates and on detected
intermediate overflow. Inputs that underflow intermediate monomials remain
outside the certified domain as documented below; they do not acquire an exact
claim merely by producing finite machine values.

Interval root isolation accepts only finite domains, finite positive target
widths, and nonzero box budgets. `newton_roots_bounded` returns structured
`RootSearchError` values for invalid requests and a `RootSearchReport` whose
`complete` bit is true only when the subdivision stack was exhausted. On work
exhaustion, unevaluated regions remain `Possible`; they are never silently
dropped or promoted to certified. The compatibility `newton_roots` entry point
uses a fixed 65,536-box ceiling and panics only for invalid parameters.

Taylor construction and arithmetic return `TaylorModelError`. The error
distinguishes variable order below one, order above 169, non-finite domain,
non-finite scalar, bounded coefficient-reservation failure, and incompatible
domains. `variable`, `constant`, `scale`, `exp`, `sin`, and the reference
Add/Sub/Mul operators are fallible; operator output is `Result`. All arithmetic
borrows its inputs and builds a fresh output, so a refusal cannot partially
mutate or escape a model. Coefficient vectors contain at most 170 entries;
polynomial multiplication remains quadratic inside that admitted envelope.
Elementary-series coefficients and Lagrange terms use outward-rounded interval
reciprocal factorials rather than treating a rounded reciprocal as exact. If
finite admitted inputs overflow an interval intermediate, the affected Taylor
remainder/evaluation normalizes to `Interval::WHOLE`; tightness is lost, but a
one-sided infinity is never emitted as a false finite-real certificate.

## Determinism class
Bit-deterministic CROSS-ISA by construction (straight-line IEEE arithmetic
+ fs-math strict functions). Evidence: FNV-64 golden hash over 500 random
DAG enclosure endpoints + affine collapses = `0x3712_a4c1_2d5e_5864`,
recorded on aarch64-apple (M4 Pro), required to match on x86-64
(Threadripper) in tests/conformance.rs. Golden-evidence policy applies.
The Taylor/root battery separately pins `0x2cdd_5d22_d018_9466`, re-frozen on
remote x86-64 after reciprocal-factorial and exponential scalars became
outward intervals, and required to match on aarch64-apple. The same test
target, outside the hash payload, includes an executable double-double
wide-domain coefficient-rounding falsifier and path-specific wide-domain
containment cases.

## Cancellation behavior
Straight-line arithmetic needs no poll points. Root isolation is a bulk search
and therefore requires an explicit box budget; the report makes exhaustion
observable and a caller may resume at the returned `Possible` boxes.

## Unsafe boundary
None. `unsafe_code` denied; no capsules.

## Feature flags
None.

## Conformance tests
`tests/eft_interval_bridge_casebook.rs`: two fixed-order, replay-complete
Casebook records exercise the public BEDROCK seams directly. The first pins
the component bits produced by `fs_math::eft` and `fs_math::dd::Dd` against
the corresponding `fs_ivl::expansion` sum/product components. The second pins
`Interval::point` evaluations of exp, ln, and tanh to the deterministic
`fs_math` result nudged outward by exactly 3, 3, and 5 ULPs. A disclosed
seeded endpoint-reference corruption must produce one byte-replayable red
record and make the Casebook merge assertion refuse it.

`tests/structured_conformance_casebook.rs`: the cheap per-crate structured PR
entrypoint. Three fixed-order records cover interval/affine containment and
zero-divisor policy, bounded Taylor/root admission and incomplete-work
receipts, and exact-predicate signs, filtered stages, SoS ties, and
antisymmetry. A disclosed predicate-reference corruption proves that one
wrong expected sign produces one replay-identical red record and is rejected
by the Casebook merge assertion.

`tests/conformance.rs`: random-DAG containment battery + golden hash, plus 600
deterministic fs-propcheck interval-pair cases proving add/multiply pointwise
containment with shrinking enabled and the fixed cases retained.
`tests/taylor_battery.rs`: functional containment, subdivision convergence,
double-double-bracketed wide-domain coefficient rounding, certified-root
honesty, Lipschitz extraction, and a cross-ISA Taylor/root golden hash.
`tests/useful_bound.rs`: E09.2 G0/G3 driver integration, including finite
Taylor classification, overflow-to-Lipschitz refusal, and a budget-exhausted
root fixture whose deterministic width trajectory is retained to the refusal.
`tests/predicates.rs`: adversarial degeneracy batteries (cocircular /
cospherical / collinear lattice configurations), dyadic and 1-ulp
perturbation classification, SoS tie determinism, measured stage-A filter
rates. In-crate predicate batteries: i128 lattice oracle agreement (2D and
3D, thousands of rounds with degeneracy hits asserted), the Kettner-class
half-ulp grid where the naive float determinant demonstrably misjudges
while the exact ladder matches integer ground truth, expansion residual
laws, and SoS totality/antisymmetry.
In-crate: dd-oracle arithmetic containment, rewrite pairs, elementary
containment vs libm oracle, trig extrema capture, zero-divisor semantics,
positive/negative finite overflow for every basic operation, actual-infinity
separation, subnormal underflow and extreme cancellation, Rump's polynomial,
set operations, affine cancellation/tightness/independence, structured
rejections. Reimplementations must pass the golden-hash case bit-for-bit.

## No-claim boundaries
- **Rigor is conditional on fs-math's declared ULP budgets** (empirically
  enforced there and independently sampled by the isolated MPFR lane for the
  elementary functions used here): if a budget were violated, elementary
  enclosures could under-cover by the violation amount. The external run is a
  deterministic finite family, not an exhaustive proof over every binary64
  input or every interior point of every interval. Basic-op enclosures are
  unconditionally rigorous.
- `TaylorModel1` is univariate and fixed-order; multivariate Taylor models and
  sparse polynomial bases are not implemented. Public construction and
  arithmetic enforce the 169-order resource/numerics envelope and use fallible
  Vec reservation. Arithmetic overflow can conservatively return a whole-line
  enclosure rather than a useful bound. The API makes no throughput claim near
  the cap and does not expose a user-selected memory budget below the fixed
  envelope.
- Useful-bound projection proves only contextual classification of the
  retained enclosure. It does not authenticate the criterion, turn a
  `Possible` root box into an existence certificate, prove the width trajectory
  converges, or establish that the suggested E09 reformulation will succeed.
- Predicates are certified for inputs whose difference monomials (degree
  ≤ 5) stay inside the normal f64 range. Non-finite coordinates and detected
  overflow fail closed; intermediate underflow remains Shewchuk's inherited
  no-claim boundary.
- `incircle`/`insphere` symbolic perturbation is NOT provided (hooks and
  2D/3D orientation SoS only). `orient3d_sos` now IS a full Edelsbrunner–Mücke
  ladder (moment-curve perturbation, exact leading-coefficient sign) — the
  earlier projection cascade was NOT antisymmetric under swaps involving the
  4th point (bead wa8i V1) and has been replaced.
- Predicate throughput: only stage-A filter rates are measured; the exact
  tail is Vec-based and deliberately unoptimized (cold path).
- No quad-double (dd suffices for current oracles; recorded on 6ys.12).
- No SIMD lanes: scalar everywhere v1 (correctness-over-lanes, per plan).
- The structured EFT/interval bridge is a fixed, cheap G0 seam check. It does
  not expose fs-la's private residual accumulator, replace the full fs-ivl
  battery, supply an independent FrankenScipy oracle, or constitute fresh
  dual-ISA/full-G5 execution evidence.
- The per-crate structured entrypoint is a representative PR subset. The broad
  random-DAG/property batteries, Taylor convergence and certified-root
  studies, adversarial predicate suites, retained cross-ISA goldens, and any
  performance or fresh dual-ISA claims remain in their dedicated lanes.
- Certificate-regime routing follows `docs/CERTIFICATE_REGIMES.md`: interval
  root isolation and Taylor enclosures support declared finite-domain events,
  local exclusions, parameter bounds, balance audits, and finite-horizon tube
  ingredients. Repeated outward propagation is not a validated reachability
  tube, and an honest wide enclosure is `NoUsefulBound` for an engineering
  decision it cannot resolve rather than a certificate for that decision.
- Affine elementary functions (exp/sin on affine forms via Chebyshev
  linearization) are not implemented — interval fallback applies.
- `Interval::midpoint` of half-infinite intervals returns 0.0 (documented
  representative-point semantics, not an enclosure claim).
