# CONTRACT: fs-dfo

> Status: PARTIAL — CMA-ES (IGO form), BIPOP restarts, Nelder–Mead,
> balanced entropic OT, NSGA-II/III, MOEA/D, hypervolume helpers and
> archives, canonical sample-CVaR delegated to `fs-robust`, and
> discrete-support Wasserstein-DRO inner sup are in force; sep/low-rank CMA, NES
> parameterizations, DE, DIRECT, TR-DFO, and fs-exec population waves
> are recorded follow-up scope.

## Purpose and layer
Derivative-free optimization engines (plan §9.3). Layer: **L4 ASCENT**.
Deps: fs-rand (keyed sampling), fs-la (eigendecomposition), fs-math, fs-obs
(canonical conformance evidence), fs-blake3 (production root-bound callback and
full-study identities), and fs-robust (canonical empirical risk algebra).
Engines are IR-agnostic (closure objectives) BY DESIGN — routing from
the fs-opt problem IR is a wiring bead once that crate stabilizes
(deliberate collision avoidance, bead 7tv.4 trail).

## Public types and semantics
- `nsga3(...)` uses the versioned
  `NSGA3_NORMALIZATION_POLICY` variant
  `frankensim-scale-equilibrated-current-generation-asf-intercepts-v1`.
  This is a FrankenSim extension, not canonical Deb--Jain normalization:
  translated objectives are equilibrated by current considered-set rank-zero
  maxima before ASF extreme selection, and ideal/extreme/intercept state is not
  retained across generations. The public descriptor binds the ASF epsilon
  (`1e-6`), span floor (`1e-30`), pivot-ratio floor (`1e-12`),
  condition-error limit (`1e-8`), residual multiplier (`512`), 64-objective
  intercept cap, current-generation candidate/ideal/extreme policies,
  row-oriented column-scaled LU policy, exact fallback, retention, and
  non-finite-input no-claim strings. Its public `replay_identity()` is the one
  domain-separated typed policy identity consumed by downstream retained
  configurations: an exhaustive, no-`..` source destructure makes every added
  policy field a compile-time identity-schema decision, and field-by-field
  mutation coverage requires every current field to move the shared root.
- `cmaes(f, x0, CmaParams, seed) -> CmaReport` — full-covariance CMA-ES
  with the standard Hansen couplings (log-rank weights, rank-µ + rank-1
  covariance updates, cumulative step-size adaptation): the natural-
  gradient/IGO parameterization. Eigendecomposition via the landed
  cyclic Jacobi with symmetrization + eigenvalue flooring (SPD
  maintenance). Stagnation stops: TolX (σ·√λmax < 1e-12·σ₀) OR TolFun
  (no f_best improvement > 1e-12 relative for 120 generations) — the
  TolFun rule is what frees budget for restart ladders (measured).
- `admit_bipop(...) -> Result<BipopAdmission, BipopError>` — callback-free
  schema-v3 authority over one proposed BIPOP run. It checks the finite,
  nonempty root point,
  strictly positive finite sigma, finite optional target, positive hard budget,
  derived seed range, dense-matrix/population storage arithmetic, population
  ladder, local budgets, scheduler ordinal, and both Philox counter domains
  before scheduler allocation or objective invocation. The receipt retains the
  actual fs-rand stream-semantics version and all admitted maxima. Its nested
  fs-la `JacobiEighAdmission` is deliberately optional: it is `Some` exactly
  when `total_budget > base_lambda`, the first boundary at which a complete CMA
  generation (and hence eigendecomposition) is reachable, and is `None` for an
  initial-callback-only schedule.
- `try_bipop_cmaes(...) -> Result<BipopReport, BipopError>` — the fallible
  execution surface. It first obtains the exact `BipopAdmission`, then runs
  large/small alternation with population doubling and deterministic
  Philox-perturbed restart starts. Per-run allocation is
  `min(lambda*250, remaining)` and CMA consumes one initial callback followed
  only by complete lambda-sized generations. Raw-input/dependency/envelope
  refusals invoke no callback; a later generated-start or generated-candidate
  refusal may follow already completed work but never exposes the affected
  non-finite decision point to the objective.
- `bipop_cmaes(...) -> BipopReport` — compatibility projection of
  `try_bipop_cmaes`: finite targets become `Some`, the historical
  `f64::NEG_INFINITY` sentinel becomes `None`, and typed refusals panic at this
  legacy boundary.
  Its schema-v3 report exposes immutable canonical `BipopRootInputs` with a
  schema-v2 root identity, one
  `BipopEvaluationRecord` for every completed production callback, a
  schema-v3 length/domain-framed, derive-key `BipopTraceIdentity`, an immutable
  schema-versioned `BipopRestartRecord` for every run, the retained hard
  `total_budget`, the exact named `best_restart`, and legacy
  best/schedule/total-evaluation projections. Each record binds lane,
  population, local cap, derived seed, start, aggregate trace interval,
  `CmaStopReason`, and its complete `CmaReport`; `validate_ledger` returns a
  typed first-invariant refusal. Local CMA best selection and global restart
  selection both use strict `f64::total_cmp`, preserving the earliest exact-bit
  winner, including signed-zero and NaN payload ordering.
- `BipopReport::validate_study_identity()` recomputes the schema-v2 full-study
  derive-key BLAKE3 over retained and recomputed root/trace receipts, every
  compatibility projection, and every ordered restart/report bit. It
  distinguishes a stale payload from a valid payload that merely names a
  different externally retained study.
- `BIPOP_ROOT_IDENTITY_SCHEMA_DECLARATION`,
  `BIPOP_TRACE_IDENTITY_SCHEMA_DECLARATION`, and
  `BIPOP_STUDY_IDENTITY_SCHEMA_DECLARATION` are the owner-local field,
  dependency, mutation, transport, and version inventories for those three
  production identities. Explicit schema-input values make domains, versions,
  algorithm constants, and canonical order executable identity inputs rather
  than documentary constants. Their source guards exhaustively destructure the
  top-level sources and every nested identity-bearing root, evaluation row,
  restart record, CMA report, and trace receipt. Their mutation lanes move each
  declared schema input and nested retained-study field independently. The
  streaming `DomainHasher`, BLAKE3 derive-key context/material flags, complete
  compression closure, and `ReplayIdentity::canonical_bytes` accessor are
  implementation dependencies, not ambient library assumptions.
- Root identity v2 deliberately inserts `root-schema-version` into the typed
  canonical preimage and changes the kind to `fs-dfo-bipop-root-v2`. Root v1
  carried a public version constant that did not move its bytes; it is therefore
  not accepted as a compatibility alias. This migration intentionally re-keys
  dependent callback-trace and full-study values.
- Callback-trace v3 and full-study v2 deliberately migrate their typed 32-byte
  roots from plain BLAKE3 to the standard BLAKE3 derive-key construction. Trace
  v2 and study v1 are not silently re-keyed or accepted as compatibility aliases:
  their old digests name a different hash mode. The derive-key context is the
  complete versioned domain, while the same domain remains in the canonical
  payload as defense-in-depth. This separates typed identities from plain
  content hashes even for identical payload bytes and separates equal payloads
  across trace and study domains at the compression-mode/key level.
- The trace v3 canonical payload remains fully specified as little-endian
  `u64(domain_len) || domain || u32(schema_version) || u64(dimension) ||
  u64(row_count) || u64(root_len) || root`, followed in row order by each
  row's `u32(schema_version)`, `u64(restart)`, checked `u64(local_offset)`,
  `u64(objective.to_bits())`, and exactly `dimension` coordinate bit patterns.
  Full-study v2 retains the v1 ordered field grammar: every field is
  `u64(label_len) || label || u64(value_len) || value`; integers and float bit
  patterns are little-endian, booleans and enum tags are single bytes, and
  ordered collections bind explicit counts and indices. The in-band study
  domain and schema version are the first canonical pair.
- The enclosing `BipopReport` deliberately remains schema v3. Its record shape,
  root-input grammar, and report-level structural invariants did not change;
  nested identity validation already follows the independently versioned
  trace/study receipts and fails closed on their legacy versions. A report bump
  would unnecessarily re-key the separate FNV/`ReplayIdentity` root through its
  bound report version. The root kind, root schema v2 grammar, and FNV digest
  path are therefore unchanged by this migration.
- `BipopReport::admit_study_identity(expected)` is the cheap callback-free
  identity plus structural-ledger gate. It deliberately does not reconstruct
  candidate generation, final CMA state, or objective semantics.
- `BipopReport::admit_study_identity_with_replay(expected, oracle)` first runs
  that cheap gate, invoking `oracle` zero times on stale payload, ledger, or
  reference refusal. It then reruns the exact retained root through production
  BIPOP and requires the independently produced complete study identity to
  match. Replay execution and completed semantic mismatch have distinct typed
  errors. The oracle is executable authority and must implement the same
  deterministic objective semantics; its panics are not caught.
- `nelder_mead(...)` — deterministic simplex polish (no randomness).
- `wasserstein_worst_case(losses, costs, n, rho) -> DroReport` —
  exact dual evaluation for a discrete-support Wasserstein-DRO inner
  supremum, with deterministic dual minimization and a recovered
  support distribution. Kink cases split mass fractionally across
  active supports instead of pretending a single argmax distribution
  realizes the dual value.
- `empirical_cvar(samples, alpha) -> Result<EmpiricalCvarReport, RobustError>`
  — a direct re-export of the canonical `fs-robust` finite-sample risk algebra.
  The report carries overflow-safe CVaR, the deterministic lower empirical
  VaR/Rockafellar–Uryasev minimizer, boundary rank, and fractional boundary
  weight; DFO does not maintain a second order-statistic implementation.

- `steer` module (bead qlvf, lane a): WORLD-FORKING steering (P9).
  `SteeredStudy` = the deterministic (population, stream-index,
  weights) triple + a base seed + the steering LINEAGE; `fork` never
  mutates the parent and records a ledger-ready `SteerEvent`;
  `advance` is a pure function of (seed, stream-index, weights) so
  every branch replays bitwise from its lineage (`fingerprint` is the
  witness). Chance constraints live in `fs_uq::chance` — fs-uq sits
  above fs-dfo through fs-bo, so the integration points that way (the
  dependency cycle the first draft hit is the layer diagram talking).

## Invariants
1. DETERMINISM: the full evolution is a pure function of the seed
   (keyed Philox sampling, `total_cmp` ranking, lowest-index
   tie-breaks) — bitwise rerun-tested and cross-ISA golden-hashed.
2. IGO INVARIANCE: strictly monotone transforms of the objective give
   BITWISE-identical trajectories (tested with exp and cube transforms;
   precondition documented: monotonicity at the resolution of sampled
   values — x³ underflow and exp saturation break injectivity near
   machine-precision convergence, measured during bring-up).
3. Translation equivariance (behavioral, tested).
4. Covariance stays SPD (symmetrize + floor at each refresh).
5. BIPOP admission uses two independent restart bounds. With base population
   `b`, the maximum zero-based restart ordinal is
   `min(total_budget - 1, 8 + ceil((8 + 249*b*(2^8 - 1))/(b + 1)) + b)`.
   Nine large rungs (`0..=8`) may execute, but rung eight terminates the
   scheduler and cannot finance another small run; only rungs `0..=7` contribute
   to the pre-terminal large-spend proof. At most `b` further launches are
   possible once no more than `b` callbacks remain, and a zero-callback
   generated-start refusal terminates immediately; both cases are accounted for
   explicitly rather than hidden behind an unbounded budget surrogate.
6. A schema-v3 receipt also bounds `n*n` dense covariance entries and their
   addressable byte span, `n*max_large_lambda` population coordinates, the
   `lambda*250` local envelope, derived-seed range, restart-normal blocks
   `2*n*max_restart_ordinal`, and the maximum per-restart candidate-normal
   blocks `2*n*lambda*generations`. These are representability/counter-uniqueness
   proofs, not promises that the allocator has enough physical memory.
7. BIPOP records are ordinal and interval-contiguous; large/small lane choice
   follows retained cumulative spend; large populations double without
   overflow; each allocated cap is exactly `min(lambda*250, remaining)`; local
   `evals = 1 + generations*lambda` and do not exceed that cap; any numeric
   callback value at or below the finite target is a terminal witness even when
   a negative NaN is the total-order best; incomplete nonterminal prefixes are
   refused; every local `CmaReport` best and the named global best are the earliest
   `total_cmp` minima of their retained callback domains; and legacy
   best/schedule/total-evaluation fields are bit-exact projections of the
   ordered ledger. Validator reachability mirrors preflight, so an early target
   cannot retroactively bypass a Jacobi authority required by the retained hard
   budget.

## Error model
Direct `cmaes` retains its compatibility panic boundary for an empty dimension
or a non-finite generated decision. Convergence failure is DATA:
`converged: false` with best-found + diagnostics. Empirical CVaR returns the
canonical structured `RobustError` and does not panic.

`admit_bipop` and `try_bipop_cmaes` return the non-exhaustive `BipopError` for
malformed roots, dependency refusal, unrepresentable storage/population/budget/
seed/counter envelopes, generated non-finite starts or candidates, internal
accounting violations, and a generated ledger that fails validation. Preflight
refusals are callback-free; an execution-time refusal does not undo already
observed callbacks, but it returns no partial `BipopReport`. `bipop_cmaes`
intentionally projects those errors to a panic for source compatibility.
Objective panics and non-finite objective outputs are not caught or normalized
in this tranche.

`BipopReport::validate_ledger` returns `BipopLedgerError { restart, invariant }`
for schema, ordinal, population/lane, checked arithmetic, interval, exact local
allocation, dependency reachability, counter range, terminal completeness,
best-selection, root/trace identity, full-study identity, or compatibility-
projection violations. It authenticates the internally retained root and trace
against their recomputed identities, but an internally self-consistent digest is
not an external root of trust and structural validation does not reconstruct
every CMA transition or objective result.

`validate_study_identity` and cheap `admit_study_identity` return
`BipopStudyAdmissionError`, distinguishing stale payload identity, invalid
canonical payload, ledger refusal, and mismatch with a separately retained
reference. `admit_study_identity_with_replay` wraps those callback-free refusals
in `BipopReplayAdmissionError::EvidenceAdmission`; only after they pass can it
return `ReplayExecution` or `SemanticMismatch`. Identity admission is not a
signature scheme: the external reference and replay oracle remain caller-
supplied authorities.

## Determinism class
Bit-deterministic per seed, cross-ISA (golden hash
`0x5441_10a6_afb1_70a1`, bumped once when the TolFun stagnation
criterion was added — semantic justification recorded; verified
identical on both reference ISAs).

## Cancellation behavior
Single-threaded v1; population waves as fs-exec sibling scopes with
cancellation draining are the recorded follow-up (G4 scope there).

## Unsafe boundary
None.

## Feature flags
None.

## Conformance tests
tests/dfo_battery.rs (benchmarks incl. condition-1e6 ellipsoid,
determinism, invariance, BIPOP schedule, NM polish, golden hash);
tests/bipop_admission.rs (G0/G3/G4) proves callback-free malformed-input,
zero-budget, conservative derived-seed, and shared fs-la Jacobi-capability
refusal; exact schema-v3 receipt fields; Jacobi `None`/`Some` reachability at
the last callback-only and first complete-generation budgets; exact scheduler,
seed, and Philox-counter boundaries; one-callback terminal behavior; live
two-block `next_normal` consumption; containment of generated non-finite CMA
candidates and restart starts before the affected callback; typed target
semantics; and complete bitwise parity with the legacy compatibility spelling.
The `cma` G0 unit lane separately proves the exact dense-covariance byte and
`isize` address boundary, plus u128 counter-product overflow and first stream
coordinate reuse;
tests/bipop_restart_ledger.rs (G0/G3 typed restart-ledger partition, local and
aggregate budget edges, causal target/stagnation/budget terminal reasons,
earliest exact/signed-zero tie selection, within-restart signed-zero and mixed
NaN/finite total ordering, finite-target witnessing independent of a negative
NaN total-order best, root-bound complete production callback replay,
legacy-projection refusal, and large-small-large schedule replay);
tests/bipop_study_identity.rs (G3/G5 production full-study identity, same-seed
reference replay, callback-free stale/reference classification, replay-backed
same-objective admission, typed replay-execution refusal, and typed completed
semantic mismatch); the `cma` G0
unit lane independently
requires refusal of a post-ninth-large suffix, allocation above `lambda*250`, a
truncated nonterminal prefix, continuation after target convergence, early-
target evasion of reachable Jacobi authority, non-finite retained start/best/
sigma state, and an otherwise exact near-`usize::MAX` terminal record whose
next-generation arithmetic wraps. Its private G3/G5 mutation lane proves that
replay rejects correctly resealed non-best sigma, noninitial/nonwinning point,
nonwinning objective, and plausible-but-false terminal-reason mutations even
when cheap admission accepts their newly named identity;
tests/probe_tmp.rs (success-rate + stagnation-stop regression; filename
is bring-up history); tests/ot_battery.rs (4 cases): Sinkhorn
marginal feasibility ≤ 1e−8 and cost symmetry ≤ 1e−8; the ε-ladder
approaching the 1D monotone-coupling CLOSED FORM monotonically
(0.11 → 0.005 at ε = 0.05 → 0.002); translation covariance
(W₂² of equal translates = t² within 0.2%); OT golden
`0x58eb_8443_224c_a689`; tests/moo_battery.rs: hypervolume vs
hand-computed 2D/3D values including dominated/out-of-reference
degenerate cases; non-dominated-sort front assignment exact;
NSGA-II on ZDT1/ZDT2 at standard budgets (pop 80 × 200 generations —
short runs measurably leave the f2-minimal arm unexplored, documented
in the test) with mean front gap ≤ 0.05 (measured ≤ 0.0008), full f1
spread, hypervolume beating scrambled-Sobol random at MATCHED
evaluations (0.87 vs 0.19 / 0.54 vs 0.00), and bitwise replay; helper
edges for crowding; knee detection hitting a synthetic elbow exactly;
canonical empirical CVaR on 2·10⁵ Gaussian samples vs the closed form
μ + σφ(z_β)/(1−β) within 0.02 (and the lower RU minimizer matching the
VaR); fractional boundaries, tied integer boundaries, mixed-sign extremes,
constant `f64::MAX`, permutations, and direct `fs-robust` parity; MOO golden
hash `0x606f_35d4_bfb8_822a` remains frozen until the canonical accumulation is
run and any bit change is deliberately justified; MC hypervolume vs exact at
m = 2/3 within
0.01 absolute and m = 6 CLOSED FORMS beyond exact reach (single point
∏(ref − p) and two-point inclusion-exclusion, within 8% at 1.2·10⁵
Sobol samples, deterministic per seed with the hit count returned as
the honesty knob); public MC-hypervolume edges match exact
hypervolume's malformed-point ignore policy, empty dimension returns
zero, and zero samples are rejected; bounded archive — dominated
inserts are no-ops, dominating inserts evict, over-capacity eviction
removes the TRUE least contributor (verified against brute-force
keep-(k−1) subset enumeration); tests/nsga3_battery.rs (12 tests):
Das–Dennis counts
(C(p+m−1, m−1) exact: 91 at (3,12), 70 at (5,4)) and on-simplex
membership; public guards reject zero-division directions, invalid
NSGA-III reference directions, empty MOEA/D weights, and zero
neighborhoods at the API boundary; DTLZ2(m=3), the m=5 matched-budget
NSGA-III/NSGA-II comparison, the matched-budget MOEA/D comparison, and bitwise
replay remain executable gates. The earlier `0.0238`/98%, `6.51` vs `4.39`, and
`2.7457` vs `2.7775` observations belong to the pre-extension maxima-normalized
lane; they are historical context, not measurements of the named
scale-equilibrated v1 variant. Its policy-bound golden schema is v2 and
binds the version and root of the single typed normalization-policy identity;
the full-study replay configuration binds that same versioned root rather than
maintaining a second policy-field encoder. Caller mutation tests require both
retained preimages to move when the shared policy root moves. The golden
is frozen as `GOLDEN_HASH_V2=0x518f_9b30_a05c_d5e4` after the complete central
normalization, study-replay, WFG-family, and WFG4 selector union passed against
one source snapshot. The former `0xd912_6c49_f1b1_6897` hash remains historical
maxima-normalized evidence and is not relabeled as current evidence.
tests/wfg4_battery.rs adds the canonical normalized WFG4 fixture at
`M=3, k=4, l=20`, following the corrected WFG toolkit's `s_multi(30,10,0.35)`,
equal-weight `r_sum`, and concave shape construction. It composes that evaluator
with production NSGA-III at a closure-counted budget of 92 × (400 + 1) = 36,892
objective evaluations per run. The battery checks the analytic `[0,0,6]`
front extreme, the output-reconstructed scaled-sphere identity over deterministic
probes, a wrong-center transformation mutant, exact evaluation accounting, front
distance, Das-Dennis direction coverage, exact 3D hypervolume, and complete
ordered-front bitwise replay. Its v2 canonical configuration identity binds units,
seed plus stream kernel/tile, algorithm constants, the versioned typed child
root produced from every field of `NSGA3_NORMALIZATION_POLICY`, budgets, gates,
versions, and capabilities. A caller mutation test requires the campaign root
to move when that child root moves.
Separate original/replay front children bind every retained decision
and objective bit; the result child binds both roots, actual transform samples,
all metrics and verdicts, and the first differing component/bit when replay fails.
The validated object receipt and four `ConformanceCase` rows retain those roots,
while machine-neutral `BenchmarkResult` rows expose hypervolume and mean front
distance for later trend wiring. The formula source is Huband et al., *A Scalable
Multi-objective Test Problem Toolkit* (corrected EMO 2005 version), cross-read
against jMetal revision `ea7e882f6b8f94b99535921674e62cda7986f20e`.
At the retained seed and fixed 36,892-evaluation budget, the first
pre-extension central run measured mean front distance
`0.23171459361058497`; that observation is not a measurement of the v2 campaign
identity. The central v2 campaign measured mean front distance
`0.21508385736544294`, worst distance `0.313057630600937`, direction coverage
`86/91`, and hypervolume `62.62495408182373`; it replayed the complete ordered
front bit-for-bit. The coarse conformance ceiling remains `0.25`. These are
fixed-fixture observations, not a convergence-rate, broad WFG-family quality,
cross-ISA, or optimizer-performance claim. The budget and orthogonal gates
remain unchanged; the v2 receipt supersedes rather than relabels the earlier
receipt.
As in the older fs-dfo aggregate batteries, impossible public-API structural
failures such as an empty front remain ordinary Rust diagnostics before aggregate
emission; gate failures after a front exists retain structured red evidence.
tests/bipop_study_replay.rs adds one G5 study-scale replay receipt for the
production `bipop_cmaes` path on a finite four-dimensional shifted-Rastrigin
fixture. Two separately executed production runs at the recorded seed must
reproduce the actual schema-v3 `BipopAdmission`, every public `BipopReport`
field (including report schema, canonical root, retained `total_budget`, trace
identity, and full-study identity), every nested `CmaReport`, the complete
ordered production callback trace, and the canonical result frame. A separate
closure-owned callback witness is bit-compared against every production
evaluation schema/global-offset/restart/local-offset/point/objective field; it is
an independent observation path rather than a substitute for the production
trace. The fixture-v6 and result-v6 identities bind `best_restart()` and every
ordered `BipopRestartRecord`: schema/ordinal, large-or-small lane, population,
allocated budget, derived seed, exact start, half-open aggregate trace interval,
causal stop reason, and complete nested report. The fixture binds every public
admission field, actual and supported BIPOP report/evaluation/restart/trace/study,
fs-rand, and Jacobi schema or semantics
versions, explicit Jacobi presence, all nested Jacobi authority fields when
present, the supported fs-rand checkpoint version and checkpoint identity
domain, both strong-hash domain identifiers, target architecture/operating
system/pointer width/endianness, and full-width `u128` counter maxima as
high/low limbs. Under `SOURCE_FILE_IDENTITY_DOMAIN`, it also binds an ordered
path plus path-domain-separated BLAKE3 for exactly sixteen `include_bytes!` inputs:
`crates/fs-dfo/Cargo.toml`, `crates/fs-dfo/tests/bipop_study_replay.rs`,
`crates/fs-dfo/src/lib.rs`,
`crates/fs-dfo/src/cma.rs`, `crates/fs-la/src/lib.rs`,
`crates/fs-la/src/eigen.rs`, `crates/fs-rand/src/lib.rs`,
`crates/fs-rand/src/philox.rs`, `crates/fs-math/src/lib.rs`,
`crates/fs-math/src/det.rs`,
`crates/fs-math/src/dd.rs`, `crates/fs-math/src/eft.rs`,
`crates/fs-math/src/payne.rs`, `crates/fs-obs/src/ident.rs`,
`crates/fs-obs/src/lib.rs`, and `crates/fs-blake3/src/lib.rs`. This is the
relevant-source snapshot claimed by the study, not an identity for the whole
workspace or build toolchain. A fixed-study KAT independently requires
dimension 4, budget 6,000, base population 8,
maximum large population 2,048, local envelope 512,000, restart ordinal 5,999,
matrix/population entries 16/8,192, restart/CMA counter maxima 47,992/47,872,
and Jacobi work 76 under the 67,108,864-element authority cap.
An algebraically distinct shifted-Rastrigin oracle independently recomputes all
recorded objective values; it deliberately shares the deterministic cosine
primitive and therefore claims algebraic, not implementation-total,
independence under its identity-bound roundoff gate. The fixture identity
records units, initial point and shift, budget, target and improvement rule,
base population, large/small restart ledger and lambda ladder, whole-generation
budget semantics, candidate/global first-stable tie rules, logical
CMA/restart stream coordinates, semantic-oracle-v6 and stream-semantics versions,
dependency versions, capabilities, and the causal input seed. The semantic
gate runs an independent trace-driven CMA shadow from the root inputs: it
reconstructs every sample, distribution transition, sigma, terminal reason,
start, seed, interval, budget/evaluation/generation count, and nested report.
It uses live `StreamCheckpoint` witnesses to prove the exact restart formula
`2*n*(records-1)` and each CMA formula
`2*n*lambda*generations`. Both the restart and per-CMA witnesses require the
current `STREAM_CHECKPOINT_VERSION`, exact stream key, admitted stream-semantics
version, and exact terminal index; the fixture binds the supported checkpoint
version and `STREAM_CHECKPOINT_IDENTITY_DOMAIN`. Every realized count must
remain within the admitted schema-v3 cap. The gate also requires contiguous
intervals, exact legacy schedule/total/best projections,
report-budget/admission-budget equality, complete terminal history, and the
earliest `total_cmp` winner rather than an existentially compatible restart.

The result-v6 identity binds that fixture plus the complete report, canonical
root preimage, production trace and trace identity, full-study identity, restart
ledger, and separately labeled closure witness. Both production identity
validation and admission against the canonical production reference run before
the study can merge. `IdentityBuilder::child` retains a compact legacy 64-bit
FNV64 fixture root, so v6 additionally embeds the complete canonical fixture
bytes and their domain-separated fixture BLAKE3. The 256-bit result digest is
therefore transitively committed to every fixture byte rather than relying on
only the legacy child projection. In addition to the compact replay roots, the
test computes two domain-separated BLAKE3 families: the replay-identity domain
hashes the exact canonical fixture and result bytes, while a distinct
event-content domain hashes the canonical fs-obs content identities of the green
Custom receipt, green `ConformanceCase` verdict, and deterministic red
corruption event. Because the result embeds the complete fixture preimage, its
strong hash transitively binds the ordered sixteen-file source snapshot as well as
the target and admission provenance.

The Custom receipt passes structural fs-obs wire validation and exposes the
actual and supported authority versions, admitted and realized stream counters,
stream-position and checkpoint identity domains, actual and supported
checkpoint versions, both strong-hash domains, source-file identity domain and
count, complete Jacobi receipt, retained budget, record count, best restart,
target tuple, versions, and explicit no-claims. The green verdict also retains
the checkpoint version, target tuple, source-file count, and an explicit
compiler-fingerprint no-claim. Mutation lanes cover the receipt's BIPOP schema,
fs-rand semantics, checkpoint version/domain, total budget, Jacobi
presence/schema, admitted CMA cap, both realized-counter fields, replay/event/
source identity domains, and target architecture. Each mutant remains
structurally wire-valid under `fs_obs::validate_line` but stales the retained
event content identity.

**None of the five 256-bit BLAKE3 values is frozen yet.** The first central run
prints one `BIPOP_STUDY_FREEZE` line containing the fixture, result, green-receipt
event, green-verdict event, and deterministic-red-event digests. That exact
five-value tuple must be captured into explicit fixed-value assertions and the
same central test rerun successfully before any member is a regression KAT.
Until that freeze-and-reverification pass, all five are strong emitted sentinels
and replay diagnostics, not immutable goldens. Because fixture, receipt, and
verdict bind architecture, operating system, pointer width, and endianness, the
forthcoming tuple is target-specific; it must not be reused as a cross-target or
cross-ISA golden without separately captured, labeled evidence for that target.
It is also specific to the ordered sixteen-file source snapshot and the executing
test binary. Neither the source hashes nor the target tuple binds compiler
version/executable, flags, profile, or broader build identity; the fixture and
receipt explicitly record that compiler fingerprint as `not-bound-no-claim`.
The ordered snapshot includes the complete test harness itself, so weakening an
oracle, mutation, admission, or merge assertion moves the fixture identity.
Target-specific expected values live in the separate
`tests/support/bipop_study_replay_goldens.rs` data module and are deliberately
excluded from that snapshot: embedding an expected fixture digest in a source
preimage of the same fixture would require a cryptographic fixed point. The
snapshotted test root owns the comparison logic and treats `None` as explicit
unfrozen/no-claim state; the external data module is reviewed regression
authority, not a self-authenticating or signed root of trust.

A disclosed corruption seed drives stable one-field red evidence over a
returned coordinate. Its deterministic red-event detail retains both compact
and strong identities for the fixture, canonical result, and corrupted result,
alongside the exact bit mutation and refusal diagnostics. Independent mutation
lanes cover trace objective, schedule, retained total evaluations, admission
budget, Jacobi dimension and presence, a separately executed adjacent-budget
production report, complete ledger substitution, causal seed, and every
green-verdict envelope/payload field. Stale, correctly resealed,
retained-reference, and semantic-self-reference admission paths all fail closed.
This battery claims exact replay only for the same ordered relevant-source
snapshot, same target tuple, and same executing test binary on that bounded
fixture. It does not claim compiler/build identity, optimizer quality,
arbitrary-input admission success, allocation recovery, retained production
callback payloads after a failed execution, other objectives/dimensions/budgets/seeds, refreshed
cross-target or cross-ISA equality, cancellation, checkpoint/resume,
authenticated or signed evidence, or performance.
tests/nelder_mead_study_replay.rs adds the corresponding standalone-family
receipt for production `nelder_mead`. It reuses the fixed two-dimensional
Rosenbrock leg from the crate golden (`x0 = [0.3, -0.2]`, simplex scale `0.2`,
soft maximum `2,000`, impossible target `-1.0`) and retains every ordered
objective callback plus all three public tuple outputs. The accounting gate
reconstructs the initial simplex, recomputes every objective bit, requires the
reported evaluation count to equal the closure count, admits only the bounded
overshoot of one already-started Nelder–Mead transition, and proves the returned
best is the global minimum of the retained evaluated trace. Separate canonical
fixture and result identities bind the units, arguments, standard coefficients,
ordering/budget semantics, dependency versions, complete callback stream, and
returned bits; independent runs must reproduce the complete result frame. A
disclosed `StreamKey` selects one returned coordinate and low mantissa bit. The
test first proves the unsealed edit fails payload validation, then reseals it,
proves exactly that bit and no callback or other output field changed, emits two
byte-identical wire-valid red `ConformanceCase` events, and catches the local
merge gate refusing the typed reference-identity mismatch. This is fixed-input
same-build evidence: the mutation seed is not an optimizer seed, the maximum is
a soft loop-entry budget rather than a strict callback cap, and the receipt adds
no optimizer-quality, broad-input, refreshed cross-ISA, cancellation,
checkpoint, authenticated-ledger, or performance claim.
tests/moead_study_replay.rs adds a full-study receipt for production `moead` on
a short four-variable ZDT1 fixture with eight ordered Das-Dennis weights and
four generations. The fixture identity binds the dimension and bounds, every
`MoeadParams` field, every weight bit, dependency and stream-semantics versions,
and the optimizer's logical stream coordinates. The result identity binds all
forty ordered objective-callback decision/objective bits and every ordered
returned-front `Individual` decision/objective bit. Its independent accounting
gate regenerates the complete initializer from the recorded stream, recomputes
ZDT1 exactly, checks callback count, dimensions and box membership, and requires
each returned individual to be an evaluated, mutually non-dominated point.
Independent production runs must reproduce the complete canonical frame. A
disclosed evidence-generator `StreamKey`, separate from the optimizer's recorded
study seed, selects one returned objective and low mantissa bit; the test proves
the stale payload is rejected, reseals it, proves that sole bit delta, emits
identical wire-valid red fs-obs evidence, and catches the local merge gate
refusing the typed retained-reference mismatch.
tests/dro_battery.rs (4 cases): the one-sample kink LP recovers the
fractional q = [0.5, 0.5] and worst-case value 0.5; tiny-scale kinks
use scale-relative recovery rather than an absolute lambda cutoff;
large radius saturates at max loss; public guards reject
empty/non-finite losses, zero sample counts, invalid radii/costs, and
rows without a zero-cost stay-put support. Each deterministic case emits
its passing aggregate outcome as a canonical fs-obs `ConformanceCase` with
seed zero after its direct assertions; earlier failures remain ordinary Rust
test diagnostics. tests/dro_oracle_battery.rs
(4 cases): closed-form endpoints and monotonicity in rho, exact
tiny-LP strong-duality oracle, robust-decision shift demo, and frozen
DRO golden hash `0xd21c_d092_b4a5_ba98`. tests/steer.rs has two logged
conformance cases plus one empty-population guard regression: a fork leaves
its parent untouched while both branches replay bitwise, opposite weights
steer sibling forks toward opposite Pareto extremes, and zero population is
rejected at construction. Aggregate results use canonical fs-obs
`ConformanceCase` events with their exact input seeds; fork and steering
measurements use separate validated fs-obs `Custom` events that also retain
those seeds. Assertions precede each passing aggregate verdict, so earlier
failures remain ordinary Rust test diagnostics.

The remaining MOO, NSGA-III, OT, and DRO-oracle batteries use the same canonical
evidence boundary. Their 22 ordinary post-assert rows retain the existing suite
and case identities and emit `ConformanceCase` events with Info/Error severity,
failure-record linting, `to_jsonl`, schema validation, print-before-terminal-
assert ordering, and exact input-seed provenance. The cases are:

- `fs-dfo-moo`: `hypervolume`, `nds`, `helper-edges`, `zdt1`, `zdt2`,
  `nsga2-replay`, `knee`, `cvar-ru`, `cvar-ties`, `mc-hv`, `mc-hv-edges`, and
  `hv-archive`.
- `fs-dfo-nsga3`: `das-dennis`, `dtlz2-m3`, `m5-vs-nsga2`, and `moead`.
- `fs-dfo-ot`: `marginals-symmetry`, `eps-ladder`, and `translation`.
- `fs-dfo-dro-oracle`: `endpoints`, `strong-duality`, and `robust-shift`.

The four frozen-hash rows are measurements rather than pre-declared passing
verdicts. `moo-golden`, `nsga3-golden`, `ot-golden`, and `dro-golden` therefore
emit validated object-shaped `Custom` companions under distinct
`<case>/measurement` scopes, with the original case identity retained as the
Custom name. Each companion records actual and expected hashes plus the exact
input roots and logical substream coordinates, then the original frozen-hash
assertion runs unchanged. This preserves the diagnostic row on a golden
mismatch without falsely emitting a passing aggregate.

Seed provenance follows the literal fixtures. MOO uses optimizer input seed 21
and Sobol comparator seed 777 for the ZDT summaries, seed 21 for replay, seed
101 with stream `(kernel=0xC7A2,tile=0)` for Gaussian CVaR, seeds 5 and 6 with
CVaR tile 1 for its golden, MC-HV roots 7, 8, 9, 10, and 42, and edge-call roots
11 and 12. NSGA-III uses seed 17 for DTLZ2, optimizer seed 23 plus MC-HV seed 99
for the many-objective comparison, seed 3 for its golden, and roots 29 and 31
for the composite MOEA/D summary. OT uses root 111 and kernel `0x0007`, with
tiles 1/2 for marginal symmetry, 3/4 for the epsilon ladder, 5 for translation,
and 6/7 for its golden. DRO-oracle uses root 131 with
`(kernel=0x0D20,tile=1)`, root 132 with kernel `0x0D21` and tiles 0..19, and
root 133 with `(kernel=0x0D22,tile=0)`; robust-shift is fixed input.

Fixed-input and multi-root aggregate summaries use seed zero, with every
subordinate root named in the typed detail. `fs_rand::StreamKey` kernel and tile
values are logical input-substream coordinates, never execution seeds. Direct
assertions and expectations remain before ordinary aggregate emission; the ZDT
loop can still emit its first completed case before a later case fails. Silent
parity, guard, and `should_panic` tests remain silent.

The core DFO and success-rate batteries complete the same migration under the
existing `fs-dfo` suite identity. Four ordinary post-assert rows use canonical
`ConformanceCase` events, failure-record linting, JSONL serialization, schema
validation, and their original case identities: `benchmarks`, `igo-invariance`,
`bipop`, and `success-rate`. The `benchmarks` composite carries aggregate seed
zero and names the CMA-ES input roots 1, 2, and 3. `igo-invariance` carries root
7, shared by the plain, exponential, and cubic objective runs. `bipop` carries
root 17. `success-rate` carries aggregate seed zero and names the five optimizer
input roots 1 through 5; its per-seed loop, majority threshold, and failed-run
stagnation bound are unchanged.

The frozen `dfo-golden` info row is not promoted to a passing aggregate before
its hash assertion. It is a validated object-shaped `Custom` measurement with
name `dfo-golden` under the distinct `dfo-golden/measurement` scope. The payload
records actual and expected hashes, aggregate input seed zero, CMA-ES input
roots 99 and 100, the fixed-input Nelder-Mead leg, and a null execution seed.
Those roots are optimizer inputs, not scheduler/execution provenance; these
tests do not create an asupersync execution context or `fs_rand::StreamKey`
substreams. The measurement remains before the original terminal hash assertion
so mismatches
retain their diagnostic row without emitting a false passing verdict. The
deterministic-evolution, translation-equivariance, and Nelder-Mead replay tests
remain assertion-only and silent.

## No-claim boundaries
- No published-ERT-table parity claims yet (in-repo BBOB-class fixtures
  only; the external COCO battery is follow-up).
- Schema-v3 `BipopAdmission` certifies the checked formulas and dependency
  authority represented in its fields for the supplied root inputs. It is not
  an authenticated root-input identity, a physical-memory-availability promise,
  a callback-panic boundary, an objective-output-finiteness guarantee, or a
  production cancellation/checkpoint capability. The versioned
  `StreamCheckpoint` in the G5 shadow is a counter-consumption witness under the
  recorded checkpoint domain, not evidence that BIPOP execution can pause and
  resume.
- A successful schema-v3 `BipopReport` retains the complete production callback
  trace and strong identities. Cheap study admission is content/structure
  authentication against a caller-retained reference; it is neither proof that
  an untrusted reference is genuine nor objective-semantic replay. The replay-
  backed method supplies the latter only relative to its caller-provided
  executable oracle and same-build deterministic semantics. Neither method is a
  signature or independent implementation oracle.
- Generated-candidate/start refusal and objective panic still publish no partial
  `BipopReport`; `Cx` cancellation with request/drain/finalize, no-half-record
  prefix evidence, pause/resume/fork identity, and interval-sheaf/schedule
  certificates remain parent bead `7tv.23` obligations.
- Root/trace/study owner declarations and authority pins cover the current
  production encoders, explicit schema inputs, nested source fields, and exact
  digest/accessor implementation closure at the source-static level. Generated
  identity-catalog/coupling regeneration and batch execution of every mutation
  lane remain proof-pending under parent `7tv.23`; source declarations alone are
  not executable proof. The declarations do not add a signature, independent
  objective oracle, whole-build identity, or proof beyond the callback-free and
  replay-relative admission boundaries above.
- The four `fs-dfo` aggregate rows attest only to their completed in-repo
  assertions at the recorded input roots. They do not promote the frozen hash
  measurement to a verdict or claim coverage of other seeds, external benchmark
  corpora, scheduler determinism, or unmeasured success probabilities.
- The standalone Nelder–Mead study receipt covers one fixed two-dimensional
  same-build run. Its mutation seed belongs only to the red evidence generator;
  the optimizer itself has no randomness. The reported maximum is a soft
  loop-entry budget and may be exceeded only while the already-started
  reflection/expansion/contraction/shrink transition drains. The fixture makes
  no strict-budget, quality, all-objective/dimension/configuration, refreshed
  cross-ISA, `Cx`, checkpoint, authenticated-ledger, or performance claim.
- The standalone MOEA/D study receipt covers one same-build four-dimensional
  ZDT1 fixture at one recorded configuration and seed. It exposes only objective
  callbacks and the returned rank-zero subset, not the final population, ideal
  history, neighborhoods, or replacement decisions. It makes no convergence,
  front-quality, hypervolume, coverage, diversity, optimizer-superiority,
  all-objective/dimension/configuration/seed, refreshed cross-ISA, `Cx`, checkpoint,
  parallelism, authenticated-ledger, external-oracle, or performance claim.
- NSGA-III ordinary normalization is the named FrankenSim
  scale-equilibrated/current-generation extension: the ideal is the
  componentwise minimum over the current environmental-selection considered
  set, current rank-zero translated maxima equilibrate the ASF, distinct ASF
  extremes form the rows of a column-scaled hyperplane system, and positive
  finite intercepts become association spans. The same rank-zero maxima floored
  at `1e-30` become the final association spans only on exact refusal for
  duplicate extremes, unsupported objective count, non-finite inputs,
  singular/ill-conditioned/tiny-pivot/high-residual solves, or invalid
  coefficients/intercepts. This is not canonical
  Deb--Jain normalization and makes no cross-generation ideal/extreme retention
  claim. On non-finite objective input the fallback is a deterministic return
  path, not authority that its spans or downstream distances are finite or
  usable.
- The WFG4 battery claims only the normalized tri-objective `k=4, l=20`
  in-repository fixture and its recorded seed/budget. It does not claim the
  complete WFG suite, external COCO parity, optimizer performance, refreshed
  cross-ISA execution, or asupersync cancellation coverage. The source cross-read
  is not an executable independent WFG/jMetal oracle comparison.
- `mc_hypervolume` is the m > 4 path; its accuracy knob is the
  sample count (standard error √(p(1−p)/n)) — no silent precision
  claim. `HvArchive` eviction uses EXACT contributions (m ≤ 4);
  MC-contribution eviction joins with its many-objective consumer.
- Module `ot`: BALANCED entropic OT only (equal masses asserted);
  unbalanced/partial transport and Sinkhorn divergences (debiasing)
  are follow-up lanes.
- Module `dro`: discrete candidate support only; continuous ambiguity
  sets, adaptive support generation, and coupled decision-dependent
  losses are follow-up scope.
- MOO scope (module `moo`): NSGA-II/III, MOEA/D, exact hypervolume
  m ≤ 4, MC hypervolume beyond m = 4, bounded hypervolume archive,
  and knee detection. Canonical sample-CVaR is re-exported from `fs-robust` at
  the crate root. Gradient-based Pareto tracing (fs-ascent
  continuation), ledger world-forking steering, and chance constraints
  are still split lanes.
- Sep-CMA/low-rank (dim > ~200), NES, DE, DIRECT, TR-DFO: not built.
- No constraint handling (fs-constraint owns kinds; integration later).
- No parallel evaluation waves yet (fs-exec bead).
