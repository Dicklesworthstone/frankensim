# CONTRACT: fs-spectral

Layer L1 spectral semantics and health monitoring. This crate owns three related,
but deliberately separate, capabilities:

1. deterministic sheaf-Laplacian gap health, confidence demotion, conditioning
   composition, and conditioning-aware route scoring;
2. RB.1a versioned spectral problem admission plus proposition-bound,
   set-valued result truth; and
3. RB.M1 versioned Maslov--Krein--Evans theorem statements, executable
   convention transforms, and explicit hypothesis/implication lattices.

The latter two capabilities are semantic and epistemic boundaries. They
classify problems or theorem statements, verify that retained evidence names
the exact proposition being used, check explicit prerequisites, and validate
result or lattice claims. They do not compute eigenpairs, prove theorems, or
manufacture scientific evidence.

## Dependencies

Production dependencies are limited to:

- `fs-blake3` for domain-separated canonical identities and the
  presented -> verified -> admitted authority typestate;
- `fs-evidence` for the existing gap-health confidence color lattice;
- `fs-exec` for the cancellation contexts polled by eigensolver-service ticks
  and bridge-lattice validation;
- `fs-la` for the deterministic Lanczos/LOBPCG backends the service
  wraps; and
- `fs-qty` for compile-time Floquet quantities and runtime-dimensioned
  normalization boundaries.

There is no FFI and no external numerical runtime.

## Existing spectral-health API

- `symmetric_eigenvalues(&[Vec<f64>])` is a deterministic small dense
  symmetric Jacobi eigensolver. It rejects empty, non-square, and
  non-symmetric inputs.
- `spectral_gap(&[f64])` returns the gap above the smallest eigenvalue,
  spectral spread, and a dimensionless health ratio.
- `GapHealthMonitor` applies a hysteresis band so marginal regions do not
  flap between `Healthy` and `Degraded`.
- `propagate` can only preserve or demote evidence color when the gap
  degrades; it never promotes confidence.
- `compose_conditioning` multiplies finite nonnegative amplification factors.
- `route` minimizes base cost plus a logarithmic conditioning penalty with
  deterministic tie behavior.

This pre-RB.1 numeric API is an unauthoritative monitoring heuristic. Its
`SpectralGap`, `Health`, and router outputs cannot construct an admitted
problem witness, method token, gap-semantics token, cluster proposition, or
`SpectralTruthV1`. In particular, its historical zero-spread behavior is not a
scientific truth input. RB.1b owns the byte-visible correction and permanent
`[a,a]` regression; RB.1a isolates it behind this explicit no-claim boundary.

## Admission API and semantics

### Product problem schema

`SpectralProblemSpecV1` keeps independent axes independent:

- real or complex scalar field;
- standard linear, generalized pencil, or exact-grade matrix polynomial
  representation;
- ordinary or descriptor role with explicit infinity policy;
- direct, monodromy/Floquet, or analytic operator-function origin;
- typed metric/form-supported structure propositions;
- complete scaling-map identities and runtime spectral dimensions;
- domain/codomain metrics, gauge convention, and zero-padding convention;
- a product set of regularity propositions;
- ordering/target convention; and
- requested candidate, partial, region, or full-finite coverage.

`validate_problem` is the only constructor for
`ValidatedSpectralProblemV1`. It canonicalizes order, returns every detected
issue that can be safely evaluated within the current bounded validation
phase, and mints no partial token on failure. A resource-cap violation returns
immediately before sorting, hashing, or pairwise inspection, so hostile
oversized input cannot force deeper work merely to accumulate secondary
diagnostics.

## Physical-domain adapter API (RB.8a)

The `adapter` module is the one-way boundary by which rotordynamics, control,
acoustic, electromagnetic, and periodic-domain crates describe an induced
spectral problem. Dependency direction remains downward: domain crates may
construct adapter descriptors using `fs-spectral` types, while `fs-spectral`
does not import, wrap, mutate, or assign scientific authority to their model
objects.

`SpectralAdapterSpecV1` binds all of the following to the exact
`SpectralProblemId` of a complete `ValidatedSpectralProblemV1`:

- source artifact and immutable physical-model version;
- source operator class plus primal/state and dual/test space identities and
  dimensions;
- unit/scaling, frame, metric, and norm crosswalks with named map artifacts;
- constraints, nullspace/gauge data, parameter/held-variable schema, and
  boundary-condition schema;
- a model-version-bound frozen linearization point or an explicit
  non-applicability artifact;
- periodic phase, Poincare section, or event-word semantics;
- retained source structure evidence when the target problem carries
  structure claims;
- exact forward fidelity and an explicit reverse-interpretation/no-claim
  boundary; and
- source-to-spectral QoI crosswalks.

Every semantic component is retained, explicitly non-applicable with a
content-addressed justification, or unresolved. `Unknown` never admits.
Descriptor targets require retained constraint data. Monodromy/Floquet targets
require a periodic source and a retained phase/section map. An alleged identity
frame map must name the same source and target frame. A frozen linearization
must name the current source-model version. Target space dimensions, metric
IDs, scaling ID, representation, descriptor role, origin, and problem identity
are cross-checked against the sealed target token.

`validate_adapter_v1` is the sole constructor for
`ValidatedSpectralAdapterV1`. It refuses lossy or ambiguous mappings, unresolved
bindings, stale linearizations, missing descriptor constraints, phase/origin
mismatches, unsupported structure, duplicate QoI interpretations, and all
target mismatches. Exact quotient fidelity must name the same retained
nullspace artifact, and neither a one-way nor a quotient adapter may claim an
exact whole-mode reverse interpretation. The QoI resource cap is enforced
before sorting or identity work. Successful validation canonicalizes QoI order and returns a
domain-separated `SpectralAdapterIdV1` plus its canonical-preimage receipt, so
schema replay is deterministic and changes to any bound model, map, frame,
convention, target, or no-claim artifact are byte-visible.

The adapter receipt proves only stable schema identity and successful local
cross-checks. Opaque source artifact, model, map, witness, and no-claim IDs are
references to upstream evidence; choosing their bytes does not prove physical
correctness, unit-transform correctness, model fidelity, or a theorem. Exact
one-way and quotient adapters deliberately carry an explicit absent-inverse
statement. No caller may reinterpret one as an inverse, lift, reconstruction,
mode-shape theorem, stability theorem, or completeness certificate. Numerical
spectral truth remains governed separately by the admission and truth APIs.

G0/G3/G5 adapter tests cover canonical QoI ordering and replay, identity
sensitivity to model/frame changes, wrong metric and scaling targets, identity
frame mismatch, lossy frames and forward maps, stale linearization, omitted
boundaries, duplicate and oversized QoIs, descriptor constraints,
periodic-phase mismatch, and source-structure retention.

## Maslov--Krein--Evans bridge statement API (RB.M1)

The `bridge` module is a theorem-statement and review boundary, not a proof
engine. It keeps `MaslovIndexObjectV1`, `KreinSpectralFlowObjectV1`,
`EvansWindingObjectV1`, and `MachineInstabilityCountObjectV1` as nonfungible
types. A label or coincident integer cannot turn one into another. A theorem
node must instead carry exact domain, versioned operator family, parameter
path/direction, count transforms, endpoint convention, multiplicity rule, and
correspondence-map identities.

`BridgeTheoremScopeV1` retains the implication lattice from the classical
finite Hamiltonian case through periodic monodromy and spatial-dynamics Evans
extensions to the bold maximal triple equality. Weaker nodes remain explicit
targets of acyclic, content-addressed projection edges. A separate machine
corollary requires a distinct physical-instability object and interpretation
map; spectral instability is never physical instability by enum relabeling.

Each node carries a product set of typed hypotheses. The validator requires
the hypotheses applicable to every count in its conclusion, including:

- symplectic path, Lagrangian Fredholm-pair, crossing-form, and endpoint
  obligations for Maslov counts;
- continuous Pontryagin path, nondegenerate Krein form, and an explicit
  neutral-signature policy for Krein flow;
- analytic Fredholm family, analytic domain/contour, essential-spectrum
  exclusion, and Evans normalization for winding counts;
- exact parameter direction and multiplicity semantics for every equality;
- monodromy, spatial-dichotomy, pairwise-correspondence, and physical
  interpretation obligations at the scopes where they apply.

A hypothesis may be witness-referenced, explicitly unresolved, or refuted.
Those are statement states, not truth constructors. Likewise,
`BridgeProofStateV1::ProofArtifactReferenced` records proof, verifier, formal
system, policy/TCB, and no-claim identities but yields only
`ReferencedNotVerified`. `ValidatedBridgeLatticeV1::scientific_authority()` is
always `ScientificCorrectnessNotProven`. This preserves the maximal coherent
conjecture without laundering an opaque digest into a theorem.
Referenced node verifiers and formal-system versions must match the lattice's
declared TCB exactly.
An unresolved neutral-Krein policy cannot be paired with a witnessed
neutral-closure hypothesis; it must remain explicitly unresolved or refuted.

Endpoint half-signatures use `DoubledSignedCountV1`, avoiding floating-point
convention ambiguity. `derive_convention_transform_v1` computes the exact
orientation sign and endpoint correction, refuses unresolved endpoint
signatures, maps reversed initial/final rules back onto canonical positive
left/right endpoints, and uses checked integer arithmetic. Periodic conclusions
that include the Maslov object require an even-dimensional symplectic state
space. Every theorem node must also preregister the falsifiers relevant to its
count product: multiple/degenerate
and tangential crossings, endpoint crossings, neutral Krein type,
non-Fredholmness, contour and essential-spectrum contact, orientation changes,
multiplicity changes, and reviewer statement mutation.

`validate_bridge_lattice_v1` checks caller budgets beneath hard schema caps
before sorting or graph work, polls an execution `Cx` initially and at bounded
node/edge/topological-sort strides, canonicalizes all set-valued inputs, rejects
duplicate or cyclic implications, and returns a domain-separated identity plus
canonical-preimage receipt. Statement version, count objects, conventions,
hypothesis/proof states, falsifiers, implication maps, TCB, and budget are all
byte-visible identity inputs. Implication edges cannot run from a lower
extension rank to a higher one. Each non-refuted maximal node must directly
project to distinct non-refuted periodic and spatial branch roots, and each
branch root must reach a non-refuted classical finite node without traversing
refuted theorem nodes or implications; machine-corollary edges cannot
substitute for that spectral coverage.

G0/G3/G5 bridge tests cover canonical replay under node/hypothesis/falsifier
permutation, orientation and endpoint transformations, unresolved endpoint
refusal, required hypothesis and falsifier gates, neutral/tangential/essential
contact preregistration, machine-count separation, implication cycles, missing
weaker scopes, lower-to-higher edges, direct projection coverage, retained
refutations, schema/convention/reviewer mutation, resource caps, cancellation,
and the fixed no-theorem authority boundary.

### Authority boundary

A generic verifier/admitter success is policy-relative and cannot directly
create a favorable proposition. Evidence producers must:

1. build the family-specific canonical proposition receipt;
2. present that receipt with an external anchor, verifier identity, and policy
   identity;
3. pass an injected `AuthorityVerifier`;
4. separately pass an injected `AuthorityAdmitter`; and
5. configure a separate `SpectralPromotionTrustRootV1` from independently
   retained full verifier and policy receipts;
6. have that root re-adjudicate both typed identities and their canonical-byte
   observations to mint an opaque `SpectralPromotionWitnessV1`; and
7. pair that opaque decision with the exact same policy-relative authority to
   obtain `AdmittedSpectralWitnessV1`.

The final pairing refuses any difference in the complete subject receipt,
anchor, verifier identity, key-policy identity, or fixed spectral promotion
context. A one-argument conversion from generic `Admitted` no longer exists.
The retained promotion audit is bounded to verifier/policy namespaces,
canonical-byte roots and lengths, and the fixed context. Its versioned witness
encoding binds those fields in addition to the pre-existing subject, anchor,
verifier, and policy identities.

Every consuming validator recomputes the expected typed proposition ID and its
independent canonical-preimage root and byte length. All must match. The audit
record retains anchor, verifier, policy, trust state, and the explicit
`ScientificCorrectnessNotProven` boundary. Policy admission and root matching
are not a theorem.

No-claim boundary: fs-blake3 v1 currently makes the promotion witness fields
opaque but allows public, copyable trust-root configuration. Therefore this
consumer hardening prevents direct and mismatched escalation, but it does not
authenticate root-instance ownership: code that can reproduce the exact root
configuration can reproduce its decision. Do not interpret a promotion token
as cryptographic caller identity or scientific correctness; a branded or
authenticated domain-owner capability is still required for that stronger
claim.

Structure and regularity receipts bind subject, scalar field, representation,
descriptor role, origin, the complete scaling context (identity, dimensions,
scale, and every map identity), domain/codomain identities and dimensions, and
the exact claim payload. Gauge receipts additionally bind the exact
`SpectralGaugeArtifactId` or `SpectralQuotientMapId`. Zero-padding receipts bind
the same witness-free gauge/reduction context, so a padding witness produced
for quotient map A cannot be replayed under quotient map B even when their
nullities and serialized counts coincide. Metric, gauge, zero-padding,
cluster, multiplicity, separation, and completeness families have distinct
canonical payload dialects. Relabeling one family or changing one semantic
axis changes the proposition identity.

The raw descriptor and canonical claim slices exposed by a validated problem
are observational views, not detachable authorization capabilities.
Authority-bearing consumers accept the complete
`ValidatedSpectralProblemV1`, preserving its problem identity and all
cross-axis validation.

`SpectralSubjectId` identifies the actual induced operator whose spectrum is
requested; metric and form IDs identify the actual supporting artifacts.
Changing that induced operator, metric, or form requires a new ID. Structure
and regularity evidence is intentionally about the induced subject and its
post-space IDs, so two reduction lineages that provably yield that same subject
may reuse it. A proposition specifically about descent through one reduction
map must instead bind `GaugeContextV1`; gauge and zero-serialization families
already do so.

### Metrics, forms, and exact structure

Metrics, symplectic forms, Krein forms, conjugations, and form-free
propositions are nonfungible typed supports even if their digest bytes happen
to coincide.

`MetricDefinitenessV1::Euclidean` is witness-free only for the unique
dimension-derived ID produced by `SpectralMetricV1::euclidean`. Caller-chosen
IDs cannot acquire positive definiteness by selecting the Euclidean variant.
Other positive-definite, indefinite, and singular metric claims require an
admitted exact proposition.

Approximate structure claims remain representable for diagnostics and future
budget-aware routing. V1 structure-preserving method tokens require a
zero-tolerance witnessed proposition. A contradiction on the selected support
fails closed. Because zero defect denotes exact structure rather than a
norm-specific error budget, an exact witnessed property conflicts with a
same-property/same-support contradiction even when the two claims name
different norm models. When an obligation accepts any symplectic or Krein
form, method admission deterministically selects one noncontradicted exact
support and records it in `AdmittedSpectralMethodClassV1`; contradictions on
unrelated forms do not erase a valid route.

A real coefficient field and conjugate-pair symmetry do not imply a real
spectrum. Real ascending/descending ordering requires theorem-closed exact
real-spectrum authority: either a noncontradicted `RealSpectrum` proposition,
standard-linear self-adjointness on an admitted positive-definite metric, or a
Hermitian-definite generalized pencil on such a metric. Indefinite-metric
self-adjointness and approximate propositions do not discharge this gate.

V1 applies a small exact-theorem closure before minting any problem token. For
a standard-linear equation on one exact support, self-adjoint implies normal,
and self-adjoint or normal excludes an admitted nonnormal proposition.
Self-adjointness on an admitted positive-definite metric additionally implies
real spectrum. An exact Hermitian-definite pencil requires that its named
metric descriptor is itself admitted positive-definite; it then implies
induced-operator self-adjointness and normality in that weight, real spectrum,
invertibility of the pencil weight, and regularity of the pencil. An admitted
invertible pencil weight implies a regular pencil, and an admitted exact-grade
invertible polynomial leading coefficient implies a regular polynomial. These
consequences are consumed, not merely conflict-checked: they can discharge
method obligations and make algebraic cardinality available without a
redundant literal regularity claim.
An explicitly contradicted consequence makes the profile inconsistent. The
positive-definite gate is essential: self-adjointness for an indefinite inner
product alone does not imply a real spectrum. Contradicting
`FiniteDimensional` is likewise inconsistent with this schema's explicit
finite metric-space dimensions. This registry is intentionally extensible:
these solid implications are a floor, not a claim that future
theorem-discovery work cannot prove stronger results.

Every complex prefix order has explicit secondary tie semantics. Projective
ordering additionally binds an admitted chart identity and the placement of
infinity. These choices participate in the problem identity; no solver may
silently substitute its own ordering convention.

V1 descriptor partial prefixes require `IncludeProjective`. When projective
infinity may be present, they additionally require
`SpectralOrderingV1::Projective`, whose chart and explicit infinity placement
make the prefix auditable. Exact theorem closure can instead exclude infinity:
an invertible generalized-pencil weight (including the consequence of a
Hermitian-definite pencil) or an exact-grade invertible polynomial leading
coefficient then permits an otherwise admissible finite ordering without
inventing an infinity location. `NoClaim` and `ExcludeWithCount` remain refused
for partial requests because V1 has no finite-only prefix scope with separate
infinity accounting. Set-valued ordering is not a prefix order. Full and
candidate set-valued results remain independent of this prefix-only rule.

Equation-specific propositions are admitted only on their native
representation: Hermitian-definite structure belongs to a generalized pencil,
gyroscopic structure to an exact quadratic polynomial, and palindromic
structure to an exact-grade matrix polynomial. A linearized representation
must carry an explicit future lineage/theorem bridge rather than reusing an
inapplicable proposition by label.

### Method-family obligations

`assess_method_class` checks schema prerequisites only; it neither chooses a
concrete implementation nor claims convergence.

- self-adjoint Lanczos: standard ordinary direct endomorphism, positive metric,
  exact self-adjoint proposition in that metric, finite-dimensional
  regularity;
- generalized self-adjoint Lanczos: ordinary direct generalized pencil,
  positive metric and an exact Hermitian-definite-pencil proposition, whose
  theorem closure supplies regular-pencil evidence;
- general Arnoldi: standard ordinary direct problem and finite-dimensional
  evidence;
- polynomial Krylov: direct exact-grade polynomial regularity; descriptor
  variants additionally require regular-descriptor evidence and an explicit
  infinity policy;
- Hamiltonian/symplectic/Krein families: typed exact form proposition,
  compatible spaces/origin, appropriate even-dimension and regularity gates;
  the Krein/J-orthogonal family additionally requires an admitted
  nondegenerate indefinite domain metric;
- monodromy Arnoldi: typed positive period, consistent multiplier/exponent
  units and branch semantics, and well-posed-monodromy evidence;
- descriptor pencil: generalized/polynomial descriptor semantics, explicit
  infinity policy, regular-descriptor evidence, plus equation-family
  regularity; and
- operator-function Krylov: ordinary standard analytic-function origin,
  non-`NoClaim` branch policy, and analytic regularity evidence.

An `Ordinary` generalized pencil must also carry proposition-bound evidence
that its weight is invertible, unless the exact Hermitian-definite-pencil
theorem already supplies that fact. An `Ordinary` matrix polynomial must carry
grade-bound evidence that its leading coefficient is invertible. Regularity of
a pencil or polynomial alone does not exclude projective roots, so it cannot
justify ordinary finite-spectrum semantics.

Domain/codomain metric dimensions always describe the operator spaces on which
the admitted problem acts. For `Quotiented`, these are already-induced
post-reduction spaces. The certified nullity belongs to the pre-reduction
lineage and is never subtracted from the declared dimension a second time.
Gauge and quotient identities participate in both evidence and problem
identity. A pre-reduction padded/omitted zero count may consequently exceed the
target-space dimension under `Quotiented`; it remains bound to the exact
gauge/reduction context and must equal the certified nullity when gap semantics
are requested.

`assess_gap_semantics` is deliberately separate from algorithm-family
admission. It mints a gap-interpretation token only when gauge/nullspace and
serialized-zero conventions are both proposition-bound and explicit. Any
declared padded/omitted count must equal the certified gauge nullity; a fixed or
quotiented problem may instead certify that its resulting sequence contains no
structural zeros. This prevents a mathematically admissible eigensolver from
silently interpreting a legacy zero-padded or zero-omitted sequence under the
wrong gap convention.

## Eigensolver service API and semantics (bead bfid)

The `service` module is the generic resumable symmetric-eigensolver
lane (ownership rule: generic operator spectra live at L1; domain
crates adapt their operators DOWNWARD to `service::SymmetricOp` — an
fs-solver `LinearOp` adapter is an L3 shim, deliberately not defined
here).

- `EigenService` wraps fs-la Lanczos/LOBPCG behind bounded,
  cancellable, resumable ticks: plain-data state, `clone()` IS the
  checkpoint, split runs replay bitwise-identically to straight runs;
  `tick(op, cx)` gates the operator dimension, polls `cx.checkpoint()`
  between backend steps, caps `steps_per_tick` at
  `MAX_STEPS_PER_TICK`, and ROLLS BACK the in-flight tick on
  cancellation or non-finite operator output so accepted state is
  never corrupted.
- `CertifiedEigenvalue` carries the Ritz value, TRUE operator
  residual, and the Weyl containment interval
  `[value − residual, value + residual]`: an EXISTENCE certificate
  (some eigenvalue lies inside), never distinctness. `gap_report`
  clusters by interval overlap (lower-bound sort, both-endpoint hull
  merge) and reports "multiplicity at the achieved resolution" plus a
  certified leading-gap lower bound that is 0 when hulls touch.
- Warm starts (`EigenService::warm`) seed from previous Ritz vectors
  for parameter continuation; wrong-size, non-finite, zero,
  rank-deficient, or overflow-sized seeds refuse typed.
- Typed refusals throughout: `InvalidQuery`, `DimensionMismatch`,
  `NonFiniteOperator` (with rollback), `InvalidSeed`,
  `SubspaceExhausted` (sticky Lanczos invariant-subspace breakdown —
  never a zero-vector re-entry), `Unconverged` (budget; resumable),
  `Cancelled` (rolled back; resumable).
- The pre-existing dense monitor path (`symmetric_eigenvalues`,
  `spectral_gap`, `GapHealthMonitor`, `propagate`, `route`) is the
  interpretation layer ABOVE this backend and is unchanged.

Service no-claims: symmetry of the operator is the caller's declared
obligation (a nonsymmetric operator voids every claim);
per-pair intervals do not certify eigenvalue distinctness or exact
multiplicity; the LOBPCG lane is PROVISIONAL pending central proof of
the rank-truncation path (the Lanczos lane is the service claim);
warm-start speedup is logged perf evidence, never an absolute claim;
no admitted-authority coupling — service outputs are numerical
results, not admission-API propositions, and gain authority only
through the admission lane's own receipts.

## Result-truth API and semantics

`SpectralTruthDraftV1` is validated only against the complete
`ValidatedSpectralProblemV1`; callers cannot bind truth to a naked problem
digest. Draft clusters expose only neutral lineage/enclosure inspection.
Favorable localization, multiplicity, internal-resolution, and defectivity
inspection exists only on the non-forgeable `ValidatedSpectralClusterV1`
instances minted inside a successful `SpectralTruthV1`.

### Set-valued clusters

An unvalidated `SpectralClusterV1` draft carries:

- stable lineage identity;
- finite real interval, finite complex box, or projective-infinity enclosure;
- candidate, estimated, or enclosed localization authority;
- independent algebraic and geometric multiplicity claims; and
- per-cluster internal state: no claim, unknown, simple, proven degenerate, or
  positively resolved, plus explicit `NoClaimUndefined` when the producer
  declines to assign internal-separation semantics. The latter does not forge
  a theorem that separation is mathematically inapplicable.

Favorable localization, multiplicity, and internal-state draft claims carry
admitted evidence, but construction alone does not establish that the evidence
belongs to the claim. Validation binds multiplicity/internal evidence to both
cluster lineage and exact enclosure, preventing a stable lineage ID from being
replayed onto changed cluster semantics. Only the validated-cluster view can
report favorable authority or infer defectivity; exact validated
algebraic/geometric multiplicities determine that inference, and numerical
value repetition never does.

The cluster lineage ID is also the canonical membership identity: changing the
represented member set requires a new ID. Internal degeneracy/resolution
receipts additionally bind witness-free algebraic and geometric multiplicity
semantics, so changing `Exact(2)` to `Bounds(2,2)`, changing either axis, or
changing membership/enclosure invalidates the retained internal witness.

### Orthogonal truth axes

Result authority, achieved coverage, scope-boundary state, and termination are
not collapsed into one total order. In particular, an estimate, a residual
bound, and a certified enclosure are distinct propositions, not levels in an
invented lattice.

Coverage is measured in algebraic cardinality:

- partial coverage distinguishes incomplete, exactly satisfied, and a repeated
  boundary cluster returned whole;
- region completeness supports a certified empty set;
- full finite coverage accounts explicitly for finite and projective/infinite
  multiplicity; and
- `NoResult` cannot retain clusters, authority, or boundary claims.

Whenever equation/descriptor regularity makes the total algebraic cardinality
known, both requested prefixes and the sum of returned exact or lower-bound
multiplicities are bounded by that total. Adding one unknown-multiplicity
cluster cannot disable this lower-bound sanity check.

Favorable partial-prefix and complete-region truth requires admitted
equation/descriptor regularity establishing a discrete spectrum. Full finite
or projective completeness additionally requires the theorem-closed known
algebraic cardinality to equal the admitted request. A request remains only a
request: candidate or `NoResult` truth can still be returned when completeness
regularity is unavailable. Directly or theorem-admitted invertibility of a
pencil weight or polynomial leading coefficient forces projective-infinity
multiplicity to zero; included or excluded nonzero infinity accounting is
refused. The same theorem boundary applies outside full accounting: favorable
projective localization, multiplicity, internal-resolution, whole-result, or
achieved-coverage claims are refused, while a raw candidate projective point
with no favorable membership claim may remain as explicitly non-authoritative
diagnostic output. Independently witnessed candidate-boundary geometry or
classification remains orthogonal and does not upgrade that point to spectral
membership; pairing it with achieved coverage does, and is therefore refused
by the same gate.

A partial cluster-closure overrun requires a matching repeated exact boundary
cluster and matching `ClusterClosed` evidence. A satisfied prefix requires
positive boundary separation. Closed-region intersection resolution cannot
exclude boundary multiplicity; open-region resolution cannot include boundary
clusters. Full-spectrum boundary truth exists only with validated full
accounting. At most one projective-infinity cluster may exist, and it requires
the descriptor's admitted include policy.

Whole-result evidence binds a canonical, evidence-free
`SpectralResultSetIdV1`. The result set includes cluster IDs, localizations,
multiplicity semantics, and internal states, so whole-result claims cannot be
replayed after any of those change.

`SpectralTruthV1` always describes the mathematical spectrum of the admitted
problem. Structural roots are not inserted or removed merely because a legacy
serialized view declared them explicitly padded or omitted; that convention
is consumed by gap interpretation and later immutable-view migration code.
Full-spectrum algebraic accounting therefore remains mathematical accounting,
not byte-layout accounting.

Region boundary policy is checked independently of the achieved-coverage
claim. Even a candidate-only result cannot retain favorable boundary evidence
whose included/excluded set contradicts `Closed`, `Open`, or
`RefuseIntersection`. Under any admitted real-spectrum requirement, including
real ordering, a certified enclosure must be real-valued or have a complex
imaginary interval containing zero; an off-axis enclosure cannot inherit
real-spectrum authority. The unique projective point at infinity is fixed by
conjugation and is compatible with an extended-real/projective spectrum; it is
refused only when independent weight/leading-coefficient invertibility
excludes infinity.

### Resource envelopes

Untrusted profiles are bounded before sorting or quadratic comparison:

- at most `MAX_STRUCTURE_CLAIMS_V1` structure claims;
- at most `MAX_REGULARITY_CLAIMS_V1` regularity claims;
- at most `MAX_SPECTRAL_CLUSTERS_V1` result clusters; and
- at most `MAX_REGION_BOUNDARY_REFERENCES_V1` region-boundary references.

Public receipt builders enforce the same collection limits before cloning or
sorting. The regularity limit is five, matching the largest compatible V1
product (for example a generalized descriptor monodromy problem carrying
finite-dimensional, regular-pencil, invertible-weight, regular-descriptor, and
well-posed-monodromy evidence). Valid at-limit fixtures reach canonical problem
identity construction; the caps are not merely early-rejection thresholds.

## Units and normalization

Floquet periods use `fs_qty::Time`; continuous branch anchors use
`fs_qty::Angle`. Spectral values cross the numerical boundary as
`QtyAny`. `SpectralScalingContextV1::normalize` requires exact runtime
dimensions and a finite positive scale; `denormalize` restores the declared
dimensions. Non-finite input, overflow, and nonzero-to-zero underflow fail
closed rather than silently erasing a spectral value. Signed zero is
canonicalized in identity-bearing numeric fields.

## Error model

Admission, truth, adapters, and bridge statements return structured reports.
Resource-limit failures occur before expensive canonicalization or graph work.
Admission/truth diagnostics use a total stable order and adjacent
deduplication; bridge diagnostics follow canonical node/hypothesis/falsifier
and implication order. Canonical identity failures are retained as typed
`CanonicalError` values. No validation function panics on untrusted input.

The legacy hysteresis constructor still panics for an inverted band, which is a
programmer configuration error.

## Determinism class

Admission, problem IDs, proposition IDs, result-set IDs, adapter and bridge
identities, issue ordering, method-support selection, cluster ordering, truth
validation, gap health, and route scoring are deterministic functions of their
inputs. Canonical sets make caller permutation irrelevant. Same-ISA
floating-point identity inputs use stable IEEE-754 bits with signed-zero
normalization. Bridge count/convention arithmetic is exact checked integer
arithmetic.

## Cancellation behavior

RB.1a validators are bounded synchronous metadata operations and do not accept
an execution `Cx`. The RB.M1 lattice validator does accept a `Cx` and polls it
before validation and at bounded node, edge, and topological-sort strides.
Numerical eigensolvers poll at their own documented step boundaries. Canonical
identity builders use the non-cancelling bounded encoder only after collection
caps have bounded their payloads. The result-set encoder's one-mebibyte field
envelope admits all 4,096 public cluster slots even when every slot uses the
largest currently legal complex-enclosure, bounded-multiplicity, and
resolved-separation semantics; the cap test exercises that worst-case shape.

## Unsafe boundary

None. Workspace `unsafe_code = "deny"` applies.

## Feature flags

None. The RB.M1 moonshot is represented only as typed statement/conjecture data;
no numerical or theorem-authority path is promoted into default execution.

## Conformance evidence

- `tests/spectral.rs`: existing deterministic Jacobi, gap, hysteresis,
  confidence propagation, conditioning, and route tests.
- `tests/admission.rs`: G0/G3/G5 battery covering authority typestate,
  wrong-preimage/anchor/verifier/policy refusal, proposition relabeling,
  exact-vs-approximate method gates, typed form nonfungibility, real-spectrum
  ordering and enclosure consistency, generalized/symplectic/Krein routing,
  indefinite-metric refusal, descriptor/Floquet/operator-function
  cross-routing, canonical Euclidean metrics, gauge/padding cross-consistency,
  unit normalization underflow/overflow, reachable structure/regularity/truth
  resource caps, theorem closure and tolerance nesting, identity permutation
  sensitivity, gauge/quotient/padding-lineage replay,
  problem/result/localization/multiplicity/internal-witness replay, certified
  empty regions, the full open/closed/refuse boundary-policy matrix,
  repeated-cluster closure overruns, ordinary and descriptor full accounting,
  regularity-gated partial/region/full truth, theorem-constrained infinity
  accounting, projective-prefix placement, algebraic/geometric capacity,
  explicit undefined-separation no-claim, impossible cluster states, and
  deterministic malformed-input reports.
- `tests/bridge.rs`: G0/G3/G5 battery covering canonical statement replay,
  exact orientation/endpoint transforms, required hypothesis/falsifier gates,
  neutral/tangential/essential-spectrum cases, machine-count separation,
  implication cycles and missing scopes, retained refutations,
  schema/convention/reviewer mutation, budgets, cancellation, and no-theorem
  authority.

The eigensolver-service battery is `tests/service.rs` (sv-001..010,
printed measurements): dense-reference accuracy falsification for both
backends and ends; degenerate-spectrum clustering at resolution;
bitwise split-run equality; warm-start continuation speedup (logged);
typed refusals incl. rank-deficient warm blocks and dimension
mismatch; Weyl-interval containment falsification; hostile-size
overflow refusals; identity-operator sticky subspace exhaustion;
bridging-interval cluster regression with non-finite filtering; and
the sparse seam via an fs-sparse CSR `SymmetricOp` against the
analytic tridiagonal-Laplacian spectrum.

## No-claim boundaries

RB.1a does not implement or claim:

- Arnoldi, Lanczos, polynomial, descriptor, monodromy, Hamiltonian,
  symplectic, Krein, or operator-function numerical kernels;
- residual, enclosure, multiplicity, separation, regularity, or completeness
  proof generation;
- scientific correctness of any externally admitted verifier or policy;
- evidence artifact storage/resolution, revocation, or distributed trust;
- concrete implementation routing, cost models, resume state, warm starts,
  cancellation, or convergence;
- sparse production sheaf-Laplacian eigensolving; or
- completeness merely because a requested scope exists.

RB.M1 additionally does not prove a Maslov--Krein--Evans equality, validate an
opaque witness/proof/checker scientifically, derive a physical instability
from a spectral crossing, remove neutral/endpoint/essential-spectrum
obstructions, or replace formal/native checking and independent adjudication.
Its validated token proves only bounded schema closure, explicit hypothesis and
falsifier coverage, acyclic lattice shape, executable convention bookkeeping,
and deterministic identity.

Those belong to RB.1b/RB.1c/RB.1d and later solver, evidence, routing, and
ledger beads. Until then, missing evidence remains an explicit refusal or
no-claim state.
