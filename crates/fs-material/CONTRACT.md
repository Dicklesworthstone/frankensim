# fs-material — CONTRACT

The constitutive-law and authority-separated identifiability kernel: materials
as mathematical objects — calibration domains, CONSISTENT tangent operators,
thermodynamic guardrails, hysteresis, uncertainty — plus source-closed,
multi-case questions whose execution and assessment identities stay distinct.

Ambition tags: elastic/hyperelastic/J2/RC-fiber laws + calibration [S];
multi-case identifiability authority schema [F]; Ogden staged [F, no-claim
below].

## Purpose and layer

Layer **L3** (FLUX support). Runtime deps: `std`, fs-ad (dual-number
energy derivatives), fs-blake3 (typed content identities), fs-evidence
(model cards, Evidence, and V&V artifacts), fs-matdb (immutable material and
constitutive-model cards), fs-qty, fs-math.
Consumers: fs-solid elasticity (tfz.13), fiber-section beams, lattice
homogenization, the P2 milestone.

## Public types and semantics

- `SmallStrainLaw` trait (Voigt 6-space, TENSOR shear components):
  `stress` / `tangent` / `update_state` / `admissibility` / `card`.
  **The tangent contract**: `tangent` is the exact derivative of the
  ALGORITHMIC stress update at the same committed state — FD-gated for
  every law in conformance (merge-gate discipline, same as adjoints).
- `IsotropicElastic` (E, ν → Lamé), `OrthotropicElastic` (engineering
  constants; construction REFUSES thermodynamically inadmissible Poisson
  sets via compliance positive-definiteness minors).
- `J2Plasticity`: radial-return mapping with linear isotropic hardening
  and the Simo–Hughes algorithmic moduli
  `C = κ I⊗I + 2μθ I_dev − 2μθ̄ n̂⊗n̂`.
- `Hyperelastic` (`NeoHookean`, `MooneyRivlin`): stored energies written
  ONCE generic over the fs-ad `Real` scalar; `piola` is the exact dual
  gradient, `tangent` the exact nested-dual Hessian (9×9). `det F ≤ 0`
  refuses structurally.
- `Uniaxial` trait + the RC flagship pair: `MenegottoPintoSteel`
  (R0/a1/a2 curvature degradation, Bauschinger via branch-state
  asymptote intersection) and `ManderConcrete` (confined envelope
  `f = f′cc·x·r/(r−1+xʳ)`, elastic unload/reload lines with residual
  strain, zero tension).
- `calibrate_bilinear`: segmented least squares recovering (E, σ_y, H)
  from monotonic data with standard-error envelopes and RMS residual.
- `evidence_stress`: wraps any law's stress in `Evidence` with the card
  attached and `in_domain` FLAGGING calibration-domain exit.
- `tensor`: Voigt helpers (deviator, contraction with shear doubling,
  von Mises, Rodrigues rotations) used by the objectivity gates.

## Invariants

1. **Tangent consistency (the merge gate)**: every law's tangent matches
   central FD of its own stress to ≤1e−5 relative across elastic branch,
   plastic branch, cyclic states, and 9×9 hyperelastic components.
2. **Frame indifference**: isotropic small-strain σ(QεQᵀ) = Qσ(ε)Qᵀ;
   hyperelastic P(QF) = Q·P(F) — randomized rotation battery.
3. **Return-map consistency**: after every J2 update, the yield function
   at the returned stress satisfies f ≤ tolerance; dissipation increments
   σ:Δεₚ ≥ 0 (associative flow), total cycle dissipation > 0.
4. **Hysteresis fixture behavior**: M-P virgin curve approaches the b·E₀
   asymptote, tangents E₀/b·E₀ at the extremes, reverse branches soften
   below the elastic line (Bauschinger), symmetric cycles dissipate;
   Mander peaks exactly at (ε_cc, f′cc) with slope 0, softens post-peak,
   unloads to the residual strain, reloads rejoining the envelope.
5. **Calibration round-trip**: synthetic bilinear data recovers E within
   1%, H within 5%, σ_y within 2%, truth inside the fitted envelope.
6. **Inadmissible parameters refuse at construction** (ν bounds,
   compliance definiteness, Ec > Esec, b ∈ [0,1), positive yield).

## Error model

`MaterialError`: `Parameters`, `State` (e.g. det F ≤ 0), `Calibration`.
Out-of-calibration-domain USE is not an error — it is flagged through
`Evidence.model.in_domain` so upstream policy decides.

## Determinism class

**D0**: pure f64 arithmetic with fs_math::det transcendentals; no
iteration counts depend on ambient state (radial return is closed-form).

## Cancellation behavior

All updates are closed-form, allocation-free; P7 by boundedness.

## Unsafe boundary

Zero `unsafe`.

## Feature flags

None.

## Conformance tests

`tests/conformance.rs` (JSON verdicts, suite `fs-material/conformance`):
mt-001 FD tangent gate (every law, cyclic states); mt-002 objectivity;
mt-003 return-map consistency + dissipation; mt-004 hysteresis fixtures;
mt-005 calibration round-trip (+ degenerate refusals); mt-006 rank-one
convexity sampling for NH/MR, constructor refusals, card completeness,
Evidence domain flagging.

## No-claim boundaries

- **Ogden is staged, not shipped**: its principal-stretch energy needs
  eigenvalue derivatives through fs-ad duals. The upstream blocker is
  now RESOLVED (bead t88x: asin/acos/atan/atan2 on `Real` + Dual chain
  rules); the staged Ogden law itself remains follow-up work here.
  Until then NH/MR are the hyperelastic set.
- **Mander cyclic rules are the simplified elastic-unload variant**
  (declared on the card): no cyclic stiffness degradation beyond the
  residual-strain rule, no tension stiffening. Full Mander cyclic rules
  are follow-up work for the fiber-beam bead.
- **J2 has no Bauschinger effect** (isotropic hardening only; kinematic/
  mixed hardening is a follow-up law, the card says so).
- **No damage/softening 3D laws yet**; failure envelopes are declarative.
- **Rank-one convexity is SAMPLED**, a necessary-condition spot test —
  not a polyconvexity proof (interval-certified convexity belongs to the
  certifier tier).
- **Calibration v0 is bilinear segmented LSQ**; nonlinear (Mander/M-P
  parameter) fitting and posterior envelopes land with fs-io CSV
  catalogs and the UQ stack.
- **Homogenization hooks**: homogenized laws register as ordinary
  `ModelCard`-carrying laws; the unit-cell pipeline itself is the
  lattice bead's scope.

## Authority-separated multi-case identifiability schema (bead I10.1)

`identifiability/authoritative.rs` is the current public I10.1 contract. It
models a campaign, not a single coupon, because complementary specimens,
protocols, environments, and observation operators can jointly break a gauge
that no individual case resolves. The schema preserves ambitious continuous,
discrete, mixed, and stratified gauge/theorem targets while refusing to report
their truth before an evidence-bearing assessment exists.

The authority chain is deliberately monotone:

1. `IdentifiabilityProblemDocument` is a canonical but **unresolved** physical
   and statistical question. Raw decode returns only this type.
2. `AdmittedIdentifiabilityProblem::resolve_and_admit` requires concrete
   Context-of-Use/material/model/case artifacts and an exact closed
   `SourceResolutionSet` for every other source. It then mints `ProblemId` and a
   separate `SourceAdmissionId`.
3. `IdentifiabilityExecutionPlan` adds coordinates, parameter actions,
   algorithms, analyzer/build/derivative semantics, seeds, budgets,
   numerical controls, arithmetic policy, initialization, stopping, and
   determinism. Every confirmatory proposition is preregistered as an exact
   `ClaimRequest` together with a sourced dimensionless error metric,
   nondimensionalization, and maximum certified error. The plan also requires
   an exact locally verified authority closure for every transitive execution
   source, and replay requires the caller to present that same closure. It mints
   `ExecutionId` without rewriting `ProblemId`.
4. `IdentifiabilityAssessment` adds product-typed propositions and explicit
   claimed-established/claimed-refuted/claimed-inconclusive/not-assessed
   conclusions with receipts. Its proposition registry must equal the exact
   preregistered execution registry; matching only coarse claim axes cannot
   substitute a new proposition. Its caller-held source authority must agree
   with both problem and execution authority on every overlapping key. It mints
   `AssessmentId` without rewriting the problem or execution. These
   deliberately named claims are not promoted to theorem tokens until a
   method-specific verifier exists.

### Current public schema

The authoritative API and all four independently domain-separated identity
stages are schema **v3**. Decoders require the exact stage version; v1/v2 bytes
are not migrated or reinterpreted implicitly.

- `SourceRef` binds a `SourceKey`, semantic `SourceKind`, expected content hash,
  exact fs-blake3 digest domain, and positive source-contract version. A hash is
  not issuer authentication. Opaque resolution must supply bytes that reproduce
  the digest while retaining the digest domain, preventing cross-domain replay.
  `SourceVerification` records the narrower fact actually established:
  `TypedArtifact`, `HashPreimage { byte_len }`, or `Unverified`; it never
  promotes byte equality into provenance or issuer trust. Synchronous generic
  hash-preimage verification refuses inputs larger than the 4 MiB canonical
  transport bound; larger sources require a future typed or streaming
  artifact-store verifier rather than an unbounded allocation disguised as
  admission.
  Typed Context/V&V/material/model artifacts additionally require their actual
  domain and schema version. `AuthorityDisposition` distinguishes byte-content
  verification, an unvalidated external trust-policy receipt, and an explicit
  unresolved state; unresolved sources cannot mint `ProblemId`.
- `StudyParameter` separates decision purpose (`Estimand`, `Nuisance`,
  `Hyperparameter`, `CalibrationControl`) from inferential treatment
  (`Estimated`, `Profiled`, `Marginalized`, source-bound `Conditioned`, or a
  DAG-checked `Derived` value). Physical quantity/domain/prior, concrete owner
  payload, population/case/specimen/field/hierarchical scope, and honest
  connectivity coverage remain coordinate-free. Case applicability is explicit:
  material-lot, field, and hierarchical scopes carry nonempty exact case sets;
  specimen scope names its one case; and only `Global` implies every campaign
  case without enumerating it.
  `ParameterPrior::Uniform` denotes a continuous uniform measure and therefore
  requires finite, positive-width support contained in the physical domain.
  Singleton support, including signed-zero aliases, refuses rather than silently
  manufacturing atomic mass. A future discrete/Dirac prior family may represent
  atoms explicitly without weakening this continuous-measure invariant.
- `JointConstraint` supports dimension-checked affine, simplex, ordered,
  external-manifold, and stochastic-coupling domains. This prevents an
  apparently valid Cartesian product of marginal domains from silently
  admitting physically impossible parameter combinations.
- `AdmissibleDomainWitness` supplies one finite value for every parameter and
  proves that the *same global point* satisfies every locally evaluable
  marginal, conditioned, affine, simplex, and ordered predicate. Derived
  parameters or opaque external/stochastic constraints require exactly one
  `OpaqueDomainMembershipClaim`. Its admitted
  `AdmissibleDomainCertificate` source has an exact kind/domain/version, and
  its binding commits to the witness, all parameter semantics, the complete
  constraint conjunction, and the exact `SourceRef`s for opaque joint
  constraints. Derived-definition source keys remain bound through the
  parameter semantics and the enclosing problem source registry. This closes
  the former gap where locally evaluable constraints could each look feasible
  while their intersection was empty. The opaque certificate remains a
  source-bound claim, not a theorem minted by fs-material.
- `StudyCaseDocument` binds case purpose, exact initial state, specimen,
  protocol, forward model, prospective or retrospective data intent,
  observations, and discrepancy semantics. `CasePhysicsSources` closes every
  embedded frame, geometry, process, preparation, load, environment, time-grid,
  and explicit-initial-state digest to an exact source kind, role-specific
  domain, contract version, and hash. `ObservationKey { case, channel }`
  prevents local channel-name aliasing across cases.
- `StudyObservation` carries exact QoI/unit/`QuantitySpec`, frame, graph port,
  observation/aggregation/sensor sources, instrument, acquisition channel,
  clock, versions, marginal noise, missingness, saturation, and prospective or
  exact raw-row semantics. At source admission, each QoI/unit is closed against
  the concrete Context of Use;
  every retrospective row is closed against the re-derived `DataLineage` and
  must match its manifest QoI, instrument, acquisition channel, and clock.
  Every `ObservationSourceRef` additionally binds the exact dataset-byte hash,
  locator domain and positive locator-contract version, immutable locator
  hash, and extraction-receipt hash. Distinct manifest IDs cannot alias one
  receipt-independent raw locator identity: changing only an extraction
  receipt cannot manufacture a second observation or cross a holdout
  partition. Receipts still move the complete typed source and artifact
  identities. Every row's dataset hash must equal the experiment authenticity
  source-byte hash. The exact instrument must occur
  exactly once in the current experiment roster, its calibration-certificate
  hash must equal the sensor source, and the observation/protocol clock must
  occur in the experiment clock topology. A blind-falsification case also needs
  an exact `BlindReleaseReceipt` accepted by
  `CalibrationSplit::blind_selection`; its authority receipt moves the source
  admission identity rather than the physical question.
- `JointNoiseModel` keeps marginals and dependence separate. Dense correlation
  requires the exact composite channel set and marginals with finite standard
  deviation; bounded/unknown noise is never assigned a fictitious Gaussian
  scale. `Independent` is not a default: it requires an exact source-bound
  `Assumption`. External kernels and honest unknown dependence remain
  representable.
- `ObservationSharingGroup` closes intentional within-case reuse of one raw row
  by naming every repeated row, its exact consumer-channel set, one
  joint-likelihood source, and a justification. The row coverage is exact:
  omitted repeats and declarations for rows that are not actually repeated both
  refuse.
- `DataReusePolicy` represents an explicit pairwise-disjoint raw-data policy or
  explicit sharing groups; there is no inferred reuse default. Admission
  compares exact raw-byte digests, manifest digests, and immutable row-source
  identities, so distinct experiment wrappers cannot hide reuse. Intentional
  reuse requires a non-overlapping `DataSharingGroup`, exact case membership,
  joint-likelihood source, and human justification. Every within- and cross-case
  sharing declaration must name the one global `DenseCorrelation` or
  `ExternalKernel` likelihood; sharing under `Independent` or `Unknown`
  dependence refuses.
- `DistributionFunctional` names location, log-scale, cross-channel
  correlation, missingness-logit, and censoring-logit targets. Derivative units
  are derived from functional and physical-parameter quantities; callers cannot
  supply contradictory units. `InfluenceRepresentation` supports direct,
  state-mediated, composite-DAG, and externally defined mathematical routes;
  computational derivative providers live only in the execution plan. No route
  carries a nonzero/proof receipt. Every influence endpoint must lie in a case
  where its parameter is active, and an influence claim's case scope must
  contain every endpoint and parameter reached through the full transitive
  composite-input DAG, not merely the final functional's endpoint.
- `GaugeDeclaration` represents continuous, discrete, mixed, stratified, and
  explicitly suspected gauges plus quotient, local-section, slice, retained,
  unresolved, or continuous-reduction-with-discrete-residual handling. A
  `GaugeSlice` binds its exact joint constraint, fixed codimension or
  stratum-indexed codimension-profile source, and an exact `GaugeSection`
  transversality/coverage source. Its constraint support must be a nonempty
  subset of the gauge members; internally understood equality/manifold
  codimensions must agree with the declaration, while inequality/order geometry
  requires the external coverage declaration. Mixed gauges can name an exact
  reduction of only their continuous component plus the induced residual
  discrete action without conflating that state with a full quotient.
  `GaugeCompositionDeclaration` is an exact hyperedge over a simultaneously
  active gauge system. It binds the complete member set, an authority-owned
  composition law, `IndependentProduct` or `Generated` semantics, effective
  algebra/orbit geometry, status, and validity cells. Every simultaneously
  active assumed system needs exactly one matching assumed hyperedge; missing,
  partial, duplicate, or scope-incompatible composition refuses instead of
  inheriting iteration order. Independent products must preserve the member
  invariants, while generated actions may declare semidirect, ordered, bracket,
  or richer interaction geometry within the union of coordinates their members
  can move. This schema is deliberately fertile ground for future groupoid,
  Lie-algebroid, sheaf, and cohomological composition theorems; the declared law
  is their substitution-resistant input, not their proof.
  Claim compatibility is fail-closed: a
  raw `Unique` fiber cannot intersect any declared gauge; `DiscreteOrbit`
  requires a purely discrete gauge containing every claimed parameter;
  `MixedOrbit` requires a mixed gauge whose continuous component has an explicit
  reduction and whose residual discrete action is source-bound;
  `QuotientUnique` requires that containment plus a non-suspected gauge and an
  explicit full quotient or full-coverage slice (a partial mixed reduction is
  insufficient); and stratified fibers require a typed
  stratification source. These checks prevent internally contradictory claim
  shapes without asserting that the declared gauge or fiber theorem is true.
- `ParameterExecutionAction` covers every physical parameter exactly once and
  must agree with its declared treatment. Built-in coordinates are checked as
  dimensionally valid full-domain bijections. Numerical rank/singular-value
  controls carry an exact `Nondimensionalization` source. All Five Explicits,
  claim requests, numerical scaling, and rank/conditioning/arithmetic choices
  are execution-identity inputs.
- `TypedIdentifiabilityClaim` is a product of information regime, local/global
  extent, fiber geometry/cardinality, mathematical quantifier, scalar/algebraic
  domain, subject, and campaign/case/stratum scope. Complex extensions,
  mixed/stratified domains, derived functionals, and stratified fibers carry
  exact admitted source definitions. Finite, orbit, quotient, positive-
  dimensional, and stratified fibers remain representable independently of
  local/global extent. These axes are not collapsed into an ordinal evidence
  color.
  Every applicable information/extent/genericity request is conjunctive rather
  than precedence-selected. Cross-axis validation refuses finite-data claims
  over prospective cases, parameters inactive throughout the requested cases,
  transitive influence endpoints outside the requested scope, and gauge-fiber
  mismatch.
  A posterior request additionally carries a `ProbabilityMeasure` joint-prior
  source and requires a declared distribution prior for every free inferential
  parameter active in any requested case. Quantifier sources are narrow:
  `QuantifierRealization` for point claims, `ReferenceMeasure` or
  `ProbabilityMeasure` for almost-everywhere claims, `QuantifierDomain` for
  universal claims, and `ProbabilityMeasure` for probability thresholds.
- `ClaimRequest` binds the complete `TypedIdentifiabilityClaim` proposition and
  a `DimensionlessErrorPolicy` into execution identity. The policy carries
  exact `DimensionlessErrorMetric` and `Nondimensionalization` sources plus a
  finite nonnegative maximum certified error. `requested_axes()` is only a
  derived planner projection; it is not independent authority or identity
  state. A positive/refuting `ClaimAssessment` must reproduce the exact
  preregistered proposition, analyzer, metric, and nondimensionalization, carry
  a separately byte-verified evidence receipt, and report a finite nonnegative
  certified error bound no larger than the preregistered maximum. Inconclusive
  and not-assessed states retain explicit reasons.
- Four owner-local identity declarations register the problem,
  source-admission, execution, and assessment surfaces independently. Each has
  its own version/domain/magic, exhaustive no-`..` field classifier, semantic
  bindings, dependency edge, mutation evidence, and coupling surface.
  `ArtifactHeader.id` remains part of exact execution/assessment transport as a
  ledger label, but is deliberately excluded from their scientific identity
  projections; all other Five Explicits remain identity-bearing.

### Current invariants

1. **No authority from transport:** problem decode re-runs structural admission,
   rejects duplicate/trailing/noncanonical encodings, and still returns only an
   unresolved document. Execution and assessment decode both require their
   exact caller-held, locally verified source sets; serialized verification
   markers cannot mint authority. Problem identity exists only after exact
   source resolution.
2. **Identity non-interference:** coordinates, algorithms, seeds, budgets,
   exact claim requests/error policies, and builds move only `ExecutionId`;
   conclusions and receipts move only `AssessmentId`; trust-policy receipts
   move `SourceAdmissionId` but not the physical `ProblemId`. A ledger
   artifact-label mutation changes exact transport but does not change
   `ExecutionId` or `AssessmentId`; execution source-authority changes do
   change `ExecutionId`.
3. **Multi-case closure:** all parameter scopes, observation endpoints, claim
   scopes, observation- and data-sharing groups, discrepancy roles,
   constraints, influences, and gauges reference the exact canonical campaign
   graph. Non-global ambiguous parameter scopes carry exact case applicability;
   influence endpoints outside parameter applicability refuse. Unreachable
   source registry entries refuse instead of perturbing identity.
4. **Model/source closure:** concrete material/model membership, parameter
   roster/quantity/nominal bounds, initial-state policy, state/protocol version,
   Context-of-Use QoI/unit, experiment QoI membership, instrument/calibration
   binding, clock-topology membership, exact split partition, blind-row source
   binding, blind-release split/hash/commitment authority, and cross-case
   raw-row reuse are re-derived from concrete artifacts rather than trusted from
   problem bytes. Shared concrete keys and shared blind releases must agree
   exactly instead of taking last-writer-wins state.
5. **No hidden defaults:** prospective/retrospective data, missingness,
   discrepancy, dependence, prior absence/not-applicability, connectivity
   absence, gauge handling, and every claim conclusion are explicit variants.
6. **Dimensional closure:** physical parameters retain exact `QuantitySpec`;
   affine joint constraints check coefficient-times-parameter dimensions;
   simplex/ordered constraints require exact common quantities; influence
   derivatives are derived; numerical controls and claim error bounds bind
   exact nondimensionalization semantics.
7. **Global domain closure:** one canonical witness covers every parameter and
   the full joint-constraint conjunction. Built-in predicates are evaluated at
   that point; opaque predicates require the exact full-conjunction membership
   binding, so individually plausible but jointly contradictory domains cannot
   mint a problem identity.
8. **Dependence closure:** independence is explicit and source-bound. Repeated
   row consumption and cross-case provenance overlap require exact sharing
   declarations, and every sharing declaration must agree with the one global
   joint-likelihood source.
9. **Confirmatory-claim closure:** execution owns a bounded nonempty registry of
   exact propositions and sourced dimensionless error policies. Assessment must
   reproduce that registry exactly; its metric and nondimensionalization must
   equal the request, and its certified error cannot exceed the preregistered
   maximum. Posterior claims bind a joint probability-measure source and require
   distribution priors for every free parameter whose applicability overlaps
   the claimed cases. Quantifier roles accept only their narrow source kinds.
10. **Graph closure:** derived parameters and composite influences must be DAGs;
   every free parameter has a declared route or an explicit no-connectivity
   reason; modeled-discrepancy parameter use must remain inside the parameter's
   exact case applicability; influence-claim support closes recursively over
   every composite input; self-correlation and incompatible
   missingness/censoring functionals refuse.
11. **Bounded deterministic transport:** maps/sets and commutative affine terms,
   symmetric correlation endpoints, sharing groups, and dense-correlation
   permutations are canonicalized; counts and text are bounded; floats retain
   exact bits with signed-zero normalization; schemas reject stale/future
   versions. The current public authority schema is v3; v1/v2 transports fail
   closed rather than being reinterpreted under the expanded domain-witness,
   case-source, row-lineage, sharing, scope, and exact-claim semantics. Generic
   synchronous source-preimage verification is bounded by the same 4 MiB
   ceiling.

### Current verification and observability

`tests/identifiability_authority.rs` and
`tests/identifiability_retrospective.rs` emit deterministic JSON diagnostics and
cover canonical problem/execution/assessment round trips, multi-case
ordering and composite endpoints, source/dataset closure, constructive
admissible-domain fixtures, source-authority separation, identity
non-interference, disconnected parameters, derived and composite cycles,
unit-invalid constraints, dense-correlation marginal rules, accidental and
declared within-/cross-case raw reuse, explicit generated and independent gauge
composition hyperedges, typed slice support/codimension/coverage, pure-discrete
versus reduced-mixed residual orbits, modeled-discrepancy case scope,
transitive composite-claim closure, treatment/action coverage, product-typed
claims, real ExperimentArtifact/CalibrationSplit
admission, partition leakage, experiment-QoI/instrument/channel/clock closure,
blind-release absence/replay/commitment/purpose/conflict cases, transitive
problem/execution/assessment authority agreement, and adversarial
transport/version inputs. `fs-evidence/tests/vv.rs` independently exercises
typed row-source injectivity, dataset anchoring, and codec/identity mutations.
Identity tests also exercise stage-domain/magic separation, caller-order
canonicalization, artifact-label exclusion, and exact stage version refusal.

The v3 verification gate additionally requires adversarial mutations for exact
claim substitution, metric/nondimensionalization substitution, certified error
above the preregistered maximum, broad or mismatched quantifier kinds, missing
joint-prior authority, a missing prior on any overlap-active free parameter,
ambiguous-scope case-set drift, and influence endpoints outside parameter or
claim applicability. Omitting any of these mutation classes is a verification
gap, not evidence that the corresponding invariant is vacuous.

The schema/codec path is deterministic D0, has bounded transport and collection
ingress, performs no I/O or parallel work, and contains no `unsafe`. Its final
aggregate structural cap is defense in depth, not yet a complete pre-allocation
work ledger for every nested public input. Downstream symbolic, numerical,
algebraic, differential, sheaf/cohomology, and theorem-discovery analyzers are
not covered by that bounded-work exemption: they must be cancellable
asupersync tile programs and finalize receipts only after drain.

### Current no-claim boundaries

- Structural admission proves closure of the declared question, not
  identifiability, observability, a nonzero sensitivity, laboratory validity, or
  model adequacy.
- A `GaugeSlice` coverage source and a residual discrete-action source are exact,
  authority-checked semantic inputs, not internal proofs of transversality,
  global orbit coverage, induced-action well-definedness, freeness, or quotient
  regularity. They make those ambitious theorem obligations preregistrable and
  substitution-resistant; a future sealed gauge verifier must discharge them.
- The admissible-domain witness is a constructive finite-f64 point, not a proof
  of robust interior, positive measure, connectedness, smooth-manifold
  regularity, or absence of disconnected or singular feasible components.
  Built-in predicates are evaluated locally and exactly under the schema's f64
  semantics. For derived/external/stochastic predicates, v3 verifies the
  *structural binding* between the witness, full conjunction, and one admitted
  certificate source; it does not parse that source or prove the external
  membership assertion true. Future [F/M] certifiers may add rational/interval/
  Taylor witnesses, independently checked decision procedures, constructive
  infeasibility certificates, and sheaf-local witnesses with explicit
  obstruction cocycles.
- Content verification proves byte equality only. `ExternalTrustReceipt`
  retains a receipt issued under an external policy that fs-material does not
  authenticate; it does not prove the issuer, experiment, or scientific
  proposition correct.
- Generic `SourceResolution::verify` is deliberately synchronous and bounded:
  opaque preimages above 4 MiB refuse. v3 does not ship a generic typed or
  streaming artifact-store verifier for arbitrary oversized opaque sources, a
  Merkle-range proof protocol, or a resumable verification receipt. Callers
  must not truncate a large artifact and treat a prefix digest as source
  admission; the existing concrete typed artifact admission paths do not imply
  that generic capability.
- Exact `CasePhysicsSources` close embedded digests to source
  kind/domain/version/hash and admitted bytes; they do not yet parse every
  geometry, process, loading, environment, or time-grid source into a common
  authority-grade semantic schema. Byte closure therefore cannot substitute
  for physical-law or laboratory validation.
- `InfluenceDeclaration` is reachability semantics. Only a separately admitted
  assessment may establish/refute nonzero influence, rank, genericity,
  globality, practical resolution, or a gauge theorem.
- A source-bound zero-discrepancy assumption is still an assumption. Missing or
  uncharacterized discrepancy never means zero.
- `DeclaredSyntheticSelfModel` additionally binds the declared producer and
  the complete forward-model `SourceRef` through a domain-separated production
  preimage. That caller-computable consistency relation is not evidence that
  the experiment actually ran that executable or model. Theorem promotion
  requires a future signed execution receipt binding the experiment hash,
  executable/build identity, producer, and exact forward-model source.
- Declared `ArtifactHeader` time and memory budgets are identity-bearing
  execution policy, but v3 does not yet estimate a document's working set or
  admit it against `memory_bytes`. Fixed transport, per-collection, and
  aggregate structural ceilings remain independent safety limits. A future
  deterministic estimator should charge nested allocation and semantic-work
  units before construction and canonical decoding, then bind that admission
  receipt into execution identity.
- Custom `Debug` implementations on direct authority/capability carriers and
  stage aggregates reduce accidental disclosure and keep aggregate output
  bounded. They are not an authorization boundary: explicit getters expose
  exact audit material, some public nested value and error diagnostics retain
  ordinary `Debug`, and v3 does not claim panic-path-wide redaction. A future
  authority-graph redaction certifier should cover every reachable public leaf,
  error path, and maximum formatted size without removing exact audit access.
- The gauge vocabulary intentionally leaves room for powerful new quotient,
  slice, stratification, groupoid, sheaf, and cohomological theorems. The v3
  fiber/gauge compatibility checks eliminate contradictory declarations; they
  do not establish freeness, properness, orbit completeness, quotient
  regularity, slice coverage, stabilizer strata, or a global fiber theorem.
  Future evidence can promote genuinely new theorems without shrinking the
  ambition of the representable claim space.
- The physical problem identity is coordinate-free but not magically invariant
  under arbitrary changes of units, likelihood, prior measure, model family,
  protocol, discretization, or data reuse. A future proved equivalence may bind
  multiple `ProblemId`s in a new theorem receipt rather than silently changing
  v3 identity semantics.
- Retrospective v3 re-derives the experiment/split lineage directly. It does not
  yet require an `AdmittedVvCase`, so full validation-plan, physical-referent,
  solution-verification, and assumptions-ledger closure is a promotion blocker,
  not a current claim. Direct experiment instrument-roster/calibration-hash and
  clock-membership checks are enforced, but synchronized-clock method/skew
  sufficiency and the complete validation-plan topology still require the full
  V&V-case authority chain. Likewise, `UnitId`/`QuantitySpec` agreement is
  declared and identity-bound but cannot be independently checked until the
  Context-of-Use schema carries dimensional unit semantics.
- Raw rows now bind dataset bytes, a versioned locator, extraction receipt, QoI,
  instrument, acquisition channel, clock, and admitted partition membership.
  Locator injectivity deliberately ignores the extraction-receipt hash, so
  reissuing or relabelling a receipt cannot duplicate one immutable raw row;
  the receipt remains artifact-identity-bearing provenance.
  They still do not carry a typed per-sample value, unit, timestamp, spatial
  coordinate/support, uncertainty contribution, parser transformation trace,
  or authenticated extraction-receipt issuer. That is stronger than a row label
  but remains weaker than a complete authority-grade measurement table.
- Reuse closure prevents an observation row or cross-case dataset from being
  silently counted as independent. It does not prove that the declared joint
  likelihood is statistically adequate, identifiable, correctly normalized,
  or robust to misspecification; `Independent` likewise remains an admitted
  assumption rather than a discovered fact.
- Material-lot, field, and hierarchical parameters now carry exact case sets,
  and applicability/endpoint checks enforce structural overlap. Those checks do
  not prove that a named lot was physically used in a specimen, that a field
  support source denotes the intended region, or that a hierarchy source has
  the claimed population semantics. Likewise, narrow quantifier source kinds
  do not by themselves prove that a realization belongs to the claimed scope,
  that a measure is normalized, or that its support equals the declared domain.
- Posterior claims bind a joint `ProbabilityMeasure` source and require
  individual distribution priors on every overlap-active free parameter. v3
  does not independently prove that the joint measure is normalized, that its
  marginals equal those parameter priors, or that its dependence/support match
  the claimed campaign scope.
- `ClaimedEstablished` and `ClaimedRefuted` are substitution-resistant content
  claims, not theorem receipts. A future sealed verifier must bind the full
  claim digest, problem/source-admission/execution identities, method/build,
  exact consumed partitioned row sources, numerical-error policy, issuer,
  checker, and trust policy before exposing an unprefixed theorem verdict.
- Execution v3 preregisters the complete proposition and sourced dimensionless
  error policy, and assessment enforces exact substitution resistance plus the
  declared certified-error ceiling. Source admission still proves only the
  metric/nondimensionalization/receipt bytes: it does not establish that those
  bytes define a mathematically valid dimensionless metric, that the analyzer
  computed the bound correctly, or that the receipt's issuer is trustworthy.
  A sealed method-specific verifier must independently check those obligations
  before the claimed bound can support a theorem token.

## Historical single-case identifiability prototype (non-authoritative)

The following retained description documents the initial single-case draft for
design archaeology. Its `StudySpecId`/`PhysicalStudyId` wrappers are
crate-private, and its prototype integration test is compile-disabled. Both are
non-authoritative: they are not the current I10.1 public contract and cannot
mint current authority.

The retained `identifiability.rs` prototype owned the admitted *subject* of
later structural, local, generic, global, and practical identifiability
analysis. Its job was to make an inverse problem closed, replayable, and
impossible to silently widen; it did not itself prove an identifiability
theorem. The authoritative multi-case schema above now owns that role.

Ambition tags: closed law/experiment schema and canonical identities [S];
downstream symbolic, numerical, algebraic, and sheaf-theoretic evidence [F/M]
remains external and is carried only through explicit receipts.

### Retained prototype types and semantics

- `MaterialModelBinding::from_cards` binds the complete immutable
  `MaterialCard`, its exact member `ConstitutiveModelCard`, the narrow canonical
  parameter block, law/state/card schema versions, parameter roster, and an
  explicitly supplied constitutive-graph content binding.
- `ParameterSpec` separates canonical physical quantity/domain/prior from the
  optimizer `ParameterCoordinate`; it also records semantic owner, population
  scope, target/nuisance/fixed class, and honest structural-observability state.
- `InitialStateBinding`, `SpecimenBinding`, `FrameBinding`, and
  `ProtocolBinding` bind initialization, geometry/process/preparation,
  handedness/orientation, load/environment/time paths, experiment clock, and
  refinement semantics.
- `DataLineage::from_vv` accepts concrete `ExperimentArtifact` and
  `CalibrationSplit` values rather than arbitrary references. It rebinds their
  exact canonical hashes, raw manifest/source/custody identities,
  preregistration, blind commitment, partition counts, parser/preprocessing,
  split-grouping identity, and V&V schema version. Only preregistered
  calibration row IDs and their immutable source identities are exposed for
  estimation; validation and blind row identities remain capability-restricted.
- `ObservationSpec` binds the measured quantity, frame, model graph node/port,
  observation operator and version, aggregation, sensor/channel model,
  calibration certificate, transfer/filter/support, clock/delay/anti-aliasing,
  marginal noise, missingness, saturation, protocol/refinement, and exact raw
  calibration rows.
- `ObservationPath` states where a parameter can enter the observation
  distribution (mean, variance, covariance, censoring, missingness, or hidden
  state). Its status distinguishes declared connectivity, symbolic nonzero,
  numerical witness, proven zero, and unresolved paths.
- `NoiseDependence` represents the dimensionless correlation `R` in
  `Sigma = D R D^T`; each channel's dimensional marginal scale remains in its
  `NoiseModel`. Channel order is explicit on input and canonicalized on
  admission.
- `GaugeClass` binds group-action, quotient, slice/inverse, stabilizer-strata,
  and evidence artifacts without conflating a gauge declaration with a proof.
- `DiscrepancyModel` keeps `NoModel`, evidence-backed `Zero`, and `Modeled`
  discrepancy distinct. Every observation needs exactly one row.
- `IdentifiabilityEvidence` has five independent fields—structural, local,
  generic, global, and practical. There is deliberately no ordering or
  automatic promotion among them.
- `IdentifiabilityStudySpec` is the opaque, canonicalized closed graph.
  `AdmittedIdentifiabilityStudy` retains two complete `StudyIdentityReceipt`
  preimages:
  - `StudySpecId` binds every field including the header artifact ID and exact
    optimizer coordinates; it is the replay identity.
  - `PhysicalStudyId` omits only the wire-only header ID and an already
    validated built-in bijective coordinate chart. Priors live in canonical
    physical coordinates, so a prior/support/version change still moves this
    identity.

### Identifiability-schema invariants

1. **Exact model membership:** a model card must be an exact content member of
   the material card. Every constitutive-owned parameter must match the model
   roster and full `QuantitySpec` dimensions; undeclared constitutive roles
   refuse.
2. **Version closure:** schema, V&V, material-card, law, state, experiment
   protocol, refinement, observation-operator, sensor, parser, prior, and
   discrepancy versions are identity-bearing. The versions duplicated in the
   Five-Explicits header must agree exactly; no migration/default is inferred.
3. **State closure:** zero state is accepted only for a model card that permits
   it. A model requiring declared state needs a nonzero exact artifact at the
   exact state-schema version.
4. **Coordinate closure:** physical and optimizer domains are finite. Identity,
   affine, and positive-log built-ins must be nonsingular, dimensionally valid,
   and map the full optimizer interval exactly onto the canonical physical
   interval. Estimated domains are nondegenerate. Signed zero is canonicalized.
5. **Observation closure:** every target/nuisance candidate has at least one
   non-proven-zero declared path, or is explicitly retained as unidentifiable
   with a nonzero witness. This is a connectivity invariant only.
6. **Nuisance and data closure:** each nuisance binds the exact admitted split;
   observation operators may consume only the split's calibration rows.
   Validation and blind rows cannot enter pre-release estimation.
7. **Unit-safe likelihood closure:** derivative `QuantitySpec` dimensions equal
   output divided by input; correlation is finite, positive semidefinite,
   exactly unit-diagonal, and covers the exact observation-channel set.
8. **Gauge closure:** gauge members are distinct, known, nonfixed parameters;
   continuous dimension is bounded by member count; action/quotient/slice and
   stabilizer artifacts are all nonzero.
9. **Explicit discrepancy closure:** missing discrepancy is never interpreted
   as zero. Modeled discrepancy requires versioned content, confounding policy,
   and a separate evidence state.
10. **Canonical transport:** exact transport is length-framed, bounded to 4 MiB,
    domain-separated, stale/future-version rejecting, duplicate rejecting, and
    trailing-byte rejecting. Decode re-runs all admission rules and then requires
    byte-for-byte re-encoding, so noncanonical order and signed-zero aliases do
    not mint alternative preimages.
11. **Identity adjudication:** receipts retain the complete canonical bytes,
    schema version, and bounded row count rather than discarding the preimage
    after hashing.

### Error, determinism, cancellation, and unsafe boundaries

`IdentifiabilityError` names malformed text/numerics, zero identities,
cardinality/duplicate/reference failures, exact version/state/nuisance/gauge/
covariance failures, V&V/material refusals, and byte offsets for canonical
transport failures. Admission never returns a partial receipt.

Schema validation and canonical encoding are deterministic D0 programs over
bounded vectors and `BTreeMap`/`BTreeSet` order. They perform no parallel work,
I/O, waiting, or unbounded iteration and therefore satisfy P7 cancellation by
boundedness. The later symbolic/numerical identifiability analyzers are not
covered by this exemption: they must run as cancellable asupersync tile programs
and emit evidence artifacts only after drain/finalize. This module contains
zero `unsafe` and has no feature flags.

### Archived compile-disabled test inventory

`tests/identifiability.rs` is guarded by `#![cfg(any())]` and therefore does not
compile or run in the current conformance gate. It is retained only as design
archaeology; the list below records what that historical prototype suite was
intended to cover. None of these archived cases is evidence for the
authoritative multi-case schema unless an enabled current test above proves the
same obligation:

- G0 exact decode/re-admit fixed points and retained exact/physical preimages;
- G0 canonical invariance under parameter, observation, path, discrepancy, and
  covariance channel-order permutations;
- G3 affine-coordinate quotient invariance while replay identity moves;
- physical prior-version identity sensitivity;
- disconnected estimates versus explicit witnessed unidentifiability;
- nuisance split mismatch and validation/blind-row leakage refusal;
- independent state/protocol/refinement/clock version refusals;
- correlation normalization and exact channel-set closure;
- valid and dangling gauge declarations;
- missing/no-model/zero discrepancy separation;
- independent structural/local/generic/global/practical evidence fields;
- stale/future schema, trailing bytes, truncation, count bombs, and
  noncanonical unit ordering; and
- signed-zero canonicalization plus NaN, reversed-domain, singular-affine, and
  invalid-log refusal.

### Identifiability no-claim boundaries

- **Schema admission is not an identifiability result.** A connected path is not
  a nonzero sensitivity; a numerical witness is not structural rank; local is
  not global; generic is not uniform; practical is not structural. Evidence
  status changes only when a separately admitted analyzer receipt says so.
- **Hashes bind bytes but do not authenticate issuers.** Neither a content hash,
  the current V&V `authenticated` metadata, nor a boolean/current calibration
  flag proves laboratory authority. Capability/issuer trust and signature
  verification remain external admission policy.
- **The constitutive graph binding is caller supplied.** `ConstitutiveGraph`
  does not yet own a canonical authority-grade semantic identity, so this module
  cannot recompute or authenticate that digest. Downstream theorem promotion
  must wait for a graph-owned identity/cross-check or state this limitation.
- **Physical quotienting is intentionally narrow.** It covers only the built-in
  coordinate charts validated here. It does not quotient gauge orbits, custom
  transforms, prior pushforwards, discretizations, frames, units, likelihoods,
  or discrepancy families. Stronger equivalence may absolutely be proved by
  future theorem-producing machinery, but it must mint a new evidence-bearing
  identity/receipt rather than silently broadening this retired prototype
  identity.
- **Gauge artifacts are declarations until proved.** A nonzero action/quotient/
  slice hash does not establish freeness, properness, orbit completeness,
  singular-stratum coverage, or a global slice theorem.
- **The historical correlation surface is cross-channel only.** Temporal
  kernels, parameter-
  dependent covariance, richer censoring/dropout likelihoods, repeated-trial
  hierarchies, derived-parameter graphs, bounded-logit/simplex/manifold charts,
  and, in this historical prototype, blind-release consumption belong to
  follow-on versioned schemas. The authoritative multi-case successor above now
  enforces blind-release admission without retroactively changing this retired
  prototype identity.
- **No physical-law validation is inferred.** Exact cards, protocols, sensor
  models, discrepancy artifacts, and evidence receipts can still encode a
  scientifically wrong claim. Their independent Gauntlet and authority checks
  remain mandatory.

## ConstitutiveGraph and law-node protocol (bead kagp)

Matter is a typed constitutive graph, not a bag of scalars. `graph.rs`
owns the seven-role decomposition as an executable protocol:

- `NodeRole` — the seven roles; TopologyBalance and BulkStorage are
  DECLARABLE but execution-refused (fs-feec/fs-rep-mesh own them); bulk
  transport, reversible coupling, interface, reaction/source, and
  internal memory are executable.
- `NodeDeclaration` / `LawNode` — every node declares ports (name +
  `fs_qty::Dims` + `TimeParity`), state slots + schema version, a
  calibration `ValidityDomain`, a differentiability class, an
  `EnergyBehavior` (including the EXPLICIT `Empirical` no-claim), and
  whether it claims a consistent tangent and/or a free energy ψ.
  `admit_node` refuses incomplete declarations and probes every claim
  (tangent shape, ψ presence for storage-claiming nodes), naming node,
  law, and failed obligation in each typed `GraphError`.
- Thermodynamic gates (test/audit surface): `check_consistent_tangent`
  (analytic vs central-FD, per entry); `check_free_energy_consistency`
  (outputs are the conjugate forces ∂ψ/∂inputs AND the tangent — the
  Hessian of ψ — is symmetric: Maxwell reciprocity);
  `check_psd_symmetric_part` (second law for force→flux blocks via
  Sylvester on the symmetric part); `check_onsager_casimir`
  (`L[i][j] = εᵢεⱼ L[j][i]` from declared port parities: even–even
  symmetric, mixed-parity antisymmetric).
- `LawRegistry` — implementations keyed by the immutable fs-matdb
  `(LawId, LawVersion)`; instantiation validates the card, checks the
  built node's identity and state-schema agreement, and admits it.
  fs-material CONSUMES card metadata, never redefines it.
- `AggregateStateSchema` — the runtime-state codec when laws coexist:
  exact round trip; version (layout-sensitive FNV fold), length, and
  count drift all refuse.
- `ConstitutiveGraph` — admitted nodes composed by typed edges (dims
  must match EXACTLY; one driver per input port), executed as ONE
  deterministic single pass in topological order (insertion-order tie
  breaks). Cycles refuse: implicit coupling belongs to the solver loop
  wrapped around the graph, never a hidden fixed point inside it.
  Execution audits declared-dissipative nodes for non-negative reported
  rates and totals the dissipation.

GRAPH NO-CLAIMS: single-pass execution is not equilibrium; the
dissipation audit checks REPORTED rates (a law that misreports is
caught only by its own gates/fixtures); free-energy and reciprocity
gates run at caller-chosen points and prove nothing globally;
objectivity/frame-indifference remains per-law fixture scope, not a
graph-level proof.
