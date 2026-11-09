<!-- fs-material contract, part 3 of 4; a higher number is newer. Index: ../CONTRACT.md. -->
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

The following description documents the initial single-case draft for design
archaeology. Its implementation and `StudySpecId`/`PhysicalStudyId` wrappers
live only in repository history; its obsolete prototype integration test remains
explicitly compile-disabled. None of these historical materials are the current
I10.1 public contract or can mint current authority.

The removed `identifiability.rs` prototype owned the admitted *subject* of
later structural, local, generic, global, and practical identifiability
analysis. Its job was to make an inverse problem closed, replayable, and
impossible to silently widen; it did not itself prove an identifiability
theorem. The authoritative multi-case schema above now owns that role.

Ambition tags: closed law/experiment schema and canonical identities [S];
downstream symbolic, numerical, algebraic, and sheaf-theoretic evidence [F/M]
remains external and is carried only through explicit receipts.

### Historical prototype types and semantics

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

`ConstitutiveGraph::linearize` uses this same value evaluator and propagates
directional derivatives through every declared consistent tangent, without
assembling a global Jacobian or advancing committed history during a probe.
`execute_cancellable` and `linearize` poll before and after each law callback;
individual law work must remain bounded. Arity, finite input/output/state,
finite tangent actions and dissipation are admitted before returning any
candidate. External directions cannot drive an internally connected input.
The four `graph_linearized` regressions compare a composed history-dependent
chain against analytic and independent central-difference derivatives, and
exercise sparse directions, cancellation and malformed-law refusal.

GRAPH NO-CLAIMS: single-pass execution is not equilibrium; the
dissipation audit checks REPORTED rates (a law that misreports is
caught only by its own gates/fixtures); free-energy and reciprocity
gates run at caller-chosen points and prove nothing globally;
objectivity/frame-indifference remains per-law fixture scope, not a
graph-level proof.
- Moist air transport uses the Wilke viscosity mixture rule and the
  Wassiljewa–Mason–Saxena conductivity rule with viscosity-based mixing
  factors. Ideal molar heat capacities determine mixture Cp/Cv; Prandtl
  number is recomputed from the resulting Cp, viscosity and conductivity.
  These dilute-gas rules, vapor Sutherland extrapolation and Eucken
  conductivity are engineering approximations without a measured error bound.
  The existing x_w = 0.15 ceiling remains; the saturation ENHANCEMENT
  FACTOR (~0.5% at 1 atm) is neglected; RH > 1
  is refused (a condensing state the ideal-gas mixture cannot
  represent — the phase no-claim made structural at this constructor's
  input; the plain try_new window no-claim stands). CO2 fraction is
  not modeled (the USSA dry-air M already carries the standard CO2).

### Material-card J2 plasticity

`plastic::source::resolve_j2_plasticity` resolves Young modulus, Poisson ratio,
yield stress and an explicitly selected hardening datum at one common physical
state. `ResolvedJ2Plasticity::try_from_material` retains an existing complete
multiphysics bundle, including additional properties and original source receipts.
Its `law()` supplies the existing J2 radial return and consistent tangent directly
to native small-strain solvers, including `fs-solid::SmallStrainTetProblem`.
Scalar, curve and sampled-surface coefficients may depend on typed fixed
pressure, composition, moisture or concentration, alongside fixed absolute
temperature. The shared external-condition admission rejects mechanical,
frequency and rate coordinates that would require an evolving physical owner.

The original hardening bases are the plastic-strain modulus H or a
monotonic uniaxial total-strain tangent Et. For the latter, H = Et E/(E-Et) is
derived after resolving both coefficients at the same state; 0 <= Et < E is
required. Perfect plasticity needs an explicit zero source. The original Et
claim and uncertainty remain intact, and no derived measured-H claim is created.
Source, quantity, domain, parameter and representability errors refuse.
Coefficients remain frozen within one load path; temperature re-resolution
starts a different law, not an energy-accounted transfer of plastic history.
The declared strain bound is limited to 0.05 and remains a caller obligation.
Creep, rate effects, softening, phase changes and finite strain are outside this law.

### Source-backed mixed J2 hardening and cyclic reversal

`J2Plasticity::new_mixed` adds an explicitly supplied nonnegative Prager
modulus `H_k` to the original associative return map. The backstress is
`a=(2/3) H_k dev(eps_p)`, reconstructed from the same retained plastic strain;
there is no second history vector or reset on load reversal. The yield function
is `q(stress-a) - sigma_y - H_iso alpha`. Its radial return and algorithmic
tangent use the same effective trial direction and denominator
`3 mu + H_iso + H_k`. `yield_function_with_state` accepts physical stress and
subtracts this backstress; the original `yield_function(stress, alpha)` is the
centered effective-stress diagnostic. Existing `new` and `H_k=0` retain the
original isotropic numerical path and model metadata.

`HardeningPropertyBasis::MixedPlasticStrainModuli` requires separate
`isotropic_hardening_modulus` and `kinematic_hardening_modulus` source claims
in Pa, with the convention above. Each retains its source identity and actual
state query. No monotonic tangent measurement identifies or silently supplies
the split. The original isotropic choices remain isotropic even if a shared
bundle contains a kinematic datum. The exact requirements are returned as a
vector because the mixed reduction has five mechanical roles.

Recoverable free energy includes `H_k dev(eps_p):dev(eps_p)/3`, reported as
`kinematic_hardening_j_m3` separately from isotropic hardening storage.
The corresponding quadratic increment belongs to backward-Euler numerical
loss. Physical dissipation remains `sigma_y delta_alpha`; selected native
linear overstress adds `eta delta_alpha^2/dt` once and includes `H_k` in the
same return denominator. Neither recoverable backstress storage nor numerical
loss becomes plastic heat.

The original `ThermoplasticSource` accepts this same explicit split, but `H_k`
must be constant at every source-line knot across its whole caloric chart.
It then adds no omitted temperature derivative to mechanical entropy. Existing
affine E/H_iso, thermal expansion, yield softening, actual finite caloric
feedback and sourced overstress retain their original owners and accounts.
Temperature-dependent H_k refuses, even if its initial and endpoint values
coincide. Native 3-D small-strain solids consume the resolved law unchanged.

Authored regressions compare an independent reverse-yield limit, distinct
cyclic response at identical monotonic tangent, nonproportional consistent
tangents, separate recoverable/physical/numerical energy, the zero-H_k limit,
source refusal, a native P1 cyclic coupon and finite caloric heating.
They do not establish executed test results or empirical cyclic calibration.
This advances MR22/MR62 without claiming nonlinear cyclic hardening, fatigue
counting, life prediction, finite-strain plasticity or broad bead completion.

### Isothermal liquid momentum and free-surface source bundle

`liquid::surface::resolve_liquid_surface_state` requires sourced density and an
explicit dynamic/kinematic viscosity basis. `LiquidCapillarity::Sourced`
additionally requires positive surface tension at the same complete query;
`Neglected` is an explicit approximation, not a zero-valued material claim.
An existing shared material bundle can be adopted without losing extra receipts.
Kinematic conversion uses the same resolved density. Source/domain/unit/selection
errors and unsupported dynamic coefficients refuse the complete reduction.

The actual consumer is `fs-lbm::material::free_surface`. No heat-capacity or
conductivity tuple is needed for this isothermal momentum model. Source conditions
and surface context remain frozen; the adapter supplies no phase diagram, ambient
gas, wetting or purity default, runtime material-state transfer, or evolving
surface-energy account. The ordinary momentum/capillarity and equilibrium-wetting
admission also accepts explicitly typed held composition, moisture and
concentration curves/surfaces, using the shared fixed-external-state policy.
The mixture remains uniform and fixed; this supplies no species transport or
mixing law. Pressure response still requires an evolving pressure owner, and
frequency/rate/strain support cannot become constant Newtonian flow. Legacy
scalar/temperature records retain their original interface selectors. The public
temperature-only support helper used by thermal, bulk-viscosity and cavitation
owners is unchanged.

`gas_field::steady::SteadyAxialGas::solve` lowers its source packs, axial volumes,
supplies and walls through `to_network` into the shared
`network::SteadyGasNetwork`. That owner also admits branched and inlet-fed
recirculating graphs, with at most 64 cells and 256 interior/boundary connections.
The existing fs-thermal `Advected` observer owns both donor-cell enthalpy and
stationary species transport; the original thermal timeline owns wall heat,
combined CFL and energy updates. Cyclic species balance uses the existing
fs-la prepared LU with bounded work and a numerical conditioning allowance;
unanchored closed loops refuse instead of inventing a composition or leakage.
`GasEnthalpyLaw` derives h(T), T(h,p), and a conservative thermal slope bound from
the original NASA-9 source. Caloric laws may differ by stationary composition
only after exact full species/fraction and common enthalpy-reference checks.
Every inlet carries its own source-mixture enthalpy. Initial EOS cell masses
are pseudo-time relaxation weights, not a physical variable-density history.
Only converged, phase-admitted common-pressure states are published, retaining
their thermochemical and material transport sources for existing flow, thermal
and acoustic consumers. Each final power gate bounds every donor/outlet caloric
closure flux separately so inverse errors cannot conceal an inlet energy error.
The axial adapter additionally checks its reconstructed inlet face, then emits
the same `AxialGasField` used by existing ducts and reeds. Balances, work, flow
and bulk-cell low-Mach admission survive this lowering. Downstream faces carry
first-order mixed-cell states; field interpolation and acoustic sampling still
need independent spatial refinement. The scope excludes pressure-drop and
branch-throat velocity prediction, transient variable-density flow, inferred
wall conductances, diffusion, reaction and condensation. Source, phase, work,
convergence and publication failures return no partial solution. This does not
close the complete M0.1-W3 qualification.

`SteadyAxialGas::solve_condensing` adds an explicit equilibrium-water route for
the same positive axial flow and supplied wall data. Each bounded temperature
proposal uses `WaterCondensationModel::partition_at_pressure_temperature` on a
one-second reference flux packet, followed by the original `Advected` observer's
stationary species and enthalpy balance. Its selective pure-water liquid outlet
carries the aligned IAPWS liquid enthalpy; its mixed gas outlet carries the
remaining source-derived composition and reduced mass rate. The original
thermal timeline supplies wall power. The per-cell scalar root does not march
a second transport model or approximate the original NASA caloric law.

The reusable constant-pressure partition retains every noncondensing species
and finite total water, with saturation capacity
`n_water = n_carrier * p_sat(T) / (p - p_sat(T))`. This is distinct from the
fixed-gas-volume capacity of `WetGasParcel`. It reuses the original gas owner's
mole-amount construction and EOS/transport resolution, plus the same IAPWS
saturation/liquid properties and caller-selected constant NASA/IAPWS caloric
alignment. Its explicit liquid enthalpy and `-dm * L_h` phase contribution do
not close a stream's sensible-enthalpy balance by themselves. The constant-p
partition neither infers a volume nor defines a physical duration or wall film.

Liquid output is a **drained** film-reservoir flux at the local equilibrium
temperature, with explicit kg/s and enthalpy W. It is not accumulating steady
storage, and downstream warm/dry gas cannot re-evaporate liquid already drained
from upstream cells. Side injections retain their own species and enthalpies.
The independent energy audit uses original NASA sensible increments plus the
shared latent gap, cancelling formation datums before arithmetic; source/audit
disagreement, loose caloric reconstruction, work exhaustion and failed total or
per-species flow budgets refuse the complete field. Latent W is diagnostic and
is never applied again beside carried liquid enthalpy.

This route publishes the SAME `AxialGasField` using explicit
`AxialGasSampling::CellwiseConstant`: station zero is the authored inlet and
state i fills (station i-1, station i]. Cached states retain their original
source evaluations and phase admission. No primitive chord between different
saturated states is manufactured, so arbitrary acoustic samples remain within
their admitted cell phases. Existing linear primitive fields and dry/network
solves are unchanged. Temperature trials, wall heat, liquid and gas flows,
source-relative power and species residuals remain inspectable per cell.

The shared 273.16--323.15 K / 20--200 kPa ideal-mixture approximation applies.
The model does not infer condensation kinetics, film motion, retained-film
storage, wall thermal evolution, pressure loss or an acoustic validity bound
from a steady residual. Authored G0/G1/G4 tests cover saturated cooling with
mass/energy closure, dry upwind limits, hot/dry side injection, live cell
refinement, arbitrary field sampling and bounded transactional refusal. These
are manufactured-source checks, not experimental breath/film validation or
completion of all M0.1-W3 requirements.

### Equilibrium substrate wetting on the liquid surface owner

`liquid::surface::resolve_liquid_wetting_state` resolves only
`equilibrium-contact-angle` [radians through liquid] from an ordered
(substrate, liquid) interface. The original interface identity, preparation
frames, source receipt, history and named environment survive. The caller
assigns that environment; this response does not resolve an ambient gas EOS.
The angle lies strictly inside (0, pi), and the same frozen-coefficient rules
as the momentum/surface bundle refuse dynamically supported data.
Missing equilibrium data never falls back to advancing/receding endpoints.
`require_liquid_state` checks the declared liquid phase/process/revision and
complete query/quantity descriptors; SurfaceSpec does not carry a bulk-card
content hash. Liquid/gas sigma remains solely in `ResolvedLiquidSurfaceState`.
The sparse free-surface consumer assigns this response to one planar domain
wall. No hysteresis evolution, moving-contact-line regularization, Marangoni
traction or measured-accuracy claim follows from that assignment.

### Incremental J2 work and physical heat

`plastic::energy` audits the original radial-return stress and history update
under the explicitly energetic interpretation of linear isotropic hardening:
elastic energy is `0.5 e:C:e`, hardening storage is `0.5 H alpha^2`, and
physical dissipation is `sigma_y delta_alpha`. Endpoint stress work also
contains the backward-Euler loss
`0.5 delta_e:C:delta_e + 0.5 H (delta_alpha)^2`; that numerical term is
reported separately and is never a physical heat source.

The mixed extension above adds its separately reported Prager storage and
quadratic increment to these accounts. Callers must retain the same elastic
coefficients and both hardening moduli across the increment.
Changing only the yield offset can release stored energy at fixed total strain.
The owner checks coefficients, tensor strain bounds, finite deviatoric plastic
history and its accumulated-strain bound, endpoint yield consistency, and an
independent absolute work residual. Proposals return immutable native stress,
history and energy accounts. These gates do not establish a measured heat
partition for an arbitrary empirical hardening curve or supply creep, phase
change, thermal expansion, finite-strain work, or parameter-work accounting.

The regression source covers independently derived pure-shear work and storage,
elastic unloading without physical heat, reversal with accumulated history,
yield reduction at fixed strain, and invalid-history refusal.

### Energetic J2 caloric feedback

`plastic::thermal::CaloricPoint` binds a virgin native J2 history to a
source-resolved, strictly increasing solid caloric chart. Yield stress and heat
capacity may vary with absolute temperature; Young's modulus E and native
energetic isotropic hardening H may be globally affine in temperature. Poisson
ratio, any selected Prager modulus, density and reference thermal geometry
remain fixed. All property roles share
one material card and source context. Each selected claim supplies one
continuously supported canonical temperature line over the full chart. Every
crossed curve knot and surface-cell edge participates in affine-slope, fixed
coefficient and positive yield/Cp admission; equal endpoints cannot hide an
interior derivative jump. Mechanically bounded strain, stress or rate
coefficients cannot be frozen while their histories evolve. Renamed temperature
coordinates require explicit absolute-temperature semantics. Variable energetic
H requires the explicit native plastic-strain-modulus basis. Uniaxial tangent
conversion retains its constant-E/Et path; variable E with exactly zero Et is
also admitted, but a generally rational derived H(T) is not treated as affine.

This reduction declares energetic isotropic hardening: `H alpha^2/2` is stored
energy; optional Prager storage is separate and also recoverable.
`sigma_y(T1) delta_alpha` is the physical plastic heat supplied to the
existing `UniformEnthalpyStep` transport owner; selected expansion adds the
signed reversible exchange below. Elastic storage, hardening storage,
boundary heat and backward-Euler numerical loss remain distinct. Every thermal
proposal starts from the same accepted mechanical/caloric history; the actual
returned temperature re-evaluates the native J2 increment before separate heat
and combined-work gates accept the candidate. A rejected or cancelled proposal
changes neither history nor time. Optional sourced isotropic thermal
expansion is supported as described below. The model excludes nonlinear or
piecewise-affine E/H slopes, variable Poisson ratio, phase change and
finite-deformation volume work; source receipts and residual gates do not
establish empirical accuracy.

### Sourced thermoelastic J2 expansion

`SourceQuery::thermal_expansion` explicitly selects either zero eigenstrain or
the canonical instantaneous linear expansion coefficient on the same card.
The selected finite scalar, curve or fixed-context surface supplies its complete
pinned temperature line, including every source knot and crossed surface cell.
Negative coefficients and zero crossings are admitted. The existing shared
piecewise-linear primitive computes `f(T)=integral(alpha dT)` without resampling;
`ThermalJ2` passes `eps-f(T)I` to the original J2 stress, consistent tangent and
history update. The reference temperature stays fixed at specimen birth.
Renamed temperature coordinates retain their explicit absolute-temperature
semantics and the original fixed nonthermal conditions.

With fixed Poisson ratio and globally affine E(T)/H(T), free energy remains
`psi=e:C(T):e/2+H(T) p²/2`. Mechanical entropy is
`s=alpha(T) tr(sigma)-(E'/E) psi_elastic-H' p²/2`; internal energy adds `T s`.
The native J2 energy account still evaluates both histories with the new law.
`coefficient_storage_change_j_m3` exposes the missing frozen-history coefficient
change `Cpsi=psi(T_new,e_old,p_old)-psi(T_old,e_old,p_old)`; the increment's
`stored_change_j_m3()` includes it exactly once. Signed reversible caloric
transfer is `tr(sigma_new) integral_Told^Tnew(alpha dT)-Cpsi-delta(T s)`, using
the direct step integral and actual pinned endpoint coefficients. With fixed
E/H and constant alpha, the original `-T_old alpha delta(tr(sigma))` arithmetic
is retained. Physical plastic heat is still `sigma_y delta_p`; numerical loss
remains separate. This closes the actual two-temperature endpoint work balance
without replacing the native return map, history, tangent or time refinement.

Source Cp and actual caloric segment slopes must stay positive after
thermomechanical coupling. Let A=max|alpha|, A'=max|alpha'| and K_max be the
maximum bulk modulus over the complete source span. The sufficient gate is
`rho c_min > T_max [9K_max A² + 3sqrt(3)K_max strain_limit A' + 6sqrt(3)|K'| strain_limit A]`.
The native energy owner bounds the actual elastic tensor by the tensor-Voigt
norm (diagonal squares plus twice shear squares), giving
`|tr(sigma)| <= 3sqrt(3)K_max strain_limit` and
`|tr(C' e)| <= 3sqrt(3)|K'| strain_limit` without assuming exact zero plastic
trace. These bound the signed terms in
`c_eff=rho c_caloric+T[alpha' tr(sigma)+2alpha tr(C' e)-9K alpha²]`, including
both sides of expansion-curve knots. Affine E/H have no second derivatives;
no unsupported plastic-history bound is needed for H. Fixed E/H retain the
previous capacity arithmetic. Failed admission publishes no material or
thermal candidate. Reference mass, geometry and Poisson ratio remain fixed;
phase change, derivative jumps in E/H and finite volume work remain outside
this model.

Focused synthetic tests cover expansion-knot integration and cooling, signed
and zero expansion, affine-modulus curves and typed source surfaces, independent
free-energy/entropy identities, data-only changes to adiabatic temperature and
native 3D loaded shape/plasticity/heat, continuation with coefficient work, and
atomic effective-capacity refusal. These fixtures are not measured material
validation; native execution of the new modulus regressions remains pending
in the configured build environment.

### Structural source bindings for thermal transport

`ThermoplasticSource::thermal_material_point` exposes the exact original card,
complete fixed context and declared temperature axis. `pin_thermal_transport`
selects a nonnegative static property at the actual accepted caloric state,
applies the same mechanical-coordinate/frequency guard as the six native
thermoplastic roles, and returns its original claim pin. The consuming
transport owner admits the property quantity. No material state, caloric law,
plastic history or transport formula is duplicated. `ThermalMaterialPoint`
admits the same typed fixed external conditions for conductivity source axes.
Its explicitly coupled additional axes keep their original physical owners;
an additional temperature cannot be frozen or introduced as a nonthermal axis.


### Declared mushy-phase volume and mobility

`phase::mushy::MushyPhaseLaw` binds the original equilibrium h chart to explicit
constant constituent densities, a permeability scale and a positive minimum
liquid-volume fraction. `PhaseVolumeFractions` converts mass fractions using
constituent specific volumes and checks the chart's bulk density. The declared
Carman–Kozeny closure supplies inverse permeability from liquid volume fraction;
it is neither a universal melting law nor a measured-data assertion. The source
state, original chart and model identity remain available. Fully liquid material
has zero added drag; states below the support floor refuse rather than inventing
solid locking. Equal constant phase densities are separately admitted for
fixed-volume consumers. Separate constituent velocities, elastic skeleton
response, pressure/composition evolution and moving boundaries remain outside
this closure.

Both phase-volume fractions are retained from their own constituent volumes.
Density contrast and nearly complete melting must not reconstruct a small solid
inventory by subtracting the liquid fraction from one. The resistance uses that
retained solid fraction even when the majority liquid fraction rounds to one;
only an actually absent solid has zero extra drag. Focused G1/G3 regressions
check independent one-kilogram phase inventories and dimensional resistance near
the liquid endpoint, including unequal-density rounding and exact pure phases.

`MushyPhaseLaw::try_from_melting_curve` consumes the original source-backed
sharp-melting chart without freezing its pure-phase densities at the fusion
temperature. Solid and liquid sensible states retain the chart's actual density;
only isothermal coexistence uses the original two endpoint densities to convert
mass fractions to additive phase volumes. The same stationary-solid branch and
Carman–Kozeny response apply. A fixed-volume consumer still requires an exactly
constant density over the entire chart, even when the fusion densities match.
G0/G1 regressions drive heating, fusion and cooling through the existing
isobaric heat owner, checking conserved mass, pressure work, phase-volume
fractions, density-dependent mobility and constant-density compatibility. This
constitutive extension supplies no moving mesh, density transport equation,
phase kinetics or experimental qualification.

### Source surfaces in hereditary thermal materials

`thermovisco::source::ResolvedThermoMaxwell` lowers the canonical typed
temperature line from scalar, curve and sampled-surface claims. Every crossed
source-cell edge, original role key, selected claim and fixed material context
is retained. Reference spring and relaxation coefficients must remain exactly
constant across the complete line; the declared modulus/time/free-strain laws
provide the supported temperature response. Source gaps and discrete-only
intervals refuse rather than supplying interpolated coverage.

`ThermoMaxwell::with_history_integration` and the corresponding selector on
`ResolvedThermoMaxwell` opt into `HistoryIntegration::ExactStrainRamp` before
creating material points. Backward Euler remains the default. Exact strain
ramps use the shared `visco::PronyRamp`, preparing once per arm for endpoint
memory, interval-mean stress, endpoint tangent and independent nonnegative heat.
The exact mode has zero constitutive algorithmic dissipation. Original source
coefficients, birth references, thermal eigenstrain, modulus-change work and
applicability limits remain attached. Coefficients are held at the proposed
endpoint temperature within each interval; changing-temperature trajectories
still require refinement, and no exact thermal path or measured accuracy is
claimed. `CaloricPoint::advance` deposits this physical heat and the same signed
reversible exchange into the existing caloric owner. Its endpoint stress and
affine-stress control modes refuse exact strain ramps because those controls
specify endpoint-load work rather than a complete compatible loading history.
Existing structural endpoint-work conjugacy gates remain binding. G0/G1/G3
tests cover scalar-Prony parity, partition-invariant isothermal decay, thermal
work and endpoint tangents, unchanged state on refusal, source selection and
prescribed-strain caloric feedback.

The thermovisco `Energy::stored_change` uses local quadratic constitutive
increments. Exact ramps reuse `PronyRamp::spring_increment`; backward Euler
uses its complementary relaxation factors for both storage and numerical loss.
The balance removes common held-history modulus work algebraically and combines
total/eigenstrain work before rounding. A weak dashpot release can therefore
remain resolved beside a much larger stationary spring energy or coefficient
work, including when absolute endpoint-energy snapshots round equally. Physical
heat, numerical loss and modulus work stay independently reported. This local
accounting does not add compensated branch-state memory or guarantee the
accumulation of sub-ULP state changes over arbitrarily many intervals. Focused
regressions cover both history integrators and a changing-modulus hold.

Hyperelastic, hereditary and J2 source bindings share a single
fixed-external-coordinate policy. Typed pressure, composition, moisture and
concentration may describe fixed external conditions. Evolving mechanical
coordinates and a second temperature coordinate cannot be frozen into a thermal
law. Frequency-dependent values refuse; thermovisco frequency applicability
alone must retain its static/relaxed endpoint. Caloric feedback admits a globally
affine common modulus multiplier and free thermal strain across the complete
retained knot union. It uses mechanical entropy `s=alpha*sigma-S_prime*Xi`
from the free-energy derivative and accounts for
internal energy separately from spring free energy. Physical dashpot loss and
signed reversible exchange, including the native thermal-eigenstrain work,
enter the original caloric transport owner; numerical loss remains separate.
Expansion requires positive total heat capacity over the complete chart and
admitted strain/history range, including instantaneous branch stiffness and
the mixed expansion/modulus derivative term. The native branch-strain envelope
is enforced on expansion-bearing endpoints. Changing reference density, phase,
or either affine derivative still requires additional terms and refuses.

### Work-mean stress control for exact hereditary heat feedback

`Point::linearize_work` returns the affine derivative of the stress used by the
original work ledger. For exact strain ramps this is the interval-mean stress,
using the original stable Prony coefficients and retained thermal/history
references; backward Euler keeps its endpoint semantics. `Step` and
`FeedbackStep` expose that work-conjugate stress even on a zero-increment hold.
`MechanicalControl::AffineMeanStress` solves an explicit affine mean-load line
at every thermal trial, then independently checks the actual mean-stress
residual and prescribed work. It requires `ExactStrainRamp`; existing endpoint
stress controls retain their refusal of an incompatible exact-ramp work law.

The original finite caloric owner continues to receive only physical dashpot
loss plus signed reversible exchange, while mechanical entropy and internal
energy close the combined balance. This enables consistent midpoint inertia
without introducing a separate constitutive or thermal solver. Coefficients
remain frozen at the proposed endpoint temperature, so changing-temperature
paths still require refinement. Source identity, strain/rate limits,
cancellation, and candidate-only publication remain binding.

`ResolvedHeatCapacityCurve::state_at_temperature` initializes the same original
single-phase enthalpy chart by its monotone piecewise-linear inverse. It retains
source identity and chart interpolation and refuses nonfinite temperatures or
extrapolation. Source-bound `caloric_point` uses this operation rather than an
independent Cp integration or arbitrary enthalpy datum.

### Shared phase-dependent conductivity

`phase::conductivity::LiquidMassFractionConductivity` owns the declared positive
multiplier linear in liquid mass fraction. It preserves the existing FEM law's
nearer-endpoint arithmetic and derivative. This explicit phenomenological model
scales a separately admitted base k(T); it does not infer a volume-mixture rule,
solid/liquid material measurements, or additional source authority. The original
fs-conduction facade and live material-flow thermal coupling consume this same
law. The immutable equilibrium h chart remains the phase-state owner.

### Shared linear magnetic permeability admission

`field_storage::linear_magnetic_permeability` exposes the existing absolute
scalar permeability admission used by `ResolvedInductor`. The linear planar
`fs-em` field owner calls this same function, so circuit storage and spatial
fields share dimensions, source support and linear lossless model restrictions.
Fixed temperature and explicitly typed pressure/composition/moisture/concentration
source snapshots are allowed, including source curves and surfaces. Both source
axes and applicability domains use the public
`ResolvedScalarProperty::is_static_at_fixed_external_conditions` classification,
backed by the same fixed-external policy as mechanical and thermal laws. Field/frequency/strain/stress dependence or bounds still require
a richer owner and refuse here. The existing circuit rebind path owns explicit
parameter work; fixed-storage winding and moving-gap owners retain their stricter
source-evolution guards. Focused field-storage tests cover pressure/composition
surfaces, typed external curves, original claim retention and dynamic-role refusal.
This adds no independent constitutive formula or magnetic heating law.

## Reversible thermal field storage

field_storage::thermal supplies source-bound reversible thermal field response for fixed-geometry linear C(T)=epsilon(T) A/d and L(T)=mu(T) N^2 A/l. The exact original material card, selected claims, typed absolute temperature and fixed external context remain attached to each refreshed native element. The existing continuous source-path owner checks the full requested interval; the shared affine-temperature derivative admission is reused by thermoplastic and field-storage consumers. Globally affine positive coefficient paths (including constant and fixed-context surface slices) supply field Helmholtz energy psi=z^2/(2K), entropy s=psi Kprime/K, internal energy U=psi+T s, and the nonpositive fixed-coordinate heat-capacity contribution -2T psi(Kprime/K)^2. A conservative full-interval bound supports the coupled owner admission of positive total caloric plus field capacity under an enforced coordinate limit. This material law does not advance a circuit or thermal state. Coupled owners must account reversible field-caloric exchange separately from Joule heat and add U, not psi, to total internal energy. Derivative jumps, field/frequency-dependent coefficients, second thermal source axes, hysteresis, phase transitions, moving geometry and changing nonthermal context refuse or require another owner. Synthetic G1 tests check potential derivatives, electric/magnetic signs, native geometry, source retention, surface paths, coefficient swaps and validity/refusal boundaries; they are not experimental material validation.

## Spatial affine magnetic thermal material

`field_storage::thermal::magnetic::ThermalMagneticMaterial` pins a complete original material snapshot, its immutable card and selected claims, an explicit absolute-temperature axis, fixed remaining context, and a positive globally affine permeability path. It shares source-path admission, Helmholtz derivatives and capacity arithmetic with lumped C(T)/L(T). For actual cell volume V and gauge-invariant |B| it returns psi=V B²/(2mu), entropy=psi mu_prime/mu, U=psi+T entropy and fixed-B heat capacity=−2T psi(mu_prime/mu)². The full-range conservative capacity bound requires the caller to enforce its B envelope and add all field contributions to the caloric capacity. This is a linear reversible fixed-geometry rung; kinked/nonaffine paths, hysteresis, saturation, transition physics and experimental validation are not inferred.

## Source-integrated directional hygroscopic strain

`state_point::integrate_orthotropic_hygroscopic_expansion` integrates three
explicit dimensionless `linear_hygroscopic_expansion_coefficient_{1,2,3}`
claims with respect to dry-basis moisture (kg water/kg dry solid). Coefficients
are engineering-strain derivatives relative to the original stress-free length.
The typed `hygroscopic_strain_reference_moisture_content` claim must equal the
reference query throughout the admitted path. Temperature is explicitly absolute
and fixed; all other query coordinates and schemas remain identical. Selected
claims are pinned and the original material-line owner admits continuous support
and exact piecewise-linear integration, including surfaces and drying. Missing
axes/directions, ambiguous sources, datum changes and domain gaps refuse.
This is a reversible additive small-strain constitutive path, not a finite-strain,
hysteretic, sorption-energy or stress-assisted transport law. Source and datum
identities accompany the complete endpoint evaluations. Tests use authored
analytic coefficients, not experimental wood validation.

## Source-integrated nonlinear moisture diffusivity

`state_point::KirchhoffMoistureDiffusivity` prepares a nonnegative diffusivity
law D(M) over a closed dry-mass moisture interval at fixed absolute temperature
and fixed remaining source coordinates. The immutable card and selected claim
remain pinned. The original typed material-line query checks every intervening
segment and surface cell; valid endpoints cannot bridge an unsupported gap.
Negative coefficients anywhere on that path refuse, while exact zero remains
a valid source value.

The chart provides D(M), the signed Kirchhoff integral and its symmetric
interval mean. Piecewise-affine source interpolation is integrated exactly up
to roundoff, without allocating during evaluation or subtracting large
cumulative potentials. Equal endpoint values use the continuous D(M) limit.
Preparation inherits the source owner's 4096-knot bound and stores a balanced
tree of nonnegative interval means in linear work and memory. Point and interval
queries require logarithmic work without allocation. Fully covered intervals
reuse their prepared means; only partial endpoint segments interpolate the
original source. Positive length-weighted combinations avoid cumulative-potential
cancellation and absolute-integral overflow. Scaled arithmetic retains a
representable contribution when its intermediate geometric weight underflows.
The retained identity binds the source, environment and complete prepared interval.
This changes evaluation cost, not the admitted physical law or source authority.

This is a constitutive primitive for homogeneous nonlinear diffusion. It lets
a consuming conservative solver propagate a wetting front through D(0)=0
without interpreting the dry state as a distinct impermeable material. It does
not average heterogeneous materials, infer sorption or hysteresis, evolve state,
or establish a continuum-error or experimental-validation claim. Authored
tests cover analytic dry-front potential, table segments, narrow intervals,
constant impermeability, surface cuts and support/negative-source refusals.

## Source-backed invertible moisture sorption

`state_point::MonotoneMoistureSorption` prepares a selected dry-basis isotherm
M(a) on a positive-width, closed relative-humidity interval in [0,1] at fixed
absolute temperature. In this dilute-vapor convention a is water activity.
The original immutable card, claim, typed source coordinates and complete
continuous path remain pinned. Original curve segments and triangulated surface
cuts supply the interpolation; source holes, negative moisture, flat/decreasing
branches, nonfinite capacities and extrapolation refuse. This single-valued
isotherm does not infer a hysteresis branch or a wet/dry basis conversion.

The chart supplies forward/inverse evaluation, one-sided positive storage
capacity, exact piecewise-affine secant capacity and signed local moisture
change. Interval changes integrate source slopes directly, avoiding subtraction
of large moisture datums. Preparation inherits the source owner's 4096-knot
limit. Point/inverse lookup is logarithmic, interval evaluation visits crossed
segments, and all prepared evaluations allocate nothing. Material, temperature,
claim and path identities let heterogeneous transport consumers preserve their
original sources. Equal activity can imply unequal material-local moisture;
this local equilibrium relation neither advances a field nor invents diffusion,
contact resistance, swelling or sorption heat. Authored tests cover nonlinear
inverse/storage values, heterogeneous partition, narrow changes, typed-domain
refusals, fixed-temperature surface cuts and unsupported gaps.

### Differential moisture sorption energy

`state_point::MoistureSorptionEnergy` admits the original selected differential
isosteric heat `q_st = h_vapor - partial_h_bound` in **J/kg water**, as a function
of dry-basis moisture M at explicitly fixed absolute temperature and source
context. This is total vapor-to-bound enthalpy, not merely the net heat relative
to liquid water. No heat of condensation, water property or material-name
default is added. Signed finite source values remain signed.

The existing `query_line_typed` owner admits every intervening curve segment
and triangulated surface cell with one pinned claim. A source gap, changed card,
wrong moisture basis, untyped temperature or extrapolation refuses. Preparation
inherits its 4096-knot limit. Integrals evaluate local segment contributions
directly with compensated accumulation, preserving narrow changes without
subtracting large cumulative potentials. Evaluation allocates nothing; a point
lookup is logarithmic and integration visits only the crossed segments. The
exact segment count is available for a consumer's work admission.

The signed integral `integral q_st(M) dM` has units J/kg dry skeleton. It is
exact for the admitted piecewise-affine source interpolation up to roundoff,
not a continuum or statistical error certificate. Under an explicit rigid-solid,
negligible bound-phase pV approximation and original ideal vapor, the consumer
obtains `partial_u_bound - u_vapor = R_water*T - q_st`. IUPAC's differential
energy/enthalpy convention is stated at
https://old.iupac.org/reports/2001/colloid_2001/manual_of_s_and_t/node31.html .
This source chart does not supply transport, hysteresis, thermal feedback,
independently authenticate compatible isotherm data, or establish entropy
production. Authored tests cover signed nonlinear integrals, narrow increments,
source pins, surface cuts, source holes and quantity/domain refusal.

## Source-backed permanent-magnet recoil

`field_storage::permanent::ResolvedRecoilMaterial` binds the reversible isotropic law `B = mu_recoil H + Br n` at one complete fixed material state. Required source properties are `recoil_magnetic_permeability` [H/m], `remanent_flux_density` [T], signed `recoil_minimum_parallel_field_strength` and `recoil_maximum_parallel_field_strength` [A/m], and `recoil_maximum_transverse_field_strength` [A/m]. The parallel interval must contain zero; all values retain their original units, selected claims, card identity, query coordinates and receipts. Generic permeability, coercivity and energy-product-derived model coefficients cannot substitute for recoil data. All coefficients and limits must be static at the selected external conditions; field/frequency dependence cannot be frozen into this linear law.

`at_query_point` refreshes the complete source-pinned bundle, including Br, slope and operating limits. It creates a new fixed-state snapshot; it does not supply live thermal parameter work or irreversible demagnetization. Specimen orientation belongs to the field owner. Computed parallel/transverse H and signed margins are checked against every source limit; those margins are model-admission estimates, not certified bounds on the numerical field. Remanence-only catalog entries remain insufficient input.

## Source-backed reversible nonlinear magnetics

`field_storage::nonlinear::ReversibleMagneticMaterial` consumes the explicit
`reversible_magnetic_field_strength` source property [A/m] along a named,
dimensioned flux-density coordinate [T]. The original matdb line query pins
one source across the whole interval from zero to the declared maximum,
including fixed-context slices of two-axis source surfaces. Every breakpoint
retains its canonical answer and receipt. Field-dependent secant permeability,
coercivity and material names are not alternative ingress paths.

Admission requires H(0)=0 and finite strictly positive segment slopes. Other
source axes and applicability bounds must be fixed external thermodynamic
conditions. Discrete-only spans, gaps, source switching, nonmonotone response
and extrapolation refuse. Runtime H evaluates the compiled canonical source
path without rebuilding material queries in each field/solver iteration.
Binary search selects the original segment in O(log knots), with no allocation
or independent permeability fit. The same piecewise-affine path supplies H/B,
dH/dB, exact segment-integrated energy `psi(B)=integral_0^B H(b) db` and
coenergy `integral_0^H B(h) dh`. Both integrals are accumulated directly so
small positive coenergy survives cancellation in `B H-psi`.
`flux_density_at_field_strength` inverts this same curve and its original
supported range; it does not introduce a separately interpolated inverse law.
The tangent at zero uses the first segment; an interior knot uses its right
derivative. Bregman loss visits only the original segments crossed by the two
field magnitudes. `bregman_density_3d_j_m3` extends the same original-slope
radial/angular decomposition to spatial fields, including rotation and crossed
saturation segments, without subtracting large stored energies. The planar
entry point keeps its reduction order. This loss is numerical dissipation,
never hysteretic or Joule heat. Six direct convexity tests cover both dimensions,
independent energy identities, rotation covariance and tiny angular losses.
These are floating-point model evaluations, not uncertainty
propagation or enclosures. Focused material regressions compare compiled
curves and fixed-temperature surface slices to canonical source queries, check
inverse and Legendre derivatives, and preserve domain refusals.

The law is reversible, isotropic and without remanence. A hysteresis loop or
initial magnetization branch does not establish that physical model. Static
temperature/material substitution creates a new law; it does not supply
magnetic heat, parameter work, hysteresis or irreversible demagnetization.
G1/G3 regressions in `fs-em/tests/nonlinear_magnetic.rs` exercise analytical
energy/tangents, material-only solved saturation, mixed-material oblique fields,
the linear limit, source-domain refusals and cancellation.

## Reversible longitudinal piezoelectric material

`piezoelectric::LongitudinalPiezoelectric` resolves density, directional
`piezoelectric_compliance_33_at_constant_field` [1/Pa], signed
`piezoelectric_strain_charge_33` [m/V = C/N], and absolute
`piezoelectric_permittivity_33_at_constant_stress` [F/m] from one original
material state. The common local direction is material axis 3, transverse
stress is free, and the electric field is longitudinal. Scalar/curve/surface
queries retain their complete selected claims and fixed external conditions;
dynamic field, stress, strain or frequency support cannot be frozen.

The strain-charge convention is `strain=s_E stress+d E`,
`D=d stress+epsilon_T E`. Coupled stability requires both `s_E>0` and
`epsilon_axial=epsilon_T-d*d/s_E>0`; positive individual diagonal coefficients
alone do not suffice. The single energy is
`w(strain,D)=strain^2/(2s_E)+(D-(d/s_E)strain)^2/(2epsilon_axial)`.
Its derivatives supply both stress and electric field. The reduced
`epsilon_axial` fixes axial strain with transverse stress still free; it is
not asserted to be the fully clamped three-dimensional permittivity.

This is a reversible, fixed-state, small-signal one-dimensional law. It does
not infer off-axis tensors, arbitrary tensor rotation, saturation, depoling,
hysteresis, leakage, dielectric heating or thermal coefficient work. Source
selection preserves evidence, not measured predictive accuracy. G1/G3 tests in
`tests/piezoelectric.rs` independently check strain-charge equations, energy
derivatives, open/short stiffness, signed/zero coupling, SI-sized synthetic
coefficients, and whole-block/source-domain refusal. Actual coupled motion,
charge and native circuit composition belong to `fs-couple::piezoelectric`.

### Full reciprocal piezoelectric tensors

`piezoelectric::tensor::TensorPiezoelectric` resolves a complete stress-charge
law from one original material card/query: all 36 `C^E` components [Pa], 18
signed `e` components [C/m²], 9 absolute `epsilon^S` components [F/m], and
density. Source names declare Mandel order `[xx,yy,zz,sqrt(2)xy,sqrt(2)yz,
sqrt(2)zx]`; the explicit nonzero source frame is shared by all three blocks.
Zero coefficients and symmetry mates must be supplied. Engineering-Voigt `e`,
strain-charge `d`, stress-free permittivity and scalar longitudinal properties
are not interchangeable with these roles. Every selected source receipt and
fixed external condition is retained; dynamic-coordinate support cannot be
silently frozen.

The complementary `resolve_strain_charge` ingress retains all 64 original
datasheet-form `s^E`, `d`, absolute `epsilon^T` and density roles. Names remain
explicitly Mandel: engineering `d` shear columns divide by sqrt(2), and each
compliance shear index divides by sqrt(2). Native scaled solves derive
`C=(s^E)^-1`, `e=d C` and `epsilon^S=epsilon^T-d C d^T` once. The complete Schur
complement must remain SPD; positive separate `s^E` and `epsilon^T` do not
suffice. Original compliance, coupling, permittivity and source receipts remain
unchanged in the retained bundle. Focused synthetic tests reconstruct both
original strain-charge equations, compare the common energy, and reject an
unstable Schur complement with individually positive source blocks.

Native, diagonally scaled Cholesky admission checks symmetric positive
definiteness of `C^E` and `epsilon^S` once. The complete storage is
`W(S,D)=S.C^E.S/2+(D-eS).(epsilon^S)^-1.(D-eS)/2`; its gradients provide stress
and electric field. Field-controlled evaluation returns the same storage at
`D=eS+epsilon^S E`, not electric enthalpy. Local fixed-D stiffness
`C^D=C^E+e^T(epsilon^S)^-1 e` is prepared once. Responses allocate no heap
storage or matrix inverses. Proper frame transformations rotate stiffness,
coupling and permittivity with their respective tensor ranks; source identity
survives and the mapped law identity includes the transformation and values.

G0/G1/G3 synthetic tests in `tests/piezoelectric_tensor.rs` check independent
stress-charge equations, charge-energy derivatives, open/short and signed
coupling limits, simultaneous specimen/tensor rotation, source completeness,
units, support and unstable-block/reflection refusal. The constitutive work is
bounded by fixed 6/3-dimensional operations. No field solve, boundary-specific
open-circuit condition, strain/field amplitude envelope, nonlinearity,
hysteresis, damping, thermal update or experimental validation is claimed by
this fixed-condition material owner; physical consumers own those operations.

`TensorPiezoelectric::longitudinal(n)` derives the existing longitudinal law
from that same complete tensor and density. For unit n in the current frame,
it assumes `T=t n⊗n` and `E=E_n n`: transverse stress and transverse electric
field vanish, while transverse strain and D may be nonzero. With
`r=mandel(n⊗n)`, prepared native elastic solves give `s=r.C^-1.r`,
`d=n.e.C^-1.r`, and `epsilon_T=n.(epsilon_S+e.C^-1.e^T).n`. The reciprocal
energy therefore retains lateral/shear relaxation and signed orientation.
The original material bundle remains unchanged, with no synthesized scalar
source claims; a separate receipt binds the tensor identity, frame and unit
direction into the longitudinal law identity. Existing rod/circuit consumers
can use this law directly. Synthetic G1/G3 tests reconstruct its response and
energy from full strain-charge source data, and check co-rotation, direction
reversal, source preservation and the scalar thermal adapter's explicit refusal.
This homogeneous constraint is not a claim about insulating lateral faces,
fringing fields, finite-wavelength rod modes or bulk charge-free waves. The
scalar affine thermal owner refuses tensor reductions until a tensor-aware
caloric law is supplied.

### Source-bound intrinsic piezoelectric relaxation

`piezoelectric::tensor::relaxation::RelaxingTensorPiezoelectric` binds
anisotropic native Maxwell branches and directed Debye polarization poles to
the original tensor card, physical query and frame. This is explicit model
selection: base `C^E` is the RELAXED stiffness and base `epsilon^S` is the
INSTANTANEOUS clamped permittivity. A single-frequency datasheet value does not
establish either limiting role. The original real coupling remains reciprocal;
no dissipative coupling coefficient, modal damping ratio or pole fit is inferred.

Mechanical sources name all 21 symmetric Mandel increment entries and a
positive time, using the existing native PSD/Gram branch admission. Signed
off-diagonal entries and rank-deficient loss remain physical. Each dielectric
pole names a positive absolute permittivity increment, positive time and all
three components of an explicitly unit direction in the original source frame.
Directions are never normalized to repair data. Frame changes rotate the
original mechanical factors and dielectric directions together with the base.

The common source owner admits at most 64 distinct auxiliary property names;
same-unit names may be shared explicitly. Retained process counts also have
separate finite limits of 64 mechanical branches and 16 dielectric poles.
Every coefficient keeps its source selection and exact units. Frequency may
restrict applicability but cannot be a frozen value-dependence; the complete
admitted band applies to every evaluation.

The positive strain/charge potential includes original base storage, native
Maxwell memory norms and directed polarization storage. With peak phasors
under `exp(-i omega t)`, mechanical increments use
`(-i omega tau)/(1-i omega tau)`; dielectric increments use
`1/(1-i omega tau)`. Direct increment accumulation preserves weak terms beside
large base tensors. Reciprocal-scaled factors avoid squaring a large rate ratio.
Recovered heat is an independent sum of positive memory/friction norms,
reported separately for mechanical and dielectric processes. Memory energy is
part of total storage, not an extra heat source. Frequency evaluation and
constitutive recovery allocate nothing; source preparation occurs at binding.

The `relaxation::transient` companion retains private native Gram and directed
polarization memories, bound to the exact source/configuration/frame and
actual factor coordinates. Zero initialization leaves the original base field
unchanged. A finite source frequency band cannot authorize arbitrary transient
histories. The scalar count is available before state allocation.

`equilibrium_memory(E)` initializes relaxed Maxwell stress and each original
Debye pole at its held-field polarization. It requires the full unrestricted
transient band, preserves native coordinates/source/frame, and refuses
unrepresentable polarization or storage. The enclosing field owner must first
solve with the relaxed permittivity and retain the resulting physical charge.
The static-only `equilibrium_response(E)` instead requires source support at DC
and may accept a finite upper frequency bound. Both recover energy from the
same represented pole coordinate, so static and retained-memory inventories
agree exactly without changing the midpoint heat or residual. Focused
regressions cover equilibrium holds, discharge heat, rotated identities,
multiple poles and representable extreme fields; these are synthetic laws.

For one coupled midpoint window, the material supplies allocation-free tangent
tensors and old-memory right-hand-side corrections. Only proposed state vectors
allocate. Endpoint memories remain transactional. Positive heat is recovered
from actual represented midpoint friction; it is never adjusted to close work.
Compensated sums and actual endpoint increments preserve weak storage changes
beside large bias energies. The enclosing spatial owner applies its explicit
energy tolerance. Midpoint is second order and passive but not L-stable; fast
unresolved modes may alternate and require time refinement.

`admit_held_temperature` admits finite caloric heating only within the exact
original selected claims. It retains the complete typed nonthermal query and
requires the owned axis to be absolute temperature. For a changed temperature,
every base and auxiliary coefficient must be independent of that axis; a flat
sample of a temperature-dependent source does not prove constancy. Each exact
selected claim retains its convex typed validity box, so checked endpoints
admit the intervening temperature interval without extrapolation. No source
switch, coefficient refresh or memory rebind occurs. Two focused source tests
exercise independent base/auxiliary domains and temperature semantics/dependence.

Five harmonic and four transient synthetic groups in
`tests/piezoelectric_tensor/relaxation.rs` and its `transient` child specify
independent coefficient limits, Maxwell decay, Debye charging, second-order
convergence, branch storage/work/heat, anisotropic null spaces, frame covariance
and source/domain/state refusals. They do not qualify a measured ceramic.
Thermal coefficient work, intrinsic coupling loss, DC leakage, domain switching,
hysteresis and depoling remain separate work.

The model conventions follow the primary Debye and piezoelectric-loss theory:
<https://doc.comsol.com/6.3/doc/com.comsol.help.acdc/acdc_ug_theory.05.36.html>
and <https://doc.comsol.com/6.4/doc/com.comsol.help.sme/sme_ug_modeling.05.109.html>.

## Explicit liquid source admission for evolving shear rate

`ResolvedLiquidSurfaceState::try_for_shear_rate_evolution` retains the original
liquid card, viscosity basis and selected claims while explicitly identifying
the temperature and physical shear-rate axes that its consumer must evolve.
It introduces no viscosity law, guessed rate or surface-energy model. The
existing frozen source constructor continues to refuse rate-dependent claims;
ordinary frozen flow ingress also refuses the evolving marker. Only an owner
that solves those physical coordinates may use this bundle. The moving-liquid
enthalpy/rheology owner in `fs-lbm` provides that integration. Focused admission
tests cover dynamic/kinematic viscosity bases, source retention and refusal of
unowned axes; caloric response with shear-dependent storage remains unsupported.

`try_for_pressure_shear_rate_evolution` additionally marks an explicitly named,
positive absolute mechanical pressure coordinate in Pascals. Pressure/stress
semantics are distinct; the canonical legacy `pressure` coordinate retains its
Pascal convention. Only declared temperature, shear and pressure response axes
are admitted. The returned `rheology_pressure_axis` is a required evolution
owner, not permission for a thermal or frozen consumer to hold pressure fixed.
Source selections and validity intervals remain pinned on every re-query.
This admission supplies no liquid compressibility or caloric pressure-work law;
consumers must independently admit density and surface-tension evolution.
The isothermal `fs-lbm` pressure/shear owner requires both coefficients invariant
in its evolving axes and uses the actual solved absolute liquid pressure for
viscosity. Tests cover source retention, explicit owner markers, unowned pressure
response, pressure/stress distinction and original source-domain refusal.

## Source-backed planar orthotropic electrical response

`conductor::planar` resolves two independent principal conductivity or
resistivity properties and an explicit fixed material-to-geometry frame.

Scalar and planar tensor snapshots may explicitly retain one solver-owned
liquid-mass-fraction coordinate using `with_liquid_mass_fraction_axis` and
`try_from_phase_material`. It must exist as a dimensional scalar in `[0,1]`.
The existing static source policy applies to every other coordinate. Refresh
preserves the declared axis, source basis and tensor frame; it does not evolve
phase or infer a solid/liquid mixture law. Coupled consumers must bind that
axis to their actual caloric liquid mass fraction. Unowned response dependence
remains refused.

Resistivity is inverted in its principal basis after source evaluation; the
existing tensor rotation computes `sigma = Q diag(sigma1,sigma2) Q^T`.
The complete material snapshot, source basis, frame identities and orientation
remain bound to the result. Missing directions never become scalar defaults.
The sealed common response also admits the existing scalar owner and a mixed
scalar/orthotropic region enum. Both feed the same native `fs-em` field solve
and its `fs-couple` temperature-feedback driver. This is a decoupled physical
plane with symmetric static Ohmic response; Hall/AC effects, inferred normal
transport and changing principal directions require additional physical laws.
Program anchor: MR43 (`frankensim-material-reality-ui2if.7.1`), remains open.

## Finite ideal-gas internal-energy capacity for coupled exchange

`ThermochemicalGasModel::maximum_isochoric_temperature_slope` bounds `dT/du`
over an explicit source-supported temperature bracket using the original NASA-9
coefficient interval bounds with bounded deterministic temperature subdivision
and each species' `cv = cp - R/M`. Every accepted subinterval proves positive
`cv`; an unresolved or unstable source range refuses. The bound is a timestep
aid; the existing source model still owns actual `u(T)` and inversion.

`ThermochemicalGasModel::with_caloric_continuity` applies the explicit
thermochemical reference policy to each retained species while preserving
composition, transport models, and source-pack identities. `GasEnthalpyLaw`
and isochoric capacity bounds may cross NASA-9 region boundaries only for
species carrying that admitted policy. Raw models retain the single-region
requirement. Each intersected source region receives its own outward-rounded
Cp/Cv bound, including the upper chart selected when the bracket ends exactly
at an internal seam. The minimum admitted capacity bounds the complete bracket;
reconciliation alone does not establish caloric stability. A fixed 13-entry
stack limits subdivision to 8,191 interval evaluations per intersected region,
with at most 16 regions per species.

For corrected models, transport enthalpy-reference identity includes the
complete retained species catalog and every explicit caloric policy,
independently of the current mixture fractions. Mixtures of the same corrected
catalog remain compatible across composition changes; raw and corrected
references, or different correction policies, remain distinguishable. Source
thermochemistry, energy inversion, thermal entropy changes, and reaction
thermochemistry therefore use the same declared caloric reference.

`IsochoricGasParcel::model` exposes that same immutable composition/source model
for coupled film queries. No constant heat capacity, enthalpy substitution or
second gas caloric law is introduced. Regressions cover source-supported
multi-region nitrogen heating, inversion, entropy and pressure response,
unchanged raw admission, and refusal of nonpositive Cv at a selected seam.

## Biased incremental longitudinal magnetostriction

`magnetostrictive::LongitudinalMagnetostrictive` resolves density, incremental
axial compliance at fixed H, the signed reciprocal piezomagnetic d33, and
incremental absolute permeability at fixed axial stress from one original
material bundle. The coefficient property names explicitly distinguish these
incremental, direction-3, stress/field conventions from an ordinary magnetic
permeability or a global magnetostriction curve. No grade name selects data.

The original query must contain an explicitly typed, finite nonzero H-bias
coordinate in A/m. This owner holds that bias and every other external condition
fixed. The shared owned-coordinate/fixed-external admission retains typed
temperature, pressure and composition and refuses unowned dynamic fields,
stress, frequency and rate dependence. Bias- or temperature-dependent source
curves select a tangent at initialization; they are not silently refreshed
along the trajectory. A scalar requires `ConstantWithinValidity`.

The reversible perturbation equations are `epsilon=s_H*sigma+d*H` and
`B=d*sigma+mu_T*H`. All four variables are increments about the already biased,
stress-free prepared reference. `s_H>0` and `mu_axial=mu_T-d*d/s_H>0` admit the
whole coupled block. With independent epsilon and B, one energy
`w=epsilon^2/(2*s_H)+(B-(d/s_H)*epsilon)^2/(2*mu_axial)` gives both stress and H.
This is incremental storage, not the absolute energy of the DC-biased device;
bias preparation and its source work remain outside the model. The constitutive
convention follows Deng and Dapino, *Magnetic flux biasing of magnetostrictive
sensors*, Smart Materials and Structures 26 (2017) 055027, Eq. (1),
<https://doi.org/10.1088/1361-665X/aa688b>. No measured constants are imported.

This is a local small-signal, longitudinal, reversible tangent with free
transverse stress. It supplies neither a full tensor, a stress-preload law,
magnetization/domain history, hysteresis, saturation, heating nor parameter work.
The dynamic consumer must declare a strain/field increment envelope; avoiding
zero-bias crossing does not prove that a measured law is linear in that envelope.
Synthetic material tests compare an independent symmetric-block inverse and
energy, work-conjugate derivatives, signed/bias/temperature data changes,
coupled instability and source-domain refusal. Motion and native circuit
composition belong to `fs-couple::magnetostrictive`.

## Source-defined linear J2 overstress and physical creep heat

`plastic::rate::LinearOverstressJ2` extends the original radial-return owner
with the explicitly declared associative Perzyna n=1 law
`alpha_dot = max(q - sigma_y - H alpha, 0)/eta`. The native tensor convention,
history, trial stress and consistent algorithmic tangent remain shared. For an
actual physical interval dt, backward Euler gives
`delta_alpha = max(f_trial,0)/(3 mu + H_iso + H_k + eta/dt)`, where `H_k=0`
for the original isotropic law and q is formed from stress minus backstress.
Positive endpoint overstress
is a legitimate committed state; the consistency gate checks
`q - sigma_y - H alpha - eta delta_alpha/dt`. Every new interval starts from
the accepted internal variables. Native viscosity and duration must be finite,
positive and representable together with the return denominator.

The energetic-hardening interpretation remains explicit:
`psi=e:C:e/2 + H_iso alpha^2/2 + H_k dev(eps_p):dev(eps_p)/3`. Physical dissipation is
`sigma_y delta_alpha + eta delta_alpha^2/dt`; the viscosity is not added to
hardening storage. The original backward-Euler numerical loss is separate and
never deposited as heat. Fixed-strain stress relaxation can therefore warm the
finite thermal store while reducing actual elastic storage. The original
rate-independent law and arithmetic remain the default when no viscosity is
selected.

`ThermoplasticSource::with_linear_overstress` selects
`plastic_overstress_viscosity` [Pa s] from the same original mechanical/caloric
material card and complete fixed external context. This role is distinct from
fluid shear viscosity. The original claim is pinned; its continuous source line
must stay finite and positive across the complete caloric interval. Piecewise
temperature dependence is admitted; unowned mechanical, rate, history and
frequency source coordinates refuse. Every thermal trial and the accepted
endpoint refresh this claim together with the other material coefficients.
Viscosity affects kinetics and dissipation, contributes no stored free energy,
and creates no additional mechanical entropy term. Existing affine E/H and
integrated expansion accounts remain intact.

`CaloricPoint` and the spatial thermal owner supply actual physical dt.
`resolve_expanded_enthalpy_for_step` and `advance_with_energy_for_step` expose
the same interval-aware source law; their duration-free counterparts refuse a
source that explicitly selected overstress. Coupling proposals and thermal
substeps cannot accumulate extra plastic time. This is small-strain linear
overstress with energetic linear hardening, not a nonlinear rate exponent,
finite-strain creep, damage, phase transformation or empirical validation.

Synthetic `tests/plastic_rate.rs` compares independent radial-return/storage/
heat expressions, a fixed-strain analytic stress-relaxation solution with
first-order refinement, all discrete tangent entries and eta/dt scaling.
The spatial consumer tests actual loaded-solid creep, finite-store heating and
endpoint temperature feedback in `fs-couple`. This advances MR22
(`frankensim-material-reality-ui2if.3.4`); its broader material/history scope
remains open.

## Constant dielectric storage admission for native charge relaxation

`field_storage::constant_electric_permittivity` admits the existing absolute
`electric_permittivity` requirement only when its selected source declares a
constant law over its admitted validity domain. A coincidentally unchanged
sample from a temperature or phase curve does not establish this constitutive
contract. Units, positivity, complete source receipt and original query
semantics remain owned by the common material resolver.

The native leaky-dielectric trajectory uses this admission to change
conductivity at the actual thermal state without silently introducing
permittivity parameter work or reversible thermal entropy. It retains the
complete common source bundle; changing material identity, source claims,
coordinate meanings, tensor frame or phase-axis binding on an existing
trajectory is refused. Variable dielectric storage needs its own work/entropy
law and is outside this entry point.

## Composition replacement at a shared heat and chemical energy endpoint

`ClosedGasParcel::replace_model_with_heat` extends the original conservative
composition boundary with one signed heat increment. At unchanged mass and
volume it reconstructs the new source-backed thermochemical model directly at
`U_end = U_start + Q` through the existing density/internal-energy inversion.
It does not require an intermediate adiabatic composition state or an
intermediate heat-only state to lie inside the material support. The original
conserved-energy replacement delegates with zero heat and keeps its existing
meaning.

`GasCompositionStep` now reports supplied heat, actual stored-energy change and
their residual. Formation energies remain part of the original species internal
energies; callers must not add their release again as an independent heat term.
All endpoint EOS, caloric, transport and acoustic properties come from the same
accepted composition, temperature and pressure. Nonfinite or unrepresentable
energy increments, source/domain failures, failed inversion and cancellation
publish no parcel. This primitive supplies a material state transition, not
reaction kinetics, species transport or open-system mass/chemical work.

## Composition, volume and mechanical work at one caloric endpoint

`ClosedGasParcel::replace_model_with_heat_work_volume` reconstructs the new
source-backed material model at unchanged mass, specified final volume and
`U_end = U_start + Q - W`. Positive W is work delivered by the gas. The original
density/internal-energy inversion resolves temperature, pressure, transport and
acoustics together, without an intermediate composition, heat-only or volume-only
state. The existing fixed-volume heat and conserved-energy methods delegate with
zero work and retain their meaning.

`GasCompositionStep` records both volumes, supplied heat, delivered work,
actual storage change and `delta U + W - Q`. The mechanical consumer must
independently close the supplied work against its actual pressure/volume
trajectory; this material primitive supplies no work quadrature or entropy
production law. Formation references remain inside species internal energies.
Nonfinite/lost exchanges, unsupported endpoints, failed energy inversion and
cancellation preserve the original parcel.

G1/G4 regression cases in the original composition module independently check
an ideal-gas isothermal expansion with equal supplied heat and delivered work,
and failed/cancelled moving-volume reconstruction. Parser checks completed;
native Rust execution was unavailable because DSR/RCH are absent and the shared
filesystem is full.


### Reference-aware hyperelastic material interface

`hyper::cards::HyperelasticMaterial` exposes `energy_density`, `first_piola`
and `consistent_tangent` with defaults delegating to the original native law.
A reference-deformation adapter must transform all three consistently and
retain the underlying elastic-domain owner. Spatial consumers can preserve the
complete material binding instead of stripping it to bare coefficients.
Both `HyperelasticMaterial` and `FiniteStrainPotential` also expose
`energy_density_increment(before, after)` for the same fixed-source,
fixed-history potential. The default uses the adapter's own endpoint energies;
it never bypasses an overridden potential through bare material coefficients.
Native and source-card hyperelastic owners use `Hyperelastic::energy_increment`;
the source-bound Cauchy tensile-limit wrapper retains that same elastic owner.
A held-reference adapter transforms both endpoints with its retained `G` and
scales the change by `det(G)`, preserving original-reference work units.

`ReferencedHyperelastic::energy_density_increment_to_reference` and
`DiscreteGradientPoint::energy_increment_to` also evaluate explicit changes of
the natural reference. Original storage is `det(G) W(F G^-1)`, including signed
work of a source-owned constant density when relaxed volume changes. Direct
determinant increments retain small natural-volume changes; principal-reference
construction retains its supplied length products, so rotating their axes alone
does no volume work. Native points require matching full/remainder energy parts.
Native `energy_increment_to` and `isochoric_energy_increment_to` admit each
endpoint only under its own original law. With a shared nonlinear shape they
carry finite coefficient changes alongside deformation: Neo-Hookean mu/lambda,
Mooney-Rivlin c10/c01, Gent mu, every Ogden mu, HGO c10/k1, and the original bulk
modulus. Gent's locking invariant, Ogden exponents/order, and HGO's fibre axes,
activation convention, dispersion and k2 define that shared shape. The existing
small-strain Ogden series shares its original invariants with the paired
amplitudes; no repeated eigensolve per coefficient or synthetic intermediate
material is introduced. When both endpoints admit the original bounded
normalized-strain expression, it also retains the absolute isochoric invariants
needed for stiffness work beneath a finite common dilation.

Changed nonlinear shapes or families retain their own endpoint remainder
potentials, with bulk evaluated separately so an unrelated large bulk preload
cannot erase the remainder change. An absent fibre family remains inactive
before exponentiation. If symmetric cross-products overflow while both own
endpoint potentials remain finite, the native increment uses their compensated
endpoint difference; it still refuses nonfinite original storage. Source and
reference events supply no history transport, heat closure or dissipative source
law. Physical source-change regressions cover all five native models and a
fixed-deformation solid stiffness event beneath an unchanged `1e24 Pa` bulk
modulus, including a held mechanical interval and reversal.

`DiscreteBulkPoint::energy_increment_to` retains both snapshots' original bulk
family, modulus, natural volume and physical mixed volume. Quadratic,
log-squared and Simo--Taylor laws propagate reference, modulus and retained-volume
changes through the same original finite-pair arithmetic. Distinct families
retain their separate original potentials. Opposing source scales fall back to
compensated endpoint products when symmetric cross-products overflow. The
operation neither resets retained volume nor alters either admitted convex
branch. Focused original-owner tests retain second-order reference/modulus work
under a `2^80 Pa` bulk preload, small work at fixed elastic strain, finite family
changes and finite storage across opposing `1e200`/`1e-200` source scales.

All three discrete bulk families test the retained total volume before scalar
inversion. That candidate is returned only when the original pressure and work
allowances admit it; perturbed pressure must still lie in the monotone range
of the original declared volume bracket. Thus interpolation roundoff in a held
pressure does not create volume motion or stored energy. Resolved loads still
use the native inverse, and the quadratic branch forms its total-volume
increment in double-double precision before rounding the endpoint. No new
pressure threshold or production tolerance is introduced. Tests distinguish
an admitted one-ulp pressure perturbation from a resolved load, strict
zero-pressure-tolerance refusal, and an out-of-branch pressure despite a loose
allowance. The original mixed-solid `1e12 Pa` bulk-preload fixture now retains
volumes bitwise with zero stored change through an unchanged-geometry interval;
its nonlinear bulk time-convergence ladder remains covered separately.

Native finite-pair arithmetic carries small changes separately from absolute
storage, using symmetric product/quotient identities, log1p/expm1 identities and
the shared double-double invariant arithmetic. It retains the original
Neo-Hookean, Mooney-Rivlin, Gent, Ogden and HGO laws, all three declared bulk
energies, and each endpoint's actual fibre activation. The original bounded
small-strain expressions are used only when both endpoints satisfy their
bounds. General Ogden powers use the exact two-eigensystem divided-difference
trace identity, including coalescing roots; the second-order AD spectral lift
is not misused as a finite increment. Both finite and derivative spectral
secants scale from the larger endpoint power, avoiding an overflowing power
ratio when the original endpoint energy, stress and tangent remain finite.
This is not stress quadrature or a
replacement potential. Both endpoint deformation domains are checked even for
coincident states; invalid or unrepresentable native storage returns NaN for
the consuming energy gate to refuse. `isochoric_energy_increment` retains the
native remainder without evaluating an unused bulk penalty. These are
floating-point identities, not interval certificates or new stability claims.

Focused tests retain small shear work beneath a `1e24 Pa` unchanged bulk
baseline for all five models, check noninteger/negative Ogden powers at finite
repeated stretches, cross fibre/chain/bulk expression boundaries, and preserve
adapter and held-reference semantics. Changing material parameters, history,
or prescribed `G` across the pair still requires its own external-work account.
The public `state_point::is_fixed_external_coordinate` classifies one coordinate;
it does not replace typed-query, source-axis, validity or path admission.


### Explicit planar dielectric source laws

`field_storage::planar` admits explicit scalar or fixed-principal-axis
in-plane permittivity. Orthotropic sources provide independently selected
`electric_permittivity_1` and `electric_permittivity_2` in absolute F/m,
with the existing named `PlanarElectricalFrame`. The shared storage source
gate retains units, selected claims and applicability; the original tensor
rotation supplies `epsilon = Q diag(epsilon1,epsilon2) Q^T`. No conductivity
orientation or omitted normal response is inferred. The physical plane must
be decoupled from the omitted normal direction.

Snapshot resolution and trajectory admission remain distinct. Retained
charge owners require every principal coefficient to be explicitly constant
throughout source validity. Equal samples do not prove constancy. Changing
dielectric coefficients or axes during a trajectory still requires an
independently accounted parameter-work/thermal-entropy model; this increment
does not supply one. All source coefficients remain in the original complete
material bundle.


### Explicit three-dimensional dielectric source laws

`field_storage::spatial::SpatialDielectricLaw` supplies the same absolute-SI
storage admission to full three-dimensional fields. Its orthotropic law needs
all three `electric_permittivity_1`, `_2`, `_3` sources and explicitly admitted
orthonormal principal directions in the mesh frame. The existing tensor rotation
forms the complete symmetric epsilon; no normal coefficient or conductor frame
is inferred. Source applicability and representable tensor positivity are
checked before returning D=epsilon E. Weighted FEEC separately checks its
assembled tensor's interval positivity.

This is a static, lossless linear snapshot. Dispersion, evolving orientation,
nonlinearity and temperature/geometry parameter work retain separate owners.
`tests/field_storage/spatial.rs` checks the uniform-capacitor limit, independent
principal-frame vector/energy projection, axis-sign invariance, and refusal of
missing, wrong-unit, nonpositive or unsupported dynamic sources and bad axes.


### Complete liquid bundles through temperature and shear evolution

`ResolvedLiquidState::try_from_material` adopts an already resolved material
bundle after checking the original four transport requirements and explicit
dynamic/kinematic viscosity basis. Additional coupled properties remain
attached. `LiquidTemperatureLaw::try_from_material` also checks the original
immutable card and temperature coordinate. Its existing exact-source refresh
therefore retains all selected claims, query meanings and receipts through
later temperatures, including an active surface-tension claim.

`ResolvedLiquidSurfaceState::try_for_shear_rate_evolution_with_capillarity`
keeps the scenario's explicit capillary choice in that same complete bundle
while preserving the original temperature/shear ownership. The old entry
continues to declare neglected capillarity. Admission of an evolving bundle
does not authorize freezing an unsupported coefficient; the receiving
runtime must own or refuse each active response. The thermal free-surface
consumer currently requires source-declared validity-wide constant sigma.
No surface-energy, Marangoni or wetting dynamics are inferred by adoption.


### Owned liquid endpoints from original caloric heating

`LiquidTemperatureLaw::prepare_caloric_endpoint` uses the original equilibrium
enthalpy chart and accepted state to prepare an owned, noncapillary
`LiquidCaloricEndpoint`. Both endpoints must remain liquid. Card/curve identity,
selected caloric claims and complete nonthermal query context remain exact.
Authored charts without receipts retain their authored assumptions.

The endpoint owns its resolved surface source and reanchored temperature law.
The ordinary reference-state path reuses initial data and resolves final
transport once; scalar heat/work trials can use the original chart without
requerying every property. `reanchored_at_temperature` preserves the complete
claim bundle and separately retained density-derived expansion metadata.
Kinematic viscosity is converted with the endpoint density.

Source/chart density disagreement is measured at both endpoints and must fit
the caller's explicit accuracy budget; neither density is rewritten.
Requested and represented specific enthalpy increments remain separate.
Fixed-pressure Cp already includes ambient p*dV. Mechanical correction work
belongs to the physical caller; projection/remap losses never become heat.
The helper changes no native state, volume or clock. The existing thermal
owner must accept those together with its mass, energy and cancellation gates.

### Relative entropy on the admitted enthalpy/phase chart

`EquilibriumEnthalpyPhaseCurve::specific_entropy_change_j_kg_k(h0,h1)`
integrates `dh/T(h)` on the SAME bounded piecewise-linear caloric chart.
Sensible intervals use their logarithmic integral, evaluated with log1p near
isothermal limits; latent plateaus contribute `delta_h/T_transition`.
Fixed-order compensated accumulation, original endpoint/source-domain
admission and explicit representability checks apply. Reversing endpoints
reverses the sign. No apparent plateau heat capacity, additional latent energy,
absolute entropy reference, chemical mixing or variable-pressure work is
invented. Its thermodynamic scope is the chart's fixed pressure/composition
convention. Coupled owners retain responsibility for boundary entropy and any
additional state/work coordinates.

Two focused regression definitions cover independent sensible-plus-latent
integrals, subinterval additivity, reversal/domain refusal and a nearly
isothermal sensible limit. Source syntax parsing and independent scalar
arithmetic passed; native Rust tests were unavailable in the implementation
environment. This is not a continuous transport or entropy-production proof.


### Source-backed planar dielectric thermodynamics

`ThermalPlanarDielectricMaterial` pins one original complete material bundle,
selected claims, nonthermal context and explicit scalar or fixed orthotropic
dielectric axes. Positive globally affine principal epsilon(T) must be supported
throughout the declared absolute-temperature interval, including interior source
knots. It reuses the original affine-path, Helmholtz derivative and capacity
helpers, and re-evaluates the full source bundle at each accepted temperature.

At fixed displacement D and volume, the law supplies
`psi=V D.epsilon(T)^-1.D/2`,
`S=-partial_T psi=V E.epsilon_T.E/2`,
`U=psi+TS` and the nonpositive fixed-D field heat capacity. The field contribution
is relative to the original zero-field caloric chart, not another heat store.
A conservative capacity bound applies over the explicit T and |D| envelope.
A thermal consumer must add ALL field contributions assigned to each caloric
receiver and keep their combined heat capacity positive.

Source continuity compares original card, exact claims, full query semantics
and actual values; a requery's stronger pinned selection receipt is not confused
with a change in physical material. Axes and composition/pressure context stay
fixed. No dispersion, nonlinear polarization, phase-dependent entropy or
electrostriction is inferred. The focused source-derivative and coupled tests
are in the fs-couple dielectric thermal test child. Final source syntax parsing,
bounded source review and independent scalar arithmetic passed; native Rust
tests were unavailable, so no executed native or empirical qualification is
asserted.


### Source-bound thermal longitudinal piezoelectric potential

`ThermalLongitudinalPiezoelectric` pins the complete original longitudinal
material bundle, card, source claims, axis quantities and fixed nonthermal
coordinates. The owned absolute-temperature interval admits affine s_E(T),
signed d(T), epsilon_T(T), constant reference density, and positive definiteness
of the full coupled block A=[[s_E,d],[d,epsilon_T]] throughout. In natural
coordinates x=[strain,D], y=A^-1*x, the same native response supplies
Psi=V*x^T*y/2, S=V*y^T*A'*y/2, U=Psi+T*S, and
C_x=-T*V*(A'*y)^T*A^-1*(A'*y).

Signed/zero coupling is valid; signed-zero reference checks use the admitted
source-path scale without replacing actual source values. A normalized
full-block eigenvalue/derivative bound supplies a conservative minimum field
heat capacity over the caller's strain/D envelope. This is an admitted model
bound, not a numerical certificate. Common affine-source and capacity routines
are shared with the existing field-storage thermodynamics; positive-coefficient
source admission retains its original behavior. No separate constitutive engine,
pyroelectric offset, expansion, depoling, hysteresis, phase transition or
dissipative ceramic law is inferred.

Focused synthetic regressions in fs-couple compare the complete analytic
potentials, thermal derivatives, signed/zero coupling and actual coupled
circuit/motion/caloric evolution. Source syntax and independent scalar references
were checked; native compilation/tests were not run in the editing environment.


### Final nozzle admission after enclosing heat or mechanical work

`GasTransportStep::admit_endpoint_nozzle_pressures` reapplies the existing
active-nozzle pressure-order check to a complete endpoint parcel map after an
enclosing physical owner supplies heat or mechanical work. Original parcel IDs
and actual positive transfers remain authoritative; the adapter samples no new
flow, spends no inventory and changes no accepted parcel. It checks pressure
ordering only and does not extend the original transport-stage entropy account
to subsequent heat/work. The multi-chamber DG integration exercises a late
heated-reservoir reversal and complete rollback. Source/API review and syntax
parsing completed; native execution remains unverified in this environment.


### Source-backed reversible isotropic thermoelastic storage

`elastic::thermal::AffineIsotropicThermoelastic` retains the original material
card, selected elastic and expansion claims, full query semantics and declared
positive temperature interval. It admits independently affine Young's modulus,
Poisson ratio and instantaneous linear expansion, with E positive and nu in
(-1, 1/2) throughout. The original pinned source integrator supplies free strain.
Bulk and shear stiffness retain their exact rational dependence on E and nu;
both thermal derivatives contribute to storage, including when expansion is zero.
Dimensionally equal semantic coefficient kinds are not interchangeable with
the native dimension-only elastic schema. External coordinates and selected
claims stay pinned through every refresh.

The original `IsotropicElastic` supplies stress. Explicit conversions between
its tensor-shear Voigt convention and the solid owner's Mandel convention
preserve work. At fixed total strain, the adapter evaluates elastic Helmholtz
energy Psi, entropy -dPsi/dT, internal energy Psi+T*S and fixed-total-strain
capacity -T*d2Psi/dT2. Total and elastic strain components must stay inside the
declared small-strain envelope. The stress-free reference is the retained
expansion source, not an inferred density or material-name temperature.

The full-interval conservative lower capacity bound uses outward-rounded
nonnegative arithmetic and the complete six-component Mandel envelope. It
includes stiffness curvature and the instantaneous expansion slope, and is
intended to be added to the actual caloric inventory by a coupled consumer. This adapter
does not create heat storage, infer current density, evolve temperature, solve
mechanical equilibrium, supply plasticity or claim finite-deformation validity.
A consumer must separately retain its reference mass, admit the actual caloric
chart and check positive combined capacity.

Focused regression definitions cover independent isotropic invariants,
pure shear work, fixed-total-strain thermal derivatives, signed modulus slopes
and expansion, simultaneous affine E/nu/alpha changes, exact relaxed/constant-source
limits, capacity bounds and source/schema/domain refusals. Variable nu with zero
expansion retains its nonzero stiffness-curvature capacity. Source paths with
kinks or nonlinear coefficients, anisotropy, phase changes and finite strain
remain outside this law. Pinned rustfmt parsing, focused diff checks and bounded
physics/API review completed; native compilation and tests for the extension
remain unexecuted because this worker exhausted shared memory/storage in its
unchanged runtime dependency.

## Non-isothermal polymer welding kinetics

`healing::kinetics::ShiftedMobility` consumes the original
`visco::temperature::TemperatureShift` and an explicitly supplied reference
mobility, using `r(T) = r_ref/a_T(T)` for both Arrhenius and WLF calibrations.
The dimensionless exposure enters the existing `ContactHistory`,
`StrengthMap` and deposited-contact observer; no second history or material
table is introduced. Seppala et al. (2017), doi:10.1039/C7SM00950J, equation 5,
motivates the reduced-time integral, not any default material coefficients.

Isothermal increments use exact rate-times-duration arithmetic. Changing
temperature follows the consumer's actual linear ramp, with allocation-free
adaptive Simpson quadrature and explicit absolute/relative dose tolerances,
evaluation count and recursion-depth ceilings. The returned error is a
numerical estimate, not a certified bound or calibration uncertainty. An
optional authored zero-mobility cutoff is split before integration; it cannot
extend the original shift domain or imply a glass-transition temperature.
Invalid paths and exhausted quadrature budgets refuse the whole increment.

Focused tests compare both models with independent high-precision integral
fixtures, including a steep WLF ramp, and exercise reference-rate orientation,
cutoff crossings, temperature-path subdivision, retained exposure and bounded
refusal. The native deposition tests carry the original finite-cell thermal
ramp into dose and cohesive strength, with atomic domain refusal. Contact
area creation, crystallization, degradation, calibration and concurrent
development of an already damaged mechanical bond remain outside this model.

## Dose-dependent opening and shear adhesion

`healing::strength::StrengthMap::new_mixed_mode` maps one original contact's
accepted mobility dose to independently supplied normal/shear peak traction and
mode-I/mode-II fracture work. The original `BilinearCohesive::new_mixed_mode`
admits every developed knot and interpolated output, retaining quadratic
initiation, frozen-onset BK toughness, irreversible damage and failed-bond
compression. The existing common-mode map and `FrozenBond` handoff remain
available to the same deposition and mechanical consumers.

All four observations increase monotonically with dose; all-zero knots can
declare an induction interval. Missing mode data, invalid envelopes, a different
mobility identity and extrapolated dose refuse. Fixed stiffnesses and BK
exponent are explicit inputs. Lookup is logarithmic in at most 1024 knots,
without numerical integration. No calibration, intimate-contact area, strength
ratio or concurrent healing/fracture is inferred. The output starts a separate
mechanical stage and cannot replace an existing damaged point.

Focused G1/G3 regressions cover pure-mode fracture work, independent BK work,
in-plane shear rotation, induction/cooling/reheating, complete-failure reclosure,
source/domain refusals and the existing deposition consumer's area partition.

## Source-defined bilinear cohesive interfaces

`cohesive::source` projects an original ordered `InterfaceSystemCard` or an
already resolved multiphysics interface bundle into the existing
`BilinearCohesive` law. Explicit `cohesive_stiffness` [Pa/m],
`cohesive_strength` [Pa], `cohesive_fracture_energy` [J/m²], and
`cohesive_compression_stiffness` [Pa/m] are required. The numerical maximum
separation is supplied separately and must include complete failure. No bulk
material name, absent property, or geometric proximity creates an adhesive bond.

Once complete failure has been accepted, opening and shear may continue beyond
that obsolete adhesive range. Only the independent normal penetration remains
bounded; every component must still be finite. The failed path evaluates the
original compression penalty without inactive adhesive energy norms, preserving
fracture loss exactly and adding no new bond work. The same rule applies to
endpoint, timed, mean-traction/tangent and temperature-transfer paths. A step
starting with surviving bonds still enforces the full separation/rate domain.
Fixed frames, original trace correspondence and the structural owner's own
displacement limits remain in force; no new sliding-contact geometry is inferred.

The full bundle, exact source selections, ordered surfaces and texture frames
remain accessible. The canonical fixed-external-condition policy applies to
both value dependence and applicability: evolving gap, damage, stress, rate or
frequency coordinates cannot be frozen into this rate-independent law. The
ordinary binding keeps external conditions fixed. An explicit temperature
successor uses `at_temperature` and the history transfer below;
there is no implicit coefficient refresh or inference of fracture heat.

The existing common-mode cohesive law still owns traction, tangents,
irreversibility, compression after failure and fracture-energy accounting.
This is not a calibrated mixed-mode adhesion model. The focused
`cohesive_source` regressions cover shared electrical/mechanical source bundles,
fixed-temperature substitution, independent fracture work and source/unit/domain
refusals; structural source-to-force tests live in `fs-solid`.

### Cohesive energy increments under retained preloads

Common-mode and frozen-onset BK trials calculate recoverable work component by
component while damage is fixed. Damage loading uses independently factored
envelope-storage and fracture increments. An unchanged compressive preload,
held orthogonal shear or large accumulated fracture inventory therefore cannot
erase a smaller representable change in another physical channel. Signed
unloading, compression after rupture and the original analytic path-mean
traction retain the same constitutive law and absolute energy gate.

At a fixed-jump temperature transfer, coefficient differences supply parameter
work and the lost intact fraction supplies new fracture release. Unchanged
compression is excluded from that release. Absolute endpoint energies remain
diagnostics and can legitimately round to the same value while an incremental
account is nonzero. The focused `cohesive::energy_tests` specify common/BK
softening and unloading, orthogonal motion, historical loss, thermal weakening
and separate physical viscosity under retained preloads. These are numerical
model checks, not new material calibration or a change to the fracture law.

### Irreversible cohesive history through temperature changes

`CohesivePoint::transfer_temperature` preserves the original jump, retained
damage and accumulated fracture loss while adopting an explicitly supplied
law from the same common- or mixed-mode family. It first changes coefficients at
fixed damage, reporting the signed change in stored surface energy as parameter
work. It then relaxes onto
the new envelope at fixed jump, recording any additional released storage as
nonnegative damage dissipation. Cooling and increased toughness cannot heal.
The inverse damage map changes internal onset/failure coordinates without
rescaling past fracture work. Subsequent loading, viscosity, mean traction and
frictional compression still use the original constitutive owners.

`ResolvedBilinearCohesive::at_temperature` retains the same ordered interface
card and exact selected claims, including explicit viscosity when selected.
`require_temperature_successor` checks the original source, property selections,
temperature meaning, all held coordinates, envelope family and separation/rate
bounds before a native owner may transfer history. The original model selection
cannot switch even if the same bundle contains both sets of coefficients.

For mixed-mode thermal transfer, `MixedModeHistory` retains the unit PHYSICAL
active-separation direction at first damage initiation and its original energy
shear fraction. The current-temperature envelope is evaluated on that direction,
never on a subsequently rotated gap. The new quadratic strength surface fixes
onset; `beta(T) = Kt |et|^2 / (Kn en^2 + Kt |et|^2)` and the original selected
`GIc(T)`, `GIIc(T)` and BK exponent determine current toughness. The explicitly
separate `envelope_shear_fraction` exposes changed energy weighting under
anisotropic stiffness changes without overwriting the physical initiation mode
or original beta. Uninitiated points retain no mode until a mechanical or
fixed-jump thermal strength crossing. Damage and old fracture work remain
irreversible under heating, cooling, unloading and changed mode toughness.
Parameter work includes both anisotropic cohesive storage and undamaged
compression; failed bonds retain compression without regained adhesion or
viscosity. Ordinary fixed-temperature mixed-mode behavior is unchanged.

This is a declared staggered endpoint model with prescribed temperature,
not a cure/fatigue law, continuous thermal-path validation, heat partition or
entropy proof. Parameter work belongs to the external thermal/material source;
additional damage is a separate physical loss. Focused G1/G3 tests cover
analytical fixed-jump work/release, unloading and changed-Gc continuation,
failed-bond compression, no healing and changed-envelope refusal. Five mixed-mode
regressions add thermal initiation, retained physical mode and original beta,
anisotropic heating/cooling, independent work quadrature through rotating reload
and rupture, pure opening/shear reductions and failed viscous compression.
Source regressions exercise the original timed three-component point with and
without physical viscosity. The existing native series joint is a pure-opening
consumer; this extension does not create a shear-equilibrium joint solver or
establish measured thermal mixed-mode calibration.

## Explicit independent isotropic elastic parameters

The native `IsotropicElastic` stress evaluation resolves the bulk and shear
modes directly from the original E/nu pair. A compensated trace and compensated
mean preserve weak volumetric strain among cancelling normal components and
avoid spurious deviatoric stress under uniform expansion. This retains both
hydrostatic response near the auxetic limit and trace-free response near
incompressibility without multiplying a dense stiffness matrix. Tensor shear
remains `sigma_xy = 2 G epsilon_xy`; the constant tangent retains its existing
representation. G0/G1/G3 tests in `elastic.rs` cover both admitted Poisson-ratio
limits, stress-unit scaling, axis permutations and the work-conjugate energy
derivative. This improves constitutive evaluation; it does not remove volumetric
locking or conditioning limits from displacement-only element formulations.

`state_point::resolve_isotropic_elastic_state_point_from_parameters` admits a
caller-selected Young/Poisson, Young/shear, or bulk/shear pair with density.
Arbitrary source property names and complete typed query conditions are retained.
The scalar role requirements reject stress, pressure, hardness and tensor-component
aliases; this adapter does not infer isotropy from a material name or directional
measurements. Positive, finite and numerically representable isotropic stiffness
is required before publishing the existing `IsotropicElasticStatePoint`.

The original resolver's `resolved()` bundle contains the selected source properties.
`derivation()` retains their roles and nominal algebraic consequences separately;
no derived constant is inserted as a measured claim. `identity()` binds the
derivation and is used by the solid element, thermal-strain and reduction ingress.
Canonical Young/Poisson states retain their previous identities. Source marginal
uncertainty remains attached; joint covariance is not assumed and uncertainty is
not propagated. Additional unselected measurements are not reconciled automatically.

`ResolvedIsotropicElasticModuli::try_from_resolved` projects that same checked
algebra from an already resolved complete bundle without requiring density.
Solid stiffness updates therefore retain the original reference mass inventory
instead of asking for a new `rho(T)`. All thermal, electrical, optical or other
properties in the supplied bundle remain present, with their exact original
receipts. Source keys and consumer domains are checked again at adoption;
finite-only resolution cannot bypass elastic stability or numerical
representability. The projection identity binds source state, parameter roles,
projection domains and nominal moduli. No additional source claim is created.

`IsotropicElasticStatePoint::{try_from_resolved,
try_from_resolved_with_parameters}` similarly adopt complete bundles for
consumers that need density. Both elastic projections expose `at_query_point`,
which refreshes every retained property through the existing shared resolver,
pins the selected claims and rebuilds the same constitutive parameterization.
An auxiliary property leaving support refuses the whole candidate, even if
the elastic constants remain available. A stricter projection-only domain
also remains enforced after refresh. Original canonical and independent-pair
resolver identities remain unchanged; additional projection-only domains use
the explicit parameter identity v2. This is endpoint evaluation, not a mass
transfer, path-support proof, constitutive-history integration or heat solve.

`IsotropicSolidStatePoint::{try_from_resolved,
try_from_resolved_with_parameters}` adds an explicitly admitted positive yield
source to that same elastic projection. The complete bundle and its auxiliary
properties survive adoption and `at_query_point`; elastic constants are never
re-resolved through a different pair. Its physical `identity()` binds elastic
and yield roles, and the ordered Hertz reduced-modulus identity carries both
solid identities into contact and specimen consumers. Canonical solid identities
remain unchanged. The reciprocal `ThermalHertzPair` differentiates the explicitly selected
Young/Poisson, Young/shear or bulk/shear source laws. It preserves the complete
solid projection, including named density and yield roles and auxiliary
properties, when querying each surface at its own temperature.
`ResolvedMaterialStatePoint::with_additional_requirements` adds disjoint source
properties at the same immutable card and typed query without re-querying any
existing property. Contact can therefore append yield while preserving every
original elastic/loss receipt, including deliberately pinned competing claims.

Focused regressions cover independent constitutive stress, source conditions,
invalid pairs, conflicting claim pins, unit rescaling, complete-bundle adoption,
multiaxis refresh and atomic auxiliary-property refusal, and the actual tetrahedral
mass/strain-energy and thermal-expansion identity path. This is a bounded MR11/MR08
parameter-conversion slice, not complete derived-property or material validation.

### Shared isotropic parameters in concentric sections

`ConcentricElasticSection::try_new_with_isotropic_parameters` consumes one
`LayerIsotropicParameters` binding per layer through the existing isotropic
projector. Density and E/nu, E/G or K/G may use arbitrary source keys; each layer
retains its complete original bundle and auxiliary properties. Derived E and G
remain nominal model consequences, never inserted measured claims. Section
mass, axial/bending stiffness, transport, strength and loss adapters consume
those same coefficients and their actual source dependencies.

`with_isotropic_shear` uses G from the original pair. `shear_property()` returns
None when shear follows that pair and Some for an explicitly named G source.
On an explicitly isotropic section, named G must select the original G requirement
of E/G or K/G; E/nu must use `with_isotropic_shear`. A separate, numerically equal
G claim cannot replace original source support or thermal derivatives. Legacy
density/E sections still accept independently named G. Source-aware accessors
retain the parent receipts and frequency/constancy restrictions.
`with_material_layers` and shear `with_material_section` preserve exact layer
geometry and parameter roles, rejecting incomplete or inadmissible candidates.
They do not perform mass transfer, history integration or geometry evolution.

Ordinary section temperature refresh re-evaluates every retained property with
its original claim pins, including nonlinear E derived from varying K/G. The
existing affine thermoelastic potential admits scalar K/G only; varying K/G
refuses because its derived E generally has additional thermal derivatives.
`ConcentricRodThermoelasticLaw::try_new_with_isotropic_shear` also admits scalar
or affine E/G, and scalar or affine E with scalar nu. Scalar reference density,
expansion, source-path support and the original strain limits remain required.
Focused G0/G1/G3 regressions specify mixed-parameter layer integrals, coupled
transport/loss response, refresh refusal and potential derivatives. They do not
establish measured-material validation or propagated joint uncertainty.

### Reversible adhesion driving equilibrium liquid wetting

`liquid::surface::adhesion::resolve_young_dupre_wetting` derives the existing
`ResolvedLiquidWettingState` from an explicit ordered interface claim
`equilibrium-work-of-adhesion` [J/m2] and the selected liquid/gas surface
tension [N/m]. The explicitly chosen ideal smooth, homogeneous, rigid-wall
model obeys `W = sigma (1 + cos(theta))`; the half-angle evaluation remains
well conditioned near zero angle. The model requires `0 < W < 2 sigma`.
Complete wetting, film formation, complete dewetting, roughness and hysteresis
require other physical owners and are refused rather than clipped.

The angle retains both original source receipts and is identified as a
Young--Dupre prediction, never inserted as a measured contact-angle claim.
The boundary checks the exact liquid card, selected tension receipt, ordered
surface roles, query point and axis semantics before the existing native LBM
wall operator consumes it. Refresh requires the original interface and the
original selected tension source at the same new physical point; no source
substitution or partial refresh is accepted. This is a frozen endpoint, not a
thermal surface-energy evolution. Directly sourced contact angles retain their
existing resolution and refresh path.

`YoungDupreLaw` also supplies the signed dry-substrate-reference surface free
energy `sigma A_lv + (sigma-W) A_sl` and unbalanced Young force per contact-line
length. These fixed-size evaluations allocate no memory. Reversible equilibrium
work is distinct from fracture toughness, peel dissipation and advancing or
receding hysteresis measurements; none are silently substituted. No droplet
pull-off force is inferred. The relation and energy convention follow
Nagy, Scientific Reports 14, 29991 (2024), equations (1)-(3),
<https://doi.org/10.1038/s41598-024-81710-6>.

Numerical regressions check the Young relation, force direction, independent
fixed-volume spherical-cap energy minimization, area scaling and domain
refusals. Native source/wall tests exercise distinct adhesion-driven capillary
responses, mass conservation, exact paired-source refresh, substitution refusal
and the unchanged directly sourced angle path.

## Mode-dependent cohesive fracture

`BilinearCohesive::new_mixed_mode` selects separate opening/shear stiffnesses,
peak strengths and fracture energies. Quadratic nominal-traction interaction
localizes the first damage-initiation point on each straight jump segment. The
undamaged shear-energy fraction at that point selects the Benzeggagh--Kenane
toughness. This implementation explicitly freezes that onset mix and envelope
for the point's subsequent history. A shared maximum energy norm then drives
irreversible scalar degradation, with undamaged compression after failure.
Pure-mode and proportional mixed-mode work equal their selected fracture energy;
nonproportional post-initiation paths retain exact analytic work/storage/loss
accounting and cannot reselect toughness or heal. The onset-mix freeze is a
declared model extension, not experimental qualification of arbitrary histories.

Admission compares the BK toughness with the quadratic initiation energy of
that same mode over the complete shear-energy fraction interval. In particular,
a tough, high-strength shear branch no longer causes a false refusal through an
unphysical combination with the opening branch's smaller toughness. The ratio
`Gc(beta) * ((1-beta)/EI0 + beta/EII0)`, with `EI0 = Tn^2/(2Kn)`, must remain
above one. Endpoint limits and the sole possible interior minimum are checked
with scaled arithmetic and a bounded deterministic derivative solve at law
construction; no mode sampling or extra work enters mechanical increments.
Unresolved near-zero softening intervals still refuse. G1/G3 regressions cover
strong mode anisotropy, exact proportional fracture work, narrow inadmissible
interior modes, reversed strength ordering and coherent energy rescaling.

The consistent algorithmic tangent includes the sensitivity of the initiation
crossing and selected toughness when a trial first initiates damage. It can be
nonsymmetric; the existing Newton/FGMRES owner consumes it without an SPD
assumption. Evaluation allocates nothing; initiation and timed damage/failure
events use analytic quadratic crossings, split at the normal compression switch.
Radius-scaled coefficients and a rationalized outward root retain initiation
fractions below `2^-64`, preventing a long increment from changing the physical
strength threshold or causing premature rupture. Retained-envelope steps
use closed-form arithmetic. Equal mode
parameters recover the original common-mode response. The scalar law getters
continue to describe pure opening; each point exposes its selected mixed envelope.

`CohesivePoint::trial_with_mean_traction` additionally returns the analytic
three-component traction average along that same straight jump path. It splits
retained-history elasticity, softening, compression and complete failure using
the original initiation/envelope history. Affine-vector/norm integrals preserve
transverse traction under nonproportional loading; a short-interval binomial
evaluation retains small tractions near failure. The mean's dot product with the
full jump change must reproduce the original independent path work within the
supplied J/m2 gate. No scalar work projection, quadrature, new damage law or
allocation is introduced. Endpoint traction and endpoint tangents remain separate
and do not differentiate the path average. Regressions compare vector means with
independent traction quadrature, shear rotations, full rupture, reversals and
straight-path subdivision.

The optional `trial_with_mean_traction_tangent` and
`trial_over_time_with_mean_traction_tangent` methods return that unchanged trial
and vector mean together with `Option<[[f64; 3]; 3]>` for
`d(mean traction)/d(endpoint jump)` [Pa/m], holding accepted history and the
physical interval fixed. The existing `fs-ad::Real`/`Dual64<3>` owner differentiates
the shared analytic path primitives and the implicit locations of damage,
failure and normal-contact crossings. Common-mode and mixed-mode evolution
include first initiation, unloading, softening and rupture. A first BK selection
reuses the original quadratic-onset sensitivity to differentiate its energy
radius and fracture energy, including both moving softening/rupture crossings
and the selected envelope. Previously selected BK histories stay frozen. The
original onset lower bound also covers a strength crossing that precedes
renewed outward energy-radius growth. Physical viscosity includes the derivative
of the integrated surviving-bond fraction,
not only its direct rate term. Original source, damage, work and dissipation
accounts are unchanged; ordinary mean-only methods execute no tangent path.

`None` explicitly leaves the consumer's numerical/semismooth fallback in charge
for an identically zero-normal path, a held onset/loading-unloading kink, or a
degenerate/unrepresentable derivative. The
endpoint tangent is never substituted, and a selected one-sided derivative is
not represented as a classical one. Allocation and work remain bounded per
point. G0/G1 regressions compare independent endpoint perturbations of the
original mean and total work, first-initiation BK gradients across integer and
fractional exponents, shear-basis covariance, delayed energy-radius growth,
retained-history symmetry, exact elastic/viscous slopes, coincident
initiation/reclosure, tiny-onset rupture and short-softening limits. An
independent fully ruptured mixed-mode solution checks the complete toughness
gradient, including its generally nonsymmetric contribution. These are numerical
constitutive tests, not adhesive calibration or
wall-clock performance claims.

`trial_over_time_with_damage_weighted_compression_tangent` also returns the
unchanged timed trial, vector mean and original damaged-area compression
`integral d(s) Kc max(-gn(s),0) ds`, with separate optional derivatives of the
mean traction and compression measure with respect to the endpoint jump.
The pressure gradient reuses the analytic pressure primitive and moving
initiation/failure/contact crossings, including first BK mode selection; it
does not freeze endpoint damage or introduce another fracture history.
Accepted history and physical time are fixed. Bond viscosity enters the
traction derivative, never the compression-pressure measure. Ordinary scalar
callers execute no derivative work; this path is bounded and allocation-free.
Held kinks or unrepresentable crossing derivatives return `None` independently
without changing admitted forces, work or state. Focused G0/G1/G3 regressions
check independent derivatives through rotating mixed modes, reclosure, complete
rupture and a small newly damaged area, plus exact intact/failed limits and
source-preserving refusal.

`resolve_mixed_mode_bilinear_cohesive` requires eight explicit ordered-interface
claims: normal/shear stiffness [Pa/m], normal/shear strength [Pa], mode-I/mode-II
fracture energy [J/m2], dimensionless BK exponent, and compression stiffness. No
missing mode is inferred. The existing `ResolvedBilinearCohesive` source bundle
and solid/P1/thermal point paths remain the owners; fixed external source
conditions and all original source-role refusals apply. The existing thermal
joint is normal-only and therefore consumes the pure-opening envelope through
the same exact point work/storage interface.

Quadratic initiation and BK interaction reference Camanho and Davila (2002),
NASA/TM-2002-211737, <https://ntrs.nasa.gov/citations/20020053651>. Tests check
pure-mode and proportional mixed fracture work, equal-mode recovery,
initiation-sensitive and retained-history tangents, nonproportional work,
straight-path subdivision, compression/shear, source-unit/domain refusals and
actual native-equilibrium load histories. Temperature/rate-dependent toughness,
fatigue, evolving mode-mix calibration and fracture heat partition remain separate.


### Explicit power-law mixed-mode fracture

`BilinearCohesive::new_power_law_mixed_mode` selects the power-law envelope
`(GI/GIc)^alpha + (GII/GIIc)^alpha = 1` with `alpha >= 1/2`, independent
of the existing BK constructor. `MixedModeCoefficients` supplies the six
mechanical coefficients; `mixed_mode_criterion()` identifies the selected
family and exponent. The legacy `mixed_mode_parameters()` getter remains
BK-specific. Generic consumers use `mixed_mode_coefficients()`.

The original quadratic initiation, frozen physical onset direction,
irreversible damage, analytic path work, vector traction, physical viscosity
and thermal transfer remain the numerical owners. Separate scaled normal and
shear amplitudes retain weak-mode toughness contributions even when the
reported shear-energy fraction rounds to one. First-initiation endpoint and
path tangents include the selected power-law toughness gradient. At
`alpha = 1/2`, pure-mode initiation has a cusp: classical derivatives
explicitly decline, while the endpoint branch uses a zero transverse
subgradient. An already accepted initiation envelope remains frozen.

Admission checks positive softening over every mode, not only the pure limits.
Writing `a = Tn^2/(2 Kn GIc)`, `b = Tt^2/(2 Kt GIIc)`, the condition is
`max(a,b) < 1` for `alpha >= 1`, or
`||(a,b)||_(alpha/(1-alpha)) < 1` for `1/2 <= alpha < 1`, with a
roundoff guard. The new power-law constructor also requires the sufficient
all-mode failure-distance bound
`2 Gmax / (min(Tn/sqrt(Kn),Tt/sqrt(Kt)) sqrt(min(Kn,Kt)))`
to fit the declared separation range. `Gmax = max(GIc,GIIc)` for
`alpha <= 1`; for `alpha > 1` it is the dual norm
`||(GIc,GIIc)||_(alpha/(alpha-1))`. This conservative admission covers
rotation after initiation; legacy BK admission is unchanged.

`resolve_power_law_mixed_mode_bilinear_cohesive` and
`CohesiveEnvelope::FrozenOnsetPowerLaw` require the original dimensionless
`cohesive_power_law_exponent` claim alongside all seven independent
mechanical interface properties. A BK exponent is not a substitute.
Temperature refresh preserves original claims, viscosity/range bounds and
the selected family; an explicit power-law exponent change reevaluates the
envelope on the retained physical onset direction without healing.
Common/BK/power-law family changes refuse even when the complete source
bundle contains every alternative.

The fracture-envelope equation follows Neves et al., Processes 10 (2022) 2730,
Eq. (5), <https://doi.org/10.3390/pr10122730>. Reported adhesive envelopes
include `alpha = 1/2`; the exponent is not a BK calibration. The tests use
manufactured coefficients and independently check proportional work,
initiation-sensitive endpoint/path gradients, weak modes, admission,
viscosity and thermal history. They do not supply a measured adhesive
calibration or validate arbitrary nonproportional experiments.

### Intact viscoelastic adhesive layers

`cohesive::relaxation::IntactViscoelasticBond` supplies small-strain, fixed-frame
normal/shear adhesive relaxation on an explicitly intact, massless thin layer.
For original thickness h, each local strain is its physical jump divided by h.
Independent effective normal and shear `GeneralizedMaxwell` spectra are required;
the normal spectrum already represents the specified joint's lateral restraint.
No Poisson-ratio assumption turns bulk E/K into interface normal data. Both shear
directions share the supplied shear spectrum. Compression remains the same intact
material; debonding, contact, friction and healing are outside this law. It cannot
be attached to a damaged `CohesivePoint` or reset its history. Temperature-dependent
coefficients use the explicit held-jump continuation below.

The original `PronyState` and exact `PronyRamp` own memory and physical heat.
`prepare(dt)` evaluates exponential ramp coefficients once per law/duration,
reusing them across trials and both shear directions. `ViscoelasticBondPoint::trial`
returns a new immutable state, endpoint traction, analytic vector mean traction,
and independent incremental work, storage and heat per original bonded area.
Its local accounts do not subtract accumulated inventories. Branch storage uses
the prepared `PronyRamp::spring_increment`, retaining weak relaxation release
even when the endpoint spring strain rounds to its accepted value. Independent
analytic held-jump checks resolve storage and heat down to `dt/tau = 1e-17`.
Held jumps relax; straight-ramp subdivision preserves the constitutive result.
Endpoint and mean stiffnesses are distinct exact derivatives at fixed accepted
memory and duration.
A prepared ramp from another law refuses even if its numeric coefficients match.

There are at most 64 branches per mode. Authored jump-norm, jump-rate and elapsed-time
bounds are enforced without truncation; those model-use bounds are not inferred
experimental qualification. Source binding requires the original ordered interface,
`cohesive_bondline_thickness` [m], and separate
`cohesive_normal_relaxation`/`cohesive_shear_relaxation` property families:
`_branch_count` [-], `_equilibrium_modulus` [Pa], and one-based indexed
`_modulus_N` [Pa]/`_time_N` [s]. Every declared branch and the complete original
query/receipt bundle are retained; missing modes, incompatible units and evolving
coefficient inputs refuse. The normal modulus is an effective normal joint response,
not an automatic alias for bulk Young or bulk compression modulus.

The existing `fs-couple::cohesive_contact::CohesiveLoadLine::trial_relaxation`
consumes the affine work-mean traction through its original area, complete 3x3
mechanical compliance, prepared LU workspace and physical/compatibility gates.
It needs one linear solve and one material trial, without finite differences or
a new motion solver. Actual free motion and the complete energy window remain
with the mechanical caller; temporal refinement still governs the approximation
of an inertial trajectory by straight jump ramps.

Constitutive background: Abaqus, Time Domain Viscoelasticity, traction-separation
elasticity, <https://docs.software.vt.edu/abaqusv2026/English/SIMACAEMATRefMap/simamat-c-timevisco.htm>.
Focused numerical tests check independent relaxation formulas, the original scalar
Maxwell owner, reversals, subdivision, shear covariance, tangents, retryable domain
refusals, cross-axis compliance and bonded-mass work/heat closure. Fixtures are
manufactured coefficients, not a measured epoxy or silicone calibration.

### Held-jump temperature continuation of adhesive relaxation

`ResolvedViscoelasticBond::at_temperature` requeries the same selected coefficient
claims on an original, dimensionally declared Kelvin axis. Its
`require_temperature_successor` guard preserves the original ordered surfaces,
card identity, complete claim plan, branch counts and indexed identities,
thickness, model-use limits, and every non-temperature query coordinate.
A coefficient curve cannot silently change geometry or add/drop/reorder memory
branches. Claim-domain refusals remain with the original material query owner.

After that source guard, `ViscoelasticBondPoint::transfer_temperature` retains each
physical Maxwell spring strain, total jump, elapsed time, past mechanical work
and physical heat. Only spring coefficients and subsequent relaxation kinetics
change. At fixed strains it accounts for signed parameter work
`h/2 (delta(E_inf) epsilon^2 + sum delta(E_j) q_j^2)` separately from mechanical
work and heat. An independent spring-force storage-change calculation uses
compensated product differences; unchanged springs contribute exact zero.
Small changed-mode work therefore survives an unrelated large stored inventory.
A change of relaxation times alone has exact zero parameter/storage work and
leaves endpoint traction unchanged. Relaxation times may cross without sorting
the declared branches or replacing their internal state.

The returned `ViscoelasticBondTemperatureStep` contains the successor point,
parameter work, storage change and their checked defect per original area.
Every failure leaves the original point retryable. The caller prepares the
successor law for its next actual mechanical interval and owns the signed
parameter-work destination. This is a coefficient continuation with held
geometry and physical strain, not a caloric, entropy, thermal-expansion,
fracture or thermorheological-simplicity calibration. No heat is inferred from
the temperature change itself. Focused tests cover memory preservation,
heating/cooling reversibility, kinetics-only continuation, branch-time crossing,
small changed work beside large storage, and native load-line continuation.

## Solid-vapor constitutive laws

`sublimation` supplies the shared `fs_thermal::sublimation::SublimationLaw`
interface. `WaterIceSublimation` evaluates the IAPWS ice-Ih pressure curve on
50..=273.16 K; `CarbonDioxideSublimation` evaluates the NIST solid-CO2 Antoine
fit on 154.26..=195.89 K, with explicit bar-to-Pa conversion. Named curves
refuse extrapolation. Pressure derivatives are analytic. Their temperature-
dependent latent heat and derivative follow Clausius-Clapeyron with ideal
vapor and negligible solid specific volume; they are not independently
measured calorimetry or a real-fluid equation of state. Original sources and
approximations are documented beside each law. Density is explicitly supplied
by the caller and held constant; no porosity, accommodation or geometry is
inferred from the material name.

`ClausiusClapeyronSublimation` admits caller-calibrated reference pressure,
reference temperature, constant latent heat, molar mass, density and a bounded
validity interval. Log-pressure evaluation avoids avoidable overflow and keeps
the exact reference anchor. All scalar laws allocate nothing and use the shared
`fs_math` deterministic transcendental kernels. Ice saturation pressure is also
available independently of density for humidity consumers. The six focused
regressions check published pressure anchors, analytic derivatives, units,
latent-heat behavior, explicit domains and representability refusals. These
laws supply constitutive data; phase inventory, transport and energy evolution
remain with their physical consumers.


## Explicit ice-relative humidity

`gas::HumidityReference` and `try_new_moist_air_with_reference` select liquid
water or ice as the relative-humidity pressure reference. Ice uses the same
bounded IAPWS curve as the sublimation law. Existing moist-air entry points
retain their liquid-water reference, including below freezing; zero humidity
retains the original dry-gas domain. Actual vapor mole fraction flows through
the original ideal-mixture, Wilke and WMS evaluator, so acoustic and transport
consumers see the same physical composition. The existing RH [0,1] and dilute
water-mole-fraction ceiling remain in force.

`resolve_moist_air_state_with_reference` retains the reference alongside both
original component source bundles. Ice pressure does not relax either card's
applicability or qualify cold transport fits. These APIs prescribe a gas state;
they do not evolve condensation or change a downstream acoustic standard's
humidity convention. Seven native shared-mixture regressions passed, including
four new reference/domain cases. The source-resolution wrapper was syntax
checked and reviewed, but was outside that native numerical closure.

### Opt-in dielectric switching and uniform-capacitor history

With `experimental-dielectric-hysteresis`,
`field_storage::dielectric::hysteresis` exposes a fixed-condition uniaxial
polarization play model in electric SI quantities. It reuses the scalar return
map and inverse of the existing magnetic numerical owner, with
`D=epsilon_b E+sum(P)`, `W=epsilon_b E²/2+sum(P²/(2chi))`, and physical
switching heat `Q=sum(Ec |delta P|)`. Exact monotone-segment work includes pinning
and saturation crossings. Each reversal must be supplied explicitly; states
retain remanence and minor-loop memory, and trials do not publish history.

`UniformHystereticCapacitor` supplies `q=A D`, `V=d E`, incremental capacitance,
work-conjugate mean voltage, and volume-scaled storage and heat. The source
adapter retains explicit scalar coefficient roles, complete original material
selection and source support. Relative permittivity, changing coefficients and
unowned frequency/field/history applicability cannot be frozen into this law.
The linear/zero-loss limits, analytic loop area, subdivision invariance,
charge inverse, tangent, geometry, source refresh and atomic refusal have
focused regressions in `tests/dielectric_hysteresis.rs`.

This is an explicitly fitted quadratic-potential specialization of an
energy-based dry-friction model, not measured ferroelectric qualification,
microscopic domain dynamics, breakdown, aging, electromechanical strain,
temperature-dependent switching or a spatial electric-field solve. Its circuit
consumer owns time integration. The broader ferroelectric/fractional/nonlocal
bead remains incomplete and this capability is default-off.

### Scalar magnetic hysteresis and uniform-core constitutive coupling

`field_storage::hysteresis::ScalarHysteresis` is an opt-in, fixed-condition
energy-based play model with explicit caller-supplied branch susceptibilities,
pinning fields and saturation polarizations. Each branch owns retained internal
polarization. The immutable law and accepted state stay separate; trial field or
flux updates do not mutate accepted memory or silently replace model parameters.
This supplies remanence, reversal/minor-loop memory and saturation independently
of the existing reversible B-H sources.

For B = mu_b H + sum(z), stored energy is
W = mu_b H²/2 + sum(z²/(2 chi)); irreversible heat is sum(Hc |delta z|).
The exact scalar constrained minimizer gives pinning/sliding/saturation behavior.
For each declared monotone segment, magnetic work is the analytical integral
H dB = delta W + Q, including crossed corners. Every field reversal must be
supplied; endpoint trapezoids do not replace this work. Positive background
permeability makes the B-driven inverse strictly monotone. Safeguarded Newton
uses the active-set tangent, has bounded work, and allocates state only for the
final candidate. Tangent values at kinks select the pinned/saturated side.

`UniformHystereticCore` exposes current- and linkage-driven candidates through
lambda=N A B and I=l H/N, including signed winding orientation. Core energy and
pinning heat scale by the actual A l volume; reported mean current is conjugate
to linkage work. Circuit trials must all use the same accepted magnetic history
until the complete coupled step accepts. This is a constitutive adapter, not an
alternate circuit time integrator or automatic spatial-field hookup.

The branch potential is quadratic with explicit saturation bounds, following
the scalar dry-friction energy construction discussed in Prigozhin et al.,
IEEE Trans. Magn. 52 (2016), doi:10.1109/TMAG.2016.2599143.
It does not reproduce the paper's general vector model or confer measured
material authority on supplied coefficients. Temperature, stress and parameters
remain fixed. Eddy/excess losses, dynamic/vector hysteresis, demagnetization
damage and thermal parameter work require separate models. Never add a loss map
for the same resolved hysteretic heat. Eleven focused tests cover analytic
segments/loop area, minor memory, monotone subdivision, inverse histories,
tangents, reversible limits, refusal and winding energy/heat conjugacy.


### Isotropic vector magnetic hysteresis

`field_storage::hysteresis::vector::VectorHysteresis` extends the original
fixed-condition branch parameters to three-dimensional polarization memory.
Each branch minimizes |z|²/(2chi) - H.z + hc|z-z_old| inside its saturation
ball. Exact interior shrinkage and a bounded scalar multiplier solve give the
constrained vector minimizer. Its consistent active-set dz/dH is symmetric
positive semidefinite; the positive background permeability makes the complete
dB/dH positive definite for the bounded, backtracked B-driven inverse.
All trials start from the same accepted memory. Branch-state storage allocates
only for the accepted field candidate, not for inverse trial responses.

For very small rotating slips around a large retained polarization, rounding
can make the analytic unsaturated return miss the original stationarity gate.
Only that exceptional case receives bounded coordinate Newton correction of
the actual retained-state residual. Each candidate must strictly improve that
residual and remain inside the original saturation ball. The branch iteration
budget and cancellation still bound all work; the final consistent tangent is
recomputed from the retained slip. No physical tolerance is relaxed, and the
original independent work audit still gates acceptance.

Reports distinguish physical pinning heat sum(hc|delta z|), stored-energy
change, endpoint H_new.delta B, and nonnegative numerical dissipation. The
numerical term is not heat; endpoint work is not an exact arbitrary rotating
path integral. Resolve direction changes and refine the history. Original
source/field/flux domains, constitutive and work residuals, bounded iterations
and cancellation gate acceptance without modifying the previous state.
An imposed B above the supported background-plus-saturation envelope returns
typed `FluxRange` with observed and maximum magnitudes [T], allowing a loading
owner to distinguish a retryable numerical boundary lift from invalid inputs.
The dielectric facade preserves this diagnosis as `DisplacementRange` [C/m²].

Eight physical tests cover the scalar state/heat limit, independent spherical
minimization, rotation covariance including tangents, an analytic rotating-loop
loss limit and history refinement, finite-difference tangents, B-driven
histories, reversible limits and atomic refusal. Eleven original scalar tests
remain passing. This is an explicit fitted, fixed-temperature/stress
quadratic-potential specialization of Prigozhin et al. (2016), not experimental
material validation, texture anisotropy, dynamic/eddy/excess loss,
demagnetization damage or automatic field/circuit/thermal integration.

The B-driven spatial interface exposes the frozen-history incremental
potential `W_new + Q` and its consistent reluctivity `dH/dB`, computed from
the same accepted return map. The potential's derivative is the endpoint H;
reversible storage alone is not the Newton merit functional. A stable relative
form uses `delta W + Q` about the original accepted energy. Focused G0/G3
regressions differentiate the potential and H with fixed previous memory in
pinned, sliding and saturated rotating states and check inverse-tangent
reciprocity. These potential values are variational objectives, never new
physical storage or additional heat accounts.

`ResolvedScalarHysteresis::isotropic_vector_law` and
`source::ResolvedVectorHysteresis::from_isotropic_source` explicitly select
isotropy using the original fitted card, coefficient roles and branch order.
The wrapper retains the original source evidence and complete query. Its
constant-fit context refresh preserves magnetic memory and refuses unsupported
conditions or changed coefficients through the existing source owner. Scalar
coupon evidence is not promoted to experimental rotating-field validation.

### Intrinsic finite-rate magnetic switching

The same vector owner optionally accepts one explicit viscosity eta [m s/H]
per branch through `with_viscosities_m_s_h`. Timed H- and B-driven candidates
minimize the original dry-friction potential plus eta|delta z|²/(2 dt), inside
the same saturation ball. The scalar saturation multiplier, full symmetric
active-set tangent and bounded inverse remain native. The implementation forms
the physical switching drive before adding viscous stiffness, avoiding the
large shifted-field cancellation of an effective-coefficient construction.
Storage always uses the original chi. The physical balance reports pinning
heat hc|delta z|, viscous heat eta|delta z|²/dt and numerical dissipation
separately. The B-driven incremental potential contains **half** the viscous
heat; the physical heat receipt contains its full value. Resolved Ohmic loss
remains a distinct owner and must not also be fitted into eta.

Finite holds can relax polarization and release heat. Every trial uses the
same accepted memory and positive finite duration; a changed kinetic profile,
clockless update, unrepresentable stiffness, failed native residual or exhausted
budget refuses atomically. An all-zero profile takes the original numerical
path. This is a first-order overdamped, isotropic, fixed-fit model requiring
time refinement, not a domain-resolved magnetic dynamics or empirical excess
loss model.

`ResolvedVectorHysteresis::from_isotropic_source_with_viscosity` selects the
additional original-card viscosity keys explicitly, retaining every existing
receipt unchanged. Repeated keys share a source coefficient while preserving
branch multiplicity. The combined material bundle retains all kinetic claims,
including zero coefficients, through complete-query refreshes and full-interval
temperature admission. Missing, negative, variable, dynamically conditioned or
dimensionally incompatible viscosity sources refuse. Synthetic native tests
cover relaxation, rate dependence, constrained rotation, tangent/potential
derivatives, physical energy accounting, zero-profile recovery and original
source validity; their authorship is not experimental material qualification.

### Compressible porous leakage and selective membrane permeation

`porous::permeability` separates intrinsic permeability [m²] from gas viscosity,
dense-membrane molar permeability [mol/(m·s·Pa)] and thickness-normalized
permeance [mol/(m²·s·Pa)]. The porous owner compiles orthotropic normal
permeability and series layer resistance once. Its isothermal ideal-gas Darcy
law conserves mass with distinct face volume rates, includes an explicitly
gas-model-bound per-layer Klinkenberg factors when calibrated, and refuses unsupported
temperature, absolute pressure, speed or material conditions. A zero normal
permeability seals the path exactly. The acoustic gas-substitution owner uses the
same intrinsic-permeability type. This is dry rigid-skeleton quasi-steady flow;
Knudsen/transition diffusion, pore storage, wet multiphase flow, pore heat storage,
Forchheimer inertia and evolving deformation are not inferred.

`thermochemical_gas::energy::parcel::transfer::species::network::transport::step_candidate_with_paths` accepts `GasTransportPaths` containing nozzles,
Maxwell–Stefan diffusion, porous barriers and dense membranes. All rates sample
the same original inventory and enter one conservative species/partial-enthalpy
exchange. Gross outgoing species allowances are shared across every mechanism;
incoming gas cannot replenish those allowances or traverse another edge during
the same explicit step. Final endpoints must preserve active total-pressure or
species partial-pressure order and pass the existing extensive entropy gate.
The original two-mechanism API delegates with empty porous/membrane slices.

Material parcel and gas-field Maxwell–Stefan transport share the original
`fs-thermochem::transport::maxwell_stefan` algebra. The default parcel entry
points retain isothermal/isobaric molar-averaged controls, original pair material
samples, absolute residual budgets and selective partial-enthalpy/entropy checks.
The shared owner polls during scaled assembly/elimination and retains its
original unscaled equations and reference-frame checks.

The explicit `GasDiffusionDriving::CompositionAndPressureMass` mode evaluates
`grad(x_i) + (x_i-Y_i) grad(ln p)` in the MASS-average frame with an independently
supplied kg/(m² s) flux-sum tolerance. Both endpoints and the mean face must lie
within the same selected binary-diffusivity claim. The original deterministic
log-ratio helper preserves small pressure differences. Reports identify the
mode, both actual flux sums, pressure gradient and pair-friction entropy rate;
only the selected frame sum is constrained. With zero net mass flux the force
contraction is the full ideal-gas chemical/pressure entropy driving, so this
mode retains the actual closed-pair/network entropy gate without an omitted
mechanical entropy term. Partial species enthalpy still crosses exactly once.
`diffusion_flux_to_with_driving`, `diffuse_with_driving`, and the fixed/adaptive
network `*_candidate_with_driving` entry points share these original owners.
All existing default APIs remain molar-averaged. Neither mode evolves bulk
momentum or supplies Soret/Dufour coupling; declared local T/p limits still apply.

`species::diffusion::network::adaptive::advance_candidate` completes a requested
positive molecular-diffusion duration on branching, cyclic or disconnected
fixed-volume networks. One original full snapshot step is compared with two
original half steps; only the two-half-step path can be accepted. Every probe
resamples actual face T/p/composition, selected pair material data and original
partial enthalpies. The method remains first order; the full/fine difference
estimates local temporal error, not a certified whole-window error. It does not
extrapolate a Richardson state, clip trace species or invent transport data.

Explicit controls bound per-parcel species-inventory error (absolute kg plus
relative scale), original-source endpoint-temperature error (K), trial-interval
size and total snapshot evaluations, including rejected coarse/fine probes.
The minimum interval is a retry floor; an exact shorter final remainder may
be attempted, and each accepted interval contains two physical half steps.
Only typed donor-withdrawal refusal, excessive local error and endpoint entropy
decrease can trigger halving. Missing coefficients, source/domain/caloric errors
and cancellation fail the whole window immediately. The existing molecular
isothermal/isobaric, molar-frame and source-support restrictions remain intact.

Every accepted half step retains the existing simultaneous gross-withdrawal,
species, original-caloric energy and entropy gates. Arrivals can move onward
only during later physical substeps. Rejected probes contribute to work counts
but never inventory, integrated mass or transported enthalpy. Canonical physical
edge/species accounts accumulate accepted gross masses and partial enthalpies;
only the last half step's detailed source receipts are retained. Compensated
accepted energy residuals use the original reference-cancelled convention, with
inventory roundoff checked separately. Final-minus-original actual total/species
mass, this energy ledger and whole-network source-qualified entropy must fit
the same absolute physical budgets, without multiplication by step count.
Late refusal publishes no partial map, including after already accepted probes.

Focused fixtures compare complete-window binary relaxation and a branched mode
with independent exponentials, tighten time accuracy, permute link order/IDs,
shift species formation datums, distinguish rejected probes from accepted
transport, and exercise source/work/accuracy/late-cancellation refusal. Closed
faces preserve parcel states bit for bit while still completing requested time.

The membrane owner uses source-selected, constant-at-fixed-temperature molar
permeability and actual species partial pressures, allowing counterpermeation.
It checks the selected claim at both gas boundaries, retains explicit zero
permeability for blocked species, and refuses missing species, stale temperature
aliases, wrong units and unsupported pressure-dependent sorption. The closure
assumes a steady negligible-storage, near-isothermal dense film with dilute Henry
sorption; competitive sorption, plasticization, swelling, transient breakthrough,
nonideal fugacity and thermal transpiration need additional constitutive laws.
Authored mixture/humidity applicability remains a source requirement.

Analytic regressions cover compressible and layered rates, anisotropy, gas
substitution, slip calibration, selective flux and reversal, common species
withdrawal, mass/enthalpy conservation, equilibrium overshoot and cancellation.
These are constitutive and numerical checks, not validation for an unnamed
porous solid or polymer. The existing measured acoustic fixtures remain active.


## External species streams through the original gas energy owner

`thermochemical_gas::energy::parcel::transfer::species::boundary` accepts
integrated species mass and complete stream energy on the original gas
thermochemical reference. `ClosedGasParcel::exchange_species_boundary` reuses
the existing species-inventory, composition, EOS, transport and relative-energy
inversion paths to return one immutable endpoint with extensive balances.
Carrier inventories, volume and external query coordinates remain intact.
Exact individual-species depletion and reentry use a retained complete source
catalog; no trace gas or replacement caloric law is inserted. Active caloric,
transport, source-pack and water-identity substitutions refuse. Candidates
borrow the immutable catalog to avoid copying source cards on each trial.

Original-owner queries expose species mass, partial pressure and partial
specific enthalpy. Opposing streams may carry energy with zero net species
mass. The external phase must debit the same signed energy on a compatible
reference; latent heat alone does not determine NASA-9 vapor enthalpy. This
primitive supplies neither an interface rate nor automatic thermal-reference
alignment. Six focused native regressions passed for influx/deposition,
depletion/reentry, reference shifts, counterflow, ordering and atomic refusal.

### Thermal entropy adapters

The shared ThermalLaw adapter for EquilibriumEnthalpyPhaseCurve delegates to
the original exact piecewise dh/T integral, including latent plateaus once.
GasEnthalpyLaw supplies the same optional entropy interface by reconstructing
both enthalpy endpoints with its original NASA/EOS law and evaluating the
cancellation-safe frozen-mixture entropy change at the retained pressure.
Source domains, bounded inversion and frozen-composition restrictions remain
in force. These adapters enable shared thermal entropy audits without a second
caloric law. Independent sensible/latent and variable-cp gas regressions cover
storage, endpoint reversal, no-ops and refused source-domain exits.

`GasEnthalpyLaw::advected_entropy_anchor` also supplies the existing stationary
species transport owner with a source-bound entropy reference. At the common
lower caloric temperature it retains original enthalpy, ideal mixing entropy,
and the source pressure-entropy increment. Pure-species standard entropies at
that temperature and standard pressure are cancelled analytically by stationary
component conservation, before a large source datum can erase mixing. The
reference binds the full source catalog, caloric policy and anchor temperature;
different pressure states remain explicit entropy contributions. The same
original Cp/T integrals extend each frozen composition away from the anchor.
Native source/transport regressions cover equal- and unequal-temperature
mixing, pressure changes, incompatible anchor temperatures and large positive
and negative species entropy-reference shifts. This closure is for stationary
nonreacting species; it does not supply reaction affinity or phase selection.

The same original gas source can bind a selective stream through thermal
`SelectiveOutletLaw`. Its declared composition, enthalpy reference and actual
outgoing enthalpy must agree with that source law. The stream uses its own
temperature and entropy, including pure-species outlets, and its boundary
enthalpy and entropy are removed once. Signed separation contributions remain
visible; an isothermal separator's entropy deficit is not supplied with
invented work or reservoir heat. Native NASA regressions cover separate stream
temperatures, an independent energy/entropy balance, incompatible composition,
cancellation and large positive and negative species entropy-reference shifts.


## Molecular sublimation with evolving gas energy

`sublimation::gas::advance_surface_with_gas` couples one original solid caloric
law to the original `ClosedGasParcel` and its retained complete source catalog.
Both temperatures, species inventories, pressure and gas properties evolve.
The explicit species binding supplies a fixed solid enthalpy reference offset;
source molar mass must agree, and h_v(T)-(h_s+offset)=L(T) must hold within the
caller's declared absolute tolerance at both accepted surface endpoints.
Incompatible caloric temperature dependence refuses rather than being hidden
by a new energy datum.

The shared molecular law supplies gross incoming and outgoing streams. Each
stream carries the original species h(T)-R_specific*T/2, with separate surface
and gas temperatures. Counterflow can therefore transfer heat at zero net mass.
The original gas boundary owner reconstructs the endpoint once. The candidate
checks actual combined mass and energy and retains signed recession; neither
input changes until the surrounding consumer accepts both endpoints together.

This is a bounded first-order explicit interval. Both gross withdrawals must
fit their initial inventories and the caller's maximum fraction, including
zero net mass with unequal temperatures; timestep refinement resolves thermal
and pressure feedback. Exact equal-temperature saturation is a stationary
turnover exception: balanced molecules may exchange repeatedly without changing
either inventory or energy. Gross stream histories and molecular pressure loads
remain active. Equality is checked on the original temperature and pressure,
before exposure multiplication; no tolerance-based equilibrium snapping occurs.
Unsupported domains, incompatible references, unrepresentable increments,
exhaustion and cancellation refuse atomically. The local Maxwellian/sticking
assumptions exclude a gas film, carrier thermal accommodation, external heat,
resolved flow and geometry changes. Multi-surface shared-chamber evolution is
separate and retains its thermostatted contract.

Six new native coupled regressions passed with the original thermochemistry,
EOS, transport and gas-energy implementation: dry-carrier sublimation,
deposition, zero-net-mass heat transfer, formation-reference invariance,
timestep convergence and atomic refusal. The six existing species-boundary
regressions also passed in the same actual-source native executable.


### Explicit equilibrium sources for tensor material memory

`visco_tensor::source::TensorElasticity` admits explicitly selected
Young/Poisson, Young/shear and bulk/shear equilibrium elasticity through the
shared isotropic material projection. Original source properties remain in the
complete resolved bundle; the native tensor Prony owner still requires every
declared branch coefficient and relaxation time. Alternative equilibrium pairs
therefore do not invent a bulk/shear relaxation spectrum. Original model tags
and canonical identities are retained, and selected source roles distinguish
the new physical laws.

The solid material parser and full-tensor plate/string consumers carry these
choices into the existing history evolution, stress, stored energy and physical
loss. Scalar Young-based spectra retain their original source assumptions.
The affine thermal tensor adapter admits independent affine bulk/shear sources
through exact stiffness derivatives, including a changing derived Poisson ratio.
The shared `TensorThermalLaw` also admits independent affine Young/shear
sources through the exact rational bulk response and its thermal curvature;
the strict affine adapter retains its common-scale restriction. Added native regressions cover source retention,
stress/work/heat, actual plate/string pressure equivalence, and the existing
projected thermal owner's free energy, entropy and internal energy. Native
workspace execution remains separate from source and syntax checks.

### Source-derived two-surface thermal Hertz compliance

`elastic::thermal::contact::ThermalHertzPair` keeps both original isotropic solid cards, complete query contexts and exact source selections. Each body has its own actual temperature. The selected Young/Poisson, Young/shear or bulk/shear coefficients must follow continuously supported affine source paths; moduli remain positive, Poisson ratios remain admissible, and yield limits remain fixed. Unsupported derivative jumps, hidden evolving axes and extrapolation refuse. Semantic absolute-temperature paths reuse the shared material admission helper, while existing dimensional-kelvin sources retain their compatibility path. Refresh preserves every original requirement and selected claim, including auxiliary properties. The original `IsotropicSolidStatePoint::reduced_modulus_against` computes harmonic compliance; its full two-temperature gradient and Hessian are analytical even when the derived Young modulus and Poisson ratio are nonlinear. The original fixed-Poisson calculation and identity remain unchanged for the canonical source model.

For a fixed-geometry native Hertz potential `psi=G(delta)*E*`, the derivatives support surface entropy `s_i=-G*dE*/dT_i`, internal energy `U=psi+sum(T_i*s_i)` and reciprocal thermal capacity. Returned source states retain their original identities and receipts. This is a constitutive capability, not a contact or thermal integrator; coupled owners must retain physical heat/work, total heat capacity, source domains and mechanical state. No geometry expansion, changing reference mass, phase change, plasticity, measured calibration or full Euler thermal closure is inferred.

Focused native G1/G3 tests compare the modulus and all first/second derivatives with independent rational compliance and finite differences, check thermodynamic internal-energy derivatives, ordered surface exchange, the constant-modulus limit, and original-card/domain/non-affine-source refusal. Additional regression sources exercise explicit independent elastic pairs, changing Poisson ratio, preserved auxiliary selections, source-domain refusal, and derivative cancellation. The general compliance Hessian need not be negative semidefinite; the coupled owner must admit the complete caloric plus contact capacity. Synthetic source fixtures establish implementation behavior only.

### Physical cohesive viscosity

`BilinearCohesive::with_viscosity` selects explicit Kelvin--Voigt traction
`(1-d) eta_i gdot_i` with independent normal/shear interface viscosities in Pa s/m.
`trial_over_time` requires positive physical time and a declared jump-rate limit;
untimed trials refuse while bonds survive. After accepted complete failure no
dashpot remains, so untimed continuation is also admitted. The original common-mode or frozen-onset
BK fracture owner still owns damage and storage. Normal viscosity acts in both
opening and compression while bonds survive, and all viscous traction vanishes
after complete failure; the separate compression penalty remains.

Endpoint total traction and its consistent, possibly nonsymmetric tangent are
separate from the original equilibrium account. Exact piecewise analytic
constant-rate path integration reports nonnegative viscous loss independently
from intrinsic fracture. Energy uses that integral, not endpoint traction times
jump. The short-interval binomial evaluation has below-roundoff truncation in
its admitted range. Longer ramps evaluate the logarithmic integral directly,
preserving nonzero loss at large failure/onset ratios without cancellation of
the leading terms. G0/G3 regressions compare full-rupture viscous work with its
independent closed form and retain the same work under reversal-path subdivision.
`trial_over_time_with_mean_traction` adds the same exactly integrated surviving
bond fraction times `eta_i gdot_i` to the equilibrium vector mean. Its three-mode
work equals intrinsic work plus the separately reported physical viscous loss;
it reuses the timed owner's integral rather than a second damage evolution.
No stored dashpot state, healing, fracture heat partition,
evolving-temperature refresh or universal adhesive calibration is inferred.

The source adapter queries intrinsic coefficients and explicit interface
viscosities in the same ordered original bundle, retains receipts, and rejects
bulk viscosity units and unowned evolving axes. Intrinsic fracture work must
exclude viscous losses to avoid double counting measured finite-rate work.


### Heterogeneous calibrated slip in layered gas barriers

`porous::permeability::DarcyGasBarrier` admits a different calibrated Klinkenberg
factor in each layer. One bounded scalar solve conserves mass through the
original ordered layers; factors are never averaged into one apparent material.
The exact common-factor fast path remains constant-time and ordinary flow
queries remain allocation-free. The shared support intersection and actual gas
model/composition admission still apply to every layer. The continuum linear
resistance adapter refuses any nonzero slip.

`mass_flow_profile[_for_model]` additionally returns the absolute layer-interface
pressures, iteration count and normalized pressure-closure residual. Sealed and
zero-area paths do not invent interior pressures. Reversing both the physical
stack and its reservoirs reverses the flow; reversing only the reservoirs can
change its magnitude when the layer factors differ. Equilibrium conductance
uses the exact small-signal series resistance. This is rigid, isothermal,
constant-composition slip flow, not a Knudsen or thermal-transpiration model.

Six new native regressions and all 19 existing porous/acoustic regressions pass,
including independent per-layer mass balance, manufactured interfaces, stack
orientation, subdivision, adjacent-float pressure differences, exact closed
limits, calibration refusal, extreme slip contrast and bounded 64-layer work.


## Joint nonisothermal multi-surface sublimation

`sublimation::gas::chamber::advance_chamber_with_gas` advances heterogeneous
surfaces against one original gas snapshot. Each surface retains its explicit
identity, species, phase law, caloric callback and fixed energy-reference offset.
Every participating species is evaluated once at the original gas state.
Canonical surface/species accumulation precedes exactly one original gas
boundary reconstruction; input ordering cannot change the accepted result.

Gross incoming streams share each species' OLD inventory allowance. Gas emitted
by another surface cannot fund an additional withdrawal in the same explicit
interval. Each solid separately bounds gross emission. `InventoryLimit` returns
a conservative retry duration instead of clipping streams. The stationary
turnover exception applies only when EVERY surface is exactly at its gas
temperature and species saturation pressure; one evolving neighbor retains all
original gross limits. Original source compatibility, reflected heat, phase
inventory, momentum and cancellation checks remain active. Common-reference
energy aggregation preserves heat transfer during counterflow without repeatedly
adding large opposing formation energies. All solids and gas pass actual total
and per-species mass balance, original caloric compatibility and energy checks
before one immutable joint candidate is returned.

The existing single-surface API shares the same physical preparation. Fixed
area, local Maxwellian vapor, original source domains, specular nonsticking
reflection and pure-solid restrictions remain. This explicit model does not
supply gas transport, mixed phases, extinction topology or external heating.
Six new native regressions cover warm/cold redistribution, different volatile
species, exact single-surface equivalence, order independence, shared gross
withdrawal, energy-reference invariance and atomic refusal. All 18 actual-source
joint/scalar/species-boundary tests passed with the original gas reconstruction.


### Reciprocal molecular loads on supported chamber surfaces

`sublimation::gas::momentum` composes force and torque impulses with the ORIGINAL
joint chamber candidate, optional sharp fusion charts and explicitly declared
Maxwell reflection. An internal observation seam consumes the same accepted
gross emission/sticking masses and original temperatures without replaying
constitutive sources. Retained gas partial pressures are cached once; canonical
compensated sums retain carrier pressure and opposing surface wrenches. The
declared unit normals, exposed-area centroids and common moment origin are
checked against the exact selected surface roster.

Emission, sticking and incident-plus-reflected nonsticking molecules contribute
separately. At equal-temperature equilibrium or a sealed surface, zero net
phase transfer still gives the complete pressure load. Carrier species exert
pressure even without phase exchange. `MaxwellReflection` explicitly declares
the diffuse fractions and shares those exact fractions with the original
reflected-heat owner; arbitrary energy accommodation is not reinterpreted as
an angular scattering law. Phase, fusion, heat and mechanical accounts escape
only after the original chamber and final cancellation gates succeed.

These are stationary supported patches in the zero-drift local Maxwellian gas
frame. Support work is zero; reciprocal gas impulses are boundary transfers,
not endpoint bulk momentum or kinetic energy in `ClosedGasParcel`. The caller
must not add ambient pressure or molecular kinetic heat a second time. Moving
surfaces, spatial gas transport, adaptive load histories and rotational guide
dynamics require their own owners. Six integration regressions cover original
candidate equivalence and source-call counts, equilibrium/carrier loads,
Maxwell heat consistency, frame/origin covariance, sharp fusion and atomic
zero-exposure/late-refusal behavior.


### Time-resolved molecular impulse on accepted adaptive histories

`sublimation::gas::momentum::history` uses the ORIGINAL chamber step-doubling
controller for both pure-solid and sharp-fusion trajectories, with optional
explicit Maxwell reflection. Each interval compares the coarse mechanical
impulse with the sum of the two fine impulses, in addition to all original
mass, temperature and constituent-inventory checks. Independent SI tolerances
apply to each surface's emitted, sticking, reflected and total normal impulse,
its moment and the complete force/moment resultants. Large carrier pressure or
opposing faces cannot hide an inaccurate recoil contribution.

Only the accepted fine intervals enter compensated history sums. Rejected and
coarse evaluations consume the original work budget but contribute no physical
mass, heat or impulse. Generic private hooks share one original controller;
existing public adaptive APIs retain zero-sized data and no-op hooks. The
complete original endpoint, stream/conservation audit, optional whole-window
fusion account and integrated loads publish together after final admission.
Storage remains O(surface count), independent of accepted-step count.

Frames, exposed solid areas and the moment origin remain fixed. This resolves
supported-surface loads without adding bulk-gas momentum or moving mechanics.
The estimator controls local temporal error, not a rigorous global quadrature
bound; original mass/species/energy budgets still cover the complete window.
Four regression groups cover source-call/original-path parity, independent
analytic thermal-trajectory quadrature with load-driven refinement, exclusion
of rejected-trial pressure, and sharp-fusion continuation plus atomic work,
zero-duration and cancellation behavior.


### Adaptive closed-chamber histories

`sublimation::gas::chamber::advance::advance_chamber_adaptive` compares one joint
explicit interval with two half intervals and accepts the finer conservative
history without extrapolation. Local absolute/relative tolerances cover every
solid mass and temperature, gas temperature and EACH retained species mass,
including trace vapor. Stable surface ordering makes states and retry decisions
independent of caller ordering on the same ISA. Typed inventory-limit refusals
reduce the whole interval; invalid source/caloric data and cancellation remain
immediate failures. Trial, minimum-step and maximum-step controls bound work.

The candidate retains only accepted endpoints and cumulative molecular streams.
Actual full-history total and per-species mass, combined energy, and each side's
stream/storage energy account must meet the original absolute balance budgets;
those tolerances are not multiplied by step count. Gas/solid inputs remain
unchanged on a late refusal. Areas and source bindings stay fixed throughout
this driver; thermal events use their separate owner. Step doubling supplies a
local time-error estimate, not a rigorous global accuracy or physical-validation
certificate. Six focused adaptive tests passed, including convergence toward a
2048-step reference, shared-inventory retries, trace-species refinement, exact
order invariance, work limits and late cancellation. The combined native
adaptive/joint/scalar/boundary suite passed 24 tests with no failures.


### Coupled molecular, Knudsen and viscous transport in porous gas paths

`thermochemical_gas::energy::parcel::transfer::species::diffusion::dusty_gas`
implements the isothermal ideal-gas dusty-gas equations in the stationary-solid
frame. One homogeneous dry pore path couples all selected binary molecular
frictions, explicit Knudsen wall frictions and the original no-slip Darcy owner.
There is no zero-molar or zero-mass-flux constraint: unequal counterdiffusion and
uphill individual-species flux can follow from the coupled molecular friction.
A separate Darcy/nozzle/diffusion edge must not count that same pore area again.

The immutable medium retains the original intrinsic-permeability barrier and
explicit pore radius, gas-filled porosity, transport tortuosity and provenance.
Effective molecular diffusivity is porosity/tortuosity times the source value;
Knudsen diffusivity uses the actual thermochemical species molar mass. The
transport tortuosity is the declared coefficient in that equation, not an
inferred squared path ratio or the acoustic high-frequency tortuosity. Added
Klinkenberg slip and heterogeneous pore layers refuse to prevent overlapping
rarefaction models or an invented unresolved interface composition.

The face composition is concentration weighted at the declared mean temperature.
Original selected binary source laws are sampled at the face and admitted at
both endpoints, including pressure-dependent diffusivity. The shared material
owner supplies mixture viscosity and thermochemistry. A symmetric positive
friction system uses the existing dense Cholesky owner, retains every species
row and rechecks the original equations after total flux assembly. Work is
bounded to 32 species; loops poll cancellation. Reports retain source states,
individual Knudsen/diffusive/viscous rates and nonnegative local dissipation.
Nonrepresentable nonzero forces refuse instead of becoming false equilibrium.

`network::transport::GasTransportPaths::dusty_gas` samples every face from the
unchanged network snapshot. Its signed species transfers share one gross
withdrawal allowance, original donor partial enthalpies, material reconstruction
and final conserved-species entropy admission with every other mechanism.
Counterflow does not replenish an outgoing allowance. Local dissipation is a
diagnostic and is never added again as gas heat. Zero area queries no inactive
sources and returns no invented face state; zero intrinsic Darcy permeability
suppresses only viscous flow when explicit conducting pore diffusivities remain.

This is a bounded central face approximation, requiring spatial and physical
time refinement. It does not resolve thermal transpiration, Soret effects,
adsorption/surface diffusion, water saturation, changing pore geometry,
pore storage, pore-size distributions, inertial Forchheimer flow or a nonlinear
steady composition profile through a finite heterogeneous laminate. Source and
geometry applicability must be supplied for the actual material.

Before an executor restart, a native executable of the actual owners passed 89
distinct focused tests, including nine new dusty-gas regressions: pure-gas
Knudsen plus compressible Darcy limits, Graham counterdiffusion, an independent
ternary uphill-flux reference, concentration weighting, source/units/closed-path
refusals, absent-side species, finite mass/energy conservation, shared nozzle
withdrawal and constitutive admission. An additional extreme-range underflow
regression and the merge onto the concurrently improved entropy API were
source-reviewed after recovery; full native rerun and workspace validation
remain unverified because the restart removed the dependency caches.


### Runnable evolving-gas ice cold trap

`cargo run -p fs-material --example nonisothermal_ice_cold_trap` advances two
10-micrometer ice films at 260/252 K through one 100 mL water-vapor chamber,
initially 256 K and 130 Pa, for 2 ms. The original NASA-9 H2O coefficients and
IAPWS ice pressure/latent law remain the thermodynamic owners. The example
explicitly derives h_s(T)=h_v(T)-L(T)-datum on 250--270 K; this compatible reduced
ice chart is not independently measured ice calorimetry. Original dilute-vapor
transport remains within its stated sub-triple continuation. Fixed film area,
917 kg/m3 density, accommodation one and a well-mixed local vapor are explicit
scenario assumptions; no chamber-design or spatial-flow validation is claimed.

The actual native example passed both default and tenfold tighter local-error
runs. Warm loss/cold deposition were 26.567540/23.881098 micrograms by default,
with gas temperature 255.821736 K; tighter results were 26.564590/23.878148
micrograms and 255.821858 K. Both runs assert whole-history original-datum mass
and energy conservation and report gross counterflux, adaptive work and residuals.
Run with argument `0.1` to tighten local temporal error targets tenfold.


### Evolving gas on the original multi-cell heat timeline

`sublimation::gas::timeline::advance_deposition_with_gas` composes original
thermal transport with simultaneous molecular exchange against one evolving
gas. Every surface retains its original cell identity, phase/caloric source,
species binding and fixed datum offset. Original heat observers check each
selected pure-solid enthalpy branch after EVERY accepted thermal substep and
birth. One joint phase interval then updates the actual masses and enthalpies
through the thermal owner's transactional boundary.

Start-of-interval exposure and exact birth/contact times remain authoritative.
The original thermal work ceiling applies to the whole requested advance.
Empty pre-birth periods jump to their next event without renewing empty work.
Actual combined gas/solid storage is compared with incoming births and external
heat, using stable initial/birth-relative energy changes and explicit datum
conversions. Late source/domain/balance failure or cancellation discards the
complete candidate, including elapsed time and consumed events.

This fixed-area heat/phase composition is first order; refine the thermal
max_dt_s to resolve splitting and explicit molecular exchange. The standalone
adaptive chamber driver does not certify this operator-splitting error. The
caller still owns pure-solid admission, supported source domains and topology;
melting, extinction, moving geometry and gas flow are not inferred. Ten new
actual-source native integration tests passed for physical heat/gas feedback,
order, births/datum shifts, exposure events, work limits, solid-branch refusal,
conservation and rollback.

### Deposition through sharp freezing and melting on the heat timeline

`sublimation::gas::timeline::advance_deposition_with_equilibrium` and its
complete-split adaptive counterpart admit explicitly selected
`GasTimelineEquilibriumSurface` bindings. Each retains the original
`SharpMeltingCaloric`; its identity, including latent energy and datum, must
match the registered thermal law. Temperature-only agreement on an isothermal
plateau cannot substitute for that identity. The original heat observer checks
the chart and a positive solid inventory after every thermal substep and birth.

The existing partitioned chamber supplies all molecular exchange, shared old-gas
withdrawal limits, optional reflected heat and original gas reconstruction.
Thermal conduction and authored sources can cross the sharp fusion plateau;
deposited vapor may melt some of the solid or concurrent heat extraction may
freeze liquid. Bulk enthalpy owns fusion energy exactly once. Whole-window
`condensed_phases` reports compare original or incoming-birth inventories with
the accepted final inventories, including conversions during heat substeps.
They are not the sum of molecular-only phase reports.

The thermal boundary commits a source-identity-bound `CondensedCellExchange`
through the same stale-cell, positive-mass, enthalpy-domain, reference-area,
work and atomic-cancellation checks as pure-solid exchange. Existing pure-solid
entry points retain their branch restrictions. The adaptive owner controls
complete heat/molecular splitting, including bulk enthalpy error at constant
temperature, and retains one work ceiling and whole-window conservation budget.

Every exposed face is explicitly admitted SOLID area; bulk phase fraction does
not determine coverage. Complete melting at any intermediate heat or molecular
endpoint refuses the entire requested advance. Nucleation on an all-liquid or
bare surface, alloy mushy ranges, changing geometry, vapor pressure work and
experimental phase-diagram validation remain outside this reduction. Focused
regressions cover analytic deposition/cooling/fusion balance, source-datum
invariance, heat-contact freezing, fusion-onset crossing, source mismatch,
complete melting, births, late cancellation and adaptive temporal refinement.

### Conservative melting completion at an exposed solid interface

`sublimation::gas::equilibrium::boundary::advance_to_melting` returns the original
fully liquid fusion endpoint when molecular exchange and optional reflected gas
heat complete melting. It reports consumed and unconsumed physical time, the
same condensed/gas endpoint and its actual phase conversion. The event equation
reuses the original source-relative extensive storage increment, including both
gross molecular streams, the original latent alignment and reflected heat.
Accepted-duration streams and heat pass the original initial-donor gates before
one shared gas reconstruction. An exact source-knot endpoint must fit one
sum-of-absolute condensed/gas energy allowance; no latent inventory is discarded.
Deposited mass may also melt, so the converted solid mass can exceed the initial
solid inventory without treating gas deposition as fusion.

Ordinary scalar and timeline entry points retain their positive-solid admission.
The completion API supplies an explicit handoff state; it does not select liquid
coverage or replace molecular transport with continuum film transport. A caller
continuing through a liquid owner must retain the original bulk chart, species
and caloric-datum offset. Any external heat split must be retried as a complete
interval at an earlier event. Frozen-rate event timing still requires physical
timestep refinement and does not establish moving-interface geometry.

`sublimation::gas::timeline::equilibrium::boundary::advance_deposition_to_melting`
extends that completion owner to one selected solid interface on the original
native heat timeline. Network heat, both molecular streams, optional reflected
gas heat and every neighboring cell share the accepted event duration. A clipped
molecular event only brackets time; the complete native/phase trial is replayed
from the unchanged accepted start. The first fully liquid native substep remains
visible even if later heat would freeze it again. Discarded trials spend the
original thermal work allowance and a bounded 64-trial event search. Initial
and scheduled births retain strict original source identity and solid inventory.

Native heat can finish melting under full coverage or with exactly zero net
molecular mass and combined molecular/reflected energy. The latter is a separate
identity phase operation commuting with native heat at the fusion endpoint:
it retains the real whole-interval old solid and vapor donors, original gross
turnover and source admission, without inventing a positive-solid post-heat seed.
Exact liquid-knot selection spends the original thermal joule allowance and is
recorded in stored energy and residual, never as an added heat source. Active
molecular completion uses the existing scalar event and its one gas reconstruction;
the timeline only converts its accepted receipts into the existing joint accounts.
Whole-window `condensed_phases` includes the final zero solid inventory, and the
receipt returns unconsumed time for an explicitly chosen liquid-interface owner.
Existing multi-surface, observed and adaptive entry points retain their strict
admission. Source-domain escape and unsupported intermediates still refuse;
physical timestep refinement and original donor limits remain necessary.

### State-dependent heat transport during deposition

`sublimation::gas::timeline::observed` supplies pure-solid and sharp-equilibrium
entry points accepting an original `DepositionObserver`. The same thermal
kernel owns conduction, radiation, interval sources, endpoint sources,
conductance overrides and timestep limits. Observer transitions also clip the
outer molecular interval. The observer begins once per requested advance;
internal heat/phase intervals continue its accepted history without resetting
source-energy accumulators. The returned observer, thermal state and finite gas
are one transactional candidate, including on late source failure.

External bath and source heat enters the combined balance exactly once through
the native thermal audit. Molecular and reflected gas heat remain internal
transfers. Heat callbacks see each thermal endpoint before molecular exchange;
the next interval sees accepted deposited or sublimated mass. The existing
`ConductionObserver::with_planar_phase_layers` can therefore update normal
thermal resistance from mass at a declared constant density and fixed
footprint. That geometry must be physically admissible for the selected phase
regime; mixed solid/liquid admission does not infer constant density.

All selected exposed areas are evaluated in one original canonical coverage
sweep per molecular interval, preserving subtraction order while avoiding a
full contact scan for every surface. This observed path is a fixed-step,
first-order heat/phase composition and requires temporal refinement; it does
not by itself provide an adaptive observer error estimate. Focused regressions cover
unchanged no-op trajectories, deposition/resistance feedback, source transition
clipping and accumulated heat, conservation and transactional source refusal.

`observed::advance_deposition_adaptive_with_gas` and
`observed::advance_deposition_adaptive_with_equilibrium` connect these same
observers to the original complete heat/phase step-doubling controller. Each
trial starts from the same committed observer; only the two accepted fine
intervals contribute physical heat and source history. `begin` runs once for
the whole request. Source transitions clip the coupled trial, including
interval-dependent observers before the first thermal cell is born. Coarse
and rejected thermal work spends the same original substep ceiling. Typed
source failures propagate without being retried as numerical error.

The local estimate covers every cell's mass, temperature and caloric increment,
plus original gas temperature and species inventories. Enthalpy error remains
visible on a sharp fusion plateau. It does not independently bound private
observer states, global time error, spatial resolution or empirical model error.
The returned observer, thermal state and gas remain one transactional candidate.
G0/G1/G3/G4 regressions exercise no-op compatibility, accepted-only pulse heat,
source transitions, growing-layer conduction refinement, equilibrium phase
inventories and rollback after late source/cancellation refusal.

`cargo run -p fs-material --example nonisothermal_ice_cold_trap --
--conducting-layers` runs a finite-vapor, two-temperature ice-layer scenario
using the original IAPWS ice and NASA gas owners, with mass-dependent conduction
to prescribed baths. It prints layer growth/retreat, temperature, gas pressure,
mass/energy closure and three temporal resolutions. Its specified dense-layer
density and conductivity are scenario inputs; it does not infer porous frost,
spatial gas transport, bare-surface nucleation or an experimentally calibrated
cold trap.
Add `--adaptive` after `--conducting-layers` to control the complete split with
explicit kg/K/J-per-kg tolerances and report accepted/rejected intervals and
executed thermal work. An optional following `0.1` tightens those local error
targets tenfold while retaining the same original physical and conservation laws.

## Sourced sharp solid/liquid melting

`phase::melting::ResolvedMeltingCurve` joins two continuously supported source
Cp/density branches and independently sourced melting temperature/fusion
enthalpy into the original `EquilibriumEnthalpyPhaseCurve`. The caller explicitly
declares one sharp equilibrium transition. Solid support must end, liquid
support must start, and the transition query must lie exactly at the sourced
melting temperature. Every query carries an explicitly typed absolute
temperature. Pressure, composition and other external context, including quantity
descriptors, must be identical throughout; an optional typed categorical source
selector is zero for solid, one for liquid, and absent from the transition query.
The constructor refuses missing density/Cp/fusion data, disconnected or discrete
support, incompatible context, nonpositive fusion energy, unrepresentable energy
increments and combined charts exceeding the original generated-knot limit.

The existing heat-capacity owner integrates sensible heat on both branches. One
fusion increment separates the solid/liquid endpoint enthalpies at the exact
same temperature; independently tabulated assigned enthalpy is not also added.
The resulting chart, its states and relative entropy remain consumable by the
original enthalpy owners. Density retains the original piecewise-linear rho(T)
on each sensible branch and the additive-specific-volume coexistence rule;
specimen mass stays invariant and the consumer must account for the resulting
volume change. Every selected Cp,
density, melting-temperature and fusion-energy receipt is retained and bound
into the curve identity. The declared interpolation allowance independently
bounds sensible interpolation and the difference between the represented fusion
plateau and its original sourced latent heat. Large enthalpy datums that distort
the plateau beyond that allowance refuse; the sourced value and receipts are
never silently replaced. This is not a combined chart-error bound, a bound on
sensible-integration roundoff, or a source-uncertainty claim. This admission does
not infer a phase diagram, nucleation, kinetics, supercooling, hysteresis, moving
geometry or a fixed-volume transport closure.


### Prescribed wet-pore gas response

`porous::wet` supplies explicitly calibrated gas relative permeability and pore
diffusivity at a prescribed liquid saturation, total porosity and gas-pore radius.
The nonwetting van Genuchten--Mualem relation uses supplied residual liquid/gas
saturations and shape parameter. Diffusion independently selects either the
Millington--Quirk 1961 relation or an explicit transport tortuosity. These are
separate mobilities: residual Darcy immobility does not invent a diffusion
connectivity threshold. Full liquid saturation returns exact zeros and no
fabricated finite gas tortuosity.

Immutable samples retain the complete calibration through shared storage.
Sampling uses bounded deterministic scalar work with no sample allocation,
retains positive transport near closure, and refuses unrepresentable positive
coefficients. Domain bounds apply before the declared endpoint extensions.
Porosity and radius evolution are prescribed; intrinsic permeability, liquid
mass, capillary pressure, mechanical work and hysteresis are not inferred.
Nine focused native tests passed with the original MaterialError and fs_math,
including independent analytic limits, exact cube-root diffusion fixtures,
geometry scaling, monotonicity, source-domain admission and numerical limits.


## Lower standard atmosphere for flight consumers

`gas::atmosphere::StandardAtmosphere1976` derives the seven hydrostatic layer
bases from the U.S. Standard Atmosphere 1976 surface state and table 4 lapse
rates. Geometric and geopotential altitude have separate query methods; the
source's effective-radius mapping and gravity agree with the pressure
integration. Each accepted sample uses the original `GasState` dry-air owner
for density, sound speed and transport. Layer-base pressure is derived once
instead of using independently rounded values. Successful construction and
queries use fixed arrays without heap allocation. Near layer boundaries,
`ln_1p` retains small increments. Repeatability is same-platform/toolchain;
no cross-ISA bit-equality claim is made for the host logarithm.

The admitted branch is 0..80000 m geometric altitude, where the source retains
constant molecular mass. Above that boundary, USSA table 8 requires additional
composition variation; this implementation refuses extrapolation. The profile
is a dry standard reference, not measured weather, wind, humidity, terrain or
a rotating/inertial Earth model. Its effective hydrostatic gravity must not be
silently substituted for a different flight-frame gravity law. Focused tests
cover source layer values, the common gas state, layer continuity, the altitude
crosswalk, hydrostatic balance and upper-domain refusal.


### Saturation in the shared gas transport owners

`PorousGasLayer::with_wet_state` retains intrinsic permeability, original axes,
source records and gas support while applying the sampled gas relative mobility.
`DarcyGasBarrier::with_wet_states` recompiles each original series layer once;
repeated saturation updates replace their preceding factors. A flooded layer
closes the pressure-driven path exactly. Existing Darcy consumers receive this
behavior through the same barrier. Dry-gas Klinkenberg calibrations refuse wet
use rather than silently extending gas-specific slip authority.

`DustyGasMedium` distinguishes reference pore geometry from its effective wet
geometry. The molecular and Knudsen terms use the independently sampled absolute
diffusivity factor; the Knudsen term also uses the current gas-pore radius.
Relative Darcy mobility is not multiplied into diffusion. Full liquid saturation
closes all gas transport before inactive gas or binary sources are queried.
The complete wet calibration remains attached to the reported medium.

These are prescribed instantaneous constitutive updates. They do not move liquid,
resize a retained gas inventory, account for displaced-gas work, infer pore
connectivity or evolve a hysteresis branch. Those operations require a coupled
storage/geometry owner. The active gas law retains its original domain and
finite-step mass, energy and entropy admission.

Eight added integration regressions cover independently derived pure/binary
flux scaling, residual-gas diffusion, exact flooding, repeated wet updates,
radius dependence, retained source identity/slip refusal and anisotropic series
layers. Both modified production files and all new tests passed Rustfmt parsing.
The nine standalone wet-law tests passed natively. The complete selected-source
Rust test typecheck also passed after restoring the real gas dependencies.
Native execution of the eight added integration tests remains blocked by full
executor storage; it is not claimed here.


### Transient heterogeneous pore-body gas storage

`network::transport::pore_domain` retains finite reservoirs and an ordered set
of actual gas inventories. Each pore body's gas volume is its authored gas-filled
porosity times bulk area times thickness. Initialization uses the original EOS;
provided inventories must already match that volume. Independently authored
homogeneous throats retain their own length, area, material and dusty-gas closure.
No averaged heterogeneous interface law or continuous laminate profile is inferred.

Each accepted substep samples every throat from one unchanged gas snapshot and
uses the existing joint species/enthalpy reconstruction. Incoming gas cannot fund
another withdrawal within that snapshot. Typed inventory and caloric-solve
refusals trigger bounded timestep halving; invalid sources and domains propagate.
A failed late substep or cancellation discards the complete requested window.
Whole-window mass, individual species, energy and entropy are rechecked from the
actual initial and final inventories using the original absolute budgets.

Reports aggregate directional species and donor-enthalpy throughput, pore uptake,
work counts and whole-window balances. Only the last substep's detailed source
report is retained. The total omitted throat void volume is measured against
stored pore-body volume and must fit an explicitly declared geometric bound.
This first-order fixed-volume model needs temporal refinement. The adsorbing
variant below connects original finite-site and solid-heat owners to this same
window loop. Moving porosity and a resolved nonlinear throat profile remain
separate physical models. Explicit skeleton conduction uses the finite-adsorbent
extension below.

Five focused regressions cover delayed forwarding, first-order convergence
against the independent three-volume diffusion eigenmodes, heterogeneous
storage/throat coefficients, bounded adaptive withdrawal and atomic refusal.
New files passed Rustfmt parsing and the selected actual-source Rust test
typecheck completed successfully. Native execution remains blocked by full
executor storage. These five tests are not claimed as run.


### Finite gas/solid adsorption with binding-energy storage

`adsorption` binds one retained gas species to a homogeneous Langmuir site
population, an original capacity/affinity/kinetic material source and two original
single-phase caloric charts. Capacity is per dry skeleton mass, affinity has
inverse-pressure units, and the bound adsorbate has an explicit enthalpy-datum
offset into the original NASA gas reference. The reference affinity fixes its
entropy datum. Temperature-dependent affinity follows the same gas/bound caloric
integrals and satisfies the corresponding van't Hoff relation; adsorption heat
is their enthalpy difference and is not added again as an external heat source.

Frozen original partial pressure and surface temperature give an analytic site
update and both integrated gross streams. Gross adsorption is limited by the
original gas inventory; total site turnover is bounded separately. Newly adsorbed
molecules may desorb within this exact frozen site solution. Original gas bulk
enthalpy travels with each direction. One finite gas reconstruction and one
bracketed, safeguarded-Newton solid caloric inversion produce a joint immutable
candidate. Both temperatures evolve; the original source domains, thermal-change
bounds, final affinity, mass/energy budgets and combined gas/site entropy must
all admit the endpoint. Cancellation or any late refusal leaves inputs intact.

The source adapter rejects stale temperature/pressure aliases and unowned
conditions. Adsorbate caloric, transport and pack identities remain bound.
Stable relative caloric operations cancel arbitrary reference energies before
thermal inversion and entropy admission. The shared
`ResolvedHeatCapacityCurve::specific_entropy_change` now integrates the original
chart's chord Cp/T directly between temperatures, preserving tiny thermal changes
even when a large enthalpy datum would round their absolute h values together.

This rigid, near-isothermal local kinetic closure uses one noncompetitive site
population and a declared bulk-enthalpy stream convention. It does not supply
free-molecular impingement, surface diffusion, liquid displacement, swelling,
hysteresis or solid heat conduction. Time refinement is required for the finite
bath's pressure/temperature feedback. It does not impose a hidden thermostat.

Eleven adsorption regressions and one new shared caloric regression cover exact
frozen kinetics and turnover, independent affinity/entropy derivatives, finite
uptake/desorption energy and entropy, equilibrium turnover, reference invariance,
piecewise caloric inversion and atomic/source refusals. All pass actual-source
Rust test typechecking together with the wet/DGM/pore-domain tests. Native
execution of these new coupled tests remains blocked by full executor storage.


### Finite adsorption and skeleton heat along a transient pore domain

`network::transport::pore_domain::adsorption::AdsorbingPoreDomain` attaches one
original finite competitive site population to each selected pore body. The
existing bounded `PoreDomain` time-window loop now also advances the original
joint gas-transport/adsorption candidate. Each interval samples original
dusty-gas throats, site affinity and kinetics on one accepted gas snapshot;
gross withdrawals share the original inventory. Bound species and the finite
dry/bound solid heat evolve with the same source calorics, and every subsequent
interval uses the updated gas composition, pressure and temperatures.

Whole-window admission reconstructs final gas plus bound-species inventories.
Native reference-cancelled gas/solid energy defects and complete-system entropy
increments accumulate with compensation under a single window allowance,
without multiplying the allowance by substep count. Reports retain gross
directional throat throughput, gross adsorption/desorption, final original
source reports and finite solid energy changes. Typed withdrawal, site turnover
and caloric range/convergence refusals may trigger the same bounded halving;
source/model refusals do not. A late refusal publishes no gas, site or clock.

Rigid geometry and fixed gas-filled pore volume remain explicit. Adsorbed site
occupancy does not become liquid saturation or an invented permeability law.
No Darcy work, adsorption heat or formation-reference energy is counted twice;
swelling, capillarity and moving pore-volume work are not inferred. Skeleton
conduction requires the explicit original-source paths described below. The original native adsorption near-isothermal and frozen-forcing
limits remain active; successful bounded stepping is not a time-error estimate.

Four focused native regressions cover the independent two-species site/finite
heat formula, changed subsequent dusty-gas throughput with persistent loading
history, typed inventory-driven refinement and atomic late failure.


### Source-backed conduction between retained finite adsorbents

`adsorption::competitive::network::step_candidate_with_conduction` and the corresponding
`AdsorbingPoreDomain::step_candidate_with_conduction` / `advance_with_conduction`
consume explicit `AdsorbentConductionLink` paths. Each side binds the exact
retained dry-skeleton caloric object and its original `ResolvedConductionLayer`.
The existing integrated k(T) owner evaluates the physical path; optional
`AdsorbentThermalContact` retains the ordered original interface card and
selected resistance claims. Interface and caloric material identities must
agree in order, including the exact card hashes.

All links sample the same accepted gas/site/solid snapshot. Equal-and-opposite
conductive heat enters the original dry-plus-bound caloric inversion together
with sorption heat, exactly once. The native gas inventory, endpoint source,
complete-system entropy and single-window closed-energy gates remain active.
A heat step that exceeds physical admission refuses without clipping its flux
or altering any participant. The original no-link APIs use this same candidate
with an empty link list.

The existing bounded pore timeline preserves accepted gas, site and thermal
history. Later intervals use the changed temperatures in original adsorption
and actual dusty-gas transport. Whole-window reports retain gross heat in both
directions, and the final substep retains original layer/interface sources.
Geometry is authored explicitly; site occupancy does not infer saturation or
permeability. No separate thermal integrator, extra Darcy heat, moving geometry
or unconditional stability is introduced. The original first-order,
frozen-forcing limits and temporal refinement requirements remain.

Six manufactured regressions cover finite heat/contact arithmetic, data-only
conductivity substitution, exact retained caloric identity, endpoint source
refusal, total-entropy rejection of oversized heat, subsequent pore-flow
feedback and continuing-window atomicity.


## Variable-caloric off-design converging-diverging nozzle

`RocketNozzleRegime` distinguishes `UnchokedSubsonic`, `ChokedSupersonic`,
`ChokedExternalCompression` and `ChokedNormalShock`. Subsonic exit pressure at
the actual exit/throat area ratio determines unchoking; sonic throat pressure
alone is not the back-pressure threshold. A valid near-stagnation subsonic
solution does not require source support at a hypothetical sonic station.

`ClosedGasParcel::normal_shock_flow` treats its unchanged source parcel as the
upstream static state. Actual mass and momentum invariants determine the
velocity/temperature relation, and a bounded variable-caloric energy solve
factors out the trivial unchanged upstream root. The accepted shock reports
original source-resolved upstream/downstream states, actual entropy production
and conservation residuals. Internal nozzle shocks are located by section
area, followed by the same downstream isentrope and subsonic pressure recovery;
axial shock location requires a separately supplied nozzle contour.

`RocketNozzleFlow` retains throat Mach, regime, optional shock report and an
optional back-pressure residual for subsonic discharge. Its existing entropy
energy residual now subtracts the separately reported physical shock entropy
rise. Total axial thrust remains `mdot*u_exit+(p_exit-p_ambient)*A_exit` with
pressure thrust counted once, including at shock-on-exit transitions. Inner
caloric stations allocate no material receipts; accepted physical stations
retain the original complete property/transport evidence. Vacuum and low back
pressure skip the unnecessary subsonic-threshold solve. Existing blowdown and
vehicle consumers use the same entry point and original accepted impulse.

The model is quasi-steady, quasi-one-dimensional, inviscid, adiabatic frozen
ideal gas on an ordinary branch with one sonic crossing. External compression
is the attached ideal-flow branch, not a claim that a real nozzle will remain
attached. No viscous separation/side-load prediction, oblique plume solution,
shock train, transient shock tracking or reaction relaxation is supplied.
Equal/higher ambient pressure still refuses forward propulsion flow.

Native actual-source checks: 38 nozzle-family tests passed, including perfect-
gas shock oracles from Mach 1.0001 through 5, variable-cp conservation and
entropy, independent internal shock areas, unchoking thresholds, pressure-thrust
continuity, restricted caloric/transport domains, cancellation and actual finite
blowdown through the shock regime. The lower standard atmosphere's three native
tests also passed. This is a scoped native build, not a full workspace result.


## Equilibrium chamber to source-preserving frozen products

`ThermochemicalGasModel::equilibrium_chamber_hp` and
`equilibrium_chamber_hp_enthalpy` bind the original neutral ideal-gas
composition owner to the existing material parcel owner. The material bank
must cover every chemistry species and retain exactly equal complete NASA-9
models, including conventions, coefficient cards, provenance and caloric
corrections. Feed amounts are independent of the bank's catalog fractions.
Chamber pressure and all other material query axes come from the typed context.

A common inlet temperature or an explicit sum of original stream enthalpies
sets the HP target. Formation enthalpy enters once. The accepted composition
uses the existing bounded generated-simplex closure, reports its maximum
roundoff correction, and refuses eligible positive traces that cannot be
represented by the frozen material owner. `with_composition` preserves active
thermodynamic models, transport selections and source-pack identities.

The final source-bearing material state is resolved once. Equilibrium mass
and its actual density determine volume, and the original `ClosedGasParcel`
is assembled without a duplicate transport evaluation. Publication checks
actual parcel mass against both feed and equilibrium amounts and actual
parcel enthalpy against both the requested target and the equilibrium result.
Original chemistry receipts, mass/energy residuals and frozen parcel are
retained together; source mismatch, transport-domain refusal, unrepresentable
trace, failed balance or cancellation publishes nothing.

The resulting storage, acoustic, caloric and nozzle composition is explicitly
frozen. The HP closure is an initial thermodynamic chamber state, not a finite-
rate combustor or a continuously fed engine. No reaction heat is added after
freezing; equilibrium nozzle relaxation, kinetics, ions, condensed phases,
real-gas behavior and flame stability remain separate physical models.

Native actual-source verification passed all 152 material thermochemical-gas
family tests, including eight new chamber/pack integration tests, the existing
storage/mixing/transfer family and off-design nozzle tests. The new cases
include an independent exothermic HP root followed by the actual nozzle,
different inlet temperatures, complete/exact source-bank admission, structural
absence versus trace underflow, original admitted-pack decode and receipt
retention, final source refusal and cancellation. Independent equation/source
review found no conservation, reference or trace-freezing blocker.

`cargo run -p fs-material --example equilibrium_rocket -- 300 3000000 100000 2`
runs a manufactured variable-cp A/B source example through the actual HP
chamber and frozen nozzle. Its native execution produced 902.699280458 K chamber
temperature with 3.47e-18 kg mass residual and 3.64e-12 J enthalpy residual.
These values validate that authored integration fixture, not real propellant
or engine performance. The full repository Cargo workspace was not run.

### Prepared reacting pressure/enthalpy material endpoints

`energy::equilibrium::pressure::EquilibriumGasEnthalpyPath` binds the original
prepared reacting PH workspace to this same exact source bank and product
freezing boundary. Typed reference T,p and original feed elements establish
an immutable INITIAL EQUILIBRIUM enthalpy reference. Queries prescribe pressure
and extensive H-H_reference; the native original-source inverse includes
chemical redistribution and reacting Cp throughout its bounded temperature
search. Source admission covers the complete declared temperature bracket.

The original feed mass is preserved exactly in every material endpoint.
Actual chemical mass and density times accepted volume are independently
checked against the caller's absolute mass allowance. Material enthalpy
increments are independently reduced from the actual frozen amounts:
`sum[(n_final-n_reference)*h_reference + n_final*integral(Cp dT)]`.
Double-double amount differences, products and sums retain small sensible
increments beside large formation datums without removing genuine chemical
formation terms. Direct endpoint freezing defects and the actual material
enthalpy residual must satisfy a separate absolute Joule allowance; rounded
absolute material enthalpy is diagnostic and never supplies the increment.

The result keeps original reacting derivatives and chemical entropy separate
from the frozen parcel's calorics, sound speed, transport and flow response.
At constant pressure with only quasistatic pV work, the prescribed enthalpy
change can equal supplied heat. Other processes require their own heat, work
and kinetic-energy balance. No reaction rate, throttle flow rate, pressure
controller, reservoir entropy or passivity claim is inferred. Failed chemistry,
transport, balance or cancellation invalidates the current query; the immutable
reference and previously owned results survive. Publication follows the final
material checks and cancellation poll.

The existing single-temperature `equilibrium_chamber_hp` now delegates to the
native inlet-aware HP entry point before using the common original-material
freezing boundary. A source chart jump between the inlet and the declared
product bracket cannot bypass the native continuity admission.

Four material regressions specify signed reacting heating/cooling against an
independent enthalpy root, 1 mJ increments beside positive and negative large
formation references, transport/final-poll failure atomicity, and original
inlet-path refusal through the material HP adapter. These integration tests
were not executed in this session because the complete material dependency
cone requires the unavailable runtime build. The native thermochemical
composition target passed all 52 tests, including seven new pressure/enthalpy
and shared-isentrope cases; that result does not assert full material or
workspace validation.


### Reacting gas storage on the original thermal timeline

`thermochemical_gas::energy::equilibrium::pressure::storage::GasEquilibriumThermalStorage` binds the original material bank and canonical native isobaric Gibbs/HP law to `ThermalLaw`. Cells store `(H-H_cold)/feed_mass` at fixed pressure and conserved elemental ratios; their actual positive cell mass may differ from the chemistry normalization. The existing thermal timeline owns every finite heat transfer and conserved joule. Enthalpy includes quasistatic pressure-volume work; reference mesh geometry does not impose a rigid-volume constraint. Reacting Cp, composition-dependent density and chemical/mixing entropy come from the same original equilibrium sources, while fully resolved material products retain separately identified frozen transport and acoustics.

A global `dT/dh` upper bound reuses the original NASA coefficient interval bound for every eligible species. The stable ideal-gas redistribution contribution to Cp is nonnegative; the admitted product/feed mass discrepancy and outward rounding are included before inversion. Uncertifiable bounds refuse. Every caloric query restarts the same numerical anchor. Direct local source heat must close the stored increment plus reported endpoint inverse residuals within the original absolute H allowance. The source inverse-temperature mean must lie between the actual endpoint reciprocal temperatures with only explicit floating-point allowance; unresolved entropy intervals refuse without clamping.

The optional `GasEquilibriumConductionLaw` attaches an explicitly declared bulk k(T) curve covering the complete caloric domain to the existing implicit, isolated-equilibrium and harmonic heat owners. Its chain rule uses reacting Cp; density-aware geometry uses `d(1/rho)/dh = (alpha/rho) dT/dh`, including changing molecular count. Nonzero derivative underflow refuses. No reacting conductivity derivative is inferred from frozen transport, and no compatible stream enthalpy datum or species-advection law is invented. Complete NASA support and cold/hot material resolution are prepared; accepted full material queries recheck original transport/source claims, including unsupported interior points.

Five focused regressions exercise analytical reacting calorics/entropy, finite original-timeline heating and cooling with independent cell mass, dissociation volume and conduction derivatives, source/transport refusal, and cancellation/retry. These new material tests were authored and independently source-reviewed but could not be executed after the native executor became unavailable.

### Sealed constant-volume equilibrium and conserved parcel chemistry

`energy::equilibrium::volume` connects the original TV/UV chemistry owner to
the same source-bearing gas material. `equilibrium_chamber_uv` takes a feed
temperature; `equilibrium_chamber_uv_energy` takes explicit total internal
energy. Both hold the declared volume and original feed mass, allowing molecular
count, pressure and temperature to change. Original formation energy is included
once; constant-pressure enthalpy and frozen heat-capacity formulas are not used
to solve the reacting endpoint.

The common product bank must contain exactly the chemistry's complete NASA-9
models. The existing composition-freezing boundary retains active transport
cards, source packs and external query axes, and refuses unrepresentable positive
traces. Final material pressure is reclosed at the exact conserved mass/volume;
the adjustment must fit the declared chemical volume/mass tolerances. Actual
mass, density, target energy, freezing energy and pressure discrepancies remain
in the result. Transport and acoustic properties are resolved at this final
temperature/pressure/composition, with no stale feed coefficients.

`ClosedGasParcel::equilibrate_adiabatic` performs the same operation atomically
on an existing rigid parcel. Its mass and volume retain their exact original
bits, and its actual stored internal energy supplies the target. The complete
bank must retain the active feed's exact caloric/transport models and source-pack
identities; it may supply initially absent products through the declared element
constraints. Source refusal, failed balances or cancellation preserve the whole
original parcel. The resulting storage and acoustic composition is frozen;
the reacting equilibrium response is returned separately. This is a bounded
thermodynamic relaxation endpoint, not a time integrator, ignition criterion,
finite-rate combustor, detonation model or experimental validation.

Focused G0/G1/G3 regressions cover an independent exothermic isomer UV root and
pressure rise, transport/acoustic consistency, reacting versus frozen heat
capacity, exact parcel mass/volume preservation, extensive scaling, source
identity, impossible energy brackets and late cancellation.



### Generated amounts and finite-rate gas reactors

`ThermochemicalGasModel::with_mole_amounts` accepts nonnegative amounts in the
bank's canonical species order and reuses the original exact-simplex and
source-preserving composition owners. Largest-amount scaling avoids overflow
while normalizing extensive inventories. Exact zeros may disappear; a positive
amount that cannot retain a positive mole fraction refuses. The resulting
immutable model retains active NASA-9, transport, water and source-pack identity.
This constructor does not change parcel mass or energy.

The L3 `fs-gas` optional `material` feature supplies the actual finite-rate
consumer through `advance_closed_parcel_adiabatic`. It integrates a declared
reversible elementary mechanism for a specified elapsed time, then uses the
existing material composition/energy replacement to publish one common gas
state. Original stored internal energy is the target; no extra reaction heat
is applied. Mass and volume retain their original bits. Source identity,
actual material balances, zero-duration identity and final cancellation are
checked before publication. Kinetics stays in L3; this crate gains no time
integrator dependency. The endpoint's acoustic composition remains frozen.

### Stable simultaneous species and additional-energy exchange

`species::network::exchange_candidate_with_heat` combines original species
streams with antisymmetric signed additional energy before one reconstruction
per affected finite gas. Heat-only nodes use the same owner. Additional energy
is separate from the existing transported partial enthalpy and may supply a
molecular energy correction; prescribing it does not establish a physical rate
law or entropy guarantee. Canonical endpoint/amount ordering and compensated
summation retain deterministic aggregation. A shared input bound covers mass
and energy records, and original gross species-withdrawal limits remain intact.

At each original node temperature T0, the thermal target uses incoming original
species enthalpy increments plus net species mass times R_i*T0, then additional
energy. The target composition is reconstructed once through the original
NASA temperature-increment owner. This cancels conserved species reference
constants before arithmetic, preserving counterflow heat when absolute formation
energies are enormous. Imported species must support T0 and the declared common
caloric bracket; source extrapolation and missing identities refuse atomically.

The local energy residual is the original-law sensible storage minus that
thermal target. The global residual sums these local closures after the paired
species/energy ledger cancels analytically. Actual total and per-species inventory
rounding is independently admitted against the intended mass changes. Displayed
absolute energy/enthalpy accounts retain their source references and are not
subtracted to judge a small caloric residual. The existing no-additional-energy
entry point delegates to this same stable owner.

Seven new regressions cover heat-only nodes, combined mass/heat, canonical order,
rollback, extreme reference shifts, zero-net-mass counterflow heat and changing
mixtures with distinct species datums. Final files passed Rustfmt parsing and
independent source/equation review. A final test typecheck was blocked before
source analysis when the executor's genuine core dependency archives disappeared;
these seven tests are not claimed as compiled or executed.


### Reciprocal free-molecular species and thermal transport

`species::effusion` samples both original Maxwellian molecular beams through a
qualified fixed thin aperture. Species i has directional molar rate
A*tau_i*p_i/sqrt(2*pi*M_i*R*T). Reciprocal transmissions are explicitly authored
or bound to original selected dimensionless constant source claims. The same
probability applies in both directions and is independent of molecular energy
and internal state. Gas caloric/transport/pack identities and active water
binding are checked directly without rebuilding a calibration composition.

Crossing energy is original partial h_i(T)-R*T/(2*M_i). The monatomic limit is
2RT per mole. Polyatomic use explicitly assumes equilibrated internal modes
independent of translational velocity. Net energy uses original relative NASA
enthalpy increments for paired counterstreams, so thermal transfer survives
large arbitrary formation-energy references. Scaled binary products retain
representable small positive rates; true underflow/overflow refuses. Species
and caller-input counts are bounded before sorting.

The authored qualification covers all mixtures/subsets of the retained original
gas catalog over its complete T/p support. Aperture diameter/minimum mean free
path <=0.1 and wall thickness/diameter <=0.01 are explicit approximation gates,
not experimental validation or rigorous error bounds. Reservoirs are stationary,
well mixed and locally Maxwellian; interface collisions, wall accommodation,
adsorption, energy filtering, bulk drift and moving-wall work are excluded.
A long accommodating pore requires its own wall-energy/transport closure.
Zero exposed area queries no inactive sources. A selected zero transmission
still admits its actual source domain at both reservoir states.

`EffusiveInterface::admit_parcel` rechecks the original source identities and
T/p qualification without sampling another beam. Both `GasTransportRates` and
`GasTransportStep` apply this admission to final coupled endpoints, including
both sides of each reciprocal transmission. An intervening heat, reaction or
mechanical update cannot publish gas outside the aperture's original support.
Zero exposed area stays inactive; selected zero transmission retains its source
domain. This does not impose pressure ordering on thermal transpiration. Three
focused network regressions and source-bound face assertions cover initial-to-
final and later coupled-domain failures, source replacement, cancellation and
geometric closure. Native execution remains blocked by missing constellation
sibling checkouts in the current environment; Rustfmt parsing passed.

`network::transport::GasTransportPaths::effusion` retains both directional
species withdrawals and their molecular energy correction. The correction uses
the actual rounded integrated beam masses and a stable temperature difference.
It joins every other mechanism's original-inventory allowance and the ONE
shared species/energy reconstruction. No pressure-order gate is applied to
thermal transpiration. The original final entropy admission remains in force.
In the pure-gas limit p_left/sqrt(T_left)=p_right/sqrt(T_right) gives zero net
mass while unequal temperatures still exchange energy. Separate apertures may
be composed with other paths; the same exposed area must not be counted twice.

Nine face regressions and five finite-network regressions cover independent
monatomic/polyatomic energy, mass dependence, thermal balance, reversal, large
reference shifts, source/geometry/domain limits, tiny-rate representability,
finite analytical temperatures, uphill-pressure mass transfer, shared nozzle
withdrawal and atomic closure/cancellation. Final files passed actual Rustfmt
parsing and independent equation/source/API review. The final combined test
typecheck was blocked by missing genuine core dependency archives; these
fourteen new tests are not claimed as compiled or executed.

The earlier wet/DGM/pore/adsorption selection passed test typechecking at commit
642daa253c4136a65787ad22eae82c20c7344438, and nine standalone wet-law tests ran
successfully. Full-workspace validation and execution of the expanded coupled
suite, including subsequent network/effusion changes, remain unverified.

### Effusion into an absorbing vacuum

`ClosedGasParcel::effusive_flux_to_vacuum` uses the same original molecular
beam and source qualification with exactly zero incident flux. Vacuum has no
invented temperature, dilute trace composition, pressure floor or finite
receiver. Each species escapes at its own molar-mass-dependent rate with the
original `h_i(T)-R*T/(2*M_i)` crossing energy. Polyatomic caloric behavior still
comes from the original NASA source, including temperature-dependent heat
capacity. Zero geometric area stays inactive; zero transmission retains source
admission.

`network::transport::effusion::vacuum::step_candidate` composes independent
vacuum outlets with the existing nozzle, diffusion, porous, membrane, dusty-gas
and reciprocal effusion paths on one original snapshot. Outgoing boundary
species share every internal path's original gross-withdrawal allowance. The
same source-relative boundary owner applies `dm_into=-dm_out` and the explicit
`+dm_out*R*T/(2*M_i)` correction before one caloric reconstruction per gas.
Thus remaining gas cools and fractionates without subtracting large formation
energies, borrowing newly arrived mass or applying nozzle enthalpy to a
molecular beam. Original source/qualification checks also apply to the actual
combined endpoints and later coupled endpoint admission.

The returned open-system account records escaped species and energy explicitly;
gas-only entropy can decrease and no vacuum reservoir temperature or closed
entropy certificate is invented. If no outlet emits, the original closed
window audit remains in force. This is an explicit first-order, fixed-volume,
locally Maxwellian model with caller-controlled temporal refinement. It does
not model pump performance, receiver filling, long accommodating pores, bulk
drift, a resolved jet or complete depletion of a finite gas inventory.

Focused G1/G3/G4 tests cover independent monatomic, diatomic and variable-heat-
capacity finite-step balances, reference invariance, binary fractionation,
convergence to analytic adiabatic blowdown, shared internal/vacuum withdrawal,
source-domain limits, canonical ordering, closure and atomic cancellation.
`effusion_blowdown` is a runnable CSV example using the existing explicitly
synthetic mixture-file loader and caller-authored aperture qualification.
Rustfmt parsing and source/API review passed; these new Rust numerical tests
and the example were not executed because required constellation sibling
checkouts are absent in the current environment.

### Adaptive finite effusion networks

`network::transport::effusion::trajectory::advance_candidate` advances a complete
requested interval through the original simultaneous gas transaction. Internal
transport, independent vacuum outlets, prescribed molecular reservoirs and
qualified drifting reservoirs all resample the actual finite gas state at every
native step. Geometry, reservoir states and authored relative velocities remain
fixed during the interval. The existing moving-vessel owners remain responsible
for motion, force and inertial energy transforms.

One full native step is compared with two half steps. Only the accepted halves
enter the endpoint and physical accounts. Dimensional local tolerances cover
temperature, pressure, every species inventory (including newly supplied
species), gross directional species streams and source-relative energy. The
second half-step's boundary energy is translated into the trial's original
partial-enthalpy chart before comparison. A large formation-energy datum cannot
hide a thermal error. This is a first-order local error estimate, not a global
accuracy certificate; tighter tolerances require temporal convergence checks.

Typed withdrawal, caloric-solve and internal-flow equilibrium/entropy overshoot
refusals can reduce the trial interval. Invalid sources, unsupported regimes,
unrepresentable arithmetic and cancellation remain hard refusals. Native calls,
accepted intervals and retained history rows are bounded. Any refusal leaves
every caller-owned parcel unchanged, including after earlier internal steps.

The complete interval checks actual original/final total and species inventories
against accepted streams. Actual endpoint internal energy is independently
checked in a fixed reference for each exact original caloric source; original
enthalpy/internal-energy increments remove formation datums before arithmetic.
The accumulated native caloric residual is retained and checked separately.
All conservation budgets apply once to the entire interval, without multiplication
by the number of steps. Open paths make no closed-system entropy claim.

With all external-boundary slices empty, this same trajectory advances closed
`GasTransportPaths`, including the native `layered_dusty_gas` interfaces and
separate molecular paths. Every accepted half step retains its original
laminate source/flux report. Layered interfaces count toward the shared
connection limit when vacuum or held-reservoir paths are composed, and inactive
layered reports count toward the retained-history bound. No second controller,
flux model or energy update is introduced. Independent equal-molar-mass binary
relaxation through parallel laminate/molecular exposures checks temporal
refinement, accepted-only throughput and fixed whole-window conservation.
Native laminate fixtures also check serial Knudsen inventory/enthalpy,
multicomponent source/reference-energy parity, shared withdrawals, later
transport of arrivals and transactional source/work/accuracy/cancellation
refusals. These regressions are authored; native execution remains unverified.

G0/G1/G3/G4 fixtures cover monatomic and diatomic analytic vacuum blowdown over
three tolerance rungs, shared gross-withdrawal retries, deterministic caloric-
reference invariance, new-species supply through two finite chambers with vacuum
fractionation, source refusal, work exhaustion and late cancellation. A
bidirectional drifting-reservoir fixture checks refinement against two native
reference resolutions with a large caloric datum. An intentionally perturbed
endpoint must fail the independent storage-energy audit even when its original
native receipts retain small closure residuals. The
`effusion_blowdown` example accepts a final `fixed|adaptive` argument; fixed mode
remains the default, and adaptive mode treats `dt_s` as its CSV reporting interval.
These are numerical checks with synthetic material sources, not measured
validation of an aperture or collision regime.

### Finite effusion receivers starting at exact vacuum

`energy::parcel::evacuated::EvacuatedGasVolume` holds only the finite volume and
auxiliary material-query axes. It has no gas temperature, composition or seed
inventory; temperature/pressure query axes are refused. Positive received
species inventories and an explicitly temperature-referenced internal energy
initialize one ordinary `ClosedGasParcel` through the original composition,
caloric, EOS and transport owners. The accepted total mass is stored directly.
Actual represented species masses, EOS density times volume and source-relative
storage are independently checked against fixed physical budgets.

`network::transport::effusion::filling::step_candidate` captures the exact
accepted species and molecular-energy corrections from the original donor's
vacuum transaction. The receiver energy relative to the incoming species'
internal energies at donor temperature is `sum_i dm_i R T_d/(2 M_i)`; the
implementation consumes the native receipt rather than reconstructing a second
beam law. A pure constant-Cv source therefore produces a first-fill temperature
`T_r/T_d = 1 + R/(2 Cv_molar)`: 4/3 for a monatomic gas and 6/5 for a diatomic
constant-capacity fixture. Mixtures fractionate according to their original
species rates and transmissions. Zero receipt preserves an exactly evacuated
receiver, including zero-area and zero-transmission cases.

The existing closed-inventory/entropy owner now admits this zero-inventory
initial receiver directly. Its original source increments and zero-amount
mixing limit give actual pair mass, energy and entropy changes without a
logarithm of zero or subtraction of arbitrary source datums. A large initial
step that overfills a small receiver can conserve energy while decreasing
entropy; it is refused with a typed `EntropyDecrease`. `filling::audit_pair`
also checks later ordinary gas endpoints against the original donor plus
vacuum, so a consumer can retain one conservation budget over many reports.

`filling::trajectory::advance_candidate` compares a full first-fill interval
with a half fill followed by a reciprocal finite-gas half step. Receiver
backflow participates immediately in state, gross species-stream and relative
energy error estimates. The original incident molar rate, receiving volume and
maximum temperature in the receiver's solve bracket also bound the initial timestep, so a
coarse trial cannot overfill a small receiver before source admission. This is
a startup time ceiling, not a pressure floor or a restriction on subsequent
thermal transpiration. Only the refined pair is accepted. Every remaining
interval uses the existing adaptive network driver. Native calls, accepted
intervals and retained rows share bounded effort; initialization history is
reserved before continuation retains more rows. The complete original-to-
final pair is independently checked against one fixed mass/energy/entropy
budget. Refusals and late cancellation publish neither endpoint.

These are lumped, thermalized, well-mixed, stationary-gas endpoints in fixed
rigid vessels. Complete capture and equilibration of the molecular beam are
authored assumptions; collisionless startup, molecular flight delay, adsorption,
pumping, wall heat accommodation and moving-wall inertia are not resolved.
Actual positive receiver states must satisfy their original material and
aperture support. A source lower-pressure limit can prevent temporal refinement;
no pressure floor or extrapolated material state is substituted. The local
step-doubling estimate is not a certified global accuracy bound.

The `effusion_receiver` CSV example starts with an exactly evacuated second
chamber, then continues reciprocal adaptive transport through ordinary gas
owners. It uses explicitly synthetic source coefficients and authored aperture
qualification and an explicitly authored synthetic minimum pressure. Existing
synthetic-loader callers keep their original pressure support. Empty receiver
thermal/composition columns are blank. Added
numerical fixtures cover molecular heating, selective mixture filling, source
reference invariance, free-expansion entropy, entropy overshoot, reciprocal
handoff, temporal refinement, source refusal and atomic cancellation; they are
not empirical validation of a physical aperture or equilibration mechanism.

## Source-bound Prony histories and exact stress control

`visco::stress_control::StressControlledMaxwell` extends the existing
`GeneralizedMaxwell`/ `PronyState` owner to exact piecewise-linear stress
ramps, holds and zero-load recovery. The symmetric positive-definite transformed
branch system is diagonalized once through the existing Jacobi eigensolver.
Work integrates stress times strain rate; dashpot heat is integrated
independently. State admission checks the constitutive stress, finite values,
positive storage and numerical energy balance before mutation. At most 64
nonnegative branches are admitted, including empty or entirely zero-weight
elastic spectra. Instantaneous stress jumps have zero viscous heat.

Zero-modulus arms retained by `fit_prony_history` keep their original order
and state coordinates. They are excluded from the energy eigensystem and add
no stiffness, stress, storage or heat. Their kinematic histories integrate the
same stress-controlled strain rate by analytic convolution, including coincident
relaxation rates; no small positive modulus or endpoint-strain approximation is
substituted. Prescribed-strain and stress control can therefore alternate on the
original fitted `PronyState`, even after its absent-arm coordinates become
nonzero. G1 regressions cover fit-to-load-control, an independent convolution
quadrature, load reversals and long holds, with no fitted-material validation
claim inferred from these synthetic checks.

Total strain integrates the same modal rate used for work, avoiding stress
subtraction divided by a weak equilibrium modulus. Whole-ramp strain bounds sum
monotone integrated initial and forcing terms, retaining finite long-hold
limits. Modal roundoff is closed through the instantaneous modulus only within
a bound derived from uncancelled reconstruction terms; strict constitutive
coherence and the independent energy ledger must still pass. This numerical
closure adds no physical load step, work or heat. Zero-duration no-ops preserve
every admitted state field bitwise.

`visco::source::ResolvedPronyReference` consumes dimensioned scalar
`prony-branch-count`, `prony-equilibrium-modulus`, and indexed
`prony-modulus-N` / `prony-relaxation-time-N` claims under the original typed
query and selection policy. It verifies and retains every source receipt and
unstated uncertainty. Optional density supplies mass only when present.
Sampled curves, missing arms, wrong dimensions, invalid spectra, unowned query
coordinates and mismatched bound original geometry refuse.

The specimen starts relaxed in incremental coordinates about the fixed source
operating point. Its nominal stress/engineering strain and work use original
area/length/volume. Strain and force controls share one persistent history and
allow holds, unloading and reversal. Positive duration and authored strain,
absolute strain-rate and elapsed-time limits are mandatory. Stress-controlled
whole-ramp bounds are conservative floating-point interval sums, not certified
outward-rounded enclosures; a large step can require subdivision. Rejected
steps preserve the entire previous state and identity. Accepted responses expose
cumulative and incremental model work and dissipation, current storage, energy
residual and a history identity. Unloading can return work; dissipation remains
nonnegative. A closed strain loop can retain branch storage.

The runnable `material_viscoelastic_reference` example requires explicit
`--time-domain-reference`, bounds input/step counts and emits output only after
the complete history is admitted. Source-fitted coefficients do not independently
validate the simulated history. No unknown preload stress/energy, nonlinear
crushing, densification, permanent set, damage, acoustic absorption or thermal
evolution is inferred. The Deng foam's dynamic calibration and the printed TPU
model-versus-measured cycle-work discrepancy remain explicit in the seed catalog.
Focused G0/G1/G3 tests cover analytic creep/recovery, independent work/heat
quadrature, subdivision, source transport, mixed controls, identity and atomic
refusals; the existing strain-ramp energy regressions remain applicable.

## Prepared thermal transport over source temperature paths

`thermal::transport::PreparedThermalTransport` binds density, caller-selected
specific heat and directional scalar conductivity from one immutable material
card to a positive, explicitly typed absolute-temperature interval. Independent
claim selections and every fixed non-temperature condition remain explicit.
Canonical `ClaimSet::query_line_typed` supplies complete support and original
source knots/receipts; the prepared union grid never smooths across a breakpoint
or silently changes the source model. The union is bounded by three canonical
line budgets. Every selected source explicitly bounds temperature; value-axis
dimensions must match the typed temperature or typed fixed static coordinates.
Second source temperatures and nonzero frequency/rate coordinates refuse.

Runtime evaluation uses one binary search, no allocation, no hashing and no
database query. It interpolates rho, c and k separately before computing rho*c,
k/(rho*c) and sqrt(k)*sqrt(rho*c). The latter avoids unnecessary intermediate
product overflow. Preparation validates source support and knot values;
unrepresentable interior results still refuse without publishing a sample or
mutating the model. Splitting source segments on the union grid can change
roundoff relative to a direct query, not the underlying piecewise-affine law.

The identity binds the original card and ordered source receipts, including fixed
conditions and selection. Original uncertainty remains on the retained answers;
derived values do not acquire an invented measured claim or uncertainty bound.
No isotropy, Cp/Cv conversion, caloric integration, phase change, evolving
moisture, covariance or coupled mechanics is inferred. Re-prepare after changing
the card, fixed conditions or selected claims. Focused integration tests cover
numerical evaluation, source selection, typed refusal and the real NASA 316 pack.

### Oriented principal-direction thermal sources

`thermal::transport::orthotropic::PreparedOrthotropicThermalTransport<D>`
composes two or three original directional paths into an explicitly declared
orthotropic law. All paths must retain the same original material card, fixed
query context and complete density/specific-heat source lines and receipts.
The caller assigns principal directions; scalar observations do not establish
material symmetry or orientation. A proper orthonormal basis is admitted
without repair and retained in the prepared identity.

Runtime evaluation reuses the original scalar interpolation and support,
then forms `K(T)=Q diag(k_i(T)) Q^T` and `alpha(T)=K(T)/(rho*c)` in SI.
Signed cross terms survive; no conductivity averaging, extrapolation or new
measured-property claim is introduced. Every original directional source
remains accessible on the prepared owner. A constant-storage view is exposed
only when every original density and specific-heat knot is unchanged over
the entire prepared path, including hidden interior excursions.

Three focused native source tests pass for rotation/receipts/support,
interior storage excursions and three-dimensional cross-plane transport.
This is a source constitutive view; consumers still admit their represented
discrete tensor, heat storage and physical evolution. In-plane D=2 does not
infer a transverse-flux reduction of an unknown 3-D material.

## Unified gas enthalpy transport

Bulk pair transfer, bulk snapshot networks and selective pair exchange share
the original species-network caloric reconstruction. Conserved species
formation datums cancel before sensible enthalpy/internal-energy increments
are evaluated; ideal-gas pressure-flow work remains explicit. Each affected
volume is reconstructed once from its simultaneous original-inventory fluxes.
This reference invariance applies to nonreacting transport with unchanged
species source conventions; it does not remove reaction formation energy.
Displayed extensive changes retain their original reference convention, while
conservation uses local datum-cancelled energy residuals. Actual total and
per-species inventory rounding remain separately bounded. Bulk-edge and public
selective-network work limits, source conflicts, gross donor allowances,
cancellation and atomic publication remain enforced.

The sequential nozzle network retains its stable-ID Lie splitting. Its energy
audit telescopes the accepted original-caloric link residuals, since each staged
endpoint becomes the next link's initial state; actual final total and species
masses are independently checked. Per-parcel energy reports accumulate the
same accepted storage increments. The global residual is not reconstructed by
subtracting large rounded formation-energy reports. This is a balance of the
chosen first-order split trajectory, not an entropy or temporal-error bound.

## Source-relative external gas and adsorption energy

`ClosedGasParcel::exchange_species_boundary_relative` accepts each species'
complete integrated energy as `dm_i*h_i(T_original) + Q_relative_i`. The external
owner supplies `Q_relative_i` directly from its original law; it is never
recovered by subtracting rounded absolute stream totals. The fixed-volume gas
endpoint solves the original final-mixture caloric increment against
`sum(dm_i*R_i*T_original + Q_relative_i)`. Pressure-flow work is included once,
without freezing heat capacity or discarding reaction energy. This API covers
nonreacting species exchange with explicitly compatible source conventions.

The complete retained catalog supports exact individual-species depletion and
reentry. Total vacuum, unavailable or mismatched sources, unrepresentable mass
changes and unsupported caloric states refuse atomically. Actual total mass,
actual species changes versus requested transfers, and reconstructed inventory
versus intended inventory are separately bounded. Zero-net-mass thermal
exchange is permitted. Absolute reporting fields retain the original energy
datum; their rounded subtraction is not the authoritative thermal residual.
The legacy absolute-energy API remains available and cannot recover precision
already lost when a caller constructed its inputs. The partial-enthalpy-change
query evaluates the original NASA-9 source; it does not certify source-region
continuity or invent a phase/latent reference convention.

Single-site and competitive adsorption pass their original-law desorption
enthalpy correction through this relative boundary. Finite solid caloric
storage, binding heat, gross site/species inventories, equilibrium-crossing
checks and combined entropy remain under their existing owners. The combined
energy audit sums the gas and solid thermal residuals before large formation
constants enter display reports. Source alignment and independently admitted
binding enthalpy remain required; a reference shift cannot manufacture source
precision or validate an empirical adsorption model.

## Source-relative molecular phase exchange and histories

Solid/vapor surface exchange evaluates original-source enthalpy differences
before formation constants are rounded into stream totals. The existing
Maxwellian molecular energy remains `h_i(T) - R_i*T/2`; its kinetic correction
and reflected heat enter the same gas reconstruction exactly once. The shared
solid update uses the admitted sublimation latent anchor together with the
actual source-alignment residual. A compatibility tolerance does not turn the
solid/vapor enthalpy gap into an exact identity or erase its integrated energy.
Mixed solid/liquid storage contributes only bulk enthalpy minus pure-solid
enthalpy; fusion already present in the original chart is never added again.
The solid update checks the caloric change actually represented by its final
stored enthalpy and independently bounds its actual mass-change roundoff.

Specular, accommodating, equilibrium-phase and joint chamber callers all pass
source-relative species energy to the original finite gas owner. Reflected
carrier heat retains its species identity and zero net mass. Gross incoming
and outgoing inventory limits, original source domains, latent alignment and
explicit exposed-solid area remain required. Reference invariance is limited
by the precision of those independently admitted sources and stored states;
an extreme datum cannot manufacture accurate latent or binding source data.

Chamber and thermal-timeline histories accumulate only accepted original
caloric residuals and storage increments. Rejected trials and discarded coarse
solutions contribute nothing. Combined, gas-stream and solid-stream defects
remain subject to the declared whole-history energy tolerance; actual total
and species mass checks remain independent. Thermal heat and incoming birth
energy retain their original owners and explicit datum offsets. Absolute
storage/stream reports retain reference energies, while their rounded
subtraction is not used to reconstruct small heat or to accept the trajectory.
These changes preserve the existing explicit kinetics, splitting, source
support and molecular assumptions; they do not introduce a new evaporation,
reaction, accommodation or real-gas model.

## Directional passive plane-stress viscosity

`PlaneStressViscosity` admits positive-semidefinite directional loss while
retaining all six source coefficients, units, identity and frequency support.
`factor_passive_resultant` normalizes the three-component work form and reuses
the tensor Prony PSD owner, including its `64*EPSILON*norm` zero-eigenvalue
roundoff allowance, with reconstruction checked at `1e-12`. Rank-deficient
stress and power use the retained Gram factor; no positive loss floor is
introduced. The established positive-definite path retains its arithmetic.
The native material tests cover axial and coupled rank-one nullspaces, rotated
work, thickness scaling, coefficient scales `1e-200..1e200`, and active-tensor
refusals. This admits declared reduced viscosity; it does not infer a bulk
constitutive law from absent data.

## Direct volumetric transport and free-expansion density

`thermal::transport::volumetric::ResolvedVolumetricThermalTransport` resolves
positive volumetric heat capacity and scalar static conductivity from one card
at an exact typed condition, retaining original answers, uncertainty and receipts.
It shares static-temperature/source-axis admission and factorized effusivity
arithmetic with prepared transport. It derives k/C_vol and sqrt(k)*sqrt(C_vol)
without inventing density, mass-specific heat or independent statistical evidence.
Refresh pins the original claims and requires the original card and coordinate
schema; a failed update leaves the old snapshot intact. A point observation does
not acquire interval support or become a caloric integral.

`thermal::transport::expansion::PreparedFreeExpansionDensity` admits a positive
reference density and continuously supported engineering free length changes.
The caller explicitly declares isotropy or fixed orthogonal principal directions,
constant mass and unconstrained free expansion. Each length ratio is normalized
at the density reference; the finite volume ratio is their product and current
density is reference density divided by that product. The model retains original
source lines and density evidence, fixed context, basis/frame, selections and
nonzero source datum residuals. CTE, infinitesimal tensor strain and unknown
temperature support cannot substitute for the declared source roles.

The union grid includes the reference temperature. Evaluation uses one binary
search and no allocation, hashing or source query. It derives stretches, volume
and density from interpolated engineering strains, rather than interpolating the
derived density. Positive stretches, volume and density must remain representable;
an extreme interior refusal publishes no sample and leaves the path unchanged.
There is no inferred mass transfer, constrained strain, phase law, evolving frame,
cross-source specimen match or propagated uncertainty.


## Source-backed fixed-volume heat capacity

`thermal::transport::isochoric::IsochoricHeatCapacityLaw` binds an explicitly
selected single-phase mass-specific Cv(T) source to the existing thermal
timeline. The source must continuously support the actual fixed mass density
and complete typed temperature path. Adding a density coordinate to an
unqualified query cannot establish source support. Fixed pressure/stress
coordinates refuse; a Cv curve evaluated along an isobar is not an isochore.
All other composition/source conditions remain fixed.

The shared capacity kernel integrates the original affine Cv exactly for
internal-energy increments and Cv/T for relative entropy. Its existing
prepared inverse provides logarithmic source lookup and analytic segment
inversion, with bounded direct fallback at unrepresentable search coordinates.
Runtime queries allocate nothing and neither query nor resample the source.
Source receipts, card identity, density, coordinates and energy datum remain
bound to the immutable law. Cell admission checks mass/volume against this
density, allowing only arithmetic roundoff. `ThermalStorageBasis::InternalEnergy`
keeps heat-only coupling available while the existing mass-advection owners
reject treating this storage as stream enthalpy.

Seven G0/G1 regression tests specify independent affine-Cv energy and entropy
integrals, heating and cooling across knots, inverse/domain boundaries, retained
source conditions, native heat-driven energy conservation, incompatible-volume
refusal and open-stream refusal without publishing a candidate. Indexed forward
queries additionally retain small relative changes and canonical source support.
This is numerical single-phase calorimetry, not an EOS, pressure prediction,
phase transition, chemical reaction, evolving strain model or experimental
material validation.

### Relative internal energy and direct-temperature conductivity

`IsochoricHeatCapacityLaw::temperature_at_relative_internal_energy` exposes
the existing preserving-energy affine-segment inverse relative to an admitted
reference temperature. Consumers must retain the supplied energy coordinate;
a sub-temperature-ULP increment may return the rounded reference temperature
without erasing the stored energy. Original source limits remain enforced.
This avoids absolute-datum cancellation without changing Cv, receipts, density,
energy identity or the separate ordinary heat-update API. Its bounded source
traversal can cross multiple original knots.

`IsochoricConductionLaw::conductivity_at_temperature_w_m_k` queries the SAME
prepared conductivity path at an already resolved temperature. It adds no
material selection, extrapolation or caloric inversion. Relative-energy
consumers therefore avoid a redundant T-to-absolute-u-to-T round trip.

`IsochoricConductionLaw::state_at_temperature` returns the same original
conductivity value and its caloric-coordinate derivatives at the exact supplied
admitted temperature: `dT/du = 1/Cv(T)` and `dk/du = (dk/dT)/Cv(T)`. It does not
round-trip through absolute energy. Interior conductivity knots retain the
original right derivative; the upper endpoint retains the left derivative.
Both caloric and conductivity support and representability refusals remain
active. `state_at_specific_internal_energy` resolves T once and delegates to
this same source path. Regression coverage compares exact state bits under
two caloric energy datums, analytic derivatives across all original segments
and knots, and immediate-outside/NaN/infinite temperature refusals.

Two focused regression definitions cover relative inversion across knots,
retained sub-ULP heat, exact support refusals and original conductivity at
knots/interiors. Pinned rustfmt parsing succeeds; native execution of these
new tests remains pending the workspace's original dependencies. The consuming
jet owns evolving mass, transport and conservation; these source methods do
not supply general stream enthalpy, pressure work or a new thermal integrator.

## Source-integrated conduction between original Cp and Cv regions

`thermal::transport::conduction::CaloricConductionLayer` retains an original
`IsobaricConductionLaw` or `IsochoricConductionLaw` and prepares its unchanged
conductivity knots once for the existing Kirchhoff integral owner. Many faces
borrow that same prepared layer. `IntegratedCaloricLink` binds the actual
ordered caloric objects, contact coverage and planar/cylindrical/spherical
geometry, then supplies the existing implicit chain or graph with the complete
source-integrated flux and analytic derivatives with respect to each original
specific h or u. The exact Cp/Cv integrals, source identity, density admission
and entropy remain with the original caloric owner.

The original two-layer kernel resolves both material-face temperatures and a
finite contact jump. Its dQ/dT derivatives chain with 1/Cp(T) or 1/Cv(T), rather
than differentiating a frozen effective conductance. The native residual's
actual secant power must agree with that same constitutive solve. Successful
face evaluations allocate nothing and query no material card; preparation is
bounded by the existing conductivity-curve knot limit. Point-only Cp/Cv paths
keep their previous admission, arithmetic and identity conventions.

Either side may use Cp or Cv; a heat transfer does not identify h with u or
supply unmodeled pressure-volume work. The layer profile is source-free,
single-phase and quasi-steady within each leg, at fixed source context.
The ordinary query retains its fixed planar/radial geometry; the explicit
density-aware planar extension below includes changing normal half-distances.
This adds no separate half-cell storage, latent account, flow,
mechanical work or time integrator. Original source-domain, contact, geometry
and cancellation refusals publish no candidate. Spatial and temporal accuracy
still require the consuming discretization's independent refinement.

Four focused native regressions pass against independent quadratic flux and
caloric references: unequal Cp/Cv source paths and analytic contact tangents;
planar/radial, equal-temperature and reversed heat flow; a large actual implicit
step with original stored-energy and entropy closure; and original-object,
geometry, source-domain and late-cancellation refusal with exact replay. These
are numerical model checks on synthetic source cards, not measured-material
validation.

## Explicit current planar geometry on original integrated caloric faces

`IntegratedCaloricLink::with_current_planar_geometry` prepares a
`CurrentPlanarCaloricLink` borrowing the original integrated face. The caller
supplies the current common area and both physical normal half-distances,
then installs those same measures on the original native observer using
`CurrentConductionGeometry`. Preparation checks the original ordered caloric
objects, reference contact identity and per-side coverage, admits the source
cells, and probes the actual current face before returning. Every native trial
checks that its current measures and immutable reference coverage match the
prepared binding. Original reference geometry, mass, h/u basis and heat storage
remain unchanged.

The wrapper calls the existing two-material Kirchhoff solve with those current
measures and the original conductivity curves. An authored areal contact
resistance remains areal; an original source contact is evaluated once per
trial at the same two energies and supplies the complete resistance. Its exact
`dQ/de` contribution remains part of the original native tangent. No material
card is re-queried, source support widened, conductivity table copied, or
independent heat inventory introduced.

`IntegratedCaloricBoundary::with_current_planar_geometry` similarly prepares
a `CurrentPlanarCaloricBoundary` for an original reservoir or exposed-ambient
channel. It rebuilds the original planar half-cell geometry and converts the
unchanged areal film resistance to total resistance using the CURRENT area,
`R_film = R''/A_current`, once during preparation. Retaining the old total
film resistance would incorrectly change that physical film. The original
bounded material-face solve, supported conductivity integral, caloric chain
rule, applied secant heat and bath entropy accounting remain with their
existing owners. A bath outside source support is still permitted only when
the actual material face remains supported.

These explicitly prepared wrappers alone opt into current mechanical geometry;
the original source adapters retain their reference-pinned admission. Each
thermal interval holds the selected mechanical measures fixed. Rebinding
prepares a later configuration without mutating the original source. This
is distinct from density-driven thickness and does not infer geometric energy
derivatives, density changes or mechanical work. Original fixed-density,
fixed-volume, pressure and composition constraints still apply and require
compatible mechanical motion. Radial equivalent distances cannot be rebound
as physical planar lengths; gray, multilayer and other pinned source owners
are not implicitly enabled.

Four focused regression definitions compare quadratic Kirchhoff/contact
oracles and film-face roots with the prepared fluxes and full energy
derivatives; exercise native chain/graph backward Euler and SDIRK2 with mixed
original Cp/Cv storage; check separate reservoir and ambient heat/entropy; and
cover stale geometry, source/reference refusal and cancellation replay.
These new native regressions are authored but unexecuted under the current
workspace dependency/build limits; pinned Rustfmt parsing is not integration
validation.

## Original reservoir response of source caloric exchange

`IsobaricHeatCapacityLaw::linear_heat_exchange` retains analytic derivatives
with respect to the initial enthalpy and reservoir temperature through every
original affine Cp interval. Its stable inlet attenuation, final-enthalpy wall
derivative and mean-temperature wall derivative do not differentiate a frozen
Cp or require material data at the reservoir temperature. The latter may lie
outside source support while all actual material temperatures remain supported.
Positive integral moments preserve weak exchange; exact equilibrium uses the
local exponential derivative limit. The original bounded exchange inverse,
signed heat, source identity, pressure and composition remain unchanged.

The same result exposes `remaining_temperature_difference_k`, the signed
reservoir-minus-outlet drive and derivative of outlet specific enthalpy with
respect to exposure, and `mean_log_exposure_derivative_k`, the signed derivative
of mean temperature with respect to log exposure. The latter equals outlet
minus mean temperature and is zero at zero exposure or equilibrium. Positive
interval moments retain it before absolute-temperature rounding, with no new
source traversal, allocation or perturbed solve. The remaining drive retains
the analytic attenuation even when the absolute outlet rounds to the reservoir.

Focused regressions differentiate an independent piecewise caloric integral
in heating and cooling across source knots, including exterior reservoirs,
and check weak, zero, equilibrium and strong exchange limits. Original affine
caloric integrals and constant-Cp Newton cooling also check both exposure
sensitivities, including NTU `1e-12`, saturated NTU `80`, and a finite mean
response when both absolute outlet and mean round to the same temperature.
These are numerical model checks, not experimental material validation.

## Integrated conductivity for volumetric and assembled storage

`CaloricConductionLayer::from_fixed_volume` retains an original
`FixedVolumeThermalLaw`, and `from_heat_storage` retains an original
`HeatStorageConductionLaw`. These use the same integrated link, film,
gray-boundary and surface-interface owners with `dT/de = rho_ref/C_vol(T)` or
`m_ref/C_total(T)`, respectively. Both the original reference-amount admission
and `DeclaredSensibleStorage` basis survive; no Cp/Cv meaning or conductivity
homogenization is inferred. Each conductivity integral retains only its
original k(T) knots, independently of the capacity grid, and many faces borrow
that one preparation. The public `caloric_storage_basis` regression checks
both storage models through native implicit film cooling against an independent
piecewise-conductivity and energy balance, including the full tangent, entropy,
supported material faces beside an out-of-span bath, original-amount refusal,
replay and cancellation. This new regression is authored but unexecuted under
the current native build limits.

## State-dependent contact on original Cp/Cv integrated faces

`IntegratedCaloricLink::new_with_contact` binds an original two-sided
`ImplicitContactLaw` to the same ordered Cp/Cv caloric objects. Its source
bounds must lie within those original energy domains. The contact supplies
the entire areal resistance; the native geometric binding declares zero
additional fixed resistance, and the existing thermal owner prevents a second
contact attachment on that integrated face.

Every nonlinear trial evaluates the contact at the same original h/u pair as
the caloric inverses. The existing Kirchhoff kernel uses that resistance and
returns its analytic sensitivity. The complete native derivative is
`dQ/de_i = (dQ/dT_i)*(dT_i/de_i) + (dQ/dR'')*(dR''/de_i)`.
A passive state-dependent resistance can make a flux derivative negative;
finite such derivatives are retained rather than replaced with the
fixed-resistance monotonicity assumption. Fixed-resistance links retain their
previous arithmetic and admission.

Preparation probes the actual source contact directly, so finite resistance
can join supported material faces even when their temperature charts do not
overlap. An intermediate perfect-contact solve is not imposed. The original
source bounds, identity, geometry, cancellation and accepted-state guarantees
remain active. The wrapper adds no allocations or independent heat storage;
the original contact source retains its own evaluation costs.
`LayerConductionResponse` temperature derivatives hold the evaluated resistance
fixed; the native link trait includes the complete contact-source chain rule.

Three additional native regressions pass: independent two-sided source-contact
derivatives including passive negative differential flux; actual implicit
mixed-Cp/Cv evolution with stored-energy and entropy closure, late cancellation
after an accepted interval and exact retry; and original source ownership,
bounded/nonrectangular support and disjoint-temperature contact admission.
All 18 focused native Cp/Cv conduction, integrated-link, boundary and contact
tests pass with the original full thermal owner and numerical solver kernels.

## Original Cp/Cv source feedback in conducting and radiating gaps

`IntegratedCaloricGrayLink::new_with_gap_resistance` binds one original
`ImplicitContactLaw` to the nonradiative part of an integrated gray gap. It
shares ordinary integrated contact's ordered caloric identity, supported energy
bounds and immutable source-trial admission. The source supplies strictly
positive finite R''(e_A,e_B) and both original energy derivatives; its reciprocal
must be representable. The declared gap's fixed HTC must be zero. Radiation,
material half-cell resistance and coatings cannot be included in this R''.
The native geometric resistance stays zero and the existing thermal owner
prevents a second contact attachment from owning the same interface.

Every initial or nonlinear trial queries the actual source once at the
original bulk h/u pair. The existing gray-gap kernel solves the two material
profiles and both surface transfers together. Its analytic resistance
sensitivity includes the resulting change in BOTH surface temperatures and
radiation: `dQ/dR'' = -Q_gap * (HTC / D)`, with D the existing coupled surface
derivative denominator. No extra surface solve or finite difference is used.
The native endpoint tangent adds `(dQ/dR'')*(dR''/de_i)` to the original
conductivity-integral and 1/Cp or 1/Cv chain. Finite source-dependent tangents
may reverse sign while strictly positive R'' retains passive signed heat flow.
Fixed-gap links keep their existing differential admission and arithmetic.

`GrayLayerGap::with_nonradiative_resistance` retains the exact authored R''.
Disabling radiation delegates to the original contact kernel with that exact
value, avoiding a reciprocal round trip through HTC. Existing fixed-HTC gaps
keep their previous fallback arithmetic. The response's temperature and
resistance derivatives are partials; the native source-bound link composes
their full original-energy chain. No additional heat store, surface-law state,
source authority or timestepper is introduced. Original source coordinates,
geometry, emissivities and external conditions remain fixed as declared.

Three new native regressions pass: independent radiative resistance and both
original caloric derivatives, including passive negative differential response
and one source query per trial; actual implicit evolution with energy/entropy
closure and late cancellation/replay; and original source/support ownership
with initial disjoint-temperature admission. The radiation-disabled test uses
a resistance whose reciprocal round trip differs by one ULP and verifies exact
recovery of the original contact flux and complete tangent. All 28 focused
native Cp/Cv and integrated exchange tests pass against the original thermal
owner and numerical kernels. Full workspace Cargo/Clippy gates remain unrun
because the declared external constellation dependency directories are absent.

## Source-integrated Cp/Cv radiation between finite bodies

`thermal::transport::conduction::IntegratedCaloricGrayLink` attaches the existing
integrated gray-gap kernel to the original finite Cp/Cv thermal chain or graph.
Both original material half-cells integrate their own k(T) curves to the two
solved material surfaces. The declared nonradiative gap HTC and gray radiation
act in parallel between those surfaces. Their reported powers are components
of ONE equal/opposite heat transfer, never additional heating channels.

Each endpoint tangent uses its own original 1/Cp or 1/Cv and exact conductivity
integral. Enthalpy and internal energy remain distinct inventories. Preparation
shares the original link's ordered caloric identity, density/mass admission,
contact coverage and geometric checks with the existing conductive adapter;
each link borrows the already prepared original source tables. No trial tables,
new thermal state or separate time integrator are allocated by this adapter.

The gray-gap source explicitly declares equal-area opaque diffuse-gray faces
with full mutual view and a nonparticipating gap. Only planar geometry is
admitted, and both physical covered areas must equal the nominal transfer area.
The native geometric resistance is zero because the gray gap owns the entire
interface; an extra authored resistance or separate contact attachment refuses.
Fixed emissivities/HTC and source-free quasi-steady single-phase legs do not
infer finite gap storage, evolving surfaces, radial view factors or phase
fronts. Each material profile stays on its own original source support. Initial
admission uses the actual gap, permitting disjoint temperature charts where a
finite supported jump exists instead of imposing an unrelated perfect contact.

Three new native regressions pass: independent nonlinear material/gray-gap
surface balances and complete mixed-Cp/Cv tangents in either direction and at
equilibrium; actual 300-second implicit evolution with conserved energy,
original entropy and late cancellation/replay; and identity, geometry,
asymmetric coverage, source-domain, disjoint-support, insulation and exact
radiation-disabled contact limits. All 25 focused native original Cp/Cv and
integrated exchange tests pass.

## Live two-surface sources on original finite Cp/Cv networks

`IntegratedCaloricSurfaceLink` binds the source-resolved material surface
kernel to the ORIGINAL ordered finite caloric cells through
`ImplicitLinkLaw`. Each original k(T) profile is prepared once and borrowed.
The original h/Cp or u/Cv storage, cell mass, caloric identity, accepted time
and entropy remain with the native chain/graph. Preparation immediately uses
the actual source and its supported face intervals; it does not require an
intermediate perfect-contact state.

One `ImplicitSurfaceExchangeLaw` evaluates both original specific energies
and both ACTUAL solved surface temperatures. This admits temperature-dependent
contact/gap transport and declared radiative exchange with cross-coupled
source response. Each energy bound lies within its own original material
chart, and exact ordered caloric objects must match. Face support may be
narrower than the material charts and exclude either bulk temperature. The
source must admit its complete coupled query; interval bounds alone do not
authorize extrapolation. Planar, cylindrical and spherical material legs retain
their original exact Kirchhoff resistance coordinates. Both actual covered
areas and the native patch area must equal the declared common interface area.
Radial faces are the two material-side traces at one shared positive interface
radius; the ordered bulk traces straddle that radius. A cylinder uses one
common axial length and a sphere its full solid angle. Changing the cylinder
length without changing the original patch therefore refuses. The same source
owns the entire interface; equivalent radial distances do not replace physical
control-volume sizes or create thermal mass. This does not infer a radiative
view factor for separated concentric faces, partial circumferential contact or
axial end effects; those layouts need their own resolved geometric model.

The complete native energy tangent combines the original material/surface
endpoint derivatives times their original 1/Cp or 1/Cv with the explicit source
partials at fixed faces, scaled by the SAME surface balance's A/D. No bulk
dk/de term is added to the differentiated conductivity integral. Signed
feedback and zero-conductance onset retain both original energy derivatives.
The same native secant applies the ONE accepted heat transfer; no independent
conduction/radiation deposit, extra surface inventory or second clock is added.

`ImplicitSurfaceExchangeState::from_conductance` applies the exact product
rule for j=G(e_A,e_B,Tsa,Tsb)*(Tsa-Tsb), including G=0 onset. Its inputs are
source-owned conductance partials; it infers no material data, heat-conversion
policy, emissivity or view factors. General enclosures continue to use the
existing radiosity owner. Trial evaluation allocates no adapter tables and
retains the accepted source result without re-query; source evaluation keeps
its own costs. Original geometry mismatch, source failure and cancellation
cannot publish a partial native thermal step.

Four new native integration tests pass: live gray/HTC exchange with crossed
energy/face dependencies against independent quadratic-k/quartic calculations
and full original-energy finite differences; actual mixed-Cv/U and Cp/H BE
evolution with energy, entropy and late cancellation/retry; independent face
support plus identity/geometry/domain refusals; and zero-conductance onset with
positive cold-side differential feedback. Together with the preceding kernel
increment, 48 selected native tests pass. Complete thermal production and
focused material sources compile; full workspace gates remain unavailable
because external constellation dependencies are absent.

## Source-resolved exchange between two material surfaces

`TwoLayerConduction::evaluate_surface_exchange` closes both original supported
conductivity integrals against one passive source at the actual ordered A/B
surface temperatures. `SurfaceExchangeState` supplies one flux density,
its passive secant, two signed surface partials and an original payload.
The payload from the accepted query is retained without another source call.
The source owns the complete interface; an additional authored contact
resistance refuses. No surface storage or distributed material source is
introduced.

The existing scalar flux driver now accepts optional independent source
temperature supports. Both ends of both face intervals are mapped through the
original Kirchhoff integrals, forming a supported lower AND upper flux bound.
A face interval may exclude its body's temperature or be an exact singleton.
Material charts can remain disjoint across a finite physical interface.
Upper/lower endpoint and interior callbacks all count against the declared
maximum iteration budget and are bracketed by cancellation checks. No callback
is extrapolated outside its source interval. Source errors retain their original
identity; interval bounds do not certify a coupled source's unsampled domain.

The same safeguarded driver retains both actual material integrals, layer heat
balance and one original A-to-B transfer. The interface derivative is
`D=1+F_A*d_A/k(Tsa)-F_B*d_B/k(Tsb)`.
Endpoint derivatives and the explicit-source response scale `A/D` differentiate
this same balance. Negative source differential slopes remain signed. An
intermediate nonpositive D requests bisection; an accepted singular or locally
unstable elimination refuses. No global uniqueness claim follows for an
arbitrary nonlinear source. Exact insulation, equilibrium and original A/B
ordering remain explicit. Existing fixed-contact and fixed-gray calls preserve
their prior arithmetic and checkpoint flow.

The kernel allocates no trial tables, introduces no new scalar solver and
resamples no conductivity chart. The original finite thermal owner retains
caloric storage, time and entropy; source evaluation retains its own costs.
General radiative enclosures remain with the existing radiosity owner.

Four new native regressions pass: independent affine-k/quartic surface balance
and derivatives in both directions; disjoint/narrow, endpoint and singleton
supports without unsupported queries; signed stable/unstable and zero limits;
and source-error/cancellation/work-budget behavior. All 44 selected original
caloric, integrated exchange and shared-surface tests pass. The complete
thermal production crate and focused material subset compile with unchanged
numerical dependencies; full workspace Cargo/Clippy remains unavailable because
external constellation dependencies are absent.

## Exterior-temperature source feedback on original Cp/Cv storage

`IntegratedCaloricGrayBoundary::with_exterior_law` binds one original
`ImplicitGrayExteriorLaw` to the actual stored energy, solved exterior surface
temperature and scheduled gas-bath temperature. Emissivity and convection can
therefore change with the exterior surface while the body retains its original
nonlinear h/Cp or u/Cv inventory. The original bulk-only source remains a
separate declared closure; its coefficients are still queried once per energy
trial. Fixed, bulk-source and exterior-source ownership are mutually exclusive.

The existing supported conductivity integral, coating jump and bounded common
surface solve remain the physical owners. Every exterior source query uses the
actual trial Ts; the accepted response carries its partial derivatives without
a second query after convergence. Source temperature support can exclude the
bulk and baths, and a coating can place the exterior outside the material
conductivity chart. The material path still requires its own support. The
source must refuse unsupported gas temperatures or coupled source conditions;
interval bounds do not authorize extrapolation or rectangular-domain inference.

For total exterior power F(e,Ts), the complete tangent is
`dQ/de = (F_Ts*Gc*dT/de + Gs*F_e)/(Gs+F_Ts)`, where Gc and Gs are the
existing material endpoint slopes after eliminating any coating resistance.
F_Ts includes both coefficient derivatives, and F_e holds Ts and Tgas fixed.
The original 1/Cp or 1/Cv appears once; an additional bulk dk/de term would
double-count the differentiated conductivity integral. Zero-emissivity/HTC
onset retains its original partials. Signed differential response is admitted
when the combined material/exterior slope is positive; a singular or locally
unstable surface elimination refuses. This local condition does not establish
global uniqueness for an arbitrary nonlinear source.

Separate bath powers, actual bath entropy temperatures, the native chain/graph
and transactional time integration remain authoritative. No extra surface
inventory, evolving flow state, finite radiative enclosure or source authority
is inferred. The adapter allocates no trial tables; the source retains its own
evaluation costs. Focused native regressions cover independent nonlinear
surface balances, full Cp/Cv tangents, distinct source support, coefficient
onset, signed response and accepted-time energy/cancellation behavior.

## Coated Cp/Cv surfaces with separate material and exterior temperatures

`IntegratedCaloricGrayBoundary` now admits the original boundary geometry's
finite nonnegative coating/contact resistance, converted from areal m² K/W to
total K/W once. The existing reservoir-layer owner solves
`Q = M integral(Tmaterial..Tcell, k(T)dT)`,
`Texterior = Tmaterial - R*Q`, and
`Q = Qconvection(Texterior) + Qradiation(Texterior)`.
The two physical temperatures are distinct: material support is enforced along
the original conductivity path, while the exterior surface and prescribed baths
may lie outside that chart. The callback's own exterior-temperature support
remains authoritative. No conductivity or surface-source extrapolation occurs.

`ReservoirLayerGeometry::evaluate_material_to_two_baths_with_resistance` reuses
the original safeguarded power root and caloric-independent conductivity
inverse. At most two existing material/film solves prepare the supported power
bracket; no nested root occurs in its interior trials. The response exposes
`material_face_temperature_k` separately from its exterior surface temperature.
Frozen material/bath and bath/bath pair secants include the same series coating.
The zero-resistance path delegates directly to the prior method and preserves
its arithmetic, callback flow and cancellation checks.

The native Cp/Cv tangent analytically eliminates the coating jump before
combining the exterior flux derivative and original bulk-source coefficient
response. It retains the original 1/Cp or 1/Cv and both conductivity-integral
endpoints. Separate bath powers, their actual entropy temperatures, the original
thermal inventory, and transactional native time integration remain owners.
The resistance represents a fixed zero-storage coating/contact; finite coating
heat inventory requires its explicitly represented thermal state or an admitted
bonded-film storage model. The declared R excludes both the convective film and
radiation. Geometry and source work constraints remain explicit.

All 36 focused native material/thermal tests pass, including five existing
shared-surface kernel tests and three new coated-boundary regressions. New
checks use independent quadratic-material/quartic-radiation balances and full
Cp/Cv/source tangents; a constant-k resistor network and exact zero-R recovery;
actual BE energy/entropy with late cancellation/replay; and hot/cold exterior
surfaces outside the material chart with supported material faces. The complete
thermal crate and original numerical kernels compiled in the focused runner;
full workspace Cargo/Clippy remains unrun because constellation dependencies
are absent.

## Source-integrated Cp/Cv convection and radiation at one surface

`thermal::transport::conduction::IntegratedCaloricGrayBoundary` couples one
original Cp/Cv material layer to distinct convective and radiative baths. It
reuses `ReservoirLayerGeometry::evaluate_material_to_two_baths` to solve the
shared material surface with its exact supported k(T) integral. The native
`ImplicitBoundaryExchangeLaw` receives the two actual signed powers and bath
temperatures, including opposing exchanges when the body receives no net heat.
The existing chain/graph retains all finite-body storage, accepted time, energy
and bath entropy. No effective bath, second half-cell or additional heat source
is inferred. The selected channel supplies the actual scheduled gas temperature;
the radiative surrounding has the explicitly bound fixed temperature.

The full original-energy tangent differentiates both conductivity-integral
endpoints and chains through the original 1/Cp or 1/Cv. Optional original
`ImplicitGraySurfaceLaw` emissivity/HTC derivatives also enter that same tangent,
including signed responses and zero-coefficient onset. This source explicitly
uses the original BULK energy, is queried once per energy trial, and holds its
coefficients fixed during the surface root. Its bounds must lie inside the
original caloric/conductivity support and it must retain the same caloric object.
No bulk dk/de term is added to the already differentiated integral.

Only the actual material path requires its conductivity source support; a bath
may lie outside it if the solved material face remains supported. Geometry is
the original positive planar area and half-distance, with an optional finite
nonnegative coating/contact resistance before the shared exterior surface.
That zero-storage resistance excludes convection and radiation. The supported
material face and exterior coating surface are distinct when resistance is
positive; conductivity support applies only to the actual material path. This remains a source-free quasi-steady
single-phase half-cell, a prescribed gas bath and a large black gray-radiation
surrounding, without finite enclosure storage or inferred view factors. Bulk
surface coefficients are not laws of a separate evolving surface temperature
or phase. The adapter allocates no trial tables and adds no source authority;
an optional source retains its own evaluation costs. Native physical-history
publication remains transactional on source failure or cancellation.

Four new native regressions pass: independent conductivity-integral/radiation
balances and Cp/Cv tangents; original coefficient feedback including zero onset;
actual backward-Euler evolution and neutral-body SDIRK2 bath entropy; and
source/geometry refusal with cancellation/replay from accepted positive time.
All 22 focused native original Cp/Cv and integrated exchange tests pass.

## Source-integrated Cp/Cv exchange with reservoirs and ambient baths

`IntegratedCaloricBoundary` reuses the prepared original Cp/Cv layer and the
existing `ReservoirLayerGeometry` material-to-bath solver. It integrates k(T)
through the planar half-cell, resolves the actual material-face temperature and
one authored film/contact resistance, and supplies the native boundary owner
with the complete outward power derivative with respect to original h or u.
Its historical h-labelled native fields do not change the storage basis.

The native reservoir or ambient channel supplies its current bath temperature,
area and energy/entropy accounts. The selected source-integrated face replaces
that channel's bulk-frozen coefficient; it does not add another authored htc,
heat source or face inventory. Original caloric-object identity and fixed
area/film are checked by the attachment path. The ordinary query also requires
the original distance; the explicit density-aware query below carries its
changing normal distance and complete tangent through the same bounded scalar
solve without material-card queries or allocation.

Material support applies to the material trace and solved material face. A
bath outside that support may exchange heat through a finite film while the
material remains supported; perfect contact refuses an unsupported prescribed
face. No source extrapolation or temperature clipping creates an endpoint.
The single-phase, fixed-context and quasi-steady half-cell assumptions of the
integrated layer remain explicit.

Three focused native regressions pass: independent affine-k film roots and
Cp/Cv energy derivatives in heating/cooling/equilibrium; a large actual
reservoir-plus-ambient step with independent caloric and bath entropy/energy
checks; and source-domain, geometry and late-cancellation refusals with exact
retry. Together with the integrated links and prior Cp/Cv conduction tests,
15 native tests pass. This is not full-workspace or experimental validation.

## Source-exact fixed-density heat capacity in implicit conduction

`thermal::transport::isochoric::conduction::IsochoricConductionLaw` binds the
original Cv internal-energy law to scalar conductivity from the same material
card and complete fixed-density source path. Conductivity must explicitly
support that density. The original cell mass/volume admission remains active.
The private scalar-conductivity preparation and interpolation kernel is shared
with `IsobaricConductionLaw`; its existing source identity convention is retained.

The native explicit and implicit conduction interfaces receive the original
temperature inverse, dT/du = 1/Cv, and dk/du = (dk/dT)/Cv. Historical h-labelled
native fields carry u on this declared `InternalEnergy` basis. Source knots
remain independent, runtime evaluation allocates nothing, and the existing
chain, graph and adaptive SDIRK2 owners retain energy/entropy and acceptance
responsibility. This introduces no density evolution, pressure-work closure,
phase change or permission to use internal energy as stream enthalpy.

Four G0/G1 tests cover exact constitutive tangents, source/density refusals,
nonlinear two-cell exchange against an independent heat root, and branched
adaptive SDIRK2 decay against an analytic network mode. All four pass in the
focused native build with the original full thermal owner and numerical solver
kernels; this is not a full workspace Cargo verification.

## Source-exact isobaric heat capacity in implicit conduction

`thermal::transport::isobaric::conduction::IsobaricConductionLaw` borrows the
original `IsobaricHeatCapacityLaw` and binds independently selected scalar
conductivity on that same card, pressure and complete typed source path. Both
sources retain every original breakpoint and receipt. Conductivity must remain
positive throughout the supported interval and explicitly support its pressure.

The adapter supplies the native `ConductivityLaw` and `ImplicitConductionLaw`
interfaces: original caloric inversion, exact dT/dh = 1/Cp(T), and analytic
dk/dh = (dk/dT)/Cp(T). Conductivity corners use the right derivative except at
the upper support endpoint. No resampled enthalpy chart, second time integrator
or runtime source query is introduced. Prepared knot lookup is logarithmic and
state evaluation allocates nothing. Native conduction and time-integration
owners retain their own physical residual, energy, entropy and acceptance gates.
Adaptive phase-error control remains unavailable because this single-phase
caloric law declares no phase fraction.

Four regressions specify separate Cp and conductivity knots, analytic tangent
checks against independent finite differences, source/card/pressure admission,
and native implicit two-cell exchange against an independently derived nonlinear
heat root with total-energy and entropy checks. All four pass in the same
focused native build. This is source-bound single-phase isobaric transport,
not a phase/EOS closure or experimental validation.

## Source-exact heat capacity in finite-rate crystallization

`thermal::transport::isobaric::crystallization::IsobaricCrystallizationPair`
binds two original phase Cp laws and a selected positive parent-minus-crystal
enthalpy datum at one fixed pressure. All sources retain one material card and
the same complete non-temperature typed context. Both Cp sources must support
the kinetic interval and reference temperature; the calorimetry reference may
lie outside that interval. The selected enthalpy gap must explicitly support
its reference temperature and pressure.

The common chart is h_parent(T) = integral_ref^T Cp_parent dT and
h_crystal(T) = integral_ref^T Cp_crystal dT - L_ref. Consequently the latent
enthalpy varies as L(T) = L_ref + integral_ref^T (Cp_parent - Cp_crystal) dT.
Preparation checks its positive represented value at every source breakpoint
and every interior quadratic minimum. Endpoint checks alone are insufficient
when the phase heat capacities differ. These are floating-point checks, not
an interval stability certificate.

The original finite-rate `CrystallizationLaw` and thermal crystallization
timeline own kinetics, mass-weighted total enthalpy, recalescence and acceptance.
No additional latent-heat source or numerical time integrator is introduced.
Each rebased branch implements `ThermalLaw` and `CrystallizationCaloric`; callers
must register the rebased parent itself on the timeline. Prepared source indexes
provide logarithmic forward and inverse caloric lookup, with direct-integral
fallback and exact clipped branch endpoints. Runtime caloric queries allocate
nothing and neither query nor resample the original sources.

Five G0/G1 regressions specify unequal-Cp latent enthalpy, arbitrary original
datum invariance, a negative interior gap despite positive endpoints, source
context/support refusals, indexed multi-knot branch round trips, and native
adiabatic recalescence with conserved total enthalpy and original-law binding.
The focused test and reduced-memory metadata builds were killed by the shared
memory limit in the genuine pinned asupersync dependency before fs-material;
these new tests have not been executed or compiler-validated in this session.
Cp and one enthalpy difference do not establish relative phase entropy, Gibbs
stability, density, pressure work, composition equivalence or a melting limit.
Those source/applicability declarations remain explicit on the original owners.

## Source-exact heating with free expansion

`thermal::transport::isobaric::expansion::IsobaricFreeExpansion` combines
the original single-phase Cp integral with the original engineering free-length
and reference-density sources at one maintained pressure. Preparation binds the
same material card, exact typed path and explicit pressure support; it never
treats Cp as Cv or infers pressure qualification from an extra query coordinate.

`advance_heat` resolves temperature through the existing analytic Cp inverse,
retains direct enthalpy and Cp/T entropy increments, and reports conserved mass,
original endpoint density/geometry, outward p Delta V and derived
Delta U = Delta H - p Delta V. It creates an immutable proposal. The specimen
is uniform and expands freely against a quasistatic load; constrained stress,
pressure evolution and mesh motion remain separate physical owners.

Original piecewise-affine principal stretches determine a cubic volume path.
Their polynomial differences are evaluated interval by interval with compensated
summation, preserving tiny expansion work when rounded endpoint densities or
volumes coincide. Preparation caches source breakpoints; a local heat step
visits only crossed intervals and performs no allocation, hashing or material
source query. Four focused regressions cover affine Cp, heating/cooling,
anisotropic source breaks, constant volume, tiny datum-independent increments,
and incompatible or pressure-unqualified source refusal.

### Original isobaric density in native integrated conduction

`IsobaricFreeExpansion` also implements the native `ThermalLaw`: temperature,
enthalpy basis and Cp/T entropy delegate to its SAME original Cp owner, while
`bulk_density_kg_m3` evaluates the retained reference density and original
free-volume ratio, rho=rho_ref/J(T). Its identity includes both original
sources. It never changes stored mass or reference mesh geometry by itself.

`isobaric::expansion::conduction::IsobaricExpansionConductionLaw` composes
that exact caloric object with its existing `IsobaricConductionLaw`. Object
identity, fixed pressure/composition and the original enthalpy datum are
retained. Preparation borrows the original tables; runtime uses the original
inverse and conductivity knots. The exact specific-volume tangent is
d(1/rho)/dh = J_T/(rho_ref*Cp(T)), with the CURRENT original Cp and the
original one-sided expansion-knot convention. Zero and negative volume
responses remain intact. No Cp-to-Cv substitution or second p dV heat occurs.

`CaloricConductionLayer::from_expanding_isobaric` carries this same source
into `IntegratedCaloricLink` and `IntegratedCaloricBoundary`. Their optional
geometry-aware native queries accept current planar half-distances and each
(dd/dh)/d, retaining original area, contact coverage, coating/film and source
ownership. The existing Kirchhoff integral and bounded surface solve remain
the flux owners. At the accepted link face,
dQ/dd_i = (dQ/dR'')/k(T_interface_i); the complete enthalpy tangent combines
this length response with original dT/dh and original contact dR''/dh.
The material-to-bath geometry term is
-Q*(d'/d)/(1+(A/d)*k(T_surface)*R_film), with R_film in K/W.
Each response uses the SAME accepted interface and original contact sample.

The native fixed-footprint reduction is a separate physical applicability
declaration. Free tangential expansion does not establish a constrained-solid
density law or fixed footprint; evolving area, 3-D strain compatibility,
stress and radial motion need their original physical owners. This extension
does not reinterpret equivalent radial resistance distances as physical
thicknesses. Cv and fixed-volume storage retain their existing constraints;
planar gray and surface-source integrated owners support the same optional
geometry tangent as described below.

Five regressions in `isobaric/expansion/tests/conduction.rs` cover current-Cp
volume tangents and source corners, complete integrated boundary/contact
finite differences with one original contact sample, an independently
eliminated heating/cooling energy equation with nonzero pressure work, and
closed graph SDIRK energy, entropy and inventory. The physical planar fixture
has normal-only free expansion and unchanged tangential stretches. These
regressions are authored and syntax checked; native execution remains blocked
by workspace storage and dependency-build resource limits.


### Expanding original Cp bodies with radiation and actual surface sources

`IntegratedCaloricGrayLink`, `IntegratedCaloricSurfaceLink` and
`IntegratedCaloricGrayBoundary` also implement the optional native planar
geometry query. They combine current density-driven half-distances with the
original h/Cp inverse and k(T) integrals. The native owner still supplies
d=m/(2*A*rho(h)) and (d'/d)=rho*d(1/rho)/dh; source heat capacities, pressure,
reference density, stored mass and reference mesh retain their original owners.

For a two-sided flux Q from A to B, the SAME accepted surface elimination gives
`Q_dA=-(Q/A)*Q_TA/k(TA)` and `Q_dB=(Q/A)*Q_TB/k(TB)`.
The endpoint partials hold source energies fixed; their density contributions
are added to original 1/Cp and source-energy terms exactly once. This includes
radiation and a separately authored nonradiative gap resistance, or a complete
original two-surface exchange law. Neither path queries the accepted source
again or repeats the nonlinear solve. Signed differential responses remain
available even when the heat/secant is passive.

For the exterior owner, C=A/d and f=1/(1+Rcoat*C*k(Tmaterial)). After eliminating
the coating, gc=f*C*k(Tbulk) and gs=f*C*k(Tmaterial). The complete tangent is
`Q_h=[F_T*(gc*T_h-Q*(d'/d)*f)+gs*F_h]/(gs+F_T)`, using the retained
Kirchhoff-integral Q and original accepted exterior partials. The existing
positive combined-stiffness admission remains authoritative for negative
exterior feedback. Bulk-only coefficients are sampled once, actual exterior
sources retain their accepted surface query, and both actual bath powers and
temperatures still reach the native entropy/heat accounts.

Area, coverage, coating resistance and independently declared gap/source
geometry remain fixed. Changing solid thickness does not infer gap motion,
new view factors, changing convection length or lateral strain constraints.
Fixed cylindrical/spherical surface-link coordinates retain their prior
support; the density-driven geometry hook refuses those radial coordinates.
There is no additional pressure-work heat: original isobaric storage already
accounts for Q=Delta H=Delta U+p*Delta V at the declared maintained pressure.

Six regressions in `isobaric/expansion/tests/conduction/gray.rs` cover both
heat-flow directions and original source-energy feedback against complete
finite differences, unchanged source-query counts, coating and signed exterior
feedback, native SDIRK closed-network mass/energy, and native backward-Euler
actual-bath entropy with independent enthalpy/pressure-work accounting. Final
Rust buffers are syntax checked; these new native regressions have not been
executed because workspace storage prevents a dependency build.

## Saturated-water boiling and condensation

`phase::vaporization::IsobaricWaterVaporization` binds the original
`fs-thermochem::water` IAPWS source to a fixed-pressure liquid/vapor enthalpy
plateau. It evolves actual vapor mass fraction, additive specific volume,
bulk density, internal energy and entropy; vapor is not represented as liquid.
The fixed total mass may expand against an external pressure load. The
returned pressure is the original source evaluation within the caller's
pressure-inversion tolerance, not a patched target pressure.

Finite heat increments independently audit Q - Delta H and
Q - Delta U - p Delta V against one explicit absolute joule budget.
Pressure work, internal energy and entropy changes use the actual represented
vapor-mass change and the original source's direct phase gaps. Exact requested
heat to either pure endpoint selects that source endpoint; excess heat still
refuses rather than extrapolating into subcooling or superheating. Critical
coalescence, unrepresentable increments and cancellation also refuse.

The existing shared-thermal `ThermalLaw` consumer uses T(h)=T_sat and the exact
plateau entropy integral Delta h/T_sat. Its geometry is reference heat-exchange
geometry, not a rigid-vessel constraint or a dynamically updated fluid mesh.
Real timeline tests cover heating/cooling with ambient loss, finite contact
between a steam condenser and a lower-pressure boiler, positive heat-transfer
entropy production, complete phase depletion/reentry and atomic refusal.
`cargo run -p fs-material --example water_boiling` runs a heating/cooling cycle
with phase inventory, volume, expansion work and independent energy checks.

### Adiabatic pressure-drop flashing

`WaterVaporizationState::flash_to_pressure` resolves the equilibrium endpoint
of a steady adiabatic throttle with no shaft work and negligible kinetic and
potential energy changes. Specific enthalpy is retained bitwise. The original
IAPWS source at the lower pressure determines actual liquid/vapor inventory,
volume, internal energy and entropy. Independent source Delta u + Delta(pv)
and signed entropy-generation checks use explicit caller budgets; p_out Delta v
is not substituted for the inlet/outlet flow-work difference.

The report retains the actual pressure drop, original pressure-inversion
residual, flow-work change and signed entropy generation. Pressure increases,
critical collapse and enthalpies requiring subcooling or superheating refuse;
no off-saturation source or valve flow-rate model is inferred. The returned
state's `law()` reuses its exact retained chart for further heating or cooling.
Tests include vapor creation from saturated 450 K liquid throttled to
101325 Pa, a heating/flash/condensation chain, conservation and source-domain,
accuracy and cancellation refusals.


## Original sources at native shear and electrical operating coordinates

`phase::transport::resolve_phase_transport_state_with_shear_rate` extends the
existing phase-first material resolution with one explicitly typed inverse-time
coordinate owned by the native flow. Temperature and phase still come from the
original enthalpy chart. Liquid viscosity and conductivity may depend on actual
shear; density and caloric storage may not acquire an unowned rate dependence.
Solid endpoints retain the existing omission of liquid-only properties.

`conductor::response::ScalarConductorResponse` retains one complete pinned scalar
conductor bundle and explicitly declares at most one electric-field magnitude
and one current-density magnitude coordinate. Native spatial consumers query
the original source at actual |E| and |J|, including all original requirements,
claim pins, coordinate quantities and external context. Conductivity and
resistivity remain distinct source bases. Every original electrical response
coordinate must be owned; frequency, strain and other unowned responses refuse.

`ScalarConductorResponse::relax_operating_point` prepares a nonlinear solver
trial by interpolating only those owned E/J magnitudes and querying the same
complete source again. Both inputs retain the original source/role/claim
identity and exactly equal non-electrical context. It never interpolates
resistivity or conductivity values. The intermediate query must independently
belong to the admitted support, including nonconvex multiaxis support; supported
endpoints alone do not authorize their connecting segment. A full step or a
source without operating coordinates requires no additional query. This helper
does not accept a physical state: the field owner must still evaluate the source
at actual E/J and satisfy its original current, conductivity and work gates.

`electrical_interface::response::ScalarElectricalInterfaceResponse` applies the
same original-source discipline to local interface voltage jump [V] and normal
current density [A/m2], oriented from the card's first surface to its second.
Signed and magnitude coordinates are distinct; neither bulk electric field nor
total face current is silently substituted. The complete ordered surface card,
property roles, selected claims and quantity conventions survive every refresh.
Temperature and explicitly declared side liquid-mass fractions remain owned by
their physical consumers. Positive finite R'' is required for this response
rung; perfect contact still requires an appropriate continuity constraint.
Optional relaxation requeries only the electrical coordinates and refuses
non-electrical context blending or unsupported source points. Uniformity,
spatial quadrature, current conservation, Joule heat and final constitutive
acceptance remain responsibilities of the native field owner.

The shared multiple-owned-coordinate predicate preserves the existing static
external-coordinate rules; it provides no evolution law itself. The electrical
response adapter is an isotropic, passive, memoryless rung, with no tensor, Hall,
hysteresis or field-storage claim. Ordinary phase metadata is retained; combined
phase-fraction/E/J source admission is not established merely by that metadata.
Actual field/current closure and finite-heat ownership belong to their original
native consumers, not this material adapter.

## Source-bound complete tensor relaxation

`visco_tensor::source::TensorPronyProperties` binds the original resolved
material state to a complete 3D Maxwell law. Equilibrium stiffness comes from
explicit isotropic parameter pairs, full orthotropic properties, or
`TensorElasticity::Mandel21` with caller-selected keys for all 21 independent
Mandel stiffness entries. The last family retains signed normal/shear and
shear/shear couplings in the original material axes; every zero is explicit.
It admits the complete positive-definite equilibrium tensor through the native
constitutive owner, without fitting an orthotropic or isotropic substitute. Each source
branch names all 21 symmetric upper-triangle Mandel entries and its relaxation
time; zero entries remain explicit data. A positive small-strain limit is also
a selected source claim. Shared names resolve once, incompatible quantity roles
refuse, and the original material state retains all selected usage receipts.

Every coefficient must be finite and constant at the admitted external state.
The law retains the intersection of all selected frequency domains, including
equilibrium, strain, time and branch data. The original tensor constructor owns
PSD admission and its roundoff policy. Missing data, frequency-sampled elastic
coefficients offered as equilibrium springs, invalid complete elasticity,
nonpositive times or strain limits, and more than 64 source branches refuse.
This binding does not condense plane stress, invent a full tensor from in-plane
coefficients, evolve temperature/moisture, or establish measured-material
authority. Full-equilibrium source regressions exercise coupled stress relaxation,
independently computed storage and heat, joint frame rotation, and affine thermal
free energy/entropy through the existing native history and projection owners.
Five original focused source tests cover full coefficients, shear convention,
orthotropy, shared names, source bands and refusals; native execution of these
new tests remains unverified in the resource-constrained authoring environment.

## Continuous volumetric storage and expansion-aware transport

`thermal::transport::volumetric::PreparedVolumetricThermalTransport` prepares
positive, continuously supported volumetric capacity C_vol(T) and directional
scalar k(T) on one fixed typed source context. It retains both original source
lines, selections, receipts, uncertainties and the immutable card identity.
Every source breakpoint enters the union grid. Point evaluation uses one binary
search, interpolates the independent inputs, then derives k/C_vol and effusivity;
no density, mass-specific heat or covariance is invented.

`storage_change_j_m3(T0,T1)` integrates the piecewise-linear C_vol chart directly
over crossed segments. `temperature_after_storage_change(T0,delta_storage)`
inverts that same monotone relation for heating and cooling using bounded knot
traversal and a scaled quadratic solve. Both use the caller's fixed volume and
caloric convention. They do not automatically represent internal energy,
enthalpy, pressure work, changing volume, latent heat or phase change. Numerical
exactness is for the selected affine model, subject to roundoff, not a measured
interpolation-error claim. Nonrepresentable updates and unsupported temperatures
refuse without mutation. A representable local inverse remains available when
the unused whole-interval integral overflows. Successful runtime calls allocate
nothing and perform no material queries.

Thermal path admission distinguishes an explicitly typed, fixed
`TemperatureDifference` test condition from a second absolute temperature.
Unknown or extra absolute-temperature response coordinates remain unsupported.

`thermal::transport::expanding::PreparedExpandingThermalTransport` composes the
original `PreparedFreeExpansionDensity` with selected Cp(T) and k(T), requiring
the same card, exact endpoints, axis and fixed context. It keeps the original
density model and its volume/reference assumptions intact. At runtime it
evaluates that model and the Cp/k grid, then derives transport; it never replaces
free-expansion density with a linear endpoint fit or a new measured claim.
Preparation checks density-only knots as well as thermal knots. Evaluation uses
at most two binary searches without allocation, hashing, queries or mutation.

Focused integration tests are `material_expanding_transport`,
`material_volumetric_curves` and `fsri_plastic_thermal_profiles` in xtask.
They cover analytic piecewise storage/inverse identities, source and numerical
refusals, original NIST expansion/Cp/k with an explicitly assumed reference
density, and real FSRI foam/plastic heating and cooling at fixed volume.

## Caloric source provenance survives the original phase-curve handoff

`EquilibriumEnthalpyPhaseCurve::source_receipts` retains the exact original
source usages from both `ResolvedHeatCapacityCurve` and
`ResolvedMeltingCurve`. Cp, density, transition temperature and latent-heat
selection/context remain available when a thermal consumer receives only the
plain curve. Existing wrapper `receipts()` methods delegate to this same
inventory. Source identity hashing and caloric interpolation are unchanged:
the receipts already contributed to the sourced chart identity.

This closes a loss of source ownership when another subsystem evolves a
coordinate held fixed while the caloric temperature line was constructed.
Consumers must inspect the actual selected claims before treating that chart as
independent of the evolving coordinate. Explicitly authored h(T,phase) charts
carry no resolved-source receipts and retain their declared law; an empty
receipt list is not empirical material validation. Tests retain Cp/density and
melting receipts through the original curve and its clone.

## IF97 water heating, boiling and superheat

`phase::water::IsobaricWater` gives one original-source, constant-pressure
enthalpy chart spanning compressed liquid, liquid/vapor coexistence and
superheated steam. It consumes only `fs-thermochem::water::if97`; original h, u,
s and v share that source's reference and are never mixed with SR1-86.
Temperature bounds must bracket saturation, remain within 273.16–1073.15 K and
admit the Region 1/2 coexistence join only through 623.15 K. This is not an ice,
Region 3, supercritical, metastable, kinetic phase-change or rigid-vessel model.

Original-source p,h and p,s inversions use safeguarded Newton/bisection, caller
residual and iteration limits, and cancellation. Coexistence resolves quality
analytically at the fixed saturation temperature. Actual source endpoints and
residuals are returned without shifting h or s. Single-phase states report no
vapor quality; mass fractions, density and specific volume remain available in
every phase. Temperature alone at saturation refuses as ambiguous.

Construction proves positive Cp over both complete sensible branches using the
source's bounded interval arithmetic. The retained global inverse slope admits
stable use of the existing `fs-thermal` timeline across latent plateaus; it is
not a sampled or locally guessed heat-capacity bound. Chart construction is
separate from repeated state queries.

The original IF97 regional fits have a finite saturation matching error:
delta_s=(h_v-h_l)/Ts-(s_v-s_l). Admission requires an explicit bound on its
magnitude. Raw source entropy and this defect remain visible. The thermal law
uses the continuous pressure-specific primitive S*=s+chi*delta_s, where chi is
zero in liquid, vapor quality in coexistence and one in vapor. It integrates
dh/T on this isobar, and must never be substituted for absolute source entropy
in a pressure-changing device. The correction is not physical entropy
generation or an added heat source.

Thermal entropy increments retain the supplied stored enthalpy difference.
The original-source ratio Delta S*/Delta h supplies the mean reciprocal
temperature between the admitted inverse endpoints; coincident temperatures
use Delta h/T. This preserves tiny sensible and latent heat contributions
when source inverse or quality rounding would erase or change them. Endpoint
accuracy still follows the declared inverse controls. Native regressions call
the shared `ThermalLaw` directly for liquid, steam, and coexistence at two
pressures, including reversed and single-ULP heat increments.

`advance_heat` retains fixed mass and reports both phase inventories, physical
volume change, delta_h, delta_u and outward p*delta_v work. Original-source
stable increments and actual rounded endpoint differences independently satisfy
the caller's unchanged absolute joule budget. Nonzero heat or phase inventory
lost to representability refuses. The shared thermal geometry remains reference
heat-exchange geometry, while physical expansion is explicitly reported.

`WaterState::isobaric_heat_capacity_response` distinguishes original finite
sensible Cp from an isothermal latent response. Strict finite-Cp queries
refuse coexistence and exact saturation endpoints; the latter explicitly
select the plateau-side generalized `dT/dh=0`. The borrowed
`phase::water::conduction::WaterConductionLaw` binds the exact registered chart
and a caller-declared bulk k(T) to the native implicit chain/graph. One original
enthalpy inversion supplies T and Cp, and analytic dk/dh uses that same state.
Original enthalpy, continuous dh/T entropy, conserved reference mass and vapor
mass fraction remain available to the thermal solver and adaptive controls.
The declared bulk k is quality-independent at saturation; it supplies no
inferred water/steam mixture conductivity, boiling coefficient, convection,
evolving geometry or rigid-vessel model. Native tests cover source Cp and
tangents, latent/omitted-sensible endpoints, source-domain refusals and
unequal-mass implicit heat exchange through boiling with an independently
reconstructed endpoint and conserved energy.

Six native material tests exercise p,h/p,s source round trips; complete
heating/condensation and phase depletion; explicit entropy-seam accounting;
domain/work/cancellation refusals; the real shared thermal history through all
three regimes; and finite steam/liquid thermal contact with conserved total
enthalpy and positive modeled entropy production. The runnable
`cargo run -p fs-material --example water_heating` heats and cools 50 g from
300 K through boiling to 700 K and back, independently reconciling every
accepted thermal interval with material calorimetry. This establishes a bounded
liquid/vapor consumer for bead `frankensim-7sga6.2.7`, not its complete ice and
all-regime scope.

## IF97 adiabatic pressure devices

`phase::water::flow` supplies steady throttles, liquid pumps, vapor compressors
and steam expanders on the same original IF97 charts. All quantities are per
kilogram of throughput. Pressure directions are explicit; kinetic/potential
changes, transient machine dynamics and flow-rate laws are outside this owner.

A throttle resolves h_out=h_in at a lower pressure, admitting compressed liquid,
wet or superheated outlets where the prepared source chart supports them.
A shaft device first resolves the original isentropic outlet from
s_out=s_in. Positive shaft work is outward: an expander gives
w=eta*(h_in-h_isentropic), while a compressor/pump gives
w=(h_in-h_isentropic)/eta. Efficiency is an explicit input in (0,1], with no
invented empirical default. The actual outlet follows from h_out=h_in-w.
For eta=1 the admitted isentropic source state is reused directly, avoiding an
unnecessary second inverse and cancellation in reconstructed target enthalpy.

Each result independently admits delta_h+w and delta_u+delta(p*v)+w against
the caller's unchanged numerical energy budget. Flow work is p_out*v_out minus
p_in*v_in, not p_out*delta_v, and is not added again to shaft work. Cross-pressure
entropy always uses original source s. Its signed change and each chart's
saturation matching defect are reported separately; the fixed-pressure thermal
entropy correction is never used as a global reference. Source-domain failures,
lost nonzero work or efficiency loss, entropy decrease outside the caller's
roundoff allowance, budget failures and cancellation publish no partial result.

Five focused native device tests pass alongside the six material tests. They
cover liquid/steam throttles, pump work against independent integration of
v*dp along a raw-source isentrope, superheated compression, wet turbine outlets,
a complete Rankine cycle, and physical/domain/cancellation refusals.
`cargo run -p fs-material --example water_rankine` exercises a 10 kPa
condenser, a 4 MPa boiler at 773.15 K, pump efficiency 0.80 and turbine efficiency
0.85. It reports states and all four process balances, including a separate
internal-energy/flow-work cycle check. Its approximately 31.22% efficiency is
the result for those illustrative inputs, not a calibrated plant prediction.
Heat-exchanger reservoir entropy and pressure losses are not supplied.

## Shared tensor strain admission in all material-history interfaces

`TensorProny::validate_mandel_strain` checks the original finite small-strain
bound in the Frobenius norm, with Mandel shear normalization. The original
Mandel `step`, tensor-shear `history().advance` and projected strain observer
use this one admission rule before constitutive history can change. Source
limits therefore survive a change of consumer or coordinate representation.
The focused regression checks shear-stress parity and identical refusal without
mutation in both time-domain interfaces. This admits supplied endpoints only;
it does not enclose strain extrema between time samples or quadrature points.

An assumed-strain projection may retain an independent compatible strain map
through `ProjectedTensorProny::with_compatible_strain_columns`. Both total and
elastic compatible strains must satisfy the same material domain as the work
map. Constitutive stiffness, Kelvin loss, Maxwell memory and conjugate forces
continue to use the original work map. Eigenstrain updates retain the additional
map, and material-memory continuation refuses a changed or dropped map.

## Volumetric storage on the shared thermal timeline

`thermal::transport::volumetric::timeline::FixedVolumeThermalLaw` borrows the
existing prepared volumetric source path and implements the original
`fs_thermal::enthalpy::ThermalLaw`. Its explicit reference mass/volume ratio
rho0 and zero-storage datum give h(T) = integral C_vol(T) dT / rho0. A compatible
cell therefore stores m*h = V0*integral C_vol(T) dT. rho0 is an authored scenario
basis, not a new measured density, temperature-dependent density or source Cp
claim. Mass and volume remain fixed. Callers must select a compatible caloric
constraint; this does not infer pressure work, phase change or moisture exchange.
Its `DeclaredSensibleStorage` basis preserves this distinction at runtime:
ordinary heat exchange remains available, while stream and phase-mass transfer
cannot treat an unspecified measured capacity as thermodynamic Cp or Cv.

Preparation retains the original source identity, includes the mass basis and
datum in the adapter identity, checks representable storage bounds and computes
an upward-rounded global temperature-slope bound from the minimum affine
capacity. Runtime storage/inversion reuse the original exact segment integrals.
`admit_cell` checks actual mass versus reference volume and supported storage;
its 64-epsilon comparison permits arithmetic roundoff only.

`FixedVolumeConductionLink::try_aligned` binds an existing ordered interface
to the two original borrowed caloric laws. Each evaluation verifies the active
link geometry, actual cell mass/volume bases and original material-law objects.
It evaluates each cell's sourced k at its current temperature, then uses the
existing series/contact conductance owner. Distances, transfer area, contact
resistance and alignment are explicit caller inputs. This is a uniform-cell,
cell-centre conductivity approximation, not an exact nonuniform half-cell
integral. Return the value through the existing implicit
`StepCoefficients::conductances_w_k`; the original timeline owns history,
sources, integration and local/global energy acceptance. Evaluation allocates
no source tables, queries no card and computes no content hashes.

Four focused G1/G3 tests in `xtask/tests/material_volumetric_timeline.rs`
cover reference-basis/datum invariance, an independently solved nonlinear
two-cell balance, original-law/geometry/domain refusals, and actual compiled
FSRI foam and density-free FRP heating on the existing timeline. The explicit
authored mass basis cancels from those fixed-volume temperature results.
An entropy interpretation is available only through the explicit opt-in below.
No experimental qualification or full spatial FEM result is added.

## Sensible entropy from the retained volumetric path (2026-10-28)

`PreparedVolumetricThermalTransport::sensible_entropy_change_j_m3_k`
integrates the retained piecewise-affine capacity as
`integral C_vol(T)/T dT` in J/(m3 K). Calling it explicitly selects a
reversible sensible-caloric interpretation at fixed work/composition
coordinates and the source reference-volume convention. An apparent
experimental capacity does not automatically become thermodynamic Cv.
No absolute entropy datum, latent heat, pressure work, chemistry, moisture
exchange or hysteresis is supplied.

Each crossed segment uses positive endpoint weights for its analytic
logarithmic integral. A bounded small-interval series avoids cancellation;
compensated accumulation retains the original source breaks. Reversal is
antisymmetric. Unsupported temperatures and unrepresentable arithmetic
refuse without changing the source path. Runtime traversal is bounded by
the original knots, allocates nothing, queries no card and hashes nothing.
The platform logarithm remains subject to floating-point roundoff; no
cross-ISA transcendental or measurement-uncertainty guarantee is added.

`FixedVolumeThermalLaw::with_sensible_entropy` opts the SAME borrowed source
law into the existing `ThermalLaw::specific_entropy_change_j_kg_k` account.
Default bindings still report entropy unavailable. Opt-in changes the
binding identity exactly once, retaining the source receipts, authored
reference mass/volume ratio and storage datum. The caller must justify
the capacity's caloric/work constraint. In exact arithmetic,
`m * delta_s = V0 * integral C_vol(T)/T dT`; the authored reference mass basis
cancels.

The adapter uses the actual stored enthalpy increment and a capacity-weighted
mean reciprocal temperature. Distinct enthalpy coordinates that round to
one temperature retain their nonzero `delta_h/T` increment. Existing
storage/inversion admission remains active: even finite endpoint
coordinates may refuse if their difference or another required intermediate
is unrepresentable. No broad all-finite-input success guarantee is made.

The original `fs_thermal` explicit and backward-Euler owners continue to
own heat flux, state acceptance, reservoir entropy, conductive generation,
unclassified source watts and temporal defects. The adapter neither
reclassifies a source as boundary heat nor counts numerical entropy as
material production. Three additional tests in
`xtask/tests/material_volumetric_timeline.rs` cover independent analytic
integrals, distinct source breaks, adjacent-temperature/enthalpy increments,
reference-basis/datum invariance, default refusal and opt-in identity, an
isolated implicit conduction balance, and actual compiled FSRI PP/mattress
foam heating. All seven timeline tests pass, including the prior four
energy/storage regressions.

## Integrated conductivity layers on the volumetric timeline (2026-10-29)

`FixedVolumeConductionLink::into_integrated_layers` explicitly selects the
existing quasi-steady, one-dimensional layer model. The returned
`FixedVolumeLayerConduction` copies each original conductivity line once,
preserving its knots and temperature support. The borrowed caloric paths retain
their source receipts, fixed reference mass basis, datum and entropy selection.
It integrates k(T) through both legs using the shared `TwoLayerConduction`
kernel; the kernel's public API and physical equations remain unchanged.

Every evaluation admits the actual active link, ordered endpoints, transfer
area, original material-law objects, mass/volume bases and stored-energy
coordinates. Each conductivity chart covers only its own bulk-to-face
temperature interval. Supports may overlap only at a perfect-contact face or
be disjoint across a finite contact jump; neither source is queried at the
opposite material's bulk temperature. The shared flux is bracketed by each
original chart's available Kirchhoff integral, and exact chart inverses plus
a safeguarded monotone solve close the contact jump. An unclosed jump at the
largest supported flux refuses the actual interface as outside material
support. Independent layer-heat and contact-jump gates remain authoritative.
Positive lengths and area, nonnegative area-specific contact
resistance, an iteration budget and numerical tolerances are explicit inputs.
The caller declares homogeneous, internally source-free legs and alignment of
the selected scalar k with heat flow. This does not establish isotropy or source
qualification for unmeasured conditions.

The response carries one signed heat transfer, a secant conductance, both
contact temperatures, endpoint derivatives and continuity diagnostics. Feed
its conductance into the original implicit `StepCoefficients` callback on every
candidate state. The secant conductance is distinct from the returned endpoint
Newton derivatives. The existing timeline owns storage, acceptance, sources,
energy conservation and optional sensible entropy. No layer storage or second
time integration is introduced. Successful runtime evaluation allocates no
memory, performs no database query and computes no content hash; work remains
bounded by the retained knots and the existing solver's iteration budget.

Three additional tests in `xtask/tests/material_volumetric_timeline.rs` cover
independent affine-conductivity integrals and derivatives, contact between two
different materials, source/geometry/domain refusals, cancellation and budget
exhaustion, and compiled FSRI polypropylene on the original implicit energy and
entropy account. All ten tests in that file pass in the focused native build.
These are bounded constitutive and timeline results; spatial convergence,
experimental qualification, phase change and dynamic moisture transport remain
outside this model.


## Exact equal-time material-history coordinates

`visco_tensor::coordinates::RelaxationCoordinates` groups displacement-driven
history rows only by bit-identical positive relaxation time. For each original
factor `A`, the existing Householder QR owner supplies `A = Q R`; complete thin
columns, including rank-deficient and zero directions, remain. With relaxed
`z = A q` or unrelaxed `z = 0` initialization, `z = Q y` is invariant and
`y' = (R q - y) / tau` preserves the original force, storage and dissipated
heat. Nearby poles are not merged and no mechanical mode or material coefficient
is fitted. Families no larger than the displacement space keep original rows.

Column scaling, orthogonality and source-entry reconstruction checks refuse
unrepresentable transforms. Explicit state/work budgets and cancellation cause
refusal rather than truncation. Immutable `RelaxationMemoryMap` slices share one
complete inverse map and recover every original material-point history for its
existing constitutive and source-domain checks. Independently authored original
histories outside this invariant subspace cannot be silently projected into it.
This is a coordinate representation of the declared original material dynamics,
not permission to lower spatial, modal or temporal convergence requirements.

## Joint phase and electrical source coordinates

`ScalarConductorResponse::new_with_liquid_mass_fraction_axis` binds the original enthalpy-derived liquid mass fraction and explicitly owned actual E/J magnitudes together. Complete source requirements, selected claims, typed coordinates and support survive phase, temperature and field refresh. Public static phase admission and generic material refresh remain strict. The existing source data supplies its admitted law; no universal phase mixture rule or higher-dimensional table is invented. Focused source regressions cover retained roles/pins, actual-field queries, phase support and unowned-coordinate refusal.

## Free dimensional change without a density requirement (2026-10-29)

`thermal::transport::expansion::PreparedFreeExpansion` consumes original finite
engineering length-change sources with an explicit reference temperature. The
new `FreeExpansionGeometryQuery` selects either an authored isotropic
interpretation or three fixed orthogonal principal directions. No density or
mass input is needed to compute the supported dimensional change.

At temperature T, each principal stretch is
`(1 + epsilon_i(T)) / (1 + epsilon_i(T_ref))`; the relative volume is their
product. Source breakpoints, nonzero fitted datum residuals, complete typed
query context and original receipts remain intact. The reference is retained as
an exact union knot. Only the independent source strains interpolate; stretches
and volume are derived from them. The identity includes the reference, frame,
original card and source receipts. Evaluation uses one binary search and
allocates nothing, queries no card and computes no hash.

`PreparedFreeExpansionDensity::try_from_expansion` optionally attaches a density
answer at that same reference and computes `rho_ref / J` under the caller's
constant-mass assumption. It requires the original immutable card and full
reference context, and reuses the prepared strain lines. The existing density
constructor, public sample accessors and v1 density identity byte order are
preserved. Existing expanding thermal transport continues to consume this
composition without separately interpolating density products.

Four additional tests in `xtask/tests/material_expansion_density.rs` exercise
independent anisotropic geometry, reference-dependent identities, optional
density composition, source/domain/frame/length refusals, and compiled NIST 304
dimensional change without an inserted density claim. All eight tests in that
suite and three existing expanding-transport tests pass in the focused native
build. Thermal coefficients retain their distinct dimensions and are refused
by this finite-strain API; mean and differential coefficients are not inferred
from one another. Constrained mechanics, changing directions, phase transitions
and uncertainty propagation remain separate model requirements.

## Rigid water vessels with evolving pressure

`phase::water::isochoric::IsochoricWater` closes equilibrium water at fixed
specific volume with internal energy as the storage coordinate. Original IF97
Regions 1, 2 and 4 supply compressed liquid, liquid/vapor coexistence and
superheated vapor over a declared 273.16–623.15 K window, with declared pressure
bounds through 100 MPa. Saturation quality follows the actual phase volumes.
Off saturation, a safeguarded pressure inverse uses the original analytic
volume derivative. A second safeguarded inverse resolves u(v,T); its derivative
is the actual isochoric heat capacity, including phase conversion and the
original Region 4 pressure derivative. It does not substitute the approximate
Clapeyron identity for that independently fitted source equation.

Every heat step keeps mass and vessel volume fixed, admits Q=delta_U and
Q=delta_H-V*delta_p independently, and reports zero boundary work. The pressure
part of enthalpy is not an additional heat source or wall work. Original-source
endpoint volumes, volume residuals, phase inventories, entropy and saturation
Gibbs mismatch remain visible; caloric values are not patched to fit budgets.
Nonzero heat or phase mass lost to representation refuses. Controls bound
volume/energy residuals, pressure/temperature iterations and source windows;
cancellation or refusal leaves the accepted state unchanged.

The temperature window is a bounded inverse search interval, not a certified
global Cv enclosure. Unsupported intermediate source queries refuse. This
owner does not claim Region 3, ice, mass flow, wall compliance, phase kinetics
or a spatially resolved vessel.

Six focused native tests pass: original p/v/u round trips; analytic two-phase
Cv against independent finite differences; sealed evaporation/condensation;
liquid pressurization and the cold-water density anomaly; heat subdivision and
small increments; domain, work, energy and cancellation refusals. The runnable
`cargo run -p fs-material --example water_rigid_vessel` heats 50 g in 5 L from
300 to 600 K and back. Pressure rises from approximately 3.54 kPa to 2.61 MPa,
with maximum complete-history energy discrepancy 5.27e-8 J and zero wall work.
These are numerical balance checks, not source-fit uncertainty estimates.
This advances the bounded phase-transition consumer in
`frankensim-7sga6.2.7`.

## Dense-water rigid-vessel calorimetry

`phase::water::dense::DenseIsochore` connects the original IF97 Region 3 source
to fixed-mass, fixed-volume heat evolution. The constructor retains one
`DenseIsochoricInterval`: original regional admission, stable derivatives and
positive Cv are established over the entire requested temperature path, with
explicit work limits and cancellation. Its enclosure and dT/du bound are reused
across state queries and heat steps; no proof is repeated inside the caloric
inverse.

A safeguarded Newton/bisection inverse uses original u(rho,T) and Cv to resolve
a heat-driven endpoint. Pressure follows the same Helmholtz equation. Original
stable increments separately provide delta_u, delta_h, delta_p and delta_s.
Accepted heat must satisfy both Q=delta_U and Q=delta_H-V*delta_p using those
increments and independently using the actual rounded source endpoints.
V*delta_p is the pressure contribution to enthalpy, while rigid-wall work is
exactly zero. No heat, entropy, caloric reference or nonzero increment is silently
clipped or reconciled. Unrepresentable heat, failed energy admission, unsupported
source states, work exhaustion and cancellation return no partial state.

This is a single-phase Region 3 equilibrium owner. It does not automatically
stitch across the Region 1/2/3 boundaries, provide high-temperature saturated
pairs, resolve a moving interface, or model vessel compliance and strength.
The separate `isochoric` owner retains the lower-temperature liquid/vapor
coexistence model. Neither substitutes its internal-energy coordinate for the
existing enthalpy timeline's declared caloric convention.

Four focused native tests pass: source u/T round trips and retained Cv bounds;
high-pressure heating/cooling with separate pressure enthalpy and zero wall
work; independent Cv/Cv-over-T integration plus subdivision and small heat;
and domain, work, energy, chart-identity and late-cancellation refusals.
The complete original water module compiles with both new volume owners, and
both vessel examples plus the existing Rankine example execute successfully.

`cargo run -p fs-material --example water_dense_vessel` heats 0.5 kg in 1 L
from 650 to approximately 750 K and cools it back. Pressure rises from about
25.584 to 78.310 MPa. Heating supplies 144.903 kJ to internal energy; enthalpy
increases by 197.629 kJ, including 52.726 kJ of V*delta_p and zero wall work.
The whole-path preparation uses 113 interval evaluations for its 649–751 K
source window, reused over 40 heat steps. Maximum complete-history numerical
energy discrepancy is 1.47e-6 J. Source uncertainty, reservoir entropy and a
container-strength prediction are outside that numerical balance claim.

## Compositional sensible heat storage (2026-10-29)

`thermal::transport::storage::PreparedHeatStorage` prepares one or more explicit
constituents at a shared temperature. `HeatStorageAmount::MassKg` multiplies a
source specific capacity; `ReferenceVolumeM3` multiplies native volumetric
capacity. Neither requires conductivity or an inferred density. Each component
retains its immutable card, exact claim selection, typed fixed conditions,
original source line, uncertainty and physical amount. The caller declares
compatible caloric/work constraints and the shared-temperature approximation.

Preparation intersects the selected temperature supports and retains every
source breakpoint inside that intersection. Constituent IDs must be unique;
sorting them before compensated summation and hashing makes input permutation
irrelevant. Capacity amounts must be positive, finite and dimensionally matched.
Preparation admits at most 64 constituents and 262,144 union knots. Runtime
capacity needs one binary search; energy, inverse temperature and optional
sensible entropy use the shared positive-affine caloric kernel. Successful
evaluation allocates nothing and performs no source query, content hashing or
numerical root iteration. Integrals traverse relative temperature intervals to
retain small changes; finite representability limits can still refuse.

Local caloric inversions now use one logarithmic knot lookup and constant
segment work. They no longer integrate every unvisited source interval before
inverting a small heat increment. Cross-knot updates retain the original
compensated whole-support admission and analytic quadratic inverse, including
the canonical cooling endpoint. The same kernel serves specific isobaric,
volumetric, tabulated and assembled capacities. G1 regressions check signed
nonlinear calorimetry, source-boundary refusal and source-query counts on
64- and 4096-knot curves; this is a work bound, not a general throughput claim.

The fixed-datum timeline laws (`IsobaricHeatCapacityLaw`,
`IsochoricHeatCapacityLaw`, `FixedVolumeThermalLaw` and `HeatStorageThermalLaw`)
prepare a compensated caloric search index in O(knots) time and storage. Repeated h-to-T
queries locate the source segment logarithmically and invert its original
affine-capacity integral analytically, including decreasing capacity. They do
not approximate the energy curve by straight lines or resample its source.
Forward T-to-energy queries reuse this same index: an intermediate source knot
between datum and target supplies a same-sign prefix plus the original local
affine-capacity integral. Cp/Cv changes to or from their datum also use this path,
including the existing convection tangent when its reference equals that datum.
There is no subtraction of large cached primitives. Other relative heat and
entropy still integrate the crossed intervals directly. Exact support endpoint
values retain the original accumulation and mass/volume conversion; a rounded
candidate outside that canonical energy span retries the original integral.
Numerically indistinguishable search knots or roundoff at a source boundary
retain direct evaluation and its original refusal semantics. Current operation
identities are Cp v3, Cv v2 and both generic storage wrappers v4; source identities
and uncertainty are unchanged. G1 tests cover both signs, interior datums, exact
support endpoints, near-datum increments and at most 16 capacity reads for
ordinary inverse queries and eight for forward queries on 64-, 4096- and
262144-knot tables. This is an algorithmic work test, not a hardware throughput
or experimental-accuracy claim.

`HeatStorageThermalLaw` binds this total storage to the existing thermal
timeline's specific-energy coordinate with an explicit total mass and datum.
Its basis is likewise `DeclaredSensibleStorage`, including after a sensible-
entropy opt-in. Both generic capacity wrappers retain this physical
interpretation in their identities; source identities and capacity curves are
unchanged.
Mass-only assemblies must match their known total mass; mixed assemblies require
additional positive mass for their volume-based constituents. The native
timeline invokes each law's `validate_cell` before birth or inventory change;
`admit_cell` also permits an explicit precheck. The timeline owns power, evolution, cancellation
and energy acceptance. Sensible entropy is opt-in and changes the law identity.

`storage::conduction::HeatStorageConductionLaw` connects this same amount-bound
caloric owner to the native implicit chain and graph solvers. The caller selects
an effective scalar conductivity source for the assembly; preparation retains
its original knots, fixed conditions and card identity across the entire storage
interval. Conductivity is not inferred from constituent capacities. Same-state
tangents use `dT/de = m_ref/C_total(T)` and `dk/de = (dk/dT)*dT/de`, preserving both
mass- and volume-based contributions, fixed reference mass and declared sensible
storage semantics. Focused public regressions cover original source tangents,
support and inventory refusals, and a branched nonlinear heat balance crossing a
conductivity breakpoint, with replay and cancellation. These new regressions are
authored but not executed in the current resource-constrained workspace.

Four focused assembly tests exercise independent mass/volume integrals,
inversion, entropy, amount scaling, source choice, permutation, domain refusals
and actual thermal timelines. Compiled MDF plus acrylic/mineral countertop works
below MDF's conductivity support; beryllium heats from its Cp source alone.
The three existing volumetric-curve tests and ten existing volumetric-timeline
tests also pass with the extracted shared caloric kernel. This addition supplies
no conductivity homogenization, phase or chemical mixing law, thermal-expansion
work, moisture evolution, absolute entropy, propagated uncertainty or automatic
justification of a uniform-temperature body.

## Dense isobaric heating, boiling and supercritical water

`phase::water::pressure::DenseIsobaricWater` binds the original Region 3
pressure/density solver to a fixed-pressure enthalpy owner. Explicit temperature
windows stay inside Region 3. Subcritical windows may retain the original
high-temperature Region 3/4 saturation pair; supercritical states have no
liquid/vapor quality or invented phase inventories. Original h, u, s and phase
volumes remain unshifted. Vapor quality is a mass fraction, and mixture volume
uses the constituent specific volumes.

Temperature, specific-enthalpy and raw-specific-entropy queries share the
retained chart. Caloric inverses use safeguarded Newton/bisection with original
Cp or Cp/T, explicit residual tolerances, source work counts, and 1..=128
outer iterations. Density inversion has its own bounded controls. Authored
caloric endpoints and latent endpoints are retained exactly. This search
window is not a global Cp certificate: unresolved interior source queries,
including near-critical resolution limits, can refuse. No `ThermalLaw`
stability bound or automatic Region 1/2/3 temperature-path stitching is claimed.

Fixed-mass heat deposition evolves H, with the independent physical balance
Q=Delta U+p_requested*Delta V. The requested mechanical pressure and actual
Helmholtz endpoint pressures remain distinct. Reports expose source pressure
residuals, endpoint h-u-p_requested*v defects, expansion work, actual phase
masses and volumes, and raw source entropy. The isobar-local entropy primitive
uses the explicit saturation defect delta_s=(h_v-h_l)/Ts-(s_v-s_l);
its correction is not an absolute entropy for pressure-changing equipment.

Stable bivariate source differences and latent lever-rule changes avoid
subtracting large reference energies for small heating steps. Acceptance gates
both stable and independently represented endpoint versions of Delta H-Q and
Delta U+p*Delta V-Q. Nonzero heat or phase mass lost to representation,
unsupported states, chart mismatch, exhausted work and early/late cancellation
return no partial state or source debit. A numerical residual budget is not a
physical source uncertainty or a license to replace source pressure.

Six focused native tests pass: 25 MPa dense-liquid/supercritical h and raw-s
round trips; 20 MPa high-temperature boiling and condensation with distinct
mass/volume fractions and one latent account; independent Cp and Cp/T integration
through the pseudocritical peak; subdivided heating/cooling, expansion work and
small/zero heat; coarse pressure closure rejected by the physical joule budget;
and domain, ambiguity, work, seam, chart, representability and cancellation
refusals. No transport coefficient, spatial boiling, moving interface, vessel
mechanics or reservoir entropy model is supplied by this equilibrium owner.

## Pressure equipment across original water source regions

`phase::water::flow::across` connects existing Region 1/2/4 Gibbs charts and the
new Region 3 Helmholtz pressure chart. Typed state/chart enums retain the exact
source owner and reuse the existing pressure-device controls and balance
report. Steady adiabatic throttles preserve original h; compressors and expanders
resolve original s at the outlet pressure, then apply the caller's explicit
isentropic efficiency. No pressure-local entropy correction is used across
pressures, and no caloric datum is shifted to manufacture regional agreement.

Every accepted device independently closes Delta h+W and
Delta u+Delta(p*v)+W, with outward shaft work positive and flow work included
once. Original h-u-p_requested*v defects at the inlet, outlet and isentropic
reference are individually limited by the energy budget so pressure errors
cannot cancel invisibly. Reports retain source entropy production, numerical
entropy residual, source-junction defects, density-source evaluation counts and
caloric iterations. Reversible devices reuse the isentropic endpoint. Numerical,
chart-domain, direction, efficiency, work and cancellation refusals return no
partial device state. The production dispatch and solves allocate nothing.

Four focused tests pass: a dense-liquid throttle reaching wet low-pressure
steam; supercritical expansion with raw entropy and declared efficiency;
hot-liquid compression across the Gibbs/dense interface with reversible return;
and invalid direction/efficiency, pressure closure, unsupported outlet, work and
cancellation refusals. The combined production water library also passes the six
new dense-calorimetry tests and all five existing pressure-device/Rankine tests:
15 integrated native tests in total.

`cargo run -p fs-material --example water_high_pressure_steam` runs a partial
steady steam train: 20 MPa water heats from 630 to 645 K through its
638.895911546 K saturation plateau, expands through an explicitly illustrative
85% efficient turbine to 10 kPa, and condenses to saturated liquid. Boiler heat
is 851.285931773 kJ/kg, including 584.286587172 kJ/kg of latent enthalpy.
Shaft output is 788.761576280 kJ/kg; the wet turbine outlet has vapor quality
0.659460798. Maximum independent stage energy residual is 1.08e-8 J/kg and
the complete open-train residual is 2.10e-9 J/kg on the native validation run.
These are numerical balances, not source or equipment measurement accuracy.

The inlet already contains substantial thermal energy. Feed-pump/preheater
return and reservoir temperatures are not modeled, so this example supplies
neither a closed-cycle efficiency nor a reservoir-entropy claim. Equipment
endpoint connections do not imply automatic temperature-path stitching of
Region 1/2/3, transient machine dynamics, pressure-drop/flow-rate laws, phase
kinetics, or resolved two-phase flow.

## Ordered material layer stacks (2026-10-29)

`thermal::transport::stack::PreparedLayerStackConduction` binds 1–64 physical
layers to their original scalar conductivity sources, thicknesses, common area
and explicit area-specific contact resistances. Each layer retains its stable
ID, immutable card identity, complete typed source context, selected receipts
and uncertainty. Physical order is preserved and contributes to identity;
repeated materials may have distinct layer IDs. No common-grid resampling or
effective-conductivity claim replaces the individual source curves.

The shared `thermal::conduction::LayerStackConduction` kernel solves one
conserved flux with exact piecewise-affine conductivity integrals. A signed
analytic Kirchhoff inverse traverses original breakpoints; safeguarded
Newton/bisection closes the final layer without a dense interface solve.
Explicit contact jumps determine the layer-face temperatures. Each original
chart covers only its own physical layer, so finite contact jumps may bridge
disjoint source intervals. A trial above a local source interval increases the
flux bracket; a trial below it decreases the bracket, without extrapolation.
Source-local Kirchhoff integrals bound the admissible flux. Two-layer stacks
delegate to the same `TwoLayerConduction` owner and its physical residual gates.

`PreparedLayerStackConduction::common_temperature_span_k()` returns the optional
intersection of the original source intervals, retaining an exact touching
point as `Some([T, T])` and disjoint sources as `None`. This describes supported
uniform temperatures, not the validity of a nonuniform steady wall. The
finite-cell transient constructor still requires a positive-width common
conductivity/storage interval; no broader caloric support is inferred.

The response exposes all ordered face temperatures, signed power, secant
conductance, both physical endpoint derivatives, and separate layer-heat and
contact-temperature residuals. Each layer's reconstructed power must meet the
requested heat tolerance. Contact representation allowances remain distinct
from that physical acceptance test. Equal-temperature and single-layer cases
have direct solutions. Successful evaluation allocates nothing, performs no
source queries or hashes, and runs checkpoints per layer and iteration under
the caller's finite iteration budget. Cancelled or failed work publishes no
response and mutates no prepared input.

### Retained material drops for thin conducting films

The two-layer and ordered-stack heat audits retain the relative temperature
displacement solved inside each original affine conductivity segment. They
reintegrate that solved displacement, including any fully crossed source
segments, before rounding the reported absolute face temperature. Rebuilding
a thin conducting film's heat from two rounded absolute temperatures can
amplify Kelvin roundoff by `k/d` and falsely reject a conservative physical
solution. The retained integral is recomputed from the solved displacement;
it is not copied from the requested common flux. Original heat tolerances,
source endpoints, contact-temperature gates and refusal of nonzero inverses
whose absolute temperature cannot change remain unchanged.

The same two-layer owner supplies bulk/contact, gray-gap and source-backed
surface-exchange closures. Ordered stacks retain these integrals for marched
layers in fixed bounded workspace; the final layer still closes the original
shooting residual against its exact cold boundary. No added source traversal,
allocation, effective conductivity or enlarged source domain is introduced.

Feed the secant conductance into the existing thermal owner's candidate-state
coefficient callback; endpoint derivatives serve a separate Newton assembly.
The original timeline continues to own storage, time integration, sources,
energy acceptance and optional sensible entropy. The caller declares fixed
planar geometry, source-free interiors and alignment of scalar conductivity
with heat flow. Unknown source orientation remains unknown. No internal layer
storage, radiation, lateral transport, phase law or inferred contact law is
supplied by this steady constitutive relation.

Four focused tests in `xtask/tests/material_layer_stack.rs` pass against the
production modules: signed inversion and analytic resistance, independent
nonlinear four-layer interfaces and derivatives, refusal/cancellation/budget
behavior, and actual compiled FSRI gypsum/mineral-wool/plywood sources. The
source integration independently constructs 0.1 W conduction and exercises
50 J transfer on the original thermal timeline with energy and sensible-entropy
closure. Existing two-layer parity and the 64-layer limit are included.

## Regional liquid/dense-water/steam caloric paths

`phase::water::regional::RegionalWater` prepares one fixed-pressure original
IF97 path across Regions 1, 3, 2 and 5, including Region 4 liquid/vapor coexistence
when present. Windows through 1073.15 K retain the original 0<p<=100 MPa domain;
windows extending above that temperature admit the original high-temperature
steam equation through 2273.15 K only at 0<p<=50 MPa. The exact critical singular
state remains unsupported. Only source pieces, saturation and junctions
intersecting the requested window are prepared. Cold compressed feedwater, high-temperature
boiling and superheated or supercritical steam can consequently share one
caloric owner. Low-pressure vapor-only windows and exact saturation endpoint
windows retain their original source admission.

Preparation retains at most five sensible pieces and four physical/source events
in bounded fixed arrays, with no heap allocation or coefficient/reference changes. Each single-phase
regional junction carries the original signed h, u, s and v differences and
temperature/pressure residuals. Callers separately budget enthalpy, entropy,
volume and boundary-temperature matching. A source matching defect is not
latent heat, a new phase, an entropy-production model, or a global equation of
state reconciliation. Original IF97 fits can produce positive gaps or negative
overlaps in h and s independently.

Temperature, h and raw-s queries share the prepared path. Caloric inverses
use original Cp or Cp/T in safeguarded Newton/bisection, bounded density
solves and a total outer-iteration budget shared across candidate pieces.
A generic caloric query refuses an unrepresented gap or multiple original
roots. Fixed-mass heat continuation retains the initial state and chooses the
first admissible root in the heating/cooling direction; it still refuses a
gap rather than inventing storage at a regional boundary. Exact zero heat
retains the state without an inverse. Numerical work and residuals remain
observable. Gibbs evaluations contribute outer iterations; the reported
source-evaluation counter counts the nested dense pressure/stiffness work.

Heat deposition evolves H and independently audits Q=Delta U+p_requested*Delta V.
Stable original-source differences are accumulated separately as sensible,
physical latent, and regional matching contributions in both enthalpy and
entropy. Reports retain phase mass and volume, endpoint pressure-closure
defects, stable differences and independent represented-endpoint residuals.
Supercritical states have no vapor quality or fictitious phase inventory, even
when their source routing reaches Region 2. Invalid budgets, exhausted work,
unrepresentable nonzero heat or phase mass, unsupported/ambiguous states,
chart mismatch and early/late cancellation produce no partial result.

Focused regressions cover heating and cooling from cold feed through dense-water
coexistence and all three source junctions into high-temperature steam;
independent h/s gap and overlap behavior, including different signs at the
Region 2/5 boundary; original Cp and Cp/T quadrature;
window-scoped preparation, low-pressure vapor and exact latent endpoints;
subdivision, tiny heat/quality and zero-heat identity; and source, matching,
work, energy and cancellation refusal. Published Region 5 states exercise both
caloric inverses, and existing throttle/shaft/cycle consumers retain the same
original h, raw s and pressure-flow-work convention at high temperature. Region 5
retains the pure-water source assumption without a chemical dissociation model.
This equilibrium pressure owner does not supply a whole-path Cp stability
certificate for `ThermalLaw`, ice,
resolved phase kinetics, spatial two-phase flow, or vessel mechanics.

## Complete high-pressure water cycles through the regional chart

`phase::water::flow::across` accepts `RegionalWater` and
`RegionalWaterState` through additional typed enum variants. The existing
throttle, compressor and expander equations retain original h and raw s,
reuse the regional inverses, propagate typed domain/gap/ambiguity refusal,
and keep the established independent enthalpy and internal-energy/flow-work
audits. Regional preparation and matching budgets remain owned by the chart.
The dispatch introduces no heap allocation, caloric datum shift, duplicate
flow work or alternative device efficiency convention.

`cargo run -p fs-material --example water_regional_rankine` runs three complete
four-device Rankine cycles: a 10 kPa saturated-liquid condenser outlet enters
an explicitly 80% efficient feed pump, a single regional boiler path heats to
823.15 K at 20 MPa or 25 MPa, or to 1500 K at 30 MPa, and an explicitly 85%
efficient turbine returns to
10 kPa, and the condenser returns exactly to the original saturated-liquid
state. Pump discharge is cold Region 1 water; each boiler traverses Region 3
and reaches Region 2 or Region 5 steam. The 20 MPa boiler includes physical
vaporization; the 25/30 MPa boilers have zero latent heat and no invented
supercritical quality.

The example reports every stage, all traversed original regional junctions, sensible/
latent/matching heat contributions, phase quality, net shaft work, thermal
efficiency, fluid entropy changes and independent energy residuals.
It uses illustrative efficiencies, no pressure-drop model, and no assigned
thermal reservoirs; it does not claim measured machinery performance or total
reservoir entropy production. The earlier hot-inlet steam-train example remains
a partial train and retains its stated limits.

Integrated regressions exercise complete cycles and cold-feed region crossings,
recover the existing 4 MPa Gibbs Rankine reference without changing its h/s datum
or device behavior, and check regional outlet-domain and cancellation refusal.
Independent evaluations of the original Gibbs equations verify high-temperature
throttling and ideal/inefficient expansion from 30 MPa/1500 K to 1 MPa, including
the original h, u+p*v and s balances.

## Adiabatic water and steam stream junctions

`WaterFlowChart::state_at_quality` and `state_at_specific_enthalpy` expose the
same original Gibbs, dense or regional state owners to finite stream devices.
Quality ingress retains the authored vapor fraction directly, including changes
too small to resolve at an absolute enthalpy datum. Charts without saturation
refuse quality requests. Enthalpy ingress preserves the original bounded
inverse, iteration counts and dense source-evaluation counts; neither ingress
adds heat, shifts source references or introduces a pressure-work term.

`phase::water::flow::across::mixing` adds `WaterFlowChart::mix_streams` for
prescribed positive water/steam mass flows on the same exact prepared pressure
chart. One original p,h inverse closes the conserved mass and enthalpy fluxes.
The resulting liquid, two-phase, vapor or supercritical state retains the
original IF97 caloric reference, with latent heat included once. Regional gaps,
ambiguous inverses and unsupported source states remain refusals. A single
inlet retains its exact original state without an inverse.

The steady adiabatic junction has no storage, shaft work, modeled pressure loss,
kinetic/potential energy or momentum balance. Its independent internal-energy
and pressure-flow-work fluxes close the same first law; flow work is already
inside enthalpy. Raw-source entropy flux supplies the signed entropy production,
without a thermal entropy-seam correction or clipping. Existing per-throughput
energy/entropy tolerances and endpoint pressure-closure gates apply before any
result is returned. The operation allocates no memory and polls cancellation
through its inlet reduction, original inverse and final acceptance.

Focused source regressions check cold-water/steam condensation against the
original saturation lever rule, nonlinear liquid temperature against a direct
Gibbs enthalpy inverse, flow scaling and port subdivision, dense supercritical
mixing without a vapor inventory, original-chart binding and transactional
early/late cancellation.

The `water_regional_rankine` example also composes a regenerative 20 MPa /
823.15 K cycle with a 1 MPa open feedwater heater. Original steam and condensate
enthalpies determine the bleed fraction needed for a saturated-liquid heater
outlet. The actual mixer, two feed pumps and two turbine stages retain their
original source states. Condenser-side flow and work carry the unextracted
mass fraction, and the complete cycle checks independent h and u+p*v balances
and closed fluid entropy at an explicit 1 kg/s boiler-throughput scale.

## Coherent thermal-fluid transport view

`fluid_transport::ThermalFluidTransport` extends the existing transport view
with its original absolute temperature and optional isobaric volumetric
expansion. Density, Cp, conductivity and viscosity remain the same immutable
checked tuple used by mechanical consumers. Original ideal-gas owners derive
beta = 1/T at fixed composition after their EOS/caloric/acoustic coherence gate.
Source-resolved liquid snapshots use their explicitly admitted temperature axis;
without a prepared density path, they require one explicitly typed
absolute-temperature coordinate. Four-property resolution alone retains unknown
expansion; a missing coefficient never means zero expansion.

`resolve_liquid_state_with_expansion` and
`LiquidTemperatureLaw::try_new_with_expansion` opt into the original card's
instantaneous isobaric `volumetric_thermal_expansion_coefficient` [1/K]. The
same immutable liquid bundle owns all five requirements and source receipts,
including the explicit dynamic/kinematic viscosity basis. Full-query and
temperature refresh pin and reevaluate the whole bundle; loss of expansion
support refuses the update without publishing a partial transport state.
Known negative and zero expansion remain distinct from missing data. Existing
thermal-fluid consumers receive this coefficient through the original liquid
type; there is no separate convection material or substituted gas law.
Independently declared density and expansion tables are not promoted into a
thermodynamically consistent EOS, a phase-stability claim or a finite-temperature
Boussinesq validity certificate. Focused `liquid_expansion` regressions cover
source retention, both viscosity bases, signed coefficients and refusal.

`LiquidTemperatureLaw::with_density_expansion` (or the corresponding
`ResolvedLiquidState` method) prepares an explicitly bounded temperature line
through the **same selected density claim**. The existing matdb line owner
admits scalar plateaus with explicit temperature support, native linear curves,
and fixed-pressure/composition slices through sampled surfaces. Every source
breakpoint and receipt is retained in shared immutable storage. At each later
temperature, `beta = -(d rho/dT)/rho` uses the local interpolation derivative
and the actual refreshed density, alongside the same Cp, conductivity and
viscosity. This follows directly from `V = mass/rho` at fixed mass. It adds no
gas-law assumption, fitted coefficient, smoothing or finite-difference step.
Tangents use original curve knots or the original triangle's analytic gradient,
so clipping a path to nearly equal evaluated densities cannot manufacture zero
expansion. Surface-cell selection reuses the canonical pinned point query;
intervals without a representable interior query explicitly refuse.
The prepared state feeds the existing thermal-fluid, natural-convection and
immersed-body consumers without separately authored expansion numbers.

The density-derived and independently declared coefficients are explicit source
choices. Binding a density path to a state that already retains an independent
expansion claim refuses, so no selected claim is silently overridden. The
liquid's `volumetric_expansion_per_k` accessor and its thermal-transport snapshot
return the same coefficient for either admitted source choice.

Source intervals must be continuously supported. Discrete samples and an
isolated density scalar cannot supply a derivative. Zero and negative
expansion are preserved; interior source knots with unequal one-sided slopes
refuse. Declared interval endpoints use the derivative from inside that
interval. Updates retain the original source and all nonthermal coordinates;
changing pressure or composition requires a newly admitted path rather than
reusing its old tangent. A density interval does not extend any other transport
claim or establish a finite-temperature convection approximation. Source-value
uncertainty is retained and is not promoted to a derivative error bound.

Focused G0/G1/G3 checks cover native curve and surface derivatives, independent
specific-volume differentiation, zero/negative expansion, endpoint and kink
semantics, missing/discrete support, and unchanged states on invalid refresh.
`fs-airflow`'s `material_cell_convection::thermal_source` regressions change only
the density source slope and verify the resulting Rayleigh number, heat-transfer
coefficient, Robin boundary and accepted backward-Euler thermal history against
independent formula/energy calculations.

An explicit authored-state constructor can declare its own temperature and
optional beta, without acquiring source-data authority.

`ResolvedThermochemicalRealGasState` supplies the same thermal view from its
original Peng-Robinson branch. Its local expansion is
`beta = -(dp/dT)_v / (v * (dp/dv)_T)` at fixed composition, derived from the
same analytic pressure derivatives that supply Cp and sound speed. Thermal
and momentum views retain identical PR density/Cp and original transport
coefficients; ideal-reference gas density or `1/T` is never substituted.
G1 tests compare against independent constant-pressure density perturbations;
G3 checks the dilute ideal limit and source-preserving dense-state substitution.
The derivative does not establish phase equilibrium, a finite-temperature
Boussinesq interval or empirical validity of the chosen critical/transport data.
`ThermochemicalPengRobinsonModel::new` accepts pair interactions in the supplied
species-binding order and permutes both matrix axes with the bindings. Reordering
three or more species therefore preserves their original pair laws and every
derived density, caloric, acoustic and thermal coefficient. G0/G1/G3 regressions
check malformed matrices, direct EOS agreement and all ternary permutations.

This is one evaluated state, not a temperature-dependent fluid law or evolving
reservoir. Consumers must keep source/domain and frozen-property assumptions,
and must verify agreement when binding mechanical and thermal snapshots.
Three focused unit definitions cover coherent gas data, authored expansion and
typed liquid temperature. Runtime qualification is reported separately.

## Source-bound Ohmic constitutive graph

`conductor::graph::MaterialOhmicNode` exposes the original scalar conductor as a public constitutive-law node with temperature-change and signed-current inputs, and signed-voltage and nonnegative-Joule-power outputs. It preserves the selected electrical property key and basis, original material bundle and claims, and continuously supported temperature slice, including fixed-context sampled surfaces. Conductivity is queried before conversion to resistance. Its piecewise-consistent tangent differentiates the native resistance and Joule-power equations; interior knots use the right interval and the final endpoint uses the left interval. Construction admits positive, representable native transport across every electrical source knot and refuses unsupported spans, unowned response coordinates, and phase-fraction evolution. The default uses fixed geometry and memoryless scalar DC transport; sourced free expansion can be bound as described below. Thermal storage, time, contacts, skin effects and phase evolution remain with their existing physical owners.

### Source-defined expansion in electrical and thermal feedback

`PreparedFreeExpansion::at_temperature_with_tangent` derives normalized
principal stretches, their exact piecewise temperature derivatives and the
volume derivative from original source segments. It uses the right interval
at interior knots and the left interval at the final endpoint, with one binary
search and no runtime allocation or source query for this geometric evaluation.
Preparation differentiates the native curve segments or the original surface
cells; it never differences clipped or merged endpoint values. Consequently a
nonzero expansion tangent survives equal rounded endpoint lengths, and an
inserted reference or another principal direction's knot creates no artificial
corner. Surface derivatives use the canonical inward tangent-cone owner even
when adjacent temperatures have no representable interior. An unrepresentable
tangent refuses its derivative request without invalidating geometry-only
evaluation. The narrow-span, signed-principal and source-surface regressions in
`tests/expansion_surface.rs` exercise these cases through the existing geometry
and exposed-area consumers.

`MaterialOhmicNode::with_free_expansion` requires the exact original material
card, reference query, full temperature interval and an explicit current axis.
It derives L=L0*s_parallel and A=A0*s_transverse1*s_transverse2 at the same
temperature as the original electrical property. The graph tangent composes
electrical resistivity/conductivity and all three geometric derivatives.
If an expansion property already belongs to the original conductor bundle,
its selected claim must remain identical. Physical cancellation of material
and geometric temperature derivatives is admitted.

The existing graph/enthalpy solver thus closes strain data, geometry,
resistance/current, Joule heat and fixed-mass temperature on one endpoint.
Five material regressions cover normalized tangents, all current axes, both
electrical source bases, derivative cancellation, fixed geometry and source
pin/refusal behavior. One coupled regression compares three strain tables over
two continuing steps against an independent finite-capacity quadratic and
terminal-work balance. Both focused test targets pass `cargo check`; execution
remains pending. This is unconstrained, quasistatic, stress-free expansion.
Source density claims, reference thermal geometry and fixed-mass caloric
storage retain their original owners. No mechanical work, moving thermal
boundary areas or constrained deformation is inferred from the new binding.

## Radial distributed heat storage on the original timeline

`thermal::transport::stack::LayerStackGeometry` opts
`PreparedLayerStackConduction::try_new_geometric` and
`PreparedLayeredThermalBody::try_new_geometric` into fixed cylindrical or
spherical shells. Existing constructors retain planar geometry and identity.
Physical thicknesses, original conductivity/capacity sources and areal contact
resistances stay distinct from resistance coordinates. Steady source-integrated
flux uses A_ref integral(dr/A(r)) and R'' A_ref/A_interface on the original
Kirchhoff/contact kernel. The inner face defines A_ref; actual ordered face areas
are exposed separately, and radial callers must not use `area_m2()` as inventory.

Each homogeneous layer is divided into equal physical-volume cells with
volume-centroid temperature radii. Actual shell volumes, fixed supplied masses,
radial interface areas and surface coverage feed the original thermal cells,
source deposition, integrated two-leg fluxes, sparse energy-coordinate Newton
and adaptive backward-Euler owner. A contact transfers the same signed joules
on both sides. Face-area sums round outward by at most one ulp only where needed
to avoid false overcoverage during original sequential subtraction.

Boundary exchange retains the original terminal-cell temperature semantics;
no exact surface temperature or extra half-cell reservoir resistance is
inferred. Radial symmetry, positive inner radius, fixed geometry and sensible
storage are required; cylinder endcaps are insulated. Spatial and time
refinement remain independent requirements. This path does not establish latent
heat, solid-centre symmetry, moving boundaries or fluid transport. Existing
planar laminate and plane-wave heating consumers explicitly refuse radial bodies.

`tests/radial_layered_storage.rs` covers analytic/contact power, nonlinear
integrals, physical inventory, native surfaces, finite-difference tangents,
implicit conservation, source/boundary energy, adaptive evolution, manufactured
spatial refinement, planar compatibility and source/geometry refusal.

## Adaptive physical-face boundary conduction

`BoundLayeredThermalBody::advance_adaptive_to_with_surface_boundaries` carries an
explicit fixed set of `LayeredSurfacePatch` values through the existing adaptive
backward-Euler owner. Every coarse and fine interval prepares the original
material/film face laws from that interval's actual authored bath temperatures,
conductances and coverage. Native boundary schedules remain the event owner.

The existing admission seam rechecks complete selected-channel coverage and the
actual material-face source support after events, including events exactly at
the target. Reservoir removal, changed coverage incompatible with the fixed
patch map, unsupported material faces and cancellation return no partial
trajectory. Accepted half steps alone supply heat and entropy; rejected/coarse
work does not enter those accounts. Existing unselected channels and empty-map
behavior retain their original semantics.

The face map preserves fixed geometry and the explicit independent normal
streamtube approximation. It adds no lateral spreading or second time
integrator. Local temporal estimates remain distinct from spatial error and
experimental validation. Focused tests compare physical-face cooling with an
independent exponential, scheduled bath/film changes with closed-form accepted
half steps, and target-event refusal with original-state preservation.

Native Rust execution of these new regressions remains unverified under the
existing dependency-build memory and disk limits; pinned formatting, source/API
review and independent scalar fixture calculations passed.

## Physical terminal-face conduction on the layered thermal timeline

`LayeredThermalStep::with_surface_boundaries` maps selected original ambient
or reservoir channels to the actual inner/outer faces of a sensible layered
body. Planar, cylindrical and spherical face-to-volume-centroid conduction
uses the existing `ReservoirLayerGeometry` conductivity-integral/film solver.
The response retains the material-face temperature, signed outward heat, secant
conductance and exact derivatives in [cell temperature, bath temperature] order.
Conductivity is integrated across its original chart, not frozen at a trace.

At most four explicit face/channel patches are prepared once. Each selected
original channel must be completely covered, and patches cannot exceed a
physical face. Partial areas represent independent normal streamtubes or angular
sectors, with no inferred lateral spreading. The original film, bath, area,
source and caloric inventories remain authored inputs; zero film conductance is
exact insulation. Unselected channels retain their original cell-temperature
semantics, and an empty map restores the original step path.

The existing boundary-coefficient ingress replaces each selected channel's
conductance. The native residual and energy/entropy owners account for its heat
once; the tridiagonal Jacobian uses the actual outward dQ/dTcell. Source support
is required throughout the material profile. A bath outside the conductivity
chart is legal only when the solved material face stays supported; perfect
contact still requires the prescribed face itself to be supported. The existing
strict all-endpoint reservoir evaluator retains its original policy.

The extension shares the original bounded scalar solver, source-integral inverse
and cancellation checks. It supplies no new convection correlation, boundary
layer model or lateral spreading law. Three unit regressions and four integration
regressions cover geometry, partial areas, nonlinear flux/tangents, source support,
native conservation and transactional refusal. Pinned formatting, source/API
review and independent scalar oracles passed. Native compilation and execution
of these new regressions remain unverified under the existing shared-memory and
storage limits.

## Distributed phase-changing layers on the original enthalpy timeline

`thermal::transport::stack::phase_transient` binds explicitly declared reference
masses and equal-volume planar, cylindrical or spherical control volumes to the
original `PreparedPhaseConductivity` and equilibrium enthalpy/phase chart. Births
specify h directly, including distinct latent states at the same temperature.
Physical interface and end-face resistance coordinates preserve actual mass,
volume and area. Optional face reservoirs retain the bulk half-cell and declared
surface resistance.

The original conduction observer and implicit/adaptive enthalpy chain own
evolution and energy/entropy accounting, including zero latent dT/dh and source
dk/dh. Geometry and reference mass remain fixed; cellwise conductivity requires
spatial refinement. A single cell can exchange with both physical end-face
baths: its inner face uses the original ambient channel and its outer face the
original reservoir channel. Each retains its own temperature, finite material
and film resistance, signed heat and bath entropy. Native exposed-area rounding
is admitted only within a bounded coverage tolerance; scaling both areal
resistances preserves the original physical face conductance and its derivative.
No boundary heat is reclassified as a generated source or counted twice.
Source-domain, original-law, inventory, topology and cancellation failures
refuse without publishing a changed timeline. This bounded MR19 transport
increment does not establish empirical source qualification or phase kinetics.

The explicit solid-cylinder and solid-sphere constructors admit a regular
zero-area centre followed by optional concentric material shells. Equal-volume
cells retain their physical volume centroids, supplied reference masses and
original phase conductivity. The centre has zero radial heat flow and refuses
an inner bath. Outer half-cell resistance is normalized at its actual physical
face; no zero-radius reference area or artificial inner surface enters the law.
The existing explicit and implicit/adaptive enthalpy owners retain all evolution,
generation, contact and energy/entropy accounting.

Ten integration regressions cover latent inventory, nonlinear multi-material
contact balances in all three geometries, physical face heating, a computed
Stefan front with spatial/time refinement, transactional refusal, source-varying
two-bath latent transport, and equal-opposite bath heat with positive conduction
entropy under explicit and adaptive evolution. Solid-centre regressions add
physical inventory/symmetry admission, nonlinear latent heating through an
outer finite film, and a computed generated-heat melt front with radial
temperature/front refinement and a single sensible-plus-latent heat account. Formatting,
parsing and independent scalar oracle checks passed. Native integration execution,
including the Stefan convergence thresholds, remains unverified because the
shared workspace's dependency build exceeded its memory/storage capacity.

## Spatial heat storage inside sourced layers (2026-10-29)

`thermal::transport::stack::transient::PreparedLayeredThermalBody` prepares
1–64 planar material layers with up to 512 spatial control volumes. Each layer
has an original conductivity query, a separately selected capacity query,
explicit whole-layer mass and a positive uniform subdivision count. Mass-based
capacity must name that whole mass; native volumetric capacity must name the
area-times-thickness volume. Preparation divides the amount among cells and
retains one per-cell storage curve per layer, shared by that layer's immutable
caloric law. Spatial refinement does not duplicate source preparations.

Capacity and conductivity retain the same immutable material card and agree
on shared typed fixed conditions. The caller still declares compatible source
populations and caloric/work constraints; common card identity does not prove
matched specimens. All capacity and conductivity paths share a positive-width
temperature interval containing the explicit common energy datum. Original
receipts, amounts, physical ordering, geometry and subdivision counts contribute
to prepared identity. No density law is inferred from the required mass input.

The bound owner's `cell_temperatures` query validates the same original cells,
borrowed material laws and contacts as conduction, then returns temperatures in
physical thickness order without recomputing flux. `prepared_body` and original
per-layer cell masses support mechanical bindings that preserve every thermal
subdivision and compare each reference mass, rather than only the global sum.

`BoundLayeredThermalBody` attaches the exact original borrowed storage laws and
checked cell/material/link IDs. Its birth events create independent cell
enthalpies, then contacts. Each internal face uses the original two-layer
conduction kernel across the adjacent half-cell distances. Contact resistance
applies only at actual material interfaces; subdivisions inside a layer have
perfect contact. The response retains each edge's heat flow, interface
temperatures, endpoint derivatives and continuity diagnostics.

Call `step_coefficients` on every candidate of the existing backward-Euler
residual owner. It checks the actual law objects, cell mass and geometry, source
layer identities and link orientation/area before returning sparse conductance
replacements. The original timeline owns energy, sources, evolution, optional
sensible entropy and acceptance. This adapter supplies no second time integrator
or nonlinear solver. Initial conductances are only initial values. Successful
edge kernels allocate nothing; candidate temperature, response and coefficient
containers use O(cells) storage. Checkpoints and explicit edge iteration budgets
bound work; cancellation/refusal changes no physical state.

The model holds geometry, mass, composition and source conditions fixed. Only
the two axial faces are represented; lateral surfaces are insulated. Boundary
heat sources/reservoirs use the existing owner's events. Each cell has one
temperature; spatial refinement and time-step refinement remain independent.
Phase change, moisture transport, changing cure, expansion work and experimental
qualification are not supplied by this sensible finite-volume assembly.

All four tests in `xtask/tests/material_layered_transient.rs` pass against the
production modules. An independently derived insulated cosine mode exercises
2/4/8/16 cells and 1/2/4/8 backward-Euler steps with separate spatial and temporal
convergence. Manufactured nonlinear heterogeneous layers check half-cell and
contact physics, exact caloric integrals, finite-difference tangents, energy
and entropy. A six-cell FSRI gypsum/mineral-wool/plywood wall retains real source
receipts, and a material substitution changes its transient from the same initial
temperature field. Binding/domain/budget/cancellation refusals preserve the
original timeline.

## Pressure-controlled ice/liquid freezing and thawing

`phase::water::freezing::FreezingWater` combines original IAPWS Ice Ih,
the released pressure-dependent Ih melting curve, and original IF97 Region 1
liquid properties. This is an ambient-pressure equilibrium owner: the melting
temperature must be in 273.15..=273.16 K, where the two source charts overlap.
Pressure cannot exceed p_melt(273.15 K), and original liquid source admission
also applies near the triple point and boiling boundary. The immutable
temperature window brackets melting, with ice no colder than 50 K and a
stable Region 1 liquid upper endpoint. A latent-only window is supported.
High-pressure freezing, supercooled liquid, other ice polymorphs, vapor,
nucleation and resolved interfaces are outside this owner.

The prepared pair retains original h, u, s, specific volumes and the actual
melting inverse diagnostics. Both Gibbs mismatch and the pressure-local fusion
entropy defect L/Tm-(s_liquid-s_ice) have explicit independent admission budgets.
Neither property reference is shifted. The defect is a source-consistency
diagnostic, not modeled entropy production.

Temperature states and explicit liquid mass fractions share one source owner.
Enthalpy inversion uses physical lever-rule storage on the fusion plateau and
safeguarded Newton/bisection with analytic Cp on each sensible branch. Exact
source knots and complete phase depletion are retained; temperature alone at
melting is ambiguous. Liquid/ice mass fractions and volume fractions remain
distinct because the phases have different densities.

Fixed-mass heat steps evolve H and independently close Delta U+p*Delta V.
Stable original-source increments and direct endpoint subtraction must both
fit the unchanged caller joule budget. Reports separate sensible and fusion
enthalpy while retaining actual phase masses, physical expansion/contraction,
boundary work, raw source entropy and the pressure-local integral of dh/T.
Latent enthalpy is counted once. The external pressure load permits expansion;
this is not a sealed rigid-vessel pressure model.

Six native tests pass: released normal-pressure melting and original-source
identity/Clapeyron consistency; complete freezing/thawing and expansion;
independent Cp/CpT quadrature; subdivided, latent-only, small and zero heat;
source overlap, matching and work refusal; and chart identity, numerical
representation and early/late cancellation. The tests compile the exact
original production sources against the new ice thermochemistry library.

## Freezing on the shared conservative thermal timeline

`phase::water::freezing::thermal::FreezingWaterThermalLaw` connects the same
ice/liquid enthalpy owner to `fs_thermal::enthalpy::ThermalLaw`. Preparation
obtains outward-rounded whole-interval Cp bounds independently for original
Ice Ih and IF97 Region 1. The maximum reciprocal lower Cp bounds dT/dh over
the entire immutable caloric path; the physical fusion plateau has exactly
zero slope. A latent-only chart prepares no sensible-branch certificate.
The retained bounds are reused across thermal steps and checkpoints.

The adapter uses the existing temperature inversion and stable pressure-local
integral of dh/T. It does not introduce a second heat solver, apparent-Cp
regularization, a sampled stability minimum, or a new energy reference.
Unresolved source bounds and cancellation yield no prepared law. Runtime
domain errors propagate through the existing transactional thermal owner.

Its entropy increment preserves the supplied stored Delta h, using the
original ice/liquid/fusion Delta S*/Delta h as mean reciprocal temperature.
Coincident inverse temperatures use Delta h/T, retaining small sensible heat
and exact latent storage despite source inverse or phase-fraction rounding.
The source-state entropy and heat-step accounts keep their original increments.
Additional `ThermalLaw` regressions cover single-ULP ice and liquid heating,
full and tiny fusion increments, reversal, domain rejection and cancellation.

`cargo run -p fs-material --example water_freeze_thaw` evolves a uniform
50 g sample at 101325 Pa through cold and warm reservoir contact. The original
sources give Tm=273.152519079769 K and fusion enthalpy 333426.514501674 J/kg.
An explicitly supplied 5 W/K contact conducts to 253.15 K until 300 s and then
293.15 K until 650 s. Temperature and phase inventory emerge from conserved
enthalpy: the sample crosses the fusion plateau, becomes cold ice, thaws, and
returns to warm liquid. The example separately reports liquid mass and volume
fractions, physical free volume, Delta U, pressure work, heat, sample/reservoir
entropy and the temporal entropy defect.

The source Gibbs mismatch is 0.016931408 J/kg and the explicit pressure-local
fusion entropy correction is 6.198518076e-5 J/(kg K); neither is hidden as
physical entropy generation. The native run's maximum independent energy
residual is 6.84e-10 J. This is arithmetic/solver closure, not experimental
prediction accuracy. The example assumes internal equilibrium and supplied
contact conductance. Reference thermal geometry remains fixed while material
volume is reported separately; no resolved front, nucleation, freezing-induced
stress, moving contact area or rigid-vessel behavior is inferred.

Ten integrated native tests pass: the six freezing-core tests plus four
thermal tests for complete branch Cp bounds and latent-only preparation,
reservoir-driven freeze/thaw with timestep-refined entropy defect, isolated
ice/water contact with conserved heat and positive thermal entropy, and
cancellation/domain atomicity. The complete original fs-thermal library,
the freezing material namespace and the executable example build and run
without numerical stubs. Full workspace/DSR verification remains unavailable
because the required sibling constellation checkouts and DSR runner are absent.

## Sparse implicit tangents for layered thermal bodies (2026-10-29)

`thermal::transport::stack::transient::implicit::LayeredThermalStep` evaluates
the existing `BackwardEulerStep` residual and attaches a derivative of that
same material-backed candidate. Its unknowns are cell energy increments in
the original owner's CellId order. For `R = delta_E - dt Q(T(E))`, the retained
operator is `J = I - dt (dQ/dT) diag(1/C)`, using original per-cell capacity and
physical endpoint derivatives of integrated conductivity/contact flux.
Secant conductance is not substituted for a Newton derivative.

The operator admits exactly its bound body cells and links. Authored constant
cell powers, convection on the original exposed areas, and reservoir exchange
remain with the original heat-rate owner. Its boundary slopes include reservoir
coverage and fixed convection/conductance. Pending event admission and target
left-limit semantics remain those of `BackwardEulerStep`; larger coupled graphs
continue to use the existing coefficient adapter and their complete Jacobian.

`with_cell_energy_j` adds signed interval-only joules in that same original
cell order. They become constant powers for one implicit solve, leaving authored
sources and accepted history unchanged. A coupled owner must solve any physical
state dependence of those inputs externally; the sparse thermal Jacobian holds
them fixed and never claims to differentiate that enclosing feedback.

Only three diagonals are retained. Tangent assembly and matrix actions require
O(cells) storage/work after material flux evaluation. Forward, transpose and
Jacobi inverse-diagonal actions allocate nothing and validate arithmetic and
cancellation before publishing their bounded output. Unequal capacities and
nonlinear conductivity generally make the operator nonsymmetric; the transpose
uses the actual opposite diagonals. Every tangent is bound to its exact trial,
with original temperatures, capacities and link responses available for clients.

The existing thermal trial still owns residuals, optional entropy and independent
local/total joule acceptance. `LayeredThermalStep::{solve,solve_with_entropy}`
adds bounded damped Newton using reusable `fs-la::tridiagonal` factors and
O(cells) Newton storage. Every proposal reevaluates original integrated k(T),
contact fluxes and storage; backtracking enforces the original source domain
and residual decrease. Work limits count rejected trials and the final optional
entropy evaluation. Only original local and aggregate joule acceptance publishes
an endpoint. Cancellation, source refusal and nonconvergence return no state.
The original residual/material evaluation retains O(cells) trial allocations.

The solve has no explicit diffusion CFL restriction, but it remains first-order
backward Euler and requires independent time and spatial refinement. It does
not silently subdivide an interval; events at the target stay pending. There
is no trajectory adjoint, condition-number certificate or unconditional
nonlinear-convergence claim. This single-interval solve does not supply
automatic time-error control.

`BoundLayeredThermalBody::advance_adaptive_to` connects the exact
integrated-conductivity/contact sparse solver to the existing `fs-thermal`
adaptive step-doubling owner. Only the accepted two-half-step trajectory
contributes heat and entropy. Local estimates use changes in original specific
storage and temperature; independent local and whole-advance joule gates still
decide physical acceptance. One residual-evaluation ceiling charges accepted,
rejected and nonconverged endpoint trials, including entropy queries. All cells
and contacts must already be active. Scheduled Boundary events apply at their
exact times, including the requested target; topology changes, unsupported phase
observables, missing sensible entropy, source-domain failures and cancellation
refuse without publishing a partial prefix. Suggested next step size is explicit.
These are local temporal estimates, not certified global error bounds, spatial
error estimates or trajectory adjoints. The single-interval
`LayeredThermalStep::solve` retains its existing pending-target-event semantics.

Focused tests in `xtask/tests/material_layered_implicit.rs` cover independent
constant coefficients, sources, convection and reservoir coverage;
nonlinear finite differences verify the full energy-coordinate matrix and the
transpose inner-product identity. Invalid/cancelled matrix actions preserve
outputs. The six-cell FSRI wall fixture now uses the production sparse solve
for its 3,000-second interval and original energy/entropy acceptance. Analytical
two-cell tests check unequal-material relaxation, conservative exchange,
workspace reuse and deterministic replay. An independent exponential decay
reference checks first-order time convergence over 8, 16 and 32 intervals.
Work-budget, source-domain, scope, event, identity, overflow and cancellation
refusal fixtures require the initial timeline to remain unchanged. These are
numerical/source-path checks, not experimental material validation.


## Source-backed spatial water freezing

`phase::water::freezing::transport` binds the unchanged original freezing chart
to independent liquid IAPWS SR6-08(2011) and empirical freshwater-ice Yen (1981)
conductivity. This composition deliberately admits standard atmospheric
pressure (101325 Pa) only, with the complete caloric temperature window inside
the transport source domains. Absent phases are not extrapolated. Melting
transport values are prepared once; sensible-phase queries retain original
source admission. Partially frozen cells require the caller to select series
or parallel layers. These directional closures use physical phase VOLUME
fractions from distinct source densities, never mass fractions or an inferred
isotropic mixture. The ice correlation supplies no pressure response or
experimental uncertainty envelope. Transport does not alter h, u, s or fusion
enthalpy.

`freezing::spatial::FreezingSlab` assembles a stationary one-dimensional plane
reference mesh on the existing `fs-thermal` timeline and geometric conduction
observer. Initial source density fixes conserved reference mass. Both faces
may have explicit temperature reservoirs with finite areal resistance;
otherwise they are insulated, as are all lateral faces. Scheduled reservoir
temperature changes use the original event owner. State-dependent internal and
half-cell boundary conductances are refreshed on every actual CFL-limited
substep. The latent plateau remains the exact original h-to-phase relation;
no apparent-Cp regularization or separate latent source is introduced.

Accepted advances independently check original-source Delta H = Q and
Delta U + p*Delta V = Q against the supplied joule tolerance. The original
thermal and reservoir entropy accounts accompany the separate temporal
entropy defect. Cell profiles expose actual enthalpy, temperature, phase mass
and volume fractions; inventory reports distinct phase masses and free volumes.
The reported equivalent frozen REFERENCE length is thickness times frozen mass
fraction, not a mechanically resolved interface or physical ice thickness.
Reference heat-transfer geometry remains frozen while free expansion and
pressure work are reported separately. This conduction model does not solve
convection, stress, rigid-vessel pressure, nucleation, supercooling, front
kinetics or a deforming boundary. It has no unconditional spatial/time accuracy
claim. Source-backed mesh/time refinement and the shared owner's independent
Stefan reference are numerical implementation evidence, not experimental
validation of every freezing specimen.

The runnable `water_freezing_front` example computes an advancing phase field
inside a 10 mm water layer and reports the resolved profile, phase inventories,
heat, internal energy, pressure work and entropy. Focused tests also cover
adiabatic equilibrium, a bath reversal, source/chart binding, energy balance,
work limits and cancellation without altering accepted state.

## Native thermoelectric constitutive graph

`thermoelectric::graph::MaterialSeebeckNode` binds the original absolute thermopower source to the public constitutive graph through three temperature-change inputs and signed current. Every value is supplied by `ResolvedSeebeckSpan`, retaining both terminal Peltier transfers, the original bulk Thomson exchange, and electrical supply. The selected claim, query quantities, fixed external context, and continuously supported temperature slice remain pinned. Its consistent tangent integrates piecewise source derivatives, including reversed edges and the isothermal limit, without differencing nearly equal thermopowers. The node declares reversible transfer with no local storage and explicitly zero dissipation; graph execution and equation-compiler execution retain that declaration, and the compiler codec preserves it. Resistance, Joule heat, Fourier conduction, caloric storage, physical time, and acceptance of complete coupled balances remain with their existing owners. No phase evolution, magnetic transport, relative-to-reference thermopower substitution, or measured authority beyond the original source is inferred.

## Radial material conduction and finite thermal feedback

`PreparedLayerConduction::implicit_link_laws` now supplies borrowed original
Kirchhoff/contact fluxes to the native implicit chain and graph. It reuses
`TwoLayerConduction`, including its analytic endpoint derivatives, with the
same caloric dT/dh and physical link/coverage/distance checks. Source-local
temperature spans, pure-phase restrictions and the admitted left derivative
at an upper source endpoint remain explicit. No database query, identity
construction or allocation occurs in a link evaluation. These are source-free
quasi-steady half-cell reconstructions; transient spatial/time refinement and
material uncertainty remain separate. Boundary faces retain their native
bulk-conductivity model unless an `ImplicitLayerBoundary` is attached. Focused
layer-network regressions exercise the actual chain/graph solve, affine-k
closed-form heat and derivatives, disjoint source
spans across contact, source endpoints and transactional refusals.

`thermal::conduction::LayerGeometry` admits planar, concentric cylindrical and
spherical source-free conduction. It lowers the exact metric
`A_contact * integral(ds/A(s))` to equivalent planar distances, retaining
physical radii and contact area. The original `TwoLayerConduction` integrates
each source k(T), solves contact temperatures and supplies endpoint derivatives.
The contact jump uses the actual area; thin-radius differences use stable
logarithm/difference forms. No source value is inferred from a material name.

`PreparedLayerConduction::new_geometric` binds those metrics to the original
active contact, source-qualified pure-phase caloric objects and per-side
coverage. `LayerConductionObserver` supplies their integrated conductances to
the existing explicit enthalpy/heat-rate/CFL owner; `resolve` supplies the same
implicit endpoint residual. One signed internal transfer changes both finite
stores, with no external heat or duplicate storage. Source/domain/event and
cancellation refusal discard the full proposed advance. Source selection and
geometric metrics are prepared once; flux evaluation performs no database query.

Explicit accepted endpoints and events are admitted by the same physical
two-layer interface solve. Each source covers only its own bulk-to-face path;
a finite contact jump may bridge disjoint temperature intervals. The endpoint
check retains the starting-state heat samples and introduces no extra heat,
receipt allocation or source query. A supported bulk temperature does not hide
an unsupported contact face on the final step. Focused G1/G3 regressions cover
disjoint hot/cold sources, reversed interface orientation, exact single-step
heat and entropy; a G0 regression refuses a final cold-face source exit while
both bulk temperatures remain supported.

The six `radial_conduction` regressions cover independent cylindrical/spherical
resistance, nonlinear conductivity integrals/tangents, thin shells and the
exterior limit, data-only rho/Cp/k changes between finite solid/liquid cells,
native time refinement/energy/entropy, implicit endpoint parity and refusal
with replay. Their integration target and the existing `layer_conduction` and
`layer_conduction_network` targets compile in the scoped Cargo check; runtime
execution is not established by that check.

This is stationary radial conduction at fixed mass/reference geometry with
quiescent liquid and explicit symmetry assumptions. The legs have no independent
thermal storage; transient spatial resolution and coefficient-freezing error
require refinement. No convection, fluid pressure/motion, confinement correction,
phase-transition propagation or experimental accuracy is claimed.

### Shared electrical geometry for native thermal networks

`MaterialOhmicNode::resolve_from_state` lowers an already resolved complete
electrical/thermal bundle through the same prepared L/A and optional free
expansion used by graph evaluation. It retains the card, full source
requirements/keys/claims, quantity schema and fixed conditions; only the
declared temperature may change. No second property query occurs. Original
source-domain admission belongs to the canonical resolved bundle, and geometry
remains bounded by the prepared interval. A focused regression checks parity
with ordinary node resolution, bundle preservation and foreign/out-of-range
refusal. The `conductor_graph` integration target compiles; this check does not
establish runtime execution of that regression.


