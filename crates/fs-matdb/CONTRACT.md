# CONTRACT: fs-matdb

> Status: ACTIVE (bead 5hmy complete: all five PRs landed). Owns the
> immutable typed material-data schema, its fail-closed insertion
> boundary, the material/constitutive card layer with supersedes
> lineage, the ordered interface-system card, and the query path where
> every answer is `Evidence<PropertySample>` + `PropertyUsageReceipt`
> with replay-verified receipt completeness. Recorded residuals for
> follow-up beads: query-time joint-uncertainty correlation refs,
> tensor/distribution payloads with frame-transform receipts (arrive
> together), wider explicit fusion policies, and the curated seed
> dataset beyond its first gas-association tranches (bead 1sxe).

The normalized-pack boundary now carries typed joint covariance/correlation
blocks and unit/basis normalization receipts. Query-time propagation of those
joint statistics remains a recorded no-claim until usage receipts select and
bind them explicitly.

## Purpose and layer

"Real material properties" are NAMED CONDITIONS with uncertainty, not
labels. fs-matdb is the L1 data layer every physics claim inherits its
weakest load-bearing material datum through: typed (dims-checked at
insertion), immutable (append-only; conflicting claims coexist),
provenance-complete (source + license + content hashes are load-bearing
fields). Layer: L1. Runtime deps: fs-qty, fs-evidence, fs-blake3 ONLY.
This crate owns no executable closures and no per-run state (L3
adapters do); it never imports L2 transforms, L3 state types, or L6
persistence.

## Public types and semantics

- `PropertyKey` — (name, `fs_qty::Dims`) pair. The dims are registered
  on first insertion; reusing a name with different dims is a refusal,
  so one vocabulary name can never quietly alias two physical
  quantities.
- `Provenance` — source citation + license + optional artifact
  `ContentHash`. Empty source or empty license REFUSES: a bare value is
  not a datum, and unlicensed data cannot enter the store.
- `ObservationDataset` — specimen/process record, method/instrument,
  observation artifact hash, covariance/censoring caveats, provenance.
  Content-addressed (`org.frankensim.fs-matdb.observation-dataset.v1`).
  Observations are what claims point at: a citation-only claim (no
  observations) is admitted but can never be Validated-class downstream —
  specimen/process match requires observations.
- `PropertyValue` — `Scalar { value, dims }` or `Curve { abscissa,
  abscissa_dims, knots, dims }` (strictly increasing finite abscissae,
  ≥ 2 knots). Tensor/distribution/model-parameter payloads arrive with
  PR-2 so frames and state schemas land together.
- `UncertaintyModel` — `Unstated` (admitted AND marked; PR-4 never lets
  it launder into a certified band), absolute `HalfWidth`, or
  `RelativeHalfWidth`, each with confidence strictly in (0, 1).
- `InterpolationPolicy` — `LinearInside`, `ConstantWithinValidity`, or
  `TabulatedOnly`. Extrapolation is never implicit.
- `PropertyClaim` — key + value + `fs_evidence::ValidityDomain` (THE
  single validity type; a second competing type is forbidden) +
  uncertainty + interpolation + observation refs + provenance.
  Content-addressed (`org.frankensim.fs-matdb.property-claim.v1`) over
  every semantic field with exact float bits.
- `NormalizedPack` — the bounded, versioned L1 artifact emitted by an
  offline compiler after raw-source and licensing policy. It carries a
  canonicalized `ClaimSet`, exact raw-envelope hash, compiler identity,
  retained redistribution decision, typed `JointStatistics`, and sorted
  `NormalizationReceipt`s. `to_bytes`/`from_bytes` own the `FSMATPK\0` v1
  binary wire format; decoding reconstructs every ordinary fs-matdb object,
  reproduces its semantic id, and then byte-reproduces the full stream.
  `from_bytes_verified` additionally requires the externally pinned whole-pack
  identity before any top-level metadata/statistics mutation can be accepted.
- `NormalizedInterfacePack` — the `FSINTPK\\0` v1 wrapper for one ordered
  `InterfaceSystemCard`. It reuses a complete `NormalizedPack` as the sole
  claim/observation/statistics/normalization payload, then binds both material
  states, both texture frames, medium, optional third body, environment, and
  named history. Decode re-admits the nested claim pack and verifies both its
  hash and the reconstructed interface-card hash before accepting canonical
  bytes. V1 carries no constitutive model cards; model-law transport requires a
  separately versioned binding and cannot be smuggled through this wrapper.
- `NormalizedModelPack` — a separate bounded `FSMODPK` v1 transport for
  immutable `ConstitutiveModelCard`s. Model cards are not laundered into
  scalar property claims: the pack retains each law/version, dimensioned
  parameter block, state convention, validity domain, source hashes, and
  provenance, plus one `ModelNormalizationReceipt` for every parameter and
  both endpoints of every validity interval. Cards encode in full-content-hash
  order; source hashes and receipt targets are strictly ordered and
  deduplicated. `from_bytes_verified` pins the whole model-pack identity before
  decoding, while each serialized card identity must independently reproduce.
  A downstream law adapter (for example the NASA-9 adapter in
  `fs-thermochem`) must still validate the exact law-specific schema.
- `NormalizedSpeciesPack` — a separate bounded `FSSPCPK` v1 transport for one
  source-declared thermochemical association. Its pack id is exactly the
  `fs_qty::chemistry::SpeciesId`; v1 retains a positive coherent-SI molar mass,
  the exact `gas`/`ideal-gas` standard-state convention, positive reference
  pressure, an opaque elemental-reference convention id, source hashes, and
  provenance. Exactly one `SpeciesNormalizationReceipt` for molar mass and one
  for reference pressure retain the source literal hash and linear SI transform.
  Canonical decode re-runs every runtime structural admission gate and
  byte-reproduces the stream; `from_bytes_verified` first pins the whole-pack
  identity.
- `JointStatistics` — a named, explicitly ordered component block for one
  observed dataset. `StatisticMember` addresses a scalar claim or one curve
  knot's abscissa/ordinate, so multiple knots never collapse into one random
  variable and one observation may carry multiple disjoint named blocks.
  Covariance and optional correlation use packed lower-triangle order.
  Matrices must be finite and positive semidefinite; covariance has a
  nonnegative diagonal, correlation is bounded to `[-1,1]` with an exact unit
  diagonal and must bit-reproduce covariance-derived correlation, every
  member's claim must cite the owning observation, and blocks for one
  observation must be member-disjoint.
- `NormalizationReceipt` — a structured `NormalizationTarget` (claim
  component, uncertainty, validity endpoint, or joint-covariance entry), hash
  of the exact source literal, six-base dimensions, and affine
  `si = source * scale + offset` basis transform. Claim/covariance targets are
  resolved and dimension-checked against retained data. Optional source/target
  frame names are retained as a pair; scalar/curve packs make no
  tensor-rotation claim. Uncertainty and covariance transforms must have exact
  positive-zero offsets because widths, fractions, and covariances cannot be
  translated. Uncertainty magnitudes and variance diagonals additionally
  require positive scales. Lower/upper receipts for one validity axis must
  agree on six-base dimensions even though the shared validity type does not
  yet carry those dimensions itself.
- `ClaimSet` — the PR-1 append-only container. `register_observation`
  and `insert_claim` are the ONLY mutations; both are fail-closed and
  idempotent by content identity. `claims_for(name)` returns EVERY
  claim for a name in insertion order — conflicting observations stay
  separate `PropertyClaim`s, and fusion is an explicit query-time
  policy (PR-4), never a map overwrite that invents a canonical value.
- `MatDbError` — total, typed refusals: `DimsMismatch`,
  `MissingLicense`, `MissingSource`, `NonFinite` (with exact bits),
  `UnusableValidity`, `InvalidUncertainty`, `MalformedCurve`,
  `UnknownObservation`, `EmptyParameterBlock`, `NonFiniteParameter`,
  `RevisionNotZero`, `SupersedesMismatch`.
- `LawId` / `ConstitutiveModelCard` (PR-2) — a law's stable (id,
  version) identity, its canonical dimensioned parameter block
  (nonempty, finite; BTreeMap = one canonical hash order),
  `StateSchemaVersion`, `InitialStatePolicy` (`ZeroInternalState` or
  `RequiresDeclaredState` — the card never implies a state it does not
  declare), the shared `ValidityDomain`, calibration source hashes, and
  load-bearing provenance. Content-addressed
  (`org.frankensim.fs-matdb.constitutive-model-card.v1`). DATA ONLY:
  the executable law-node protocol is L3 fs-material (bead kagp).
  `canonical_parameters_hash` separately addresses only the sorted parameter
  names, exact value bits, and six-base dimensions under
  `org.frankensim.fs-matdb.canonical-parameter-block.v1`; it is never a model
  identity by itself and must travel with law/version/state-schema and
  implementation-contract identities. Its preimage is the identity version as
  little-endian `u32`, the parameter count as little-endian `u64`, then each
  BTreeMap-ordered name, exact little-endian `f64` bits, and six signed
  dimension bytes, with every part framed by a little-endian `u64` byte count.
  Minting first runs the ordinary card
  admission gates, so empty or non-finite parameter blocks cannot acquire an
  authoritative canonical hash.
- `MaterialStateId` / `MaterialCard` (PR-2) — a NAMED MATERIAL STATE
  (chemistry + phase + temper/process + revision) carrying its claim
  set, its model cards, the by-key and by-law indexes, and explicit
  lineage. Constructors are `assemble` (revision 0 only —
  `RevisionNotZero` otherwise) and `supersede` (revision exactly +1,
  `supersedes` bound to the predecessor's content hash; the
  predecessor is untouched and stays retrievable). No mutable access
  exists after construction. The card hash
  (`org.frankensim.fs-matdb.material-card.v1`) binds the id, schema
  version, lineage link, every claim/observation content id, and every
  model-card hash — so it binds the full transitive content.
- `SurfaceSpec` / `SystemContext` / `InterfaceSystemCard` (PR-3) — an
  ORDERED interface system: surface A (material state + opaque
  texture-frame id; blank refuses), surface B, and the system context
  (medium, optional third body, environment, NAMED history state —
  each blank member refuses). Friction, wetting, contact conductance,
  wear, and adhesion are claimed against the SYSTEM, never against an
  unordered bulk pair: `(a, b)` and `(b, a)` hash differently, history
  is identity-bearing, and wetting is a solid–liquid–gas system (the
  liquid is the medium, the gas is the environment). Content-addressed
  (`org.frankensim.fs-matdb.interface-system-card.v1`) over both
  ordered surfaces, the full context, and the transitive
  claim/observation/model identities.

## Invariants

- APPEND-ONLY: no API removes or mutates a stored claim or observation.
  Same content re-inserts idempotently under the same id; different
  content for the same key coexists under a distinct id.
- DIMS AT THE DOOR: a claim whose payload dims disagree with its key —
  or whose key name is already registered with other dims — never
  enters.
- PROVENANCE IS LOAD-BEARING: no source ⇒ refusal; no license ⇒
  refusal. A citation alone does not make a value Validated: that
  requires specimen/process-matched observations (enforced at the PR-4
  color boundary; the linkage is stored here).
- VALIDITY IS THE SHARED TYPE: `fs_evidence::ValidityDomain`, with its
  intersection composition law, is the only validity representation.
- CONTENT IDENTITY: claim/observation ids are domain-separated BLAKE3
  over length-framed semantic fields with exact IEEE-754 bits.
- PACK CANONICALITY: observations and claims encode in semantic-id order;
  joint blocks encode in `(observation id, block id)` order; normalization
  receipts encode in structured-target order. Duplicate unordered entries
  refuse. Portable packs require finite values and canonical positive zero.
- MODEL-PACK CANONICALITY: cards encode in full semantic-id order and require
  nonempty sorted source hashes plus an exact provenance artifact. Every
  normalized parameter and validity endpoint has exactly one structurally
  linked receipt; dangling, missing, duplicate, parameter-dimension-mismatched,
  or endpoint-transform-incoherent receipts refuse. The model pack does not
  infer a law schema from parameter names; validity-axis dimensions remain
  compiler-declared because the shared validity type does not store them.
- SPECIES-PACK CANONICALITY: one pack carries exactly one species association;
  its id, positive molar mass, exact v1 phase/EOS, reference pressure, elemental
  convention, sorted unique source hashes, and provenance are identity-bearing.
  Both numeric fields have exactly one dimension-matched, positive-scale,
  positive-zero-offset receipt targeting the explicit six-base SI basis.
- JOINT STATISTICS STAY JOINT: covariance/correlation are never collapsed into
  nominal values or caveat text. Member order and every lower-triangle entry
  are identity-bearing normalized bytes.

## Error model

Total functions; no panics in library paths. Ordinary claim/card/query
refusals use `MatDbError`; normalized artifact refusals use `PackError` with
stable field/resource/byte-offset or semantic-identity context. Non-finite
ordinary-data refusals carry exact bits.

## Determinism class

Fully deterministic: pure data structures, `BTreeMap`/content-id ordering,
content hashes over canonical byte encodings, and exact IEEE-754 transport.
Pack covariance admission derives dimensionless correlations from exact stored
variances for pairwise/source-correlation checks, then gates each covariance
and correlation matrix with one fixed-order, outward-rounded interval LDLT
pass. Zero variance requires an exact zero row/column; a pivot is accepted only
with a positive lower bound or an exact structural-zero proof. Negative,
overflowed, or rounding-ambiguous pivots refuse rather than tolerating a false
PSD claim. The gate is deliberately incomplete for ill-conditioned or
rank-deficient valid matrices: refusal is not evidence of indefiniteness. This
is an admission check, not a solver or a re-estimation of source statistics.
Query interpolation retains its separately declared evaluator semantics.

## Cancellation behavior

All operations are synchronous and contain no I/O or solve. Pack decode is
explicitly bounded by byte, collection, per-block-member, and cumulative cubic
PSD-work budgets. Reference admission uses ordered lookup rather than nested
linear scans. Declared outer/member/observation counts are checked against
minimum remaining payload bytes before proportional semantic allocation. The
PR-4 query path polls at claim granularity if selection ever becomes
super-linear.

## Unsafe boundary

`#![forbid(unsafe_code)]` inherited from workspace lints; no unsafe.

## Feature flags

None.

## Conformance tests

`tests/conformance.rs`: scalar and curve round-trips with dims checks
at insertion; the same-name/different-dims registry refusal;
license-missing and source-missing refusals; NaN scalar / NaN knot /
unordered knots / short-curve refusals; NaN validity-axis refusal;
invalid-confidence and negative-half-width refusals; dangling
observation-reference refusal; conflicting-claims preservation (two
different densities for one key coexist under distinct ids, both
retrievable, insertion order preserved); content-hash stability
(re-insertion is idempotent) and sensitivity (any semantic field change
moves the id); observation registration idempotency.

`tests/cards.rs` (PR-2): genesis assembly with by-key/by-law indexes;
nonzero-genesis-revision, empty-parameter-block, NaN-parameter, and
unlicensed-model refusals; supersession chain 0→1→2 with predecessor
hashes bound and predecessors immutable; model-card hash field
sensitivity (parameter value/version/state-policy/validity/sources);
material-card hash binding claims, models, and the named-state id.

`tests/interface.rs` (PR-3): surface-order and history hash
sensitivity; three-phase wetting with advancing/receding hysteresis as
coexisting claims; unnamed-texture-frame and blank
medium/history/third-body refusals.

`tests/query.rs` (PR-4 + PR-5): honest evidence slices (Estimate band
from stated uncertainty; Unstated → numerical no-claim PROVEN never to
certify); complete point-sensitive receipts; extrapolation, unknown
property, and non-finite point refusals; explicit-fusion ambiguity
refusal and the observation-backed preference policy; curve
interpolation inside knots with exact-hit tagging and beyond-data
refusal (fail-closed ordering: validity containment gates before
evaluation, so `MissingQueryAxis` is only reachable through
unconstrained validity); and the RECEIPT-COMPLETENESS MUTATION BATTERY —
`ClaimSet::verify_receipt` replays the query from the receipt's own
fields, 11 per-field mutations all refuse typed
(`ReceiptMismatch { field }` / `UnknownPolicyTag` /
`EvaluatorVersionDrift` / the replay's own refusals) and every mutation
moves the receipt content hash.

`tests/pack.rs`: fixed v1 canonical byte-length/hash golden, exact
byte/semantic round-trip, and externally pinned whole-pack verification;
permutation-invariant block/receipt ordering; named
curve-knot members and multiple disjoint blocks per observation; typed
covariance/correlation preservation and aggregate PSD-work admission;
malformed shape, negative diagonal, exact-negative/rounded-zero PSD regression,
invalid/inconsistent correlation, overlapping joint blocks, unknown/uncited
reference, negative zero, unlinked/dimension-mismatched/translated statistical
normalization target, negative statistical scale, contradictory validity-axis
dimensions, partial frame receipt, untrusted-count preflight,
truncation, trailing-byte, and semantic-id tamper refusals.

`tests/model_pack.rs`: G0 canonical two-card permutation and exact binary
round-trip; G3 complete receipt coverage, parameter-dimension linkage,
validity-endpoint transform coherence, portable positive-zero/source/provenance
gates, nested card-identity tampering, whole-pack tampering, and trailing-byte
refusals.

`tests/species_pack.rs`: G0 source/receipt permutation canonicalization and
verified exact binary round-trip; G3 phase/EOS/positive-value/provenance gates,
complete dimension-linked receipt coverage, pack/species identity binding,
whole-pack tampering, untrusted-length preflight, and trailing-byte refusals.

`xtask/tests/matdb_pack_cli.rs`: G3 compilation of the committed methane,
nitrogen, oxygen, argon, carbon-dioxide, water-vapor, and carbon-monoxide seed
manifests twice into byte-identical, identity-verified species packs, including
retained NASA source/license/standard-state associations and separately bounded
NIST displayed-precision agreement checks. The material-pack battery likewise
compiles each committed bulk tranche twice and verifies exact identity,
condition flags, source/license linkage, no-claim exclusions, and bounded
comparison evidence, including NASA-CR-115153's inhibited water/ethylene-glycol
coolant specification and the conflict-preserving N0602-001 nitrile O-ring
compatibility tranche. The NGYC N42 battery additionally pins supplier, sintered
grade, and coating identity while proving that incompatible SI/CGS energy-product
prints survive as distinct claims. The NACA TN 2680 battery compiles the
supplier-minimum-purity 2,2,4-trimethylpentane tranche and verifies all fifteen
apparatus-bound Table I maximum-flame-speed rows without promoting the report's
empirical fit to a bulk-material law. The FACE G CDTRF-G 2023 v1 battery checks
five exact volume fractions, their unit sum, and two conflict-preserving
calculated-RON prints for the same named formulation.
The NASA UAM winding-insulation battery compiles the source-linked MW-16C
polyimide-wire thermal-class pair, the selected 0.08 mm Nomex 410 slot liner,
and the CoolTherm EP-2000 actual/omitted cure steps as three separate packs. It
also proves that assembly-level PDIV, modeled stress, dielectric, conductivity,
and lifetime-law claims are not laundered into the bulk constituents.
The NASA-CR-195445 rotary-coating battery compiles the source-linked PS-200
feedstock composition and the OMC Test 3/Test 6 before/after RMS finishes while
pinning the aluminum-alloy configuration, incomplete process identity, crack,
local breakthrough, and foreshortened-run caveats. It refuses to infer a wear
rate, friction law, coating thickness, transferable durability, or generic
Wankel-housing authority.

## No-claim boundaries

- Storing a claim asserts NOTHING about its truth: fs-matdb records who
  said what, where it holds, and how uncertain it is — evidence colors
  and certification live in fs-evidence and are assigned at the PR-4
  query boundary, never at insertion.
- `Unstated` uncertainty is a marked absence, not zero uncertainty.
- Interface/system properties (friction, wetting, contact conductance)
  are NOT expressible as bulk claims — they wait for
  InterfaceSystemCard (PR-3), because they are system+history
  properties, never unordered pair constants.
- Total magnetic moment is integrated over an actual body and is
  deliberately NOT a storable geometry-free material scalar; store
  magnetization/specific-moment curves instead.
- No complete seed dataset ships in this crate. The repository's first
  `data/matdb/seed-v1` tranches contain raw, offline-compiled gas-species
  associations, six exact-temperature Aluminum 6061-T6 scalar claims, and
  four exact-temperature OFHC Copper scalar claims. The first PTFE/Teflon
  tranche adds four polynomial-derived exact-temperature thermal claims while
  retaining the source's missing grade and process identity. A NASA LaRC PEEK
  plate tranche adds four exact-point conductivity claims, four exact-point
  specific-heat claims, and one source-temperature-unspecified density. These
  sit beside fourteen exact-temperature AISI 4140 mechanical-property claims
  and three small-sample AISI 1045 tensile claims. The first pinned AISI 52100
  bearing-steel tranche adds six actual-composition claims, five Rockwell C
  scale readings, and four exact retained-austenite fractions. The first
  pinned AISI 9310 gear-steel tranche adds nine nominal-composition claims,
  two conflicting case Rockwell C statements, one core Rockwell C statement,
  and one carburized case-depth claim. Two named NASA/NAPC gear-oil packs add
  nine explicitly unit-bearing flash-point, pour-point, specific-gravity, and
  total-acid-number claims while refusing the source's unitless viscosity
  numbers. The named Rheolube 2000 bearing-grease pack adds one NLGI
  named-scale reading, one exact-temperature density, and one time/temperature-
  bound oil-separation fraction. The named Pennzane SHF X-2000 bearing-oil pack
  adds three exact-temperature kinematic-viscosity values, one viscosity-index
  scale reading, flash and pour points, and one exact-temperature density. The
  first pinned gray-iron tranche adds nine reported composition/carbon-equivalent
  values, four quantified microstructure fields, and two low-precision
  graph-digitized room-state properties. The NASA-CR-115153 coolant tranche
  adds six formulation-bound endpoints, one source-condition-unspecified
  density, one source-condition-unspecified conductivity, and three exact-point
  transcriptions of an approximate heat-capacity law. The N0602-001 nitrile
  tranche adds one TGA semi-volatile claim, two aromatic-content-keyed
  absorbed-fuel claims, three partitioning statistics, two conflicting printed
  swell slopes, one regression intercept, and one coefficient of determination.
  The first supplier-pinned N42 tranche adds one remanence, one coercivity, and
  two conflict-preserving maximum-energy-product claims for NGYC sintered,
  nickel-coated cubes. A separate Jinshan N42 tranche adds exact `25 degC` and
  `120 degC` remanence, intrinsic-coercivity, and maximum-energy-product
  endpoints plus the source's two interval-average temperature coefficients
  for its pristine, wire-cut commercial sintered state. The NACA TN 2680
  iso-octane tranche adds one
  supplier-minimum-purity claim and fifteen atmospheric Bunsen-flame maximum-
  speed observations with explicit temperature, oxidizer, flow, burner, and
  missing-condition metadata. The FACE G CDTRF-G 2023 v1 tranche adds five
  composition claims on the source's volume-fraction basis and two separately
  provenance-linked calculated-RON claims (`94` and `93.9`) for the identical
  formulation.
  The gas associations do not define air or exhaust mixture compositions,
  humidity, or combustion completeness. The Aluminum claims do not define a
  continuous constitutive curve or a general-purpose design card; their NIST
  polynomial-fit errors lack the confidence metadata needed by the current
  statistical uncertainty variants and therefore remain observation caveats
  with explicitly `Unstated` runtime uncertainty. OFHC thermal conductivity
  pins RRR=100, while its specific-heat observation preserves the source's
  unstated RRR rather than laundering that condition across properties. The
  combined NIST source scope also does not select between UNS C10100 and
  C10200. The PTFE/Teflon claims bind only NIST's source label and exact
  temperature points; no resin grade, crystallinity, filler, processing state,
  density, continuous constitutive law, or application-specific authority is
  inferred. The PEEK claims bind one NASA THERMIC plate and the source's
  thermal-model inputs, retain the stricter 300-525 K range, and expose missing
  grade/process identity plus density test temperature; no continuous curve,
  seal, tribology, or lifetime authority is inferred. The 4140 claims bind
  NASA's QQ-S-624 heat 137M186, one-inch bar, normalize/harden/oil-quench/
  temper schedule, and Rockwell C33 condition; they are not generic values for
  the grade or its separately reported C44 branch. The 1045 claims bind one
  cold-drawn bar and ASTM E8 specimen series;
  their Student-t intervals are derived from three printed replicates under an
  explicit iid-normal assumption. Because that source omits test temperature,
  each claim requires the fail-closed
  `source_test_temperature_known = 0` query axis and supplies no temperature
  validity interval. The flag acknowledges missing metadata; it does not make
  the values temperature-independent. No hardness or joint covariance is
  admitted from that source. The 52100 claims bind one NASA
  consumable-vacuum-melted ingot, its reported chemistry, a common
  austenitize/oil-quench/first-temper spine, and five separately keyed
  second-temper states. Rockwell C is a named empirical scale reading in
  dimensionless storage, not a ratio quantity. The `<2%` retained-austenite
  result stays censored in an observation, the predictive equation's
  `+/-1`-point accuracy is not reused as table-measurement uncertainty, and no
  hot-hardness curve or long-time stability is claimed. The 9310 claims bind
  one NASA CVM heat and the complete carburize/quench/subzero/temper/grind/
  stress-relief schedule. The source's detailed C58 and summary-level C60 case
  statements remain separate provenance-linked claims with `Unstated`
  uncertainty; no average or preferred value is invented. Lubricant-dependent
  surface-fatigue lives remain excluded as gear-system properties. The
  NASA/NAPC oil packs retain public code, batch, basestock-class, and
  specification context, but the source says formulation chemistry is
  proprietary and omits the viscosity unit. No viscosity claim, formulation
  identity, or gear-system fatigue/EHD result is laundered into those bulk-oil
  cards. The Rheolube 2000 pack binds its Pennzane base-oil association,
  thickener identity, and three typical source values, but refuses approximate,
  unitless, extrapolated, or test-system rows as exact bulk claims. It supplies
  no generic bearing-life or service-envelope authority. The Pennzane oil pack
  similarly binds one named MAC base fluid and seven typical properties while
  refusing the unrepresentable Celsius-interval coefficient, system wear scar,
  and ambiguous extracted vapor-pressure exponents. It supplies no formulated-
  oil, continuous-curve, or bearing-life authority. The gray-
  iron claims bind Wang et al.'s S2-S charge, chemistry, Sr-FeSi inoculation,
  EN-1561 Type II mould, and fully pearlitic/type-A-graphite state. Its Figure 8
  strength and conductivity centers retain only graph-supported precision;
  one-standard-deviation bars without confidence metadata remain observation
  caveats, not runtime half-widths. Exact test temperature is absent and must
  be acknowledged with `source_test_temperature_known = 0`. The NASA-CR-115153
  coolant pack binds one inhibited water/ethylene-glycol specification rather
  than a generic mixture: printed formulation endpoints remain separate bounds,
  and no midpoint or ethylene-glycol balance is inferred. Density and thermal
  conductivity retain missing-temperature and missing-pressure flags; the
  report's unstated BTU convention remains explicit for conductivity and three
  exact-point SI evaluations of its approximate heat-capacity relation. Those
  points do not authorize a continuous law, interpolation, a freeze/boil
  envelope, corrosion/compatibility behavior, or service-life authority. The
  N0602-001 claims bind one source-coded nitrile O-ring and the reported JP-8/FT
  matrix, not the nitrile family. Compound formulation, supplier, cure,
  hardness, lot, dimensions, exposure temperature/duration, raw points, and
  coefficient uncertainty remain absent. The summary-row and plot slopes
  (`0.451` and `0.463`) coexist without fusion; the approximate prediction-
  interval overlap remains observation-only, and the regression intercept is
  not a certified shrinkage value. No generic seal, ozone/aging, permeability,
  compression-set, compatibility, or service-life authority follows. The NGYC
  N42 claims bind the supplier-named sintered grade and nickel-coated
  cube family from one CC-BY-4.0 article. The source does not say whether the
  values were measured or supplier-nominal and omits test temperature, method,
  lot, chemistry, detailed process, intrinsic coercivity, recoil permeability,
  demagnetization curves, and temperature coefficients. Its printed
  `318.3 kJ/m^3` and `42 MGOe` maximum-energy-product representations normalize
  to different SI values and therefore remain separate claims; neither is a
  resolved design allowable. The Jinshan N42 claims remain a separate material
  state and may not be fused with NGYC: they bind a different supplier, the
  source-described `10 mm x 10 mm x 6 mm` pristine wire-cut specimen, and a
  NIM 6500C measurement campaign. The source's rounded endpoints and printed
  coefficients are retained independently because they do not reproduce one
  another exactly. The plotted `80 degC` curve is not digitized, and missing
  production lot, chemistry, original sinter schedule, coating, magnetization
  protocol, temperature-control details, field sweep, repeats, uncertainty,
  recoil permeability, and irreversible-loss boundary remain fail-closed.
  Neither the two endpoint states nor the interval-average coefficients are a
  continuous law or extrapolation authority. The NACA TN 2680 claims bind a
  supplier-claimed `99.6 mol%` minimum-purity 2,2,4-trimethylpentane fuel and
  fifteen Table I
  maximum-flame-speed rows to the reported atmospheric Bunsen-flame apparatus.
  Supplier, lot, exact assay, impurity composition, row-level pressure, exact
  maximizing equivalence ratio, raw images, dispersion, and confidence metadata
  remain absent. The observations are not a pure-fluid property card, gasoline
  surrogate specification, reaction mechanism, transferable burner law, or
  authority for density, viscosity, surface tension, heat capacity, latent heat,
  vapor pressure, octane rating, or the report's empirical fit. The FACE G
  CDTRF-G 2023 v1 claims bind one named surrogate's five published volume
  fractions and two internally inconsistent calculated-RON prints. They are not
  an assay of FACE G, a fungible gasoline recipe, pure-component cards, or a
  combustion mechanism. Component suppliers, lots, purities, preparation and
  mixing state, volume-contraction treatment, CFR-engine RON, and statistical
  uncertainty remain absent; no density, viscosity, vapor pressure,
  distillation, heat capacity, latent heat, flame speed, ignition delay, storage,
  compatibility, or emissions authority follows. Bead 1sxe still owns the
  N42 recoil-limit evidence and the remaining curated material/property
  and interface-system dataset. The WO 2018/125520 Formulation 8 claims bind
  one patent-table comparator with four named source-era commercial products,
  exact as-added mass fractions, and eight printed 5W-30 property results. They
  are not authority for the SAE grade generally, current trademarked products,
  component interchangeability, detailed additive chemistry, a continuous
  viscosity law, aging, tribology, wear, oxidation, deposits, emissions,
  engine life, or service intervals. Component lots, manufacturing dates,
  final blending protocol, method editions, repeats, dispersion, and confidence
  metadata remain absent. Permission to redistribute patent text is explicitly
  not a patent-practice or trademark license.
  The NASA UAM winding-insulation claims bind three constituents from one
  NASA electric-aircraft motorette campaign. The MW-16C `240 degC` and
  `20000 h` values are cross-bound as one NEMA/ASTM thermal-class basis, not an
  Arrhenius or service-life law. The Nomex 410 claim is only the `0.08 mm`
  selected slot-liner thickness. The CoolTherm EP-2000 process state completed
  the `180 degC` post-cure but explicitly omitted the manufacturer-recommended
  `210 degC` final step; it is not a fully cured generic epoxy card. Wire
  vendor/lot, raw classification evidence, Nomex conditioning/metrology, epoxy
  lot, full cure schedule, degree of cure, uncertainty, and confidence remain
  absent. Campaign PDIV, hot-spot, and stress results depend on the assembled
  geometry and history and therefore remain excluded from the bulk packs, as
  do dielectric strength, thermal conductivity, fatigue, adhesion, moisture,
  activation-energy, and service-life claims.
  The NASA-CR-195445 rotary-coating claims bind the patent-linked PS-200 powder
  composition and four RMS finish observations to the report's exact air-cooled
  OMC Test 3 and Test 6 configurations. Feedstock mass fractions are not
  as-sprayed phase fractions. The aluminum alloy, SX-331 chemistry, coating
  thickness/profile, spray parameters, finish metrology, repeats, and
  uncertainty remain absent. The cracked Test 3 TBC and Test 6's uneven honing,
  local PS-200 breakthrough, leakage, hot spots, and `1.5 h` early stop are
  retained no-claim boundaries; no wear rate, friction coefficient, thermal
  property, transferable durability, generic Wankel identity, or service-life
  authority follows. Source-text redistribution grants no patent-practice or
  trademark rights.
  No equilibrium computation happens here (fs-thermochem consumes phase data;
  this crate only stores it).
- The L1 pack codec does not parse handbooks, CSV, NASA tables, license text,
  or other raw formats and does not decide whether terms permit
  redistribution; those are L6/offline compiler responsibilities. A nonblank
  retained decision is provenance, not legal advice.
- The model-pack codec proves transport integrity and generic card admission,
  not that an arbitrary card is NASA-9, Arrhenius kinetics, or any other named
  physical law. Executable consumers retain that law-specific validation
  obligation; no kinetics executor or model-card/species linkage is claimed.
- A species pack preserves an explicit source association; it does not
  authenticate chemical identity, validate the molar mass or elemental
  convention against an independent authority, link a NASA/kinetics card, or
  supply thermodynamic, kinetic, equilibrium, or transport evaluation.
- Joint statistics are preserved and validated but are not yet selected or
  propagated by `PropertyUsageReceipt`; no query result may claim correlated
  uncertainty until that later authority surface binds the exact block.
- Frame names in normalization receipts record provenance only. Scalar/curve
  payloads do not carry tensor components, rotation matrices, or a claim that
  a frame conversion was physically valid.
- `ValidityDomain` does not yet retain axis dimensions. A validity-bound
  normalization target proves that the claim/axis/endpoint exists, but its
  six-base dimensions remain compiler-supplied provenance until the shared
  validity schema grows a typed axis registry.
- A pack may contain already-normalized SI values with no transform receipt;
  the L1 codec therefore proves every present receipt is linked, not that the
  receipt set exhausts every numeric field. Source-format policy owns that
  completeness check in the offline compiler.
