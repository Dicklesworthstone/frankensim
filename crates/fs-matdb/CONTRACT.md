# CONTRACT: fs-matdb

> Status: ACTIVE (bead 5hmy, PR-1 of 5). Owns the immutable typed
> material-data schema and its fail-closed insertion boundary. The query
> path with `Evidence<PropertySample>` + `PropertyUsageReceipt` is PR-4;
> MaterialCard/ConstitutiveModelCard are PR-2; InterfaceSystemCard is
> PR-3.

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
- `ClaimSet` — the PR-1 append-only container. `register_observation`
  and `insert_claim` are the ONLY mutations; both are fail-closed and
  idempotent by content identity. `claims_for(name)` returns EVERY
  claim for a name in insertion order — conflicting observations stay
  separate `PropertyClaim`s, and fusion is an explicit query-time
  policy (PR-4), never a map overwrite that invents a canonical value.
- `MatDbError` — total, typed refusals: `DimsMismatch`,
  `MissingLicense`, `MissingSource`, `NonFinite` (with exact bits),
  `UnusableValidity`, `InvalidUncertainty`, `MalformedCurve`,
  `UnknownObservation`.

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

## Error model

Total functions; no panics in library paths; every refusal is a typed
`MatDbError` naming the gate and the offending field (non-finite
refusals carry exact bits).

## Determinism class

Fully deterministic: pure data structures, `BTreeMap` iteration order,
content hashes over canonical byte encodings. No floating-point
arithmetic is performed on stored values in PR-1 (validation compares
and hashes only), so there is no rounding class to declare yet;
interpolation arithmetic arrives with the PR-4 evaluator and will
declare its class then.

## Cancellation behavior

Not applicable in PR-1: all operations are bounded, small, and
synchronous (no solves, no I/O). The PR-4 query path polls at claim
granularity if selection ever becomes super-linear.

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
- No seed data ships in this crate (bead 1sxe owns the curated
  dataset); no equilibrium computation happens here (fs-thermochem
  consumes phase data; this crate only stores it).
