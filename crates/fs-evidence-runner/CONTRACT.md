# CONTRACT: fs-evidence-runner

## Purpose and layer

Layer: **TOOL**.

`fs-evidence-runner` owns the frozen Runner V2 orchestration vocabulary and its
bounded, pure validation rules. The base schema makes command intent, terminal
state, diagnostics, limits, budgets, typed values, logical publication
selection, capability policy, and presented identity references explicit
before lifecycle, process, storage, verification, or admission implementations
are allowed to consume them.

This initial slice is declarations and validation only. It does not execute a
case, emit lifecycle records, parse hostile bytes, access a filesystem, spawn a
process, persist an object, publish a bundle, verify scientific evidence, or
mint authority.

As a TOOL crate, it may orchestrate lower layers. No product-layer crate may
depend on it, and no type in this crate becomes a lower-layer authority merely
because it is accepted by a validator.

## Product generation and wire version

The public product-generation name is `RunnerSpecV2`. The API generation is
exactly `2`.

The first frozen wire schema version is exactly `1`. It has an explicit
no-predecessor policy: there is no wire version zero, legacy decoder,
compatibility alias, implicit migration, or best-effort reinterpretation.
API generation and wire version are distinct typed fields and cannot be
substituted for one another.

Every canonical identity schema introduced by this Runner V2 family uses wire
version `1` and a `.v1` domain even when its nominal Rust type ends in `V2`.
Domains follow exactly:

```text
org.frankensim.fs-evidence-runner.<schema-kebab-name>.v1
```

The schema name table is closed and contains no aliases. A later wire-shape
change must rotate its version and domain; an API-generation change alone does
not silently move wire bytes.

## Public types and semantics

This crate owns the following base-schema families and their pure validators.
Later leaves may consume or seal registrations for them, but may not redeclare
their Rust or wire types.

### Terminal states and refusal reasons

`ProofExitV2` is a closed `u16` catalog:

| Code | Variant |
| ---: | --- |
| 0 | `Pass` |
| 10 | `Failed` |
| 11 | `Refused` |
| 12 | `NoData` |
| 13 | `Stale` |
| 14 | `EnvironmentInvalid` |
| 15 | `Blocked` |
| 16 | `Unsupported` |
| 17 | `NotRun` |
| 18 | `Cancelled` |
| 19 | `TimedOut` |
| 64 | `Usage` |
| 70 | `InternalError` |

Only `Refused` carries a `RefusedReasonV2`. Its closed `u16` catalog is:

1. `InvalidEvidence`
2. `NonCanonicalEvidence`
3. `EvidenceIdentityMismatch`
4. `EvidenceTampered`
5. `LimitExceeded`
6. `UnsafeArtifactPlacement`
7. `ArtifactCollision`
8. `LifecycleViolation`
9. `PolicyRefused`
10. `AuthorityBoundaryViolation`
11. `MigrationRefused`

Unknown codes refuse. Family-specific reasons use separately registered,
bounded namespaces and cannot extend either catalog.

### Commands, profiles, dispositions, roots, and paths

`RunnerCommandV2` is exact:

| Code | Variant |
| ---: | --- |
| 0 | `List` |
| 1 | `Check` |
| 2 | `SelfTest` |
| 3 | `Run` |
| 4 | `Negative` |
| 5 | `Replay` |

`RunProfileV2` is `1 Smoke` or `2 Full`.

`ArtifactDispositionV2` is `1 LifecycleOnlyNoBundle` or
`2 DurableBundleRequired`. `List` carries typed absence rather than a third
disposition.

`RootClassV2` is `1 InputArtifactRoot`, `2 OutputArtifactRoot`, or
`3 Other` carrying a registered policy ID.

`PlatformPathProfileV2` is:

1. `PosixDescriptorRelativeV1`
2. `WindowsHandleRelativeV1`
3. `ContentStoreObjectKeyV1`

An unavailable platform profile produces `Unsupported`; it is not silently
substituted.

Command applicability is table-driven:

- `List` has no family, mode, profile, budgets, publication selection, or
  disposition.
- `Check` uses its sealed preflight manifest, accepts no caller family, mode,
  or profile, and is `LifecycleOnlyNoBundle`.
- `SelfTest` uses its sealed bounded internal manifest, accepts no caller
  family, mode, or profile, and is `LifecycleOnlyNoBundle`.
- `Run` requires one explicit registered family mode and one explicit `Smoke`
  or `Full` profile. There is no default mode.
- `Negative` names one immutable negative case whose sealed manifest supplies
  family, mode, and profile.
- `Replay` names one immutable source manifest whose lineage supplies family,
  mode, and profile.
- `Run`, `Negative`, and `Replay` require `DurableBundleRequired`.

Omitted, ambiguous, duplicated, or command-inapplicable selectors are `Usage`.

### Lifecycle record vocabulary and state-bearing roles

`LifecycleRecordKindV2` is exact:

1. `RunStart`
2. `CaseStart`
3. `FamilyRow`
4. `CaseTerminal`
5. `RunSummary`
6. `RunTerminal`

`StateBearingRecordRoleV2` is exact:

1. `PreRunDiagnostic`
2. `ExecutedCaseTerminal`
3. `SuppressedCaseTerminal`
4. `RunTerminal`

The role/state matrix is closed:

- `PreRunDiagnostic` permits `Usage`, `Refused`, `NoData`, `Stale`,
  `EnvironmentInvalid`, `Blocked`, `Unsupported`, `Cancelled`, `TimedOut`, or
  controlled `InternalError`. It requires exactly one matching base diagnostic
  and has no lifecycle record.
- `ExecutedCaseTerminal` permits every non-`Usage` state except `NotRun`.
- `SuppressedCaseTerminal` permits only `NotRun`.
- `RunTerminal` permits every non-`Usage` state except `NotRun`.
- `Pass` forbids a diagnostic.
- `Failed` requires `case.conformance_mismatch`.
- `NotRun` requires `runner.not_run`.
- Every other permitted non-pass state requires its matching base diagnostic.
- `RefusedReasonV2` is present exactly when the state is `Refused`.
- Active `Cancelled` and `TimedOut` require their correctly typed proven-drain
  roots.
- Active `InternalError` requires a correctly typed controlled-drain root.

This base slice validates the matrix. Lifecycle slot emission, precedence, and
immediate terminal sequencing are downstream responsibilities.

### Diagnostics and repairs

`DiagnosticCodeV2` is a closed base catalog:

| Code | Stable name |
| ---: | --- |
| 1 | `case.conformance_mismatch` |
| 2 | `runner.not_run` |
| 3 | `runner.refused` |
| 4 | `runner.no_data` |
| 5 | `runner.stale` |
| 6 | `runner.environment_invalid` |
| 7 | `runner.blocked` |
| 8 | `runner.unsupported` |
| 9 | `runner.cancelled` |
| 10 | `runner.timed_out` |
| 11 | `runner.usage` |
| 12 | `runner.internal_error` |

Code 1 covers the first conformance divergence in either the comparison or
effect lane. `comparison.mismatch`, family aliases, and lane-specific
substitutes are not valid base codes.

`RetryabilityV2` is:

| Code | Variant |
| ---: | --- |
| 0 | `Never` |
| 1 | `SameInvocation` |
| 2 | `AfterInputChange` |
| 3 | `AfterEnvironmentChange` |
| 4 | `AfterPrerequisiteChange` |

`RepairActionKindV2` is:

1. `ChangeArguments`
2. `SupplyEvidence`
3. `RegenerateCanonicalEvidence`
4. `RefreshEvidence`
5. `ReduceResourceDemand`
6. `ChooseSafeArtifactDestination`
7. `RestoreLifecycle`
8. `UpdatePolicyOrCapability`
9. `RegisterMigration`
10. `RetrySameInvocation`
11. `ContactOwner`
12. `InspectRetainedArtifact`

`ActionableDiagnosticV2` has these exact ordered fields:

1. base or registered diagnostic code
2. optional registered namespace
3. retryability
4. optional typed expected value
5. optional typed observed value
6. owner
7. zero through sixteen ordered prerequisites
8. no-claim scope
9. one through sixteen contiguous, ranked `RepairActionV2` values

`RepairActionV2` has exact fields `rank: u8`, kind, validated target token,
optional typed expected value, optional typed replacement value, validated
owner token, and optional display hint. Ranks are unique and contiguous from
one. A display hint is single-line UTF-8 of at most 256 bytes, contains no NUL
or control character, and is never executable. There is no shell, script,
callback, URI-launch, or command-string field. Rendering is generated from the
structured fields.

A repair action has a maximum canonical size of 1024 bytes. A complete
actionable diagnostic has a maximum canonical size of 8192 bytes and must fit
its enclosing record, case, run, stdout, and stderr grants. A large
`Text`/`OpaqueBytes`/detail value that would exceed the total is retained as a
bounded `DiagnosticLog` artifact and represented inline only by its typed
content root and length plus an `InspectRetainedArtifact` repair. Validation
may not produce a field-valid but unencodable mandatory diagnostic.

`DecisionDetailNamespaceV2` is a sealed `u16` registry. Family details can
refine a diagnostic only inside their registered namespace; they cannot add
terminal states, refusal reasons, precedence, or authority.

### NotRun basis

`NotRunCauseV2` has exactly three `u16` cases:

1. `PriorCancelled`, carrying `CancelledStopRootV2`
2. `PriorTimedOut`, carrying `TimedOutStopRootV2`
3. `PriorControlledInternalError`, carrying `DrainedInternalErrorRootV2`

Every value also binds the lowest remaining manifest ordinal as a `u32`.
The causal root must have the exact nominal type for its cause. There is no
generic unattempted, filtered, skipped, profile-filter, or catch-all cause.
This slice validates the basis and ordinal arithmetic only; it does not emit
suppressed lifecycle slots.

### Typed values and absence

`TypedValueV2` has exactly these variants:

`I8`, `I16`, `I32`, `I64`, `I128`, `U8`, `U16`, `U32`, `U64`, `U128`,
`Rational`, `Decimal`, `F32Bits`, `F64Bits`, `Digest`, `Quantity`, `Token`,
`Text`, `RelativePath`, and `OpaqueBytes`.

Their exact outer `u16` tags are:

| Tag | Variant | Tag | Variant |
| ---: | --- | ---: | --- |
| 1 | `I8` | 11 | `Rational` |
| 2 | `I16` | 12 | `Decimal` |
| 3 | `I32` | 13 | `F32Bits` |
| 4 | `I64` | 14 | `F64Bits` |
| 5 | `I128` | 15 | `Digest` |
| 6 | `U8` | 16 | `Quantity` |
| 7 | `U16` | 17 | `Token` |
| 8 | `U32` | 18 | `Text` |
| 9 | `U64` | 19 | `RelativePath` |
| 10 | `U128` | 20 | `OpaqueBytes` |

`NumericValueV2`, used only inside `Quantity`, is deliberately nonrecursive
and uses the same tags 1 through 14 for the corresponding numeric variants.

Absence is a separate Option-like wire sum. It is never a `TypedValueV2`
variant, an empty value, or an all-zero sentinel. Its tags are exactly
`0 Absent` and `1 Present`.

Canonical scalar rules are:

- A rational is exactly an `i128` numerator and a nonzero positive `u128`
  denominator. Sign appears only in the numerator and the pair is reduced to
  lowest terms.
- A decimal is exactly an `i128` coefficient times ten to the negative `i32`
  scale. Scale is restricted to `-6144..=6144`. Zero is coefficient zero and
  scale zero; nonzero values remove every removable trailing decimal zero.
  There is no negative zero.
- `F32Bits` and `F64Bits` preserve the exact raw bits, including signed zero,
  subnormals, infinities, and NaN payloads. They do not expose a silent total
  numerical ordering.
- A quantity is one numeric typed value plus one `UnitV2`.
- `UnitV2` has a normalized `i128/u128` rational scale and seven `i16` SI
  base-dimension exponents in this exact order: length, mass, time, electric
  current, thermodynamic temperature, amount of substance, and luminous
  intensity.
- Stable tokens are nonempty lowercase ASCII, at most 128 bytes, with nonempty
  alphanumeric segments separated only by `.`, `_`, or `-`.
- Text is bounded UTF-8.
- Logical bundle-relative paths use `/`, are at most 240 bytes and 32
  components, and reject empty, `.`, `..`, absolute, drive, NUL, and backslash
  components.
- `OpaqueBytes` is capped at 8192 bytes.
- Integer decimal renderings have at most 39 digits excluding a sign.
- Raw presented bytes cannot enter token, text, unit, path, or other validated
  APIs without a checked constructor.

### Digests, logical units, artifact roles, and axes

`DigestRoleV2` is exact:

`Spec`, `Invocation`, `Run`, `Source`, `Build`, `Toolchain`, `CaseManifest`,
`ArtifactEncoded`, `ArtifactContent`, `StoredObject`, `ArtifactInventory`,
`LifecycleLog`, `RunSummary`, `RunTerminal`, `BundleManifest`,
`DurablePublication`, `Seal`, `PublishedBundleReceipt`, `Policy`,
`CandidateBytes`, `CandidateSchema`, `SourceClosure`, `ClaimScope`,
`ProducerManifest`, and `RegisteredFamilyDomain`.

Their exact `u16` tags are:

| Tag | Role | Tag | Role |
| ---: | --- | ---: | --- |
| 1 | `Spec` | 14 | `RunTerminal` |
| 2 | `Invocation` | 15 | `BundleManifest` |
| 3 | `Run` | 16 | `DurablePublication` |
| 4 | `Source` | 17 | `Seal` |
| 5 | `Build` | 18 | `PublishedBundleReceipt` |
| 6 | `Toolchain` | 19 | `Policy` |
| 7 | `CaseManifest` | 20 | `CandidateBytes` |
| 8 | `ArtifactEncoded` | 21 | `CandidateSchema` |
| 9 | `ArtifactContent` | 22 | `SourceClosure` |
| 10 | `StoredObject` | 23 | `ClaimScope` |
| 11 | `ArtifactInventory` | 24 | `ProducerManifest` |
| 12 | `LifecycleLog` | 25 | `RegisteredFamilyDomain` |
| 13 | `RunSummary` |  |  |

A digest binds its role, registered domain, and exactly 32 bytes. Present
all-zero bytes are syntactically valid and remain distinct from absence. This
contract does not require or claim discovery of a BLAKE3 preimage whose digest
is zero.

`LogicalUnitV2` is exact:

`EncodedBytes`, `ExpandedBytes`, `StoredBytes`, `LogicalBytes`, `Count`,
`Records`, `Rows`, `Elements`, `Samples`, `Iterations`, `Operations`, `Cycles`,
`Nanoseconds`, `Seconds`, `Dimensionless`, and `RegisteredUnit`.

Their exact outer `u16` tags are:

| Tag | Unit | Tag | Unit |
| ---: | --- | ---: | --- |
| 1 | `EncodedBytes` | 9 | `Samples` |
| 2 | `ExpandedBytes` | 10 | `Iterations` |
| 3 | `StoredBytes` | 11 | `Operations` |
| 4 | `LogicalBytes` | 12 | `Cycles` |
| 5 | `Count` | 13 | `Nanoseconds` |
| 6 | `Records` | 14 | `Seconds` |
| 7 | `Rows` | 15 | `Dimensionless` |
| 8 | `Elements` | 16 | `RegisteredUnit(u16)` |

The tag-16 payload is one nonzero sealed registered-unit ID. Tag 16 without an
ID, an ID supplied to a fixed variant, or ID zero refuses.

`ArtifactRoleV2` is:

1. `Observation`
2. `ComparisonDetail`
3. `EffectDetail`
4. `DiagnosticLog`
5. `FamilyEvidence`
6. `PerformanceEvidence`
7. `ReplaySupport`
8. `RegisteredFamilyRole`, carrying a sealed family role ID

`LogicalExtentAxisV2` is:

1. `Payload`
2. `Records`
3. `Rows`
4. `Elements`
5. `Samples`
6. `Iterations`
7. `Operations`
8. `Cycles`
9. `Duration`
10. `RegisteredAxis`

Registered family domains, root policies, diagnostics, roles, units, schemas,
and executable descriptors are bounded sealed data. They are not callbacks or
runtime schema implementations.

### Wire widths

Unless a field explicitly says otherwise:

- closed catalog discriminants and registry IDs are `u16`;
- case ordinals and counts are `u32`;
- encoded, stored, expanded, and byte-budget fields are `u64`;
- logical extents and logical-work quantities are `u128`;
- rational components are `i128/u128`;
- decimal components are `i128/i32`;
- unit exponents are seven ordered `i16` values.

Widths, signs, and canonical endian/frame rules are part of wire V1. A
width-changing or sign-changing reinterpretation is not compatible wire V1.

## Runner limits

`RunnerLimitsV2` is the sole owner of the base ceilings below. Families may
tighten a base ceiling down to its documented structural minimum, but never
widen it.

### Invocation, lifecycle, and child output

- 64 argv tokens
- 8192 bytes per argv token
- 64 KiB aggregate argv
- 16 KiB per lifecycle record
- 256 records and 256 KiB per case
- 254 family rows per case
- 256 cases per invocation
- 4096 records and 4 MiB per lifecycle document
- 5 MiB complete canonical command-result stdout
- 4 MiB per child stdout
- 16 MiB Smoke or 128 MiB Full combined child stdout
- 64 KiB per child stderr
- 256 KiB combined child stderr
- 1 MiB discarded overflow drain per child stream
- 1 MiB per manifest
- nesting depth 32
- 256 comparison or effect nodes
- 512 expression edges
- 4096 memoized evaluation visits

### Values, collections, diagnostics, and registries

- 256 map entries
- 4096 generic array items
- 32 path segments
- 39 integer digits
- 16 bytes per rational numerator or denominator
- 16 bytes per decimal coefficient
- absolute decimal scale 6144
- 16 logical extents per artifact
- 256 observation keys per case
- 8 KiB text
- 128-byte stable tokens
- 240-byte logical bundle-relative paths
- 8192 opaque bytes per value
- 32 diagnostics per case and 256 per run
- 16 prerequisites and 16 repairs per diagnostic
- 64 modes per family
- 256 extension diagnostics per family
- 64 artifact roles per family
- 64 root policies per family
- 64 registered units per family
- 64 digest domains per family
- 64 decision-detail namespaces per family
- 64 output classes per family
- 64 extension schemas per family
- 64 executable descriptors per family
- 65,536 retained unknown-extension bytes per document

### Artifacts and publication

- 256 artifacts
- 64 MiB encoded per artifact
- 64 MiB expanded per artifact
- 64 MiB stored per Posix or Windows artifact
- 64 MiB plus exactly at most 4 KiB canonical envelope overhead per
  ContentStore artifact
- 64 MiB Smoke or 512 MiB Full encoded bundle total
- 64 MiB Smoke or 512 MiB Full expanded bundle total
- 65 MiB Smoke or 513 MiB Full artifact-stored aggregate
- 8 MiB system-publication stored total in either profile
- 73 MiB Smoke or 521 MiB Full whole-publication stored total

For Posix and Windows, stored length equals encoded length. For ContentStore,
stored length equals encoded length plus the exact canonical envelope
non-payload length. Sums include complete stored objects, use checked `u64`
arithmetic, and refuse overflow before allocation.

`system_publication_stored_bytes` is the checked sum of the six logical
`runner/` objects later frozen by `BundleLayoutV1`: lifecycle log, run
terminal, artifact inventory, bundle manifest, publication intent, and seal.
Artifact inventory entries remain artifact accounting. Command results,
physical transaction residue, store generations, commit-service metadata, and
recovery locators are not logical publication objects.

`publication_stored_bytes` equals checked
`artifact_stored_bytes + system_publication_stored_bytes`. No object may be
omitted or counted twice. Writers, verifiers, replay, and deterministic logs
report artifact-stored, system-stored, and whole-publication totals separately.

### Structural feasibility

A sealed executable family has at least one mode, one case, one executable
descriptor, one comparison root, and one effect root. Expression-node and edge
caps admit the declared roots. Per-case lifecycle capacity admits at least
`CaseStart` and `CaseTerminal`.

Run lifecycle feasibility is checked as:

```text
3 + sum_for_each_case(2 + declared_maximum_family_rows)
```

The leading three records are `RunStart`, `RunSummary`, and `RunTerminal`.
This result must fit the declared run record capacity and the base maximum of
4096. The independent 256-case, 256-record-per-case, and 254-family-row
ceilings are not blindly multiplied. Zero is legal only for genuinely optional
extension, artifact, unknown-data, or family-row capacity. A jointly
infeasible limit vector refuses before proportional allocation.

Every cap carries its exact encoded, stored, expanded, logical, node, item, or
byte unit. No downstream parser may invent or silently derive another bound.

## Runner budgets

`RunnerBudgetsV2` has this frozen field order:

1. `wall_time_ns`
2. `max_resident_bytes`
3. `max_child_processes`
4. `max_parallel_children`
5. `logical_work_limit`
6. `logical_work_unit`
7. `lifecycle_encoded_bytes`
8. `command_result_stdout_bytes`
9. `combined_child_stdout_bytes`
10. `combined_child_stderr_bytes`
11. `artifact_encoded_bytes`
12. `artifact_stored_bytes`
13. `artifact_expanded_bytes`
14. `system_publication_stored_bytes`
15. `publication_stored_bytes`
16. `stop_observation_ns`
17. `drain_ns`
18. `finalize_ns`

All fields are semantic invocation inputs. Wall time, resident memory,
stop-observation, drain, finalize, and command-result stdout grants are
nonzero. Total children may be zero only for
`LifecycleOnlyNoBundle`; parallel children never exceed total children.
The checked sum of stop-observation, drain, and finalize allowances never
exceeds wall time.

Every output grant is at or below `RunnerLimitsV2`. Logical work is an exact
`u128` paired with one `LogicalUnitV2`.

Smoke admits at most:

- 900 seconds wall time
- 16 GiB resident memory
- 32 parallel children
- 256 total children

Full admits at most:

- 86,400 seconds wall time
- 128 GiB resident memory
- 64 parallel children
- 256 total children

These are admission ceilings, not performance claims.

The command-result grant caps the complete canonical `List`,
`LifecycleOnlyCommandResultV2`, or `DurableCommandResultV2` frame, including
framing, embedded lifecycle document, and published-bundle receipt when
present. The base cap is 5 MiB; lifecycle remains at most 4 MiB, and
`RunnerCatalogV2` and each embedded receipt are at most 1 MiB. Joint
feasibility is checked. CLI stderr is zero on success or one
`ActionableDiagnosticV2` in a canonical failure frame of at most 16 KiB.
Durable lifecycle staging consumes the same lifecycle and resident-memory
budgets and cannot create an unbudgeted second copy.

Zero system-publication grant is legal only for `List` or
`LifecycleOnlyNoBundle` paths that write no durable publication.

## Publication selection and root-capability policy

`PublicationProtocolV2` is exact:

1. `PosixDescriptorRenameAndDirectorySyncV1`
2. `WindowsHandleReplaceAndDirectoryFlushV1`
3. `ContentStoreAtomicCommitV1`

`DestinationAdmissionModeV2` is exactly `1 Absent` or
`2 PreExistingEmpty`. ContentStore code 2 additionally requires an exclusive,
leased-empty namespace.

`LogicalBundlePathV1` and `ContentStoreObjectKeyV1` are distinct validated
types. Neither contains a physical capability prefix, transaction key,
generation, credential, attempt locator, or acquired handle.

`PublicationTargetV2` is:

- `PosixRelative(LogicalBundlePathV1)`
- `WindowsRelative(LogicalBundlePathV1)`
- `ContentStoreLogicalKey(ContentStoreObjectKeyV1)`

`PublicationSelectionV2` binds, in order:

1. `PlatformPathProfileV2`
2. `PublicationProtocolV2`
3. `DestinationAdmissionModeV2`
4. `PublicationTargetV2`

Durable commands require a selection. `List`, `Check`, and `SelfTest` carry
typed absence. A selection is semantic intent, not evidence that a destination
exists, is empty, is safe, was acquired, or was durably committed.

`RootCapabilityAccessV2` is `1 ReadOnlyInput` or `2 DurableOutput`.

`RootCapabilityRightV2` is:

1. `Traverse`
2. `ReadObject`
3. `Enumerate`
4. `CreateObject`
5. `PopulateEmptyDestination`
6. `SyncObject`
7. `SyncContainer`
8. `AcquireExclusiveLease`
9. `QueryGeneration`
10. `CommitCompareAndSwap`

`RootCapabilityPolicyV2` has exact ordered fields:

1. root class
2. platform path profile
3. access
4. canonical nonempty right set
5. `freshness_policy_id: u16`
6. `revocation_policy_id: u16`
7. `overlap_policy_id: u16`
8. no-claim scope

Its domain is:

```text
org.frankensim.fs-evidence-runner.root-capability-policy.v1
```

Its semantic root constructor is private.

The exact least-privilege right sets are:

- Posix or Windows read-only input: `Traverse`, `ReadObject`, `Enumerate`.
- Posix or Windows durable output: `Traverse`, `Enumerate`, `CreateObject`,
  `SyncObject`, `SyncContainer`; add `PopulateEmptyDestination` only for
  `PreExistingEmpty`.
- ContentStore read-only input: `ReadObject`, `Enumerate`,
  `QueryGeneration`.
- ContentStore durable output with `Absent`: `CreateObject`,
  `QueryGeneration`, `CommitCompareAndSwap`.
- ContentStore durable output with `PreExistingEmpty`: the `Absent` set plus
  `Enumerate` and `AcquireExclusiveLease`.

A policy cannot omit or add a right for its protocol/access/mode cell. A
broader physical capability may be accepted later only through an affine,
narrowed view exposing exactly the policy set. Neither the broader capability
nor its physical details enter semantic identity or deterministic logs.

Invocation and bundle schemas later bind an ordered, duplicate-rejected set of
root-capability policy roots sorted by embedded root class. `Run` and
`Negative` require one durable output policy. `Replay` requires one read-only
input and one durable output policy.

Opaque policy IDs are declarations, not evidence. Policy validation therefore
has three distinct stages:

1. intrinsic field, canonical-right, and nonzero-ID validation;
2. exact least-privilege validation against `PublicationSelectionV2`;
3. registration validation against a bounded, duplicate-rejected
   `RootPolicyRegistryProjectionV2`.

The projection exact-set lists registered freshness and revocation policy IDs
and registered overlap policy IDs. The only base overlap relation is
`1 RequireInputOutputDisjoint`. Replay input and output policies must carry the
same registered nonzero overlap-policy ID and that registration must declare
`RequireInputOutputDisjoint`. Distinct IDs, merely nonzero IDs, or absent
registry data refuse.

This establishes only declaration consistency. The downstream capability
broker/provider must separately prove that the two acquired resources are
actually disjoint before any write. Equal, aliased, substituted, or
unadjudicable resources refuse without effects. Other roots require a
registered `RootClassV2` policy ID.

Publication selection does not duplicate budgets or capability policy.
Descriptors, handles, broker slots, credentials, namespace prefixes,
generations, and attempts remain nonsemantic observations acquired and checked
by the downstream capability owner.

## Presented nominal roots

This base schema owns private-field, role-specific presented reference wrappers
and their `.v1` descriptors for:

- `SourceIdentityRootV2` — `source-identity`
- `BuildIdentityRootV2` — `build-identity`
- `ToolchainIdentityRootV2` — `toolchain-identity`
- `CaseManifestRootV2` — `case-manifest`
- `ArtifactEncodedRootV2` — `artifact-encoded`
- `ArtifactContentRootV2` — `artifact-content`
- `StoredObjectRootV2` — `stored-object`
- `ArtifactInventoryRootV2` — `artifact-inventory`
- `LifecycleLogRootV2` — `lifecycle-log`
- `RunSummaryRootV2` — `run-summary`
- `RunTerminalRecordRootV2` — `run-terminal-record`
- `BundleManifestRootV2` — `bundle-manifest`
- `PresentedPublicationCommitRefV2` — `presented-publication-commit-ref`
- `DurablePublicationIdentityV2` — `durable-publication-identity`
- `SealRootV2` — `seal`
- `PublishedBundleReceiptRootV2` — `published-bundle-receipt`
- `AuthorityScopeRootV2` — `authority-scope`
- `ExternalMutationSetRootV2` — `external-mutation-set`
- `ArtifactSetRootV2` — `artifact-set`
- `ResourceIdentityRootV2` — `resource-identity`
- `RunnerLimitsSchemaRootV2` — `runner-limits-schema`
- `RunnerLimitsRootV2` — `runner-limits`
- `RunnerBudgetsSchemaRootV2` — `runner-budgets-schema`
- `RunnerBudgetsRootV2` — `runner-budgets`
- `RootCapabilityPolicyRootV2` — `root-capability-policy`
- `NoClaimScopeRootV1` — `no-claim-scope`
- `CancelledStopRootV2` — `cancelled-stop`
- `TimedOutStopRootV2` — `timed-out-stop`
- `DrainedInternalErrorRootV2` — `drained-internal-error`

A checked parser validates the registered domain, expected role, exact 32-byte
shape, and canonical lowercase textual form. A presented wrapper proves no
existence, byte possession, content equivalence, lifecycle completion,
durability, verification, admission, or authority.

There is no generic root wrapper, cross-role conversion, semantic identity
constructor, lifecycle-root constructor, durability constructor, or seal
constructor in this base slice. The limits, budgets, and root-capability
modules own their private semantic root constructors. The three stop wrappers
are presented-only here; lifecycle phase 24.1.2 remains their sole semantic
stop/drain constructor owner.

## Invariants

- Every closed catalog has one exact code-to-variant mapping. Unknown codes
  refuse; aliases and best-effort decoding do not exist.
- API generation 2 and wire version 1 remain separately typed and
  non-interchangeable.
- Wire V1 has no predecessor or implicit migration.
- Canonical values have exactly one admitted representation except where raw
  IEEE bits intentionally preserve distinct encodings.
- Presence is explicit and never inferred from zero, empty bytes, or an
  all-zero digest.
- Counts, lengths, arithmetic relations, and structural minima are checked
  before proportional allocation.
- Every family limit is equal to or tighter than the base limit and remains
  jointly feasible.
- State, reason, record role, diagnostic, and causal-root relationships obey
  the closed Cartesian matrix.
- A base diagnostic cannot be replaced by a family code. A family namespace
  cannot add a terminal state or refusal reason.
- Repairs are structured, bounded, ranked, non-executable data.
- Logical paths and object keys remain distinct from physical capability
  locators.
- Publication selection and root-capability policy are semantic inputs, not
  proof that physical capabilities or durable effects exist.
- Presented identities are nominal, role-specific references and cannot mint
  semantic or authority-bearing identities.
- Canonical identity input changes move the corresponding identity; variation
  in explicitly nonsemantic physical acquisition does not.
- Pure validation does not write, spawn, retain ambient state, consult an
  environment-dependent default, or emit lifecycle or publication artifacts.

## Error model

All constructors and validators are fail-closed and return typed, deterministic
errors. Validation never repairs, truncates, aliases, coerces, defaults, or
normalizes malformed external input silently.

Errors distinguish at least:

- unknown or disallowed catalog and registry codes;
- wrong role, domain, width, sign, or nominal root type;
- noncanonical rational, decimal, token, text, path, unit, or collection data;
- wrong state/reason/role/diagnostic combinations;
- missing, extra, duplicated, misordered, or incompatible command selectors;
- limit, budget, structural-minimum, joint-feasibility, and checked-arithmetic
  failures;
- wrong or noncontiguous repair ranks and incompatible replacement types;
- invalid publication protocol/profile/mode/target combinations;
- missing, extra, duplicated, or substituted capability rights;
- unsupported platform profiles;
- attempts to cross an authority boundary.

An actionable failure names the precise dimension, typed unit, expected bound
or value, observed value, owner, ordered prerequisites, no-claim scope, and
ranked structured repairs when those data are applicable. An error is
diagnostic data only; it is not evidence of scientific verification,
admission, promotion, or durable execution.

Pure validators do not panic for caller-controlled values. Integer and size
arithmetic is checked. Refusal happens before allocation proportional to a
rejected count or length.

## Determinism class

All base-schema validation is deterministic for the same complete typed input:

- catalogs, registries, and ordered field tables have fixed order;
- sets use their declared canonical ordering and reject duplicates;
- error precedence and the first reported divergent dimension are stable;
- rendering is generated from structured data and has canonical field order;
- no locale, wall clock, scheduler, process identifier, absolute path,
  filesystem enumeration, environment variable, or randomized hash order
  influences a result;
- validation results are independent of worker count because this slice has no
  parallel execution.

Determinism of a pure schema result is not a performance, scientific, process,
storage, lifecycle, or durability claim.

## Cancellation behavior

Cancellation is intentionally not part of this base slice. These APIs perform
only bounded, in-memory validation under `RunnerLimitsV2`; they do not wait on
I/O, spawn work, acquire a lease, or enter a lifecycle operation. Once invoked,
a pure validator runs to completion and returns success or a typed refusal.

This non-applicability must not be generalized to later lifecycle, child
process, filesystem, ContentStore, verification, or publication work. Those
owners must use explicit asupersync scopes, bounded stop observation,
request-drain-finalize sequencing, and typed drained roots.

## Unsafe boundary

The crate is pure safe Rust and forbids unsafe code. It has no FFI, dynamic
loading, build script, proc macro, platform syscall, memory mapping, raw
descriptor, raw handle, or credential implementation.

## Feature flags

None. No feature may widen authority, relax canonical validation, enable a
mock production path, or alter wire semantics.

## Dependency and source policy

Phase-1 normal code has exactly one direct dependency:

```toml
fs-blake3 = { path = "../fs-blake3" }
```

It is used only for canonical hashing and nominal identity construction. This
phase has no build dependencies, development dependencies, optional
dependencies, renamed dependencies, target-specific dependencies, ambient
registry sources, or feature routes.

The closed eventual direct normal-dependency allowlist for the complete TOOL
crate is exactly:

- `fs-blake3`
- `asupersync`
- `fsqlite` from the pinned sibling source with exactly its `async-api`
  feature

The latter two rows belong exclusively to the later capability, cancellation,
process, and persistent-ContentStore phase. Their future addition cannot
invalidate the Phase-1 proof that base-schema modules import only `fs-blake3`.
No build, proc-macro, optional, renamed, target-split, registry-substituted, or
unowned dependency route is admitted. Direct and transitive feature/source
identity must unify with the pinned constellation.

The crate imports, reexports, and constructs no promotion,
verification-admission, scientific-authority, or rjoq type.

Source-closure evidence exact-set binds every file and generated input owned or
consumed by a leaf: implementation modules, this contract and rustdoc,
schema/domain/descriptor and constructor tables, independent oracles, unit and
property fixtures, mutation fixtures, in-process E2E projection rows, log
schemas, direct dependency/source/feature identities, and build configuration.
Missing, extra, stale, reordered, ambient, unpinned, duplicate-owner, or
mixed-snapshot inputs move the closure or refuse. A retained binary or root
cannot substitute for absent source.

## Retention and logging

The production base library prints nothing and performs no ambient logging.
Its deterministic test harness and retained conformance evidence record:

- API, wire, source, build, toolchain, target, and feature roots;
- exact catalog and cap counts;
- state/reason/role/diagnostic matrix cells;
- typed expected and observed values;
- causal roots and manifest ordinals;
- owner, prerequisites, no-claim scope, and repairs;
- minimized relative case names and relative artifact paths;
- separately typed encoded, stored, expanded, logical, cycle, duration, and
  performance quantities where those are semantic inputs.

Logs do not retain absolute paths, PIDs, wall-clock timestamps, scheduler
latency, environment secrets, credentials, physical capability locators, raw
bulk payloads, or undeclared timing. Raw rejected payloads are represented by
bounded typed metadata or retained relative artifact references.

Coverage manifests expose nonzero eligible, passed, failed, and unsupported
counts as applicable and keep focused Phase-1 proof distinct from repository
policy, DSR, release-built E2E, and downstream authority proof.

## Conformance tests

The base slice requires all of the following independent evidence classes:

1. handwritten literal-oracle tests for every catalog, code, version,
   predecessor rule, role, unit, limit, and ordered schema field;
2. exhaustive Cartesian state/reason/role/diagnostic and command-applicability
   matrices;
3. zero, one-below, exact, one-above, extrema, checked-overflow, and structural
   feasibility boundary tests;
4. property and metamorphic tests for canonical equality, one-field identity
   movement, canonical ordering, and limit tightening;
5. public-API tests, compile-fail tests, and doctests proving nominality,
   privacy, raw-versus-validated separation, and the absence of coercion or
   authority surfaces;
6. independent wire-width, schema, domain, endian/frame, mutation, and unknown
   code tests;
7. no-mock integration tests using real public constructors and validators;
8. deterministic structured-logging tests;
9. exact source-closure and coverage-manifest tests;
10. in-process, source-closed Runner V2 base E2E projection tests.

Compile-fail cases have a corresponding safe positive constructor or
read-only-accessor example. The no-mock fixture emits no lifecycle, bundle,
scientific receipt, admission receipt, or promotion witness.

Focused crate tests establish only base-schema construction and validation.
Repository-wide dependency policy, contract scanning, documentation facts,
source/claim inventories, generated grammar, DSR quality, release-built E2E,
and downstream authority remain separate proof obligations.

## In-process E2E projection handoff

This leaf owns a source-closed `RunnerV2BaseE2eProjectionV1` test manifest.
It contains exact journey-keyed projections for:

- `scripts/ci/e2e_evidence_runner_publication_state_v2.sh`
- `scripts/ci/e2e_evidence_runner_publication_v2.sh`
- `scripts/ci/e2e_evidence_verifier_v1.sh`
- `scripts/ci/canonical_evidence_runner_v2.sh`
- `scripts/ci/verify_runner_rjoq_handoff_v1.sh`

The rows cover every base literal, limit, budget, unit, capability-policy,
one-field identity mutation, exact diagnostic, and no-claim boundary consumed
by each journey. This leaf runs every eligible projection row **in process**
through real public constructors and validators with no mocks. It emits
deterministic detailed logs and reports eligible, passed, failed, unsupported,
projection-E2E, and logging counts.

This leaf does not create, edit, invoke, or claim release-built execution of
those scripts. Each sole downstream script owner must exact-set consume its
corresponding immutable projection root, reject missing, extra, duplicate,
stale, or unmapped rows, and execute every eligible row. Final program closure
requires all five release-built projections; their execution is not a
prerequisite for closing this pure base-schema leaf.

## Downstream ownership exclusions

The following responsibilities are explicitly outside this base slice:

- RunnerSpec semantic construction and sealed family/registry membership:
  `frankensim-epic-foundations-huq.24.1.1.3` and its designated registry
  leaves.
- Invocation and run semantic identity constructors:
  `frankensim-epic-foundations-huq.24.1.1.2`.
- Artifact registration, artifact inventory, and bundle-manifest semantic
  schemas: `frankensim-epic-foundations-huq.24.1.1.4` and their designated
  later owners.
- Lifecycle record emission, slot sequencing, precedence, lifecycle roots,
  cancellation draining, and terminal publication:
  `frankensim-epic-foundations-huq.24.1.2`.
- Hostile-byte framing, parsing, codec enforcement, and migration:
  `frankensim-epic-foundations-huq.24.2.1`.
- Physical capability acquisition and affine narrowing, process execution,
  filesystem operations, persistent ContentStore operations, durability,
  recovery, publication commit, seal construction, and published receipt
  construction: `frankensim-epic-foundations-huq.24.2.2`.
- Verifier execution and independent evidence checking:
  `frankensim-epic-foundations-huq.24.3`.
- Authority coherence, admission, promotion, and rjoq handoff:
  `frankensim-epic-foundations-huq.24.5`.
- TOOL dependency policy, specialized contract checks, source/claim scans,
  documentation facts, count gates, generated retained inventories, and CI/DSR
  wiring: `frankensim-epic-foundations-huq.24.1.3.1`.
- Root-free ownership registries, generated grammar, coverage joins, and exact
  registry inventory: `frankensim-epic-foundations-huq.24.1.3.2`.

The sole downstream owners of the five release-built scripts are respectively
`frankensim-epic-foundations-huq.24.2.2.2`,
`frankensim-epic-foundations-huq.24.2.2.3`,
`frankensim-epic-foundations-huq.24.3.3.3`,
`frankensim-epic-foundations-huq.24.4.1.3`, and
`frankensim-epic-foundations-huq.24.5.3.1`. This crate provides immutable
projection inputs to those owners; it does not duplicate their script
ownership.

## No-claim boundaries

Acceptance by this crate means only that base Runner V2 data are
well-formed, canonical, within their declared limits and budgets, and
consistent with the closed pure-validation matrices.

It does **not** claim:

- that a source, build, toolchain, case, artifact, object, log, manifest,
  publication, receipt, policy, or resource exists;
- that presented bytes have the asserted meaning or match retained content;
- that a case ran, a lifecycle completed, cancellation drained, or a process
  terminated;
- that a path or logical key maps to a safe physical destination;
- that a capability was acquired, remained fresh, was unrevoked, or enforced
  a policy;
- that an artifact or bundle was stored, synchronized, committed, sealed,
  recoverable, durable, or published;
- that evidence is scientifically correct, validated, independently verified,
  admissible, or sufficient for promotion;
- that an admission, waiver, release, rjoq, or other authority boundary was
  crossed;
- that any performance, resource-use, platform, or timing target was met;
- that repository-wide policy, DSR, or release-built E2E proof passed.

Presented roots, logical publication choices, policy roots, diagnostics,
coverage manifests, and in-process E2E projection results remain
non-authoritative data. Scientific, durability, verification, admission,
promotion, release, and rjoq authority can arise only from their separately
owned, explicitly capability-bearing protocols and retained evidence.
