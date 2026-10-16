# CONTRACT: fs-evidence-runner

## Purpose and layer

Layer: **TOOL**.

`fs-evidence-runner` owns the frozen Runner V2 orchestration vocabulary and its
bounded, pure validation rules. The base schema makes command intent, terminal
state, diagnostics, limits, budgets, typed values, logical publication
selection, capability policy, and presented identity references explicit
before lifecycle, process, storage, verification, or admission implementations
are allowed to consume them.

This base family is declarations and bounded, in-memory validation, including
the explicitly named pure Stage-A local evaluators. It does not execute an
external, lifecycle, or release-built case, emit lifecycle records, parse
hostile bytes, access a filesystem, spawn a process, persist an object, publish
a bundle, verify scientific evidence, or mint authority.

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

Non-wire semantic projections also carry an explicit, independently versioned
domain. Their first frozen form uses `.v1`, but that suffix describes the
projection schema rather than claiming that the projection itself is a wire
frame. In particular, the initial extension registry, logical extent,
constructor-owner handoff, and root-free evaluator-member guard projections
all use distinct `.v1` domains.

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
`CommandSelectorPresenceV2` represents all five caller-selectable dimensions
with exact `Absent`, `Singular`, `Duplicate`, or `Ambiguous` cardinality.
`validate_command_selector_presence_v2` checks the frozen order family, mode,
profile, negative case, then replay source and returns the first exact
`CommandSelectorUsageV2`: `RunnerUsage`, selector field, command-specific
`Absent` or `Singular` expectation, observed cardinality, and stable owner. It
does not parse tokens, infer defaults, or treat a sealed manifest-derived value
as a caller selector.

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

### Fixed artifact codecs and logical extents

`ArtifactCodecIdV2` is a closed declaration-only catalog:

| Code | Variant |
| ---: | --- |
| 0 | `Identity` |
| 1 | `ZstdFrameV1` |

The fixed codec rows consume no family registry slot and expose no encoder,
decoder, frame parser, checksum verifier, migration, or artifact-inventory
behavior. Unknown codec codes refuse.

`LogicalExtentV2` has private fields and the exact semantic order:

1. `axis`, including its tag and optional registered ID;
2. `value: u128`;
3. `unit`, including its tag and optional registered ID.

Its first non-wire semantic domain is
`org.frankensim.fs-evidence-runner.logical-extent-projection.v1`.
`LogicalExtentFieldV1` freezes the same three-row order. The semantic root
binds every field and distinguishes zero, one, and `u128::MAX`. A bare
syntactic registered axis or unit cannot enter the registry-free constructor.

The exact base axis/unit cells are:

| Axis | Admitted unit |
| --- | --- |
| `Payload` | `LogicalBytes` |
| `Records` | `Records` |
| `Rows` | `Rows` |
| `Elements` | `Elements` |
| `Samples` | `Samples` |
| `Iterations` | `Iterations` |
| `Operations` | `Operations` |
| `Cycles` | `Cycles` |
| `Duration` | `Nanoseconds` |
| `Duration` | `Seconds` |

`Nanoseconds` is the canonical Duration unit. One second scales to exactly
`1_000_000_000/1` nanoseconds. Cross-axis units refuse.

### Base extension registry and exact conversions

`BaseExtensionRegistryProjectionV2` contains exactly three independent typed
categories:

1. registered artifact-role descriptors;
2. registered logical-unit descriptors;
3. registered logical-extent-axis descriptors.

Its first non-wire domain is
`org.frankensim.fs-evidence-runner.base-extension-registry-projection.v1`.
Every projection binds the exact `RunnerLimitsV2` semantic root. IDs are
nonzero and namespace-local, so the same numeric ID may appear once in each
different category. Duplicate IDs within one category refuse. Names are
globally namespaced and unique across all three categories. Every descriptor
binds its name, owner, no-claim scope, and category-specific fields.

Registered axes have a nonempty canonical allowed-unit set. The canonical unit
occurs exactly once with scale `1/1`; every scale is positive, reduced, and
exact. Registered unit references must resolve in the same exact registry.
The axis constructor enforces an absolute 4,096-row pre-allocation ceiling,
while the enclosing projection additionally enforces its admitted generic
array ceiling. Registry categories independently admit zero, one, and their
exact 64-row base ceilings and refuse 65. Exact reconstruction reports
category-specific missing, extra, or mutated data and is permutation-invariant
within each category.

Registry-aware extent construction proves only membership in the supplied
projection. Conversion uses exact rational arithmetic and succeeds only when
the final extent value is an integral `u128`. Unavailable units, cross-axis
substitution, fractional results, and arithmetic overflow refuse.
`normalized_unit_scale_ratio_v2` and `convert_rational_quantity_v2` separately
convert dimension-compatible `UnitV2` values through exact normalized rational
scales; dimension mismatch, noncanonical input, and checked-product overflow
refuse. None of these declarations proves execution, measurement, physical
units, scientific validity, or authority.

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

The canonical schema contains exactly 71 fields in this order. `S/F` denotes
the sole Smoke/Full profile difference; every other value is profile-equal.
`RunnerLimitDescriptorV2` additionally freezes each row's semantic unit,
tightenability, and structural-minimum rule.

| Code | Field | Width | Base value |
| ---: | --- | --- | ---: |
| 1 | `argv_tokens` | `u32` | 64 |
| 2 | `argv_token_bytes` | `u64` | 8 KiB |
| 3 | `argv_aggregate_bytes` | `u64` | 64 KiB |
| 4 | `lifecycle_record_encoded_bytes` | `u64` | 16 KiB |
| 5 | `case_lifecycle_records` | `u32` | 256 |
| 6 | `case_lifecycle_encoded_bytes` | `u64` | 256 KiB |
| 7 | `family_rows_per_case` | `u32` | 254 |
| 8 | `invocation_cases` | `u32` | 256 |
| 9 | `lifecycle_document_records` | `u32` | 4096 |
| 10 | `lifecycle_document_encoded_bytes` | `u64` | 4 MiB |
| 11 | `command_result_stdout_bytes` | `u64` | 5 MiB |
| 12 | `child_stdout_bytes` | `u64` | 4 MiB |
| 13 | `combined_child_stdout_bytes` | `u64` | 16/128 MiB S/F |
| 14 | `child_stderr_bytes` | `u64` | 64 KiB |
| 15 | `combined_child_stderr_bytes` | `u64` | 256 KiB |
| 16 | `manifest_encoded_bytes` | `u64` | 1 MiB |
| 17 | `nesting_depth` | `u32` | 32 |
| 18 | `comparison_nodes` | `u32` | 256 |
| 19 | `effect_nodes` | `u32` | 256 |
| 20 | `text_bytes` | `u64` | 8 KiB |
| 21 | `stable_token_bytes` | `u64` | 128 |
| 22 | `bundle_relative_path_bytes` | `u64` | 240 |
| 23 | `diagnostics_per_case` | `u32` | 32 |
| 24 | `diagnostics_per_run` | `u32` | 256 |
| 25 | `prerequisites_per_diagnostic` | `u32` | 16 |
| 26 | `repairs_per_diagnostic` | `u32` | 16 |
| 27 | `artifacts` | `u32` | 256 |
| 28 | `artifact_encoded_bytes` | `u64` | 64 MiB |
| 29 | `artifact_expanded_bytes` | `u64` | 64 MiB |
| 30 | `artifact_stored_bytes` | `u64` | 64 MiB + 4 KiB |
| 31 | `bundle_encoded_bytes` | `u64` | 64/512 MiB S/F |
| 32 | `bundle_expanded_bytes` | `u64` | 64/512 MiB S/F |
| 33 | `artifact_stored_aggregate_bytes` | `u64` | 65/513 MiB S/F |
| 34 | `system_publication_stored_bytes` | `u64` | 8 MiB |
| 35 | `publication_stored_bytes` | `u64` | 73/521 MiB S/F |
| 36 | `child_stream_discard_bytes` | `u64` | 1 MiB |
| 37 | `modes_per_family` | `u32` | 64 |
| 38 | `extension_diagnostics_per_family` | `u32` | 256 |
| 39 | `artifact_roles_per_family` | `u32` | 64 |
| 40 | `root_policies_per_family` | `u32` | 64 |
| 41 | `registered_units_per_family` | `u32` | 64 |
| 42 | `digest_domains_per_family` | `u32` | 64 |
| 43 | `extension_schemas_per_family` | `u32` | 64 |
| 44 | `executable_descriptors_per_family` | `u32` | 64 |
| 45 | `map_entries` | `u32` | 256 |
| 46 | `generic_array_items` | `u32` | 4096 |
| 47 | `path_segments` | `u32` | 32 |
| 48 | `integer_digits` | `u32` | 39 |
| 49 | `rational_component_bytes` | `u64` | 16 |
| 50 | `decimal_coefficient_bytes` | `u64` | 16 |
| 51 | `decimal_absolute_scale` | `u32` | 6144 |
| 52 | `logical_extents_per_artifact` | `u32` | 16 |
| 53 | `observation_keys_per_case` | `u32` | 256 |
| 54 | `decision_detail_namespaces` | `u32` | 64 |
| 55 | `output_classes` | `u32` | 64 |
| 56 | `opaque_value_bytes` | `u64` | 8192 |
| 57 | `retained_unknown_extension_bytes` | `u64` | 65,536 |
| 58 | `expression_edges` | `u32` | 512 |
| 59 | `memoized_evaluation_visits` | `u32` | 4096 |
| 60 | `repair_action_encoded_bytes` | `u64` | 1024 |
| 61 | `actionable_diagnostic_encoded_bytes` | `u64` | 8192 |
| 62 | `failure_stderr_encoded_bytes` | `u64` | 16,384 |
| 63 | `runner_catalog_encoded_bytes` | `u64` | 1 MiB |
| 64 | `published_bundle_receipt_encoded_bytes` | `u64` | 1 MiB |
| 65 | `content_store_envelope_non_payload_bytes` | `u64` | 4096 |
| 66 | `registered_extent_axes_per_family` | `u32` | 64 |
| 67 | `registered_observation_keys_per_family` | `u32` | 4096 |
| 68 | `registered_authority_scopes_per_family` | `u32` | 64 |
| 69 | `registered_external_root_classes_per_family` | `u32` | 64 |
| 70 | `registered_evaluation_units_per_family` | `u32` | 64 |
| 71 | `registered_resource_identities_per_family` | `u32` | 256 |

Codes 66 through 71 are independent profile-equal categories, may each be
tightened to zero, and cannot borrow unused capacity from another category.
Code 53 remains the per-case observation-key cap; code 67 bounds a family's
registered observation-key union. Code 55 remains output classes. Fixed codec
rows consume no registry slot.

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

`RunnerLimitsV2::admit_family` remains the compatibility fail-first boundary.
`validate_family_complete` and `admit_family_complete` additionally expose one
fixed-size report with at most one refusal per limit field. Its public iterator
is always in ascending field order, while
`compatibility_first_violation` identifies the exact refusal the fail-first
boundary would return. Same-field precedence is individual fixed/base/minimum,
then declared-minimum order/width/value in presented order, then frozen joint
nested/artifact/case feasibility order, then executable-family shape order.
The report checks every independently rejected field across those phases,
retains no caller-sized diagnostic collection, and never allocates in
proportion to the declared-minimum or per-case-row slices.

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

Relational overflow diagnostics retain the exact mathematical sum as `U128`
rather than substituting `u64::MAX` or dropping an operand. Timeout-sum
refusals target the stop-observation, drain, or finalize inputs; a
whole-publication sum overflow targets either the artifact-stored or
system-publication-stored summand. An equation mismatch may instead target the
presented `publication_stored_bytes` total. Every repair target is bounded
non-executable data.

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

## Runner V2 base schema-impact declaration

`schema_impact` owns the result-free, declaration-only schema-impact contract
for the Runner V2 base leaf. It does not decode untrusted wire data, inspect
retained artifacts, execute a migration, observe a downstream runtime, or mint
compatibility merely because two schemas have similar shapes. Hostile bounded
decoding remains owned by
`frankensim-epic-foundations-huq.24.2.1.1`; this leaf presents only values
constructed through its closed Rust API.

The contract has independently typed, closed catalogs for frame version,
authority state, slot use, manifest relation, registry kind, authority surface,
field wire kind, field layout, impact disposition, and migration policy. A
field descriptor binds its one-based ordinal, nonzero field code, stable field
name and semantic-type ID, wire kind, structural layout, optional reciprocal
field code, and optional version-slot code. A frame descriptor binds its Rust
schema name, frame version, domain, raw terminal-version magic, exact ordered
field descriptors, separately typed API generation and Runner wire version,
predecessor policy, and optional nominal-root role. Version slots, registry
fragments, impact rows, and the complete manifest each have distinct canonical
domains and nominal roots. There is no raw-digest constructor, cross-role root
conversion, positional substitution between frozen-base and leaf-extension
values, or caller-selected manifest ordinal.

The production declaration is exactly
`runner_v2_base_schema_impact_manifest_v1()`. Its three canonically ordered
rows describe `CanonicalSchemaFieldDescriptorV1`,
`CanonicalSchemaFrameDescriptorV1`, and
`CanonicalSchemaVersionSlotDescriptorV1`. They use the immutable FrozenBase
registry and no leaf-extension fragment. Separately, the common manifest API
admits a leaf extension only from a crate-owned static declaration that binds
its owner, stable fragment ID, exact FrozenBase root, exact compatible source
member, and complete descriptor sequence. Being listed by a caller is not
source closure. Mixed snapshots; unknown, missing, extra, duplicate, reordered,
substituted, or conflicting members; and an extension contributing no used
role refuse before a manifest is constructed.

Each impact row makes one of six dispositions explicit:
`NewV1NoPredecessor`, `UnchangedV1`, `MigratedV1ToV2`,
`DecodeOnlyLegacyV1`, `RetiredV1`, or
`InapplicableNoCanonicalFrame`. Its migration policy, historical/current frame
bindings, optional legacy container, predecessor relation, and authority
surfaces must satisfy the complete closed matrix. Decode-only, retired, and
inapplicable rows—and every other row without an authoritative current
frame—expose an exactly empty authority surface. No disposition or policy
silently inherits authority. The manifest derives canonical row ordinals,
validates exact row-role coverage and unique roots, resolves reciprocal
parent/child slots and legacy containers, and checks the bounded acyclic
dependency graph. A compatibility-evidence edge is rejected whenever its
fixed-point reachability includes a row with a nonempty authority surface.

All AC60 cardinality maxima and canonical-size constants are inclusive. The
canonical-size constants are defensive safety guards, not a claim that the
current nested V1 grammar can construct an admitted component whose encoding
reaches the guard:

| Surface | Inclusive maximum |
| --- | ---: |
| stable token or registered text | 128 bytes |
| source path | 240 bytes |
| fields per frame | 256 |
| authority surfaces per row | 6 |
| predecessor, parent, or child slots | 256 |
| impact rows per manifest | 256 |
| dependency edges per manifest | 512 |
| registered roles | 64 |
| extension registry fragments | 256 |
| canonical field descriptor | 1,024 bytes |
| canonical version slot | 2,048 bytes |
| canonical frame descriptor | 262,144 bytes |
| canonical registry fragment | 65,536 bytes |
| canonical impact row | 1,048,576 bytes |
| canonical schema-impact manifest | 1,048,576 bytes |

The independently derived tight reachability bounds for the present V1 grammar
are:

| V1 component | Tight grammar bound |
| --- | ---: |
| field descriptor | 298 bytes |
| version-slot descriptor | 558 bytes |
| frame descriptor | 77,631 bytes |
| LeafExtension registry fragment | 27,417 bytes |
| impact row | 477,318 bytes plus source-path length, therefore at most 477,558 bytes |
| schema-impact manifest | 119,685 bytes |

These grammar bounds do not replace or relax the inclusive safety guards. The
immutable FrozenBase registry retains its exact frozen sequence and independent
oracle; the variable LeafExtension bound above follows its one-through-64 role
grammar. Every proportional count, edge aggregate, and canonical length uses
checked arithmetic. Production-valid boundary proof constructs the largest
semantically valid source-frozen value and every exact cardinality maximum.
Manifest admission independently recomputes the FrozenBase registry root and
every bound LeafExtension registry root from the exact kind, owner, fragment
ID, inherited FrozenBase root, and descriptor content before duplicate,
ordering, resolver, graph, or authority joins. A stored registry root is never
accepted as its own identity oracle.

When a jointly maximal combination cannot satisfy an earlier invariant, the
test proves the exact named refusal instead of padding, weakening, or mutating
the component to manufacture a guard-length value.

Guard reachability is proved separately at the shared canonical-frame seam. A
nonallocating count-only preflight accepts an exactly guard-length synthetic
frame and refuses guard plus one and checked arithmetic overflow before output
allocation. That seam proof is not evidence that a production component can
reach its guard. Private-field mutation, corrupted nested canonical bytes,
one-blob framing under a component magic, and reserved or padding fields are
not production-valid boundary fixtures.

For admitted values, canonical materialization performs the count-only pass
first, then makes exactly one output `Vec` allocation at the preflighted
capacity. The second pass must produce exactly the preflighted length without
growing beyond it. The component computes its root while borrowing that
buffer, then consumes the frame and transfers the private `Vec` by ownership.
There is no full-buffer clone, shrink-box conversion, second materialization
allocation, mutable or public byte exposure, or authority widening. The
complete composite-size preflight occurs before any
caller-cardinality-sized validation index; a divergent second pass cannot
produce a frame.

The schema-impact test log uses six distinct semantic partitions: positive,
expected refusal, expected failure, mutation, unsupported, and inapplicable.
Expected counts are result-free declarations; observed matched and mismatched
counts remain separate. A complete report must reconcile every expected case
exactly once, retain matched counts for each partition, reject missing, extra,
duplicate, reordered, skipped, or unlogged cases, and retain only bounded
typed diagnostics plus the first divergence. Logs and reports bind the
compatible-source snapshot root, schema-impact manifest root, close-repair
manifest root, and closed log-schema root.

The schema-owned production getter
`runner_v2_base_schema_impact_log_case_manifest_v1()` first reconstructs the
admitted production schema-impact manifest, then uses the crate-private
`source_frozen_schema_impact_log_case_manifest_v1` translator. The translator
exact-joins each declared case to one admitted manifest entry, the compiled
source-member root, the kind-checked FrozenBase or LeafExtension registry
identity, row owner and no-claim, row root, relation, derived local ordinal,
and exact predecessor and parent/child-slot counts. LeafExtension log context
retains its source-frozen owner and fragment ID; FrozenBase context admits
neither. Context, expected-case, and case-manifest construction remain
crate-private, so a public event cannot invent different source or registry
context.

This getter and translator declare expected partitions and expected-result
roots only. They execute no case, observe no runtime, produce no matched or
terminal count, inspect no retained artifact, perform no migration or hostile
decode, and mint no compatibility, scientific, admission, close-decision, or
other authority. Logging cannot create source-frozen authority merely by
repeating a root; its result-free case manifest is derived only from a
successfully validated schema-impact manifest.

The Phase-1 downstream handoff separately exposes the complete immutable
`RunnerV2PhaseOneExpectedLedgerV1`; the six closed partition names or their
aggregate counts are not a substitute for that ledger. Every
`RunnerV2PhaseOneExpectedLedgerCellV1` binds one contiguous ledger ordinal,
the exact source ordinal, distinct stable ledger-cell ID, source case ID,
source class and slash-separated path, close group and facet, execution scope,
partition, expected closed decision, exactly one expected-result root or
registered reason, semantic journey, route ID and owner, compatible-source
snapshot, source close-manifest and source-closure roots, optional historical
downstream-contribution root, Five Explicits root, source cell root, and
no-claim. The expected-result root is only a source-declared
expected-decision contract. It is never an observed or output root, and a
downstream executor must derive its actual root from the real outcome rather
than copying the declaration.

Ledger admission reconstructs every row from the exact full-set close
manifest, rejects any missing, extra, duplicate, reordered, reclassified,
wrong-result, wrong-reason, wrong-owner, wrong-route, stale, or mixed-snapshot
row, and independently reconciles all six per-partition counts and ordered
ledger-cell/source-case ID sets. Its root binds the independent six-name/code
vocabulary root, every complete ordered row, all ordered per-partition ID
sets, and zero required unexpected mismatches, execution failures, and
unexplained skips. A vocabulary-only, count-only, subset, or metadata-only
root cannot satisfy the ledger contract.

`RunnerV2PhaseOneContractContributionV2` binds that complete ledger root into
the downstream case-manifest declaration root, immutable payload root,
result-free contribution, distinct Deferred envelope, and aggregate
projection root. It still performs no downstream file read, process launch,
case execution, result matching, artifact retention, or authority grant. The
production-field mutation proof changes one real schema-field declaration
and verifies identity propagation through the enclosing frame, impact row,
impact manifest, schema-impact projection, payload, result-free contribution,
Deferred envelope, and aggregate while every unaffected sibling remains
byte-identical. Independent literal encoders, exact-set/order mutations, and
decision/root mismatch tests prevent any production table, advertised total,
or downstream result from acting as its own oracle.

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

The pure compiled source closure exact-set binds every file and generated input
owned or consumed by this leaf: implementation modules, this contract and
rustdoc, schema/domain/descriptor and constructor tables, independent oracles,
unit and property fixtures, mutation fixtures, in-process E2E projection rows,
log schemas, the current direct-dependency declaration, and build
configuration. Every entry carries one exact workspace-relative path, typed
owner, typed source route, declarative expected-source identity, common
snapshot policy and identity, encoded length, and content root. Missing, extra,
stale, reordered, ambient, duplicate-owner, wrong-route, wrong-identity,
wrong-length, or mixed-snapshot presented inputs move the closure or refuse.
The in-process closure corpus is exactly one positive reconstruction plus
fourteen deliberate refusal cases. Every refusal compares the complete
presented `ConstructionErrorV2` kind, field, expected text, and observed text;
an arbitrary error is never counted as a match.
The declarative identity and compile-time content root are not live revision or
supply-chain proof. Ambient live-tree discovery, lock/constellation pin
verification, generated-input discovery, and retained live evidence belong
solely to `frankensim-epic-foundations-huq.24.1.3.1`. A retained binary or root
cannot substitute for absent source.

## Retention and logging

The production base library prints nothing and performs no ambient logging.
Its deterministic test harness returns canonical conformance log data that
records:

- API, wire, source, build, toolchain, target, and feature roots;
- exact catalog and cap counts;
- state/reason/role/diagnostic matrix cells;
- typed expected and observed values;
- causal roots and manifest ordinals;
- diagnostic owner, prerequisite and repair counts, and no-claim scope; exact
  ordered diagnostic/repair values remain in the corresponding typed row or
  detail manifest instead of being flattened into open log fields;
- immutable manifest roots separately from caller-context execution roots;
- expected and observed detail-manifest roots, exact detail-cell counts, and a
  bounded typed first divergence;
- artifact-stored, system-publication-stored, and whole-publication-stored
  values as three distinct fields with exact stored-byte units and checked-sum
  reconciliation;
- minimized relative journey/case names and the exact typed downstream-script
  mapping;
- the closed logical-unit token and the exact count/stored-byte quantities
  admitted by the field catalog. Richer encoded, expanded, logical-work,
  cycle, duration, or performance inputs remain in their typed row/detail
  manifests rather than appearing as undeclared log fields.

Every logged source, build, toolchain, no-claim, cancelled, timed-out, and
controlled-internal-error digest validates both its closed role and its exact
registered nominal domain. Sharing `RunTerminal` or another broad digest role
does not permit cross-field substitution.

Logs do not retain absolute paths, PIDs, wall-clock timestamps, scheduler
latency, environment secrets, credentials, physical capability locators, raw
bulk payloads, or undeclared timing. Raw rejected payloads are represented by
bounded typed metadata or downstream retained relative artifact references.
The pure logging envelope admits at most 64 closed fields per event, 4,096
events per log, 1 MiB of canonical bytes per event, and 64 MiB for the complete
canonical log. These are inclusive size guards, not permission to violate the
closed event matrix or the base-E2E log's exact three-argument reproduction
tuple. The leaf-close writer separately uses the structured component-root
descriptor defined below. The largest legal publication terminal has exactly
63 fields and a valid 4,096-event state machine is admitted; 65 fields, 4,097
events, any structurally invalid value below a guard, and every checked
arithmetic overflow refuse before canonicalization, reconciliation, or any
further retained allocation.
This pure in-process leaf does not fabricate a retained-artifact claim by
pointing at one of its source files. Its report carries an explicit typed
absent-retained-artifact claim, and every log event emitted by
`run_base_e2e_projection_v1` has no retained-artifact path. The public event
constructor can separately represent a validated relative retained-artifact
path for downstream producers; that capability does not turn a source or
script mapping into retained evidence.

The leaf-close bounded writer is a separate complete-document state machine.
Before it admits any detail it reserves exactly the declared terminal-reserve
bytes for one canonical terminal event. Detail classes are closed and typed:
cell, stage, diagnostic, and first divergence. Their exact ordered manifest
binds every class and digest before admission. The writer preserves every
admitted prefix and seals permanently on the first detail that cannot fit
inside the remaining detail budget.

A normal stream ends in one complete `Complete` terminal. The first
over-budget detail instead emits one complete `LogBudgetExceeded` terminal
containing all of:

- rejected detail class and zero-based ordinal;
- rejected detail digest and exact omitted-suffix count;
- total byte budget and reserved terminal bytes;
- first-divergence stage;
- exact resource-return and drain outcomes;
- diagnostic owner and exact repair-manifest root;
- no-claim scope; and
- the structured minimized-reproduction descriptor.

The terminal itself must fit the bytes reserved before detail admission. A
missing, extra, duplicate, reordered, digest-substituted, or truncated detail;
post-terminal write; checked arithmetic overflow; absent terminal; incomplete
terminal; partial canonical document; overwritten prefix; silent truncation;
or overflow reported as success refuses. Canonical bytes and roots are
materialized only after `finish` has admitted the exact complete stream.
Normal and overflow documents replay byte-for-byte and root-for-root.
Mutating any terminal field moves both terminal and document identities.
Rejected caller-controlled data are represented only by stable class, ordinal,
digest, bounded counts, and registered redaction metadata; no error, Debug
surface, terminal, log, or reproduction string echoes the rejected value.

Coverage and execution reports use one exact, non-contradictory vocabulary:
eligible positive cells, matched positive cells, expected-refusal or mutation
cells, matched expected refusals, explicitly unsupported cells with typed
reasons, and unexpected mismatches. Green requires the two matched counts to
equal their eligible counts, exactly one result for every selected ID, and zero
unexpected mismatches. A nonzero expected-refusal corpus is required evidence;
it is not a request for a green run to contain a failure. Focused Phase-1 proof
stays distinct from repository policy, DSR, release-built E2E, and downstream
authority proof.

A structurally complete red log is still a valid, inspectable log document:
its summaries must exactly reconcile every failed partition. A row-red log has
at least one failed case terminal, and each such terminal carries the bounded
first-divergence identity/root. That root names either the exact typed
expected/observed detail divergence or, when the row contract itself fails
before a detail comparison exists, the exact typed row-contract divergence. A
source-closure-only red log may have no failed case terminal; its final
projection summary instead reports a nonzero exact source-closure failure
partition. `is_green` is the explicit conjunction of row and source-closure
success. Red evidence is never discarded merely because it does not satisfy
that predicate; malformed, incomplete, duplicated, or arithmetically
inconsistent logs still refuse construction. Logging validates the divergence
root's closed shape and placement; the projection producer must additionally
prove that its emitted value equals the exact retained divergence descriptor
root.

## Semantic seed and Five Explicits

`value::SeedMaterialV2` is the sole semantic workload-seed payload owned by
this leaf. It is exactly 32 bytes. Its canonical command-line grammar is
exactly `--seed none` or `--seed seed-256:` followed by 64 lowercase
hexadecimal digits. Omission, duplication, uppercase, alternate prefixes,
signs, separators, whitespace, Unicode, malformed lengths, and implicit or
ambient seed sources refuse. Seed data are reproducibility inputs only; they
mint no scientific, verification, admission, or publication authority.
Every `SeedErrorV2` variant has a unique stable code and name, the exact
`fs-evidence-runner/value` owner, typed retryability, one declarative
non-executable repair kind and target, and the same closed no-claim boundary.
Display, Debug, diagnostics, canonical refusal data, logs, and reproduction
text retain only bounded lengths, offsets, codes, and redacted classes; they
never echo rejected operands or seed material.

The three source-declared case policies remain distinct:

- `NoRandomness` carries no material and requires the registered
  no-randomness inapplicability record.
- `FixedManifest` obtains exact material only from the sealed case manifest
  and binds that manifest identity.
- `InvocationDerived` requires an explicitly provided, nonzero invocation
  seed and binds the exact stable case identity plus one registered derivation
  domain. Generator and minimizer versions are separate exact inputs.

No property, metamorphic, mutation, or fuzz label implies randomness. The
current base corpus is deterministic; any future randomized cell must be added
to the exact source oracle together with its complete provenance. An
applicable cell with absent source-backed material or provenance refuses
instead of deriving material from its case ID, facet, clock, process, worker,
thread, map order, or operating system.

Every stable close cell exact-joins one independent Five Explicits oracle row.
That row names separately rooted profiles for semantic numeric inputs,
semantic numeric grants, expected numeric observations, seed policy and
provenance, all seven hard/soft budget axes, declared version requirements,
the capability contract, and the no-claim scope. Shared component profiles are
permitted only when every cell names them explicitly, their resolved values
are bound into that cell root, and missing, extra, duplicate, reordered,
unknown, or unused rows or profiles refuse. Facet, source path, ordinal, and
execution scope are classification data and cannot synthesize semantic
explicit values.

The semantic numeric surface includes only named values crossing the case or
evidence contract boundary. Incidental implementation counters and source
ordinals are not semantic inputs. Numeric values use the complete closed
`NumericValueV2` union and bind either a canonical physical `UnitV2` or a
logical `LogicalUnitV2`; registered logical units additionally bind their
registry identity. A genuinely empty semantic surface is represented by an
explicit exact-empty declaration, never by an omitted profile.

Budgets retain seven axis codes and source-declared hard and soft values:
u64 nanoseconds for time, u64 logical bytes for memory, u128 plus an exact
logical unit for work, u32 for processes, and u64 encoded bytes for artifact,
output, and log limits. Soft values are never inferred from hard values.
Total-child and parallel-child constraints are separate process-shape fields;
neither is a soft process budget. Construction validates widths, units,
`soft <= hard`, governing Runner ceilings, checked aggregate feasibility, and
canonical-frame feasibility.

Version declarations preserve API generation, wire version and predecessor
policy, schema identity, source-closure identity, any source-known build and
toolchain identities, exact target or platform-matrix identity, exact
profile/configuration identity, and canonical feature-set identity. A hash of
the cell ID, class, path, facet, or scope is not a source, build, toolchain,
target, profile, or feature identity. Execution evidence supplies the actual
runtime-owned version fields and exact-joins them to the declaration.

Capability contracts use registered nominal capability IDs and bind the
capability registry/policy identity. A Bead owner, driver, or script token is
never a capability. Declarations retain exact required and permitted sets.
Execution evidence separately retains actual required, granted, observed,
returned, and revoked sets and validates membership, canonical ordering,
uniqueness, required/granted and observed/granted containment, and terminal
return/revocation reconciliation. A contribution-only row carries a typed
deferred-observation contract; it cannot represent unobserved capabilities as
an observed empty set.

Source-known explicits move the close-cell and close-manifest roots.
Runtime-resolved explicits move the presented-evidence, detailed-log,
minimized-reproduction, aggregate-execution, and report roots. Logs expose the
Five Explicits root and each component root directly so the first divergent
component is inspectable without reverse engineering an enclosing hash.
Reproduction data retain applicable identities and safe resolved values but
never raw forbidden values. Result-free declarations and immutable downstream
contributions remain distinct from execution proof.

`coverage::BaseCoverageManifestV1` is the sole source-authoritative,
result-free AC38 coverage inventory. The older projection-local coverage
inventory (`projection::BaseCoverageInventoryV1` and its local class/case
types) is compatibility-only auxiliary data and cannot replace, widen, or
prove the AC38 manifest. Likewise, compatibility accessors named only
`eligible`, `passed`, `failed`, or `root` cannot replace the exact partition
accessors or the explicitly named `manifest_root` and `execution_root`.
Each coverage `source_path` is a result-free mapping to the workspace-relative
source owner or designated external harness. A downstream harness path may be
declared before that downstream-owned script exists; the path is neither
retained evidence nor a claim that the harness ran. The in-process coverage
selection is derived only from the manifest's exact `ProjectionE2e`,
`RuntimeLogging`, and `SourceClosure` declarations. Runtime observation keys
are exact-set compared with that independent selection, so a missing, extra,
or non-local observation refuses instead of silently narrowing or widening
the checked scope.

## Runner V2 Stage-A base-values declaration and rootless local evaluator

`runner_v2::work_packages::base_values` owns the source-authoritative Stage-A
declaration and pure domain evaluator for foundational work package
`frankensim-epic-foundations-huq.24.1.1.1.1`, under the semantic ownership of
parent `frankensim-epic-foundations-huq.24.1.1.1`. The public declaration
surface is `declare_24_1_1_1_1_v1`,
`RunnerV2BaseValuesStageADeclarationV1`, and
`RunnerV2StageADeclarationRootV1`. The rootless handoff and its read-only
inspection types are public, but their constructors and
`evaluate_24_1_1_1_1_cell_v1` remain crate-private. Stage A does not expose a
public `run_*` function. Work package
`frankensim-epic-foundations-huq.24.1.1.1.7` is the sole owner of the future
public `run_24_1_1_1_1_v1` wrapper, fresh evaluator invocation, and
attempt-specific execution evidence.

The declaration is canonical result-free contract data. It exact-binds the
package ID, stable cells, independent expected-oracle rows, parent-projection
rows, inspectable fixture and operation declarations, per-cell result-free case
manifests, retained-domain obligations, auxiliary mutation obligations,
declaration-side Five Explicits, local route, deferred common requirements,
future broad-source requirements, the child owner-source fragment, the current
dependency-source closure, schema-impact deferral, rootless AC58 fragment,
shard and resume inapplicability, and no-claim. Its root contains no evaluator
invocation result or runtime actual.

### Exact operation, fixture, and oracle inventory

The new evaluator corpus contains exactly `71 * 12 + 15 = 867` source-ordered
cells. Every one of the 71 `RunnerLimitsV2` fields has a distinct cell for
every boundary kind in this exact order:

1. `Zero`
2. `One`
3. `StructuralMinimum`
4. `OneBelowStructuralMinimum`
5. `SmokeCeiling`
6. `SmokeTightened`
7. `SmokeOneOver`
8. `FullCeiling`
9. `FullTightened`
10. `FullOneOver`
11. `RepresentationalMaximum`
12. `CheckedRepresentationalOverflowRefusal`

A semantically undefined combination remains present with a typed
`Inapplicable` expectation. Coincident numeric values do not collapse cell
IDs. The canonical tightening rule, exact overflow operation, expected
outcome/reason/partition, and independent oracle are declaration data; none is
inferred from the evaluator result or from the production descriptor table
under test.

Each independent oracle row binds the cell ID, outcome, reason, partition,
ordered numeric cardinality, and every numeric name, heterogeneous value, and
unit. Diagnostic presence is explicit; when present, the root also binds its
code, owner, retryability, ordered prerequisite cardinality and identities,
and ordered repair cardinality, rank, kind, and target. Root-sensitivity tests
mutate each of those components independently, including numeric order,
diagnostic presence, prerequisite identity, and repair cardinality. The fresh
867-row exact join compares every nested field rather than only the outer
outcome. Its frozen aggregate is 422 accepted, 389 refused, and 56
inapplicable rows; 2,617 numeric observations; 445 diagnostics; 388 repairs;
and 56 prerequisites. The 852 limit rows account for 411 accepted, 388
refused, and 53 inapplicable rows, with refusal reasons partitioned as 201
above-ceiling, 20 fixed-representation, 59 below-minimum, 37 joint-feasibility,
and 71 checked-overflow cases.

The remaining fifteen cells are, in exact source order:

1. typed absence distinct from present zero;
2. named binary32 IEEE total order;
3. named binary64 IEEE total order;
4. exact `None` capability contract;
5. exact deferred-common-requirement membership;
6. reordered deferred-common-requirement refusal;
7. exact future-source membership;
8. rootless AC58 classification;
9. exact owner-source fragment;
10. exact local route;
11. diagnostic redaction and forbidden-value no-echo;
12. structured reproduction declaration;
13. compile-fail implicit IEEE ordering surface;
14. shard inapplicability; and
15. resume inapplicability.

The redaction cell consumes an owned sentinel value through the real redaction
boundary and then drops it. Its checked observed field is exactly the bounded
redaction placeholder, while the sentinel must be absent from the structured
observation, `Display`, `Debug`, and complete raw-cell `Debug` projections.
This is a non-echo proof over a presented sensitive value, not a test that
starts with an already-redacted placeholder.

Every cell publicly exposes its stable ordinal and ID, verification group,
typed operation, source-ordered companion normalization, independent oracle
root, and result-free case-manifest root. The shared limit fixture is likewise
inspectable: it has exactly one executable case, zero family rows for that
case, an explicitly present-empty declared-minimum list, and checked lifecycle
minimum `3 + (2 + 0) = 5`. A tightened limit cell carries an exact
source-ordered list with at most one row per companion field and the exact
same-width value needed to preserve joint feasibility. The per-cell case
manifest binds the exact operation, companion values, shared fixture, and
declaration-side Five Explicits; it is never an outcome or invocation record.

The declaration also preserves the complete ordered pre-Runner-V2 checked-value
obligation catalog outside these 867 new cells. Its eight closed facets are
numeric literals; units; tokens/text/paths; catalogs and nominal identities;
property/metamorphic behavior; mutation/refusal behavior; API/compile-fail
behavior; and fault/resource/no-mock integration. These retained obligations
preserve all existing integer, rational, decimal, IEEE, unit, token, text,
opaque-byte, path, catalog, nominal-root, normalization, determinism,
mutation, redaction, and local-integration functionality; they are not
silently counted as new evaluator cells. A separate exact 71-row auxiliary
mutation catalog binds one stable, field-ordinal-matched wrong-primitive-width
refusal obligation for every limit field. It is also outside the 867-cell
count.

### Direct Five Explicits and capability-none declaration

Stage A directly binds its source, build inputs, toolchain, target/profile,
schema inventory, and feature declaration. No ordinal, cell ID, class, path,
facet, scope, owner, expected partition, or other classification axis may
synthesize any of those identities. The exact dependency-source bytes form
the source identity. Build identity separately binds the workspace and crate
manifests, `Cargo.lock`, `constellation.lock`, this contract, and the
`fs-blake3` manifest. Toolchain identity binds `rust-toolchain.toml`. Target is
`TargetIndependentPureValidation`, profile is `CrateTest`, and the feature
declaration is an explicitly rooted exact-empty set validated against the
absence of a crate feature table and the absence of every implicit Cargo
feature introduced by `optional = true` in normal, build, nested dependency,
or target-specific dependency tables. Comments, quoted text, `optional =
false`, and unrelated metadata do not create false feature declarations.

Numeric inputs, numeric grants, and expected numeric observations are each
explicitly present and exactly empty. Seed is explicitly inapplicable under
`NoRandomnessByContract`. The exact `LocalSourceValidation` hard/soft budget
rows are:

| Axis | Hard | Soft | Unit |
| --- | ---: | ---: | --- |
| time | 60,000,000,000 | 45,000,000,000 | nanoseconds |
| memory | 536,870,912 | 402,653,184 | logical bytes |
| logical work | 1,000,000 | 750,000 | operations |
| processes | 1 | 0 | count |
| artifacts | 67,108,864 | 50,331,648 | encoded bytes |
| output | 5,242,880 | 4,194,304 | encoded bytes |
| logs | 67,108,864 | 50,331,648 | encoded bytes |

Every Stage-A cell binds the frozen
`BaseCoverageCloseCapabilityRegistryV1`,
`BaseCoverageCloseCapabilityProfileRegistryV1`, and
`BaseCoverageCloseCapabilityContractV1` for profile
`fs-evidence-runner.close-capability.none.v1`. Required and permitted
capability sets are both exactly empty. A Bead, owner, route, driver, script,
or physical resource is never a capability. Granted, observed, returned,
revoked, and resource-reconciliation values belong to later runtime evidence
and are structurally absent here.

### Deferred common contracts, routes, and projections

The exact ordered 31-row `RunnerV2CommonContractRequirementV1` catalog
preregisters every canonical, execution, and retention slot later owned by
work packages `.4`, `.5`, and `.6`. Each row binds its stable slot ID, API and
wire versions, no-predecessor rule, semantic and realization owners, future
nominal role/domain, exact nonempty subset of the canonical/execution/retention
planes, fulfillment stage, `.7` resolution owner, no-claim, and an uninhabited
typed future root that must be `Absent`. Missing, extra, duplicate, reordered,
unknown-version, wrong-owner, wrong-role/domain/plane/stage, wildcard,
range-only, predicted, classification-derived, copied, raw-hash, or
prematurely realized requirements refuse.

The child declares exactly one route:
`runner-v2.route.24-1-1-1-1.local.work-package.v1`, class `LocalOnly`,
future public entry point
`fs_evidence_runner::runner_v2::work_packages::run_24_1_1_1_1_v1`,
execution owner `.7`, capability profile `None`, and typed-absent external
driver. Route counts are exactly `LocalInProcess = 1`, `ExecutionOwned = 0`,
`ContributionOnly = 0`, and `Inapplicable = 0`. This is a route declaration,
not evidence that the wrapper exists or ran.

Every one of the 867 stable cells has exactly one ordered result-free
parent-projection row. It binds the cell, consumer route and owner, dispatcher,
expected partition, per-cell case-manifest root, no-claim, and the exact future
POSIX and native-Windows E2E paths:

- `scripts/ci/runner_v2_base_work_packages_e2e.sh`
- `scripts/ci/runner_v2_base_work_packages_e2e.ps1`

Wildcards, subset selection, transitive route inference, hidden skips, and
downstream status substitution refuse. Stage A does not claim either script
exists or executed.

### Source fragments and typed future source absence

The content-rooted owner-source fragment and the dependency-source closure are
distinct exact sets. The owner fragment has exactly these two child-owned
members:

- `crates/fs-evidence-runner/src/runner_v2/handoff.rs`
- `crates/fs-evidence-runner/src/runner_v2/work_packages/base_values.rs`

The current dependency-source closure has exactly sixteen content-rooted
members whose bytes can change Stage-A semantics:

- `crates/fs-blake3/src/lib.rs`
- `crates/fs-evidence-runner/src/lib.rs`
- `crates/fs-evidence-runner/src/canonical.rs`
- `crates/fs-evidence-runner/src/catalog.rs`
- `crates/fs-evidence-runner/src/construction.rs`
- `crates/fs-evidence-runner/src/coverage.rs`
- `crates/fs-evidence-runner/src/identity.rs`
- `crates/fs-evidence-runner/src/limits.rs`
- `crates/fs-evidence-runner/src/path.rs`
- `crates/fs-evidence-runner/src/projection.rs`
- `crates/fs-evidence-runner/src/schema_impact.rs`
- `crates/fs-evidence-runner/src/value.rs`
- `crates/fs-evidence-runner/src/runner_v2.rs`
- `crates/fs-evidence-runner/src/runner_v2/handoff.rs`
- `crates/fs-evidence-runner/src/runner_v2/work_packages.rs`
- `crates/fs-evidence-runner/src/runner_v2/work_packages/base_values.rs`

The two-member ownership statement cannot stand in for the sixteen-member
semantic dependency closure, and neither can stand in for the complete broad
crate closure or live revision/supply-chain proof.

Every current owner and dependency member is declared from one
workspace-relative source literal that feeds both its logical path and its
compile-time byte inclusion. Exact-set validation independently recomputes the
expected member root from separately included expected bytes and the nominal
owner/dependency hash domain; it does not merely trust the presented root.
Negative conformance tests refuse missing, extra, duplicate, reordered, and
wrong-path members as well as wrong, swapped, correctly resealed altered, and
wrong-domain content roots.

The existing broad source inventory has 27 members. Stage A freezes these
thirteen final ordinals and paths as typed-absent future content requirements,
yielding an eventual 40-member broad inventory:

1. `crates/fs-evidence-runner/src/runner_v2.rs`
2. `crates/fs-evidence-runner/src/runner_v2/handoff.rs`
3. `crates/fs-evidence-runner/src/runner_v2/work_packages.rs`
4. `crates/fs-evidence-runner/src/runner_v2/work_packages/base_values.rs`
5. `crates/fs-evidence-runner/src/runner_v2/work_packages/diagnostics.rs`
6. `crates/fs-evidence-runner/src/runner_v2/work_packages/schema_registry.rs`
7. `crates/fs-evidence-runner/src/runner_v2/work_packages/runtime_evidence.rs`
8. `crates/fs-evidence-runner/src/runner_v2/work_packages/routes.rs`
9. `crates/fs-evidence-runner/src/runner_v2/work_packages/detailed_logging.rs`
10. `crates/fs-evidence-runner/src/runner_v2/work_packages/execution.rs`
11. `crates/fs-evidence-runner/tests/runner_v2_base_work_packages.rs`
12. `scripts/ci/runner_v2_base_work_packages_e2e.sh`
13. `scripts/ci/runner_v2_base_work_packages_e2e.ps1`

Every future content root is `TypedOptionV1::Absent` in Stage A, including for
a path whose local source happens to exist. The final ordinal/path inventory
is a requirement for `.7` to resolve on one compatible snapshot, not a present
40-member source-closure claim.

### AC58 ownership and the rootless handoff

The owned Stage-A schema inventory is exactly 43 distinct names. Exactly 42
canonical names are exact-set in `RunnerV2SchemaImpactDeferralV1` for resolution
by the dedicated `.3` schema-registry owner, including the
`runner-v2-raw-outcome-reason-contract-v1` compatibility contract. The
forty-third name is the separately classified rootless handoff described below.
The future `.3` manifest root is typed absent. Stage A does not fabricate a
`CanonicalSchemaImpactRowV1`, historical/current frame, legacy container,
migration result, or schema-impact manifest root for those schemas.

The sole local AC58 exception is the exact
`runner-v2-local-work-package-handoff-v1` semantic type. Its lightweight
`RunnerV2RootlessAc58FragmentV1` declares
`InapplicableNoCanonicalFrame`, `NoSchemaPredecessor`, an explicitly
present-empty authority surface, and the Stage-A no-claim. This fragment is
bound by the Stage-A declaration but is not a canonical schema-impact row and
does not give the handoff canonical bytes or an identity.

`RunnerV2LocalWorkPackageHandoffV1` is a bounded, source-ordered,
noncanonical, non-authoritative Rust value. The admitted complete report has
at most 2,048 cells; each cell has at most eight safe typed numeric
observations; each diagnostic has at most eight strictly ordered
prerequisites and four contiguously ranked non-executable repairs. Accepted
cells carry no diagnostic, and every refused, failed, unsupported, or
inapplicable cell carries exactly one. Cell IDs and observation names are
strictly ordered and unique.

The handoff contains only package and cell identity, raw checked outcome and
reason, safe typed numeric observations, and bounded structured diagnostics.
It contains no expected result, oracle root, executable repair text, canonical
bytes or root, attempt identity, actual Five Explicits, AC57 disposition or
envelope, canonical route/log/reproduction/receipt/telemetry root, execution
partition, retained-artifact claim, physical capability/resource value, or
authority. Public accessors permit inspection only; compile-fail tests freeze
the absence of caller construction, root conversion, attempt access, and
canonical-content substitution.

### Stage-A proof and no-claim boundary

`declare_24_1_1_1_1_v1` never invokes the evaluator. Unit, literal, boundary,
property/metamorphic, state/model, mutation, API/compile-fail,
fault/resource/cancellation-model, safe-field/redaction, reproduction,
source-closure, and no-mock local-integration tests invoke the real
crate-private evaluator afresh and compare its complete source-ordered handoff
with independent declaration-side expectations. Those reports and verdicts
are ephemeral verification; they do not become `.7` execution evidence.

The common JSONL event, terminal reservation, redaction, first-divergence,
reproduction, relative-artifact, raw-audit, telemetry, and operator-view
contracts are typed future requirements for `.6`. Stage A emits no JSONL,
terminal log, reproduction instance, retained artifact, receipt, telemetry, or
operator view. Sharding and resume are explicitly inapplicable because the
bounded local evaluator is a complete single pass and each invocation
recomputes the entire package.

Acceptance proves only that the canonical Stage-A declaration is complete,
that the pure evaluator obeys its checked domain contract when freshly
invoked, and that the rootless handoff preserves its structural no-authority
boundary. It proves no external or release-built execution, actual runtime
Five Explicits, AC57 observation, acquired capability, resource
reconciliation, retained or durable artifact, common log/reproduction/receipt
instance, telemetry, scientific validity, admission, promotion, DSR success,
release E2E result, or downstream authority.

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
9. exact source-closure and coverage-manifest tests, with unit, boundary,
   property/metamorphic, schema/descriptor, mutation, compile-fail/doctest,
   no-mock integration, projection, logging, and source-closure classes
   reported separately;
10. in-process, source-closed Runner V2 base E2E projection tests.

The frozen base coverage inventory contains exactly 78 source-classified
compile-fail contracts. Ten of them independently reject implicit `cmp`, `<`,
`BTreeSet`, slice `sort`, and `sort_by_key` use for each of the binary32 and
binary64 wrapper types. Each fence carries one expected compiler error code,
and the inventory test exact-matches the complete error-code distribution. The
frozen base total is therefore exactly `217 + 78 + 29 = 324`; adding the exact
121 close-extension cells yields exactly 445 full source-manifest and close
cells.

Command-applicability evidence enumerates all `4^5` selector-cardinality
vectors for each of the six commands, accepts exactly one vector per command,
and compares the complete Usage refusal. Profile-bound evidence separately
checks both `Smoke` and `Full` exact ceilings and one-over refusals; passing one
profile cannot stand in for the other.

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
- `scripts/ci/e2e_evidence_verifier_v2.sh`
- `scripts/ci/canonical_evidence_runner_v2.sh`
- `scripts/ci/verify_runner_rjoq_handoff_v1.sh`

The exact routes are:

| Journey | Downstream owner | Driver | Script or native route | Immutable case manifest |
| --- | --- | --- | --- | --- |
| publication-state-v2 | `frankensim-epic-foundations-huq.24.2.2.2` | `e2e-evidence-runner-publication-state-v2-driver` | `scripts/ci/e2e_evidence_runner_publication_state_v2.sh` | `scripts/ci/manifests/evidence_runner_publication_state_v2_cases.v1.json` |
| publication-v2 | `frankensim-epic-foundations-huq.24.2.2.3.3` | `e2e-evidence-runner-publication-v2-driver` | `scripts/ci/e2e_evidence_runner_publication_v2.sh` | `scripts/ci/manifests/evidence_runner_publication_v2_cases.v1.json` |
| verifier-v2 | `frankensim-epic-foundations-huq.24.3.3.3.3` | `e2e-evidence-verifier-v2-driver` | `scripts/ci/e2e_evidence_verifier_v2.sh` | `scripts/ci/manifests/evidence_verifier_v2_cases.v1.json` |
| canonical-runner-v2 | `frankensim-epic-foundations-huq.24.4.1.4` | `canonical-evidence-runner-v2-e2e-driver` | `scripts/ci/canonical_evidence_runner_v2.sh` | `scripts/ci/manifests/canonical_evidence_runner_v2_cases.v1.json` |
| rjoq-handoff-v1 | `frankensim-epic-foundations-huq.24.5.3.1` | `verify-runner-rjoq-handoff-v1-driver` | `scripts/ci/verify_runner_rjoq_handoff_v1.sh` | `scripts/ci/manifests/runner_rjoq_handoff_verifier_v1_cases.v1.json` |

The five row lists are handwritten and journey-specific. Every row binds its
journey, exact downstream owner, driver, script or native route, immutable
manifest path and manifest root, consumption rationale, bounded fixture or
closed subcase-manifest reference, exact expected decision and detail,
semantic-cell partition, unit, no-claim scope, compiled source-closure root,
and closed log-schema root. The immutable context-free journey-manifest root
is distinct from the caller-context-bound journey-execution root.

Each row exposes an exact ordered, bounded, public detail-cell manifest.
Callers can inspect every stable cell ID, semantic ordinal, expected decision,
and closed typed expected/refusal/unsupported payload without reverse
engineering an opaque hash. Detail manifests reject missing, extra, duplicate,
reordered, stale, or substituted cells. A containing row may bind the existing
bounded `RegisteredDecisionDetailProjectionV2`; that reference includes its
registered namespace, local code, content root, encoded length, registry root,
and projection root, but neither extends `ActionableDiagnosticV2` nor mints
comparison, verification, scientific, or admission authority.

The public handoff includes the typed row, detail-cell, detail-manifest,
presented-result, checked-result, comparison-report, execution-report,
aggregate-report, harness-context, and explicit retained-artifact-absence
vocabularies. Expected and observed detail-manifest roots remain distinct, as
do immutable `manifest_root`, comparison-only `comparison_root`, private
row-witness roots, and caller-context-bound `execution_root`; a generic `root`
compatibility accessor does not erase those distinctions.
External observers construct cells through the bounded checked
`BaseE2eDetailCellV1::new` path and close complete ordered slices through
`BaseE2eDecisionDetailManifestV1::from_cells`. The row-based
`BaseE2ePresentedRowResultV1::new_with_observed_detail_cells` constructor
recomputes and checks kind, ordinal range, order, unique IDs and roots, manifest
count and root, matched-cell count, and the exact first divergent stable ID.
That public constructor is comparison-only. Exact caller-presented cells and
aggregate counts, including a positive-only row with an empty detail manifest,
may produce a green equality result only inside
`BaseE2eJourneyComparisonReportV1`. The report has its own comparison domain
and `comparison_root`; it has no execution witness, no `execution_root` API,
and no conversion into `BaseE2eJourneyExecutionReportV1`. The deprecated
`join_base_e2e_journey_results_v1` name is only a source-migration alias for
the comparison API; its retained harness argument cannot change that type or
mint execution evidence.
When the typed cells are exact but another semantic row partition is red, the
only admitted failure ID is the closed `row.contract` sentinel, and the
row-contract divergence root binds that exact value. The constructor retains
only the first typed observed divergence after validation. The older
descriptor-only constructors are deprecated compatibility paths because an
opaque root/count cannot establish typed cell membership or order. Their
root-bound observation mode remains explicitly unverified: they may be
compared and inspected as mismatches, but they can never satisfy a green
typed-detail comparison or expose cached cells as caller-observed cells.

Only `run_base_e2e_journey_v1` and the five-journey aggregate can construct
`BaseE2eJourneyExecutionReportV1`. Their private executed-row finalizer consumes
`BaseE2eCaseExecutionV1` directly, reconstructs and reconciles every decision,
semantic partition, detail manifest, matched-detail count, and first
divergence from the retained in-process cells, and never binds caller-shaped
`BaseE2ePresentedRowResultV1` data. It then mints one private witness per row
under the distinct source-closed in-process execution class and domain. The
witness binds the journey manifest root and code, one-based ordered row
ordinal, row ID and kind, semantic manifest root, journey mapping root,
harness-context root, all actual observed partitions, detail
root/count/matched values, first divergence, and the witness-absent comparison
row root. The checked execution-row root retains the exact witness root, and
the journey-execution root retains every ordered row-result and witness root.
Journey finalization rejects missing, extra, duplicate, reordered,
reconstructed-mismatch, or substituted witnesses. Comparison rows retain
explicit witness absence and cannot enter this finalizer.
`Accept`, `Refuse`, and `Unsupported` cells have disjoint payload rules;
limit/budget repair ranks are bounded to 1 through 16 and their owner/target
labels must be canonical `StableTokenV2` values.
The log validator rejects duplicate journey manifest roots, duplicate journey
execution roots, any manifest/execution cross-substitution, and either
aggregate root reusing a journey execution root. The projection layer additionally
proves exact equality between every logged scoped execution root and its
journey report, and between the logged aggregate execution root and the
aggregate report.

The rows cover every base literal, limit, budget, unit, capability-policy,
one-field identity mutation, exact diagnostic, and no-claim boundary consumed
by each journey. This leaf runs every eligible projection row **in process**
through real constructors and validators with no mocks. Observations
are not cached across rows or journeys: each reported row count therefore
binds a distinct execution, while repeated semantic cases in different
journeys are deliberately rerun. Refusal and
unsupported observations compare exact typed error kind or code, field, owner,
expected and observed values, unit, repair data, and unsupported reason.
Aggregated matrices bind every ordered stable subcase ID and detail into a
closed manifest with an exact count. Caller-presented row results are
exact-compared by journey, row ID, semantic root, detail root, and order but
never enter execution aggregation. In-process executions instead pass through
the private witness finalizer before aggregation.
Catalog-size log fields are frozen expected-inventory metadata, not claims
about how far one observed execution progressed. `checked_cells` is the exact
nonzero observed execution prefix: it equals the inventory only for a complete
matrix run and remains shorter when a red run stops early. Failures while
constructing command applicability/setup, the base budget, or the capability
registry are explicitly assigned semantic ordinal 1. Matrix and aggregation
helpers reject zero progress instead of padding it to one, so no unobserved
cell can be fabricated merely to make a red result structurally nonzero.
Expected catalog, limit, budget, diagnostic, and command rows come from
handwritten source-closed oracle tables independent of the production lookup
tables under test, and every oracle-table root is bound into its semantic
manifest. Literal-table encoders accept their explicit handwritten row slices,
bind independently matched nominal identities without calling production
catalog accessors, and are mutation-tested through both the table root and the
containing semantic row root. Mixed positive/refusal matrices retain typed
per-subcase progress, so an early refusal-cell failure cannot be mislabeled as
a positive-cell failure. The aggregate accepts only the execution-report type,
retains the five ordered scoped execution roots, and binds them into one
caller-context aggregate execution root; a comparison report is structurally
ineligible. It does not discard execution roots after computing counts or
substitute the context-free projection manifest root.
Deterministic detailed logs reconcile the positive, expected-refusal,
unsupported, and unexpected-mismatch partitions plus separately nonzero
projection-E2E and logging counts. Publication-storage log events expose the
three independently projected accounting cells—artifact-stored,
system-publication-stored, and whole-publication-stored—without claiming that
any bytes were physically retained; the closed log schema requires exact
stored-byte typing and checked `artifact + system = publication`
reconciliation.

This leaf does not create, edit, invoke, or claim release-built execution of
those scripts. Each sole downstream script owner must exact-set consume its
corresponding immutable journey-manifest root, separately record and verify its
caller-context journey-execution root, reject missing, extra, duplicate,
reordered, stale, unmapped, cross-journey, context-substituted, or
manifest/execution-root-substituted inputs, and execute every eligible row.
Final program closure requires all five release-built projections; their
execution is not a prerequisite for closing this pure base-schema leaf.
None of these public descriptors, roots, reports, or in-process results claims
live source provenance, retained runtime evidence, downstream script
execution, durable publication, scientific verification, or admission.

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
- Root-free ownership registries, generated grammar, cross-leaf and external
  coverage joins, and exact registry inventory:
  `frankensim-epic-foundations-huq.24.1.3.2`. This leaf retains ownership of
  its own in-process manifest/subcase/result/log joins.

The sole downstream owners of the five release-built scripts are respectively
`frankensim-epic-foundations-huq.24.2.2.2`,
`frankensim-epic-foundations-huq.24.2.2.3.3`,
`frankensim-epic-foundations-huq.24.3.3.3.3`,
`frankensim-epic-foundations-huq.24.4.1.4`, and
`frankensim-epic-foundations-huq.24.5.3.1`. This crate provides immutable
projection inputs to those owners; it does not duplicate their script
ownership.

## No-claim boundaries

Acceptance by this crate means only that base Runner V2 data are
well-formed, canonical, within their declared limits and budgets, and
consistent with the closed pure-validation matrices.

It does **not** claim:

- that any nominal source, build, toolchain, case, artifact, object,
  externally retained log or manifest, publication, receipt, acquired policy,
  or resource exists outside the bounded in-memory data constructed here;
- that presented bytes have the asserted meaning or match retained content;
- that any downstream or release-built case ran; the pure in-process
  projection-row validation performed here does not imply that a lifecycle
  completed, cancellation drained, or a process terminated;
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

Presented roots, compiled source identities, logical publication choices,
policy roots, diagnostics, coverage manifests, and in-process E2E projection
results remain non-authoritative data. Scientific, live-source, durability,
verification, admission, promotion, release, and rjoq authority can arise only
from their separately owned, explicitly capability-bearing protocols and
retained evidence.
