//! Closed Runner V2 wire catalogs.
//!
//! This module owns names and unsigned 16-bit tags only. It performs no byte
//! parsing, registration lookup, lifecycle work, or authority-bearing
//! admission. Unknown tags are refused rather than preserved or guessed.

use core::fmt;
use core::num::NonZeroU16;

/// The public Runner product generation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct RunnerApiGeneration(u16);

impl RunnerApiGeneration {
    /// The sole admitted API generation for this crate.
    pub const V2: Self = Self(2);

    /// Every admitted API generation, in increasing numeric order.
    pub const ALL: [Self; 1] = [Self::V2];

    /// Decode the exact API-generation number.
    ///
    /// # Errors
    ///
    /// Returns [`UnknownCatalogCode`] for every value other than two.
    pub const fn from_code(code: u16) -> Result<Self, UnknownCatalogCode> {
        match code {
            2 => Ok(Self::V2),
            _ => Err(UnknownCatalogCode::new("RunnerApiGeneration", code)),
        }
    }

    /// Exact numeric generation.
    #[must_use]
    pub const fn code(self) -> u16 {
        self.0
    }

    /// Frozen public product-generation name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        let _ = self;
        "RunnerSpecV2"
    }
}

/// The frozen Runner wire-schema version.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct RunnerWireVersion(u16);

impl RunnerWireVersion {
    /// The first and sole admitted wire version.
    pub const V1: Self = Self(1);

    /// Every admitted wire version, in increasing numeric order.
    pub const ALL: [Self; 1] = [Self::V1];

    /// Decode the exact wire version.
    ///
    /// # Errors
    ///
    /// Returns [`UnknownCatalogCode`] for every value other than one.
    pub const fn from_code(code: u16) -> Result<Self, UnknownCatalogCode> {
        match code {
            1 => Ok(Self::V1),
            _ => Err(UnknownCatalogCode::new("RunnerWireVersion", code)),
        }
    }

    /// Exact numeric wire version.
    #[must_use]
    pub const fn code(self) -> u16 {
        self.0
    }

    /// Stable wire-version name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        let _ = self;
        "runner-wire-v1"
    }

    /// Frozen predecessor rule for this wire version.
    #[must_use]
    pub const fn predecessor_policy(self) -> WirePredecessorPolicyV1 {
        let _ = self;
        WirePredecessorPolicyV1::NoPredecessor
    }
}

/// Wire V1's closed predecessor policy.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum WirePredecessorPolicyV1 {
    /// Wire V1 has no predecessor, migration alias, or legacy decoder.
    NoPredecessor,
}

impl WirePredecessorPolicyV1 {
    /// Every admitted predecessor policy.
    pub const ALL: [Self; 1] = [Self::NoPredecessor];

    /// Stable policy name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::NoPredecessor => "no-predecessor",
        }
    }

    /// Predecessor version, absent by construction for wire V1.
    #[must_use]
    pub const fn predecessor(self) -> Option<RunnerWireVersion> {
        match self {
            Self::NoPredecessor => None,
        }
    }
}

/// Frozen API-generation value for `RunnerSpecV2`.
pub const RUNNER_SPEC_V2_API_GENERATION: RunnerApiGeneration = RunnerApiGeneration::V2;

/// Frozen first Runner wire-schema version.
pub const RUNNER_V2_WIRE_VERSION: RunnerWireVersion = RunnerWireVersion::V1;

/// Frozen no-predecessor rule for Runner wire V1.
pub const RUNNER_V2_PREDECESSOR_POLICY: WirePredecessorPolicyV1 =
    WirePredecessorPolicyV1::NoPredecessor;

/// A closed catalog rejected an unknown unsigned 16-bit code.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct UnknownCatalogCode {
    catalog: &'static str,
    code: u16,
}

impl UnknownCatalogCode {
    const fn new(catalog: &'static str, code: u16) -> Self {
        Self { catalog, code }
    }

    /// Rust type name of the closed catalog that refused the code.
    #[must_use]
    pub const fn catalog(self) -> &'static str {
        self.catalog
    }

    /// Unknown code that was refused.
    #[must_use]
    pub const fn code(self) -> u16 {
        self.code
    }
}

impl fmt::Display for UnknownCatalogCode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} refuses unknown u16 code {}", self.catalog, self.code)
    }
}

impl std::error::Error for UnknownCatalogCode {}

macro_rules! count_variants {
    ($($variant:ident),+ $(,)?) => {
        <[()]>::len(&[$(count_variants!(@one $variant)),+])
    };
    (@one $variant:ident) => { () };
}

macro_rules! closed_u16_catalog {
    (
        $(#[$enum_meta:meta])*
        pub enum $catalog:ident {
            $($variant:ident = $code:literal => $stable_name:literal),+ $(,)?
        }
    ) => {
        $(#[$enum_meta])*
        #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
        #[repr(u16)]
        pub enum $catalog {
            $(
                #[doc = concat!("Stable catalog value `", $stable_name, "`.")]
                $variant = $code,
            )+
        }

        impl $catalog {
            /// Every catalog value, in frozen wire order.
            pub const ALL: [Self; count_variants!($($variant),+)] = [
                $(Self::$variant,)+
            ];

            /// Decode one exact unsigned 16-bit catalog code.
            ///
            /// # Errors
            ///
            /// Returns [`UnknownCatalogCode`] when `code` is not one of the
            /// frozen discriminants.
            pub const fn from_code(code: u16) -> Result<Self, UnknownCatalogCode> {
                match code {
                    $($code => Ok(Self::$variant),)+
                    _ => Err(UnknownCatalogCode::new(stringify!($catalog), code)),
                }
            }

            /// Exact unsigned 16-bit wire code.
            #[must_use]
            pub const fn code(self) -> u16 {
                self as u16
            }

            /// Stable lowercase catalog name.
            #[must_use]
            pub const fn name(self) -> &'static str {
                match self {
                    $(Self::$variant => $stable_name,)+
                }
            }
        }
    };
}

closed_u16_catalog! {
    /// Closed Runner V2 terminal-state catalog.
    pub enum ProofExitV2 {
        Pass = 0 => "pass",
        Failed = 10 => "failed",
        Refused = 11 => "refused",
        NoData = 12 => "no-data",
        Stale = 13 => "stale",
        EnvironmentInvalid = 14 => "environment-invalid",
        Blocked = 15 => "blocked",
        Unsupported = 16 => "unsupported",
        NotRun = 17 => "not-run",
        Cancelled = 18 => "cancelled",
        TimedOut = 19 => "timed-out",
        Usage = 64 => "usage",
        InternalError = 70 => "internal-error",
    }
}

closed_u16_catalog! {
    /// Closed reasons carried only by [`ProofExitV2::Refused`].
    pub enum RefusedReasonV2 {
        InvalidEvidence = 1 => "invalid-evidence",
        NonCanonicalEvidence = 2 => "non-canonical-evidence",
        EvidenceIdentityMismatch = 3 => "evidence-identity-mismatch",
        EvidenceTampered = 4 => "evidence-tampered",
        LimitExceeded = 5 => "limit-exceeded",
        UnsafeArtifactPlacement = 6 => "unsafe-artifact-placement",
        ArtifactCollision = 7 => "artifact-collision",
        LifecycleViolation = 8 => "lifecycle-violation",
        PolicyRefused = 9 => "policy-refused",
        AuthorityBoundaryViolation = 10 => "authority-boundary-violation",
        MigrationRefused = 11 => "migration-refused",
    }
}

closed_u16_catalog! {
    /// Closed Runner V2 command catalog.
    pub enum RunnerCommandV2 {
        List = 0 => "list",
        Check = 1 => "check",
        SelfTest = 2 => "self-test",
        Run = 3 => "run",
        Negative = 4 => "negative",
        Replay = 5 => "replay",
    }
}

closed_u16_catalog! {
    /// Closed Runner execution-profile catalog.
    pub enum RunProfileV2 {
        Smoke = 1 => "smoke",
        Full = 2 => "full",
    }
}

closed_u16_catalog! {
    /// Closed artifact-disposition catalog.
    pub enum ArtifactDispositionV2 {
        LifecycleOnlyNoBundle = 1 => "lifecycle-only-no-bundle",
        DurableBundleRequired = 2 => "durable-bundle-required",
    }
}

closed_u16_catalog! {
    /// Closed logical platform-path profiles.
    pub enum PlatformPathProfileV2 {
        PosixDescriptorRelativeV1 = 1 => "posix-descriptor-relative-v1",
        WindowsHandleRelativeV1 = 2 => "windows-handle-relative-v1",
        ContentStoreObjectKeyV1 = 3 => "content-store-object-key-v1",
    }
}

closed_u16_catalog! {
    /// Closed lifecycle-record vocabulary.
    pub enum LifecycleRecordKindV2 {
        RunStart = 1 => "run-start",
        CaseStart = 2 => "case-start",
        FamilyRow = 3 => "family-row",
        CaseTerminal = 4 => "case-terminal",
        RunSummary = 5 => "run-summary",
        RunTerminal = 6 => "run-terminal",
    }
}

closed_u16_catalog! {
    /// Closed roles whose records carry a terminal state.
    pub enum StateBearingRecordRoleV2 {
        PreRunDiagnostic = 1 => "pre-run-diagnostic",
        ExecutedCaseTerminal = 2 => "executed-case-terminal",
        SuppressedCaseTerminal = 3 => "suppressed-case-terminal",
        RunTerminal = 4 => "run-terminal",
    }
}

closed_u16_catalog! {
    /// Closed base diagnostic-code catalog.
    pub enum DiagnosticCodeV2 {
        CaseConformanceMismatch = 1 => "case.conformance_mismatch",
        RunnerNotRun = 2 => "runner.not_run",
        RunnerRefused = 3 => "runner.refused",
        RunnerNoData = 4 => "runner.no_data",
        RunnerStale = 5 => "runner.stale",
        RunnerEnvironmentInvalid = 6 => "runner.environment_invalid",
        RunnerBlocked = 7 => "runner.blocked",
        RunnerUnsupported = 8 => "runner.unsupported",
        RunnerCancelled = 9 => "runner.cancelled",
        RunnerTimedOut = 10 => "runner.timed_out",
        RunnerUsage = 11 => "runner.usage",
        RunnerInternalError = 12 => "runner.internal_error",
    }
}

closed_u16_catalog! {
    /// Closed diagnostic retryability catalog.
    pub enum RetryabilityV2 {
        Never = 0 => "never",
        SameInvocation = 1 => "same-invocation",
        AfterInputChange = 2 => "after-input-change",
        AfterEnvironmentChange = 3 => "after-environment-change",
        AfterPrerequisiteChange = 4 => "after-prerequisite-change",
    }
}

closed_u16_catalog! {
    /// Closed structured repair-action kinds.
    pub enum RepairActionKindV2 {
        ChangeArguments = 1 => "change-arguments",
        SupplyEvidence = 2 => "supply-evidence",
        RegenerateCanonicalEvidence = 3 => "regenerate-canonical-evidence",
        RefreshEvidence = 4 => "refresh-evidence",
        ReduceResourceDemand = 5 => "reduce-resource-demand",
        ChooseSafeArtifactDestination = 6 => "choose-safe-artifact-destination",
        RestoreLifecycle = 7 => "restore-lifecycle",
        UpdatePolicyOrCapability = 8 => "update-policy-or-capability",
        RegisterMigration = 9 => "register-migration",
        RetrySameInvocation = 10 => "retry-same-invocation",
        ContactOwner = 11 => "contact-owner",
        InspectRetainedArtifact = 12 => "inspect-retained-artifact",
    }
}

closed_u16_catalog! {
    /// Closed causal codes for `NotRun` slots.
    pub enum NotRunCauseCodeV2 {
        PriorCancelled = 1 => "prior-cancelled",
        PriorTimedOut = 2 => "prior-timed-out",
        PriorControlledInternalError = 3 => "prior-controlled-internal-error",
    }
}

closed_u16_catalog! {
    /// Closed outer tags for Runner V2 typed values.
    pub enum TypedValueTagV2 {
        I8 = 1 => "i8",
        I16 = 2 => "i16",
        I32 = 3 => "i32",
        I64 = 4 => "i64",
        I128 = 5 => "i128",
        U8 = 6 => "u8",
        U16 = 7 => "u16",
        U32 = 8 => "u32",
        U64 = 9 => "u64",
        U128 = 10 => "u128",
        Rational = 11 => "rational",
        Decimal = 12 => "decimal",
        F32Bits = 13 => "f32-bits",
        F64Bits = 14 => "f64-bits",
        Digest = 15 => "digest",
        Quantity = 16 => "quantity",
        Token = 17 => "token",
        Text = 18 => "text",
        RelativePath = 19 => "relative-path",
        OpaqueBytes = 20 => "opaque-bytes",
    }
}

closed_u16_catalog! {
    /// Closed option-sum tags; absence never borrows a payload sentinel.
    pub enum TypedOptionTagV1 {
        Absent = 0 => "absent",
        Present = 1 => "present",
    }
}

closed_u16_catalog! {
    /// Closed digest-role tags.
    pub enum DigestRoleV2 {
        Spec = 1 => "spec",
        Invocation = 2 => "invocation",
        Run = 3 => "run",
        Source = 4 => "source",
        Build = 5 => "build",
        Toolchain = 6 => "toolchain",
        CaseManifest = 7 => "case-manifest",
        ArtifactEncoded = 8 => "artifact-encoded",
        ArtifactContent = 9 => "artifact-content",
        StoredObject = 10 => "stored-object",
        ArtifactInventory = 11 => "artifact-inventory",
        LifecycleLog = 12 => "lifecycle-log",
        RunSummary = 13 => "run-summary",
        RunTerminal = 14 => "run-terminal",
        BundleManifest = 15 => "bundle-manifest",
        DurablePublication = 16 => "durable-publication",
        Seal = 17 => "seal",
        PublishedBundleReceipt = 18 => "published-bundle-receipt",
        Policy = 19 => "policy",
        CandidateBytes = 20 => "candidate-bytes",
        CandidateSchema = 21 => "candidate-schema",
        SourceClosure = 22 => "source-closure",
        ClaimScope = 23 => "claim-scope",
        ProducerManifest = 24 => "producer-manifest",
        RegisteredFamilyDomain = 25 => "registered-family-domain",
    }
}

closed_u16_catalog! {
    /// Closed logical publication protocols.
    pub enum PublicationProtocolV2 {
        PosixDescriptorRenameAndDirectorySyncV1 = 1
            => "posix-descriptor-rename-and-directory-sync-v1",
        WindowsHandleReplaceAndDirectoryFlushV1 = 2
            => "windows-handle-replace-and-directory-flush-v1",
        ContentStoreAtomicCommitV1 = 3 => "content-store-atomic-commit-v1",
    }
}

closed_u16_catalog! {
    /// Closed destination-admission modes.
    pub enum DestinationAdmissionModeV2 {
        Absent = 1 => "absent",
        PreExistingEmpty = 2 => "pre-existing-empty",
    }
}

closed_u16_catalog! {
    /// Closed semantic root-capability access modes.
    pub enum RootCapabilityAccessV2 {
        ReadOnlyInput = 1 => "read-only-input",
        DurableOutput = 2 => "durable-output",
    }
}

closed_u16_catalog! {
    /// Closed semantic root-capability rights.
    pub enum RootCapabilityRightV2 {
        Traverse = 1 => "traverse",
        ReadObject = 2 => "read-object",
        Enumerate = 3 => "enumerate",
        CreateObject = 4 => "create-object",
        PopulateEmptyDestination = 5 => "populate-empty-destination",
        SyncObject = 6 => "sync-object",
        SyncContainer = 7 => "sync-container",
        AcquireExclusiveLease = 8 => "acquire-exclusive-lease",
        QueryGeneration = 9 => "query-generation",
        CommitCompareAndSwap = 10 => "commit-compare-and-swap",
    }
}

closed_u16_catalog! {
    /// Closed overlap relation carried by the base root-policy registry.
    pub enum OverlapPolicyRelationV2 {
        RequireInputOutputDisjoint = 1 => "require-input-output-disjoint",
    }
}

/// One outer-tag descriptor for a registered-payload catalog.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct TaggedCatalogDescriptorV2 {
    tag: u16,
    name: &'static str,
    requires_registered_id: bool,
}

impl TaggedCatalogDescriptorV2 {
    const fn new(tag: u16, name: &'static str, requires_registered_id: bool) -> Self {
        Self {
            tag,
            name,
            requires_registered_id,
        }
    }

    /// Exact unsigned 16-bit outer tag.
    #[must_use]
    pub const fn tag(self) -> u16 {
        self.tag
    }

    /// Stable lowercase variant name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        self.name
    }

    /// Whether the variant carries one nonzero registered identifier.
    #[must_use]
    pub const fn requires_registered_id(self) -> bool {
        self.requires_registered_id
    }
}

/// A tagged catalog value could not be constructed canonically.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CatalogValueError {
    /// The outer unsigned 16-bit tag is not in the closed catalog.
    UnknownCode(UnknownCatalogCode),
    /// A registered-payload variant omitted its identifier.
    MissingRegisteredId {
        /// Rust type name of the owning catalog.
        catalog: &'static str,
        /// Outer registered-payload tag.
        tag: u16,
    },
    /// A fixed variant was incorrectly accompanied by a registered identifier.
    UnexpectedRegisteredId {
        /// Rust type name of the owning catalog.
        catalog: &'static str,
        /// Fixed outer tag.
        tag: u16,
        /// Unexpected identifier.
        registered_id: u16,
    },
    /// Registered identifiers are nonzero.
    ZeroRegisteredId {
        /// Rust type name of the owning catalog.
        catalog: &'static str,
        /// Outer registered-payload tag.
        tag: u16,
    },
}

impl fmt::Display for CatalogValueError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnknownCode(error) => fmt::Display::fmt(error, f),
            Self::MissingRegisteredId { catalog, tag } => {
                write!(f, "{catalog} tag {tag} requires one registered id")
            }
            Self::UnexpectedRegisteredId {
                catalog,
                tag,
                registered_id,
            } => write!(
                f,
                "{catalog} fixed tag {tag} forbids registered id {registered_id}"
            ),
            Self::ZeroRegisteredId { catalog, tag } => {
                write!(f, "{catalog} tag {tag} refuses registered id zero")
            }
        }
    }
}

impl std::error::Error for CatalogValueError {}

macro_rules! registered_payload_catalog {
    (
        $(#[$enum_meta:meta])*
        pub enum $catalog:ident {
            $($fixed_variant:ident = $fixed_tag:literal => $fixed_name:literal,)+
            @ $registered_variant:ident = $registered_tag:literal
                => $registered_name:literal $(,)?
        }
    ) => {
        $(#[$enum_meta])*
        #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub enum $catalog {
            $(
                #[doc = concat!("Stable catalog value `", $fixed_name, "`.")]
                $fixed_variant,
            )+
            #[doc = concat!(
                "Stable registered-payload catalog value `",
                $registered_name,
                "`."
            )]
            $registered_variant(NonZeroU16),
        }

        impl $catalog {
            /// Every outer tag, in frozen wire order.
            pub const ALL: [
                TaggedCatalogDescriptorV2;
                count_variants!($($fixed_variant),+) + 1
            ] = [
                $(
                    TaggedCatalogDescriptorV2::new(
                        $fixed_tag,
                        $fixed_name,
                        false,
                    ),
                )+
                TaggedCatalogDescriptorV2::new(
                    $registered_tag,
                    $registered_name,
                    true,
                ),
            ];

            /// Decode a fixed variant without a registered payload.
            ///
            /// # Errors
            ///
            /// Registered-payload tags return
            /// [`CatalogValueError::MissingRegisteredId`]; unknown tags return
            /// [`CatalogValueError::UnknownCode`].
            pub const fn from_code(code: u16) -> Result<Self, CatalogValueError> {
                Self::from_tag(code, None)
            }

            /// Decode an outer tag and its optional registered identifier.
            ///
            /// # Errors
            ///
            /// Refuses unknown tags, missing or zero registered identifiers,
            /// and identifiers attached to fixed variants.
            pub const fn from_tag(
                tag: u16,
                registered_id: Option<u16>,
            ) -> Result<Self, CatalogValueError> {
                match tag {
                    $(
                        $fixed_tag => match registered_id {
                            None => Ok(Self::$fixed_variant),
                            Some(registered_id) => Err(
                                CatalogValueError::UnexpectedRegisteredId {
                                    catalog: stringify!($catalog),
                                    tag,
                                    registered_id,
                                },
                            ),
                        },
                    )+
                    $registered_tag => match registered_id {
                        None => Err(CatalogValueError::MissingRegisteredId {
                            catalog: stringify!($catalog),
                            tag,
                        }),
                        Some(0) => Err(CatalogValueError::ZeroRegisteredId {
                            catalog: stringify!($catalog),
                            tag,
                        }),
                        Some(registered_id) => {
                            match NonZeroU16::new(registered_id) {
                                Some(id) => Ok(Self::$registered_variant(id)),
                                None => Err(CatalogValueError::ZeroRegisteredId {
                                    catalog: stringify!($catalog),
                                    tag,
                                }),
                            }
                        }
                    },
                    _ => Err(CatalogValueError::UnknownCode(
                        UnknownCatalogCode::new(stringify!($catalog), tag),
                    )),
                }
            }

            /// Exact unsigned 16-bit outer tag.
            #[must_use]
            pub const fn tag(self) -> u16 {
                match self {
                    $(Self::$fixed_variant => $fixed_tag,)+
                    Self::$registered_variant(_) => $registered_tag,
                }
            }

            /// Stable lowercase variant name.
            #[must_use]
            pub const fn name(self) -> &'static str {
                match self {
                    $(Self::$fixed_variant => $fixed_name,)+
                    Self::$registered_variant(_) => $registered_name,
                }
            }

            /// Registered payload identifier, present only on its designated
            /// variant.
            #[must_use]
            pub const fn registered_id(self) -> Option<u16> {
                match self {
                    $(Self::$fixed_variant => None,)+
                    Self::$registered_variant(id) => Some(id.get()),
                }
            }
        }
    };
}

registered_payload_catalog! {
    /// Root class plus the registered policy id required by `Other`.
    pub enum RootClassV2 {
        InputArtifactRoot = 1 => "input-artifact-root",
        OutputArtifactRoot = 2 => "output-artifact-root",
        @ Other = 3 => "other",
    }
}

registered_payload_catalog! {
    /// Logical unit plus the registered unit id required by `RegisteredUnit`.
    pub enum LogicalUnitV2 {
        EncodedBytes = 1 => "encoded-bytes",
        ExpandedBytes = 2 => "expanded-bytes",
        StoredBytes = 3 => "stored-bytes",
        LogicalBytes = 4 => "logical-bytes",
        Count = 5 => "count",
        Records = 6 => "records",
        Rows = 7 => "rows",
        Elements = 8 => "elements",
        Samples = 9 => "samples",
        Iterations = 10 => "iterations",
        Operations = 11 => "operations",
        Cycles = 12 => "cycles",
        Nanoseconds = 13 => "nanoseconds",
        Seconds = 14 => "seconds",
        Dimensionless = 15 => "dimensionless",
        @ RegisteredUnit = 16 => "registered-unit",
    }
}

registered_payload_catalog! {
    /// Artifact role plus the family role id required by the extension case.
    pub enum ArtifactRoleV2 {
        Observation = 1 => "observation",
        ComparisonDetail = 2 => "comparison-detail",
        EffectDetail = 3 => "effect-detail",
        DiagnosticLog = 4 => "diagnostic-log",
        FamilyEvidence = 5 => "family-evidence",
        PerformanceEvidence = 6 => "performance-evidence",
        ReplaySupport = 7 => "replay-support",
        @ RegisteredFamilyRole = 8 => "registered-family-role",
    }
}

registered_payload_catalog! {
    /// Logical extent axis plus the registered axis id required by extensions.
    pub enum LogicalExtentAxisV2 {
        Payload = 1 => "payload",
        Records = 2 => "records",
        Rows = 3 => "rows",
        Elements = 4 => "elements",
        Samples = 5 => "samples",
        Iterations = 6 => "iterations",
        Operations = 7 => "operations",
        Cycles = 8 => "cycles",
        Duration = 9 => "duration",
        @ RegisteredAxis = 10 => "registered-axis",
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;

    use super::*;

    fn assert_fieldless_catalog<T>(
        catalog: &'static str,
        actual: &[T],
        expected: &[(T, u16, &'static str)],
        code: impl Fn(T) -> u16,
        name: impl Fn(T) -> &'static str,
        decode: impl Fn(u16) -> Result<T, UnknownCatalogCode>,
        unknown_codes: &[u16],
    ) where
        T: Copy + fmt::Debug + Eq + Ord,
    {
        assert_eq!(
            core::mem::size_of::<T>(),
            core::mem::size_of::<u16>(),
            "{catalog} must retain an exact u16 discriminant"
        );
        assert_eq!(actual.len(), expected.len(), "{catalog} count drift");
        for (index, ((expected_value, expected_code, expected_name), actual_value)) in
            expected.iter().zip(actual).enumerate()
        {
            assert_eq!(
                actual_value, expected_value,
                "{catalog} order drift at {index}"
            );
            assert_eq!(
                code(*actual_value),
                *expected_code,
                "{catalog} code drift at {index}"
            );
            assert_eq!(
                name(*actual_value),
                *expected_name,
                "{catalog} name drift at {index}"
            );
            assert_eq!(
                decode(*expected_code),
                Ok(*expected_value),
                "{catalog} roundtrip drift at {index}"
            );
        }

        let codes = expected
            .iter()
            .map(|(_, code, _)| *code)
            .collect::<BTreeSet<_>>();
        let names = expected
            .iter()
            .map(|(_, _, name)| *name)
            .collect::<BTreeSet<_>>();
        assert_eq!(codes.len(), expected.len(), "{catalog} duplicate code");
        assert_eq!(names.len(), expected.len(), "{catalog} duplicate name");

        for &unknown in unknown_codes {
            let error = decode(unknown).expect_err("unknown code must refuse");
            assert_eq!(error.catalog(), catalog);
            assert_eq!(error.code(), unknown);
        }
    }

    fn assert_tagged_catalog<T>(
        catalog: &'static str,
        actual: &[TaggedCatalogDescriptorV2],
        expected: &[(u16, &'static str, bool)],
        from_code: impl Fn(u16) -> Result<T, CatalogValueError>,
        from_tag: impl Fn(u16, Option<u16>) -> Result<T, CatalogValueError>,
        tag: impl Fn(T) -> u16,
        name: impl Fn(T) -> &'static str,
        registered_id: impl Fn(T) -> Option<u16>,
        fixed_values: &[(T, u16, &'static str)],
        registered_tag: u16,
        registered_name: &'static str,
        unknown: u16,
    ) where
        T: Copy + fmt::Debug + Eq,
    {
        assert_eq!(actual.len(), expected.len(), "{catalog} count drift");
        for (index, (descriptor, (expected_tag, expected_name, requires_id))) in
            actual.iter().zip(expected).enumerate()
        {
            assert_eq!(
                descriptor.tag(),
                *expected_tag,
                "{catalog} tag drift at {index}"
            );
            assert_eq!(
                descriptor.name(),
                *expected_name,
                "{catalog} name drift at {index}"
            );
            assert_eq!(
                descriptor.requires_registered_id(),
                *requires_id,
                "{catalog} payload rule drift at {index}"
            );
        }

        for &(value, expected_tag, expected_name) in fixed_values {
            assert_eq!(tag(value), expected_tag);
            assert_eq!(name(value), expected_name);
            assert_eq!(registered_id(value), None);
            assert_eq!(from_code(expected_tag), Ok(value));
            assert_eq!(from_tag(expected_tag, None), Ok(value));
            assert!(matches!(
                from_tag(expected_tag, Some(7)),
                Err(CatalogValueError::UnexpectedRegisteredId {
                    catalog: observed_catalog,
                    tag: observed_tag,
                    registered_id: 7,
                }) if observed_catalog == catalog && observed_tag == expected_tag
            ));
        }

        assert!(matches!(
            from_code(registered_tag),
            Err(CatalogValueError::MissingRegisteredId {
                catalog: observed_catalog,
                tag: observed_tag,
            }) if observed_catalog == catalog && observed_tag == registered_tag
        ));
        assert!(matches!(
            from_tag(registered_tag, Some(0)),
            Err(CatalogValueError::ZeroRegisteredId {
                catalog: observed_catalog,
                tag: observed_tag,
            }) if observed_catalog == catalog && observed_tag == registered_tag
        ));
        let registered = from_tag(registered_tag, Some(u16::MAX))
            .expect("nonzero registered id must be retained");
        assert_eq!(tag(registered), registered_tag);
        assert_eq!(name(registered), registered_name);
        assert_eq!(registered_id(registered), Some(u16::MAX));
        assert!(matches!(
            from_tag(unknown, Some(7)),
            Err(CatalogValueError::UnknownCode(error))
                if error.catalog() == catalog && error.code() == unknown
        ));

        let tags = expected
            .iter()
            .map(|(tag, _, _)| *tag)
            .collect::<BTreeSet<_>>();
        let names = expected
            .iter()
            .map(|(_, name, _)| *name)
            .collect::<BTreeSet<_>>();
        assert_eq!(tags.len(), expected.len(), "{catalog} duplicate tag");
        assert_eq!(names.len(), expected.len(), "{catalog} duplicate name");
    }

    #[test]
    fn api_wire_and_predecessor_literal_oracle() {
        let expected_api = [(RunnerApiGeneration::V2, 2, "RunnerSpecV2")];
        assert_eq!(RunnerApiGeneration::ALL, [expected_api[0].0]);
        assert_eq!(RUNNER_SPEC_V2_API_GENERATION, expected_api[0].0);
        assert_eq!(expected_api[0].0.code(), expected_api[0].1);
        assert_eq!(expected_api[0].0.name(), expected_api[0].2);
        assert_eq!(
            RunnerApiGeneration::from_code(expected_api[0].1),
            Ok(expected_api[0].0)
        );
        for unknown in [0, 1, 3, u16::MAX] {
            let error = RunnerApiGeneration::from_code(unknown)
                .expect_err("unknown API generation must refuse");
            assert_eq!(error.catalog(), "RunnerApiGeneration");
            assert_eq!(error.code(), unknown);
        }

        let expected_wire = [(RunnerWireVersion::V1, 1, "runner-wire-v1")];
        assert_eq!(RunnerWireVersion::ALL, [expected_wire[0].0]);
        assert_eq!(RUNNER_V2_WIRE_VERSION, expected_wire[0].0);
        assert_eq!(expected_wire[0].0.code(), expected_wire[0].1);
        assert_eq!(expected_wire[0].0.name(), expected_wire[0].2);
        assert_eq!(
            RunnerWireVersion::from_code(expected_wire[0].1),
            Ok(expected_wire[0].0)
        );
        for unknown in [0, 2, u16::MAX] {
            let error = RunnerWireVersion::from_code(unknown)
                .expect_err("unknown wire version must refuse");
            assert_eq!(error.catalog(), "RunnerWireVersion");
            assert_eq!(error.code(), unknown);
        }

        let expected_predecessor = [(
            WirePredecessorPolicyV1::NoPredecessor,
            "no-predecessor",
            None,
        )];
        assert_eq!(WirePredecessorPolicyV1::ALL, [expected_predecessor[0].0]);
        assert_eq!(RUNNER_V2_PREDECESSOR_POLICY, expected_predecessor[0].0);
        assert_eq!(expected_predecessor[0].0.name(), expected_predecessor[0].1);
        assert_eq!(
            expected_predecessor[0].0.predecessor(),
            expected_predecessor[0].2
        );
        assert_eq!(
            RunnerWireVersion::V1.predecessor_policy(),
            WirePredecessorPolicyV1::NoPredecessor
        );
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn terminal_command_and_diagnostic_literal_oracles() {
        let proof_exit = [
            (ProofExitV2::Pass, 0, "pass"),
            (ProofExitV2::Failed, 10, "failed"),
            (ProofExitV2::Refused, 11, "refused"),
            (ProofExitV2::NoData, 12, "no-data"),
            (ProofExitV2::Stale, 13, "stale"),
            (ProofExitV2::EnvironmentInvalid, 14, "environment-invalid"),
            (ProofExitV2::Blocked, 15, "blocked"),
            (ProofExitV2::Unsupported, 16, "unsupported"),
            (ProofExitV2::NotRun, 17, "not-run"),
            (ProofExitV2::Cancelled, 18, "cancelled"),
            (ProofExitV2::TimedOut, 19, "timed-out"),
            (ProofExitV2::Usage, 64, "usage"),
            (ProofExitV2::InternalError, 70, "internal-error"),
        ];
        assert_fieldless_catalog(
            "ProofExitV2",
            &ProofExitV2::ALL,
            &proof_exit,
            ProofExitV2::code,
            ProofExitV2::name,
            ProofExitV2::from_code,
            &[1, 9, 20, 63, 65, 69, 71, u16::MAX],
        );

        let refused_reason = [
            (RefusedReasonV2::InvalidEvidence, 1, "invalid-evidence"),
            (
                RefusedReasonV2::NonCanonicalEvidence,
                2,
                "non-canonical-evidence",
            ),
            (
                RefusedReasonV2::EvidenceIdentityMismatch,
                3,
                "evidence-identity-mismatch",
            ),
            (RefusedReasonV2::EvidenceTampered, 4, "evidence-tampered"),
            (RefusedReasonV2::LimitExceeded, 5, "limit-exceeded"),
            (
                RefusedReasonV2::UnsafeArtifactPlacement,
                6,
                "unsafe-artifact-placement",
            ),
            (RefusedReasonV2::ArtifactCollision, 7, "artifact-collision"),
            (
                RefusedReasonV2::LifecycleViolation,
                8,
                "lifecycle-violation",
            ),
            (RefusedReasonV2::PolicyRefused, 9, "policy-refused"),
            (
                RefusedReasonV2::AuthorityBoundaryViolation,
                10,
                "authority-boundary-violation",
            ),
            (RefusedReasonV2::MigrationRefused, 11, "migration-refused"),
        ];
        assert_fieldless_catalog(
            "RefusedReasonV2",
            &RefusedReasonV2::ALL,
            &refused_reason,
            RefusedReasonV2::code,
            RefusedReasonV2::name,
            RefusedReasonV2::from_code,
            &[0, 12, u16::MAX],
        );

        let command = [
            (RunnerCommandV2::List, 0, "list"),
            (RunnerCommandV2::Check, 1, "check"),
            (RunnerCommandV2::SelfTest, 2, "self-test"),
            (RunnerCommandV2::Run, 3, "run"),
            (RunnerCommandV2::Negative, 4, "negative"),
            (RunnerCommandV2::Replay, 5, "replay"),
        ];
        assert_fieldless_catalog(
            "RunnerCommandV2",
            &RunnerCommandV2::ALL,
            &command,
            RunnerCommandV2::code,
            RunnerCommandV2::name,
            RunnerCommandV2::from_code,
            &[6, u16::MAX],
        );

        let run_profile = [
            (RunProfileV2::Smoke, 1, "smoke"),
            (RunProfileV2::Full, 2, "full"),
        ];
        assert_fieldless_catalog(
            "RunProfileV2",
            &RunProfileV2::ALL,
            &run_profile,
            RunProfileV2::code,
            RunProfileV2::name,
            RunProfileV2::from_code,
            &[0, 3, u16::MAX],
        );

        let disposition = [
            (
                ArtifactDispositionV2::LifecycleOnlyNoBundle,
                1,
                "lifecycle-only-no-bundle",
            ),
            (
                ArtifactDispositionV2::DurableBundleRequired,
                2,
                "durable-bundle-required",
            ),
        ];
        assert_fieldless_catalog(
            "ArtifactDispositionV2",
            &ArtifactDispositionV2::ALL,
            &disposition,
            ArtifactDispositionV2::code,
            ArtifactDispositionV2::name,
            ArtifactDispositionV2::from_code,
            &[0, 3, u16::MAX],
        );

        let path_profile = [
            (
                PlatformPathProfileV2::PosixDescriptorRelativeV1,
                1,
                "posix-descriptor-relative-v1",
            ),
            (
                PlatformPathProfileV2::WindowsHandleRelativeV1,
                2,
                "windows-handle-relative-v1",
            ),
            (
                PlatformPathProfileV2::ContentStoreObjectKeyV1,
                3,
                "content-store-object-key-v1",
            ),
        ];
        assert_fieldless_catalog(
            "PlatformPathProfileV2",
            &PlatformPathProfileV2::ALL,
            &path_profile,
            PlatformPathProfileV2::code,
            PlatformPathProfileV2::name,
            PlatformPathProfileV2::from_code,
            &[0, 4, u16::MAX],
        );

        let record_kind = [
            (LifecycleRecordKindV2::RunStart, 1, "run-start"),
            (LifecycleRecordKindV2::CaseStart, 2, "case-start"),
            (LifecycleRecordKindV2::FamilyRow, 3, "family-row"),
            (LifecycleRecordKindV2::CaseTerminal, 4, "case-terminal"),
            (LifecycleRecordKindV2::RunSummary, 5, "run-summary"),
            (LifecycleRecordKindV2::RunTerminal, 6, "run-terminal"),
        ];
        assert_fieldless_catalog(
            "LifecycleRecordKindV2",
            &LifecycleRecordKindV2::ALL,
            &record_kind,
            LifecycleRecordKindV2::code,
            LifecycleRecordKindV2::name,
            LifecycleRecordKindV2::from_code,
            &[0, 7, u16::MAX],
        );

        let record_role = [
            (
                StateBearingRecordRoleV2::PreRunDiagnostic,
                1,
                "pre-run-diagnostic",
            ),
            (
                StateBearingRecordRoleV2::ExecutedCaseTerminal,
                2,
                "executed-case-terminal",
            ),
            (
                StateBearingRecordRoleV2::SuppressedCaseTerminal,
                3,
                "suppressed-case-terminal",
            ),
            (StateBearingRecordRoleV2::RunTerminal, 4, "run-terminal"),
        ];
        assert_fieldless_catalog(
            "StateBearingRecordRoleV2",
            &StateBearingRecordRoleV2::ALL,
            &record_role,
            StateBearingRecordRoleV2::code,
            StateBearingRecordRoleV2::name,
            StateBearingRecordRoleV2::from_code,
            &[0, 5, u16::MAX],
        );

        let diagnostic = [
            (
                DiagnosticCodeV2::CaseConformanceMismatch,
                1,
                "case.conformance_mismatch",
            ),
            (DiagnosticCodeV2::RunnerNotRun, 2, "runner.not_run"),
            (DiagnosticCodeV2::RunnerRefused, 3, "runner.refused"),
            (DiagnosticCodeV2::RunnerNoData, 4, "runner.no_data"),
            (DiagnosticCodeV2::RunnerStale, 5, "runner.stale"),
            (
                DiagnosticCodeV2::RunnerEnvironmentInvalid,
                6,
                "runner.environment_invalid",
            ),
            (DiagnosticCodeV2::RunnerBlocked, 7, "runner.blocked"),
            (DiagnosticCodeV2::RunnerUnsupported, 8, "runner.unsupported"),
            (DiagnosticCodeV2::RunnerCancelled, 9, "runner.cancelled"),
            (DiagnosticCodeV2::RunnerTimedOut, 10, "runner.timed_out"),
            (DiagnosticCodeV2::RunnerUsage, 11, "runner.usage"),
            (
                DiagnosticCodeV2::RunnerInternalError,
                12,
                "runner.internal_error",
            ),
        ];
        assert_fieldless_catalog(
            "DiagnosticCodeV2",
            &DiagnosticCodeV2::ALL,
            &diagnostic,
            DiagnosticCodeV2::code,
            DiagnosticCodeV2::name,
            DiagnosticCodeV2::from_code,
            &[0, 13, u16::MAX],
        );

        let retryability = [
            (RetryabilityV2::Never, 0, "never"),
            (RetryabilityV2::SameInvocation, 1, "same-invocation"),
            (RetryabilityV2::AfterInputChange, 2, "after-input-change"),
            (
                RetryabilityV2::AfterEnvironmentChange,
                3,
                "after-environment-change",
            ),
            (
                RetryabilityV2::AfterPrerequisiteChange,
                4,
                "after-prerequisite-change",
            ),
        ];
        assert_fieldless_catalog(
            "RetryabilityV2",
            &RetryabilityV2::ALL,
            &retryability,
            RetryabilityV2::code,
            RetryabilityV2::name,
            RetryabilityV2::from_code,
            &[5, u16::MAX],
        );

        let repair = [
            (RepairActionKindV2::ChangeArguments, 1, "change-arguments"),
            (RepairActionKindV2::SupplyEvidence, 2, "supply-evidence"),
            (
                RepairActionKindV2::RegenerateCanonicalEvidence,
                3,
                "regenerate-canonical-evidence",
            ),
            (RepairActionKindV2::RefreshEvidence, 4, "refresh-evidence"),
            (
                RepairActionKindV2::ReduceResourceDemand,
                5,
                "reduce-resource-demand",
            ),
            (
                RepairActionKindV2::ChooseSafeArtifactDestination,
                6,
                "choose-safe-artifact-destination",
            ),
            (RepairActionKindV2::RestoreLifecycle, 7, "restore-lifecycle"),
            (
                RepairActionKindV2::UpdatePolicyOrCapability,
                8,
                "update-policy-or-capability",
            ),
            (
                RepairActionKindV2::RegisterMigration,
                9,
                "register-migration",
            ),
            (
                RepairActionKindV2::RetrySameInvocation,
                10,
                "retry-same-invocation",
            ),
            (RepairActionKindV2::ContactOwner, 11, "contact-owner"),
            (
                RepairActionKindV2::InspectRetainedArtifact,
                12,
                "inspect-retained-artifact",
            ),
        ];
        assert_fieldless_catalog(
            "RepairActionKindV2",
            &RepairActionKindV2::ALL,
            &repair,
            RepairActionKindV2::code,
            RepairActionKindV2::name,
            RepairActionKindV2::from_code,
            &[0, 13, u16::MAX],
        );

        let not_run_cause = [
            (NotRunCauseCodeV2::PriorCancelled, 1, "prior-cancelled"),
            (NotRunCauseCodeV2::PriorTimedOut, 2, "prior-timed-out"),
            (
                NotRunCauseCodeV2::PriorControlledInternalError,
                3,
                "prior-controlled-internal-error",
            ),
        ];
        assert_fieldless_catalog(
            "NotRunCauseCodeV2",
            &NotRunCauseCodeV2::ALL,
            &not_run_cause,
            NotRunCauseCodeV2::code,
            NotRunCauseCodeV2::name,
            NotRunCauseCodeV2::from_code,
            &[0, 4, u16::MAX],
        );
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn value_digest_and_role_literal_oracles() {
        let typed_value = [
            (TypedValueTagV2::I8, 1, "i8"),
            (TypedValueTagV2::I16, 2, "i16"),
            (TypedValueTagV2::I32, 3, "i32"),
            (TypedValueTagV2::I64, 4, "i64"),
            (TypedValueTagV2::I128, 5, "i128"),
            (TypedValueTagV2::U8, 6, "u8"),
            (TypedValueTagV2::U16, 7, "u16"),
            (TypedValueTagV2::U32, 8, "u32"),
            (TypedValueTagV2::U64, 9, "u64"),
            (TypedValueTagV2::U128, 10, "u128"),
            (TypedValueTagV2::Rational, 11, "rational"),
            (TypedValueTagV2::Decimal, 12, "decimal"),
            (TypedValueTagV2::F32Bits, 13, "f32-bits"),
            (TypedValueTagV2::F64Bits, 14, "f64-bits"),
            (TypedValueTagV2::Digest, 15, "digest"),
            (TypedValueTagV2::Quantity, 16, "quantity"),
            (TypedValueTagV2::Token, 17, "token"),
            (TypedValueTagV2::Text, 18, "text"),
            (TypedValueTagV2::RelativePath, 19, "relative-path"),
            (TypedValueTagV2::OpaqueBytes, 20, "opaque-bytes"),
        ];
        assert_fieldless_catalog(
            "TypedValueTagV2",
            &TypedValueTagV2::ALL,
            &typed_value,
            TypedValueTagV2::code,
            TypedValueTagV2::name,
            TypedValueTagV2::from_code,
            &[0, 21, u16::MAX],
        );

        let option = [
            (TypedOptionTagV1::Absent, 0, "absent"),
            (TypedOptionTagV1::Present, 1, "present"),
        ];
        assert_fieldless_catalog(
            "TypedOptionTagV1",
            &TypedOptionTagV1::ALL,
            &option,
            TypedOptionTagV1::code,
            TypedOptionTagV1::name,
            TypedOptionTagV1::from_code,
            &[2, u16::MAX],
        );

        let digest_role = [
            (DigestRoleV2::Spec, 1, "spec"),
            (DigestRoleV2::Invocation, 2, "invocation"),
            (DigestRoleV2::Run, 3, "run"),
            (DigestRoleV2::Source, 4, "source"),
            (DigestRoleV2::Build, 5, "build"),
            (DigestRoleV2::Toolchain, 6, "toolchain"),
            (DigestRoleV2::CaseManifest, 7, "case-manifest"),
            (DigestRoleV2::ArtifactEncoded, 8, "artifact-encoded"),
            (DigestRoleV2::ArtifactContent, 9, "artifact-content"),
            (DigestRoleV2::StoredObject, 10, "stored-object"),
            (DigestRoleV2::ArtifactInventory, 11, "artifact-inventory"),
            (DigestRoleV2::LifecycleLog, 12, "lifecycle-log"),
            (DigestRoleV2::RunSummary, 13, "run-summary"),
            (DigestRoleV2::RunTerminal, 14, "run-terminal"),
            (DigestRoleV2::BundleManifest, 15, "bundle-manifest"),
            (DigestRoleV2::DurablePublication, 16, "durable-publication"),
            (DigestRoleV2::Seal, 17, "seal"),
            (
                DigestRoleV2::PublishedBundleReceipt,
                18,
                "published-bundle-receipt",
            ),
            (DigestRoleV2::Policy, 19, "policy"),
            (DigestRoleV2::CandidateBytes, 20, "candidate-bytes"),
            (DigestRoleV2::CandidateSchema, 21, "candidate-schema"),
            (DigestRoleV2::SourceClosure, 22, "source-closure"),
            (DigestRoleV2::ClaimScope, 23, "claim-scope"),
            (DigestRoleV2::ProducerManifest, 24, "producer-manifest"),
            (
                DigestRoleV2::RegisteredFamilyDomain,
                25,
                "registered-family-domain",
            ),
        ];
        assert_fieldless_catalog(
            "DigestRoleV2",
            &DigestRoleV2::ALL,
            &digest_role,
            DigestRoleV2::code,
            DigestRoleV2::name,
            DigestRoleV2::from_code,
            &[0, 26, u16::MAX],
        );
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn publication_capability_and_extension_literal_oracles() {
        let protocol = [
            (
                PublicationProtocolV2::PosixDescriptorRenameAndDirectorySyncV1,
                1,
                "posix-descriptor-rename-and-directory-sync-v1",
            ),
            (
                PublicationProtocolV2::WindowsHandleReplaceAndDirectoryFlushV1,
                2,
                "windows-handle-replace-and-directory-flush-v1",
            ),
            (
                PublicationProtocolV2::ContentStoreAtomicCommitV1,
                3,
                "content-store-atomic-commit-v1",
            ),
        ];
        assert_fieldless_catalog(
            "PublicationProtocolV2",
            &PublicationProtocolV2::ALL,
            &protocol,
            PublicationProtocolV2::code,
            PublicationProtocolV2::name,
            PublicationProtocolV2::from_code,
            &[0, 4, u16::MAX],
        );

        let destination = [
            (DestinationAdmissionModeV2::Absent, 1, "absent"),
            (
                DestinationAdmissionModeV2::PreExistingEmpty,
                2,
                "pre-existing-empty",
            ),
        ];
        assert_fieldless_catalog(
            "DestinationAdmissionModeV2",
            &DestinationAdmissionModeV2::ALL,
            &destination,
            DestinationAdmissionModeV2::code,
            DestinationAdmissionModeV2::name,
            DestinationAdmissionModeV2::from_code,
            &[0, 3, u16::MAX],
        );

        let access = [
            (RootCapabilityAccessV2::ReadOnlyInput, 1, "read-only-input"),
            (RootCapabilityAccessV2::DurableOutput, 2, "durable-output"),
        ];
        assert_fieldless_catalog(
            "RootCapabilityAccessV2",
            &RootCapabilityAccessV2::ALL,
            &access,
            RootCapabilityAccessV2::code,
            RootCapabilityAccessV2::name,
            RootCapabilityAccessV2::from_code,
            &[0, 3, u16::MAX],
        );

        let rights = [
            (RootCapabilityRightV2::Traverse, 1, "traverse"),
            (RootCapabilityRightV2::ReadObject, 2, "read-object"),
            (RootCapabilityRightV2::Enumerate, 3, "enumerate"),
            (RootCapabilityRightV2::CreateObject, 4, "create-object"),
            (
                RootCapabilityRightV2::PopulateEmptyDestination,
                5,
                "populate-empty-destination",
            ),
            (RootCapabilityRightV2::SyncObject, 6, "sync-object"),
            (RootCapabilityRightV2::SyncContainer, 7, "sync-container"),
            (
                RootCapabilityRightV2::AcquireExclusiveLease,
                8,
                "acquire-exclusive-lease",
            ),
            (
                RootCapabilityRightV2::QueryGeneration,
                9,
                "query-generation",
            ),
            (
                RootCapabilityRightV2::CommitCompareAndSwap,
                10,
                "commit-compare-and-swap",
            ),
        ];
        assert_fieldless_catalog(
            "RootCapabilityRightV2",
            &RootCapabilityRightV2::ALL,
            &rights,
            RootCapabilityRightV2::code,
            RootCapabilityRightV2::name,
            RootCapabilityRightV2::from_code,
            &[0, 11, u16::MAX],
        );

        let overlap = [(
            OverlapPolicyRelationV2::RequireInputOutputDisjoint,
            1,
            "require-input-output-disjoint",
        )];
        assert_fieldless_catalog(
            "OverlapPolicyRelationV2",
            &OverlapPolicyRelationV2::ALL,
            &overlap,
            OverlapPolicyRelationV2::code,
            OverlapPolicyRelationV2::name,
            OverlapPolicyRelationV2::from_code,
            &[0, 2, u16::MAX],
        );

        let root_class_descriptors = [
            (1, "input-artifact-root", false),
            (2, "output-artifact-root", false),
            (3, "other", true),
        ];
        let root_class_fixed = [
            (RootClassV2::InputArtifactRoot, 1, "input-artifact-root"),
            (RootClassV2::OutputArtifactRoot, 2, "output-artifact-root"),
        ];
        assert_tagged_catalog(
            "RootClassV2",
            &RootClassV2::ALL,
            &root_class_descriptors,
            RootClassV2::from_code,
            RootClassV2::from_tag,
            RootClassV2::tag,
            RootClassV2::name,
            RootClassV2::registered_id,
            &root_class_fixed,
            3,
            "other",
            4,
        );

        let logical_unit_descriptors = [
            (1, "encoded-bytes", false),
            (2, "expanded-bytes", false),
            (3, "stored-bytes", false),
            (4, "logical-bytes", false),
            (5, "count", false),
            (6, "records", false),
            (7, "rows", false),
            (8, "elements", false),
            (9, "samples", false),
            (10, "iterations", false),
            (11, "operations", false),
            (12, "cycles", false),
            (13, "nanoseconds", false),
            (14, "seconds", false),
            (15, "dimensionless", false),
            (16, "registered-unit", true),
        ];
        let logical_unit_fixed = [
            (LogicalUnitV2::EncodedBytes, 1, "encoded-bytes"),
            (LogicalUnitV2::ExpandedBytes, 2, "expanded-bytes"),
            (LogicalUnitV2::StoredBytes, 3, "stored-bytes"),
            (LogicalUnitV2::LogicalBytes, 4, "logical-bytes"),
            (LogicalUnitV2::Count, 5, "count"),
            (LogicalUnitV2::Records, 6, "records"),
            (LogicalUnitV2::Rows, 7, "rows"),
            (LogicalUnitV2::Elements, 8, "elements"),
            (LogicalUnitV2::Samples, 9, "samples"),
            (LogicalUnitV2::Iterations, 10, "iterations"),
            (LogicalUnitV2::Operations, 11, "operations"),
            (LogicalUnitV2::Cycles, 12, "cycles"),
            (LogicalUnitV2::Nanoseconds, 13, "nanoseconds"),
            (LogicalUnitV2::Seconds, 14, "seconds"),
            (LogicalUnitV2::Dimensionless, 15, "dimensionless"),
        ];
        assert_tagged_catalog(
            "LogicalUnitV2",
            &LogicalUnitV2::ALL,
            &logical_unit_descriptors,
            LogicalUnitV2::from_code,
            LogicalUnitV2::from_tag,
            LogicalUnitV2::tag,
            LogicalUnitV2::name,
            LogicalUnitV2::registered_id,
            &logical_unit_fixed,
            16,
            "registered-unit",
            17,
        );

        let artifact_role_descriptors = [
            (1, "observation", false),
            (2, "comparison-detail", false),
            (3, "effect-detail", false),
            (4, "diagnostic-log", false),
            (5, "family-evidence", false),
            (6, "performance-evidence", false),
            (7, "replay-support", false),
            (8, "registered-family-role", true),
        ];
        let artifact_role_fixed = [
            (ArtifactRoleV2::Observation, 1, "observation"),
            (ArtifactRoleV2::ComparisonDetail, 2, "comparison-detail"),
            (ArtifactRoleV2::EffectDetail, 3, "effect-detail"),
            (ArtifactRoleV2::DiagnosticLog, 4, "diagnostic-log"),
            (ArtifactRoleV2::FamilyEvidence, 5, "family-evidence"),
            (
                ArtifactRoleV2::PerformanceEvidence,
                6,
                "performance-evidence",
            ),
            (ArtifactRoleV2::ReplaySupport, 7, "replay-support"),
        ];
        assert_tagged_catalog(
            "ArtifactRoleV2",
            &ArtifactRoleV2::ALL,
            &artifact_role_descriptors,
            ArtifactRoleV2::from_code,
            ArtifactRoleV2::from_tag,
            ArtifactRoleV2::tag,
            ArtifactRoleV2::name,
            ArtifactRoleV2::registered_id,
            &artifact_role_fixed,
            8,
            "registered-family-role",
            9,
        );

        let axis_descriptors = [
            (1, "payload", false),
            (2, "records", false),
            (3, "rows", false),
            (4, "elements", false),
            (5, "samples", false),
            (6, "iterations", false),
            (7, "operations", false),
            (8, "cycles", false),
            (9, "duration", false),
            (10, "registered-axis", true),
        ];
        let axis_fixed = [
            (LogicalExtentAxisV2::Payload, 1, "payload"),
            (LogicalExtentAxisV2::Records, 2, "records"),
            (LogicalExtentAxisV2::Rows, 3, "rows"),
            (LogicalExtentAxisV2::Elements, 4, "elements"),
            (LogicalExtentAxisV2::Samples, 5, "samples"),
            (LogicalExtentAxisV2::Iterations, 6, "iterations"),
            (LogicalExtentAxisV2::Operations, 7, "operations"),
            (LogicalExtentAxisV2::Cycles, 8, "cycles"),
            (LogicalExtentAxisV2::Duration, 9, "duration"),
        ];
        assert_tagged_catalog(
            "LogicalExtentAxisV2",
            &LogicalExtentAxisV2::ALL,
            &axis_descriptors,
            LogicalExtentAxisV2::from_code,
            LogicalExtentAxisV2::from_tag,
            LogicalExtentAxisV2::tag,
            LogicalExtentAxisV2::name,
            LogicalExtentAxisV2::registered_id,
            &axis_fixed,
            10,
            "registered-axis",
            11,
        );
    }
}
