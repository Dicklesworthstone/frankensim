//! Frozen Runner V2 invocation budgets.
//!
//! Construction is deliberately two-stage. `RunnerBudgetsV2::try_new` checks
//! intrinsic relations that do not require command context. `admit` then checks
//! the selected profile, artifact disposition, and admitted Runner V2 limits.
//! Neither stage performs I/O or grants execution, storage, or admission
//! authority.

use crate::catalog::{
    ArtifactDispositionV2, DigestRoleV2, LogicalUnitV2, RepairActionKindV2, RunProfileV2,
};
use crate::extension::BaseExtensionRegistryProjectionV2;
use crate::identity::{DigestValueV2, RunnerBudgetsRootV2};
use crate::limits::RunnerLimitsV2;
use fs_blake3::{ContentHash, hash_domain};

const GIB: u64 = 1024 * 1024 * 1024;
const SECOND_NS: u64 = 1_000_000_000;

/// Canonical identity domain for the ordered budget-value projection.
pub const RUNNER_BUDGETS_PROJECTION_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.runner-budgets.v1";

/// Canonical domain binding a registered-unit budget to the exact extension
/// registry that established membership.
pub const RUNNER_BUDGETS_EXTENSION_REGISTRY_BINDING_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.runner-budgets-extension-registry-binding.v1";

/// Exact number of fields in the Runner V2 budget schema.
pub const RUNNER_BUDGET_FIELD_COUNT_V2: usize = 18;

/// One exact field in the ordered Runner V2 budget schema.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum RunnerBudgetFieldV2 {
    /// Total wall-clock allowance in nanoseconds.
    WallTimeNs = 1,
    /// Maximum resident memory in logical bytes.
    MaxResidentBytes = 2,
    /// Maximum number of child processes over the invocation.
    MaxChildProcesses = 3,
    /// Maximum number of child processes active concurrently.
    MaxParallelChildren = 4,
    /// Exact logical-work allowance.
    LogicalWorkLimit = 5,
    /// Unit that interprets the logical-work allowance.
    LogicalWorkUnit = 6,
    /// Aggregate lifecycle-document encoded bytes.
    LifecycleEncodedBytes = 7,
    /// Atomic command-result stdout encoded bytes.
    CommandResultStdoutBytes = 8,
    /// Aggregate child stdout encoded bytes.
    CombinedChildStdoutBytes = 9,
    /// Aggregate child stderr encoded bytes.
    CombinedChildStderrBytes = 10,
    /// Aggregate encoded bytes of published artifacts.
    ArtifactEncodedBytes = 11,
    /// Aggregate stored bytes of published artifacts.
    ArtifactStoredBytes = 12,
    /// Aggregate expanded bytes of published artifacts.
    ArtifactExpandedBytes = 13,
    /// Aggregate stored bytes of the six system publication objects.
    SystemPublicationStoredBytes = 14,
    /// Whole-publication stored bytes.
    PublicationStoredBytes = 15,
    /// Stop-observation timeout in nanoseconds.
    StopObservationNs = 16,
    /// Cancellation-drain timeout in nanoseconds.
    DrainNs = 17,
    /// Finalization timeout in nanoseconds.
    FinalizeNs = 18,
}

impl RunnerBudgetFieldV2 {
    /// Exact canonical field order.
    pub const ALL: [Self; RUNNER_BUDGET_FIELD_COUNT_V2] = [
        Self::WallTimeNs,
        Self::MaxResidentBytes,
        Self::MaxChildProcesses,
        Self::MaxParallelChildren,
        Self::LogicalWorkLimit,
        Self::LogicalWorkUnit,
        Self::LifecycleEncodedBytes,
        Self::CommandResultStdoutBytes,
        Self::CombinedChildStdoutBytes,
        Self::CombinedChildStderrBytes,
        Self::ArtifactEncodedBytes,
        Self::ArtifactStoredBytes,
        Self::ArtifactExpandedBytes,
        Self::SystemPublicationStoredBytes,
        Self::PublicationStoredBytes,
        Self::StopObservationNs,
        Self::DrainNs,
        Self::FinalizeNs,
    ];

    /// Exact one-based schema ordinal.
    #[must_use]
    pub const fn ordinal(self) -> u16 {
        self as u16
    }

    /// Resolve an exact one-based schema ordinal.
    #[must_use]
    pub const fn from_ordinal(ordinal: u16) -> Option<Self> {
        match ordinal {
            1 => Some(Self::WallTimeNs),
            2 => Some(Self::MaxResidentBytes),
            3 => Some(Self::MaxChildProcesses),
            4 => Some(Self::MaxParallelChildren),
            5 => Some(Self::LogicalWorkLimit),
            6 => Some(Self::LogicalWorkUnit),
            7 => Some(Self::LifecycleEncodedBytes),
            8 => Some(Self::CommandResultStdoutBytes),
            9 => Some(Self::CombinedChildStdoutBytes),
            10 => Some(Self::CombinedChildStderrBytes),
            11 => Some(Self::ArtifactEncodedBytes),
            12 => Some(Self::ArtifactStoredBytes),
            13 => Some(Self::ArtifactExpandedBytes),
            14 => Some(Self::SystemPublicationStoredBytes),
            15 => Some(Self::PublicationStoredBytes),
            16 => Some(Self::StopObservationNs),
            17 => Some(Self::DrainNs),
            18 => Some(Self::FinalizeNs),
            _ => None,
        }
    }

    /// Static descriptor for this exact field.
    #[must_use]
    pub const fn descriptor(self) -> &'static RunnerBudgetDescriptorV2 {
        &RUNNER_BUDGET_DESCRIPTORS_V2[(self as usize) - 1]
    }
}

/// Frozen primitive or tagged-sum width of one budget field.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum RunnerBudgetWidthV2 {
    /// Unsigned 32-bit integer.
    U32,
    /// Unsigned 64-bit integer.
    U64,
    /// Unsigned 128-bit integer.
    U128,
    /// Closed logical-unit tagged sum with an optional registered identifier.
    LogicalUnitTaggedSum,
}

/// Unit carried by a budget field and by its deterministic refusals.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum RunnerBudgetUnitV2 {
    /// Nanoseconds.
    Nanoseconds,
    /// Logical bytes before storage-envelope interpretation.
    LogicalBytes,
    /// Discrete object or process count.
    Count,
    /// Exact logical-work magnitude.
    LogicalWork,
    /// Closed tag interpreting logical work.
    LogicalWorkUnit,
    /// Canonically encoded bytes.
    EncodedBytes,
    /// Bytes charged to storage.
    StoredBytes,
    /// Bytes after deterministic expansion.
    ExpandedBytes,
}

/// Static descriptor for one ordered budget field.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct RunnerBudgetDescriptorV2 {
    /// Typed field identity.
    pub field: RunnerBudgetFieldV2,
    /// Exact one-based canonical position.
    pub ordinal: u16,
    /// Stable snake-case field name.
    pub name: &'static str,
    /// Frozen primitive or tagged-sum width.
    pub width: RunnerBudgetWidthV2,
    /// Semantic unit.
    pub unit: RunnerBudgetUnitV2,
}

/// Exact ordered Runner V2 budget descriptors.
pub const RUNNER_BUDGET_DESCRIPTORS_V2: [RunnerBudgetDescriptorV2; RUNNER_BUDGET_FIELD_COUNT_V2] = [
    RunnerBudgetDescriptorV2 {
        field: RunnerBudgetFieldV2::WallTimeNs,
        ordinal: 1,
        name: "wall_time_ns",
        width: RunnerBudgetWidthV2::U64,
        unit: RunnerBudgetUnitV2::Nanoseconds,
    },
    RunnerBudgetDescriptorV2 {
        field: RunnerBudgetFieldV2::MaxResidentBytes,
        ordinal: 2,
        name: "max_resident_bytes",
        width: RunnerBudgetWidthV2::U64,
        unit: RunnerBudgetUnitV2::LogicalBytes,
    },
    RunnerBudgetDescriptorV2 {
        field: RunnerBudgetFieldV2::MaxChildProcesses,
        ordinal: 3,
        name: "max_child_processes",
        width: RunnerBudgetWidthV2::U32,
        unit: RunnerBudgetUnitV2::Count,
    },
    RunnerBudgetDescriptorV2 {
        field: RunnerBudgetFieldV2::MaxParallelChildren,
        ordinal: 4,
        name: "max_parallel_children",
        width: RunnerBudgetWidthV2::U32,
        unit: RunnerBudgetUnitV2::Count,
    },
    RunnerBudgetDescriptorV2 {
        field: RunnerBudgetFieldV2::LogicalWorkLimit,
        ordinal: 5,
        name: "logical_work_limit",
        width: RunnerBudgetWidthV2::U128,
        unit: RunnerBudgetUnitV2::LogicalWork,
    },
    RunnerBudgetDescriptorV2 {
        field: RunnerBudgetFieldV2::LogicalWorkUnit,
        ordinal: 6,
        name: "logical_work_unit",
        width: RunnerBudgetWidthV2::LogicalUnitTaggedSum,
        unit: RunnerBudgetUnitV2::LogicalWorkUnit,
    },
    RunnerBudgetDescriptorV2 {
        field: RunnerBudgetFieldV2::LifecycleEncodedBytes,
        ordinal: 7,
        name: "lifecycle_encoded_bytes",
        width: RunnerBudgetWidthV2::U64,
        unit: RunnerBudgetUnitV2::EncodedBytes,
    },
    RunnerBudgetDescriptorV2 {
        field: RunnerBudgetFieldV2::CommandResultStdoutBytes,
        ordinal: 8,
        name: "command_result_stdout_bytes",
        width: RunnerBudgetWidthV2::U64,
        unit: RunnerBudgetUnitV2::EncodedBytes,
    },
    RunnerBudgetDescriptorV2 {
        field: RunnerBudgetFieldV2::CombinedChildStdoutBytes,
        ordinal: 9,
        name: "combined_child_stdout_bytes",
        width: RunnerBudgetWidthV2::U64,
        unit: RunnerBudgetUnitV2::EncodedBytes,
    },
    RunnerBudgetDescriptorV2 {
        field: RunnerBudgetFieldV2::CombinedChildStderrBytes,
        ordinal: 10,
        name: "combined_child_stderr_bytes",
        width: RunnerBudgetWidthV2::U64,
        unit: RunnerBudgetUnitV2::EncodedBytes,
    },
    RunnerBudgetDescriptorV2 {
        field: RunnerBudgetFieldV2::ArtifactEncodedBytes,
        ordinal: 11,
        name: "artifact_encoded_bytes",
        width: RunnerBudgetWidthV2::U64,
        unit: RunnerBudgetUnitV2::EncodedBytes,
    },
    RunnerBudgetDescriptorV2 {
        field: RunnerBudgetFieldV2::ArtifactStoredBytes,
        ordinal: 12,
        name: "artifact_stored_bytes",
        width: RunnerBudgetWidthV2::U64,
        unit: RunnerBudgetUnitV2::StoredBytes,
    },
    RunnerBudgetDescriptorV2 {
        field: RunnerBudgetFieldV2::ArtifactExpandedBytes,
        ordinal: 13,
        name: "artifact_expanded_bytes",
        width: RunnerBudgetWidthV2::U64,
        unit: RunnerBudgetUnitV2::ExpandedBytes,
    },
    RunnerBudgetDescriptorV2 {
        field: RunnerBudgetFieldV2::SystemPublicationStoredBytes,
        ordinal: 14,
        name: "system_publication_stored_bytes",
        width: RunnerBudgetWidthV2::U64,
        unit: RunnerBudgetUnitV2::StoredBytes,
    },
    RunnerBudgetDescriptorV2 {
        field: RunnerBudgetFieldV2::PublicationStoredBytes,
        ordinal: 15,
        name: "publication_stored_bytes",
        width: RunnerBudgetWidthV2::U64,
        unit: RunnerBudgetUnitV2::StoredBytes,
    },
    RunnerBudgetDescriptorV2 {
        field: RunnerBudgetFieldV2::StopObservationNs,
        ordinal: 16,
        name: "stop_observation_ns",
        width: RunnerBudgetWidthV2::U64,
        unit: RunnerBudgetUnitV2::Nanoseconds,
    },
    RunnerBudgetDescriptorV2 {
        field: RunnerBudgetFieldV2::DrainNs,
        ordinal: 17,
        name: "drain_ns",
        width: RunnerBudgetWidthV2::U64,
        unit: RunnerBudgetUnitV2::Nanoseconds,
    },
    RunnerBudgetDescriptorV2 {
        field: RunnerBudgetFieldV2::FinalizeNs,
        ordinal: 18,
        name: "finalize_ns",
        width: RunnerBudgetWidthV2::U64,
        unit: RunnerBudgetUnitV2::Nanoseconds,
    },
];

/// Mutable, unadmitted input for intrinsic budget construction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RunnerBudgetsCandidateV2 {
    /// Candidate total wall-clock allowance in nanoseconds.
    pub wall_time_ns: u64,
    /// Candidate maximum resident memory in logical bytes.
    pub max_resident_bytes: u64,
    /// Candidate total child-process count.
    pub max_child_processes: u32,
    /// Candidate concurrent child-process count.
    pub max_parallel_children: u32,
    /// Candidate exact logical-work allowance.
    pub logical_work_limit: u128,
    /// Candidate logical-work unit.
    pub logical_work_unit: LogicalUnitV2,
    /// Candidate lifecycle-document encoded-byte allowance.
    pub lifecycle_encoded_bytes: u64,
    /// Candidate atomic command-result stdout allowance.
    pub command_result_stdout_bytes: u64,
    /// Candidate aggregate child-stdout allowance.
    pub combined_child_stdout_bytes: u64,
    /// Candidate aggregate child-stderr allowance.
    pub combined_child_stderr_bytes: u64,
    /// Candidate artifact encoded-byte allowance.
    pub artifact_encoded_bytes: u64,
    /// Candidate artifact stored-byte allowance.
    pub artifact_stored_bytes: u64,
    /// Candidate artifact expanded-byte allowance.
    pub artifact_expanded_bytes: u64,
    /// Candidate system-object stored-byte allowance.
    pub system_publication_stored_bytes: u64,
    /// Candidate whole-publication stored-byte allowance.
    pub publication_stored_bytes: u64,
    /// Candidate stop-observation timeout in nanoseconds.
    pub stop_observation_ns: u64,
    /// Candidate cancellation-drain timeout in nanoseconds.
    pub drain_ns: u64,
    /// Candidate finalization timeout in nanoseconds.
    pub finalize_ns: u64,
}

/// Intrinsically valid, immutable Runner V2 budgets.
///
/// Admitted values are inspected through typed, read-only accessors:
///
/// ```
/// use fs_evidence_runner::RunnerBudgetsV2;
///
/// fn wall_time_ns(budgets: &RunnerBudgetsV2) -> u64 {
///     budgets.wall_time_ns()
/// }
/// ```
///
/// A caller cannot post-mutate or widen an intrinsically validated grant:
///
/// ```compile_fail
/// use fs_evidence_runner::RunnerBudgetsV2;
///
/// fn widen_wall_time(budgets: &mut RunnerBudgetsV2) {
///     budgets.wall_time_ns = u64::MAX;
/// }
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RunnerBudgetsV2 {
    wall_time_ns: u64,
    max_resident_bytes: u64,
    max_child_processes: u32,
    max_parallel_children: u32,
    logical_work_limit: u128,
    logical_work_unit: LogicalUnitV2,
    lifecycle_encoded_bytes: u64,
    command_result_stdout_bytes: u64,
    combined_child_stdout_bytes: u64,
    combined_child_stderr_bytes: u64,
    artifact_encoded_bytes: u64,
    artifact_stored_bytes: u64,
    artifact_expanded_bytes: u64,
    system_publication_stored_bytes: u64,
    publication_stored_bytes: u64,
    stop_observation_ns: u64,
    drain_ns: u64,
    finalize_ns: u64,
}

/// Intrinsically valid budgets whose registered logical-work unit is bound to
/// the exact extension registry that established membership.
///
/// The nested budget retains the frozen 18-field schema unchanged. The
/// separate binding root prevents the same syntactic registered ID from being
/// replayed under a different descriptor set.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RegistryBoundRunnerBudgetsV2 {
    budgets: RunnerBudgetsV2,
    extension_registry_root: ContentHash,
    binding_root: ContentHash,
}

/// Contextually admitted budgets. This is still non-authoritative schema data:
/// it proves only profile/disposition/limit consistency.
///
/// The nested immutable budget remains available through a read-only
/// accessor:
///
/// ```
/// use fs_evidence_runner::AdmittedRunnerBudgetsV2;
///
/// fn wall_time_ns(admitted: &AdmittedRunnerBudgetsV2) -> u64 {
///     admitted.budgets().wall_time_ns()
/// }
/// ```
///
/// Contextual admission does not expose a cap-widening mutation path:
///
/// ```compile_fail
/// use fs_evidence_runner::AdmittedRunnerBudgetsV2;
///
/// fn widen_admitted_wall_time(admitted: &mut AdmittedRunnerBudgetsV2) {
///     admitted.budgets.wall_time_ns = u64::MAX;
/// }
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AdmittedRunnerBudgetsV2 {
    budgets: RunnerBudgetsV2,
    profile: RunProfileV2,
    disposition: ArtifactDispositionV2,
}

/// Contextually admitted registered-unit budgets retaining their exact
/// extension-registry membership binding.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RegistryBoundAdmittedRunnerBudgetsV2 {
    admitted: AdmittedRunnerBudgetsV2,
    extension_registry_root: ContentHash,
    binding_root: ContentHash,
}

/// Heterogeneous value retained by a budget refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum RunnerBudgetValueV2 {
    /// Unsigned 32-bit value.
    U32(u32),
    /// Unsigned 64-bit value.
    U64(u64),
    /// Unsigned 128-bit value.
    U128(u128),
    /// Logical-unit tag and its optional registered identifier.
    LogicalUnit {
        /// Frozen closed-catalog or registered-unit tag.
        tag: u16,
        /// Registered identifier, present only for the registered-unit tag.
        registered_id: Option<u16>,
    },
}

/// Exact expectation retained by a budget refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum RunnerBudgetExpectationV2 {
    /// The value must be nonzero.
    NonZero,
    /// The value must be zero.
    Zero,
    /// The value must not exceed the carried ceiling.
    AtMost(RunnerBudgetValueV2),
    /// The value must meet or exceed the carried floor.
    AtLeast(RunnerBudgetValueV2),
    /// The value must equal the carried value.
    Exactly(RunnerBudgetValueV2),
    /// A registered logical unit must resolve in the exact bound extension
    /// registry.
    RegisteredInExtensionRegistry,
}

/// Deterministic class of a budget refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum RunnerBudgetViolationKindV2 {
    /// An intrinsically required allowance was zero.
    Zero,
    /// Concurrent children exceeded total children.
    ParallelChildrenExceedTotal,
    /// Checked addition of timeout components overflowed.
    TimeoutSumOverflow,
    /// Timeout components exceeded wall time.
    TimeoutSumExceedsWall,
    /// A profile-owned ceiling was exceeded.
    ProfileCeilingExceeded,
    /// An admitted Runner limit was exceeded.
    LimitExceeded,
    /// The selected context requires zero.
    ContextualZeroRequired,
    /// The selected context requires a nonzero value.
    ContextualNonZeroRequired,
    /// Atomic stdout could not contain the lifecycle document.
    CommandResultCannotContainLifecycle,
    /// Stored artifact bytes were below encoded artifact bytes.
    ArtifactStoredBelowEncoded,
    /// Checked whole-publication addition overflowed.
    PublicationSumOverflow,
    /// Whole-publication bytes disagreed with the exact component sum.
    PublicationEquationMismatch,
    /// A syntactic registered logical-work unit was not resolved through an
    /// exact extension registry.
    UnregisteredLogicalWorkUnit,
}

/// Precise, bounded budget refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct RunnerBudgetViolationV2 {
    kind: RunnerBudgetViolationKindV2,
    field: RunnerBudgetFieldV2,
    unit: RunnerBudgetUnitV2,
    expected: RunnerBudgetExpectationV2,
    observed: RunnerBudgetValueV2,
}

impl RunnerBudgetViolationV2 {
    fn new(
        kind: RunnerBudgetViolationKindV2,
        field: RunnerBudgetFieldV2,
        expected: RunnerBudgetExpectationV2,
        observed: RunnerBudgetValueV2,
    ) -> Self {
        Self {
            kind,
            field,
            unit: field.descriptor().unit,
            expected,
            observed,
        }
    }

    /// Stable refusal class.
    #[must_use]
    pub const fn kind(&self) -> RunnerBudgetViolationKindV2 {
        self.kind
    }

    /// Budget field that refused admission.
    #[must_use]
    pub const fn field(&self) -> RunnerBudgetFieldV2 {
        self.field
    }

    /// Semantic unit of the expected and observed values.
    #[must_use]
    pub const fn unit(&self) -> RunnerBudgetUnitV2 {
        self.unit
    }

    /// Exact predicate required for admission.
    #[must_use]
    pub const fn expected(&self) -> RunnerBudgetExpectationV2 {
        self.expected
    }

    /// Exact value that violated the predicate.
    #[must_use]
    pub const fn observed(&self) -> RunnerBudgetValueV2 {
        self.observed
    }

    /// Stable owner used by structured diagnostics.
    #[must_use]
    pub const fn owner(&self) -> &'static str {
        "fs-evidence-runner.runner-budgets"
    }

    /// One-based rank of the primary bounded repair recommendation.
    #[must_use]
    pub const fn repair_rank(&self) -> u8 {
        1
    }

    /// Closed non-executable repair class appropriate to this refusal.
    #[must_use]
    pub const fn repair_kind(&self) -> RepairActionKindV2 {
        match self.expected {
            RunnerBudgetExpectationV2::AtMost(_) => RepairActionKindV2::ReduceResourceDemand,
            RunnerBudgetExpectationV2::RegisteredInExtensionRegistry => {
                RepairActionKindV2::UpdatePolicyOrCapability
            }
            RunnerBudgetExpectationV2::NonZero
            | RunnerBudgetExpectationV2::Zero
            | RunnerBudgetExpectationV2::AtLeast(_)
            | RunnerBudgetExpectationV2::Exactly(_) => RepairActionKindV2::ChangeArguments,
        }
    }

    /// Stable structured repair target; this is data, not an executable command.
    #[must_use]
    pub const fn repair_target(&self) -> &'static str {
        match self.kind {
            RunnerBudgetViolationKindV2::TimeoutSumOverflow
            | RunnerBudgetViolationKindV2::TimeoutSumExceedsWall => {
                "stop_observation_ns-or-drain_ns-or-finalize_ns"
            }
            RunnerBudgetViolationKindV2::PublicationSumOverflow => {
                "artifact_stored_bytes-or-system_publication_stored_bytes"
            }
            _ => self.field.descriptor().name,
        }
    }
}

impl RunnerBudgetsV2 {
    /// Construct an intrinsically valid, immutable budget value using only
    /// fixed logical-work units.
    ///
    /// A syntactic `LogicalUnitV2::RegisteredUnit` cannot establish registry
    /// membership on this registry-free path.
    pub fn try_new(candidate: RunnerBudgetsCandidateV2) -> Result<Self, RunnerBudgetViolationV2> {
        Self::try_new_intrinsic(candidate, None)
    }

    /// Construct intrinsically valid budgets while resolving a registered
    /// logical-work unit through one exact extension-registry projection.
    pub fn try_new_with_extension_registry(
        candidate: RunnerBudgetsCandidateV2,
        registry: &BaseExtensionRegistryProjectionV2,
    ) -> Result<RegistryBoundRunnerBudgetsV2, RunnerBudgetViolationV2> {
        let budgets = Self::try_new_intrinsic(candidate, Some(registry))?;
        let extension_registry_root = *registry.root();
        let binding_root =
            runner_budgets_extension_registry_binding_root(&budgets, extension_registry_root);
        Ok(RegistryBoundRunnerBudgetsV2 {
            budgets,
            extension_registry_root,
            binding_root,
        })
    }

    fn try_new_intrinsic(
        candidate: RunnerBudgetsCandidateV2,
        registry: Option<&BaseExtensionRegistryProjectionV2>,
    ) -> Result<Self, RunnerBudgetViolationV2> {
        if let Some(registered_id) = candidate.logical_work_unit.registered_id()
            && registry.is_none_or(|registry| registry.logical_unit(registered_id).is_err())
        {
            return Err(RunnerBudgetViolationV2::new(
                RunnerBudgetViolationKindV2::UnregisteredLogicalWorkUnit,
                RunnerBudgetFieldV2::LogicalWorkUnit,
                RunnerBudgetExpectationV2::RegisteredInExtensionRegistry,
                RunnerBudgetValueV2::LogicalUnit {
                    tag: candidate.logical_work_unit.tag(),
                    registered_id: Some(registered_id),
                },
            ));
        }

        for (field, value) in [
            (RunnerBudgetFieldV2::WallTimeNs, candidate.wall_time_ns),
            (
                RunnerBudgetFieldV2::MaxResidentBytes,
                candidate.max_resident_bytes,
            ),
            (
                RunnerBudgetFieldV2::CommandResultStdoutBytes,
                candidate.command_result_stdout_bytes,
            ),
            (
                RunnerBudgetFieldV2::StopObservationNs,
                candidate.stop_observation_ns,
            ),
            (RunnerBudgetFieldV2::DrainNs, candidate.drain_ns),
            (RunnerBudgetFieldV2::FinalizeNs, candidate.finalize_ns),
        ] {
            if value == 0 {
                return Err(RunnerBudgetViolationV2::new(
                    RunnerBudgetViolationKindV2::Zero,
                    field,
                    RunnerBudgetExpectationV2::NonZero,
                    RunnerBudgetValueV2::U64(0),
                ));
            }
        }

        if candidate.max_parallel_children > candidate.max_child_processes {
            return Err(RunnerBudgetViolationV2::new(
                RunnerBudgetViolationKindV2::ParallelChildrenExceedTotal,
                RunnerBudgetFieldV2::MaxParallelChildren,
                RunnerBudgetExpectationV2::AtMost(RunnerBudgetValueV2::U32(
                    candidate.max_child_processes,
                )),
                RunnerBudgetValueV2::U32(candidate.max_parallel_children),
            ));
        }

        let timeout_total_u128 = u128::from(candidate.stop_observation_ns)
            + u128::from(candidate.drain_ns)
            + u128::from(candidate.finalize_ns);
        let timeout_total = u64::try_from(timeout_total_u128).map_err(|_| {
            RunnerBudgetViolationV2::new(
                RunnerBudgetViolationKindV2::TimeoutSumOverflow,
                RunnerBudgetFieldV2::FinalizeNs,
                RunnerBudgetExpectationV2::AtMost(RunnerBudgetValueV2::U128(u128::from(
                    candidate.wall_time_ns,
                ))),
                RunnerBudgetValueV2::U128(timeout_total_u128),
            )
        })?;
        if timeout_total > candidate.wall_time_ns {
            return Err(RunnerBudgetViolationV2::new(
                RunnerBudgetViolationKindV2::TimeoutSumExceedsWall,
                RunnerBudgetFieldV2::FinalizeNs,
                RunnerBudgetExpectationV2::AtMost(RunnerBudgetValueV2::U64(candidate.wall_time_ns)),
                RunnerBudgetValueV2::U64(timeout_total),
            ));
        }

        Ok(Self {
            wall_time_ns: candidate.wall_time_ns,
            max_resident_bytes: candidate.max_resident_bytes,
            max_child_processes: candidate.max_child_processes,
            max_parallel_children: candidate.max_parallel_children,
            logical_work_limit: candidate.logical_work_limit,
            logical_work_unit: candidate.logical_work_unit,
            lifecycle_encoded_bytes: candidate.lifecycle_encoded_bytes,
            command_result_stdout_bytes: candidate.command_result_stdout_bytes,
            combined_child_stdout_bytes: candidate.combined_child_stdout_bytes,
            combined_child_stderr_bytes: candidate.combined_child_stderr_bytes,
            artifact_encoded_bytes: candidate.artifact_encoded_bytes,
            artifact_stored_bytes: candidate.artifact_stored_bytes,
            artifact_expanded_bytes: candidate.artifact_expanded_bytes,
            system_publication_stored_bytes: candidate.system_publication_stored_bytes,
            publication_stored_bytes: candidate.publication_stored_bytes,
            stop_observation_ns: candidate.stop_observation_ns,
            drain_ns: candidate.drain_ns,
            finalize_ns: candidate.finalize_ns,
        })
    }

    /// Contextually admit budgets against profile, disposition, and immutable
    /// limits.
    pub fn admit(
        self,
        profile: RunProfileV2,
        disposition: ArtifactDispositionV2,
        limits: &RunnerLimitsV2,
    ) -> Result<AdmittedRunnerBudgetsV2, RunnerBudgetViolationV2> {
        self.validate_profile_and_process_context(profile, disposition)?;
        self.validate_limit_grants(limits)?;
        self.validate_internal_relations()?;
        self.validate_disposition_and_publication(disposition)?;
        Ok(AdmittedRunnerBudgetsV2 {
            budgets: self,
            profile,
            disposition,
        })
    }

    fn validate_profile_and_process_context(
        &self,
        profile: RunProfileV2,
        disposition: ArtifactDispositionV2,
    ) -> Result<(), RunnerBudgetViolationV2> {
        let (max_wall_time_ns, max_resident_bytes, max_parallel, max_children) =
            profile_ceilings(profile);
        require_at_most_u64(
            RunnerBudgetViolationKindV2::ProfileCeilingExceeded,
            RunnerBudgetFieldV2::WallTimeNs,
            self.wall_time_ns,
            max_wall_time_ns,
        )?;
        require_at_most_u64(
            RunnerBudgetViolationKindV2::ProfileCeilingExceeded,
            RunnerBudgetFieldV2::MaxResidentBytes,
            self.max_resident_bytes,
            max_resident_bytes,
        )?;
        require_at_most_u32(
            RunnerBudgetViolationKindV2::ProfileCeilingExceeded,
            RunnerBudgetFieldV2::MaxChildProcesses,
            self.max_child_processes,
            max_children,
        )?;
        require_at_most_u32(
            RunnerBudgetViolationKindV2::ProfileCeilingExceeded,
            RunnerBudgetFieldV2::MaxParallelChildren,
            self.max_parallel_children,
            max_parallel,
        )?;

        if disposition == ArtifactDispositionV2::DurableBundleRequired
            && self.max_child_processes == 0
        {
            return Err(RunnerBudgetViolationV2::new(
                RunnerBudgetViolationKindV2::ContextualNonZeroRequired,
                RunnerBudgetFieldV2::MaxChildProcesses,
                RunnerBudgetExpectationV2::NonZero,
                RunnerBudgetValueV2::U32(0),
            ));
        }
        Ok(())
    }

    fn validate_limit_grants(
        &self,
        limits: &RunnerLimitsV2,
    ) -> Result<(), RunnerBudgetViolationV2> {
        require_limit(
            RunnerBudgetFieldV2::LifecycleEncodedBytes,
            self.lifecycle_encoded_bytes,
            limits.lifecycle_document_encoded_bytes(),
        )?;
        require_limit(
            RunnerBudgetFieldV2::CommandResultStdoutBytes,
            self.command_result_stdout_bytes,
            limits.command_result_stdout_bytes(),
        )?;
        require_limit(
            RunnerBudgetFieldV2::CombinedChildStdoutBytes,
            self.combined_child_stdout_bytes,
            limits.combined_child_stdout_bytes(),
        )?;
        require_limit(
            RunnerBudgetFieldV2::CombinedChildStderrBytes,
            self.combined_child_stderr_bytes,
            limits.combined_child_stderr_bytes(),
        )?;
        require_limit(
            RunnerBudgetFieldV2::ArtifactEncodedBytes,
            self.artifact_encoded_bytes,
            limits.bundle_encoded_bytes(),
        )?;
        require_limit(
            RunnerBudgetFieldV2::ArtifactStoredBytes,
            self.artifact_stored_bytes,
            limits.artifact_stored_aggregate_bytes(),
        )?;
        require_limit(
            RunnerBudgetFieldV2::ArtifactExpandedBytes,
            self.artifact_expanded_bytes,
            limits.bundle_expanded_bytes(),
        )?;
        require_limit(
            RunnerBudgetFieldV2::SystemPublicationStoredBytes,
            self.system_publication_stored_bytes,
            limits.system_publication_stored_bytes(),
        )?;
        require_limit(
            RunnerBudgetFieldV2::PublicationStoredBytes,
            self.publication_stored_bytes,
            limits.publication_stored_bytes(),
        )?;
        Ok(())
    }

    fn validate_internal_relations(&self) -> Result<(), RunnerBudgetViolationV2> {
        if self.command_result_stdout_bytes < self.lifecycle_encoded_bytes {
            return Err(RunnerBudgetViolationV2::new(
                RunnerBudgetViolationKindV2::CommandResultCannotContainLifecycle,
                RunnerBudgetFieldV2::CommandResultStdoutBytes,
                RunnerBudgetExpectationV2::AtLeast(RunnerBudgetValueV2::U64(
                    self.lifecycle_encoded_bytes,
                )),
                RunnerBudgetValueV2::U64(self.command_result_stdout_bytes),
            ));
        }
        if self.artifact_stored_bytes < self.artifact_encoded_bytes {
            return Err(RunnerBudgetViolationV2::new(
                RunnerBudgetViolationKindV2::ArtifactStoredBelowEncoded,
                RunnerBudgetFieldV2::ArtifactStoredBytes,
                RunnerBudgetExpectationV2::AtLeast(RunnerBudgetValueV2::U64(
                    self.artifact_encoded_bytes,
                )),
                RunnerBudgetValueV2::U64(self.artifact_stored_bytes),
            ));
        }
        Ok(())
    }

    fn validate_disposition_and_publication(
        &self,
        disposition: ArtifactDispositionV2,
    ) -> Result<(), RunnerBudgetViolationV2> {
        match disposition {
            ArtifactDispositionV2::LifecycleOnlyNoBundle => {
                for (field, value) in [
                    (
                        RunnerBudgetFieldV2::ArtifactEncodedBytes,
                        self.artifact_encoded_bytes,
                    ),
                    (
                        RunnerBudgetFieldV2::ArtifactStoredBytes,
                        self.artifact_stored_bytes,
                    ),
                    (
                        RunnerBudgetFieldV2::ArtifactExpandedBytes,
                        self.artifact_expanded_bytes,
                    ),
                    (
                        RunnerBudgetFieldV2::SystemPublicationStoredBytes,
                        self.system_publication_stored_bytes,
                    ),
                    (
                        RunnerBudgetFieldV2::PublicationStoredBytes,
                        self.publication_stored_bytes,
                    ),
                ] {
                    if value != 0 {
                        return Err(RunnerBudgetViolationV2::new(
                            RunnerBudgetViolationKindV2::ContextualZeroRequired,
                            field,
                            RunnerBudgetExpectationV2::Zero,
                            RunnerBudgetValueV2::U64(value),
                        ));
                    }
                }
            }
            ArtifactDispositionV2::DurableBundleRequired => {
                if self.system_publication_stored_bytes == 0 {
                    return Err(RunnerBudgetViolationV2::new(
                        RunnerBudgetViolationKindV2::ContextualNonZeroRequired,
                        RunnerBudgetFieldV2::SystemPublicationStoredBytes,
                        RunnerBudgetExpectationV2::NonZero,
                        RunnerBudgetValueV2::U64(0),
                    ));
                }
            }
        }

        let computed_publication = checked_publication_stored_bytes(
            self.artifact_stored_bytes,
            self.system_publication_stored_bytes,
        )?;
        if self.publication_stored_bytes != computed_publication {
            return Err(RunnerBudgetViolationV2::new(
                RunnerBudgetViolationKindV2::PublicationEquationMismatch,
                RunnerBudgetFieldV2::PublicationStoredBytes,
                RunnerBudgetExpectationV2::Exactly(RunnerBudgetValueV2::U64(computed_publication)),
                RunnerBudgetValueV2::U64(self.publication_stored_bytes),
            ));
        }
        Ok(())
    }

    /// Exact ordered value projection for tests, source closure, and the
    /// nominal RunnerBudgets root owner. Integers are big-endian and
    /// `LogicalUnitV2::RegisteredUnit` alone appends its nonzero ID.
    #[must_use]
    pub fn canonical_projection(&self) -> Vec<u8> {
        let mut bytes = Vec::with_capacity(180);
        bytes.extend_from_slice(b"FSRUNNER-BUDGETS\x01");
        bytes.extend_from_slice(&self.wall_time_ns.to_be_bytes());
        bytes.extend_from_slice(&self.max_resident_bytes.to_be_bytes());
        bytes.extend_from_slice(&self.max_child_processes.to_be_bytes());
        bytes.extend_from_slice(&self.max_parallel_children.to_be_bytes());
        bytes.extend_from_slice(&self.logical_work_limit.to_be_bytes());
        bytes.extend_from_slice(&self.logical_work_unit.tag().to_be_bytes());
        if let Some(registered_id) = self.logical_work_unit.registered_id() {
            bytes.extend_from_slice(&registered_id.to_be_bytes());
        }
        bytes.extend_from_slice(&self.lifecycle_encoded_bytes.to_be_bytes());
        bytes.extend_from_slice(&self.command_result_stdout_bytes.to_be_bytes());
        bytes.extend_from_slice(&self.combined_child_stdout_bytes.to_be_bytes());
        bytes.extend_from_slice(&self.combined_child_stderr_bytes.to_be_bytes());
        bytes.extend_from_slice(&self.artifact_encoded_bytes.to_be_bytes());
        bytes.extend_from_slice(&self.artifact_stored_bytes.to_be_bytes());
        bytes.extend_from_slice(&self.artifact_expanded_bytes.to_be_bytes());
        bytes.extend_from_slice(&self.system_publication_stored_bytes.to_be_bytes());
        bytes.extend_from_slice(&self.publication_stored_bytes.to_be_bytes());
        bytes.extend_from_slice(&self.stop_observation_ns.to_be_bytes());
        bytes.extend_from_slice(&self.drain_ns.to_be_bytes());
        bytes.extend_from_slice(&self.finalize_ns.to_be_bytes());
        bytes
    }

    /// Domain-separated root of `canonical_projection`.
    #[must_use]
    pub fn canonical_projection_root(&self) -> ContentHash {
        hash_domain(
            RUNNER_BUDGETS_PROJECTION_DOMAIN_V1,
            &self.canonical_projection(),
        )
    }

    /// Nominal semantic identity of this exact intrinsically valid budget
    /// vector.
    ///
    /// This is non-authoritative schema identity; contextual admission still
    /// requires [`Self::admit`], and neither form grants execution resources.
    #[must_use]
    pub fn semantic_root(&self) -> RunnerBudgetsRootV2 {
        let content = self.canonical_projection_root();
        let digest = DigestValueV2::from_array(
            DigestRoleV2::Policy,
            RunnerBudgetsRootV2::DESCRIPTOR.domain_witness(),
            *content.as_bytes(),
        );
        RunnerBudgetsRootV2::from_digest(digest)
            .expect("the private budgets constructor fixes the nominal role and domain")
    }

    /// Read one budget value through the exact heterogeneous field catalog.
    #[must_use]
    pub const fn value(&self, field: RunnerBudgetFieldV2) -> RunnerBudgetValueV2 {
        match field {
            RunnerBudgetFieldV2::WallTimeNs => RunnerBudgetValueV2::U64(self.wall_time_ns),
            RunnerBudgetFieldV2::MaxResidentBytes => {
                RunnerBudgetValueV2::U64(self.max_resident_bytes)
            }
            RunnerBudgetFieldV2::MaxChildProcesses => {
                RunnerBudgetValueV2::U32(self.max_child_processes)
            }
            RunnerBudgetFieldV2::MaxParallelChildren => {
                RunnerBudgetValueV2::U32(self.max_parallel_children)
            }
            RunnerBudgetFieldV2::LogicalWorkLimit => {
                RunnerBudgetValueV2::U128(self.logical_work_limit)
            }
            RunnerBudgetFieldV2::LogicalWorkUnit => RunnerBudgetValueV2::LogicalUnit {
                tag: self.logical_work_unit.tag(),
                registered_id: self.logical_work_unit.registered_id(),
            },
            RunnerBudgetFieldV2::LifecycleEncodedBytes => {
                RunnerBudgetValueV2::U64(self.lifecycle_encoded_bytes)
            }
            RunnerBudgetFieldV2::CommandResultStdoutBytes => {
                RunnerBudgetValueV2::U64(self.command_result_stdout_bytes)
            }
            RunnerBudgetFieldV2::CombinedChildStdoutBytes => {
                RunnerBudgetValueV2::U64(self.combined_child_stdout_bytes)
            }
            RunnerBudgetFieldV2::CombinedChildStderrBytes => {
                RunnerBudgetValueV2::U64(self.combined_child_stderr_bytes)
            }
            RunnerBudgetFieldV2::ArtifactEncodedBytes => {
                RunnerBudgetValueV2::U64(self.artifact_encoded_bytes)
            }
            RunnerBudgetFieldV2::ArtifactStoredBytes => {
                RunnerBudgetValueV2::U64(self.artifact_stored_bytes)
            }
            RunnerBudgetFieldV2::ArtifactExpandedBytes => {
                RunnerBudgetValueV2::U64(self.artifact_expanded_bytes)
            }
            RunnerBudgetFieldV2::SystemPublicationStoredBytes => {
                RunnerBudgetValueV2::U64(self.system_publication_stored_bytes)
            }
            RunnerBudgetFieldV2::PublicationStoredBytes => {
                RunnerBudgetValueV2::U64(self.publication_stored_bytes)
            }
            RunnerBudgetFieldV2::StopObservationNs => {
                RunnerBudgetValueV2::U64(self.stop_observation_ns)
            }
            RunnerBudgetFieldV2::DrainNs => RunnerBudgetValueV2::U64(self.drain_ns),
            RunnerBudgetFieldV2::FinalizeNs => RunnerBudgetValueV2::U64(self.finalize_ns),
        }
    }

    /// Recover an explicitly unadmitted candidate.
    #[must_use]
    pub const fn to_candidate(self) -> RunnerBudgetsCandidateV2 {
        RunnerBudgetsCandidateV2 {
            wall_time_ns: self.wall_time_ns,
            max_resident_bytes: self.max_resident_bytes,
            max_child_processes: self.max_child_processes,
            max_parallel_children: self.max_parallel_children,
            logical_work_limit: self.logical_work_limit,
            logical_work_unit: self.logical_work_unit,
            lifecycle_encoded_bytes: self.lifecycle_encoded_bytes,
            command_result_stdout_bytes: self.command_result_stdout_bytes,
            combined_child_stdout_bytes: self.combined_child_stdout_bytes,
            combined_child_stderr_bytes: self.combined_child_stderr_bytes,
            artifact_encoded_bytes: self.artifact_encoded_bytes,
            artifact_stored_bytes: self.artifact_stored_bytes,
            artifact_expanded_bytes: self.artifact_expanded_bytes,
            system_publication_stored_bytes: self.system_publication_stored_bytes,
            publication_stored_bytes: self.publication_stored_bytes,
            stop_observation_ns: self.stop_observation_ns,
            drain_ns: self.drain_ns,
            finalize_ns: self.finalize_ns,
        }
    }

    /// Admitted wall-clock allowance in nanoseconds.
    #[must_use]
    pub const fn wall_time_ns(&self) -> u64 {
        self.wall_time_ns
    }

    /// Admitted maximum resident memory in logical bytes.
    #[must_use]
    pub const fn max_resident_bytes(&self) -> u64 {
        self.max_resident_bytes
    }

    /// Admitted total child-process count.
    #[must_use]
    pub const fn max_child_processes(&self) -> u32 {
        self.max_child_processes
    }

    /// Admitted parallel child-process count.
    #[must_use]
    pub const fn max_parallel_children(&self) -> u32 {
        self.max_parallel_children
    }

    /// Admitted exact logical-work allowance.
    #[must_use]
    pub const fn logical_work_limit(&self) -> u128 {
        self.logical_work_limit
    }

    /// Unit paired with the logical-work allowance.
    #[must_use]
    pub const fn logical_work_unit(&self) -> LogicalUnitV2 {
        self.logical_work_unit
    }

    /// Admitted lifecycle encoded-byte allowance.
    #[must_use]
    pub const fn lifecycle_encoded_bytes(&self) -> u64 {
        self.lifecycle_encoded_bytes
    }

    /// Admitted atomic command-result stdout allowance.
    #[must_use]
    pub const fn command_result_stdout_bytes(&self) -> u64 {
        self.command_result_stdout_bytes
    }

    /// Admitted aggregate child-stdout allowance.
    #[must_use]
    pub const fn combined_child_stdout_bytes(&self) -> u64 {
        self.combined_child_stdout_bytes
    }

    /// Admitted aggregate child-stderr allowance.
    #[must_use]
    pub const fn combined_child_stderr_bytes(&self) -> u64 {
        self.combined_child_stderr_bytes
    }

    /// Admitted artifact encoded-byte allowance.
    #[must_use]
    pub const fn artifact_encoded_bytes(&self) -> u64 {
        self.artifact_encoded_bytes
    }

    /// Admitted artifact stored-byte allowance.
    #[must_use]
    pub const fn artifact_stored_bytes(&self) -> u64 {
        self.artifact_stored_bytes
    }

    /// Admitted artifact expanded-byte allowance.
    #[must_use]
    pub const fn artifact_expanded_bytes(&self) -> u64 {
        self.artifact_expanded_bytes
    }

    /// Admitted stored-byte allowance for the six system objects.
    #[must_use]
    pub const fn system_publication_stored_bytes(&self) -> u64 {
        self.system_publication_stored_bytes
    }

    /// Admitted whole-publication stored-byte allowance.
    #[must_use]
    pub const fn publication_stored_bytes(&self) -> u64 {
        self.publication_stored_bytes
    }

    /// Admitted stop-observation allowance in nanoseconds.
    #[must_use]
    pub const fn stop_observation_ns(&self) -> u64 {
        self.stop_observation_ns
    }

    /// Admitted drain allowance in nanoseconds.
    #[must_use]
    pub const fn drain_ns(&self) -> u64 {
        self.drain_ns
    }

    /// Admitted finalization allowance in nanoseconds.
    #[must_use]
    pub const fn finalize_ns(&self) -> u64 {
        self.finalize_ns
    }
}

impl RegistryBoundRunnerBudgetsV2 {
    /// Read-only access to the unchanged frozen 18-field budget value.
    #[must_use]
    pub const fn budgets(&self) -> &RunnerBudgetsV2 {
        &self.budgets
    }

    /// Exact extension-registry projection root that established membership.
    #[must_use]
    pub const fn extension_registry_root(&self) -> ContentHash {
        self.extension_registry_root
    }

    /// Domain-separated semantic binding of budget and registry roots.
    #[must_use]
    pub const fn semantic_binding_root(&self) -> ContentHash {
        self.binding_root
    }

    /// Contextually admit the bound budgets without discarding the exact
    /// registry-membership witness.
    pub fn admit(
        self,
        profile: RunProfileV2,
        disposition: ArtifactDispositionV2,
        limits: &RunnerLimitsV2,
    ) -> Result<RegistryBoundAdmittedRunnerBudgetsV2, RunnerBudgetViolationV2> {
        Ok(RegistryBoundAdmittedRunnerBudgetsV2 {
            admitted: self.budgets.admit(profile, disposition, limits)?,
            extension_registry_root: self.extension_registry_root,
            binding_root: self.binding_root,
        })
    }
}

impl RegistryBoundAdmittedRunnerBudgetsV2 {
    /// Read-only contextually admitted budget.
    #[must_use]
    pub const fn admitted(&self) -> &AdmittedRunnerBudgetsV2 {
        &self.admitted
    }

    /// Exact extension-registry projection root that established membership.
    #[must_use]
    pub const fn extension_registry_root(&self) -> ContentHash {
        self.extension_registry_root
    }

    /// Domain-separated semantic budget/registry binding.
    #[must_use]
    pub const fn semantic_binding_root(&self) -> ContentHash {
        self.binding_root
    }
}

fn runner_budgets_extension_registry_binding_root(
    budgets: &RunnerBudgetsV2,
    extension_registry_root: ContentHash,
) -> ContentHash {
    let mut bytes = Vec::with_capacity(64);
    bytes.extend_from_slice(budgets.canonical_projection_root().as_bytes());
    bytes.extend_from_slice(extension_registry_root.as_bytes());
    hash_domain(RUNNER_BUDGETS_EXTENSION_REGISTRY_BINDING_DOMAIN_V1, &bytes)
}

impl AdmittedRunnerBudgetsV2 {
    /// Intrinsically valid budgets admitted by this contextual witness.
    #[must_use]
    pub const fn budgets(&self) -> &RunnerBudgetsV2 {
        &self.budgets
    }

    /// Profile whose ceilings were applied.
    #[must_use]
    pub const fn profile(&self) -> RunProfileV2 {
        self.profile
    }

    /// Artifact disposition whose contextual rules were applied.
    #[must_use]
    pub const fn disposition(&self) -> ArtifactDispositionV2 {
        self.disposition
    }
}

fn profile_ceilings(profile: RunProfileV2) -> (u64, u64, u32, u32) {
    match profile {
        RunProfileV2::Smoke => (900 * SECOND_NS, 16 * GIB, 32, 256),
        RunProfileV2::Full => (86_400 * SECOND_NS, 128 * GIB, 64, 256),
    }
}

fn require_at_most_u32(
    kind: RunnerBudgetViolationKindV2,
    field: RunnerBudgetFieldV2,
    observed: u32,
    ceiling: u32,
) -> Result<(), RunnerBudgetViolationV2> {
    if observed <= ceiling {
        Ok(())
    } else {
        Err(RunnerBudgetViolationV2::new(
            kind,
            field,
            RunnerBudgetExpectationV2::AtMost(RunnerBudgetValueV2::U32(ceiling)),
            RunnerBudgetValueV2::U32(observed),
        ))
    }
}

fn require_at_most_u64(
    kind: RunnerBudgetViolationKindV2,
    field: RunnerBudgetFieldV2,
    observed: u64,
    ceiling: u64,
) -> Result<(), RunnerBudgetViolationV2> {
    if observed <= ceiling {
        Ok(())
    } else {
        Err(RunnerBudgetViolationV2::new(
            kind,
            field,
            RunnerBudgetExpectationV2::AtMost(RunnerBudgetValueV2::U64(ceiling)),
            RunnerBudgetValueV2::U64(observed),
        ))
    }
}

fn require_limit(
    field: RunnerBudgetFieldV2,
    observed: u64,
    ceiling: u64,
) -> Result<(), RunnerBudgetViolationV2> {
    require_at_most_u64(
        RunnerBudgetViolationKindV2::LimitExceeded,
        field,
        observed,
        ceiling,
    )
}

/// Checked whole-publication grant equation used by contextual admission.
pub fn checked_publication_stored_bytes(
    artifact_stored_bytes: u64,
    system_publication_stored_bytes: u64,
) -> Result<u64, RunnerBudgetViolationV2> {
    let exact_sum = u128::from(artifact_stored_bytes) + u128::from(system_publication_stored_bytes);
    u64::try_from(exact_sum).map_err(|_| {
        RunnerBudgetViolationV2::new(
            RunnerBudgetViolationKindV2::PublicationSumOverflow,
            RunnerBudgetFieldV2::PublicationStoredBytes,
            RunnerBudgetExpectationV2::AtMost(RunnerBudgetValueV2::U128(u128::from(u64::MAX))),
            RunnerBudgetValueV2::U128(exact_sum),
        )
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::catalog::ArtifactRoleV2;
    use crate::extension::{
        BaseExtensionRegistryProjectionV2, RegisteredArtifactRoleDescriptorV2,
        RegisteredLogicalUnitDescriptorV2,
    };
    use crate::identity::NoClaimScopeRootV1;
    use crate::value::StableTokenV2;

    const TEST_MIB: u64 = 1024 * 1024;

    const EXPECTED_NAMES: [&str; RUNNER_BUDGET_FIELD_COUNT_V2] = [
        "wall_time_ns",
        "max_resident_bytes",
        "max_child_processes",
        "max_parallel_children",
        "logical_work_limit",
        "logical_work_unit",
        "lifecycle_encoded_bytes",
        "command_result_stdout_bytes",
        "combined_child_stdout_bytes",
        "combined_child_stderr_bytes",
        "artifact_encoded_bytes",
        "artifact_stored_bytes",
        "artifact_expanded_bytes",
        "system_publication_stored_bytes",
        "publication_stored_bytes",
        "stop_observation_ns",
        "drain_ns",
        "finalize_ns",
    ];

    const EXPECTED_WIDTHS: [RunnerBudgetWidthV2; RUNNER_BUDGET_FIELD_COUNT_V2] = [
        RunnerBudgetWidthV2::U64,
        RunnerBudgetWidthV2::U64,
        RunnerBudgetWidthV2::U32,
        RunnerBudgetWidthV2::U32,
        RunnerBudgetWidthV2::U128,
        RunnerBudgetWidthV2::LogicalUnitTaggedSum,
        RunnerBudgetWidthV2::U64,
        RunnerBudgetWidthV2::U64,
        RunnerBudgetWidthV2::U64,
        RunnerBudgetWidthV2::U64,
        RunnerBudgetWidthV2::U64,
        RunnerBudgetWidthV2::U64,
        RunnerBudgetWidthV2::U64,
        RunnerBudgetWidthV2::U64,
        RunnerBudgetWidthV2::U64,
        RunnerBudgetWidthV2::U64,
        RunnerBudgetWidthV2::U64,
        RunnerBudgetWidthV2::U64,
    ];

    const EXPECTED_UNITS: [RunnerBudgetUnitV2; RUNNER_BUDGET_FIELD_COUNT_V2] = [
        RunnerBudgetUnitV2::Nanoseconds,
        RunnerBudgetUnitV2::LogicalBytes,
        RunnerBudgetUnitV2::Count,
        RunnerBudgetUnitV2::Count,
        RunnerBudgetUnitV2::LogicalWork,
        RunnerBudgetUnitV2::LogicalWorkUnit,
        RunnerBudgetUnitV2::EncodedBytes,
        RunnerBudgetUnitV2::EncodedBytes,
        RunnerBudgetUnitV2::EncodedBytes,
        RunnerBudgetUnitV2::EncodedBytes,
        RunnerBudgetUnitV2::EncodedBytes,
        RunnerBudgetUnitV2::StoredBytes,
        RunnerBudgetUnitV2::ExpandedBytes,
        RunnerBudgetUnitV2::StoredBytes,
        RunnerBudgetUnitV2::StoredBytes,
        RunnerBudgetUnitV2::Nanoseconds,
        RunnerBudgetUnitV2::Nanoseconds,
        RunnerBudgetUnitV2::Nanoseconds,
    ];

    const EXPECTED_DURABLE_VALUES: [RunnerBudgetValueV2; RUNNER_BUDGET_FIELD_COUNT_V2] = [
        RunnerBudgetValueV2::U64(100_000_000_000),
        RunnerBudgetValueV2::U64(1_073_741_824),
        RunnerBudgetValueV2::U32(8),
        RunnerBudgetValueV2::U32(4),
        RunnerBudgetValueV2::U128(1_000),
        RunnerBudgetValueV2::LogicalUnit {
            tag: 11,
            registered_id: None,
        },
        RunnerBudgetValueV2::U64(1_000),
        RunnerBudgetValueV2::U64(4_000),
        RunnerBudgetValueV2::U64(2_000),
        RunnerBudgetValueV2::U64(1_000),
        RunnerBudgetValueV2::U64(100),
        RunnerBudgetValueV2::U64(104),
        RunnerBudgetValueV2::U64(200),
        RunnerBudgetValueV2::U64(72),
        RunnerBudgetValueV2::U64(176),
        RunnerBudgetValueV2::U64(1_000_000_000),
        RunnerBudgetValueV2::U64(1_000_000_000),
        RunnerBudgetValueV2::U64(1_000_000_000),
    ];

    fn durable_candidate() -> RunnerBudgetsCandidateV2 {
        RunnerBudgetsCandidateV2 {
            wall_time_ns: 100 * SECOND_NS,
            max_resident_bytes: GIB,
            max_child_processes: 8,
            max_parallel_children: 4,
            logical_work_limit: 1_000,
            logical_work_unit: LogicalUnitV2::Operations,
            lifecycle_encoded_bytes: 1000,
            command_result_stdout_bytes: 4000,
            combined_child_stdout_bytes: 2000,
            combined_child_stderr_bytes: 1000,
            artifact_encoded_bytes: 100,
            artifact_stored_bytes: 104,
            artifact_expanded_bytes: 200,
            system_publication_stored_bytes: 72,
            publication_stored_bytes: 176,
            stop_observation_ns: SECOND_NS,
            drain_ns: SECOND_NS,
            finalize_ns: SECOND_NS,
        }
    }

    fn lifecycle_candidate() -> RunnerBudgetsCandidateV2 {
        RunnerBudgetsCandidateV2 {
            artifact_encoded_bytes: 0,
            artifact_stored_bytes: 0,
            artifact_expanded_bytes: 0,
            system_publication_stored_bytes: 0,
            publication_stored_bytes: 0,
            max_child_processes: 0,
            max_parallel_children: 0,
            ..durable_candidate()
        }
    }

    #[allow(
        clippy::too_many_arguments,
        reason = "the helper compares every independent field of the frozen structured refusal tuple"
    )]
    fn assert_refusal_tuple(
        error: RunnerBudgetViolationV2,
        kind: RunnerBudgetViolationKindV2,
        field: RunnerBudgetFieldV2,
        unit: RunnerBudgetUnitV2,
        expected: RunnerBudgetExpectationV2,
        observed: RunnerBudgetValueV2,
        repair_kind: RepairActionKindV2,
        repair_target: &str,
    ) {
        assert_eq!(error.kind(), kind);
        assert_eq!(error.field(), field);
        assert_eq!(error.unit(), unit);
        assert_eq!(error.expected(), expected);
        assert_eq!(error.observed(), observed);
        assert_eq!(error.owner(), "fs-evidence-runner.runner-budgets");
        assert_eq!(error.repair_rank(), 1);
        assert_eq!(error.repair_kind(), repair_kind);
        assert_eq!(error.repair_target(), repair_target);
    }

    #[test]
    fn independent_literal_oracle_freezes_all_18_fields_and_widths() {
        assert_eq!(RunnerBudgetFieldV2::ALL.len(), 18);
        for index in 0..RUNNER_BUDGET_FIELD_COUNT_V2 {
            let field = RunnerBudgetFieldV2::ALL[index];
            let descriptor = RUNNER_BUDGET_DESCRIPTORS_V2[index];
            assert_eq!(field.ordinal(), u16::try_from(index + 1).unwrap());
            assert_eq!(
                RunnerBudgetFieldV2::from_ordinal(field.ordinal()),
                Some(field)
            );
            assert_eq!(descriptor.field, field);
            assert_eq!(descriptor.name, EXPECTED_NAMES[index]);
            assert_eq!(descriptor.width, EXPECTED_WIDTHS[index]);
            assert_eq!(descriptor.unit, EXPECTED_UNITS[index]);
        }
        assert_eq!(RunnerBudgetFieldV2::from_ordinal(0), None);
        assert_eq!(RunnerBudgetFieldV2::from_ordinal(19), None);

        let budgets = RunnerBudgetsV2::try_new(durable_candidate()).expect("literal fixture");
        for (field, expected) in RunnerBudgetFieldV2::ALL
            .into_iter()
            .zip(EXPECTED_DURABLE_VALUES)
        {
            assert_eq!(
                budgets.value(field),
                expected,
                "{}",
                field.descriptor().name
            );
        }
    }

    #[test]
    #[allow(
        clippy::too_many_lines,
        reason = "the test keeps every intrinsic nonzero relation and timeout-arithmetic boundary in one literal refusal oracle"
    )]
    fn intrinsic_nonzero_relations_and_timeout_arithmetic_refuse_precisely() {
        for (field, unit, repair_target) in [
            (
                RunnerBudgetFieldV2::WallTimeNs,
                RunnerBudgetUnitV2::Nanoseconds,
                "wall_time_ns",
            ),
            (
                RunnerBudgetFieldV2::MaxResidentBytes,
                RunnerBudgetUnitV2::LogicalBytes,
                "max_resident_bytes",
            ),
            (
                RunnerBudgetFieldV2::CommandResultStdoutBytes,
                RunnerBudgetUnitV2::EncodedBytes,
                "command_result_stdout_bytes",
            ),
            (
                RunnerBudgetFieldV2::StopObservationNs,
                RunnerBudgetUnitV2::Nanoseconds,
                "stop_observation_ns",
            ),
            (
                RunnerBudgetFieldV2::DrainNs,
                RunnerBudgetUnitV2::Nanoseconds,
                "drain_ns",
            ),
            (
                RunnerBudgetFieldV2::FinalizeNs,
                RunnerBudgetUnitV2::Nanoseconds,
                "finalize_ns",
            ),
        ] {
            let mut candidate = durable_candidate();
            match field {
                RunnerBudgetFieldV2::WallTimeNs => candidate.wall_time_ns = 0,
                RunnerBudgetFieldV2::MaxResidentBytes => candidate.max_resident_bytes = 0,
                RunnerBudgetFieldV2::CommandResultStdoutBytes => {
                    candidate.command_result_stdout_bytes = 0;
                }
                RunnerBudgetFieldV2::StopObservationNs => candidate.stop_observation_ns = 0,
                RunnerBudgetFieldV2::DrainNs => candidate.drain_ns = 0,
                RunnerBudgetFieldV2::FinalizeNs => candidate.finalize_ns = 0,
                _ => unreachable!(),
            }
            let error = RunnerBudgetsV2::try_new(candidate).unwrap_err();
            assert_refusal_tuple(
                error,
                RunnerBudgetViolationKindV2::Zero,
                field,
                unit,
                RunnerBudgetExpectationV2::NonZero,
                RunnerBudgetValueV2::U64(0),
                RepairActionKindV2::ChangeArguments,
                repair_target,
            );
        }

        let mut candidate = durable_candidate();
        candidate.max_parallel_children = 9;
        let error = RunnerBudgetsV2::try_new(candidate).unwrap_err();
        assert_refusal_tuple(
            error,
            RunnerBudgetViolationKindV2::ParallelChildrenExceedTotal,
            RunnerBudgetFieldV2::MaxParallelChildren,
            RunnerBudgetUnitV2::Count,
            RunnerBudgetExpectationV2::AtMost(RunnerBudgetValueV2::U32(8)),
            RunnerBudgetValueV2::U32(9),
            RepairActionKindV2::ReduceResourceDemand,
            "max_parallel_children",
        );

        let mut candidate = durable_candidate();
        candidate.stop_observation_ns = u64::MAX - 1;
        candidate.drain_ns = 1;
        candidate.finalize_ns = 1;
        let error = RunnerBudgetsV2::try_new(candidate).unwrap_err();
        assert_refusal_tuple(
            error,
            RunnerBudgetViolationKindV2::TimeoutSumOverflow,
            RunnerBudgetFieldV2::FinalizeNs,
            RunnerBudgetUnitV2::Nanoseconds,
            RunnerBudgetExpectationV2::AtMost(RunnerBudgetValueV2::U128(u128::from(
                candidate.wall_time_ns,
            ))),
            RunnerBudgetValueV2::U128(u128::from(u64::MAX) + 1),
            RepairActionKindV2::ReduceResourceDemand,
            "stop_observation_ns-or-drain_ns-or-finalize_ns",
        );

        let mut candidate = durable_candidate();
        candidate.wall_time_ns = 3 * SECOND_NS - 1;
        let error = RunnerBudgetsV2::try_new(candidate).unwrap_err();
        assert_refusal_tuple(
            error,
            RunnerBudgetViolationKindV2::TimeoutSumExceedsWall,
            RunnerBudgetFieldV2::FinalizeNs,
            RunnerBudgetUnitV2::Nanoseconds,
            RunnerBudgetExpectationV2::AtMost(RunnerBudgetValueV2::U64(3 * SECOND_NS - 1)),
            RunnerBudgetValueV2::U64(3 * SECOND_NS),
            RepairActionKindV2::ReduceResourceDemand,
            "stop_observation_ns-or-drain_ns-or-finalize_ns",
        );

        let mut exact = durable_candidate();
        exact.wall_time_ns = 3 * SECOND_NS;
        RunnerBudgetsV2::try_new(exact).expect("timeout sum equal to wall time");
    }

    #[test]
    #[allow(
        clippy::too_many_lines,
        reason = "the test exhaustively couples both profiles with all exact and one-over contextual ceilings"
    )]
    fn profile_boundaries_are_exact_and_one_over_refuses() {
        for (profile, wall_time_ns, resident_bytes, parallel_children, total_children) in [
            (RunProfileV2::Smoke, 900 * SECOND_NS, 16 * GIB, 32, 256),
            (RunProfileV2::Full, 86_400 * SECOND_NS, 128 * GIB, 64, 256),
        ] {
            let limits = RunnerLimitsV2::base(profile);
            let mut candidate = durable_candidate();
            candidate.wall_time_ns = wall_time_ns;
            candidate.max_resident_bytes = resident_bytes;
            candidate.max_parallel_children = parallel_children;
            candidate.max_child_processes = total_children;
            let exact = RunnerBudgetsV2::try_new(candidate)
                .expect("intrinsically valid exact profile vector")
                .admit(
                    profile,
                    ArtifactDispositionV2::DurableBundleRequired,
                    &limits,
                )
                .expect("all four exact profile ceilings");
            assert_eq!(exact.profile(), profile);
            assert_eq!(
                exact.disposition(),
                ArtifactDispositionV2::DurableBundleRequired
            );
            assert_eq!(exact.budgets().wall_time_ns(), wall_time_ns);
            assert_eq!(exact.budgets().max_resident_bytes(), resident_bytes);
            assert_eq!(exact.budgets().max_parallel_children(), parallel_children);
            assert_eq!(exact.budgets().max_child_processes(), total_children);

            for field in [
                RunnerBudgetFieldV2::WallTimeNs,
                RunnerBudgetFieldV2::MaxResidentBytes,
                RunnerBudgetFieldV2::MaxChildProcesses,
                RunnerBudgetFieldV2::MaxParallelChildren,
            ] {
                let mut candidate = durable_candidate();
                match field {
                    RunnerBudgetFieldV2::WallTimeNs => {
                        candidate.wall_time_ns = wall_time_ns + 1;
                    }
                    RunnerBudgetFieldV2::MaxResidentBytes => {
                        candidate.max_resident_bytes = resident_bytes + 1;
                    }
                    RunnerBudgetFieldV2::MaxChildProcesses => {
                        candidate.max_child_processes = total_children + 1;
                    }
                    RunnerBudgetFieldV2::MaxParallelChildren => {
                        candidate.max_child_processes = total_children;
                        candidate.max_parallel_children = parallel_children + 1;
                    }
                    _ => unreachable!(),
                }
                let error = RunnerBudgetsV2::try_new(candidate)
                    .expect("one-over profile vector remains intrinsically valid")
                    .admit(
                        profile,
                        ArtifactDispositionV2::DurableBundleRequired,
                        &limits,
                    )
                    .expect_err("one-over profile ceiling");
                assert_eq!(
                    error.kind(),
                    RunnerBudgetViolationKindV2::ProfileCeilingExceeded
                );
                assert_eq!(error.field(), field);
                assert_eq!(error.unit(), field.descriptor().unit);
                let (expected, observed) = match field {
                    RunnerBudgetFieldV2::WallTimeNs => (
                        RunnerBudgetExpectationV2::AtMost(RunnerBudgetValueV2::U64(wall_time_ns)),
                        RunnerBudgetValueV2::U64(wall_time_ns + 1),
                    ),
                    RunnerBudgetFieldV2::MaxResidentBytes => (
                        RunnerBudgetExpectationV2::AtMost(RunnerBudgetValueV2::U64(resident_bytes)),
                        RunnerBudgetValueV2::U64(resident_bytes + 1),
                    ),
                    RunnerBudgetFieldV2::MaxChildProcesses => (
                        RunnerBudgetExpectationV2::AtMost(RunnerBudgetValueV2::U32(total_children)),
                        RunnerBudgetValueV2::U32(total_children + 1),
                    ),
                    RunnerBudgetFieldV2::MaxParallelChildren => (
                        RunnerBudgetExpectationV2::AtMost(RunnerBudgetValueV2::U32(
                            parallel_children,
                        )),
                        RunnerBudgetValueV2::U32(parallel_children + 1),
                    ),
                    _ => unreachable!(),
                };
                assert_eq!(error.expected(), expected);
                assert_eq!(error.observed(), observed);
                assert_eq!(error.owner(), "fs-evidence-runner.runner-budgets");
            }
        }

        let shared = RunnerBudgetsV2::try_new(durable_candidate())
            .expect("one intrinsic vector valid for both profiles");
        let smoke = shared
            .admit(
                RunProfileV2::Smoke,
                ArtifactDispositionV2::DurableBundleRequired,
                &RunnerLimitsV2::base(RunProfileV2::Smoke),
            )
            .expect("shared vector is Smoke-admissible");
        let full = shared
            .admit(
                RunProfileV2::Full,
                ArtifactDispositionV2::DurableBundleRequired,
                &RunnerLimitsV2::base(RunProfileV2::Full),
            )
            .expect("shared vector is Full-admissible");
        assert_eq!(
            smoke.budgets().semantic_root(),
            full.budgets().semantic_root(),
            "intrinsic budget identity is context-free"
        );
        assert_ne!(
            smoke, full,
            "profile remains an identity-sensitive field of contextual admission"
        );
        assert_eq!(smoke.profile(), RunProfileV2::Smoke);
        assert_eq!(full.profile(), RunProfileV2::Full);
    }

    #[test]
    #[allow(
        clippy::too_many_lines,
        reason = "the test is the single exhaustive oracle for lifecycle-only zero rules and durable publication algebra"
    )]
    fn disposition_zero_rules_and_publication_equation_are_exact() {
        let limits = RunnerLimitsV2::base(RunProfileV2::Smoke);
        let admitted = RunnerBudgetsV2::try_new(lifecycle_candidate())
            .unwrap()
            .admit(
                RunProfileV2::Smoke,
                ArtifactDispositionV2::LifecycleOnlyNoBundle,
                &limits,
            )
            .unwrap();
        assert_eq!(
            admitted.disposition(),
            ArtifactDispositionV2::LifecycleOnlyNoBundle
        );

        for (field, unit, repair_target) in [
            (
                RunnerBudgetFieldV2::ArtifactEncodedBytes,
                RunnerBudgetUnitV2::EncodedBytes,
                "artifact_encoded_bytes",
            ),
            (
                RunnerBudgetFieldV2::ArtifactStoredBytes,
                RunnerBudgetUnitV2::StoredBytes,
                "artifact_stored_bytes",
            ),
            (
                RunnerBudgetFieldV2::ArtifactExpandedBytes,
                RunnerBudgetUnitV2::ExpandedBytes,
                "artifact_expanded_bytes",
            ),
            (
                RunnerBudgetFieldV2::SystemPublicationStoredBytes,
                RunnerBudgetUnitV2::StoredBytes,
                "system_publication_stored_bytes",
            ),
            (
                RunnerBudgetFieldV2::PublicationStoredBytes,
                RunnerBudgetUnitV2::StoredBytes,
                "publication_stored_bytes",
            ),
        ] {
            let mut candidate = lifecycle_candidate();
            match field {
                RunnerBudgetFieldV2::ArtifactEncodedBytes => {
                    candidate.artifact_encoded_bytes = 1;
                    candidate.artifact_stored_bytes = 1;
                    candidate.publication_stored_bytes = 1;
                }
                RunnerBudgetFieldV2::ArtifactStoredBytes => {
                    candidate.artifact_stored_bytes = 1;
                    candidate.publication_stored_bytes = 1;
                }
                RunnerBudgetFieldV2::ArtifactExpandedBytes => {
                    candidate.artifact_expanded_bytes = 1;
                }
                RunnerBudgetFieldV2::SystemPublicationStoredBytes => {
                    candidate.system_publication_stored_bytes = 1;
                    candidate.publication_stored_bytes = 1;
                }
                RunnerBudgetFieldV2::PublicationStoredBytes => {
                    candidate.publication_stored_bytes = 1;
                }
                _ => unreachable!(),
            }
            let error = RunnerBudgetsV2::try_new(candidate)
                .unwrap()
                .admit(
                    RunProfileV2::Smoke,
                    ArtifactDispositionV2::LifecycleOnlyNoBundle,
                    &limits,
                )
                .unwrap_err();
            assert_refusal_tuple(
                error,
                RunnerBudgetViolationKindV2::ContextualZeroRequired,
                field,
                unit,
                RunnerBudgetExpectationV2::Zero,
                RunnerBudgetValueV2::U64(1),
                RepairActionKindV2::ChangeArguments,
                repair_target,
            );
        }

        let mut candidate = durable_candidate();
        candidate.command_result_stdout_bytes = candidate.lifecycle_encoded_bytes - 1;
        let error = RunnerBudgetsV2::try_new(candidate)
            .unwrap()
            .admit(
                RunProfileV2::Smoke,
                ArtifactDispositionV2::DurableBundleRequired,
                &limits,
            )
            .unwrap_err();
        assert_refusal_tuple(
            error,
            RunnerBudgetViolationKindV2::CommandResultCannotContainLifecycle,
            RunnerBudgetFieldV2::CommandResultStdoutBytes,
            RunnerBudgetUnitV2::EncodedBytes,
            RunnerBudgetExpectationV2::AtLeast(RunnerBudgetValueV2::U64(1_000)),
            RunnerBudgetValueV2::U64(999),
            RepairActionKindV2::ChangeArguments,
            "command_result_stdout_bytes",
        );

        let mut candidate = durable_candidate();
        candidate.artifact_stored_bytes = candidate.artifact_encoded_bytes - 1;
        let error = RunnerBudgetsV2::try_new(candidate)
            .unwrap()
            .admit(
                RunProfileV2::Smoke,
                ArtifactDispositionV2::DurableBundleRequired,
                &limits,
            )
            .unwrap_err();
        assert_refusal_tuple(
            error,
            RunnerBudgetViolationKindV2::ArtifactStoredBelowEncoded,
            RunnerBudgetFieldV2::ArtifactStoredBytes,
            RunnerBudgetUnitV2::StoredBytes,
            RunnerBudgetExpectationV2::AtLeast(RunnerBudgetValueV2::U64(100)),
            RunnerBudgetValueV2::U64(99),
            RepairActionKindV2::ChangeArguments,
            "artifact_stored_bytes",
        );

        let mut candidate = durable_candidate();
        candidate.publication_stored_bytes -= 1;
        let error = RunnerBudgetsV2::try_new(candidate)
            .unwrap()
            .admit(
                RunProfileV2::Smoke,
                ArtifactDispositionV2::DurableBundleRequired,
                &limits,
            )
            .unwrap_err();
        assert_refusal_tuple(
            error,
            RunnerBudgetViolationKindV2::PublicationEquationMismatch,
            RunnerBudgetFieldV2::PublicationStoredBytes,
            RunnerBudgetUnitV2::StoredBytes,
            RunnerBudgetExpectationV2::Exactly(RunnerBudgetValueV2::U64(176)),
            RunnerBudgetValueV2::U64(175),
            RepairActionKindV2::ChangeArguments,
            "publication_stored_bytes",
        );

        let mut candidate = durable_candidate();
        candidate.max_child_processes = 0;
        candidate.max_parallel_children = 0;
        let error = RunnerBudgetsV2::try_new(candidate)
            .unwrap()
            .admit(
                RunProfileV2::Smoke,
                ArtifactDispositionV2::DurableBundleRequired,
                &limits,
            )
            .unwrap_err();
        assert_refusal_tuple(
            error,
            RunnerBudgetViolationKindV2::ContextualNonZeroRequired,
            RunnerBudgetFieldV2::MaxChildProcesses,
            RunnerBudgetUnitV2::Count,
            RunnerBudgetExpectationV2::NonZero,
            RunnerBudgetValueV2::U32(0),
            RepairActionKindV2::ChangeArguments,
            "max_child_processes",
        );

        let mut candidate = durable_candidate();
        candidate.system_publication_stored_bytes = 0;
        candidate.publication_stored_bytes = candidate.artifact_stored_bytes;
        let error = RunnerBudgetsV2::try_new(candidate)
            .unwrap()
            .admit(
                RunProfileV2::Smoke,
                ArtifactDispositionV2::DurableBundleRequired,
                &limits,
            )
            .unwrap_err();
        assert_refusal_tuple(
            error,
            RunnerBudgetViolationKindV2::ContextualNonZeroRequired,
            RunnerBudgetFieldV2::SystemPublicationStoredBytes,
            RunnerBudgetUnitV2::StoredBytes,
            RunnerBudgetExpectationV2::NonZero,
            RunnerBudgetValueV2::U64(0),
            RepairActionKindV2::ChangeArguments,
            "system_publication_stored_bytes",
        );
    }

    #[test]
    #[allow(
        clippy::too_many_lines,
        reason = "the test intentionally enumerates every independently admitted output grant and its exact limit refusal tuple"
    )]
    fn every_output_grant_is_checked_against_the_admitted_limits() {
        let limits = RunnerLimitsV2::base(RunProfileV2::Smoke);
        for field in [
            RunnerBudgetFieldV2::LifecycleEncodedBytes,
            RunnerBudgetFieldV2::CommandResultStdoutBytes,
            RunnerBudgetFieldV2::CombinedChildStdoutBytes,
            RunnerBudgetFieldV2::CombinedChildStderrBytes,
            RunnerBudgetFieldV2::ArtifactEncodedBytes,
            RunnerBudgetFieldV2::ArtifactStoredBytes,
            RunnerBudgetFieldV2::ArtifactExpandedBytes,
            RunnerBudgetFieldV2::SystemPublicationStoredBytes,
            RunnerBudgetFieldV2::PublicationStoredBytes,
        ] {
            let mut candidate = durable_candidate();
            let (ceiling, unit, repair_target) = match field {
                RunnerBudgetFieldV2::LifecycleEncodedBytes => {
                    let ceiling = limits.lifecycle_document_encoded_bytes();
                    candidate.lifecycle_encoded_bytes = ceiling + 1;
                    candidate.command_result_stdout_bytes = candidate.lifecycle_encoded_bytes;
                    (
                        ceiling,
                        RunnerBudgetUnitV2::EncodedBytes,
                        "lifecycle_encoded_bytes",
                    )
                }
                RunnerBudgetFieldV2::CommandResultStdoutBytes => {
                    let ceiling = limits.command_result_stdout_bytes();
                    candidate.command_result_stdout_bytes = ceiling + 1;
                    (
                        ceiling,
                        RunnerBudgetUnitV2::EncodedBytes,
                        "command_result_stdout_bytes",
                    )
                }
                RunnerBudgetFieldV2::CombinedChildStdoutBytes => {
                    let ceiling = limits.combined_child_stdout_bytes();
                    candidate.combined_child_stdout_bytes = ceiling + 1;
                    (
                        ceiling,
                        RunnerBudgetUnitV2::EncodedBytes,
                        "combined_child_stdout_bytes",
                    )
                }
                RunnerBudgetFieldV2::CombinedChildStderrBytes => {
                    let ceiling = limits.combined_child_stderr_bytes();
                    candidate.combined_child_stderr_bytes = ceiling + 1;
                    (
                        ceiling,
                        RunnerBudgetUnitV2::EncodedBytes,
                        "combined_child_stderr_bytes",
                    )
                }
                RunnerBudgetFieldV2::ArtifactEncodedBytes => {
                    let ceiling = limits.bundle_encoded_bytes();
                    candidate.artifact_encoded_bytes = ceiling + 1;
                    candidate.artifact_stored_bytes = candidate.artifact_encoded_bytes;
                    candidate.publication_stored_bytes =
                        candidate.artifact_stored_bytes + candidate.system_publication_stored_bytes;
                    (
                        ceiling,
                        RunnerBudgetUnitV2::EncodedBytes,
                        "artifact_encoded_bytes",
                    )
                }
                RunnerBudgetFieldV2::ArtifactStoredBytes => {
                    let ceiling = limits.artifact_stored_aggregate_bytes();
                    candidate.artifact_stored_bytes = ceiling + 1;
                    candidate.publication_stored_bytes =
                        candidate.artifact_stored_bytes + candidate.system_publication_stored_bytes;
                    (
                        ceiling,
                        RunnerBudgetUnitV2::StoredBytes,
                        "artifact_stored_bytes",
                    )
                }
                RunnerBudgetFieldV2::ArtifactExpandedBytes => {
                    let ceiling = limits.bundle_expanded_bytes();
                    candidate.artifact_expanded_bytes = ceiling + 1;
                    (
                        ceiling,
                        RunnerBudgetUnitV2::ExpandedBytes,
                        "artifact_expanded_bytes",
                    )
                }
                RunnerBudgetFieldV2::SystemPublicationStoredBytes => {
                    let ceiling = limits.system_publication_stored_bytes();
                    candidate.system_publication_stored_bytes = ceiling + 1;
                    candidate.publication_stored_bytes =
                        candidate.artifact_stored_bytes + candidate.system_publication_stored_bytes;
                    (
                        ceiling,
                        RunnerBudgetUnitV2::StoredBytes,
                        "system_publication_stored_bytes",
                    )
                }
                RunnerBudgetFieldV2::PublicationStoredBytes => {
                    let ceiling = limits.publication_stored_bytes();
                    candidate.publication_stored_bytes = ceiling + 1;
                    (
                        ceiling,
                        RunnerBudgetUnitV2::StoredBytes,
                        "publication_stored_bytes",
                    )
                }
                _ => unreachable!(),
            };
            let error = RunnerBudgetsV2::try_new(candidate)
                .unwrap()
                .admit(
                    RunProfileV2::Smoke,
                    ArtifactDispositionV2::DurableBundleRequired,
                    &limits,
                )
                .unwrap_err();
            let expected_ceiling = match field {
                RunnerBudgetFieldV2::LifecycleEncodedBytes => 4 * TEST_MIB,
                RunnerBudgetFieldV2::CommandResultStdoutBytes => 5 * TEST_MIB,
                RunnerBudgetFieldV2::CombinedChildStdoutBytes => 16 * TEST_MIB,
                RunnerBudgetFieldV2::CombinedChildStderrBytes => TEST_MIB / 4,
                RunnerBudgetFieldV2::ArtifactEncodedBytes
                | RunnerBudgetFieldV2::ArtifactExpandedBytes => 64 * TEST_MIB,
                RunnerBudgetFieldV2::ArtifactStoredBytes => 65 * TEST_MIB,
                RunnerBudgetFieldV2::SystemPublicationStoredBytes => 8 * TEST_MIB,
                RunnerBudgetFieldV2::PublicationStoredBytes => 73 * TEST_MIB,
                _ => unreachable!(),
            };
            assert_eq!(
                ceiling, expected_ceiling,
                "{repair_target} Smoke ceiling drifted from the literal refusal oracle",
            );
            assert_refusal_tuple(
                error,
                RunnerBudgetViolationKindV2::LimitExceeded,
                field,
                unit,
                RunnerBudgetExpectationV2::AtMost(RunnerBudgetValueV2::U64(expected_ceiling)),
                RunnerBudgetValueV2::U64(expected_ceiling + 1),
                RepairActionKindV2::ReduceResourceDemand,
                repair_target,
            );
        }
    }

    #[test]
    fn logical_work_is_exact_u128_and_registered_units_require_exact_registry_membership() {
        let mut candidate = durable_candidate();
        candidate.logical_work_limit = u128::MAX;
        candidate.logical_work_unit = LogicalUnitV2::from_tag(16, Some(65_535)).unwrap();
        let error = RunnerBudgetsV2::try_new(candidate).expect_err("bare registered unit");
        assert_refusal_tuple(
            error,
            RunnerBudgetViolationKindV2::UnregisteredLogicalWorkUnit,
            RunnerBudgetFieldV2::LogicalWorkUnit,
            RunnerBudgetUnitV2::LogicalWorkUnit,
            RunnerBudgetExpectationV2::RegisteredInExtensionRegistry,
            RunnerBudgetValueV2::LogicalUnit {
                tag: 16,
                registered_id: Some(65_535),
            },
            RepairActionKindV2::UpdatePolicyOrCapability,
            "logical_work_unit",
        );

        let no_claim = NoClaimScopeRootV1::parse_presented(
            DigestRoleV2::ClaimScope,
            NoClaimScopeRootV1::DESCRIPTOR.domain(),
            &"55".repeat(32),
        )
        .unwrap();
        let name = StableTokenV2::new("org.example.unit.work").unwrap();
        let owner = StableTokenV2::new("org.example.owner").unwrap();
        let same_numeric_wrong_category = RegisteredArtifactRoleDescriptorV2::new(
            ArtifactRoleV2::from_tag(8, Some(65_535)).unwrap(),
            StableTokenV2::new("org.example.role.work").unwrap(),
            owner.clone(),
            no_claim.clone(),
        )
        .unwrap();
        let limits = RunnerLimitsV2::base(RunProfileV2::Smoke);
        let wrong_category = BaseExtensionRegistryProjectionV2::try_new(
            &limits,
            &[same_numeric_wrong_category],
            &[],
            &[],
        )
        .unwrap();
        assert_eq!(
            RunnerBudgetsV2::try_new_with_extension_registry(candidate, &wrong_category)
                .expect_err("same numeric ID in a role namespace")
                .kind(),
            RunnerBudgetViolationKindV2::UnregisteredLogicalWorkUnit
        );

        let descriptor = RegisteredLogicalUnitDescriptorV2::new(
            candidate.logical_work_unit,
            name,
            owner,
            no_claim,
        )
        .unwrap();
        let registry =
            BaseExtensionRegistryProjectionV2::try_new(&limits, &[], &[descriptor], &[]).unwrap();
        let bound = RunnerBudgetsV2::try_new_with_extension_registry(candidate, &registry).unwrap();
        assert_eq!(bound.budgets().logical_work_limit(), u128::MAX);
        assert_eq!(
            bound.budgets().logical_work_unit(),
            LogicalUnitV2::from_tag(16, Some(65_535)).unwrap()
        );
        assert_eq!(bound.extension_registry_root(), *registry.root());
        let projection = bound.budgets().canonical_projection();
        assert!(
            projection
                .windows(4)
                .any(|window| window == [0, 16, 255, 255])
        );

        let alternate_descriptor = RegisteredLogicalUnitDescriptorV2::new(
            candidate.logical_work_unit,
            StableTokenV2::new("org.example.unit.work").unwrap(),
            StableTokenV2::new("org.example.alternate-owner").unwrap(),
            NoClaimScopeRootV1::parse_presented(
                DigestRoleV2::ClaimScope,
                NoClaimScopeRootV1::DESCRIPTOR.domain(),
                &"55".repeat(32),
            )
            .unwrap(),
        )
        .unwrap();
        let alternate_registry =
            BaseExtensionRegistryProjectionV2::try_new(&limits, &[], &[alternate_descriptor], &[])
                .unwrap();
        let alternate_bound =
            RunnerBudgetsV2::try_new_with_extension_registry(candidate, &alternate_registry)
                .unwrap();
        assert_eq!(
            bound.budgets().semantic_root(),
            alternate_bound.budgets().semantic_root(),
            "the frozen 18-field budget identity remains unchanged"
        );
        assert_ne!(
            bound.semantic_binding_root(),
            alternate_bound.semantic_binding_root(),
            "registry descriptor drift must move the enclosing membership binding"
        );
        let admitted = bound
            .admit(
                RunProfileV2::Smoke,
                ArtifactDispositionV2::DurableBundleRequired,
                &limits,
            )
            .expect("registry-bound budget admission");
        assert_eq!(admitted.extension_registry_root(), *registry.root());
    }

    #[test]
    fn every_one_field_mutation_moves_the_canonical_projection_and_root() {
        let base = RunnerBudgetsV2::try_new(durable_candidate()).unwrap();
        let base_projection = base.canonical_projection();
        let base_root = base.canonical_projection_root();
        let base_semantic_root = base.semantic_root();
        assert_eq!(base_semantic_root.role(), DigestRoleV2::Policy);
        assert_eq!(
            base_semantic_root.domain(),
            RunnerBudgetsRootV2::DESCRIPTOR.domain()
        );
        assert_eq!(base_semantic_root.bytes(), base_root.as_bytes());
        for field in RunnerBudgetFieldV2::ALL {
            let mut candidate = durable_candidate();
            match field {
                RunnerBudgetFieldV2::WallTimeNs => candidate.wall_time_ns += 1,
                RunnerBudgetFieldV2::MaxResidentBytes => candidate.max_resident_bytes += 1,
                RunnerBudgetFieldV2::MaxChildProcesses => candidate.max_child_processes += 1,
                RunnerBudgetFieldV2::MaxParallelChildren => {
                    candidate.max_parallel_children += 1;
                }
                RunnerBudgetFieldV2::LogicalWorkLimit => candidate.logical_work_limit += 1,
                RunnerBudgetFieldV2::LogicalWorkUnit => {
                    candidate.logical_work_unit = LogicalUnitV2::Cycles;
                }
                RunnerBudgetFieldV2::LifecycleEncodedBytes => {
                    candidate.lifecycle_encoded_bytes += 1;
                }
                RunnerBudgetFieldV2::CommandResultStdoutBytes => {
                    candidate.command_result_stdout_bytes += 1;
                }
                RunnerBudgetFieldV2::CombinedChildStdoutBytes => {
                    candidate.combined_child_stdout_bytes += 1;
                }
                RunnerBudgetFieldV2::CombinedChildStderrBytes => {
                    candidate.combined_child_stderr_bytes += 1;
                }
                RunnerBudgetFieldV2::ArtifactEncodedBytes => {
                    candidate.artifact_encoded_bytes += 1;
                }
                RunnerBudgetFieldV2::ArtifactStoredBytes => {
                    candidate.artifact_stored_bytes += 1;
                }
                RunnerBudgetFieldV2::ArtifactExpandedBytes => {
                    candidate.artifact_expanded_bytes += 1;
                }
                RunnerBudgetFieldV2::SystemPublicationStoredBytes => {
                    candidate.system_publication_stored_bytes += 1;
                }
                RunnerBudgetFieldV2::PublicationStoredBytes => {
                    candidate.publication_stored_bytes += 1;
                }
                RunnerBudgetFieldV2::StopObservationNs => {
                    candidate.stop_observation_ns += 1;
                }
                RunnerBudgetFieldV2::DrainNs => candidate.drain_ns += 1,
                RunnerBudgetFieldV2::FinalizeNs => candidate.finalize_ns += 1,
            }
            let mutated = RunnerBudgetsV2::try_new(candidate).unwrap();
            let changed_fields = RunnerBudgetFieldV2::ALL
                .into_iter()
                .filter(|candidate_field| {
                    mutated.value(*candidate_field) != base.value(*candidate_field)
                })
                .collect::<Vec<_>>();
            assert_eq!(
                changed_fields,
                vec![field],
                "{} mutation must change exactly its own typed budget field",
                field.descriptor().name
            );
            assert_ne!(
                mutated.canonical_projection(),
                base_projection,
                "{}",
                field.descriptor().name
            );
            assert_ne!(
                mutated.canonical_projection_root(),
                base_root,
                "{}",
                field.descriptor().name
            );
            assert_ne!(
                mutated.semantic_root().bytes(),
                base_semantic_root.bytes(),
                "{}",
                field.descriptor().name
            );
        }
    }

    #[test]
    fn precise_refusal_carries_unit_expectation_observation_and_owner() {
        let mut candidate = durable_candidate();
        candidate.max_parallel_children = candidate.max_child_processes + 1;
        let error = RunnerBudgetsV2::try_new(candidate).unwrap_err();
        assert_eq!(error.field(), RunnerBudgetFieldV2::MaxParallelChildren);
        assert_eq!(error.unit(), RunnerBudgetUnitV2::Count);
        assert_eq!(
            error.expected(),
            RunnerBudgetExpectationV2::AtMost(RunnerBudgetValueV2::U32(8))
        );
        assert_eq!(error.observed(), RunnerBudgetValueV2::U32(9));
        assert_eq!(error.owner(), "fs-evidence-runner.runner-budgets");
        assert_eq!(error.repair_rank(), 1);
        assert_eq!(
            error.repair_kind(),
            RepairActionKindV2::ReduceResourceDemand
        );
        assert_eq!(error.repair_target(), "max_parallel_children");
    }

    #[test]
    fn publication_sum_overflow_is_typed_even_without_concrete_storage() {
        let limits = RunnerLimitsV2::base(RunProfileV2::Full);
        let mut candidate = durable_candidate();
        candidate.artifact_encoded_bytes = 1;
        candidate.artifact_stored_bytes = u64::MAX;
        candidate.system_publication_stored_bytes = 1;
        candidate.publication_stored_bytes = 0;
        let budgets = RunnerBudgetsV2::try_new(candidate).unwrap();
        // The normal contextual limit check is intentionally earlier and
        // rejects this unbounded grant before the sum is evaluated.
        let error = budgets
            .admit(
                RunProfileV2::Full,
                ArtifactDispositionV2::DurableBundleRequired,
                &limits,
            )
            .unwrap_err();
        assert_eq!(error.kind(), RunnerBudgetViolationKindV2::LimitExceeded);

        let error = checked_publication_stored_bytes(u64::MAX, 1).unwrap_err();
        assert_eq!(
            error.kind(),
            RunnerBudgetViolationKindV2::PublicationSumOverflow
        );
        assert_eq!(error.field(), RunnerBudgetFieldV2::PublicationStoredBytes);
        assert_eq!(error.unit(), RunnerBudgetUnitV2::StoredBytes);
        assert_eq!(
            error.expected(),
            RunnerBudgetExpectationV2::AtMost(RunnerBudgetValueV2::U128(u128::from(u64::MAX,)))
        );
        assert_eq!(
            error.observed(),
            RunnerBudgetValueV2::U128(u128::from(u64::MAX) + 1)
        );
        assert_eq!(error.owner(), "fs-evidence-runner.runner-budgets");
        assert_eq!(error.repair_rank(), 1);
        assert_eq!(
            error.repair_kind(),
            RepairActionKindV2::ReduceResourceDemand
        );
        assert_eq!(
            error.repair_target(),
            "artifact_stored_bytes-or-system_publication_stored_bytes"
        );
    }

    #[test]
    fn publication_accounting_accepts_zero_one_exact_and_maximum_boundaries() {
        assert_eq!(checked_publication_stored_bytes(0, 0).unwrap(), 0);
        assert_eq!(checked_publication_stored_bytes(1, 0).unwrap(), 1);
        assert_eq!(checked_publication_stored_bytes(0, 1).unwrap(), 1);
        assert_eq!(
            checked_publication_stored_bytes(u64::MAX - 1, 1).unwrap(),
            u64::MAX
        );
        assert_eq!(
            checked_publication_stored_bytes(u64::MAX, 0).unwrap(),
            u64::MAX
        );
        let error = checked_publication_stored_bytes(u64::MAX, 1).unwrap_err();
        assert_eq!(
            error.kind(),
            RunnerBudgetViolationKindV2::PublicationSumOverflow
        );
    }
}
