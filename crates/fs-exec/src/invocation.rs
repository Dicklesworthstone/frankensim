//! Affine, invocation-scoped resource accounting.
//!
//! `asupersync::Budget` is a copyable propagation envelope.  It is not a
//! spend ledger.  This module supplies the separate, non-`Clone` authority a
//! composed scientific invocation needs when sibling calls must not recreate
//! the ambient poll or cost allowance.

use crate::{CancelGate, Cx};
use fs_alloc::{LeaseCharge, OperationMemoryLease};
use fs_blake3::{ContentHash, DomainHasher};
#[cfg(test)]
use fs_blake3::hash_domain;

pub use asupersync::time::{TimeSource, VirtualClock, WallClock};
pub use asupersync::types::Time;

/// Version of the canonical invocation-accounting receipt.
///
/// Version 2 binds the exact finalizer partition and report commitment for
/// finalizable children. Version-1 receipts remain historical evidence; this
/// producer/verifier deliberately does not silently reinterpret them.
pub const INVOCATION_RECEIPT_VERSION: u32 = 2;
/// Maximum child rows accepted by the producer and semantic verifier.
pub const INVOCATION_RECEIPT_MAX_CHILDREN: u64 = 1_048_576;

const CHILD_ID_DOMAIN: &str = "frankensim.fs-exec.invocation-child.v2";
const CHILD_RECEIPT_DOMAIN: &str = "frankensim.fs-exec.invocation-child-receipt.v2";
const INVOCATION_RECEIPT_DOMAIN: &str = "frankensim.fs-exec.invocation-receipt.v2";
const FINALIZATION_REPORT_DOMAIN: &str = "frankensim.fs-exec.finalization-report.v2";
const FINALIZED_CHILD_RECEIPT_DOMAIN: &str = "frankensim.fs-exec.finalized-child-receipt.v1";
const INVOCATION_ERROR_DOMAIN: &str = "frankensim.fs-exec.invocation-error.v1";

/// Version of the post-cancel cleanup/finalization report.
pub const FINALIZATION_REPORT_VERSION: u32 = 2;

macro_rules! resource_unit {
    ($name:ident, $repr:ty, $docs:literal) => {
        #[doc = $docs]
        #[derive(Debug, Clone, Copy, Default, PartialEq, Eq, PartialOrd, Ord)]
        pub struct $name($repr);

        impl $name {
            /// Construct a typed quantity.
            #[must_use]
            pub const fn new(value: $repr) -> Self {
                Self(value)
            }

            /// Raw diagnostic value in this type's declared unit.
            #[must_use]
            pub const fn get(self) -> $repr {
                self.0
            }
        }
    };
}

resource_unit!(WorkUnits, u128, "Declared logical work units.");
resource_unit!(PollUnits, u32, "Cancellation/deadline poll opportunities.");
resource_unit!(CostUnits, u64, "Abstract monetary or energy cost units.");
resource_unit!(EvaluationUnits, u64, "Scientific evaluation count.");
resource_unit!(MemoryBytes, u64, "Concurrent live memory bytes.");
resource_unit!(OutputBytes, u64, "Retained publication capacity in bytes.");

/// Resources isolated from ordinary scientific work and reserved for terminal
/// cleanup. Finalization is deliberately limited to logical work and polls:
/// it may release already-admitted storage and seal or abort publication, but
/// it cannot allocate new scientific memory, mint evaluations, or spend output
/// capacity as compute.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct FinalizationResources {
    work: WorkUnits,
    polls: PollUnits,
}

impl FinalizationResources {
    /// Construct a typed cleanup reserve.
    #[must_use]
    pub const fn new(work: WorkUnits, polls: PollUnits) -> Self {
        Self { work, polls }
    }

    /// Reserved cleanup work.
    #[must_use]
    pub const fn work(self) -> WorkUnits {
        self.work
    }

    /// Reserved cleanup polls.
    #[must_use]
    pub const fn polls(self) -> PollUnits {
        self.polls
    }

    /// Resource-vector form used when a caller builds the complete root
    /// preflight. Only work and polls are nonzero.
    #[must_use]
    pub const fn as_invocation_resources(self) -> InvocationResources {
        InvocationResources::new(
            self.work,
            self.polls,
            CostUnits::new(0),
            EvaluationUnits::new(0),
            MemoryBytes::new(0),
            OutputBytes::new(0),
        )
    }

    fn checked_sub(self, used: Self) -> Result<Self, InvocationError> {
        Ok(Self {
            work: WorkUnits::new(self.work.get().checked_sub(used.work.get()).ok_or(
                InvocationError::ResourceExceeded {
                    resource: "finalization-work",
                    requested: used.work.get(),
                    available: self.work.get(),
                },
            )?),
            polls: PollUnits::new(self.polls.get().checked_sub(used.polls.get()).ok_or(
                InvocationError::ResourceExceeded {
                    resource: "finalization-polls",
                    requested: u128::from(used.polls.get()),
                    available: u128::from(self.polls.get()),
                },
            )?),
        })
    }
}

/// Dimensioned affine capacities.  Deliberately no generic numeric indexing or
/// cross-kind conversion exists.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct InvocationResources {
    work: WorkUnits,
    polls: PollUnits,
    cost: CostUnits,
    evaluations: EvaluationUnits,
    memory: MemoryBytes,
    output: OutputBytes,
}

impl InvocationResources {
    /// Construct one dimensionally explicit resource vector.
    #[must_use]
    pub const fn new(
        work: WorkUnits,
        polls: PollUnits,
        cost: CostUnits,
        evaluations: EvaluationUnits,
        memory: MemoryBytes,
        output: OutputBytes,
    ) -> Self {
        Self {
            work,
            polls,
            cost,
            evaluations,
            memory,
            output,
        }
    }

    /// Declared logical work.
    #[must_use]
    pub const fn work(self) -> WorkUnits {
        self.work
    }

    /// Poll opportunities.
    #[must_use]
    pub const fn polls(self) -> PollUnits {
        self.polls
    }

    /// Cost allowance.
    #[must_use]
    pub const fn cost(self) -> CostUnits {
        self.cost
    }

    /// Evaluation allowance.
    #[must_use]
    pub const fn evaluations(self) -> EvaluationUnits {
        self.evaluations
    }

    /// Concurrent-memory ceiling.
    #[must_use]
    pub const fn memory(self) -> MemoryBytes {
        self.memory
    }

    /// Retained-output capacity.
    #[must_use]
    pub const fn output(self) -> OutputBytes {
        self.output
    }

    /// Dimension-preserving checked subtraction.
    ///
    /// # Errors
    /// Refuses the first insufficient dimension in canonical resource order.
    pub fn checked_sub(self, requested: Self) -> Result<Self, InvocationError> {
        Ok(Self {
            work: WorkUnits(
                self.work
                    .0
                    .checked_sub(requested.work.0)
                    .ok_or_else(|| exceeded("work", requested.work.0, self.work.0))?,
            ),
            polls: PollUnits(self.polls.0.checked_sub(requested.polls.0).ok_or_else(|| {
                exceeded(
                    "polls",
                    u128::from(requested.polls.0),
                    u128::from(self.polls.0),
                )
            })?),
            cost: CostUnits(self.cost.0.checked_sub(requested.cost.0).ok_or_else(|| {
                exceeded(
                    "cost",
                    u128::from(requested.cost.0),
                    u128::from(self.cost.0),
                )
            })?),
            evaluations: EvaluationUnits(
                self.evaluations
                    .0
                    .checked_sub(requested.evaluations.0)
                    .ok_or_else(|| {
                        exceeded(
                            "evaluations",
                            u128::from(requested.evaluations.0),
                            u128::from(self.evaluations.0),
                        )
                    })?,
            ),
            memory: MemoryBytes(self.memory.0.checked_sub(requested.memory.0).ok_or_else(
                || {
                    exceeded(
                        "memory-bytes",
                        u128::from(requested.memory.0),
                        u128::from(self.memory.0),
                    )
                },
            )?),
            output: OutputBytes(self.output.0.checked_sub(requested.output.0).ok_or_else(
                || {
                    exceeded(
                        "output-bytes",
                        u128::from(requested.output.0),
                        u128::from(self.output.0),
                    )
                },
            )?),
        })
    }

    /// Dimension-preserving checked addition.
    ///
    /// # Errors
    /// Refuses representational overflow without changing either operand.
    pub fn checked_add(self, returned: Self) -> Result<Self, InvocationError> {
        Ok(Self {
            work: WorkUnits(
                self.work
                    .0
                    .checked_add(returned.work.0)
                    .ok_or(InvocationError::ArithmeticOverflow { resource: "work" })?,
            ),
            polls: PollUnits(
                self.polls
                    .0
                    .checked_add(returned.polls.0)
                    .ok_or(InvocationError::ArithmeticOverflow { resource: "polls" })?,
            ),
            cost: CostUnits(
                self.cost
                    .0
                    .checked_add(returned.cost.0)
                    .ok_or(InvocationError::ArithmeticOverflow { resource: "cost" })?,
            ),
            evaluations: EvaluationUnits(
                self.evaluations
                    .0
                    .checked_add(returned.evaluations.0)
                    .ok_or(InvocationError::ArithmeticOverflow {
                        resource: "evaluations",
                    })?,
            ),
            memory: MemoryBytes(self.memory.0.checked_add(returned.memory.0).ok_or(
                InvocationError::ArithmeticOverflow {
                    resource: "memory-bytes",
                },
            )?),
            output: OutputBytes(self.output.0.checked_add(returned.output.0).ok_or(
                InvocationError::ArithmeticOverflow {
                    resource: "output-bytes",
                },
            )?),
        })
    }
}

fn exceeded(resource: &'static str, requested: u128, available: u128) -> InvocationError {
    InvocationError::ResourceExceeded {
        resource,
        requested,
        available,
    }
}

/// Fixed admission envelope. Accuracy and capability are immutable identities,
/// not spendable counters.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InvocationLimits {
    resources: InvocationResources,
    deadline: Option<Time>,
    accuracy_obligation: ContentHash,
    capability_scope: ContentHash,
}

impl InvocationLimits {
    /// Construct a complete invocation envelope.
    #[must_use]
    pub const fn new(
        resources: InvocationResources,
        deadline: Option<Time>,
        accuracy_obligation: ContentHash,
        capability_scope: ContentHash,
    ) -> Self {
        Self {
            resources,
            deadline,
            accuracy_obligation,
            capability_scope,
        }
    }

    /// Affine capacity dimensions.
    #[must_use]
    pub const fn resources(&self) -> InvocationResources {
        self.resources
    }

    /// Absolute logical deadline.
    #[must_use]
    pub const fn deadline(&self) -> Option<Time> {
        self.deadline
    }

    /// Immutable accuracy/tolerance obligation identity.
    #[must_use]
    pub const fn accuracy_obligation(&self) -> ContentHash {
        self.accuracy_obligation
    }

    /// Immutable capability-authority scope identity.
    #[must_use]
    pub const fn capability_scope(&self) -> ContentHash {
        self.capability_scope
    }
}

/// Explicit semantic plan identity bound into strong invocation evidence.
///
/// `schema_root` identifies the plan grammar/encoder, `schema_version` selects
/// its exact admitted version, and `plan_root` identifies the concrete checked
/// operation plan. The binding is data only; it does not claim that the plan's
/// scientific mathematics is correct.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct InvocationPlanBinding {
    schema_root: ContentHash,
    schema_version: u32,
    plan_root: ContentHash,
}

impl InvocationPlanBinding {
    /// Construct one explicit checked-plan binding.
    #[must_use]
    pub const fn new(
        schema_root: ContentHash,
        schema_version: u32,
        plan_root: ContentHash,
    ) -> Self {
        Self {
            schema_root,
            schema_version,
            plan_root,
        }
    }

    /// Identity of the plan schema/encoder.
    #[must_use]
    pub const fn schema_root(self) -> ContentHash {
        self.schema_root
    }

    /// Exact plan-schema version.
    #[must_use]
    pub const fn schema_version(self) -> u32 {
        self.schema_version
    }

    /// Concrete checked-operation plan root.
    #[must_use]
    pub const fn plan_root(self) -> ContentHash {
        self.plan_root
    }
}

/// Opaque one-shot root-admission token.
///
/// The token is deliberately neither `Clone` nor `Copy`. Constructing it
/// validates the complete typed preflight against the caller's envelope;
/// [`Self::begin`] consumes it exactly once, so nested stages receive only
/// affine child leases and cannot reissue the admitted root authority.
#[derive(Debug)]
pub struct InvocationAdmission {
    invocation_id: ContentHash,
    plan_binding: Option<InvocationPlanBinding>,
    limits: InvocationLimits,
    required: InvocationResources,
}

/// One-shot admission issuer for one scientific invocation.
///
/// The issuer is consumed when it seals a plan, so a coordinator must create
/// a distinct issuer for a distinct invocation and cannot remint the same
/// invocation from one authority object.
#[derive(Debug, Default)]
pub struct InvocationAdmitter {
    _private: (),
}

impl InvocationAdmitter {
    /// Create one unused invocation issuer.
    #[must_use]
    pub const fn new() -> Self {
        Self { _private: () }
    }

    /// Seal a complete preflight into a one-use admission token and consume
    /// this issuer.
    ///
    /// # Errors
    /// Refuses the first insufficient resource in canonical dimensional order.
    pub fn admit(
        self,
        invocation_id: ContentHash,
        limits: InvocationLimits,
        required: InvocationResources,
    ) -> Result<InvocationAdmission, InvocationError> {
        self.admit_inner(invocation_id, None, limits, required)
    }

    /// Seal a complete preflight with an explicit checked-operation plan
    /// binding. Strong invocation-wide publication accepts only this path;
    /// [`Self::admit`] remains the clearly unbound compatibility entry.
    ///
    /// # Errors
    /// Refuses the first insufficient resource in canonical dimensional order.
    pub fn admit_bound(
        self,
        invocation_id: ContentHash,
        plan_binding: InvocationPlanBinding,
        limits: InvocationLimits,
        required: InvocationResources,
    ) -> Result<InvocationAdmission, InvocationError> {
        self.admit_inner(invocation_id, Some(plan_binding), limits, required)
    }

    fn admit_inner(
        self,
        invocation_id: ContentHash,
        plan_binding: Option<InvocationPlanBinding>,
        limits: InvocationLimits,
        required: InvocationResources,
    ) -> Result<InvocationAdmission, InvocationError> {
        limits.resources.checked_sub(required)?;
        Ok(InvocationAdmission {
            invocation_id,
            plan_binding,
            limits,
            required,
        })
    }
}

impl InvocationAdmission {
    /// Consume this admission and mint the sole root spend authority.
    ///
    /// # Errors
    /// Refuses an already-reached absolute deadline.
    pub fn begin<'clock>(
        self,
        cx: &'clock Cx<'_>,
        clock: &'clock dyn TimeSource,
    ) -> Result<InvocationBudget<'clock>, InvocationError> {
        InvocationBudget::admit(self, cx.cancel_gate(), cx.lease(), clock)
    }
}

/// Terminal disposition retained by a child or invocation receipt.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InvocationDisposition {
    /// Authority closed without a latched error. Unused capacity is returned;
    /// a caller claiming an exact plan must separately verify exact spend.
    Completed,
    /// Cancellation was requested or observed and the operation drained.
    Cancelled,
    /// A typed admission/runtime fault refused publication.
    Refused,
}

/// Fail-closed affine-ledger refusal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum InvocationError {
    /// One typed capacity was insufficient.
    ResourceExceeded {
        /// Stable resource name.
        resource: &'static str,
        /// Requested units.
        requested: u128,
        /// Available units.
        available: u128,
    },
    /// Checked accounting overflowed.
    ArithmeticOverflow {
        /// Stable resource name.
        resource: &'static str,
    },
    /// The absolute deadline was reached or passed.
    DeadlineExpired {
        /// Stable observing phase.
        phase: &'static str,
        /// Absolute deadline.
        deadline_ns: u64,
        /// Clock observation.
        observed_ns: u64,
    },
    /// Cancellation was observed after spending one poll opportunity.
    Cancelled {
        /// Stable observing phase.
        phase: &'static str,
    },
    /// The backing operation-memory lease refused a reservation.
    MemoryRefused {
        /// Stable allocation site.
        what: &'static str,
        /// Requested bytes.
        requested: u64,
        /// Bytes live at refusal.
        used: u64,
        /// Enforced limit.
        limit: u64,
    },
    /// A scientific phase explicitly refused its domain result.
    ExplicitRefusal {
        /// Stable refusing phase.
        phase: &'static str,
        /// Content identity of the structured domain refusal.
        reason: ContentHash,
    },
    /// A child phase label must be non-empty before identity derivation.
    EmptyPhase,
    /// A lease was used after terminal disposition.
    InactiveChild,
    /// A child cannot close while a nested child remains live.
    LiveNestedChildren {
        /// Number of unfinished descendants immediately below it.
        count: u64,
    },
    /// A child cannot close while memory reservations remain live.
    LiveMemoryReservations {
        /// Bytes still held.
        bytes: u64,
    },
    /// Root finalization found an unfinished child.
    UnfinishedChild {
        /// Deterministic child identity.
        child: ContentHash,
    },
    /// A finalizable child attempted the legacy close path and would have
    /// bypassed its reserved cleanup authority.
    FinalizationRequired,
    /// A nested child attempted to receive output authority owned by an
    /// enclosing invocation-wide transactional publication.
    TransactionalOutputScopeViolation {
        /// Finalizable ancestor that exclusively owns publication.
        ancestor: ContentHash,
        /// Rejected nested child phase.
        phase: &'static str,
        /// Nonzero output capacity requested by the nested child.
        requested: u64,
    },
    /// A finalizer operation was attempted before publication was prepared.
    PublicationNotPrepared,
    /// Publication was already committed or aborted for this child.
    PublicationAlreadySealed,
    /// Successful publication is forbidden after cancellation or refusal.
    PublicationForbidden,
    /// Finalization cannot seal until one mandatory terminal step completes.
    FinalizationIncomplete {
        /// Stable missing-step label.
        step: &'static str,
    },
    /// A finalization report could not be joined to the named child receipt.
    FinalizationReceiptMismatch {
        /// Stable mismatch label.
        invariant: &'static str,
    },
    /// The root invocation was already sealed and may only replay evidence.
    InvocationAlreadyFinalized,
    /// The backing lease observed an impossible release invariant.
    MemoryReleaseInvariant,
    /// Heap capacity for immutable evidence could not be reserved.
    ///
    /// This is a retryable producer-control refusal, not scientific failure
    /// evidence. It is never latched into a child or invocation receipt.
    EvidenceAllocationRefused {
        /// Stable evidence collection being prepared.
        what: &'static str,
        /// Number of typed items whose capacity was requested.
        requested_items: u64,
    },
}

impl core::fmt::Display for InvocationError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::ResourceExceeded {
                resource,
                requested,
                available,
            } => write!(
                formatter,
                "invocation resource `{resource}` requested {requested} units with {available} available"
            ),
            Self::ArithmeticOverflow { resource } => {
                write!(formatter, "invocation `{resource}` accounting overflowed")
            }
            Self::DeadlineExpired {
                phase,
                deadline_ns,
                observed_ns,
            } => write!(
                formatter,
                "invocation deadline {deadline_ns} ns expired during {phase} at {observed_ns} ns"
            ),
            Self::Cancelled { phase } => {
                write!(formatter, "invocation cancelled during {phase}")
            }
            Self::MemoryRefused {
                what,
                requested,
                used,
                limit,
            } => write!(
                formatter,
                "invocation memory refused {requested} B for `{what}` with {used}/{limit} B live"
            ),
            Self::ExplicitRefusal { phase, reason } => {
                write!(
                    formatter,
                    "invocation phase `{phase}` refused result {reason}"
                )
            }
            Self::EmptyPhase => formatter.write_str("invocation child phase must be non-empty"),
            Self::InactiveChild => formatter.write_str("invocation child is no longer active"),
            Self::LiveNestedChildren { count } => write!(
                formatter,
                "invocation child still owns {count} unfinished nested child lease(s)"
            ),
            Self::LiveMemoryReservations { bytes } => write!(
                formatter,
                "invocation child still owns {bytes} B of live memory reservations"
            ),
            Self::UnfinishedChild { child } => {
                write!(formatter, "invocation child {child} was not finalized")
            }
            Self::FinalizationRequired => formatter.write_str(
                "invocation child owns a reserved finalizer and must enter finalization before close",
            ),
            Self::TransactionalOutputScopeViolation {
                ancestor,
                phase,
                requested,
            } => write!(
                formatter,
                "invocation child phase `{phase}` requested {requested} output byte(s) inside transactional publication owned by {ancestor}"
            ),
            Self::PublicationNotPrepared => {
                formatter.write_str("invocation publication was not prepared")
            }
            Self::PublicationAlreadySealed => {
                formatter.write_str("invocation publication is already sealed")
            }
            Self::PublicationForbidden => formatter.write_str(
                "invocation success publication is forbidden after cancellation or refusal",
            ),
            Self::FinalizationIncomplete { step } => {
                write!(formatter, "invocation finalization is missing `{step}`")
            }
            Self::FinalizationReceiptMismatch { invariant } => write!(
                formatter,
                "invocation finalization receipt violated `{invariant}`"
            ),
            Self::InvocationAlreadyFinalized => {
                formatter.write_str("invocation was already finalized")
            }
            Self::MemoryReleaseInvariant => {
                formatter.write_str("invocation backing memory lease violated release accounting")
            }
            Self::EvidenceAllocationRefused {
                what,
                requested_items,
            } => write!(
                formatter,
                "invocation evidence `{what}` could not reserve capacity for {requested_items} item(s)"
            ),
        }
    }
}

impl core::error::Error for InvocationError {}

fn error_disposition(error: &InvocationError) -> InvocationDisposition {
    match error {
        InvocationError::DeadlineExpired { .. } | InvocationError::Cancelled { .. } => {
            InvocationDisposition::Cancelled
        }
        _ => InvocationDisposition::Refused,
    }
}

/// First backing-lease refusal retained in the immutable invocation receipt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InvocationMemoryRefusal {
    what: &'static str,
    requested: u64,
    used: u64,
    limit: u64,
}

impl InvocationMemoryRefusal {
    /// Stable allocation site.
    #[must_use]
    pub const fn what(&self) -> &'static str {
        self.what
    }

    /// Bytes requested by the refused reservation.
    #[must_use]
    pub const fn requested_bytes(&self) -> u64 {
        self.requested
    }

    /// Bytes live at refusal.
    #[must_use]
    pub const fn used_bytes(&self) -> u64 {
        self.used
    }

    /// Enforced backing limit.
    #[must_use]
    pub const fn limit_bytes(&self) -> u64 {
        self.limit
    }
}

/// Structured reason an immutable receipt failed semantic verification.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ReceiptSemanticError {
    /// The receipt schema is not understood by this verifier.
    UnsupportedVersion {
        /// Encountered schema version.
        found: u32,
    },
    /// The canonical root does not bind the retained fields.
    RootMismatch,
    /// A child receipt violated one named invariant.
    Child {
        /// Deterministic child ordinal.
        ordinal: u64,
        /// Stable invariant name.
        invariant: &'static str,
    },
    /// The invocation-level receipt violated one named invariant.
    Invocation {
        /// Stable invariant name.
        invariant: &'static str,
    },
    /// The receipt exceeds the verifier's explicit row-work envelope.
    WorkLimitExceeded {
        /// Child rows presented.
        children: u64,
        /// Maximum child rows admitted by this schema.
        limit: u64,
    },
    /// Scratch capacity for bounded semantic verification could not be
    /// reserved.
    AllocationRefused {
        /// Stable scratch collection name.
        what: &'static str,
        /// Typed items whose capacity was requested.
        requested_items: u64,
    },
}

impl core::fmt::Display for ReceiptSemanticError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::UnsupportedVersion { found } => {
                write!(formatter, "unsupported invocation receipt version {found}")
            }
            Self::RootMismatch => formatter.write_str("invocation receipt root mismatch"),
            Self::Child { ordinal, invariant } => write!(
                formatter,
                "invocation child ordinal {ordinal} violated `{invariant}`"
            ),
            Self::Invocation { invariant } => {
                write!(formatter, "invocation receipt violated `{invariant}`")
            }
            Self::WorkLimitExceeded { children, limit } => write!(
                formatter,
                "invocation receipt has {children} child rows with verifier limit {limit}"
            ),
            Self::AllocationRefused {
                what,
                requested_items,
            } => write!(
                formatter,
                "invocation receipt verifier could not reserve `{what}` for {requested_items} item(s)"
            ),
        }
    }
}

impl core::error::Error for ReceiptSemanticError {}

#[derive(Debug)]
struct ChildState {
    id: ContentHash,
    plan_binding: Option<InvocationPlanBinding>,
    parent: Option<usize>,
    ordinal: u64,
    phase: &'static str,
    granted: InvocationResources,
    remaining: InvocationResources,
    direct_consumed: InvocationResources,
    memory_current: u64,
    subtree_memory_current: u64,
    direct_memory_peak: u64,
    memory_peak: u64,
    memory_requested: u128,
    memory_released: u128,
    output_retained: u64,
    live_children: u64,
    finalization_required: bool,
    finalization_started: bool,
    finalization_granted: FinalizationResources,
    finalization_remaining: FinalizationResources,
    finalization_publication_scope: Option<InvocationPublicationScope>,
    finalization_publication: FinalizationPublication,
    finalization_report_root: Option<ContentHash>,
    failure: Option<InvocationError>,
    failure_inherited: bool,
    disposition: Option<InvocationDisposition>,
}

/// Recoverable authority class for an unfinished invocation child.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UnfinishedChildAuthority {
    /// Legacy/ordinary child budget.
    Ordinary,
    /// Finalizable child whose scientific handle never transitioned.
    FinalizableScientific,
    /// Child already transitioned into cleanup/finalization.
    Finalizer,
}

/// Deterministic read-only descriptor for abandoned child discovery.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UnfinishedChild {
    id: ContentHash,
    parent: Option<ContentHash>,
    ordinal: u64,
    phase: &'static str,
    authority: UnfinishedChildAuthority,
    live_children: u64,
    live_memory_bytes: u64,
    publication_scope: Option<InvocationPublicationScope>,
    publication: FinalizationPublication,
}

impl UnfinishedChild {
    /// Exact child identity accepted by the matching recovery API.
    #[must_use]
    pub const fn id(&self) -> ContentHash {
        self.id
    }

    /// Parent identity, or `None` for a root-level child.
    #[must_use]
    pub const fn parent(&self) -> Option<ContentHash> {
        self.parent
    }

    /// Deterministic global creation ordinal.
    #[must_use]
    pub const fn ordinal(&self) -> u64 {
        self.ordinal
    }

    /// Stable phase label.
    #[must_use]
    pub const fn phase(&self) -> &'static str {
        self.phase
    }

    /// Recovery API/type required for this child.
    #[must_use]
    pub const fn authority(&self) -> UnfinishedChildAuthority {
        self.authority
    }

    /// Direct unfinished descendants.
    #[must_use]
    pub const fn live_children(&self) -> u64 {
        self.live_children
    }

    /// Direct or descendant live logical memory.
    #[must_use]
    pub const fn live_memory_bytes(&self) -> u64 {
        self.live_memory_bytes
    }

    /// Current finalizer publication state.
    #[must_use]
    pub const fn publication(&self) -> FinalizationPublication {
        self.publication
    }

    /// Selected authority scope, if a terminal publication path was chosen.
    #[must_use]
    pub const fn publication_scope(&self) -> Option<InvocationPublicationScope> {
        self.publication_scope
    }
}

/// Immutable accounting evidence for one child lease.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ChildReceipt {
    id: ContentHash,
    parent: Option<ContentHash>,
    ordinal: u64,
    phase: &'static str,
    granted: InvocationResources,
    consumed: InvocationResources,
    direct_consumed: InvocationResources,
    returned: InvocationResources,
    direct_memory_peak: u64,
    memory_peak: u64,
    memory_requested: u128,
    memory_released: u128,
    output_retained: u64,
    finalization: Option<ChildFinalizationEvidence>,
    failure: Option<InvocationError>,
    failure_inherited: bool,
    disposition: InvocationDisposition,
    root: ContentHash,
}

impl ChildReceipt {
    /// Deterministic child identity.
    #[must_use]
    pub const fn id(&self) -> ContentHash {
        self.id
    }

    /// Parent child identity, or `None` for a root-level phase.
    #[must_use]
    pub const fn parent(&self) -> Option<ContentHash> {
        self.parent
    }

    /// Global deterministic issue ordinal.
    #[must_use]
    pub const fn ordinal(&self) -> u64 {
        self.ordinal
    }

    /// Stable phase name.
    #[must_use]
    pub const fn phase(&self) -> &'static str {
        self.phase
    }

    /// Original transferred allowance.
    #[must_use]
    pub const fn granted(&self) -> InvocationResources {
        self.granted
    }

    /// Permanently spent consumables and retained output.
    #[must_use]
    pub const fn consumed(&self) -> InvocationResources {
        self.consumed
    }

    /// Resources spent directly by this phase, excluding descendants.
    #[must_use]
    pub const fn direct_consumed(&self) -> InvocationResources {
        self.direct_consumed
    }

    /// Unused capacities returned exactly once.
    #[must_use]
    pub const fn returned(&self) -> InvocationResources {
        self.returned
    }

    /// Peak concurrent memory under this child.
    #[must_use]
    pub const fn memory_peak_bytes(&self) -> u64 {
        self.memory_peak
    }

    /// Peak bytes reserved directly by this phase, excluding descendants.
    #[must_use]
    pub const fn direct_memory_peak_bytes(&self) -> u64 {
        self.direct_memory_peak
    }

    /// Cumulative bytes directly reserved by this child.
    #[must_use]
    pub const fn memory_requested_bytes(&self) -> u128 {
        self.memory_requested
    }

    /// Cumulative direct reservations released by this child.
    #[must_use]
    pub const fn memory_released_bytes(&self) -> u128 {
        self.memory_released
    }

    /// Retained output bytes.
    #[must_use]
    pub const fn output_retained_bytes(&self) -> u64 {
        self.output_retained
    }

    /// Exact finalizer partition and report commitment, when this child used
    /// the post-cancel finalizer protocol.
    #[must_use]
    pub const fn finalization(&self) -> Option<&ChildFinalizationEvidence> {
        self.finalization.as_ref()
    }

    /// First latched runtime refusal, when this child did not complete.
    #[must_use]
    pub const fn failure(&self) -> Option<&InvocationError> {
        self.failure.as_ref()
    }

    /// Whether this child did not originate the failure: it was propagated
    /// from a descendant or copied from a pre-existing root failure.
    #[must_use]
    pub const fn failure_inherited(&self) -> bool {
        self.failure_inherited
    }

    /// Terminal child disposition.
    #[must_use]
    pub const fn disposition(&self) -> InvocationDisposition {
        self.disposition
    }

    /// Canonical child-receipt root.
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }
}

/// Terminal publication state sealed by a child finalizer.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FinalizationPublication {
    /// The finalizer has not yet decided whether output may be committed.
    Pending,
    /// A successful child passed its final terminal poll and owns a one-use
    /// prepared-publication token.
    Prepared,
    /// Output was deliberately left unchanged.
    Aborted,
    /// Output capacity was retained exactly once after the final poll.
    Committed {
        /// Retained output bytes.
        bytes: OutputBytes,
    },
}

/// Authority scope of a terminal finalizer publication decision.
///
/// `ChildLocal` preserves independently durable child output, but cannot prove
/// that the enclosing invocation completed successfully. `InvocationAtomic`
/// is reserved for the root two-phase protocol and is accepted as a strong
/// claim only when joined through a `FinalizedInvocationReceipt`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InvocationPublicationScope {
    /// The child decision may outlive a later root cancellation or refusal.
    ChildLocal,
    /// Destination mutation and root terminal evidence share one commit point.
    InvocationAtomic,
}

/// Exact finalizer partition committed into a version-2 child receipt.
///
/// These fields are producer-derived from the closed child state rather than
/// reconstructed from aggregate scientific consumption. A verifier can
/// therefore require equality with the corresponding finalization report
/// instead of accepting any cleanup claim that merely fits beneath generic
/// child totals.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ChildFinalizationEvidence {
    scientific_granted: InvocationResources,
    scientific_direct_consumed: InvocationResources,
    scientific_returned: InvocationResources,
    granted: FinalizationResources,
    consumed: FinalizationResources,
    returned: FinalizationResources,
    publication_scope: InvocationPublicationScope,
    publication: FinalizationPublication,
    report_root: ContentHash,
}

impl ChildFinalizationEvidence {
    /// Ordinary scientific allowance, excluding the isolated finalizer.
    #[must_use]
    pub const fn scientific_granted(&self) -> InvocationResources {
        self.scientific_granted
    }

    /// Ordinary resources spent directly by this child, excluding both
    /// descendants and finalizer consumption.
    #[must_use]
    pub const fn scientific_direct_consumed(&self) -> InvocationResources {
        self.scientific_direct_consumed
    }

    /// Ordinary unused capacity returned by the child.
    #[must_use]
    pub const fn scientific_returned(&self) -> InvocationResources {
        self.scientific_returned
    }

    /// Cleanup reserve isolated at child admission.
    #[must_use]
    pub const fn granted(&self) -> FinalizationResources {
        self.granted
    }

    /// Cleanup resources consumed by the finalizer.
    #[must_use]
    pub const fn consumed(&self) -> FinalizationResources {
        self.consumed
    }

    /// Cleanup resources returned when the child closed.
    #[must_use]
    pub const fn returned(&self) -> FinalizationResources {
        self.returned
    }

    /// Authority scope of the terminal publication decision.
    #[must_use]
    pub const fn publication_scope(&self) -> InvocationPublicationScope {
        self.publication_scope
    }

    /// Exact terminal publication state.
    #[must_use]
    pub const fn publication(&self) -> FinalizationPublication {
        self.publication
    }

    /// Canonical finalization-report root committed by the producer.
    #[must_use]
    pub const fn report_root(&self) -> ContentHash {
        self.report_root
    }
}

/// Immutable post-cancel cleanup evidence. It attests only the invocation
/// authority and resource steps performed by [`ChildFinalizer`]; a scheduler
/// or operation must separately supply its real worker-drain witness.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FinalizationReport {
    version: u32,
    plan_binding: Option<InvocationPlanBinding>,
    child: ContentHash,
    granted: FinalizationResources,
    consumed: FinalizationResources,
    returned: FinalizationResources,
    publication_scope: InvocationPublicationScope,
    publication: FinalizationPublication,
    failure: Option<InvocationError>,
    disposition: InvocationDisposition,
    root: ContentHash,
}

impl FinalizationReport {
    /// Report schema version.
    #[must_use]
    pub const fn version(&self) -> u32 {
        self.version
    }

    /// Explicit checked-plan binding, or `None` for the unbound compatibility
    /// admission path.
    #[must_use]
    pub const fn plan_binding(&self) -> Option<InvocationPlanBinding> {
        self.plan_binding
    }

    /// Exact invocation-child identity finalized by this report.
    #[must_use]
    pub const fn child(&self) -> ContentHash {
        self.child
    }

    /// Cleanup reserve isolated at child admission.
    #[must_use]
    pub const fn granted(&self) -> FinalizationResources {
        self.granted
    }

    /// Cleanup resources consumed after the finalizer transition.
    #[must_use]
    pub const fn consumed(&self) -> FinalizationResources {
        self.consumed
    }

    /// Unused cleanup capacity returned exactly once.
    #[must_use]
    pub const fn returned(&self) -> FinalizationResources {
        self.returned
    }

    /// Authority scope of this terminal publication decision.
    #[must_use]
    pub const fn publication_scope(&self) -> InvocationPublicationScope {
        self.publication_scope
    }

    /// Sealed publication decision.
    #[must_use]
    pub const fn publication(&self) -> FinalizationPublication {
        self.publication
    }

    /// First invocation failure, preserved unchanged through cleanup.
    #[must_use]
    pub const fn failure(&self) -> Option<&InvocationError> {
        self.failure.as_ref()
    }

    /// Terminal child disposition.
    #[must_use]
    pub const fn disposition(&self) -> InvocationDisposition {
        self.disposition
    }

    /// Canonical finalization-report root.
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }

    /// Verify the report-local resource and root invariants.
    #[must_use]
    pub fn verifies_integrity(&self) -> bool {
        finalization_report_semantics(self).is_ok()
    }

    /// Join this report to the exact semantically verified immutable child
    /// receipt. The returned receipt cannot be constructed from caller-authored
    /// success booleans.
    ///
    /// # Errors
    /// Refuses a malformed invocation receipt or any child/root/disposition,
    /// failure, resource, or publication mismatch.
    pub fn join(
        &self,
        invocation: &InvocationReceipt,
    ) -> Result<FinalizedChildReceipt, InvocationError> {
        invocation.verify_semantics().map_err(|_| {
            InvocationError::FinalizationReceiptMismatch {
                invariant: "invocation-receipt",
            }
        })?;
        finalization_report_semantics(self)?;
        if self.plan_binding != invocation.plan_binding {
            return Err(InvocationError::FinalizationReceiptMismatch {
                invariant: "plan-binding",
            });
        }
        let child = invocation
            .children()
            .iter()
            .find(|candidate| candidate.id() == self.child)
            .ok_or(InvocationError::FinalizationReceiptMismatch {
                invariant: "child-exists",
            })?;
        let child_finalization =
            child
                .finalization()
                .ok_or(InvocationError::FinalizationReceiptMismatch {
                    invariant: "child-finalization-present",
                })?;
        if child_finalization.granted() != self.granted
            || child_finalization.consumed() != self.consumed
            || child_finalization.returned() != self.returned
            || child_finalization.publication_scope() != self.publication_scope
            || child_finalization.publication() != self.publication
            || child_finalization.report_root() != self.root
        {
            return Err(InvocationError::FinalizationReceiptMismatch {
                invariant: "child-finalization-equality",
            });
        }
        if child.disposition() != self.disposition {
            return Err(InvocationError::FinalizationReceiptMismatch {
                invariant: "disposition",
            });
        }
        if child.failure() != self.failure.as_ref() {
            return Err(InvocationError::FinalizationReceiptMismatch {
                invariant: "first-failure",
            });
        }
        let granted = child.granted();
        if granted.work().get() < self.granted.work().get()
            || granted.polls().get() < self.granted.polls().get()
        {
            return Err(InvocationError::FinalizationReceiptMismatch {
                invariant: "cleanup-grant",
            });
        }
        let direct = child.direct_consumed();
        if direct.work().get() < self.consumed.work().get()
            || direct.polls().get() < self.consumed.polls().get()
        {
            return Err(InvocationError::FinalizationReceiptMismatch {
                invariant: "cleanup-consumption",
            });
        }
        match self.publication {
            FinalizationPublication::Committed { bytes }
                if child.output_retained_bytes() != bytes.get() =>
            {
                return Err(InvocationError::FinalizationReceiptMismatch {
                    invariant: "committed-output",
                });
            }
            FinalizationPublication::Aborted if child.output_retained_bytes() != 0 => {
                return Err(InvocationError::FinalizationReceiptMismatch {
                    invariant: "aborted-output",
                });
            }
            FinalizationPublication::Pending | FinalizationPublication::Prepared => {
                return Err(InvocationError::FinalizationReceiptMismatch {
                    invariant: "publication-terminal",
                });
            }
            FinalizationPublication::Committed { .. } | FinalizationPublication::Aborted => {}
        }
        let mut joined = FinalizedChildReceipt {
            invocation_root: invocation.root(),
            child: child.clone(),
            finalization: self.clone(),
            root: ContentHash([0; 32]),
        };
        joined.root = finalized_child_receipt_root(&joined);
        Ok(joined)
    }
}

/// Receipt that cryptographically joins sparse/domain finalization evidence to
/// one verified immutable invocation child.
///
/// ```compile_fail
/// fn forge_join(
///     invocation_root: fs_blake3::ContentHash,
///     child: fs_exec::ChildReceipt,
///     finalization: fs_exec::FinalizationReport,
///     root: fs_blake3::ContentHash,
/// ) -> fs_exec::FinalizedChildReceipt {
///     fs_exec::FinalizedChildReceipt {
///         invocation_root,
///         child,
///         finalization,
///         root,
///     }
/// }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FinalizedChildReceipt {
    invocation_root: ContentHash,
    child: ChildReceipt,
    finalization: FinalizationReport,
    root: ContentHash,
}

impl FinalizedChildReceipt {
    /// Parent invocation receipt root.
    #[must_use]
    pub const fn invocation_root(&self) -> ContentHash {
        self.invocation_root
    }

    /// Exact generic child receipt.
    #[must_use]
    pub const fn child(&self) -> &ChildReceipt {
        &self.child
    }

    /// Post-cancel finalization report.
    #[must_use]
    pub const fn finalization(&self) -> &FinalizationReport {
        &self.finalization
    }

    /// Canonical joined-receipt root.
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }

    /// Verify local component roots and their exact
    /// identity/disposition/failure binding.
    ///
    /// This structural check cannot authenticate the stored parent invocation
    /// root by itself. Use [`Self::verify_against`] when validating transported
    /// evidence against its parent receipt.
    #[must_use]
    pub fn verifies_integrity(&self) -> bool {
        self.child.root() == child_receipt_root(&self.child)
            && self.finalization.verifies_integrity()
            && self.child.finalization().is_some_and(|evidence| {
                evidence.granted() == self.finalization.granted()
                    && evidence.consumed() == self.finalization.consumed()
                    && evidence.returned() == self.finalization.returned()
                    && evidence.publication_scope() == self.finalization.publication_scope()
                    && evidence.publication() == self.finalization.publication()
                    && evidence.report_root() == self.finalization.root()
            })
            && self.child.id() == self.finalization.child()
            && self.child.disposition() == self.finalization.disposition()
            && self.child.failure() == self.finalization.failure()
            && self.root == finalized_child_receipt_root(self)
    }

    /// Verify this joined receipt against the complete parent invocation.
    ///
    /// # Errors
    /// Refuses a malformed or substituted invocation, a missing/different
    /// child, any finalizer mismatch, or a rehashed joined receipt whose
    /// components do not exactly equal a fresh verified join.
    pub fn verify_against(&self, invocation: &InvocationReceipt) -> Result<(), InvocationError> {
        if invocation.root() != self.invocation_root {
            return Err(InvocationError::FinalizationReceiptMismatch {
                invariant: "invocation-root",
            });
        }
        let expected = self.finalization.join(invocation)?;
        if expected != *self {
            return Err(InvocationError::FinalizationReceiptMismatch {
                invariant: "joined-receipt",
            });
        }
        Ok(())
    }
}

/// Immutable terminal receipt.  This value is cloneable evidence; it contains
/// no live resource authority.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InvocationReceipt {
    version: u32,
    invocation_id: ContentHash,
    plan_binding: Option<InvocationPlanBinding>,
    limits: InvocationLimits,
    required: InvocationResources,
    remaining: InvocationResources,
    children: Vec<ChildReceipt>,
    last_deadline_observation: Option<Time>,
    memory_peak: u64,
    memory_requested: u128,
    memory_released: u128,
    memory_refusals: u128,
    memory_first_refusal: Option<InvocationMemoryRefusal>,
    output_retained: u64,
    failure: Option<InvocationError>,
    failure_origin: Option<ContentHash>,
    disposition: InvocationDisposition,
    root: ContentHash,
}

impl InvocationReceipt {
    /// Receipt schema version.
    #[must_use]
    pub const fn version(&self) -> u32 {
        self.version
    }

    /// Stable invocation identity.
    #[must_use]
    pub const fn invocation_id(&self) -> ContentHash {
        self.invocation_id
    }

    /// Explicit checked-plan binding, or `None` for an unbound compatibility
    /// invocation.
    #[must_use]
    pub const fn plan_binding(&self) -> Option<InvocationPlanBinding> {
        self.plan_binding
    }

    /// Admitted limits and immutable obligations.
    #[must_use]
    pub const fn limits(&self) -> &InvocationLimits {
        &self.limits
    }

    /// Preflight requirement.
    #[must_use]
    pub const fn required(&self) -> InvocationResources {
        self.required
    }

    /// Unspent capacity at terminal finalization.
    #[must_use]
    pub const fn remaining(&self) -> InvocationResources {
        self.remaining
    }

    /// Ordered child receipts.
    #[must_use]
    pub fn children(&self) -> &[ChildReceipt] {
        &self.children
    }

    /// Peak backing-memory live set.
    #[must_use]
    pub const fn memory_peak_bytes(&self) -> u64 {
        self.memory_peak
    }

    /// Cumulative bytes admitted by the backing memory lease.
    #[must_use]
    pub const fn memory_requested_bytes(&self) -> u128 {
        self.memory_requested
    }

    /// Cumulative bytes released by completed RAII reservations.
    #[must_use]
    pub const fn memory_released_bytes(&self) -> u128 {
        self.memory_released
    }

    /// Count of backing memory refusals retained by this transaction.
    #[must_use]
    pub const fn memory_refusals(&self) -> u128 {
        self.memory_refusals
    }

    /// First backing memory refusal, when any occurred.
    #[must_use]
    pub const fn memory_first_refusal(&self) -> Option<&InvocationMemoryRefusal> {
        self.memory_first_refusal.as_ref()
    }

    /// Last logical-clock observation made for deadline enforcement.
    #[must_use]
    pub const fn last_deadline_observation(&self) -> Option<Time> {
        self.last_deadline_observation
    }

    /// Retained output bytes.
    #[must_use]
    pub const fn output_retained_bytes(&self) -> u64 {
        self.output_retained
    }

    /// First latched transaction failure, when terminal disposition is not
    /// completed.
    #[must_use]
    pub const fn failure(&self) -> Option<&InvocationError> {
        self.failure.as_ref()
    }

    /// Child that first latched the root failure, or `None` when the failure
    /// originated at root admission/finalization.
    #[must_use]
    pub const fn failure_origin(&self) -> Option<ContentHash> {
        self.failure_origin
    }

    /// Terminal disposition.
    #[must_use]
    pub const fn disposition(&self) -> InvocationDisposition {
        self.disposition
    }

    /// Canonical accounting root.
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }

    /// Recompute the canonical accounting root and all typed conservation,
    /// topology, memory, output, and disposition invariants.
    #[must_use]
    pub fn verifies_integrity(&self) -> bool {
        self.verify_semantics().is_ok()
    }

    /// Verify the canonical root and the complete affine receipt semantics.
    ///
    /// # Errors
    /// Returns the first invariant failure in deterministic schema order.
    pub fn verify_semantics(&self) -> Result<(), ReceiptSemanticError> {
        verify_receipt_semantics(self)
    }
}

/// Non-cloneable root invocation authority.
pub struct InvocationBudget<'clock> {
    invocation_id: ContentHash,
    plan_binding: Option<InvocationPlanBinding>,
    limits: InvocationLimits,
    required: InvocationResources,
    remaining: InvocationResources,
    clock: &'clock dyn TimeSource,
    cancel_gate: &'clock CancelGate,
    last_deadline_observation: Option<Time>,
    _ambient_memory: Option<LeaseCharge>,
    backing_memory: OperationMemoryLease,
    children: Vec<ChildState>,
    next_ordinal: u64,
    failure: Option<InvocationError>,
    failure_origin: Option<ContentHash>,
    sealed_receipt: Option<InvocationReceipt>,
}

impl core::fmt::Debug for InvocationBudget<'_> {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        formatter
            .debug_struct("InvocationBudget")
            .field("invocation_id", &self.invocation_id)
            .field("limits", &self.limits)
            .field("required", &self.required)
            .field("remaining", &self.remaining)
            .field("children", &self.children.len())
            .finish_non_exhaustive()
    }
}

impl<'clock> InvocationBudget<'clock> {
    /// Admit a complete plan before any child can spend authority.
    ///
    /// # Errors
    /// Refuses the first resource in fixed work, poll, cost, evaluation,
    /// memory, output order, then an already-expired deadline.
    fn admit(
        admission: InvocationAdmission,
        cancel_gate: &'clock CancelGate,
        ambient_memory: Option<&OperationMemoryLease>,
        clock: &'clock dyn TimeSource,
    ) -> Result<Self, InvocationError> {
        let InvocationAdmission {
            invocation_id,
            plan_binding,
            limits,
            required,
        } = admission;
        let last_deadline_observation = if let Some(deadline) = limits.deadline {
            let now = clock.now();
            if now >= deadline {
                return Err(InvocationError::DeadlineExpired {
                    phase: "invocation-admission",
                    deadline_ns: deadline.as_nanos(),
                    observed_ns: now.as_nanos(),
                });
            }
            Some(now)
        } else {
            None
        };
        let ambient_memory = ambient_memory
            .map(|lease| lease.reserve("invocation-root-memory", required.memory.0))
            .transpose()
            .map_err(|refusal| InvocationError::MemoryRefused {
                what: refusal.what,
                requested: refusal.requested_bytes,
                used: refusal.used_bytes,
                limit: refusal.limit_bytes,
            })?;
        Ok(Self {
            invocation_id,
            plan_binding,
            limits,
            required,
            remaining: required,
            clock,
            cancel_gate,
            last_deadline_observation,
            _ambient_memory: ambient_memory,
            backing_memory: OperationMemoryLease::bounded(required.memory.0),
            children: Vec::new(),
            next_ordinal: 0,
            failure: None,
            failure_origin: None,
            sealed_receipt: None,
        })
    }

    /// Transfer an exact affine allowance to one sequential child.
    ///
    /// # Errors
    /// Refuses an empty phase, insufficient capacity, or ordinal overflow
    /// before mutation.
    pub fn split_child<'budget>(
        &'budget mut self,
        phase: &'static str,
        grant: InvocationResources,
    ) -> Result<ChildBudget<'budget, 'clock>, InvocationError> {
        let node = self.open_child(None, phase, grant, FinalizationResources::default(), false)?;
        Ok(ChildBudget { owner: self, node })
    }

    /// Transfer ordinary scientific resources plus an isolated terminal
    /// cleanup reserve to one sequential child.
    ///
    /// The returned [`FinalizableChildBudget`] cannot spend the finalization
    /// reserve. It must transition through
    /// [`FinalizableChildBudget::begin_finalization`] before closing,
    /// including on cancellation or refusal.
    ///
    /// # Errors
    /// Refuses an empty phase, insufficient combined capacity, or checked
    /// resource overflow before mutation.
    pub fn split_finalizable_child<'budget>(
        &'budget mut self,
        phase: &'static str,
        grant: InvocationResources,
        finalization: FinalizationResources,
    ) -> Result<FinalizableChildBudget<'budget, 'clock>, InvocationError> {
        let node = self.open_child(None, phase, grant, finalization, true)?;
        Ok(FinalizableChildBudget {
            child: ChildBudget { owner: self, node },
        })
    }

    fn open_child(
        &mut self,
        parent: Option<usize>,
        phase: &'static str,
        grant: InvocationResources,
        finalization: FinalizationResources,
        finalization_required: bool,
    ) -> Result<usize, InvocationError> {
        if self.sealed_receipt.is_some() {
            return Err(InvocationError::InvocationAlreadyFinalized);
        }
        if let Some(error) = parent
            .and_then(|index| self.children.get(index))
            .and_then(|state| state.failure.clone())
            .or_else(|| self.failure.clone())
        {
            return Err(error);
        }
        if phase.is_empty() {
            return Err(InvocationError::EmptyPhase);
        }
        let child_count = u64::try_from(self.children.len()).map_err(|_| {
            InvocationError::ArithmeticOverflow {
                resource: "child-count",
            }
        })?;
        if child_count >= INVOCATION_RECEIPT_MAX_CHILDREN {
            let error = InvocationError::ResourceExceeded {
                resource: "child-count",
                requested: u128::from(child_count) + 1,
                available: u128::from(INVOCATION_RECEIPT_MAX_CHILDREN),
            };
            self.latch_failure(parent, error.clone());
            return Err(error);
        }
        let ordinal = self.next_ordinal;
        let next_ordinal = match ordinal.checked_add(1) {
            Some(next) => next,
            None => {
                let error = InvocationError::ArithmeticOverflow {
                    resource: "child-ordinal",
                };
                self.latch_failure(parent, error.clone());
                return Err(error);
            }
        };
        let total_grant = match grant.checked_add(finalization.as_invocation_resources()) {
            Ok(total) => total,
            Err(error) => {
                self.latch_failure(parent, error.clone());
                return Err(error);
            }
        };
        let mut ancestor = parent;
        let mut transactional_ancestor = None;
        while let Some(index) = ancestor {
            let state = self
                .children
                .get(index)
                .ok_or(InvocationError::InactiveChild)?;
            if state.finalization_required {
                transactional_ancestor = Some(state.id);
                break;
            }
            ancestor = state.parent;
        }
        if let Some(ancestor) = transactional_ancestor.filter(|_| total_grant.output.get() != 0) {
            let error = InvocationError::TransactionalOutputScopeViolation {
                ancestor,
                phase,
                requested: total_grant.output.get(),
            };
            self.latch_failure(parent, error.clone());
            return Err(error);
        }
        let available = match parent {
            Some(index) => {
                let state = self
                    .children
                    .get(index)
                    .ok_or(InvocationError::InactiveChild)?;
                if state.disposition.is_some() {
                    return Err(InvocationError::InactiveChild);
                }
                state.remaining
            }
            None => self.remaining,
        };
        if let Some(index) = parent {
            let direct_live = self.children[index].memory_current;
            let allocatable = available
                .memory
                .0
                .checked_sub(direct_live)
                .ok_or(InvocationError::MemoryReleaseInvariant)?;
            if total_grant.memory.0 > allocatable {
                let error = exceeded(
                    "memory-bytes",
                    u128::from(total_grant.memory.0),
                    u128::from(allocatable),
                );
                self.latch_failure(parent, error.clone());
                return Err(error);
            }
        }
        let remaining = match available.checked_sub(total_grant) {
            Ok(remaining) => remaining,
            Err(error) => {
                self.latch_failure(parent, error.clone());
                return Err(error);
            }
        };
        let live_children = match parent {
            Some(index) => match self.children[index].live_children.checked_add(1) {
                Some(count) => count,
                None => {
                    let error = InvocationError::ArithmeticOverflow {
                        resource: "live-children",
                    };
                    self.latch_failure(parent, error.clone());
                    return Err(error);
                }
            },
            None => 0,
        };
        let parent_id = parent.map(|index| self.children[index].id);
        let id = child_id(
            self.invocation_id,
            self.plan_binding,
            parent_id,
            ordinal,
            phase,
            grant,
            finalization,
        );
        match parent {
            Some(index) => {
                let parent_state = &mut self.children[index];
                parent_state.remaining = remaining;
                parent_state.live_children = live_children;
            }
            None => self.remaining = remaining,
        }
        self.next_ordinal = next_ordinal;
        let node = self.children.len();
        self.children.push(ChildState {
            id,
            plan_binding: self.plan_binding,
            parent,
            ordinal,
            phase,
            granted: total_grant,
            remaining: grant,
            direct_consumed: InvocationResources::default(),
            memory_current: 0,
            subtree_memory_current: 0,
            direct_memory_peak: 0,
            memory_peak: 0,
            memory_requested: 0,
            memory_released: 0,
            output_retained: 0,
            live_children: 0,
            finalization_required,
            finalization_started: false,
            finalization_granted: finalization,
            finalization_remaining: finalization,
            finalization_publication_scope: None,
            finalization_publication: FinalizationPublication::Pending,
            finalization_report_root: None,
            failure: None,
            failure_inherited: false,
            disposition: None,
        });
        Ok(node)
    }

    fn latch_failure(&mut self, mut node: Option<usize>, error: InvocationError) {
        let origin = node.map(|index| self.children[index].id);
        let origin_node = node;
        let root_already_failed = self.failure.is_some();
        while let Some(index) = node {
            let state = &mut self.children[index];
            if state.failure.is_none() {
                state.failure = Some(error.clone());
                state.failure_inherited = root_already_failed || Some(index) != origin_node;
            }
            node = state.parent;
        }
        if self.failure.is_none() {
            self.failure = Some(error);
            self.failure_origin = origin;
        }
    }

    fn close_child(&mut self, node: usize) -> Result<InvocationDisposition, InvocationError> {
        let (parent, returned, disposition) = {
            let state = self
                .children
                .get(node)
                .ok_or(InvocationError::InactiveChild)?;
            if state.disposition.is_some() {
                return Err(InvocationError::InactiveChild);
            }
            if state.live_children != 0 {
                return Err(InvocationError::LiveNestedChildren {
                    count: state.live_children,
                });
            }
            if state.memory_current != 0 {
                return Err(InvocationError::LiveMemoryReservations {
                    bytes: state.memory_current,
                });
            }
            if state.subtree_memory_current != 0 {
                return Err(InvocationError::MemoryReleaseInvariant);
            }
            let returned = state
                .remaining
                .checked_add(state.finalization_remaining.as_invocation_resources())?;
            (
                state.parent,
                returned,
                state
                    .failure
                    .as_ref()
                    .map_or(InvocationDisposition::Completed, error_disposition),
            )
        };
        match parent {
            Some(index) => {
                let parent_state = &mut self.children[index];
                parent_state.remaining = parent_state.remaining.checked_add(returned)?;
                parent_state.live_children = parent_state.live_children.checked_sub(1).ok_or(
                    InvocationError::ArithmeticOverflow {
                        resource: "live-children",
                    },
                )?;
            }
            None => self.remaining = self.remaining.checked_add(returned)?,
        }
        self.children[node].disposition = Some(disposition);
        Ok(disposition)
    }

    fn current_failure(&self, node: Option<usize>) -> Option<InvocationError> {
        node.and_then(|index| self.children[index].failure.clone())
            .or_else(|| self.failure.clone())
    }

    fn observe_deadline(
        &mut self,
        node: Option<usize>,
        phase: &'static str,
    ) -> Result<(), InvocationError> {
        if let Some(error) = self.current_failure(node) {
            return Err(error);
        }
        let Some(deadline) = self.limits.deadline else {
            return Ok(());
        };
        let now = self.clock.now();
        self.last_deadline_observation = Some(now);
        if now < deadline {
            return Ok(());
        }
        let error = InvocationError::DeadlineExpired {
            phase,
            deadline_ns: deadline.as_nanos(),
            observed_ns: now.as_nanos(),
        };
        self.cancel_gate.request();
        self.latch_failure(node, error.clone());
        Err(error)
    }

    fn observe_cancellation(
        &mut self,
        node: Option<usize>,
        phase: &'static str,
    ) -> Result<(), InvocationError> {
        if let Some(error) = self.current_failure(node) {
            return Err(error);
        }
        if !self.cancel_gate.is_requested() {
            return Ok(());
        }
        let error = InvocationError::Cancelled { phase };
        self.latch_failure(node, error.clone());
        Err(error)
    }

    fn observe_terminal(
        &mut self,
        node: Option<usize>,
        phase: &'static str,
    ) -> Result<(), InvocationError> {
        self.observe_deadline(node, phase)?;
        self.observe_cancellation(node, phase)
    }

    fn unfinished_child_descriptor(&self, node: usize) -> UnfinishedChild {
        let state = &self.children[node];
        let authority = if state.finalization_required {
            if state.finalization_started {
                UnfinishedChildAuthority::Finalizer
            } else {
                UnfinishedChildAuthority::FinalizableScientific
            }
        } else {
            UnfinishedChildAuthority::Ordinary
        };
        UnfinishedChild {
            id: state.id,
            parent: state.parent.map(|parent| self.children[parent].id),
            ordinal: state.ordinal,
            phase: state.phase,
            authority,
            live_children: state.live_children,
            live_memory_bytes: state.memory_current.max(state.subtree_memory_current),
            publication_scope: state.finalization_publication_scope,
            publication: state.finalization_publication,
        }
    }

    /// All unfinished children in deterministic parent-before-child creation
    /// order. This lets a caller recover after a panic even when no child id
    /// was copied before unwind.
    #[must_use]
    pub fn unfinished_children(&self) -> Vec<UnfinishedChild> {
        self.children
            .iter()
            .enumerate()
            .filter(|(_, state)| state.disposition.is_none())
            .map(|(node, _)| self.unfinished_child_descriptor(node))
            .collect()
    }

    /// Deterministic deepest recoverable child.
    ///
    /// Append-only topology places descendants after parents. Reverse-order
    /// leaf selection therefore exposes the next inside-out recovery target
    /// instead of repeatedly returning a blocked ancestor.
    #[must_use]
    pub fn next_unfinished_child(&self) -> Option<UnfinishedChild> {
        self.children
            .iter()
            .enumerate()
            .rev()
            .find(|(_, state)| state.disposition.is_none() && state.live_children == 0)
            .or_else(|| {
                self.children
                    .iter()
                    .enumerate()
                    .rev()
                    .find(|(_, state)| state.disposition.is_none())
            })
            .map(|(node, _)| self.unfinished_child_descriptor(node))
    }

    /// Replay immutable finalization evidence for an already sealed child.
    ///
    /// This is evidence replay only: it does not re-open authority, return
    /// resources, or repeat cleanup/publication. The recomputed report must
    /// match the report root committed into the child state.
    ///
    /// # Errors
    /// Refuses an unknown child, a legacy/unsealed child, or any stored
    /// finalization commitment mismatch.
    pub fn replay_child_finalization(
        &self,
        child: ContentHash,
    ) -> Result<FinalizationReport, InvocationError> {
        let state = self.children.iter().find(|state| state.id == child).ok_or(
            InvocationError::FinalizationReceiptMismatch {
                invariant: "child-exists",
            },
        )?;
        let report = finalization_report_from_state(state)?;
        finalization_report_semantics(&report)?;
        if state.finalization_report_root != Some(report.root()) {
            return Err(InvocationError::FinalizationReceiptMismatch {
                invariant: "stored-finalization-root",
            });
        }
        Ok(report)
    }

    /// Recover fail-closed finalizer authority after a caught unwind or an
    /// otherwise abandoned child handle.
    ///
    /// Mutable access to the root proves the prior child/finalizer borrow no
    /// longer exists, so this cannot duplicate live authority. An unfinished,
    /// non-committed child first inherits any existing ancestor failure,
    /// observes deadline/cancellation, or latches the supplied unwind refusal.
    /// A child whose atomic publication already committed is recovered only
    /// for evidence sealing/replay and is not retroactively relabeled.
    ///
    /// # Errors
    /// Refuses an empty phase, unknown/legacy child identity, or malformed
    /// already-closed state.
    pub fn recover_child_finalizer<'budget>(
        &'budget mut self,
        child: ContentHash,
        phase: &'static str,
        reason: ContentHash,
    ) -> Result<ChildFinalizer<'budget, 'clock>, InvocationError> {
        if self.sealed_receipt.is_some() {
            return Err(InvocationError::InvocationAlreadyFinalized);
        }
        if phase.is_empty() {
            return Err(InvocationError::EmptyPhase);
        }
        let node = self
            .children
            .iter()
            .position(|state| state.id == child)
            .ok_or(InvocationError::FinalizationReceiptMismatch {
                invariant: "child-exists",
            })?;
        if !self.children[node].finalization_required {
            return Err(InvocationError::FinalizationReceiptMismatch {
                invariant: "child-is-finalizable",
            });
        }
        if self.children[node].disposition.is_some() {
            let report = finalization_report_from_state(&self.children[node])?;
            if self.children[node].finalization_report_root != Some(report.root()) {
                return Err(InvocationError::FinalizationReceiptMismatch {
                    invariant: "stored-finalization-root",
                });
            }
            return Ok(ChildFinalizer { owner: self, node });
        }
        self.children[node].finalization_started = true;
        if !matches!(
            self.children[node].finalization_publication,
            FinalizationPublication::Committed { .. }
        ) {
            if self.children[node].failure.is_none() {
                if let Some(failure) = self.current_failure(Some(node)) {
                    self.latch_failure(Some(node), failure);
                } else {
                    let _ = self.observe_terminal(Some(node), phase);
                }
            }
            if self.children[node].failure.is_none() {
                self.latch_failure(
                    Some(node),
                    InvocationError::ExplicitRefusal { phase, reason },
                );
            }
        }
        Ok(ChildFinalizer { owner: self, node })
    }

    /// Recover an abandoned ordinary child fail-closed.
    ///
    /// Mutable root access proves the previous `ChildBudget` borrow no longer
    /// exists. The recovered child inherits the existing first failure,
    /// observes deadline/cancellation in normal precedence order, or latches
    /// the supplied unwind refusal before any authority is returned.
    ///
    /// # Errors
    /// Refuses an empty phase, unknown/closed child, or a finalizable child
    /// (which must use [`Self::recover_child_finalizer`]).
    pub fn recover_child_budget<'budget>(
        &'budget mut self,
        child: ContentHash,
        phase: &'static str,
        reason: ContentHash,
    ) -> Result<ChildBudget<'budget, 'clock>, InvocationError> {
        if self.sealed_receipt.is_some() {
            return Err(InvocationError::InvocationAlreadyFinalized);
        }
        if phase.is_empty() {
            return Err(InvocationError::EmptyPhase);
        }
        let node = self
            .children
            .iter()
            .position(|state| state.id == child)
            .ok_or(InvocationError::FinalizationReceiptMismatch {
                invariant: "child-exists",
            })?;
        if self.children[node].finalization_required {
            return Err(InvocationError::FinalizationRequired);
        }
        if self.children[node].disposition.is_some() {
            return Err(InvocationError::InactiveChild);
        }
        if self.children[node].failure.is_none() {
            if let Some(failure) = self.current_failure(Some(node)) {
                self.latch_failure(Some(node), failure);
            } else {
                let _ = self.observe_terminal(Some(node), phase);
            }
        }
        if self.children[node].failure.is_none() {
            self.latch_failure(
                Some(node),
                InvocationError::ExplicitRefusal { phase, reason },
            );
        }
        Ok(ChildBudget { owner: self, node })
    }

    /// Seal or exactly replay a terminal immutable receipt.
    ///
    /// # Errors
    /// Refuses unfinished children or a backing-memory invariant violation
    /// without consuming root authority, so abandoned children can be
    /// recovered and closed before retrying.
    pub fn finish(&mut self) -> Result<InvocationReceipt, InvocationError> {
        if let Some(receipt) = &self.sealed_receipt {
            return try_clone_invocation_receipt(receipt, "sealed-invocation-replay");
        }
        if let Some(unfinished) = self.next_unfinished_child() {
            return Err(InvocationError::UnfinishedChild {
                child: unfinished.id(),
            });
        }
        let _ = self.observe_terminal(None, "invocation-finalize");
        let memory = self.backing_memory.receipt();
        let (memory_requested, memory_released) =
            self.children
                .iter()
                .try_fold((0_u128, 0_u128), |(requested, released), state| {
                    Ok::<_, InvocationError>((
                        requested.checked_add(state.memory_requested).ok_or(
                            InvocationError::ArithmeticOverflow {
                                resource: "memory-requested",
                            },
                        )?,
                        released.checked_add(state.memory_released).ok_or(
                            InvocationError::ArithmeticOverflow {
                                resource: "memory-released",
                            },
                        )?,
                    ))
                })?;
        if memory.used_bytes != 0
            || memory.release_invariant_violations != 0
            || memory_requested != memory_released
            || memory_requested != memory.requested_bytes
        {
            return Err(InvocationError::MemoryReleaseInvariant);
        }
        let mut children = try_evidence_vec("invocation-child-receipts", self.children.len())?;
        for state in &self.children {
            children.push(child_receipt(&self.children, state)?);
        }
        let output_retained = children.iter().try_fold(0_u64, |sum, child| {
            sum.checked_add(child.output_retained)
                .ok_or(InvocationError::ArithmeticOverflow {
                    resource: "output-retained",
                })
        })?;
        let memory_first_refusal =
            memory
                .first_refusal
                .as_ref()
                .map(|refusal| InvocationMemoryRefusal {
                    what: refusal.what,
                    requested: refusal.requested_bytes,
                    used: refusal.used_bytes,
                    limit: refusal.limit_bytes,
                });
        let disposition = self
            .failure
            .as_ref()
            .map_or(InvocationDisposition::Completed, error_disposition);
        let mut receipt = InvocationReceipt {
            version: INVOCATION_RECEIPT_VERSION,
            invocation_id: self.invocation_id,
            plan_binding: self.plan_binding,
            limits: self.limits.clone(),
            required: self.required,
            remaining: self.remaining,
            children,
            last_deadline_observation: self.last_deadline_observation,
            memory_peak: memory.peak_bytes,
            memory_requested,
            memory_released,
            memory_refusals: memory.refusals,
            memory_first_refusal,
            output_retained,
            failure: self.failure.clone(),
            failure_origin: self.failure_origin,
            disposition,
            root: ContentHash([0; 32]),
        };
        receipt.root = invocation_receipt_root(&receipt);
        receipt
            .verify_semantics()
            .map_err(|_| InvocationError::FinalizationReceiptMismatch {
                invariant: "producer-invocation-receipt",
            })?;
        let returned = try_clone_invocation_receipt(&receipt, "invocation-receipt-return")?;
        // The ambient parent charge covered the root capacity until this exact
        // terminal cut. Release it once; replay uses the cached receipt.
        self._ambient_memory.take();
        self.sealed_receipt = Some(receipt);
        Ok(returned)
    }
}

/// Non-cloneable affine child authority. `finish` consumes it, returning unused
/// capacities exactly once to its parent.
pub struct ChildBudget<'budget, 'clock> {
    owner: &'budget mut InvocationBudget<'clock>,
    node: usize,
}

impl core::fmt::Debug for ChildBudget<'_, '_> {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        formatter
            .debug_struct("ChildBudget")
            .field("id", &self.owner.children[self.node].id)
            .field("phase", &self.owner.children[self.node].phase)
            .field("remaining", &self.owner.children[self.node].remaining)
            .finish_non_exhaustive()
    }
}

impl<'budget, 'clock> ChildBudget<'budget, 'clock> {
    /// Deterministic child identity.
    #[must_use]
    pub fn id(&self) -> ContentHash {
        self.owner.children[self.node].id
    }

    /// Remaining typed capacity.
    #[must_use]
    pub fn remaining(&self) -> InvocationResources {
        self.owner.children[self.node].remaining
    }

    /// Split a nested affine child from this child's remaining capacity.
    ///
    /// # Errors
    /// Refuses an empty phase or insufficient capacity before mutation.
    pub fn split_child<'child>(
        &'child mut self,
        phase: &'static str,
        grant: InvocationResources,
    ) -> Result<ChildBudget<'child, 'clock>, InvocationError> {
        let node = self.owner.open_child(
            Some(self.node),
            phase,
            grant,
            FinalizationResources::default(),
            false,
        )?;
        Ok(ChildBudget {
            owner: &mut *self.owner,
            node,
        })
    }

    /// Split a nested child with an isolated terminal cleanup reserve.
    ///
    /// # Errors
    /// Refuses empty phase, checked overflow, or insufficient parent
    /// scientific capacity before mutation.
    pub fn split_finalizable_child<'child>(
        &'child mut self,
        phase: &'static str,
        grant: InvocationResources,
        finalization: FinalizationResources,
    ) -> Result<FinalizableChildBudget<'child, 'clock>, InvocationError> {
        let node = self
            .owner
            .open_child(Some(self.node), phase, grant, finalization, true)?;
        Ok(FinalizableChildBudget {
            child: ChildBudget {
                owner: &mut *self.owner,
                node,
            },
        })
    }

    fn ensure_active(&self) -> Result<(), InvocationError> {
        if self.owner.children[self.node].disposition.is_some() {
            Err(InvocationError::InactiveChild)
        } else if let Some(error) = self.owner.children[self.node].failure.clone() {
            Err(error)
        } else {
            Ok(())
        }
    }

    fn latch(&mut self, error: InvocationError) -> InvocationError {
        self.owner.latch_failure(Some(self.node), error.clone());
        error
    }

    /// Spend declared logical work.
    ///
    /// # Errors
    /// Refuses over-consumption or stale authority.
    pub fn charge_work(&mut self, amount: WorkUnits) -> Result<(), InvocationError> {
        self.ensure_active()?;
        let state = &self.owner.children[self.node];
        let remaining = match state.remaining.work.0.checked_sub(amount.0) {
            Some(remaining) => remaining,
            None => {
                let available = state.remaining.work.0;
                return Err(self.latch(exceeded("work", amount.0, available)));
            }
        };
        let direct = match state.direct_consumed.work.0.checked_add(amount.0) {
            Some(direct) => direct,
            None => {
                return Err(self.latch(InvocationError::ArithmeticOverflow { resource: "work" }));
            }
        };
        let state = &mut self.owner.children[self.node];
        state.remaining.work.0 = remaining;
        state.direct_consumed.work.0 = direct;
        Ok(())
    }

    /// Spend abstract cost.
    ///
    /// # Errors
    /// Refuses over-consumption or stale authority.
    pub fn charge_cost(&mut self, amount: CostUnits) -> Result<(), InvocationError> {
        self.ensure_active()?;
        let state = &self.owner.children[self.node];
        let remaining = match state.remaining.cost.0.checked_sub(amount.0) {
            Some(remaining) => remaining,
            None => {
                let available = state.remaining.cost.0;
                return Err(self.latch(exceeded(
                    "cost",
                    u128::from(amount.0),
                    u128::from(available),
                )));
            }
        };
        let direct = match state.direct_consumed.cost.0.checked_add(amount.0) {
            Some(direct) => direct,
            None => {
                return Err(self.latch(InvocationError::ArithmeticOverflow { resource: "cost" }));
            }
        };
        let state = &mut self.owner.children[self.node];
        state.remaining.cost.0 = remaining;
        state.direct_consumed.cost.0 = direct;
        Ok(())
    }

    /// Spend scientific evaluations.
    ///
    /// # Errors
    /// Refuses over-consumption or stale authority.
    pub fn charge_evaluations(&mut self, amount: EvaluationUnits) -> Result<(), InvocationError> {
        self.ensure_active()?;
        let state = &self.owner.children[self.node];
        let remaining = match state.remaining.evaluations.0.checked_sub(amount.0) {
            Some(remaining) => remaining,
            None => {
                let available = state.remaining.evaluations.0;
                return Err(self.latch(exceeded(
                    "evaluations",
                    u128::from(amount.0),
                    u128::from(available),
                )));
            }
        };
        let direct = match state.direct_consumed.evaluations.0.checked_add(amount.0) {
            Some(direct) => direct,
            None => {
                return Err(self.latch(InvocationError::ArithmeticOverflow {
                    resource: "evaluations",
                }));
            }
        };
        let state = &mut self.owner.children[self.node];
        state.remaining.evaluations.0 = remaining;
        state.direct_consumed.evaluations.0 = direct;
        Ok(())
    }

    /// Check deadline, spend one poll, then observe cancellation in that fixed
    /// order.
    ///
    /// # Errors
    /// Refuses expired deadline, exhausted poll allowance, or cancellation.
    pub fn poll(&mut self, phase: &'static str) -> Result<(), InvocationError> {
        self.ensure_active()?;
        self.owner.observe_deadline(Some(self.node), phase)?;
        let state = &self.owner.children[self.node];
        let Some(remaining) = state.remaining.polls.0.checked_sub(1) else {
            return Err(self.latch(exceeded("polls", 1, 0)));
        };
        let Some(direct) = state.direct_consumed.polls.0.checked_add(1) else {
            return Err(self.latch(InvocationError::ArithmeticOverflow { resource: "polls" }));
        };
        {
            let state = &mut self.owner.children[self.node];
            state.remaining.polls.0 = remaining;
            state.direct_consumed.polls.0 = direct;
        }
        if self.owner.cancel_gate.is_requested() {
            return Err(self.latch(InvocationError::Cancelled { phase }));
        }
        Ok(())
    }

    /// Reserve concurrent memory through both the child sub-cap and the root
    /// operation-memory lease. The returned guard releases on drop/unwind.
    ///
    /// # Errors
    /// Refuses a child-cap or backing-lease overrun before allocation.
    pub fn reserve_memory<'child>(
        &'child mut self,
        what: &'static str,
        bytes: MemoryBytes,
    ) -> Result<InvocationMemoryReservation<'child, 'budget, 'clock>, InvocationError> {
        self.ensure_active()?;
        let state = &self.owner.children[self.node];
        let next = match state.memory_current.checked_add(bytes.0) {
            Some(next) => next,
            None => {
                return Err(self.latch(InvocationError::ArithmeticOverflow {
                    resource: "memory-bytes",
                }));
            }
        };
        let next_requested = match state.memory_requested.checked_add(u128::from(bytes.0)) {
            Some(next) => next,
            None => {
                return Err(self.latch(InvocationError::ArithmeticOverflow {
                    resource: "memory-requested",
                }));
            }
        };
        if next > state.remaining.memory.0 {
            let available = state.remaining.memory.0;
            return Err(self.latch(exceeded(
                "memory-bytes",
                u128::from(next),
                u128::from(available),
            )));
        }
        let mut ancestor = Some(self.node);
        while let Some(index) = ancestor {
            let state = &self.owner.children[index];
            if state.subtree_memory_current.checked_add(bytes.0).is_none() {
                return Err(self.latch(InvocationError::ArithmeticOverflow {
                    resource: "subtree-memory-bytes",
                }));
            }
            ancestor = state.parent;
        }
        let charge = match self.owner.backing_memory.reserve(what, bytes.0) {
            Ok(charge) => charge,
            Err(refusal) => {
                let error = InvocationError::MemoryRefused {
                    what: refusal.what,
                    requested: refusal.requested_bytes,
                    used: refusal.used_bytes,
                    limit: refusal.limit_bytes,
                };
                return Err(self.latch(error));
            }
        };
        let mut ancestor = Some(self.node);
        while let Some(index) = ancestor {
            let state = &mut self.owner.children[index];
            state.subtree_memory_current = state
                .subtree_memory_current
                .checked_add(bytes.0)
                .expect("subtree memory was preflighted");
            state.memory_peak = state.memory_peak.max(state.subtree_memory_current);
            ancestor = state.parent;
        }
        let state = &mut self.owner.children[self.node];
        state.memory_current = next;
        state.direct_memory_peak = state.direct_memory_peak.max(next);
        state.memory_requested = next_requested;
        Ok(InvocationMemoryReservation {
            child: self,
            bytes: bytes.0,
            _charge: charge,
        })
    }

    /// Permanently retain publication capacity.
    ///
    /// # Errors
    /// Refuses output overrun or stale authority.
    pub fn publish_output(&mut self, bytes: OutputBytes) -> Result<(), InvocationError> {
        self.ensure_active()?;
        self.owner
            .observe_terminal(Some(self.node), "child-publication")?;
        let state = &self.owner.children[self.node];
        let remaining = match state.remaining.output.0.checked_sub(bytes.0) {
            Some(remaining) => remaining,
            None => {
                let available = state.remaining.output.0;
                return Err(self.latch(exceeded(
                    "output-bytes",
                    u128::from(bytes.0),
                    u128::from(available),
                )));
            }
        };
        let retained = match state.output_retained.checked_add(bytes.0) {
            Some(retained) => retained,
            None => {
                return Err(self.latch(InvocationError::ArithmeticOverflow {
                    resource: "output-retained",
                }));
            }
        };
        let direct = match state.direct_consumed.output.0.checked_add(bytes.0) {
            Some(direct) => direct,
            None => {
                return Err(self.latch(InvocationError::ArithmeticOverflow {
                    resource: "output-bytes",
                }));
            }
        };
        let state = &mut self.owner.children[self.node];
        state.remaining.output.0 = remaining;
        state.output_retained = retained;
        state.direct_consumed.output.0 = direct;
        Ok(())
    }

    /// Latch a structured scientific refusal so terminal receipts cannot
    /// misrepresent a domain error as successful completion.
    pub fn refuse(&mut self, phase: &'static str, reason: ContentHash) -> InvocationError {
        self.latch(InvocationError::ExplicitRefusal { phase, reason })
    }

    /// Return unused authority exactly once and retain terminal disposition.
    ///
    /// # Errors
    /// Refuses live nested children or live memory reservations.
    pub fn finish(self) -> Result<InvocationDisposition, InvocationError> {
        if self.owner.children[self.node].finalization_required {
            return Err(InvocationError::FinalizationRequired);
        }
        let _ = self
            .owner
            .observe_terminal(Some(self.node), "child-finalize");
        self.owner.close_child(self.node)
    }
}

/// Scientific-work authority for a child whose cleanup resources were
/// isolated at admission. This type intentionally has no `finish` or
/// `publish_output` method: callers must consume it through
/// [`Self::begin_finalization`], and publication is available only from the
/// resulting [`ChildFinalizer`].
///
/// ```compile_fail
/// fn bypass_finalizer(child: fs_exec::FinalizableChildBudget<'_, '_>) {
///     let _ = child.finish();
/// }
/// ```
///
/// ```compile_fail
/// fn reissue_finalizer(child: fs_exec::FinalizableChildBudget<'_, '_>) {
///     let _first = child.begin_finalization();
///     let _second = child.begin_finalization();
/// }
/// ```
pub struct FinalizableChildBudget<'budget, 'clock> {
    child: ChildBudget<'budget, 'clock>,
}

impl core::fmt::Debug for FinalizableChildBudget<'_, '_> {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        formatter
            .debug_struct("FinalizableChildBudget")
            .field("id", &self.child.id())
            .field("remaining", &self.child.remaining())
            .finish_non_exhaustive()
    }
}

impl<'budget, 'clock> FinalizableChildBudget<'budget, 'clock> {
    /// Deterministic child identity.
    #[must_use]
    pub fn id(&self) -> ContentHash {
        self.child.id()
    }

    /// Remaining ordinary scientific resources. The cleanup reserve is not
    /// visible or spendable through this authority.
    #[must_use]
    pub fn remaining(&self) -> InvocationResources {
        self.child.remaining()
    }

    /// Split an ordinary nested child from scientific capacity.
    pub fn split_child<'child>(
        &'child mut self,
        phase: &'static str,
        grant: InvocationResources,
    ) -> Result<ChildBudget<'child, 'clock>, InvocationError> {
        self.child.split_child(phase, grant)
    }

    /// Split a nested child with its own isolated cleanup reserve.
    pub fn split_finalizable_child<'child>(
        &'child mut self,
        phase: &'static str,
        grant: InvocationResources,
        finalization: FinalizationResources,
    ) -> Result<FinalizableChildBudget<'child, 'clock>, InvocationError> {
        self.child
            .split_finalizable_child(phase, grant, finalization)
    }

    /// Spend ordinary scientific work.
    pub fn charge_work(&mut self, amount: WorkUnits) -> Result<(), InvocationError> {
        self.child.charge_work(amount)
    }

    /// Spend ordinary cost capacity.
    pub fn charge_cost(&mut self, amount: CostUnits) -> Result<(), InvocationError> {
        self.child.charge_cost(amount)
    }

    /// Spend scientific evaluation capacity.
    pub fn charge_evaluations(&mut self, amount: EvaluationUnits) -> Result<(), InvocationError> {
        self.child.charge_evaluations(amount)
    }

    /// Poll during scientific work.
    pub fn poll(&mut self, phase: &'static str) -> Result<(), InvocationError> {
        self.child.poll(phase)
    }

    /// Reserve already-admitted ordinary memory. The reservation must drop
    /// before this authority can transition into finalization.
    pub fn reserve_memory<'child>(
        &'child mut self,
        what: &'static str,
        bytes: MemoryBytes,
    ) -> Result<InvocationMemoryReservation<'child, 'budget, 'clock>, InvocationError> {
        self.child.reserve_memory(what, bytes)
    }

    /// Latch one structured scientific refusal.
    pub fn refuse(&mut self, phase: &'static str, reason: ContentHash) -> InvocationError {
        self.child.refuse(phase, reason)
    }

    /// Consume scientific authority and enter the only post-terminal state
    /// that can spend the isolated cleanup reserve.
    #[must_use]
    pub fn begin_finalization(self) -> ChildFinalizer<'budget, 'clock> {
        let ChildBudget { owner, node } = self.child;
        let state = &mut owner.children[node];
        debug_assert!(state.finalization_required);
        debug_assert!(!state.finalization_started);
        state.finalization_started = true;
        let _ = owner.observe_terminal(Some(node), "child-finalization-begin");
        ChildFinalizer { owner, node }
    }
}

impl InvocationPoll for FinalizableChildBudget<'_, '_> {
    fn invocation_poll(&mut self, phase: &'static str) -> Result<(), InvocationError> {
        self.poll(phase)
    }

    fn invocation_polls_remaining(&self) -> PollUnits {
        self.remaining().polls()
    }
}

/// Result of a cleanup poll. Terminal state is reported without disabling the
/// finalizer; cleanup authority remains usable until it seals its report.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FinalizationObservation {
    /// No cancellation, deadline, or earlier refusal is latched.
    Clear,
    /// The first terminal cause is latched; cleanup must continue and success
    /// publication is forbidden.
    Terminal(InvocationDisposition),
}

/// One-use token proving that a successful child spent its final terminal poll
/// and was still non-terminal at that boundary.
#[derive(Debug, PartialEq, Eq)]
pub struct PreparedPublication {
    child: ContentHash,
}

/// Failed atomic publication together with the still-unpublished staged value.
///
/// Returning the staged value makes the unchanged-destination guarantee
/// observable without requiring rollback code or a fallible caller callback.
#[derive(Debug)]
pub struct PublicationCommitError<T> {
    error: InvocationError,
    staged: T,
}

impl<T> PublicationCommitError<T> {
    /// Structured reason publication was refused.
    #[must_use]
    pub const fn error(&self) -> &InvocationError {
        &self.error
    }

    /// Recover both the refusal and the value that never reached the
    /// destination.
    #[must_use]
    pub fn into_parts(self) -> (InvocationError, T) {
        (self.error, self.staged)
    }
}

/// Affine terminal authority. It may spend only the isolated cleanup reserve,
/// abort or commit publication once, and seal the child. It cannot restart
/// scientific work, allocate memory, or mint nested children.
///
/// ```compile_fail
/// fn duplicate_finalizer(finalizer: fs_exec::ChildFinalizer<'_, '_>) {
///     let _second = finalizer.clone();
/// }
/// ```
///
/// ```compile_fail
/// fn spend_cleanup_as_science(finalizer: &mut fs_exec::ChildFinalizer<'_, '_>) {
///     finalizer.charge_work(fs_exec::WorkUnits::new(1));
/// }
/// ```
pub struct ChildFinalizer<'budget, 'clock> {
    owner: &'budget mut InvocationBudget<'clock>,
    node: usize,
}

impl core::fmt::Debug for ChildFinalizer<'_, '_> {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        let state = &self.owner.children[self.node];
        formatter
            .debug_struct("ChildFinalizer")
            .field("child", &state.id)
            .field("remaining", &state.finalization_remaining)
            .field("publication", &state.finalization_publication)
            .finish_non_exhaustive()
    }
}

impl ChildFinalizer<'_, '_> {
    fn ensure_cleanup_open(&self) -> Result<(), InvocationError> {
        match self.owner.children[self.node].finalization_publication {
            FinalizationPublication::Pending => Ok(()),
            FinalizationPublication::Prepared
            | FinalizationPublication::Aborted
            | FinalizationPublication::Committed { .. } => {
                Err(InvocationError::PublicationAlreadySealed)
            }
        }
    }

    fn ensure_local_resources_closed(&self) -> Result<(), InvocationError> {
        let state = &self.owner.children[self.node];
        if state.live_children != 0 {
            return Err(InvocationError::LiveNestedChildren {
                count: state.live_children,
            });
        }
        if state.memory_current != 0 || state.subtree_memory_current != 0 {
            return Err(InvocationError::LiveMemoryReservations {
                bytes: state.subtree_memory_current.max(state.memory_current),
            });
        }
        Ok(())
    }

    /// Exact child identity being finalized.
    #[must_use]
    pub fn child_id(&self) -> ContentHash {
        self.owner.children[self.node].id
    }

    /// Remaining cleanup-only capacity.
    #[must_use]
    pub fn remaining(&self) -> FinalizationResources {
        self.owner.children[self.node].finalization_remaining
    }

    /// Spend declared cleanup work even after cancellation/refusal is latched.
    ///
    /// # Errors
    /// Refuses only stale authority or cleanup-reserve exhaustion. It never
    /// clears or replaces the first scientific terminal cause.
    pub fn charge_cleanup_work(&mut self, amount: WorkUnits) -> Result<(), InvocationError> {
        if self.owner.children[self.node].disposition.is_some() {
            return Err(InvocationError::InactiveChild);
        }
        self.ensure_cleanup_open()?;
        let state = &self.owner.children[self.node];
        let remaining = match state
            .finalization_remaining
            .work
            .get()
            .checked_sub(amount.get())
        {
            Some(remaining) => remaining,
            None => {
                let error = InvocationError::ResourceExceeded {
                    resource: "finalization-work",
                    requested: amount.get(),
                    available: state.finalization_remaining.work.get(),
                };
                self.owner.latch_failure(Some(self.node), error.clone());
                return Err(error);
            }
        };
        let consumed = match state.direct_consumed.work.get().checked_add(amount.get()) {
            Some(consumed) => consumed,
            None => {
                let error = InvocationError::ArithmeticOverflow {
                    resource: "finalization-work",
                };
                self.owner.latch_failure(Some(self.node), error.clone());
                return Err(error);
            }
        };
        let state = &mut self.owner.children[self.node];
        state.finalization_remaining.work = WorkUnits::new(remaining);
        state.direct_consumed.work = WorkUnits::new(consumed);
        Ok(())
    }

    /// Spend one cleanup poll while preserving finalizer usability after a
    /// terminal observation.
    ///
    /// Deadline and cancellation are observed only when no earlier failure is
    /// latched, so cleanup can never overwrite the first terminal cause.
    pub fn poll_cleanup(
        &mut self,
        phase: &'static str,
    ) -> Result<FinalizationObservation, InvocationError> {
        if self.owner.children[self.node].disposition.is_some() {
            return Err(InvocationError::InactiveChild);
        }
        self.ensure_cleanup_open()?;
        if self.owner.current_failure(Some(self.node)).is_none() {
            if let Some(deadline) = self.owner.limits.deadline {
                let now = self.owner.clock.now();
                self.owner.last_deadline_observation = Some(now);
                if now >= deadline {
                    self.owner.cancel_gate.request();
                    self.owner.latch_failure(
                        Some(self.node),
                        InvocationError::DeadlineExpired {
                            phase,
                            deadline_ns: deadline.as_nanos(),
                            observed_ns: now.as_nanos(),
                        },
                    );
                }
            }
        }
        let state = &self.owner.children[self.node];
        let remaining = match state.finalization_remaining.polls.get().checked_sub(1) {
            Some(remaining) => remaining,
            None => {
                let error = InvocationError::ResourceExceeded {
                    resource: "finalization-polls",
                    requested: 1,
                    available: 0,
                };
                self.owner.latch_failure(Some(self.node), error.clone());
                return Err(error);
            }
        };
        let consumed = match state.direct_consumed.polls.get().checked_add(1) {
            Some(consumed) => consumed,
            None => {
                let error = InvocationError::ArithmeticOverflow {
                    resource: "finalization-polls",
                };
                self.owner.latch_failure(Some(self.node), error.clone());
                return Err(error);
            }
        };
        {
            let state = &mut self.owner.children[self.node];
            state.finalization_remaining.polls = PollUnits::new(remaining);
            state.direct_consumed.polls = PollUnits::new(consumed);
        }
        if self.owner.current_failure(Some(self.node)).is_none()
            && self.owner.cancel_gate.is_requested()
        {
            self.owner
                .latch_failure(Some(self.node), InvocationError::Cancelled { phase });
        }
        Ok(self
            .owner
            .current_failure(Some(self.node))
            .as_ref()
            .map_or(FinalizationObservation::Clear, |failure| {
                FinalizationObservation::Terminal(error_disposition(failure))
            }))
    }

    /// Spend the mandatory final pre-publication poll and mint a one-use token
    /// only while the child is still successful.
    pub fn prepare_publication(&mut self) -> Result<PreparedPublication, InvocationError> {
        match self.owner.children[self.node].finalization_publication {
            FinalizationPublication::Pending => {}
            FinalizationPublication::Prepared
            | FinalizationPublication::Aborted
            | FinalizationPublication::Committed { .. } => {
                return Err(InvocationError::PublicationAlreadySealed);
            }
        }
        self.ensure_local_resources_closed()?;
        if matches!(
            self.poll_cleanup("child-finalization-pre-publication")?,
            FinalizationObservation::Terminal(_)
        ) {
            return Err(InvocationError::PublicationForbidden);
        }
        self.owner.children[self.node].finalization_publication = FinalizationPublication::Prepared;
        Ok(PreparedPublication {
            child: self.owner.children[self.node].id,
        })
    }

    /// Commit retained output with explicitly child-local durability.
    ///
    /// Deadline/cancellation is checked before and immediately after the
    /// infallible swap. A request that wins either check converts the prepared
    /// state to `Aborted`; the second check restores the old destination and
    /// returns the staged value. The successful post-swap check is the
    /// publication linearization point.
    ///
    /// `declared_bytes` is logical capacity accounting supplied by the caller;
    /// this generic layer does not inspect `T` and makes no claim that it
    /// equals an allocation size, serialized length, or content identity.
    ///
    /// A successful return does **not** seal the enclosing invocation. Later
    /// root cancellation/refusal may therefore coexist with this durable
    /// child output. Consumers requiring unchanged destination on any
    /// invocation failure must use the root invocation-atomic protocol.
    pub fn commit_child_local_publication<T>(
        &mut self,
        prepared: PreparedPublication,
        declared_bytes: OutputBytes,
        destination: &mut T,
        staged: T,
    ) -> Result<T, PublicationCommitError<T>> {
        self.commit_publication_inner(prepared, declared_bytes, destination, staged, || {})
    }

    /// Compatibility spelling for [`Self::commit_child_local_publication`].
    ///
    /// This method has child-local, not invocation-atomic, semantics.
    pub fn commit_publication<T>(
        &mut self,
        prepared: PreparedPublication,
        declared_bytes: OutputBytes,
        destination: &mut T,
        staged: T,
    ) -> Result<T, PublicationCommitError<T>> {
        self.commit_child_local_publication(prepared, declared_bytes, destination, staged)
    }

    fn commit_publication_inner<T>(
        &mut self,
        prepared: PreparedPublication,
        declared_bytes: OutputBytes,
        destination: &mut T,
        staged: T,
        after_swap: impl FnOnce(),
    ) -> Result<T, PublicationCommitError<T>> {
        if prepared.child != self.owner.children[self.node].id {
            return Err(PublicationCommitError {
                error: InvocationError::FinalizationReceiptMismatch {
                    invariant: "prepared-child",
                },
                staged,
            });
        }
        if self.owner.children[self.node].finalization_publication
            != FinalizationPublication::Prepared
        {
            return Err(PublicationCommitError {
                error: InvocationError::PublicationNotPrepared,
                staged,
            });
        }
        if let Err(error) = self.ensure_local_resources_closed() {
            // A recovered prepared finalizer can observe abandoned descendants
            // that were inaccessible while the original borrow was live.
            // Re-open preparation so the caller can drain/release and either
            // spend another final poll or abort without stranding authority.
            self.owner.children[self.node].finalization_publication =
                FinalizationPublication::Pending;
            return Err(PublicationCommitError { error, staged });
        }
        if self.owner.current_failure(Some(self.node)).is_none() {
            let _ = self
                .owner
                .observe_terminal(Some(self.node), "child-finalization-commit");
        }
        if self.owner.current_failure(Some(self.node)).is_some() {
            self.owner.children[self.node].finalization_publication_scope =
                Some(InvocationPublicationScope::ChildLocal);
            self.owner.children[self.node].finalization_publication =
                FinalizationPublication::Aborted;
            return Err(PublicationCommitError {
                error: InvocationError::PublicationForbidden,
                staged,
            });
        }
        let state = &self.owner.children[self.node];
        let remaining = match state
            .remaining
            .output
            .get()
            .checked_sub(declared_bytes.get())
        {
            Some(remaining) => remaining,
            None => {
                let error = InvocationError::ResourceExceeded {
                    resource: "output-bytes",
                    requested: u128::from(declared_bytes.get()),
                    available: u128::from(state.remaining.output.get()),
                };
                self.owner.latch_failure(Some(self.node), error.clone());
                self.owner.children[self.node].finalization_publication_scope =
                    Some(InvocationPublicationScope::ChildLocal);
                self.owner.children[self.node].finalization_publication =
                    FinalizationPublication::Aborted;
                return Err(PublicationCommitError { error, staged });
            }
        };
        let retained = match state.output_retained.checked_add(declared_bytes.get()) {
            Some(retained) => retained,
            None => {
                let error = InvocationError::ArithmeticOverflow {
                    resource: "output-retained",
                };
                self.owner.latch_failure(Some(self.node), error.clone());
                self.owner.children[self.node].finalization_publication_scope =
                    Some(InvocationPublicationScope::ChildLocal);
                self.owner.children[self.node].finalization_publication =
                    FinalizationPublication::Aborted;
                return Err(PublicationCommitError { error, staged });
            }
        };
        let consumed = match state
            .direct_consumed
            .output
            .get()
            .checked_add(declared_bytes.get())
        {
            Some(consumed) => consumed,
            None => {
                let error = InvocationError::ArithmeticOverflow {
                    resource: "output-bytes",
                };
                self.owner.latch_failure(Some(self.node), error.clone());
                self.owner.children[self.node].finalization_publication_scope =
                    Some(InvocationPublicationScope::ChildLocal);
                self.owner.children[self.node].finalization_publication =
                    FinalizationPublication::Aborted;
                return Err(PublicationCommitError { error, staged });
            }
        };
        // `mem::replace` is the only real destination mutation supported by
        // this protocol. It is infallible and runs after every refusal check,
        // so accounting cannot say Committed without the destination changing,
        // and a refusal always returns `staged` with the destination untouched.
        let replaced = core::mem::replace(destination, staged);
        after_swap();
        if self.owner.current_failure(Some(self.node)).is_none() {
            let _ = self
                .owner
                .observe_terminal(Some(self.node), "child-finalization-commit");
        }
        if self.owner.current_failure(Some(self.node)).is_some() {
            let unpublished = core::mem::replace(destination, replaced);
            self.owner.children[self.node].finalization_publication_scope =
                Some(InvocationPublicationScope::ChildLocal);
            self.owner.children[self.node].finalization_publication =
                FinalizationPublication::Aborted;
            return Err(PublicationCommitError {
                error: InvocationError::PublicationForbidden,
                staged: unpublished,
            });
        }
        let state = &mut self.owner.children[self.node];
        state.remaining.output = OutputBytes::new(remaining);
        state.output_retained = retained;
        state.direct_consumed.output = OutputBytes::new(consumed);
        state.finalization_publication_scope = Some(InvocationPublicationScope::ChildLocal);
        state.finalization_publication = FinalizationPublication::Committed {
            bytes: declared_bytes,
        };
        Ok(replaced)
    }

    /// Seal a child-local unchanged-destination decision for cancellation,
    /// refusal, or a successful no-output operation.
    pub fn abort_child_local_publication(&mut self) -> Result<(), InvocationError> {
        match self.owner.children[self.node].finalization_publication {
            FinalizationPublication::Pending | FinalizationPublication::Prepared => {
                self.owner.children[self.node].finalization_publication_scope =
                    Some(InvocationPublicationScope::ChildLocal);
                self.owner.children[self.node].finalization_publication =
                    FinalizationPublication::Aborted;
                Ok(())
            }
            FinalizationPublication::Aborted | FinalizationPublication::Committed { .. } => {
                Err(InvocationError::PublicationAlreadySealed)
            }
        }
    }

    /// Compatibility spelling for [`Self::abort_child_local_publication`].
    ///
    /// This method seals child-local, not invocation-atomic, evidence.
    pub fn abort_publication(&mut self) -> Result<(), InvocationError> {
        self.abort_child_local_publication()
    }

    /// Return unused cleanup capacity, close the child exactly once, and mint
    /// immutable finalization evidence.
    ///
    /// # Errors
    /// Refuses live nested authority/memory or an unsealed publication state.
    pub fn finish(&mut self) -> Result<FinalizationReport, InvocationError> {
        if self.owner.children[self.node].disposition.is_some() {
            let report = finalization_report_from_state(&self.owner.children[self.node])?;
            finalization_report_semantics(&report)?;
            if self.owner.children[self.node].finalization_report_root != Some(report.root()) {
                return Err(InvocationError::FinalizationReceiptMismatch {
                    invariant: "stored-finalization-root",
                });
            }
            return Ok(report);
        }
        self.ensure_local_resources_closed()?;
        let state = &self.owner.children[self.node];
        if matches!(
            state.finalization_publication,
            FinalizationPublication::Pending | FinalizationPublication::Prepared
        ) {
            return Err(InvocationError::FinalizationIncomplete {
                step: "publication-seal",
            });
        }
        if state.finalization_publication == FinalizationPublication::Aborted {
            let _ = self
                .owner
                .observe_terminal(Some(self.node), "child-finalization-seal");
        }
        let state = &self.owner.children[self.node];
        let disposition = state
            .failure
            .as_ref()
            .map_or(InvocationDisposition::Completed, error_disposition);
        let report = prospective_finalization_report(state, disposition)?;
        finalization_report_semantics(&report)?;
        let closed = self.owner.close_child(self.node)?;
        debug_assert_eq!(closed, disposition);
        self.owner.children[self.node].finalization_report_root = Some(report.root());
        Ok(report)
    }
}

/// Small object-safe poll seam for lower-layer progress engines.
pub trait InvocationPoll {
    /// Observe deadline/cancellation while consuming one affine poll.
    fn invocation_poll(&mut self, phase: &'static str) -> Result<(), InvocationError>;

    /// Remaining poll opportunities.
    fn invocation_polls_remaining(&self) -> PollUnits;
}

impl InvocationPoll for ChildBudget<'_, '_> {
    fn invocation_poll(&mut self, phase: &'static str) -> Result<(), InvocationError> {
        self.poll(phase)
    }

    fn invocation_polls_remaining(&self) -> PollUnits {
        self.remaining().polls()
    }
}

/// RAII memory reservation. Scientific code continues spending through
/// [`Self::budget`] while the allocation charge remains live.
pub struct InvocationMemoryReservation<'child, 'budget, 'clock> {
    child: &'child mut ChildBudget<'budget, 'clock>,
    bytes: u64,
    _charge: LeaseCharge,
}

impl<'budget, 'clock> InvocationMemoryReservation<'_, 'budget, 'clock> {
    /// Continue using the same child authority while this memory is live.
    pub fn budget(&mut self) -> &mut ChildBudget<'budget, 'clock> {
        self.child
    }

    /// Reserved bytes.
    #[must_use]
    pub const fn bytes(&self) -> MemoryBytes {
        MemoryBytes(self.bytes)
    }
}

impl Drop for InvocationMemoryReservation<'_, '_, '_> {
    fn drop(&mut self) {
        let node = self.child.node;
        let mut violation = false;
        {
            let state = &mut self.child.owner.children[node];
            match (
                state.memory_current.checked_sub(self.bytes),
                state.memory_released.checked_add(u128::from(self.bytes)),
            ) {
                (Some(current), Some(released)) => {
                    state.memory_current = current;
                    state.memory_released = released;
                }
                _ => {
                    state.memory_current = u64::MAX;
                    state.memory_released = u128::MAX;
                    violation = true;
                }
            }
        }
        let mut ancestor = Some(node);
        while let Some(index) = ancestor {
            let state = &mut self.child.owner.children[index];
            match state.subtree_memory_current.checked_sub(self.bytes) {
                Some(current) => state.subtree_memory_current = current,
                None => {
                    state.subtree_memory_current = u64::MAX;
                    violation = true;
                }
            }
            ancestor = state.parent;
        }
        if violation {
            self.child
                .owner
                .latch_failure(Some(node), InvocationError::MemoryReleaseInvariant);
        }
    }
}

fn try_evidence_vec<T>(
    what: &'static str,
    requested_items: usize,
) -> Result<Vec<T>, InvocationError> {
    let requested_items_u64 = u64::try_from(requested_items).unwrap_or(u64::MAX);
    let mut values = Vec::new();
    values.try_reserve_exact(requested_items).map_err(|_| {
        InvocationError::EvidenceAllocationRefused {
            what,
            requested_items: requested_items_u64,
        }
    })?;
    Ok(values)
}

fn try_clone_invocation_receipt(
    receipt: &InvocationReceipt,
    what: &'static str,
) -> Result<InvocationReceipt, InvocationError> {
    let mut children = try_evidence_vec(what, receipt.children.len())?;
    children.extend(receipt.children.iter().cloned());
    Ok(InvocationReceipt {
        version: receipt.version,
        invocation_id: receipt.invocation_id,
        plan_binding: receipt.plan_binding,
        limits: receipt.limits.clone(),
        required: receipt.required,
        remaining: receipt.remaining,
        children,
        last_deadline_observation: receipt.last_deadline_observation,
        memory_peak: receipt.memory_peak,
        memory_requested: receipt.memory_requested,
        memory_released: receipt.memory_released,
        memory_refusals: receipt.memory_refusals,
        memory_first_refusal: receipt.memory_first_refusal.clone(),
        output_retained: receipt.output_retained,
        failure: receipt.failure.clone(),
        failure_origin: receipt.failure_origin,
        disposition: receipt.disposition,
        root: receipt.root,
    })
}

fn child_receipt(
    states: &[ChildState],
    state: &ChildState,
) -> Result<ChildReceipt, InvocationError> {
    let returned = state
        .remaining
        .checked_add(state.finalization_remaining.as_invocation_resources())?;
    let consumed = state.granted.checked_sub(returned)?;
    let parent = state.parent.map(|index| states[index].id);
    let finalization = if state.finalization_required {
        let report_root =
            state
                .finalization_report_root
                .ok_or(InvocationError::FinalizationIncomplete {
                    step: "report-commitment",
                })?;
        let publication_scope = state.finalization_publication_scope.ok_or(
            InvocationError::FinalizationIncomplete {
                step: "publication-scope",
            },
        )?;
        let finalizer_consumed = state
            .finalization_granted
            .checked_sub(state.finalization_remaining)?;
        let finalizer_consumed_resources = finalizer_consumed.as_invocation_resources();
        Some(ChildFinalizationEvidence {
            scientific_granted: state
                .granted
                .checked_sub(state.finalization_granted.as_invocation_resources())?,
            scientific_direct_consumed: state
                .direct_consumed
                .checked_sub(finalizer_consumed_resources)?,
            scientific_returned: state.remaining,
            granted: state.finalization_granted,
            consumed: finalizer_consumed,
            returned: state.finalization_remaining,
            publication_scope,
            publication: state.finalization_publication,
            report_root,
        })
    } else {
        None
    };
    let mut receipt = ChildReceipt {
        id: state.id,
        parent,
        ordinal: state.ordinal,
        phase: state.phase,
        granted: state.granted,
        consumed,
        direct_consumed: state.direct_consumed,
        returned,
        direct_memory_peak: state.direct_memory_peak,
        memory_peak: state.memory_peak,
        memory_requested: state.memory_requested,
        memory_released: state.memory_released,
        output_retained: state.output_retained,
        finalization,
        failure: state.failure.clone(),
        failure_inherited: state.failure_inherited,
        disposition: state.disposition.ok_or(InvocationError::InactiveChild)?,
        root: ContentHash([0; 32]),
    };
    receipt.root = child_receipt_root(&receipt);
    Ok(receipt)
}

fn child_id(
    invocation: ContentHash,
    plan_binding: Option<InvocationPlanBinding>,
    parent: Option<ContentHash>,
    ordinal: u64,
    phase: &str,
    scientific_grant: InvocationResources,
    finalization_grant: FinalizationResources,
) -> ContentHash {
    let mut hasher = DomainHasher::new(CHILD_ID_DOMAIN);
    hash_field(&mut hasher, "invocation", invocation.as_bytes());
    hash_plan_binding(&mut hasher, plan_binding);
    hash_field(&mut hasher, "parent-present", &[u8::from(parent.is_some())]);
    if let Some(parent) = parent {
        hash_field(&mut hasher, "parent", parent.as_bytes());
    }
    hash_field(&mut hasher, "ordinal", &ordinal.to_le_bytes());
    hash_field(&mut hasher, "phase", phase.as_bytes());
    hash_resources(
        &mut hasher,
        [
            "scientific-grant.work",
            "scientific-grant.polls",
            "scientific-grant.cost",
            "scientific-grant.evaluations",
            "scientific-grant.memory",
            "scientific-grant.output",
        ],
        scientific_grant,
    );
    hash_finalization_resources(
        &mut hasher,
        ["finalization-grant.work", "finalization-grant.polls"],
        finalization_grant,
    );
    hasher.finalize()
}

fn hash_plan_binding(hasher: &mut DomainHasher, binding: Option<InvocationPlanBinding>) {
    hash_field(
        hasher,
        "plan-binding-present",
        &[u8::from(binding.is_some())],
    );
    if let Some(binding) = binding {
        hash_field(
            hasher,
            "plan-binding.schema-root",
            binding.schema_root.as_bytes(),
        );
        hash_field(
            hasher,
            "plan-binding.schema-version",
            &binding.schema_version.to_le_bytes(),
        );
        hash_field(
            hasher,
            "plan-binding.plan-root",
            binding.plan_root.as_bytes(),
        );
    }
}

fn hash_field(hasher: &mut DomainHasher, label: &str, value: &[u8]) {
    hasher.update(&(label.len() as u64).to_le_bytes());
    hasher.update(label.as_bytes());
    hasher.update(&(value.len() as u64).to_le_bytes());
    hasher.update(value);
}

fn hash_resources(hasher: &mut DomainHasher, labels: [&str; 6], resources: InvocationResources) {
    hash_field(hasher, labels[0], &resources.work.0.to_le_bytes());
    hash_field(hasher, labels[1], &resources.polls.0.to_le_bytes());
    hash_field(hasher, labels[2], &resources.cost.0.to_le_bytes());
    hash_field(hasher, labels[3], &resources.evaluations.0.to_le_bytes());
    hash_field(hasher, labels[4], &resources.memory.0.to_le_bytes());
    hash_field(hasher, labels[5], &resources.output.0.to_le_bytes());
}

fn hash_finalization_resources(
    hasher: &mut DomainHasher,
    labels: [&str; 2],
    resources: FinalizationResources,
) {
    hash_field(hasher, labels[0], &resources.work.get().to_le_bytes());
    hash_field(hasher, labels[1], &resources.polls.get().to_le_bytes());
}

fn hash_finalization_publication(hasher: &mut DomainHasher, publication: FinalizationPublication) {
    let (tag, retained) = match publication {
        FinalizationPublication::Pending => (0_u8, None),
        FinalizationPublication::Prepared => (1_u8, None),
        FinalizationPublication::Aborted => (2_u8, None),
        FinalizationPublication::Committed { bytes } => (3_u8, Some(bytes.get())),
    };
    hash_field(hasher, "publication.tag", &[tag]);
    hash_field(
        hasher,
        "publication.bytes-present",
        &[u8::from(retained.is_some())],
    );
    if let Some(retained) = retained {
        hash_field(hasher, "publication.bytes", &retained.to_le_bytes());
    }
}

fn hash_publication_scope(hasher: &mut DomainHasher, scope: InvocationPublicationScope) {
    let tag = match scope {
        InvocationPublicationScope::ChildLocal => 0_u8,
        InvocationPublicationScope::InvocationAtomic => 1_u8,
    };
    hash_field(hasher, "publication-scope", &[tag]);
}

fn encode_disposition(disposition: InvocationDisposition) -> u8 {
    match disposition {
        InvocationDisposition::Completed => 0,
        InvocationDisposition::Cancelled => 1,
        InvocationDisposition::Refused => 2,
    }
}

#[allow(clippy::too_many_lines)]
fn invocation_error_root(error: &InvocationError) -> ContentHash {
    let mut hasher = DomainHasher::new(INVOCATION_ERROR_DOMAIN);
    let tag = match error {
        InvocationError::ResourceExceeded { .. } => 0,
        InvocationError::ArithmeticOverflow { .. } => 1,
        InvocationError::DeadlineExpired { .. } => 2,
        InvocationError::Cancelled { .. } => 3,
        InvocationError::MemoryRefused { .. } => 4,
        InvocationError::ExplicitRefusal { .. } => 5,
        InvocationError::InactiveChild => 6,
        InvocationError::LiveNestedChildren { .. } => 7,
        InvocationError::LiveMemoryReservations { .. } => 8,
        InvocationError::UnfinishedChild { .. } => 9,
        InvocationError::MemoryReleaseInvariant => 10,
        InvocationError::EmptyPhase => 11,
        InvocationError::FinalizationRequired => 12,
        InvocationError::TransactionalOutputScopeViolation { .. } => 13,
        InvocationError::PublicationNotPrepared => 14,
        InvocationError::PublicationAlreadySealed => 15,
        InvocationError::PublicationForbidden => 16,
        InvocationError::FinalizationIncomplete { .. } => 17,
        InvocationError::FinalizationReceiptMismatch { .. } => 18,
        InvocationError::InvocationAlreadyFinalized => 19,
        InvocationError::EvidenceAllocationRefused { .. } => 20,
    };
    hash_field(&mut hasher, "tag", &[tag]);
    match error {
        InvocationError::ResourceExceeded {
            resource,
            requested,
            available,
        } => {
            hash_field(&mut hasher, "resource", resource.as_bytes());
            hash_field(&mut hasher, "requested", &requested.to_le_bytes());
            hash_field(&mut hasher, "available", &available.to_le_bytes());
        }
        InvocationError::ArithmeticOverflow { resource } => {
            hash_field(&mut hasher, "resource", resource.as_bytes());
        }
        InvocationError::DeadlineExpired {
            phase,
            deadline_ns,
            observed_ns,
        } => {
            hash_field(&mut hasher, "phase", phase.as_bytes());
            hash_field(&mut hasher, "deadline-nanos", &deadline_ns.to_le_bytes());
            hash_field(&mut hasher, "observed-nanos", &observed_ns.to_le_bytes());
        }
        InvocationError::Cancelled { phase } => {
            hash_field(&mut hasher, "phase", phase.as_bytes());
        }
        InvocationError::MemoryRefused {
            what,
            requested,
            used,
            limit,
        } => {
            hash_field(&mut hasher, "what", what.as_bytes());
            hash_field(&mut hasher, "requested", &requested.to_le_bytes());
            hash_field(&mut hasher, "used", &used.to_le_bytes());
            hash_field(&mut hasher, "limit", &limit.to_le_bytes());
        }
        InvocationError::ExplicitRefusal { phase, reason } => {
            hash_field(&mut hasher, "phase", phase.as_bytes());
            hash_field(&mut hasher, "reason", reason.as_bytes());
        }
        InvocationError::TransactionalOutputScopeViolation {
            ancestor,
            phase,
            requested,
        } => {
            hash_field(&mut hasher, "ancestor", ancestor.as_bytes());
            hash_field(&mut hasher, "phase", phase.as_bytes());
            hash_field(&mut hasher, "requested", &requested.to_le_bytes());
        }
        InvocationError::LiveNestedChildren { count } => {
            hash_field(&mut hasher, "count", &count.to_le_bytes());
        }
        InvocationError::LiveMemoryReservations { bytes } => {
            hash_field(&mut hasher, "bytes", &bytes.to_le_bytes());
        }
        InvocationError::UnfinishedChild { child } => {
            hash_field(&mut hasher, "child", child.as_bytes());
        }
        InvocationError::FinalizationIncomplete { step } => {
            hash_field(&mut hasher, "step", step.as_bytes());
        }
        InvocationError::FinalizationReceiptMismatch { invariant } => {
            hash_field(&mut hasher, "invariant", invariant.as_bytes());
        }
        InvocationError::EvidenceAllocationRefused {
            what,
            requested_items,
        } => {
            hash_field(&mut hasher, "what", what.as_bytes());
            hash_field(
                &mut hasher,
                "requested-items",
                &requested_items.to_le_bytes(),
            );
        }
        InvocationError::EmptyPhase
        | InvocationError::InactiveChild
        | InvocationError::FinalizationRequired
        | InvocationError::PublicationNotPrepared
        | InvocationError::PublicationAlreadySealed
        | InvocationError::PublicationForbidden
        | InvocationError::InvocationAlreadyFinalized
        | InvocationError::MemoryReleaseInvariant => {}
    }
    hasher.finalize()
}

fn child_receipt_root(receipt: &ChildReceipt) -> ContentHash {
    let mut hasher = DomainHasher::new(CHILD_RECEIPT_DOMAIN);
    hash_field(&mut hasher, "id", receipt.id.as_bytes());
    hash_field(
        &mut hasher,
        "parent-present",
        &[u8::from(receipt.parent.is_some())],
    );
    if let Some(parent) = receipt.parent {
        hash_field(&mut hasher, "parent", parent.as_bytes());
    }
    hash_field(&mut hasher, "ordinal", &receipt.ordinal.to_le_bytes());
    hash_field(&mut hasher, "phase", receipt.phase.as_bytes());
    hash_resources(
        &mut hasher,
        [
            "granted.work",
            "granted.polls",
            "granted.cost",
            "granted.evaluations",
            "granted.memory",
            "granted.output",
        ],
        receipt.granted,
    );
    hash_resources(
        &mut hasher,
        [
            "consumed.work",
            "consumed.polls",
            "consumed.cost",
            "consumed.evaluations",
            "consumed.memory",
            "consumed.output",
        ],
        receipt.consumed,
    );
    hash_resources(
        &mut hasher,
        [
            "direct-consumed.work",
            "direct-consumed.polls",
            "direct-consumed.cost",
            "direct-consumed.evaluations",
            "direct-consumed.memory",
            "direct-consumed.output",
        ],
        receipt.direct_consumed,
    );
    hash_resources(
        &mut hasher,
        [
            "returned.work",
            "returned.polls",
            "returned.cost",
            "returned.evaluations",
            "returned.memory",
            "returned.output",
        ],
        receipt.returned,
    );
    hash_field(
        &mut hasher,
        "direct-memory-peak",
        &receipt.direct_memory_peak.to_le_bytes(),
    );
    hash_field(
        &mut hasher,
        "memory-peak",
        &receipt.memory_peak.to_le_bytes(),
    );
    hash_field(
        &mut hasher,
        "memory-requested",
        &receipt.memory_requested.to_le_bytes(),
    );
    hash_field(
        &mut hasher,
        "memory-released",
        &receipt.memory_released.to_le_bytes(),
    );
    hash_field(
        &mut hasher,
        "output-retained",
        &receipt.output_retained.to_le_bytes(),
    );
    hash_field(
        &mut hasher,
        "finalization-present",
        &[u8::from(receipt.finalization.is_some())],
    );
    if let Some(finalization) = &receipt.finalization {
        hash_resources(
            &mut hasher,
            [
                "finalization.scientific-granted.work",
                "finalization.scientific-granted.polls",
                "finalization.scientific-granted.cost",
                "finalization.scientific-granted.evaluations",
                "finalization.scientific-granted.memory",
                "finalization.scientific-granted.output",
            ],
            finalization.scientific_granted,
        );
        hash_resources(
            &mut hasher,
            [
                "finalization.scientific-direct-consumed.work",
                "finalization.scientific-direct-consumed.polls",
                "finalization.scientific-direct-consumed.cost",
                "finalization.scientific-direct-consumed.evaluations",
                "finalization.scientific-direct-consumed.memory",
                "finalization.scientific-direct-consumed.output",
            ],
            finalization.scientific_direct_consumed,
        );
        hash_resources(
            &mut hasher,
            [
                "finalization.scientific-returned.work",
                "finalization.scientific-returned.polls",
                "finalization.scientific-returned.cost",
                "finalization.scientific-returned.evaluations",
                "finalization.scientific-returned.memory",
                "finalization.scientific-returned.output",
            ],
            finalization.scientific_returned,
        );
        hash_finalization_resources(
            &mut hasher,
            ["finalization.granted.work", "finalization.granted.polls"],
            finalization.granted,
        );
        hash_finalization_resources(
            &mut hasher,
            ["finalization.consumed.work", "finalization.consumed.polls"],
            finalization.consumed,
        );
        hash_finalization_resources(
            &mut hasher,
            ["finalization.returned.work", "finalization.returned.polls"],
            finalization.returned,
        );
        hash_publication_scope(&mut hasher, finalization.publication_scope);
        hash_finalization_publication(&mut hasher, finalization.publication);
        hash_field(
            &mut hasher,
            "finalization.report-root",
            finalization.report_root.as_bytes(),
        );
    }
    hash_field(
        &mut hasher,
        "failure-inherited",
        &[u8::from(receipt.failure_inherited)],
    );
    hash_field(
        &mut hasher,
        "failure-present",
        &[u8::from(receipt.failure.is_some())],
    );
    if let Some(failure) = &receipt.failure {
        hash_field(
            &mut hasher,
            "failure-root",
            invocation_error_root(failure).as_bytes(),
        );
    }
    hash_field(
        &mut hasher,
        "disposition",
        &[encode_disposition(receipt.disposition)],
    );
    hasher.finalize()
}

fn finalization_report_root(report: &FinalizationReport) -> ContentHash {
    let mut hasher = DomainHasher::new(FINALIZATION_REPORT_DOMAIN);
    hash_field(&mut hasher, "version", &report.version.to_le_bytes());
    hash_plan_binding(&mut hasher, report.plan_binding);
    hash_field(&mut hasher, "child", report.child.as_bytes());
    hash_finalization_resources(
        &mut hasher,
        ["granted.work", "granted.polls"],
        report.granted,
    );
    hash_finalization_resources(
        &mut hasher,
        ["consumed.work", "consumed.polls"],
        report.consumed,
    );
    hash_finalization_resources(
        &mut hasher,
        ["returned.work", "returned.polls"],
        report.returned,
    );
    hash_publication_scope(&mut hasher, report.publication_scope);
    hash_finalization_publication(&mut hasher, report.publication);
    hash_field(
        &mut hasher,
        "failure-present",
        &[u8::from(report.failure.is_some())],
    );
    if let Some(failure) = &report.failure {
        hash_field(
            &mut hasher,
            "failure-root",
            invocation_error_root(failure).as_bytes(),
        );
    }
    hash_field(
        &mut hasher,
        "disposition",
        &[encode_disposition(report.disposition)],
    );
    hasher.finalize()
}

fn finalization_report_from_state(
    state: &ChildState,
) -> Result<FinalizationReport, InvocationError> {
    let disposition = state
        .disposition
        .ok_or(InvocationError::FinalizationIncomplete {
            step: "child-close",
        })?;
    prospective_finalization_report(state, disposition)
}

fn prospective_finalization_report(
    state: &ChildState,
    disposition: InvocationDisposition,
) -> Result<FinalizationReport, InvocationError> {
    if !state.finalization_required || !state.finalization_started {
        return Err(InvocationError::FinalizationIncomplete {
            step: "finalizer-start",
        });
    }
    if matches!(
        state.finalization_publication,
        FinalizationPublication::Pending | FinalizationPublication::Prepared
    ) {
        return Err(InvocationError::FinalizationIncomplete {
            step: "publication-seal",
        });
    }
    let granted = state.finalization_granted;
    let returned = state.finalization_remaining;
    let consumed = granted.checked_sub(returned)?;
    let publication_scope =
        state
            .finalization_publication_scope
            .ok_or(InvocationError::FinalizationIncomplete {
                step: "publication-scope",
            })?;
    let mut report = FinalizationReport {
        version: FINALIZATION_REPORT_VERSION,
        plan_binding: state.plan_binding,
        child: state.id,
        granted,
        consumed,
        returned,
        publication_scope,
        publication: state.finalization_publication,
        failure: state.failure.clone(),
        disposition,
        root: ContentHash([0; 32]),
    };
    report.root = finalization_report_root(&report);
    Ok(report)
}

fn finalization_report_semantics(report: &FinalizationReport) -> Result<(), InvocationError> {
    if report.version != FINALIZATION_REPORT_VERSION {
        return Err(InvocationError::FinalizationReceiptMismatch {
            invariant: "version",
        });
    }
    if report.root != finalization_report_root(report) {
        return Err(InvocationError::FinalizationReceiptMismatch {
            invariant: "report-root",
        });
    }
    let consumed = report.granted.checked_sub(report.returned).map_err(|_| {
        InvocationError::FinalizationReceiptMismatch {
            invariant: "resource-conservation",
        }
    })?;
    if consumed != report.consumed {
        return Err(InvocationError::FinalizationReceiptMismatch {
            invariant: "consumed-definition",
        });
    }
    if matches!(
        report.publication,
        FinalizationPublication::Pending | FinalizationPublication::Prepared
    ) {
        return Err(InvocationError::FinalizationReceiptMismatch {
            invariant: "publication-terminal",
        });
    }
    let expected = report
        .failure
        .as_ref()
        .map_or(InvocationDisposition::Completed, error_disposition);
    if report.disposition != expected {
        return Err(InvocationError::FinalizationReceiptMismatch {
            invariant: "disposition",
        });
    }
    if report.failure.is_some()
        && matches!(
            report.publication,
            FinalizationPublication::Committed { .. }
        )
    {
        return Err(InvocationError::FinalizationReceiptMismatch {
            invariant: "failed-publication",
        });
    }
    if report
        .failure
        .as_ref()
        .is_some_and(|failure| !failure_evidence_is_valid(failure))
    {
        return Err(InvocationError::FinalizationReceiptMismatch {
            invariant: "failure-evidence",
        });
    }
    Ok(())
}

fn finalized_child_receipt_root(receipt: &FinalizedChildReceipt) -> ContentHash {
    let mut hasher = DomainHasher::new(FINALIZED_CHILD_RECEIPT_DOMAIN);
    hash_field(
        &mut hasher,
        "invocation-root",
        receipt.invocation_root.as_bytes(),
    );
    hash_field(&mut hasher, "child-root", receipt.child.root().as_bytes());
    hash_field(
        &mut hasher,
        "finalization-root",
        receipt.finalization.root().as_bytes(),
    );
    hasher.finalize()
}

fn child_semantic_error(child: &ChildReceipt, invariant: &'static str) -> ReceiptSemanticError {
    ReceiptSemanticError::Child {
        ordinal: child.ordinal,
        invariant,
    }
}

fn invocation_semantic_error(invariant: &'static str) -> ReceiptSemanticError {
    ReceiptSemanticError::Invocation { invariant }
}

fn verify_receipt_child_count(children: usize) -> Result<u64, ReceiptSemanticError> {
    let children =
        u64::try_from(children).map_err(|_| ReceiptSemanticError::WorkLimitExceeded {
            children: u64::MAX,
            limit: INVOCATION_RECEIPT_MAX_CHILDREN,
        })?;
    if children > INVOCATION_RECEIPT_MAX_CHILDREN {
        return Err(ReceiptSemanticError::WorkLimitExceeded {
            children,
            limit: INVOCATION_RECEIPT_MAX_CHILDREN,
        });
    }
    Ok(children)
}

fn try_verifier_vec<T>(
    what: &'static str,
    requested_items: usize,
) -> Result<Vec<T>, ReceiptSemanticError> {
    let requested_items_u64 = u64::try_from(requested_items).unwrap_or(u64::MAX);
    let mut values = Vec::new();
    values.try_reserve_exact(requested_items).map_err(|_| {
        ReceiptSemanticError::AllocationRefused {
            what,
            requested_items: requested_items_u64,
        }
    })?;
    Ok(values)
}

struct ReceiptTopology {
    sorted_by_id: Vec<usize>,
    parents: Vec<Option<usize>>,
    nested_consumed: Vec<InvocationResources>,
    nested_granted: Vec<InvocationResources>,
    nested_returned: Vec<InvocationResources>,
    immediate_child_memory_peak: Vec<u64>,
    subtree_output_retained: Vec<u64>,
    descendant_has_output_grant: Vec<bool>,
}

impl ReceiptTopology {
    fn build(children: &[ChildReceipt]) -> Result<Self, ReceiptSemanticError> {
        verify_receipt_child_count(children.len())?;

        let mut sorted_by_id = try_verifier_vec("child-id-index", children.len())?;
        for index in 0..children.len() {
            sorted_by_id.push(index);
        }
        sorted_by_id.sort_unstable_by(|left, right| {
            children[*left]
                .id
                .as_bytes()
                .cmp(children[*right].id.as_bytes())
        });
        for pair in sorted_by_id.windows(2) {
            if children[pair[0]].id == children[pair[1]].id {
                let duplicate = pair[0].max(pair[1]);
                return Err(child_semantic_error(&children[duplicate], "unique-id"));
            }
        }

        let mut parents = try_verifier_vec("child-parent-index", children.len())?;
        for (index, child) in children.iter().enumerate() {
            let ordinal = u64::try_from(index)
                .map_err(|_| child_semantic_error(child, "ordinal-representable"))?;
            if child.ordinal != ordinal {
                return Err(child_semantic_error(child, "ordinal-order"));
            }
            if child.phase.is_empty() {
                return Err(child_semantic_error(child, "non-empty-phase"));
            }
            let parent = child
                .parent
                .map(|parent_id| {
                    Self::find_in_index(children, &sorted_by_id, parent_id)
                        .filter(|parent| *parent < index)
                        .ok_or_else(|| child_semantic_error(child, "parent-precedes-child"))
                })
                .transpose()?;
            parents.push(parent);
        }

        let mut nested_consumed = try_verifier_vec("nested-consumed-aggregates", children.len())?;
        let mut nested_granted = try_verifier_vec("nested-granted-aggregates", children.len())?;
        let mut nested_returned = try_verifier_vec("nested-returned-aggregates", children.len())?;
        let mut immediate_child_memory_peak =
            try_verifier_vec("child-memory-peak-aggregates", children.len())?;
        let mut subtree_output_retained =
            try_verifier_vec("subtree-output-aggregates", children.len())?;
        let mut descendant_has_output_grant =
            try_verifier_vec("descendant-output-grant-aggregates", children.len())?;
        for child in children {
            nested_consumed.push(InvocationResources::default());
            nested_granted.push(InvocationResources::default());
            nested_returned.push(InvocationResources::default());
            immediate_child_memory_peak.push(0);
            subtree_output_retained.push(child.output_retained);
            descendant_has_output_grant.push(false);
        }

        for index in (0..children.len()).rev() {
            let Some(parent) = parents[index] else {
                continue;
            };
            let parent_child = &children[parent];
            let child_subtree_output = subtree_output_retained[index];
            let child_subtree_has_output_grant =
                children[index].granted.output.get() != 0 || descendant_has_output_grant[index];
            nested_consumed[parent] = nested_consumed[parent]
                .checked_add(children[index].consumed)
                .map_err(|_| child_semantic_error(parent_child, "nested-consumption-sum"))?;
            nested_granted[parent] = nested_granted[parent]
                .checked_add(children[index].granted)
                .map_err(|_| child_semantic_error(parent_child, "nested-grant-sum"))?;
            nested_returned[parent] = nested_returned[parent]
                .checked_add(children[index].returned)
                .map_err(|_| child_semantic_error(parent_child, "nested-return-sum"))?;
            immediate_child_memory_peak[parent] =
                immediate_child_memory_peak[parent].max(children[index].memory_peak);
            subtree_output_retained[parent] = subtree_output_retained[parent]
                .checked_add(child_subtree_output)
                .ok_or_else(|| {
                    child_semantic_error(parent_child, "finalizer-subtree-output-sum")
                })?;
            descendant_has_output_grant[parent] |= child_subtree_has_output_grant;
        }

        Ok(Self {
            sorted_by_id,
            parents,
            nested_consumed,
            nested_granted,
            nested_returned,
            immediate_child_memory_peak,
            subtree_output_retained,
            descendant_has_output_grant,
        })
    }

    fn find_in_index(
        children: &[ChildReceipt],
        sorted_by_id: &[usize],
        id: ContentHash,
    ) -> Option<usize> {
        sorted_by_id
            .binary_search_by(|index| children[*index].id.as_bytes().cmp(id.as_bytes()))
            .ok()
            .map(|position| sorted_by_id[position])
    }

    fn index_of(&self, children: &[ChildReceipt], id: ContentHash) -> Option<usize> {
        Self::find_in_index(children, &self.sorted_by_id, id)
    }

    fn descends_from(&self, candidate: usize, ancestor: usize) -> bool {
        let mut parent = self.parents[candidate];
        for _ in 0..self.parents.len() {
            let Some(index) = parent else {
                return false;
            };
            if index == ancestor {
                return true;
            }
            parent = self.parents[index];
        }
        false
    }
}

fn failure_evidence_is_valid(error: &InvocationError) -> bool {
    match error {
        InvocationError::ResourceExceeded {
            resource,
            requested,
            available,
        } => {
            matches!(
                *resource,
                "work"
                    | "polls"
                    | "cost"
                    | "evaluations"
                    | "memory-bytes"
                    | "output-bytes"
                    | "child-count"
                    | "finalization-work"
                    | "finalization-polls"
            ) && requested > available
        }
        InvocationError::DeadlineExpired {
            deadline_ns,
            observed_ns,
            ..
        } => observed_ns >= deadline_ns,
        InvocationError::MemoryRefused {
            requested,
            used,
            limit,
            ..
        } => {
            *requested != 0
                && used <= limit
                && used
                    .checked_add(*requested)
                    .is_none_or(|total| total > *limit)
        }
        InvocationError::ArithmeticOverflow { resource } => matches!(
            *resource,
            "work"
                | "polls"
                | "cost"
                | "evaluations"
                | "memory-bytes"
                | "output-bytes"
                | "child-count"
                | "child-ordinal"
                | "live-children"
                | "memory-requested"
                | "memory-released"
                | "subtree-memory-bytes"
                | "output-retained"
                | "finalization-work"
                | "finalization-polls"
        ),
        InvocationError::TransactionalOutputScopeViolation {
            phase, requested, ..
        } => !phase.is_empty() && *requested != 0,
        InvocationError::EmptyPhase
        | InvocationError::InactiveChild
        | InvocationError::LiveNestedChildren { .. }
        | InvocationError::LiveMemoryReservations { .. }
        | InvocationError::UnfinishedChild { .. }
        | InvocationError::FinalizationRequired
        | InvocationError::PublicationNotPrepared
        | InvocationError::PublicationAlreadySealed
        | InvocationError::PublicationForbidden
        | InvocationError::FinalizationIncomplete { .. }
        | InvocationError::FinalizationReceiptMismatch { .. }
        | InvocationError::InvocationAlreadyFinalized
        | InvocationError::MemoryReleaseInvariant
        | InvocationError::EvidenceAllocationRefused { .. } => false,
        InvocationError::Cancelled { .. } | InvocationError::ExplicitRefusal { .. } => true,
    }
}

fn memory_refusal_matches_failure(
    refusal: &InvocationMemoryRefusal,
    failure: &InvocationError,
) -> bool {
    matches!(
        failure,
        InvocationError::MemoryRefused {
            what,
            requested,
            used,
            limit,
        } if *what == refusal.what
            && *requested == refusal.requested
            && *used == refusal.used
            && *limit == refusal.limit
    )
}

fn failure_requires_child_origin(failure: &InvocationError) -> bool {
    matches!(
        failure,
        InvocationError::MemoryRefused { .. }
            | InvocationError::ExplicitRefusal { .. }
            | InvocationError::TransactionalOutputScopeViolation { .. }
    ) || matches!(
        failure,
        InvocationError::ResourceExceeded { resource, .. }
            | InvocationError::ArithmeticOverflow { resource }
            if matches!(*resource, "finalization-work" | "finalization-polls")
    )
}

fn failure_requires_finalizable_origin(failure: &InvocationError) -> bool {
    matches!(
        failure,
        InvocationError::ResourceExceeded { resource, .. }
            | InvocationError::ArithmeticOverflow { resource }
            if matches!(*resource, "finalization-work" | "finalization-polls")
    )
}

fn resource_capacity_for_failure(
    receipt: &InvocationReceipt,
    origin: Option<&ChildReceipt>,
    resource: &str,
) -> Option<u128> {
    let scientific = origin.map_or(receipt.required, |child| {
        child
            .finalization
            .as_ref()
            .map_or(child.granted, |finalization| {
                finalization.scientific_granted
            })
    });
    match resource {
        "work" => Some(scientific.work.get()),
        "polls" => Some(u128::from(scientific.polls.get())),
        "cost" => Some(u128::from(scientific.cost.get())),
        "evaluations" => Some(u128::from(scientific.evaluations.get())),
        "memory-bytes" => Some(u128::from(scientific.memory.get())),
        "output-bytes" => Some(u128::from(scientific.output.get())),
        "child-count" => Some(u128::from(INVOCATION_RECEIPT_MAX_CHILDREN)),
        "finalization-work" => origin
            .and_then(|child| child.finalization.as_ref())
            .map(|finalization| finalization.granted.work.get()),
        "finalization-polls" => origin
            .and_then(|child| child.finalization.as_ref())
            .map(|finalization| u128::from(finalization.granted.polls.get())),
        _ => None,
    }
}

fn verify_deadline_semantics(receipt: &InvocationReceipt) -> Result<(), ReceiptSemanticError> {
    match (receipt.limits.deadline, receipt.last_deadline_observation) {
        (None, None) => {
            if matches!(
                &receipt.failure,
                Some(InvocationError::DeadlineExpired { .. })
            ) {
                return Err(invocation_semantic_error("deadline-without-limit"));
            }
        }
        (Some(deadline), Some(observed)) => {
            if let Some(InvocationError::DeadlineExpired {
                deadline_ns,
                observed_ns,
                ..
            }) = &receipt.failure
            {
                if *deadline_ns != deadline.as_nanos() || *observed_ns != observed.as_nanos() {
                    return Err(invocation_semantic_error("deadline-failure-observation"));
                }
            } else if observed >= deadline {
                return Err(invocation_semantic_error(
                    "nondeadline-observation-before-limit",
                ));
            }
        }
        _ => return Err(invocation_semantic_error("deadline-observation-presence")),
    }
    Ok(())
}

fn verify_failure_propagation(
    receipt: &InvocationReceipt,
    topology: &ReceiptTopology,
) -> Result<(), ReceiptSemanticError> {
    if let Some(failure) = &receipt.failure {
        let origin_index = receipt
            .failure_origin
            .and_then(|origin| topology.index_of(&receipt.children, origin));
        let origin = origin_index.map(|index| &receipt.children[index]);
        if failure_requires_child_origin(failure) && receipt.failure_origin.is_none() {
            return Err(invocation_semantic_error("failure-origin-required"));
        }
        if failure_requires_finalizable_origin(failure) {
            if origin.is_some_and(|child| child.finalization.is_none()) {
                return Err(invocation_semantic_error(
                    "finalization-failure-origin-kind",
                ));
            }
        }
        if let InvocationError::TransactionalOutputScopeViolation { ancestor, .. } = failure {
            let ancestor_index = topology
                .index_of(&receipt.children, *ancestor)
                .ok_or_else(|| invocation_semantic_error("transactional-output-ancestor-exists"))?;
            let ancestor_child = &receipt.children[ancestor_index];
            if ancestor_child.finalization.is_none() {
                return Err(child_semantic_error(
                    ancestor_child,
                    "transactional-output-ancestor-finalizable",
                ));
            }
            if origin_index.is_some_and(|origin| {
                origin != ancestor_index && !topology.descends_from(origin, ancestor_index)
            }) {
                return Err(invocation_semantic_error(
                    "transactional-output-origin-descends",
                ));
            }
        }
        if let InvocationError::ResourceExceeded {
            resource,
            available,
            ..
        } = failure
        {
            if resource_capacity_for_failure(receipt, origin, resource)
                .is_some_and(|capacity| *available > capacity)
            {
                return Err(origin.map_or_else(
                    || invocation_semantic_error("failure-available-within-origin-grant"),
                    |child| child_semantic_error(child, "failure-available-within-origin-grant"),
                ));
            }
        }
    }
    match (&receipt.failure, receipt.failure_origin) {
        (None | Some(_), None) => {}
        (None, Some(_)) => {
            return Err(invocation_semantic_error("failure-origin-without-failure"));
        }
        (Some(failure), Some(origin)) => {
            let origin = topology
                .index_of(&receipt.children, origin)
                .map(|index| &receipt.children[index])
                .ok_or_else(|| invocation_semantic_error("failure-origin-exists"))?;
            if origin.failure.as_ref() != Some(failure) {
                return Err(child_semantic_error(origin, "failure-origin-matches"));
            }
        }
    }
    for (index, child) in receipt.children.iter().enumerate() {
        let Some(failure) = &child.failure else {
            if child.failure_inherited {
                return Err(child_semantic_error(
                    child,
                    "inherited-flag-without-failure",
                ));
            }
            continue;
        };
        let is_origin = receipt.failure_origin == Some(child.id);
        if is_origin == child.failure_inherited {
            return Err(child_semantic_error(child, "failure-origin-marker"));
        }
        if receipt.failure.as_ref() != Some(failure) {
            return Err(child_semantic_error(child, "failure-propagates-to-root"));
        }
        if let Some(parent) = topology.parents[index] {
            let parent = &receipt.children[parent];
            if parent.failure.as_ref() != Some(failure) {
                return Err(child_semantic_error(child, "failure-propagates-to-parent"));
            }
        }
    }
    Ok(())
}

#[allow(clippy::too_many_lines)]
fn verify_receipt_semantics(receipt: &InvocationReceipt) -> Result<(), ReceiptSemanticError> {
    if receipt.version != INVOCATION_RECEIPT_VERSION {
        return Err(ReceiptSemanticError::UnsupportedVersion {
            found: receipt.version,
        });
    }
    verify_receipt_child_count(receipt.children.len())?;
    if invocation_receipt_root(receipt) != receipt.root {
        return Err(ReceiptSemanticError::RootMismatch);
    }
    receipt
        .limits
        .resources
        .checked_sub(receipt.required)
        .map_err(|_| invocation_semantic_error("required-within-limits"))?;
    verify_deadline_semantics(receipt)?;
    let topology = ReceiptTopology::build(&receipt.children)?;

    for (index, child) in receipt.children.iter().enumerate() {
        let (scientific_grant, finalization_grant) = child.finalization.as_ref().map_or(
            (child.granted, FinalizationResources::default()),
            |finalization| (finalization.scientific_granted, finalization.granted),
        );
        if child_id(
            receipt.invocation_id,
            receipt.plan_binding,
            child.parent,
            child.ordinal,
            child.phase,
            scientific_grant,
            finalization_grant,
        ) != child.id
        {
            return Err(child_semantic_error(child, "derived-id"));
        }
        if child_receipt_root(child) != child.root {
            return Err(child_semantic_error(child, "receipt-root"));
        }
        let consumed = child
            .granted
            .checked_sub(child.returned)
            .map_err(|_| child_semantic_error(child, "granted-returned"))?;
        if consumed != child.consumed {
            return Err(child_semantic_error(child, "consumed-definition"));
        }
        let expected_consumed = child
            .direct_consumed
            .checked_add(topology.nested_consumed[index])
            .map_err(|_| child_semantic_error(child, "direct-plus-nested-consumption"))?;
        if expected_consumed != child.consumed {
            return Err(child_semantic_error(child, "subtree-conservation"));
        }
        let mut replay = child
            .granted
            .checked_sub(topology.nested_granted[index])
            .and_then(|available| available.checked_add(topology.nested_returned[index]))
            .map_err(|_| child_semantic_error(child, "nested-affine-transfer"))?;
        replay = replay
            .checked_sub(child.direct_consumed)
            .map_err(|_| child_semantic_error(child, "direct-affine-spend"))?;
        if replay != child.returned {
            return Err(child_semantic_error(child, "returned-conservation"));
        }
        if child.direct_consumed.memory != MemoryBytes::new(0)
            || child.consumed.memory != MemoryBytes::new(0)
            || child.returned.memory != child.granted.memory
        {
            return Err(child_semantic_error(child, "memory-is-reusable-capacity"));
        }
        let descendant_memory_peak = topology.immediate_child_memory_peak[index];
        let minimum_memory_peak = child.direct_memory_peak.max(descendant_memory_peak);
        let maximum_memory_peak = child
            .direct_memory_peak
            .checked_add(descendant_memory_peak)
            .ok_or_else(|| child_semantic_error(child, "memory-peak-bound"))?;
        if child.memory_requested != child.memory_released
            || child.memory_peak < minimum_memory_peak
            || child.memory_peak > maximum_memory_peak
            || child.memory_peak > child.granted.memory.0
            || u128::from(child.direct_memory_peak) > child.memory_requested
            || (child.memory_requested == 0) != (child.direct_memory_peak == 0)
        {
            return Err(child_semantic_error(child, "memory-receipt"));
        }
        if child.direct_consumed.output.0 != child.output_retained {
            return Err(child_semantic_error(child, "direct-output-retention"));
        }
        if let Some(finalization) = &child.finalization {
            if topology.subtree_output_retained[index] != child.output_retained {
                return Err(child_semantic_error(child, "finalizer-descendant-output"));
            }
            if topology.descendant_has_output_grant[index] {
                return Err(child_semantic_error(
                    child,
                    "finalizer-descendant-output-grant",
                ));
            }
            let finalizer_consumed = finalization
                .granted
                .checked_sub(finalization.returned)
                .map_err(|_| child_semantic_error(child, "finalizer-granted-returned"))?;
            if finalizer_consumed != finalization.consumed {
                return Err(child_semantic_error(child, "finalizer-consumed-definition"));
            }
            let replay_granted = finalization
                .scientific_granted
                .checked_add(finalization.granted.as_invocation_resources())
                .map_err(|_| child_semantic_error(child, "finalizer-grant-partition"))?;
            let replay_direct = finalization
                .scientific_direct_consumed
                .checked_add(finalization.consumed.as_invocation_resources())
                .map_err(|_| child_semantic_error(child, "finalizer-direct-partition"))?;
            let replay_returned = finalization
                .scientific_returned
                .checked_add(finalization.returned.as_invocation_resources())
                .map_err(|_| child_semantic_error(child, "finalizer-return-partition"))?;
            if replay_granted != child.granted {
                return Err(child_semantic_error(child, "finalizer-grant-partition"));
            }
            if replay_direct != child.direct_consumed {
                return Err(child_semantic_error(child, "finalizer-direct-partition"));
            }
            if replay_returned != child.returned {
                return Err(child_semantic_error(child, "finalizer-return-partition"));
            }
            let bound_report = FinalizationReport {
                version: FINALIZATION_REPORT_VERSION,
                plan_binding: receipt.plan_binding,
                child: child.id,
                granted: finalization.granted,
                consumed: finalization.consumed,
                returned: finalization.returned,
                publication_scope: finalization.publication_scope,
                publication: finalization.publication,
                failure: child.failure.clone(),
                disposition: child.disposition,
                root: finalization.report_root,
            };
            if finalization_report_semantics(&bound_report).is_err() {
                return Err(child_semantic_error(child, "finalizer-report-root"));
            }
            match finalization.publication {
                FinalizationPublication::Pending | FinalizationPublication::Prepared => {
                    return Err(child_semantic_error(
                        child,
                        "finalizer-publication-terminal",
                    ));
                }
                FinalizationPublication::Aborted if child.output_retained != 0 => {
                    return Err(child_semantic_error(child, "finalizer-aborted-output"));
                }
                FinalizationPublication::Committed { bytes }
                    if child.output_retained != bytes.get() =>
                {
                    return Err(child_semantic_error(child, "finalizer-committed-output"));
                }
                FinalizationPublication::Committed { .. } if child.failure.is_some() => {
                    return Err(child_semantic_error(child, "finalizer-failed-publication"));
                }
                FinalizationPublication::Aborted | FinalizationPublication::Committed { .. } => {}
            }
        }
        let expected_disposition = child
            .failure
            .as_ref()
            .map_or(InvocationDisposition::Completed, error_disposition);
        if child.disposition != expected_disposition {
            return Err(child_semantic_error(child, "derived-disposition"));
        }
        if child
            .failure
            .as_ref()
            .is_some_and(|failure| !failure_evidence_is_valid(failure))
        {
            return Err(child_semantic_error(child, "failure-evidence"));
        }
    }

    let mut replay = receipt.required;
    for child in receipt
        .children
        .iter()
        .filter(|candidate| candidate.parent.is_none())
    {
        replay = replay
            .checked_sub(child.granted)
            .and_then(|available| available.checked_add(child.returned))
            .map_err(|_| invocation_semantic_error("root-affine-transfer"))?;
    }
    if replay != receipt.remaining {
        return Err(invocation_semantic_error("root-conservation"));
    }

    let (memory_requested, memory_released, output_retained) = receipt.children.iter().try_fold(
        (0_u128, 0_u128, 0_u64),
        |(requested, released, output), child| {
            Ok::<_, ReceiptSemanticError>((
                requested
                    .checked_add(child.memory_requested)
                    .ok_or_else(|| invocation_semantic_error("memory-requested-sum"))?,
                released
                    .checked_add(child.memory_released)
                    .ok_or_else(|| invocation_semantic_error("memory-released-sum"))?,
                output
                    .checked_add(child.output_retained)
                    .ok_or_else(|| invocation_semantic_error("output-retained-sum"))?,
            ))
        },
    )?;
    let memory_peak = receipt
        .children
        .iter()
        .filter(|child| child.parent.is_none())
        .map(|child| child.memory_peak)
        .max()
        .unwrap_or(0);
    if memory_requested != receipt.memory_requested
        || memory_released != receipt.memory_released
        || receipt.memory_requested != receipt.memory_released
        || memory_peak != receipt.memory_peak
        || receipt.memory_peak > receipt.required.memory.0
        || receipt.remaining.memory != receipt.required.memory
    {
        return Err(invocation_semantic_error("root-memory-receipt"));
    }
    if output_retained != receipt.output_retained
        || receipt
            .required
            .output
            .0
            .checked_sub(receipt.remaining.output.0)
            != Some(receipt.output_retained)
    {
        return Err(invocation_semantic_error("root-output-receipt"));
    }
    if receipt.memory_refusals > 1
        || (receipt.memory_refusals == 0) != receipt.memory_first_refusal.is_none()
    {
        return Err(invocation_semantic_error("memory-refusal-evidence"));
    }
    match (&receipt.memory_first_refusal, &receipt.failure) {
        (Some(refusal), Some(failure))
            if memory_refusal_matches_failure(refusal, failure)
                && refusal.limit == receipt.required.memory.0
                && refusal.used <= receipt.memory_peak
                && receipt
                    .children
                    .iter()
                    .any(|child| child.failure.as_ref() == Some(failure)) => {}
        (None, Some(InvocationError::MemoryRefused { .. })) | (Some(_), _) => {
            return Err(invocation_semantic_error("memory-refusal-first-fault"));
        }
        _ => {}
    }
    let expected_disposition = receipt
        .failure
        .as_ref()
        .map_or(InvocationDisposition::Completed, error_disposition);
    if receipt.disposition != expected_disposition {
        return Err(invocation_semantic_error("root-derived-disposition"));
    }
    if receipt
        .children
        .iter()
        .any(|child| child.disposition != InvocationDisposition::Completed)
        && receipt.failure.is_none()
    {
        return Err(invocation_semantic_error("child-failure-propagates"));
    }
    if receipt
        .failure
        .as_ref()
        .is_some_and(|failure| !failure_evidence_is_valid(failure))
    {
        return Err(invocation_semantic_error("root-failure-evidence"));
    }
    verify_failure_propagation(receipt, &topology)?;
    Ok(())
}

fn invocation_receipt_root(receipt: &InvocationReceipt) -> ContentHash {
    let mut hasher = DomainHasher::new(INVOCATION_RECEIPT_DOMAIN);
    hash_field(&mut hasher, "version", &receipt.version.to_le_bytes());
    hash_field(
        &mut hasher,
        "invocation-id",
        receipt.invocation_id.as_bytes(),
    );
    hash_plan_binding(&mut hasher, receipt.plan_binding);
    hash_resources(
        &mut hasher,
        [
            "limits.work",
            "limits.polls",
            "limits.cost",
            "limits.evaluations",
            "limits.memory",
            "limits.output",
        ],
        receipt.limits.resources,
    );
    hash_field(
        &mut hasher,
        "deadline-present",
        &[u8::from(receipt.limits.deadline.is_some())],
    );
    if let Some(deadline) = receipt.limits.deadline {
        hash_field(
            &mut hasher,
            "deadline-nanos",
            &deadline.as_nanos().to_le_bytes(),
        );
    }
    hash_field(
        &mut hasher,
        "accuracy-obligation",
        receipt.limits.accuracy_obligation.as_bytes(),
    );
    hash_field(
        &mut hasher,
        "capability-scope",
        receipt.limits.capability_scope.as_bytes(),
    );
    hash_resources(
        &mut hasher,
        [
            "required.work",
            "required.polls",
            "required.cost",
            "required.evaluations",
            "required.memory",
            "required.output",
        ],
        receipt.required,
    );
    hash_resources(
        &mut hasher,
        [
            "remaining.work",
            "remaining.polls",
            "remaining.cost",
            "remaining.evaluations",
            "remaining.memory",
            "remaining.output",
        ],
        receipt.remaining,
    );
    hash_field(
        &mut hasher,
        "last-deadline-observation-present",
        &[u8::from(receipt.last_deadline_observation.is_some())],
    );
    if let Some(observed) = receipt.last_deadline_observation {
        hash_field(
            &mut hasher,
            "last-deadline-observation-nanos",
            &observed.as_nanos().to_le_bytes(),
        );
    }
    hash_field(
        &mut hasher,
        "memory-peak",
        &receipt.memory_peak.to_le_bytes(),
    );
    hash_field(
        &mut hasher,
        "memory-requested",
        &receipt.memory_requested.to_le_bytes(),
    );
    hash_field(
        &mut hasher,
        "memory-released",
        &receipt.memory_released.to_le_bytes(),
    );
    hash_field(
        &mut hasher,
        "memory-refusals",
        &receipt.memory_refusals.to_le_bytes(),
    );
    hash_field(
        &mut hasher,
        "memory-first-refusal-present",
        &[u8::from(receipt.memory_first_refusal.is_some())],
    );
    if let Some(refusal) = &receipt.memory_first_refusal {
        hash_field(
            &mut hasher,
            "memory-first-refusal.what",
            refusal.what.as_bytes(),
        );
        hash_field(
            &mut hasher,
            "memory-first-refusal.requested",
            &refusal.requested.to_le_bytes(),
        );
        hash_field(
            &mut hasher,
            "memory-first-refusal.used",
            &refusal.used.to_le_bytes(),
        );
        hash_field(
            &mut hasher,
            "memory-first-refusal.limit",
            &refusal.limit.to_le_bytes(),
        );
    }
    hash_field(
        &mut hasher,
        "output-retained",
        &receipt.output_retained.to_le_bytes(),
    );
    hash_field(
        &mut hasher,
        "failure-origin-present",
        &[u8::from(receipt.failure_origin.is_some())],
    );
    if let Some(origin) = receipt.failure_origin {
        hash_field(&mut hasher, "failure-origin", origin.as_bytes());
    }
    hash_field(
        &mut hasher,
        "failure-present",
        &[u8::from(receipt.failure.is_some())],
    );
    if let Some(failure) = &receipt.failure {
        hash_field(
            &mut hasher,
            "failure-root",
            invocation_error_root(failure).as_bytes(),
        );
    }
    hash_field(
        &mut hasher,
        "disposition",
        &[encode_disposition(receipt.disposition)],
    );
    hash_field(
        &mut hasher,
        "child-count",
        &(receipt.children.len() as u64).to_le_bytes(),
    );
    for child in &receipt.children {
        hash_field(&mut hasher, "child.id", child.id.as_bytes());
        hash_field(&mut hasher, "child.root", child.root.as_bytes());
    }
    hasher.finalize()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Budget, CancelGate, ExecMode, StreamKey};
    use fs_alloc::{ArenaConfig, ArenaPool};

    fn resources(value: u64) -> InvocationResources {
        InvocationResources::new(
            WorkUnits::new(u128::from(value)),
            PollUnits::new(value as u32),
            CostUnits::new(value),
            EvaluationUnits::new(value),
            MemoryBytes::new(value),
            OutputBytes::new(value),
        )
    }

    fn resource_vector(values: [u64; 6]) -> InvocationResources {
        let [work, polls, cost, evaluations, memory, output] = values;
        InvocationResources::new(
            WorkUnits::new(u128::from(work)),
            PollUnits::new(u32::try_from(polls).expect("test poll value fits u32")),
            CostUnits::new(cost),
            EvaluationUnits::new(evaluations),
            MemoryBytes::new(memory),
            OutputBytes::new(output),
        )
    }

    fn identities() -> (ContentHash, ContentHash, ContentHash) {
        (
            hash_domain("test.invocation", b"id"),
            hash_domain("test.accuracy", b"obligation"),
            hash_domain("test.capability", b"scope"),
        )
    }

    fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
        let gate = CancelGate::new_clock_free();
        let pool = ArenaPool::new(ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 1,
                    kernel_id: 2,
                    tile: 3,
                    iteration: 4,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            f(&cx)
        })
    }

    fn with_gate_cx<R>(f: impl FnOnce(&CancelGate, &Cx<'_>) -> R) -> R {
        let gate = CancelGate::new_clock_free();
        let pool = ArenaPool::new(ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 1,
                    kernel_id: 2,
                    tile: 3,
                    iteration: 4,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            f(&gate, &cx)
        })
    }

    fn with_leased_cx<R>(lease: &OperationMemoryLease, f: impl FnOnce(&Cx<'_>) -> R) -> R {
        let gate = CancelGate::new_clock_free();
        let refusals = crate::cx::RefusalSink::default();
        let pool = ArenaPool::new(ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new_with_refusal_sink(
                &gate,
                arena,
                StreamKey {
                    seed: 1,
                    kernel_id: 2,
                    tile: 3,
                    iteration: 4,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
                &refusals,
                lease,
            );
            f(&cx)
        })
    }

    #[test]
    fn affine_children_conserve_each_dimension_and_memory_releases() {
        let clock = VirtualClock::new();
        let (id, accuracy, capability) = identities();
        let limits = InvocationLimits::new(resources(20), None, accuracy, capability);
        let receipt = with_cx(|cx| {
            let admission = InvocationAdmitter::new()
                .admit(id, limits, resources(10))
                .unwrap();
            let mut root = admission.begin(cx, &clock).unwrap();
            {
                let mut child = root.split_child("phase-a", resources(10)).unwrap();
                child.charge_work(WorkUnits::new(7)).unwrap();
                child.charge_cost(CostUnits::new(5)).unwrap();
                child.charge_evaluations(EvaluationUnits::new(2)).unwrap();
                child.poll("phase-a.poll").unwrap();
                {
                    let mut memory = child
                        .reserve_memory("invocation-test", MemoryBytes::new(8))
                        .unwrap();
                    memory.budget().publish_output(OutputBytes::new(3)).unwrap();
                }
                assert_eq!(child.finish().unwrap(), InvocationDisposition::Completed);
            }
            root.finish().unwrap()
        });
        assert!(receipt.verifies_integrity());
        assert_eq!(receipt.children().len(), 1);
        assert_eq!(receipt.children()[0].consumed().work(), WorkUnits::new(7));
        assert_eq!(receipt.children()[0].memory_peak_bytes(), 8);
        assert_eq!(receipt.output_retained_bytes(), 3);
    }

    #[test]
    fn admission_and_deadline_refusals_are_exact_and_ordered() {
        let clock = VirtualClock::starting_at(Time::from_nanos(5));
        let (id, accuracy, capability) = identities();
        let limits = InvocationLimits::new(
            resources(9),
            Some(Time::from_nanos(10)),
            accuracy,
            capability,
        );
        assert!(matches!(
            InvocationAdmitter::new().admit(id, limits, resources(10)),
            Err(InvocationError::ResourceExceeded {
                resource: "work",
                requested: 10,
                available: 9
            })
        ));
        let expired = InvocationLimits::new(
            resources(10),
            Some(Time::from_nanos(5)),
            accuracy,
            capability,
        );
        with_cx(|cx| {
            let admission = InvocationAdmitter::new()
                .admit(id, expired, resources(10))
                .unwrap();
            assert!(matches!(
                admission.begin(cx, &clock),
                Err(InvocationError::DeadlineExpired {
                    phase: "invocation-admission",
                    deadline_ns: 5,
                    observed_ns: 5
                })
            ));
        });
    }

    #[test]
    fn admission_refuses_one_below_in_every_resource_dimension() {
        let (id, accuracy, capability) = identities();
        let required = resources(10);
        let one_below = [
            (resource_vector([9, 10, 10, 10, 10, 10]), "work"),
            (resource_vector([10, 9, 10, 10, 10, 10]), "polls"),
            (resource_vector([10, 10, 9, 10, 10, 10]), "cost"),
            (resource_vector([10, 10, 10, 9, 10, 10]), "evaluations"),
            (resource_vector([10, 10, 10, 10, 9, 10]), "memory-bytes"),
            (resource_vector([10, 10, 10, 10, 10, 9]), "output-bytes"),
        ];
        for (available, resource) in one_below {
            let limits = InvocationLimits::new(available, None, accuracy, capability);
            assert!(matches!(
                InvocationAdmitter::new().admit(id, limits, required),
                Err(InvocationError::ResourceExceeded {
                    resource: observed,
                    requested: 10,
                    available: 9,
                }) if observed == resource
            ));
        }
    }

    #[test]
    fn child_runtime_refuses_overrun_in_every_resource_dimension() {
        for resource in [
            "work",
            "polls",
            "cost",
            "evaluations",
            "memory-bytes",
            "output-bytes",
        ] {
            let clock = VirtualClock::new();
            let (id, accuracy, capability) = identities();
            let receipt = with_cx(|cx| {
                let admission = InvocationAdmitter::new()
                    .admit(
                        id,
                        InvocationLimits::new(resources(4), None, accuracy, capability),
                        resources(4),
                    )
                    .unwrap();
                let mut root = admission.begin(cx, &clock).unwrap();
                let mut child = root.split_child("overrun", resources(4)).unwrap();
                let failure = match resource {
                    "work" => child.charge_work(WorkUnits::new(5)),
                    "polls" => {
                        let mut result = Ok(());
                        for _ in 0..5 {
                            result = child.poll("overrun.poll");
                            if result.is_err() {
                                break;
                            }
                        }
                        result
                    }
                    "cost" => child.charge_cost(CostUnits::new(5)),
                    "evaluations" => child.charge_evaluations(EvaluationUnits::new(5)),
                    "memory-bytes" => child
                        .reserve_memory("overrun-memory", MemoryBytes::new(5))
                        .map(drop),
                    "output-bytes" => child.publish_output(OutputBytes::new(5)),
                    _ => unreachable!(),
                };
                assert!(matches!(
                    failure,
                    Err(InvocationError::ResourceExceeded {
                        resource: observed,
                        requested,
                        available,
                    }) if observed == resource && requested > available
                ));
                assert_eq!(child.finish().unwrap(), InvocationDisposition::Refused);
                root.finish().unwrap()
            });
            assert_eq!(receipt.disposition(), InvocationDisposition::Refused);
            assert!(receipt.verifies_integrity());
        }
    }

    #[test]
    fn empty_child_phase_is_rejected_before_identity_or_ordinal_mutation() {
        let clock = VirtualClock::new();
        let (id, accuracy, capability) = identities();
        let receipt = with_cx(|cx| {
            let admission = InvocationAdmitter::new()
                .admit(
                    id,
                    InvocationLimits::new(resources(4), None, accuracy, capability),
                    resources(4),
                )
                .unwrap();
            let mut root = admission.begin(cx, &clock).unwrap();
            assert!(matches!(
                root.split_child("", resources(4)),
                Err(InvocationError::EmptyPhase)
            ));
            let child = root.split_child("valid", resources(4)).unwrap();
            assert_eq!(child.finish().unwrap(), InvocationDisposition::Completed);
            root.finish().unwrap()
        });
        assert_eq!(receipt.children().len(), 1);
        assert_eq!(receipt.children()[0].ordinal(), 0);
        assert_eq!(receipt.children()[0].phase(), "valid");
        assert!(receipt.verifies_integrity());
    }

    #[test]
    fn root_memory_is_reserved_once_against_the_ambient_operation_lease() {
        let clock = VirtualClock::new();
        let (id, accuracy, capability) = identities();
        let occupied_lease = OperationMemoryLease::bounded(10);
        let occupied = occupied_lease.reserve("existing-operation", 4).unwrap();
        with_leased_cx(&occupied_lease, |cx| {
            let admission = InvocationAdmitter::new()
                .admit(
                    id,
                    InvocationLimits::new(resources(10), None, accuracy, capability),
                    resources(8),
                )
                .unwrap();
            assert!(matches!(
                admission.begin(cx, &clock),
                Err(InvocationError::MemoryRefused {
                    what: "invocation-root-memory",
                    requested: 8,
                    used: 4,
                    limit: 10,
                })
            ));
        });
        drop(occupied);

        let admitted_lease = OperationMemoryLease::bounded(10);
        with_leased_cx(&admitted_lease, |cx| {
            let admission = InvocationAdmitter::new()
                .admit(
                    id,
                    InvocationLimits::new(resources(10), None, accuracy, capability),
                    resources(8),
                )
                .unwrap();
            let mut root = admission.begin(cx, &clock).unwrap();
            assert_eq!(admitted_lease.receipt().used_bytes, 8);
            let receipt = root.finish().unwrap();
            assert!(receipt.verifies_integrity());
            assert_eq!(admitted_lease.receipt().used_bytes, 0);
        });
    }

    #[test]
    fn nested_child_ids_and_receipts_replay_deterministically() {
        fn run() -> InvocationReceipt {
            let clock = VirtualClock::new();
            let (id, accuracy, capability) = identities();
            let limits = InvocationLimits::new(resources(12), None, accuracy, capability);
            with_cx(|cx| {
                let admission = InvocationAdmitter::new()
                    .admit(id, limits, resources(12))
                    .unwrap();
                let mut root = admission.begin(cx, &clock).unwrap();
                {
                    let mut parent = root.split_child("parent", resources(12)).unwrap();
                    {
                        let nested = parent.split_child("nested", resources(4)).unwrap();
                        assert_eq!(nested.finish().unwrap(), InvocationDisposition::Completed);
                    }
                    assert_eq!(parent.finish().unwrap(), InvocationDisposition::Completed);
                }
                root.finish().unwrap()
            })
        }
        let first = run();
        let second = run();
        assert_eq!(first, second);
        assert!(first.verifies_integrity());
    }

    #[test]
    fn receipt_work_and_allocation_limits_refuse_structurally() {
        assert!(matches!(
            verify_receipt_child_count(
                usize::try_from(INVOCATION_RECEIPT_MAX_CHILDREN + 1)
                    .expect("configured child limit fits this target"),
            ),
            Err(ReceiptSemanticError::WorkLimitExceeded {
                children,
                limit: INVOCATION_RECEIPT_MAX_CHILDREN,
            }) if children == INVOCATION_RECEIPT_MAX_CHILDREN + 1
        ));
        assert!(matches!(
            try_verifier_vec::<u8>("verifier-test", usize::MAX),
            Err(ReceiptSemanticError::AllocationRefused {
                what: "verifier-test",
                requested_items,
            }) if requested_items == u64::try_from(usize::MAX).unwrap_or(u64::MAX)
        ));
        assert!(matches!(
            try_evidence_vec::<u8>("producer-test", usize::MAX),
            Err(InvocationError::EvidenceAllocationRefused {
                what: "producer-test",
                requested_items,
            }) if requested_items == u64::try_from(usize::MAX).unwrap_or(u64::MAX)
        ));
    }

    #[test]
    fn semantic_verifier_handles_deep_and_wide_receipts_with_one_index_pass() {
        fn synthetic_receipt(children_count: usize, deep: bool) -> InvocationReceipt {
            let invocation_id = hash_domain("test.verifier-shape.invocation", b"bounded");
            let accuracy = hash_domain("test.verifier-shape.accuracy", b"exact");
            let capability = hash_domain("test.verifier-shape.capability", b"unit");
            let resources = InvocationResources::default();
            let mut children: Vec<ChildReceipt> = Vec::with_capacity(children_count);
            for index in 0..children_count {
                let parent = deep
                    .then(|| index.checked_sub(1).map(|parent| children[parent].id))
                    .flatten();
                let ordinal = u64::try_from(index).expect("test child ordinal is representable");
                let id = child_id(
                    invocation_id,
                    None,
                    parent,
                    ordinal,
                    "bounded-verifier-shape",
                    resources,
                    FinalizationResources::default(),
                );
                let mut child = ChildReceipt {
                    id,
                    parent,
                    ordinal,
                    phase: "bounded-verifier-shape",
                    granted: resources,
                    consumed: resources,
                    direct_consumed: resources,
                    returned: resources,
                    direct_memory_peak: 0,
                    memory_peak: 0,
                    memory_requested: 0,
                    memory_released: 0,
                    output_retained: 0,
                    finalization: None,
                    failure: None,
                    failure_inherited: false,
                    disposition: InvocationDisposition::Completed,
                    root: ContentHash([0; 32]),
                };
                child.root = child_receipt_root(&child);
                children.push(child);
            }
            let limits = InvocationLimits::new(resources, None, accuracy, capability);
            let mut receipt = InvocationReceipt {
                version: INVOCATION_RECEIPT_VERSION,
                invocation_id,
                plan_binding: None,
                limits,
                required: resources,
                remaining: resources,
                children,
                last_deadline_observation: None,
                memory_peak: 0,
                memory_requested: 0,
                memory_released: 0,
                memory_refusals: 0,
                memory_first_refusal: None,
                output_retained: 0,
                failure: None,
                failure_origin: None,
                disposition: InvocationDisposition::Completed,
                root: ContentHash([0; 32]),
            };
            receipt.root = invocation_receipt_root(&receipt);
            receipt
        }

        let deep = synthetic_receipt(4_096, true);
        let wide = synthetic_receipt(4_096, false);
        assert!(deep.verify_semantics().is_ok());
        assert!(wide.verify_semantics().is_ok());
    }

    #[test]
    fn first_fault_latches_and_derives_refused_receipts() {
        let clock = VirtualClock::new();
        let (id, accuracy, capability) = identities();
        let receipt = with_cx(|cx| {
            let admission = InvocationAdmitter::new()
                .admit(
                    id,
                    InvocationLimits::new(resources(10), None, accuracy, capability),
                    resources(10),
                )
                .unwrap();
            let mut root = admission.begin(cx, &clock).unwrap();
            let mut child = root.split_child("overrun", resources(10)).unwrap();
            assert!(matches!(
                child.charge_work(WorkUnits::new(11)),
                Err(InvocationError::ResourceExceeded {
                    resource: "work",
                    requested: 11,
                    available: 10,
                })
            ));
            assert_eq!(child.finish().unwrap(), InvocationDisposition::Refused);
            root.finish().unwrap()
        });
        assert_eq!(receipt.disposition(), InvocationDisposition::Refused);
        assert!(receipt.failure().is_some());
        assert!(receipt.verifies_integrity());
    }

    #[test]
    fn cancellation_after_one_poll_drains_and_cannot_complete() {
        let clock = VirtualClock::new();
        let (id, accuracy, capability) = identities();
        let receipt = with_gate_cx(|gate, cx| {
            let admission = InvocationAdmitter::new()
                .admit(
                    id,
                    InvocationLimits::new(resources(4), None, accuracy, capability),
                    resources(4),
                )
                .unwrap();
            let mut root = admission.begin(cx, &clock).unwrap();
            let mut child = root.split_child("cancelled", resources(4)).unwrap();
            gate.request();
            assert!(matches!(
                child.poll("cancelled.poll"),
                Err(InvocationError::Cancelled {
                    phase: "cancelled.poll"
                })
            ));
            assert_eq!(child.finish().unwrap(), InvocationDisposition::Cancelled);
            root.finish().unwrap()
        });
        assert_eq!(receipt.disposition(), InvocationDisposition::Cancelled);
        assert!(receipt.verifies_integrity());
    }

    #[test]
    fn nested_memory_peak_counts_concurrent_parent_and_child_reservations() {
        let clock = VirtualClock::new();
        let (id, accuracy, capability) = identities();
        let receipt = with_cx(|cx| {
            let admission = InvocationAdmitter::new()
                .admit(
                    id,
                    InvocationLimits::new(resources(8), None, accuracy, capability),
                    resources(8),
                )
                .unwrap();
            let mut root = admission.begin(cx, &clock).unwrap();
            {
                let mut parent = root.split_child("parent", resources(8)).unwrap();
                {
                    let mut parent_memory = parent
                        .reserve_memory("parent-memory", MemoryBytes::new(4))
                        .unwrap();
                    {
                        let mut nested = parent_memory
                            .budget()
                            .split_child("nested", resources(4))
                            .unwrap();
                        {
                            let _nested_memory = nested
                                .reserve_memory("nested-memory", MemoryBytes::new(4))
                                .unwrap();
                        }
                        assert_eq!(nested.finish().unwrap(), InvocationDisposition::Completed);
                    }
                }
                assert_eq!(parent.finish().unwrap(), InvocationDisposition::Completed);
            }
            root.finish().unwrap()
        });
        assert_eq!(receipt.memory_peak_bytes(), 8);
        assert_eq!(receipt.memory_requested_bytes(), 8);
        assert_eq!(receipt.memory_released_bytes(), 8);
        assert_eq!(receipt.children()[0].direct_memory_peak_bytes(), 4);
        assert_eq!(receipt.children()[0].memory_peak_bytes(), 8);
        assert!(receipt.verifies_integrity());
    }

    #[test]
    fn semantic_receipt_preserves_cumulative_memory_beyond_u64() {
        let clock = VirtualClock::new();
        let (id, accuracy, capability) = identities();
        let memory_only = resource_vector([0, 0, 0, 0, 1, 0]);
        let receipt = with_cx(|cx| {
            let admission = InvocationAdmitter::new()
                .admit(
                    id,
                    InvocationLimits::new(memory_only, None, accuracy, capability),
                    memory_only,
                )
                .unwrap();
            let mut root = admission.begin(cx, &clock).unwrap();
            {
                let mut child = root.split_child("cumulative-memory", memory_only).unwrap();
                {
                    let _reservation = child
                        .reserve_memory("one-live-byte", MemoryBytes::new(1))
                        .unwrap();
                }
                assert_eq!(child.finish().unwrap(), InvocationDisposition::Completed);
            }
            root.finish().unwrap()
        });
        assert!(receipt.verifies_integrity());

        // This is the exact terminal shape produced by sequentially reusing a
        // one-byte live cap u64::MAX + 1 times. The mutation avoids an
        // infeasible loop while proving the v2 schema neither narrows nor
        // truncates cumulative evidence.
        let cumulative = u128::from(u64::MAX) + 1;
        let mut widened = receipt;
        widened.children[0].memory_requested = cumulative;
        widened.children[0].memory_released = cumulative;
        widened.children[0].root = child_receipt_root(&widened.children[0]);
        widened.memory_requested = cumulative;
        widened.memory_released = cumulative;
        widened.root = invocation_receipt_root(&widened);
        assert_eq!(widened.memory_requested_bytes(), cumulative);
        assert_eq!(widened.memory_released_bytes(), cumulative);
        assert!(widened.verifies_integrity());
    }

    #[test]
    fn semantic_verifier_rejects_rehashed_descendant_memory_peak_underclaim() {
        let clock = VirtualClock::new();
        let (id, accuracy, capability) = identities();
        let receipt = with_cx(|cx| {
            let admission = InvocationAdmitter::new()
                .admit(
                    id,
                    InvocationLimits::new(resources(8), None, accuracy, capability),
                    resources(8),
                )
                .unwrap();
            let mut root = admission.begin(cx, &clock).unwrap();
            {
                let mut parent = root.split_child("parent", resources(8)).unwrap();
                {
                    let mut nested = parent.split_child("nested", resources(8)).unwrap();
                    {
                        let _memory = nested
                            .reserve_memory("nested-memory", MemoryBytes::new(4))
                            .unwrap();
                    }
                    assert_eq!(nested.finish().unwrap(), InvocationDisposition::Completed);
                }
                assert_eq!(parent.finish().unwrap(), InvocationDisposition::Completed);
            }
            root.finish().unwrap()
        });
        assert!(receipt.verifies_integrity());

        let mut forged = receipt;
        forged.children[0].memory_peak = 0;
        forged.children[0].root = child_receipt_root(&forged.children[0]);
        forged.memory_peak = 0;
        forged.root = invocation_receipt_root(&forged);
        assert!(matches!(
            forged.verify_semantics(),
            Err(ReceiptSemanticError::Child {
                invariant: "memory-receipt",
                ..
            })
        ));
    }

    #[test]
    fn semantic_verifier_rejects_rehashed_unmarked_sibling_failures() {
        let clock = VirtualClock::new();
        let (id, accuracy, capability) = identities();
        let receipt = with_cx(|cx| {
            let admission = InvocationAdmitter::new()
                .admit(
                    id,
                    InvocationLimits::new(resources(8), None, accuracy, capability),
                    resources(8),
                )
                .unwrap();
            let mut root = admission.begin(cx, &clock).unwrap();
            for phase in ["first", "second"] {
                let child = root.split_child(phase, resources(4)).unwrap();
                assert_eq!(child.finish().unwrap(), InvocationDisposition::Completed);
            }
            root.finish().unwrap()
        });
        assert!(receipt.verifies_integrity());

        let failure = InvocationError::ExplicitRefusal {
            phase: "forged-siblings",
            reason: hash_domain("test.forged-siblings", b"same-failure"),
        };
        let mut forged = receipt;
        for child in &mut forged.children {
            child.failure = Some(failure.clone());
            child.disposition = InvocationDisposition::Refused;
            child.root = child_receipt_root(child);
        }
        forged.failure = Some(failure);
        // Mark children[0] as the origin so the forgery passes the
        // invocation-level "failure-origin-required" and
        // "failure-origin-matches" gates; the UNMARKED sibling
        // (children[1]: carries the failure, neither origin nor
        // inherited) must then trip "failure-origin-marker".
        forged.failure_origin = Some(forged.children[0].id);
        forged.disposition = InvocationDisposition::Refused;
        forged.root = invocation_receipt_root(&forged);
        assert!(matches!(
            forged.verify_semantics(),
            Err(ReceiptSemanticError::Child {
                invariant: "failure-origin-marker",
                ..
            })
        ));
    }

    #[test]
    fn semantic_verifier_rejects_rehashed_deadline_and_first_fault_forgery() {
        let clock = VirtualClock::new();
        let (id, accuracy, capability) = identities();
        let receipt = with_cx(|cx| {
            let admission = InvocationAdmitter::new()
                .admit(
                    id,
                    InvocationLimits::new(
                        resources(4),
                        Some(Time::from_nanos(10)),
                        accuracy,
                        capability,
                    ),
                    resources(4),
                )
                .unwrap();
            let mut root = admission.begin(cx, &clock).unwrap();
            let mut child = root.split_child("refusal", resources(4)).unwrap();
            assert!(matches!(
                child.charge_work(WorkUnits::new(5)),
                Err(InvocationError::ResourceExceeded { .. })
            ));
            assert_eq!(child.finish().unwrap(), InvocationDisposition::Refused);
            root.finish().unwrap()
        });
        assert!(receipt.verifies_integrity());

        let mut forged_deadline = receipt.clone();
        forged_deadline.last_deadline_observation = Some(Time::from_nanos(10));
        forged_deadline.root = invocation_receipt_root(&forged_deadline);
        assert!(matches!(
            forged_deadline.verify_semantics(),
            Err(ReceiptSemanticError::Invocation {
                invariant: "nondeadline-observation-before-limit"
            })
        ));

        let mut forged_impossible = receipt.clone();
        forged_impossible.failure = Some(InvocationError::ResourceExceeded {
            resource: "work",
            requested: 4,
            available: 4,
        });
        forged_impossible.children[0].failure = forged_impossible.failure.clone();
        forged_impossible.children[0].root = child_receipt_root(&forged_impossible.children[0]);
        forged_impossible.root = invocation_receipt_root(&forged_impossible);
        assert!(matches!(
            forged_impossible.verify_semantics(),
            Err(ReceiptSemanticError::Child {
                invariant: "failure-evidence",
                ..
            })
        ));

        let mut forged_resource_label = receipt.clone();
        forged_resource_label.failure = Some(InvocationError::ResourceExceeded {
            resource: "invented",
            requested: 5,
            available: 4,
        });
        forged_resource_label.children[0].failure = forged_resource_label.failure.clone();
        forged_resource_label.children[0].root =
            child_receipt_root(&forged_resource_label.children[0]);
        forged_resource_label.root = invocation_receipt_root(&forged_resource_label);
        assert!(matches!(
            forged_resource_label.verify_semantics(),
            Err(ReceiptSemanticError::Child {
                invariant: "failure-evidence",
                ..
            })
        ));

        let mut forged_overflow_label = receipt.clone();
        forged_overflow_label.failure = Some(InvocationError::ArithmeticOverflow {
            resource: "invented-overflow",
        });
        forged_overflow_label.children[0].failure = forged_overflow_label.failure.clone();
        forged_overflow_label.children[0].root =
            child_receipt_root(&forged_overflow_label.children[0]);
        forged_overflow_label.root = invocation_receipt_root(&forged_overflow_label);
        assert!(matches!(
            forged_overflow_label.verify_semantics(),
            Err(ReceiptSemanticError::Child {
                invariant: "failure-evidence",
                ..
            })
        ));

        let forged_memory_failure = InvocationError::MemoryRefused {
            what: "forged-memory",
            requested: 5,
            used: 0,
            limit: 4,
        };
        let mut forged_memory_count = receipt.clone();
        forged_memory_count.memory_refusals = 2;
        forged_memory_count.memory_first_refusal = Some(InvocationMemoryRefusal {
            what: "forged-memory",
            requested: 5,
            used: 0,
            limit: 4,
        });
        forged_memory_count.failure = Some(forged_memory_failure.clone());
        forged_memory_count.children[0].failure = Some(forged_memory_failure);
        forged_memory_count.children[0].root = child_receipt_root(&forged_memory_count.children[0]);
        forged_memory_count.root = invocation_receipt_root(&forged_memory_count);
        assert!(matches!(
            forged_memory_count.verify_semantics(),
            Err(ReceiptSemanticError::Invocation {
                invariant: "memory-refusal-evidence"
            })
        ));

        let mut forged_unsealable = receipt.clone();
        forged_unsealable.failure = Some(InvocationError::MemoryReleaseInvariant);
        forged_unsealable.children[0].failure = forged_unsealable.failure.clone();
        forged_unsealable.children[0].root = child_receipt_root(&forged_unsealable.children[0]);
        forged_unsealable.root = invocation_receipt_root(&forged_unsealable);
        assert!(matches!(
            forged_unsealable.verify_semantics(),
            Err(ReceiptSemanticError::Child {
                invariant: "failure-evidence",
                ..
            })
        ));

        let mut forged_root_refusal = receipt.clone();
        forged_root_refusal.failure = Some(InvocationError::ExplicitRefusal {
            phase: "forged-root-refusal",
            reason: hash_domain("test.forged-root-refusal", b"no-child-origin"),
        });
        forged_root_refusal.failure_origin = None;
        forged_root_refusal.disposition = InvocationDisposition::Refused;
        forged_root_refusal.root = invocation_receipt_root(&forged_root_refusal);
        assert!(matches!(
            forged_root_refusal.verify_semantics(),
            Err(ReceiptSemanticError::Invocation {
                invariant: "failure-origin-required"
            })
        ));

        let mut forged_finalizer_origin = receipt.clone();
        let finalizer_failure = InvocationError::ResourceExceeded {
            resource: "finalization-work",
            requested: 2,
            available: 1,
        };
        forged_finalizer_origin.failure = Some(finalizer_failure.clone());
        forged_finalizer_origin.failure_origin = Some(forged_finalizer_origin.children[0].id);
        forged_finalizer_origin.disposition = InvocationDisposition::Refused;
        forged_finalizer_origin.children[0].failure = Some(finalizer_failure);
        forged_finalizer_origin.children[0].failure_inherited = false;
        forged_finalizer_origin.children[0].disposition = InvocationDisposition::Refused;
        forged_finalizer_origin.children[0].root =
            child_receipt_root(&forged_finalizer_origin.children[0]);
        forged_finalizer_origin.root = invocation_receipt_root(&forged_finalizer_origin);
        assert!(matches!(
            forged_finalizer_origin.verify_semantics(),
            Err(ReceiptSemanticError::Invocation {
                invariant: "finalization-failure-origin-kind"
            })
        ));

        let mut forged_fault = receipt;
        // Sever the origin link and mark the child as inherited: with the
        // origin consistent, the child loop reaches the propagation check
        // and the child's foreign failure trips "failure-propagates-to-root"
        // instead of the earlier "failure-origin-matches"/"failure-origin-marker"
        // gates. ResourceExceeded("work") does not require a child origin.
        forged_fault.failure_origin = None;
        forged_fault.children[0].failure = Some(InvocationError::ExplicitRefusal {
            phase: "forged",
            reason: hash_domain("test.forged-fault", b"different"),
        });
        forged_fault.children[0].failure_inherited = true;
        forged_fault.children[0].root = child_receipt_root(&forged_fault.children[0]);
        forged_fault.root = invocation_receipt_root(&forged_fault);
        assert!(matches!(
            forged_fault.verify_semantics(),
            Err(ReceiptSemanticError::Child {
                invariant: "failure-propagates-to-root",
                ..
            })
        ));
    }

    #[test]
    fn semantic_verifier_rejects_descendant_transactional_output_authority() {
        let clock = VirtualClock::new();
        let (id, accuracy, capability) = identities();
        let scientific = resource_vector([0, 0, 0, 0, 0, 1]);
        let required = scientific;
        let receipt = with_cx(|cx| {
            let admission = InvocationAdmitter::new()
                .admit(
                    id,
                    InvocationLimits::new(required, None, accuracy, capability),
                    required,
                )
                .unwrap();
            let mut root = admission.begin(cx, &clock).unwrap();
            let mut outer = root
                .split_finalizable_child(
                    "transaction-owner",
                    scientific,
                    FinalizationResources::default(),
                )
                .unwrap();
            let nested = outer
                .split_child("zero-output-descendant", InvocationResources::default())
                .unwrap();
            assert_eq!(nested.finish().unwrap(), InvocationDisposition::Completed);
            let mut finalizer = outer.begin_finalization();
            finalizer.abort_publication().unwrap();
            finalizer.finish().unwrap();
            root.finish().unwrap()
        });
        assert!(receipt.verifies_integrity());

        let mut unused_grant = receipt;
        let nested_grant = resource_vector([0, 0, 0, 0, 0, 1]);
        let nested = &mut unused_grant.children[1];
        nested.granted = nested_grant;
        nested.returned = nested_grant;
        nested.id = child_id(
            unused_grant.invocation_id,
            unused_grant.plan_binding,
            nested.parent,
            nested.ordinal,
            nested.phase,
            nested_grant,
            FinalizationResources::default(),
        );
        nested.root = child_receipt_root(nested);
        unused_grant.root = invocation_receipt_root(&unused_grant);
        assert!(matches!(
            unused_grant.verify_semantics(),
            Err(ReceiptSemanticError::Child {
                invariant: "finalizer-descendant-output-grant",
                ..
            })
        ));

        let mut retained = unused_grant;
        retained.children[1].returned.output = OutputBytes::new(0);
        retained.children[1].consumed.output = OutputBytes::new(1);
        retained.children[1].direct_consumed.output = OutputBytes::new(1);
        retained.children[1].output_retained = 1;
        retained.children[1].root = child_receipt_root(&retained.children[1]);
        retained.children[0].returned.output = OutputBytes::new(0);
        retained.children[0].consumed.output = OutputBytes::new(1);
        retained.children[0]
            .finalization
            .as_mut()
            .unwrap()
            .scientific_returned
            .output = OutputBytes::new(0);
        retained.children[0].root = child_receipt_root(&retained.children[0]);
        retained.remaining.output = OutputBytes::new(0);
        retained.output_retained = 1;
        retained.root = invocation_receipt_root(&retained);
        assert!(matches!(
            retained.verify_semantics(),
            Err(ReceiptSemanticError::Child {
                invariant: "finalizer-descendant-output",
                ..
            })
        ));
    }

    #[test]
    fn request_winning_immediately_after_swap_rolls_destination_back() {
        let clock = VirtualClock::new();
        let scientific = InvocationResources::new(
            WorkUnits::new(0),
            PollUnits::new(0),
            CostUnits::new(0),
            EvaluationUnits::new(0),
            MemoryBytes::new(0),
            OutputBytes::new(1),
        );
        let finalization = FinalizationResources::new(WorkUnits::new(0), PollUnits::new(1));
        let required = scientific
            .checked_add(finalization.as_invocation_resources())
            .unwrap();
        let (id, accuracy, capability) = identities();
        with_gate_cx(|gate, cx| {
            let admission = InvocationAdmitter::new()
                .admit(
                    id,
                    InvocationLimits::new(required, None, accuracy, capability),
                    required,
                )
                .unwrap();
            let mut root = admission.begin(cx, &clock).unwrap();
            let child = root
                .split_finalizable_child("post-swap-race", scientific, finalization)
                .unwrap();
            let mut finalizer = child.begin_finalization();
            let prepared = finalizer.prepare_publication().unwrap();
            let mut destination = 7_u64;
            let error = finalizer
                .commit_publication_inner(
                    prepared,
                    OutputBytes::new(1),
                    &mut destination,
                    11_u64,
                    || gate.request(),
                )
                .unwrap_err();
            assert_eq!(error.error(), &InvocationError::PublicationForbidden);
            assert_eq!(error.into_parts().1, 11);
            assert_eq!(destination, 7, "the post-swap check must roll back");
            let report = finalizer.finish().unwrap();
            assert_eq!(report.publication(), FinalizationPublication::Aborted);
            assert_eq!(report.disposition(), InvocationDisposition::Cancelled);
            let receipt = root.finish().unwrap();
            assert_eq!(receipt.output_retained_bytes(), 0);
            assert!(receipt.verifies_integrity());
        });
    }

    #[test]
    fn finalizer_mutants_cannot_cross_report_child_or_invocation_boundaries() {
        let clock = VirtualClock::new();
        let scientific = InvocationResources::new(
            WorkUnits::new(4),
            PollUnits::new(0),
            CostUnits::new(0),
            EvaluationUnits::new(0),
            MemoryBytes::new(0),
            OutputBytes::new(0),
        );
        let finalization = FinalizationResources::new(WorkUnits::new(2), PollUnits::new(0));
        let required = scientific
            .checked_add(finalization.as_invocation_resources())
            .unwrap();
        let (id, accuracy, capability) = identities();
        let (report, receipt) = with_cx(|cx| {
            let admission = InvocationAdmitter::new()
                .admit(
                    id,
                    InvocationLimits::new(required, None, accuracy, capability),
                    required,
                )
                .unwrap();
            let mut root = admission.begin(cx, &clock).unwrap();
            let mut child = root
                .split_finalizable_child("mutation-target", scientific, finalization)
                .unwrap();
            child.charge_work(WorkUnits::new(1)).unwrap();
            let mut finalizer = child.begin_finalization();
            finalizer.charge_cleanup_work(WorkUnits::new(2)).unwrap();
            finalizer.abort_publication().unwrap();
            let report = finalizer.finish().unwrap();
            let receipt = root.finish().unwrap();
            (report, receipt)
        });
        assert!(report.verifies_integrity());
        assert!(receipt.verifies_integrity());

        // Mutant 1 removes one cleanup charge while keeping report-local
        // conservation/root self-consistent. Exact child evidence must kill
        // the substitution at join.
        let mut removed_charge = report.clone();
        removed_charge.consumed =
            FinalizationResources::new(WorkUnits::new(1), removed_charge.consumed.polls());
        removed_charge.returned =
            FinalizationResources::new(WorkUnits::new(1), removed_charge.returned.polls());
        removed_charge.root = finalization_report_root(&removed_charge);
        assert!(removed_charge.verifies_integrity());
        assert!(matches!(
            removed_charge.join(&receipt),
            Err(InvocationError::FinalizationReceiptMismatch {
                invariant: "child-finalization-equality"
            })
        ));

        // Mutant 2 relabels child-local evidence as invocation-atomic while
        // keeping the report root self-consistent. The independently retained
        // child partition must reject that authority upgrade.
        let mut upgraded_scope = report.clone();
        upgraded_scope.publication_scope = InvocationPublicationScope::InvocationAtomic;
        upgraded_scope.root = finalization_report_root(&upgraded_scope);
        assert!(upgraded_scope.verifies_integrity());
        assert!(matches!(
            upgraded_scope.join(&receipt),
            Err(InvocationError::FinalizationReceiptMismatch {
                invariant: "child-finalization-equality"
            })
        ));

        // Mutant 3 replaces both finalizer fields and its commitment but does
        // not rewrite the independently recorded scientific partition. The
        // child verifier rejects the shifted direct-spend boundary even after
        // every enclosing hash is recomputed.
        let mut shifted_partition = receipt.clone();
        let child = &mut shifted_partition.children[0];
        let evidence = child.finalization.as_mut().unwrap();
        evidence.consumed = removed_charge.consumed;
        evidence.returned = removed_charge.returned;
        evidence.report_root = removed_charge.root;
        child.root = child_receipt_root(child);
        shifted_partition.root = invocation_receipt_root(&shifted_partition);
        assert!(matches!(
            shifted_partition.verify_semantics(),
            Err(ReceiptSemanticError::Child {
                invariant: "finalizer-direct-partition",
                ..
            })
        ));

        // Mutant 4 claims an availability at the finalizer fault that exceeds
        // the exact reserve bound into the same origin child. Rehashing the
        // report, child, and invocation cannot make that producer state
        // reachable.
        let mut impossible_availability = receipt.clone();
        let failure = InvocationError::ResourceExceeded {
            resource: "finalization-work",
            requested: 101,
            available: 100,
        };
        let mut forged_report = report.clone();
        forged_report.failure = Some(failure.clone());
        forged_report.disposition = InvocationDisposition::Refused;
        forged_report.root = finalization_report_root(&forged_report);
        impossible_availability.failure = Some(failure.clone());
        impossible_availability.failure_origin = Some(impossible_availability.children[0].id);
        impossible_availability.disposition = InvocationDisposition::Refused;
        impossible_availability.children[0].failure = Some(failure);
        impossible_availability.children[0].failure_inherited = false;
        impossible_availability.children[0].disposition = InvocationDisposition::Refused;
        impossible_availability.children[0]
            .finalization
            .as_mut()
            .unwrap()
            .report_root = forged_report.root;
        impossible_availability.children[0].root =
            child_receipt_root(&impossible_availability.children[0]);
        impossible_availability.root = invocation_receipt_root(&impossible_availability);
        assert!(matches!(
            impossible_availability.verify_semantics(),
            Err(ReceiptSemanticError::Child {
                invariant: "failure-available-within-origin-grant",
                ..
            })
        ));

        // Mutant 4 substitutes only a report commitment and rehashes every
        // enclosing receipt. Reconstructing the report from exact child fields
        // rejects the arbitrary commitment.
        let mut wrong_report_root = receipt;
        wrong_report_root.children[0]
            .finalization
            .as_mut()
            .unwrap()
            .report_root = hash_domain("test.finalizer-mutant", b"wrong-report");
        wrong_report_root.children[0].root = child_receipt_root(&wrong_report_root.children[0]);
        wrong_report_root.root = invocation_receipt_root(&wrong_report_root);
        assert!(matches!(
            wrong_report_root.verify_semantics(),
            Err(ReceiptSemanticError::Child {
                invariant: "finalizer-report-root",
                ..
            })
        ));
    }
}
