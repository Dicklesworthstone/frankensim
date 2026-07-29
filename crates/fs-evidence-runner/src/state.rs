//! Pure Runner V2 role/state matrix validation.
//!
//! Validation order is frozen: role/state compatibility, refusal-reason
//! presence, exact base diagnostic, then presented drain-root kind. The result
//! describes a valid cell only; it emits no lifecycle record and proves no
//! drain, execution, or authority fact.

use core::fmt;

use crate::catalog::{
    DiagnosticCodeV2, LifecycleRecordKindV2, NotRunCauseCodeV2, ProofExitV2, RefusedReasonV2,
    StateBearingRecordRoleV2,
};
use crate::identity::{CancelledStopRootV2, DrainedInternalErrorRootV2, TimedOutStopRootV2};

const MAX_NOT_RUN_MANIFEST_CASES_V2: u32 = 256;

/// Presented nominal drain-root kind required by an active terminal state.
///
/// These are descriptive kinds, not root constructors or proof that draining
/// occurred. Nominal wrapper types replace this projection when their owning
/// module is available.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum PresentedDrainRootKindV2 {
    /// A cancellation request was observed and the active work drained.
    CancelledStopRoot,
    /// A timeout was observed and the active work drained.
    TimedOutStopRoot,
    /// A controlled internal error retained a completed drain basis.
    DrainedInternalErrorRoot,
}

impl PresentedDrainRootKindV2 {
    /// Stable descriptive name; this is not a wire discriminant.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::CancelledStopRoot => "cancelled-stop-root",
            Self::TimedOutStopRoot => "timed-out-stop-root",
            Self::DrainedInternalErrorRoot => "drained-internal-error-root",
        }
    }
}

/// Exact typed cause for one contiguous `NotRun` suffix.
///
/// Each variant accepts only its nominal causal-root type. There is no generic
/// digest/root payload, profile-filter case, or catch-all case.
///
/// A checked nominal cancellation root enters only its matching cause
/// variant:
///
/// ```
/// use fs_evidence_runner::identity::CancelledStopRootV2;
/// use fs_evidence_runner::state::NotRunCauseV2;
///
/// let root = CancelledStopRootV2::parse_presented(
///     CancelledStopRootV2::DESCRIPTOR.role(),
///     CancelledStopRootV2::DESCRIPTOR.domain(),
///     &"55".repeat(32),
/// )
/// .unwrap();
/// let cause = NotRunCauseV2::PriorCancelled(root);
/// assert_eq!(cause.name(), "prior-cancelled");
/// ```
///
/// Swapping nominal roots does not type-check:
///
/// ```compile_fail,E0308
/// use fs_evidence_runner::identity::TimedOutStopRootV2;
/// use fs_evidence_runner::state::{NotRunCauseV2, NotRunCauseV2::PriorCancelled};
///
/// fn swap_root(root: TimedOutStopRootV2) -> NotRunCauseV2 {
///     PriorCancelled(root)
/// }
/// ```
///
/// A generic digest cannot replace the nominal causal root:
///
/// ```compile_fail,E0308
/// use fs_evidence_runner::identity::DigestValueV2;
/// use fs_evidence_runner::state::NotRunCauseV2;
///
/// fn erase_nominality(root: DigestValueV2) -> NotRunCauseV2 {
///     NotRunCauseV2::PriorCancelled(root)
/// }
/// ```
///
/// Every variant requires its causal root:
///
/// ```compile_fail,E0308
/// use fs_evidence_runner::state::NotRunCauseV2;
///
/// let _: NotRunCauseV2 = NotRunCauseV2::PriorCancelled;
/// ```
///
/// No profile-filter or generic suppression variant exists:
///
/// ```compile_fail,E0599
/// use fs_evidence_runner::state::NotRunCauseV2;
///
/// let _ = NotRunCauseV2::PriorProfileFilter;
/// ```
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
#[repr(u16)]
pub enum NotRunCauseV2 {
    /// A prior accepted cancellation drained before the remaining suffix.
    PriorCancelled(CancelledStopRootV2) = 1,
    /// A prior accepted timeout drained before the remaining suffix.
    PriorTimedOut(TimedOutStopRootV2) = 2,
    /// A prior controlled internal error drained before the remaining suffix.
    PriorControlledInternalError(DrainedInternalErrorRootV2) = 3,
}

impl NotRunCauseV2 {
    /// Exact closed catalog code carried by this typed cause.
    #[must_use]
    pub const fn cause_code(&self) -> NotRunCauseCodeV2 {
        match self {
            Self::PriorCancelled(_) => NotRunCauseCodeV2::PriorCancelled,
            Self::PriorTimedOut(_) => NotRunCauseCodeV2::PriorTimedOut,
            Self::PriorControlledInternalError(_) => {
                NotRunCauseCodeV2::PriorControlledInternalError
            }
        }
    }

    /// Exact unsigned 16-bit wire code.
    #[must_use]
    pub const fn code(&self) -> u16 {
        self.cause_code().code()
    }

    /// Stable lowercase cause name.
    #[must_use]
    pub const fn name(&self) -> &'static str {
        self.cause_code().name()
    }

    /// Cancellation root, present only for [`Self::PriorCancelled`].
    #[must_use]
    pub const fn prior_cancelled_root(&self) -> Option<&CancelledStopRootV2> {
        match self {
            Self::PriorCancelled(root) => Some(root),
            Self::PriorTimedOut(_) | Self::PriorControlledInternalError(_) => None,
        }
    }

    /// Timeout root, present only for [`Self::PriorTimedOut`].
    #[must_use]
    pub const fn prior_timed_out_root(&self) -> Option<&TimedOutStopRootV2> {
        match self {
            Self::PriorTimedOut(root) => Some(root),
            Self::PriorCancelled(_) | Self::PriorControlledInternalError(_) => None,
        }
    }

    /// Controlled-error root, present only for
    /// [`Self::PriorControlledInternalError`].
    #[must_use]
    pub const fn prior_controlled_internal_error_root(
        &self,
    ) -> Option<&DrainedInternalErrorRootV2> {
        match self {
            Self::PriorControlledInternalError(root) => Some(root),
            Self::PriorCancelled(_) | Self::PriorTimedOut(_) => None,
        }
    }
}

/// Validated basis for the first slot in one contiguous `NotRun` suffix.
///
/// Its only semantic fields are the exact typed cause and the lowest remaining
/// manifest ordinal. The manifest count is constructor context, not duplicated
/// semantic data.
///
/// Fields are private, so validation cannot be bypassed with a struct literal:
///
/// ```
/// use fs_evidence_runner::identity::CancelledStopRootV2;
/// use fs_evidence_runner::state::{NotRunBasisV2, NotRunCauseV2};
///
/// let root = CancelledStopRootV2::parse_presented(
///     CancelledStopRootV2::DESCRIPTOR.role(),
///     CancelledStopRootV2::DESCRIPTOR.domain(),
///     &"66".repeat(32),
/// )
/// .unwrap();
/// let basis = NotRunBasisV2::new(NotRunCauseV2::PriorCancelled(root), 2, 5).unwrap();
/// assert_eq!(basis.lowest_remaining_manifest_ordinal(), 2);
/// assert_eq!(basis.remaining_case_count(5).unwrap(), 3);
/// ```
///
/// ```compile_fail,E0451
/// use fs_evidence_runner::state::{NotRunBasisV2, NotRunCauseV2};
///
/// fn forge(cause: NotRunCauseV2) -> NotRunBasisV2 {
///     NotRunBasisV2 {
///         cause,
///         lowest_remaining_manifest_ordinal: 0,
///     }
/// }
/// ```
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct NotRunBasisV2 {
    cause: NotRunCauseV2,
    lowest_remaining_manifest_ordinal: u32,
}

impl NotRunBasisV2 {
    /// Validate a typed cause and the lowest ordinal of a nonempty remaining
    /// manifest suffix.
    ///
    /// `ordered_case_count` is checked against the frozen 256-case base
    /// ceiling before ordinal arithmetic. No row vector is allocated or
    /// inspected.
    ///
    /// # Errors
    ///
    /// Returns [`NotRunBasisErrorV2`] when the manifest is empty, exceeds the
    /// frozen case ceiling, or does not contain the supplied ordinal.
    pub const fn new(
        cause: NotRunCauseV2,
        lowest_remaining_manifest_ordinal: u32,
        ordered_case_count: u32,
    ) -> Result<Self, NotRunBasisErrorV2> {
        match validate_not_run_ordinal_v2(lowest_remaining_manifest_ordinal, ordered_case_count) {
            Ok(()) => {}
            Err(error) => return Err(error),
        }
        Ok(Self {
            cause,
            lowest_remaining_manifest_ordinal,
        })
    }

    /// Exact typed cause for every slot in the remaining suffix.
    #[must_use]
    pub const fn cause(&self) -> &NotRunCauseV2 {
        &self.cause
    }

    /// Lowest manifest ordinal suppressed by this cause.
    #[must_use]
    pub const fn lowest_remaining_manifest_ordinal(&self) -> u32 {
        self.lowest_remaining_manifest_ordinal
    }

    /// Terminal state required for every slot covered by this basis.
    #[must_use]
    pub const fn state(&self) -> ProofExitV2 {
        ProofExitV2::NotRun
    }

    /// Exact base diagnostic required for every slot covered by this basis.
    #[must_use]
    pub const fn diagnostic(&self) -> DiagnosticCodeV2 {
        DiagnosticCodeV2::RunnerNotRun
    }

    /// Compute the number of remaining manifest slots without allocating rows.
    ///
    /// # Errors
    ///
    /// Revalidates the supplied manifest count against the frozen ceiling and
    /// this basis's lowest ordinal before subtraction.
    pub const fn remaining_case_count(
        &self,
        ordered_case_count: u32,
    ) -> Result<u32, NotRunBasisErrorV2> {
        match validate_not_run_ordinal_v2(
            self.lowest_remaining_manifest_ordinal,
            ordered_case_count,
        ) {
            Ok(()) => {}
            Err(error) => return Err(error),
        }
        Ok(ordered_case_count - self.lowest_remaining_manifest_ordinal)
    }
}

/// Deterministic refusal from [`NotRunBasisV2`] ordinal validation.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum NotRunBasisErrorV2 {
    /// A zero-case manifest has no remaining manifest ordinal.
    EmptyManifest,
    /// The manifest case count exceeds the frozen base ceiling.
    ManifestCaseCountExceedsMaximum {
        /// Supplied manifest case count.
        observed: u32,
        /// Frozen maximum manifest case count.
        maximum: u32,
    },
    /// The alleged lowest remaining ordinal is outside the manifest.
    LowestRemainingOrdinalOutOfRange {
        /// Supplied ordinal.
        observed: u32,
        /// Supplied nonzero, in-limit manifest case count.
        ordered_case_count: u32,
    },
}

impl fmt::Display for NotRunBasisErrorV2 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyManifest => write!(f, "a zero-case manifest has no NotRun suffix"),
            Self::ManifestCaseCountExceedsMaximum { observed, maximum } => write!(
                f,
                "manifest case count {observed} exceeds frozen maximum {maximum}"
            ),
            Self::LowestRemainingOrdinalOutOfRange {
                observed,
                ordered_case_count,
            } => write!(
                f,
                "lowest remaining ordinal {observed} is outside case count {ordered_case_count}"
            ),
        }
    }
}

impl std::error::Error for NotRunBasisErrorV2 {}

const fn validate_not_run_ordinal_v2(
    lowest_remaining_manifest_ordinal: u32,
    ordered_case_count: u32,
) -> Result<(), NotRunBasisErrorV2> {
    if ordered_case_count == 0 {
        return Err(NotRunBasisErrorV2::EmptyManifest);
    }
    if ordered_case_count > MAX_NOT_RUN_MANIFEST_CASES_V2 {
        return Err(NotRunBasisErrorV2::ManifestCaseCountExceedsMaximum {
            observed: ordered_case_count,
            maximum: MAX_NOT_RUN_MANIFEST_CASES_V2,
        });
    }
    if lowest_remaining_manifest_ordinal >= ordered_case_count {
        return Err(NotRunBasisErrorV2::LowestRemainingOrdinalOutOfRange {
            observed: lowest_remaining_manifest_ordinal,
            ordered_case_count,
        });
    }
    Ok(())
}

/// One state cell that passed the complete closed matrix.
///
/// Fields are private so callers cannot manufacture a validated terminal
/// combination without [`validate_state_cell_v2`].
///
/// ```
/// use fs_evidence_runner::catalog::{ProofExitV2, StateBearingRecordRoleV2};
/// use fs_evidence_runner::state::validate_state_cell_v2;
///
/// let cell = validate_state_cell_v2(
///     StateBearingRecordRoleV2::ExecutedCaseTerminal,
///     ProofExitV2::Pass,
///     None,
///     None,
///     None,
/// )
/// .unwrap();
/// assert_eq!(cell.state(), ProofExitV2::Pass);
/// ```
///
/// An unvalidated candidate cannot convert directly into a validated terminal
/// cell:
///
/// ```compile_fail,E0277
/// use fs_evidence_runner::catalog::{ProofExitV2, StateBearingRecordRoleV2};
/// use fs_evidence_runner::state::{StateValidationInputV2, ValidatedStateCellV2};
///
/// let candidate = StateValidationInputV2::new(
///     StateBearingRecordRoleV2::ExecutedCaseTerminal,
///     ProofExitV2::Pass,
///     None,
///     None,
///     None,
/// );
/// let _terminal: ValidatedStateCellV2 = candidate.into();
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ValidatedStateCellV2 {
    role: StateBearingRecordRoleV2,
    state: ProofExitV2,
    refused_reason: Option<RefusedReasonV2>,
    diagnostic: Option<DiagnosticCodeV2>,
    drain_basis: Option<PresentedDrainRootKindV2>,
}

/// Complete unvalidated input to the closed state matrix.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StateValidationInputV2 {
    role: StateBearingRecordRoleV2,
    state: ProofExitV2,
    refused_reason: Option<RefusedReasonV2>,
    diagnostic: Option<DiagnosticCodeV2>,
    drain_basis: Option<PresentedDrainRootKindV2>,
}

impl StateValidationInputV2 {
    /// Assemble one candidate matrix cell without validating it.
    #[must_use]
    pub const fn new(
        role: StateBearingRecordRoleV2,
        state: ProofExitV2,
        refused_reason: Option<RefusedReasonV2>,
        diagnostic: Option<DiagnosticCodeV2>,
        drain_basis: Option<PresentedDrainRootKindV2>,
    ) -> Self {
        Self {
            role,
            state,
            refused_reason,
            diagnostic,
            drain_basis,
        }
    }

    /// Candidate state-bearing role.
    #[must_use]
    pub const fn role(self) -> StateBearingRecordRoleV2 {
        self.role
    }

    /// Candidate terminal state.
    #[must_use]
    pub const fn state(self) -> ProofExitV2 {
        self.state
    }

    /// Candidate refusal reason.
    #[must_use]
    pub const fn refused_reason(self) -> Option<RefusedReasonV2> {
        self.refused_reason
    }

    /// Candidate base diagnostic.
    #[must_use]
    pub const fn diagnostic(self) -> Option<DiagnosticCodeV2> {
        self.diagnostic
    }

    /// Candidate presented drain-root kind.
    #[must_use]
    pub const fn drain_basis(self) -> Option<PresentedDrainRootKindV2> {
        self.drain_basis
    }
}

impl ValidatedStateCellV2 {
    /// State-bearing role admitted by the matrix.
    #[must_use]
    pub const fn role(self) -> StateBearingRecordRoleV2 {
        self.role
    }

    /// Terminal state admitted for the role.
    #[must_use]
    pub const fn state(self) -> ProofExitV2 {
        self.state
    }

    /// Refusal reason, present exactly for [`ProofExitV2::Refused`].
    #[must_use]
    pub const fn refused_reason(self) -> Option<RefusedReasonV2> {
        self.refused_reason
    }

    /// Exact base diagnostic, absent exactly for [`ProofExitV2::Pass`].
    #[must_use]
    pub const fn diagnostic(self) -> Option<DiagnosticCodeV2> {
        self.diagnostic
    }

    /// Required presented drain-root kind for active stop states.
    #[must_use]
    pub const fn drain_basis(self) -> Option<PresentedDrainRootKindV2> {
        self.drain_basis
    }

    /// Lifecycle record kind implied by the role.
    ///
    /// Pre-run diagnostics deliberately return `None`: even pre-run
    /// cancellation, timeout, and controlled-internal-error cells have no
    /// lifecycle-record basis.
    #[must_use]
    pub const fn lifecycle_record_kind(self) -> Option<LifecycleRecordKindV2> {
        lifecycle_record_kind_for_role(self.role)
    }
}

/// Deterministic refusal from the closed role/state matrix.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StateValidationErrorV2 {
    /// The state is not legal for the selected state-bearing role.
    StateNotAllowedForRole {
        /// Selected role.
        role: StateBearingRecordRoleV2,
        /// Disallowed state.
        state: ProofExitV2,
    },
    /// `Refused` omitted its one required closed refusal reason.
    MissingRefusedReason,
    /// A non-refused state attempted to carry a refusal reason.
    UnexpectedRefusedReason {
        /// State that forbids the reason.
        state: ProofExitV2,
        /// Unexpected reason.
        observed: RefusedReasonV2,
    },
    /// A non-pass state omitted its exact base diagnostic.
    MissingDiagnostic {
        /// Required diagnostic.
        expected: DiagnosticCodeV2,
    },
    /// `Pass` attempted to carry a diagnostic.
    UnexpectedDiagnostic {
        /// Unexpected diagnostic.
        observed: DiagnosticCodeV2,
    },
    /// A non-pass state carried a base diagnostic belonging to another state.
    WrongDiagnostic {
        /// Required diagnostic.
        expected: DiagnosticCodeV2,
        /// Supplied diagnostic.
        observed: DiagnosticCodeV2,
    },
    /// An active stop state omitted its typed presented drain-root kind.
    MissingDrainBasis {
        /// Required presented root kind.
        expected: PresentedDrainRootKindV2,
    },
    /// A state with no active drain requirement attempted to carry a basis.
    UnexpectedDrainBasis {
        /// State that forbids an active drain basis in this role.
        state: ProofExitV2,
        /// Unexpected presented root kind.
        observed: PresentedDrainRootKindV2,
    },
    /// An active stop state carried the wrong nominal drain-root kind.
    WrongDrainBasis {
        /// Required presented root kind.
        expected: PresentedDrainRootKindV2,
        /// Supplied presented root kind.
        observed: PresentedDrainRootKindV2,
    },
}

impl fmt::Display for StateValidationErrorV2 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::StateNotAllowedForRole { role, state } => write!(
                f,
                "state {} is not allowed for role {}",
                state.name(),
                role.name()
            ),
            Self::MissingRefusedReason => {
                write!(f, "state refused requires exactly one refusal reason")
            }
            Self::UnexpectedRefusedReason { state, observed } => write!(
                f,
                "state {} forbids refusal reason {}",
                state.name(),
                observed.name()
            ),
            Self::MissingDiagnostic { expected } => {
                write!(f, "state requires diagnostic {}", expected.name())
            }
            Self::UnexpectedDiagnostic { observed } => {
                write!(f, "state pass forbids diagnostic {}", observed.name())
            }
            Self::WrongDiagnostic { expected, observed } => write!(
                f,
                "wrong diagnostic {}; expected {}",
                observed.name(),
                expected.name()
            ),
            Self::MissingDrainBasis { expected } => {
                write!(f, "active stop state requires {}", expected.name())
            }
            Self::UnexpectedDrainBasis { state, observed } => write!(
                f,
                "state {} forbids drain basis {} in this role",
                state.name(),
                observed.name()
            ),
            Self::WrongDrainBasis { expected, observed } => write!(
                f,
                "wrong drain basis {}; expected {}",
                observed.name(),
                expected.name()
            ),
        }
    }
}

impl std::error::Error for StateValidationErrorV2 {}

/// Validate one complete role/state/reason/diagnostic/drain-basis cell.
///
/// Error precedence is part of the API and is evaluated in this order:
///
/// 1. role/state compatibility;
/// 2. refusal-reason presence or absence;
/// 3. exact diagnostic presence and code;
/// 4. exact active drain-root kind.
///
/// # Errors
///
/// Returns [`StateValidationErrorV2`] at the first failed rule.
pub const fn validate_state_cell_v2(
    role: StateBearingRecordRoleV2,
    state: ProofExitV2,
    refused_reason: Option<RefusedReasonV2>,
    diagnostic: Option<DiagnosticCodeV2>,
    drain_basis: Option<PresentedDrainRootKindV2>,
) -> Result<ValidatedStateCellV2, StateValidationErrorV2> {
    if !role_allows_state(role, state) {
        return Err(StateValidationErrorV2::StateNotAllowedForRole { role, state });
    }

    match (state, refused_reason) {
        (ProofExitV2::Refused, None) => {
            return Err(StateValidationErrorV2::MissingRefusedReason);
        }
        (ProofExitV2::Refused, Some(_)) | (_, None) => {}
        (_, Some(observed)) => {
            return Err(StateValidationErrorV2::UnexpectedRefusedReason { state, observed });
        }
    }

    match (diagnostic_for_state(state), diagnostic) {
        (None, Some(observed)) => {
            return Err(StateValidationErrorV2::UnexpectedDiagnostic { observed });
        }
        (Some(expected), None) => {
            return Err(StateValidationErrorV2::MissingDiagnostic { expected });
        }
        (Some(expected), Some(observed)) if expected.code() != observed.code() => {
            return Err(StateValidationErrorV2::WrongDiagnostic { expected, observed });
        }
        (None, None) | (Some(_), Some(_)) => {}
    }

    match (drain_basis_for(role, state), drain_basis) {
        (None, Some(observed)) => {
            return Err(StateValidationErrorV2::UnexpectedDrainBasis { state, observed });
        }
        (Some(expected), None) => {
            return Err(StateValidationErrorV2::MissingDrainBasis { expected });
        }
        (Some(expected), Some(observed)) if !same_drain_kind(expected, observed) => {
            return Err(StateValidationErrorV2::WrongDrainBasis { expected, observed });
        }
        (None, None) | (Some(_), Some(_)) => {}
    }

    Ok(ValidatedStateCellV2 {
        role,
        state,
        refused_reason,
        diagnostic,
        drain_basis,
    })
}

const fn same_drain_kind(left: PresentedDrainRootKindV2, right: PresentedDrainRootKindV2) -> bool {
    matches!(
        (left, right),
        (
            PresentedDrainRootKindV2::CancelledStopRoot,
            PresentedDrainRootKindV2::CancelledStopRoot
        ) | (
            PresentedDrainRootKindV2::TimedOutStopRoot,
            PresentedDrainRootKindV2::TimedOutStopRoot
        ) | (
            PresentedDrainRootKindV2::DrainedInternalErrorRoot,
            PresentedDrainRootKindV2::DrainedInternalErrorRoot
        )
    )
}

/// Validate one assembled [`StateValidationInputV2`].
///
/// # Errors
///
/// Returns [`StateValidationErrorV2`] at the first failed matrix rule.
pub const fn validate_state_v2(
    input: StateValidationInputV2,
) -> Result<ValidatedStateCellV2, StateValidationErrorV2> {
    validate_state_cell_v2(
        input.role,
        input.state,
        input.refused_reason,
        input.diagnostic,
        input.drain_basis,
    )
}

/// Whether a state belongs to the role's frozen terminal-state set.
#[must_use]
pub const fn role_allows_state(role: StateBearingRecordRoleV2, state: ProofExitV2) -> bool {
    match role {
        StateBearingRecordRoleV2::PreRunDiagnostic => matches!(
            state,
            ProofExitV2::Usage
                | ProofExitV2::Refused
                | ProofExitV2::NoData
                | ProofExitV2::Stale
                | ProofExitV2::EnvironmentInvalid
                | ProofExitV2::Blocked
                | ProofExitV2::Unsupported
                | ProofExitV2::Cancelled
                | ProofExitV2::TimedOut
                | ProofExitV2::InternalError
        ),
        StateBearingRecordRoleV2::ExecutedCaseTerminal | StateBearingRecordRoleV2::RunTerminal => {
            !matches!(state, ProofExitV2::Usage | ProofExitV2::NotRun)
        }
        StateBearingRecordRoleV2::SuppressedCaseTerminal => {
            matches!(state, ProofExitV2::NotRun)
        }
    }
}

/// Exact base diagnostic required by a state.
///
/// Only `Pass` maps to absence.
#[must_use]
pub const fn diagnostic_for_state(state: ProofExitV2) -> Option<DiagnosticCodeV2> {
    match state {
        ProofExitV2::Pass => None,
        ProofExitV2::Failed => Some(DiagnosticCodeV2::CaseConformanceMismatch),
        ProofExitV2::Refused => Some(DiagnosticCodeV2::RunnerRefused),
        ProofExitV2::NoData => Some(DiagnosticCodeV2::RunnerNoData),
        ProofExitV2::Stale => Some(DiagnosticCodeV2::RunnerStale),
        ProofExitV2::EnvironmentInvalid => Some(DiagnosticCodeV2::RunnerEnvironmentInvalid),
        ProofExitV2::Blocked => Some(DiagnosticCodeV2::RunnerBlocked),
        ProofExitV2::Unsupported => Some(DiagnosticCodeV2::RunnerUnsupported),
        ProofExitV2::NotRun => Some(DiagnosticCodeV2::RunnerNotRun),
        ProofExitV2::Cancelled => Some(DiagnosticCodeV2::RunnerCancelled),
        ProofExitV2::TimedOut => Some(DiagnosticCodeV2::RunnerTimedOut),
        ProofExitV2::Usage => Some(DiagnosticCodeV2::RunnerUsage),
        ProofExitV2::InternalError => Some(DiagnosticCodeV2::RunnerInternalError),
    }
}

/// Active presented drain-root kind required by one role/state pair.
///
/// Pre-run stop diagnostics intentionally return `None`; they have not
/// entered the lifecycle and cannot present an active lifecycle drain root.
#[must_use]
pub const fn drain_basis_for(
    role: StateBearingRecordRoleV2,
    state: ProofExitV2,
) -> Option<PresentedDrainRootKindV2> {
    if matches!(role, StateBearingRecordRoleV2::PreRunDiagnostic) {
        return None;
    }

    match state {
        ProofExitV2::Cancelled => Some(PresentedDrainRootKindV2::CancelledStopRoot),
        ProofExitV2::TimedOut => Some(PresentedDrainRootKindV2::TimedOutStopRoot),
        ProofExitV2::InternalError => Some(PresentedDrainRootKindV2::DrainedInternalErrorRoot),
        _ => None,
    }
}

/// Lifecycle-record kind implied by a state-bearing role.
#[must_use]
pub const fn lifecycle_record_kind_for_role(
    role: StateBearingRecordRoleV2,
) -> Option<LifecycleRecordKindV2> {
    match role {
        StateBearingRecordRoleV2::PreRunDiagnostic => None,
        StateBearingRecordRoleV2::ExecutedCaseTerminal
        | StateBearingRecordRoleV2::SuppressedCaseTerminal => {
            Some(LifecycleRecordKindV2::CaseTerminal)
        }
        StateBearingRecordRoleV2::RunTerminal => Some(LifecycleRecordKindV2::RunTerminal),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::catalog::RunProfileV2;
    use crate::limits::RunnerLimitsCandidateV2;

    const STATES: [ProofExitV2; 13] = [
        ProofExitV2::Pass,
        ProofExitV2::Failed,
        ProofExitV2::Refused,
        ProofExitV2::NoData,
        ProofExitV2::Stale,
        ProofExitV2::EnvironmentInvalid,
        ProofExitV2::Blocked,
        ProofExitV2::Unsupported,
        ProofExitV2::NotRun,
        ProofExitV2::Cancelled,
        ProofExitV2::TimedOut,
        ProofExitV2::Usage,
        ProofExitV2::InternalError,
    ];

    const ROLES: [StateBearingRecordRoleV2; 4] = [
        StateBearingRecordRoleV2::PreRunDiagnostic,
        StateBearingRecordRoleV2::ExecutedCaseTerminal,
        StateBearingRecordRoleV2::SuppressedCaseTerminal,
        StateBearingRecordRoleV2::RunTerminal,
    ];

    const REASONS: [Option<RefusedReasonV2>; 12] = [
        None,
        Some(RefusedReasonV2::InvalidEvidence),
        Some(RefusedReasonV2::NonCanonicalEvidence),
        Some(RefusedReasonV2::EvidenceIdentityMismatch),
        Some(RefusedReasonV2::EvidenceTampered),
        Some(RefusedReasonV2::LimitExceeded),
        Some(RefusedReasonV2::UnsafeArtifactPlacement),
        Some(RefusedReasonV2::ArtifactCollision),
        Some(RefusedReasonV2::LifecycleViolation),
        Some(RefusedReasonV2::PolicyRefused),
        Some(RefusedReasonV2::AuthorityBoundaryViolation),
        Some(RefusedReasonV2::MigrationRefused),
    ];

    const DIAGNOSTICS: [Option<DiagnosticCodeV2>; 13] = [
        None,
        Some(DiagnosticCodeV2::CaseConformanceMismatch),
        Some(DiagnosticCodeV2::RunnerNotRun),
        Some(DiagnosticCodeV2::RunnerRefused),
        Some(DiagnosticCodeV2::RunnerNoData),
        Some(DiagnosticCodeV2::RunnerStale),
        Some(DiagnosticCodeV2::RunnerEnvironmentInvalid),
        Some(DiagnosticCodeV2::RunnerBlocked),
        Some(DiagnosticCodeV2::RunnerUnsupported),
        Some(DiagnosticCodeV2::RunnerCancelled),
        Some(DiagnosticCodeV2::RunnerTimedOut),
        Some(DiagnosticCodeV2::RunnerUsage),
        Some(DiagnosticCodeV2::RunnerInternalError),
    ];

    const DRAIN_BASES: [Option<PresentedDrainRootKindV2>; 4] = [
        None,
        Some(PresentedDrainRootKindV2::CancelledStopRoot),
        Some(PresentedDrainRootKindV2::TimedOutStopRoot),
        Some(PresentedDrainRootKindV2::DrainedInternalErrorRoot),
    ];

    const ZERO_DIGEST_HEX: &str =
        "0000000000000000000000000000000000000000000000000000000000000000";
    const ONE_DIGEST_HEX: &str = "0101010101010101010101010101010101010101010101010101010101010101";
    const TWO_DIGEST_HEX: &str = "0202020202020202020202020202020202020202020202020202020202020202";

    fn cancelled_root() -> CancelledStopRootV2 {
        CancelledStopRootV2::parse_presented(
            CancelledStopRootV2::DESCRIPTOR.role(),
            CancelledStopRootV2::DESCRIPTOR.domain(),
            ZERO_DIGEST_HEX,
        )
        .expect("literal cancellation root must parse")
    }

    fn timed_out_root() -> TimedOutStopRootV2 {
        TimedOutStopRootV2::parse_presented(
            TimedOutStopRootV2::DESCRIPTOR.role(),
            TimedOutStopRootV2::DESCRIPTOR.domain(),
            ONE_DIGEST_HEX,
        )
        .expect("literal timeout root must parse")
    }

    fn controlled_error_root() -> DrainedInternalErrorRootV2 {
        DrainedInternalErrorRootV2::parse_presented(
            DrainedInternalErrorRootV2::DESCRIPTOR.role(),
            DrainedInternalErrorRootV2::DESCRIPTOR.domain(),
            TWO_DIGEST_HEX,
        )
        .expect("literal controlled-error root must parse")
    }

    #[test]
    fn not_run_causes_have_exact_codes_names_and_nominal_accessors() {
        let cancelled_root = cancelled_root();
        let timed_out_root = timed_out_root();
        let controlled_error_root = controlled_error_root();
        let cancelled = NotRunCauseV2::PriorCancelled(cancelled_root.clone());
        let timed_out = NotRunCauseV2::PriorTimedOut(timed_out_root.clone());
        let controlled = NotRunCauseV2::PriorControlledInternalError(controlled_error_root.clone());

        let expected = [
            (
                &cancelled,
                NotRunCauseCodeV2::PriorCancelled,
                1,
                "prior-cancelled",
            ),
            (
                &timed_out,
                NotRunCauseCodeV2::PriorTimedOut,
                2,
                "prior-timed-out",
            ),
            (
                &controlled,
                NotRunCauseCodeV2::PriorControlledInternalError,
                3,
                "prior-controlled-internal-error",
            ),
        ];
        for (cause, expected_code, expected_wire_code, expected_name) in expected {
            assert_eq!(cause.cause_code(), expected_code);
            assert_eq!(cause.code(), expected_wire_code);
            assert_eq!(cause.name(), expected_name);
        }

        assert_eq!(
            cancelled.prior_cancelled_root(),
            Some(&cancelled_root),
            "cancellation retains only its nominal root"
        );
        assert_eq!(cancelled.prior_timed_out_root(), None);
        assert_eq!(cancelled.prior_controlled_internal_error_root(), None);

        assert_eq!(timed_out.prior_cancelled_root(), None);
        assert_eq!(timed_out.prior_timed_out_root(), Some(&timed_out_root));
        assert_eq!(timed_out.prior_controlled_internal_error_root(), None);

        assert_eq!(controlled.prior_cancelled_root(), None);
        assert_eq!(controlled.prior_timed_out_root(), None);
        assert_eq!(
            controlled.prior_controlled_internal_error_root(),
            Some(&controlled_error_root)
        );

        let _: fn(CancelledStopRootV2) -> NotRunCauseV2 = NotRunCauseV2::PriorCancelled;
        let _: fn(TimedOutStopRootV2) -> NotRunCauseV2 = NotRunCauseV2::PriorTimedOut;
        let _: fn(DrainedInternalErrorRootV2) -> NotRunCauseV2 =
            NotRunCauseV2::PriorControlledInternalError;

        assert_eq!(
            NotRunCauseCodeV2::ALL,
            [
                NotRunCauseCodeV2::PriorCancelled,
                NotRunCauseCodeV2::PriorTimedOut,
                NotRunCauseCodeV2::PriorControlledInternalError,
            ]
        );
        for unknown in [0, 4, u16::MAX] {
            assert!(
                NotRunCauseCodeV2::from_code(unknown).is_err(),
                "unknown cause code {unknown} must refuse"
            );
        }
    }

    #[test]
    fn not_run_basis_validates_manifest_boundaries_and_exact_diagnostic() {
        assert_eq!(
            MAX_NOT_RUN_MANIFEST_CASES_V2,
            RunnerLimitsCandidateV2::base(RunProfileV2::Smoke).invocation_cases,
            "local NotRun ceiling drifted from Smoke base limits"
        );
        assert_eq!(
            MAX_NOT_RUN_MANIFEST_CASES_V2,
            RunnerLimitsCandidateV2::base(RunProfileV2::Full).invocation_cases,
            "local NotRun ceiling drifted from Full base limits"
        );

        let first = NotRunBasisV2::new(NotRunCauseV2::PriorCancelled(cancelled_root()), 0, 1)
            .expect("ordinal zero is the first slot in a one-case manifest");
        assert_eq!(first.lowest_remaining_manifest_ordinal(), 0);
        assert_eq!(first.remaining_case_count(1), Ok(1));
        assert_eq!(first.state(), ProofExitV2::NotRun);
        assert_eq!(first.diagnostic(), DiagnosticCodeV2::RunnerNotRun);
        assert_eq!(
            first.cause().cause_code(),
            NotRunCauseCodeV2::PriorCancelled
        );

        let last = NotRunBasisV2::new(
            NotRunCauseV2::PriorTimedOut(timed_out_root()),
            255,
            MAX_NOT_RUN_MANIFEST_CASES_V2,
        )
        .expect("ordinal 255 is the last slot in the maximum-size manifest");
        assert_eq!(last.lowest_remaining_manifest_ordinal(), 255);
        assert_eq!(
            last.remaining_case_count(MAX_NOT_RUN_MANIFEST_CASES_V2),
            Ok(1)
        );

        assert_eq!(
            NotRunBasisV2::new(
                NotRunCauseV2::PriorControlledInternalError(controlled_error_root()),
                0,
                0,
            ),
            Err(NotRunBasisErrorV2::EmptyManifest)
        );
        assert_eq!(
            NotRunBasisV2::new(
                NotRunCauseV2::PriorControlledInternalError(controlled_error_root()),
                u32::MAX,
                MAX_NOT_RUN_MANIFEST_CASES_V2 + 1,
            ),
            Err(NotRunBasisErrorV2::ManifestCaseCountExceedsMaximum {
                observed: 257,
                maximum: 256,
            }),
            "count-ceiling error precedes ordinal validation"
        );
        assert_eq!(
            NotRunBasisV2::new(
                NotRunCauseV2::PriorControlledInternalError(controlled_error_root()),
                MAX_NOT_RUN_MANIFEST_CASES_V2,
                MAX_NOT_RUN_MANIFEST_CASES_V2,
            ),
            Err(NotRunBasisErrorV2::LowestRemainingOrdinalOutOfRange {
                observed: 256,
                ordered_case_count: 256,
            })
        );
    }

    #[test]
    fn not_run_remaining_suffix_arithmetic_is_exhaustive_and_allocation_free() {
        let cause = NotRunCauseV2::PriorCancelled(cancelled_root());
        for lowest in 0..MAX_NOT_RUN_MANIFEST_CASES_V2 {
            let basis = NotRunBasisV2::new(cause.clone(), lowest, MAX_NOT_RUN_MANIFEST_CASES_V2)
                .expect("every in-range lowest ordinal must validate");
            assert_eq!(
                basis.remaining_case_count(MAX_NOT_RUN_MANIFEST_CASES_V2),
                Ok(MAX_NOT_RUN_MANIFEST_CASES_V2 - lowest),
                "remaining suffix arithmetic drift at ordinal {lowest}"
            );
        }

        let zero_family_row_basis = NotRunBasisV2::new(cause, 0, 256)
            .expect("suppression is independent of family-row allocation");
        assert_eq!(zero_family_row_basis.remaining_case_count(256), Ok(256));
        assert_eq!(
            zero_family_row_basis.remaining_case_count(0),
            Err(NotRunBasisErrorV2::EmptyManifest)
        );
        assert_eq!(
            zero_family_row_basis.remaining_case_count(257),
            Err(NotRunBasisErrorV2::ManifestCaseCountExceedsMaximum {
                observed: 257,
                maximum: 256,
            })
        );
    }

    fn oracle_role_allows(role: StateBearingRecordRoleV2, state: ProofExitV2) -> bool {
        match role {
            StateBearingRecordRoleV2::PreRunDiagnostic => matches!(
                state,
                ProofExitV2::Usage
                    | ProofExitV2::Refused
                    | ProofExitV2::NoData
                    | ProofExitV2::Stale
                    | ProofExitV2::EnvironmentInvalid
                    | ProofExitV2::Blocked
                    | ProofExitV2::Unsupported
                    | ProofExitV2::Cancelled
                    | ProofExitV2::TimedOut
                    | ProofExitV2::InternalError
            ),
            StateBearingRecordRoleV2::ExecutedCaseTerminal
            | StateBearingRecordRoleV2::RunTerminal => {
                !matches!(state, ProofExitV2::Usage | ProofExitV2::NotRun)
            }
            StateBearingRecordRoleV2::SuppressedCaseTerminal => state == ProofExitV2::NotRun,
        }
    }

    fn oracle_diagnostic(state: ProofExitV2) -> Option<DiagnosticCodeV2> {
        match state {
            ProofExitV2::Pass => None,
            ProofExitV2::Failed => Some(DiagnosticCodeV2::CaseConformanceMismatch),
            ProofExitV2::Refused => Some(DiagnosticCodeV2::RunnerRefused),
            ProofExitV2::NoData => Some(DiagnosticCodeV2::RunnerNoData),
            ProofExitV2::Stale => Some(DiagnosticCodeV2::RunnerStale),
            ProofExitV2::EnvironmentInvalid => Some(DiagnosticCodeV2::RunnerEnvironmentInvalid),
            ProofExitV2::Blocked => Some(DiagnosticCodeV2::RunnerBlocked),
            ProofExitV2::Unsupported => Some(DiagnosticCodeV2::RunnerUnsupported),
            ProofExitV2::NotRun => Some(DiagnosticCodeV2::RunnerNotRun),
            ProofExitV2::Cancelled => Some(DiagnosticCodeV2::RunnerCancelled),
            ProofExitV2::TimedOut => Some(DiagnosticCodeV2::RunnerTimedOut),
            ProofExitV2::Usage => Some(DiagnosticCodeV2::RunnerUsage),
            ProofExitV2::InternalError => Some(DiagnosticCodeV2::RunnerInternalError),
        }
    }

    fn oracle_drain(
        role: StateBearingRecordRoleV2,
        state: ProofExitV2,
    ) -> Option<PresentedDrainRootKindV2> {
        if role == StateBearingRecordRoleV2::PreRunDiagnostic {
            return None;
        }
        match state {
            ProofExitV2::Cancelled => Some(PresentedDrainRootKindV2::CancelledStopRoot),
            ProofExitV2::TimedOut => Some(PresentedDrainRootKindV2::TimedOutStopRoot),
            ProofExitV2::InternalError => Some(PresentedDrainRootKindV2::DrainedInternalErrorRoot),
            _ => None,
        }
    }

    fn oracle_record_kind(role: StateBearingRecordRoleV2) -> Option<LifecycleRecordKindV2> {
        match role {
            StateBearingRecordRoleV2::PreRunDiagnostic => None,
            StateBearingRecordRoleV2::ExecutedCaseTerminal
            | StateBearingRecordRoleV2::SuppressedCaseTerminal => {
                Some(LifecycleRecordKindV2::CaseTerminal)
            }
            StateBearingRecordRoleV2::RunTerminal => Some(LifecycleRecordKindV2::RunTerminal),
        }
    }

    fn oracle(
        role: StateBearingRecordRoleV2,
        state: ProofExitV2,
        reason: Option<RefusedReasonV2>,
        diagnostic: Option<DiagnosticCodeV2>,
        drain: Option<PresentedDrainRootKindV2>,
    ) -> Result<(), StateValidationErrorV2> {
        if !oracle_role_allows(role, state) {
            return Err(StateValidationErrorV2::StateNotAllowedForRole { role, state });
        }

        match (state, reason) {
            (ProofExitV2::Refused, None) => {
                return Err(StateValidationErrorV2::MissingRefusedReason);
            }
            (ProofExitV2::Refused, Some(_)) | (_, None) => {}
            (_, Some(observed)) => {
                return Err(StateValidationErrorV2::UnexpectedRefusedReason { state, observed });
            }
        }

        match (oracle_diagnostic(state), diagnostic) {
            (None, Some(observed)) => {
                return Err(StateValidationErrorV2::UnexpectedDiagnostic { observed });
            }
            (Some(expected), None) => {
                return Err(StateValidationErrorV2::MissingDiagnostic { expected });
            }
            (Some(expected), Some(observed)) if expected != observed => {
                return Err(StateValidationErrorV2::WrongDiagnostic { expected, observed });
            }
            (None, None) | (Some(_), Some(_)) => {}
        }

        match (oracle_drain(role, state), drain) {
            (None, Some(observed)) => {
                return Err(StateValidationErrorV2::UnexpectedDrainBasis { state, observed });
            }
            (Some(expected), None) => {
                return Err(StateValidationErrorV2::MissingDrainBasis { expected });
            }
            (Some(expected), Some(observed)) if expected != observed => {
                return Err(StateValidationErrorV2::WrongDrainBasis { expected, observed });
            }
            (None, None) | (Some(_), Some(_)) => {}
        }

        Ok(())
    }

    #[test]
    fn exact_diagnostic_mapping_literal_oracle() {
        let expected = [
            (ProofExitV2::Pass, None),
            (
                ProofExitV2::Failed,
                Some(DiagnosticCodeV2::CaseConformanceMismatch),
            ),
            (ProofExitV2::Refused, Some(DiagnosticCodeV2::RunnerRefused)),
            (ProofExitV2::NoData, Some(DiagnosticCodeV2::RunnerNoData)),
            (ProofExitV2::Stale, Some(DiagnosticCodeV2::RunnerStale)),
            (
                ProofExitV2::EnvironmentInvalid,
                Some(DiagnosticCodeV2::RunnerEnvironmentInvalid),
            ),
            (ProofExitV2::Blocked, Some(DiagnosticCodeV2::RunnerBlocked)),
            (
                ProofExitV2::Unsupported,
                Some(DiagnosticCodeV2::RunnerUnsupported),
            ),
            (ProofExitV2::NotRun, Some(DiagnosticCodeV2::RunnerNotRun)),
            (
                ProofExitV2::Cancelled,
                Some(DiagnosticCodeV2::RunnerCancelled),
            ),
            (
                ProofExitV2::TimedOut,
                Some(DiagnosticCodeV2::RunnerTimedOut),
            ),
            (ProofExitV2::Usage, Some(DiagnosticCodeV2::RunnerUsage)),
            (
                ProofExitV2::InternalError,
                Some(DiagnosticCodeV2::RunnerInternalError),
            ),
        ];
        for (state, diagnostic) in expected {
            assert_eq!(diagnostic_for_state(state), diagnostic, "{state:?}");
        }
        assert_eq!(
            diagnostic_for_state(ProofExitV2::Failed)
                .expect("failed diagnostic")
                .name(),
            "case.conformance_mismatch"
        );
        assert_ne!(
            DiagnosticCodeV2::CaseConformanceMismatch.name(),
            "comparison.mismatch"
        );
    }

    #[test]
    fn role_state_sets_and_lifecycle_basis_are_exact() {
        let pre_run = [
            ProofExitV2::Refused,
            ProofExitV2::NoData,
            ProofExitV2::Stale,
            ProofExitV2::EnvironmentInvalid,
            ProofExitV2::Blocked,
            ProofExitV2::Unsupported,
            ProofExitV2::Cancelled,
            ProofExitV2::TimedOut,
            ProofExitV2::Usage,
            ProofExitV2::InternalError,
        ];
        let active_terminal = [
            ProofExitV2::Pass,
            ProofExitV2::Failed,
            ProofExitV2::Refused,
            ProofExitV2::NoData,
            ProofExitV2::Stale,
            ProofExitV2::EnvironmentInvalid,
            ProofExitV2::Blocked,
            ProofExitV2::Unsupported,
            ProofExitV2::Cancelled,
            ProofExitV2::TimedOut,
            ProofExitV2::InternalError,
        ];
        for state in STATES {
            assert_eq!(
                role_allows_state(StateBearingRecordRoleV2::PreRunDiagnostic, state),
                pre_run.contains(&state)
            );
            assert_eq!(
                role_allows_state(StateBearingRecordRoleV2::ExecutedCaseTerminal, state),
                active_terminal.contains(&state)
            );
            assert_eq!(
                role_allows_state(StateBearingRecordRoleV2::RunTerminal, state),
                active_terminal.contains(&state)
            );
            assert_eq!(
                role_allows_state(StateBearingRecordRoleV2::SuppressedCaseTerminal, state),
                state == ProofExitV2::NotRun
            );
        }

        assert_eq!(
            lifecycle_record_kind_for_role(StateBearingRecordRoleV2::PreRunDiagnostic),
            None
        );
        assert_eq!(
            lifecycle_record_kind_for_role(StateBearingRecordRoleV2::ExecutedCaseTerminal),
            Some(LifecycleRecordKindV2::CaseTerminal)
        );
        assert_eq!(
            lifecycle_record_kind_for_role(StateBearingRecordRoleV2::SuppressedCaseTerminal),
            Some(LifecycleRecordKindV2::CaseTerminal)
        );
        assert_eq!(
            lifecycle_record_kind_for_role(StateBearingRecordRoleV2::RunTerminal),
            Some(LifecycleRecordKindV2::RunTerminal)
        );

        for state in [
            ProofExitV2::Cancelled,
            ProofExitV2::TimedOut,
            ProofExitV2::InternalError,
        ] {
            assert_eq!(
                drain_basis_for(StateBearingRecordRoleV2::PreRunDiagnostic, state),
                None,
                "pre-run stop diagnostics have no lifecycle drain basis"
            );
        }
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn exhaustive_cartesian_matrix_and_error_precedence() {
        let mut visited = 0usize;
        let mut accepted = 0usize;

        for role in ROLES {
            for state in STATES {
                for reason in REASONS {
                    for diagnostic in DIAGNOSTICS {
                        for drain in DRAIN_BASES {
                            visited += 1;
                            let expected = oracle(role, state, reason, diagnostic, drain);
                            let actual =
                                validate_state_cell_v2(role, state, reason, diagnostic, drain);
                            match (expected, actual) {
                                (Err(expected), Err(actual)) => {
                                    assert_eq!(
                                        actual, expected,
                                        "precedence drift for role={role:?} state={state:?} \
                                         reason={reason:?} diagnostic={diagnostic:?} \
                                         drain={drain:?}"
                                    );
                                }
                                (Ok(()), Ok(cell)) => {
                                    accepted += 1;
                                    assert_eq!(cell.role(), role);
                                    assert_eq!(cell.state(), state);
                                    assert_eq!(cell.refused_reason(), reason);
                                    assert_eq!(cell.diagnostic(), diagnostic);
                                    assert_eq!(cell.drain_basis(), drain);
                                    assert_eq!(
                                        cell.lifecycle_record_kind(),
                                        oracle_record_kind(role)
                                    );
                                }
                                (expected, actual) => panic!(
                                    "matrix drift for role={role:?} state={state:?} \
                                     reason={reason:?} diagnostic={diagnostic:?} \
                                     drain={drain:?}: expected={expected:?} actual={actual:?}"
                                ),
                            }
                        }
                    }
                }
            }
        }

        assert_eq!(visited, 13 * 12 * 4 * 13 * 4);
        assert_eq!(accepted, 63, "exactly 63 reason-expanded cells are valid");
    }

    #[test]
    fn active_stop_states_require_the_exact_nominal_basis() {
        let cases = [
            (
                ProofExitV2::Cancelled,
                DiagnosticCodeV2::RunnerCancelled,
                PresentedDrainRootKindV2::CancelledStopRoot,
                PresentedDrainRootKindV2::TimedOutStopRoot,
            ),
            (
                ProofExitV2::TimedOut,
                DiagnosticCodeV2::RunnerTimedOut,
                PresentedDrainRootKindV2::TimedOutStopRoot,
                PresentedDrainRootKindV2::DrainedInternalErrorRoot,
            ),
            (
                ProofExitV2::InternalError,
                DiagnosticCodeV2::RunnerInternalError,
                PresentedDrainRootKindV2::DrainedInternalErrorRoot,
                PresentedDrainRootKindV2::CancelledStopRoot,
            ),
        ];

        for (state, diagnostic, expected, wrong) in cases {
            assert!(matches!(
                validate_state_cell_v2(
                    StateBearingRecordRoleV2::RunTerminal,
                    state,
                    None,
                    Some(diagnostic),
                    None,
                ),
                Err(StateValidationErrorV2::MissingDrainBasis {
                    expected: observed,
                }) if observed == expected
            ));
            assert!(matches!(
                validate_state_cell_v2(
                    StateBearingRecordRoleV2::RunTerminal,
                    state,
                    None,
                    Some(diagnostic),
                    Some(wrong),
                ),
                Err(StateValidationErrorV2::WrongDrainBasis {
                    expected: observed_expected,
                    observed,
                }) if observed_expected == expected && observed == wrong
            ));
            let cell = validate_state_cell_v2(
                StateBearingRecordRoleV2::RunTerminal,
                state,
                None,
                Some(diagnostic),
                Some(expected),
            )
            .expect("exact nominal drain basis must validate");
            assert_eq!(cell.drain_basis(), Some(expected));
        }
    }
}
