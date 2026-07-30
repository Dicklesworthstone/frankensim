//! Bounded, rootless handoff from one local Stage-A evaluator.
//!
//! The handoff is intentionally noncanonical and non-authoritative. It carries
//! only stable cell identity, raw checked outcomes, safe typed numeric
//! observations, and bounded structured diagnostics. Expected outcomes,
//! oracle roots, attempts, AC57 data, runtime Five Explicits, receipts,
//! telemetry, retained artifacts, and authority are structurally absent.

use crate::catalog::{DiagnosticCodeV2, LogicalUnitV2, RepairActionKindV2, RetryabilityV2};
use crate::construction::{ConstructionErrorKindV2, ConstructionErrorV2};
use crate::limits::{RunnerLimitUnitV2, RunnerLimitValueV2};
use crate::value::{NumericValueV2, StableTokenV2, UnitV2};

/// Maximum cells in one complete foundational work-package handoff.
pub const RUNNER_V2_LOCAL_HANDOFF_MAX_CELLS_V1: usize = 2_048;
/// Maximum safe numeric observations attached to one raw cell.
pub const RUNNER_V2_LOCAL_HANDOFF_MAX_NUMERIC_OBSERVATIONS_V1: usize = 8;
/// Maximum prerequisites retained by one raw diagnostic.
pub const RUNNER_V2_LOCAL_HANDOFF_MAX_PREREQUISITES_V1: usize = 8;
/// Maximum non-executable repairs retained by one raw diagnostic.
pub const RUNNER_V2_LOCAL_HANDOFF_MAX_REPAIRS_V1: usize = 4;

/// Raw outcome produced by a real local checked evaluator.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u16)]
pub enum RunnerV2RawOutcomeKindV1 {
    /// The checked operation accepted the presented domain value.
    Accepted = 1,
    /// The checked operation refused the presented domain value.
    Refused = 2,
    /// The checked operation reached one explicit modeled failure.
    Failed = 3,
    /// The closed implementation does not support this registered cell.
    Unsupported = 4,
    /// The registered facet has no runtime operation in this pure evaluator.
    Inapplicable = 5,
}

impl RunnerV2RawOutcomeKindV1 {
    /// Exact stable code.
    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }
}

/// Bounded, closed reason attached to a raw local outcome.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u16)]
pub enum RunnerV2RawReasonV1 {
    /// The checked operation accepted its exact input.
    ExactCheckedValue = 1,
    /// A value fell below its structural minimum.
    BelowStructuralMinimum = 2,
    /// A value exceeded its registered profile ceiling.
    AboveProfileCeiling = 3,
    /// A fixed representation property was changed.
    FixedRepresentationChanged = 4,
    /// A heterogeneous value used the wrong primitive width.
    WrongPrimitiveWidth = 5,
    /// Checked representational arithmetic overflowed.
    CheckedRepresentationalOverflow = 6,
    /// Individually valid fields violated a joint equation.
    JointFeasibilityViolation = 7,
    /// A literal, code, tag, role, or nominal domain was unknown.
    UnknownClosedValue = 8,
    /// A value was malformed or noncanonical.
    MalformedOrNoncanonical = 9,
    /// A required value was absent.
    RequiredValueAbsent = 10,
    /// A value was present where typed absence was required.
    UnexpectedValuePresent = 11,
    /// A set or ordered declaration had the wrong exact membership.
    ExactMembershipMismatch = 12,
    /// A pure compile-time or declaration facet has no runtime operation.
    PureDeclarationFacet = 13,
    /// Cancellation is inapplicable to one bounded, allocation-checked cell.
    CancellationInapplicable = 14,
    /// Sharding is inapplicable to this deterministic local evaluator.
    ShardInapplicable = 15,
    /// Resume is inapplicable; a new invocation evaluates the full package.
    ResumeInapplicable = 16,
    /// A checked source declaration or source fragment was incompatible.
    SourceDeclarationMismatch = 17,
    /// An implementation invariant failed after checked construction.
    InternalInvariantFailure = 18,
    /// The closed evaluator has no implementation for a registered value.
    UnsupportedClosedValue = 19,
}

impl RunnerV2RawReasonV1 {
    /// Exact stable code.
    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }
}

/// Safe typed numeric value retained in a raw handoff.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RunnerV2SafeNumericValueV1 {
    /// One value from the closed Runner numeric union.
    Numeric(NumericValueV2),
    /// One heterogeneous Runner limit value.
    Limit(RunnerLimitValueV2),
    /// One exact bounded count.
    Count(u64),
}

/// Required unit carried by every safe numeric observation.
///
/// Limit units remain separate from the general logical-unit catalog so a
/// limit observation cannot silently lose the field-specific semantic unit
/// used by limit admission.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RunnerV2SafeNumericUnitV1 {
    /// One canonical physical unit with exact scale and dimensions.
    Physical(UnitV2),
    /// One unit from the closed Runner logical-unit catalog.
    Logical(LogicalUnitV2),
    /// One exact unit from the heterogeneous Runner limit catalog.
    Limit(RunnerLimitUnitV2),
}

/// One named, unit-bearing safe numeric observation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunnerV2SafeNumericObservationV1 {
    name: StableTokenV2,
    value: RunnerV2SafeNumericValueV1,
    unit: RunnerV2SafeNumericUnitV1,
}

impl RunnerV2SafeNumericObservationV1 {
    pub(crate) fn limit(
        name: StableTokenV2,
        value: RunnerLimitValueV2,
        unit: RunnerLimitUnitV2,
    ) -> Self {
        Self {
            name,
            value: RunnerV2SafeNumericValueV1::Limit(value),
            unit: RunnerV2SafeNumericUnitV1::Limit(unit),
        }
    }

    pub(crate) fn numeric(name: StableTokenV2, value: NumericValueV2, unit: LogicalUnitV2) -> Self {
        Self {
            name,
            value: RunnerV2SafeNumericValueV1::Numeric(value),
            unit: RunnerV2SafeNumericUnitV1::Logical(unit),
        }
    }

    pub(crate) fn physical(name: StableTokenV2, value: NumericValueV2, unit: UnitV2) -> Self {
        Self {
            name,
            value: RunnerV2SafeNumericValueV1::Numeric(value),
            unit: RunnerV2SafeNumericUnitV1::Physical(unit),
        }
    }

    pub(crate) fn count(name: StableTokenV2, value: u64) -> Self {
        Self {
            name,
            value: RunnerV2SafeNumericValueV1::Count(value),
            unit: RunnerV2SafeNumericUnitV1::Logical(LogicalUnitV2::Count),
        }
    }

    /// Stable observation name.
    #[must_use]
    pub const fn name(&self) -> &StableTokenV2 {
        &self.name
    }

    /// Safe typed value.
    #[must_use]
    pub const fn value(&self) -> &RunnerV2SafeNumericValueV1 {
        &self.value
    }

    /// Exact required unit.
    #[must_use]
    pub const fn unit(&self) -> RunnerV2SafeNumericUnitV1 {
        self.unit
    }
}

/// One ranked, non-executable repair descriptor.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunnerV2RawRepairV1 {
    rank: u8,
    kind: RepairActionKindV2,
    target: StableTokenV2,
}

impl RunnerV2RawRepairV1 {
    pub(crate) const fn new(rank: u8, kind: RepairActionKindV2, target: StableTokenV2) -> Self {
        Self { rank, kind, target }
    }

    /// One-based contiguous rank.
    #[must_use]
    pub const fn rank(&self) -> u8 {
        self.rank
    }

    /// Closed non-executable repair class.
    #[must_use]
    pub const fn kind(&self) -> RepairActionKindV2 {
        self.kind
    }

    /// Stable semantic repair target.
    #[must_use]
    pub const fn target(&self) -> &StableTokenV2 {
        &self.target
    }
}

/// One bounded structured diagnostic produced by a raw evaluator cell.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunnerV2RawDiagnosticV1 {
    code: DiagnosticCodeV2,
    owner: StableTokenV2,
    retryability: RetryabilityV2,
    prerequisites: Box<[StableTokenV2]>,
    repairs: Box<[RunnerV2RawRepairV1]>,
}

impl RunnerV2RawDiagnosticV1 {
    pub(crate) fn new(
        code: DiagnosticCodeV2,
        owner: StableTokenV2,
        retryability: RetryabilityV2,
        prerequisites: Vec<StableTokenV2>,
        repairs: Vec<RunnerV2RawRepairV1>,
    ) -> Result<Self, ConstructionErrorV2> {
        if prerequisites.len() > RUNNER_V2_LOCAL_HANDOFF_MAX_PREREQUISITES_V1 {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                "runner_v2.handoff.diagnostic.prerequisites",
                "at most eight stable prerequisites",
                prerequisites.len(),
            ));
        }
        if repairs.len() > RUNNER_V2_LOCAL_HANDOFF_MAX_REPAIRS_V1 {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                "runner_v2.handoff.diagnostic.repairs",
                "at most four ranked non-executable repairs",
                repairs.len(),
            ));
        }
        validate_strict_token_order("runner_v2.handoff.diagnostic.prerequisites", &prerequisites)?;
        for (index, repair) in repairs.iter().enumerate() {
            let expected_rank = u8::try_from(index + 1).map_err(|_| {
                ConstructionErrorV2::new(
                    ConstructionErrorKindV2::ArithmeticOverflow,
                    "runner_v2.handoff.diagnostic.repair_rank",
                    "a one-based rank representable as u8",
                    index,
                )
            })?;
            if repair.rank != expected_rank {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::OutOfRange,
                    "runner_v2.handoff.diagnostic.repair_rank",
                    "contiguous one-based ranks in presented order",
                    repair.rank,
                ));
            }
        }
        for (index, repair) in repairs.iter().enumerate() {
            if repairs[..index]
                .iter()
                .any(|prior| prior.kind == repair.kind && prior.target == repair.target)
            {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Duplicate,
                    "runner_v2.handoff.diagnostic.repairs",
                    "distinct ranked repair kind and target alternatives",
                    index,
                ));
            }
        }
        Ok(Self {
            code,
            owner,
            retryability,
            prerequisites: prerequisites.into_boxed_slice(),
            repairs: repairs.into_boxed_slice(),
        })
    }

    /// Closed diagnostic code.
    #[must_use]
    pub const fn code(&self) -> DiagnosticCodeV2 {
        self.code
    }

    /// Stable source owner.
    #[must_use]
    pub const fn owner(&self) -> &StableTokenV2 {
        &self.owner
    }

    /// Closed retry classification.
    #[must_use]
    pub const fn retryability(&self) -> RetryabilityV2 {
        self.retryability
    }

    /// Exact ordered prerequisites.
    #[must_use]
    pub fn prerequisites(&self) -> &[StableTokenV2] {
        &self.prerequisites
    }

    /// Exact ranked non-executable repairs.
    #[must_use]
    pub fn repairs(&self) -> &[RunnerV2RawRepairV1] {
        &self.repairs
    }
}

/// One raw local cell observation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunnerV2RawCellObservationV1 {
    cell_id: StableTokenV2,
    outcome: RunnerV2RawOutcomeKindV1,
    reason: RunnerV2RawReasonV1,
    numeric: Box<[RunnerV2SafeNumericObservationV1]>,
    diagnostic: Option<RunnerV2RawDiagnosticV1>,
}

impl RunnerV2RawCellObservationV1 {
    pub(crate) fn new(
        cell_id: StableTokenV2,
        outcome: RunnerV2RawOutcomeKindV1,
        reason: RunnerV2RawReasonV1,
        numeric: Vec<RunnerV2SafeNumericObservationV1>,
        diagnostic: Option<RunnerV2RawDiagnosticV1>,
    ) -> Result<Self, ConstructionErrorV2> {
        if numeric.len() > RUNNER_V2_LOCAL_HANDOFF_MAX_NUMERIC_OBSERVATIONS_V1 {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                "runner_v2.handoff.cell.numeric",
                "at most eight safe numeric observations",
                numeric.len(),
            ));
        }
        validate_observation_order(&numeric)?;
        if matches!(outcome, RunnerV2RawOutcomeKindV1::Accepted) && diagnostic.is_some() {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Unexpected,
                "runner_v2.handoff.cell.diagnostic",
                "no diagnostic on an accepted raw cell",
                true,
            ));
        }
        if !matches!(outcome, RunnerV2RawOutcomeKindV1::Accepted) && diagnostic.is_none() {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Missing,
                "runner_v2.handoff.cell.diagnostic",
                "one bounded diagnostic for every non-accepted raw cell",
                false,
            ));
        }
        let inapplicable_reason = matches!(
            reason,
            RunnerV2RawReasonV1::PureDeclarationFacet
                | RunnerV2RawReasonV1::CancellationInapplicable
                | RunnerV2RawReasonV1::ShardInapplicable
                | RunnerV2RawReasonV1::ResumeInapplicable
        );
        let reason_compatible = match outcome {
            RunnerV2RawOutcomeKindV1::Accepted => {
                matches!(reason, RunnerV2RawReasonV1::ExactCheckedValue)
            }
            RunnerV2RawOutcomeKindV1::Inapplicable => inapplicable_reason,
            RunnerV2RawOutcomeKindV1::Refused => matches!(
                reason,
                RunnerV2RawReasonV1::BelowStructuralMinimum
                    | RunnerV2RawReasonV1::AboveProfileCeiling
                    | RunnerV2RawReasonV1::FixedRepresentationChanged
                    | RunnerV2RawReasonV1::WrongPrimitiveWidth
                    | RunnerV2RawReasonV1::CheckedRepresentationalOverflow
                    | RunnerV2RawReasonV1::JointFeasibilityViolation
                    | RunnerV2RawReasonV1::UnknownClosedValue
                    | RunnerV2RawReasonV1::MalformedOrNoncanonical
                    | RunnerV2RawReasonV1::RequiredValueAbsent
                    | RunnerV2RawReasonV1::UnexpectedValuePresent
                    | RunnerV2RawReasonV1::ExactMembershipMismatch
                    | RunnerV2RawReasonV1::SourceDeclarationMismatch
            ),
            RunnerV2RawOutcomeKindV1::Failed => {
                matches!(reason, RunnerV2RawReasonV1::InternalInvariantFailure)
            }
            RunnerV2RawOutcomeKindV1::Unsupported => {
                matches!(reason, RunnerV2RawReasonV1::UnsupportedClosedValue)
            }
        };
        if !reason_compatible {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Incompatible,
                "runner_v2.handoff.cell.outcome_reason",
                "the closed outcome and reason compatibility matrix",
                reason.code(),
            ));
        }
        if let Some(diagnostic) = &diagnostic {
            let (expected_code, expected_retryability) = match outcome {
                RunnerV2RawOutcomeKindV1::Accepted => {
                    return Err(ConstructionErrorV2::new(
                        ConstructionErrorKindV2::Unexpected,
                        "runner_v2.handoff.cell.diagnostic",
                        "no diagnostic for an accepted raw cell",
                        diagnostic.code.code(),
                    ));
                }
                RunnerV2RawOutcomeKindV1::Refused => (
                    DiagnosticCodeV2::RunnerRefused,
                    RetryabilityV2::AfterInputChange,
                ),
                RunnerV2RawOutcomeKindV1::Failed => {
                    (DiagnosticCodeV2::RunnerInternalError, RetryabilityV2::Never)
                }
                RunnerV2RawOutcomeKindV1::Unsupported => (
                    DiagnosticCodeV2::RunnerUnsupported,
                    RetryabilityV2::AfterPrerequisiteChange,
                ),
                RunnerV2RawOutcomeKindV1::Inapplicable => {
                    (DiagnosticCodeV2::RunnerNotRun, RetryabilityV2::Never)
                }
            };
            if diagnostic.code != expected_code {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Incompatible,
                    "runner_v2.handoff.cell.diagnostic_code",
                    "the exact diagnostic code for the raw outcome",
                    diagnostic.code.code(),
                ));
            }
            if diagnostic.retryability != expected_retryability {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Incompatible,
                    "runner_v2.handoff.cell.retryability",
                    "the exact retryability for the raw outcome",
                    diagnostic.retryability.code(),
                ));
            }
            if matches!(
                outcome,
                RunnerV2RawOutcomeKindV1::Inapplicable | RunnerV2RawOutcomeKindV1::Unsupported
            ) && diagnostic.prerequisites.is_empty()
            {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Missing,
                    "runner_v2.handoff.cell.applicability_prerequisite",
                    "one or more stable applicability prerequisites",
                    0_usize,
                ));
            }
            if matches!(outcome, RunnerV2RawOutcomeKindV1::Inapplicable)
                && !diagnostic.repairs.is_empty()
            {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Unexpected,
                    "runner_v2.handoff.cell.inapplicable_repairs",
                    "no repair for a permanently inapplicable local facet",
                    diagnostic.repairs.len(),
                ));
            }
        }
        Ok(Self {
            cell_id,
            outcome,
            reason,
            numeric: numeric.into_boxed_slice(),
            diagnostic,
        })
    }

    /// Stable source-owned cell identity.
    #[must_use]
    pub const fn cell_id(&self) -> &StableTokenV2 {
        &self.cell_id
    }

    /// Raw checked outcome.
    #[must_use]
    pub const fn outcome(&self) -> RunnerV2RawOutcomeKindV1 {
        self.outcome
    }

    /// Closed raw reason.
    #[must_use]
    pub const fn reason(&self) -> RunnerV2RawReasonV1 {
        self.reason
    }

    /// Safe typed numeric observations.
    #[must_use]
    pub fn numeric(&self) -> &[RunnerV2SafeNumericObservationV1] {
        &self.numeric
    }

    /// Bounded diagnostic for a non-accepted outcome.
    #[must_use]
    pub const fn diagnostic(&self) -> Option<&RunnerV2RawDiagnosticV1> {
        self.diagnostic.as_ref()
    }
}

/// Complete rootless raw report returned by one local evaluator.
///
/// The type exposes only inspection of the complete source-ordered raw rows:
///
/// ```
/// use fs_evidence_runner::runner_v2::RunnerV2LocalWorkPackageHandoffV1;
///
/// fn inspect(report: &RunnerV2LocalWorkPackageHandoffV1) {
///     let _ = report.package_id();
///     let _ = report.cells();
/// }
/// ```
///
/// Its fields are private:
///
/// ```compile_fail,E0451
/// use fs_evidence_runner::runner_v2::RunnerV2LocalWorkPackageHandoffV1;
///
/// let _ = RunnerV2LocalWorkPackageHandoffV1 {
///     package_id: panic!(),
///     cells: Vec::new().into_boxed_slice(),
/// };
/// ```
///
/// It has no canonical or nominal root:
///
/// ```compile_fail,E0599
/// use fs_evidence_runner::runner_v2::RunnerV2LocalWorkPackageHandoffV1;
///
/// fn forbidden(report: &RunnerV2LocalWorkPackageHandoffV1) {
///     let _ = report.root();
/// }
/// ```
///
/// It has no attempt or AC57 accessor:
///
/// ```compile_fail,E0599
/// use fs_evidence_runner::runner_v2::RunnerV2LocalWorkPackageHandoffV1;
///
/// fn forbidden(report: &RunnerV2LocalWorkPackageHandoffV1) {
///     let _ = report.attempt_identity();
/// }
/// ```
///
/// It cannot substitute for canonical content:
///
/// ```compile_fail,E0308
/// use fs_blake3::ContentHash;
/// use fs_evidence_runner::runner_v2::RunnerV2LocalWorkPackageHandoffV1;
///
/// fn canonical(_: ContentHash) {}
/// fn forbidden(report: RunnerV2LocalWorkPackageHandoffV1) {
///     canonical(report);
/// }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunnerV2LocalWorkPackageHandoffV1 {
    package_id: StableTokenV2,
    cells: Box<[RunnerV2RawCellObservationV1]>,
}

impl RunnerV2LocalWorkPackageHandoffV1 {
    pub(crate) fn new(
        package_id: StableTokenV2,
        declared_cell_ids: &[StableTokenV2],
        cells: Vec<RunnerV2RawCellObservationV1>,
    ) -> Result<Self, ConstructionErrorV2> {
        if declared_cell_ids.is_empty() {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Missing,
                "runner_v2.handoff.declared_cells",
                "one complete nonempty source-ordered declaration",
                0_usize,
            ));
        }
        if declared_cell_ids.len() > RUNNER_V2_LOCAL_HANDOFF_MAX_CELLS_V1 {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                "runner_v2.handoff.declared_cells",
                "at most 2048 declared cells",
                declared_cell_ids.len(),
            ));
        }
        validate_unique_tokens("runner_v2.handoff.declared_cells", declared_cell_ids)?;
        if cells.len() < declared_cell_ids.len() {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Missing,
                "runner_v2.handoff.cells",
                "exactly one raw cell for every declared cell",
                cells.len(),
            ));
        }
        if cells.len() > declared_cell_ids.len() {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Unexpected,
                "runner_v2.handoff.cells",
                "no raw cell beyond the complete declaration",
                cells.len(),
            ));
        }
        for (index, (declared, actual)) in declared_cell_ids.iter().zip(&cells).enumerate() {
            if declared != &actual.cell_id {
                return Err(ConstructionErrorV2::new(
                    if declared_cell_ids
                        .iter()
                        .any(|candidate| candidate == &actual.cell_id)
                    {
                        ConstructionErrorKindV2::OutOfOrder
                    } else {
                        ConstructionErrorKindV2::Incompatible
                    },
                    "runner_v2.handoff.cells",
                    "the complete declaration-side source order",
                    index,
                ));
            }
        }
        Ok(Self {
            package_id,
            cells: cells.into_boxed_slice(),
        })
    }

    /// Stable work-package identity.
    #[must_use]
    pub const fn package_id(&self) -> &StableTokenV2 {
        &self.package_id
    }

    /// Complete source-ordered raw cell observations.
    #[must_use]
    pub fn cells(&self) -> &[RunnerV2RawCellObservationV1] {
        &self.cells
    }

    /// Validate exact source order against one declaration-side cell sequence.
    ///
    /// This is a structural comparison only. It does not compare expected
    /// results, create a canonical root, or mint execution evidence.
    pub fn validate_exact_cell_order(
        &self,
        expected: &[StableTokenV2],
    ) -> Result<(), ConstructionErrorV2> {
        if expected.len() < self.cells.len() {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Unexpected,
                "runner_v2.handoff.expected_cells",
                "exactly one expected ID for every raw cell",
                expected.len(),
            ));
        }
        if expected.len() > self.cells.len() {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Missing,
                "runner_v2.handoff.expected_cells",
                "exactly one expected ID for every raw cell",
                expected.len(),
            ));
        }
        for (index, (actual, expected_id)) in self.cells.iter().zip(expected).enumerate() {
            if actual.cell_id != *expected_id {
                return Err(ConstructionErrorV2::new(
                    if expected
                        .iter()
                        .any(|candidate| candidate == &actual.cell_id)
                    {
                        ConstructionErrorKindV2::OutOfOrder
                    } else {
                        ConstructionErrorKindV2::Incompatible
                    },
                    "runner_v2.handoff.expected_cells",
                    "the exact declaration-side stable cell order",
                    index,
                ));
            }
        }
        Ok(())
    }
}

fn validate_observation_order(
    observations: &[RunnerV2SafeNumericObservationV1],
) -> Result<(), ConstructionErrorV2> {
    for pair in observations.windows(2) {
        if pair[0].name.as_str() >= pair[1].name.as_str() {
            return Err(ConstructionErrorV2::new(
                if pair[0].name == pair[1].name {
                    ConstructionErrorKindV2::Duplicate
                } else {
                    ConstructionErrorKindV2::OutOfOrder
                },
                "runner_v2.handoff.cell.numeric",
                "strict stable-name order with no duplicates",
                observations.len(),
            ));
        }
    }
    Ok(())
}

fn validate_strict_token_order(
    field: &'static str,
    values: &[StableTokenV2],
) -> Result<(), ConstructionErrorV2> {
    for pair in values.windows(2) {
        if pair[0].as_str() >= pair[1].as_str() {
            return Err(ConstructionErrorV2::new(
                if pair[0] == pair[1] {
                    ConstructionErrorKindV2::Duplicate
                } else {
                    ConstructionErrorKindV2::OutOfOrder
                },
                field,
                "strict stable-token order with no duplicates",
                values.len(),
            ));
        }
    }
    Ok(())
}

fn validate_unique_tokens(
    field: &'static str,
    values: &[StableTokenV2],
) -> Result<(), ConstructionErrorV2> {
    for (index, value) in values.iter().enumerate() {
        if values[..index].contains(value) {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Duplicate,
                field,
                "one occurrence of every stable token",
                index,
            ));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::catalog::{DiagnosticCodeV2, RepairActionKindV2, RetryabilityV2};

    const OUTCOME_CODES: [(RunnerV2RawOutcomeKindV1, u16); 5] = [
        (RunnerV2RawOutcomeKindV1::Accepted, 1),
        (RunnerV2RawOutcomeKindV1::Refused, 2),
        (RunnerV2RawOutcomeKindV1::Failed, 3),
        (RunnerV2RawOutcomeKindV1::Unsupported, 4),
        (RunnerV2RawOutcomeKindV1::Inapplicable, 5),
    ];

    const REASON_CODES: [(RunnerV2RawReasonV1, u16); 19] = [
        (RunnerV2RawReasonV1::ExactCheckedValue, 1),
        (RunnerV2RawReasonV1::BelowStructuralMinimum, 2),
        (RunnerV2RawReasonV1::AboveProfileCeiling, 3),
        (RunnerV2RawReasonV1::FixedRepresentationChanged, 4),
        (RunnerV2RawReasonV1::WrongPrimitiveWidth, 5),
        (RunnerV2RawReasonV1::CheckedRepresentationalOverflow, 6),
        (RunnerV2RawReasonV1::JointFeasibilityViolation, 7),
        (RunnerV2RawReasonV1::UnknownClosedValue, 8),
        (RunnerV2RawReasonV1::MalformedOrNoncanonical, 9),
        (RunnerV2RawReasonV1::RequiredValueAbsent, 10),
        (RunnerV2RawReasonV1::UnexpectedValuePresent, 11),
        (RunnerV2RawReasonV1::ExactMembershipMismatch, 12),
        (RunnerV2RawReasonV1::PureDeclarationFacet, 13),
        (RunnerV2RawReasonV1::CancellationInapplicable, 14),
        (RunnerV2RawReasonV1::ShardInapplicable, 15),
        (RunnerV2RawReasonV1::ResumeInapplicable, 16),
        (RunnerV2RawReasonV1::SourceDeclarationMismatch, 17),
        (RunnerV2RawReasonV1::InternalInvariantFailure, 18),
        (RunnerV2RawReasonV1::UnsupportedClosedValue, 19),
    ];

    fn token(value: &str) -> StableTokenV2 {
        StableTokenV2::new(value).expect("test token")
    }

    fn ordered_tokens(prefix: &str, count: usize) -> Vec<StableTokenV2> {
        (0..count)
            .map(|index| token(&format!("{prefix}-{index:04}")))
            .collect()
    }

    fn ordered_numeric(count: usize) -> Vec<RunnerV2SafeNumericObservationV1> {
        (0..count)
            .map(|index| {
                RunnerV2SafeNumericObservationV1::count(
                    token(&format!("numeric-{index:02}")),
                    u64::try_from(index).expect("bounded test index"),
                )
            })
            .collect()
    }

    fn ranked_repairs(count: usize) -> Vec<RunnerV2RawRepairV1> {
        (0..count)
            .map(|index| {
                RunnerV2RawRepairV1::new(
                    u8::try_from(index + 1).expect("bounded test rank"),
                    RepairActionKindV2::ChangeArguments,
                    token(&format!("repair-target-{index:02}")),
                )
            })
            .collect()
    }

    fn accepted_with_numeric(
        id: StableTokenV2,
        numeric: Vec<RunnerV2SafeNumericObservationV1>,
    ) -> RunnerV2RawCellObservationV1 {
        RunnerV2RawCellObservationV1::new(
            id,
            RunnerV2RawOutcomeKindV1::Accepted,
            RunnerV2RawReasonV1::ExactCheckedValue,
            numeric,
            None,
        )
        .expect("accepted raw cell")
    }

    fn accepted_token(id: StableTokenV2) -> RunnerV2RawCellObservationV1 {
        accepted_with_numeric(id, vec![])
    }

    fn accepted(id: &str) -> RunnerV2RawCellObservationV1 {
        accepted_token(token(id))
    }

    fn diagnostic_with(
        code: DiagnosticCodeV2,
        retryability: RetryabilityV2,
        prerequisites: Vec<StableTokenV2>,
        repairs: Vec<RunnerV2RawRepairV1>,
    ) -> RunnerV2RawDiagnosticV1 {
        RunnerV2RawDiagnosticV1::new(
            code,
            token("fs-evidence-runner.runner-v2.handoff-tests"),
            retryability,
            prerequisites,
            repairs,
        )
        .expect("valid test diagnostic")
    }

    fn diagnostic_for_outcome(
        outcome: RunnerV2RawOutcomeKindV1,
    ) -> Option<RunnerV2RawDiagnosticV1> {
        let (code, retryability) = match outcome {
            RunnerV2RawOutcomeKindV1::Accepted => return None,
            RunnerV2RawOutcomeKindV1::Refused => (
                DiagnosticCodeV2::RunnerRefused,
                RetryabilityV2::AfterInputChange,
            ),
            RunnerV2RawOutcomeKindV1::Failed => {
                (DiagnosticCodeV2::RunnerInternalError, RetryabilityV2::Never)
            }
            RunnerV2RawOutcomeKindV1::Unsupported => (
                DiagnosticCodeV2::RunnerUnsupported,
                RetryabilityV2::AfterPrerequisiteChange,
            ),
            RunnerV2RawOutcomeKindV1::Inapplicable => {
                (DiagnosticCodeV2::RunnerNotRun, RetryabilityV2::Never)
            }
        };
        let prerequisites = if matches!(
            outcome,
            RunnerV2RawOutcomeKindV1::Unsupported | RunnerV2RawOutcomeKindV1::Inapplicable
        ) {
            vec![token("registered-applicability-prerequisite")]
        } else {
            vec![]
        };
        Some(diagnostic_with(code, retryability, prerequisites, vec![]))
    }

    fn reason_is_compatible(
        outcome: RunnerV2RawOutcomeKindV1,
        reason: RunnerV2RawReasonV1,
    ) -> bool {
        match outcome {
            RunnerV2RawOutcomeKindV1::Accepted => {
                matches!(reason, RunnerV2RawReasonV1::ExactCheckedValue)
            }
            RunnerV2RawOutcomeKindV1::Refused => matches!(
                reason,
                RunnerV2RawReasonV1::BelowStructuralMinimum
                    | RunnerV2RawReasonV1::AboveProfileCeiling
                    | RunnerV2RawReasonV1::FixedRepresentationChanged
                    | RunnerV2RawReasonV1::WrongPrimitiveWidth
                    | RunnerV2RawReasonV1::CheckedRepresentationalOverflow
                    | RunnerV2RawReasonV1::JointFeasibilityViolation
                    | RunnerV2RawReasonV1::UnknownClosedValue
                    | RunnerV2RawReasonV1::MalformedOrNoncanonical
                    | RunnerV2RawReasonV1::RequiredValueAbsent
                    | RunnerV2RawReasonV1::UnexpectedValuePresent
                    | RunnerV2RawReasonV1::ExactMembershipMismatch
                    | RunnerV2RawReasonV1::SourceDeclarationMismatch
            ),
            RunnerV2RawOutcomeKindV1::Failed => {
                matches!(reason, RunnerV2RawReasonV1::InternalInvariantFailure)
            }
            RunnerV2RawOutcomeKindV1::Unsupported => {
                matches!(reason, RunnerV2RawReasonV1::UnsupportedClosedValue)
            }
            RunnerV2RawOutcomeKindV1::Inapplicable => matches!(
                reason,
                RunnerV2RawReasonV1::PureDeclarationFacet
                    | RunnerV2RawReasonV1::CancellationInapplicable
                    | RunnerV2RawReasonV1::ShardInapplicable
                    | RunnerV2RawReasonV1::ResumeInapplicable
            ),
        }
    }

    fn semantic_reason(outcome: RunnerV2RawOutcomeKindV1) -> RunnerV2RawReasonV1 {
        match outcome {
            RunnerV2RawOutcomeKindV1::Accepted => RunnerV2RawReasonV1::ExactCheckedValue,
            RunnerV2RawOutcomeKindV1::Refused => RunnerV2RawReasonV1::MalformedOrNoncanonical,
            RunnerV2RawOutcomeKindV1::Failed => RunnerV2RawReasonV1::InternalInvariantFailure,
            RunnerV2RawOutcomeKindV1::Unsupported => RunnerV2RawReasonV1::UnsupportedClosedValue,
            RunnerV2RawOutcomeKindV1::Inapplicable => RunnerV2RawReasonV1::PureDeclarationFacet,
        }
    }

    fn assert_error(
        error: ConstructionErrorV2,
        kind: ConstructionErrorKindV2,
        field: &'static str,
        observed: &str,
    ) {
        assert_eq!(error.kind(), kind, "{error}");
        assert_eq!(error.field(), field, "{error}");
        assert_eq!(error.observed(), observed, "{error}");
    }

    fn refused(id: &str) -> RunnerV2RawCellObservationV1 {
        let diagnostic = diagnostic_with(
            DiagnosticCodeV2::RunnerRefused,
            RetryabilityV2::AfterInputChange,
            vec![token("present-canonical-input")],
            vec![RunnerV2RawRepairV1::new(
                1,
                RepairActionKindV2::ChangeArguments,
                token("canonical-input"),
            )],
        );
        RunnerV2RawCellObservationV1::new(
            token(id),
            RunnerV2RawOutcomeKindV1::Refused,
            RunnerV2RawReasonV1::MalformedOrNoncanonical,
            vec![],
            Some(diagnostic),
        )
        .expect("refused raw cell")
    }

    #[test]
    fn outcome_and_reason_codes_and_all_ninety_five_semantic_pairs_are_exact() {
        for (outcome, expected_code) in OUTCOME_CODES {
            assert_eq!(outcome.code(), expected_code);
        }
        for (reason, expected_code) in REASON_CODES {
            assert_eq!(reason.code(), expected_code);
        }

        for (outcome, _) in OUTCOME_CODES {
            for (reason, _) in REASON_CODES {
                let result = RunnerV2RawCellObservationV1::new(
                    token("semantic-matrix-cell"),
                    outcome,
                    reason,
                    vec![],
                    diagnostic_for_outcome(outcome),
                );
                if reason_is_compatible(outcome, reason) {
                    let cell = result.expect("compatible outcome and reason pair");
                    assert_eq!(cell.cell_id(), &token("semantic-matrix-cell"));
                    assert_eq!(cell.outcome(), outcome);
                    assert_eq!(cell.reason(), reason);
                    assert!(cell.numeric().is_empty());
                    assert_eq!(
                        cell.diagnostic().is_some(),
                        !matches!(outcome, RunnerV2RawOutcomeKindV1::Accepted)
                    );
                } else {
                    assert_error(
                        result.expect_err("incompatible outcome and reason pair"),
                        ConstructionErrorKindV2::Incompatible,
                        "runner_v2.handoff.cell.outcome_reason",
                        &reason.code().to_string(),
                    );
                }
            }
        }
    }

    #[test]
    fn diagnostic_code_and_retryability_matrices_are_exhaustive() {
        let cases = [
            (
                RunnerV2RawOutcomeKindV1::Refused,
                DiagnosticCodeV2::RunnerRefused,
                RetryabilityV2::AfterInputChange,
            ),
            (
                RunnerV2RawOutcomeKindV1::Failed,
                DiagnosticCodeV2::RunnerInternalError,
                RetryabilityV2::Never,
            ),
            (
                RunnerV2RawOutcomeKindV1::Unsupported,
                DiagnosticCodeV2::RunnerUnsupported,
                RetryabilityV2::AfterPrerequisiteChange,
            ),
            (
                RunnerV2RawOutcomeKindV1::Inapplicable,
                DiagnosticCodeV2::RunnerNotRun,
                RetryabilityV2::Never,
            ),
        ];

        for (outcome, expected_code, expected_retryability) in cases {
            let reason = semantic_reason(outcome);
            let prerequisites = || {
                if matches!(
                    outcome,
                    RunnerV2RawOutcomeKindV1::Unsupported | RunnerV2RawOutcomeKindV1::Inapplicable
                ) {
                    vec![token("registered-applicability-prerequisite")]
                } else {
                    vec![]
                }
            };

            for code in DiagnosticCodeV2::ALL {
                let result = RunnerV2RawCellObservationV1::new(
                    token("diagnostic-code-matrix-cell"),
                    outcome,
                    reason,
                    vec![],
                    Some(diagnostic_with(
                        code,
                        expected_retryability,
                        prerequisites(),
                        vec![],
                    )),
                );
                if code == expected_code {
                    assert_eq!(
                        result
                            .expect("exact diagnostic code")
                            .diagnostic()
                            .expect("non-accepted diagnostic")
                            .code(),
                        code
                    );
                } else {
                    assert_error(
                        result.expect_err("incompatible diagnostic code"),
                        ConstructionErrorKindV2::Incompatible,
                        "runner_v2.handoff.cell.diagnostic_code",
                        &code.code().to_string(),
                    );
                }
            }

            for retryability in RetryabilityV2::ALL {
                let result = RunnerV2RawCellObservationV1::new(
                    token("retryability-matrix-cell"),
                    outcome,
                    reason,
                    vec![],
                    Some(diagnostic_with(
                        expected_code,
                        retryability,
                        prerequisites(),
                        vec![],
                    )),
                );
                if retryability == expected_retryability {
                    assert_eq!(
                        result
                            .expect("exact retryability")
                            .diagnostic()
                            .expect("non-accepted diagnostic")
                            .retryability(),
                        retryability
                    );
                } else {
                    assert_error(
                        result.expect_err("incompatible retryability"),
                        ConstructionErrorKindV2::Incompatible,
                        "runner_v2.handoff.cell.retryability",
                        &retryability.code().to_string(),
                    );
                }
            }

            assert_error(
                RunnerV2RawCellObservationV1::new(
                    token("missing-diagnostic-cell"),
                    outcome,
                    reason,
                    vec![],
                    None,
                )
                .expect_err("every non-accepted outcome requires a diagnostic"),
                ConstructionErrorKindV2::Missing,
                "runner_v2.handoff.cell.diagnostic",
                "false",
            );
        }

        assert_error(
            RunnerV2RawCellObservationV1::new(
                token("unexpected-diagnostic-cell"),
                RunnerV2RawOutcomeKindV1::Accepted,
                RunnerV2RawReasonV1::ExactCheckedValue,
                vec![],
                diagnostic_for_outcome(RunnerV2RawOutcomeKindV1::Refused),
            )
            .expect_err("accepted outcomes cannot carry diagnostics"),
            ConstructionErrorKindV2::Unexpected,
            "runner_v2.handoff.cell.diagnostic",
            "true",
        );
    }

    #[test]
    fn prerequisite_and_repair_semantics_are_outcome_exact() {
        for (outcome, code, retryability) in [
            (
                RunnerV2RawOutcomeKindV1::Refused,
                DiagnosticCodeV2::RunnerRefused,
                RetryabilityV2::AfterInputChange,
            ),
            (
                RunnerV2RawOutcomeKindV1::Failed,
                DiagnosticCodeV2::RunnerInternalError,
                RetryabilityV2::Never,
            ),
            (
                RunnerV2RawOutcomeKindV1::Unsupported,
                DiagnosticCodeV2::RunnerUnsupported,
                RetryabilityV2::AfterPrerequisiteChange,
            ),
        ] {
            let prerequisites = if matches!(outcome, RunnerV2RawOutcomeKindV1::Unsupported) {
                vec![token("registered-applicability-prerequisite")]
            } else {
                vec![]
            };
            let cell = RunnerV2RawCellObservationV1::new(
                token("repair-permitted-cell"),
                outcome,
                semantic_reason(outcome),
                vec![],
                Some(diagnostic_with(
                    code,
                    retryability,
                    prerequisites,
                    ranked_repairs(1),
                )),
            )
            .expect("repairs are permitted for refused, failed, and unsupported outcomes");
            assert_eq!(
                cell.diagnostic()
                    .expect("non-accepted diagnostic")
                    .repairs()
                    .len(),
                1
            );
        }

        for outcome in [
            RunnerV2RawOutcomeKindV1::Unsupported,
            RunnerV2RawOutcomeKindV1::Inapplicable,
        ] {
            let (code, retryability) = match outcome {
                RunnerV2RawOutcomeKindV1::Unsupported => (
                    DiagnosticCodeV2::RunnerUnsupported,
                    RetryabilityV2::AfterPrerequisiteChange,
                ),
                RunnerV2RawOutcomeKindV1::Inapplicable => {
                    (DiagnosticCodeV2::RunnerNotRun, RetryabilityV2::Never)
                }
                _ => unreachable!("loop contains only applicability outcomes"),
            };
            assert_error(
                RunnerV2RawCellObservationV1::new(
                    token("missing-applicability-prerequisite-cell"),
                    outcome,
                    semantic_reason(outcome),
                    vec![],
                    Some(diagnostic_with(code, retryability, vec![], vec![])),
                )
                .expect_err("applicability outcomes require a prerequisite"),
                ConstructionErrorKindV2::Missing,
                "runner_v2.handoff.cell.applicability_prerequisite",
                "0",
            );
        }

        assert_error(
            RunnerV2RawCellObservationV1::new(
                token("inapplicable-repair-cell"),
                RunnerV2RawOutcomeKindV1::Inapplicable,
                RunnerV2RawReasonV1::PureDeclarationFacet,
                vec![],
                Some(diagnostic_with(
                    DiagnosticCodeV2::RunnerNotRun,
                    RetryabilityV2::Never,
                    vec![token("registered-applicability-prerequisite")],
                    ranked_repairs(1),
                )),
            )
            .expect_err("permanently inapplicable facets have no repair"),
            ConstructionErrorKindV2::Unexpected,
            "runner_v2.handoff.cell.inapplicable_repairs",
            "1",
        );
    }

    #[test]
    fn all_collection_bounds_accept_zero_and_max_and_refuse_one_over() {
        assert_error(
            RunnerV2LocalWorkPackageHandoffV1::new(token("empty-declaration-package"), &[], vec![])
                .expect_err("a handoff requires a nonempty declaration"),
            ConstructionErrorKindV2::Missing,
            "runner_v2.handoff.declared_cells",
            "0",
        );
        assert_error(
            RunnerV2LocalWorkPackageHandoffV1::new(
                token("empty-cell-package"),
                &[token("declared-cell")],
                vec![],
            )
            .expect_err("a nonempty declaration requires all cells"),
            ConstructionErrorKindV2::Missing,
            "runner_v2.handoff.cells",
            "0",
        );

        let maximum_declared = ordered_tokens("maximum-cell", RUNNER_V2_LOCAL_HANDOFF_MAX_CELLS_V1);
        let maximum_cells = maximum_declared
            .iter()
            .cloned()
            .map(accepted_token)
            .collect();
        let maximum_handoff = RunnerV2LocalWorkPackageHandoffV1::new(
            token("maximum-cell-package"),
            &maximum_declared,
            maximum_cells,
        )
        .expect("the exact cell maximum is admitted");
        assert_eq!(
            maximum_handoff.cells().len(),
            RUNNER_V2_LOCAL_HANDOFF_MAX_CELLS_V1
        );

        let too_many_declared = ordered_tokens(
            "declared-overflow-cell",
            RUNNER_V2_LOCAL_HANDOFF_MAX_CELLS_V1 + 1,
        );
        assert_error(
            RunnerV2LocalWorkPackageHandoffV1::new(
                token("declared-overflow-package"),
                &too_many_declared,
                vec![],
            )
            .expect_err("one-over declared cells are refused before allocation"),
            ConstructionErrorKindV2::TooLarge,
            "runner_v2.handoff.declared_cells",
            &(RUNNER_V2_LOCAL_HANDOFF_MAX_CELLS_V1 + 1).to_string(),
        );

        let too_many_cells =
            ordered_tokens("maximum-cell", RUNNER_V2_LOCAL_HANDOFF_MAX_CELLS_V1 + 1)
                .into_iter()
                .map(accepted_token)
                .collect();
        assert_error(
            RunnerV2LocalWorkPackageHandoffV1::new(
                token("cell-overflow-package"),
                &maximum_declared,
                too_many_cells,
            )
            .expect_err("one-over raw cells are unexpected"),
            ConstructionErrorKindV2::Unexpected,
            "runner_v2.handoff.cells",
            &(RUNNER_V2_LOCAL_HANDOFF_MAX_CELLS_V1 + 1).to_string(),
        );

        assert!(
            accepted_with_numeric(token("empty-numeric-cell"), vec![])
                .numeric()
                .is_empty()
        );
        let maximum_numeric = accepted_with_numeric(
            token("maximum-numeric-cell"),
            ordered_numeric(RUNNER_V2_LOCAL_HANDOFF_MAX_NUMERIC_OBSERVATIONS_V1),
        );
        assert_eq!(
            maximum_numeric.numeric().len(),
            RUNNER_V2_LOCAL_HANDOFF_MAX_NUMERIC_OBSERVATIONS_V1
        );
        assert_error(
            RunnerV2RawCellObservationV1::new(
                token("numeric-overflow-cell"),
                RunnerV2RawOutcomeKindV1::Accepted,
                RunnerV2RawReasonV1::ExactCheckedValue,
                ordered_numeric(RUNNER_V2_LOCAL_HANDOFF_MAX_NUMERIC_OBSERVATIONS_V1 + 1),
                None,
            )
            .expect_err("one-over numeric observations are refused"),
            ConstructionErrorKindV2::TooLarge,
            "runner_v2.handoff.cell.numeric",
            &(RUNNER_V2_LOCAL_HANDOFF_MAX_NUMERIC_OBSERVATIONS_V1 + 1).to_string(),
        );

        let empty_diagnostic = diagnostic_with(
            DiagnosticCodeV2::RunnerRefused,
            RetryabilityV2::AfterInputChange,
            vec![],
            vec![],
        );
        assert!(empty_diagnostic.prerequisites().is_empty());
        assert!(empty_diagnostic.repairs().is_empty());

        let maximum_diagnostic = diagnostic_with(
            DiagnosticCodeV2::RunnerRefused,
            RetryabilityV2::AfterInputChange,
            ordered_tokens("prerequisite", RUNNER_V2_LOCAL_HANDOFF_MAX_PREREQUISITES_V1),
            ranked_repairs(RUNNER_V2_LOCAL_HANDOFF_MAX_REPAIRS_V1),
        );
        assert_eq!(
            maximum_diagnostic.prerequisites().len(),
            RUNNER_V2_LOCAL_HANDOFF_MAX_PREREQUISITES_V1
        );
        assert_eq!(
            maximum_diagnostic.repairs().len(),
            RUNNER_V2_LOCAL_HANDOFF_MAX_REPAIRS_V1
        );

        assert_error(
            RunnerV2RawDiagnosticV1::new(
                DiagnosticCodeV2::RunnerRefused,
                token("prerequisite-overflow-owner"),
                RetryabilityV2::AfterInputChange,
                ordered_tokens(
                    "prerequisite",
                    RUNNER_V2_LOCAL_HANDOFF_MAX_PREREQUISITES_V1 + 1,
                ),
                vec![],
            )
            .expect_err("one-over prerequisites are refused"),
            ConstructionErrorKindV2::TooLarge,
            "runner_v2.handoff.diagnostic.prerequisites",
            &(RUNNER_V2_LOCAL_HANDOFF_MAX_PREREQUISITES_V1 + 1).to_string(),
        );
        assert_error(
            RunnerV2RawDiagnosticV1::new(
                DiagnosticCodeV2::RunnerRefused,
                token("repair-overflow-owner"),
                RetryabilityV2::AfterInputChange,
                vec![],
                ranked_repairs(RUNNER_V2_LOCAL_HANDOFF_MAX_REPAIRS_V1 + 1),
            )
            .expect_err("one-over repairs are refused"),
            ConstructionErrorKindV2::TooLarge,
            "runner_v2.handoff.diagnostic.repairs",
            &(RUNNER_V2_LOCAL_HANDOFF_MAX_REPAIRS_V1 + 1).to_string(),
        );
    }

    #[test]
    fn constructor_and_revalidation_exact_join_errors_report_first_mismatch_indexes() {
        let declared = [token("cell-z"), token("cell-a"), token("cell-m")];
        let handoff = RunnerV2LocalWorkPackageHandoffV1::new(
            token("runner-v2.work-package.24-1-1-1-1.v1"),
            &declared,
            vec![accepted("cell-z"), refused("cell-a"), accepted("cell-m")],
        )
        .expect("source order is authoritative and need not be lexical");
        assert_eq!(
            handoff.package_id(),
            &token("runner-v2.work-package.24-1-1-1-1.v1")
        );
        assert_eq!(
            handoff
                .cells()
                .iter()
                .map(RunnerV2RawCellObservationV1::cell_id)
                .collect::<Vec<_>>(),
            declared.iter().collect::<Vec<_>>()
        );
        handoff
            .validate_exact_cell_order(&declared)
            .expect("exact order");

        assert_error(
            RunnerV2LocalWorkPackageHandoffV1::new(
                token("missing-cell-package"),
                &declared,
                vec![accepted("cell-z"), refused("cell-a")],
            )
            .expect_err("a short cell join is missing"),
            ConstructionErrorKindV2::Missing,
            "runner_v2.handoff.cells",
            "2",
        );
        assert_error(
            RunnerV2LocalWorkPackageHandoffV1::new(
                token("unexpected-cell-package"),
                &declared,
                vec![
                    accepted("cell-z"),
                    refused("cell-a"),
                    accepted("cell-m"),
                    accepted("cell-extra"),
                ],
            )
            .expect_err("a long cell join is unexpected"),
            ConstructionErrorKindV2::Unexpected,
            "runner_v2.handoff.cells",
            "4",
        );
        assert_error(
            RunnerV2LocalWorkPackageHandoffV1::new(
                token("out-of-order-cell-package"),
                &declared,
                vec![accepted("cell-z"), accepted("cell-m"), refused("cell-a")],
            )
            .expect_err("a declared ID in the wrong position is out of order"),
            ConstructionErrorKindV2::OutOfOrder,
            "runner_v2.handoff.cells",
            "1",
        );
        assert_error(
            RunnerV2LocalWorkPackageHandoffV1::new(
                token("incompatible-cell-package"),
                &declared,
                vec![
                    accepted("cell-z"),
                    accepted("cell-unknown"),
                    accepted("cell-m"),
                ],
            )
            .expect_err("an undeclared ID is incompatible"),
            ConstructionErrorKindV2::Incompatible,
            "runner_v2.handoff.cells",
            "1",
        );
        assert_error(
            RunnerV2LocalWorkPackageHandoffV1::new(
                token("duplicate-declaration-package"),
                &[token("cell-z"), token("cell-a"), token("cell-z")],
                vec![accepted("cell-z"), accepted("cell-a"), accepted("cell-z")],
            )
            .expect_err("declaration IDs are globally unique"),
            ConstructionErrorKindV2::Duplicate,
            "runner_v2.handoff.declared_cells",
            "2",
        );

        assert_error(
            handoff
                .validate_exact_cell_order(&[token("cell-z"), token("cell-a")])
                .expect_err("a short expected join is unexpected relative to raw cells"),
            ConstructionErrorKindV2::Unexpected,
            "runner_v2.handoff.expected_cells",
            "2",
        );
        assert_error(
            handoff
                .validate_exact_cell_order(&[
                    token("cell-z"),
                    token("cell-a"),
                    token("cell-m"),
                    token("cell-extra"),
                ])
                .expect_err("a long expected join is missing raw cells"),
            ConstructionErrorKindV2::Missing,
            "runner_v2.handoff.expected_cells",
            "4",
        );
        assert_error(
            handoff
                .validate_exact_cell_order(&[token("cell-z"), token("cell-m"), token("cell-a")])
                .expect_err("a known expected ID in the wrong position is out of order"),
            ConstructionErrorKindV2::OutOfOrder,
            "runner_v2.handoff.expected_cells",
            "1",
        );
        assert_error(
            handoff
                .validate_exact_cell_order(&[
                    token("cell-z"),
                    token("cell-unknown"),
                    token("cell-m"),
                ])
                .expect_err("an expected set omitting the actual ID is incompatible"),
            ConstructionErrorKindV2::Incompatible,
            "runner_v2.handoff.expected_cells",
            "1",
        );
    }

    #[test]
    fn safe_numeric_union_units_and_accessors_round_trip_without_erasure() {
        let physical_unit = UnitV2::from_parts(1, 1, [1, 0, 0, 0, 0, 0, 0]).expect("physical unit");
        let cell = accepted_with_numeric(
            token("safe-numeric-cell"),
            vec![
                RunnerV2SafeNumericObservationV1::count(token("count"), 7),
                RunnerV2SafeNumericObservationV1::limit(
                    token("limit"),
                    RunnerLimitValueV2::U64(1_024),
                    RunnerLimitUnitV2::EncodedBytes,
                ),
                RunnerV2SafeNumericObservationV1::numeric(
                    token("logical"),
                    NumericValueV2::U32(5),
                    LogicalUnitV2::Dimensionless,
                ),
                RunnerV2SafeNumericObservationV1::physical(
                    token("physical"),
                    NumericValueV2::I64(-3),
                    physical_unit,
                ),
            ],
        );
        let numeric = cell.numeric();

        assert_eq!(numeric[0].name(), &token("count"));
        assert_eq!(numeric[0].value(), &RunnerV2SafeNumericValueV1::Count(7));
        assert_eq!(
            numeric[0].unit(),
            RunnerV2SafeNumericUnitV1::Logical(LogicalUnitV2::Count)
        );
        assert_eq!(numeric[1].name(), &token("limit"));
        assert_eq!(
            numeric[1].value(),
            &RunnerV2SafeNumericValueV1::Limit(RunnerLimitValueV2::U64(1_024))
        );
        assert_eq!(
            numeric[1].unit(),
            RunnerV2SafeNumericUnitV1::Limit(RunnerLimitUnitV2::EncodedBytes)
        );
        assert_eq!(numeric[2].name(), &token("logical"));
        assert_eq!(
            numeric[2].value(),
            &RunnerV2SafeNumericValueV1::Numeric(NumericValueV2::U32(5))
        );
        assert_eq!(
            numeric[2].unit(),
            RunnerV2SafeNumericUnitV1::Logical(LogicalUnitV2::Dimensionless)
        );
        assert_eq!(numeric[3].name(), &token("physical"));
        assert_eq!(
            numeric[3].value(),
            &RunnerV2SafeNumericValueV1::Numeric(NumericValueV2::I64(-3))
        );
        assert_eq!(
            numeric[3].unit(),
            RunnerV2SafeNumericUnitV1::Physical(physical_unit)
        );
    }

    #[test]
    fn numeric_and_prerequisite_ordering_refuses_duplicates_and_descents() {
        for (numeric, kind) in [
            (
                vec![
                    RunnerV2SafeNumericObservationV1::count(token("numeric-a"), 1),
                    RunnerV2SafeNumericObservationV1::count(token("numeric-a"), 2),
                ],
                ConstructionErrorKindV2::Duplicate,
            ),
            (
                vec![
                    RunnerV2SafeNumericObservationV1::count(token("numeric-b"), 1),
                    RunnerV2SafeNumericObservationV1::count(token("numeric-a"), 2),
                ],
                ConstructionErrorKindV2::OutOfOrder,
            ),
        ] {
            assert_error(
                RunnerV2RawCellObservationV1::new(
                    token("numeric-order-cell"),
                    RunnerV2RawOutcomeKindV1::Accepted,
                    RunnerV2RawReasonV1::ExactCheckedValue,
                    numeric,
                    None,
                )
                .expect_err("numeric names require strict order and uniqueness"),
                kind,
                "runner_v2.handoff.cell.numeric",
                "2",
            );
        }

        for (prerequisites, kind) in [
            (
                vec![token("prerequisite-a"), token("prerequisite-a")],
                ConstructionErrorKindV2::Duplicate,
            ),
            (
                vec![token("prerequisite-b"), token("prerequisite-a")],
                ConstructionErrorKindV2::OutOfOrder,
            ),
        ] {
            assert_error(
                RunnerV2RawDiagnosticV1::new(
                    DiagnosticCodeV2::RunnerRefused,
                    token("prerequisite-order-owner"),
                    RetryabilityV2::AfterInputChange,
                    prerequisites,
                    vec![],
                )
                .expect_err("prerequisites require strict order and uniqueness"),
                kind,
                "runner_v2.handoff.diagnostic.prerequisites",
                "2",
            );
        }
    }

    #[test]
    fn repair_ranks_duplicate_identity_and_accessors_are_exact() {
        for (ranks, observed) in [
            (&[0_u8][..], "0"),
            (&[2_u8][..], "2"),
            (&[1_u8, 1][..], "1"),
            (&[1_u8, 3][..], "3"),
            (&[2_u8, 1][..], "2"),
            (&[1_u8, 3, 2][..], "3"),
        ] {
            let repairs = ranks
                .iter()
                .enumerate()
                .map(|(index, rank)| {
                    RunnerV2RawRepairV1::new(
                        *rank,
                        RepairActionKindV2::ChangeArguments,
                        token(&format!("rank-target-{index:02}")),
                    )
                })
                .collect();
            assert_error(
                RunnerV2RawDiagnosticV1::new(
                    DiagnosticCodeV2::RunnerRefused,
                    token("repair-rank-owner"),
                    RetryabilityV2::AfterInputChange,
                    vec![],
                    repairs,
                )
                .expect_err("repair ranks must be contiguous and one based"),
                ConstructionErrorKindV2::OutOfRange,
                "runner_v2.handoff.diagnostic.repair_rank",
                observed,
            );
        }

        assert_error(
            RunnerV2RawDiagnosticV1::new(
                DiagnosticCodeV2::RunnerRefused,
                token("duplicate-repair-owner"),
                RetryabilityV2::AfterInputChange,
                vec![],
                vec![
                    RunnerV2RawRepairV1::new(
                        1,
                        RepairActionKindV2::ChangeArguments,
                        token("same-target"),
                    ),
                    RunnerV2RawRepairV1::new(
                        2,
                        RepairActionKindV2::ChangeArguments,
                        token("same-target"),
                    ),
                ],
            )
            .expect_err("the same repair kind and target is a duplicate"),
            ConstructionErrorKindV2::Duplicate,
            "runner_v2.handoff.diagnostic.repairs",
            "1",
        );

        let distinct = RunnerV2RawDiagnosticV1::new(
            DiagnosticCodeV2::RunnerRefused,
            token("distinct-repair-owner"),
            RetryabilityV2::AfterInputChange,
            vec![],
            vec![
                RunnerV2RawRepairV1::new(
                    1,
                    RepairActionKindV2::ChangeArguments,
                    token("first-target"),
                ),
                RunnerV2RawRepairV1::new(
                    2,
                    RepairActionKindV2::ChangeArguments,
                    token("second-target"),
                ),
                RunnerV2RawRepairV1::new(
                    3,
                    RepairActionKindV2::SupplyEvidence,
                    token("second-target"),
                ),
            ],
        )
        .expect("same kind with a new target and new kind with same target are distinct");
        assert_eq!(distinct.code(), DiagnosticCodeV2::RunnerRefused);
        assert_eq!(distinct.owner(), &token("distinct-repair-owner"));
        assert_eq!(distinct.retryability(), RetryabilityV2::AfterInputChange);
        assert!(distinct.prerequisites().is_empty());
        assert_eq!(distinct.repairs().len(), 3);
        assert_eq!(distinct.repairs()[0].rank(), 1);
        assert_eq!(
            distinct.repairs()[0].kind(),
            RepairActionKindV2::ChangeArguments
        );
        assert_eq!(distinct.repairs()[0].target(), &token("first-target"));
        assert_eq!(distinct.repairs()[1].rank(), 2);
        assert_eq!(distinct.repairs()[1].target(), &token("second-target"));
        assert_eq!(distinct.repairs()[2].rank(), 3);
        assert_eq!(
            distinct.repairs()[2].kind(),
            RepairActionKindV2::SupplyEvidence
        );
        assert_eq!(distinct.repairs()[2].target(), &token("second-target"));
    }

    #[test]
    fn repair_action_catalog_is_fully_usable_as_non_executable_data() {
        for (index, kind) in RepairActionKindV2::ALL.into_iter().enumerate() {
            let repair = RunnerV2RawRepairV1::new(1, kind, token("catalog-repair-target"));
            let diagnostic = diagnostic_with(
                DiagnosticCodeV2::RunnerRefused,
                RetryabilityV2::AfterInputChange,
                vec![],
                vec![repair],
            );
            assert_eq!(diagnostic.repairs()[0].rank(), 1);
            assert_eq!(diagnostic.repairs()[0].kind(), kind);
            assert_eq!(
                diagnostic.repairs()[0].target(),
                &token("catalog-repair-target")
            );
            assert_eq!(
                kind.code(),
                u16::try_from(index + 1).expect("catalog index")
            );
        }
    }
}
