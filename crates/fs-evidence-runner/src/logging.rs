//! Deterministic, bounded base-projection logs returned as typed data.

use crate::construction::{ConstructionErrorKindV2, ConstructionErrorV2};
use crate::path::LogicalBundlePathV1;
use crate::value::{StableTokenV2, TypedValueV2};
use std::collections::BTreeSet;

/// Maximum typed detail fields in one base E2E log event.
pub const BASE_E2E_LOG_FIELDS_MAX_V1: usize = 64;
/// Maximum symbolic reproduction arguments in one base E2E log event.
pub const BASE_E2E_REPRO_ARGS_MAX_V1: usize = 32;

/// Closed event kind for base-projection execution.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum BaseE2eLogKindV1 {
    /// A journey is about to evaluate its frozen row set.
    JourneyStart,
    /// One frozen row reached a deterministic terminal decision.
    CaseTerminal,
    /// A journey emitted its exact eligible/pass/fail/unsupported counts.
    JourneySummary,
    /// The complete five-journey projection emitted its aggregate counts.
    ProjectionSummary,
}

/// Closed outcome vocabulary for projection rows.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum BaseE2eOutcomeV1 {
    /// Constructor or validator agreed with its independent expected result.
    Passed,
    /// Constructor or validator disagreed with its independent expected result.
    Failed,
    /// A deliberately platform-dependent cell is not locally adjudicable.
    Unsupported,
    /// Start/summary event with no case outcome.
    NotApplicable,
}

/// One canonical named detail value.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseE2eLogFieldV1 {
    name: StableTokenV2,
    value: TypedValueV2,
}

impl BaseE2eLogFieldV1 {
    /// Bind a validated stable field name to one typed value.
    #[must_use]
    pub const fn new(name: StableTokenV2, value: TypedValueV2) -> Self {
        Self { name, value }
    }

    /// Stable field name.
    #[must_use]
    pub const fn name(&self) -> &StableTokenV2 {
        &self.name
    }

    /// Typed field value.
    #[must_use]
    pub const fn value(&self) -> &TypedValueV2 {
        &self.value
    }
}

/// A reproduction argument that cannot contain an ambient absolute path,
/// credential, shell fragment, or live process selector.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum SymbolicReproductionArgV1 {
    /// Symbolic workspace root supplied by the downstream harness.
    WorkspaceRoot,
    /// Symbolic source snapshot root supplied by the downstream harness.
    SourceSnapshot,
    /// Exact validated semantic argument.
    Literal(StableTokenV2),
}

/// One deterministic, bounded event.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseE2eLogEventV1 {
    logical_sequence: u32,
    journey: StableTokenV2,
    case: Option<StableTokenV2>,
    kind: BaseE2eLogKindV1,
    outcome: BaseE2eOutcomeV1,
    fields: Box<[BaseE2eLogFieldV1]>,
    relative_artifact: Option<LogicalBundlePathV1>,
    reproduction: Box<[SymbolicReproductionArgV1]>,
}

impl BaseE2eLogEventV1 {
    /// Validate one event and canonicalize the nonsemantic detail-field set.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        logical_sequence: u32,
        journey: StableTokenV2,
        case: Option<StableTokenV2>,
        kind: BaseE2eLogKindV1,
        outcome: BaseE2eOutcomeV1,
        mut fields: Vec<BaseE2eLogFieldV1>,
        relative_artifact: Option<LogicalBundlePathV1>,
        reproduction: Vec<SymbolicReproductionArgV1>,
    ) -> Result<Self, ConstructionErrorV2> {
        if fields.len() > BASE_E2E_LOG_FIELDS_MAX_V1 {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                "base_e2e_log.fields",
                "at most 64 typed fields",
                fields.len(),
            ));
        }
        if reproduction.len() > BASE_E2E_REPRO_ARGS_MAX_V1 {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                "base_e2e_log.reproduction",
                "at most 32 symbolic arguments",
                reproduction.len(),
            ));
        }
        let mut seen = BTreeSet::new();
        for field in &fields {
            if is_forbidden_ambient_name(field.name.as_str()) {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Incompatible,
                    "base_e2e_log.field_name",
                    "a declared semantic field, never ambient or sensitive process state",
                    field.name.as_str(),
                ));
            }
            if !seen.insert(field.name.as_str()) {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Duplicate,
                    "base_e2e_log.fields",
                    "unique stable field names",
                    field.name.as_str(),
                ));
            }
        }
        fields.sort_by(|left, right| left.name.cmp(&right.name));
        for argument in &reproduction {
            if let SymbolicReproductionArgV1::Literal(value) = argument
                && is_forbidden_ambient_name(value.as_str())
            {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Incompatible,
                    "base_e2e_log.reproduction",
                    "symbolic roots or declared non-sensitive semantic literals",
                    value.as_str(),
                ));
            }
        }

        let shape_valid = match kind {
            BaseE2eLogKindV1::JourneyStart => {
                case.is_none() && outcome == BaseE2eOutcomeV1::NotApplicable
            }
            BaseE2eLogKindV1::CaseTerminal => {
                case.is_some()
                    && matches!(
                        outcome,
                        BaseE2eOutcomeV1::Passed
                            | BaseE2eOutcomeV1::Failed
                            | BaseE2eOutcomeV1::Unsupported
                    )
            }
            BaseE2eLogKindV1::JourneySummary | BaseE2eLogKindV1::ProjectionSummary => {
                case.is_none() && outcome == BaseE2eOutcomeV1::NotApplicable
            }
        };
        if !shape_valid {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Incompatible,
                "base_e2e_log.kind_case_outcome",
                "the closed event-shape matrix",
                format_args!("{kind:?}/{}/{outcome:?}", case.is_some()),
            ));
        }

        Ok(Self {
            logical_sequence,
            journey,
            case,
            kind,
            outcome,
            fields: fields.into_boxed_slice(),
            relative_artifact,
            reproduction: reproduction.into_boxed_slice(),
        })
    }

    /// Globally contiguous deterministic sequence.
    #[must_use]
    pub const fn logical_sequence(&self) -> u32 {
        self.logical_sequence
    }

    /// Stable journey key.
    #[must_use]
    pub const fn journey(&self) -> &StableTokenV2 {
        &self.journey
    }

    /// Stable case key when this is a case event.
    #[must_use]
    pub const fn case(&self) -> Option<&StableTokenV2> {
        self.case.as_ref()
    }

    /// Event kind.
    #[must_use]
    pub const fn kind(&self) -> BaseE2eLogKindV1 {
        self.kind
    }

    /// Case outcome, or `NotApplicable` for start/summary events.
    #[must_use]
    pub const fn outcome(&self) -> BaseE2eOutcomeV1 {
        self.outcome
    }

    /// Canonically ordered typed detail fields.
    #[must_use]
    pub fn fields(&self) -> &[BaseE2eLogFieldV1] {
        &self.fields
    }

    /// Optional retained artifact path, always logical and relative.
    #[must_use]
    pub const fn relative_artifact(&self) -> Option<&LogicalBundlePathV1> {
        self.relative_artifact.as_ref()
    }

    /// Symbolic, non-shell reproduction arguments.
    #[must_use]
    pub fn reproduction(&self) -> &[SymbolicReproductionArgV1] {
        &self.reproduction
    }
}

fn is_forbidden_ambient_name(value: &str) -> bool {
    matches!(
        value,
        "pid"
            | "process-id"
            | "wall-clock"
            | "wall-clock-timestamp"
            | "timestamp"
            | "scheduler-latency"
            | "environment-secret"
            | "environment-value"
            | "secret"
            | "credential"
            | "physical-path"
            | "absolute-path"
            | "ambient-path"
            | "raw-payload"
    )
}

/// A validated contiguous event document.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseE2eLogV1 {
    events: Box<[BaseE2eLogEventV1]>,
}

impl BaseE2eLogV1 {
    /// Validate nonempty, zero-based contiguous sequence and a final aggregate
    /// summary.
    pub fn new(events: Vec<BaseE2eLogEventV1>) -> Result<Self, ConstructionErrorV2> {
        if events.is_empty() {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Missing,
                "base_e2e_log.events",
                "at least one event",
                0,
            ));
        }
        for (expected, event) in events.iter().enumerate() {
            let expected = u32::try_from(expected).map_err(|_| {
                ConstructionErrorV2::new(
                    ConstructionErrorKindV2::TooLarge,
                    "base_e2e_log.logical_sequence",
                    "event ordinal representable as u32",
                    expected,
                )
            })?;
            if event.logical_sequence != expected {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::OutOfOrder,
                    "base_e2e_log.logical_sequence",
                    "zero-based contiguous sequence",
                    event.logical_sequence,
                ));
            }
        }
        if events.last().map(BaseE2eLogEventV1::kind) != Some(BaseE2eLogKindV1::ProjectionSummary) {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Missing,
                "base_e2e_log.projection_summary",
                "exactly one final projection summary",
                "absent or nonterminal",
            ));
        }
        if events[..events.len() - 1]
            .iter()
            .any(|event| event.kind == BaseE2eLogKindV1::ProjectionSummary)
        {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Duplicate,
                "base_e2e_log.projection_summary",
                "only the final event is the projection summary",
                "early duplicate",
            ));
        }
        Ok(Self {
            events: events.into_boxed_slice(),
        })
    }

    /// Contiguous typed events.
    #[must_use]
    pub fn events(&self) -> &[BaseE2eLogEventV1] {
        &self.events
    }
}

#[cfg(test)]
mod tests {
    use super::{
        BASE_E2E_LOG_FIELDS_MAX_V1, BASE_E2E_REPRO_ARGS_MAX_V1, BaseE2eLogEventV1,
        BaseE2eLogFieldV1, BaseE2eLogKindV1, BaseE2eLogV1, BaseE2eOutcomeV1,
        SymbolicReproductionArgV1,
    };
    use crate::catalog::DigestRoleV2;
    use crate::identity::SourceIdentityRootV2;
    use crate::path::LogicalBundlePathV1;
    use crate::value::{NumericValueV2, QuantityV2, StableTokenV2, TypedValueV2, UnitV2};

    fn token(value: &str) -> StableTokenV2 {
        StableTokenV2::new(value).expect("fixture token")
    }

    fn event(
        sequence: u32,
        kind: BaseE2eLogKindV1,
        outcome: BaseE2eOutcomeV1,
        case: Option<&str>,
        fields: Vec<BaseE2eLogFieldV1>,
        reproduction: Vec<SymbolicReproductionArgV1>,
    ) -> Result<BaseE2eLogEventV1, crate::ConstructionErrorV2> {
        BaseE2eLogEventV1::new(
            sequence,
            token("publication-state"),
            case.map(token),
            kind,
            outcome,
            fields,
            Some(LogicalBundlePathV1::new("evidence/case.log").expect("relative artifact")),
            reproduction,
        )
    }

    #[test]
    fn field_order_is_canonical_and_duplicates_refuse() {
        let event = BaseE2eLogEventV1::new(
            0,
            token("publication-state"),
            None,
            BaseE2eLogKindV1::JourneyStart,
            BaseE2eOutcomeV1::NotApplicable,
            vec![
                BaseE2eLogFieldV1::new(token("z"), TypedValueV2::U8(2)),
                BaseE2eLogFieldV1::new(token("a"), TypedValueV2::U8(1)),
            ],
            None,
            vec![SymbolicReproductionArgV1::WorkspaceRoot],
        )
        .expect("valid event");
        assert_eq!(event.fields()[0].name().as_str(), "a");
        assert_eq!(event.fields()[1].name().as_str(), "z");

        assert!(
            BaseE2eLogEventV1::new(
                0,
                token("publication-state"),
                None,
                BaseE2eLogKindV1::JourneyStart,
                BaseE2eOutcomeV1::NotApplicable,
                vec![
                    BaseE2eLogFieldV1::new(token("a"), TypedValueV2::U8(1)),
                    BaseE2eLogFieldV1::new(token("a"), TypedValueV2::U8(2)),
                ],
                None,
                Vec::new(),
            )
            .is_err()
        );
    }

    #[test]
    fn log_sequence_and_terminal_summary_are_exact() {
        let start = BaseE2eLogEventV1::new(
            0,
            token("publication-state"),
            None,
            BaseE2eLogKindV1::JourneyStart,
            BaseE2eOutcomeV1::NotApplicable,
            Vec::new(),
            None,
            Vec::new(),
        )
        .expect("start");
        let summary = BaseE2eLogEventV1::new(
            1,
            token("all"),
            None,
            BaseE2eLogKindV1::ProjectionSummary,
            BaseE2eOutcomeV1::NotApplicable,
            Vec::new(),
            None,
            Vec::new(),
        )
        .expect("summary");
        assert!(BaseE2eLogV1::new(vec![start.clone(), summary]).is_ok());
        assert!(BaseE2eLogV1::new(vec![start]).is_err());
    }

    #[test]
    fn event_shape_matrix_and_collection_bounds_are_exact() {
        let shapes = [
            (
                BaseE2eLogKindV1::JourneyStart,
                BaseE2eOutcomeV1::NotApplicable,
                None,
                true,
            ),
            (
                BaseE2eLogKindV1::CaseTerminal,
                BaseE2eOutcomeV1::Passed,
                Some("case"),
                true,
            ),
            (
                BaseE2eLogKindV1::CaseTerminal,
                BaseE2eOutcomeV1::Failed,
                Some("case"),
                true,
            ),
            (
                BaseE2eLogKindV1::CaseTerminal,
                BaseE2eOutcomeV1::Unsupported,
                Some("case"),
                true,
            ),
            (
                BaseE2eLogKindV1::JourneySummary,
                BaseE2eOutcomeV1::NotApplicable,
                None,
                true,
            ),
            (
                BaseE2eLogKindV1::ProjectionSummary,
                BaseE2eOutcomeV1::NotApplicable,
                None,
                true,
            ),
            (
                BaseE2eLogKindV1::JourneyStart,
                BaseE2eOutcomeV1::Passed,
                None,
                false,
            ),
            (
                BaseE2eLogKindV1::CaseTerminal,
                BaseE2eOutcomeV1::NotApplicable,
                Some("case"),
                false,
            ),
            (
                BaseE2eLogKindV1::JourneySummary,
                BaseE2eOutcomeV1::NotApplicable,
                Some("case"),
                false,
            ),
        ];
        for (kind, outcome, case, expected) in shapes {
            assert_eq!(
                event(0, kind, outcome, case, Vec::new(), Vec::new()).is_ok(),
                expected,
                "{kind:?}/{outcome:?}/{case:?}"
            );
        }

        let exact_fields = (0..BASE_E2E_LOG_FIELDS_MAX_V1)
            .map(|index| {
                BaseE2eLogFieldV1::new(
                    token(&format!("field.{index}")),
                    TypedValueV2::U32(u32::try_from(index).expect("bounded index")),
                )
            })
            .collect();
        assert!(
            event(
                0,
                BaseE2eLogKindV1::JourneyStart,
                BaseE2eOutcomeV1::NotApplicable,
                None,
                exact_fields,
                Vec::new()
            )
            .is_ok()
        );
        let one_over_fields = (0..=BASE_E2E_LOG_FIELDS_MAX_V1)
            .map(|index| {
                BaseE2eLogFieldV1::new(
                    token(&format!("field.{index}")),
                    TypedValueV2::U32(u32::try_from(index).expect("bounded index")),
                )
            })
            .collect();
        assert!(
            event(
                0,
                BaseE2eLogKindV1::JourneyStart,
                BaseE2eOutcomeV1::NotApplicable,
                None,
                one_over_fields,
                Vec::new()
            )
            .is_err()
        );
        let exact_reproduction = (0..BASE_E2E_REPRO_ARGS_MAX_V1)
            .map(|index| SymbolicReproductionArgV1::Literal(token(&format!("case.{index}"))))
            .collect();
        assert!(
            event(
                0,
                BaseE2eLogKindV1::JourneyStart,
                BaseE2eOutcomeV1::NotApplicable,
                None,
                Vec::new(),
                exact_reproduction
            )
            .is_ok()
        );
        let one_over_reproduction = (0..=BASE_E2E_REPRO_ARGS_MAX_V1)
            .map(|index| SymbolicReproductionArgV1::Literal(token(&format!("case.{index}"))))
            .collect();
        assert!(
            event(
                0,
                BaseE2eLogKindV1::JourneyStart,
                BaseE2eOutcomeV1::NotApplicable,
                None,
                Vec::new(),
                one_over_reproduction
            )
            .is_err()
        );
    }

    #[test]
    fn detailed_semantic_fields_are_typed_deterministic_and_nonambient() {
        let source = SourceIdentityRootV2::parse_presented(
            DigestRoleV2::Source,
            SourceIdentityRootV2::DESCRIPTOR.domain(),
            &"11".repeat(32),
        )
        .expect("fixture source root");
        let seconds = UnitV2::from_parts(1, 1, [0, 0, 1, 0, 0, 0, 0]).expect("seconds unit");
        let fields = vec![
            BaseE2eLogFieldV1::new(token("api-generation"), TypedValueV2::U16(2)),
            BaseE2eLogFieldV1::new(
                token("source-root"),
                TypedValueV2::Digest(source.digest().clone()),
            ),
            BaseE2eLogFieldV1::new(token("case-ordinal"), TypedValueV2::U32(7)),
            BaseE2eLogFieldV1::new(
                token("duration"),
                TypedValueV2::Quantity(QuantityV2::new(NumericValueV2::U64(3), seconds)),
            ),
            BaseE2eLogFieldV1::new(token("owner"), TypedValueV2::Token(token("runner.owner"))),
        ];
        let reproduction = vec![
            SymbolicReproductionArgV1::WorkspaceRoot,
            SymbolicReproductionArgV1::SourceSnapshot,
            SymbolicReproductionArgV1::Literal(token("case")),
        ];
        let first = event(
            0,
            BaseE2eLogKindV1::CaseTerminal,
            BaseE2eOutcomeV1::Passed,
            Some("case"),
            fields.clone(),
            reproduction.clone(),
        )
        .expect("typed deterministic event");
        let second = event(
            0,
            BaseE2eLogKindV1::CaseTerminal,
            BaseE2eOutcomeV1::Passed,
            Some("case"),
            fields,
            reproduction,
        )
        .expect("same deterministic event");
        assert_eq!(first, second);
        assert_eq!(
            first
                .relative_artifact()
                .expect("relative artifact")
                .as_str(),
            "evidence/case.log"
        );

        for forbidden in [
            "pid",
            "process-id",
            "wall-clock-timestamp",
            "timestamp",
            "scheduler-latency",
            "environment-secret",
            "credential",
            "physical-path",
            "absolute-path",
            "raw-payload",
        ] {
            assert!(
                event(
                    0,
                    BaseE2eLogKindV1::JourneyStart,
                    BaseE2eOutcomeV1::NotApplicable,
                    None,
                    vec![BaseE2eLogFieldV1::new(
                        token(forbidden),
                        TypedValueV2::U8(1)
                    )],
                    Vec::new()
                )
                .is_err(),
                "{forbidden}"
            );
            assert!(
                event(
                    0,
                    BaseE2eLogKindV1::JourneyStart,
                    BaseE2eOutcomeV1::NotApplicable,
                    None,
                    Vec::new(),
                    vec![SymbolicReproductionArgV1::Literal(token(forbidden))]
                )
                .is_err(),
                "reproduction {forbidden}"
            );
        }
        assert!(LogicalBundlePathV1::new("/absolute/evidence.log").is_err());
    }
}
