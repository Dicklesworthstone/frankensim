//! Expansion-program risks from the new-domains charter.
//!
//! This register is deliberately separate from the addendum's R1--R10
//! scientific/engineering register in the crate root.  PR-001--PR-012 govern
//! the expansion program itself: ownership, leading indicators, quantitative
//! trip points, mitigations, contingencies, and phase-exit review gates.

use core::fmt::Write as _;

use crate::json_escape;

/// Stable identifier for one expansion-program risk.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ProgramRiskId {
    /// Validated-flow wrapping and certificate cost explosion.
    Pr001,
    /// Contact candidate and certificate explosion.
    Pr002,
    /// Hiptmair--Xu robustness does not materialize.
    Pr003,
    /// Chemistry stiffness overwhelms the integrator lanes.
    Pr004,
    /// Material-data licensing or experimental scarcity.
    Pr005,
    /// Event-bearing hybrid adjoints are invalid.
    Pr006,
    /// Same-layer Cargo cycles emerge from the expanded crate atlas.
    Pr007,
    /// Schema or hash migration errors during the six-base transition.
    Pr008,
    /// Theorem-tool trusted-computing-base failure.
    Pr009,
    /// Single-node scale ceiling conflicts with accelerator policy.
    Pr010,
    /// Interoperability adoption fails.
    Pr011,
    /// Scientific receipts are misused as safety/regulatory certificates.
    Pr012,
}

impl ProgramRiskId {
    /// Every program risk in canonical order.
    pub const ALL: [Self; 12] = [
        Self::Pr001,
        Self::Pr002,
        Self::Pr003,
        Self::Pr004,
        Self::Pr005,
        Self::Pr006,
        Self::Pr007,
        Self::Pr008,
        Self::Pr009,
        Self::Pr010,
        Self::Pr011,
        Self::Pr012,
    ];

    /// Stable externally visible code.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::Pr001 => "PR-001",
            Self::Pr002 => "PR-002",
            Self::Pr003 => "PR-003",
            Self::Pr004 => "PR-004",
            Self::Pr005 => "PR-005",
            Self::Pr006 => "PR-006",
            Self::Pr007 => "PR-007",
            Self::Pr008 => "PR-008",
            Self::Pr009 => "PR-009",
            Self::Pr010 => "PR-010",
            Self::Pr011 => "PR-011",
            Self::Pr012 => "PR-012",
        }
    }
}

/// Five-point governance rating used for likelihood and impact.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum RiskRating {
    /// Very low.
    VeryLow = 1,
    /// Low.
    Low = 2,
    /// Moderate.
    Moderate = 3,
    /// High.
    High = 4,
    /// Very high.
    VeryHigh = 5,
}

impl RiskRating {
    /// Numeric 1--5 representation used by the canonical artifact.
    #[must_use]
    pub const fn score(self) -> u8 {
        self as u8
    }
}

/// Phase exit at which a program risk must be reviewed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReviewGate {
    /// E0a constitutional/schema exit.
    E0a,
    /// E0b architecture/manifest exit.
    E0b,
    /// E0c data-readiness exit.
    E0c,
    /// E0d theorem-foundry exit.
    E0d,
    /// E2 validated geometry/contact exit.
    E2,
    /// E4 field-solver exit.
    E4,
    /// E5 reacting-flow exit.
    E5,
    /// E6 scale-qualification exit.
    E6,
    /// E7 workflow/assurance exit.
    E7,
}

impl ReviewGate {
    /// Stable phase-exit spelling.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::E0a => "E0a",
            Self::E0b => "E0b",
            Self::E0c => "E0c",
            Self::E0d => "E0d",
            Self::E2 => "E2",
            Self::E4 => "E4",
            Self::E5 => "E5",
            Self::E6 => "E6",
            Self::E7 => "E7",
        }
    }
}

/// Comparator for a numeric leading-indicator trip point.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TriggerComparison {
    /// Trip when the observation is at least the threshold.
    GreaterThanOrEqual,
    /// Trip when the observation is strictly below the threshold.
    LessThan,
}

/// Admissible numeric domain for a leading-indicator aggregate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ObservationDomain {
    /// Any finite value greater than or equal to zero.
    NonNegativeReal,
    /// A finite fraction in the closed interval `[0, 1]`.
    Fraction,
    /// A finite, non-negative, exactly integral count.
    NonNegativeInteger,
}

impl ObservationDomain {
    /// Stable artifact spelling.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::NonNegativeReal => "non-negative-real",
            Self::Fraction => "fraction-0-to-1",
            Self::NonNegativeInteger => "non-negative-integer",
        }
    }

    fn admits(self, value: f64) -> bool {
        match self {
            Self::NonNegativeReal => value >= 0.0,
            Self::Fraction => (0.0..=1.0).contains(&value),
            Self::NonNegativeInteger => value >= 0.0 && value.fract() == 0.0,
        }
    }
}

impl TriggerComparison {
    /// Stable operator spelling for artifacts and diagnostics.
    #[must_use]
    pub const fn symbol(self) -> &'static str {
        match self {
            Self::GreaterThanOrEqual => ">=",
            Self::LessThan => "<",
        }
    }

    /// Apply the comparator to one finite observation.
    #[must_use]
    pub fn is_triggered(self, observed: f64, threshold: f64) -> bool {
        match self {
            Self::GreaterThanOrEqual => observed >= threshold,
            Self::LessThan => observed < threshold,
        }
    }
}

/// Quantitative trip point for one leading indicator.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct QuantitativeTrigger {
    /// Comparison applied to the aggregate observation.
    pub comparison: TriggerComparison,
    /// Finite numeric threshold.
    pub threshold: f64,
    /// Explicit unit for the threshold and observation.
    pub unit: &'static str,
    /// Admissible numeric domain for the aggregate.
    pub domain: ObservationDomain,
    /// Minimum number of raw samples supporting the aggregate observation.
    pub min_samples: usize,
}

impl QuantitativeTrigger {
    /// Judge one explicitly unit-tagged aggregate observation.
    ///
    /// Missing observations are handled by [`assess_program_risks`]. This
    /// method fails closed on unit mismatch, non-finite values, and an
    /// insufficient sample count before applying the comparator.
    #[must_use]
    pub fn assess(self, observed: f64, samples: usize, unit: &str) -> AssessmentStatus {
        if unit != self.unit {
            AssessmentStatus::UnitMismatch
        } else if !observed.is_finite() {
            AssessmentStatus::NonFinite
        } else if !self.domain.admits(observed) {
            AssessmentStatus::OutOfRange
        } else if samples < self.min_samples {
            AssessmentStatus::UnderSampled
        } else if self.comparison.is_triggered(observed, self.threshold) {
            AssessmentStatus::Triggered
        } else {
            AssessmentStatus::Clear
        }
    }
}

/// Named owner for a program-risk workstream.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ProgramRiskOwner {
    /// Stable workstream role.
    pub role: &'static str,
    /// Bead that owns mitigation and contingency work.
    pub bead_id: &'static str,
}

/// One complete expansion-program risk row.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ProgramRisk {
    /// Stable id.
    pub id: ProgramRiskId,
    /// Short human-readable name.
    pub name: &'static str,
    /// Named workstream and owning Bead.
    pub owner: ProgramRiskOwner,
    /// Initial likelihood rating.
    pub likelihood: RiskRating,
    /// Initial impact rating.
    pub impact: RiskRating,
    /// Machine-facing leading-indicator name.
    pub leading_indicator: &'static str,
    /// Numeric trigger, unit, and sampling floor.
    pub trigger: QuantitativeTrigger,
    /// Action intended to reduce likelihood or impact before the trigger.
    pub mitigation: &'static str,
    /// Kill, refusal, fallback, or escalation action after the trigger.
    pub contingency: &'static str,
    /// Likelihood rating after the mitigation is applied.
    pub residual_likelihood: RiskRating,
    /// Impact rating after the mitigation is applied.
    pub residual_impact: RiskRating,
    /// Phase exit that must review this row.
    pub review_gate: ReviewGate,
}

const PROGRAM_RISKS: [ProgramRisk; 12] = [
    ProgramRisk {
        id: ProgramRiskId::Pr001,
        name: "Validated-flow certificate cost explosion",
        owner: ProgramRiskOwner {
            role: "validated-step",
            bead_id: "frankensim-ext-time-validated-step-ow2o",
        },
        likelihood: RiskRating::High,
        impact: RiskRating::High,
        leading_indicator: "certified_flow_wall_time_ratio_p95",
        trigger: QuantitativeTrigger {
            comparison: TriggerComparison::GreaterThanOrEqual,
            threshold: 10.0,
            unit: "ratio-to-estimated-baseline",
            domain: ObservationDomain::NonNegativeReal,
            min_samples: 5,
        },
        mitigation: "Reduce wrapping through representation changes and adaptive refinement before promoting the validated lane.",
        contingency: "Disable the validated lane by default for the affected regime and return Estimated or Unknown evidence color.",
        residual_likelihood: RiskRating::Low,
        residual_impact: RiskRating::Moderate,
        review_gate: ReviewGate::E2,
    },
    ProgramRisk {
        id: ProgramRiskId::Pr002,
        name: "Contact candidate and certificate explosion",
        owner: ProgramRiskOwner {
            role: "contact-detection",
            bead_id: "frankensim-ext-contact-detection-ccd-tqag",
        },
        likelihood: RiskRating::High,
        impact: RiskRating::High,
        leading_indicator: "contact_candidate_certificates_per_accepted_pair_p95",
        trigger: QuantitativeTrigger {
            comparison: TriggerComparison::GreaterThanOrEqual,
            threshold: 100.0,
            unit: "candidate-certificates-per-accepted-pair",
            domain: ObservationDomain::NonNegativeReal,
            min_samples: 10,
        },
        mitigation: "Partition and refine the broad phase while retaining every conservative candidate needed by the certificate contract.",
        contingency: "Cap or refuse unresolved contacts with an explicit no-claim result; never silently prune candidates.",
        residual_likelihood: RiskRating::Low,
        residual_impact: RiskRating::Moderate,
        review_gate: ReviewGate::E2,
    },
    ProgramRisk {
        id: ProgramRiskId::Pr003,
        name: "HX preconditioner robustness failure",
        owner: ProgramRiskOwner {
            role: "hx-preconditioner",
            bead_id: "frankensim-ext-solver-hx-preconditioner-12l2",
        },
        likelihood: RiskRating::Moderate,
        impact: RiskRating::VeryHigh,
        leading_indicator: "hx_krylov_iterations_worst",
        trigger: QuantitativeTrigger {
            comparison: TriggerComparison::GreaterThanOrEqual,
            threshold: 200.0,
            unit: "iterations",
            domain: ObservationDomain::NonNegativeInteger,
            min_samples: 12,
        },
        mitigation: "Tune auxiliary and coarse spaces across the full coefficient/topology robustness matrix.",
        contingency: "Restrict admitted contrasts or topologies and route to a policy-admitted direct or AMG fallback.",
        residual_likelihood: RiskRating::Low,
        residual_impact: RiskRating::High,
        review_gate: ReviewGate::E4,
    },
    ProgramRisk {
        id: ProgramRiskId::Pr004,
        name: "Chemistry stiffness overwhelms integrator lanes",
        owner: ProgramRiskOwner {
            role: "chemistry-ladder",
            bead_id: "frankensim-ext-gas-chemistry-ladder-paqh",
        },
        likelihood: RiskRating::High,
        impact: RiskRating::High,
        leading_indicator: "chemistry_substeps_per_flow_step_p95",
        trigger: QuantitativeTrigger {
            comparison: TriggerComparison::GreaterThanOrEqual,
            threshold: 1_000.0,
            unit: "chemistry-substeps-per-flow-step",
            domain: ObservationDomain::NonNegativeReal,
            min_samples: 100,
        },
        mitigation: "Select IMEX, global, or skeletal chemistry according to measured stiffness and error receipts.",
        contingency: "Kill the detailed chemistry lane for the affected regime and report the admitted reduced model explicitly.",
        residual_likelihood: RiskRating::Low,
        residual_impact: RiskRating::Moderate,
        review_gate: ReviewGate::E5,
    },
    ProgramRisk {
        id: ProgramRiskId::Pr005,
        name: "Material-data licensing and scarcity",
        owner: ProgramRiskOwner {
            role: "material-dataset",
            bead_id: "frankensim-ext-matdb-seed-dataset-1sxe",
        },
        likelihood: RiskRating::High,
        impact: RiskRating::High,
        leading_indicator: "redistributable_required_claim_coverage",
        trigger: QuantitativeTrigger {
            comparison: TriggerComparison::LessThan,
            threshold: 0.70,
            unit: "fraction",
            domain: ObservationDomain::Fraction,
            min_samples: 10,
        },
        mitigation: "Source, cite, license, and uncertainty-qualify the required material claims before dataset promotion.",
        contingency: "Cut scope or quarantine non-redistributable local data behind explicit no-claim boundaries.",
        residual_likelihood: RiskRating::Moderate,
        residual_impact: RiskRating::Moderate,
        review_gate: ReviewGate::E0c,
    },
    ProgramRisk {
        id: ProgramRiskId::Pr006,
        name: "Hybrid-adjoint invalidity",
        owner: ProgramRiskOwner {
            role: "adjoint-composition",
            bead_id: "frankensim-ext-adjoint-composition-easb",
        },
        likelihood: RiskRating::Moderate,
        impact: RiskRating::VeryHigh,
        leading_indicator: "event_bearing_adjoint_fd_pass_fraction",
        trigger: QuantitativeTrigger {
            comparison: TriggerComparison::LessThan,
            threshold: 0.80,
            unit: "fraction",
            domain: ObservationDomain::Fraction,
            min_samples: 10,
        },
        mitigation: "Repair event derivatives and validate composed adjoints against independent finite-difference fixtures.",
        contingency: "Disable the classical hybrid adjoint and route to generalized or derivative-free optimization.",
        residual_likelihood: RiskRating::Low,
        residual_impact: RiskRating::High,
        review_gate: ReviewGate::E2,
    },
    ProgramRisk {
        id: ProgramRiskId::Pr007,
        name: "Same-layer Cargo cycle",
        owner: ProgramRiskOwner {
            role: "manifest-fixture",
            bead_id: "frankensim-ext-manifest-fixture-r56j",
        },
        likelihood: RiskRating::Moderate,
        impact: RiskRating::VeryHigh,
        leading_indicator: "same_layer_cargo_cycle_count",
        trigger: QuantitativeTrigger {
            comparison: TriggerComparison::GreaterThanOrEqual,
            threshold: 1.0,
            unit: "cycles",
            domain: ObservationDomain::NonNegativeInteger,
            min_samples: 1,
        },
        mitigation: "Validate every proposed edge against a complete manifest fixture before changing the workspace graph.",
        contingency: "Move ownership or the shared API to a neutral lower layer and block the manifest phase exit.",
        residual_likelihood: RiskRating::VeryLow,
        residual_impact: RiskRating::High,
        review_gate: ReviewGate::E0b,
    },
    ProgramRisk {
        id: ProgramRiskId::Pr008,
        name: "Six-base schema or hash migration error",
        owner: ProgramRiskOwner {
            role: "ledger-migration",
            bead_id: "frankensim-ext-ledger-package-migration-h61n",
        },
        likelihood: RiskRating::Moderate,
        impact: RiskRating::VeryHigh,
        leading_indicator: "migration_replay_or_hash_mismatch_count",
        trigger: QuantitativeTrigger {
            comparison: TriggerComparison::GreaterThanOrEqual,
            threshold: 1.0,
            unit: "mismatches",
            domain: ObservationDomain::NonNegativeInteger,
            min_samples: 1,
        },
        mitigation: "Retain old identities and exercise known-version replay plus old/new crosswalks over the migration corpus.",
        contingency: "Halt cutover, preserve the old schema authority, and repair the migrator before any identity promotion.",
        residual_likelihood: RiskRating::VeryLow,
        residual_impact: RiskRating::VeryHigh,
        review_gate: ReviewGate::E0a,
    },
    ProgramRisk {
        id: ProgramRiskId::Pr009,
        name: "Theorem-tool trusted-computing-base failure",
        owner: ProgramRiskOwner {
            role: "theorem-foundry",
            bead_id: "frankensim-ext-theorem-foundry-infra-zxob",
        },
        likelihood: RiskRating::Low,
        impact: RiskRating::VeryHigh,
        leading_indicator: "formal_kernel_native_checker_disagreement_count",
        trigger: QuantitativeTrigger {
            comparison: TriggerComparison::GreaterThanOrEqual,
            threshold: 1.0,
            unit: "disagreements",
            domain: ObservationDomain::NonNegativeInteger,
            min_samples: 1,
        },
        mitigation: "Cross-check the formal kernel and native checker over a retained adversarial theorem corpus.",
        contingency: "Quarantine affected theorem cards, freeze the kernel/version, and block downstream promotion.",
        residual_likelihood: RiskRating::VeryLow,
        residual_impact: RiskRating::VeryHigh,
        review_gate: ReviewGate::E0d,
    },
    ProgramRisk {
        id: ProgramRiskId::Pr010,
        name: "Scale ceiling conflicts with accelerator policy",
        owner: ProgramRiskOwner {
            role: "scale-qualification",
            bead_id: "frankensim-ext-scale-qualification-0h2j",
        },
        likelihood: RiskRating::High,
        impact: RiskRating::High,
        leading_indicator: "flagship_budget_misses_without_admitted_accelerator",
        trigger: QuantitativeTrigger {
            comparison: TriggerComparison::GreaterThanOrEqual,
            threshold: 2.0,
            unit: "flagship-decks",
            domain: ObservationDomain::NonNegativeInteger,
            min_samples: 2,
        },
        mitigation: "Reduce fidelity or optimize the native substrate against matched-QoI scale receipts.",
        contingency: "Mark the capability unavailable and escalate any dependency-policy change through constitutional review.",
        residual_likelihood: RiskRating::Moderate,
        residual_impact: RiskRating::High,
        review_gate: ReviewGate::E6,
    },
    ProgramRisk {
        id: ProgramRiskId::Pr011,
        name: "Workflow interoperability adoption failure",
        owner: ProgramRiskOwner {
            role: "workflow-interop",
            bead_id: "frankensim-ext-workflow-interop-lz8f",
        },
        likelihood: RiskRating::Moderate,
        impact: RiskRating::High,
        leading_indicator: "partner_trial_acceptance_fraction",
        trigger: QuantitativeTrigger {
            comparison: TriggerComparison::LessThan,
            threshold: 0.25,
            unit: "fraction",
            domain: ObservationDomain::Fraction,
            min_samples: 8,
        },
        mitigation: "Rank partner refusal causes and address only the highest-evidence native workflow gaps.",
        contingency: "Stop broad adapter work and narrow the supported boundary to export-only or native workflows.",
        residual_likelihood: RiskRating::Low,
        residual_impact: RiskRating::Moderate,
        review_gate: ReviewGate::E7,
    },
    ProgramRisk {
        id: ProgramRiskId::Pr012,
        name: "Scientific receipts misrepresented as safety certificates",
        owner: ProgramRiskOwner {
            role: "safety-assurance",
            bead_id: "frankensim-ext-safety-emc-assurance-te0w",
        },
        likelihood: RiskRating::Moderate,
        impact: RiskRating::VeryHigh,
        leading_indicator: "safety_or_regulatory_misrepresentation_count",
        trigger: QuantitativeTrigger {
            comparison: TriggerComparison::GreaterThanOrEqual,
            threshold: 1.0,
            unit: "exported-reports",
            domain: ObservationDomain::NonNegativeInteger,
            min_samples: 1,
        },
        mitigation: "Keep scientific evidence color and no-claim language explicit in every exported assurance surface.",
        contingency: "Block export, relabel the artifact, and require independent conformance, authorization, and legal review.",
        residual_likelihood: RiskRating::VeryLow,
        residual_impact: RiskRating::VeryHigh,
        review_gate: ReviewGate::E7,
    },
];

/// Canonical PR-001--PR-012 register.
#[must_use]
pub fn program_risks() -> &'static [ProgramRisk] {
    &PROGRAM_RISKS
}

/// Look up one canonical row.
#[must_use]
pub fn program_risk(id: ProgramRiskId) -> &'static ProgramRisk {
    let index = ProgramRiskId::ALL
        .iter()
        .position(|candidate| *candidate == id)
        .expect("every ProgramRiskId is registered");
    &PROGRAM_RISKS[index]
}

/// One aggregate leading-indicator observation supplied at session end.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ProgramRiskObservation<'a> {
    /// Risk whose leading indicator was measured.
    pub id: ProgramRiskId,
    /// Aggregate value in the row's declared unit.
    pub value: f64,
    /// Explicit unit. It must exactly match the canonical row's unit.
    pub unit: &'a str,
    /// Number of raw samples supporting the aggregate.
    pub samples: usize,
}

/// Fail-closed disposition of one register row.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AssessmentStatus {
    /// One finite, sufficiently sampled observation did not trip the trigger.
    Clear,
    /// One finite, sufficiently sampled observation tripped the trigger.
    Triggered,
    /// No observation was supplied.
    Missing,
    /// More than one aggregate observation was supplied for the same id.
    Duplicate,
    /// The one supplied aggregate was NaN or infinite.
    NonFinite,
    /// The supplied unit did not exactly match the register row.
    UnitMismatch,
    /// The value was outside the row's typed numeric domain.
    OutOfRange,
    /// The aggregate did not meet the declared sample floor.
    UnderSampled,
}

impl AssessmentStatus {
    /// Stable artifact spelling.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::Clear => "clear",
            Self::Triggered => "triggered",
            Self::Missing => "missing",
            Self::Duplicate => "duplicate",
            Self::NonFinite => "non-finite",
            Self::UnitMismatch => "unit-mismatch",
            Self::OutOfRange => "out-of-range",
            Self::UnderSampled => "under-sampled",
        }
    }

    /// Only a measured, sufficiently sampled, untripped row is green.
    #[must_use]
    pub const fn is_alert(self) -> bool {
        !matches!(self, Self::Clear)
    }
}

/// Assessment of one canonical program-risk row.
#[derive(Debug, Clone, PartialEq)]
pub struct ProgramRiskAssessmentRow {
    /// Canonical id.
    pub id: ProgramRiskId,
    /// Fail-closed status.
    pub status: AssessmentStatus,
    /// Finite value when exactly one observation was supplied.
    pub observed_value: Option<f64>,
    /// Bounded UTF-8-safe preview of the supplied unit for one observation.
    pub observed_unit: Option<String>,
    /// Exact byte length of the supplied unit for one observation.
    pub observed_unit_bytes: Option<usize>,
    /// Sample count when exactly one observation was supplied.
    pub samples: Option<usize>,
    /// Number of aggregate observations supplied for this id.
    pub observation_count: usize,
}

/// Canonically ordered assessment of all twelve program risks.
#[derive(Debug, Clone, PartialEq)]
pub struct ProgramRiskAssessment {
    rows: Vec<ProgramRiskAssessmentRow>,
}

/// Maximum retained bytes from a caller-supplied observation unit.
pub const MAX_OBSERVED_UNIT_PREVIEW_BYTES: usize = 64;

fn bounded_unit_preview(unit: &str) -> String {
    let mut end = 0usize;
    for (offset, character) in unit.char_indices() {
        let next = offset + character.len_utf8();
        if next > MAX_OBSERVED_UNIT_PREVIEW_BYTES {
            break;
        }
        end = next;
    }
    unit[..end].to_owned()
}

impl ProgramRiskAssessment {
    /// All rows in PR-001--PR-012 order.
    #[must_use]
    pub fn rows(&self) -> &[ProgramRiskAssessmentRow] {
        &self.rows
    }

    /// Every non-green row. Missing or malformed evidence is an alert.
    pub fn alerts(&self) -> impl Iterator<Item = &ProgramRiskAssessmentRow> {
        self.rows.iter().filter(|row| row.status.is_alert())
    }

    /// Number of non-green rows.
    #[must_use]
    pub fn alert_count(&self) -> usize {
        self.alerts().count()
    }

    /// True only when all twelve rows have one finite, sufficiently sampled,
    /// untripped observation.
    #[must_use]
    pub fn all_clear(&self) -> bool {
        !self.rows.is_empty() && self.alert_count() == 0
    }

    /// Deterministic JSON assessment in canonical risk order.
    #[must_use]
    pub fn to_json(&self) -> String {
        let mut out =
            String::from("{\"schema\":\"frankensim.program-risk-assessment.v1\",\"all_clear\":");
        out.push_str(if self.all_clear() { "true" } else { "false" });
        write!(out, ",\"alert_count\":{},\"rows\":[", self.alert_count())
            .expect("writing to a String is infallible");
        for (index, row) in self.rows.iter().enumerate() {
            if index > 0 {
                out.push(',');
            }
            let risk = program_risk(row.id);
            write!(
                out,
                "{{\"id\":\"{}\",\"status\":\"{}\",\"observed_value\":{},\"observed_unit\":{},\"observed_unit_bytes\":{},\"samples\":{},\"observation_count\":{},\"trigger\":{{\"comparison\":\"{}\",\"threshold\":{},\"unit\":\"{}\",\"domain\":\"{}\",\"min_samples\":{}}}}}",
                row.id.code(),
                row.status.code(),
                row.observed_value.map_or_else(
                    || "null".to_owned(),
                    |value| value.to_string()
                ),
                row.observed_unit.as_ref().map_or_else(
                    || "null".to_owned(),
                    |unit| format!("\"{}\"", json_escape(unit))
                ),
                row.observed_unit_bytes
                    .map_or_else(|| "null".to_owned(), |bytes| bytes.to_string()),
                row.samples
                    .map_or_else(|| "null".to_owned(), |samples| samples.to_string()),
                row.observation_count,
                risk.trigger.comparison.symbol(),
                risk.trigger.threshold,
                json_escape(risk.trigger.unit),
                risk.trigger.domain.code(),
                risk.trigger.min_samples,
            )
            .expect("writing to a String is infallible");
        }
        out.push_str("]}");
        out
    }
}

/// Assess all twelve canonical rows.
///
/// Input order never affects the result. Missing, duplicate, non-finite, and
/// under-sampled observations are alerts rather than implicit green rows.
#[must_use]
pub fn assess_program_risks(observations: &[ProgramRiskObservation<'_>]) -> ProgramRiskAssessment {
    let mut rows = Vec::with_capacity(PROGRAM_RISKS.len());
    for risk in PROGRAM_RISKS {
        let mut matches = observations
            .iter()
            .filter(|observation| observation.id == risk.id);
        let first = matches.next();
        let extras = matches.count();
        let observation_count = usize::from(first.is_some()).saturating_add(extras);
        let (status, observed_value, observed_unit, observed_unit_bytes, samples) =
            match (first, extras) {
                (None, _) => (AssessmentStatus::Missing, None, None, None, None),
                (Some(_), 1..) => (AssessmentStatus::Duplicate, None, None, None, None),
                (Some(observation), 0) => {
                    let status = risk.trigger.assess(
                        observation.value,
                        observation.samples,
                        observation.unit,
                    );
                    (
                        status,
                        observation.value.is_finite().then_some(observation.value),
                        Some(bounded_unit_preview(observation.unit)),
                        Some(observation.unit.len()),
                        Some(observation.samples),
                    )
                }
            };
        rows.push(ProgramRiskAssessmentRow {
            id: risk.id,
            status,
            observed_value,
            observed_unit,
            observed_unit_bytes,
            samples,
            observation_count,
        });
    }
    ProgramRiskAssessment { rows }
}

/// Deterministic, schema-versioned register artifact.
#[must_use]
pub fn program_risk_register_json() -> String {
    let mut out = String::from("{\"schema\":\"frankensim.program-risk-register.v1\",\"risks\":[");
    for (index, risk) in PROGRAM_RISKS.iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        write!(
            out,
            "{{\"id\":\"{}\",\"name\":\"{}\",\"owner\":{{\"role\":\"{}\",\"bead_id\":\"{}\"}},\"likelihood\":{},\"impact\":{},\"leading_indicator\":\"{}\",\"trigger\":{{\"comparison\":\"{}\",\"threshold\":{},\"unit\":\"{}\",\"domain\":\"{}\",\"min_samples\":{}}},\"mitigation\":\"{}\",\"contingency\":\"{}\",\"residual_risk\":{{\"likelihood\":{},\"impact\":{}}},\"review_gate\":\"{}\"}}",
            risk.id.code(),
            json_escape(risk.name),
            json_escape(risk.owner.role),
            json_escape(risk.owner.bead_id),
            risk.likelihood.score(),
            risk.impact.score(),
            json_escape(risk.leading_indicator),
            risk.trigger.comparison.symbol(),
            risk.trigger.threshold,
            json_escape(risk.trigger.unit),
            risk.trigger.domain.code(),
            risk.trigger.min_samples,
            json_escape(risk.mitigation),
            json_escape(risk.contingency),
            risk.residual_likelihood.score(),
            risk.residual_impact.score(),
            risk.review_gate.code(),
        )
        .expect("writing to a String is infallible");
    }
    out.push_str("]}");
    out
}
