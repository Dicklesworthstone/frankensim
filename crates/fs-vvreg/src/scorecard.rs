//! Public V&V scorecard: outcome metrics as one deterministic artifact.
//!
//! The scorecard is a pure projection of (a) the registered validation
//! corpus and (b) caller-supplied ledgered run results plus adversarial
//! assessments. It reports outcome metrics per (QoI, regime) cell —
//! reference counts by portfolio axis, prediction-error rows with reference
//! uncertainty shown next to the model error, interval coverage,
//! false-acceptance counts, and regime-limitation statements — and renders
//! byte-identical Markdown and canonical JSON.
//!
//! Honesty rules:
//! - A cell without executed evidence renders `NO-DATA`, never zero. Zero
//!   executed adversarial challenges is `NO-DATA`, not a clean bill.
//! - Empirical interval coverage is `NO-DATA` until the e07 coverage
//!   machinery is live; nominal coverage is never extrapolated.
//! - Supplying run records or assessments grants no authority: the
//!   scorecard is a diagnostic projection, and every claim cap declared by
//!   the corpus remains in force. Pass/fail counts here are envelope
//!   arithmetic, not evidence colours.
//! - A caller-built corpus is labelled `caller-built` so a synthetic
//!   registry can never masquerade as the public seeded corpus.

use core::fmt;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;

use fs_blake3::{ContentHash, hash_domain};

use crate::adversarial::{
    ADVERSARIAL_CASE_SCHEMA_VERSION, AdversarialAssessment, AdversarialEvidence,
    AdversarialRegistry, AdversarialScorecardError, DominantUncertainty, HonestyVerdict,
};
use crate::corpus::{
    CORPUS_SCHEMA_VERSION, ContextRange, CorpusEnvelope, CorpusRegistry, EvidenceLevel,
    LEVEL_C_COOLING_QOIS, MAX_CORPUS_TEXT_BYTES,
};
use crate::portfolio::EvidenceAxis;

/// Canonical scorecard wire and identity schema.
pub const SCORECARD_SCHEMA_VERSION: u32 = 1;
/// Maximum ledgered run records accepted by one scorecard build.
pub const MAX_SCORECARD_RUN_RECORDS: usize = 4_096;

const SCORECARD_IDENTITY_DOMAIN: &str = "org.frankensim.fs-vvreg.vv-scorecard.v1";
const RUN_RECORD_IDENTITY_DOMAIN: &str = "org.frankensim.fs-vvreg.scorecard-run-record.v1";
const COVERAGE_NO_DATA_REASON: &str = "empirical interval-coverage machinery (e07) is not live; nominal coverage is never extrapolated";
const EMPTY_REGIME_LABEL: &str = "(no declared regime restriction)";
const UNREGISTERED_REGIME_LABEL: &str = "(no dataset registered)";

/// Portfolio axes counted as external reference coverage. Numerical
/// verification, transferability, and independent reproduction are portfolio
/// dimensions of the program itself, not independent external references.
pub const EXTERNAL_AXES: [EvidenceAxis; 4] = [
    EvidenceAxis::CrossCodeAgreement,
    EvidenceAxis::ControlledExperimentalValidation,
    EvidenceAxis::BlindPredictiveValidation,
    EvidenceAxis::FieldMonitoring,
];

/// Fail-closed scorecard build refusal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ScorecardError {
    /// Too many run records were supplied.
    ResourceLimit {
        /// Maximum accepted count.
        limit: usize,
        /// Supplied count.
        observed: usize,
    },
    /// A run-record field is invalid.
    InvalidRunField {
        /// Stable field name.
        field: &'static str,
        /// Stable reason.
        reason: &'static str,
    },
    /// A run record names a dataset the corpus does not contain.
    UnknownDataset {
        /// Requested dataset id.
        dataset_id: String,
    },
    /// A run record names a metric its dataset does not declare.
    UnknownMetric {
        /// Dataset id.
        dataset_id: String,
        /// Requested metric.
        metric: String,
    },
    /// An assessment belongs to another adversarial registry.
    ForeignAssessment {
        /// Assessment case id.
        case_id: String,
    },
    /// More than one assessment was supplied for one adversarial case.
    DuplicateAssessment {
        /// Duplicate case id.
        case_id: String,
    },
}

impl fmt::Display for ScorecardError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ResourceLimit { limit, observed } => {
                write!(
                    formatter,
                    "run-record count {observed} exceeds limit {limit}"
                )
            }
            Self::InvalidRunField { field, reason } => {
                write!(formatter, "run-record field `{field}` {reason}")
            }
            Self::UnknownDataset { dataset_id } => {
                write!(formatter, "run record names unknown dataset `{dataset_id}`")
            }
            Self::UnknownMetric { dataset_id, metric } => write!(
                formatter,
                "run record names metric `{metric}` not declared by dataset `{dataset_id}`"
            ),
            Self::ForeignAssessment { case_id } => write!(
                formatter,
                "assessment for `{case_id}` belongs to another adversarial registry"
            ),
            Self::DuplicateAssessment { case_id } => {
                write!(formatter, "duplicate assessment for `{case_id}`")
            }
        }
    }
}

impl std::error::Error for ScorecardError {}

impl From<AdversarialScorecardError> for ScorecardError {
    fn from(error: AdversarialScorecardError) -> Self {
        match error {
            AdversarialScorecardError::ForeignAssessment { case_id } => {
                Self::ForeignAssessment { case_id }
            }
            AdversarialScorecardError::DuplicateAssessment { case_id } => {
                Self::DuplicateAssessment { case_id }
            }
        }
    }
}

/// Reference uncertainty attached to one ledgered run comparison.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ReferenceUncertainty {
    /// Symmetric absolute half-width in the metric's dimensions.
    Bounded {
        /// Non-negative half-width.
        half_width: f64,
    },
    /// The retained reference states no quantitative uncertainty.
    Unstated,
}

/// One ledgered run result compared against one corpus reference metric.
///
/// This is an input record, not authority: it binds a solver run (through
/// its ledger identity) to a corpus dataset metric so the scorecard can
/// report the comparison. Admission validates arithmetic and text bounds;
/// dataset and metric existence are checked against the corpus at build.
#[derive(Debug, Clone, PartialEq)]
pub struct ScorecardRunRecord {
    dataset_id: String,
    metric: String,
    predicted: f64,
    reference: f64,
    reference_uncertainty: ReferenceUncertainty,
    run_identity: ContentHash,
}

impl ScorecardRunRecord {
    /// Admit one run record with fail-closed field validation.
    ///
    /// # Errors
    /// Refuses blank/oversized/control-character ids, non-finite values,
    /// and negative or non-finite uncertainty half-widths.
    pub fn try_new(
        dataset_id: &str,
        metric: &str,
        predicted: f64,
        reference: f64,
        reference_uncertainty: ReferenceUncertainty,
        run_identity: ContentHash,
    ) -> Result<Self, ScorecardError> {
        validate_run_text("dataset_id", dataset_id)?;
        validate_run_text("metric", metric)?;
        if !predicted.is_finite() {
            return Err(ScorecardError::InvalidRunField {
                field: "predicted",
                reason: "must be finite",
            });
        }
        if !reference.is_finite() {
            return Err(ScorecardError::InvalidRunField {
                field: "reference",
                reason: "must be finite",
            });
        }
        if let ReferenceUncertainty::Bounded { half_width } = reference_uncertainty
            && (!half_width.is_finite() || half_width < 0.0)
        {
            return Err(ScorecardError::InvalidRunField {
                field: "reference_uncertainty",
                reason: "half-width must be finite and non-negative",
            });
        }
        Ok(Self {
            dataset_id: dataset_id.to_string(),
            metric: metric.to_string(),
            predicted,
            reference,
            reference_uncertainty,
            run_identity,
        })
    }

    /// Compared dataset id.
    #[must_use]
    pub fn dataset_id(&self) -> &str {
        &self.dataset_id
    }

    /// Compared corpus metric.
    #[must_use]
    pub fn metric(&self) -> &str {
        &self.metric
    }

    /// Model prediction in the metric's dimensions.
    #[must_use]
    pub const fn predicted(&self) -> f64 {
        self.predicted
    }

    /// External reference value in the metric's dimensions.
    #[must_use]
    pub const fn reference(&self) -> f64 {
        self.reference
    }

    /// Declared reference uncertainty.
    #[must_use]
    pub const fn reference_uncertainty(&self) -> ReferenceUncertainty {
        self.reference_uncertainty
    }

    /// Ledger identity binding the producing run.
    #[must_use]
    pub const fn run_identity(&self) -> ContentHash {
        self.run_identity
    }

    /// Signed prediction error (`predicted - reference`).
    #[must_use]
    pub const fn signed_error(&self) -> f64 {
        self.predicted - self.reference
    }
}

fn validate_run_text(field: &'static str, value: &str) -> Result<(), ScorecardError> {
    if value.trim().is_empty() {
        return Err(ScorecardError::InvalidRunField {
            field,
            reason: "is blank",
        });
    }
    if value.len() > MAX_CORPUS_TEXT_BYTES {
        return Err(ScorecardError::InvalidRunField {
            field,
            reason: "exceeds the byte limit",
        });
    }
    if value.chars().any(char::is_control) {
        return Err(ScorecardError::InvalidRunField {
            field,
            reason: "contains a control character",
        });
    }
    Ok(())
}

/// Envelope arithmetic verdict for one run comparison. This is not an
/// evidence colour: an in-envelope comparison proves envelope arithmetic
/// only, and an `Unpinned` corpus envelope supports no verdict at all.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EnvelopeCheck {
    /// The prediction satisfies the pinned envelope.
    Pass,
    /// The prediction violates the pinned envelope.
    Fail,
    /// The dataset declares no defensible scalar envelope for this metric.
    Unpinned,
}

impl EnvelopeCheck {
    /// Stable scorecard spelling.
    #[must_use]
    pub const fn slug(self) -> &'static str {
        match self {
            Self::Pass => "pass",
            Self::Fail => "fail",
            Self::Unpinned => "unpinned",
        }
    }
}

/// One rendered run-comparison row inside a cell.
#[derive(Debug, Clone, PartialEq)]
pub struct CellRunRow {
    dataset_id: String,
    run_identity: ContentHash,
    predicted: f64,
    reference: f64,
    reference_uncertainty: ReferenceUncertainty,
    signed_error: f64,
    envelope: EnvelopeCheck,
}

impl CellRunRow {
    /// Compared dataset id.
    #[must_use]
    pub fn dataset_id(&self) -> &str {
        &self.dataset_id
    }

    /// Ledger identity binding the producing run.
    #[must_use]
    pub const fn run_identity(&self) -> ContentHash {
        self.run_identity
    }

    /// Model prediction.
    #[must_use]
    pub const fn predicted(&self) -> f64 {
        self.predicted
    }

    /// External reference value.
    #[must_use]
    pub const fn reference(&self) -> f64 {
        self.reference
    }

    /// Declared reference uncertainty.
    #[must_use]
    pub const fn reference_uncertainty(&self) -> ReferenceUncertainty {
        self.reference_uncertainty
    }

    /// Signed prediction error.
    #[must_use]
    pub const fn signed_error(&self) -> f64 {
        self.signed_error
    }

    /// Envelope arithmetic verdict.
    #[must_use]
    pub const fn envelope(&self) -> EnvelopeCheck {
        self.envelope
    }
}

/// False-acceptance report for one cell. Zero executed challenges is
/// `NoData`, never a zero count.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FalseAcceptanceCell {
    /// No adversarial challenge bound to this cell's datasets was executed.
    NoData,
    /// At least one bound challenge was executed.
    Counted {
        /// Executed (non-`NO-DATA`) bound challenges.
        executed: usize,
        /// Executed challenges that were false acceptances.
        false_acceptances: usize,
    },
}

/// One deterministic (QoI, regime) outcome cell.
#[derive(Debug, Clone, PartialEq)]
pub struct ScorecardCell {
    qoi: String,
    regime: String,
    regime_statement: String,
    dataset_ids: Vec<String>,
    axis_counts: [usize; EvidenceAxis::ALL.len()],
    external_datasets: usize,
    runs: Vec<CellRunRow>,
    envelope_pass: usize,
    envelope_fail: usize,
    envelope_unpinned: usize,
    false_acceptance: FalseAcceptanceCell,
}

impl ScorecardCell {
    /// Stable QoI/metric identifier.
    #[must_use]
    pub fn qoi(&self) -> &str {
        &self.qoi
    }

    /// Canonical regime label derived from the acceptance-record validity
    /// ranges.
    #[must_use]
    pub fn regime(&self) -> &str {
        &self.regime
    }

    /// Regime-limitation statement derived from the declared validity
    /// domain.
    #[must_use]
    pub fn regime_statement(&self) -> &str {
        &self.regime_statement
    }

    /// Sorted ids of corpus datasets declaring this (QoI, regime) cell.
    #[must_use]
    pub fn dataset_ids(&self) -> &[String] {
        &self.dataset_ids
    }

    /// Number of cell datasets supplying the named portfolio coordinate.
    #[must_use]
    pub const fn datasets_on_axis(&self, axis: EvidenceAxis) -> usize {
        self.axis_counts[axis.index()]
    }

    /// Number of cell datasets carrying at least one external axis (see
    /// [`EXTERNAL_AXES`]).
    #[must_use]
    pub const fn external_datasets(&self) -> usize {
        self.external_datasets
    }

    /// Sorted run-comparison rows.
    #[must_use]
    pub fn runs(&self) -> &[CellRunRow] {
        &self.runs
    }

    /// In-envelope run count (envelope arithmetic, not evidence colour).
    #[must_use]
    pub const fn envelope_pass(&self) -> usize {
        self.envelope_pass
    }

    /// Out-of-envelope run count.
    #[must_use]
    pub const fn envelope_fail(&self) -> usize {
        self.envelope_fail
    }

    /// Runs whose dataset metric declares no pinned envelope.
    #[must_use]
    pub const fn envelope_unpinned(&self) -> usize {
        self.envelope_unpinned
    }

    /// False-acceptance report for this cell.
    #[must_use]
    pub const fn false_acceptance(&self) -> FalseAcceptanceCell {
        self.false_acceptance
    }
}

/// Cell accumulator used during scorecard assembly.
struct CellDraft {
    dataset_ids: BTreeSet<String>,
    axis_counts: [usize; EvidenceAxis::ALL.len()],
    external_datasets: usize,
    runs: Vec<CellRunRow>,
}

/// One adversarial case row retained for the JSON render.
#[derive(Debug, Clone, PartialEq)]
struct AdversarialCaseRow {
    case_id: &'static str,
    regime: &'static str,
    attacked: &'static str,
    evidence: String,
    verdict: &'static str,
    dominant: &'static str,
    false_acceptance: Option<bool>,
    limitation: &'static str,
}

/// Deterministic public V&V scorecard.
///
/// Build one with [`build_scorecard`]; render it with
/// [`VvScorecard::render_markdown`] and [`VvScorecard::render_json`]. Both
/// renders are pure functions of the typed content, and the identity binds
/// the typed content rather than rendered bytes.
#[derive(Debug, Clone, PartialEq)]
pub struct VvScorecard {
    corpus_seeded: bool,
    corpus_digest: ContentHash,
    adversarial_registry_identity: ContentHash,
    dataset_total: usize,
    level_counts: [usize; 5],
    axis_totals: [usize; EvidenceAxis::ALL.len()],
    run_record_count: usize,
    supplied_assessment_count: usize,
    executed_assessment_count: usize,
    false_acceptance_total: usize,
    adversarial_case_count: usize,
    cells: Vec<ScorecardCell>,
    known_gaps: Vec<String>,
    adversarial_rows: Vec<AdversarialCaseRow>,
    adversarial_markdown: String,
    identity: ContentHash,
}

const fn level_index(level: EvidenceLevel) -> usize {
    match level {
        EvidenceLevel::Analytic => 0,
        EvidenceLevel::CrossCode => 1,
        EvidenceLevel::PublishedExperiment => 2,
        EvidenceLevel::Blind => 3,
        EvidenceLevel::Field => 4,
    }
}

const LEVEL_CODES: [&str; 5] = ["A", "B", "C", "D", "E"];

fn regime_label(ranges: &[ContextRange]) -> String {
    if ranges.is_empty() {
        return EMPTY_REGIME_LABEL.to_string();
    }
    let mut parts = Vec::with_capacity(ranges.len());
    for range in ranges {
        if range.lo.dims == range.hi.dims {
            parts.push(format!(
                "{} in [{}, {}] {}",
                range.name,
                range.lo.value,
                range.hi.value,
                range.lo.dims.unit_string()
            ));
        } else {
            parts.push(format!("{} in [{}, {}]", range.name, range.lo, range.hi));
        }
    }
    parts.join("; ")
}

fn envelope_check(envelope: &CorpusEnvelope, predicted: f64, reference: f64) -> EnvelopeCheck {
    match envelope {
        CorpusEnvelope::Tolerance { atol, rtol } => {
            if (predicted - reference).abs() <= atol + rtol * reference.abs() {
                EnvelopeCheck::Pass
            } else {
                EnvelopeCheck::Fail
            }
        }
        CorpusEnvelope::Interval { lo, hi } => {
            if predicted >= *lo && predicted <= *hi {
                EnvelopeCheck::Pass
            } else {
                EnvelopeCheck::Fail
            }
        }
        CorpusEnvelope::Unpinned { .. } => EnvelopeCheck::Unpinned,
    }
}

/// Build the scorecard as a pure projection of the supplied corpus, the
/// adversarial registry, ledgered run records, and executed assessments.
///
/// # Errors
/// Refuses run records naming unknown datasets or undeclared metrics,
/// resource-limit violations, and assessments that are foreign to or
/// duplicated within the supplied adversarial registry. Refusal is total:
/// no partial scorecard is produced.
#[allow(clippy::too_many_lines)] // The build is one linear, auditable assembly pass.
pub fn build_scorecard(
    corpus: &CorpusRegistry,
    adversarial: &AdversarialRegistry,
    run_records: &[ScorecardRunRecord],
    assessments: &[AdversarialAssessment],
) -> Result<VvScorecard, ScorecardError> {
    if run_records.len() > MAX_SCORECARD_RUN_RECORDS {
        return Err(ScorecardError::ResourceLimit {
            limit: MAX_SCORECARD_RUN_RECORDS,
            observed: run_records.len(),
        });
    }
    // The adversarial renderer performs the registry-binding and
    // per-case-uniqueness validation this build also relies on.
    let adversarial_markdown = adversarial.render_regime_limitations(assessments)?;

    let mut by_case: BTreeMap<&str, &AdversarialAssessment> = BTreeMap::new();
    for assessment in assessments {
        by_case.insert(assessment.case_id(), assessment);
    }

    let mut level_counts = [0_usize; 5];
    let mut axis_totals = [0_usize; EvidenceAxis::ALL.len()];
    for dataset in corpus.datasets() {
        level_counts[level_index(dataset.evidence_level())] += 1;
        for &axis in dataset.evidence_level().portfolio_axes() {
            axis_totals[axis.index()] += 1;
        }
    }

    // The regime label is value-injective: f64 `Display` is shortest
    // round-trip, so distinct bounds always yield distinct labels and the
    // (metric, label) key cannot merge different regimes.
    let mut cells: BTreeMap<(String, String), CellDraft> = BTreeMap::new();
    for dataset in corpus.datasets() {
        for record in dataset.acceptance_envelopes() {
            let key = (record.metric.clone(), regime_label(&record.regime));
            let draft = cells.entry(key).or_insert_with(|| CellDraft {
                dataset_ids: BTreeSet::new(),
                axis_counts: [0; EvidenceAxis::ALL.len()],
                external_datasets: 0,
                runs: Vec::new(),
            });
            if draft.dataset_ids.insert(dataset.id().to_string()) {
                for &axis in dataset.evidence_level().portfolio_axes() {
                    draft.axis_counts[axis.index()] += 1;
                }
                if dataset
                    .evidence_level()
                    .portfolio_axes()
                    .iter()
                    .any(|axis| EXTERNAL_AXES.contains(axis))
                {
                    draft.external_datasets += 1;
                }
            }
        }
    }

    for record in run_records {
        let Some(dataset) = corpus.dataset(record.dataset_id()) else {
            return Err(ScorecardError::UnknownDataset {
                dataset_id: record.dataset_id().to_string(),
            });
        };
        let Some(acceptance) = dataset
            .acceptance_envelopes()
            .iter()
            .find(|acceptance| acceptance.metric == record.metric())
        else {
            return Err(ScorecardError::UnknownMetric {
                dataset_id: record.dataset_id().to_string(),
                metric: record.metric().to_string(),
            });
        };
        let key = (acceptance.metric.clone(), regime_label(&acceptance.regime));
        let envelope = envelope_check(&acceptance.envelope, record.predicted(), record.reference());
        let draft = cells
            .get_mut(&key)
            .expect("every dataset acceptance record registered its cell");
        draft.runs.push(CellRunRow {
            dataset_id: record.dataset_id().to_string(),
            run_identity: record.run_identity(),
            predicted: record.predicted(),
            reference: record.reference(),
            reference_uncertainty: record.reference_uncertainty(),
            signed_error: record.signed_error(),
            envelope,
        });
    }

    // Attribute executed adversarial challenges to the cells of their
    // retained corpus datasets. Planned-evidence cases have no retained
    // dataset and therefore no cell; they stay visible in the adversarial
    // section as NO-DATA.
    let mut executed_by_dataset: BTreeMap<&str, (usize, usize)> = BTreeMap::new();
    let mut executed_assessment_count = 0_usize;
    let mut false_acceptance_total = 0_usize;
    for assessment in assessments {
        if assessment.verdict() != HonestyVerdict::NoData {
            executed_assessment_count += 1;
            if assessment.is_false_acceptance() {
                false_acceptance_total += 1;
            }
            let case = adversarial
                .cases()
                .iter()
                .find(|case| case.id == assessment.case_id())
                .expect("render_regime_limitations validated case membership");
            if let AdversarialEvidence::Retained { dataset_id } = case.evidence {
                let entry = executed_by_dataset.entry(dataset_id).or_insert((0, 0));
                entry.0 += 1;
                if assessment.is_false_acceptance() {
                    entry.1 += 1;
                }
            }
        }
    }

    let mut finished = Vec::with_capacity(cells.len());
    let mut known_gaps = Vec::new();
    for ((qoi, regime), mut draft) in cells {
        draft.runs.sort_by(|a, b| {
            (
                a.dataset_id.as_str(),
                a.run_identity.0,
                a.predicted.to_bits(),
                a.reference.to_bits(),
            )
                .cmp(&(
                    b.dataset_id.as_str(),
                    b.run_identity.0,
                    b.predicted.to_bits(),
                    b.reference.to_bits(),
                ))
        });
        let mut envelope_pass = 0;
        let mut envelope_fail = 0;
        let mut envelope_unpinned = 0;
        for run in &draft.runs {
            match run.envelope {
                EnvelopeCheck::Pass => envelope_pass += 1,
                EnvelopeCheck::Fail => envelope_fail += 1,
                EnvelopeCheck::Unpinned => envelope_unpinned += 1,
            }
        }
        let mut executed = 0;
        let mut false_acceptances = 0;
        for dataset_id in &draft.dataset_ids {
            if let Some((cell_executed, cell_false)) = executed_by_dataset.get(dataset_id.as_str())
            {
                executed += cell_executed;
                false_acceptances += cell_false;
            }
        }
        let false_acceptance = if executed == 0 {
            FalseAcceptanceCell::NoData
        } else {
            FalseAcceptanceCell::Counted {
                executed,
                false_acceptances,
            }
        };
        let regime_statement = if regime == EMPTY_REGIME_LABEL {
            format!(
                "The corpus declares no regime restriction for `{qoi}` beyond each dataset's context of use; absence of a restriction is not evidence of transferability."
            )
        } else {
            format!(
                "`{qoi}` references apply only within {regime}; outside this declared context the corpus asserts no claim."
            )
        };
        if draft.external_datasets == 0 {
            known_gaps.push(format!("qoi={qoi} regime={regime} external_datasets=0"));
        }
        finished.push(ScorecardCell {
            qoi,
            regime,
            regime_statement,
            dataset_ids: draft.dataset_ids.into_iter().collect(),
            axis_counts: draft.axis_counts,
            external_datasets: draft.external_datasets,
            runs: draft.runs,
            envelope_pass,
            envelope_fail,
            envelope_unpinned,
            false_acceptance,
        });
    }
    for qoi in LEVEL_C_COOLING_QOIS {
        if !finished.iter().any(|cell| cell.qoi == *qoi) {
            known_gaps.push(format!(
                "qoi={qoi} regime={UNREGISTERED_REGIME_LABEL} external_datasets=0"
            ));
        }
    }
    known_gaps.sort();

    let mut adversarial_rows = Vec::with_capacity(adversarial.cases().len());
    for case in adversarial.cases() {
        let assessment = by_case.get(case.id).copied();
        let evidence = match case.evidence {
            AdversarialEvidence::Retained { dataset_id } => format!("retained:{dataset_id}"),
            AdversarialEvidence::Planned { tracking_bead, .. } => {
                format!("NO-DATA:{tracking_bead}")
            }
        };
        adversarial_rows.push(AdversarialCaseRow {
            case_id: case.id,
            regime: case.regime,
            attacked: case.attacked_assumption.slug(),
            evidence,
            verdict: assessment.map_or("NO-DATA", |value| value.verdict().slug()),
            dominant: assessment
                .and_then(|value| value.outcome().dominant())
                .map_or("NO-DATA", DominantUncertainty::slug),
            false_acceptance: assessment
                .filter(|value| value.verdict() != HonestyVerdict::NoData)
                .map(AdversarialAssessment::is_false_acceptance),
            limitation: case.regime_limitation,
        });
    }

    let mut scorecard = VvScorecard {
        corpus_seeded: corpus.is_seeded(),
        corpus_digest: corpus.digest(),
        adversarial_registry_identity: adversarial.identity(),
        dataset_total: corpus.datasets().len(),
        level_counts,
        axis_totals,
        run_record_count: run_records.len(),
        supplied_assessment_count: assessments.len(),
        executed_assessment_count,
        false_acceptance_total,
        adversarial_case_count: adversarial.cases().len(),
        cells: finished,
        known_gaps,
        adversarial_rows,
        adversarial_markdown,
        identity: ContentHash([0; 32]),
    };
    scorecard.identity = scorecard.compute_identity(assessments);
    Ok(scorecard)
}

fn push_len_bytes(bytes: &mut Vec<u8>, len: usize) {
    bytes.extend_from_slice(&u64::try_from(len).unwrap_or(u64::MAX).to_le_bytes());
}

fn push_text_bytes(bytes: &mut Vec<u8>, value: &str) {
    push_len_bytes(bytes, value.len());
    bytes.extend_from_slice(value.as_bytes());
}

fn hex(hash: ContentHash) -> String {
    const DIGITS: &[u8; 16] = b"0123456789abcdef";
    let mut output = String::with_capacity(64);
    for byte in hash.0 {
        output.push(char::from(DIGITS[usize::from(byte >> 4)]));
        output.push(char::from(DIGITS[usize::from(byte & 0x0f)]));
    }
    output
}

impl VvScorecard {
    /// Whether the projected corpus is the seeded workspace corpus.
    #[must_use]
    pub const fn corpus_seeded(&self) -> bool {
        self.corpus_seeded
    }

    /// Projected corpus registry digest.
    #[must_use]
    pub const fn corpus_digest(&self) -> ContentHash {
        self.corpus_digest
    }

    /// Sorted (QoI, regime) outcome cells.
    #[must_use]
    pub fn cells(&self) -> &[ScorecardCell] {
        &self.cells
    }

    /// Sorted loud gap list: cells (and scoped QoIs without any dataset)
    /// with zero external reference coverage.
    #[must_use]
    pub fn known_gaps(&self) -> &[String] {
        &self.known_gaps
    }

    /// Executed (non-`NO-DATA`) assessment count.
    #[must_use]
    pub const fn executed_assessments(&self) -> usize {
        self.executed_assessment_count
    }

    /// Total false acceptances over executed assessments.
    #[must_use]
    pub const fn false_acceptance_total(&self) -> usize {
        self.false_acceptance_total
    }

    /// Canonical scorecard identity over the typed content.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    fn compute_identity(&self, assessments: &[AdversarialAssessment]) -> ContentHash {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&SCORECARD_SCHEMA_VERSION.to_le_bytes());
        bytes.extend_from_slice(&CORPUS_SCHEMA_VERSION.to_le_bytes());
        bytes.extend_from_slice(&ADVERSARIAL_CASE_SCHEMA_VERSION.to_le_bytes());
        bytes.push(u8::from(self.corpus_seeded));
        bytes.extend_from_slice(&self.corpus_digest.0);
        bytes.extend_from_slice(&self.adversarial_registry_identity.0);
        push_len_bytes(&mut bytes, self.dataset_total);
        for count in self.level_counts {
            push_len_bytes(&mut bytes, count);
        }
        for count in self.axis_totals {
            push_len_bytes(&mut bytes, count);
        }
        push_len_bytes(&mut bytes, self.run_record_count);
        push_len_bytes(&mut bytes, self.supplied_assessment_count);
        push_len_bytes(&mut bytes, self.executed_assessment_count);
        push_len_bytes(&mut bytes, self.false_acceptance_total);
        push_len_bytes(&mut bytes, self.cells.len());
        for cell in &self.cells {
            push_text_bytes(&mut bytes, &cell.qoi);
            push_text_bytes(&mut bytes, &cell.regime);
            push_text_bytes(&mut bytes, &cell.regime_statement);
            push_len_bytes(&mut bytes, cell.dataset_ids.len());
            for id in &cell.dataset_ids {
                push_text_bytes(&mut bytes, id);
            }
            for count in cell.axis_counts {
                push_len_bytes(&mut bytes, count);
            }
            push_len_bytes(&mut bytes, cell.external_datasets);
            push_len_bytes(&mut bytes, cell.runs.len());
            for run in &cell.runs {
                push_text_bytes(&mut bytes, &run.dataset_id);
                bytes.extend_from_slice(&run.run_identity.0);
                bytes.extend_from_slice(&run.predicted.to_bits().to_le_bytes());
                bytes.extend_from_slice(&run.reference.to_bits().to_le_bytes());
                match run.reference_uncertainty {
                    ReferenceUncertainty::Bounded { half_width } => {
                        bytes.push(1);
                        bytes.extend_from_slice(&half_width.to_bits().to_le_bytes());
                    }
                    ReferenceUncertainty::Unstated => bytes.push(2),
                }
                bytes.push(match run.envelope {
                    EnvelopeCheck::Pass => 1,
                    EnvelopeCheck::Fail => 2,
                    EnvelopeCheck::Unpinned => 3,
                });
            }
            match cell.false_acceptance {
                FalseAcceptanceCell::NoData => bytes.push(1),
                FalseAcceptanceCell::Counted {
                    executed,
                    false_acceptances,
                } => {
                    bytes.push(2);
                    push_len_bytes(&mut bytes, executed);
                    push_len_bytes(&mut bytes, false_acceptances);
                }
            }
        }
        push_len_bytes(&mut bytes, self.known_gaps.len());
        for gap in &self.known_gaps {
            push_text_bytes(&mut bytes, gap);
        }
        push_len_bytes(&mut bytes, assessments.len());
        let mut assessment_identities = assessments
            .iter()
            .map(|assessment| assessment.identity().0)
            .collect::<Vec<_>>();
        assessment_identities.sort_unstable();
        for identity in assessment_identities {
            bytes.extend_from_slice(&identity);
        }
        hash_domain(SCORECARD_IDENTITY_DOMAIN, &bytes)
    }

    fn render_uncertainty(uncertainty: ReferenceUncertainty) -> String {
        match uncertainty {
            ReferenceUncertainty::Bounded { half_width } => format!("+/-{half_width}"),
            ReferenceUncertainty::Unstated => "uncertainty unstated".to_string(),
        }
    }

    /// Render the publishable Markdown scorecard. Byte-identical for equal
    /// scorecard content.
    #[must_use]
    #[allow(clippy::too_many_lines)] // One linear document assembly.
    pub fn render_markdown(&self) -> String {
        let mut out = String::from("# FrankenSim public V&V scorecard\n\n");
        let _ = writeln!(out, "schema: {SCORECARD_SCHEMA_VERSION}");
        let _ = writeln!(out, "corpus_schema: {CORPUS_SCHEMA_VERSION}");
        let _ = writeln!(
            out,
            "corpus_authority: {}",
            if self.corpus_seeded {
                "seeded"
            } else {
                "caller-built (NOT the public corpus)"
            }
        );
        let _ = writeln!(out, "corpus_digest: {}", hex(self.corpus_digest));
        let _ = writeln!(
            out,
            "adversarial_registry: {}",
            hex(self.adversarial_registry_identity)
        );
        let _ = writeln!(out, "datasets: {}", self.dataset_total);
        let levels = LEVEL_CODES
            .iter()
            .zip(self.level_counts)
            .map(|(code, count)| format!("{code}={count}"))
            .collect::<Vec<_>>()
            .join(" ");
        let _ = writeln!(out, "datasets_by_level: {levels}");
        let axes = EvidenceAxis::ALL
            .iter()
            .map(|axis| format!("{}={}", axis.slug(), self.axis_totals[axis.index()]))
            .collect::<Vec<_>>()
            .join(" ");
        let _ = writeln!(out, "datasets_by_axis: {axes}");
        let _ = writeln!(out, "ledgered_run_records: {}", self.run_record_count);
        let _ = writeln!(
            out,
            "executed_adversarial_challenges: {}/{}",
            self.executed_assessment_count, self.adversarial_case_count
        );
        if self.executed_assessment_count == 0 {
            out.push_str("false_acceptance_total: NO-DATA (0 executed challenges)\n");
        } else {
            let _ = writeln!(
                out,
                "false_acceptance_total: {} of {} executed",
                self.false_acceptance_total, self.executed_assessment_count
            );
        }
        let _ = writeln!(
            out,
            "interval_coverage: NO-DATA ({COVERAGE_NO_DATA_REASON})"
        );
        let external = EXTERNAL_AXES
            .iter()
            .map(|axis| axis.slug())
            .collect::<Vec<_>>()
            .join(", ");
        let _ = writeln!(out, "external_axes: {external}");
        out.push_str(
            "\nThis scorecard is a deterministic projection of the registered validation \
             corpus and the supplied ledgered run results. It grants no authority: a cell \
             with data reports outcome arithmetic; a cell without data reports NO-DATA, \
             never zero. Corpus claim caps remain in force regardless of anything shown \
             here.\n\n",
        );

        out.push_str("## Known gaps\n\n");
        if self.known_gaps.is_empty() {
            out.push_str("- none\n");
        } else {
            for gap in &self.known_gaps {
                let _ = writeln!(out, "- {gap}");
            }
        }
        out.push('\n');

        out.push_str("## Per-QoI/regime cells\n\n");
        out.push_str(
            "Axis order in the `axes` column: numerical-verification / cross-code-agreement / \
             controlled-experimental-validation / blind-predictive-validation / \
             field-monitoring / transferability-across-regimes / independent-reproduction.\n\n",
        );
        out.push_str(
            "| qoi | regime | refs | axes | external | prediction error | envelope | coverage | false acceptance |\n",
        );
        out.push_str("| --- | --- | --- | --- | --- | --- | --- | --- | --- |\n");
        for cell in &self.cells {
            let axes = cell
                .axis_counts
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join("/");
            let error = if cell.runs.is_empty() {
                "NO-DATA".to_string()
            } else {
                let mut min = f64::INFINITY;
                let mut max = f64::NEG_INFINITY;
                let mut max_abs = 0.0_f64;
                for run in &cell.runs {
                    min = min.min(run.signed_error);
                    max = max.max(run.signed_error);
                    max_abs = max_abs.max(run.signed_error.abs());
                }
                format!(
                    "n={} min={min} max={max} max_abs={max_abs}",
                    cell.runs.len()
                )
            };
            let envelope = if cell.runs.is_empty() {
                "NO-DATA".to_string()
            } else {
                format!(
                    "pass={} fail={} unpinned={}",
                    cell.envelope_pass, cell.envelope_fail, cell.envelope_unpinned
                )
            };
            let false_acceptance = match cell.false_acceptance {
                FalseAcceptanceCell::NoData => "NO-DATA".to_string(),
                FalseAcceptanceCell::Counted {
                    executed,
                    false_acceptances,
                } => format!("{false_acceptances} of {executed} executed"),
            };
            let _ = writeln!(
                out,
                "| {} | {} | {} | {} | {} | {} | {} | NO-DATA | {} |",
                cell.qoi,
                cell.regime,
                cell.dataset_ids.len(),
                axes,
                cell.external_datasets,
                error,
                envelope,
                false_acceptance
            );
        }
        out.push('\n');

        if self.run_record_count != 0 {
            out.push_str("### Run detail\n\n");
            out.push_str(
                "| qoi | regime | dataset | predicted | reference | reference uncertainty | signed error | envelope | run identity |\n",
            );
            out.push_str("| --- | --- | --- | --- | --- | --- | --- | --- | --- |\n");
            for cell in &self.cells {
                for run in &cell.runs {
                    let _ = writeln!(
                        out,
                        "| {} | {} | {} | {} | {} | {} | {} | {} | {} |",
                        cell.qoi,
                        cell.regime,
                        run.dataset_id,
                        run.predicted,
                        run.reference,
                        Self::render_uncertainty(run.reference_uncertainty),
                        run.signed_error,
                        run.envelope.slug(),
                        hex(run.run_identity)
                    );
                }
            }
            out.push('\n');
        }

        out.push_str("## Regime limitations\n\n");
        for cell in &self.cells {
            let _ = writeln!(out, "- {}", cell.regime_statement);
        }
        out.push('\n');
        out.push_str(&self.adversarial_markdown);
        out.push('\n');
        let _ = writeln!(out, "identity: {}", hex(self.identity));
        out
    }

    fn push_json_f64(out: &mut String, value: f64) {
        out.push_str("{\"display\":");
        crate::push_json_str(out, &format!("{value}"));
        out.push_str(",\"bits\":");
        crate::push_f64_bits(out, value);
        out.push('}');
    }

    /// Render the canonical machine-readable JSON scorecard for the e16
    /// dashboard. Byte-identical for equal scorecard content.
    #[must_use]
    #[allow(clippy::too_many_lines)] // One linear canonical serialization.
    pub fn render_json(&self) -> String {
        let mut out = String::from("{\"schema\":");
        let _ = write!(out, "{SCORECARD_SCHEMA_VERSION}");
        let _ = write!(out, ",\"corpus_schema\":{CORPUS_SCHEMA_VERSION}");
        out.push_str(",\"corpus_authority\":");
        crate::push_json_str(
            &mut out,
            if self.corpus_seeded {
                "seeded"
            } else {
                "caller-built"
            },
        );
        out.push_str(",\"corpus_digest\":");
        crate::push_json_str(&mut out, &hex(self.corpus_digest));
        out.push_str(",\"adversarial_registry\":");
        crate::push_json_str(&mut out, &hex(self.adversarial_registry_identity));
        let _ = write!(out, ",\"dataset_total\":{}", self.dataset_total);
        out.push_str(",\"datasets_by_level\":{");
        for (index, (code, count)) in LEVEL_CODES.iter().zip(self.level_counts).enumerate() {
            if index != 0 {
                out.push(',');
            }
            crate::push_json_str(&mut out, code);
            let _ = write!(out, ":{count}");
        }
        out.push_str("},\"datasets_by_axis\":{");
        for (index, axis) in EvidenceAxis::ALL.iter().enumerate() {
            if index != 0 {
                out.push(',');
            }
            crate::push_json_str(&mut out, axis.slug());
            let _ = write!(out, ":{}", self.axis_totals[axis.index()]);
        }
        out.push('}');
        let _ = write!(out, ",\"ledgered_run_records\":{}", self.run_record_count);
        let _ = write!(
            out,
            ",\"supplied_assessments\":{}",
            self.supplied_assessment_count
        );
        let _ = write!(
            out,
            ",\"executed_assessments\":{}",
            self.executed_assessment_count
        );
        out.push_str(",\"false_acceptance_total\":");
        if self.executed_assessment_count == 0 {
            out.push_str("{\"status\":\"no-data\"}");
        } else {
            let _ = write!(
                out,
                "{{\"status\":\"counted\",\"count\":{}}}",
                self.false_acceptance_total
            );
        }
        out.push_str(",\"interval_coverage\":{\"status\":\"no-data\",\"reason\":");
        crate::push_json_str(&mut out, COVERAGE_NO_DATA_REASON);
        out.push_str("},\"external_axes\":[");
        for (index, axis) in EXTERNAL_AXES.iter().enumerate() {
            if index != 0 {
                out.push(',');
            }
            crate::push_json_str(&mut out, axis.slug());
        }
        out.push_str("],\"known_gaps\":[");
        for (index, gap) in self.known_gaps.iter().enumerate() {
            if index != 0 {
                out.push(',');
            }
            crate::push_json_str(&mut out, gap);
        }
        out.push_str("],\"cells\":[");
        for (cell_index, cell) in self.cells.iter().enumerate() {
            if cell_index != 0 {
                out.push(',');
            }
            out.push_str("{\"qoi\":");
            crate::push_json_str(&mut out, &cell.qoi);
            out.push_str(",\"regime\":");
            crate::push_json_str(&mut out, &cell.regime);
            out.push_str(",\"regime_statement\":");
            crate::push_json_str(&mut out, &cell.regime_statement);
            out.push_str(",\"datasets\":[");
            for (index, id) in cell.dataset_ids.iter().enumerate() {
                if index != 0 {
                    out.push(',');
                }
                crate::push_json_str(&mut out, id);
            }
            out.push_str("],\"axis_counts\":{");
            for (index, axis) in EvidenceAxis::ALL.iter().enumerate() {
                if index != 0 {
                    out.push(',');
                }
                crate::push_json_str(&mut out, axis.slug());
                let _ = write!(out, ":{}", cell.axis_counts[axis.index()]);
            }
            out.push('}');
            let _ = write!(out, ",\"external_datasets\":{}", cell.external_datasets);
            out.push_str(",\"runs\":[");
            for (index, run) in cell.runs.iter().enumerate() {
                if index != 0 {
                    out.push(',');
                }
                out.push_str("{\"dataset\":");
                crate::push_json_str(&mut out, &run.dataset_id);
                out.push_str(",\"run_identity\":");
                crate::push_json_str(&mut out, &hex(run.run_identity));
                out.push_str(",\"predicted\":");
                Self::push_json_f64(&mut out, run.predicted);
                out.push_str(",\"reference\":");
                Self::push_json_f64(&mut out, run.reference);
                out.push_str(",\"reference_uncertainty\":");
                match run.reference_uncertainty {
                    ReferenceUncertainty::Bounded { half_width } => {
                        out.push_str("{\"status\":\"bounded\",\"half_width\":");
                        Self::push_json_f64(&mut out, half_width);
                        out.push('}');
                    }
                    ReferenceUncertainty::Unstated => {
                        out.push_str("{\"status\":\"unstated\"}");
                    }
                }
                out.push_str(",\"signed_error\":");
                Self::push_json_f64(&mut out, run.signed_error);
                out.push_str(",\"envelope\":");
                crate::push_json_str(&mut out, run.envelope.slug());
                out.push('}');
            }
            out.push_str("],\"envelope\":");
            if cell.runs.is_empty() {
                out.push_str("{\"status\":\"no-data\"}");
            } else {
                let _ = write!(
                    out,
                    "{{\"status\":\"counted\",\"pass\":{},\"fail\":{},\"unpinned\":{}}}",
                    cell.envelope_pass, cell.envelope_fail, cell.envelope_unpinned
                );
            }
            out.push_str(",\"interval_coverage\":{\"status\":\"no-data\",\"reason\":");
            crate::push_json_str(&mut out, COVERAGE_NO_DATA_REASON);
            out.push_str("},\"false_acceptance\":");
            match cell.false_acceptance {
                FalseAcceptanceCell::NoData => out.push_str("{\"status\":\"no-data\"}"),
                FalseAcceptanceCell::Counted {
                    executed,
                    false_acceptances,
                } => {
                    let _ = write!(
                        out,
                        "{{\"status\":\"counted\",\"executed\":{executed},\"count\":{false_acceptances}}}"
                    );
                }
            }
            out.push('}');
        }
        out.push_str("],\"adversarial\":{\"schema\":");
        let _ = write!(out, "{ADVERSARIAL_CASE_SCHEMA_VERSION}");
        out.push_str(",\"registry\":");
        crate::push_json_str(&mut out, &hex(self.adversarial_registry_identity));
        out.push_str(",\"false_acceptance_count\":");
        if self.executed_assessment_count == 0 {
            out.push_str("{\"status\":\"no-data\"}");
        } else {
            let _ = write!(
                out,
                "{{\"status\":\"counted\",\"count\":{}}}",
                self.false_acceptance_total
            );
        }
        out.push_str(",\"cases\":[");
        for (index, row) in self.adversarial_rows.iter().enumerate() {
            if index != 0 {
                out.push(',');
            }
            out.push_str("{\"id\":");
            crate::push_json_str(&mut out, row.case_id);
            out.push_str(",\"regime\":");
            crate::push_json_str(&mut out, row.regime);
            out.push_str(",\"attacked_assumption\":");
            crate::push_json_str(&mut out, row.attacked);
            out.push_str(",\"evidence\":");
            crate::push_json_str(&mut out, &row.evidence);
            out.push_str(",\"verdict\":");
            crate::push_json_str(&mut out, row.verdict);
            out.push_str(",\"dominant\":");
            crate::push_json_str(&mut out, row.dominant);
            out.push_str(",\"false_acceptance\":");
            match row.false_acceptance {
                Some(true) => out.push_str("true"),
                Some(false) => out.push_str("false"),
                None => out.push_str("null"),
            }
            out.push_str(",\"limitation\":");
            crate::push_json_str(&mut out, row.limitation);
            out.push('}');
        }
        out.push_str("]},\"identity\":");
        crate::push_json_str(&mut out, &hex(self.identity));
        out.push('}');
        out
    }
}

/// Deterministic content identity for one run record, usable by callers
/// that persist run comparisons before scorecard assembly.
#[must_use]
pub fn run_record_identity(record: &ScorecardRunRecord) -> ContentHash {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(&SCORECARD_SCHEMA_VERSION.to_le_bytes());
    push_text_bytes(&mut bytes, record.dataset_id());
    push_text_bytes(&mut bytes, record.metric());
    bytes.extend_from_slice(&record.predicted().to_bits().to_le_bytes());
    bytes.extend_from_slice(&record.reference().to_bits().to_le_bytes());
    match record.reference_uncertainty() {
        ReferenceUncertainty::Bounded { half_width } => {
            bytes.push(1);
            bytes.extend_from_slice(&half_width.to_bits().to_le_bytes());
        }
        ReferenceUncertainty::Unstated => bytes.push(2),
    }
    bytes.extend_from_slice(&record.run_identity().0);
    hash_domain(RUN_RECORD_IDENTITY_DOMAIN, &bytes)
}
