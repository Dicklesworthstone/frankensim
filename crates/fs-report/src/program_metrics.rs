//! Program metrics dashboard: outcome metrics as one deterministic artifact.
//!
//! The dashboard is a pure projection of caller-supplied metric rows and a
//! recorded generation history. It reports the program's OUTCOME metrics —
//! how much of the claim surface is anchored to external references, how much
//! of the adversarial suite has actually been executed, how far registered
//! capabilities have matured — and renders byte-identical Markdown and
//! canonical JSON so its history is diffable and its claims are auditable.
//!
//! Honesty rules:
//! - `NO-DATA` and a measured zero are DIFFERENT states and are never
//!   conflated. `NO-DATA` means no measurement machinery exists, so any
//!   number would be invented. A measured zero means the population is
//!   enumerable and the answer is genuinely none — rendering that as
//!   `NO-DATA` would hide a known-bad result, which is the opposite failure.
//! - A ratio's denominator is [`NonZeroU64`]: an empty population can never
//!   masquerade as `0%`.
//! - Every row must state what it does NOT capture. A row without an honesty
//!   statement is refused, because an unqualified number is precisely the
//!   failure mode this dashboard exists to prevent.
//! - Every measured row must cite at least one data source. A measurement
//!   nobody can trace back to an artifact is not evidence.
//! - Trend is computed against a RECORDED prior generation, never against the
//!   dashboard's own current render, so generation is idempotent and the
//!   artifact cannot become its own trend basis.
//! - Supplying rows or history grants no authority. The dashboard reports
//!   movement and coverage; it never upgrades a claim, and a metric moving in
//!   the better direction is not evidence that any underlying claim is true.
//! - Kernel throughput, crate counts, and test-file counts are deliberately
//!   absent from the metric set: they are diagnostics, not outcomes. See
//!   [`DEMOTED_DIAGNOSTICS`].

use core::fmt::{self, Write as _};
use core::num::NonZeroU64;
use std::collections::{BTreeMap, BTreeSet};

use fs_blake3::{ContentHash, hash_domain};
use fs_vvreg::adversarial::AdversarialRegistry;
use fs_vvreg::corpus::CorpusRegistry;
use fs_vvreg::portfolio::EvidenceAxis;
use fs_vvreg::scorecard::{EXTERNAL_AXES, VvScorecard};

/// Schema version of the rendered dashboard artifacts.
pub const DASHBOARD_SCHEMA_VERSION: u32 = 1;

/// Schema version of one recorded history generation.
pub const HISTORY_SCHEMA_VERSION: u32 = 1;

/// Upper bound on metric rows in one dashboard.
pub const MAX_METRIC_ROWS: usize = 512;

/// Upper bound on recorded generations admitted from a history file.
pub const MAX_HISTORY_GENERATIONS: usize = 4_096;

/// Upper bound on any caller-supplied text field, in bytes.
pub const MAX_METRIC_TEXT_BYTES: usize = 4_096;

const DASHBOARD_IDENTITY_DOMAIN: &str = "org.frankensim.fs-report.program-metrics.v1";
const SOURCE_IDENTITY_DOMAIN: &str = "org.frankensim.fs-report.program-metrics-sources.v1";

/// Measurements deliberately excluded from the metric set, with the reason.
///
/// The program optimizes what it measures, so it must measure outcomes. These
/// are legitimate signals but they are diagnostics: they move without the
/// program getting better at predicting reality, and each already has its own
/// lane. Listing them in the artifact explains an absence instead of hiding it.
pub const DEMOTED_DIAGNOSTICS: [(&str, &str); 4] = [
    (
        "kernel throughput",
        "a performance diagnostic owned by the roofline lane; a faster kernel is not a better \
         prediction",
    ),
    (
        "crate count",
        "inventory, not capability; the capability maturity registry is the outcome measure",
    ),
    (
        "integration-test file count",
        "inventory, not proof; check-docs already pins it and a test file is not an outcome",
    ),
    (
        "open issue counts",
        "the beads store churns on every unrelated issue edit, which would make this checked \
         artifact stale for every agent in the repository",
    ),
];

/// Fail-closed dashboard build refusal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DashboardError {
    /// A caller-supplied collection exceeded its admitted bound.
    ResourceLimit {
        /// The admitted bound.
        limit: usize,
        /// The observed size.
        observed: usize,
    },
    /// A metric row field failed validation.
    InvalidMetricField {
        /// Stable machine name of the offending field.
        field: &'static str,
        /// Why the value was refused.
        reason: &'static str,
    },
    /// Two rows declared the same metric identifier.
    DuplicateMetric {
        /// The repeated identifier.
        metric: String,
    },
    /// A measured row cited no data source.
    MissingSourceCitation {
        /// The offending identifier.
        metric: String,
    },
    /// Recorded generations were not strictly increasing.
    NonMonotonicHistory {
        /// The generation number that did not advance.
        generation: u64,
    },
    /// A required source registry was empty, so no metric could be projected.
    ///
    /// An empty registry is a READ FAILURE, not a program with nothing in it.
    /// Rendering it as zeros would invent the most flattering possible answer.
    EmptySourceRegistry {
        /// Stable machine name of the registry.
        source: &'static str,
    },
    /// A source registry carried a value the projection cannot interpret.
    UninterpretableSource {
        /// Stable machine name of the registry.
        source: &'static str,
        /// The offending value.
        value: String,
    },
    /// The projection was handed a corpus that is not the seeded public one.
    UnseededCorpus,
}

impl fmt::Display for DashboardError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ResourceLimit { limit, observed } => write!(
                formatter,
                "dashboard input exceeds the admitted bound: limit {limit}, observed {observed}"
            ),
            Self::InvalidMetricField { field, reason } => {
                write!(formatter, "metric field `{field}` is invalid: {reason}")
            }
            Self::DuplicateMetric { metric } => {
                write!(formatter, "metric `{metric}` is declared more than once")
            }
            Self::MissingSourceCitation { metric } => write!(
                formatter,
                "measured metric `{metric}` cites no data source; an untraceable measurement is \
                 not evidence"
            ),
            Self::NonMonotonicHistory { generation } => write!(
                formatter,
                "recorded generation {generation} does not advance the history; generations must \
                 strictly increase"
            ),
            Self::EmptySourceRegistry { source } => write!(
                formatter,
                "source registry `{source}` is empty; an unreadable registry is a read failure, \
                 not a program with nothing in it, and is never rendered as zeros"
            ),
            Self::UninterpretableSource { source, value } => write!(
                formatter,
                "source registry `{source}` carries uninterpretable value `{value}`"
            ),
            Self::UnseededCorpus => write!(
                formatter,
                "the supplied corpus is caller-built; a synthetic registry may not masquerade as \
                 the public program state"
            ),
        }
    }
}

impl std::error::Error for DashboardError {}

/// A known measurement.
///
/// Both shapes carry exact integers. The dashboard performs no floating-point
/// arithmetic: percentages are rendered from integer basis points, so a render
/// cannot drift with rounding mode or build mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MetricObservation {
    /// A part of an enumerable population.
    Ratio {
        /// Members of the population satisfying the metric.
        numerator: u64,
        /// Total enumerable population. Nonzero by construction.
        denominator: NonZeroU64,
    },
    /// An absolute count with no meaningful denominator.
    Count {
        /// The counted value.
        value: u64,
    },
}

impl MetricObservation {
    /// A ratio over a nonzero population.
    #[must_use]
    pub const fn ratio(numerator: u64, denominator: NonZeroU64) -> Self {
        Self::Ratio {
            numerator,
            denominator,
        }
    }

    /// An absolute count.
    #[must_use]
    pub const fn count(value: u64) -> Self {
        Self::Count { value }
    }

    /// Hundredths of a percent, or `None` for a count.
    ///
    /// Integer arithmetic throughout: `numerator * 10_000 / denominator` in
    /// `u128`, so no float ever reaches a render.
    #[must_use]
    pub fn basis_points(self) -> Option<u64> {
        match self {
            Self::Ratio {
                numerator,
                denominator,
            } => {
                let scaled = u128::from(numerator) * 10_000 / u128::from(denominator.get());
                Some(u64::try_from(scaled).unwrap_or(u64::MAX))
            }
            Self::Count { .. } => None,
        }
    }

    /// Ordering against another observation of the SAME shape.
    ///
    /// Comparing a ratio with a count returns `None`: the metric changed
    /// shape, so no movement is claimed rather than an invented one. Ratios
    /// are compared by exact `u128` cross-multiplication, never by division.
    #[must_use]
    pub fn compare(self, other: Self) -> Option<core::cmp::Ordering> {
        match (self, other) {
            (
                Self::Ratio {
                    numerator: left_numerator,
                    denominator: left_denominator,
                },
                Self::Ratio {
                    numerator: right_numerator,
                    denominator: right_denominator,
                },
            ) => {
                let left = u128::from(left_numerator) * u128::from(right_denominator.get());
                let right = u128::from(right_numerator) * u128::from(left_denominator.get());
                Some(left.cmp(&right))
            }
            (Self::Count { value: left }, Self::Count { value: right }) => Some(left.cmp(&right)),
            _ => None,
        }
    }

    /// Human-readable value cell.
    #[must_use]
    pub fn render(self) -> String {
        match self {
            Self::Ratio {
                numerator,
                denominator,
            } => {
                let basis = self.basis_points().unwrap_or(0);
                format!(
                    "{numerator} of {denominator} ({}.{:02}%)",
                    basis / 100,
                    basis % 100
                )
            }
            Self::Count { value } => format!("{value}"),
        }
    }
}

/// One metric's value, distinguishing "unmeasured" from "measured as zero".
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MetricCell {
    /// No measurement machinery exists, so no number is invented.
    NoData {
        /// Why the metric cannot be measured yet.
        reason: String,
        /// Tracked work that would make it live, when one is identified.
        unblocked_by: Option<String>,
    },
    /// A real measurement, including a genuine zero.
    Measured(MetricObservation),
}

impl MetricCell {
    /// Whether this cell carries a measurement.
    #[must_use]
    pub const fn is_measured(&self) -> bool {
        matches!(self, Self::Measured(_))
    }

    /// The observation, or `None` when the metric is `NO-DATA`.
    #[must_use]
    pub const fn observation(&self) -> Option<MetricObservation> {
        match self {
            Self::Measured(observation) => Some(*observation),
            Self::NoData { .. } => None,
        }
    }

    /// Stable machine status token.
    #[must_use]
    pub const fn status(&self) -> &'static str {
        match self {
            Self::NoData { .. } => "no-data",
            Self::Measured(MetricObservation::Ratio { .. }) => "ratio",
            Self::Measured(MetricObservation::Count { .. }) => "count",
        }
    }
}

/// Which direction of movement counts as improvement.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MetricDirection {
    /// A larger value is better.
    HigherIsBetter,
    /// A smaller value is better.
    LowerIsBetter,
    /// Movement is reported without a better/worse judgement.
    Neutral,
}

impl MetricDirection {
    /// Stable machine name.
    #[must_use]
    pub const fn slug(self) -> &'static str {
        match self {
            Self::HigherIsBetter => "higher-is-better",
            Self::LowerIsBetter => "lower-is-better",
            Self::Neutral => "neutral",
        }
    }
}

/// Metric grouping, in render order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum MetricFamily {
    /// Does the program predict reality and turn inputs into decisions?
    Outcome,
    /// How anchored is the claim surface to external evidence?
    Portfolio,
    /// Maturity, honesty debt, and promotion discipline.
    Governance,
}

impl MetricFamily {
    /// Canonical family order.
    pub const ALL: [Self; 3] = [Self::Outcome, Self::Portfolio, Self::Governance];

    /// Stable machine name.
    #[must_use]
    pub const fn slug(self) -> &'static str {
        match self {
            Self::Outcome => "outcome",
            Self::Portfolio => "portfolio",
            Self::Governance => "governance",
        }
    }

    /// Human-readable section title.
    #[must_use]
    pub const fn title(self) -> &'static str {
        match self {
            Self::Outcome => "Outcome metrics",
            Self::Portfolio => "Evidence portfolio metrics",
            Self::Governance => "Governance metrics",
        }
    }

    /// Stable zero-based position in [`Self::ALL`].
    #[must_use]
    pub const fn index(self) -> usize {
        match self {
            Self::Outcome => 0,
            Self::Portfolio => 1,
            Self::Governance => 2,
        }
    }
}

/// Movement of one metric against the most recent recorded generation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Trend {
    /// No recorded generation carries this metric.
    NoPrior,
    /// The recorded value is identical.
    Unchanged,
    /// Moved in this metric's better direction.
    Improved,
    /// Moved in this metric's worse direction.
    Worsened,
    /// Moved, but the metric declares no better/worse direction.
    Changed,
    /// The prior generation had no data and this one measures it.
    NewlyMeasured,
    /// The prior generation measured it and this one does not.
    LostData,
    /// The observation changed shape, so no movement is claimed.
    NotComparable,
}

impl Trend {
    /// Stable machine name.
    #[must_use]
    pub const fn slug(self) -> &'static str {
        match self {
            Self::NoPrior => "no-prior",
            Self::Unchanged => "unchanged",
            Self::Improved => "improved",
            Self::Worsened => "worsened",
            Self::Changed => "changed",
            Self::NewlyMeasured => "newly-measured",
            Self::LostData => "lost-data",
            Self::NotComparable => "not-comparable",
        }
    }

    /// Human-readable trend cell.
    #[must_use]
    pub const fn render(self) -> &'static str {
        match self {
            Self::NoPrior => "no prior generation",
            Self::Unchanged => "unchanged",
            Self::Improved => "improved",
            Self::Worsened => "WORSENED",
            Self::Changed => "changed",
            Self::NewlyMeasured => "newly measured",
            Self::LostData => "LOST DATA",
            Self::NotComparable => "not comparable",
        }
    }
}

fn validate_text(field: &'static str, value: &str) -> Result<(), DashboardError> {
    if value.trim().is_empty() {
        return Err(DashboardError::InvalidMetricField {
            field,
            reason: "must not be blank",
        });
    }
    if value.len() > MAX_METRIC_TEXT_BYTES {
        return Err(DashboardError::InvalidMetricField {
            field,
            reason: "exceeds the admitted text bound",
        });
    }
    if value.chars().any(char::is_control) {
        return Err(DashboardError::InvalidMetricField {
            field,
            reason: "must not contain control characters",
        });
    }
    Ok(())
}

/// One dashboard row: a metric, its value, its provenance, and its limits.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MetricRow {
    id: String,
    title: String,
    family: MetricFamily,
    direction: MetricDirection,
    cell: MetricCell,
    sources: Vec<String>,
    honesty: String,
}

impl MetricRow {
    /// Admit one metric row.
    ///
    /// Refusals, all total: a blank or control-bearing identifier, title,
    /// honesty statement, source reference, or `NO-DATA` reason; a measured
    /// row citing no source.
    ///
    /// The honesty statement is mandatory for every row, measured or not. A
    /// number without a stated boundary is the failure this artifact exists
    /// to prevent, so the type system refuses to carry one.
    ///
    /// # Errors
    ///
    /// [`DashboardError::InvalidMetricField`] for any failed field check, and
    /// [`DashboardError::MissingSourceCitation`] for an uncited measurement.
    pub fn try_new(
        id: &str,
        title: &str,
        family: MetricFamily,
        direction: MetricDirection,
        cell: MetricCell,
        sources: &[&str],
        honesty: &str,
    ) -> Result<Self, DashboardError> {
        validate_text("id", id)?;
        validate_text("title", title)?;
        validate_text("honesty", honesty)?;
        if let MetricCell::NoData {
            reason,
            unblocked_by,
        } = &cell
        {
            validate_text("no_data_reason", reason)?;
            if let Some(bead) = unblocked_by {
                validate_text("unblocked_by", bead)?;
            }
        }
        for source in sources {
            validate_text("source", source)?;
        }
        if cell.is_measured() && sources.is_empty() {
            return Err(DashboardError::MissingSourceCitation {
                metric: id.to_string(),
            });
        }
        let mut ordered: Vec<String> = sources
            .iter()
            .map(|source| (*source).to_string())
            .collect::<BTreeSet<String>>()
            .into_iter()
            .collect();
        ordered.shrink_to_fit();
        Ok(Self {
            id: id.to_string(),
            title: title.to_string(),
            family,
            direction,
            cell,
            sources: ordered,
            honesty: honesty.to_string(),
        })
    }

    /// Stable machine identifier.
    #[must_use]
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Human-readable metric name.
    #[must_use]
    pub fn title(&self) -> &str {
        &self.title
    }

    /// Render grouping.
    #[must_use]
    pub const fn family(&self) -> MetricFamily {
        self.family
    }

    /// Which direction counts as improvement.
    #[must_use]
    pub const fn direction(&self) -> MetricDirection {
        self.direction
    }

    /// The measured value, or the typed `NO-DATA` state.
    #[must_use]
    pub const fn cell(&self) -> &MetricCell {
        &self.cell
    }

    /// Sorted, deduplicated data-source references.
    #[must_use]
    pub fn sources(&self) -> &[String] {
        &self.sources
    }

    /// What this metric does NOT capture.
    #[must_use]
    pub fn honesty(&self) -> &str {
        &self.honesty
    }
}

/// One recorded dashboard generation, used only as a trend basis.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HistoryGeneration {
    generation: u64,
    label: String,
    source_identity: ContentHash,
    metrics: BTreeMap<String, MetricCell>,
}

impl HistoryGeneration {
    /// Admit one recorded generation.
    ///
    /// # Errors
    ///
    /// [`DashboardError::InvalidMetricField`] for a blank or oversized label,
    /// and [`DashboardError::ResourceLimit`] for an oversized metric set.
    pub fn try_new(
        generation: u64,
        label: &str,
        source_identity: ContentHash,
        metrics: BTreeMap<String, MetricCell>,
    ) -> Result<Self, DashboardError> {
        validate_text("history_label", label)?;
        if metrics.len() > MAX_METRIC_ROWS {
            return Err(DashboardError::ResourceLimit {
                limit: MAX_METRIC_ROWS,
                observed: metrics.len(),
            });
        }
        Ok(Self {
            generation,
            label: label.to_string(),
            source_identity,
            metrics,
        })
    }

    /// Monotonic generation number.
    #[must_use]
    pub const fn generation(&self) -> u64 {
        self.generation
    }

    /// Caller-supplied label for the recorded state.
    #[must_use]
    pub fn label(&self) -> &str {
        &self.label
    }

    /// Identity of the source values that produced this generation.
    #[must_use]
    pub const fn source_identity(&self) -> ContentHash {
        self.source_identity
    }

    /// Recorded value for one metric, if the generation carried it.
    #[must_use]
    pub fn metric(&self, id: &str) -> Option<&MetricCell> {
        self.metrics.get(id)
    }
}

/// Recorded generations in strictly increasing order.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct MetricHistory {
    generations: Vec<HistoryGeneration>,
}

impl MetricHistory {
    /// A history with no recorded generation. Every trend is `NoPrior`.
    #[must_use]
    pub const fn empty() -> Self {
        Self {
            generations: Vec::new(),
        }
    }

    /// Admit recorded generations.
    ///
    /// # Errors
    ///
    /// [`DashboardError::ResourceLimit`] beyond [`MAX_HISTORY_GENERATIONS`],
    /// and [`DashboardError::NonMonotonicHistory`] when a generation number
    /// fails to advance — a history that can reorder is not a trend basis.
    pub fn try_new(generations: Vec<HistoryGeneration>) -> Result<Self, DashboardError> {
        if generations.len() > MAX_HISTORY_GENERATIONS {
            return Err(DashboardError::ResourceLimit {
                limit: MAX_HISTORY_GENERATIONS,
                observed: generations.len(),
            });
        }
        for pair in generations.windows(2) {
            if pair[1].generation <= pair[0].generation {
                return Err(DashboardError::NonMonotonicHistory {
                    generation: pair[1].generation,
                });
            }
        }
        Ok(Self { generations })
    }

    /// The most recently recorded generation.
    #[must_use]
    pub fn latest(&self) -> Option<&HistoryGeneration> {
        self.generations.last()
    }

    /// Number of recorded generations.
    #[must_use]
    pub fn len(&self) -> usize {
        self.generations.len()
    }

    /// Whether no generation has been recorded.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.generations.is_empty()
    }

    /// The generation number a new recording should carry.
    #[must_use]
    pub fn next_generation(&self) -> u64 {
        self.latest()
            .map_or(1, |generation| generation.generation.saturating_add(1))
    }
}

fn trend_for(
    current: &MetricCell,
    prior: Option<&MetricCell>,
    direction: MetricDirection,
) -> Trend {
    let Some(prior) = prior else {
        return Trend::NoPrior;
    };
    match (prior.observation(), current.observation()) {
        (None, None) => Trend::Unchanged,
        (None, Some(_)) => Trend::NewlyMeasured,
        (Some(_), None) => Trend::LostData,
        (Some(before), Some(after)) => match after.compare(before) {
            None => Trend::NotComparable,
            Some(core::cmp::Ordering::Equal) => Trend::Unchanged,
            Some(ordering) => match direction {
                MetricDirection::Neutral => Trend::Changed,
                MetricDirection::HigherIsBetter => {
                    if ordering == core::cmp::Ordering::Greater {
                        Trend::Improved
                    } else {
                        Trend::Worsened
                    }
                }
                MetricDirection::LowerIsBetter => {
                    if ordering == core::cmp::Ordering::Less {
                        Trend::Improved
                    } else {
                        Trend::Worsened
                    }
                }
            },
        },
    }
}

/// A rendered row: an admitted metric plus its computed trend.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DashboardRow {
    row: MetricRow,
    trend: Trend,
}

impl DashboardRow {
    /// The admitted metric.
    #[must_use]
    pub const fn metric(&self) -> &MetricRow {
        &self.row
    }

    /// Movement against the most recent recorded generation.
    #[must_use]
    pub const fn trend(&self) -> Trend {
        self.trend
    }
}

/// The deterministic program metrics dashboard.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProgramDashboard {
    rows: Vec<DashboardRow>,
    measured: usize,
    no_data: usize,
    prior_generation: Option<u64>,
    prior_label: String,
    sources_identical_to_prior: bool,
    source_identity: ContentHash,
    identity: ContentHash,
}

/// Project metric rows and recorded history into the dashboard.
///
/// The projection is pure and total in its refusals: no partial dashboard is
/// produced. Rows are ordered by `(family, id)` so caller order cannot move
/// the artifact.
///
/// # Errors
///
/// [`DashboardError::ResourceLimit`] beyond [`MAX_METRIC_ROWS`] and
/// [`DashboardError::DuplicateMetric`] for a repeated identifier.
pub fn build_dashboard(
    rows: &[MetricRow],
    history: &MetricHistory,
) -> Result<ProgramDashboard, DashboardError> {
    if rows.len() > MAX_METRIC_ROWS {
        return Err(DashboardError::ResourceLimit {
            limit: MAX_METRIC_ROWS,
            observed: rows.len(),
        });
    }
    let mut ordered: BTreeMap<(usize, String), &MetricRow> = BTreeMap::new();
    for row in rows {
        let key = (row.family.index(), row.id.clone());
        if ordered.insert(key, row).is_some() {
            return Err(DashboardError::DuplicateMetric {
                metric: row.id.clone(),
            });
        }
    }

    let latest = history.latest();
    let mut built = Vec::with_capacity(ordered.len());
    let mut measured = 0_usize;
    let mut no_data = 0_usize;
    for row in ordered.into_values() {
        let prior = latest.and_then(|generation| generation.metric(&row.id));
        let trend = trend_for(&row.cell, prior, row.direction);
        if row.cell.is_measured() {
            measured += 1;
        } else {
            no_data += 1;
        }
        built.push(DashboardRow {
            row: row.clone(),
            trend,
        });
    }

    let source_identity = compute_source_identity(&built);
    let sources_identical_to_prior =
        latest.is_some_and(|generation| generation.source_identity == source_identity);
    let mut dashboard = ProgramDashboard {
        rows: built,
        measured,
        no_data,
        prior_generation: latest.map(HistoryGeneration::generation),
        prior_label: latest.map_or_else(String::new, |generation| generation.label.clone()),
        sources_identical_to_prior,
        source_identity,
        identity: ContentHash([0; 32]),
    };
    dashboard.identity = dashboard.compute_identity();
    Ok(dashboard)
}

fn push_len_bytes(bytes: &mut Vec<u8>, len: usize) {
    bytes.extend_from_slice(&u64::try_from(len).unwrap_or(u64::MAX).to_le_bytes());
}

fn push_text_bytes(bytes: &mut Vec<u8>, value: &str) {
    push_len_bytes(bytes, value.len());
    bytes.extend_from_slice(value.as_bytes());
}

fn push_cell_bytes(bytes: &mut Vec<u8>, cell: &MetricCell) {
    match cell {
        MetricCell::NoData {
            reason,
            unblocked_by,
        } => {
            bytes.push(1);
            push_text_bytes(bytes, reason);
            match unblocked_by {
                Some(bead) => {
                    bytes.push(1);
                    push_text_bytes(bytes, bead);
                }
                None => bytes.push(0),
            }
        }
        MetricCell::Measured(MetricObservation::Ratio {
            numerator,
            denominator,
        }) => {
            bytes.push(2);
            bytes.extend_from_slice(&numerator.to_le_bytes());
            bytes.extend_from_slice(&denominator.get().to_le_bytes());
        }
        MetricCell::Measured(MetricObservation::Count { value }) => {
            bytes.push(3);
            bytes.extend_from_slice(&value.to_le_bytes());
        }
    }
}

/// Identity of the metric VALUES alone, excluding trend and history.
///
/// Two generations sharing this identity were computed from identical source
/// values, which is what makes "unchanged" a provable statement rather than a
/// coincidence of equal renders.
fn compute_source_identity(rows: &[DashboardRow]) -> ContentHash {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(&DASHBOARD_SCHEMA_VERSION.to_le_bytes());
    push_len_bytes(&mut bytes, rows.len());
    for entry in rows {
        push_text_bytes(&mut bytes, &entry.row.id);
        push_cell_bytes(&mut bytes, &entry.row.cell);
    }
    hash_domain(SOURCE_IDENTITY_DOMAIN, &bytes)
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

fn escape_json_into(out: &mut String, value: &str) {
    for character in value.chars() {
        match character {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            control if (control as u32) < 0x20 => {
                let code = control as u32;
                let _ = write!(out, "\\u{code:04x}");
            }
            other => out.push(other),
        }
    }
}

fn push_json_str(out: &mut String, value: &str) {
    out.push('"');
    escape_json_into(out, value);
    out.push('"');
}

fn push_cell_json(out: &mut String, cell: &MetricCell) {
    match cell {
        MetricCell::NoData {
            reason,
            unblocked_by,
        } => {
            out.push_str("{\"status\":\"no-data\",\"reason\":");
            push_json_str(out, reason);
            out.push_str(",\"unblocked_by\":");
            match unblocked_by {
                Some(bead) => push_json_str(out, bead),
                None => out.push_str("null"),
            }
            out.push('}');
        }
        MetricCell::Measured(observation) => match observation {
            MetricObservation::Ratio {
                numerator,
                denominator,
            } => {
                let basis = observation.basis_points().unwrap_or(0);
                let _ = write!(
                    out,
                    "{{\"status\":\"ratio\",\"numerator\":{numerator},\"denominator\":{},\
                     \"basis_points\":{basis}}}",
                    denominator.get()
                );
            }
            MetricObservation::Count { value } => {
                let _ = write!(out, "{{\"status\":\"count\",\"value\":{value}}}");
            }
        },
    }
}

impl ProgramDashboard {
    /// Rows in render order: `(family, id)`.
    #[must_use]
    pub fn rows(&self) -> &[DashboardRow] {
        &self.rows
    }

    /// Rows carrying a measurement.
    #[must_use]
    pub const fn measured_count(&self) -> usize {
        self.measured
    }

    /// Rows rendering loud `NO-DATA`.
    #[must_use]
    pub const fn no_data_count(&self) -> usize {
        self.no_data
    }

    /// Generation number of the trend basis, if any was recorded.
    #[must_use]
    pub const fn prior_generation(&self) -> Option<u64> {
        self.prior_generation
    }

    /// Whether the trend basis was computed from identical source values.
    #[must_use]
    pub const fn sources_identical_to_prior(&self) -> bool {
        self.sources_identical_to_prior
    }

    /// Identity of the metric values alone.
    #[must_use]
    pub const fn source_identity(&self) -> ContentHash {
        self.source_identity
    }

    /// Identity of the whole projection, including computed trends.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    fn compute_identity(&self) -> ContentHash {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&DASHBOARD_SCHEMA_VERSION.to_le_bytes());
        bytes.extend_from_slice(&self.source_identity.0);
        bytes.push(u8::from(self.sources_identical_to_prior));
        match self.prior_generation {
            Some(generation) => {
                bytes.push(1);
                bytes.extend_from_slice(&generation.to_le_bytes());
            }
            None => bytes.push(0),
        }
        push_text_bytes(&mut bytes, &self.prior_label);
        push_len_bytes(&mut bytes, self.rows.len());
        for entry in &self.rows {
            push_text_bytes(&mut bytes, &entry.row.id);
            push_text_bytes(&mut bytes, &entry.row.title);
            push_text_bytes(&mut bytes, entry.row.family.slug());
            push_text_bytes(&mut bytes, entry.row.direction.slug());
            push_cell_bytes(&mut bytes, &entry.row.cell);
            push_len_bytes(&mut bytes, entry.row.sources.len());
            for source in &entry.row.sources {
                push_text_bytes(&mut bytes, source);
            }
            push_text_bytes(&mut bytes, &entry.row.honesty);
            push_text_bytes(&mut bytes, entry.trend.slug());
        }
        hash_domain(DASHBOARD_IDENTITY_DOMAIN, &bytes)
    }

    /// One JSONL line recording this dashboard as a future trend basis.
    ///
    /// The recorded line carries only the metric VALUES and the source
    /// identity — never the trends, which are derived. Recording is a
    /// deliberate, separate act so that generating the dashboard can never
    /// mutate its own trend basis.
    #[must_use]
    pub fn render_history_line(&self, generation: u64, label: &str) -> String {
        let mut out = String::new();
        let _ = write!(
            out,
            "{{\"schema\":{HISTORY_SCHEMA_VERSION},\"generation\":{generation},\"label\":"
        );
        push_json_str(&mut out, label);
        out.push_str(",\"source_identity\":");
        push_json_str(&mut out, &hex(self.source_identity));
        out.push_str(",\"metrics\":{");
        for (index, entry) in self.rows.iter().enumerate() {
            if index != 0 {
                out.push(',');
            }
            push_json_str(&mut out, &entry.row.id);
            out.push(':');
            push_cell_json(&mut out, &entry.row.cell);
        }
        out.push_str("}}\n");
        out
    }

    /// Deterministic human-readable dashboard.
    #[allow(clippy::too_many_lines)] // One linear document; splitting it would hide the render order.
    #[must_use]
    pub fn render_markdown(&self) -> String {
        let mut out = String::new();
        out.push_str("# FrankenSim program metrics dashboard\n\n");
        let _ = writeln!(out, "schema: {DASHBOARD_SCHEMA_VERSION}");
        let _ = writeln!(out, "metrics: {}", self.rows.len());
        let _ = writeln!(out, "measured: {}", self.measured);
        let _ = writeln!(out, "no_data: {}", self.no_data);
        match self.prior_generation {
            Some(generation) => {
                let _ = writeln!(
                    out,
                    "trend_basis: recorded generation {generation} ({})",
                    self.prior_label
                );
                let _ = writeln!(
                    out,
                    "sources_identical_to_prior: {}",
                    if self.sources_identical_to_prior {
                        "yes (every unchanged cell is provably unchanged, not coincidentally equal)"
                    } else {
                        "no"
                    }
                );
            }
            None => {
                out.push_str(
                    "trend_basis: NO-DATA (no generation recorded yet; every trend cell reads \
                     `no prior generation`)\n",
                );
            }
        }
        let _ = writeln!(out, "source_identity: {}", hex(self.source_identity));
        out.push('\n');
        out.push_str(
            "This dashboard measures OUTCOMES. A `NO-DATA` row means no measurement machinery \
             exists yet, so no number is invented; a measured `0` means the population is \
             enumerable and the answer is genuinely none. The two are never conflated, and a \
             measured zero is deliberately left visible rather than hidden behind `NO-DATA`.\n\n",
        );

        for family in MetricFamily::ALL {
            let rows: Vec<&DashboardRow> = self
                .rows
                .iter()
                .filter(|entry| entry.row.family == family)
                .collect();
            if rows.is_empty() {
                continue;
            }
            let _ = writeln!(out, "## {}", family.title());
            out.push('\n');
            out.push_str("| metric | value | trend | direction | sources |\n");
            out.push_str("| --- | --- | --- | --- | --- |\n");
            for entry in rows {
                let value = match &entry.row.cell {
                    MetricCell::NoData { unblocked_by, .. } => match unblocked_by {
                        Some(bead) => format!("NO-DATA (needs {bead})"),
                        None => "NO-DATA".to_string(),
                    },
                    MetricCell::Measured(observation) => observation.render(),
                };
                let sources = if entry.row.sources.is_empty() {
                    "-".to_string()
                } else {
                    entry.row.sources.join("; ")
                };
                let _ = writeln!(
                    out,
                    "| {} | {} | {} | {} | {} |",
                    entry.row.title,
                    value,
                    entry.trend.render(),
                    entry.row.direction.slug(),
                    sources
                );
            }
            out.push('\n');
        }

        out.push_str("## What each metric does not capture\n\n");
        for entry in &self.rows {
            let _ = writeln!(out, "- `{}` — {}", entry.row.id, entry.row.honesty);
        }
        out.push('\n');

        out.push_str("## Why metrics are missing\n\n");
        let mut gaps = 0_usize;
        for entry in &self.rows {
            if let MetricCell::NoData {
                reason,
                unblocked_by,
            } = &entry.row.cell
            {
                gaps += 1;
                match unblocked_by {
                    Some(bead) => {
                        let _ = writeln!(out, "- `{}` — {reason} (tracked: {bead})", entry.row.id);
                    }
                    None => {
                        let _ = writeln!(out, "- `{}` — {reason}", entry.row.id);
                    }
                }
            }
        }
        if gaps == 0 {
            out.push_str("- none; every registered metric is measured\n");
        }
        out.push('\n');

        out.push_str("## Deliberately excluded\n\n");
        out.push_str(
            "These are legitimate signals that are NOT outcome metrics. They move without the \
             program getting better at predicting reality, and each has its own lane.\n\n",
        );
        for (name, reason) in DEMOTED_DIAGNOSTICS {
            let _ = writeln!(out, "- {name} — {reason}");
        }
        out.push('\n');

        let _ = writeln!(out, "identity: {}", hex(self.identity));
        out
    }

    /// Deterministic machine-readable dashboard.
    #[allow(clippy::too_many_lines)] // Field order is the schema; splitting it would obscure it.
    #[must_use]
    pub fn render_json(&self) -> String {
        let mut out = String::new();
        let _ = write!(
            out,
            "{{\"schema\":{DASHBOARD_SCHEMA_VERSION},\"metrics\":{},\"measured\":{},\
             \"no_data\":{}",
            self.rows.len(),
            self.measured,
            self.no_data
        );
        out.push_str(",\"trend_basis\":");
        match self.prior_generation {
            Some(generation) => {
                let _ = write!(out, "{{\"status\":\"recorded\",\"generation\":{generation}");
                out.push_str(",\"label\":");
                push_json_str(&mut out, &self.prior_label);
                let _ = write!(
                    out,
                    ",\"sources_identical\":{}}}",
                    self.sources_identical_to_prior
                );
            }
            None => out.push_str("{\"status\":\"no-data\",\"reason\":\"no generation recorded\"}"),
        }
        out.push_str(",\"source_identity\":");
        push_json_str(&mut out, &hex(self.source_identity));
        out.push_str(",\"rows\":[");
        for (index, entry) in self.rows.iter().enumerate() {
            if index != 0 {
                out.push(',');
            }
            out.push_str("{\"id\":");
            push_json_str(&mut out, &entry.row.id);
            out.push_str(",\"title\":");
            push_json_str(&mut out, &entry.row.title);
            out.push_str(",\"family\":");
            push_json_str(&mut out, entry.row.family.slug());
            out.push_str(",\"direction\":");
            push_json_str(&mut out, entry.row.direction.slug());
            out.push_str(",\"value\":");
            push_cell_json(&mut out, &entry.row.cell);
            out.push_str(",\"trend\":");
            push_json_str(&mut out, entry.trend.slug());
            out.push_str(",\"sources\":[");
            for (position, source) in entry.row.sources.iter().enumerate() {
                if position != 0 {
                    out.push(',');
                }
                push_json_str(&mut out, source);
            }
            out.push_str("],\"does_not_capture\":");
            push_json_str(&mut out, &entry.row.honesty);
            out.push('}');
        }
        out.push_str("],\"excluded\":[");
        for (index, (name, reason)) in DEMOTED_DIAGNOSTICS.iter().enumerate() {
            if index != 0 {
                out.push(',');
            }
            out.push_str("{\"measurement\":");
            push_json_str(&mut out, name);
            out.push_str(",\"reason\":");
            push_json_str(&mut out, reason);
            out.push('}');
        }
        out.push_str("],\"identity\":");
        push_json_str(&mut out, &hex(self.identity));
        out.push('}');
        out
    }
}

// ---------------------------------------------------------------------------
// FrankenSim projection
//
// The concrete metric set for THIS program, projected from registry values.
// This section takes typed registries and returns rows: it performs no I/O and
// has no filesystem, clock, or network access, so it structurally cannot
// fabricate a measurement that its inputs do not already contain.
// ---------------------------------------------------------------------------

/// Lowest maturity level counted as "numerically verified or better".
const VERIFIED_LEVEL: u8 = 2;

/// Lowest maturity level counted as "integrated workflow or better".
const INTEGRATED_LEVEL: u8 = 3;

fn level_rank(level: &str) -> Option<u8> {
    match level {
        "L1" => Some(1),
        "L2" => Some(2),
        "L3" => Some(3),
        "L4" => Some(4),
        "L5" => Some(5),
        _ => None,
    }
}

/// A ratio, or `NO-DATA` when the population is empty.
///
/// This is the structural guarantee stated in the module docs: an empty
/// denominator cannot become `0%`, because the only way to express a ratio is
/// with a nonzero denominator.
fn ratio_or_no_data(
    numerator: usize,
    denominator: usize,
    empty_reason: &str,
    unblocked_by: Option<&str>,
) -> MetricCell {
    match NonZeroU64::new(u64::try_from(denominator).unwrap_or(u64::MAX)) {
        Some(total) => MetricCell::Measured(MetricObservation::ratio(
            u64::try_from(numerator).unwrap_or(u64::MAX),
            total,
        )),
        None => MetricCell::NoData {
            reason: empty_reason.to_string(),
            unblocked_by: unblocked_by.map(str::to_string),
        },
    }
}

fn no_data(reason: &str, unblocked_by: Option<&str>) -> MetricCell {
    MetricCell::NoData {
        reason: reason.to_string(),
        unblocked_by: unblocked_by.map(str::to_string),
    }
}

/// Registry values the FrankenSim projection reads.
///
/// Deliberately a typed struct rather than paths: the projection cannot reach
/// past what a caller has already loaded and validated.
#[derive(Debug, Clone, Copy)]
pub struct ProgramSources {
    /// Whether the corpus supplied is the seeded public one.
    pub corpus_seeded: bool,
    /// Registered validation datasets.
    pub corpus_datasets: usize,
    /// Datasets supplying at least one EXTERNAL portfolio axis.
    pub external_datasets: usize,
    /// Datasets supplying the blind-predictive axis.
    pub blind_predictive_datasets: usize,
    /// Datasets supplying the independent-reproduction axis.
    pub independent_reproduction_datasets: usize,
    /// Declared (QoI, regime) claim cells.
    pub claim_cells: usize,
    /// Claim cells carrying at least one external reference.
    pub externally_anchored_cells: usize,
    /// Claim cells with no external reference at all.
    pub unanchored_cells: usize,
    /// Registered adversarial challenges.
    pub adversarial_cases: usize,
    /// Adversarial challenges actually executed.
    pub executed_assessments: usize,
    /// Registered capabilities.
    pub capabilities: usize,
    /// Capabilities at L2 or above.
    pub capabilities_verified: usize,
    /// Capabilities at L3 or above.
    pub capabilities_integrated: usize,
}

impl ProgramSources {
    /// Count capability levels, refusing an empty or uninterpretable registry.
    ///
    /// # Errors
    ///
    /// [`DashboardError::EmptySourceRegistry`] when the registry is empty —
    /// the read-failure trap — and [`DashboardError::UninterpretableSource`]
    /// for a level outside `L1..=L5`.
    pub fn count_levels(
        levels: &BTreeMap<String, String>,
    ) -> Result<(usize, usize), DashboardError> {
        if levels.is_empty() {
            return Err(DashboardError::EmptySourceRegistry {
                source: "capability-maturity",
            });
        }
        let mut verified = 0_usize;
        let mut integrated = 0_usize;
        for level in levels.values() {
            let Some(rank) = level_rank(level) else {
                return Err(DashboardError::UninterpretableSource {
                    source: "capability-maturity",
                    value: level.clone(),
                });
            };
            if rank >= VERIFIED_LEVEL {
                verified += 1;
            }
            if rank >= INTEGRATED_LEVEL {
                integrated += 1;
            }
        }
        Ok((verified, integrated))
    }

    /// Count the FrankenSim program sources straight out of the live registries.
    ///
    /// This is the only place registry shapes are turned into metric inputs, so
    /// the artifact lane and any harness that renders it read identically — a
    /// second implementation could drift, and a drifting dashboard is worse than
    /// no dashboard.
    ///
    /// # Errors
    ///
    /// Propagates [`ProgramSources::count_levels`] refusals: an empty or
    /// uninterpretable capability registry is a read failure, never zeros.
    pub fn from_registries(
        corpus: &CorpusRegistry,
        adversarial: &AdversarialRegistry,
        scorecard: &VvScorecard,
        capability_levels: &BTreeMap<String, String>,
    ) -> Result<Self, DashboardError> {
        let (capabilities_verified, capabilities_integrated) =
            Self::count_levels(capability_levels)?;
        let datasets = corpus.datasets();
        let on_axis = |axis: EvidenceAxis| {
            datasets
                .iter()
                .filter(|dataset| dataset.evidence_level().portfolio_axes().contains(&axis))
                .count()
        };
        let external_datasets = datasets
            .iter()
            .filter(|dataset| {
                let axes = dataset.evidence_level().portfolio_axes();
                EXTERNAL_AXES.iter().any(|axis| axes.contains(axis))
            })
            .count();
        let cells = scorecard.cells();
        let externally_anchored_cells = cells
            .iter()
            .filter(|cell| cell.external_datasets() != 0)
            .count();
        Ok(Self {
            corpus_seeded: corpus.is_seeded(),
            corpus_datasets: datasets.len(),
            external_datasets,
            blind_predictive_datasets: on_axis(EvidenceAxis::BlindPredictiveValidation),
            independent_reproduction_datasets: on_axis(EvidenceAxis::IndependentReproduction),
            claim_cells: cells.len(),
            externally_anchored_cells,
            unanchored_cells: cells.len() - externally_anchored_cells,
            adversarial_cases: adversarial.cases().len(),
            executed_assessments: scorecard.executed_assessments(),
            capabilities: capability_levels.len(),
            capabilities_verified,
            capabilities_integrated,
        })
    }
}

const CORPUS_SOURCE: &str = "fs_vvreg::corpus seeded validation registry";
const SCORECARD_SOURCE: &str = "vv-scorecard.json (fs_vvreg::scorecard)";
const ADVERSARIAL_SOURCE: &str = "fs_vvreg::adversarial registry";
const MATURITY_SOURCE: &str = "capability-maturity.json";

/// Project the FrankenSim program metric set from registry values.
///
/// Every row is either a real measurement with a cited source, or a loud
/// `NO-DATA` naming the tracked work that would make it live. A measured zero
/// is reported as zero: several rows below are genuinely zero today, and the
/// point of the artifact is that those stay visible.
///
/// # Errors
///
/// [`DashboardError::UnseededCorpus`] for a caller-built corpus, plus any
/// refusal from [`MetricRow::try_new`].
#[allow(clippy::too_many_lines)] // The metric set is a flat declaration; splitting it would scatter the policy.
pub fn frankensim_rows(sources: ProgramSources) -> Result<Vec<MetricRow>, DashboardError> {
    if !sources.corpus_seeded {
        return Err(DashboardError::UnseededCorpus);
    }
    if sources.corpus_datasets == 0 {
        return Err(DashboardError::EmptySourceRegistry {
            source: "vv-corpus",
        });
    }

    let rows = vec![
        // ---- Outcome ----
        MetricRow::try_new(
            "decision-turnaround",
            "Time from dirty CAD input to a defensible decision artifact",
            MetricFamily::Outcome,
            MetricDirection::LowerIsBetter,
            no_data(
                "no end-to-end acceptance lane records stage timings, so there is nothing to \
                 measure without inventing it",
                Some("f85xj.6.11"),
            ),
            &[],
            "when live this will time OUR examples on OUR hardware, which is not the same as a \
             real user's setup, data, or interruptions",
        )?,
        MetricRow::try_new(
            "blind-prediction-error",
            "Prediction error against blind held-out experimental references",
            MetricFamily::Outcome,
            MetricDirection::LowerIsBetter,
            no_data(
                "no ledgered run-result store exists and no Level-D blind reference is admitted, \
                 so no model-versus-reality error can be computed",
                Some("f85xj.7.5"),
            ),
            &[],
            "error against a reference bounds nothing on its own: the reference's own uncertainty \
             and the regime it was measured in both constrain what the number means",
        )?,
        MetricRow::try_new(
            "empirical-interval-coverage",
            "Empirical coverage of predicted uncertainty intervals",
            MetricFamily::Outcome,
            MetricDirection::HigherIsBetter,
            no_data(
                "the empirical coverage machinery is not live; nominal coverage is never \
                 extrapolated into an empirical claim",
                Some("f85xj.7.2"),
            ),
            &[],
            "coverage measured on the calibration population says nothing about coverage under \
             the distribution shift a real design study introduces",
        )?,
        MetricRow::try_new(
            "false-acceptance-rate",
            "Rate at which adversarial challenges were wrongly accepted",
            MetricFamily::Outcome,
            MetricDirection::LowerIsBetter,
            no_data(
                "zero registered adversarial challenges have been executed, so the rate's \
                 denominator is empty; a rate over zero trials is unrepresentable, not zero",
                Some("f85xj.7.5"),
            ),
            &[],
            "a false-acceptance rate only covers the failure modes someone thought to write a \
             challenge for; it is silent about unimagined ones",
        )?,
        MetricRow::try_new(
            "error-budget-completeness",
            "Share of user-facing outputs carrying a complete error budget",
            MetricFamily::Outcome,
            MetricDirection::HigherIsBetter,
            no_data(
                "no audit enumerates user-facing outputs against their budget terms",
                Some("f85xj.8.7"),
            ),
            &[],
            "a complete budget means every term the model KNOWS about is present; it cannot count \
             terms nobody has identified yet",
        )?,
        MetricRow::try_new(
            "import-admission-rate",
            "Supplier CAD import: clean, repaired, and refused rates",
            MetricFamily::Outcome,
            MetricDirection::HigherIsBetter,
            no_data(
                "no retained real supplier CAD corpus exists, and rates measured on fixtures we \
                 authored would be self-graded",
                Some("f85xj.11.6"),
            ),
            &[],
            "import success measures admission, not fidelity: a file that imports cleanly can \
             still carry geometry that means something different downstream",
        )?,
        MetricRow::try_new(
            "surrogate-escalation-correctness",
            "Correctness of certify-or-escalate decisions by learned components",
            MetricFamily::Outcome,
            MetricDirection::HigherIsBetter,
            no_data(
                "no decisive-metrics instrumentation exists for learned components",
                Some("f85xj.14.2"),
            ),
            &[],
            "escalation correctness is measured against cases where the truth is known, which are \
             systematically the easier ones",
        )?,
        MetricRow::try_new(
            "time-to-explain",
            "Time to explain a surprising result through ledger lineage",
            MetricFamily::Outcome,
            MetricDirection::LowerIsBetter,
            no_data(
                "the ledger exposes no explain-query session instrumentation, so no timing \
                 surface exists to read",
                None,
            ),
            &[],
            "even when live this measures the tool's explain path, not whether the explanation \
             actually convinced the engineer reading it",
        )?,
        MetricRow::try_new(
            "decision-changes-from-omitted-uncertainty",
            "Compliance verdicts flipped by adding a previously omitted budget term",
            MetricFamily::Outcome,
            MetricDirection::Neutral,
            no_data(
                "no ledger diff of decision verdicts across budget-term introductions is \
                 retained; this is the single best evidence that the error-budget program has \
                 decision value, and it is not yet collected",
                Some("f85xj.8.7"),
            ),
            &[],
            "a flipped verdict shows the term mattered on re-run; it does not show the new \
             verdict is correct, only that the earlier one was underdetermined",
        )?,
        MetricRow::try_new(
            "user-study-measurements",
            "Setup time, diagnosis time, and decision quality from real user sessions",
            MetricFamily::Outcome,
            MetricDirection::HigherIsBetter,
            no_data(
                "no user-study measurement exists; the nearest current proxies are quickstart \
                 timings and the external-reproduction friction log, neither of which is a user \
                 study",
                Some("f85xj.7.6"),
            ),
            &[],
            "proxies measured on people who already know the system are the opposite of the \
             population this metric is supposed to describe",
        )?,
        MetricRow::try_new(
            "cross-machine-reproducibility",
            "Share of claims replayed bitwise on an independent machine",
            MetricFamily::Outcome,
            MetricDirection::HigherIsBetter,
            no_data(
                "cross-ISA determinism is proven per-artifact by golden couplings and per-host by \
                 perf baselines, but no lane aggregates them into a program-level replay rate",
                None,
            ),
            &[],
            "bitwise replay on a second machine proves determinism, not correctness: two machines \
             can reproduce the same wrong number exactly",
        )?,
        MetricRow::try_new(
            "independent-reproduction",
            "Datasets reproduced by an independent team or implementation lineage",
            MetricFamily::Outcome,
            MetricDirection::HigherIsBetter,
            ratio_or_no_data(
                sources.independent_reproduction_datasets,
                sources.corpus_datasets,
                "no validation dataset is registered",
                Some("f85xj.7.6"),
            ),
            &[CORPUS_SOURCE],
            "this counts datasets DECLARING the independent-reproduction axis; the declaration is \
             a registry fact, and the current value is a genuine zero rather than an unmeasured \
             one",
        )?,
        // ---- Portfolio ----
        MetricRow::try_new(
            "externally-anchored-claim-cells",
            "Claim cells carrying at least one external reference",
            MetricFamily::Portfolio,
            MetricDirection::HigherIsBetter,
            ratio_or_no_data(
                sources.externally_anchored_cells,
                sources.claim_cells,
                "no (QoI, regime) claim cell is declared",
                None,
            ),
            &[SCORECARD_SOURCE],
            "anchoring counts REFERENCES attached to a cell, not agreement with them; a cell can \
             be anchored and still predict the reference badly",
        )?,
        MetricRow::try_new(
            "external-reference-datasets",
            "Validation datasets supplying an external evidence axis",
            MetricFamily::Portfolio,
            MetricDirection::HigherIsBetter,
            ratio_or_no_data(
                sources.external_datasets,
                sources.corpus_datasets,
                "no validation dataset is registered",
                None,
            ),
            &[CORPUS_SOURCE],
            "external means cross-code, controlled-experimental, blind-predictive, or field \
             monitoring; our own numerical verification is deliberately excluded from the \
             numerator because agreeing with ourselves is not external evidence",
        )?,
        MetricRow::try_new(
            "blind-predictive-datasets",
            "Validation datasets on the blind-predictive axis",
            MetricFamily::Portfolio,
            MetricDirection::HigherIsBetter,
            ratio_or_no_data(
                sources.blind_predictive_datasets,
                sources.corpus_datasets,
                "no validation dataset is registered",
                Some("f85xj.4.9"),
            ),
            &[CORPUS_SOURCE],
            "the honesty exam: a prediction made before the reference was unblinded. The current \
             value is a real zero, which is the single most important number on this dashboard",
        )?,
        MetricRow::try_new(
            "adversarial-suite-execution",
            "Registered adversarial challenges actually executed",
            MetricFamily::Portfolio,
            MetricDirection::HigherIsBetter,
            ratio_or_no_data(
                sources.executed_assessments,
                sources.adversarial_cases,
                "no adversarial challenge is registered",
                Some("f85xj.7.5"),
            ),
            &[ADVERSARIAL_SOURCE, SCORECARD_SOURCE],
            "execution is not survival: running a challenge says nothing about whether the \
             program passed it, and an unexecuted suite is a registry of good intentions",
        )?,
        MetricRow::try_new(
            "unanchored-claim-cells",
            "Claim cells with no external reference at all",
            MetricFamily::Portfolio,
            MetricDirection::LowerIsBetter,
            MetricCell::Measured(MetricObservation::count(
                u64::try_from(sources.unanchored_cells).unwrap_or(u64::MAX),
            )),
            &[SCORECARD_SOURCE],
            "an absolute count, so it grows when the program declares new claim cells; a rising \
             number can mean expanding scope rather than decaying evidence",
        )?,
        // ---- Governance ----
        MetricRow::try_new(
            "capabilities-at-l2-plus",
            "Registered capabilities at L2 (numerically verified) or above",
            MetricFamily::Governance,
            MetricDirection::HigherIsBetter,
            ratio_or_no_data(
                sources.capabilities_verified,
                sources.capabilities,
                "no capability is registered",
                None,
            ),
            &[MATURITY_SOURCE],
            "registry levels are declarations backed by cited evidence, not independent audits; \
             the registry's own maturity is L1",
        )?,
        MetricRow::try_new(
            "capabilities-at-l3-plus",
            "Registered capabilities at L3 (integrated workflow) or above",
            MetricFamily::Governance,
            MetricDirection::HigherIsBetter,
            ratio_or_no_data(
                sources.capabilities_integrated,
                sources.capabilities,
                "no capability is registered",
                None,
            ),
            &[MATURITY_SOURCE],
            "L3 requires an admitted end-to-end integration claim; the current value is a real \
             zero, and no crate count or test count can move it",
        )?,
    ];
    Ok(rows)
}
