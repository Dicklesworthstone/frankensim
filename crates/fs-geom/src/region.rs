//! `Region`: one abstract shape presented through many charts, with
//! agreement as a CHECKABLE PROPOSITION (plan §7.1) — sampled, seeded,
//! deterministic, and localized when it fails.

use crate::{Chart, Point3};
use fs_evidence::{NumericalCertificate, NumericalKind, ProvenanceHash};
use fs_exec::{Cancelled, Cx};
use std::fmt::Write as _;
use std::sync::Arc;

/// One chart entry: the presentation plus how it was obtained.
#[derive(Clone)]
pub struct RegionChart {
    /// The chart.
    pub chart: Arc<dyn Chart>,
    /// Content-address of the operation that produced this presentation.
    pub provenance: ProvenanceHash,
}

/// An abstract region: one or more charts + provenance. The FIRST chart is
/// the deterministic primary (routing beads add cost-aware selection).
#[derive(Clone)]
pub struct Region {
    charts: Vec<RegionChart>,
}

impl Region {
    /// A region presented by one chart.
    #[must_use]
    pub fn from_chart(chart: Arc<dyn Chart>, provenance: ProvenanceHash) -> Self {
        Region {
            charts: vec![RegionChart { chart, provenance }],
        }
    }

    /// Add another presentation of the SAME abstract region.
    #[must_use]
    pub fn with_chart(mut self, chart: Arc<dyn Chart>, provenance: ProvenanceHash) -> Self {
        self.charts.push(RegionChart { chart, provenance });
        self
    }

    /// The presentations, in insertion order.
    #[must_use]
    pub fn charts(&self) -> &[RegionChart] {
        &self.charts
    }

    /// The deterministic primary presentation.
    #[must_use]
    pub fn primary(&self) -> &RegionChart {
        &self.charts[0]
    }

    /// Check pairwise chart agreement on seeded sample points over the
    /// union support box. Two charts agree at `x` when their declared
    /// signed-distance intervals overlap after `config.tolerance_abs` slack;
    /// a disagreement therefore means at least one chart's geometry OR error
    /// declaration is wrong. Missing, malformed, or non-finite evidence makes
    /// the proposition [`AgreementStatus::Unknown`], never agreed.
    /// Deterministic: same seed, same report. Polls cancellation per
    /// sample (Decalogue P7).
    ///
    /// # Errors
    /// [`Cancelled`] when the context's gate was requested mid-check.
    #[allow(clippy::too_many_lines)] // One ordered seeded sampling and reduction protocol.
    pub fn check_agreement(
        &self,
        config: &AgreementConfig,
        cx: &Cx<'_>,
    ) -> Result<AgreementReport, Cancelled> {
        let mut report = AgreementReport::empty();
        if config.samples == 0 {
            report.record_unknown(
                AgreementUnknown::global(AgreementUnknownReason::ZeroSamples),
                config.max_diagnostics,
            );
        }
        if self.charts.len() < 2 {
            report.record_unknown(
                AgreementUnknown::global(AgreementUnknownReason::InsufficientCharts {
                    found: self.charts.len(),
                }),
                config.max_diagnostics,
            );
        }
        if !config.tolerance_abs.is_finite() || config.tolerance_abs < 0.0 {
            report.record_unknown(
                AgreementUnknown::global(AgreementUnknownReason::InvalidTolerance {
                    bits: config.tolerance_abs.to_bits(),
                }),
                config.max_diagnostics,
            );
        }
        if report.unknown_count > 0 {
            return Ok(report);
        }

        let mut support: Option<crate::Aabb> = None;
        for rc in &self.charts {
            let chart_support = rc.chart.support();
            if !valid_support(chart_support) {
                report.record_unknown(
                    AgreementUnknown::chart(
                        rc.chart.name(),
                        AgreementUnknownReason::InvalidSupport,
                    ),
                    config.max_diagnostics,
                );
                continue;
            }
            support = Some(match support {
                Some(accumulated) => accumulated.union(&chart_support),
                None => chart_support,
            });
        }
        if report.unknown_count > 0 {
            return Ok(report);
        }
        let Some(support) = support else {
            report.record_unknown(
                AgreementUnknown::global(AgreementUnknownReason::InvalidSupport),
                config.max_diagnostics,
            );
            return Ok(report);
        };
        let span = (support.max.x - support.min.x)
            .max(support.max.y - support.min.y)
            .max(support.max.z - support.min.z);
        let box_ = support.inflate(0.1 * span.max(1e-3));
        if !valid_support(box_) {
            report.record_unknown(
                AgreementUnknown::global(AgreementUnknownReason::InvalidSupport),
                config.max_diagnostics,
            );
            return Ok(report);
        }

        let mut state = config.seed | 1;
        let mut unit = move || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((state >> 11) as f64) / (1u64 << 53) as f64
        };
        for _ in 0..config.samples {
            cx.checkpoint()?;
            let p = Point3::new(
                box_.min.x + (box_.max.x - box_.min.x) * unit(),
                box_.min.y + (box_.max.y - box_.min.y) * unit(),
                box_.min.z + (box_.max.z - box_.min.z) * unit(),
            );
            if !finite_point(p) {
                report.record_unknown(
                    AgreementUnknown::at(p, AgreementUnknownReason::NonFiniteSamplePoint),
                    config.max_diagnostics,
                );
                continue;
            }

            // Evaluate each chart once. Besides avoiding repeated work, this
            // prevents a stateful or otherwise malformed chart from returning
            // different evidence for different pairs at the same point.
            let mut samples = Vec::with_capacity(self.charts.len());
            for rc in &self.charts {
                let sample = rc.chart.eval(p, cx);
                match validate_sample(&sample) {
                    Ok(()) => samples.push(Some(sample)),
                    Err(reason) => {
                        report.record_unknown(
                            AgreementUnknown::chart_at(rc.chart.name(), p, reason),
                            config.max_diagnostics,
                        );
                        samples.push(None);
                    }
                }
            }

            for i in 0..self.charts.len() {
                for j in (i + 1)..self.charts.len() {
                    let (Some(a), Some(b)) = (samples[i], samples[j]) else {
                        continue;
                    };
                    let Some((gap, allowed, excess)) = comparison(&a, &b, config.tolerance_abs)
                    else {
                        report.record_unknown(
                            AgreementUnknown::pair_at(
                                self.charts[i].chart.name(),
                                self.charts[j].chart.name(),
                                p,
                                AgreementUnknownReason::NonFiniteComparison,
                            ),
                            config.max_diagnostics,
                        );
                        continue;
                    };
                    let evidence = a.error.kind.max(b.error.kind);
                    report.checked += 1;
                    report.worst_excess = Some(
                        report
                            .worst_excess
                            .map_or(excess, |current| current.max(excess)),
                    );
                    report.weakest_evidence = Some(
                        report
                            .weakest_evidence
                            .map_or(evidence, |current| current.max(evidence)),
                    );
                    if excess > 0.0 {
                        report.disagreement_count += 1;
                        report.strongest_counterexample_evidence = Some(
                            report
                                .strongest_counterexample_evidence
                                .map_or(evidence, |current| current.min(evidence)),
                        );
                        if report.disagreements.len() < config.max_diagnostics {
                            report.disagreements.push(Disagreement {
                                at: p,
                                chart_a: self.charts[i].chart.name(),
                                chart_b: self.charts[j].chart.name(),
                                sd_a: a.signed_distance,
                                sd_b: b.signed_distance,
                                gap,
                                allowed,
                                evidence,
                            });
                        }
                    }
                }
            }
        }
        report.status = if report.disagreement_count > 0 {
            AgreementStatus::Disagreed
        } else if report.unknown_count > 0 || report.checked == 0 {
            AgreementStatus::Unknown
        } else {
            AgreementStatus::Agreed
        };
        Ok(report)
    }
}

fn finite_point(p: Point3) -> bool {
    p.x.is_finite() && p.y.is_finite() && p.z.is_finite()
}

fn valid_support(support: crate::Aabb) -> bool {
    finite_point(support.min)
        && finite_point(support.max)
        && support.min.x <= support.max.x
        && support.min.y <= support.max.y
        && support.min.z <= support.max.z
        && (support.max.x - support.min.x).is_finite()
        && (support.max.y - support.min.y).is_finite()
        && (support.max.z - support.min.z).is_finite()
}

#[allow(clippy::float_cmp)] // An Exact certificate must have a numerically singleton interval.
fn validate_sample(sample: &crate::ChartSample) -> Result<(), AgreementUnknownReason> {
    if !sample.signed_distance.is_finite() {
        return Err(AgreementUnknownReason::NonFiniteSignedDistance {
            bits: sample.signed_distance.to_bits(),
        });
    }
    if let Some(gradient) = sample.gradient {
        for (component, value) in [("x", gradient.x), ("y", gradient.y), ("z", gradient.z)] {
            if !value.is_finite() {
                return Err(AgreementUnknownReason::NonFiniteGradient {
                    component,
                    bits: value.to_bits(),
                });
            }
        }
    }
    if let Some(lipschitz) = sample.lipschitz
        && (!lipschitz.is_finite() || lipschitz < 0.0)
    {
        return Err(AgreementUnknownReason::InvalidLipschitz {
            bits: lipschitz.to_bits(),
        });
    }
    if sample.error.kind == NumericalKind::NoClaim {
        return Err(AgreementUnknownReason::NoClaim);
    }
    let NumericalCertificate { kind, lo, hi } = sample.error;
    if !lo.is_finite()
        || !hi.is_finite()
        || lo > hi
        || !(hi - lo).is_finite()
        || (kind == NumericalKind::Exact && lo != hi)
    {
        return Err(AgreementUnknownReason::MalformedCertificate {
            lo_bits: lo.to_bits(),
            hi_bits: hi.to_bits(),
        });
    }
    if sample.signed_distance < lo || sample.signed_distance > hi {
        return Err(AgreementUnknownReason::ValueOutsideCertificate {
            value_bits: sample.signed_distance.to_bits(),
            lo_bits: lo.to_bits(),
            hi_bits: hi.to_bits(),
        });
    }
    Ok(())
}

fn comparison(
    a: &crate::ChartSample,
    b: &crate::ChartSample,
    tolerance_abs: f64,
) -> Option<(f64, f64, f64)> {
    let gap = (a.signed_distance - b.signed_distance).abs();
    let allowed = if a.signed_distance <= b.signed_distance {
        (a.error.hi - a.signed_distance) + (b.signed_distance - b.error.lo) + tolerance_abs
    } else {
        (a.signed_distance - a.error.lo) + (b.error.hi - b.signed_distance) + tolerance_abs
    };
    let excess = gap - allowed;
    (gap.is_finite() && allowed.is_finite() && excess.is_finite()).then_some((gap, allowed, excess))
}

/// Agreement-check configuration (seeded — the check is replayable).
#[derive(Debug, Clone)]
pub struct AgreementConfig {
    /// Sample-point count.
    pub samples: u64,
    /// Replay seed.
    pub seed: u64,
    /// Extra absolute slack beyond the charts' declared errors.
    pub tolerance_abs: f64,
    /// Cap on localized diagnostics kept (first-K, deterministic).
    pub max_diagnostics: usize,
}

impl Default for AgreementConfig {
    fn default() -> Self {
        AgreementConfig {
            samples: 2_000,
            seed: 0x9E0_A62E,
            tolerance_abs: 1e-9,
            max_diagnostics: 8,
        }
    }
}

/// One localized disagreement: the point, the charts, the numbers — the
/// diagnostic an agent needs to decide WHICH presentation to distrust.
#[derive(Debug, Clone, PartialEq)]
pub struct Disagreement {
    /// Where the charts disagree.
    pub at: Point3,
    /// First chart's name.
    pub chart_a: &'static str,
    /// Second chart's name.
    pub chart_b: &'static str,
    /// First chart's signed distance.
    pub sd_a: f64,
    /// Second chart's signed distance.
    pub sd_b: f64,
    /// The observed gap.
    pub gap: f64,
    /// The gap the declared errors would have allowed.
    pub allowed: f64,
    /// Weakest certificate class used by this counterexample.
    pub evidence: NumericalKind,
}

/// A sampled agreement check's three-valued verdict.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AgreementStatus {
    /// Every requested pairwise comparison was valid and no counterexample
    /// was found within the declared intervals and configured slack.
    Agreed,
    /// At least one valid comparison produced a counterexample.
    Disagreed,
    /// The requested proposition could not be evaluated completely.
    Unknown,
}

impl AgreementStatus {
    fn name(self) -> &'static str {
        match self {
            AgreementStatus::Agreed => "agreed",
            AgreementStatus::Disagreed => "disagreed",
            AgreementStatus::Unknown => "unknown",
        }
    }
}

/// Why an agreement proposition could not be evaluated.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AgreementUnknownReason {
    /// Agreement requires at least two presentations. This check does not
    /// establish that their implementations or provenance are independent.
    InsufficientCharts {
        /// Number of presentations available.
        found: usize,
    },
    /// An empty sample set provides no evidence.
    ZeroSamples,
    /// Slack must be finite and non-negative; the exact input bits are kept.
    InvalidTolerance {
        /// Exact IEEE-754 bits supplied by the caller.
        bits: u64,
    },
    /// A chart's support was non-finite, inverted, or too wide for finite
    /// sampling arithmetic.
    InvalidSupport,
    /// Support interpolation produced a non-finite query point.
    NonFiniteSamplePoint,
    /// A chart returned a non-finite signed distance.
    NonFiniteSignedDistance {
        /// Exact IEEE-754 bits returned by the chart.
        bits: u64,
    },
    /// A chart returned a non-finite gradient component.
    NonFiniteGradient {
        /// Component name (`x`, `y`, or `z`).
        component: &'static str,
        /// Exact IEEE-754 bits returned by the chart.
        bits: u64,
    },
    /// A claimed Lipschitz bound was negative or non-finite.
    InvalidLipschitz {
        /// Exact IEEE-754 bits returned by the chart.
        bits: u64,
    },
    /// A chart explicitly made no numerical claim.
    NoClaim,
    /// Certificate bounds were non-finite, inverted, impossibly wide, or an
    /// `Exact` certificate had nonzero width.
    MalformedCertificate {
        /// Exact IEEE-754 bits of the lower bound.
        lo_bits: u64,
        /// Exact IEEE-754 bits of the upper bound.
        hi_bits: u64,
    },
    /// The reported signed distance was outside its own declared interval.
    ValueOutsideCertificate {
        /// Exact IEEE-754 bits of the reported value.
        value_bits: u64,
        /// Exact IEEE-754 bits of the lower bound.
        lo_bits: u64,
        /// Exact IEEE-754 bits of the upper bound.
        hi_bits: u64,
    },
    /// Otherwise-valid finite values overflowed the comparison arithmetic.
    NonFiniteComparison,
}

impl AgreementUnknownReason {
    fn code(&self) -> &'static str {
        match self {
            Self::InsufficientCharts { .. } => "insufficient-charts",
            Self::ZeroSamples => "zero-samples",
            Self::InvalidTolerance { .. } => "invalid-tolerance",
            Self::InvalidSupport => "invalid-support",
            Self::NonFiniteSamplePoint => "non-finite-sample-point",
            Self::NonFiniteSignedDistance { .. } => "non-finite-signed-distance",
            Self::NonFiniteGradient { .. } => "non-finite-gradient",
            Self::InvalidLipschitz { .. } => "invalid-lipschitz",
            Self::NoClaim => "no-claim",
            Self::MalformedCertificate { .. } => "malformed-certificate",
            Self::ValueOutsideCertificate { .. } => "value-outside-certificate",
            Self::NonFiniteComparison => "non-finite-comparison",
        }
    }
}

/// A localized reason the agreement proposition is unknown.
#[derive(Debug, Clone, PartialEq)]
pub struct AgreementUnknown {
    /// Query point, when the failure arose during sampling.
    pub at: Option<Point3>,
    /// First implicated chart, when any.
    pub chart_a: Option<&'static str>,
    /// Second implicated chart for pairwise failures.
    pub chart_b: Option<&'static str>,
    /// Structured failure class and exact non-finite input bits where useful.
    pub reason: AgreementUnknownReason,
}

impl AgreementUnknown {
    fn global(reason: AgreementUnknownReason) -> Self {
        Self {
            at: None,
            chart_a: None,
            chart_b: None,
            reason,
        }
    }

    fn at(at: Point3, reason: AgreementUnknownReason) -> Self {
        Self {
            at: Some(at),
            ..Self::global(reason)
        }
    }

    fn chart(chart: &'static str, reason: AgreementUnknownReason) -> Self {
        Self {
            chart_a: Some(chart),
            ..Self::global(reason)
        }
    }

    fn chart_at(chart: &'static str, at: Point3, reason: AgreementUnknownReason) -> Self {
        Self {
            at: Some(at),
            chart_a: Some(chart),
            chart_b: None,
            reason,
        }
    }

    fn pair_at(
        chart_a: &'static str,
        chart_b: &'static str,
        at: Point3,
        reason: AgreementUnknownReason,
    ) -> Self {
        Self {
            at: Some(at),
            chart_a: Some(chart_a),
            chart_b: Some(chart_b),
            reason,
        }
    }
}

/// The agreement verdict with numerical strength and localized diagnostics.
#[derive(Debug, Clone, PartialEq)]
pub struct AgreementReport {
    /// Three-valued result. Invalid or absent evidence is never agreement.
    pub status: AgreementStatus,
    /// Pairwise comparisons performed.
    pub checked: u64,
    /// Worst observed `gap - allowed`; negative values preserve the smallest
    /// observed safety margin. `None` means no valid comparison occurred.
    pub worst_excess: Option<f64>,
    /// Weakest certificate class used by any valid comparison.
    pub weakest_evidence: Option<NumericalKind>,
    /// Strongest certificate class among observed counterexamples.
    pub strongest_counterexample_evidence: Option<NumericalKind>,
    /// Total counterexamples, including those beyond the diagnostic cap.
    pub disagreement_count: u64,
    /// Total unknown observations, including those beyond the diagnostic cap.
    pub unknown_count: u64,
    /// First-K localized disagreements.
    pub disagreements: Vec<Disagreement>,
    /// First-K localized reasons the proposition was not fully evaluable.
    pub unknowns: Vec<AgreementUnknown>,
}

impl AgreementReport {
    fn empty() -> Self {
        Self {
            status: AgreementStatus::Unknown,
            checked: 0,
            worst_excess: None,
            weakest_evidence: None,
            strongest_counterexample_evidence: None,
            disagreement_count: 0,
            unknown_count: 0,
            disagreements: Vec::new(),
            unknowns: Vec::new(),
        }
    }

    fn record_unknown(&mut self, unknown: AgreementUnknown, max_diagnostics: usize) {
        self.unknown_count += 1;
        if self.unknowns.len() < max_diagnostics {
            self.unknowns.push(unknown);
        }
    }

    /// True only for a fully evaluated check with no counterexamples.
    #[must_use]
    pub fn is_agreed(&self) -> bool {
        self.status == AgreementStatus::Agreed
    }

    /// Canonical JSON (deterministic — replayable evidence for the
    /// explain() payload).
    #[must_use]
    pub fn to_json(&self) -> String {
        let mut s = String::with_capacity(128);
        let _ = write!(
            s,
            "{{\"checked\":{},\"status\":\"{}\",\"worst_excess\":",
            self.checked,
            self.status.name()
        );
        write_optional_f64(&mut s, self.worst_excess);
        s.push_str(",\"weakest_evidence\":");
        write_optional_kind(&mut s, self.weakest_evidence);
        s.push_str(",\"strongest_counterexample_evidence\":");
        write_optional_kind(&mut s, self.strongest_counterexample_evidence);
        let _ = write!(
            s,
            ",\"disagreement_count\":{},\"unknown_count\":{},\"disagreements\":[",
            self.disagreement_count, self.unknown_count
        );
        for (i, d) in self.disagreements.iter().enumerate() {
            if i > 0 {
                s.push(',');
            }
            s.push_str("{\"at\":[");
            write_json_f64(&mut s, d.at.x);
            s.push(',');
            write_json_f64(&mut s, d.at.y);
            s.push(',');
            write_json_f64(&mut s, d.at.z);
            s.push_str("],\"charts\":[");
            write_json_string(&mut s, d.chart_a);
            s.push(',');
            write_json_string(&mut s, d.chart_b);
            let _ = write!(
                s,
                "],\"sd\":[{},{}],\"gap\":{},\"allowed\":{},\"evidence\":\"{}\"}}",
                d.sd_a,
                d.sd_b,
                d.gap,
                d.allowed,
                kind_name(d.evidence)
            );
        }
        s.push_str("],\"unknowns\":[");
        for (i, unknown) in self.unknowns.iter().enumerate() {
            if i > 0 {
                s.push(',');
            }
            write_unknown(&mut s, unknown);
        }
        s.push_str("]}");
        s
    }
}

fn kind_name(kind: NumericalKind) -> &'static str {
    match kind {
        NumericalKind::Exact => "exact",
        NumericalKind::Enclosure => "enclosure",
        NumericalKind::Estimate => "estimate",
        NumericalKind::NoClaim => "no-claim",
    }
}

fn write_optional_kind(out: &mut String, kind: Option<NumericalKind>) {
    if let Some(kind) = kind {
        let _ = write!(out, "\"{}\"", kind_name(kind));
    } else {
        out.push_str("null");
    }
}

fn write_optional_f64(out: &mut String, value: Option<f64>) {
    if let Some(value) = value {
        debug_assert!(value.is_finite());
        let _ = write!(out, "{value}");
    } else {
        out.push_str("null");
    }
}

fn write_json_f64(out: &mut String, value: f64) {
    if value.is_finite() {
        let _ = write!(out, "{value}");
    } else {
        let _ = write!(out, "\"bits:{:016x}\"", value.to_bits());
    }
}

fn write_json_string(out: &mut String, value: &str) {
    out.push('"');
    for ch in value.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\u{08}' => out.push_str("\\b"),
            '\u{0c}' => out.push_str("\\f"),
            c if c <= '\u{1f}' => {
                let _ = write!(out, "\\u{:04x}", u32::from(c));
            }
            c => out.push(c),
        }
    }
    out.push('"');
}

fn write_unknown(out: &mut String, unknown: &AgreementUnknown) {
    out.push_str("{\"reason\":");
    write_json_string(out, unknown.reason.code());
    if let Some(at) = unknown.at {
        out.push_str(",\"at\":[");
        write_json_f64(out, at.x);
        out.push(',');
        write_json_f64(out, at.y);
        out.push(',');
        write_json_f64(out, at.z);
        out.push(']');
    }
    if let Some(chart) = unknown.chart_a {
        out.push_str(",\"chart_a\":");
        write_json_string(out, chart);
    }
    if let Some(chart) = unknown.chart_b {
        out.push_str(",\"chart_b\":");
        write_json_string(out, chart);
    }
    match &unknown.reason {
        AgreementUnknownReason::InsufficientCharts { found } => {
            let _ = write!(out, ",\"found\":{found}");
        }
        AgreementUnknownReason::InvalidTolerance { bits }
        | AgreementUnknownReason::NonFiniteSignedDistance { bits }
        | AgreementUnknownReason::InvalidLipschitz { bits } => {
            let _ = write!(out, ",\"bits\":\"{bits:016x}\"");
        }
        AgreementUnknownReason::NonFiniteGradient { component, bits } => {
            let _ = write!(
                out,
                ",\"component\":\"{component}\",\"bits\":\"{bits:016x}\""
            );
        }
        AgreementUnknownReason::MalformedCertificate { lo_bits, hi_bits } => {
            let _ = write!(
                out,
                ",\"lo_bits\":\"{lo_bits:016x}\",\"hi_bits\":\"{hi_bits:016x}\""
            );
        }
        AgreementUnknownReason::ValueOutsideCertificate {
            value_bits,
            lo_bits,
            hi_bits,
        } => {
            let _ = write!(
                out,
                ",\"value_bits\":\"{value_bits:016x}\",\"lo_bits\":\"{lo_bits:016x}\",\"hi_bits\":\"{hi_bits:016x}\""
            );
        }
        AgreementUnknownReason::ZeroSamples
        | AgreementUnknownReason::InvalidSupport
        | AgreementUnknownReason::NonFiniteSamplePoint
        | AgreementUnknownReason::NoClaim
        | AgreementUnknownReason::NonFiniteComparison => {}
    }
    out.push('}');
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fixtures::{LyingSphereChart, SphereChart};
    use asupersync::types::Budget;
    use fs_exec::{CancelGate, ExecMode, StreamKey};

    fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                gate,
                arena,
                StreamKey {
                    seed: 7,
                    kernel_id: 7,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            f(&cx)
        })
    }

    fn sphere() -> SphereChart {
        SphereChart {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.5,
        }
    }

    #[test]
    fn identical_charts_agree_and_reports_are_deterministic() {
        let region = Region::from_chart(Arc::new(sphere()), ProvenanceHash::of_bytes(b"a"))
            .with_chart(Arc::new(sphere()), ProvenanceHash::of_bytes(b"b"));
        let gate = CancelGate::new();
        let (r1, r2) = with_cx(&gate, |cx| {
            let cfg = AgreementConfig::default();
            (
                region.check_agreement(&cfg, cx).expect("not cancelled"),
                region.check_agreement(&cfg, cx).expect("not cancelled"),
            )
        });
        assert_eq!(r1.status, AgreementStatus::Agreed);
        assert!(r1.disagreements.is_empty() && r1.unknowns.is_empty());
        assert_eq!(r1.weakest_evidence, Some(NumericalKind::Exact));
        assert!(r1.worst_excess.is_some_and(|excess| excess < 0.0));
        assert_eq!(r1.to_json(), r2.to_json(), "seeded check replays (G5)");
    }

    #[test]
    fn lying_chart_is_detected_with_localized_diagnostics() {
        let region = Region::from_chart(Arc::new(sphere()), ProvenanceHash::of_bytes(b"honest"))
            .with_chart(
                Arc::new(LyingSphereChart {
                    sphere: sphere(),
                    bias: 0.05,
                }),
                ProvenanceHash::of_bytes(b"liar"),
            );
        let gate = CancelGate::new();
        let report = with_cx(&gate, |cx| {
            region
                .check_agreement(&AgreementConfig::default(), cx)
                .expect("not cancelled")
        });
        assert_eq!(
            report.status,
            AgreementStatus::Disagreed,
            "the bias must be caught"
        );
        assert_eq!(
            report.strongest_counterexample_evidence,
            Some(NumericalKind::Exact)
        );
        let d = &report.disagreements[0];
        assert!(d.gap > d.allowed);
        assert!(
            d.chart_a == "fixture/lying-sphere" || d.chart_b == "fixture/lying-sphere",
            "diagnostic names the liar: {d:?}"
        );
        assert!((d.gap - 0.05).abs() < 1e-9, "gap localizes the bias size");
    }

    #[test]
    fn agreement_check_is_cancellable() {
        let region = Region::from_chart(Arc::new(sphere()), ProvenanceHash::of_bytes(b"a"))
            .with_chart(Arc::new(sphere()), ProvenanceHash::of_bytes(b"b"));
        let gate = CancelGate::new();
        gate.request();
        let outcome = with_cx(&gate, |cx| {
            region.check_agreement(&AgreementConfig::default(), cx)
        });
        assert_eq!(outcome, Err(Cancelled), "pre-requested gate cancels");
    }

    #[derive(Clone, Copy)]
    struct ProbeChart {
        name: &'static str,
        sample: crate::ChartSample,
        support: crate::Aabb,
    }

    impl Chart for ProbeChart {
        fn eval(&self, _x: Point3, _cx: &Cx<'_>) -> crate::ChartSample {
            self.sample
        }

        fn support(&self) -> crate::Aabb {
            self.support
        }

        fn name(&self) -> &'static str {
            self.name
        }
    }

    fn probe(sample: crate::ChartSample) -> ProbeChart {
        ProbeChart {
            name: "probe",
            sample,
            support: crate::Aabb::new(Point3::new(-2.0, -2.0, -2.0), Point3::new(2.0, 2.0, 2.0)),
        }
    }

    fn region_with_probe(probe: ProbeChart) -> Region {
        Region::from_chart(Arc::new(sphere()), ProvenanceHash::of_bytes(b"honest"))
            .with_chart(Arc::new(probe), ProvenanceHash::of_bytes(b"probe"))
    }

    #[test]
    fn zero_samples_and_one_chart_are_unknown_not_vacuously_agreed() {
        let gate = CancelGate::new();
        let zero = with_cx(&gate, |cx| {
            let region = Region::from_chart(Arc::new(sphere()), ProvenanceHash::of_bytes(b"a"))
                .with_chart(Arc::new(sphere()), ProvenanceHash::of_bytes(b"b"));
            region
                .check_agreement(
                    &AgreementConfig {
                        samples: 0,
                        ..AgreementConfig::default()
                    },
                    cx,
                )
                .expect("not cancelled")
        });
        assert_eq!(zero.status, AgreementStatus::Unknown);
        assert_eq!(zero.checked, 0);
        assert!(matches!(
            zero.unknowns[0].reason,
            AgreementUnknownReason::ZeroSamples
        ));

        let one = with_cx(&gate, |cx| {
            Region::from_chart(Arc::new(sphere()), ProvenanceHash::of_bytes(b"only"))
                .check_agreement(&AgreementConfig::default(), cx)
                .expect("not cancelled")
        });
        assert_eq!(one.status, AgreementStatus::Unknown);
        assert!(matches!(
            one.unknowns[0].reason,
            AgreementUnknownReason::InsufficientCharts { found: 1 }
        ));
    }

    #[test]
    fn non_finite_tolerances_are_unknown_without_sampling() {
        let gate = CancelGate::new();
        for tolerance_abs in [f64::NAN, f64::INFINITY] {
            let report = with_cx(&gate, |cx| {
                let region = Region::from_chart(Arc::new(sphere()), ProvenanceHash::of_bytes(b"a"))
                    .with_chart(Arc::new(sphere()), ProvenanceHash::of_bytes(b"b"));
                region
                    .check_agreement(
                        &AgreementConfig {
                            tolerance_abs,
                            ..AgreementConfig::default()
                        },
                        cx,
                    )
                    .expect("not cancelled")
            });
            assert_eq!(report.status, AgreementStatus::Unknown);
            assert_eq!(report.checked, 0);
            assert!(matches!(
                report.unknowns[0].reason,
                AgreementUnknownReason::InvalidTolerance { .. }
            ));
        }
    }

    #[test]
    fn non_finite_output_and_no_claim_are_unknown() {
        let gate = CancelGate::new();
        let nan_report = with_cx(&gate, |cx| {
            region_with_probe(probe(crate::ChartSample {
                signed_distance: f64::NAN,
                gradient: None,
                lipschitz: Some(1.0),
                error: NumericalCertificate::exact(f64::NAN),
            }))
            .check_agreement(
                &AgreementConfig {
                    samples: 1,
                    ..AgreementConfig::default()
                },
                cx,
            )
            .expect("not cancelled")
        });
        assert_eq!(nan_report.status, AgreementStatus::Unknown);
        assert_eq!(nan_report.checked, 0);
        assert_eq!(nan_report.unknown_count, 1);
        assert!(matches!(
            nan_report.unknowns[0].reason,
            AgreementUnknownReason::NonFiniteSignedDistance { .. }
        ));
        assert!(!nan_report.to_json().contains(":NaN"));

        let no_claim_report = with_cx(&gate, |cx| {
            region_with_probe(probe(crate::ChartSample {
                signed_distance: 0.0,
                gradient: None,
                lipschitz: Some(1.0),
                error: NumericalCertificate::no_claim(),
            }))
            .check_agreement(
                &AgreementConfig {
                    samples: 1,
                    ..AgreementConfig::default()
                },
                cx,
            )
            .expect("not cancelled")
        });
        assert_eq!(no_claim_report.status, AgreementStatus::Unknown);
        assert_eq!(no_claim_report.checked, 0);
        assert_eq!(no_claim_report.unknown_count, 1);
        assert!(matches!(
            no_claim_report.unknowns[0].reason,
            AgreementUnknownReason::NoClaim
        ));
    }

    #[test]
    fn malformed_certificate_and_support_are_unknown() {
        let gate = CancelGate::new();
        let malformed = probe(crate::ChartSample {
            signed_distance: 0.0,
            gradient: None,
            lipschitz: Some(1.0),
            error: NumericalCertificate {
                kind: NumericalKind::Exact,
                lo: -1.0,
                hi: 1.0,
            },
        });
        let report = with_cx(&gate, |cx| {
            region_with_probe(malformed)
                .check_agreement(
                    &AgreementConfig {
                        samples: 1,
                        ..AgreementConfig::default()
                    },
                    cx,
                )
                .expect("not cancelled")
        });
        assert_eq!(report.status, AgreementStatus::Unknown);
        assert!(matches!(
            report.unknowns[0].reason,
            AgreementUnknownReason::MalformedCertificate { .. }
        ));

        let invalid_support = ProbeChart {
            support: crate::Aabb {
                min: Point3::new(f64::NAN, -1.0, -1.0),
                max: Point3::new(1.0, 1.0, 1.0),
            },
            ..malformed
        };
        let report = with_cx(&gate, |cx| {
            region_with_probe(invalid_support)
                .check_agreement(&AgreementConfig::default(), cx)
                .expect("not cancelled")
        });
        assert_eq!(report.status, AgreementStatus::Unknown);
        assert!(matches!(
            report.unknowns[0].reason,
            AgreementUnknownReason::InvalidSupport
        ));
    }

    #[test]
    fn exact_disagreement_survives_a_zero_diagnostic_cap() {
        let region = Region::from_chart(Arc::new(sphere()), ProvenanceHash::of_bytes(b"honest"))
            .with_chart(
                Arc::new(LyingSphereChart {
                    sphere: sphere(),
                    bias: 0.05,
                }),
                ProvenanceHash::of_bytes(b"liar"),
            );
        let gate = CancelGate::new();
        let report = with_cx(&gate, |cx| {
            region
                .check_agreement(
                    &AgreementConfig {
                        samples: 1,
                        max_diagnostics: 0,
                        ..AgreementConfig::default()
                    },
                    cx,
                )
                .expect("not cancelled")
        });
        assert_eq!(report.status, AgreementStatus::Disagreed);
        assert_eq!(report.disagreement_count, 1);
        assert!(report.disagreements.is_empty());
        assert_eq!(
            report.strongest_counterexample_evidence,
            Some(NumericalKind::Exact)
        );
    }

    #[test]
    fn agreement_report_preserves_the_weakest_certificate_class() {
        let exact = probe(crate::ChartSample {
            signed_distance: 0.0,
            gradient: None,
            lipschitz: Some(1.0),
            error: NumericalCertificate::exact(0.0),
        });
        let estimated = ProbeChart {
            name: "estimated-probe",
            sample: crate::ChartSample {
                error: NumericalCertificate::estimate(-0.1, 0.1),
                ..exact.sample
            },
            ..exact
        };
        let region = Region::from_chart(Arc::new(exact), ProvenanceHash::of_bytes(b"exact"))
            .with_chart(Arc::new(estimated), ProvenanceHash::of_bytes(b"estimated"));
        let gate = CancelGate::new();
        let report = with_cx(&gate, |cx| {
            region
                .check_agreement(
                    &AgreementConfig {
                        samples: 1,
                        ..AgreementConfig::default()
                    },
                    cx,
                )
                .expect("not cancelled")
        });
        assert_eq!(report.status, AgreementStatus::Agreed);
        assert_eq!(report.weakest_evidence, Some(NumericalKind::Estimate));
        assert_eq!(report.strongest_counterexample_evidence, None);
    }
}
