//! `Region`: one abstract shape presented through many charts, with
//! agreement as a CHECKABLE PROPOSITION (plan §7.1) — sampled, seeded,
//! deterministic, and localized when it fails.

use crate::{Chart, Point3};
use fs_evidence::ProvenanceHash;
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
    /// union support box. Two charts agree at `x` when their signed
    /// distances differ by no more than the SUM of their declared error
    /// half-widths plus `config.tolerance_abs` — a disagreement therefore
    /// means at least one chart's geometry OR error declaration is wrong.
    /// Deterministic: same seed, same report. Polls cancellation per
    /// sample (Decalogue P7).
    ///
    /// # Errors
    /// [`Cancelled`] when the context's gate was requested mid-check.
    pub fn check_agreement(
        &self,
        config: &AgreementConfig,
        cx: &Cx<'_>,
    ) -> Result<AgreementReport, Cancelled> {
        let mut support = self.charts[0].chart.support();
        for rc in &self.charts[1..] {
            support = support.union(&rc.chart.support());
        }
        // Sample the inflated box so boundary behavior is exercised too.
        let box_ = support.inflate(0.1 * (support.max.x - support.min.x).abs().max(1e-3));
        let mut state = config.seed | 1;
        let mut unit = move || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((state >> 11) as f64) / (1u64 << 53) as f64
        };
        let mut report = AgreementReport {
            checked: 0,
            worst_gap: 0.0,
            agreed: true,
            disagreements: Vec::new(),
        };
        for _ in 0..config.samples {
            cx.checkpoint()?;
            let p = Point3::new(
                box_.min.x + (box_.max.x - box_.min.x) * unit(),
                box_.min.y + (box_.max.y - box_.min.y) * unit(),
                box_.min.z + (box_.max.z - box_.min.z) * unit(),
            );
            for i in 0..self.charts.len() {
                for j in (i + 1)..self.charts.len() {
                    let a = self.charts[i].chart.eval(p, cx);
                    let b = self.charts[j].chart.eval(p, cx);
                    let gap = (a.signed_distance - b.signed_distance).abs();
                    let allowed = half_width(&a.error)
                        + half_width(&b.error)
                        + config.tolerance_abs;
                    report.worst_gap = report.worst_gap.max(gap - allowed);
                    report.checked += 1;
                    if gap > allowed {
                        report.agreed = false;
                        if report.disagreements.len() < config.max_diagnostics {
                            report.disagreements.push(Disagreement {
                                at: p,
                                chart_a: self.charts[i].chart.name(),
                                chart_b: self.charts[j].chart.name(),
                                sd_a: a.signed_distance,
                                sd_b: b.signed_distance,
                                gap,
                                allowed,
                            });
                        }
                    }
                }
            }
        }
        Ok(report)
    }
}

fn half_width(cert: &fs_evidence::NumericalCertificate) -> f64 {
    if cert.hi.is_finite() && cert.lo.is_finite() {
        0.5 * (cert.hi - cert.lo)
    } else {
        f64::INFINITY
    }
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
}

/// The agreement verdict with localized diagnostics.
#[derive(Debug, Clone, PartialEq)]
pub struct AgreementReport {
    /// Pairwise comparisons performed.
    pub checked: u64,
    /// Worst observed `gap - allowed` (≤ 0 when everything agreed).
    pub worst_gap: f64,
    /// True when every comparison stayed within its allowance.
    pub agreed: bool,
    /// First-K localized disagreements.
    pub disagreements: Vec<Disagreement>,
}

impl AgreementReport {
    /// Canonical JSON (deterministic — replayable evidence for the
    /// explain() payload).
    #[must_use]
    pub fn to_json(&self) -> String {
        let mut s = String::with_capacity(128);
        let _ = write!(
            s,
            "{{\"checked\":{},\"agreed\":{},\"worst_gap\":{},\"disagreements\":[",
            self.checked, self.agreed, self.worst_gap
        );
        for (i, d) in self.disagreements.iter().enumerate() {
            if i > 0 {
                s.push(',');
            }
            let _ = write!(
                s,
                "{{\"at\":[{},{},{}],\"charts\":[\"{}\",\"{}\"],\"sd\":[{},{}],\
                 \"gap\":{},\"allowed\":{}}}",
                d.at.x, d.at.y, d.at.z, d.chart_a, d.chart_b, d.sd_a, d.sd_b, d.gap, d.allowed
            );
        }
        s.push_str("]}");
        s
    }
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
        let region = Region::from_chart(
            Arc::new(sphere()),
            ProvenanceHash::of_bytes(b"a"),
        )
        .with_chart(Arc::new(sphere()), ProvenanceHash::of_bytes(b"b"));
        let gate = CancelGate::new();
        let (r1, r2) = with_cx(&gate, |cx| {
            let cfg = AgreementConfig::default();
            (
                region.check_agreement(&cfg, cx).expect("not cancelled"),
                region.check_agreement(&cfg, cx).expect("not cancelled"),
            )
        });
        assert!(r1.agreed && r1.disagreements.is_empty());
        assert_eq!(r1.to_json(), r2.to_json(), "seeded check replays (G5)");
    }

    #[test]
    fn lying_chart_is_detected_with_localized_diagnostics() {
        let region = Region::from_chart(
            Arc::new(sphere()),
            ProvenanceHash::of_bytes(b"honest"),
        )
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
        assert!(!report.agreed, "the bias must be caught");
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
        let region = Region::from_chart(
            Arc::new(sphere()),
            ProvenanceHash::of_bytes(b"a"),
        )
        .with_chart(Arc::new(sphere()), ProvenanceHash::of_bytes(b"b"));
        let gate = CancelGate::new();
        gate.request();
        let outcome = with_cx(&gate, |cx| {
            region.check_agreement(&AgreementConfig::default(), cx)
        });
        assert_eq!(outcome, Err(Cancelled), "pre-requested gate cancels");
    }
}
