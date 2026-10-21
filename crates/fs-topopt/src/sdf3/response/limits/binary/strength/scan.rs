//! Explicit threshold enumeration with original response AND physical-stress gates.
use super::*;

#[derive(Debug, Clone, Copy)]
pub struct BinaryStressScanOptions3 {
    /// Every policy except endpoint.threshold, which the schedule supplies.
    pub stress: BinaryResponseStressOptions3,
    pub max_candidates: usize,
    /// Numerical solid fraction cap, independent of the gray source volume.
    pub volume_cap: f64,
    /// Applied to the joint aggregate AND each original response/stress residual.
    pub feasibility_tolerance: f64,
}
impl Default for BinaryStressScanOptions3 {
    fn default() -> Self {
        Self { stress: BinaryResponseStressOptions3::default(), max_candidates: 64,
            volume_cap: 1.0, feasibility_tolerance: 0.0 }
    }
}
/// Complete candidate summaries; only the incumbent retains displacement fields.
#[derive(Debug, Clone, PartialEq)]
pub struct BinaryStressCandidate3 {
    pub threshold: f64,
    pub volume_fraction: f64,
    pub constraint: f64,
    pub worst_response_residual: f64,
    pub stress_aggregate: f64,
    pub sampled_physical_max: f64,
    pub stress_residual: f64,
    pub response_feasible: bool,
    pub joint_feasible: bool,
    pub volume_feasible: bool,
}
/// Least-volume feasible COMPLETED endpoint among this caller's schedule.
/// Equal-volume ties retain the first. This is not a global binary optimizer.
#[derive(Debug, Clone)]
pub struct BinaryStressScan3 {
    pub candidates: Vec<BinaryStressCandidate3>,
    pub best: Option<BinaryResponseStressEvaluation3>,
    pub best_index: Option<usize>,
    /// A late stop preserves completed evidence, not the partial current trial.
    /// Ok with no best means that this schedule found no joint-feasible endpoint.
    pub outcome: Result<(), ResponseStressError3>,
    /// Cumulative work including rejected, cancelled and incomplete candidates.
    pub work: SolveWork,
}
impl<O: AdaptiveSdf3Elasticity> CutDensityStudy3<O> {
    /// Project the source once, then solve each actual binary experiment family
    /// and test ALL physical limits in caller order. Force windows and stress
    /// need not be monotone with threshold, so no bisection or feasibility
    /// interpolation is used. Duplicate masks are deliberately solved again.
    ///
    /// Policy and the complete schedule are admitted before callbacks or solves.
    /// Later failures retain a completed prefix and its least-volume feasible
    /// field. No mask replaces the source material. A candidate adds only its
    /// normal response preparation/primal solves plus bounded stress sampling.
    /// Sampling has a per-candidate whole-family cap; solve work is cumulative.
    pub fn scan_binary_response_stress_limits(&mut self, rho: &[f64],
        cases: &[ResponseLimitCase3<'_>], thresholds: &[f64], options: BinaryStressScanOptions3,
        control: &mut SolveControl<'_>) -> Result<BinaryStressScan3, ResponseStressError3> {
        control.checkpoint("response-binary-stress-scan-start")?;
        if thresholds.is_empty() || thresholds.len() > options.max_candidates
            || !options.volume_cap.is_finite() || !(0.0..=1.0).contains(&options.volume_cap)
            || !options.feasibility_tolerance.is_finite() || options.feasibility_tolerance < 0.0 {
            return Err(ResponseStressError3::Invalid("invalid binary stress threshold schedule or acceptance policy"));
        }
        for &threshold in thresholds {
            let mut policy = options.stress; policy.endpoint.threshold = threshold; policy.validate()?;
        }
        admit_limits(self, cases, options.stress.endpoint.limits)?;
        if rho.len() != self.cells() || rho.iter().any(|r| !r.is_finite() || !(0.0..=1.0).contains(r)) {
            return Err(ResponseStressError3::Invalid("invalid binary stress source controls"));
        }
        let source = self.design(rho, control)?;
        let mut report = BinaryStressScan3 { candidates: Vec::new(), best: None, best_index: None,
            outcome: Ok(()), work: control.work() };
        report.outcome = (|| {
            for (index, &threshold) in thresholds.iter().enumerate() {
                control.checkpoint("response-binary-stress-scan-candidate")?;
                let mut policy = options.stress; policy.endpoint.threshold = threshold;
                let response = self.binary_from_projected(&source.projected, source.volume,
                    cases, policy.endpoint, control)?;
                let endpoint = self.binary_physical_stress(response, policy, control)?;
                let joint_feasible = endpoint.feasible(options.feasibility_tolerance);
                let volume_feasible = endpoint.response.volume_fraction <= options.volume_cap;
                let summary = BinaryStressCandidate3 { threshold, volume_fraction: endpoint.response.volume_fraction,
                    constraint: endpoint.constraint, worst_response_residual: endpoint.response.worst_residual,
                    stress_aggregate: endpoint.stress_aggregate, sampled_physical_max: endpoint.sampled_physical_max,
                    stress_residual: endpoint.stress_residual,
                    response_feasible: endpoint.response.feasible(options.feasibility_tolerance),
                    joint_feasible, volume_feasible };
                control.checkpoint("response-binary-stress-scan-retain")?;
                report.candidates.push(summary);
                if joint_feasible && volume_feasible && report.best.as_ref()
                    .is_none_or(|best| endpoint.response.volume_fraction < best.response.volume_fraction) {
                    report.best = Some(endpoint); report.best_index = Some(index);
                }
                control.checkpoint("response-binary-stress-scan-accepted")?;
            }
            control.checkpoint("response-binary-stress-scan-publish")?;
            Ok(())
        })();
        report.work = control.work();
        Ok(report)
    }
}
