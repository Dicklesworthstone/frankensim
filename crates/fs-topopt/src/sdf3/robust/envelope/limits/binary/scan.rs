//! Bounded threshold search; a gray feasible design is never a fallback endpoint.
use super::*;

#[derive(Debug, Clone, Copy)]
pub struct BinaryComplianceScanOptions3 {
    pub robustness: RobustOptions3,
    pub void_stiffness: f64,
    pub max_load_cases: usize,
    pub max_candidates: usize,
    /// Independent upper limit on BINARY DILATED material volume, in [0,1].
    pub dilated_volume_limit: f64,
    /// Numerical tolerance of the smooth aggregate only; no physical cap slack.
    pub smooth_tolerance: f64,
}
impl Default for BinaryComplianceScanOptions3 {
    fn default() -> Self {
        Self { robustness: RobustOptions3::default(), void_stiffness: 1e-3,
            max_load_cases: 64, max_candidates: 64, dilated_volume_limit: 1.0,
            smooth_tolerance: 1e-6 }
    }
}
impl BinaryComplianceScanOptions3 {
    // Shared admission for independent-case and continuous-load scans. This
    // does not choose their feasibility rule; each result owns that decision.
    pub(in crate::sdf3::robust::envelope) fn validate(self, params: SimpParams,
        thresholds: &[f64]) -> Result<(), BinaryComplianceError3> {
        if thresholds.is_empty() || thresholds.len() > self.max_candidates
            || !self.dilated_volume_limit.is_finite() || !(0.0..=1.0).contains(&self.dilated_volume_limit)
            || !self.smooth_tolerance.is_finite() || !(0.0..1.0).contains(&self.smooth_tolerance) {
            return Err(invalid("invalid binary candidate schedule, material limit or smooth tolerance"));
        }
        for &threshold in thresholds { validate_options(params, self.endpoint(threshold))?; }
        Ok(())
    }
    pub(in crate::sdf3::robust::envelope) fn endpoint(self, threshold: f64) -> BinaryComplianceOptions3 {
        BinaryComplianceOptions3 { threshold, robustness: self.robustness,
            void_stiffness: self.void_stiffness, max_load_cases: self.max_load_cases }
    }
}
/// Every threshold at which any of the three binary masks can change.
/// Because material uses `projected >= threshold`, each open interval between
/// consecutive positive projected values has the same masks as its upper end.
/// Adding 1.0 covers the final interval. This enumerates MATERIAL states only;
/// feasibility still requires a fresh physical solve for every threshold.
pub(in crate::sdf3::robust::envelope) fn critical_thresholds(
    source: &[Design3; 3], max_candidates: usize,
) -> Result<Vec<f64>, BinaryComplianceError3> {
    if max_candidates == 0 { return Err(invalid("binary candidate allowance must be positive")); }
    let mut thresholds = Vec::new();
    for design in source {
        for &value in &design.projected {
            if !value.is_finite() || !(0.0..=1.0).contains(&value) {
                return Err(invalid("binary source projection left [0,1]"));
            }
            if value > 0.0 { thresholds.push(value); }
        }
    }
    thresholds.push(1.0);
    thresholds.sort_by(f64::total_cmp);
    thresholds.dedup_by(|a, b| *a == *b);
    if thresholds.len() > max_candidates {
        return Err(invalid("distinct binary material states exceed candidate allowance"));
    }
    Ok(thresholds)
}

#[derive(Debug, Clone, PartialEq)]
pub struct BinaryComplianceCandidate3 {
    pub index: usize,
    pub threshold: f64,
    pub eroded_volume: f64,
    pub nominal_volume: f64,
    pub dilated_volume: f64,
    pub case_compliances: Vec<f64>,
    pub envelope: f64,
    pub maximum_ratio: f64,
    pub physical_caps_satisfied: bool,
    /// Direct physical cap checks, smooth constraint and dilated material cap.
    pub feasible: bool,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BinaryComplianceScanStop3 { ScheduleComplete, EvaluationStopped }
#[derive(Debug, Clone)]
pub struct BinaryComplianceScan3 {
    pub options: BinaryComplianceScanOptions3,
    pub candidates: Vec<BinaryComplianceCandidate3>,
    /// Least dilated material among COMPLETED, feasible endpoints. Equal-volume
    /// ties retain the first schedule entry; the threshold order is authored.
    pub best: Option<BinaryComplianceEvaluation3>,
    pub selected_index: Option<usize>,
    pub stop: BinaryComplianceScanStop3,
    pub error: Option<BinaryComplianceError3>,
    pub work: SolveWork,
}
impl<O: Sdf3Elasticity> CutDensityStudy3<O> {
    /// Threshold one retained three-field projection family. Solve all eroded
    /// loads independently for EVERY candidate; no feasibility monotonicity,
    /// interpolation, gray-field reuse, or optimizer restart is assumed.
    /// Keep only the current candidate and best endpoint's displacement family.
    ///
    /// Entire schedule/policy admission precedes the source filter. All solves
    /// and preparation share control. Late interruption retains the completed
    /// summaries and best feasible endpoint, with stopped status and spent work.
    /// No feasible endpoint means exactly None, not the original gray design or
    /// a proof that all thresholds (or geometries) are infeasible.
    #[allow(clippy::too_many_arguments)]
    pub fn scan_binary_compliance_caps(&mut self, rho: &[f64], forces: &[&[f64]], caps: &ComplianceCaps3,
        thresholds: &[f64], options: BinaryComplianceScanOptions3, control: &mut SolveControl<'_>)
        -> Result<BinaryComplianceScan3, BinaryComplianceError3> {
        options.validate(self.params(), thresholds)?;
        admit(self, rho, forces, caps, options.endpoint(thresholds[0]))?;
        let mut report = BinaryComplianceScan3 { options, candidates: Vec::new(), best: None,
            selected_index: None, stop: BinaryComplianceScanStop3::ScheduleComplete,
            error: None, work: control.work() };
        let attempt = (|| -> Result<(), BinaryComplianceError3> {
            let source = self.binary_cap_source(rho, options.robustness, control)?;
            for (index, &threshold) in thresholds.iter().enumerate() {
                control.checkpoint("sdf3-binary-cap-scan-candidate")?;
                let endpoint = self.binary_caps_from_source(&source, forces, caps, options.endpoint(threshold), control)?;
                let feasible = endpoint.feasible(options.smooth_tolerance)
                    && endpoint.dilated.volume_fraction <= options.dilated_volume_limit;
                report.candidates.push(BinaryComplianceCandidate3 { index, threshold,
                    eroded_volume: endpoint.eroded.volume_fraction,
                    nominal_volume: endpoint.nominal.volume_fraction, dilated_volume: endpoint.dilated.volume_fraction,
                    case_compliances: endpoint.case_compliances.clone(), envelope: endpoint.envelope.value,
                    maximum_ratio: endpoint.envelope.maximum,
                    physical_caps_satisfied: caps.satisfied_by(&endpoint.case_compliances), feasible });
                if feasible && report.best.as_ref().is_none_or(|best|
                    endpoint.dilated.volume_fraction < best.dilated.volume_fraction) {
                    report.best = Some(endpoint); report.selected_index = Some(index);
                }
                // Store paired fields before a late cancellation can escape.
                control.checkpoint("sdf3-binary-cap-scan-accepted")?;
            }
            Ok(())
        })();
        if let Err(error) = attempt { report.error = Some(error); report.stop = BinaryComplianceScanStop3::EvaluationStopped; }
        report.work = control.work();
        Ok(report)
    }

    /// Exhaustively evaluate every distinct three-field binary material state
    /// admitted by the current gray source, up to `max_candidates`. The source
    /// is filtered/projected once; every distinct threshold still receives fresh
    /// binary-eroded load solves. Refusal after source construction leaves spent
    /// work observable through `control`; no sampled subset is substituted.
    #[allow(clippy::too_many_arguments)]
    pub fn scan_all_binary_compliance_caps(&mut self, rho: &[f64], forces: &[&[f64]],
        caps: &ComplianceCaps3, options: BinaryComplianceScanOptions3,
        control: &mut SolveControl<'_>) -> Result<BinaryComplianceScan3, BinaryComplianceError3> {
        options.validate(self.params(), &[1.0])?;
        admit(self, rho, forces, caps, options.endpoint(1.0))?;
        let source = self.binary_cap_source(rho, options.robustness, control)?;
        let thresholds = critical_thresholds(&source, options.max_candidates)?;
        let mut report = BinaryComplianceScan3 { options, candidates: Vec::new(), best: None,
            selected_index: None, stop: BinaryComplianceScanStop3::ScheduleComplete,
            error: None, work: control.work() };
        let attempt = (|| -> Result<(), BinaryComplianceError3> {
            for (index, &threshold) in thresholds.iter().enumerate() {
                control.checkpoint("sdf3-binary-cap-scan-candidate")?;
                let endpoint = self.binary_caps_from_source(&source, forces, caps,
                    options.endpoint(threshold), control)?;
                let feasible = endpoint.feasible(options.smooth_tolerance)
                    && endpoint.dilated.volume_fraction <= options.dilated_volume_limit;
                report.candidates.push(BinaryComplianceCandidate3 { index, threshold,
                    eroded_volume: endpoint.eroded.volume_fraction,
                    nominal_volume: endpoint.nominal.volume_fraction,
                    dilated_volume: endpoint.dilated.volume_fraction,
                    case_compliances: endpoint.case_compliances.clone(),
                    envelope: endpoint.envelope.value, maximum_ratio: endpoint.envelope.maximum,
                    physical_caps_satisfied: caps.satisfied_by(&endpoint.case_compliances), feasible });
                if feasible && report.best.as_ref().is_none_or(|best|
                    endpoint.dilated.volume_fraction < best.dilated.volume_fraction) {
                    report.best = Some(endpoint); report.selected_index = Some(index);
                }
                control.checkpoint("sdf3-binary-cap-scan-accepted")?;
            }
            Ok(())
        })();
        if let Err(error) = attempt {
            report.error = Some(error); report.stop = BinaryComplianceScanStop3::EvaluationStopped;
        }
        report.work = control.work();
        Ok(report)
    }
}

#[cfg(test)]
#[path = "scan_tests.rs"]
mod tests;
