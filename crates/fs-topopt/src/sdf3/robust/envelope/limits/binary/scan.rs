//! Bounded threshold search: gray feasibility is never a fallback endpoint.
use super::*;

#[derive(Debug, Clone, Copy)]
pub struct BinaryComplianceScanOptions3 {
    pub robustness: RobustOptions3,
    pub void_stiffness: f64,
    pub max_load_cases: usize,
    pub max_candidates: usize,
    pub dilated_volume_limit: f64,
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
    fn endpoint(self, threshold: f64) -> BinaryComplianceOptions3 {
        BinaryComplianceOptions3 { threshold, robustness: self.robustness,
            void_stiffness: self.void_stiffness, max_load_cases: self.max_load_cases }
    }
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
    pub feasible: bool,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BinaryComplianceScanStop3 { ScheduleComplete, EvaluationStopped }
#[derive(Debug, Clone)]
pub struct BinaryComplianceScan3 {
    pub options: BinaryComplianceScanOptions3,
    pub candidates: Vec<BinaryComplianceCandidate3>,
    pub best: Option<BinaryComplianceEvaluation3>,
    pub selected_index: Option<usize>,
    pub stop: BinaryComplianceScanStop3,
    pub error: Option<BinaryComplianceError3>,
    pub work: SolveWork,
}
impl<O: Sdf3Elasticity> CutDensityStudy3<O> {
    /// Filter once, then re-solve each independent load at each actual binary
    /// eroded stiffness. No monotonicity, gray-field reuse or interpolation.
    /// Late cancellation retains completed candidates and the best endpoint.
    #[allow(clippy::too_many_arguments)]
    pub fn scan_binary_compliance_caps(&mut self, rho: &[f64], forces: &[&[f64]], caps: &ComplianceCaps3,
        thresholds: &[f64], options: BinaryComplianceScanOptions3, control: &mut SolveControl<'_>)
        -> Result<BinaryComplianceScan3, BinaryComplianceError3> {
        if thresholds.is_empty() || thresholds.len() > options.max_candidates
            || !options.dilated_volume_limit.is_finite() || !(0.0..=1.0).contains(&options.dilated_volume_limit)
            || !options.smooth_tolerance.is_finite() || !(0.0..1.0).contains(&options.smooth_tolerance) {
            return Err(invalid("invalid binary candidate schedule, material limit or smooth tolerance"));
        }
        for &threshold in thresholds { validate_options(self.params(), options.endpoint(threshold))?; }
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
                control.checkpoint("sdf3-binary-cap-scan-accepted")?;
            }
            Ok(())
        })();
        if let Err(error) = attempt {
            report.error = Some(error);
            report.stop = BinaryComplianceScanStop3::EvaluationStopped;
        }
        report.work = control.work();
        Ok(report)
    }
}
