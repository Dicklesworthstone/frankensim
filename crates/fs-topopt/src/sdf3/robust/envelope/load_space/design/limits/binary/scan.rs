//! Select solid/void material against the continuous load set, not its columns.
use super::*;
use crate::sdf3::robust::envelope::limits::binary::admit;
pub use crate::sdf3::robust::envelope::limits::binary::{
    BinaryComplianceScanOptions3 as BinaryLoadSpaceScanOptions3,
    BinaryComplianceScanStop3 as BinaryLoadSpaceScanStop3,
};

#[derive(Debug, Clone, PartialEq)]
pub struct BinaryLoadSpaceCandidate3 {
    pub index: usize,
    pub threshold: f64,
    pub eroded_volume: f64,
    pub nominal_volume: f64,
    pub dilated_volume: f64,
    /// Diagnostics only; passing these separately cannot authorize this row.
    pub basis_compliances: Vec<f64>,
    pub worst_compliance: f64,
    pub maximum: f64,
    pub smooth_maximum: f64,
    pub load_set_feasible: bool,
    /// The full load-space requirement AND the separate dilated material cap.
    pub feasible: bool,
}

#[derive(Debug, Clone)]
pub struct BinaryLoadSpaceScan3 {
    pub cap: LoadSpaceCap3,
    pub options: BinaryLoadSpaceScanOptions3,
    /// Only fully solved and spectrally checked candidates, in authored order.
    pub candidates: Vec<BinaryLoadSpaceCandidate3>,
    /// Least dilated volume among completed feasible entries. Equal-volume
    /// ties keep the first entry, not a later threshold or a gray design.
    pub best: Option<BinaryLoadSpaceEvaluation3>,
    pub selected_index: Option<usize>,
    pub stop: BinaryLoadSpaceScanStop3,
    pub error: Option<BinaryComplianceError3>,
    pub work: SolveWork,
}

impl<O: Sdf3Elasticity> CutDensityStudy3<O> {
    /// Admit the entire schedule before work, filter/project the raw design
    /// ONCE, then independently solve all binary-eroded basis loads and the
    /// full cross-compliance spectrum at EVERY threshold. The original binary
    /// solver and spectral admission own these computations; a finite-case
    /// envelope or the gray maximizing direction is never the selection gate.
    ///
    /// All thresholds share the fixed cap, material model, void stiffness,
    /// source projection and physical work budget. No binary mask is refiltered.
    /// Candidates are not skipped on failure or ranked by assumed monotonicity.
    /// Keep only the current and best displacement families, plus summaries.
    ///
    /// An interrupted scan retains its completed prefix and best endpoint, but
    /// is not ScheduleComplete. No feasible endpoint means None, not a global
    /// infeasibility proof. Incoming nominal state is unchanged on returned
    /// outcomes. This is a fixed-grid positive-void-stiffness material search,
    /// not geometric extraction, a gradient, or a continuum certificate.
    #[allow(clippy::too_many_arguments)]
    pub fn scan_binary_load_space_cap(&mut self, rho: &[f64], basis: &[&[f64]],
        cap: &LoadSpaceCap3, thresholds: &[f64], options: BinaryLoadSpaceScanOptions3,
        control: &mut SolveControl<'_>) -> Result<BinaryLoadSpaceScan3, BinaryComplianceError3> {
        self.require_homogeneous_load_space().map_err(BinaryComplianceError3::Invalid)?;
        options.validate(self.params(), thresholds)?;
        let cases = basis_caps(cap, basis.len())?;
        admit(self, rho, basis, &cases, options.endpoint(thresholds[0]))?;
        let mut report = BinaryLoadSpaceScan3 { cap: cap.clone(), options,
            candidates: Vec::new(), best: None, selected_index: None,
            stop: BinaryLoadSpaceScanStop3::ScheduleComplete, error: None, work: control.work() };
        let attempt = (|| -> Result<(), BinaryComplianceError3> {
            let source = self.binary_cap_source(rho, options.robustness, control)?;
            for (index, &threshold) in thresholds.iter().enumerate() {
                control.checkpoint("sdf3-binary-load-space-scan-candidate")?;
                let physical = self.binary_caps_from_source(&source, basis, &cases,
                    options.endpoint(threshold), control)?;
                let endpoint = self.binary_load_space_from_physical(physical, basis, cap, control)?;
                let load_set_feasible = endpoint.feasible(options.smooth_tolerance);
                let feasible = load_set_feasible
                    && endpoint.physical.dilated.volume_fraction <= options.dilated_volume_limit;
                report.candidates.push(BinaryLoadSpaceCandidate3 { index, threshold,
                    eroded_volume: endpoint.physical.eroded.volume_fraction,
                    nominal_volume: endpoint.physical.nominal.volume_fraction,
                    dilated_volume: endpoint.physical.dilated.volume_fraction,
                    basis_compliances: endpoint.physical.case_compliances.clone(),
                    worst_compliance: endpoint.worst_compliance,
                    maximum: endpoint.spectrum.maximum,
                    smooth_maximum: endpoint.spectrum.smooth_maximum,
                    load_set_feasible, feasible });
                if feasible && report.best.as_ref().is_none_or(|best|
                    endpoint.physical.dilated.volume_fraction < best.physical.dilated.volume_fraction) {
                    report.best = Some(endpoint);
                    report.selected_index = Some(index);
                }
                // Retain the complete paired fields before polling a late stop.
                control.checkpoint("sdf3-binary-load-space-scan-accepted")?;
            }
            Ok(())
        })();
        if let Err(error) = attempt {
            report.error = Some(error);
            report.stop = BinaryLoadSpaceScanStop3::EvaluationStopped;
        }
        report.work = control.work();
        Ok(report)
    }

    /// Exhaustively evaluate every distinct three-field binary material state
    /// generated by this gray source, bounded by `max_candidates`. This closes
    /// the gap where a coarse authored threshold schedule can skip a feasible
    /// lower-material endpoint. The source projection is built once, while each
    /// distinct state still receives fresh full-load-space physics and spectrum.
    #[allow(clippy::too_many_arguments)]
    pub fn scan_all_binary_load_space_cap(&mut self, rho: &[f64], basis: &[&[f64]],
        cap: &LoadSpaceCap3, options: BinaryLoadSpaceScanOptions3,
        control: &mut SolveControl<'_>) -> Result<BinaryLoadSpaceScan3, BinaryComplianceError3> {
        self.require_homogeneous_load_space().map_err(BinaryComplianceError3::Invalid)?;
        options.validate(self.params(), &[1.0])?;
        let cases = basis_caps(cap, basis.len())?;
        admit(self, rho, basis, &cases, options.endpoint(1.0))?;
        let source = self.binary_cap_source(rho, options.robustness, control)?;
        let thresholds = crate::sdf3::robust::envelope::limits::binary::critical_thresholds(
            &source, options.max_candidates)?;
        let mut report = BinaryLoadSpaceScan3 { cap: cap.clone(), options,
            candidates: Vec::new(), best: None, selected_index: None,
            stop: BinaryLoadSpaceScanStop3::ScheduleComplete, error: None, work: control.work() };
        let attempt = (|| -> Result<(), BinaryComplianceError3> {
            for (index, &threshold) in thresholds.iter().enumerate() {
                control.checkpoint("sdf3-binary-load-space-scan-candidate")?;
                let physical = self.binary_caps_from_source(&source, basis, &cases,
                    options.endpoint(threshold), control)?;
                let endpoint = self.binary_load_space_from_physical(physical, basis, cap, control)?;
                let load_set_feasible = endpoint.feasible(options.smooth_tolerance);
                let feasible = load_set_feasible
                    && endpoint.physical.dilated.volume_fraction <= options.dilated_volume_limit;
                report.candidates.push(BinaryLoadSpaceCandidate3 { index, threshold,
                    eroded_volume: endpoint.physical.eroded.volume_fraction,
                    nominal_volume: endpoint.physical.nominal.volume_fraction,
                    dilated_volume: endpoint.physical.dilated.volume_fraction,
                    basis_compliances: endpoint.physical.case_compliances.clone(),
                    worst_compliance: endpoint.worst_compliance,
                    maximum: endpoint.spectrum.maximum,
                    smooth_maximum: endpoint.spectrum.smooth_maximum,
                    load_set_feasible, feasible });
                if feasible && report.best.as_ref().is_none_or(|best|
                    endpoint.physical.dilated.volume_fraction < best.physical.dilated.volume_fraction) {
                    report.best = Some(endpoint); report.selected_index = Some(index);
                }
                control.checkpoint("sdf3-binary-load-space-scan-accepted")?;
            }
            Ok(())
        })();
        if let Err(error) = attempt {
            report.error = Some(error); report.stop = BinaryLoadSpaceScanStop3::EvaluationStopped;
        }
        report.work = control.work();
        Ok(report)
    }
}

#[cfg(test)]
#[path = "scan_tests.rs"]
mod tests;
