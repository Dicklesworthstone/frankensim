//! Re-solved solid/void endpoints of all three robust material realizations.
//!
//! Compliance caps apply to the actual binary ERODED stiffness, not the gray
//! solution. DILATED cut volume is the material objective. This is a fixed-grid
//! ersatz-material model: authored dead loads, supports and cut quadrature
//! remain in voids; an attached material body load follows the binary occupancy.
//! It is not geometry extraction, a zero-void solve or a manufacturing proof.
use super::ComplianceCaps3;
use crate::sdf3::{CutDensityStudy3, Design3, Sdf3Elasticity, checked_solve_preconditioned, failure};
use crate::sdf3::robust::RobustOptions3;
use crate::sdf3::robust::envelope::EnvelopeValue3;
use crate::{EvaluationStop, SimpParams, SolveControl, SolveWork};
mod scan;
pub use scan::{BinaryComplianceCandidate3, BinaryComplianceScan3, BinaryComplianceScanOptions3, BinaryComplianceScanStop3};
pub(in crate::sdf3::robust::envelope) use scan::critical_thresholds;


#[derive(Debug, Clone, Copy)]
pub struct BinaryComplianceOptions3 {
    pub robustness: RobustOptions3,
    /// Solid iff the corresponding projected physical density >= threshold.
    /// The SAME threshold is used for all realizations; ties are never split.
    pub threshold: f64,
    /// Explicit reference-stiffness multiplier for every void cell, in (0,1).
    /// This need not equal the gray model's e_min; changing it changes physics.
    pub void_stiffness: f64,
    pub max_load_cases: usize,
}
impl Default for BinaryComplianceOptions3 {
    fn default() -> Self {
        Self { robustness: RobustOptions3::default(), threshold: 0.5,
            void_stiffness: 1e-3, max_load_cases: 64 }
    }
}
#[derive(Debug, Clone, PartialEq)]
pub enum BinaryComplianceError3 {
    Invalid(&'static str),
    Evaluation(EvaluationStop),
}
impl From<EvaluationStop> for BinaryComplianceError3 {
    fn from(e: EvaluationStop) -> Self { Self::Evaluation(e) }
}
impl std::fmt::Display for BinaryComplianceError3 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "binary compliance check stopped: {self:?}")
    }
}
impl std::error::Error for BinaryComplianceError3 {}

/// In the source study's active-cell order. Numerical cut volumes are used;
/// the positive void stiffness does not count as solid material volume.
#[derive(Debug, Clone, PartialEq)]
pub struct BinaryMaterial3 {
    pub solid: Vec<bool>,
    pub volume_fraction: f64,
    pub source_volume_fraction: f64,
}
#[derive(Debug, Clone)]
pub struct BinaryComplianceEvaluation3 {
    /// Gray model which produced the three projected fields, not binary SIMP.
    pub params: SimpParams,
    pub options: BinaryComplianceOptions3,
    pub caps: ComplianceCaps3,
    pub eroded: BinaryMaterial3,
    pub nominal: BinaryMaterial3,
    pub dilated: BinaryMaterial3,
    /// Fresh binary-eroded solves, in cap/load order. Never averaged or dropped.
    pub case_compliances: Vec<f64>,
    pub displacements: Vec<Vec<f64>>,
    pub envelope: EnvelopeValue3,
    pub work: SolveWork,
}
impl BinaryComplianceEvaluation3 {
    /// Exact multipliers used for the returned eroded displacement fields.
    /// These are NOT raw densities to pass back through the gray filter.
    #[must_use]
    pub fn eroded_stiffness_scales(&self) -> Vec<f64> {
        self.eroded.solid.iter().map(|&s| if s { 1.0 } else { self.options.void_stiffness }).collect()
    }
    /// Every physical cap must pass directly, without numerical slack. The
    /// separate tolerance belongs only to the conservative smooth aggregate.
    #[must_use]
    pub fn feasible(&self, smooth_tolerance: f64) -> bool {
        smooth_tolerance.is_finite() && (0.0..1.0).contains(&smooth_tolerance)
            && self.envelope.value.is_finite() && self.envelope.value - 1.0 <= smooth_tolerance
            && self.caps.satisfied_by(&self.case_compliances)
    }
}
fn invalid(message: &'static str) -> BinaryComplianceError3 { BinaryComplianceError3::Invalid(message) }
fn validate_options(params: SimpParams, options: BinaryComplianceOptions3) -> Result<(), BinaryComplianceError3> {
    if !options.robustness.valid(params)
        || !options.threshold.is_finite() || options.threshold <= 0.0 || options.threshold > 1.0
        || !options.void_stiffness.is_finite() || options.void_stiffness <= 0.0 || options.void_stiffness >= 1.0
        || options.max_load_cases == 0 {
        return Err(invalid("invalid binary threshold, stiffness, robustness or load allowance"));
    }
    Ok(())
}
pub(in crate::sdf3::robust::envelope) fn admit<O: Sdf3Elasticity>(study: &CutDensityStudy3<O>, rho: &[f64], forces: &[&[f64]],
    caps: &ComplianceCaps3, options: BinaryComplianceOptions3) -> Result<(), BinaryComplianceError3> {
    validate_options(study.params(), options)?;
    if forces.len() != caps.caps().len() || forces.len() > options.max_load_cases
        || forces.iter().any(|f| f.len() != study.operator().n() || f.iter().any(|v| !v.is_finite()))
        || rho.len() != study.cells() || rho.iter().any(|r| !r.is_finite() || !(0.0..=1.0).contains(r)) {
        return Err(invalid("one finite load per cap and raw density per cell required"));
    }
    Ok(())
}
fn material(source: &Design3, mass: &[f64], threshold: f64, control: &mut SolveControl<'_>)
    -> Result<BinaryMaterial3, EvaluationStop> {
    let mut solid = Vec::with_capacity(mass.len());
    let mut volume_fraction = 0.0;
    for (i, (&r, &m)) in source.projected.iter().zip(mass).enumerate() {
        if i % 256 == 0 { control.checkpoint("sdf3-binary-cap-mask")?; }
        // Physical overrides were applied by the original projection owner:
        // exact 0 stays void and exact 1 stays solid for every admitted threshold.
        let present = r >= threshold;
        solid.push(present);
        if present { volume_fraction += m; }
    }
    if !volume_fraction.is_finite() { return Err(failure("sdf3-binary-cap-volume")); }
    Ok(BinaryMaterial3 { solid, volume_fraction, source_volume_fraction: source.volume })
}

impl<O: Sdf3Elasticity> CutDensityStudy3<O> {
    /// Filter once and threshold all three physical projections, then solve
    /// EVERY load on binary eroded material with the explicit void stiffness.
    /// No binary mask is refiltered and no gray displacement is reused. The
    /// existing prepared solver and true-residual gate perform all state solves.
    /// There is no adjoint or derivative claim for this discontinuous mapping.
    ///
    /// Homogeneous supports and fixed authored force vectors are required.
    /// An attached material body load follows the exact binary eroded mask;
    /// void stiffness contributes no body force. Authored dead loads and
    /// supports on void cells remain in place. The source nominal scales
    /// and parameters are preserved on success and every returned error. All
    /// filtering, preparation and load solves share the supplied control.
    pub fn evaluate_binary_compliance_caps(&mut self, rho: &[f64], forces: &[&[f64]],
        caps: &ComplianceCaps3, options: BinaryComplianceOptions3, control: &mut SolveControl<'_>)
        -> Result<BinaryComplianceEvaluation3, BinaryComplianceError3> {
        admit(self, rho, forces, caps, options)?;
        let source = self.binary_cap_source(rho, options.robustness, control)?;
        self.binary_caps_from_source(&source, forces, caps, options, control)
    }

    pub(in crate::sdf3::robust::envelope) fn binary_cap_source(&self, rho: &[f64], robustness: RobustOptions3, control: &mut SolveControl<'_>)
        -> Result<[Design3; 3], EvaluationStop> {
        control.checkpoint("sdf3-binary-cap-start")?;
        let filtered = self.filter_apply(rho, control)?;
        Ok([
            self.design_from_filtered(&filtered, self.params.eta + robustness.threshold_offset, control)?,
            self.design_from_filtered(&filtered, self.params.eta, control)?,
            self.design_from_filtered(&filtered, self.params.eta - robustness.threshold_offset, control)?,
        ])
    }

    pub(in crate::sdf3::robust::envelope) fn binary_caps_from_source(&mut self, source: &[Design3; 3], forces: &[&[f64]],
        caps: &ComplianceCaps3, options: BinaryComplianceOptions3, control: &mut SolveControl<'_>)
        -> Result<BinaryComplianceEvaluation3, BinaryComplianceError3> {
        let eroded = material(&source[0], &self.mass, options.threshold, control)?;
        let nominal = material(&source[1], &self.mass, options.threshold, control)?;
        let dilated = material(&source[2], &self.mass, options.threshold, control)?;
        // Refuse an ordering defect rather than changing any physical mask.
        if eroded.solid.iter().zip(&nominal.solid).zip(&dilated.solid)
            .any(|((&e, &n), &d)| (e && !n) || (n && !d)) {
            return Err(failure("sdf3-binary-cap-ordering").into());
        }
        let mut scales = Vec::with_capacity(self.cells());
        for (i, &s) in eroded.solid.iter().enumerate() {
            if i % 256 == 0 { control.checkpoint("sdf3-binary-cap-scales")?; }
            scales.push(if s { 1.0 } else { options.void_stiffness });
        }
        let body = if self.has_material_body_load() {
            let mut occupancy = Vec::with_capacity(self.cells());
            for (i, &solid) in eroded.solid.iter().enumerate() {
                if i % 256 == 0 { control.checkpoint("sdf3-binary-cap-body-load")?; }
                occupancy.push(if solid { 1.0 } else { 0.0 });
            }
            self.material_body_force(&occupancy, control)?
        } else { None };
        let previous = self.operator.scales().to_vec();
        self.operator.set_scales(&scales).map_err(|_| failure("sdf3-binary-cap-scales"))?;
        let result = (|| {
            let prepared = self.operator.prepare_elasticity(control)?;
            let mut case_compliances = Vec::with_capacity(forces.len());
            let mut displacements = Vec::with_capacity(forces.len());
            for force in forces {
                control.checkpoint("sdf3-binary-cap-elasticity")?;
                let rhs = self.external_force(force, body.as_deref(), control)?;
                let u = checked_solve_preconditioned(&self.operator, &prepared, &rhs,
                    1e-11, "sdf3-binary-cap-elasticity", control)?;
                let c: f64 = rhs.iter().zip(&u).map(|(f, u)| f * u).sum();
                if !c.is_finite() || c < 0.0 { return Err(failure("sdf3-binary-cap-compliance").into()); }
                case_compliances.push(c); displacements.push(u);
            }
            let envelope = caps.envelope().aggregate(&case_compliances)?;
            control.checkpoint("sdf3-binary-cap-publish")?;
            Ok(BinaryComplianceEvaluation3 { params: self.params, options, caps: caps.clone(),
                eroded, nominal, dilated, case_compliances, displacements, envelope, work: control.work() })
        })();
        self.operator.set_scales(&previous).expect("incoming scales remain admitted");
        result
    }
}

#[cfg(test)]
#[path = "binary/tests.rs"]
mod tests;
