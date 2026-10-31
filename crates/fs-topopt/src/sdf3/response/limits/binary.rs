//! Threshold a projected density field, then re-solve its actual response limits.
//!
//! This is a fixed-grid binary MATERIAL check, not a body-fitted geometry export.
//! A false cell retains the declared positive ersatz stiffness; quadrature,
//! authored dead loads, support patches and observation laws remain. Configured
//! native material body force follows binary occupancy and vanishes in void.
//! Feasibility depends on that stiffness floor and is not a continuum, zero-void,
//! connectedness, minimum-feature-size or manufacturing certificate.
use super::*;
use crate::sdf3::PhysicalRegion3;

mod scan;
pub use scan::{BinaryCandidate3, BinaryScan3, BinaryScanOptions3};
pub mod strength;

/// Independent endpoint policy; changing void stiffness changes the physics.
#[derive(Debug, Clone, Copy)]
pub struct BinaryEndpointOptions3 {
    /// Solid when projected density >= threshold. All ties move together; no
    /// index-based tie splitting to manufacture an exact material volume.
    /// Must lie in (0,1], so exact void/solid physical regions stay meaningful.
    pub threshold: f64,
    /// Reference-material multiplier for void cells, strictly between 0 and 1.
    /// Never silently substitute the gray design's regularization floor.
    pub void_stiffness: f64,
    /// Same original signed inequalities and smooth upper-maximum policy.
    pub limits: ResponseLimitOptions3,
}
impl Default for BinaryEndpointOptions3 {
    fn default() -> Self {
        Self { threshold: 0.5, void_stiffness: 1e-3, limits: ResponseLimitOptions3::default() }
    }
}
impl BinaryEndpointOptions3 {
    fn validate(self) -> Result<(), ResponseError3> {
        if !self.threshold.is_finite() || self.threshold <= 0.0 || self.threshold > 1.0
            || !self.void_stiffness.is_finite() || self.void_stiffness <= 0.0 || self.void_stiffness >= 1.0 {
            return Err(ResponseError3::Invalid("binary threshold must be in (0,1] and void stiffness in (0,1)"));
        }
        Ok(())
    }
}

/// Complete independently solved endpoint, distinct from a differentiable raw
/// design evaluation. No gradient, adjoint, or gray-design feasibility is reused.
#[derive(Debug, Clone)]
pub struct BinaryLimitEvaluation3 {
    pub threshold: f64,
    pub void_stiffness: f64,
    /// Same active-cell order as the source study; protected regions override.
    pub solid: Vec<bool>,
    /// Original projected numerical volume, before thresholding.
    pub source_volume_fraction: f64,
    /// Numerical cut-volume fraction of solid cells; ersatz volume is excluded.
    pub volume_fraction: f64,
    pub constraint: f64,
    pub worst_residual: f64,
    /// Input experiment/limit order; no force or observation is dropped.
    pub measurements: Vec<Vec<f64>>,
    pub residuals: Vec<Vec<f64>>,
    /// Fresh fields on the same constrained DOFs, for the binary material.
    pub displacements: Vec<Vec<f64>>,
    /// Cumulative shared work, not just the solves for this endpoint.
    pub work: SolveWork,
}
impl BinaryLimitEvaluation3 {
    /// Actual stiffness multipliers used for these fields, not densities to
    /// feed back through the gray pipeline's filter and projection.
    #[must_use]
    pub fn stiffness_scales(&self) -> Vec<f64> {
        self.solid.iter().map(|&s| if s { 1.0 } else { self.void_stiffness }).collect()
    }
    /// Both the smooth upper aggregate and EVERY original signed limit must
    /// pass. A tolerance does not turn a numerical endpoint into a certificate.
    #[must_use]
    pub fn feasible(&self, tolerance: f64) -> bool {
        tolerance.is_finite() && tolerance >= 0.0
            && self.constraint.is_finite() && self.constraint <= tolerance
            && self.worst_residual.is_finite() && self.worst_residual <= tolerance
            && !self.residuals.is_empty() && self.residuals.iter().all(|r|
                !r.is_empty() && r.iter().all(|v| v.is_finite() && *v <= tolerance))
    }
}

impl<O: AdaptiveSdf3Elasticity> CutDensityStudy3<O> {
    /// Filter/project the source ONCE, threshold the resulting physical field,
    /// then prepare and solve fresh independent equilibria. The binary field is
    /// NOT refiltered (which would reintroduce gray material). Solid cells have
    /// unit reference stiffness; void cells use options.void_stiffness exactly.
    ///
    /// Prescribed-motion liftings and reaction functionals are recomputed at the
    /// endpoint stiffness through the SAME primal/load owner as the differentiable
    /// evaluator. Authored dead forces and support laws stay in place, including
    /// on weak cells; configured material body forces follow binary occupancy.
    /// No load-removal rule or stiffness floor is guessed.
    /// Protected physical solid/void cells are retained. Both success and returned
    /// errors restore incoming material. The report never installs a binary mask
    /// as raw optimizer controls or claims a valid manufactured/continuum shape.
    pub fn evaluate_binary_response_limits(&mut self, rho: &[f64], cases: &[ResponseLimitCase3<'_>],
        options: BinaryEndpointOptions3, control: &mut SolveControl<'_>)
        -> Result<BinaryLimitEvaluation3, ResponseError3> {
        control.checkpoint("response-binary-start")?;
        options.validate()?;
        admit_limits(self, cases, options.limits)?;
        if rho.len() != self.cells() || rho.iter().any(|r| !r.is_finite() || !(0.0..=1.0).contains(r)) {
            return Err(ResponseError3::Invalid("one source density in [0,1] per active cell required"));
        }
        let source = self.design(rho, control)?;
        self.binary_from_projected(&source.projected, source.volume, cases, options, control)
    }

    fn binary_from_projected(&mut self, projected: &[f64], source_volume: f64,
        cases: &[ResponseLimitCase3<'_>], options: BinaryEndpointOptions3, control: &mut SolveControl<'_>)
        -> Result<BinaryLimitEvaluation3, ResponseError3> {
        control.checkpoint("response-binary-mask")?;
        let mut solid = Vec::with_capacity(self.cells());
        let mut scales = Vec::with_capacity(self.cells());
        let mut volume = 0.0;
        for (i, (&r, &mass)) in projected.iter().zip(&self.mass).enumerate() {
            if i % 256 == 0 { control.checkpoint("response-binary-mask")?; }
            let present = match self.physical_regions.as_ref().map(|r| r[i]) {
                Some(PhysicalRegion3::Solid) => true,
                Some(PhysicalRegion3::Void) => false,
                _ => r >= options.threshold,
            };
            solid.push(present);
            scales.push(if present { 1.0 } else { options.void_stiffness });
            if present { volume = finite(volume+mass)?; }
        }
        let previous = self.operator.scales().to_vec();
        self.operator.set_scales(&scales)?;
        let result = (|| {
            let prepared = self.operator.prepare_elasticity(control)?;
            // Binary physical occupancy has no mass in ersatz void, even
            // though its stiffness is deliberately positive.
            let occupancy: Vec<f64> = solid.iter().map(|&present| if present { 1.0 } else { 0.0 }).collect();
            let physical::LimitPrimals3 { displacements, measurements, residuals, all } =
                self.response_limit_primals(cases, &occupancy, &prepared, control)?;
            let (constraint, worst_residual, _) = aggregate(&all, options.limits.sharpness)?;
            control.checkpoint("response-binary-publish")?;
            Ok(BinaryLimitEvaluation3 { threshold: options.threshold, void_stiffness: options.void_stiffness,
                solid, source_volume_fraction: source_volume, volume_fraction: volume,
                constraint, worst_residual, measurements, residuals, displacements, work: control.work() })
        })();
        self.operator.set_scales(&previous).expect("previous admitted scales remain valid");
        result
    }
}
