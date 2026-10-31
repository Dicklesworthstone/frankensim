//! Simultaneous load uncertainty: f = sum_i q_i f_i, ||q||_2 <= 1.
//!
//! The largest eigenvalue of F^T K_eroded^-1 F is the worst compliance over
//! this entire load set, not just its basis cases. The smooth spectral maximum
//! has an exact discrete derivative even at repeated eigenvalues. All units,
//! amplitudes and correlations are authored through the columns of F; these
//! columns are NOT probabilities or independent mutually exclusive scenarios.
pub mod design;
mod audit;
use super::{RobustEvaluation3, RobustOptions3};
use crate::gradient_check::{audit_design_map, GradientCheckOptions, GradientSample, MultiLoadGradientCheck};
use crate::pipeline::LoadCase;
use crate::sdf3::{failure, CutDensityStudy3, Sdf3Elasticity};
use crate::{EvaluationStop, SolveControl, SolveWork};

/// Bound the existing dense eigensolver's non-interruptible sweep. Larger load
/// spaces need a separately budgeted spectral implementation, not this path.
pub const MAX_LOAD_SPACE_DIMENSION: usize = 16;
const SPECTRAL_TOLERANCE: f64 = 1e-9;

#[derive(Debug, Clone, PartialEq)]
pub struct LoadSpace3 {
    dimension: usize,
    reference: f64,
    temperature: f64,
}
impl LoadSpace3 {
    /// `reference` is one fixed positive compliance scale for the whole set.
    /// Geometry, basis columns, reference and temperature stay fixed in a run.
    pub fn new(dimension: usize, reference: f64, temperature: f64) -> Result<Self, &'static str> {
        if dimension == 0 || dimension > MAX_LOAD_SPACE_DIMENSION
            || !reference.is_finite() || reference <= 0.0 || !(1.0 / reference).is_finite()
            || !temperature.is_finite() || temperature <= 0.0
            || !(temperature * fs_math::det::ln(dimension as f64)).is_finite() {
            return Err("invalid load-space dimension, reference or smoothing temperature");
        }
        Ok(Self { dimension, reference, temperature })
    }
    #[must_use] pub fn dimension(&self) -> usize { self.dimension }
    #[must_use] pub fn reference(&self) -> f64 { self.reference }
    #[must_use] pub fn temperature(&self) -> f64 { self.temperature }
    /// Algebraic bias in exact arithmetic, NOT a floating-point enclosure.
    #[must_use] pub fn bias_bound(&self) -> f64 { self.temperature * fs_math::det::ln(self.dimension as f64) }
}

#[derive(Debug, Clone)]
pub struct LoadSpaceSpectrum3 {
    /// Symmetric, row-major normalized compliance Gram matrix. Reciprocity is
    /// checked BEFORE averaging the two numerically solved cross terms.
    pub gram: Vec<f64>,
    pub eigenvalues: Vec<f64>,
    /// Columns in row-major order; ascending eigenvalue order.
    pub eigenvectors: Vec<f64>,
    pub weights: Vec<f64>,
    pub maximum: f64,
    pub smooth_maximum: f64,
    pub reciprocity_defect: f64,
    pub eigenpair_residual: f64,
    pub orthogonality_defect: f64,
}

/// Admit the numerical spectrum; no eigenvalue clipping or diagonal loading.
fn spectrum(mut gram: Vec<f64>, policy: &LoadSpace3) -> Result<LoadSpaceSpectrum3, EvaluationStop> {
    let n = policy.dimension;
    if gram.len() != n*n || gram.iter().any(|v| !v.is_finite()) { return Err(failure("sdf3-load-space-gram")); }
    let scale = gram.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
    let mut reciprocity_defect = 0.0_f64;
    for i in 0..n { for j in i+1..n {
        let (a, b) = (gram[i*n+j], gram[j*n+i]);
        let defect = if scale == 0.0 { 0.0 } else { (a/scale - b/scale).abs() };
        reciprocity_defect = reciprocity_defect.max(defect);
        gram[i*n+j] = 0.5*a + 0.5*b;
        gram[j*n+i] = gram[i*n+j];
    }}
    if reciprocity_defect > SPECTRAL_TOLERANCE { return Err(failure("sdf3-load-space-reciprocity")); }
    let normalized: Vec<_> = gram.iter().map(|g| if scale == 0.0 { 0.0 } else { g/scale }).collect();
    let (values, mut vectors) = fs_la::eigen::jacobi_eigh(&normalized, n);
    if values.iter().chain(&vectors).any(|v| !v.is_finite())
        || values[0] < -SPECTRAL_TOLERANCE || values[n-1] < 0.0 {
        return Err(failure("sdf3-load-space-spectrum"));
    }
    // Fix each eigenvector's sign by its first largest-magnitude coefficient.
    // Within a repeated eigenspace the worst direction need not be unique.
    for k in 0..n {
        let mut pivot = 0;
        for i in 1..n { if vectors[i*n+k].abs() > vectors[pivot*n+k].abs() { pivot = i; } }
        if vectors[pivot*n+k] < 0.0 { for i in 0..n { vectors[i*n+k] = -vectors[i*n+k]; } }
    }
    let mut eigenpair_residual = 0.0_f64;
    let mut orthogonality_defect = 0.0_f64;
    for i in 0..n { for k in 0..n {
        let av: f64 = (0..n).map(|j| normalized[i*n+j]*vectors[j*n+k]).sum();
        eigenpair_residual = eigenpair_residual.max((av - values[k]*vectors[i*n+k]).abs());
        let vv: f64 = (0..n).map(|j| vectors[j*n+i]*vectors[j*n+k]).sum();
        orthogonality_defect = orthogonality_defect.max((vv - if i == k { 1.0 } else { 0.0 }).abs());
    }}
    if eigenpair_residual > SPECTRAL_TOLERANCE || orthogonality_defect > SPECTRAL_TOLERANCE {
        return Err(failure("sdf3-load-space-eigenpair"));
    }
    let eigenvalues: Vec<_> = values.iter().map(|v| v*scale).collect();
    let maximum = eigenvalues[n-1];
    let mut weights: Vec<_> = eigenvalues.iter()
        .map(|v| fs_math::det::exp((v-maximum)/policy.temperature)).collect();
    let total: f64 = weights.iter().sum();
    let smooth_maximum = maximum + policy.temperature*fs_math::det::ln(total);
    if eigenvalues.iter().any(|v| !v.is_finite()) || !total.is_finite() || total < 1.0 || !smooth_maximum.is_finite() {
        return Err(failure("sdf3-load-space-range"));
    }
    for w in &mut weights { *w /= total; }
    Ok(LoadSpaceSpectrum3 { gram, eigenvalues, eigenvectors: vectors, weights, maximum,
        smooth_maximum, reciprocity_defect, eigenpair_residual, orthogonality_defect })
}

#[derive(Debug, Clone)]
pub struct LoadSpaceEvaluation3 {
    pub policy: LoadSpace3,
    /// Independent basis solutions and three material fields. Its embedded
    /// mean compliance/gradient is a DIAGNOSTIC, not the spectral objective.
    pub physical: RobustEvaluation3,
    pub spectrum: LoadSpaceSpectrum3,
    /// Derivative of smooth_maximum through the eroded SIMP and filter map.
    pub gradient: Vec<f64>,
    /// Unit-norm coefficients; an eigenvalue tie gives one maximizing direction.
    pub worst_coefficients: Vec<f64>,
    pub worst_force: Vec<f64>,
    /// Linear combination of solved basis fields, not an extra linear solve.
    pub worst_displacement: Vec<f64>,
    /// Actual f^T u of that combination in physical compliance units.
    pub worst_compliance: f64,
    pub work: SolveWork,
}

impl<O: Sdf3Elasticity> CutDensityStudy3<O> {
    // A common material force is an affine offset, not another basis column.
    // The existing Gram eigenproblem solves only a zero-centered load set.
    fn require_homogeneous_load_space(&self) -> Result<(), &'static str> {
        if self.has_material_body_load() {
            return Err("continuous load-space envelopes do not support a material body-load offset");
        }
        Ok(())
    }

    /// Solve each basis column once on eroded material and retain ALL cross
    /// compliances. Independent maxima/averages miss adverse combinations.
    /// The exact maximum is reported alongside a differentiable spectral
    /// log-sum-exp; the latter, not a selected eigenvector, drives gradients.
    ///
    /// Homogeneous supports, density-independent linear dead loads and fixed
    /// geometry only. A zero-centered load ellipsoid is defined by F, not a
    /// distribution; no nonzero nominal-load offset is implied. Singular or
    /// dependent basis columns are allowed, but adding a duplicate CHANGES
    /// the set of attainable forces. All returned outcomes preserve scales.
    pub fn evaluate_robust_load_space(&mut self, rho: &[f64], basis: &[&[f64]],
        robustness: RobustOptions3, policy: &LoadSpace3, control: &mut SolveControl<'_>)
        -> Result<LoadSpaceEvaluation3, EvaluationStop> {
        self.require_homogeneous_load_space().map_err(|_| failure("sdf3-load-space-material-body-load"))?;
        let m = policy.dimension;
        let dofs = self.operator.n();
        if basis.len() != m || basis.iter().any(|f| f.len() != dofs || f.iter().any(|v| !v.is_finite())) {
            return Err(failure("sdf3-load-space-basis"));
        }
        control.checkpoint("sdf3-load-space-start")?;
        let loads: Vec<_> = basis.iter().map(|&force| LoadCase { force, weight: 1.0/m as f64 }).collect();
        let physical = self.evaluate_robust(rho, &loads, robustness, control)?;
        let states = &physical.eroded.objective.displacements;
        let mut gram = vec![0.0; m*m];
        for i in 0..m { for j in 0..m {
            control.checkpoint("sdf3-load-space-cross-work")?;
            let mut work = 0.0;
            for d in 0..dofs {
                if d % 256 == 0 { control.checkpoint("sdf3-load-space-cross-work")?; }
                if !self.operator.fixed()[d/3] { work += basis[i][d]*states[j][d]; }
            }
            gram[i*m+j] = work/policy.reference;
        }}
        control.checkpoint("sdf3-load-space-spectrum")?;
        let spectrum = spectrum(gram, policy)?;
        control.checkpoint("sdf3-load-space-spectrum")?;
        let filtered = self.filter_apply(rho, control)?;
        let design = self.design_from_filtered(&filtered, self.params.eta + robustness.threshold_offset, control)?;
        let mut local = vec![0.0; self.cells()];
        let mut combined = vec![0.0; dofs];
        for k in 0..m {
            for d in 0..dofs {
                if d % 256 == 0 { control.checkpoint("sdf3-load-space-adjoint")?; }
                combined[d] = (0..m).map(|i| spectrum.eigenvectors[i*m+k]*states[i][d]).sum();
            }
            let energies = self.operator.scale_quadratic_forms(&combined).map_err(|_| failure("sdf3-load-space-energy"))?;
            for (g, e) in local.iter_mut().zip(energies) { *g -= (e/policy.reference)*spectrum.weights[k]; }
        }
        // k=m-1 leaves combined at the maximizing displacement.
        let worst_coefficients: Vec<_> = (0..m).map(|i| spectrum.eigenvectors[i*m+m-1]).collect();
        let mut worst_force = vec![0.0; dofs];
        let mut worst_compliance = 0.0;
        for d in 0..dofs {
            if d % 256 == 0 { control.checkpoint("sdf3-load-space-witness")?; }
            worst_force[d] = (0..m).map(|i| worst_coefficients[i]*basis[i][d]).sum();
            if !self.operator.fixed()[d/3] { worst_compliance += worst_force[d]*combined[d]; }
        }
        let p = self.params;
        for ((g, &r), &slope) in local.iter_mut().zip(&design.projected).zip(&design.slope) {
            let power = if p.penal == 1.0 { 1.0 } else if r == 0.0 { 0.0 } else { fs_math::det::pow(r, p.penal-1.0) };
            *g *= (1.0-p.e_min)*p.penal*power*slope;
        }
        if local.iter().chain(&combined).chain(&worst_force).any(|v| !v.is_finite())
            || !worst_compliance.is_finite() || worst_compliance < 0.0 {
            return Err(failure("sdf3-load-space-gradient"));
        }
        let gradient = self.pullback(&local, control)?;
        control.checkpoint("sdf3-load-space-publish")?;
        Ok(LoadSpaceEvaluation3 { policy: policy.clone(), physical, spectrum, gradient,
            worst_coefficients, worst_force, worst_displacement: combined, worst_compliance, work: control.work() })
    }
}

#[derive(Debug, Clone)]
pub struct LoadSpaceGradientCheck3 {
    pub policy: LoadSpace3,
    pub robustness: RobustOptions3,
    /// Compliance quantities denote the smooth normalized spectral objective.
    pub numerical: MultiLoadGradientCheck,
}
pub fn controlled_load_space_gradient_check<O: Sdf3Elasticity>(
    study: &mut CutDensityStudy3<O>, rho: &[f64], basis: &[&[f64]], robustness: RobustOptions3,
    policy: &LoadSpace3, options: GradientCheckOptions, control: &mut SolveControl<'_>,
) -> Result<LoadSpaceGradientCheck3, EvaluationStop> {
    study.require_homogeneous_load_space().map_err(|_| failure("sdf3-load-space-material-body-load"))?;
    options.assert_valid();
    robustness.assert_valid(study.params());
    let mut baseline_compliance = 0.0;
    // Remove only the DESIGN-INDEPENDENT temperature*ln(dimension) offset
    // before differencing. Adding it first can erase small but nonzero load
    // sensitivities under loose caps. The analytic gradient and gate stay exact.
    let mut numerical = audit_design_map(study.params(), rho, options, control, |point, baseline, control| {
        let e = study.evaluate_robust_load_space(point, basis, robustness, policy, control)?;
        if baseline { baseline_compliance = e.spectrum.smooth_maximum; }
        Ok(GradientSample { compliance: audit::centered_maximum(&e.spectrum, policy)?, compliance_gradient: e.gradient,
            volume: e.physical.dilated_volume_fraction, volume_gradient: e.physical.dilated_volume_gradient })
    })?;
    // Keep the public baseline in the original objective's units and convention.
    // All directional differences are invariant to the constant removed above.
    numerical.baseline_compliance = baseline_compliance;
    Ok(LoadSpaceGradientCheck3 { policy: policy.clone(), robustness, numerical })
}

#[cfg(test)]
#[path = "load_space/tests.rs"]
mod tests;
