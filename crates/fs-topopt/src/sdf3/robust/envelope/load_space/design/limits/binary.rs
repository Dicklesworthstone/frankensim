//! Recheck a solid/void endpoint against the ENTIRE simultaneous-load set.
//! The gray maximizing direction and separately passing basis cases are not
//! substitutes for the spectrum of the freshly solved binary material.
mod scan;
pub use scan::{BinaryLoadSpaceCandidate3, BinaryLoadSpaceScan3,
    BinaryLoadSpaceScanOptions3, BinaryLoadSpaceScanStop3};
use super::LoadSpaceCap3;
use super::super::super::{LoadSpaceSpectrum3, SPECTRAL_TOLERANCE, spectrum};
use crate::sdf3::robust::envelope::limits::{ComplianceCaps3,
    binary::BinaryComplianceEvaluation3};
pub use crate::sdf3::robust::envelope::limits::binary::{
    BinaryComplianceOptions3, BinaryComplianceError3,
};
use crate::sdf3::{CutDensityStudy3, Sdf3Elasticity, failure};
use crate::{SolveControl, SolveWork};

/// Complete discrete binary-material assessment, not a geometry certificate.
#[derive(Debug, Clone)]
pub struct BinaryLoadSpaceEvaluation3 {
    pub cap: LoadSpaceCap3,
    /// Three masks and the independent binary-eroded basis solves. Its finite
    /// case envelope is DIAGNOSTIC: use this result's spectrum and feasible().
    pub physical: BinaryComplianceEvaluation3,
    pub spectrum: LoadSpaceSpectrum3,
    pub worst_coefficients: Vec<f64>,
    pub worst_force: Vec<f64>,
    pub worst_displacement: Vec<f64>,
    pub worst_compliance: f64,
    /// Difference between normalized force work and the numerical eigenvalue.
    pub witness_relative_defect: f64,
    pub work: SolveWork,
}
impl BinaryLoadSpaceEvaluation3 {
    /// Shared physical-cap admission with the gray load-space optimizer. The
    /// smoothing tolerance never permits either physical work or the numerical
    /// spectral maximum to exceed the requirement. No volume cap is implied.
    #[must_use]
    pub fn feasible(&self, smooth_tolerance: f64) -> bool {
        self.cap.admits_spectrum(self.cap.policy(), &self.spectrum,
            self.worst_compliance, smooth_tolerance)
    }
}

fn basis_caps(cap: &LoadSpaceCap3, count: usize) -> Result<ComplianceCaps3, BinaryComplianceError3> {
    if count != cap.policy().dimension() {
        return Err(BinaryComplianceError3::Invalid("binary load-space basis count differs from its cap"));
    }
    ComplianceCaps3::new(vec![cap.cap(); count], cap.policy().temperature(), count)
        .map_err(BinaryComplianceError3::Invalid)
}

impl<O: Sdf3Elasticity> CutDensityStudy3<O> {
    /// Filter/project/threshold once using the existing binary material owner,
    /// solve EVERY basis column on the actual binary eroded stiffness, then
    /// compute all cross compliances and the worst simultaneous combination.
    /// Neither a gray displacement nor a previously selected direction is used.
    ///
    /// The same load-space spectrum owner checks reciprocity, positive
    /// semidefiniteness, orthogonality and eigenpair residuals. A direct f^T u
    /// consistency check also precedes publication. No extra elasticity solve,
    /// derivative claim, eigenvalue clipping or artificial stiffness is added.
    ///
    /// Homogeneous supports and fixed dead-load columns only. The declared
    /// positive void stiffness, loads in voids and fixed cut geometry remain
    /// exactly those of evaluate_binary_compliance_caps. All operations share
    /// control; returned stops leave the incoming study unchanged, with no
    /// partially assessed family published. At most 16 columns are admitted
    /// by LoadSpaceCap3. This is not zero-void or continuum certification.
    pub fn evaluate_binary_load_space_cap(&mut self, rho: &[f64], basis: &[&[f64]],
        cap: &LoadSpaceCap3, options: BinaryComplianceOptions3, control: &mut SolveControl<'_>)
        -> Result<BinaryLoadSpaceEvaluation3, BinaryComplianceError3> {
        self.require_homogeneous_load_space().map_err(BinaryComplianceError3::Invalid)?;
        let cases = basis_caps(cap, basis.len())?;
        // These original per-column caps let the existing owner perform the
        // actual binary solves. Their finite-case feasibility is never promoted
        // to a pass for the continuous set; the spectrum below owns that check.
        let physical = self.evaluate_binary_compliance_caps(rho, basis, &cases, options, control)?;
        self.binary_load_space_from_physical(physical, basis, cap, control)
    }

    // Shared by the single endpoint and scanner. Only already-admitted,
    // freshly solved binary fields enter here; there is no public cache ingress.
    fn binary_load_space_from_physical(&self, physical: BinaryComplianceEvaluation3,
        basis: &[&[f64]], cap: &LoadSpaceCap3, control: &mut SolveControl<'_>)
        -> Result<BinaryLoadSpaceEvaluation3, BinaryComplianceError3> {
        self.require_homogeneous_load_space().map_err(BinaryComplianceError3::Invalid)?;
        let policy = cap.policy();
        let n = policy.dimension();
        let states = &physical.displacements;
        let dofs = self.operator().n();
        let fixed = self.operator().fixed();
        let mut gram = vec![0.0; n*n];
        for i in 0..n { for j in 0..n {
            let mut work = 0.0;
            for d in 0..dofs {
                if d % 256 == 0 { control.checkpoint("sdf3-binary-load-space-cross-work")?; }
                if !fixed[d/3] { work += basis[i][d]*states[j][d]; }
            }
            gram[i*n+j] = work / cap.cap();
        }}
        control.checkpoint("sdf3-binary-load-space-spectrum")?;
        let spectrum = spectrum(gram, policy)?;
        control.checkpoint("sdf3-binary-load-space-spectrum")?;
        let worst_coefficients: Vec<_> = (0..n).map(|i| spectrum.eigenvectors[i*n+n-1]).collect();
        let mut worst_force = vec![0.0; dofs];
        let mut worst_displacement = vec![0.0; dofs];
        let mut worst_compliance = 0.0;
        for d in 0..dofs {
            if d % 256 == 0 { control.checkpoint("sdf3-binary-load-space-witness")?; }
            worst_force[d] = (0..n).map(|i| worst_coefficients[i]*basis[i][d]).sum();
            worst_displacement[d] = (0..n).map(|i| worst_coefficients[i]*states[i][d]).sum();
            if !fixed[d/3] { worst_compliance += worst_force[d]*worst_displacement[d]; }
        }
        let normalized_work = worst_compliance / cap.cap();
        let scale = normalized_work.abs().max(spectrum.maximum.abs());
        let witness_relative_defect = if scale == 0.0 { 0.0 }
            else { (normalized_work/scale-spectrum.maximum/scale).abs() };
        if !worst_compliance.is_finite() || worst_compliance < 0.0 || !normalized_work.is_finite()
            || worst_force.iter().chain(&worst_displacement).any(|v| !v.is_finite())
            || !witness_relative_defect.is_finite() || witness_relative_defect > SPECTRAL_TOLERANCE {
            return Err(failure("sdf3-binary-load-space-witness").into());
        }
        control.checkpoint("sdf3-binary-load-space-publish")?;
        Ok(BinaryLoadSpaceEvaluation3 { cap: cap.clone(), physical, spectrum,
            worst_coefficients, worst_force, worst_displacement, worst_compliance,
            witness_relative_defect, work: control.work() })
    }
}

#[cfg(test)]
#[path = "binary_tests.rs"]
mod tests;
