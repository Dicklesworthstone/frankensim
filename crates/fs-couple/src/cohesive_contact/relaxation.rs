//! Intact relaxing bonds on the original three-component mechanical load line.
//! Reuses its reference area, full body response, prepared LU and physical gates.
//! The affine mean-traction equation needs one solve and one material trial.
//! Actual mechanics still own free motion, accepted time and the complete work
//! window; a straight ramp is not an exact continuous inertial trajectory.
use super::{CohesiveContactError, CohesiveLoadLine, checked, norm_inf};
use fs_material::cohesive::relaxation::{
    ViscoelasticBondPoint, ViscoelasticBondRamp, ViscoelasticBondStep,
};

#[derive(Debug, Clone)]
pub struct ViscoelasticBondLoadStep<'a> {
    pub material: ViscoelasticBondStep<'a>,
    pub mean_force_n: [f64; 3],
    pub compatibility_residual_m: [f64; 3],
    pub compatibility_tolerance_m: f64,
}
impl CohesiveLoadLine {
    /// Solve g=g_free-S*A*mean_traction(g) with exact hereditary means and
    /// the same original mechanical compliance. Neither endpoint force nor a
    /// guessed damping coefficient replaces this work-conjugate traction.
    /// Cancellation/refusal leaves the caller's accepted Prony histories intact.
    pub fn trial_relaxation<'a>(
        &mut self,
        point: &ViscoelasticBondPoint<'a>,
        prepared: &ViscoelasticBondRamp<'_>,
        free_jump: [f64; 3],
        mut keep_going: impl FnMut() -> bool,
    ) -> Result<ViscoelasticBondLoadStep<'a>, CohesiveContactError> {
        if !keep_going() {
            return Err(CohesiveContactError::Cancelled);
        }
        if free_jump.iter().any(|v| !v.is_finite()) {
            return Err(CohesiveContactError::Invalid(
                "finite free bond jump required",
            ));
        }
        let old = point.jump_m();
        let held = point.held_mean_traction(prepared)?;
        let stiffness = prepared.mean_stiffness_pa_m();
        let mut matrix = [0.0; 9];
        let mut rhs = [0.0; 3];
        for i in 0..3 {
            rhs[i] = checked(free_jump[i] - old[i])?;
            for j in 0..3 {
                let response = checked(self.compliance_m_n[i][j] * self.area_m2)?;
                matrix[i * 3 + j] =
                    checked((if i == j { 1.0 } else { 0.0 }) + response * stiffness[j])?;
                rhs[i] = checked(rhs[i] - response * held[j])?;
            }
        }
        let mut increment = [0.0; 3];
        self.solve.solve_into(&matrix, &rhs, &mut increment)?;
        if !keep_going() {
            return Err(CohesiveContactError::Cancelled);
        }
        if norm_inf(increment) > self.config.maximum_jump_increment_m {
            return Err(CohesiveContactError::Range(
                "per-step viscoelastic bond jump",
            ));
        }
        let jump = core::array::from_fn(|i| old[i] + increment[i]);
        let material = point.trial(prepared, jump, self.config.material_energy_tolerance_j_m2)?;
        let force = material.mean_traction_pa.map(|v| v * self.area_m2);
        if force.iter().any(|v| !v.is_finite()) || norm_inf(force) > self.config.maximum_force_n {
            return Err(CohesiveContactError::Range("viscoelastic bond held force"));
        }
        let response = core::array::from_fn::<_, 3, _>(|i| {
            (0..3)
                .map(|j| self.compliance_m_n[i][j] * force[j])
                .sum::<f64>()
        });
        let mut residual = [0.0; 3];
        for i in 0..3 {
            residual[i] = checked(jump[i] - free_jump[i] + response[i])?;
        }
        let free_increment = core::array::from_fn(|i| free_jump[i] - old[i]);
        let scale = norm_inf(increment)
            .max(norm_inf(free_increment))
            .max(norm_inf(response));
        let tolerance = checked(
            self.config.jump_absolute_tolerance_m + self.config.jump_relative_tolerance * scale,
        )?;
        if norm_inf(residual) > tolerance {
            return Err(CohesiveContactError::NotConverged {
                residual_m: norm_inf(residual),
                tolerance_m: tolerance,
                iterations: 1,
            });
        }
        if !keep_going() {
            return Err(CohesiveContactError::Cancelled);
        }
        Ok(ViscoelasticBondLoadStep {
            material,
            mean_force_n: force,
            compatibility_residual_m: residual,
            compatibility_tolerance_m: tolerance,
        })
    }
}

#[cfg(test)]
mod tests;
