//! Full-kinematic material contractions shared by accepted endpoint owners.
use super::super::material::{MaterialGradientBudget3, MaterialScaleGradient3};
use super::super::{TransientAdjointError3, check, validate_values};
use crate::linear3::TetThermalStrainField;
use crate::loads3::transient::{TetTransientRun3, TransientError3, response};
use fs_exec::Cx;

impl TetTransientRun3<'_, '_, '_, '_> {
    pub(in crate::loads3::transient::adjoint) fn admit_material_gradient(
        &self,
        force_cotangent: &[f64],
        budget: MaterialGradientBudget3,
        cx: &Cx<'_>,
    ) -> Result<(), TransientAdjointError3> {
        check(cx)?;
        validate_values(force_cotangent, self.free_dofs().len(), cx)?;
        let full = self
            .problem
            .nodes_m
            .len()
            .checked_mul(3)
            .ok_or(TransientAdjointError3::Budget)?;
        if full > budget.maximum_full_dofs
            || self.problem.tetrahedra.len() > budget.maximum_elements
        {
            return Err(TransientAdjointError3::Budget);
        }
        Ok(())
    }

    /// Contract -lambda' d(M a + (eta M + zeta K) v + K q - fT), with
    /// lambda zero on constrained equations but ACTUAL full endpoint motion.
    /// The initialization mass solve uses the same residual and contraction.
    /// Unlike a reduced DynamicState, this carries every support cross term.
    ///
    /// The internal caller owns the association with this run's exact accepted
    /// endpoint, frozen material/thermal data and fixed step/control schedule.
    pub(in crate::loads3::transient::adjoint) fn material_scale_pullback_nodal_at(
        &self,
        motion: &[Vec<[f64; 3]>; 3],
        force_cotangent: &[f64],
        budget: MaterialGradientBudget3,
        cx: &Cx<'_>,
    ) -> Result<MaterialScaleGradient3, TransientAdjointError3> {
        self.admit_material_gradient(force_cotangent, budget, cx)?;
        let nodes = self.problem.nodes_m.len();
        for field in motion {
            if field.len() != nodes {
                return Err(TransientAdjointError3::Invalid(
                    "material pullback needs full endpoint q/v/a for every node",
                ));
            }
            for chunk in field.chunks(256) {
                check(cx)?;
                if chunk.iter().flatten().any(|value| !value.is_finite()) {
                    return Err(TransientAdjointError3::Invalid(
                        "nonfinite material pullback kinematics",
                    ));
                }
            }
        }
        let mut lambda = vec![[0.0; 3]; nodes];
        for (i, &dof) in self.free_dofs().iter().enumerate() {
            if i % 256 == 0 {
                check(cx)?;
            }
            lambda[dof / 3][dof % 3] = force_cotangent[i];
        }
        let thermal = self.thermal_strains.as_ref().map(|field| match field {
            response::FrozenThermal3::Uniform(value) => TetThermalStrainField::Uniform(value),
            response::FrozenThermal3::PerElement(values) => {
                TetThermalStrainField::PerElement(values)
            }
        });
        let (element_scales, rayleigh) = self
            .problem
            .transient_material_contraction(
                motion,
                &lambda,
                thermal,
                self.config.mass_damping_per_s,
                self.config.stiffness_damping_s,
                cx,
            )
            .map_err(TransientError3::from)?;
        check(cx)?;
        Ok(MaterialScaleGradient3 {
            element_scales,
            rayleigh,
        })
    }
}
