//! Accepted-state material-ambient adjoints through the original full h owner.

use super::*;
use crate::MaterialPointwiseAmbientRadiation;
use crate::radiation::pointwise::material::adjoint::MaterialRadiationParameterGradient;

/// One physical transpose supplies thermal history/source and original
/// emissivity/bath controls. No nonlinear iteration is differentiated.
#[derive(Debug)]
pub struct MaterialAmbientEnthalpyGradient {
    pub thermal: EnthalpyStepGradient,
    pub radiation: MaterialRadiationParameterGradient,
}

impl<'m> EnthalpyBackwardEuler<'m, '_> {
    /// Recheck an unchanged material-emissivity endpoint and retain its complete
    /// enthalpy tangent, including phase feedback when dT/dh is zero.
    ///
    /// The original mesh, storage charts, boundary and source receipts must
    /// agree. Caloric/source corners and unsupported two-sided neighborhoods
    /// refuse with typed errors; the forward Newton's one-sided slope does not
    /// establish an adjoint. Scalar values and original curve ordinates are
    /// available as controls at fixed abscissae, chart and source conditions.
    #[allow(clippy::too_many_arguments)]
    pub fn linearize_accepted_with_material_ambient_radiation(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        radiation: &MaterialPointwiseAmbientRadiation<'_>,
        accepted: EnthalpyStepSolution,
        physical_target_j: f64,
        max_radiation_entries: usize,
    ) -> Result<EnthalpyStepLinearization<'m>, EnthalpyAdjointError> {
        poll(cx, 0)?;
        if !physical_target_j.is_finite()
            || physical_target_j < 0.0
            || !std::ptr::eq(self.mesh, problem.mesh)
            || !radiation.matches(problem.mesh, problem.boundary)
            || !radiation.matches_curves(cx, |vertex| self.curve_for_vertex(vertex))?
            || accepted.pointwise_enclosure_radiation.is_some()
            || !accepted.pointwise_radiation.as_ref().is_some_and(|report| {
                report
                    .patches
                    .iter()
                    .all(|patch| !patch.material_states.is_empty())
            })
        {
            return Err(EnthalpyAdjointError::InvalidInput(
                "material ambient adjoint requires the original physical target, exact storage/boundary and accepted source-state report",
            ));
        }
        if radiation.adjoint_workspace_entries()? > max_radiation_entries {
            return Err(EnthalpyAdjointError::InvalidInput(
                "material radiation adjoint exceeds its declared entry allowance",
            ));
        }
        self.linearize_accepted_model_with_target(
            cx,
            problem,
            interfaces,
            old_h,
            dt_s,
            config,
            accepted,
            Some(physical_target_j),
            Some(PointwiseBoundary::MaterialAmbient(radiation)),
            Some(max_radiation_entries),
        )
    }
}

impl EnthalpyStepLinearization<'_> {
    /// Contract original material emissivity controls using outward-load
    /// weights. State feedback already belongs to the complete h transpose.
    pub fn material_ambient_parameter_pullback(
        &self,
        cx: &Cx<'_>,
        outward_load_weights: &[f64],
    ) -> Result<MaterialRadiationParameterGradient, EnthalpyAdjointError> {
        vector(cx, outward_load_weights, self.masses.len())?;
        match &self.pointwise_radiation {
            Some(PointwiseBoundaryAssembly::Ambient(assembly)) => {
                Ok(assembly.material_parameter_pullback(cx, outward_load_weights)?)
            }
            _ => Err(EnthalpyAdjointError::InvalidInput(
                "original material ambient controls require their retained complete radiation tangent",
            )),
        }
    }

    /// Solve one complete transpose and return thermal plus original-source
    /// controls. Direct objective dependence on controls belongs to the caller.
    pub fn pullback_with_material_ambient_controls(
        &self,
        cx: &Cx<'_>,
        seed: &[f64],
        linear: LinearConfig,
    ) -> Result<MaterialAmbientEnthalpyGradient, EnthalpyAdjointError> {
        let thermal = self.pullback(cx, seed, linear)?;
        let mut weights = self.nodal_load_pullback(cx, &thermal.adjoint)?;
        for (vertex, weight) in weights.iter_mut().enumerate() {
            if vertex % ASSEMBLY_TILE == 0 {
                poll(cx, vertex)?;
            }
            *weight = -*weight;
        }
        let radiation = self.material_ambient_parameter_pullback(cx, &weights)?;
        Ok(MaterialAmbientEnthalpyGradient { thermal, radiation })
    }
}
