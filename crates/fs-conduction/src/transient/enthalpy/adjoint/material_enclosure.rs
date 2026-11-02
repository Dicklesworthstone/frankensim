//! Accepted material enclosure controls on the original complete h transpose.

use super::*;
use crate::MaterialPointwiseEnclosureRadiation;
use crate::radiation::pointwise_enclosure::material::adjoint::MaterialEnclosureParameterGradient;

#[derive(Debug)]
pub struct MaterialEnclosureEnthalpyGradient {
    pub thermal: EnthalpyStepGradient,
    pub radiation: MaterialEnclosureParameterGradient,
}

impl<'m> EnthalpyBackwardEuler<'m, '_> {
    /// Recheck the unchanged physical endpoint and retain its complete h
    /// tangent, including original material emissivity and reflected irradiation.
    /// Source/phase corners with moving coordinates refuse typed two-sided
    /// admission. Geometry, chart and source knot abscissae remain fixed.
    #[allow(clippy::too_many_arguments)]
    pub fn linearize_accepted_with_material_enclosure_radiation(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        radiation: &MaterialPointwiseEnclosureRadiation<'_>,
        accepted: EnthalpyStepSolution,
        physical_target_j: f64,
        max_radiation_entries: usize,
    ) -> Result<EnthalpyStepLinearization<'m>, EnthalpyAdjointError> {
        poll(cx, 0)?;
        if !physical_target_j.is_finite() || physical_target_j < 0.0
            || !std::ptr::eq(self.mesh, problem.mesh)
            || !radiation.matches(problem.mesh, problem.boundary)
            || !radiation.matches_curves(cx, |vertex| self.curve_for_vertex(vertex))?
            || accepted.pointwise_radiation.is_some()
            || !accepted.pointwise_enclosure_radiation.as_ref().is_some_and(|report| {
                report.patches.iter().all(|patch| !patch.material_states.is_empty())
            })
        {
            return Err(EnthalpyAdjointError::InvalidInput(
                "material enclosure adjoint requires the original physical target, exact storage/boundary and accepted source-state report",
            ));
        }
        if radiation.adjoint_workspace_entries()? > max_radiation_entries {
            return Err(EnthalpyAdjointError::InvalidInput(
                "material enclosure adjoint exceeds its declared entry allowance",
            ));
        }
        self.linearize_accepted_model_with_target(cx, problem, interfaces, old_h, dt_s,
            config, accepted, Some(physical_target_j),
            Some(PointwiseBoundary::MaterialEnclosure(radiation)), Some(max_radiation_entries))
    }
}

impl EnthalpyStepLinearization<'_> {
    /// Contract original emissivity controls at fixed h using outward-load
    /// weights. Reflected state feedback belongs to the complete h transpose.
    pub fn material_enclosure_parameter_pullback(
        &self,
        cx: &Cx<'_>,
        outward_load_weights: &[f64],
    ) -> Result<MaterialEnclosureParameterGradient, EnthalpyAdjointError> {
        vector(cx, outward_load_weights, self.masses.len())?;
        match &self.pointwise_radiation {
            Some(PointwiseBoundaryAssembly::MaterialEnclosure(assembly)) => {
                Ok(assembly.material_parameter_pullback(cx, outward_load_weights)?)
            }
            _ => Err(EnthalpyAdjointError::InvalidInput(
                "original material enclosure controls require their retained complete radiation tangent",
            )),
        }
    }

    /// One complete transpose returns history/source and original emissivity/
    /// open-bath controls. Direct objective control dependence is caller-owned.
    pub fn pullback_with_material_enclosure_controls(
        &self,
        cx: &Cx<'_>,
        seed: &[f64],
        linear: LinearConfig,
    ) -> Result<MaterialEnclosureEnthalpyGradient, EnthalpyAdjointError> {
        let thermal = self.pullback(cx, seed, linear)?;
        let mut weights = self.nodal_load_pullback(cx, &thermal.adjoint)?;
        for (vertex, weight) in weights.iter_mut().enumerate() {
            if vertex % ASSEMBLY_TILE == 0 { poll(cx, vertex)?; }
            *weight = -*weight;
        }
        let radiation = self.material_enclosure_parameter_pullback(cx, &weights)?;
        Ok(MaterialEnclosureEnthalpyGradient { thermal, radiation })
    }
}
