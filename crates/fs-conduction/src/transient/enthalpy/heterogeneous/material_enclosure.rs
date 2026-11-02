//! Original heterogeneous storage with evolving material reflection.
use super::*;

#[derive(Debug)]
pub struct HeterogeneousMaterialEnclosureEnthalpyGradient {
    pub thermal: HeterogeneousEnthalpyStepGradient,
    pub radiation: crate::radiation::pointwise_enclosure::material::adjoint::MaterialEnclosureParameterGradient,
}

impl<'m> HeterogeneousEnthalpyBackwardEuler<'m, '_> {
    /// Bind original material reflection and heterogeneous storage controls at
    /// the unchanged accepted endpoint, with strict source/chart admission.
    #[allow(clippy::too_many_arguments)]
    pub fn linearize_accepted_with_material_enclosure_radiation(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>, old_h: &[f64], dt_s: f64,
        config: EnthalpyStepConfig, radiation: &crate::MaterialPointwiseEnclosureRadiation<'_>,
        accepted: EnthalpyStepSolution, physical_target_j: f64, max_radiation_entries: usize,
    ) -> Result<HeterogeneousEnthalpyStepLinearization<'m>, EnthalpyAdjointError> {
        let inner = self.inner.linearize_accepted_with_material_enclosure_radiation(
            cx, problem, interfaces, old_h, dt_s, config, radiation, accepted,
            physical_target_j, max_radiation_entries)?;
        self.bind_linearization(cx, old_h, inner)
    }
}

impl HeterogeneousEnthalpyStepLinearization<'_> {
    /// One physical transpose supplies reference-density, source/history and
    /// original emissivity/open-bath controls, including reflection feedback.
    pub fn pullback_with_material_enclosure_controls(
        &self, cx: &Cx<'_>, seed: &[f64], linear: LinearConfig,
    ) -> Result<HeterogeneousMaterialEnclosureEnthalpyGradient, EnthalpyAdjointError> {
        let thermal = self.pullback(cx, seed, linear)?;
        let mut weights = self.inner.nodal_load_pullback(cx, &thermal.transport.adjoint)?;
        for (vertex, weight) in weights.iter_mut().enumerate() {
            if vertex % ASSEMBLY_TILE == 0 { poll(cx, vertex)?; }
            *weight = -*weight;
        }
        let radiation = self.inner.material_enclosure_parameter_pullback(cx, &weights)?;
        Ok(HeterogeneousMaterialEnclosureEnthalpyGradient { thermal, radiation })
    }
}

impl HeterogeneousEnthalpyBackwardEuler<'_, '_> {
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_material_enclosure_radiation(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>, old_h: &[f64], dt_s: f64,
        config: EnthalpyStepConfig, radiation: &crate::MaterialPointwiseEnclosureRadiation<'_>,
    ) -> Result<EnthalpyStepSolution, EnthalpyError> {
        self.inner.advance_with_material_enclosure_radiation(cx, problem, interfaces, old_h, dt_s, config, radiation)
    }

    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_endpoint_with_material_enclosure_radiation(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>, old_h: &[f64], candidate_h: &[f64],
        dt_s: f64, config: EnthalpyStepConfig, radiation: &crate::MaterialPointwiseEnclosureRadiation<'_>,
    ) -> Result<super::super::EnthalpyEndpointEvaluation, EnthalpyError> {
        self.inner.evaluate_endpoint_with_material_enclosure_radiation(
            cx, problem, interfaces, old_h, candidate_h, dt_s, config, radiation)
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn advance_with_temperature_contacts_and_material_enclosure_radiation(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>,
        contacts: &crate::interface::temperature::TemperatureContacts<'_, '_>,
        old_h: &[f64], dt_s: f64,
        config: crate::interface::temperature::enthalpy::ContactEnthalpyNewtonConfig,
        radiation: &crate::MaterialPointwiseEnclosureRadiation<'_>,
    ) -> Result<(EnthalpyStepSolution, crate::interface::temperature::TemperatureContactEvaluation), EnthalpyError> {
        self.inner.advance_with_temperature_contacts(cx, problem, contacts, old_h, dt_s, config,
            Some(super::super::PointwiseBoundary::MaterialEnclosure(radiation)))
    }
}
