//! Original heterogeneous storage with source-dependent ambient emissivity.
use super::*;

impl HeterogeneousEnthalpyBackwardEuler<'_, '_> {
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_material_ambient_radiation(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>, old_h: &[f64], dt_s: f64,
        config: EnthalpyStepConfig, radiation: &crate::MaterialPointwiseAmbientRadiation<'_>,
    ) -> Result<EnthalpyStepSolution, EnthalpyError> {
        self.inner.advance_with_material_ambient_radiation(
            cx, problem, interfaces, old_h, dt_s, config, radiation,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_endpoint_with_material_ambient_radiation(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>, old_h: &[f64], candidate_h: &[f64],
        dt_s: f64, config: EnthalpyStepConfig, radiation: &crate::MaterialPointwiseAmbientRadiation<'_>,
    ) -> Result<super::super::EnthalpyEndpointEvaluation, EnthalpyError> {
        self.inner.evaluate_endpoint_with_material_ambient_radiation(
            cx, problem, interfaces, old_h, candidate_h, dt_s, config, radiation,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn advance_with_temperature_contacts_and_material_ambient_radiation(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>,
        contacts: &crate::interface::temperature::TemperatureContacts<'_, '_>,
        old_h: &[f64], dt_s: f64,
        config: crate::interface::temperature::enthalpy::ContactEnthalpyNewtonConfig,
        radiation: &crate::MaterialPointwiseAmbientRadiation<'_>,
    ) -> Result<(EnthalpyStepSolution, crate::interface::temperature::TemperatureContactEvaluation), EnthalpyError> {
        self.inner.advance_with_temperature_contacts(cx, problem, contacts, old_h, dt_s, config,
            Some(super::super::PointwiseBoundary::MaterialAmbient(radiation)))
    }
}
