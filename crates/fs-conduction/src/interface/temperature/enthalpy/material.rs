//! Material emissivity and source contact in the original coupled h solve.
use super::*;

impl TemperatureContacts<'_, '_> {
    /// Solve contact resistance, nodal material emissivity and enthalpy together.
    /// The direct emissivity/phase feedback remains active on latent plateaus.
    /// Every surface retains its exact chart, pinned source and external bath;
    /// radiation and original convection keep separate physical heat accounts.
    /// Material reports are retained in `conduction.pointwise_radiation`.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_enthalpy_newton_with_material_ambient_radiation(
        &self, cx: &Cx<'_>, engine: &HeterogeneousEnthalpyBackwardEuler<'_, '_>,
        problem: ConductionProblem<'_>, old_h: &[f64], dt_s: f64,
        config: ContactEnthalpyNewtonConfig,
        radiation: &crate::MaterialPointwiseAmbientRadiation<'_>,
    ) -> Result<ContactEnthalpyStep, EnthalpyError> {
        let (solved, actual) = engine.advance_with_temperature_contacts_and_material_ambient_radiation(
            cx, problem, self, old_h, dt_s, config, radiation,
        )?;
        let physical = engine.evaluate_endpoint_with_material_ambient_radiation(
            cx, problem, Some(actual.frozen_interfaces()), old_h,
            &solved.specific_enthalpy_j_kg, dt_s, config.step, radiation,
        )?;
        accept_coupled_endpoint(cx, solved, actual, physical, config)
    }
}
