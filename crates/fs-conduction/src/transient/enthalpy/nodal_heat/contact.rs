//! Held signed nodal heat with the original live contact Newton operator.
//!
//! Every nonlinear trial queries the original contact sources at its actual
//! temperature and phase traces. The native tangent retains their derivatives,
//! nonmatching quadrature and optional pointwise radiation. Geometry, pressure,
//! gaps and other nonthermal source coordinates stay held by the contact owner.
//! This is backward Euler on the original reference inventory, not a mechanical
//! contact or moving-boundary model.

use super::super::heterogeneous::EnthalpyRadiation;
use super::*;
use crate::interface::temperature::TemperatureContacts;
use crate::interface::temperature::enthalpy::ContactEnthalpyNewtonConfig;

impl EnthalpyBackwardEuler<'_, '_> {
    /// Solve original temperature/phase-dependent contacts and optional native
    /// pointwise radiation with signed interval heat in the same Newton solve.
    ///
    /// The heat is held during this solve. An endpoint-dependent physical
    /// owner must recompute it and audit the resulting endpoint before advancing
    /// its history. No frozen-contact outer iteration or second tangent is used.
    /// Both the contact residual norm and the original local/global joule gates
    /// must pass. The endpoint retains the actual original contact receipts.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_nodal_heat_and_temperature_contacts(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        contacts: &TemperatureContacts<'_, '_>,
        old_specific_h: &[f64],
        dt_s: f64,
        config: ContactEnthalpyNewtonConfig,
        nodal_heat_j: &[Dd],
        radiation: EnthalpyRadiation<'_>,
    ) -> Result<NodalHeatStep, EnthalpyError> {
        poll(cx, 0)?;
        let heat = admit_heat(cx, self.masses.len(), nodal_heat_j)?;
        let radiation = super::radiation::pointwise(self, cx, problem, radiation)?;
        let (h, newton, actual) = self.solve_with_temperature_contacts(
            cx,
            problem,
            contacts,
            old_specific_h,
            dt_s,
            config,
            radiation,
            Some(&heat),
            None,
        )?;
        let context = StepContext {
            storage: self,
            cx,
            problem,
            interfaces: Some(actual.frozen_interfaces()),
            old: old_specific_h,
            dt: dt_s,
            pointwise_radiation: radiation,
        };
        let mut endpoint = evaluate(&context, &h, &heat)?;
        if endpoint.residual_norm_j > config.residual_tolerance_j {
            return Err(EnthalpyError::EndpointResidual {
                residual_j: endpoint.residual_norm_j,
                tolerance_j: config.residual_tolerance_j,
            });
        }
        endpoint.temperature_contacts = Some(actual);
        self.finish_nodal_heat(cx, h, newton, endpoint, config.step)
    }

    /// Re-query the original contacts and radiation at a supplied enthalpy
    /// endpoint with its freshly computed signed heat. This takes no time step.
    ///
    /// Exact mesh, material-caloric identities, prescribed states and source
    /// domains remain mandatory. The caller must admit the contact residual
    /// norm, maximum free-row defect and independent global energy balance.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_endpoint_with_nodal_heat_and_temperature_contacts(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        contacts: &TemperatureContacts<'_, '_>,
        old_specific_h: &[f64],
        candidate_specific_h: &[f64],
        dt_s: f64,
        config: ContactEnthalpyNewtonConfig,
        nodal_heat_j: &[Dd],
        radiation: EnthalpyRadiation<'_>,
    ) -> Result<NodalHeatEndpoint, EnthalpyError> {
        self.admit_step(cx, problem, old_specific_h, dt_s, config.step)?;
        config.validate()?;
        if candidate_specific_h.len() != self.masses.len() {
            return Err(EnthalpyError::InvalidInput(
                "candidate enthalpy must match the original reference inventory",
            ));
        }
        let heat = admit_heat(cx, self.masses.len(), nodal_heat_j)?;
        let radiation = super::radiation::pointwise(self, cx, problem, radiation)?;
        contacts.check_mesh(problem.mesh)?;
        contacts.check_material_curves(cx, |vertex| Some(self.curve_for_vertex(vertex)))?;
        let mut context = StepContext {
            storage: self,
            cx,
            problem,
            interfaces: None,
            old: old_specific_h,
            dt: dt_s,
            pointwise_radiation: radiation,
        };
        context.check_prescribed_values(candidate_specific_h)?;
        let temperature = context.temperatures(candidate_specific_h)?;
        let actual = contacts.evaluate_enthalpy_with_query_refusal(
            cx,
            problem.boundary,
            candidate_specific_h,
            &temperature,
            &mut |_| {},
        )?;
        context.interfaces = Some(actual.frozen_interfaces());
        let mut endpoint = evaluate(&context, candidate_specific_h, &heat)?;
        endpoint.temperature_contacts = Some(actual);
        Ok(endpoint)
    }
}
