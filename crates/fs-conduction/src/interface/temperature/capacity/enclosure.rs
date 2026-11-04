//! Live capacity/contact responses inside the existing closed-enclosure owner.

use fs_exec::Cx;

use crate::radiation::{EnclosureExchangeReport, GrayDiffuseEnclosure};
use crate::transient::backward_euler::{
    BackwardEuler, EndpointEvaluation, NonlinearStepConfig, StepConfig, StepSolution,
    enclosure::EnclosureStepResponse,
};
use crate::{AmbientRadiationConfig, ConductionError, ConductionProblem, ThermalBoundary};

use super::super::{
    ContactIterationConfig, TemperatureContactEvaluation, TemperatureContacts, invalid, poll,
};

/// One unchanged capacity field accepted against live contact and closed exchange.
#[derive(Debug)]
pub struct ContactCapacityEnclosureStep {
    /// Last numerical response, with its applied contact and enclosure boundary.
    /// Its inner residual is not the complete physical endpoint residual.
    pub conduction: StepSolution,
    /// Contact operator used by that last numerical response.
    pub applied_contacts: TemperatureContactEvaluation,
    /// Pinned sources freshly queried at the returned field.
    pub physical_contacts: TemperatureContactEvaluation,
    /// Fresh full residual and storage balance, with original convection and
    /// Neumann heat separated from internal reflected exchange.
    pub physical: EndpointEvaluation,
    /// Frozen reflected boundary used by the final numerical response.
    pub combined_boundary: ThermalBoundary,
    /// Fresh radiosity and applied heat evidence. Iterations counts enclosure
    /// callbacks; solid/Krylov work includes all nested contact responses.
    pub radiation: EnclosureExchangeReport,
    /// Original absolute full-residual allowance from the contact policy, J.
    pub physical_residual_tolerance_j: f64,
    /// Every thermal response, including discarded contact/enclosure trials.
    pub thermal_responses: usize,
    /// Newton updates across every thermal response.
    pub nonlinear_updates: usize,
    /// Rejected Newton trials across every thermal response.
    pub nonlinear_backtracks: usize,
    /// Krylov iterations across every thermal response.
    pub krylov_iterations: usize,
}

impl TemperatureContacts<'_, '_> {
    /// Couple source-surface contact fixed points to closed area-mean gray
    /// enclosure exchange. Each trial retains the original full old field;
    /// only freshly audited complete contact/radiation balance may be accepted.
    /// The shared enclosure owner preserves authored convection and Neumann
    /// loads, prescribed-node storage reactions and internal radiative heat.
    /// No pointwise/open exchange or monolithic radiation tangent is supplied.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_capacity_with_enclosure_radiation(
        &self,
        cx: &Cx<'_>,
        engine: &BackwardEuler<'_>,
        problem: ConductionProblem<'_>,
        old: &[f64],
        dt_s: f64,
        step_config: StepConfig,
        nonlinear: Option<NonlinearStepConfig>,
        coupling: ContactIterationConfig,
        enclosure: &GrayDiffuseEnclosure,
        radiation_config: AmbientRadiationConfig,
    ) -> Result<ContactCapacityEnclosureStep, ConductionError> {
        self.capacity_enclosure(
            cx, engine, problem, old, dt_s, step_config, nonlinear, coupling,
            None, enclosure, radiation_config,
        )
    }

    /// Use the complete bulk/contact Newton response inside that same closed
    /// enclosure iteration. Continuous curves retain generalized primal
    /// tangents and source-domain backtracking under the declared face cap.
    /// Reflection remains a bounded outer fixed point with independent fresh
    /// nodal and energy gates; it is not a contact/radiation sensitivity owner.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_capacity_newton_with_enclosure_radiation(
        &self,
        cx: &Cx<'_>,
        engine: &BackwardEuler<'_>,
        problem: ConductionProblem<'_>,
        old: &[f64],
        dt_s: f64,
        step_config: StepConfig,
        nonlinear: NonlinearStepConfig,
        coupling: ContactIterationConfig,
        max_face_pairs: usize,
        enclosure: &GrayDiffuseEnclosure,
        radiation_config: AmbientRadiationConfig,
    ) -> Result<ContactCapacityEnclosureStep, ConductionError> {
        self.capacity_enclosure(
            cx, engine, problem, old, dt_s, step_config, Some(nonlinear), coupling,
            Some(max_face_pairs), enclosure, radiation_config,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn capacity_enclosure(
        &self,
        cx: &Cx<'_>,
        engine: &BackwardEuler<'_>,
        problem: ConductionProblem<'_>,
        old: &[f64],
        dt_s: f64,
        step_config: StepConfig,
        nonlinear: Option<NonlinearStepConfig>,
        coupling: ContactIterationConfig,
        max_face_pairs: Option<usize>,
        enclosure: &GrayDiffuseEnclosure,
        radiation_config: AmbientRadiationConfig,
    ) -> Result<ContactCapacityEnclosureStep, ConductionError> {
        poll(cx, 0)?;
        self.check_mesh(problem.mesh)?;
        coupling.admit(step_config.linear.max_iterations)?;
        if max_face_pairs == Some(0) {
            return Err(invalid("<capacity-contact-enclosure>", "Newton contact-face cap must be positive"));
        }
        let responses_per_trial = if max_face_pairs.is_some() { 1 } else { coupling.max_iterations };
        let accepted = engine.advance_enclosure_with_responses(
            cx, problem, old, dt_s, step_config, nonlinear, enclosure, radiation_config,
            Some(coupling.residual_tolerance_j), responses_per_trial,
            |combined, inner_step, policy, inner_allowance_j| {
                let inner_coupling = ContactIterationConfig {
                    residual_tolerance_j: coupling.residual_tolerance_j.min(inner_allowance_j),
                    ..coupling
                };
                let contact = match (max_face_pairs, policy) {
                    (Some(limit), Some(policy)) => self.advance_capacity_newton(
                        cx, engine, combined, old, dt_s, inner_step, policy, inner_coupling, limit,
                    )?,
                    (None, policy) => self.advance_capacity(
                        cx, engine, combined, old, dt_s, inner_step, policy, inner_coupling,
                    )?,
                    _ => return Err(invalid("<capacity-contact-enclosure>", "contact Newton lost its nonlinear policy")),
                };
                Ok(EnclosureStepResponse {
                    conduction: contact.conduction,
                    nonlinear: None,
                    retained: contact.applied_contacts,
                    thermal_responses: contact.iterations,
                    nonlinear_updates: contact.nonlinear_updates,
                    nonlinear_backtracks: contact.nonlinear_backtracks,
                    krylov_iterations: contact.krylov_iterations,
                })
            },
            |combined, candidate| {
                let contacts = self.evaluate(cx, combined.boundary, candidate)?;
                let physical = engine.evaluate_prescribed_endpoint(
                    cx, combined, Some(contacts.frozen_interfaces()), old, candidate,
                    dt_s, step_config, None,
                )?;
                Ok((physical, contacts))
            },
        )?;
        poll(cx, accepted.thermal_responses)?;
        Ok(ContactCapacityEnclosureStep {
            conduction: accepted.response.conduction,
            applied_contacts: accepted.response.retained,
            physical_contacts: accepted.physical_state,
            physical: accepted.physical,
            combined_boundary: accepted.combined_boundary,
            radiation: accepted.radiation,
            physical_residual_tolerance_j: accepted.physical_residual_tolerance_j,
            thermal_responses: accepted.thermal_responses,
            nonlinear_updates: accepted.nonlinear_updates,
            nonlinear_backtracks: accepted.nonlinear_backtracks,
            krylov_iterations: accepted.krylov_iterations,
        })
    }
}
