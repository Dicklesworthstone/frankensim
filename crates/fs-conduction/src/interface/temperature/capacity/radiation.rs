//! Original capacity/contact responses inside the shared ambient patch iterator.

use fs_exec::Cx;

use crate::radiation::{AmbientRadiationTrial, solve_ambient_radiation_with};
use crate::transient::backward_euler::{
    BackwardEuler, EndpointEvaluation, NonlinearStepConfig, StepConfig, StepSolution,
};
use crate::{
    AmbientRadiationConfig, AmbientRadiationPatch, AmbientRadiationReport,
    ConductionError, ConductionProblem, ThermalBoundary,
};

use super::super::{
    ContactIterationConfig, TemperatureContactEvaluation, TemperatureContacts,
    finite, invalid, poll,
};
use super::ContactCapacityStep;

/// One unchanged capacity field accepted against live contact and ambient radiation.
#[derive(Debug)]
pub struct ContactCapacityRadiationStep {
    /// Last numerical response, with its applied contact and radiative secant.
    /// Its inner residual is not the complete physical endpoint residual.
    pub conduction: StepSolution,
    /// Contact operator used by that last numerical response.
    pub applied_contacts: TemperatureContactEvaluation,
    /// Pinned sources freshly queried at the returned field.
    pub physical_contacts: TemperatureContactEvaluation,
    /// Complete physical residual, storage and reservoir heat, with original
    /// convection separated from radiation and prescribed storage counted once.
    pub physical: EndpointEvaluation,
    /// Exact frozen radiative boundary used by the final numerical response.
    pub combined_boundary: ThermalBoundary,
    /// Patch heat evidence; iterations counts radiation callbacks, while its
    /// solid/Krylov work includes every nested contact response.
    pub radiation: AmbientRadiationReport,
    /// Fixed absolute full-residual gate supplied by the contact policy, J.
    pub physical_residual_tolerance_j: f64,
    /// All thermal responses, including discarded contact/radiation trials.
    pub thermal_responses: usize,
    /// Newton updates across every thermal response.
    pub nonlinear_updates: usize,
    /// Rejected Newton trials across every thermal response.
    pub nonlinear_backtracks: usize,
    /// Krylov iterations across every thermal response.
    pub krylov_iterations: usize,
}

struct Trial {
    contact: ContactCapacityStep,
    physical: Option<(TemperatureContactEvaluation, EndpointEvaluation)>,
}

impl TemperatureContacts<'_, '_> {
    /// Couple the bounded source-contact fixed point to area-mean ambient gray
    /// radiation. Source surfaces retain their original interpolation without
    /// requiring a unique tangent. Every contact and radiative trial uses the
    /// same full old temperature field, including prescribed-temperature jumps.
    ///
    /// The shared radiation owner supplies secants and raw temperature/watt
    /// gates. A fresh endpoint evaluation additionally checks the complete
    /// contact/radiation nodal residual and both physical and reported heat
    /// balances. Internal contact heat never enters external reservoir input.
    /// No final replacement solve, radiation derivative, continuum bound or
    /// arbitrary-law convergence claim is supplied by this nested fixed point.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_capacity_with_ambient_radiation(
        &self,
        cx: &Cx<'_>,
        engine: &BackwardEuler<'_>,
        problem: ConductionProblem<'_>,
        old: &[f64],
        dt_s: f64,
        step_config: StepConfig,
        nonlinear: Option<NonlinearStepConfig>,
        coupling: ContactIterationConfig,
        patches: &[AmbientRadiationPatch],
        radiation_config: AmbientRadiationConfig,
    ) -> Result<ContactCapacityRadiationStep, ConductionError> {
        self.capacity_ambient(
            cx, engine, problem, old, dt_s, step_config, nonlinear, coupling,
            None, patches, radiation_config,
        )
    }

    /// Use the complete bulk/contact Newton response inside the same ambient
    /// radiation iteration. Continuous curves keep their generalized primal
    /// tangents and source-domain backtracking; radiation remains a bounded
    /// outer fixed point. The explicit face cap bounds each contact tangent.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_capacity_newton_with_ambient_radiation(
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
        patches: &[AmbientRadiationPatch],
        radiation_config: AmbientRadiationConfig,
    ) -> Result<ContactCapacityRadiationStep, ConductionError> {
        self.capacity_ambient(
            cx, engine, problem, old, dt_s, step_config, Some(nonlinear), coupling,
            Some(max_face_pairs), patches, radiation_config,
        )
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn capacity_ambient(
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
        patches: &[AmbientRadiationPatch],
        radiation_config: AmbientRadiationConfig,
    ) -> Result<ContactCapacityRadiationStep, ConductionError> {
        poll(cx, 0)?;
        self.check_mesh(problem.mesh)?;
        radiation_config.validate()?;
        coupling.admit(step_config.linear.max_iterations)?;
        if old.len() != problem.mesh.vertex_count() || max_face_pairs == Some(0) {
            return Err(invalid("<capacity-contact-radiation>",
                "complete original history and a positive Newton contact-face cap are required"));
        }
        if let Some(policy) = nonlinear { policy.validate()?; }
        let responses_per_trial = if max_face_pairs.is_some() { 1 } else { coupling.max_iterations };
        let max_responses = product(responses_per_trial, radiation_config.max_iterations)?;
        let max_krylov = product(max_responses, step_config.linear.max_iterations)?;
        let max_updates = product(max_responses, nonlinear.map_or(0, |policy| policy.max_iterations))?;
        let max_backtracks = product(max_updates,
            nonlinear.map_or(0, |policy| policy.line_search.max_backtracks))?;
        if let Some(policy) = nonlinear {
            let trials = policy.line_search.max_backtracks.checked_add(1)
                .ok_or_else(|| invalid("<capacity-contact-radiation>", "Newton trial allowance overflows"))?;
            product(max_updates, trials)?;
        }

        let mut predictor = old.to_vec();
        for (slot, &(vertex, prescribed)) in problem.boundary.dirichlet().iter().enumerate() {
            if slot % 256 == 0 { poll(cx, slot)?; }
            *predictor.get_mut(vertex).ok_or_else(|| invalid(
                "<capacity-contact-radiation>", "prescribed vertex lies outside original history",
            ))? = prescribed;
        }
        // This fixed predictor retains old free and new prescribed values.
        // Its contact sources are unchanged by later radiative Robin secants.
        let predictor_contacts = self.evaluate(cx, problem.boundary, &predictor)?;
        engine.evaluate_prescribed_endpoint(
            cx, problem, Some(predictor_contacts.frozen_interfaces()), old, &predictor,
            dt_s, step_config, Some(patches),
        )?;
        let target = coupling.residual_tolerance_j;
        let inner_target = finite(target * 0.25)?;
        let numerical_allowance = finite(inner_target * 0.5)?;
        if numerical_allowance <= 0.0 {
            return Err(invalid("<capacity-contact-radiation>", "inner residual allowance is not representable"));
        }
        let inner_coupling = ContactIterationConfig {
            residual_tolerance_j: inner_target,
            ..coupling
        };
        let mut thermal_responses = 0;
        let mut nonlinear_updates = 0;
        let mut nonlinear_backtracks = 0;
        let mut krylov_iterations = 0;
        let accepted = solve_ambient_radiation_with::<_, ConductionError>(
            cx, problem, patches, radiation_config,
            |combined_problem| {
                poll(cx, thermal_responses)?;
                let initial = engine.evaluate_prescribed_endpoint(
                    cx, combined_problem, Some(predictor_contacts.frozen_interfaces()),
                    old, &predictor, dt_s, step_config, None,
                )?;
                let mut inner_step = step_config;
                let inner_nonlinear = if let Some(mut policy) = nonlinear {
                    // Capacity Newton adds these terms; each receives half
                    // the contact allowance, never a widened physical gate.
                    policy.residual_atol_j = policy.residual_atol_j.min(numerical_allowance);
                    policy.residual_rtol = tighter_relative(
                        policy.residual_rtol, numerical_allowance, initial.residual_norm_j,
                    )?;
                    Some(policy)
                } else {
                    inner_step.linear.tolerance = tighter_relative(
                        step_config.linear.tolerance, numerical_allowance, initial.residual_norm_j,
                    )?;
                    None
                };
                let contact = match (max_face_pairs, inner_nonlinear) {
                    (Some(limit), Some(policy)) => self.advance_capacity_newton(
                        cx, engine, combined_problem, old, dt_s, inner_step,
                        policy, inner_coupling, limit,
                    )?,
                    (None, policy) => self.advance_capacity(
                        cx, engine, combined_problem, old, dt_s, inner_step, policy, inner_coupling,
                    )?,
                    _ => return Err(invalid("<capacity-contact-radiation>", "contact Newton lost its nonlinear policy")),
                };
                count(&mut thermal_responses, contact.iterations, max_responses)?;
                count(&mut nonlinear_updates, contact.nonlinear_updates, max_updates)?;
                count(&mut nonlinear_backtracks, contact.nonlinear_backtracks, max_backtracks)?;
                count(&mut krylov_iterations, contact.krylov_iterations, max_krylov)?;
                Ok(AmbientRadiationTrial {
                    robin_fluxes: contact.conduction.robin_fluxes.clone(),
                    robin_out_w: contact.conduction.robin_out_w,
                    solid_iterations: contact.nonlinear_updates,
                    krylov_iterations: contact.krylov_iterations,
                    solution: Trial { contact, physical: None },
                })
            },
            |trial, rows, convective_out_w| {
                let temperature = &trial.contact.conduction.temperature;
                let contacts = self.evaluate(cx, problem.boundary, temperature)?;
                let physical = engine.evaluate_prescribed_endpoint(
                    cx, problem, Some(contacts.frozen_interfaces()), old, temperature,
                    dt_s, step_config, Some(patches),
                )?;
                // Also audit the heat rows returned by the shared radiation
                // owner. Contact is internal; the actual bath reaction includes
                // prescribed-node capacity exactly once in either account.
                let radiation_out_w = rows.iter().try_fold(0.0,
                    |total, row| finite(total + row.nonlinear_heat_w))?;
                let reported_net_input_w = finite(
                    physical.source_w + physical.dirichlet_in_w - physical.neumann_out_w
                        - convective_out_w - radiation_out_w,
                )?;
                let reported_energy_residual_j = finite(
                    physical.stored_energy_change_j - dt_s * reported_net_input_w,
                )?;
                poll(cx, rows.len())?;
                if physical.residual_norm_j > target
                    || physical.energy_residual_j.abs() > step_config.energy_tolerance_j
                    || reported_energy_residual_j.abs() > step_config.energy_tolerance_j
                {
                    return Ok(false);
                }
                trial.physical = Some((contacts, physical));
                Ok(true)
            },
        )?;
        let (physical_contacts, physical) = accepted.solution.physical.ok_or_else(||
            invalid("<capacity-contact-radiation>", "accepted radiation lost its physical contact audit"))?;
        poll(cx, thermal_responses)?;
        Ok(ContactCapacityRadiationStep {
            conduction: accepted.solution.contact.conduction,
            applied_contacts: accepted.solution.contact.applied_contacts,
            physical_contacts,
            physical,
            combined_boundary: accepted.combined_boundary,
            radiation: accepted.radiation,
            physical_residual_tolerance_j: target,
            thermal_responses,
            nonlinear_updates,
            nonlinear_backtracks,
            krylov_iterations,
        })
    }
}

fn tighter_relative(requested: f64, allowance_j: f64, initial_norm_j: f64)
    -> Result<f64, ConductionError>
{
    let tolerance = if initial_norm_j == 0.0 { requested }
        else { requested.min(allowance_j / initial_norm_j) };
    if tolerance.is_finite() && tolerance > 0.0 { Ok(tolerance) }
    else { Err(invalid("<capacity-contact-radiation>", "inner relative tolerance is not representable")) }
}

fn product(a: usize, b: usize) -> Result<usize, ConductionError> {
    a.checked_mul(b).ok_or_else(|| invalid(
        "<capacity-contact-radiation>", "combined radiation/contact work allowance overflows",
    ))
}

fn count(total: &mut usize, added: usize, cap: usize) -> Result<(), ConductionError> {
    *total = total.checked_add(added).filter(|&value| value <= cap)
        .ok_or_else(|| invalid("<capacity-contact-radiation>", "counted numerical work exceeds its allowance"))?;
    Ok(())
}
