//! Temperature-dependent contact with primary enthalpy and immutable history.

use fs_exec::Cx;
use crate::{AmbientRadiationConfig, AmbientRadiationPatch, AmbientRadiationReport,
    ConductionProblem, ThermalBoundary};
use crate::transient::enthalpy::{
    EnthalpyEndpointEvaluation, EnthalpyError, EnthalpyStepConfig, EnthalpyStepSolution,
};
use crate::transient::enthalpy::heterogeneous::HeterogeneousEnthalpyBackwardEuler;
use super::{ContactIterationConfig, TemperatureContactEvaluation, TemperatureContacts, invalid, poll};

/// One accepted enthalpy field and the original and re-evaluated contact laws.
#[derive(Debug)]
pub struct ContactEnthalpyStep {
    /// Last counted nonlinear response. Advance history with its specific h,
    /// never by attempting to invert its possibly isothermal temperature field.
    pub conduction: EnthalpyStepSolution,
    /// Source states and operator actually used by the final inner response.
    pub applied_contacts: TemperatureContactEvaluation,
    /// Fresh source states at the unchanged, chart-resolved endpoint temperatures.
    pub physical_contacts: TemperatureContactEvaluation,
    /// Full nodal enthalpy residual and independent physical energy audit.
    pub physical: EnthalpyEndpointEvaluation,
    /// All contact trials, including discarded nonlinear responses.
    pub iterations: usize,
    /// Sum of actual Newton iterations over every contact trial.
    pub nonlinear_iterations: usize,
    /// Sum of actual Krylov iterations over every contact trial.
    pub krylov_iterations: usize,
    /// Final inner radiation report; its counters cover that contact trial.
    pub radiation: Option<AmbientRadiationReport>,
    /// Exact secant boundary used by the final radiative enthalpy response.
    /// The separate `physical` audit uses the actual endpoint radiation law.
    pub combined_boundary: Option<ThermalBoundary>,
    /// All radiation trials across every contact iteration, not just the last.
    pub radiation_iterations: usize,
}

impl TemperatureContacts<'_, '_> {
    /// Advance heterogeneous enthalpy with source-backed R(T_interface).
    ///
    /// Every contact trial calls the existing spatial enthalpy solver from the
    /// SAME old h. That owner retains heterogeneous reference masses, distinct
    /// charts, nonlinear conductivity and optional phase-conductivity feedback.
    /// Only the contact-query temperature guess is relaxed. Neither h nor the
    /// final field is mixed or reconstructed from temperature.
    ///
    /// Requery every face at the candidate's actual temperature and reassemble
    /// the complete nodal joule residual before acceptance. A latent plateau or
    /// a closed total-energy sum cannot conceal a local enthalpy imbalance.
    /// The original contact residual target and energy tolerance never widen.
    /// Inner Newton tolerances must be tight enough to satisfy those gates.
    ///
    /// Prescribed Neumann/Robin data and fixed contact geometry are supported.
    /// Dirichlet enthalpy, radiation/air feedback and nonlinear-contact adjoints
    /// are not supplied by this entry point. No state is returned on failure.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_enthalpy(
        &self,
        cx: &Cx<'_>,
        engine: &HeterogeneousEnthalpyBackwardEuler<'_, '_>,
        problem: ConductionProblem<'_>,
        old_h: &[f64],
        dt_s: f64,
        step_config: EnthalpyStepConfig,
        coupling: ContactIterationConfig,
    ) -> Result<ContactEnthalpyStep, EnthalpyError> {
        self.advance_enthalpy_inner(cx, engine, problem, old_h, dt_s, step_config, coupling, None)
    }

    /// Compose source-backed contact with the existing nonlinear gray-radiation
    /// owner. Every contact/radiation/Newton trial keeps the same physical h.
    /// The final independent audit re-evaluates BOTH R(T_interface) and radiation;
    /// radiation is external heat, while contact transfer remains internal.
    /// Prescribed convection retains its own heat account and may share an
    /// exterior patch with radiation under the existing composition rules.
    ///
    /// Work bounds cover the complete product of contact trials, radiation
    /// trials, Newton attempts and Krylov cycles. A failed inner solve, source
    /// query or complete physical gate returns no field or altered history.
    /// This does not add airflow feedback or a nonlinear-contact adjoint.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_enthalpy_with_ambient_radiation(
        &self,
        cx: &Cx<'_>,
        engine: &HeterogeneousEnthalpyBackwardEuler<'_, '_>,
        problem: ConductionProblem<'_>,
        old_h: &[f64],
        dt_s: f64,
        step_config: EnthalpyStepConfig,
        coupling: ContactIterationConfig,
        patches: &[AmbientRadiationPatch],
        radiation_config: AmbientRadiationConfig,
    ) -> Result<ContactEnthalpyStep, EnthalpyError> {
        self.advance_enthalpy_inner(cx, engine, problem, old_h, dt_s, step_config,
            coupling, Some((patches, radiation_config)))
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn advance_enthalpy_inner(
        &self,
        cx: &Cx<'_>,
        engine: &HeterogeneousEnthalpyBackwardEuler<'_, '_>,
        problem: ConductionProblem<'_>,
        old_h: &[f64],
        dt_s: f64,
        step_config: EnthalpyStepConfig,
        coupling: ContactIterationConfig,
        radiation: Option<(&[AmbientRadiationPatch], AmbientRadiationConfig)>,
    ) -> Result<ContactEnthalpyStep, EnthalpyError> {
        poll(cx, 0)?;
        self.check_mesh(problem.mesh)?;
        if !problem.boundary.dirichlet().is_empty() {
            return Err(EnthalpyError::UnsupportedDirichlet);
        }
        step_config.newton.validate().map_err(EnthalpyError::Newton)?;
        let radiation_factor = match radiation {
            Some((patches, config)) => {
                config.validate()?;
                if patches.is_empty() {
                    return Err(EnthalpyError::InvalidInput("contact radiation requires declared patches"));
                }
                config.max_iterations
            },
            None => 1,
        };
        let max_responses = coupling.max_iterations.checked_mul(radiation_factor)
            .ok_or(EnthalpyError::InvalidInput("contact radiation work allowance overflows"))?;
        let per_response = radiation_factor.checked_mul(step_config.max_newton_iterations)
            .and_then(|n| n.checked_mul(step_config.newton.linear_restart))
            .and_then(|n| n.checked_mul(step_config.newton.max_linear_cycles))
            .ok_or(EnthalpyError::InvalidInput("contact enthalpy inner work allowance overflows"))?;
        let max_krylov = coupling.admit(per_response)?;
        let max_newton = max_responses.checked_mul(step_config.max_newton_iterations)
            .ok_or(EnthalpyError::InvalidInput("contact enthalpy Newton allowance overflows"))?;
        if old_h.len() != problem.mesh.vertex_count()
            || engine.vertex_material_ids().len() != old_h.len()
        {
            return Err(EnthalpyError::InvalidInput("contact enthalpy history and reference mesh differ"));
        }
        let mut guess = Vec::new();
        guess.try_reserve_exact(old_h.len())
            .map_err(|_| EnthalpyError::InvalidInput("contact enthalpy query field allocation refused"))?;
        for (vertex, &h) in old_h.iter().enumerate() {
            if vertex % 256 == 0 { poll(cx, vertex)?; }
            let curve = engine.phase_curve_at(vertex)
                .ok_or(EnthalpyError::InvalidInput("contact enthalpy vertex has no phase chart"))?;
            let state = curve.state_at_specific_enthalpy(h)
                .map_err(|source| EnthalpyError::Phase { vertex, source })?;
            guess.push(state.temperature_k());
        }
        let mut newton = 0usize;
        let mut krylov = 0usize;
        let mut radiation_iterations = 0usize;
        let mut last_residual = f64::INFINITY;
        for iteration in 1..=coupling.max_iterations {
            poll(cx, iteration)?;
            let applied = self.evaluate(cx, problem.boundary, &guess)?;
            let (solved, report, boundary, trial_newton, trial_krylov) = match radiation {
                Some((patches, config)) => {
                    let response = engine.advance_with_ambient_radiation(cx, problem,
                        Some(applied.frozen_interfaces()), old_h, dt_s, step_config, patches, config)?;
                    radiation_iterations = radiation_iterations.checked_add(response.radiation.iterations)
                        .filter(|&n| n <= max_responses)
                        .ok_or(EnthalpyError::InvalidInput("contact radiation work exceeds allowance"))?;
                    let work = (response.radiation.solid_iterations, response.radiation.krylov_iterations);
                    (response.conduction, Some(response.radiation), Some(response.combined_boundary), work.0, work.1)
                },
                None => {
                    let solved = engine.advance(cx, problem, Some(applied.frozen_interfaces()),
                        old_h, dt_s, step_config)?;
                    let mut work = 0usize;
                    for row in &solved.newton.history {
                        poll(cx, work)?;
                        work = work.checked_add(row.linear_iterations)
                            .ok_or(EnthalpyError::InvalidInput("contact enthalpy inner work overflows"))?;
                    }
                    let iterations = solved.newton.iterations;
                    (solved, None, None, iterations, work)
                },
            };
            newton = newton.checked_add(trial_newton)
                .filter(|&n| n <= max_newton)
                .ok_or(EnthalpyError::InvalidInput("contact enthalpy Newton work exceeds allowance"))?;
            krylov = krylov.checked_add(trial_krylov)
                .filter(|&n| n <= max_krylov)
                .ok_or(EnthalpyError::InvalidInput("contact enthalpy Krylov work exceeds allowance"))?;
            let actual = self.evaluate(cx, problem.boundary, &solved.temperature)?;
            let physical = engine.evaluate_endpoint(
                cx, problem, Some(actual.frozen_interfaces()), old_h,
                &solved.specific_enthalpy_j_kg, dt_s, step_config, radiation.map(|(patches, _)| patches),
            )?;
            poll(cx, iteration)?;
            last_residual = physical.residual_norm_j;
            if last_residual.is_finite() && last_residual <= coupling.residual_tolerance_j
                && physical.energy_residual_j.is_finite()
                && physical.energy_residual_j.abs() <= step_config.energy_tolerance_j
            {
                return Ok(ContactEnthalpyStep {
                    conduction: solved, applied_contacts: applied, physical_contacts: actual,
                    physical, iterations: iteration, nonlinear_iterations: newton,
                    krylov_iterations: krylov,
                    radiation: report, combined_boundary: boundary, radiation_iterations,
                });
            }
            coupling.relax(cx, &mut guess, &solved.temperature)?;
        }
        Err(invalid("<contact-enthalpy-step>", format!(
            "temperature-dependent contact exhausted {} enthalpy responses; last physical residual {last_residual:e} J, target {:e} J",
            coupling.max_iterations, coupling.residual_tolerance_j,
        )).into())
    }
}


/// Rechecked heterogeneous enthalpy endpoint with the complete source-backed
/// contact Jacobian expressed in specific-enthalpy coordinates.
pub struct ContactEnthalpyLinearization<'m> {
    step: crate::transient::enthalpy::heterogeneous::HeterogeneousEnthalpyStepLinearization<'m>,
    contact: super::tangent::ContactTemperatureLinearization,
}

/// One complete physical transpose and named whole-law resistance controls.
#[derive(Debug)]
pub struct ContactEnthalpyGradient {
    /// Existing history, source and reference-density derivatives.
    pub storage: crate::transient::enthalpy::heterogeneous::HeterogeneousEnthalpyStepGradient,
    /// Independent multipliers of each complete matching resistance law at one.
    pub resistance_multipliers: Vec<crate::interface::ContactResistanceGradient>,
}

impl ContactEnthalpyLinearization<'_> {
    #[must_use]
    pub fn transport(
        &self,
    ) -> &crate::transient::enthalpy::heterogeneous::HeterogeneousEnthalpyStepLinearization<'_> {
        &self.step
    }

    #[must_use]
    pub fn contact(&self) -> &super::tangent::ContactTemperatureLinearization {
        &self.contact
    }

    pub fn temperature_pullback(
        &self,
        cx: &Cx<'_>,
        seed: &[f64],
    ) -> Result<Vec<f64>, crate::transient::enthalpy::adjoint::EnthalpyAdjointError> {
        self.step.transport().temperature_pullback(cx, seed)
    }

    pub fn pullback_with_resistance_controls(
        &self,
        cx: &Cx<'_>,
        seed: &[f64],
        linear: crate::LinearConfig,
        max_face_pairs: usize,
    ) -> Result<ContactEnthalpyGradient, crate::transient::enthalpy::adjoint::EnthalpyAdjointError> {
        poll(cx, 0)?;
        let storage = self.step.pullback(cx, seed, linear)?;
        let load_dual = self.step.transport().nodal_load_pullback(cx, &storage.transport.adjoint)?;
        let resistance_multipliers = self.contact.physical_contacts().frozen_interfaces()
            .matching_resistance_scale_pullback(
                cx, self.contact.temperature(), &load_dual, max_face_pairs,
            )?;
        Ok(ContactEnthalpyGradient { storage, resistance_multipliers })
    }
}

impl TemperatureContacts<'_, '_> {
    /// Revalidate and differentiate an accepted non-radiative enthalpy/contact
    /// endpoint without another primal solve. The outer contact iteration's
    /// original physical residual target remains authoritative.
    #[allow(clippy::too_many_arguments)]
    pub fn linearize_accepted_enthalpy<'m>(
        &self,
        cx: &Cx<'_>,
        engine: &HeterogeneousEnthalpyBackwardEuler<'m, '_>,
        problem: ConductionProblem<'_>,
        old_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        physical_residual_tolerance_j: f64,
        accepted: ContactEnthalpyStep,
        max_face_pairs: usize,
    ) -> Result<ContactEnthalpyLinearization<'m>, crate::transient::enthalpy::adjoint::EnthalpyAdjointError> {
        poll(cx, 0)?;
        if accepted.radiation.is_some() || accepted.combined_boundary.is_some() {
            return Err(crate::transient::enthalpy::adjoint::EnthalpyAdjointError::InvalidInput(
                "radiative temperature-contact enthalpy needs the radiative contact adjoint owner",
            ));
        }
        if !physical_residual_tolerance_j.is_finite() || physical_residual_tolerance_j <= 0.0 {
            return Err(crate::transient::enthalpy::adjoint::EnthalpyAdjointError::InvalidInput(
                "positive finite physical contact residual target required",
            ));
        }
        let contact = self.linearize_temperature(
            cx, problem.boundary, &accepted.conduction.temperature, max_face_pairs,
        )?;
        if contact.physical_contacts().states() != accepted.physical_contacts.states() {
            return Err(crate::transient::enthalpy::adjoint::EnthalpyAdjointError::InvalidInput(
                "rebound contact source state differs from the accepted physical endpoint",
            ));
        }
        let step = engine.linearize_accepted_with_target(
            cx, problem, Some(contact.physical_contacts().frozen_interfaces()),
            old_h, dt_s, config, accepted.conduction, physical_residual_tolerance_j,
        )?.with_contact_temperature_feedback(cx, &contact)?;
        Ok(ContactEnthalpyLinearization { step, contact })
    }
}


pub struct ContactEnthalpyRadiationLinearization<'m> {
    step: crate::transient::enthalpy::heterogeneous::HeterogeneousEnthalpyRadiationStepLinearization<'m>,
    contact: super::tangent::ContactTemperatureLinearization,
}

#[derive(Debug)]
pub struct ContactEnthalpyRadiationGradient {
    pub radiation: crate::transient::enthalpy::heterogeneous::HeterogeneousEnthalpyRadiationStepGradient,
    pub resistance_multipliers: Vec<crate::interface::ContactResistanceGradient>,
}

impl ContactEnthalpyRadiationLinearization<'_> {
    #[must_use]
    pub fn transport(
        &self,
    ) -> &crate::transient::enthalpy::heterogeneous::HeterogeneousEnthalpyRadiationStepLinearization<'_> {
        &self.step
    }

    #[must_use]
    pub fn contact(&self) -> &super::tangent::ContactTemperatureLinearization {
        &self.contact
    }

    pub fn temperature_pullback(
        &self,
        cx: &Cx<'_>,
        seed: &[f64],
    ) -> Result<Vec<f64>, crate::transient::enthalpy::adjoint::EnthalpyAdjointError> {
        self.step.transport().temperature_pullback(cx, seed)
    }
}


impl ContactEnthalpyRadiationLinearization<'_> {
    pub fn pullback_with_resistance_controls(
        &self,
        cx: &Cx<'_>,
        seed: &[f64],
        linear: crate::LinearConfig,
        max_face_pairs: usize,
    ) -> Result<ContactEnthalpyRadiationGradient, crate::transient::enthalpy::adjoint::EnthalpyAdjointError> {
        poll(cx, 0)?;
        let radiation = self.step.pullback(cx, seed, linear)?;
        let load_dual = self.step.transport().nodal_load_pullback(cx, &radiation.storage.transport.adjoint)?;
        let resistance_multipliers = self.contact.physical_contacts().frozen_interfaces()
            .matching_resistance_scale_pullback(
                cx, self.contact.temperature(), &load_dual, max_face_pairs,
            )?;
        Ok(ContactEnthalpyRadiationGradient { radiation, resistance_multipliers })
    }
}


impl TemperatureContacts<'_, '_> {
    #[allow(clippy::too_many_arguments)]
    pub fn linearize_accepted_enthalpy_with_ambient_radiation<'m>(
        &self,
        cx: &Cx<'_>,
        engine: &HeterogeneousEnthalpyBackwardEuler<'m, '_>,
        problem: ConductionProblem<'_>,
        old_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        physical_residual_tolerance_j: f64,
        patches: &[AmbientRadiationPatch],
        accepted: ContactEnthalpyStep,
        max_feedback_entries: usize,
        max_face_pairs: usize,
    ) -> Result<ContactEnthalpyRadiationLinearization<'m>, crate::transient::enthalpy::adjoint::EnthalpyAdjointError> {
        poll(cx, 0)?;
        if accepted.radiation.is_none() || accepted.combined_boundary.is_none() {
            return Err(crate::transient::enthalpy::adjoint::EnthalpyAdjointError::InvalidInput(
                "radiative contact adjoint requires an accepted radiative-contact endpoint",
            ));
        }
        if !physical_residual_tolerance_j.is_finite() || physical_residual_tolerance_j <= 0.0 {
            return Err(crate::transient::enthalpy::adjoint::EnthalpyAdjointError::InvalidInput(
                "positive finite physical contact/radiation residual target required",
            ));
        }
        let contact = self.linearize_temperature(
            cx, problem.boundary, &accepted.conduction.temperature, max_face_pairs,
        )?;
        if contact.physical_contacts().states() != accepted.physical_contacts.states() {
            return Err(crate::transient::enthalpy::adjoint::EnthalpyAdjointError::InvalidInput(
                "rebound radiative-contact source state differs from accepted physical endpoint",
            ));
        }
        let step = engine
            .linearize_accepted_parameters_with_ambient_radiation_target(
                cx,
                problem,
                Some(contact.physical_contacts().frozen_interfaces()),
                old_h,
                dt_s,
                config,
                patches,
                accepted.conduction,
                max_feedback_entries,
                physical_residual_tolerance_j,
            )?
            .with_contact_temperature_feedback(cx, &contact)?;
        Ok(ContactEnthalpyRadiationLinearization { step, contact })
    }
}
