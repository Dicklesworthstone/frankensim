//! Nonlinear backward Euler with immutable physical history.
//!
//! Solve F(T) = C(T - T_old) + dt (A(T)T - b) = 0. The tangent is
//! C + dt J_steady(T), including the existing element K'(T) contribution.
//! Newton trials never update T_old, and every merit evaluation reassembles
//! the actual constitutive law. A small update is NOT a convergence test.

use super::*;
use crate::assemble::{AssembledSystem, assemble_jacobian_with_optional_interfaces};
use crate::solve::LineSearch;
use fs_solver::{FgmresState, LinearOp};
use crate::PointwiseAmbientRadiation;
use super::pointwise::RadiativeDirection;
use crate::interface::temperature::{TemperatureContactEvaluation, TemperatureContacts};

/// Explicit nonlinear work, stopping and globalization policy for one endpoint.
/// Heat capacity remains temperature independent. The ordinary endpoint keeps
/// its supplied contact operator fixed; the source-backed contact entry point
/// includes the original temperature-dependent contact law in the same solve.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NonlinearStepConfig {
    /// Maximum accepted Newton updates, not a convergence condition.
    pub max_iterations: usize,
    /// Relative tolerance against the INITIAL free transient residual norm.
    /// The scale is fixed during this solve; it is not the absolute-temperature
    /// load C*T_old, which could conceal a poor temperature correction.
    pub residual_rtol: f64,
    /// Absolute Euclidean free-residual floor in joules.
    pub residual_atol_j: f64,
    /// Armijo policy on the actual transient residual. Invalid material-range
    /// trials consume this same budget instead of extrapolating the curve.
    pub line_search: LineSearch,
}

impl Default for NonlinearStepConfig {
    fn default() -> Self {
        Self {
            max_iterations: 32,
            residual_rtol: 1.0e-10,
            residual_atol_j: 1.0e-10,
            line_search: LineSearch::default(),
        }
    }
}

impl NonlinearStepConfig {
    /// Admit the policy before allocation or assembly.
    ///
    /// # Errors
    /// Refuses zero update budget, non-finite/negative tolerances, relative
    /// tolerance outside (0,1), or invalid Armijo/shrink factors.
    pub fn validate(&self) -> Result<(), ConductionError> {
        if self.max_iterations == 0
            || !(self.residual_rtol.is_finite() && self.residual_rtol > 0.0 && self.residual_rtol < 1.0)
            || !(self.residual_atol_j.is_finite() && self.residual_atol_j >= 0.0)
            || !(self.line_search.armijo_c.is_finite()
                && self.line_search.armijo_c > 0.0 && self.line_search.armijo_c < 1.0)
            || !(self.line_search.shrink.is_finite()
                && self.line_search.shrink > 0.0 && self.line_search.shrink < 1.0)
        {
            return Err(invalid("nonlinear endpoint requires a positive Newton budget, rtol in (0,1), nonnegative joule tolerance and Armijo/shrink in (0,1)"));
        }
        Ok(())
    }
}

/// A residual- AND energy-admitted endpoint with its nonlinear solve evidence.
#[derive(Debug, Clone)]
pub struct NonlinearStepSolution {
    /// Same endpoint, flux and discrete-energy vocabulary as the linear step.
    /// Its Krylov count is the sum over all Newton corrections; its relative
    /// residual is the worst recomputed INNER linear-solve relative residual.
    pub step: StepSolution,
    /// Accepted Newton updates (zero if the initial state already passed).
    pub nonlinear_iterations: usize,
    /// Rejected Newton trials, including material-range refusals.
    pub backtracks: usize,
    /// Initial free residual norm, joules. With prescribed endpoint data,
    /// the predictor has old FREE temperatures and NEW fixed temperatures;
    /// physical history in the capacity term remains the full old field.
    pub initial_residual_j: f64,
    /// Recomputed norm of F(T_new) on the published field, joules.
    pub residual_j: f64,
    /// Fixed atol + rtol * initial_residual_j gate, joules. `residual_j`
    /// exceeds it only when a contact Newton stopped at the field's
    /// resolution (`advance_nonlinear_prescribed_with_temperature_contacts`).
    pub threshold_j: f64,
}

/// A nonlinear backward-Euler endpoint whose source-backed contact law was
/// evaluated and differentiated at every accepted Newton field.
#[derive(Debug)]
pub struct ContactNonlinearStepSolution {
    /// Original nonlinear endpoint evidence and physical thermal state.
    pub step: NonlinearStepSolution,
    /// Contact source states at the exact accepted endpoint.
    pub contacts: TemperatureContactEvaluation,
}

pub(super) struct Evaluation {
    pub(super) system: AssembledSystem,
    pub(super) residual: Vec<f64>,
    pub(super) norm: f64,
}

impl BackwardEuler<'_> {
    /// Advance one nonlinear endpoint from immutable old temperatures.
    ///
    /// Uses the same mesh, lumped capacity, heterogeneous material laws,
    /// boundaries and contact operators as `advance`. The NEW endpoint's k(T)
    /// is evaluated on every trial; conductivity is never frozen at T_old.
    /// Prescribed temperatures must already equal history. Constant models are
    /// also admitted, making the linear path an independent comparison seam.
    ///
    /// `config.linear.max_iterations` is a TOTAL inner-iteration cap across
    /// all Newton corrections of this endpoint. Restarts are capped at 64 and
    /// shortened for the remaining budget. Cancellation is checked at Newton,
    /// line-search and Krylov-cycle boundaries, as well as in assembly.
    ///
    /// # Errors
    /// Invalid data, sampled-material extrapolation, exhausted Newton/Krylov or
    /// line-search budgets, non-finite arithmetic, failed joule balance, or
    /// cancellation return no new physical state. The caller's history is never
    /// modified, even after several accepted internal Newton trials.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_nonlinear(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>, old: &[f64], dt_s: f64,
        config: StepConfig, nonlinear: NonlinearStepConfig,
    ) -> Result<NonlinearStepSolution, ConductionError> {
        self.advance_nonlinear_with_boundary(cx, problem, interfaces, old, dt_s,
            config, nonlinear, false, None)
    }

    /// Advance k(T) conduction with explicitly prescribed ENDPOINT temperatures.
    ///
    /// The full `old` field remains immutable physical history. New Dirichlet
    /// values come from `problem.boundary` and are installed in the predictor
    /// before the first constitutive/residual evaluation. Newton and line-search
    /// trials change only free nodes, so conductivity always sees the actual
    /// endpoint field, including its new boundary data.
    ///
    /// Storage and prescribed reaction include the boundary-node temperature
    /// changes. As with [`Self::advance_prescribed`], the caller supplies the
    /// time grid and endpoint data; no continuous profile or instantaneous heat
    /// impulse is inferred. Capacity and contact resistance remain independent
    /// of temperature. Existing residual, energy and total-work gates apply.
    ///
    /// # Errors
    /// The refusals from [`Self::advance_nonlinear`], except that changed fixed
    /// temperatures are explicitly admitted. An out-of-range predictor refuses
    /// rather than extrapolating its material law. Failure returns no new field.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_nonlinear_prescribed(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>, old: &[f64], dt_s: f64,
        config: StepConfig, nonlinear: NonlinearStepConfig,
    ) -> Result<NonlinearStepSolution, ConductionError> {
        self.advance_nonlinear_with_boundary(cx, problem, interfaces, old, dt_s,
            config, nonlinear, true, None)
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn advance_nonlinear_with_boundary(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>, old: &[f64], dt_s: f64,
        config: StepConfig, nonlinear: NonlinearStepConfig, prescribed: bool,
        pointwise: Option<&PointwiseAmbientRadiation<'_>>,
    ) -> Result<NonlinearStepSolution, ConductionError> {
        poll(cx, 0)?;
        nonlinear.validate()?;
        if config.linear.restart == 0 {
            return Err(invalid("nonlinear endpoint requires a positive FGMRES restart"));
        }
        let dofs = if prescribed {
            self.admit_endpoint_step(cx, problem, old, dt_s, config)?
        } else {
            self.admit_step(cx, problem, old, dt_s, config)?
        };
        if pointwise.is_some_and(|radiation| !radiation.matches(problem.mesh, problem.boundary)) {
            return Err(invalid("pointwise radiation must retain the exact capacity mesh and physical boundary"));
        }
        let mut temperature = old.to_vec();
        if prescribed {
            for (slot, &vertex) in dofs.fixed().iter().enumerate() {
                if slot % 512 == 0 { poll(cx, slot)?; }
                temperature[vertex] = dofs.prescribed()[vertex];
            }
        }
        let mut ev = self.evaluate_nonlinear_endpoint(cx, problem, interfaces, old,
            &temperature, dt_s, &dofs, pointwise)?;
        let initial_residual_j = ev.norm;
        let threshold_j = finite(nonlinear.residual_atol_j + nonlinear.residual_rtol * initial_residual_j)?;
        let mut updates = 0;
        let mut work = 0;
        let mut backtracks = 0_usize;
        let mut worst_linear_residual = 0.0_f64;
        loop {
            poll(cx, updates)?;
            if ev.norm <= threshold_j {
                let residual_j = ev.norm;
                let radiation = pointwise.map(|owner| owner.load_and_report(cx, &temperature)).transpose()?;
                let step = self.finish_step_with_external_heat(cx, problem, interfaces, old, dt_s, config,
                    &dofs, &ev.system, temperature, worst_linear_residual, work,
                    radiation.as_ref().map(|(load, report)| (load.as_slice(), report.outward_heat_w)))?;
                return Ok(NonlinearStepSolution { step, nonlinear_iterations: updates,
                    backtracks, initial_residual_j, residual_j, threshold_j });
            }
            if updates == nonlinear.max_iterations {
                return Err(ConductionError::NotConverged {
                    iterations: updates, residual: ev.norm, threshold: threshold_j,
                });
            }
            if work == config.linear.max_iterations {
                return Err(ConductionError::LinearSolveFailed {
                    iteration: updates, krylov_iterations: work,
                    true_relative_residual: 1.0, tolerance: config.linear.tolerance,
                });
            }
            let jacobian = assemble_jacobian_with_optional_interfaces(cx, self.mesh,
                problem.boundary, problem.material, &temperature, interfaces, problem.element_materials)?;
            let full = axpy_csr(&self.capacity, 1.0, &jacobian, dt_s);
            // Newton increments on fixed nodes are ZERO: no absolute lift.
            let (matrix, _) = reduce_matrix_and_lift(&full, &dofs);
            // Precondition with the SPD Picard block, not with a nonsymmetric
            // K'(T) Jacobian masquerading as an SPD operator.
            let picard = axpy_csr(&self.capacity, 1.0, &ev.system.operator, dt_s);
            let (picard, _) = reduce_matrix_and_lift(&picard, &dofs);
            let rhs: Vec<_> = ev.residual.iter().map(|&r| -r).collect();
            let remaining = LinearConfig {
                max_iterations: config.linear.max_iterations - work, ..config.linear
            };
            let radiation = pointwise.map(|owner| owner.assemble(cx, &temperature)).transpose()?;
            let (direction, relative, iters) = solve_direction(cx, &matrix, &picard, &rhs, remaining, updates,
                radiation.as_ref().map(|radiation| (radiation, &dofs, dt_s)))?;
            work += iters;
            worst_linear_residual = worst_linear_residual.max(relative);
            let mut alpha = 1.0;
            let mut rejected = 0;
            loop {
                poll(cx, updates)?;
                let mut trial = temperature.clone();
                for (index, &vertex) in dofs.free().iter().enumerate() {
                    if index % 512 == 0 { poll(cx, index)?; }
                    trial[vertex] = temperature[vertex] + alpha * direction[index];
                }
                let candidate = if trial.iter().all(|t| t.is_finite()
                    && (pointwise.is_none() || *t > 0.0)) {
                    match self.evaluate_nonlinear_endpoint(cx, problem, interfaces, old,
                        &trial, dt_s, &dofs, pointwise) {
                        Ok(candidate) => Some(candidate),
                        Err(ConductionError::OutsideTemperatureSpan { .. }) => None,
                        Err(ConductionError::MaterialQuery { ref property, .. })
                            if pointwise.is_some() && property == crate::SURFACE_EMISSIVITY_PROPERTY => None,
                        Err(error) => return Err(error),
                    }
                } else { None };
                if let Some(candidate) = candidate {
                    let sufficient = (1.0 - nonlinear.line_search.armijo_c * alpha) * ev.norm;
                    if candidate.norm <= threshold_j || candidate.norm <= sufficient {
                        temperature = trial;
                        ev = candidate;
                        break;
                    }
                }
                if rejected == nonlinear.line_search.max_backtracks {
                    return Err(ConductionError::LineSearchFailed {
                        iteration: updates, backtracks: rejected, smallest_step: alpha,
                    });
                }
                rejected += 1;
                backtracks = backtracks.checked_add(1)
                    .ok_or_else(|| invalid("nonlinear backtrack count overflow"))?;
                alpha *= nonlinear.line_search.shrink;
            }
            updates += 1;
        }
    }

    /// Advance one prescribed endpoint with source-backed contact feedback
    /// inside the same Newton solve as storage and nonlinear bulk conductivity.
    ///
    /// Every residual trial freshly evaluates the pinned contact law at the
    /// trial temperature. The tangent adds only dK_contact/dT*T to the ordinary
    /// bulk-plus-frozen-contact Jacobian, so the physical contact block is
    /// counted exactly once. Source-domain trial misses backtrack through the
    /// existing Armijo policy; no claim switching or extrapolation is allowed.
    /// Continuous source curves use their original primal tangent at knots and
    /// inward domain endpoints. Surface laws retain their unique-gradient
    /// requirement; accepted-endpoint derivatives and adjoints still require
    /// the separate strict two-sided contact linearization.
    ///
    /// The accepted endpoint is returned only after the original nonlinear
    /// residual and energy gates pass, with one exception for the residual
    /// gate: a final line search rejected with its smallest trial equal to the
    /// field bit for bit. The Newton step then moves no free temperature (the
    /// field's floating-point resolution along the Newton direction), so that
    /// field goes to the unchanged energy gate, and its `residual_j` stays
    /// above `threshold_j` in the report. Every other stall refuses. Callers
    /// keep their own physical gates. Geometry, pressure and all nonthermal
    /// source coordinates remain fixed for the step.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn advance_nonlinear_prescribed_with_temperature_contacts(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        contacts: &TemperatureContacts<'_, '_>,
        old: &[f64],
        dt_s: f64,
        config: StepConfig,
        nonlinear: NonlinearStepConfig,
        max_face_pairs: usize,
    ) -> Result<ContactNonlinearStepSolution, ConductionError> {
        poll(cx, 0)?;
        contacts.check_mesh(problem.mesh)?;
        nonlinear.validate()?;
        if max_face_pairs == 0 || config.linear.restart == 0 {
            return Err(invalid("temperature-contact Newton requires positive face and FGMRES restart budgets"));
        }
        let dofs = self.admit_endpoint_step(cx, problem, old, dt_s, config)?;
        let mut temperature = old.to_vec();
        for (slot, &vertex) in dofs.fixed().iter().enumerate() {
            if slot % 512 == 0 { poll(cx, slot)?; }
            temperature[vertex] = dofs.prescribed()[vertex];
        }
        let mut applied = contacts.evaluate(cx, problem.boundary, &temperature)?;
        let mut ev = self.evaluate_endpoint(
            cx, problem, Some(applied.frozen_interfaces()), old, &temperature, dt_s, &dofs,
        )?;
        let initial_residual_j = ev.norm;
        let threshold_j = finite(
            nonlinear.residual_atol_j + nonlinear.residual_rtol * initial_residual_j,
        )?;
        let mut updates = 0usize;
        let mut work = 0usize;
        let mut backtracks = 0usize;
        let mut worst_linear_residual = 0.0_f64;
        let mut at_resolution = false;
        loop {
            poll(cx, updates)?;
            if ev.norm <= threshold_j || at_resolution {
                let residual_j = ev.norm;
                let step = self.finish_step(
                    cx, problem, Some(applied.frozen_interfaces()), old, dt_s, config,
                    &dofs, &ev.system, temperature, worst_linear_residual, work,
                )?;
                return Ok(ContactNonlinearStepSolution {
                    step: NonlinearStepSolution {
                        step,
                        nonlinear_iterations: updates,
                        backtracks,
                        initial_residual_j,
                        residual_j,
                        threshold_j,
                    },
                    contacts: applied,
                });
            }
            if updates == nonlinear.max_iterations {
                return Err(ConductionError::NotConverged {
                    iterations: updates, residual: ev.norm, threshold: threshold_j,
                });
            }
            if work == config.linear.max_iterations {
                return Err(ConductionError::LinearSolveFailed {
                    iteration: updates, krylov_iterations: work,
                    true_relative_residual: 1.0, tolerance: config.linear.tolerance,
                });
            }

            let contact = contacts.linearize_temperature_for_newton(
                cx, problem.boundary, &temperature, max_face_pairs,
            )?;
            let bulk = assemble_jacobian_with_optional_interfaces(
                cx, self.mesh, problem.boundary, problem.material, &temperature,
                Some(contact.physical_contacts().frozen_interfaces()),
                problem.element_materials,
            )?;
            // bulk already contains K_contact. feedback is only dK_contact/dT*T.
            let complete = axpy_csr(&bulk, 1.0, contact.feedback(), 1.0);
            let full = axpy_csr(&self.capacity, 1.0, &complete, dt_s);
            let (matrix, _) = reduce_matrix_and_lift(&full, &dofs);
            let picard = axpy_csr(&self.capacity, 1.0, &ev.system.operator, dt_s);
            let (picard, _) = reduce_matrix_and_lift(&picard, &dofs);
            let rhs: Vec<_> = ev.residual.iter().map(|&r| -r).collect();
            let remaining = LinearConfig {
                max_iterations: config.linear.max_iterations - work, ..config.linear
            };
            let (direction, relative, iters) =
                solve_direction(cx, &matrix, &picard, &rhs, remaining, updates, None)?;
            work += iters;
            worst_linear_residual = worst_linear_residual.max(relative);

            let mut alpha = 1.0;
            let mut rejected = 0usize;
            loop {
                poll(cx, updates)?;
                let mut trial = temperature.clone();
                for (index, &vertex) in dofs.free().iter().enumerate() {
                    if index % 512 == 0 { poll(cx, index)?; }
                    trial[vertex] = temperature[vertex] + alpha * direction[index];
                }
                let candidate = if trial.iter().all(|t| t.is_finite() && *t > 0.0) {
                    let mut source_domain = false;
                    let contact_trial = match contacts.evaluate_with_query_refusal(
                        cx, problem.boundary, &trial, &mut |error| {
                            source_domain = matches!(
                                error,
                                fs_matdb::MatDbError::NoClaimInDomain { .. }
                                    | fs_matdb::MatDbError::PinnedClaimOutOfDomain { .. }
                                    | fs_matdb::MatDbError::OutsideKnotSpan { .. }
                                    | fs_matdb::MatDbError::OutsideSurfaceSupport { .. }
                            );
                        },
                    ) {
                        Ok(value) => Some(value),
                        Err(_) if source_domain => None,
                        Err(error) => return Err(error),
                    };
                    match contact_trial {
                        Some(contact_trial) => match self.evaluate_endpoint(
                            cx, problem, Some(contact_trial.frozen_interfaces()), old,
                            &trial, dt_s, &dofs,
                        ) {
                            Ok(candidate) => Some((candidate, contact_trial)),
                            Err(ConductionError::OutsideTemperatureSpan { .. }) => None,
                            Err(error) => return Err(error),
                        },
                        None => None,
                    }
                } else {
                    None
                };
                if let Some((candidate, contact_trial)) = candidate {
                    let sufficient =
                        (1.0 - nonlinear.line_search.armijo_c * alpha) * ev.norm;
                    if candidate.norm <= threshold_j || candidate.norm <= sufficient {
                        temperature = trial;
                        ev = candidate;
                        applied = contact_trial;
                        break;
                    }
                }
                if rejected == nonlinear.line_search.max_backtracks {
                    // The line search evaluates its smallest trial last. When
                    // that trial is the field to the bit, the Newton step moves
                    // no free temperature: the field is at its floating-point
                    // resolution along the Newton direction, where a target
                    // below the residual's rounding floor cannot be met. That
                    // field goes to the unchanged energy gate in finish_step.
                    at_resolution = dofs.free().iter()
                        .all(|&vertex| trial[vertex].to_bits() == temperature[vertex].to_bits());
                    if at_resolution {
                        break;
                    }
                    return Err(ConductionError::LineSearchFailed {
                        iteration: updates, backtracks: rejected, smallest_step: alpha,
                    });
                }
                rejected += 1;
                backtracks = backtracks.checked_add(1)
                    .ok_or_else(|| invalid("nonlinear contact backtrack count overflow"))?;
                alpha *= nonlinear.line_search.shrink;
            }
            if !at_resolution {
                updates += 1;
            }
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn evaluate_endpoint(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>, old: &[f64], temperature: &[f64],
        dt: f64, dofs: &DofMap,
    ) -> Result<Evaluation, ConductionError> {
        poll(cx, 0)?;
        self.validate_capacity_temperature(cx, temperature)?;
        let system = assemble_operator_scaled_with_interfaces(cx, self.mesh, problem.boundary,
            problem.material, problem.source, temperature, None, interfaces, problem.element_materials)?;
        let delta = temperature.iter().zip(old).map(|(t, old)| finite(t - old))
            .collect::<Result<Vec<_>, _>>()?;
        let mut storage = vec![0.0; temperature.len()];
        self.capacity.spmv(&delta, &mut storage);
        poll(cx, 0)?;
        let residual = dofs.free().iter().map(|&v|
            finite(storage[v] + dt * finite(system.heat_residual_row(v, temperature))?))
            .collect::<Result<Vec<_>, _>>()?;
        poll(cx, 0)?;
        let norm = finite(norm2(&residual))?;
        Ok(Evaluation { system, residual, norm })
    }
}

fn solve_direction(
    cx: &Cx<'_>, matrix: &Csr, picard: &Csr, rhs: &[f64],
    config: LinearConfig, iteration: usize,
    radiation: Option<(&crate::PointwiseRadiationAssembly, &DofMap, f64)>,
) -> Result<(Vec<f64>, f64, usize), ConductionError> {
    poll(cx, iteration)?;
    let scale = rhs.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
    if scale == 0.0 { return Ok((vec![0.0; rhs.len()], 0.0, 0)); }
    let normalized: Vec<_> = rhs.iter().map(|v| v / scale).collect();
    let op = RadiativeDirection::new(cx, CsrOp::general(matrix.clone()), radiation);
    let pre = spd_preconditioner(picard);
    let restart = config.restart.min(rhs.len()).min(64);
    let mut state = FgmresState::new(&normalized, restart);
    while state.rel_residual() >= config.tolerance && state.iters < config.max_iterations {
        poll(cx, state.iters)?;
        let before = state.iters;
        state.restart = restart.min(config.max_iterations - before);
        state.run(&op, &pre, &normalized, config.tolerance, 1);
        op.check()?;
        if state.iters == before { break; }
    }
    poll(cx, state.iters)?;
    let mut applied = vec![0.0; rhs.len()];
    op.apply(&state.x, &mut applied);
    op.check()?;
    let residual = normalized.iter().zip(applied).map(|(b, a)| finite(b - a))
        .collect::<Result<Vec<_>, _>>()?;
    let relative = finite(norm2(&residual) / norm2(&normalized))?;
    if relative >= config.tolerance {
        return Err(ConductionError::LinearSolveFailed { iteration,
            krylov_iterations: state.iters, true_relative_residual: relative,
            tolerance: config.tolerance });
    }
    let direction = state.x.iter().map(|v| finite(v * scale)).collect::<Result<_, _>>()?;
    Ok((direction, relative, state.iters))
}
