//! Steady source-backed contact coupling through the original spatial solver.

/// Revalidated steady tangents and adjoints of the complete contact law.
pub mod linearization;
pub use linearization::{ContactSteadyGradient, ContactSteadyLinearization};
/// Coupled sparse Newton iteration for strong source-law temperature feedback.
pub mod newton;
pub use newton::{ContactSteadyNewtonSolution, SteadyContactNewtonConfig};

use fs_exec::Cx;
use fs_solver::norm2;
use crate::{ConductionError, ConductionProblem, ConductionSolution, DofMap, EnergyBalance,
    InitialGuess, RobinFlux, SolveConfig, solve_with_interfaces};
use crate::assemble::{assemble_operator_scaled_with_interfaces, full_residual};
use crate::solve::energy_balance;
use super::{TemperatureContactEvaluation, TemperatureContacts, finite, invalid, poll, positive};

/// Bounded contact iteration with physical steady tolerances, in watts.
#[derive(Debug, Clone, Copy)]
pub struct SteadyContactConfig {
    /// Maximum complete spatial thermal responses, including discarded trials.
    pub max_iterations: usize,
    /// Convex update fraction of the next source-query field, in (0,1].
    pub relaxation: f64,
    /// Absolute Euclidean free-node residual target, W.
    pub residual_tolerance_w: f64,
    /// Absolute whole-domain heat-balance closure allowance, W.
    pub energy_tolerance_w: f64,
}

/// One accepted field, its counted solve, and freshly evaluated physical contacts.
#[derive(Debug)]
pub struct ContactSteadySolution {
    /// Last counted spatial response. Its report uses applied_contacts.
    /// The physical source law and its acceptance evidence are retained below.
    pub conduction: ConductionSolution,
    /// Contact source states actually used by that spatial thermal response.
    pub applied_contacts: TemperatureContactEvaluation,
    /// Source states freshly queried at the unchanged returned temperatures.
    pub physical_contacts: TemperatureContactEvaluation,
    /// Physical (vertex, outward residual W) pairs in ascending free-node order.
    pub free_residual_w: Vec<(usize, f64)>,
    /// Euclidean norm of the physical free-node residual, W.
    pub residual_norm_w: f64,
    /// Physical source/Neumann/Robin/Dirichlet balance from the existing owner.
    pub energy: EnergyBalance,
    /// Physical Robin heat reports from the existing boundary integration owner.
    pub robin_fluxes: Vec<RobinFlux>,
    /// (vertex, reservoir heat entering W) in ascending Dirichlet-node order.
    pub dirichlet_reactions_w: Vec<(usize, f64)>,
    /// Complete thermal responses, including all discarded contact trials.
    pub iterations: usize,
    /// Sum of nonlinear updates in every counted spatial response.
    pub nonlinear_iterations: usize,
    /// Sum of actual Krylov iterations in every counted spatial response.
    pub krylov_iterations: usize,
}

impl SteadyContactConfig {
    fn admit(self, solver: &SolveConfig) -> Result<(usize, usize), ConductionError> {
        if self.max_iterations == 0 || solver.stop.max_iterations == 0
            || solver.linear.max_iterations == 0 || solver.linear.restart == 0
            || !(self.relaxation.is_finite() && self.relaxation > 0.0 && self.relaxation <= 1.0)
            || !(self.residual_tolerance_w.is_finite() && self.residual_tolerance_w > 0.0)
            || !(self.energy_tolerance_w.is_finite() && self.energy_tolerance_w > 0.0)
        {
            return Err(invalid("<steady-contact>", "positive steady work and watt residual/heat tolerances and relaxation in (0,1] are required"));
        }
        let nonlinear = self.max_iterations.checked_mul(solver.stop.max_iterations)
            .ok_or_else(|| invalid("<steady-contact>", "combined nonlinear allowance overflows"))?;
        // The original FGMRES owner rounds its limit to a full restart cycle.
        // This also bounds the original CG path without changing either solver.
        let per_linear = solver.linear.max_iterations.div_ceil(solver.linear.restart)
            .checked_mul(solver.linear.restart)
            .ok_or_else(|| invalid("<steady-contact>", "Krylov cycle allowance overflows"))?;
        let krylov = nonlinear.checked_mul(per_linear)
            .ok_or_else(|| invalid("<steady-contact>", "combined Krylov allowance overflows"))?;
        Ok((nonlinear, krylov))
    }
}

impl TemperatureContacts<'_, '_> {
    /// Solve steady conduction with fresh local R(T_A,T_B) or declared mean laws.
    ///
    /// Matching laws retain their declared face-mean or pointwise convention;
    /// nonmatching laws retain common-refinement quadrature. Every trial delegates the spatial
    /// problem, heterogeneous material assignment and nonlinear bulk k(T) to
    /// the original steady solver. Only the next query field is relaxed.
    ///
    /// The initial argument is a complete positive Kelvin field. Prescribed
    /// boundary values replace its fixed entries before the first source query;
    /// solver.initial is superseded by this explicit field and subsequent
    /// relaxed warm starts. Pure-Neumann and fully prescribed limitations
    /// remain those of the existing steady solver.
    ///
    /// Acceptance reassembles the full physical operator at the unchanged
    /// candidate, using freshly queried contacts and candidate bulk response.
    /// Both free-node residual and independent heat-accounting gates must pass
    /// their original watt tolerances. No global cancellation, inner step-size
    /// stop, hidden final solve or enlarged tolerance can admit a false balance.
    ///
    /// Source support and selection remain pinned. Any source, cancellation,
    /// inner-solve or coupling refusal returns no field and mutates no caller
    /// state. This bounded fixed point does not assert uniqueness or convergence
    /// for arbitrary positive source laws. Use linearize_steady_endpoint for
    /// a separately revalidated complete steady tangent and adjoint.
    pub fn solve_steady(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>, initial: &[f64],
        mut solver: SolveConfig, coupling: SteadyContactConfig,
    ) -> Result<ContactSteadySolution, ConductionError> {
        poll(cx, 0)?;
        self.check_mesh(problem.mesh)?;
        let (max_nonlinear, max_krylov) = coupling.admit(&solver)?;
        if initial.len() != problem.mesh.vertex_count() {
            return Err(invalid("<steady-contact>", "initial field and contact mesh differ"));
        }
        for (i, &temperature) in initial.iter().enumerate() {
            if i % 256 == 0 { poll(cx, i)?; }
            positive(temperature)?;
        }
        let dofs = DofMap::new(problem.boundary, problem.mesh.vertex_count())?;
        let mut guess = dofs.scatter(&dofs.free().iter().map(|&i| initial[i]).collect::<Vec<_>>());
        // Avoid cloning an unrelated caller-supplied initial vector every pass.
        solver.initial = InitialGuess::Free(Vec::new());
        let mut nonlinear_iterations = 0usize;
        let mut krylov_iterations = 0usize;
        let mut last_residual = f64::INFINITY;
        let mut last_closure = f64::INFINITY;
        for iteration in 1..=coupling.max_iterations {
            poll(cx, iteration)?;
            let applied = self.evaluate(cx, problem.boundary, &guess)?;
            let mut trial = solver.clone();
            trial.initial = InitialGuess::Free(dofs.free().iter().map(|&i| guess[i]).collect());
            let conduction = solve_with_interfaces(cx, problem, applied.frozen_interfaces(), trial)?;
            nonlinear_iterations = nonlinear_iterations.checked_add(conduction.report.iterations)
                .filter(|&n| n <= max_nonlinear)
                .ok_or_else(|| invalid("<steady-contact>", "counted nonlinear work exceeds its combined allowance"))?;
            for record in &conduction.report.linear {
                krylov_iterations = krylov_iterations.checked_add(record.iterations)
                    .filter(|&n| n <= max_krylov)
                    .ok_or_else(|| invalid("<steady-contact>", "counted Krylov work exceeds its combined allowance"))?;
            }
            let physical = self.evaluate(cx, problem.boundary, &conduction.temperature)?;
            let system = assemble_operator_scaled_with_interfaces(
                cx, problem.mesh, problem.boundary, problem.material, problem.source,
                &conduction.temperature, None, Some(physical.frozen_interfaces()), problem.element_materials,
            )?;
            let residual = full_residual(&system, &conduction.temperature);
            for (i, &value) in residual.iter().enumerate() {
                if i % 256 == 0 { poll(cx, i)?; }
                finite(value)?;
            }
            let free: Vec<_> = dofs.free().iter().map(|&i| residual[i]).collect();
            last_residual = finite(norm2(&free))?;
            let (energy, robin_fluxes) = energy_balance(
                problem.mesh, problem.boundary, problem.source, &system, &dofs, &conduction.temperature,
            );
            for value in [energy.source_w, energy.neumann_out_w, energy.robin_out_w,
                energy.dirichlet_in_w, energy.closure_w, energy.scale_w] { finite(value)?; }
            last_closure = energy.closure_w;
            poll(cx, iteration)?;
            if last_residual <= coupling.residual_tolerance_w
                && last_closure.abs() <= coupling.energy_tolerance_w
            {
                return Ok(ContactSteadySolution {
                    conduction, applied_contacts: applied, physical_contacts: physical,
                    free_residual_w: dofs.free().iter().copied().zip(free).collect(),
                    residual_norm_w: last_residual, energy, robin_fluxes,
                    dirichlet_reactions_w: dofs.fixed().iter().map(|&i| (i, residual[i])).collect(),
                    iterations: iteration, nonlinear_iterations, krylov_iterations,
                });
            }
            // Prescribed entries keep their exact boundary values; only free
            // temperatures belong to the relaxed constitutive query field.
            for (slot, &vertex) in dofs.free().iter().enumerate() {
                if slot % 256 == 0 { poll(cx, slot)?; }
                guess[vertex] = positive((1.0-coupling.relaxation)*guess[vertex]
                    + coupling.relaxation*conduction.temperature[vertex])?;
            }
        }
        Err(invalid("<steady-contact>", format!(
            "source contact exhausted {} thermal responses; physical residual {last_residual:e} W (target {:e}), heat closure {last_closure:e} W (allowance {:e})",
            coupling.max_iterations, coupling.residual_tolerance_w, coupling.energy_tolerance_w,
        )))
    }
}
