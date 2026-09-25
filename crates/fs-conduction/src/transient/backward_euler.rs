//! Immutable-history backward-Euler steps for partitioned thermal coupling.
//!
//! [`BackwardEuler::advance`] solves a linear correction from the SAME old
//! field. [`BackwardEuler::advance_nonlinear`] solves the endpoint k(T) law
//! with Newton/FGMRES. A caller may iterate Robin references without advancing
//! physical time repeatedly. Only its accepted field becomes the next history.
//!
//! Both paths share exact P1 capacity, steady assembly, contact operators,
//! boundary integrals and discrete energy accounting. Heat capacity and contact
//! resistance are temperature independent. [`StepLinearization`] adds discrete
//! endpoint/history derivatives; fluid storage is not introduced.

mod nonlinear;
pub use nonlinear::{NonlinearStepConfig, NonlinearStepSolution};
mod adjoint;
pub use adjoint::StepLinearization;

use fs_exec::Cx;
use fs_solver::{CheckedCgConfig, CheckedCgError, CsrOp, checked_cg, norm2};
use fs_sparse::Csr;

use crate::assemble::{AssembledSystem, assemble_operator_scaled_with_interfaces, reduce_matrix_and_lift, DofMap};
use crate::solve::{energy_balance, spd_preconditioner};
use crate::{ConductionError, ConductionMesh, ConductionProblem, InterfaceFlux, LinearConfig,
    RobinFlux, ThermalInterfaces};
use super::{assemble_capacitance, assemble_capacitance_from, axpy_csr, VolumetricHeatCapacity};

/// Work and absolute discrete energy-closure budgets for one solid step.
#[derive(Debug, Clone, Copy)]
pub struct StepConfig {
    /// Krylov tolerance/budget. For the nonlinear entry point the iteration
    /// cap covers ALL Newton corrections, not each correction separately.
    /// The linear path permits at most two residual-defect repairs within
    /// this same iteration cap; neither the tolerance nor the energy gate widens.
    pub linear: LinearConfig,
    /// Maximum absolute storage-minus-net-input mismatch in joules.
    pub energy_tolerance_j: f64,
}

/// Accepted backward-Euler solid response; no steady-state closure is implied.
#[derive(Debug, Clone)]
pub struct StepSolution {
    /// Published full nodal temperature at the end of the step, kelvin.
    pub temperature: Vec<f64>,
    /// Endpoint external Robin exchanges; suitable for existing air coupling.
    pub robin_fluxes: Vec<RobinFlux>,
    /// Endpoint internal contact exchanges, never counted as external power.
    pub contact_fluxes: Vec<InterfaceFlux>,
    /// Endpoint volumetric generation integrated with the production source rule.
    pub source_w: f64,
    /// Endpoint outward Neumann heat, watts.
    pub neumann_out_w: f64,
    /// Endpoint outward Robin heat, watts.
    pub robin_out_w: f64,
    /// Endpoint prescribed-temperature reaction, including consistent capacity.
    pub dirichlet_in_w: f64,
    /// `1^T C (T_new - T_old)` of the actually published temperatures, joules.
    pub stored_energy_change_j: f64,
    /// Storage minus dt times net external input, joules; checked independently.
    pub energy_residual_j: f64,
    /// Recomputed relative residual of the normalized correction solve (or
    /// the worst such inner residual over nonlinear Newton corrections).
    /// Final absolute-temperature rounding is checked by the separate energy gate.
    pub relative_residual: f64,
    /// Total Krylov iterations for this solid response.
    pub krylov_iterations: usize,
}

/// Capacity bound to one exact mesh. Reused across time steps and coupling trials.
/// History is supplied by immutable borrow on each call, never modified in place.
#[derive(Debug)]
pub struct BackwardEuler<'m> {
    mesh: &'m ConductionMesh,
    capacity: Csr,
}

impl<'m> BackwardEuler<'m> {
    /// Prepare one uniform, explicitly declared volumetric heat capacity.
    pub fn uniform(cx: &Cx<'_>, mesh: &'m ConductionMesh, capacity: VolumetricHeatCapacity)
        -> Result<Self, ConductionError>
    {
        poll(cx, 0)?;
        let capacity = assemble_capacitance(cx, mesh, capacity)?;
        poll(cx, 0)?;
        Ok(Self { mesh, capacity })
    }

    /// Prepare one declared volumetric heat capacity per tetrahedron, in mesh order.
    pub fn per_element(cx: &Cx<'_>, mesh: &'m ConductionMesh, capacities: &[VolumetricHeatCapacity])
        -> Result<Self, ConductionError>
    {
        poll(cx, 0)?;
        if capacities.len() != mesh.element_count() {
            return Err(ConductionError::FieldLength { field: "element heat capacity",
                expected: mesh.element_count(), found: capacities.len() });
        }
        let capacity = assemble_capacitance_from(cx, mesh, |e| capacities[e].value_j_per_m3_k())?;
        poll(cx, 0)?;
        Ok(Self { mesh, capacity })
    }

    /// Evaluate one constant-conductivity endpoint from immutable history.
    ///
    /// Existing Dirichlet values must already match history: an instantaneous
    /// prescribed-temperature jump is not silently assigned an energy impulse.
    /// Pure-Neumann transient problems are allowed because capacity anchors them.
    /// Invalid fields, nonlinear materials, nonconvergence and failed energy
    /// closure return no new state. Use [`Self::advance_nonlinear`] with an
    /// explicit policy for k(T). Constant matching-face contact is supported.
    pub fn advance(&self, cx: &Cx<'_>, problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>, old: &[f64], dt_s: f64, config: StepConfig)
        -> Result<StepSolution, ConductionError>
    {
        let dofs = self.admit_step(cx, problem, old, dt_s, config)?;
        for e in 0..self.mesh.element_count() {
            if e % 512 == 0 { poll(cx, e)?; }
            let model = match problem.element_materials {
                Some(materials) => materials.model_for(e)?, None => problem.material,
            };
            if model.is_temperature_dependent() {
                return Err(invalid("linear backward Euler requires temperature-independent conductivity; use advance_nonlinear with an explicit nonlinear policy"));
            }
        }
        let system = assemble_operator_scaled_with_interfaces(cx, self.mesh, problem.boundary,
            problem.material, problem.source, old, None, interfaces, problem.element_materials)?;
        let lhs = axpy_csr(&self.capacity, 1.0, &system.operator, dt_s);
        // The unknown is a CORRECTION. Its fixed entries are zero; the absolute
        // Dirichlet lift must not be added again to this right-hand side.
        let (matrix, _) = reduce_matrix_and_lift(&lhs, &dofs);
        let mut applied = vec![0.0; self.mesh.vertex_count()];
        system.operator.spmv(old, &mut applied);
        poll(cx, 0)?;
        let rhs: Vec<f64> = dofs.free().iter().map(|&v| finite(dt_s * (system.load[v] - applied[v])))
            .collect::<Result<_, _>>()?;
        let (correction, relative_residual, krylov_iterations) = solve(cx, &matrix, &rhs, config.linear)?;
        let mut temperature = old.to_vec();
        for (slot, &v) in dofs.free().iter().enumerate() {
            if slot % 512 == 0 { poll(cx, slot)?; }
            temperature[v] = finite(old[v] + correction[slot])?;
        }
        self.finish_step(cx, problem, interfaces, old, dt_s, config, &dofs,
            &system, temperature, relative_residual, krylov_iterations)
    }

    fn admit_step(&self, cx: &Cx<'_>, problem: ConductionProblem<'_>,
        old: &[f64], dt_s: f64, config: StepConfig) -> Result<DofMap, ConductionError>
    {
        poll(cx, 0)?;
        if !std::ptr::eq(self.mesh, problem.mesh) {
            return Err(invalid("capacity and spatial problem must use the same mesh"));
        }
        if !(dt_s.is_finite() && dt_s > 0.0 && config.energy_tolerance_j.is_finite()
            && config.energy_tolerance_j > 0.0 && config.linear.tolerance.is_finite()
            && config.linear.tolerance > 0.0 && config.linear.tolerance < 1.0
            && config.linear.max_iterations > 0)
        { return Err(invalid("positive finite step, energy tolerance and Krylov budget required")); }
        let n = self.mesh.vertex_count();
        if old.len() != n {
            return Err(ConductionError::FieldLength { field: "previous temperature", expected: n, found: old.len() });
        }
        for (i, &t) in old.iter().enumerate() { if i % 512 == 0 { poll(cx, i)?; } finite(t)?; }
        if let Some(materials) = problem.element_materials { materials.validate_for(self.mesh)?; }
        let dofs = DofMap::new(problem.boundary, n)?;
        for &v in dofs.fixed() {
            if old[v] != dofs.prescribed()[v] {
                return Err(invalid("history must match constant Dirichlet values; boundary jumps need an explicit impulse model"));
            }
        }
        Ok(dofs)
    }

    #[allow(clippy::too_many_arguments)]
    fn finish_step(&self, cx: &Cx<'_>, problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>, old: &[f64], dt_s: f64, config: StepConfig,
        dofs: &DofMap, system: &AssembledSystem, temperature: Vec<f64>,
        relative_residual: f64, krylov_iterations: usize) -> Result<StepSolution, ConductionError>
    {
        // Account using the rounded temperature field that the caller actually
        // receives, not an unpublished high-accuracy correction vector.
        let delta: Vec<f64> = temperature.iter().zip(old).map(|(a,b)| a-b).collect();
        let mut applied = vec![0.0; temperature.len()];
        self.capacity.spmv(&delta, &mut applied);
        poll(cx, 0)?;
        let stored_energy_change_j = sum(applied.iter().copied())?;
        let (energy, robin_fluxes) = energy_balance(self.mesh, problem.boundary, problem.source,
            system, dofs, &temperature);
        let dirichlet_in_w = finite(energy.dirichlet_in_w
            + sum(dofs.fixed().iter().map(|&v| applied[v]))? / dt_s)?;
        let net = finite(energy.source_w + dirichlet_in_w - energy.neumann_out_w - energy.robin_out_w)?;
        let energy_residual_j = finite(stored_energy_change_j - dt_s * net)?;
        if energy_residual_j.abs() > config.energy_tolerance_j {
            return Err(invalid(&format!("transient energy residual {energy_residual_j} J exceeds {} J", config.energy_tolerance_j)));
        }
        let contact_fluxes = match interfaces {
            Some(interfaces) => interfaces.fluxes(&temperature)?, None => Vec::new(),
        };
        poll(cx, krylov_iterations)?;
        Ok(StepSolution { temperature, robin_fluxes, contact_fluxes,
            source_w: energy.source_w, neumann_out_w: energy.neumann_out_w,
            robin_out_w: energy.robin_out_w, dirichlet_in_w, stored_energy_change_j,
            energy_residual_j, relative_residual, krylov_iterations })
    }
}

fn solve(cx: &Cx<'_>, matrix: &Csr, rhs: &[f64], config: LinearConfig)
    -> Result<(Vec<f64>, f64, usize), ConductionError>
{
    poll(cx, 0)?;
    let scale = rhs.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
    if scale == 0.0 { return Ok((vec![0.0; rhs.len()], 0.0, 0)); }
    let normalized: Vec<f64> = rhs.iter().map(|v| v/scale).collect();
    let op = CsrOp::symmetric(matrix.clone());
    let pre = spd_preconditioner(matrix);
    // A recursive CG stop can miss the recomputed gate by roundoff even
    // with work left (the 319.6 K transient-fan sizing regression). Repair
    // that actual defect, rather than accepting it or weakening the tolerance.
    let checked = checked_cg(&op, &pre, &normalized, CheckedCgConfig {
        tolerance: config.tolerance,
        max_iterations: config.max_iterations,
        max_corrections: 2,
    }, |iterations| poll(cx, iterations)).map_err(|error| match error {
        CheckedCgError::Interrupted(reason) => reason,
        CheckedCgError::InvalidInput(what) => invalid(what),
    })?;
    let relative = finite(checked.report.rel_residual)?;
    let iterations = checked.report.iters;
    if !checked.report.converged_euclidean() {
        return Err(ConductionError::LinearSolveFailed { iteration: 0, krylov_iterations: iterations,
            true_relative_residual: relative, tolerance: config.tolerance });
    }
    let solution = checked.x.iter().map(|v| finite(v*scale)).collect::<Result<_,_>>()?;
    poll(cx, iterations)?;
    Ok((solution, relative, iterations))
}
fn sum(values: impl IntoIterator<Item=f64>) -> Result<f64, ConductionError> {
    values.into_iter().try_fold(0.0, |s,v| finite(s+v))
}
fn finite(value: f64) -> Result<f64, ConductionError> {
    if value.is_finite() { Ok(value) } else {
        Err(ConductionError::NonFinite { field: "backward Euler", bits: value.to_bits() })
    }
}
fn invalid(what: &str) -> ConductionError {
    ConductionError::Config { parameter: "backward Euler", what: what.to_string() }
}
fn poll(cx: &Cx<'_>, at: usize) -> Result<(), ConductionError> {
    cx.checkpoint().map_err(|_| ConductionError::Cancelled { stage: "backward-euler-coupling", at })
}

#[cfg(test)]
mod tests;
