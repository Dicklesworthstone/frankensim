//! Nonsymmetric material Jacobian for the existing Robin response API.
//! The residual, K'(T) assembly, preconditioner and Krylov algorithm all come
//! from the production conduction/solver stack; no second nonlinear model.

use fs_solver::{FgmresState, LinearOp};
use fs_sparse::{Csr, ops};

use super::{ConductionError, ConductionProblem, Cx, DofMap, LinearConfig,
    ThermalInterfaces, checked, failed, invalid, poll, true_residual};
use crate::assemble::{assemble_jacobian_with_optional_interfaces, element_temperature,
    reduce_matrix_and_lift};

pub(super) struct TangentSystem {
    transpose: Csr,
    pub(super) smooth: bool,
}

pub(super) fn prepare(
    cx: &Cx<'_>, problem: ConductionProblem<'_>, interfaces: Option<&ThermalInterfaces>,
    temperature: &[f64], dofs: &DofMap,
) -> Result<(Csr, TangentSystem), ConductionError> {
    let mut smooth = true;
    for element in 0..problem.mesh.element_count() {
        poll(cx, element)?;
        let model = match problem.element_materials {
            Some(assigned) => assigned.model_for(element)?,
            None => problem.material,
        };
        if !model.is_temperature_dependent() { continue; }
        let t = checked(element_temperature(problem.mesh, element, temperature))?;
        let slope = model.tensor_derivative_at(t)?;
        for row in slope { for value in row { checked(value)?; } }
        // The producer's piecewise-linear tables select the right segment at
        // a knot. Neighbouring representable temperatures expose a change of
        // slope or a validity endpoint without finite-differencing the PDE.
        // Equal slopes across a redundant knot remain differentiable.
        for neighbour in [t.next_down(), t.next_up()] {
            match model.tensor_derivative_at(neighbour) {
                Ok(other) => {
                    for row in other { for value in row { checked(value)?; } }
                    smooth &= slope == other;
                }
                Err(ConductionError::OutsideTemperatureSpan { .. }) => smooth = false,
                Err(error) => return Err(error),
            }
        }
    }
    let full = assemble_jacobian_with_optional_interfaces(cx, problem.mesh,
        problem.boundary, problem.material, temperature, interfaces, problem.element_materials)?;
    // Prescribed temperatures have zero perturbation: discard the PRIMAL lift.
    let (jacobian, _) = reduce_matrix_and_lift(&full, dofs);
    poll(cx, 0)?;
    let transpose = ops::transpose(&jacobian);
    poll(cx, 0)?;
    Ok((jacobian, TangentSystem { transpose, smooth }))
}

/// Borrow both orientations, avoiding a matrix/transpose clone per derivative.
struct Oriented<'a> { forward: &'a Csr, reverse: &'a Csr }
impl LinearOp for Oriented<'_> {
    fn n(&self) -> usize { self.forward.nrows() }
    fn apply(&self, x: &[f64], y: &mut [f64]) { self.forward.spmv(x, y); }
    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) { self.reverse.spmv(x, y); }
}

impl TangentSystem {
    pub(super) fn solve(
        &self, cx: &Cx<'_>, jacobian: &Csr, dofs: &DofMap, config: LinearConfig,
        full_rhs: &[f64], transposed: bool,
    ) -> Result<(Vec<f64>, f64, usize), ConductionError> {
        poll(cx, 0)?;
        if !self.smooth {
            return Err(invalid("material tangent lies at a slope discontinuity or sampled validity endpoint; no unique two-sided derivative is returned"));
        }
        if config.restart == 0 || config.max_iterations == 0 {
            return Err(invalid("positive nonlinear tangent restart and Krylov iteration budget required"));
        }
        let rhs = dofs.gather(full_rhs);
        let scale = rhs.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
        if scale == 0.0 { return Ok((vec![0.0; full_rhs.len()], 0.0, 0)); }
        let normalized: Vec<f64> = rhs.iter().map(|v| v / scale).collect();
        let operator = if transposed {
            Oriented { forward: &self.transpose, reverse: jacobian }
        } else { Oriented { forward: jacobian, reverse: &self.transpose } };
        let preconditioner = crate::solve::spd_preconditioner(operator.forward);
        // At most 64 Arnoldi columns between checkpoint boundaries. This is
        // a ceiling on the caller's restart, not a silent larger work budget.
        let restart = config.restart.min(dofs.n()).min(64);
        let mut state = FgmresState::new(&normalized, restart);
        while state.rel_residual() >= config.tolerance && state.iters < config.max_iterations {
            poll(cx, state.iters)?;
            let before = state.iters;
            // FGMRES retains only the iterate across cycles; shortening the
            // final cycle enforces the exact INNER-iteration budget.
            state.restart = restart.min(config.max_iterations - before);
            state.run(&operator, &preconditioner, &normalized, config.tolerance, 1);
            if state.iters == before { break; }
        }
        poll(cx, state.iters)?;
        let residual = true_residual(operator.forward, &state.x, &normalized)?;
        if residual >= config.tolerance { return Err(failed(state.iters, residual, config)); }
        let mut full = vec![0.0; full_rhs.len()];
        for (index, &vertex) in dofs.free().iter().enumerate() {
            poll(cx, index)?;
            full[vertex] = checked(state.x[index] * scale)?;
        }
        Ok((full, residual, state.iters))
    }
}

#[cfg(test)]
mod tests;
