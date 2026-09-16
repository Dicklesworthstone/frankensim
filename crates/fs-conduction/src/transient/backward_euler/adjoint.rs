//! Discrete endpoint derivatives for immutable-history backward Euler.
//!
//! R = C(T-new - T-old)/dt + A(T-new)T-new - b. Its endpoint Jacobian is
//! C/dt + J-steady, NOT a steady-state tangent and NOT J evaluated at T-old.
//! A coupled air adjoint can use the same RobinResponse, then carry its total
//! nodal-load multiplier backwards as (C/dt)^T lambda. Time grids, geometry,
//! prescribed temperatures, material laws and contact resistance stay fixed.

use super::*;
use crate::adjoint::robin::RobinResponse;
use crate::assemble::assemble_jacobian_with_optional_interfaces;

/// A newly solved transient endpoint and its owned response operator.
/// Dereferencing exposes the same Robin tangent/adjoint operations as the
/// steady producer, without representing this state as a steady solution.
pub struct StepLinearization<'a> {
    primal: StepSolution,
    response: RobinResponse,
    capacity: &'a Csr,
    dofs: DofMap,
    old: Vec<f64>,
    dt_s: f64,
    source_load_w: Vec<f64>,
}

impl std::ops::Deref for StepLinearization<'_> {
    type Target = RobinResponse;
    fn deref(&self) -> &Self::Target { &self.response }
}

impl BackwardEuler<'_> {
    /// Solve one real endpoint and prepare its discrete response operator.
    /// Use the same old field, dt, references and policies as the forward run
    /// when reconstructing an accepted endpoint during a reverse traversal.
    ///
    /// No physical history is mutated. The primal obeys the existing residual,
    /// energy, material-range and cancellation gates. The derivative has the
    /// caller's linear budget. A material kink retains the primal but refuses
    /// tangent/adjoint requests, as on the steady Robin path.
    ///
    /// # Errors
    /// All errors from advance/advance_nonlinear and Robin port admission.
    #[allow(clippy::too_many_arguments)]
    pub fn linearize_step<'a>(
        &'a self, cx: &Cx<'_>, problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>, old: &[f64], dt_s: f64,
        config: StepConfig, nonlinear: Option<NonlinearStepConfig>, regions: &[&str],
    ) -> Result<StepLinearization<'a>, ConductionError> {
        let primal = match nonlinear {
            Some(policy) => self.advance_nonlinear(cx, problem, interfaces, old, dt_s, config, policy)?.step,
            None => self.advance(cx, problem, interfaces, old, dt_s, config)?,
        };
        let dofs = DofMap::new(problem.boundary, self.mesh.vertex_count())?;
        let jacobian = assemble_jacobian_with_optional_interfaces(cx, self.mesh,
            problem.boundary, problem.material, &primal.temperature, interfaces, problem.element_materials)?;
        let inverse_dt = finite(1.0 / dt_s)?;
        let full = axpy_csr(&self.capacity, inverse_dt, &jacobian, 1.0);
        let (matrix, _) = reduce_matrix_and_lift(&full, &dofs);
        let response = RobinResponse::for_endpoint(cx, problem, &primal, matrix,
            dofs.clone(), config.linear, regions)?;
        // Same P1 source integration and accumulation order as assemble.rs;
        // isolate the body load rather than subtracting large Robin loads.
        let mut source_load_w = vec![0.0; self.mesh.vertex_count()];
        for element in 0..self.mesh.element_count() {
            if element % 512 == 0 { poll(cx, element)?; }
            let tet = self.mesh.complex().tets[element];
            let volume = self.mesh.element_volume(element);
            for a in 0..4 {
                let mut load = 0.0;
                for (b, &vertex) in tet.iter().enumerate() {
                    let weight = if a == b { volume / 10.0 } else { volume / 20.0 };
                    load = weight.mul_add(problem.source.at(vertex as usize), load);
                }
                let vertex = tet[a] as usize;
                source_load_w[vertex] = finite(source_load_w[vertex] + load)?;
            }
        }
        poll(cx, 0)?;
        Ok(StepLinearization { primal, response, capacity: &self.capacity, dofs,
            old: old.to_vec(), dt_s, source_load_w })
    }
}

impl StepLinearization<'_> {
    /// The actual residual- and energy-checked transient primal.
    #[must_use]
    pub const fn primal(&self) -> &StepSolution { &self.primal }

    fn multiplier(&self, cx: &Cx<'_>, values: &[f64]) -> Result<Vec<f64>, ConductionError> {
        if values.len() != self.old.len() {
            return Err(ConductionError::FieldLength { field: "transient adjoint multiplier",
                expected: self.old.len(), found: values.len() });
        }
        let mut result = values.to_vec();
        for (i, &value) in values.iter().enumerate() {
            if i % 512 == 0 { poll(cx, i)?; }
            finite(value)?;
        }
        for &vertex in self.dofs.fixed() { result[vertex] = 0.0; }
        Ok(result)
    }

    /// Pull a TOTAL endpoint load multiplier back to the previous temperature.
    /// In a coupled problem pass the converged coupled nodal-load adjoint, not
    /// a frozen-reference solid adjoint. Prescribed-temperature controls are zero.
    pub fn previous_temperature_pullback(&self, cx: &Cx<'_>, nodal_load_adjoint: &[f64])
        -> Result<Vec<f64>, ConductionError> {
        let lambda = self.multiplier(cx, nodal_load_adjoint)?;
        let mut result = vec![0.0; lambda.len()];
        self.capacity.spmv(&lambda, &mut result);
        for (i, value) in result.iter_mut().enumerate() {
            if i % 512 == 0 { poll(cx, i)?; }
            *value = finite(*value / self.dt_s)?;
        }
        for &vertex in self.dofs.fixed() { result[vertex] = 0.0; }
        poll(cx, 0)?;
        Ok(result)
    }

    /// Derivative with respect to a multiplier on THIS endpoint's entire
    /// volumetric source at multiplier one. Zero source gives zero, not a
    /// division by its amplitude. This is not a Neumann or Robin control.
    pub fn source_multiplier_pullback(&self, cx: &Cx<'_>, nodal_load_adjoint: &[f64])
        -> Result<f64, ConductionError> {
        let lambda = self.multiplier(cx, nodal_load_adjoint)?;
        let mut value = 0.0;
        for (i, (&weight, &load)) in lambda.iter().zip(&self.source_load_w).enumerate() {
            if i % 512 == 0 { poll(cx, i)?; }
            value = finite(value + finite(weight * load)?)?;
        }
        Ok(value)
    }

    /// Derivative with respect to a common multiplier on the entire capacity
    /// matrix at multiplier one, holding the OLD state independently fixed.
    /// A trajectory sums this local contraction while propagating history.
    pub fn capacity_multiplier_pullback(&self, cx: &Cx<'_>, nodal_load_adjoint: &[f64])
        -> Result<f64, ConductionError> {
        let lambda = self.multiplier(cx, nodal_load_adjoint)?;
        let delta = self.primal.temperature.iter().zip(&self.old)
            .map(|(new, old)| finite(new - old)).collect::<Result<Vec<_>, _>>()?;
        let mut storage = vec![0.0; delta.len()];
        self.capacity.spmv(&delta, &mut storage);
        let mut value = 0.0;
        for (i, (&weight, &change)) in lambda.iter().zip(&storage).enumerate() {
            if i % 512 == 0 { poll(cx, i)?; }
            value = finite(value - finite(weight * finite(change / self.dt_s)?)?)?;
        }
        Ok(value)
    }
}
