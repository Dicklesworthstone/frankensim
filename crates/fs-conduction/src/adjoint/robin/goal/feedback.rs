//! Residual-checked solid seam for a consumer-owned Robin feedback law.

use super::super::{
    ConductionError, ConductionProblem, Cx, DofMap, LinearConfig, RobinResponse, ThermalInterfaces,
    add, admit_linear, assemble_operator_scaled_with_interfaces, bind_ports, checked, failed,
    invalid, nonlinear, poll, reduce, temperature_dependent, true_residual, vector,
};
use super::DiscreteGoalComparison;

/// Exact reference changes and transpose feedback supplied by the owner of
/// the external law r(T). This is numerical data, not a certified derivative.
/// `fs-airflow` constructs these from its actual exponential air-path model.
#[derive(Debug, Clone)]
pub struct RobinGoalFeedback {
    /// r(T_reference) minus the references bound in the solid boundary.
    pub reference_shift_k: Vec<f64>,
    /// r(T_reference) minus r(T_approximate), in selected port order.
    pub reference_difference_k: Vec<f64>,
    /// A^T B^T lambda where A maps wall means to external references and B
    /// maps references to solid nodal loads. One weight per wall mean.
    pub wall_adjoint: Vec<f64>,
}

/// Owned reference tangent and two actual solid residuals, prepared without
/// another primal solve. Both fields use the same mesh, material law, loads,
/// prescribed temperatures, and fixed contact resistance. A consumer may
/// close a Robin-reference feedback law through the response pullback and
/// then verify its full reduced transpose equation with `compare`.
pub struct RobinGoalLinearization {
    response: RobinResponse,
    reference_rhs: Vec<f64>,
    reference_residual: Vec<f64>,
    residual_difference: Vec<f64>,
    temperature_difference: Vec<f64>,
    tolerance: f64,
}

impl RobinGoalLinearization {
    /// Admit both fields and bind the actual material/contact tangent and
    /// selected uniform Robin ports. The supplied reference must first pass
    /// the solid residual gate at its declared boundary references. Inner
    /// pullbacks use one tenth of the requested tolerance, leaving room for
    /// the separately checked complete feedback equation.
    ///
    /// # Errors
    /// The same material, field, budget, prescribed-value, contact, reference
    /// residual and interruption refusals as `compare_discrete_goal`.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        linear: LinearConfig,
        reference_temperature: &[f64],
        approximate_temperature: &[f64],
        regions: &[&str],
    ) -> Result<Self, ConductionError> {
        poll(cx, 0)?;
        admit_linear(linear)?;
        let n = problem.mesh.vertex_count();
        vector(cx, reference_temperature, n)?;
        vector(cx, approximate_temperature, n)?;
        let dofs = DofMap::new(problem.boundary, n)?;
        if dofs.fixed().is_empty() && !problem.boundary.has_robin() {
            return Err(ConductionError::SingularPureNeumann);
        }
        for &vertex in dofs.fixed() {
            poll(cx, vertex)?;
            if reference_temperature[vertex] != dofs.prescribed()[vertex]
                || approximate_temperature[vertex] != dofs.prescribed()[vertex]
            {
                return Err(invalid(
                    "Robin goal comparison requires unchanged prescribed temperatures",
                ));
            }
        }
        let ports = bind_ports(cx, problem, regions)?;
        let assemble = |temperature| {
            assemble_operator_scaled_with_interfaces(
                cx,
                problem.mesh,
                problem.boundary,
                problem.material,
                problem.source,
                temperature,
                None,
                interfaces,
                problem.element_materials,
            )
        };
        let reference_system = assemble(reference_temperature)?;
        let (matrix, reference_rhs) = reduce(&reference_system, &dofs);
        let relative = true_residual(&matrix, &dofs.gather(reference_temperature), &reference_rhs)?;
        if relative >= linear.tolerance {
            return Err(failed(0, relative, linear));
        }
        let reference_residual =
            crate::assemble::residual(&reference_system, &dofs, reference_temperature);
        let approximate_system = assemble(approximate_temperature)?;
        let approximate_residual =
            crate::assemble::residual(&approximate_system, &dofs, approximate_temperature);
        let residual_difference = reference_residual
            .iter()
            .zip(&approximate_residual)
            .map(|(a, b)| checked(a - b))
            .collect::<Result<Vec<_>, _>>()?;
        let temperature_difference = reference_temperature
            .iter()
            .zip(approximate_temperature)
            .map(|(a, b)| checked(a - b))
            .collect::<Result<Vec<_>, _>>()?;
        let (_, robin_fluxes) = crate::solve::energy_balance(
            problem.mesh,
            problem.boundary,
            problem.source,
            &reference_system,
            &dofs,
            reference_temperature,
        );
        let (matrix, nonlinear) = if temperature_dependent(cx, problem)? {
            if linear.restart == 0 {
                return Err(invalid("nonlinear Robin goal needs positive restart"));
            }
            let (matrix, tangent) =
                nonlinear::prepare(cx, problem, interfaces, reference_temperature, &dofs)?;
            (matrix, Some(tangent))
        } else {
            (matrix, None)
        };
        let inner = LinearConfig {
            tolerance: linear.tolerance * 0.1,
            ..linear
        };
        admit_linear(inner)?;
        let response = RobinResponse {
            temperature: reference_temperature.to_vec(),
            robin_fluxes,
            matrix,
            dofs,
            ports,
            linear: inner,
            nonlinear,
        };
        poll(cx, 0)?;
        Ok(Self {
            response,
            reference_rhs,
            reference_residual,
            residual_difference,
            temperature_difference,
            tolerance: linear.tolerance,
        })
    }

    /// Checked reference response used by the external implicit solve.
    #[must_use]
    pub const fn response(&self) -> &RobinResponse {
        &self.response
    }

    /// Verify the complete reduced primal and transpose equations, then
    /// accumulate lambda_i times the corrected residual difference. With
    /// R_c(T)=R_s(T)-B[r(T)-r_bound], the transpose is J_s^T-W^T A^T B^T.
    /// The residual is normalized against the ORIGINAL goal, so a large
    /// feedback term cannot make an inaccurate coupled adjoint appear good.
    /// Fixed-node adjoints must be zero. A consumer must supply the actual
    /// feedback law and transpose; this seam does not establish their truth.
    ///
    /// # Errors
    /// Refuses malformed/nonfinite feedback, a nonsmooth tangent, failed
    /// coupled primal or dual gate, changed fixed-node adjoints, cancellation,
    /// or unrepresentable arithmetic. No partial comparison is published.
    pub fn compare(
        &self,
        cx: &Cx<'_>,
        nodal_weights: &[f64],
        adjoint: &[f64],
        feedback: &RobinGoalFeedback,
        dual_iterations: usize,
    ) -> Result<DiscreteGoalComparison, ConductionError> {
        let response = &self.response;
        let n = response.temperature.len();
        vector(cx, nodal_weights, n)?;
        vector(cx, adjoint, n)?;
        if !response.has_smooth_material_tangent() {
            return Err(invalid("Robin goal has no unique smooth material tangent"));
        }
        for values in [
            &feedback.reference_shift_k,
            &feedback.reference_difference_k,
            &feedback.wall_adjoint,
        ] {
            vector(cx, values, response.ports.len())?;
        }
        for &vertex in response.dofs.fixed() {
            if adjoint[vertex] != 0.0 {
                return Err(invalid("prescribed-node adjoint must be zero"));
            }
        }
        let mut reference_load = vec![0.0; n];
        let mut difference_load = vec![0.0; n];
        let mut wall_adjoint = vec![0.0; n];
        for (i, port) in response.ports.iter().enumerate() {
            for (vertices, area) in &port.faces {
                poll(cx, i)?;
                for &vertex in vertices {
                    let b = port.htc_w_m2_k * (area / 3.0);
                    add(
                        &mut reference_load[vertex],
                        b * feedback.reference_shift_k[i],
                    )?;
                    add(
                        &mut difference_load[vertex],
                        b * feedback.reference_difference_k[i],
                    )?;
                    add(
                        &mut wall_adjoint[vertex],
                        (area / port.area_m2 / 3.0) * feedback.wall_adjoint[i],
                    )?;
                }
            }
        }
        let mut coupled_residual = Vec::with_capacity(response.dofs.n());
        let mut coupled_rhs = Vec::with_capacity(response.dofs.n());
        for (i, &vertex) in response.dofs.free().iter().enumerate() {
            poll(cx, i)?;
            coupled_residual.push(checked(
                self.reference_residual[i] - reference_load[vertex],
            )?);
            coupled_rhs.push(checked(self.reference_rhs[i] + reference_load[vertex])?);
        }
        let primal_relative_residual = residual_ratio(&coupled_residual, &coupled_rhs)?;
        let config = LinearConfig {
            tolerance: self.tolerance,
            ..response.linear
        };
        if primal_relative_residual >= self.tolerance {
            return Err(failed(0, primal_relative_residual, config));
        }
        let weights = response.dofs.gather(nodal_weights);
        let scale = weights
            .iter()
            .map(|v| v.abs())
            .fold(0.0_f64, f64::max)
            .max(f64::MIN_POSITIVE);
        let mut transpose = vec![0.0; response.dofs.n()];
        for (i, &vertex) in response.dofs.free().iter().enumerate() {
            poll(cx, i)?;
            let (columns, values) = response.matrix.row(i);
            for (&column, &value) in columns.iter().zip(values) {
                add(&mut transpose[column], value * (adjoint[vertex] / scale))?;
            }
        }
        let mut dual_residual = Vec::with_capacity(weights.len());
        let mut normalized_weights = Vec::with_capacity(weights.len());
        for (i, &vertex) in response.dofs.free().iter().enumerate() {
            poll(cx, i)?;
            let weight = checked(weights[i] / scale)?;
            normalized_weights.push(weight);
            dual_residual.push(checked(
                transpose[i] - wall_adjoint[vertex] / scale - weight,
            )?);
        }
        let dual_relative_residual = residual_ratio(&dual_residual, &normalized_weights)?;
        if dual_relative_residual >= self.tolerance {
            return Err(failed(dual_iterations, dual_relative_residual, config));
        }
        let mut nodal_contributions = vec![0.0; n];
        let mut signed_residual_change = 0.0;
        let mut signed_goal_change = 0.0;
        for (i, &vertex) in response.dofs.free().iter().enumerate() {
            poll(cx, i)?;
            let difference = checked(self.residual_difference[i] - difference_load[vertex])?;
            nodal_contributions[vertex] = checked(adjoint[vertex] * difference)?;
            add(&mut signed_residual_change, nodal_contributions[vertex])?;
        }
        for (i, &weight) in nodal_weights.iter().enumerate() {
            poll(cx, i)?;
            add(
                &mut signed_goal_change,
                weight * self.temperature_difference[i],
            )?;
        }
        let linearization_remainder = checked(signed_goal_change - signed_residual_change)?;
        poll(cx, dual_iterations)?;
        Ok(DiscreteGoalComparison {
            signed_residual_change,
            linearization_remainder,
            signed_goal_change,
            nodal_contributions,
            primal_relative_residual,
            dual_relative_residual,
            dual_iterations,
            uses_nonlinear_jacobian: response.uses_nonlinear_jacobian(),
        })
    }
}

fn residual_ratio(residual: &[f64], rhs: &[f64]) -> Result<f64, ConductionError> {
    let scale = residual
        .iter()
        .chain(rhs)
        .map(|v| v.abs())
        .fold(0.0_f64, f64::max);
    if scale == 0.0 {
        return Ok(0.0);
    }
    let r = residual
        .iter()
        .map(|v| checked(v / scale))
        .collect::<Result<Vec<_>, _>>()?;
    let b = rhs
        .iter()
        .map(|v| checked(v / scale))
        .collect::<Result<Vec<_>, _>>()?;
    checked(fs_solver::norm2(&r) / fs_solver::norm2(&b).max(f64::MIN_POSITIVE))
}
