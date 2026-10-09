//! Two-field goal estimates for actual area-mean-dependent Robin laws.
//! Both residuals are reassembled at their own law values. The reference
//! adjoint includes coefficient AND reference feedback and the observed
//! finite linearization remainder is retained, never called a bound.

use super::super::{
    ConductionError, ConductionProblem, Cx, DofMap, LinearConfig, RobinResponse,
    ThermalInterfaces, add, admit_linear, assemble_operator_scaled_with_interfaces,
    checked, failed, invalid, poll, reduce, temperature_dependent, true_residual, vector,
};
use super::super::goal::DiscreteGoalComparison;

mod radiation;

impl RobinResponse {
    /// Compare two full fields under the same area-mean-dependent Robin laws.
    ///
    /// `reference_points[j]` is `[h, r, dh/dmean, dr/dmean]` evaluated at the
    /// REFERENCE field's area-mean temperature on `regions[j]`. The units are
    /// W/(m^2 K), K, W/(m^2 K^2), and K/K. `approximate_values[j]` is `[h, r]`
    /// evaluated at the APPROXIMATE field's mean, not at the reference mean.
    /// The law owner is responsible for these evaluations, domains and slopes;
    /// supplying data here does not certify their physical or derivative truth.
    ///
    /// Only the selected uniform Robin rows are rebound. Geometry, source,
    /// prescribed temperatures, material assignments and matching contacts stay
    /// unchanged. The complete reference residual is checked even for a zero
    /// goal. One genuine transposed tangent solve includes smooth k(T), contact,
    /// dh/dmean and dr/dmean, retaining the consistent face mass matrix. No
    /// primal solve and no differentiation of solver iterations is performed.
    ///
    /// The signed nodal contributions are lambda_i [R(T_ref)-R(T_approx)]_i.
    /// Their sum plus the returned finite linearization remainder equals the
    /// directly observed linear goal difference, to rounding. Approximate law
    /// values need no derivative and may cross smooth segments. This is an
    /// Estimated discrete comparison, not a continuum discretization bound,
    /// effectivity certificate, or a bound on omitted physics. The existing
    /// `uses_nonlinear_jacobian` field continues to describe the material
    /// tangent; Robin feedback is always included by this method.
    ///
    /// The two factors per selected region must fit `max_feedback_entries`;
    /// at most 64 regions are admitted. This bounds feedback storage, not the
    /// entire assembly or allocator. Existing Krylov and cancellation gates
    /// apply. Any failure leaves the input fields and boundaries unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn compare_mean_robin_goal_at(
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        linear: LinearConfig,
        reference_temperature: &[f64],
        approximate_temperature: &[f64],
        regions: &[&str],
        reference_points: &[[f64; 4]],
        approximate_values: &[[f64; 2]],
        nodal_weights: &[f64],
        max_feedback_entries: usize,
    ) -> Result<DiscreteGoalComparison, ConductionError> {
        Self::compare_mean_robin_goal_inner(
            cx, problem, interfaces, linear, reference_temperature, approximate_temperature,
            regions, reference_points, approximate_values, None, nodal_weights, max_feedback_entries,
        )
    }

    /// Include reference feedback between different mean-temperature regions.
    /// `reference_feedback[i*m+j]` is the additional derivative of reference i
    /// with respect to mean j, in K/K. Its diagonal adds to the local slopes;
    /// never supply a derivative twice. Both fields retain their OWN actual
    /// law values. This is the same estimated two-field comparison as
    /// `compare_mean_robin_goal_at`, with the complete coupled transpose and
    /// four retained factor vectors per region instead of two. Every physical
    /// residual, material-domain, prescribed-value and Krylov gate still applies.
    #[allow(clippy::too_many_arguments)]
    pub fn compare_mean_robin_goal_with_reference_feedback_at(
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        linear: LinearConfig,
        reference_temperature: &[f64],
        approximate_temperature: &[f64],
        regions: &[&str],
        reference_points: &[[f64; 4]],
        approximate_values: &[[f64; 2]],
        reference_feedback: &[f64],
        nodal_weights: &[f64],
        max_feedback_entries: usize,
    ) -> Result<DiscreteGoalComparison, ConductionError> {
        Self::compare_mean_robin_goal_inner(
            cx, problem, interfaces, linear, reference_temperature, approximate_temperature,
            regions, reference_points, approximate_values, Some(reference_feedback), nodal_weights, max_feedback_entries,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn compare_mean_robin_goal_inner(
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        linear: LinearConfig,
        reference_temperature: &[f64],
        approximate_temperature: &[f64],
        regions: &[&str],
        reference_points: &[[f64; 4]],
        approximate_values: &[[f64; 2]],
        reference_feedback: Option<&[f64]>,
        nodal_weights: &[f64],
        max_feedback_entries: usize,
    ) -> Result<DiscreteGoalComparison, ConductionError> {
        poll(cx, 0)?;
        admit_linear(linear)?;
        let n = problem.mesh.vertex_count();
        let factors = if reference_feedback.is_some() { 4 } else { 2 };
        if regions.len() > 64 || regions.len() != reference_points.len()
            || regions.len() != approximate_values.len()
            || n.checked_mul(factors).and_then(|v| v.checked_mul(regions.len()))
                .is_none_or(|v| v > max_feedback_entries)
        {
            return Err(invalid("mean-Robin goal needs matching law rows within the 64-region and factor-entry budgets"));
        }
        if let Some(feedback) = reference_feedback {
            vector(cx, feedback, regions.len()*regions.len())?;
        }
        for field in [reference_temperature, approximate_temperature, nodal_weights] {
            vector(cx, field, n)?;
        }
        let mut reference_rows = Vec::with_capacity(regions.len());
        let mut approximate_rows = Vec::with_capacity(regions.len());
        let mut slopes = Vec::with_capacity(regions.len());
        let mut reference_slopes = Vec::with_capacity(regions.len());
        for (j, &name) in regions.iter().enumerate() {
            poll(cx, j)?;
            vector(cx, &reference_points[j], 4)?;
            vector(cx, &approximate_values[j], 2)?;
            let index = problem.boundary.region_names().iter().position(|v| v == name)
                .ok_or_else(|| invalid("unknown mean-Robin goal region"))?;
            let [h, r, dh, dr] = reference_points[j];
            let [ha, ra] = approximate_values[j];
            reference_rows.push((index, h, r));
            approximate_rows.push((index, ha, ra));
            slopes.push(dh);
            reference_slopes.push(dr);
        }
        let reference_boundary = problem.boundary.with_uniform_robin_replacements(&reference_rows)?;
        let approximate_boundary = problem.boundary.with_uniform_robin_replacements(&approximate_rows)?;
        let reference_problem = ConductionProblem { boundary: &reference_boundary, ..problem };
        let dofs = DofMap::new(&reference_boundary, n)?;
        for &v in dofs.fixed() {
            poll(cx, v)?;
            if reference_temperature[v] != dofs.prescribed()[v]
                || approximate_temperature[v] != dofs.prescribed()[v]
            {
                return Err(invalid("mean-Robin goal comparison requires unchanged prescribed temperatures in both fields"));
            }
        }
        let assemble = |boundary, temperature| {
            assemble_operator_scaled_with_interfaces(cx, problem.mesh, boundary,
                problem.material, problem.source, temperature, None, interfaces,
                problem.element_materials)
        };
        let reference_system = assemble(&reference_boundary, reference_temperature)?;
        let (matrix, rhs) = reduce(&reference_system, &dofs);
        let primal_relative_residual =
            true_residual(&matrix, &dofs.gather(reference_temperature), &rhs)?;
        if primal_relative_residual >= linear.tolerance {
            return Err(failed(0, primal_relative_residual, linear));
        }
        let reference_residual =
            crate::assemble::residual(&reference_system, &dofs, reference_temperature);
        let approximate_system = assemble(&approximate_boundary, approximate_temperature)?;
        let approximate_residual =
            crate::assemble::residual(&approximate_system, &dofs, approximate_temperature);
        // Use the already verified complete, nonsymmetric adjoint path. The
        // fixed-state partial and mean feedback are not counted a second time.
        let gradient = if let Some(feedback) = reference_feedback {
            Self::pullback_mean_robin_with_reference_feedback(cx, reference_problem,
                interfaces, linear, reference_temperature, regions, &slopes,
                &reference_slopes, feedback, nodal_weights, max_feedback_entries)?
        } else {
            Self::pullback_mean_robin(cx, reference_problem, interfaces, linear,
                reference_temperature, regions, &slopes, &reference_slopes, nodal_weights,
                max_feedback_entries)?
        };
        let mut nodal_contributions = vec![0.0; n];
        let mut signed_residual_change = 0.0;
        for (i, &v) in dofs.free().iter().enumerate() {
            poll(cx, i)?;
            let difference = checked(reference_residual[i] - approximate_residual[i])?;
            let contribution = checked(gradient.nodal_load[v] * difference)?;
            nodal_contributions[v] = contribution;
            add(&mut signed_residual_change, contribution)?;
        }
        let mut signed_goal_change = 0.0;
        for v in 0..n {
            poll(cx, v)?;
            let difference = checked(reference_temperature[v] - approximate_temperature[v])?;
            add(&mut signed_goal_change, checked(nodal_weights[v] * difference)?)?;
        }
        let linearization_remainder = checked(signed_goal_change - signed_residual_change)?;
        let uses_nonlinear_jacobian = temperature_dependent(cx, reference_problem)?;
        poll(cx, gradient.iterations)?;
        Ok(DiscreteGoalComparison {
            signed_residual_change, linearization_remainder, signed_goal_change,
            nodal_contributions, primal_relative_residual,
            dual_relative_residual: gradient.relative_residual,
            dual_iterations: gradient.iterations, uses_nonlinear_jacobian,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{ConductionMesh, ConductivityModel, ConductivityTable, ScalarField,
        ThermalBoundary, ThermalBoundaryBuilder, ThermalBc, SolveConfig, InitialGuess};
    use crate::fixtures::{box_grid, on_box_face};
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

    fn with_gate<T>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> T) -> T {
        ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(gate, arena,
            StreamKey { seed: 51, kernel_id: 820, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic)))
    }
    fn config() -> SolveConfig {
        let mut c = SolveConfig::default();
        c.initial = InitialGuess::Uniform(320.0);
        c.linear.tolerance = 1e-10;
        c.stop.residual_rtol = 1e-12;
        c.stop.step_atol = 0.0;
        c
    }
    fn wall_mean(mesh: &ConductionMesh, temperature: &[f64]) -> f64 {
        let mut area = 0.0;
        let mut integral = 0.0;
        for face in mesh.boundary().iter().filter(|f| on_box_face(f.centroid[0], 0.1)) {
            area += face.area;
            integral += face.area / 3.0 * face.vertices.iter()
                .map(|&v| temperature[v as usize]).sum::<f64>();
        }
        integral / area
    }
    fn point(mean: f64) -> [f64; 4] {
        [40.0 + 0.8 * (mean - 313.0), 293.0 + 0.35 * (mean - 293.0), 0.8, 0.35]
    }
    fn boundary(mesh: &ConductionMesh, h: f64, r: f64) -> ThermalBoundary {
        ThermalBoundaryBuilder::new(mesh)
            .region("hot", |f| on_box_face(f.centroid[0], 0.0),
                ThermalBc::dirichlet(340.0).unwrap()).unwrap()
            .region("wall", |f| on_box_face(f.centroid[0], 0.1),
                ThermalBc::robin(h, r).unwrap()).unwrap()
            .adiabatic_remainder().finish().unwrap()
    }
    fn physical(cx: &Cx<'_>) -> (ConductionMesh, ConductivityModel, ScalarField, ThermalBoundary, Vec<f64>) {
        let (complex, positions) = box_grid([3, 2, 2], [0.1, 0.04, 0.03]);
        let mesh = ConductionMesh::new(complex, positions).unwrap();
        let material = ConductivityModel::isotropic(ConductivityTable::declared_curve(
            vec![(250.0, 3.0), (450.0, 13.0)]).unwrap());
        let source = ScalarField::Nodal(mesh.positions().iter()
            .map(|p| 1000.0 + 50000.0 * p[1]).collect());
        let mut h = 40.0;
        let mut r = 293.0;
        for _ in 0..160 {
            let bc = boundary(&mesh, h, r);
            let problem = ConductionProblem { mesh: &mesh, material: &material,
                source: &source, boundary: &bc, element_materials: None };
            let solution = crate::solve(cx, problem, config()).unwrap();
            let [next_h, next_r, _, _] = point(wall_mean(&mesh, &solution.temperature));
            if (h-next_h).abs() < 1e-13*h && (r-next_r).abs() < 1e-13*r {
                return (mesh, material, source, bc, solution.temperature);
            }
            h = 0.5*(h+next_h);
            r = 0.5*(r+next_r);
        }
        panic!("physical mean-Robin fixed point did not converge");
    }

    #[test]
    fn total_goal_remainder_is_second_order_and_both_law_states_are_used() {
        with_gate(&CancelGate::new_clock_free(), |cx| {
            let (mesh, material, source, bc, reference) = physical(cx);
            let problem = ConductionProblem { mesh: &mesh, material: &material,
                source: &source, boundary: &bc, element_materials: None };
            let n = mesh.vertex_count();
            let mut weights = vec![0.0; n]; weights[n-1] = 1.0;
            let p = point(wall_mean(&mesh, &reference));
            let mut remainders = Vec::new();
            for step in [1.0, 0.5, 0.25] {
                let approximate: Vec<_> = mesh.positions().iter().zip(&reference)
                    .map(|(x, t)| t - step * x[0] / 0.1 * (1.0 + x[1] / 0.04)).collect();
                let q = point(wall_mean(&mesh, &approximate));
                let report = RobinResponse::compare_mean_robin_goal_at(cx, problem, None,
                    config().linear, &reference, &approximate, &["wall"], &[p],
                    &[[q[0], q[1]]], &weights, 2*n).unwrap();
                let expected = reference[n-1] - approximate[n-1];
                assert!((report.signed_goal_change - expected).abs() < 1e-12);
                assert!((report.signed_residual_change + report.linearization_remainder
                    - expected).abs() < 1e-12);
                assert!((report.nodal_contributions.iter().sum::<f64>()
                    - report.signed_residual_change).abs() < 1e-12);
                assert!(report.primal_relative_residual < config().linear.tolerance);
                assert!(report.dual_relative_residual < config().linear.tolerance);
                assert!(report.uses_nonlinear_jacobian);
                for &v in DofMap::new(&bc, n).unwrap().fixed() {
                    assert_eq!(report.nodal_contributions[v], 0.0);
                }
                let wrong = RobinResponse::compare_mean_robin_goal_at(cx, problem, None,
                    config().linear, &reference, &approximate, &["wall"], &[p],
                    &[[p[0], p[1]]], &weights, 2*n).unwrap();
                assert!((wrong.signed_residual_change-report.signed_residual_change).abs()
                    > 1e-3 * step, "fixture must distinguish a frozen approximate boundary");
                remainders.push(report.linearization_remainder.abs());
            }
            assert!(remainders[0] > 1e-6);
            for pair in remainders.windows(2) {
                assert!(pair[0] / pair[1] > 3.7 && pair[0] / pair[1] < 4.3,
                    "a total tangent leaves a quadratic remainder: {remainders:?}");
            }
        });
    }

    #[test]
    fn identical_zero_and_scaled_goals_preserve_the_checked_equations() {
        with_gate(&CancelGate::new_clock_free(), |cx| {
            let (mesh, material, source, bc, reference) = physical(cx);
            let problem = ConductionProblem { mesh: &mesh, material: &material,
                source: &source, boundary: &bc, element_materials: None };
            let n = mesh.vertex_count();
            let p = point(wall_mean(&mesh, &reference));
            let q = [[p[0], p[1]]];
            let mut weights = vec![0.0; n]; weights[n-1] = 1.0;
            let identical = RobinResponse::compare_mean_robin_goal_at(cx, problem, None,
                config().linear, &reference, &reference, &["wall"], &[p], &q, &weights, 2*n).unwrap();
            assert_eq!(identical.signed_goal_change, 0.0);
            assert_eq!(identical.signed_residual_change, 0.0);
            assert_eq!(identical.linearization_remainder, 0.0);
            let approximate: Vec<_> = mesh.positions().iter().zip(&reference)
                .map(|(x, t)| t - x[0] / 0.1).collect();
            let a = point(wall_mean(&mesh, &approximate));
            let q = [[a[0], a[1]]];
            let base = RobinResponse::compare_mean_robin_goal_at(cx, problem, None,
                config().linear, &reference, &approximate, &["wall"], &[p], &q, &weights, 2*n).unwrap();
            for scale in [-3.0, 1e-9, 1e9] {
                weights[n-1] = scale;
                let result = RobinResponse::compare_mean_robin_goal_at(cx, problem, None,
                    config().linear, &reference, &approximate, &["wall"], &[p], &q, &weights, 2*n).unwrap();
                assert!((result.signed_residual_change / scale - base.signed_residual_change).abs() < 1e-8);
                assert!((result.linearization_remainder / scale - base.linearization_remainder).abs() < 1e-8);
            }
            let zero = RobinResponse::compare_mean_robin_goal_at(cx, problem, None,
                config().linear, &reference, &approximate, &["wall"], &[p], &q, &vec![0.0;n], 2*n).unwrap();
            assert_eq!(zero.dual_iterations, 0);
            assert_eq!(zero.signed_residual_change, 0.0);
            assert_eq!(zero.linearization_remainder, 0.0);
        });
    }

    #[test]
    fn malformed_unconverged_over_budget_and_cancelled_inputs_refuse() {
        with_gate(&CancelGate::new_clock_free(), |cx| {
            let (mesh, material, source, bc, reference) = physical(cx);
            let problem = ConductionProblem { mesh: &mesh, material: &material,
                source: &source, boundary: &bc, element_materials: None };
            let n = mesh.vertex_count();
            let p = point(wall_mean(&mesh, &reference));
            let q = [[p[0], p[1]]];
            let mut weights = vec![0.0; n]; weights[n-1] = 1.0;
            let compare = |cx: &Cx<'_>, linear, r: &[f64], a: &[f64],
                names: &[&str], points: &[[f64;4]], values: &[[f64;2]], w: &[f64], budget| {
                RobinResponse::compare_mean_robin_goal_at(cx, problem, None, linear,
                    r, a, names, points, values, w, budget)
            };
            assert!(compare(cx, config().linear, &reference, &reference, &["wall"], &[p], &q, &weights, 2*n-1).is_err());
            assert!(compare(cx, config().linear, &reference, &reference, &["wall"], &[], &q, &weights, 2*n).is_err());
            assert!(compare(cx, config().linear, &reference, &reference, &["absent"], &[p], &q, &weights, 2*n).is_err());
            assert!(compare(cx, config().linear, &reference, &reference, &["wall", "wall"], &[p,p], &[q[0],q[0]], &weights, 4*n).is_err());
            for values in [[0.0, p[1]], [-1.0, p[1]], [f64::NAN, p[1]], [p[0], f64::INFINITY]] {
                assert!(compare(cx, config().linear, &reference, &reference, &["wall"], &[p], &[values], &weights, 2*n).is_err());
            }
            let mut changed = reference.clone(); changed[0] += 1.0;
            assert!(compare(cx, config().linear, &reference, &changed, &["wall"], &[p], &q, &weights, 2*n).is_err());
            let mut stale = p; stale[1] += 10.0;
            assert!(compare(cx, config().linear, &reference, &reference, &["wall"], &[stale], &q, &vec![0.0;n], 2*n).is_err(),
                "even a zero goal must reject an unconverged physical reference");
            let short = LinearConfig { max_iterations: 1, restart: 1, ..config().linear };
            assert!(matches!(compare(cx, short, &reference, &reference, &["wall"], &[p], &q, &weights, 2*n),
                Err(ConductionError::LinearSolveFailed { .. })));
            let gate = CancelGate::new_clock_free(); gate.request();
            with_gate(&gate, |cancelled| assert!(matches!(compare(cancelled, config().linear,
                &reference, &reference, &["wall"], &[p], &q, &weights, 2*n),
                Err(ConductionError::Cancelled { .. }))));
        });
    }
}
