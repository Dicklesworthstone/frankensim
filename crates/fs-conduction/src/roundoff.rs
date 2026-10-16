//! A-priori roundoff bound for one nodal temperature of a steady solve.
//!
//! The computed field `T` solves a floating-point version of `A(T) T = b`:
//! assembly rounds every entry of `A` and `b`, and the residual gate is
//! itself evaluated in floating point. Each assembled row is a sum of at
//! most `k` rounded products, so the discrete equations the solver actually
//! satisfied differ from the exact ones by a row perturbation bounded
//! componentwise by `γ_k (|b| + |A||T|)` (Higham, *Accuracy and Stability
//! of Numerical Algorithms*, §3.1, `γ_k = k u / (1 − k u)`).
//!
//! That perturbation reaches the nodal QoI `T_v` through the adjoint
//! `λ = A⁻ᵀ e_v`, so `|δT_v| ≤ γ_k Σ_i |λ_i| (|b_i| + (|A||T|)_i)`.
//!
//! `k` counts the row's non-zeros plus a per-element operation allowance for
//! every tetrahedron that contributes to the row. The bound uses the
//! ASSEMBLED magnitudes. Element contributions that cancel inside one entry
//! are not tracked individually, so this is an Estimated componentwise
//! bound, not an interval certificate. The solver's own stopping error is a
//! separate term (solver-algebraic) and is not included here.

use fs_exec::Cx;
use fs_solver::{CgState, CsrOp, norm2};

use crate::assemble::{DofMap, assemble_operator_scaled_with_interfaces, reduce_matrix_and_lift};
use crate::{ConductionError, ConductionProblem, LinearConfig, ThermalInterfaces};

/// Floating-point operations allowed per tetrahedron contribution to one
/// assembled entry (gradients, the 3×3 tensor product, the volume scale and
/// the Robin/source quadrature), generously rounded up.
pub const OPERATIONS_PER_ELEMENT_CONTRIBUTION: usize = 64;

/// The bound and the evidence that produced it.
#[derive(Debug, Clone, PartialEq)]
pub struct NodalRoundoffBound {
    /// Bound on the nodal temperature's roundoff error, K. Zero on a
    /// prescribed (Dirichlet) vertex, whose value is declared data.
    pub half_width_k: f64,
    /// `γ_k` used for every row.
    pub gamma: f64,
    /// The largest per-row operation count `k`.
    pub operations: usize,
    /// Recomputed relative residual of the adjoint solve.
    pub adjoint_relative_residual: f64,
    /// Conjugate-gradient iterations of the adjoint solve.
    pub adjoint_iterations: usize,
}

/// Bound the roundoff error of `temperature[vertex]` for the exact problem
/// (boundary, materials, contacts) that produced `temperature`.
///
/// # Errors
/// A mis-sized field or out-of-range vertex, assembly refusals (including a
/// temperature outside a material's sampled span), an adjoint solve that
/// misses `linear.tolerance`, non-finite arithmetic, or cancellation.
pub fn nodal_roundoff_bound(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    interfaces: Option<&ThermalInterfaces>,
    temperature: &[f64],
    vertex: usize,
    linear: LinearConfig,
) -> Result<NodalRoundoffBound, ConductionError> {
    let n = problem.mesh.vertex_count();
    if temperature.len() != n {
        return Err(ConductionError::FieldLength {
            field: "roundoff temperature",
            expected: n,
            found: temperature.len(),
        });
    }
    if vertex >= n {
        return Err(ConductionError::FieldLength {
            field: "roundoff QoI vertex",
            expected: n,
            found: vertex,
        });
    }
    let system = assemble_operator_scaled_with_interfaces(
        cx,
        problem.mesh,
        problem.boundary,
        problem.material,
        problem.source,
        temperature,
        None,
        interfaces,
        problem.element_materials,
    )?;
    let dofs = DofMap::new(problem.boundary, n)?;
    let mut elements_at = vec![0usize; n];
    for tet in &problem.mesh.complex().tets {
        for &v in tet {
            elements_at[v as usize] += 1;
        }
    }
    let unit = f64::EPSILON / 2.0;
    let mut operations = 0usize;
    for &v in dofs.free() {
        let row = system.operator.row(v).0.len()
            + OPERATIONS_PER_ELEMENT_CONTRIBUTION * elements_at[v].max(1);
        operations = operations.max(row);
    }
    #[allow(clippy::cast_precision_loss)] // operation counts are far below 2^53
    let ku = operations as f64 * unit;
    if !(ku < 0.5) {
        return Err(ConductionError::Config {
            parameter: "roundoff operations",
            what: format!("k*u = {ku} leaves gamma_k undefined"),
        });
    }
    let gamma = ku / (1.0 - ku);
    let Some(slot) = dofs.slot_of(vertex) else {
        return Ok(NodalRoundoffBound {
            half_width_k: 0.0,
            gamma,
            operations,
            adjoint_relative_residual: 0.0,
            adjoint_iterations: 0,
        });
    };
    let (matrix, _) = reduce_matrix_and_lift(&system.operator, &dofs);
    let mut rhs = vec![0.0; dofs.n()];
    rhs[slot] = 1.0;
    let op = CsrOp::symmetric(matrix.clone());
    let pre = crate::solve::spd_preconditioner(&matrix);
    let mut state = CgState::new(&op, &pre, &rhs);
    while state.rel_residual() >= linear.tolerance && state.iters < linear.max_iterations {
        cx.checkpoint().map_err(|_| ConductionError::Cancelled {
            stage: "roundoff-adjoint",
            at: state.iters,
        })?;
        let before = state.iters;
        let batch = (linear.max_iterations - before).min(16);
        state.run(&op, &pre, linear.tolerance, batch);
        if state.iters == before {
            break;
        }
    }
    let mut applied = vec![0.0; rhs.len()];
    matrix.spmv(&state.x, &mut applied);
    let defect: Vec<f64> = applied.iter().zip(&rhs).map(|(a, b)| a - b).collect();
    let adjoint_relative_residual = norm2(&defect);
    if !adjoint_relative_residual.is_finite() || adjoint_relative_residual >= linear.tolerance {
        return Err(ConductionError::Config {
            parameter: "roundoff adjoint",
            what: format!(
                "adjoint residual {adjoint_relative_residual:e} missed tolerance {:e} after {} iterations",
                linear.tolerance, state.iters
            ),
        });
    }
    let mut half_width_k = 0.0_f64;
    for (i, &v) in dofs.free().iter().enumerate() {
        if i % 4096 == 0 {
            cx.checkpoint().map_err(|_| ConductionError::Cancelled {
                stage: "roundoff-rows",
                at: i,
            })?;
        }
        let (columns, values) = system.operator.row(v);
        let magnitude = columns
            .iter()
            .zip(values)
            .fold(system.load[v].abs(), |sum, (&c, &a)| a.abs().mul_add(temperature[c].abs(), sum));
        half_width_k = state.x[i].abs().mul_add(magnitude, half_width_k);
    }
    half_width_k *= gamma;
    if !half_width_k.is_finite() {
        return Err(ConductionError::Config {
            parameter: "roundoff bound",
            what: "the componentwise bound is not finite".to_string(),
        });
    }
    Ok(NodalRoundoffBound {
        half_width_k,
        gamma,
        operations,
        adjoint_relative_residual,
        adjoint_iterations: state.iters,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fixtures::{box_grid, on_box_face};
    use crate::{
        ConductionMesh, ConductivityModel, ScalarField, SolveConfig, ThermalBc,
        ThermalBoundaryBuilder,
    };
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
    use fs_rep_mesh::TetComplex;

    fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            f(&Cx::new(
                &CancelGate::new(),
                arena,
                StreamKey { seed: 7, kernel_id: 907, tile: 0, iteration: 0 },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ))
        })
    }

    fn solve_hottest(cx: &Cx<'_>, reverse: bool) -> (f64, NodalRoundoffBound) {
        let (complex, positions) = box_grid([6, 4, 4], [0.06, 0.04, 0.04]);
        let mut tets = complex.tets.clone();
        if reverse {
            tets.reverse();
        }
        let complex = TetComplex::from_tets(complex.vertex_count, tets);
        let mesh = ConductionMesh::new(complex, positions).unwrap();
        let material = ConductivityModel::isotropic_declared(167.0).unwrap();
        let source = ScalarField::Uniform(2.0e6);
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .region("sink", |f| on_box_face(f.centroid[0], 0.0), ThermalBc::robin(900.0, 300.0).unwrap())
            .unwrap()
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let mut config = SolveConfig::default();
        config.linear.tolerance = 1e-13;
        config.stop.residual_rtol = 1e-12;
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &material,
            element_materials: None,
            source: &source,
        };
        let solution = crate::solve::solve(cx, problem, config.clone()).unwrap();
        let (vertex, &peak) = solution
            .temperature
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.total_cmp(b.1))
            .unwrap();
        let bound =
            nodal_roundoff_bound(cx, problem, None, &solution.temperature, vertex, config.linear)
                .unwrap();
        (peak, bound)
    }

    #[test]
    fn bound_is_tiny_positive_and_covers_an_assembly_order_change() {
        with_cx(|cx| {
            let (forward, bound) = solve_hottest(cx, false);
            let (reversed, reversed_bound) = solve_hottest(cx, true);
            assert!(bound.half_width_k > 0.0 && bound.half_width_k < 1e-6 * forward, "{bound:?}");
            // Reversing the element order changes every summation order in
            // assembly; the peak moves only by rounding, which each run's bound
            // must cover (the two runs' errors add at most).
            let moved = (forward - reversed).abs();
            assert!(
                moved <= bound.half_width_k + reversed_bound.half_width_k,
                "order change moved the peak by {moved:e} K, bounds {:e} + {:e}",
                bound.half_width_k,
                reversed_bound.half_width_k
            );
            println!(
                "{{\"peak_k\":{forward},\"order_shift_k\":{moved:e},\"bound_k\":{:e},\"gamma\":{:e},\"k\":{}}}",
                bound.half_width_k, bound.gamma, bound.operations
            );
        });
    }

    #[test]
    fn a_prescribed_vertex_has_no_roundoff_and_bad_inputs_refuse() {
        with_cx(|cx| {
            let (complex, positions) = box_grid([2, 2, 2], [0.02, 0.02, 0.02]);
            let mesh = ConductionMesh::new(complex, positions).unwrap();
            let material = ConductivityModel::isotropic_declared(10.0).unwrap();
            let source = ScalarField::Uniform(1.0e5);
            let boundary = ThermalBoundaryBuilder::new(&mesh)
                .region("cold", |f| on_box_face(f.centroid[0], 0.0), ThermalBc::dirichlet(300.0).unwrap())
                .unwrap()
                .adiabatic_remainder()
                .finish()
                .unwrap();
            let problem = ConductionProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material,
                element_materials: None,
                source: &source,
            };
            let config = SolveConfig::default();
            let solution = crate::solve::solve(cx, problem, config.clone()).unwrap();
            let fixed = (0..mesh.vertex_count())
                .find(|&v| mesh.positions()[v][0] == 0.0)
                .unwrap();
            let bound =
                nodal_roundoff_bound(cx, problem, None, &solution.temperature, fixed, config.linear)
                    .unwrap();
            assert_eq!(bound.half_width_k, 0.0);
            assert!(
                nodal_roundoff_bound(cx, problem, None, &solution.temperature[1..], 0, config.linear)
                    .is_err()
            );
            assert!(
                nodal_roundoff_bound(
                    cx,
                    problem,
                    None,
                    &solution.temperature,
                    mesh.vertex_count(),
                    config.linear
                )
                .is_err()
            );
        });
    }
}
