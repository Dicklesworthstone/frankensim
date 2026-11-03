//! G1/G3: independent simplex pressure projection and one retained factor.

use super::*;
use crate::linear3::TetAssemblyBudget;
use crate::linear3::tet10::{TET10_EDGES, promote_p1_to_p2};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_solver::CheckedCgConfig;
use std::cell::Cell;

const IDENTITY: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];

fn with_cx(run: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 71,
                kernel_id: 12,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}

fn controls() -> MixedImpulseControls {
    MixedImpulseControls {
        balance: SolidImpulseControls {
            solve: CheckedCgConfig {
                tolerance: 1e-11,
                max_iterations: 64,
                max_corrections: 2,
            },
            maximum_workspace_components: 4096,
            linear_tolerance_kg_m_s: 1e-10,
            angular_tolerance_kg_m2_s: 1e-10,
            work_tolerance_j: 1e-10,
        },
        velocity_scale_m_s: 1.0,
        pressure_impulse_scale_pa_s: 1.0,
        work_scale_j: 1.0,
        maximum_assembly_bytes: 1 << 20,
        maximum_factor_bytes: 1 << 20,
        maximum_relative_volume_rate_per_s: 1e-9,
    }
}

fn specimen(cx: &Cx<'_>) -> (Vec<[f64; 3]>, PreparedMixedMesh, ReferenceMass) {
    let (nodes, elements) = promote_p1_to_p2(
        &[[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        &[[0, 1, 2, 3]],
        None,
        cx,
    )
    .unwrap();
    let mesh = PreparedMixedMesh::new(&nodes, &elements, TetAssemblyBudget::standard(), 4, 64, cx)
        .unwrap();
    let mass = ReferenceMass::from_tet10_quadrature(
        mesh.elements.iter().map(|element| {
            (
                element.connectivity,
                element.points.iter().map(|q| (q.shapes, q.weight_m3)),
            )
        }),
        &[6.0], // unit total mass on the unit reference tetrahedron
        30,
        1 << 20,
        cx,
    )
    .unwrap();
    (nodes, mesh, mass)
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-10 * expected.abs().max(1.0),
        "{actual:.17e} != {expected:.17e}"
    );
}

fn point_load(cell: [usize; 10], weights: [f64; 10], axis: usize) -> Vec<(usize, f64)> {
    cell.into_iter()
        .zip(weights)
        .map(|(node, w)| (3 * node + axis, w))
        .collect()
}

#[test]
fn one_factor_reuses_exact_constrained_point_response_and_original_event_account() {
    with_cx(|cx| {
        let (nodes, mesh, mass) = specimen(cx);
        let samples = Cell::new(0);
        let response = ConstrainedImpulseResponse::prepare(
            &mesh,
            &mass,
            &[],
            None,
            controls(),
            |_, _| {
                samples.set(samples.get() + 1);
                Ok(IDENTITY)
            },
            cx,
        )
        .unwrap();
        assert_eq!(samples.get(), 64);
        assert_eq!(
            response.kkt_inertia(),
            Inertia {
                positive: 30,
                negative: 4
            }
        );
        let center = [
            -0.125, -0.125, -0.125, -0.125, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
        ];
        let mut vertex = [0.0; 10];
        vertex[1] = 1.0;
        let cell = mesh.elements[0].connectivity;
        let loads = [
            point_load(cell, center, 0),
            point_load(cell, center, 1),
            point_load(cell, vertex, 0),
        ];
        // Independent rational integration of N_a*N_b and L_p*grad(N_a),
        // followed by exact elimination of the 30+4 saddle system. These
        // entries include pressure-induced cross-axis and cross-point response.
        let exact = [
            [17735.0 / 6528.0, -1505.0 / 13056.0, -1895.0 / 408.0],
            [-1505.0 / 13056.0, 17735.0 / 6528.0, 35.0 / 816.0],
            [-1895.0 / 408.0, 35.0 / 816.0, 3670.0 / 51.0],
        ];
        for (column, load) in loads.iter().enumerate() {
            let (velocity, solve) = response.free_velocity_response(load, cx).unwrap();
            assert!(solve.converged_euclidean());
            for (row, trace) in loads.iter().enumerate() {
                let value: f64 = trace.iter().map(|&(d, w)| w * velocity[d / 3][d % 3]).sum();
                near(value, exact[row][column]);
            }
        }
        // Mass-only center mobility is 95/32. It is a different impact law.
        assert!(exact[0][0] < 95.0 / 32.0);
        let (increment, _) = response.free_velocity_response(&loads[0], cx).unwrap();
        let before = vec![[0.25, -0.125, 0.0625]; nodes.len()];
        let (after, report) = response
            .impulse_candidate(&nodes, &before, &loads[0], [0.0; 3], cx)
            .unwrap();
        for ((old, change), next) in before.iter().zip(&increment).zip(&after) {
            for axis in 0..3 {
                near(next[axis], old[axis] + change[axis]);
            }
        }
        near(report.external_work_j, 0.25 + 0.5 * exact[0][0]);
        near(report.momentum.kinetic_change_j, report.external_work_j);
        near(report.pressure_work_j, 0.0);
        assert_eq!(report.kkt_inertia, Some(response.kkt_inertia()));
        assert_eq!(samples.get(), 64); // no endpoint material/cofactor reevaluation
    });
}

#[test]
fn exact_pressure_modes_have_no_mobility_and_incompatible_background_refuses() {
    with_cx(|cx| {
        let (nodes, mesh, mass) = specimen(cx);
        let response = ConstrainedImpulseResponse::prepare(
            &mesh,
            &mass,
            &[],
            None,
            controls(),
            |_, _| Ok(IDENTITY),
            cx,
        )
        .unwrap();
        let gradients = [[-1.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let multipliers = [0.5, -0.25, 0.75, -1.0];
        let mut loads = Vec::new();
        for (a, &node) in mesh.elements[0].connectivity.iter().enumerate() {
            for axis in 0..3 {
                // Exact reference integral B_pa = integral L_p grad(N_a).
                // No production quadrature or constraint entries form this load.
                let value = (0..4)
                    .map(|p| {
                        let derivative = if a < 4 {
                            (if p == a { 3.0 } else { -1.0 }) * gradients[a][axis] / 120.0
                        } else {
                            let [i, j] = TET10_EDGES[a - 4];
                            ((if p == i { 2.0 } else { 1.0 }) * gradients[j][axis]
                                + (if p == j { 2.0 } else { 1.0 }) * gradients[i][axis])
                                / 30.0
                        };
                        multipliers[p] * derivative
                    })
                    .sum();
                loads.push((3 * node + axis, value));
            }
        }
        let (increment, _) = response.free_velocity_response(&loads, cx).unwrap();
        for &v in increment.iter().flatten() {
            near(v, 0.0);
        }
        let before = vec![[0.125, -0.25, 0.5]; nodes.len()];
        let (after, report) = response
            .impulse_candidate(&nodes, &before, &loads, [0.0; 3], cx)
            .unwrap();
        for (v, expected) in after.iter().zip(&before) {
            for axis in 0..3 {
                near(v[axis], expected[axis]);
            }
        }
        for (&actual, expected) in report.pressure_impulses_pa_s.iter().zip(multipliers) {
            near(actual, expected);
        }
        near(report.external_work_j, 0.0);
        near(report.pressure_work_j, 0.0);
        let expansion: Vec<_> = nodes.iter().map(|x| x.map(|v| 0.02 * v)).collect();
        assert!(response.validate_velocity(&expansion, cx).is_err());
        assert!(
            response
                .impulse_candidate(&nodes, &expansion, &loads, [0.0; 3], cx)
                .is_err()
        );
        // A failed state candidate cannot disturb the factored response.
        assert_eq!(
            response.free_velocity_response(&loads, cx).unwrap().0,
            increment
        );
    });
}
