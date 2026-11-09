//! G1/G3 implicit truss design derivatives, including imposed-state load terms.

use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_truss::elastic::sensitivity::displacement_functional;
use fs_truss::elastic::{
    ElasticControls, ElasticError, ElasticImposedState, solve, solve_with_imposed,
};
use fs_truss::{GroundLimits, GroundStructure, LayoutCase, TrussConstructionError};

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let gate = CancelGate::new();
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 1,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate);
    });
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-20),
        "{actual:e} != {expected:e}",
    );
}

fn star(cx: &Cx<'_>, angle: f64, force: [f64; 2]) -> (GroundStructure, LayoutCase) {
    let (s, c) = angle.sin_cos();
    let rotate = |p: [f64; 2]| [c * p[0] - s * p[1], s * p[0] + c * p[1]];
    let ground = GroundStructure::try_from_parts(
        &[[-1.0, 0.0], [0.0, -1.0], [-1.0, -1.0], [0.0, 0.0]].map(rotate),
        &[(0, 3), (1, 3), (2, 3)],
        &[1.0, 1.0, 2.0_f64.sqrt()],
        GroundLimits::default(),
        cx,
    )
    .unwrap();
    let case = LayoutCase::try_new(
        vec![[true; 2], [true; 2], [true; 2], [false; 2]],
        vec![[0.0; 2], [0.0; 2], [0.0; 2], rotate(force)],
        4,
    )
    .unwrap();
    (ground, case)
}

fn objective(
    ground: &GroundStructure,
    case: &LayoutCase,
    areas: &[f64],
    young: &[f64],
    imposed: ElasticImposedState<'_>,
    weights: &[[f64; 2]],
    cx: &Cx<'_>,
) -> f64 {
    let response = solve_with_imposed(
        ground,
        case,
        areas,
        young,
        imposed,
        ElasticControls::default(),
        cx,
    )
    .unwrap();
    weights
        .iter()
        .zip(&response.displacements_m)
        .map(|(g, u)| g[0] * u[0] + g[1] * u[1])
        .sum()
}

#[test]
fn g1_settlement_compliance_has_the_correct_positive_area_derivative() {
    with_cx(|cx, _| {
        let ground = GroundStructure::try_from_parts(
            &[[0.0, 0.0], [2.0, 0.0], [3.0, 0.0]],
            &[(0, 1), (1, 2)],
            &[2.0, 1.0],
            GroundLimits::default(),
            cx,
        )
        .unwrap();
        let case = LayoutCase::try_new(
            vec![[true, true], [false, true], [true, true]],
            vec![[2.0, 0.0], [10.0, 0.0], [-3.0, 0.0]],
            3,
        )
        .unwrap();
        let result = displacement_functional(
            &ground,
            &case,
            &[2.0, 3.0],
            &[4e6, 2e6],
            ElasticImposedState {
                member_stress_free_strain: Some(&[1e-6, -2e-6]),
                prescribed_displacements_m: Some(&[[-2e-6, 0.0], [0.0; 2], [8e-6, 0.0]]),
            },
            case.loads(),
            ElasticControls::default(),
            cx,
        )
        .unwrap();
        // K_ff = 10e6; u1=7e-6, lambda1=1e-6; N=[28,18].
        // The right member's elastic and adjoint extensions have opposite signs.
        near(result.value, 42e-6, 1e-12);
        let left = result.members[0].unwrap();
        let right = result.members[1].unwrap();
        near(left.area_derivative, -14e-6, 1e-12);
        near(right.area_derivative, 6e-6, 1e-12);
        near(left.young_modulus_derivative, -7e-12, 1e-12);
        near(right.young_modulus_derivative, 9e-12, 1e-12);
        near(left.stress_free_strain_derivative, 8.0, 1e-12);
        near(right.stress_free_strain_derivative, -6.0, 1e-12);
        near(result.prescribed_displacement_gradient[0][0], 6.0, 1e-12);
        near(result.prescribed_displacement_gradient[2][0], 3.0, 1e-12);
        near(result.load_gradient[1][0], 1e-6, 1e-12);
        assert_eq!(result.load_gradient[0], [0.0; 2]);
        assert_eq!(result.prescribed_displacement_gradient[1][0], 0.0);
    });
}

#[test]
#[allow(clippy::too_many_lines)] // One fixture, independent differences of each physical input.
fn g3_redundant_thermal_and_settlement_gradients_match_centered_differences() {
    with_cx(|cx, _| {
        for angle in [0.0_f64, 0.37] {
            let (ground, case) = star(cx, angle, [-38.0, -63.0]);
            let (s, c) = angle.sin_cos();
            let rotate = |p: [f64; 2]| [c * p[0] - s * p[1], s * p[0] + c * p[1]];
            let areas = [1.0; 3];
            let young = [2e6, 3e6, 6e6 * 2.0_f64.sqrt()];
            let strains = [2e-6, -1e-6, 3e-6];
            let motion = [[1e-6, -2e-6], [-3e-6, 1e-6], [2e-6, 4e-6], [0.0; 2]].map(rotate);
            let weights = [[1.0, -2.0], [0.0; 2], [0.0; 2], [4.0, -3.0]].map(rotate);
            let state = ElasticImposedState {
                member_stress_free_strain: Some(&strains),
                prescribed_displacements_m: Some(&motion),
            };
            let result = displacement_functional(
                &ground,
                &case,
                &areas,
                &young,
                state,
                &weights,
                ElasticControls::default(),
                cx,
            )
            .unwrap();
            near(result.value, 46e-6, 1e-12);
            for (index, expected) in [-44e-6 / 7.0, -27e-6, 12e-6].into_iter().enumerate() {
                near(
                    result.members[index].unwrap().area_derivative,
                    expected,
                    1e-11,
                );
            }
            for index in 0..3 {
                let derivative = result.members[index].unwrap();
                let mut plus = areas;
                let mut minus = areas;
                let step = 2e-5 * areas[index];
                plus[index] += step;
                minus[index] -= step;
                let difference = (objective(&ground, &case, &plus, &young, state, &weights, cx)
                    - objective(&ground, &case, &minus, &young, state, &weights, cx))
                    / (2.0 * step);
                near(derivative.area_derivative, difference, 2e-7);

                let mut plus = young;
                let mut minus = young;
                let step = 2e-5 * young[index];
                plus[index] += step;
                minus[index] -= step;
                let difference = (objective(&ground, &case, &areas, &plus, state, &weights, cx)
                    - objective(&ground, &case, &areas, &minus, state, &weights, cx))
                    / (2.0 * step);
                near(derivative.young_modulus_derivative, difference, 2e-7);

                let mut plus = strains;
                let mut minus = strains;
                let step = 1e-8;
                plus[index] += step;
                minus[index] -= step;
                let difference = (objective(
                    &ground,
                    &case,
                    &areas,
                    &young,
                    ElasticImposedState {
                        member_stress_free_strain: Some(&plus),
                        ..state
                    },
                    &weights,
                    cx,
                ) - objective(
                    &ground,
                    &case,
                    &areas,
                    &young,
                    ElasticImposedState {
                        member_stress_free_strain: Some(&minus),
                        ..state
                    },
                    &weights,
                    cx,
                )) / (2.0 * step);
                near(derivative.stress_free_strain_derivative, difference, 2e-7);
            }
            for node in 0..3 {
                for axis in 0..2 {
                    let mut plus = motion;
                    let mut minus = motion;
                    let step = 1e-8;
                    plus[node][axis] += step;
                    minus[node][axis] -= step;
                    let difference = (objective(
                        &ground,
                        &case,
                        &areas,
                        &young,
                        ElasticImposedState {
                            prescribed_displacements_m: Some(&plus),
                            ..state
                        },
                        &weights,
                        cx,
                    ) - objective(
                        &ground,
                        &case,
                        &areas,
                        &young,
                        ElasticImposedState {
                            prescribed_displacements_m: Some(&minus),
                            ..state
                        },
                        &weights,
                        cx,
                    )) / (2.0 * step);
                    near(
                        result.prescribed_displacement_gradient[node][axis],
                        difference,
                        2e-7,
                    );
                }
            }
            for axis in 0..2 {
                let mut plus = case.loads().to_vec();
                let mut minus = plus.clone();
                let step = 1e-3;
                plus[3][axis] += step;
                minus[3][axis] -= step;
                let plus = LayoutCase::try_new(case.supported().to_vec(), plus, 4).unwrap();
                let minus = LayoutCase::try_new(case.supported().to_vec(), minus, 4).unwrap();
                let difference = (objective(&ground, &plus, &areas, &young, state, &weights, cx)
                    - objective(&ground, &minus, &areas, &young, state, &weights, cx))
                    / (2.0 * step);
                near(result.load_gradient[3][axis], difference, 2e-7);
            }
            assert!(result.adjoint_relative_equilibrium_residual < 1e-12);
        }
    });
}

#[test]
fn g0_removed_members_have_no_topology_activation_derivative() {
    with_cx(|cx, _| {
        let (ground, case) = star(cx, 0.0, [11.0, 15.0]);
        let result = displacement_functional(
            &ground,
            &case,
            &[1.0, 1.0, 0.0],
            &[2e6, 3e6, 1e9],
            ElasticImposedState::default(),
            case.loads(),
            ElasticControls::default(),
            cx,
        )
        .unwrap();
        assert!(result.members[0].is_some());
        assert!(result.members[1].is_some());
        assert_eq!(result.members[2], None);
        near(result.members[0].unwrap().area_derivative, -60.5e-6, 1e-12);
        near(result.members[1].unwrap().area_derivative, -75e-6, 1e-12);
    });
}

#[test]
fn g4_invalid_objectives_unsupported_adjoint_shared_budget_and_cancellation_refuse() {
    with_cx(|cx, gate| {
        let (ground, case) = star(cx, 0.0, [11.0, 15.0]);
        let areas = [1.0; 3];
        let young = [2e6, 3e6, 6e6 * 2.0_f64.sqrt()];
        let state = ElasticImposedState::default();
        let controls = ElasticControls::default();
        assert!(
            displacement_functional(
                &ground,
                &case,
                &areas,
                &young,
                state,
                &[[0.0; 2]; 3],
                controls,
                cx,
            )
            .is_err()
        );
        assert!(
            displacement_functional(
                &ground,
                &case,
                &areas,
                &young,
                state,
                &[[0.0; 2], [0.0; 2], [0.0; 2], [f64::NAN, 0.0]],
                controls,
                cx,
            )
            .is_err()
        );

        let primal = solve(&ground, &case, &areas, &young, controls, cx).unwrap();
        let too_small = ElasticControls {
            maximum_iterations: primal.iterations,
            ..controls
        };
        assert!(matches!(
            displacement_functional(
                &ground,
                &case,
                &areas,
                &young,
                state,
                case.loads(),
                too_small,
                cx,
            ),
            Err(ElasticError::NotConverged { .. })
        ));

        let (_, unloaded_y) = star(cx, 0.0, [11.0, 0.0]);
        assert!(matches!(
            displacement_functional(
                &ground,
                &unloaded_y,
                &[1.0, 0.0, 0.0],
                &young,
                state,
                &[[0.0; 2], [0.0; 2], [0.0; 2], [0.0, 1.0]],
                controls,
                cx,
            ),
            Err(ElasticError::UnsupportedLoad { node: 3, axis: 1 })
        ));
        gate.request();
        assert!(matches!(
            displacement_functional(
                &ground,
                &case,
                &areas,
                &young,
                state,
                case.loads(),
                controls,
                cx,
            ),
            Err(ElasticError::Construction(
                TrussConstructionError::Cancelled { .. }
            ))
        ));
    });
}
