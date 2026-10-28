use super::*;
use crate::cohesive_contact::CohesiveLoadLine;
use fs_material::cohesive::{BilinearCohesive, CohesiveViscosity, mixed_mode::MixedModeParameters};

fn config() -> CohesivePatchConfig {
    CohesivePatchConfig {
        max_points: 8,
        contact: CohesiveContactConfig {
            max_iterations: 32,
            max_material_evaluations: 3000,
            max_backtracks: 16,
            jump_absolute_tolerance_m: 1e-14,
            jump_relative_tolerance: 1e-11,
            maximum_jump_increment_m: 0.02,
            maximum_force_n: 100.0,
            material_energy_tolerance_j_m2: 1e-8,
        },
    }
}
fn mixed() -> BilinearCohesive {
    BilinearCohesive::new_mixed_mode(
        MixedModeParameters {
            normal_stiffness_pa_m: 1e9,
            tangential_stiffness_pa_m: 4e8,
            normal_strength_pa: 2e5,
            tangential_strength_pa: 1.6e5,
            mode_i_fracture_energy_j_m2: 100.0,
            mode_ii_fracture_energy_j_m2: 400.0,
            bk_exponent: 2.0,
        },
        2e9,
        0.02,
    )
    .unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}

#[test]
fn heterogeneous_full_cross_compliance_recovers_exact_elastic_solution_with_linear_material_work() {
    let laws = [
        mixed(),
        BilinearCohesive::new(7e8, 3e5, 300.0, 1e9, 0.02).unwrap(),
        mixed(),
    ];
    let old = [[5e-5, -4e-5, 3e-5], [4e-5, 3e-5, -2e-5], [3e-5, 2e-5, 1e-5]];
    let points: Vec<_> = laws
        .iter()
        .zip(old)
        .map(|(law, g)| law.bonded_point().trial(g, 1e-8).unwrap().point)
        .collect();
    let target = [
        [1e-4, 8e-5, -6e-5],
        [8e-5, -7e-5, 4e-5],
        [7e-5, -5e-5, -3e-5],
    ];
    let areas = vec![2e-6, 3e-6, 5e-6];
    let stiffness = [[1e9, 4e8, 4e8], [7e8; 3], [1e9, 4e8, 4e8]];
    let expected: [[f64; 3]; 3] = core::array::from_fn(|p| {
        core::array::from_fn(|i| 0.5 * areas[p] * stiffness[p][i] * (old[p][i] + target[p][i]))
    });
    // Strictly diagonally dominant finite compliance, with every inter-site
    // and inter-axis term present. The second map must retain its asymmetry.
    for asymmetric in [false, true] {
        let compliance: Vec<_> = (0..9)
            .flat_map(|i| {
                (0..9).map(move |j| {
                    if i == j {
                        6e-4 + f64::from(i) * 1e-5
                    } else {
                        1e-5 * f64::from(1 + (i + j) % 3)
                            * if asymmetric && i > j { 0.7 } else { 1.0 }
                    }
                })
            })
            .collect();
        let free: Vec<_> = (0..3)
            .map(|p| {
                core::array::from_fn(|i| {
                    target[p][i]
                        + (0..9)
                            .map(|j| compliance[(3 * p + i) * 9 + j] * expected[j / 3][j % 3])
                            .sum::<f64>()
                })
            })
            .collect();
        let mut solver =
            CohesivePatchLoadLine::new(areas.clone(), compliance.clone(), config()).unwrap();
        let result = solver.trial(&points, &free, 0.01, || true).unwrap();
        for p in 0..3 {
            for i in 0..3 {
                close(result.jump_m[p][i], target[p][i], 3e-14);
                close(result.mean_force_n[p][i], expected[p][i], 3e-10);
            }
        }
        // At most two Newton corrections: N initial + 2*(6N Jacobian + N line).
        // Perturbing/reintegrating the whole patch for every column exceeds it.
        assert!(result.material_evaluations <= 15 * points.len());
        assert_eq!(solver.compliance_m_per_n(), compliance);
        assert_eq!(solver.areas_m2(), areas);
    }
}

#[test]
fn failed_edge_redistributes_load_to_an_intact_neighbor_instead_of_solving_sites_independently() {
    let law = BilinearCohesive::new(1e9, 2e5, 100.0, 2e9, 0.02).unwrap();
    let intact = law.bonded_point();
    let broken = intact
        .trial([0.002, 0.0, 0.0], 1e-8)
        .unwrap()
        .point
        .trial([0.0; 3], 1e-8)
        .unwrap()
        .point;
    assert!(broken.failed());
    let mut s = vec![0.0; 36];
    // Positive translational and rotational fixture compliance for two sites:
    // s=1e-3, c=8e-4. Their local responses cannot be separated.
    for axis in 0..3 {
        s[axis * 6 + axis] = 1e-3;
        s[(3 + axis) * 6 + 3 + axis] = 1e-3;
        s[axis * 6 + 3 + axis] = 8e-4;
        s[(3 + axis) * 6 + axis] = 8e-4;
    }
    let mut solver = CohesivePatchLoadLine::new(vec![1e-6; 2], s, config()).unwrap();
    let free = [[1e-4, 0.0, 0.0]; 2];
    let both = solver
        .trial(&[intact; 2], &free, 0.01, || true)
        .unwrap()
        .mean_force_n[1][0];
    close(both, 1e-4 * 500.0 / (1.0 + 500.0 * 1.8e-3), 1e-12);
    let result = solver
        .trial(&[broken, intact], &free, 0.01, || true)
        .unwrap();
    assert_eq!(result.mean_force_n[0], [0.0; 3]);
    close(
        result.mean_force_n[1][0],
        1e-4 * 500.0 / (1.0 + 500.0 * 1e-3),
        1e-12,
    );
    assert!(result.mean_force_n[1][0] > 1.25 * both);
    close(
        result.jump_m[0][0],
        1e-4 - 8e-4 * result.mean_force_n[1][0],
        1e-14,
    );
}

#[test]
fn a_single_patch_site_recovers_the_existing_load_line_including_physical_viscosity() {
    let law = mixed().with_viscosity(CohesiveViscosity::new(2e6, 7e5, 100.0).unwrap());
    let compliance = [
        [2e-6, 4e-7, -3e-7],
        [4e-7, 3e-6, 5e-7],
        [-3e-7, 5e-7, 1.5e-6],
    ];
    let mut one = CohesiveLoadLine::new(2e-5, compliance, config().contact).unwrap();
    let mut patch = CohesivePatchLoadLine::new(
        vec![2e-5],
        compliance.into_iter().flatten().collect(),
        config(),
    )
    .unwrap();
    let mut point = law.bonded_point();
    for free in [
        [3e-4, 6e-4, -2e-4],
        [-2e-4, -1e-3, 3e-4],
        [5e-5, 2e-5, 0.0],
        [0.001, 0.012, -0.003],
    ] {
        let expected = one.trial(point, free, 0.02, || true).unwrap();
        let actual = patch.trial(&[point], &[free], 0.02, || true).unwrap();
        for i in 0..3 {
            close(
                actual.jump_m[0][i],
                expected.material.equilibrium.point.jump_m()[i],
                2e-12,
            );
            close(actual.mean_force_n[0][i], expected.mean_force_n[i], 1e-8);
        }
        point = expected.material.equilibrium.point;
    }
    assert!(point.failed());
}

#[test]
fn finite_two_body_lever_arms_peel_progressively_and_conserve_force_torque_and_energy() {
    const X: [f64; 3] = [-0.02, 0.0, 0.02];
    let ma = [2.0, 2.0, 2.0, 0.003];
    let mb = [3.0, 3.0, 3.0, 0.004];
    let dt = 1.0 / 5000.0;
    let areas = [1e-6 / 3.0; 3];
    let relative = [0.008, 0.01, -0.005, -2.0];
    let jump = |a: [f64; 4], b: [f64; 4], x: f64| {
        [b[0] - a[0] - x * (b[3] - a[3]), b[1] - a[1], b[2] - a[2]]
    };
    for law in [
        mixed(),
        mixed().with_viscosity(CohesiveViscosity::new(2e6, 7e5, 100.0).unwrap()),
    ] {
        let mut compliance = vec![0.0; 81];
        for p in 0..3 {
            for q in 0..3 {
                for i in 0..3 {
                    compliance[(3 * p + i) * 9 + 3 * q + i] = 0.5
                        * dt
                        * dt
                        * (1.0 / ma[i]
                            + 1.0 / mb[i]
                            + if i == 0 {
                                X[p] * X[q] * (1.0 / ma[3] + 1.0 / mb[3])
                            } else {
                                0.0
                            });
                }
            }
        }
        let mut solver = CohesivePatchLoadLine::new(areas.to_vec(), compliance, config()).unwrap();
        let (mut qa, mut qb) = ([0.0; 4], [0.0; 4]);
        let mut va = core::array::from_fn::<_, 4, _>(|i| -mb[i] / (ma[i] + mb[i]) * relative[i]);
        let mut vb = core::array::from_fn::<_, 4, _>(|i| ma[i] / (ma[i] + mb[i]) * relative[i]);
        let kinetic = |a: [f64; 4], b: [f64; 4]| {
            (0..4)
                .map(|i| 0.5 * (ma[i] * a[i] * a[i] + mb[i] * b[i] * b[i]))
                .sum::<f64>()
        };
        let initial = kinetic(va, vb);
        let mut points = [law.bonded_point(); 3];
        let mut dissipation = 0.0;
        let mut viscous = 0.0;
        let mut edge_ahead = false;
        let mut first_failure = [None; 3];
        for sample in 0..500 {
            let free_a = core::array::from_fn(|i| qa[i] + dt * va[i]);
            let free_b = core::array::from_fn(|i| qb[i] + dt * vb[i]);
            let free = X.map(|x| jump(free_a, free_b, x));
            let result = solver.trial(&points, &free, dt, || true).unwrap();
            let force: [f64; 3] =
                core::array::from_fn(|i| result.mean_force_n.iter().map(|f| f[i]).sum());
            let moment: f64 = result
                .mean_force_n
                .iter()
                .zip(X)
                .map(|(f, x)| x * f[0])
                .sum();
            let qa_force = [force[0], force[1], force[2], -moment];
            let before_a = va;
            let before_b = vb;
            for i in 0..4 {
                qa[i] += dt * va[i] + 0.5 * dt * dt * qa_force[i] / ma[i];
                qb[i] += dt * vb[i] - 0.5 * dt * dt * qa_force[i] / mb[i];
                va[i] += dt * qa_force[i] / ma[i];
                vb[i] -= dt * qa_force[i] / mb[i];
                close(ma[i] * va[i] + mb[i] * vb[i], 0.0, 3e-14);
                close(ma[i] * (va[i] - before_a[i]), dt * qa_force[i], 3e-14);
                close(mb[i] * (vb[i] - before_b[i]), -dt * qa_force[i], 3e-14);
            }
            let mut work = 0.0;
            for p in 0..3 {
                let actual_jump = jump(qa, qb, X[p]);
                let (actual, mean) = points[p]
                    .trial_over_time_with_mean_traction(actual_jump, dt, 1e-8)
                    .unwrap();
                for i in 0..3 {
                    close(actual_jump[i], result.jump_m[p][i], 5e-13);
                    close(areas[p] * mean[i], result.mean_force_n[p][i], 2e-9);
                }
                work += dot(
                    result.mean_force_n[p],
                    core::array::from_fn(|i| actual_jump[i] - points[p].jump_m()[i]),
                );
                dissipation += areas[p]
                    * (actual.equilibrium.energy.dissipated + actual.viscous_dissipated_j_m2);
                viscous += areas[p] * actual.viscous_dissipated_j_m2;
                assert!(actual.equilibrium.point.damage() >= points[p].damage());
                points[p] = actual.equilibrium.point;
                if points[p].failed() && first_failure[p].is_none() {
                    first_failure[p] = Some(sample);
                }
            }
            close(
                kinetic(va, vb) - kinetic(before_a, before_b) + work,
                0.0,
                3e-15,
            );
            let stored: f64 = points
                .iter()
                .zip(areas)
                .map(|(p, a)| a * p.stored_energy_j_m2())
                .sum();
            close(kinetic(va, vb) + stored + dissipation, initial, 3e-11);
            edge_ahead |= sample < 250
                && points[2].damage() > points[1].damage() + 0.02
                && points[2].damage() > points[0].damage() + 0.05;
        }
        assert!(edge_ahead);
        assert!(
            first_failure[2].is_some(),
            "the tensile edge must separate completely: {first_failure:?}"
        );
        assert!(first_failure[1].is_none_or(|s| s > first_failure[2].unwrap()));
        assert!(first_failure[0].is_none_or(|s| s > first_failure[2].unwrap()));
        assert_eq!(viscous > 0.0, law.viscosity().is_some());
    }
}

#[test]
fn cancellation_partial_point_refusal_and_global_evaluation_budget_leave_all_histories_unchanged() {
    let ordinary = mixed();
    let rate = mixed().with_viscosity(CohesiveViscosity::new(2e6, 7e5, 0.02).unwrap());
    let points = [ordinary.bonded_point(), rate.bonded_point()];
    let saved = points;
    let mut owner = CohesivePatchLoadLine::new(vec![1e-6; 2], vec![0.0; 36], config()).unwrap();
    let free = [[1e-5, 2e-5, 0.0]; 2];
    assert_eq!(
        owner.trial(&points, &free, 0.01, || false),
        Err(CohesivePatchError::Cancelled)
    );
    let mut polls = 0;
    assert_eq!(
        owner.trial(&points, &free, 0.01, || {
            polls += 1;
            polls < 8
        }),
        Err(CohesivePatchError::Cancelled)
    );
    let mut short = config();
    short.contact.max_material_evaluations = 2;
    let mut limited = CohesivePatchLoadLine::new(vec![1e-6; 2], vec![0.0; 36], short).unwrap();
    assert!(matches!(
        limited.trial(&points, &free, 0.01, || true),
        Err(CohesivePatchError::Budget(_))
    ));
    assert!(
        owner
            .trial(&points, &[[0.001, 0.0, 0.0]; 2], 0.01, || true)
            .is_err()
    );
    assert_eq!(points, saved);
    let expected = owner
        .trial(&points, &free, 0.01, || true)
        .unwrap()
        .jump_m
        .to_vec();
    let replay = owner.trial(&saved, &free, 0.01, || true).unwrap();
    assert_eq!(replay.jump_m, expected);
    assert!(replay.material_evaluations <= config().contact.max_material_evaluations);
    assert!(
        CohesivePatchLoadLine::new(vec![1e-6; MAX_COHESIVE_PATCH_POINTS + 1], vec![], config())
            .is_err()
    );
    assert!(CohesivePatchLoadLine::new(vec![1e-6; 2], vec![0.0; 9], config()).is_err());
}
