//! G1/G3: physical quadratic inertia, retained nonlinear state and event work.

use super::*;
use fs_solid::hyper3::mixed::conservative::Tet10ImpulsePoint;
use fs_solid::hyper3::mixed::impulse::{MixedImpulseControls, SolidImpulseControls};
use fs_solver::CheckedCgConfig;

fn controls() -> MixedImpulseControls {
    MixedImpulseControls {
        balance: SolidImpulseControls {
            solve: CheckedCgConfig {
                tolerance: 1e-13,
                max_iterations: 64,
                max_corrections: 2,
            },
            maximum_workspace_components: 4096,
            linear_tolerance_kg_m_s: 1e-11,
            angular_tolerance_kg_m2_s: 1e-11,
            work_tolerance_j: 1e-10,
        },
        // No pressure constraint or sparse KKT factor is constructed here.
        velocity_scale_m_s: 0.0,
        pressure_impulse_scale_pa_s: 0.0,
        work_scale_j: 0.0,
        maximum_assembly_bytes: 0,
        maximum_factor_bytes: 0,
        maximum_relative_volume_rate_per_s: 0.0,
    }
}

#[test]
fn full_and_original_bulk_impulses_keep_state_and_continue_the_original_work_account() {
    with_cx(false, |cx| {
        let (nodes, elements) = mesh(cx);
        let material = law(false, KAPPA);
        let densities = [DENSITY];
        let (u0, v0) = initial_state(&nodes);
        for full in [false, true] {
            let mut source = definition(
                &nodes,
                &elements,
                MixedLawField::Uniform(&material),
                &densities,
            );
            let pressure = initial_pressure(false);
            if full {
                source.formulation = ConservativeFormulation::FullMaterial;
            }
            let solid = source
                .prepare_with_pressure(&u0, &v0, if full { &[] } else { &pressure }, cx)
                .unwrap();
            let before = state_bits(&solid);
            let delta: Vec<[f64; 3]> = nodes
                .iter()
                .zip(&u0)
                .map(|(x, u)| {
                    let spin = cross([0.03, -0.02, 0.04], core::array::from_fn(|c| x[c] + u[c]));
                    core::array::from_fn(|c| [0.02, -0.01, 0.03][c] + spin[c])
                })
                .collect();
            // Exact independent P2 simplex integrals, including signed corner
            // couplings. No assembled production matrix supplies the oracle.
            let impulses: Vec<_> = elements[0]
                .iter()
                .enumerate()
                .flat_map(|(a, &node)| {
                    let delta = &delta;
                    let element = &elements[0];
                    (0..3).map(move |c| {
                        (
                            3 * node + c,
                            (0..10)
                                .map(|b| {
                                    DENSITY * VOLUME * mass_fraction(a, b) * delta[element[b]][c]
                                })
                                .sum::<f64>(),
                        )
                    })
                })
                .collect();
            let event = solid
                .impulse_candidate(&[], &impulses, None, [0.0; 3], controls(), cx)
                .unwrap();
            let expected_v: Vec<[f64; 3]> = v0
                .iter()
                .zip(&delta)
                .map(|(v, d)| core::array::from_fn(|c| v[c] + d[c]))
                .collect();
            let old = independent_motion(&nodes, &elements[0], &u0, &v0);
            let expected = independent_motion(&nodes, &elements[0], &u0, &expected_v);
            for (actual, wanted) in event.solid.velocity_mps().iter().zip(&expected_v) {
                for c in 0..3 {
                    near(actual[c], wanted[c], 2e-12);
                }
            }
            for c in 0..3 {
                near(
                    event.impulse.momentum.after.linear_kg_m_s[c],
                    expected.1[c],
                    2e-12,
                );
                near(
                    event.impulse.momentum.after.angular_kg_m2_s[c],
                    expected.2[c],
                    2e-12,
                );
                near(
                    event.impulse.external.linear_kg_m_s[c],
                    expected.1[c] - old.1[c],
                    2e-12,
                );
                near(
                    event.impulse.external.angular_kg_m2_s[c],
                    expected.2[c] - old.2[c],
                    2e-12,
                );
            }
            near(
                event.impulse.momentum.kinetic_change_j,
                expected.0 - old.0,
                2e-12,
            );
            near(event.impulse.external_work_j, expected.0 - old.0, 2e-12);
            assert!(event.impulse.momentum.support_impulses_kg_m_s.is_empty());
            assert_eq!(event.impulse.kkt_inertia, None);
            assert!(event.impulse.pressure_nodes.is_empty());
            assert!(event.impulse.pressure_impulses_pa_s.is_empty());
            assert!(event.impulse.volume_rate_before_m3_s.is_empty());
            assert!(event.impulse.volume_rate_after_m3_s.is_empty());
            assert_eq!(event.impulse.pressure_work_j, 0.0);
            assert_eq!(event.impulse.pressure.linear_kg_m_s, [0.0; 3]);
            assert_eq!(event.impulse.pressure.angular_kg_m2_s, [0.0; 3]);
            assert_eq!(event.solid.displacement_m(), solid.displacement_m());
            assert_eq!(
                event.solid.algorithmic_pressure_pa(),
                solid.algorithmic_pressure_pa()
            );
            assert_eq!(event.solid.bulk_volume_ratios(), solid.bulk_volume_ratios());
            assert_eq!(event.solid.stored_energy_j(), solid.stored_energy_j());
            assert_eq!(
                event.solid.bulk_stored_energy_j(),
                solid.bulk_stored_energy_j()
            );
            assert_eq!(
                event.solid.reference_density_kg_m3(),
                solid.reference_density_kg_m3()
            );
            assert_eq!(event.solid.time_s(), solid.time_s());
            assert_eq!(state_bits(&solid), before);
            near(
                event.impulse.momentum.work_residual_j,
                0.0,
                controls().balance.work_tolerance_j,
            );

            // The event is part of this trajectory's original work history.
            // A native nonlinear step must not see it as unexplained energy.
            let next = event
                .solid
                .advance(0.002, &[], &[], [0.0; 3], settings(), WORK_ALLOWANCE, cx)
                .unwrap();
            near(next.energy.cumulative_balance_error_j, 0.0, WORK_ALLOWANCE);
            near(
                next.solid.total_energy_j() - solid.total_energy_j(),
                event.impulse.external_work_j,
                WORK_ALLOWANCE,
            );
            assert_ne!(next.solid.displacement_m(), event.solid.displacement_m());
        }
    });
}

#[test]
fn signed_quadratic_support_impulses_retain_moving_support_work() {
    with_cx(false, |cx| {
        let (nodes, elements) = mesh(cx);
        let material = law(false, KAPPA);
        let densities = [DENSITY];
        let u0 = vec![[0.0; 3]; nodes.len()];
        let v0 = vec![[0.0, 0.1, 0.0]; nodes.len()];
        let source = ConservativeTetDefinition {
            formulation: ConservativeFormulation::FullMaterial,
            ..definition(
                &nodes,
                &elements,
                MixedLawField::Uniform(&material),
                &densities,
            )
        };
        let solid = source.prepare(&u0, &v0, cx).unwrap();
        let node = elements[0][4];
        let free = 3 * node + 1;
        let held: Vec<_> = (0..3 * nodes.len()).filter(|d| *d != free).collect();
        let event = solid
            .impulse_candidate(
                &held,
                &[(free, 0.012), (free, 0.018)],
                None,
                [0.0; 3],
                controls(),
                cx,
            )
            .unwrap();
        let delta = 0.03 / (DENSITY * VOLUME * mass_fraction(4, 4));
        let mut expected_v = v0.clone();
        expected_v[node][1] += delta;
        let mut support_work = 0.0;
        for &(dof, impulse) in &event.impulse.momentum.support_impulses_kg_m_s {
            let a = elements[0].iter().position(|&n| n == dof / 3).unwrap();
            let expected = if dof % 3 == 1 {
                DENSITY * VOLUME * mass_fraction(a, 4) * delta
            } else {
                0.0
            };
            near(impulse, expected, 2e-13);
            support_work += v0[dof / 3][dof % 3] * expected;
        }
        assert!(support_work > 0.0);
        near(event.impulse.momentum.support_work_j, support_work, 2e-13);
        near(
            event.impulse.external_work_j,
            0.03 * (0.1 + 0.5 * delta),
            2e-13,
        );
        let expected = independent_motion(&nodes, &elements[0], &u0, &expected_v);
        near(event.solid.kinetic_energy_j(), expected.0, 2e-13);
        for dof in held {
            assert_eq!(
                event.solid.velocity_mps()[dof / 3][dof % 3],
                v0[dof / 3][dof % 3]
            );
        }
        near(
            event.impulse.momentum.work_residual_j,
            0.0,
            controls().balance.work_tolerance_j,
        );
    });
}

#[test]
fn impulse_refusals_keep_the_accepted_solid_and_pressure_gauges_require_constraints() {
    with_cx(false, |cx| {
        let (nodes, elements) = mesh(cx);
        let material = law(false, KAPPA);
        let densities = [DENSITY];
        let zero = vec![[0.0; 3]; nodes.len()];
        let source = definition(
            &nodes,
            &elements,
            MixedLawField::Uniform(&material),
            &densities,
        );
        let solid = source.prepare(&zero, &zero, cx).unwrap();
        let before = state_bits(&solid);
        let mut limited = controls();
        limited.balance.maximum_workspace_components = 0;
        assert!(
            solid
                .impulse_candidate(&[], &[(0, 0.1)], None, [0.0; 3], limited, cx)
                .is_err()
        );
        assert!(
            solid
                .impulse_candidate(&[1, 1], &[], None, [0.0; 3], controls(), cx)
                .is_err()
        );
        assert!(
            solid
                .impulse_candidate(&[], &[(0, f64::NAN)], None, [0.0; 3], controls(), cx)
                .is_err()
        );
        with_cx(true, |stopped| {
            assert!(
                solid
                    .impulse_candidate(&[], &[(0, 0.1)], None, [0.0; 3], controls(), stopped)
                    .is_err()
            );
        });
        assert_eq!(state_bits(&solid), before);
        assert!(
            solid
                .impulse_candidate(&[], &[], Some(0), [0.0; 3], controls(), cx)
                .is_err()
        );
        assert_eq!(state_bits(&solid), before);
    });
}

#[test]
fn point_response_uses_each_original_formulation_and_continues_its_cumulative_work() {
    with_cx(false, |cx| {
        let (nodes, elements) = mesh(cx);
        let material = law(false, KAPPA);
        let density = [DENSITY];
        let zero = vec![[0.0; 3]; nodes.len()];
        let velocity = vec![[0.02, -0.03, 0.01]; nodes.len()];
        let center = Tet10ImpulsePoint {
            tetrahedron: 0,
            barycentric: [0.25; 4],
        };
        let corner = Tet10ImpulsePoint {
            tetrahedron: 0,
            barycentric: [0.0, 1.0, 0.0, 0.0],
        };
        let rows = [
            (center, [1.0, 0.0, 0.0]),
            (center, [0.0, 1.0, 0.0]),
            (corner, [1.0, 0.0, 0.0]),
        ];
        let mut source = definition(
            &nodes,
            &elements,
            MixedLawField::Uniform(&material),
            &density,
        );
        let original = source.formulation;
        for formulation in [
            ConservativeFormulation::FullMaterial,
            original,
            ConservativeFormulation::ExactIncompressible {
                maximum_relative_volume_residual: 1e-9,
                maximum_initial_relative_volume_rate_per_s: 1e-9,
                maximum_point_volume_work_defect: 1e-12,
            },
        ] {
            source.formulation = formulation;
            let solid = source.prepare(&zero, &velocity, cx).unwrap();
            let before = state_bits(&solid);
            let exact_mode = matches!(
                formulation,
                ConservativeFormulation::ExactIncompressible { .. }
            );
            let mut control = controls();
            if exact_mode {
                control.balance.solve.tolerance = 1e-11;
                control.velocity_scale_m_s = 1.0;
                control.pressure_impulse_scale_pa_s = 1.0;
                control.work_scale_j = 1.0;
                control.maximum_assembly_bytes = 1 << 20;
                control.maximum_factor_bytes = 1 << 20;
                control.maximum_relative_volume_rate_per_s = 1e-9;
            }
            let response = solid
                .point_impulse_response(&[], None, control, cx)
                .unwrap();
            let mobility = response.point_mobility_matrix_per_kg(&rows, 3, cx).unwrap();
            // Independent rational P2 mass and P1 volume integrals. The exact
            // pressure response has physical cross-axis terms; the original
            // finite-compliance inertia has independent Cartesian components.
            let expected = if exact_mode {
                [
                    17735.0 / 6528.0,
                    -1505.0 / 13056.0,
                    -1895.0 / 408.0,
                    -1505.0 / 13056.0,
                    17735.0 / 6528.0,
                    35.0 / 816.0,
                    -1895.0 / 408.0,
                    35.0 / 816.0,
                    3670.0 / 51.0,
                ]
            } else {
                [
                    95.0 / 32.0,
                    0.0,
                    -55.0 / 8.0,
                    0.0,
                    95.0 / 32.0,
                    0.0,
                    -55.0 / 8.0,
                    0.0,
                    100.0,
                ]
            };
            for (&actual, wanted) in mobility.matrix_per_kg.iter().zip(expected) {
                near(actual, wanted, 2e-9);
            }
            assert_eq!(response.kkt_inertia().is_some(), exact_mode);
            assert!(response.point_mobility_matrix_per_kg(&rows, 2, cx).is_err());
            let coefficients = [0.002, -0.001, -0.001];
            let event = response
                .point_impulses_candidate(
                    &[(center, [0.002, -0.001, 0.0]), (corner, [-0.001, 0.0, 0.0])],
                    [0.0; 3],
                    cx,
                )
                .unwrap();
            let quadratic_work: f64 = (0..3)
                .flat_map(|i| {
                    (0..3)
                        .map(move |j| 0.5 * coefficients[i] * expected[3 * i + j] * coefficients[j])
                })
                .sum();
            let work = 0.02 * 0.001 + (-0.03) * (-0.001) + quadratic_work;
            near(event.impulse.external_work_j, work, 2e-11);
            near(event.impulse.momentum.kinetic_change_j, work, 2e-11);
            near(event.impulse.pressure_work_j, 0.0, 2e-11);
            assert_eq!(event.impulse.kkt_inertia, response.kkt_inertia());
            assert_eq!(event.solid.displacement_m(), solid.displacement_m());
            assert_eq!(
                event.solid.algorithmic_pressure_pa(),
                solid.algorithmic_pressure_pa()
            );
            assert_eq!(event.solid.bulk_volume_ratios(), solid.bulk_volume_ratios());
            assert_eq!(event.solid.stored_energy_j(), solid.stored_energy_j());
            assert_eq!(event.solid.time_s(), solid.time_s());
            assert_eq!(state_bits(&solid), before);
            for (c, expected) in [0.001, -0.001, 0.0].into_iter().enumerate() {
                near(event.impulse.external.linear_kg_m_s[c], expected, 2e-12);
            }
            let next = event
                .solid
                .advance(0.0002, &[], &[], [0.0; 3], settings(), WORK_ALLOWANCE, cx)
                .unwrap();
            near(next.energy.cumulative_balance_error_j, 0.0, WORK_ALLOWANCE);
            near(
                next.solid.total_energy_j() - solid.total_energy_j(),
                work,
                WORK_ALLOWANCE,
            );
            assert_ne!(next.solid.displacement_m(), event.solid.displacement_m());
        }
    });
}

#[test]
fn curved_point_impulse_keeps_signed_moving_support_work_and_material_storage() {
    with_cx(false, |cx| {
        let (nodes, elements) = mesh(cx);
        let material = law(false, KAPPA);
        let densities = [DENSITY];
        // Exact P2 shear with det(F)=1. Its curved physical point differs from
        // the affine corner map while the original reference mass stays fixed.
        let displacement: Vec<_> = nodes
            .iter()
            .map(|x| [0.2 * x[1] * x[2], 0.0, 0.0])
            .collect();
        let velocity = vec![[0.0, 0.1, 0.0]; nodes.len()];
        let solid = definition(
            &nodes,
            &elements,
            MixedLawField::Uniform(&material),
            &densities,
        )
        .prepare(&displacement, &velocity, cx)
        .unwrap();
        let point = Tet10ImpulsePoint {
            tetrahedron: 0,
            barycentric: [0.25; 4],
        };
        let before = solid.point_kinematics(point, cx).unwrap();
        for (c, wanted) in [0.2625, 0.25, 0.25].into_iter().enumerate() {
            near(before.position_m[c], wanted, 2e-14);
        }
        let free = 3 * elements[0][4] + 1;
        let held: Vec<_> = (0..3 * nodes.len()).filter(|d| *d != free).collect();
        let response = solid
            .point_impulse_response(&held, None, controls(), cx)
            .unwrap();
        let mobility = response
            .point_directional_inverse_mass_per_kg(point, [0.0, 1.0, 0.0], cx)
            .unwrap();
        near(mobility.inverse_mass_per_kg, 105.0 / 128.0, 2e-12);
        let event = response
            .point_impulses_candidate(&[(point, [0.0, 0.12, 0.0])], [0.0; 3], cx)
            .unwrap();
        let delta = 0.03 / (DENSITY * VOLUME * mass_fraction(4, 4));
        near(
            event
                .solid
                .point_kinematics(point, cx)
                .unwrap()
                .velocity_mps[1],
            0.1 + 0.25 * delta,
            2e-13,
        );
        near(
            event.impulse.external_work_j,
            0.12 * (0.1 + 0.125 * delta),
            2e-13,
        );
        near(event.impulse.momentum.support_work_j, -0.004125, 2e-13);
        for (c, wanted) in [-0.03, 0.0, 0.0315].into_iter().enumerate() {
            near(event.impulse.external.angular_kg_m2_s[c], wanted, 2e-13);
        }
        for dof in held {
            assert_eq!(
                event.solid.velocity_mps()[dof / 3][dof % 3],
                velocity[dof / 3][dof % 3]
            );
        }
        assert_eq!(event.solid.stored_energy_j(), solid.stored_energy_j());
        assert_eq!(event.solid.bulk_volume_ratios(), solid.bulk_volume_ratios());
        assert_eq!(
            event.solid.algorithmic_pressure_pa(),
            solid.algorithmic_pressure_pa()
        );
        assert_eq!(event.solid.displacement_m(), displacement);
    });
}
