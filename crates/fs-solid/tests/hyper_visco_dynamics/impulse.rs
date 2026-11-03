//! Instantaneous momentum transfer retains the original relaxed Maxwell body.
use super::*;
use fs_solid::hyper3::visco::dynamics::{SolidImpulseControls, TetImpulsePoint};
use fs_solver::CheckedCgConfig;

const RHO: [f64; 4] = [30.0, 150.0, 30.0, 30.0];
const EXTENSION: f64 = 1e-4;
const HOLD: f64 = 0.007;

fn impulse_controls() -> SolidImpulseControls {
    SolidImpulseControls {
        solve: CheckedCgConfig {
            tolerance: 1e-13,
            max_iterations: 64,
            max_corrections: 2,
        },
        maximum_workspace_components: 4096,
        linear_tolerance_kg_m_s: 1e-11,
        angular_tolerance_kg_m2_s: 1e-11,
        work_tolerance_j: 1e-11,
    }
}

fn relaxed<'a>(law: &'a ViscoHyperelastic, cx: &Cx<'_>) -> FiniteViscoTetSolid<'a> {
    let mut solid = FiniteViscoTetSolid::try_new(
        &NODES,
        &TETS,
        ViscoelasticLawField::Uniform(law),
        TetAssemblyBudget::default(),
        cx,
    )
    .unwrap();
    let fixed: Vec<_> = (0..12)
        .map(|dof| (dof, if dof == 3 { EXTENSION } else { 0.0 }))
        .collect();
    for dt in [0.0, HOLD] {
        solid = solid
            .advance(
                dt,
                &fixed,
                &[],
                [0.0; 3],
                settings(),
                controls().material,
                cx,
            )
            .unwrap()
            .solid;
    }
    assert!(solid.element_states()[0].dissipated > 0.0);
    solid
}

// Independent simplex integral: rho*N_i*N_j is cubic, so the mean density
// and diagonal/lumped mass cannot provide this impulse oracle.
fn mass(i: usize, j: usize) -> f64 {
    VOLUME * (RHO.iter().sum::<f64>() + RHO[i] + RHO[j]) / if i == j { 60.0 } else { 120.0 }
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn inventory(x: &[[f64; 3]; 4], v: &[[f64; 3]; 4], origin: [f64; 3]) -> ([f64; 3], [f64; 3], f64) {
    let (mut linear, mut angular, mut kinetic) = ([0.0; 3], [0.0; 3], 0.0);
    for i in 0..4 {
        let p = std::array::from_fn(|c| (0..4).map(|j| mass(i, j) * v[j][c]).sum::<f64>());
        let torque = cross(std::array::from_fn(|c| x[i][c] - origin[c]), p);
        for c in 0..3 {
            linear[c] += p[c];
            angular[c] += torque[c];
            kinetic += 0.5 * v[i][c] * p[c];
        }
    }
    (linear, angular, kinetic)
}

#[test]
fn graded_kick_keeps_relaxed_memory_and_accounts_moving_supports_before_native_advance() {
    with_cx(|cx| {
        let law = material(true);
        let old_v = [[0.0, 0.1, 0.0]; 4];
        let motion =
            FiniteViscoTetMotion::try_new_with_nodal_density(relaxed(&law, cx), &RHO, &old_v, cx)
                .unwrap();
        let old_u = motion.solid().displacement_m().to_vec();
        let old_states = motion.solid().element_states().to_vec();
        let old_stress = motion.solid().element_first_piola_pa().to_vec();
        let old_storage = motion.solid().stored_energy_j();
        let x = std::array::from_fn(|i| std::array::from_fn(|c| NODES[i][c] + old_u[i][c]));
        let held: Vec<_> = (0..12).filter(|&dof| dof != 4).collect();
        let origin = [0.4, -0.3, 0.2];
        // Duplicate load entries are one physical impulse, not two events.
        let event = motion
            .impulse_candidate(
                &held,
                &[(4, 0.012), (4, 0.018)],
                origin,
                impulse_controls(),
                cx,
            )
            .unwrap();
        let mut expected_v = old_v;
        expected_v[1][1] += 0.02; // J/M_11 = 0.03/1.5, not 0.03/mean-density mass.
        for i in 0..4 {
            for c in 0..3 {
                close(event.motion.velocity_mps()[i][c], expected_v[i][c], 1e-13);
            }
        }
        let mut support = [0.0; 12];
        for &(dof, value) in &event.impulse.support_impulses_kg_m_s {
            support[dof] += value;
        }
        for i in 0..4 {
            for c in 0..3 {
                let expected = if i != 1 && c == 1 {
                    mass(i, 1) * 0.02
                } else {
                    0.0
                };
                close(support[3 * i + c], expected, 1e-13);
            }
        }
        let before = inventory(&x, &old_v, origin);
        let after = inventory(&x, &expected_v, origin);
        let applied_torque = cross(
            std::array::from_fn(|c| x[1][c] - origin[c]),
            [0.0, 0.03, 0.0],
        );
        for c in 0..3 {
            close(event.impulse.before.linear_kg_m_s[c], before.0[c], 1e-13);
            close(event.impulse.after.linear_kg_m_s[c], after.0[c], 1e-13);
            close(event.impulse.before.angular_kg_m2_s[c], before.1[c], 1e-13);
            close(event.impulse.after.angular_kg_m2_s[c], after.1[c], 1e-13);
            close(
                event.impulse.applied.angular_kg_m2_s[c],
                applied_torque[c],
                1e-13,
            );
            close(event.impulse.linear_residual_kg_m_s[c], 0.0, 1e-13);
            close(event.impulse.angular_residual_kg_m2_s[c], 0.0, 1e-13);
        }
        close(event.impulse.applied_work_j, 0.03 * 0.11, 1e-14);
        close(
            event.impulse.support_work_j,
            0.1 * 3.0 * mass(0, 1) * 0.02,
            1e-14,
        );
        close(event.impulse.kinetic_change_j, after.2 - before.2, 1e-14);
        close(event.impulse.work_residual_j, 0.0, 1e-14);
        assert!(event.impulse.support_work_j > 0.003);

        // One simultaneous resultant must survive cancellation regardless of
        // where its small contribution lies in the supplied contact loads.
        // The oracle remains the exact graded simplex mass above, including
        // its off-diagonal reaction at every moving held node.
        let large = 2.0_f64.powi(54);
        for loads in [
            [(4, large), (4, 0.03), (4, -large)],
            [(4, -large), (4, 0.03), (4, large)],
            [(4, large), (4, -large), (4, 0.03)],
        ] {
            let simultaneous = motion
                .impulse_candidate(&held, &loads, origin, impulse_controls(), cx)
                .unwrap();
            for i in 0..4 {
                for c in 0..3 {
                    close(
                        simultaneous.motion.velocity_mps()[i][c],
                        expected_v[i][c],
                        1e-13,
                    );
                }
            }
            for &(dof, reaction) in &simultaneous.impulse.support_impulses_kg_m_s {
                let expected = if dof % 3 == 1 {
                    mass(dof / 3, 1) * 0.02
                } else {
                    0.0
                };
                close(reaction, expected, 1e-13);
            }
            for c in 0..3 {
                let applied = [0.0, 0.03, 0.0][c];
                close(
                    simultaneous.impulse.applied.linear_kg_m_s[c],
                    applied,
                    1e-13,
                );
                close(
                    simultaneous.impulse.applied.angular_kg_m2_s[c],
                    applied_torque[c],
                    1e-13,
                );
                close(
                    simultaneous.impulse.support.linear_kg_m_s[c],
                    after.0[c] - before.0[c] - applied,
                    1e-13,
                );
                close(
                    simultaneous.impulse.support.angular_kg_m2_s[c],
                    after.1[c] - before.1[c] - applied_torque[c],
                    1e-13,
                );
            }
            close(
                simultaneous.impulse.kinetic_change_j,
                after.2 - before.2,
                1e-14,
            );
            close(simultaneous.impulse.applied_work_j, 0.03 * 0.11, 1e-14);
            close(
                simultaneous.impulse.support_work_j,
                0.1 * 3.0 * mass(0, 1) * 0.02,
                1e-14,
            );
            assert_eq!(simultaneous.motion.solid().element_states(), old_states);
            assert_eq!(simultaneous.motion.solid().displacement_m(), old_u);
            assert_eq!(simultaneous.motion.solid().time_s(), HOLD);
        }
        assert_eq!(event.motion.solid().displacement_m(), old_u);
        assert_eq!(event.motion.solid().element_states(), old_states);
        assert_eq!(event.motion.solid().element_first_piola_pa(), old_stress);
        assert_eq!(event.motion.solid().stored_energy_j(), old_storage);
        assert_eq!(event.motion.solid().time_s(), HOLD);
        assert_eq!(event.motion.reference_mass_kg(), motion.reference_mass_kg());
        assert_eq!(
            event.motion.reference_nodal_density_kg_m3(),
            Some(RHO.as_slice())
        );
        assert_eq!(motion.velocity_mps(), old_v);

        // One genuinely free shear coordinate continues through the native
        // nonlinear dynamic owner; held coordinates translate at their old v.
        let dt = 0.001;
        let fixed: Vec<_> = held
            .iter()
            .map(|&dof| (dof, old_u[dof / 3][dof % 3] + dt * old_v[dof / 3][dof % 3]))
            .collect();
        let next = event
            .motion
            .advance(dt, &fixed, &[], [0.0; 3], settings(), controls(), cx)
            .unwrap();
        let actual = &next.motion.solid().element_states()[0];
        let response = law
            .step_response(&old_states[0], &actual.deformation_gradient, dt)
            .unwrap();
        assert_eq!(*actual, response.state);
        assert_ne!(
            actual.deformation_gradient,
            old_states[0].deformation_gradient
        );
        close(
            next.material_energy.dissipated_j,
            VOLUME * response.dissipated_pa,
            1e-15,
        );
        close(next.motion.solid().time_s(), HOLD + dt, 1e-15);
        close(
            next.balance_error_j,
            0.0,
            controls().material.equilibrium_work_tolerance_j,
        );
        assert_eq!(
            next.motion.reference_nodal_density_kg_m3(),
            Some(RHO.as_slice())
        );
    });
}

#[test]
fn free_rigid_spin_uses_consistent_angular_inertia_and_failed_events_preserve_history() {
    with_cx(|cx| {
        let law = material(true);
        let old_v = [[0.12, -0.03, 0.01]; 4];
        let motion =
            FiniteViscoTetMotion::try_new_with_nodal_density(relaxed(&law, cx), &RHO, &old_v, cx)
                .unwrap();
        let old_u = motion.solid().displacement_m().to_vec();
        let old_states = motion.solid().element_states().to_vec();
        let origin = [-0.2, 0.3, 0.4];
        let x: [[f64; 3]; 4] =
            std::array::from_fn(|i| std::array::from_fn(|c| NODES[i][c] + old_u[i][c]));
        let delta: [[f64; 3]; 4] = std::array::from_fn(|i| {
            let spin = cross(
                [0.03, -0.02, 0.04],
                std::array::from_fn(|c| x[i][c] - origin[c]),
            );
            std::array::from_fn(|c| [0.01, -0.02, 0.03][c] + spin[c])
        });
        let impulses: Vec<_> = (0..12)
            .map(|dof| {
                (
                    dof,
                    (0..4).map(|j| mass(dof / 3, j) * delta[j][dof % 3]).sum(),
                )
            })
            .collect();
        let event = motion
            .impulse_candidate(&[], &impulses, origin, impulse_controls(), cx)
            .unwrap();
        let expected = std::array::from_fn(|i| std::array::from_fn(|c| old_v[i][c] + delta[i][c]));
        let expected_inventory = inventory(&x, &expected, origin);
        for i in 0..4 {
            for c in 0..3 {
                close(event.motion.velocity_mps()[i][c], expected[i][c], 1e-12);
            }
        }
        for c in 0..3 {
            close(
                event.impulse.after.linear_kg_m_s[c],
                expected_inventory.0[c],
                1e-12,
            );
            close(
                event.impulse.after.angular_kg_m2_s[c],
                expected_inventory.1[c],
                1e-12,
            );
            close(event.impulse.support.linear_kg_m_s[c], 0.0, 0.0);
        }
        close(event.impulse.kinetic_after_j, expected_inventory.2, 1e-13);
        assert!(event.impulse.support_impulses_kg_m_s.is_empty());
        assert_eq!(event.motion.solid().element_states(), old_states);
        let mut limited = impulse_controls();
        limited.maximum_workspace_components = 0;
        assert!(
            motion
                .impulse_candidate(&[], &impulses, origin, limited, cx)
                .is_err()
        );
        assert!(
            motion
                .impulse_candidate(&[1, 1], &impulses, origin, impulse_controls(), cx)
                .is_err()
        );
        let cancel = CancelGate::new_clock_free();
        cancel.request();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let stopped = Cx::new(
                &cancel,
                arena,
                StreamKey {
                    seed: 42,
                    kernel_id: 17,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            assert!(
                motion
                    .impulse_candidate(&[], &impulses, origin, impulse_controls(), &stopped)
                    .is_err()
            );
        });
        assert_eq!(motion.solid().displacement_m(), old_u);
        assert_eq!(motion.solid().element_states(), old_states);
        assert_eq!(motion.solid().time_s(), HOLD);
        assert_eq!(motion.velocity_mps(), old_v);
        let retry = motion
            .impulse_candidate(&[], &impulses, origin, impulse_controls(), cx)
            .unwrap();
        assert_eq!(retry.motion.velocity_mps(), event.motion.velocity_mps());
        assert_eq!(
            retry.motion.solid().element_states(),
            event.motion.solid().element_states()
        );

        // Equal/opposite collinear impulses have zero net force, torque and
        // rigid-translation work. At this baseline velocity both local kicks
        // are below one ULP, so only the represented free-DOF balance detects
        // their loss. A global-only momentum gate would incorrectly accept.
        let fast_velocity = [[1e16, 0.0, 0.0]; 4];
        let fast = FiniteViscoTetMotion::try_new_with_nodal_density(
            relaxed(&law, cx),
            &RHO,
            &fast_velocity,
            cx,
        )
        .unwrap();
        assert!(
            fast.impulse_candidate(
                &[],
                &[(0, 0.01), (3, -0.01)],
                [0.0; 3],
                impulse_controls(),
                cx,
            )
            .is_err()
        );
        assert_eq!(fast.velocity_mps(), fast_velocity);
    });
}

#[test]
fn interior_point_trace_and_supported_mobility_share_the_exact_graded_mass() {
    with_cx(|cx| {
        let law = material(true);
        let old_v: [[f64; 3]; 4] = std::array::from_fn(|i| {
            let i = i as f64;
            [0.02 + 0.03 * i, -0.01 + 0.02 * i, 0.04 - 0.01 * i]
        });
        let motion =
            FiniteViscoTetMotion::try_new_with_nodal_density(relaxed(&law, cx), &RHO, &old_v, cx)
                .unwrap();
        let point = TetImpulsePoint {
            tetrahedron: 0,
            weights: [0.1, 0.2, 0.3, 0.4],
        };
        let direction: [f64; 3] = [0.6, 0.8, 0.0];
        let impulse = 0.03;
        let applied = direction.map(|n| impulse * n);
        let origin = [-0.2, 0.3, 0.4];
        let x: [[f64; 3]; 4] = std::array::from_fn(|i| {
            std::array::from_fn(|c| NODES[i][c] + motion.solid().displacement_m()[i][c])
        });
        let before = motion.point_kinematics(point, cx).unwrap();
        let position = [0.2 * (1.0 + EXTENSION), 0.3, 0.4];
        let velocity = [0.08, 0.03, 0.02];
        for c in 0..3 {
            close(before.position_m[c], position[c], 1e-15);
            close(before.velocity_mps[c], velocity[c], 1e-15);
        }
        // Each Cartesian mass block is A/12, with diagonal [10,18,10,10],
        // off-diagonal 7 at the dense node and 5 between the other nodes.
        // Exact inversion gives M^-1*w = [-370,-160,482,908]/1775.
        // In particular, a point impulse produces negative nodal responses:
        // neither lumping nor substituting the cell mean can pass this oracle.
        let inverse_weights = [-370.0, -160.0, 482.0, 908.0].map(|x| x / 1775.0);
        let free_mobility = 2194.0 / 8875.0;
        let supported: Vec<_> = (0..12).filter(|&dof| dof != 4).collect();
        for held in [&[][..], supported.as_slice()] {
            let restrained = !held.is_empty();
            let expected_mobility = if restrained {
                (0.2 * direction[1]).powi(2) / 1.5
            } else {
                free_mobility
            };
            let mobility = motion
                .point_directional_inverse_mass_per_kg(
                    held,
                    point,
                    direction,
                    impulse_controls(),
                    cx,
                )
                .unwrap();
            assert!(mobility.solve.converged_euclidean());
            close(mobility.inverse_mass_per_kg, expected_mobility, 1e-13);
            let event = motion
                .point_impulse_candidate(held, point, applied, origin, impulse_controls(), cx)
                .unwrap();
            let expected: [[f64; 3]; 4] = std::array::from_fn(|i| {
                std::array::from_fn(|c| {
                    old_v[i][c]
                        + if restrained {
                            if i == 1 && c == 1 {
                                impulse * 0.2 * direction[1] / 1.5
                            } else {
                                0.0
                            }
                        } else {
                            inverse_weights[i] * applied[c]
                        }
                })
            });
            for i in 0..4 {
                for c in 0..3 {
                    close(event.motion.velocity_mps()[i][c], expected[i][c], 1e-13);
                }
            }
            let after = event.motion.point_kinematics(point, cx).unwrap();
            assert_eq!(after.position_m, before.position_m);
            let directional_change: f64 = (0..3)
                .map(|c| direction[c] * (after.velocity_mps[c] - before.velocity_mps[c]))
                .sum();
            close(directional_change, impulse * expected_mobility, 1e-14);
            let point_work: f64 = (0..3)
                .map(|c| {
                    applied[c]
                        * (before.velocity_mps[c]
                            + 0.5 * (after.velocity_mps[c] - before.velocity_mps[c]))
                })
                .sum();
            close(event.impulse.applied_work_j, point_work, 1e-14);
            let torque = cross(std::array::from_fn(|c| position[c] - origin[c]), applied);
            for c in 0..3 {
                close(event.impulse.applied.linear_kg_m_s[c], applied[c], 1e-14);
                close(event.impulse.applied.angular_kg_m2_s[c], torque[c], 1e-14);
            }
            let initial_inventory = inventory(&x, &old_v, origin);
            let final_inventory = inventory(&x, &expected, origin);
            close(
                event.impulse.kinetic_change_j,
                final_inventory.2 - initial_inventory.2,
                1e-14,
            );
            let expected_support_work: f64 = if restrained {
                let delta = impulse * 0.2 * direction[1] / 1.5;
                held.iter()
                    .map(|&dof| {
                        let i = dof / 3;
                        let c = dof % 3;
                        let dp = if c == 1 { mass(i, 1) * delta } else { 0.0 };
                        old_v[i][c] * (dp - point.weights[i] * applied[c])
                    })
                    .sum()
            } else {
                0.0
            };
            close(event.impulse.support_work_j, expected_support_work, 1e-14);
            close(
                event.impulse.kinetic_change_j,
                point_work + expected_support_work,
                1e-14,
            );
            assert_eq!(
                event.motion.solid().element_states(),
                motion.solid().element_states()
            );
            assert_eq!(event.motion.solid().time_s(), HOLD);
        }
        // A fully restrained interior trace has zero mobility without changing
        // its nonzero accepted material-point velocity.
        let all: Vec<_> = (0..12).collect();
        close(
            motion
                .point_directional_inverse_mass_per_kg(
                    &all,
                    point,
                    direction,
                    impulse_controls(),
                    cx,
                )
                .unwrap()
                .inverse_mass_per_kg,
            0.0,
            0.0,
        );
        for bad in [
            TetImpulsePoint {
                tetrahedron: 1,
                ..point
            },
            TetImpulsePoint {
                weights: [0.1, 0.2, 0.3, 0.3],
                ..point
            },
            TetImpulsePoint {
                weights: [-0.1, 0.2, 0.3, 0.6],
                ..point
            },
            TetImpulsePoint {
                weights: [f64::NAN, 0.2, 0.3, 0.5],
                ..point
            },
        ] {
            assert!(motion.point_kinematics(bad, cx).is_err());
            assert!(
                motion
                    .point_directional_inverse_mass_per_kg(
                        &[],
                        bad,
                        direction,
                        impulse_controls(),
                        cx,
                    )
                    .is_err()
            );
            assert!(
                motion
                    .point_impulse_candidate(&[], bad, applied, origin, impulse_controls(), cx,)
                    .is_err()
            );
        }
        for bad_direction in [[0.0; 3], [2.0, 0.0, 0.0]] {
            assert!(
                motion
                    .point_directional_inverse_mass_per_kg(
                        &[],
                        point,
                        bad_direction,
                        impulse_controls(),
                        cx,
                    )
                    .is_err()
            );
        }
        assert_eq!(motion.velocity_mps(), old_v);
        assert_eq!(motion.point_kinematics(point, cx).unwrap(), before);
    });
}
