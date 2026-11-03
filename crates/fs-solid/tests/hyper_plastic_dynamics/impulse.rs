//! Graded-mass impulse events on an already plastically set finite solid.
use super::*;
use fs_solid::hyper3::plastic::dynamics::SolidImpulseControls;
use fs_solver::CheckedCgConfig;

const RHO: [f64; 4] = [3.0, 15.0, 3.0, 3.0];
const MASS_SHARES: [f64; 4] = [0.225, 0.325, 0.225, 0.225];

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

fn plastically_set<'a>(material: &'a FiniteJ2, cx: &Cx<'_>) -> FinitePlasticTetSolid<'a> {
    let mut solid = FinitePlasticTetSolid::try_new(
        &NODES,
        &TETS,
        FinitePlasticLawField::Uniform(material),
        TetAssemblyBudget::standard(),
        cx,
    )
    .unwrap();
    // Commuting isochoric return, followed by unloading to its stress-free
    // plastic shape. The retained alpha/hardening energy must survive impact.
    let a = 0.20;
    let plastic = (3.0 * MU * a - YIELD) / (3.0 * MU + HARDENING);
    for stretch in [a, plastic] {
        let diagonal = [
            stretch.exp(),
            (-0.5 * stretch).exp(),
            (-0.5 * stretch).exp(),
        ];
        let fixed: Vec<_> = (0..12)
            .map(|dof| (dof, (diagonal[dof % 3] - 1.0) * NODES[dof / 3][dof % 3]))
            .collect();
        solid = solid
            .advance(&fixed, &[], [0.0; 3], settings(), controls().material, cx)
            .unwrap()
            .solid;
    }
    close(solid.element_states()[0].alpha, plastic, 1e-12);
    assert!(solid.element_states()[0].alpha > 0.18);
    assert!(solid.stored_energy_j() > 2.0);
    solid
}

// Exact integral of the original P1 density times two P1 shape functions.
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

#[test]
fn free_impulse_translates_the_plastically_set_graded_body_through_native_motion() {
    with_cx(|cx| {
        let material = law();
        let v0 = [0.2, -0.1, 0.05];
        let delta = [0.03, -0.02, 0.04];
        let expected_v: [f64; 3] = std::array::from_fn(|c| v0[c] + delta[c]);
        let motion = FinitePlasticTetMotion::try_new_with_nodal_density(
            plastically_set(&material, cx),
            &RHO,
            &[v0; 4],
            cx,
        )
        .unwrap();
        let old_u = motion.solid().displacement_m().to_vec();
        let old_states = motion.solid().element_states().to_vec();
        let old_f = motion.solid().element_deformation().to_vec();
        let old_stress = motion.solid().element_first_piola_pa().to_vec();
        let storage = motion.solid().stored_energy_j();
        let origin = [-0.2, 0.3, 0.4];
        let first_moment: [f64; 3] = std::array::from_fn(|c| {
            (0..4)
                .map(|i| MASS_SHARES[i] * (NODES[i][c] + old_u[i][c] - origin[c]))
                .sum()
        });
        let impulses: Vec<_> = (0..12)
            .map(|dof| (dof, MASS_SHARES[dof / 3] * delta[dof % 3]))
            .collect();
        let event = motion
            .impulse_candidate(&[], &impulses, origin, impulse_controls(), cx)
            .unwrap();
        for velocity in event.motion.velocity_mps() {
            for c in 0..3 {
                close(velocity[c], expected_v[c], 1e-13);
            }
        }
        let before_angular = cross(first_moment, v0);
        let after_angular = cross(first_moment, expected_v);
        let applied_angular = cross(first_moment, delta);
        for c in 0..3 {
            close(event.impulse.before.linear_kg_m_s[c], v0[c], 1e-13);
            close(event.impulse.after.linear_kg_m_s[c], expected_v[c], 1e-13);
            close(event.impulse.applied.linear_kg_m_s[c], delta[c], 1e-13);
            close(
                event.impulse.before.angular_kg_m2_s[c],
                before_angular[c],
                1e-13,
            );
            close(
                event.impulse.after.angular_kg_m2_s[c],
                after_angular[c],
                1e-13,
            );
            close(
                event.impulse.applied.angular_kg_m2_s[c],
                applied_angular[c],
                1e-13,
            );
            close(event.impulse.linear_residual_kg_m_s[c], 0.0, 1e-13);
            close(event.impulse.angular_residual_kg_m2_s[c], 0.0, 1e-13);
        }
        close(
            event.impulse.kinetic_change_j,
            0.5 * (norm_squared(expected_v) - norm_squared(v0)),
            1e-14,
        );
        close(
            event.impulse.kinetic_change_j,
            event.impulse.applied_work_j,
            1e-14,
        );
        close(event.impulse.support_work_j, 0.0, 0.0);
        assert!(event.impulse.support_impulses_kg_m_s.is_empty());
        assert_eq!(event.motion.solid().displacement_m(), old_u);
        assert_eq!(event.motion.solid().element_deformation(), old_f);
        assert_eq!(event.motion.solid().element_states(), old_states);
        assert_eq!(event.motion.solid().element_first_piola_pa(), old_stress);
        assert_eq!(event.motion.solid().stored_energy_j(), storage);
        assert_eq!(event.motion.reference_mass_kg(), motion.reference_mass_kg());
        assert_eq!(
            event.motion.reference_nodal_density_kg_m3(),
            Some(RHO.as_slice())
        );
        assert_eq!(motion.velocity_mps(), [v0; 4]);

        // The accepted impulse becomes the original velocity of the next
        // native interval. Translation must not replay/reset the plastic return.
        let dt = 0.02;
        let next = event
            .motion
            .advance(&[], &[], [0.0; 3], dt, settings(), controls(), cx)
            .unwrap();
        assert!(next.equilibrium.reactions_n.is_empty());
        assert_eq!(next.motion.solid().element_states(), old_states);
        for i in 0..4 {
            for c in 0..3 {
                close(
                    next.motion.solid().displacement_m()[i][c],
                    old_u[i][c] + dt * expected_v[c],
                    1e-10,
                );
                close(next.motion.velocity_mps()[i][c], expected_v[c], 1e-8);
            }
        }
        close(next.material_energy.stored_before_j, storage, 1e-9);
        close(next.material_energy.stored_after_j, storage, 1e-9);
        close(next.material_energy.plastic_dissipation_j, 0.0, 0.0);
        close(next.inertia.kinetic_change_j, 0.0, 1e-9);
        close(next.balance_error_j, 0.0, 1e-9);
        assert_eq!(
            next.motion.reference_nodal_density_kg_m3(),
            Some(RHO.as_slice())
        );
    });
}

#[test]
fn deformed_plastic_body_reports_moving_support_work_and_preserves_state_on_refusal() {
    with_cx(|cx| {
        let material = law();
        let old_v = [[0.0, 0.1, 0.0]; 4];
        let motion = FinitePlasticTetMotion::try_new_with_nodal_density(
            plastically_set(&material, cx),
            &RHO,
            &old_v,
            cx,
        )
        .unwrap();
        let old_u = motion.solid().displacement_m().to_vec();
        let old_states = motion.solid().element_states().to_vec();
        let old_f = motion.solid().element_deformation().to_vec();
        let old_stress = motion.solid().element_first_piola_pa().to_vec();
        let storage = motion.solid().stored_energy_j();
        let held: Vec<_> = (0..12).filter(|&dof| dof != 4).collect();
        let origin = [0.4, -0.3, 0.2];
        let impulses = [(4, 0.03)];
        let event = motion
            .impulse_candidate(&held, &impulses, origin, impulse_controls(), cx)
            .unwrap();
        let delta_v = 0.2; // 0.03 / exact graded M_11=0.15 kg.
        let mut support = [0.0; 12];
        for &(dof, value) in &event.impulse.support_impulses_kg_m_s {
            support[dof] += value;
        }
        let mut total_torque = [0.0; 3];
        for i in 0..4 {
            let impulse = [0.0, mass(i, 1) * delta_v, 0.0];
            let torque = cross(
                std::array::from_fn(|c| NODES[i][c] + old_u[i][c] - origin[c]),
                impulse,
            );
            for c in 0..3 {
                let expected_v = old_v[i][c] + if i == 1 && c == 1 { delta_v } else { 0.0 };
                close(event.motion.velocity_mps()[i][c], expected_v, 1e-13);
                close(
                    support[3 * i + c],
                    if i == 1 { 0.0 } else { impulse[c] },
                    1e-13,
                );
                total_torque[c] += torque[c];
            }
        }
        let applied_torque = cross(
            std::array::from_fn(|c| NODES[1][c] + old_u[1][c] - origin[c]),
            [0.0, 0.03, 0.0],
        );
        for c in 0..3 {
            close(
                event.impulse.applied.angular_kg_m2_s[c],
                applied_torque[c],
                1e-13,
            );
            close(
                event.impulse.support.angular_kg_m2_s[c],
                total_torque[c] - applied_torque[c],
                1e-13,
            );
            close(event.impulse.angular_residual_kg_m2_s[c], 0.0, 1e-13);
        }
        close(
            event.impulse.applied_work_j,
            0.03 * (0.1 + 0.5 * delta_v),
            1e-14,
        );
        close(
            event.impulse.support_work_j,
            0.1 * 3.0 * mass(0, 1) * delta_v,
            1e-14,
        );
        close(
            event.impulse.kinetic_change_j,
            0.1 * MASS_SHARES[1] * delta_v + 0.5 * mass(1, 1) * delta_v * delta_v,
            1e-14,
        );
        close(event.impulse.work_residual_j, 0.0, 1e-14);
        assert!(event.impulse.support_work_j > 0.003);

        // The original graded mass must receive the net simultaneous load,
        // even when large opposite impulses surround a much smaller one.
        // Neither plastic history nor the moving supports may absorb its loss.
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
                    let expected_v = old_v[i][c] + if i == 1 && c == 1 { delta_v } else { 0.0 };
                    close(simultaneous.motion.velocity_mps()[i][c], expected_v, 1e-13);
                }
            }
            for &(dof, reaction) in &simultaneous.impulse.support_impulses_kg_m_s {
                let expected = if dof % 3 == 1 {
                    mass(dof / 3, 1) * delta_v
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
                    if c == 1 {
                        3.0 * mass(0, 1) * delta_v
                    } else {
                        0.0
                    },
                    1e-13,
                );
                close(
                    simultaneous.impulse.support.angular_kg_m2_s[c],
                    total_torque[c] - applied_torque[c],
                    1e-13,
                );
            }
            close(
                simultaneous.impulse.kinetic_change_j,
                0.1 * MASS_SHARES[1] * delta_v + 0.5 * mass(1, 1) * delta_v * delta_v,
                1e-14,
            );
            close(
                simultaneous.impulse.applied_work_j,
                0.03 * (0.1 + 0.5 * delta_v),
                1e-14,
            );
            close(
                simultaneous.impulse.support_work_j,
                0.1 * 3.0 * mass(0, 1) * delta_v,
                1e-14,
            );
            assert_eq!(simultaneous.motion.solid().element_states(), old_states);
            assert_eq!(simultaneous.motion.solid().element_deformation(), old_f);
            assert_eq!(simultaneous.motion.solid().displacement_m(), old_u);
        }
        assert_eq!(event.motion.solid().displacement_m(), old_u);
        assert_eq!(event.motion.solid().element_deformation(), old_f);
        assert_eq!(event.motion.solid().element_states(), old_states);
        assert_eq!(event.motion.solid().element_first_piola_pa(), old_stress);
        assert_eq!(event.motion.solid().stored_energy_j(), storage);
        assert_eq!(
            event.motion.reference_nodal_density_kg_m3(),
            Some(RHO.as_slice())
        );
        let mut limited = impulse_controls();
        limited.maximum_workspace_components = 0;
        assert!(
            motion
                .impulse_candidate(&held, &impulses, origin, limited, cx)
                .is_err()
        );
        assert!(
            motion
                .impulse_candidate(&held, &[(4, f64::NAN)], origin, impulse_controls(), cx)
                .is_err()
        );
        let cancel = CancelGate::new_clock_free();
        cancel.request();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let stopped = Cx::new(
                &cancel,
                arena,
                StreamKey {
                    seed: 3,
                    kernel_id: 17,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            assert!(
                motion
                    .impulse_candidate(&held, &impulses, origin, impulse_controls(), &stopped)
                    .is_err()
            );
        });
        assert_eq!(motion.velocity_mps(), old_v);
        assert_eq!(motion.solid().element_states(), old_states);
        assert_eq!(motion.solid().displacement_m(), old_u);
        assert_eq!(motion.solid().element_deformation(), old_f);
        assert_eq!(motion.solid().element_first_piola_pa(), old_stress);
        assert_eq!(motion.solid().stored_energy_j(), storage);
        let retry = motion
            .impulse_candidate(&held, &impulses, origin, impulse_controls(), cx)
            .unwrap();
        assert_eq!(retry.motion.velocity_mps(), event.motion.velocity_mps());
        assert_eq!(
            retry.motion.solid().element_states(),
            event.motion.solid().element_states()
        );
    });
}
