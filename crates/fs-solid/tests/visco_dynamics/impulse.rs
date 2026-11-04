use super::*;

fn full_system(law: &TetViscoMaterial, free: &[usize], cx: &Cx<'_>) -> ViscoElements {
    let fixed: Vec<_> = (0..12).filter(|d| !free.contains(d)).collect();
    ViscoElements::new_with_full_reference_mass(
        SolidMeshRef::Tet4 {
            nodes_m: &NODES,
            elements: &CELLS,
        },
        TetViscoField::Uniform(law),
        &fixed,
        TetAssemblyBudget::standard(),
        ViscoStateBudget {
            max_internal_variables: 1000,
        },
        cx,
    )
    .unwrap()
}

#[test]
fn native_impulse_preserves_prony_history_and_exact_free_body_momenta() {
    with_cx(|cx, _| {
        let law = coupled();
        let free: Vec<_> = (0..12).collect();
        let solid = system(&law, &free, cx);
        let prepared = solid.dynamic_operator(0.002, config(), cx).unwrap();
        let mut history = solid.state(cx).unwrap();
        let mut preload = [0.0; 12];
        preload[3] = 1e-4;
        preload[7] = -2e-5;
        history.step_internal_force(&preload, 0.01, cx).unwrap();
        let initial_v: Vec<_> = (0..4).flat_map(|_| [0.01, -0.02, 0.03]).collect();
        let mut state = prepared.start(history, &initial_v, cx).unwrap();
        state.advance(&[0.0; 12], cx).unwrap();
        let before = state.clone();
        assert!(before.material_history().energy().dissipated_j > 0.0);
        assert_eq!(state.reference_nodes_m(), &NODES);
        assert_eq!(solid.reference_nodes_m(), &NODES);
        let j = [0.4, -0.25, 0.15];
        let mut applied = [0.0; 12];
        applied[3..6].copy_from_slice(&j);
        let response = state.momentum_response(&applied, cx).unwrap();
        let mass = law.density_kg_m3 / 6.0;
        // Exact P1 inverse: M=(m/20)*(I+11^T), hence a node-1
        // impulse gives 16 J/m there and -4 J/m at each other node.
        for node in 0..4 {
            for axis in 0..3 {
                let expected = if node == 1 { 16.0 } else { -4.0 } * j[axis] / mass;
                close(response[3 * node + axis], expected, 3e-12);
            }
        }
        assert_eq!(state.velocity_mps(), before.velocity_mps());
        let report = state.apply_momentum_impulse(&applied, 2e-12, cx).unwrap();
        for d in 0..12 {
            close(
                state.velocity_mps()[d] - before.velocity_mps()[d],
                response[d],
                3e-12,
            );
        }
        for origin in [[0.0; 3], [0.7, -0.2, 0.4]] {
            let mut linear = [0.0; 3];
            let mut angular = [0.0; 3];
            for (node, x) in NODES.iter().enumerate() {
                let arm: [f64; 3] = core::array::from_fn(|a| {
                    x[a] + before.material_history().displacement_m()[3 * node + a] - origin[a]
                });
                let p: [f64; 3] =
                    core::array::from_fn(|a| report.momentum_change_kg_m_s[3 * node + a]);
                for a in 0..3 {
                    linear[a] += p[a];
                    let b = (a + 1) % 3;
                    let c = (a + 2) % 3;
                    angular[a] += arm[b] * p[c] - arm[c] * p[b];
                }
            }
            let arm: [f64; 3] = core::array::from_fn(|a| {
                NODES[1][a] + before.material_history().displacement_m()[3 + a] - origin[a]
            });
            for a in 0..3 {
                let b = (a + 1) % 3;
                let c = (a + 2) % 3;
                assert!((linear[a] - j[a]).abs() < 2e-12);
                assert!((angular[a] - (arm[b] * j[c] - arm[c] * j[b])).abs() < 2e-12);
            }
        }
        assert!(report.support_impulse_kg_m_s.iter().all(|j| *j == 0.0));
        assert_eq!(report.support_work_j, 0.0);
        close(report.external_work_j, report.kinetic_change_j, 3e-12);
        close(
            state.external_work_j() - before.external_work_j(),
            report.external_work_j,
            3e-12,
        );
        assert_eq!(state.time_s(), before.time_s());
        assert_eq!(state.steps(), before.steps());
        assert_eq!(state.eigenstrain_work_j(), before.eigenstrain_work_j());
        assert_eq!(state.parameter_work_j(), before.parameter_work_j());
        assert_eq!(
            state.material_history().energy(),
            before.material_history().energy()
        );
        assert_eq!(
            state.material_history().displacement_m(),
            before.material_history().displacement_m()
        );
        // Hidden Maxwell memories must produce the identical subsequent
        // constitutive response under the same prescribed strain ramp.
        let mut old_memory = before.material_history().clone();
        let mut new_memory = state.material_history().clone();
        let a = old_memory.step_internal_force(&preload, 0.003, cx).unwrap();
        let b = new_memory.step_internal_force(&preload, 0.003, cx).unwrap();
        assert_eq!(a.full_force_n, b.full_force_n);
        assert_eq!(a.full_mean_force_n, b.full_mean_force_n);
        assert_eq!(a.energy, b.energy);
        let step = state.advance(&[0.0; 12], cx).unwrap();
        assert_eq!(state.steps(), before.steps() + 1);
        close(step.time_s, before.time_s() + 0.002, 1e-15);
        assert!(step.constitutive.step_dissipated_j >= 0.0);
    });
}

#[test]
fn moving_support_impulses_retain_cross_mass_reactions_and_signed_work() {
    with_cx(|cx, _| {
        let law = coupled();
        let solid = full_system(&law, &[3], cx);
        let prepared = solid.dynamic_operator(0.001, config(), cx).unwrap();
        let mut v = [0.0; 12];
        v[0] = 0.7;
        v[3] = -0.02;
        v[6] = -0.1;
        let mut state = prepared
            .start_with_support_motion(solid.state(cx).unwrap(), &v, cx)
            .unwrap();
        let mut j = [0.0; 12];
        j[3] = 2.0;
        j[0] = -0.3;
        let report = state.apply_momentum_impulse(&j, 1e-12, cx).unwrap();
        let delta = 2.0 / (law.density_kg_m3 / 60.0);
        close(state.velocity_mps()[3], v[3] + delta, 1e-14);
        for d in (0..12).filter(|d| *d != 3) {
            assert_eq!(state.velocity_mps()[d], v[d]);
        }
        // P1 cross mass is half its diagonal: each supported x row
        // receives one kg m/s when the free row receives two.
        close(report.support_impulse_kg_m_s[0], 1.3, 1e-14);
        close(report.support_impulse_kg_m_s[6], 1.0, 1e-14);
        close(report.support_impulse_kg_m_s[9], 1.0, 1e-14);
        assert_eq!(report.support_impulse_kg_m_s[3], 0.0);
        close(report.external_work_j, -0.13, 1e-13);
        close(report.support_work_j, 0.81, 1e-13);
        close(report.kinetic_change_j, 0.68, 1e-13);
        close(state.support_work_j(), 0.81, 1e-13);
        close(state.external_work_j(), -0.13, 1e-13);
        assert_eq!(state.time_s(), 0.0);
        assert_eq!(state.steps(), 0);
        let mut support = v;
        support[3] = 0.0;
        state
            .advance_with_support_velocity(&[0.0; 12], &support, cx)
            .unwrap();

        // A fully supported body has no mass solve: the actuator opposes
        // every applied impulse and its work cancels applied impulse work.
        let held = full_system(&law, &[], cx);
        let operator = held.dynamic_operator(0.001, config(), cx).unwrap();
        let mut fixed = operator
            .start_with_support_motion(held.state(cx).unwrap(), &v, cx)
            .unwrap();
        let report = fixed.apply_momentum_impulse(&j, 0.0, cx).unwrap();
        assert_eq!(report.iterations, 0);
        assert_eq!(fixed.velocity_mps(), &v);
        assert_eq!(report.kinetic_change_j, 0.0);
        assert_eq!(report.external_work_j, -report.support_work_j);
        for d in 0..12 {
            assert_eq!(report.support_impulse_kg_m_s[d], -j[d]);
        }
    });
}

#[test]
fn impulse_refusals_preserve_all_state_and_reject_unrepresented_velocity_changes() {
    with_cx(|cx, gate| {
        let law = coupled();
        let solid = system(&law, &(0..12).collect::<Vec<_>>(), cx);
        let prepared = solid.dynamic_operator(0.002, config(), cx).unwrap();
        let mut state = prepared
            .start(solid.state(cx).unwrap(), &[1e16; 12], cx)
            .unwrap();
        let before = state.clone();
        let mut impulse = [0.0; 12];
        impulse[3] = 0.01;
        assert!(matches!(
            state.apply_momentum_impulse(&impulse, 1e-12, cx),
            Err(ViscoElementError::ImpulseMomentum { .. })
        ));
        assert!(state.apply_momentum_impulse(&[0.0; 11], 1e-12, cx).is_err());
        assert!(
            state
                .apply_momentum_impulse(&[f64::NAN; 12], 1e-12, cx)
                .is_err()
        );
        assert!(state.apply_momentum_impulse(&[0.0; 12], -1.0, cx).is_err());
        let legacy = system(&law, &[3], cx);
        let operator = legacy.dynamic_operator(0.002, config(), cx).unwrap();
        let mut reduced_only = operator
            .start(legacy.state(cx).unwrap(), &[0.0; 12], cx)
            .unwrap();
        assert!(
            reduced_only
                .apply_momentum_impulse(&impulse, 1e-12, cx)
                .is_err()
        );
        let supported = full_system(&law, &[3], cx);
        let operator = supported.dynamic_operator(0.002, config(), cx).unwrap();
        let mut supported_state = operator
            .start(supported.state(cx).unwrap(), &[0.0; 12], cx)
            .unwrap();
        let mut opposing = [0.0; 12];
        opposing[0] = 1e20;
        opposing[3] = 2.0;
        // The exact reaction at node 0 is 1-1e20, which f64 cannot
        // represent to this tolerance. Its missing unit must be refused.
        assert!(matches!(
            supported_state.apply_momentum_impulse(&opposing, 1e-12, cx),
            Err(ViscoElementError::ImpulseMomentum { .. })
        ));
        assert_eq!(supported_state.velocity_mps(), &[0.0; 12]);
        gate.request();
        assert!(state.apply_momentum_impulse(&[0.0; 12], 1e-12, cx).is_err());
        assert_eq!(state.velocity_mps(), before.velocity_mps());
        assert_eq!(state.kinetic_energy_j(), before.kinetic_energy_j());
        assert_eq!(
            state.material_history().energy(),
            before.material_history().energy()
        );
        assert_eq!(
            state.material_history().displacement_m(),
            before.material_history().displacement_m()
        );
        assert_eq!(state.external_work_j(), before.external_work_j());
        assert_eq!(state.support_work_j(), before.support_work_j());
        assert_eq!(state.time_s(), before.time_s());
        assert_eq!(state.steps(), before.steps());
    });
}

#[test]
fn opposing_impulses_retain_kinetic_work_below_the_velocity_midpoint_ulp() {
    with_cx(|cx, _| {
        let mut law = coupled();
        law.density_kg_m3 = 120.0;
        let solid = system(&law, &(0..12).collect::<Vec<_>>(), cx);
        let prepared = solid.dynamic_operator(0.002, config(), cx).unwrap();
        let speed = 2.0_f64.powi(53);
        let mut state = prepared
            .start(solid.state(cx).unwrap(), &[speed; 12], cx)
            .unwrap();
        let mut impulse = [0.0; 12];
        impulse[0] = 2.0;
        impulse[3] = -2.0;
        // rho=120 on the unit tetrahedron gives each scalar mass block
        // I + 11^T. These opposing impulses therefore give exact changes
        // (+2,-2,0,0), with zero net linear impulse and Delta K=4 J.
        // The first velocity midpoint is 2^53+1, unrepresentable as f64.
        let report = state.apply_momentum_impulse(&impulse, 1e-12, cx).unwrap();
        assert_eq!(state.velocity_mps()[0], speed + 2.0);
        assert_eq!(state.velocity_mps()[3], speed - 2.0);
        close(report.kinetic_change_j, 4.0, 1e-14);
        close(report.external_work_j, 4.0, 1e-14);
        close(state.external_work_j(), 4.0, 1e-14);
        assert_eq!(report.support_work_j, 0.0);
        assert_eq!(report.balance_residual_j, 0.0);
        assert_eq!(state.time_s(), 0.0);
        assert_eq!(state.steps(), 0);
        assert_eq!(state.material_history().energy().dissipated_j, 0.0);
    });
}

#[test]
fn prepared_impulse_discards_or_publishes_without_restarting_the_original_trajectory() {
    with_cx(|cx, _| {
        let law = coupled();
        let solid = system(&law, &(0..12).collect::<Vec<_>>(), cx);
        let prepared = solid.dynamic_operator(0.002, config(), cx).unwrap();
        let mut history = solid.state(cx).unwrap();
        let mut preload = [0.0; 12];
        preload[3] = 1e-4;
        history.step_internal_force(&preload, 0.01, cx).unwrap();
        let mut state = prepared.start(history, &[0.0; 12], cx).unwrap();
        state.advance(&[0.0; 12], cx).unwrap();
        let before_velocity = state.velocity_mps().to_vec();
        let before_displacement = state.material_history().displacement_m().to_vec();
        let before_energy = state.material_history().energy();
        let before_time = state.time_s();
        let before_steps = state.steps();
        let before_kinetic = state.kinetic_energy_j();
        let before_work = state.external_work_j();
        let mut impulse = [0.0; 12];
        impulse[3] = 0.25;
        {
            let proposal = state.prepare_momentum_impulse(&impulse, 1e-12, cx).unwrap();
            assert_eq!(
                proposal.original().velocity_mps(),
                before_velocity.as_slice()
            );
            assert_eq!(
                proposal.original().material_history().energy(),
                before_energy
            );
            assert_eq!(proposal.original().time_s(), before_time);
            assert_eq!(proposal.original().steps(), before_steps);
            assert_eq!(proposal.report().kinetic_before_j, before_kinetic);
            assert_ne!(proposal.velocity_mps(), before_velocity.as_slice());
        }
        assert_eq!(state.velocity_mps(), before_velocity.as_slice());
        assert_eq!(state.kinetic_energy_j(), before_kinetic);
        assert_eq!(state.external_work_j(), before_work);
        assert_eq!(
            state.material_history().displacement_m(),
            before_displacement.as_slice()
        );
        assert_eq!(state.material_history().energy(), before_energy);
        let proposal = state.prepare_momentum_impulse(&impulse, 1e-12, cx).unwrap();
        let expected_velocity = proposal.velocity_mps().to_vec();
        let expected_kinetic = proposal.report().kinetic_after_j;
        let report = proposal.publish();
        assert_eq!(state.velocity_mps(), expected_velocity.as_slice());
        assert_eq!(state.kinetic_energy_j(), expected_kinetic);
        assert_eq!(
            state.external_work_j(),
            before_work + report.external_work_j
        );
        assert_eq!(
            state.material_history().displacement_m(),
            before_displacement.as_slice()
        );
        assert_eq!(state.material_history().energy(), before_energy);
        assert_eq!(state.time_s(), before_time);
        assert_eq!(state.steps(), before_steps);
    });
}
