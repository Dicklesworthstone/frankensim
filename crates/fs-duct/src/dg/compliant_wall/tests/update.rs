use super::*;

#[test]
fn live_wall_transfer_follows_the_analytic_oscillator_and_keeps_one_work_account() {
    let old = wall(0.0);
    let new = WallPin {
        surface_density: 0.5,
        stiffness_per_area: 4e6,
        resistance: 0.0,
    };
    let mut dg = duct(DgBoundary::Closed).with_compliant_wall(old).unwrap();
    let dt = 2e-6;
    let original_propagator = dg.propagator(dt, 1e-13).unwrap();
    let mut state = dg.project_initial_state(|_| 10.0, |_| 0.0);
    for _ in 0..30 {
        original_propagator.step(&mut state, &[]).unwrap();
    }
    original_propagator.sync_nodal(&mut state).unwrap();
    let before = state.clone();
    let motion = dg.compliant_wall_state(&state).unwrap();
    let update = dg
        .rebind_compliant_wall_regions(&mut state, &[new; 2], 1.0, &mut || false)
        .unwrap();
    let transferred = dg.compliant_wall_state(&state).unwrap();
    assert_eq!(state.pressure(), before.pressure());
    assert_eq!(state.volume_flow(), before.volume_flow());
    assert_eq!(state.step_count(), 30);
    assert_eq!(state.initial_energy_j(), before.initial_energy_j());
    assert_eq!(
        state.cumulative_dissipated_j(),
        before.cumulative_dissipated_j()
    );
    for i in 0..dg.field_dofs() {
        close(
            transferred.displacement_m[i],
            motion.displacement_m[i],
            1e-6,
        );
        close(
            new.surface_density * transferred.velocity_m_s[i],
            old.surface_density * motion.velocity_m_s[i],
            1e-3,
        );
    }
    let surface = core::f64::consts::TAU * RADIUS * LENGTH;
    let expected_work = 0.5
        * surface
        * ((new.stiffness_per_area - old.stiffness_per_area) * motion.displacement_m[0].powi(2)
            + (old.surface_density.powi(2) / new.surface_density - old.surface_density)
                * motion.velocity_m_s[0].powi(2));
    close(
        update.parameter_work_j,
        expected_work,
        before.stored_energy_j(),
    );
    close(
        update.energy_after_j - update.energy_before_j,
        expected_work,
        before.stored_energy_j(),
    );
    assert!(update.parameter_work_j < 0.0);
    assert_eq!(state.cumulative_parameter_work_j(), update.parameter_work_j);
    let saved = state.clone();
    assert!(original_propagator.step(&mut state, &[]).is_err());
    assert_eq!(state, saved);

    // For a spatially uniform sealed cylinder, C=p+beta*w is invariant and
    // sigma*w''+(K+beta)*w=C. Use its independent initial-value solution.
    let beta = 2.0 * dg.gas().density * dg.gas().sound_speed.powi(2) / RADIUS;
    let w0 = transferred.displacement_m[0];
    let v0 = transferred.velocity_m_s[0];
    let c = state.pressure()[0] + beta * w0;
    let equilibrium = c / (new.stiffness_per_area + beta);
    let omega = ((new.stiffness_per_area + beta) / new.surface_density).sqrt();
    let propagation = dg.propagator(dt, 1e-13).unwrap();
    for step in 1..=30 {
        let frame = propagation.step(&mut state, &[]).unwrap();
        propagation.sync_nodal(&mut state).unwrap();
        let (sin, cos) = (omega * dt * step as f64).sin_cos();
        let expected_w = equilibrium + (w0 - equilibrium) * cos + v0 / omega * sin;
        let expected_v = -(w0 - equilibrium) * omega * sin + v0 * cos;
        let expected_p = c - beta * expected_w;
        let motion = dg.compliant_wall_state(&state).unwrap();
        for i in 0..dg.field_dofs() {
            close(state.pressure()[i], expected_p, 10.0);
            close(motion.displacement_m[i], expected_w, 1e-6);
            close(motion.velocity_m_s[i], expected_v, 1e-3);
        }
        assert_eq!(
            frame.ledger.cumulative_parameter_work_j,
            update.parameter_work_j
        );
        assert_eq!(frame.ledger.cumulative_port_work_j, 0.0);
        close(
            frame.ledger.energy_balance_residual_j,
            0.0,
            before.initial_energy_j(),
        );
        close(
            state.stored_energy_j(),
            update.energy_after_j,
            before.initial_energy_j(),
        );
    }
}

#[test]
fn regional_update_retains_wave_and_diffusion_history_and_measures_physical_work() {
    let old = [
        wall(30.0),
        WallPin {
            surface_density: 0.7,
            stiffness_per_area: 4e7,
            resistance: 80.0,
        },
    ];
    let new = [
        WallPin {
            surface_density: 0.3,
            stiffness_per_area: 2e7,
            resistance: 40.0,
        },
        WallPin {
            surface_density: 0.5,
            stiffness_per_area: 3e7,
            resistance: 120.0,
        },
    ];
    let mut dg = duct(DgBoundary::Port(PortId(0)))
        .with_viscothermal(8, [50.0, 2000.0], 0.01)
        .unwrap()
        .with_compliant_wall_regions(&old)
        .unwrap();
    let mut state = dg.project_initial_state(|x| 10.0 + x, |x| 1e-7 * (1.0 + x));
    // A legitimate nonzero initial memory makes an accidental cold restart
    // observable without requiring another dynamics solver to create it.
    for (i, value) in state.tilde[2 * dg.field_dofs()..].iter_mut().enumerate() {
        *value = 1e-6 * ((i * 17 % 23) as f64 - 11.0);
    }
    state.initial_energy_j = state.stored_energy_j();
    let original = dg.clone();
    let before = state.clone();
    let old_motion = dg.compliant_wall_state(&state).unwrap();
    let prefix = dg.compliant_wall.as_ref().unwrap().start;
    let update = dg
        .rebind_compliant_wall_regions(&mut state, &new, 1.0, &mut || false)
        .unwrap();
    let motion = dg.compliant_wall_state(&state).unwrap();
    assert_eq!(&state.tilde[..prefix], &before.tilde[..prefix]);
    assert_eq!(dg.ports, original.ports);
    assert_eq!(dg.viscothermal_model(), original.viscothermal_model());
    let n = dg.total_dofs();
    for row in 0..prefix {
        assert_eq!(
            &dg.l_tilde[row * n..row * n + prefix],
            &original.l_tilde[row * n..row * n + prefix]
        );
    }
    let local = dg.mesh().order + 1;
    for i in 0..dg.field_dofs() {
        let e = i / local;
        close(motion.displacement_m[i], old_motion.displacement_m[i], 1e-4);
        close(
            new[e].surface_density * motion.velocity_m_s[i],
            old[e].surface_density * old_motion.velocity_m_s[i],
            1.0,
        );
    }
    let (points, weights) = gauss_legendre_rule(8);
    let mut expected_work = 0.0;
    for (e, element) in dg.elements.iter().enumerate() {
        for (&point, &weight) in points.iter().zip(&weights) {
            let shape = dg.lgl.eval_basis(point);
            let interp = |values: &[f64]| -> f64 {
                (0..local).map(|i| shape[i] * values[e * local + i]).sum()
            };
            let w = interp(&old_motion.displacement_m);
            let momentum = old[e].surface_density * interp(&old_motion.velocity_m_s);
            expected_work += 0.5
                * core::f64::consts::TAU
                * RADIUS
                * 0.5
                * element.length
                * weight
                * ((new[e].stiffness_per_area - old[e].stiffness_per_area) * w * w
                    + (1.0 / new[e].surface_density - 1.0 / old[e].surface_density)
                        * momentum
                        * momentum);
        }
    }
    close(
        update.parameter_work_j,
        expected_work,
        before.stored_energy_j(),
    );
    let accepted = state.clone();
    let next = new.map(|pin| WallPin {
        resistance: 4.0 * pin.resistance,
        ..pin
    });
    let viscous = dg
        .rebind_compliant_wall_regions(&mut state, &next, 0.0, &mut || false)
        .unwrap();
    assert_eq!(viscous.parameter_work_j, 0.0);
    assert_eq!(state.tilde, accepted.tilde);
    assert_eq!(state.cumulative_parameter_work_j(), update.parameter_work_j);
    close(
        dg.compliant_wall_state(&state).unwrap().dissipated_power_w,
        4.0 * motion.dissipated_power_w,
        motion.dissipated_power_w,
    );
}

#[test]
fn work_budget_foreign_state_history_and_late_cancellation_refuse_without_publication() {
    let original = duct(DgBoundary::Closed)
        .with_compliant_wall(wall(0.0))
        .unwrap();
    let mut initial = original.project_initial_state(|_| 10.0, |_| 0.0);
    let propagation = original.propagator(2e-6, 1e-13).unwrap();
    for _ in 0..10 {
        propagation.step(&mut initial, &[]).unwrap();
    }
    let replacement = [WallPin {
        stiffness_per_area: 2e7,
        ..wall(0.0)
    }; 2];
    for budget in [0.0, f64::NAN, -1.0] {
        let mut dg = original.clone();
        let mut state = initial.clone();
        assert!(
            dg.rebind_compliant_wall_regions(&mut state, &replacement, budget, &mut || false)
                .is_err()
        );
        assert_eq!(dg.operator_identity(), original.operator_identity());
        assert_eq!(state, initial);
    }
    let mut dg = original.clone();
    let mut foreign = duct(DgBoundary::Closed).zero_state();
    let saved = foreign.clone();
    assert!(
        dg.rebind_compliant_wall_regions(&mut foreign, &replacement, 1.0, &mut || false)
            .is_err()
    );
    assert_eq!(foreign, saved);
    assert_eq!(dg.operator_identity(), original.operator_identity());
    let mut polls = 0;
    let mut reference = original.clone();
    reference
        .rebind_compliant_wall_regions(&mut initial.clone(), &replacement, 1.0, &mut || {
            polls += 1;
            false
        })
        .unwrap();
    let mut state = initial.clone();
    let mut calls = 0;
    assert!(matches!(
        dg.rebind_compliant_wall_regions(&mut state, &replacement, 1.0, &mut || {
            calls += 1;
            calls == polls
        }),
        Err(DuctError::Cancelled { .. })
    ));
    assert_eq!(state, initial);
    assert_eq!(dg.operator_identity(), original.operator_identity());
    let material = fs_material::visco::GeneralizedMaxwell::new(1e6, vec![(5e5, 0.001)]).unwrap();
    let uniform = duct(DgBoundary::Closed)
        .with_compliant_wall_relaxation(wall(0.0), &material)
        .unwrap();
    let regional = duct(DgBoundary::Closed)
        .with_compliant_wall_relaxation_regions(
            &[
                wall(0.0),
                WallPin {
                    surface_density: 0.7,
                    ..wall(0.0)
                },
            ],
            &[None, Some(&material)],
        )
        .unwrap();
    for mut hereditary in [uniform, regional] {
        let mut memory = hereditary.zero_state();
        for branch in &hereditary.compliant_wall.as_ref().unwrap().relaxation {
            memory.tilde[branch.start..branch.start + branch.dofs].fill(1e-6);
        }
        memory.initial_energy_j = memory.stored_energy_j();
        assert!(
            hereditary
                .compliant_wall_state(&memory)
                .unwrap()
                .relaxation_stored_energy_j
                > 0.0
        );
        let saved = memory.clone();
        let identity = hereditary.operator_identity();
        assert!(
            hereditary
                .rebind_compliant_wall_regions(&mut memory, &replacement, 1.0, &mut || false)
                .is_err()
        );
        assert_eq!(memory, saved);
        assert_eq!(hereditary.operator_identity(), identity);
    }
}
