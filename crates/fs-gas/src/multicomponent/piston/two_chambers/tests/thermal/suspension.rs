//! Dissipation has one mechanical debit and one declared thermal destination.

use super::*;

fn suspended(
    left: &mut MixtureEuler1d,
    right: &mut MixtureEuler1d,
    wall: &mut SharedPistonThermal<'_>,
    piston: &mut Piston,
    duration: f64,
    fraction: f64,
    cancelled: impl FnMut() -> bool,
) -> Result<SuspendedPistonStep<SharedPistonThermalStepReport>, SharedPistonThermalError> {
    left.step_with_shared_suspended_piston_and_heat(
        right,
        wall,
        piston,
        Time::new(duration),
        StepControl::default(),
        heat_controls(),
        0.02,
        PistonSuspension::new(8.0, 2.5, 0.95).unwrap(),
        fraction,
        controls(),
        cancelled,
    )
}

#[test]
fn g1_damping_heat_enters_one_wall_and_drives_the_original_pressure_field() {
    let mut endpoints = Vec::new();
    for fraction in [0.0, 0.35, 1.0] {
        let (mut left, mut right) = fields([2, 2], 1e18);
        let g = [9.0, 0.2];
        let mut wall = finite(&left, &right, 1.4, g);
        let mut piston = Piston::new(0.4, 0.35).unwrap();
        let initial_energy = left.totals().energy
            + right.totals().energy
            + 0.6 * 1.4
            + 0.5 * 0.4 * 0.35_f64.powi(2)
            + 0.5 * 8.0 * 0.05_f64.powi(2);
        let (mut external_loss, mut drive) = (0.0, 0.0);
        for _ in 0..16 {
            let x0 = left.mesh().positions()[2];
            let v0 = piston.velocity_m_s();
            let t0 = [
                end_temperature(&left, 1, 0.0),
                end_temperature(&right, 0, 0.0),
            ];
            let tw0 = wall_temperature(&wall);
            let capacity = [
                2.5 * left.cell(1).unwrap().density() * left.mesh().volumes()[1],
                1.75 * right.cell(0).unwrap().density() * right.mesh().volumes()[0],
            ];
            let report = suspended(
                &mut left,
                &mut right,
                &mut wall,
                &mut piston,
                0.005,
                fraction,
                || false,
            )
            .unwrap();
            let p = &report.piston;
            let dt = p.gas.transport.left.elapsed.value();
            let w = 0.5 * (v0 + piston.velocity_m_s());
            // Independent linear midpoint mechanics, forced by the actual two
            // resolved gas pressure impulses, with authored k, c and anchor.
            let expected_w = (2.0 * 0.4 * v0 + p.gas_impulse_n_s + dt * (0.02 - 8.0 * (x0 - 0.95)))
                / (2.0 * 0.4 + 2.5 * dt + 0.5 * 8.0 * dt * dt);
            close(w, expected_w, 4e-11);
            let damping = 2.5 * w * w * dt;
            let returned = fraction * damping;
            close(report.damping_energy_j, damping, 3e-14);
            close(p.gas.heat.damping_heat_to_wall_j, returned, 3e-14);
            close(
                p.gas.heat.damping_loss_external_j,
                (1.0 - fraction) * damping,
                3e-14,
            );
            close(p.gas.heat.damping_heat_residual_j, 0.0, 3e-14);
            // Each H is checked against a separate three-capacity algebraic
            // solution. H2 starts AFTER deposition and uses moved gas amounts.
            let (tw1, q1) = linear_heat(t0, capacity, tw0, g, 0.5 * dt);
            let actual_q2 = [
                p.gas.after.left.energy - p.gas.transport.left.after.energy,
                p.gas.after.right.energy - p.gas.transport.right.after.energy,
            ];
            let tm = [
                end_temperature(&left, 1, actual_q2[0]),
                end_temperature(&right, 0, actual_q2[1]),
            ];
            let moved_capacity = [
                2.5 * left.cell(1).unwrap().density() * left.mesh().volumes()[1],
                1.75 * right.cell(0).unwrap().density() * right.mesh().volumes()[0],
            ];
            let (tw2, q2) = linear_heat(tm, moved_capacity, tw1 + returned / 0.6, g, 0.5 * dt);
            for (i, (before, middle)) in [
                (p.gas.before.left, p.gas.transport.left),
                (p.gas.before.right, p.gas.transport.right),
            ]
            .into_iter()
            .enumerate()
            {
                close(middle.before.energy - before.energy, q1[i], 3e-12);
                close(actual_q2[i], q2[i], 3e-12);
            }
            close(wall_temperature(&wall), tw2, 1e-11);
            close(p.gas.heat.wall_energy_change_j, 0.6 * (tw2 - tw0), 1e-11);
            assert!(p.admitted(controls()));
            assert!(report.total_energy_residual_j.abs() < controls().energy_tolerance_j);
            close(p.gas.heat.energy_residual_j, 0.0, 1e-12);
            external_loss += (1.0 - fraction) * damping;
            drive += 0.02 * (left.mesh().positions()[2] - x0);
        }
        let final_energy = left.totals().energy
            + right.totals().energy
            + 0.6 * wall_temperature(&wall)
            + 0.5 * 0.4 * piston.velocity_m_s().powi(2)
            + 0.5 * 8.0 * (left.mesh().positions()[2] - 0.95).powi(2);
        close(
            final_energy - initial_energy + external_loss - drive,
            0.0,
            2e-10,
        );
        endpoints.push([
            wall_temperature(&wall),
            pressure(&left, 1),
            piston.velocity_m_s(),
        ]);
    }
    // Returned dissipation is a real thermodynamic source; it changes pressure
    // and the following pressure/inertia root rather than only its receipt.
    for component in 0..3 {
        assert!(
            endpoints[1][component] > endpoints[0][component] + 1e-6,
            "{endpoints:?}"
        );
        assert!(
            endpoints[2][component] > endpoints[1][component] + 1e-6,
            "{endpoints:?}"
        );
    }
}

#[test]
fn g0_insulated_and_zero_suspension_limits_retain_native_trajectories() {
    let (mut left, mut right) = fields([2, 2], 0.0);
    let (mut native_left, mut native_right) = fields([2, 2], 0.0);
    let mut wall = finite(&left, &right, 3.0, [0.0; 2]);
    let mut piston = Piston::new(0.4, 0.2).unwrap();
    let mut native_piston = piston;
    for _ in 0..4 {
        let native = native_left
            .step_with_shared_suspended_piston(
                &mut native_right,
                &mut native_piston,
                Time::new(0.01),
                StepControl::default(),
                0.02,
                PistonSuspension::new(8.0, 2.5, 0.95).unwrap(),
                controls(),
                || false,
            )
            .unwrap();
        let actual = suspended(
            &mut left,
            &mut right,
            &mut wall,
            &mut piston,
            0.01,
            0.0,
            || false,
        )
        .unwrap();
        assert_eq!(actual.piston.gas.transport, native.piston.gas);
        assert_eq!(actual.damping_energy_j, native.damping_energy_j);
        assert_eq!(actual.spring_energy_j, native.spring_energy_j);
        assert_eq!(left.states, native_left.states);
        assert_eq!(right.states, native_right.states);
        assert_eq!(left.mesh(), native_left.mesh());
        assert_eq!(right.mesh(), native_right.mesh());
        assert_eq!(piston, native_piston);
        assert_eq!(wall_temperature(&wall), 3.0);
    }
    let (mut left, mut right) = fields([2, 2], 0.0);
    let (mut native_left, mut native_right) = fields([2, 2], 0.0);
    let mut wall = finite(&left, &right, 3.0, [1.7, 0.4]);
    let mut native_wall = finite(&native_left, &native_right, 3.0, [1.7, 0.4]);
    let mut piston = Piston::new(0.4, 0.2).unwrap();
    let mut native_piston = piston;
    let native = advance(
        &mut native_left,
        &mut native_right,
        &mut native_wall,
        &mut native_piston,
        0.01,
        || false,
    )
    .unwrap();
    let actual = left
        .step_with_shared_suspended_piston_and_heat(
            &mut right,
            &mut wall,
            &mut piston,
            Time::new(0.01),
            StepControl::default(),
            heat_controls(),
            0.0,
            PistonSuspension::new(0.0, 0.0, 10.0).unwrap(),
            1.0,
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(actual.piston, native);
    assert_eq!(left.states, native_left.states);
    assert_eq!(right.states, native_right.states);
    assert_eq!(wall.storage(), native_wall.storage());
    assert_eq!(piston, native_piston);
}

#[test]
fn g1_damper_deposition_uses_original_nonlinear_wall_energy_and_cancels_its_datum() {
    let mut endpoints = Vec::new();
    for datum in [0.0, 1024.0] {
        let (mut left, mut right) = fields([2, 2], datum * 1e12);
        let mut wall = SharedPistonThermal::new(
            &left,
            &right,
            source_storage(datum, ThermalStorageBasis::InternalEnergy),
            [0.0; 2],
            SharedWallControl::default(),
            || false,
        )
        .unwrap();
        let mut piston = Piston::new(0.4, 0.7).unwrap();
        let mut returned = 0.0;
        for _ in 0..6 {
            let initial_t = wall_temperature(&wall);
            let report = left
                .step_with_shared_suspended_piston_and_heat(
                    &mut right,
                    &mut wall,
                    &mut piston,
                    Time::new(0.01),
                    StepControl::default(),
                    heat_controls(),
                    0.0,
                    PistonSuspension::new(4.0, 20.0, 1.0).unwrap(),
                    0.75,
                    controls(),
                    || false,
                )
                .unwrap();
            let heat = report.piston.gas.heat;
            let expected_q = 0.75 * report.damping_energy_j;
            returned += expected_q;
            close(
                wall_temperature(&wall),
                (initial_t.powi(2) + expected_q / 0.08).sqrt(),
                5e-12,
            );
            close(heat.wall_energy_change_j, expected_q, 1e-12);
            assert_eq!(heat.heat_to_gas_j, [0.0; 2]);
            assert!(report.total_energy_residual_j.abs() < controls().energy_tolerance_j);
        }
        let final_t = wall_temperature(&wall);
        close(final_t, (9.0 + returned / 0.08).sqrt(), 5e-12);
        let frozen_capacity_result = 3.0 + returned / (0.16 * 3.0);
        assert!((final_t - frozen_capacity_result).abs() > 1e-4);
        endpoints.push((left.states, right.states, piston, final_t));
    }
    assert_eq!(endpoints[0].0, endpoints[1].0);
    assert_eq!(endpoints[0].1, endpoints[1].1);
    assert_eq!(endpoints[0].2, endpoints[1].2);
    close(endpoints[0].3, endpoints[1].3, 5e-12);
}

#[test]
fn g4_shortening_and_late_refusal_never_deposit_rejected_damping_heat() {
    let setup = || {
        let (mut left, mut right) = fields([12, 20], 0.0);
        let mut wall = finite(&left, &right, 10.0, [100.0; 2]);
        let mut piston = Piston::new(0.4, 0.2).unwrap();
        suspended(
            &mut left,
            &mut right,
            &mut wall,
            &mut piston,
            0.001,
            0.6,
            || false,
        )
        .unwrap();
        (left, right, wall, piston)
    };
    let (mut left, mut right, mut wall, mut piston) = setup();
    let mut polls = 0;
    let accepted = suspended(
        &mut left,
        &mut right,
        &mut wall,
        &mut piston,
        1.0,
        0.6,
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    assert!(accepted.piston.gas.restarts > 0);
    let (mut replay_left, mut replay_right, mut replay_wall, mut replay_piston) = setup();
    let replay = suspended(
        &mut replay_left,
        &mut replay_right,
        &mut replay_wall,
        &mut replay_piston,
        accepted.piston.gas.transport.left.elapsed.value(),
        0.6,
        || false,
    )
    .unwrap();
    assert_eq!(left.states, replay_left.states);
    assert_eq!(right.states, replay_right.states);
    assert_eq!(wall.storage(), replay_wall.storage());
    assert_eq!(piston, replay_piston);
    assert_eq!(accepted.piston.gas.heat, replay.piston.gas.heat);
    let (mut failed_left, mut failed_right, mut failed_wall, mut failed_piston) = setup();
    let before = [failed_left.states.clone(), failed_right.states.clone()];
    let meshes = [failed_left.mesh().clone(), failed_right.mesh().clone()];
    let original_piston = failed_piston;
    let old_storage = failed_wall.storage().clone();
    let old_report = failed_wall.last_report();
    let old_time = failed_left.time();
    let mut calls = 0;
    assert_eq!(
        suspended(
            &mut failed_left,
            &mut failed_right,
            &mut failed_wall,
            &mut failed_piston,
            1.0,
            0.6,
            || {
                calls += 1;
                calls == polls
            },
        )
        .unwrap_err(),
        SharedPistonThermalError::Cancelled
    );
    assert_eq!(failed_left.states, before[0]);
    assert_eq!(failed_right.states, before[1]);
    assert_eq!(failed_left.mesh(), &meshes[0]);
    assert_eq!(failed_right.mesh(), &meshes[1]);
    assert_eq!(failed_piston, original_piston);
    assert_eq!(failed_wall.storage(), &old_storage);
    assert_eq!(failed_wall.last_report(), old_report);
    assert_eq!(failed_left.time(), old_time);
    assert_eq!(failed_right.time(), old_time);
    let retried = suspended(
        &mut failed_left,
        &mut failed_right,
        &mut failed_wall,
        &mut failed_piston,
        1.0,
        0.6,
        || false,
    )
    .unwrap();
    assert_eq!(retried, accepted);
    assert_eq!(failed_left.states, left.states);
    assert_eq!(failed_right.states, right.states);
    assert_eq!(failed_wall.storage(), wall.storage());

    let (mut left, mut right) = fields([2, 2], 0.0);
    let mut full = finite(&left, &right, 3.0, [1.7, 0.4]);
    let mut piston = Piston::new(0.4, 0.2).unwrap();
    let result = suspended(
        &mut left,
        &mut right,
        &mut full,
        &mut piston,
        0.01,
        0.6,
        || false,
    )
    .unwrap();
    let (mut left, mut right) = fields([2, 2], 0.0);
    let mut wall = SharedPistonThermal::new(
        &left,
        &right,
        finite(&left, &right, 3.0, [1.7, 0.4]).storage().clone(),
        [1.7, 0.4],
        SharedWallControl {
            maximum_source_evaluations: result.piston.gas.heat.source_evaluations - 1,
        },
        || false,
    )
    .unwrap();
    let before = [left.states.clone(), right.states.clone()];
    let mut piston = Piston::new(0.4, 0.2).unwrap();
    let original = piston;
    assert!(
        left.step_with_shared_suspended_piston_and_heat(
            &mut right,
            &mut wall,
            &mut piston,
            Time::new(0.01),
            StepControl {
                maximum_halvings: 0,
                ..StepControl::default()
            },
            heat_controls(),
            0.02,
            PistonSuspension::new(8.0, 2.5, 0.95).unwrap(),
            0.6,
            controls(),
            || false,
        )
        .is_err()
    );
    assert_eq!(left.states, before[0]);
    assert_eq!(right.states, before[1]);
    assert_eq!(wall_temperature(&wall), 3.0);
    assert_eq!(wall.last_report(), None);
    assert_eq!(piston, original);
    for fraction in [-0.1, 1.1, f64::NAN] {
        assert!(
            suspended(
                &mut left,
                &mut right,
                &mut wall,
                &mut piston,
                0.01,
                fraction,
                || false
            )
            .is_err()
        );
        assert_eq!(left.states, before[0]);
        assert_eq!(right.states, before[1]);
        assert_eq!(wall_temperature(&wall), 3.0);
        assert_eq!(piston, original);
    }
}

#[test]
fn g0_material_contacts_and_returned_damper_heat_keep_separate_energy_accounts() {
    let (mut left, mut right) = contacts();
    let mut wall = finite(&left, &right, 3.0, [1.7, 0.4]);
    let mut piston = Piston::new(0.4, 0.25).unwrap();
    let mut workspace = [64, 80].map(|n| {
        ContactWorkspace::new(
            n,
            ContactWorkspace::required_workspace_bytes(n).unwrap(),
            || false,
        )
        .unwrap()
    });
    let equilibrium = [ContactEquilibriumControl {
        relative_tolerance: 1e-10,
        maximum_energy_defect_j: 1.0,
    }; 2];
    let mut defect = 0.0;
    for _ in 0..8 {
        let [a, b] = &mut workspace;
        let report = left
            .step_with_shared_suspended_piston_heat_and_contact_equilibrium(
                &mut right,
                &mut wall,
                &mut piston,
                Time::new(0.001),
                StepControl::default(),
                heat_controls(),
                0.02,
                PistonSuspension::new(8.0, 0.8, 0.95).unwrap(),
                0.75,
                controls(),
                equilibrium,
                [a, b],
                || false,
            )
            .unwrap();
        let p = &report.piston;
        let correction = sum([
            p.gas.transport.left.correction.signed_energy_j,
            p.gas.transport.right.correction.signed_energy_j,
        ]);
        defect += correction.abs();
        let residual = sum([
            p.gas.after.left.energy - p.gas.before.left.energy,
            p.gas.after.right.energy - p.gas.before.right.energy,
            p.gas.heat.wall_energy_change_j,
            p.kinetic_energy_change_j,
            report.spring_energy_change_j,
            0.25 * report.damping_energy_j,
            -report.applied_force_work_j,
            -correction,
        ]);
        close(residual, report.total_energy_residual_j, 1e-14);
        assert!(residual.abs() < controls().energy_tolerance_j);
        // Counting the returned heat again as a loss must fail this same gate.
        assert!(
            (residual + p.gas.heat.damping_heat_to_wall_j).abs() > controls().energy_tolerance_j
        );
        close(p.gas.heat.energy_residual_j, 0.0, 1e-12);
        assert!(p.admitted(controls()));
    }
    assert!(defect > 1e-9);
    for (side, flow) in [&left, &right].into_iter().enumerate() {
        for cell in flow.cell_count() / 2 - 3..flow.cell_count() / 2 + 3 {
            close(pressure(flow, cell), [1.2, 0.8][side], 3e-10);
        }
    }
    assert!((pressure(&left, 63) - 1.2).abs() > 1e-4);
    let before = [left.states.clone(), right.states.clone()];
    let storage = wall.storage().clone();
    let last = wall.last_report();
    let original = piston;
    let [a, b] = &mut workspace;
    let mut refusing = equilibrium;
    refusing[1].maximum_energy_defect_j = 0.0;
    assert!(
        left.step_with_shared_suspended_piston_heat_and_contact_equilibrium(
            &mut right,
            &mut wall,
            &mut piston,
            Time::new(0.001),
            StepControl {
                maximum_halvings: 0,
                ..StepControl::default()
            },
            heat_controls(),
            0.02,
            PistonSuspension::new(8.0, 0.8, 0.95).unwrap(),
            0.75,
            controls(),
            refusing,
            [a, b],
            || false,
        )
        .is_err()
    );
    assert_eq!(left.states, before[0]);
    assert_eq!(right.states, before[1]);
    assert_eq!(wall.storage(), &storage);
    assert_eq!(wall.last_report(), last);
    assert_eq!(piston, original);
}

#[test]
fn g2_pressure_motion_and_returned_dissipation_converge_under_time_refinement() {
    let mut endpoints = Vec::new();
    for steps in [4, 8, 16, 64] {
        let (mut left, mut right) = fields([2, 2], 0.0);
        let mut wall = finite(&left, &right, 3.0, [30.0, 7.0]);
        let mut piston = Piston::new(0.4, 0.35).unwrap();
        let dt = 0.04 / steps as f64;
        for _ in 0..steps {
            let report = suspended(
                &mut left,
                &mut right,
                &mut wall,
                &mut piston,
                dt,
                0.6,
                || false,
            )
            .unwrap();
            assert_eq!(report.piston.gas.transport.left.elapsed.value(), dt);
            assert_eq!(report.piston.gas.restarts, 0);
        }
        endpoints.push([
            piston.velocity_m_s(),
            pressure(&left, 1),
            pressure(&right, 0),
            wall_temperature(&wall),
        ]);
    }
    // Fixed spatial resolution: a converging split time integrator, with no
    // second-order or simultaneous thermomechanical convergence claim.
    for component in 0..4 {
        let error = |i: usize| (endpoints[i][component] - endpoints[3][component]).abs();
        assert!(
            error(1) < 0.8 * error(0),
            "component {component}: {endpoints:?}"
        );
        assert!(
            error(2) < 0.8 * error(1),
            "component {component}: {endpoints:?}"
        );
    }
}
