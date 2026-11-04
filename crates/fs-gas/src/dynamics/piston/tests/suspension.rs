use super::*;

#[test]
fn zero_suspension_keeps_native_motion_and_preload_balances_pressure() {
    let mut original = field(16, 0.4, 0.0);
    let mut suspended = field(16, 0.4, 0.0);
    let mut a = Piston::new(0.4, 0.07).unwrap();
    let mut b = a;
    let free = PistonSuspension::new(0.0, 0.0, 12.0).unwrap();
    let load = PistonLoad {
        exterior_pressure_pa: 1.4,
        external_force_n: -0.1,
    };
    for _ in 0..12 {
        let expected = original
            .step_with_piston(
                &mut a,
                Time::new(0.003),
                StepControl::default(),
                load,
                controls(),
                || false,
            )
            .unwrap();
        let actual = suspended
            .step_with_suspended_piston(
                &mut b,
                Time::new(0.003),
                StepControl::default(),
                load,
                free,
                controls(),
                || false,
            )
            .unwrap();
        assert_eq!(actual.piston, expected);
        assert_eq!(a, b);
        assert_eq!(original.states(), suspended.states());
        assert_eq!(original.mesh(), suspended.mesh());
        assert_eq!(actual.spring_energy_j, [0.0; 2]);
        assert_eq!(actual.damping_energy_j, 0.0);
    }

    let mut balanced = field(16, 0.4, 0.0);
    let before = balanced.states().to_vec();
    let mut piston = Piston::new(0.4, 0.0).unwrap();
    // p*A=0.7 N, spring=-0.5 N and independent drive=-0.2 N.
    let spring = PistonSuspension::new(4.0, 0.3, 0.875).unwrap();
    let report = balanced
        .step_with_suspended_piston(
            &mut piston,
            Time::new(0.003),
            StepControl::default(),
            PistonLoad {
                exterior_pressure_pa: 0.0,
                external_force_n: -0.2,
            },
            spring,
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(piston.velocity_m_s(), 0.0);
    assert_eq!(report.piston.displacement_m, 0.0);
    assert_eq!(report.spring_energy_j, [0.03125; 2]);
    assert_eq!(report.spring_energy_change_j, 0.0);
    assert_eq!(report.damping_energy_j, 0.0);
    for (actual, expected) in balanced.states().iter().zip(before) {
        close(actual.density, expected.density, 1e-14);
        close(actual.energy, expected.energy, 1e-14);
    }
}

#[test]
fn pressure_spring_and_damper_share_midpoint_and_close_complete_trajectory() {
    let mut gas = field(24, 0.4, 0.0);
    let mut piston = Piston::new(0.7, 0.12).unwrap();
    let spring = PistonSuspension::new(8.0, 1.2, 0.92).unwrap();
    let load = PistonLoad {
        exterior_pressure_pa: 1.2,
        external_force_n: 0.05,
    };
    let initial = gas.totals();
    let initial_velocity = piston.velocity_m_s();
    let initial_energy = initial.energy
        + 0.5 * piston.mass_kg() * initial_velocity.powi(2)
        + 0.5 * 8.0 * (1.0_f64 - 0.92).powi(2);
    let mut heat = 0.0;
    let mut work = 0.0;
    let mut impulse = 0.0;
    let mut support = 0.0;
    for _ in 0..150 {
        let x0 = gas.mesh().positions()[24];
        let v0 = piston.velocity_m_s();
        let step = gas
            .step_with_suspended_piston(
                &mut piston,
                Time::new(0.004),
                StepControl::default(),
                load,
                spring,
                controls(),
                || false,
            )
            .unwrap();
        let report = &step.piston;
        let dt = report.gas.elapsed.value();
        let x1 = gas.mesh().positions()[24];
        let midpoint = 0.5 * (v0 + piston.velocity_m_s());
        // Independently solve the discrete linear mechanical equation with the
        // native field's integrated pressure impulse as its forcing term.
        let expected = (2.0 * 0.7 * v0
            + report.gas_impulse_n_s
            + dt * (0.05 - AREA * 1.2 - 8.0 * (x0 - 0.92)))
            / (2.0 * 0.7 + 1.2 * dt + 0.5 * 8.0 * dt * dt);
        close(midpoint, expected, 3e-11);
        close(
            step.spring_energy_change_j,
            0.5 * 8.0 * ((x1 - 0.92).powi(2) - (x0 - 0.92).powi(2)),
            2e-15,
        );
        close(step.damping_energy_j, 1.2 * dt * midpoint * midpoint, 2e-16);
        assert!(step.damping_energy_j >= 0.0);
        assert!(step.total_energy_residual_j.abs() <= controls().energy_tolerance_j);
        assert!(step.mechanical_energy_residual_j.abs() <= controls().energy_tolerance_j);
        close(step.applied_force_work_j, 0.05 * (x1 - x0), 1e-16);
        heat += step.damping_energy_j;
        work += step.applied_force_work_j - report.exterior_pressure_work_j;
        impulse += dt * (0.05 - AREA * 1.2) + step.suspension_impulse_n_s;
        support += report.fixed_wall_impulse_n_s;
    }
    let final_state = gas.totals();
    let x = gas.mesh().positions()[24];
    let final_energy = final_state.energy
        + 0.5 * 0.7 * piston.velocity_m_s().powi(2)
        + 0.5 * 8.0 * (x - 0.92).powi(2);
    close(final_energy - initial_energy + heat, work, 2e-9);
    close(final_state.mass, initial.mass, 2e-12);
    close(
        final_state.momentum - initial.momentum + 0.7 * (piston.velocity_m_s() - initial_velocity),
        impulse - support,
        2e-9,
    );
    assert!(heat > 1e-4);
    assert!((x - 1.0).abs() > 1e-3);
}

#[test]
fn suspension_refusal_and_final_cancellation_keep_accepted_state() {
    for (k, c, rest) in [
        (-1.0, 0.0, 1.0),
        (0.0, -1.0, 1.0),
        (f64::NAN, 1.0, 1.0),
        (1.0, 1.0, f64::INFINITY),
    ] {
        assert!(PistonSuspension::new(k, c, rest).is_err());
    }
    let spring = PistonSuspension::new(8.0, 1.2, 0.92).unwrap();
    let load = PistonLoad {
        exterior_pressure_pa: 1.2,
        external_force_n: 0.05,
    };
    let mut good = field(16, 0.4, 0.0);
    let mut moving = Piston::new(0.7, 0.12).unwrap();
    let mut polls = 0;
    good.step_with_suspended_piston(
        &mut moving,
        Time::new(0.004),
        StepControl::default(),
        load,
        spring,
        controls(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    for cancellation_at in [1, polls] {
        let mut gas = field(16, 0.4, 0.0);
        let mut piston = Piston::new(0.7, 0.12).unwrap();
        let before_piston = piston;
        let before = gas.states().to_vec();
        let mesh = gas.mesh().clone();
        let mut observed = 0;
        assert_eq!(
            gas.step_with_suspended_piston(
                &mut piston,
                Time::new(0.004),
                StepControl::default(),
                load,
                spring,
                controls(),
                || {
                    observed += 1;
                    observed == cancellation_at
                },
            ),
            Err(PistonError::Cancelled)
        );
        assert_eq!(gas.states(), before);
        assert_eq!(gas.mesh(), &mesh);
        assert_eq!(gas.time().value(), 0.0);
        assert_eq!(piston, before_piston);
    }
}
