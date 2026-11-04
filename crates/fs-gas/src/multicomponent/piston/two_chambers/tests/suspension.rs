use super::*;

fn pressure_pair(reference: f64) -> (MixtureEuler1d, MixtureEuler1d) {
    // conserved() takes pressure, independently of each source's gas constant.
    (
        chamber(left_gas(reference), 0.0, 1.0, 12, AREA, 2.0),
        chamber(right_gas(reference), 1.0, 1.0, 20, AREA, 1.0),
    )
}

#[test]
fn g1_preloaded_shared_equilibrium_and_zero_suspension_retain_native_fields() {
    let (mut left, mut right) = pressure_pair(0.0);
    let mut piston = Piston::new(0.4, 0.0).unwrap();
    let spring = PistonSuspension::new(4.0, 0.3, 0.875).unwrap();
    // The resolved pressure difference supplies 0.7 N, opposed by 0.5 N
    // spring preload and 0.2 N independently authored force.
    let report = left
        .step_with_shared_suspended_piston(
            &mut right,
            &mut piston,
            Time::new(0.001),
            StepControl::default(),
            -0.2,
            spring,
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(piston.velocity_m_s(), 0.0);
    assert_eq!(report.piston.displacement_m, 0.0);
    assert_eq!(report.spring_energy_j, [0.03125; 2]);
    assert_eq!(report.damping_energy_j, 0.0);
    let dt = report.piston.gas.left.elapsed.value();
    close(report.piston.gas_impulse_n_s, 0.7 * dt, 2e-15);
    close(report.suspension_impulse_n_s, -0.5 * dt, 2e-15);
    assert_eq!(left.mesh().positions()[left.cell_count()], 1.0);
    assert_eq!(right.mesh().positions()[0], 1.0);

    let (mut original_left, mut original_right) = pressure_pair(0.0);
    let (mut suspended_left, mut suspended_right) = pressure_pair(0.0);
    let mut a = Piston::new(0.4, 0.07).unwrap();
    let mut b = a;
    let free = PistonSuspension::new(0.0, 0.0, 8.0).unwrap();
    for _ in 0..4 {
        let expected = step(&mut original_left, &mut original_right, &mut a, || false).unwrap();
        let actual = suspended_left
            .step_with_shared_suspended_piston(
                &mut suspended_right,
                &mut b,
                Time::new(0.001),
                StepControl::default(),
                0.0,
                free,
                controls(),
                || false,
            )
            .unwrap();
        assert_eq!(actual.piston, expected);
        assert_eq!(original_left.states, suspended_left.states);
        assert_eq!(original_right.states, suspended_right.states);
        assert_eq!(original_left.mesh(), suspended_left.mesh());
        assert_eq!(original_right.mesh(), suspended_right.mesh());
        assert_eq!(a, b);
    }
}

#[test]
fn g1_both_pressures_close_suspended_motion_storage_species_and_independent_gauges() {
    let (mut left, mut right) = pressure_pair(0.0);
    let (mut shifted_left, mut shifted_right) = pressure_pair(1e18);
    let mut piston = Piston::new(0.4, 0.03).unwrap();
    let mut shifted_piston = piston;
    let spring = PistonSuspension::new(8.0, 0.8, 0.95).unwrap();
    let initial = [left.totals(), right.totals()];
    let initial_energy = initial[0].energy
        + initial[1].energy
        + 0.5 * 0.4 * 0.03_f64.powi(2)
        + 0.5 * 8.0 * 0.05_f64.powi(2);
    let mut loss = 0.0;
    let mut drive_work = 0.0;
    let mut impulse = 0.0;
    let mut support = 0.0;
    for _ in 0..60 {
        let x0 = left.mesh().positions()[left.cell_count()];
        let v0 = piston.velocity_m_s();
        let report = left
            .step_with_shared_suspended_piston(
                &mut right,
                &mut piston,
                Time::new(0.002),
                StepControl::default(),
                0.02,
                spring,
                controls(),
                || false,
            )
            .unwrap();
        let shifted = shifted_left
            .step_with_shared_suspended_piston(
                &mut shifted_right,
                &mut shifted_piston,
                Time::new(0.002),
                StepControl::default(),
                0.02,
                spring,
                controls(),
                || false,
            )
            .unwrap();
        assert_eq!(report, shifted);
        assert_eq!(left.states, shifted_left.states);
        assert_eq!(right.states, shifted_right.states);
        assert_eq!(piston, shifted_piston);
        let mechanical = &report.piston;
        let dt = mechanical.gas.left.elapsed.value();
        let x1 = left.mesh().positions()[left.cell_count()];
        let w = 0.5 * (v0 + piston.velocity_m_s());
        // Solve the discrete linear mechanics independently, using the actual
        // sum of the two source pressure impulses as the forcing term.
        let expected =
            (2.0 * 0.4 * v0 + mechanical.gas_impulse_n_s + dt * (0.02 - 8.0 * (x0 - 0.95)))
                / (2.0 * 0.4 + 0.8 * dt + 0.5 * 8.0 * dt * dt);
        close(w, expected, 3e-11);
        close(
            report.spring_energy_change_j,
            0.5 * 8.0 * ((x1 - 0.95).powi(2) - (x0 - 0.95).powi(2)),
            2e-15,
        );
        close(report.damping_energy_j, 0.8 * dt * w * w, 2e-16);
        assert_eq!(mechanical.exterior_pressure_work_j, 0.0);
        assert!(report.total_energy_residual_j.abs() <= controls().energy_tolerance_j);
        assert_eq!(x1, right.mesh().positions()[0]);
        assert_eq!(left.time(), right.time());
        loss += report.damping_energy_j;
        drive_work += 0.02 * (x1 - x0);
        impulse += 0.02 * dt + report.suspension_impulse_n_s;
        support += mechanical.fixed_wall_impulse_n_s;
    }
    let after = [left.totals(), right.totals()];
    let x = left.mesh().positions()[left.cell_count()];
    let final_energy = after[0].energy
        + after[1].energy
        + 0.5 * 0.4 * piston.velocity_m_s().powi(2)
        + 0.5 * 8.0 * (x - 0.95).powi(2);
    close(final_energy - initial_energy + loss, drive_work, 2e-9);
    close(
        after[0].momentum - initial[0].momentum + after[1].momentum - initial[1].momentum
            + 0.4 * (piston.velocity_m_s() - 0.03)
            + support,
        impulse,
        2e-10,
    );
    for (j, flow) in [&left, &right].into_iter().enumerate() {
        for k in 0..flow.gas().species_count() {
            close(
                after[j].species_masses[k],
                initial[j].species_masses[k],
                3e-13,
            );
        }
    }
    for i in 0..left.cell_count() {
        assert_eq!(left.cell(i).unwrap().partial_densities[1], 0.0);
    }
    assert!(loss > 1e-5);
    assert!(x > 1.001);
    assert!(pressure(&left, left.cell_count() - 1) < pressure(&left, 0));
    assert!(pressure(&right, 0) > pressure(&right, right.cell_count() - 1));
}
