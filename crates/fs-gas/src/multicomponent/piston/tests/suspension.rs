use super::*;

#[test]
fn suspended_mixture_retains_species_formation_gauges_and_complete_work() {
    let mut gas = field(24, [0.0; 3]);
    let mut shifted = field(24, [1e18, -7e17, 3e17]);
    let mut piston = Piston::new(0.7, 0.08).unwrap();
    let mut other = piston;
    let before = gas.totals();
    let spring = PistonSuspension::new(8.0, 1.2, 0.92).unwrap();
    let load = PistonLoad {
        exterior_pressure_pa: 1.3,
        external_force_n: 0.05,
    };
    let initial_energy =
        before.energy + 0.5 * 0.7 * 0.08_f64.powi(2) + 0.5 * 8.0 * (1.0_f64 - 0.92).powi(2);
    let mut heat = 0.0;
    let mut work = 0.0;
    for _ in 0..40 {
        let step = gas
            .step_with_suspended_piston(
                &mut piston,
                Time::new(0.003),
                StepControl::default(),
                load,
                spring,
                controls(),
                || false,
            )
            .unwrap();
        let reference = shifted
            .step_with_suspended_piston(
                &mut other,
                Time::new(0.003),
                StepControl::default(),
                load,
                spring,
                controls(),
                || false,
            )
            .unwrap();
        assert_eq!(step, reference);
        assert_eq!(piston, other);
        assert_eq!(gas.states, shifted.states);
        assert_eq!(gas.mesh(), shifted.mesh());
        for k in 0..3 {
            assert_eq!(step.piston.gas.boundary_change.species_masses[k], 0.0);
            close(
                step.piston.gas.conservation_residual.species_masses[k],
                0.0,
                4e-13,
            );
        }
        for state in gas.states.chunks_exact(5) {
            assert_eq!(state[2], 0.0);
        }
        heat += step.damping_energy_j;
        work += step.applied_force_work_j - step.piston.exterior_pressure_work_j;
        assert!(step.total_energy_residual_j.abs() <= controls().energy_tolerance_j);
    }
    let after = gas.totals();
    let x = gas.mesh().positions()[24];
    let final_energy =
        after.energy + 0.5 * 0.7 * piston.velocity_m_s().powi(2) + 0.5 * 8.0 * (x - 0.92).powi(2);
    close(final_energy - initial_energy + heat, work, 6e-10);
    for k in 0..3 {
        close(after.species_masses[k], before.species_masses[k], 3e-12);
    }
    assert!(heat > 0.0);
}

#[test]
fn zero_mixture_suspension_is_exact_and_cancellation_is_transactional() {
    let mut gas = field(16, [0.0; 3]);
    let mut reference = field(16, [0.0; 3]);
    let mut piston = Piston::new(0.4, -0.04).unwrap();
    let mut other = piston;
    let load = PistonLoad {
        exterior_pressure_pa: 1.5,
        external_force_n: -0.1,
    };
    let zero = PistonSuspension::new(0.0, 0.0, -3.0).unwrap();
    let suspended = gas
        .step_with_suspended_piston(
            &mut piston,
            Time::new(0.004),
            StepControl::default(),
            load,
            zero,
            controls(),
            || false,
        )
        .unwrap();
    let original = reference
        .step_with_piston(
            &mut other,
            Time::new(0.004),
            StepControl::default(),
            load,
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(suspended.piston, original);
    assert_eq!(piston, other);
    assert_eq!(gas.states, reference.states);
    assert_eq!(gas.mesh(), reference.mesh());

    let spring = PistonSuspension::new(4.0, 0.3, 1.0).unwrap();
    let before = gas.states.clone();
    let mesh = gas.mesh().clone();
    let time = gas.time();
    let accepted = piston;
    assert_eq!(
        gas.step_with_suspended_piston(
            &mut piston,
            Time::new(0.004),
            StepControl::default(),
            load,
            spring,
            controls(),
            || true,
        ),
        Err(PistonError::Cancelled)
    );
    assert_eq!(gas.states, before);
    assert_eq!(gas.mesh(), &mesh);
    assert_eq!(gas.time(), time);
    assert_eq!(piston, accepted);
}
