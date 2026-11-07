//! Independent scalar kinetics and P1 field balances exercise the real native
//! time owner, including its separate resolved-conduction and thermal accounts.

use super::*;

fn rate_problem(
    depth: f64,
    sigma: f64,
    viscosity: Option<f64>,
) -> (PlanarHystereticMagnetostatics, Vec<EddyConductor>) {
    let (magnetic, conductor) = sources_with_viscosity(MU, true, sigma, viscosity);
    problem_from_sources(paired_squares(depth), magnetic, conductor)
}

fn balanced_current(drive_h: f64) -> [f64; 8] {
    let mut currents = [6.0 * drive_h; 8];
    currents[4..].fill(-6.0 * drive_h);
    currents
}

// One-dimensional backward-Euler minimization of the ORIGINAL quadratic
// branch storage plus dry friction and eta*(z-z_old)^2/(2*dt). This closed
// form never invokes the production vector constitutive solver or tangent.
fn kinetic_play(old_z: f64, h: f64, eta: f64, dt: f64) -> f64 {
    let stiffness = CHI.recip() + eta / dt;
    let drive = h + eta / dt * old_z;
    old_z
        .clamp((drive - PINNING) / stiffness, (drive + PINNING) / stiffness)
        .clamp(-SATURATION, SATURATION)
}

fn kinetic_relaxation(old_z: f64, old_b: f64, eta: f64, sigma: f64, dt: f64) -> (f64, f64, f64) {
    // Opposite center potentials give zero mean dA/dt in the connected
    // conductor. Its exact projected P1 balance is H + sigma*delta(B)/(24dt)=0.
    let mut lo = -30.0;
    let mut hi = 30.0;
    for _ in 0..90 {
        let h = 0.5 * (lo + hi);
        let b = MU * h + kinetic_play(old_z, h, eta, dt);
        if h + sigma / (24.0 * dt) * (b - old_b) < 0.0 {
            lo = h;
        } else {
            hi = h;
        }
    }
    let h = 0.5 * (lo + hi);
    let z = kinetic_play(old_z, h, eta, dt);
    (h, z, MU * h + z)
}

#[test]
fn intrinsic_switching_and_resolved_eddy_diffusion_close_independent_physical_accounts() {
    let depth = 0.3;
    let eta = 0.6;
    for (sigma, gauge, sign) in [(18.0, 0.0, 1.0), (36.0, 100.0, -1.0)] {
        let (owner, conductors) = rate_problem(depth, sigma, Some(eta));
        let initial = owner.demagnetized(gauge, &mut || false).unwrap();
        let mut state = PlanarHystereticEddy::new(
            &initial,
            conductors,
            &[],
            &boundary(6, gauge),
            0.0,
            controls(),
            &mut || false,
        )
        .unwrap();
        let (mut old_h, mut old_z, mut old_b) = (0.0, 0.0, 0.0);
        let mut total_viscous = 0.0;
        let mut total_pinning = 0.0;
        let mut total_joule = 0.0;
        for (index, dt) in [1.0, 0.001, 0.01, 0.1, 1.0, 10.0, 100.0]
            .into_iter()
            .enumerate()
        {
            // Prepare a known saturated state by a balanced physical source,
            // then release that source without replacing its native memory.
            let (h, z, b, drive_h) = if index == 0 {
                let h = sign * 12.0;
                let z = kinetic_play(old_z, h, eta, dt);
                let b = MU * h + z;
                (h, z, b, h + sigma / (24.0 * dt) * (b - old_b))
            } else {
                let (h, z, b) = kinetic_relaxation(old_z, old_b, eta, sigma, dt);
                (h, z, b, 0.0)
            };
            let next = state
                .step_to_candidate(
                    state.time_s() + dt,
                    state.conductors().to_vec(),
                    &balanced_current(drive_h),
                    controls(),
                    &mut || false,
                )
                .unwrap();
            let step = next.step();
            let stored_before = 2.0 * depth * storage(old_h, old_z);
            let stored_after = 2.0 * depth * storage(h, z);
            let pinning = 2.0 * depth * PINNING * (z - old_z).abs();
            let viscous = 2.0 * depth * eta * (z - old_z).powi(2) / dt;
            let joule = sigma * depth * (b - old_b).powi(2) / (12.0 * dt);
            let field_work = 2.0 * depth * h * (b - old_b);
            let numerical = field_work - (stored_after - stored_before) - pinning - viscous;
            near(next.state().potential_wb_m()[6] - gauge, 0.5 * b);
            near(next.state().potential_wb_m()[7] - gauge, -0.5 * b);
            near(step.stored_before_j, stored_before);
            near(step.stored_after_j, stored_after);
            near(step.magnetic.pinning_heat_j, pinning);
            near(step.magnetic.viscous_heat_j, viscous);
            near(step.magnetic.irreversible_heat_j, pinning + viscous);
            near(step.hysteresis_heat_j, pinning + viscous);
            near(step.joule_heat_j, joule);
            near(step.heat_j, pinning + viscous + joule);
            near(step.numerical_dissipation_j, numerical);
            near(step.supplied_work_j, 2.0 * depth * drive_h * (b - old_b));
            near(step.magnetic.field_work_j, field_work);
            near(step.energy_residual_j, 0.0);
            assert!(step.numerical_dissipation_j >= 0.0);
            for cell in 0..8 {
                let direction = if cell < 4 { 1.0 } else { -1.0 };
                let memory = &next.state().magnetic_state().cell_states()[cell];
                for axis in 0..2 {
                    near(
                        memory.field_a_m()[axis],
                        direction * h * DIRECTIONS[cell % 4][axis],
                    );
                    near(
                        memory.polarizations_t()[0][axis],
                        direction * z * DIRECTIONS[cell % 4][axis],
                    );
                    near(
                        next.state().magnetic_state().magnetic_flux_density_t()[cell][axis],
                        direction * b * DIRECTIONS[cell % 4][axis],
                    );
                }
                near(step.magnetic.cells[cell].pinning_heat_j, pinning / 8.0);
                near(step.magnetic.cells[cell].viscous_heat_j, viscous / 8.0);
                near(
                    step.magnetic.cells[cell].irreversible_heat_j,
                    (pinning + viscous) / 8.0,
                );
                near(step.cells[cell].heat_j, joule / 8.0);
            }
            total_pinning += pinning;
            total_viscous += viscous;
            total_joule += joule;
            (state, _) = next.into_parts();
            (old_h, old_z, old_b) = (h, z, b);
        }
        assert!(total_pinning > 0.0 && total_viscous > 0.0 && total_joule > 0.0);
        assert!(sign * old_z > 0.39 && old_h.abs() < 1e-4);
    }
}

#[test]
fn intrinsic_switching_rate_changes_the_field_without_any_resolved_eddy_conductor() {
    let eta = 0.6;
    let depth = 0.3;
    let (owner, _) = rate_problem(depth, 18.0, Some(eta));
    let initial = owner.demagnetized(0.0, &mut || false).unwrap();
    let state = PlanarHystereticEddy::new(
        &initial,
        vec![],
        &(0..8).collect::<Vec<_>>(),
        &boundary(6, 0.0),
        0.0,
        controls(),
        &mut || false,
    )
    .unwrap();
    let mut switched = Vec::new();
    for dt in [0.01, 1.0] {
        let h = 6.0;
        let z = kinetic_play(0.0, h, eta, dt);
        let next = state
            .step_to_candidate(dt, vec![], &balanced_current(h), controls(), &mut || false)
            .unwrap();
        near(next.state().potential_wb_m()[6], 0.5 * (MU * h + z));
        near(
            next.step().magnetic.pinning_heat_j,
            2.0 * depth * PINNING * z,
        );
        near(
            next.step().magnetic.viscous_heat_j,
            2.0 * depth * eta * z * z / dt,
        );
        assert_eq!(next.step().joule_heat_j, 0.0);
        near(next.step().heat_j, next.step().hysteresis_heat_j);
        near(next.step().energy_residual_j, 0.0);
        switched.push(next.state().magnetic_state().cell_states()[0].polarizations_t()[0][0]);
    }
    assert!(switched[1] > 10.0 * switched[0]);
    assert_eq!(state.time_s(), 0.0);
    assert_eq!(state.magnetic_state().cell_states(), initial.cell_states());
}

#[test]
fn explicitly_zero_intrinsic_viscosity_recovers_existing_eddy_history_bitwise() {
    let (legacy, legacy_conductors) = rate_problem(0.3, 18.0, None);
    let (zero, zero_conductors) = rate_problem(0.3, 18.0, Some(0.0));
    let legacy_initial = legacy.demagnetized(0.0, &mut || false).unwrap();
    let zero_initial = zero.demagnetized(0.0, &mut || false).unwrap();
    let mut legacy_state = PlanarHystereticEddy::new(
        &legacy_initial,
        legacy_conductors,
        &[],
        &boundary(6, 0.0),
        0.0,
        controls(),
        &mut || false,
    )
    .unwrap();
    let mut zero_state = PlanarHystereticEddy::new(
        &zero_initial,
        zero_conductors,
        &[],
        &boundary(6, 0.0),
        0.0,
        controls(),
        &mut || false,
    )
    .unwrap();
    for (target, drive) in [(0.1, 15.0), (0.4, 0.0), (0.8, -12.0), (2.0, 4.0)] {
        let baseline = legacy_state
            .step_to_candidate(
                target,
                legacy_state.conductors().to_vec(),
                &balanced_current(drive),
                controls(),
                &mut || false,
            )
            .unwrap();
        let with_zero = zero_state
            .step_to_candidate(
                target,
                zero_state.conductors().to_vec(),
                &balanced_current(drive),
                controls(),
                &mut || false,
            )
            .unwrap();
        assert_eq!(
            baseline.state().potential_wb_m(),
            with_zero.state().potential_wb_m()
        );
        assert_eq!(
            baseline.state().energy_j().to_bits(),
            with_zero.state().energy_j().to_bits()
        );
        for (a, b) in baseline
            .state()
            .magnetic_state()
            .cell_states()
            .iter()
            .zip(with_zero.state().magnetic_state().cell_states())
        {
            assert_eq!(a.field_a_m(), b.field_a_m());
            assert_eq!(a.flux_density_t(), b.flux_density_t());
            assert_eq!(a.polarizations_t(), b.polarizations_t());
        }
        let a = baseline.step();
        let b = with_zero.step();
        for (a, b) in [
            (a.magnetic.pinning_heat_j, b.magnetic.pinning_heat_j),
            (a.hysteresis_heat_j, b.hysteresis_heat_j),
            (a.joule_heat_j, b.joule_heat_j),
            (a.heat_j, b.heat_j),
            (a.supplied_work_j, b.supplied_work_j),
            (a.numerical_dissipation_j, b.numerical_dissipation_j),
        ] {
            assert_eq!(a.to_bits(), b.to_bits());
        }
        assert_eq!(b.magnetic.viscous_heat_j, 0.0);
        (legacy_state, _) = baseline.into_parts();
        (zero_state, _) = with_zero.into_parts();
    }
}

#[test]
fn clockless_failed_and_cancelled_kinetic_trials_leave_history_available_for_retry() {
    let (owner, conductors) = rate_problem(0.3, 18.0, Some(0.6));
    let initial = owner.demagnetized(0.0, &mut || false).unwrap();
    let current = balanced_current(12.0);
    assert!(
        owner
            .candidate(
                &initial,
                &current,
                &boundary(6, 0.0),
                controls().magnetic,
                &mut || false,
            )
            .is_err()
    );
    let state = PlanarHystereticEddy::new(
        &initial,
        conductors,
        &[],
        &boundary(6, 0.0),
        0.0,
        controls(),
        &mut || false,
    )
    .unwrap();
    let mut limited = controls();
    limited.magnetic.nonlinear.max_newton_iterations = 0;
    assert!(
        state
            .step_to_candidate(
                0.1,
                state.conductors().to_vec(),
                &current,
                limited,
                &mut || false,
            )
            .is_err()
    );
    let mut polls = 0;
    let cancelled = state.step_to_candidate(
        0.1,
        state.conductors().to_vec(),
        &current,
        controls(),
        &mut || {
            polls += 1;
            polls >= 24
        },
    );
    assert!(matches!(cancelled, Err(HystereticEddyError::Cancelled)));
    assert_eq!(state.time_s(), 0.0);
    assert_eq!(state.potential_wb_m(), initial.potential_wb_m());
    assert_eq!(state.magnetic_state().cell_states(), initial.cell_states());
    let retry = state
        .step_to_candidate(
            0.1,
            state.conductors().to_vec(),
            &current,
            controls(),
            &mut || false,
        )
        .unwrap();
    assert!(retry.step().magnetic.viscous_heat_j > 0.0);
    assert_eq!(state.magnetic_state().cell_states(), initial.cell_states());
}
