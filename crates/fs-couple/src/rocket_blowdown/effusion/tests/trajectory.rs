//! Continuous kinetic-theory limits and accepted-only adaptive physical windows.
mod exchange;
use super::super::trajectory::{
    EffusionAdaptiveSolve, EffusionTrajectoryCandidate, EffusionTrajectoryError,
};
use super::*;
use fs_material::thermochemical_gas::energy::parcel::transfer::species::network::{
    SpeciesNetworkError, transport::GasTransportError,
};

fn quiet() -> EffusionRocketEnvironment {
    EffusionRocketEnvironment {
        inertial_frame: FRAME,
        gravity_world_m_s2: Vec3::ZERO,
        external_force_world_n: Vec3::ZERO,
    }
}
fn adaptive(scale: f64) -> EffusionAdaptiveSolve {
    EffusionAdaptiveSolve {
        initial_step_s: 1.0,
        minimum_step_s: 1e-5,
        maximum_step_s: 1.0,
        maximum_attempts: 1024,
        maximum_accepted_intervals: 512,
        relative_tolerance: 1e-3 * scale,
        temperature_tolerance_k: 0.05 * scale,
        pressure_tolerance_pa: 5.0 * scale,
        species_mass_tolerance_kg: 1e-10 * scale,
        position_tolerance_m: 1e-4 * scale,
        velocity_tolerance_m_s: 1e-3 * scale,
    }
}
fn evolve(
    before: &GasRocketState,
    outlets: &[EffusionOutlet<'_>],
    duration: f64,
    adaptive: EffusionAdaptiveSolve,
) -> EffusionTrajectoryCandidate {
    before
        .advance_effusion_adaptive(
            quiet(),
            outlets,
            duration,
            endpoint(),
            controls(),
            adaptive,
            || false,
        )
        .unwrap()
}
fn species_masses(state: &GasRocketState) -> Vec<f64> {
    let mixture = state.reservoir().model().thermodynamics();
    let molar: Vec<_> = mixture
        .components()
        .iter()
        .map(|row| row.molar_mass())
        .collect();
    mixture
        .mole_composition()
        .to_mass_fractions(&molar)
        .unwrap()
        .fractions()
        .iter()
        .map(|fraction| fraction * state.reservoir().mass_kg())
        .collect()
}

// The exact continuous monatomic cooling solution gives p(t) and m(t).
// Independently integrate its recoil acceleration -A*p/(2*(m_dry+m_gas));
// no discrete gas step, sampled exhaust velocity or controller is reused.
fn continuous_monatomic_velocity(before: &GasRocketState, duration: f64, lambda: f64) -> f64 {
    const PANELS: usize = 4096;
    let h = duration / PANELS as f64;
    let mut sum = 0.0;
    for i in 0..=PANELS {
        let z = 1.0 + lambda * (i as f64 * h) / 6.0;
        let pressure = P0 * z.powi(-8);
        let gas_mass = before.reservoir().mass_kg() * z.powi(-6);
        let acceleration = -0.5 * AREA * pressure / (before.vehicle().dry_mass_kg + gas_mass);
        let weight = if i == 0 || i == PANELS {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        sum += weight * acceleration;
    }
    before.vehicle().velocity_m_s.x + h * sum / 3.0
}

#[test]
fn tighter_local_controls_improve_exact_monatomic_cooling_and_continuous_recoil() {
    let before = initial(PURE, Vec3::ZERO);
    let law = interface(before.reservoir().model());
    let outlets = [outlet(&law, Vec3::new(1.0, 0.0, 0.0))];
    let duration = 4.0;
    let loose = evolve(&before, &outlets, duration, adaptive(1.0));
    let tight = evolve(&before, &outlets, duration, adaptive(0.01));
    let lambda = (AREA / VOLUME) * (R * T0 / (2.0 * core::f64::consts::PI * 0.028)).sqrt();
    let z = 1.0 + lambda * duration / 6.0;
    let exact_t = T0 * z.powi(-2);
    let exact_p = P0 * z.powi(-8);
    let exact_mass = before.reservoir().mass_kg() * z.powi(-6);
    let exact_v = continuous_monatomic_velocity(&before, duration, lambda);
    for (coarse, fine, exact) in [
        (
            loose.endpoint().reservoir().state().state().temperature,
            tight.endpoint().reservoir().state().state().temperature,
            exact_t,
        ),
        (
            loose.endpoint().reservoir().state().state().pressure,
            tight.endpoint().reservoir().state().state().pressure,
            exact_p,
        ),
        (
            loose.endpoint().reservoir().mass_kg(),
            tight.endpoint().reservoir().mass_kg(),
            exact_mass,
        ),
        (
            loose.endpoint().vehicle().velocity_m_s.x,
            tight.endpoint().vehicle().velocity_m_s.x,
            exact_v,
        ),
    ] {
        assert!(
            (fine - exact).abs() < 0.3 * (coarse - exact).abs(),
            "tight={fine:.17e}, loose={coarse:.17e}, continuous={exact:.17e}"
        );
    }
    assert!(tight.step().accepted_intervals > loose.step().accepted_intervals);
    assert!(tight.step().maximum_accepted_error <= 1.0);
    assert_eq!(tight.step().duration_s, duration);
    assert!(tight.endpoint().reservoir().state().state().temperature < T0);
}

#[test]
fn adaptive_mixture_resolves_fractionation_against_fine_original_source_steps() {
    let before = initial(MIXTURE, Vec3::ZERO);
    let law = interface(before.reservoir().model());
    let outlets = [outlet(&law, Vec3::new(1.0, 0.0, 0.0))];
    let mut reference = before.clone();
    for _ in 0..256 {
        reference = reference
            .advance_effusion(FRAME, &outlets, 1.0 / 256.0, endpoint(), controls(), || {
                false
            })
            .unwrap()
            .into_parts()
            .0;
    }
    let loose = evolve(&before, &outlets, 1.0, adaptive(1.0));
    let tight = evolve(&before, &outlets, 1.0, adaptive(0.01));
    let reference_masses = species_masses(&reference);
    let discrepancy = |state: &GasRocketState| {
        let actual = state.reservoir().state().state();
        let target = reference.reservoir().state().state();
        (actual.temperature - target.temperature).abs() / T0
            + (actual.pressure - target.pressure).abs() / P0
            + species_masses(state)
                .iter()
                .zip(&reference_masses)
                .map(|(a, b)| (a - b).abs())
                .sum::<f64>()
                / reference.reservoir().mass_kg()
            + state
                .vehicle()
                .velocity_m_s
                .sub(reference.vehicle().velocity_m_s)
                .norm_squared()
                .sqrt()
    };
    assert!(discrepancy(tight.endpoint()) < 0.3 * discrepancy(loose.endpoint()));
    let fractions = tight
        .endpoint()
        .reservoir()
        .model()
        .thermodynamics()
        .mole_composition()
        .fractions();
    let reference_fractions = reference
        .reservoir()
        .model()
        .thermodynamics()
        .mole_composition()
        .fractions();
    assert!(fractions[1] < 0.5); // The original light B species escapes faster.
    near(fractions[1], reference_fractions[1], 1e-5);
    near(
        tight.endpoint().reservoir().state().state().temperature,
        reference.reservoir().state().state().temperature,
        0.005,
    );
}

#[test]
fn accepted_physical_totals_close_the_window_without_an_unrepresentable_final_tail() {
    let before = initial(MIXTURE, Vec3::new(2.0, -1.0, 0.5));
    let saved = before.clone();
    let law = interface(before.reservoir().model());
    // Eight naive additions of 0.0125 leave a tiny remainder below the physical
    // gas owner's representable withdrawal. The controller must split its last
    // pair and reach 0.1, rather than discard time or relax a physical budget.
    let schedule = EffusionAdaptiveSolve {
        initial_step_s: 0.0125,
        maximum_step_s: 0.0125,
        minimum_step_s: 1e-6,
        ..adaptive(1000.0)
    };
    let candidate = before
        .advance_effusion_adaptive(
            loaded_environment(),
            &[outlet(&law, Vec3::new(0.6, 0.8, 0.0))],
            0.1,
            endpoint(),
            controls(),
            schedule,
            || false,
        )
        .unwrap();
    let step = candidate.step();
    assert_eq!(step.duration_s, 0.1);
    assert_eq!(step.rejected_intervals, 0);
    assert_eq!(step.attempted_intervals, step.accepted_intervals);
    assert_eq!(step.accepted_steps.len(), 2 * step.accepted_intervals);
    let (mut time, mut escaped, mut energy, mut work) = (0.0, 0.0, 0.0, 0.0);
    let (mut momentum, mut impulse) = (Vec3::ZERO, Vec3::ZERO);
    let (mut mass_residual, mut momentum_residual, mut energy_residual) = (0.0, Vec3::ZERO, 0.0);
    let (mut gravity_impulse, mut gravity_work) = (Vec3::ZERO, 0.0);
    for accepted in &step.accepted_steps {
        time += accepted.duration_s;
        escaped += accepted.expelled_mass_kg;
        momentum = momentum.add(accepted.exhaust_momentum_world_kg_m_s);
        energy += accepted.exhaust_energy_world_j;
        impulse = impulse.add(accepted.exhaust.external_impulse_n_s);
        work += accepted.exhaust.external_work_j;
        mass_residual += accepted.gas.exchange.total_mass_residual_kg;
        momentum_residual = momentum_residual.add(accepted.momentum_residual_kg_m_s);
        energy_residual += accepted.energy_residual_j;
        gravity_impulse = gravity_impulse.add(accepted.gravity_impulse_uncertainty_n_s);
        gravity_work += accepted.gravity_work_uncertainty_j;
    }
    near(time, 0.1, 1e-15);
    near(step.total_escaped_mass_kg, escaped, 1e-18);
    vector_near(step.total_escaped_momentum_world_kg_m_s, momentum, 2e-15);
    near(step.total_escaped_energy_world_j, energy, 2e-12);
    vector_near(step.total_external_impulse_n_s, impulse, 2e-15);
    near(step.total_external_work_j, work, 2e-12);
    near(step.total_mass_residual_kg, mass_residual, 1e-25);
    vector_near(
        step.total_momentum_residual_kg_m_s,
        momentum_residual,
        1e-25,
    );
    near(step.total_energy_residual_j, energy_residual, 1e-18);
    vector_near(
        step.total_gravity_impulse_uncertainty_n_s,
        gravity_impulse,
        1e-25,
    );
    near(step.total_gravity_work_uncertainty_j, gravity_work, 1e-25);
    let after = candidate.endpoint();
    near(
        before.reservoir().mass_kg() - after.reservoir().mass_kg(),
        escaped,
        1e-13,
    );
    let physical_momentum = after
        .vehicle()
        .velocity_m_s
        .scale(mass(after))
        .sub(before.vehicle().velocity_m_s.scale(mass(&before)))
        .add(momentum)
        .sub(impulse);
    vector_near(
        physical_momentum,
        step.total_momentum_residual_kg_m_s,
        2e-12,
    );
    near(
        internal(after) - internal(&before) + kinetic(after) - kinetic(&before) + energy - work,
        step.total_energy_residual_j,
        2e-8,
    );
    assert_eq!(before, saved);
}

#[test]
fn adaptive_closed_loaded_coast_retains_gas_and_matches_constant_mass_ballistics() {
    let before = initial(PURE, Vec3::new(3.0, 2.0, -1.0));
    let foreign = initial("A 1 0.028 3.5 0 1.7e-5 0.025", Vec3::ZERO);
    let law = interface(foreign.reservoir().model());
    let closed = EffusionOutlet {
        area_m2: 0.0,
        ..outlet(&law, Vec3::new(1.0, 0.0, 0.0))
    };
    let environment = loaded_environment();
    let candidate = before
        .advance_effusion_adaptive(
            environment,
            &[closed],
            1.0,
            endpoint(),
            controls(),
            EffusionAdaptiveSolve {
                initial_step_s: 0.125,
                maximum_step_s: 0.5,
                ..adaptive(1e-4)
            },
            || false,
        )
        .unwrap();
    let total_force = environment
        .external_force_world_n
        .add(environment.gravity_world_m_s2.scale(mass(&before)));
    let acceleration = total_force.scale(1.0 / mass(&before));
    let displacement = before.vehicle().velocity_m_s.add(acceleration.scale(0.5));
    vector_near(
        candidate.endpoint().vehicle().velocity_m_s,
        before.vehicle().velocity_m_s.add(acceleration),
        2e-12,
    );
    vector_near(
        candidate.endpoint().vehicle().position_m,
        before.vehicle().position_m.add(displacement),
        2e-12,
    );
    assert_eq!(candidate.endpoint().reservoir(), before.reservoir());
    assert_eq!(candidate.step().total_escaped_mass_kg, 0.0);
    assert_eq!(
        candidate.step().total_escaped_momentum_world_kg_m_s,
        Vec3::ZERO
    );
    assert_eq!(candidate.step().total_escaped_energy_world_j, 0.0);
    vector_near(
        candidate.step().total_external_impulse_n_s,
        total_force,
        2e-14,
    );
    near(
        candidate.step().total_external_work_j,
        total_force.dot(displacement),
        2e-12,
    );
}

#[test]
fn only_typed_withdrawal_retries_and_minimum_or_work_limits_publish_no_endpoint() {
    let before = initial(PURE, Vec3::ZERO);
    let saved = before.clone();
    let law = interface(before.reservoir().model());
    let outlets = [outlet(&law, Vec3::new(1.0, 0.0, 0.0))];
    let mut physical = controls();
    physical.transport.exchange.maximum_species_fraction = 0.001;
    assert!(matches!(
        before.advance_effusion(FRAME, &outlets, 1.0, endpoint(), physical, || false),
        Err(EffusionRocketError::Gas(GasTransportError::Exchange(
            SpeciesNetworkError::Withdrawal { .. }
        )))
    ));
    let candidate = before
        .advance_effusion_adaptive(
            quiet(),
            &outlets,
            1.0,
            endpoint(),
            physical,
            adaptive(1000.0),
            || false,
        )
        .unwrap();
    assert!(candidate.step().rejected_intervals > 0);
    assert_eq!(
        candidate.step().attempted_intervals,
        candidate.step().accepted_intervals + candidate.step().rejected_intervals
    );
    assert!(candidate.step().accepted_steps[0].duration_s < 0.05);
    assert!(matches!(
        before.advance_effusion_adaptive(
            quiet(),
            &outlets,
            1.0,
            endpoint(),
            physical,
            EffusionAdaptiveSolve {
                minimum_step_s: 0.125,
                ..adaptive(1000.0)
            },
            || false,
        ),
        Err(EffusionTrajectoryError::Step {
            error: EffusionRocketError::Gas(GasTransportError::Exchange(
                SpeciesNetworkError::Withdrawal { .. }
            )),
            ..
        })
    ));
    assert!(matches!(before.advance_effusion_adaptive(
        quiet(), &outlets, 1.0, endpoint(), controls(),
        EffusionAdaptiveSolve { minimum_step_s: 1.0, ..adaptive(1e-10) }, || false,
    ), Err(EffusionTrajectoryError::MinimumStep { normalized_error, .. }) if normalized_error > 1.0));
    assert!(matches!(
        before.advance_effusion_adaptive(
            quiet(),
            &outlets,
            1.0,
            endpoint(),
            physical,
            EffusionAdaptiveSolve {
                maximum_attempts: 1,
                ..adaptive(1000.0)
            },
            || false,
        ),
        Err(EffusionTrajectoryError::WorkBudget {
            attempted_intervals: 1,
            ..
        })
    ));
    assert!(matches!(
        before.advance_effusion_adaptive(
            quiet(),
            &[],
            0.5,
            endpoint(),
            controls(),
            EffusionAdaptiveSolve {
                initial_step_s: 0.25,
                maximum_step_s: 0.25,
                maximum_accepted_intervals: 1,
                ..adaptive(1.0)
            },
            || false,
        ),
        Err(EffusionTrajectoryError::WorkBudget {
            accepted_intervals: 1,
            ..
        })
    ));
    assert_eq!(before, saved);
}

#[test]
fn foreign_sources_frames_and_late_cancellation_are_atomic_and_never_retried_as_error() {
    let before = initial(PURE, Vec3::ZERO);
    let saved = before.clone();
    let foreign = initial("A 1 0.028 3.5 0 1.7e-5 0.025", Vec3::ZERO);
    let other_law = interface(foreign.reservoir().model());
    let error = before
        .advance_effusion_adaptive(
            quiet(),
            &[outlet(&other_law, Vec3::new(1.0, 0.0, 0.0))],
            1.0,
            endpoint(),
            controls(),
            adaptive(1.0),
            || false,
        )
        .unwrap_err();
    assert!(matches!(error, EffusionTrajectoryError::Step {
        elapsed_s, duration_s, ..
    } if elapsed_s == 0.0 && duration_s == 1.0));
    assert!(matches!(
        before.advance_effusion_adaptive(
            EffusionRocketEnvironment {
                inertial_frame: FRAME + 1,
                ..quiet()
            },
            &[],
            1.0,
            endpoint(),
            controls(),
            adaptive(1.0),
            || false,
        ),
        Err(EffusionTrajectoryError::Step {
            error: EffusionRocketError::Invalid(_),
            ..
        })
    ));
    let law = interface(before.reservoir().model());
    let outlets = [outlet(&law, Vec3::new(1.0, 0.0, 0.0))];
    let schedule = EffusionAdaptiveSolve {
        initial_step_s: 0.05,
        maximum_step_s: 0.05,
        ..adaptive(1.0)
    };
    let mut polls = 0;
    before
        .advance_effusion_adaptive(
            quiet(),
            &outlets,
            0.1,
            endpoint(),
            controls(),
            schedule,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    for stop in [1, polls - 1, polls] {
        let mut count = 0;
        assert!(
            before
                .advance_effusion_adaptive(
                    quiet(),
                    &outlets,
                    0.1,
                    endpoint(),
                    controls(),
                    schedule,
                    || {
                        count += 1;
                        count == stop
                    },
                )
                .is_err()
        );
        assert_eq!(before, saved);
    }
}
