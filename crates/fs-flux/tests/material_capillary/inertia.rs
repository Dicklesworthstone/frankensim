//! Actual original material/interface integration for inertial capillary motion.
use super::contact_line::{limits, sourced_wetting};
use super::*;
use fs_flux::capillary::ContactLineMotion;
use fs_flux::capillary::inertia::{InertialCapillaryError, InertialCapillarySolve};
use fs_flux::capillary::material::inertia::advance::{
    CapillaryInertiaAdvanceControls, CapillaryInertiaAdvanceError,
};
use fs_flux::capillary::material::inertia::{
    CapillaryInertiaParameters, InertialCapillaryBoundary, InertialCapillaryIntervalControls,
    InertialMeniscusState, MaterialInertiaError,
};

fn parameters() -> CapillaryInertiaParameters {
    CapillaryInertiaParameters {
        radius_m: 0.001,
        tube_length_m: 0.1,
        liquid_added_length_m: 0.0,
        displaced_added_length_m: 0.0,
        liquid_entry_loss_coefficient: 0.0,
        liquid_exit_loss_coefficient: 0.0,
        displaced_entry_loss_coefficient: 0.0,
        displaced_exit_loss_coefficient: 0.0,
        maximum_velocity_m_s: 0.01,
    }
}
fn state(velocity: f64) -> InertialMeniscusState {
    InertialMeniscusState {
        wet_length_m: 0.02,
        velocity_m_s: velocity,
    }
}
fn controls() -> InertialCapillaryIntervalControls {
    InertialCapillaryIntervalControls {
        mechanics: InertialCapillarySolve {
            maximum_iterations: 80,
            pressure_absolute_tolerance_pa: 1e-8,
            velocity_absolute_tolerance_m_s: 1e-12,
            time_absolute_tolerance_s: 1e-12,
            energy_absolute_tolerance_j: 1e-18,
            energy_relative_tolerance: 1e-9,
        },
        maximum_stopping_events: 16,
    }
}
fn adaptive() -> CapillaryInertiaAdvanceControls {
    CapillaryInertiaAdvanceControls {
        initial_step_s: 0.02,
        minimum_step_s: 1e-10,
        maximum_step_s: 0.02,
        maximum_attempts: 4096,
        relative_tolerance: 1e-6,
        length_absolute_tolerance_m: 1e-11,
        velocity_absolute_tolerance_m_s: 1e-11,
        event_time_absolute_tolerance_s: 1e-8,
    }
}
fn near(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "{actual} != {expected}"
    );
}
fn displaced(rho: f64, mu: f64) -> DisplacedFixture {
    DisplacedFixture {
        label: "original frozen displaced-fluid fixture",
        transport: FluidTransportState::try_new(rho, 1000.0, 0.5, mu).unwrap(),
        reads: Cell::new(0),
        refuse: false,
    }
}
fn balanced(duration: f64) -> InertialCapillaryBoundary {
    InertialCapillaryBoundary {
        applied_pressure_pa: -2.0 * 0.05 * 1.0_f64.cos() / parameters().radius_m,
        opposing_gravity_m_s2: 0.0,
        duration_s: duration,
    }
}

#[test]
fn source_owned_added_inertance_gives_finite_dry_startup_without_regularization() {
    let (liquid, wetting) = fixture(1000.0, 0.05, 0.05, 1.0);
    let tube = MaterialCapillaryTube::new(&liquid, &wetting).unwrap();
    let dry = InertialMeniscusState {
        wet_length_m: 0.0,
        velocity_m_s: 0.0,
    };
    assert!(matches!(
        tube.with_inertia(parameters(), dry),
        Err(MaterialInertiaError::Invalid(_))
    ));
    let parameters = CapillaryInertiaParameters {
        liquid_added_length_m: 0.002,
        liquid_entry_loss_coefficient: 1.0,
        liquid_exit_loss_coefficient: 1.0,
        ..parameters()
    };
    let original = tube.with_inertia(parameters, dry).unwrap();
    let dt = 1e-7;
    let candidate = original
        .step_candidate(
            InertialCapillaryBoundary {
                applied_pressure_pa: 0.0,
                opposing_gravity_m_s2: 0.0,
                duration_s: dt,
            },
            controls(),
            || false,
        )
        .unwrap();
    let acceleration = 2.0 * 0.05 * 1.0_f64.cos()
        / (parameters.radius_m * 1000.0 * parameters.liquid_added_length_m);
    near(
        candidate.step().state.velocity_m_s / dt,
        acceleration,
        0.0,
        1e-4,
    );
    near(
        candidate.step().state.wet_length_m / dt.powi(2),
        0.5 * acceleration,
        0.0,
        1e-4,
    );
    assert!(candidate.step().kinetic_energy_j > 0.0);
    assert!(candidate.step().reservoir_dissipation_j > 0.0);
    assert!(candidate.step().energy_balance_residual_j.abs() < 1e-18);
    assert_eq!(original.state(), dry);
    assert!(std::ptr::eq(candidate.endpoint().tube().liquid(), &liquid));
    assert!(std::ptr::eq(
        candidate.endpoint().tube().wetting(),
        &wetting
    ));
}

#[test]
fn both_original_fluids_determine_inertance_and_directional_reservoir_losses_once() {
    let (liquid, wetting) = fixture(1000.0, 0.05, 0.05, 1.0);
    let displaced = displaced(800.0, 0.02);
    let p = CapillaryInertiaParameters {
        liquid_added_length_m: 0.002,
        displaced_added_length_m: 0.003,
        liquid_entry_loss_coefficient: 0.6,
        liquid_exit_loss_coefficient: 1.1,
        displaced_entry_loss_coefficient: 0.7,
        displaced_exit_loss_coefficient: 1.4,
        ..parameters()
    };
    let area = PI * p.radius_m.powi(2);
    for velocity in [0.0004, -0.0004] {
        let owner = MaterialCapillaryTube::new(&liquid, &wetting)
            .unwrap()
            .with_inertia_two_fluid(p, state(velocity), &displaced)
            .unwrap();
        let reads = displaced.reads.get();
        let input = InertialCapillaryBoundary {
            applied_pressure_pa: 0.0,
            opposing_gravity_m_s2: 2.0,
            duration_s: 1e-4,
        };
        let candidate = owner.step_candidate(input, controls(), || false).unwrap();
        let step = candidate.step();
        let u = 0.5 * (velocity + step.state.velocity_m_s);
        let x = 0.5 * (state(velocity).wet_length_m + step.state.wet_length_m);
        let coefficient = if velocity > 0.0 {
            0.5 * area * (1000.0 * 0.6 + 800.0 * 1.4)
        } else {
            0.5 * area * (1000.0 * 1.1 + 800.0 * 0.7)
        };
        near(
            step.liquid_viscous_dissipation_j,
            8.0 * PI * 0.05 * x * u * u * input.duration_s,
            1e-22,
            1e-8,
        );
        near(
            step.displaced_viscous_dissipation_j,
            8.0 * PI * 0.02 * (p.tube_length_m - x) * u * u * input.duration_s,
            1e-22,
            1e-8,
        );
        near(
            step.reservoir_dissipation_j,
            coefficient * u.abs().powi(3) * input.duration_s,
            1e-24,
            1e-8,
        );
        near(
            step.displaced_mass_kg,
            800.0 * step.transferred_volume_m3,
            1e-24,
            1e-12,
        );
        near(
            step.transferred_liquid_mass_kg,
            1000.0 * step.transferred_volume_m3,
            1e-24,
            1e-12,
        );
        near(
            step.outlet_potential_energy_j,
            800.0 * 2.0 * p.tube_length_m * step.transferred_volume_m3,
            1e-22,
            1e-10,
        );
        assert!(step.energy_balance_residual_j.abs() < 1e-18);
        assert_eq!(displaced.reads.get(), reads);
        assert!(std::ptr::eq(
            candidate.endpoint().displaced_fluid().unwrap(),
            &displaced as &dyn FluidTransport
        ));
        assert_eq!(owner.state(), state(velocity));
    }
    assert_eq!(displaced.reads.get(), 2);
}

#[test]
fn a_moving_contact_line_coasts_to_a_resolved_stop_then_pins_without_zeno() {
    let (liquid, _) = fixture(1000.0, 0.05, 0.05, 1.0);
    let wetting = sourced_wetting(&liquid, 1.0, 2.0, 0.8, 1.2);
    let displaced = displaced(1000.0, 0.05);
    let p = parameters();
    let owner = MaterialCapillaryTube::new(&liquid, &wetting)
        .unwrap()
        .with_contact_line(limits())
        .unwrap()
        .with_inertia_two_fluid(p, state(0.001), &displaced)
        .unwrap();
    let candidate = owner
        .advance_candidate(balanced(0.02), controls(), adaptive(), || false)
        .unwrap();
    let step = candidate.step();
    let area = PI * p.radius_m.powi(2);
    let mass = area * 1000.0 * p.tube_length_m;
    let bulk = 8.0 * PI * 0.05 * p.tube_length_m;
    let line = 2.0 * PI * p.radius_m * 2.0;
    let resistance = bulk + line;
    let threshold_force = area * (2.0 * 0.05 / p.radius_m) * (1.0_f64.cos() - 1.2_f64.cos());
    let a = threshold_force / resistance;
    let tau = mass / resistance;
    let stop_time = tau * (0.001 / a).ln_1p();
    let displacement = tau * (0.001 - a * (0.001 / a).ln_1p());
    let pinning_loss = threshold_force * displacement;
    let kinetic0 = 0.5 * mass * 0.001_f64.powi(2);
    near(
        step.state.wet_length_m - state(0.001).wet_length_m,
        displacement,
        2e-10,
        0.0,
    );
    near(step.first_stopping_time_s.unwrap(), stop_time, 2e-7, 0.0);
    near(step.pinning_dissipation_j, pinning_loss, 2e-14, 0.0);
    near(
        step.contact_line_viscous_dissipation_j,
        line / resistance * (kinetic0 - pinning_loss),
        2e-14,
        0.0,
    );
    near(
        step.liquid_viscous_dissipation_j + step.displaced_viscous_dissipation_j,
        bulk / resistance * (kinetic0 - pinning_loss),
        2e-14,
        0.0,
    );
    near(step.total_dissipation_j, kinetic0, 1e-18, 0.0);
    assert_eq!(step.elapsed_time_s, 0.02);
    assert_eq!(step.state.velocity_m_s, 0.0);
    assert_eq!(step.final_motion, ContactLineMotion::Pinned);
    assert_eq!(step.stopping_events, 1);
    assert!(candidate.report().rejected_intervals > 0);
    assert!(candidate.report().maximum_normalized_local_error <= 1.0);
    assert_eq!(owner.state(), state(0.001));
    assert_eq!(displaced.reads.get(), 1);
}

#[test]
fn viscous_coasting_has_no_fabricated_finite_stop_or_reflected_velocity() {
    let (liquid, wetting) = fixture(1000.0, 0.05, 0.05, 1.0);
    let displaced = displaced(1000.0, 0.05);
    let p = parameters();
    let owner = MaterialCapillaryTube::new(&liquid, &wetting)
        .unwrap()
        .with_inertia_two_fluid(p, state(0.001), &displaced)
        .unwrap();
    let candidate = owner
        .advance_candidate(balanced(0.02), controls(), adaptive(), || false)
        .unwrap();
    let tau = 1000.0 * p.radius_m.powi(2) / (8.0 * 0.05);
    let decay = (-0.02 / tau).exp();
    near(
        candidate.step().state.velocity_m_s,
        0.001 * decay,
        1e-9,
        0.0,
    );
    near(
        candidate.step().state.wet_length_m,
        state(0.001).wet_length_m + tau * 0.001 * (1.0 - decay),
        2e-10,
        0.0,
    );
    assert!(candidate.step().state.velocity_m_s > 0.0);
    assert_eq!(candidate.step().stopping_events, 0);
    assert_eq!(candidate.step().first_stopping_time_s, None);
    assert!(candidate.step().energy_balance_residual_j.abs() < 1e-18);
}

#[test]
fn cancellation_work_and_invalid_end_models_refuse_without_publishing_state() {
    let (liquid, _) = fixture(1000.0, 0.05, 0.05, 1.0);
    let wetting = sourced_wetting(&liquid, 1.0, 2.0, 0.8, 1.2);
    let tube = MaterialCapillaryTube::new(&liquid, &wetting)
        .unwrap()
        .with_contact_line(limits())
        .unwrap();
    assert!(matches!(
        tube.with_inertia(
            CapillaryInertiaParameters {
                displaced_added_length_m: 0.01,
                ..parameters()
            },
            state(0.0)
        ),
        Err(MaterialInertiaError::Material(
            MaterialCapillaryError::Unsupported(_)
        ))
    ));
    assert!(matches!(
        tube.with_inertia(parameters(), state(1.0)),
        Err(MaterialInertiaError::Invalid(_))
    ));
    let displaced = displaced(1000.0, 0.05);
    let owner = tube
        .with_inertia_two_fluid(parameters(), state(0.001), &displaced)
        .unwrap();
    assert!(matches!(
        owner.step_candidate(balanced(0.02), controls(), || true),
        Err(MaterialInertiaError::Cancelled)
    ));
    assert!(matches!(
        owner.step_candidate(
            balanced(0.02),
            InertialCapillaryIntervalControls {
                maximum_stopping_events: 0,
                ..controls()
            },
            || false
        ),
        Err(MaterialInertiaError::StoppingBudget { .. })
    ));
    assert!(matches!(
        owner.advance_candidate(
            balanced(0.02),
            controls(),
            CapillaryInertiaAdvanceControls {
                maximum_attempts: 1,
                ..adaptive()
            },
            || false
        ),
        Err(CapillaryInertiaAdvanceError::Attempts { .. })
    ));
    let mut polls = 0;
    assert!(matches!(
        owner.advance_candidate(balanced(0.02), controls(), adaptive(), || {
            polls += 1;
            polls > 20
        }),
        Err(CapillaryInertiaAdvanceError::Cancelled)
    ));
    assert_eq!(owner.state(), state(0.001));
    assert_eq!(displaced.reads.get(), 1);
}

#[test]
fn zero_time_and_a_short_final_tail_preserve_the_complete_requested_interval() {
    let (liquid, _) = fixture(1000.0, 0.05, 0.05, 1.0);
    let wetting = sourced_wetting(&liquid, 1.0, 2.0, 0.8, 1.2);
    let owner = MaterialCapillaryTube::new(&liquid, &wetting)
        .unwrap()
        .with_contact_line(limits())
        .unwrap()
        .with_inertia(parameters(), state(0.0))
        .unwrap();
    let c = CapillaryInertiaAdvanceControls {
        initial_step_s: 0.003,
        minimum_step_s: 0.003,
        maximum_step_s: 0.003,
        ..adaptive()
    };
    for duration in [0.0, 0.001, 0.01] {
        let candidate = owner
            .advance_candidate(balanced(duration), controls(), c, || false)
            .unwrap();
        assert_eq!(candidate.step().elapsed_time_s, duration);
        assert_eq!(candidate.step().state, state(0.0));
        assert_eq!(candidate.step().transferred_volume_m3, 0.0);
        assert_eq!(candidate.step().total_dissipation_j, 0.0);
        assert_eq!(candidate.step().kinetic_energy_change_j, 0.0);
        assert_eq!(candidate.step().energy_balance_residual_j, 0.0);
    }
    let near_end = MaterialCapillaryTube::new(&liquid, &wetting)
        .unwrap()
        .with_inertia(
            parameters(),
            InertialMeniscusState {
                wet_length_m: 0.1 - 1e-9,
                velocity_m_s: 0.001,
            },
        )
        .unwrap();
    assert!(matches!(
        near_end.step_candidate(
            InertialCapillaryBoundary {
                applied_pressure_pa: 100.0,
                opposing_gravity_m_s2: 0.0,
                duration_s: 1.0,
            },
            controls(),
            || false
        ),
        Err(MaterialInertiaError::Mechanics(
            InertialCapillaryError::StepTooLarge(_)
        ))
    ));
}
