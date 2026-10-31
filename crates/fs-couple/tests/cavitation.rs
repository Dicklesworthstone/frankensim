//! G1/G4/G5 integration of original water, conserved film inventory, bubble
//! dynamics and causal acoustic observation. No measured-accuracy claim.

use fs_aeroac::monopole::{CompactMonopoleObserver, MonopoleLimits};
use fs_couple::cavitation::bubble::{MaterialBubble, MaterialBubbleError};
use fs_couple::cavitation::{
    CavitationFilmEnvironment, CavitationSolveConfig, MaterialCavitatingFilm,
    MaterialCavitationFilmError,
};
use fs_flux::cavitation::{BubbleGas, BubbleLimits, BubbleNumerics, BubbleState, CavitationError};
use fs_material::cavitation::{
    CavitationProperties, CavitationState, SaturatedVaporState, SaturatedWaterCavitation,
};
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_qty::Temperature;
use fs_thermochem::water::if97::If97SaturationPair;
use fs_tribo::resistive_film::{FilmCell, FilmChannel, FilmLimits, GapPort};

const AMBIENT_PA: f64 = 101_325.0;
const AREA_M2: f64 = 1e-4;
const INITIAL_GAP_M: f64 = 20e-6;
const CHANNEL_GAP_M: f64 = 10e-6;
const RADIUS_M: f64 = 100e-6;

fn water(temperature_k: f64) -> SaturatedWaterCavitation {
    SaturatedWaterCavitation::try_new(
        If97SaturationPair::at_temperature(Temperature::new(temperature_k)).unwrap(),
        0.072, // Prescribed test interface; no temperature-dependent sigma fit is claimed.
    )
    .unwrap()
}
fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-10 * expected.abs().max(1e-15),
        "{actual} != {expected}"
    );
}
fn film<F: CavitationProperties>(source: &F) -> MaterialCavitatingFilm<'_, F> {
    try_film(source, 0.01).unwrap()
}
fn try_film<F: CavitationProperties>(
    source: &F,
    maximum_vapor_density_ratio: f64,
) -> Result<MaterialCavitatingFilm<'_, F>, MaterialCavitationFilmError> {
    let gap = |reference_m, closure| GapPort {
        reference_m,
        closure: vec![closure],
    };
    MaterialCavitatingFilm::new(
        source,
        vec![FilmCell {
            area_m2: AREA_M2,
            gap: gap(INITIAL_GAP_M, 1.0),
        }],
        vec![FilmChannel {
            from: 0,
            to: None,
            width_m: 1e-3,
            length_m: 1e-2,
            gap: gap(CHANNEL_GAP_M, 0.0),
        }],
        1,
        FilmLimits {
            minimum_cell_gap_m: 1e-6,
            maximum_gap_m: 1e-3,
            maximum_pressure_pa: 2e5,
        },
        CavitationFilmEnvironment {
            ambient_pressure_pa: AMBIENT_PA,
            ambient_liquid_fraction: 1.0,
            maximum_vapor_density_ratio,
            solve: CavitationSolveConfig {
                maximum_iterations: 32,
                absolute_mass_tolerance_kg: 1e-17,
                relative_mass_tolerance: 1e-11,
            },
        },
    )
}

#[test]
fn film_rejects_dense_vapor_and_requires_an_explicit_positive_density_ratio_limit() {
    let liquid = CavitationState::try_new(300.0, 1000.0, 0.001, 0.072, 2000.0, 1500.0).unwrap();
    let dense = liquid
        .with_saturated_vapor(SaturatedVaporState::try_new(50.0, 2.2e6).unwrap())
        .unwrap();
    let dilute = liquid
        .with_saturated_vapor(SaturatedVaporState::try_new(1.0, 2.2e6).unwrap())
        .unwrap();
    assert!(matches!(
        try_film(&dense, 0.01),
        Err(MaterialCavitationFilmError::Film(_))
    ));
    assert!(try_film(&dilute, 0.01).is_ok());
    // Missing optional vapor data is a declared assumption, not invented data.
    assert!(try_film(&liquid, 0.01).is_ok());
    for invalid in [0.0, -0.01, 0.1001, f64::NAN, f64::INFINITY] {
        assert!(try_film(&dilute, invalid).is_err());
    }
}
fn bubble(source: &SaturatedWaterCavitation) -> MaterialBubble<'_, SaturatedWaterCavitation> {
    MaterialBubble::new(
        source,
        BubbleGas::equilibrium(source.state(), RADIUS_M, AMBIENT_PA, 1.4).unwrap(),
        BubbleLimits {
            minimum_radius_m: 0.5 * RADIUS_M,
            maximum_radius_m: 2.0 * RADIUS_M,
            maximum_wall_mach_number: 0.01,
            maximum_vapor_density_ratio: 0.01,
        },
        BubbleState::try_new(0.0, RADIUS_M, 0.0).unwrap(),
        AMBIENT_PA,
    )
    .unwrap()
}
fn numerics() -> BubbleNumerics {
    BubbleNumerics {
        relative_tolerance: 1e-8,
        radius_absolute_tolerance_m: 1e-13,
        velocity_absolute_tolerance_m_s: 1e-8,
        initial_step_s: 1e-8,
        minimum_step_s: 1e-12,
        maximum_step_s: 1e-7,
        maximum_attempts: 10_000,
    }
}
fn acoustic_limits() -> MonopoleLimits {
    MonopoleLimits {
        maximum_ka: 0.1,
        minimum_distance_radii: 10.0,
        maximum_pressure_ratio: 0.01,
    }
}
fn pressure(t: f64) -> f64 {
    AMBIENT_PA + 200.0 * (2.0 * std::f64::consts::PI * 2000.0 * t).sin()
}

#[test]
fn original_water_controls_film_viscosity_pressure_force_and_conserved_mass() {
    let mut inflows = Vec::new();
    for temperature in [293.15, 333.15] {
        let source = water(temperature);
        let film = film(&source);
        let material = source.state();
        assert!(std::ptr::eq(film.source(), &source));
        assert_eq!(film.material(), material);
        assert_eq!(
            film.film().geometry().viscosity_pa_s(),
            material.dynamic_viscosity_pa_s()
        );
        assert_eq!(
            film.film().liquid().density_kg_m3,
            material.liquid_density_kg_m3()
        );
        assert_eq!(
            film.film().liquid().cavitation_pressure_pa,
            material.saturation_pressure_pa()
        );
        let initial = film.initial_state(&[0.0]).unwrap();
        let saved = initial.clone();
        let dt = 0.01;
        let step = film
            .advance(&initial, &[-INITIAL_GAP_M], &[0.0], dt, || false)
            .unwrap();
        // Independent single passage Poiseuille law and swept-volume inventory.
        let conductance =
            1e-3 * CHANNEL_GAP_M.powi(3) / (12.0 * material.dynamic_viscosity_pa_s() * 1e-2);
        let inflow = conductance * (AMBIENT_PA - material.saturation_pressure_pa());
        let volume_before = AREA_M2 * INITIAL_GAP_M;
        near(
            step.state.liquid_fraction()[0],
            (volume_before + dt * inflow) / (2.0 * volume_before),
        );
        assert_eq!(
            step.state.pressure_pa(),
            &[material.saturation_pressure_pa()]
        );
        near(step.channel_volume_flow_m3_s[0], -inflow);
        near(
            step.generalized_force[0],
            AREA_M2 * (material.saturation_pressure_pa() - AMBIENT_PA),
        );
        near(
            step.report.liquid_mass_before_kg,
            material.liquid_density_kg_m3() * volume_before,
        );
        near(
            step.report.boundary_mass_out_kg,
            -material.liquid_density_kg_m3() * dt * inflow,
        );
        assert!(step.report.mass_residual_kg.abs() < 1e-18);
        assert_eq!(initial, saved);
        let refill = film
            .advance(&step.state, &[0.0], &[0.0], dt, || false)
            .unwrap();
        assert_eq!(refill.state.liquid_fraction(), &[1.0]);
        near(
            refill.state.pressure_pa()[0],
            2.0 * AMBIENT_PA - material.saturation_pressure_pa(),
        );
        near(
            refill.report.liquid_mass_after_kg,
            step.report.liquid_mass_before_kg,
        );
        inflows.push(inflow);
    }
    assert!(inflows[1] > inflows[0]); // Actual water transport changes the flow.
}

#[test]
fn film_history_replays_identical_content_flux_and_receipts() {
    let source = water(293.15);
    let film = film(&source);
    let initial = film.initial_state(&[0.0]).unwrap();
    let mut a = initial.clone();
    let mut b = initial;
    for q in [
        -INITIAL_GAP_M,
        -0.5 * INITIAL_GAP_M,
        0.0,
        -INITIAL_GAP_M,
        0.0,
    ] {
        let first = film.advance(&a, &[q], &[0.0], 0.01, || false).unwrap();
        let second = film.advance(&b, &[q], &[0.0], 0.01, || false).unwrap();
        assert_eq!(first, second);
        a = first.state;
        b = second.state;
    }
}

#[test]
fn one_material_and_clock_drive_repeatable_bubble_and_causal_pressure() {
    let source = water(293.15);
    let mut first = bubble(&source);
    let mut second = bubble(&source);
    let observer = first
        .monopole_observer(0.1, 100_000.0, acoustic_limits())
        .unwrap();
    assert!(std::ptr::eq(first.source(), &source));
    assert_eq!(first.model().liquid(), source.state());
    assert_eq!(observer.medium(), source.acoustic_state().unwrap());
    for endpoint in [10e-6, 20e-6, 30e-6] {
        let a = first
            .advance_observed(endpoint, &numerics(), &pressure, &observer, &mut || false)
            .unwrap();
        let b = second
            .advance_observed(endpoint, &numerics(), &pressure, &observer, &mut || false)
            .unwrap();
        assert_eq!(a, b);
        assert_eq!(first.state(), second.state());
        assert_eq!(a.acoustic.source_time_s, first.state().time_s());
        near(
            a.acoustic.arrival_time_s - a.acoustic.source_time_s,
            0.1 / source.state().liquid_sound_speed_m_s(),
        );
        let volume_acceleration = a.bubble.observation.volume_acceleration_m3_s2;
        near(
            a.acoustic.pressure_pa,
            source.state().liquid_density_kg_m3() * volume_acceleration
                / (4.0 * std::f64::consts::PI * 0.1),
        );
    }
}

#[test]
fn foreign_medium_and_late_observer_refusal_cannot_publish_bubble_candidate() {
    let source = water(293.15);
    let mut bubble = bubble(&source);
    let accepted = bubble.state();
    let foreign = CompactMonopoleObserver::new(
        FluidAcousticState::try_new(500.0, 1000.0).unwrap(),
        0.1,
        100_000.0,
        acoustic_limits(),
    )
    .unwrap();
    assert_eq!(
        bubble.advance_observed(10e-6, &numerics(), &pressure, &foreign, &mut || false),
        Err(MaterialBubbleError::ObserverMediumMismatch)
    );
    assert_eq!(bubble.state(), accepted);
    // Observer construction has no source radius. Its domain fails only after
    // the complete native bubble interval has produced a valid candidate.
    let too_close = bubble
        .monopole_observer(5.0 * RADIUS_M, 100_000.0, acoustic_limits())
        .unwrap();
    assert!(matches!(
        bubble.advance_observed(10e-6, &numerics(), &pressure, &too_close, &mut || false),
        Err(MaterialBubbleError::Acoustic(_))
    ));
    assert_eq!(bubble.state(), accepted);
}

#[test]
fn final_publication_cancellation_preserves_state_and_retry_replays() {
    let source = water(293.15);
    let mut reference = bubble(&source);
    let observer = reference
        .monopole_observer(0.1, 100_000.0, acoustic_limits())
        .unwrap();
    let mut total_polls = 0;
    let expected = reference
        .advance_observed(10e-6, &numerics(), &pressure, &observer, &mut || {
            total_polls += 1;
            false
        })
        .unwrap();
    assert!(total_polls > 10);
    let mut interrupted = bubble(&source);
    let accepted = interrupted.state();
    let mut polls = 0;
    let result =
        interrupted.advance_observed(10e-6, &numerics(), &pressure, &observer, &mut || {
            polls += 1;
            polls == total_polls
        });
    assert_eq!(
        result,
        Err(MaterialBubbleError::Dynamics(CavitationError::Cancelled))
    );
    assert_eq!(polls, total_polls);
    assert_eq!(interrupted.state(), accepted);
    let retried = interrupted
        .advance_observed(10e-6, &numerics(), &pressure, &observer, &mut || false)
        .unwrap();
    assert_eq!(retried, expected);
    assert_eq!(interrupted.state(), reference.state());
}
