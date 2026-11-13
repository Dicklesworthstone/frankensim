//! G1 physical volume-source identities and original-equation reuse; G4 atomic
//! cancellation. This is a source export, not a translating-bubble model.

#![cfg(feature = "cavitation")]

use fs_aeroac::curle3d::moving::{
    MovingCompactLimits, MovingCompactObserver, SourceKinematics, UniformSourceMotion,
};
use fs_aeroac::monopole::MonopoleLimits;
use fs_couple::cavitation::bubble::{MaterialBubble, MaterialBubbleError};
use fs_couple::pcm_wav::observation::moving::{
    CompactSourceHistory, MovingCompactRenderer, MovingPressureBudget, MovingPressureConfig,
};
use fs_flux::cavitation::{
    BubbleGas, BubbleLimits, BubbleNumerics, BubblePressureSample, BubbleState, CavitationError,
    KellerMiksisBubble, KellerMiksisLimits,
};
use fs_material::cavitation::CavitationState;
use fs_material::fluid_acoustics::FluidAcoustics;

const RADIUS: f64 = 100e-6;
const AMBIENT: f64 = 101_325.0;
const END_TIME: f64 = 50e-6;

fn liquid() -> CavitationState {
    CavitationState::try_new(293.15, 998.2, 0.001, 0.072, 2339.0, 1482.0).unwrap()
}

fn limits() -> BubbleLimits {
    BubbleLimits {
        minimum_radius_m: 0.5 * RADIUS,
        maximum_radius_m: 2.0 * RADIUS,
        maximum_wall_mach_number: 0.01,
        maximum_vapor_density_ratio: 0.01,
    }
}

fn numerics() -> BubbleNumerics {
    BubbleNumerics {
        relative_tolerance: 1e-9,
        radius_absolute_tolerance_m: 1e-14,
        velocity_absolute_tolerance_m_s: 1e-9,
        initial_step_s: 1e-8,
        minimum_step_s: 1e-12,
        maximum_step_s: 1e-7,
        maximum_attempts: 20_000,
    }
}

fn drive(time_s: f64) -> BubblePressureSample {
    let omega = 2.0 * core::f64::consts::PI * 2000.0;
    BubblePressureSample {
        absolute_pressure_pa: AMBIENT + 200.0 * (omega * time_s).sin(),
        time_derivative_pa_s: 200.0 * omega * (omega * time_s).cos(),
    }
}

fn rp(source: &CavitationState) -> MaterialBubble<'_, CavitationState> {
    MaterialBubble::new(
        source,
        BubbleGas::equilibrium(*source, RADIUS, AMBIENT, 1.4).unwrap(),
        limits(),
        BubbleState::try_new(0.0, RADIUS, 0.02).unwrap(),
        AMBIENT,
    )
    .unwrap()
}

fn km(source: &CavitationState) -> MaterialBubble<'_, CavitationState, KellerMiksisBubble> {
    MaterialBubble::new_compressible(
        source,
        BubbleGas::equilibrium(*source, RADIUS, AMBIENT, 1.4).unwrap(),
        KellerMiksisLimits {
            bubble: limits(),
            maximum_pressure_ratio: 0.01,
        },
        BubbleState::try_new(0.0, RADIUS, 0.0).unwrap(),
        drive(0.0),
    )
    .unwrap()
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-12 * expected.abs().max(1e-20),
        "{actual} != {expected}"
    );
}

#[test]
fn read_only_source_is_the_native_volume_derivative_on_the_same_clock() {
    let liquid = liquid();
    let bubble = rp(&liquid);
    let saved = bubble.state();
    let native = bubble.observe(AMBIENT).unwrap();
    let sample = bubble.compact_source_sample(AMBIENT).unwrap();
    let radius = saved.radius_m();
    let velocity = saved.wall_velocity_m_s();
    near(
        sample.strength.volume_flow_m3_s,
        4.0 * core::f64::consts::PI * radius * radius * velocity,
    );
    near(
        sample.strength.volume_flow_rate_m3_s2,
        4.0 * core::f64::consts::PI
            * radius
            * (2.0 * velocity * velocity + radius * native.wall_acceleration_m_s2),
    );
    assert_eq!(sample.source_time_s, saved.time_s());
    assert_eq!(sample.strength.volume_flow_m3_s, native.volume_rate_m3_s);
    assert_eq!(
        sample.strength.volume_flow_rate_m3_s2,
        native.volume_acceleration_m3_s2
    );
    assert_eq!(sample.strength.force_on_fluid_n, [0.0; 3]);
    assert_eq!(sample.strength.force_rate_on_fluid_n_s, [0.0; 3]);
    assert_eq!(bubble.state(), saved);
    assert!(bubble.compact_source_sample(f64::NAN).is_err());
    assert_eq!(bubble.state(), saved);
}

#[test]
fn rayleigh_plesset_export_preserves_native_dynamics_and_stationary_pressure() {
    let liquid = liquid();
    let mut bubble = rp(&liquid);
    let pressure = |time| drive(time).absolute_pressure_pa;
    let mut native_state = bubble.state();
    let native = bubble
        .model()
        .advance_controlled(
            &mut native_state,
            END_TIME,
            &numerics(),
            &pressure,
            &mut || false,
        )
        .unwrap();
    let exported = bubble
        .advance_compact_source(END_TIME, &numerics(), &pressure, &mut || false)
        .unwrap();
    assert_eq!(exported.bubble, native);
    assert_eq!(bubble.state(), native_state);
    assert_eq!(exported.source.source_time_s, END_TIME);
    assert_eq!(
        exported.source,
        bubble.compact_source_sample(pressure(END_TIME)).unwrap()
    );
    assert!(core::ptr::eq(bubble.source(), &liquid));
    assert_eq!(bubble.material(), liquid);

    let distance = 0.1;
    let old = bubble
        .monopole_observer(
            distance,
            100_000.0,
            MonopoleLimits {
                maximum_ka: 0.1,
                minimum_distance_radii: 10.0,
                maximum_pressure_ratio: 0.01,
            },
        )
        .unwrap()
        .sample(
            END_TIME,
            native_state.radius_m(),
            native.observation.volume_acceleration_m3_s2,
        )
        .unwrap();
    let observer = MovingCompactObserver::new(
        liquid.acoustic_state().unwrap(),
        [distance, 0.0, 0.0],
        MovingCompactLimits {
            maximum_mach: 0.0,
            minimum_distance_m: 0.01,
        },
    )
    .unwrap();
    let pressure = observer
        .prepare_emission(END_TIME, SourceKinematics::default())
        .unwrap()
        .evaluate(exported.source.strength)
        .unwrap();
    near(pressure.pressure_pa, old.pressure_pa);
    near(pressure.receiver_time_s, old.arrival_time_s);
    assert_eq!(pressure.dipole_pressure_pa, 0.0);
}

#[test]
fn keller_miksis_export_retains_the_analytic_drive_derivative_and_native_history() {
    let liquid = liquid();
    let mut bubble = km(&liquid);
    let initial = bubble.compact_source_sample(drive(0.0)).unwrap();
    // The initial compressibility response is nonzero even when wall speed
    // and instantaneous pressure perturbation are both zero.
    let acceleration = -drive(0.0).time_derivative_pa_s
        / (liquid.liquid_density_kg_m3() * liquid.liquid_sound_speed_m_s()
            + 4.0 * liquid.dynamic_viscosity_pa_s() / RADIUS);
    assert_eq!(initial.strength.volume_flow_m3_s, 0.0);
    near(
        initial.strength.volume_flow_rate_m3_s2,
        4.0 * core::f64::consts::PI * RADIUS * RADIUS * acceleration,
    );
    assert!(initial.strength.volume_flow_rate_m3_s2 < 0.0);

    let mut native_state = bubble.state();
    let native = bubble
        .model()
        .advance_controlled(
            &mut native_state,
            END_TIME,
            &numerics(),
            &drive,
            &mut || false,
        )
        .unwrap();
    let exported = bubble
        .advance_compact_source(END_TIME, &numerics(), &drive, &mut || false)
        .unwrap();
    assert_eq!(bubble.state(), native_state);
    assert_eq!(exported.bubble, native);
    assert_eq!(
        exported.source.strength.volume_flow_m3_s,
        native.observation.volume_rate_m3_s
    );
    assert_eq!(
        exported.source.strength.volume_flow_rate_m3_s2,
        native.observation.volume_acceleration_m3_s2
    );
    assert_eq!(
        exported.source,
        bubble.compact_source_sample(drive(END_TIME)).unwrap()
    );
}

#[test]
fn final_export_cancellation_is_atomic_and_retry_replays_the_physical_endpoint() {
    let liquid = liquid();
    let mut baseline = km(&liquid);
    let mut polls = 0usize;
    let expected = baseline
        .advance_compact_source(END_TIME, &numerics(), &drive, &mut || {
            polls += 1;
            false
        })
        .unwrap();
    let mut bubble = km(&liquid);
    let saved = bubble.state();
    let mut calls = 0usize;
    let result = bubble.advance_compact_source(END_TIME, &numerics(), &drive, &mut || {
        calls += 1;
        calls == polls
    });
    assert_eq!(
        result,
        Err(MaterialBubbleError::Dynamics(CavitationError::Cancelled))
    );
    assert_eq!(bubble.state(), saved);
    let retried = bubble
        .advance_compact_source(END_TIME, &numerics(), &drive, &mut || false)
        .unwrap();
    assert_eq!(retried, expected);
    assert_eq!(bubble.state(), baseline.state());
}

#[test]
fn accepted_bubble_history_renders_physical_microphone_pressure_without_source_evolution() {
    let liquid = liquid();
    let mut bubble = km(&liquid);
    let initial = bubble.observe(drive(0.0)).unwrap();
    let mut history = vec![bubble.compact_source_sample(drive(0.0)).unwrap()];
    let mut native_accelerations = vec![initial.volume_acceleration_m3_s2];
    for index in 1..=32_u32 {
        let time = f64::from(index) / 1_000_000.0;
        let step = bubble
            .advance_compact_source(time, &numerics(), &drive, &mut || false)
            .unwrap();
        assert_eq!(step.source.source_time_s, time);
        history.push(step.source);
        native_accelerations.push(step.bubble.observation.volume_acceleration_m3_s2);
    }
    let saved_state = bubble.state();
    let saved_history = history.clone();
    let distance = 0.1;
    let medium = bubble.material().acoustic_state().unwrap();
    let delay = distance / medium.sound_speed_m_s();
    let observer = MovingCompactObserver::new(
        medium,
        [distance, 0.0, 0.0],
        MovingCompactLimits {
            maximum_mach: 0.0,
            minimum_distance_m: 0.01,
        },
    )
    .unwrap();
    let sources = [CompactSourceHistory {
        motion: UniformSourceMotion::new(0.0, [0.0; 3], [0.0; 3]).unwrap(),
        samples: &history,
        source_bandwidth_hz: 100_000.0,
    }];
    // Observe interior source endpoints: retained samples on both sides keep
    // roundoff at a retarded endpoint inside the actual accepted history.
    let mut renderer = MovingCompactRenderer::new(
        observer,
        &sources,
        MovingPressureConfig {
            first_receiver_time_s: delay + history[1].source_time_s,
            sample_rate_hz: 1_000_000,
            samples: 31,
            max_block: 31,
            receiver_bandwidth_hz: 100_000.0,
            maximum_time_residual_s: 1e-15,
        },
        MovingPressureBudget {
            maximum_sources: 1,
            maximum_source_samples: 33,
            maximum_output_samples: 31,
            maximum_payload_bytes: 16_384,
            maximum_source_receiver_pairs: 31,
        },
        &mut || false,
    )
    .unwrap();
    let expected: Vec<f64> = native_accelerations[1..32]
        .iter()
        .map(|acceleration| {
            medium.density_kg_m3() * acceleration / (4.0 * core::f64::consts::PI * distance)
        })
        .collect();
    let scale = expected.iter().map(|p| p.abs()).fold(0.0_f64, f64::max);
    assert!(scale > 0.0);
    let mut pressure = [0.0; 31];
    renderer.block(&mut pressure).unwrap();
    for (index, (&actual, &expected)) in pressure.iter().zip(&expected).enumerate() {
        assert!(
            (actual - expected).abs() <= 2e-10 * scale,
            "sample {index}: {actual} Pa != native {expected} Pa"
        );
        near(
            renderer.receiver_time_at(index as u64).unwrap() - delay,
            history[index + 1].source_time_s,
        );
    }
    assert_eq!(renderer.samples_rendered(), 31);
    assert_eq!(bubble.state(), saved_state);
    assert_eq!(history, saved_history);
}
