//! G1/G3 physical moving fields, original clocks and coherent directional audio.

#[path = "moving_monopole_field/paths.rs"]
mod paths;

use fs_aeroac::curle3d::moving::{
    CompactSourceStrength, MovingCompactLimits, MovingCompactObserver, TimedCompactSource,
    UniformSourceMotion,
};
use fs_bem::near_field::FirstOrder;
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_couple::pcm_wav::observation::moving::field::{CompactFieldSample, MovingMonopoleRenderer};
use fs_couple::pcm_wav::observation::moving::stationary::{
    StationaryCompactRenderer, StationaryCompactSample,
};
use fs_couple::pcm_wav::observation::moving::{
    CompactSourceHistory, MovingCompactRenderer, MovingPressureBudget, MovingPressureConfig,
};
use fs_material::fluid_acoustics::FluidAcousticState;
use std::f64::consts::PI;

fn medium() -> FluidAcousticState {
    FluidAcousticState::try_new(1.2, 300.0).unwrap()
}
fn observer(position: [f64; 3]) -> MovingCompactObserver {
    MovingCompactObserver::new(
        medium(),
        position,
        MovingCompactLimits {
            maximum_mach: 0.9,
            minimum_distance_m: 0.01,
        },
    )
    .unwrap()
}
fn budget() -> MovingPressureBudget {
    MovingPressureBudget {
        maximum_sources: 8,
        maximum_source_samples: 2048,
        maximum_output_samples: 4096,
        maximum_payload_bytes: 131_072,
        maximum_source_receiver_pairs: 65_536,
    }
}
fn config() -> MovingPressureConfig {
    MovingPressureConfig {
        first_receiver_time_s: 0.05,
        sample_rate_hz: 1000,
        samples: 81,
        max_block: 81,
        receiver_bandwidth_hz: 10.0,
        maximum_time_residual_s: 1e-12,
    }
}
// One original polynomial is exactly represented on arbitrary Hermite knots.
fn strength(t: f64, scale: f64) -> CompactSourceStrength {
    CompactSourceStrength {
        volume_flow_m3_s: scale * (0.001 + 0.002 * t - 0.003 * t * t + 0.004 * t.powi(3)),
        volume_flow_rate_m3_s2: scale * (0.002 - 0.006 * t + 0.012 * t * t),
        ..Default::default()
    }
}
fn history(scale: f64) -> [TimedCompactSource; 5] {
    [-0.25, -0.1, 0.05, 0.2, 0.5].map(|source_time_s| TimedCompactSource {
        source_time_s,
        strength: strength(source_time_s, scale),
    })
}
fn source(samples: &[TimedCompactSource], velocity: f64) -> CompactSourceHistory<'_> {
    CompactSourceHistory {
        motion: UniformSourceMotion::new(0.0, [0.0; 3], [velocity, 0.0, 0.0]).unwrap(),
        samples,
        source_bandwidth_hz: 1.0,
    }
}
fn renderer<'a>(
    sources: &[CompactSourceHistory<'a>],
    microphone: FirstOrder,
) -> MovingMonopoleRenderer<'a> {
    MovingMonopoleRenderer::new(
        observer([12.0, 0.0, 0.0]),
        sources,
        microphone,
        config(),
        budget(),
        &mut || false,
    )
    .unwrap()
}
fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-12),
        "actual={actual:e}, expected={expected:e}"
    );
}

// Independent collinear Green-potential reduction, with source derivatives
// on tau. Receiver speed changes tau, never the Eulerian radiation gain.
fn collinear(time: f64, source_speed: f64, receiver_speed: f64) -> (f64, f64) {
    let x = 12.0 + receiver_speed * time;
    let kappa = 1.0 - source_speed / 300.0;
    let tau = (time - x / 300.0) / kappa;
    let radius = x - source_speed * tau;
    let s = strength(tau, 1.0);
    let pressure = 1.2
        * (s.volume_flow_rate_m3_s2 / radius
            + s.volume_flow_m3_s * source_speed / (radius * radius))
        / (4.0 * PI * kappa * kappa);
    let velocity = (s.volume_flow_m3_s / (radius * radius)
        + s.volume_flow_rate_m3_s2 / (300.0 * radius))
        / (4.0 * PI * kappa * kappa);
    (pressure, velocity)
}

#[test]
fn moving_fields_retain_induction_source_clock_and_ideal_proximity_response() {
    let samples = history(1.0);
    let microphone = FirstOrder::new(0.5, [-1.0, 0.0, 0.0]).unwrap();
    for speed in [-60.0, 0.0, 60.0] {
        let sources = [source(&samples, speed)];
        let mut fields = renderer(&sources, microphone);
        let mut output = [CompactFieldSample::default(); 81];
        fields.field_block(&mut output).unwrap();
        let mut pressure_owner = MovingCompactRenderer::new(
            observer([12.0, 0.0, 0.0]),
            &sources,
            config(),
            budget(),
            &mut || false,
        )
        .unwrap();
        let mut pressure = [0.0; 81];
        pressure_owner.block(&mut pressure).unwrap();
        for (index, field) in output.iter().enumerate() {
            let time = fields.receiver_time_at(index as u64).unwrap();
            let (p, u) = collinear(time, speed, 0.0);
            assert_eq!(field.receiver_time_s, time);
            close(field.pressure_pa, p, 2e-12);
            close(field.particle_velocity_m_s[0], u, 2e-12);
            assert_eq!(&field.particle_velocity_m_s[1..], &[0.0; 2]);
            close(field.instantaneous_intensity_w_m2[0], p * u, 4e-12);
            close(field.microphone_pressure_pa, 0.5 * (p + 360.0 * u), 2e-12);
            assert_eq!(field.pressure_pa, pressure[index]);
            assert!((u - p / 360.0).abs() > 0.5 * u.abs());
        }
        assert_eq!(fields.maximum_retarded_iterations(), 0);
    }
}

#[test]
fn coherent_fields_interfere_before_intensity_or_directional_pickup() {
    let positive = history(1.0);
    let negative = history(-1.0);
    let a = source(&positive, 60.0);
    let b = source(&negative, 60.0);
    let mic = FirstOrder::new(0.3, [-1.0, 0.0, 0.0]).unwrap();
    let mut one = renderer(&[a], mic);
    let mut two = renderer(&[a, a], mic);
    let mut cancelling = renderer(&[a, b], mic);
    let mut single = [CompactFieldSample::default(); 81];
    let mut doubled = single;
    let mut null = single;
    one.field_block(&mut single).unwrap();
    two.field_block(&mut doubled).unwrap();
    cancelling.field_block(&mut null).unwrap();
    for ((single, double), zero) in single.into_iter().zip(doubled).zip(null) {
        close(double.pressure_pa, 2.0 * single.pressure_pa, 2e-12);
        close(
            double.microphone_pressure_pa,
            2.0 * single.microphone_pressure_pa,
            2e-12,
        );
        for axis in 0..3 {
            close(
                double.particle_velocity_m_s[axis],
                2.0 * single.particle_velocity_m_s[axis],
                2e-12,
            );
            close(
                double.instantaneous_intensity_w_m2[axis],
                4.0 * single.instantaneous_intensity_w_m2[axis],
                2e-12,
            );
        }
        assert_eq!(zero.pressure_pa, 0.0);
        assert_eq!(zero.particle_velocity_m_s, [0.0; 3]);
        assert_eq!(zero.instantaneous_intensity_w_m2, [0.0; 3]);
        assert_eq!(zero.microphone_pressure_pa, 0.0);
    }
}

#[test]
fn stationary_public_fields_and_moving_zero_speed_limit_agree() {
    let samples = history(1.0);
    let sources = [source(&samples, 0.0)];
    let mic = FirstOrder::new(0.5, [-1.0, 0.0, 0.0]).unwrap();
    let mut stationary = StationaryCompactRenderer::new(
        observer([12.0, 0.0, 0.0]),
        &sources,
        &[[0.0; 3]],
        mic,
        config(),
        budget(),
        &mut || false,
    )
    .unwrap();
    let mut moving = renderer(&sources, mic);
    // The preserved stationary alias is the same physical field record.
    let mut expected = [StationaryCompactSample::default(); 81];
    let mut observed = [CompactFieldSample::default(); 81];
    stationary.field_block(&mut expected).unwrap();
    moving.field_block(&mut observed).unwrap();
    for (actual, expected) in observed.iter().zip(expected) {
        close(actual.pressure_pa, expected.pressure_pa, 2e-12);
        close(
            actual.microphone_pressure_pa,
            expected.microphone_pressure_pa,
            2e-12,
        );
        for axis in 0..3 {
            close(
                actual.particle_velocity_m_s[axis],
                expected.particle_velocity_m_s[axis],
                2e-12,
            );
            close(
                actual.instantaneous_intensity_w_m2[axis],
                expected.instantaneous_intensity_w_m2[axis],
                2e-12,
            );
        }
    }
}

#[test]
fn field_and_audio_blocks_share_one_partition_invariant_transactional_clock() {
    let samples = history(1.0);
    let sources = [source(&samples, 60.0)];
    let mic = FirstOrder::new(0.5, [-1.0, 0.0, 0.0]).unwrap();
    let mut whole = renderer(&sources, mic);
    let mut expected = [CompactFieldSample::default(); 81];
    whole.field_block(&mut expected).unwrap();
    let mut split = renderer(&sources, mic);
    let mut fields = [CompactFieldSample::default(); 7];
    split.field_block(&mut fields).unwrap();
    assert_eq!(fields, expected[..7]);
    let mut microphone = [0.0; 74];
    for block in microphone.chunks_mut(11) {
        PressureRenderer::block(&mut split, block).unwrap();
    }
    assert_eq!(
        microphone,
        expected[7..]
            .iter()
            .map(|f| f.microphone_pressure_pa)
            .collect::<Vec<_>>()
            .as_slice()
    );
    assert_eq!(split.samples_rendered(), 81);
    assert!(split.block(&mut [0.0]).is_err());
    assert!(split.validate_sample_rate(999).is_err());
    assert_eq!(split.samples_rendered(), 81);

    let mut retry = renderer(&sources, mic);
    let sentinel = CompactFieldSample {
        pressure_pa: 123.0,
        ..Default::default()
    };
    let mut output = [sentinel; 81];
    let mut polls = 0;
    assert!(
        retry
            .field_block_controlled(&mut output, &mut || {
                polls += 1;
                polls == 20
            })
            .is_err()
    );
    assert_eq!(output, [sentinel; 81]);
    assert_eq!(retry.samples_rendered(), 0);
    retry.field_block(&mut output).unwrap();
    assert_eq!(output, expected);
}

#[test]
fn force_histories_and_unbudgeted_field_storage_refuse_before_output() {
    let mut samples = history(1.0);
    for derivative_only in [false, true] {
        if derivative_only {
            samples[3].strength.force_rate_on_fluid_n_s = [0.0, 1e-30, 0.0];
            samples[3].strength.force_on_fluid_n = [0.0; 3];
        } else {
            samples[3].strength.force_on_fluid_n = [1e-30, 0.0, 0.0];
        }
        let sources = [source(&samples, 60.0)];
        assert!(
            MovingMonopoleRenderer::new(
                observer([12.0, 0.0, 0.0]),
                &sources,
                FirstOrder::default(),
                config(),
                budget(),
                &mut || false,
            )
            .is_err()
        );
    }
    let samples = history(1.0);
    let sources = [source(&samples, 60.0)];
    assert!(
        MovingMonopoleRenderer::new(
            observer([12.0, 0.0, 0.0]),
            &sources,
            FirstOrder::default(),
            config(),
            MovingPressureBudget {
                maximum_payload_bytes: config().max_block * size_of::<CompactFieldSample>(),
                ..budget()
            },
            &mut || false,
        )
        .is_err()
    );
}
