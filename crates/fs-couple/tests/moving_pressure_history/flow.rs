//! G1/G3/G4: fixed laboratory source histories retain the original flow path.

use super::{budget, close, config, history, strength};
use fs_aeroac::curle3d::moving::MovingCompactLimits;
use fs_aeroac::curle3d::moving::flow::{UniformFlowCompactObserver, UniformFlowReceiver};
use fs_couple::pcm_wav::observation::moving::{
    CompactMotion, CompactSourceHistory, MovingCompactRenderer,
};
use fs_material::fluid_acoustics::FluidAcousticState;

fn receiver(rho: f64, flow: [f64; 3], position: [f64; 3]) -> UniformFlowReceiver {
    UniformFlowReceiver::new(
        FluidAcousticState::try_new(rho, 300.0).unwrap(),
        flow,
        position,
        MovingCompactLimits {
            maximum_mach: 0.9,
            minimum_distance_m: 0.01,
        },
    )
    .unwrap()
}

fn wind() -> UniformFlowReceiver {
    receiver(1.2, [180.0, 0.0, 0.0], [-12.0, 0.0, 0.0])
}

#[test]
fn original_cubic_histories_sum_in_the_same_flow_without_doppler_inflation() {
    let common = wind();
    let samples = history(&[-0.1, 0.0, 0.1, 0.2], 0.0, 1.0);
    let second = history(&[-0.1, 0.0, 0.1, 0.2], 0.0, -0.3);
    let sources = [
        CompactSourceHistory {
            motion: common.bind_source([0.0; 3]).unwrap(),
            samples: &samples,
            source_bandwidth_hz: 1.0,
        },
        CompactSourceHistory {
            motion: common.bind_source([1.0, -1.0, 0.5]).unwrap(),
            samples: &second,
            source_bandwidth_hz: 1.0,
        },
    ];
    let mut config = config();
    config.receiver_bandwidth_hz = 1.0;
    let build = || {
        MovingCompactRenderer::with_uniform_flow(&sources, config, budget(), &mut || false).unwrap()
    };
    let mut renderer = build();
    // The rest-frame stationary-point factor is >1 upstream. It must never
    // inflate the actual bandwidth of these fixed laboratory endpoints.
    assert!(
        sources[0]
            .motion
            .prepare_receiver(0.05, 1e-12)
            .unwrap()
            .pressure_in_fluid_frame()
            .source_time_rate()
            > 2.0
    );
    assert_eq!(*renderer.observer(), common);
    assert_eq!(renderer.maximum_retarded_iterations(), 0);
    let microphone = renderer.observer_at(config.first_receiver_time_s).unwrap();
    assert_eq!(microphone.medium(), common.medium());
    close(
        microphone.observer_position_m()[0],
        -12.0 - 180.0 * config.first_receiver_time_s,
    );
    let mut output = [0.0; 81];
    renderer.block(&mut output).unwrap();
    for (i, actual) in output.iter().enumerate() {
        let time = renderer.receiver_time_at(i as u64).unwrap();
        let mut expected = 0.0;
        for (source, scale) in sources.iter().zip([1.0, -0.3]) {
            let sample = source.motion.prepare_receiver(time, 1e-12).unwrap();
            let pressure = sample.pressure_in_fluid_frame();
            assert_eq!(pressure.source_time_s(), sample.source_time_s());
            assert_eq!(pressure.receiver_time_s(), time);
            expected += sample
                .evaluate(strength(sample.source_time_s(), scale))
                .unwrap()
                .pressure_pa;
        }
        close(*actual, expected);
    }
    let mut partitioned = build();
    let mut chunks = [0.0; 81];
    for chunk in chunks.chunks_mut(7) {
        partitioned.block(chunk).unwrap();
    }
    assert_eq!(output.map(f64::to_bits), chunks.map(f64::to_bits));
}

#[test]
fn shared_receiver_and_finite_history_bandwidth_are_required() {
    let samples = history(&[-0.1, 0.0, 0.1, 0.2], 0.0, 1.0);
    let original = CompactSourceHistory {
        motion: wind().bind_source([0.0; 3]).unwrap(),
        samples: &samples,
        source_bandwidth_hz: 1.0,
    };
    for other in [
        receiver(1.3, [180.0, 0.0, 0.0], [-12.0, 0.0, 0.0]),
        receiver(1.2, [179.0, 0.0, 0.0], [-12.0, 0.0, 0.0]),
        receiver(1.2, [180.0, 0.0, 0.0], [-11.0, 0.0, 0.0]),
    ] {
        let incompatible = CompactSourceHistory {
            motion: other.bind_source([1.0; 3]).unwrap(),
            ..original
        };
        assert!(
            MovingCompactRenderer::with_uniform_flow(
                &[original, incompatible],
                config(),
                budget(),
                &mut || false
            )
            .is_err()
        );
    }
    let mut narrow = config();
    narrow.receiver_bandwidth_hz = 0.999;
    assert!(
        MovingCompactRenderer::with_uniform_flow(&[original], narrow, budget(), &mut || false)
            .is_err()
    );
    let absent_prehistory = history(&[0.0, 0.1, 0.2], 0.0, 1.0);
    let missing = CompactSourceHistory {
        samples: &absent_prehistory,
        ..original
    };
    assert!(
        MovingCompactRenderer::with_uniform_flow(&[missing], config(), budget(), &mut || false)
            .is_err()
    );
    let unresolved = CompactSourceHistory {
        source_bandwidth_hz: 5.0,
        ..original
    };
    assert!(
        MovingCompactRenderer::with_uniform_flow(&[unresolved], config(), budget(), &mut || false)
            .is_err()
    );
}

#[test]
fn exact_history_endpoints_use_cached_geometry_at_large_epochs() {
    let motion = wind().bind_source([0.0; 3]).unwrap();
    for time in [f64::NAN, f64::INFINITY, -f64::INFINITY, 1e30] {
        assert!(motion.prepare_emission(time).is_err());
    }
    for epoch in [0.0, 1e12] {
        let times = [epoch + 0.1, epoch + 0.2];
        let samples = history(&times, -epoch, 1.0);
        let source = CompactSourceHistory {
            motion,
            samples: &samples,
            source_bandwidth_hz: 1.0,
        };
        let endpoint = motion.prepare_emission(times[1]).unwrap();
        let mut config = config();
        config.first_receiver_time_s = endpoint.receiver_time_s();
        config.samples = 1;
        config.maximum_time_residual_s = 1e-4;
        let native = motion
            .prepare_endpoint(&wind(), times[1], config.first_receiver_time_s)
            .unwrap();
        assert_eq!(native.source_time_s(), times[1]);
        assert_eq!(native.receiver_time_s(), endpoint.receiver_time_s());
        assert_eq!(
            native.distance_m(),
            motion
                .prepare_emission(0.0)
                .unwrap()
                .pressure_in_fluid_frame()
                .distance_m()
        );
        let expected = endpoint.evaluate(samples[1].strength).unwrap().pressure_pa;
        let mut renderer =
            MovingCompactRenderer::with_uniform_flow(&[source], config, budget(), &mut || false)
                .unwrap();
        let mut output = [0.0];
        renderer.block(&mut output).unwrap();
        assert_eq!(output[0].to_bits(), expected.to_bits());
        // A wider tolerance must not turn an actually later arrival into
        // extrapolated source history or a tolerance-sized endpoint clamp.
        config.first_receiver_time_s = config.first_receiver_time_s.next_up();
        assert!(
            MovingCompactRenderer::with_uniform_flow(&[source], config, budget(), &mut || false)
                .is_err()
        );
    }
}

#[test]
fn cancelled_flow_callback_preserves_output_and_retries_bitwise() {
    let samples = history(&[-0.1, 0.0, 0.1, 0.2], 0.0, 1.0);
    let source = CompactSourceHistory::<UniformFlowCompactObserver> {
        motion: wind().bind_source([0.0; 3]).unwrap(),
        samples: &samples,
        source_bandwidth_hz: 1.0,
    };
    let build = || {
        MovingCompactRenderer::with_uniform_flow(&[source], config(), budget(), &mut || false)
            .unwrap()
    };
    let mut renderer = build();
    let mut output = [123.0; 81];
    let mut polls = 0;
    assert!(
        renderer
            .block_controlled(&mut output, &mut || {
                polls += 1;
                polls == 6
            })
            .is_err()
    );
    assert_eq!(output, [123.0; 81]);
    assert_eq!(renderer.samples_rendered(), 0);
    renderer.block(&mut output).unwrap();
    let mut expected = [0.0; 81];
    build().block(&mut expected).unwrap();
    assert_eq!(output.map(f64::to_bits), expected.map(f64::to_bits));
}
