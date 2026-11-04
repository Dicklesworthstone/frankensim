//! Accepted moving-source histories in flow use the original receiver renderer.
use fs_aeroac::curle3d::moving::flow::UniformFlowFrame;
use fs_aeroac::curle3d::moving::trajectory::{
    MovingReceiverBudget, MovingReceiverConfig, PreparedKinematicTrajectory,
    PreparedMovingReceiverSource, TimedSourceKinematics, TrajectoryBudget,
};
use fs_aeroac::curle3d::moving::{CompactSourceStrength, SourceKinematics, TimedCompactSource};
use fs_couple::pcm_wav::observation::moving::{
    CompactSourceHistory, MovingCompactRenderer, MovingPressureBudget, MovingPressureConfig,
};
use fs_material::fluid_acoustics::FluidAcousticState;

fn path(times: &[f64], kind: usize) -> PreparedKinematicTrajectory {
    let samples: Vec<_> = times
        .iter()
        .map(|&t| TimedSourceKinematics {
            source_time_s: t,
            kinematics: match kind {
                0 => SourceKinematics {
                    position_m: [20.0 * t + 2.0 * t * t, 0.5, 0.0],
                    velocity_m_s: [20.0 + 4.0 * t, 0.0, 0.0],
                    acceleration_m_s2: [4.0, 0.0, 0.0],
                },
                1 => SourceKinematics {
                    position_m: [-1.0 + 10.0 * t, -0.5 + 2.0 * t * t, 1.0],
                    velocity_m_s: [10.0, 4.0 * t, 0.0],
                    acceleration_m_s2: [0.0, 4.0, 0.0],
                },
                _ => SourceKinematics {
                    position_m: [8.0 + 4.0 * t, 3.0 - 2.0 * t, -1.0],
                    velocity_m_s: [4.0, -2.0, 0.0],
                    acceleration_m_s2: [0.0; 3],
                },
            },
        })
        .collect();
    PreparedKinematicTrajectory::in_uniform_flow(
        FluidAcousticState::try_new(1.2, 300.0).unwrap(),
        &samples,
        UniformFlowFrame::new([70.0, -20.0, 15.0], 0.0, [1.0, -2.0, 0.5]).unwrap(),
        0.8,
        TrajectoryBudget {
            maximum_segments: 8,
            maximum_payload_bytes: 131_072,
            maximum_bound_nodes: 16_384,
            maximum_subdivision_depth: 20,
            maximum_root_iterations: 96,
        },
        &mut || false,
    )
    .unwrap()
}

fn pair<'a>(
    source: &'a PreparedKinematicTrajectory,
    receiver: &'a PreparedKinematicTrajectory,
) -> PreparedMovingReceiverSource<'a> {
    PreparedMovingReceiverSource::new(
        source,
        receiver,
        MovingReceiverConfig {
            receiver_time_interval_s: [0.0, 0.16],
            minimum_distance_m: 0.1,
            maximum_time_residual_s: 1e-12,
        },
        MovingReceiverBudget {
            maximum_segment_pairs: 64,
            maximum_bound_nodes: 65_536,
            maximum_subdivision_depth: 24,
            maximum_root_iterations: 96,
        },
        &mut || false,
    )
    .unwrap()
}

fn strength(t: f64, scale: f64) -> CompactSourceStrength {
    CompactSourceStrength {
        volume_flow_m3_s: scale * (0.001 + 0.002 * t + 0.003 * t * t - 0.001 * t * t * t),
        volume_flow_rate_m3_s2: scale * (0.002 + 0.006 * t - 0.003 * t * t),
        force_on_fluid_n: [
            scale * (2.0 + 3.0 * t - 2.0 * t * t + t * t * t),
            scale * (-1.0 + 2.0 * t * t),
            scale * (0.5 - 0.4 * t),
        ],
        force_rate_on_fluid_n_s: [
            scale * (3.0 - 4.0 * t + 3.0 * t * t),
            scale * 4.0 * t,
            scale * (-0.4),
        ],
    }
}

#[test]
fn moving_flow_sources_render_coherently_on_the_original_clock_and_callback_path() {
    let receiver = path(&[0.0, 0.1, 0.2], 2);
    let a = path(&[-0.3, 0.0, 0.4], 0);
    let b = path(&[-0.3, 0.0, 0.4], 1);
    let pairs = [pair(&a, &receiver), pair(&b, &receiver)];
    let histories: Vec<Vec<_>> = [1.0, -0.4]
        .iter()
        .map(|&scale| {
            [-0.3, 0.0, 0.4]
                .iter()
                .map(|&time| TimedCompactSource {
                    source_time_s: time,
                    strength: strength(time, scale),
                })
                .collect()
        })
        .collect();
    let sources: Vec<_> = pairs
        .iter()
        .zip(&histories)
        .map(|(motion, samples)| CompactSourceHistory {
            motion,
            samples,
            source_bandwidth_hz: 1.0,
        })
        .collect();
    let config = MovingPressureConfig {
        first_receiver_time_s: 0.04,
        sample_rate_hz: 1000,
        samples: 101,
        max_block: 101,
        receiver_bandwidth_hz: 10.0,
        maximum_time_residual_s: 1e-12,
    };
    let budget = MovingPressureBudget {
        maximum_sources: 4,
        maximum_source_samples: 128,
        maximum_output_samples: 4096,
        maximum_payload_bytes: 262_144,
        maximum_source_receiver_pairs: 1_000_000,
    };
    let build = || {
        MovingCompactRenderer::with_moving_receiver(&sources, config, budget, &mut || false)
            .unwrap()
    };
    let mut renderer = build();
    assert!(core::ptr::eq(*renderer.observer(), &receiver));
    let mut output = [0.0; 101];
    renderer.block(&mut output).unwrap();
    for (i, &pressure) in output.iter().enumerate() {
        let time = renderer.receiver_time_at(i as u64).unwrap();
        let mut expected = 0.0;
        for (pair, scale) in pairs.iter().zip([1.0, -0.4]) {
            let sample = pair.prepare_receiver(time, &mut || false).unwrap();
            expected += sample
                .pressure
                .evaluate(strength(sample.pressure.source_time_s(), scale))
                .unwrap()
                .pressure_pa;
        }
        assert!((pressure - expected).abs() < 2e-12 * expected.abs().max(0.01));
    }
    let mut partitioned = build();
    let mut small_blocks = [0.0; 101];
    for block in small_blocks.chunks_mut(9) {
        partitioned.block(block).unwrap();
    }
    assert_eq!(output.map(f64::to_bits), small_blocks.map(f64::to_bits));

    let mut cancelled = build();
    let mut pending = [123.0; 101];
    let mut polls = 0;
    assert!(
        cancelled
            .block_controlled(&mut pending, &mut || {
                polls += 1;
                polls == 12
            })
            .is_err()
    );
    assert_eq!(pending, [123.0; 101]);
    assert_eq!(cancelled.samples_rendered(), 0);
    cancelled.block(&mut pending).unwrap();
    assert_eq!(output.map(f64::to_bits), pending.map(f64::to_bits));
}
