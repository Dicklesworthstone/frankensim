//! G3/G4 original accelerating stress histories through the existing audio path.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::moving::trajectory::{
    PreparedSourceTrajectory, TimedSourceKinematics, TrajectoryBudget,
};
use fs_aeroac::curle3d::moving::{
    MovingCompactLimits, MovingCompactObserver, SourceKinematics, UniformSourceMotion,
};
use fs_aeroac::curle3d::quadrupole::transient::TransientQuadrupoleSource;
use fs_couple::pcm_wav::decimate::Decimator;
use fs_couple::pcm_wav::encode_pcm16_wav;
use fs_couple::pcm_wav::observation::DecimatedRenderer;
use fs_couple::pcm_wav::observation::moving::{CompactSourceMotion, MovingPressureConfig};
use fs_couple::pcm_wav::observation::quadrupole::QuadrupolePressureBudget;
use fs_couple::pcm_wav::observation::quadrupole::moving::{
    MovingQuadrupoleRenderer, MovingQuadrupoleSourceWindow, QuadrupoleMotion,
};
use fs_couple::pcm_wav::stream::{Pcm16WavStream, render_pressure_pcm16};
use fs_exec::CancelGate;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::dd::Dd;
use std::{cell::Cell, io::Cursor};

fn observer() -> MovingCompactObserver {
    MovingCompactObserver::new(
        FluidAcousticState::try_new(1.2, 340.0).unwrap(),
        [10.0, 2.0, 1.0],
        MovingCompactLimits {
            maximum_mach: 0.8,
            minimum_distance_m: 0.05,
        },
    )
    .unwrap()
}

fn trajectory() -> PreparedSourceTrajectory {
    let knots = [-0.1, 0.0, 0.1].map(|t: f64| {
        let v = [34.0, 0.0, 0.0];
        let a = [120.0, -40.0, 30.0];
        let j = [900.0, 150.0, -300.0];
        TimedSourceKinematics {
            source_time_s: t,
            kinematics: SourceKinematics {
                position_m: core::array::from_fn(|i| {
                    v[i] * t + a[i] * t * t / 2.0 + j[i] * t * t * t / 6.0
                }),
                velocity_m_s: core::array::from_fn(|i| v[i] + a[i] * t + j[i] * t * t / 2.0),
                acceleration_m_s2: core::array::from_fn(|i| a[i] + j[i] * t),
            },
        }
    });
    PreparedSourceTrajectory::new(
        observer(),
        &knots,
        TrajectoryBudget {
            maximum_segments: 2,
            maximum_payload_bytes: 100_000,
            maximum_bound_nodes: 8192,
            maximum_subdivision_depth: 20,
            maximum_root_iterations: 96,
        },
        &mut || false,
    )
    .unwrap()
}

fn config() -> MovingPressureConfig {
    MovingPressureConfig {
        first_receiver_time_s: 0.05,
        sample_rate_hz: 8000,
        samples: 128,
        max_block: 64,
        receiver_bandwidth_hz: 500.0,
        maximum_time_residual_s: 1e-12,
    }
}
fn budget() -> QuadrupolePressureBudget {
    QuadrupolePressureBudget {
        maximum_sources: 4,
        maximum_output_samples: 1000,
        maximum_source_path_samples: 4000,
        maximum_scratch_bytes: 4096,
    }
}
fn strength(index: usize, time: f64) -> Result<TransientQuadrupoleSource, AeroacError> {
    let omega = 2.0 * std::f64::consts::PI * 400.0;
    let (s, c) = (omega * time + index as f64 * 0.3).sin_cos();
    let q = [0.007, -0.003, 0.002, 0.004, 0.001, -0.005];
    Ok(TransientQuadrupoleSource {
        stress_moment_n_m: q.map(|v| v * c),
        stress_rate_n_m_s: q.map(|v| -omega * v * s),
        stress_acceleration_n_m_s2: q.map(|v| -omega * omega * v * c),
        ..Default::default()
    })
}
type Reader = fn(usize, f64) -> Result<TransientQuadrupoleSource, AeroacError>;
fn build<'a, M: QuadrupoleMotion<Observer = MovingCompactObserver> + 'a>(
    sources: &'a [MovingQuadrupoleSourceWindow<M>],
    cfg: MovingPressureConfig,
) -> MovingQuadrupoleRenderer<'a, Reader, M> {
    MovingQuadrupoleRenderer::new(
        observer(),
        sources,
        strength as Reader,
        cfg,
        budget(),
        &mut || false,
    )
    .unwrap()
}

#[test]
fn mixed_original_motion_owners_share_source_clocks_and_chunk_invariant_output() {
    let curved = trajectory();
    let uniform = UniformSourceMotion::new(0.0, [0.2, 0.1, -0.1], [20.0, -5.0, 0.0]).unwrap();
    let sources = [
        MovingQuadrupoleSourceWindow {
            motion: CompactSourceMotion::Trajectory(&curved),
            source_time_interval_s: [-0.1, 0.1],
            source_bandwidth_hz: 400.0,
        },
        MovingQuadrupoleSourceWindow {
            motion: CompactSourceMotion::Uniform(uniform),
            source_time_interval_s: [-0.1, 0.1],
            source_bandwidth_hz: 400.0,
        },
    ];
    let cfg = config();
    let mut one = build(&sources, cfg);
    assert_eq!(one.maximum_retarded_iterations(), 96 * cfg.samples);
    let expected: Vec<_> = (0..cfg.samples)
        .map(|index| {
            let time = one.receiver_time_at(index).unwrap();
            let a = curved
                .prepare_quadrupole_receiver(time, 1e-12, &mut || false)
                .unwrap();
            let b = observer()
                .prepare_uniform_quadrupole(time, &uniform, 1e-12)
                .unwrap();
            (Dd::from_f64(
                a.pressure_pa(strength(0, a.source_time_s()).unwrap())
                    .unwrap(),
            ) + Dd::from_f64(
                b.pressure_pa(strength(1, b.source_time_s()).unwrap())
                    .unwrap(),
            ))
            .to_f64()
        })
        .collect();
    let mut a = vec![0.0; cfg.samples as usize];
    for part in a.chunks_mut(64) {
        one.block(part).unwrap();
    }
    let mut many = build(&sources, cfg);
    let mut b = vec![0.0; cfg.samples as usize];
    for part in b.chunks_mut(7) {
        many.block(part).unwrap();
    }
    assert_eq!(a, expected);
    assert_eq!(b, expected);
    assert_eq!(many.samples_rendered(), cfg.samples);
}

#[test]
fn trajectory_history_observer_and_doppler_admission_precede_streaming() {
    let curved = trajectory();
    let sources = [MovingQuadrupoleSourceWindow {
        motion: &curved,
        source_time_interval_s: [-0.1, 0.1],
        source_bandwidth_hz: 400.0,
    }];
    let other =
        MovingCompactObserver::new(observer().medium(), [11.0, 2.0, 1.0], observer().limits())
            .unwrap();
    assert!(
        MovingQuadrupoleRenderer::new(
            other,
            &sources,
            strength as Reader,
            config(),
            budget(),
            &mut || false
        )
        .is_err()
    );
    for cfg in [
        MovingPressureConfig {
            receiver_bandwidth_hz: 400.0,
            ..config()
        },
        MovingPressureConfig {
            first_receiver_time_s: 0.2,
            ..config()
        },
    ] {
        assert!(
            MovingQuadrupoleRenderer::new(
                observer(),
                &sources,
                strength as Reader,
                cfg,
                budget(),
                &mut || false
            )
            .is_err()
        );
    }
    let missing = [MovingQuadrupoleSourceWindow {
        source_time_interval_s: [0.045, 0.1],
        ..sources[0]
    }];
    assert!(
        MovingQuadrupoleRenderer::new(
            observer(),
            &missing,
            strength as Reader,
            config(),
            budget(),
            &mut || false
        )
        .is_err()
    );
    assert!(
        MovingQuadrupoleRenderer::new(
            observer(),
            &sources,
            strength as Reader,
            config(),
            QuadrupolePressureBudget {
                maximum_source_path_samples: 127,
                ..budget()
            },
            &mut || false
        )
        .is_err()
    );
}

#[test]
fn failed_or_cancelled_accelerating_blocks_leave_output_and_clock_unchanged() {
    let curved = trajectory();
    let sources = [MovingQuadrupoleSourceWindow {
        motion: &curved,
        source_time_interval_s: [-0.1, 0.1],
        source_bandwidth_hz: 400.0,
    }];
    let reads = Cell::new(0);
    let fail = Cell::new(true);
    let reader = |i, t| {
        reads.set(reads.get() + 1);
        if fail.get() && reads.get() == 4 {
            return Err(AeroacError::InvalidParameter {
                what: "missing original stress event",
            });
        }
        strength(i, t)
    };
    let mut renderer = MovingQuadrupoleRenderer::new(
        observer(),
        &sources,
        reader,
        config(),
        budget(),
        &mut || false,
    )
    .unwrap();
    let mut output = [77.0; 8];
    assert!(renderer.block(&mut output).is_err());
    assert_eq!(output, [77.0; 8]);
    assert_eq!(renderer.samples_rendered(), 0);
    fail.set(false);
    let mut checks = 0;
    assert!(
        renderer
            .block_controlled(&mut output, &mut || {
                checks += 1;
                checks == 5
            })
            .is_err()
    );
    assert_eq!(output, [77.0; 8]);
    assert_eq!(renderer.samples_rendered(), 0);
    renderer.block(&mut output).unwrap();
    let mut expected = [0.0; 8];
    build(&sources, config()).block(&mut expected).unwrap();
    assert_eq!(output, expected);
}

#[test]
fn accelerating_pressure_reaches_the_original_decimator_and_pcm_wav_writer() {
    let curved = trajectory();
    let sources = [MovingQuadrupoleSourceWindow {
        motion: &curved,
        source_time_interval_s: [-0.1, 0.1],
        source_bandwidth_hz: 400.0,
    }];
    let cfg = MovingPressureConfig {
        sample_rate_hz: 16000,
        ..config()
    };
    let mut raw = build(&sources, cfg);
    let mut pressure = [0.0; 128];
    for part in pressure.chunks_mut(64) {
        raw.block(part).unwrap();
    }
    let mut filter = Decimator::new(2, 1).unwrap();
    let expected: Vec<_> = pressure
        .chunks_exact(2)
        .map(|part| {
            let value = filter.preview(part).unwrap()[0];
            filter.commit();
            value
        })
        .collect();
    let (expected_wav, clips) = encode_pcm16_wav(&expected, 8000, 0.01).unwrap();
    let mut renderer = DecimatedRenderer::new(build(&sources, cfg), 16000, 8000, 13).unwrap();
    let mut stream = Pcm16WavStream::new(Cursor::new(Vec::new()), 8000, 0.01, 13).unwrap();
    render_pressure_pcm16(
        &mut renderer,
        &mut stream,
        &CancelGate::new(),
        &mut [0.0; 13],
        64,
    )
    .unwrap();
    let (sink, summary) = stream.finish().unwrap();
    assert_eq!(sink.into_inner(), expected_wav);
    assert_eq!(summary.samples, 64);
    assert_eq!(summary.clipped_samples, clips as u64);
}
