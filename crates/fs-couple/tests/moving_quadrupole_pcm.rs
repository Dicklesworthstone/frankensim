//! Retarded stress queries, continuous admission and the original WAV pipeline.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::moving::{MovingCompactLimits, MovingCompactObserver, UniformSourceMotion};
use fs_aeroac::curle3d::quadrupole::transient::TransientQuadrupoleSource;
use fs_couple::pcm_wav::decimate::Decimator;
use fs_couple::pcm_wav::encode_pcm16_wav;
use fs_couple::pcm_wav::observation::DecimatedRenderer;
use fs_couple::pcm_wav::observation::moving::MovingPressureConfig;
use fs_couple::pcm_wav::observation::quadrupole::QuadrupolePressureBudget;
use fs_couple::pcm_wav::observation::quadrupole::moving::{
    MovingQuadrupoleRenderer, MovingQuadrupoleSourceWindow,
};
use fs_couple::pcm_wav::stream::{Pcm16WavStream, render_pressure_pcm16};
use fs_exec::CancelGate;
use fs_material::fluid_acoustics::FluidAcousticState;
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
fn source() -> MovingQuadrupoleSourceWindow {
    MovingQuadrupoleSourceWindow {
        motion: UniformSourceMotion::new(0.0, [0.0; 3], [34.0, 0.0, 0.0]).unwrap(),
        source_time_interval_s: [-1.0, 1.0],
        source_bandwidth_hz: 400.0,
    }
}
fn strength(_index: usize, t: f64) -> Result<TransientQuadrupoleSource, AeroacError> {
    let omega = 2.0 * std::f64::consts::PI * 400.0;
    let (s, c) = (omega * t).sin_cos();
    let q = [0.007, -0.003, 0.002, 0.004, 0.001, -0.005];
    Ok(TransientQuadrupoleSource {
        stress_moment_n_m: q.map(|v| v * c),
        stress_rate_n_m_s: q.map(|v| -omega * v * s),
        stress_acceleration_n_m_s2: q.map(|v| -omega * omega * v * c),
        ..Default::default()
    })
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
type Reader = fn(usize, f64) -> Result<TransientQuadrupoleSource, AeroacError>;
fn build(
    sources: &[MovingQuadrupoleSourceWindow],
    cfg: MovingPressureConfig,
) -> MovingQuadrupoleRenderer<'_, Reader> {
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
fn original_retarded_source_clock_and_chunk_invariant_pressure() {
    let sources = [source()];
    let cfg = config();
    let mut one = build(&sources, cfg);
    let expected: Vec<_> = (0..cfg.samples)
        .map(|i| {
            let t = one.receiver_time_at(i).unwrap();
            let p = observer()
                .prepare_uniform_quadrupole(t, &sources[0].motion, 1e-12)
                .unwrap();
            p.pressure_pa(strength(0, p.source_time_s()).unwrap())
                .unwrap()
        })
        .collect();
    let mut a = vec![0.0; cfg.samples as usize];
    for chunk in a.chunks_mut(64) {
        one.block(chunk).unwrap();
    }
    let mut many = build(&sources, cfg);
    let mut b = vec![0.0; cfg.samples as usize];
    for chunk in b.chunks_mut(7) {
        many.block(chunk).unwrap();
    }
    assert_eq!(a, expected);
    assert_eq!(b, expected);
    assert_eq!(many.samples_rendered(), 128);
    assert!(many.validate_sample_count(1).is_err());
    assert!(many.validate_sample_rate(16000).is_err());
}

#[test]
fn signed_physical_sources_cancel_before_audio_conversion() {
    let sources = [source(), source()];
    let reader = |i, t| {
        let mut q = strength(i, t)?;
        if i == 1 {
            q.stress_moment_n_m = q.stress_moment_n_m.map(|v| -v);
            q.stress_rate_n_m_s = q.stress_rate_n_m_s.map(|v| -v);
            q.stress_acceleration_n_m_s2 = q.stress_acceleration_n_m_s2.map(|v| -v);
        }
        Ok(q)
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
    let mut output = [9.0; 32];
    renderer.block(&mut output).unwrap();
    assert_eq!(output, [0.0; 32]);
}

#[test]
fn source_read_failure_and_cancellation_do_not_publish_or_advance() {
    let sources = [source()];
    let failing = Cell::new(true);
    let reads = Cell::new(0);
    let reader = |i, t| {
        reads.set(reads.get() + 1);
        if failing.get() && reads.get() == 4 {
            return Err(AeroacError::InvalidParameter {
                what: "fixture missing accepted stress",
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
    let mut output = [99.0; 8];
    assert!(renderer.block(&mut output).is_err());
    assert_eq!(output, [99.0; 8]);
    assert_eq!(renderer.samples_rendered(), 0);
    failing.set(false);
    let mut checks = 0;
    assert!(
        renderer
            .block_controlled(&mut output, &mut || {
                checks += 1;
                checks == 5
            })
            .is_err()
    );
    assert_eq!(output, [99.0; 8]);
    assert_eq!(renderer.samples_rendered(), 0);
    renderer.block(&mut output).unwrap();
    let mut expected = [0.0; 8];
    build(&sources, config()).block(&mut expected).unwrap();
    assert_eq!(output, expected);
}

#[test]
fn whole_window_admission_covers_doppler_support_work_and_close_passes() {
    let original = [source()];
    let refused = |sources: &[MovingQuadrupoleSourceWindow], cfg, bounds| {
        assert!(
            MovingQuadrupoleRenderer::new(
                observer(),
                sources,
                strength as Reader,
                cfg,
                bounds,
                &mut || false
            )
            .is_err()
        );
    };
    // This tone fits below Nyquist but exceeds an unshifted bandwidth claim.
    refused(
        &original,
        MovingPressureConfig {
            receiver_bandwidth_hz: 400.0,
            ..config()
        },
        budget(),
    );
    refused(
        &original,
        MovingPressureConfig {
            receiver_bandwidth_hz: 4000.0,
            ..config()
        },
        budget(),
    );
    refused(
        &[MovingQuadrupoleSourceWindow {
            source_time_interval_s: [0.05, 1.0],
            ..source()
        }],
        config(),
        budget(),
    );
    refused(
        &[MovingQuadrupoleSourceWindow {
            source_time_interval_s: [-1.0, 0.02],
            ..source()
        }],
        config(),
        budget(),
    );
    refused(
        &original,
        config(),
        QuadrupolePressureBudget {
            maximum_source_path_samples: 127,
            ..budget()
        },
    );
    refused(
        &original,
        config(),
        QuadrupolePressureBudget {
            maximum_scratch_bytes: 511,
            ..budget()
        },
    );
    refused(
        &original,
        config(),
        QuadrupolePressureBudget {
            maximum_sources: 0,
            ..budget()
        },
    );
    refused(
        &original,
        MovingPressureConfig {
            first_receiver_time_s: 1e16,
            ..config()
        },
        budget(),
    );

    // Both endpoint ranges are allowed; the true closest approach lies between.
    let pass = [MovingQuadrupoleSourceWindow {
        motion: UniformSourceMotion::new(0.0, [9.0, 2.01, 1.0], [10.0, 0.0, 0.0]).unwrap(),
        source_bandwidth_hz: 0.0,
        ..source()
    }];
    refused(
        &pass,
        MovingPressureConfig {
            sample_rate_hz: 1000,
            samples: 101,
            receiver_bandwidth_hz: 100.0,
            ..config()
        },
        budget(),
    );
}

#[test]
fn moving_stress_uses_original_antialias_decimator_and_wav_writer() {
    let sources = [source()];
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
            let p = filter.preview(part).unwrap()[0];
            filter.commit();
            p
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
    assert_eq!(renderer.source().samples_rendered(), 128);
    let (sink, summary) = stream.finish().unwrap();
    assert_eq!(sink.into_inner(), expected_wav);
    assert_eq!(summary.samples, 64);
    assert_eq!(summary.clipped_samples, clips as u64);
}
