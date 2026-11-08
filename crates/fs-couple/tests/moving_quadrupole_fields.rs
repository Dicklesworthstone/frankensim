//! Complete moving stress fields retain their causal near field and audio clock.

use fs_aeroac::curle3d::moving::trajectory::{
    MovingReceiverBudget, MovingReceiverConfig, MovingReflectionObserver,
    PreparedKinematicTrajectory, PreparedMovingReceiverSource,
    PreparedReflectedMovingReceiverSource, ReflectedTrajectoryBudget, TimedSourceKinematics,
    TrajectoryBudget,
};
use fs_aeroac::curle3d::moving::{
    MovingCompactLimits, MovingCompactObserver, SourceKinematics, UniformSourceMotion,
};
use fs_aeroac::curle3d::quadrupole::recorded::QuadrupoleSampleClock;
use fs_aeroac::curle3d::quadrupole::transient::history::{
    PreparedQuadrupoleHistory, QuadrupoleHistoryBudget, QuadrupoleHistoryEnds,
};
use fs_aeroac::curle3d::quadrupole::transient::{
    PreparedTransientQuadrupoleObserver, RetardedQuadrupoleField, TransientQuadrupoleSource,
};
use fs_aeroac::curle3d::{IdealAcousticPlane, IdealPlaneBoundary};
use fs_bem::near_field::FirstOrder;
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_couple::pcm_wav::observation::moving::MovingPressureConfig;
use fs_couple::pcm_wav::observation::moving::force::ForceHistoryQuadrature;
use fs_couple::pcm_wav::observation::quadrupole::moving::field::{
    MovingQuadrupoleFieldBudget, MovingQuadrupoleFieldRenderer,
};
use fs_couple::pcm_wav::observation::quadrupole::moving::{
    MovingQuadrupoleRenderer, MovingQuadrupoleSourceWindow,
};
use fs_couple::pcm_wav::observation::quadrupole::{
    QuadrupoleFieldSample, QuadrupolePressureBudget,
};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::dd::Dd;

const SAMPLES: usize = 17;
const TENSOR: [f64; 6] = [0.7, -0.2, 0.4, 0.3, -0.1, 0.23];

fn medium() -> FluidAcousticState {
    FluidAcousticState::try_new(1.2, 300.0).unwrap()
}

fn observer(position: [f64; 3]) -> MovingCompactObserver {
    MovingCompactObserver::new(
        medium(),
        position,
        MovingCompactLimits {
            maximum_mach: 0.8,
            minimum_distance_m: 0.01,
        },
    )
    .unwrap()
}

fn microphone() -> FirstOrder {
    FirstOrder::new(0.25, [0.0, 1.0, 0.0]).unwrap()
}

fn config() -> MovingPressureConfig {
    MovingPressureConfig {
        first_receiver_time_s: 0.10,
        sample_rate_hz: 1000,
        samples: SAMPLES as u64,
        max_block: SAMPLES,
        receiver_bandwidth_hz: 20.0,
        maximum_time_residual_s: 1e-12,
    }
}

fn budget() -> MovingQuadrupoleFieldBudget {
    MovingQuadrupoleFieldBudget {
        pressure: QuadrupolePressureBudget {
            maximum_sources: 2,
            maximum_output_samples: 128,
            maximum_source_path_samples: 512,
            maximum_scratch_bytes: 1_000_000,
        },
        maximum_history_bytes: 1_000_000,
        quadrature: ForceHistoryQuadrature {
            absolute_tolerance_n_s_m3: 1e-13,
            maximum_subdivision_depth: 20,
            maximum_evaluations: 1_000_000,
        },
    }
}

fn window<M>(motion: M) -> MovingQuadrupoleSourceWindow<M> {
    MovingQuadrupoleSourceWindow {
        motion,
        source_time_interval_s: [0.0, 0.3],
        source_bandwidth_hz: 1.0,
    }
}

// A clamped cubic reconstructs this original polynomial exactly. Its integral
// and derivatives below are an independent stationary-field oracle.
fn stress(t: f64, tensor: [f64; 6]) -> TransientQuadrupoleSource {
    TransientQuadrupoleSource {
        stress_integral_n_m_s: tensor.map(|q| q * (t.powi(3) / 3.0 + 0.075 * t.powi(4))),
        stress_moment_n_m: tensor.map(|q| q * (t * t + 0.3 * t.powi(3))),
        stress_rate_n_m_s: tensor.map(|q| q * (2.0 * t + 0.9 * t * t)),
        stress_acceleration_n_m_s2: tensor.map(|q| q * (2.0 + 1.8 * t)),
    }
}

fn history_with_onset(
    position: [f64; 3],
    tensor: [f64; 6],
    onset: TransientQuadrupoleSource,
) -> PreparedQuadrupoleHistory {
    let last = stress(0.3, tensor);
    PreparedQuadrupoleHistory::new(
        &[position],
        QuadrupoleSampleClock::new(0.0, 0.3, 0.0).unwrap(),
        &[0.0, 0.3],
        &[onset.stress_moment_n_m, last.stress_moment_n_m],
        &[QuadrupoleHistoryEnds {
            initial_rate_n_m_s: onset.stress_rate_n_m_s,
            final_rate_n_m_s: last.stress_rate_n_m_s,
            initial_integral_n_m_s: onset.stress_integral_n_m_s,
        }],
        QuadrupoleHistoryBudget {
            maximum_sources: 1,
            maximum_samples_per_source: 2,
            maximum_workspace_bytes: 16_384,
        },
        || false,
    )
    .unwrap()
}

fn history(position: [f64; 3], tensor: [f64; 6]) -> PreparedQuadrupoleHistory {
    history_with_onset(position, tensor, TransientQuadrupoleSource::default())
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-14 + 3e-9 * expected.abs(),
        "actual={actual:e}, expected={expected:e}"
    );
}

fn check_field(actual: QuadrupoleFieldSample, expected: RetardedQuadrupoleField) {
    close(actual.pressure_pa, expected.pressure_pa);
    for axis in 0..3 {
        close(
            actual.particle_velocity_m_s[axis],
            expected.particle_velocity_m_s[axis],
        );
        close(
            actual.instantaneous_intensity_w_m2[axis],
            expected.instantaneous_intensity_w_m2[axis],
        );
        assert_eq!(
            actual.instantaneous_intensity_w_m2[axis],
            actual.pressure_pa * actual.particle_velocity_m_s[axis]
        );
    }
    close(
        actual.microphone_pressure_pa,
        0.25 * expected.pressure_pa - 0.75 * 360.0 * expected.particle_velocity_m_s[1],
    );
}

#[test]
fn stationary_stress_recovers_transverse_near_fields_and_directional_audio() {
    let position = [0.0; 3];
    let receiver = [3.0, 0.0, 0.0];
    let sources = [window(
        UniformSourceMotion::new(0.0, position, [0.0; 3]).unwrap(),
    )];
    let oracle =
        PreparedTransientQuadrupoleObserver::new(medium(), receiver, &[position], 1).unwrap();
    for tensor in [[0.0, 0.0, 0.0, 1.0, 0.0, 0.0], TENSOR] {
        let history = history(position, tensor);
        let mut renderer = MovingQuadrupoleFieldRenderer::new_with_zero_stress_prehistory(
            observer(receiver),
            &sources,
            &history,
            microphone(),
            config(),
            budget(),
            &mut || false,
        )
        .unwrap();
        let mut fields = [QuadrupoleFieldSample::default(); SAMPLES];
        renderer.field_block(&mut fields).unwrap();
        for field in fields {
            let expected = oracle
                .evaluate_at_time(field.receiver_time_s, |_, t| Ok(stress(t, tensor)))
                .unwrap();
            check_field(field, expected);
            if tensor[0] == 0.0 {
                assert_eq!(field.pressure_pa, 0.0);
                assert_eq!(field.instantaneous_intensity_w_m2, [0.0; 3]);
                assert!(field.particle_velocity_m_s[1].abs() > 1e-7);
                assert!(field.microphone_pressure_pa.abs() > 1e-5);
            }
        }
    }
}

#[test]
fn blocks_audio_and_cancelled_retries_preserve_one_clock_and_charge_actual_work() {
    let position = [0.2, 0.1, 0.3];
    let sources = [window(
        UniformSourceMotion::new(0.0, position, [12.0, -3.0, 1.0]).unwrap(),
    )];
    let history = history(position, TENSOR);
    let make = |max_block| {
        MovingQuadrupoleFieldRenderer::new_with_zero_stress_prehistory(
            observer([5.0, 3.0, 2.0]),
            &sources,
            &history,
            microphone(),
            MovingPressureConfig {
                max_block,
                ..config()
            },
            budget(),
            &mut || false,
        )
        .unwrap()
    };
    let mut whole = make(SAMPLES);
    let mut expected = [QuadrupoleFieldSample::default(); SAMPLES];
    let mut complete_polls = 0;
    whole
        .field_block_controlled(&mut expected, &mut || {
            complete_polls += 1;
            false
        })
        .unwrap();
    let mut split = make(5);
    let mut actual = [QuadrupoleFieldSample::default(); SAMPLES];
    for chunk in actual.chunks_mut(5) {
        split.field_block(chunk).unwrap();
    }
    assert_eq!(actual, expected);
    assert_eq!(split.max_block_len(), 5);
    assert_eq!(
        split.quadrature_evaluations(),
        whole.quadrature_evaluations()
    );
    let mut audio = make(5);
    let mut samples = [0.0; SAMPLES];
    for chunk in samples.chunks_mut(5) {
        PressureRenderer::block(&mut audio, chunk).unwrap();
    }
    for (index, field) in expected.iter().enumerate() {
        assert_eq!(samples[index], field.microphone_pressure_pa);
        assert_eq!(
            field.receiver_time_s,
            audio.receiver_time_at(index as u64).unwrap()
        );
    }
    assert_eq!(audio.samples_rendered(), SAMPLES as u64);
    let mut past_end = [123.0];
    assert!(audio.block(&mut past_end).is_err());
    assert_eq!(past_end, [123.0]);
    assert_eq!(audio.samples_rendered(), SAMPLES as u64);

    let sentinel = QuadrupoleFieldSample {
        pressure_pa: 123.0,
        ..Default::default()
    };
    let mut output = [sentinel; SAMPLES];
    let mut retry = make(SAMPLES);
    let before = retry.quadrature_evaluations();
    let mut polls = 0;
    assert!(
        retry
            .field_block_controlled(&mut output, &mut || {
                polls += 1;
                polls == complete_polls / 2
            })
            .is_err()
    );
    assert_eq!(output, [sentinel; SAMPLES]);
    assert_eq!(retry.samples_rendered(), 0);
    let failed_work = retry.quadrature_evaluations() - before;
    assert!(failed_work > 0);
    retry.field_block(&mut output).unwrap();
    assert_eq!(output, expected);
    assert_eq!(
        retry.quadrature_evaluations(),
        whole.quadrature_evaluations() + failed_work
    );
}

#[test]
fn onset_binding_complete_past_and_finite_budgets_are_admission_gates() {
    let sources = [window(
        UniformSourceMotion::new(0.0, [0.0; 3], [0.0; 3]).unwrap(),
    )];
    let rejects = |history: &PreparedQuadrupoleHistory, bounds| {
        MovingQuadrupoleFieldRenderer::new_with_zero_stress_prehistory(
            observer([3.0, 0.0, 0.0]),
            &sources,
            history,
            microphone(),
            config(),
            bounds,
            &mut || false,
        )
        .is_err()
    };
    for component in 0..3 {
        let mut onset = TransientQuadrupoleSource::default();
        match component {
            0 => onset.stress_moment_n_m[0] = 1.0,
            1 => onset.stress_rate_n_m_s[0] = 1.0,
            _ => onset.stress_integral_n_m_s[0] = 1.0,
        }
        assert!(rejects(
            &history_with_onset([0.0; 3], TENSOR, onset),
            budget()
        ));
    }
    assert!(rejects(&history([0.1, 0.0, 0.0], TENSOR), budget()));
    let accepted = history([0.0; 3], TENSOR);
    assert!(rejects(
        &accepted,
        MovingQuadrupoleFieldBudget {
            maximum_history_bytes: 0,
            ..budget()
        }
    ));
    assert!(rejects(
        &accepted,
        MovingQuadrupoleFieldBudget {
            pressure: QuadrupolePressureBudget {
                maximum_source_path_samples: SAMPLES as u64 - 1,
                ..budget().pressure
            },
            ..budget()
        }
    ));

    // The first stress interval spans the output window: four evaluations
    // cannot finish its first adaptive Simpson integral, even across retries.
    let mut limited = MovingQuadrupoleFieldRenderer::new_with_zero_stress_prehistory(
        observer([3.0, 0.0, 0.0]),
        &sources,
        &accepted,
        microphone(),
        config(),
        MovingQuadrupoleFieldBudget {
            quadrature: ForceHistoryQuadrature {
                maximum_evaluations: 4,
                ..budget().quadrature
            },
            ..budget()
        },
        &mut || false,
    )
    .unwrap();
    let mut output = [123.0];
    for _ in 0..2 {
        assert!(limited.block(&mut output).is_err());
        assert_eq!(output, [123.0]);
        assert_eq!(limited.samples_rendered(), 0);
        assert_eq!(limited.quadrature_evaluations(), 4);
    }

    // Present retarded events are clear, but the source crossed the receiver
    // at t=0.05. The complete near-field past must not be silently discarded.
    let crossing = [window(
        UniformSourceMotion::new(0.0, [0.0; 3], [20.0, 0.0, 0.0]).unwrap(),
    )];
    let late = MovingPressureConfig {
        first_receiver_time_s: 0.20,
        ..config()
    };
    assert!(
        MovingQuadrupoleRenderer::new(
            observer([1.0, 0.0, 0.0]),
            &crossing,
            |_, t| accepted.source_at_clock_dd(0, Dd::from_f64(t)),
            late,
            budget().pressure,
            &mut || false,
        )
        .is_ok()
    );
    assert!(
        MovingQuadrupoleFieldRenderer::new_with_zero_stress_prehistory(
            observer([1.0, 0.0, 0.0]),
            &crossing,
            &accepted,
            microphone(),
            late,
            budget(),
            &mut || false,
        )
        .is_err()
    );
}

fn stationary(position_m: [f64; 3]) -> SourceKinematics {
    SourceKinematics {
        position_m,
        velocity_m_s: [0.0; 3],
        acceleration_m_s2: [0.0; 3],
    }
}

fn receiver(t: f64) -> SourceKinematics {
    SourceKinematics {
        position_m: [7.0 - 4.0 * t, 3.0 + 2.0 * t, 1.0 + 0.5 * t * t],
        velocity_m_s: [-4.0, 2.0, t],
        acceleration_m_s2: [0.0, 0.0, 1.0],
    }
}

fn path(motion: impl Fn(f64) -> SourceKinematics) -> PreparedKinematicTrajectory {
    let knots = [0.0, 0.04, 0.11, 0.3].map(|source_time_s| TimedSourceKinematics {
        source_time_s,
        kinematics: motion(source_time_s),
    });
    PreparedKinematicTrajectory::new(
        medium(),
        &knots,
        0.8,
        TrajectoryBudget {
            maximum_segments: 8,
            maximum_payload_bytes: 131_072,
            maximum_bound_nodes: 8192,
            maximum_subdivision_depth: 20,
            maximum_root_iterations: 96,
        },
        &mut || false,
    )
    .unwrap()
}

fn pair_config() -> MovingReceiverConfig {
    MovingReceiverConfig {
        receiver_time_interval_s: [0.08, 0.18],
        minimum_distance_m: 0.01,
        maximum_time_residual_s: 1e-12,
    }
}

fn pair_budget() -> MovingReceiverBudget {
    MovingReceiverBudget {
        maximum_segment_pairs: 128,
        maximum_bound_nodes: 32_768,
        maximum_subdivision_depth: 20,
        maximum_root_iterations: 96,
    }
}

fn image_budget() -> ReflectedTrajectoryBudget {
    ReflectedTrajectoryBudget {
        maximum_payload_bytes: 1_000_000,
        maximum_bound_nodes: 8192,
        maximum_subdivision_depth: 20,
    }
}

#[test]
fn moving_receiver_samples_the_current_spatial_field_without_extra_amplitude_doppler() {
    let position = [2.0, 1.0, 0.5];
    let source = path(|_| stationary(position));
    let receiver = path(receiver);
    let motion = PreparedMovingReceiverSource::new(
        &source,
        &receiver,
        pair_config(),
        pair_budget(),
        &mut || false,
    )
    .unwrap();
    let sources = [window(&motion)];
    let history = history(position, TENSOR);
    let mut renderer = MovingQuadrupoleFieldRenderer::new_with_zero_stress_prehistory(
        &receiver,
        &sources,
        &history,
        microphone(),
        config(),
        budget(),
        &mut || false,
    )
    .unwrap();
    let mut original = MovingQuadrupoleRenderer::new(
        &receiver,
        &sources,
        |_, t| history.source_at_clock_dd(0, Dd::from_f64(t)),
        config(),
        budget().pressure,
        &mut || false,
    )
    .unwrap();
    let mut pressure = [0.0; SAMPLES];
    original.block(&mut pressure).unwrap();
    let mut fields = [QuadrupoleFieldSample::default(); SAMPLES];
    renderer.field_block(&mut fields).unwrap();
    for (index, field) in fields.into_iter().enumerate() {
        assert_eq!(field.pressure_pa, pressure[index]);
        let current = receiver
            .kinematics_at(field.receiver_time_s)
            .unwrap()
            .position_m;
        let oracle =
            PreparedTransientQuadrupoleObserver::new(medium(), current, &[position], 1).unwrap();
        check_field(
            field,
            oracle
                .evaluate_at_time(field.receiver_time_s, |_, t| Ok(stress(t, TENSOR)))
                .unwrap(),
        );
    }
}

#[test]
fn oblique_rigid_and_pressure_release_images_transform_velocity_and_stress_once() {
    let position = [2.0, 1.0, 0.5];
    let source = path(|_| stationary(position));
    let receiver = path(receiver);
    let history = history(position, TENSOR);
    for boundary in [
        IdealPlaneBoundary::Rigid,
        IdealPlaneBoundary::PressureRelease,
    ] {
        let plane = IdealAcousticPlane::new([0.0; 3], [0.6, 0.8, 0.0], boundary).unwrap();
        let walls =
            MovingReflectionObserver::new(&receiver, &[plane], 2, image_budget(), &mut || false)
                .unwrap();
        let motion = PreparedReflectedMovingReceiverSource::new(
            &source,
            &walls,
            pair_config(),
            image_budget(),
            pair_budget(),
            &mut || false,
        )
        .unwrap();
        let sources = [window(&motion)];
        let mut renderer = MovingQuadrupoleFieldRenderer::new_with_zero_stress_prehistory(
            &walls,
            &sources,
            &history,
            microphone(),
            config(),
            budget(),
            &mut || false,
        )
        .unwrap();
        assert_eq!(renderer.source_count(), 1);
        assert_eq!(renderer.source_path_count(), 2);
        let mut fields = [QuadrupoleFieldSample::default(); SAMPLES];
        renderer.field_block(&mut fields).unwrap();
        for field in fields {
            let current = receiver
                .kinematics_at(field.receiver_time_s)
                .unwrap()
                .position_m;
            let oracle =
                PreparedTransientQuadrupoleObserver::new(medium(), current, &[position], 1)
                    .unwrap()
                    .with_orthogonal_planes(&[plane], 2)
                    .unwrap();
            check_field(
                field,
                oracle
                    .evaluate_at_time(field.receiver_time_s, |_, t| Ok(stress(t, TENSOR)))
                    .unwrap(),
            );
        }
    }
}
