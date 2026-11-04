//! Independently moving physical receivers through the native frame/PCM path.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::moving::trajectory::{
    MovingReceiverBudget, MovingReceiverConfig, PreparedKinematicTrajectory,
    PreparedMovingReceiverSource, TimedSourceKinematics, TrajectoryBudget,
};
use fs_aeroac::curle3d::moving::{CompactSourceStrength, SourceKinematics, TimedCompactSource};
use fs_couple::pcm_wav::encode_pcm16_wav_interleaved;
use fs_couple::pcm_wav::observation::PressureFrameRenderer;
use fs_couple::pcm_wav::observation::moving::{
    CompactSourceHistory, MovingCompactRenderer, MovingPressureBudget, MovingPressureConfig,
    frames::{MovingPressureFrameBudget, MovingPressureFrameRenderer},
};
use fs_couple::pcm_wav::stream::{Pcm16WavStream, WavFrameProgress, render_pressure_frames_pcm16};
use fs_couple::render::RenderError;
use fs_exec::CancelGate;
use fs_material::fluid_acoustics::FluidAcousticState;
use std::{f64::consts::PI, io::Cursor};

type Renderer<'a> = MovingCompactRenderer<'a, &'a PreparedMovingReceiverSource<'a>>;
type Frames<'a> = MovingPressureFrameRenderer<'static, 'a, &'a PreparedMovingReceiverSource<'a>, 2>;

fn path(origin: f64, velocity: f64) -> PreparedKinematicTrajectory {
    let times = if origin == 0.0 { [-0.2, 0.4] } else { [0.0, 0.16] };
    let samples = times.map(|source_time_s| TimedSourceKinematics {
        source_time_s,
        kinematics: SourceKinematics {
            position_m: [origin + velocity * source_time_s, 0.0, 0.0],
            velocity_m_s: [velocity, 0.0, 0.0],
            acceleration_m_s2: [0.0; 3],
        },
    });
    PreparedKinematicTrajectory::new(
        FluidAcousticState::try_new(1.2, 300.0).unwrap(), &samples,
        if velocity == 0.0 { 0.0 } else { 0.4 },
        TrajectoryBudget { maximum_segments: 8, maximum_payload_bytes: 131_072,
            maximum_bound_nodes: 4096, maximum_subdivision_depth: 16,
            maximum_root_iterations: 64 }, &mut || false,
    ).unwrap()
}
fn pair<'a>(source: &'a PreparedKinematicTrajectory, receiver: &'a PreparedKinematicTrajectory)
    -> PreparedMovingReceiverSource<'a>
{
    PreparedMovingReceiverSource::new(source, receiver,
        MovingReceiverConfig { receiver_time_interval_s: [0.0, 0.16],
            minimum_distance_m: 0.1, maximum_time_residual_s: 1e-12 },
        MovingReceiverBudget { maximum_segment_pairs: 64, maximum_bound_nodes: 4096,
            maximum_subdivision_depth: 16, maximum_root_iterations: 64 }, &mut || false,
    ).unwrap()
}
fn strength(t: f64) -> CompactSourceStrength {
    CompactSourceStrength {
        volume_flow_m3_s: 0.001 + 0.002*t + 0.003*t*t - 0.001*t*t*t,
        volume_flow_rate_m3_s2: 0.002 + 0.006*t - 0.003*t*t,
        force_on_fluid_n: [2.0 + 3.0*t - 2.0*t*t + t*t*t, 0.0, 0.0],
        force_rate_on_fluid_n_s: [3.0 - 4.0*t + 3.0*t*t, 0.0, 0.0],
    }
}
fn history() -> [TimedCompactSource; 4] {
    [-0.2, 0.0, 0.2, 0.4].map(|source_time_s| TimedCompactSource {
        source_time_s, strength: strength(source_time_s),
    })
}
fn config() -> MovingPressureConfig {
    MovingPressureConfig { first_receiver_time_s: 0.04, sample_rate_hz: 1000,
        samples: 101, max_block: 101, receiver_bandwidth_hz: 3.0,
        maximum_time_residual_s: 1e-12 }
}
fn renderer<'a>(motion: &'a PreparedMovingReceiverSource<'a>, samples: &'a [TimedCompactSource],
    config: MovingPressureConfig) -> Result<Renderer<'a>, RenderError>
{
    MovingCompactRenderer::with_moving_receiver(&[CompactSourceHistory {
        motion, samples, source_bandwidth_hz: 1.0,
    }], config, MovingPressureBudget { maximum_sources: 4, maximum_source_samples: 128,
        maximum_output_samples: 4096, maximum_payload_bytes: 131_072,
        maximum_source_receiver_pairs: 16_384 }, &mut || false)
}
fn budget() -> MovingPressureFrameBudget {
    MovingPressureFrameBudget { maximum_output_values: 8192,
        maximum_source_receiver_pairs: 32_768, maximum_retarded_iterations: 2_097_152,
        maximum_payload_bytes: 262_144 }
}
fn frames<'a>(left: &'a PreparedMovingReceiverSource<'a>, right: &'a PreparedMovingReceiverSource<'a>,
    samples: &'a [TimedCompactSource]) -> Frames<'a>
{
    MovingPressureFrameRenderer::new(["left", "right"],
        [renderer(left, samples, config()).unwrap(), renderer(right, samples, config()).unwrap()],
        budget(), &mut || false).unwrap()
}
fn bits(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
        expected.iter().map(|x| x.to_bits()).collect::<Vec<_>>());
}

#[test]
fn distinct_receiver_velocities_keep_their_own_retarded_pressure_and_mono_bits() {
    let source = path(0.0, 0.0);
    let (left_path, right_path) = (path(12.0, -30.0), path(20.0, 60.0));
    let (left, right) = (pair(&source, &left_path), pair(&source, &right_path));
    let samples = history();
    let mut expected = [[0.0; 101]; 2];
    for (channel, motion) in expected.iter_mut().zip([&left, &right]) {
        renderer(motion, &samples, config()).unwrap().block(channel).unwrap();
    }
    for block in [1, 7, 31, 101] {
        let mut array = frames(&left, &right, &samples);
        let mut output = [0.0; 202];
        for part in output.chunks_mut(2*block) { array.frame_block(part).unwrap(); }
        for (i, frame) in output.chunks_exact(2).enumerate() {
            let t = array.receiver_time_at(i as u64).unwrap();
            for (channel, (origin, velocity)) in [(12.0, -30.0), (20.0, 60.0)].into_iter().enumerate() {
                assert_eq!(frame[channel].to_bits(), expected[channel][i].to_bits());
                let radius = origin + velocity*t;
                let tau = t - radius/300.0;
                let s = strength(tau);
                // Stationary emitter: exact retarded Qdot/r, F/r^2 and Fdot/(c*r).
                // Receiver velocity changes the sampled emission clock, NOT gain.
                let analytic = (1.2*s.volume_flow_rate_m3_s2/radius
                    + s.force_on_fluid_n[0]/(radius*radius)
                    + s.force_rate_on_fluid_n_s[0]/(300.0*radius))/(4.0*PI);
                assert!((frame[channel]-analytic).abs() < 2e-12, "{i} {channel}");
            }
        }
        assert_eq!(array.frames_rendered(), 101);
        assert_eq!(array.remaining_frames(), 0);
        assert_eq!(array.channel_ids(), &["left", "right"]);
        assert!(array.channels().iter().all(|r| r.samples_rendered() == 101));
        assert_ne!(output[0].to_bits(), output[1].to_bits());
    }
}

#[test]
fn late_channel_cancellation_and_final_checkpoint_roll_back_every_cursor() {
    let source = path(0.0, 0.0);
    let (lp, rp) = (path(12.0, -30.0), path(20.0, 60.0));
    let (left, right) = (pair(&source, &lp), pair(&source, &rp));
    let samples = history();
    let mut left_polls = 0;
    renderer(&left, &samples, config()).unwrap().block_controlled(&mut [0.0; 9],
        &mut || { left_polls += 1; false }).unwrap();
    let mut expected = [0.0; 18];
    let mut all_polls = 0;
    frames(&left, &right, &samples).frame_block_controlled(&mut expected,
        &mut || { all_polls += 1; false }).unwrap();
    for cancel_at in [1, left_polls+2, all_polls] {
        let mut array = frames(&left, &right, &samples);
        let mut output = [123.0; 18];
        let mut polls = 0;
        assert!(matches!(array.frame_block_controlled(&mut output,
            &mut || { polls += 1; polls == cancel_at }),
            Err(RenderError::Acoustic(AeroacError::Cancelled { .. }))));
        assert_eq!(polls, cancel_at);
        assert_eq!(output, [123.0; 18]);
        assert!(array.channels().iter().all(|r| r.samples_rendered() == 0));
        array.frame_block(&mut output).unwrap();
        bits(&output, &expected);
        assert!(array.channels().iter().all(|r| r.samples_rendered() == 9));
    }
}

#[test]
fn aggregate_budgets_and_unaligned_clocks_refuse_instead_of_silent_resampling() {
    let source = path(0.0, 0.0);
    let (lp, rp) = (path(12.0, -30.0), path(20.0, 60.0));
    let (left, right) = (pair(&source, &lp), pair(&source, &rp));
    let samples = history();
    let measured = frames(&left, &right, &samples);
    let exact = MovingPressureFrameBudget { maximum_output_values: 202,
        maximum_source_receiver_pairs: measured.maximum_source_receiver_pairs(),
        maximum_retarded_iterations: measured.maximum_retarded_iterations(),
        maximum_payload_bytes: measured.retained_payload_bytes() };
    let build = |ids: [&'static str; 2], b, right_config| MovingPressureFrameRenderer::new(ids,
        [renderer(&left, &samples, config()).unwrap(), renderer(&right, &samples, right_config).unwrap()],
        b, &mut || false);
    assert!(build(["left", "right"], exact, config()).is_ok());
    for b in [MovingPressureFrameBudget { maximum_output_values: 201, ..exact },
        MovingPressureFrameBudget { maximum_source_receiver_pairs: exact.maximum_source_receiver_pairs-1, ..exact },
        MovingPressureFrameBudget { maximum_retarded_iterations: exact.maximum_retarded_iterations-1, ..exact },
        MovingPressureFrameBudget { maximum_payload_bytes: exact.maximum_payload_bytes-1, ..exact }]
    { assert!(build(["left", "right"], b, config()).is_err()); }
    for c in [MovingPressureConfig { first_receiver_time_s: 0.041, ..config() },
        MovingPressureConfig { sample_rate_hz: 1001, ..config() },
        MovingPressureConfig { samples: 100, ..config() }]
    { assert!(build(["left", "right"], budget(), c).is_err()); }
    for ids in [["same", "same"], [" ", "right"], ["left\n", "right"]] {
        assert!(build(ids, budget(), config()).is_err());
    }
    let mut started = renderer(&left, &samples, config()).unwrap();
    started.block(&mut [0.0; 1]).unwrap();
    assert!(MovingPressureFrameRenderer::new(["left", "right"],
        [started, renderer(&right, &samples, config()).unwrap()], budget(), &mut || false).is_err());
    // The approaching mic requires the original 1.1 source-clock expansion.
    assert!(renderer(&left, &samples, MovingPressureConfig { receiver_bandwidth_hz: 1.0, ..config() }).is_err());
}

#[test]
fn aligned_prefix_continues_and_complete_frame_shape_is_enforced() {
    let source = path(0.0, 0.0);
    let (lp, rp) = (path(12.0, -30.0), path(20.0, 60.0));
    let (left, right) = (pair(&source, &lp), pair(&source, &rp));
    let samples = history();
    let mut expected = [0.0; 30];
    frames(&left, &right, &samples).frame_block(&mut expected).unwrap();
    let mut observers = [renderer(&left, &samples, config()).unwrap(),
        renderer(&right, &samples, MovingPressureConfig { max_block: 7, ..config() }).unwrap()];
    for r in &mut observers { r.block(&mut [0.0; 3]).unwrap(); }
    let mut array = MovingPressureFrameRenderer::new(["left", "right"], observers,
        budget(), &mut || false).unwrap();
    assert_eq!(array.frames_rendered(), 3);
    assert_eq!(array.max_frame_block_len(), 7);
    for length in [0, 3, 16] {
        let mut invalid = vec![123.0; length];
        assert!(array.frame_block(&mut invalid).is_err());
        assert_eq!(invalid, vec![123.0; length]);
        assert_eq!(array.frames_rendered(), 3);
    }
    assert!(array.validate_frame_rate(999).is_err());
    assert!(array.validate_frame_count(99).is_err());
    let mut output = [0.0; 14];
    array.frame_block(&mut output).unwrap();
    bits(&output, &expected[6..20]);
    let mut originals = array.into_channels();
    for (channel, r) in originals.iter_mut().enumerate() {
        assert_eq!(r.samples_rendered(), 10);
        let mut next = [0.0; 5];
        r.block(&mut next).unwrap();
        for i in 0..5 { assert_eq!(next[i].to_bits(), expected[20+2*i+channel].to_bits()); }
    }
}

#[test]
fn native_multichannel_pcm_preserves_physical_order_and_pressure_calibration() {
    let source = path(0.0, 0.0);
    let (lp, rp) = (path(12.0, -30.0), path(20.0, 60.0));
    let (left, right) = (pair(&source, &lp), pair(&source, &rp));
    let samples = history();
    let mut reference = [0.0; 202];
    frames(&left, &right, &samples).frame_block(&mut reference).unwrap();
    let (expected, clips) = encode_pcm16_wav_interleaved(&reference, 1000, 2, 0.02).unwrap();
    for block in [1, 7, 32] {
        let mut array = frames(&left, &right, &samples);
        let mut stream = Pcm16WavStream::new_interleaved(Cursor::new(Vec::new()), 1000, 2, 0.02, block).unwrap();
        assert_eq!(render_pressure_frames_pcm16(&mut array, &mut stream,
            &CancelGate::new_clock_free(), &mut vec![0.0; 2*block], 101).unwrap(),
            WavFrameProgress::Completed { frames: 101 });
        let (sink, summary) = stream.finish().unwrap();
        assert_eq!(sink.into_inner(), expected);
        assert_eq!((summary.channels, summary.frames, summary.samples, summary.clipped_samples),
            (2, 101, 202, clips as u64));
    }
}

#[test]
fn late_physical_pressure_failure_does_not_consume_the_healthy_microphone() {
    let source = path(0.0, 0.0);
    let (lp, rp) = (path(12.0, -30.0), path(0.11, 0.0));
    let (left, right) = (pair(&source, &lp), pair(&source, &rp));
    let samples = [-0.2, 0.0, 0.2, 0.4].map(|source_time_s| TimedCompactSource {
        source_time_s,
        strength: CompactSourceStrength { force_on_fluid_n: [1e308, 0.0, 0.0],
            ..Default::default() },
    });
    let mut array = frames(&left, &right, &samples);
    let mut output = [123.0; 8];
    assert!(array.frame_block(&mut output).is_err());
    assert_eq!(output, [123.0; 8]);
    assert!(array.channels().iter().all(|r| r.samples_rendered() == 0));
    let mut originals = array.into_channels();
    originals[0].block(&mut [0.0; 4]).unwrap();
    assert_eq!(originals[0].samples_rendered(), 4);
    assert_eq!(originals[1].samples_rendered(), 0);
}
