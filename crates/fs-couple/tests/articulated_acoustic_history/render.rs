//! Complete native articulated frames through retarded pressure and PCM.

use super::*;
use fs_aeroac::curle3d::moving::{CompactSourceStrength, SourceKinematics};
use fs_aeroac::curle3d::moving::trajectory::{
    MovingReceiverBudget, MovingReceiverConfig, PreparedKinematicTrajectory,
    PreparedMovingReceiverSource, TrajectoryBudget,
};
use fs_couple::pcm_wav::encode_pcm16_wav;
use fs_couple::pcm_wav::observation::moving::{
    CompactSourceHistory, MovingCompactRenderer, MovingPressureBudget, MovingPressureConfig,
};
use fs_couple::pcm_wav::observation::moving::rigid::articulated::{
    LinkAcousticUpdate, append_link_frame,
};
use fs_couple::pcm_wav::stream::{Pcm16WavStream, render_pressure_pcm16};
use fs_exec::CancelGate;
use fs_material::fluid_acoustics::FluidAcousticState;
use std::{f64::consts::PI, io::Cursor};

fn update<'a>(history: &'a mut RigidAcousticHistory, motion: &'a PreparedLinkPointMotion)
    -> LinkAcousticUpdate<'a> {
    LinkAcousticUpdate { history, motion, strength: BodyFrameCompactSourceStrength::default() }
}

#[test]
fn complete_link_frame_rejects_a_late_source_or_cancellation_without_partial_histories() {
    let motion = quiet_motion();
    let mut a = history(3);
    let mut b = history(2);
    let mut scratch = [None; 2];
    let strength = BodyFrameCompactSourceStrength::default();
    append_link_frame(0.0, &mut [update(&mut a, &motion), update(&mut b, &motion)], &mut scratch, &mut || false).unwrap();
    let accepted_a = (a.kinematics().to_vec(), a.sources().to_vec());
    let accepted_b = (b.kinematics().to_vec(), b.sources().to_vec());
    for cancel_at in 1..=4 {
        let mut polls = 0;
        assert!(append_link_frame(0.1, &mut [update(&mut a, &motion), update(&mut b, &motion)], &mut scratch,
            &mut || { polls += 1; polls == cancel_at }).is_err());
        assert_eq!(polls, cancel_at);
        assert_eq!(a.kinematics(), accepted_a.0.as_slice());
        assert_eq!(a.sources(), accepted_a.1.as_slice());
        assert_eq!(b.kinematics(), accepted_b.0.as_slice());
        assert_eq!(b.sources(), accepted_b.1.as_slice());
    }
    let bad = LinkAcousticUpdate { history: &mut b, motion: &motion,
        strength: BodyFrameCompactSourceStrength { volume_flow_m3_s: f64::NAN, ..strength } };
    assert!(append_link_frame(0.1, &mut [update(&mut a, &motion), bad], &mut scratch, &mut || false).is_err());
    assert_eq!((a.len(), b.len()), (1, 1));
    // Retry overwrites scratch left by every earlier failure.
    append_link_frame(0.1, &mut [update(&mut a, &motion), update(&mut b, &motion)], &mut scratch, &mut || false).unwrap();
    assert!(scratch.iter().all(Option::is_none));
    assert_eq!((a.len(), b.len()), (2, 2));
    // Last history is now full; the first still has a free slot.
    assert!(append_link_frame(0.2, &mut [update(&mut a, &motion), update(&mut b, &motion)], &mut scratch, &mut || false).is_err());
    assert_eq!((a.len(), b.len()), (2, 2));
    assert_eq!(a.sources()[1].source_time_s, b.sources()[1].source_time_s);
}

#[test]
fn frame_shape_is_checked_and_old_scratch_never_becomes_a_new_event() {
    let motion = quiet_motion();
    let mut a = history(2);
    let mut b = history(2);
    let stale = a.attachment().sample_link(-20.0, &motion,
        BodyFrameCompactSourceStrength::default()).unwrap();
    let mut scratch = [Some(stale); 2];
    assert!(append_link_frame(1.0, &mut [update(&mut a, &motion), update(&mut b, &motion)], &mut scratch[..1], &mut || false).is_err());
    assert!(a.is_empty() && b.is_empty());
    assert!(append_link_frame(1.0, &mut [], &mut [], &mut || false).is_err());
    append_link_frame(1.0, &mut [update(&mut a, &motion), update(&mut b, &motion)], &mut scratch, &mut || false).unwrap();
    for h in [&a, &b] {
        assert_eq!(h.sources()[0].source_time_s, 1.0);
        assert_eq!(h.kinematics()[0].source_time_s, 1.0);
    }
}

const ARMS: [[f64; 3]; 3] = [[0.25, 0.15, -0.1], [-0.4, 0.3, 0.2], [8.0, 3.0, 1.0]];
fn source(t: f64, scale: f64) -> BodyFrameCompactSourceStrength {
    BodyFrameCompactSourceStrength {
        volume_flow_m3_s: scale * (0.0001 + 0.0002*t + 0.0001*t*t),
        volume_flow_rate_m3_s2: scale * (0.0002 + 0.0002*t),
        force_on_fluid_body_n: [scale*(1.0+0.2*t), scale*(0.4-0.3*t), scale*(-0.2+0.1*t*t)],
        force_rate_on_fluid_body_n_s: [scale*0.2, -scale*0.3, scale*0.2*t],
    }
}
fn record(intervals: usize) -> [RigidAcousticHistory; 3] {
    let model = slider();
    let gravity = Vec3::new(0.0, 0.0, -9.81);
    let mut histories = ARMS.map(|a| RigidAcousticHistory::new(
        attachment(Vec3::new(a[0], a[1], a[2])),
        RigidHistoryBudget { maximum_samples: intervals+1, maximum_payload_bytes: 1_048_576 }).unwrap());
    let mut scratch = [None; 3];
    for i in 0..=intervals {
        let t = -0.1 + 0.4 * i as f64 / intervals as f64;
        let q = [0.2+2.0*t+0.4*t*t, 0.8+0.3*t+0.2*t*t];
        let rates = [2.0+0.8*t, 0.3+0.4*t];
        let base = BaseState::prescribed(
            Se3::from_parts(So3::identity(), Vec3::new(0.7*t, -0.2*t, 0.1*t*t)).unwrap(),
            Twist::new(Vec3::default(), Vec3::new(0.7, -0.2, 0.2*t)),
            Twist::new(Vec3::default(), Vec3::new(0.0, 0.0, 0.2)));
        let kin = forward_kinematics(&model, base, &q, &rates).unwrap();
        let dyns = inverse_dynamics(&model, base, &q, &rates, &[0.8, 0.4], gravity,
            &[Wrench::default(); 2]).unwrap();
        let root = dyns.prepare_point_motion(&kin, 0, gravity).unwrap();
        let sliding = dyns.prepare_point_motion(&kin, 1, gravity).unwrap();
        let [a, b, microphone] = &mut histories;
        append_link_frame(t, &mut [
            LinkAcousticUpdate { history: a, motion: &sliding, strength: source(t, 1.0) },
            LinkAcousticUpdate { history: b, motion: &root, strength: source(t, -0.4) },
            LinkAcousticUpdate { history: microphone, motion: &root, strength: Default::default() },
        ], &mut scratch, &mut || false).unwrap();
    }
    histories
}
fn path(h: &RigidAcousticHistory) -> PreparedKinematicTrajectory {
    PreparedKinematicTrajectory::new(FluidAcousticState::try_new(1.2, 300.0).unwrap(),
        h.kinematics(), 0.2, TrajectoryBudget {
            maximum_segments: 128, maximum_payload_bytes: 1_048_576,
            maximum_bound_nodes: 65_536, maximum_subdivision_depth: 20,
            maximum_root_iterations: 96,
        }, &mut || false).unwrap()
}
fn pair<'a>(s: &'a PreparedKinematicTrajectory, r: &'a PreparedKinematicTrajectory)
    -> PreparedMovingReceiverSource<'a> {
    PreparedMovingReceiverSource::new(s, r, MovingReceiverConfig {
        receiver_time_interval_s: [0.04, 0.141], minimum_distance_m: 0.2,
        maximum_time_residual_s: 1e-12,
    }, MovingReceiverBudget {
        maximum_segment_pairs: 16_384, maximum_bound_nodes: 262_144,
        maximum_subdivision_depth: 24, maximum_root_iterations: 96,
    }, &mut || false).unwrap()
}

fn observe(intervals: usize, exercise_output: bool) -> [f64; 101] {
    let histories = record(intervals);
    let paths = histories.each_ref().map(path);
    let pairs = [pair(&paths[0], &paths[2]), pair(&paths[1], &paths[2])];
    let tracks = [0, 1].map(|i| CompactSourceHistory {
        motion: &pairs[i], samples: histories[i].sources(), source_bandwidth_hz: 4.0,
    });
    let config = MovingPressureConfig {
        first_receiver_time_s: 0.04, sample_rate_hz: 1000, samples: 101, max_block: 101,
        receiver_bandwidth_hz: 12.0, maximum_time_residual_s: 1e-12,
    };
    let budget = MovingPressureBudget {
        maximum_sources: 2, maximum_source_samples: 300, maximum_output_samples: 101,
        maximum_payload_bytes: 1_048_576, maximum_source_receiver_pairs: 202,
    };
    let build = || MovingCompactRenderer::with_moving_receiver(&tracks, config, budget, &mut || false).unwrap();
    let mut output = [0.0; 101];
    build().block(&mut output).unwrap();
    if exercise_output {
        let mut split = build();
        let mut partitioned = [0.0; 101];
        for block in partitioned.chunks_mut(13) { split.block(block).unwrap(); }
        assert_eq!(output.map(f64::to_bits), partitioned.map(f64::to_bits));
        let mut cancelled = build();
        let mut pending = [123.0; 101];
        let mut polls = 0;
        assert!(cancelled.block_controlled(&mut pending,
            &mut || { polls += 1; polls == 12 }).is_err());
        assert_eq!(pending, [123.0; 101]);
        assert_eq!(cancelled.samples_rendered(), 0);
        cancelled.block(&mut pending).unwrap();
        assert_eq!(output.map(f64::to_bits), pending.map(f64::to_bits));
        let (expected_wav, clips) = encode_pcm16_wav(&output, 1000, 0.02).unwrap();
        let mut stream = Pcm16WavStream::new(Cursor::new(Vec::new()), 1000, 0.02, 13).unwrap();
        render_pressure_pcm16(&mut build(), &mut stream, &CancelGate::new(), &mut [0.0; 13], 101).unwrap();
        let (sink, summary) = stream.finish().unwrap();
        assert_eq!(sink.into_inner(), expected_wav);
        assert_eq!(summary.clipped_samples, clips as u64);
        assert_eq!(summary.samples, 101);
    }
    output
}

// Analytic Cartesian mechanism, independent of native recursion, quintics,
// motion preparation and production retarded-root/pressure implementation.
fn rotate(t: f64, v: [f64; 3]) -> [f64; 3] {
    let (s, c) = (0.2+2.0*t+0.4*t*t).sin_cos();
    [c*v[0]-s*v[1], s*v[0]+c*v[1], v[2]]
}
fn point(t: f64, kind: usize) -> SourceKinematics {
    let a = ARMS[kind];
    let (r, dr, ddr) = if kind == 0 { (0.8+0.3*t+0.2*t*t, 0.3+0.4*t, 0.4) }
        else { (0.0, 0.0, 0.0) };
    let w = 2.0+0.8*t;
    let x = r+a[0];
    let p = rotate(t, [x, a[1], a[2]]);
    let v = rotate(t, [dr-w*a[1], w*x, 0.0]);
    let acc = rotate(t, [ddr-0.8*a[1]-w*w*x, 2.0*w*dr+0.8*x-w*w*a[1], 0.0]);
    SourceKinematics { position_m: [p[0]+0.7*t, p[1]-0.2*t, p[2]+0.1*t*t],
        velocity_m_s: [v[0]+0.7, v[1]-0.2, v[2]+0.2*t],
        acceleration_m_s2: [acc[0], acc[1], acc[2]+0.2] }
}
fn strength_world(t: f64, scale: f64) -> CompactSourceStrength {
    let s = source(t, scale);
    let f = s.force_on_fluid_body_n;
    let df = s.force_rate_on_fluid_body_n_s;
    let w = 2.0+0.8*t;
    CompactSourceStrength { volume_flow_m3_s: s.volume_flow_m3_s,
        volume_flow_rate_m3_s2: s.volume_flow_rate_m3_s2,
        force_on_fluid_n: rotate(t, f),
        force_rate_on_fluid_n_s: rotate(t, [df[0]-w*f[1], df[1]+w*f[0], df[2]]) }
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 { (0..3).map(|i| a[i]*b[i]).sum() }
fn norm(a: [f64; 3]) -> f64 { dot(a,a).sqrt() }
fn reference_pressure(time: f64) -> f64 {
    let receiver = point(time, 2).position_m;
    [1.0, -0.4].into_iter().enumerate().map(|(kind, scale)| {
        let mut lo = -0.1;
        let mut hi = 0.3;
        for _ in 0..96 {
            let tau = 0.5*(lo+hi);
            let source = point(tau, kind).position_m;
            let r = norm(core::array::from_fn(|i| receiver[i]-source[i]));
            if tau+r/300.0 < time { lo=tau; } else { hi=tau; }
        }
        let tau=0.5*(lo+hi);
        let motion=point(tau,kind);
        let s=strength_world(tau,scale);
        let offset=core::array::from_fn(|i|receiver[i]-motion.position_m[i]);
        let r=norm(offset);
        let n=offset.map(|x|x/r);
        let m=motion.velocity_m_s.map(|v|v/300.0);
        let mr=dot(n,m);
        let m2=dot(m,m);
        let ar=dot(n,motion.acceleration_m_s2);
        let k=1.0-mr;
        let q=1.2*(s.volume_flow_rate_m3_s2/(r*k*k)
            +s.volume_flow_m3_s*(r*ar/300.0+300.0*(mr-m2))/(r*r*k.powi(3)));
        let f=dot(n,s.force_rate_on_fluid_n_s)/(300.0*r*k*k)
            -dot(s.force_on_fluid_n,m)/(r*r*k*k)
            +dot(n,s.force_on_fluid_n)*(1.0-m2+r*ar/90_000.0)/(r*r*k.powi(3));
        (q+f)/(4.0*PI)
    }).sum()
}

#[test]
fn joint_driven_sources_and_microphone_reach_coherent_retarded_pressure_and_pcm() {
    let output=observe(80,true);
    for (i,value) in output.into_iter().enumerate() {
        let expected=reference_pressure(0.04+i as f64/1000.0);
        assert!((value-expected).abs()<2e-12, "{value:e} != {expected:e}");
    }
}

#[test]
fn refining_accepted_mechanical_events_reduces_pressure_error_without_changing_the_output_clock() {
    let coarse=observe(40,false);
    let fine=observe(80,false);
    let mut coarse_error=0.0_f64;
    let mut fine_error=0.0_f64;
    for i in 0..101 {
        let expected=reference_pressure(0.04+i as f64/1000.0);
        coarse_error=coarse_error.max((coarse[i]-expected).abs());
        fine_error=fine_error.max((fine[i]-expected).abs());
    }
    // Cubic source interpolation has third-order derivative error; its
    // derivative, not only its value, participates in acoustic radiation.
    assert!(coarse_error<2e-11 && fine_error<0.3*coarse_error,
        "coarse={coarse_error:e}, fine={fine_error:e}");
}
