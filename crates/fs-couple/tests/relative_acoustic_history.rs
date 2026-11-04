//! Moving local source coordinates reach the original causal pressure renderer.
use fs_aeroac::curle3d::moving::{MovingCompactLimits, MovingCompactObserver};
use fs_aeroac::curle3d::moving::trajectory::{PreparedSourceTrajectory, TrajectoryBudget};
use fs_couple::pcm_wav::observation::moving::{
    CompactSourceHistory, MovingCompactRenderer, MovingPressureBudget, MovingPressureConfig,
};
use fs_couple::pcm_wav::observation::moving::rigid::{
    BodyFrameCompactSourceStrength, RigidAcousticAttachment, RigidAcousticHistory,
    RigidHistoryBudget, RelativeAcousticMotion,
};
use fs_couple::pcm_wav::observation::moving::rigid::articulated::{
    LinkAcousticUpdate, append_link_frame_with_relative_motion,
};
use fs_ga::{Mat3, Se3, Vec3 as GVec, Wrench as GWrench};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_mbd::articulated::{
    ArticulatedModel, BaseState, JointModel, Link, PreparedLinkPointMotion, SpatialInertia,
    forward_kinematics, inverse_dynamics,
};
use fs_mbd::{Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState,
    UnitQuaternion, Vec3, Wrench};

const OMEGA: f64 = 5.0;
fn link(time: f64) -> PreparedLinkPointMotion {
    let inertia = SpatialInertia::new(2.0, GVec::default(), Mat3 {
        m: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
    }).unwrap();
    let model = ArticulatedModel::new(vec![Link::new("carrier", None, Se3::identity(),
        JointModel::revolute(GVec::new(0.0, 0.0, 1.0), None).unwrap(), inertia)]).unwrap();
    let base = BaseState::stationary(Se3::identity());
    let kin = forward_kinematics(&model, base, &[OMEGA*time], &[OMEGA]).unwrap();
    inverse_dynamics(&model, base, &[OMEGA*time], &[OMEGA], &[0.0], GVec::default(),
        &[GWrench::default()]).unwrap().prepare_point_motion(&kin, 0, GVec::default()).unwrap()
}
fn body(time: f64) -> fs_mbd::PreparedRigidPointMotion {
    let properties = MassProperties::new(2.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap();
    let state = RigidBodyState::new(Pose::new(Vec3::ZERO,
        UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), OMEGA*time).unwrap()).unwrap(),
        Vec3::ZERO, Vec3::new(0.0, 0.0, OMEGA)).unwrap();
    RigidBodyIntegrator::new(Gravity::ZERO).prepare_point_motion(state, properties, Wrench::ZERO).unwrap()
}
fn relative(t: f64) -> RelativeAcousticMotion {
    RelativeAcousticMotion {
        displacement_body_m: Vec3::new(0.4*t + 0.1*t*t, 0.0, 0.0),
        velocity_body_m_s: Vec3::new(0.4 + 0.2*t, 0.0, 0.0),
        acceleration_body_m_s2: Vec3::new(0.2, 0.0, 0.0),
    }
}
fn strength() -> BodyFrameCompactSourceStrength {
    BodyFrameCompactSourceStrength { force_on_fluid_body_n: [2.0, 1.0, 0.0], ..Default::default() }
}
fn history(capacity: usize) -> RigidAcousticHistory {
    RigidAcousticHistory::new(RigidAcousticAttachment::new(Vec3::new(0.5, 0.0, 0.0)).unwrap(),
        RigidHistoryBudget { maximum_samples: capacity, maximum_payload_bytes: 1_048_576 }).unwrap()
}
fn analytic(t: f64) -> ([f64; 3], [f64; 3], [f64; 3], [f64; 3], [f64; 3]) {
    let r = 0.5 + 0.4*t + 0.1*t*t;
    let d = 0.4 + 0.2*t;
    let (s, c) = (OMEGA*t).sin_cos();
    ([r*c, r*s, 0.0], [d*c-OMEGA*r*s, d*s+OMEGA*r*c, 0.0],
     [(0.2-OMEGA*OMEGA*r)*c-2.0*OMEGA*d*s,
      (0.2-OMEGA*OMEGA*r)*s+2.0*OMEGA*d*c, 0.0],
     [2.0*c-s, 2.0*s+c, 0.0], [-OMEGA*(2.0*s+c), OMEGA*(2.0*c-s), 0.0])
}
fn close(a: [f64; 3], b: [f64; 3]) {
    for (x, y) in a.into_iter().zip(b) { assert!((x-y).abs() < 2e-12, "{x:e} != {y:e}"); }
}

#[test]
fn native_carriers_export_consistent_migrating_motion_and_force_derivatives() {
    let mut h = history(3);
    for t in [-0.1, 0.0, 0.2] {
        let carrier = body(t);
        let articulated = link(t);
        let a = h.attachment();
        let x = a.sample_relative(t, &carrier, relative(t), strength()).unwrap();
        let y = a.sample_link_relative(t, &articulated, relative(t), strength()).unwrap();
        let (p,v,acc,f,df) = analytic(t);
        for sample in [x,y] {
            close(sample.kinematics.kinematics.position_m, p);
            close(sample.kinematics.kinematics.velocity_m_s, v);
            close(sample.kinematics.kinematics.acceleration_m_s2, acc);
            close(sample.source.strength.force_on_fluid_n, f);
            close(sample.source.strength.force_rate_on_fluid_n_s, df);
        }
        assert_eq!(a.sample_relative(t, &carrier, Default::default(), strength()).unwrap(),
            a.sample(t, &carrier, strength()).unwrap());
        assert_eq!(a.sample_link_relative(t, &articulated, Default::default(), strength()).unwrap(),
            a.sample_link(t, &articulated, strength()).unwrap());
        h.append_relative(t, &carrier, relative(t), strength(), &mut || false).unwrap();
        assert_eq!(h.kinematics().last().unwrap(), &x.kinematics);
    }
}

fn dot(a: [f64;3], b: [f64;3]) -> f64 { (0..3).map(|i| a[i]*b[i]).sum() }
fn reference_pressure(time: f64) -> f64 {
    let observer = [8.0,3.0,1.0];
    let mut lo = -0.1; let mut hi = 0.3;
    for _ in 0..80 {
        let tau = 0.5*(lo+hi);
        let p = analytic(tau).0;
        let r = core::array::from_fn(|i| observer[i]-p[i]);
        if tau+dot(r,r).sqrt()/300.0 < time { lo=tau; } else { hi=tau; }
    }
    let (p,v,a,f,df) = analytic(0.5*(lo+hi));
    let r: [f64;3] = core::array::from_fn(|i| observer[i]-p[i]);
    let distance=dot(r,r).sqrt(); let n=r.map(|v| v/distance); let m=v.map(|v| v/300.0);
    let k=1.0-dot(n,m);
    (dot(n,df)/(300.0*distance*k*k)-dot(f,m)/(distance*distance*k*k)
        +dot(n,f)*(1.0-dot(m,m)+distance*dot(n,a)/90_000.0)/(distance*distance*k*k*k))
        /(4.0*std::f64::consts::PI)
}

#[test]
fn migrating_source_reaches_native_retarded_pressure_and_partition_independent_callbacks() {
    let mut h = history(81);
    let pointer = h.kinematics().as_ptr();
    for i in 0..=80 {
        let t = -0.1 + 0.005*f64::from(i);
        h.append_link_relative(t, &link(t), relative(t), strength(), &mut || false).unwrap();
    }
    assert_eq!(pointer, h.kinematics().as_ptr());
    let observer = MovingCompactObserver::new(FluidAcousticState::try_new(1.2,300.0).unwrap(),
        [8.0,3.0,1.0], MovingCompactLimits { maximum_mach: 0.2, minimum_distance_m: 0.1 }).unwrap();
    let trajectory = PreparedSourceTrajectory::new(observer, h.kinematics(), TrajectoryBudget {
        maximum_segments: 128, maximum_payload_bytes: 1_048_576,
        maximum_bound_nodes: 16_384, maximum_subdivision_depth: 20, maximum_root_iterations: 96,
    }, &mut || false).unwrap();
    let tracks = [CompactSourceHistory { motion: &trajectory, samples: h.sources(), source_bandwidth_hz: 2.0 }];
    let config = MovingPressureConfig { first_receiver_time_s: 0.02, sample_rate_hz: 1000,
        samples: 81, max_block: 81, receiver_bandwidth_hz: 10.0, maximum_time_residual_s: 1e-12 };
    let budget = MovingPressureBudget { maximum_sources: 1, maximum_source_samples: 128,
        maximum_output_samples: 128, maximum_payload_bytes: 1_048_576, maximum_source_receiver_pairs: 128 };
    let build = || MovingCompactRenderer::new(observer, &tracks, config, budget, &mut || false).unwrap();
    let mut renderer = build(); let mut output = [0.0;81]; renderer.block(&mut output).unwrap();
    for (i, &p) in output.iter().enumerate() {
        let t=renderer.receiver_time_at(i as u64).unwrap();
        assert!((p-reference_pressure(t)).abs() < 1e-9);
    }
    let mut partitioned=build(); let mut small=[0.0;81];
    for block in small.chunks_mut(7) { partitioned.block(block).unwrap(); }
    assert_eq!(output.map(f64::to_bits),small.map(f64::to_bits));
    let mut retry=build(); let mut pending=[123.0;81]; let mut polls=0;
    assert!(retry.block_controlled(&mut pending,&mut || {polls+=1;polls==8}).is_err());
    assert_eq!(pending,[123.0;81]); assert_eq!(retry.samples_rendered(),0);
    retry.block(&mut pending).unwrap();
    assert_eq!(output.map(f64::to_bits),pending.map(f64::to_bits));
}

#[test]
fn relative_frames_keep_all_histories_unchanged_on_late_failure_and_cancellation() {
    let carrier=link(0.0); let mut a=history(2); let mut b=history(2); let mut scratch=[None;2];
    let zero=RelativeAcousticMotion::default();
    let invalid=RelativeAcousticMotion { velocity_body_m_s: Vec3::new(f64::NAN,0.0,0.0), ..zero };
    assert!(append_link_frame_with_relative_motion(0.0,&mut updates(&mut a,&mut b,&carrier),
        &[relative(0.0),invalid],&mut scratch,&mut || false).is_err());
    assert!(a.is_empty() && b.is_empty()); assert!(scratch[0].is_some());
    assert!(append_link_frame_with_relative_motion(0.0,&mut updates(&mut a,&mut b,&carrier),
        &[zero],&mut scratch,&mut || false).is_err());
    let mut polls=0;
    assert!(append_link_frame_with_relative_motion(0.0,&mut updates(&mut a,&mut b,&carrier),
        &[relative(0.0),zero],&mut scratch,&mut || {polls+=1;polls==4}).is_err());
    assert!(a.is_empty() && b.is_empty());
    append_link_frame_with_relative_motion(0.0,&mut updates(&mut a,&mut b,&carrier),
        &[zero,relative(0.0)],&mut scratch,&mut || false).unwrap();
    assert_eq!(a.kinematics()[0],a.attachment().kinematics_link(0.0,&carrier).unwrap());
    assert_eq!(b.kinematics()[0],b.attachment().kinematics_link_relative(0.0,&carrier,relative(0.0)).unwrap());
    assert_eq!(a.len(),1); assert_eq!(b.len(),1); assert!(scratch.iter().all(Option::is_none));
}

fn updates<'a>(a: &'a mut RigidAcousticHistory, b: &'a mut RigidAcousticHistory,
    carrier: &'a PreparedLinkPointMotion) -> [LinkAcousticUpdate<'a>; 2] {
    [LinkAcousticUpdate { history: a, motion: carrier, strength: strength() },
     LinkAcousticUpdate { history: b, motion: carrier, strength: strength() }]
}
