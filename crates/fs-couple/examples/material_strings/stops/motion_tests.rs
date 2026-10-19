//! File geometry/motion -> original material performer -> physical pressure/PCM.
use super::*;
use crate::stops::Profile;
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::performance::motion::{
    StopMotionKnot, StopMotionProgram, ObstacleMotionLimits,
};

const STATIC: &str = "frankensim-string-stops-v1\nstop support 0 0.43 positive fixed 0 5000 1 0.2 10 0.001\n";
const MOTION: &str = "motion 0.001 0.1 3\nknot 12 -0.000003\nknot 36 -0.000003\nknot 72 -0.000001\n";
fn program() -> StopMotionProgram {
    StopMotionProgram {
        knots: [(12, -3e-6), (36, -3e-6), (72, -1e-6)].into_iter()
            .map(|(sample, value)| StopMotionKnot { sample, offsets_m: vec![value] }).collect(),
        limits: ObstacleMotionLimits { maximum_abs_offset_m: 0.001, maximum_speed_m_s: 0.1 },
        max_knots: 128, max_coordinate_values: 4096,
    }
}
fn build(cx: &Cx<'_>, gate: &CancelGate, spec: input::Spec, profile: &str, rho: f64) -> Ready {
    spec.compile_with_profiles(cx, gate, None, Some(Profile::parse(profile.as_bytes()).unwrap()),
        |name| Ok(pack(name, if name == "body.fsmatpk" { 7800.0 } else { rho }, 200_000.0, None))).unwrap()
}

#[test]
fn file_motion_matches_direct_binding_and_original_pressure_pcm() {
    scoped(|cx, gate| {
        let text = format!("{STATIC}{MOTION}");
        let mut from_file = build(cx, gate, spec(), &text, 7800.0);
        assert_eq!(from_file.stop_profile.as_ref().unwrap().bytes, text.as_bytes());
        assert_eq!(from_file.performance.stop_motion_program(), Some(&program()));
        let mut manual = build(cx, gate, spec(), STATIC, 7800.0);
        manual.performance = manual.performance.with_stop_motion(program(), gate).unwrap();
        assert_eq!(manual.performance.preload_energy_j(), from_file.performance.preload_energy_j());
        let expected = pressure(&mut manual, 1);
        assert!(expected.iter().any(|p| p.abs() > 1e-12));
        assert_eq!(pressure(&mut from_file, 37), expected);
        assert_eq!(from_file.performance.bank().system().obstacle_offsets_m().unwrap(), &[-1e-6]);
        assert_eq!(from_file.performance.bank().system().last_frame().unwrap().obstacle_work_j, 0.0,
            "held final geometry must not keep supplying work");
        let mut wav = build(cx, gate, spec(), &text, 7800.0);
        let (bytes, summary) = wav.write(Cursor::new(Vec::new()), gate).unwrap();
        let (encoded, clipped) = fs_couple::pcm_wav::encode_pcm16_wav(&expected, wav.spec.rate, wav.spec.full_scale).unwrap();
        assert_eq!(bytes.into_inner(), encoded);
        assert_eq!(summary.clipped_samples, clipped as u64);
        assert_eq!(summary.samples, wav.spec.samples);
    });
}

#[test]
fn zero_trajectory_preserves_static_profiles_and_zero_motion_work() {
    scoped(|cx, gate| {
        let mut fixed = build(cx, gate, spec(), STATIC, 7800.0);
        let text = format!("{STATIC}motion 0.001 0.1 1\nknot 256 0\n");
        let mut zero = build(cx, gate, spec(), &text, 7800.0);
        assert!(fixed.performance.stop_motion_program().is_none());
        assert_eq!(fixed.performance.preload_energy_j(), zero.performance.preload_energy_j());
        assert_eq!(pressure(&mut fixed, 128), pressure(&mut zero, 37));
        assert_eq!(zero.performance.bank().system().obstacle_offsets_m().unwrap(), &[0.0]);
        assert_eq!(zero.performance.bank().system().last_frame().unwrap().obstacle_work_j, 0.0);
    });
}

#[test]
fn independent_normal_paths_keep_two_plane_material_maps_and_change_physical_output() {
    scoped(|cx, gate| {
        let input = include_str!("../polarized-pluck.material-strings").replace("48000 4096", "48000 256");
        let profile = "frankensim-string-stops-v1\n\
stop first 0 0.43 normal 1 0 fixed 0 5000 1 0.2 10 0.001\n\
stop second 0 0.31 normal 0 1 fixed 0.00002 3000 1 0.1 10 0.001\n\
motion 0.001 0.1 3\nknot 12 -0.000003 0\nknot 36 -0.000003 -0.000002\nknot 72 -0.000001 -0.000002\n";
        let make = |rho| build(cx, gate, input::Spec::parse(input.as_bytes()).unwrap(), profile, rho);
        let mut a = make(7800.0); let mut b = make(15600.0);
        assert_eq!(a.performance.bank().system().components()[1].modes(), b.performance.bank().system().components()[1].modes());
        let n = a.spec.strings[0].bending_modes;
        let x = &a.performance.bank().system().obstacle(0).unwrap().left.shapes;
        let y = &a.performance.bank().system().obstacle(1).unwrap().left.shapes;
        assert!(x[n..].iter().all(|v| *v == 0.0));
        assert!(y[..n].iter().chain(&y[2*n..]).all(|v| *v == 0.0));
        let mut first = vec![0.0; 12]; a.performance.block(&mut first).unwrap();
        assert_eq!(a.performance.bank().system().obstacle_offsets_m().unwrap(), &[-3e-6, 0.0]);
        let mut rest = vec![0.0; 244];
        for block in rest.chunks_mut(37) { a.performance.block(block).unwrap(); }
        first.extend(rest);
        assert_eq!(a.performance.bank().system().obstacle_offsets_m().unwrap(), &[-1e-6, -2e-6]);
        let heavy = pressure(&mut b, 128);
        assert!(first.iter().zip(heavy).any(|(a, b)| (a-b).abs() > 1e-12));
        assert_eq!(pressure(&mut make(7800.0), 1), first);
    });
}

#[test]
fn malformed_or_out_of_window_motion_refuses_before_material_io() {
    for tail in ["motion 0.001 0.1 0\n", "motion 0.001 0.1 129\n",
        "motion 0.001 0.1 1\nknot 12\n", "motion 0.001 0.1 1\nknot 12 0 0\n",
        "motion 0.001 0.1 1\nknot 12 NaN\n", "motion 0.001 0.1 1\nknot 12 0\nextra\n"] {
        assert!(Profile::parse(format!("{STATIC}{tail}").as_bytes()).is_err());
    }
    assert!(Profile::parse(b"frankensim-string-stops-v1\nmotion 0.001 0.1 1\nknot 12 0\n").is_err());
    scoped(|cx, gate| {
        for tail in ["motion 0.001 0.1 1\nknot 0 0\n", "motion 0.001 0.1 1\nknot 257 0\n",
            "motion 0.001 0.1 2\nknot 12 0\nknot 12 0\n", "motion 0.001 0.1 1\nknot 12 0.002\n",
            "motion 0.001 0.00001 1\nknot 12 -0.000003\n"] {
            let profile = Profile::parse(format!("{STATIC}{tail}").as_bytes()).unwrap();
            let mut calls = 0;
            let result = spec().compile_with_profiles(cx, gate, None, Some(profile), |name| {
                calls += 1; Ok(pack(name, 7800.0, 200_000.0, None))
            });
            assert!(result.is_err()); assert_eq!(calls, 0);
        }
        for (stops, samples, rate) in [(0, 256, 48000), (1, 0, 48000), (1, 256, 0)] {
            assert!(program().validate(stops, samples, rate, gate).is_err());
        }
    });
}
