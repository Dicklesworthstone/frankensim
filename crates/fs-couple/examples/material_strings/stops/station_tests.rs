//! Native geometry data -> original material/contact/DG performer -> real PCM.
use super::*;
use crate::stops::Profile;

fn input() -> input::Spec {
    // Keep the source example's two planes, torsion, bonded section, complete
    // cavity and causal fluid owner. A smaller retained basis bounds this test;
    // neither this basis nor the example is a measured fidelity certificate.
    let text = include_str!("../migrating-finger.material-strings")
        .replace("48000 4096", "48000 256").replace("20 4 1 0.8", "20 1 1 0.8");
    input::Spec::parse(text.as_bytes()).unwrap()
}
fn profile() -> String {
    include_str!("../migrating-finger.stops")
        .replace(" 480 ", " 32 ").replace(" 1920 ", " 128 ").replace(" 2400 ", " 160 ")
        .replace(" 2880 ", " 192 ").replace(" 3360 ", " 224 ")
        .replace("0.2 0.48", "0.2 0.4312").replace("0.2 0.46", "0.2 0.4304")
}
fn build(cx: &Cx<'_>, gate: &CancelGate, text: &str) -> Ready {
    input().compile_with_profiles(cx, gate, None, Some(Profile::parse(text.as_bytes()).unwrap()),
        |name| Ok(pack(name, 7800.0, 200_000.0, None))).unwrap()
}

#[test]
fn native_moving_finger_uses_original_physics_clock_and_pcm() {
    scoped(|cx, gate| {
        let declared = input::Spec::parse(include_bytes!("../migrating-finger.material-strings")).unwrap();
        assert_eq!(declared.strings[0].bending_modes, 4);
        Profile::parse(include_bytes!("../migrating-finger.stops")).unwrap().validate_performance(&declared, gate).unwrap();
        let text = profile(); let source = Profile::parse(text.as_bytes()).unwrap();
        let fixed = text.split_once("station-motion").unwrap().0;
        let mut from_file = build(cx, gate, &text);
        let mut direct = build(cx, gate, fixed);
        direct.performance = direct.performance.with_stop_station_motion(source.station_program().unwrap().clone(), gate).unwrap();
        assert_eq!(from_file.performance.stop_station_program(), source.station_program());
        assert_eq!(from_file.stop_profile.as_ref().unwrap().bytes, text.as_bytes());
        let expected = pressure(&mut direct, 1);
        assert!(expected.iter().any(|p| p.abs() > 1e-15));
        assert_eq!(pressure(&mut from_file, 37), expected);
        let bank = from_file.performance.bank();
        assert_eq!(bank.stop_station_fractions().unwrap(), &[0.2, 0.4304]);
        assert_eq!(bank.system().obstacle_offsets_m().unwrap(), &[0.00003, 0.00003]);
        assert_eq!(bank.source().stops()[1].location.station_fraction, 0.43);
        assert_eq!(bank.system().last_frame().unwrap().obstacle_work_j, 0.0);
        assert!(bank.system().last_frame().unwrap().nonlinear.is_some());
        assert_eq!(bank.source().bank().geometric_stretching().unwrap().potentials()[0].transverse_planes(), 2);
        assert!(bank.source().bank().cylinder_gas_loading().is_some());
        assert!(bank.source().bank().material_cavity().is_some());
        assert_eq!(bank.system().components().iter().map(|m| m.modes().len()).collect::<Vec<_>>(), vec![3, 2, 1, 12]);
        assert!(bank.system().components()[3].states().iter().any(|s| s.velocity_m_sqrt_kg_per_s != 0.0));
        let mut wav = build(cx, gate, &text);
        let (output, summary) = wav.write(Cursor::new(Vec::<u8>::new()), gate).unwrap();
        let encoded = fs_couple::pcm_wav::encode_pcm16_wav(&expected, wav.spec.rate, wav.spec.full_scale).unwrap();
        assert_eq!(output.into_inner(), encoded.0);
        assert_eq!(summary.clipped_samples, encoded.1 as u64);
        assert_eq!(summary.samples, 256);
    });
}

#[test]
fn constant_station_path_retains_exact_fixed_stop_pressure() {
    scoped(|cx, gate| {
        let moving = profile();
        let constant = moving.replace("0.4312", "0.43").replace("0.4304", "0.43");
        let fixed = moving.split_once("station-motion").unwrap().0;
        let mut original = build(cx, gate, fixed);
        let mut admitted = build(cx, gate, &constant);
        assert_eq!(pressure(&mut original, 128), pressure(&mut admitted, 37));
        assert!(admitted.performance.bank().system().last_obstacle_frames().unwrap()
            .iter().all(|c| c.attachment_work_j == 0.0));
    });
}

#[test]
fn incomplete_unphysical_and_moving_bow_normal_paths_refuse_before_pack_io() {
    let declarations = "frankensim-string-stops-v1\nstop finger 0 0.43 positive fixed 0 5000 1 0.2 10 0.001\n";
    for tail in ["station-motion 1 0.001 1048576 0\n", "station-motion 1 0.001 1048576 129\n",
        "station-motion 1 0.001 0 1\nstation-knot 12 0.43\n",
        "station-motion 1 0.001 1048577 1\nstation-knot 12 0.43\n",
        "station-motion 1 0.001 1048576 1\nstation-knot 12\n",
        "station-motion 1 0.001 1048576 1\nstation-knot 12 0.43 0.43\n",
        "station-motion 1 0.001 1048576 1\nstation-knot 12 NaN\n",
        "station-motion 1 0.001 1048576 1\nstation-knot 12 0.43\nextra\n"] {
        assert!(Profile::parse(format!("{declarations}{tail}").as_bytes()).is_err());
    }
    let text = profile();
    let tail = text.split_once("station-motion").unwrap().1;
    assert!(Profile::parse(format!("{text}station-motion{tail}").as_bytes()).is_err());
    // Both optional geometry blocks share the same declared stops and clock;
    // their ordering does not introduce another coordinate interpolation route.
    let (before_normal, normal_and_stations) = text.split_once("motion ").unwrap();
    let (normal, stations) = normal_and_stations.split_once("station-motion").unwrap();
    let reordered = format!("{before_normal}station-motion{stations}motion {normal}");
    assert_eq!(Profile::parse(reordered.as_bytes()).unwrap().station_program(), Profile::parse(text.as_bytes()).unwrap().station_program());
    scoped(|cx, gate| {
        for invalid in [
            text.replace("station-knot 32", "station-knot 0"),
            text.replace("station-knot 192", "station-knot 257"),
            text.replace("0.4312", "1"),
            text.replace("0.4312", "0.7"),
            text.replace("station-motion 1 0.00003", "station-motion 1 0.000000001"),
            text.replace("station-knot 32 0.2", "station-knot 32 0.20001"),
        ] {
            let source = Profile::parse(invalid.as_bytes()).unwrap();
            let mut calls = 0;
            let result = input().compile_with_profiles(cx, gate, None, Some(source), |name| {
                calls += 1; Ok(pack(name, 7800.0, 200_000.0, None))
            });
            assert!(result.is_err()); assert_eq!(calls, 0);
        }
    });
}
