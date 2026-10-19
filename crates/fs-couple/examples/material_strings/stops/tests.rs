//! Synthetic profile -> actual sources -> joint contact/preload -> pressure/PCM.
use super::*;
use crate::stops::{Profile, profile_paths};

const PROFILE: &str = "frankensim-string-stops-v1\nstop support 0 0.43 positive fixed 0 5000 1 0.2 10 0.001\n";
fn stopped(cx: &Cx<'_>, gate: &CancelGate, text: &str, rho: f64) -> Ready {
    spec().compile_with_profiles(cx, gate, None, Some(Profile::parse(text.as_bytes()).unwrap()),
        |name| Ok(pack(name, if name == "wire.fsmatpk" { rho } else { 7800.0 }, 200_000.0, None))).unwrap()
}

#[test]
fn profile_drives_contact_preload_release_and_original_pressure_pcm() {
    scoped(|cx, gate| {
        let mut a = stopped(cx, gate, PROFILE, 7800.0);
        assert_eq!(a.stop_profile.as_ref().unwrap().bytes, PROFILE.as_bytes());
        assert_eq!(a.performance.bank().system().obstacle_count(), 1);
        assert!(a.performance.preload_energy_j().unwrap() > 0.0);
        assert!(a.performance.bank().source().bank().material_cavity().is_some());
        assert!(a.performance.bank().system().last_frame().is_none());
        let normal = a.performance.bank().system().obstacle(0).unwrap();
        assert!(normal.right.shapes.iter().all(|x| *x == 0.0));
        assert!(normal.left.shapes[a.spec.strings[0].bending_modes..].iter().all(|x| *x == 0.0));
        let expected = pressure(&mut a, 37);
        assert!(expected.iter().any(|p| p.abs() > 1e-12));
        let mut b = stopped(cx, gate, PROFILE, 7800.0);
        let (bytes, summary) = b.write(Cursor::new(Vec::new()), gate).unwrap();
        let (reference, clips) = fs_couple::pcm_wav::encode_pcm16_wav(&expected, b.spec.rate, b.spec.full_scale).unwrap();
        assert_eq!(bytes.into_inner(), reference);
        assert_eq!(summary.clipped_samples, clips as u64);
        assert_eq!(summary.samples, b.spec.samples);
        assert_eq!(pressure(&mut stopped(cx, gate, PROFILE, 7800.0), 1), expected);
    });
}

#[test]
fn data_only_stop_and_wire_changes_reach_the_same_body_output() {
    scoped(|cx, gate| {
        let mut a = stopped(cx, gate, PROFILE, 7800.0);
        let mut b = stopped(cx, gate, &PROFILE.replace("5000", "50000"), 7800.0);
        for (ma, mb) in a.performance.bank().system().components().iter().zip(b.performance.bank().system().components()) {
            assert_eq!(ma.modes(), mb.modes(), "stop data must not retune the material basis");
        }
        let pa = pressure(&mut a, 37); let pb = pressure(&mut b, 128);
        assert!(pa.iter().zip(&pb).any(|(a,b)| (a-b).abs() > 1e-12));
        let mut c = stopped(cx, gate, PROFILE, 15600.0);
        let n = &a.performance.bank().system().obstacle(0).unwrap().left.shapes;
        let heavy = &c.performance.bank().system().obstacle(0).unwrap().left.shapes;
        for (x,y) in n[..a.spec.strings[0].bending_modes].iter().zip(heavy) {
            assert!((y/x - 0.5_f64.sqrt()).abs() < 1e-12);
        }
        assert_eq!(a.performance.bank().system().components()[1].modes(), c.performance.bank().system().components()[1].modes());
        let pc = pressure(&mut c, 128);
        assert!(pa.iter().zip(pc).any(|(a,c)| (a-c).abs() > 1e-12));
    });
}

#[test]
fn receiver_is_the_actual_retained_bridge_footprint_not_a_new_body() {
    scoped(|cx, gate| {
        let text = PROFILE.replace("fixed", "bridge:0");
        let a = stopped(cx, gate, &text, 7800.0);
        let bank = a.performance.bank();
        let bridge = &bank.source().bank().strings()[0].connections[0].body;
        let normal = bank.system().obstacle(0).unwrap();
        assert_eq!(normal.right.component, bank.cards().len() + bridge.component);
        assert_eq!(normal.right.shapes, bridge.shapes);
        assert_eq!(bank.source().stops()[0].body.as_ref().unwrap(), bridge);
        assert_eq!(bank.system().components().len(), compile(cx, gate, spec(), 7800.0).performance.bank().system().components().len());
    });
}

#[test]
fn legacy_output_and_explicit_profile_options_keep_their_meanings() {
    scoped(|cx, gate| {
        let mut a = compile(cx, gate, spec(), 7800.0);
        let mut b = spec().compile_with_profiles(cx, gate, None, None,
            |name| Ok(pack(name, 7800.0, 200_000.0, None))).unwrap();
        assert!(b.stop_profile.is_none());
        assert_eq!(pressure(&mut a, 128), pressure(&mut b, 37));
    });
    let args: Vec<std::ffi::OsString> = ["--stop-profile", "a", "--body-profile", "b"].map(Into::into).to_vec();
    assert_eq!(profile_paths(&args).unwrap(), (Some(Path::new("b")), Some(Path::new("a"))));
    for words in [vec!["--stop-profile"], vec!["--unknown", "a"], vec!["--stop-profile", "a", "--stop-profile", "b"]] {
        let args: Vec<std::ffi::OsString> = words.into_iter().map(Into::into).collect();
        assert!(profile_paths(&args).is_err());
    }
}

#[test]
fn malformed_profiles_and_unknown_sources_refuse_before_material_loading() {
    for (from, to) in [("0.43", "1.0"), ("positive", "sideways"), ("fixed", "bridge:no"),
        ("0 5000", "NaN 5000"), ("5000 1", "5000 0.5"), ("0.2 10", "-0.2 10"),
        ("10 0.001", "0 0.001"), ("support", "bad/name")] {
        assert!(Profile::parse(PROFILE.replace(from, to).as_bytes()).is_err(), "accepted {to}");
    }
    assert!(Profile::parse(format!("{PROFILE}{}", PROFILE.lines().nth(1).unwrap()).as_bytes()).is_err());
    assert!(Profile::parse(&vec![b' '; crate::stops::MAX_BYTES + 1]).is_err());
    scoped(|cx, gate| {
        for text in [PROFILE.replace("support 0", "support 9"), PROFILE.replace("fixed", "bridge:9")] {
            let profile = Profile::parse(text.as_bytes()).unwrap();
            let mut calls = 0;
            let result = spec().compile_with_profiles(cx, gate, None, Some(profile), |name| {
                calls += 1; Ok(pack(name, 7800.0, 200_000.0, None))
            });
            assert!(result.is_err()); assert_eq!(calls, 0);
        }
    });
}

#[test]
fn stopped_layered_two_plane_plucks_keep_the_current_physical_basis() {
    scoped(|cx, gate| {
        let input = include_str!("../polarized-pluck.material-strings")
            .replace("48000 4096", "48000 128");
        let make = || input::Spec::parse(input.as_bytes()).unwrap().compile_with_profiles(
            cx, gate, None, Some(Profile::parse(PROFILE.as_bytes()).unwrap()),
            |name| Ok(pack(name, 7800.0, 200_000.0, None))).unwrap();
        let mut a = make();
        let n = a.spec.strings[0].bending_modes;
        let normal = &a.performance.bank().system().obstacle(0).unwrap().left.shapes;
        assert_eq!(normal.len(), 2*n + a.spec.strings[0].twist_modes);
        assert!(normal[2*n..].iter().all(|x| *x == 0.0));
        // The fixture declares one physical transverse direction (0.6, 0.8).
        // Contact follows that direction in both mass-normalized bending planes.
        for (first, second) in normal[..n].iter().zip(&normal[n..2*n]) {
            assert!((second/first - 0.8/0.6).abs() < 1e-12);
        }
        assert!(a.performance.bank().system().port(0).unwrap().left.shapes[2*n..]
            .iter().any(|x| *x != 0.0), "bow surface must retain material torsion");
        assert!(a.performance.preload_energy_j().unwrap() > 0.0);
        let expected = pressure(&mut a, 1);
        assert!(expected.iter().any(|p| p.abs() > 1e-12));
        assert_eq!(pressure(&mut make(), 37), expected);
    });
}
