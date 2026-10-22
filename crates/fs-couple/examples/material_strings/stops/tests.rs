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

const DIRECTED: &str = include_str!("../directional-pluck.material-strings");
const NORMAL: &str = include_str!("../directional-fingerboard.stops");

#[test]
fn independent_force_and_stop_directions_round_trip_to_pressure_and_pcm() {
    scoped(|cx, gate| {
        let make = |profile: &str| input::Spec::parse(DIRECTED.replace("48000 4096", "48000 128").as_bytes())
            .unwrap().compile_with_profiles(cx, gate, None, Some(Profile::parse(profile.as_bytes()).unwrap()),
                |name| Ok(pack(name, 7800.0, 200_000.0, None))).unwrap();
        let mut a = make(NORMAL);
        let n = a.spec.strings[0].bending_modes;
        let bank = a.performance.bank();
        let normal = &bank.system().obstacle(0).unwrap().left.shapes;
        let bow = &bank.system().port(0).unwrap().left.shapes;
        assert!(normal[..n].iter().all(|x| *x == 0.0));
        assert!(normal[n..2*n].iter().any(|x| *x != 0.0));
        assert!(normal[2*n..].iter().all(|x| *x == 0.0));
        assert!(bow[n..2*n].iter().all(|x| *x == 0.0));
        assert!(bow[2*n..].iter().any(|x| *x != 0.0));
        assert_eq!(a.stop_profile.as_ref().unwrap().bytes, NORMAL.as_bytes());
        let expected = pressure(&mut a, 37);
        assert!(expected.iter().any(|p| p.abs() > 1e-12));
        let mut b = make(NORMAL);
        let (bytes, _) = b.write(Cursor::new(Vec::new()), gate).unwrap();
        let (reference, _) = fs_couple::pcm_wav::encode_pcm16_wav(&expected, b.spec.rate, b.spec.full_scale).unwrap();
        assert_eq!(bytes.into_inner(), reference);
        assert_eq!(pressure(&mut make(NORMAL), 1), expected);
        // Geometry-only intervention: no material modes or observer gains changed.
        let other = pressure(&mut make(&NORMAL.replace("normal 0 1", "normal 1 0")), 37);
        assert!(expected.iter().zip(other).any(|(a,b)| (a-b).abs() > 1e-12));
    });
}

#[test]
fn malformed_or_unrepresented_directions_do_not_become_projected_forces() {
    for replacement in ["along 0 0", "along 2 0", "along NaN 1", "along 0", "along 0 1 7"] {
        assert!(input::Spec::parse(DIRECTED.replace("along 0 1", replacement).as_bytes()).is_err());
    }
    let one_plane = DIRECTED.replace("polarizations two 1 0 0.6 0.8", "polarizations single");
    assert!(input::Spec::parse(one_plane.as_bytes()).is_err());
    assert!(input::Spec::parse(one_plane.replace("along 0 1", "along -1 0").as_bytes()).is_ok());
    for replacement in ["normal 0 0", "normal 2 0", "normal 0 NaN", "normal 0", "normal 0 1 7"] {
        assert!(Profile::parse(NORMAL.replace("normal 0 1", replacement).as_bytes()).is_err());
    }
}

#[test]
fn explicit_receiver_normals_require_a_real_collinear_bridge_footprint() {
    scoped(|cx, gate| {
        let input = DIRECTED.replace("48000 4096", "48000 128");
        let construct = |input: &str, normal: &str| input::Spec::parse(input.as_bytes()).unwrap()
            .compile_with_profiles(cx, gate, None, Some(Profile::parse(normal.as_bytes()).unwrap()),
                |name| Ok(pack(name, 7800.0, 200_000.0, None)));
        let profile = NORMAL.replace("fixed", "bridge:0");
        assert!(construct(&input, &profile).is_err(), "oblique scalar receiver was silently repurposed");
        let aligned = input.replace("polarizations two 1 0 0.6 0.8", "polarizations two 1 0 0 1");
        let a = construct(&aligned, &profile).unwrap();
        let bridge = &a.performance.bank().source().bank().strings()[0].connections[0].body.shapes;
        assert_eq!(&a.performance.bank().system().obstacle(0).unwrap().right.shapes, bridge);
        let b = construct(&aligned, &profile.replace("normal 0 1", "normal 0 -1")).unwrap();
        let receiver = &b.performance.bank().system().obstacle(0).unwrap().right.shapes;
        for (x,y) in receiver.iter().zip(bridge) { assert_eq!(*x, -*y); }
    });
}


const BOW_SURFACE: &str = "frankensim-string-stops-v1\nbow-surface 0 bridge:0\n";

fn surface_input() -> String {
    include_str!("../pluck.material-strings").replace("48000 4096", "48000 128")
        .replace("initial\tbow\t0e0\t0e0", "initial\tbow\t1e-3\t1e-1")
}

fn surface_ready(cx: &Cx<'_>, gate: &CancelGate, profile: Option<&str>, body_rho: f64) -> Ready {
    input::Spec::parse(surface_input().as_bytes()).unwrap().compile_with_profiles(cx, gate, None,
        profile.map(|text| Profile::parse(text.as_bytes()).unwrap()),
        |name| Ok(pack(name, if name == "body.fsmatpk" { body_rho } else { 7800.0 }, 200_000.0, None))).unwrap()
}

#[test]
fn data_selected_bow_surface_returns_reaction_to_original_material_body_and_pressure() {
    scoped(|cx, gate| {
        let mut a = surface_ready(cx, gate, Some(BOW_SURFACE), 7800.0);
        let mut heavy = surface_ready(cx, gate, Some(BOW_SURFACE), 15600.0);
        let mut prescribed = surface_ready(cx, gate, None, 7800.0);
        let bank = a.performance.bank();
        assert_eq!(bank.system().obstacle_count(), 0, "surfaces-only profile adds no normal law");
        assert_eq!(a.stop_profile.as_ref().unwrap().bytes, BOW_SURFACE.as_bytes());
        let surface = &bank.source().body_surfaces()[0];
        let bridge = &bank.source().bank().strings()[0].connections[0].body;
        assert_eq!(surface.bow, 0);
        assert_eq!(&surface.body, bridge);
        let right = bank.system().port(0).unwrap().right.as_ref().unwrap();
        assert_eq!(right.component, bank.cards().len() + bridge.component);
        assert_eq!(right.shapes, bridge.shapes);
        for (light, dense) in right.shapes.iter().zip(
            &heavy.performance.bank().system().port(0).unwrap().right.as_ref().unwrap().shapes) {
            assert!((dense - light*0.5_f64.sqrt()).abs() < 1e-12*light.abs().max(1.0));
        }
        assert_eq!(bank.system().port(0).unwrap().left, heavy.performance.bank().system().port(0).unwrap().left);
        for (a, b) in bank.system().components().iter().zip(prescribed.performance.bank().system().components()) {
            assert_eq!(a.modes(), b.modes(), "binding must not retune the original body/string modes");
        }
        assert_eq!(bank.system().acceleration_pressure_weights(),
            prescribed.performance.bank().system().acceleration_pressure_weights());
        assert!(a.performance.preload_energy_j().unwrap() > 0.0);
        assert_eq!(a.performance.preload_energy_j(), prescribed.performance.preload_energy_j());
        let (mut out, mut dense_out) = (Vec::new(), Vec::new());
        let mut traction_difference = 0.0;
        let mut traction_magnitude = 0.0;
        for _ in 0..a.spec.samples {
            let (mut sample, mut dense_sample) = ([0.0], [0.0]);
            a.performance.block(&mut sample).unwrap();
            heavy.performance.block(&mut dense_sample).unwrap();
            out.push(sample[0]); dense_out.push(dense_sample[0]);
            let frame = a.performance.bank().system().last_frame().unwrap();
            let dense = heavy.performance.bank().system().last_frame().unwrap();
            traction_difference += (frame.ports[0].traction_n - dense.ports[0].traction_n).abs();
            traction_magnitude += frame.ports[0].traction_n.abs();
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
            assert!(dense.energy_residual_j.abs() <= dense.energy_tolerance_j);
        }
        assert!(out.iter().all(|p| p.is_finite()) && out.iter().any(|p| p.abs() > 1e-12));
        assert!(out.iter().zip(dense_out).any(|(a,b)| (a-b).abs() > 1e-12));
        assert!(traction_difference > 1e-12 && traction_magnitude > 1e-12);
        assert!(out.iter().zip(pressure(&mut prescribed, 37)).any(|(a,b)| (a-b).abs() > 1e-12));
        assert_eq!(out, pressure(&mut surface_ready(cx, gate, Some(BOW_SURFACE), 7800.0), 37));
        let combined = format!("{PROFILE}bow-surface 0 bridge:0\n");
        let mut stopped = surface_ready(cx, gate, Some(&combined), 7800.0);
        assert_eq!(stopped.performance.bank().system().obstacle_count(), 1);
        assert_eq!(stopped.performance.bank().source().body_surfaces().len(), 1);
        assert!(pressure(&mut stopped, 37).iter().any(|p| p.abs() > 1e-12));
    });
}

#[test]
fn bow_surface_uses_the_signed_collinear_receiver_and_refuses_unsupported_data_before_io() {
    for row in ["bow-surface 0", "bow-surface 0 fixed", "bow-surface 0 bridge:no",
        "bow-surface 16 bridge:0", "bow-surface 0 bridge:0 extra",
        "bow-surface 0 bridge:0\nbow-surface 0 bridge:0"] {
        assert!(Profile::parse(format!("frankensim-string-stops-v1\n{row}\n").as_bytes()).is_err());
    }
    // Distinct physical contacts may share the same receiving footprint.
    assert!(Profile::parse(format!("{BOW_SURFACE}bow-surface 1 bridge:0\n").as_bytes()).is_ok());
    assert!(Profile::parse(format!("{BOW_SURFACE}motion 0.001 1 1\nknot 1 0\n").as_bytes()).is_err());
    assert!(Profile::parse(format!("{PROFILE}motion 0.001 1 1\nknot 1 0\nbow-surface 0 bridge:0\n").as_bytes()).is_err());
    scoped(|cx, gate| {
        for (input, profile) in [
            (surface_input(), BOW_SURFACE.replace("surface 0", "surface 1")),
            (surface_input(), BOW_SURFACE.replace("bridge:0", "bridge:1")),
            (surface_input().replace("friction stribeck 0.8 0.4 0.1 0", "friction stribeck 0.8 0.4 0.1 0\nbow_ribbon 0.01 3"), BOW_SURFACE.into()),
            (DIRECTED.replace("48000 4096", "48000 128"), BOW_SURFACE.into()),
        ] {
            let mut calls = 0;
            let result = input::Spec::parse(input.as_bytes()).unwrap().compile_with_profiles(cx, gate, None,
                Some(Profile::parse(profile.as_bytes()).unwrap()), |name| {
                    calls += 1; Ok(pack(name, 7800.0, 200_000.0, None))
                });
            assert!(result.is_err()); assert_eq!(calls, 0);
        }
        for (direction, sign) in [("1 0", 1.0), ("-1 0", -1.0)] {
            let text = DIRECTED.replace("48000 4096", "48000 128")
                .replace("polarizations two 1 0 0.6 0.8", &format!("polarizations two 1 0 {direction}"));
            let a = input::Spec::parse(text.as_bytes()).unwrap().compile_with_profiles(cx, gate, None,
                Some(Profile::parse(BOW_SURFACE.as_bytes()).unwrap()),
                |name| Ok(pack(name, 7800.0, 200_000.0, None))).unwrap();
            let bank = a.performance.bank();
            let bridge = &bank.source().bank().strings()[0].connections[0].body.shapes;
            for (actual, original) in bank.system().port(0).unwrap().right.as_ref().unwrap().shapes.iter().zip(bridge) {
                assert_eq!(*actual, sign*original);
            }
        }
    });
}
