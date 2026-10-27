//! G0/G3/G5: native source -> continuing material bank -> physical PCM, with
//! independent direct-library replay and callback-size equivalence.
use super::*;

const PATH: &str = "bow_station_motion 2 0.00003 1048576 2\n\
    bow_station_knot 64 0.202\n\
    bow_station_knot 128 0.201\n";

fn source(ribbon: bool, motion: &str) -> String {
    let mut text = String::from_utf8(spec().bytes).unwrap();
    let track = GestureTrack {
        id: "bow".into(),
        target: GestureTarget::BowStroke { string: 0 },
        initial: GestureValue::Bow {
            velocity_m_per_s: 0.02,
            normal_force_n: 0.001,
            station: 0.2,
        },
        events: vec![],
    };
    let schedule = GestureSchedule::try_new(1000, vec![track]).unwrap();
    text.truncate(text.find("\nschedule\n").unwrap());
    text.push_str("\nschedule\n");
    text.push_str(std::str::from_utf8(&schedule.to_canonical_bytes()).unwrap());
    let ribbon = if ribbon { "bow_ribbon 0.01 3\n" } else { "" };
    text.replace("forces preload", &format!("{ribbon}{motion}forces preload"))
}

#[test]
fn native_bow_position_path_matches_the_library_and_streamed_pcm() {
    scoped(|cx, gate| {
        for ribbon in [false, true] {
            let moving = input::Spec::parse(source(ribbon, PATH).as_bytes()).unwrap();
            let program = moving.bow_station_program(gate).unwrap().unwrap();
            let mut from_file = compile(cx, gate, moving, 7800.0);
            let fixed = input::Spec::parse(source(ribbon, "").as_bytes()).unwrap();
            let mut direct = compile(cx, gate, fixed.clone(), 7800.0);
            let original_modes = direct.performance.bank().system().components()[0]
                .modes()
                .to_vec();
            direct.performance = direct
                .performance
                .with_bow_station_motion(program.clone(), gate)
                .unwrap();
            let expected = pressure(&mut direct, 1);
            assert!(expected.iter().any(|value| value.abs() > 1e-15));
            assert_eq!(pressure(&mut from_file, 37), expected);
            let bank = from_file.performance.bank();
            assert_eq!(
                bank.bow_station_fractions().unwrap(),
                program.knots.last().unwrap().station_fractions
            );
            assert_eq!(bank.source().bows()[0].location.station_fraction, 0.2);
            assert_eq!(bank.system().components()[0].modes(), original_modes);
            assert!(bank.source().bank().material_cavity().is_some());
            let mut held = compile(cx, gate, fixed, 7800.0);
            assert_ne!(
                pressure(&mut held, 37),
                expected,
                "moving contact must reach physical pressure"
            );
            let mut wav = compile(
                cx,
                gate,
                input::Spec::parse(source(ribbon, PATH).as_bytes()).unwrap(),
                7800.0,
            );
            let (output, summary) = wav.write(Cursor::new(Vec::<u8>::new()), gate).unwrap();
            let encoded =
                fs_couple::pcm_wav::encode_pcm16_wav(&expected, wav.spec.rate, wav.spec.full_scale)
                    .unwrap();
            assert_eq!(output.into_inner(), encoded.0);
            assert_eq!(summary.samples, 256);
            assert_eq!(summary.clipped_samples, encoded.1 as u64);
        }
    });
}

#[test]
fn native_constant_bow_position_is_the_original_hold() {
    scoped(|cx, gate| {
        let hold = PATH.replace("0.202", "0.2").replace("0.201", "0.2");
        for ribbon in [false, true] {
            let mut original = compile(
                cx,
                gate,
                input::Spec::parse(source(ribbon, "").as_bytes()).unwrap(),
                7800.0,
            );
            let mut moving = compile(
                cx,
                gate,
                input::Spec::parse(source(ribbon, &hold).as_bytes()).unwrap(),
                7800.0,
            );
            assert_eq!(pressure(&mut moving, 17), pressure(&mut original, 128));
        }
    });
}

#[test]
fn malformed_or_unphysical_bow_paths_refuse_at_source_admission() {
    for bad in [
        "bow_station_motion 2 0.00003 1048576 0\n".to_string(),
        PATH.replace("1048576 2", "1048576 129"),
        PATH.replace("1048576 2", "0 2"),
        PATH.replace("1048576 2", "1048577 2"),
        PATH.replace("motion 2 ", "motion NaN "),
        PATH.replace("0.00003", "0.000000001"),
        PATH.replace("knot 64", "knot 0"),
        PATH.replace("knot 128", "knot 64"),
        PATH.replace("knot 128", "knot 257"),
        PATH.replace("64 0.202", "64"),
        PATH.replace("64 0.202", "64 0.202 0.4"),
        PATH.replace("0.202", "NaN"),
        PATH.replace("0.202", "1"),
        PATH.replace("0.202", "0.4"),
        format!("{PATH}{PATH}"),
    ] {
        assert!(
            input::Spec::parse(source(false, &bad).as_bytes()).is_err(),
            "admitted {bad}"
        );
    }
    // The centre is inside the string, but the complete finite ribbon is not.
    let crossing = "bow_station_motion 100000 1 1048576 1\nbow_station_knot 64 0.995\n";
    assert!(input::Spec::parse(source(true, crossing).as_bytes()).is_err());
    input::Spec::parse(include_bytes!("../migrating-bow.material-strings")).unwrap();
}

#[test]
fn bow_geometry_is_rechecked_before_material_io() {
    scoped(|cx, gate| {
        let original = input::Spec::parse(source(false, PATH).as_bytes()).unwrap();
        for thermal in [false, true] {
            let mut changed = original.clone();
            if thermal {
                changed.contact_thermal = true;
            } else {
                changed.samples = 32;
            }
            let mut reads = 0;
            let result = changed.compile(cx, gate, |name| {
                reads += 1;
                Ok(pack(name, 7800.0, 200_000.0, None))
            });
            assert!(result.is_err());
            assert_eq!(reads, 0);
        }
    });
}
