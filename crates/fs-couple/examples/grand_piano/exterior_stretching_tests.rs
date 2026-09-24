//! The selected string constitutive must reach played mechanics AND exterior Pa.
use super::*;
use crate::{board, performance::Performance};

fn specification(courses: &[Course], nonlinear: bool) -> String {
    let mut text = format!("{}\n", linear::string_stretching::HEADER);
    for c in courses {
        if nonlinear { text.push_str(&format!("stretch,{},150000,0.2\n", c.midi)); }
        else { text.push_str(&format!("linear,{}\n", c.midi)); }
    }
    text
}
fn controls(courses: &[Course], text: &str) -> Controls {
    Controls::from_texts(courses, None, None, Some("estimated")).unwrap()
        .with_string_stretching(text, courses).unwrap()
}
fn same(a: &engine::Instrument, b: &engine::Instrument) {
    assert_eq!(a.bank.q, b.bank.q); assert_eq!(a.bank.v, b.bank.v);
    assert_eq!(a.energy_j(), b.energy_j());
    assert_eq!(a.radiation_energy_j(), b.radiation_energy_j());
    assert_eq!(a.accounting.input_work_j, b.accounting.input_work_j);
    assert_eq!(a.accounting.dissipated_j(), b.accounting.dissipated_j());
}

#[test]
fn complete_stretching_is_explicit_and_never_silently_used_as_linear_harmonic_data() {
    let parse = |args: &[&str]| Options::parse(&args.iter().map(|s| s.to_string()).collect::<Vec<_>>());
    let options = parse(&["--string-stretching", "strings.fsps", "score.mid", "--modes", "48"]).unwrap();
    assert_eq!(options.string_stretching.as_deref(), Some("strings.fsps"));
    assert_eq!(options.midi.as_deref(), Some("score.mid"));
    assert!(Options::harmonic(&["--string-stretching".into(), "strings.fsps".into()]).is_err());
    for args in [vec!["--string-stretching"], vec!["--string-stretching", ""],
        vec!["--string-stretching", "  "], vec!["--string-stretching", "--modes", "12"],
        vec!["--string-stretching", "a", "--string-stretching", "b"]] {
        assert!(parse(&args).is_err());
    }
    let courses = steinway_scale::courses().unwrap(); let courses = &courses[48..50];
    for text in [String::new(), format!("{}\n", linear::string_stretching::HEADER),
        specification(&courses[..1], true), specification(courses, true).replace("150000", "NaN"),
        specification(courses, true).replace("0.2", "0.31")] {
        assert!(Controls::from_texts(courses, None, None, None).unwrap()
            .with_string_stretching(&text, courses).is_err());
    }
    // Even an unplayed second key requires an explicit constitutive choice.
    assert!(Controls::from_texts(courses, None, None, None).unwrap()
        .instrument(courses.to_vec(), &board::demonstration(), &options).is_err());
    assert!(Controls::load(&Options { string_stretching: Some("missing-string-material.fsps".into()),
        ..Options::default() }, courses).is_err());
    assert!(controls(courses, &specification(courses, true))
        .with_string_stretching(&specification(courses, true), courses).is_err());
}

#[test]
fn explicit_linear_selection_preserves_jack_material_and_pedal_motion_bit_for_bit() {
    let courses = vec![steinway_scale::courses().unwrap()[48]];
    let options = Options { modes: 12, ..Options::default() };
    let mut selected = controls(&courses, &specification(&courses, false))
        .instrument(courses.clone(), &board::demonstration(), &options).unwrap();
    let mut original = Controls::from_texts(&courses, None, None, Some("estimated")).unwrap()
        .instrument(courses.clone(), &board::demonstration(), &options).unwrap();
    assert!(!selected.bank.has_string_stretching());
    for p in [&mut selected, &mut original] { p.jack_on(69, 70., 0.007).unwrap(); }
    for i in 0..2400 {
        if i == 1200 { for p in [&mut selected, &mut original] {
            p.set_sustain(0.4).unwrap(); p.note_off(69).unwrap();
        }}
        assert_eq!(selected.step().unwrap(), original.step().unwrap());
    }
    same(&selected, &original);
}

#[test]
fn nonlinear_selection_reaches_loaded_and_one_way_stereo_without_a_second_mechanical_clock() {
    for loaded in [false, true] {
        let (board_text, courses, obj, spec) = crate::tests::small_source_inputs();
        let material = specification(&courses, true);
        let options = Options { modes: 12, ..Options::default() };
        let selected = controls(&courses, &material);
        let direct = Controls::from_texts(&courses, None, None, Some("estimated")).unwrap();
        let linear = Controls::from_texts(&courses, None, None, Some("estimated")).unwrap();
        let mut scene = crate::prepare_controlled(&board_text, courses.clone(), &obj, spec, &options, selected).unwrap();
        let mut manual = crate::prepare_controlled(&board_text, courses.clone(), &obj, crate::tests::small_source_inputs().3, &options, direct).unwrap();
        let mut reference = crate::prepare_controlled(&board_text, courses.clone(), &obj, crate::tests::small_source_inputs().3, &options, linear).unwrap();
        manual.piano.configure_string_stretching(
            &linear::string_stretching::Specification::read(&material, &courses).unwrap()).unwrap();
        // Rest frequencies/boundary coordinates are NOT retuned by an amplitude law.
        assert_eq!(scene.piano.bank.modes.iter().map(|m| m.omega).collect::<Vec<_>>(),
            reference.piano.bank.modes.iter().map(|m| m.omega).collect::<Vec<_>>());
        let (baked, _, _) = crate::bake(&mut scene, loaded).unwrap();
        crate::bake(&mut manual, loaded).unwrap(); crate::bake(&mut reference, loaded).unwrap();
        let score = || Performance::read("sample,event,key,value\n0,note_on,69,1\n1200,sustain,0,0.4\n1200,note_off,69,0\n",
            &[69], 2400).unwrap();
        let output = crate::exterior_audio::render(&mut scene.piano, score(), 2400, &baked, 2.).unwrap();
        let mut program = score(); let mut linear_program = score();
        let mut maximum_extension = 0.0_f64;
        for sample in 0..2400 {
            program.dispatch(sample, &mut manual.piano).unwrap(); manual.piano.step().unwrap();
            linear_program.dispatch(sample, &mut reference.piano).unwrap(); reference.piano.step().unwrap();
            for i in 0..manual.piano.bank.strings.len() {
                maximum_extension = maximum_extension.max(manual.piano.bank.string_stretching_observation(i)
                    .unwrap().stretching_energy_j);
            }
        }
        same(&scene.piano, &manual.piano);
        assert!(maximum_extension > 0.);
        assert_ne!(scene.piano.bank.q, reference.piano.bank.q);
        assert!(scene.piano.accounting.felt_loss_j > 0. && scene.piano.accounting.damper_loss_j > 0.);
        assert!((scene.piano.energy_j() + scene.piano.accounting.dissipated_j()
            - scene.piano.accounting.input_work_j).abs() < 1e-7);
        assert!(output.peak_pa > 1e-14); assert_eq!(&output.wav[..4], b"RIFF");
        assert_eq!(u16::from_le_bytes([output.wav[22], output.wav[23]]), 2);
        let data = output.wav.windows(4).position(|w| w == b"data").unwrap() + 8;
        for frame in output.wav[data..].chunks_exact(4) { assert_eq!(&frame[..2], &frame[2..]); }
        assert!(options.report(&scene.piano).contains("nonlinear string extension: true"));
    }
}
