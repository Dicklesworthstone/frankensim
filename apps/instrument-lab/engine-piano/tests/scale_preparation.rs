//! G0/G3: preparation follows the offline physical owner, then actual live
//! mechanics consumes those courses. These are not measured-instrument tests.
use fs_live_piano::{Config, Description, LivePiano, PianoAudio, preparation::{self, Tuning},
    refmods::{geometry::{self, Course}, linear::BoardMode, stretch_tuning}};

fn description() -> Description {
    let mut bridge = [0.0; 88];
    bridge.fill(0.08);
    Description {
        courses: geometry::demonstration_scale().unwrap(),
        board: vec![BoardMode { frequency_hz: 180.0, damping_ratio: 0.015, bridge, volume: 0.01 }],
        damping: Some(vec![2.0 * 0.015 * std::f64::consts::TAU * 180.0]),
        surface: vec![fs_live_piano::refmods::board_geometry::SurfaceSample {
            position_m: [0.0, 0.0, 0.0], area_m2: 0.01, mode_shape: vec![1.0],
        }],
        board_provenance: "manufactured one-coordinate board".into(),
        board_mass_kg: 1.0, board_area_m2: 0.01, exact_low_modes: 1,
    }
}
fn tuning() -> Tuning {
    Tuning::new(442.0, [0.0, 1.0, 1.0, 0.0], [1.0, 0.0, 0.0, 0.0], None).unwrap()
}

#[test]
fn absent_selection_keeps_original_bytes_and_allocation() {
    let mut desc = description();
    let before = desc.courses.clone(); let pointer = desc.courses.as_ptr();
    let board = desc.board_bytes();
    let report = preparation::prepare(&mut desc, None, None).unwrap();
    assert_eq!(desc.courses, before); assert_eq!(desc.courses.as_ptr(), pointer);
    assert_eq!(desc.board_bytes(), board);
    assert!(report.octaves.is_empty()); assert!(!report.imported_scale);
    assert_eq!(report.json(), "{\"importedScale\":false,\"courses\":88,\"tuning\":null}");
}

#[test]
fn live_preparation_is_the_shared_solver_with_only_physical_tension_changed() {
    for weights in [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 1.0, 0.0]] {
        let mut desc = description(); let original = desc.courses.clone();
        let board = desc.board_bytes();
        let policy = Tuning::new(442.0, weights, [1.0, 0.0, 0.0, 0.0], None).unwrap();
        let mut expected = original.clone();
        let reference = stretch_tuning::partials::tune(&mut expected, 442.0, policy.profile()).unwrap();
        let report = preparation::prepare(&mut desc, None, Some(policy)).unwrap();
        assert_eq!(desc.courses, expected); assert_eq!(desc.board_bytes(), board);
        assert_eq!(report.octaves.len(), 76);
        for (a, b) in report.octaves.iter().zip(reference) {
            assert_eq!(a.lower_midi, b.lower_midi);
            assert_eq!(a.signed_beats_hz, b.signed_beats_hz);
            assert_eq!(a.weighted_rms_hz.to_bits(), b.weighted_rms_hz.to_bits());
        }
        for (old, new) in original.iter().zip(&desc.courses) {
            assert_eq!(*old, Course { tension_n: old.tension_n, ..*new });
        }
        let a4 = desc.courses.iter().find(|c| c.midi == 69).unwrap();
        assert!((a4.partial_hz(1, a4.tension_n) - 442.0).abs() < 1e-9);
        preparation::prepare(&mut desc, None, Some(policy)).unwrap();
        assert_eq!(desc.courses, expected); // no accumulated retuning
    }
}

#[test]
fn imported_csv_preserves_author_physics_unless_tuning_is_explicit() {
    let mut desc = description(); let board = desc.board_bytes();
    let mut source = desc.courses.clone();
    for c in &mut source {
        c.tension_n *= 1.02; c.flexural_rigidity_nm2 *= 1.1;
        c.hammer_mass_kg *= 1.05; c.felt_thickness_m *= 0.95;
    }
    let csv = geometry::write_scale(&source);
    let report = preparation::prepare(&mut desc, Some(&csv), None).unwrap();
    assert!(report.imported_scale); assert!(report.tuning.is_none());
    assert_eq!(desc.courses, source); assert_eq!(desc.board_bytes(), board);
    let mut expected = source;
    stretch_tuning::partials::tune(&mut expected, 442.0, tuning().profile()).unwrap();
    preparation::prepare(&mut desc, Some(&csv), Some(tuning())).unwrap();
    assert_eq!(desc.courses, expected); assert_eq!(desc.board_bytes(), board);
}

#[test]
fn failed_import_or_late_tuning_retains_the_whole_original_description() {
    let mut desc = description(); let before = desc.courses.clone(); let board = desc.board_bytes();
    for csv in [String::new(), "not a scale".into(), "x".repeat(preparation::MAX_SCALE_BYTES + 1),
        geometry::write_scale(&[before[0], before[0]])] {
        assert!(preparation::prepare(&mut desc, Some(&csv), Some(tuning())).is_err());
        assert_eq!(desc.courses, before); assert_eq!(desc.board_bytes(), board);
    }
    let source: Vec<_> = before.iter().copied().filter(|c| c.midi != 81).collect();
    let csv = geometry::write_scale(&source);
    assert!(preparation::prepare(&mut desc, Some(&csv), Some(tuning())).is_err());
    let band = Tuning::new(442.0, [1.0,0.0,0.0,0.0], [1.0,0.0,0.0,0.0], Some(1.0)).unwrap();
    assert!(preparation::prepare(&mut desc, Some(&geometry::write_scale(&before)), Some(band)).is_err());
    assert_eq!(desc.courses, before); assert_eq!(desc.board_bytes(), board);
}

#[test]
fn typed_and_wasm_controls_refuse_nonphysical_or_ambiguous_policies() {
    assert!(Tuning::from_values(&[]).unwrap().is_none());
    let packet = [442.0, 0.0,1.0,1.0,0.0, 1.0,0.0,0.0,0.0, 0.0];
    let policy = Tuning::from_values(&packet).unwrap().unwrap();
    assert_eq!(policy.concert_a_hz(), 442.0);
    assert_eq!(policy.profile().bass.values(), [0.0,1.0,1.0,0.0]);
    assert!(policy.profile().max_partial_hz.is_none());
    for (index, value) in [(0, f64::NAN), (0, 451.0), (1, -1.0), (5, f64::INFINITY), (9, -1.0)] {
        let mut bad = packet; bad[index] = value;
        assert!(Tuning::from_values(&bad).is_err());
    }
    for length in [1, 9, 11] { assert!(Tuning::from_values(&vec![1.0; length]).is_err()); }
    assert!(Tuning::new(440.0, [0.0; 4], [1.0; 4], None).is_err());
}

#[test]
fn prepared_tensions_reach_actual_live_partials_and_high_rate_pressure() {
    let mut desc = description();
    desc.courses.retain(|c| [57,69,81].contains(&c.midi));
    for c in &mut desc.courses { c.unison = 1; c.detune_cents = 0.0; c.duplex_length_m = 0.0; }
    let mut expected = description();
    expected.courses = desc.courses.clone();
    stretch_tuning::partials::tune(&mut expected.courses, 442.0, tuning().profile()).unwrap();
    preparation::prepare(&mut desc, None, Some(tuning())).unwrap();
    let mut cfg = Config::live(48_000);
    cfg.max_partials = 12; cfg.course_mode_budget = 12; cfg.max_dynamic_modes = 36;
    cfg.max_courses = 3; cfg.duplex = false; cfg.partial_band_hz = 6_000.0;
    let native = LivePiano::new(&desc, cfg.clone()).unwrap();
    for c in &desc.courses {
        assert!((native.partials_hz(c.midi)[0] / c.partial_hz(1, c.tension_n) - 1.0).abs() < 1e-12);
    }
    for ratio in [1,2,4] {
        let mut actual = PianoAudio::new(&desc, cfg.clone(), ratio).unwrap();
        let mut reference = PianoAudio::new(&expected, cfg.clone(), ratio).unwrap();
        actual.note_on(81, 1.0).unwrap(); reference.note_on(81, 1.0).unwrap();
        let mut a = vec![0.0_f32; 768]; let mut b = a.clone();
        actual.render(&mut a).unwrap();
        for block in b.chunks_mut(37) { reference.render(block).unwrap(); }
        assert_eq!(a, b); assert!(a.iter().any(|p| *p != 0.0));
        assert_eq!(actual.piano().counters.contact_failures, 0);
    }
}
