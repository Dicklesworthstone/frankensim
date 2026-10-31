//! Geometric string extension on the live kernel (`Config::string_stretching`,
//! `src/stretch.rs`): the reference's own strain, storage and file, with one
//! kick per step. The wider comparison is `examples/validate.rs --string-waves`.
use fs_live_piano::refmods::{engine, geometry, linear::string_stretching::{Specification, HEADER}, steinway_scale};
use fs_live_piano::stretch::plain_steel_waves_spec;
use fs_live_piano::{Config, Description, LivePiano};

const RATE: u32 = 48_000;

fn scale() -> Vec<geometry::Course> {
    geometry::read_scale(fs_live_piano::kernel::SCALE_CSV).unwrap()
}

fn single(midi: u8) -> (geometry::Course, Description) {
    let course = scale()[usize::from(midi - 21)];
    (course, Description::with_courses(vec![course], 800.0).unwrap())
}

/// `validate.rs`'s isolated-key configuration.
fn config(spec: Option<String>) -> Config {
    let mut cfg = Config::live(RATE);
    cfg.max_partials = 256;
    cfg.course_mode_budget = usize::MAX / 4;
    cfg.static_attachment = false;
    cfg.damper_points = 1;
    cfg.string_stretching = spec;
    cfg
}

fn bridge_velocities(piano: &mut LivePiano, midi: u8, steps: usize) -> Vec<f64> {
    (0..steps).map(|_| {
        piano.step();
        piano.bridge_velocity(midi).unwrap()
    }).collect()
}

fn struck(desc: &Description, spec: Option<String>, midi: u8, velocity: f64) -> LivePiano {
    let mut piano = LivePiano::new(desc, config(spec)).unwrap();
    piano.note_on(midi, velocity).unwrap();
    piano.apply_pending_notes();
    piano
}

#[test]
fn plain_steel_rows_follow_the_scale_and_wound_strings_stay_linear() {
    let scale = scale();
    let spec = Specification::read(&plain_steel_waves_spec(&scale, 10_000.0, 0.01, 0.05), &scale).unwrap();
    // The graded A3 file (axial-a3.fspx) has EA 166505 N from the same rule (d = 1.030 mm).
    let a3 = spec.courses[&57].unwrap().axial_rigidity_n;
    assert!((a3 / 166_505.0 - 1.0).abs() < 1e-3, "A3 EA {a3}");
    // The wound bass (core about a third of the mass) stays linear; MIDI 41-108 are plain steel.
    let waves: Vec<u8> = spec.waves.keys().copied().collect();
    assert_eq!(waves, (41..=108).collect::<Vec<u8>>());
    assert!(spec.courses[&40].is_none());
}

#[test]
fn an_all_linear_file_keeps_the_linear_strings_bit_for_bit() {
    let (_, desc) = single(57);
    let mut plain = struck(&desc, None, 57, 4.5);
    let mut linear = struck(&desc, Some(format!("{HEADER}\nlinear,57\n")), 57, 4.5);
    for _ in 0..4800 {
        assert_eq!(plain.step().to_bits(), linear.step().to_bits());
    }
}

#[test]
fn a3_waves_follow_the_reference_engine_and_close_the_ledger() {
    let (course, desc) = single(57);
    let text = plain_steel_waves_spec(&[course], 10_000.0, 0.01, 0.05);
    let material = steinway_scale::hammer_material(&course).unwrap();
    let mut reference = engine::Instrument::new_with_string_damping(vec![course], &desc.board, RATE, 1, 256, true,
        vec![material], Some(engine::ShankGeometry::published()), None, None, false).unwrap();
    if let Some(c) = &desc.damping {
        reference.configure_bare_board_damping(c).unwrap();
    }
    // The reference's waves coupling is a per-mode row; the live engine's is one scalar. To compare the same
    // equations the reference is given the engine's own labelled scalar times its vertical bridge row.
    let mut material_spec = Specification::read(&text, &[course]).unwrap();
    let rows = material_spec.waves.keys().map(|&key| {
        (key, desc.board.iter().map(|m| fs_live_piano::stretch::LIVE_BRIDGE_AXIAL_RATIO * m.bridge[usize::from(key - 21)]).collect())
    }).collect();
    material_spec.set_bridge_axial_rows(rows).unwrap();
    reference.configure_string_stretching(&material_spec).unwrap();
    reference.note_on(57, 4.5).unwrap();
    let mut live = struck(&desc, Some(text), 57, 4.5);
    let mut linear = struck(&desc, None, 57, 4.5);
    let row = reference.bank.strings[0].bridge.clone();
    let n = reference.bank.modes.len();
    let (mut peak, mut to_reference, mut to_linear) = (0.0_f64, 0.0_f64, 0.0_f64);
    for _ in 0..4800 {
        reference.step().unwrap();
        let want: f64 = row.iter().zip(&reference.bank.v[n..]).map(|(g, v)| g * v).sum();
        live.step();
        linear.step();
        let got = live.bridge_velocity(57).unwrap();
        peak = peak.max(want.abs());
        to_reference = to_reference.max((got - want).abs());
        to_linear = to_linear.max((linear.bridge_velocity(57).unwrap() - want).abs());
    }
    // Two methods for the same equations (measured 9.7e-4 over 0.2 s, against the
    // reference's own 1.2e-3 between 48 and 192 kHz), versus the waves' own effect.
    assert!(to_reference < 2e-3 * peak, "live vs reference {:.3e}", to_reference / peak);
    assert!(to_linear > 5e-2 * peak, "waves changed the motion by only {:.3e}", to_linear / peak);
    let ledger = live.close_ledger();
    assert!(ledger.residual_j.abs() < 1e-12 * ledger.input_j, "residual {:e} J", ledger.residual_j);
    // The method's own work stays far below the strike's (measured -9e-6 J of 0.12 J).
    assert!(ledger.coupling_j.abs() < 1e-3 * ledger.input_j, "coupling {:e} J", ledger.coupling_j);
    assert_eq!(live.counters.stretch_out_of_domain, 0);
}

#[test]
fn a_clone_continues_the_stretching_strings_bit_for_bit() {
    let desc = Description::embedded().unwrap();
    let mut cfg = Config::live(RATE);
    cfg.string_stretching = Some(plain_steel_waves_spec(&desc.courses, 10_000.0, 0.01, 0.05));
    let mut piano = LivePiano::new(&desc, cfg).unwrap();
    piano.note_on(45, 5.0).unwrap();
    piano.note_on(57, 4.0).unwrap();
    for _ in 0..2000 {
        piano.step();
    }
    let mut copy = piano.clone();
    for _ in 0..4000 {
        assert_eq!(piano.step().to_bits(), copy.step().to_bits());
    }
    assert_eq!(piano.close_ledger().residual_j.to_bits(), copy.close_ledger().residual_j.to_bits());
}
