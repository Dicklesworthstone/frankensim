//! Geometric string extension on the live kernel (`Config::string_stretching`,
//! `src/stretch.rs`): the reference's own strain, storage and file, with one
//! kick per step. The wider comparison is `examples/validate.rs --string-waves`.
use fs_live_piano::refmods::{engine, geometry, linear::string_stretching::{Specification, HEADER}, steinway_scale};
use fs_live_piano::stretch::plain_steel_waves_spec;
use fs_live_piano::{Config, Description, LivePiano};
pub use fs_live_piano::refmods;

// Compile the actual private adapter here so its force pullbacks have a
// direct G0 test without exposing implementation state as a public API.
#[path = "../src/stretch.rs"]
#[allow(dead_code)]
mod stretch_adapter;

const RATE: u32 = 48_000;

/// G0: both independent termination forces are gradients of the actual
/// owner's stored energy, including when no axial oscillator fits the band.
#[test]
fn independent_transverse_and_axial_reactions_match_stored_energy_gradients() {
    use refmods::linear::string_stretching::{Material, WaveMaterial};
    for band in [100.0, 10_000.0] {
        let mut model = stretch_adapter::Stretch::new(
            Material { axial_rigidity_n: 150_000.0, maximum_slope: 0.2 },
            Some(WaveMaterial { bandwidth_hz: band, mode_cap: 8, damping_ratio: 0.01, maximum_strain: 0.01 }),
            1.0, 0.006, 500.0, 4, RATE).unwrap();
        assert_eq!(model.axial_modes() > 0, band > 100.0);
        let (q, beta) = ([2e-5, -1e-5, 4e-6, 3e-6], [0.02, -0.01, 0.006, -0.004]);
        let (b, axial, h) = (3e-6, -2e-6, 1e-9);
        let kick = model.step(&q, &beta, b, axial, 1.0 / f64::from(RATE), &mut [0.0; 4]);
        let transverse = (model.stored(&q, &beta, b + h, axial) - model.stored(&q, &beta, b - h, axial)) / (2.0 * h);
        let longitudinal = (model.stored(&q, &beta, b, axial + h) - model.stored(&q, &beta, b, axial - h)) / (2.0 * h);
        for (got, want) in [(kick.reaction, transverse), (kick.axial_reaction, longitudinal)] {
            assert!((got - want).abs() < 1e-7 * want.abs().max(1e-8), "band {band}: reaction {got:e}, energy gradient {want:e}");
        }
        // A signed two-coordinate row and a separate vertical row pull both
        // gradients back to the same board energy, including opposite signs.
        for (vertical, axial_row) in [(0.4, -0.07), (-0.2, 0.11)] {
            let gradient = (model.stored(&q, &beta, b + h * vertical, axial + h * axial_row)
                - model.stored(&q, &beta, b - h * vertical, axial - h * axial_row)) / (2.0 * h);
            let force = vertical * kick.reaction + axial_row * kick.axial_reaction;
            assert!((force - gradient).abs() < 1e-7 * gradient.abs().max(1e-8));
        }
    }
}

#[test]
fn stretch_without_waves_ignores_an_axial_port_exactly() {
    use refmods::linear::string_stretching::Material;
    let mut a = stretch_adapter::Stretch::new(
        Material { axial_rigidity_n: 150_000.0, maximum_slope: 0.2 },
        None, 1.0, 0.006, 500.0, 4, RATE).unwrap();
    let mut b = a.clone();
    let (q, beta) = ([2e-5, -1e-5, 4e-6, 3e-6], [0.02, -0.01, 0.006, -0.004]);
    let (mut va, mut vb) = ([0.0; 4], [0.0; 4]);
    let ka = a.step(&q, &beta, 3e-6, 0.0, 1.0 / f64::from(RATE), &mut va);
    let kb = b.step(&q, &beta, 3e-6, 0.25, 1.0 / f64::from(RATE), &mut vb);
    assert_eq!(va.map(f64::to_bits), vb.map(f64::to_bits));
    assert_eq!(ka.reaction.to_bits(), kb.reaction.to_bits());
    assert_eq!(ka.axial_reaction, 0.0);
    assert_eq!(kb.axial_reaction, 0.0);
    assert_eq!(a.stored(&q, &beta, 3e-6, 0.0).to_bits(), b.stored(&q, &beta, 3e-6, 0.25).to_bits());
}

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
    let axial = &desc.bridge_axial_rows[&57];
    assert_eq!(axial.len(), desc.board.len());
    assert!(axial.iter().any(|a| a.abs() > 1e-8));
    assert!(axial.iter().zip(&desc.board).any(|(a, mode)|
        (*a - 0.06 * mode.bridge[36]).abs() > 1e-8), "geometry must replace the former authored scalar");
    let material = steinway_scale::hammer_material(&course).unwrap();
    let mut reference = engine::Instrument::new_with_string_damping(vec![course], &desc.board, RATE, 1, 256, true,
        vec![material], Some(engine::ShankGeometry::published()), None, None, false).unwrap();
    if let Some(c) = &desc.damping {
        reference.configure_bare_board_damping(c).unwrap();
    }
    // The reference and live engine consume the same geometry-derived axial
    // termination motion, including the bridge's retained rotations.
    let mut material_spec = Specification::read(&text, &[course]).unwrap();
    let rows = material_spec.waves.keys().map(|&key| {
        (key, desc.bridge_axial_rows[&key].clone())
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
    // Compare the two integration methods against the waves' own effect.
    // The acceptance bounds are unchanged from the original scalar-port test.
    assert!(to_reference < 2e-3 * peak, "live vs reference {:.3e}", to_reference / peak);
    assert!(to_linear > 5e-2 * peak, "waves changed the motion by only {:.3e}", to_linear / peak);
    let ledger = live.close_ledger();
    assert!(ledger.residual_j.abs() < 1e-12 * ledger.input_j, "residual {:e} J", ledger.residual_j);
    // The method's own work must stay far below the strike's.
    assert!(ledger.coupling_j.abs() < 1e-3 * ledger.input_j, "coupling {:e} J", ledger.coupling_j);
    assert_eq!(live.counters.stretch_out_of_domain, 0);
}

#[test]
fn a_clone_continues_the_stretching_strings_bit_for_bit() {
    let courses = scale().into_iter().filter(|c| [45, 57].contains(&c.midi)).collect();
    let desc = Description::with_courses(courses, 800.0).unwrap();
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

/// G0: existing board files keep their original linear interpretation. They
/// cannot invent a longitudinal bridge port, while `stretch` needs none.
#[test]
fn legacy_board_keeps_linear_and_stretch_but_refuses_waves_without_geometry() {
    let mut desc = Description::embedded().unwrap();
    desc.courses.retain(|c| c.midi == 57);
    assert!(desc.bridge_axial_rows.is_empty());
    assert!(LivePiano::new(&desc, config(None)).is_ok());
    assert!(LivePiano::new(&desc, config(Some(format!("{HEADER}\nstretch,57,166505,0.2\n")))).is_ok());
    let waves = plain_steel_waves_spec(&desc.courses, 10_000.0, 0.01, 0.05);
    let error = LivePiano::new(&desc, config(Some(waves))).err().unwrap();
    assert!(error.contains("key 57") && error.contains("geometry-derived axial bridge row"), "{error}");
}

/// G0/G3: a serialized board carries signed modal axial motion exactly, and
/// the existing exact-prefix truncation carries the same coordinates along.
#[test]
fn axial_board_rows_round_trip_and_truncate_with_their_modal_basis() {
    let mut desc = Description::embedded().unwrap();
    let r = desc.board.len();
    let row: Vec<f64> = (0..r).map(|i| (i as f64 - 7.0) * 0.0125).collect();
    desc.bridge_axial_rows.insert(57, row.clone());
    desc.bridge_axial_rows.insert(60, vec![0.0; r]);
    let bytes = desc.board_bytes();
    assert_eq!(&bytes[..8], b"FSLPBRD2");
    let back = Description::from_board_bytes(desc.courses.clone(), &bytes).unwrap();
    assert_eq!(back.bridge_axial_rows, desc.bridge_axial_rows);
    assert_eq!(back.board_bytes(), bytes);
    let prefix = back.truncated(4).unwrap();
    assert_eq!(prefix.bridge_axial_rows[&57], row[..4]);
    assert_eq!(prefix.bridge_axial_rows[&60], vec![0.0; 4]);
    assert!(Description::from_board_bytes(desc.courses.clone(), &bytes[..bytes.len() - 1]).is_err());
    desc.bridge_axial_rows.get_mut(&57).unwrap().pop();
    let error = LivePiano::new(&desc, config(None)).err().unwrap();
    assert!(error.contains("one finite coefficient per board coordinate"), "{error}");
}

/// G3: an independent axial port changes the actual rendered mechanics. A
/// zero geometric row is accepted as zero, never replaced with a heuristic.
#[test]
fn geometric_axial_motion_reaches_sounding_strings_and_clone_rebuilds() {
    let courses = scale().into_iter().filter(|c| [57, 60].contains(&c.midi)).collect();
    let desc = Description::with_courses(courses, 800.0).unwrap();
    let text = plain_steel_waves_spec(&desc.courses, 10_000.0, 0.01, 0.05);
    let bytes = desc.board_bytes();
    let mut zero = Description::from_board_bytes(desc.courses.clone(), &bytes).unwrap();
    for row in zero.bridge_axial_rows.values_mut() { row.fill(0.0); }
    let mut linked = struck(&desc, Some(text.clone()), 57, 4.5);
    let mut unlinked = struck(&zero, Some(text), 57, 4.5);
    let mut difference = 0.0_f64;
    for _ in 0..1200 {
        linked.step(); unlinked.step();
        difference = difference.max((linked.bridge_velocity(57).unwrap() - unlinked.bridge_velocity(57).unwrap()).abs());
    }
    assert!(difference > 1e-10, "axial geometry did not reach the sounding string: {difference:e}");
    let mut copy = linked.clone();
    for piano in [&mut linked, &mut copy] {
        piano.note_on(60, 2.0).unwrap();
        piano.apply_pending_notes();
    }
    for _ in 0..2400 {
        assert_eq!(linked.step().to_bits(), copy.step().to_bits());
    }
    for piano in [&mut linked, &mut copy] {
        assert!(piano.take_error().is_none());
        let ledger = piano.close_ledger();
        assert!(ledger.residual_j.abs() < 1e-11 * ledger.input_j, "residual {:e} J", ledger.residual_j);
    }
}
