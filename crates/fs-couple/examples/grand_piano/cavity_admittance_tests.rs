//! The supplied sealed enclosure through the force-driven harmonic frontend.
use super::*;
use fs_math::c64::C64;
use super::cavity_exterior_tests::{card, inputs};

fn fresh_path() -> std::path::PathBuf {
    static NEXT: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    let stamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    std::env::temp_dir().join(format!("frankensim-cavity-admittance-{}-{stamp}-{}.fspc",
        std::process::id(), NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed)))
}
fn write_card(text: &str) -> String {
    let path = fresh_path();
    std::fs::OpenOptions::new().create_new(true).write(true).open(&path).unwrap()
        .write_all(text.as_bytes()).unwrap();
    path.to_str().unwrap().to_owned()
}
fn table(csv: &str) -> (Vec<&str>, Vec<Vec<&str>>) {
    let mut lines = csv.lines().filter(|line| !line.starts_with('#') && !line.is_empty());
    let header = lines.next().unwrap().split(',').collect();
    let rows = lines.map(|line| line.split(',').collect()).collect();
    (header, rows)
}
fn value(header: &[&str], row: &[&str], name: &str) -> f64 {
    let index = header.iter().position(|&field| field == name).unwrap();
    row[index].parse().unwrap()
}
fn complex(header: &[&str], row: &[&str], one_way: bool) -> C64 {
    let (real, imag) = if one_way { ("one_way_real", "one_way_imag") } else { ("real", "imag") };
    C64::new(value(header, row, real), value(header, row, imag))
}
fn close(actual: C64, expected: C64) {
    assert!((actual - expected).abs() <= 1e-10 * expected.abs().max(f64::MIN_POSITIVE),
        "harmonic frontend changed owner result: {actual:?} versus {expected:?}");
}

#[test]
fn supplied_cavity_changes_harmonic_mobility_pressure_and_complete_air_power() {
    let (board, courses, obj, spec) = inputs();
    let path = write_card(card());
    let (options, damping) = admittance_options(&[
        "--cavity".into(), path, "--modes".into(), "12".into()]).unwrap();
    assert!(damping);
    let plain = playback::Options { modes: 12, ..playback::Options::default() };
    let baseline = admittance_controlled(&board, &courses, &obj, &spec, 69, &plain, true).unwrap();
    let selected = admittance_controlled(&board, &courses, &obj, &spec, 69, &options, true).unwrap();
    assert!(selected.contains("authored sealed acoustic integration fixture"));
    let (old_header, old_rows) = table(&baseline);
    let (header, rows) = table(&selected);
    assert_eq!(old_header.len(), 15); assert_eq!(header.len(), 16);
    assert_eq!(header.iter().copied().filter(|field| *field != "cavity_w").collect::<Vec<_>>(), old_header);
    assert_eq!(rows.len(), spec.frequencies * (courses.len() + spec.receivers.len()));
    assert_eq!(rows.len(), old_rows.len());
    let mut changed_bridge = [false; 2]; let mut changed_pressure = [false; 2];
    let mut cavity_dissipates = false; let mut radiates = false;
    for (row, old) in rows.iter().zip(&old_rows) {
        assert_eq!(row.len(), 16); assert_eq!(old.len(), 15); assert_eq!(&row[..3], &old[..3]);
        assert!(row.iter().enumerate().filter(|(i, _)| *i != 1)
            .all(|(_, text)| text.parse::<f64>().unwrap().is_finite()));
        for (slot, one_way) in [false, true].into_iter().enumerate() {
            let actual = complex(&header, row, one_way); let previous = complex(&old_header, old, one_way);
            let changed = (actual - previous).abs() > 1e-6 * previous.abs();
            if row[1] == "bridge" { changed_bridge[slot] |= changed; }
            if row[1] == "receiver" { changed_pressure[slot] |= changed; }
        }
        let input = value(&header, row, "input_w");
        let wood = value(&header, row, "wood_w"); let strings = value(&header, row, "string_w");
        let radiation = value(&header, row, "radiation_w"); let air = value(&header, row, "cavity_w");
        let scale = input.abs() + wood.abs() + strings.abs() + radiation.abs() + air.abs();
        assert!((input - wood - strings - radiation - air).abs() < 1e-12 + 1e-7 * scale);
        assert!(value(&header, row, "power_defect_w").abs() < 1e-12 + 1e-7 * scale);
        assert!(value(&header, row, "backward_error") < 1e-9);
        assert!(air >= 0.); cavity_dissipates |= air > 0.; radiates |= radiation > 0.;
    }
    assert_eq!(changed_bridge, [true; 2], "cavity reaction must remain in both mobility comparisons");
    assert_eq!(changed_pressure, [true; 2], "both receiver comparisons must observe cavity-loaded motion");
    assert!(cavity_dissipates && radiates);
}

#[test]
fn conservative_structure_keeps_cavity_loss_and_matches_manual_geometric_owner() {
    let (source, courses, obj, spec) = inputs();
    let path = write_card(card());
    let (options, damping) = admittance_options(&[
        "--cavity".into(), path, "--modes".into(), "12".into(), "--lossless-structure".into()]).unwrap();
    assert!(!damping);
    let csv = admittance_controlled(&source, &courses, &obj, &spec, 69, &options, damping).unwrap();
    let (header, rows) = table(&csv);
    let mut dissipates = false;
    for row in &rows {
        assert_eq!(value(&header, row, "wood_w"), 0.);
        assert_eq!(value(&header, row, "string_w"), 0.);
        let input = value(&header, row, "input_w");
        let radiation = value(&header, row, "radiation_w"); let air = value(&header, row, "cavity_w");
        assert!((input - radiation - air).abs() < 1e-12 + 1e-7 * input.abs());
        dissipates |= air > 0.;
    }
    assert!(dissipates, "lossless structure must not erase the separately supplied acoustic drag");
    // Construct the same geometric projection and continuous owner independently
    // of CLI admission, then use the actual BEM load at one in-band frequency.
    let keys: Vec<_> = courses.iter().map(|course| course.midi).collect();
    let board = prepare_board_motion(&source, &keys, spec.board_band_hz, &options).unwrap();
    let mut model = bridge_response::BridgeResponse::new(&courses, &board.modes,
        RATE * options.substeps as u32, 0.45 * f64::from(RATE), options.modes, damping).unwrap();
    if let Some(c) = &board.physical_damping { model.configure_bare_board_damping(c).unwrap(); }
    let projection = cavity::Specification::read(card()).unwrap().project(&board).unwrap();
    model.configure_cavity(&projection.loaded(model.bank()).unwrap()).unwrap();
    let boundary = Boundary::from_obj(&obj, &spec, board.motion.as_ref().unwrap()).unwrap()
        .loaded(model.bank()).unwrap();
    let omega = spec.omega()[spec.frequencies / 2]; let hz = omega / std::f64::consts::TAU;
    let load = exterior_loading::sample(&boundary, &spec, omega).unwrap();
    for one_way in [false, true] {
        let solved = model.solve(hz, 69, C64::ONE, if one_way { None } else { Some(&load.impedance) }).unwrap();
        let pressure = load.pressure(&solved, omega).unwrap();
        assert!(solved.cavity_loss_w > 0.);
        let chosen: Vec<_> = rows.iter().filter(|row| value(&header, row, "frequency_hz") == hz).collect();
        assert_eq!(chosen.len(), courses.len() + spec.receivers.len());
        for row in chosen {
            let index = row[2].parse::<usize>().unwrap();
            let expected = if row[1] == "bridge" {
                solved.bridge_velocity[courses.iter().position(|course| usize::from(course.midi) == index).unwrap()]
            } else { pressure[index] };
            close(complex(&header, row, one_way), expected);
            if !one_way {
                assert!((value(&header, row, "cavity_w") - solved.cavity_loss_w).abs()
                    < 1e-12 * solved.cavity_loss_w);
            }
        }
    }
    assert!(model.bank().q.iter().chain(&model.bank().v).all(|value| *value == 0.));
}

#[test]
fn harmonic_cavity_rejects_missing_cards_incompatible_exteriors_and_prescribed_response() {
    let (board, courses, obj, spec) = inputs();
    let missing = fresh_path(); assert!(!missing.exists());
    let options = playback::Options { cavity: Some(missing.to_str().unwrap().to_owned()),
        modes: 12, ..playback::Options::default() };
    let error = admittance_controlled(&board, &courses, &obj, &spec, 69, &options, true).unwrap_err();
    assert!(error.contains(missing.to_str().unwrap()), "{error}");
    let invalid = playback::Options { cavity: Some(write_card("not a cavity specification\n")),
        ..options.clone() };
    // Source-card failure precedes structural preparation; it must not be
    // swallowed into an uncoupled or differently sourced board calculation.
    let error = admittance_controlled("invalid structural source", &courses, &obj, &spec,
        69, &invalid, true).unwrap_err();
    assert!(error.contains(cavity::HEADER), "{error}");
    let path = write_card(card());
    let valid = playback::Options { cavity: Some(path.clone()), ..options };
    let error = admittance_controlled_body(&board, &courses, None, &spec, 69, &valid, true, false).unwrap_err();
    assert!(error.contains("outer enclosure OBJ"), "{error}");
    let (bare_board, _, bare_skin, bare_spec) = tests::small_source_inputs();
    let error = admittance_controlled(&bare_board, &courses, &bare_skin, &bare_spec, 69, &valid, true).unwrap_err();
    assert!(error.contains("moving board underside"), "{error}");
    let output = fresh_path();
    let error = run(&["response".into(), "missing.fsb".into(), "missing.csv".into(),
        "missing.obj".into(), "missing.fspe".into(), output.to_str().unwrap().into(),
        "--cavity".into(), path]).unwrap_err();
    assert!(error.contains("cavity"), "{error}");
    assert!(!output.exists());
}
