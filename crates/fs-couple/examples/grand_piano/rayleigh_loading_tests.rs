use super::*;
use super::super::{geometry::Course, linear::BoardMode, steinway_scale};
use std::f64::consts::PI;

fn fixture() -> (Bank, Vec<SurfaceSample>) {
    let course = steinway_scale::courses().unwrap().into_iter().find(|c| c.midi == 69).unwrap();
    let course = Course { unison: 1, duplex_length_m: 0.0, ..course };
    let mut bridge = [0.0; 88]; bridge[69 - 21] = 2.0;
    let radius = 0.15;
    let board = [BoardMode { frequency_hz: 250.0, damping_ratio: 0.0,
        bridge, volume: PI * radius * radius }];
    let bank = Bank::new(&[course], &board, 192_000, 8_000.0, 20, false).unwrap();
    let rings = 4;
    let mut surface = Vec::new();
    for m in 0..rings {
        let r0 = radius * m as f64 / rings as f64;
        let r1 = radius * (m + 1) as f64 / rings as f64;
        let rc = f64::midpoint(r0, r1);
        let sectors = 6 * (m + 1);
        for sector in 0..sectors {
            let theta = TAU * (sector as f64 + 0.5) / sectors as f64;
            surface.push(SurfaceSample { position_m: [rc * det::cos(theta), rc * det::sin(theta), 0.0],
                area_m2: PI * (r1 * r1 - r0 * r0) / sectors as f64,
                mode_shape: vec![1.0] });
        }
    }
    (bank, surface)
}

fn medium() -> Medium { Medium { density: 1.2, sound_speed: 343.0 } }

#[test]
fn unsupported_complete_port_count_refuses_before_the_rayleigh_sweep() {
    let count = MAX_PORTS + 1;
    assert!(count <= MAX_MODAL_RAYLEIGH_FIELDS,
        "fixture must fit the kernel while exceeding the passive realization");
    let course = steinway_scale::courses().unwrap().into_iter().find(|c| c.midi == 69).unwrap();
    let course = Course { unison: 1, duplex_length_m: 0.0, ..course };
    let board: Vec<_> = (0..count).map(|i| BoardMode {
        frequency_hz: 250.0 + i as f64, damping_ratio: 0.0,
        bridge: [0.0; 88], volume: 0.0,
    }).collect();
    let bank = Bank::new(&[course], &board, 192_000, 8_000.0, 20, false).unwrap();
    let surface = [SurfaceSample { position_m: [0.0, 0.0, 0.0],
        area_m2: 0.001, mode_shape: vec![1.0; count] }];
    let error = prepare(&surface, &bank, &Options::parse("1,40,500,65").unwrap(),
        medium(), 48_000, None).err().expect("unsupported full image must refuse");
    assert_eq!(error, format!("Rayleigh loading retains 1..{MAX_PORTS} complete board modes and 1..{} area samples; reduce the supplied physical image explicitly", fs_bem::helmholtz::MAX_DENSE_PANELS));
}

#[test]
fn complete_loaded_projection_reaches_the_actual_passive_fit() {
    let (bank, surface) = fixture();
    let options = Options::parse("1,40,500,65").unwrap();
    let result = prepare(&surface, &bank, &options, medium(), 48_000, None).unwrap();
    let (_, reference) = bank.board_pressure_basis().unwrap();
    let loaded_shape = bank.project_board_shape(&[1.0]).unwrap()[0];
    eprintln!("piston loaded shape={loaded_shape:.17e}, reference_hz={:.17e}", det::sqrt(reference[0])/TAU);
    assert!((loaded_shape - 1.0).abs() > 1e-5,
        "fixture must exercise string mass loading of the board shape");
    assert!(result.report.reference_max_hz > 250.0);
    assert_eq!(result.report.surface_samples, surface.len());
    assert_eq!(result.report.board_modes, 1);
    assert_eq!(result.report.reference_min_hz.to_bits(), (det::sqrt(reference[0]) / TAU).to_bits());
    assert!(result.report.peak_error < 0.05 && result.report.rms_error < 0.02);
    assert!(result.report.resistance_peak_error < 0.10 && result.report.resistance_rms_error < 0.05);
    assert!(result.model.poles.iter().any(|pole| pole.coupling[0] != 0.0));
    // Use the unchanged independent piston owner, after the bank's actual
    // loaded map. Compare a held-out interior frequency as well as band edges.
    let shape = bank.project_board_shape(&[1.0]).unwrap()[0];
    for hz in [40.0, 193.0, 500.0] {
        let air = medium();
        let original = fs_bem::helmholtz::baffled_piston_impedance(0.15, TAU * hz / air.sound_speed, air, 4)
            .unwrap().scale(PI * 0.15 * 0.15 * shape * shape);
        let fitted = result.model.impedance(TAU * hz).unwrap()[0];
        assert!((fitted - original).abs() / original.abs() < 0.04,
            "loaded piston fit at {hz}Hz: fitted={fitted:?}, original={original:?}");
    }
    eprintln!("{}", result.report);
}

#[test]
fn two_half_spaces_add_reaction_without_changing_the_source_coordinates() {
    let (bank, surface) = fixture();
    let one = Options::parse("1,40,500,65").unwrap();
    let two = Options { sides: 2, ..one };
    let a = prepare(&surface, &bank, &one, medium(), 48_000, None).unwrap();
    let b = prepare(&surface, &bank, &two, medium(), 48_000, None).unwrap();
    for hz in [40.0, 193.0, 500.0] {
        let first = a.model.impedance(TAU * hz).unwrap()[0];
        let second = b.model.impedance(TAU * hz).unwrap()[0];
        assert!((second - first.scale(2.0)).abs() / second.abs() < 2e-12);
    }
    assert_eq!(a.report.reference_min_hz.to_bits(), b.report.reference_min_hz.to_bits());
    assert_eq!(a.report.surface_samples, b.report.surface_samples);
}

#[test]
fn loading_requires_explicit_sides_resolved_clocks_and_complete_reference_band() {
    for text in ["", "40,500,65", "0,40,500,65", "3,40,500,65", "1,0,500,65",
        "1,500,40,65", "1,NaN,500,65", "1,40,inf,65", "1,40,500,64", "1,40,500,259"] {
        assert!(Options::parse(text).is_err(), "accepted {text}");
    }
    let options = Options::parse("1,40,500,65").unwrap();
    assert!(options.validate_clock(48_000, 192_000).is_ok());
    assert!(options.validate_clock(48_000, 96_001).is_err());
    assert!(Options::parse("1,40,22000,65").unwrap().validate_clock(48_000, 192_000).is_err());
    assert!(Options::parse("1,40,11000,65").unwrap().validate_clock(48_000, 192_000).is_err());
    let (bank, surface) = fixture();
    let reference_hz = det::sqrt(bank.board_pressure_basis().unwrap().1[0]) / TAU;
    for invalid in [Options { min_hz: reference_hz * 1.01, ..options },
        Options { max_hz: reference_hz * 0.99, ..options }] {
        assert!(prepare(&surface, &bank, &invalid, medium(), 48_000, None).is_err());
    }
}

#[test]
fn nonplanar_or_spatially_unresolved_sources_refuse_before_fitting() {
    let (bank, mut surface) = fixture();
    let options = Options::parse("1,40,500,65").unwrap();
    surface[1].position_m[2] = 0.001;
    assert!(prepare(&surface, &bank, &options, medium(), 48_000, None).is_err());
    // A large XY translation cannot hide a physically nonplanar source.
    for cell in &mut surface { cell.position_m[0] += 1e12; }
    assert!(prepare(&surface, &bank, &options, medium(), 48_000, None).is_err());
    let (_, mut coarse) = fixture();
    coarse[0].area_m2 = 2.0;
    assert!(prepare(&coarse, &bank, &options, medium(), 48_000, None).is_err());
}

#[test]
fn physical_radiation_changes_string_board_decay_and_closes_the_combined_energy() {
    use super::super::engine::radiation::Prepared;
    // Independent square fixture: two physical strings and two signed board
    // fields. These are explicit test parameters, not measured piano data.
    let courses: Vec<_> = steinway_scale::courses().unwrap().into_iter()
        .filter(|c| [60, 69].contains(&c.midi))
        .map(|c| Course { unison: 1, detune_cents: 0.0, duplex_length_m: 0.0, ..c }).collect();
    let mut first = [0.0; 88]; let mut second = [0.0; 88];
    first[60 - 21] = 2.0; second[60 - 21] = -1.0;
    first[69 - 21] = -1.0; second[69 - 21] = 2.0;
    let modes = [
        BoardMode { frequency_hz: 240.0, damping_ratio: 0.0, bridge: first,
            volume: 8.0 * 0.48 * 0.48 / (PI * PI) },
        BoardMode { frequency_hz: 390.0, damping_ratio: 0.0, bridge: second, volume: 0.0 },
    ];
    let mut bank = Bank::new(&courses, &modes, 192_000, 10_000.0, 24, false).unwrap();
    assert!(bank.board_pressure_basis().unwrap().0[1].abs() > 0.01);
    let cells = 12;
    let surface: Vec<_> = (0..cells).flat_map(|i| (0..cells).map(move |j| {
        let x = (i as f64 + 0.5) / cells as f64;
        let y = (j as f64 + 0.5) / cells as f64;
        SurfaceSample { position_m: [0.48 * x, 0.48 * y, 0.0],
            area_m2: (0.48 / cells as f64).powi(2),
            mode_shape: vec![2.0 * (PI * x).sin() * (PI * y).sin(),
                2.0 * (2.0 * PI * x).sin() * (PI * y).sin()] }
    })).collect();
    let load = prepare(&surface, &bank, &Options::parse("1,40,500,65").unwrap(), medium(), 48_000, None).unwrap();
    let mut air = Prepared::new(&load.model, bank.rate, bank.board_count).unwrap();
    let offset = bank.modes.len();
    bank.v[offset] = 0.03; bank.v[offset + 1] = -0.02;
    let initial = bank.energy();
    let mut previous = initial;
    let zero = vec![0.0; bank.contact_strings.len()];
    let mut radiated = 0.0;
    let mut modal_roundoff = 0.0;
    let mut maximum_air = 0.0_f64;
    let mut maximum_balance = 0.0_f64;
    for _ in 0..bank.rate / 2 {
        radiated += air.before(&mut bank.v[offset..]).unwrap();
        bank.predict(); bank.finish(&zero);
        modal_roundoff += bank.last_modal_loss_j;
        bank.commit();
        radiated += air.after(&mut bank.v[offset..]).unwrap();
        maximum_air = maximum_air.max(air.energy());
        let total = bank.energy() + air.energy();
        assert!(total.is_finite() && total >= 0.0);
        assert!(total <= previous + 1e-10 * initial, "unforced air/board energy increased");
        maximum_balance = maximum_balance.max((total + radiated + modal_roundoff - initial).abs() / initial);
        previous = total;
    }
    assert!(maximum_balance < 1e-8, "combined work balance {maximum_balance:e}");
    assert!(modal_roundoff.abs() < 1e-8 * initial, "numerical drift cannot count as radiation loss");
    assert!(radiated > 0.01 * initial, "missing radiative decay");
    assert!(maximum_air > 1e-4 * initial, "missing reactive acoustic storage");
    assert!(bank.energy() < 0.95 * initial, "receiver-only output cannot replace physical reaction");
    eprintln!("Rayleigh reaction: initial={initial:e} J, radiated={radiated:e} J, mechanics={:e} J, air={:e} J, max relative ledger={maximum_balance:e}", bank.energy(), air.energy());
}

#[test]
fn projected_rayleigh_loading_builds_port_basis() {
    let (bank, surface) = fixture();
    let options = Options::parse("1,40,500,65").unwrap();
    let result = prepare(&surface, &bank, &options, medium(), 48_000, Some(1)).unwrap();
    assert!(result.basis.is_some());
    assert_eq!(result.report.acoustic_ports, 1);
    assert!(result.report.surface_rms_error.is_some());
}

#[test]
fn explicit_projection_reaches_the_existing_runtime_with_all_33_board_coordinates() {
    let count = MAX_PORTS + 1;
    let course = steinway_scale::courses().unwrap().into_iter().find(|c| c.midi == 69).unwrap();
    let course = Course { unison: 1, duplex_length_m: 0.0, ..course };
    let board: Vec<_> = (0..count).map(|i| BoardMode {
        frequency_hz: 250.0 + i as f64, damping_ratio: 0.0,
        bridge: [0.0; 88], volume: 0.0,
    }).collect();
    let bank = Bank::new(&[course], &board, 192_000, 8_000.0, 20, false).unwrap();
    let surface = [SurfaceSample { position_m: [0.0, 0.0, 0.0],
        area_m2: 0.001, mode_shape: vec![1.0; count] }];
    let options = Options::parse("1,40,500,65").unwrap();
    for invalid in [0, MAX_PORTS + 1] {
        assert!(prepare(&surface, &bank, &options, medium(), 48_000, Some(invalid)).is_err());
    }
    let load = prepare(&surface, &bank, &options, medium(), 48_000, Some(1)).unwrap();
    let basis = load.basis.as_ref().unwrap();
    assert_eq!(bank.board_count, count);
    assert_eq!(load.report.board_modes, count);
    assert_eq!(basis.board_ports(), count);
    assert_eq!(basis.ports(), 1);
    assert_eq!(load.model.ports, 1);
    assert_eq!(load.report.acoustic_ports, 1);
    assert!(load.report.surface_rms_error.unwrap() < 1e-12);
    assert!(load.report.peak_error < 0.05 && load.report.rms_error < 0.02);
    assert!(load.report.resistance_peak_error < 0.10 && load.report.resistance_rms_error < 0.05);
    let mut air = radiation::Prepared::new_projected(&load.model, bank.rate, basis).unwrap();
    let mut velocity = vec![0.01; count];
    let initial = 0.5 * velocity.iter().map(|v| v * v).sum::<f64>();
    let mut loss = 0.0;
    for _ in 0..480 {
        loss += air.before(&mut velocity).unwrap() + air.after(&mut velocity).unwrap();
    }
    let final_energy = 0.5 * velocity.iter().map(|v| v * v).sum::<f64>() + air.energy();
    assert!((final_energy + loss - initial).abs() < 2e-13);
    assert!(air.energy() > 0.0 && loss > 0.0);
    assert_eq!(velocity.len(), count);
}
