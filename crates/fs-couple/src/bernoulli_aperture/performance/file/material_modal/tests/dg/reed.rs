//! The typed, sourced lumped reed of the DG-valve file (`reed` record).
use super::{pack, render, valve, BladePropagator, CancelGate, LoadedAperture,
    MaterialModalDgValvePerformance, PlateValveInputError, PressureRenderer};

const REED: &str = "reed 0.000357 0.0056 17720 2.2e-05 3464 0.114 test:oboe-description-reed-rows";
// The lay the reed-performance file gives a massive reed: 1e7 x width, exponent 2, loss
// chi = c / (k_lay H^2) with c = 2 zeta sqrt(k m) = 0.0629 N s/m, H = 0.000357 m: about 8.82 s/m.
const LAY: &str = "contact 56000 2 8.82";

fn card(reed: &str) -> String {
    let text = super::input()
        .replace("valve 0.0004", "valve 0.000357")
        .replace("slit_length 0.001", "slit_length none")
        .replace("contact 100000000 2 0.5", LAY);
    let (head, rest) = text.split_once("loss tensor").expect("the shipped card has a loss row");
    let (_, tail) = rest.split_once("compile_limits").expect("and a compile_limits row");
    format!("{head}{reed}\ncompile_limits{tail}")
}

fn refused(text: &str) -> (usize, &'static str) {
    match MaterialModalDgValvePerformance::from_bytes_with_sources(
        text.as_bytes(), 37, &CancelGate::new_clock_free(), None,
        |_| Ok(pack(900.0, false)),
    ) {
        Err(PlateValveInputError::Input { line, what }) => (line, what),
        Err(other) => panic!("refused by another owner than the input record: {other:?}"),
        Ok(_) => panic!("a refused reed record was accepted"),
    }
}

#[test]
fn the_typed_reed_is_one_modal_mode_with_the_stated_numbers_and_no_plate() {
    let source = super::build(&card(REED), 37);
    let info = source.info();
    assert_eq!((info.blades, info.modes, info.nodes, info.triangles), (1, 1, 0, 0));
    assert!(source.material_records().is_empty() && source.solid_record().is_none());
    let record = source.reed_record().expect("the typed reed is retained");
    assert_eq!(record.source(), "test:oboe-description-reed-rows");
    assert_eq!(
        (record.rest_opening_m(), record.width_m(), record.closing_pressure_pa(),
            record.mass_kg(), record.stiffness_n_m(), record.damping_ratio()),
        (0.000357, 0.0056, 17720.0, 2.2e-5, 3464.0, 0.114),
    );
    let aperture = valve(&source);
    let spec = aperture.spec();
    let load = source.renderer().system().load();
    assert_eq!(spec.modes.len(), 1);
    let mode = spec.modes[0];
    assert_eq!((mode.mass_kg, mode.stiffness_n_m), (2.2e-5, 3464.0));
    // A = (k / P_c) H and c = 2 zeta sqrt(k m), bit for bit as the live reed (`reed_bore`) derives them.
    assert_eq!(mode.pressure_area_m2.to_bits(), (3464.0_f64 / 17720.0 * 0.000357).to_bits());
    assert_eq!(mode.damping_n_s_m.to_bits(),
        (2.0 * 0.114 * fs_math::det::sqrt(3464.0 * 2.2e-5)).to_bits());
    assert_eq!((mode.gap_projection, mode.force_projection), (1.0, [0.0, 0.0]));
    assert_eq!((spec.rest_opening_m, spec.width_m), (0.000357, 0.0056));
    assert_eq!(spec.density_kg_m3.to_bits(), load.duct().inlet_gas().density.to_bits());
    assert_eq!(spec.impedance_pa_s_m3.to_bits(), load.impedance_pa_s_m3().to_bits());
    assert_eq!(aperture.blade_propagator(), BladePropagator::ExactZoh);
    let rest = aperture.modal_states();
    assert_eq!(rest.len(), 1);
    assert_eq!((rest[0].displacement_m, rest[0].velocity_m_s), (0.0, 0.0));
}

#[test]
fn every_unstated_or_contradictory_reed_field_is_refused_for_its_own_reason() {
    let positive = "a typed reed needs positive opening, width, closing pressure, mass and stiffness; none is derived";
    for (reed, expected) in [
        ("reed 0.000357 0.0056 17720 0 3464 0.114 t", positive),
        ("reed 0.000357 0.0056 17720 2.2e-05 0 0.114 t", positive),
        ("reed 0.000357 0.0056 0 2.2e-05 3464 0.114 t", positive),
        ("reed 0.000357 0 17720 2.2e-05 3464 0.114 t", positive),
        ("reed 0.000357 0.0056 17720 2.2e-05 3464 -0.1 t", "a typed reed needs a nonnegative damping ratio"),
        ("reed 0.0004 0.0056 17720 2.2e-05 3464 0.114 t", "reed HEIGHT must equal the valve REST_GAP"),
        ("reed 0.000357 0.0056 17720 2.2e-05 3464 0.114", "missing field"),
        ("reed 0.000357 0.0056 17720 2.2e-05 3464 0.114 t u", "unexpected extra field"),
        ("reed 0.000357 0.0056 17720 2.2e-05 3464 0.114 sourcé",
            "a typed reed needs a printable source token of at most 256 bytes"),
        ("reed 0.000357 0.0056 1e-320 2.2e-05 3464 0.114 t",
            "a typed reed's pressure area k H / P_c and damping must be finite"),
    ] {
        let (line, what) = refused(&card(reed));
        assert_eq!(what, expected, "{reed}");
        assert!(line > 0, "{reed}");
    }
    let long = format!("reed 0.000357 0.0056 17720 2.2e-05 3464 0.114 {}", "a".repeat(257));
    assert_eq!(refused(&card(&long)).1,
        "a typed reed needs a printable source token of at most 256 bytes");
    let at_the_limit = format!("reed 0.000357 0.0056 17720 2.2e-05 3464 0.114 {}", "a".repeat(256));
    assert!(super::build(&card(&at_the_limit), 37).reed_record().is_some());
}

#[test]
fn a_typed_reed_selects_no_slit_and_owns_no_plate_record() {
    let slit = card(REED).replace("slit_length none", "slit_length 0.001");
    assert_eq!(refused(&slit).1, "a typed reed carries no viscous slit: declare slit_length none");
    let plate = "a typed reed has no plate: wetted surfaces, pressure shapes, closure, source forces and source pressures are plate records";
    for record in ["wetted_surfaces 0", "pressure_shapes 0", "closure", "source_forces 0", "source_pressures 0"] {
        let text = card(REED).replace("\ncompile_limits", &format!("\n{record}\ncompile_limits"));
        assert_eq!(refused(&text).1, plate, "{record}");
    }
    // `blades`, `solid` and a plate `loss` row are the other families: the reed record does not mix.
    let other = "unexpected record kind or order";
    for family in ["blades 1", "solid tet4 test:source"] {
        let text = card(REED).replace("\ncompile_limits", &format!("\n{family}\ncompile_limits"));
        assert_eq!(refused(&text).1, other, "{family}");
    }
    let loss_first = card(&format!("loss tensor eta11 eta22 eta66 eta12 eta16 eta26 1 1000000\n{REED}"));
    assert_eq!(refused(&loss_first).1, other);
}

#[test]
fn the_typed_reed_renders_finite_and_callback_invariant_on_the_dg_bore() {
    let text = card(REED);
    let (reference, rendered) = render(super::build(&text, 37));
    assert_eq!(reference.len(), 512);
    assert!(reference.iter().all(|p| p.is_finite()) && reference.iter().any(|p| p.abs() > 1e-15));
    assert_eq!(rendered.samples_rendered(), 512);
    for block in [1, 512] {
        let (other, _) = render(super::build(&text, block));
        assert_eq!(reference, other);
    }
}

#[test]
fn the_typed_reed_stays_inside_its_bounds_at_the_recipes_blowing_pressure() {
    // The posted Db4 page's operating pressure, 7400 Pa, ramped over 20 ms as the page's gesture does,
    // on the small shipped bore: the reed must move without leaving the card's physical bounds.
    let text = card(REED)
        .replace("audio 96000 96000 512 20", "audio 96000 96000 2400 20")
        .replace("events\t2\nevent\t0e0\t1e-3\tpressure-pa\t5e0\nevent\t4e-3\t1e-3\tpressure-pa\t0e0",
            "events\t1\nevent\t0e0\t2e-2\tpressure-pa\t7.4e3");
    assert!(text.contains("pressure-pa\t7.4e3"), "the gesture suffix was not replaced");
    let (pressure, rendered) = render(super::build(&text, 64));
    assert_eq!(pressure.len(), 2400);
    assert!(pressure.iter().all(|p| p.is_finite()));
    let LoadedAperture::Modal(reed) = rendered.system().aperture() else {
        panic!("a typed reed is a modal aperture");
    };
    let end = reed.modal_states()[0];
    assert!(end.displacement_m != 0.0 && end.displacement_m.abs() < 1e-3, "{end:?}");
}

fn bits(pressure: &[f64]) -> Vec<u64> { pressure.iter().map(|p| p.to_bits()).collect() }

fn with_row(row: &str) -> String {
    card(REED).replace("\ncompile_limits", &format!("\n{row}\ncompile_limits"))
}

#[test]
fn the_contraction_row_scales_only_the_jet_width_and_one_is_the_original_flow_bit_for_bit() {
    let reference = bits(&render(super::build(&card(REED), 37)).0);
    // No row and C_c = 1 are the same flow, bit for bit.
    assert_eq!(bits(&render(super::build(&with_row("contraction 1"), 37)).0), reference);
    // C_c < 1 changes the render ...
    let contracted = bits(&render(super::build(&with_row("contraction 0.611"), 37)).0);
    assert_ne!(contracted, reference);
    // ... and is exactly the physical width scaled by C_c in the flow law: the render of the card whose WIDTH is
    // C_c * 0.0056 with no contraction row, bit for bit (pressure area, swept volume, contact and limits do not
    // see the jet width).
    let scaled = REED.replace("0.000357 0.0056 ", &format!("0.000357 {} ", 0.611_f64 * 0.0056));
    assert_ne!(scaled, REED);
    assert_eq!(bits(&render(super::build(&card(&scaled), 37)).0), contracted);
    let source = super::build(&with_row("contraction 0.611"), 37);
    assert_eq!(source.reed_record().unwrap().contraction(), Some(0.611));
    assert_eq!(valve(&source).jet_contraction(), Some(0.611));
    assert_eq!(super::build(&card(REED), 37).reed_record().unwrap().contraction(), None);
}

#[test]
fn a_contraction_outside_the_unit_interval_or_misplaced_is_refused_for_its_own_reason() {
    let range = "a typed reed's contraction needs a coefficient C_c in (0, 1]";
    for (row, expected) in [
        ("contraction 0", range),
        ("contraction -0.5", range),
        ("contraction 1.0000001", range),
        ("contraction", "missing field"),
        ("contraction 0.6 0.7", "unexpected extra field"),
    ] {
        assert_eq!(refused(&with_row(row)).1, expected, "{row}");
    }
    // Twice, or before the reed, or without it, is not the grammar.
    let other = "unexpected record kind or order";
    assert_eq!(refused(&with_row("contraction 0.6\ncontraction 0.6")).1, other);
    assert_eq!(refused(&card(&format!("contraction 0.6\n{REED}"))).1, other);
}

#[test]
fn the_constriction_noise_row_is_a_seeded_series_pressure_the_callbacks_cannot_change() {
    let quiet = bits(&render(super::build(&card(REED), 37)).0);
    assert!(super::build(&card(REED), 37).constriction_noise().is_none());
    let text = with_row("constriction_noise 1 0.00087");
    let source = super::build(&text, 37);
    let card_read = source.constriction_noise().copied().expect("the row is retained");
    assert_eq!((card_read.seed, card_read.duct_radius_m), (1, 0.00087));
    let loud = bits(&render(source).0);
    assert_ne!(loud, quiet);
    // The noise is silent until the last sample's bore flow was forward, so the render starts as the quiet one.
    let first = loud.iter().zip(&quiet).position(|(a, b)| a != b).expect("the renders differ");
    assert!(first >= 1, "the noise spoke at the first sample, before any flow");
    // One seed is one stream on the mechanical clock whatever the callback partition ...
    for block in [1, 512] {
        assert_eq!(bits(&render(super::build(&text, block)).0), loud, "block {block}");
    }
    // ... another seed or another duct radius is another stream.
    assert_ne!(bits(&render(super::build(&with_row("constriction_noise 2 0.00087"), 37)).0), loud);
    assert_ne!(bits(&render(super::build(&with_row("constriction_noise 1 0.0012"), 37)).0), loud);
}

#[test]
fn a_constriction_noise_row_is_refused_for_its_own_reason_and_only_after_the_reed_and_its_contraction() {
    let radius = "constriction_noise duct radius must be finite and positive";
    for (row, expected) in [
        ("constriction_noise 1 0", radius),
        ("constriction_noise 1 -0.001", radius),
        ("constriction_noise", "missing field"),
        ("constriction_noise 1", "missing field"),
        ("constriction_noise 1.5 0.001", "invalid numeric field"),
        ("constriction_noise -1 0.001", "invalid numeric field"),
        ("constriction_noise 1 0.001 7", "unexpected extra field"),
    ] {
        assert_eq!(refused(&with_row(row)).1, expected, "{row}");
    }
    let other = "unexpected record kind or order";
    // The contraction row is read directly after the reed; the noise after it. Reversed or repeated is not the grammar.
    assert_eq!(refused(&with_row("constriction_noise 1 0.001\ncontraction 0.611")).1, other);
    assert_eq!(refused(&with_row("constriction_noise 1 0.001\nconstriction_noise 1 0.001")).1, other);
    let both = with_row("contraction 0.611\nconstriction_noise 1 0.00087");
    let source = super::build(&both, 37);
    assert_eq!(source.reed_record().unwrap().contraction(), Some(0.611));
    assert_eq!(source.constriction_noise().unwrap().seed, 1);
}

#[test]
fn the_typed_reed_with_the_live_pages_two_rows_stays_inside_its_bounds_at_the_recipes_pressure() {
    let text = with_row("contraction 0.611\nconstriction_noise 1 0.00087")
        .replace("audio 96000 96000 512 20", "audio 96000 96000 2400 20")
        .replace("events\t2\nevent\t0e0\t1e-3\tpressure-pa\t5e0\nevent\t4e-3\t1e-3\tpressure-pa\t0e0",
            "events\t1\nevent\t0e0\t2e-2\tpressure-pa\t7.4e3");
    assert!(text.contains("pressure-pa\t7.4e3"), "the gesture suffix was not replaced");
    let (pressure, rendered) = render(super::build(&text, 64));
    assert_eq!(pressure.len(), 2400);
    assert!(pressure.iter().all(|p| p.is_finite()));
    let LoadedAperture::Modal(reed) = rendered.system().aperture() else {
        panic!("a typed reed is a modal aperture");
    };
    let end = reed.modal_states()[0];
    assert!(end.displacement_m != 0.0 && end.displacement_m.abs() < 1e-3, "{end:?}");
}

#[test]
fn the_aperture_takes_a_contraction_once_inside_the_unit_interval_and_never_with_a_viscous_slit() {
    let source = super::build(&card(REED), 37);
    let gas = *source.renderer().system().load().duct().inlet_gas();
    let slit = crate::bernoulli_aperture::viscous::ViscousSlit::from_gas(&gas, 0.001).unwrap();
    let mut viscous = valve(&source).clone();
    viscous.set_viscous_slit(slit).unwrap();
    assert!(viscous.set_jet_contraction(0.611).is_err());
    let mut once = valve(&source).clone();
    once.set_jet_contraction(0.611).unwrap();
    assert!(once.set_jet_contraction(0.5).is_err());
    assert!(once.set_viscous_slit(slit).is_err());
    assert_eq!(once.jet_contraction(), Some(0.611));
    let mut bad = valve(&source).clone();
    for coefficient in [0.0, -0.1, 1.0000001, f64::NAN, f64::INFINITY] {
        assert!(bad.set_jet_contraction(coefficient).is_err(), "{coefficient}");
    }
    assert_eq!(bad.jet_contraction(), None);
}
