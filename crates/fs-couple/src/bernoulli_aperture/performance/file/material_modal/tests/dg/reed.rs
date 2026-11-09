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

// The live oboe pages' staple row: C_P 0.8 and S_c 4 mm^2, one measured reed's numbers (Almeida, Vergez and Causse 2007).
const STAPLE: &str = "staple_recovery 0.8 0.000004";

#[test]
fn the_staple_row_recovers_nothing_at_zero_and_a_real_recovery_changes_the_render() {
    let reference = bits(&render(super::build(&card(REED), 37)).0);
    // C_P = 0 recovers nothing: the original solve, bit for bit.
    assert_eq!(bits(&render(super::build(&with_row("staple_recovery 0 0.000004"), 37)).0), reference);
    let recovered = render(super::build(&with_row(STAPLE), 37)).0;
    assert!(recovered.iter().all(|p| p.is_finite()));
    assert_ne!(bits(&recovered), reference);
    let source = super::build(&with_row(STAPLE), 37);
    assert_eq!(source.reed_record().unwrap().staple_recovery(), Some((0.8, 4e-6)));
    assert!(valve(&source).staple_recovery().is_some());
    assert!(super::build(&card(REED), 37).reed_record().unwrap().staple_recovery().is_none());
    // Callback partition cannot change it.
    for block in [1, 512] {
        assert_eq!(bits(&render(super::build(&with_row(STAPLE), block)).0), bits(&recovered), "block {block}");
    }
}

#[test]
fn a_staple_row_is_refused_for_its_own_reason_and_only_directly_after_the_reed() {
    let range = "a typed reed's staple_recovery needs C_P in [0, 1] and a positive constriction area in m^2";
    for (row, expected) in [
        ("staple_recovery -0.1 0.000004", range),
        ("staple_recovery 1.0000001 0.000004", range),
        ("staple_recovery 0.8 0", range),
        ("staple_recovery 0.8 -0.000004", range),
        ("staple_recovery 0.8", "missing field"),
        ("staple_recovery", "missing field"),
        ("staple_recovery 0.8 0.000004 1", "unexpected extra field"),
        ("staple_recovery nan 0.000004", "physical scalars must be finite"),
    ] {
        assert_eq!(refused(&with_row(row)).1, expected, "{row}");
    }
    // Twice, or after the contraction, is not the grammar; the contraction and noise rows follow it.
    let other = "unexpected record kind or order";
    assert_eq!(refused(&with_row(&format!("{STAPLE}\n{STAPLE}"))).1, other);
    assert_eq!(refused(&with_row(&format!("contraction 0.611\n{STAPLE}"))).1, other);
    let all = with_row(&format!("{STAPLE}\ncontraction 0.611\nconstriction_noise 7 0.00087"));
    let source = super::build(&all, 37);
    assert_eq!(source.reed_record().unwrap().staple_recovery(), Some((0.8, 4e-6)));
    assert_eq!(source.reed_record().unwrap().contraction(), Some(0.611));
    assert_eq!(source.constriction_noise().unwrap().seed, 7);
}

#[test]
fn the_aperture_takes_a_staple_recovery_once_and_composes_viscosity_in_either_order() {
    use crate::bernoulli_aperture::viscous::{StapleRecovery, ViscousSlit};
    let source = super::build(&card(REED), 37);
    let gas = *source.renderer().system().load().duct().inlet_gas();
    let slit = ViscousSlit::from_gas(&gas, 0.001).unwrap();
    let recovery = StapleRecovery::new(0.8, 4e-6).unwrap();
    let mut viscous = valve(&source).clone();
    viscous.set_viscous_slit(slit).unwrap();
    viscous.set_staple_recovery(recovery).unwrap();
    let mut once = valve(&source).clone();
    once.set_staple_recovery(recovery).unwrap();
    assert!(once.set_staple_recovery(recovery).is_err());
    once.set_viscous_slit(slit).unwrap();
    assert_eq!((viscous.viscous_slit(), viscous.staple_recovery()), (Some(slit), Some(recovery)));
    assert_eq!((once.viscous_slit(), once.staple_recovery()), (Some(slit), Some(recovery)));
    assert!(once.set_jet_contraction(0.611).is_err());
    // A contraction and a staple compose, in either order.
    let mut contracted = valve(&source).clone();
    contracted.set_staple_recovery(recovery).unwrap();
    contracted.set_jet_contraction(0.611).unwrap();
    let mut other_order = valve(&source).clone();
    other_order.set_jet_contraction(0.611).unwrap();
    other_order.set_staple_recovery(recovery).unwrap();
    assert_eq!((contracted.jet_contraction(), contracted.staple_recovery()), (Some(0.611), Some(recovery)));
    assert_eq!((other_order.jet_contraction(), other_order.staple_recovery()), (Some(0.611), Some(recovery)));
    // The separated-entrance law still has no condensed modal owner.
    let separated = ViscousSlit::separated_jet_from_gas(&gas, 0.001).unwrap();
    assert!(once.set_viscous_slit(separated).is_err());
    assert_eq!(once.viscous_slit(), Some(slit));
    let mut separated_first = valve(&source).clone();
    separated_first.set_viscous_slit(separated).unwrap();
    assert!(separated_first.set_staple_recovery(recovery).is_err());
    assert_eq!(separated_first.staple_recovery(), None);
}

#[test]
fn the_staple_domain_refuses_at_the_step_not_silently() {
    // C_P (S / S_c)^2 >= 1 on forward flow is the live domain refusal. A constriction of 0.5 mm^2 against the reed's
    // 2 mm^2 rest area makes r = 4 and C_P r^2 = 12.8 at the very first forward step.
    let text = with_row("staple_recovery 0.8 0.0000005")
        .replace("audio 96000 96000 512 20", "audio 96000 96000 600 20")
        .replace("events\t2\nevent\t0e0\t1e-3\tpressure-pa\t5e0\nevent\t4e-3\t1e-3\tpressure-pa\t0e0",
            "events\t1\nevent\t0e0\t2e-3\tpressure-pa\t7.4e3");
    let source = super::build(&text, 64);
    let block = source.renderer().config().max_block;
    let mut renderer = source.into_renderer();
    let mut out = vec![0.0; 600];
    // The callback is declared at 64 samples: hand it blocks of that size, as `render` does, and keep the first refusal.
    let error = out
        .chunks_mut(block)
        .find_map(|chunk| renderer.block(chunk).err())
        .expect("the staple domain must refuse a wide jet, not render it");
    assert!(format!("{error:?}").contains("staple pressure recovery left its domain"), "{error:?}");
}

#[test]
fn the_live_pages_three_rows_render_inside_the_bounds_with_closed_work_and_flow_balances() {
    // The step refuses an unresolved work or flow balance, so a completed render has closed both at every sample.
    let text = with_row(&format!("{STAPLE}\ncontraction 0.611\nconstriction_noise 7 0.00087"))
        .replace("audio 96000 96000 512 20", "audio 96000 96000 2400 20")
        .replace("events\t2\nevent\t0e0\t1e-3\tpressure-pa\t5e0\nevent\t4e-3\t1e-3\tpressure-pa\t0e0",
            "events\t1\nevent\t0e0\t2e-2\tpressure-pa\t7.4e3");
    assert!(text.contains("pressure-pa\t7.4e3"), "the gesture suffix was not replaced");
    let (pressure, rendered) = render(super::build(&text, 64));
    assert_eq!(pressure.len(), 2400);
    assert!(pressure.iter().all(|p| p.is_finite()));
    let frame = &rendered.last_frame().expect("an accepted frame").aperture;
    let flows = frame.jet_flow_m3_s.abs() + frame.swept_flow_m3_s.abs() + frame.bore_flow_m3_s.abs();
    assert!(frame.flow_residual_m3_s.abs() <= 1e-9 * flows.max(1e-30), "{frame:?}");
    assert!(frame.dissipated_energy_j.is_finite());
    let LoadedAperture::Modal(reed) = rendered.system().aperture() else {
        panic!("a typed reed is a modal aperture");
    };
    let end = reed.modal_states()[0];
    assert!(end.displacement_m != 0.0 && end.displacement_m.abs() < 1e-3, "{end:?}");
}

// The oboe pages' tongue: a 2 mm x 5.6 mm strip (1.12e-5 m^2), the reed's own face k H / P_c = 6.979e-5 m^2 as the effective
// area, 0.2 N, held until 10 ms and let go over 5 ms. The 25 ms card below runs past the release.
const TONGUE: &str = "tongue_contact 1.12e-5 6.979e-5 0.2 0.01 0.005";

fn long(text: String) -> String {
    let text = text
        .replace("audio 96000 96000 512 20", "audio 96000 96000 2400 20")
        .replace("events\t2\nevent\t0e0\t1e-3\tpressure-pa\t5e0\nevent\t4e-3\t1e-3\tpressure-pa\t0e0",
            "events\t1\nevent\t0e0\t2e-2\tpressure-pa\t7.4e3");
    assert!(text.contains("pressure-pa\t7.4e3"), "the gesture suffix was not replaced");
    text
}

/// The refusal of a card by any owner, as text.
fn refused_by_anyone(text: &str) -> String {
    match MaterialModalDgValvePerformance::from_bytes_with_sources(
        text.as_bytes(), 37, &CancelGate::new_clock_free(), None,
        |_| Ok(pack(900.0, false)),
    ) {
        Err(error) => format!("{error:?}"),
        Ok(_) => panic!("a refused tongue row was accepted"),
    }
}

/// The reed's displacement [m] after its first `n` samples (one block of `n`).
fn displacement_after(text: &str, n: usize) -> f64 {
    let mut renderer = super::build(text, n).into_renderer();
    let mut out = vec![0.0; n];
    renderer.block(&mut out).unwrap();
    let LoadedAperture::Modal(reed) = renderer.system().aperture() else {
        panic!("a typed reed is a modal aperture");
    };
    reed.modal_states()[0].displacement_m
}

#[test]
fn the_tongue_row_opens_one_force_port_on_the_reed_and_its_absence_opens_none() {
    let plain = super::build(&long(card(REED)), 37);
    assert!(plain.reed_record().unwrap().tongue_contact().is_none());
    assert_eq!(valve(&plain).spec().modes[0].force_projection, [0.0, 0.0]);
    let source = super::build(&long(with_row(TONGUE)), 37);
    assert!(source.reed_record().unwrap().tongue_contact().is_some());
    assert_eq!(valve(&source).spec().modes[0].force_projection, [-1.0, 0.0]);
    // The row changes nothing of the mode itself.
    let (a, b) = (valve(&plain).spec().modes[0], valve(&source).spec().modes[0]);
    assert_eq!(
        (a.mass_kg.to_bits(), a.stiffness_n_m.to_bits(), a.damping_n_s_m.to_bits(), a.pressure_area_m2.to_bits(),
            a.gap_projection.to_bits()),
        (b.mass_kg.to_bits(), b.stiffness_n_m.to_bits(), b.damping_n_s_m.to_bits(), b.pressure_area_m2.to_bits(),
            b.gap_projection.to_bits()),
    );
    let (_, without) = render(plain);
    assert!(without.applied_plate_force_controls(0).is_empty());
    let (_, with) = render(source);
    let events = with.applied_plate_force_controls(0);
    assert!(events.len() > 2 && events.windows(2).all(|w| w[0].sample < w[1].sample));
    assert!(events.iter().all(|e| e.port == 0 && e.force_n >= 0.0), "a closing force is positive newtons on the port");
    assert_eq!(events[0].sample, 0);
    // The tongue starts holding at least its own pull (modal factor < 1 times 0.2 N, plus the pressure it must carry)
    // and lets go completely by the end of the release, 15 ms = 1440 mechanical samples.
    assert!(events[0].force_n > 0.1, "{:?}", events[0]);
    let last = events.last().unwrap();
    assert!(last.force_n == 0.0 && last.sample >= 1439, "{last:?}");
    assert!(with.plate_force_ports(0).is_empty());
}

#[test]
fn the_tongue_presses_the_reed_closed_until_it_lets_go_and_the_callbacks_cannot_change_it() {
    let control = long(card(REED));
    let held = long(with_row(TONGUE));
    // 8 mechanical samples (83 us) into the note the control's reed has barely felt the blowing pressure (about 2 mN); the
    // tongue's opening load is about 0.44 N (its 0.2 N times the modal factor, plus the pressure it must carry until 10 ms),
    // which has already pushed the reed toward its lay.
    let open = displacement_after(&control, 8);
    let pressed = displacement_after(&held, 8);
    assert!(pressed < open - 2e-5, "tongue {pressed} m, control {open} m");
    let (plain, _) = render(super::build(&control, 64));
    let (tongued, rendered) = render(super::build(&held, 64));
    assert!(tongued.iter().all(|p| p.is_finite()) && tongued.iter().any(|p| p.abs() > 1e-15));
    assert_ne!(bits(&plain), bits(&tongued));
    for block in [1, 512] {
        assert_eq!(bits(&render(super::build(&held, block)).0), bits(&tongued), "block {block}");
    }
    let LoadedAperture::Modal(reed) = rendered.system().aperture() else {
        panic!("a typed reed is a modal aperture");
    };
    let end = reed.modal_states()[0];
    assert!(end.displacement_m.abs() < 1e-3, "{end:?}");
}

#[test]
fn a_tongue_row_is_refused_for_its_own_reason_and_only_after_the_contraction() {
    let positive = "a typed reed's tongue_contact needs a positive area, effective area, force and release time";
    for (row, expected) in [
        ("tongue_contact 0 6.979e-5 0.2 0.01 0.005", positive),
        ("tongue_contact 1.12e-5 -6.979e-5 0.2 0.01 0.005", positive),
        ("tongue_contact 1.12e-5 6.979e-5 0 0.01 0.005", positive),
        ("tongue_contact 1.12e-5 6.979e-5 0.2 0.01 0", positive),
        ("tongue_contact 1.12e-4 6.979e-5 0.2 0.01 0.005", "tongue_contact area exceeds the vibrating reed blade area"),
        ("tongue_contact 1.12e-5 6.979e-5 0.2 -0.01 0.005", "tongue_contact requires a nonnegative release_start_s"),
        ("tongue_contact 1.12e-5 6.979e-5 0.2 0.01", "missing field"),
        ("tongue_contact", "missing field"),
        ("tongue_contact 1.12e-5 6.979e-5 0.2 0.01 0.005 1", "unexpected extra field"),
        ("tongue_contact nan 6.979e-5 0.2 0.01 0.005", "physical scalars must be finite"),
    ] {
        assert_eq!(refused(&long(with_row(row))).1, expected, "{row}");
    }
    // A release after the last observed control tick is the performance's refusal, not the card grammar's.
    let late = refused_by_anyone(&long(with_row("tongue_contact 1.12e-5 6.979e-5 0.2 0.05 0.005")));
    assert!(late.contains("tongue_contact release starts after the last observed control tick"), "{late}");
    // The row follows the contraction and precedes the noise: reversed, repeated or before the reed is not the grammar.
    let other = "unexpected record kind or order";
    assert_eq!(refused(&long(with_row(&format!("{TONGUE}\ncontraction 0.611")))).1, other);
    assert_eq!(refused(&long(with_row(&format!("{TONGUE}\n{STAPLE}")))).1, other);
    assert_eq!(refused(&long(with_row(&format!("{TONGUE}\n{TONGUE}")))).1, other);
    assert_eq!(refused(&long(with_row(&format!("constriction_noise 7 0.00087\n{TONGUE}")))).1, other);
    assert_eq!(refused(&long(card(&format!("{TONGUE}\n{REED}")))).1, other);
    let all = long(with_row(&format!("{STAPLE}\ncontraction 0.611\n{TONGUE}\nconstriction_noise 7 0.00087")));
    let source = super::build(&all, 37);
    let record = source.reed_record().unwrap();
    assert_eq!(record.staple_recovery(), Some((0.8, 4e-6)));
    assert_eq!(record.contraction(), Some(0.611));
    assert!(record.tongue_contact().is_some());
    assert_eq!(source.constriction_noise().unwrap().seed, 7);
}

#[test]
fn the_tongue_with_the_live_pages_other_rows_stays_inside_the_bounds_with_closed_balances() {
    let text = long(with_row(&format!("{STAPLE}\ncontraction 0.611\n{TONGUE}\nconstriction_noise 7 0.00087")));
    let (pressure, rendered) = render(super::build(&text, 64));
    assert_eq!(pressure.len(), 2400);
    assert!(pressure.iter().all(|p| p.is_finite()));
    let frame = &rendered.last_frame().expect("an accepted frame").aperture;
    let flows = frame.jet_flow_m3_s.abs() + frame.swept_flow_m3_s.abs() + frame.bore_flow_m3_s.abs();
    assert!(frame.flow_residual_m3_s.abs() <= 1e-9 * flows.max(1e-30), "{frame:?}");
    let LoadedAperture::Modal(reed) = rendered.system().aperture() else {
        panic!("a typed reed is a modal aperture");
    };
    let end = reed.modal_states()[0];
    assert!(end.displacement_m != 0.0 && end.displacement_m.abs() < 1e-3, "{end:?}");
    // The same card without the tongue is another note.
    let plain = long(with_row(&format!("{STAPLE}\ncontraction 0.611\nconstriction_noise 7 0.00087")));
    assert_ne!(bits(&render(super::build(&plain, 64)).0), bits(&pressure));
}

#[test]
fn a_tongue_is_a_typed_reed_force_and_a_plate_source_refuses_it() {
    use crate::render::schedule::reed::TongueContactInput;
    let tongue = TongueContactInput::new(1.12e-5, 6.979e-5, 0.2, 0.001, 0.005).unwrap();
    let plate = super::build(&super::input(), 37).into_renderer();
    match plate.with_reed_closing_force(&tongue, 6.979e-5) {
        Err(error) => assert!(format!("{error:?}").contains("a tongue is a typed-reed force"), "{error:?}"),
        Ok(_) => panic!("a plate source accepted a tongue"),
    }
}

#[test]
fn the_dg_tongue_force_events_are_the_shared_series_on_the_same_numbers() {
    use crate::render::schedule::reed::TongueContactInput;
    use fs_scenario::gesture::GestureSchedule;
    // The long card's own gesture suffix: 48 kHz control ticks, one 20 ms ramp to 7.4 kPa on the track `mouth`.
    let schedule = GestureSchedule::from_canonical_bytes(
        "frankensim-gesture-schedule-v1\ncontrol_rate_hz\t48000\ntracks\t1\ntrack\tmouth\tblowing-pressure\n\
         initial\tpressure-pa\t0e0\nevents\t1\nevent\t0e0\t2e-2\tpressure-pa\t7.4e3\n"
            .as_bytes(),
    )
    .expect("the long card's gesture suffix");
    let source = super::build(&long(with_row(TONGUE)), 37);
    let record = source.reed_record().unwrap();
    let face = record.pressure_area_m2();
    let tongue: TongueContactInput = *record.tongue_contact().unwrap();
    // 2400 audio samples at 96 kHz against 48 kHz control ticks: the last observed tick is (2400 - 1) / 2 = 1199.
    let series = tongue.force_series(&schedule, "mouth", face, 1199).unwrap();
    assert!(series.len() > 2 && series.last().unwrap().1 == 0.0);
    let (_, rendered) = render(source);
    let events = rendered.applied_plate_force_controls(0);
    assert_eq!(events.len(), series.len());
    for (event, (tick, force_n)) in events.iter().zip(&series) {
        assert_eq!(event.sample, tick * 2, "two audio samples per control tick");
        assert_eq!(event.port, 0);
        assert_eq!(event.force_n.to_bits(), force_n.to_bits());
    }
}
