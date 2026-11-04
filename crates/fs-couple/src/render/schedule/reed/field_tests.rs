use super::*;
use crate::driving_point::regional::field::tests::field;

const INPUT_FIXTURE: &str = include_str!("../../../../examples/reed-flanged-baffled.performance");
/// The fixture's physics without its sidecar declaration: these library tests
/// add their own observation records, and the library never refuses an
/// undeclared law (ruling 6532 puts that in the sidecar writer).
static INPUT: std::sync::LazyLock<String> =
    std::sync::LazyLock::new(|| INPUT_FIXTURE.replace("observer_law stokes-kirchhoff\n", ""));

fn source() -> String {
    INPUT
        .replace("audio 48000 12000 2", "audio 48000 4801 2")
        .replace("ambient mixture 293.15 100000 mole 2 N2 0.7 O2 0.3", "ambient 300 101325 0")
        .replace("reed 0.0004 0.013 6000 0 0 0.35", "reed 0.0004 0.013 6000 0.00001 100 0.35")
        .replace("segments 2\ncylinder 0.0015 0.2\ncylinder 0.0015 0.3 gas mixture 293.15 100000 mole 3 N2 0.4 O2 0.2 CO2 0.4",
            "segments 1\ncylinder 0.0015 0.57")
}
fn load(text: &str, temperature: f64, cells: usize, block: usize) -> ReedPerformance {
    ReedPerformance::from_bytes_with_gas_field(text.as_bytes(), block,
        ObservationPolicy::Diagnostic, None, field(temperature), cells).unwrap()
}
fn voice(p: &ReedPerformance) -> &ReedBoreVoice {
    let RenderVoice::ReedBore(voice) = &p.renderer.context().voices[0] else { panic!("reed voice") };
    voice
}
fn render(text: &str, temperature: f64, block: usize) -> (Vec<f64>, [u64; 4]) {
    let mut p = load(text, temperature, 8, block);
    let mut pressure = vec![0.0; p.info.samples as usize];
    for chunk in pressure.chunks_mut(block) { p.pressure_block(chunk).unwrap(); }
    assert_eq!(p.remaining_samples(), 0);
    assert!(pressure.iter().all(|p| p.is_finite()));
    let v = voice(&p);
    (pressure, [v.reed_y, v.reed_v, v.p_plus_prev, v.line.incoming()].map(f64::to_bits))
}

#[test]
fn exact_field_boundaries_reach_reed_feedback_and_physical_microphone() {
    let p = load(&source(), 280.0, 8, 37);
    let (field, cells) = p.gas_field().unwrap();
    assert_eq!(cells, 8);
    assert_eq!(p.segment_gases().len(), 8);
    let inlet = *field.state_at(0.0).unwrap().state();
    let outlet = *field.state_at(field.length_m()).unwrap().state();
    assert_eq!(voice(&p).rho, inlet.density);
    assert_eq!(voice(&p).inlet_gas.sound_speed, inlet.sound_speed);
    assert_ne!(p.segment_gases()[0].temperature, inlet.temperature);
    assert_eq!(p.outlet_receiver().unwrap().gas, outlet);
    assert_ne!(p.segment_gases().last().unwrap(), &outlet);
    assert_eq!(p.info.segments, 1, "source geometry count is not the refinement count");
}

#[test]
fn downstream_field_changes_massive_reed_motion_and_microphone_with_fixed_exterior() {
    let text = source().replace("termination flanged", "exterior 300 101325 0\ntermination flanged");
    let a = load(&text, 300.0, 8, 37);
    let b = load(&text, 280.0, 8, 37);
    assert_eq!(voice(&a).rho, voice(&b).rho);
    assert_eq!(voice(&a).zc, voice(&b).zc);
    assert_eq!(a.outlet_receiver().unwrap().gas, b.outlet_receiver().unwrap().gas);
    assert_eq!(a.info.compiled_controls, b.info.compiled_controls);
    let (pa, ma) = render(&text, 300.0, 37);
    let (pb, mb) = render(&text, 280.0, 37);
    let norm: f64 = pa.iter().map(|p| p*p).sum();
    let difference: f64 = pa.iter().zip(&pb).map(|(p,q)| (p-q)*(p-q)).sum();
    assert!(norm > 0.0 && difference > 1e-8*norm);
    assert_ne!(&ma[..2], &mb[..2], "downstream field must react back on the reed");
    assert_eq!((pb, mb), render(&text, 280.0, 113));
}

#[test]
fn scheduled_field_keeps_original_compact_voice_at_every_sample() {
    let text = source().replace("baffled_outlet 0 0 0.5 4 16 1000", "listener_m 1");
    let mut p = load(&text, 285.0, 8, 37);
    let geometry = Duct { segments: vec![Segment::Cylinder { radius: 0.0015, length: 0.57 }] };
    let reed = BeatingReed {
        rest_opening_m: 0.0004, width_m: 0.013, closing_pressure_pa: 6000.0,
        mass_kg: 0.00001, stiffness_n_m: 100.0, damping_ratio: 0.35,
        blowing_pressure_pa: 0.0, attack_s: 0.0,
    };
    let mut direct = ReedBoreVoice::new_with_gas_field(&geometry, &field(285.0), 8,
        reed, Termination::FlangedOpen, PlateBank::default(), 1.0, 48000, 4801, None).unwrap();
    let schedule = GestureSchedule::from_canonical_bytes(
        text.split_once("\nschedule\n").unwrap().1.as_bytes()).unwrap();
    let ambient = *p.ambient_gas();
    let mut observer = p.reference_observer(&ambient);
    for n in 0..4801 {
        direct.set_blowing_pressure(schedule.sample("blow", n*700/48000).unwrap());
        let mut a = [0.0]; let mut b = [0.0];
        p.block(&mut a).unwrap();
        direct.step_block(&mut b).unwrap();
        let parts = direct.last_observer_parts();
        assert_eq!(a[0].to_bits(), (parts.internal_pa + observer.process(parts.radiated_pa)).to_bits());
        assert_eq!(voice(&p).reed_y.to_bits(), direct.reed_y.to_bits());
        assert_eq!(voice(&p).reed_v.to_bits(), direct.reed_v.to_bits());
    }
}

#[test]
fn uniform_field_keeps_original_segment_count_and_refuses_conflicting_authority() {
    let text = source();
    let p = load(&text, 300.0, 64, 37);
    assert_eq!(p.segment_gases().len(), 1);
    for bad in [
        text.replace("cylinder 0.0015 0.57", "cylinder 0.0015 0.57 gas 300 0"),
        text.replace("ambient 300 101325 0", "ambient 300 100000 0"),
        text.replace("cylinder 0.0015 0.57", "cylinder 0.0015 0.56"),
    ] {
        assert!(ReedPerformance::from_bytes_with_gas_field(bad.as_bytes(), 37,
            ObservationPolicy::Diagnostic, None, field(285.0), 8).is_err());
    }
    for cells in [0, 65, usize::MAX] {
        assert!(ReedPerformance::from_bytes_with_gas_field(text.as_bytes(), 37,
            ObservationPolicy::Diagnostic, None, field(285.0), cells).is_err());
    }
}

#[test]
fn continuous_field_uses_the_scheduled_source_and_same_physical_microphone() {
    let text = source().replace("termination flanged", "exterior 300 101325 0\ntermination flanged");
    let make = |temperature| ReedPerformance::from_bytes_with_field_resolution(
        text.as_bytes(), 113, ObservationPolicy::Diagnostic, None, field(temperature),
        GasFieldResolution::Continuous(fs_duct::RiccatiSteps(16)),
    ).unwrap();
    let (mut a, mut b) = (make(300.0), make(280.0));
    assert_eq!(a.gas_field_resolution(), b.gas_field_resolution());
    assert_eq!(b.gas_field_resolution(), Some(
        GasFieldResolution::Continuous(fs_duct::RiccatiSteps(16)),
    ));
    assert_eq!(b.segment_gases().len(), 1);
    assert_eq!(voice(&a).rho, voice(&b).rho);
    assert_eq!(voice(&a).zc, voice(&b).zc);
    assert_eq!(a.outlet_receiver().unwrap().gas, b.outlet_receiver().unwrap().gas);
    let mut pa = vec![0.0; a.info.samples as usize];
    let mut pb = vec![0.0; b.info.samples as usize];
    for (a_block, b_block) in pa.chunks_mut(113).zip(pb.chunks_mut(113)) {
        a.pressure_block(a_block).unwrap();
        b.pressure_block(b_block).unwrap();
    }
    assert_eq!(a.remaining_samples(), 0);
    assert_eq!(b.remaining_samples(), 0);
    let norm: f64 = pa.iter().map(|p| p * p).sum();
    let difference: f64 = pa.iter().zip(&pb).map(|(p, q)| (p - q).powi(2)).sum();
    assert!(pa.iter().chain(&pb).all(|p| p.is_finite()));
    assert!(norm > 0.0 && difference > 1e-8 * norm);
    assert_ne!(voice(&a).reed_y.to_bits(), voice(&b).reed_y.to_bits());
}

// ---- `open_hole_interaction tmmi EPSILON` (opt-in radiation coupling of the open holes and the bell) ----

/// The fixture bore with two open holes 47.5 mm apart; the lengths sum exactly to the field length.
fn holed(text: &str) -> String {
    text.replace("termination flanged", "exterior 300 101325 0\ntermination flanged")
        .replace("segments 1\ncylinder 0.0015 0.57",
            "segments 5\ncylinder 0.0015 0.285\nhole 0.0013 0.0005 0.0015 1\ncylinder 0.0015 0.0475\n\
             hole 0.0013 0.0005 0.0015 1\ncylinder 0.0015 0.2375")
}
fn with_row(text: &str, row: &str) -> String {
    text.replace("compile_limits", &format!("{row}\ncompile_limits"))
}
fn continuous_load(text: &str, outlet_temperature: f64) -> Result<ReedPerformance, ReedPerformanceError> {
    ReedPerformance::from_bytes_with_field_resolution(
        text.as_bytes(), 113, ObservationPolicy::Diagnostic, None, field(outlet_temperature),
        GasFieldResolution::Continuous(fs_duct::RiccatiSteps(16)),
    )
}
fn pressure_of(text: &str) -> Vec<f64> {
    let mut p = continuous_load(text, 280.0).unwrap();
    let mut out = vec![0.0; p.info.samples as usize];
    for chunk in out.chunks_mut(113) { p.pressure_block(chunk).unwrap(); }
    assert!(out.iter().all(|x| x.is_finite()));
    out
}
fn distance(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| (x - y).powi(2)).sum::<f64>().sqrt()
}

#[test]
fn the_open_hole_interaction_row_is_off_by_default_and_reaches_the_rendered_pressure() {
    let base = holed(&source());
    let plain = pressure_of(&base);
    assert_eq!(plain, pressure_of(&base));
    let norm = plain.iter().map(|x| x * x).sum::<f64>().sqrt();
    let coupled = pressure_of(&with_row(&base, "open_hole_interaction tmmi 2"));
    let weak = pressure_of(&with_row(&base, "open_hole_interaction tmmi 1e15"));
    let (dc, dw) = (distance(&plain, &coupled), distance(&plain, &weak));
    assert!(norm > 0.0 && dc > 1.0e-9 * norm, "coupled differs: {dc} of {norm}");
    // Vanishing coupling reproduces the independent holes far more closely than eps = 2 does.
    assert!(dw < 0.1 * dc, "decoupled {dw} vs coupled {dc}");
}

#[test]
fn the_open_hole_interaction_row_is_validated_and_needs_a_continuous_binding() {
    let base = holed(&source());
    for bad in [
        "open_hole_interaction tmmi 0", "open_hole_interaction tmmi -2", "open_hole_interaction tmmi nan",
        "open_hole_interaction tmmi inf", "open_hole_interaction tmmi", "open_hole_interaction tmmi 2 3",
        "open_hole_interaction foo 2", "open_hole_interaction",
    ] {
        assert!(continuous_load(&with_row(&base, bad), 280.0).is_err(), "{bad}");
    }
    let declared = with_row(&base, "open_hole_interaction tmmi 2");
    assert!(continuous_load(&declared, 280.0).is_ok());
    // A midpoint binding has no continuous operator to couple; a uniform field gives the regional
    // wrapper no continuous binding either.
    assert!(ReedPerformance::from_bytes_with_gas_field(declared.as_bytes(), 113,
        ObservationPolicy::Diagnostic, None, field(280.0), 8).is_err());
    assert!(continuous_load(&declared, 300.0).is_err());
}

// ---- `tonehole_lengths lefebvre-scavone-2012` (opt-in printed thick-wall hole lengths) ----

#[test]
fn the_tonehole_lengths_row_is_off_by_default_and_reaches_the_rendered_pressure() {
    let base = holed(&source());
    let plain = pressure_of(&base);
    let norm = plain.iter().map(|x| x * x).sum::<f64>().sqrt();
    let printed = pressure_of(&with_row(&base, "tonehole_lengths lefebvre-scavone-2012"));
    let both = pressure_of(&with_row(&with_row(&base, "open_hole_interaction tmmi 2"),
        "tonehole_lengths lefebvre-scavone-2012"));
    let coupled = pressure_of(&with_row(&base, "open_hole_interaction tmmi 2"));
    assert!(norm > 0.0 && distance(&plain, &printed) > 1.0e-9 * norm);
    assert!(distance(&printed, &both) > 1.0e-9 * norm, "the coupling still acts with the printed lengths");
    assert!(distance(&coupled, &both) > 1.0e-9 * norm, "the printed lengths still act with the coupling");
    assert_eq!(plain, pressure_of(&base));
}

#[test]
fn the_tonehole_lengths_row_is_validated_ordered_and_needs_a_continuous_binding() {
    let base = holed(&source());
    for bad in [
        "tonehole_lengths", "tonehole_lengths foo", "tonehole_lengths none",
        "tonehole_lengths lefebvre-scavone-2012 extra",
    ] {
        assert!(continuous_load(&with_row(&base, bad), 280.0).is_err(), "{bad}");
    }
    let declared = with_row(&base, "tonehole_lengths lefebvre-scavone-2012");
    assert!(continuous_load(&declared, 280.0).is_ok());
    // After the interaction row it is fine; before it, the header order refuses.
    let ordered = with_row(&with_row(&base, "open_hole_interaction tmmi 2"), "tonehole_lengths lefebvre-scavone-2012");
    assert!(continuous_load(&ordered, 280.0).is_ok());
    let swapped = with_row(&with_row(&base, "tonehole_lengths lefebvre-scavone-2012"), "open_hole_interaction tmmi 2");
    assert!(continuous_load(&swapped, 280.0).is_err());
    assert!(ReedPerformance::from_bytes_with_gas_field(declared.as_bytes(), 113,
        ObservationPolicy::Diagnostic, None, field(280.0), 8).is_err());
    assert!(continuous_load(&declared, 300.0).is_err());
}

// ---- `constriction_noise SEED DUCT_RADIUS_M` (the brass lip row's duct-constriction noise on the reed channel) ----

#[test]
fn the_constriction_noise_row_is_off_by_default_seeded_and_reaches_the_rendered_pressure() {
    let base = holed(&source());
    let plain = pressure_of(&base);
    let norm = plain.iter().map(|x| x * x).sum::<f64>().sqrt();
    let row = "constriction_noise 7 0.0015";
    let noisy = pressure_of(&with_row(&base, row));
    assert!(norm > 0.0 && distance(&plain, &noisy) > 1.0e-9 * norm, "the row must reach the pressure");
    let rms = |v: &[f64]| (v.iter().map(|x| x * x).sum::<f64>() / v.len() as f64).sqrt();
    assert!(rms(&noisy) > 0.3 * rms(&plain) && rms(&noisy) < 3.0 * rms(&plain),
        "a series source of K times the blowing pressure leaves the level near the plain note");
    assert_eq!(noisy, pressure_of(&with_row(&base, row)), "one seed, one stream");
    assert_ne!(noisy, pressure_of(&with_row(&base, "constriction_noise 8 0.0015")), "another seed, another stream");
    assert_eq!(plain, pressure_of(&base), "no row, the original voice");
}

#[test]
fn the_constriction_noise_row_is_validated_ordered_and_once_only() {
    let base = holed(&source());
    for bad in [
        "constriction_noise", "constriction_noise 7", "constriction_noise 7 0", "constriction_noise 7 -0.0015",
        "constriction_noise 7 nan", "constriction_noise 7 inf", "constriction_noise -3 0.0015",
        "constriction_noise 7 0.0015 extra", "constriction_noise x 0.0015",
    ] {
        assert!(continuous_load(&with_row(&base, bad), 280.0).is_err(), "{bad}");
    }
    let declared = with_row(&base, "constriction_noise 7 0.0015");
    assert!(continuous_load(&declared, 280.0).is_ok());
    // After the hole rows it is fine; before them, the header order refuses.
    let ordered = with_row(&with_row(&with_row(&base, "open_hole_interaction tmmi 2"), "tonehole_lengths lefebvre-scavone-2012"),
        "constriction_noise 7 0.0015");
    assert!(continuous_load(&ordered, 280.0).is_ok());
    let swapped = with_row(&with_row(&base, "constriction_noise 7 0.0015"), "open_hole_interaction tmmi 2");
    assert!(continuous_load(&swapped, 280.0).is_err());
    // Once only.
    assert!(continuous_load(&with_row(&declared, "constriction_noise 8 0.0015"), 280.0).is_err());
}
