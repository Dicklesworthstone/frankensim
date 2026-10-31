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
