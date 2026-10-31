use super::*;
use fs_material::gas_species::{GasSpecies, resolve_gas_mixture};
use fs_qty::semantic::{Composition, CompositionBasis};

const INPUT: &str = include_str!("../../../../examples/reed-regional-mixture.performance");
const DEFAULT: &str = "mixture 293.15 100000 mole 2 N2 0.7 O2 0.3";
const ENRICHED: &str = "mixture 293.15 100000 mole 3 N2 0.4 O2 0.2 CO2 0.4";

fn load(text: &str, block: usize) -> ReedPerformance {
    ReedPerformance::from_bytes(text.as_bytes(), block).unwrap()
}
fn voice(p: &ReedPerformance) -> &ReedBoreVoice {
    let RenderVoice::ReedBore(voice) = &p.renderer.context().voices[0] else { panic!("expected reed") };
    voice
}
fn motion(p: &ReedBoreVoice) -> [u64; 5] {
    [p.reed_y, p.reed_v, p.p_plus_prev, p.f_jet_prev, p.line.incoming()].map(f64::to_bits)
}
fn render(text: &str, block: usize) -> (Vec<f64>, [u64; 5]) {
    let mut p = load(text, block);
    let mut out = vec![0.0; p.info().samples as usize];
    for part in out.chunks_mut(block) { p.block(part).unwrap(); }
    assert_eq!(p.remaining_samples(), 0);
    assert!(out.iter().all(|p| p.is_finite()));
    (out, motion(voice(&p)))
}
fn geometry() -> Duct {
    Duct { segments: vec![Segment::Cylinder { radius: 0.0015, length: 0.2 },
        Segment::Cylinder { radius: 0.0015, length: 0.3 }] }
}
fn reed(mass: f64) -> BeatingReed {
    BeatingReed { rest_opening_m: 0.0004, width_m: 0.013, closing_pressure_pa: 6000.0,
        mass_kg: mass, stiffness_n_m: if mass > 0.0 { 100.0 } else { 0.0 }, damping_ratio: 0.35,
        blowing_pressure_pa: 0.0, attack_s: 0.0 }
}
fn source(text: &str, mass: f64) -> String {
    text.replace("reed 0.0004 0.013 6000 0 0 0.35", &format!(
        "reed 0.0004 0.013 6000 {mass} {} 0.35", reed(mass).stiffness_n_m))
}

#[test]
fn regional_input_matches_direct_material_and_voice_construction_at_every_sample() {
    let make = |species: &[GasSpecies], x| resolve_gas_mixture(species,
        &Composition::new(CompositionBasis::MoleFraction, x).unwrap(), 293.15, 100000.0).unwrap();
    let gases = [make(&[GasSpecies::Nitrogen, GasSpecies::Oxygen], vec![0.7,0.3]),
        make(&[GasSpecies::Nitrogen, GasSpecies::Oxygen, GasSpecies::CarbonDioxide], vec![0.4,0.2,0.4])];
    for mass in [0.0, 1e-5] {
        let text = source(INPUT, mass);
        let mut actual = load(&text, 37);
        assert_eq!(actual.segment_gases(), &gases);
        assert_eq!(*actual.ambient_gas(), gases[0]);
        let mut direct = ReedBoreVoice::new_with_segment_gases(&geometry(), &gases, reed(mass),
            Termination::UnflangedOpen, PlateBank::default(), 1.0, 48000, 4801, None).unwrap();
        let schedule = GestureSchedule::from_canonical_bytes(text.split_once("\nschedule\n").unwrap().1.as_bytes()).unwrap();
        let mut observer = actual.reference_observer(&gases[0]);
        for n in 0..4801 {
            direct.set_blowing_pressure(schedule.sample("blow", n*700/48000).unwrap());
            let mut a = [0.0]; let mut b = [0.0];
            actual.block(&mut a).unwrap(); direct.step_block(&mut b).unwrap();
            let parts = direct.last_observer_parts();
            assert_eq!(a[0].to_bits(), (parts.internal_pa + observer.process(parts.radiated_pa)).to_bits());
            assert_eq!(motion(voice(&actual)), motion(&direct));
        }
        assert_eq!(render(&text, 37), render(&text, 113));
    }
}

#[test]
fn uniform_explicit_regions_preserve_the_original_reed_initialization_and_runtime() {
    let uniform = INPUT.replace(&format!(" gas {ENRICHED}"), "")
        .replace(DEFAULT, "293.15 100000 0")
        .replace("cylinder 0.0015 0.2", "cylinder 0.0015 0.2 gas 293.15 0")
        .replace("cylinder 0.0015 0.3", "cylinder 0.0015 0.3 gas 293.15 0");
    for mass in [0.0, 1e-5] {
        let text = source(&uniform, mass);
        let mut actual = load(&text, 37);
        let gas = GasState::try_new_moist_air(293.15, 100000.0, 0.0).unwrap();
        let mut original = ReedBoreVoice::new(&geometry(), &gas, reed(mass),
            Termination::UnflangedOpen, PlateBank::default(), 1.0, 48000, 4801, None).unwrap();
        assert_eq!(motion(voice(&actual)), motion(&original));
        let schedule = GestureSchedule::from_canonical_bytes(text.split_once("\nschedule\n").unwrap().1.as_bytes()).unwrap();
        let mut observer = actual.reference_observer(&gas);
        for n in 0..4801 {
            original.set_blowing_pressure(schedule.sample("blow", n*700/48000).unwrap());
            let mut a = [0.0]; let mut b = [0.0];
            actual.block(&mut a).unwrap(); original.step_block(&mut b).unwrap();
            let parts = original.last_observer_parts();
            assert_eq!(a[0].to_bits(), (parts.internal_pa + observer.process(parts.radiated_pa)).to_bits());
            assert_eq!(motion(voice(&actual)), motion(&original));
        }
    }
}

#[test]
fn downstream_only_gas_change_reacts_back_on_massive_reed_motion() {
    let text = source(INPUT, 1e-5);
    let uniform = text.replace(&format!(" gas {ENRICHED}"), "");
    let mut a = load(&text, 37); let mut b = load(&uniform, 37);
    assert_eq!(voice(&a).rho, voice(&b).rho);
    assert_eq!(voice(&a).zc, voice(&b).zc);
    assert_eq!(a.info().compiled_controls, b.info().compiled_controls);
    let mut changed_motion = false;
    let (mut energy, mut difference) = (0.0, 0.0);
    for _ in 0..4801 {
        let mut pa = [0.0]; let mut pb = [0.0];
        a.block(&mut pa).unwrap(); b.block(&mut pb).unwrap();
        energy += pb[0]*pb[0]; difference += (pa[0]-pb[0])*(pa[0]-pb[0]);
        changed_motion |= voice(&a).reed_y != voice(&b).reed_y || voice(&a).reed_v != voice(&b).reed_v;
    }
    assert!(changed_motion, "downstream medium must change mechanics, not just observer gain");
    assert!(energy > 0.0 && difference > 1e-8*energy);
}

#[test]
fn actual_terminal_medium_not_inlet_medium_decides_the_radiation_band() {
    let text = INPUT.replace(DEFAULT, "293.15 100000 0")
        .replace(&format!(" gas {ENRICHED}"), "")
        .replace("cylinder 0.0015 0.2", "cylinder 0.0022 0.2 gas mixture 293.15 100000 mole 1 CO2 1")
        .replace("cylinder 0.0015 0.3", "cylinder 0.0022 0.3");
    let actual = load(&text, 37);
    assert_ne!(actual.segment_gases()[0], *actual.ambient_gas());
    assert_eq!(actual.segment_gases()[1], *actual.ambient_gas(), "omission inherits ambient, not the preceding region");
    assert_eq!(voice(&actual).rho, actual.segment_gases()[0].density);
    assert_eq!(voice(&actual).inlet_gas.sound_speed, actual.segment_gases()[0].sound_speed);
    let radius = 0.0022;
    assert_eq!(voice(&actual).zc, actual.segment_gases()[0].density*actual.segment_gases()[0].sound_speed
        /(core::f64::consts::PI*radius*radius));
    // A uniform CO2 trial would falsely refuse this valid two-region duct.
    let wrong = text.replace("ambient 293.15 100000 0", "ambient mixture 293.15 100000 mole 1 CO2 1");
    assert!(matches!(ReedPerformance::from_bytes(wrong.as_bytes(), 37), Err(ReedPerformanceError::Model(_))));
}

#[test]
fn regional_gas_syntax_keeps_shared_material_semantics_and_original_error_lines() {
    let text = INPUT.replace(&format!("gas {ENRICHED}"), "gas 303.15 0.5");
    assert_eq!(load(&text, 37).segment_gases()[1], GasState::try_new_moist_air(303.15, 100000.0, 0.5).unwrap());
    let text = INPUT.replace(ENRICHED, "mixture 300 100000 mass 2 N2 0.8 Ar 0.2");
    let expected = resolve_gas_mixture(&[GasSpecies::Nitrogen, GasSpecies::Argon],
        &Composition::new(CompositionBasis::MassFraction, vec![0.8,0.2]).unwrap(), 300.0, 100000.0).unwrap();
    assert_eq!(load(&text, 37).segment_gases()[1], expected);
    let text = INPUT.replace(ENRICHED, "mixture 300 100000 mole 1 Xe 1");
    let Err(error) = ReedPerformance::from_bytes(text.as_bytes(), 37) else { panic!("unsupported species accepted") };
    assert!(matches!(error, ReedPerformanceError::Gas { line: 9, .. }));
    assert!(std::error::Error::source(&error).is_some());
    for record in ["mixture 300 50000 mole 1 N2 1", "mixture 300 100000 mole 2 N2 0.5 nitrogen 0.5",
        "mixture 300 100000 mole 1 N2 1 extra", "293.15", "mixture"] {
        let text = INPUT.replace(ENRICHED, record);
        assert!(ReedPerformance::from_bytes(text.as_bytes(), 37).is_err(), "accepted {record}");
    }
    let text = INPUT.replace("segments 2\ncylinder 0.0015 0.2", "segments 4\ncylinder 0.0015 0.2\nhole 0.0004 0.001 0.0015 1")
        .replace("cylinder 0.0015 0.3 gas", "hole 0.0004 0.001 0.0015 1 gas")
        .replace("\ncompile_limits", "\ncylinder 0.0015 0.3\ncompile_limits");
    assert!(matches!(ReedPerformance::from_bytes(text.as_bytes(), 37), Err(ReedPerformanceError::Model(_))));
}
