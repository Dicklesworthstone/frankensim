use super::*;
use crate::pcm_wav::observation::PressureRenderer;
use fs_material::gas_species::input::read_gas_state;

const SOURCE: &str = include_str!("../../../../../../examples/plate-valve-regional-gas.performance");
const CO2: &str = "mixture 293.15 101325 mole 3 N2 0.4 O2 0.2 CO2 0.4";
fn with_exterior(text: &str, record: &str) -> String {
    text.replacen("\nnetwork ", &format!("\nexterior {record}\nnetwork "), 1)
}
fn load(text: &str, block: usize) -> PlateValvePerformance {
    let resolved = text.replace("audio 96000 9602", "audio 192000 19204");
    PlateValvePerformance::from_bytes(resolved.as_bytes(), block, &CancelGate::new()).unwrap()
}
fn network(p: &PlateValvePerformance) -> &ApertureNetwork {
    let CoupledAperture::Network(n) = p.renderer().system() else { panic!("network source") }; n
}
fn pressure(p: PlateValvePerformance, block: usize) -> Vec<f64> {
    let mut result = vec![0.0; p.info().samples as usize];
    let mut source = p.into_renderer();
    for chunk in result.chunks_mut(block) { source.block(chunk).unwrap(); }
    assert!(result.iter().all(|v| v.is_finite())); result
}

#[test]
fn exterior_changes_radiation_and_plate_feedback_not_internal_constitutives() {
    let text = with_exterior(SOURCE, CO2);
    let a = load(SOURCE, 37); let b = load(&text, 37);
    let gas = read_gas_state(&mut CO2.split_ascii_whitespace(), None).unwrap();
    assert_eq!(b.exterior_gas(), Some(&gas));
    assert_eq!(a.ambient_gas(), b.ambient_gas());
    assert_eq!(network(&a).section_gases(), network(&b).section_gases());
    assert_eq!(network(&a).spec().sections, network(&b).spec().sections);
    // Includes the cavity, junction and every material-derived loss node.
    for (index, (old, new)) in network(&a).spec().nodes.iter().zip(&network(&b).spec().nodes).enumerate() {
        if index != 2 { assert_eq!(old, new, "interior node {index}"); }
    }
    let expected = BaffledRadiationLoad::new(0.007, gas.density, gas.sound_speed,
        1.0/192000.0, 1000.0, &CancelGate::new()).unwrap();
    assert_eq!(b.radiation_loads()[0], (2, expected));
    assert_ne!(a.radiation_loads()[0].1, expected);
    let (original, changed) = (a.renderer().system().aperture().spec(), b.renderer().system().aperture().spec());
    assert_eq!(original.density_kg_m3, changed.density_kg_m3);
    assert_eq!(original.impedance_pa_s_m3, changed.impedance_pa_s_m3);
    assert_eq!(original.mass_kg, changed.mass_kg);
    assert_eq!(original.stiffness_n_m, changed.stiffness_n_m);
    let (mut a, mut b) = (a.into_renderer(), b.into_renderer());
    let (mut sound, mut motion) = (false, false);
    for _ in 0..2048 {
        let (mut pa, mut pb) = ([0.0], [0.0]);
        a.block(&mut pa).unwrap(); b.block(&mut pb).unwrap();
        sound |= pa[0].to_bits() != pb[0].to_bits();
        motion |= a.system().aperture().state() != b.system().aperture().state();
    }
    assert!(sound && motion, "exterior must affect loaded physics, not only an audio gain");
}

#[test]
fn matched_exterior_and_callback_partitions_preserve_the_original_regional_sound() {
    let matched = with_exterior(SOURCE, "313.15 101325 0.5");
    let reference = pressure(load(SOURCE, 37), 37);
    assert_eq!(reference, pressure(load(&matched, 113), 113));
    let changed = with_exterior(SOURCE, CO2);
    assert_eq!(pressure(load(&changed, 37), 37), pressure(load(&changed, 113), 113));
    let internal = changed.replace("observation network-baffled 2 0 0 0.2 8 32 1000", "observation inlet");
    let p = load(&internal, 37);
    assert!(p.exterior_gas().is_some());
    assert!(p.renderer().baffled_receiver().is_none());
    assert_eq!(p.radiation_loads(), load(&changed, 37).radiation_loads());
}

#[test]
fn every_radiating_endpoint_uses_the_common_exterior_without_averaging_gases() {
    let text = with_exterior(SOURCE, CO2)
        .replace("duct_node cavity 0.00001 0.0015 0.02 200000", "duct_node baffled-low-ka 1000");
    let p = load(&text, 37); let gas = p.exterior_gas().unwrap();
    assert_eq!(p.radiation_loads().len(), 2);
    for &(node, load) in p.radiation_loads() {
        let section = network(&p).spec().sections.iter().find(|s| s.nodes.contains(&node)).unwrap();
        assert_eq!(load, BaffledRadiationLoad::new(section.radius_m, gas.density,
            gas.sound_speed, 1.0/192000.0, 1000.0, &CancelGate::new()).unwrap());
    }
    assert_ne!(network(&p).section_gases().unwrap()[0], *gas);
}

#[test]
fn radiation_admission_uses_the_actual_exterior_not_a_false_local_gas_trial() {
    let gas = |name| read_gas_state(&mut format!("mixture 293.15 101325 mole 1 {name} 1")
        .split_ascii_whitespace(), None).unwrap();
    let (slow, fast) = (gas("CO2"), gas("N2"));
    let make = |air: GasState| DuctInput::Tube { spec: UniformTubeSpec {
        length_m: 0.25, radius_m: 0.02, terminal_reflection: 0.0,
        max_length_error_m: 0.003, max_wave_memory_bytes: 1048576, sound_speed_m_s: air.sound_speed,
    }, radiation_band: Some(1100.0), loss: None };
    let gate = CancelGate::new(); let dt = 1.0/96000.0;
    assert!(make(slow).prepare(&slow, dt, ApertureObservation::Inlet, &gate).is_err());
    let actual = make(slow).prepare_in_exterior(&slow, dt, ApertureObservation::Inlet, Some(&fast), &gate).unwrap();
    assert_eq!(actual.gases.as_ref().unwrap(), &[slow]);
    assert_eq!(actual.radiation[0].1, BaffledRadiationLoad::new(0.02, fast.density,
        fast.sound_speed, dt, 1100.0, &gate).unwrap());
    assert!(make(fast).prepare_in_exterior(&fast, dt, ApertureObservation::Inlet, Some(&slow), &gate).is_err());
}

#[test]
fn exterior_records_retain_material_errors_and_refuse_ambiguous_boundaries() {
    for record in ["mixture 293.15 101325 mole 1 He 1", "mixture 300 101325 mole 1 H2O 1"] {
        let result = PlateValvePerformance::from_bytes(with_exterior(SOURCE, record).as_bytes(), 37, &CancelGate::new());
        assert!(matches!(result, Err(PlateValveInputError::Gas { line: 4, .. })));
    }
    for text in [with_exterior(SOURCE, "293.15 50000 0"),
        with_exterior(SOURCE, CO2).replace("duct_node baffled-low-ka 1000", "duct_node reflection -0.9"),
        with_exterior(SOURCE, CO2).replace("\nnetwork ", "\nexterior 293.15 101325 0\nnetwork "),
        with_exterior(SOURCE, "mixture 300 101325 mole 2 N2 0.5 N2 0.5"),
        with_exterior(SOURCE, "mixture 300 101325 mole 1 N2 1 extra")]
    {
        assert!(PlateValvePerformance::from_bytes(text.as_bytes(), 37, &CancelGate::new()).is_err());
    }
    let gate = CancelGate::new(); gate.request();
    assert!(PlateValvePerformance::from_bytes(with_exterior(SOURCE, CO2).as_bytes(), 37, &gate).is_err());
}
