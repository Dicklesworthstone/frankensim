use super::*;
use fs_exec::CancelGate;
use fs_material::gas_species::input::read_gas_state;

const SOURCE: &str = include_str!("../../../examples/plate-valve-regional-gas.performance");

fn source() -> AperturePerformance {
    file::PlateValvePerformance::from_bytes(SOURCE.as_bytes(), 113, &CancelGate::new())
        .unwrap().into_renderer()
}
fn gas(pressure: f64) -> GasState {
    read_gas_state(&mut format!("mixture 293.15 {pressure} mole 3 N2 0.4 O2 0.2 CO2 0.4")
        .split_ascii_whitespace(), None).unwrap()
}
fn local(p: &AperturePerformance) -> GasState {
    let (CoupledAperture::Network(n), ApertureObservation::NetworkBaffled {node, ..})
        = (&p.system, p.observation) else { panic!("regional fixture has an outlet") };
    let (index, _) = n.spec().sections.iter().enumerate().find(|(_, s)| s.nodes.contains(&node)).unwrap();
    n.section_gases().unwrap()[index]
}
fn in_gas(exterior: GasState) -> AperturePerformance {
    let old = source();
    AperturePerformance::new_with_exterior_gas(old.system, old.observation,
        old.schedule, old.config, exterior).unwrap()
}

#[test]
fn independent_receiver_matches_existing_rayleigh_without_changing_mechanics() {
    let mut reference = source();
    let exterior = gas(local(&reference).pressure);
    let mut observed = in_gas(exterior);
    let (CoupledAperture::Network(n), ApertureObservation::NetworkBaffled {node, receiver})
        = (&reference.system, reference.observation) else { unreachable!() };
    let section = n.spec().sections.iter().find(|s| s.nodes.contains(&node)).unwrap();
    let medium = RayleighMedium { density: exterior.density, sound_speed: exterior.sound_speed };
    let mut mic = BaffledPressure::circular_outlet(section.radius_m,
        reference.config.sample_rate_hz, receiver, medium).unwrap();
    assert_eq!(observed.baffled_medium(), Some(medium));
    assert_eq!(observed.exterior_gas(), Some(&exterior));
    let mut changed = false;
    for _ in 0..2048 {
        let (mut a, mut b) = ([0.0], [0.0]);
        reference.block(&mut a).unwrap();
        observed.block(&mut b).unwrap();
        let (CoupledAperture::Network(n), CoupledAperture::Network(m))
            = (&reference.system, &observed.system) else { unreachable!() };
        let flow = n.node_frame(node).unwrap().net_flow_into_node_m3_s;
        assert_eq!(flow.to_bits(), m.node_frame(node).unwrap().net_flow_into_node_m3_s.to_bits());
        assert_eq!(b[0].to_bits(), mic.step(&[flow]).unwrap().to_bits());
        changed |= a[0].to_bits() != b[0].to_bits();
    }
    assert!(changed, "exterior must change receiver pressure, not internal flow");
    assert_eq!(reference.applied_controls(), observed.applied_controls());
}

#[test]
fn explicit_matched_exterior_retains_legacy_samples_and_partitioning() {
    let mut old = source();
    let exterior = local(&old);
    let mut matched = in_gas(exterior);
    let mut partitioned = in_gas(exterior);
    let mut expected = vec![0.0; 2048];
    let mut actual = expected.clone();
    let mut split = expected.clone();
    for block in expected.chunks_mut(113) { old.block(block).unwrap(); }
    for block in actual.chunks_mut(37) { matched.block(block).unwrap(); }
    for block in split.chunks_mut(1) { partitioned.block(block).unwrap(); }
    assert_eq!(expected, actual);
    assert_eq!(actual, split);
}

#[test]
fn exterior_admission_refuses_unknown_static_pressure_and_wrong_observation() {
    let pressure = local(&source()).pressure;
    for exterior in [gas(pressure / 2.0), GasState { density: f64::NAN, ..gas(pressure) }] {
        let p = source();
        assert!(AperturePerformance::new_with_exterior_gas(p.system, p.observation,
            p.schedule, p.config, exterior).is_err());
    }
    let p = source();
    assert!(AperturePerformance::new_with_exterior_gas(p.system, ApertureObservation::Inlet,
        p.schedule, p.config, gas(pressure)).is_err());
    let text = SOURCE.replace(" gas 313.15 0.5", "").replace(" gas 283.15 0", "");
    let p = file::PlateValvePerformance::from_bytes(text.as_bytes(), 37, &CancelGate::new())
        .unwrap().into_renderer();
    assert!(AperturePerformance::new_with_exterior_gas(p.system, p.observation,
        p.schedule, p.config, gas(pressure)).is_err());
}
