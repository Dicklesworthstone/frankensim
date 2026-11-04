use super::*;

fn branch(gas: &GasRecord, radius_m: f64) -> DgBranch {
    DgBranch {
        port: PORTS[0],
        duct: DgDuct::new(
            &Duct { segments: vec![Segment::Cylinder { radius: radius_m, length: 0.012 }] },
            gas.state(), DgMesh { elements: 1, order: 2 },
            DgBoundary::Port(PortId(20)), DgBoundary::Port(PortId(21)),
        ).unwrap(),
    }
}

#[test]
fn computed_chimney_geometry_is_the_actual_outer_area_and_material_source() {
    let gas = exterior("293.15 101325 0");
    let other = exterior("mixture 303.15 101325 mole 2 N2 0.6 CO2 0.4");
    let main = duct(&gas);
    let branches = vec![branch(&other, 0.0015)];
    admit_branches(&main, &branches, 1.0 / 8000.0, 64).unwrap();
    admit_openings(&main, &branches, &gas, &openings()).unwrap();
    let mut wrong = openings();
    wrong[0].1.radius_m = 0.001; // A branch outlet cannot borrow a different circle.
    assert!(admit_openings(&main, &branches, &gas, &wrong).is_err());
    assert!(admit_branches(&main, &branches, 1.0 / 8000.0, 17).is_err());
    let wrong_pressure = vec![branch(&exterior("293.15 90000 0"), 0.0015)];
    assert!(admit_branches(&main, &wrong_pressure, 1.0 / 8000.0, 64).is_err());
    let mut unknown = branch(&gas, 0.0015);
    unknown.port = PortId(91);
    assert!(admit_branches(&main, &[unknown], 1.0 / 8000.0, 64).is_err());
    assert!(admit_branches(
        &main, &[branch(&gas, 0.0015), branch(&gas, 0.0015)], 1.0 / 8000.0, 64,
    ).is_err());
}

#[test]
fn finite_branch_outlet_flows_drive_the_same_mutual_source_and_microphone() {
    let gas = exterior("293.15 101325 0");
    let chimney_gas = exterior("mixture 303.15 101325 mole 2 N2 0.6 CO2 0.4");
    let controls = ComputedBaffledInput::parse(
        "20 65 4 8 12 20 0.01 67108864 2 100 1e-8 64 1048576 257",
    ).unwrap();
    let openings = openings();
    let mut load = DgPortLoad::new_with_computed_baffled_branches(
        duct(&gas), 8000, 1e-11, 64, &gas, controls, &openings,
        vec![branch(&chimney_gas, openings[0].1.radius_m)],
        &CancelGate::new_clock_free(),
    ).unwrap();
    assert_eq!(load.terminal_ports(), &PORTS);
    assert_eq!(load.branches().len(), 1);
    assert_eq!(load.branches()[0].port(), PORTS[0]);
    assert_eq!(load.branches()[0].duct().outlet_gas(), chimney_gas.state());
    assert_eq!(load.terminal_gas(PORTS[0]).unwrap(), *chimney_gas.state());
    assert_eq!(load.terminal_area_m2(PORTS[0]).unwrap(),
        Some(core::f64::consts::PI * 0.0015 * 0.0015));
    assert_eq!(load.radiation_record().unwrap().gas_record(), &gas);
    let config = CircularOutletReceiver {
        position_m: [0.1, 0.05, 1.0], radial_rings: 8, angular_points: 32,
        maximum_frequency_hz: 1000.0,
    };
    let mut combined = load.baffled_receiver(8000, config).unwrap();
    let mut individual: Vec<_> = openings.iter().map(|(_, aperture)|
        BaffledPressure::circular_openings_multirate_differentiating(
            &[*aperture], 8000, 8000, config, gas.state().into(),
        ).unwrap()
    ).collect();
    let mut energy = [0.0; 2];
    let mut maximum_chimney_storage: f64 = 0.0;
    for sample in 0..512 {
        load.preview_step(if sample < 8 { 1e-9 } else { 0.0 }).unwrap().accept();
        let frames = load.last_terminal_frames();
        assert_eq!([frames[0].port, frames[1].port], PORTS);
        let flows = [frames[0].outward_flow_m3_s, frames[1].outward_flow_m3_s];
        let pressure = combined.step_trace(&flows).unwrap();
        let a = individual[0].step_trace(&flows[..1]).unwrap();
        let b = individual[1].step_trace(&flows[1..]).unwrap();
        assert!((pressure - a - b).abs() <= 2e-12 * (a.abs() + b.abs()).max(1e-300));
        energy[0] += a * a;
        energy[1] += b * b;
        maximum_chimney_storage = maximum_chimney_storage.max(
            load.branches()[0].state().stored_energy_j(),
        );
        assert_eq!(load.branches()[0].state().step_count(), load.state().step_count());
    }
    assert!(maximum_chimney_storage > 0.0);
    assert!(energy.iter().all(|value| *value > 0.0));
}
