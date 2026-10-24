use super::*;

fn chimney_input(gas: &str) -> String {
    baffled_tube_input(RADIATION_EXTERIOR)
        .replace("1e-11 128", "1e-11 256")
        .replace("\ndg_radiation", &format!(
            "\ndg_side_openings 1\ndg_side_chimney 7 1 0.01 0 0.0006 0.012 1 2\ndg_chimney_gas {gas}\ndg_radiation",
        ))
}

#[test]
fn finite_chimney_retains_primitive_material_geometry_and_original_wall_owner() {
    let text = chimney_input("mixture 303.15 101325 mole 2 N2 0.6 CO2 0.4");
    let (source, records) = read_receipted(&text, None).unwrap();
    assert_eq!(records.len(), 2);
    assert_eq!(records[0], source.chimneys[0].gas);
    assert_eq!(records[1], source.radiation.as_ref().unwrap().1);
    assert_ne!(records[0].state().sound_speed, ambient().sound_speed);
    let branches = source.prepare_chimneys(&CancelGate::new_clock_free()).unwrap();
    assert_eq!(branches.len(), 1);
    assert_eq!(branches[0].port, PortId(7));
    let actual = &branches[0].duct;
    assert_eq!(actual.inlet_gas(), records[0].state());
    assert_eq!(actual.outlet_gas(), records[0].state());
    assert_eq!(actual.total_length(), 0.012);
    assert_eq!(actual.total_dofs(), 54);
    assert_eq!(actual.inlet_area_m2(), core::f64::consts::PI * 0.0006 * 0.0006);
    assert!(actual.viscothermal_model().is_some());
    let expected = DgDuct::new(
        &Duct { segments: vec![Segment::Cylinder { radius: 0.0006, length: 0.012 }] },
        records[0].state(), DgMesh { elements: 1, order: 2 },
        DgBoundary::Port(PortId(0)), DgBoundary::Port(PortId(1)),
    ).unwrap().with_viscothermal(8, [50.0, 2000.0], 0.01).unwrap();
    assert_eq!(actual.operator_identity(), expected.operator_identity());
    let other = read(&chimney_input("293.15 101325 0")).unwrap()
        .prepare_chimneys(&CancelGate::new_clock_free()).unwrap();
    assert_ne!(actual.operator_identity(), other[0].duct.operator_identity());
}

#[test]
fn chimney_thermochemical_receipts_use_the_shared_species_sources() {
    use fs_material::gas_species::input::ThermochemicalPhase;
    let a = radiation_gas_support::packs("dga", 0.028, 3.5, 200.0);
    let b = radiation_gas_support::packs("dgb", 0.044, 4.5, 200.0);
    let bindings = [a.binding(), b.binding()];
    let sources = ThermochemicalSources {
        species: &bindings, water_species: None, phase: ThermochemicalPhase::FrozenGas,
    };
    let text = chimney_input("thermochemical 298.15 101325 mole 2 dga 0.3 dgb 0.7");
    let (input, records) = read_receipted(&text, Some(&sources)).unwrap();
    assert!(matches!(&records[0], GasRecord::Thermochemical(_)));
    let branches = input.prepare_chimneys(&CancelGate::new_clock_free()).unwrap();
    assert_eq!(branches[0].duct.inlet_gas(), records[0].state());
    assert_ne!(branches[0].duct.inlet_gas(), &ambient());
}

#[test]
fn mixed_flush_and_finite_openings_keep_declared_channel_order_and_total_memory() {
    let text = chimney_input("293.15 101325 0")
        .replace("baffled 2 2", "baffled 3 2")
        .replace("dg_side_openings 1", "dg_side_openings 2")
        .replace("\ndg_radiation", "\ndg_side_opening 9 2 -0.01 0 0.0005\ndg_radiation");
    let source = read(&text).unwrap();
    assert_eq!(source.side_openings.iter().map(|(port, _)| port.id).collect::<Vec<_>>(),
        vec![PortId(7), PortId(9)]);
    assert_eq!(source.chimneys.len(), 1);
    assert_eq!(source.chimneys[0].port, PortId(7));
    // 3*(2+1)*(2+2*8) main states plus 1*(2+1)*(2+2*8) branch states.
    assert!(read(&text.replace("1e-11 256", "1e-11 216")).is_ok());
    assert!(read(&text.replace("1e-11 256", "1e-11 215")).is_err());
}

#[test]
fn finite_chimneys_refuse_incomplete_physics_and_unresolved_state_budgets() {
    let valid = chimney_input("293.15 101325 0");
    for text in [
        valid.replace("0.012 1 2", "0 1 2"),
        valid.replace("0.012 1 2", "0.012 0 2"),
        valid.replace("0.012 1 2", "0.012 1 0"),
        valid.replace("0.012 1 2", "0.012 1"),
        valid.replace("0.012 1 2", "0.012 1 2 extra"),
        valid.replace("dg_chimney_gas", "gas"),
        valid.replace("\ndg_chimney_gas 293.15 101325 0", ""),
        valid.replace("dg_chimney_gas 293.15 101325", "dg_chimney_gas 293.15 90000"),
        valid.replace("dg_chimney_gas 293.15", "dg_chimney_gas NaN"),
        valid.replace("1e-11 256", "1e-11 161"),
    ] {
        assert!(read_receipted(&text, None).is_err(), "admitted {text}");
    }
    assert!(read(&valid.replace("1e-11 256", "1e-11 162")).is_ok());
    let source = read(&valid).unwrap();
    let gate = CancelGate::new_clock_free();
    gate.request();
    assert!(source.prepare_chimneys(&gate).is_err());
}
