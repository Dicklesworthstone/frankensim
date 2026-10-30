use super::*;
use fs_material::gas::GasSpec;

fn ambient() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}

fn read(text: &str) -> Result<DgSourceInput, DgSourceError> {
    let mut lines = text.lines().enumerate();
    let (line, header) = lines.next().unwrap();
    let mut fields = header.split_ascii_whitespace();
    let kind = fields.next().unwrap();
    let source = DgSourceInput::read(
        kind, line + 1, &mut fields,
        &mut || {
            let (line, text) = lines.next().ok_or_else(|| bad(0, "missing source row"))?;
            Ok::<_, DgSourceError>((line + 1, text.split_ascii_whitespace()))
        },
        &ambient(), None, &mut Vec::new(), &CancelGate::new_clock_free(),
    )?;
    assert!(lines.next().is_none());
    Ok(source)
}

const TUBE: &str = "dg_tube 0.05 0.001 open 2 2 8 50 2000 0.01 1e-11 128";
// Twelve wall poles: at 50-2000 Hz and 1 %, eight resolve the 1 mm end but
// refuse from about 1.6 mm, inside the 2 mm cone (.362 item 13).
const BORE: &str = "dg_bore 2 open 2 2 12 50 2000 0.01 1e-11 224 4\ndg_segment cylinder 0.001 0.025\ndg_segment cone 0.001 0.002 0.025";

#[test]
fn g0_uniform_source_retains_original_geometry_gas_and_exact_mass_wall_owner() {
    let load = read(TUBE).unwrap().prepare(
        &ambient(), 1.0 / 96000.0, &CancelGate::new_clock_free(),
    ).unwrap();
    assert_eq!(load.duct().inlet_gas(), &ambient());
    assert_eq!(load.duct().outlet_gas(), &ambient());
    assert_eq!(load.duct().total_length(), 0.05);
    assert_eq!(load.duct().total_dofs(), 108);
    assert!(load.duct().viscothermal_model().is_some());
    assert!(load.duct().spatial_viscothermal_admission().is_none());
    assert_eq!(load.state().step_count(), 0);
}

#[test]
fn g3_authored_bore_uses_one_existing_spatial_wall_operator() {
    let load = read(BORE).unwrap().prepare(
        &ambient(), 1.0 / 96000.0, &CancelGate::new_clock_free(),
    ).unwrap();
    let duct = load.duct();
    let wall = duct.spatial_viscothermal_admission().unwrap();
    assert_eq!(wall.resolution.quadrature_points, 4);
    assert_eq!(wall.resolution.poles_per_field, 12);
    assert_eq!(wall.memory_dofs, 192);
    assert_eq!(duct.total_dofs(), 204);
    assert_eq!(duct.total_length(), 0.05);
    assert_eq!(duct.inlet_area_m2(), core::f64::consts::PI * 0.001 * 0.001);
    assert_eq!(duct.outlet_area_m2(), core::f64::consts::PI * 0.002 * 0.002);
    assert_eq!(duct.viscothermal_band_rad_s(), Some(wall.resolution.band_rad_s));
    assert!(duct.viscothermal_model().is_none());
    assert!(duct.has_viscothermal_memory());
    let uniform = read(TUBE).unwrap().prepare(
        &ambient(), 1.0 / 96000.0, &CancelGate::new_clock_free(),
    ).unwrap();
    assert_ne!(duct.operator_identity(), uniform.duct().operator_identity());
    assert_ne!(load.impedance_pa_s_m3(), uniform.impedance_pa_s_m3());
}

#[test]
fn g0_source_refuses_incomplete_geometry_unbounded_memory_and_unresolved_band() {
    for text in [
        BORE.replace("2 open 2", "2 open 1"),
        BORE.replace("224 4", "203 4"),
        BORE.replace("224 4", "224 2"),
        BORE.replace("224 4", "2048 65"),
        BORE.replace("dg_segment cone", "cone"),
        BORE.replace("0.001 0.002", "-0.001 0.002"),
        BORE.replace("50 2000", "2000 50"),
        BORE.replace("224 4", "224 4 gas_field 2 mole 1 missing"),
        TUBE.replace("128", "107"),
    ] {
        assert!(read(&text).is_err());
    }
    let source = read(&BORE.replace("50 2000", "50 50000")).unwrap();
    assert!(source.prepare(&ambient(), 1.0 / 96000.0, &CancelGate::new_clock_free()).is_err());
    // Eight poles cannot resolve the cone: the refusal names the wall law's
    // own reason and where along the bore it refused.
    let source = read(&BORE.replace(" 12 50 2000", " 8 50 2000")).unwrap();
    match source.prepare(&ambient(), 1.0 / 96000.0, &CancelGate::new_clock_free()) {
        Err(DgSourceError::Physics(AcousticRealizeError::Duct(
            fs_duct::DuctError::WallDiffusionUnresolved { radius_m, x_m, .. },
        ))) => assert!(radius_m > 0.0015 && radius_m <= 0.002 && x_m > 0.025, "{radius_m} {x_m}"),
        Err(other) => panic!("eight poles must refuse at the cone: {other:?}"),
        Ok(_) => panic!("eight poles must refuse at the cone"),
    }
    let source = read(BORE).unwrap();
    let gate = CancelGate::new_clock_free();
    gate.request();
    assert!(source.prepare(&ambient(), 1.0 / 96000.0, &gate).is_err());
}

#[path = "../../render/schedule/reed/field_input/tests/support.rs"]
mod radiation_gas_support;

// Collect the exact source receipts without preparing any DG propagator or
// fitting a radiation operator. Full-source tests reject an unconsumed row.
fn read_receipted(
    text: &str,
    sources: Option<&ThermochemicalSources<'_>>,
) -> Result<(DgSourceInput, Vec<GasRecord>), DgSourceError> {
    let mut lines = text.lines().enumerate();
    let (line, header) = lines.next().ok_or_else(|| bad(1, "missing source header"))?;
    let mut fields = header.split_ascii_whitespace();
    let kind = fields.next().ok_or_else(|| bad(line + 1, "missing source kind"))?;
    let mut records = Vec::new();
    let source = DgSourceInput::read(
        kind, line + 1, &mut fields,
        &mut || {
            let (line, text) = lines.next().ok_or_else(|| bad(0, "missing source row"))?;
            Ok::<_, DgSourceError>((line + 1, text.split_ascii_whitespace()))
        },
        &ambient(), sources, &mut records, &CancelGate::new_clock_free(),
    )?;
    if let Some((line, _)) = lines.next() {
        return Err(bad(line + 1, "unexpected source test row"));
    }
    Ok((source, records))
}

const COMPUTED_CONTROLS: &str =
    "20 65 4 8 12 20 0.1 67108864 2 100 1e-8 64 1048576 257";
const RADIATION_EXTERIOR: &str = "mixture 293.15 101325 mole 2 N2 0.8 CO2 0.2";

#[test]
fn authored_flush_side_openings_keep_original_faces_ids_and_independent_areas() {
    let text = baffled_tube_input(RADIATION_EXTERIOR).replace("\ndg_radiation", concat!(
        "\ndg_side_openings 1\ndg_side_opening 7 1 0.01 -0.002 0.0004",
        "\ndg_radiation",
    ));
    let (source, records) = read_receipted(&text, None).unwrap();
    let (port, aperture) = source.side_openings[0];
    assert_eq!(port, DgSidePort { id: PortId(7), downstream_element: 1 });
    assert_eq!(aperture.center_m, [0.01, -0.002]);
    assert_eq!(aperture.radius_m, 0.0004);
    assert_eq!(source.mesh, DgMesh { elements: 2, order: 2 });
    assert_eq!(records.len(), 1);
    assert_eq!(source.radiation.as_ref().unwrap().1, records[0]);
    assert!(read(&baffled_tube_input(RADIATION_EXTERIOR)).unwrap().side_openings.is_empty());
    for invalid in [
        text.replace("dg_side_openings 1", "dg_side_openings 0"),
        text.replace("dg_side_openings 1", "dg_side_openings 8"),
        text.replace("dg_side_opening 7 1", "dg_side_opening 0 1"),
        text.replace("dg_side_opening 7 1", "dg_side_opening 1 1"),
        text.replace("dg_side_opening 7 1", "dg_side_opening 7 0"),
        text.replace("dg_side_opening 7 1", "dg_side_opening 7 2"),
        text.replace("0.0004\n", "0\n"),
        text.replace("0.01 -0.002", "NaN -0.002"),
        text.replace("dg_side_openings 1", "dg_side_openings 2")
            .replace("\ndg_radiation", "\ndg_side_opening 8 1 0.02 0 0.0004\ndg_radiation"),
        text.replace("dg_side_opening 7 1", "dg_side_opening 7 1 extra"),
    ] {
        assert!(read_receipted(&invalid, None).is_err(), "admitted {invalid}");
    }
}

fn baffled_tube_input(exterior: &str) -> String {
    format!("{}\ndg_radiation computed-baffled {COMPUTED_CONTROLS}\ndg_exterior {exterior}",
        TUBE.replace(" open ", " baffled "))
}

#[test]
fn g0_baffled_input_retains_exact_geometry_controls_and_original_gas_records() {
    for exterior in ["293.15 101325 0", RADIATION_EXTERIOR] {
        let mut fields = exterior.split_ascii_whitespace();
        let expected = fs_material::gas_species::input::read_gas_record(
            &mut fields, None, None,
        ).unwrap();
        assert!(fields.next().is_none());
        let (source, records) = read_receipted(&baffled_tube_input(exterior), None).unwrap();
        assert!(matches!(source.outlet, DgBoundary::Port(PortId(1))));
        assert!(source.field.is_none());
        assert!(matches!(source.geometry.segments.as_slice(),
            [Segment::Cylinder { radius, length }] if *radius == 0.001 && *length == 0.05));
        let (input, actual) = source.radiation.as_ref().unwrap();
        assert_eq!(actual, &expected);
        assert_eq!(records, vec![expected]);
        assert_eq!(actual.state().pressure.to_bits(), ambient().pressure.to_bits());
        if exterior == RADIATION_EXTERIOR {
            assert_ne!(actual.state().sound_speed.to_bits(), ambient().sound_speed.to_bits());
        }
        assert_eq!(input.minimum_frequency_hz(), 20.0);
        let options = input.options();
        assert_eq!(options.frequency_samples, 65);
        assert_eq!(options.orders, [4, 8, 12]);
        assert_eq!(options.fit_iterations, 20);
        assert_eq!(options.maximum_relative_error, 0.1);
        assert_eq!(options.max_fit_memory_bytes, 67_108_864);
        assert_eq!(options.causality.subtractions, 2);
        assert_eq!(options.causality.tail_bound_ratio, 100.0);
        assert_eq!(options.causality.relative_data_tolerance, 1e-8);
        let limits = input.limits();
        assert_eq!(limits.max_states, 64);
        assert_eq!(limits.max_memory_bytes, 1_048_576);
        assert_eq!(limits.sampled_passivity_points, 257);
    }
}

#[test]
fn g0_baffled_records_follow_geometry_and_field_and_keep_thermochemical_receipts() {
    use fs_material::gas_species::input::{ThermochemicalPhase, read_gas_record};
    let a = radiation_gas_support::packs("dga", 0.028, 3.5, 200.0);
    let b = radiation_gas_support::packs("dgb", 0.044, 4.5, 200.0);
    let bindings = [a.binding(), b.binding()];
    let sources = ThermochemicalSources {
        species: &bindings, water_species: None, phase: ThermochemicalPhase::FrozenGas,
    };
    let (header, geometry) = BORE.split_once('\n').unwrap();
    let header = format!("{} gas_field 2 mole 2 dgb dga",
        header.replace(" open ", " baffled "));
    let stations = "gas_station 0 293.15 0.25 0.75\ngas_station 0.05 303.15 0.5 0.5";
    let radiation = format!("dg_radiation computed-baffled {COMPUTED_CONTROLS}");
    let exterior = "thermochemical 298.15 101325 mole 2 dga 0.2 dgb 0.8";
    let text = format!("{header}\n{geometry}\n{stations}\n{radiation}\ndg_exterior {exterior}");
    let (source, records) = read_receipted(&text, Some(&sources)).unwrap();
    assert!(matches!(source.geometry.segments.as_slice(),
        [Segment::Cylinder { radius, length: first },
         Segment::Cone { inlet_radius, outlet_radius, length: last }]
            if *radius == 0.001 && *first == 0.025 && *inlet_radius == 0.001
                && *outlet_radius == 0.002 && *last == 0.025));
    assert!(matches!(source.outlet, DgBoundary::Port(PortId(1))));
    let field = source.field.as_ref().unwrap();
    assert_eq!(field.stations_m(), &[0.0, 0.05]);
    assert_eq!(records.len(), 3);
    for (index, position) in [0.0, 0.05].into_iter().enumerate() {
        assert_eq!(records[index], GasRecord::Thermochemical(field.state_at(position).unwrap()));
    }
    let mut fields = exterior.split_ascii_whitespace();
    let expected = read_gas_record(&mut fields, None, Some(&sources)).unwrap();
    assert!(fields.next().is_none());
    assert!(matches!(&expected, GasRecord::Thermochemical(_)));
    assert_eq!(&records[2], &expected);
    let (_, actual) = source.radiation.as_ref().unwrap();
    assert_eq!(actual, &expected);
    assert_ne!(actual.state().sound_speed.to_bits(),
        field.state_at(0.05).unwrap().state().sound_speed.to_bits());
    for wrong_order in [
        format!("{header}\n{radiation}\n{geometry}\n{stations}\ndg_exterior {exterior}"),
        format!("{header}\n{geometry}\n{radiation}\n{stations}\ndg_exterior {exterior}"),
        format!("{header}\n{geometry}\n{stations}\ndg_exterior {exterior}\n{radiation}"),
    ] {
        assert!(read_receipted(&wrong_order, Some(&sources)).is_err(), "{wrong_order}");
    }
}

#[test]
fn g0_baffled_input_refuses_missing_malformed_extra_fields_and_static_pressure_jumps() {
    let valid = baffled_tube_input(RADIATION_EXTERIOR);
    for text in [
        valid.split_once("\ndg_radiation").unwrap().0.to_owned(),
        valid.split_once("\ndg_exterior").unwrap().0.to_owned(),
        valid.replace("dg_radiation", "radiation"),
        valid.replace("computed-baffled", "artifact"),
        valid.replace("computed-baffled 20", "computed-baffled NaN"),
        valid.replace(" 257\n", "\n"),
        valid.replace(" 257\n", " 257 999\n"),
        valid.replace(" 4 8 12 ", " 4 4 12 "),
        valid.replace("dg_exterior", "exterior"),
        valid.replace(RADIATION_EXTERIOR, "293.15 101325"),
        valid.replace("mole 2", "volume 2"),
        valid.replace("CO2 0.2", "missing 0.2"),
        valid.replace(RADIATION_EXTERIOR, "thermochemical 293.15 101325 mole 1 missing 1"),
        format!("{valid} extra"),
        valid.replace("101325", "90000"),
    ] {
        assert!(read_receipted(&text, None).is_err(), "admitted {text}");
    }
    let error = read_receipted(
        &valid.replace("computed-baffled 20", "computed-baffled NaN"), None,
    ).err().unwrap();
    assert!(matches!(error, DgSourceError::Radiation { line: 2, .. }));
    let error = read_receipted(&valid.replace("101325", "90000"), None).err().unwrap();
    assert!(matches!(error, DgSourceError::Input { line: 3, .. }));
}

#[path = "tests/chimneys.rs"]
mod chimneys;

mod terminal_channels;

fn branched_steady_input() -> String {
    format!(concat!(
        "dg_bore 1 baffled 2 1 12 50 2000 0.01 1e-11 228 3\n",
        "dg_segment cylinder 0.001 0.05\n",
        "dg_side_openings 1\n",
        "dg_side_chimney 7 1 0.01 0 0.0004 0.01 1 1\n",
        "dg_chimney_gas 293.15 101325 0\n",
        "dg_gas_network mole 2 dga dgb\n",
        "gas_inlet 0.000001 330 0.75 0.25\n",
        "gas_bore 2 0.0000005\n",
        "gas_cell 0.025 0.0005 290 300\n",
        "gas_cell 0.05 0.0005 290 300\n",
        "gas_side 7 0.025 0.0000005 1\n",
        "gas_cell 0.01 0.00025 290 300\n",
        "gas_solve 250 360 0.0000001 64 0.01 5000 200000 0.00000001 0.00000001 0.05\n",
        "dg_radiation computed-baffled {}\n",
        "dg_exterior 293.15 101325 0"
    ), COMPUTED_CONTROLS)
}

#[test]
fn g3_branched_steady_gas_changes_bore_cooling_and_original_chimney_operator() {
    use fs_material::gas_species::input::ThermochemicalPhase;
    let a = radiation_gas_support::packs("dga", 0.028, 3.5, 200.0);
    let b = radiation_gas_support::packs("dgb", 0.044, 4.5, 200.0);
    let bindings = [a.binding(), b.binding()];
    let sources = ThermochemicalSources {
        species: &bindings, water_species: None, phase: ThermochemicalPhase::FrozenGas,
    };
    let text = branched_steady_input();
    let (source, records) = read_receipted(&text, Some(&sources)).unwrap();
    let network = source.steady_gas_network().unwrap();
    let solution = network.solution();
    assert_eq!(solution.mass_flow_kg_s(), &[1e-6, 5e-7, 5e-7]);
    assert!(solution.species().maximum_residual_kg_s() < 1e-20);
    assert!(solution.maximum_power_residual_bound_w() <= 1e-8);
    assert!(solution.maximum_energy_residual_j() <= 1e-8);
    assert_eq!(network.problem().cells.len(), 3);
    assert_eq!(network.problem().flows.len(), 2);
    assert_eq!(network.junctions(), &[(7, 0.025)]);
    let main = network.main_field();
    let branch = network.branches()[0].field();
    assert_eq!(main.stations_m(), &[0.0, 0.025, 0.05]);
    assert_eq!(branch.stations_m(), &[0.0, 0.01]);
    assert_eq!(branch.state_at(0.0).unwrap(), main.state_at(0.025).unwrap());
    assert!(branch.state_at(0.01).unwrap().state().temperature
        < branch.state_at(0.0).unwrap().state().temperature);
    assert_eq!(records.len(), 7); // Old homogeneous declaration, five field knots, exterior.
    let branches = source.prepare_chimneys(&CancelGate::new_clock_free()).unwrap();
    let duct = &branches[0].duct;
    assert_eq!(duct.total_dofs(), 76);
    assert_eq!(duct.inlet_gas(), branch.state_at(0.0).unwrap().state());
    assert_eq!(duct.outlet_gas(), branch.state_at(0.01).unwrap().state());
    assert!(duct.spatial_viscothermal_admission().is_some());
    assert!(duct.viscothermal_model().is_none());
    let main_duct = DgDuct::new_with_gas_field(&source.geometry, main, source.mesh,
        DgBoundary::Port(PortId(0)), source.outlet).unwrap()
        .with_side_ports(&[source.side_openings[0].0]).unwrap();
    assert_eq!(main_duct.side_port_source(PortId(7)).unwrap().station_m, 0.025);
    assert_eq!(main_duct.side_port_source(PortId(7)).unwrap().gas, *duct.inlet_gas());
    // Same inlet and walls: sending less gas out the side leaves more warm
    // breath downstream in the main bore, while the side chimney cools more.
    let other_text = text.replace("gas_bore 2 0.0000005", "gas_bore 2 0.00000075")
        .replace("gas_side 7 0.025 0.0000005", "gas_side 7 0.025 0.00000025");
    let (other, _) = read_receipted(&other_text, Some(&sources)).unwrap();
    let other_network = other.steady_gas_network().unwrap();
    assert!(other_network.main_field().state_at(0.05).unwrap().state().temperature
        > main.state_at(0.05).unwrap().state().temperature + 0.1);
    assert!(other_network.branches()[0].field().state_at(0.01).unwrap().state().temperature
        < branch.state_at(0.01).unwrap().state().temperature - 0.1);
    let other_branches = other.prepare_chimneys(&CancelGate::new_clock_free()).unwrap();
    assert_ne!(other_branches[0].duct.operator_identity(), duct.operator_identity());
}

#[test]
fn g0_branched_steady_gas_refuses_unbalanced_unlocated_or_unrepresented_paths() {
    use fs_material::gas_species::input::ThermochemicalPhase;
    let a = radiation_gas_support::packs("dga", 0.028, 3.5, 200.0);
    let b = radiation_gas_support::packs("dgb", 0.044, 4.5, 200.0);
    let bindings = [a.binding(), b.binding()];
    let sources = ThermochemicalSources {
        species: &bindings, water_species: None, phase: ThermochemicalPhase::FrozenGas,
    };
    let valid = branched_steady_input();
    for invalid in [
        valid.replace("gas_inlet 0.000001", "gas_inlet 0.000002"),
        valid.replace("gas_bore 2 0.0000005", "gas_bore 2 0"),
        valid.replace("gas_side 7", "gas_side 9"),
        valid.replace("gas_side 7 0.025", "gas_side 7 0.02"),
        valid.replace("gas_side 7 0.025", "gas_side 7 0.05"),
        valid.replace("gas_side 7 0.025 0.0000005 1", "gas_side 7 0.025 0.0000005 0"),
        valid.replace("gas_bore 2", "gas_bore 65"),
        valid.replace("1e-11 228", "1e-11 227"),
        valid.replace("\ndg_gas_network", "\ndg_terminal_channels 1 4\ndg_terminal_slit 7 0.001 0.0004 0.001 4 0.1\ndg_gas_network"),
        valid.replace("\ndg_gas_network", "\ndg_terminal_channels 1 0\ndg_terminal_closed 1\ndg_gas_network"),
    ] {
        assert!(read_receipted(&invalid, Some(&sources)).is_err(), "admitted {invalid}");
    }
    // A thermal boundary does not silently relocate an existing acoustic port.
    let displaced = valid.replace("gas_cell 0.025", "gas_cell 0.02")
        .replace("gas_side 7 0.025", "gas_side 7 0.02")
        .replace("dg_side_chimney 7 1", "dg_side_chimney 7 2")
        .replace("baffled 2 1", "baffled 3 1").replace("1e-11 228", "1e-11 256")
        .replace(" 12 50", " 8 50");
    let (source, _) = read_receipted(&displaced, Some(&sources)).unwrap();
    assert!(matches!(source.prepare(&ambient(), 1.0 / 96000.0,
        &CancelGate::new_clock_free()), Err(DgSourceError::Input { what, .. })
        if what.contains("thermal side junction")));
}

