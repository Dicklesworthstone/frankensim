//! Source data -> original stationary gas -> original DG/flue pressure.
use fs_couple::dg_port::input::{DgSourceError, DgSourceInput};
use fs_couple::flue_loop::dg::performance::file::FlueDgFilePerformance;
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_exec::CancelGate;
use fs_material::gas::{GasSpec, GasState};
use fs_material::gas_field::AxialGasField;
use fs_material::gas_species::input::{GasRecord, ThermochemicalPhase, ThermochemicalSources};

#[path = "../src/render/schedule/reed/field_input/tests/support.rs"]
mod support;

const GEOMETRY: &str = "dg_segment cylinder 0.001 0.025\ndg_segment cone 0.001 0.002 0.025";
const SOURCE: &str = concat!(
    "dg_bore 2 open 2 2 12 50 2000 0.01 1e-11 256 4 gas_steady 2 1 mole 2 N2 CO2\n",
    "dg_segment cylinder 0.001 0.025\ndg_segment cone 0.001 0.002 0.025\n",
    "gas_inlet 0.000001 330 1 0\n",
    "gas_cell 0.025 0.00005 290 300\ngas_cell 0.05 0.00005 290 300\n",
    "gas_injection 1 0.0000005 295 0 1\n",
    "gas_solve 275 350 0.0000001 64 0.01 2000 100000 0.00000001 0.00000001 0.05",
);

fn ambient() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}
fn read(
    text: &str, sources: Option<&ThermochemicalSources<'_>>, cancel_at: Option<&str>,
) -> Result<(DgSourceInput, Vec<GasRecord>), DgSourceError> {
    let mut lines = text.lines().enumerate();
    let (line, header) = lines.next().unwrap();
    let mut fields = header.split_ascii_whitespace();
    let kind = fields.next().unwrap();
    let mut records = Vec::new();
    let gate = CancelGate::new_clock_free();
    let source = DgSourceInput::read(
        kind, line + 1, &mut fields,
        &mut || {
            let (line, text) = lines.next().ok_or(DgSourceError::Input {
                line: 0, what: "missing test source row",
            })?;
            if cancel_at.is_some_and(|key| text.starts_with(key)) { gate.request(); }
            Ok::<_, DgSourceError>((line + 1, text.split_ascii_whitespace()))
        },
        &ambient(), sources, &mut records, &gate,
    )?;
    assert!(lines.next().is_none(), "the complete source must be consumed");
    Ok((source, records))
}
fn prescribed(field: &AxialGasField) -> String {
    // The material owner canonicalizes these exact source IDs as CO2,N2.
    let mut text = format!(
        "dg_bore 2 open 2 2 12 50 2000 0.01 1e-11 256 4 gas_field {} mole 2 CO2 N2\n{GEOMETRY}",
        field.stations_m().len(),
    );
    for (x, state) in field.stations_m().iter().zip(field.primitive_states()) {
        let fractions = state.composition.fractions();
        text.push_str(&format!("\ngas_station {x} {} {} {}", state.temperature_k, fractions[0], fractions[1]));
    }
    text
}
fn performance(source: &str) -> String {
    let example = include_str!("../examples/flue-dg.performance");
    let header = example.split_once("\ngesture\n").unwrap().0
        .replace("output 2048", "output 512")
        .replace("dg_tube 0.05 0.001 open 1 2 8 50 2000 0.01 1e-11 256", source);
    format!("{header}\ngesture\nfrankensim-gesture-schedule-v1\ncontrol_rate_hz\t1000\ntracks\t1\ntrack\tblow\tblowing-pressure\ninitial\tpressure-pa\t2e1\nevents\t0\n")
}
fn render(text: &str, sources: &ThermochemicalSources<'_>, block: usize) -> Vec<f64> {
    let source = FlueDgFilePerformance::from_bytes_with_gas_sources(
        text.as_bytes(), block, &CancelGate::new_clock_free(), Some(sources),
    ).unwrap();
    let mut renderer = source.into_renderer();
    let mut pressure = vec![0.0; 512];
    for chunk in pressure.chunks_mut(block) { renderer.block(chunk).unwrap(); }
    assert!(pressure.iter().all(|p| p.is_finite()));
    assert_eq!(renderer.samples_rendered(), 512);
    pressure
}

#[test]
fn generated_geometry_sources_and_dg_operator_equal_the_original_prescribed_field() {
    let a = support::packs("N2", 0.028, 3.5, 200.0);
    let b = support::packs("CO2", 0.044, 4.5, 200.0);
    let bindings = [a.binding(), b.binding()];
    let sources = ThermochemicalSources { species: &bindings, water_species: None, phase: ThermochemicalPhase::FrozenGas };
    let (generated, records) = read(SOURCE, Some(&sources), None).unwrap();
    let gas = generated.steady_gas().unwrap();
    let cells = &gas.problem().cells;
    let pi = core::f64::consts::PI;
    let expected_volume = [pi * 0.001_f64.powi(2) * 0.025,
        pi * 0.025 * (0.001_f64.powi(2) + 0.001 * 0.002 + 0.002_f64.powi(2)) / 3.0];
    let expected_wall = [2.0 * pi * 0.001 * 0.025, pi * 0.003 * 0.025_f64.hypot(0.001)];
    for (index, cell) in cells.iter().enumerate() {
        assert!((cell.volume_m3 / expected_volume[index] - 1.0).abs() < 2e-14);
        assert!((cell.wall_area_m2 / expected_wall[index] - 1.0).abs() < 2e-14);
    }
    assert_eq!(gas.solution().mass_flow_kg_s(), &[1e-6, 1.5e-6]);
    assert!(gas.solution().maximum_power_residual_bound_w() <= gas.controls().steady_power_tolerance_w);
    let field = gas.solution().field().clone();
    assert_eq!(records.len(), field.stations_m().len());
    for (record, &x) in records.iter().zip(field.stations_m()) {
        assert_eq!(record, &GasRecord::Thermochemical(field.state_at(x).unwrap()));
    }
    let (direct, direct_records) = read(&prescribed(&field), Some(&sources), None).unwrap();
    assert!(direct.steady_gas().is_none());
    assert_eq!(records, direct_records);
    let gate = CancelGate::new_clock_free();
    let load = generated.prepare(&ambient(), 1.0 / 96000.0, &gate).unwrap();
    let expected = direct.prepare(&ambient(), 1.0 / 96000.0, &gate).unwrap();
    assert_eq!(load.duct().gas_field(), Some(&field));
    assert_eq!(load.duct().operator_identity(), expected.duct().operator_identity());
    assert_eq!(load.impedance_pa_s_m3(), expected.impedance_pa_s_m3());
    assert_ne!(load.duct().inlet_gas(), load.duct().outlet_gas());
}

#[test]
fn generated_gas_changes_actual_flue_pressure_without_changing_callback_clock() {
    let a = support::packs("N2", 0.028, 3.5, 200.0);
    let b = support::packs("CO2", 0.044, 4.5, 200.0);
    let bindings = [a.binding(), b.binding()];
    let sources = ThermochemicalSources { species: &bindings, water_species: None, phase: ThermochemicalPhase::FrozenGas };
    let (input, _) = read(SOURCE, Some(&sources), None).unwrap();
    let field = input.steady_gas().unwrap().solution().field();
    let generated = performance(SOURCE);
    let reference = render(&generated, &sources, 37);
    assert!(reference.iter().any(|p| p.abs() > 1e-6));
    assert_eq!(reference, render(&performance(&prescribed(field)), &sources, 37));
    assert_eq!(reference, render(&generated, &sources, 1));
    assert_eq!(reference, render(&generated, &sources, 512));
    for changed in [
        SOURCE.replace("gas_cell 0.05 0.00005 290", "gas_cell 0.05 0.00005 320"),
        SOURCE.replace("gas_injection 1 0.0000005 295 0 1", "gas_injection 1 0.0000005 295 1 0"),
    ] {
        let (other, _) = read(&changed, Some(&sources), None).unwrap();
        let other = other.steady_gas().unwrap().solution().field();
        assert_eq!(field.state_at(0.0).unwrap(), other.state_at(0.0).unwrap());
        assert_ne!(field.state_at(0.05).unwrap(), other.state_at(0.05).unwrap());
        assert!(reference.iter().zip(render(&performance(&changed), &sources, 37))
            .any(|(a, b)| (a - b).abs() > 1e-10));
    }
}

#[test]
fn generated_dg_gas_refuses_bad_sources_unmodeled_flow_and_late_cancellation() {
    let a = support::packs("N2", 0.028, 3.5, 200.0);
    let b = support::packs("CO2", 0.044, 4.5, 200.0);
    let bindings = [a.binding(), b.binding()];
    let sources = ThermochemicalSources { species: &bindings, water_species: None, phase: ThermochemicalPhase::FrozenGas };
    assert!(read(SOURCE, None, None).is_err());
    for input in [
        SOURCE.replace("gas_steady 2 1", "gas_steady 0 1"),
        SOURCE.replace("gas_cell 0.05", "gas_cell 0.049"),
        SOURCE.replace("gas_cell 0.025", "gas_cell 0.051"),
        SOURCE.replace("gas_injection 1", "gas_injection 2"),
        SOURCE.replace("mole 2 N2 CO2", "mole 2 N2 N2"),
        SOURCE.replace("0.01 2000 100000", "0.01 1 100000"),
        SOURCE.replace("gas_inlet 0.000001", "gas_inlet 0.01"),
        SOURCE.replace("gas_inlet 0.000001 330", "gas_inlet 0.000001 1001"),
        SOURCE.replace("256 4 gas_steady", "256 4 gas_field 3 mole 2 N2 CO2 gas_steady"),
    ] {
        assert!(read(&input, Some(&sources), None).is_err(), "{input}");
    }
    for suffix in ["dg_side_openings 1", "dg_terminal_channels 1 4"] {
        let text = format!("{}\n{suffix}", SOURCE.replace("2 open 2", "2 baffled 2"));
        let error = read(&text, Some(&sources), None).err().expect("unsupported mean flow must refuse");
        assert!(matches!(error, DgSourceError::Input { .. }));
        assert!(error.to_string().contains("gas_steady"));
    }
    let cancelled = read(SOURCE, Some(&sources), Some("gas_solve")).err().unwrap();
    assert!(matches!(cancelled, DgSourceError::SteadyGas {
        source: fs_material::gas_field::steady::SteadyGasError::Thermal(
            fs_thermal::enthalpy::DepositionError::Cancelled), ..
    }));
    let retry = read(SOURCE, Some(&sources), None).unwrap().0;
    assert!(retry.steady_gas().is_some());
}
