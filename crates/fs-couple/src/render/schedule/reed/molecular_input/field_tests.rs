//! Continuous source fields and molecular laws through the actual performance.
use super::super::{GasFieldResolution, ObservationPolicy, ReedPerformance, field_input};
use super::*;
use fs_material::gas_species::input::{ThermochemicalPhase, ThermochemicalSources};

#[path = "../field_input/tests/support.rs"]
mod support;
use support::packs;

const FIELD: &str = "gas_field continuous 8 3 mole 2 A B\ngas_station 0 300 0.8 0.2\ngas_station 0.075 330 0.7 0.3\ngas_station 0.15 310 0.6 0.4\n";

fn declaration(a: &str, b: &str) -> String {
    format!(
        "molecular declared 2 {a} 1 vibration 1 1300 1 power 2.0265 293.15 0 200 500 {b} 1 vibration 1 1900 1 power 32.424 293.15 0 200 500\nmolecular_resolution 8 4096\nmolecular_wall none\n"
    )
}

fn input() -> String {
    let source = include_str!("../../../../../examples/reed-molecular-synthetic.performance")
        .replace("audio 48000 9606", "audio 8000 1601")
        .replace(
            "ambient mixture 293.15 101325 mole 2 N2 0.7 O2 0.3",
            "ambient 293.15 101325 0",
        );
    let start = source.find("molecular declared ").unwrap();
    let end = source.find("compile_limits ").unwrap();
    format!(
        "{}{}{}",
        &source[..start],
        declaration("A", "B"),
        &source[end..]
    )
}

fn load(text: &str, sources: &ThermochemicalSources<'_>, block: usize) -> ReedPerformance {
    ReedPerformance::from_bytes_with_gas_sources(
        text.as_bytes(),
        block,
        ObservationPolicy::Diagnostic,
        Some(sources),
    )
    .unwrap()
}

fn render(mut source: ReedPerformance, block: usize) -> Vec<f64> {
    let mut pressure = vec![0.0; source.info().samples as usize];
    for part in pressure.chunks_mut(block) {
        source.pressure_block(part).unwrap();
    }
    assert_eq!(source.remaining_samples(), 0);
    assert!(pressure.iter().all(|p| p.is_finite()));
    pressure
}

#[test]
fn continuous_molecular_field_uses_its_own_sources_for_feedback_and_microphone() {
    let (a, b) = (packs("A", 0.028, 3.5, 200.0), packs("B", 0.044, 4.5, 200.0));
    let bindings = [a.binding(), b.binding()];
    let sources = ThermochemicalSources {
        species: &bindings,
        water_species: None,
        phase: ThermochemicalPhase::FrozenGas,
    };
    let text = input().replace(
        "molecular declared ",
        &format!("{FIELD}molecular declared "),
    );
    let performance = load(&text, &sources, 37);
    let field = performance.gas_field().unwrap().0.clone();
    let resolution = GasFieldResolution::Continuous(fs_duct::RiccatiSteps(8));
    assert_eq!(performance.gas_field_resolution(), Some(resolution));
    let point = performance.molecular_driving_point().unwrap();
    assert_eq!(point.material().gas_field(), Some(&field));
    assert_eq!(point.inlet_gas(), field.state_at(0.0).unwrap().state());
    assert_ne!(point.inlet_gas(), performance.ambient_gas());
    assert_eq!(point.exterior_gas(), field.state_at(0.15).unwrap().state());
    assert_eq!(
        performance.outlet_receiver().unwrap().gas,
        *point.exterior_gas()
    );
    assert_eq!(performance.segment_gases(), &[*point.inlet_gas()]);
    assert_eq!(performance.gas_records().len(), 4);
    assert_eq!(field.model().source_packs().len(), 2);
    // An injected field takes the same molecular route and physical samples.
    // Its source field retains receipts, while gas_records keeps its existing
    // authored-record meaning. The explicit lossless observer is identical.
    let injected = ReedPerformance::from_bytes_with_field_resolution(
        input().as_bytes(),
        113,
        ObservationPolicy::Diagnostic,
        Some(&sources),
        field,
        resolution,
    )
    .unwrap();
    assert_eq!(injected.gas_records().len(), 1);
    let baseline = render(performance, 37);
    assert_eq!(baseline, render(injected, 113));
    assert!(baseline.iter().any(|p| *p != 0.0));
    let faster = text.replace("power 2.0265", "power 0.20265");
    let changed = render(load(&faster, &sources, 37), 37);
    let energy: f64 = baseline.iter().map(|p| p * p).sum();
    let difference: f64 = baseline
        .iter()
        .zip(&changed)
        .map(|(p, q)| (p - q).powi(2))
        .sum();
    assert!(
        difference > 1e-8 * energy,
        "field molecular rates must reach physical Pa"
    );
}

#[test]
fn molecular_fields_refuse_midpoint_images_and_reservoirs_missing_downstream() {
    let (a, b) = (packs("A", 0.028, 3.5, 200.0), packs("B", 0.044, 4.5, 200.0));
    let bindings = [a.binding(), b.binding()];
    let sources = ThermochemicalSources {
        species: &bindings,
        water_species: None,
        phase: ThermochemicalPhase::FrozenGas,
    };
    let text = input().replace(
        "molecular declared ",
        &format!("{FIELD}molecular declared "),
    );
    let midpoint = text.replace("gas_field continuous 8", "gas_field midpoint 8");
    assert!(
        matches!(ReedPerformance::from_bytes_with_gas_sources(midpoint.as_bytes(), 37,
        ObservationPolicy::Diagnostic, Some(&sources)),
        Err(ReedPerformanceError::Input { what, .. }) if what.contains("explicit continuous"))
    );
    // B is absent at the inlet but present downstream. Retaining only the
    // evaluated inlet state would silently erase its molecular law there.
    let field_text = FIELD.replace("300 0.8 0.2", "300 1 0");
    let missing = declaration("A", "B")
        .replace("declared 2", "declared 1")
        .replace(" B 1 vibration 1 1900 1 power 32.424 293.15 0 200 500", "");
    let without_law = input().replace(&declaration("A", "B"), &missing);
    let missing = without_law.replace(
        "molecular declared ",
        &format!("{field_text}molecular declared "),
    );
    assert!(matches!(
        ReedPerformance::from_bytes_with_gas_sources(
            missing.as_bytes(),
            37,
            ObservationPolicy::Diagnostic,
            Some(&sources)
        ),
        Err(ReedPerformanceError::Model(_))
    ));
    let out_of_domain = text
        .replace("gas_station 0.075 330", "gas_station 0.075 390")
        .replace("0 200 500", "0 200 340");
    assert!(matches!(
        ReedPerformance::from_bytes_with_gas_sources(
            out_of_domain.as_bytes(),
            37,
            ObservationPolicy::Diagnostic,
            Some(&sources)
        ),
        Err(ReedPerformanceError::Model(_))
    ));
    // Parser admission itself is attached to the actual field inlet, even
    // though the legacy ambient has neither declared source species.
    let mut records = Vec::new();
    let mut reader = Reader {
        lines: FIELD.lines(),
        line: 0,
    };
    let (field, _) =
        field_input::read(&mut reader, Some(&sources), 101325.0, &mut records).unwrap();
    let ambient = GasRecord::Legacy(GasState::try_new_moist_air(293.15, 101325.0, 0.0).unwrap());
    let law = declaration("A", "B");
    let mut reader = Reader {
        lines: law.lines(),
        line: 0,
    };
    let parsed = read(&mut reader, &ambient, Some(&field)).unwrap();
    assert_eq!(parsed.state.gas(), field.state_at(0.0).unwrap().state());
    assert_eq!(parsed.law.species().len(), 2);
}

#[test]
fn steady_generated_gas_and_molecular_reservoirs_share_the_existing_performance() {
    let (a, b) = (
        packs("N2", 0.028, 3.5, 200.0),
        packs("CO2", 0.044, 4.5, 200.0),
    );
    let bindings = [a.binding(), b.binding()];
    let sources = ThermochemicalSources {
        species: &bindings,
        water_species: None,
        phase: ThermochemicalPhase::FrozenGas,
    };
    let text = include_str!("../../../../../examples/reed-steady-gas.performance")
        .replace("audio 48000 12000", "audio 8000 1601")
        .replace("observer_law stokes-kirchhoff", "observer_law none")
        .replace(
            "compile_limits ",
            &format!("{}compile_limits ", declaration("N2", "CO2")),
        );
    let performance = load(&text, &sources, 37);
    let solution = performance.steady_gas().unwrap().solution();
    let field = solution.field();
    let point = performance.molecular_driving_point().unwrap();
    assert_eq!(performance.gas_field().unwrap().0, field);
    assert_eq!(point.material().gas_field(), Some(field));
    assert_eq!(point.inlet_gas(), field.state_at(0.0).unwrap().state());
    assert_eq!(
        point.exterior_gas(),
        field.state_at(field.length_m()).unwrap().state()
    );
    assert_eq!(
        performance.outlet_receiver().unwrap().gas,
        *point.exterior_gas()
    );
    assert_eq!(field.model().source_packs().len(), 2);
    assert!(field.uniform_state().is_none());
    assert!(
        solution.maximum_energy_residual_j()
            <= performance
                .steady_gas()
                .unwrap()
                .controls()
                .energy_tolerance_j
    );
    assert_eq!(performance.renderer().samples_rendered(), 0);
}
