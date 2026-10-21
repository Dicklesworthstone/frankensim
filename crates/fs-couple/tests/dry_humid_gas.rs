//! A common source registry must admit a dry exterior and humid instrument gas.
//! Synthetic thermodynamic/transport packs exercise wiring, not measured air data.
use fs_couple::pcm_wav::observation::{ObservationPolicy, PressureRenderer};
use fs_couple::render::schedule::reed::ReedPerformance;
use fs_material::gas_species::input::{
    GasInputError, GasRecord, ThermochemicalPhase, ThermochemicalSources, read_gas_record,
};

// Reuse the existing synthetic NASA-9/transport pack fixture, not a new gas law.
#[path = "../src/render/schedule/reed/field_input/tests/support.rs"]
mod support;
use support::packs;

fn base() -> String {
    include_str!("../examples/reed-flanged-baffled.performance")
        .replace("audio 48000 12000 2", "audio 48000 4801 2")
        .replace("ambient mixture 293.15 100000 mole 2 N2 0.7 O2 0.3", "ambient 300 100000 0")
        .replace("reed 0.0004 0.013 6000 0 0 0.35", "reed 0.0004 0.013 6000 0.00001 100 0.35")
        .replace("segments 2\ncylinder 0.0015 0.2\ncylinder 0.0015 0.3 gas mixture 293.15 100000 mole 3 N2 0.4 O2 0.2 CO2 0.4",
            "segments 1\ncylinder 0.0015 0.57")
}
fn load(text: &str, sources: &ThermochemicalSources<'_>, block: usize) -> ReedPerformance {
    ReedPerformance::from_bytes_with_gas_sources(
        text.as_bytes(), block, ObservationPolicy::Diagnostic, Some(sources),
    ).unwrap()
}
fn render(mut performance: ReedPerformance, block: usize) -> Vec<f64> {
    let mut pressure = vec![0.0; performance.info().samples as usize];
    for chunk in pressure.chunks_mut(block) { performance.block(chunk).unwrap(); }
    assert!(pressure.iter().all(|p| p.is_finite()));
    pressure
}

#[test]
fn known_dry_records_keep_exact_source_receipts_without_a_fictitious_water_component() {
    let (dry, vapor) = (packs("A", 0.028, 3.5, 200.0), packs("B", 0.018, 4.0, 200.0));
    let bindings = [dry.binding(), vapor.binding()];
    // B, not a name heuristic such as H2O, is the explicitly bound water identity.
    let water = vapor.species.association().species();
    let checked = ThermochemicalSources { species: &bindings, water_species: Some(water),
        phase: ThermochemicalPhase::LiquidWaterLimit };
    let frozen = ThermochemicalSources { phase: ThermochemicalPhase::FrozenGas, ..checked };
    for basis in ["mole", "mass"] {
        // A known-dry gas does not need a liquid-water saturation extrapolation
        // at 250 K. Its own NASA/transport domains still govern admissibility.
        for temperature in [250.0, 300.0] {
            let text = format!("thermochemical {temperature} 100000 {basis} 1 A 1 suffix");
            let mut fields = text.split_ascii_whitespace();
            let actual = read_gas_record(&mut fields, Some(200000.0), Some(&checked)).unwrap();
            let reference = read_gas_record(&mut text.split_ascii_whitespace(),
                Some(200000.0), Some(&frozen)).unwrap();
            assert_eq!(actual, reference);
            assert_eq!(fields.collect::<Vec<_>>(), ["suffix"]);
            let GasRecord::Thermochemical(gas) = actual else { panic!("original source receipt"); };
            assert_eq!(gas.state().pressure, 100000.0);
            assert_eq!(gas.state().water_mole_fraction, 0.0);
            assert_eq!(gas.components().len(), 1);
            assert_eq!(gas.source_packs().len(), 1);
            assert_eq!(gas.source_packs()[0].species_pack, dry.species.content_hash());
        }
    }
}

#[test]
fn dry_admission_does_not_disable_water_identity_saturation_or_source_domain_checks() {
    let (dry, vapor) = (packs("A", 0.028, 3.5, 200.0), packs("B", 0.018, 4.0, 200.0));
    let bindings = [dry.binding(), vapor.binding()];
    let checked = ThermochemicalSources { species: &bindings,
        water_species: Some(vapor.species.association().species()),
        phase: ThermochemicalPhase::LiquidWaterLimit };
    let unknown = ThermochemicalSources { water_species: None, ..checked };
    for (text, sources) in [
        ("thermochemical 300 100000 mole 2 A 0.5 B 0.5", &checked),
        ("thermochemical 199 100000 mole 1 A 1", &checked),
        ("thermochemical 300 100000 mole 1 A 1", &unknown),
        ("thermochemical 300 100000 mole 2 A 0.99 B 0.01", &unknown),
    ] {
        let mut fields = text.split_ascii_whitespace();
        assert!(matches!(read_gas_record(&mut fields, None, Some(sources)),
            Err(GasInputError::Material(_))));
        assert_eq!(fields.collect::<Vec<_>>(), text.split_ascii_whitespace().collect::<Vec<_>>());
    }
    let wet = read_gas_record(&mut "thermochemical 300 100000 mole 2 A 0.99 B 0.01"
        .split_ascii_whitespace(), None, Some(&checked)).unwrap();
    let GasRecord::Thermochemical(wet) = wet else { panic!("original source receipt"); };
    assert_eq!(wet.source_packs().len(), 2);
    assert!((wet.state().water_mole_fraction-0.01).abs() < 1e-15);
}

#[test]
fn dry_exterior_and_humid_bore_share_one_policy_in_actual_reed_pressure_rendering() {
    let (dry, vapor) = (packs("A", 0.028, 3.5, 200.0), packs("B", 0.018, 4.0, 200.0));
    let bindings = [dry.binding(), vapor.binding()];
    let checked = ThermochemicalSources { species: &bindings,
        water_species: Some(vapor.species.association().species()),
        phase: ThermochemicalPhase::LiquidWaterLimit };
    let frozen = ThermochemicalSources { phase: ThermochemicalPhase::FrozenGas, ..checked };
    let text = base().replace("ambient 300 100000 0",
        "ambient thermochemical 300 100000 mole 2 A 0.99 B 0.01")
        .replace("termination flanged",
        "exterior thermochemical 300 100000 mole 1 A 1\ntermination flanged");
    let performance = load(&text, &checked, 37);
    let reference = load(&text, &frozen, 113);
    assert_eq!(performance.gas_records(), reference.gas_records());
    assert_eq!(performance.outlet_receiver().unwrap().gas.water_mole_fraction, 0.0);
    let actual = render(performance, 37);
    assert!(actual.iter().any(|p| *p != 0.0));
    assert_eq!(actual, render(reference, 113));
    let all_dry = text.replace("mole 2 A 0.99 B 0.01", "mole 1 A 1");
    let dry_output = render(load(&all_dry, &checked, 37), 37);
    assert_ne!(actual, dry_output, "bore composition must change actual rendered pressure");
}
