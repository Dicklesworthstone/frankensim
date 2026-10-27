//! Native and live callers must consume the same canonical material sources.
//! Synthetic gases test source plumbing and physical coupling, not gas accuracy.
#[path = "../src/bin/music_render/thermochemical.rs"]
mod thermochemical;
#[path = "../src/bin/music_render/json.rs"]
mod json;
use json::json_string;
#[path = "../src/render/schedule/reed/field_input/tests/support.rs"]
mod support;

use fs_couple::pcm_wav::observation::{ObservationPolicy, PressureRenderer};
use fs_couple::render::schedule::reed::ReedPerformance;
use fs_material::gas_species::input::{
    GasRecord, ThermochemicalPhase, read_gas_record,
    packs::{SourcePackBytes, read_source_packs},
};
use std::collections::BTreeMap;
use std::path::Path;

fn options() -> thermochemical::Options {
    thermochemical::Options::extract(&[
        "--gas-component", "gas/species", "gas/nasa9", "gas/transport",
        "--gas-water", "none", "--gas-phase", "frozen",
    ].map(String::from)).unwrap().1
}

fn files(packs: &support::Packs, legacy: bool) -> BTreeMap<String, Vec<u8>> {
    [
        ("gas/species".into(), packs.species.to_bytes()),
        ("gas/nasa9".into(), packs.nasa9.to_bytes()),
        ("gas/transport".into(), if legacy {
            packs.transport.claims_pack().to_bytes()
        } else { packs.transport.to_bytes() }),
    ].into_iter().collect()
}
fn load(files: &BTreeMap<String, Vec<u8>>) -> Result<thermochemical::LoadedSources, String> {
    options().load_with(|path: &Path| files.get(path.to_str().unwrap()).cloned()
        .ok_or_else(|| format!("unattached gas source {}", path.display())))
}

const GAS: &str = "thermochemical 300 101325 mole 1 custom-alpha 1";

#[test]
fn canonical_transport_retains_the_live_gas_state_and_all_original_pack_identities() {
    let packs = support::packs("custom-alpha", 0.028, 3.5, 250.0);
    let source_files = files(&packs, false);
    let native = load(&source_files).unwrap();
    let bindings = native.bindings();
    let registry = native.registry(&bindings).unwrap();
    let record = read_gas_record(&mut GAS.split_ascii_whitespace(), None, Some(&registry)).unwrap();
    let live = read_source_packs(GAS, &[SourcePackBytes {
        species: &source_files["gas/species"], nasa9: &source_files["gas/nasa9"],
        transport: &source_files["gas/transport"],
    }], None, ThermochemicalPhase::FrozenGas).unwrap();
    assert_eq!(record, live);
    let GasRecord::Thermochemical(gas) = &record else { panic!("source lost") };
    assert_eq!(gas.source_packs()[0].transport_pack, packs.transport.content_hash());
    assert_eq!(bindings[0].transport.as_ref().unwrap().0, &packs.transport);

    // Historical raw claims retain their old wrapper, with identical physics.
    let legacy = load(&files(&packs, true)).unwrap();
    let legacy_bindings = legacy.bindings();
    let legacy_registry = legacy.registry(&legacy_bindings).unwrap();
    let old = read_gas_record(&mut GAS.split_ascii_whitespace(), None, Some(&legacy_registry)).unwrap();
    assert_eq!(record.state(), old.state());
    assert_ne!(record, old, "the explicit material-card identity must not be rewritten");
}

fn render(mass: f64, block: usize) -> Vec<f64> {
    let packs = support::packs("custom-alpha", mass, 3.5, 250.0);
    let loaded = load(&files(&packs, false)).unwrap();
    let bindings = loaded.bindings();
    let registry = loaded.registry(&bindings).unwrap();
    let text = include_str!("../examples/reed-duct.performance")
        .replace("ambient 293.15 101325 0", &format!("ambient {GAS}"));
    let mut voice = ReedPerformance::from_bytes_with_gas_sources(
        text.as_bytes(), block, ObservationPolicy::Diagnostic, Some(&registry),
    ).unwrap();
    let mut pressure = vec![0.0; voice.info().samples as usize];
    for chunk in pressure.chunks_mut(block) { voice.block(chunk).unwrap(); }
    assert!(pressure.iter().all(|value| value.is_finite()));
    assert!(pressure.iter().any(|value| value.abs() > 1e-6));
    pressure
}

#[test]
fn changing_only_gas_source_data_changes_the_native_reed_pressure_not_just_metadata() {
    let baseline = render(0.028, 37);
    assert_eq!(baseline, render(0.028, 512));
    assert_ne!(baseline, render(0.032, 37));
}

#[test]
fn missing_corrupt_and_oversized_sources_never_fall_back_to_air_or_filesystem() {
    let packs = support::packs("custom-alpha", 0.028, 3.5, 250.0);
    let mut source_files = files(&packs, false);
    source_files.remove("gas/transport");
    assert!(load(&source_files).err().unwrap().contains("unattached gas source"));
    source_files.insert("gas/transport".into(), b"FSMCDPK\0corrupt".to_vec());
    let error = load(&source_files).err().unwrap();
    assert!(error.contains("material-card pack"), "{error}");
    assert!(!error.contains("claim pack"), "{error}");
    let error = options().load_with(|_| Ok(vec![0; 16 * 1024 * 1024 + 1])).err().unwrap();
    assert!(error.contains("16 MiB"), "{error}");
    assert!(options().load_with(|_| Ok(vec![])).is_err());
}
