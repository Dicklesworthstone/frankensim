//! Exercises the real file -> pack decoder -> material resolver -> spatial
//! electrical/thermal solver -> command output path. Synthetic data only.

mod directional;

use super::*;
use fs_blake3::hash_bytes;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, NormalizedPack, ObservationDataset,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, UncertaintyModel,
};
use fs_material::conductor::{ELECTRICAL_RESISTIVITY_DIMS, ELECTRICAL_RESISTIVITY_PROPERTY};
use fs_material::state_point::{THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY};
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

/// Registers the one raw observation a compiled pack must retain. For these
/// fixtures it is the authored analytic law itself, not a measurement.
pub(super) fn authored_observation(claims: &mut ClaimSet, source: &str) {
    claims.register_observation(ObservationDataset {
        specimen: "authored synthetic fixture".into(), method: "analytic law, authored for tests".into(),
        artifact: hash_bytes(source.as_bytes()), caveats: "not measured; test data only".into(),
        provenance: Provenance { source: source.into(), license: "CC0: authored synthetic test data".into(),
            artifact: None },
    }).unwrap();
}

fn material_bytes(rho: f64, k: f64, a: f64, b: f64, maximum_k: f64) -> Vec<u8> {
    let mut claims = ClaimSet::new();
    for (name, dims, base, slope) in [
        (ELECTRICAL_RESISTIVITY_PROPERTY, ELECTRICAL_RESISTIVITY_DIMS, rho, a),
        (THERMAL_CONDUCTIVITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS, k, b),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Curve {
                abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                knots: vec![(290.0, base * (1.0 - 10.0 * slope)),
                    (maximum_k, base * (1.0 + slope * (maximum_k - 300.0)))], dims,
            },
            validity: ValidityDomain::unconstrained().with("T", 290.0, maximum_k)
                .with("pressure", 101_325.0, 101_325.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: Vec::new(),
            provenance: Provenance {
                source: "authored analytic resistor-network fixture, not measured metal".into(),
                license: "CC0: authored synthetic test data".into(), artifact: None,
            },
        }).unwrap();
    }
    authored_observation(&mut claims, "authored analytic resistor-network fixture, not measured metal");
    let pack = NormalizedPack::new(
        "synthetic-electrothermal-fixture", "test-fixture-v1", hash_bytes(b"synthetic analytic data"),
        "CC0: authored synthetic test data", claims, Vec::new(), Vec::new(),
    ).unwrap();
    NormalizedMaterialCardPack::new(MaterialStateId {
        chemistry: "synthetic fixture".into(), phase: "solid".into(),
        process: "authored analytic response".into(), revision: 0,
    }, pack).unwrap().to_bytes()
}

fn scene() -> String {
    r#"{
      "schema":"frankensim.electrothermal.v1",
      "cells":[3,1,1], "cell_size_m":1.0,
      "materials":[{"pack":"material.fsmcdpk", "conducting":true,
        "temperature_axis":"T", "temperature_kind":"legacy",
        "state":[{"name":"pressure", "value":"101325 Pa", "kind":"legacy"}]}],
      "material_indices":[0,0,0], "initial_temperature_k":300.0,
      "thermal_faces":[{"temperature_k":300.0},{"temperature_k":300.0},
        "adiabatic","adiabatic","adiabatic","adiabatic"],
      "electrodes":[
        {"lo":[0,0,0],"hi":[1,1,1],"voltage_v":0.0},
        {"lo":[2,0,0],"hi":[3,1,1],"voltage_v":3.0}],
      "numerics":{"max_passes":100,"max_thermal_passes":100,
        "temperature_tolerance_k":1e-9,"power_tolerance_w":1e-9,
        "relative_tolerance":1e-10,"conductivity_relative_tolerance":1e-10,
        "energy_tolerance":1e-12,"energy_max_iterations":2000}
    }"#.into()
}

fn fixture(bytes: &[u8], scene: &str) -> PathBuf {
    static NEXT: AtomicU64 = AtomicU64::new(0);
    let id = NEXT.fetch_add(1, Ordering::Relaxed);
    let timestamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let dir = std::env::temp_dir().join(format!("frankensim-electrothermal-{}-{timestamp}-{id}", std::process::id()));
    std::fs::create_dir(&dir).unwrap();
    std::fs::write(dir.join("material.fsmcdpk"), bytes).unwrap();
    let path = dir.join("scene.json");
    std::fs::write(&path, scene).unwrap();
    path
}

fn execute(path: &Path) -> (CommandOutput, J) {
    let output = run(&[path.as_os_str().to_owned()], true);
    assert_eq!(output.exit_code, exit::SUCCESS, "{}", output.stderr);
    assert!(output.stderr.is_empty());
    let result = J::parse(&output.stdout).unwrap();
    (output, result)
}
fn component(value: &J, name: &str, index: usize) -> f64 {
    value.get(name).unwrap().as_array().unwrap()[index].as_f64().unwrap()
}

#[test]
fn g1_g3_file_only_material_swap_changes_current_heating_and_temperature() {
    let text = scene();
    for (rho, k, expected_t, expected_p) in [(2.0, 1.0, 303.375, 4.5), (4.0, 2.0, 300.84375, 2.25)] {
        let bytes = material_bytes(rho, k, 0.0, 0.0, 350.0);
        let path = fixture(&bytes, &text);
        let (output, result) = execute(&path);
        assert!((component(&result, "temperature_k", 1) - expected_t).abs() < 1e-8);
        assert!((result.f64_field("total_joule_w").unwrap() - expected_p).abs() < 1e-9);
        assert!((result.f64_field("electrode_power_w").unwrap() - expected_p).abs() < 1e-9);
        assert!(result.f64_field("thermal_balance_relative_residual").unwrap() < 1e-9);
        let sources = result.get("source_states").unwrap().as_array().unwrap();
        assert_eq!(sources.len(), 3);
        assert!(sources[1].str_field("electrical_state").is_some());
        assert_ne!(sources[0].str_field("thermal_state"), sources[1].str_field("thermal_state"));
        // Replay is byte-identical and does not rewrite the source scene/pack.
        assert_eq!(output, run(&[path.as_os_str().to_owned()], true));
        assert_eq!(std::fs::read(&path).unwrap(), text.as_bytes());
        assert_eq!(std::fs::read(path.parent().unwrap().join("material.fsmcdpk")).unwrap(), bytes);
    }
}

#[test]
fn g1_both_temperature_dependent_properties_close_the_same_spatial_problem() {
    let path = fixture(&material_bytes(2.0, 1.0, 0.02, 0.03, 350.0), &scene());
    let (_, result) = execute(&path);
    // Independent two-node symmetric resistor-network reduction:
    // P=9/[2(1+.02*d)], x=P/[4(1+.03*x)],
    // d=P/[2(1+.03*x)]+P/[4(1+.03*d)].
    let expected = [301.028_953_367_629_24, 303.030_231_392_828_7, 301.028_953_367_629_24];
    for (cell, expected_t) in expected.into_iter().enumerate() {
        let temperature = component(&result, "temperature_k", cell);
        assert!((temperature - expected_t).abs() < 1e-7);
        let resolved = component(&result, "resolved_conductivity_w_m_k", cell);
        assert!((resolved - (1.0 + 0.03 * (temperature - 300.0))).abs() < 1e-12);
    }
    assert!((result.f64_field("total_joule_w").unwrap() - 4.242_862_874_447_63).abs() < 1e-8);
    assert!(result.f64_field("conductivity_relative_defect").unwrap() <= 1e-10);
    assert!(result.f64_field("passes").unwrap() > 1.0);
    assert!(result.f64_field("thermal_solves").unwrap() > result.f64_field("passes").unwrap());
}

#[test]
fn g0_source_support_refusal_publishes_no_partial_result() {
    let path = fixture(&material_bytes(2.0, 1.0, 0.0, 0.0, 300.5), &scene());
    let output = run(&[path.as_os_str().to_owned()], true);
    assert_eq!(output.exit_code, exit::REFUSED);
    assert!(output.stdout.is_empty());
    assert!(J::parse(&output.stderr).unwrap().str_field("message").unwrap().contains("material.thermal_conductivity"));
}

#[test]
fn g0_malformed_scene_or_pack_cannot_enter_the_solver() {
    let bytes = material_bytes(2.0, 1.0, 0.0, 0.0, 350.0);
    for bad in [
        scene().replace("\"material_indices\":[0,0,0]", "\"material_indices\":[0,0]"),
        scene().replace("\"material_indices\":[0,0,0]", "\"material_indices\":[0,0,1]"),
        scene().replace("\"max_passes\":100", "\"max_passes\":0"),
        scene().replace("\"cell_size_m\":1.0", "\"cell_size_m\":0.0"),
        scene().replace("\"schema\":", "\"unexpected\":true,\"schema\":"),
        scene().replace("\"conducting\":true", "\"conducting\":\"true\""),
    ] {
        let path = fixture(&bytes, &bad);
        let output = run(&[path.as_os_str().to_owned()], true);
        assert_eq!(output.exit_code, exit::INPUT, "{}", output.stderr);
        assert!(output.stdout.is_empty());
    }
    let path = fixture(b"not a canonical material pack", &scene());
    let output = run(&[path.as_os_str().to_owned()], true);
    assert_eq!(output.exit_code, exit::INPUT);
    assert!(output.stdout.is_empty());
}

#[test]
fn g0_explicit_coordinate_conventions_are_not_silently_erased() {
    let bytes = material_bytes(2.0, 1.0, 0.0, 0.0, 350.0);
    // A semantic query cannot match the fixture's legacy source axis simply
    // because both happen to be kelvin. The material owner must refuse it.
    let text = scene().replace("\"temperature_kind\":\"legacy\"", "\"temperature_kind\":\"absolute-temperature\"");
    let path = fixture(&bytes, &text);
    let output = run(&[path.as_os_str().to_owned()], true);
    assert_eq!(output.exit_code, exit::REFUSED);
    assert!(output.stdout.is_empty());
    assert!(quantity("temperature-difference", Temperature::DIMS).is_ok());
    assert!(quantity("pressure", Temperature::DIMS).is_err());
}
