//! End-to-end ordered interface packs -> spatial current -> heat -> temperature.
//! All numerical properties are synthetic analytic fixtures, not measured metals.

use super::*;
use fs_blake3::hash_bytes;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, NormalizedInterfacePack,
    NormalizedPack, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    SurfaceSpec, SystemContext, UncertaintyModel,
};
use fs_material::conductor::{ELECTRICAL_RESISTIVITY_DIMS, ELECTRICAL_RESISTIVITY_PROPERTY};
use fs_material::electrical_interface::{
    ELECTRICAL_CONTACT_RESISTANCE_DIMS as RESISTANCE_DIMS,
    ELECTRICAL_CONTACT_RESISTANCE_PROPERTY as PROPERTY,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS, SPECIFIC_HEAT_CAPACITY_PROPERTY,
    THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY,
};
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

const SCENE: &str = r#"{
  "schema":"frankensim.electrothermal.v1", "cells":[4,1,1], "cell_size_m":1.0,
  "materials":[
    {"pack":"a.fsmcdpk","conducting":true,"temperature_axis":"T","temperature_kind":"legacy",
     "state":[{"name":"pressure","value":"101325 Pa","kind":"legacy"}]},
    {"pack":"b.fsmcdpk","conducting":true,"temperature_axis":"T","temperature_kind":"legacy",
     "state":[{"name":"pressure","value":"101325 Pa","kind":"legacy"}]}],
  "material_indices":[0,0,1,1], "initial_temperature_k":300.0,
  "thermal_faces":[{"temperature_k":300.0},{"temperature_k":300.0},
    "adiabatic","adiabatic","adiabatic","adiabatic"],
  "electrodes":[{"lo":[0,0,0],"hi":[1,1,1],"voltage_v":0.0},
    {"lo":[3,0,0],"hi":[4,1,1],"voltage_v":3.0}],
  "electrical_interfaces":[{"pack":"joint.fsintpk","faces":[[1,2]],
    "temperature_axes":["Ta","Tb"],"temperature_kind":"legacy",
    "state":[{"name":"normal_pressure","value":"100000 Pa","kind":"legacy"}],
    "heat_fraction_to_a":0.25}],
  "numerics":{"max_passes":100,"max_thermal_passes":100,
    "temperature_tolerance_k":1e-9,"power_tolerance_w":1e-9,"relative_tolerance":1e-10,
    "conductivity_relative_tolerance":1e-10,"energy_tolerance":1e-12,"energy_max_iterations":2000}
}"#;

fn material_id(side: &str) -> MaterialStateId {
    MaterialStateId { chemistry: format!("synthetic contact solid {side}"),
        phase: "solid".into(), process: "analytic test fixture".into(), revision: 0 }
}

fn provenance() -> Provenance {
    Provenance { source: "authored analytic contact fixture, not measurements".into(),
        license: "CC0: synthetic test data".into(), artifact: None }
}

fn pack(mut claims: ClaimSet, name: &str) -> NormalizedPack {
    super::tests::authored_observation(&mut claims, "authored analytic contact fixture, not measurements");
    NormalizedPack::new(name, "contact-test-v1", hash_bytes(b"authored contact fixture"),
        "CC0: synthetic test data", claims, Vec::new(), Vec::new()).unwrap()
}

fn bulk(side: &str) -> Vec<u8> {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (ELECTRICAL_RESISTIVITY_PROPERTY, ELECTRICAL_RESISTIVITY_DIMS, 2.0),
        (THERMAL_CONDUCTIVITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS, 1.0),
        (DENSITY_PROPERTY, fs_qty::Density::DIMS, 1.0),
        (SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS, 1.0),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims), value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 290.0, 350.0)
                .with("pressure", 101_325.0, 101_325.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: Vec::new(), provenance: provenance(),
        }).unwrap();
    }
    NormalizedMaterialCardPack::new(material_id(side), pack(claims, side)).unwrap().to_bytes()
}

fn joint(rc: f64, slope: f64, upper_b: f64, reversed: bool, dims: Dims) -> Vec<u8> {
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(PROPERTY, dims),
        value: PropertyValue::Curve { abscissa: "Ta".into(), abscissa_dims: Temperature::DIMS,
            knots: vec![(290.0, rc * (1.0 - 10.0*slope)), (350.0, rc * (1.0 + 50.0*slope))], dims },
        validity: ValidityDomain::unconstrained().with("Ta", 290.0, 350.0)
            .with("Tb", 290.0, upper_b).with("normal_pressure", 100_000.0, 200_000.0),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
        observations: Vec::new(), provenance: provenance(),
    }).unwrap();
    NormalizedInterfacePack::new(
        SurfaceSpec { material: material_id(if reversed {"B"} else {"A"}), texture_frame: "declared-A-cell-state".into() },
        SurfaceSpec { material: material_id(if reversed {"A"} else {"B"}), texture_frame: "declared-B-cell-state".into() },
        SystemContext { medium: "synthetic joint".into(), third_body: None,
            environment: "two adjacent cell-temperature samples".into(), history: "fixed contact fixture".into() },
        pack(claims, "joint"),
    ).unwrap().to_bytes()
}

fn fixture(text: &str, joint: &[u8]) -> PathBuf {
    static NEXT: AtomicU64 = AtomicU64::new(0);
    let timestamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let dir = std::env::temp_dir().join(format!("frankensim-contact-{}-{timestamp}-{}",
        std::process::id(), NEXT.fetch_add(1, Ordering::Relaxed)));
    std::fs::create_dir(&dir).unwrap();
    std::fs::write(dir.join("a.fsmcdpk"), bulk("A")).unwrap();
    std::fs::write(dir.join("b.fsmcdpk"), bulk("B")).unwrap();
    std::fs::write(dir.join("joint.fsintpk"), joint).unwrap();
    let path = dir.join("scene.json");
    std::fs::write(&path, text).unwrap();
    path
}

fn execute(path: &Path) -> (CommandOutput, J) {
    let output = run(&[path.as_os_str().to_owned()], true);
    assert_eq!(output.exit_code, exit::SUCCESS, "{}", output.stderr);
    let parsed = J::parse(&output.stdout).unwrap();
    (output, parsed)
}
fn component(value: &J, key: &str, index: usize) -> f64 {
    value.get(key).unwrap().as_array().unwrap()[index].as_f64().unwrap()
}
fn interface(value: &J) -> &J {
    &value.get("electrical_interfaces").unwrap().as_array().unwrap()[0]
}

#[test]
fn g1_g3_interface_pack_only_swaps_change_current_local_heat_and_temperature() {
    for rc in [0.0, 2.0, 6.0] {
        let bytes = joint(rc, 0.0, 350.0, false, RESISTANCE_DIMS);
        let path = fixture(SCENE, &bytes);
        let (output, result) = execute(&path);
        let current = 3.0 / (4.0 + rc);
        let qa = current * current * (2.0 + 0.25 * rc);
        let qb = current * current * (2.0 + 0.75 * rc);
        let expected_rises = [
            (5.0*qa+3.0*qb)/16.0, (15.0*qa+9.0*qb)/16.0,
            (9.0*qa+15.0*qb)/16.0, (3.0*qa+5.0*qb)/16.0,
        ]; // exact inverse of the four-cell unit-conductivity thermal matrix
        for (cell, rise) in expected_rises.into_iter().enumerate() {
            assert!((component(&result, "temperature_k", cell) - 300.0 - rise).abs() < 1e-8);
        }
        assert!((component(&result, "joule_w", 1) - qa).abs() < 1e-10);
        assert!((component(&result, "joule_w", 2) - qb).abs() < 1e-10);
        assert!((result.f64_field("total_joule_w").unwrap() - 3.0*current).abs() < 1e-10);
        assert!((result.f64_field("electrode_power_w").unwrap() - 3.0*current).abs() < 1e-10);
        assert!(result.f64_field("thermal_balance_relative_residual").unwrap() < 1e-9);
        assert_eq!(interface(&result).f64_field("resistance_ohm_m2"), Some(rc));
        assert_eq!(interface(&result).str_field("pack_identity").unwrap(),
            NormalizedInterfacePack::from_bytes(&bytes).unwrap().content_hash().to_hex());
        assert_eq!(output, run(&[path.as_os_str().to_owned()], true));
        assert_eq!(std::fs::read(&path).unwrap(), SCENE.as_bytes());
        assert_eq!(std::fs::read(path.parent().unwrap().join("joint.fsintpk")).unwrap(), bytes);
    }
}

#[test]
fn g1_temperature_dependent_joint_is_requeried_at_the_returned_state() {
    let path = fixture(SCENE, &joint(2.0, 0.02, 350.0, false, RESISTANCE_DIMS));
    let problem = load(&path).unwrap();
    let initial = interfaces::resolve(&problem.electrical_interfaces, &[300.0;4], &CancelGate::new()).unwrap();
    let initial = J::parse(&interfaces::render(&problem, &initial)).unwrap();
    let (_, result) = execute(&path);
    // Independent nonlinear network: R=2*(1+.02*(Ta-300)), I=3/(4+R).
    let expected = [300.356_610_279_130_73, 301.069_830_837_392_2,
        301.164_235_343_321_8, 300.388_078_447_773_9];
    for (cell, t) in expected.into_iter().enumerate() {
        assert!((component(&result, "temperature_k", cell) - t).abs() < 1e-8);
    }
    let rc = interface(&result).f64_field("resistance_ohm_m2").unwrap();
    assert!((rc - 2.0*(1.0+0.02*(component(&result, "temperature_k", 1)-300.0))).abs() < 1e-12);
    assert!((result.f64_field("total_joule_w").unwrap() - 1.489_377_453_809_321).abs() < 1e-9);
    assert_ne!(initial.as_array().unwrap()[0].str_field("state_identity"),
        interface(&result).str_field("state_identity"));
    assert!(result.f64_field("passes").unwrap() > 2.0);
}

#[test]
fn g0_both_cell_temperatures_and_pressure_must_remain_in_interface_support() {
    for (text, upper_b) in [
        (SCENE.to_owned(), 300.2),
        (SCENE.replace("100000 Pa", "99999 Pa"), 350.0),
    ] {
        let path = fixture(&text, &joint(2.0, 0.0, upper_b, false, RESISTANCE_DIMS));
        let output = run(&[path.as_os_str().to_owned()], true);
        assert_eq!(output.exit_code, exit::REFUSED, "{}", output.stderr);
        assert!(output.stdout.is_empty());
        let diagnostic = J::parse(&output.stderr).unwrap();
        assert!(diagnostic.str_field("message").unwrap().contains("interface.electrical_resistance"));
    }
}

#[test]
fn g0_order_units_duplicate_faces_and_temperature_conventions_are_not_inferred() {
    let good = joint(2.0, 0.0, 350.0, false, RESISTANCE_DIMS);
    for text in [
        SCENE.replace("\"faces\":[[1,2]]", "\"faces\":[[2,1]]"),
        SCENE.replace("\"faces\":[[1,2]]", "\"faces\":[[1,2],[2,1]]"),
        SCENE.replace("\"faces\":[[1,2]]", "\"faces\":[[0,2]]"),
        SCENE.replace("\"heat_fraction_to_a\":0.25", "\"heat_fraction_to_a\":1.25"),
        SCENE.replace("\"temperature_axes\":[\"Ta\",\"Tb\"],\"temperature_kind\":\"legacy\"",
            "\"temperature_axes\":[\"Ta\",\"Tb\"],\"temperature_kind\":\"temperature-difference\""),
    ] {
        let path = fixture(&text, &good);
        let output = run(&[path.as_os_str().to_owned()], true);
        assert_eq!(output.exit_code, exit::INPUT, "{}", output.stderr);
        assert!(output.stdout.is_empty());
    }
    for (bytes, expected_exit) in [
        (joint(2.0, 0.0, 350.0, true, RESISTANCE_DIMS), exit::INPUT),
        (joint(2.0, 0.0, 350.0, false, ELECTRICAL_RESISTIVITY_DIMS), exit::REFUSED),
    ] {
        let path = fixture(SCENE, &bytes);
        let output = run(&[path.as_os_str().to_owned()], true);
        assert_eq!(output.exit_code, expected_exit, "{}", output.stderr);
        assert!(output.stdout.is_empty());
    }
}

#[test]
fn g1_transient_joint_heat_is_stored_once_and_survives_drive_switch_off() {
    let text = SCENE.replace("{\"temperature_k\":300.0}", "\"adiabatic\"");
    let path = fixture(&text, &joint(2.0, 0.0, 350.0, false, RESISTANCE_DIMS));
    let time = path.parent().unwrap().join("time.json");
    std::fs::write(&time, r#"{
        "schema":"frankensim.electrothermal.transient.v1","time_step_s":0.1,"steps":8,
        "temperature_min_k":290.0,"temperature_max_k":350.0,
        "enthalpy_interpolation_error_j_kg":0.001,"energy_absolute_tolerance_j":1e-8,
        "energy_relative_tolerance":1e-9,"schedule":[[1,1,1],[5,0,1]]
    }"#).unwrap();
    let args = [path.as_os_str().to_owned(), OsString::from("--transient"), time.as_os_str().to_owned()];
    let output = run(&args, true);
    assert_eq!(output.exit_code, exit::SUCCESS, "{}", output.stderr);
    let result = J::parse(&output.stdout).unwrap();
    let records = result.get("records").unwrap().as_array().unwrap();
    assert_eq!(records.len(), 8);
    let mut stored = 0.0;
    for (step, record) in records.iter().enumerate() {
        let power_j = if step < 4 {0.15} else {0.0}; // V²/(R_bulk+R_contact) * dt
        assert!((record.f64_field("source_j").unwrap() - power_j).abs() < 1e-10);
        assert!(record.f64_field("closure_j").unwrap().abs() < 1e-8);
        assert!(record.f64_field("physical_outflow_j").unwrap().abs() < 1e-10);
        stored += record.f64_field("stored_energy_change_j").unwrap();
    }
    assert!((stored - 0.6).abs() < 1e-8);
    assert_eq!(result.f64_field("total_joule_w"), Some(0.0));
    let expected = [300.093_712_110_723_74, 300.185_077_438_680_4,
        300.208_554_453_966_5, 300.112_655_996_629_37];
    for (cell, t) in expected.into_iter().enumerate() {
        assert!((component(&result, "temperature_k", cell) - t).abs() < 1e-7);
    }
    assert_eq!(interface(&result).f64_field("resistance_ohm_m2"), Some(2.0));
    assert_eq!(output, run(&args, true));
}
