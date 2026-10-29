//! Actual file/pack/query/solver tests. All material values are synthetic.
use super::*;
use fs_blake3::hash_bytes;
use fs_conduction::interface::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as CONTACT_DIMS,
    AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as CONTACT_PROPERTY,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, NormalizedInterfacePack, NormalizedPack,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, SurfaceSpec, SystemContext,
    UncertaintyModel,
};
use fs_material::conductor::{ELECTRICAL_RESISTIVITY_DIMS, ELECTRICAL_RESISTIVITY_PROPERTY};
use fs_material::electrical_interface::{ELECTRICAL_CONTACT_RESISTANCE_DIMS, ELECTRICAL_CONTACT_RESISTANCE_PROPERTY};
use fs_material::state_point::{
    DENSITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS, SPECIFIC_HEAT_CAPACITY_PROPERTY,
    THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY,
};
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

const SCENE: &str = r#"{
  "schema":"frankensim.electrothermal.v1","cells":[3,1,1],"cell_size_m":1.0,
  "materials":[{"pack":"bulk.fsmcdpk","conducting":false,
    "temperature_axis":"T","temperature_kind":"legacy",
    "state":[{"name":"pressure","value":"101325 Pa","kind":"legacy"}]}],
  "material_indices":[0,0,0],"initial_temperature_k":300.0,
  "thermal_faces":[{"heat_flux_w_m2":4.0},{"temperature_k":300.0},
    "adiabatic","adiabatic","adiabatic","adiabatic"],
  "electrodes":[{"lo":[0,0,0],"hi":[1,1,1],"voltage_v":0.0}],
  "thermal_interfaces":[{"pack":"joint.fsintpk","faces":[[0,1]],
    "temperature_axes":["Ta","Tb"],"temperature_kind":"legacy",
    "state":[{"name":"normal_pressure","value":"100000 Pa","kind":"legacy"}]}],
  "numerics":{"max_passes":100,"max_thermal_passes":100,
    "temperature_tolerance_k":1e-9,"power_tolerance_w":1e-9,"relative_tolerance":1e-10,
    "conductivity_relative_tolerance":1e-10,"thermal_contact_relative_tolerance":1e-10,
    "energy_tolerance":1e-12,"energy_max_iterations":2000}
}"#;
const TIME: &str = r#"{
  "schema":"frankensim.electrothermal.transient.v1","time_step_s":0.1,"steps":8,
  "temperature_min_k":290.0,"temperature_max_k":340.0,
  "enthalpy_interpolation_error_j_kg":1e-4,"energy_absolute_tolerance_j":1e-8,
  "energy_relative_tolerance":1e-9,"schedule":[[1,0.0,1.0],[5,0.0,0.0]]
}"#;

fn material_id() -> MaterialStateId {
    MaterialStateId { chemistry: "synthetic thermal-joint solid".into(), phase: "solid".into(),
        process: "authored analytic fixture".into(), revision: 0 }
}
fn provenance() -> Provenance {
    Provenance { source: "synthetic analytic thermal-joint fixture, not measured material".into(),
        license: "CC0: authored test data".into(), artifact: None }
}
fn pack(name: &str, mut claims: ClaimSet) -> NormalizedPack {
    super::tests::authored_observation(&mut claims, "synthetic analytic thermal-joint fixture, not measured material");
    NormalizedPack::new(name, "thermal-joint-test-v1", hash_bytes(b"synthetic fixture"),
        "CC0: authored test data", claims, Vec::new(), Vec::new()).unwrap()
}
fn bulk(conducting: bool) -> Vec<u8> {
    let mut properties = vec![
        (THERMAL_CONDUCTIVITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS, 1.0),
        (DENSITY_PROPERTY, fs_qty::Density::DIMS, 1.0),
        (SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS, 2.0),
    ];
    if conducting { properties.push((ELECTRICAL_RESISTIVITY_PROPERTY, ELECTRICAL_RESISTIVITY_DIMS, 2.0)); }
    let mut claims = ClaimSet::new();
    for (name, dims, value) in properties {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims), value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 290.0, 400.0)
                .with("pressure", 101_325.0, 101_325.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: Vec::new(), provenance: provenance(),
        }).unwrap();
    }
    NormalizedMaterialCardPack::new(material_id(), pack("thermal-joint-bulk", claims)).unwrap().to_bytes()
}
fn joint(r: f64, slope: f64, upper: f64, dims: Dims, electrical: bool) -> Vec<u8> {
    let validity = ValidityDomain::unconstrained().with("Ta", 290.0, upper)
        .with("Tb", 290.0, upper).with("normal_pressure", 100_000.0, 100_000.0);
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(CONTACT_PROPERTY, dims),
        value: PropertyValue::Curve { abscissa: "Ta".into(), abscissa_dims: Temperature::DIMS,
            knots: vec![(290.0, r*(1.0-10.0*slope)), (upper, r*(1.0+slope*(upper-300.0)))], dims },
        validity: validity.clone(), uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: Vec::new(), provenance: provenance(),
    }).unwrap();
    if electrical {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(ELECTRICAL_CONTACT_RESISTANCE_PROPERTY, ELECTRICAL_CONTACT_RESISTANCE_DIMS),
            value: PropertyValue::Scalar { value: 2.0, dims: ELECTRICAL_CONTACT_RESISTANCE_DIMS },
            validity, uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: Vec::new(), provenance: provenance(),
        }).unwrap();
    }
    NormalizedInterfacePack::new(
        SurfaceSpec { material: material_id(), texture_frame: "authored side A".into() },
        SurfaceSpec { material: material_id(), texture_frame: "authored side B".into() },
        SystemContext { medium: "synthetic zero-thickness joint".into(), third_body: None,
            environment: "analytic fixture".into(), history: "fixed authored state".into() },
        pack("thermal-joint-interface", claims),
    ).unwrap().to_bytes()
}
fn fixture(scene: &str, bulk: &[u8], joint: &[u8]) -> PathBuf {
    static NEXT: AtomicU64 = AtomicU64::new(0);
    let timestamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let dir = std::env::temp_dir().join(format!("frankensim-thermal-joint-{}-{timestamp}-{}",
        std::process::id(), NEXT.fetch_add(1, Ordering::Relaxed)));
    std::fs::create_dir(&dir).unwrap();
    std::fs::write(dir.join("bulk.fsmcdpk"), bulk).unwrap();
    std::fs::write(dir.join("joint.fsintpk"), joint).unwrap();
    std::fs::write(dir.join("time.json"), TIME).unwrap();
    let path = dir.join("scene.json");
    std::fs::write(&path, scene).unwrap();
    path
}
fn success(args: &[OsString]) -> (CommandOutput, J) {
    let output = run(args, true);
    assert_eq!(output.exit_code, exit::SUCCESS, "{}", output.stderr);
    assert!(output.stderr.is_empty());
    let value = J::parse(&output.stdout).unwrap();
    (output, value)
}
fn temperature(value: &J, cell: usize) -> f64 {
    value.get("temperature_k").unwrap().as_array().unwrap()[cell].as_f64().unwrap()
}
fn first_joint(value: &J) -> &J {
    &value.get("thermal_interfaces").unwrap().as_array().unwrap()[0]
}

#[test]
fn g1_g3_thermal_only_interface_swaps_change_heat_transfer_without_electrical_claims() {
    let material = bulk(false);
    let mut identities = Vec::new();
    for (r, peak) in [(1.0, 314.0), (3.0, 322.0)] {
        let interface = joint(r, 0.0, 350.0, CONTACT_DIMS, false);
        let path = fixture(SCENE, &material, &interface);
        let args = [path.as_os_str().to_owned()];
        let (output, value) = success(&args);
        for (cell, expected) in [peak, 306.0, 302.0].into_iter().enumerate() {
            assert!((temperature(&value, cell)-expected).abs() < 1e-8);
        }
        assert_eq!(value.f64_field("total_joule_w"), Some(0.0));
        assert!(value.f64_field("thermal_balance_relative_residual").unwrap() < 1e-9);
        assert_eq!(first_joint(&value).f64_field("resistance_m2_k_w"), Some(r));
        identities.push(first_joint(&value).str_field("pack_identity").unwrap().to_owned());
        assert_eq!(output, run(&args, true));
        assert_eq!(std::fs::read(&path).unwrap(), SCENE.as_bytes());
        assert_eq!(std::fs::read(path.parent().unwrap().join("joint.fsintpk")).unwrap(), interface);
        assert_eq!(std::fs::read(path.parent().unwrap().join("bulk.fsmcdpk")).unwrap(), material);
    }
    assert_ne!(identities[0], identities[1]);
}

#[test]
fn g1_nonlinear_joint_is_resolved_at_the_accepted_temperature_not_the_initial_state() {
    let path = fixture(SCENE, &bulk(false), &joint(1.0, 0.02, 350.0, CONTACT_DIMS, false));
    let (_, value) = success(&[path.as_os_str().to_owned()]);
    // d = 4*(2.5 + 1 + 0.02*d), so d = 14/0.92. Other cells stay unchanged.
    assert!((temperature(&value, 0) - (300.0 + 14.0/0.92)).abs() < 1e-8);
    assert!((temperature(&value, 1) - 306.0).abs() < 1e-8);
    let contact = first_joint(&value);
    let resolved = contact.f64_field("resistance_m2_k_w").unwrap();
    let used = contact.f64_field("solved_resistance_m2_k_w").unwrap();
    assert!((resolved - (1.0 + 0.02*(temperature(&value, 0)-300.0))).abs() < 1e-12);
    assert!((resolved-used).abs()/resolved < 1e-10);
    assert!(value.f64_field("thermal_contact_relative_defect").unwrap() <= 1e-10);
    assert!(contact.str_field("source_receipt_identity").is_some());
}

#[test]
fn g1_electrical_and_thermal_responses_share_one_face_without_double_counting() {
    let electrical = r#""electrical_interfaces":[{"pack":"joint.fsintpk","faces":[[1,2]],
      "temperature_axes":["Ta","Tb"],"temperature_kind":"legacy","heat_fraction_to_a":0.25,
      "state":[{"name":"normal_pressure","value":"100000 Pa","kind":"legacy"}]}],"numerics":"#;
    let scene = SCENE.replace("\"cells\":[3,1,1]", "\"cells\":[4,1,1]")
        .replace("\"material_indices\":[0,0,0]", "\"material_indices\":[0,0,0,0]")
        .replace("\"conducting\":false", "\"conducting\":true")
        .replace("\"faces\":[[0,1]]", "\"faces\":[[1,2]]")
        .replace("{\"heat_flux_w_m2\":4.0}", "{\"temperature_k\":300.0}")
        .replace("\"voltage_v\":0.0}]", "\"voltage_v\":0.0},{\"lo\":[3,0,0],\"hi\":[4,1,1],\"voltage_v\":3.0}]")
        .replace("\"numerics\":", electrical);
    let path = fixture(&scene, &bulk(true), &joint(2.0, 0.0, 350.0, CONTACT_DIMS, true));
    let (_, value) = success(&[path.as_os_str().to_owned()]);
    // Two 2-ohm bulk cells plus a 2-ohm contact carry 0.5 A. Electrical contact
    // heat is 0.5 W, split 0.125/0.375 W; thermal resistance is independently 2.
    for (cell, expected) in [300.34375, 301.03125, 301.21875, 300.40625].into_iter().enumerate() {
        assert!((temperature(&value, cell)-expected).abs() < 1e-8);
    }
    assert!((value.f64_field("total_joule_w").unwrap()-1.5).abs() < 1e-9);
    assert!((value.f64_field("electrode_power_w").unwrap()-1.5).abs() < 1e-9);
    let electric = &value.get("electrical_interfaces").unwrap().as_array().unwrap()[0];
    assert_eq!(first_joint(&value).get("cells"), electric.get("cells"));
    assert!(value.f64_field("thermal_balance_relative_residual").unwrap() < 1e-9);
}

#[test]
fn g1_g3_transient_joint_swap_changes_redistribution_not_total_stored_energy() {
    let scene = SCENE.replace("{\"heat_flux_w_m2\":4.0},{\"temperature_k\":300.0}", "\"adiabatic\",\"adiabatic\"")
        .replace("\"numerics\":", "\"power_w\":[2.0,0.0,0.0],\"numerics\":");
    let mut peaks = Vec::new();
    for r in [1.0, 3.0] {
        let path = fixture(&scene, &bulk(false), &joint(r, 0.0, 350.0, CONTACT_DIMS, false));
        let time = path.parent().unwrap().join("time.json");
        let args = [path.as_os_str().to_owned(), OsString::from("--transient"), time.into_os_string()];
        let (output, value) = success(&args);
        let stored: f64 = (0..3).map(|c| 2.0*(temperature(&value,c)-300.0)).sum();
        assert!((stored-0.8).abs() < 1e-8);
        let records = value.get("records").unwrap().as_array().unwrap();
        assert_eq!(records.len(), 8);
        for (step, record) in records.iter().enumerate() {
            assert!(record.f64_field("closure_j").unwrap().abs() < 1e-8);
            assert!((record.f64_field("source_j").unwrap() - if step < 4 { 0.2 } else { 0.0 }).abs() < 1e-12);
        }
        assert_eq!(first_joint(&value).f64_field("resistance_m2_k_w"), Some(r));
        assert_eq!(output, run(&args, true));
        peaks.push(temperature(&value,0));
    }
    assert!(peaks[1] > peaks[0]);
}

#[test]
fn g0_missing_controls_invalid_units_duplicate_faces_and_unsupported_final_states_refuse() {
    for text in [
        SCENE.replace(",\"thermal_contact_relative_tolerance\":1e-10", ""),
        SCENE.replace("\"thermal_contact_relative_tolerance\":1e-10", "\"thermal_contact_relative_tolerance\":1.0"),
        SCENE.replace("\"faces\":[[0,1]]", "\"faces\":[[0,1],[1,0]]"),
        SCENE.replace("\"temperature_axes\":[\"Ta\",\"Tb\"]", "\"temperature_axes\":[\"Ta\",\"Ta\"]"),
    ] {
        let path = fixture(&text, &bulk(false), &joint(1.0,0.0,350.0,CONTACT_DIMS,false));
        let output = run(&[path.into_os_string()],true);
        assert_eq!(output.exit_code,exit::INPUT,"{}",output.stderr);
        assert!(output.stdout.is_empty());
    }
    for (upper, dims) in [(310.0, CONTACT_DIMS), (350.0, ELECTRICAL_CONTACT_RESISTANCE_DIMS)] {
        let path = fixture(SCENE, &bulk(false), &joint(1.0,0.0,upper,dims,false));
        let output = run(&[path.into_os_string()],true);
        assert_eq!(output.exit_code,exit::REFUSED,"{}",output.stderr);
        assert!(output.stdout.is_empty());
        assert!(output.stderr.contains("interface.thermal_resistance"));
    }
}
