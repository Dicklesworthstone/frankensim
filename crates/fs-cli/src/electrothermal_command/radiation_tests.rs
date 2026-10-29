//! Actual file -> canonical card -> spatial electrical/thermal/radiation path.
//! Synthetic coefficients only; no measured-material qualification is asserted.
use super::*;
use fs_blake3::hash_bytes;
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialStateId, NormalizedPack,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, UncertaintyModel};
use fs_conduction::radiation::{EMISSIVITY_DIMS, SURFACE_EMISSIVITY_PROPERTY, STEFAN_BOLTZMANN_W_M2_K4 as SIGMA};
use fs_material::conductor::{ELECTRICAL_RESISTIVITY_DIMS, ELECTRICAL_RESISTIVITY_PROPERTY};
use fs_material::state_point::{DENSITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS,
    SPECIFIC_HEAT_CAPACITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY};
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

const SCENE: &str = r#"{
 "schema":"frankensim.electrothermal.v1","cells":[3,1,1],"cell_size_m":1,
 "materials":[{"pack":"material.fsmcdpk","conducting":true,"temperature_axis":"T",
 "temperature_kind":"legacy","state":[{"name":"pressure","value":"101325 Pa","kind":"legacy"}]}],
 "material_indices":[0,0,0],"initial_temperature_k":300,
 "thermal_faces":["adiabatic","adiabatic","adiabatic","adiabatic","adiabatic","adiabatic"],
 "electrodes":[{"lo":[0,0,0],"hi":[1,1,1],"voltage_v":0},
 {"lo":[2,0,0],"hi":[3,1,1],"voltage_v":30}],
 "numerics":{"max_passes":100,"max_thermal_passes":250,"temperature_tolerance_k":1e-9,
 "power_tolerance_w":1e-9,"relative_tolerance":1e-10,"conductivity_relative_tolerance":1e-10,
 "energy_tolerance":1e-12,"energy_max_iterations":2000},
 "radiation":{"max_iterations":250,"temperature_tolerance_k":1e-9,"power_tolerance_w":1e-9,
 "relative_tolerance":1e-10,"relaxation":0.5,"surfaces":[{"name":"exterior",
 "material_index":0,"faces":[[1,3]],"surroundings_k":300,"initial_surface_temperature_k":450}]}
}"#;

fn material(eps: f64, slope: f64, maximum_surface: f64) -> Vec<u8> {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (ELECTRICAL_RESISTIVITY_PROPERTY, ELECTRICAL_RESISTIVITY_DIMS, 2.0),
        (THERMAL_CONDUCTIVITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS, 1.0),
        (DENSITY_PROPERTY, fs_qty::Density::DIMS, 1.0),
        (SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS, 4.0),
        (SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS, eps),
    ] {
        let optical = name == SURFACE_EMISSIVITY_PROPERTY;
        let maximum = if optical { maximum_surface } else { 1200.0 };
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: if optical { PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                knots: vec![(250.0, value - 50.0*slope), (maximum, value + (maximum-300.0)*slope)], dims } }
                else { PropertyValue::Scalar { value, dims } },
            validity: ValidityDomain::unconstrained().with("T", 250.0, maximum).with("pressure", 101_325.0, 101_325.0),
            interpolation: if optical { InterpolationPolicy::LinearInside } else { InterpolationPolicy::ConstantWithinValidity },
            uncertainty: UncertaintyModel::Unstated, observations: Vec::new(),
            provenance: Provenance { source: "authored synthetic radiation fixture".into(),
                license: "CC0: authored test data".into(), artifact: None },
        }).unwrap();
    }
    super::tests::authored_observation(&mut claims, "authored synthetic radiation fixture");
    let pack = NormalizedPack::new("radiation-fixture", "test-v1", hash_bytes(b"authored radiation fixture"),
        "CC0: authored test data", claims, Vec::new(), Vec::new()).unwrap();
    NormalizedMaterialCardPack::new(MaterialStateId { chemistry: "synthetic radiating conductor".into(),
        phase: "solid".into(), process: "analytic fixture".into(), revision: 0 }, pack).unwrap().to_bytes()
}
fn fixture(bytes: &[u8], text: &str) -> PathBuf {
    static NEXT: AtomicU64 = AtomicU64::new(0);
    let stamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let directory = std::env::temp_dir().join(format!("frankensim-radiation-{}-{stamp}-{}",
        std::process::id(), NEXT.fetch_add(1, Ordering::Relaxed)));
    std::fs::create_dir(&directory).unwrap();
    std::fs::write(directory.join("material.fsmcdpk"), bytes).unwrap();
    let path = directory.join("scene.json"); std::fs::write(&path, text).unwrap(); path
}
fn execute(path: &Path) -> J {
    let output = run(&[path.as_os_str().to_owned()], true);
    assert_eq!(output.exit_code, exit::SUCCESS, "{}", output.stderr);
    assert_eq!(output, run(&[path.as_os_str().to_owned()], true));
    J::parse(&output.stdout).unwrap()
}
fn surface(result: &J) -> &J { &result.get("radiation").unwrap().as_array().unwrap()[0] }
fn cell(result: &J) -> f64 { result.get("temperature_k").unwrap().as_array().unwrap()[1].as_f64().unwrap() }

#[test]
fn g1_g3_emissivity_only_swaps_change_surface_and_bulk_temperature_not_authored_power() {
    for eps in [0.2, 0.8] {
        let bytes = material(eps, 0.0, 1200.0);
        let path = fixture(&bytes, SCENE);
        let result = execute(&path);
        let expected = (300.0_f64.powi(4) + 450.0 / (eps * SIGMA)).sqrt().sqrt();
        assert!((surface(&result).f64_field("surface_temperature_k").unwrap() - expected).abs() < 1e-7);
        assert!((cell(&result) - expected - 225.0).abs() < 1e-7);
        assert!((result.f64_field("total_joule_w").unwrap() - 450.0).abs() < 1e-8);
        assert!((surface(&result).f64_field("radiated_w").unwrap() - 450.0).abs() < 1e-7);
        assert!(surface(&result).f64_field("heat_defect_w").unwrap() < 1e-7);
        assert_eq!(std::fs::read(path.parent().unwrap().join("material.fsmcdpk")).unwrap(), bytes);
    }
}

#[test]
fn g1_temperature_dependent_emissivity_and_convection_use_one_surface_and_half_cell() {
    let text = SCENE.replace("\"surroundings_k\":300", "\"convection\":{\"h_w_m2_k\":2,\"ambient_k\":280},\"surroundings_k\":300");
    let path = fixture(&material(0.2, 0.0005, 1200.0), &text);
    let result = execute(&path);
    let face = surface(&result);
    let t = face.f64_field("surface_temperature_k").unwrap();
    assert!((t - 393.6478510454975).abs() < 1e-7);
    assert!((cell(&result) - t - 225.0).abs() < 1e-7);
    assert!((face.f64_field("emissivity").unwrap() - (0.2 + 0.0005*(t-300.0))).abs() < 1e-12);
    assert!((face.f64_field("radiated_w").unwrap() + face.f64_field("convective_w").unwrap() - 450.0).abs() < 1e-7);
}

#[test]
fn g0_emissivity_queries_use_supported_surface_not_hotter_cell_and_refuse_true_domain_exit() {
    let path = fixture(&material(0.2, 0.0, 500.0), SCENE);
    let result = execute(&path);
    assert!(cell(&result) > 500.0);
    assert!(surface(&result).f64_field("surface_temperature_k").unwrap() < 500.0);
    let text = SCENE.replace("\"initial_surface_temperature_k\":450", "\"initial_surface_temperature_k\":300");
    let path = fixture(&material(0.2, 0.0, 300.5), &text);
    let output = run(&[path.as_os_str().to_owned()], true);
    assert_eq!(output.exit_code, exit::REFUSED);
    assert!(output.stdout.is_empty());
    assert!(output.stderr.contains("radiation.material"));
    for text in [SCENE.replace("[[1,3]]", "[[1,3],[1,3]]"), SCENE.replace("[[1,3]]", "[[1,1]]"),
        SCENE.replace("\"relaxation\":0.5", "\"relaxation\":0")] {
        let path = fixture(&material(0.2, 0.0, 1200.0), &text);
        let output = run(&[path.as_os_str().to_owned()], true);
        assert_eq!(output.exit_code, exit::INPUT); assert!(output.stdout.is_empty());
    }
}

#[test]
fn g1_transient_radiative_cooling_equals_lost_caloric_energy() {
    let text = SCENE.replace("\"initial_temperature_k\":300", "\"initial_temperature_k\":400");
    let mut outgoing = Vec::new();
    for eps in [0.2, 0.8] {
        let path = fixture(&material(eps, 0.0, 1200.0), &text);
        let time = path.parent().unwrap().join("time.json");
        std::fs::write(&time, r#"{"schema":"frankensim.electrothermal.transient.v1",
            "time_step_s":0.1,"steps":5,"temperature_min_k":250,"temperature_max_k":1200,
            "enthalpy_interpolation_error_j_kg":1e-8,"energy_absolute_tolerance_j":1e-7,
            "energy_relative_tolerance":1e-9,"schedule":[[1,0,0]]}"#).unwrap();
        let output = run(&[path.as_os_str().to_owned(), "--transient".into(), time.into_os_string()], true);
        assert_eq!(output.exit_code, exit::SUCCESS, "{}", output.stderr);
        let result = J::parse(&output.stdout).unwrap();
        let storage: f64 = result.get("temperature_k").unwrap().as_array().unwrap().iter()
            .map(|t| 4.0 * (t.as_f64().unwrap() - 400.0)).sum();
        let records = result.get("records").unwrap().as_array().unwrap();
        let heat: f64 = records.iter().map(|r| r.f64_field("physical_outflow_j").unwrap()).sum();
        assert!((storage + heat).abs() < 1e-6); assert!(heat > 0.0);
        for record in records {
            assert!(record.f64_field("closure_j").unwrap().abs() < 1e-7);
            assert_eq!(record.f64_field("source_j").unwrap(), 0.0);
        }
        outgoing.push(heat);
    }
    assert!(outgoing[1] > outgoing[0]);
}
