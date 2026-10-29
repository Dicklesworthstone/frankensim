//! File-to-solver transient tests. All material values are authored analytic
//! fixtures, not experimental metal data. No test replaces the physical owners.

use super::*;
use fs_blake3::hash_bytes;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, NormalizedPack, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, UncertaintyModel,
};
use fs_material::conductor::{ELECTRICAL_RESISTIVITY_DIMS, ELECTRICAL_RESISTIVITY_PROPERTY};
use fs_material::state_point::{
    DENSITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS, SPECIFIC_HEAT_CAPACITY_PROPERTY,
    THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY,
};
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

fn pack(density: f64, nonlinear: bool, include_capacity: bool) -> Vec<u8> {
    let mut claims = ClaimSet::new();
    for (name, dims, value, slope) in [
        (ELECTRICAL_RESISTIVITY_PROPERTY, ELECTRICAL_RESISTIVITY_DIMS, 2.0, if nonlinear { 0.03 } else { 0.0 }),
        (THERMAL_CONDUCTIVITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS, 1.0, if nonlinear { 0.02 } else { 0.0 }),
        (DENSITY_PROPERTY, fs_qty::Density::DIMS, density, if nonlinear { 0.01 } else { 0.0 }),
        (SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS, 3.0, if nonlinear { 0.015 } else { 0.0 }),
    ] {
        if name == SPECIFIC_HEAT_CAPACITY_PROPERTY && !include_capacity { continue; }
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Curve {
                abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                knots: vec![(290.0, value * (1.0 - 10.0 * slope)), (350.0, value * (1.0 + 50.0 * slope))], dims,
            },
            validity: ValidityDomain::unconstrained().with("T", 290.0, 350.0)
                .with("pressure", 101325.0, 101325.0),
            uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
            observations: Vec::new(), provenance: Provenance {
                source: "authored transient resistor-network fixture".into(),
                license: "CC0: authored test data".into(), artifact: None,
            },
        }).unwrap();
    }
    super::tests::authored_observation(&mut claims, "authored transient resistor-network fixture");
    let claims = NormalizedPack::new("synthetic-transient", "test-v1", hash_bytes(b"authored thermal storage"),
        "CC0: authored test data", claims, Vec::new(), Vec::new()).unwrap();
    NormalizedMaterialCardPack::new(MaterialStateId {
        chemistry: "synthetic".into(), phase: "solid".into(), process: "analytic".into(), revision: 0,
    }, claims).unwrap().to_bytes()
}

fn scene() -> &'static str {
    r#"{"schema":"frankensim.electrothermal.v1","cells":[3,1,1],"cell_size_m":1.0,
      "materials":[{"pack":"material.fsmcdpk","conducting":true,"temperature_axis":"T","temperature_kind":"legacy",
        "state":[{"name":"pressure","value":"101325 Pa","kind":"legacy"}]}],
      "material_indices":[0,0,0],"initial_temperature_k":300.0,
      "thermal_faces":["adiabatic","adiabatic","adiabatic","adiabatic","adiabatic","adiabatic"],
      "electrodes":[{"lo":[0,0,0],"hi":[1,1,1],"voltage_v":0.0},{"lo":[2,0,0],"hi":[3,1,1],"voltage_v":3.0}],
      "numerics":{"max_passes":100,"max_thermal_passes":100,"temperature_tolerance_k":1e-9,"power_tolerance_w":1e-9,
        "relative_tolerance":1e-10,"conductivity_relative_tolerance":1e-10,"energy_tolerance":1e-12,"energy_max_iterations":2000}}"#
}

fn request() -> String {
    r#"{"schema":"frankensim.electrothermal.transient.v1","time_step_s":0.1,"steps":8,
      "temperature_min_k":290.0,"temperature_max_k":350.0,"enthalpy_interpolation_error_j_kg":0.00001,
      "energy_absolute_tolerance_j":0.0000001,"energy_relative_tolerance":0.000000001,
      "schedule":[[1,1.0,1.0],[5,0.0,1.0]]}"#.into()
}

fn fixture(bytes: &[u8], request: &str) -> (PathBuf, PathBuf) {
    static NEXT: AtomicU64 = AtomicU64::new(0);
    let id = NEXT.fetch_add(1, Ordering::Relaxed);
    let timestamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let dir = std::env::temp_dir().join(format!("frankensim-thermal-time-{}-{timestamp}-{id}", std::process::id()));
    std::fs::create_dir(&dir).unwrap();
    std::fs::write(dir.join("material.fsmcdpk"), bytes).unwrap();
    let (scene_path, time_path) = (dir.join("scene.json"), dir.join("time.json"));
    std::fs::write(&scene_path, scene()).unwrap();
    std::fs::write(&time_path, request).unwrap();
    (scene_path, time_path)
}

fn execute(paths: &(PathBuf, PathBuf)) -> CommandOutput {
    run(&[paths.0.as_os_str().to_owned(), "--transient".into(), paths.1.as_os_str().to_owned()], true)
}
fn numbers(result: &J, key: &str) -> Vec<f64> {
    result.get(key).unwrap().as_array().unwrap().iter().map(|n| n.as_f64().unwrap()).collect()
}

#[test]
fn g1_g3_transient_file_path_preserves_mass_and_energy_across_drive_switch_off() {
    let mut peaks = Vec::new();
    for mass in [2.0, 4.0] {
        let bytes = pack(mass, false, true);
        let text = request();
        let paths = fixture(&bytes, &text);
        let output = execute(&paths);
        assert_eq!(output.exit_code, exit::SUCCESS, "{}", output.stderr);
        let result = J::parse(&output.stdout).unwrap();
        assert_eq!(numbers(&result, "cell_mass_kg"), vec![mass; 3]);
        let temperatures = numbers(&result, "temperature_k");
        let stored: f64 = temperatures.iter().map(|t| mass * 3.0 * (t - 300.0)).sum();
        // Four powered 0.1-s steps at 3 V across a 2-ohm unit conductor.
        assert!((stored - 1.8).abs() < 1e-6);
        assert_eq!(result.f64_field("total_joule_w"), Some(0.0));
        let records = result.get("records").unwrap().as_array().unwrap();
        assert_eq!(records.len(), 8);
        assert!(records[7].f64_field("max_temperature_k").unwrap() < records[3].f64_field("max_temperature_k").unwrap());
        for row in records {
            assert!(row.f64_field("closure_j").unwrap().abs() < 1e-7);
            assert_eq!(row.f64_field("physical_outflow_j"), Some(0.0));
        }
        assert_eq!(output, execute(&paths));
        assert_eq!(std::fs::read_to_string(&paths.0).unwrap(), scene());
        assert_eq!(std::fs::read_to_string(&paths.1).unwrap(), text);
        assert_eq!(std::fs::read(paths.0.parent().unwrap().join("material.fsmcdpk")).unwrap(), bytes);
        peaks.push(temperatures[1]);
    }
    assert!(peaks[1] < peaks[0]);
}

#[test]
fn g1_variable_capacity_uses_integrated_enthalpy_not_frozen_startup_cp() {
    let paths = fixture(&pack(2.0, true, true), &request());
    let output = execute(&paths);
    assert_eq!(output.exit_code, exit::SUCCESS, "{}", output.stderr);
    let result = J::parse(&output.stdout).unwrap();
    let t = numbers(&result, "temperature_k");
    let stored: f64 = t.iter().map(|t| {
        let d = t - 300.0;
        2.0 * 3.0 * (d + 0.5 * 0.015 * d * d)
    }).sum();
    let source: f64 = result.get("records").unwrap().as_array().unwrap().iter()
        .map(|r| r.f64_field("source_j").unwrap()).sum();
    // Continuous source Cp integral differs from the admitted chord chart by
    // at most the declared interpolation allowance, not by frozen-Cp error.
    // 1425 generated knots stay within the material owner's 4096-knot cap.
    assert!((stored - source).abs() < 2e-4);
    assert_eq!(numbers(&result, "cell_mass_kg"), vec![2.0; 3]);
    let materials = result.get("materials").unwrap().as_array().unwrap();
    assert!(materials[0].str_field("storage_chart_identity").is_some());
    assert!(materials[0].f64_field("storage_receipts").unwrap() >= 4.0);
}

#[test]
fn g0_missing_capacity_and_late_storage_domain_failure_publish_no_trajectory() {
    for (bytes, request) in [
        (pack(2.0, false, false), request()),
        (pack(2.0, false, true), request().replace("\"temperature_max_k\":350.0", "\"temperature_max_k\":300.001")),
    ] {
        let paths = fixture(&bytes, &request);
        let output = execute(&paths);
        assert_eq!(output.exit_code, exit::REFUSED, "{}", output.stderr);
        assert!(output.stdout.is_empty());
        assert!(J::parse(&output.stderr).unwrap().str_field("message").unwrap().contains("material.storage"));
    }
}

#[test]
fn g0_transient_time_and_schedule_admission_precede_execution() {
    for request in [request().replace("\"steps\":8", "\"steps\":0"),
        request().replace("[5,0.0,1.0]", "[1,0.0,1.0]"),
        request().replace("[1,1.0,1.0]", "[2,1.0,1.0]"),
        request().replace("\"time_step_s\":0.1", "\"time_step_s\":0.0")]
    {
        let paths = fixture(&pack(2.0, false, true), &request);
        let output = execute(&paths);
        assert_eq!(output.exit_code, exit::INPUT);
        assert!(output.stdout.is_empty());
    }
}

#[test]
fn g4_cancelled_transient_does_not_modify_the_input_problem() {
    let paths = fixture(&pack(2.0, false, true), &request());
    let problem = load(&paths.0).unwrap();
    let request = source::transient::load(&paths.1, problem.domain.cell_count()).unwrap();
    let before = problem.initial_temperature.clone();
    let electric = problem.electric.clone();
    let gate = CancelGate::new();
    gate.request();
    assert!(matches!(source::transient::march(&problem, &request, &gate), Err(ChtError::Cancelled)));
    assert_eq!(problem.initial_temperature, before);
    assert_eq!(problem.electric, electric);
}
