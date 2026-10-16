//! G0 count/conservation checks and G3 order replay through the real binary.
use super::*;

fn generated(input: &mut J) {
    let enclosure = enclosure(input);
    for key in ["view_factors", "row_sum_tolerance", "reciprocity_relative_tolerance", "evidence"] {
        remove(enclosure, key);
    }
    put(enclosure, "generate_view_factors", J::parse(r#"{
        "seed":"41",
        "samples_per_replicate":512,
        "replicates":4,
        "max_triangles":256,
        "max_ray_triangle_tests":"1000000",
        "projection":{
            "max_iterations":10000,
            "row_sum_tolerance":1e-12,
            "max_abs_adjustment":0.25
        }
    }"#).unwrap());
}

fn generation(input: &mut J) -> &mut J {
    member(enclosure(input), "generate_view_factors")
}

/// Twenty-six conforming cubes surround one omitted central cube.
/// Its twelve inward-facing boundary triangles form a closed radiating cavity;
/// every other boundary triangle is opaque and thermally adiabatic.
fn cavity() -> J {
    let mut input = J::parse(BASE).unwrap();
    let mut solid = J::parse(r#"{
        "conductivity_w_m_k":10,
        "source_w_m3":0,
        "adiabatic_remainder":true,
        "surfaces":[
            {"name":"emitter","htc_w_m2_k":10,
             "faces":[[21,25,41],[21,37,41],[21,22,38],[21,37,38],[21,22,26],[21,25,26]]},
            {"name":"receiver","htc_w_m2_k":25,
             "faces":[[22,26,42],[22,38,42],[25,26,42],[25,41,42],[37,38,42],[37,41,42]]}
        ]
    }"#).unwrap();
    let coordinates = [-1.0, 0.0, 1.0, 2.0];
    let index = |x: usize, y: usize, z: usize| x + 4 * y + 16 * z;
    let mut vertices = Vec::new();
    for z in 0..4 {
        for y in 0..4 {
            for x in 0..4 {
                vertices.push(J::Array(vec![num(coordinates[x]), num(coordinates[y]), num(coordinates[z])]));
            }
        }
    }
    let mut tetrahedra = Vec::new();
    for z in 0..3 {
        for y in 0..3 {
            for x in 0..3 {
                if [x,y,z] == [1,1,1] { continue; }
                let v = [
                    index(x,y,z), index(x+1,y,z), index(x,y+1,z), index(x+1,y+1,z),
                    index(x,y,z+1), index(x+1,y,z+1), index(x,y+1,z+1), index(x+1,y+1,z+1),
                ];
                for local in [[0,1,3,7],[0,1,5,7],[0,2,3,7],[0,2,6,7],[0,4,5,7],[0,4,6,7]] {
                    tetrahedra.push(J::Array(local.iter().map(|&i| num(v[i] as f64)).collect()));
                }
            }
        }
    }
    put(&mut solid, "vertices_m", J::Array(vertices));
    put(&mut solid, "tetrahedra", J::Array(tetrahedra));
    put(&mut input, "solid", solid);
    generated(&mut input);
    input
}

fn matrix(value: &J) -> Vec<Vec<f64>> {
    value.as_array().unwrap().iter().map(|row|
        row.as_array().unwrap().iter().map(|value| value.as_f64().unwrap()).collect()
    ).collect()
}

fn count(value: &J) -> u64 {
    value.as_str().expect("retained u64 counts are decimal strings").parse().unwrap()
}

#[test]
fn g0_generated_cavity_retains_raw_counts_errors_and_projection() {
    let result = run(&cavity());
    assert_eq!(result.str_field("authority"), Some("nominal-estimate"));
    near(n(result.get("objective").unwrap(), "value_k"), 300.0, 1e-5);
    let radiation = result.get("radiation").unwrap();
    let evidence = radiation.get("view_factor_evidence").unwrap();
    assert_eq!(evidence.str_field("kind"), Some("geometry-qmc"));
    assert_eq!(evidence.str_field("seed"), Some("41"));
    assert_eq!(evidence.str_field("count_encoding"), Some("decimal-u64-string"));
    assert_eq!(evidence.str_field("mesh_identity").unwrap().len(), 64);
    assert_eq!(evidence.str_field("estimate_identity").unwrap().len(), 64);
    assert_eq!(evidence.get("surface_names"), Some(&J::parse(r#"["emitter","receiver"]"#).unwrap()));
    let samples = count(evidence.get("samples_per_surface").unwrap());
    assert_eq!(samples, 2048);
    assert!(count(evidence.get("ray_triangle_tests").unwrap()) <= 1_000_000);
    let hits = evidence.get("hit_counts").unwrap().as_array().unwrap();
    let raw = matrix(evidence.get("raw_view_factors").unwrap());
    let errors = matrix(evidence.get("empirical_standard_errors").unwrap());
    let projected = matrix(radiation.get("view_factors").unwrap());
    let areas = evidence.get("areas_m2").unwrap().as_array().unwrap();
    let mut adjustment = 0.0_f64;
    for i in 0..2 {
        near(areas[i].as_f64().unwrap(), 3.0, 1e-14);
        let hits = hits[i].as_array().unwrap();
        let mut accounted = hits.iter().map(count).sum::<u64>();
        for field in ["escaped_counts", "unassigned_counts", "backface_counts"] {
            let losses = count(&evidence.get(field).unwrap().as_array().unwrap()[i]);
            assert_eq!(losses, 0, "{field}");
            accounted += losses;
        }
        assert_eq!(accounted, samples);
        near(projected[i].iter().sum(), 1.0, 2e-12);
        for j in 0..2 {
            near(raw[i][j], count(&hits[j]) as f64 / samples as f64, 1e-15);
            assert!(errors[i][j].is_finite() && errors[i][j] >= 0.0);
            near(areas[i].as_f64().unwrap() * projected[i][j],
                areas[j].as_f64().unwrap() * projected[j][i], 1e-11);
            adjustment = adjustment.max((raw[i][j] - projected[i][j]).abs());
        }
    }
    let projection = evidence.get("projection").unwrap();
    near(n(projection, "max_abs_adjustment"), adjustment, 1e-15);
    assert!(adjustment <= n(projection, "max_abs_adjustment_budget"));
    assert!(n(projection, "iterations") >= 1.0);
    assert!(evidence.str_field("scope").unwrap().contains("zero sampled losses do not prove"));
    assert!(evidence.str_field("scope").unwrap().contains("not calibrated confidence bounds"));
    near(n(radiation, "radiative_out_w"), 0.0, 1e-7);
}

#[test]
fn g3_generated_patch_order_replays_and_seed_changes_identity() {
    let input = cavity();
    let original = output(&input);
    let original_json = success(&original);
    let mut reordered = input.clone();
    rows(member(enclosure(&mut reordered), "surfaces")).reverse();
    let replay = output(&reordered);
    success(&replay);
    assert_eq!(original.stdout, replay.stdout);
    put(generation(&mut reordered), "seed", J::Str("42".into()));
    let reseeded = run(&reordered);
    assert_ne!(
        original_json.path(&["radiation","view_factor_evidence","estimate_identity"]),
        reseeded.path(&["radiation","view_factor_evidence","estimate_identity"]),
    );
}

#[test]
fn g0_generated_open_geometry_refuses_instead_of_normalizing_escape() {
    let mut gap = J::parse(BASE).unwrap();
    generated(&mut gap);
    let result = output(&gap);
    assert!(!result.status.success());
    assert!(result.stdout.is_empty());
    let diagnostic = J::parse(std::str::from_utf8(&result.stderr).unwrap()).unwrap();
    assert!(diagnostic.str_field("message").unwrap().contains("sampled rays leave"));
}

#[test]
fn g0_generated_input_ownership_and_work_budgets_refuse() {
    let input = cavity();
    for (field, value) in [
        ("view_factors", J::parse("[[0,1],[1,0]]").unwrap()),
        ("row_sum_tolerance", num(1e-12)),
        ("reciprocity_relative_tolerance", num(1e-12)),
        ("evidence", J::parse(r#"{"kind":"analytic","geometry":"caller label"}"#).unwrap()),
    ] {
        let mut invalid = input.clone();
        put(enclosure(&mut invalid), field, value);
        let result = output(&invalid);
        assert!(!result.status.success(), "{field}");
        assert!(result.stdout.is_empty());
    }
    for (field, value) in [
        ("seed", num(41.0)),
        ("samples_per_replicate", num(3.0)),
        ("replicates", num(1.0)),
        ("unknown_sampler_knob", num(1.0)),
    ] {
        let mut invalid = input.clone();
        put(generation(&mut invalid), field, value);
        let result = output(&invalid);
        assert!(!result.status.success(), "{field}");
        assert!(result.stdout.is_empty());
    }
    let mut invalid_projection = input.clone();
    put(member(generation(&mut invalid_projection), "projection"), "unknown_projection_knob", num(1.0));
    let result = output(&invalid_projection);
    assert!(!result.status.success());
    assert!(result.stdout.is_empty());
    for (field, value) in [
        ("max_triangles", num(1.0)),
        ("max_ray_triangle_tests", J::Str("1".into())),
    ] {
        let mut exhausted = input.clone();
        put(generation(&mut exhausted), field, value);
        let result = output(&exhausted);
        assert_eq!(result.status.code(), Some(6), "{field}: {}", String::from_utf8_lossy(&result.stderr));
        assert!(result.stdout.is_empty());
    }
}

#[test]
fn g0_generated_heated_example_closes_solid_and_air_energy() {
    let input = J::parse(include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
        "/../../examples/cooling-network/enclosure-radiation-generated.json"))).unwrap();
    let result = run(&input);
    assert_eq!(result.str_field("authority"), Some("nominal-estimate"));
    let radiation = result.get("radiation").unwrap();
    assert_eq!(radiation.path(&["view_factor_evidence", "kind"]).unwrap().as_str(), Some("geometry-qmc"));
    assert!(n(result.get("objective").unwrap(), "value_k") > 300.0);
    near(n(&result, "source_w"), 5.2, 1e-10);
    near(n(&result, "robin_out_w"), 5.2, 1e-7);
    let wall_heat: f64 = result.get("walls").unwrap().as_array().unwrap()
        .iter().map(|wall| n(wall, "outward_heat_w")).sum();
    near(wall_heat, 5.2, 1e-7);
    let air = input.get("air").unwrap();
    let capacity = n(air, "density_kg_m3") * n(air, "specific_heat_j_kg_k");
    let air_heat: f64 = result.get("branches").unwrap().as_array().unwrap()
        .iter().map(|branch| capacity * n(branch, "flow_m3_s")
            * (n(branch, "outlet_k") - n(branch, "inlet_k"))).sum();
    near(air_heat, 5.2, 1e-7);
    near(n(&result, "air_heat_imbalance_w"), 0.0, 1e-7);
    near(n(radiation, "radiative_out_w"), 0.0, 1e-7);
}
