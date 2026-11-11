use super::*;
use fs_viz::vtu::VtuChecker;

const NODES: &str = "[[1,2,3],[2,2,3],[1,4,3],[2,4,3],[1,2,6],[2,2,6],[1,4,6],[2,4,6]]";
const CORNERS: &str = "[[0,1,2,3,4,5,6,7]]";

fn displacement(offset: f64) -> Vec<f64> {
    (0..24).map(|i| offset + f64::from(i) / 1000.0).collect()
}

fn fields(load_space: bool) -> String {
    let density = if load_space {
        "\"raw_density\":0.7,\"eroded_density\":0.4,\"nominal_density\":0.6,\"dilated_density\":0.8"
    } else {
        "\"raw_density\":0.7,\"projected_density\":0.6"
    };
    let displacement = if load_space {
        format!(
            "\"displacement_layout\":\"basis column then node-major xyz, on eroded material\",\"basis_displacements_m\":[{:?},{:?}],\"worst_displacement_m\":{:?}",
            displacement(0.0),
            displacement(0.1),
            displacement(0.2)
        )
    } else {
        format!(
            "\"displacement_layout\":\"load-case then node-major xyz\",\"displacements_m\":[{:?},{:?}],\"cell_relaxed_max_pa\":[[120.5],[80.0]]",
            displacement(0.0),
            displacement(0.1)
        )
    };
    format!(
        "\"cells\":[{{{density},\"cut_volume_m3\":2.5}}],\"physical_nodes_m\":{NODES},\"cell_nodes\":{CORNERS},{displacement}"
    )
}

fn adaptive() -> String {
    format!(
        "{{\"schema\":\"sdf3-study-design.v1\",\"state\":\"solved\",{}}}",
        fields(false)
    )
}

fn incumbent(load_space: bool) -> String {
    let (schema, selection) = if load_space {
        ("sdf3-load-space-design.v1", "least-material-feasible")
    } else {
        ("sdf3-stress-design.v1", "best-feasible")
    };
    // A different last accepted point must never leak into selected output.
    format!(
        "{{\"schema\":{schema:?},\"state\":\"solved\",\"selected_design\":{selection:?},\"selected_feasible\":true,\"selected\":{{{}}},\"last_accepted\":{{{}}}}}",
        fields(load_space),
        fields(load_space).replace("0.7", "0.9")
    )
}

fn affine_strength_design() -> String {
    let selected = format!(
        "{{{},\"offset_displacement_m\":{:?},\"strength\":{{\"schema\":\"sdf3-affine-strength-fields.v1\",\"cell_physical_envelope_pa\":[120.5],\"cell_relaxed_envelope_pa\":[180.5],\"cell_sample_count\":[8]}}}}",
        fields(true),
        displacement(-0.1),
    );
    let accepted = selected
        .replace("[120.5]", "[999.0]")
        .replace("[180.5]", "[1999.0]");
    format!(
        "{{\"schema\":\"sdf3-affine-load-space-design.v1\",\"strength_constrained\":true,\"selected_design\":\"least-material-feasible\",\"selected_feasible\":true,\"selected\":{selected},\"last_accepted\":{accepted}}}"
    )
}

const RESPONSE_WITNESSES: &str = r#"[{"index":0,"name":"tip_lower","kind":"mean-displacement","unit":"m","sense":"at-least","coefficients":[-1,0]},{"index":1,"name":"support-upper","kind":"support-force","unit":"N","sense":"at-most","coefficients":[0.6,0.8]},{"index":2,"name":"moment-upper","kind":"support-moment","unit":"N*m","sense":"at-most","coefficients":[0,-1]}]"#;

fn affine_response_design(witnesses: &str) -> String {
    affine_strength_design()
        .replace("\"strength_constrained\":true", "\"strength_constrained\":true,\"response_constrained\":true")
        .replace("\"strength\":{", &format!(
            "\"responses\":{{\"schema\":\"sdf3-affine-response-fields.v1\",\"witnesses\":{witnesses}}},\"strength\":{{",
        ))
}

fn grid(driver: &str, design: &str) -> UnstructuredGrid {
    let Projection::Available(xml) = render(driver, design.as_bytes()).unwrap() else {
        panic!("a solved feasible design must have a field projection");
    };
    VtuChecker::parse_ascii(&xml).unwrap()
}

fn values<'a>(grid: &'a UnstructuredGrid, name: &str) -> &'a [f64] {
    let array = grid.arrays.iter().find(|a| a.name == name).unwrap();
    let DataValues::Float64(values) = &array.values else {
        panic!("f64 field");
    };
    values
}

#[test]
fn adaptive_export_preserves_physical_coordinates_voxel_order_and_case_vectors() {
    let grid = grid(SDF3_DRIVER, &adaptive());
    assert_eq!(grid.points[0], [1.0, 2.0, 3.0]);
    assert_eq!(grid.points[7], [2.0, 4.0, 6.0]);
    assert_eq!(grid.cells_connectivity, vec![0, 1, 2, 3, 4, 5, 6, 7]);
    assert_eq!(grid.cells_types, vec![CellType::Voxel as u8]);
    assert_eq!(values(&grid, "raw_density"), [0.7]);
    assert_eq!(values(&grid, "projected_density"), [0.6]);
    assert_eq!(values(&grid, "cut_volume_m3"), [2.5]);
    assert_eq!(values(&grid, "displacement_case_0"), displacement(0.0));
    assert_eq!(values(&grid, "displacement_case_1"), displacement(0.1));
    assert_eq!(grid.arrays.len(), 5);
    let vector = grid
        .arrays
        .iter()
        .find(|a| a.name == "displacement_case_1")
        .unwrap();
    assert_eq!(vector.components, 3);
    assert_eq!(vector.association, DataAssociation::PointData);
}

#[test]
fn stress_export_uses_selected_incumbent_and_sampled_stress() {
    let grid = grid(STRESS3_DRIVER, &incumbent(false));
    assert_eq!(values(&grid, "raw_density"), [0.7]);
    assert_eq!(values(&grid, "sampled_relaxed_max_stress"), [120.5]);
    assert_eq!(grid.arrays.len(), 6);
}

#[test]
fn stress_export_keeps_each_cells_maximum_over_the_load_cases() {
    // A second voxel beside the first in x; the first cell's larger stress
    // sits in case 0, the second cell's in case 1.
    let wide = |offset: f64| {
        let values: Vec<f64> = (0..36).map(|i| offset + f64::from(i) / 1000.0).collect();
        format!("{values:?}")
    };
    let two_cells = incumbent(false)
        .replace(NODES, "[[1,2,3],[2,2,3],[1,4,3],[2,4,3],[1,2,6],[2,2,6],[1,4,6],[2,4,6],[3,2,3],[3,4,3],[3,2,6],[3,4,6]]")
        .replace(CORNERS, "[[0,1,2,3,4,5,6,7],[1,8,3,9,5,10,7,11]]")
        .replace(
            "\"cut_volume_m3\":2.5}]",
            "\"cut_volume_m3\":2.5},{\"raw_density\":0.5,\"projected_density\":0.4,\"cut_volume_m3\":1.5}]",
        )
        .replace(&format!("{:?}", displacement(0.0)), &wide(0.0))
        .replace(&format!("{:?}", displacement(0.1)), &wide(0.1));
    let stress = |cases: &str| {
        two_cells.replace(
            "\"cell_relaxed_max_pa\":[[120.5],[80.0]]",
            &format!("\"cell_relaxed_max_pa\":{cases}"),
        )
    };
    let grid = grid(STRESS3_DRIVER, &stress("[[120.5,10.0],[80.0,50.0]]"));
    assert_eq!(values(&grid, "raw_density"), [0.7, 0.5]);
    assert_eq!(values(&grid, "sampled_relaxed_max_stress"), [120.5, 50.0]);
    for malformed in [
        "[120.5,50.0]",
        "[[120.5,10.0]]",
        "[[120.5,10.0],[80.0]]",
        "[[120.5,10.0],[80.0,-1.0]]",
    ] {
        assert!(
            render(STRESS3_DRIVER, stress(malformed).as_bytes()).is_err(),
            "accepted malformed sampled stress: {malformed}"
        );
    }
}

#[test]
fn affine_export_keeps_common_self_weight_separate_from_uncertain_basis_fields() {
    let affine = incumbent(true)
        .replace(
            "sdf3-load-space-design.v1",
            "sdf3-affine-load-space-design.v1",
        )
        .replace(
            "\"worst_displacement_m\":",
            &format!(
                "\"offset_displacement_m\":{:?},\"worst_displacement_m\":",
                displacement(-0.1)
            ),
        );
    let grid = grid(AFFINE_LOAD_SPACE3_DRIVER, &affine);
    assert_eq!(values(&grid, "offset_displacement"), displacement(-0.1));
    assert_eq!(values(&grid, "basis_displacement_0"), displacement(0.0));
    assert_eq!(values(&grid, "worst_displacement"), displacement(0.2));
    assert_eq!(grid.arrays.len(), 9);
    let missing = incumbent(true).replace(
        "sdf3-load-space-design.v1",
        "sdf3-affine-load-space-design.v1",
    );
    assert!(render(AFFINE_LOAD_SPACE3_DRIVER, missing.as_bytes()).is_err());
    assert!(render(LOAD_SPACE3_DRIVER, affine.as_bytes()).is_err());
}

#[test]
fn g3_affine_compliance_export_preserves_selected_fields_and_volume_goal() {
    let source = affine_strength_design().replace("least-material-feasible", "least-compliance-feasible");
    let goal = r#""goal":{"kind":"worst-compliance","volume_limit_fraction":0.85},"selected_design""#;
    let design = source.replace("\"selected_design\"", goal);
    let projected = grid(AFFINE_LOAD_SPACE3_DRIVER, &design);
    assert_eq!(values(&projected, "raw_density"), [0.7]);
    assert_eq!(values(&projected, "offset_displacement"), displacement(-0.1));
    assert_eq!(values(&projected, "sampled_physical_stress_envelope"), [120.5]);
    assert!(render(AFFINE_LOAD_SPACE3_DRIVER, source.as_bytes()).is_err());
    assert!(render(AFFINE_LOAD_SPACE3_DRIVER,
        design.replace("worst-compliance", "minimum-material").as_bytes()).is_err());
    assert!(render(AFFINE_LOAD_SPACE3_DRIVER,
        design.replace("0.85", "0.0").as_bytes()).is_err());
}

#[test]
fn affine_strength_export_preserves_selected_pointwise_envelopes_and_sample_counts() {
    let design = affine_strength_design();
    let grid = grid(AFFINE_LOAD_SPACE3_DRIVER, &design);
    assert_eq!(values(&grid, "sampled_physical_stress_envelope"), [120.5]);
    assert_eq!(values(&grid, "sampled_relaxed_stress_envelope"), [180.5]);
    assert_eq!(values(&grid, "worst_displacement"), displacement(0.2));
    assert_eq!(values(&grid, "offset_displacement"), displacement(-0.1));
    assert_eq!(grid.arrays.len(), 12);
    let counts = grid
        .arrays
        .iter()
        .find(|a| a.name == "stress_sample_count")
        .unwrap();
    assert_eq!(counts.association, DataAssociation::CellData);
    assert_eq!(counts.values, DataValues::Int64(vec![8]));
    for name in [
        "sampled_physical_stress_envelope",
        "sampled_relaxed_stress_envelope",
    ] {
        let stress = grid.arrays.iter().find(|a| a.name == name).unwrap();
        assert_eq!(stress.unit.as_deref(), Some("Pa"));
        assert_eq!(stress.association, DataAssociation::CellData);
    }
}

#[test]
fn affine_strength_export_refuses_missing_or_invalid_sampled_evidence() {
    let original = affine_strength_design();
    for malformed in [
        original.replace(
            "\"strength_constrained\":true",
            "\"strength_constrained\":false",
        ),
        original.replace("\"strength\":", "\"discarded_strength\":"),
        original.replace(
            "sdf3-affine-strength-fields.v1",
            "unsupported-strength-fields",
        ),
        original.replace(
            "\"cell_physical_envelope_pa\":[120.5]",
            "\"cell_physical_envelope_pa\":[]",
        ),
        original.replace("[120.5]", "[-1.0]"),
        original.replace("[180.5]", "[1e999]"),
        original.replace("\"cell_sample_count\":[8]", "\"cell_sample_count\":[-1]"),
        original.replace("\"cell_sample_count\":[8]", "\"cell_sample_count\":[8.5]"),
        original.replace("\"cell_sample_count\":[8]", "\"cell_sample_count\":[0]"),
    ] {
        assert!(
            render(AFFINE_LOAD_SPACE3_DRIVER, malformed.as_bytes()).is_err(),
            "accepted malformed affine strength fields: {malformed}"
        );
    }
}

#[test]
fn affine_response_export_keeps_every_signed_bound_adversary_and_common_offset() {
    let original = affine_response_design(RESPONSE_WITNESSES);
    let (selected, accepted) = original.split_once("\"last_accepted\":").unwrap();
    let design = format!(
        "{selected}\"last_accepted\":{}",
        accepted.replace("[-1,0]", "[0,1]")
    );
    let grid = grid(AFFINE_LOAD_SPACE3_DRIVER, &design);
    let offset = displacement(-0.1);
    let first = displacement(0.0);
    let second = displacement(0.1);
    for (name, coefficients) in [
        ("response_displacement_0_tip_lower", [-1.0, 0.0]),
        ("response_displacement_1_support-upper", [0.6, 0.8]),
        ("response_displacement_2_moment-upper", [0.0, -1.0]),
    ] {
        let expected: Vec<_> = (0..offset.len())
            .map(|i| (offset[i] + coefficients[0] * first[i]) + coefficients[1] * second[i])
            .collect();
        assert_eq!(values(&grid, name), expected);
        let array = grid.arrays.iter().find(|a| a.name == name).unwrap();
        assert_eq!(array.association, DataAssociation::PointData);
        assert_eq!(array.components, 3);
        assert_eq!(array.unit.as_deref(), Some("m"));
    }
    assert_eq!(grid.arrays.len(), 15);
    assert_eq!(values(&grid, "sampled_physical_stress_envelope"), [120.5]);
    assert_eq!(values(&grid, "worst_displacement"), displacement(0.2));
}

#[test]
fn affine_response_export_refuses_missing_reordered_or_invalid_adversaries() {
    let original = affine_response_design(RESPONSE_WITNESSES);
    for malformed in [
        original.replace(
            "\"response_constrained\":true",
            "\"response_constrained\":false",
        ),
        original.replace("\"responses\":", "\"missing_responses\":"),
        original.replace("sdf3-affine-response-fields.v1", "unknown-response-fields"),
        original.replace(RESPONSE_WITNESSES, "[]"),
        original.replace("\"index\":1", "\"index\":0"),
        original.replace("\"index\":1", "\"index\":1.0"),
        original.replace("support-upper", "tip_lower"),
        original.replace("tip_lower", "tip/lower"),
        original.replace("\"unit\":\"m\"", "\"unit\":\"N\""),
        original.replace("at-least", "absolute"),
        original.replace("[0.6,0.8]", "[0.6]"),
        original.replace("[0.6,0.8]", "[1,1]"),
        original.replace("[0.6,0.8]", "[1e999,0]"),
    ] {
        assert!(
            render(AFFINE_LOAD_SPACE3_DRIVER, malformed.as_bytes()).is_err(),
            "accepted invalid response witness fields: {malformed}"
        );
    }
    let oversized = (0..65).map(|index| format!(
        "{{\"index\":{index},\"name\":\"limit-{index}\",\"kind\":\"mean-displacement\",\"unit\":\"m\",\"sense\":\"at-most\",\"coefficients\":[1,0]}}"
    )).collect::<Vec<_>>().join(",");
    assert!(
        render(
            AFFINE_LOAD_SPACE3_DRIVER,
            affine_response_design(&format!("[{oversized}]")).as_bytes()
        )
        .is_err()
    );
    let huge = format!("{:?}", vec![f64::MAX; 24]);
    let overflow = original
        .replace("[-1,0]", "[1,0]")
        .replace(&format!("{:?}", displacement(-0.1)), &huge)
        .replace(&format!("{:?}", displacement(0.0)), &huge);
    assert!(render(AFFINE_LOAD_SPACE3_DRIVER, overflow.as_bytes()).is_err());
}

#[test]
fn affine_response_export_admits_the_complete_expansion_before_allocating_fields() {
    let design = J::parse(&affine_response_design(RESPONSE_WITNESSES)).unwrap();
    let fields = design.get("selected").unwrap();
    let mut grid = grid(AFFINE_LOAD_SPACE3_DRIVER, &affine_strength_design());
    let count = grid.arrays.len();
    // Three bound fields, each with 24 displacement scalars and two coefficients.
    let mut insufficient = 77;
    assert!(affine_responses(&design, fields, &mut grid, &mut insufficient).is_err());
    assert_eq!(grid.arrays.len(), count);
    let mut exact = 78;
    affine_responses(&design, fields, &mut grid, &mut exact).unwrap();
    assert_eq!(exact, 0);
    assert_eq!(grid.arrays.len(), count + 3);
}

#[test]
fn load_space_export_preserves_all_robust_densities_and_worst_combination() {
    let design = incumbent(true);
    let grid = grid(LOAD_SPACE3_DRIVER, &design);
    for (name, value) in [
        ("raw_density", 0.7),
        ("eroded_density", 0.4),
        ("nominal_density", 0.6),
        ("dilated_density", 0.8),
    ] {
        assert_eq!(values(&grid, name), [value]);
    }
    assert_eq!(values(&grid, "basis_displacement_0"), displacement(0.0));
    assert_eq!(values(&grid, "basis_displacement_1"), displacement(0.1));
    assert_eq!(values(&grid, "worst_displacement"), displacement(0.2));
    assert_eq!(grid.arrays.len(), 8);
    let Projection::Available(first) = render(LOAD_SPACE3_DRIVER, design.as_bytes()).unwrap()
    else {
        panic!()
    };
    let Projection::Available(second) = render(LOAD_SPACE3_DRIVER, design.as_bytes()).unwrap()
    else {
        panic!()
    };
    assert_eq!(first, second);
    assert!(first.contains("m^3"));
    assert!(first.contains("Unit=\"m\""));
}

#[test]
fn no_solved_or_feasible_incumbent_means_no_invented_field() {
    let empty = b"{\"schema\":\"sdf3-study-design.v1\",\"state\":\"no-solved-baseline\",\"cells\":[],\"displacements_m\":[]}";
    assert!(matches!(
        render(SDF3_DRIVER, empty).unwrap(),
        Projection::Unavailable("no-solved-baseline")
    ));
    for (driver, load_space) in [(STRESS3_DRIVER, false), (LOAD_SPACE3_DRIVER, true)] {
        let infeasible = incumbent(load_space)
            .replace("\"selected_feasible\":true", "\"selected_feasible\":false");
        assert!(matches!(
            render(driver, infeasible.as_bytes()).unwrap(),
            Projection::Unavailable("selected-design-not-feasible")
        ));
    }
}

#[test]
fn malformed_shape_nonfinite_values_and_wrong_driver_refuse() {
    let original = adaptive();
    for malformed in [
        original.replace(CORNERS, "[[0,1,3,2,4,5,7,6]]"),
        original.replace(CORNERS, "[[0,1,2,3,4,5,6,8]]"),
        original.replace(CORNERS, "[[0,1,2,3,4,5,6,7.0]]"),
        original.replace(CORNERS, "[]"),
        original.replace("[1,2,3]", "[1,2]"),
        original.replace("\"raw_density\":0.7", "\"raw_density\":1.1"),
        original.replace("\"cut_volume_m3\":2.5", "\"cut_volume_m3\":-1"),
        original.replace("\"cut_volume_m3\":2.5", "\"cut_volume_m3\":1e999"),
        original.replace("load-case then node-major xyz", "component-major xyz"),
        original.replace(&format!("{:?}", displacement(0.0)), "[0,0,0]"),
    ] {
        assert!(
            render(SDF3_DRIVER, malformed.as_bytes()).is_err(),
            "accepted malformed design: {malformed}"
        );
    }
    assert!(render(LOAD_SPACE3_DRIVER, original.as_bytes()).is_err());
    let missing_stress = incumbent(false).replace(
        "\"cell_relaxed_max_pa\":[[120.5],[80.0]]",
        "\"cell_relaxed_max_pa\":[]",
    );
    assert!(render(STRESS3_DRIVER, missing_stress.as_bytes()).is_err());
}

#[test]
fn field_case_envelope_refuses_unbounded_expansion() {
    let vector = format!("{:?}", displacement(0.0));
    let original = adaptive();
    let malformed = original.replace(
        &format!("[{:?},{:?}]", displacement(0.0), displacement(0.1)),
        &format!("[{}]", vec![vector; 65].join(",")),
    );
    assert!(render(SDF3_DRIVER, malformed.as_bytes()).is_err());
}
