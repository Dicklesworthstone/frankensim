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
            "\"displacement_layout\":\"load-case then node-major xyz\",\"displacements_m\":[{:?},{:?}],\"cell_relaxed_max_pa\":[120.5]",
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
        "\"cell_relaxed_max_pa\":[120.5]",
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
