//! Report projects the selected retained design through the shared VTU writer.
use super::*;
use fs_viz::vtu::{CellType, DataValues, UnstructuredGrid, VtuChecker};

fn values<'a>(grid: &'a UnstructuredGrid, name: &str) -> &'a [f64] {
    let array = grid.arrays.iter().find(|a| a.name == name).unwrap();
    let DataValues::Float64(values) = &array.values else {
        panic!("f64 study field");
    };
    values
}

fn numbers(value: &J) -> Vec<f64> {
    value
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_f64().unwrap())
        .collect()
}

#[test]
fn report_exports_selected_robust_fields_from_ledger_after_source_is_relocated() {
    let source = include_str!("../../../../examples/marquee/bracket-3d-load-space.fsim")
        .replace(":max-updates 80", ":max-updates 2")
        .replace(":compliance-limit-j 10.0", ":compliance-limit-j 10000.0");
    let dir = scratch("load-space-vtu");
    let (ledger, result) = run(&dir, "robust", &source, 6);
    let design = data(&ledger, &result, "design");
    assert_eq!(design.get("selected_feasible"), Some(&J::Bool(true)));
    let selected = design.get("selected").unwrap();
    fs::rename(dir.join("robust.fsim"), dir.join("relocated.fsim")).unwrap();
    let pointer = result.str_field("run_id").unwrap();
    let exported = document(
        &command("report")
            .arg(pointer)
            .arg(&ledger)
            .output()
            .unwrap(),
        0,
    );
    assert_eq!(exported.str_field("field_vtu_status"), Some("available"));
    assert_eq!(
        exported.get("field_vtu_source_design"),
        result.path(&["receipt", "design"])
    );
    assert!(
        exported
            .str_field("field_vtu_scope")
            .unwrap()
            .contains("background")
    );
    let path = exported.str_field("field_vtu").unwrap();
    let xml = fs::read_to_string(path).unwrap();
    let grid = VtuChecker::parse_ascii(&xml).unwrap();
    let nodes = selected
        .get("physical_nodes_m")
        .and_then(J::as_array)
        .unwrap();
    assert_eq!(grid.points.len(), nodes.len());
    for (actual, expected) in grid.points.iter().zip(nodes) {
        assert_eq!(actual.as_slice(), numbers(expected));
    }
    let cells = selected.get("cells").and_then(J::as_array).unwrap();
    assert_eq!(grid.cells_types, vec![CellType::Voxel as u8; cells.len()]);
    let connectivity: Vec<_> = selected
        .get("cell_nodes")
        .and_then(J::as_array)
        .unwrap()
        .iter()
        .flat_map(|c| c.as_array().unwrap())
        .map(|v| v.number_raw().unwrap().parse::<i64>().unwrap())
        .collect();
    assert_eq!(grid.cells_connectivity, connectivity);
    for name in [
        "raw_density",
        "eroded_density",
        "nominal_density",
        "dilated_density",
        "cut_volume_m3",
    ] {
        let expected: Vec<_> = cells.iter().map(|c| number(c, name)).collect();
        assert_eq!(values(&grid, name), expected);
    }
    let basis = selected
        .get("basis_displacements_m")
        .and_then(J::as_array)
        .unwrap();
    for (i, expected) in basis.iter().enumerate() {
        assert_eq!(
            values(&grid, &format!("basis_displacement_{i}")),
            numbers(expected)
        );
    }
    assert_eq!(
        values(&grid, "worst_displacement"),
        numbers(selected.get("worst_displacement_m").unwrap())
    );
    let repeated = document(
        &command("report")
            .arg(pointer)
            .arg(&ledger)
            .output()
            .unwrap(),
        0,
    );
    assert_eq!(
        fs::read_to_string(repeated.str_field("field_vtu").unwrap()).unwrap(),
        xml
    );
}
