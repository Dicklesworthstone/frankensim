//! Integration tests for VTU and XDMF visualization export.
//!
//! Bead: `frankensim-extreal-program-f85xj.6.8`

use fs_viz::{
    CellType, DataArray, ExportFormat, FIELD_REGISTRY, UnstructuredGrid, VtuChecker, VtuError,
    VtuWriter, XdmfWriter, find_field_by_name, validate_output_request,
};

#[test]
fn test_vtu_basic_tetra_ascii_roundtrip() {
    let mut grid = UnstructuredGrid::new();
    let p0 = grid.add_point(0.0, 0.0, 0.0);
    let p1 = grid.add_point(1.0, 0.0, 0.0);
    let p2 = grid.add_point(0.0, 1.0, 0.0);
    let p3 = grid.add_point(0.0, 0.0, 1.0);

    grid.add_tetra(p0, p1, p2, p3);

    let temp_vals = vec![300.0, 310.0, 320.0, 350.0];
    grid.add_array(DataArray::new_point_scalar("Temperature", temp_vals).with_unit("K"));

    let region_vals = vec![1i32];
    grid.add_array(DataArray::new_cell_int32("RegionId", region_vals));

    let xml = VtuWriter::write_ascii(&grid).expect("write VTU");
    assert!(xml.contains("<VTKFile type=\"UnstructuredGrid\""));
    assert!(xml.contains("NumberOfPoints=\"4\""));
    assert!(xml.contains("NumberOfCells=\"1\""));
    assert!(xml.contains("Name=\"Temperature\""));
    assert!(xml.contains("Name=\"RegionId\""));

    let report = VtuChecker::check(&xml).expect("check VTU");
    assert_eq!(report.num_points, 4);
    assert_eq!(report.num_cells, 1);
    assert_eq!(report.array_count, 2);

    let temp_extrema = report
        .array_extrema
        .iter()
        .find(|(name, _)| name == "Temperature")
        .expect("temperature in report");
    let arr_eq_bits = |a: &[f64], b: &[f64]| {
        a.len() == b.len()
            && a.iter()
                .zip(b.iter())
                .all(|(x, y)| x.to_bits() == y.to_bits())
    };
    assert!(arr_eq_bits(&temp_extrema.1, &[300.0, 350.0]));

    // Independent readback
    let parsed_grid = VtuChecker::parse_ascii(&xml).expect("parse grid");
    assert_eq!(parsed_grid.points, grid.points);
    assert_eq!(parsed_grid.cells_connectivity, grid.cells_connectivity);
    assert_eq!(parsed_grid.cells_types, grid.cells_types);
}

#[test]
fn test_vtu_determinism_same_input() {
    let mut grid1 = UnstructuredGrid::new();
    let p0 = grid1.add_point(0.0, 0.0, 0.0);
    let p1 = grid1.add_point(1.0, 0.0, 0.0);
    let p2 = grid1.add_point(0.0, 1.0, 0.0);
    let p3 = grid1.add_point(0.0, 0.0, 1.0);
    grid1.add_tetra(p0, p1, p2, p3);
    grid1.add_array(DataArray::new_point_scalar(
        "Temperature",
        vec![300.0, 310.0, 320.0, 350.0],
    ));

    let mut grid2 = UnstructuredGrid::new();
    let q0 = grid2.add_point(0.0, 0.0, 0.0);
    let q1 = grid2.add_point(1.0, 0.0, 0.0);
    let q2 = grid2.add_point(0.0, 1.0, 0.0);
    let q3 = grid2.add_point(0.0, 0.0, 1.0);
    grid2.add_tetra(q0, q1, q2, q3);
    grid2.add_array(DataArray::new_point_scalar(
        "Temperature",
        vec![300.0, 310.0, 320.0, 350.0],
    ));

    let xml1 = VtuWriter::write_ascii(&grid1).unwrap();
    let xml2 = VtuWriter::write_ascii(&grid2).unwrap();
    assert_eq!(xml1, xml2, "exact byte-level determinism");

    let hash1 = VtuWriter::content_hash(&grid1).unwrap();
    let hash2 = VtuWriter::content_hash(&grid2).unwrap();
    assert_eq!(hash1, hash2, "content hash identical");
}

#[test]
fn test_vtu_validation_failures() {
    // 1. Empty grid
    let empty_grid = UnstructuredGrid::new();
    assert_eq!(empty_grid.validate(), Err(VtuError::EmptyGrid));

    // 2. Non-finite point
    let mut bad_point_grid = UnstructuredGrid::new();
    bad_point_grid.add_point(f64::NAN, 0.0, 0.0);
    assert!(matches!(
        bad_point_grid.validate(),
        Err(VtuError::NonFinitePointCoordinate { .. })
    ));

    // 3. Out of bounds connectivity
    let mut bad_conn_grid = UnstructuredGrid::new();
    bad_conn_grid.add_point(0.0, 0.0, 0.0);
    bad_conn_grid.add_cell(CellType::Vertex, &[5]); // index 5 > 0
    assert!(matches!(
        bad_conn_grid.validate(),
        Err(VtuError::InvalidCellIndex { .. })
    ));

    // 4. Array length mismatch
    let mut bad_arr_grid = UnstructuredGrid::new();
    let p0 = bad_arr_grid.add_point(0.0, 0.0, 0.0);
    let p1 = bad_arr_grid.add_point(1.0, 0.0, 0.0);
    bad_arr_grid.add_cell(CellType::Line, &[p0, p1]);
    bad_arr_grid.add_array(DataArray::new_point_scalar("Temperature", vec![300.0])); // 1 item, expected 2
    assert!(matches!(
        bad_arr_grid.validate(),
        Err(VtuError::ArrayLengthMismatch { .. })
    ));

    // 5. Non-finite field value
    let mut nan_field_grid = UnstructuredGrid::new();
    let p0 = nan_field_grid.add_point(0.0, 0.0, 0.0);
    let p1 = nan_field_grid.add_point(1.0, 0.0, 0.0);
    nan_field_grid.add_cell(CellType::Line, &[p0, p1]);
    nan_field_grid.add_array(DataArray::new_point_scalar(
        "Temperature",
        vec![300.0, f64::INFINITY],
    ));
    assert!(matches!(
        nan_field_grid.validate(),
        Err(VtuError::NonFiniteFieldValue { .. })
    ));
}

#[test]
fn test_xdmf_export_and_binary_companion() {
    let mut grid = UnstructuredGrid::new();
    let p0 = grid.add_point(0.0, 0.0, 0.0);
    let p1 = grid.add_point(1.0, 0.0, 0.0);
    let p2 = grid.add_point(0.0, 1.0, 0.0);
    let p3 = grid.add_point(0.0, 0.0, 1.0);
    grid.add_tetra(p0, p1, p2, p3);
    grid.add_array(DataArray::new_point_scalar(
        "Temperature",
        vec![300.0, 310.0, 320.0, 350.0],
    ));

    let (xmf, bin) = XdmfWriter::write_xdmf_with_binary(&grid, "mesh.bin").expect("write XDMF");
    assert!(xmf.contains("<Xdmf Version=\"3.0\">"));
    assert!(xmf.contains("<Topology TopologyType=\"Tetrahedron\""));
    assert!(xmf.contains("mesh.bin"));
    assert!(!bin.is_empty());
    assert_eq!(bin.len(), (4 * 8) + (4 * 3 * 8) + (4 * 8)); // 4 conn + 4*3 coords + 4 temp
}

#[test]
fn test_field_registry() {
    assert!(!FIELD_REGISTRY.is_empty());

    let temp_desc = find_field_by_name("Temperature").expect("Temperature found");
    assert_eq!(temp_desc.semantic_id, "field.thermal.temperature");
    assert_eq!(temp_desc.default_unit, "K");

    let flux_desc = find_field_by_name("HeatFlux").expect("HeatFlux found");
    assert_eq!(flux_desc.components, 3);

    assert!(validate_output_request("Temperature", ExportFormat::Vtu).is_ok());
    assert!(validate_output_request("UnknownField123", ExportFormat::Vtu).is_err());
}
