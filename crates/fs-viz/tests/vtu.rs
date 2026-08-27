//! Integration tests for VTU and XDMF visualization export.
//!
//! Bead: `frankensim-extreal-program-f85xj.6.8`

use fs_viz::{
    CellType, DataArray, DataAssociation, DataValues, ExportFormat, FIELD_REGISTRY,
    UnstructuredGrid, VtuChecker, VtuError, VtuWriter, XdmfWriter, find_field_by_name,
    validate_output_request,
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
fn test_vtu_validation_rejects_connectivity_length_mismatches_without_panicking() {
    let mut out_of_bounds = UnstructuredGrid::new();
    out_of_bounds.add_point(0.0, 0.0, 0.0);
    out_of_bounds.cells_connectivity = vec![0];
    out_of_bounds.cells_offsets = vec![2];
    out_of_bounds.cells_types = vec![CellType::Line.type_id()];
    assert_eq!(
        out_of_bounds.validate(),
        Err(VtuError::ConnectivityOffsetOutOfBounds {
            cell: 0,
            offset: 2,
            connectivity_len: 1,
        })
    );

    let mut trailing = UnstructuredGrid::new();
    trailing.add_point(0.0, 0.0, 0.0);
    trailing.cells_connectivity = vec![0, 0];
    trailing.cells_offsets = vec![1];
    trailing.cells_types = vec![CellType::Vertex.type_id()];
    assert_eq!(
        trailing.validate(),
        Err(VtuError::ConnectivityLengthMismatch {
            referenced: 1,
            available: 2,
        })
    );
}

#[test]
fn test_vtu_validation_rejects_array_item_count_overflow() {
    let mut grid = UnstructuredGrid::new();
    grid.add_point(0.0, 0.0, 0.0);
    grid.add_point(1.0, 0.0, 0.0);
    grid.add_array(DataArray {
        name: "overflow".to_string(),
        association: DataAssociation::PointData,
        components: usize::MAX,
        unit: None,
        values: DataValues::Float64(Vec::new()),
    });
    assert_eq!(
        grid.validate(),
        Err(VtuError::ArrayItemCountOverflow {
            array: "overflow".to_string(),
            entities: 2,
            components: usize::MAX,
        })
    );
}

#[test]
fn test_vtu_checker_preserves_every_writer_data_type() {
    let mut grid = UnstructuredGrid::new();
    let p0 = grid.add_point(0.0, 0.0, 0.0);
    let p1 = grid.add_point(1.0, 0.0, 0.0);
    let p2 = grid.add_point(0.0, 1.0, 0.0);
    let p3 = grid.add_point(0.0, 0.0, 1.0);
    grid.add_tetra(p0, p1, p2, p3);
    grid.add_array(DataArray {
        name: "f32".to_string(),
        association: DataAssociation::PointData,
        components: 1,
        unit: None,
        values: DataValues::Float32(vec![1.25, 2.5, 3.75, 5.0]),
    });
    grid.add_array(DataArray {
        name: "i64".to_string(),
        association: DataAssociation::PointData,
        components: 1,
        unit: None,
        values: DataValues::Int64(vec![i64::MIN, -1, 0, i64::MAX]),
    });
    grid.add_array(DataArray {
        name: "u8".to_string(),
        association: DataAssociation::PointData,
        components: 1,
        unit: None,
        values: DataValues::UInt8(vec![0, 1, 254, 255]),
    });

    let xml = VtuWriter::write_ascii(&grid).expect("write VTU");
    let parsed = VtuChecker::parse_ascii(&xml).expect("parse every emitted VTU data type");
    assert_eq!(parsed.arrays, grid.arrays);
}

#[test]
fn test_vtu_xml_attributes_are_escaped_and_round_trip() {
    let mut grid = UnstructuredGrid::new();
    let p0 = grid.add_point(0.0, 0.0, 0.0);
    grid.add_cell(CellType::Vertex, &[p0]);
    grid.add_array(
        DataArray::new_point_scalar("temperature <&\"'>", vec![300.0]).with_unit("K <&\"'>"),
    );

    let xml = VtuWriter::write_ascii(&grid).expect("write escaped VTU");
    assert!(xml.contains("Name=\"temperature &lt;&amp;&quot;&apos;&gt;\""));
    assert!(xml.contains("Unit=\"K &lt;&amp;&quot;&apos;&gt;\""));
    let parsed = VtuChecker::parse_ascii(&xml).expect("parse escaped VTU");
    assert_eq!(parsed.arrays, grid.arrays);
}

#[test]
fn test_vtu_checker_rejects_tampered_piece_counts() {
    let mut grid = UnstructuredGrid::new();
    let p0 = grid.add_point(0.0, 0.0, 0.0);
    grid.add_cell(CellType::Vertex, &[p0]);
    let xml = VtuWriter::write_ascii(&grid).expect("write VTU");

    let bad_points = xml.replacen("NumberOfPoints=\"1\"", "NumberOfPoints=\"2\"", 1);
    assert!(matches!(
        VtuChecker::check(&bad_points),
        Err(VtuError::ParseError { detail }) if detail.contains("declares 2 points")
    ));

    let bad_cells = xml.replacen("NumberOfCells=\"1\"", "NumberOfCells=\"2\"", 1);
    assert!(matches!(
        VtuChecker::check(&bad_cells),
        Err(VtuError::ParseError { detail }) if detail.contains("declares 2 cells")
    ));
}

#[test]
fn test_vtu_checker_requires_exact_cell_array_names() {
    let mut grid = UnstructuredGrid::new();
    let p0 = grid.add_point(0.0, 0.0, 0.0);
    grid.add_cell(CellType::Vertex, &[p0]);
    let xml = VtuWriter::write_ascii(&grid).expect("write VTU");
    let tampered = xml.replacen("Name=\"connectivity\"", "AliasName=\"connectivity\"", 1);

    assert!(matches!(
        VtuChecker::check(&tampered),
        Err(VtuError::ConnectivityOffsetOutOfBounds {
            cell: 0,
            offset: 1,
            connectivity_len: 0,
        })
    ));
}

#[test]
fn test_vtu_checker_rejects_invalid_component_metadata() {
    let mut grid = UnstructuredGrid::new();
    let p0 = grid.add_point(0.0, 0.0, 0.0);
    grid.add_cell(CellType::Vertex, &[p0]);
    grid.add_array(DataArray::new_point_scalar("scalar", vec![1.0]));
    let xml = VtuWriter::write_ascii(&grid).expect("write VTU");
    let tampered = xml.replacen(
        "NumberOfComponents=\"1\"",
        "NumberOfComponents=\"not-a-number\"",
        1,
    );

    assert!(matches!(
        VtuChecker::check(&tampered),
        Err(VtuError::ParseError { detail }) if detail.contains("invalid NumberOfComponents")
    ));
}

#[test]
fn test_vtu_checker_rejects_false_core_array_metadata() {
    let mut grid = UnstructuredGrid::new();
    let p0 = grid.add_point(0.0, 0.0, 0.0);
    grid.add_cell(CellType::Vertex, &[p0]);
    let xml = VtuWriter::write_ascii(&grid).expect("write VTU");

    for (tampered, expected_detail) in [
        (
            xml.replacen(
                "type=\"Float64\" Name=\"Points\"",
                "type=\"Float32\" Name=\"Points\"",
                1,
            ),
            "Points DataArray `type` must be `Float64`",
        ),
        (
            xml.replacen(
                "Name=\"Points\" NumberOfComponents=\"3\"",
                "Name=\"Points\" NumberOfComponents=\"2\"",
                1,
            ),
            "Points DataArray `NumberOfComponents` must be `3`",
        ),
        (
            xml.replacen(
                "type=\"Int64\" Name=\"connectivity\" format=\"ascii\"",
                "type=\"Int32\" Name=\"connectivity\" format=\"ascii\"",
                1,
            ),
            "Cells connectivity DataArray `type` must be `Int64`",
        ),
        (
            xml.replacen(
                "type=\"UInt8\" Name=\"types\" format=\"ascii\"",
                "type=\"UInt8\" Name=\"types\" format=\"binary\"",
                1,
            ),
            "Cells types DataArray `format` must be `ascii`",
        ),
    ] {
        assert_ne!(tampered, xml, "test fixture must mutate writer output");
        assert!(matches!(
            VtuChecker::check(&tampered),
            Err(VtuError::ParseError { detail }) if detail.contains(expected_detail)
        ));
    }
}

#[test]
fn test_vtu_checker_rejects_non_ascii_field_encoding() {
    let mut grid = UnstructuredGrid::new();
    let p0 = grid.add_point(0.0, 0.0, 0.0);
    grid.add_cell(CellType::Vertex, &[p0]);
    grid.add_array(DataArray::new_point_scalar("scalar", vec![1.0]));
    let xml = VtuWriter::write_ascii(&grid).expect("write VTU");
    let tampered = xml.replacen(
        "Name=\"scalar\" NumberOfComponents=\"1\" format=\"ascii\"",
        "Name=\"scalar\" NumberOfComponents=\"1\" format=\"binary\"",
        1,
    );
    assert_ne!(tampered, xml, "test fixture must mutate writer output");

    assert!(matches!(
        VtuChecker::check(&tampered),
        Err(VtuError::ParseError { detail })
            if detail.contains("array `scalar` DataArray `format` must be `ascii`")
    ));
}

#[test]
fn test_vtu_duplicate_names_are_scoped_to_their_association() {
    let mut grid = UnstructuredGrid::new();
    let p0 = grid.add_point(0.0, 0.0, 0.0);
    grid.add_cell(CellType::Vertex, &[p0]);
    grid.add_array(DataArray::new_point_scalar("id", vec![1.0]));
    grid.add_array(DataArray::new_cell_scalar("id", vec![2.0]));

    let xml = VtuWriter::write_ascii(&grid).expect("same name is valid across associations");
    let parsed = VtuChecker::parse_ascii(&xml).expect("round-trip association-scoped names");
    assert_eq!(parsed.arrays, grid.arrays);
}

#[test]
fn test_vtu_rejects_xml_forbidden_characters() {
    let mut grid = UnstructuredGrid::new();
    let p0 = grid.add_point(0.0, 0.0, 0.0);
    grid.add_cell(CellType::Vertex, &[p0]);
    grid.add_array(DataArray::new_point_scalar("bad\0name", vec![1.0]));

    assert_eq!(
        VtuWriter::write_ascii(&grid),
        Err(VtuError::InvalidXmlCharacter {
            context: "data array name",
            index: 3,
            codepoint: 0,
        })
    );
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
    assert!(xmf.contains("NumberType=\"Float\" Precision=\"8\""));
    assert!(xmf.contains("Endian=\"Little\""));
    assert!(!bin.is_empty());
    assert_eq!(bin.len(), (4 * 8) + (4 * 3 * 8) + (4 * 8)); // 4 conn + 4*3 coords + 4 temp
}

#[test]
fn test_xdmf_field_metadata_matches_binary_value_types() {
    let mut grid = UnstructuredGrid::new();
    let p0 = grid.add_point(0.0, 0.0, 0.0);
    let p1 = grid.add_point(1.0, 0.0, 0.0);
    let p2 = grid.add_point(0.0, 1.0, 0.0);
    let p3 = grid.add_point(0.0, 0.0, 1.0);
    grid.add_tetra(p0, p1, p2, p3);
    for (name, values) in [
        ("f32", DataValues::Float32(vec![1.0, 2.0, 3.0, 4.0])),
        ("i32", DataValues::Int32(vec![1, 2, 3, 4])),
        ("i64", DataValues::Int64(vec![1, 2, 3, 4])),
        ("u8", DataValues::UInt8(vec![1, 2, 3, 4])),
    ] {
        grid.add_array(DataArray {
            name: name.to_string(),
            association: DataAssociation::PointData,
            components: 1,
            unit: None,
            values,
        });
    }

    let (xmf, _) = XdmfWriter::write_xdmf_with_binary(&grid, "typed.bin").expect("write XDMF");
    for (name, number_type, precision) in [
        ("f32", "Float", 4),
        ("i32", "Int", 4),
        ("i64", "Int", 8),
        ("u8", "UInt", 1),
    ] {
        let attribute_start = xmf
            .find(&format!("<Attribute Name=\"{name}\""))
            .expect("attribute exists");
        let attribute = &xmf[attribute_start..];
        let attribute_end = attribute.find("</Attribute>").expect("attribute closes");
        let attribute = &attribute[..attribute_end];
        assert!(
            attribute.contains(&format!(
                "NumberType=\"{number_type}\" Precision=\"{precision}\""
            )),
            "wrong XDMF type metadata for {name}: {attribute}"
        );
    }
}

#[test]
fn test_xdmf_escapes_attribute_names_and_binary_file_text() {
    let mut grid = UnstructuredGrid::new();
    let p0 = grid.add_point(0.0, 0.0, 0.0);
    grid.add_cell(CellType::Vertex, &[p0]);
    grid.add_array(DataArray::new_point_scalar("field <&\"'>", vec![1.0]));

    let (xmf, _) = XdmfWriter::write_xdmf_with_binary(&grid, "mesh<&>\tline\nreturn\r.bin")
        .expect("write escaped XDMF");
    assert!(xmf.contains("Name=\"field &lt;&amp;&quot;&apos;&gt;\""));
    assert!(xmf.contains("mesh&lt;&amp;&gt;&#x9;line&#xA;return&#xD;.bin"));
    assert!(!xmf.contains("mesh<&>\tline\nreturn\r.bin"));
}

#[test]
fn test_xdmf_mixed_topology_includes_cell_codes() {
    let mut grid = UnstructuredGrid::new();
    let p0 = grid.add_point(0.0, 0.0, 0.0);
    let p1 = grid.add_point(1.0, 0.0, 0.0);
    let p2 = grid.add_point(0.0, 1.0, 0.0);
    let p3 = grid.add_point(0.0, 0.0, 1.0);
    grid.add_tetra(p0, p1, p2, p3);
    grid.add_triangle(p0, p1, p2);

    let (xmf, bin) =
        XdmfWriter::write_xdmf_with_binary(&grid, "mixed.bin").expect("write mixed XDMF topology");
    assert!(xmf.contains("TopologyType=\"Mixed\" NumberOfElements=\"2\""));
    assert!(xmf.contains("Dimensions=\"9\""));

    let topology: Vec<u64> = bin[..9 * 8]
        .chunks_exact(8)
        .map(|bytes| u64::from_le_bytes(bytes.try_into().expect("eight bytes")))
        .collect();
    assert_eq!(topology, vec![6, 0, 1, 2, 3, 4, 0, 1, 2]);
}

#[test]
fn test_xdmf_refuses_cell_types_without_lossless_mapping() {
    let mut grid = UnstructuredGrid::new();
    let p0 = grid.add_point(0.0, 0.0, 0.0);
    let p1 = grid.add_point(1.0, 0.0, 0.0);
    let p2 = grid.add_point(0.0, 1.0, 0.0);
    grid.add_cell(CellType::TriangleStrip, &[p0, p1, p2]);

    assert_eq!(
        XdmfWriter::write_xdmf_with_binary(&grid, "strip.bin"),
        Err(VtuError::UnsupportedXdmfCellType {
            cell: 0,
            type_id: CellType::TriangleStrip.type_id(),
        })
    );
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
