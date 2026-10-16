//! Conformance tests for CAE ecosystem interchange tiers (bead frankensim-extreal-program-f85xj.11.5).

use fs_io::{
    CAE_CAPABILITY_MATRIX, CaeDirection, CaeQuarantineStatus, GmshElement, GmshElementType,
    GmshMesh, GmshNode, TabularColumn, TabularDataset, export_arrow_ipc_stream, export_tabular_csv,
    parse_abaqus_inp, parse_gmsh, parse_nastran_bdf, write_gmsh_msh2,
};

#[test]
fn test_gmsh_roundtrip_msh2() {
    let mut mesh = GmshMesh::new("2.2");
    mesh.physical_names.insert((2, 100), "inlet".to_string());
    mesh.physical_names
        .insert((3, 200), "solid_domain".to_string());

    mesh.nodes.push(GmshNode {
        id: 1,
        coords: [0.0, 0.0, 0.0],
    });
    mesh.nodes.push(GmshNode {
        id: 2,
        coords: [1.0, 0.0, 0.0],
    });
    mesh.nodes.push(GmshNode {
        id: 3,
        coords: [0.0, 1.0, 0.0],
    });
    mesh.nodes.push(GmshNode {
        id: 4,
        coords: [0.0, 0.0, 1.0],
    });

    mesh.elements.push(GmshElement {
        id: 1,
        element_type: GmshElementType::Triangle,
        tags: vec![100, 1],
        node_ids: vec![1, 2, 3],
    });
    mesh.elements.push(GmshElement {
        id: 2,
        element_type: GmshElementType::Tetrahedron,
        tags: vec![200, 2],
        node_ids: vec![1, 2, 3, 4],
    });

    let msh_ascii = write_gmsh_msh2(&mesh);
    assert!(msh_ascii.contains("$MeshFormat"));
    assert!(msh_ascii.contains("$Nodes"));
    assert!(msh_ascii.contains("$Elements"));

    let (quarantined, receipt) = parse_gmsh(&msh_ascii).expect("parse gmsh");
    assert_eq!(receipt.node_count, 4);
    assert_eq!(receipt.element_count, 2);
    assert_eq!(receipt.physical_group_count, 2);

    let parsed_mesh = quarantined.into_inner();
    assert_eq!(parsed_mesh.nodes.len(), 4);
    assert_eq!(parsed_mesh.elements.len(), 2);
    assert_eq!(
        parsed_mesh.elements[1].element_type,
        GmshElementType::Tetrahedron
    );
}

#[test]
fn test_abaqus_inp_thermal_extraction() {
    let inp = r#"
** Simple thermal test
*HEADING
Model: Heatsink
*NODE, NSET=ALL_NODES
1, 0.0, 0.0, 0.0
2, 0.1, 0.0, 0.0
3, 0.0, 0.1, 0.0
4, 0.0, 0.0, 0.1
*ELEMENT, TYPE=DC3D4, ELSET=SOLID_PART
1, 1, 2, 3, 4
*NSET, NSET=BASE_SURFACE
1, 2, 3
*MATERIAL, NAME=ALUMINUM
*CONDUCTIVITY
167.0
*SPECIFIC HEAT
896.0
*DENSITY
2700.0
*BOUNDARY
BASE_SURFACE, 11, 11, 353.15
*DYNAMIC, EXPLICIT
0.01, 1.0
"#;

    let (quarantined, receipt) = parse_abaqus_inp(inp).expect("parse inp");
    assert_eq!(receipt.node_count, 4);
    assert_eq!(receipt.element_count, 1);
    assert_eq!(receipt.material_count, 1);
    assert_eq!(receipt.unsupported_card_count, 2); // *HEADING, *DYNAMIC

    let model = quarantined.into_inner();
    assert_eq!(model.materials.len(), 1);
    assert_eq!(model.materials[0].name, "ALUMINUM");
    assert_eq!(model.materials[0].thermal_conductivity, Some(167.0));
    assert_eq!(model.materials[0].specific_heat, Some(896.0));
    assert_eq!(model.materials[0].density, Some(2700.0));
    assert_eq!(model.boundary_conditions.len(), 1);
}

#[test]
fn test_nastran_bdf_thermal_extraction() {
    let bdf = r#"
$ NASTRAN BDF Thermal Deck
GRID,1,,0.0,0.0,0.0
GRID,2,,0.05,0.0,0.0
GRID,3,,0.0,0.05,0.0
GRID,4,,0.0,0.0,0.05
CTETRA,101,1,1,2,3,4
MAT4,1,200.0,900.0,2700.0
TEMP,1,340.0
FORCE,1,1,0,100.0,0.0,0.0,1.0
"#;

    let (quarantined, receipt) = parse_nastran_bdf(bdf).expect("parse bdf");
    assert_eq!(receipt.node_count, 4);
    assert_eq!(receipt.element_count, 1);
    assert_eq!(receipt.material_count, 1);
    assert_eq!(receipt.unsupported_card_count, 1); // FORCE

    let model = quarantined.into_inner();
    assert_eq!(model.materials[0].thermal_conductivity, Some(200.0));
    assert_eq!(model.boundary_conditions.len(), 1);
}

/// G0 regression for bead frankensim-egdbd: receipt identity must bind the
/// admitted values, not merely aggregate node/element/card counts.
#[test]
fn cae_receipts_distinguish_equal_count_models() {
    let abaqus_a = "*NODE\n1,0,0,0\n*MATERIAL,NAME=M\n*CONDUCTIVITY\n10\n";
    let abaqus_b = "*NODE\n1,0,0,0\n*MATERIAL,NAME=M\n*CONDUCTIVITY\n20\n";
    let (abaqus_a_model, abaqus_a_receipt) = parse_abaqus_inp(abaqus_a).expect("Abaqus A");
    let (abaqus_b_model, abaqus_b_receipt) = parse_abaqus_inp(abaqus_b).expect("Abaqus B");
    assert_ne!(abaqus_a_receipt.content_hash, abaqus_b_receipt.content_hash);
    assert_ne!(
        abaqus_a_model.source_receipt.source_hash,
        abaqus_b_model.source_receipt.source_hash
    );
    assert_eq!(
        abaqus_a_model.source_receipt.source_hash,
        fs_obs::fnv1a64(abaqus_a.as_bytes())
    );
    let (abaqus_a_again, abaqus_a_receipt_again) =
        parse_abaqus_inp(abaqus_a).expect("Abaqus A replay");
    assert_eq!(
        abaqus_a_receipt.content_hash,
        abaqus_a_receipt_again.content_hash
    );
    assert_eq!(
        abaqus_a_model.source_receipt.source_hash,
        abaqus_a_again.source_receipt.source_hash
    );

    let bdf_a = "GRID,1,,0,0,0\nGRID,2,,1,0,0\nCTETRA,1,1,1,1,1,1\n";
    let bdf_b = "GRID,1,,0,0,0\nGRID,2,,1,0,0\nCTETRA,1,1,1,1,1,2\n";
    let (bdf_a_model, bdf_a_receipt) = parse_nastran_bdf(bdf_a).expect("BDF A");
    let (bdf_b_model, bdf_b_receipt) = parse_nastran_bdf(bdf_b).expect("BDF B");
    assert_ne!(bdf_a_receipt.content_hash, bdf_b_receipt.content_hash);
    assert_ne!(
        bdf_a_model.source_receipt.source_hash,
        bdf_b_model.source_receipt.source_hash
    );
    assert_eq!(
        bdf_a_model.source_receipt.source_hash,
        fs_obs::fnv1a64(bdf_a.as_bytes())
    );
    let (bdf_a_again, bdf_a_receipt_again) = parse_nastran_bdf(bdf_a).expect("BDF A replay");
    assert_eq!(bdf_a_receipt.content_hash, bdf_a_receipt_again.content_hash);
    assert_eq!(
        bdf_a_model.source_receipt.source_hash,
        bdf_a_again.source_receipt.source_hash
    );

    let mut gmsh_a = GmshMesh::new("2.2");
    gmsh_a.nodes.push(GmshNode {
        id: 1,
        coords: [0.0, 0.0, 0.0],
    });
    gmsh_a.elements.push(GmshElement {
        id: 1,
        element_type: GmshElementType::Point,
        tags: Vec::new(),
        node_ids: vec![1],
    });
    let mut gmsh_b = gmsh_a.clone();
    gmsh_b.nodes[0].coords[0] = 1.0;
    let gmsh_a_text = write_gmsh_msh2(&gmsh_a);
    let gmsh_b_text = write_gmsh_msh2(&gmsh_b);
    let (gmsh_a_model, gmsh_a_receipt) = parse_gmsh(&gmsh_a_text).expect("Gmsh A");
    let (gmsh_b_model, gmsh_b_receipt) = parse_gmsh(&gmsh_b_text).expect("Gmsh B");
    assert_ne!(gmsh_a_receipt.content_hash, gmsh_b_receipt.content_hash);
    assert_ne!(
        gmsh_a_model.source_receipt.source_hash,
        gmsh_b_model.source_receipt.source_hash
    );
    assert_eq!(
        gmsh_a_model.source_receipt.source_hash,
        fs_obs::fnv1a64(gmsh_a_text.as_bytes())
    );

    let (gmsh_a_again, gmsh_a_receipt_again) = parse_gmsh(&gmsh_a_text).expect("Gmsh A replay");
    assert_eq!(
        gmsh_a_receipt.content_hash,
        gmsh_a_receipt_again.content_hash
    );
    assert_eq!(
        gmsh_a_model.source_receipt.source_hash,
        gmsh_a_again.source_receipt.source_hash
    );
}

#[test]
fn test_tabular_csv_and_arrow_ipc_export() {
    let temp_col = TabularColumn::new_f64("temperature", "K", vec![300.0, 320.5, 345.2, 380.1]);
    let flux_col =
        TabularColumn::new_f64("heat_flux_z", "W/m^2", vec![1000.0, 1500.0, 2200.0, 3100.0]);
    let node_col = TabularColumn::new_i64("node_id", "id", vec![1, 2, 3, 4]);

    let dataset = TabularDataset::new("thermal_results", vec![node_col, temp_col, flux_col])
        .expect("valid dataset");

    // CSV export
    let (csv_text, csv_receipt) = export_tabular_csv(&dataset);
    assert!(csv_text.contains("temperature,heat_flux_z"));
    assert!(csv_text.contains("K,W/m^2"));
    assert_eq!(csv_receipt.row_count, 4);
    assert_eq!(csv_receipt.column_count, 3);

    // Arrow IPC export
    let (arrow_bytes, arrow_receipt) = export_arrow_ipc_stream(&dataset);
    assert!(arrow_bytes.starts_with(b"ARROW1\0\0"));
    assert_eq!(arrow_receipt.row_count, 4);
    assert_eq!(arrow_receipt.column_count, 3);
}

#[test]
fn test_cae_capability_matrix_completeness() {
    assert!(CAE_CAPABILITY_MATRIX.len() >= 8);

    let gmsh_entry = CAE_CAPABILITY_MATRIX
        .iter()
        .find(|e| e.format_id == "Gmsh-MSH")
        .expect("gmsh in matrix");
    assert_eq!(
        gmsh_entry.quarantine_status,
        CaeQuarantineStatus::NativeCertified
    );
    assert_eq!(gmsh_entry.direction, CaeDirection::Bidirectional);

    let exo_entry = CAE_CAPABILITY_MATRIX
        .iter()
        .find(|e| e.format_id == "Exodus-II")
        .expect("exodus in matrix");
    assert_eq!(
        exo_entry.quarantine_status,
        CaeQuarantineStatus::QuarantinedAdapter
    );
}
