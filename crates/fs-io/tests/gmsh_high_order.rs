//! IRP M1.3-W4: second-order Gmsh solids (ids 11, 17, 12) and their node
//! orders. The Gmsh reference coordinates below are transcribed from the Gmsh
//! manual's node-ordering figures (https://gmsh.info/doc/texinfo/#Node-ordering),
//! independently of the permutation tables they pin. The imported-mesh test
//! reads a committed fixture made offline by an out-of-process gmsh adapter.

use fs_io::gmsh::{
    GMSH_HEX20_TO_CANONICAL, GMSH_HEX27_TO_CANONICAL, GMSH_TET10_TO_CANONICAL, GmshElement,
    GmshElementType, GmshLimits, GmshMesh, GmshNode, GmshPromotionRefusal, GmshSolidDeclaration,
    parse_gmsh, parse_gmsh_with_limits, promote_gmsh_solid, write_gmsh_msh2,
};
use fs_io::{IoError, LengthUnit, ReferenceState};
use fs_solid::linear3::hex::{
    GMSH20_TO_CANONICAL, GMSH27_TO_CANONICAL, HEX_EDGES, HEX27_REFERENCE_NODES,
};
use fs_solid::linear3::tet10::{QuadraticTetProblem, TET10_EDGES, TetQuadrature};
use fs_solid::linear3::{TetAssemblyBudget, TetElasticMaterial, TetMaterialField};

/// Gmsh order: vertices, then edges (0,1),(1,2),(2,0),(3,0),(3,2),(3,1).
const GMSH_TET10_REFERENCE: [[f64; 3]; 10] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [0.5, 0.0, 0.0],
    [0.5, 0.5, 0.0],
    [0.0, 0.5, 0.0],
    [0.0, 0.0, 0.5],
    [0.0, 0.5, 0.5],
    [0.5, 0.0, 0.5],
];

/// Gmsh order: corners, edges (0,1),(0,3),(0,4),(1,2),(1,5),(2,3),(2,6),
/// (3,7),(4,5),(4,7),(5,6),(6,7), faces z-,y-,x-,x+,y+,z+, centre.
const GMSH_HEX27_REFERENCE: [[f64; 3]; 27] = [
    [-1.0, -1.0, -1.0],
    [1.0, -1.0, -1.0],
    [1.0, 1.0, -1.0],
    [-1.0, 1.0, -1.0],
    [-1.0, -1.0, 1.0],
    [1.0, -1.0, 1.0],
    [1.0, 1.0, 1.0],
    [-1.0, 1.0, 1.0],
    [0.0, -1.0, -1.0],
    [-1.0, 0.0, -1.0],
    [-1.0, -1.0, 0.0],
    [1.0, 0.0, -1.0],
    [1.0, -1.0, 0.0],
    [0.0, 1.0, -1.0],
    [1.0, 1.0, 0.0],
    [-1.0, 1.0, 0.0],
    [0.0, -1.0, 1.0],
    [-1.0, 0.0, 1.0],
    [1.0, 0.0, 1.0],
    [0.0, 1.0, 1.0],
    [0.0, 0.0, -1.0],
    [0.0, -1.0, 0.0],
    [-1.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [0.0, 0.0, 0.0],
];

/// A sheared, non-dyadic affine map, so coordinates need all 17 digits.
fn affine(p: [f64; 3]) -> [f64; 3] {
    const A: [[f64; 3]; 3] = [
        [0.0731, 0.0113, -0.0047],
        [0.0029, 0.0517, 0.0131],
        [-0.0061, 0.0089, 0.0293],
    ];
    const B: [f64; 3] = [0.1234567, -0.0456789, 0.3141592];
    core::array::from_fn(|i| B[i] + A[i][0] * p[0] + A[i][1] * p[1] + A[i][2] * p[2])
}

fn reference(element_type: GmshElementType) -> Vec<[f64; 3]> {
    match element_type {
        GmshElementType::Tetrahedron10 => GMSH_TET10_REFERENCE.to_vec(),
        GmshElementType::Hexahedron20 => GMSH_HEX27_REFERENCE[..20].to_vec(),
        GmshElementType::Hexahedron27 => GMSH_HEX27_REFERENCE.to_vec(),
        other => panic!("no reference element for {other:?}"),
    }
}

/// One affine element per second-order type, each on its own nodes.
fn high_order_mesh() -> GmshMesh {
    let mut mesh = GmshMesh::new("2.2");
    mesh.physical_names.insert((3, 7), "solid".to_string());
    let mut next_node = 1u64;
    for (id, element_type) in [
        GmshElementType::Tetrahedron10,
        GmshElementType::Hexahedron20,
        GmshElementType::Hexahedron27,
    ]
    .into_iter()
    .enumerate()
    {
        let mut node_ids = Vec::new();
        for p in reference(element_type) {
            mesh.nodes.push(GmshNode {
                id: next_node,
                coords: affine(p),
            });
            node_ids.push(next_node);
            next_node += 1;
        }
        mesh.elements.push(GmshElement {
            id: id as u64 + 1,
            element_type,
            tags: vec![7, 1],
            node_ids,
        });
    }
    mesh
}

fn coordinate_bits(mesh: &GmshMesh) -> Vec<(u64, [u64; 3])> {
    mesh.nodes
        .iter()
        .map(|n| (n.id, n.coords.map(f64::to_bits)))
        .collect()
}

#[test]
fn g0_tet10_hex20_hex27_round_trip() {
    let mesh = high_order_mesh();
    let text = write_gmsh_msh2(&mesh);
    let (quarantined, receipt) = parse_gmsh(&text).expect("second-order solids parse");
    let parsed = quarantined.into_inner();
    assert_eq!(parsed, mesh);
    assert_eq!(coordinate_bits(&parsed), coordinate_bits(&mesh));
    assert_eq!(write_gmsh_msh2(&parsed), text, "rewrite is byte-identical");
    assert_eq!(
        receipt.element_census.into_iter().collect::<Vec<_>>(),
        vec![(11, 1), (12, 1), (17, 1)]
    );
    assert_eq!(receipt.node_count, 10 + 20 + 27);

    // Gmsh's own MSH 2.2 writer prints %.16g; the parse is format-agnostic
    // and the canonical rewrite is then a fixed point.
    let gmsh_text = "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n10\n\
        1 0 0 0\n2 1 0 0\n3 0 1 0\n4 0 0 1\n5 0.5 0 0\n6 0.5 0.5 0\n\
        7 0 0.5 0\n8 0 0 0.5\n9 0 0.5 0.5\n10 0.5 0 0.5\n$EndNodes\n\
        $Elements\n1\n1 11 2 7 1 1 2 3 4 5 6 7 8 9 10\n$EndElements\n";
    let (quarantined, _) = parse_gmsh(gmsh_text).expect("Gmsh-formatted TET10 parses");
    let parsed = quarantined.into_inner();
    for (node, expected) in parsed.nodes.iter().zip(GMSH_TET10_REFERENCE) {
        assert_eq!(node.coords.map(f64::to_bits), expected.map(f64::to_bits));
    }
    assert_eq!(
        parsed.elements[0].element_type,
        GmshElementType::Tetrahedron10
    );
    let canonical = write_gmsh_msh2(&parsed);
    let (again, _) = parse_gmsh(&canonical).expect("canonical rewrite parses");
    assert_eq!(write_gmsh_msh2(&again.into_inner()), canonical);
}

fn close(a: [f64; 3], b: [f64; 3]) -> bool {
    (0..3).all(|k| (a[k] - b[k]).abs() <= 1e-12)
}

fn mean(points: &[[f64; 3]]) -> [f64; 3] {
    core::array::from_fn(|k| points.iter().map(|p| p[k]).sum::<f64>() / points.len() as f64)
}

/// Canonical-order coordinates of each parsed element.
fn canonical_coordinates(mesh: &GmshMesh) -> Vec<(GmshElementType, Vec<[f64; 3]>)> {
    mesh.elements
        .iter()
        .map(|element| {
            let ids = element.canonical_node_ids().expect("canonical order");
            let coords = ids
                .iter()
                .map(|id| {
                    mesh.nodes
                        .iter()
                        .find(|n| n.id == *id)
                        .expect("node listed")
                        .coords
                })
                .collect();
            (element.element_type, coords)
        })
        .collect()
}

#[test]
fn g0_midside_nodes_at_midpoints_pin_ordering() {
    assert_eq!(GMSH_HEX20_TO_CANONICAL, GMSH20_TO_CANONICAL);
    assert_eq!(GMSH_HEX27_TO_CANONICAL, GMSH27_TO_CANONICAL);
    for map in [
        &GMSH_TET10_TO_CANONICAL[..],
        &GMSH_HEX20_TO_CANONICAL[..],
        &GMSH_HEX27_TO_CANONICAL[..],
    ] {
        let mut seen = vec![false; map.len()];
        for &c in map {
            assert!(!seen[c], "permutation repeats {c}");
            seen[c] = true;
        }
    }
    let (quarantined, _) = parse_gmsh(&write_gmsh_msh2(&high_order_mesh())).expect("parse");
    let mesh = quarantined.into_inner();
    for (element_type, x) in canonical_coordinates(&mesh) {
        match element_type {
            GmshElementType::Tetrahedron10 => {
                for (e, [a, b]) in TET10_EDGES.into_iter().enumerate() {
                    assert!(
                        close(x[4 + e], mean(&[x[a], x[b]])),
                        "TET10 canonical node {} is not the midpoint of ({a},{b})",
                        4 + e
                    );
                }
            }
            GmshElementType::Hexahedron20 | GmshElementType::Hexahedron27 => {
                for (a, q) in HEX27_REFERENCE_NODES.iter().take(x.len()).enumerate() {
                    let expected = affine(q.map(f64::from));
                    assert!(
                        close(x[a], expected),
                        "{element_type:?} canonical node {a} is not at reference {q:?}"
                    );
                }
                for (e, [a, b]) in HEX_EDGES.into_iter().enumerate() {
                    let node = HEX27_REFERENCE_NODES
                        .iter()
                        .position(|q| {
                            (0..3).all(|k| {
                                2 * i32::from(q[k])
                                    == i32::from(HEX27_REFERENCE_NODES[a][k])
                                        + i32::from(HEX27_REFERENCE_NODES[b][k])
                            })
                        })
                        .expect("edge node exists");
                    assert!(close(x[node], mean(&[x[a], x[b]])), "hex edge {e}");
                }
                if element_type == GmshElementType::Hexahedron27 {
                    for face in 20..26 {
                        let q = HEX27_REFERENCE_NODES[face];
                        let axis = q.iter().position(|v| *v != 0).expect("face axis");
                        let corners: Vec<[f64; 3]> = (0..8)
                            .filter(|&c| HEX27_REFERENCE_NODES[c][axis] == q[axis])
                            .map(|c| x[c])
                            .collect();
                        assert_eq!(corners.len(), 4);
                        assert!(close(x[face], mean(&corners)), "hex face node {face}");
                    }
                    assert!(close(x[26], mean(&x[..8])), "hex body node");
                }
            }
            other => panic!("unexpected {other:?}"),
        }
    }
    // Negative control: reading Gmsh order as canonical misplaces the TET10
    // midside nodes of edges (1,3) and (2,3), so the check above discriminates.
    let tet = &mesh.elements[0];
    let raw: Vec<[f64; 3]> = tet
        .node_ids
        .iter()
        .map(|id| mesh.nodes[(*id - 1) as usize].coords)
        .collect();
    let misplaced: Vec<usize> = TET10_EDGES
        .iter()
        .enumerate()
        .filter(|(e, [a, b])| !close(raw[4 + e], mean(&[raw[*a], raw[*b]])))
        .map(|(e, _)| 4 + e)
        .collect();
    assert_eq!(misplaced, vec![8, 9]);
}

fn single_element_text(type_id: u32, nodes: usize) -> String {
    let mut text = String::from("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n");
    text.push_str(&format!("{nodes}\n"));
    for i in 1..=nodes {
        text.push_str(&format!("{i} {} 0 0\n", i as f64 * 0.1));
    }
    text.push_str("$EndNodes\n$Elements\n1\n");
    text.push_str(&format!("1 {type_id} 2 1 1"));
    for i in 1..=nodes {
        text.push_str(&format!(" {i}"));
    }
    text.push_str("\n$EndElements\n");
    text
}

#[test]
fn g0_unknown_type_refuses() {
    // Second-order surfaces, prisms, pyramids, serendipity faces, third-order
    // tets and nonsense ids: none has an implemented node order.
    for type_id in [0, 8, 9, 10, 13, 14, 16, 18, 19, 29, 92, 255] {
        let result = parse_gmsh(&single_element_text(type_id, 30));
        assert!(
            matches!(result, Err(IoError::Unsupported { .. })),
            "type {type_id} must refuse, got {result:?}"
        );
    }
    for (type_id, nodes) in [(11, 10), (17, 20), (12, 27)] {
        let (q, _) = parse_gmsh(&single_element_text(type_id, nodes)).expect("supported");
        assert_eq!(q.into_inner().elements[0].node_ids.len(), nodes);
        // One node short is malformed, not silently padded.
        let short = single_element_text(type_id, nodes - 1);
        assert!(matches!(parse_gmsh(&short), Err(IoError::Malformed { .. })));
    }
    let surface = GmshElement {
        id: 1,
        element_type: GmshElementType::Triangle,
        tags: vec![],
        node_ids: vec![1, 2, 3],
    };
    assert!(matches!(
        surface.canonical_node_ids(),
        Err(IoError::Unsupported { .. })
    ));
    let truncated = GmshElement {
        id: 2,
        element_type: GmshElementType::Hexahedron27,
        tags: vec![],
        node_ids: (1..=20).collect(),
    };
    assert!(matches!(
        truncated.canonical_node_ids(),
        Err(IoError::Malformed { .. })
    ));
}

#[test]
fn g0_element_cap_and_cap_plus_one() {
    let text = write_gmsh_msh2(&high_order_mesh());
    let at_cap = GmshLimits {
        max_nodes: 57,
        max_elements: 3,
        max_bytes: text.len(),
    };
    let (q, receipt) = parse_gmsh_with_limits(&text, &at_cap).expect("exactly at every cap");
    assert_eq!((receipt.node_count, receipt.element_count), (57, 3));
    assert_eq!(q.into_inner().elements.len(), 3);
    for limits in [
        GmshLimits {
            max_elements: 2,
            ..at_cap.clone()
        },
        GmshLimits {
            max_nodes: 56,
            ..at_cap.clone()
        },
        GmshLimits {
            max_bytes: text.len() - 1,
            ..at_cap.clone()
        },
    ] {
        assert!(
            matches!(
                parse_gmsh_with_limits(&text, &limits),
                Err(IoError::ResourceBound { .. })
            ),
            "cap + 1 must refuse under {limits:?}"
        );
    }
}

// ---- The committed gmsh/OpenCASCADE fixture (crates/fs-io/tests/fixtures/gmsh,
// produced offline by make_holed_block.py; the receipt names tool and version).

const HOLED_BLOCK_TET10: &str = include_str!("fixtures/gmsh/holed_block_tet10.msh");
const ADAPTER_RECEIPT: &str = include_str!("fixtures/gmsh/holed_block.receipt.txt");

fn adapter_receipt() -> std::collections::BTreeMap<&'static str, &'static str> {
    ADAPTER_RECEIPT
        .lines()
        .filter_map(|line| line.split_once('='))
        .collect()
}

fn receipt_number(key: &str) -> f64 {
    adapter_receipt()[key]
        .parse()
        .unwrap_or_else(|_| panic!("receipt {key} is a number"))
}

fn millimetres() -> GmshSolidDeclaration {
    GmshSolidDeclaration {
        unit: Some(LengthUnit {
            label: "mm".to_string(),
            metres_per_unit: 1e-3,
        }),
        reference_state: Some(ReferenceState::Unloaded),
    }
}

fn with_cx<R>(f: impl FnOnce(&fs_exec::Cx<'_>) -> R) -> R {
    let gate = fs_exec::CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = fs_exec::Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 0x63,
                kernel_id: 63,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn tet_min_dihedral_deg(p: [[f64; 3]; 4]) -> f64 {
    let sub = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let cross = |a: [f64; 3], b: [f64; 3]| {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    };
    let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let faces = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]];
    let normals: Vec<[f64; 3]> = faces
        .iter()
        .enumerate()
        .map(|(k, f)| {
            let n = cross(sub(p[f[1]], p[f[0]]), sub(p[f[2]], p[f[0]]));
            if dot(n, sub(p[k], p[f[0]])) > 0.0 {
                n.map(|v| -v)
            } else {
                n
            }
        })
        .collect();
    let mut best = 180.0f64;
    for k in 0..4 {
        for l in k + 1..4 {
            let c = dot(normals[k], normals[l])
                / (dot(normals[k], normals[k]) * dot(normals[l], normals[l])).sqrt();
            best = best.min(180.0 - c.clamp(-1.0, 1.0).acos().to_degrees());
        }
    }
    best
}

#[test]
fn g1_imported_mesh_promotes_and_assembles() {
    let receipt = adapter_receipt();
    assert_eq!(receipt["adapter.tool"], "gmsh");
    assert_eq!(receipt["adapter.kernel"], "OpenCASCADE");
    assert!(!receipt["adapter.version"].is_empty());
    assert!(receipt["authority"].starts_with("Estimate"));
    assert_eq!(receipt["unit"], "mm");
    assert_eq!(
        HOLED_BLOCK_TET10.len(),
        receipt["tet10.bytes"].parse::<usize>().unwrap()
    );
    let elements: usize = receipt["tet10.elements"].parse().unwrap();
    let nodes: usize = receipt["tet10.nodes"].parse().unwrap();

    // Quarantined at parse: the census is the receipt's.
    let (quarantined, census) = parse_gmsh(HOLED_BLOCK_TET10).expect("adapter mesh parses");
    assert_eq!(
        census.element_census.iter().collect::<Vec<_>>(),
        vec![(&11u8, &elements)]
    );
    let (solid, promoted) = promote_gmsh_solid(quarantined, &millimetres()).expect("promotes");
    assert_eq!(promoted.element_type, 11);
    assert_eq!((promoted.elements, promoted.nodes), (elements, nodes));
    assert_eq!(promoted.reference_state, ReferenceState::Unloaded);
    assert_eq!(
        solid.region_names.get(&1).map(String::as_str),
        Some("block")
    );

    // Assembled by TET10: the element volumes are gmsh's curved volumes
    // (its own Jacobians on an exact rule), and the mass is rho V.
    let tet10: Vec<[usize; 10]> = solid
        .elements
        .iter()
        .map(|e| e.as_slice().try_into().expect("ten nodes"))
        .collect();
    let density = 1000.0;
    let material = TetElasticMaterial::try_new(density, 2.0e9, 0.3, [0x63; 32]).unwrap();
    let assembly = with_cx(|cx| {
        QuadraticTetProblem {
            nodes_m: &solid.nodes_m,
            elements: &tet10,
            materials: TetMaterialField::Uniform(&material),
            fixed_dofs: &[],
            budget: TetAssemblyBudget::standard(),
            quadrature: TetQuadrature {
                stiffness_points_per_dir: 3,
                mass_points_per_dir: 4,
            },
        }
        .assemble(cx)
    })
    .expect("TET10 assembles the imported mesh");
    let mm3 = 1e-9;
    let rel = |a: f64, b: f64| ((a - b) / b).abs();
    let volume: f64 = assembly.element_volumes_m3.iter().sum();
    assert!(rel(volume, receipt_number("tet10.volume_mm3") * mm3) <= 1e-12);
    let smallest = assembly
        .element_volumes_m3
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min);
    let largest = assembly
        .element_volumes_m3
        .iter()
        .copied()
        .fold(0.0, f64::max);
    assert!(
        rel(
            smallest,
            receipt_number("tet10.element_volume_min_mm3") * mm3
        ) <= 1e-12
    );
    assert!(
        rel(
            largest,
            receipt_number("tet10.element_volume_max_mm3") * mm3
        ) <= 1e-12
    );
    assert!(rel(assembly.total_mass_kg, density * volume) <= 1e-12);
    assert!(assembly.minimum_scaled_jacobian > 0.0);
    // The curved volume is closer to the CAD solid than the straight one.
    let cad = receipt_number("cad.volume_mm3") * mm3;
    let straight = receipt_number("tet10.straight_sided_volume_mm3") * mm3;
    assert!((volume - cad).abs() < 0.01 * (straight - cad).abs());
    // Quality census: corner-tet dihedrals in metres match the adapter's.
    let dihedral = tet10
        .iter()
        .map(|e| tet_min_dihedral_deg(core::array::from_fn(|k| solid.nodes_m[e[k]])))
        .fold(180.0, f64::min);
    assert!((dihedral - receipt_number("tet10.min_corner_dihedral_deg")).abs() <= 1e-9);
    eprintln!(
        "{{\"suite\":\"irp-s1\",\"bead\":\"M1.3-W4\",\"test\":\"g1_imported_mesh_promotes_and_assembles\",\
         \"inputs\":{{\"source\":\"gmsh-adapter\",\"adapter\":{{\"tool\":\"{}\",\"version\":\"{}\"}},\
         \"fixture\":\"crates/fs-io/tests/fixtures/gmsh/holed_block_tet10.msh\"}},\
         \"observed\":{{\"elements\":{elements},\"nodes\":{nodes},\"volume_m3\":{volume:e},\
         \"receipt_volume_m3\":{:e},\"cad_volume_m3\":{cad:e},\"min_corner_dihedral_deg\":{dihedral},\
         \"min_scaled_jacobian\":{:e}}},\"tolerance\":1e-12,\"verdict\":\"pass\"}}",
        receipt["adapter.tool"],
        receipt["adapter.version"],
        receipt_number("tet10.volume_mm3") * mm3,
        assembly.minimum_scaled_jacobian
    );
}

#[test]
fn g0_promotion_refuses_undeclared_and_unsupported() {
    let promote = |text: &str, d: &GmshSolidDeclaration| {
        let (q, _) = parse_gmsh(text).expect("parses");
        promote_gmsh_solid(q, d).map(|_| ())
    };
    let mut no_unit = millimetres();
    no_unit.unit = None;
    assert_eq!(
        promote(HOLED_BLOCK_TET10, &no_unit),
        Err(GmshPromotionRefusal::UnitMissing)
    );
    let mut zero_unit = millimetres();
    zero_unit.unit = Some(LengthUnit {
        label: "mm".to_string(),
        metres_per_unit: 0.0,
    });
    assert_eq!(
        promote(HOLED_BLOCK_TET10, &zero_unit),
        Err(GmshPromotionRefusal::UnitMissing)
    );
    let mut no_state = millimetres();
    no_state.reference_state = None;
    assert_eq!(
        promote(HOLED_BLOCK_TET10, &no_state),
        Err(GmshPromotionRefusal::ReferenceStateMissing)
    );

    let mesh = high_order_mesh();
    // Three volume families in one file.
    assert!(matches!(
        promote(&write_gmsh_msh2(&mesh), &millimetres()),
        Err(GmshPromotionRefusal::MixedVolumeTypes {
            first: 11,
            second: 17
        })
    ));
    let mut tet_only = mesh.clone();
    tet_only.elements.truncate(1);
    tet_only.nodes.truncate(10);
    let (q, _) = parse_gmsh(&write_gmsh_msh2(&tet_only)).unwrap();
    let (solid, receipt) = promote_gmsh_solid(q, &millimetres()).expect("one TET10");
    assert_eq!(receipt.unreferenced_nodes, 0);
    // Canonical order and unit applied: the midside of canonical edge (1,3)
    // is node 8, in metres.
    let x = |i: usize| solid.nodes_m[solid.elements[0][i]];
    for k in 0..3 {
        assert!((x(8)[k] - (x(1)[k] + x(3)[k]) / 2.0).abs() <= 1e-15);
        assert!((x(0)[k] - affine(GMSH_TET10_REFERENCE[0])[k] * 1e-3).abs() <= 1e-18);
    }
    let mut prism = tet_only.clone();
    prism.elements.push(GmshElement {
        id: 9,
        element_type: GmshElementType::Prism,
        tags: vec![7, 1],
        node_ids: vec![1, 2, 3, 4, 5, 6],
    });
    assert!(matches!(
        promote(&write_gmsh_msh2(&prism), &millimetres()),
        Err(GmshPromotionRefusal::NoCanonicalOrder {
            element: 9,
            type_id: 6
        })
    ));
    let mut surface_only = tet_only.clone();
    surface_only.elements = vec![GmshElement {
        id: 1,
        element_type: GmshElementType::Triangle,
        tags: vec![7, 1],
        node_ids: vec![1, 2, 3],
    }];
    assert_eq!(
        promote(&write_gmsh_msh2(&surface_only), &millimetres()),
        Err(GmshPromotionRefusal::NoVolumeElements)
    );
    let mut dangling = tet_only.clone();
    dangling.elements[0].node_ids[9] = 99;
    assert!(matches!(
        promote(&write_gmsh_msh2(&dangling), &millimetres()),
        Err(GmshPromotionRefusal::UnknownNode { node: 99, .. })
    ));
    let mut repeated = tet_only.clone();
    repeated.elements[0].node_ids[9] = 1;
    assert!(matches!(
        promote(&write_gmsh_msh2(&repeated), &millimetres()),
        Err(GmshPromotionRefusal::RepeatedNode { .. })
    ));
    let mut duplicate = tet_only.clone();
    duplicate.nodes[1].id = 1;
    assert!(matches!(
        promote(&write_gmsh_msh2(&duplicate), &millimetres()),
        Err(GmshPromotionRefusal::DuplicateNodeId { id: 1 })
    ));
}
