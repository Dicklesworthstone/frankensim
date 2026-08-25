//! No-mock multi-region pipeline proof (bead frankensim-s93ej.1).
//!
//! The committed two-solid reference project
//! (`data/reference-project/multi-region-interface.fsim`) parses through the
//! real wire grammar, both meshes admit through the real STL reader, every
//! declared assignment resolves through the real project assignment surface,
//! the production PLC mesher volumetricizes and independently audits the
//! shared-interface pair, and the real conduction consumer opens the audited
//! labeled complex.

use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_io::read_stl;
use fs_mesh::{RegionId, RegionKind, RegionSpec, UnverifiedPlc, VolumetricPolicy, volumetricize};
use fs_project::ImportedMeshLibrary;

const DATA: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../data/reference-project");

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x7E7,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

/// Independent signed tet volume from retained positions (auditor formula).
fn tet_volume(positions: &[[f64; 3]], tet: &[u32; 4]) -> f64 {
    let p = |i: u32| positions[i as usize];
    let (a, b, c, d) = (p(tet[0]), p(tet[1]), p(tet[2]), p(tet[3]));
    let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let ac = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
    let ad = [d[0] - a[0], d[1] - a[1], d[2] - a[2]];
    let cross = [
        ab[1] * ac[2] - ab[2] * ac[1],
        ab[2] * ac[0] - ab[0] * ac[2],
        ab[0] * ac[1] - ab[1] * ac[0],
    ];
    (cross[0] * ad[0] + cross[1] * ad[1] + cross[2] * ad[2]) / 6.0
}

#[test]
fn multi_region_fixture_resolves_volumetricizes_and_opens() {
    // -- project parses through the real wire grammar --------------------
    let project_path = format!("{DATA}/multi-region-interface.fsim");
    let src = std::fs::read_to_string(&project_path).expect("committed fixture present");
    let decoded = fs_project::parse_sexpr(&src).expect("canonical parse");
    let artifacts = decoded
        .spec
        .geometry
        .clone()
        .expect("multi-region fixture declares geometry rows");
    assert_eq!(artifacts.len(), 2, "two solids share one interface");

    // -- both meshes admit through the real STL reader -------------------
    // Raw geometry is snapshotted before the Soups move into the imported
    // library: the PLC below needs positions/triangles independent of
    // library ownership, and Soups are not promised Clone.
    let mut library = ImportedMeshLibrary::new();
    let mut raw = Vec::new();
    for artifact in &artifacts {
        let mesh_path = format!("{DATA}/multi-region-{role}.stl", role = artifact.role);
        let bytes = std::fs::read(&mesh_path)
            .unwrap_or_else(|e| panic!("committed mesh missing: {mesh_path}: {e}"));
        let soup: fs_rep_mesh::Soup = fs_io::stl::read_stl(&bytes).expect("stl admits");
        assert_eq!(soup.triangles.len(), 12, "unit cube fixture topology");
        raw.push((
            soup.positions.clone(),
            soup.triangles.clone(),
        ));
    }
    for (artifact, (positions, triangles)) in artifacts.iter().zip(raw.iter()) {
        let soup = fs_rep_mesh::Soup {
            positions: positions.clone(),
            triangles: triangles.clone(),
        };
        library.insert(artifact, soup, "m", Vec::new());
    }

    with_cx(|cx| {
        // -- assignments resolve against the imported library ------------
        let resolution = fs_project::resolve_geometry_assignments(
            &decoded.spec,
            &library,
            fs_io::AssignmentLimits::default(),
            cx,
        );
        assert!(
            resolution.admissible(),
            "assignment violations: {:?}",
            resolution.violations
        );
        assert_eq!(resolution.artifacts.len(), 2);

        // -- production PLC: merged table, two solid region specs --------
        let cold = &raw[0];
        let hot = &raw[1];
        assert_eq!(cold.0.len(), 8, "cube vertex count");
        assert_eq!(hot.0.len(), 8, "cube vertex count");
        let mut vertices = cold.0.clone();
        let vertex_offset = vertices.len() as u32;
        vertices.extend(hot.0.iter().copied());

        let regions = vec![
            RegionSpec {
                id: RegionId(1),
                kind: RegionKind::Solid,
                seed: [0.5, 0.5, 0.5],
                triangles: cold.1.clone(),
            },
            RegionSpec {
                id: RegionId(2),
                kind: RegionKind::Solid,
                seed: [1.5, 0.5, 0.5],
                triangles: hot.1
                    .iter()
                    .map(|t| [t[0] + vertex_offset, t[1] + vertex_offset, t[2] + vertex_offset])
                    .collect(),
            },
        ];

        // -- production mesher: the crate's own audited entry point ------
        let plc = UnverifiedPlc::new(vertices, regions);
        let policy = VolumetricPolicy::fixture_default("m");
        let audited = fs_mesh::volumetricize(plc, policy, cx)
            .expect("production mesher volumetricizes and audits");

        // -- independent volume recomputation from retained bytes --------
        let labeled = audited.labeled();
        let positions = labeled.positions();
        let mut per_region = std::collections::BTreeMap::new();
        for (tet, region) in labeled.tets().iter().zip(labeled.region_of_tet()) {
            *per_region.entry(region.0).or_insert(0.0) += tet_volume(positions, tet);
        }
        assert_eq!(per_region.len(), 2, "both declared regions retain tets");
        for id in [1u32, 2] {
            let volume = per_region[&id];
            assert!(
                (volume - 1.0).abs() < 1e-9,
                "region {id} volume {volume} != analytic unit-cube 1.0"
            );
        }
        assert!(
            labeled.region_of_tet().len() == labeled.tets().len(),
            "every retained tet carries exactly one region label"
        );

        // -- the real conduction consumer opens the audited complex ------
        let complex =
            fs_rep_mesh::TetComplex::from_tets(positions.len(), labeled.tets().to_vec());
        assert_eq!(complex.vertex_count, positions.len());
        let mesh = fs_conduction::ConductionMesh::new(complex, positions.clone())
            .expect("conduction consumer opens the audited labeled complex");
        let _ = mesh;
    });
}
