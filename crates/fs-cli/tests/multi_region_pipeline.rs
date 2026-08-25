//! No-mock multi-region pipeline proof (bead frankensim-s93ej.1).
//!
//! The committed two-solid reference project
//! (`data/reference-project/multi-region-interface.fsim`) parses through the
//! real wire grammar, both meshes admit through the real STL reader, every
//! declared assignment resolves through the real project assignment surface,
//! the production PLC mesher volumetricizes and independently audits the
//! shared-interface pair. The real conduction consumer then derives
//! region-owned P1 traces from that same audited labeled complex and the real
//! project adapter resolves every coincident trace back to the declared,
//! oriented interface and retained source faces.

use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::Point3;
use fs_mesh::{RegionId, RegionKind, RegionSpec, UnverifiedPlc, VolumetricPolicy};
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

/// Exact-bit vertex welding across the two committed solids, mirroring
/// the proven adjacent-solids idiom in fs-mesh's own volumetric battery:
/// the shared x=1 plane must become ONE vertex ring, not duplicates.
fn weld(parts: &[Vec<[f64; 3]>]) -> (Vec<[f64; 3]>, Vec<Vec<u32>>) {
    let mut verts = Vec::new();
    let mut index = std::collections::BTreeMap::new();
    let mut remaps = Vec::new();
    for part in parts {
        let mut remap = Vec::with_capacity(part.len());
        for p in part {
            let key = [p[0].to_bits(), p[1].to_bits(), p[2].to_bits()];
            let id = if let Some(&existing) = index.get(&key) {
                existing
            } else {
                let id = u32::try_from(verts.len()).expect("vertex count");
                verts.push(*p);
                index.insert(key, id);
                id
            };
            remap.push(id);
        }
        remaps.push(remap);
    }
    (verts, remaps)
}

/// Retention + independent recheck: the audited labeled complex is
/// serialized to canonical JSONL with a domain-separated content root,
/// dropped to disk, reloaded from those bytes alone, and re-verified by
/// recomputation - the retained artifact must stand on its own without
/// the producing process in the room.
#[test]
fn multi_region_artifacts_survive_retention_and_recheck() {
    use fs_blake3::hash_domain;

    const RETENTION_DOMAIN: &str = "org.frankensim.fs-cli.tests.multi-region-retention.v1";
    let project_path = format!("{DATA}/multi-region-interface.fsim");
    let src = std::fs::read_to_string(&project_path).expect("fixture present");
    let decoded = fs_project::parse_sexpr(&src).expect("parse");
    let artifacts = decoded.spec.geometry.clone().expect("geometry rows");

    let mut raw = Vec::new();
    for artifact in &artifacts {
        let mesh_path = format!("{DATA}/multi-region-{role}.stl", role = artifact.role);
        let soup: fs_rep_mesh::Soup =
            fs_io::stl::read_stl(&std::fs::read(&mesh_path).expect("mesh bytes")).expect("stl");
        raw.push((
            soup.positions
                .iter()
                .map(|p| [p.x, p.y, p.z])
                .collect::<Vec<[f64; 3]>>(),
            soup.triangles.clone(),
        ));
    }

    with_cx(|cx| {
        let mut library = ImportedMeshLibrary::new();
        for (artifact, (positions, triangles)) in artifacts.iter().zip(raw.iter()) {
            let soup = fs_rep_mesh::Soup {
                positions: positions
                    .iter()
                    .map(|p| Point3 {
                        x: p[0],
                        y: p[1],
                        z: p[2],
                    })
                    .collect(),
                triangles: triangles.clone(),
            };
            library.insert(artifact, soup, "m", Vec::new());
        }
        let resolution = fs_project::resolve_geometry_assignments(
            &decoded.spec,
            &library,
            fs_io::AssignmentLimits::default(),
            cx,
        );
        assert!(resolution.admissible());

        let (vertices, remaps) = weld(&[raw[0].0.clone(), raw[1].0.clone()]);
        let remap_tris = |tris: &[[u32; 3]], remap: &[u32]| -> Vec<[u32; 3]> {
            tris.iter()
                .map(|t| {
                    [
                        remap[t[0] as usize],
                        remap[t[1] as usize],
                        remap[t[2] as usize],
                    ]
                })
                .collect()
        };
        let regions = vec![
            RegionSpec {
                id: RegionId(1),
                kind: RegionKind::Solid,
                seed: [0.5, 0.5, 0.5],
                triangles: remap_tris(&raw[0].1, &remaps[0]),
            },
            RegionSpec {
                id: RegionId(2),
                kind: RegionKind::Solid,
                seed: [1.5, 0.5, 0.5],
                triangles: remap_tris(&raw[1].1, &remaps[1]),
            },
        ];
        let plc = UnverifiedPlc::new(vertices, regions);
        let audited = fs_mesh::volumetricize(plc, VolumetricPolicy::fixture_default("m"), cx)
            .expect("production mesher");
        let labeled = audited.labeled();

        // -- retain: canonical JSONL + content root ----------------------
        let mut lines = Vec::new();
        lines.push(format!(
            "{{\"kind\":\"header\",\"length_unit\":\"{}\"}}",
            "m"
        ));
        for (tet, region) in labeled.tets().iter().zip(labeled.region_of_tet()) {
            lines.push(format!(
                "{{\"kind\":\"tet\",\"v\":[{},{},{},{}],\"region\":{}}}",
                tet[0], tet[1], tet[2], tet[3], region.0
            ));
        }
        for p in labeled.positions() {
            lines.push(format!(
                "{{\"kind\":\"position\",\"p\":[{:?},{:?},{:?}]}}",
                p[0], p[1], p[2]
            ));
        }
        let payload = lines.join("\n");
        let root = hash_domain(RETENTION_DOMAIN, payload.as_bytes());
        let root_hex: String = root.0.iter().map(|b| format!("{b:02x}")).collect();
        let dir = std::env::temp_dir().join("mri-retention");
        std::fs::create_dir_all(&dir).expect("artifact dir");
        let artifact = dir.join("labeled_complex.jsonl");
        std::fs::write(&artifact, &payload).expect("retain");

        // -- independent recheck from disk bytes alone -------------------
        let reloaded = std::fs::read_to_string(&artifact).expect("retained artifact reads");
        assert_eq!(
            hash_domain(RETENTION_DOMAIN, reloaded.as_bytes()).0,
            root.0,
            "content root survives the round trip"
        );
        let mut reload_positions = Vec::<[f64; 3]>::new();
        let mut reload_tets = Vec::<([u32; 4], u32)>::new();
        // Deterministic parse: split on known key markers rather than a
        // full JSON parser (this test has no serde dependency).
        for line in reloaded.lines() {
            let grab = |key: &str| -> Vec<u32> {
                let k = format!("\"{key}\":[");
                let i = line.find(&k).expect("key present") + k.len();
                let end = line[i..].find(']').expect("array close") + i;
                line[i..end]
                    .split(',')
                    .filter_map(|t| t.trim().parse::<u32>().ok())
                    .collect()
            };
            if line.contains("\"kind\":\"tet\"") {
                let v = grab("v");
                let rk = "\"region\":";
                let ri = line.find(rk).expect("region key") + rk.len();
                let rend = line[ri..]
                    .find(|c: char| !c.is_ascii_digit())
                    .unwrap_or(line.len() - ri);
                let r: u32 = line[ri..ri + rend].parse().expect("region id");
                reload_tets.push(([v[0], v[1], v[2], v[3]], r));
            } else if line.contains("\"kind\":\"position\"") {
                let kf = "\"p\":[";
                let i = line.find(kf).expect("p key") + kf.len();
                let end = line[i..].find(']').expect("close") + i;
                let xyz: Vec<f64> = line[i..end]
                    .split(',')
                    .filter_map(|t| t.trim().parse::<f64>().ok())
                    .collect();
                reload_positions.push([xyz[0], xyz[1], xyz[2]]);
            }
        }
        assert_eq!(reload_positions.len(), labeled.positions().len());
        assert_eq!(reload_tets.len(), labeled.tets().len());

        // -- recompute volumes purely from the reloaded artifact ---------
        let mut per_region = std::collections::BTreeMap::new();
        for (tet, region) in &reload_tets {
            *per_region.entry(*region).or_insert(0.0) += tet_volume(&reload_positions, tet).abs();
        }
        assert_eq!(per_region.len(), 2);
        for id in [1u32, 2] {
            assert!(
                (per_region[&id] - 1.0).abs() < 1e-9,
                "reloaded region {id} volume {} != 1.0",
                per_region[&id]
            );
        }
    });
}

#[test]
fn multi_region_fixture_resolves_volumetricizes_and_binds_contact_traces() {
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
    let mut soups = Vec::new();
    let mut raw = Vec::new();
    for artifact in &artifacts {
        let mesh_path = format!("{DATA}/multi-region-{role}.stl", role = artifact.role);
        let bytes = std::fs::read(&mesh_path)
            .unwrap_or_else(|e| panic!("committed mesh missing: {mesh_path}: {e}"));
        let soup: fs_rep_mesh::Soup = fs_io::stl::read_stl(&bytes).expect("stl admits");
        assert_eq!(soup.triangles.len(), 12, "unit cube fixture topology");
        raw.push((
            soup.positions
                .iter()
                .map(|p| [p.x, p.y, p.z])
                .collect::<Vec<[f64; 3]>>(),
            soup.triangles.clone(),
        ));
        soups.push(soup);
    }
    for (artifact, soup) in artifacts.iter().zip(soups) {
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

        // -- production PLC: exact-welded table, two solid region specs --
        let (vertices, remaps) = weld(&[raw[0].0.clone(), raw[1].0.clone()]);
        let remap_tris = |tris: &[[u32; 3]], remap: &[u32]| -> Vec<[u32; 3]> {
            tris.iter()
                .map(|t| {
                    [
                        remap[t[0] as usize],
                        remap[t[1] as usize],
                        remap[t[2] as usize],
                    ]
                })
                .collect()
        };
        let regions = vec![
            RegionSpec {
                id: RegionId(1),
                kind: RegionKind::Solid,
                seed: [0.5, 0.5, 0.5],
                triangles: remap_tris(&raw[0].1, &remaps[0]),
            },
            RegionSpec {
                id: RegionId(2),
                kind: RegionKind::Solid,
                seed: [1.5, 0.5, 0.5],
                triangles: remap_tris(&raw[1].1, &remaps[1]),
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
        // The producer's signed-tet convention fixes an orientation this
        // file does not own; the independent check is MAGNITUDE against
        // the analytic unit cube, plus agreement with the auditor's
        // producer-side accumulation.
        let witness = audited.witness().per_region_producer.clone();
        for id in [1u32, 2] {
            let volume = per_region[&id].abs();
            assert!(
                (volume - 1.0).abs() < 1e-9,
                "region {id} volume {volume} != analytic unit-cube 1.0"
            );
            let producer = witness
                .iter()
                .find(|(rid, _)| rid.0 == id)
                .map(|(_, v)| *v)
                .unwrap_or(f64::NAN);
            assert!(
                (producer.abs() - volume).abs() < 1e-9,
                "witness/auditor disagreement on region {id}: {producer} vs {volume}"
            );
        }
        assert!(
            labeled.region_of_tet().len() == labeled.tets().len(),
            "every retained tet carries exactly one region label"
        );

        // -- real contact topology + declared interface lowering ----------
        let complex = fs_rep_mesh::TetComplex::from_tets(positions.len(), labeled.tets().to_vec());
        assert_eq!(complex.vertex_count, positions.len());
        let labels = labeled
            .region_of_tet()
            .iter()
            .map(|region| region.0)
            .collect::<Vec<_>>();
        let mesh =
            fs_conduction::ConductionMesh::new_region_owned(complex, positions.to_vec(), &labels)
                .expect("conduction consumer derives region-owned traces from the audited volume");
        let candidates = fs_conduction::ThermalInterfaces::coincident_face_pairs(&mesh)
            .expect("region-owned mesh has valid matching-P1 candidates");
        assert_eq!(
            candidates.len(),
            2,
            "the square joint has two trace triangles"
        );

        let interfaces = fs_project::resolve_conduction_interface_pairs(
            &decoded.spec,
            &library,
            fs_io::AssignmentLimits::default(),
            fs_project::ConductionInterfaceLimits::DEFAULT,
            &mesh,
            cx,
        );
        assert!(
            interfaces.admissible(),
            "interface-lowering violations: {:?}",
            interfaces.violations
        );
        assert_eq!(interfaces.pairs.len(), candidates.len());
        for pair in interfaces.pairs {
            assert_eq!(pair.interface, "cold-hot-joint");
            assert_eq!(pair.from_region, "cold");
            assert_eq!(pair.to_region, "hot");
            assert!(
                !pair.interface_sources.is_empty(),
                "the declared interface selector must retain its source face"
            );
        }
    });
}
