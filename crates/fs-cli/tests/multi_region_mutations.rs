//! Mutation catalog for the multi-region no-mock journey (bead
//! frankensim-s93ej.1, tail item c). Each case perturbs the committed
//! two-solid fixture along one axis named in the acceptance matrix and
//! pins the expected behavior: benign transforms preserve per-region
//! analytic volumes, corrupt geometry refuses with typed errors.

use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_mesh::volumetric::{RegionId, RegionKind, RegionSpec, UnverifiedPlc, VolumetricPolicy};
use fs_project::ImportedMeshLibrary;

const DATA: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../data/reference-project");

struct RawMesh {
    positions: Vec<[f64; 3]>,
    triangles: Vec<[u32; 3]>,
}

fn load_committed() -> Vec<RawMesh> {
    let project_path = format!("{DATA}/multi-region-interface.fsim");
    let src = std::fs::read_to_string(&project_path).expect("fixture present");
    let decoded = fs_project::parse_sexpr(&src).expect("parse");
    let artifacts = decoded.spec.geometry.clone().expect("geometry rows");
    let mut raw = Vec::new();
    for artifact in &artifacts {
        let mesh_path = format!("{DATA}/multi-region-{role}.stl", role = artifact.role);
        let bytes = std::fs::read(&mesh_path).expect("mesh bytes");
        let soup: fs_rep_mesh::Soup = fs_io::stl::read_stl(&bytes).expect("stl admits");
        raw.push(RawMesh {
            positions: soup
                .positions
                .iter()
                .map(|p| [p.x, p.y, p.z])
                .collect::<Vec<[f64; 3]>>(),
            triangles: soup.triangles.clone(),
        });
    }
    assert_eq!(raw.len(), 2);
    raw
}

/// Exact-bit welding of the two solids' vertex tables (shared x=1 ring
/// becomes one vertex set), mirroring fs-mesh's adjacent-solids idiom.
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

fn region_specs(
    raw: &[RawMesh],
    remaps: &[Vec<u32>],
    seeds: &[[f64; 3]],
) -> Vec<RegionSpec> {
    raw.iter()
        .enumerate()
        .map(|(i, mesh)| RegionSpec {
            id: RegionId(i as u32 + 1),
            kind: RegionKind::Solid,
            seed: seeds[i],
            triangles: remap_tris(&mesh.triangles, &remaps[i]),
        })
        .collect()
}

fn remap_tris(tris: &[[u32; 3]], remap: &[u32]) -> Vec<[u32; 3]> {
    tris.iter()
        .map(|t| {
            [
                remap[t[0] as usize],
                remap[t[1] as usize],
                remap[t[2] as usize],
            ]
        })
        .collect()
}

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

/// Volumetricize and return per-region |volume| recomputed independently
/// from the retained labeled complex.
fn volumes_for(raw: &[RawMesh], seeds: &[[f64; 3]]) -> std::collections::BTreeMap<u32, f64> {
    with_cx(|cx| {
        let (vertices, remaps) = weld(&[raw[0].positions.clone(), raw[1].positions.clone()]);
        let regions = region_specs(raw, &remaps, seeds);
        let plc = UnverifiedPlc::new(vertices, regions);
        let audited =
            fs_mesh::volumetricize(plc, VolumetricPolicy::fixture_default("m"), cx)
                .expect("production mesher handles the mutated fixture");
        let labeled = audited.labeled();
        let positions = labeled.positions();
        let mut per_region = std::collections::BTreeMap::new();
        for (tet, region) in labeled.tets().iter().zip(labeled.region_of_tet()) {
            *per_region.entry(region.0).or_insert(0.0) +=
                signed_tet_volume(positions, tet).abs();
        }
        per_region
    })
}

fn signed_tet_volume(positions: &[[f64; 3]], tet: &[u32; 4]) -> f64 {
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

fn assert_unit_cubes(volumes: &std::collections::BTreeMap<u32, f64>) {
    assert_eq!(volumes.len(), 2, "both regions retain tets");
    for id in [1u32, 2] {
        let v = volumes[&id];
        assert!(
            (v - 1.0).abs() < 1e-9,
            "region {id} volume {v} != analytic 1.0"
        );
    }
}

#[test]
fn triangle_permutation_is_volume_invariant() {
    let mut raw = load_committed();
    for mesh in &mut raw {
        mesh.triangles.reverse();
    }
    let volumes = volumes_for(&raw, &[[0.5, 0.5, 0.5], [1.5, 0.5, 0.5]]);
    assert_unit_cubes(&volumes);
}

#[test]
fn vertex_permutation_with_remap_is_volume_invariant() {
    let mut raw = load_committed();
    for mesh in &mut raw {
        // Reverse vertex order and remap every triangle index to match:
        // identical geometry, different encoding.
        let n = mesh.positions.len();
        let remap: std::collections::BTreeMap<usize, usize> =
            (0..n).map(|i| (i, n - 1 - i)).collect();
        mesh.triangles = mesh
            .triangles
            .iter()
            .map(|t| {
                [
                    remap[&(t[0] as usize)] as u32,
                    remap[&(t[1] as usize)] as u32,
                    remap[&(t[2] as usize)] as u32,
                ]
            })
            .collect();
        mesh.positions.reverse();
    }
    let volumes = volumes_for(&raw, &[[0.5, 0.5, 0.5], [1.5, 0.5, 0.5]]);
    assert_unit_cubes(&volumes);
}

#[test]
fn translation_is_volume_invariant() {
    let mut raw = load_committed();
    for mesh in &mut raw {
        for p in &mut mesh.positions {
            p[0] += 10.0;
            p[1] += 7.0;
        }
    }
    // Seeds travel with the solids.
    let volumes = volumes_for(&raw, &[[10.5, 7.5, 0.5], [11.5, 7.5, 0.5]]);
    assert_unit_cubes(&volumes);
}

#[test]
fn quarter_rotation_is_volume_invariant() {
    // Exact integer 90-degree rotation about z: (x,y,z) -> (-y,x,z).
    let mut raw = load_committed();
    for mesh in &mut raw {
        for p in &mut mesh.positions {
            *p = [-p[1], p[0], p[2]];
        }
    }
    // Seeds rotate identically: cold (0.5,0.5)->(-0.5,0.5); hot
    // (1.5,0.5)->(-0.5,1.5).
    let volumes = volumes_for(&raw, &[[-0.5, 0.5, 0.5], [-0.5, 1.5, 0.5]]);
    assert_unit_cubes(&volumes);
}

#[test]
fn unit_rescaling_scales_volumes_by_the_cube() {
    let mut raw = load_committed();
    for mesh in &mut raw {
        for p in &mut mesh.positions {
            p[0] *= 2.0;
            p[1] *= 2.0;
            p[2] *= 2.0;
        }
    }
    let volumes = volumes_for(&raw, &[[1.0, 1.0, 1.0], [3.0, 1.0, 1.0]]);
    assert_eq!(volumes.len(), 2);
    for id in [1u32, 2] {
        assert!(
            (volumes[&id] - 8.0).abs() < 1e-9,
            "region {id} volume {} != 8.0",
            volumes[&id]
        );
    }
}

#[test]
fn single_reversed_triangle_refuses_as_not_closed() {
    let mut raw = load_committed();
    let t = raw[0].triangles[0];
    raw[0].triangles[0] = [t[0], t[2], t[1]];
    with_cx(|cx| {
        let (vertices, remaps) = weld(&[raw[0].positions.clone(), raw[1].positions.clone()]);
        let regions = region_specs(&raw, &remaps, &[[0.5, 0.5, 0.5], [1.5, 0.5, 0.5]]);
        let plc = UnverifiedPlc::new(vertices, regions);
        let outcome =
            fs_mesh::volumetricize(plc, VolumetricPolicy::fixture_default("m"), cx);
        match outcome {
            Ok(_) => panic!("reversed face must not volumetricize"),
            Err(e) => {
                let rendered = format!("{e:?}");
                assert!(
                    rendered.contains("NotClosedManifold")
                        || rendered.contains("UnrecoveredConstraint")
                        || rendered.contains("Audit"),
                    "typed refusal expected, got {rendered}"
                );
            }
        }
    });
}

#[test]
fn deleted_triangle_refuses_as_open_surface() {
    let mut raw = load_committed();
    raw[0].triangles.pop();
    with_cx(|cx| {
        let (vertices, remaps) = weld(&[raw[0].positions.clone(), raw[1].positions.clone()]);
        let regions = region_specs(&raw, &remaps, &[[0.5, 0.5, 0.5], [1.5, 0.5, 0.5]]);
        let plc = UnverifiedPlc::new(vertices, regions);
        let outcome =
            fs_mesh::volumetricize(plc, VolumetricPolicy::fixture_default("m"), cx);
        assert!(
            outcome.is_err(),
            "an open solid must refuse, never fill silently"
        );
    });
}
