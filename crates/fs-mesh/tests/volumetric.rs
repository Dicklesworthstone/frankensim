//! Corpus for constrained multi-region PLC volumetricization (s93ej.1).
//!
//! Independent checks: winding-backed audit inside `volumetricize`,
//! analytic box volumes, region exclusivity, and carve of convex-hull
//! fill that is not a declared solid.

use asupersync::types::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_mesh::{
    RegionId, RegionKind, RegionSpec, UnverifiedPlc, VolumetricError, VolumetricPolicy,
    box_triangles, box_vertices, volumetricize,
};
use std::collections::BTreeMap;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x593E_1001,
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

fn policy() -> VolumetricPolicy {
    VolumetricPolicy::fixture_default("m")
}

fn weld(parts: &[Vec<[f64; 3]>]) -> (Vec<[f64; 3]>, Vec<Vec<u32>>) {
    let mut verts = Vec::new();
    let mut index: BTreeMap<[u64; 3], u32> = BTreeMap::new();
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

fn solid_box(
    id: u32,
    seed: [f64; 3],
    x0: f64,
    x1: f64,
    y0: f64,
    y1: f64,
    z0: f64,
    z1: f64,
) -> (Vec<[f64; 3]>, RegionSpec) {
    let verts = box_vertices(x0, x1, y0, y1, z0, z1);
    let spec = RegionSpec {
        id: RegionId(id),
        kind: RegionKind::Solid,
        seed,
        triangles: box_triangles(0),
    };
    (verts, spec)
}

#[test]
fn unit_cube_is_one_solid_of_volume_one() {
    with_cx(|cx| {
        let (verts, spec) = solid_box(1, [0.5, 0.5, 0.5], 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
        let first = volumetricize(
            UnverifiedPlc::new(verts.clone(), vec![spec.clone()]),
            policy(),
            cx,
        )
        .expect("cube");
        let second =
            volumetricize(UnverifiedPlc::new(verts, vec![spec]), policy(), cx).expect("replay");
        assert_eq!(first.labeled().tets(), second.labeled().tets());
        assert_eq!(
            first.labeled().region_of_tet(),
            second.labeled().region_of_tet()
        );
        assert!(first.labeled().tets().len() >= 5);
        assert!(
            first
                .labeled()
                .region_of_tet()
                .iter()
                .all(|r| *r == RegionId(1))
        );
        let vol = first.witness().per_region_surface[0].1;
        assert!((vol - 1.0).abs() < 1e-12, "surface volume {vol}");
        let retained = first.witness().per_region_producer[0].1;
        assert!((retained - 1.0).abs() < 1e-12, "retained volume {retained}");
        assert!(first.witness().excluded_exterior.abs() < 1e-12);
    });
}

fn l_vertices() -> Vec<[f64; 3]> {
    // 8 bottom + 8 top of the [0,2]×[0,1]×[0,1] ∪ [0,1]×[1,2]×[0,1] solid.
    vec![
        [0.0, 0.0, 0.0], // 0
        [1.0, 0.0, 0.0], // 1
        [2.0, 0.0, 0.0], // 2
        [0.0, 1.0, 0.0], // 3
        [1.0, 1.0, 0.0], // 4
        [2.0, 1.0, 0.0], // 5
        [0.0, 2.0, 0.0], // 6
        [1.0, 2.0, 0.0], // 7
        [0.0, 0.0, 1.0], // 8
        [1.0, 0.0, 1.0], // 9
        [2.0, 0.0, 1.0], // 10
        [0.0, 1.0, 1.0], // 11
        [1.0, 1.0, 1.0], // 12
        [2.0, 1.0, 1.0], // 13
        [0.0, 2.0, 1.0], // 14
        [1.0, 2.0, 1.0], // 15
    ]
}

fn quad(a: u32, b: u32, c: u32, d: u32) -> [[u32; 3]; 2] {
    [[a, b, c], [a, c, d]]
}

fn l_triangles() -> Vec<[u32; 3]> {
    let mut tris = Vec::new();
    // Bottom, outward −z.
    tris.extend_from_slice(&quad(0, 3, 4, 1));
    tris.extend_from_slice(&quad(1, 4, 5, 2));
    tris.extend_from_slice(&quad(3, 6, 7, 4));
    // Top, outward +z.
    tris.extend_from_slice(&quad(8, 9, 12, 11));
    tris.extend_from_slice(&quad(9, 10, 13, 12));
    tris.extend_from_slice(&quad(11, 12, 15, 14));
    // Vertical perimeter, outward.
    tris.extend_from_slice(&quad(0, 1, 9, 8)); // y = 0
    tris.extend_from_slice(&quad(1, 2, 10, 9));
    tris.extend_from_slice(&quad(2, 5, 13, 10)); // x = 2
    tris.extend_from_slice(&quad(5, 4, 12, 13)); // y = 1, x ∈ [1, 2]
    tris.extend_from_slice(&quad(4, 7, 15, 12)); // x = 1, y ∈ [1, 2]
    tris.extend_from_slice(&quad(7, 6, 14, 15)); // y = 2
    tris.extend_from_slice(&quad(6, 3, 11, 14)); // x = 0, y ∈ [1, 2]
    tris.extend_from_slice(&quad(3, 0, 8, 11)); // x = 0, y ∈ [0, 1]
    tris
}

#[test]
fn l_shape_carves_the_convex_hull_notch() {
    with_cx(|cx| {
        let verts = l_vertices();
        let tris = l_triangles();
        let spec = RegionSpec {
            id: RegionId(1),
            kind: RegionKind::Solid,
            seed: [0.5, 0.5, 0.5],
            triangles: tris,
        };
        let out = volumetricize(UnverifiedPlc::new(verts, vec![spec]), policy(), cx).expect("L");
        let vol = out.witness().per_region_surface[0].1;
        assert!((vol - 3.0).abs() < 1e-9, "L surface volume {vol}");
        let retained = out.witness().per_region_producer[0].1;
        assert!((retained - 3.0).abs() < 1e-9, "L retained {retained}");
        // Convex hull of these 16 vertices is a pentagonal prism: the
        // missing corner (2,2) is not present, so the carved notch is
        // the right triangle (1,1)-(2,1)-(1,2) extruded through z, volume 1/2.
        assert!(
            (out.witness().excluded_exterior - 0.5).abs() < 1e-9,
            "notch must be carved, exterior {}",
            out.witness().excluded_exterior
        );
        assert!(
            out.labeled()
                .region_of_tet()
                .iter()
                .all(|r| *r == RegionId(1))
        );
    });
}

#[test]
fn adjacent_solids_keep_the_interface_and_two_volumes() {
    with_cx(|cx| {
        let left = box_vertices(0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
        let right = box_vertices(1.0, 2.0, 0.0, 1.0, 0.0, 1.0);
        let (verts, remaps) = weld(&[left, right]);
        let regions = vec![
            RegionSpec {
                id: RegionId(1),
                kind: RegionKind::Solid,
                seed: [0.5, 0.5, 0.5],
                triangles: remap_tris(&box_triangles(0), &remaps[0]),
            },
            RegionSpec {
                id: RegionId(2),
                kind: RegionKind::Solid,
                seed: [1.5, 0.5, 0.5],
                triangles: remap_tris(&box_triangles(0), &remaps[1]),
            },
        ];
        let out =
            volumetricize(UnverifiedPlc::new(verts, regions), policy(), cx).expect("adjacent");
        let mut by_region = BTreeMap::new();
        for &(id, vol) in &out.witness().per_region_producer {
            by_region.insert(id, vol);
        }
        assert!((by_region[&RegionId(1)] - 1.0).abs() < 1e-9);
        assert!((by_region[&RegionId(2)] - 1.0).abs() < 1e-9);
        let n1 = out
            .labeled()
            .region_of_tet()
            .iter()
            .filter(|r| **r == RegionId(1))
            .count();
        let n2 = out
            .labeled()
            .region_of_tet()
            .iter()
            .filter(|r| **r == RegionId(2))
            .count();
        assert!(n1 > 0 && n2 > 0);
        assert!(!out.labeled().source_faces().is_empty());
    });
}

/// Unimodular shear `(x,y,z) → (x+z, y, z)`. Volume is preserved and
/// `x = const` faces stay exactly planar, but the shared interface is
/// no longer axis-aligned. A 3-4-5 rotation of the same brick was
/// measured to break f64 coplanarity and storm the Steiner budget.
fn shear_xz(p: [f64; 3]) -> [f64; 3] {
    [p[0] + p[2], p[1], p[2]]
}

#[test]
fn rotated_adjacent_solids_keep_the_interface_and_two_volumes() {
    with_cx(|cx| {
        let left: Vec<[f64; 3]> = box_vertices(0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
            .into_iter()
            .map(shear_xz)
            .collect();
        let right: Vec<[f64; 3]> = box_vertices(1.0, 2.0, 0.0, 1.0, 0.0, 1.0)
            .into_iter()
            .map(shear_xz)
            .collect();
        let (verts, remaps) = weld(&[left, right]);
        let regions = vec![
            RegionSpec {
                id: RegionId(1),
                kind: RegionKind::Solid,
                seed: shear_xz([0.5, 0.5, 0.5]),
                triangles: remap_tris(&box_triangles(0), &remaps[0]),
            },
            RegionSpec {
                id: RegionId(2),
                kind: RegionKind::Solid,
                seed: shear_xz([1.5, 0.5, 0.5]),
                triangles: remap_tris(&box_triangles(0), &remaps[1]),
            },
        ];
        let out = volumetricize(UnverifiedPlc::new(verts, regions), policy(), cx)
            .expect("sheared adjacent");
        let mut by_region = BTreeMap::new();
        for &(id, vol) in &out.witness().per_region_producer {
            by_region.insert(id, vol);
        }
        assert!(
            (by_region[&RegionId(1)] - 1.0).abs() < 1e-9,
            "sheared left volume {}",
            by_region[&RegionId(1)]
        );
        assert!(
            (by_region[&RegionId(2)] - 1.0).abs() < 1e-9,
            "sheared right volume {}",
            by_region[&RegionId(2)]
        );
        let n1 = out
            .labeled()
            .region_of_tet()
            .iter()
            .filter(|r| **r == RegionId(1))
            .count();
        let n2 = out
            .labeled()
            .region_of_tet()
            .iter()
            .filter(|r| **r == RegionId(2))
            .count();
        assert!(n1 > 0 && n2 > 0, "rotated regions empty: n1={n1} n2={n2}");
        assert!(!out.labeled().source_faces().is_empty());
    });
}

#[test]
fn nested_cavity_is_discarded_from_the_solid() {
    with_cx(|cx| {
        let outer = box_vertices(0.0, 3.0, 0.0, 3.0, 0.0, 3.0);
        let inner = box_vertices(1.0, 2.0, 1.0, 2.0, 1.0, 2.0);
        let (verts, remaps) = weld(&[outer, inner]);
        let inner_tris = remap_tris(&box_triangles(0), &remaps[1]);
        let mut solid_tris = remap_tris(&box_triangles(0), &remaps[0]);
        for tri in &inner_tris {
            solid_tris.push([tri[0], tri[2], tri[1]]);
        }
        let regions = vec![
            RegionSpec {
                id: RegionId(1),
                kind: RegionKind::Solid,
                seed: [0.5, 0.5, 0.5],
                triangles: solid_tris,
            },
            RegionSpec {
                id: RegionId(2),
                kind: RegionKind::Cavity,
                seed: [1.5, 1.5, 1.5],
                triangles: inner_tris,
            },
        ];
        let out = volumetricize(UnverifiedPlc::new(verts, regions), policy(), cx).expect("cavity");
        assert_eq!(out.witness().per_region_producer.len(), 1);
        let retained = out.witness().per_region_producer[0].1;
        // Solid surface is the shell (outer minus inner) so volume 26.
        assert!((retained - 26.0).abs() < 1e-8, "shell volume {retained}");
        assert!(out.witness().excluded_cavity > 0.5);
        assert!(
            out.labeled()
                .region_of_tet()
                .iter()
                .all(|r| *r == RegionId(1))
        );
    });
}

#[test]
fn disjoint_solids_carve_the_gap() {
    with_cx(|cx| {
        let a = box_vertices(0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
        let b = box_vertices(2.0, 3.0, 0.0, 1.0, 0.0, 1.0);
        let (verts, remaps) = weld(&[a, b]);
        let regions = vec![
            RegionSpec {
                id: RegionId(1),
                kind: RegionKind::Solid,
                seed: [0.5, 0.5, 0.5],
                triangles: remap_tris(&box_triangles(0), &remaps[0]),
            },
            RegionSpec {
                id: RegionId(2),
                kind: RegionKind::Solid,
                seed: [2.5, 0.5, 0.5],
                triangles: remap_tris(&box_triangles(0), &remaps[1]),
            },
        ];
        let out =
            volumetricize(UnverifiedPlc::new(verts, regions), policy(), cx).expect("disjoint");
        let total: f64 = out
            .witness()
            .per_region_producer
            .iter()
            .map(|(_, v)| *v)
            .sum();
        assert!((total - 2.0).abs() < 1e-9, "two unit cubes {total}");
        assert!(
            out.witness().excluded_exterior > 0.5,
            "gap must be exterior"
        );
    });
}

#[test]
fn open_surface_is_refused() {
    with_cx(|cx| {
        let verts = box_vertices(0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
        let mut tris = box_triangles(0);
        tris.pop();
        let spec = RegionSpec {
            id: RegionId(1),
            kind: RegionKind::Solid,
            seed: [0.5, 0.5, 0.5],
            triangles: tris,
        };
        let err =
            volumetricize(UnverifiedPlc::new(verts, vec![spec]), policy(), cx).expect_err("open");
        assert!(matches!(err, VolumetricError::NotClosedManifold { .. }));
    });
}

#[test]
fn seed_on_a_face_is_refused() {
    with_cx(|cx| {
        let (verts, mut spec) = solid_box(1, [0.5, 0.5, 0.0], 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
        spec.seed = [0.5, 0.5, 0.0];
        let err = volumetricize(UnverifiedPlc::new(verts, vec![spec]), policy(), cx)
            .expect_err("boundary seed");
        assert!(matches!(
            err,
            VolumetricError::SeedOnBoundary { .. } | VolumetricError::SeedNotLocated { .. }
        ));
    });
}

#[test]
fn missing_cavity_seed_refuses_the_enclosed_chamber() {
    with_cx(|cx| {
        let outer = box_vertices(0.0, 3.0, 0.0, 3.0, 0.0, 3.0);
        let inner = box_vertices(1.0, 2.0, 1.0, 2.0, 1.0, 2.0);
        let (verts, remaps) = weld(&[outer, inner]);
        // Solid surface is the Boolean difference: outer plus reversed inner.
        let mut tris = remap_tris(&box_triangles(0), &remaps[0]);
        for tri in remap_tris(&box_triangles(0), &remaps[1]) {
            tris.push([tri[0], tri[2], tri[1]]);
        }
        let spec = RegionSpec {
            id: RegionId(1),
            kind: RegionKind::Solid,
            seed: [0.5, 0.5, 0.5],
            triangles: tris,
        };
        let err = volumetricize(UnverifiedPlc::new(verts, vec![spec]), policy(), cx)
            .expect_err("unlabeled cavity");
        assert!(matches!(err, VolumetricError::UnlabeledChamber));
    });
}

#[test]
fn duplicate_region_and_nonfinite_are_refused() {
    with_cx(|cx| {
        let (verts, spec) = solid_box(1, [0.5, 0.5, 0.5], 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
        let mut clone = spec.clone();
        clone.seed = [0.25, 0.25, 0.25];
        let err = volumetricize(
            UnverifiedPlc::new(verts.clone(), vec![spec.clone(), clone]),
            policy(),
            cx,
        )
        .expect_err("duplicate");
        assert!(matches!(err, VolumetricError::DuplicateRegion { .. }));

        let mut bad = verts;
        bad[0][0] = f64::NAN;
        let err =
            volumetricize(UnverifiedPlc::new(bad, vec![spec]), policy(), cx).expect_err("nan");
        assert!(matches!(err, VolumetricError::NonFinite { .. }));
    });
}

#[test]
fn permutation_of_region_triangles_is_deterministic() {
    with_cx(|cx| {
        let (verts, mut spec) = solid_box(1, [0.5, 0.5, 0.5], 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
        let a = volumetricize(
            UnverifiedPlc::new(verts.clone(), vec![spec.clone()]),
            policy(),
            cx,
        )
        .expect("a");
        spec.triangles.reverse();
        for tri in &mut spec.triangles {
            *tri = [tri[0], tri[2], tri[1]];
        }
        // Reversing every triangle flips the surface; flip back by
        // reversing the triangle list only, keeping orientation.
        let (verts2, mut spec2) = solid_box(1, [0.5, 0.5, 0.5], 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
        spec2.triangles.rotate_left(3);
        spec2.triangles.swap(0, 5);
        let b = volumetricize(UnverifiedPlc::new(verts2, vec![spec2]), policy(), cx).expect("b");
        assert_eq!(a.labeled().tets(), b.labeled().tets());
        let _ = verts;
        let _ = spec;
    });
}

#[test]
fn axis_permutation_is_a_rigid_transform_of_volume_one() {
    with_cx(|cx| {
        let (verts, mut spec) = solid_box(1, [0.5, 0.5, 0.5], 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
        let verts: Vec<[f64; 3]> = verts.iter().map(|p| [p[1], p[2], p[0]]).collect();
        spec.seed = [0.5, 0.5, 0.5];
        let out = volumetricize(UnverifiedPlc::new(verts, vec![spec]), policy(), cx)
            .expect("permuted axes");
        let vol = out.witness().per_region_surface[0].1;
        assert!((vol - 1.0).abs() < 1e-12, "permuted volume {vol}");
    });
}

#[test]
fn unit_rescaling_scales_volume_by_the_cube() {
    with_cx(|cx| {
        let (verts, mut spec) = solid_box(1, [0.5, 0.5, 0.5], 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
        let verts: Vec<[f64; 3]> = verts
            .iter()
            .map(|p| [2.0 * p[0], 2.0 * p[1], 2.0 * p[2]])
            .collect();
        spec.seed = [1.0, 1.0, 1.0];
        let out =
            volumetricize(UnverifiedPlc::new(verts, vec![spec]), policy(), cx).expect("scaled");
        let vol = out.witness().per_region_producer[0].1;
        assert!((vol - 8.0).abs() < 1e-12, "scaled volume {vol}");
    });
}

#[test]
fn thin_slab_keeps_its_analytic_volume() {
    with_cx(|cx| {
        let (verts, spec) = solid_box(1, [0.5, 0.5, 0.005], 0.0, 1.0, 0.0, 1.0, 0.0, 0.01);
        let out = volumetricize(UnverifiedPlc::new(verts, vec![spec]), policy(), cx).expect("slab");
        let vol = out.witness().per_region_producer[0].1;
        assert!((vol - 0.01).abs() < 1e-12, "slab volume {vol}");
    });
}

#[test]
fn overlapping_solids_refuse_the_unlabeled_intersection() {
    with_cx(|cx| {
        let a = box_vertices(0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
        let b = box_vertices(0.5, 1.5, 0.0, 1.0, 0.0, 1.0);
        let (verts, remaps) = weld(&[a, b]);
        let regions = vec![
            RegionSpec {
                id: RegionId(1),
                kind: RegionKind::Solid,
                seed: [0.25, 0.5, 0.5],
                triangles: remap_tris(&box_triangles(0), &remaps[0]),
            },
            RegionSpec {
                id: RegionId(2),
                kind: RegionKind::Solid,
                seed: [1.25, 0.5, 0.5],
                triangles: remap_tris(&box_triangles(0), &remaps[1]),
            },
        ];
        let err =
            volumetricize(UnverifiedPlc::new(verts, regions), policy(), cx).expect_err("overlap");
        assert!(
            matches!(
                err,
                VolumetricError::UnlabeledChamber
                    | VolumetricError::AmbiguousChamber { .. }
                    | VolumetricError::Audit { .. }
            ),
            "overlap refused as {err:?}"
        );
    });
}

#[test]
fn vertex_budget_refuses_before_meshing() {
    with_cx(|cx| {
        let (verts, spec) = solid_box(1, [0.5, 0.5, 0.5], 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
        let mut pol = policy();
        pol.max_vertices = 4;
        let err = volumetricize(UnverifiedPlc::new(verts, vec![spec]), pol, cx).expect_err("cap");
        assert!(matches!(
            err,
            VolumetricError::Budget {
                what: "vertices",
                ..
            }
        ));
    });
}

#[test]
fn pre_cancelled_context_refuses_without_an_audited_mesh() {
    let gate = CancelGate::new();
    gate.request();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x593E_1001,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let (verts, spec) = solid_box(1, [0.5, 0.5, 0.5], 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
        let err = volumetricize(UnverifiedPlc::new(verts, vec![spec]), policy(), &cx)
            .expect_err("cancelled");
        assert!(
            matches!(err, VolumetricError::Mesh(fs_mesh::MeshError::Cancelled)),
            "cancelled as {err:?}"
        );
    });
}
