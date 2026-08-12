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

fn permutation_sign(sorted: [u32; 3], tri: [u32; 3]) -> i32 {
    // Even permutation of the sorted triple → +1, odd → −1.
    let even = [
        [sorted[0], sorted[1], sorted[2]],
        [sorted[1], sorted[2], sorted[0]],
        [sorted[2], sorted[0], sorted[1]],
    ];
    if even.contains(&tri) { 1 } else { -1 }
}

/// Merge closed box surfaces and drop opposite coincident pairs
/// (the cancelled shared faces of a Boolean union).
fn union_boxes(parts: &[Vec<[f64; 3]>]) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    let (verts, remaps) = weld(parts);
    let mut net: BTreeMap<[u32; 3], i32> = BTreeMap::new();
    let mut sample: BTreeMap<[u32; 3], [u32; 3]> = BTreeMap::new();
    for remap in &remaps {
        for tri in remap_tris(&box_triangles(0), remap) {
            let mut key = tri;
            key.sort_unstable();
            *net.entry(key).or_insert(0) += permutation_sign(key, tri);
            sample.entry(key).or_insert(tri);
        }
    }
    let tris = net
        .into_iter()
        .filter_map(|(key, n)| {
            if n == 0 {
                None
            } else if n > 0 {
                Some(sample[&key])
            } else {
                let t = sample[&key];
                Some([t[0], t[2], t[1]])
            }
        })
        .collect();
    (verts, tris)
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

#[test]
fn l_shape_carves_the_convex_hull_notch() {
    with_cx(|cx| {
        let (verts, tris) = union_boxes(&[
            box_vertices(0.0, 1.0, 0.0, 1.0, 0.0, 1.0),
            box_vertices(1.0, 2.0, 0.0, 1.0, 0.0, 1.0),
            box_vertices(0.0, 1.0, 1.0, 2.0, 0.0, 1.0),
        ]);
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
        assert!(
            out.witness().excluded_exterior > 0.5,
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
