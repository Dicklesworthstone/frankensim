//! Constrained refinement between recovery and carving (bridge plan B2,
//! WORK 1), MEASURED on the two-fin comb 2026-09-02:
//!
//! * interior-only refinement (walls protected by their equatorial spheres,
//!   `split_walls: false`, the default) inserts NOTHING on this thin body:
//!   all 171 tets above the radius-edge target have circumcenters inside some
//!   wall's equatorial sphere, so the mesh, its volume and its walls are
//!   untouched and the evidence says `no-progress` with 171 encroach-skips —
//!   wall splitting is the whole game for thin bodies;
//! * incremental re-recovery seeded from a previous correspondence is a
//!   no-op when nothing changed (zero Steiner points, every facet satisfied),
//!   which is the invariant the split path builds on.
//!
//! The split path itself (`split_walls: true`) is opt-in and not yet claimed:
//! 23 splits left 13 of 60 facets unrecovered within the recovery budget.

use asupersync::types::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::Point3;
use fs_mesh::{
    RecoveryOptions, RefinementOptions, RegionId, RegionKind, RegionSpec, UnverifiedPlc,
    VolumetricPolicy, delaunay, recover_facets, recover_facets_with_points, recover_segments,
    volumetricize,
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
                seed: 0xC0_3C,
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

const BASE_X: f64 = 0.080;
const BASE_Y: f64 = 0.060;
const BASE_Z: f64 = 0.005;
const FIN_W: f64 = 0.006;
const FIN_H: f64 = 0.020;

struct Builder {
    verts: Vec<[f64; 3]>,
    index: BTreeMap<[u64; 3], u32>,
    tris: Vec<[u32; 3]>,
}

impl Builder {
    fn vid(&mut self, p: [f64; 3]) -> u32 {
        let key = [p[0].to_bits(), p[1].to_bits(), p[2].to_bits()];
        if let Some(&id) = self.index.get(&key) {
            return id;
        }
        let id = u32::try_from(self.verts.len()).expect("vertex count");
        self.verts.push(p);
        self.index.insert(key, id);
        id
    }
    fn quad(&mut self, corners: [[f64; 3]; 4], outward: [f64; 3]) {
        let [a, b, c, d] = corners;
        let u = sub(b, a);
        let v = sub(c, a);
        let n = cross(u, v);
        let (a, b, c, d) = if dot(n, outward) > 0.0 {
            (a, b, c, d)
        } else {
            (a, d, c, b)
        };
        let (ia, ib, ic, id) = (self.vid(a), self.vid(b), self.vid(c), self.vid(d));
        self.tris.push([ia, ib, ic]);
        self.tris.push([ia, ic, id]);
    }
}

/// The same grid-triangulated comb prism as `comb_prism.rs`.
fn comb(fins: usize) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    let fin_x: Vec<(f64, f64)> = (0..fins)
        .map(|k| {
            let x0 = 0.008 + 0.018 * k as f64;
            (x0, x0 + FIN_W)
        })
        .collect();
    let mut xs = vec![0.0];
    for &(x0, x1) in &fin_x {
        xs.push(x0);
        xs.push(x1);
    }
    xs.push(BASE_X);
    let top = BASE_Z + FIN_H;
    let mut b = Builder {
        verts: Vec::new(),
        index: BTreeMap::new(),
        tris: Vec::new(),
    };
    let is_fin = |i: usize| i % 2 == 1;
    for i in 0..xs.len() - 1 {
        let (x0, x1) = (xs[i], xs[i + 1]);
        b.quad(
            [
                [x0, 0.0, 0.0],
                [x1, 0.0, 0.0],
                [x1, BASE_Y, 0.0],
                [x0, BASE_Y, 0.0],
            ],
            [0.0, 0.0, -1.0],
        );
        let z_top = if is_fin(i) { top } else { BASE_Z };
        b.quad(
            [
                [x0, 0.0, z_top],
                [x1, 0.0, z_top],
                [x1, BASE_Y, z_top],
                [x0, BASE_Y, z_top],
            ],
            [0.0, 0.0, 1.0],
        );
        for (y, ny) in [(0.0, -1.0), (BASE_Y, 1.0)] {
            b.quad(
                [[x0, y, 0.0], [x1, y, 0.0], [x1, y, BASE_Z], [x0, y, BASE_Z]],
                [0.0, ny, 0.0],
            );
            if is_fin(i) {
                b.quad(
                    [[x0, y, BASE_Z], [x1, y, BASE_Z], [x1, y, top], [x0, y, top]],
                    [0.0, ny, 0.0],
                );
            }
        }
        if is_fin(i) {
            b.quad(
                [
                    [x0, 0.0, BASE_Z],
                    [x0, BASE_Y, BASE_Z],
                    [x0, BASE_Y, top],
                    [x0, 0.0, top],
                ],
                [-1.0, 0.0, 0.0],
            );
            b.quad(
                [
                    [x1, 0.0, BASE_Z],
                    [x1, BASE_Y, BASE_Z],
                    [x1, BASE_Y, top],
                    [x1, 0.0, top],
                ],
                [1.0, 0.0, 0.0],
            );
        }
    }
    b.quad(
        [
            [0.0, 0.0, 0.0],
            [0.0, BASE_Y, 0.0],
            [0.0, BASE_Y, BASE_Z],
            [0.0, 0.0, BASE_Z],
        ],
        [-1.0, 0.0, 0.0],
    );
    b.quad(
        [
            [BASE_X, 0.0, 0.0],
            [BASE_X, BASE_Y, 0.0],
            [BASE_X, BASE_Y, BASE_Z],
            [BASE_X, 0.0, BASE_Z],
        ],
        [1.0, 0.0, 0.0],
    );
    (b.verts, b.tris)
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

/// Independent minimum dihedral (degrees) from inward face normals.
fn min_dihedral_deg(p: &[[f64; 3]; 4]) -> f64 {
    let faces = [[0, 1, 2, 3], [0, 1, 3, 2], [0, 2, 3, 1], [1, 2, 3, 0]];
    let mut normals = Vec::with_capacity(4);
    for f in faces {
        let n = cross(sub(p[f[1]], p[f[0]]), sub(p[f[2]], p[f[0]]));
        let l = norm(n).max(1e-300);
        let n = [n[0] / l, n[1] / l, n[2] / l];
        let n = if dot(n, sub(p[f[3]], p[f[0]])) > 0.0 {
            n
        } else {
            [-n[0], -n[1], -n[2]]
        };
        normals.push(n);
    }
    let mut worst = 180.0f64;
    for i in 0..4 {
        for j in (i + 1)..4 {
            let c = (-dot(normals[i], normals[j])).clamp(-1.0, 1.0);
            worst = worst.min(c.acos().to_degrees());
        }
    }
    worst
}

/// Independent circumradius over shortest edge (Cayley–Menger free form).

fn region_volumes(labeled: &fs_mesh::LabeledTetComplex) -> BTreeMap<u32, f64> {
    let pts = labeled.positions();
    let mut volumes = BTreeMap::new();
    for (tet, region) in labeled.tets().iter().zip(labeled.region_of_tet()) {
        let p = [
            pts[tet[0] as usize],
            pts[tet[1] as usize],
            pts[tet[2] as usize],
            pts[tet[3] as usize],
        ];
        let v = dot(sub(p[1], p[0]), cross(sub(p[2], p[0]), sub(p[3], p[0]))).abs() / 6.0;
        *volumes.entry(region.0).or_insert(0.0) += v;
    }
    volumes
}

fn comb_policy(verts: usize, refinement: Option<RefinementOptions>) -> VolumetricPolicy {
    VolumetricPolicy {
        length_unit: "m".to_string(),
        recovery: RecoveryOptions::default(),
        max_vertices: verts,
        max_tets: 4_000_000,
        refinement,
    }
}

#[test]
fn interior_only_refinement_leaves_the_thin_comb_untouched_and_says_so() {
    let (verts, tris) = comb(2);
    let spec = |tris: Vec<[u32; 3]>| RegionSpec {
        id: RegionId(1),
        kind: RegionKind::Solid,
        seed: [0.04, 0.03, 0.0025],
        triangles: tris,
    };
    let base = with_cx(|cx| {
        volumetricize(
            UnverifiedPlc::new(verts.clone(), vec![spec(tris.clone())]),
            comb_policy(verts.len(), None),
            cx,
        )
    })
    .expect("unrefined comb volumetricizes");
    let refined = with_cx(|cx| {
        volumetricize(
            UnverifiedPlc::new(verts.clone(), vec![spec(tris)]),
            comb_policy(verts.len(), Some(RefinementOptions::default())),
            cx,
        )
    })
    .expect("refinement never breaks a wall, so the refined comb volumetricizes");
    let base_census = base.labeled().quality();
    let census = refined.labeled().quality();
    let evidence = refined.labeled().recovery().refinement;
    eprintln!(
        "base: tets {} min_dihedral {:.3} max_radius_edge {:.3}; refined: tets {} min_dihedral {:.3} max_radius_edge {:.3}; evidence {}",
        base_census.tets,
        base_census.min_dihedral_deg,
        base_census.max_radius_edge,
        census.tets,
        census.min_dihedral_deg,
        census.max_radius_edge,
        evidence.to_json()
    );
    assert_eq!(base.labeled().recovery().refinement.stop, "off");
    // MEASURED: every offender encroaches a wall; nothing is inserted and the
    // evidence discloses exactly that.
    assert_eq!(evidence.rounds, 1, "{}", evidence.to_json());
    assert_eq!(evidence.steiner_inserted, 0, "{}", evidence.to_json());
    assert_eq!(evidence.stop, "no-progress", "{}", evidence.to_json());
    assert!(evidence.offenders_remaining > 100, "{}", evidence.to_json());
    assert_eq!(
        evidence.encroach_skipped,
        evidence.offenders_remaining,
        "every offender was blocked by wall protection: {}",
        evidence.to_json()
    );
    assert_eq!(evidence.walls_split, 0);
    assert_eq!(census.tets, base_census.tets);
    assert_eq!(census.vertices, base_census.vertices);
    assert_eq!(
        region_volumes(refined.labeled()),
        region_volumes(base.labeled())
    );
    assert!((evidence.worst_after - census.max_radius_edge).abs() < 1e-9 * census.max_radius_edge);
}

/// The invariant the wall-split path builds on: re-recovering from the
/// previous correspondence, with nothing inserted in between, changes
/// nothing and inserts nothing.
#[test]
fn incremental_re_recovery_from_the_previous_tiling_is_a_no_op() {
    let (verts, tris) = comb(2);
    let mut segs = std::collections::BTreeSet::new();
    let mut unique = std::collections::BTreeSet::new();
    for t in &tris {
        for (u, v) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            segs.insert(if u < v { [u, v] } else { [v, u] });
        }
        let mut k = *t;
        k.sort_unstable();
        unique.insert(k);
    }
    let segs: Vec<[u32; 2]> = segs.into_iter().collect();
    let loops: Vec<Vec<u32>> = unique.iter().map(|f| vec![f[0], f[1], f[2]]).collect();
    let opts = RecoveryOptions::default();
    with_cx(|cx| {
        let points: Vec<Point3> = verts
            .iter()
            .map(|p| Point3::new(p[0], p[1], p[2]))
            .collect();
        let mut tetra = delaunay(&points, cx).expect("delaunay");
        let (seg_stats, _) = recover_segments(&mut tetra, &segs, opts, cx).expect("segments");
        assert_eq!(seg_stats.unrecovered, 0);
        let (first, table) = recover_facets(&mut tetra, &loops, opts, cx).expect("facets");
        assert_eq!(first.unrecovered, 0);
        assert!(
            first.steiner_inserted > 0,
            "the comb needs facet Steiner points"
        );
        let (again, table2) =
            recover_facets_with_points(&mut tetra, &loops, &[], &table.rows, opts, cx)
                .expect("incremental re-recovery");
        assert_eq!(again.unrecovered, 0, "{}", again.to_json());
        assert_eq!(
            again.steiner_inserted,
            0,
            "seeded from the previous tiling, nothing needs inserting: {}",
            again.to_json()
        );
        assert_eq!(again.recovered, loops.len() as u64);
        let rows: std::collections::BTreeSet<_> = table.rows.iter().collect();
        let rows2: std::collections::BTreeSet<_> = table2.rows.iter().collect();
        assert_eq!(rows, rows2, "the correspondence is reproduced row for row");
    });
}

/// PROBE (ignored): the opt-in wall-split path under `FS_MESH_TRACE_RECOVERY=1`.
#[test]
#[ignore = "diagnostic for the wall-split increment; run explicitly with the trace variable"]
fn probe_wall_splitting_on_the_comb() {
    let (verts, tris) = comb(2);
    let spec = RegionSpec {
        id: RegionId(1),
        kind: RegionKind::Solid,
        seed: [0.04, 0.03, 0.0025],
        triangles: tris,
    };
    let opts = RefinementOptions {
        split_walls: true,
        max_rounds: 1,
        ..RefinementOptions::default()
    };
    let outcome = with_cx(|cx| {
        volumetricize(
            UnverifiedPlc::new(verts.clone(), vec![spec]),
            comb_policy(verts.len(), Some(opts)),
            cx,
        )
    });
    match outcome {
        Ok(audited) => eprintln!(
            "split path: {}",
            audited.labeled().recovery().refinement.to_json()
        ),
        Err(error) => eprintln!("split path refused: {error:?}"),
    }
}
