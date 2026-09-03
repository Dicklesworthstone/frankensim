//! Finned-heatsink comb prism through the production `volumetricize`
//! path — the shape behind `examples/heatsink-fan/heatsink.stl`.
//!
//! Base 80 x 60 x 5 mm, N fins 6 x 60 x 20 mm at x = 8 + 18k mm, ONE closed
//! manifold shell, grid triangulated (every facet an axis-aligned rectangle
//! split on one diagonal; no T-junctions, no slivers). Every rectangle's
//! four corners are co-circular, so the Delaunay kernel's tie-break picks a
//! diagonal per rectangle and the facet's own triangulation is not the one
//! the mesh holds for roughly a third of them; before facet recovery
//! accepted any coplanar tiling, this body refused with unrecovered facets
//! at every Steiner budget and round cap (MEASURED: one fin, 3 of 36 facets
//! at depth 12 and at depth 24). The tests pin that the default policy now
//! volumetricizes one, two and four fins with the analytic volume under the
//! independent audit, and bound the Steiner spend so a regression back to
//! blind bisection (hundreds of points, then refusal) cannot pass quietly.

use asupersync::types::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::Point3;
use fs_mesh::{
    RecoveryOptions, RegionId, RegionKind, RegionSpec, UnverifiedPlc, VolumetricPolicy, delaunay,
    recover_facets, recover_segments, volumetricize,
};
use std::collections::{BTreeMap, BTreeSet};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0xC0_3B,
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
const SEED: [f64; 3] = [0.04, 0.03, 0.0025];

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
    /// Rectangle a-b-c-d (cyclic), oriented so its normal points along `outward`.
    fn quad(&mut self, corners: [[f64; 3]; 4], outward: [f64; 3]) {
        let [a, b, c, d] = corners;
        let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let v = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        let n = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ];
        let dot = n[0] * outward[0] + n[1] * outward[1] + n[2] * outward[2];
        assert!(dot != 0.0, "degenerate quad");
        let (a, b, c, d) = if dot > 0.0 {
            (a, b, c, d)
        } else {
            (a, d, c, b)
        };
        let (ia, ib, ic, id) = (self.vid(a), self.vid(b), self.vid(c), self.vid(d));
        self.tris.push([ia, ib, ic]);
        self.tris.push([ia, ic, id]);
    }
}

/// Comb prism with `fins` fins: (vertices, triangles, analytic volume).
/// The same construction as `examples/heatsink-fan/generate_heatsink_stl.py`.
fn comb(fins: usize) -> (Vec<[f64; 3]>, Vec<[u32; 3]>, f64) {
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
    let is_fin = |i: usize| i % 2 == 1; // intervals alternate gap, fin, gap, ...
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
    // Closed manifold: every directed edge once, its reverse once.
    let mut edges: BTreeMap<(u32, u32), u32> = BTreeMap::new();
    for t in &b.tris {
        for e in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            *edges.entry(e).or_insert(0) += 1;
        }
    }
    for (&(u, v), &k) in &edges {
        assert_eq!(k, 1, "edge ({u},{v}) used {k} times");
        assert!(edges.contains_key(&(v, u)), "edge ({u},{v}) has no twin");
    }
    let volume = BASE_X * BASE_Y * BASE_Z + fins as f64 * FIN_W * BASE_Y * FIN_H;
    (b.verts, b.tris, volume)
}

fn policy(vertices: usize) -> VolumetricPolicy {
    VolumetricPolicy {
        length_unit: "m".to_string(),
        recovery: RecoveryOptions::default(),
        max_vertices: vertices,
        max_tets: 4_000_000,
        refinement: None,
    }
}

fn rel_close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol * b.abs()
}

/// One, two and four fins volumetricize under the DEFAULT recovery policy,
/// pass the independent winding/volume audit, and reproduce the analytic
/// volume from both the producer's and the auditor's tet-volume formula.
#[test]
fn comb_prism_volumetricizes_with_the_analytic_volume_under_the_default_policy() {
    for fins in [1usize, 2, 4] {
        let (verts, tris, analytic) = comb(fins);
        let spec = RegionSpec {
            id: RegionId(1),
            kind: RegionKind::Solid,
            seed: SEED,
            triangles: tris.clone(),
        };
        let audited = with_cx(|cx| {
            volumetricize(
                UnverifiedPlc::new(verts.clone(), vec![spec]),
                policy(verts.len()),
                cx,
            )
        })
        .unwrap_or_else(|error| panic!("{fins}-fin comb refused: {error}"));
        let witness = audited.witness();
        for (label, rows) in [
            ("producer", &witness.per_region_producer),
            ("auditor", &witness.per_region_auditor),
            ("surface", &witness.per_region_surface),
        ] {
            assert_eq!(rows.len(), 1, "{fins}-fin {label}: one region");
            assert_eq!(rows[0].0, RegionId(1));
            assert!(
                rel_close(rows[0].1, analytic, 1e-9),
                "{fins}-fin {label} volume {} vs analytic {analytic}",
                rows[0].1
            );
        }
        assert!(
            !audited.labeled().tets().is_empty(),
            "{fins}-fin comb retained no tets"
        );
    }
}

/// The staged replica of `AdmittedPlc::recover` pins the mechanism: every
/// segment and every facet of the four-fin comb is recovered under the
/// default caps, and the Steiner spend stays far below the budget (MEASURED
/// 139 facet points, 8 rounds; the pre-tiling code spent >1000 on two fins
/// and still refused).
#[test]
fn comb_prism_facets_are_all_recovered_with_a_bounded_steiner_spend() {
    let (verts, tris, _) = comb(4);
    let mut segs: BTreeSet<[u32; 2]> = BTreeSet::new();
    let mut unique: BTreeSet<[u32; 3]> = BTreeSet::new();
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
        assert_eq!(
            seg_stats.unrecovered,
            0,
            "segments: {}",
            seg_stats.to_json()
        );
        let (facet_stats, table) = recover_facets(&mut tetra, &loops, opts, cx).expect("facets");
        assert_eq!(
            facet_stats.unrecovered,
            0,
            "facets: {}",
            facet_stats.to_json()
        );
        assert_eq!(facet_stats.recovered, loops.len() as u64);
        assert!(
            facet_stats.steiner_inserted <= 600,
            "facet Steiner spend regressed: {}",
            facet_stats.to_json()
        );
        // Every facet has rows, and every recorded row is a live mesh face.
        let mut faces = BTreeSet::new();
        for tet in tetra.tets() {
            for skip in 0..4 {
                let mut f: Vec<u32> = (0..4).filter(|&i| i != skip).map(|i| tet[i]).collect();
                if f.contains(&fs_mesh::GHOST) {
                    continue;
                }
                f.sort_unstable();
                faces.insert([f[0], f[1], f[2]]);
            }
        }
        let mut per_facet = vec![0usize; loops.len()];
        for (face, fid) in &table.rows {
            assert!(
                faces.contains(face),
                "recorded row {face:?} is not a mesh face"
            );
            per_facet[*fid as usize] += 1;
        }
        assert!(per_facet.iter().all(|&n| n > 0), "a facet has no rows");
        // Tile areas reproduce every facet's area (no gap, no overlap).
        let pts: Vec<[f64; 3]> = tetra.points().iter().map(|p| [p.x, p.y, p.z]).collect();
        let area = |a: [f64; 3], b: [f64; 3], c: [f64; 3]| {
            let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
            let v = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
            let n = [
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0],
            ];
            0.5 * (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt()
        };
        for (fid, f) in unique.iter().enumerate() {
            let want = area(pts[f[0] as usize], pts[f[1] as usize], pts[f[2] as usize]);
            let got: f64 = table
                .rows
                .iter()
                .filter(|(_, id)| *id as usize == fid)
                .map(|(face, _)| {
                    area(
                        pts[face[0] as usize],
                        pts[face[1] as usize],
                        pts[face[2] as usize],
                    )
                })
                .sum();
            assert!(
                rel_close(got, want, 1e-9),
                "facet {fid} tiles sum to {got} but the facet area is {want}"
            );
        }
    });
}
