//! Uniform 1→8 refinement of the recovered comb (bridge plan B2, the
//! h-ladder control): each rung multiplies the tets by eight, keeps every
//! region's volume, splits every recovered wall four ways in place (each
//! refined source face is still a boundary face of exactly one tet, with its
//! parent facet inherited), mints no flat tets, and keeps the smallest
//! dihedral bounded (shortest-diagonal octahedron split, Liu–Joe 1996). The
//! numbers asserted were MEASURED on 2026-09-02: two-fin comb, 425 base tets
//! (min dihedral 6.86°, max radius-edge 45.0) → 3,400 → 27,200 tets with min
//! dihedral 4.738° at rung 1 and exactly 4.738° again at rung 2 (the interior
//! class appears once, then persists) and radius-edge unchanged (scale
//! invariant: uniform refinement cannot improve it). The assertion keeps the
//! 0.5× bound. The generation-stability test on a single tet lives in
//! `src/uniform.rs`.

use asupersync::types::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_mesh::{
    RecoveryOptions, RegionId, RegionKind, RegionSpec, UnverifiedPlc, VolumetricPolicy,
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

/// Per-region volume (absolute) and the orientation census: every retained
/// tet must carry the SAME orientation sign (whichever convention the
/// producer uses), and refinement must keep it — children take the parent's
/// exact sign.
fn per_region_volume(labeled: &fs_mesh::LabeledTetComplex) -> (BTreeMap<u32, f64>, i64) {
    let pts = labeled.positions();
    let mut volumes = BTreeMap::new();
    let mut sign_sum = 0i64;
    for (tet, region) in labeled.tets().iter().zip(labeled.region_of_tet()) {
        let p = [
            pts[tet[0] as usize],
            pts[tet[1] as usize],
            pts[tet[2] as usize],
            pts[tet[3] as usize],
        ];
        let signed = dot(sub(p[1], p[0]), cross(sub(p[2], p[0]), sub(p[3], p[0]))) / 6.0;
        assert!(signed != 0.0, "no retained tet is flat");
        sign_sum += if signed > 0.0 { 1 } else { -1 };
        *volumes.entry(region.0).or_insert(0.0) += signed.abs();
    }
    let n = i64::try_from(labeled.tets().len()).expect("tet count");
    assert!(
        sign_sum == n || sign_sum == -n,
        "every retained tet carries one orientation sign; sum {sign_sum} of {n}"
    );
    (volumes, sign_sum.signum())
}

fn sorted_face(face: [u32; 3]) -> [u32; 3] {
    let mut f = face;
    f.sort_unstable();
    f
}

#[test]
fn uniform_refinement_of_the_comb_keeps_volume_walls_and_quality_over_two_rungs() {
    let (verts, tris) = comb(2);
    let spec = RegionSpec {
        id: RegionId(1),
        kind: RegionKind::Solid,
        seed: [0.04, 0.03, 0.0025],
        triangles: tris,
    };
    let policy = VolumetricPolicy {
        length_unit: "m".to_string(),
        recovery: RecoveryOptions::default(),
        max_vertices: verts.len(),
        max_tets: 4_000_000,
    };
    let audited = with_cx(|cx| volumetricize(UnverifiedPlc::new(verts, vec![spec]), policy, cx))
        .expect("comb volumetricizes");
    let base = audited.labeled().clone();
    let (base_volume, base_sign) = per_region_volume(&base);
    let base_census = base.quality();
    let base_sources = base.source_faces().len();
    let mut rung = base.clone();
    for k in 1..=2usize {
        rung = rung.refine_uniform();
        let expected_tets = base.tets().len() * 8usize.pow(k as u32);
        assert_eq!(rung.tets().len(), expected_tets, "rung {k} tet count");
        assert_eq!(
            rung.region_of_tet().len(),
            expected_tets,
            "rung {k} label count"
        );
        assert_eq!(
            rung.source_faces().len(),
            base_sources * 4usize.pow(k as u32),
            "rung {k} source-face count"
        );
        // Volume: per region, to rounding of the midpoints.
        let (volume, sign) = per_region_volume(&rung);
        assert_eq!(sign, base_sign, "rung {k} keeps the base orientation sign");
        for (region, &v0) in &base_volume {
            let v = volume[region];
            assert!(
                ((v - v0) / v0).abs() < 1e-12,
                "rung {k} region {region} volume {v} vs base {v0}"
            );
        }
        // Walls: every refined source face is a face of exactly one tet
        // (a boundary face), and its parent facet id is one the base used.
        let mut face_count: BTreeMap<[u32; 3], u32> = BTreeMap::new();
        for tet in rung.tets() {
            for f in [
                [tet[0], tet[1], tet[2]],
                [tet[0], tet[1], tet[3]],
                [tet[0], tet[2], tet[3]],
                [tet[1], tet[2], tet[3]],
            ] {
                *face_count.entry(sorted_face(f)).or_insert(0) += 1;
            }
        }
        let base_parents: std::collections::BTreeSet<u32> =
            base.source_faces().iter().map(|(_, p)| *p).collect();
        for (face, parent) in rung.source_faces() {
            assert_eq!(
                face_count.get(&sorted_face(*face)).copied(),
                Some(1),
                "rung {k} source face {face:?} must be a boundary face of exactly one tet"
            );
            assert!(
                base_parents.contains(parent),
                "rung {k} parent facet {parent} unknown"
            );
        }
        // Quality: no flat tets, dihedrals bounded. MEASURED 2026-09-02 on the
        // two-fin comb: base 6.86°, rung 1 4.738°, rung 2 4.738° (see the header).
        let census = rung.quality();
        assert_eq!(census.flat_tets, 0, "rung {k} minted flat tets");
        assert!(
            census.min_dihedral_deg >= 0.5 * base_census.min_dihedral_deg,
            "rung {k} min dihedral {} vs base {}",
            census.min_dihedral_deg,
            base_census.min_dihedral_deg
        );
        eprintln!(
            "rung {k}: tets {} vertices {} min_dihedral {:.3} max_radius_edge {:.3} (base {:.3} / {:.3})",
            census.tets,
            census.vertices,
            census.min_dihedral_deg,
            census.max_radius_edge,
            base_census.min_dihedral_deg,
            base_census.max_radius_edge
        );
    }
}
