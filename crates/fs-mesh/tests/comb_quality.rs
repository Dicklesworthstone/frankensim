//! The quality census of a retained complex agrees with an independent
//! computation (bridge plan B3/A12): `LabeledTetComplex::quality` is what
//! the conduction receipt discloses, so it is checked here against a
//! second implementation of dihedral, radius-edge and volume on the
//! four-fin comb prism. It also pins the flat-tet repair (bridge plan B2a):
//! before it, the comb carried coplanar quadruples at the fin roots (3 of
//! 220 tets on one fin, 9 of 722 on four, dihedral 0.000°) that a P1 solve
//! would turn into equal-temperature constraints; after edge removal there
//! are none and the smallest dihedral is 6.4° / 6.1° (MEASURED 2026-09-02;
//! the assertion keeps 2x headroom). Radius-edge is still disclosed, not
//! enforced — refinement is B2c.

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
fn radius_edge(p: &[[f64; 3]; 4]) -> Option<f64> {
    let a = sub(p[1], p[0]);
    let b = sub(p[2], p[0]);
    let c = sub(p[3], p[0]);
    let vol6 = dot(a, cross(b, c)).abs();
    if vol6 == 0.0 {
        return None;
    }
    let t = [
        cross(b, c)[0] * dot(a, a) + cross(c, a)[0] * dot(b, b) + cross(a, b)[0] * dot(c, c),
        cross(b, c)[1] * dot(a, a) + cross(c, a)[1] * dot(b, b) + cross(a, b)[1] * dot(c, c),
        cross(b, c)[2] * dot(a, a) + cross(c, a)[2] * dot(b, b) + cross(a, b)[2] * dot(c, c),
    ];
    let r = norm(t) / (2.0 * vol6);
    let mut emin = f64::INFINITY;
    for i in 0..4 {
        for j in (i + 1)..4 {
            emin = emin.min(norm(sub(p[i], p[j])));
        }
    }
    Some(r / emin)
}

#[test]
fn quality_census_matches_an_independent_computation_on_the_comb() {
    for fins in [1usize, 4] {
        let (verts, tris) = comb(fins);
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
        let audited =
            with_cx(|cx| volumetricize(UnverifiedPlc::new(verts, vec![spec]), policy, cx))
                .expect("comb volumetricizes");
        let labeled = audited.labeled();
        let census = labeled.quality();
        let pts = labeled.positions();
        let mut volumes = Vec::new();
        let mut min_dih = 180.0f64;
        let mut slivers = 0u32;
        for t in labeled.tets() {
            let p = [
                pts[t[0] as usize],
                pts[t[1] as usize],
                pts[t[2] as usize],
                pts[t[3] as usize],
            ];
            volumes.push(dot(sub(p[1], p[0]), cross(sub(p[2], p[0]), sub(p[3], p[0]))).abs() / 6.0);
            let d = min_dihedral_deg(&p);
            min_dih = min_dih.min(d);
            if d < 5.0 {
                slivers += 1;
            }
        }
        let largest = volumes.iter().copied().fold(0.0, f64::max);
        let mut flats = 0u32;
        let mut max_re = 0.0f64;
        for (t, v) in labeled.tets().iter().zip(&volumes) {
            if *v <= 1e-9 * largest {
                flats += 1;
                continue;
            }
            let p = [
                pts[t[0] as usize],
                pts[t[1] as usize],
                pts[t[2] as usize],
                pts[t[3] as usize],
            ];
            if let Some(re) = radius_edge(&p) {
                max_re = max_re.max(re);
            }
        }
        println!(
            "census fins={fins}: {} | repair: {} | independent: min_dihedral={min_dih:.6} max_radius_edge={max_re:.6} slivers={slivers} flats={flats}",
            census.to_json(),
            labeled.flat_repair().to_json()
        );
        // The repair pass ran before the audit: whatever flat tets remain are
        // exactly the ones it could not remove, and every removal kept the
        // volume (the audit that produced `audited` re-checked it).
        assert_eq!(census.flat_tets, labeled.flat_repair().unrepaired);
        assert_eq!(
            labeled.flat_repair().found,
            labeled.flat_repair().repaired + labeled.flat_repair().unrepaired
        );
        assert!(
            labeled.flat_repair().found > 0,
            "the comb is the flat-tet fixture; if the kernel stops producing them, retire this pin"
        );
        assert_eq!(
            census.flat_tets,
            0,
            "flat tets remain after repair: {}",
            census.to_json()
        );
        assert_eq!(census.slivers_below_5deg, 0, "{}", census.to_json());
        assert!(
            census.min_dihedral_deg >= 3.0,
            "min dihedral regressed below 3 deg (measured 6.4 / 6.1 on 2026-09-02): {}",
            census.to_json()
        );
        assert_eq!(census.tets, labeled.tets().len());
        assert_eq!(census.vertices, pts.len());
        assert_eq!(census.flat_tets, flats, "flat-tet count");
        assert_eq!(census.slivers_below_5deg, slivers, "sliver count");
        assert!(
            (census.min_dihedral_deg - min_dih).abs() <= 1e-6,
            "min dihedral {} vs {min_dih}",
            census.min_dihedral_deg
        );
        assert!(
            (census.max_radius_edge - max_re).abs() <= 1e-6 * max_re.max(1.0),
            "max radius-edge {} vs {max_re}",
            census.max_radius_edge
        );
        // Recovery evidence is retained and names the policy that ran.
        let evidence = labeled.recovery();
        assert_eq!(evidence.options, RecoveryOptions::default());
        assert_eq!(evidence.facets.unrecovered, 0);
        assert!(evidence.to_json().contains("\"max_steiner\":4000"));
    }
}
