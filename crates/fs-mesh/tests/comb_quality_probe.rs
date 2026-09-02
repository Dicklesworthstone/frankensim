//! PROBE (untracked, bridge plan B2 groundwork): measure the tet quality the
//! production `volumetricize` path leaves on the four-fin comb — min dihedral,
//! max radius-edge, sliver count — so the B2 quality floor is pinned from a
//! measurement, not a textbook default. Prints; asserts only that it ran.

use asupersync::types::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_mesh::{RecoveryOptions, RegionId, RegionKind, RegionSpec, UnverifiedPlc, VolumetricPolicy, volumetricize};
use std::collections::BTreeMap;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey { seed: 0xC0_3C, kernel_id: 1, tile: 0, iteration: 0 },
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

struct Builder { verts: Vec<[f64; 3]>, index: BTreeMap<[u64; 3], u32>, tris: Vec<[u32; 3]> }
impl Builder {
    fn vid(&mut self, p: [f64; 3]) -> u32 {
        let key = [p[0].to_bits(), p[1].to_bits(), p[2].to_bits()];
        if let Some(&id) = self.index.get(&key) { return id; }
        let id = u32::try_from(self.verts.len()).expect("vertex count");
        self.verts.push(p); self.index.insert(key, id); id
    }
    fn quad(&mut self, corners: [[f64; 3]; 4], outward: [f64; 3]) {
        let [a, b, c, d] = corners;
        let u = [b[0]-a[0], b[1]-a[1], b[2]-a[2]];
        let v = [c[0]-a[0], c[1]-a[1], c[2]-a[2]];
        let n = [u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0]];
        let dot = n[0]*outward[0] + n[1]*outward[1] + n[2]*outward[2];
        let (a, b, c, d) = if dot > 0.0 { (a, b, c, d) } else { (a, d, c, b) };
        let (ia, ib, ic, id) = (self.vid(a), self.vid(b), self.vid(c), self.vid(d));
        self.tris.push([ia, ib, ic]); self.tris.push([ia, ic, id]);
    }
}

fn comb(fins: usize) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    let fin_x: Vec<(f64, f64)> = (0..fins).map(|k| { let x0 = 0.008 + 0.018 * k as f64; (x0, x0 + FIN_W) }).collect();
    let mut xs = vec![0.0];
    for &(x0, x1) in &fin_x { xs.push(x0); xs.push(x1); }
    xs.push(BASE_X);
    let top = BASE_Z + FIN_H;
    let mut b = Builder { verts: Vec::new(), index: BTreeMap::new(), tris: Vec::new() };
    let is_fin = |i: usize| i % 2 == 1;
    for i in 0..xs.len() - 1 {
        let (x0, x1) = (xs[i], xs[i + 1]);
        b.quad([[x0,0.0,0.0],[x1,0.0,0.0],[x1,BASE_Y,0.0],[x0,BASE_Y,0.0]], [0.0,0.0,-1.0]);
        let z_top = if is_fin(i) { top } else { BASE_Z };
        b.quad([[x0,0.0,z_top],[x1,0.0,z_top],[x1,BASE_Y,z_top],[x0,BASE_Y,z_top]], [0.0,0.0,1.0]);
        for (y, ny) in [(0.0, -1.0), (BASE_Y, 1.0)] {
            b.quad([[x0,y,0.0],[x1,y,0.0],[x1,y,BASE_Z],[x0,y,BASE_Z]], [0.0,ny,0.0]);
            if is_fin(i) { b.quad([[x0,y,BASE_Z],[x1,y,BASE_Z],[x1,y,top],[x0,y,top]], [0.0,ny,0.0]); }
        }
        if is_fin(i) {
            b.quad([[x0,0.0,BASE_Z],[x0,BASE_Y,BASE_Z],[x0,BASE_Y,top],[x0,0.0,top]], [-1.0,0.0,0.0]);
            b.quad([[x1,0.0,BASE_Z],[x1,BASE_Y,BASE_Z],[x1,BASE_Y,top],[x1,0.0,top]], [1.0,0.0,0.0]);
        }
    }
    b.quad([[0.0,0.0,0.0],[0.0,BASE_Y,0.0],[0.0,BASE_Y,BASE_Z],[0.0,0.0,BASE_Z]], [-1.0,0.0,0.0]);
    b.quad([[BASE_X,0.0,0.0],[BASE_X,BASE_Y,0.0],[BASE_X,BASE_Y,BASE_Z],[BASE_X,0.0,BASE_Z]], [1.0,0.0,0.0]);
    (b.verts, b.tris)
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] { [a[0]-b[0], a[1]-b[1], a[2]-b[2]] }
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] { [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]] }
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 { a[0]*b[0]+a[1]*b[1]+a[2]*b[2] }
fn norm(a: [f64; 3]) -> f64 { dot(a, a).sqrt() }

/// Minimum dihedral angle (degrees) of a tet from its face normals.
fn min_dihedral_deg(p: &[[f64; 3]; 4]) -> f64 {
    let faces = [[0, 1, 2, 3], [0, 1, 3, 2], [0, 2, 3, 1], [1, 2, 3, 0]];
    let mut normals = Vec::with_capacity(4);
    for f in faces {
        let n = cross(sub(p[f[1]], p[f[0]]), sub(p[f[2]], p[f[0]]));
        let n = { let l = norm(n); [n[0]/l, n[1]/l, n[2]/l] };
        // orient inward: point toward the opposite vertex
        let to_opp = sub(p[f[3]], p[f[0]]);
        let n = if dot(n, to_opp) > 0.0 { n } else { [-n[0], -n[1], -n[2]] };
        normals.push(n);
    }
    let mut worst = 180.0f64;
    for i in 0..4 { for j in (i+1)..4 {
        let c = (-dot(normals[i], normals[j])).clamp(-1.0, 1.0);
        worst = worst.min(c.acos().to_degrees());
    } }
    worst
}

/// Radius-edge ratio: circumradius / shortest edge.
fn radius_edge(p: &[[f64; 3]; 4]) -> f64 {
    let a = sub(p[1], p[0]); let b = sub(p[2], p[0]); let c = sub(p[3], p[0]);
    let vol6 = dot(a, cross(b, c));
    let num = {
        let t = [
            cross(b, c)[0]*dot(a,a) + cross(c, a)[0]*dot(b,b) + cross(a, b)[0]*dot(c,c),
            cross(b, c)[1]*dot(a,a) + cross(c, a)[1]*dot(b,b) + cross(a, b)[1]*dot(c,c),
            cross(b, c)[2]*dot(a,a) + cross(c, a)[2]*dot(b,b) + cross(a, b)[2]*dot(c,c),
        ];
        norm(t)
    };
    let r = num / (2.0 * vol6.abs());
    let mut emin = f64::INFINITY;
    for i in 0..4 { for j in (i+1)..4 { emin = emin.min(norm(sub(p[i], p[j]))); } }
    r / emin
}

#[test]
fn comb_quality_probe() {
    for fins in [1usize, 4] {
        let (verts, tris) = comb(fins);
        let spec = RegionSpec { id: RegionId(1), kind: RegionKind::Solid, seed: [0.04, 0.03, 0.0025], triangles: tris };
        let policy = VolumetricPolicy { length_unit: "m".to_string(), recovery: RecoveryOptions::default(), max_vertices: verts.len(), max_tets: 4_000_000 };
        let audited = with_cx(|cx| volumetricize(UnverifiedPlc::new(verts, vec![spec]), policy, cx)).expect("comb volumetricizes");
        let l = audited.labeled();
        let pts = l.positions();
        let mut min_dih = 180.0f64; let mut max_re = 0.0f64; let mut slivers5 = 0; let mut slivers1 = 0; let mut bad_re = 0;
        let mut hist = [0usize; 6]; // <1, <5, <10, <20, <40, >=40 degrees
        let mut min_vol = f64::INFINITY; let mut max_vol = 0.0f64;
        for t in l.tets() {
            let p = [pts[t[0] as usize], pts[t[1] as usize], pts[t[2] as usize], pts[t[3] as usize]];
            let d = min_dihedral_deg(&p); let re = radius_edge(&p);
            let vol = dot(sub(p[1],p[0]), cross(sub(p[2],p[0]), sub(p[3],p[0]))).abs() / 6.0;
            min_vol = min_vol.min(vol); max_vol = max_vol.max(vol);
            min_dih = min_dih.min(d); max_re = max_re.max(re);
            if d < 5.0 { slivers5 += 1; } if d < 1.0 { slivers1 += 1; } if re > 2.0 { bad_re += 1; }
            let bin = if d < 1.0 {0} else if d < 5.0 {1} else if d < 10.0 {2} else if d < 20.0 {3} else if d < 40.0 {4} else {5};
            hist[bin] += 1;
        }
        println!(
            "QUALITY fins={fins} tets={} verts={} min_dihedral_deg={min_dih:.3} max_radius_edge={max_re:.2} slivers(<5deg)={slivers5} slivers(<1deg)={slivers1} radius_edge>2={bad_re} dihedral_hist[<1,<5,<10,<20,<40,>=40]={hist:?} vol_min={min_vol:.3e} vol_max={max_vol:.3e}",
            l.tets().len(), pts.len()
        );
    }
}
