//! Body corpus for facet recovery (bridge plan B4 / C3): bodies that are not
//! the axis-aligned comb must volumetricize under the default policy with
//! the exact analytic volume and a clean independent audit. Each body here is
//! a falsifier of the recovery machinery; a refusal is a finding, not a
//! tolerance question. MEASURED 2026-09-02 after the fixes this corpus forced:
//!
//! * the two-fin comb rotated out of axis alignment (0.61 rad about z, 0.37
//!   rad about x, translated; non-dyadic coordinates, every facet oblique):
//!   247 tets, 99 vertices, 60/60 facets recovered with 43 Steiner points,
//!   min dihedral 2.5°, max radius-edge 34, exact volume, no flat tets. It
//!   found three defects: the kernel's coplanar-ghost apex rounding onto the
//!   plane (a swallowed hull vertex), the flat-tet repair's blanket refusal
//!   of wall-touching flats, and a boundary sliver that only a boundary drop
//!   removes (CONTRACT items 18 and 19);
//! * a plate with a rectangular through-hole (genus one: inner walls whose
//!   outward normal points into the hole, annular top and bottom): 180 tets,
//!   72 vertices, 5 Steiner points, exact volume;
//! * the tracked four-fin heatsink shell and its rotated twin
//!   (`examples/heatsink-fan/heatsink.stl`, `heatsink-rotated.stl`): the
//!   twin found three more defects — facet tilings that did not meet along a
//!   shared segment (a hole the seed flood walked through), a segment
//!   midpoint the mesh kept beside its own edge, and coplanar face sets that
//!   were double covers (CONTRACT items 20 and 21). 379 tets, 154 vertices,
//!   108/108 facets, exact volume, no flat tet after items 22-23;
//! * a vented enclosure (box shell, interior void reaching the outside
//!   through a slot, thin walls): 188 facets, 1137 tets, 368 vertices,
//!   exact volume, smallest dihedral 2.0°.

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
                seed: 0xB0_D1,
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

struct Builder {
    verts: Vec<[f64; 3]>,
    index: BTreeMap<[u64; 3], u32>,
    tris: Vec<[u32; 3]>,
}

impl Builder {
    fn new() -> Self {
        Builder {
            verts: Vec::new(),
            index: BTreeMap::new(),
            tris: Vec::new(),
        }
    }
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
    /// A quad split into two triangles, wound so the normal points along
    /// `outward`.
    fn quad(&mut self, corners: [[f64; 3]; 4], outward: [f64; 3]) {
        let [a, b, c, d] = corners;
        let n = cross(sub(b, a), sub(c, a));
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

const BASE_X: f64 = 0.080;
const BASE_Y: f64 = 0.060;
const BASE_Z: f64 = 0.005;
const FIN_W: f64 = 0.006;
const FIN_H: f64 = 0.020;

/// The grid-triangulated comb prism of `comb_prism.rs`.
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
    let mut b = Builder::new();
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

fn comb_volume(fins: usize) -> f64 {
    BASE_X * BASE_Y * BASE_Z + fins as f64 * FIN_W * BASE_Y * FIN_H
}

/// Rotate by `angle_z` about z, then `angle_x` about x, then translate.
fn rotate(p: [f64; 3], angle_z: f64, angle_x: f64, shift: [f64; 3]) -> [f64; 3] {
    let (sz, cz) = angle_z.sin_cos();
    let (sx, cx) = angle_x.sin_cos();
    let q = [cz * p[0] - sz * p[1], sz * p[0] + cz * p[1], p[2]];
    let r = [q[0], cx * q[1] - sx * q[2], sx * q[1] + cx * q[2]];
    [r[0] + shift[0], r[1] + shift[1], r[2] + shift[2]]
}

/// Plate `X×Y×Z` with a rectangular through-hole `[hx0,hx1]×[hy0,hy1]`.
fn plate_with_hole() -> (Vec<[f64; 3]>, Vec<[u32; 3]>, f64, [f64; 3]) {
    let (x, y, z) = (0.080, 0.060, 0.005);
    let (hx0, hx1, hy0, hy1) = (0.030, 0.050, 0.020, 0.040);
    let mut b = Builder::new();
    // Outer walls, split where the annular grid puts vertices on them.
    for (xx, nx) in [(0.0, -1.0), (x, 1.0)] {
        for (y0, y1) in [(0.0, hy0), (hy0, hy1), (hy1, y)] {
            b.quad(
                [[xx, y0, 0.0], [xx, y1, 0.0], [xx, y1, z], [xx, y0, z]],
                [nx, 0.0, 0.0],
            );
        }
    }
    for (yy, ny) in [(0.0, -1.0), (y, 1.0)] {
        for (x0, x1) in [(0.0, hx0), (hx0, hx1), (hx1, x)] {
            b.quad(
                [[x0, yy, 0.0], [x1, yy, 0.0], [x1, yy, z], [x0, yy, z]],
                [0.0, ny, 0.0],
            );
        }
    }
    // Hole walls: outward (from the solid) points INTO the hole.
    b.quad(
        [
            [hx0, hy0, 0.0],
            [hx0, hy1, 0.0],
            [hx0, hy1, z],
            [hx0, hy0, z],
        ],
        [1.0, 0.0, 0.0],
    );
    b.quad(
        [
            [hx1, hy0, 0.0],
            [hx1, hy1, 0.0],
            [hx1, hy1, z],
            [hx1, hy0, z],
        ],
        [-1.0, 0.0, 0.0],
    );
    b.quad(
        [
            [hx0, hy0, 0.0],
            [hx1, hy0, 0.0],
            [hx1, hy0, z],
            [hx0, hy0, z],
        ],
        [0.0, 1.0, 0.0],
    );
    b.quad(
        [
            [hx0, hy1, 0.0],
            [hx1, hy1, 0.0],
            [hx1, hy1, z],
            [hx0, hy1, z],
        ],
        [0.0, -1.0, 0.0],
    );
    // Annular top and bottom as a 3x3 grid minus the centre cell, so every
    // edge meets its neighbours at shared vertices (full-height side strips
    // would leave T-junctions at the hole corners).
    for (zz, nz) in [(0.0, -1.0), (z, 1.0)] {
        for (x0, x1, y0, y1) in [
            (0.0, hx0, 0.0, hy0),
            (0.0, hx0, hy0, hy1),
            (0.0, hx0, hy1, y),
            (hx1, x, 0.0, hy0),
            (hx1, x, hy0, hy1),
            (hx1, x, hy1, y),
            (hx0, hx1, 0.0, hy0),
            (hx0, hx1, hy1, y),
        ] {
            b.quad(
                [[x0, y0, zz], [x1, y0, zz], [x1, y1, zz], [x0, y1, zz]],
                [0.0, 0.0, nz],
            );
        }
    }
    let volume = x * y * z - (hx1 - hx0) * (hy1 - hy0) * z;
    (b.verts, b.tris, volume, [0.010, 0.010, 0.0025])
}

/// A closed box shell of wall thickness `t` with one full-height slot vent
/// through the `x = bx` wall: the solid is the shell, and its interior void
/// reaches the outside through the vent, so nothing is a sealed cavity and
/// every wall is a thin feature. Returns points, triangles, the exact
/// analytic volume and a seed inside the `x = 0` wall.
fn vented_enclosure() -> (Vec<[f64; 3]>, Vec<[u32; 3]>, f64, [f64; 3]) {
    let (bx, by, bz, t) = (0.060f64, 0.040f64, 0.030f64, 0.005f64);
    let vent = (0.014f64, 0.026f64, t, bz - t);
    // Breakpoints shared by EVERY face of both boxes: a face that carries a
    // vent edge subdivides its neighbours' shared edges too, or the shell has
    // T-junctions and is not a closed manifold.
    let sorted = |mut v: Vec<f64>| {
        v.sort_by(f64::total_cmp);
        v.dedup();
        v
    };
    let xs = sorted(vec![0.0, t, bx - t, bx]);
    let ys = sorted(vec![0.0, t, vent.0, vent.1, by - t, by]);
    let zs = sorted(vec![0.0, t, vent.2, vent.3, bz - t, bz]);
    let clip = |v: &[f64], lo: f64, hi: f64| -> Vec<f64> {
        v.iter()
            .copied()
            .filter(|x| *x >= lo - 1e-12 && *x <= hi + 1e-12)
            .collect()
    };
    let mut b = Builder::new();
    let grid_face = |b: &mut Builder,
                     axis: usize,
                     fixed: f64,
                     outward: [f64; 3],
                     us: &[f64],
                     vs: &[f64],
                     skip_vent: bool| {
        for wu in us.windows(2) {
            for wv in vs.windows(2) {
                if skip_vent
                    && wu[0] >= vent.0 - 1e-12
                    && wu[1] <= vent.1 + 1e-12
                    && wv[0] >= vent.2 - 1e-12
                    && wv[1] <= vent.3 + 1e-12
                {
                    continue;
                }
                let at = |u: f64, v: f64| match axis {
                    0 => [fixed, u, v],
                    1 => [u, fixed, v],
                    _ => [u, v, fixed],
                };
                b.quad(
                    [
                        at(wu[0], wv[0]),
                        at(wu[1], wv[0]),
                        at(wu[1], wv[1]),
                        at(wu[0], wv[1]),
                    ],
                    outward,
                );
            }
        }
    };
    grid_face(&mut b, 0, 0.0, [-1.0, 0.0, 0.0], &ys, &zs, false);
    grid_face(&mut b, 0, bx, [1.0, 0.0, 0.0], &ys, &zs, true);
    grid_face(&mut b, 1, 0.0, [0.0, -1.0, 0.0], &xs, &zs, false);
    grid_face(&mut b, 1, by, [0.0, 1.0, 0.0], &xs, &zs, false);
    grid_face(&mut b, 2, 0.0, [0.0, 0.0, -1.0], &xs, &ys, false);
    grid_face(&mut b, 2, bz, [0.0, 0.0, 1.0], &xs, &ys, false);
    let (xi, yi, zi) = (
        clip(&xs, t, bx - t),
        clip(&ys, t, by - t),
        clip(&zs, t, bz - t),
    );
    // Outward from the SOLID points INTO the void.
    grid_face(&mut b, 0, t, [1.0, 0.0, 0.0], &yi, &zi, false);
    grid_face(&mut b, 0, bx - t, [-1.0, 0.0, 0.0], &yi, &zi, true);
    grid_face(&mut b, 1, t, [0.0, 1.0, 0.0], &xi, &zi, false);
    grid_face(&mut b, 1, by - t, [0.0, -1.0, 0.0], &xi, &zi, false);
    grid_face(&mut b, 2, t, [0.0, 0.0, 1.0], &xi, &yi, false);
    grid_face(&mut b, 2, bz - t, [0.0, 0.0, -1.0], &xi, &yi, false);
    // The vent tube's four walls, from x = bx - t to x = bx.
    let (vy0, vy1, vz0, vz1) = vent;
    b.quad(
        [
            [bx - t, vy0, vz0],
            [bx, vy0, vz0],
            [bx, vy0, vz1],
            [bx - t, vy0, vz1],
        ],
        [0.0, 1.0, 0.0],
    );
    b.quad(
        [
            [bx - t, vy1, vz0],
            [bx, vy1, vz0],
            [bx, vy1, vz1],
            [bx - t, vy1, vz1],
        ],
        [0.0, -1.0, 0.0],
    );
    b.quad(
        [
            [bx - t, vy0, vz0],
            [bx, vy0, vz0],
            [bx, vy1, vz0],
            [bx - t, vy1, vz0],
        ],
        [0.0, 0.0, 1.0],
    );
    b.quad(
        [
            [bx - t, vy0, vz1],
            [bx, vy0, vz1],
            [bx, vy1, vz1],
            [bx - t, vy1, vz1],
        ],
        [0.0, 0.0, -1.0],
    );
    let void = (bx - 2.0 * t) * (by - 2.0 * t) * (bz - 2.0 * t);
    let tube = (vy1 - vy0) * (vz1 - vz0) * t;
    (
        b.verts,
        b.tris,
        bx * by * bz - void - tube,
        [t / 2.0, by / 2.0, bz / 2.0],
    )
}

fn solve_body(
    verts: Vec<[f64; 3]>,
    tris: Vec<[u32; 3]>,
    seed: [f64; 3],
) -> Result<fs_mesh::AuditedLabeledTetComplex, fs_mesh::VolumetricError> {
    let spec = RegionSpec {
        id: RegionId(1),
        kind: RegionKind::Solid,
        seed,
        triangles: tris,
    };
    let policy = VolumetricPolicy {
        length_unit: "m".to_string(),
        recovery: RecoveryOptions::default(),
        max_vertices: verts.len(),
        max_tets: 4_000_000,
        refinement: None,
    };
    with_cx(|cx| volumetricize(UnverifiedPlc::new(verts, vec![spec]), policy, cx))
}

fn region_volume(audited: &fs_mesh::AuditedLabeledTetComplex) -> f64 {
    audited
        .witness()
        .per_region_auditor
        .iter()
        .map(|(_, v)| *v)
        .sum()
}

#[test]
fn the_rotated_two_fin_comb_volumetricizes_with_the_exact_volume() {
    let (verts, tris) = comb(2);
    let angle_z = 0.61;
    let angle_x = 0.37;
    let shift = [0.1, 0.2, 0.05];
    let rotated: Vec<[f64; 3]> = verts
        .iter()
        .map(|p| rotate(*p, angle_z, angle_x, shift))
        .collect();
    let seed = rotate([0.04, 0.03, 0.0025], angle_z, angle_x, shift);
    let audited = solve_body(rotated, tris, seed).expect("the rotated comb volumetricizes");
    let volume = region_volume(&audited);
    let expected = comb_volume(2);
    assert!(
        ((volume - expected) / expected).abs() < 1e-9,
        "rotated comb volume {volume} vs analytic {expected}"
    );
    let census = audited.labeled().quality();
    assert_eq!(census.flat_tets, 0, "{}", census.to_json());
    eprintln!(
        "rotated comb: tets {} vertices {} min_dihedral {:.3} max_radius_edge {:.3}; recovery {}",
        census.tets,
        census.vertices,
        census.min_dihedral_deg,
        census.max_radius_edge,
        audited.labeled().recovery().facets.to_json()
    );
}

#[test]
fn a_plate_with_a_through_hole_volumetricizes_with_the_exact_volume() {
    let (verts, tris, expected, seed) = plate_with_hole();
    let audited = solve_body(verts, tris, seed).expect("the plate with a hole volumetricizes");
    let volume = region_volume(&audited);
    assert!(
        ((volume - expected) / expected).abs() < 1e-9,
        "plate volume {volume} vs analytic {expected}"
    );
    let census = audited.labeled().quality();
    assert_eq!(census.flat_tets, 0, "{}", census.to_json());
    eprintln!(
        "plate with hole: tets {} vertices {} min_dihedral {:.3} max_radius_edge {:.3}; recovery {}",
        census.tets,
        census.vertices,
        census.min_dihedral_deg,
        census.max_radius_edge,
        audited.labeled().recovery().facets.to_json()
    );
}

/// Regression for the coplanar-ghost apex defect (fs-mesh CONTRACT item 18):
/// the Delaunay of the rotated comb's point set keeps every input vertex and
/// passes the FULL exact audit. Before the fix it silently dropped hull vertex
/// 9 and carried 227 local-Delaunay violations, which downstream surfaced as
/// 14 unrecoverable segments.
#[test]
fn the_rotated_comb_point_set_keeps_every_input_vertex_and_audits_clean() {
    use fs_geom::Point3;
    let (verts, _) = comb(2);
    for (label, angle_z, angle_x) in [("axis-aligned", 0.0, 0.0), ("rotated", 0.61, 0.37)] {
        let pts: Vec<Point3> = verts
            .iter()
            .map(|p| {
                let q = rotate(*p, angle_z, angle_x, [0.1, 0.2, 0.05]);
                Point3::new(q[0], q[1], q[2])
            })
            .collect();
        let tetra = with_cx(|cx| fs_mesh::delaunay(&pts, cx)).expect("delaunay");
        let mut incidence = vec![0usize; pts.len()];
        for tet in tetra.tets() {
            for &v in &tet {
                if (v as usize) < pts.len() {
                    incidence[v as usize] += 1;
                }
            }
        }
        let absent: Vec<usize> = (0..pts.len()).filter(|&i| incidence[i] == 0).collect();
        let audit = tetra.audit(true);
        eprintln!(
            "{label}: {} input vertices, {} tets, absent vertices {:?}; exact audit violations {}: {:?}",
            pts.len(),
            tetra.tets().len(),
            absent,
            audit.violations.len(),
            audit.violations.iter().take(3).collect::<Vec<_>>()
        );
        for &i in absent.iter().take(4) {
            eprintln!("  absent {i}: {:?}", verts[i]);
        }
        assert!(
            absent.is_empty(),
            "{label}: every input vertex must be a vertex of the Delaunay tetrahedralization; absent {absent:?}"
        );
    }
}

/// Minimal ASCII-STL reader for the tracked example shells (no parser crate
/// in this crate's dev-dependencies): vertices deduplicated by bits, one
/// triangle per facet, winding as written.
fn read_ascii_stl(path: &std::path::Path) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    let text = std::fs::read_to_string(path).expect("tracked STL is readable");
    let mut b = Builder::new();
    let mut current: Vec<[f64; 3]> = Vec::new();
    for line in text.lines() {
        let mut words = line.split_whitespace();
        if words.next() == Some("vertex") {
            let mut v = [0.0f64; 3];
            for k in 0..3 {
                v[k] = words
                    .next()
                    .and_then(|w| w.parse::<f64>().ok())
                    .expect("vertex coordinate");
            }
            current.push(v);
            if current.len() == 3 {
                let ids = [b.vid(current[0]), b.vid(current[1]), b.vid(current[2])];
                b.tris.push(ids);
                current.clear();
            }
        }
    }
    (b.verts, b.tris)
}

/// A vented enclosure: a box shell whose interior void reaches the outside
/// through a slot, so the flood must keep the void OUTSIDE the solid without
/// any of it being a sealed cavity, and every wall is a thin feature between
/// two constraint surfaces. MEASURED 2026-09-03: 188 facets recovered with a
/// single Steiner point, 1137 tets, 368 vertices, EXACT volume, no flat tet,
/// smallest dihedral 2.0°.
#[test]
fn a_vented_enclosure_volumetricizes_with_the_exact_volume() {
    let (verts, tris, expected, seed) = vented_enclosure();
    assert_eq!(tris.len(), 188, "the tracked coarse shell");
    let audited = solve_body(verts, tris, seed).expect("the vented enclosure volumetricizes");
    let volume = region_volume(&audited);
    assert!(
        ((volume - expected) / expected).abs() < 1e-9,
        "volume {volume} vs analytic {expected}"
    );
    let facets = audited.labeled().recovery().facets;
    assert_eq!(facets.recovered, 188, "{}", facets.to_json());
    let census = audited.labeled().quality();
    assert_eq!(
        (census.tets, census.vertices, census.flat_tets),
        (1137, 368, 0),
        "{}",
        census.to_json()
    );
    assert!(
        census.min_dihedral_deg > 1.0,
        "above the conduction floor: {}",
        census.to_json()
    );
    eprintln!(
        "vented enclosure: tets {} vertices {} min_dihedral {:.3} max_radius_edge {:.3}; recovery {}",
        census.tets,
        census.vertices,
        census.min_dihedral_deg,
        census.max_radius_edge,
        facets.to_json()
    );
}

/// The tracked heatsink shell (generated by `examples/heatsink-fan/
/// generate_heatsink_stl.py`, 108 facets, 52.8 cm³) volumetricizes with the
/// exact volume from its tracked bytes.
#[test]
fn the_tracked_heatsink_shell_volumetricizes_with_the_exact_volume() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/heatsink-fan");
    let expected = 0.080 * 0.060 * 0.005 + 4.0 * 0.006 * 0.060 * 0.020;
    let (verts, tris) = read_ascii_stl(&root.join("heatsink.stl"));
    assert_eq!(tris.len(), 108, "one closed 108-facet shell");
    let audited =
        solve_body(verts, tris, [0.04, 0.03, 0.0025]).expect("the tracked shell volumetricizes");
    let volume = region_volume(&audited);
    assert!(
        ((volume - expected) / expected).abs() < 1e-9,
        "volume {volume} vs analytic {expected}"
    );
    let census = audited.labeled().quality();
    assert_eq!(census.flat_tets, 0, "{}", census.to_json());
    eprintln!(
        "heatsink.stl: tets {} vertices {} min_dihedral {:.3} max_radius_edge {:.3}; recovery {}",
        census.tets,
        census.vertices,
        census.min_dihedral_deg,
        census.max_radius_edge,
        audited.labeled().recovery().facets.to_json()
    );
}

/// The rotated twin of the tracked shell (`generate_heatsink_stl.py OUT
/// --rotate 35 21 --shift 0.1 0.2 0.05`, tracked as `heatsink-rotated.stl`):
/// four fins, every facet oblique, round-trip-exact coordinates, the seed the
/// tracked `heatsink-fan-rotated.fsim` declares. MEASURED 2026-09-03 after the
/// fixes it forced (CONTRACT items 20–23): 379 tets, 154 vertices (98 Steiner
/// points), 108/108 facets, exact volume, no zero-volume tet left — the two
/// coplanar-cluster flats the flip repair could not clear go in two rounds
/// of Steiner perturbation, and the smallest dihedral is 1.64°. The counts
/// are pinned so a change in recovery is a visible event, not a silent drift.
#[test]
fn the_rotated_heatsink_shell_volumetricizes_with_the_exact_volume() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/heatsink-fan");
    let expected = 0.080 * 0.060 * 0.005 + 4.0 * 0.006 * 0.060 * 0.020;
    let (verts, tris) = read_ascii_stl(&root.join("heatsink-rotated.stl"));
    assert_eq!(tris.len(), 108, "one closed 108-facet shell");
    let seed = [0.1155587886810283, 0.2434655989354321, 0.06936274267667861];
    let audited = solve_body(verts, tris, seed).expect("the rotated shell volumetricizes");
    let volume = region_volume(&audited);
    assert!(
        ((volume - expected) / expected).abs() < 1e-9,
        "volume {volume} vs analytic {expected}"
    );
    let facets = audited.labeled().recovery().facets;
    assert_eq!(facets.recovered, 108, "{}", facets.to_json());
    let census = audited.labeled().quality();
    assert_eq!(
        (census.tets, census.vertices, census.flat_tets),
        (379, 154, 0),
        "{}",
        census.to_json()
    );
    assert!(
        census.min_dihedral_deg > 1.0,
        "above the conduction floor: {}",
        census.to_json()
    );
    eprintln!(
        "heatsink-rotated.stl: tets {} vertices {} min_dihedral {:.3} max_radius_edge {:.3}; recovery {}; repair {}",
        census.tets,
        census.vertices,
        census.min_dihedral_deg,
        census.max_radius_edge,
        facets.to_json(),
        audited.labeled().flat_repair().to_json()
    );
}
