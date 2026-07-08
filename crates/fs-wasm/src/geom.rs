//! geom.rs — Tier-2 geometry demos for the WebGL showcase:
//!
//! * [`marching_cubes`] — a real isosurface polygonizer over a `res³` scalar
//!   field, returning a Three.js-ready triangle soup with gradient normals.
//! * [`sdf_volume`]     — a real F-rep / SDF CSG tree sampled into a `res³`
//!   signed-distance grid for WebGL volume raymarching.
//!
//! The isosurface is extracted with **marching tetrahedra** (each grid cube is
//! split into 6 tetrahedra sharing the main diagonal, then each tetrahedron is
//! polygonized generically from its four corner signs). This is the
//! table-light, trap-proof cousin of classic marching cubes: it needs no
//! 4096-entry hand-transcribed case table (which cannot be validated at a
//! glance), yet extracts the identical level set. Vertex normals come from the
//! analytic field gradient (central differences), so shading is correct
//! regardless of triangle winding.
//!
//! Every input is clamped and the triangle budget is capped: nothing traps.

use fs_math::det;

/* ----------------------------------------------------------------------- */
/*  Scalar fields (signed: negative inside, gradient points outward)        */
/* ----------------------------------------------------------------------- */

/// The chosen scalar field `φ(g)` sampled in the normalized cube `g∈[-1,1]³`.
/// `kind`: 0 = Gyroid TPMS, 1 = metaballs (four fixed centres), 2 = torus.
/// Internally each field maps `g` to its natural domain; the returned geometry
/// therefore always lives in `[-1,1]³` for the viz.
fn field_phi(kind: u32, gx: f64, gy: f64, gz: f64) -> f64 {
    match kind {
        0 => {
            // Gyroid: sin x cos y + sin y cos z + sin z cos x, ~1.5 periods.
            let s = std::f64::consts::PI * 1.5;
            let (x, y, z) = (gx * s, gy * s, gz * s);
            det::sin(x) * det::cos(y) + det::sin(y) * det::cos(z) + det::sin(z) * det::cos(x)
        }
        1 => {
            // Metaballs: potential of four centres; φ = 1 − Σ r²/dist² (inside < 0).
            let p = [gx * 1.4, gy * 1.4, gz * 1.4];
            const C: [[f64; 4]; 4] = [
                [0.5, 0.3, 0.0, 0.30],
                [-0.5, -0.2, 0.2, 0.32],
                [0.0, 0.5, -0.4, 0.26],
                [-0.2, -0.5, -0.3, 0.24],
            ];
            let mut sum = 0.0;
            for c in C.iter() {
                let dx = p[0] - c[0];
                let dy = p[1] - c[1];
                let dz = p[2] - c[2];
                sum += c[3] * c[3] / (dx * dx + dy * dy + dz * dz + 1e-4);
            }
            1.0 - sum
        }
        _ => {
            // Torus in the x–z plane (topologically genus-1, "torus-ish").
            let big_r = 0.55f64;
            let small_r = 0.24f64;
            let q = (gx * gx + gz * gz).sqrt() - big_r;
            q * q + gy * gy - small_r * small_r
        }
    }
}

/// Analytic outward normal `normalize(∇φ)` via central differences.
fn field_normal(kind: u32, g: [f64; 3]) -> [f64; 3] {
    let e = 1.0e-3f64;
    let dx = field_phi(kind, g[0] + e, g[1], g[2]) - field_phi(kind, g[0] - e, g[1], g[2]);
    let dy = field_phi(kind, g[0], g[1] + e, g[2]) - field_phi(kind, g[0], g[1] - e, g[2]);
    let dz = field_phi(kind, g[0], g[1], g[2] + e) - field_phi(kind, g[0], g[1], g[2] - e);
    let len = (dx * dx + dy * dy + dz * dz).sqrt();
    if len > 1e-12 {
        [dx / len, dy / len, dz / len]
    } else {
        [0.0, 0.0, 1.0]
    }
}

/// Linear crossing point on the segment `pa→pb` where `φ` reaches `lev`.
fn interp_edge(pa: [f64; 3], pb: [f64; 3], va: f64, vb: f64, lev: f64) -> [f64; 3] {
    let denom = vb - va;
    let mut t = if denom.abs() > 1e-12 { (lev - va) / denom } else { 0.5 };
    t = t.clamp(0.0, 1.0);
    [
        pa[0] + t * (pb[0] - pa[0]),
        pa[1] + t * (pb[1] - pa[1]),
        pa[2] + t * (pb[2] - pa[2]),
    ]
}

/// Push one triangle (three positions) plus its three gradient normals into
/// `out`, in the layout `v0(3) v1(3) v2(3) n0(3) n1(3) n2(3)`.
fn emit_tri(out: &mut Vec<f64>, kind: u32, a: [f64; 3], b: [f64; 3], c: [f64; 3]) {
    for p in [a, b, c] {
        out.push(p[0]);
        out.push(p[1]);
        out.push(p[2]);
    }
    for p in [a, b, c] {
        let nrm = field_normal(kind, p);
        out.push(nrm[0]);
        out.push(nrm[1]);
        out.push(nrm[2]);
    }
}

/// Polygonize one tetrahedron generically from its four corner signs relative
/// to `lev`. Returns the number of triangles emitted (0, 1, or 2).
fn polygonise_tet(p: &[[f64; 3]; 4], val: &[f64; 4], lev: f64, kind: u32, out: &mut Vec<f64>) -> usize {
    let mut inside = [0usize; 4];
    let mut n_in = 0usize;
    let mut outside = [0usize; 4];
    let mut n_out = 0usize;
    for i in 0..4 {
        if val[i] < lev {
            inside[n_in] = i;
            n_in += 1;
        } else {
            outside[n_out] = i;
            n_out += 1;
        }
    }
    let cross = |a: usize, b: usize| interp_edge(p[a], p[b], val[a], val[b], lev);
    match n_in {
        1 => {
            let a = inside[0];
            emit_tri(
                out,
                kind,
                cross(a, outside[0]),
                cross(a, outside[1]),
                cross(a, outside[2]),
            );
            1
        }
        3 => {
            let b = outside[0];
            emit_tri(
                out,
                kind,
                cross(b, inside[0]),
                cross(b, inside[1]),
                cross(b, inside[2]),
            );
            1
        }
        2 => {
            let (a, b) = (inside[0], inside[1]);
            let (c, d) = (outside[0], outside[1]);
            let e_ac = cross(a, c);
            let e_ad = cross(a, d);
            let e_bd = cross(b, d);
            let e_bc = cross(b, c);
            emit_tri(out, kind, e_ac, e_ad, e_bd);
            emit_tri(out, kind, e_ac, e_bd, e_bc);
            2
        }
        _ => 0, // all in or all out — no surface
    }
}

/// Marching-tetrahedra isosurface of a real `res³` scalar field.
///
/// `kind`: 0 = Gyroid TPMS `sin x cos y + …`, 1 = metaballs (four centres),
/// 2 = a torus. `iso` shifts the extraction level (clamped to `[-1,1]`;
/// default `0` gives the natural surface).
///
/// Output layout (a flat `Float64Array` for Three.js):
/// - `[0]` = `triCount`, the number of triangles emitted.
/// - then `18 * triCount` values: per triangle `v0(3) v1(3) v2(3) n0(3) n1(3)
///   n2(3)` — three vertex positions (in the `[-1,1]³` cube) followed by three
///   per-vertex gradient normals. Total length `1 + 18*triCount`.
///
/// `res` clamped to `[8,48]`; extraction stops once `triCount` reaches ~60000
/// (report the actual count via `[0]` / the return length).
pub fn marching_cubes(res_in: usize, kind_in: u32, iso_in: f64) -> Vec<f64> {
    let res = res_in.clamp(8, 48);
    let iso = iso_in.clamp(-1.0, 1.0);
    let kind = if kind_in > 2 { 0 } else { kind_in };
    let tri_cap = 60000usize;

    let coord = |i: usize| -1.0 + 2.0 * i as f64 / (res as f64 - 1.0);
    let vidx = |i: usize, j: usize, k: usize| (k * res + j) * res + i;
    let mut vals = vec![0.0f64; res * res * res];
    for k in 0..res {
        let z = coord(k);
        for j in 0..res {
            let y = coord(j);
            for i in 0..res {
                vals[vidx(i, j, k)] = field_phi(kind, coord(i), y, z);
            }
        }
    }

    // Corner index = i + 2j + 4k (matching the tet decomposition below).
    const CORNER: [[usize; 3]; 8] = [
        [0, 0, 0],
        [1, 0, 0],
        [0, 1, 0],
        [1, 1, 0],
        [0, 0, 1],
        [1, 0, 1],
        [0, 1, 1],
        [1, 1, 1],
    ];
    // Six tetrahedra sharing the cube's 0–7 main diagonal.
    const TETS: [[usize; 4]; 6] = [
        [0, 7, 1, 3],
        [0, 7, 3, 2],
        [0, 7, 2, 6],
        [0, 7, 6, 4],
        [0, 7, 4, 5],
        [0, 7, 5, 1],
    ];

    let mut out: Vec<f64> = Vec::new();
    out.push(0.0); // triCount placeholder
    let mut tri_count = 0usize;

    let cells = res.saturating_sub(1);
    'outer: for k in 0..cells {
        for j in 0..cells {
            for i in 0..cells {
                let mut cpos = [[0.0f64; 3]; 8];
                let mut cval = [0.0f64; 8];
                for c in 0..8 {
                    let ci = i + CORNER[c][0];
                    let cj = j + CORNER[c][1];
                    let ck = k + CORNER[c][2];
                    cpos[c] = [coord(ci), coord(cj), coord(ck)];
                    cval[c] = vals[vidx(ci, cj, ck)];
                }
                for tet in TETS.iter() {
                    let pp = [cpos[tet[0]], cpos[tet[1]], cpos[tet[2]], cpos[tet[3]]];
                    let vv = [cval[tet[0]], cval[tet[1]], cval[tet[2]], cval[tet[3]]];
                    tri_count += polygonise_tet(&pp, &vv, iso, kind, &mut out);
                    if tri_count >= tri_cap {
                        break 'outer;
                    }
                }
            }
        }
    }
    out[0] = tri_count as f64;
    out
}

/* ----------------------------------------------------------------------- */
/*  sdf_volume — an F-rep / SDF CSG tree sampled into a signed-distance grid */
/* ----------------------------------------------------------------------- */

fn sd_sphere(p: [f64; 3], r: f64) -> f64 {
    (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt() - r
}

fn sd_box(p: [f64; 3], b: [f64; 3]) -> f64 {
    let q = [p[0].abs() - b[0], p[1].abs() - b[1], p[2].abs() - b[2]];
    let qx = q[0].max(0.0);
    let qy = q[1].max(0.0);
    let qz = q[2].max(0.0);
    (qx * qx + qy * qy + qz * qz).sqrt() + q[0].max(q[1]).max(q[2]).min(0.0)
}

fn sd_torus(p: [f64; 3], big: f64, small: f64) -> f64 {
    let qx = (p[0] * p[0] + p[2] * p[2]).sqrt() - big;
    (qx * qx + p[1] * p[1]).sqrt() - small
}

/// Polynomial smooth-min (opSmoothUnion).
fn smin(a: f64, b: f64, k: f64) -> f64 {
    let h = (0.5 + 0.5 * (b - a) / k).clamp(0.0, 1.0);
    b * (1.0 - h) + a * h - k * h * (1.0 - h)
}

/// Smooth-max (used for smooth intersection / subtraction).
fn smax(a: f64, b: f64, k: f64) -> f64 {
    -smin(-a, -b, k)
}

fn sdf_scene(kind: u32, x: f64, y: f64, z: f64, t: f64) -> f64 {
    let ang = 2.0 * std::f64::consts::PI * t;
    match kind {
        0 => {
            // Smooth union of a sphere and a box, orbiting/morphing with t.
            let cx = 0.4 * det::cos(ang);
            let cy = 0.4 * det::sin(ang);
            let s = sd_sphere([x - cx, y - cy, z], 0.45);
            let b = sd_box([x + 0.5 * cx, y, z - 0.2 * det::sin(ang)], [0.35, 0.35, 0.35]);
            smin(s, b, 0.25)
        }
        1 => {
            // Box with an orbiting spherical bite carved out (subtraction).
            let b = sd_box([x, y, z], [0.55, 0.55, 0.55]);
            let sph = sd_sphere([x - 0.5 * det::cos(ang), y, z - 0.5 * det::sin(ang)], 0.4);
            smax(b, -sph, 0.1)
        }
        _ => {
            // Intersection of a torus and a pulsing sphere → morphing shell.
            let tor = sd_torus([x, y, z], 0.55, 0.25);
            let sph = sd_sphere([x, y, z], 0.55 + 0.25 * det::sin(ang));
            smax(tor, sph, 0.1)
        }
    }
}

/// Sample a real F-rep / SDF CSG scene into a `res³` signed-distance grid,
/// row-major with **x fastest** (`index = x + res*(y + res*z)`), for WebGL
/// volume raymarching. The scene is animated/morphed by `t∈[0,1]`.
///
/// `kind`: 0 = smooth-union of a moving sphere+box, 1 = box minus an orbiting
/// sphere, 2 = torus ∩ pulsing sphere. Coordinates span `[-1,1]³`; values are
/// signed distances (negative inside the solid), roughly in `[-2,2]`.
///
/// Output length `res*res*res`. `res` clamped to `[8,64]`, `t` to `[0,1]`.
pub fn sdf_volume(res_in: usize, kind_in: u32, t_in: f64) -> Vec<f64> {
    let res = res_in.clamp(8, 64);
    let t = t_in.clamp(0.0, 1.0);
    let kind = if kind_in > 2 { 0 } else { kind_in };
    let coord = |i: usize| -1.0 + 2.0 * i as f64 / (res as f64 - 1.0);
    let mut out = vec![0.0f64; res * res * res];
    for zc in 0..res {
        let z = coord(zc);
        for yc in 0..res {
            let y = coord(yc);
            for xc in 0..res {
                out[xc + res * (yc + res * zc)] = sdf_scene(kind, coord(xc), y, z, t);
            }
        }
    }
    out
}
