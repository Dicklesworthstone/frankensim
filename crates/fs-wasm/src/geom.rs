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
use fs_viz::Grid3;

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
/// `res` is clamped to `[8,48]` and the complete surface is limited to 60,000
/// triangles. Invalid/non-finite input, allocation refusal, or a surface over
/// that budget returns the bounded failure sentinel `[0]`; no partial surface
/// is serialized.
pub fn marching_cubes(res_in: usize, kind_in: u32, iso_in: f64) -> Vec<f64> {
    let res = res_in.clamp(8, 48);
    let iso = iso_in.clamp(-1.0, 1.0);
    let kind = if kind_in > 2 { 0 } else { kind_in };
    let tri_cap = 60000usize;
    let node_limit = res * res * res;
    let Ok(grid) = Grid3::from_fn([res; 3], [-1.0; 3], [1.0; 3], node_limit, |point| {
        field_phi(kind, point[0], point[1], point[2])
    }) else {
        return vec![0.0];
    };
    let Ok(mesh) = grid.isosurface(iso, tri_cap) else {
        return vec![0.0];
    };
    let Some(output_len) = mesh
        .triangles()
        .len()
        .checked_mul(18)
        .and_then(|payload| payload.checked_add(1))
    else {
        return vec![0.0];
    };
    let mut out = Vec::new();
    if out.try_reserve_exact(output_len).is_err() {
        return vec![0.0];
    }
    out.push(mesh.triangles().len() as f64);
    for triangle in mesh.triangles() {
        let points = (*triangle).map(|index| mesh.vertices()[index as usize]);
        for point in points {
            out.extend_from_slice(&point);
        }
        for point in points {
            out.extend_from_slice(&field_normal(kind, point));
        }
    }
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
            let b = sd_box(
                [x + 0.5 * cx, y, z - 0.2 * det::sin(ang)],
                [0.35, 0.35, 0.35],
            );
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

#[cfg(test)]
mod tests {
    use super::marching_cubes;

    #[test]
    fn shared_isosurface_serializes_the_browser_triangle_layout() {
        let output = marching_cubes(16, 2, 0.0);
        let triangle_count = output[0] as usize;
        assert!(triangle_count > 0);
        assert_eq!(output.len(), 1 + 18 * triangle_count);
        assert!(output.iter().all(|value| value.is_finite()));
        for triangle in output[1..].chunks_exact(18) {
            for normal in triangle[9..].chunks_exact(3) {
                let length = normal[0]
                    .mul_add(
                        normal[0],
                        normal[1].mul_add(normal[1], normal[2] * normal[2]),
                    )
                    .sqrt();
                assert!((length - 1.0).abs() < 1e-12);
            }
        }
        assert_eq!(output, marching_cubes(16, 2, 0.0));
    }

    #[test]
    fn shared_isosurface_maps_nonfinite_input_to_the_bounded_sentinel() {
        assert_eq!(marching_cubes(16, 0, f64::NAN), vec![0.0]);
    }
}
