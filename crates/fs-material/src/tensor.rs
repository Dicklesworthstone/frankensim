//! Voigt/tensor helpers: invariants, deviators, rotations. Convention:
//! `[xx, yy, zz, xy, yz, zx]` with TENSOR shear components (σ_xy stored
//! once, no engineering factor-2 — the factor lives in the contraction
//! helpers where it belongs).

use crate::Voigt;

/// Trace of a Voigt tensor.
#[must_use]
pub fn trace(v: &Voigt) -> f64 {
    v[0] + v[1] + v[2]
}

/// Deviatoric part.
#[must_use]
pub fn deviator(v: &Voigt) -> Voigt {
    let m = trace(v) / 3.0;
    [v[0] - m, v[1] - m, v[2] - m, v[3], v[4], v[5]]
}

/// Double contraction `a : b` (tensor components ⇒ shear terms twice).
#[must_use]
pub fn contract(a: &Voigt, b: &Voigt) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + 2.0 * (a[3] * b[3] + a[4] * b[4] + a[5] * b[5])
}

/// Von Mises equivalent stress `√(3/2 s:s)` of the deviator.
#[must_use]
pub fn von_mises(v: &Voigt) -> f64 {
    let s = deviator(v);
    fs_math::det::sqrt(1.5 * contract(&s, &s))
}

/// Rotate a Voigt tensor by a 3×3 rotation matrix: `Q ε Qᵀ`.
#[must_use]
pub fn rotate(v: &Voigt, q: &[[f64; 3]; 3]) -> Voigt {
    // Expand to full tensor, rotate, re-pack.
    let t = [[v[0], v[3], v[5]], [v[3], v[1], v[4]], [v[5], v[4], v[2]]];
    let mut qt = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let mut s = 0.0;
            for k in 0..3 {
                for l in 0..3 {
                    s += q[i][k] * t[k][l] * q[j][l];
                }
            }
            qt[i][j] = s;
        }
    }
    [qt[0][0], qt[1][1], qt[2][2], qt[0][1], qt[1][2], qt[0][2]]
}

/// Rotation matrix about a unit axis by an angle (Rodrigues; the
/// objectivity-test helper).
#[must_use]
pub fn rotation(axis: [f64; 3], angle: f64) -> [[f64; 3]; 3] {
    let (s, c) = (fs_math::det::sin(angle), fs_math::det::cos(angle));
    let [x, y, z] = axis;
    let one_c = 1.0 - c;
    [
        [
            c + x * x * one_c,
            x * y * one_c - z * s,
            x * z * one_c + y * s,
        ],
        [
            y * x * one_c + z * s,
            c + y * y * one_c,
            y * z * one_c - x * s,
        ],
        [
            z * x * one_c - y * s,
            z * y * one_c + x * s,
            c + z * z * one_c,
        ],
    ]
}

/// 3×3 matrix product (deformation-gradient helpers).
#[must_use]
pub fn matmul3(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut out = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let mut s = 0.0;
            for (k, bk) in b.iter().enumerate() {
                s += a[i][k] * bk[j];
            }
            out[i][j] = s;
        }
    }
    out
}
