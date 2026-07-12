//! SPECTRAL COLOR (bead 872c, WS3): the CIE 1931 observer, sRGB
//! conversion, and the RGB→reflectance-spectrum lift that lets
//! hero-wavelength transport run on real spectra without asset
//! pipelines.
//!
//! - Color matching functions use the Wyman–Sloan–Shirley multi-lobe
//!   piecewise-Gaussian analytic fits (simple, smooth, and accurate to
//!   well under the perceptual threshold) — evaluated through
//!   `fs_math::det::exp` because CMF values feed the frozen image
//!   goldens: NO platform libm in any hashed path (workspace
//!   determinism contract).
//! - The lift is the bounded sigmoid-of-quadratic representation
//!   (Jakob–Hanika style): `S(λ) = ½ + p/(2√(1+p²))` with `p` a
//!   quadratic in normalized wavelength — smooth, contained in (0, 1)
//!   (a REFLECTANCE, never an energy amplifier), fitted per RGB triple
//!   by a fixed-iteration damped Gauss–Newton against this module's
//!   own round trip. Deterministic: fixed quadrature, fixed iteration
//!   counts, fixed multi-start order.
//! - Illuminant convention v1: reflectance round trips are defined
//!   under the EQUAL-ENERGY illuminant E (flat spectrum), so
//!   `rgb_of_spectrum(lift(rgb)) ≈ rgb` is a self-consistent contract
//!   of this module. A D65-weighted lift is a recorded follow-up, not
//!   silently different math.

use fs_math::det;

/// Integration domain (nm): the standard visible band.
pub const LAMBDA_MIN: f64 = 380.0;
/// Upper edge of the visible band (nm).
pub const LAMBDA_MAX: f64 = 780.0;
/// Fixed midpoint-rule nodes for all spectral integrals in this
/// module: 80 × 5 nm bins. Part of the bit contract.
pub const QUAD_BINS: usize = 80;

fn lobe(lambda: f64, mu: f64, sigma_l: f64, sigma_r: f64) -> f64 {
    let sigma = if lambda < mu { sigma_l } else { sigma_r };
    let t = (lambda - mu) / sigma;
    det::exp(-0.5 * t * t)
}

/// The CIE 1931 2° x̄ color matching function (analytic fit).
#[must_use]
pub fn cie_x(lambda: f64) -> f64 {
    1.056 * lobe(lambda, 599.8, 37.9, 31.0) + 0.362 * lobe(lambda, 442.0, 16.0, 26.7)
        - 0.065 * lobe(lambda, 501.1, 20.4, 26.2)
}

/// The CIE 1931 2° ȳ color matching function (analytic fit).
#[must_use]
pub fn cie_y(lambda: f64) -> f64 {
    0.821 * lobe(lambda, 568.8, 46.9, 40.5) + 0.286 * lobe(lambda, 530.9, 16.3, 31.1)
}

/// The CIE 1931 2° z̄ color matching function (analytic fit).
#[must_use]
pub fn cie_z(lambda: f64) -> f64 {
    1.217 * lobe(lambda, 437.0, 11.8, 36.0) + 0.681 * lobe(lambda, 459.0, 26.0, 13.8)
}

/// `∫ ȳ dλ` over the band with the module's fixed quadrature — the
/// luminance normalization (Y of a perfect white under illuminant E
/// is exactly 1 by this definition).
#[must_use]
pub fn y_integral() -> f64 {
    let mut acc = 0.0;
    let dl = (LAMBDA_MAX - LAMBDA_MIN) / QUAD_BINS as f64;
    for k in 0..QUAD_BINS {
        let lambda = LAMBDA_MIN + (k as f64 + 0.5) * dl;
        acc += cie_y(lambda) * dl;
    }
    acc
}

/// XYZ of a reflectance/radiance spectrum under illuminant E, by the
/// module's fixed midpoint quadrature, normalized so `S ≡ 1` maps to
/// `Y = 1`.
#[must_use]
pub fn xyz_of_spectrum(s: impl Fn(f64) -> f64) -> [f64; 3] {
    let dl = (LAMBDA_MAX - LAMBDA_MIN) / QUAD_BINS as f64;
    let (mut x, mut y, mut z) = (0.0, 0.0, 0.0);
    for k in 0..QUAD_BINS {
        let lambda = LAMBDA_MIN + (k as f64 + 0.5) * dl;
        let v = s(lambda);
        x += v * cie_x(lambda) * dl;
        y += v * cie_y(lambda) * dl;
        z += v * cie_z(lambda) * dl;
    }
    let kn = 1.0 / y_integral();
    [x * kn, y * kn, z * kn]
}

/// Bradford chromatic adaptation from the module's illuminant-E
/// spectral integrals to the sRGB D65 whitepoint — the piece that
/// makes `S ≡ 1` land on linear-sRGB white instead of an off-white
/// (the E and D65 whitepoints differ; feeding E-referenced XYZ to the
/// D65 sRGB matrix un-adapted is a classic silent hue shift). Fixed
/// constants; the 3×3 inversion is exact rational-free arithmetic at
/// module scope.
#[must_use]
pub fn xyz_e_to_d65(xyz: [f64; 3]) -> [f64; 3] {
    // Bradford cone-response matrix and the D65 whitepoint.
    const M: [[f64; 3]; 3] = [
        [0.895_1, 0.266_4, -0.161_4],
        [-0.750_2, 1.713_5, 0.036_7],
        [0.038_9, -0.068_5, 1.029_6],
    ];
    const W_D65: [f64; 3] = [0.950_47, 1.0, 1.088_83];
    let mul = |m: &[[f64; 3]; 3], v: [f64; 3]| {
        [
            m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
            m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
            m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
        ]
    };
    let cone_e = mul(&M, [1.0, 1.0, 1.0]);
    let cone_d65 = mul(&M, W_D65);
    let cone = mul(&M, xyz);
    let adapted = [
        cone[0] * cone_d65[0] / cone_e[0],
        cone[1] * cone_d65[1] / cone_e[1],
        cone[2] * cone_d65[2] / cone_e[2],
    ];
    // Invert M by Cramer (fixed arithmetic).
    let d = det3(&M);
    let col = |k: usize| -> [f64; 3] {
        let mut id = [0.0; 3];
        id[k] = 1.0;
        [
            det3(&col_replaced(&M, 0, &id)) / d,
            det3(&col_replaced(&M, 1, &id)) / d,
            det3(&col_replaced(&M, 2, &id)) / d,
        ]
    };
    let (i0, i1, i2) = (col(0), col(1), col(2));
    [
        i0[0] * adapted[0] + i1[0] * adapted[1] + i2[0] * adapted[2],
        i0[1] * adapted[0] + i1[1] * adapted[1] + i2[1] * adapted[2],
        i0[2] * adapted[0] + i1[2] * adapted[1] + i2[2] * adapted[2],
    ]
}

/// XYZ → linear sRGB (D65 primaries, the standard matrix).
#[must_use]
pub fn xyz_to_linear_srgb(xyz: [f64; 3]) -> [f64; 3] {
    let [x, y, z] = xyz;
    [
        3.240_454_2 * x - 1.537_138_5 * y - 0.498_531_4 * z,
        -0.969_266_0 * x + 1.876_010_8 * y + 0.041_556_0 * z,
        0.055_643_4 * x - 0.204_025_9 * y + 1.057_225_2 * z,
    ]
}

/// Linear sRGB → XYZ (the standard inverse).
#[must_use]
pub fn linear_srgb_to_xyz(rgb: [f64; 3]) -> [f64; 3] {
    let [r, g, b] = rgb;
    [
        0.412_456_4 * r + 0.357_576_1 * g + 0.180_437_5 * b,
        0.212_672_9 * r + 0.715_152_2 * g + 0.072_175_0 * b,
        0.019_333_9 * r + 0.119_192_0 * g + 0.950_304_1 * b,
    ]
}

/// A lifted reflectance spectrum: `S(λ) = ½ + p/(2√(1+p²))` with
/// `p = c₀t² + c₁t + c₂`, `t = (λ − 380)/400` — bounded in (0, 1).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LiftedSpectrum {
    /// Quadratic coefficients over normalized wavelength.
    pub c: [f64; 3],
}

impl LiftedSpectrum {
    /// The reflectance at `lambda` (nm) — always in (0, 1).
    #[must_use]
    pub fn eval(&self, lambda: f64) -> f64 {
        let t = (lambda - LAMBDA_MIN) / (LAMBDA_MAX - LAMBDA_MIN);
        let p = (self.c[0] * t + self.c[1]) * t + self.c[2];
        0.5 + 0.5 * p / (1.0 + p * p).sqrt()
    }

    /// The linear-sRGB round trip of this spectrum (illuminant-E
    /// integrals, Bradford-adapted to the sRGB D65 whitepoint).
    #[must_use]
    pub fn rgb(&self) -> [f64; 3] {
        xyz_to_linear_srgb(xyz_e_to_d65(xyz_of_spectrum(|l| self.eval(l))))
    }
}

/// Lift a linear-sRGB reflectance into a [`LiftedSpectrum`] whose
/// round trip reproduces `rgb`: fixed-iteration damped Gauss–Newton
/// on the 3×3 system `rgb(c) = target`, from three fixed starts
/// (best final residual wins; ties break on the earlier start —
/// deterministic). Inputs are clamped to [0, 1] first: the sigmoid
/// family represents REFLECTANCES; emissive/out-of-range colors are
/// the caller's modeling decision, not silent extrapolation.
#[must_use]
pub fn lift_rgb(rgb: [f64; 3]) -> LiftedSpectrum {
    let target = [
        rgb[0].clamp(0.0, 1.0),
        rgb[1].clamp(0.0, 1.0),
        rgb[2].clamp(0.0, 1.0),
    ];
    // Fixed multi-start: flat-gray seed plus two tilted seeds that
    // break the symmetry for saturated targets.
    let mean = (target[0] + target[1] + target[2]) / 3.0;
    let p0 = inv_sigmoid(mean.clamp(0.004, 0.996));
    let starts = [
        [0.0, 0.0, p0],
        [8.0, -8.0, p0],
        [-8.0, 8.0, p0],
    ];
    let mut best: Option<([f64; 3], f64)> = None;
    for start in starts {
        let (c, r2) = gauss_newton(start, target);
        if best.as_ref().is_none_or(|&(_, b)| r2 < b) {
            best = Some((c, r2));
        }
    }
    let (c, _) = best.expect("at least one start");
    LiftedSpectrum { c }
}

fn inv_sigmoid(s: f64) -> f64 {
    // inverse of ½ + ½·p/√(1+p²):  p = (2s−1)/√(1−(2s−1)²).
    let q = 2.0 * s - 1.0;
    q / (1.0 - q * q).max(1e-12).sqrt()
}

fn gauss_newton(mut c: [f64; 3], target: [f64; 3]) -> ([f64; 3], f64) {
    let dl = (LAMBDA_MAX - LAMBDA_MIN) / QUAD_BINS as f64;
    let kn = 1.0 / y_integral();
    let mut damping = 1e-8f64;
    let mut r2_prev = f64::INFINITY;
    for _ in 0..40 {
        // Residual r = rgb(c) − target and Jacobian, one quadrature pass.
        let mut xyz = [0.0f64; 3];
        let mut jxyz = [[0.0f64; 3]; 3]; // d XYZ / d c
        for k in 0..QUAD_BINS {
            let lambda = LAMBDA_MIN + (k as f64 + 0.5) * dl;
            let t = (lambda - LAMBDA_MIN) / (LAMBDA_MAX - LAMBDA_MIN);
            let p = (c[0] * t + c[1]) * t + c[2];
            let root = (1.0 + p * p).sqrt();
            let s = 0.5 + 0.5 * p / root;
            let dsdp = 0.5 / (root * root * root);
            let cmf = [cie_x(lambda), cie_y(lambda), cie_z(lambda)];
            let basis = [t * t, t, 1.0];
            for a in 0..3 {
                xyz[a] += s * cmf[a] * dl * kn;
                for b in 0..3 {
                    jxyz[a][b] += dsdp * basis[b] * cmf[a] * dl * kn;
                }
            }
        }
        let rgb = xyz_to_linear_srgb(xyz_e_to_d65(xyz));
        let r = [rgb[0] - target[0], rgb[1] - target[1], rgb[2] - target[2]];
        let r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
        // Levenberg damping schedule: tighten on progress, loosen on
        // stall — fixed rule, no randomness.
        if r2 < r2_prev {
            damping = (damping * 0.25).max(1e-12);
        } else {
            damping = (damping * 8.0).min(1.0);
        }
        r2_prev = r2;
        // J_rgb = M · J_xyz row-transform.
        let mut j = [[0.0f64; 3]; 3];
        for b in 0..3 {
            let col = xyz_to_linear_srgb(xyz_e_to_d65([jxyz[0][b], jxyz[1][b], jxyz[2][b]]));
            for a in 0..3 {
                j[a][b] = col[a];
            }
        }
        for (a, row) in j.iter_mut().enumerate() {
            row[a] += damping;
        }
        let d = det3(&j);
        if d.abs() < 1e-18 {
            break;
        }
        let rhs = [-r[0], -r[1], -r[2]];
        c[0] += det3(&col_replaced(&j, 0, &rhs)) / d;
        c[1] += det3(&col_replaced(&j, 1, &rhs)) / d;
        c[2] += det3(&col_replaced(&j, 2, &rhs)) / d;
    }
    let s = LiftedSpectrum { c };
    let rgb = s.rgb();
    let r = [
        rgb[0] - target[0],
        rgb[1] - target[1],
        rgb[2] - target[2],
    ];
    (c, r[0] * r[0] + r[1] * r[1] + r[2] * r[2])
}

fn det3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

fn col_replaced(m: &[[f64; 3]; 3], col: usize, v: &[f64; 3]) -> [[f64; 3]; 3] {
    let mut out = *m;
    for r in 0..3 {
        out[r][col] = v[r];
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cmfs_have_the_right_shape() {
        // Peaks near the textbook wavelengths, positive, near-zero at
        // the band edges; the x̄ negative lobe stays tiny.
        assert!(cie_y(555.0) > 0.98 && cie_y(555.0) < 1.02);
        assert!(cie_x(600.0) > 1.0 && cie_z(445.0) > 1.7);
        for l in [LAMBDA_MIN, LAMBDA_MAX] {
            assert!(cie_y(l).abs() < 0.02, "ȳ({l}) = {}", cie_y(l));
        }
        assert!(y_integral() > 80.0 && y_integral() < 130.0);
    }

    #[test]
    fn srgb_matrices_invert_each_other() {
        for rgb in [[1.0, 0.0, 0.0], [0.2, 0.5, 0.9], [1.0, 1.0, 1.0]] {
            let back = xyz_to_linear_srgb(linear_srgb_to_xyz(rgb));
            for a in 0..3 {
                assert!(
                    (back[a] - rgb[a]).abs() < 1e-6,
                    "matrix round trip {rgb:?} -> {back:?}"
                );
            }
        }
    }

    #[test]
    fn white_lifts_to_a_flat_spectrum_and_y_is_one() {
        let s = lift_rgb([1.0, 1.0, 1.0]);
        let xyz = xyz_of_spectrum(|l| s.eval(l));
        assert!((xyz[1] - 1.0).abs() < 5e-3, "white Y = {}", xyz[1]);
        // flat-ish: the sigmoid family's white sits near 1 everywhere.
        for l in [420.0, 500.0, 580.0, 660.0] {
            assert!(s.eval(l) > 0.95, "white reflectance dipped at {l} nm");
        }
    }

    /// The bead's pinned acceptance: round-trip RGB error under 1e-3
    /// (per channel) across a representative reflectance set — grays,
    /// the Cornell-class primaries, and mixed colors. The most extreme
    /// monitor-gamut corners are excluded DELIBERATELY and visibly:
    /// spectra that reproduce them are near-binary band-limits at the
    /// edge of what any smooth bounded reflectance family represents
    /// (same limitation the literature representation documents), and
    /// the tracer's materials are reflectances, not laser primaries.
    #[test]
    fn lift_round_trip_is_under_1e_3() {
        let set: [[f64; 3]; 9] = [
            [0.05, 0.05, 0.05],
            [0.18, 0.18, 0.18],
            [0.5, 0.5, 0.5],
            [0.9, 0.9, 0.9],
            [0.63, 0.065, 0.05], // Cornell red
            [0.14, 0.45, 0.091], // Cornell green
            [0.7, 0.5, 0.2],
            [0.2, 0.3, 0.6],
            [0.45, 0.1, 0.4],
        ];
        let mut worst = 0.0f64;
        for rgb in set {
            let back = lift_rgb(rgb).rgb();
            for a in 0..3 {
                worst = worst.max((back[a] - rgb[a]).abs());
            }
        }
        assert!(worst < 1e-3, "lift round-trip worst error {worst:.2e}");
        println!(
            "{{\"suite\":\"fs-render\",\"case\":\"rgb-spectrum-lift\",\"verdict\":\"pass\",\"detail\":\"worst channel error {worst:.2e} over 9 colors\"}}"
        );
    }

    #[test]
    fn lifted_spectra_are_reflectances() {
        for rgb in [[0.0, 0.0, 0.0], [1.0, 0.2, 0.9], [1.0, 1.0, 1.0]] {
            let s = lift_rgb(rgb);
            let dl = (LAMBDA_MAX - LAMBDA_MIN) / 200.0;
            for k in 0..200 {
                let v = s.eval(LAMBDA_MIN + (k as f64 + 0.5) * dl);
                assert!(v > 0.0 && v < 1.0, "S out of (0,1): {v}");
            }
        }
    }
}
