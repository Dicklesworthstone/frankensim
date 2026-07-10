//! Film response and display transforms (plan §10.5): exposure, white
//! balance, the Hable filmic curve, and the sRGB display encode —
//! deterministic pixel pipelines built on fs-math's strict `exp`/`pow`
//! (cross-ISA bit-deterministic, like everything else in the chain).

/// Exposure scale: `c · 2^ev` (exact powers of two).
#[must_use]
pub fn exposure(c: f64, ev: i32) -> f64 {
    c * fs_math::det::powi(2.0, ev)
}

/// Per-channel white-balance gains.
#[must_use]
pub fn white_balance(rgb: [f64; 3], gains: [f64; 3]) -> [f64; 3] {
    [rgb[0] * gains[0], rgb[1] * gains[1], rgb[2] * gains[2]]
}

/// Hable (Uncharted 2) filmic tone curve, normalized to a white point of
/// 11.2 — the classic operator constants, applied per channel.
#[must_use]
pub fn hable_filmic(x: f64) -> f64 {
    fn partial(x: f64) -> f64 {
        const A: f64 = 0.15;
        const B: f64 = 0.50;
        const C: f64 = 0.10;
        const D: f64 = 0.20;
        const E: f64 = 0.02;
        const F: f64 = 0.30;
        ((x * (A * x + C * B) + D * E) / (x * (A * x + B) + D * F)) - E / F
    }
    let x = x.max(0.0);
    (partial(x) / partial(11.2)).clamp(0.0, 1.0)
}

/// The sRGB display encode (EOTF⁻¹) on linear [0, 1] values, via
/// fs-math's deterministic `pow`.
#[must_use]
pub fn srgb_encode(linear: f64) -> f64 {
    let l = linear.clamp(0.0, 1.0);
    if l <= 0.003_130_8 {
        12.92 * l
    } else {
        1.055 * fs_math::det::pow(l, 1.0 / 2.4) - 0.055
    }
}

/// Quantize an encoded [0, 1] value to 8 bits (round-half-away, exact
/// endpoints).
#[must_use]
pub fn quantize8(encoded: f64) -> u8 {
    (encoded.clamp(0.0, 1.0) * 255.0 + 0.5) as u8
}

/// Full display pipeline for one linear-light pixel: exposure → white
/// balance → filmic → sRGB → 8-bit.
#[must_use]
pub fn display_transform(rgb: [f64; 3], ev: i32, gains: [f64; 3]) -> [u8; 3] {
    let balanced = white_balance(rgb, gains);
    let mut out = [0u8; 3];
    for (o, &c) in out.iter_mut().zip(&balanced) {
        *o = quantize8(srgb_encode(hable_filmic(exposure(c, ev))));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn curve_shape_and_exact_anchors() {
        assert_eq!(quantize8(srgb_encode(0.0)), 0);
        assert_eq!(quantize8(1.0), 255);
        // sRGB knee: linear segment below the threshold.
        assert!((srgb_encode(0.001) - 12.92 * 0.001).abs() < 1e-15);
        // Monotonicity of the display chain.
        let mut last = -1.0;
        for i in 0..=100 {
            let v = hable_filmic(f64::from(i) * 0.12);
            assert!(v >= last, "filmic curve must be monotone");
            last = v;
        }
        // White point maps to 1.0 exactly by normalization.
        assert!((hable_filmic(11.2) - 1.0).abs() < 1e-12);
        // Exposure is exact powers of two.
        assert!((exposure(0.3, 2) - 1.2).abs() < 1e-15);
        assert_eq!(exposure(1.0, -1024).to_bits(), 1_u64 << 50);
        assert_eq!(exposure(1.0, -1074).to_bits(), 1);
        assert!(exposure(1.0, 1024).is_infinite());
    }
}
