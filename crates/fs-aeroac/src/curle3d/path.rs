//! Source-to-receiver geometry for cached harmonic Green functions.
//!
//! Rounding x-y, its norm, or k*r before trigonometric reduction discards a
//! short path difference on a long carrier path. Compensating the subsequent
//! source sum cannot recover that lost interference phase. Keep the original
//! binary64 coordinates in DD through distance and phase preparation; only
//! cached field coefficients and public distance metadata round to binary64.

use fs_math::{c64::C64, dd::Dd, det};

#[derive(Clone, Copy)]
pub(super) struct HarmonicPath {
    pub(super) radius_m: f64,
    pub(super) direction: [f64; 3],
    radius: Dd,
}

impl HarmonicPath {
    pub(super) fn between(observer: [f64; 3], source: [f64; 3]) -> Option<Self> {
        if observer.iter().chain(&source).any(|x| !x.is_finite()) {
            return None;
        }
        let displacement = core::array::from_fn::<_, 3, _>(|axis| {
            Dd::from_f64(observer[axis]) - Dd::from_f64(source[axis])
        });
        if displacement.iter().any(|&x| !finite(x)) {
            return None;
        }
        // Scale before squaring, as in hypot: neither large nor small admitted
        // separations should overflow/underflow merely while finding the norm.
        let scale = displacement
            .iter()
            .fold(0.0_f64, |scale, x| scale.max(x.hi.abs()));
        if scale == 0.0 {
            return None;
        }
        let scale = Dd::from_f64(scale);
        let scaled = displacement.map(|x| x / scale);
        let radius = scaled
            .iter()
            .fold(Dd::ZERO, |sum, &x| sum + x * x)
            .sqrt()
            * scale;
        if !finite(radius) || radius.hi <= 0.0 {
            return None;
        }
        Some(Self {
            radius_m: radius.to_f64(),
            direction: displacement.map(|x| (x / radius).to_f64()),
            radius,
        })
    }

    /// exp(i Re(k) r), without rounding the many-wavelength carrier product.
    /// The original material's supplied k remains authoritative; this improves
    /// arithmetic, not the material model or the precision of authored inputs.
    pub(super) fn phase(self, wavenumber: f64) -> Option<C64> {
        let angle = self.radius * Dd::from_f64(wavenumber);
        if !finite(angle) {
            return None;
        }
        let (sin, cos) = det::sin_cos_dd(angle);
        (finite(sin) && finite(cos)).then(|| C64::new(cos.to_f64(), sin.to_f64()))
    }
}

fn finite(value: Dd) -> bool {
    value.hi.is_finite() && value.lo.is_finite()
}
