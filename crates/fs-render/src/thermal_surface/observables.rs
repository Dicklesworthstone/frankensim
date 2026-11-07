//! Absolute visible measurements of an accepted diffuse thermal source.
//!
//! The source retains its original SI Planck/emissivity law. CIE responses use
//! the calibrated renderer's published observer, constant inside each 5 nm
//! bin. Source integration splits at those bin boundaries AND emissivity knots,
//! so a narrow accepted emission feature cannot fall between sampling bins.
//! This does not recover observer variation inside a bin or establish a
//! certified quadrature/material error. No white normalization, display
//! exposure, thermal-state advance or additional cooling debit is involved.

use core::f64::consts::PI;

use fs_blake3::ContentHash;
use fs_math::eft::two_sum;

use super::{DiffuseThermalSurface, ThermalSurfaceError};
use crate::photometry::{observer, source::LUMINOUS_EFFICACY_LM_W};
use crate::spectral::{LAMBDA_MAX, LAMBDA_MIN, QUAD_BINS};

/// Gross emission over the renderer's 380–780 nm vacuum-wavelength band.
///
/// Exitance is per unit area of ONE emitting face. Both faces of a sheet need
/// their own area in a power calculation. These values are observations of the
/// source, not net heat exchange with its surroundings. A finite-band source
/// does not acquire a bolometric emissivity or cooling law from this report.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalVisibleObservables {
    surface_identity: ContentHash,
    visible_exitance_w_m2: f64,
    xyz: [f64; 3],
    luminous_exitance_lm_m2: f64,
}

impl ThermalVisibleObservables {
    /// Accepted source, temperature and emissivity-law identity.
    #[must_use]
    pub const fn surface_identity(self) -> ContentHash {
        self.surface_identity
    }

    /// Published observer and its retained 5 nm discretization.
    #[must_use]
    pub const fn observer_id(self) -> &'static str {
        observer::OBSERVER_ID
    }

    /// Integrated visible radiant exitance, `pi * integral(L_lambda d_lambda)`,
    /// in W/m² of one emitting face. This omits infrared/ultraviolet radiation.
    #[must_use]
    pub const fn visible_exitance_w_m2(self) -> f64 {
        self.visible_exitance_w_m2
    }

    /// Absolute `683 * integral(L_lambda * [x_bar, y_bar, z_bar] d_lambda)`.
    /// Y is luminance in cd/m²; X and Z retain the same absolute scaling.
    /// No adaptation to a display white or normalization by source power occurs.
    #[must_use]
    pub const fn xyz(self) -> [f64; 3] {
        self.xyz
    }

    /// Photopic luminance in cd/m², exactly the Y component of absolute XYZ.
    #[must_use]
    pub const fn luminance_cd_m2(self) -> f64 {
        self.xyz[1]
    }

    /// Lambertian luminous exitance `pi * luminance`, lm/m² per emitting face.
    #[must_use]
    pub const fn luminous_exitance_lm_m2(self) -> f64 {
        self.luminous_exitance_lm_m2
    }

    /// CIE xy chromaticity, or `None` when all represented XYZ values are zero.
    /// Scaling before summation keeps the ratio defined for large finite XYZ.
    #[must_use]
    pub fn chromaticity_xy(self) -> Option<[f64; 2]> {
        let scale = self.xyz.into_iter().fold(0.0_f64, f64::max);
        if scale == 0.0 {
            return None;
        }
        let [x, y, z] = self.xyz.map(|value| value / scale);
        let total = x + y + z;
        Some([x / total, y / total])
    }
}

impl DiffuseThermalSurface<'_> {
    /// Measure the original source with the calibrated visible observer.
    ///
    /// A spectral emissivity law must cover the complete 380–780 nm observer
    /// band; missing wavelengths refuse even for a nominally dark source.
    /// Four-point Gauss–Legendre integration on each smooth source interval
    /// uses at most 376 pointwise source evaluations, no heap allocation and
    /// no copied spectral table. The observer remains constant within each
    /// original 5 nm bin; this report does not change tracer packet semantics.
    /// Floating-point overflow refuses instead of publishing infinite light.
    pub fn visible_observables(self) -> Result<ThermalVisibleObservables, ThermalSurfaceError> {
        let knots = match self {
            Self::Gray(_) => &[][..],
            Self::Spectral(surface) => {
                let knots = surface.emissivity_knots();
                if knots[0].0 > LAMBDA_MIN || knots[knots.len() - 1].0 < LAMBDA_MAX {
                    return Err(ThermalSurfaceError::Wavelength);
                }
                knots
            }
        };
        let mut sums = [0.0; 4];
        let mut residuals = [0.0; 4];
        let mut knot = knots.partition_point(|&(wavelength, _)| wavelength <= LAMBDA_MIN);
        let bin_width = (LAMBDA_MAX - LAMBDA_MIN) / QUAD_BINS as f64;
        // Nodes and weights on [0,1], ordered by wavelength for fixed replay.
        const QUADRATURE: [(f64, f64); 4] = [
            (0.069_431_844_202_973_71, 0.173_927_422_568_726_92),
            (0.330_009_478_207_571_87, 0.326_072_577_431_273_05),
            (0.669_990_521_792_428_1, 0.326_072_577_431_273_05),
            (0.930_568_155_797_026_2, 0.173_927_422_568_726_92),
        ];
        for bin in 0..QUAD_BINS {
            let mut lower = LAMBDA_MIN + bin as f64 * bin_width;
            let bin_upper = LAMBDA_MIN + (bin + 1) as f64 * bin_width;
            let [x, y, z] = observer::xyz_at_bin(bin);
            while lower < bin_upper {
                while knot < knots.len() && knots[knot].0 <= lower {
                    knot += 1;
                }
                let upper = knots.get(knot).map_or(bin_upper, |k| k.0.min(bin_upper));
                let width = upper - lower;
                for (fraction, weight) in QUADRATURE {
                    let wavelength = lower + fraction * width;
                    let radiance = self.radiance_w_m2_sr_nm(wavelength)?;
                    let contribution = radiance * (weight * width);
                    for (channel, response) in [1.0, x, y, z].into_iter().enumerate() {
                        let (sum, residual) = two_sum(sums[channel], contribution * response);
                        sums[channel] = sum;
                        residuals[channel] += residual;
                    }
                }
                lower = upper;
            }
        }
        let integral: [f64; 4] = core::array::from_fn(|i| sums[i] + residuals[i]);
        let visible_exitance_w_m2 = PI * integral[0];
        let xyz =
            [integral[1], integral[2], integral[3]].map(|value| LUMINOUS_EFFICACY_LM_W * value);
        let luminous_exitance_lm_m2 = PI * xyz[1];
        if [visible_exitance_w_m2, luminous_exitance_lm_m2]
            .into_iter()
            .chain(xyz)
            .any(|value| !value.is_finite() || value < 0.0)
        {
            return Err(ThermalSurfaceError::Unrepresentable);
        }
        Ok(ThermalVisibleObservables {
            surface_identity: self.identity(),
            visible_exitance_w_m2,
            xyz,
            luminous_exitance_lm_m2,
        })
    }
}

#[cfg(test)]
mod tests;
