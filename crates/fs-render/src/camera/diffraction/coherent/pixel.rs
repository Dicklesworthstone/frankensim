//! Positive integration of the original coherent field over a finite pixel.
//!
//! The pupil is contained in the unit disk in its own Fourier coordinates.
//! Its intensity spectrum therefore has support |u_x|, |u_y| <= 2, including
//! all source cross terms and the existing defocus phase. With the original
//! medium-dependent diffraction scale s, a subcell of tangent width h spans
//! at most 2*h/s radians of intensity phase. Subdivision uses that support,
//! not sampled intensities which could mistake dark fringes for a zero field.
//! Four-point Gauss rules in each direction then integrate positive powers.
//! This controls resolution; it is not a certified error bound for either
//! quadrature or the original scalar paraxial pupil model.

use super::{
    CoherentDiffractionError, Dd, PreparedCoherentDiffraction, admit_coordinate, admit_count, poll,
    rounded,
};

// Gauss-Legendre nodes and positive weights on [0, 1]. No endpoints are
// sampled; weights add to one before scaling by the physical pixel measure.
const GAUSS: [(f64, f64); 4] = [
    (0.069_431_844_202_973_71, 0.173_927_422_568_726_93),
    (0.330_009_478_207_571_87, 0.326_072_577_431_273_07),
    (0.669_990_521_792_428_1, 0.326_072_577_431_273_07),
    (0.930_568_155_797_026_2, 0.173_927_422_568_726_93),
];

/// Rectangular detector footprint in the original camera's tangent coordinates.
/// Bounds follow the same pupil axes as `PreparedCoherentDiffraction::observe`.
/// For a perpendicular image plane at distance d, map each physical lateral
/// bound x to x/d: irradiance is the tangent density divided by d², so the
/// integrated power is unchanged. No solid-angle cosine Jacobian applies.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoherentDiffractionPixel {
    minimum_tangent: [f64; 2],
    maximum_tangent: [f64; 2],
    extent_tangent: [f64; 2],
    area_tangent: f64,
}

impl CoherentDiffractionPixel {
    /// Admit a nonempty rectangle wholly inside the original paraxial domain.
    pub fn try_new(
        minimum_tangent: [f64; 2],
        maximum_tangent: [f64; 2],
    ) -> Result<Self, CoherentDiffractionError> {
        for axis in 0..2 {
            if !minimum_tangent[axis].is_finite()
                || !maximum_tangent[axis].is_finite()
                || minimum_tangent[axis] >= maximum_tangent[axis]
            {
                return Err(CoherentDiffractionError::InvalidPixel);
            }
        }
        // A rectangle is convex; its farthest corner determines admission.
        admit_coordinate(std::array::from_fn(|axis| {
            minimum_tangent[axis].abs().max(maximum_tangent[axis].abs())
        }))?;
        let extent_tangent =
            std::array::from_fn(|axis| maximum_tangent[axis] - minimum_tangent[axis]);
        let area_tangent = extent_tangent[0] * extent_tangent[1];
        if !area_tangent.is_finite() || area_tangent <= 0.0 {
            return Err(CoherentDiffractionError::Unrepresentable);
        }
        Ok(Self {
            minimum_tangent,
            maximum_tangent,
            extent_tangent,
            area_tangent,
        })
    }

    /// Lower tangent-coordinate bounds on the original pupil axes.
    #[must_use]
    pub const fn minimum_tangent(self) -> [f64; 2] {
        self.minimum_tangent
    }

    /// Upper tangent-coordinate bounds on the original pupil axes.
    #[must_use]
    pub const fn maximum_tangent(self) -> [f64; 2] {
        self.maximum_tangent
    }

    /// Pixel measure per d(tan(theta_x)) d(tan(theta_y)), not steradians.
    #[must_use]
    pub const fn area_tangent(self) -> f64 {
        self.area_tangent
    }
}

/// Resolution and exact logical work bounds checked before any field query.
/// Each subcell uses 16 positive quadrature samples; every sample evaluates
/// each authored pupil wave once, including zero amplitudes. Factor reductions
/// retain the original preparation's independently bounded factor count.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoherentDiffractionPixelLimits {
    /// Largest intensity phase span across one subcell on either axis, in
    /// `(0, 1]` radians. Smaller values refine the integral. This resolution
    /// setting is not a tolerance or a certified numerical error bound.
    pub maximum_phase_span_rad: f64,
    /// Maximum total quadrature samples for this pixel.
    pub maximum_samples: usize,
    /// Maximum original pupil-transfer evaluations across all samples.
    pub maximum_pupil_evaluations: usize,
}

/// Collected monochromatic power from the original mean and covariance banks.
/// This is an integral in watts, not a centre sample or a pixel-average density.
/// No phase or mean field is assigned to an intensity-integrating detector.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoherentDiffractionPixelPower {
    /// Collected power of the deterministic coherent mean.
    pub mean_power_w: f64,
    /// Collected power of independent covariance columns, excluding the mean.
    pub fluctuation_power_w: f64,
    /// Ensemble-mean collected power, without isolated-source renormalization.
    pub power_w: f64,
    /// Subcells on the two original pupil axes; zero for an empty source bank.
    pub subdivisions: [usize; 2],
    /// Number of positive quadrature samples used.
    pub quadrature_samples: usize,
    /// Number of original pupil-transfer evaluations used.
    pub pupil_evaluations: usize,
}

impl PreparedCoherentDiffraction<'_, '_> {
    /// Integrate coherent and partially coherent power over a finite pixel.
    ///
    /// The original medium wavelength, signed pupil phase and physical
    /// aperture set both the propagated fields and the integration resolution.
    /// Positive Gauss weights multiply powers AFTER each source bank's fields
    /// interfere. Powers from distinct covariance columns add independently.
    /// The complete work count is admitted first. The fixed traversal uses
    /// constant scratch, allocates nothing and polls for every sample and wave.
    /// Cancellation publishes no partial result and leaves this preparation
    /// reusable. Refine `maximum_phase_span_rad` to assess quadrature error;
    /// scalar/paraxial and pupil-model errors are separate.
    pub fn integrate_pixel(
        &self,
        pixel: CoherentDiffractionPixel,
        limits: CoherentDiffractionPixelLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CoherentDiffractionPixelPower, CoherentDiffractionError> {
        poll(&mut cancelled)?;
        let (subdivisions, quadrature_samples, pupil_evaluations) =
            self.pixel_work(pixel, limits)?;
        let mut power = [Dd::ZERO; 2];
        if quadrature_samples != 0 {
            let step: [Dd; 2] = std::array::from_fn(|axis| {
                Dd::from_f64(pixel.extent_tangent[axis]) / Dd::from_f64(subdivisions[axis] as f64)
            });
            let cell_area = Dd::from_f64(pixel.area_tangent)
                / Dd::from_f64(subdivisions[0] as f64)
                / Dd::from_f64(subdivisions[1] as f64);
            for y in 0..subdivisions[1] {
                for x in 0..subdivisions[0] {
                    for &(ny, wy) in &GAUSS {
                        for &(nx, wx) in &GAUSS {
                            poll(&mut cancelled)?;
                            let at = std::array::from_fn(|axis| {
                                let (index, node) = if axis == 0 { (x, nx) } else { (y, ny) };
                                (Dd::from_f64(pixel.minimum_tangent[axis])
                                    + step[axis]
                                        * (Dd::from_f64(index as f64) + Dd::from_f64(node)))
                                .to_f64()
                            });
                            let (_, mean, fluctuations) = self.field_moments(at, &mut cancelled)?;
                            let weight = cell_area * Dd::from_f64(wx) * Dd::from_f64(wy);
                            power[0] = power[0] + weight * mean;
                            power[1] = power[1] + weight * fluctuations;
                        }
                    }
                }
            }
        }
        let result = CoherentDiffractionPixelPower {
            mean_power_w: rounded(power[0])?,
            fluctuation_power_w: rounded(power[1])?,
            power_w: rounded(power[0] + power[1])?,
            subdivisions,
            quadrature_samples,
            pupil_evaluations,
        };
        poll(&mut cancelled)?;
        Ok(result)
    }

    fn pixel_work(
        &self,
        pixel: CoherentDiffractionPixel,
        limits: CoherentDiffractionPixelLimits,
    ) -> Result<([usize; 2], usize, usize), CoherentDiffractionError> {
        if !limits.maximum_phase_span_rad.is_finite()
            || limits.maximum_phase_span_rad <= 0.0
            || limits.maximum_phase_span_rad > 1.0
        {
            return Err(CoherentDiffractionError::InvalidPixelResolution);
        }
        if self.waves == 0 {
            return Ok(([0; 2], 0, 0));
        }
        let over_samples = || CoherentDiffractionError::Limit {
            quantity: "pixel quadrature samples",
            requested: usize::MAX,
            maximum: limits.maximum_samples,
        };
        let maximum_step = self.pupil.scale * (0.5 * limits.maximum_phase_span_rad);
        if maximum_step <= 0.0 {
            return Err(CoherentDiffractionError::Unrepresentable);
        }
        let mut subdivisions = [0; 2];
        for (axis, count) in subdivisions.iter_mut().enumerate() {
            // ceil is an exact elementary rounding operation, independent of
            // the platform's transcendental library.
            let needed = (pixel.extent_tangent[axis] / maximum_step).ceil().max(1.0);
            if !needed.is_finite() || needed >= usize::MAX as f64 {
                return Err(over_samples());
            }
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            {
                *count = needed as usize;
            }
        }
        let samples = subdivisions[0]
            .checked_mul(subdivisions[1])
            .and_then(|cells| cells.checked_mul(16))
            .ok_or_else(over_samples)?;
        admit_count("pixel quadrature samples", samples, limits.maximum_samples)?;
        let evaluations =
            samples
                .checked_mul(self.waves)
                .ok_or(CoherentDiffractionError::Limit {
                    quantity: "pixel pupil evaluations",
                    requested: usize::MAX,
                    maximum: limits.maximum_pupil_evaluations,
                })?;
        admit_count(
            "pixel pupil evaluations",
            evaluations,
            limits.maximum_pupil_evaluations,
        )?;
        Ok((subdivisions, samples, evaluations))
    }
}
