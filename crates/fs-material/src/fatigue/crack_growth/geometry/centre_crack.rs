//! Elastic secant geometry for a centred through crack under remote tension.
//!
//! AFGROW DTD Handbook §11.3, Table 11.3.8 (p. 11.3.11), gives
//! `K = sigma_gross * sqrt(pi*a) * sqrt(sec(pi*a/W))`:
//! <https://www.afgrow.net/applications/DTDHandbook/pdfs/Sec11_3.pdf>.
//! Here `W` is the full plate width and `a` the physical half-crack length.
//! NASA TN D-5140, symbols pp. 2–3 and p. 6, supplies these conventions and
//! reports agreement of the secant approximation with Isida within 0.3% for
//! `2*a/W <= 0.8`: <https://ntrs.nasa.gov/citations/19690012030>.
//!
//! That ratio is the chosen supported formula range, not a material fit or an
//! error bound for interpolation of a sampled table. Applicability requires a
//! linear-elastic plate long enough in the loading direction for remote uniform
//! gross-section tension to apply. No numerical grip-distance/aspect criterion,
//! finite-height correction, bending, plastic-zone correction, or 3D crack-front
//! solution is inferred. The NASA report's separate effective-length plasticity
//! correction is deliberately not part of this elastic geometry.

use fs_math::det;

use super::super::{CrackGrowthError, log1p_over_x, positive_log_ratio};
use super::{CrackGeometryPoint, CrackGeometryTableContext, TabulatedCrackGeometry, poll};

/// Finite-width elastic centre-crack solution, with physical half-length `a`.
/// The input stress is the remote gross-section stress, not net-ligament stress.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FiniteWidthCentreCrack {
    plate_width_m: f64,
}

impl FiniteWidthCentreCrack {
    /// Declare the full plate width. The caller establishes the remote-tension
    /// plate assumptions documented by this module; width alone cannot do so.
    pub fn new(plate_width_m: f64) -> Result<Self, CrackGrowthError> {
        if !plate_width_m.is_finite() || plate_width_m <= 0.0 {
            return Err(CrackGrowthError::Invalid(
                "centre-crack full plate width must be finite and positive",
            ));
        }
        let maximum = 0.4 * plate_width_m;
        if maximum == 0.0 || maximum >= 0.5 * plate_width_m {
            return Err(CrackGrowthError::Arithmetic);
        }
        Ok(Self { plate_width_m })
    }

    pub fn plate_width_m(&self) -> f64 {
        self.plate_width_m
    }

    /// Largest supported physical half-length: full crack / full width = 0.8.
    pub fn maximum_half_length_m(&self) -> f64 {
        0.4 * self.plate_width_m
    }

    /// Evaluate `Y = 1/sqrt(cos(pi*a/W))` without table interpolation.
    /// A strictly positive physical half-length is required.
    pub fn factor_at(&self, half_length_m: f64) -> Result<f64, CrackGrowthError> {
        let maximum = self.maximum_half_length_m();
        if !half_length_m.is_finite() || half_length_m <= 0.0 || half_length_m > maximum {
            return Err(CrackGrowthError::CrackLengthOutsideDomain {
                length_m: half_length_m,
                minimum: 0.0,
                maximum,
            });
        }
        let angle = std::f64::consts::PI * (half_length_m / self.plate_width_m);
        let factor = 1.0 / det::sqrt(det::cos(angle));
        if !factor.is_finite() || factor <= 0.0 {
            return Err(CrackGrowthError::Arithmetic);
        }
        Ok(factor)
    }

    /// Sample exactly the declared ordered half-lengths. No knots are inserted,
    /// sorted, discarded, or extrapolated. The resulting table uses its usual
    /// log-log interpolation; callers choose resolution for their application.
    pub fn tabulate(
        &self,
        half_lengths_m: Vec<f64>,
        maximum_points: usize,
    ) -> Result<TabulatedCrackGeometry, CrackGrowthError> {
        self.tabulate_with_cancel(half_lengths_m, maximum_points, &mut || false)
    }

    /// Bounded construction which publishes no partial table on cancellation.
    pub fn tabulate_with_cancel(
        &self,
        half_lengths_m: Vec<f64>,
        maximum_points: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<TabulatedCrackGeometry, CrackGrowthError> {
        poll(cancelled)?;
        if half_lengths_m.len() > maximum_points {
            return Err(CrackGrowthError::Budget);
        }
        let mut points = Vec::with_capacity(half_lengths_m.len());
        for half_length_m in half_lengths_m {
            poll(cancelled)?;
            points.push(CrackGeometryPoint {
                crack_length_m: half_length_m,
                geometry_factor: self.factor_at(half_length_m)?,
            });
        }
        let context = CrackGeometryTableContext {
            source_reference: format!(
                "Finite-width centre-through-crack elastic secant model; full plate width W={:.17e} m; \
                 AFGROW DTD Handbook Table 11.3.8, p. 11.3.11 \
                 (https://www.afgrow.net/applications/DTDHandbook/pdfs/Sec11_3.pdf); \
                 supported 2a/W <= 0.8 from NASA TN D-5140, p. 6 \
                 (https://ntrs.nasa.gov/citations/19690012030)",
                self.plate_width_m
            ),
            crack_length_convention:
                "physical half-length a of a centred through crack; full crack length 2a, metres"
                    .into(),
            nominal_stress_basis:
                "remote uniform gross-section tensile stress normal to the crack; long elastic plate"
                    .into(),
        };
        let mut table =
            TabulatedCrackGeometry::new_with_cancel(context, points, maximum_points, cancelled)?;
        table.centre_crack = Some(*self);
        Ok(table)
    }

    /// Choose samples for a declared relative geometry-factor interpolation
    /// tolerance. Both endpoints are retained exactly; their interval must lie
    /// in this model's support. Work and storage are bounded by `maximum_points`,
    /// independent of the number of fatigue cycles subsequently integrated.
    ///
    /// In log length `x=ln(a)`, `f(x)=ln(Y)` is convex and has increasing
    /// curvature `f''=(theta*tan(theta)+theta^2*sec(theta)^2)/2`, where
    /// `theta=pi*a/W`. On a log interval of width `h`, the interpolation error
    /// is at most `h^2*f''(right)/8`. Refinement continues until that bound is
    /// at most `ln(1+tolerance)` for every interval. In real arithmetic the
    /// table therefore overestimates Y by no more than the requested fraction.
    ///
    /// This is an analytical interpolation bound evaluated in ordinary floating
    /// point, not an interval-arithmetic certificate. It does not bound the
    /// handbook approximation, material calibration, crack-growth life error,
    /// or applicability of this elastic geometry to a particular component.
    pub fn tabulate_to_tolerance(
        &self,
        minimum_half_length_m: f64,
        maximum_half_length_m: f64,
        relative_factor_tolerance: f64,
        maximum_points: usize,
    ) -> Result<TabulatedCrackGeometry, CrackGrowthError> {
        self.tabulate_to_tolerance_with_cancel(
            minimum_half_length_m,
            maximum_half_length_m,
            relative_factor_tolerance,
            maximum_points,
            &mut || false,
        )
    }

    /// Cancel between bounded refinements and during final table construction.
    /// Exhausted point budgets and indistinguishable midpoints refuse without
    /// returning a partially resolved table or relaxing the requested tolerance.
    pub fn tabulate_to_tolerance_with_cancel(
        &self,
        minimum_half_length_m: f64,
        maximum_half_length_m: f64,
        relative_factor_tolerance: f64,
        maximum_points: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<TabulatedCrackGeometry, CrackGrowthError> {
        poll(cancelled)?;
        self.factor_at(minimum_half_length_m)?;
        self.factor_at(maximum_half_length_m)?;
        if minimum_half_length_m >= maximum_half_length_m
            || !relative_factor_tolerance.is_finite()
            || relative_factor_tolerance <= 0.0
        {
            return Err(CrackGrowthError::Invalid(
                "declare an increasing half-length interval and positive finite relative Y tolerance",
            ));
        }
        if maximum_points < 2 {
            return Err(CrackGrowthError::Budget);
        }
        let log_tolerance = relative_factor_tolerance * log1p_over_x(relative_factor_tolerance);
        let mut lengths = vec![minimum_half_length_m];
        let mut pending = vec![maximum_half_length_m];
        let mut point_count = 2;
        let mut left = minimum_half_length_m;
        // Depth-first, left-to-right traversal: every split adds one unique
        // point, so both stored endpoints and visited intervals are O(budget).
        while let Some(right) = pending.pop() {
            poll(cancelled)?;
            let h = positive_log_ratio(left, right);
            let theta = std::f64::consts::PI * (right / self.plate_width_m);
            let (sine, cosine) = det::sin_cos(theta);
            let scaled_theta = theta / cosine;
            let curvature = 0.5 * (scaled_theta * sine + scaled_theta * scaled_theta);
            let bound = (h * h / 8.0) * curvature;
            if !h.is_finite() || h <= 0.0 || !bound.is_finite() || bound < 0.0 {
                return Err(CrackGrowthError::Arithmetic);
            }
            if bound <= log_tolerance {
                lengths.push(right);
                left = right;
                continue;
            }
            if point_count >= maximum_points {
                return Err(CrackGrowthError::Budget);
            }
            let midpoint = if h < 0.5 {
                left + left * det::expm1(0.5 * h)
            } else {
                det::exp(0.5 * det::ln(left) + 0.5 * det::ln(right))
            };
            if !midpoint.is_finite() || midpoint <= left || midpoint >= right {
                return Err(CrackGrowthError::Arithmetic);
            }
            point_count += 1;
            pending.push(right);
            pending.push(midpoint);
        }
        let mut table = self.tabulate_with_cancel(lengths, maximum_points, cancelled)?;
        table.context.source_reference.push_str(&format!(
            "; log-length curvature-bound refinement with requested relative Y tolerance={relative_factor_tolerance:.17e}"
        ));
        poll(cancelled)?;
        Ok(table)
    }
}
