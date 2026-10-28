//! Bounded single-liquid Cox contact-line closure with an inviscid outer gas.
//!
//! Chan et al., JFM 900 A8 (2020), equations (1.1)-(1.3):
//! https://doi.org/10.1017/jfm.2020.499
//! g(theta_app) = g(theta_micro) + Ca ln(L/l),
//! g(theta) = integral_0^theta (u - sin(u) cos(u))/(2 sin(u)) du.
//!
//! The full angular function avoids applying the small-angle cubic law to
//! hydrophobic surfaces. Effective cutoff l and observation scale L are
//! explicit model inputs; l is not silently identified with a measured slip
//! length. Newtonian creeping flow, negligible gas viscosity and separated
//! scales remain assumptions. This closure does not resolve wall slip,
//! nucleation, film entrainment, or a conservative surface-energy transfer.

use std::sync::Arc;

use fs_exec::CancelGate;

use super::{Error, checkpoint};

/// Explicit physical scales and admitted model domain. No universal defaults.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoxVoinovControls {
    /// Effective microscopic matching cutoff [m], strictly positive.
    pub microscopic_cutoff_m: f64,
    /// Scale at which the apparent boundary angle is imposed [m], greater than l.
    pub observation_length_m: f64,
    /// Admitted |mu U / sigma|, in (0, 0.01]. A model bound, not calibration.
    pub maximum_capillary_number: f64,
    /// Lower admitted microscopic AND apparent angle [rad], at least 1e-4.
    pub minimum_angle_radians: f64,
    /// Upper admitted microscopic AND apparent angle [rad], at most 3.
    pub maximum_angle_radians: f64,
}

/// Declared hydrodynamic model and shared, fixed-order angular quadrature table.
/// Clones share the table; accepted VOF steps do not rebuild or integrate it.
#[derive(Clone, Debug)]
pub struct CoxVoinovLaw {
    controls: CoxVoinovControls,
    logarithmic_scale_ratio: f64,
    table: Arc<AngularTable>,
}

/// Physical coefficients actually used in one model prediction.
/// The endpoint material receipt remains a STATIC angle measurement.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoxVoinovSample {
    /// Signed Ca = mu U / sigma; positive denotes advancing liquid.
    pub signed_capillary_number: f64,
    /// Current source viscosity [Pa s], including thermal evolution.
    pub dynamic_viscosity_pa_s: f64,
    /// The native owner's admitted constant surface tension [N/m].
    pub surface_tension_n_m: f64,
}

const INTERVALS: usize = 2048;
const MAX_ANGLE: f64 = 3.0;
const STEP: f64 = MAX_ANGLE / INTERVALS as f64;

#[derive(Debug)]
struct AngularTable {
    values: Vec<f64>,
    slopes: Vec<f64>,
}

impl CoxVoinovLaw {
    /// Admit the model and prepare its bounded angular table once.
    /// Cancellation leaves no installed boundary or runtime history.
    pub fn try_new(controls: CoxVoinovControls, gate: &CancelGate) -> Result<Self, Error> {
        checkpoint(gate)?;
        if !controls.microscopic_cutoff_m.is_finite()
            || controls.microscopic_cutoff_m <= 0.0
            || !controls.observation_length_m.is_finite()
            || controls.observation_length_m <= controls.microscopic_cutoff_m
            || !controls.maximum_capillary_number.is_finite()
            || controls.maximum_capillary_number <= 0.0
            || controls.maximum_capillary_number > 0.01
            || !controls.minimum_angle_radians.is_finite()
            || controls.minimum_angle_radians < 1e-4
            || !controls.maximum_angle_radians.is_finite()
            || controls.maximum_angle_radians > MAX_ANGLE
            || controls.minimum_angle_radians >= controls.maximum_angle_radians
        {
            return Err(Error::Invalid(
                "invalid Cox scales, low-capillary-number or angular domain",
            ));
        }
        // Difference of logs avoids overflow of the dimensional ratio.
        let logarithmic_scale_ratio = fs_math::det::ln(controls.observation_length_m)
            - fs_math::det::ln(controls.microscopic_cutoff_m);
        if !logarithmic_scale_ratio.is_finite() || logarithmic_scale_ratio <= 0.0 {
            return Err(Error::Invalid("Cox scale separation is unrepresentable"));
        }
        let table = Arc::new(AngularTable::try_new(gate)?);
        Ok(Self {
            controls,
            logarithmic_scale_ratio,
            table,
        })
    }

    /// Complete caller-declared model inputs.
    #[must_use]
    pub const fn controls(&self) -> &CoxVoinovControls {
        &self.controls
    }

    pub(super) fn admit_angle(&self, angle: f64) -> Result<(), Error> {
        if !angle.is_finite()
            || angle < self.controls.minimum_angle_radians
            || angle > self.controls.maximum_angle_radians
        {
            return Err(Error::Invalid(
                "Cox contact angle leaves its declared angular domain",
            ));
        }
        Ok(())
    }

    pub(super) fn admit_capillary_number(&self, capillary_number: f64) -> Result<(), Error> {
        if !capillary_number.is_finite()
            || capillary_number.abs() > self.controls.maximum_capillary_number
        {
            return Err(Error::Invalid(
                "Cox contact speed leaves its low-capillary-number domain",
            ));
        }
        Ok(())
    }

    /// Apparent angle from a microscopic endpoint and SIGNED capillary number.
    /// No clipping at zero/pi or extrapolation through a film-transition limit.
    /// The fixed table uses cubic Hermite interpolation and safeguarded inversion;
    /// there is no per-call allocation, quadrature or transcendental evaluation.
    pub fn apparent_angle_radians(
        &self,
        microscopic_angle_radians: f64,
        signed_capillary_number: f64,
    ) -> Result<f64, Error> {
        self.admit_angle(microscopic_angle_radians)?;
        self.admit_capillary_number(signed_capillary_number)?;
        if signed_capillary_number == 0.0 {
            return Ok(microscopic_angle_radians);
        }
        let target = self.table.value(microscopic_angle_radians)
            + signed_capillary_number * self.logarithmic_scale_ratio;
        let lower = self.table.value(self.controls.minimum_angle_radians);
        let upper = self.table.value(self.controls.maximum_angle_radians);
        if !target.is_finite() || target < lower || target > upper {
            return Err(Error::Invalid(
                "Cox moving contact line has no angle in its admitted domain",
            ));
        }
        let angle = self.table.inverse(target);
        // The final roundoff-sized bracket must not widen the declared domain.
        Ok(angle.clamp(
            self.controls.minimum_angle_radians,
            self.controls.maximum_angle_radians,
        ))
    }
}

impl AngularTable {
    fn try_new(gate: &CancelGate) -> Result<Self, Error> {
        let mut values = Vec::with_capacity(INTERVALS + 1);
        let mut slopes = Vec::with_capacity(INTERVALS + 1);
        values.push(0.0);
        slopes.push(0.0);
        let (mut sum, mut correction) = (0.0, 0.0);
        // Four-point Gauss-Legendre quadrature on each fixed interval.
        for i in 0..INTERVALS {
            if i % 64 == 0 {
                checkpoint(gate)?;
            }
            let midpoint = (i as f64 + 0.5) * STEP;
            let half = 0.5 * STEP;
            let integral = half
                * (0.34785484513745385
                    * (angular_slope(midpoint - 0.8611363115940526 * half)
                        + angular_slope(midpoint + 0.8611363115940526 * half))
                    + 0.6521451548625461
                        * (angular_slope(midpoint - 0.3399810435848563 * half)
                            + angular_slope(midpoint + 0.3399810435848563 * half)));
            let adjusted = integral - correction;
            let next = sum + adjusted;
            correction = (next - sum) - adjusted;
            sum = next;
            values.push(sum);
            slopes.push(angular_slope((i + 1) as f64 * STEP));
        }
        checkpoint(gate)?;
        Ok(Self { values, slopes })
    }

    fn interpolate(&self, i: usize, t: f64) -> f64 {
        let y0 = self.values[i];
        let delta = self.values[i + 1] - y0;
        let m0 = STEP * self.slopes[i];
        let m1 = STEP * self.slopes[i + 1];
        y0 + t * (m0 + t * (3.0 * delta - 2.0 * m0 - m1 + t * (m0 + m1 - 2.0 * delta)))
    }

    fn value(&self, angle: f64) -> f64 {
        let index = ((angle / STEP) as usize).min(INTERVALS - 1);
        self.interpolate(index, (angle - index as f64 * STEP) / STEP)
    }

    fn inverse(&self, target: f64) -> f64 {
        let index = self
            .values
            .partition_point(|&value| value < target)
            .saturating_sub(1)
            .min(INTERVALS - 1);
        let (mut lower, mut upper) = (0.0, 1.0);
        // Initial binary search locates a single cell. Thirty-six bounded
        // arithmetic-only bisections then resolve angle to about 2.2e-14 rad.
        for _ in 0..36 {
            let middle = f64::midpoint(lower, upper);
            if self.interpolate(index, middle) < target {
                lower = middle;
            } else {
                upper = middle;
            }
        }
        (index as f64 + f64::midpoint(lower, upper)) * STEP
    }
}

fn angular_slope(angle: f64) -> f64 {
    if angle < 0.01 {
        // Removable singularity: avoid cancellation in u - sin(u) cos(u).
        let square = angle * angle;
        square * (1.0 / 3.0 + square * (-1.0 / 90.0 + square * (13.0 / 7560.0 + square / 10800.0)))
    } else {
        let (sin, cos) = fs_math::det::sin_cos(angle);
        (angle - sin * cos) / (2.0 * sin)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn controls() -> CoxVoinovControls {
        CoxVoinovControls {
            microscopic_cutoff_m: 1e-9,
            observation_length_m: 1e-3,
            maximum_capillary_number: 0.01,
            minimum_angle_radians: 1e-4,
            maximum_angle_radians: 3.0,
        }
    }

    // Independent composite Simpson quadrature, using platform trig rather
    // than the production fixed Gauss/Hermite table.
    fn reference_g(angle: f64) -> f64 {
        let n = 16_384;
        let h = angle / n as f64;
        let slope = |u: f64| {
            if u == 0.0 {
                0.0
            } else {
                (u - u.sin() * u.cos()) / (2.0 * u.sin())
            }
        };
        let mut result = slope(angle);
        for i in 1..n {
            result += (if i % 2 == 0 { 2.0 } else { 4.0 }) * slope(i as f64 * h);
        }
        result * h / 3.0
    }

    #[test]
    fn g1_full_angular_law_matches_independent_quadrature_and_direction() {
        let law = CoxVoinovLaw::try_new(controls(), &CancelGate::new_clock_free()).unwrap();
        for microscopic in [0.3, 1.2, 2.4, 2.9] {
            assert_eq!(
                law.apparent_angle_radians(microscopic, 0.0).unwrap(),
                microscopic
            );
            for capillary in [-0.0001, 0.0001] {
                let apparent = law.apparent_angle_radians(microscopic, capillary).unwrap();
                assert_eq!(apparent > microscopic, capillary > 0.0);
                let expected = capillary * (1e6_f64).ln();
                let actual = reference_g(apparent) - reference_g(microscopic);
                assert!(
                    (actual - expected).abs() < 2e-10,
                    "{microscopic} {capillary}: {actual}"
                );
            }
        }
        let small: f64 = 1e-3;
        assert!((law.table.value(small) / (small.powi(3) / 9.0) - 1.0).abs() < 1e-7);
        // At a hydrophobic angle the full angular relation is observably
        // different from the unrestricted theta^3 approximation.
        let microscopic: f64 = 2.9;
        let apparent = law.apparent_angle_radians(microscopic, 0.001).unwrap();
        let cubic = (microscopic.powi(3) + 9.0 * 0.001 * (1e6_f64).ln()).cbrt();
        assert!((apparent - cubic).abs() > 1e-4);
    }

    #[test]
    fn g3_scale_similarity_and_monotone_inversion() {
        let gate = CancelGate::new_clock_free();
        let a = CoxVoinovLaw::try_new(controls(), &gate).unwrap();
        let b = CoxVoinovLaw::try_new(
            CoxVoinovControls {
                microscopic_cutoff_m: 1e-6,
                observation_length_m: 1.0,
                ..controls()
            },
            &gate,
        )
        .unwrap();
        let mut previous = 0.0;
        for i in 0..101 {
            let capillary = -0.005 + 0.0001 * i as f64;
            let angle = a.apparent_angle_radians(1.5, capillary).unwrap();
            assert!(angle > previous);
            previous = angle;
            assert!((angle - b.apparent_angle_radians(1.5, capillary).unwrap()).abs() < 1e-12);
        }
    }

    #[test]
    fn g0_g4_invalid_scales_high_ca_dewetting_and_cancellation_refuse() {
        let gate = CancelGate::new_clock_free();
        let a = CoxVoinovLaw::try_new(controls(), &gate).unwrap();
        assert!(a.apparent_angle_radians(1.5, 0.0101).is_err());
        assert!(a.apparent_angle_radians(0.1, -0.001).is_err());
        assert!(a.apparent_angle_radians(3.0, 0.001).is_err());
        assert!(a.apparent_angle_radians(f64::NAN, 0.0).is_err());
        for invalid in [
            CoxVoinovControls {
                microscopic_cutoff_m: 0.0,
                ..controls()
            },
            CoxVoinovControls {
                observation_length_m: 1e-10,
                ..controls()
            },
            CoxVoinovControls {
                maximum_capillary_number: 0.1,
                ..controls()
            },
            CoxVoinovControls {
                maximum_angle_radians: core::f64::consts::PI,
                ..controls()
            },
        ] {
            assert!(CoxVoinovLaw::try_new(invalid, &gate).is_err());
        }
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(CoxVoinovLaw::try_new(controls(), &cancelled).is_err());
    }
}
