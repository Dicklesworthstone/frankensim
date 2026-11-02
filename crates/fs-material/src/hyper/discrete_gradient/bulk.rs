//! Condensed original bulk storage for conservative mixed elasticity.
//!
//! With the accepted local volume `a` fixed, interval pressure is the scalar
//! discrete derivative `p = (U(b)-U(a))/(b-a)`. Its derivative with respect to
//! `b` is `integral_0^1 t U''(a+t*(b-a)) dt`, so the inverse is unique on a
//! strictly convex volume bracket. Quadratic bulk is inverted analytically;
//! log-squared bulk retains its original energy on `0 < J < exp(1)`.
//!
//! This pressure is an interval work multiplier, not endpoint `U'(b)`. A mixed
//! spatial owner retains `a` at each quadrature point and enforces its weak
//! endpoint constraint against `b`. The same P1 interval pressure then pairs
//! the geometric volume change with the original bulk-energy change. No
//! quadratic replacement or pressure-midpoint approximation is introduced.
//! Initial retained volumes belong to a compatible mixed endpoint state; they
//! must not be silently reset to geometric det(F) at every accepted interval.
//!
//! Scalar work is explicitly bounded. Cancellation is polled by the spatial
//! owner at its quadrature boundary; one query performs at most 128 iterations.

use fs_math::det;

use super::super::mixed::{MixedFiniteStrainPotential, original_volume, pressure_potential};
use super::super::{Hyperelastic, Volumetric};
use crate::MaterialError;

/// Hard upper bound on one local scalar solve; the caller selects its own cap.
pub const MAX_BULK_SECANT_ITERATIONS: usize = 128;

/// Explicit original-volume branch, scalar work cap and local error allowances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DiscreteBulkSettings {
    pub minimum_volume_ratio: f64,
    pub maximum_volume_ratio: f64,
    pub maximum_iterations: usize,
    /// Absolute residual of the scalar pressure equation [Pa].
    pub absolute_pressure_tolerance_pa: f64,
    /// Additional allowance times max(|requested pressure|, |secant pressure|).
    pub relative_pressure_tolerance: f64,
    /// Absolute defect of pressure work versus original bulk-energy change [Pa].
    pub absolute_work_tolerance_pa: f64,
    /// Additional allowance times max(|U_before|, |U_after|, |pressure work|).
    pub relative_work_tolerance: f64,
}

impl DiscreteBulkSettings {
    fn validate(self, volume: Volumetric) -> Result<(), MaterialError> {
        if !self.minimum_volume_ratio.is_finite()
            || self.minimum_volume_ratio <= 0.0
            || !self.maximum_volume_ratio.is_finite()
            || self.maximum_volume_ratio <= self.minimum_volume_ratio
            || self.maximum_iterations == 0
            || self.maximum_iterations > MAX_BULK_SECANT_ITERATIONS
            || ![
                self.absolute_pressure_tolerance_pa,
                self.relative_pressure_tolerance,
                self.absolute_work_tolerance_pa,
                self.relative_work_tolerance,
            ]
            .iter()
            .all(|v| v.is_finite() && *v >= 0.0)
            || (matches!(volume, Volumetric::LogSquared { .. })
                && det::ln(self.maximum_volume_ratio) >= 1.0)
        {
            return Err(MaterialError::Parameters {
                what: "discrete bulk requires a positive ordered convex volume bracket, 1..=128 iterations, and finite nonnegative allowances".into(),
            });
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DiscreteBulkResponse {
    pub volume_ratio_after: f64,
    /// Original physical endpoint bulk storage per reference volume [Pa].
    pub energy_after_pa: f64,
    /// Original energy difference, evaluated without subtracting nearby squares.
    pub energy_change_pa: f64,
    /// Original endpoint derivative U'(J_after), not the interval multiplier.
    pub endpoint_pressure_pa: f64,
    /// dJ_after / d(interval pressure), positive on the admitted branch [Pa^-1].
    pub compliance_per_pa: f64,
    /// Computed secant pressure minus the requested interval pressure [Pa].
    pub pressure_residual_pa: f64,
    pub work_pa: f64,
    pub work_defect_pa: f64,
    /// Zero for the analytic quadratic inverse or an unchanged exact endpoint.
    pub iterations: usize,
}

/// Immutable accepted quadrature volume and the original, unmodified bulk law.
#[derive(Clone, Copy, Debug)]
pub struct DiscreteBulkPoint {
    volume: Volumetric,
    accepted_volume_ratio: f64,
    energy_before_pa: f64,
    settings: DiscreteBulkSettings,
}

fn state_error(what: impl Into<String>) -> MaterialError {
    MaterialError::State {
        what: format!("conservative original bulk: {}", what.into()),
    }
}

// The finite series resolves coalescing logarithms without platform libm or
// subtraction of two separately rounded logarithms. |r|<=1/8 implies |s|<=1/15;
// the omitted atanh tail after s^19 is below 1e-25 relative.
fn small_log_ratio(r: f64) -> f64 {
    let s = r / (2.0 + r);
    let mut power = s;
    let mut sum = s;
    for n in 1..10 {
        power *= s * s;
        sum += power / f64::from(2 * n + 1);
    }
    2.0 * sum
}

// Normalized log-squared secant and its derivative with respect to b.
// The near-pair series is
// (1/a) sum_n (-1)^n (ln(a)-H_n) r^n/(n+1), r=(b-a)/a.
// Its differentiated series retains U''(a)/2 at equality without subtracting
// nearly equal stresses. The omitted derivative is O(r^20), below roundoff
// on this branch; the independently evaluated work identity is still gated.
fn log_secant(a: f64, b: f64) -> (f64, f64, f64) {
    let delta = b - a;
    let l0 = det::ln(a);
    let l1 = det::ln(b);
    let r = delta / a;
    if r.abs() <= 0.125 {
        let mut coefficients = [0.0; 21];
        let mut harmonic = 0.0;
        for (n, coefficient) in coefficients.iter_mut().enumerate() {
            if n > 0 {
                harmonic += 1.0 / n as f64;
            }
            let sign = if n % 2 == 0 { 1.0 } else { -1.0 };
            *coefficient = sign * (l0 - harmonic) / (n + 1) as f64;
        }
        let mut value = coefficients[20];
        let mut derivative = 0.0_f64;
        for coefficient in coefficients[..20].iter().rev() {
            derivative = derivative.mul_add(r, value);
            value = value.mul_add(r, *coefficient);
        }
        let dl = small_log_ratio(r);
        (value / a, derivative / a / a, (l0 + 0.5 * dl) * dl)
    } else {
        let change = 0.5 * (l1 + l0) * (l1 - l0);
        let value = change / delta;
        (value, (l1 / b - value) / delta, change)
    }
}

impl DiscreteBulkPoint {
    /// Initialize a compatible mixed bulk state from ORIGINAL endpoint pressure.
    /// The spatial owner still admits its weak geometric-volume constraint.
    pub fn from_endpoint_pressure(
        law: &Hyperelastic,
        endpoint_pressure_pa: f64,
        settings: DiscreteBulkSettings,
    ) -> Result<Self, MaterialError> {
        let volume = original_volume(law)?;
        let initial = pressure_potential(volume, endpoint_pressure_pa)?;
        Self::try_new(law, initial.volume_ratio, settings)
    }

    pub fn try_new(
        law: &Hyperelastic,
        accepted_volume_ratio: f64,
        settings: DiscreteBulkSettings,
    ) -> Result<Self, MaterialError> {
        law.validate_mixed()?;
        let volume = original_volume(law)?;
        settings.validate(volume)?;
        let point = Self {
            volume,
            accepted_volume_ratio,
            energy_before_pa: 0.0,
            settings,
        };
        point.rebind(accepted_volume_ratio)
    }

    pub fn accepted_volume_ratio(&self) -> f64 {
        self.accepted_volume_ratio
    }

    pub fn energy_before_pa(&self) -> f64 {
        self.energy_before_pa
    }

    /// Original endpoint U'(J_before), useful as a physical pressure predictor.
    pub fn endpoint_pressure_pa(&self) -> f64 {
        self.endpoint_pressure(self.accepted_volume_ratio)
    }

    /// Reuse the already admitted original law and controls at a new endpoint.
    /// This creates a candidate; the accepted point is never modified.
    pub fn rebind(&self, volume_ratio: f64) -> Result<Self, MaterialError> {
        self.check_volume(volume_ratio)?;
        let energy_before_pa = self.volume.energy(volume_ratio);
        if !energy_before_pa.is_finite() || !self.endpoint_pressure(volume_ratio).is_finite() {
            return Err(state_error(
                "endpoint bulk storage or pressure is unrepresentable",
            ));
        }
        Ok(Self {
            accepted_volume_ratio: volume_ratio,
            energy_before_pa,
            ..*self
        })
    }

    fn check_volume(&self, value: f64) -> Result<(), MaterialError> {
        if !value.is_finite()
            || value < self.settings.minimum_volume_ratio
            || value > self.settings.maximum_volume_ratio
        {
            Err(state_error(
                "volume is outside the explicitly admitted branch",
            ))
        } else {
            Ok(())
        }
    }

    fn endpoint_pressure(&self, j: f64) -> f64 {
        self.volume.kappa()
            * match self.volume {
                Volumetric::Quadratic { .. } => j - 1.0,
                Volumetric::LogSquared { .. } => det::ln(j) / j,
            }
    }

    // All entries normalized by kappa: secant, d(secant)/dJ_after and Delta U.
    fn normalized_pair(&self, b: f64) -> (f64, f64, f64) {
        let a = self.accepted_volume_ratio;
        match self.volume {
            Volumetric::Quadratic { .. } => {
                let secant = 0.5 * (a - 1.0) + 0.5 * (b - 1.0);
                (secant, 0.5, secant * (b - a))
            }
            Volumetric::LogSquared { .. } => log_secant(a, b),
        }
    }

    fn candidate(
        &self,
        pressure_pa: f64,
        b: f64,
        iterations: usize,
    ) -> Result<(DiscreteBulkResponse, bool), MaterialError> {
        self.check_volume(b)?;
        let (secant, derivative, change) = self.normalized_pair(b);
        let kappa = self.volume.kappa();
        let secant_pa = kappa * secant;
        let response = DiscreteBulkResponse {
            volume_ratio_after: b,
            energy_after_pa: self.volume.energy(b),
            energy_change_pa: kappa * change,
            endpoint_pressure_pa: self.endpoint_pressure(b),
            compliance_per_pa: (1.0 / kappa) / derivative,
            pressure_residual_pa: secant_pa - pressure_pa,
            work_pa: pressure_pa * (b - self.accepted_volume_ratio),
            work_defect_pa: 0.0,
            iterations,
        };
        let response = DiscreteBulkResponse {
            work_defect_pa: response.work_pa - response.energy_change_pa,
            ..response
        };
        let pressure_allowance = self.settings.absolute_pressure_tolerance_pa
            + self.settings.relative_pressure_tolerance * pressure_pa.abs().max(secant_pa.abs());
        let work_allowance = self.settings.absolute_work_tolerance_pa
            + self.settings.relative_work_tolerance
                * self
                    .energy_before_pa
                    .abs()
                    .max(response.energy_after_pa.abs())
                    .max(response.work_pa.abs());
        if ![
            response.energy_after_pa,
            response.energy_change_pa,
            response.endpoint_pressure_pa,
            response.compliance_per_pa,
            response.pressure_residual_pa,
            response.work_pa,
            response.work_defect_pa,
            pressure_allowance,
            work_allowance,
        ]
        .iter()
        .all(|v| v.is_finite())
            || !derivative.is_finite()
            || derivative <= 0.0
            || response.compliance_per_pa <= 0.0
        {
            return Err(state_error(
                "secant, positive compliance or local work is unrepresentable",
            ));
        }
        let admitted = response.pressure_residual_pa.abs() <= pressure_allowance
            && response.work_defect_pa.abs() <= work_allowance;
        Ok((response, admitted))
    }

    /// Condense the endpoint volume against one finite interval pressure.
    /// The branch and work limits are checked without mutating this point.
    pub fn response(&self, pressure_pa: f64) -> Result<DiscreteBulkResponse, MaterialError> {
        if !pressure_pa.is_finite() {
            return Err(state_error("interval pressure must be finite"));
        }
        if matches!(self.volume, Volumetric::Quadratic { .. }) {
            let b = 1.0
                + (2.0 * (pressure_pa / self.volume.kappa()) - (self.accepted_volume_ratio - 1.0));
            let (result, admitted) = self.candidate(pressure_pa, b, 0)?;
            return if admitted {
                Ok(result)
            } else {
                Err(state_error(
                    "analytic quadratic inverse exceeds pressure or work allowance",
                ))
            };
        }
        if pressure_pa == self.endpoint_pressure_pa() {
            let (result, admitted) = self.candidate(pressure_pa, self.accepted_volume_ratio, 0)?;
            if admitted {
                return Ok(result);
            }
        }
        let (mut lo, mut hi) = (
            self.settings.minimum_volume_ratio,
            self.settings.maximum_volume_ratio,
        );
        let (low, low_ok) = self.candidate(pressure_pa, lo, 0)?;
        let (high, high_ok) = self.candidate(pressure_pa, hi, 0)?;
        if low.pressure_residual_pa > 0.0 || high.pressure_residual_pa < 0.0 {
            return Err(state_error(
                "interval pressure has no endpoint in the admitted volume bracket",
            ));
        }
        if low_ok {
            return Ok(low);
        }
        if high_ok {
            return Ok(high);
        }
        let (_, derivative, _) = self.normalized_pair(self.accepted_volume_ratio);
        let predictor = self.accepted_volume_ratio
            + (pressure_pa / self.volume.kappa()
                - self.endpoint_pressure_pa() / self.volume.kappa())
                / derivative;
        let mut b = if predictor > lo && predictor < hi {
            predictor
        } else {
            0.5 * lo + 0.5 * hi
        };
        for iteration in 1..=self.settings.maximum_iterations {
            let (result, admitted) = self.candidate(pressure_pa, b, iteration)?;
            if admitted {
                return Ok(result);
            }
            if result.pressure_residual_pa < 0.0 {
                lo = b;
            } else {
                hi = b;
            }
            let next = b - result.pressure_residual_pa * result.compliance_per_pa;
            b = if next > lo && next < hi {
                next
            } else {
                0.5 * lo + 0.5 * hi
            };
        }
        Err(state_error(
            "scalar inverse exhausted its pressure/work iteration budget",
        ))
    }
}
