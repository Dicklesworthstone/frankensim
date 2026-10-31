//! Finite-gap concentric Couette shear for a monotone generalized-Newtonian liquid.
//!
//! Angular momentum gives `tau(r) = tau_i (ri/r)^2`, and no slip requires
//! `abs(omega) = integral[ri,ro] gamma(tau(r)) dr/r`. A source viscosity line
//! is piecewise affine in physical shear rate: `mu = a + b gamma`. Therefore
//! `tau = gamma mu` has an analytic monotone inverse, and each contribution to
//! `2 abs(omega) = integral gamma d(tau)/tau` is integrated analytically.
//! Only the bounded scalar inner-wall stress solve remains; no radial mesh,
//! sampled effective viscosity, fitted power law or per-iteration allocation
//! is introduced. The original source owner supplies the complete viscosity
//! line and refreshes it at accepted thermal endpoints.
//!
//! This is steady, laminar, purely circumferential shear between concentric,
//! no-slip cylinders. Pressure support, end leakage, fluid inertia/storage,
//! yield surfaces, memory, thermal gradients and Taylor instabilities are
//! outside this owner. Its caller must admit the quasisteady creeping regime.

use crate::{TriboError, positive_finite};
use fs_math::det;

/// Explicit storage and evaluation bound on a source viscosity line.
pub const MAX_JOURNAL_VISCOSITY_KNOTS: usize = 256;

/// Actual concentric cylinder dimensions in SI units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JournalGeometry {
    pub inner_radius_m: f64,
    pub outer_radius_m: f64,
    pub length_m: f64,
}

impl JournalGeometry {
    /// Exact Newtonian geometric factor: torque = -factor * mu * omega.
    pub fn resistance_factor_m3(self) -> Result<f64, TriboError> {
        for value in [self.inner_radius_m, self.outer_radius_m, self.length_m] {
            positive_finite(value, "journal_geometry")?;
        }
        if self.outer_radius_m <= self.inner_radius_m {
            return invalid("journal_open_radial_gap");
        }
        // Factor the difference of squares, including for narrow clearances.
        let factor = 4.0 * core::f64::consts::PI * self.length_m
            * self.inner_radius_m.powi(2) * self.outer_radius_m.powi(2)
            / ((self.outer_radius_m - self.inner_radius_m)
                * (self.outer_radius_m + self.inner_radius_m));
        positive_finite(factor, "journal_geometry_resistance")?;
        Ok(factor)
    }
}

/// One endpoint of the original piecewise-affine dynamic-viscosity law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JournalViscosityKnot {
    pub shear_rate_per_s: f64,
    pub dynamic_viscosity_pa_s: f64,
}

/// Bounded scalar torque solve; residuals are physical angular rates [1/s].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JournalShearSolve {
    pub maximum_iterations: usize,
    pub absolute_angular_rate_tolerance_per_s: f64,
    pub relative_angular_rate_tolerance: f64,
}

/// One work-conjugate instantaneous native response, not integrated heat.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JournalShearResponse {
    /// Signed torque opposing the supplied relative angular velocity [N m].
    pub resisting_torque_nm: f64,
    pub dissipated_power_w: f64,
    /// Physical `sqrt(2 D:D)` at [inner, outer] walls, both nonnegative [1/s].
    pub shear_rate_per_s: [f64; 2],
    /// Signed reconstructed angular-speed magnitude minus the requested one.
    pub angular_rate_defect_per_s: f64,
    pub iterations: usize,
}

/// Compiled exact source line, frozen at one thermal/geometry configuration.
#[derive(Clone, Debug)]
pub struct GeneralizedNewtonianJournal {
    geometry: JournalGeometry,
    knots: Vec<JournalViscosityKnot>,
    stresses: Vec<f64>,
    slopes: Vec<f64>,
    solve: JournalShearSolve,
    factor: f64,
    constant_viscosity: Option<f64>,
    minimum_viscosity: f64,
}

fn invalid<T>(field: &'static str) -> Result<T, TriboError> {
    Err(TriboError::InvalidInput { field })
}

impl GeneralizedNewtonianJournal {
    /// Admit positive finite viscosity including zero shear, with a strictly
    /// increasing stress law over every COMPLETE source segment. No extension
    /// beyond the last supplied rate is implied. Shear thinning is admitted
    /// when d(gamma*mu)/d(gamma) remains positive at both segment endpoints.
    pub fn new(geometry: JournalGeometry, knots: Vec<JournalViscosityKnot>,
        solve: JournalShearSolve) -> Result<Self, TriboError> {
        let factor = geometry.resistance_factor_m3()?;
        if !(2..=MAX_JOURNAL_VISCOSITY_KNOTS).contains(&knots.len())
            || knots[0].shear_rate_per_s != 0.0
            || !(1..=128).contains(&solve.maximum_iterations)
            || !solve.absolute_angular_rate_tolerance_per_s.is_finite()
            || solve.absolute_angular_rate_tolerance_per_s <= 0.0
            || !solve.relative_angular_rate_tolerance.is_finite()
            || !(0.0..1.0).contains(&solve.relative_angular_rate_tolerance)
        { return invalid("journal_rheology_controls_or_zero_shear_support"); }
        let mut stresses = Vec::with_capacity(knots.len());
        let mut slopes = Vec::with_capacity(knots.len() - 1);
        let mut minimum_viscosity = f64::INFINITY;
        for (i, knot) in knots.iter().enumerate() {
            positive_finite(knot.dynamic_viscosity_pa_s, "journal_source_viscosity")?;
            if !knot.shear_rate_per_s.is_finite() || knot.shear_rate_per_s < 0.0 {
                return invalid("journal_source_shear_rate");
            }
            let stress = knot.shear_rate_per_s * knot.dynamic_viscosity_pa_s;
            if !stress.is_finite() { return invalid("journal_source_stress_overflow"); }
            minimum_viscosity = minimum_viscosity.min(knot.dynamic_viscosity_pa_s);
            if i > 0 {
                let previous = knots[i - 1];
                let delta = knot.shear_rate_per_s - previous.shear_rate_per_s;
                if delta <= 0.0 || stress <= stresses[i - 1] {
                    return invalid("journal_source_monotone_stress");
                }
                let slope = (knot.dynamic_viscosity_pa_s - previous.dynamic_viscosity_pa_s) / delta;
                for endpoint in [previous, *knot] {
                    positive_finite(endpoint.dynamic_viscosity_pa_s
                        + slope * endpoint.shear_rate_per_s, "journal_source_stress_tangent")?;
                }
                slopes.push(slope);
            }
            stresses.push(stress);
        }
        let constant_viscosity = knots.iter().all(|p|
            p.dynamic_viscosity_pa_s == knots[0].dynamic_viscosity_pa_s)
            .then_some(knots[0].dynamic_viscosity_pa_s);
        Ok(Self { geometry, knots, stresses, slopes, solve, factor,
            constant_viscosity, minimum_viscosity })
    }

    /// Minimum over the complete admitted rate interval, for conservative Re checks.
    pub const fn minimum_viscosity_pa_s(&self) -> f64 { self.minimum_viscosity }
    pub fn knots(&self) -> &[JournalViscosityKnot] { &self.knots }

    /// The same admitted affine source line, useful for coefficient-change gates.
    pub fn viscosity_at_shear_rate(&self, rate: f64) -> Result<f64, TriboError> {
        if !rate.is_finite() || rate < 0.0 {
            return invalid("journal_source_shear_rate");
        }
        let hi = self.knots.partition_point(|point| point.shear_rate_per_s < rate);
        if hi == self.knots.len() { return invalid("journal_source_shear_rate_exceeded"); }
        if rate == self.knots[hi].shear_rate_per_s { return Ok(self.knots[hi].dynamic_viscosity_pa_s); }
        let lo = hi - 1;
        let value = self.knots[lo].dynamic_viscosity_pa_s
            + self.slopes[lo] * (rate - self.knots[lo].shear_rate_per_s);
        positive_finite(value, "journal_interpolated_source_viscosity")?;
        Ok(value)
    }

    /// Evaluate the actual nonlinear torque and its SAME mechanical power.
    /// The immutable line is shared by every bounded bisection evaluation.
    pub fn response(&self, omega: f64) -> Result<JournalShearResponse, TriboError> {
        if !omega.is_finite() { return invalid("journal_angular_velocity"); }
        if omega == 0.0 {
            return Ok(JournalShearResponse { resisting_torque_nm: 0.0,
                dissipated_power_w: 0.0, shear_rate_per_s: [0.0; 2],
                angular_rate_defect_per_s: 0.0, iterations: 0 });
        }
        let g = self.geometry;
        let ratio = (g.inner_radius_m / g.outer_radius_m).powi(2);
        let speed = omega.abs();
        if let Some(mu) = self.constant_viscosity {
            let inner = 2.0 * g.outer_radius_m.powi(2) * speed
                / ((g.outer_radius_m - g.inner_radius_m) * (g.outer_radius_m + g.inner_radius_m));
            if !inner.is_finite() || inner > self.knots.last().expect("admitted line").shear_rate_per_s {
                return invalid("journal_source_shear_rate_exceeded");
            }
            return self.finish(omega, inner * mu, [inner, inner * ratio], 0.0, 0,
                Some(-self.factor * mu * omega));
        }
        let mut lower = 0.0;
        let mut upper = *self.stresses.last().expect("admitted line");
        let (_, maximum) = self.angular_speed(upper, ratio)?;
        if speed > maximum { return invalid("journal_source_shear_rate_exceeded"); }
        let tolerance = self.solve.absolute_angular_rate_tolerance_per_s
            + self.solve.relative_angular_rate_tolerance * speed;
        for iteration in 1..=self.solve.maximum_iterations {
            let stress = lower + 0.5 * (upper - lower);
            let (rates, actual) = self.angular_speed(stress, ratio)?;
            let defect = actual - speed;
            if defect.abs() <= tolerance {
                return self.finish(omega, stress, rates, defect, iteration, None);
            }
            if actual < speed { lower = stress; } else { upper = stress; }
        }
        invalid("journal_angular_rate_solve_exhausted")
    }

    fn finish(&self, omega: f64, stress: f64, rates: [f64; 2], defect: f64,
        iterations: usize, exact_torque: Option<f64>) -> Result<JournalShearResponse, TriboError> {
        let torque = exact_torque.unwrap_or_else(|| -omega.signum() * 2.0
            * core::f64::consts::PI * self.geometry.length_m
            * self.geometry.inner_radius_m.powi(2) * stress);
        let power = -torque * omega;
        if !torque.is_finite() || !power.is_finite() || power <= 0.0 {
            return invalid("journal_torque_or_power_unrepresentable");
        }
        Ok(JournalShearResponse { resisting_torque_nm: torque, dissipated_power_w: power,
            shear_rate_per_s: rates, angular_rate_defect_per_s: defect, iterations })
    }

    fn shear_rate(&self, stress: f64) -> Result<f64, TriboError> {
        if stress == 0.0 { return Ok(0.0); }
        let hi = self.stresses.partition_point(|value| *value < stress);
        if hi == self.stresses.len() { return invalid("journal_source_stress_exceeded"); }
        if stress == self.stresses[hi] { return Ok(self.knots[hi].shear_rate_per_s); }
        let lo = hi - 1;
        let p = self.knots[lo];
        let slope = self.slopes[lo];
        let tangent = p.dynamic_viscosity_pa_s + slope * p.shear_rate_per_s;
        let delta_stress = stress - self.stresses[lo];
        let discriminant = tangent * tangent + 4.0 * slope * delta_stress;
        let delta = 2.0 * delta_stress / (tangent + det::sqrt(discriminant));
        let rate = p.shear_rate_per_s + delta;
        if !rate.is_finite() || rate < p.shear_rate_per_s || rate > self.knots[hi].shear_rate_per_s {
            return invalid("journal_source_inverse_unrepresentable");
        }
        Ok(rate)
    }

    fn angular_speed(&self, inner_stress: f64, ratio: f64)
        -> Result<([f64; 2], f64), TriboError> {
        let inner = self.shear_rate(inner_stress)?;
        let outer = self.shear_rate(inner_stress * ratio)?;
        let mut integral = 0.0;
        for (i, pair) in self.knots.windows(2).enumerate() {
            let a = outer.max(pair[0].shear_rate_per_s);
            let b = inner.min(pair[1].shear_rate_per_s);
            if b <= a { continue; }
            let delta = b - a;
            let slope = self.slopes[i];
            let mu_a = pair[0].dynamic_viscosity_pa_s
                + slope * (a - pair[0].shear_rate_per_s);
            let x = slope * delta / mu_a;
            let ratio_log = log1p_over_x(x);
            // Equivalent to 2Δγ-(intercept/slope) log(mu_b/mu_a),
            // without dividing by a tiny constitutive slope.
            let value = delta * (2.0 - ratio_log) + a * x * ratio_log;
            if !value.is_finite() || value < 0.0 {
                return invalid("journal_source_integral_unrepresentable");
            }
            integral += value;
        }
        let speed = 0.5 * integral;
        if !speed.is_finite() { return invalid("journal_angular_rate_overflow"); }
        Ok(([inner, outer], speed))
    }
}

// log(1+x)/x, using a fixed convergent series near zero so constant-viscosity
// limits retain relative accuracy. The outer logarithm uses the shared det owner.
fn log1p_over_x(x: f64) -> f64 {
    if x.abs() < 0.125 {
        let mut sum = 1.0;
        let mut power = 1.0;
        for n in 1..=20 {
            power *= -x;
            sum += power / f64::from(n + 1);
        }
        sum
    } else {
        det::ln(1.0 + x) / x
    }
}

#[cfg(test)]
mod tests;
