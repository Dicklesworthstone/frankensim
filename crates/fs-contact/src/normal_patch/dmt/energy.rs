//! Stable exact-piece work and the native discrete-gradient energy numerator.
use super::{DmtError, SourceDmtSphereContact};

impl SourceDmtSphereContact {
    /// Conservative work for prescribed signed approach, J. Integrates the
    /// original force analytically, avoiding cancellation between endpoint
    /// potentials for small motion. No intrinsic hysteresis or heat is added.
    pub fn reversible_work_between(&self, from_m: f64, to_m: f64) -> Result<f64, DmtError> {
        self.evaluate(from_m)?;
        self.evaluate(to_m)?;
        if from_m == to_m {
            return Ok(0.0);
        }
        let (lo, hi, sign) = if from_m < to_m {
            (from_m, to_m, 1.0)
        } else {
            (to_m, from_m, -1.0)
        };
        let a = lo.max(-self.interaction_range_m);
        let b = hi.min(0.0);
        let gap_work = if b > a {
            -self.tail_stiffness_n_per_m
                * (b - a)
                * (0.5 * (a + self.interaction_range_m) + 0.5 * (b + self.interaction_range_m))
        } else {
            0.0
        };
        let a = lo.max(0.0);
        let b = hi.max(0.0);
        let contact_work = if b > a {
            (mean_hertz_force(self.hertz_coefficient, a, b) - self.pull_off_magnitude_n) * (b - a)
        } else {
            0.0
        };
        finite(sign * (gap_work + contact_work), "DMT reversible work")
    }

    /// `U(to)-U(from)-F(midpoint)*(to-from)`, evaluated without subtracting
    /// nearly equal energies. The midpoint is supplied by the original state
    /// owner; its floating-point affine-map offset is not assumed to vanish.
    /// This uses the same single-valued potential and accepts its two corners.
    pub fn midpoint_energy_remainder(
        &self,
        from_m: f64,
        to_m: f64,
        midpoint_m: f64,
    ) -> Result<f64, DmtError> {
        self.evaluate(from_m)?;
        self.evaluate(to_m)?;
        self.evaluate(midpoint_m)?;
        if from_m == to_m {
            return Ok(0.0);
        }
        let dx = to_m - from_m;
        let (center, low) = two_sum(0.5 * from_m, 0.5 * to_m);
        let hertz = if from_m > 0.0 && to_m > 0.0 && (2.5 * (0.5 * dx / center)).abs() <= 0.5 {
            // Odd binomial terms beyond the linear term of delta^(5/2).
            // |z| <= .2 makes terms through 33 sufficient for f64 rounding.
            let z = 0.5 * dx / center;
            let mut choose = 1.0;
            let mut power = 1.0;
            let mut sum = 0.0;
            for n in 1..=33 {
                choose *= (3.5 - n as f64) / n as f64;
                power *= z;
                if n >= 3 && n % 2 == 1 {
                    sum += choose * power;
                }
            }
            let remainder = 0.8 * self.hertz_coefficient * center * center.sqrt() * center * sum;
            let force_offset = hertz_force_difference(self.hertz_coefficient, center, midpoint_m)
                + 1.5 * self.hertz_coefficient * center.sqrt() * low;
            remainder + force_offset * dx
        } else {
            let lo = from_m.min(to_m).max(0.0);
            let hi = from_m.max(to_m).max(0.0);
            let work = if hi > lo {
                mean_hertz_force(self.hertz_coefficient, lo, hi)
                    * (hi - lo)
                    * if dx > 0.0 { 1.0 } else { -1.0 }
            } else {
                0.0
            };
            work - hertz_force(self.hertz_coefficient, midpoint_m) * dx
        };
        // The attraction potential is S/2 [delta_+^2-(delta+z0)_+^2].
        // Each fully active quadratic uses its centered remainder, so the
        // large constant attraction cannot hide a tiny Hertz cubic numerator.
        let attraction = self.tail_stiffness_n_per_m
            * (hinge_remainder(from_m, to_m, midpoint_m, 0.0, center, low)
                - hinge_remainder(
                    from_m,
                    to_m,
                    midpoint_m,
                    self.interaction_range_m,
                    center,
                    low,
                ));
        finite(hertz + attraction, "DMT midpoint energy remainder")
    }
}

fn mean_hertz_force(k: f64, a: f64, b: f64) -> f64 {
    // 2 K/5 * (b^(5/2)-a^(5/2))/(b-a), factored in sqrt(a/b).
    let r = (a / b).sqrt();
    let polynomial = 1.0 + r * (1.0 + r * (1.0 + r * (1.0 + r)));
    (k * b * b.sqrt()) * (0.4 * polynomial / (1.0 + r))
}
fn hertz_force(k: f64, x: f64) -> f64 {
    let x = x.max(0.0);
    k * x * x.sqrt()
}
fn hertz_force_difference(k: f64, a: f64, b: f64) -> f64 {
    let a = a.max(0.0);
    let b = b.max(0.0);
    let hi = a.max(b);
    if hi == 0.0 {
        return 0.0;
    }
    let r = (a.min(b) / hi).sqrt();
    k * hi.sqrt() * ((1.0 + r + r * r) / (1.0 + r)) * (a - b)
}
fn hinge_remainder(a: f64, b: f64, m: f64, shift: f64, center: f64, low: f64) -> f64 {
    let dx = b - a;
    if a >= -shift && b >= -shift && m >= -shift {
        // The shift cancels analytically, including affine midpoint roundoff.
        ((center - m) + low) * dx
    } else {
        let ap = (a + shift).max(0.0);
        let bp = (b + shift).max(0.0);
        let mp = (m + shift).max(0.0);
        0.5 * (bp - ap) * (bp + ap) - mp * dx
    }
}
fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let sum = a + b;
    let v = sum - a;
    (sum, (a - (sum - v)) + (b - v))
}
fn finite(value: f64, field: &'static str) -> Result<f64, DmtError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DmtError::Unrepresentable { field })
    }
}
