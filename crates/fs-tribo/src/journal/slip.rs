//! Finite wall slip for the original finite-gap generalized-Newtonian journal.
//!
//! A Navier slip length multiplies the actual cylindrical shear rate,
//! `abs(du_theta/dr - u_theta/r)`, not `abs(du_theta/dr)`. With inner-wall
//! stress `tau_i` and `q=(ri/ro)^2`, the relative cylinder speed is
//! `abs(omega) = integral gamma(tau(r)) dr/r + bi*gamma(tau_i)/ri
//!              + bo*gamma(q*tau_i)/ro`.
//! The original source line supplies both wall viscosities and their stress
//! tangents; no frozen apparent viscosity or radial discretization is needed.
//!
//! Bulk viscous heating and each sliding-interface loss are separate outputs.
//! This owner does not decide which body receives the interface heat, infer a
//! slip length from roughness, or extend the original laminar regime admission.

use super::{GeneralizedNewtonianJournal, JournalShearResponse, invalid};
use crate::TriboError;

/// Prescribed scalar Navier slip lengths on the actual cylindrical walls [m].
///
/// Zero imposes no slip. Lengths remain fixed during an evaluation, while the
/// interfacial mobilities `b / mu(gamma_wall)` follow the actual source law.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct JournalWallSlipLengths {
    pub inner_m: f64,
    pub outer_m: f64,
}

/// Torque, wall slip and the three physically distinct heat contributions.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JournalSlipResponse {
    /// Original signed torque and mechanical dissipation. Its total power
    /// includes the bulk and both interfaces, using the requested wall speed.
    pub shear: JournalShearResponse,
    pub bulk_dissipated_power_w: f64,
    /// Nonnegative heat production at [inner, outer] interfaces [W].
    pub wall_dissipated_power_w: [f64; 2],
    /// Fluid velocity minus wall velocity, at [inner, outer] walls [m/s].
    pub wall_slip_velocity_m_per_s: [f64; 2],
    /// Reconstructed bulk/interface heat minus the conjugate mechanical work.
    /// This retains the bounded scalar solve defect; it is not heat to add.
    pub power_defect_w: f64,
}

/// The original source viscosity line with explicitly supplied wall slip.
#[derive(Clone, Debug)]
pub struct NavierSlipJournal {
    journal: GeneralizedNewtonianJournal,
    slip: JournalWallSlipLengths,
    angular_slip: [f64; 2],
    ratio: f64,
    drop: f64,
}

impl NavierSlipJournal {
    pub fn new(
        journal: GeneralizedNewtonianJournal,
        slip: JournalWallSlipLengths,
    ) -> Result<Self, TriboError> {
        let g = journal.geometry;
        let angular_slip = [slip.inner_m / g.inner_radius_m, slip.outer_m / g.outer_radius_m];
        if [slip.inner_m, slip.outer_m, angular_slip[0], angular_slip[1]]
            .iter()
            .any(|value| !value.is_finite() || *value < 0.0)
        {
            return invalid("journal_wall_slip_length");
        }
        let radius_ratio = g.inner_radius_m / g.outer_radius_m;
        let ratio = radius_ratio.powi(2);
        // Preserve the physical stress width through representable thin gaps.
        let drop = if radius_ratio < 0.5 {
            1.0 - ratio
        } else {
            (g.outer_radius_m - g.inner_radius_m) / g.outer_radius_m * (1.0 + radius_ratio)
        };
        Ok(Self { journal, slip, angular_slip, ratio, drop })
    }

    pub const fn journal(&self) -> &GeneralizedNewtonianJournal { &self.journal }
    pub const fn wall_slip_lengths(&self) -> JournalWallSlipLengths { self.slip }

    /// Evaluate exact source-line torque with curved-wall Navier slip.
    ///
    /// The speed is the signed INNER-minus-OUTER angular velocity. Uniform
    /// rotation of both cylinders contributes no strain or dissipation. The
    /// safeguarded solve allocates nothing and uses analytic stress tangents.
    pub fn response(&self, omega: f64) -> Result<JournalSlipResponse, TriboError> {
        if !omega.is_finite() { return invalid("journal_angular_velocity"); }
        if omega == 0.0 || self.angular_slip == [0.0; 2] {
            let shear = self.journal.response(omega)?;
            return self.finish(shear, omega.abs() + shear.angular_rate_defect_per_s);
        }
        let speed = omega.abs();
        if let Some(mu) = self.journal.constant_viscosity {
            let inner = speed / (0.5 * self.drop + self.angular_slip[0]
                + self.angular_slip[1] * self.ratio);
            if !inner.is_finite()
                || inner > self.journal.knots.last().expect("admitted line").shear_rate_per_s
            {
                return invalid("journal_source_shear_rate_exceeded");
            }
            let rates = [inner, inner * self.ratio];
            let bulk = 0.5 * self.drop * inner;
            let actual = bulk + self.angular_slip[0] * rates[0] + self.angular_slip[1] * rates[1];
            let shear = self.journal.finish(omega, mu * inner, rates, actual - speed, 0, None)?;
            return self.finish(shear, bulk);
        }
        let mut lower = 0.0;
        let mut upper = *self.journal.stresses.last().expect("admitted line");
        let (_, maximum, _, _) = self.angular_speed(upper)?;
        if speed > maximum { return invalid("journal_source_shear_rate_exceeded"); }
        let tolerance = self.journal.solve.absolute_angular_rate_tolerance_per_s
            + self.journal.solve.relative_angular_rate_tolerance * speed;
        let seed_rate = (speed / (0.5 * self.drop + self.angular_slip[0]
            + self.angular_slip[1] * self.ratio))
            .min(self.journal.knots.last().expect("admitted line").shear_rate_per_s);
        let mut stress = seed_rate * self.journal.viscosity_at_shear_rate(seed_rate)?;
        let mut previous_step = upper - lower;
        for iteration in 1..=self.journal.solve.maximum_iterations {
            let (rates, actual, tangent, bulk) = self.angular_speed(stress)?;
            let defect = actual - speed;
            if defect.abs() <= tolerance {
                let shear = self.journal.finish(omega, stress, rates, defect, iteration, None)?;
                return self.finish(shear, bulk);
            }
            if actual < speed { lower = stress; } else { upper = stress; }
            let candidate = stress - defect / tangent;
            let next = if tangent > 0.0 && tangent.is_finite()
                && candidate > lower && candidate < upper
                && (candidate - stress).abs() <= 0.5 * previous_step
            {
                candidate
            } else {
                lower + 0.5 * (upper - lower)
            };
            previous_step = (next - stress).abs();
            stress = next;
        }
        invalid("journal_angular_rate_solve_exhausted")
    }

    fn angular_speed(&self, stress: f64) -> Result<([f64; 2], f64, f64, f64), TriboError> {
        let (rates, bulk, mut tangent) = self.journal.angular_speed(stress, self.ratio, self.drop)?;
        let speed = bulk + self.angular_slip[0] * rates[0] + self.angular_slip[1] * rates[1];
        for (wall, ratio) in [(0, 1.0), (1, self.ratio)] {
            let rate = rates[wall];
            // At a source knot take the preceding segment's one-sided tangent;
            // bracketing remains authoritative across constitutive kinks.
            let hi = self.journal.knots.partition_point(|p| p.shear_rate_per_s < rate).max(1);
            let mu = self.journal.viscosity_at_shear_rate(rate)?;
            let stress_tangent = mu + rate * self.journal.slopes[hi - 1];
            tangent += self.angular_slip[wall] * ratio / stress_tangent;
        }
        if !speed.is_finite() { return invalid("journal_angular_rate_overflow"); }
        Ok((rates, speed, tangent, bulk))
    }

    fn finish(
        &self,
        shear: JournalShearResponse,
        bulk_speed: f64,
    ) -> Result<JournalSlipResponse, TriboError> {
        let torque = shear.resisting_torque_nm.abs();
        let rates = shear.shear_rate_per_s;
        let bulk = torque * bulk_speed;
        let wall = [
            torque * (self.angular_slip[0] * rates[0]),
            torque * (self.angular_slip[1] * rates[1]),
        ];
        let sign = -shear.resisting_torque_nm.signum();
        let slip = [
            -sign * (self.slip.inner_m * rates[0]),
            sign * (self.slip.outer_m * rates[1]),
        ];
        let defect = bulk + wall[0] + wall[1] - shear.dissipated_power_w;
        if [bulk, wall[0], wall[1], slip[0], slip[1], defect].iter().any(|v| !v.is_finite()) {
            return invalid("journal_wall_slip_power_unrepresentable");
        }
        Ok(JournalSlipResponse {
            shear,
            bulk_dissipated_power_w: bulk,
            wall_dissipated_power_w: wall,
            wall_slip_velocity_m_per_s: slip,
            power_defect_w: defect,
        })
    }
}
