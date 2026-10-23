//! A declared compliant drive: a massless spring-dashpot in series between a
//! friction port's prescribed surface motion and the contact.
//!
//! Pitteroff & Woodhouse 1998 ("Finite width bowing, Part I", Acta Acustica
//! 84:543-562) model bow hair as a massless longitudinal spring-dashpot (their
//! eq. 1, used as eq. 8): the longitudinal force on the string equals the hair
//! force `f = -(k eta + c d/dt eta)`, with `eta` the hair stretch relative to
//! the rigid bow stick. The hair has no mass, so the force at the stick end
//! equals the force at the contact.
//!
//! The spring force is taken at the midpoint stretch (the trapezoidal rule)
//! and the dashpot at the backward-Euler velocity, so over a step the force
//! `f = k (eta + eta0)/2 + c (eta - eta0)/dt` does work `f (eta - eta0) =
//! 1/2 k (eta^2 - eta0^2) + c (eta - eta0)^2/dt` EXACTLY: the elastic energy
//! `1/2 k eta^2` and the dashpot loss `c (d eta/dt)^2 dt` close the discrete
//! ledger with no spurious numerical energy (a plain backward-Euler spring
//! would leave a `1/2 k (eta - eta0)^2` residual). The element then adds a
//! series compliance `1/(k/2 + c/dt)` to the port's load line and shifts the
//! free slip by the relaxation of the stored stretch, so the string no longer
//! follows the stick rigidly. For a finite-width contact each cell carries its
//! own such element, so the contact-point velocities differ across the width
//! and differential slip EMERGES (P&W Fig. 9: the rigid bow's singular edge
//! forces spread into a smooth interior distribution). Both energies enter the
//! joint ledger, so the per-step balance still closes.
//!
//! This is the generic "declared compliant drive" of plan item M4.4. When the
//! driven-body interface (M4.4-W1) lands, this series element re-homes onto it;
//! until then it is attached directly to the existing friction ports.

use super::*;

/// One port's series spring-dashpot drive and its committed stretch state.
/// The stretch is the displacement of the rigid stick ahead of the contact
/// point; it is carried across steps as the hair's memory.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SeriesCompliantDrive {
    stiffness_n_m: f64,
    damping_kg_s: f64,
    eta_m: f64,
}

impl SeriesCompliantDrive {
    /// A drive with positive finite stiffness and nonnegative finite damping,
    /// starting from a relaxed (zero) stretch.
    ///
    /// # Errors
    /// Non-physical stiffness or damping.
    pub fn new(stiffness_n_m: f64, damping_kg_s: f64) -> Result<Self, ModalFrictionError> {
        if !stiffness_n_m.is_finite()
            || stiffness_n_m <= 0.0
            || !damping_kg_s.is_finite()
            || damping_kg_s < 0.0
        {
            return Err(invalid(
                "a series compliant drive needs positive finite stiffness and nonnegative finite damping",
            )
            .into());
        }
        Ok(Self {
            stiffness_n_m,
            damping_kg_s,
            eta_m: 0.0,
        })
    }

    /// Series spring constant `k` [N/m].
    #[must_use]
    pub const fn stiffness_n_m(&self) -> f64 {
        self.stiffness_n_m
    }
    /// Series dashpot constant `c` [kg/s].
    #[must_use]
    pub const fn damping_kg_s(&self) -> f64 {
        self.damping_kg_s
    }
    /// Committed hair stretch `eta0` [m].
    #[must_use]
    pub const fn stretch_m(&self) -> f64 {
        self.eta_m
    }

    /// Effective one-step series denominator `k/2 + c/dt` [N/m] (midpoint
    /// spring, backward-Euler dashpot).
    #[must_use]
    fn denominator(&self, dt: f64) -> f64 {
        0.5 * self.stiffness_n_m + self.damping_kg_s / dt
    }

    /// Series compliance `1/(k/2 + c/dt)` [m/N] added to the port's self
    /// compliance before the load line is built.
    #[must_use]
    pub fn series_compliance_m_n(&self, dt: f64) -> f64 {
        1.0 / self.denominator(dt)
    }

    /// Free-slip shift `+k eta0/(k/2 + c/dt)` [m] from the stored stretch: the
    /// drive the stored stretch supplies before any contact force.
    #[must_use]
    pub fn free_slip_shift_m(&self, dt: f64) -> f64 {
        self.stiffness_n_m * self.eta_m / self.denominator(dt)
    }

    /// Stretch after a step given the solved contact force `f`:
    /// `eta = eta0 + (f - k eta0)/(k/2 + c/dt)`, the midpoint-rule update.
    #[must_use]
    pub fn staged_stretch_m(&self, contact_force_n: f64, dt: f64) -> f64 {
        self.eta_m + (contact_force_n - self.stiffness_n_m * self.eta_m) / self.denominator(dt)
    }

    /// The reduction of the friction-interface slip relative to a rigid drive,
    /// `delta eta = eta - eta0` [m], for the solved contact force `f`. The true
    /// interface slip is `slip_rigid - stretch_increment`.
    #[must_use]
    pub fn stretch_increment_m(&self, contact_force_n: f64, dt: f64) -> f64 {
        self.staged_stretch_m(contact_force_n, dt) - self.eta_m
    }

    /// Stored elastic energy `1/2 k eta^2` [J] at a given stretch.
    #[must_use]
    pub fn spring_energy_j(&self, stretch_m: f64) -> f64 {
        0.5 * self.stiffness_n_m * stretch_m * stretch_m
    }

    /// Dashpot dissipation `c (eta - eta0)^2 / dt` [J] over one step reaching
    /// `stretch_m`. Nonnegative.
    #[must_use]
    pub fn dashpot_loss_j(&self, stretch_m: f64, dt: f64) -> f64 {
        let increment = stretch_m - self.eta_m;
        self.damping_kg_s * increment * increment / dt
    }

    /// Commit the staged stretch as the drive's memory for the next step.
    pub fn commit(&mut self, stretch_m: f64) {
        self.eta_m = stretch_m;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stiff_drive_approaches_the_rigid_limit() {
        // As k -> infinity the series compliance and the stored-stretch effects
        // vanish, so a stiff hair recovers the rigid drive.
        let dt = 1.0 / 48_000.0;
        let drive = SeriesCompliantDrive::new(1.0e12, 10.0).unwrap();
        assert!(drive.series_compliance_m_n(dt) < 1e-11);
        assert!(drive.free_slip_shift_m(dt).abs() <= 0.0); // eta0 = 0
        // A 1 N contact force barely stretches a stiff hair.
        assert!(drive.staged_stretch_m(1.0, dt).abs() < 1e-11);
    }

    #[test]
    fn soft_drive_has_large_compliance_and_stretches() {
        let dt = 1.0 / 48_000.0;
        let drive = SeriesCompliantDrive::new(1.0e4, 0.0).unwrap();
        // Midpoint spring, undamped: compliance 1/(k/2) and eta = 2 f/k.
        assert!((drive.series_compliance_m_n(dt) - 2.0e-4).abs() <= 1e-12);
        assert!((drive.staged_stretch_m(2.0, dt) - 4.0e-4).abs() <= 1e-12);
        assert!((drive.spring_energy_j(4.0e-4) - 0.5 * 1.0e4 * 1.6e-7).abs() <= 1e-18);
        assert!(drive.dashpot_loss_j(4.0e-4, dt).abs() <= 0.0);
    }

    #[test]
    fn stretch_increment_matches_series_relation() {
        // delta eta = (f - k eta0)/(k/2 + c/dt): the slip reduction identity.
        let dt = 1.0 / 96_000.0;
        let mut drive = SeriesCompliantDrive::new(5.0e4, 8.0).unwrap();
        drive.commit(3.0e-5);
        let f = 1.7;
        let expected = (f - 5.0e4 * 3.0e-5) / (0.5 * 5.0e4 + 8.0 / dt);
        assert!((drive.stretch_increment_m(f, dt) - expected).abs() <= 1e-18);
    }

    #[test]
    fn midpoint_rule_closes_the_discrete_energy_identity() {
        // The whole reason for the midpoint spring: the work the end-of-step
        // force does over the step equals the stored-energy change plus the
        // dashpot loss EXACTLY, f (eta - eta0) = 1/2 k (eta^2 - eta0^2) + loss,
        // for any stretch history and force. A backward-Euler spring would miss
        // this by 1/2 k (eta - eta0)^2.
        let dt = 1.0 / 48_000.0;
        for (k, c, eta0, f) in [
            (6.0e4, 10.0, 0.0, 0.7),
            (6.0e4, 10.0, 2.0e-5, -1.3),
            (1.0e4, 0.0, -4.0e-5, 0.9),
            (2.0e4, 8.0, 5.0e-5, 0.0),
        ] {
            let mut drive = SeriesCompliantDrive::new(k, c).unwrap();
            drive.commit(eta0);
            let eta = drive.staged_stretch_m(f, dt);
            let work = f * (eta - eta0);
            let stored = drive.spring_energy_j(eta) - drive.spring_energy_j(eta0);
            let loss = drive.dashpot_loss_j(eta, dt);
            // The identity is exact in real arithmetic; the only slack is the
            // rounding of the separately summed stored and loss terms.
            let scale = work.abs().max(stored.abs()).max(loss.abs()).max(1e-15);
            assert!(
                (work - stored - loss).abs() <= 1e-12 * scale,
                "k={k} c={c}: work {work} != stored {stored} + loss {loss}"
            );
        }
    }

    #[test]
    fn dashpot_loss_is_nonnegative_and_damping_free_stores_only() {
        let dt = 1.0 / 48_000.0;
        let damped = SeriesCompliantDrive::new(5.0e4, 12.0).unwrap();
        for f in [-3.0, -0.2, 0.0, 0.5, 4.0] {
            let eta = damped.staged_stretch_m(f, dt);
            assert!(damped.dashpot_loss_j(eta, dt) >= 0.0);
        }
    }

    #[test]
    fn invalid_parameters_refuse() {
        assert!(SeriesCompliantDrive::new(0.0, 1.0).is_err());
        assert!(SeriesCompliantDrive::new(-1.0, 1.0).is_err());
        assert!(SeriesCompliantDrive::new(1.0e4, -1.0).is_err());
        assert!(SeriesCompliantDrive::new(f64::NAN, 1.0).is_err());
        assert!(SeriesCompliantDrive::new(1.0e4, f64::INFINITY).is_err());
    }
}
