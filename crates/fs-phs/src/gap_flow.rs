//! Rung-A gap-flow law for a flow-controlling deformable aperture.
//!
//! IRP M3.1-W3 (DEF-13). The legacy [`crate::bernoulli_volume_flow`] is the
//! inviscid limit `U = w h sgn(dp) sqrt(2|dp|/rho)`: it is singular at
//! closure, has no laminar channel term, and gives an infinite
//! `dU/d(dp)` at zero drop. This owner keeps that law as its declared
//! inviscid rung and adds the two terms a real slit has:
//!
//! ```text
//! dp = rho U |U| / (2 (C_c w h)^2) + 12 mu L_c U / (w h^3)
//! ```
//!
//! with `C_c` the contraction coefficient, `L_c` the channel length and
//! `mu` the dynamic viscosity. The quasi-steady solve is analytic and
//! regular at `dp = 0` and at `h = 0` whenever the laminar term is
//! present. The unsteady term `d/dt(I(h) U)` is exposed through
//! [`GapFlowLaw::inertance`]; a caller that keeps an inertance state adds
//! `I dU/dt + U dI/dt` as a pressure source, while a caller that does not
//! declares the quasi-steady rung.
//!
//! The law is the same port for a single reed, a paired reed, a modal
//! aperture, a brass lip and a vocal fold. It is not instrument-shaped.

/// The rung-A quasi-steady gap-flow law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GapFlowLaw {
    /// Contraction coefficient `C_c` in `(0, 1]`.
    pub contraction: f64,
    /// Laminar channel length `L_c` [m], `>= 0`.
    pub channel_length_m: f64,
    /// Unsteady effective length `L_eff` [m], `>= 0`.
    pub effective_length_m: f64,
    /// Dynamic viscosity `mu` [Pa s], `>= 0`.
    pub dynamic_viscosity_pa_s: f64,
}

/// Admission refusal for a gap-flow evaluation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GapFlowError {
    /// A width, opening, density or pressure is not finite, or a width or
    /// density is not positive.
    InvalidInput,
    /// A coefficient is outside its declared domain.
    InvalidLaw,
}

impl core::fmt::Display for GapFlowError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidInput => f.write_str("gap flow: non-finite or non-positive width, opening, density or pressure"),
            Self::InvalidLaw => f.write_str("gap flow: contraction, channel length, effective length or viscosity outside its domain"),
        }
    }
}

impl std::error::Error for GapFlowError {}

impl Default for GapFlowLaw {
    /// The inviscid, zero-length limit: bit-for-bit
    /// [`crate::bernoulli_volume_flow`] wherever that law is defined.
    fn default() -> Self {
        Self {
            contraction: 1.0,
            channel_length_m: 0.0,
            effective_length_m: 0.0,
            dynamic_viscosity_pa_s: 0.0,
        }
    }
}

impl GapFlowLaw {
    /// Admit the law's own coefficients.
    ///
    /// # Errors
    /// [`GapFlowError::InvalidLaw`] when a coefficient is non-finite or
    /// outside its domain.
    pub fn admit(&self) -> Result<(), GapFlowError> {
        let valid = self.contraction.is_finite()
            && self.contraction > 0.0
            && self.contraction <= 1.0
            && self.channel_length_m.is_finite()
            && self.channel_length_m >= 0.0
            && self.effective_length_m.is_finite()
            && self.effective_length_m >= 0.0
            && self.dynamic_viscosity_pa_s.is_finite()
            && self.dynamic_viscosity_pa_s >= 0.0;
        if valid {
            Ok(())
        } else {
            Err(GapFlowError::InvalidLaw)
        }
    }

    /// The unsteady inertance `I(h) = rho L_eff / (w h)` [kg/m^2].
    ///
    /// Returns `None` at `h = 0` (the slit is closed and carries no
    /// unsteady mass) and for invalid inputs.
    #[must_use]
    pub fn inertance(&self, width_m: f64, opening_m: f64, density: f64) -> Option<f64> {
        if width_m > 0.0
            && opening_m > 0.0
            && density > 0.0
            && width_m.is_finite()
            && opening_m.is_finite()
            && density.is_finite()
        {
            Some(density * self.effective_length_m / (width_m * opening_m))
        } else {
            None
        }
    }

    /// The quasi-steady pressure drop for a signed volume flow.
    ///
    /// `dp = A U|U| + B U`, with `A = rho/(2(C_c w h)^2)` and
    /// `B = 12 mu L_c/(w h^3)`. It is the inverse of
    /// [`Self::volume_flow`] and returns the same sign as `flow`.
    #[must_use]
    pub fn pressure_drop(&self, width_m: f64, opening_m: f64, flow_m3_s: f64, density: f64) -> f64 {
        if !(width_m > 0.0) || !(opening_m > 0.0) || !(density > 0.0) || !flow_m3_s.is_finite() {
            return 0.0;
        }
        let a = density / (2.0 * (self.contraction * width_m * opening_m).powi(2));
        let b = 12.0 * self.dynamic_viscosity_pa_s * self.channel_length_m
            / (width_m * opening_m.powi(3));
        a * flow_m3_s * flow_m3_s.abs() + b * flow_m3_s
    }

    /// The unsteady rung-A pressure drop,
    /// `dp = A U|U| + B U + I(h) dU/dt + U dI/dt`.
    ///
    /// `I dU/dt + U dI/dt` is the exact derivative of `I(h) U` for a
    /// constant effective length; `dI/dt = -I (dh/dt)/h`. A closed slit
    /// (`h = 0`) has no unsteady mass, so the quasi-steady value is
    /// returned.
    #[must_use]
    pub fn pressure_drop_unsteady(
        &self,
        width_m: f64,
        opening_m: f64,
        flow_m3_s: f64,
        flow_rate_m3_s2: f64,
        opening_rate_m_s: f64,
        density: f64,
    ) -> f64 {
        let quasi = self.pressure_drop(width_m, opening_m, flow_m3_s, density);
        let Some(inertance) = self.inertance(width_m, opening_m, density) else {
            return quasi;
        };
        if !flow_rate_m3_s2.is_finite() || !opening_rate_m_s.is_finite() {
            return quasi;
        }
        quasi + inertance * flow_rate_m3_s2 - inertance * opening_rate_m_s * flow_m3_s / opening_m
    }

    /// The signed quasi-steady volume flow through a slit.
    ///
    /// The inviscid limit is evaluated with the legacy owner's exact
    /// arithmetic, so `GapFlowLaw::default()` and
    /// [`crate::bernoulli_volume_flow`] agree bit for bit. With the
    /// laminar term present the solve is
    /// `U = (-B + sqrt(B^2 + 4 A |dp|)) / (2 A)`.
    #[must_use]
    pub fn volume_flow(&self, width_m: f64, opening_m: f64, dp_pa: f64, density: f64) -> f64 {
        match self.volume_flow_and_jacobian(width_m, opening_m, dp_pa, density) {
            Ok((flow, _, _)) => flow,
            Err(_) => 0.0,
        }
    }

    /// Volume flow and its analytic partials `(U, dU/d(dp), dU/dh)`.
    ///
    /// The derivatives are those of the quasi-steady law; at `h = 0` the
    /// aperture is closed and both are zero. The pure-Bernoulli limit is
    /// intentionally singular at `dp = 0` (`dU/d(dp) = +inf`), which is
    /// exactly the defect DEF-13 names; a law with a positive laminar
    /// coefficient is regular.
    ///
    /// # Errors
    /// [`GapFlowError::InvalidInput`] for a non-finite or non-positive
    /// input, or [`GapFlowError::InvalidLaw`] for an invalid coefficient.
    pub fn volume_flow_and_jacobian(
        &self,
        width_m: f64,
        opening_m: f64,
        dp_pa: f64,
        density: f64,
    ) -> Result<(f64, f64, f64), GapFlowError> {
        self.admit()?;
        let valid = width_m.is_finite()
            && width_m > 0.0
            && opening_m.is_finite()
            && opening_m >= 0.0
            && dp_pa.is_finite()
            && density.is_finite()
            && density > 0.0;
        if !valid {
            return Err(GapFlowError::InvalidInput);
        }
        if opening_m == 0.0 {
            return Ok((0.0, 0.0, 0.0));
        }
        let b = 12.0 * self.dynamic_viscosity_pa_s * self.channel_length_m
            / (width_m * opening_m.powi(3));
        if b == 0.0 && self.contraction == 1.0 {
            // The declared inviscid rung: keep the legacy owner's arithmetic
            // and its exact dead-zone behavior, bit for bit.
            let flow = crate::bernoulli_volume_flow(width_m, opening_m, dp_pa, density);
            let derivative = if dp_pa.abs() < 1.0e-12 {
                f64::INFINITY
            } else {
                width_m * opening_m * (2.0 / (density * dp_pa.abs())).sqrt() * 0.5
            };
            return Ok((flow, derivative, flow / opening_m));
        }
        let a = density / (2.0 * (self.contraction * width_m * opening_m).powi(2));
        let q = dp_pa.abs();
        let disc = (b * b + 4.0 * a * q).sqrt();
        // Rationalized root: U = 2q/(sqrt(B^2+4Aq)+B), which avoids the
        // cancellation of (disc - B) when the laminar term dominates.
        let denom = disc + b;
        let u_mag = 2.0 * q / denom;
        let flow = if dp_pa < 0.0 { -u_mag } else { u_mag };
        // dU/d|dp| = 1 / sqrt(B^2 + 4 A q)
        let d_dp = 1.0 / disc;
        // A = a0 h^-2, B = b0 h^-3, so dA/dh = -2A/h and dB/dh = -3B/h.
        let d_disc = (b * (-3.0 * b / opening_m) + 2.0 * q * (-2.0 * a / opening_m)) / disc;
        let d_denom = d_disc - 3.0 * b / opening_m;
        let d_h = -2.0 * q * d_denom / (denom * denom);
        let d_h = if dp_pa < 0.0 { -d_h } else { d_h };
        Ok((flow, d_dp, d_h))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn law() -> GapFlowLaw {
        GapFlowLaw {
            contraction: 0.61,
            channel_length_m: 0.002,
            effective_length_m: 0.003,
            dynamic_viscosity_pa_s: 1.81e-5,
        }
    }

    #[test]
    fn bernoulli_limit_is_bitwise() {
        let law = GapFlowLaw::default();
        for opening in [1e-6, 1e-4, 4e-4, 1e-3] {
            for dp in [-1e3, -1e-3, -1e-12, 0.0, 1e-12, 1e-3, 1e3] {
                assert_eq!(
                    law.volume_flow(0.013, opening, dp, 1.2).to_bits(),
                    crate::bernoulli_volume_flow(0.013, opening, dp, 1.2).to_bits(),
                    "opening {opening}, dp {dp}"
                );
            }
        }
    }

    #[test]
    fn no_dead_zone_and_regular_jacobian() {
        let law = law();
        for opening in [1e-9, 1e-6, 4e-4] {
            let (u0, d0, h0) = law
                .volume_flow_and_jacobian(0.013, opening, 0.0, 1.2)
                .unwrap();
            assert_eq!(u0, 0.0);
            assert!(d0.is_finite() && d0 > 0.0, "dU/ddp {d0} at h {opening}");
            assert!(h0.is_finite());
            for dp in [-1e3, -1e-3, -1e-12, 1e-12, 1e-3, 1e3] {
                let (u, du, dh) = law
                    .volume_flow_and_jacobian(0.013, opening, dp, 1.2)
                    .unwrap();
                assert!(u.is_finite() && du.is_finite() && dh.is_finite());
                assert_eq!(u.signum(), dp.signum(), "the flow follows the drop");
            }
        }
        let (u, _, _) = law
            .volume_flow_and_jacobian(0.013, 0.0, 100.0, 1.2)
            .unwrap();
        assert_eq!(u, 0.0, "a closed slit carries nothing");
    }

    #[test]
    fn poiseuille_limit_is_exact() {
        let law = GapFlowLaw {
            contraction: 1.0,
            channel_length_m: 0.002,
            effective_length_m: 0.0,
            dynamic_viscosity_pa_s: 1.81e-5,
        };
        let (w, h, dp): (f64, f64, f64) = (0.013, 4e-4, 1e-6);
        let expected =
            w * h.powi(3) * dp / (12.0 * law.dynamic_viscosity_pa_s * law.channel_length_m);
        let actual = law.volume_flow(w, h, dp, 1.2);
        assert!(
            (actual - expected).abs() <= 1e-12 * expected.abs().max(1.0),
            "{actual} vs {expected}"
        );
    }

    #[test]
    fn jacobian_matches_central_differences() {
        let law = law();
        let (w, h, rho) = (0.013, 4e-4, 1.2);
        for dp in [-1e3, -1e-3, 0.0, 1e-3, 1e3] {
            let (_, du, dh) = law.volume_flow_and_jacobian(w, h, dp, rho).unwrap();
            // A relative step keeps the central difference above the
            // subtraction's round-off while staying inside the local
            // quadratic regime.
            let e_dp = 1e-6 * dp.abs().max(1.0);
            let e_h = 1e-6 * h;
            let fd_dp = (law.volume_flow(w, h, dp + e_dp, rho)
                - law.volume_flow(w, h, dp - e_dp, rho))
                / (2.0 * e_dp);
            let fd_h = (law.volume_flow(w, h + e_h, dp, rho)
                - law.volume_flow(w, h - e_h, dp, rho))
                / (2.0 * e_h);
            assert!(
                (du - fd_dp).abs() <= 1e-5 * fd_dp.abs().max(1e-12),
                "dU/ddp {du} vs {fd_dp}"
            );
            assert!(
                (dh - fd_h).abs() <= 1e-5 * fd_h.abs().max(1e-12),
                "dU/dh {dh} vs {fd_h}"
            );
        }
    }

    #[test]
    fn pressure_drop_inverts_volume_flow() {
        let law = law();
        for dp in [-1e3, -1e-6, 0.0, 1e-6, 1e3] {
            let u = law.volume_flow(0.013, 4e-4, dp, 1.2);
            let back = law.pressure_drop(0.013, 4e-4, u, 1.2);
            assert!(
                (back - dp).abs() <= 1e-10 * dp.abs().max(1.0),
                "{back} vs {dp}"
            );
        }
    }

    #[test]
    fn inertance_is_rho_l_over_wh() {
        let law = law();
        let got = law.inertance(0.013, 4e-4, 1.2).unwrap();
        assert!((got - 1.2 * 0.003 / (0.013 * 4e-4)).abs() <= 1e-15);
        assert_eq!(law.inertance(0.013, 0.0, 1.2), None);
    }

    #[test]
    fn unsteady_pressure_is_the_exact_d_dt_of_inertance_times_flow() {
        let law = law();
        let (w, h, rho) = (0.013, 4e-4, 1.2);
        let u = 2.3e-5;
        let du = 1.7e-2;
        let dh = 3.1e-3;
        let inertance = law.inertance(w, h, rho).unwrap();
        let expected = law.pressure_drop(w, h, u, rho) + inertance * du - inertance * dh * u / h;
        let actual = law.pressure_drop_unsteady(w, h, u, du, dh, rho);
        assert!(
            (actual - expected).abs() <= 1e-12 * expected.abs().max(1.0),
            "{actual} vs {expected}"
        );
        // A closed slit has no unsteady mass; only the (zero) quasi-steady term remains.
        assert_eq!(law.pressure_drop_unsteady(w, 0.0, u, du, dh, rho), 0.0);
    }
}
