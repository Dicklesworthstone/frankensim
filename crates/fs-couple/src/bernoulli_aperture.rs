//! Coupling-layer view of the fs-phs Bernoulli valve.
//!
//! The constitutive law lives in [`fs_phs::bernoulli_volume_flow`] and
//! [`fs_phs::quasistatic_aperture_opening`]. This type is only the
//! geometry + closing-pressure data a coupling step carries.

/// Quasistatic aperture: `y = H max(0, 1 − Δp/P_c)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BernoulliAperture {
    /// Rest opening `H` [m].
    pub rest_opening_m: f64,
    /// Slit width [m].
    pub width_m: f64,
    /// Pressure drop that just closes the slit [Pa].
    pub closing_pressure_pa: f64,
}

impl BernoulliAperture {
    /// Opening height at pressure drop `dp = p_upstream − p_downstream`.
    #[must_use]
    pub fn opening_m(self, dp: f64) -> f64 {
        fs_phs::quasistatic_aperture_opening(self.rest_opening_m, self.closing_pressure_pa, dp)
    }

    /// Volume flow [m³/s] through the slit, `U = w y sgn(Δp) √(2|Δp|/ρ)`.
    #[must_use]
    pub fn volume_flow(self, dp: f64, density: f64) -> f64 {
        fs_phs::bernoulli_volume_flow(self.width_m, self.opening_m(dp), dp, density)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn closes_at_the_named_pressure_and_flows_as_sqrt_dp() {
        let a = BernoulliAperture {
            rest_opening_m: 4.0e-4,
            width_m: 0.01,
            closing_pressure_pa: 1_000.0,
        };
        assert_eq!(a.opening_m(1_000.0), 0.0);
        assert_eq!(a.volume_flow(1_200.0, 1.2), 0.0);
        assert!((a.opening_m(0.0) - 4.0e-4).abs() < 1.0e-16);
        let u1 = a.volume_flow(100.0, 1.2);
        let u4 = a.volume_flow(400.0, 1.2);
        // y halves when Δp goes 100 → 400, √Δp doubles, so U is unchanged
        // at those two points? y(100)=0.9 H, y(400)=0.6 H, √4=2,
        // U4/U1 = (0.6/0.9)*2 = 1.333.
        assert!((u4 / u1 - 4.0_f64.sqrt() * 0.6 / 0.9).abs() < 1.0e-12);
        assert!(u1 > 0.0);
    }
}
