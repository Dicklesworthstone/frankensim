//! Ordinary-water/liquid-vapor surface tension from IAPWS R1-76(2014).
//!
//! The original ITS-90 correlation is sigma = B*tau^mu*(1+b*tau), with
//! tau = 1-T/Tc, B = 235.8 mN/m, b = -0.625 and mu = 1.256. It describes
//! ordinary liquid water in equilibrium with its PURE vapor. This evaluator
//! admits the inclusive triple-to-critical interval, 273.16--647.096 K;
//! the release's optional supercooled extrapolation is not exposed.
//!
//! Values and analytic temperature derivatives have coherent SI units. The
//! derivative follows the saturation curve, not an arbitrary pressure path.
//! At the critical endpoint both the value and its limiting first derivative
//! are zero; this does not create a distinguishable liquid/vapor interface.
//! No other gas pairing, solute/surfactant, wetting, curvature correction,
//! surface entropy or heat supply is inferred. The correlation is separate
//! from both saturation caloric sources and cannot change their phase datums.
//!
//! Attribution: International Association for the Properties of Water and Steam,
//! Revised Release on Surface Tension of Ordinary Water Substance (June 2014),
//! equation and Table 1: <https://iapws.org/technical-guidance/release/Surf-H2O.download>.
//! Source uncertainty is distinct from numerical evaluation error. Evaluation
//! is allocation-free and uses the repository's deterministic elementary math.

use core::fmt;
use fs_math::det;
use fs_qty::{Qty, SurfaceTension, Temperature};

/// Original source release and temperature scale.
pub const SOURCE_RELEASE: &str = "IAPWS R1-76(2014), ITS-90";
/// Official equation, coefficients and rounded reference table.
pub const SOURCE_URL: &str = "https://iapws.org/technical-guidance/release/Surf-H2O.download";
/// Lower inclusive source temperature [K].
pub const TRIPLE_POINT_TEMPERATURE_K: f64 = 273.16;
/// Upper inclusive source temperature [K]; no phase split exists here.
pub const CRITICAL_TEMPERATURE_K: f64 = 647.096;

/// Surface-tension derivative along the saturation curve [N/(m K)].
pub type SurfaceTensionTemperatureDerivative = Qty<0, 1, -2, -1, 0, 0>;

/// Source-domain or finite-representation refusal; no clipping/extrapolation.
#[derive(Clone, Debug, PartialEq)]
pub enum WaterSurfaceTensionError {
    TemperatureOutOfRange { temperature_k: f64 },
    Unrepresentable,
}
impl fmt::Display for WaterSurfaceTensionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "IAPWS water/liquid-vapor surface tension: {self:?}")
    }
}
impl std::error::Error for WaterSurfaceTensionError {}

/// Immutable original-source coefficient and analytic tangent at one temperature.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterSurfaceTensionState {
    temperature_k: f64,
    surface_tension_n_m: f64,
    derivative_n_m_k: f64,
}
impl WaterSurfaceTensionState {
    /// Evaluate the original R1-76(2014) liquid/pure-vapor correlation.
    ///
    /// One fractional power supplies both sigma and d(sigma)/dT. Forming tau
    /// from Tc-T preserves the distance of adjacent temperatures from criticality.
    pub fn at_temperature(temperature: Temperature) -> Result<Self, WaterSurfaceTensionError> {
        let t = temperature.value();
        if !t.is_finite() || !(TRIPLE_POINT_TEMPERATURE_K..=CRITICAL_TEMPERATURE_K).contains(&t) {
            return Err(WaterSurfaceTensionError::TemperatureOutOfRange { temperature_k: t });
        }
        const B_N_M: f64 = 0.2358;
        const B: f64 = -0.625;
        const MU: f64 = 1.256;
        let tau = (CRITICAL_TEMPERATURE_K - t) / CRITICAL_TEMPERATURE_K;
        let power = det::pow(tau, MU - 1.0);
        let sigma = B_N_M * tau * power * (1.0 + B * tau);
        let derivative = if tau == 0.0 {
            0.0
        } else {
            -(B_N_M / CRITICAL_TEMPERATURE_K) * power * (MU + B * (MU + 1.0) * tau)
        };
        if !sigma.is_finite() || sigma < 0.0 || !derivative.is_finite() || derivative > 0.0 {
            return Err(WaterSurfaceTensionError::Unrepresentable);
        }
        Ok(Self {
            temperature_k: t,
            surface_tension_n_m: sigma,
            derivative_n_m_k: derivative,
        })
    }

    pub const fn temperature(self) -> Temperature {
        Temperature::new(self.temperature_k)
    }
    /// Original water/pure-vapor coefficient [N/m], zero only at criticality.
    pub const fn surface_tension(self) -> SurfaceTension {
        SurfaceTension::new(self.surface_tension_n_m)
    }
    /// Analytic saturation-path derivative [N/(m K)], with its critical limit.
    pub const fn temperature_derivative(self) -> SurfaceTensionTemperatureDerivative {
        SurfaceTensionTemperatureDerivative::new(self.derivative_n_m_k)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::water::{SaturatedWaterState, if97::If97SaturationPair};

    fn state(t: f64) -> WaterSurfaceTensionState {
        WaterSurfaceTensionState::at_temperature(Temperature::new(t)).unwrap()
    }

    #[test]
    fn g1_original_iapws_table_values_use_newtons_per_metre() {
        // Literal Table 1 columns (1),(4), degrees C and mN/m. The tolerance
        // is half the displayed 0.01 mN/m precision, not material uncertainty.
        let table = [
            (0.01, 75.65),
            (5.0, 74.94),
            (20.0, 72.74),
            (25.0, 71.97),
            (50.0, 67.94),
            (100.0, 58.91),
            (150.0, 48.74),
            (200.0, 37.67),
            (250.0, 26.04),
            (300.0, 14.36),
            (350.0, 3.67),
            (360.0, 1.88),
            (370.0, 0.39),
        ];
        for (celsius, millinewtons_per_metre) in table {
            // Decimal 0.01 + 273.15 rounds below the exact admitted binary64
            // triple-point constant; select the source endpoint explicitly.
            let t = if celsius == 0.01 {
                TRIPLE_POINT_TEMPERATURE_K
            } else {
                celsius + 273.15
            };
            let surface = state(t);
            let coefficient: SurfaceTension = surface.surface_tension();
            assert!(
                (coefficient.value() - 1e-3 * millinewtons_per_metre).abs() <= 5.0001e-6,
                "Table 1 at {celsius} C: {} N/m",
                coefficient.value(),
            );
        }
    }

    #[test]
    fn g1_analytic_surface_tangent_matches_independent_temperature_differences() {
        for t in [274.0, 298.15, 373.15, 473.15, 573.15, 623.15, 646.0] {
            let h = 0.001;
            let fd = (state(t - 2.0 * h).surface_tension().value()
                - 8.0 * state(t - h).surface_tension().value()
                + 8.0 * state(t + h).surface_tension().value()
                - state(t + 2.0 * h).surface_tension().value())
                / (12.0 * h);
            let slope: SurfaceTensionTemperatureDerivative = state(t).temperature_derivative();
            assert!(slope.value() < 0.0);
            assert!((slope.value() - fd).abs() <= 1e-8 * slope.value().abs());
        }
    }

    #[test]
    fn g0_surface_domain_and_critical_limit_preserve_phase_coalescence() {
        for t in [
            f64::NAN,
            f64::INFINITY,
            f64::NEG_INFINITY,
            TRIPLE_POINT_TEMPERATURE_K.next_down(),
            CRITICAL_TEMPERATURE_K.next_up(),
        ] {
            assert!(matches!(
                WaterSurfaceTensionState::at_temperature(Temperature::new(t)),
                Err(WaterSurfaceTensionError::TemperatureOutOfRange { .. })
            ));
        }
        let critical = state(CRITICAL_TEMPERATURE_K);
        assert_eq!(critical.surface_tension().value(), 0.0);
        assert_eq!(critical.temperature_derivative().value(), 0.0);
        let adjacent = state(CRITICAL_TEMPERATURE_K.next_down());
        assert!(adjacent.surface_tension().value() > 0.0);
        assert!(adjacent.temperature_derivative().value() < 0.0);
        assert!(adjacent.surface_tension().value().is_finite());
        assert!(adjacent.temperature_derivative().value().is_finite());
        let coalesced =
            SaturatedWaterState::at_temperature(Temperature::new(CRITICAL_TEMPERATURE_K)).unwrap();
        assert_eq!(coalesced.specific_volume_gap().value(), 0.0);
        assert_eq!(coalesced.surface_tension().unwrap(), critical);
        assert!(
            If97SaturationPair::at_temperature(Temperature::new(CRITICAL_TEMPERATURE_K)).is_err()
        );
    }

    #[test]
    fn g0_existing_saturation_owners_query_the_same_interface_temperature() {
        for t in [TRIPLE_POINT_TEMPERATURE_K, 298.15, 373.15, 623.15] {
            let temperature = Temperature::new(t);
            let supplementary = SaturatedWaterState::at_temperature(temperature).unwrap();
            let industrial = If97SaturationPair::at_temperature(temperature).unwrap();
            let direct = state(t);
            assert_eq!(supplementary.surface_tension().unwrap(), direct);
            assert_eq!(industrial.surface_tension().unwrap(), direct);
            assert_eq!(direct.temperature(), temperature);
        }
        // IF97's optional sub-triple liquid domain does not widen this source.
        let sub_triple = If97SaturationPair::at_temperature(Temperature::new(273.15)).unwrap();
        assert!(matches!(
            sub_triple.surface_tension(),
            Err(WaterSurfaceTensionError::TemperatureOutOfRange { .. })
        ));
    }
}
