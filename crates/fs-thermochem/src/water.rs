//! Saturation properties of ordinary water from IAPWS SR1-86(1992).
//!
//! Implements the ITS-90 correlations in equations (1)--(9), over the release's
//! inclusive triple-to-critical interval. Pressure, phase densities, enthalpy
//! and entropy share this one source; internal energy is h-p/rho. This is not
//! an off-saturation equation of state, a metastability model, a transport law,
//! or an evaporation/condensation time integrator.
//!
//! Source: IAPWS, Revised Supplementary Release on Saturation Properties of
//! Ordinary Water Substance (September 1992),
//! <https://iapws.org/public/documents/6dGkr/Supp-sat.pdf>.
//! The release refers property uncertainties to the IAPWS skeleton tables;
//! numerical inversion tolerances below are not material uncertainty bounds.
//! The u_liquid=s_liquid=0 triple-point reference is retained to the precision
//! of the published coefficients, without a separate numerical datum shift.

use crate::{MassSpecificEnergyQuantityV1, MassSpecificThermalQuantityV1};
use core::fmt;
use fs_math::det;
use fs_qty::{Density, Pressure, Qty, Temperature};

/// Fixed source release and temperature scale.
pub const SOURCE_RELEASE: &str = "IAPWS SR1-86(1992), ITS-90";
/// Source equations and reference table used by this evaluator.
pub const SOURCE_URL: &str = "https://iapws.org/public/documents/6dGkr/Supp-sat.pdf";
/// Version of the deterministic source evaluation and common phase-gap tree.
pub const EVALUATOR_VERSION: u32 = 1;
/// Lower inclusive source temperature [K], the liquid/ice-I/vapor triple point.
pub const TRIPLE_POINT_TEMPERATURE_K: f64 = 273.16;
/// Upper inclusive source temperature [K]. Two-phase quality is undefined here.
pub const CRITICAL_TEMPERATURE_K: f64 = 647.096;
/// Critical pressure [Pa] from release section 2.
pub const CRITICAL_PRESSURE_PA: f64 = 22.064e6;
/// Critical mass density [kg/m3] from release section 2.
pub const CRITICAL_DENSITY_KG_M3: f64 = 322.0;
/// Hard work ceiling for the scalar inverse, in addition to the caller's budget.
pub const MAX_PRESSURE_ITERATIONS: usize = 256;

/// Pressure slope along saturation [Pa/K], not an isochoric EOS derivative.
pub type SaturationPressureDerivative = Qty<-1, 1, -2, -1, 0, 0>;
/// Mass-specific volume [m3/kg].
pub type MassSpecificVolume = Qty<3, -1, 0, 0, 0, 0>;

/// Source-domain, numerical-work, or cancellation refusal; no state is mutated.
#[derive(Clone, Debug, PartialEq)]
pub enum WaterSaturationError {
    /// Nonfinite temperature or a temperature outside the inclusive source range.
    TemperatureOutOfRange { temperature_k: f64 },
    /// Nonfinite pressure or pressure outside the evaluated endpoint pressures.
    PressureOutOfRange {
        pressure_pa: f64,
        minimum_pa: f64,
        maximum_pa: f64,
    },
    /// Invalid caller accuracy or work bound.
    InvalidSolve(&'static str),
    /// A finite, ordered saturation pair could not be represented.
    Unrepresentable,
    /// Work or representable temperature resolution exhausted before acceptance.
    NotConverged {
        iterations: usize,
        pressure_residual_pa: f64,
    },
    /// Caller cancellation observed before publishing a result.
    Cancelled,
}
impl fmt::Display for WaterSaturationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "IAPWS saturation: {self:?}")
    }
}
impl std::error::Error for WaterSaturationError {}

/// Explicit numerical accuracy and bounded work for T(p) on the same source.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterSaturationSolve {
    /// Absolute residual of the ORIGINAL pressure correlation [Pa].
    pub pressure_tolerance_pa: f64,
    /// Positive and no larger than MAX_PRESSURE_ITERATIONS.
    pub maximum_iterations: usize,
}

/// Accepted inverse state and its actual original-source residual.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterSaturationSolution {
    state: SaturatedWaterState,
    iterations: usize,
    pressure_residual_pa: f64,
}
impl WaterSaturationSolution {
    /// Original source state accepted by the inverse.
    pub const fn state(&self) -> &SaturatedWaterState {
        &self.state
    }
    /// Consume the solve receipt and retain its immutable state.
    pub const fn into_state(self) -> SaturatedWaterState {
        self.state
    }
    /// Pressure evaluations after the endpoint check; zero at exact endpoints.
    pub const fn iterations(&self) -> usize {
        self.iterations
    }
    /// Signed evaluated-minus-requested pressure [Pa].
    pub const fn pressure_residual(&self) -> Pressure {
        Pressure::new(self.pressure_residual_pa)
    }
}

/// Validated immutable liquid/vapor saturation pair, with common caloric datum.
/// At the critical endpoint both phases coalesce and every latent gap is zero.
/// Consumers requiring distinguishable phase inventories must refuse that state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SaturatedWaterState {
    temperature_k: f64,
    pressure_pa: f64,
    pressure_derivative_pa_k: f64,
    liquid_density_kg_m3: f64,
    vapor_density_kg_m3: f64,
    liquid_enthalpy_j_kg: f64,
    liquid_entropy_j_kg_k: f64,
    liquid_internal_energy_j_kg: f64,
    specific_volume_gap_m3_kg: f64,
    latent_enthalpy_j_kg: f64,
    latent_entropy_j_kg_k: f64,
    latent_internal_energy_j_kg: f64,
}

impl SaturatedWaterState {
    /// Evaluate equations (1)--(9) without allocation or extrapolation.
    /// Fractional powers use the deterministic math owner. The explicit phase
    /// gaps cancel the common auxiliary functions before evaluation; their
    /// direct values can differ by rounding from subtracting reported totals.
    pub fn at_temperature(temperature: Temperature) -> Result<Self, WaterSaturationError> {
        let t = temperature.value();
        if !t.is_finite() || !(TRIPLE_POINT_TEMPERATURE_K..=CRITICAL_TEMPERATURE_K).contains(&t) {
            return Err(WaterSaturationError::TemperatureOutOfRange { temperature_k: t });
        }
        let tau = (CRITICAL_TEMPERATURE_K - t) / CRITICAL_TEMPERATURE_K;
        let theta = t / CRITICAL_TEMPERATURE_K;
        let (pressure, dp) = pressure_and_derivative(t);

        // One sixth-root supplies every rational exponent in equations (2),(3).
        let q = det::pow(tau, 1.0 / 6.0);
        let q2 = q * q;
        let q4 = q2 * q2;
        let q8 = q4 * q4;
        let q16 = q8 * q8;
        let q32 = q16 * q16;
        let q64 = q32 * q32;
        let q128 = q64 * q64;
        let liquid_excess = 1.99274064 * q2 + 1.09965342 * q4
            - 0.510839303 * q8 * q2
            - 1.75493479 * q32
            - 45.5170352 * q64 * q16 * q4 * q2
            - 6.74694450e5 * q128 * q64 * q16 * q8 * q4;
        let vapor_log_ratio = -2.03150240 * q2
            - 2.68302940 * q4
            - 5.38626492 * q8
            - 17.2991605 * q16 * q2
            - 44.7586581 * q32 * q4 * q
            - 63.9201063 * q64 * q4 * q2 * q;
        let liquid_density = CRITICAL_DENSITY_KG_M3 * (1.0 + liquid_excess);
        let vapor_density = CRITICAL_DENSITY_KG_M3 * det::exp(vapor_log_ratio);
        // v_v-v_l = [expm1(-log(rho_v/rho_c)) + B/(1+B)]/rho_c.
        // This keeps the small positive gap near critical coalescence.
        let volume_gap = (det::expm1(-vapor_log_ratio) + liquid_excess / (1.0 + liquid_excess))
            / CRITICAL_DENSITY_KG_M3;

        let root_theta = det::sqrt(theta);
        let theta2 = theta * theta;
        let theta4 = theta2 * theta2;
        let theta_4_5 = theta4 * root_theta;
        let theta_54_5 = det::powi(theta, 54) * root_theta;
        let theta_neg19 = det::powi(theta, -19);
        // Keep every published digit of d_alpha: liquid h near the triple point
        // depends on cancellation in this source-defined reference convention.
        let alpha = 1000.0
            * (-1135.905627715 - 5.65134998e-8 * theta_neg19
                + 2690.66631 * theta
                + 127.287297 * theta_4_5
                - 135.003439 * theta4 * theta
                + 0.981825814 * theta_54_5);
        let phi = (1000.0 / CRITICAL_TEMPERATURE_K)
            * (2319.5246
                + (19.0 / 20.0) * (-5.65134998e-8) * theta_neg19 / theta
                + 2690.66631 * det::ln(theta)
                + (9.0 / 7.0) * 127.287297 * theta_4_5 / theta
                + (5.0 / 4.0) * (-135.003439) * theta4
                + (109.0 / 107.0) * 0.981825814 * theta_54_5 / theta);
        let liquid_enthalpy = alpha + t * dp / liquid_density;
        let liquid_entropy = phi + dp / liquid_density;
        let liquid_internal_energy = liquid_enthalpy - pressure / liquid_density;
        let latent_entropy = dp * volume_gap;
        let latent_enthalpy = t * latent_entropy;
        let latent_internal_energy = (t * dp - pressure) * volume_gap;
        if ![
            pressure,
            dp,
            liquid_density,
            vapor_density,
            liquid_enthalpy,
            liquid_entropy,
            liquid_internal_energy,
            volume_gap,
            latent_entropy,
            latent_enthalpy,
            latent_internal_energy,
            liquid_enthalpy + latent_enthalpy,
            liquid_entropy + latent_entropy,
            liquid_internal_energy + latent_internal_energy,
        ]
        .iter()
        .all(|value| value.is_finite())
            || pressure <= 0.0
            || dp <= 0.0
            || vapor_density <= 0.0
            || liquid_density < vapor_density
            || volume_gap < 0.0
            || latent_enthalpy < 0.0
        {
            return Err(WaterSaturationError::Unrepresentable);
        }
        Ok(Self {
            temperature_k: t,
            pressure_pa: pressure,
            pressure_derivative_pa_k: dp,
            liquid_density_kg_m3: liquid_density,
            vapor_density_kg_m3: vapor_density,
            liquid_enthalpy_j_kg: liquid_enthalpy,
            liquid_entropy_j_kg_k: liquid_entropy,
            liquid_internal_energy_j_kg: liquid_internal_energy,
            specific_volume_gap_m3_kg: volume_gap,
            latent_enthalpy_j_kg: latent_enthalpy,
            latent_entropy_j_kg_k: latent_entropy,
            latent_internal_energy_j_kg: latent_internal_energy,
        })
    }

    /// Invert the original monotone saturation pressure with safeguarded Newton
    /// steps and bisection fallback. The endpoint remains an original T-based
    /// source evaluation; the requested pressure is not patched into its fields.
    /// The lower pressure limit is the evaluated release at 273.16 K, rather
    /// than its rounded table entry. Every iteration polls cancellation.
    pub fn at_pressure(
        pressure: Pressure,
        controls: WaterSaturationSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<WaterSaturationSolution, WaterSaturationError> {
        if !controls.pressure_tolerance_pa.is_finite()
            || controls.pressure_tolerance_pa <= 0.0
            || controls.maximum_iterations == 0
            || controls.maximum_iterations > MAX_PRESSURE_ITERATIONS
        {
            return Err(WaterSaturationError::InvalidSolve(
                "positive finite pressure tolerance and bounded iterations required",
            ));
        }
        poll(&mut cancelled)?;
        let target = pressure.value();
        let minimum = pressure_and_derivative(TRIPLE_POINT_TEMPERATURE_K).0;
        if !target.is_finite() || target < minimum || target > CRITICAL_PRESSURE_PA {
            return Err(WaterSaturationError::PressureOutOfRange {
                pressure_pa: target,
                minimum_pa: minimum,
                maximum_pa: CRITICAL_PRESSURE_PA,
            });
        }
        if target == minimum || target == CRITICAL_PRESSURE_PA {
            let t = if target == minimum {
                TRIPLE_POINT_TEMPERATURE_K
            } else {
                CRITICAL_TEMPERATURE_K
            };
            let state = Self::at_temperature(Temperature::new(t))?;
            poll(&mut cancelled)?;
            return Ok(WaterSaturationSolution {
                state,
                iterations: 0,
                pressure_residual_pa: 0.0,
            });
        }
        let mut lower = TRIPLE_POINT_TEMPERATURE_K;
        let mut upper = CRITICAL_TEMPERATURE_K;
        // A reciprocal-temperature seed interpolates only the source endpoints.
        // It changes convergence speed, never the accepted thermodynamic law.
        let fraction = det::ln(target / minimum) / det::ln(CRITICAL_PRESSURE_PA / minimum);
        let mut t = 1.0 / ((1.0 - fraction) / lower + fraction / upper);
        let mut last_residual = f64::INFINITY;
        for iterations in 1..=controls.maximum_iterations {
            poll(&mut cancelled)?;
            let (evaluated, slope) = pressure_and_derivative(t);
            let residual = evaluated - target;
            last_residual = residual;
            if residual.abs() <= controls.pressure_tolerance_pa {
                let state = Self::at_temperature(Temperature::new(t))?;
                poll(&mut cancelled)?;
                return Ok(WaterSaturationSolution {
                    state,
                    iterations,
                    pressure_residual_pa: residual,
                });
            }
            if residual < 0.0 {
                lower = t;
            } else {
                upper = t;
            }
            let newton = t - residual / slope;
            let next = if newton.is_finite() && newton > lower && newton < upper && newton != t {
                newton
            } else {
                lower + 0.5 * (upper - lower)
            };
            if next <= lower || next >= upper {
                return Err(WaterSaturationError::NotConverged {
                    iterations,
                    pressure_residual_pa: residual,
                });
            }
            t = next;
        }
        Err(WaterSaturationError::NotConverged {
            iterations: controls.maximum_iterations,
            pressure_residual_pa: last_residual,
        })
    }

    /// Saturation temperature on the ITS-90 scale [K].
    pub const fn temperature(self) -> Temperature {
        Temperature::new(self.temperature_k)
    }
    /// Evaluated equilibrium vapor pressure [Pa].
    pub const fn pressure(self) -> Pressure {
        Pressure::new(self.pressure_pa)
    }
    /// Analytic derivative of the source pressure correlation [Pa/K].
    pub const fn pressure_derivative(self) -> SaturationPressureDerivative {
        SaturationPressureDerivative::new(self.pressure_derivative_pa_k)
    }
    /// Saturated liquid density [kg/m3].
    pub const fn liquid_density(self) -> Density {
        Density::new(self.liquid_density_kg_m3)
    }
    /// Saturated vapor density [kg/m3].
    pub const fn vapor_density(self) -> Density {
        Density::new(self.vapor_density_kg_m3)
    }
    /// Saturated liquid enthalpy in the release's common phase datum [J/kg].
    pub const fn liquid_specific_enthalpy(self) -> MassSpecificEnergyQuantityV1 {
        MassSpecificEnergyQuantityV1::new(self.liquid_enthalpy_j_kg)
    }
    /// Saturated vapor enthalpy, liquid enthalpy plus the direct latent gap [J/kg].
    pub fn vapor_specific_enthalpy(self) -> MassSpecificEnergyQuantityV1 {
        MassSpecificEnergyQuantityV1::new(self.liquid_enthalpy_j_kg + self.latent_enthalpy_j_kg)
    }
    /// Saturated liquid entropy in the release's common phase datum [J/(kg K)].
    pub const fn liquid_specific_entropy(self) -> MassSpecificThermalQuantityV1 {
        MassSpecificThermalQuantityV1::new(self.liquid_entropy_j_kg_k)
    }
    /// Saturated vapor entropy, liquid entropy plus the direct latent gap [J/(kg K)].
    pub fn vapor_specific_entropy(self) -> MassSpecificThermalQuantityV1 {
        MassSpecificThermalQuantityV1::new(self.liquid_entropy_j_kg_k + self.latent_entropy_j_kg_k)
    }
    /// Saturated liquid internal energy, h_l-p/rho_l [J/kg].
    pub const fn liquid_specific_internal_energy(self) -> MassSpecificEnergyQuantityV1 {
        MassSpecificEnergyQuantityV1::new(self.liquid_internal_energy_j_kg)
    }
    /// Saturated vapor internal energy, u_l plus the direct internal-energy gap [J/kg].
    pub fn vapor_specific_internal_energy(self) -> MassSpecificEnergyQuantityV1 {
        MassSpecificEnergyQuantityV1::new(
            self.liquid_internal_energy_j_kg + self.latent_internal_energy_j_kg,
        )
    }
    /// Positive vapor-minus-liquid specific volume, evaluated before cancellation [m3/kg].
    pub const fn specific_volume_gap(self) -> MassSpecificVolume {
        MassSpecificVolume::new(self.specific_volume_gap_m3_kg)
    }
    /// Vaporization enthalpy T(dp/dT)(v_v-v_l), zero at criticality [J/kg].
    pub const fn latent_enthalpy(self) -> MassSpecificEnergyQuantityV1 {
        MassSpecificEnergyQuantityV1::new(self.latent_enthalpy_j_kg)
    }
    /// Vaporization entropy (dp/dT)(v_v-v_l), zero at criticality [J/(kg K)].
    pub const fn latent_entropy(self) -> MassSpecificThermalQuantityV1 {
        MassSpecificThermalQuantityV1::new(self.latent_entropy_j_kg_k)
    }
    /// Vaporization internal-energy change (T dp/dT-p)(v_v-v_l) [J/kg].
    pub const fn latent_internal_energy(self) -> MassSpecificEnergyQuantityV1 {
        MassSpecificEnergyQuantityV1::new(self.latent_internal_energy_j_kg)
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), WaterSaturationError> {
    if cancelled() {
        Err(WaterSaturationError::Cancelled)
    } else {
        Ok(())
    }
}

// Equation (1) and its analytic derivative. Only admitted temperatures enter.
fn pressure_and_derivative(t: f64) -> (f64, f64) {
    let tau = (CRITICAL_TEMPERATURE_K - t) / CRITICAL_TEMPERATURE_K;
    let root = det::sqrt(tau);
    let tau2 = tau * tau;
    let tau3 = tau2 * tau;
    let tau4 = tau2 * tau2;
    let polynomial = -7.85951783 * tau + 1.84408259 * tau * root - 11.7866497 * tau3
        + 22.6807411 * tau3 * root
        - 15.9618719 * tau4
        + 1.80122502 * tau3 * tau4 * root;
    let derivative_tau = -7.85951783 + 1.5 * 1.84408259 * root - 3.0 * 11.7866497 * tau2
        + 3.5 * 22.6807411 * tau2 * root
        - 4.0 * 15.9618719 * tau3
        + 7.5 * 1.80122502 * tau3 * tau3 * root;
    let log_ratio = (CRITICAL_TEMPERATURE_K / t) * polynomial;
    let pressure = CRITICAL_PRESSURE_PA * det::exp(log_ratio);
    (pressure, -(pressure / t) * (log_ratio + derivative_tau))
}

#[cfg(test)]
mod tests;
