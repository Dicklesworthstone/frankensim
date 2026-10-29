//! Original IF97 Region 3 dense-water Helmholtz formulation in coherent SI.
//!
//! Equation (28) and Tables 30–32 of IAPWS R7-97(2012) supply every coefficient
//! and derivative. Attribution: International Association for the Properties
//! of Water and Steam, <https://iapws.org/technical-guidance/release/IF97-Rev.download>.
//! The source uses density and temperature directly. The separate [`pressure`]
//! owner provides bounded pressure inversion and paired saturation endpoints;
//! the source evaluator itself supplies no mixture or metastable continuation.
//!
//! Admission requires the original B23/100 MPa pressure bounds, 623.15–863.15 K,
//! positive Cv and isothermal stiffness, and the stable side of Region 4 below
//! Tc. A conservative density search ceiling lies beyond the 100 MPa isobar
//! throughout this temperature interval (verified with outward intervals).
//! It excludes remote extrapolated roots of the Helmholtz polynomial; positive
//! local compressibility alone would not establish source-domain membership.
//! The exact critical state refuses because the returned Cp/compressibility
//! would be singular. No neighborhood is smoothed or replaced by ideal gas.
//!
//! Source-fit uncertainty remains that of the IAPWS release. In particular,
//! Region 3 and Region 4 are independently fitted original equations: a raw
//! Region 3 state is not moved onto Region 4 by modifying its pressure or energy.

use super::{
    CRITICAL_PRESSURE_PA, CRITICAL_TEMPERATURE_K, IsothermalCompressibility, MAXIMUM_PRESSURE_PA,
    R_SPECIFIC_J_KG_K, SpecificVolumePressureDerivative, SpecificVolumeTemperatureDerivative,
    ThermalExpansion, boundary23_pressure, saturation_pressure,
};
use crate::water::{MassSpecificVolume, SaturationPressureDerivative};
use crate::{MassSpecificEnergyQuantityV1, MassSpecificThermalQuantityV1};
use core::fmt;
use fs_math::det;
use fs_qty::{Density, Pressure, Qty, Temperature, Velocity};

/// Stable fixed-density caloric differences from the same original coefficients.
pub mod increments;

/// Stable changes between original dense-water states at arbitrary endpoints.
pub mod difference;
/// Bounded pressure inversion and original Region 3/4 saturation endpoints.
pub mod pressure;

/// Original critical density and Helmholtz density reduction [kg/m3], Eq. (4).
pub const CRITICAL_DENSITY_KG_M3: f64 = 322.0;
/// Inclusive temperature at the Region 1/3 join [K].
pub const MINIMUM_TEMPERATURE_K: f64 = 623.15;
/// Inclusive temperature at the B23/100 MPa intersection [K].
pub const MAXIMUM_TEMPERATURE_K: f64 = 863.15;
/// Version of the fixed coefficient-power and scaled-derivative operation tree.
pub const EVALUATOR_VERSION: u32 = 1;

// This is a conservative root-search envelope, not a new physical density
// limit. p(2.5*rho_c,T)>100 MPa over the entire Region 3 temperature interval.
const DENSITY_CEILING_KG_M3: f64 = 2.5 * CRITICAL_DENSITY_KG_M3;

/// (dp/drho)_T in Pa/(kg/m3), equivalently m2/s2.
pub type IsothermalPressureDensityDerivative = Qty<2, 0, -2, 0, 0, 0>;
/// (du/drho)_T in (J/kg)/(kg/m3).
pub type InternalEnergyDensityDerivative = Qty<5, -1, -2, 0, 0, 0>;

/// Stable single-phase identity; equality to Tc has no liquid/vapor distinction.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum If97DensePhase {
    Liquid,
    Vapor,
    /// Both temperature and pressure strictly exceed the source critical values.
    Supercritical,
    /// The critical-temperature isotherm away from its singular critical point.
    CriticalIsotherm,
}

/// Source-domain, stability or finite-representation refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum If97DenseError {
    InvalidInput(&'static str),
    OutsideRegion {
        density_kg_m3: f64,
        temperature_k: f64,
        /// None when the input is outside the finite density/temperature envelope.
        pressure_pa: Option<f64>,
    },
    /// A bulk density on the wrong stable side of saturation is not a mixture.
    UnsupportedPhase {
        density_kg_m3: f64,
        temperature_k: f64,
        pressure_pa: f64,
        saturation_pressure_pa: f64,
    },
    Unstable {
        isothermal_stiffness_m2_s2: f64,
        cv_j_kg_k: f64,
    },
    CriticalSingularity,
    Unrepresentable(&'static str),
}
impl fmt::Display for If97DenseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "IAPWS-IF97 Region 3: {self:?}")
    }
}
impl std::error::Error for If97DenseError {}

/// Immutable original density/temperature state; no caloric reference shifts.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct If97DenseState {
    phase: If97DensePhase,
    density_kg_m3: f64,
    temperature_k: f64,
    pressure_pa: f64,
    internal_energy_j_kg: f64,
    enthalpy_j_kg: f64,
    entropy_j_kg_k: f64,
    helmholtz_energy_j_kg: f64,
    gibbs_energy_j_kg: f64,
    cv_j_kg_k: f64,
    cp_j_kg_k: f64,
    pressure_density_derivative: f64,
    pressure_temperature_derivative: f64,
    sound_speed_m_s: f64,
}

impl If97DenseState {
    /// Evaluate the original stable Region 3 source without iteration/allocation.
    pub fn at_density_temperature(
        density: Density,
        temperature: Temperature,
    ) -> Result<Self, If97DenseError> {
        let (rho, t) = (density.value(), temperature.value());
        if !rho.is_finite() || rho <= 0.0 || !t.is_finite() {
            return Err(If97DenseError::InvalidInput(
                "positive finite density and finite temperature required",
            ));
        }
        if !(MINIMUM_TEMPERATURE_K..=MAXIMUM_TEMPERATURE_K).contains(&t)
            || rho > DENSITY_CEILING_KG_M3
        {
            return Err(If97DenseError::OutsideRegion {
                density_kg_m3: rho,
                temperature_k: t,
                pressure_pa: None,
            });
        }
        if rho == CRITICAL_DENSITY_KG_M3 && t == CRITICAL_TEMPERATURE_K {
            return Err(If97DenseError::CriticalSingularity);
        }
        let f = helmholtz(rho / CRITICAL_DENSITY_KG_M3, CRITICAL_TEMPERATURE_K / t);
        let rt = R_SPECIFIC_J_KG_K * t;
        let p = rho * rt * f.d;
        let cv = -R_SPECIFIC_J_KG_K * f.tt;
        let stiffness = rt * (2.0 * f.d + f.dd);
        let coupling = f.d - f.dt;
        let dp_dt = rho * R_SPECIFIC_J_KG_K * coupling;
        let p23 = boundary23_pressure(temperature)
            .map_err(|_| If97DenseError::InvalidInput("B23 temperature"))?
            .value();
        if !p.is_finite() || p < p23 || p > MAXIMUM_PRESSURE_PA {
            return Err(If97DenseError::OutsideRegion {
                density_kg_m3: rho,
                temperature_k: t,
                pressure_pa: Some(p),
            });
        }
        if !cv.is_finite() || !stiffness.is_finite() || cv <= 0.0 || stiffness <= 0.0 {
            return Err(If97DenseError::Unstable {
                isothermal_stiffness_m2_s2: stiffness,
                cv_j_kg_k: cv,
            });
        }
        let phase = if t < CRITICAL_TEMPERATURE_K {
            let ps = saturation_pressure(temperature)
                .map_err(|_| If97DenseError::InvalidInput("Region 4 temperature"))?
                .value();
            if (rho > CRITICAL_DENSITY_KG_M3 && p < ps)
                || (rho < CRITICAL_DENSITY_KG_M3 && p > ps)
                || rho == CRITICAL_DENSITY_KG_M3
            {
                return Err(If97DenseError::UnsupportedPhase {
                    density_kg_m3: rho,
                    temperature_k: t,
                    pressure_pa: p,
                    saturation_pressure_pa: ps,
                });
            }
            if rho > CRITICAL_DENSITY_KG_M3 {
                If97DensePhase::Liquid
            } else {
                If97DensePhase::Vapor
            }
        } else if t == CRITICAL_TEMPERATURE_K {
            If97DensePhase::CriticalIsotherm
        } else if p > CRITICAL_PRESSURE_PA {
            If97DensePhase::Supercritical
        } else {
            If97DensePhase::Vapor
        };
        let u = rt * f.t;
        let h = rt * (f.t + f.d);
        let s = R_SPECIFIC_J_KG_K * (f.t - f.value);
        let free = rt * f.value;
        let gibbs = rt * (f.value + f.d);
        let cp = cv + (t / (rho * rho)) * (dp_dt * dp_dt) / stiffness;
        let sound2 = stiffness + (t / (rho * rho)) * (dp_dt * dp_dt) / cv;
        if ![u, h, s, free, gibbs, cp, sound2, dp_dt]
            .into_iter()
            .all(f64::is_finite)
            || cp <= 0.0
            || sound2 <= 0.0
        {
            return Err(If97DenseError::Unrepresentable(
                "original Helmholtz properties",
            ));
        }
        Ok(Self {
            phase,
            density_kg_m3: rho,
            temperature_k: t,
            pressure_pa: p,
            internal_energy_j_kg: u,
            enthalpy_j_kg: h,
            entropy_j_kg_k: s,
            helmholtz_energy_j_kg: free,
            gibbs_energy_j_kg: gibbs,
            cv_j_kg_k: cv,
            cp_j_kg_k: cp,
            pressure_density_derivative: stiffness,
            pressure_temperature_derivative: dp_dt,
            sound_speed_m_s: det::sqrt(sound2),
        })
    }
    pub const fn phase(self) -> If97DensePhase {
        self.phase
    }
    pub const fn density(self) -> Density {
        Density::new(self.density_kg_m3)
    }
    pub const fn temperature(self) -> Temperature {
        Temperature::new(self.temperature_k)
    }
    pub const fn pressure(self) -> Pressure {
        Pressure::new(self.pressure_pa)
    }
    pub fn specific_volume(self) -> MassSpecificVolume {
        MassSpecificVolume::new(1.0 / self.density_kg_m3)
    }
    pub const fn specific_internal_energy(self) -> MassSpecificEnergyQuantityV1 {
        MassSpecificEnergyQuantityV1::new(self.internal_energy_j_kg)
    }
    pub const fn specific_enthalpy(self) -> MassSpecificEnergyQuantityV1 {
        MassSpecificEnergyQuantityV1::new(self.enthalpy_j_kg)
    }
    pub const fn specific_entropy(self) -> MassSpecificThermalQuantityV1 {
        MassSpecificThermalQuantityV1::new(self.entropy_j_kg_k)
    }
    pub const fn specific_helmholtz_energy(self) -> MassSpecificEnergyQuantityV1 {
        MassSpecificEnergyQuantityV1::new(self.helmholtz_energy_j_kg)
    }
    pub const fn specific_gibbs_energy(self) -> MassSpecificEnergyQuantityV1 {
        MassSpecificEnergyQuantityV1::new(self.gibbs_energy_j_kg)
    }
    pub const fn cv(self) -> MassSpecificThermalQuantityV1 {
        MassSpecificThermalQuantityV1::new(self.cv_j_kg_k)
    }
    pub const fn cp(self) -> MassSpecificThermalQuantityV1 {
        MassSpecificThermalQuantityV1::new(self.cp_j_kg_k)
    }
    pub const fn sound_speed(self) -> Velocity {
        Velocity::new(self.sound_speed_m_s)
    }
    pub const fn dp_drho_at_temperature(self) -> IsothermalPressureDensityDerivative {
        IsothermalPressureDensityDerivative::new(self.pressure_density_derivative)
    }
    pub const fn dp_dt_at_density(self) -> SaturationPressureDerivative {
        SaturationPressureDerivative::new(self.pressure_temperature_derivative)
    }
    pub fn du_drho_at_temperature(self) -> InternalEnergyDensityDerivative {
        InternalEnergyDensityDerivative::new(
            (self.pressure_pa - self.temperature_k * self.pressure_temperature_derivative)
                / (self.density_kg_m3 * self.density_kg_m3),
        )
    }
    /// (dh/dT)_rho=Cv+(dp/dT)_rho/rho, distinct from isobaric Cp.
    pub fn dh_dt_at_density(self) -> MassSpecificThermalQuantityV1 {
        MassSpecificThermalQuantityV1::new(
            self.cv_j_kg_k + self.pressure_temperature_derivative / self.density_kg_m3,
        )
    }
    pub fn dv_dt_at_pressure(self) -> SpecificVolumeTemperatureDerivative {
        SpecificVolumeTemperatureDerivative::new(
            self.pressure_temperature_derivative
                / (self.density_kg_m3 * self.density_kg_m3 * self.pressure_density_derivative),
        )
    }
    pub fn dv_dp_at_temperature(self) -> SpecificVolumePressureDerivative {
        SpecificVolumePressureDerivative::new(
            -1.0 / (self.density_kg_m3 * self.density_kg_m3 * self.pressure_density_derivative),
        )
    }
    pub fn thermal_expansion(self) -> ThermalExpansion {
        ThermalExpansion::new(
            self.pressure_temperature_derivative
                / (self.density_kg_m3 * self.pressure_density_derivative),
        )
    }
    pub fn isothermal_compressibility(self) -> IsothermalCompressibility {
        IsothermalCompressibility::new(
            1.0 / (self.density_kg_m3 * self.pressure_density_derivative),
        )
    }
}

// Scaled derivatives avoid dividing each of the 39 polynomial terms by density
// or temperature. Each cached nonnegative power is formed once per source query.
#[derive(Default)]
struct Helmholtz {
    value: f64,
    d: f64,  // delta*phi_delta
    dd: f64, // delta^2*phi_delta_delta
    t: f64,  // tau*phi_tau
    tt: f64, // tau^2*phi_tau_tau
    dt: f64, // delta*tau*phi_delta_tau
}
fn helmholtz(delta: f64, tau: f64) -> Helmholtz {
    let mut density_power = [1.0; 12];
    let mut temperature_power = [1.0; 27];
    for i in 1..density_power.len() {
        density_power[i] = density_power[i - 1] * delta;
    }
    for j in 1..temperature_power.len() {
        temperature_power[j] = temperature_power[j - 1] * tau;
    }
    let mut out = Helmholtz {
        value: LOG_COEFFICIENT * det::ln(delta),
        d: LOG_COEFFICIENT,
        dd: -LOG_COEFFICIENT,
        ..Helmholtz::default()
    };
    for &(i, j, n) in &TERMS {
        let term = n * density_power[i] * temperature_power[j];
        let di = i as f64;
        let tj = j as f64;
        out.value += term;
        out.d += di * term;
        out.dd += di * (di - 1.0) * term;
        out.t += tj * term;
        out.tt += tj * (tj - 1.0) * term;
        out.dt += di * tj * term;
    }
    out
}

// IAPWS R7-97(2012), Table 30. Term 1 multiplies ln(delta); terms 2–40
// multiply delta^I*tau^J in the released table order.
const LOG_COEFFICIENT: f64 = 1.0658070028513;
const TERMS: [(usize, usize, f64); 39] = [
    (0, 0, -15.732845290239),
    (0, 1, 20.944396974307),
    (0, 2, -7.6867707878716),
    (0, 7, 2.6185947787954),
    (0, 10, -2.8080781148620),
    (0, 12, 1.2053369696517),
    (0, 23, -0.0084566812812502),
    (1, 2, -1.2654315477714),
    (1, 6, -1.1524407806681),
    (1, 15, 0.88521043984318),
    (1, 17, -0.64207765181607),
    (2, 0, 0.38493460186671),
    (2, 2, -0.85214708824206),
    (2, 6, 4.8972281541877),
    (2, 7, -3.0502617256965),
    (2, 22, 0.039420536879154),
    (2, 26, 0.12558408424308),
    (3, 0, -0.27999329698710),
    (3, 2, 1.3899799569460),
    (3, 4, -2.0189915023570),
    (3, 16, -0.0082147637173963),
    (3, 26, -0.47596035734923),
    (4, 0, 0.043984074473500),
    (4, 2, -0.44476435428739),
    (4, 4, 0.90572070719733),
    (4, 26, 0.70522450087967),
    (5, 1, 0.10770512626332),
    (5, 3, -0.32913623258954),
    (5, 26, -0.50871062041158),
    (6, 0, -0.022175400873096),
    (6, 2, 0.094260751665092),
    (6, 26, 0.16436278447961),
    (7, 2, -0.013503372241348),
    (8, 26, -0.014834345352472),
    (9, 2, 0.00057922953628084),
    (9, 26, 0.0032308904703711),
    (10, 0, 0.000080964802996215),
    (10, 1, -0.00016557679795037),
    (11, 26, -0.000044923899061815),
];

#[cfg(test)]
mod tests;
