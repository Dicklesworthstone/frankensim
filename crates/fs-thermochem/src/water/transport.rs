//! Pressure-dependent industrial water/steam transport on the original IF97 EOS.
//!
//! Viscosity uses IAPWS R12-08, equations (10)--(12), Tables 1--2 and section 3:
//! the industrial formulation explicitly sets the viscosity critical factor to
//! one. Conductivity uses IAPWS R15-11, equations (15)--(26), Tables 1--3 and 6,
//! including its INDUSTRIAL critical enhancement. This is not the distinct
//! IAPWS-95 general/scientific formulation near the critical point.
//!
//! The immutable IF97 owner supplies density, Cp, Cv and isothermal response at
//! the same pressure/temperature. Its original region, phase, caloric datum,
//! source matching defects and dense-root diagnostics remain available. No
//! reference-pressure liquid fit or dilute-vapor surrogate is substituted.
//! The scope is the existing stable IF97 Regions 1--3, 273.16--1073.15 K and
//! at most 100 MPa. Saturation sides must be explicitly selected by that owner;
//! a two-phase bulk quality is not assigned single-phase transport. The exact
//! critical singularity, ice, metastable continuations and Region 5 remain out
//! of scope. The releases' source uncertainty is separate from numerical error.
//!
//! Attribution: International Association for the Properties of Water and Steam.
//! Publication with attribution is permitted by both original releases:
//! <https://iapws.org/documents/release/viscosity.download> (2008), and
//! <https://iapws.org/documents/release/ThCond.download> (2011).

use super::if97::region3::pressure::DensePressureSolve;
use super::if97::regions::{If97FluidError, If97FluidState};
use core::fmt;
use fs_math::det;
use fs_qty::{Pressure, Qty, Temperature};

pub const VISCOSITY_SOURCE: &str = "IAPWS R12-08 (2008), industrial section 3";
pub const CONDUCTIVITY_SOURCE: &str = "IAPWS R15-11 (2011), industrial section 3";
pub const MINIMUM_TEMPERATURE_K: f64 = 273.16;
pub const MAXIMUM_TEMPERATURE_K: f64 = 1073.15;
pub const MAXIMUM_PRESSURE_PA: f64 = 100.0e6;
const REFERENCE_TEMPERATURE_K: f64 = 647.096;
const REFERENCE_DENSITY_KG_M3: f64 = 322.0;
const REFERENCE_PRESSURE_PA: f64 = 22.064e6;
// R15-11's own reference value, not IF97's different thermodynamic R.
const CONDUCTIVITY_GAS_CONSTANT_J_KG_K: f64 = 461.51805;

/// Dynamic shear viscosity [Pa s]. No bulk-viscosity coefficient is inferred.
pub type DynamicViscosity = Qty<-1, 1, -1, 0, 0, 0>;
/// Thermal conductivity [W/(m K)].
pub type ThermalConductivity = Qty<1, 1, -3, -1, 0, 0>;

#[derive(Clone, Debug, PartialEq)]
pub enum WaterTransportError {
    Thermodynamics(If97FluidError),
    OutsideDomain,
    Unrepresentable,
    Cancelled,
}
impl From<If97FluidError> for WaterTransportError {
    fn from(error: If97FluidError) -> Self {
        match error {
            If97FluidError::Cancelled => Self::Cancelled,
            other => Self::Thermodynamics(other),
        }
    }
}
impl fmt::Display for WaterTransportError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "IAPWS industrial water transport: {self:?}")
    }
}
impl std::error::Error for WaterTransportError {}

/// Original EOS plus its named industrial transport laws at the SAME state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct If97TransportState {
    thermodynamics: If97FluidState,
    viscosity_pa_s: f64,
    background_conductivity_w_m_k: f64,
    critical_conductivity_w_m_k: f64,
    conductivity_w_m_k: f64,
}

impl If97TransportState {
    /// Bind an already accepted original regional or explicit saturation-side
    /// state. This performs no second EOS solve or phase selection.
    pub fn from_state(thermodynamics: If97FluidState) -> Result<Self, WaterTransportError> {
        let t = thermodynamics.temperature().value();
        let p = thermodynamics.pressure().value();
        if !t.is_finite()
            || !(MINIMUM_TEMPERATURE_K..=MAXIMUM_TEMPERATURE_K).contains(&t)
            || !p.is_finite()
            || p <= 0.0
            || p > MAXIMUM_PRESSURE_PA
        {
            return Err(WaterTransportError::OutsideDomain);
        }
        let rho = thermodynamics.density().value();
        let cp = thermodynamics.cp().value();
        let cv = thermodynamics.cv().value();
        let beta = thermodynamics.isothermal_compressibility().value();
        if ![rho, cp, cv, beta].into_iter().all(positive) || cv > cp {
            return Err(WaterTransportError::Unrepresentable);
        }
        let viscosity_pa_s = viscosity(t, rho);
        let background_conductivity_w_m_k = background_conductivity(t, rho);
        let critical_conductivity_w_m_k =
            critical_conductivity(t, rho, cp, cv, beta, viscosity_pa_s);
        let conductivity_w_m_k = background_conductivity_w_m_k + critical_conductivity_w_m_k;
        if ![
            viscosity_pa_s,
            background_conductivity_w_m_k,
            conductivity_w_m_k,
        ]
        .into_iter()
        .all(positive)
            || !critical_conductivity_w_m_k.is_finite()
            || critical_conductivity_w_m_k < 0.0
        {
            return Err(WaterTransportError::Unrepresentable);
        }
        Ok(Self {
            thermodynamics,
            viscosity_pa_s,
            background_conductivity_w_m_k,
            critical_conductivity_w_m_k,
            conductivity_w_m_k,
        })
    }

    /// Resolve one pressure/temperature through the existing stable IF97 router
    /// and then attach transport. Work bounds/cancellation stay on that router;
    /// refused states publish neither a thermodynamic nor transport candidate.
    pub fn at_pressure_temperature(
        pressure: Pressure,
        temperature: Temperature,
        controls: DensePressureSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, WaterTransportError> {
        if cancelled() {
            return Err(WaterTransportError::Cancelled);
        }
        let source = If97FluidState::at_pressure_temperature(
            pressure,
            temperature,
            controls,
            &mut cancelled,
        )?;
        let result = Self::from_state(source)?;
        if cancelled() {
            return Err(WaterTransportError::Cancelled);
        }
        Ok(result)
    }

    /// Original source state, including region, phase and dense-solve diagnostics.
    pub const fn thermodynamics(&self) -> &If97FluidState {
        &self.thermodynamics
    }
    /// Industrial shear viscosity, with the source's critical factor set to one.
    pub const fn dynamic_viscosity(&self) -> DynamicViscosity {
        DynamicViscosity::new(self.viscosity_pa_s)
    }
    /// Total industrial conductivity, including the actual critical enhancement.
    pub const fn thermal_conductivity(&self) -> ThermalConductivity {
        ThermalConductivity::new(self.conductivity_w_m_k)
    }
    /// Dilute-times-finite-density contribution, retained separately for diagnosis.
    pub const fn background_thermal_conductivity(&self) -> ThermalConductivity {
        ThermalConductivity::new(self.background_conductivity_w_m_k)
    }
    /// Additive industrial critical contribution; this is not a heat source.
    pub const fn critical_thermal_conductivity(&self) -> ThermalConductivity {
        ThermalConductivity::new(self.critical_conductivity_w_m_k)
    }
}

fn positive(v: f64) -> bool {
    v.is_finite() && v > 0.0
}
fn polynomial<const N: usize>(x: f64, coefficients: &[f64; N]) -> f64 {
    coefficients
        .iter()
        .rev()
        .fold(0.0, |value, &coefficient| value * x + coefficient)
}

// R12-08 Tables 1 and 2. Rows are i, columns j; omitted terms are zero.
const MU_DILUTE: [f64; 4] = [1.67752, 2.20462, 0.6366564, -0.241605];
const MU_DENSE: [[f64; 7]; 6] = [
    [
        0.520094, 0.222531, -0.281378, 0.161913, -0.0325372, 0.0, 0.0,
    ],
    [0.0850895, 0.999115, -0.906851, 0.257399, 0.0, 0.0, 0.0],
    [-1.08374, 1.88797, -0.772479, 0.0, 0.0, 0.0, 0.0],
    [
        -0.289555,
        1.26613,
        -0.489837,
        0.0,
        0.0698452,
        0.0,
        -0.00435673,
    ],
    [0.0, 0.0, -0.257040, 0.0, 0.0, 0.00872102, 0.0],
    [0.0, 0.120573, 0.0, 0.0, 0.0, 0.0, -0.000593264],
];
fn viscosity(t: f64, rho: f64) -> f64 {
    let theta = t / REFERENCE_TEMPERATURE_K;
    let delta = rho / REFERENCE_DENSITY_KG_M3;
    let dilute = 100.0 * det::sqrt(theta) / polynomial(theta.recip(), &MU_DILUTE);
    let residual = MU_DENSE.iter().rev().fold(0.0, |value, row| {
        value * (theta.recip() - 1.0) + polynomial(delta - 1.0, row)
    });
    1e-6 * dilute * det::exp(delta * residual)
}

// R15-11 Tables 1 and 2. Same i/j convention as the original release.
const K_DILUTE: [f64; 5] = [
    2.443221e-3,
    1.323095e-2,
    6.770357e-3,
    -3.454586e-3,
    4.096266e-4,
];
const K_DENSE: [[f64; 6]; 5] = [
    [
        1.60397357,
        -0.646013523,
        0.111443906,
        0.102997357,
        -0.0504123634,
        0.00609859258,
    ],
    [
        2.33771842,
        -2.78843778,
        1.53616167,
        -0.463045512,
        0.0832827019,
        -0.00719201245,
    ],
    [
        2.19650529,
        -4.54580785,
        3.55777244,
        -1.40944978,
        0.275418278,
        -0.0205938816,
    ],
    [
        -1.21051378,
        1.60812989,
        -0.621178141,
        0.0716373224,
        0.0,
        0.0,
    ],
    [
        -2.7203370,
        4.57586331,
        -3.18369245,
        1.1168348,
        -0.19268305,
        0.012913842,
    ],
];
fn background_conductivity(t: f64, rho: f64) -> f64 {
    let theta = t / REFERENCE_TEMPERATURE_K;
    let delta = rho / REFERENCE_DENSITY_KG_M3;
    let dilute = det::sqrt(theta) / polynomial(theta.recip(), &K_DILUTE);
    let residual = K_DENSE.iter().rev().fold(0.0, |value, row| {
        value * (theta.recip() - 1.0) + polynomial(delta - 1.0, row)
    });
    1e-3 * dilute * det::exp(delta * residual)
}

// R15-11 Table 6. Each row is one density interval j and contains A_0j..A_5j.
const REFERENCE_RESPONSE: [[f64; 6]; 5] = [
    [
        6.53786807199516,
        -5.61149954923348,
        3.39624167361325,
        -2.27492629730878,
        10.2631854662709,
        1.97815050331519,
    ],
    [
        6.52717759281799,
        -6.30816983387575,
        8.08379285492595,
        -9.82240510197603,
        12.1358413791395,
        -5.54349664571295,
    ],
    [
        5.35500529896124,
        -3.96415689925446,
        8.91990208918795,
        -12.0338729505790,
        9.19494865194302,
        -2.16866274479712,
    ],
    [
        1.55225959906681,
        0.464621290821181,
        8.93237374861479,
        -11.0321960061126,
        6.16780999933360,
        -0.965458722086812,
    ],
    [
        1.11999926419994,
        0.595748562571649,
        9.88952565078920,
        -10.3255051147040,
        4.66861294457414,
        -0.503243546373828,
    ],
];
fn critical_conductivity(t: f64, rho: f64, cp: f64, cv: f64, beta: f64, mu: f64) -> f64 {
    let theta = t / REFERENCE_TEMPERATURE_K;
    let delta = rho / REFERENCE_DENSITY_KG_M3;
    let j = [0.310559006, 0.776397516, 1.242236025, 1.863354037]
        .partition_point(|&boundary| delta > boundary);
    let reference = polynomial(delta, &REFERENCE_RESPONSE[j]).recip();
    // zeta = d(rho/rho*)/d(p/p*) = p* delta beta_T. R15-11 p12,
    // footnote 2 caps these two auxiliary conductivity inputs at 1e13. It
    // does not replace the thermodynamic owner's original Cp or compressibility.
    let zeta = (REFERENCE_PRESSURE_PA * delta * beta).min(1e13);
    let cp_dimensionless = (cp / CONDUCTIVITY_GAS_CONSTANT_J_KG_K).min(1e13);
    let excess = delta * (zeta - reference * 1.5 / theta);
    // The zero-excess rule is explicitly part of source equation (23).
    if excess <= 0.0 {
        return 0.0;
    }
    let correlation_length_nm = 0.13 * det::pow(excess / 0.06, 0.630 / 1.239);
    let y = correlation_length_nm / 0.40;
    if y < 1.2e-7 {
        return 0.0;
    }
    let cv_over_cp = cv / cp;
    let decay = (y.recip() + (y / delta) * (y / delta) / 3.0).recip();
    let z = 2.0 / (core::f64::consts::PI * y)
        * ((1.0 - cv_over_cp) * det::atan(y) + cv_over_cp * y + det::expm1(-decay));
    1e-3 * 177.8514 * delta * cp_dimensionless * theta / (mu / 1e-6) * z
}

#[cfg(test)]
mod tests;
