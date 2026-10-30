//! Original IAPWS-IF97 liquid and steam Gibbs formulations in coherent SI.
//!
//! Implements Regions 1 and 2, the Region 4 saturation equations, and the B23
//! boundary from IAPWS R7-97(2012), equations (5)--(7), (15)--(17), (29)--(31).
//! Coefficients are transcribed from that release's Tables 1, 2, 10, 11 and 34.
//! Attribution: International Association for the Properties of Water and Steam,
//! <https://iapws.org/technical-guidance/release/IF97-Rev.download>.
//!
//! This is a separate source from the parent module's SR1-86 saturation model.
//! No coefficients, caloric datums or endpoints are mixed between the releases.
//! The separate `region3` module supplies the original dense-water Helmholtz
//! source. Region 5, ice and the separate metastable-vapor equation are not
//! implemented. Gibbs-region requests never extrapolate into those domains. IF97's
//! stated 273.15--273.16 K extension includes states metastable relative to ice.
//! Source uncertainties remain those in release section 12, not floating-point
//! error bars. The raw Gibbs energies of the two saturation branches have the
//! finite source consistency error described in section 10.2; no latent heat,
//! entropy, or phase reference is silently adjusted to force equality.

use super::{MassSpecificVolume, SaturationPressureDerivative};
use crate::{MassSpecificEnergyQuantityV1, MassSpecificThermalQuantityV1};
use core::fmt;
use fs_math::det;
use fs_qty::{Density, Pressure, Qty, Temperature, Velocity};

pub mod bounds;
pub mod increments;
pub mod region3;
/// Original-source regional routing, saturation pairs and boundary diagnostics.
pub mod regions;

/// Authoritative coefficient release and temperature scale.
pub const SOURCE_RELEASE: &str = "IAPWS R7-97(2012), IF97, ITS-90";
/// Official original equation, coefficient and verification-table document.
pub const SOURCE_URL: &str = "https://iapws.org/technical-guidance/release/IF97-Rev.download";
/// Deterministic polynomial/derivative and saturation-anchor operation tree.
pub const EVALUATOR_VERSION: u32 = 1;
/// Specific gas constant fixed by this source [J/(kg K)], equation (1).
pub const R_SPECIFIC_J_KG_K: f64 = 461.526;
/// Inclusive lower engineering temperature [K].
pub const MINIMUM_TEMPERATURE_K: f64 = 273.15;
/// Inclusive Region 1 upper temperature [K].
pub const REGION1_MAXIMUM_TEMPERATURE_K: f64 = 623.15;
/// Region 2 upper temperature [K]; hotter Region 5 is not implemented.
pub const REGION2_MAXIMUM_TEMPERATURE_K: f64 = 1073.15;
/// Source critical temperature [K]. The exact critical state has singular response.
pub const CRITICAL_TEMPERATURE_K: f64 = 647.096;
/// Source critical pressure [Pa], equation (3).
pub const CRITICAL_PRESSURE_PA: f64 = 22.064e6;
/// Upper pressure of both implemented regions [Pa].
pub const MAXIMUM_PRESSURE_PA: f64 = 100e6;
pub(super) const R1_PRESSURE_STAR_PA: f64 = 16.53e6;
pub(super) const R1_TEMPERATURE_STAR_K: f64 = 1386.0;
pub(super) const R2_PRESSURE_STAR_PA: f64 = 1e6;
pub(super) const R2_TEMPERATURE_STAR_K: f64 = 540.0;

/// Supported original Gibbs equation; a saturation point requires an explicit side.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum If97Region {
    /// Equation (7): compressed or saturated liquid through 623.15 K.
    Region1,
    /// Equations (15)--(17): vapor and the admitted high-temperature fluid.
    Region2,
}

/// Source branch and its saturation/supercritical status; no ice stability claim.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum If97Phase {
    /// Region 1 away from its saturation boundary.
    Liquid,
    /// Explicit Region 1 side of Region 4.
    SaturatedLiquid,
    /// Region 2 outside the saturated or supercritical classifications below.
    Vapor,
    /// Explicit Region 2 side of Region 4, through 623.15 K.
    SaturatedVapor,
    /// Region 2 with both T and p strictly above the critical constants.
    Supercritical,
}

/// Domain or finite-representation refusal from the original source evaluator.
#[derive(Clone, Debug, PartialEq)]
pub enum If97Error {
    /// Invalid/nonfinite input or a saturation query outside its source interval.
    InvalidInput(&'static str),
    /// Requested equation does not cover this point, including Region 3/5 points.
    OutsideRegion {
        region: If97Region,
        pressure_pa: f64,
        temperature_k: f64,
    },
    /// Saturation requires Region 3 above the implemented two-branch interval.
    UnsupportedSaturation { temperature_k: f64 },
    /// A required property or derivative cannot be represented as a finite value.
    Unrepresentable(&'static str),
}
impl fmt::Display for If97Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "IAPWS-IF97: {self:?}")
    }
}
impl std::error::Error for If97Error {}

/// Specific-volume derivative at fixed pressure [m3/(kg K)].
pub type SpecificVolumeTemperatureDerivative = Qty<3, -1, 0, -1, 0, 0>;
/// Specific-volume derivative at fixed temperature [m3/(kg Pa)].
pub type SpecificVolumePressureDerivative = Qty<4, -2, 2, 0, 0, 0>;
/// Volumetric thermal expansion coefficient [1/K]. Water may have negative expansion.
pub type ThermalExpansion = Qty<0, 0, 0, -1, 0, 0>;
/// Isothermal compressibility [1/Pa].
pub type IsothermalCompressibility = Qty<1, -1, 2, 0, 0, 0>;

/// Immutable admitted p/T source state, including original Gibbs derivatives.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct If97State {
    region: If97Region,
    phase: If97Phase,
    pressure_pa: f64,
    temperature_k: f64,
    volume_m3_kg: f64,
    enthalpy_j_kg: f64,
    internal_energy_j_kg: f64,
    entropy_j_kg_k: f64,
    gibbs_energy_j_kg: f64,
    cp_j_kg_k: f64,
    cv_j_kg_k: f64,
    volume_temperature_derivative: f64,
    volume_pressure_derivative: f64,
    sound_speed_m_s: f64,
}
impl If97State {
    /// Evaluate the explicitly requested original region without allocation.
    /// Both canonical Region 4 p(T) and T(p) anchors are recognized as saturated
    /// points; no tolerance band admits unrelated off-domain p/T combinations.
    pub fn at_pressure_temperature(
        pressure: Pressure,
        temperature: Temperature,
        region: If97Region,
    ) -> Result<Self, If97Error> {
        let (p, t) = (pressure.value(), temperature.value());
        admit(p, t, region)?;
        let saturated = t <= REGION1_MAXIMUM_TEMPERATURE_K && on_saturation(p, t);
        let phase = match (region, saturated) {
            (If97Region::Region1, true) => If97Phase::SaturatedLiquid,
            (If97Region::Region1, false) => If97Phase::Liquid,
            (If97Region::Region2, true) => If97Phase::SaturatedVapor,
            (If97Region::Region2, false)
                if t > CRITICAL_TEMPERATURE_K && p > CRITICAL_PRESSURE_PA =>
            {
                If97Phase::Supercritical
            }
            (If97Region::Region2, false) => If97Phase::Vapor,
        };
        evaluate(p, t, region, phase)
    }

    /// Actual original Gibbs equation used.
    pub const fn region(self) -> If97Region {
        self.region
    }
    /// Explicit phase/branch status of the admitted point.
    pub const fn phase(self) -> If97Phase {
        self.phase
    }
    /// Absolute source pressure [Pa].
    pub const fn pressure(self) -> Pressure {
        Pressure::new(self.pressure_pa)
    }
    /// Absolute ITS-90 source temperature [K].
    pub const fn temperature(self) -> Temperature {
        Temperature::new(self.temperature_k)
    }
    /// Specific volume [m3/kg].
    pub const fn specific_volume(self) -> MassSpecificVolume {
        MassSpecificVolume::new(self.volume_m3_kg)
    }
    /// Density [kg/m3], reciprocal of the original Gibbs specific volume.
    pub fn density(self) -> Density {
        Density::new(1.0 / self.volume_m3_kg)
    }
    /// Original specific enthalpy [J/kg], without reference reconciliation.
    pub const fn specific_enthalpy(self) -> MassSpecificEnergyQuantityV1 {
        MassSpecificEnergyQuantityV1::new(self.enthalpy_j_kg)
    }
    /// Original specific internal energy [J/kg].
    pub const fn specific_internal_energy(self) -> MassSpecificEnergyQuantityV1 {
        MassSpecificEnergyQuantityV1::new(self.internal_energy_j_kg)
    }
    /// Original specific entropy [J/(kg K)].
    pub const fn specific_entropy(self) -> MassSpecificThermalQuantityV1 {
        MassSpecificThermalQuantityV1::new(self.entropy_j_kg_k)
    }
    /// Original specific Gibbs energy evaluated directly as RT gamma [J/kg].
    pub const fn specific_gibbs_energy(self) -> MassSpecificEnergyQuantityV1 {
        MassSpecificEnergyQuantityV1::new(self.gibbs_energy_j_kg)
    }
    /// Isobaric heat capacity (dh/dT)_p [J/(kg K)].
    pub const fn cp(self) -> MassSpecificThermalQuantityV1 {
        MassSpecificThermalQuantityV1::new(self.cp_j_kg_k)
    }
    /// Isochoric heat capacity (du/dT)_v [J/(kg K)].
    pub const fn cv(self) -> MassSpecificThermalQuantityV1 {
        MassSpecificThermalQuantityV1::new(self.cv_j_kg_k)
    }
    /// Analytic (dv/dT)_p [m3/(kg K)].
    pub const fn dv_dt_at_pressure(self) -> SpecificVolumeTemperatureDerivative {
        SpecificVolumeTemperatureDerivative::new(self.volume_temperature_derivative)
    }
    /// Analytic (dv/dp)_T [m3/(kg Pa)], negative on these stable fluid branches.
    pub const fn dv_dp_at_temperature(self) -> SpecificVolumePressureDerivative {
        SpecificVolumePressureDerivative::new(self.volume_pressure_derivative)
    }
    /// Analytic (dh/dp)_T = v-T(dv/dT)_p [m3/kg].
    pub fn dh_dp_at_temperature(self) -> MassSpecificVolume {
        MassSpecificVolume::new(
            self.volume_m3_kg - self.temperature_k * self.volume_temperature_derivative,
        )
    }
    /// Analytic (ds/dp)_T = -(dv/dT)_p [m3/(kg K)].
    pub fn ds_dp_at_temperature(self) -> SpecificVolumeTemperatureDerivative {
        SpecificVolumeTemperatureDerivative::new(-self.volume_temperature_derivative)
    }
    /// Volumetric thermal expansion (dv/dT)_p/v [1/K].
    pub fn thermal_expansion(self) -> ThermalExpansion {
        ThermalExpansion::new(self.volume_temperature_derivative / self.volume_m3_kg)
    }
    /// Isothermal compressibility -(dv/dp)_T/v [1/Pa].
    pub fn isothermal_compressibility(self) -> IsothermalCompressibility {
        IsothermalCompressibility::new(-self.volume_pressure_derivative / self.volume_m3_kg)
    }
    /// Isentropic sound speed from the same Gibbs second derivatives [m/s].
    pub const fn sound_speed(self) -> Velocity {
        Velocity::new(self.sound_speed_m_s)
    }
}

/// Both original Gibbs branches at one canonical Region 4 saturation point.
/// Their raw latent h/s and Gibbs mismatch remain visible to the phase owner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct If97SaturationPair {
    liquid: If97State,
    vapor: If97State,
}
impl If97SaturationPair {
    /// Retain the requested pressure bit-for-bit and evaluate both branches at
    /// its one Region 4 temperature. Above p_s(623.15 K), Region 3 is required.
    pub fn at_pressure(pressure: Pressure) -> Result<Self, If97Error> {
        let t = saturation_temperature(pressure)?.value();
        // Admission uses the requested pressure as well as the inverse's
        // rounded temperature: one p ULP above this join can invert below it.
        if pressure.value() > saturation_pressure_raw(REGION1_MAXIMUM_TEMPERATURE_K) {
            return Err(If97Error::UnsupportedSaturation { temperature_k: t });
        }
        Self::evaluate(pressure.value(), t)
    }
    /// Retain the requested temperature and use its original Region 4 pressure.
    pub fn at_temperature(temperature: Temperature) -> Result<Self, If97Error> {
        let p = saturation_pressure(temperature)?.value();
        Self::evaluate(p, temperature.value())
    }
    fn evaluate(p: f64, t: f64) -> Result<Self, If97Error> {
        if t > REGION1_MAXIMUM_TEMPERATURE_K {
            return Err(If97Error::UnsupportedSaturation { temperature_k: t });
        }
        Ok(Self {
            liquid: evaluate(p, t, If97Region::Region1, If97Phase::SaturatedLiquid)?,
            vapor: evaluate(p, t, If97Region::Region2, If97Phase::SaturatedVapor)?,
        })
    }
    /// Original saturated liquid source state.
    pub const fn liquid(&self) -> &If97State {
        &self.liquid
    }
    /// Original saturated vapor source state.
    pub const fn vapor(&self) -> &If97State {
        &self.vapor
    }
    /// Shared canonical saturation pressure [Pa].
    pub const fn pressure(self) -> Pressure {
        self.liquid.pressure()
    }
    /// Shared canonical saturation temperature [K].
    pub const fn temperature(self) -> Temperature {
        self.liquid.temperature()
    }
    /// Raw vapor-minus-liquid enthalpy [J/kg], not replaced by T delta_s.
    pub fn latent_enthalpy(self) -> MassSpecificEnergyQuantityV1 {
        MassSpecificEnergyQuantityV1::new(self.vapor.enthalpy_j_kg - self.liquid.enthalpy_j_kg)
    }
    /// Raw vapor-minus-liquid entropy [J/(kg K)].
    pub fn latent_entropy(self) -> MassSpecificThermalQuantityV1 {
        MassSpecificThermalQuantityV1::new(self.vapor.entropy_j_kg_k - self.liquid.entropy_j_kg_k)
    }
    /// Original finite Gibbs inconsistency g_v-g_l [J/kg]; never silently zeroed.
    pub fn gibbs_energy_gap(self) -> MassSpecificEnergyQuantityV1 {
        MassSpecificEnergyQuantityV1::new(
            self.vapor.gibbs_energy_j_kg - self.liquid.gibbs_energy_j_kg,
        )
    }
}

/// Original Region 4 pressure, equation (30), over 273.15..=647.096 K.
/// The rounded released coefficients are retained, including endpoint roundoff.
pub fn saturation_pressure(temperature: Temperature) -> Result<Pressure, If97Error> {
    let t = temperature.value();
    if !t.is_finite() || !(MINIMUM_TEMPERATURE_K..=CRITICAL_TEMPERATURE_K).contains(&t) {
        return Err(If97Error::InvalidInput(
            "Region 4 temperature must be in 273.15..=647.096 K",
        ));
    }
    let p = saturation_pressure_raw(t);
    if !p.is_finite() || p <= 0.0 {
        return Err(If97Error::Unrepresentable("saturation pressure"));
    }
    Ok(Pressure::new(p))
}

/// Original Region 4 analytic inverse, equation (31), without iteration.
/// Pressure limits are evaluated equation-(30) endpoints. Those endpoints and
/// the 623.15 K join map exactly back to their known source temperatures, so a
/// forward/inverse roundtrip cannot move a closed source endpoint out of range.
pub fn saturation_temperature(pressure: Pressure) -> Result<Temperature, If97Error> {
    let p = pressure.value();
    let minimum = saturation_pressure_raw(MINIMUM_TEMPERATURE_K);
    let maximum = saturation_pressure_raw(CRITICAL_TEMPERATURE_K);
    if !p.is_finite() || p < minimum || p > maximum {
        return Err(If97Error::InvalidInput(
            "pressure is outside the evaluated Region 4 endpoints",
        ));
    }
    let t = if p == minimum {
        MINIMUM_TEMPERATURE_K
    } else if p == maximum {
        CRITICAL_TEMPERATURE_K
    } else if p == saturation_pressure_raw(REGION1_MAXIMUM_TEMPERATURE_K) {
        REGION1_MAXIMUM_TEMPERATURE_K
    } else {
        saturation_temperature_raw(p)
    };
    if !t.is_finite() || !(MINIMUM_TEMPERATURE_K..=CRITICAL_TEMPERATURE_K).contains(&t) {
        return Err(If97Error::Unrepresentable("saturation temperature"));
    }
    Ok(Temperature::new(t))
}

/// Analytic Region 4 pressure slope [Pa/K], from implicit equation (29).
pub fn saturation_pressure_derivative(
    temperature: Temperature,
) -> Result<SaturationPressureDerivative, If97Error> {
    let p = saturation_pressure(temperature)?.value();
    let t = temperature.value();
    let theta = t + R4[8] / (t - R4[9]);
    let beta = det::sqrt(det::sqrt(p / 1e6));
    let a = theta * theta + R4[0] * theta + R4[1];
    let b = R4[2] * theta * theta + R4[3] * theta + R4[4];
    let numerator = (2.0 * theta + R4[0]) * beta * beta
        + (2.0 * R4[2] * theta + R4[3]) * beta
        + 2.0 * R4[5] * theta
        + R4[6];
    let slope = -4e6 * beta * beta * beta * numerator / (2.0 * a * beta + b)
        * (1.0 - R4[8] / ((t - R4[9]) * (t - R4[9])));
    if !slope.is_finite() || slope <= 0.0 {
        return Err(If97Error::Unrepresentable("saturation pressure slope"));
    }
    Ok(SaturationPressureDerivative::new(slope))
}

/// B23 dividing pressure from equation (5), only for 623.15..=863.15 K.
pub fn boundary23_pressure(temperature: Temperature) -> Result<Pressure, If97Error> {
    let t = temperature.value();
    if !t.is_finite() || !(623.15..=863.15).contains(&t) {
        return Err(If97Error::InvalidInput(
            "B23 temperature must be in 623.15..=863.15 K",
        ));
    }
    Ok(Pressure::new(boundary23_pressure_raw(t)))
}
pub(super) fn boundary23_pressure_raw(t: f64) -> f64 {
    1e6 * (B23_COEFFICIENTS[0] - B23_COEFFICIENTS[1] * t + B23_COEFFICIENTS[2] * t * t)
}

fn on_saturation(p: f64, t: f64) -> bool {
    p == saturation_pressure_raw(t)
        || saturation_temperature(Pressure::new(p)).is_ok_and(|source| source.value() == t)
}
pub(super) fn admit(p: f64, t: f64, region: If97Region) -> Result<(), If97Error> {
    if !p.is_finite() || p <= 0.0 || !t.is_finite() {
        return Err(If97Error::InvalidInput(
            "positive finite pressure and finite temperature required",
        ));
    }
    let inside = p <= MAXIMUM_PRESSURE_PA
        && match region {
            If97Region::Region1 => {
                (MINIMUM_TEMPERATURE_K..=REGION1_MAXIMUM_TEMPERATURE_K).contains(&t)
                    && (p >= saturation_pressure_raw(t) || on_saturation(p, t))
            }
            If97Region::Region2 => {
                (MINIMUM_TEMPERATURE_K..=REGION2_MAXIMUM_TEMPERATURE_K).contains(&t)
                    && if t <= REGION1_MAXIMUM_TEMPERATURE_K {
                        p <= saturation_pressure_raw(t) || on_saturation(p, t)
                    } else if t <= 863.15 {
                        p <= boundary23_pressure_raw(t)
                    } else {
                        true
                    }
            }
        };
    if inside {
        Ok(())
    } else {
        Err(If97Error::OutsideRegion {
            region,
            pressure_pa: p,
            temperature_k: t,
        })
    }
}

fn evaluate(p: f64, t: f64, region: If97Region, phase: If97Phase) -> Result<If97State, If97Error> {
    let (tau, gamma) = match region {
        If97Region::Region1 => {
            let tau = R1_TEMPERATURE_STAR_K / t;
            (tau, region1(p / R1_PRESSURE_STAR_PA, tau))
        }
        If97Region::Region2 => {
            let tau = R2_TEMPERATURE_STAR_K / t;
            (tau, region2(p / R2_PRESSURE_STAR_PA, tau))
        }
    };
    let rt = R_SPECIFIC_J_KG_K * t;
    let volume = (rt / p) * gamma.p;
    let h = rt * tau * gamma.t;
    let u = rt * (tau * gamma.t - gamma.p);
    let s = R_SPECIFIC_J_KG_K * (tau * gamma.t - gamma.value);
    let cp_r = -tau * tau * gamma.tt;
    let coupling = gamma.p - tau * gamma.pt;
    let cv_r = cp_r + coupling * coupling / gamma.pp;
    let dv_dt = (R_SPECIFIC_J_KG_K / p) * coupling;
    let dv_dp = ((rt / p) / p) * gamma.pp;
    let sound2 = rt * gamma.p * gamma.p / (-coupling * coupling / cp_r - gamma.pp);
    let (cp, cv, gibbs) = (
        R_SPECIFIC_J_KG_K * cp_r,
        R_SPECIFIC_J_KG_K * cv_r,
        rt * gamma.value,
    );
    if ![
        volume,
        h,
        u,
        s,
        cp,
        cv,
        dv_dt,
        dv_dp,
        sound2,
        gibbs,
        1.0 / volume,
        dv_dt / volume,
        -dv_dp / volume,
        volume - t * dv_dt,
    ]
    .iter()
    .all(|x| x.is_finite())
        || volume <= 0.0
        || cp <= 0.0
        || cv <= 0.0
        || dv_dp >= 0.0
        || sound2 <= 0.0
    {
        return Err(If97Error::Unrepresentable(
            "Gibbs properties or stable second derivatives",
        ));
    }
    Ok(If97State {
        region,
        phase,
        pressure_pa: p,
        temperature_k: t,
        volume_m3_kg: volume,
        enthalpy_j_kg: h,
        internal_energy_j_kg: u,
        entropy_j_kg_k: s,
        gibbs_energy_j_kg: gibbs,
        cp_j_kg_k: cp,
        cv_j_kg_k: cv,
        volume_temperature_derivative: dv_dt,
        volume_pressure_derivative: dv_dp,
        sound_speed_m_s: det::sqrt(sound2),
    })
}

// Pressure derivatives are scaled as pi*gamma_pi, pi^2*gamma_pipi and
// pi*gamma_pitau. This retains the Region 2 ideal-gas limit without forming
// enormous 1/pi^2 intermediates at low pressure.
struct Gibbs {
    value: f64,
    t: f64,
    tt: f64,
    p: f64,
    pp: f64,
    pt: f64,
}
#[derive(Default, Clone, Copy)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, x: f64) {
        let x = x - self.correction;
        let next = self.value + x;
        self.correction = (next - self.value) - x;
        self.value = next;
    }
}
fn powers<const N: usize>(x: f64, zero: usize) -> [f64; N] {
    let mut result = [1.0; N];
    for i in zero + 1..N {
        result[i] = result[i - 1] * x;
    }
    let inverse = 1.0 / x;
    for i in (0..zero).rev() {
        result[i] = result[i + 1] * inverse;
    }
    result
}
fn region1(pi: f64, tau: f64) -> Gibbs {
    let x = powers::<33>(7.1 - pi, 0);
    let y = powers::<61>(tau - 1.222, 43);
    let mut sums = [Sum::default(); 6];
    for &GibbsTerm { i, j, n } in &R1_TERMS {
        let (ix, jy) = (i as usize, (j + 43) as usize);
        let ni = n * f64::from(i);
        let nj = n * f64::from(j);
        sums[0].add(n * x[ix] * y[jy]);
        if j != 0 {
            sums[1].add(nj * x[ix] * y[jy - 1]);
        }
        if j != 0 && j != 1 {
            sums[2].add(nj * f64::from(j - 1) * x[ix] * y[jy - 2]);
        }
        if i != 0 {
            sums[3].add(-ni * x[ix - 1] * y[jy]);
            if j != 0 {
                sums[5].add(-ni * f64::from(j) * x[ix - 1] * y[jy - 1]);
            }
        }
        if i >= 2 {
            sums[4].add(ni * f64::from(i - 1) * x[ix - 2] * y[jy]);
        }
    }
    Gibbs {
        value: sums[0].value,
        t: sums[1].value,
        tt: sums[2].value,
        p: pi * sums[3].value,
        pp: pi * pi * sums[4].value,
        pt: pi * sums[5].value,
    }
}
fn region2(pi: f64, tau: f64) -> Gibbs {
    let x = powers::<25>(pi, 0);
    let y = powers::<59>(tau - 0.5, 0);
    let ideal = powers::<11>(tau, 7);
    let mut sums = [Sum::default(); 6];
    sums[0].add(det::ln(pi));
    sums[3].add(1.0);
    sums[4].add(-1.0);
    for &IdealTerm { j, n } in &R2_IDEAL {
        let k = (j + 7) as usize;
        sums[0].add(n * ideal[k]);
        if j != 0 {
            sums[1].add(n * f64::from(j) * ideal[k - 1]);
        }
        if j != 0 && j != 1 {
            sums[2].add(n * f64::from(j * (j - 1)) * ideal[k - 2]);
        }
    }
    for &GibbsTerm { i, j, n } in &R2_RESIDUAL {
        let (ix, jy) = (i as usize, j as usize);
        let coefficient = n * x[ix];
        let term = coefficient * y[jy];
        sums[0].add(term);
        sums[3].add(f64::from(i) * term);
        sums[4].add(f64::from(i * (i - 1)) * term);
        if j != 0 {
            let derivative = coefficient * f64::from(j) * y[jy - 1];
            sums[1].add(derivative);
            sums[5].add(f64::from(i) * derivative);
        }
        if j >= 2 {
            sums[2].add(coefficient * f64::from(j * (j - 1)) * y[jy - 2]);
        }
    }
    Gibbs {
        value: sums[0].value,
        t: sums[1].value,
        tt: sums[2].value,
        p: sums[3].value,
        pp: sums[4].value,
        pt: sums[5].value,
    }
}

fn saturation_pressure_raw(t: f64) -> f64 {
    let theta = t + R4[8] / (t - R4[9]);
    let a = theta * theta + R4[0] * theta + R4[1];
    let b = R4[2] * theta * theta + R4[3] * theta + R4[4];
    let c = R4[5] * theta * theta + R4[6] * theta + R4[7];
    let beta = 2.0 * c / (-b + det::sqrt(b * b - 4.0 * a * c));
    let square = beta * beta;
    1e6 * square * square
}
fn saturation_temperature_raw(p: f64) -> f64 {
    let beta = det::sqrt(det::sqrt(p / 1e6));
    let e = beta * beta + R4[2] * beta + R4[5];
    let f = R4[0] * beta * beta + R4[3] * beta + R4[6];
    let g = R4[1] * beta * beta + R4[4] * beta + R4[7];
    let d = 2.0 * g / (-f - det::sqrt(f * f - 4.0 * e * g));
    // Algebraic cancellation of 4*n10*d before rounding the discriminant.
    0.5 * (R4[9] + d - det::sqrt((R4[9] - d) * (R4[9] - d) - 4.0 * R4[8]))
}

#[derive(Clone, Copy)]
pub(super) struct GibbsTerm {
    pub(super) i: i32,
    pub(super) j: i32,
    pub(super) n: f64,
}
#[derive(Clone, Copy)]
pub(super) struct IdealTerm {
    pub(super) j: i32,
    pub(super) n: f64,
}
const fn g(i: i32, j: i32, n: f64) -> GibbsTerm {
    GibbsTerm { i, j, n }
}
const fn ideal(j: i32, n: f64) -> IdealTerm {
    IdealTerm { j, n }
}

// Table 2, equation (7). Every published coefficient digit is retained.
pub(super) const R1_TERMS: [GibbsTerm; 34] = [
    g(0, -2, 0.14632971213167),
    g(0, -1, -0.84548187169114),
    g(0, 0, -0.37563603672040e1),
    g(0, 1, 0.33855169168385e1),
    g(0, 2, -0.95791963387872),
    g(0, 3, 0.15772038513228),
    g(0, 4, -0.16616417199501e-1),
    g(0, 5, 0.81214629983568e-3),
    g(1, -9, 0.28319080123804e-3),
    g(1, -7, -0.60706301565874e-3),
    g(1, -1, -0.18990068218419e-1),
    g(1, 0, -0.32529748770505e-1),
    g(1, 1, -0.21841717175414e-1),
    g(1, 3, -0.52838357969930e-4),
    g(2, -3, -0.47184321073267e-3),
    g(2, 0, -0.30001780793026e-3),
    g(2, 1, 0.47661393906987e-4),
    g(2, 3, -0.44141845330846e-5),
    g(2, 17, -0.72694996297594e-15),
    g(3, -4, -0.31679644845054e-4),
    g(3, 0, -0.28270797985312e-5),
    g(3, 6, -0.85205128120103e-9),
    g(4, -5, -0.22425281908000e-5),
    g(4, -2, -0.65171222895601e-6),
    g(4, 10, -0.14341729937924e-12),
    g(5, -8, -0.40516996860117e-6),
    g(8, -11, -0.12734301741641e-8),
    g(8, -6, -0.17424871230634e-9),
    g(21, -29, -0.68762131295531e-18),
    g(23, -31, 0.14478307828521e-19),
    g(29, -38, 0.26335781662795e-22),
    g(30, -39, -0.11947622640071e-22),
    g(31, -40, 0.18228094581404e-23),
    g(32, -41, -0.93537087292458e-25),
];
// Table 10, stable-vapor equation (16); NOT the metastable substitutions.
pub(super) const R2_IDEAL: [IdealTerm; 9] = [
    ideal(0, -0.96927686500217e1),
    ideal(1, 0.10086655968018e2),
    ideal(-5, -0.56087911283020e-2),
    ideal(-4, 0.71452738081455e-1),
    ideal(-3, -0.40710498223928),
    ideal(-2, 0.14240819171444e1),
    ideal(-1, -0.43839511319450e1),
    ideal(2, -0.28408632460772),
    ideal(3, 0.21268463753307e-1),
];
// Table 11, stable-vapor residual equation (17).
pub(super) const R2_RESIDUAL: [GibbsTerm; 43] = [
    g(1, 0, -0.17731742473213e-2),
    g(1, 1, -0.17834862292358e-1),
    g(1, 2, -0.45996013696365e-1),
    g(1, 3, -0.57581259083432e-1),
    g(1, 6, -0.50325278727930e-1),
    g(2, 1, -0.33032641670203e-4),
    g(2, 2, -0.18948987516315e-3),
    g(2, 4, -0.39392777243355e-2),
    g(2, 7, -0.43797295650573e-1),
    g(2, 36, -0.26674547914087e-4),
    g(3, 0, 0.20481737692309e-7),
    g(3, 1, 0.43870667284435e-6),
    g(3, 3, -0.32277677238570e-4),
    g(3, 6, -0.15033924542148e-2),
    g(3, 35, -0.40668253562649e-1),
    g(4, 1, -0.78847309559367e-9),
    g(4, 2, 0.12790717852285e-7),
    g(4, 3, 0.48225372718507e-6),
    g(5, 7, 0.22922076337661e-5),
    g(6, 3, -0.16714766451061e-10),
    g(6, 16, -0.21171472321355e-2),
    g(6, 35, -0.23895741934104e2),
    g(7, 0, -0.59059564324270e-17),
    g(7, 11, -0.12621808899101e-5),
    g(7, 25, -0.38946842435739e-1),
    g(8, 8, 0.11256211360459e-10),
    g(8, 36, -0.82311340897998e1),
    g(9, 13, 0.19809712802088e-7),
    g(10, 4, 0.10406965210174e-18),
    g(10, 10, -0.10234747095929e-12),
    g(10, 14, -0.10018179379511e-8),
    g(16, 29, -0.80882908646985e-10),
    g(16, 50, 0.10693031879409),
    g(18, 57, -0.33662250574171),
    g(20, 20, 0.89185845355421e-24),
    g(20, 35, 0.30629316876232e-12),
    g(20, 48, -0.42002467698208e-5),
    g(21, 21, -0.59056029685639e-25),
    g(22, 53, 0.37826947613457e-5),
    g(23, 39, -0.12768608934681e-14),
    g(24, 26, 0.73087610595061e-28),
    g(24, 40, 0.55414715350778e-16),
    g(24, 58, -0.94369707241210e-6),
];
// Table 1, equation (5), shared with original-source interval admission.
pub(super) const B23_COEFFICIENTS: [f64; 3] =
    [348.05185628969, 1.1671859879975, 0.0010192970039326];

// Table 34, equations (29)--(31).
const R4: [f64; 10] = [
    0.11670521452767e4,
    -0.72421316703206e6,
    -0.17073846940092e2,
    0.12020824702470e5,
    -0.32325550322333e7,
    0.14915108613530e2,
    -0.48232657361591e4,
    0.40511340542057e6,
    -0.23855557567849,
    0.65017534844798e3,
];

#[cfg(test)]
mod tests;
