//! Original IF97 fluid routing and explicit regional matching defects.
//!
//! Regions 1, 2 and 3 retain their own equations, references and admission.
//! Routing supplies no interpolated EOS, latent-heat adjustment or continuity
//! correction. In particular, a single-phase regional boundary is not a phase
//! transition. Its original h/u/s/v mismatch remains visible to the consumer.
//! IAPWS R7-97(2012), sections 2, 4 and 10 describe these regional boundaries
//! and the finite matching errors: <https://iapws.org/technical-guidance/release/IF97-Rev.download>.
//!
//! The stable engineering ingress starts at the canonical 273.16 K anchor;
//! ice, metastable continuations, Region 5 and the exact critical singularity
//! remain unsupported. Near-critical source solves may also refuse resolution.

use super::region3::pressure::{
    DensePressureError, DensePressureSolution, DensePressureSolve, If97DenseSaturationPair,
};
use super::region3::{If97DensePhase, If97DenseState};
use super::{
    CRITICAL_PRESSURE_PA, CRITICAL_TEMPERATURE_K, If97Error, If97Phase, If97Region,
    If97SaturationPair, If97State, IsothermalCompressibility, MAXIMUM_PRESSURE_PA,
    REGION1_MAXIMUM_TEMPERATURE_K, REGION2_MAXIMUM_TEMPERATURE_K, SpecificVolumePressureDerivative,
    SpecificVolumeTemperatureDerivative, ThermalExpansion, boundary23_pressure_raw, on_saturation,
    saturation_pressure,
};
use crate::water::MassSpecificVolume;
use crate::{MassSpecificEnergyQuantityV1, MassSpecificThermalQuantityV1};
use core::fmt;
use fs_qty::{Density, Pressure, Temperature, Velocity};

pub const MINIMUM_STABLE_TEMPERATURE_K: f64 = 273.16;

/// Original source equation, independent of liquid/vapor phase identity.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum If97FluidRegion {
    Region1,
    Region2,
    Region3,
}

/// Original phase classification; a saturation pair additionally names its sides.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum If97FluidPhase {
    Liquid,
    SaturatedLiquid,
    Vapor,
    SaturatedVapor,
    Supercritical,
    CriticalIsotherm,
}

#[derive(Clone, Debug, PartialEq)]
pub enum If97FluidError {
    Invalid(&'static str),
    OutsideDomain,
    Gibbs(If97Error),
    Dense(DensePressureError),
    AmbiguousSaturation,
    CriticalSingularity,
    BoundaryNotConverged {
        iterations: usize,
        temperature_bounds_k: [f64; 2],
    },
    Cancelled,
}
impl From<If97Error> for If97FluidError {
    fn from(error: If97Error) -> Self {
        Self::Gibbs(error)
    }
}
impl From<DensePressureError> for If97FluidError {
    fn from(error: DensePressureError) -> Self {
        match error {
            DensePressureError::Cancelled => Self::Cancelled,
            DensePressureError::AmbiguousSaturation => Self::AmbiguousSaturation,
            DensePressureError::CriticalSingularity => Self::CriticalSingularity,
            other => Self::Dense(other),
        }
    }
}
impl fmt::Display for If97FluidError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "original IF97 fluid regions: {self:?}")
    }
}
impl std::error::Error for If97FluidError {}

/// Unchanged source state plus the dense solver's requested-pressure diagnostics.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum If97FluidState {
    Gibbs(If97State),
    Dense(DensePressureSolution),
}

macro_rules! source_property {
    ($name:ident, $quantity:ty) => {
        pub fn $name(self) -> $quantity {
            match self {
                Self::Gibbs(s) => s.$name(),
                Self::Dense(s) => s.state.$name(),
            }
        }
    };
}

impl If97FluidState {
    /// Select the original stable source. Canonical p(T) AND T(p) saturation
    /// anchors require an explicit pair; neither side is chosen by roundoff.
    pub fn at_pressure_temperature(
        pressure: Pressure,
        temperature: Temperature,
        controls: DensePressureSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, If97FluidError> {
        poll(&mut cancelled)?;
        validate_point(pressure.value(), temperature.value())?;
        let (p, t) = (pressure.value(), temperature.value());
        if t <= CRITICAL_TEMPERATURE_K && on_saturation(p, t) {
            return Err(If97FluidError::AmbiguousSaturation);
        }
        let region = if t <= REGION1_MAXIMUM_TEMPERATURE_K {
            if p > saturation_pressure(temperature)?.value() {
                If97FluidRegion::Region1
            } else {
                If97FluidRegion::Region2
            }
        } else if t <= 863.15 && p > boundary23_pressure_raw(t) {
            If97FluidRegion::Region3
        } else {
            If97FluidRegion::Region2
        };
        Self::at_pressure_temperature_in_region(pressure, temperature, region, controls, cancelled)
    }

    /// Explicit source selection for prepared sensible branches and boundary
    /// endpoints. Each original equation still enforces its own closed domain.
    /// Region 1/2 saturation sides are explicit; Region 3 needs a paired solve.
    pub fn at_pressure_temperature_in_region(
        pressure: Pressure,
        temperature: Temperature,
        region: If97FluidRegion,
        controls: DensePressureSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, If97FluidError> {
        poll(&mut cancelled)?;
        validate_point(pressure.value(), temperature.value())?;
        let state = match region {
            If97FluidRegion::Region1 | If97FluidRegion::Region2 => {
                let source = if region == If97FluidRegion::Region1 {
                    If97Region::Region1
                } else {
                    If97Region::Region2
                };
                Self::Gibbs(If97State::at_pressure_temperature(
                    pressure,
                    temperature,
                    source,
                )?)
            }
            If97FluidRegion::Region3 => {
                if temperature.value() < CRITICAL_TEMPERATURE_K
                    && on_saturation(pressure.value(), temperature.value())
                {
                    return Err(If97FluidError::AmbiguousSaturation);
                }
                Self::Dense(If97DenseState::at_pressure_temperature(
                    pressure,
                    temperature,
                    controls,
                    &mut cancelled,
                )?)
            }
        };
        poll(&mut cancelled)?;
        Ok(state)
    }

    pub fn region(self) -> If97FluidRegion {
        match self {
            Self::Gibbs(s) => match s.region() {
                If97Region::Region1 => If97FluidRegion::Region1,
                If97Region::Region2 => If97FluidRegion::Region2,
            },
            Self::Dense(_) => If97FluidRegion::Region3,
        }
    }
    pub fn phase(self) -> If97FluidPhase {
        match self {
            Self::Gibbs(s) => match s.phase() {
                If97Phase::Liquid => If97FluidPhase::Liquid,
                If97Phase::SaturatedLiquid => If97FluidPhase::SaturatedLiquid,
                If97Phase::Vapor => If97FluidPhase::Vapor,
                If97Phase::SaturatedVapor => If97FluidPhase::SaturatedVapor,
                If97Phase::Supercritical => If97FluidPhase::Supercritical,
            },
            Self::Dense(s) => match s.state.phase() {
                If97DensePhase::Liquid => If97FluidPhase::Liquid,
                If97DensePhase::Vapor => If97FluidPhase::Vapor,
                If97DensePhase::Supercritical => If97FluidPhase::Supercritical,
                If97DensePhase::CriticalIsotherm => If97FluidPhase::CriticalIsotherm,
            },
        }
    }
    pub fn requested_pressure(self) -> Pressure {
        match self {
            Self::Gibbs(s) => s.pressure(),
            Self::Dense(s) => Pressure::new(s.requested_pressure_pa),
        }
    }
    pub const fn pressure_residual_pa(self) -> f64 {
        match self {
            Self::Gibbs(_) => 0.0,
            Self::Dense(s) => s.pressure_residual_pa,
        }
    }
    /// Dense branch/root iterations; Gibbs evaluation performs no density solve.
    pub const fn iterations(self) -> usize {
        match self {
            Self::Gibbs(_) => 0,
            Self::Dense(s) => s.iterations,
        }
    }
    /// Original dense pressure/stiffness queries; zero for the Gibbs sources.
    pub const fn source_evaluations(self) -> usize {
        match self {
            Self::Gibbs(_) => 0,
            Self::Dense(s) => s.evaluations,
        }
    }

    source_property!(pressure, Pressure);
    source_property!(temperature, Temperature);
    source_property!(density, Density);
    source_property!(specific_volume, MassSpecificVolume);
    source_property!(specific_enthalpy, MassSpecificEnergyQuantityV1);
    source_property!(specific_internal_energy, MassSpecificEnergyQuantityV1);
    source_property!(specific_entropy, MassSpecificThermalQuantityV1);
    source_property!(specific_gibbs_energy, MassSpecificEnergyQuantityV1);
    source_property!(cp, MassSpecificThermalQuantityV1);
    source_property!(cv, MassSpecificThermalQuantityV1);
    source_property!(sound_speed, Velocity);
    source_property!(dv_dt_at_pressure, SpecificVolumeTemperatureDerivative);
    source_property!(dv_dp_at_temperature, SpecificVolumePressureDerivative);
    source_property!(thermal_expansion, ThermalExpansion);
    source_property!(isothermal_compressibility, IsothermalCompressibility);

    pub fn dh_dp_at_temperature(self) -> MassSpecificVolume {
        MassSpecificVolume::new(
            self.specific_volume().value()
                - self.temperature().value() * self.dv_dt_at_pressure().value(),
        )
    }
    pub fn ds_dp_at_temperature(self) -> SpecificVolumeTemperatureDerivative {
        SpecificVolumeTemperatureDerivative::new(-self.dv_dt_at_pressure().value())
    }
}

/// Original saturation endpoints; there is no bulk-quality or mixture EOS here.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum If97FluidSaturationPair {
    Gibbs(If97SaturationPair),
    Dense(If97DenseSaturationPair),
}
impl If97FluidSaturationPair {
    pub fn at_temperature(
        temperature: Temperature,
        controls: DensePressureSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, If97FluidError> {
        poll(&mut cancelled)?;
        let t = temperature.value();
        if !t.is_finite() || !(MINIMUM_STABLE_TEMPERATURE_K..CRITICAL_TEMPERATURE_K).contains(&t) {
            return Err(If97FluidError::OutsideDomain);
        }
        let pair = if t <= REGION1_MAXIMUM_TEMPERATURE_K {
            Self::Gibbs(If97SaturationPair::at_temperature(temperature)?)
        } else {
            Self::Dense(If97DenseSaturationPair::at_temperature(
                temperature,
                controls,
                &mut cancelled,
            )?)
        };
        poll(&mut cancelled)?;
        Ok(pair)
    }
    pub fn at_pressure(
        pressure: Pressure,
        controls: DensePressureSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, If97FluidError> {
        poll(&mut cancelled)?;
        let p = pressure.value();
        let minimum = saturation_pressure(Temperature::new(MINIMUM_STABLE_TEMPERATURE_K))?.value();
        if !p.is_finite() || p < minimum || p >= CRITICAL_PRESSURE_PA {
            return Err(If97FluidError::OutsideDomain);
        }
        // Select by the requested pressure: the rounded inverse can land on the
        // other side of 623.15 K for the immediately neighboring p value.
        let join = saturation_pressure(Temperature::new(REGION1_MAXIMUM_TEMPERATURE_K))?.value();
        // Preserve exact forward anchors through their original constructor.
        // The independently fitted T(p) equation can otherwise round below the
        // minimum temperature or onto the wrong side of the regional join.
        let pair = if p == minimum {
            Self::Gibbs(If97SaturationPair::at_temperature(Temperature::new(
                MINIMUM_STABLE_TEMPERATURE_K,
            ))?)
        } else if p == join {
            Self::Gibbs(If97SaturationPair::at_temperature(Temperature::new(
                REGION1_MAXIMUM_TEMPERATURE_K,
            ))?)
        } else if p < join {
            Self::Gibbs(If97SaturationPair::at_pressure(pressure)?)
        } else {
            Self::Dense(If97DenseSaturationPair::at_pressure(
                pressure,
                controls,
                &mut cancelled,
            )?)
        };
        if pair.temperature().value() < MINIMUM_STABLE_TEMPERATURE_K {
            return Err(If97FluidError::OutsideDomain);
        }
        poll(&mut cancelled)?;
        Ok(pair)
    }
    pub fn liquid(self) -> If97FluidState {
        match self {
            Self::Gibbs(pair) => If97FluidState::Gibbs(*pair.liquid()),
            Self::Dense(pair) => If97FluidState::Dense(pair.liquid_solution()),
        }
    }
    pub fn vapor(self) -> If97FluidState {
        match self {
            Self::Gibbs(pair) => If97FluidState::Gibbs(*pair.vapor()),
            Self::Dense(pair) => If97FluidState::Dense(pair.vapor_solution()),
        }
    }
    pub fn pressure(self) -> Pressure {
        match self {
            Self::Gibbs(pair) => pair.pressure(),
            Self::Dense(pair) => pair.pressure(),
        }
    }
    pub fn temperature(self) -> Temperature {
        match self {
            Self::Gibbs(pair) => pair.temperature(),
            Self::Dense(pair) => pair.temperature(),
        }
    }
    pub fn saturation_pressure_residual_pa(self) -> f64 {
        match self {
            Self::Gibbs(pair) => {
                super::saturation_pressure_raw(pair.temperature().value()) - pair.pressure().value()
            }
            Self::Dense(pair) => pair.saturation_pressure_residual_pa(),
        }
    }
    pub fn latent_enthalpy(self) -> MassSpecificEnergyQuantityV1 {
        match self {
            Self::Gibbs(pair) => pair.latent_enthalpy(),
            Self::Dense(pair) => pair.latent_enthalpy(),
        }
    }
    pub fn latent_entropy(self) -> MassSpecificThermalQuantityV1 {
        match self {
            Self::Gibbs(pair) => pair.latent_entropy(),
            Self::Dense(pair) => pair.latent_entropy(),
        }
    }
    pub fn gibbs_energy_gap(self) -> MassSpecificEnergyQuantityV1 {
        match self {
            Self::Gibbs(pair) => pair.gibbs_energy_gap(),
            Self::Dense(pair) => pair.gibbs_energy_gap(),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum If97JunctionKind {
    LiquidDense,
    DenseVapor,
}

/// Separate bounded work for the numerical B23 temperature bracket and density.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct If97JunctionSolve {
    pub pressure: DensePressureSolve,
    pub temperature_tolerance_k: f64,
    pub maximum_iterations: usize,
}

/// Two admitted source endpoints in increasing-temperature regional order.
///
/// B23 endpoints can have slightly different temperatures because the original
/// forward polynomial does not necessarily represent the requested pressure at
/// a single binary64 temperature. The width is reported and bounded; this is a
/// numerical source-admission bracket, not a certified real-arithmetic root.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct If97RegionalJunction {
    kind: If97JunctionKind,
    lower: If97FluidState,
    upper: If97FluidState,
    boundary_iterations: usize,
}
impl If97RegionalJunction {
    pub fn new(
        pressure: Pressure,
        kind: If97JunctionKind,
        controls: If97JunctionSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, If97FluidError> {
        poll(&mut cancelled)?;
        if !controls.temperature_tolerance_k.is_finite()
            || controls.temperature_tolerance_k <= 0.0
            || !(1..=128).contains(&controls.maximum_iterations)
        {
            return Err(If97FluidError::Invalid(
                "positive finite boundary temperature width and 1..=128 iterations required",
            ));
        }
        let p = pressure.value();
        let join = saturation_pressure(Temperature::new(REGION1_MAXIMUM_TEMPERATURE_K))?.value();
        if !p.is_finite()
            || p <= join.max(boundary23_pressure_raw(REGION1_MAXIMUM_TEMPERATURE_K))
            || p > MAXIMUM_PRESSURE_PA
        {
            return Err(If97FluidError::OutsideDomain);
        }
        let (temperatures, boundary_iterations, lower_region, upper_region) = match kind {
            If97JunctionKind::LiquidDense => (
                [REGION1_MAXIMUM_TEMPERATURE_K; 2],
                0,
                If97FluidRegion::Region1,
                If97FluidRegion::Region3,
            ),
            If97JunctionKind::DenseVapor => {
                let (bounds, iterations) = boundary_temperatures(p, controls, &mut cancelled)?;
                (
                    bounds,
                    iterations,
                    If97FluidRegion::Region3,
                    If97FluidRegion::Region2,
                )
            }
        };
        let lower = If97FluidState::at_pressure_temperature_in_region(
            pressure,
            Temperature::new(temperatures[0]),
            lower_region,
            controls.pressure,
            &mut cancelled,
        )?;
        let upper = If97FluidState::at_pressure_temperature_in_region(
            pressure,
            Temperature::new(temperatures[1]),
            upper_region,
            controls.pressure,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(Self {
            kind,
            lower,
            upper,
            boundary_iterations,
        })
    }
    pub const fn kind(self) -> If97JunctionKind {
        self.kind
    }
    pub const fn lower(self) -> If97FluidState {
        self.lower
    }
    pub const fn upper(self) -> If97FluidState {
        self.upper
    }
    pub fn pressure(self) -> Pressure {
        self.lower.requested_pressure()
    }
    pub fn temperature_bounds_k(self) -> [f64; 2] {
        [
            self.lower.temperature().value(),
            self.upper.temperature().value(),
        ]
    }
    pub const fn boundary_iterations(self) -> usize {
        self.boundary_iterations
    }
    pub fn source_evaluations(self) -> usize {
        self.lower.source_evaluations() + self.upper.source_evaluations()
    }
    pub fn pressure_residuals_pa(self) -> [f64; 2] {
        [
            self.lower.pressure_residual_pa(),
            self.upper.pressure_residual_pa(),
        ]
    }
    /// Original upper-minus-lower enthalpy, never an added latent heat.
    pub fn enthalpy_change(self) -> MassSpecificEnergyQuantityV1 {
        MassSpecificEnergyQuantityV1::new(
            self.upper.specific_enthalpy().value() - self.lower.specific_enthalpy().value(),
        )
    }
    pub fn internal_energy_change(self) -> MassSpecificEnergyQuantityV1 {
        MassSpecificEnergyQuantityV1::new(
            self.upper.specific_internal_energy().value()
                - self.lower.specific_internal_energy().value(),
        )
    }
    pub fn entropy_change(self) -> MassSpecificThermalQuantityV1 {
        MassSpecificThermalQuantityV1::new(
            self.upper.specific_entropy().value() - self.lower.specific_entropy().value(),
        )
    }
    pub fn volume_change(self) -> MassSpecificVolume {
        MassSpecificVolume::new(
            self.upper.specific_volume().value() - self.lower.specific_volume().value(),
        )
    }
    pub fn gibbs_energy_change(self) -> MassSpecificEnergyQuantityV1 {
        MassSpecificEnergyQuantityV1::new(
            self.upper.specific_gibbs_energy().value() - self.lower.specific_gibbs_energy().value(),
        )
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), If97FluidError> {
    if cancelled() {
        Err(If97FluidError::Cancelled)
    } else {
        Ok(())
    }
}
fn validate_point(p: f64, t: f64) -> Result<(), If97FluidError> {
    if !p.is_finite() || p <= 0.0 || !t.is_finite() {
        return Err(If97FluidError::Invalid(
            "positive finite pressure and finite temperature required",
        ));
    }
    if p > MAXIMUM_PRESSURE_PA
        || !(MINIMUM_STABLE_TEMPERATURE_K..=REGION2_MAXIMUM_TEMPERATURE_K).contains(&t)
    {
        return Err(If97FluidError::OutsideDomain);
    }
    if p == CRITICAL_PRESSURE_PA && t == CRITICAL_TEMPERATURE_K {
        return Err(If97FluidError::CriticalSingularity);
    }
    Ok(())
}
fn boundary_temperatures(
    p: f64,
    controls: If97JunctionSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<([f64; 2], usize), If97FluidError> {
    let mut lower = REGION1_MAXIMUM_TEMPERATURE_K;
    let mut upper = 863.15_f64;
    if boundary23_pressure_raw(upper) < p {
        // The Region 3 upper-T boundary and the Region 2 >863.15 K rule supply
        // admitted one-sided endpoints even if rounded B23 is just below 100 MPa.
        let bounds = [upper, upper.next_up()];
        if bounds[1] - bounds[0] <= controls.temperature_tolerance_k {
            return Ok((bounds, 0));
        }
        return Err(If97FluidError::BoundaryNotConverged {
            iterations: 0,
            temperature_bounds_k: bounds,
        });
    }
    for iteration in 1..=controls.maximum_iterations {
        poll(cancelled)?;
        let middle = lower + 0.5 * (upper - lower);
        if middle == lower || middle == upper {
            return Err(If97FluidError::BoundaryNotConverged {
                iterations: iteration - 1,
                temperature_bounds_k: [lower, upper],
            });
        }
        if boundary23_pressure_raw(middle) <= p {
            lower = middle;
        } else {
            upper = middle;
        }
        if upper - lower <= controls.temperature_tolerance_k {
            poll(cancelled)?;
            return Ok(([lower, upper], iteration));
        }
    }
    Err(If97FluidError::BoundaryNotConverged {
        iterations: controls.maximum_iterations,
        temperature_bounds_k: [lower, upper],
    })
}

#[cfg(test)]
mod tests;
