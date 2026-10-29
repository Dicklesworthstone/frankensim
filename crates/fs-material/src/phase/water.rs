//! Original IF97 liquid heating, boiling and superheated steam on one isobar.
//!
//! Regions 1, 2 and 4 share their released h/u/s reference. Saturation is limited
//! to 623.15 K: Region 3, ice, metastability and transport kinetics are absent.
//! Raw source properties are never shifted. Because the regional Gibbs fits
//! have a finite saturation matching defect, the fixed-pressure thermal entropy
//! primitive is explicitly S*=s+chi*delta_s, where delta_s=(h_v-h_l)/Ts-(s_v-s_l).
//! Chi is zero in liquid, vapor quality in coexistence, and one in vapor.
//! S* integrates dh/T on THIS chart; it is not an absolute entropy for devices
//! changing pressure. Such devices must use the original source entropy.
//!
//! The existing thermal timeline supplies heat transfer. Its geometry remains
//! reference heat-exchange geometry; constant-pressure expansion is reported
//! separately and must not be interpreted as a rigid sealed vessel.

/// Adiabatic throttles and shaft devices using original entropy across pressures.
pub mod flow;

use fs_qty::{Pressure, Temperature};
use fs_thermochem::water::if97::bounds::{
    If97HeatCapacityBound, If97HeatCapacityBoundError, If97HeatCapacityBoundSolve,
    isobaric_heat_capacity_bound,
};
use fs_thermochem::water::if97::increments::isobaric_difference;
use fs_thermochem::water::if97::{If97Error, If97Region, If97SaturationPair, If97State};

/// Numerical work and explicit source-junction admission for one pressure chart.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterChartControls {
    pub minimum_temperature_k: f64,
    pub maximum_temperature_k: f64,
    /// Absolute original-source p,h inversion residual [J/kg].
    pub enthalpy_tolerance_j_kg: f64,
    /// Maximum original regional state evaluations in each inverse, 1..=128.
    pub maximum_iterations: usize,
    /// Admitted magnitude of the complete saturation entropy matching defect.
    /// This is source-fit inconsistency, not invented entropy production.
    pub maximum_entropy_seam_j_kg_k: f64,
    /// Independent bounded interval proof of positive Cp on each sensible branch.
    pub heat_capacity: If97HeatCapacityBoundSolve,
}

/// Independent numerical controls for inversion of RAW source entropy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterEntropySolve {
    pub entropy_tolerance_j_kg_k: f64,
    pub maximum_iterations: usize,
}

/// Liquid/vapor phase identity; vapor is never represented as a solid/liquid phase.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WaterPhase {
    CompressedLiquid,
    SaturatedLiquid,
    LiquidVapor,
    SaturatedVapor,
    SuperheatedVapor,
}

/// Coordinate of an original-source inverse; residual units follow the coordinate.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WaterCoordinate {
    Enthalpy,
    Entropy,
}

/// Explicit bounded-domain, source, representability or physical-budget refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum WaterError {
    Source(If97Error),
    HeatCapacity(If97HeatCapacityBoundError),
    Invalid(&'static str),
    OutsideDomain,
    AmbiguousSaturation,
    ChartMismatch,
    EntropySeam {
        defect_j_kg_k: f64,
        maximum_j_kg_k: f64,
    },
    NotConverged {
        coordinate: WaterCoordinate,
        iterations: usize,
        residual: f64,
    },
    EnergyBalance {
        maximum_residual_j: f64,
        tolerance_j: f64,
    },
    Unrepresentable(&'static str),
    Cancelled,
}
impl core::fmt::Display for WaterError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "IF97 isobaric water: {self:?}")
    }
}
impl std::error::Error for WaterError {}
impl From<If97Error> for WaterError {
    fn from(e: If97Error) -> Self {
        Self::Source(e)
    }
}
impl From<If97HeatCapacityBoundError> for WaterError {
    fn from(e: If97HeatCapacityBoundError) -> Self {
        if e == If97HeatCapacityBoundError::Cancelled {
            Self::Cancelled
        } else {
            Self::HeatCapacity(e)
        }
    }
}

/// Immutable chart with its original saturation pair and whole-interval Cp proof.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsobaricWater {
    saturation: If97SaturationPair,
    minimum: If97State,
    maximum: If97State,
    controls: WaterChartControls,
    bounds: [Option<If97HeatCapacityBound>; 2],
    entropy_seam_j_kg_k: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum Point {
    Single(If97State),
    Saturated { quality: f64 },
}

/// Complete validated source state and its retained original pressure chart.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterState {
    chart: IsobaricWater,
    point: Point,
}
impl WaterState {
    pub const fn law(&self) -> &IsobaricWater {
        &self.chart
    }
    pub fn pressure_pa(self) -> f64 {
        self.chart.saturation.pressure().value()
    }
    pub fn temperature_k(self) -> f64 {
        match self.point {
            Point::Single(s) => s.temperature().value(),
            Point::Saturated { .. } => self.chart.saturation.temperature().value(),
        }
    }
    fn property(self, value: impl Fn(If97State) -> f64) -> f64 {
        match self.point {
            Point::Single(s) => value(s),
            Point::Saturated { quality } => {
                let a = value(*self.chart.saturation.liquid());
                let b = value(*self.chart.saturation.vapor());
                if quality == 0.0 {
                    a
                } else if quality == 1.0 {
                    b
                } else {
                    a + quality * (b - a)
                }
            }
        }
    }
    pub fn specific_enthalpy_j_kg(self) -> f64 {
        self.property(|s| s.specific_enthalpy().value())
    }
    pub fn specific_internal_energy_j_kg(self) -> f64 {
        self.property(|s| s.specific_internal_energy().value())
    }
    /// Raw IF97 source entropy, including its finite regional matching defect.
    pub fn specific_entropy_j_kg_k(self) -> f64 {
        self.property(|s| s.specific_entropy().value())
    }
    pub fn specific_volume_m3_kg(self) -> f64 {
        self.property(|s| s.specific_volume().value())
    }
    pub fn bulk_density_kg_m3(self) -> f64 {
        1.0 / self.specific_volume_m3_kg()
    }
    /// Real vapor inventory fraction, also 0/1 in the single-phase branches.
    pub fn vapor_mass_fraction(self) -> f64 {
        match self.point {
            Point::Saturated { quality } => quality,
            Point::Single(s) => {
                if s.region() == If97Region::Region1 {
                    0.0
                } else {
                    1.0
                }
            }
        }
    }
    pub fn liquid_mass_fraction(self) -> f64 {
        1.0 - self.vapor_mass_fraction()
    }
    /// Thermodynamic quality exists only at saturation, including pure endpoints.
    pub const fn vapor_quality(self) -> Option<f64> {
        match self.point {
            Point::Saturated { quality } => Some(quality),
            Point::Single(_) => None,
        }
    }
    pub fn phase(self) -> WaterPhase {
        match self.point {
            Point::Single(s) => {
                if s.region() == If97Region::Region1 {
                    WaterPhase::CompressedLiquid
                } else {
                    WaterPhase::SuperheatedVapor
                }
            }
            Point::Saturated { quality } => {
                if quality == 0.0 {
                    WaterPhase::SaturatedLiquid
                } else if quality == 1.0 {
                    WaterPhase::SaturatedVapor
                } else {
                    WaterPhase::LiquidVapor
                }
            }
        }
    }
    /// Explicit pressure-specific source matching correction, never physical heat.
    pub fn thermal_entropy_correction_j_kg_k(self) -> f64 {
        self.vapor_mass_fraction() * self.chart.entropy_seam_j_kg_k
    }
    /// Continuous dh/T primitive on this isobar only; not cross-pressure entropy.
    pub fn thermal_specific_entropy_j_kg_k(self) -> f64 {
        self.specific_entropy_j_kg_k() + self.thermal_entropy_correction_j_kg_k()
    }
}

/// Accepted p,h inverse and actual unpatched source residual.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterEnthalpySolution {
    state: WaterState,
    iterations: usize,
    residual: f64,
}
impl WaterEnthalpySolution {
    pub const fn state(&self) -> &WaterState {
        &self.state
    }
    pub const fn into_state(self) -> WaterState {
        self.state
    }
    pub const fn iterations(&self) -> usize {
        self.iterations
    }
    pub const fn enthalpy_residual_j_kg(&self) -> f64 {
        self.residual
    }
}
/// Accepted p,s inverse of RAW source entropy, for pressure-changing devices.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterEntropySolution {
    state: WaterState,
    iterations: usize,
    residual: f64,
}
impl WaterEntropySolution {
    pub const fn state(&self) -> &WaterState {
        &self.state
    }
    pub const fn into_state(self) -> WaterState {
        self.state
    }
    pub const fn iterations(&self) -> usize {
        self.iterations
    }
    pub const fn entropy_residual_j_kg_k(&self) -> f64 {
        self.residual
    }
}

/// Proposed fixed-mass heat exchange against the chart's constant pressure load.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterHeatStep {
    pub initial: WaterState,
    pub state: WaterState,
    pub mass_kg: f64,
    pub supplied_heat_j: f64,
    /// Separately rounded [liquid,vapor] inventories; total mass remains fixed.
    pub initial_phase_masses_kg: [f64; 2],
    pub phase_masses_kg: [f64; 2],
    /// Separately rounded initial/final volumes; their subtraction can lose expansion.
    pub volumes_m3: [f64; 2],
    /// Original-source volume increment evaluated before absolute-value rounding.
    pub volume_change_m3: f64,
    pub enthalpy_change_j: f64,
    pub internal_energy_change_j: f64,
    pub boundary_work_j: f64,
    pub enthalpy_residual_j: f64,
    pub energy_residual_j: f64,
    /// Actual rounded endpoint m*(h1-h0)-Q, independently admitted.
    pub represented_enthalpy_residual_j: f64,
    /// Actual rounded endpoint m*(u1-u0)+p*m*(v1-v0)-Q, independently admitted.
    pub represented_energy_residual_j: f64,
    pub source_entropy_change_j_k: f64,
    pub thermal_entropy_change_j_k: f64,
    pub entropy_seam_correction_change_j_k: f64,
    pub inversion_iterations: usize,
    pub inversion_residual_j_kg: f64,
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), WaterError> {
    if cancelled() {
        Err(WaterError::Cancelled)
    } else {
        Ok(())
    }
}
fn inverse_controls(tolerance: f64, iterations: usize) -> Result<(), WaterError> {
    if !tolerance.is_finite() || tolerance < 0.0 || !(1..=128).contains(&iterations) {
        Err(WaterError::Invalid(
            "finite nonnegative inverse tolerance and 1..=128 iterations required",
        ))
    } else {
        Ok(())
    }
}
impl WaterCoordinate {
    fn value(self, s: WaterState) -> f64 {
        match self {
            Self::Enthalpy => s.specific_enthalpy_j_kg(),
            Self::Entropy => s.specific_entropy_j_kg_k(),
        }
    }
}

impl IsobaricWater {
    /// Admit a whole liquid/coexistence/vapor chart using source interval bounds.
    /// Tmin/Tmax must bracket saturation and lie within [273.16,1073.15] K.
    /// Equal saturation endpoints permit omission of either sensible branch.
    pub fn new(
        pressure: Pressure,
        controls: WaterChartControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, WaterError> {
        poll(&mut cancelled)?;
        inverse_controls(
            controls.enthalpy_tolerance_j_kg,
            controls.maximum_iterations,
        )?;
        if !controls.minimum_temperature_k.is_finite()
            || !controls.maximum_temperature_k.is_finite()
            || controls.minimum_temperature_k < 273.16
            || controls.maximum_temperature_k > 1073.15
            || controls.minimum_temperature_k > controls.maximum_temperature_k
            || !controls.maximum_entropy_seam_j_kg_k.is_finite()
            || controls.maximum_entropy_seam_j_kg_k < 0.0
        {
            return Err(WaterError::Invalid(
                "bounded temperatures and nonnegative entropy seam budget required",
            ));
        }
        let saturation = If97SaturationPair::at_pressure(pressure)?;
        let ts = saturation.temperature().value();
        if ts > 623.15 || controls.minimum_temperature_k > ts || controls.maximum_temperature_k < ts
        {
            return Err(WaterError::OutsideDomain);
        }
        let minimum = If97State::at_pressure_temperature(
            pressure,
            Temperature::new(controls.minimum_temperature_k),
            If97Region::Region1,
        )?;
        let maximum = If97State::at_pressure_temperature(
            pressure,
            Temperature::new(controls.maximum_temperature_k),
            If97Region::Region2,
        )?;
        let (liquid, vapor) = (saturation.liquid(), saturation.vapor());
        let latent = vapor.specific_enthalpy().value() - liquid.specific_enthalpy().value();
        let ds = vapor.specific_entropy().value() - liquid.specific_entropy().value();
        let defect = latent / ts - ds;
        if !latent.is_finite()
            || latent <= 0.0
            || !ds.is_finite()
            || ds <= 0.0
            || !defect.is_finite()
        {
            return Err(WaterError::Unrepresentable("saturation interval"));
        }
        if (controls.minimum_temperature_k < ts
            && (minimum.specific_enthalpy().value() >= liquid.specific_enthalpy().value()
                || minimum.specific_entropy().value() >= liquid.specific_entropy().value()))
            || (controls.maximum_temperature_k > ts
                && (maximum.specific_enthalpy().value() <= vapor.specific_enthalpy().value()
                    || maximum.specific_entropy().value() <= vapor.specific_entropy().value()))
        {
            return Err(WaterError::Unrepresentable(
                "sensible branch coordinate endpoints",
            ));
        }
        if defect.abs() > controls.maximum_entropy_seam_j_kg_k {
            return Err(WaterError::EntropySeam {
                defect_j_kg_k: defect,
                maximum_j_kg_k: controls.maximum_entropy_seam_j_kg_k,
            });
        }
        let mut bounds = [None, None];
        for (i, region, lo, hi) in [
            (0, If97Region::Region1, controls.minimum_temperature_k, ts),
            (1, If97Region::Region2, ts, controls.maximum_temperature_k),
        ] {
            if lo < hi {
                bounds[i] = Some(isobaric_heat_capacity_bound(
                    pressure,
                    region,
                    Temperature::new(lo),
                    Temperature::new(hi),
                    controls.heat_capacity,
                    &mut cancelled,
                )?);
            }
        }
        poll(&mut cancelled)?;
        Ok(Self {
            saturation,
            minimum,
            maximum,
            controls,
            bounds,
            entropy_seam_j_kg_k: defect,
        })
    }
    pub const fn saturation(&self) -> &If97SaturationPair {
        &self.saturation
    }
    pub fn pressure_pa(&self) -> f64 {
        self.saturation.pressure().value()
    }
    pub fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        [
            self.minimum.specific_enthalpy().value(),
            self.maximum.specific_enthalpy().value(),
        ]
    }
    pub fn entropy_bounds_j_kg_k(&self) -> [f64; 2] {
        [
            self.minimum.specific_entropy().value(),
            self.maximum.specific_entropy().value(),
        ]
    }
    pub const fn saturation_entropy_defect_j_kg_k(&self) -> f64 {
        self.entropy_seam_j_kg_k
    }
    pub const fn heat_capacity_bounds(&self) -> &[Option<If97HeatCapacityBound>; 2] {
        &self.bounds
    }
    /// Proven global inverse slope from outward-rounded whole-branch Cp bounds.
    pub fn maximum_temperature_slope(&self) -> f64 {
        self.bounds
            .iter()
            .flatten()
            .map(If97HeatCapacityBound::maximum_temperature_slope)
            .fold(0.0, f64::max)
    }
    fn saturated(&self, quality: f64) -> WaterState {
        WaterState {
            chart: *self,
            point: Point::Saturated { quality },
        }
    }
    fn sourced(&self, state: If97State) -> WaterState {
        if state.temperature().value() == self.saturation.temperature().value() {
            self.saturated(if state.region() == If97Region::Region1 {
                0.0
            } else {
                1.0
            })
        } else {
            WaterState {
                chart: *self,
                point: Point::Single(state),
            }
        }
    }
    /// Resolve an authored single-phase temperature; Ts requires a quality or h/s.
    pub fn state_at_temperature(
        &self,
        temperature: Temperature,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<WaterState, WaterError> {
        poll(&mut cancelled)?;
        let t = temperature.value();
        let ts = self.saturation.temperature().value();
        if !t.is_finite()
            || t < self.controls.minimum_temperature_k
            || t > self.controls.maximum_temperature_k
        {
            return Err(WaterError::OutsideDomain);
        }
        if t == ts {
            return Err(WaterError::AmbiguousSaturation);
        }
        let state = If97State::at_pressure_temperature(
            self.saturation.pressure(),
            temperature,
            if t < ts {
                If97Region::Region1
            } else {
                If97Region::Region2
            },
        )?;
        poll(&mut cancelled)?;
        Ok(self.sourced(state))
    }
    pub fn state_at_specific_enthalpy(
        &self,
        h: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<WaterEnthalpySolution, WaterError> {
        let (state, iterations, residual) = self.invert(
            h,
            WaterCoordinate::Enthalpy,
            self.controls.enthalpy_tolerance_j_kg,
            self.controls.maximum_iterations,
            &mut cancelled,
        )?;
        Ok(WaterEnthalpySolution {
            state,
            iterations,
            residual,
        })
    }
    /// Invert raw source entropy, including the raw saturation gap, never S*.
    pub fn state_at_specific_entropy(
        &self,
        s: f64,
        controls: WaterEntropySolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<WaterEntropySolution, WaterError> {
        let (state, iterations, residual) = self.invert(
            s,
            WaterCoordinate::Entropy,
            controls.entropy_tolerance_j_kg_k,
            controls.maximum_iterations,
            &mut cancelled,
        )?;
        Ok(WaterEntropySolution {
            state,
            iterations,
            residual,
        })
    }
    fn invert(
        &self,
        target: f64,
        coordinate: WaterCoordinate,
        tolerance: f64,
        maximum_iterations: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(WaterState, usize, f64), WaterError> {
        poll(cancelled)?;
        inverse_controls(tolerance, maximum_iterations)?;
        let (minimum, maximum) = (self.sourced(self.minimum), self.sourced(self.maximum));
        let (liquid, vapor) = (self.saturated(0.0), self.saturated(1.0));
        let value = |state| coordinate.value(state);
        if !target.is_finite() || target < value(minimum) || target > value(maximum) {
            return Err(WaterError::OutsideDomain);
        }
        for state in [minimum, liquid, vapor, maximum] {
            if target == value(state) {
                poll(cancelled)?;
                return Ok((state, 0, 0.0));
            }
        }
        if target > value(liquid) && target < value(vapor) {
            let quality = (target - value(liquid)) / (value(vapor) - value(liquid));
            if !quality.is_finite() || quality <= 0.0 || quality >= 1.0 {
                return Err(WaterError::Unrepresentable("coexistence quality"));
            }
            let state = self.saturated(quality);
            let residual = value(state) - target;
            if residual.abs() > tolerance {
                return Err(WaterError::NotConverged {
                    coordinate,
                    iterations: 0,
                    residual,
                });
            }
            poll(cancelled)?;
            return Ok((state, 0, residual));
        }
        let (mut lower, mut upper, region) = if target < value(liquid) {
            (minimum, liquid, If97Region::Region1)
        } else {
            (vapor, maximum, If97Region::Region2)
        };
        let mut t = lower.temperature_k()
            + (upper.temperature_k() - lower.temperature_k())
                * ((target - value(lower)) / (value(upper) - value(lower)));
        let mut residual = f64::INFINITY;
        let mut evaluated = 0;
        for iterations in 1..=maximum_iterations {
            poll(cancelled)?;
            if t <= lower.temperature_k() || t >= upper.temperature_k() {
                t = lower.temperature_k() + 0.5 * (upper.temperature_k() - lower.temperature_k());
            }
            if t <= lower.temperature_k() || t >= upper.temperature_k() {
                break;
            }
            let source = If97State::at_pressure_temperature(
                self.saturation.pressure(),
                Temperature::new(t),
                region,
            )?;
            evaluated = iterations;
            let state = self.sourced(source);
            residual = value(state) - target;
            if residual.abs() <= tolerance {
                poll(cancelled)?;
                return Ok((state, iterations, residual));
            }
            let derivative = source.cp().value()
                / if coordinate == WaterCoordinate::Entropy {
                    t
                } else {
                    1.0
                };
            if !residual.is_finite() || !derivative.is_finite() || derivative <= 0.0 {
                return Err(WaterError::Unrepresentable("original caloric inverse"));
            }
            if residual < 0.0 {
                lower = state;
            } else {
                upper = state;
            }
            let newton = t - residual / derivative;
            t = if newton > lower.temperature_k() && newton < upper.temperature_k() && newton != t {
                newton
            } else {
                lower.temperature_k() + 0.5 * (upper.temperature_k() - lower.temperature_k())
            };
        }
        for state in [lower, upper] {
            let end_residual = value(state) - target;
            if end_residual.abs() < residual.abs() {
                residual = end_residual;
            }
            if end_residual.abs() <= tolerance {
                poll(cancelled)?;
                return Ok((state, evaluated, end_residual));
            }
        }
        Err(WaterError::NotConverged {
            coordinate,
            iterations: evaluated,
            residual,
        })
    }

    // [delta_h, delta_u, delta_v, raw_delta_s, thermal_delta_s, seam_correction].
    // Clamp only branch-integration coordinates, never states, heat or phase mass.
    fn changes(&self, before: WaterState, after: WaterState) -> Result<[f64; 6], WaterError> {
        let ts = self.saturation.temperature().value();
        let mut result = [0.0; 6];
        for (region, a, b) in [
            (
                If97Region::Region1,
                before.temperature_k().min(ts),
                after.temperature_k().min(ts),
            ),
            (
                If97Region::Region2,
                before.temperature_k().max(ts),
                after.temperature_k().max(ts),
            ),
        ] {
            if a == b {
                continue;
            }
            let d = isobaric_difference(
                self.saturation.pressure(),
                region,
                Temperature::new(a),
                Temperature::new(b),
            )?;
            result[0] += d.specific_enthalpy_change().value();
            result[1] += d.specific_internal_energy_change().value();
            result[2] += d.specific_volume_change().value();
            result[3] += d.specific_entropy_change().value();
        }
        result[4] = result[3];
        let dq = after.vapor_mass_fraction() - before.vapor_mass_fraction();
        let (l, v) = (self.saturation.liquid(), self.saturation.vapor());
        let dh = dq * (v.specific_enthalpy().value() - l.specific_enthalpy().value());
        result[0] += dh;
        result[1] +=
            dq * (v.specific_internal_energy().value() - l.specific_internal_energy().value());
        result[2] += dq * (v.specific_volume().value() - l.specific_volume().value());
        result[3] += dq * (v.specific_entropy().value() - l.specific_entropy().value());
        result[4] += dh / ts;
        result[5] = dq * self.entropy_seam_j_kg_k;
        if !result.into_iter().all(f64::is_finite) {
            return Err(WaterError::Unrepresentable("source increments"));
        }
        Ok(result)
    }
    /// Fixed-pressure dh/T integral, with the source seam explicitly accounted.
    /// Evaluated at actual inverse endpoints; each h residual is bounded by the
    /// declared inverse tolerance. No source-fit uncertainty is claimed.
    pub fn specific_entropy_change_j_kg_k(
        &self,
        initial_h: f64,
        final_h: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64, WaterError> {
        let a = self
            .state_at_specific_enthalpy(initial_h, &mut cancelled)?
            .into_state();
        let b = self
            .state_at_specific_enthalpy(final_h, &mut cancelled)?
            .into_state();
        let change = self.changes(a, b)?[4];
        poll(&mut cancelled)?;
        Ok(change)
    }
    /// Fixed mass and fixed pressure, with independent original and represented
    /// endpoint balances. Both sensible branches and latent heat occur once.
    /// All four energy defects use the caller's unchanged absolute J budget.
    pub fn advance_heat(
        &self,
        current: WaterState,
        mass_kg: f64,
        heat_j: f64,
        energy_tolerance_j: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<WaterHeatStep, WaterError> {
        poll(&mut cancelled)?;
        if current.chart != *self {
            return Err(WaterError::ChartMismatch);
        }
        if !mass_kg.is_finite()
            || mass_kg <= 0.0
            || !heat_j.is_finite()
            || !energy_tolerance_j.is_finite()
            || energy_tolerance_j < 0.0
        {
            return Err(WaterError::Invalid(
                "positive mass, finite heat and nonnegative J budget required",
            ));
        }
        let h0 = current.specific_enthalpy_j_kg();
        let dq = heat_j / mass_kg;
        if !dq.is_finite() || (heat_j != 0.0 && dq == 0.0) {
            return Err(WaterError::Unrepresentable("specific heat"));
        }
        let [lo, hi] = self.enthalpy_bounds_j_kg();
        if heat_j < mass_kg * (lo - h0) || heat_j > mass_kg * (hi - h0) {
            return Err(WaterError::OutsideDomain);
        }
        let mut target = h0 + dq;
        // Exact requests for actual source knots retain complete phase depletion
        // and sensible-domain endpoints; excess heat is never clipped.
        if heat_j != 0.0 {
            for h in [
                lo,
                self.saturated(0.0).specific_enthalpy_j_kg(),
                self.saturated(1.0).specific_enthalpy_j_kg(),
                hi,
            ] {
                if heat_j == mass_kg * (h - h0) {
                    target = h;
                    break;
                }
            }
            if target == h0 {
                return Err(WaterError::Unrepresentable(
                    "heat lost in enthalpy coordinate",
                ));
            }
        }
        let solved = if heat_j == 0.0 {
            WaterEnthalpySolution {
                state: current,
                iterations: 0,
                residual: 0.0,
            }
        } else {
            let tolerance = self
                .controls
                .enthalpy_tolerance_j_kg
                .min(energy_tolerance_j / mass_kg);
            let (state, iterations, residual) = self.invert(
                target,
                WaterCoordinate::Enthalpy,
                tolerance,
                self.controls.maximum_iterations,
                &mut cancelled,
            )?;
            WaterEnthalpySolution {
                state,
                iterations,
                residual,
            }
        };
        let state = solved.state;
        let d = self.changes(current, state)?;
        let enthalpy = mass_kg * d[0];
        let internal = mass_kg * d[1];
        let volume_change = mass_kg * d[2];
        let work = self.pressure_pa() * volume_change;
        let represented_h = mass_kg * (state.specific_enthalpy_j_kg() - h0);
        let represented_u = mass_kg
            * (state.specific_internal_energy_j_kg() - current.specific_internal_energy_j_kg());
        let represented_work = self.pressure_pa()
            * (mass_kg * (state.specific_volume_m3_kg() - current.specific_volume_m3_kg()));
        let residuals = [
            enthalpy - heat_j,
            internal + work - heat_j,
            represented_h - heat_j,
            represented_u + represented_work - heat_j,
        ];
        let volumes = [
            mass_kg * current.specific_volume_m3_kg(),
            mass_kg * state.specific_volume_m3_kg(),
        ];
        if ![enthalpy, internal, volume_change, work]
            .into_iter()
            .all(f64::is_finite)
            || !residuals.into_iter().all(f64::is_finite)
            || !volumes.into_iter().all(|v| v.is_finite() && v > 0.0)
            || (heat_j != 0.0 && (enthalpy == 0.0 || represented_h == 0.0))
        {
            return Err(WaterError::Unrepresentable(
                "physical heat or endpoint state",
            ));
        }
        let maximum_residual = residuals.into_iter().map(f64::abs).fold(0.0, f64::max);
        if maximum_residual > energy_tolerance_j {
            return Err(WaterError::EnergyBalance {
                maximum_residual_j: maximum_residual,
                tolerance_j: energy_tolerance_j,
            });
        }
        let masses = |s: WaterState| {
            [
                mass_kg * s.liquid_mass_fraction(),
                mass_kg * s.vapor_mass_fraction(),
            ]
        };
        for s in [current, state] {
            let inventory = masses(s);
            if (s.liquid_mass_fraction() > 0.0 && inventory[0] == 0.0)
                || (s.vapor_mass_fraction() > 0.0 && inventory[1] == 0.0)
            {
                return Err(WaterError::Unrepresentable("phase inventory"));
            }
        }
        poll(&mut cancelled)?;
        Ok(WaterHeatStep {
            initial: current,
            state,
            mass_kg,
            supplied_heat_j: heat_j,
            initial_phase_masses_kg: masses(current),
            phase_masses_kg: masses(state),
            volumes_m3: volumes,
            volume_change_m3: volume_change,
            enthalpy_change_j: enthalpy,
            internal_energy_change_j: internal,
            boundary_work_j: work,
            enthalpy_residual_j: residuals[0],
            energy_residual_j: residuals[1],
            represented_enthalpy_residual_j: residuals[2],
            represented_energy_residual_j: residuals[3],
            source_entropy_change_j_k: mass_kg * d[3],
            thermal_entropy_change_j_k: mass_kg * d[4],
            entropy_seam_correction_change_j_k: mass_kg * d[5],
            inversion_iterations: solved.iterations,
            inversion_residual_j_kg: solved.residual,
        })
    }
}

#[cfg(test)]
mod tests;
