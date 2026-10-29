//! Pressure-controlled dense-water heating and high-temperature boiling.
//!
//! The original IF97 Region 3 pressure/density solve supplies every single-phase
//! state; Region 4 supplies coexistence temperature. Enthalpy, internal energy,
//! entropy and volume retain their source reference. A fixed-mass pressure load
//! absorbs p*Delta V work, so supplied heat changes H rather than U alone.
//!
//! The requested pressure and the actual Helmholtz pressure are distinct finite
//! quantities. Their closure defect is exposed and admitted in the independent
//! U+pV heat balance. No endpoint pressure, enthalpy or entropy is patched to make
//! a balance pass. Subcritical coexistence uses mass-based quality and the actual
//! liquid/vapor volumes. Supercritical states have no liquid/vapor inventory.
//!
//! A temperature window is a search domain, not a certified global Cp enclosure.
//! Individual source solves can refuse unresolved near-critical states. There is
//! no Region 1/2/3 stitching, spatial flow, transport model or thermal-network
//! stability claim. The optional entropy primitive S*=s+chi*delta_s is local to
//! this isobar; pressure-changing equipment must use the original source s.

use super::{WaterCoordinate, WaterEntropySolve};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::water::if97::region3::difference::difference;
use fs_thermochem::water::if97::region3::pressure::{
    DensePressureError, DensePressureSolve, If97DenseSaturationPair,
};
use fs_thermochem::water::if97::region3::{If97DensePhase, If97DenseState};
use fs_thermochem::water::if97::{
    CRITICAL_PRESSURE_PA, CRITICAL_TEMPERATURE_K, MAXIMUM_PRESSURE_PA, boundary23_pressure,
    saturation_temperature,
};

/// Bounded temperature/enthalpy work and explicit original-source matching limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DenseIsobaricControls {
    pub minimum_temperature_k: f64,
    pub maximum_temperature_k: f64,
    pub enthalpy_tolerance_j_kg: f64,
    pub maximum_iterations: usize,
    pub maximum_entropy_seam_j_kg_k: f64,
    pub pressure: DensePressureSolve,
}

#[derive(Clone, Debug, PartialEq)]
pub enum DenseIsobaricError {
    Source(DensePressureError),
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
impl From<DensePressureError> for DenseIsobaricError {
    fn from(error: DensePressureError) -> Self {
        if error == DensePressureError::Cancelled {
            Self::Cancelled
        } else {
            Self::Source(error)
        }
    }
}
impl core::fmt::Display for DenseIsobaricError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "IF97 dense-water isobar: {self:?}")
    }
}
impl std::error::Error for DenseIsobaricError {}

/// Phase identity never assigns a vapor quality to a supercritical fluid.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DenseIsobaricPhase {
    Liquid,
    SaturatedLiquid,
    LiquidVapor,
    SaturatedVapor,
    Vapor,
    Supercritical,
    CriticalIsotherm,
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum Point {
    Single(If97DenseState),
    Saturated { quality: f64 },
}

/// Prepared Region 3 isobar with optional original high-temperature saturation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DenseIsobaricWater {
    pressure: Pressure,
    controls: DenseIsobaricControls,
    saturation: Option<If97DenseSaturationPair>,
    endpoints: [Point; 2],
    entropy_seam_j_kg_k: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DenseIsobaricState {
    chart: DenseIsobaricWater,
    point: Point,
}
impl DenseIsobaricState {
    pub const fn law(&self) -> &DenseIsobaricWater {
        &self.chart
    }
    pub fn pressure_pa(self) -> f64 {
        self.chart.pressure.value()
    }
    pub fn temperature_k(self) -> f64 {
        match self.point {
            Point::Single(state) => state.temperature().value(),
            Point::Saturated { .. } => self.chart.saturation.unwrap().temperature().value(),
        }
    }
    pub const fn single_phase_source(self) -> Option<If97DenseState> {
        match self.point {
            Point::Single(state) => Some(state),
            Point::Saturated { .. } => None,
        }
    }
    fn property(self, value: impl Fn(If97DenseState) -> f64) -> f64 {
        match self.point {
            Point::Single(state) => value(state),
            Point::Saturated { quality } => {
                let pair = self.chart.saturation.unwrap();
                let liquid = value(*pair.liquid());
                let vapor = value(*pair.vapor());
                if quality == 0.0 {
                    liquid
                } else if quality == 1.0 {
                    vapor
                } else {
                    liquid + quality * (vapor - liquid)
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
    /// Original source entropy; no pressure-specific matching adjustment.
    pub fn specific_entropy_j_kg_k(self) -> f64 {
        self.property(|s| s.specific_entropy().value())
    }
    pub fn specific_volume_m3_kg(self) -> f64 {
        self.property(|s| s.specific_volume().value())
    }
    pub fn bulk_density_kg_m3(self) -> f64 {
        1.0 / self.specific_volume_m3_kg()
    }
    pub fn phase(self) -> DenseIsobaricPhase {
        match self.point {
            Point::Single(s) => match s.phase() {
                If97DensePhase::Liquid => DenseIsobaricPhase::Liquid,
                If97DensePhase::Vapor => DenseIsobaricPhase::Vapor,
                If97DensePhase::Supercritical => DenseIsobaricPhase::Supercritical,
                If97DensePhase::CriticalIsotherm => DenseIsobaricPhase::CriticalIsotherm,
            },
            Point::Saturated { quality: 0.0 } => DenseIsobaricPhase::SaturatedLiquid,
            Point::Saturated { quality: 1.0 } => DenseIsobaricPhase::SaturatedVapor,
            Point::Saturated { .. } => DenseIsobaricPhase::LiquidVapor,
        }
    }
    pub const fn vapor_quality(self) -> Option<f64> {
        match self.point {
            Point::Saturated { quality } => Some(quality),
            Point::Single(_) => None,
        }
    }
    /// None where the source has no distinct liquid/vapor phase inventories.
    pub fn vapor_mass_fraction(self) -> Option<f64> {
        match self.point {
            Point::Saturated { quality } => Some(quality),
            Point::Single(s) => match s.phase() {
                If97DensePhase::Liquid => Some(0.0),
                If97DensePhase::Vapor => Some(1.0),
                If97DensePhase::Supercritical | If97DensePhase::CriticalIsotherm => None,
            },
        }
    }
    pub fn liquid_mass_fraction(self) -> Option<f64> {
        self.vapor_mass_fraction().map(|x| 1.0 - x)
    }
    pub fn vapor_volume_fraction(self) -> Option<f64> {
        match self.point {
            Point::Saturated { quality } => Some(
                quality
                    * self
                        .chart
                        .saturation
                        .unwrap()
                        .vapor()
                        .specific_volume()
                        .value()
                    / self.specific_volume_m3_kg(),
            ),
            Point::Single(_) => self.vapor_mass_fraction(),
        }
    }
    /// Actual original-EOS pressure errors [liquid,vapor], or a repeated single value.
    pub fn source_pressure_residuals_pa(self) -> [f64; 2] {
        match self.point {
            Point::Single(s) => [s.pressure().value() - self.pressure_pa(); 2],
            Point::Saturated { .. } => {
                let pair = self.chart.saturation.unwrap();
                [
                    pair.liquid().pressure().value() - self.pressure_pa(),
                    pair.vapor().pressure().value() - self.pressure_pa(),
                ]
            }
        }
    }
    /// h-u-p_requested*v; includes unpatched source pressure closure and rounding.
    pub fn pressure_closure_defect_j_kg(self) -> f64 {
        self.specific_enthalpy_j_kg()
            - self.specific_internal_energy_j_kg()
            - self.pressure_pa() * self.specific_volume_m3_kg()
    }
    fn thermal_entropy_coordinate(self) -> f64 {
        let Some(pair) = self.chart.saturation else {
            return 0.0;
        };
        match self.point {
            Point::Saturated { quality } => quality,
            Point::Single(s) => {
                if s.temperature().value() > pair.temperature().value() {
                    1.0
                } else {
                    0.0
                }
            }
        }
    }
    pub fn thermal_entropy_correction_j_kg_k(self) -> f64 {
        self.thermal_entropy_coordinate() * self.chart.entropy_seam_j_kg_k
    }
    /// Isobar-local dh/T primitive; never substitute it for cross-pressure entropy.
    pub fn thermal_specific_entropy_j_kg_k(self) -> f64 {
        self.specific_entropy_j_kg_k() + self.thermal_entropy_correction_j_kg_k()
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DenseIsobaricSolution {
    state: DenseIsobaricState,
    iterations: usize,
    source_evaluations: usize,
    residual: f64,
}
impl DenseIsobaricSolution {
    pub const fn state(&self) -> &DenseIsobaricState {
        &self.state
    }
    pub const fn into_state(self) -> DenseIsobaricState {
        self.state
    }
    pub const fn iterations(&self) -> usize {
        self.iterations
    }
    pub const fn source_evaluations(&self) -> usize {
        self.source_evaluations
    }
    pub const fn enthalpy_residual_j_kg(&self) -> f64 {
        self.residual
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DenseIsobaricEntropySolution {
    state: DenseIsobaricState,
    iterations: usize,
    source_evaluations: usize,
    residual: f64,
}
impl DenseIsobaricEntropySolution {
    pub const fn state(&self) -> &DenseIsobaricState {
        &self.state
    }
    pub const fn into_state(self) -> DenseIsobaricState {
        self.state
    }
    pub const fn iterations(&self) -> usize {
        self.iterations
    }
    pub const fn source_evaluations(&self) -> usize {
        self.source_evaluations
    }
    pub const fn entropy_residual_j_kg_k(&self) -> f64 {
        self.residual
    }
}

/// Fixed-mass calorimetry against the requested constant-pressure mechanical load.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DenseIsobaricHeatStep {
    pub initial: DenseIsobaricState,
    pub state: DenseIsobaricState,
    pub mass_kg: f64,
    pub supplied_heat_j: f64,
    /// None for supercritical inventory, otherwise [liquid,vapor] masses.
    pub initial_phase_masses_kg: Option<[f64; 2]>,
    pub phase_masses_kg: Option<[f64; 2]>,
    pub volumes_m3: [f64; 2],
    pub volume_change_m3: f64,
    pub enthalpy_change_j: f64,
    pub internal_energy_change_j: f64,
    /// Positive outward mechanical work, p_requested*Delta V.
    pub boundary_work_j: f64,
    pub enthalpy_residual_j: f64,
    pub energy_residual_j: f64,
    pub represented_enthalpy_residual_j: f64,
    pub represented_energy_residual_j: f64,
    /// Endpoint h-u-p_requested*v defects, each multiplied by the retained mass.
    pub pressure_closure_defects_j: [f64; 2],
    pub source_entropy_change_j_k: f64,
    pub thermal_entropy_change_j_k: f64,
    pub entropy_seam_correction_change_j_k: f64,
    pub inversion_iterations: usize,
    pub source_evaluations: usize,
    pub inversion_residual_j_kg: f64,
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), DenseIsobaricError> {
    if cancelled() {
        Err(DenseIsobaricError::Cancelled)
    } else {
        Ok(())
    }
}
fn inverse_controls(tolerance: f64, iterations: usize) -> Result<(), DenseIsobaricError> {
    if !tolerance.is_finite() || tolerance < 0.0 || !(1..=128).contains(&iterations) {
        Err(DenseIsobaricError::Invalid(
            "finite nonnegative inverse tolerance and 1..=128 iterations required",
        ))
    } else {
        Ok(())
    }
}
fn value(coordinate: WaterCoordinate, state: DenseIsobaricState) -> f64 {
    match coordinate {
        WaterCoordinate::Enthalpy => state.specific_enthalpy_j_kg(),
        WaterCoordinate::Entropy => state.specific_entropy_j_kg_k(),
    }
}

impl DenseIsobaricWater {
    /// Prepare only original Region 3. A subcritical window may include coexistence.
    /// The exact critical point and every source-region transition remain excluded.
    pub fn new(
        pressure: Pressure,
        controls: DenseIsobaricControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DenseIsobaricError> {
        poll(&mut cancelled)?;
        inverse_controls(
            controls.enthalpy_tolerance_j_kg,
            controls.maximum_iterations,
        )?;
        let (p, lo, hi) = (
            pressure.value(),
            controls.minimum_temperature_k,
            controls.maximum_temperature_k,
        );
        if !p.is_finite()
            || p <= 0.0
            || p > MAXIMUM_PRESSURE_PA
            || !lo.is_finite()
            || !hi.is_finite()
            || lo < 623.15
            || hi > 863.15
            || lo > hi
            || !controls.maximum_entropy_seam_j_kg_k.is_finite()
            || controls.maximum_entropy_seam_j_kg_k < 0.0
        {
            return Err(DenseIsobaricError::Invalid(
                "finite Region 3 pressure/temperature window and nonnegative entropy seam budget required",
            ));
        }
        let p23 = boundary23_pressure(Temperature::new(hi))
            .map_err(|_| DenseIsobaricError::OutsideDomain)?
            .value();
        if p < p23
            || (p == CRITICAL_PRESSURE_PA
                && lo <= CRITICAL_TEMPERATURE_K
                && hi >= CRITICAL_TEMPERATURE_K)
        {
            return Err(DenseIsobaricError::OutsideDomain);
        }
        let saturation = if p < CRITICAL_PRESSURE_PA {
            let ts = saturation_temperature(pressure)
                .map_err(|_| DenseIsobaricError::OutsideDomain)?
                .value();
            if lo <= ts && ts <= hi {
                Some(If97DenseSaturationPair::at_pressure(
                    pressure,
                    controls.pressure,
                    &mut cancelled,
                )?)
            } else {
                None
            }
        } else {
            None
        };
        let mut entropy_seam = 0.0;
        if let Some(pair) = saturation {
            let dh = pair.latent_enthalpy().value();
            let ds = pair.latent_entropy().value();
            entropy_seam = dh / pair.temperature().value() - ds;
            if !dh.is_finite()
                || dh <= 0.0
                || !ds.is_finite()
                || ds <= 0.0
                || !entropy_seam.is_finite()
            {
                return Err(DenseIsobaricError::Unrepresentable(
                    "saturation coordinates",
                ));
            }
            if entropy_seam.abs() > controls.maximum_entropy_seam_j_kg_k {
                return Err(DenseIsobaricError::EntropySeam {
                    defect_j_kg_k: entropy_seam,
                    maximum_j_kg_k: controls.maximum_entropy_seam_j_kg_k,
                });
            }
        }
        let mut endpoint = |t, quality| -> Result<Point, DenseIsobaricError> {
            if saturation.is_some_and(|pair| t == pair.temperature().value()) {
                Ok(Point::Saturated { quality })
            } else {
                let solved = If97DenseState::at_pressure_temperature(
                    pressure,
                    Temperature::new(t),
                    controls.pressure,
                    &mut cancelled,
                )?;
                Ok(Point::Single(solved.state))
            }
        };
        let endpoints = [endpoint(lo, 0.0)?, endpoint(hi, 1.0)?];
        let chart = Self {
            pressure,
            controls,
            saturation,
            endpoints,
            entropy_seam_j_kg_k: entropy_seam,
        };
        let [minimum, maximum] = endpoints.map(|point| chart.wrap(point));
        if minimum.specific_enthalpy_j_kg() >= maximum.specific_enthalpy_j_kg()
            || minimum.specific_entropy_j_kg_k() >= maximum.specific_entropy_j_kg_k()
        {
            return Err(DenseIsobaricError::Unrepresentable(
                "ordered caloric endpoints",
            ));
        }
        if let Some(pair) = saturation {
            for (source, endpoint, is_lower) in [
                (*pair.liquid(), minimum, true),
                (*pair.vapor(), maximum, false),
            ] {
                if endpoint.temperature_k() != pair.temperature().value()
                    && ((is_lower
                        && (endpoint.specific_enthalpy_j_kg()
                            >= source.specific_enthalpy().value()
                            || endpoint.specific_entropy_j_kg_k()
                                >= source.specific_entropy().value()))
                        || (!is_lower
                            && (endpoint.specific_enthalpy_j_kg()
                                <= source.specific_enthalpy().value()
                                || endpoint.specific_entropy_j_kg_k()
                                    <= source.specific_entropy().value())))
                {
                    return Err(DenseIsobaricError::Unrepresentable(
                        "sensible branch endpoints",
                    ));
                }
            }
        }
        poll(&mut cancelled)?;
        Ok(chart)
    }
    pub const fn controls(&self) -> &DenseIsobaricControls {
        &self.controls
    }
    pub fn pressure_pa(&self) -> f64 {
        self.pressure.value()
    }
    pub const fn saturation(&self) -> Option<&If97DenseSaturationPair> {
        self.saturation.as_ref()
    }
    pub const fn saturation_entropy_defect_j_kg_k(&self) -> f64 {
        self.entropy_seam_j_kg_k
    }
    pub fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        self.endpoints
            .map(|point| self.wrap(point).specific_enthalpy_j_kg())
    }
    pub fn entropy_bounds_j_kg_k(&self) -> [f64; 2] {
        self.endpoints
            .map(|point| self.wrap(point).specific_entropy_j_kg_k())
    }
    fn wrap(&self, point: Point) -> DenseIsobaricState {
        DenseIsobaricState {
            chart: *self,
            point,
        }
    }
    pub fn state_at_temperature(
        &self,
        temperature: Temperature,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DenseIsobaricState, DenseIsobaricError> {
        poll(&mut cancelled)?;
        let t = temperature.value();
        if !t.is_finite()
            || t < self.controls.minimum_temperature_k
            || t > self.controls.maximum_temperature_k
        {
            return Err(DenseIsobaricError::OutsideDomain);
        }
        if self
            .saturation
            .is_some_and(|pair| t == pair.temperature().value())
        {
            return Err(DenseIsobaricError::AmbiguousSaturation);
        }
        for point in self.endpoints {
            if let Point::Single(state) = point {
                if t == state.temperature().value() {
                    poll(&mut cancelled)?;
                    return Ok(self.wrap(point));
                }
            }
        }
        let solved = If97DenseState::at_pressure_temperature(
            self.pressure,
            temperature,
            self.controls.pressure,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(self.wrap(Point::Single(solved.state)))
    }
    pub fn state_at_quality(
        &self,
        quality: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DenseIsobaricState, DenseIsobaricError> {
        poll(&mut cancelled)?;
        if self.saturation.is_none() {
            return Err(DenseIsobaricError::OutsideDomain);
        }
        if !quality.is_finite() || !(0.0..=1.0).contains(&quality) {
            return Err(DenseIsobaricError::Invalid(
                "quality must be finite and in [0,1]",
            ));
        }
        poll(&mut cancelled)?;
        Ok(self.wrap(Point::Saturated { quality }))
    }
    pub fn state_at_specific_enthalpy(
        &self,
        target: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DenseIsobaricSolution, DenseIsobaricError> {
        let (state, iterations, source_evaluations, residual) = self.invert(
            target,
            WaterCoordinate::Enthalpy,
            self.controls.enthalpy_tolerance_j_kg,
            self.controls.maximum_iterations,
            &mut cancelled,
        )?;
        Ok(DenseIsobaricSolution {
            state,
            iterations,
            source_evaluations,
            residual,
        })
    }
    pub fn state_at_specific_entropy(
        &self,
        target: f64,
        controls: WaterEntropySolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DenseIsobaricEntropySolution, DenseIsobaricError> {
        let (state, iterations, source_evaluations, residual) = self.invert(
            target,
            WaterCoordinate::Entropy,
            controls.entropy_tolerance_j_kg_k,
            controls.maximum_iterations,
            &mut cancelled,
        )?;
        Ok(DenseIsobaricEntropySolution {
            state,
            iterations,
            source_evaluations,
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
    ) -> Result<(DenseIsobaricState, usize, usize, f64), DenseIsobaricError> {
        poll(cancelled)?;
        inverse_controls(tolerance, maximum_iterations)?;
        let [mut lower, mut upper] = self.endpoints.map(|point| self.wrap(point));
        let read = |state| value(coordinate, state);
        if !target.is_finite() || target < read(lower) || target > read(upper) {
            return Err(DenseIsobaricError::OutsideDomain);
        }
        for state in [lower, upper] {
            if target == read(state) {
                poll(cancelled)?;
                return Ok((state, 0, 0, 0.0));
            }
        }
        if self.saturation.is_some() {
            let liquid = self.wrap(Point::Saturated { quality: 0.0 });
            let vapor = self.wrap(Point::Saturated { quality: 1.0 });
            for state in [liquid, vapor] {
                if target == read(state) {
                    poll(cancelled)?;
                    return Ok((state, 0, 0, 0.0));
                }
            }
            if target > read(liquid) && target < read(vapor) {
                let quality = (target - read(liquid)) / (read(vapor) - read(liquid));
                if !quality.is_finite() || quality <= 0.0 || quality >= 1.0 {
                    return Err(DenseIsobaricError::Unrepresentable("coexistence inventory"));
                }
                let state = self.wrap(Point::Saturated { quality });
                let residual = read(state) - target;
                if residual.abs() > tolerance {
                    return Err(DenseIsobaricError::NotConverged {
                        coordinate,
                        iterations: 0,
                        residual,
                    });
                }
                poll(cancelled)?;
                return Ok((state, 0, 0, residual));
            }
            if target < read(liquid) {
                upper = liquid;
            } else {
                lower = vapor;
            }
        }
        let mut t = lower.temperature_k()
            + (upper.temperature_k() - lower.temperature_k())
                * ((target - read(lower)) / (read(upper) - read(lower)));
        let mut best_residual = f64::INFINITY;
        let mut evaluated = 0;
        let mut source_evaluations = 0;
        for iterations in 1..=maximum_iterations {
            poll(cancelled)?;
            if t <= lower.temperature_k() || t >= upper.temperature_k() {
                t = lower.temperature_k() + 0.5 * (upper.temperature_k() - lower.temperature_k());
            }
            if t <= lower.temperature_k() || t >= upper.temperature_k() {
                break;
            }
            let solved = If97DenseState::at_pressure_temperature(
                self.pressure,
                Temperature::new(t),
                self.controls.pressure,
                &mut *cancelled,
            )?;
            let state = self.wrap(Point::Single(solved.state));
            evaluated = iterations;
            source_evaluations += solved.evaluations;
            let residual = read(state) - target;
            if residual.abs() < best_residual.abs() {
                best_residual = residual;
            }
            if residual.abs() <= tolerance {
                poll(cancelled)?;
                return Ok((state, iterations, source_evaluations, residual));
            }
            let derivative = solved.state.cp().value()
                / if coordinate == WaterCoordinate::Entropy {
                    t
                } else {
                    1.0
                };
            if !residual.is_finite() || !derivative.is_finite() || derivative <= 0.0 {
                return Err(DenseIsobaricError::Unrepresentable(
                    "original caloric inverse",
                ));
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
            let residual = read(state) - target;
            if residual.abs() < best_residual.abs() {
                best_residual = residual;
            }
            if residual.abs() <= tolerance {
                poll(cancelled)?;
                return Ok((state, evaluated, source_evaluations, residual));
            }
        }
        Err(DenseIsobaricError::NotConverged {
            coordinate,
            iterations: evaluated,
            residual: best_residual,
        })
    }

    // [Delta h, Delta u, Delta v, Delta s] before absolute-property subtraction.
    // The original polynomial difference applies across valid source endpoints;
    // quality contributes the lever-rule difference without an invented Cp.
    fn changes(&self, before: DenseIsobaricState, after: DenseIsobaricState) -> [f64; 4] {
        let values = |a, b| {
            let d = difference(a, b);
            [
                d.specific_enthalpy_change().value(),
                d.specific_internal_energy_change().value(),
                d.specific_volume_change().value(),
                d.specific_entropy_change().value(),
            ]
        };
        match (before.point, after.point) {
            (Point::Single(a), Point::Single(b)) => values(a, b),
            (Point::Saturated { quality: a }, Point::Saturated { quality: b }) => {
                let pair = self.saturation.unwrap();
                values(*pair.liquid(), *pair.vapor()).map(|v| (b - a) * v)
            }
            (Point::Single(a), Point::Saturated { quality }) => {
                let pair = self.saturation.unwrap();
                let latent = values(*pair.liquid(), *pair.vapor());
                if a.temperature().value() > pair.temperature().value() {
                    let base = values(a, *pair.vapor());
                    core::array::from_fn(|i| base[i] - (1.0 - quality) * latent[i])
                } else {
                    let base = values(a, *pair.liquid());
                    core::array::from_fn(|i| base[i] + quality * latent[i])
                }
            }
            (Point::Saturated { quality }, Point::Single(b)) => {
                let pair = self.saturation.unwrap();
                let latent = values(*pair.liquid(), *pair.vapor());
                if b.temperature().value() > pair.temperature().value() {
                    let base = values(*pair.vapor(), b);
                    core::array::from_fn(|i| base[i] + (1.0 - quality) * latent[i])
                } else {
                    let base = values(*pair.liquid(), b);
                    core::array::from_fn(|i| base[i] - quality * latent[i])
                }
            }
        }
    }

    /// Deposit Q into H, with independent U and constant-pressure boundary work.
    /// All source and physical residuals pass before returning the proposed state.
    pub fn advance_heat(
        &self,
        current: DenseIsobaricState,
        mass_kg: f64,
        heat_j: f64,
        energy_tolerance_j: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DenseIsobaricHeatStep, DenseIsobaricError> {
        poll(&mut cancelled)?;
        if current.chart != *self {
            return Err(DenseIsobaricError::ChartMismatch);
        }
        if !mass_kg.is_finite()
            || mass_kg <= 0.0
            || !heat_j.is_finite()
            || !energy_tolerance_j.is_finite()
            || energy_tolerance_j < 0.0
        {
            return Err(DenseIsobaricError::Invalid(
                "positive finite mass, finite heat and nonnegative J budget required",
            ));
        }
        let h0 = current.specific_enthalpy_j_kg();
        let dh = heat_j / mass_kg;
        if !dh.is_finite() || (heat_j != 0.0 && dh == 0.0) {
            return Err(DenseIsobaricError::Unrepresentable("specific heat"));
        }
        let [lo, hi] = self.enthalpy_bounds_j_kg();
        if heat_j < mass_kg * (lo - h0) || heat_j > mass_kg * (hi - h0) {
            return Err(DenseIsobaricError::OutsideDomain);
        }
        let mut target = h0 + dh;
        if heat_j != 0.0 {
            for knot in [
                Some(lo),
                Some(hi),
                self.saturation
                    .map(|s| s.liquid().specific_enthalpy().value()),
                self.saturation
                    .map(|s| s.vapor().specific_enthalpy().value()),
            ]
            .into_iter()
            .flatten()
            {
                if heat_j == mass_kg * (knot - h0) {
                    target = knot;
                }
            }
            if target == h0 {
                return Err(DenseIsobaricError::Unrepresentable(
                    "heat lost in H coordinate",
                ));
            }
        }
        let (state, iterations, source_evaluations, residual) = if heat_j == 0.0 {
            (current, 0, 0, 0.0)
        } else {
            self.invert(
                target,
                WaterCoordinate::Enthalpy,
                self.controls
                    .enthalpy_tolerance_j_kg
                    .min(energy_tolerance_j / mass_kg),
                self.controls.maximum_iterations,
                &mut cancelled,
            )?
        };
        let changes = self.changes(current, state);
        let enthalpy = mass_kg * changes[0];
        let internal = mass_kg * changes[1];
        let volume_change = mass_kg * changes[2];
        let work = self.pressure_pa() * volume_change;
        let source_entropy = mass_kg * changes[3];
        let entropy_correction = mass_kg
            * (state.thermal_entropy_coordinate() - current.thermal_entropy_coordinate())
            * self.entropy_seam_j_kg_k;
        let volumes = [current, state].map(|s| mass_kg * s.specific_volume_m3_kg());
        let closure = [current, state].map(|s| mass_kg * s.pressure_closure_defect_j_kg());
        let represented_h = mass_kg * (state.specific_enthalpy_j_kg() - h0);
        let represented_u = mass_kg
            * (state.specific_internal_energy_j_kg() - current.specific_internal_energy_j_kg());
        let represented_work = self.pressure_pa()
            * mass_kg
            * (state.specific_volume_m3_kg() - current.specific_volume_m3_kg());
        let residuals = [
            enthalpy - heat_j,
            internal + work - heat_j,
            represented_h - heat_j,
            represented_u + represented_work - heat_j,
        ];
        let phase_masses = |state: DenseIsobaricState| {
            state
                .vapor_mass_fraction()
                .map(|x| [mass_kg * (1.0 - x), mass_kg * x])
        };
        let inventories = [phase_masses(current), phase_masses(state)];
        if ![
            enthalpy,
            internal,
            volume_change,
            work,
            source_entropy,
            entropy_correction,
        ]
        .into_iter()
        .all(f64::is_finite)
            || !volumes.into_iter().all(|v| v.is_finite() && v > 0.0)
            || !closure.into_iter().all(f64::is_finite)
            || !residuals.into_iter().all(f64::is_finite)
            || inventories
                .iter()
                .flatten()
                .any(|m| !m.into_iter().all(|x| x.is_finite() && *x >= 0.0))
            || (heat_j != 0.0
                && (enthalpy == 0.0
                    || represented_h == 0.0
                    || (heat_j > 0.0) != (enthalpy > 0.0)
                    || (heat_j > 0.0) != (represented_h > 0.0)))
        {
            return Err(DenseIsobaricError::Unrepresentable(
                "physical heat, inventory or endpoint state",
            ));
        }
        for (state, inventory) in [current, state].into_iter().zip(inventories) {
            if let (Some(x), Some(masses)) = (state.vapor_mass_fraction(), inventory) {
                if (x > 0.0 && masses[1] == 0.0) || (x < 1.0 && masses[0] == 0.0) {
                    return Err(DenseIsobaricError::Unrepresentable(
                        "phase mass lost to representation",
                    ));
                }
            }
        }
        let maximum_residual = residuals.into_iter().map(f64::abs).fold(0.0, f64::max);
        if maximum_residual > energy_tolerance_j {
            return Err(DenseIsobaricError::EnergyBalance {
                maximum_residual_j: maximum_residual,
                tolerance_j: energy_tolerance_j,
            });
        }
        poll(&mut cancelled)?;
        Ok(DenseIsobaricHeatStep {
            initial: current,
            state,
            mass_kg,
            supplied_heat_j: heat_j,
            initial_phase_masses_kg: inventories[0],
            phase_masses_kg: inventories[1],
            volumes_m3: volumes,
            volume_change_m3: volume_change,
            enthalpy_change_j: enthalpy,
            internal_energy_change_j: internal,
            boundary_work_j: work,
            enthalpy_residual_j: residuals[0],
            energy_residual_j: residuals[1],
            represented_enthalpy_residual_j: residuals[2],
            represented_energy_residual_j: residuals[3],
            pressure_closure_defects_j: closure,
            source_entropy_change_j_k: source_entropy,
            thermal_entropy_change_j_k: source_entropy + entropy_correction,
            entropy_seam_correction_change_j_k: entropy_correction,
            inversion_iterations: iterations,
            source_evaluations,
            inversion_residual_j_kg: residual,
        })
    }
}

#[cfg(test)]
mod tests;
