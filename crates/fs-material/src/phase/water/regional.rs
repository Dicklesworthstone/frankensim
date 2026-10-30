//! Original-source water calorimetry across the IF97 liquid, dense and steam regions.
//!
//! A prepared isobar retains at most four sensible pieces, its physical
//! liquid/vapor plateau, and the two possible single-phase regional junctions.
//! The original equations do not match exactly at those junctions. Their finite
//! h/s/v changes are admitted by explicit budgets and reported separately from
//! sensible and latent changes; no reference shift or artificial phase is used.
//! A caloric gap has no fabricated inverse. A globally ambiguous inverse
//! refuses; heating/cooling can select the first admitted root in its direction
//! from the retained state. Every heat step still passes its physical J budget.
//!
//! This is an equilibrium endpoint/path owner, not a global Cp certificate,
//! continuously reconciled EOS, thermal-network stability proof or flow model.

use super::{WaterCoordinate, WaterEntropySolve};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::water::if97::increments::isobaric_difference;
use fs_thermochem::water::if97::region3::difference::difference as dense_difference;
use fs_thermochem::water::if97::region3::pressure::DensePressureSolve;
use fs_thermochem::water::if97::regions::{
    If97FluidError, If97FluidPhase, If97FluidSaturationPair, If97FluidState, If97JunctionKind,
    If97JunctionSolve, If97RegionalJunction,
};
use fs_thermochem::water::if97::{
    CRITICAL_PRESSURE_PA, CRITICAL_TEMPERATURE_K, MAXIMUM_PRESSURE_PA, boundary23_pressure,
    saturation_pressure, saturation_temperature,
};

/// Original-source matching allowances, separate from caloric inverse accuracy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RegionalJunctionControls {
    pub maximum_enthalpy_jump_j_kg: f64,
    pub maximum_entropy_jump_j_kg_k: f64,
    pub maximum_volume_jump_m3_kg: f64,
    pub maximum_temperature_gap_k: f64,
    pub maximum_iterations: usize,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RegionalWaterControls {
    pub minimum_temperature_k: f64,
    pub maximum_temperature_k: f64,
    pub enthalpy_tolerance_j_kg: f64,
    /// Total sensible-state queries across all candidate pieces in one inverse.
    pub maximum_iterations: usize,
    pub maximum_saturation_entropy_defect_j_kg_k: f64,
    pub pressure: DensePressureSolve,
    pub junctions: RegionalJunctionControls,
}

#[derive(Clone, Debug, PartialEq)]
pub enum RegionalWaterError {
    Source(If97FluidError),
    Invalid(&'static str),
    OutsideDomain,
    AmbiguousSaturation,
    ChartMismatch,
    JunctionBudget {
        kind: If97JunctionKind,
    },
    SaturationEntropyDefect {
        defect_j_kg_k: f64,
        maximum_j_kg_k: f64,
    },
    TemperatureGap {
        bounds_k: [f64; 2],
    },
    CaloricGap {
        coordinate: WaterCoordinate,
        target: f64,
        enclosing_values: [f64; 2],
    },
    AmbiguousCaloricState {
        coordinate: WaterCoordinate,
        target: f64,
        temperature_bounds_k: [f64; 2],
        candidates: usize,
    },
    NoDirectionalRoot,
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
impl From<If97FluidError> for RegionalWaterError {
    fn from(error: If97FluidError) -> Self {
        match error {
            If97FluidError::Cancelled => Self::Cancelled,
            If97FluidError::AmbiguousSaturation => Self::AmbiguousSaturation,
            other => Self::Source(other),
        }
    }
}
impl core::fmt::Display for RegionalWaterError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "original IF97 regional isobar: {self:?}")
    }
}
impl std::error::Error for RegionalWaterError {}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RegionalWaterPhase {
    Liquid,
    SaturatedLiquid,
    LiquidVapor,
    SaturatedVapor,
    Vapor,
    Supercritical,
    CriticalIsotherm,
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct Piece {
    lower: If97FluidState,
    upper: If97FluidState,
}
#[derive(Clone, Copy, Debug, PartialEq)]
enum Event {
    Saturation(If97FluidSaturationPair),
    Junction(If97RegionalJunction),
}
impl Event {
    fn lower(self) -> If97FluidState {
        match self {
            Self::Saturation(s) => s.liquid(),
            Self::Junction(j) => j.lower(),
        }
    }
    fn upper(self) -> If97FluidState {
        match self {
            Self::Saturation(s) => s.vapor(),
            Self::Junction(j) => j.upper(),
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq)]
enum Point {
    Single {
        source: If97FluidState,
        piece: usize,
    },
    Saturated {
        quality: f64,
    },
}

/// Immutable regional partition; successful queries allocate no heap memory.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RegionalWater {
    pressure: Pressure,
    controls: RegionalWaterControls,
    pieces: [Option<Piece>; 4],
    events: [Option<Event>; 3],
    event_count: usize,
    saturation_entropy_defect_j_kg_k: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RegionalWaterState {
    chart: RegionalWater,
    point: Point,
}
impl RegionalWaterState {
    pub const fn law(&self) -> &RegionalWater {
        &self.chart
    }
    pub fn pressure_pa(self) -> f64 {
        self.chart.pressure.value()
    }
    pub fn temperature_k(self) -> f64 {
        match self.point {
            Point::Single { source, .. } => source.temperature().value(),
            Point::Saturated { .. } => self.chart.saturation().unwrap().temperature().value(),
        }
    }
    pub const fn single_phase_source(self) -> Option<If97FluidState> {
        match self.point {
            Point::Single { source, .. } => Some(source),
            _ => None,
        }
    }
    fn property(self, read: impl Fn(If97FluidState) -> f64) -> f64 {
        match self.point {
            Point::Single { source, .. } => read(source),
            Point::Saturated { quality } => {
                let pair = self.chart.saturation().unwrap();
                let (a, b) = (read(pair.liquid()), read(pair.vapor()));
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
    pub fn specific_entropy_j_kg_k(self) -> f64 {
        self.property(|s| s.specific_entropy().value())
    }
    pub fn specific_volume_m3_kg(self) -> f64 {
        self.property(|s| s.specific_volume().value())
    }
    pub fn bulk_density_kg_m3(self) -> f64 {
        1.0 / self.specific_volume_m3_kg()
    }
    pub fn phase(self) -> RegionalWaterPhase {
        match self.point {
            Point::Saturated { quality: 0.0 } => RegionalWaterPhase::SaturatedLiquid,
            Point::Saturated { quality: 1.0 } => RegionalWaterPhase::SaturatedVapor,
            Point::Saturated { .. } => RegionalWaterPhase::LiquidVapor,
            Point::Single { source, .. } => match source.phase() {
                If97FluidPhase::Liquid => RegionalWaterPhase::Liquid,
                If97FluidPhase::SaturatedLiquid => RegionalWaterPhase::SaturatedLiquid,
                If97FluidPhase::Vapor => RegionalWaterPhase::Vapor,
                If97FluidPhase::SaturatedVapor => RegionalWaterPhase::SaturatedVapor,
                If97FluidPhase::Supercritical => RegionalWaterPhase::Supercritical,
                If97FluidPhase::CriticalIsotherm => RegionalWaterPhase::CriticalIsotherm,
            },
        }
    }
    pub const fn vapor_quality(self) -> Option<f64> {
        match self.point {
            Point::Saturated { quality } => Some(quality),
            _ => None,
        }
    }
    pub fn vapor_mass_fraction(self) -> Option<f64> {
        match self.phase() {
            RegionalWaterPhase::Liquid | RegionalWaterPhase::SaturatedLiquid => Some(0.0),
            RegionalWaterPhase::Vapor | RegionalWaterPhase::SaturatedVapor => Some(1.0),
            RegionalWaterPhase::LiquidVapor => self.vapor_quality(),
            RegionalWaterPhase::Supercritical | RegionalWaterPhase::CriticalIsotherm => None,
        }
    }
    pub fn liquid_mass_fraction(self) -> Option<f64> {
        self.vapor_mass_fraction().map(|x| 1.0 - x)
    }
    pub fn vapor_volume_fraction(self) -> Option<f64> {
        if let Some(x) = self.vapor_quality() {
            Some(
                x * self
                    .chart
                    .saturation()
                    .unwrap()
                    .vapor()
                    .specific_volume()
                    .value()
                    / self.specific_volume_m3_kg(),
            )
        } else {
            self.vapor_mass_fraction()
        }
    }
    pub fn source_pressure_residuals_pa(self) -> [f64; 2] {
        match self.point {
            Point::Single { source, .. } => [source.pressure().value() - self.pressure_pa(); 2],
            Point::Saturated { .. } => {
                let s = self.chart.saturation().unwrap();
                [
                    s.liquid().pressure().value() - self.pressure_pa(),
                    s.vapor().pressure().value() - self.pressure_pa(),
                ]
            }
        }
    }
    pub fn pressure_closure_defect_j_kg(self) -> f64 {
        self.specific_enthalpy_j_kg()
            - self.specific_internal_energy_j_kg()
            - self.pressure_pa() * self.specific_volume_m3_kg()
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RegionalWaterSolution {
    state: RegionalWaterState,
    iterations: usize,
    dense_evaluations: usize,
    residual: f64,
    candidate_count: usize,
}
impl RegionalWaterSolution {
    pub const fn state(&self) -> &RegionalWaterState {
        &self.state
    }
    pub const fn into_state(self) -> RegionalWaterState {
        self.state
    }
    pub const fn iterations(self) -> usize {
        self.iterations
    }
    /// Original dense pressure/stiffness queries; Gibbs evaluation adds zero.
    pub const fn source_evaluations(self) -> usize {
        self.dense_evaluations
    }
    pub const fn dense_source_evaluations(self) -> usize {
        self.dense_evaluations
    }
    pub const fn enthalpy_residual_j_kg(self) -> f64 {
        self.residual
    }
    pub const fn candidate_count(self) -> usize {
        self.candidate_count
    }
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RegionalWaterEntropySolution {
    solution: RegionalWaterSolution,
}
impl RegionalWaterEntropySolution {
    pub const fn state(&self) -> &RegionalWaterState {
        &self.solution.state
    }
    pub const fn into_state(self) -> RegionalWaterState {
        self.solution.state
    }
    pub const fn iterations(self) -> usize {
        self.solution.iterations
    }
    pub const fn source_evaluations(self) -> usize {
        self.solution.dense_evaluations
    }
    pub const fn dense_source_evaluations(self) -> usize {
        self.solution.dense_evaluations
    }
    pub const fn entropy_residual_j_kg_k(self) -> f64 {
        self.solution.residual
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), RegionalWaterError> {
    if cancelled() {
        Err(RegionalWaterError::Cancelled)
    } else {
        Ok(())
    }
}
fn inverse_controls(tolerance: f64, iterations: usize) -> Result<(), RegionalWaterError> {
    if !tolerance.is_finite() || tolerance < 0.0 || !(1..=128).contains(&iterations) {
        Err(RegionalWaterError::Invalid(
            "finite nonnegative tolerance and 1..=128 caloric iterations required",
        ))
    } else {
        Ok(())
    }
}
fn value(coordinate: WaterCoordinate, state: RegionalWaterState) -> f64 {
    match coordinate {
        WaterCoordinate::Enthalpy => state.specific_enthalpy_j_kg(),
        WaterCoordinate::Entropy => state.specific_entropy_j_kg_k(),
    }
}

impl RegionalWater {
    pub fn new(
        pressure: Pressure,
        controls: RegionalWaterControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, RegionalWaterError> {
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
        let j = controls.junctions;
        if !p.is_finite()
            || p <= 0.0
            || p > MAXIMUM_PRESSURE_PA
            || !lo.is_finite()
            || !hi.is_finite()
            || lo < 273.16
            || hi > 1073.15
            || lo > hi
            || ![
                controls.maximum_saturation_entropy_defect_j_kg_k,
                j.maximum_enthalpy_jump_j_kg,
                j.maximum_entropy_jump_j_kg_k,
                j.maximum_volume_jump_m3_kg,
                j.maximum_temperature_gap_k,
            ]
            .into_iter()
            .all(|x| x.is_finite() && x >= 0.0)
            || !(1..=256).contains(&j.maximum_iterations)
            || !controls.pressure.pressure_tolerance_pa.is_finite()
            || controls.pressure.pressure_tolerance_pa < 0.0
            || !controls.pressure.density_tolerance_kg_m3.is_finite()
            || controls.pressure.density_tolerance_kg_m3 <= 0.0
            || !(1..=256).contains(&controls.pressure.maximum_iterations)
        {
            return Err(RegionalWaterError::Invalid(
                "finite original-source window and nonnegative matching budgets required",
            ));
        }
        if p == CRITICAL_PRESSURE_PA && lo <= CRITICAL_TEMPERATURE_K && hi >= CRITICAL_TEMPERATURE_K
        {
            return Err(RegionalWaterError::OutsideDomain);
        }
        let mut chart = Self {
            pressure,
            controls,
            pieces: [None; 4],
            events: [None; 3],
            event_count: 0,
            saturation_entropy_defect_j_kg_k: 0.0,
        };
        // Very low-pressure vapor can be valid although no admitted liquid/vapor
        // pair exists in this temperature window. Never extrapolate saturation.
        if p < CRITICAL_PRESSURE_PA {
            // Recognize exact p(T) endpoints before asking the independently
            // rounded T(p) equation whether the pair lies inside the window.
            let anchor = [lo, hi].into_iter().find(|t| {
                *t < CRITICAL_TEMPERATURE_K
                    && saturation_pressure(Temperature::new(*t))
                        .is_ok_and(|source| source.value() == p)
            });
            let ts = anchor.or_else(|| saturation_temperature(pressure).ok().map(|t| t.value()));
            if let Some(ts) = ts {
                if lo <= ts && ts <= hi {
                    let pair = if let Some(anchor) = anchor {
                        If97FluidSaturationPair::at_temperature(
                            Temperature::new(anchor),
                            controls.pressure,
                            &mut cancelled,
                        )?
                    } else {
                        If97FluidSaturationPair::at_pressure(
                            pressure,
                            controls.pressure,
                            &mut cancelled,
                        )?
                    };
                    let defect = pair.latent_enthalpy().value() / pair.temperature().value()
                        - pair.latent_entropy().value();
                    if !defect.is_finite()
                        || defect.abs() > controls.maximum_saturation_entropy_defect_j_kg_k
                    {
                        return Err(RegionalWaterError::SaturationEntropyDefect {
                            defect_j_kg_k: defect,
                            maximum_j_kg_k: controls.maximum_saturation_entropy_defect_j_kg_k,
                        });
                    }
                    chart.saturation_entropy_defect_j_kg_k = defect;
                    chart.push_event(Event::Saturation(pair))?;
                }
            }
        }
        let join_pressure = saturation_pressure(Temperature::new(623.15))
            .map_err(If97FluidError::from)?
            .value();
        let source_junction = If97JunctionSolve {
            pressure: controls.pressure,
            temperature_tolerance_k: j.maximum_temperature_gap_k,
            maximum_iterations: j.maximum_iterations,
        };
        if p > join_pressure && lo <= 623.15 && hi >= 623.15 {
            chart.admit_junction(If97RegionalJunction::new(
                pressure,
                If97JunctionKind::LiquidDense,
                source_junction,
                &mut cancelled,
            )?)?;
        }
        // Prepare B23 only if the requested window can actually cross it.
        let reaches_b23 = hi >= 863.15
            || (hi > 623.15
                && p <= boundary23_pressure(Temperature::new(hi))
                    .map_err(If97FluidError::from)?
                    .value());
        let starts_before_b23 = lo <= 623.15
            || (lo <= 863.15
                && p >= boundary23_pressure(Temperature::new(lo))
                    .map_err(If97FluidError::from)?
                    .value());
        if p > join_pressure && reaches_b23 && starts_before_b23 {
            let junction = If97RegionalJunction::new(
                pressure,
                If97JunctionKind::DenseVapor,
                source_junction,
                &mut cancelled,
            )?;
            let [a, b] = junction.temperature_bounds_k();
            if a < lo || b > hi {
                return Err(RegionalWaterError::TemperatureGap { bounds_k: [a, b] });
            }
            chart.admit_junction(junction)?;
        }
        if lo == hi && chart.saturation().is_none() {
            return Err(RegionalWaterError::Invalid(
                "a single temperature requires a physical latent plateau",
            ));
        }
        let first = if chart.event_count > 0
            && chart.events[0].unwrap().lower().temperature().value() == lo
        {
            chart.events[0].unwrap().lower()
        } else {
            If97FluidState::at_pressure_temperature(
                pressure,
                Temperature::new(lo),
                controls.pressure,
                &mut cancelled,
            )?
        };
        let last = if chart.event_count > 0
            && chart.events[chart.event_count - 1]
                .unwrap()
                .upper()
                .temperature()
                .value()
                == hi
        {
            chart.events[chart.event_count - 1].unwrap().upper()
        } else {
            If97FluidState::at_pressure_temperature(
                pressure,
                Temperature::new(hi),
                controls.pressure,
                &mut cancelled,
            )?
        };
        let mut left = first;
        for i in 0..=chart.event_count {
            let right = if i < chart.event_count {
                chart.events[i].unwrap().lower()
            } else {
                last
            };
            if left.region() != right.region()
                || left.temperature().value() > right.temperature().value()
                || left.specific_enthalpy().value() > right.specific_enthalpy().value()
                || left.specific_entropy().value() > right.specific_entropy().value()
                || ![left.cp().value(), right.cp().value()]
                    .into_iter()
                    .all(|x| x.is_finite() && x > 0.0)
            {
                return Err(RegionalWaterError::Unrepresentable(
                    "ordered same-region sensible piece",
                ));
            }
            chart.pieces[i] = Some(Piece {
                lower: left,
                upper: right,
            });
            if i < chart.event_count {
                left = chart.events[i].unwrap().upper();
            }
        }
        poll(&mut cancelled)?;
        Ok(chart)
    }
    fn push_event(&mut self, event: Event) -> Result<(), RegionalWaterError> {
        if self.event_count == self.events.len() {
            return Err(RegionalWaterError::Unrepresentable(
                "regional partition capacity",
            ));
        }
        let mut i = self.event_count;
        while i > 0
            && self.events[i - 1].unwrap().lower().temperature().value()
                > event.lower().temperature().value()
        {
            self.events[i] = self.events[i - 1];
            i -= 1;
        }
        self.events[i] = Some(event);
        self.event_count += 1;
        Ok(())
    }
    fn admit_junction(&mut self, junction: If97RegionalJunction) -> Result<(), RegionalWaterError> {
        let c = self.controls.junctions;
        let [a, b] = junction.temperature_bounds_k();
        if junction.enthalpy_change().value().abs() > c.maximum_enthalpy_jump_j_kg
            || junction.entropy_change().value().abs() > c.maximum_entropy_jump_j_kg_k
            || junction.volume_change().value().abs() > c.maximum_volume_jump_m3_kg
            || b - a > c.maximum_temperature_gap_k
        {
            return Err(RegionalWaterError::JunctionBudget {
                kind: junction.kind(),
            });
        }
        self.push_event(Event::Junction(junction))
    }
    pub const fn controls(&self) -> &RegionalWaterControls {
        &self.controls
    }
    pub fn pressure_pa(&self) -> f64 {
        self.pressure.value()
    }
    pub fn saturation(&self) -> Option<&If97FluidSaturationPair> {
        self.events.iter().flatten().find_map(|e| match e {
            Event::Saturation(s) => Some(s),
            _ => None,
        })
    }
    pub fn junctions(&self) -> impl Iterator<Item = &If97RegionalJunction> {
        self.events.iter().flatten().filter_map(|e| match e {
            Event::Junction(j) => Some(j),
            _ => None,
        })
    }
    pub const fn saturation_entropy_defect_j_kg_k(&self) -> f64 {
        self.saturation_entropy_defect_j_kg_k
    }
    fn wrap(&self, point: Point) -> RegionalWaterState {
        RegionalWaterState {
            chart: *self,
            point,
        }
    }
    fn wrap_source(&self, source: If97FluidState, piece: usize) -> RegionalWaterState {
        if let Some(s) = self.saturation() {
            if source == s.liquid() {
                return self.wrap(Point::Saturated { quality: 0.0 });
            }
            if source == s.vapor() {
                return self.wrap(Point::Saturated { quality: 1.0 });
            }
        }
        self.wrap(Point::Single { source, piece })
    }
    fn bounds(&self, coordinate: WaterCoordinate) -> [f64; 2] {
        let mut bounds = [f64::INFINITY, f64::NEG_INFINITY];
        for (i, piece) in self
            .pieces
            .iter()
            .enumerate()
            .filter_map(|(i, p)| p.map(|p| (i, p)))
        {
            for source in [piece.lower, piece.upper] {
                let x = value(coordinate, self.wrap_source(source, i));
                bounds[0] = bounds[0].min(x);
                bounds[1] = bounds[1].max(x);
            }
        }
        bounds
    }
    pub fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        self.bounds(WaterCoordinate::Enthalpy)
    }
    pub fn entropy_bounds_j_kg_k(&self) -> [f64; 2] {
        self.bounds(WaterCoordinate::Entropy)
    }
    pub fn state_at_temperature(
        &self,
        temperature: Temperature,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RegionalWaterState, RegionalWaterError> {
        poll(&mut cancelled)?;
        let t = temperature.value();
        if !t.is_finite()
            || t < self.controls.minimum_temperature_k
            || t > self.controls.maximum_temperature_k
        {
            return Err(RegionalWaterError::OutsideDomain);
        }
        if self
            .saturation()
            .is_some_and(|s| s.temperature().value() == t)
        {
            return Err(RegionalWaterError::AmbiguousSaturation);
        }
        for (i, piece) in self
            .pieces
            .iter()
            .enumerate()
            .filter_map(|(i, p)| p.map(|p| (i, p)))
        {
            let (lo, hi) = (
                piece.lower.temperature().value(),
                piece.upper.temperature().value(),
            );
            if lo <= t && t <= hi {
                let source = if t == lo {
                    piece.lower
                } else if t == hi {
                    piece.upper
                } else {
                    If97FluidState::at_pressure_temperature_in_region(
                        self.pressure,
                        temperature,
                        piece.lower.region(),
                        self.controls.pressure,
                        &mut cancelled,
                    )?
                };
                poll(&mut cancelled)?;
                return Ok(self.wrap_source(source, i));
            }
        }
        let bounds = self
            .junctions()
            .map(|j| j.temperature_bounds_k())
            .find(|b| b[0] < t && t < b[1])
            .unwrap_or([t, t]);
        Err(RegionalWaterError::TemperatureGap { bounds_k: bounds })
    }
    pub fn state_at_quality(
        &self,
        quality: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RegionalWaterState, RegionalWaterError> {
        poll(&mut cancelled)?;
        if self.saturation().is_none() {
            return Err(RegionalWaterError::OutsideDomain);
        }
        if !quality.is_finite() || !(0.0..=1.0).contains(&quality) {
            return Err(RegionalWaterError::Invalid(
                "quality must be finite in [0,1]",
            ));
        }
        poll(&mut cancelled)?;
        Ok(self.wrap(Point::Saturated { quality }))
    }

    pub fn state_at_specific_enthalpy(
        &self,
        target: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RegionalWaterSolution, RegionalWaterError> {
        self.invert(
            target,
            WaterCoordinate::Enthalpy,
            self.controls.enthalpy_tolerance_j_kg,
            self.controls.maximum_iterations,
            None,
            &mut cancelled,
        )
    }
    pub fn state_at_specific_entropy(
        &self,
        target: f64,
        controls: WaterEntropySolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RegionalWaterEntropySolution, RegionalWaterError> {
        Ok(RegionalWaterEntropySolution {
            solution: self.invert(
                target,
                WaterCoordinate::Entropy,
                controls.entropy_tolerance_j_kg_k,
                controls.maximum_iterations,
                None,
                &mut cancelled,
            )?,
        })
    }

    // An integer path segment and its unscaled local coordinate avoid losing
    // tiny representable vapor fractions by adding them to an O(1) path index.
    fn key(&self, point: Point) -> (usize, f64) {
        match point {
            Point::Single { source, piece } => (2 * piece, source.temperature().value()),
            Point::Saturated { quality } => {
                let i = self
                    .events
                    .iter()
                    .position(|e| matches!(e, Some(Event::Saturation(_))))
                    .unwrap();
                (2 * i + 1, quality)
            }
        }
    }
    fn invert(
        &self,
        target: f64,
        coordinate: WaterCoordinate,
        tolerance: f64,
        maximum_iterations: usize,
        continuation: Option<(Point, bool)>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<RegionalWaterSolution, RegionalWaterError> {
        poll(cancelled)?;
        inverse_controls(tolerance, maximum_iterations)?;
        let bounds = self.bounds(coordinate);
        if !target.is_finite() || target < bounds[0] || target > bounds[1] {
            return Err(RegionalWaterError::OutsideDomain);
        }
        let mut candidates: [Option<(RegionalWaterState, f64)>; 5] = [None; 5];
        let mut count = 0;
        let mut iterations = 0;
        let mut density_queries = 0;
        let mut enclosing = [f64::NEG_INFINITY, f64::INFINITY];
        let mut record_endpoint = |state: RegionalWaterState| {
            let x = value(coordinate, state);
            if x <= target {
                enclosing[0] = enclosing[0].max(x);
            }
            if x >= target {
                enclosing[1] = enclosing[1].min(x);
            }
        };
        if let Some(pair) = self.saturation() {
            let lower = self.wrap(Point::Saturated { quality: 0.0 });
            let upper = self.wrap(Point::Saturated { quality: 1.0 });
            record_endpoint(lower);
            record_endpoint(upper);
            let (a, b) = (value(coordinate, lower), value(coordinate, upper));
            if a <= target && target <= b {
                let q = if target == a {
                    0.0
                } else if target == b {
                    1.0
                } else {
                    (target - a) / (b - a)
                };
                if !q.is_finite() || (target > a && q == 0.0) || (target < b && q == 1.0) {
                    return Err(RegionalWaterError::Unrepresentable("coexistence inventory"));
                }
                let state = self.wrap(Point::Saturated { quality: q });
                let residual = value(coordinate, state) - target;
                if residual.abs() > tolerance {
                    return Err(RegionalWaterError::NotConverged {
                        coordinate,
                        iterations: 0,
                        residual,
                    });
                }
                // `pair` retention establishes that this is a physical plateau.
                debug_assert!(pair.latent_enthalpy().value() > 0.0);
                candidates[count] = Some((state, residual));
                count += 1;
            }
        }
        for (index, piece) in self
            .pieces
            .iter()
            .enumerate()
            .filter_map(|(i, p)| p.map(|p| (i, p)))
        {
            let lower = self.wrap_source(piece.lower, index);
            let upper = self.wrap_source(piece.upper, index);
            record_endpoint(lower);
            record_endpoint(upper);
            if target < value(coordinate, lower) || target > value(coordinate, upper) {
                continue;
            }
            let (state, residual, used, queries) = self.solve_piece(
                index,
                target,
                coordinate,
                tolerance,
                maximum_iterations - iterations,
                cancelled,
            )?;
            iterations += used;
            density_queries += queries;
            if candidates[..count]
                .iter()
                .flatten()
                .any(|(s, _)| s.point == state.point)
            {
                continue;
            }
            candidates[count] = Some((state, residual));
            count += 1;
        }
        if count == 0 {
            return Err(RegionalWaterError::CaloricGap {
                coordinate,
                target,
                enclosing_values: enclosing,
            });
        }
        let selected = if let Some((initial, heating)) = continuation {
            let start = self.key(initial);
            let mut selected: Option<(RegionalWaterState, f64)> = None;
            for &(state, residual) in candidates[..count].iter().flatten() {
                let key = self.key(state.point);
                if (heating && key > start) || (!heating && key < start) {
                    let replace = selected.is_none_or(|(old, _)| {
                        if heating {
                            key < self.key(old.point)
                        } else {
                            key > self.key(old.point)
                        }
                    });
                    if replace {
                        selected = Some((state, residual));
                    }
                }
            }
            selected.ok_or(RegionalWaterError::NoDirectionalRoot)?
        } else if count == 1 {
            candidates[0].unwrap()
        } else {
            let mut ts = [f64::INFINITY, f64::NEG_INFINITY];
            for (state, _) in candidates[..count].iter().flatten() {
                ts[0] = ts[0].min(state.temperature_k());
                ts[1] = ts[1].max(state.temperature_k());
            }
            return Err(RegionalWaterError::AmbiguousCaloricState {
                coordinate,
                target,
                temperature_bounds_k: ts,
                candidates: count,
            });
        };
        poll(cancelled)?;
        Ok(RegionalWaterSolution {
            state: selected.0,
            residual: selected.1,
            iterations,
            dense_evaluations: density_queries,
            candidate_count: count,
        })
    }
    fn solve_piece(
        &self,
        index: usize,
        target: f64,
        coordinate: WaterCoordinate,
        tolerance: f64,
        maximum_iterations: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(RegionalWaterState, f64, usize, usize), RegionalWaterError> {
        let piece = self.pieces[index].unwrap();
        let (mut a, mut b) = (
            self.wrap_source(piece.lower, index),
            self.wrap_source(piece.upper, index),
        );
        for state in [a, b] {
            if value(coordinate, state) == target {
                return Ok((state, 0.0, 0, 0));
            }
        }
        let mut t = a.temperature_k()
            + (b.temperature_k() - a.temperature_k())
                * ((target - value(coordinate, a)) / (value(coordinate, b) - value(coordinate, a)));
        let mut queries = 0;
        let mut used = 0;
        let mut best = f64::INFINITY;
        for iteration in 1..=maximum_iterations {
            poll(cancelled)?;
            if t <= a.temperature_k() || t >= b.temperature_k() {
                t = a.temperature_k() + 0.5 * (b.temperature_k() - a.temperature_k());
            }
            if t <= a.temperature_k() || t >= b.temperature_k() {
                break;
            }
            let source = If97FluidState::at_pressure_temperature_in_region(
                self.pressure,
                Temperature::new(t),
                piece.lower.region(),
                self.controls.pressure,
                &mut *cancelled,
            )?;
            let state = self.wrap_source(source, index);
            queries += source.source_evaluations();
            used = iteration;
            let residual = value(coordinate, state) - target;
            if residual.abs() < best.abs() {
                best = residual;
            }
            if residual.abs() <= tolerance {
                return Ok((state, residual, used, queries));
            }
            let derivative = source.cp().value()
                / if coordinate == WaterCoordinate::Entropy {
                    t
                } else {
                    1.0
                };
            if !residual.is_finite() || !derivative.is_finite() || derivative <= 0.0 {
                return Err(RegionalWaterError::Unrepresentable(
                    "regional caloric inverse",
                ));
            }
            if residual < 0.0 {
                a = state;
            } else {
                b = state;
            }
            let newton = t - residual / derivative;
            t = if newton > a.temperature_k() && newton < b.temperature_k() && newton != t {
                newton
            } else {
                a.temperature_k() + 0.5 * (b.temperature_k() - a.temperature_k())
            };
        }
        for state in [a, b] {
            let residual = value(coordinate, state) - target;
            if residual.abs() < best.abs() {
                best = residual;
            }
            if residual.abs() <= tolerance {
                return Ok((state, residual, used, queries));
            }
        }
        Err(RegionalWaterError::NotConverged {
            coordinate,
            iterations: used,
            residual: best,
        })
    }

    fn changes(&self, before: Point, after: Point) -> Result<Changes, RegionalWaterError> {
        let (a, b, sign) = if self.key(before) <= self.key(after) {
            (before, after, 1.0)
        } else {
            (after, before, -1.0)
        };
        let (ka, kb) = (self.key(a).0, self.key(b).0);
        let mut changes = Changes::default();
        for (index, piece) in self
            .pieces
            .iter()
            .enumerate()
            .filter_map(|(i, p)| p.map(|p| (i, p)))
        {
            let tag = 2 * index;
            if kb < tag || ka > tag {
                continue;
            }
            let lower = if ka == tag {
                match a {
                    Point::Single { source, .. } => source,
                    _ => piece.lower,
                }
            } else {
                piece.lower
            };
            let upper = if kb == tag {
                match b {
                    Point::Single { source, .. } => source,
                    _ => piece.upper,
                }
            } else {
                piece.upper
            };
            let delta = source_changes(self.pressure, lower, upper)?;
            add(&mut changes.sensible, delta, sign);
        }
        for (index, event) in self
            .events
            .iter()
            .enumerate()
            .filter_map(|(i, e)| e.map(|e| (i, e)))
        {
            let tag = 2 * index + 1;
            if kb < tag || ka > tag {
                continue;
            }
            match event {
                Event::Saturation(pair) => {
                    let qa = if ka == tag {
                        match a {
                            Point::Saturated { quality } => quality,
                            _ => 0.0,
                        }
                    } else {
                        0.0
                    };
                    let qb = if kb == tag {
                        match b {
                            Point::Saturated { quality } => quality,
                            _ => 1.0,
                        }
                    } else {
                        1.0
                    };
                    add(
                        &mut changes.latent,
                        endpoint_changes(pair.liquid(), pair.vapor()),
                        sign * (qb - qa),
                    );
                }
                Event::Junction(junction) => {
                    add(
                        &mut changes.matching,
                        endpoint_changes(junction.lower(), junction.upper()),
                        sign,
                    );
                }
            }
        }
        Ok(changes)
    }
}

// Components are [h,u,v,s], all specific; the original matching differences
// are a source representation discrepancy, not a new thermodynamic phase.
#[derive(Default)]
struct Changes {
    sensible: [f64; 4],
    latent: [f64; 4],
    matching: [f64; 4],
}
impl Changes {
    fn total(&self) -> [f64; 4] {
        core::array::from_fn(|i| self.sensible[i] + self.latent[i] + self.matching[i])
    }
}
fn add(total: &mut [f64; 4], change: [f64; 4], scale: f64) {
    for (a, b) in total.iter_mut().zip(change) {
        *a += scale * b;
    }
}
fn endpoint_changes(a: If97FluidState, b: If97FluidState) -> [f64; 4] {
    if let (If97FluidState::Dense(a), If97FluidState::Dense(b)) = (a, b) {
        let d = dense_difference(a.state, b.state);
        [
            d.specific_enthalpy_change().value(),
            d.specific_internal_energy_change().value(),
            d.specific_volume_change().value(),
            d.specific_entropy_change().value(),
        ]
    } else {
        [
            b.specific_enthalpy().value() - a.specific_enthalpy().value(),
            b.specific_internal_energy().value() - a.specific_internal_energy().value(),
            b.specific_volume().value() - a.specific_volume().value(),
            b.specific_entropy().value() - a.specific_entropy().value(),
        ]
    }
}
fn source_changes(
    p: Pressure,
    a: If97FluidState,
    b: If97FluidState,
) -> Result<[f64; 4], RegionalWaterError> {
    match (a, b) {
        (If97FluidState::Gibbs(a), If97FluidState::Gibbs(b)) if a.region() == b.region() => {
            let d = isobaric_difference(p, a.region(), a.temperature(), b.temperature())
                .map_err(If97FluidError::from)?;
            Ok([
                d.specific_enthalpy_change().value(),
                d.specific_internal_energy_change().value(),
                d.specific_volume_change().value(),
                d.specific_entropy_change().value(),
            ])
        }
        (If97FluidState::Dense(_), If97FluidState::Dense(_)) => Ok(endpoint_changes(a, b)),
        _ => Err(RegionalWaterError::Unrepresentable(
            "sensible difference across different equations",
        )),
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RegionalWaterHeatStep {
    pub initial: RegionalWaterState,
    pub state: RegionalWaterState,
    pub mass_kg: f64,
    pub supplied_heat_j: f64,
    pub initial_phase_masses_kg: Option<[f64; 2]>,
    pub phase_masses_kg: Option<[f64; 2]>,
    pub volumes_m3: [f64; 2],
    pub volume_change_m3: f64,
    pub enthalpy_change_j: f64,
    pub internal_energy_change_j: f64,
    pub boundary_work_j: f64,
    pub sensible_enthalpy_change_j: f64,
    pub latent_enthalpy_change_j: f64,
    pub regional_enthalpy_matching_change_j: f64,
    pub source_entropy_change_j_k: f64,
    pub sensible_entropy_change_j_k: f64,
    pub latent_entropy_change_j_k: f64,
    pub regional_entropy_matching_change_j_k: f64,
    pub enthalpy_residual_j: f64,
    pub energy_residual_j: f64,
    pub represented_enthalpy_residual_j: f64,
    pub represented_energy_residual_j: f64,
    pub pressure_closure_defects_j: [f64; 2],
    pub inversion_iterations: usize,
    pub source_evaluations: usize,
    pub inversion_residual_j_kg: f64,
    /// Number of original caloric roots before directional continuation selected one.
    pub inversion_candidates: usize,
}

impl RegionalWater {
    /// Deposit heat into original H and independently close U+p_requested*V.
    ///
    /// Caloric overlap uses the first root in the heating/cooling direction.
    /// Gaps still refuse: no artificial isothermal storage is inserted at a
    /// single-phase regional boundary. Matching changes remain explicit in the
    /// report and are not identified as physical latent heat or entropy creation.
    pub fn advance_heat(
        &self,
        current: RegionalWaterState,
        mass_kg: f64,
        heat_j: f64,
        energy_tolerance_j: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RegionalWaterHeatStep, RegionalWaterError> {
        poll(&mut cancelled)?;
        if current.chart != *self {
            return Err(RegionalWaterError::ChartMismatch);
        }
        if !mass_kg.is_finite()
            || mass_kg <= 0.0
            || !heat_j.is_finite()
            || !energy_tolerance_j.is_finite()
            || energy_tolerance_j < 0.0
        {
            return Err(RegionalWaterError::Invalid(
                "positive finite mass, finite heat and nonnegative J budget required",
            ));
        }
        let h0 = current.specific_enthalpy_j_kg();
        let dh = heat_j / mass_kg;
        if !dh.is_finite() || (heat_j != 0.0 && dh == 0.0) {
            return Err(RegionalWaterError::Unrepresentable("specific heat"));
        }
        let mut target = h0 + dh;
        // Retain exactly authored source knots without silently clipping a target.
        if heat_j != 0.0 {
            for (index, piece) in self
                .pieces
                .iter()
                .enumerate()
                .filter_map(|(i, p)| p.map(|p| (i, p)))
            {
                for source in [piece.lower, piece.upper] {
                    let h = self.wrap_source(source, index).specific_enthalpy_j_kg();
                    if heat_j == mass_kg * (h - h0) {
                        target = h;
                    }
                }
            }
            if !target.is_finite() || target == h0 {
                return Err(RegionalWaterError::Unrepresentable(
                    "heat lost in H coordinate",
                ));
            }
        }
        let solved = if heat_j == 0.0 {
            RegionalWaterSolution {
                state: current,
                iterations: 0,
                dense_evaluations: 0,
                residual: 0.0,
                candidate_count: 1,
            }
        } else {
            self.invert(
                target,
                WaterCoordinate::Enthalpy,
                self.controls
                    .enthalpy_tolerance_j_kg
                    .min(energy_tolerance_j / mass_kg),
                self.controls.maximum_iterations,
                Some((current.point, heat_j > 0.0)),
                &mut cancelled,
            )?
        };
        let state = solved.state;
        let changes = self.changes(current.point, state.point)?;
        let total = changes.total().map(|x| mass_kg * x);
        let [enthalpy, internal, volume_change, entropy] = total;
        let work = self.pressure_pa() * volume_change;
        let raw_h = mass_kg * (state.specific_enthalpy_j_kg() - h0);
        let raw_u = mass_kg
            * (state.specific_internal_energy_j_kg() - current.specific_internal_energy_j_kg());
        let raw_work = self.pressure_pa()
            * mass_kg
            * (state.specific_volume_m3_kg() - current.specific_volume_m3_kg());
        let residuals = [
            enthalpy - heat_j,
            internal + work - heat_j,
            raw_h - heat_j,
            raw_u + raw_work - heat_j,
        ];
        let volumes = [current, state].map(|s| mass_kg * s.specific_volume_m3_kg());
        let closures = [current, state].map(|s| mass_kg * s.pressure_closure_defect_j_kg());
        let inventories = [current, state].map(|s| {
            s.vapor_mass_fraction()
                .map(|x| [mass_kg * (1.0 - x), mass_kg * x])
        });
        if !total.into_iter().all(f64::is_finite)
            || !work.is_finite()
            || !residuals.into_iter().all(f64::is_finite)
            || !closures.into_iter().all(f64::is_finite)
            || !volumes.into_iter().all(|x| x.is_finite() && x > 0.0)
            || inventories
                .iter()
                .flatten()
                .any(|m| !m.into_iter().all(|x| x.is_finite() && *x >= 0.0))
            || (heat_j != 0.0
                && (enthalpy == 0.0
                    || raw_h == 0.0
                    || (enthalpy > 0.0) != (heat_j > 0.0)
                    || (raw_h > 0.0) != (heat_j > 0.0)))
        {
            return Err(RegionalWaterError::Unrepresentable(
                "physical heat, volume or retained energy",
            ));
        }
        for (state, inventory) in [current, state].into_iter().zip(inventories) {
            if let (Some(x), Some(masses)) = (state.vapor_mass_fraction(), inventory) {
                if (x > 0.0 && masses[1] == 0.0) || (x < 1.0 && masses[0] == 0.0) {
                    return Err(RegionalWaterError::Unrepresentable(
                        "phase mass lost to representation",
                    ));
                }
            }
        }
        let maximum_residual = residuals.into_iter().map(f64::abs).fold(0.0, f64::max);
        if maximum_residual > energy_tolerance_j {
            return Err(RegionalWaterError::EnergyBalance {
                maximum_residual_j: maximum_residual,
                tolerance_j: energy_tolerance_j,
            });
        }
        poll(&mut cancelled)?;
        Ok(RegionalWaterHeatStep {
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
            sensible_enthalpy_change_j: mass_kg * changes.sensible[0],
            latent_enthalpy_change_j: mass_kg * changes.latent[0],
            regional_enthalpy_matching_change_j: mass_kg * changes.matching[0],
            source_entropy_change_j_k: entropy,
            sensible_entropy_change_j_k: mass_kg * changes.sensible[3],
            latent_entropy_change_j_k: mass_kg * changes.latent[3],
            regional_entropy_matching_change_j_k: mass_kg * changes.matching[3],
            enthalpy_residual_j: residuals[0],
            energy_residual_j: residuals[1],
            represented_enthalpy_residual_j: residuals[2],
            represented_energy_residual_j: residuals[3],
            pressure_closure_defects_j: closures,
            inversion_iterations: solved.iterations,
            source_evaluations: solved.dense_evaluations,
            inversion_residual_j_kg: solved.residual,
            inversion_candidates: solved.candidate_count,
        })
    }
}

#[cfg(test)]
mod tests;
