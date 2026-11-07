//! Equilibrium water in a sealed, rigid vessel using original IF97 calorics.
//!
//! Internal energy is the storage coordinate: heat changes U, while the fixed
//! mass and volume determine pressure and phase inventory. This owner uses the
//! original Regions 1, 2 and 4 in its default mode, with 273.16 <= T <= 623.15 K.
//! It resolves p(v,T) on a stable branch or the saturation mixture selected by v.
//! The explicit regional-coexistence mode retains original Region 1/2 or Region 3
//! saturation pairs for 273.16 <= T < 647.096 K and refuses departure from their
//! liquid/vapor volume envelope. Its window cannot cross the 623.15 K regional
//! fit seam. It uses the SAME internal-energy inverse and
//! heat accounts; no pressure, latent heat, EOS or phase-transfer rate is fitted.
//! Neither mode applies a constant-pressure entropy correction. Ice, metastable
//! states, wall compliance and phase kinetics require other owners.
//!
//! The temperature bracket is a bounded search window, not a certified global
//! heat-capacity enclosure. Unsupported intermediate states refuse. Pressure
//! inversions retain the actual source volume and its residual; fixed-volume
//! heat acceptance independently checks both Delta U=Q and
//! Delta H=Delta U+v_fixed*Delta p. No endpoint property is patched to force a
//! balance. Original IF97 saturation Gibbs mismatch remains visible. Regional
//! phase-pressure residuals also leave a raw h-u-p_pair*v defect; this is exposed
//! separately and is never repaired or classified as physical heat.

use super::WaterPhase;
use fs_qty::{Pressure, Temperature};
use fs_thermochem::water::MassSpecificVolume;
use fs_thermochem::water::if97::region3::pressure::DensePressureSolve;
use fs_thermochem::water::if97::regions::{
    If97FluidError, If97FluidSaturationPair, If97FluidState,
};
use fs_thermochem::water::if97::{
    CRITICAL_TEMPERATURE_K, If97Error, If97Phase, If97Region, If97SaturationPair, If97State,
    R_SPECIFIC_J_KG_K, saturation_pressure_derivative,
};

/// Explicit source bounds, residual tolerances and bounded nonlinear work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsochoricWaterControls {
    pub minimum_temperature_k: f64,
    pub maximum_temperature_k: f64,
    pub minimum_pressure_pa: f64,
    pub maximum_pressure_pa: f64,
    /// Absolute original-source volume residual, not a relative density error.
    pub volume_tolerance_m3_kg: f64,
    pub internal_energy_tolerance_j_kg: f64,
    /// Maximum Newton/bisection iterations in each p(v,T) inverse, 1..=128.
    /// The saturation pair and at most two pressure-bracket states are additional.
    pub maximum_pressure_iterations: usize,
    /// Maximum temperature evaluations in each u(v,T) inverse, 1..=128.
    pub maximum_temperature_iterations: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum IsochoricWaterError {
    Invalid(&'static str),
    Source(If97Error),
    RegionalSource(If97FluidError),
    OutsideDomain,
    VolumeOutsideDomain {
        temperature_k: f64,
        requested_m3_kg: f64,
        admitted_m3_kg: [f64; 2],
    },
    VolumeNotConverged {
        iterations: usize,
        residual_m3_kg: f64,
    },
    EnergyNotConverged {
        iterations: usize,
        residual_j_kg: f64,
    },
    ChartMismatch,
    EnergyBalance {
        maximum_residual_j: f64,
        tolerance_j: f64,
    },
    Unrepresentable(&'static str),
    Cancelled,
}
impl From<If97Error> for IsochoricWaterError {
    fn from(error: If97Error) -> Self {
        Self::Source(error)
    }
}
impl From<If97FluidError> for IsochoricWaterError {
    fn from(error: If97FluidError) -> Self {
        if error == If97FluidError::Cancelled {
            Self::Cancelled
        } else {
            Self::RegionalSource(error)
        }
    }
}
impl core::fmt::Display for IsochoricWaterError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "IF97 rigid water vessel: {self:?}")
    }
}
impl std::error::Error for IsochoricWaterError {}

#[derive(Clone, Copy, Debug, PartialEq)]
enum Source {
    Single(If97State),
    Mixture {
        pair: If97FluidSaturationPair,
        quality: f64,
    },
}
impl Source {
    fn property(self, read: impl Fn(If97FluidState) -> f64) -> f64 {
        match self {
            Self::Single(s) => read(If97FluidState::Gibbs(s)),
            Self::Mixture { pair, quality } => {
                let l = read(pair.liquid());
                let v = read(pair.vapor());
                if quality == 0.0 {
                    l
                } else if quality == 1.0 {
                    v
                } else {
                    l + quality * (v - l)
                }
            }
        }
    }
    fn temperature(self) -> f64 {
        self.property(|s| s.temperature().value())
    }
    fn pressure(self) -> f64 {
        match self {
            Self::Single(s) => s.pressure().value(),
            // Retain the pair pressure separately from finite EOS-root residuals.
            Self::Mixture { pair, .. } => pair.pressure().value(),
        }
    }
    fn volume(self) -> f64 {
        self.property(|s| s.specific_volume().value())
    }
    fn u(self) -> f64 {
        self.property(|s| s.specific_internal_energy().value())
    }
    fn h(self) -> f64 {
        self.property(|s| s.specific_enthalpy().value())
    }
    fn s(self) -> f64 {
        self.property(|s| s.specific_entropy().value())
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct Point {
    source: Source,
    cv_j_kg_k: f64,
    pressure_iterations: usize,
}

/// Prepared fixed specific volume and bounded original-source inverse window.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsochoricWater {
    volume_m3_kg: f64,
    controls: IsochoricWaterControls,
    regional_coexistence: Option<DensePressureSolve>,
    minimum: Point,
    maximum: Point,
}

/// Original equilibrium endpoint, including the retained fixed-volume owner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsochoricWaterState {
    chart: IsochoricWater,
    point: Point,
}
impl IsochoricWaterState {
    pub const fn law(&self) -> &IsochoricWater {
        &self.chart
    }
    /// Original retained single-phase IF97 source; no re-evaluation or fit.
    pub const fn single_phase_source(&self) -> Option<If97State> {
        match self.point.source {
            Source::Single(state) => Some(state),
            Source::Mixture { .. } => None,
        }
    }
    /// Original retained coexistence endpoints, including exact saturation ends.
    pub const fn saturation_pair(&self) -> Option<If97SaturationPair> {
        match self.point.source {
            Source::Mixture {
                pair: If97FluidSaturationPair::Gibbs(pair),
                ..
            } => Some(pair),
            _ => None,
        }
    }
    /// Original coexistence pair from either mode, including dense roots.
    pub const fn regional_saturation_pair(&self) -> Option<If97FluidSaturationPair> {
        match self.point.source {
            Source::Mixture { pair, .. } => Some(pair),
            Source::Single(_) => None,
        }
    }
    /// Actual phase EOS pressures minus the declared pair pressure [Pa].
    /// Single-phase Gibbs states have no paired pressure solve and return zero.
    pub fn source_pressure_residuals_pa(self) -> [f64; 2] {
        match self.regional_saturation_pair() {
            Some(pair) => [
                pair.liquid().pressure().value() - pair.pressure().value(),
                pair.vapor().pressure().value() - pair.pressure().value(),
            ],
            None => [0.0; 2],
        }
    }
    /// RAW h-u-p_pair*v [J/kg], including phase pressure-root and rounding defects.
    /// This remains an inspectable source defect, not heat or an entropy source.
    pub fn pressure_closure_defect_j_kg(self) -> f64 {
        self.specific_enthalpy_j_kg()
            - self.specific_internal_energy_j_kg()
            - self.pressure_pa() * self.specific_volume_m3_kg()
    }
    pub fn temperature_k(self) -> f64 {
        self.point.source.temperature()
    }
    pub fn pressure_pa(self) -> f64 {
        self.point.source.pressure()
    }
    pub fn specific_internal_energy_j_kg(self) -> f64 {
        self.point.source.u()
    }
    pub fn specific_enthalpy_j_kg(self) -> f64 {
        self.point.source.h()
    }
    /// Unmodified original source entropy; no pressure-dependent datum shift.
    pub fn specific_entropy_j_kg_k(self) -> f64 {
        self.point.source.s()
    }
    /// Actual source volume, never replaced by the requested fixed volume.
    pub fn specific_volume_m3_kg(self) -> f64 {
        self.point.source.volume()
    }
    pub const fn fixed_specific_volume_m3_kg(self) -> f64 {
        self.chart.volume_m3_kg
    }
    pub fn volume_residual_m3_kg(self) -> f64 {
        self.specific_volume_m3_kg() - self.chart.volume_m3_kg
    }
    pub fn bulk_density_kg_m3(self) -> f64 {
        1.0 / self.specific_volume_m3_kg()
    }
    /// Analytic dU/dT at fixed v. On coexistence this includes phase conversion.
    pub const fn heat_capacity_at_volume_j_kg_k(self) -> f64 {
        self.point.cv_j_kg_k
    }
    pub const fn pressure_iterations(self) -> usize {
        self.point.pressure_iterations
    }
    pub fn vapor_quality(self) -> Option<f64> {
        match self.point.source {
            Source::Mixture { quality, .. } => Some(quality),
            Source::Single(s) => match s.phase() {
                If97Phase::SaturatedLiquid => Some(0.0),
                If97Phase::SaturatedVapor => Some(1.0),
                _ => None,
            },
        }
    }
    pub fn vapor_mass_fraction(self) -> f64 {
        match self.point.source {
            Source::Mixture { quality, .. } => quality,
            Source::Single(s) => {
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
    pub fn phase(self) -> WaterPhase {
        match self.vapor_quality() {
            Some(x) if x == 0.0 => WaterPhase::SaturatedLiquid,
            Some(x) if x == 1.0 => WaterPhase::SaturatedVapor,
            Some(_) => WaterPhase::LiquidVapor,
            None if self.vapor_mass_fraction() == 0.0 => WaterPhase::CompressedLiquid,
            None => WaterPhase::SuperheatedVapor,
        }
    }
    /// Original g_v-g_l at a coexistence state; not artificial entropy production.
    pub fn saturation_gibbs_gap_j_kg(self) -> Option<f64> {
        match self.point.source {
            Source::Mixture { pair, .. } => Some(pair.gibbs_energy_gap().value()),
            Source::Single(_) => None,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsochoricWaterSolution {
    state: IsochoricWaterState,
    temperature_iterations: usize,
    pressure_iterations: usize,
    residual_j_kg: f64,
}
impl IsochoricWaterSolution {
    pub const fn state(&self) -> &IsochoricWaterState {
        &self.state
    }
    pub const fn into_state(self) -> IsochoricWaterState {
        self.state
    }
    pub const fn temperature_iterations(&self) -> usize {
        self.temperature_iterations
    }
    /// Total pressure iterations of this inverse; prepared endpoints are reused.
    pub const fn pressure_iterations(&self) -> usize {
        self.pressure_iterations
    }
    pub const fn internal_energy_residual_j_kg(&self) -> f64 {
        self.residual_j_kg
    }
}

/// Accepted fixed-mass, fixed-volume heat exchange with independent balances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsochoricWaterHeatStep {
    pub initial: IsochoricWaterState,
    pub state: IsochoricWaterState,
    pub mass_kg: f64,
    pub fixed_volume_m3: f64,
    pub supplied_heat_j: f64,
    pub internal_energy_change_j: f64,
    pub enthalpy_change_j: f64,
    /// V_fixed * (p_final-p_initial), a change of H-U, not boundary work.
    pub pressure_enthalpy_change_j: f64,
    /// Exactly zero: the boundary is rigid, not a constant-pressure load.
    pub boundary_work_j: f64,
    pub internal_energy_residual_j: f64,
    /// Delta H - V_fixed Delta p - Q, independently checked in joules.
    pub enthalpy_residual_j: f64,
    /// Actual source endpoint volumes minus the fixed vessel volume.
    pub endpoint_volume_residuals_m3: [f64; 2],
    /// RAW M*(h-u-p_pair*v) at initial and final endpoints [J].
    /// Their change remains inside enthalpy_residual_j; no source value is repaired.
    pub endpoint_pressure_closure_defects_j: [f64; 2],
    pub initial_phase_masses_kg: [f64; 2],
    pub phase_masses_kg: [f64; 2],
    /// Fluid entropy change only; wall/reservoir entropy is unspecified.
    pub source_entropy_change_j_k: f64,
    pub temperature_iterations: usize,
    pub pressure_iterations: usize,
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), IsochoricWaterError> {
    if cancelled() {
        Err(IsochoricWaterError::Cancelled)
    } else {
        Ok(())
    }
}

impl IsochoricWaterControls {
    fn validate(self, regional_coexistence: bool) -> Result<(), IsochoricWaterError> {
        if !self.minimum_temperature_k.is_finite()
            || !self.maximum_temperature_k.is_finite()
            || self.minimum_temperature_k < 273.16
            || (!regional_coexistence && self.maximum_temperature_k > 623.15)
            || (regional_coexistence && self.maximum_temperature_k >= CRITICAL_TEMPERATURE_K)
            || (regional_coexistence
                && self.minimum_temperature_k <= 623.15
                && self.maximum_temperature_k > 623.15)
            || self.minimum_temperature_k >= self.maximum_temperature_k
            || !self.minimum_pressure_pa.is_finite()
            || !self.maximum_pressure_pa.is_finite()
            || self.minimum_pressure_pa <= 0.0
            || self.maximum_pressure_pa > 100e6
            || self.minimum_pressure_pa >= self.maximum_pressure_pa
            || !self.volume_tolerance_m3_kg.is_finite()
            || self.volume_tolerance_m3_kg < 0.0
            || !self.internal_energy_tolerance_j_kg.is_finite()
            || self.internal_energy_tolerance_j_kg < 0.0
            || !(1..=128).contains(&self.maximum_pressure_iterations)
            || !(1..=128).contains(&self.maximum_temperature_iterations)
        {
            Err(IsochoricWaterError::Invalid(
                "ordered source-domain bounds, nonnegative residual budgets and 1..=128 iterations required",
            ))
        } else {
            Ok(())
        }
    }
}

impl IsochoricWater {
    pub fn new(
        specific_volume: MassSpecificVolume,
        controls: IsochoricWaterControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, IsochoricWaterError> {
        Self::prepare(specific_volume, controls, None, cancelled)
    }

    /// Original regional saturation pairs on a bounded coexistence isochore.
    ///
    /// Every evaluated temperature must satisfy 273.16 <= T < 647.096 K,
    /// the declared pressure window, and v_liquid <= v_fixed <= v_vapor.
    /// Exact phase endpoints retain the pair; heating or rebinding beyond its
    /// envelope refuses instead of entering an unmodeled single-phase branch.
    /// The window must stay entirely on one side of the 623.15 K IF97 fit seam:
    /// maximum <= 623.15 K or minimum > 623.15 K. No globally monotone caloric
    /// inverse is claimed across that source discontinuity.
    /// The supplied dense pressure controls and all source residuals are retained.
    pub fn new_regional_coexistence(
        specific_volume: MassSpecificVolume,
        controls: IsochoricWaterControls,
        pressure: DensePressureSolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, IsochoricWaterError> {
        Self::prepare(specific_volume, controls, Some(pressure), cancelled)
    }

    /// Rebuild a fixed-volume chart without changing its original source mode.
    /// Explicit replacement budgets allow a finite cell to bind its new mass.
    pub fn rebind_specific_volume(
        &self,
        specific_volume: MassSpecificVolume,
        controls: IsochoricWaterControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, IsochoricWaterError> {
        Self::prepare(
            specific_volume,
            controls,
            self.regional_coexistence,
            cancelled,
        )
    }

    /// Dense-root controls only for the explicitly selected coexistence mode.
    pub const fn regional_coexistence_controls(&self) -> Option<DensePressureSolve> {
        self.regional_coexistence
    }

    fn prepare(
        specific_volume: MassSpecificVolume,
        controls: IsochoricWaterControls,
        regional_coexistence: Option<DensePressureSolve>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, IsochoricWaterError> {
        poll(&mut cancelled)?;
        controls.validate(regional_coexistence.is_some())?;
        if let Some(pressure) = regional_coexistence {
            if !pressure.pressure_tolerance_pa.is_finite()
                || pressure.pressure_tolerance_pa < 0.0
                || !pressure.density_tolerance_kg_m3.is_finite()
                || pressure.density_tolerance_kg_m3 <= 0.0
                || !(1..=256).contains(&pressure.maximum_iterations)
            {
                return Err(IsochoricWaterError::Invalid(
                    "regional coexistence requires bounded original dense pressure controls",
                ));
            }
        }
        let v = specific_volume.value();
        if !v.is_finite() || v <= 0.0 || !(1.0 / v).is_finite() {
            return Err(IsochoricWaterError::Invalid(
                "positive representable specific volume required",
            ));
        }
        let minimum = resolve(
            v,
            controls,
            regional_coexistence,
            controls.minimum_temperature_k,
            &mut cancelled,
        )?;
        let maximum = resolve(
            v,
            controls,
            regional_coexistence,
            controls.maximum_temperature_k,
            &mut cancelled,
        )?;
        if minimum.source.u() >= maximum.source.u() {
            return Err(IsochoricWaterError::Unrepresentable(
                "ordered internal-energy endpoints",
            ));
        }
        poll(&mut cancelled)?;
        Ok(Self {
            volume_m3_kg: v,
            controls,
            regional_coexistence,
            minimum,
            maximum,
        })
    }
    pub const fn fixed_specific_volume_m3_kg(&self) -> f64 {
        self.volume_m3_kg
    }
    pub const fn controls(&self) -> &IsochoricWaterControls {
        &self.controls
    }
    pub fn internal_energy_bounds_j_kg(&self) -> [f64; 2] {
        [self.minimum.source.u(), self.maximum.source.u()]
    }
    pub fn state_at_temperature(
        &self,
        temperature: Temperature,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<IsochoricWaterState, IsochoricWaterError> {
        poll(&mut cancelled)?;
        let t = temperature.value();
        if !t.is_finite()
            || t < self.controls.minimum_temperature_k
            || t > self.controls.maximum_temperature_k
        {
            return Err(IsochoricWaterError::OutsideDomain);
        }
        let point = if t == self.controls.minimum_temperature_k {
            self.minimum
        } else if t == self.controls.maximum_temperature_k {
            self.maximum
        } else {
            resolve(
                self.volume_m3_kg,
                self.controls,
                self.regional_coexistence,
                t,
                &mut cancelled,
            )?
        };
        poll(&mut cancelled)?;
        Ok(IsochoricWaterState {
            chart: *self,
            point,
        })
    }
    pub fn state_at_specific_internal_energy(
        &self,
        internal_energy_j_kg: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<IsochoricWaterSolution, IsochoricWaterError> {
        self.invert(
            internal_energy_j_kg,
            self.controls.internal_energy_tolerance_j_kg,
            &mut cancelled,
        )
    }
    fn invert(
        &self,
        target: f64,
        tolerance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<IsochoricWaterSolution, IsochoricWaterError> {
        poll(cancelled)?;
        let (mut lower, mut upper) = (self.minimum, self.maximum);
        if !target.is_finite() || target < lower.source.u() || target > upper.source.u() {
            return Err(IsochoricWaterError::OutsideDomain);
        }
        let answer = |point, temperature_iterations, pressure_iterations, residual_j_kg| {
            IsochoricWaterSolution {
                state: IsochoricWaterState {
                    chart: *self,
                    point,
                },
                temperature_iterations,
                pressure_iterations,
                residual_j_kg,
            }
        };
        for p in [lower, upper] {
            if target == p.source.u() {
                poll(cancelled)?;
                return Ok(answer(p, 0, 0, 0.0));
            }
        }
        let mut t = lower.source.temperature()
            + (upper.source.temperature() - lower.source.temperature())
                * ((target - lower.source.u()) / (upper.source.u() - lower.source.u()));
        let mut pressure_iterations = 0;
        let mut temperature_iterations = 0;
        let mut residual = f64::INFINITY;
        for i in 1..=self.controls.maximum_temperature_iterations {
            poll(cancelled)?;
            if t <= lower.source.temperature() || t >= upper.source.temperature() {
                t = lower.source.temperature()
                    + 0.5 * (upper.source.temperature() - lower.source.temperature());
            }
            if t <= lower.source.temperature() || t >= upper.source.temperature() {
                break;
            }
            let point = resolve(
                self.volume_m3_kg,
                self.controls,
                self.regional_coexistence,
                t,
                cancelled,
            )?;
            pressure_iterations += point.pressure_iterations;
            temperature_iterations = i;
            residual = point.source.u() - target;
            if residual.abs() <= tolerance {
                poll(cancelled)?;
                return Ok(answer(point, i, pressure_iterations, residual));
            }
            if residual < 0.0 {
                lower = point;
            } else {
                upper = point;
            }
            let newton = t - residual / point.cv_j_kg_k;
            t = if newton > lower.source.temperature()
                && newton < upper.source.temperature()
                && newton != t
            {
                newton
            } else {
                lower.source.temperature()
                    + 0.5 * (upper.source.temperature() - lower.source.temperature())
            };
        }
        for point in [lower, upper] {
            let r = point.source.u() - target;
            if r.abs() < residual.abs() {
                residual = r;
            }
            if r.abs() <= tolerance {
                poll(cancelled)?;
                return Ok(answer(
                    point,
                    temperature_iterations,
                    pressure_iterations,
                    r,
                ));
            }
        }
        Err(IsochoricWaterError::EnergyNotConverged {
            iterations: temperature_iterations,
            residual_j_kg: residual,
        })
    }

    /// Quasistatic heat exchange in a sealed rigid vessel. No mass/shaft flows.
    pub fn advance_heat(
        &self,
        current: IsochoricWaterState,
        mass_kg: f64,
        supplied_heat_j: f64,
        energy_tolerance_j: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<IsochoricWaterHeatStep, IsochoricWaterError> {
        poll(&mut cancelled)?;
        if current.chart != *self {
            return Err(IsochoricWaterError::ChartMismatch);
        }
        if !mass_kg.is_finite()
            || mass_kg <= 0.0
            || !supplied_heat_j.is_finite()
            || !energy_tolerance_j.is_finite()
            || energy_tolerance_j < 0.0
        {
            return Err(IsochoricWaterError::Invalid(
                "positive mass, finite heat and nonnegative joule budget required",
            ));
        }
        let u0 = current.specific_internal_energy_j_kg();
        let dq = supplied_heat_j / mass_kg;
        if !dq.is_finite() || (supplied_heat_j != 0.0 && dq == 0.0) {
            return Err(IsochoricWaterError::Unrepresentable("specific heat"));
        }
        let [lo, hi] = self.internal_energy_bounds_j_kg();
        if supplied_heat_j < mass_kg * (lo - u0) || supplied_heat_j > mass_kg * (hi - u0) {
            return Err(IsochoricWaterError::OutsideDomain);
        }
        let mut target = u0 + dq;
        if supplied_heat_j != 0.0 {
            for u in [lo, hi] {
                if supplied_heat_j == mass_kg * (u - u0) {
                    target = u;
                    break;
                }
            }
            if target == u0 {
                return Err(IsochoricWaterError::Unrepresentable(
                    "heat lost in internal-energy coordinate",
                ));
            }
        }
        let solution = if supplied_heat_j == 0.0 {
            IsochoricWaterSolution {
                state: current,
                temperature_iterations: 0,
                pressure_iterations: 0,
                residual_j_kg: 0.0,
            }
        } else {
            self.invert(
                target,
                self.controls
                    .internal_energy_tolerance_j_kg
                    .min(energy_tolerance_j / mass_kg),
                &mut cancelled,
            )?
        };
        let state = solution.state;
        let volume = mass_kg * self.volume_m3_kg;
        let internal = mass_kg * (state.specific_internal_energy_j_kg() - u0);
        let enthalpy =
            mass_kg * (state.specific_enthalpy_j_kg() - current.specific_enthalpy_j_kg());
        let pressure_enthalpy = volume * (state.pressure_pa() - current.pressure_pa());
        let residuals = [
            internal - supplied_heat_j,
            enthalpy - pressure_enthalpy - supplied_heat_j,
        ];
        let volume_residuals = [
            mass_kg * current.volume_residual_m3_kg(),
            mass_kg * state.volume_residual_m3_kg(),
        ];
        let entropy =
            mass_kg * (state.specific_entropy_j_kg_k() - current.specific_entropy_j_kg_k());
        let pressure_closure_defects = [
            mass_kg * current.pressure_closure_defect_j_kg(),
            mass_kg * state.pressure_closure_defect_j_kg(),
        ];
        if !volume.is_finite()
            || volume <= 0.0
            || ![internal, enthalpy, pressure_enthalpy, entropy]
                .into_iter()
                .all(f64::is_finite)
            || !residuals
                .into_iter()
                .chain(volume_residuals)
                .chain(pressure_closure_defects)
                .all(f64::is_finite)
            || (supplied_heat_j != 0.0 && internal == 0.0)
        {
            return Err(IsochoricWaterError::Unrepresentable(
                "rigid-vessel heat or endpoint state",
            ));
        }
        let maximum_residual = residuals.into_iter().map(f64::abs).fold(0.0, f64::max);
        if maximum_residual > energy_tolerance_j {
            return Err(IsochoricWaterError::EnergyBalance {
                maximum_residual_j: maximum_residual,
                tolerance_j: energy_tolerance_j,
            });
        }
        let masses = |s: IsochoricWaterState| {
            [
                mass_kg * s.liquid_mass_fraction(),
                mass_kg * s.vapor_mass_fraction(),
            ]
        };
        for s in [current, state] {
            let m = masses(s);
            if (s.liquid_mass_fraction() > 0.0 && m[0] == 0.0)
                || (s.vapor_mass_fraction() > 0.0 && m[1] == 0.0)
            {
                return Err(IsochoricWaterError::Unrepresentable("phase inventory"));
            }
        }
        poll(&mut cancelled)?;
        Ok(IsochoricWaterHeatStep {
            initial: current,
            state,
            mass_kg,
            fixed_volume_m3: volume,
            supplied_heat_j,
            internal_energy_change_j: internal,
            enthalpy_change_j: enthalpy,
            pressure_enthalpy_change_j: pressure_enthalpy,
            boundary_work_j: 0.0,
            internal_energy_residual_j: residuals[0],
            enthalpy_residual_j: residuals[1],
            endpoint_volume_residuals_m3: volume_residuals,
            endpoint_pressure_closure_defects_j: pressure_closure_defects,
            initial_phase_masses_kg: masses(current),
            phase_masses_kg: masses(state),
            source_entropy_change_j_k: entropy,
            temperature_iterations: solution.temperature_iterations,
            pressure_iterations: solution.pressure_iterations,
        })
    }
}

fn resolve(
    volume: f64,
    controls: IsochoricWaterControls,
    regional_coexistence: Option<DensePressureSolve>,
    temperature: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Point, IsochoricWaterError> {
    poll(cancelled)?;
    if let Some(pressure) = regional_coexistence {
        let pair = If97FluidSaturationPair::at_temperature(
            Temperature::new(temperature),
            pressure,
            &mut *cancelled,
        )?;
        return coexistence_point(volume, controls, temperature, pair, cancelled);
    }
    let pair = If97SaturationPair::at_temperature(Temperature::new(temperature))?;
    let (l, v) = (*pair.liquid(), *pair.vapor());
    let (vl, vv, ps) = (
        l.specific_volume().value(),
        v.specific_volume().value(),
        pair.pressure().value(),
    );
    if volume >= vl && volume <= vv {
        return coexistence_point(
            volume,
            controls,
            temperature,
            If97FluidSaturationPair::Gibbs(pair),
            cancelled,
        );
    }
    let region = if volume < vl {
        If97Region::Region1
    } else {
        If97Region::Region2
    };
    let (low_p, high_p) = if region == If97Region::Region1 {
        (
            ps.max(controls.minimum_pressure_pa),
            controls.maximum_pressure_pa,
        )
    } else {
        (
            controls.minimum_pressure_pa,
            ps.min(controls.maximum_pressure_pa),
        )
    };
    if low_p > high_p {
        return Err(IsochoricWaterError::OutsideDomain);
    }
    let evaluate = |p| {
        If97State::at_pressure_temperature(Pressure::new(p), Temperature::new(temperature), region)
    };
    let mut lower = if low_p == ps { l } else { evaluate(low_p)? };
    let mut upper = if high_p == ps { v } else { evaluate(high_p)? };
    // A saturated endpoint must use the selected Gibbs branch, not its other phase.
    if region == If97Region::Region1 && high_p == ps {
        upper = l;
    }
    if region == If97Region::Region2 && low_p == ps {
        lower = v;
    }
    let low_v = lower.specific_volume().value();
    let high_v = upper.specific_volume().value();
    if volume > low_v || volume < high_v {
        return Err(IsochoricWaterError::VolumeOutsideDomain {
            temperature_k: temperature,
            requested_m3_kg: volume,
            admitted_m3_kg: [high_v, low_v],
        });
    }
    let accepted = |s: If97State, n| Point {
        source: Source::Single(s),
        cv_j_kg_k: s.cv().value(),
        pressure_iterations: n,
    };
    for s in [lower, upper] {
        if (s.specific_volume().value() - volume).abs() <= controls.volume_tolerance_m3_kg {
            poll(cancelled)?;
            return Ok(accepted(s, 0));
        }
    }
    // Cheap physical seeds only start Newton; accepted properties always come
    // from the original equation with independently checked volume residual.
    let mut p = if region == If97Region::Region2 {
        R_SPECIFIC_J_KG_K * temperature / volume
    } else {
        lower.pressure().value() + (volume - low_v) / lower.dv_dp_at_temperature().value()
    };
    let mut residual = f64::INFINITY;
    let mut evaluated = 0;
    for i in 1..=controls.maximum_pressure_iterations {
        poll(cancelled)?;
        let (lo, hi) = (lower.pressure().value(), upper.pressure().value());
        if !p.is_finite() || p <= lo || p >= hi {
            p = lo + 0.5 * (hi - lo);
        }
        if p <= lo || p >= hi {
            break;
        }
        let state = evaluate(p)?;
        evaluated = i;
        residual = state.specific_volume().value() - volume;
        if residual.abs() <= controls.volume_tolerance_m3_kg {
            poll(cancelled)?;
            return Ok(accepted(state, i));
        }
        if residual > 0.0 {
            lower = state;
        } else {
            upper = state;
        }
        let newton = p - residual / state.dv_dp_at_temperature().value();
        p = if newton > lower.pressure().value() && newton < upper.pressure().value() && newton != p
        {
            newton
        } else {
            lower.pressure().value() + 0.5 * (upper.pressure().value() - lower.pressure().value())
        };
    }
    for s in [lower, upper] {
        let r = s.specific_volume().value() - volume;
        if r.abs() < residual.abs() {
            residual = r;
        }
        if r.abs() <= controls.volume_tolerance_m3_kg {
            poll(cancelled)?;
            return Ok(accepted(s, evaluated));
        }
    }
    Err(IsochoricWaterError::VolumeNotConverged {
        iterations: evaluated,
        residual_m3_kg: residual,
    })
}

fn coexistence_point(
    volume: f64,
    controls: IsochoricWaterControls,
    temperature: f64,
    pair: If97FluidSaturationPair,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Point, IsochoricWaterError> {
    let (l, v) = (pair.liquid(), pair.vapor());
    let (vl, vv, ps) = (
        l.specific_volume().value(),
        v.specific_volume().value(),
        pair.pressure().value(),
    );
    if ps < controls.minimum_pressure_pa || ps > controls.maximum_pressure_pa {
        return Err(IsochoricWaterError::OutsideDomain);
    }
    if volume < vl || volume > vv {
        return Err(IsochoricWaterError::VolumeOutsideDomain {
            temperature_k: temperature,
            requested_m3_kg: volume,
            admitted_m3_kg: [vl, vv],
        });
    }
    let quality = if volume == vl {
        0.0
    } else if volume == vv {
        1.0
    } else {
        (volume - vl) / (vv - vl)
    };
    if !quality.is_finite()
        || !(0.0..=1.0).contains(&quality)
        || (volume > vl && volume < vv && (quality == 0.0 || quality == 1.0))
    {
        return Err(IsochoricWaterError::Unrepresentable("saturation quality"));
    }
    let slope = saturation_pressure_derivative(Temperature::new(temperature))?.value();
    let along = |s: If97FluidState| {
        let vt = s.dv_dt_at_pressure().value();
        let vp = s.dv_dp_at_temperature().value();
        // Differentiate ORIGINAL u(T,p) at the phase's actual EOS pressure.
        // Dense root residuals do not authorize replacing p by ps in
        // du = Cp*dT - p*dv - T*(dv/dT)_p*dp.
        let p = s.pressure().value();
        (
            vt + vp * slope,
            s.cp().value() - p * vt + (-temperature * vt - p * vp) * slope,
        )
    };
    let (vl_t, ul_t) = along(l);
    let (vv_t, uv_t) = along(v);
    let quality_t = -(vl_t + quality * (vv_t - vl_t)) / (vv - vl);
    let cv = ul_t
        + quality * (uv_t - ul_t)
        + quality_t * (v.specific_internal_energy().value() - l.specific_internal_energy().value());
    let source = Source::Mixture { pair, quality };
    if !cv.is_finite()
        || cv <= 0.0
        || (source.volume() - volume).abs() > controls.volume_tolerance_m3_kg
    {
        return Err(IsochoricWaterError::Unrepresentable(
            "stable coexistence isochore or represented volume",
        ));
    }
    poll(cancelled)?;
    Ok(Point {
        source,
        cv_j_kg_k: cv,
        pressure_iterations: l.iterations() + v.iterations(),
    })
}

#[cfg(test)]
mod regional_tests;
#[cfg(test)]
mod tests;
