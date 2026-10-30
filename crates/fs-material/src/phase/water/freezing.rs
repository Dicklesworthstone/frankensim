//! Pressure-controlled freezing and thawing on original ice and liquid sources.
//!
//! The IAPWS Ice Ih Gibbs equation and IF97 Region 1 retain their original
//! enthalpy, internal-energy and entropy references. The released melting curve
//! selects the coexistence temperature. Only their admitted overlap is used:
//! melting must lie in 273.15..=273.16 K, so this is an ambient-pressure
//! ice/liquid owner, not a high-pressure-ice or supercooled-liquid model.
//!
//! Liquid fraction is a mass fraction. Phase volumes use the actual distinct
//! densities, including expansion on freezing. Heat changes H; the separate
//! pressure load exchanges p*Delta V. This is not a sealed rigid vessel.
//! The isobar-local entropy primitive integrates dh/T with its original fusion
//! enthalpy. Its explicitly reported source matching correction must never be
//! substituted for absolute entropy in pressure-changing equipment.

use fs_qty::{Pressure, Temperature};
use fs_thermochem::water::ice::{
    IceIhError, IceIhState, IceMeltingSolution, IceMeltingSolve, MINIMUM_TEMPERATURE_K,
    TRIPLE_POINT_PRESSURE_PA, isobaric_difference as ice_difference, melting_pressure,
    melting_temperature,
};
use fs_thermochem::water::if97::increments::isobaric_difference as liquid_difference;
use fs_thermochem::water::if97::{If97Error, If97Region, If97State};

/// Explicit source matching and bounded caloric/melting work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreezingWaterControls {
    pub minimum_temperature_k: f64,
    pub maximum_temperature_k: f64,
    pub enthalpy_tolerance_j_kg: f64,
    pub maximum_iterations: usize,
    pub maximum_entropy_seam_j_kg_k: f64,
    pub maximum_gibbs_mismatch_j_kg: f64,
    pub melting: IceMeltingSolve,
}

#[derive(Clone, Debug, PartialEq)]
pub enum FreezingWaterError {
    IceSource(IceIhError),
    LiquidSource(If97Error),
    Invalid(&'static str),
    OutsideDomain,
    AmbiguousMelting,
    ChartMismatch,
    EntropySeam {
        defect_j_kg_k: f64,
        maximum_j_kg_k: f64,
    },
    GibbsMismatch {
        defect_j_kg: f64,
        maximum_j_kg: f64,
    },
    NotConverged {
        iterations: usize,
        residual_j_kg: f64,
    },
    EnergyBalance {
        maximum_residual_j: f64,
        tolerance_j: f64,
    },
    Unrepresentable(&'static str),
    Cancelled,
}
impl From<IceIhError> for FreezingWaterError {
    fn from(error: IceIhError) -> Self {
        if error == IceIhError::Cancelled {
            Self::Cancelled
        } else {
            Self::IceSource(error)
        }
    }
}
impl From<If97Error> for FreezingWaterError {
    fn from(error: If97Error) -> Self {
        Self::LiquidSource(error)
    }
}
impl core::fmt::Display for FreezingWaterError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "ice/liquid-water isobar: {self:?}")
    }
}
impl std::error::Error for FreezingWaterError {}

/// Retained original source pair; no Gibbs or entropy datum reconciliation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreezingMeltingPair {
    melting: IceMeltingSolution,
    ice: IceIhState,
    liquid: If97State,
    entropy_defect_j_kg_k: f64,
    gibbs_mismatch_j_kg: f64,
}
impl FreezingMeltingPair {
    pub const fn ice(&self) -> &IceIhState {
        &self.ice
    }
    pub const fn liquid(&self) -> &If97State {
        &self.liquid
    }
    pub fn temperature(&self) -> Temperature {
        self.ice.temperature()
    }
    /// Requested mechanical pressure shared by both original Gibbs sources.
    pub fn pressure(&self) -> Pressure {
        self.ice.pressure()
    }
    pub const fn melting_solution(&self) -> &IceMeltingSolution {
        &self.melting
    }
    pub fn latent_enthalpy_j_kg(&self) -> f64 {
        self.liquid.specific_enthalpy().value() - self.ice.specific_enthalpy().value()
    }
    pub fn latent_entropy_j_kg_k(&self) -> f64 {
        self.liquid.specific_entropy().value() - self.ice.specific_entropy().value()
    }
    /// L/Tm minus the original source entropy difference; not entropy production.
    pub const fn entropy_defect_j_kg_k(&self) -> f64 {
        self.entropy_defect_j_kg_k
    }
    pub const fn gibbs_mismatch_j_kg(&self) -> f64 {
        self.gibbs_mismatch_j_kg
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FreezingWaterPhase {
    Ice,
    MeltingIce,
    IceLiquid,
    MeltingLiquid,
    Liquid,
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum Point {
    Ice(IceIhState),
    Melting { liquid_fraction: f64 },
    Liquid(If97State),
}

/// Immutable pressure chart whose temperature window brackets melting.
/// Whole-interval thermal stability bounds are prepared by the thermal adapter.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreezingWater {
    controls: FreezingWaterControls,
    melting: FreezingMeltingPair,
    minimum: IceIhState,
    maximum: If97State,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreezingWaterState {
    chart: FreezingWater,
    point: Point,
}
impl FreezingWaterState {
    pub const fn law(&self) -> &FreezingWater {
        &self.chart
    }
    pub fn pressure_pa(self) -> f64 {
        self.chart.pressure_pa()
    }
    pub fn temperature_k(self) -> f64 {
        match self.point {
            Point::Ice(s) => s.temperature().value(),
            Point::Liquid(s) => s.temperature().value(),
            Point::Melting { .. } => self.chart.melting.temperature().value(),
        }
    }
    pub const fn ice_source(self) -> Option<IceIhState> {
        match self.point {
            Point::Ice(source) => Some(source),
            _ => None,
        }
    }
    pub const fn liquid_source(self) -> Option<If97State> {
        match self.point {
            Point::Liquid(source) => Some(source),
            _ => None,
        }
    }
    fn properties(self) -> [f64; 4] {
        match self.point {
            Point::Ice(s) => ice_properties(s),
            Point::Liquid(s) => liquid_properties(s),
            Point::Melting { liquid_fraction: x } => {
                let a = ice_properties(self.chart.melting.ice);
                let b = liquid_properties(self.chart.melting.liquid);
                if x == 0.0 {
                    a
                } else if x == 1.0 {
                    b
                } else {
                    core::array::from_fn(|i| a[i] + x * (b[i] - a[i]))
                }
            }
        }
    }
    pub fn specific_enthalpy_j_kg(self) -> f64 {
        self.properties()[0]
    }
    pub fn specific_internal_energy_j_kg(self) -> f64 {
        self.properties()[1]
    }
    pub fn specific_volume_m3_kg(self) -> f64 {
        self.properties()[2]
    }
    /// Original source entropy, without the isobar-local fusion correction.
    pub fn specific_entropy_j_kg_k(self) -> f64 {
        self.properties()[3]
    }
    pub fn bulk_density_kg_m3(self) -> f64 {
        1.0 / self.specific_volume_m3_kg()
    }
    pub const fn liquid_mass_fraction(self) -> f64 {
        match self.point {
            Point::Ice(_) => 0.0,
            Point::Melting { liquid_fraction } => liquid_fraction,
            Point::Liquid(_) => 1.0,
        }
    }
    pub fn ice_mass_fraction(self) -> f64 {
        1.0 - self.liquid_mass_fraction()
    }
    pub fn liquid_volume_fraction(self) -> f64 {
        match self.point {
            Point::Melting { liquid_fraction } => {
                liquid_fraction * self.chart.melting.liquid.specific_volume().value()
                    / self.specific_volume_m3_kg()
            }
            _ => self.liquid_mass_fraction(),
        }
    }
    pub fn ice_volume_fraction(self) -> f64 {
        match self.point {
            Point::Melting { liquid_fraction } => {
                (1.0 - liquid_fraction) * self.chart.melting.ice.specific_volume().value()
                    / self.specific_volume_m3_kg()
            }
            _ => self.ice_mass_fraction(),
        }
    }
    pub const fn phase(self) -> FreezingWaterPhase {
        match self.point {
            Point::Ice(_) => FreezingWaterPhase::Ice,
            Point::Liquid(_) => FreezingWaterPhase::Liquid,
            Point::Melting {
                liquid_fraction: 0.0,
            } => FreezingWaterPhase::MeltingIce,
            Point::Melting {
                liquid_fraction: 1.0,
            } => FreezingWaterPhase::MeltingLiquid,
            Point::Melting { .. } => FreezingWaterPhase::IceLiquid,
        }
    }
    pub fn thermal_entropy_correction_j_kg_k(self) -> f64 {
        self.liquid_mass_fraction() * self.chart.melting.entropy_defect_j_kg_k
    }
    /// Pressure-local dh/T primitive, not an absolute entropy for pressure devices.
    pub fn thermal_specific_entropy_j_kg_k(self) -> f64 {
        self.specific_entropy_j_kg_k() + self.thermal_entropy_correction_j_kg_k()
    }
    pub fn pressure_closure_defect_j_kg(self) -> f64 {
        self.specific_enthalpy_j_kg()
            - self.specific_internal_energy_j_kg()
            - self.pressure_pa() * self.specific_volume_m3_kg()
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreezingWaterEnthalpySolution {
    state: FreezingWaterState,
    iterations: usize,
    residual_j_kg: f64,
}
impl FreezingWaterEnthalpySolution {
    pub const fn state(&self) -> &FreezingWaterState {
        &self.state
    }
    pub const fn into_state(self) -> FreezingWaterState {
        self.state
    }
    pub const fn iterations(self) -> usize {
        self.iterations
    }
    pub const fn enthalpy_residual_j_kg(self) -> f64 {
        self.residual_j_kg
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), FreezingWaterError> {
    if cancelled() {
        Err(FreezingWaterError::Cancelled)
    } else {
        Ok(())
    }
}
fn ice_properties(s: IceIhState) -> [f64; 4] {
    [
        s.specific_enthalpy().value(),
        s.specific_internal_energy().value(),
        s.specific_volume().value(),
        s.specific_entropy().value(),
    ]
}
fn liquid_properties(s: If97State) -> [f64; 4] {
    [
        s.specific_enthalpy().value(),
        s.specific_internal_energy().value(),
        s.specific_volume().value(),
        s.specific_entropy().value(),
    ]
}

impl FreezingWater {
    /// Admit the original ice/liquid overlap without extrapolating IF97 below
    /// 273.15 K. Tmin<=Tm<=Tmax; equal endpoints admit a latent-only chart.
    pub fn new(
        pressure: Pressure,
        controls: FreezingWaterControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, FreezingWaterError> {
        poll(&mut cancelled)?;
        if !controls.minimum_temperature_k.is_finite()
            || !controls.maximum_temperature_k.is_finite()
            || controls.minimum_temperature_k < MINIMUM_TEMPERATURE_K
            || controls.minimum_temperature_k > controls.maximum_temperature_k
            || !controls.enthalpy_tolerance_j_kg.is_finite()
            || controls.enthalpy_tolerance_j_kg < 0.0
            || !(1..=128).contains(&controls.maximum_iterations)
            || !controls.maximum_entropy_seam_j_kg_k.is_finite()
            || controls.maximum_entropy_seam_j_kg_k < 0.0
            || !controls.maximum_gibbs_mismatch_j_kg.is_finite()
            || controls.maximum_gibbs_mismatch_j_kg < 0.0
        {
            return Err(FreezingWaterError::Invalid(
                "ordered source temperatures, nonnegative matching/inverse budgets and 1..=128 iterations required",
            ));
        }
        let p = pressure.value();
        let maximum_pressure = melting_pressure(Temperature::new(273.15))?.value();
        if !p.is_finite() || p < TRIPLE_POINT_PRESSURE_PA || p > maximum_pressure {
            return Err(FreezingWaterError::OutsideDomain);
        }
        let melting = melting_temperature(pressure, controls.melting, &mut cancelled)?;
        let tm = melting.temperature().value();
        if !(273.15..=273.16).contains(&tm)
            || controls.minimum_temperature_k > tm
            || controls.maximum_temperature_k < tm
        {
            return Err(FreezingWaterError::OutsideDomain);
        }
        let ice = IceIhState::at_pressure_temperature(pressure, melting.temperature())?;
        let liquid = If97State::at_pressure_temperature(
            pressure,
            melting.temperature(),
            If97Region::Region1,
        )?;
        let minimum = IceIhState::at_pressure_temperature(
            pressure,
            Temperature::new(controls.minimum_temperature_k),
        )?;
        let maximum = If97State::at_pressure_temperature(
            pressure,
            Temperature::new(controls.maximum_temperature_k),
            If97Region::Region1,
        )?;
        let latent = liquid.specific_enthalpy().value() - ice.specific_enthalpy().value();
        let ds = liquid.specific_entropy().value() - ice.specific_entropy().value();
        let entropy_defect = latent / tm - ds;
        let gibbs_mismatch =
            liquid.specific_gibbs_energy().value() - ice.specific_gibbs_energy().value();
        if !latent.is_finite()
            || latent <= 0.0
            || !ds.is_finite()
            || ds <= 0.0
            || !entropy_defect.is_finite()
            || !gibbs_mismatch.is_finite()
            || ice.specific_volume().value() <= liquid.specific_volume().value()
            || (controls.minimum_temperature_k < tm
                && minimum.specific_enthalpy().value() >= ice.specific_enthalpy().value())
            || (controls.maximum_temperature_k > tm
                && maximum.specific_enthalpy().value() <= liquid.specific_enthalpy().value())
        {
            return Err(FreezingWaterError::Unrepresentable(
                "original fusion/source intervals",
            ));
        }
        if entropy_defect.abs() > controls.maximum_entropy_seam_j_kg_k {
            return Err(FreezingWaterError::EntropySeam {
                defect_j_kg_k: entropy_defect,
                maximum_j_kg_k: controls.maximum_entropy_seam_j_kg_k,
            });
        }
        if gibbs_mismatch.abs() > controls.maximum_gibbs_mismatch_j_kg {
            return Err(FreezingWaterError::GibbsMismatch {
                defect_j_kg: gibbs_mismatch,
                maximum_j_kg: controls.maximum_gibbs_mismatch_j_kg,
            });
        }
        poll(&mut cancelled)?;
        Ok(Self {
            controls,
            melting: FreezingMeltingPair {
                melting,
                ice,
                liquid,
                entropy_defect_j_kg_k: entropy_defect,
                gibbs_mismatch_j_kg: gibbs_mismatch,
            },
            minimum,
            maximum,
        })
    }
    pub const fn controls(&self) -> &FreezingWaterControls {
        &self.controls
    }
    pub const fn melting(&self) -> &FreezingMeltingPair {
        &self.melting
    }
    pub fn pressure_pa(&self) -> f64 {
        self.melting.pressure().value()
    }
    pub fn temperature_bounds_k(&self) -> [f64; 2] {
        [
            self.controls.minimum_temperature_k,
            self.controls.maximum_temperature_k,
        ]
    }
    pub fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        [
            self.minimum.specific_enthalpy().value(),
            self.maximum.specific_enthalpy().value(),
        ]
    }
    pub fn melting_entropy_defect_j_kg_k(&self) -> f64 {
        self.melting.entropy_defect_j_kg_k
    }
    fn wrap(&self, point: Point) -> FreezingWaterState {
        FreezingWaterState {
            chart: *self,
            point,
        }
    }
    fn minimum_state(&self) -> FreezingWaterState {
        if self.minimum.temperature() == self.melting.temperature() {
            self.wrap(Point::Melting {
                liquid_fraction: 0.0,
            })
        } else {
            self.wrap(Point::Ice(self.minimum))
        }
    }
    fn maximum_state(&self) -> FreezingWaterState {
        if self.maximum.temperature() == self.melting.temperature() {
            self.wrap(Point::Melting {
                liquid_fraction: 1.0,
            })
        } else {
            self.wrap(Point::Liquid(self.maximum))
        }
    }
    pub fn state_at_temperature(
        &self,
        temperature: Temperature,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FreezingWaterState, FreezingWaterError> {
        poll(&mut cancelled)?;
        let t = temperature.value();
        if !t.is_finite()
            || t < self.controls.minimum_temperature_k
            || t > self.controls.maximum_temperature_k
        {
            return Err(FreezingWaterError::OutsideDomain);
        }
        let tm = self.melting.temperature().value();
        let point = if t < tm {
            Point::Ice(IceIhState::at_pressure_temperature(
                self.melting.pressure(),
                temperature,
            )?)
        } else if t > tm {
            Point::Liquid(If97State::at_pressure_temperature(
                self.melting.pressure(),
                temperature,
                If97Region::Region1,
            )?)
        } else {
            return Err(FreezingWaterError::AmbiguousMelting);
        };
        poll(&mut cancelled)?;
        Ok(self.wrap(point))
    }
    /// Explicit liquid MASS fraction on the physical fusion plateau.
    pub fn state_at_liquid_fraction(
        &self,
        liquid_fraction: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FreezingWaterState, FreezingWaterError> {
        poll(&mut cancelled)?;
        if !liquid_fraction.is_finite() || !(0.0..=1.0).contains(&liquid_fraction) {
            return Err(FreezingWaterError::Invalid(
                "liquid mass fraction must be in [0,1]",
            ));
        }
        poll(&mut cancelled)?;
        Ok(self.wrap(Point::Melting { liquid_fraction }))
    }
    pub fn state_at_specific_enthalpy(
        &self,
        h: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FreezingWaterEnthalpySolution, FreezingWaterError> {
        self.invert(h, self.controls.enthalpy_tolerance_j_kg, &mut cancelled)
    }
    fn invert(
        &self,
        target: f64,
        tolerance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FreezingWaterEnthalpySolution, FreezingWaterError> {
        poll(cancelled)?;
        let (minimum, maximum) = (self.minimum_state(), self.maximum_state());
        let (ice, liquid) = (
            self.wrap(Point::Melting {
                liquid_fraction: 0.0,
            }),
            self.wrap(Point::Melting {
                liquid_fraction: 1.0,
            }),
        );
        let h = |s: FreezingWaterState| s.specific_enthalpy_j_kg();
        if !target.is_finite() || target < h(minimum) || target > h(maximum) {
            return Err(FreezingWaterError::OutsideDomain);
        }
        let solved = |state, iterations, residual_j_kg| FreezingWaterEnthalpySolution {
            state,
            iterations,
            residual_j_kg,
        };
        for state in [minimum, ice, liquid, maximum] {
            if target == h(state) {
                poll(cancelled)?;
                return Ok(solved(state, 0, 0.0));
            }
        }
        if target > h(ice) && target < h(liquid) {
            let x = (target - h(ice)) / self.melting.latent_enthalpy_j_kg();
            if !x.is_finite() || x <= 0.0 || x >= 1.0 {
                return Err(FreezingWaterError::Unrepresentable("fusion inventory"));
            }
            let state = self.wrap(Point::Melting { liquid_fraction: x });
            let residual = h(state) - target;
            if residual.abs() > tolerance {
                return Err(FreezingWaterError::NotConverged {
                    iterations: 0,
                    residual_j_kg: residual,
                });
            }
            poll(cancelled)?;
            return Ok(solved(state, 0, residual));
        }
        let ice_branch = target < h(ice);
        let (mut a, mut b) = if ice_branch {
            (minimum, ice)
        } else {
            (liquid, maximum)
        };
        let mut t = a.temperature_k()
            + (b.temperature_k() - a.temperature_k()) * ((target - h(a)) / (h(b) - h(a)));
        let (mut best, mut used) = (f64::INFINITY, 0);
        for iteration in 1..=self.controls.maximum_iterations {
            poll(cancelled)?;
            if t <= a.temperature_k() || t >= b.temperature_k() {
                t = a.temperature_k() + 0.5 * (b.temperature_k() - a.temperature_k());
            }
            if t <= a.temperature_k() || t >= b.temperature_k() {
                break;
            }
            let (point, derivative) = if ice_branch {
                let s = IceIhState::at_pressure_temperature(
                    self.melting.pressure(),
                    Temperature::new(t),
                )?;
                (Point::Ice(s), s.cp().value())
            } else {
                let s = If97State::at_pressure_temperature(
                    self.melting.pressure(),
                    Temperature::new(t),
                    If97Region::Region1,
                )?;
                (Point::Liquid(s), s.cp().value())
            };
            let state = self.wrap(point);
            let residual = h(state) - target;
            used = iteration;
            if residual.abs() < best.abs() {
                best = residual;
            }
            if residual.abs() <= tolerance {
                poll(cancelled)?;
                return Ok(solved(state, used, residual));
            }
            if !residual.is_finite() || !derivative.is_finite() || derivative <= 0.0 {
                return Err(FreezingWaterError::Unrepresentable("caloric inverse"));
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
            let residual = h(state) - target;
            if residual.abs() < best.abs() {
                best = residual;
            }
            if residual.abs() <= tolerance {
                poll(cancelled)?;
                return Ok(solved(state, used, residual));
            }
        }
        Err(FreezingWaterError::NotConverged {
            iterations: used,
            residual_j_kg: best,
        })
    }

    fn changes(
        &self,
        a: FreezingWaterState,
        b: FreezingWaterState,
    ) -> Result<Changes, FreezingWaterError> {
        let tm = self.melting.temperature().value();
        let ice = ice_difference(
            self.melting.pressure(),
            Temperature::new(a.temperature_k().min(tm)),
            Temperature::new(b.temperature_k().min(tm)),
        )?;
        let liquid = liquid_difference(
            self.melting.pressure(),
            If97Region::Region1,
            Temperature::new(a.temperature_k().max(tm)),
            Temperature::new(b.temperature_k().max(tm)),
        )?;
        let ice = [
            ice.specific_enthalpy_change().value(),
            ice.specific_internal_energy_change().value(),
            ice.specific_volume_change().value(),
            ice.specific_entropy_change().value(),
        ];
        let liquid = [
            liquid.specific_enthalpy_change().value(),
            liquid.specific_internal_energy_change().value(),
            liquid.specific_volume_change().value(),
            liquid.specific_entropy_change().value(),
        ];
        let x = b.liquid_mass_fraction() - a.liquid_mass_fraction();
        let (solid, fluid) = (
            ice_properties(self.melting.ice),
            liquid_properties(self.melting.liquid),
        );
        Ok(Changes {
            sensible: core::array::from_fn(|i| ice[i] + liquid[i]),
            fusion: core::array::from_fn(|i| x * (fluid[i] - solid[i])),
        })
    }
    /// Stable pressure-local integral of dh/T, including latent enthalpy once.
    /// The returned quantity uses this chart's fixed pressure and composition.
    pub fn thermal_entropy_change_j_kg_k(
        &self,
        initial: FreezingWaterState,
        final_state: FreezingWaterState,
    ) -> Result<f64, FreezingWaterError> {
        if initial.chart != *self || final_state.chart != *self {
            return Err(FreezingWaterError::ChartMismatch);
        }
        let changes = self.changes(initial, final_state)?;
        let entropy = changes.sensible[3] + changes.fusion[0] / self.melting.temperature().value();
        if entropy.is_finite() {
            Ok(entropy)
        } else {
            Err(FreezingWaterError::Unrepresentable(
                "thermal entropy change",
            ))
        }
    }
    /// ThermalLaw-facing relative entropy integral on the admitted enthalpy path.
    pub fn specific_entropy_change_j_kg_k(
        &self,
        initial_h_j_kg: f64,
        final_h_j_kg: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64, FreezingWaterError> {
        let a = self
            .state_at_specific_enthalpy(initial_h_j_kg, &mut cancelled)?
            .into_state();
        let b = self
            .state_at_specific_enthalpy(final_h_j_kg, &mut cancelled)?
            .into_state();
        let entropy = self.thermal_entropy_change_j_kg_k(a, b)?;
        poll(&mut cancelled)?;
        Ok(entropy)
    }
}

struct Changes {
    sensible: [f64; 4],
    fusion: [f64; 4],
}
impl Changes {
    fn total(&self) -> [f64; 4] {
        core::array::from_fn(|i| self.sensible[i] + self.fusion[i])
    }
}

/// Heat and physical inventory changes for an accepted fixed-mass pressure step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreezingWaterHeatStep {
    pub initial: FreezingWaterState,
    pub state: FreezingWaterState,
    pub mass_kg: f64,
    pub supplied_heat_j: f64,
    /// Original phase masses [ice, liquid].
    pub initial_phase_masses_kg: [f64; 2],
    pub phase_masses_kg: [f64; 2],
    pub volumes_m3: [f64; 2],
    pub volume_change_m3: f64,
    pub enthalpy_change_j: f64,
    pub internal_energy_change_j: f64,
    pub boundary_work_j: f64,
    pub sensible_enthalpy_change_j: f64,
    pub fusion_enthalpy_change_j: f64,
    pub source_entropy_change_j_k: f64,
    pub thermal_entropy_change_j_k: f64,
    pub entropy_seam_correction_change_j_k: f64,
    pub enthalpy_residual_j: f64,
    pub energy_residual_j: f64,
    pub represented_enthalpy_residual_j: f64,
    pub represented_energy_residual_j: f64,
    pub pressure_closure_defects_j: [f64; 2],
    pub inversion_iterations: usize,
    pub inversion_residual_j_kg: f64,
}

impl FreezingWater {
    /// Deposit heat into H and independently close Delta U+p*Delta V.
    /// Original phase energies, volumes and source matching defects are retained.
    pub fn advance_heat(
        &self,
        current: FreezingWaterState,
        mass_kg: f64,
        heat_j: f64,
        energy_tolerance_j: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FreezingWaterHeatStep, FreezingWaterError> {
        poll(&mut cancelled)?;
        if current.chart != *self {
            return Err(FreezingWaterError::ChartMismatch);
        }
        if !mass_kg.is_finite()
            || mass_kg <= 0.0
            || !heat_j.is_finite()
            || !energy_tolerance_j.is_finite()
            || energy_tolerance_j < 0.0
        {
            return Err(FreezingWaterError::Invalid(
                "positive finite mass, finite heat and nonnegative J budget required",
            ));
        }
        let h0 = current.specific_enthalpy_j_kg();
        let dh = heat_j / mass_kg;
        if !dh.is_finite() || (heat_j != 0.0 && dh == 0.0) {
            return Err(FreezingWaterError::Unrepresentable("specific heat"));
        }
        let [lo, hi] = self.enthalpy_bounds_j_kg();
        if heat_j < mass_kg * (lo - h0) || heat_j > mass_kg * (hi - h0) {
            return Err(FreezingWaterError::OutsideDomain);
        }
        let mut target = h0 + dh;
        if heat_j != 0.0 {
            for knot in [
                lo,
                hi,
                self.melting.ice.specific_enthalpy().value(),
                self.melting.liquid.specific_enthalpy().value(),
            ] {
                if heat_j == mass_kg * (knot - h0) {
                    target = knot;
                }
            }
            if !target.is_finite() || target == h0 {
                return Err(FreezingWaterError::Unrepresentable(
                    "heat lost in H coordinate",
                ));
            }
        }
        let solved = if heat_j == 0.0 {
            FreezingWaterEnthalpySolution {
                state: current,
                iterations: 0,
                residual_j_kg: 0.0,
            }
        } else {
            self.invert(
                target,
                self.controls
                    .enthalpy_tolerance_j_kg
                    .min(energy_tolerance_j / mass_kg),
                &mut cancelled,
            )?
        };
        let state = solved.state;
        let changes = self.changes(current, state)?;
        let [enthalpy, internal, volume_change, source_entropy] =
            changes.total().map(|x| mass_kg * x);
        let work = self.pressure_pa() * volume_change;
        let thermal_entropy = mass_kg
            * (changes.sensible[3] + changes.fusion[0] / self.melting.temperature().value());
        let correction = mass_kg
            * (state.liquid_mass_fraction() - current.liquid_mass_fraction())
            * self.melting.entropy_defect_j_kg_k;
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
        let closure = [current, state].map(|s| mass_kg * s.pressure_closure_defect_j_kg());
        let inventories = [current, state].map(|s| {
            [
                mass_kg * s.ice_mass_fraction(),
                mass_kg * s.liquid_mass_fraction(),
            ]
        });
        if ![
            enthalpy,
            internal,
            volume_change,
            work,
            source_entropy,
            thermal_entropy,
            correction,
        ]
        .into_iter()
        .all(f64::is_finite)
            || !volumes.into_iter().all(|v| v.is_finite() && v > 0.0)
            || !closure.into_iter().all(f64::is_finite)
            || !residuals.into_iter().all(f64::is_finite)
            || !inventories
                .into_iter()
                .flatten()
                .all(|m| m.is_finite() && m >= 0.0)
            || (heat_j != 0.0
                && (enthalpy == 0.0
                    || raw_h == 0.0
                    || (heat_j > 0.0) != (enthalpy > 0.0)
                    || (heat_j > 0.0) != (raw_h > 0.0)))
        {
            return Err(FreezingWaterError::Unrepresentable(
                "physical heat, inventory or endpoint state",
            ));
        }
        for (s, inventory) in [current, state].into_iter().zip(inventories) {
            if (s.ice_mass_fraction() > 0.0 && inventory[0] == 0.0)
                || (s.liquid_mass_fraction() > 0.0 && inventory[1] == 0.0)
            {
                return Err(FreezingWaterError::Unrepresentable(
                    "phase mass lost to representation",
                ));
            }
        }
        let maximum_residual = residuals.into_iter().map(f64::abs).fold(0.0, f64::max);
        if maximum_residual > energy_tolerance_j {
            return Err(FreezingWaterError::EnergyBalance {
                maximum_residual_j: maximum_residual,
                tolerance_j: energy_tolerance_j,
            });
        }
        poll(&mut cancelled)?;
        Ok(FreezingWaterHeatStep {
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
            fusion_enthalpy_change_j: mass_kg * changes.fusion[0],
            source_entropy_change_j_k: source_entropy,
            thermal_entropy_change_j_k: thermal_entropy,
            entropy_seam_correction_change_j_k: correction,
            enthalpy_residual_j: residuals[0],
            energy_residual_j: residuals[1],
            represented_enthalpy_residual_j: residuals[2],
            represented_energy_residual_j: residuals[3],
            pressure_closure_defects_j: closure,
            inversion_iterations: solved.iterations,
            inversion_residual_j_kg: solved.residual_j_kg,
        })
    }
}

#[cfg(test)]
mod tests;
