//! Fixed-volume heating of original IF97 Region 3 dense water.
//!
//! The retained source interval proves positive Cv and original regional
//! admission throughout the temperature path. Heat changes internal energy;
//! pressure follows the original density/temperature Helmholtz equation.
//! Enthalpy changes by Delta U + V Delta p although the rigid boundary does no
//! mechanical work. Source increments preserve small changes before absolute
//! caloric subtraction, and actual rounded endpoints are independently checked.
//! No Region 1/2/3 stitching, phase mixture, wall compliance or kinetics is
//! inferred. Use the separate isochoric water owner for low-temperature boiling.

/// Rigid dense-water volumes on the existing conservative thermal timeline.
pub mod thermal;

use fs_qty::{Density, Temperature};
use fs_thermochem::water::if97::region3::increments::{
    DenseIsochoricInterval, DenseIsochoricIntervalError, DenseIsochoricIntervalSolve,
};
use fs_thermochem::water::if97::region3::{If97DenseError, If97DenseState};

/// Explicit temperature window, source-domain proof and caloric inverse work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DenseWaterControls {
    pub minimum_temperature_k: f64,
    pub maximum_temperature_k: f64,
    pub internal_energy_tolerance_j_kg: f64,
    pub maximum_iterations: usize,
    pub interval: DenseIsochoricIntervalSolve,
}

#[derive(Clone, Debug, PartialEq)]
pub enum DenseWaterError {
    Source(If97DenseError),
    Interval(DenseIsochoricIntervalError),
    Invalid(&'static str),
    OutsideDomain,
    ChartMismatch,
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
impl From<If97DenseError> for DenseWaterError {
    fn from(value: If97DenseError) -> Self {
        Self::Source(value)
    }
}
impl From<DenseIsochoricIntervalError> for DenseWaterError {
    fn from(value: DenseIsochoricIntervalError) -> Self {
        if value == DenseIsochoricIntervalError::Cancelled {
            Self::Cancelled
        } else {
            Self::Interval(value)
        }
    }
}
impl core::fmt::Display for DenseWaterError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "IF97 dense-water isochore: {self:?}")
    }
}
impl std::error::Error for DenseWaterError {}

/// One fixed-density source chart, prepared once and reused across heat steps.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DenseIsochore {
    density: Density,
    controls: DenseWaterControls,
    interval: DenseIsochoricInterval,
    endpoints: [If97DenseState; 2],
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DenseWaterState {
    chart: DenseIsochore,
    source: If97DenseState,
}
impl DenseWaterState {
    pub const fn law(&self) -> &DenseIsochore {
        &self.chart
    }
    pub const fn source(self) -> If97DenseState {
        self.source
    }
    pub fn temperature_k(self) -> f64 {
        self.source.temperature().value()
    }
    pub fn pressure_pa(self) -> f64 {
        self.source.pressure().value()
    }
    pub fn density_kg_m3(self) -> f64 {
        self.source.density().value()
    }
    pub fn specific_volume_m3_kg(self) -> f64 {
        self.source.specific_volume().value()
    }
    pub fn specific_internal_energy_j_kg(self) -> f64 {
        self.source.specific_internal_energy().value()
    }
    pub fn specific_enthalpy_j_kg(self) -> f64 {
        self.source.specific_enthalpy().value()
    }
    pub fn specific_entropy_j_kg_k(self) -> f64 {
        self.source.specific_entropy().value()
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DenseWaterSolution {
    state: DenseWaterState,
    iterations: usize,
    residual_j_kg: f64,
}
impl DenseWaterSolution {
    pub const fn state(&self) -> &DenseWaterState {
        &self.state
    }
    pub const fn into_state(self) -> DenseWaterState {
        self.state
    }
    pub const fn iterations(&self) -> usize {
        self.iterations
    }
    pub const fn internal_energy_residual_j_kg(&self) -> f64 {
        self.residual_j_kg
    }
}

/// Proposed fixed-mass, fixed-volume heat exchange; no partial state publication.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DenseWaterHeatStep {
    pub initial: DenseWaterState,
    pub state: DenseWaterState,
    pub mass_kg: f64,
    pub volume_m3: f64,
    pub supplied_heat_j: f64,
    pub internal_energy_change_j: f64,
    pub enthalpy_change_j: f64,
    /// V*Delta p changes enthalpy; it is not outward rigid-wall work.
    pub pressure_enthalpy_change_j: f64,
    pub boundary_work_j: f64,
    pub entropy_change_j_k: f64,
    pub energy_residual_j: f64,
    pub enthalpy_residual_j: f64,
    pub represented_energy_residual_j: f64,
    pub represented_enthalpy_residual_j: f64,
    pub inversion_iterations: usize,
    pub inversion_residual_j_kg: f64,
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), DenseWaterError> {
    if cancelled() {
        Err(DenseWaterError::Cancelled)
    } else {
        Ok(())
    }
}

impl DenseIsochore {
    pub fn new(
        density: Density,
        controls: DenseWaterControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DenseWaterError> {
        poll(&mut cancelled)?;
        if !controls.internal_energy_tolerance_j_kg.is_finite()
            || controls.internal_energy_tolerance_j_kg < 0.0
            || !(1..=128).contains(&controls.maximum_iterations)
            || !controls.minimum_temperature_k.is_finite()
            || !controls.maximum_temperature_k.is_finite()
            || controls.minimum_temperature_k >= controls.maximum_temperature_k
        {
            return Err(DenseWaterError::Invalid(
                "ordered finite temperatures, nonnegative energy tolerance and 1..=128 iterations required",
            ));
        }
        let interval = DenseIsochoricInterval::new(
            density,
            Temperature::new(controls.minimum_temperature_k),
            Temperature::new(controls.maximum_temperature_k),
            controls.interval,
            &mut cancelled,
        )?;
        let endpoints = [
            controls.minimum_temperature_k,
            controls.maximum_temperature_k,
        ]
        .map(|t| If97DenseState::at_density_temperature(density, Temperature::new(t)));
        let endpoints = [endpoints[0].clone()?, endpoints[1].clone()?];
        if endpoints[0].specific_internal_energy().value()
            >= endpoints[1].specific_internal_energy().value()
        {
            return Err(DenseWaterError::Unrepresentable(
                "distinct caloric endpoints",
            ));
        }
        poll(&mut cancelled)?;
        Ok(Self {
            density,
            controls,
            interval,
            endpoints,
        })
    }
    pub const fn density(&self) -> Density {
        self.density
    }
    pub const fn source_interval(&self) -> &DenseIsochoricInterval {
        &self.interval
    }
    pub fn internal_energy_bounds_j_kg(&self) -> [f64; 2] {
        self.endpoints.map(|s| s.specific_internal_energy().value())
    }
    /// Whole-path bound on dT/du, obtained from the retained positive Cv enclosure.
    pub fn maximum_temperature_slope(&self) -> f64 {
        (1.0 / self.interval.minimum_heat_capacity_j_kg_k()).next_up()
    }
    fn wrap(&self, source: If97DenseState) -> DenseWaterState {
        DenseWaterState {
            chart: *self,
            source,
        }
    }
    pub fn state_at_temperature(
        &self,
        temperature: Temperature,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DenseWaterState, DenseWaterError> {
        poll(&mut cancelled)?;
        let t = temperature.value();
        if !t.is_finite()
            || t < self.controls.minimum_temperature_k
            || t > self.controls.maximum_temperature_k
        {
            return Err(DenseWaterError::OutsideDomain);
        }
        let source = If97DenseState::at_density_temperature(self.density, temperature)?;
        poll(&mut cancelled)?;
        Ok(self.wrap(source))
    }
    pub fn state_at_specific_internal_energy(
        &self,
        target_j_kg: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DenseWaterSolution, DenseWaterError> {
        self.invert(
            target_j_kg,
            self.controls.internal_energy_tolerance_j_kg,
            &mut cancelled,
        )
    }
    fn invert(
        &self,
        target: f64,
        tolerance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DenseWaterSolution, DenseWaterError> {
        poll(cancelled)?;
        let [mut lo, mut hi] = self.endpoints;
        let u = |s: If97DenseState| s.specific_internal_energy().value();
        if !target.is_finite() || target < u(lo) || target > u(hi) {
            return Err(DenseWaterError::OutsideDomain);
        }
        for source in [lo, hi] {
            if target == u(source) {
                poll(cancelled)?;
                return Ok(DenseWaterSolution {
                    state: self.wrap(source),
                    iterations: 0,
                    residual_j_kg: 0.0,
                });
            }
        }
        let mut t = lo.temperature().value()
            + (hi.temperature().value() - lo.temperature().value())
                * ((target - u(lo)) / (u(hi) - u(lo)));
        let mut evaluations = 0;
        let mut best_residual = f64::INFINITY;
        for iteration in 1..=self.controls.maximum_iterations {
            poll(cancelled)?;
            let (lower_t, upper_t) = (lo.temperature().value(), hi.temperature().value());
            if t <= lower_t || t >= upper_t {
                t = lower_t + 0.5 * (upper_t - lower_t);
            }
            if t <= lower_t || t >= upper_t {
                break;
            }
            let source = If97DenseState::at_density_temperature(self.density, Temperature::new(t))?;
            evaluations = iteration;
            let residual = u(source) - target;
            if residual.abs() < best_residual.abs() {
                best_residual = residual;
            }
            if residual.abs() <= tolerance {
                poll(cancelled)?;
                return Ok(DenseWaterSolution {
                    state: self.wrap(source),
                    iterations: iteration,
                    residual_j_kg: residual,
                });
            }
            if residual < 0.0 {
                lo = source;
            } else {
                hi = source;
            }
            let newton = t - residual / source.cv().value();
            t = if newton > lo.temperature().value()
                && newton < hi.temperature().value()
                && newton != t
            {
                newton
            } else {
                lo.temperature().value()
                    + 0.5 * (hi.temperature().value() - lo.temperature().value())
            };
        }
        for source in [lo, hi] {
            let residual = u(source) - target;
            if residual.abs() < best_residual.abs() {
                best_residual = residual;
            }
            if residual.abs() <= tolerance {
                poll(cancelled)?;
                return Ok(DenseWaterSolution {
                    state: self.wrap(source),
                    iterations: evaluations,
                    residual_j_kg: residual,
                });
            }
        }
        Err(DenseWaterError::NotConverged {
            iterations: evaluations,
            residual_j_kg: best_residual,
        })
    }
    /// Heat against a rigid boundary; U is storage, and H retains its pV account.
    pub fn advance_heat(
        &self,
        current: DenseWaterState,
        mass_kg: f64,
        heat_j: f64,
        energy_tolerance_j: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DenseWaterHeatStep, DenseWaterError> {
        poll(&mut cancelled)?;
        if current.chart != *self {
            return Err(DenseWaterError::ChartMismatch);
        }
        if !mass_kg.is_finite()
            || mass_kg <= 0.0
            || !heat_j.is_finite()
            || !energy_tolerance_j.is_finite()
            || energy_tolerance_j < 0.0
        {
            return Err(DenseWaterError::Invalid(
                "positive mass, finite heat and nonnegative J budget required",
            ));
        }
        let u0 = current.specific_internal_energy_j_kg();
        let du = heat_j / mass_kg;
        if !du.is_finite() || (heat_j != 0.0 && du == 0.0) {
            return Err(DenseWaterError::Unrepresentable("specific heat"));
        }
        let [lo, hi] = self.internal_energy_bounds_j_kg();
        if heat_j < mass_kg * (lo - u0) || heat_j > mass_kg * (hi - u0) {
            return Err(DenseWaterError::OutsideDomain);
        }
        let mut target = u0 + du;
        if heat_j != 0.0 {
            for knot in [lo, hi] {
                if heat_j == mass_kg * (knot - u0) {
                    target = knot;
                }
            }
            if target == u0 {
                return Err(DenseWaterError::Unrepresentable(
                    "heat lost in U coordinate",
                ));
            }
        }
        let solution = if heat_j == 0.0 {
            DenseWaterSolution {
                state: current,
                iterations: 0,
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
        let d = self
            .interval
            .difference(current.source.temperature(), state.source.temperature())?;
        let internal = mass_kg * d.specific_internal_energy_change().value();
        let enthalpy = mass_kg * d.specific_enthalpy_change().value();
        let volume = mass_kg / self.density.value();
        let pressure_enthalpy = volume * d.pressure_change().value();
        let entropy = mass_kg * d.specific_entropy_change().value();
        let represented_u = mass_kg * (state.specific_internal_energy_j_kg() - u0);
        let represented_h =
            mass_kg * (state.specific_enthalpy_j_kg() - current.specific_enthalpy_j_kg());
        let represented_pv = volume * (state.pressure_pa() - current.pressure_pa());
        let residuals = [
            internal - heat_j,
            enthalpy - pressure_enthalpy - heat_j,
            represented_u - heat_j,
            represented_h - represented_pv - heat_j,
        ];
        if !volume.is_finite()
            || volume <= 0.0
            || ![internal, enthalpy, pressure_enthalpy, entropy]
                .into_iter()
                .all(f64::is_finite)
            || !residuals.into_iter().all(f64::is_finite)
            || (heat_j != 0.0 && (internal == 0.0 || represented_u == 0.0))
        {
            return Err(DenseWaterError::Unrepresentable(
                "physical heat or endpoint state",
            ));
        }
        let maximum_residual = residuals.into_iter().map(f64::abs).fold(0.0, f64::max);
        if maximum_residual > energy_tolerance_j {
            return Err(DenseWaterError::EnergyBalance {
                maximum_residual_j: maximum_residual,
                tolerance_j: energy_tolerance_j,
            });
        }
        poll(&mut cancelled)?;
        Ok(DenseWaterHeatStep {
            initial: current,
            state,
            mass_kg,
            volume_m3: volume,
            supplied_heat_j: heat_j,
            internal_energy_change_j: internal,
            enthalpy_change_j: enthalpy,
            pressure_enthalpy_change_j: pressure_enthalpy,
            boundary_work_j: 0.0,
            entropy_change_j_k: entropy,
            energy_residual_j: residuals[0],
            enthalpy_residual_j: residuals[1],
            represented_energy_residual_j: residuals[2],
            represented_enthalpy_residual_j: residuals[3],
            inversion_iterations: solution.iterations,
            inversion_residual_j_kg: solution.residual_j_kg,
        })
    }
}

#[cfg(test)]
mod tests;
