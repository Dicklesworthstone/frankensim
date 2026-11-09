//! Scenario temperature instruments on the original accepted enthalpy timeline.
//!
//! The explicit support maps scenario state slots to persistent native cell IDs;
//! neither map iteration order nor a nearest-cell guess defines a probe. Each
//! native thermal substep drives the original compiled mount and sensor dynamics
//! with an affine reconstruction of its actual old/new caloric temperatures.
//! This reconstruction is an observation approximation, not an assertion that a
//! nonlinear thermal trajectory is affine. Refine the native thermal steps.
//!
//! The observer contributes no heat, conductance or timestep restriction. It can
//! be composed with physical owners through `fs_thermal::observers::CoupledObservers`.
//! Native candidate cloning publishes heat and instrument history together;
//! cancellation or a rejected physical step leaves both original states intact.

use std::collections::{BTreeMap, BTreeSet};

use fs_scenario::sensor::transient::{SensorState, SensorTransientError};
use fs_scenario::{CompiledSensorOperator, EntityId, SensorQuantity};
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionObserver, DepositionSimulation, EnergyAudit, Event,
    StepConfig,
};

/// Authored correspondence from one scenario entity to original thermal cells.
/// Slot `i` of the compiled sensor reads `cells[i]`. The geometry owner supplies
/// this correspondence and any patch weights; cell centres do not establish a
/// surface footprint or a continuum point-evaluation certificate.
#[derive(Clone, Copy, Debug)]
pub struct ThermalSensorSupport<'a> {
    /// Current entity represented by the mapped thermal field.
    pub entity: EntityId,
    /// Distinct, active persistent cell IDs in scenario state-slot order.
    pub cells: &'a [CellId],
}

/// A refusal preserves the original thermal and instrument histories.
#[derive(Clone, Debug, PartialEq)]
pub enum ThermalSensorError {
    /// Quantity, entity, support or accepted-history disagreement.
    Invalid(&'static str),
    /// The original caloric/thermal/cancellation owner refused.
    Thermal(DepositionError),
    /// The original compiled instrument or dynamic model refused.
    Sensor(SensorTransientError),
}

impl From<DepositionError> for ThermalSensorError {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}

impl core::fmt::Display for ThermalSensorError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(what) => write!(f, "thermal scenario sensor: {what}"),
            Self::Thermal(error) => write!(f, "thermal scenario sensor: {error}"),
            Self::Sensor(error) => write!(f, "thermal scenario sensor: {error}"),
        }
    }
}

impl std::error::Error for ThermalSensorError {}

/// One passive scenario instrument attached to accepted finite-cell heat.
///
/// Mutable data are deeply cloned: original cell correspondence, temperatures
/// at the last accepted endpoint and the instrument's scalar dynamic state.
/// Storage is O(the compiled state dimension), with no retained field history.
#[derive(Clone, Debug, PartialEq)]
pub struct TemperatureSensorObserver<'a> {
    cells: Vec<CellId>,
    accepted_temperatures_k: Vec<f64>,
    state: SensorState<'a>,
}

impl<'a> TemperatureSensorObserver<'a> {
    /// Bind a compiled temperature instrument at the native accepted clock.
    /// Initial instrument reading is explicit; it need not be in equilibrium
    /// for a first-order sensor. All support cells must already be active.
    pub fn new(
        operator: &'a CompiledSensorOperator,
        support: ThermalSensorSupport<'_>,
        thermal: &DepositionSimulation<'_>,
        initial_reading_k: f64,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ThermalSensorError> {
        checkpoint()?;
        if operator.quantity() != SensorQuantity::Temperature {
            return Err(ThermalSensorError::Invalid(
                "instrument must measure temperature",
            ));
        }
        if operator.entity() != support.entity {
            return Err(ThermalSensorError::Invalid(
                "sensor and thermal support entities differ",
            ));
        }
        if support.cells.len() != operator.operator().len() {
            return Err(ThermalSensorError::Invalid(
                "cell map must match the compiled state dimension",
            ));
        }
        let mut unique = BTreeSet::new();
        for (index, &id) in support.cells.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            if !unique.insert(id) {
                return Err(ThermalSensorError::Invalid(
                    "duplicate cell in sensor support",
                ));
            }
        }
        let temperatures = sample_cells(support.cells, thermal, &mut checkpoint)?;
        let mut cells = Vec::new();
        cells
            .try_reserve_exact(support.cells.len())
            .map_err(|_| ThermalSensorError::Invalid("sensor cell-map allocation refused"))?;
        cells.extend_from_slice(support.cells);
        let mut failure = None;
        let state = SensorState::new(
            operator,
            thermal.time_s(),
            &temperatures,
            initial_reading_k,
            &mut || match checkpoint() {
                Ok(()) => false,
                Err(error) => {
                    failure = Some(error);
                    true
                }
            },
        );
        if let Some(error) = failure {
            return Err(error.into());
        }
        let state = state.map_err(ThermalSensorError::Sensor)?;
        checkpoint()?;
        Ok(Self {
            cells,
            accepted_temperatures_k: temperatures,
            state,
        })
    }

    /// Original compiled instrument, accepted physical time, response and target.
    #[must_use]
    pub const fn state(&self) -> &SensorState<'a> {
        &self.state
    }

    /// Stable field-slot correspondence retained at construction.
    #[must_use]
    pub fn cells(&self) -> &[CellId] {
        &self.cells
    }

    /// Advance native heat and this instrument in one existing transaction.
    /// For additional physical owners, pass a `CoupledObservers` containing this
    /// observer to the same native `advance_observed` entry point instead.
    pub fn advance_candidate<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        config: StepConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'m>, Self, EnergyAudit), ThermalSensorError> {
        self.check_accepted(thermal, &mut checkpoint)?;
        thermal.advance_observed(self, target_s, config, checkpoint)
    }

    /// Observe an already solved native thermal candidate without solving heat
    /// again. This connects implicit and multiphysics candidates to the same
    /// instrument response used by the explicit substep observer.
    ///
    /// Both thermal states remain borrowed and unpublished. Accept the returned
    /// sensor together with the physical candidate only after this succeeds.
    /// The caller declares an affine temperature reconstruction across this
    /// single native interval and supplies the duration used by that solver.
    /// Unresolved events or strongly curved temperatures require smaller
    /// physical intervals; endpoint observation does not recover hidden steps.
    pub fn observe_endpoint_candidate(
        &self,
        accepted_thermal: &DepositionSimulation<'_>,
        endpoint_thermal: &DepositionSimulation<'_>,
        duration_s: f64,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ThermalSensorError> {
        self.check_accepted(accepted_thermal, &mut checkpoint)?;
        let (state, temperatures) =
            self.response_candidate(duration_s, endpoint_thermal, &mut checkpoint)?;
        let mut cells = Vec::new();
        cells
            .try_reserve_exact(self.cells.len())
            .map_err(|_| ThermalSensorError::Invalid("sensor cell-map allocation refused"))?;
        cells.extend_from_slice(&self.cells);
        checkpoint()?;
        Ok(Self {
            cells,
            accepted_temperatures_k: temperatures,
            state,
        })
    }

    fn response_candidate(
        &self,
        duration_s: f64,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(SensorState<'a>, Vec<f64>), ThermalSensorError> {
        checkpoint()?;
        let start_s = self.state.time_s();
        if !duration_s.is_finite()
            || duration_s <= 0.0
            || !(start_s + duration_s == thermal.time_s()
                || thermal.time_s() - start_s == duration_s)
        {
            return Err(ThermalSensorError::Invalid(
                "sensor interval differs from accepted thermal interval",
            ));
        }
        let temperatures = sample_cells(&self.cells, thermal, checkpoint)?;
        let mut failure = None;
        let candidate = self.state.step_affine_interval_candidate(
            thermal.time_s(),
            duration_s,
            &temperatures,
            &mut || match checkpoint() {
                Ok(()) => false,
                Err(error) => {
                    failure = Some(error);
                    true
                }
            },
        );
        if let Some(error) = failure {
            return Err(error.into());
        }
        let candidate = candidate.map_err(ThermalSensorError::Sensor)?;
        checkpoint()?;
        Ok((candidate, temperatures))
    }

    fn check_accepted(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), ThermalSensorError> {
        checkpoint()?;
        if self.state.time_s() != thermal.time_s() {
            return Err(ThermalSensorError::Invalid(
                "thermal and instrument clocks differ",
            ));
        }
        for (index, (&id, &old)) in self
            .cells
            .iter()
            .zip(&self.accepted_temperatures_k)
            .enumerate()
        {
            if index % 256 == 0 {
                checkpoint()?;
            }
            if thermal.temperature_k(id)? != old {
                return Err(ThermalSensorError::Invalid(
                    "thermal field differs from accepted sensor history",
                ));
            }
        }
        checkpoint()?;
        Ok(())
    }
}

impl DepositionObserver for TemperatureSensorObserver<'_> {
    type Error = ThermalSensorError;

    fn next_transition_time_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<f64>, Self::Error> {
        self.check_accepted(thermal, checkpoint)?;
        Ok(None)
    }

    fn event(&mut self, _: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        // Original births, boundary changes and contact changes do not jump
        // temperatures of already bound cells. They must not reset sensor lag.
        self.check_accepted(thermal, &mut || Ok(()))
    }

    fn step(
        &mut self,
        start_s: f64,
        duration_s: f64,
        old_temperatures_k: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        checkpoint()?;
        if self.state.time_s() != start_s {
            return Err(ThermalSensorError::Invalid(
                "sensor interval differs from accepted thermal interval",
            ));
        }
        for (index, (&id, &old)) in self
            .cells
            .iter()
            .zip(&self.accepted_temperatures_k)
            .enumerate()
        {
            if index % 256 == 0 {
                checkpoint()?;
            }
            if old_temperatures_k.get(&id) != Some(&old) {
                return Err(ThermalSensorError::Invalid(
                    "stale starting thermal field for sensor",
                ));
            }
        }
        let (candidate, temperatures) = self.response_candidate(duration_s, thermal, checkpoint)?;
        checkpoint()?;
        self.state = candidate;
        self.accepted_temperatures_k = temperatures;
        Ok(())
    }
}

fn sample_cells(
    cells: &[CellId],
    thermal: &DepositionSimulation<'_>,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<Vec<f64>, ThermalSensorError> {
    checkpoint()?;
    let mut values = Vec::new();
    values
        .try_reserve_exact(cells.len())
        .map_err(|_| ThermalSensorError::Invalid("sensor temperature allocation refused"))?;
    for (index, &id) in cells.iter().enumerate() {
        if index % 256 == 0 {
            checkpoint()?;
        }
        values.push(thermal.temperature_k(id)?);
    }
    checkpoint()?;
    Ok(values)
}
