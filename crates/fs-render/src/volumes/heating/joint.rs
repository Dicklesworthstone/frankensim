//! Incident absorption and thermal emission on one native enthalpy trajectory.
//!
//! Illumination and thermal radiation may use different wavelength packets and
//! source-supported material tables. They read the SAME native start state and
//! receive the SAME accepted substep, so conduction, absorption, reabsorption,
//! bath exchange and latent storage evolve together. Spectral weights remain
//! with each physical owner; visible source data never imply infrared support.
//!
//! The original `CoupledObservers` owns source composition, conservative
//! transfers, timestep restrictions and event ordering. This adapter introduces
//! no second thermal integrator or radiative constitutive law. It is explicit
//! coupling; time and photon refinement are needed to assess nonlinear feedback.

use super::exchange::feedback::{RadiativeFeedbackError, TemperatureRadiativeExchange};
use super::feedback::OpticalFeedbackError;
use super::illumination::TemperatureIlluminationHeating;
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionObserver, DepositionSimulation, Event, HeatTransfer,
};
use fs_thermal::observers::{CoupledObservers, ObserverCouplingError};
use std::collections::BTreeMap;

/// Complete joint trajectories with separate temporal and photon assessments.
pub mod refinement;

/// Preserve the failing illumination, radiation or native thermal owner.
pub type JointOpticalError = ObserverCouplingError<OpticalFeedbackError, RadiativeFeedbackError>;

/// A jointly evolved optical source pair with independent spectral packets.
///
/// Optical grids and mappings are supplied explicitly by each original owner;
/// they need not have the same resolution or wavelength count. Both must refer
/// to active cells of the supplied native thermal state. Consistency of the
/// underlying specimen/source data is a modeling responsibility, not inferred
/// from equal densities, equal wavelengths, or material names.
#[derive(Clone, Debug)]
pub struct JointOpticalHeating<'g, const NI: usize, const NR: usize> {
    pub(crate) observers: CoupledObservers<
        TemperatureIlluminationHeating<'g, NI>,
        TemperatureRadiativeExchange<'g, NR>,
    >,
}

impl<'g, const NI: usize, const NR: usize> JointOpticalHeating<'g, NI, NR> {
    /// Join admitted optical owners at their common accepted native clock.
    /// No photons are sampled and no thermal state is advanced by admission.
    /// Retained component histories remain with their original owners.
    pub fn try_new(
        illumination: TemperatureIlluminationHeating<'g, NI>,
        radiation: TemperatureRadiativeExchange<'g, NR>,
        thermal: &DepositionSimulation<'_>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, JointOpticalError> {
        checkpoint()?;
        if thermal.cells().is_empty()
            || illumination.time_s() != thermal.time_s()
            || radiation.time_s() != thermal.time_s()
        {
            return Err(DepositionError::InvalidInput(
                "joint optical owners require one active native thermal state and clock",
            )
            .into());
        }
        illumination
            .optical_field()
            .validate_temperatures(thermal, &mut checkpoint)
            .map_err(ObserverCouplingError::First)?;
        radiation
            .field()
            .validate_temperatures(thermal, &mut checkpoint)
            .map_err(|error| ObserverCouplingError::Second(error.into()))?;
        checkpoint()?;
        Ok(Self {
            observers: CoupledObservers::new(illumination, radiation),
        })
    }

    /// Accepted illumination history and its physical incident/escape accounts.
    pub fn illumination(&self) -> &TemperatureIlluminationHeating<'g, NI> {
        &self.observers.first
    }

    /// Accepted volume/bath exchange history and conservative internal transfers.
    pub fn radiation(&self) -> &TemperatureRadiativeExchange<'g, NR> {
        &self.observers.second
    }

    /// Common accepted time, in seconds, after a successful native advance.
    pub fn time_s(&self) -> f64 {
        self.illumination().time_s()
    }

    /// Recover the original physical owners at the accepted common endpoint.
    pub fn into_parts(
        self,
    ) -> (
        TemperatureIlluminationHeating<'g, NI>,
        TemperatureRadiativeExchange<'g, NR>,
    ) {
        self.observers.into_parts()
    }
}

impl<const NI: usize, const NR: usize> DepositionObserver for JointOpticalHeating<'_, NI, NR> {
    type Error = JointOpticalError;

    fn begin(&mut self) -> Result<(), Self::Error> {
        self.observers.begin()
    }

    fn requires_step_interval(&self) -> bool {
        self.observers.requires_step_interval()
    }

    fn next_transition_time_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<f64>, Self::Error> {
        self.observers.next_transition_time_s(thermal, checkpoint)
    }

    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.observers.sources(thermal, checkpoint)
    }

    fn step_sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        start_s: f64,
        duration_s: f64,
        end_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.observers
            .step_sources(thermal, start_s, duration_s, end_s, checkpoint)
    }

    fn conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<u64, f64>, Self::Error> {
        self.observers.conductances(thermal, checkpoint)
    }

    fn reservoir_conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.observers.reservoir_conductances(thermal, checkpoint)
    }

    fn heat_transfers(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<HeatTransfer>, Self::Error> {
        self.observers.heat_transfers(thermal, checkpoint)
    }

    fn limit_timestep_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        proposed_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, Self::Error> {
        self.observers
            .limit_timestep_s(thermal, proposed_s, checkpoint)
    }

    fn event(
        &mut self,
        event: &Event,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), Self::Error> {
        self.observers.event(event, thermal)
    }

    fn step(
        &mut self,
        start_s: f64,
        duration_s: f64,
        old: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.observers
            .step(start_s, duration_s, old, thermal, checkpoint)
    }
}
