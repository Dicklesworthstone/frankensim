//! Add the frame's physical/reversible heat to the original transport observer.
//! Every native event, interval-source and transport hook remains with its owner.

use super::{Cx, FrameEnthalpyError, ThermoplasticFrame};
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionObserver, DepositionSimulation, Event, HeatTransfer,
    collect_source_rates,
};
use std::collections::BTreeMap;

#[derive(Clone)]
pub(super) struct FrameHeat<'s, 'm, 'cx, O> {
    pub(super) frame: &'s ThermoplasticFrame<'m>,
    pub(super) inner: O,
    pub(super) cx: &'s Cx<'cx>,
    pub(super) watts: BTreeMap<CellId, f64>,
}

impl<O: DepositionObserver> DepositionObserver for FrameHeat<'_, '_, '_, O> {
    type Error = FrameEnthalpyError<O::Error>;

    fn begin(&mut self) -> Result<(), Self::Error> {
        self.inner.begin().map_err(FrameEnthalpyError::Observer)
    }

    fn requires_step_interval(&self) -> bool {
        self.inner.requires_step_interval()
    }

    fn next_transition_time_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<f64>, Self::Error> {
        self.inner
            .next_transition_time_s(thermal, check)
            .map_err(FrameEnthalpyError::Observer)
    }

    fn step_sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        start_s: f64,
        duration_s: f64,
        end_s: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        collect_source_rates(thermal, |emit| {
            self.visit_step_sources(thermal, start_s, duration_s, end_s, emit, check)
        })
    }

    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        collect_source_rates(thermal, |emit| self.visit_sources(thermal, emit, check))
    }

    fn visit_sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        emit: &mut dyn FnMut(CellId, f64) -> Result<(), DepositionError>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.frame.check_thermal(thermal, false, self.cx)?;
        self.inner
            .visit_sources(thermal, emit, check)
            .map_err(FrameEnthalpyError::Observer)?;
        for (&cell, &watts) in &self.watts {
            check()?;
            emit(cell, watts)?;
        }
        Ok(())
    }

    fn visit_step_sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        start: f64,
        dt: f64,
        end: f64,
        emit: &mut dyn FnMut(CellId, f64) -> Result<(), DepositionError>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.inner
            .visit_step_sources(thermal, start, dt, end, emit, check)
            .map_err(FrameEnthalpyError::Observer)
    }

    fn conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<u64, f64>, Self::Error> {
        self.inner
            .conductances(thermal, check)
            .map_err(FrameEnthalpyError::Observer)
    }

    fn reservoir_conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.inner
            .reservoir_conductances(thermal, check)
            .map_err(FrameEnthalpyError::Observer)
    }

    fn ambient_conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.inner
            .ambient_conductances(thermal, check)
            .map_err(FrameEnthalpyError::Observer)
    }

    fn heat_transfers(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<HeatTransfer>, Self::Error> {
        self.inner
            .heat_transfers(thermal, check)
            .map_err(FrameEnthalpyError::Observer)
    }

    fn limit_timestep_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        proposed_s: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, Self::Error> {
        self.inner
            .limit_timestep_s(thermal, proposed_s, check)
            .map_err(FrameEnthalpyError::Observer)
    }

    fn event(
        &mut self,
        event: &Event,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), Self::Error> {
        self.inner
            .event(event, thermal)
            .map_err(FrameEnthalpyError::Observer)?;
        self.frame.check_thermal(thermal, false, self.cx)
    }

    fn step(
        &mut self,
        start_s: f64,
        duration_s: f64,
        old_temperatures: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.inner
            .step(start_s, duration_s, old_temperatures, thermal, check)
            .map_err(FrameEnthalpyError::Observer)?;
        self.frame.check_thermal(thermal, false, self.cx)
    }
}
