//! Delegate the complete native observer protocol; emit reversible power once.
use std::collections::BTreeMap;
use fs_thermal::enthalpy::{Event, HeatTransfer, collect_source_rates};
use super::*;

#[derive(Clone)]
pub(super) struct ReversibleHeat<'a,'m,'cx,O> {
    bindings:&'a Bindings<'m>,
    pub(super) inner:O,
    cx:&'a Cx<'cx>,
    watts:Vec<f64>,
    pub(super) applied_j:Vec<Dd>,
}
impl<'a,'m,'cx,O:DepositionObserver> ReversibleHeat<'a,'m,'cx,O> {
    pub(super) fn new(bindings:&'a Bindings<'m>,inner:O,heat:&[f64],dt:f64,cx:&'a Cx<'cx>)
        ->Result<Self,CutThermoelasticError3<O::Error>> {
        let mut watts=Vec::with_capacity(heat.len());
        for &q in heat {
            poll(cx)?;
            let power=finite(q/dt)?;
            if q!=0.0 && power==0.0 {
                return Err(CutThermoelasticError3::Invalid("reversible heat underflows its power conversion"));
            }
            watts.push(power);
        }
        Ok(Self {bindings,inner,cx,watts,applied_j:vec![Dd::ZERO;heat.len()]})
    }
}
impl<O: DepositionObserver> DepositionObserver for ReversibleHeat<'_, '_, '_, O> {
    type Error=CutThermoelasticError3<O::Error>;
    fn begin(&mut self)->Result<(),Self::Error> {
        self.inner.begin().map_err(CutThermoelasticError3::Observer)
    }
    fn requires_step_interval(&self)->bool {self.inner.requires_step_interval()}
    fn next_transition_time_s(&mut self,thermal:&DepositionSimulation<'_>,
        check:&mut dyn FnMut()->Result<(),DepositionError>)->Result<Option<f64>,Self::Error> {
        self.inner.next_transition_time_s(thermal,check).map_err(CutThermoelasticError3::Observer)
    }
    fn step_sources(&mut self,thermal:&DepositionSimulation<'_>,start:f64,dt:f64,end:f64,
        check:&mut dyn FnMut()->Result<(),DepositionError>)->Result<BTreeMap<CellId,f64>,Self::Error> {
        collect_source_rates(thermal, |emit| {
            self.visit_step_sources(thermal, start, dt, end, emit, check)
        })
    }
    fn sources(&mut self,thermal:&DepositionSimulation<'_>,
        checkpoint:&mut dyn FnMut()->Result<(),DepositionError>)->Result<BTreeMap<CellId,f64>,Self::Error> {
        collect_source_rates(thermal, |emit| self.visit_sources(thermal, emit, checkpoint))
    }
    fn visit_sources(&mut self, thermal: &DepositionSimulation<'_>,
        emit: &mut dyn FnMut(CellId, f64) -> Result<(), DepositionError>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>) -> Result<(), Self::Error> {
        self.bindings.admit_heat(thermal, self.cx)?;
        self.inner.visit_sources(thermal, emit, checkpoint)
            .map_err(CutThermoelasticError3::Observer)?;
        for (bound, &power) in self.bindings.cells.iter().zip(&self.watts) {
            checkpoint()?;
            if power!=0.0 {emit(bound.cell, power)?;}
        }
        Ok(())
    }
    fn visit_step_sources(&mut self, thermal: &DepositionSimulation<'_>,
        start: f64, dt: f64, end: f64,
        emit: &mut dyn FnMut(CellId, f64) -> Result<(), DepositionError>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>) -> Result<(), Self::Error> {
        self.inner.visit_step_sources(thermal, start, dt, end, emit, checkpoint)
            .map_err(CutThermoelasticError3::Observer)
    }
    fn conductances(&mut self,thermal:&DepositionSimulation<'_>,
        check:&mut dyn FnMut()->Result<(),DepositionError>)->Result<BTreeMap<u64,f64>,Self::Error> {
        self.inner.conductances(thermal,check).map_err(CutThermoelasticError3::Observer)
    }
    fn ambient_conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.inner
            .ambient_conductances(thermal, check)
            .map_err(CutThermoelasticError3::Observer)
    }
    fn reservoir_conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.inner
            .reservoir_conductances(thermal, check)
            .map_err(CutThermoelasticError3::Observer)
    }
    fn heat_transfers(&mut self,thermal:&DepositionSimulation<'_>,
        check:&mut dyn FnMut()->Result<(),DepositionError>)->Result<Vec<HeatTransfer>,Self::Error> {
        self.inner.heat_transfers(thermal,check).map_err(CutThermoelasticError3::Observer)
    }
    fn limit_timestep_s(&mut self,thermal:&DepositionSimulation<'_>,proposed:f64,
        check:&mut dyn FnMut()->Result<(),DepositionError>)->Result<f64,Self::Error> {
        self.inner.limit_timestep_s(thermal,proposed,check).map_err(CutThermoelasticError3::Observer)
    }
    fn event(&mut self,event:&Event,thermal:&DepositionSimulation<'_>)->Result<(),Self::Error> {
        if !matches!(event, Event::Boundary(_, _)) {
            return Err(CutThermoelasticError3::Invalid("cut-solid interval cannot change thermal inventory or contact topology"));
        }
        self.inner.event(event,thermal).map_err(CutThermoelasticError3::Observer)?;
        self.bindings.admit_heat(thermal,self.cx)
    }
    fn step(&mut self,start:f64,dt:f64,old:&BTreeMap<CellId,f64>,thermal:&DepositionSimulation<'_>,
        check:&mut dyn FnMut()->Result<(),DepositionError>)->Result<(),Self::Error> {
        self.inner.step(start,dt,old,thermal,check).map_err(CutThermoelasticError3::Observer)?;
        self.bindings.admit_heat(thermal,self.cx)?;
        for (i, (applied, &power)) in self.applied_j.iter_mut().zip(&self.watts).enumerate() {
            if i%64==0 {check()?;}
            *applied=*applied+Dd::from_f64(dt)*Dd::from_f64(power);
            finite(applied.to_f64())?;
        }
        Ok(())
    }
}
