//! Spatial self-heating followed by actual scheduled, stress-free tip births.
//! The thermal timeline owns event order; the original bridge owns tip geometry.
//! No birth or history update occurs during nonlinear trials. A failed final
//! birth discards the interval and its complete event batch. This bounded lane
//! inherits the existing single-layer, clamped-seed, perfect-end-bond model.

use core::convert::Infallible;
use fs_solid::bridge::growth::{BirthReport, CellBirth};
use fs_thermal::enthalpy::{CellId, DepositionObserver, Event, StepConfig};
use super::*;

/// Finite electrical storage and material births accepted in one transaction.
pub mod electrical;

/// One future thermal birth and its explicit unstrained projected axial length.
#[derive(Clone, Copy, Debug)]
pub struct TipBirth<'m> {
    pub binding: BridgeThermalBinding<'m>,
    pub length_m: f64,
}

#[derive(Debug)]
pub enum GrowthError<E> {
    Interval(BridgeEnthalpyError<E>),
    Birth(BridgeEnthalpyError<Infallible>),
}
impl<E: core::fmt::Debug> core::fmt::Display for GrowthError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "growing bridge/enthalpy: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for GrowthError<E> {}

/// Zero-duration event accounts, separate from the preceding solved interval.
#[derive(Debug)]
pub struct BirthAudit {
    pub thermal: EnergyAudit,
    /// Actual append reports in the thermal schedule's order, never ID order.
    pub mechanical: Vec<BirthReport>,
    pub energy_residual_j: f64,
}
pub struct GrowthReport<S> {
    pub interval: BridgeEnthalpyReport<S>,
    pub births: BirthAudit,
    /// A new tip has not yet equilibrated under the enlarged load field.
    pub birth_equilibrium_pending: bool,
    pub total_energy_residual_j: f64,
}

// Record events only. Appends run afterward with the caller's cancellation gate,
// rather than hiding an entire beam clone in an uninterruptible observer event.
#[derive(Clone, Default)]
struct BirthEvents { ids: Vec<(CellId, u64)> }
impl DepositionObserver for BirthEvents {
    type Error = DepositionError;
    fn begin(&mut self) -> Result<(), Self::Error> { self.ids.clear(); Ok(()) }
    fn event(&mut self, event: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        if let Event::Birth(cell) = event { self.ids.push((cell.id, cell.body_id)); }
        Ok(())
    }
    fn step(&mut self, _: f64, _: f64, _: &BTreeMap<CellId,f64>,
        _: &DepositionSimulation<'_>, _: &mut dyn FnMut()->Result<(),DepositionError>)
        -> Result<(),Self::Error> {
        Err(DepositionError::InvalidInput("birth synchronization cannot advance time"))
    }
}
fn poll(cancelled: &mut impl FnMut()->bool) -> Result<(),DepositionError> {
    if cancelled() { Err(DepositionError::Cancelled) } else { Ok(()) }
}

/// An accepted beam, heat field and bounded future material registry.
/// Requires an already born single-layer clamped seed; arbitrary nozzle motion,
/// moving contact geometry, multi-layer tip growth and impact are not inferred.
#[derive(Clone)]
pub struct GrowingBridge<'t,'m> {
    thermal: DepositionSimulation<'t>,
    bridge: solid::Bridge<'m>,
    bindings: Vec<BridgeThermalBinding<'m>>,
    future: BTreeMap<CellId,TipBirth<'m>>,
    body: u64,
    pending: bool,
}
impl<'t,'m> GrowingBridge<'t,'m> {
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(thermal:DepositionSimulation<'t>, bridge:solid::Bridge<'m>,
        body:u64, bindings:Vec<BridgeThermalBinding<'m>>, future:Vec<TipBirth<'m>>,
        maximum_registered_cells:usize, config:solid::Config, event_energy_tolerance_j:f64,
        mut cancelled:impl FnMut()->bool)
        -> Result<(Self,BirthAudit),BridgeEnthalpyError<Infallible>> {
        poll(&mut cancelled)?;
        let count=bindings.len().checked_add(future.len()).ok_or_else(||invalid("registry size overflow"))?;
        if maximum_registered_cells==0 || count>maximum_registered_cells {
            return Err(invalid("spatial birth registry exceeds its explicit bound"));
        }
        // The original growth owner admits the seed and small-rotation controls.
        bridge.append_cantilever_controlled(&[],config,&mut cancelled)?;
        let mut ids=BTreeSet::new();
        for b in &bindings {
            poll(&mut cancelled)?;
            if !ids.insert(b.axial_cell) { return Err(invalid("tip growth needs one layer per axial cell")); }
        }
        let mut registered=BTreeMap::new();
        for birth in future {
            poll(&mut cancelled)?;
            if !birth.length_m.is_finite() || birth.length_m<=0.0
                || !ids.insert(birth.binding.axial_cell)
                || thermal.cells().contains_key(&birth.binding.section.cell)
                || registered.insert(birth.binding.section.cell,birth).is_some() {
                return Err(invalid("future axial/thermal IDs must be new and lengths positive"));
            }
        }
        let state=Self {thermal,bridge,bindings,future:registered,body,pending:true};
        state.admit(&mut cancelled)?;
        state.synchronize(config,event_energy_tolerance_j,&mut cancelled)
    }
    pub fn thermal(&self)->&DepositionSimulation<'t> { &self.thermal }
    pub fn bridge(&self)->&solid::Bridge<'m> { &self.bridge }
    pub fn bindings(&self)->&[BridgeThermalBinding<'m>] { &self.bindings }
    pub fn remaining_births(&self)->usize { self.future.len() }
    pub fn birth_equilibrium_pending(&self)->bool { self.pending }

    fn admit(&self, cancelled:&mut impl FnMut()->bool)->Result<(),BridgeEnthalpyError<Infallible>> {
        super::admit_bindings::<Infallible>(&self.thermal,&self.bridge,&self.bindings,||poll(cancelled))?;
        let mut cells=BTreeSet::new();
        for binding in &self.bindings {
            poll(cancelled)?;
            if self.thermal.cells()[&binding.section.cell].body_id!=self.body {
                return Err(invalid("bound spatial material belongs to a different body"));
            }
            cells.insert(binding.section.cell);
        }
        for (id,cell) in self.thermal.cells() {
            poll(cancelled)?;
            if cell.body_id==self.body && !cells.contains(id) {
                return Err(invalid("selected thermal body contains unrepresented bridge material"));
            }
        }
        poll(cancelled)?;
        Ok(())
    }

    fn synchronize(&self, config:solid::Config, tolerance:f64, cancelled:&mut impl FnMut()->bool)
        -> Result<(Self,BirthAudit),BridgeEnthalpyError<Infallible>> {
        poll(cancelled)?;
        let (thermal,events,audit)=self.thermal.advance_observed_limited(&BirthEvents::default(),
            self.thermal.time_s(),StepConfig {max_dt_s:1.0,energy_tolerance_j:tolerance},0,
            ||poll(cancelled))?;
        let mut next=self.clone(); next.thermal=thermal;
        let mut reports=Vec::new();
        for (id,body) in events.ids {
            poll(cancelled)?;
            if body!=self.body {
                if next.future.contains_key(&id) { return Err(invalid("registered tip birth has the wrong body")); }
                continue;
            }
            let birth=next.future.remove(&id).ok_or_else(||invalid("spatial thermal birth lacks a tip binding"))?;
            let cell=&next.thermal.cells()[&id];
            let delivered=CellBirth {id:birth.binding.axial_cell,length_m:birth.length_m,
                mass_kg:cell.mass_kg,temperature_k:next.thermal.temperature_k(id)?,
                law:birth.binding.section.material.law()};
            let (bridge,report)=next.bridge.append_cantilever_controlled(
                std::slice::from_ref(&delivered),config,cancelled)?;
            next.bridge=bridge;
            next.bindings.push(birth.binding);
            next.pending=true;
            reports.push(report);
        }
        next.admit(cancelled)?;
        let residual=finite(next.bridge.stored_energy_j()-self.bridge.stored_energy_j()
            +audit.stored_energy_change_j-audit.incoming_energy_j-audit.external_heat_j)?;
        if residual.abs()>tolerance { return Err(invalid("spatial birth storage and incoming enthalpy do not close")); }
        poll(cancelled)?;
        Ok((next,BirthAudit {thermal:audit,mechanical:reports,energy_residual_j:residual}))
    }

    fn finish<S,E>(&self, candidate:BridgeEnthalpyCandidate<'t,'m,S>,config:BridgeEnthalpyConfig,
        cancelled:&mut impl FnMut()->bool)->Result<(Self,GrowthReport<S>),GrowthError<E>> {
        let (thermal,bridge,interval)=candidate.into_parts();
        let staged=Self {thermal,bridge,bindings:self.bindings.clone(),future:self.future.clone(),
            body:self.body,pending:false};
        let mechanics=solid::Config {section:config.feedback.mechanics,maximum_deflection_m:config.maximum_deflection_m};
        let (next,births)=staged.synchronize(mechanics,config.feedback.thermal_energy_tolerance_j,cancelled)
            .map_err(GrowthError::Birth)?;
        let total=next.bridge.stored_energy_j()-self.bridge.stored_energy_j()
            +next.thermal.energy_j()-self.thermal.energy_j()
            +interval.mechanical.energy_j.numerical_dissipation-interval.mechanical.load_work_j
            -interval.non_mechanical_heat_j-births.thermal.incoming_energy_j-births.thermal.external_heat_j;
        if !total.is_finite() || total.abs()>config.feedback.total_energy_tolerance_j {
            return Err(GrowthError::Birth(invalid("complete spatial interval/birth energy does not close")));
        }
        poll(cancelled).map_err(|e|GrowthError::Birth(e.into()))?;
        let pending=next.pending;
        Ok((next,GrowthReport {interval,births,birth_equilibrium_pending:pending,total_energy_residual_j:total}))
    }

    /// Caller splits at scheduled events and load discontinuities. All old-cell
    /// loads/energy guesses describe the pre-event interval. A failed last birth
    /// returns neither the new heat field nor the new mechanical history.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate<S,E>(&self,target_s:f64,gravity_m_s2:f64,loads_n_m:&[f64],
        config:BridgeEnthalpyConfig,guess_j:&[f64],
        resolve:impl Fn(&DepositionSimulation<'t>)->Result<ThermalCoefficients<S>,E>,
        mut cancelled:impl FnMut()->bool)->Result<(Self,GrowthReport<S>),GrowthError<E>> {
        let candidate=super::step_candidate(&self.thermal,&self.bridge,&self.bindings,target_s,
            gravity_m_s2,loads_n_m,config,guess_j,resolve,&mut cancelled).map_err(GrowthError::Interval)?;
        self.finish(candidate,config,&mut cancelled)
    }

    /// Source-pinned conduction through the interval followed by the same atomic
    /// births. Use PreparedMaterialNetwork::extend_thermal on the accepted new
    /// topology to keep old source selections; no contacts are inferred here.
    #[allow(clippy::too_many_arguments)]
    pub fn step_transport_candidate(&self,transport:&PreparedMaterialNetwork<'_,'_>,
        target_s:f64,gravity_m_s2:f64,loads_n_m:&[f64],config:BridgeEnthalpyConfig,guess_j:&[f64],
        mut cancelled:impl FnMut()->bool)
        -> Result<(Self,GrowthReport<MaterialNetworkSample>),GrowthError<TransportError>> {
        let candidate=super::step_transport_candidate(&self.thermal,&self.bridge,&self.bindings,transport,
            target_s,gravity_m_s2,loads_n_m,config,guess_j,&mut cancelled).map_err(GrowthError::Interval)?;
        self.finish(candidate,config,&mut cancelled)
    }
}
