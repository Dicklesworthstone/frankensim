//! Actual contact work enters the existing multi-cell enthalpy network once.
//! The original network owns conduction, boundary losses, phase and subcycling.
//! Its observer and the contact-source refresh publish with mechanics, or not at all.
use std::collections::{BTreeMap, BTreeSet};
use super::*;
use super::feedback::MaterialFeedbackFrame;
use super::transport::{CoupledThermalStep, ImpactThermalTrial};
use super::super::MaterialContact;
use fs_contact::normal_patch::source::SourcePointContact;
use fs_material::phase::UniformEnthalpyStep;
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionObserver,
    DepositionSimulation, EnergyAudit, Event, StepConfig};

/// One reservoir is the SAME finite thermal cell, not another heat destination.
/// All reservoirs and all active cells must be bound exactly once. Several
/// mechanical contacts may already deliver their work to that one reservoir.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReservoirCell {
    pub reservoir: String,
    pub cell: CellId,
}

/// Accepted per-cell transport observations. Net exchange is reconstructed from
/// the original network's actual enthalpy change minus its integrated contact
/// source, not from a second implementation of its conduction/boundary laws.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NetworkCellHeat {
    pub cell: CellId,
    pub contact_heat_j: f64,
    pub net_exchange_heat_j: f64,
}

/// The original network audit and the existing combined contact/heat account.
#[derive(Debug)]
pub struct NetworkImpactFrame {
    pub impact: MaterialFeedbackFrame<NetworkCellHeat>,
    /// Includes contact heat as a source to the thermal subsystem. It is internal
    /// to the combined system and is subtracted in `non_contact_heat_j`.
    pub network: EnergyAudit,
    /// Authored/observer source and environment heat, excluding contact work.
    pub non_contact_heat_j: f64,
}

impl<'m> HeatedImpactSystem<'m> {
    /// Advance contact mechanics, every bound thermal cell and a thermal observer
    /// as ONE transaction. The observer may resolve original material-dependent
    /// conductances and additional sources on every actual thermal substep.
    ///
    /// This uses `step_with_material_feedback` and the original network's
    /// `advance_observed`: no pairwise splitting or second heat integrator.
    /// `resolve` is the existing unloaded contact-source callback; return None
    /// to retain frozen mechanics. Thermal contact work is held as a power over
    /// this mechanical sample and integrated on the network's actual CFL grid.
    /// It must not also be supplied by the observer or authored boundary sources.
    ///
    /// The network must have processed its initial events. Its time must equal
    /// the mechanical sample time; every cell must match the reservoir's original
    /// borrowed caloric object, mass, reference geometry and accepted enthalpy.
    /// No pending event may execute inside this operation, even at the endpoint:
    /// topology/birth/boundary events require an explicit coupled event boundary.
    /// Existing stateless phase curves can change phase; mechanical material
    /// feedback retains its separate supported-phase gates.
    ///
    /// Thermal and source time accuracy still require max_dt refinement; the
    /// conduction CFL is not a bound on nonlinear observer source stiffness.
    /// The thermal energy allowance cannot enlarge this system's heat allowance.
    /// As with the existing observer API, Clone must isolate mutable physical
    /// history. External side effects and panics cannot be rolled back.
    pub fn step_with_thermal_network<'n, O: DepositionObserver>(
        &mut self, thermal: &mut DepositionSimulation<'n>, observer: &mut O,
        bindings: &[ReservoirCell], config: StepConfig, forces: &[f64], gate: &CancelGate,
        resolve: impl FnMut(&MaterialContact, &[CaloricReservoir<'m>])
            -> Result<Option<SourcePointContact>, ImpactError>,
    ) -> Result<NetworkImpactFrame, ImpactError>
    where O::Error: core::fmt::Display {
        if config.energy_tolerance_j > self.heat_tolerance_j
            || !config.energy_tolerance_j.is_finite() || config.energy_tolerance_j < 0.0 {
            return Err(invalid("network transport cannot enlarge the contact heat allowance"));
        }
        let mut staged = None;
        let impact = self.step_with_material_feedback(forces, gate, resolve, |trial| {
            let (next, next_observer, audit, proposals) = propose(
                thermal, observer, bindings, config, &trial, gate)?;
            staged = Some((next, next_observer, audit));
            Ok(proposals)
        })?;
        // Original source, thermal and mechanical gates have all succeeded.
        // Publication is infallible; no callback runs after this point.
        let (next, next_observer, network) = staged.expect("accepted thermal callback ran once");
        *thermal = next;
        *observer = next_observer;
        Ok(NetworkImpactFrame {
            non_contact_heat_j: impact.thermal.boundary_heat_j, impact, network,
        })
    }
}

type Proposed<'n, O> = (DepositionSimulation<'n>, O, EnergyAudit,
    CoupledThermalStep<NetworkCellHeat>);
fn propose<'n, O: DepositionObserver>(
    thermal: &DepositionSimulation<'n>, observer: &O, bindings: &[ReservoirCell],
    config: StepConfig, trial: &ImpactThermalTrial<'_>, gate: &CancelGate,
) -> Result<Proposed<'n, O>, ImpactError>
where O::Error: core::fmt::Display {
    let inputs = trial.reservoirs();
    if bindings.len() != inputs.len() || thermal.cells().len() != inputs.len() {
        return Err(invalid("thermal network needs complete one-to-one reservoir/cell coverage"));
    }
    let start = (trial.sample() - 1) as f64 * trial.duration_s();
    let target = trial.sample() as f64 * trial.duration_s();
    let duration = target - start;
    if thermal.time_s() != start || !positive(duration) || !target.is_finite() {
        return Err(invalid("thermal network and mechanics do not share the accepted sample clock"));
    }
    let mut names = BTreeMap::new();
    let mut used = BTreeSet::new();
    for binding in bindings {
        if names.insert(binding.reservoir.as_str(), binding.cell).is_some()
            || !used.insert(binding.cell) {
            return Err(invalid("thermal network reservoir/cell binding is duplicated"));
        }
    }
    let mut ids = Vec::with_capacity(inputs.len());
    let mut powers = BTreeMap::new();
    for input in inputs {
        let id = *names.get(input.key).ok_or_else(||invalid("unbound thermal reservoir"))?;
        let cell = thermal.cells().get(&id).ok_or_else(||invalid("inactive bound thermal cell"))?;
        let i = input.thermal;
        if !thermal.uses_material_law(cell.material_id, i.curve)
            || cell.mass_kg != i.mass_kg || cell.geometry.volume_m3 != i.volume_m3
            || cell.geometry.surface_area_m2 != i.surface_area_m2
            || cell.specific_h_j_kg != i.initial.specific_enthalpy_j_kg() {
            return Err(invalid("network cell differs from its original caloric reservoir state or geometry"));
        }
        let power = i.internal_heat_j / duration;
        if !power.is_finite() || (i.internal_heat_j != 0.0 && power == 0.0) {
            return Err(invalid("contact heat cannot be represented on the thermal sample clock"));
        }
        powers.insert(id, power);
        ids.push(id);
    }
    let source = ContactHeating { observer: observer.clone(), powers,
        received: ids.iter().map(|&id|(id, (0.0, 0.0))).collect() };
    let checkpoint = || if gate.is_requested() { Err(DepositionError::Cancelled) } else { Ok(()) };
    let (next, source, audit) = thermal.advance_observed(&source, target, config, checkpoint)
        .map_err(|e|match e {
            NetworkError::Thermal(DepositionError::Cancelled) => ImpactError::Cancelled,
            other => ImpactError::Owner(other.to_string()),
        })?;
    let mut total_contact = 0.0;
    let mut transfers = Vec::with_capacity(inputs.len());
    for (&id, input) in ids.iter().zip(inputs) {
        let i = input.thermal;
        let h = next.cells()[&id].specific_h_j_kg;
        let state = i.curve.state_at_specific_enthalpy(h)
            .map_err(|e|ImpactError::Owner(e.to_string()))?;
        let contact = source.received[&id].0;
        let exchange = i.mass_kg * (h - i.initial.specific_enthalpy_j_kg()) - contact;
        if !exchange.is_finite() || (contact-i.internal_heat_j).abs() > config.energy_tolerance_j {
            return Err(invalid("thermal substeps did not deliver the original contact work"));
        }
        total_contact += contact;
        transfers.push(UniformEnthalpyStep { state, external_heat_j: exchange,
            energy_residual_tolerance_j: config.energy_tolerance_j,
            report: NetworkCellHeat {cell:id, contact_heat_j:contact, net_exchange_heat_j:exchange} });
    }
    let boundary_heat_j = audit.external_heat_j - total_contact;
    if !boundary_heat_j.is_finite() || audit.events != 0 || audit.incoming_mass_kg != 0.0 {
        return Err(invalid("network transport changed the coupled material topology"));
    }
    Ok((next, source.observer, audit, CoupledThermalStep { transfers, boundary_heat_j }))
}

#[derive(Clone)]
struct ContactHeating<O> {
    observer: O,
    powers: BTreeMap<CellId, f64>,
    received: BTreeMap<CellId, (f64, f64)>,
}
#[derive(Debug)]
enum NetworkError<E> { Thermal(DepositionError), Observer(E) }
impl<E> From<DepositionError> for NetworkError<E> {
    fn from(error: DepositionError) -> Self { Self::Thermal(error) }
}
impl<E: core::fmt::Display> core::fmt::Display for NetworkError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {Self::Thermal(e)=>write!(f,"thermal network: {e}"),
            Self::Observer(e)=>write!(f,"thermal observer: {e}")}
    }
}
impl<O: DepositionObserver> DepositionObserver for ContactHeating<O> {
    type Error = NetworkError<O::Error>;
    fn begin(&mut self) -> Result<(),Self::Error> {
        self.observer.begin().map_err(NetworkError::Observer)
    }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(),Self::Error> {
        Err(DepositionError::InvalidInput("process thermal events at a separate coupled boundary").into())
    }
    fn sources(&mut self, thermal:&DepositionSimulation<'_>, checkpoint:&mut dyn FnMut()->Result<(),DepositionError>)
        ->Result<BTreeMap<CellId,f64>,Self::Error> {
        let mut sources=self.observer.sources(thermal,checkpoint).map_err(NetworkError::Observer)?;
        for (&id,&power) in &self.powers {
            checkpoint()?;
            if power != 0.0 {
                let value=sources.entry(id).or_insert(0.0);
                *value+=power;
                if !value.is_finite() {return Err(DepositionError::Arithmetic.into());}
            }
        }
        Ok(sources)
    }
    fn conductances(&mut self, thermal:&DepositionSimulation<'_>, checkpoint:&mut dyn FnMut()->Result<(),DepositionError>)
        ->Result<BTreeMap<u64,f64>,Self::Error> {
        self.observer.conductances(thermal,checkpoint).map_err(NetworkError::Observer)
    }
    fn step(&mut self,start:f64,dt:f64,old:&BTreeMap<CellId,f64>,thermal:&DepositionSimulation<'_>,
        checkpoint:&mut dyn FnMut()->Result<(),DepositionError>)->Result<(),Self::Error> {
        for (&id,&power) in &self.powers {
            checkpoint()?;
            let (sum,correction)=self.received.get_mut(&id).expect("bound receiving cell");
            let y=power*dt-*correction;
            let next=*sum+y;
            *correction=(next-*sum)-y;
            *sum=next;
            if !sum.is_finite() || !correction.is_finite() {return Err(DepositionError::Arithmetic.into());}
        }
        self.observer.step(start,dt,old,thermal,checkpoint).map_err(NetworkError::Observer)
    }
}
