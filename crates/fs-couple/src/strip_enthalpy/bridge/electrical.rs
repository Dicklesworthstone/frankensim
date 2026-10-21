//! One enthalpy residual for finite electrical storage and dissipative beams.
//!
//! The original circuit/enthalpy owner supplies the nonlinear solve, Kirchhoff
//! update and branch-local Joule deposition. Its coefficient callback evaluates
//! the original beam from the same accepted history and adds ONLY physical
//! Maxwell heat. No thermal step is nested in another thermal step. All three
//! endpoints publish together after independent combined work/storage admission.
//! Fixed C/L, electrical paths and thermal geometry; no field force or reciprocal
//! thermoelastic expansion is inferred. Endpoint R(T), k(T), relaxation and BE
//! mechanics make the coupled method first-order in general.

use core::convert::Infallible;
use super::*;
use crate::circuit_enthalpy::{self, CircuitEnthalpyError, CircuitEnthalpyInput,
    CircuitEnthalpySolve, CircuitThermalCoefficients, ResistorCell};
use crate::circuit_enthalpy::contacts::{ContactHeatCells, ContactThermalCoefficients};
use crate::circuit_enthalpy::contacts::material::{
    ElectricalContactNetworkError, PreparedElectricalContactNetwork,
};
use crate::material_circuit::{MaterialCircuit, MaterialCircuitSolve, MaterialCircuitStep};
use crate::circuit_enthalpy::material::MaterialNetworkError;
use fs_material::thermovisco::source::SourceError;
use fs_thermal::enthalpy::CellId;

#[derive(Clone, Copy, Debug)]
pub struct ElectricalConfig {
    pub bridge: BridgeEnthalpyConfig,
    pub circuit: MaterialCircuitSolve,
    /// Includes every electrical solve and final fresh endpoint admission.
    pub maximum_circuit_trials: usize,
}

/// Original material and mechanics refusals survive the shared numerical owner.
#[derive(Debug)]
pub enum EvaluationError {
    Material(MaterialNetworkError),
    Contact(ElectricalContactNetworkError),
    Context { cell: CellId, error: SourceError },
    Mechanical(solid::Error),
    Thermal(DepositionError),
    Invalid(&'static str),
    StructuralBudget,
}
impl From<DepositionError> for EvaluationError {
    fn from(e: DepositionError) -> Self { Self::Thermal(e) }
}
impl core::fmt::Display for EvaluationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "electrical/beam source: {self:?}")
    }
}
impl std::error::Error for EvaluationError {}

#[derive(Debug)]
pub enum ElectricalError {
    Binding(BridgeEnthalpyError<Infallible>),
    Coupled(CircuitEnthalpyError<EvaluationError>),
    Invalid(&'static str),
    Cancelled,
}
impl core::fmt::Display for ElectricalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "electrical/beam/enthalpy: {self:?}")
    }
}
impl std::error::Error for ElectricalError {}

pub struct ElectricalReport {
    pub electrical: MaterialCircuitStep,
    /// Original explicit A/B contact heat fractions; no default partition.
    /// Contact losses and source snapshots are in `electrical.contacts`.
    pub contact_heat_destinations: Vec<ContactHeatCells>,
    pub mechanical: solid::Report,
    pub layer_energy_j: BTreeMap<(u64,u64), Energy>,
    /// Includes both physical Joule and Maxwell heat; neither is added again.
    pub thermal: EnergyAudit,
    /// Original complete transport/electrical material samples at the endpoint.
    pub material: MaterialNetworkSample,
    pub nonlinear: NewtonReport,
    pub circuit_trials: usize,
    pub structural_trials: usize,
    pub maximum_cell_residual_j: f64,
    pub deposited_joule_j: f64,
    pub deposited_mechanical_heat_j: f64,
    /// Independently evaluated boundaries/heaters, excluding both internal losses.
    pub other_external_heat_j: f64,
    /// Delta electrical + elastic + thermal + numerical - supply - load - heat.
    pub total_energy_residual_j: f64,
}

pub struct ElectricalCandidate<'t,'m> {
    thermal: DepositionSimulation<'t>,
    bridge: solid::Bridge<'m>,
    circuit: MaterialCircuit,
    report: ElectricalReport,
}
impl<'t,'m> ElectricalCandidate<'t,'m> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> { &self.thermal }
    pub fn bridge(&self) -> &solid::Bridge<'m> { &self.bridge }
    pub fn circuit(&self) -> &MaterialCircuit { &self.circuit }
    pub fn report(&self) -> &ElectricalReport { &self.report }
    pub fn into_parts(self) -> (DepositionSimulation<'t>,solid::Bridge<'m>,MaterialCircuit,ElectricalReport) {
        (self.thermal,self.bridge,self.circuit,self.report)
    }
}
struct Trial<'m> {
    bridge: solid::Bridge<'m>,
    mechanical: solid::Report,
    local: BTreeMap<(u64,u64),Energy>,
    material: MaterialNetworkSample,
    heat_j: f64,
    other_j: f64,
}

/// Simultaneously solve material R(T), k(T), relaxation and all thermal cells.
///
/// Heat routes come from the prepared network itself, not a second mapping.
/// Every circuit and beam proposal starts from the original charges/fluxes and
/// Maxwell histories. Only the final accepted trial's three states are returned.
/// Authored thermal boundaries remain active. Target events remain pending;
/// process them explicitly after acceptance. Electrical topology changes and
/// series electrical contacts need their own state/heat-route owners and are
/// not silently admitted here. The original circuit driver enforces that limit.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate<'t,'m>(
    thermal: &DepositionSimulation<'t>, bridge: &solid::Bridge<'m>,
    bindings: &[BridgeThermalBinding<'_>], circuit: &MaterialCircuit,
    transport: &PreparedMaterialNetwork<'_,'_>, ports: &[f64],
    target_s: f64, gravity_m_s2: f64, loads_n_m: &[f64],
    config: ElectricalConfig, guess_j: &[f64], cancelled: impl FnMut()->bool,
) -> Result<ElectricalCandidate<'t,'m>,ElectricalError> {
    let cancel=RefCell::new(cancelled);
    let destinations=transport.resistor_destinations();
    driven_candidate(thermal,bridge,bindings,circuit,&destinations,&[],ports,
        target_s,gravity_m_s2,loads_n_m,config,guess_j,|endpoint| {
            let bulk=transport.resolve(endpoint,||(*cancel.borrow_mut())())
                .map_err(EvaluationError::Material)?;
            Ok(ContactThermalCoefficients {bulk,contacts:Vec::new()})
        },||(*cancel.borrow_mut())())
}

/// The same joint electrical/beam/thermal solve, including interface Joule heat.
///
/// Contact resistance is evaluated at both actual endpoint material states by
/// the existing ordered-interface owner. Bulk and two-sided contact heat routes
/// come from this SAME prepared network. Physical contact and Maxwell losses
/// are each deposited once by their original owners. Retain `transport` across
/// accepted steps to preserve bulk and interface source selections.
///
/// The thermal link provides fixed area and A/B order; its explicitly supplied
/// heat fraction is not a predicted effusivity partition. No adhesion, pressure
/// evolution, opening, arc, field force or contact storage is inferred here.
#[allow(clippy::too_many_arguments)]
pub fn step_with_contacts_candidate<'t,'m>(
    thermal: &DepositionSimulation<'t>, bridge: &solid::Bridge<'m>,
    bindings: &[BridgeThermalBinding<'_>], circuit: &MaterialCircuit,
    transport: &PreparedElectricalContactNetwork<'_,'_>, ports: &[f64],
    target_s: f64, gravity_m_s2: f64, loads_n_m: &[f64],
    config: ElectricalConfig, guess_j: &[f64], cancelled: impl FnMut()->bool,
) -> Result<ElectricalCandidate<'t,'m>,ElectricalError> {
    let cancel=RefCell::new(cancelled);
    if (*cancel.borrow_mut())() { return Err(ElectricalError::Cancelled); }
    let destinations=transport.resistor_destinations();
    let contacts=transport.heat_destinations(thermal).map_err(|error|
        ElectricalError::Coupled(CircuitEnthalpyError::Source(EvaluationError::Contact(error))))?;
    driven_candidate(thermal,bridge,bindings,circuit,&destinations,&contacts,ports,
        target_s,gravity_m_s2,loads_n_m,config,guess_j,
        |endpoint|transport.resolve(endpoint,||(*cancel.borrow_mut())()).map_err(EvaluationError::Contact),
        ||(*cancel.borrow_mut())())
}

// Both entry points use the original contact-capable circuit driver. With no
// contacts that driver calls the original bulk-only circuit update unchanged.
// Share beam admission, physical heat generation and combined energy accounting.
#[allow(clippy::too_many_arguments)]
fn driven_candidate<'t,'m>(
    thermal: &DepositionSimulation<'t>, bridge: &solid::Bridge<'m>,
    bindings: &[BridgeThermalBinding<'_>], circuit: &MaterialCircuit,
    destinations: &[ResistorCell], contact_routes: &[ContactHeatCells], ports: &[f64],
    target_s: f64, gravity_m_s2: f64, loads_n_m: &[f64],
    config: ElectricalConfig, guess_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'t>)
        -> Result<ContactThermalCoefficients<MaterialNetworkSample>,EvaluationError>,
    cancelled: impl FnMut()->bool,
) -> Result<ElectricalCandidate<'t,'m>,ElectricalError> {
    let cancel=RefCell::new(cancelled);
    let checkpoint=|| if (*cancel.borrow_mut())() { Err(DepositionError::Cancelled) } else { Ok(()) };
    checkpoint().map_err(|_|ElectricalError::Cancelled)?;
    let c=config.bridge.feedback;
    if c.maximum_structural_trials==0 || bridge.cells().len()>c.maximum_cells {
        return Err(ElectricalError::Invalid("bounded structural work and cells required"));
    }
    let bound=super::admit_bindings::<Infallible>(thermal,bridge,bindings,checkpoint)
        .map_err(ElectricalError::Binding)?;
    let by_heat:BTreeMap<_,_>=bindings.iter().map(|b|(b.section.cell,b)).collect();
    let routes:BTreeMap<_,_>=destinations.iter().map(|r|(r.branch,r.cell)).collect();
    let dt=target_s-thermal.time_s();
    let mechanics=solid::Config {section:c.mechanics,maximum_deflection_m:config.bridge.maximum_deflection_m};
    let controls=CircuitEnthalpySolve {maximum_cells:c.maximum_cells,
        maximum_circuit_trials:config.maximum_circuit_trials,maximum_newton_iterations:c.maximum_newton_iterations,
        difference_step_j:c.difference_step_j,cell_energy_tolerance_j:c.cell_energy_tolerance_j,
        thermal_energy_tolerance_j:c.thermal_energy_tolerance_j,total_energy_tolerance_j:c.total_energy_tolerance_j,
        circuit:config.circuit,newton:c.newton};
    let trials=Cell::new(0usize);
    let candidate=circuit_enthalpy::contacts::step_candidate(thermal,target_s,
        CircuitEnthalpyInput {circuit,ports,destinations},contact_routes,controls,guess_j,|endpoint| {
        checkpoint()?;
        if trials.get()>=c.maximum_structural_trials { return Err(EvaluationError::StructuralBudget); }
        trials.set(trials.get()+1);
        let ContactThermalCoefficients {bulk:mut material,contacts}=resolve(endpoint)?;
        // Conductive and electrical-only cells both need the original mechanical
        // source context; absence of k does not bypass electrical admission.
        for sample in &material.source.cells {
            checkpoint()?;
            if let (Some(binding),Some(state))=(by_heat.get(&sample.cell),&sample.conductivity) {
                binding.section.material.admit_state_point(state)
                    .map_err(|error|EvaluationError::Context {cell:sample.cell,error})?;
            }
        }
        for resistor in &material.resistors {
            checkpoint()?;
            let cell=*routes.get(&resistor.branch)
                .ok_or(EvaluationError::Invalid("material resistor lacks its original heat route"))?;
            if let Some(binding)=by_heat.get(&cell) {
                binding.section.material.admit_state_point(resistor.conductor.material())
                    .map_err(|error|EvaluationError::Context {cell,error})?;
            }
        }
        let other_j=finite(dt*endpoint.external_step_power_w(&material.sources_w,
            &material.conductances_w_k,checkpoint)?)?;
        let mut temperatures=BTreeMap::new();
        for (&id,binding) in &bound {
            checkpoint()?;
            temperatures.insert(id,endpoint.temperature_k(binding.section.cell)?);
        }
        let (next,mechanical,local)=bridge.advance_with_cell_energy(target_s,&temperatures,
            gravity_m_s2,loads_n_m,mechanics,&mut ||(*cancel.borrow_mut())())
            .map_err(EvaluationError::Mechanical)?;
        let mut heat_j=0.0;
        for (id,energy) in &local {
            checkpoint()?;
            if !energy.dissipated.is_finite() || energy.dissipated<0.0 {
                return Err(EvaluationError::Invalid("nonphysical mechanical heat"));
            }
            let power=finite(energy.dissipated/dt)?;
            if energy.dissipated!=0.0 && power==0.0 {
                return Err(EvaluationError::Invalid("mechanical heat underflowed to zero power"));
            }
            heat_j=finite(heat_j+finite(dt*power)?)?;
            if power!=0.0 {
                let source=material.sources_w.entry(bound[id].section.cell).or_insert(0.0);
                *source=finite(*source+power)?;
            }
        }
        if (heat_j-mechanical.energy_j.dissipated).abs()>c.total_energy_tolerance_j {
            return Err(EvaluationError::Invalid("local and total physical mechanical heat differ"));
        }
        checkpoint()?;
        Ok(ContactThermalCoefficients {contacts,bulk:CircuitThermalCoefficients {
            resistors:material.resistors,sources_w:material.sources_w,
            conductances_w_k:material.conductances_w_k,source:Trial {bridge:next,mechanical,local,
                material:material.source,heat_j,other_j}}})
    },||(*cancel.borrow_mut())()).map_err(ElectricalError::Coupled)?;
    let (thermal,next_circuit,r)=candidate.into_parts();
    let contact_heat_destinations=r.source.destinations;
    let trial=r.source.material;
    let total=(r.electrical.stored_after_j-r.electrical.stored_before_j)
        +(trial.bridge.stored_energy_j()-bridge.stored_energy_j())+r.thermal.stored_energy_change_j
        +trial.mechanical.energy_j.numerical_dissipation-r.electrical.native.supplied
        -trial.mechanical.load_work_j-trial.other_j;
    let heat_defect=r.thermal.external_heat_j-r.deposited_joule_j-trial.heat_j-trial.other_j;
    if !total.is_finite() || !heat_defect.is_finite() || total.abs()>c.total_energy_tolerance_j
        || heat_defect.abs()>c.total_energy_tolerance_j {
        return Err(ElectricalError::Invalid("electrical, mechanical and thermal work/storage do not close"));
    }
    checkpoint().map_err(|_|ElectricalError::Cancelled)?;
    Ok(ElectricalCandidate {thermal,bridge:trial.bridge,circuit:next_circuit,
        report:ElectricalReport {electrical:r.electrical,contact_heat_destinations,mechanical:trial.mechanical,
            layer_energy_j:trial.local,thermal:r.thermal,material:trial.material,nonlinear:r.nonlinear,
            circuit_trials:r.circuit_trials,structural_trials:trials.get(),
            maximum_cell_residual_j:r.maximum_cell_residual_j,deposited_joule_j:r.deposited_joule_j,
            deposited_mechanical_heat_j:trial.heat_j,other_external_heat_j:trial.other_j,
            total_energy_residual_j:total}})
}
