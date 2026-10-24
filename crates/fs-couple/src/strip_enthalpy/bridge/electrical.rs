//! One enthalpy residual for finite electrical storage and dissipative beams.
//!
//! The original circuit/enthalpy owner supplies the nonlinear solve, Kirchhoff
//! update and branch-local Joule deposition. Its coefficient callback evaluates
//! the original beam from the same accepted history and adds physical Maxwell
//! loss plus signed reversible exchange. No thermal step is nested in another. All three
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
use crate::circuit_enthalpy::contacts::heat::{HeatNetworkError, HeatNetworkSample, PreparedHeatNetwork};
use fs_conduction::transient::deposition::radiation::RadiationCellSources;
use fs_exec::Cx;

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
    Heat(HeatNetworkError),
    Context { cell: CellId, error: SourceError },
    Mechanical(solid::Error),
    Thermodynamics(StripEnthalpyError<Infallible>),
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

pub struct ElectricalReport<S = MaterialNetworkSample> {
    pub electrical: MaterialCircuitStep,
    /// Original explicit A/B contact heat fractions; no default partition.
    /// Contact losses and source snapshots are in `electrical.contacts`.
    pub contact_heat_destinations: Vec<ContactHeatCells>,
    pub mechanical: solid::Report,
    pub layer_energy_j: BTreeMap<(u64,u64), Energy>,
    /// Includes Joule/dashpot heat and reversible exchange; none is added again.
    pub thermal: EnergyAudit,
    /// Original complete transport/electrical material samples at the endpoint.
    pub material: S,
    pub nonlinear: NewtonReport,
    pub circuit_trials: usize,
    pub structural_trials: usize,
    pub maximum_cell_residual_j: f64,
    pub deposited_joule_j: f64,
    /// Signed dashpot heat plus reversible material exchange, already deposited.
    pub deposited_mechanical_heat_j: f64,
    pub reversible_heat_j: f64,
    pub mechanical_internal_energy_change_j: f64,
    pub mechanical_entropy_energy_change_j: f64,
    /// Roundoff-net energy of closed radiative exchange, not external supply.
    /// It is retained in the thermal source account but excluded from external heat.
    pub closed_radiation_net_j: f64,
    /// Independently evaluated boundaries/heaters, excluding internal losses and radiation.
    pub other_external_heat_j: f64,
    /// Delta electrical + mechanical internal + caloric + numerical - supply - load - heat.
    pub total_energy_residual_j: f64,
}

pub struct ElectricalCandidate<'t,'m,S = MaterialNetworkSample> {
    thermal: DepositionSimulation<'t>,
    bridge: solid::Bridge<'m>,
    circuit: MaterialCircuit,
    report: ElectricalReport<S>,
}
impl<'t,'m,S> ElectricalCandidate<'t,'m,S> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> { &self.thermal }
    pub fn bridge(&self) -> &solid::Bridge<'m> { &self.bridge }
    pub fn circuit(&self) -> &MaterialCircuit { &self.circuit }
    pub fn report(&self) -> &ElectricalReport<S> { &self.report }
    pub fn into_parts(self) -> (DepositionSimulation<'t>,solid::Bridge<'m>,MaterialCircuit,ElectricalReport<S>) {
        (self.thermal,self.bridge,self.circuit,self.report)
    }
}
struct Trial<'m,S> {
    bridge: solid::Bridge<'m>,
    mechanical: solid::Report,
    local: BTreeMap<(u64,u64),Energy>,
    material: S,
    heat_j: f64,
    reversible_heat_j: f64,
    mechanical_internal_energy_change_j: f64,
    mechanical_entropy_energy_change_j: f64,
    other_j: f64,
    radiation_net_j: f64,
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

/// Include source-backed thermal contact resistance and closed gray-diffuse
/// radiation in the SAME electrical/beam/enthalpy residual. The prepared network
/// owns material selection, contact series resistance, radiosity and heat routes;
/// no heat step is performed before or after the shared nonlinear solve.
///
/// Retain this prepared network across accepted intervals. Its original endpoint
/// reports survive in `report.material`. Radiation is internal exchange: its
/// roundoff-net power is reported, never treated as external supply or silently
/// redistributed. Fixed reference areas/view factors are not deformed-beam geometry.
#[allow(clippy::too_many_arguments)]
pub fn step_heat_network_candidate<'t,'m>(
    cx: &Cx<'_>, thermal: &DepositionSimulation<'t>, bridge: &solid::Bridge<'m>,
    bindings: &[BridgeThermalBinding<'_>], circuit: &MaterialCircuit,
    transport: &PreparedHeatNetwork<'_,'_>, ports: &[f64],
    target_s: f64, gravity_m_s2: f64, loads_n_m: &[f64],
    config: ElectricalConfig, guess_j: &[f64],
) -> Result<ElectricalCandidate<'t,'m,HeatNetworkSample>,ElectricalError> {
    cx.checkpoint().map_err(|_|ElectricalError::Cancelled)?;
    let destinations=transport.resistor_destinations();
    let contacts=transport.heat_destinations(thermal).map_err(|error|
        ElectricalError::Coupled(CircuitEnthalpyError::Source(EvaluationError::Heat(error))))?;
    let result=driven_candidate(thermal,bridge,bindings,circuit,&destinations,&contacts,ports,
        target_s,gravity_m_s2,loads_n_m,config,guess_j,
        |endpoint|transport.resolve(cx,endpoint).map_err(EvaluationError::Heat),
        ||cx.checkpoint().is_err());
    match result {
        Err(ElectricalError::Coupled(CircuitEnthalpyError::Source(
            EvaluationError::Heat(HeatNetworkError::Cancelled))))=>Err(ElectricalError::Cancelled),
        other=>other,
    }
}

// Keep concrete source reports; only original source owners enter this private
// adapter. No caller can forge a radiative loss or replace its physical solver.
trait SourceSample {
    fn bulk(&self) -> &MaterialNetworkSample;
    fn radiation(&self) -> Option<&RadiationCellSources> { None }
}
impl SourceSample for MaterialNetworkSample {
    fn bulk(&self) -> &MaterialNetworkSample { self }
}
impl SourceSample for HeatNetworkSample {
    fn bulk(&self) -> &MaterialNetworkSample { &self.material }
    fn radiation(&self) -> Option<&RadiationCellSources> { self.radiation.as_ref() }
}

// All entry points use the original contact-capable circuit driver. With no
// contacts that driver calls the original bulk-only circuit update unchanged.
// Share beam admission, physical heat generation and combined energy accounting.
#[allow(clippy::too_many_arguments)]
fn driven_candidate<'t,'m,S: SourceSample>(
    thermal: &DepositionSimulation<'t>, bridge: &solid::Bridge<'m>,
    bindings: &[BridgeThermalBinding<'_>], circuit: &MaterialCircuit,
    destinations: &[ResistorCell], contact_routes: &[ContactHeatCells], ports: &[f64],
    target_s: f64, gravity_m_s2: f64, loads_n_m: &[f64],
    config: ElectricalConfig, guess_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'t>)
        -> Result<ContactThermalCoefficients<S>,EvaluationError>,
    cancelled: impl FnMut()->bool,
) -> Result<ElectricalCandidate<'t,'m,S>,ElectricalError> {
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
        for sample in &material.source.bulk().cells {
            checkpoint()?;
            if let (Some(binding),Some(state))=(by_heat.get(&sample.cell),&sample.conductivity) {
                binding.section.material.admit_state_point(state)
                    .map_err(|error|EvaluationError::Context {cell:sample.cell,error})?;
            }
        }
        // Radiation-only material also needs the mechanical source conditions;
        // absence of a resistor or conductive link cannot bypass that admission.
        if let Some(radiation)=material.source.radiation() {
            for sample in radiation.surfaces() {
                checkpoint()?;
                if let Some(binding)=by_heat.get(&sample.cell) {
                    binding.section.material.admit_state_point(&sample.material)
                        .map_err(|error|EvaluationError::Context {cell:sample.cell,error})?;
                }
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
        let radiation_net_j=finite(dt*material.source.radiation().map_or(0.0,
            RadiationCellSources::net_heat_w))?;
        let other_j=finite(finite(dt*endpoint.external_step_power_w(&material.sources_w,
            &material.conductances_w_k,checkpoint)?)?-radiation_net_j)?;
        let mut temperatures=BTreeMap::new();
        for (&id,binding) in &bound {
            checkpoint()?;
            temperatures.insert(id,endpoint.temperature_k(binding.section.cell)?);
        }
        let (next,mechanical,local)=bridge.advance_with_cell_energy(target_s,&temperatures,
            gravity_m_s2,loads_n_m,mechanics,&mut ||(*cancel.borrow_mut())())
            .map_err(EvaluationError::Mechanical)?;
        let exchange=bridge_exchange(thermal,endpoint,bridge,&next,&bound,&local,checkpoint)
            .map_err(EvaluationError::Thermodynamics)?;
        let mut heat_j=0.0;
        for (id,account) in &exchange.layers {
            checkpoint()?;
            heat_j=finite(heat_j+deposit_heat_j(&mut material.sources_w,bound[id].section.cell,
                account.total_material_heat_j,dt,c.total_energy_tolerance_j)
                .map_err(EvaluationError::Thermodynamics)?)?;
        }
        if (exchange.dissipated_j-mechanical.energy_j.dissipated).abs()>c.total_energy_tolerance_j
            || (heat_j-exchange.dissipated_j-exchange.reversible_heat_j).abs()>c.total_energy_tolerance_j {
            return Err(EvaluationError::Invalid("spatial dashpot/reversible heat differs from total"));
        }
        checkpoint()?;
        Ok(ContactThermalCoefficients {contacts,bulk:CircuitThermalCoefficients {
            resistors:material.resistors,sources_w:material.sources_w,
            conductances_w_k:material.conductances_w_k,source:Trial {bridge:next,mechanical,local,
                material:material.source,heat_j,other_j,radiation_net_j,
                reversible_heat_j:exchange.reversible_heat_j,
                mechanical_internal_energy_change_j:exchange.mechanical_internal_energy_change_j,
                mechanical_entropy_energy_change_j:exchange.mechanical_entropy_energy_change_j}}})
    },||(*cancel.borrow_mut())()).map_err(ElectricalError::Coupled)?;
    let (thermal,next_circuit,r)=candidate.into_parts();
    let contact_heat_destinations=r.source.destinations;
    let trial=r.source.material;
    let total=(r.electrical.stored_after_j-r.electrical.stored_before_j)
        +trial.mechanical_internal_energy_change_j+r.thermal.stored_energy_change_j
        +trial.mechanical.energy_j.numerical_dissipation-r.electrical.native.supplied
        -trial.mechanical.load_work_j-trial.other_j;
    let heat_defect=r.thermal.external_heat_j-r.deposited_joule_j-trial.heat_j
        -trial.other_j-trial.radiation_net_j;
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
            reversible_heat_j:trial.reversible_heat_j,
            mechanical_internal_energy_change_j:trial.mechanical_internal_energy_change_j,
            mechanical_entropy_energy_change_j:trial.mechanical_entropy_energy_change_j,
            closed_radiation_net_j:trial.radiation_net_j,
            total_energy_residual_j:total}})
}
