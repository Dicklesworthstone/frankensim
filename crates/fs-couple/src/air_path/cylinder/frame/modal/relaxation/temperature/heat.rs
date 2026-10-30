//! Finite source-layer heat and rate feedback on the same native trajectory.
//! Exact frozen-start-rate motion supplies original positive Maxwell heat;
//! the existing enthalpy residual advances caloric storage/conduction. Endpoint
//! rate preparation and BOTH histories publish together. This is a first-order
//! thermal split, not an endpoint-power heat estimate or a new integrator.
use super::{FrameRelaxationLayerTemperature, PreparedLayerTemperatures};
use super::shared::{admit_bindings, PreparedSharedTemperatures};
use core::ops::Range;
use super::super::{FrameRelaxationBinding, FrameRelaxationError};
use crate::circuit_enthalpy::{CircuitEnthalpyError, CircuitEnthalpySolve, DrivenCoefficients,
    DriverAccount, ThermalDriver, cell_driver_candidate_entropy_option};
use crate::circuit_enthalpy::material::{MaterialNetwork, MaterialNetworkError,
    MaterialNetworkSample, PreparedMaterialNetwork};
use crate::frame_vibration::thermal::{FrameLayerCell, FrameThermalSolve, FrameHeatEntropy};
use crate::material_circuit::MaterialCircuitSolve;
use crate::render::schedule::force::coupled::{CoupledModalSystem,CoupledModalFrame,ModalCouplingError};
use crate::render::schedule::force::coupled::coefficients::heat::{MaterialHeatBudget,MaterialHeatIntegral};
use crate::render::schedule::force::coupled::coefficients::rates::MaterialRelaxationTimeReport;
use fs_exec::CancelGate;
use fs_solver::NewtonReport;
use fs_thermal::enthalpy::{DepositedCell,DepositionSimulation,EnergyAudit};
use std::collections::{BTreeMap,BTreeSet};

#[derive(Debug)]
pub enum RelaxationHeatSourceError {
    Material(FrameRelaxationError),
    Network(MaterialNetworkError),
    Native(ModalCouplingError),
    Port(fs_ports::CoupleError),
}
impl core::fmt::Display for RelaxationHeatSourceError {
    fn fmt(&self,f:&mut core::fmt::Formatter<'_>)->core::fmt::Result {write!(f,"source relaxation heat: {self:?}")}
}
impl std::error::Error for RelaxationHeatSourceError {}
pub type RelaxationHeatError=CircuitEnthalpyError<RelaxationHeatSourceError>;
type Result<T>=core::result::Result<T,RelaxationHeatError>;
impl From<ModalCouplingError> for RelaxationHeatError {
    fn from(e:ModalCouplingError)->Self {match e {
        ModalCouplingError::Cancelled=>Self::Cancelled,
        other=>Self::Source(RelaxationHeatSourceError::Native(other)),
    }}
}
impl From<FrameRelaxationError> for RelaxationHeatError {
    fn from(e:FrameRelaxationError)->Self {match e {
        FrameRelaxationError::Cancelled=>Self::Cancelled,
        other=>Self::Source(RelaxationHeatSourceError::Material(other)),
    }}
}
impl From<fs_ports::CoupleError> for RelaxationHeatError {
    fn from(e:fs_ports::CoupleError)->Self {Self::Source(RelaxationHeatSourceError::Port(e))}
}
fn network(e:MaterialNetworkError)->RelaxationHeatError {
    match e {
        MaterialNetworkError::Thermal(fs_thermal::enthalpy::DepositionError::Cancelled)=>RelaxationHeatError::Cancelled,
        other=>RelaxationHeatError::Source(RelaxationHeatSourceError::Network(other)),
    }
}
fn invalid(s:&'static str)->RelaxationHeatError {RelaxationHeatError::Invalid(s)}
fn check(gate:&CancelGate)->Result<()> {if gate.is_requested(){Err(RelaxationHeatError::Cancelled)}else{Ok(())}}
fn finite(x:f64)->Result<f64> {if x.is_finite(){Ok(x)}else{Err(invalid("nonfinite relaxation/heat arithmetic"))}}
fn same_positive(a:f64,b:f64)->bool {a.is_finite()&&b.is_finite()&&a>0.0&&b>0.0&&(a-b).abs()<=1e-10*b}
fn controls(c:FrameThermalSolve)->CircuitEnthalpySolve {
    CircuitEnthalpySolve {maximum_cells:c.maximum_cells,maximum_circuit_trials:c.maximum_material_trials,
        maximum_newton_iterations:c.maximum_newton_iterations,difference_step_j:c.difference_step_j,
        cell_energy_tolerance_j:c.cell_energy_tolerance_j,thermal_energy_tolerance_j:c.thermal_energy_tolerance_j,
        total_energy_tolerance_j:c.total_energy_tolerance_j,newton:c.newton,
        circuit:MaterialCircuitSolve {energy_tolerance_j:0.0,power_tolerance_w:0.0}}
}

struct Receiver { destination:FrameLayerCell, original:DepositedCell, source_cell:usize, group:usize }

/// Integrated heat of the ORIGINAL layer, not a fitted sharing fraction.
#[derive(Clone,Copy,Debug,PartialEq)]
pub struct RelaxationLayerHeat {pub destination:FrameLayerCell,pub heat_j:f64}

pub struct RelaxationHeatReport {
    pub mechanical:CoupledModalFrame,
    pub thermal:EnergyAudit,
    pub source:MaterialNetworkSample,
    pub layers:Vec<RelaxationLayerHeat>,
    /// Ranges in `layers`, one per original material binding in declaration
    /// order. Element/layer indices are LOCAL to the corresponding source frame.
    pub material_layer_ranges:Vec<Range<usize>>,
    pub rates:MaterialRelaxationTimeReport,
    pub deposited_material_heat_j:f64,
    /// Other branches, fluid/connection losses and component damping, NOT heat
    /// assigned to this solid. Every unbound Maxwell group is independently integrated.
    pub external_mechanical_loss_j:f64,
    pub external_thermal_heat_j:f64,
    pub total_energy_residual_j:f64,
    pub entropy:FrameHeatEntropy,
    pub nonlinear:NewtonReport,
    pub maximum_cell_residual_j:f64,
    pub material_trials:usize,
}

/// Source-defined frames' finite layers inside ONE complete mechanical assembly.
/// Each binding keeps its own layer inventory; shared conduction uses the same
/// original thermal network. Unbound material memories retain explicit loss
/// sinks. Whole-layer thermal lumping is a caller-selected physical model;
/// internal gradients, expansion, changing storage moduli/mass and phase change
/// are not inferred. The original source rate-only admission remains binding.
pub struct FrameRelaxationThermalSystem<'m,'f,'n,'a,F> {
    materials:Vec<FrameRelaxationBinding<'m,'f,F>>,
    layer_ranges:Vec<Range<usize>>,
    mechanical:CoupledModalSystem,
    source:&'n MaterialNetwork<'a>,
    prepared:PreparedMaterialNetwork<'n,'a>,
    receivers:Vec<Receiver>,
    groups:Vec<core::ops::Range<usize>>,
    outside_groups:Vec<usize>,
    heat_budget:MaterialHeatBudget,
    integral:Option<MaterialHeatIntegral>,
    origin_s:f64,
}

impl<'m,'f,'n,'a,F> FrameRelaxationThermalSystem<'m,'f,'n,'a,F> {
    /// Bind each original material layer to exactly one whole solid caloric cell.
    /// All source cards, contexts, density, physical layer volume and mass must
    /// match. Distinct layers may not silently share a cell. Initial q/v, fluid
    /// momenta and Maxwell histories are retained; no zero-state restart occurs.
    #[allow(clippy::too_many_arguments)]
    pub fn new(material:FrameRelaxationBinding<'m,'f,F>,mechanical:CoupledModalSystem,
        thermal:&DepositionSimulation<'_>,source:&'n MaterialNetwork<'a>,destinations:&[FrameLayerCell],
        heat_budget:MaterialHeatBudget,gate:&CancelGate)->Result<Self> {
        Self::new_shared(vec![(material,destinations.to_vec())],mechanical,thermal,source,heat_budget,gate)
    }

    /// Join independently sourced frames to the SAME finite-heat transaction.
    /// Bindings may arrive in any declaration order and retain local member IDs.
    /// Their native memory blocks must be disjoint, and every layer needs its
    /// own original whole-cell inventory. One thermal network may conduct across
    /// bodies; no separate stepping, temperature averaging or heat fractions.
    ///
    /// All selected material histories, including separated blocks, enter ONE
    /// heat observation and ONE native rate rebuild at the caloric endpoint.
    /// Unselected prefix/interior/suffix histories remain independent sinks.
    /// A common `F` can be a borrowed `&dyn FluidTransport` for distinct original
    /// fluids. No source identities or transport coefficients are erased.
    pub fn new_shared(bindings:Vec<(FrameRelaxationBinding<'m,'f,F>,Vec<FrameLayerCell>)>,
        mechanical:CoupledModalSystem,thermal:&DepositionSimulation<'_>,source:&'n MaterialNetwork<'a>,
        heat_budget:MaterialHeatBudget,gate:&CancelGate)->Result<Self> {
        check(gate)?;
        if mechanical.samples_rendered()!=0 || !source.resistors.is_empty() {
            return Err(invalid("relaxation heat requires a cold shared clock and a thermal-only source network"));
        }
        // Bound all layer inventories before allocating receivers or querying
        // the material network. The complete observer also admits foreign gaps.
        let count=bindings.iter().try_fold(0usize,|sum,(binding,destinations)| {
            let count=binding.members().iter().try_fold(0usize,|n,m|
                n.checked_add(m.source.elastic().layers().len()).ok_or_else(||invalid("shared layer count overflow")))?;
            if count!=destinations.len() {return Err(invalid("one whole-cell destination per original relaxation layer required"));}
            sum.checked_add(count).ok_or_else(||invalid("shared layer count overflow"))
        })?;
        if bindings.is_empty()||count>heat_budget.max_groups {
            return Err(invalid("shared layer heat exceeds its nonempty source-group budget"));
        }
        let (materials,destinations):(Vec<_>,Vec<_>)=bindings.into_iter().unzip();
        admit_bindings(&materials,&mechanical,gate)?;
        let prepared=source.prepare(thermal,||gate.is_requested()).map_err(network)?;
        let mut used=BTreeSet::new();let mut receivers=Vec::with_capacity(count);
        let mut layer_ranges=Vec::with_capacity(materials.len());
        let mut selected=Vec::with_capacity(count);
        for (material,destinations) in materials.iter().zip(destinations) {
            check(gate)?;
            material.states(&mechanical,gate)?;
            let routes:BTreeMap<_,_>=destinations.iter().map(|d|((d.element,d.layer),d.cell)).collect();
            if routes.len()!=destinations.len() {return Err(invalid("duplicate layer heat destination"));}
            let layer_start=receivers.len();
            let range=material.memory_range();let mut cursor=range.start;
            for member in material.members() {
                let (_,_,length)=material.source().frame().reference_strain_traces(member.element)
                    .map_err(|e|FrameRelaxationError::Frame(e.into()))?;
                for (layer,(solid,geometry)) in member.source.elastic().layers().iter()
                    .zip(member.source.elastic().layer_integrals()).enumerate() {
                    check(gate)?;
                    let cell=*routes.get(&(member.element,layer)).ok_or_else(||invalid("missing original layer cell"))?;
                    if !used.insert(cell){return Err(invalid("each relaxation layer across ALL frames needs its own whole caloric inventory"));}
                    let index=source.cells.iter().position(|s|s.cell==cell).ok_or_else(||invalid("missing original caloric source"))?;
                    let caloric=&source.cells[index];
                    let original=thermal.cells().get(&cell).ok_or_else(||invalid("inactive layer heat cell"))?;
                    let density=geometry.line_density_kg_m()/geometry.area_m2();
                    let initial_temperature=solid.material.query_point().iter()
                        .find(|(axis,_)|axis==caloric.point.temperature_axis).map(|(_,t)|*t)
                        .ok_or_else(||invalid("source layer lacks its caloric temperature axis"))?;
                    if solid.material.card_identity()!=caloric.point.material.content_hash()
                        || solid.material.query()!=*caloric.point.context
                        || original.material_id!=caloric.material_id
                        || !same_positive(original.mass_kg,geometry.line_density_kg_m()*length)
                        || !same_positive(original.geometry.volume_m3,geometry.area_m2()*length)
                        || caloric.curve.knots().iter().any(|k|k.liquid_mass_fraction!=0.0||!same_positive(k.bulk_density_kg_m3,density))
                        || thermal.temperature_k(cell)?.to_bits()!=initial_temperature.to_bits() {
                        return Err(invalid("relaxation heat must retain original card/context and whole fixed-solid mass/volume/caloric law"));
                    }
                    let first=cursor;
                    while material.coordinates().get(cursor-range.start).is_some_and(|c|c.element==member.element&&c.layer==layer) {cursor+=1;}
                    selected.push((first..cursor,receivers.len()));
                    receivers.push(Receiver {destination:FrameLayerCell {element:member.element,layer,cell},
                        original:original.clone(),source_cell:index,group:0});
                }
            }
            if cursor!=range.end {return Err(invalid("layer heat did not cover every original material memory"));}
            layer_ranges.push(layer_start..receivers.len());
        }
        // Physical reporting stays in declaration order; heat ranges are in
        // native memory order. Include every FOREIGN INTERIOR GAP, not only a
        // prefix and suffix. No per-frame allocation of a global loss total.
        selected.sort_unstable_by_key(|(r,i)|(r.start,r.end,*i));
        let mut groups=Vec::new();let mut outside_groups=Vec::new();let mut cursor=0;
        for (range,receiver) in selected {
            check(gate)?;
            if range.start<cursor {return Err(invalid("selected layer heat ranges overlap"));}
            if cursor<range.start {outside_groups.push(groups.len());groups.push(cursor..range.start);}
            receivers[receiver].group=groups.len();cursor=range.end;groups.push(range);
        }
        if cursor<mechanical.material_memory().len() {
            outside_groups.push(groups.len());groups.push(cursor..mechanical.material_memory().len());
        }
        if groups.len()>heat_budget.max_groups {return Err(invalid("complete selected and foreign heat groups exceed the budget"));}
        let mut owner=Self {materials,layer_ranges,mechanical,source,prepared,receivers,groups,outside_groups,heat_budget,
            integral:None,origin_s:thermal.time_s()};
        // Same pinned rate-only source admission, including validity-wide
        // CONSTANT storage. All bodies are queried before any publication.
        propose_shared(&owner.materials,&owner.mechanical,source,&owner.receivers,
            &owner.layer_ranges,thermal,gate)?;
        owner.integral=Some(owner.mechanical.prepare_material_heat(&owner.groups,heat_budget,gate)?);
        check(gate)?;Ok(owner)
    }

    pub fn mechanical(&self)->&CoupledModalSystem {&self.mechanical}
    /// First declared binding; existing single-frame callers retain their meaning.
    pub fn material(&self)->&FrameRelaxationBinding<'m,'f,F> {&self.materials[0]}
    pub fn materials(&self)->&[FrameRelaxationBinding<'m,'f,F>] {&self.materials}
    /// One original binding's receiver/report range, in declaration order.
    pub fn material_layer_ranges(&self)->&[Range<usize>] {&self.layer_ranges}

    /// One exact frozen-rate mechanical sample, one original enthalpy solve,
    /// then original rates at the receiving endpoint temperatures. Physical heat
    /// is integrated per layer; it is never an energy residual or endpoint power.
    /// No partial mechanics, caloric state, source snapshot or clock publishes.
    /// Changing rates has no storage jump. Refine this first-order split in time.
    pub fn step(&mut self,thermal:&mut DepositionSimulation<'_>,external:&[f64],
        solve:FrameThermalSolve,energy_guess_j:&[f64],gate:&CancelGate)->Result<RelaxationHeatReport> {
        check(gate)?;
        let sample=self.mechanical.samples_rendered();
        let next=sample.checked_add(1).ok_or_else(||invalid("relaxation heat clock overflow"))?;
        let dt=self.mechanical.sample_period_s();
        let start=finite(self.origin_s+sample as f64*dt)?;
        let target=finite(self.origin_s+next as f64*dt)?;
        if thermal.time_s().to_bits()!=start.to_bits()||target<=start {
            return Err(invalid("relaxation and caloric histories must share their original advancing clock"));
        }
        for (receiver,layer) in self.receivers.iter().zip(self.materials.iter().flat_map(|b|b.members()).flat_map(|m|m.source.elastic().layers())) {
            let axis=self.source.cells[receiver.source_cell].point.temperature_axis;
            let previous=layer.material.query_point().iter().find(|(name,_)|name==axis).map(|(_,t)|*t)
                .ok_or_else(||invalid("accepted material temperature axis disappeared"))?;
            if thermal.temperature_k(receiver.destination.cell)?.to_bits()!=previous.to_bits() {
                return Err(invalid("finite heat and relaxation rates were changed outside joint publication"));
            }
        }
        if self.integral.as_ref().map_or(Ok(false),|i|i.is_current(&self.mechanical,gate))?==false {
            self.integral=Some(self.mechanical.prepare_material_heat(&self.groups,self.heat_budget,gate)?);
        }
        let heat=self.integral.as_ref().ok_or_else(||invalid("material heat integral is absent"))?
            .evaluate(&self.mechanical,external,gate)?;
        let total=heat.iter().try_fold(0.0,|s,h|finite(s+h))?;
        let selected=self.receivers.iter().try_fold(0.0,|s,r|finite(s+heat[r.group]))?;
        let foreign=self.outside_groups.iter().try_fold(0.0,|s,&i|finite(s+heat[i]))?;
        let before=self.mechanical.total_energy_j()?;
        let (materials,source,prepared,receivers,layer_ranges)=
            (&self.materials,self.source,&self.prepared,&self.receivers,&self.layer_ranges);
        let (candidate,mut report,source_update,native_update)=self.mechanical.step_with_admission(external,gate,
            |mechanical,frame,_transfers| -> Result<_> {
                if finite(total-frame.material_dissipation_j)?.abs()>solve.total_energy_tolerance_j {
                    return Err(invalid("independent source-layer heat does not equal native total material loss"));
                }
                let outside=finite(foreign+frame.component_dissipation_j+frame.connection_dissipation_j+frame.radiation_dissipation_j)?;
                let layers:Vec<_>=receivers.iter().map(|r|RelaxationLayerHeat {destination:r.destination,heat_j:heat[r.group]}).collect();
                let (candidate,(),audit,entropy)=cell_driver_candidate_entropy_option(thermal,target,controls(solve),energy_guess_j,
                    |trial,_remaining| -> Result<_> {
                        check(gate)?;
                        let coefficients=prepared.resolve(trial,||gate.is_requested()).map_err(network)?;
                        let update=propose_shared(materials,mechanical,source,receivers,layer_ranges,trial,gate)?;
                        Ok(DrivenCoefficients {circuit:Driver {frame:frame.clone(),before,heat:selected,outside},
                            sources_w:coefficients.sources_w,conductances_w_k:coefficients.conductances_w_k,
                            source:(coefficients.source,update),circuit_trials:1})
                    },
                    |_driver,duration,checkpoint| {
                        let mut powers=BTreeMap::new();
                        for layer in &layers {
                            checkpoint()?;
                            let power=finite(layer.heat_j/duration)?;
                            if layer.heat_j>0.0&&power==0.0 {return Err(invalid("layer heat power underflow"));}
                            powers.insert(layer.destination.cell,power);
                        }
                        Ok(powers)
                    },true,||gate.is_requested())?;
                let (source,update)=audit.source;
                let native=mechanical.prepare_staged_material_relaxation_times(&update.updates,gate)?;
                let entropy=entropy.ok_or_else(||invalid("missing original material heat entropy account"))?;
                let report=RelaxationHeatReport {mechanical:audit.electrical,thermal:audit.thermal,source,layers,
                    material_layer_ranges:layer_ranges.to_vec(),rates:MaterialRelaxationTimeReport::default(),deposited_material_heat_j:audit.deposited_joule_j,
                    external_mechanical_loss_j:outside,external_thermal_heat_j:audit.non_electrical_heat_j,
                    total_energy_residual_j:audit.total_energy_residual_j,entropy:entropy.into(),nonlinear:audit.nonlinear,
                    maximum_cell_residual_j:audit.maximum_cell_residual_j,material_trials:audit.circuit_trials};
                Ok((candidate,report,update,native))
            })?;
        // The native admission made the final cancellation check. No fallible
        // work may split publication of mechanics, source conditions and heat.
        report.rates=self.mechanical.install_material_relaxation_times(native_update);
        source_update.install(&mut self.materials);
        *thermal=candidate;
        if !report.rates.changes.is_empty(){self.integral=None;}
        Ok(report)
    }
}

fn propose_shared<F>(materials:&[FrameRelaxationBinding<'_, '_, F>],mechanical:&CoupledModalSystem,
    source:&MaterialNetwork<'_>,receivers:&[Receiver],ranges:&[Range<usize>],
    thermal:&DepositionSimulation<'_>,gate:&CancelGate)->Result<PreparedSharedTemperatures> {
    let mut proposals=Vec::with_capacity(materials.len());
    for (material,range) in materials.iter().zip(ranges) {
        check(gate)?;
        proposals.push(propose(material,mechanical,source,&receivers[range.clone()],thermal,gate)?);
    }
    PreparedSharedTemperatures::new(proposals,mechanical.material_memory().len(),gate).map_err(Into::into)
}

fn propose<F>(material:&FrameRelaxationBinding<'_, '_, F>,mechanical:&CoupledModalSystem,
    source:&MaterialNetwork<'_>,receivers:&[Receiver],thermal:&DepositionSimulation<'_>,gate:&CancelGate)
    ->Result<PreparedLayerTemperatures> {
    let mut inputs=Vec::with_capacity(receivers.len());
    for receiver in receivers {
        check(gate)?;
        let cell=thermal.cells().get(&receiver.destination.cell).ok_or_else(||invalid("layer heat cell disappeared"))?;
        let original=&receiver.original;
        if cell.id!=original.id||cell.material_id!=original.material_id||cell.body_id!=original.body_id
            ||cell.segment_id!=original.segment_id||cell.mass_kg.to_bits()!=original.mass_kg.to_bits()||cell.geometry!=original.geometry {
            return Err(invalid("original layer caloric inventory or geometry changed"));
        }
        let point=&source.cells[receiver.source_cell].point;
        inputs.push(FrameRelaxationLayerTemperature {element:receiver.destination.element,layer:receiver.destination.layer,
            card:point.material,temperature_axis:point.temperature_axis,temperature_k:thermal.temperature_k(cell.id)?});
    }
    material.prepare_layer_temperatures(mechanical,&inputs,gate).map_err(Into::into)
}
struct Driver {frame:CoupledModalFrame,before:f64,heat:f64,outside:f64}
impl ThermalDriver for Driver {
    type Model=();type Step=CoupledModalFrame;
    fn account(&self)->DriverAccount<'_> {DriverAccount {resistors:&[],contacts:&[],heat_j:self.heat,
        numerical_dissipation_j:self.frame.numerical_dissipation_j,
        stored_change_j:self.frame.total_energy_j()-self.before,supplied_j:self.frame.external_work_j-self.outside}}
    fn into_parts(self)->((),CoupledModalFrame){((),self.frame)}
}

#[cfg(test)]
mod tests;
