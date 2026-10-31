//! Source-bound breakthrough of one pure permeant in a dense laminate.
//!
//! The original material cards explicitly supply sorbed-phase diffusivity D
//! [m2/s] and concentration-based Henry solubility S [mol/(m3 Pa)]. P = D S
//! shares the steady membrane coefficient owner. No generic material defaults,
//! inverse Henry convention, intrinsic porous k or STP-volume conversion enter.
//!
//! This bounded rung holds temperature, geometry and all other authored source
//! coordinates fixed. Pure-gas reservoir pressures are prescribed; their molar
//! transfers and the membrane's finite inventory are separate in every report.
//! The thermostat and reservoirs are external, with no invented closed-system
//! energy balance. Mixed-gas competition, swelling and nonlinear sorption need
//! a different constitutive model, not silent extrapolation of these sources.

use std::sync::Arc;
use super::{MembranePermeabilityMaterial, MembraneInterfaceMaterial, MembraneDiffusivity,
    MolarHenrySolubility, MolarPermeance, MEMBRANE_DIFFUSIVITY_DIMS,
    MOLAR_HENRY_SOLUBILITY_DIMS, MOLAR_PERMEANCE_DIMS, MOLAR_INTERFACE_PERMEANCE_PROPERTY,
    GasEnergyError, updated_point, source_error};
use super::series::admit_henry;
use crate::state_point::{ResolvedMaterialStatePoint, ResolvedInterfaceStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
    resolve_interface_state_point};
use fs_matdb::QueryPoint;
use fs_thermochem::SpeciesId;

pub mod slab;
pub use slab::{MembraneDiffusion, MembraneDiffusionLayer, MembraneDiffusionState,
    MembraneDiffusionWorkspace, MembraneDiffusionSolve, MembraneDiffusionFlux,
    MembraneDiffusionStep, MembraneDiffusionCandidate, MembraneDiffusionError,
    MAX_TRANSIENT_MEMBRANE_LAYERS, MAX_TRANSIENT_MEMBRANE_CELLS};

pub const MEMBRANE_DIFFUSIVITY_PROPERTY:&str="gas_membrane_diffusivity";
pub const MOLAR_HENRY_SOLUBILITY_PROPERTY:&str="gas_molar_henry_solubility";

/// Geometry and the existing species/material binding; this consumer resolves
/// BOTH D and S from that original card and selection, not a second schema.
#[derive(Clone,Debug)]
pub struct MembraneDiffusionLayerSource<'a> {
    pub thickness_m:f64,
    pub cells:usize,
    pub material:MembranePermeabilityMaterial<'a>,
}

#[derive(Clone,Debug,PartialEq)]
pub enum ResolvedMembraneDiffusionError {
    Invalid(&'static str),
    Source(GasEnergyError),
    Diffusion(MembraneDiffusionError),
    Cancelled,
}
impl core::fmt::Display for ResolvedMembraneDiffusionError {
    fn fmt(&self,f:&mut core::fmt::Formatter<'_>)->core::fmt::Result {
        write!(f,"source-bound transient membrane: {self:?}")
    }
}
impl std::error::Error for ResolvedMembraneDiffusionError {}
impl From<GasEnergyError> for ResolvedMembraneDiffusionError {
    fn from(error:GasEnergyError)->Self {
        if error==GasEnergyError::Cancelled {Self::Cancelled}else{Self::Source(error)}
    }
}
impl From<MembraneDiffusionError> for ResolvedMembraneDiffusionError {
    fn from(error:MembraneDiffusionError)->Self {
        if error==MembraneDiffusionError::Cancelled {Self::Cancelled}else{Self::Diffusion(error)}
    }
}
use ResolvedMembraneDiffusionError as Error;
fn poll(cancelled:&mut impl FnMut()->bool)->Result<(),Error> {
    if cancelled() {Err(Error::Cancelled)}else{Ok(())}
}

#[derive(Clone,Debug)]
struct LayerBinding<'a> {
    binding:MembranePermeabilityMaterial<'a>,
    original:ResolvedMaterialStatePoint,
}
#[derive(Clone,Debug)]
struct InterfaceBinding<'a> {
    binding:MembraneInterfaceMaterial<'a>,
    original:ResolvedInterfaceStatePoint,
}
#[derive(Clone,Debug)]
struct Bindings<'a> {
    species:SpeciesId,
    temperature_k:f64,
    layers:Vec<LayerBinding<'a>>,
    interfaces:Vec<Option<InterfaceBinding<'a>>>,
}

/// Exact original-claim receipts at minimum and maximum actual Henry pressure.
/// For a step, these include both its initial and final cells and local faces.
/// Constant-at-T laws and interval applicability cover the intervening pressure
/// interval. Each interface range includes both of its actual ordered sides.
#[derive(Clone,Debug,PartialEq)]
pub struct MembraneDiffusionSourceSupport {
    pub layers:Vec<[ResolvedMaterialStatePoint;2]>,
    pub interfaces:Vec<Option<[ResolvedInterfaceStatePoint;2]>>,
}

#[derive(Clone,Debug)]
pub struct ResolvedMembraneDiffusion<'a> {
    bindings:Arc<Bindings<'a>>,
    state:MembraneDiffusionState,
    support:MembraneDiffusionSourceSupport,
}
#[derive(Clone,Debug)]
pub struct ResolvedMembraneDiffusionCandidate<'a> {
    state:ResolvedMembraneDiffusion<'a>,
    step:MembraneDiffusionStep,
}
impl<'a> ResolvedMembraneDiffusionCandidate<'a> {
    pub const fn state(&self)->&ResolvedMembraneDiffusion<'a> {&self.state}
    pub const fn step(&self)->&MembraneDiffusionStep {&self.step}
    pub fn into_parts(self)->(ResolvedMembraneDiffusion<'a>,MembraneDiffusionStep) {(self.state,self.step)}
}

fn selection_point(context:&QueryPoint,temperature_k:f64)->Result<QueryPoint,Error> {
    let pressure=context.axes().get("pressure").copied()
        .filter(|p|p.is_finite()&&*p>=0.0)
        .ok_or(Error::Invalid("each membrane source requires an explicit nonnegative canonical selection pressure"))?;
    Ok(updated_point(context,temperature_k,pressure)?)
}

impl<'a> ResolvedMembraneDiffusion<'a> {
    /// One pure permeant, isothermal planar layers and a declared initial
    /// loading per cell. None is equilibrium interlayer contact; Some uses the
    /// existing ordered-interface binding and molar pressure-based permeance.
    /// Initial reservoirs are supplied so original source support includes the
    /// actual initial faces, including a degassed film's inlet boundary jump.
    pub fn new(temperature_k:f64,area_m2:f64,layers:&[MembraneDiffusionLayerSource<'a>],
        interfaces:&[Option<MembraneInterfaceMaterial<'a>>],initial_partial_pressures_pa:&[f64],
        initial_reservoirs_pa:[f64;2],mut cancelled:impl FnMut()->bool)->Result<Self,Error>
    {
        poll(&mut cancelled)?;
        if !temperature_k.is_finite()||temperature_k<=0.0||!area_m2.is_finite()||area_m2<=0.0
            ||layers.is_empty()||layers.len()>MAX_TRANSIENT_MEMBRANE_LAYERS
            ||interfaces.len()!=layers.len()-1 {
            return Err(Error::Invalid("positive temperature/area and a bounded nonempty ordered laminate required"));
        }
        let count=layers.iter().try_fold(0usize,|count,layer|count.checked_add(layer.cells));
        if count!=Some(initial_partial_pressures_pa.len())||initial_partial_pressures_pa.len()>MAX_TRANSIENT_MEMBRANE_CELLS {
            return Err(Error::Invalid("bounded initial sorbed loading must cover every authored cell"));
        }
        let species=layers[0].material.species.clone();
        let mut bound_layers=Vec::with_capacity(layers.len());
        let mut coefficients=Vec::with_capacity(layers.len());
        let d_requirement=ScalarPropertyRequirement::try_new(MEMBRANE_DIFFUSIVITY_PROPERTY,
            MEMBRANE_DIFFUSIVITY_DIMS,ScalarAdmissibility::NonNegative).map_err(source_error)?;
        let s_requirement=ScalarPropertyRequirement::try_new(MOLAR_HENRY_SOLUBILITY_PROPERTY,
            MOLAR_HENRY_SOLUBILITY_DIMS,ScalarAdmissibility::StrictlyPositive).map_err(source_error)?;
        for layer in layers {
            poll(&mut cancelled)?;
            if layer.material.species!=species||layer.cells==0||!layer.thickness_m.is_finite()||layer.thickness_m<=0.0 {
                return Err(Error::Invalid("each transient layer must bind the same pure species and positive geometry"));
            }
            let point=selection_point(layer.material.context,temperature_k)?;
            let original=resolve_material_state_point(layer.material.material,&point,
                &[d_requirement.clone(),s_requirement.clone()],layer.material.selection.clone()).map_err(source_error)?;
            let d=original.property(MEMBRANE_DIFFUSIVITY_PROPERTY)
                .ok_or(Error::Invalid("resolved dense-phase diffusivity is missing"))?;
            let s=original.property(MOLAR_HENRY_SOLUBILITY_PROPERTY)
                .ok_or(Error::Invalid("resolved molar Henry solubility is missing"))?;
            admit_henry(d)?;admit_henry(s)?;
            coefficients.push(MembraneDiffusionLayer {thickness_m:layer.thickness_m,cells:layer.cells,
                diffusivity:MembraneDiffusivity::try_new(d.value_si()).map_err(GasEnergyError::from)?,
                solubility:MolarHenrySolubility::try_new(s.value_si()).map_err(GasEnergyError::from)?});
            bound_layers.push(LayerBinding {binding:layer.material.clone(),original});
        }
        let requirement=ScalarPropertyRequirement::try_new(MOLAR_INTERFACE_PERMEANCE_PROPERTY,
            MOLAR_PERMEANCE_DIMS,ScalarAdmissibility::NonNegative).map_err(source_error)?;
        let mut bound_interfaces=Vec::with_capacity(interfaces.len());
        let mut permeances=Vec::with_capacity(interfaces.len());
        for (i,binding) in interfaces.iter().enumerate() {
            poll(&mut cancelled)?;
            let Some(binding)=binding else {bound_interfaces.push(None);permeances.push(None);continue;};
            if binding.species!=species||&binding.interface.surface_a().material!=layers[i].material.material.id()
                ||&binding.interface.surface_b().material!=layers[i+1].material.material.id() {
                return Err(Error::Invalid("transient interlayer source must bind the same species and ordered adjacent material states"));
            }
            let point=selection_point(binding.context,temperature_k)?;
            let original=resolve_interface_state_point(binding.interface,&point,
                &[requirement.clone()],binding.selection.clone()).map_err(source_error)?;
            let property=original.property(MOLAR_INTERFACE_PERMEANCE_PROPERTY)
                .ok_or(Error::Invalid("resolved molar interface permeance is missing"))?;
            admit_henry(property)?;
            permeances.push(Some(MolarPermeance::try_new(property.value_si()).map_err(GasEnergyError::from)?));
            bound_interfaces.push(Some(InterfaceBinding {binding:binding.clone(),original}));
        }
        let model=MembraneDiffusion::new(area_m2,&coefficients,&permeances,&mut cancelled)?;
        let state=model.initial_state(initial_partial_pressures_pa,&mut cancelled)?;
        let flux=state.sample_fluxes(initial_reservoirs_pa,&mut cancelled)?;
        let bindings=Arc::new(Bindings {species,temperature_k,layers:bound_layers,interfaces:bound_interfaces});
        let sample=(&state,flux.face_partial_pressures_pa.as_slice());
        let support=bindings.admit([sample,sample],&mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(Self {bindings,state,support})
    }
    pub fn species(&self)->&SpeciesId {&self.bindings.species}
    pub fn temperature_k(&self)->f64 {self.bindings.temperature_k}
    /// The immutable numerical state, including actual retained sorbed amount.
    pub const fn state(&self)->&MembraneDiffusionState {&self.state}
    pub const fn source_support(&self)->&MembraneDiffusionSourceSupport {&self.support}
    /// Preserve the exact selected original claims while admitting actual old
    /// and new cell/face pressures. A failed source check or cancellation leaves
    /// this owner untouched; numerical workspace changes carry no inventory.
    pub fn step_candidate(&self,reservoirs_pa:[f64;2],duration_s:f64,controls:MembraneDiffusionSolve,
        work:&mut MembraneDiffusionWorkspace,mut cancelled:impl FnMut()->bool)
        ->Result<ResolvedMembraneDiffusionCandidate<'a>,Error>
    {
        poll(&mut cancelled)?;
        let old=self.state.sample_fluxes(reservoirs_pa,&mut cancelled)?;
        let candidate=self.state.step_candidate(reservoirs_pa,duration_s,controls,work,&mut cancelled)?;
        let (state,step)=candidate.into_parts();
        let support=self.bindings.admit([(&self.state,&old.face_partial_pressures_pa),
            (&state,&step.face_partial_pressures_pa)],&mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(ResolvedMembraneDiffusionCandidate {state:Self {bindings:self.bindings.clone(),state,support},step})
    }
}

impl Bindings<'_> {
    fn admit(&self,samples:[(&MembraneDiffusionState,&[[Option<f64>;2]]);2],
        cancelled:&mut impl FnMut()->bool)->Result<MembraneDiffusionSourceSupport,Error>
    {
        let mut layers=Vec::with_capacity(self.layers.len());
        for (i,binding) in self.layers.iter().enumerate() {
            poll(cancelled)?;
            let mut lo=f64::INFINITY;let mut hi:f64=0.0;
            for (state,faces) in samples {
                let range=&state.model().layer_cell_ranges()[i];
                for (local,&q) in state.partial_pressures_pa()[range.clone()].iter().enumerate() {
                    if local%128==0 {poll(cancelled)?;}
                    lo=lo.min(q);hi=hi.max(q);
                }
                // Within-layer traces are convex combinations of cells, so
                // only these exterior layer traces can extend its pressure range.
                for q in [faces[range.start][1],faces[range.end][0]].into_iter().flatten() {
                    lo=lo.min(q);hi=hi.max(q);
                }
            }
            let low=binding.original.at_query_point(binding.binding.material,
                &updated_point(binding.binding.context,self.temperature_k,lo)?).map_err(source_error)?;
            poll(cancelled)?;
            let high=if hi==lo {low.clone()}else {binding.original.at_query_point(binding.binding.material,
                &updated_point(binding.binding.context,self.temperature_k,hi)?).map_err(source_error)?};
            layers.push([low,high]);
        }
        let mut interfaces=Vec::with_capacity(self.interfaces.len());
        for (i,binding) in self.interfaces.iter().enumerate() {
            poll(cancelled)?;
            let Some(binding)=binding else {interfaces.push(None);continue;};
            let mut lo=f64::INFINITY;let mut hi:f64=0.0;
            for (state,faces) in samples {
                let face=state.model().layer_cell_ranges()[i].end;
                for q in faces[face] {
                    let q=q.ok_or(Error::Invalid("coincident seals leave an interface trace unknown; cannot qualify that interface source"))?;
                    lo=lo.min(q);hi=hi.max(q);
                }
            }
            let low=binding.original.at_query_point(binding.binding.interface,
                &updated_point(binding.binding.context,self.temperature_k,lo)?).map_err(source_error)?;
            poll(cancelled)?;
            let high=if hi==lo {low.clone()}else {binding.original.at_query_point(binding.binding.interface,
                &updated_point(binding.binding.context,self.temperature_k,hi)?).map_err(source_error)?};
            interfaces.push(Some([low,high]));
        }
        poll(cancelled)?;
        Ok(MembraneDiffusionSourceSupport {layers,interfaces})
    }
}

#[cfg(test)]
#[path="transient/tests.rs"]
mod tests;
