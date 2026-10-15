//! Thermal state -> source contact -> motion -> contact heat, without resets.
//! Constitutive coefficients are refreshed only while a contact is unloaded and
//! held during compression. Continuous loaded thermoelastic work is not modeled.
use super::*;
use super::super::{MaterialContact,refresh::ContactSourceUpdate};
use super::transport::{CoupledThermalStep,ImpactThermalTrial,TransportedImpactFrame};
use fs_contact::normal_patch::source::{SourcePointContact,state::PointContactMaterialCards};
use fs_matdb::QueryPoint;

/// Accepted coupled tick and explicit sample-and-hold constitutive work counts.
#[derive(Debug)]
pub struct MaterialFeedbackFrame<R> {
    pub thermal:TransportedImpactFrame<R>,
    /// Contacts whose complete source changed at the start of this tick.
    pub refreshed_contacts:usize,
    /// Compressed contacts held at their previous admitted source state.
    pub held_loaded_contacts:usize,
}

impl<'a> HeatedImpactSystem<'a> {
    /// Refresh unloaded contacts from accepted thermal data, then advance the
    /// existing coupled mechanical/thermal step. Both operations publish together.
    ///
    /// `resolve` sees the original contact and all current finite stores. Return
    /// None for an unchanged/unbound contact. It is called only when indentation
    /// is exactly zero. Loaded contacts retain their admitted source until a later
    /// separated sample; the report makes that hold explicit. This is a declared
    /// sample-and-hold approximation, not continuous loaded thermal feedback.
    /// Callers must bound within-contact changes and material validity separately.
    ///
    /// The general thermal callback is the SAME joint carrier path as
    /// `step_with_coupled_transport`, so conduction/cooling and internal contact
    /// heat keep their original owners and energy accounting. Source changes
    /// preserve geometry, inertia, state, memory and time. Failed numerical,
    /// source or thermal work (including late cancellation) also restores the
    /// previous contact sources and loss-path selection. No heat, momentum or
    /// clock is reset or repaired. Callback side effects/panics cannot be undone;
    /// callbacks must propose rather than mutate external physical owners.
    pub fn step_with_material_feedback<R>(&mut self,forces:&[f64],gate:&CancelGate,
        mut resolve:impl FnMut(&MaterialContact,&[CaloricReservoir<'a>])
            ->Result<Option<SourcePointContact>,ImpactError>,
        transport:impl FnMut(ImpactThermalTrial<'_>)->Result<CoupledThermalStep<R>,ImpactError>,
    )->Result<MaterialFeedbackFrame<R>,ImpactError> {
        if gate.is_requested() {return Err(ImpactError::Cancelled);}
        let mut updates=Vec::new();let mut held_loaded_contacts=0;
        for (i,c) in self.mechanics.material_contacts().iter().enumerate() {
            if gate.is_requested() {return Err(ImpactError::Cancelled);}
            let x=self.mechanics.state();
            let (d,_)=c.kinematics(x,x).map_err(|source|ImpactError::MaterialContact {contact:i,source})?;
            if d>0.0 {held_loaded_contacts+=1;continue;}
            if let Some(source)=resolve(c,&self.stores)? {
                if &source!=c.source() {updates.push(ContactSourceUpdate {contact:c.key().into(),source});}
            }
        }
        let refreshed_contacts=updates.len();
        let staged=self.mechanics.stage_contact_refresh(updates,gate)?;
        let previous=staged.map(|next|self.mechanics.swap_contact_refresh(next));
        let result=self.step_with_coupled_transport(forces,gate,transport);
        match result {
            Ok(thermal)=>Ok(MaterialFeedbackFrame {thermal,refreshed_contacts,held_loaded_contacts}),
            Err(error)=>{
                if let Some(previous)=previous {self.mechanics.swap_contact_refresh(previous);}
                Err(error)
            }
        }
    }
}

/// Resolve a same-temperature contact from the exact two thermal/material owners.
///
/// Both reservoirs must be fully solid and their charts must belong to the
/// matching ordered bulk cards. Both temperatures must be equal: the existing
/// isotropic contact binder requires a common physical state, and this adapter
/// never averages unequal bodies. `reference` supplies every nonthermal source
/// coordinate and its exact schema; only T is replaced with accepted temperature.
///
/// This fixed-geometry/fixed-inertia path refuses changed density and a mismatch
/// between the caloric density and the bulk-contact density. Thermal expansion,
/// changed inertial geometry and phase-change mechanics require their own owners.
/// Source cards and interface identity must match the already bound contact;
/// changing a material assignment is a separate explicit reconfiguration.
/// The same reusable source resolver also serves nonthermal generic consumers.
pub fn resolve_isothermal_contact(contact:&MaterialContact,cards:&PointContactMaterialCards<'_>,
    reference:&QueryPoint,stores:[&CaloricReservoir<'_>;2])->Result<SourcePointContact,ImpactError> {
    let old=contact.source().model().normal_state().elastic();
    let originals=[old.surface_a_state(),old.surface_b_state()];
    let materials=[cards.surface_a,cards.surface_b];
    if cards.interface.content_hash()!=old.interface().card_identity() {
        return Err(invalid("thermal refresh must retain the original interface card"));
    }
    for i in 0..2 {
        let state=stores[i].state();
        if !state.is_fully_solid() || stores[i].curve().material_card_identity()!=materials[i].content_hash()
            || originals[i].resolved().card_identity()!=materials[i].content_hash() {
            return Err(invalid("thermal contact requires both original solid material/chart bindings"));
        }
    }
    if stores[0].state().temperature_k()!=stores[1].state().temperature_k() {
        return Err(invalid("isothermal contact refresh cannot average different body temperatures"));
    }
    if reference.axes().len()!=old.state_point().len()
        || !reference.axes().contains_key("T")
        || old.state_point().iter().any(|(axis,value)|axis!="T" && reference.axes().get(axis)!=Some(value))
        || reference.axis_quantities()!=old.surface_a_state().resolved().axis_quantities() {
        return Err(invalid("thermal refresh must retain every nonthermal coordinate and axis schema"));
    }
    let point=reference.clone().with("T",stores[0].state().temperature_k())
        .map_err(|e|ImpactError::Owner(e.to_string()))?;
    let source=cards.resolve(&point,contact.source().geometry(),contact.source().envelope())
        .map_err(|e|ImpactError::Owner(e.to_string()))?;
    let new=source.model().normal_state().elastic();
    for (i,state) in [new.surface_a_state(),new.surface_b_state()].into_iter().enumerate() {
        if state.density_kg_m3()!=originals[i].density_kg_m3()
            || state.density_kg_m3()!=stores[i].state().bulk_density_kg_m3() {
            return Err(invalid("thermal contact refresh requires unchanged consistent bulk/caloric density"));
        }
    }
    Ok(source)
}
