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
/// Both temperatures must be equal; this convenience API never averages bodies.
/// All source, density, phase and context admission is shared with
/// [`resolve_contact_at_temperatures`]. The common temperature is also the
/// interface-law temperature, preserving the original isothermal behavior.
pub fn resolve_isothermal_contact(contact:&MaterialContact,cards:&PointContactMaterialCards<'_>,
    reference:&QueryPoint,stores:[&CaloricReservoir<'_>;2])->Result<SourcePointContact,ImpactError> {
    if stores[0].state().temperature_k()!=stores[1].state().temperature_k() {
        return Err(invalid("isothermal contact refresh cannot average different body temperatures"));
    }
    resolve_contact_at_temperatures(contact,cards,reference,stores,stores[0].state().temperature_k())
}

/// Resolve each contact material at its own accepted reservoir temperature.
///
/// The ordered stores must belong to the exact original bulk cards and remain
/// fully solid with unchanged consistent density. `reference` preserves all
/// nonthermal coordinates and schemas. The explicit `interface_temperature_k`
/// is used ONLY for the interface-property/model point; it does not replace
/// either bulk temperature. It may be a declared surface trace, another thermal
/// owner's interface state, or a prescribed condition. The caller must supply
/// and justify that choice. No arithmetic/effusivity average is silently chosen.
///
/// This returns the same SourcePointContact as the generic material resolver.
/// Use it inside `step_with_material_feedback` to retain atomic source, motion,
/// heat and clock acceptance. That runtime refreshes only unloaded contacts and
/// holds the source throughout compression. This function does not implement
/// loaded thermoelastic work, thermal expansion, plasticity or moving melt.
/// Changing card assignments or body inertia is a separate reconfiguration.
pub fn resolve_contact_at_temperatures(contact:&MaterialContact,cards:&PointContactMaterialCards<'_>,
    reference:&QueryPoint,stores:[&CaloricReservoir<'_>;2],interface_temperature_k:f64)
    ->Result<SourcePointContact,ImpactError> {
    let old=contact.source().model().normal_state().elastic();
    let originals=[old.surface_a_state(),old.surface_b_state()];
    let materials=[cards.surface_a,cards.surface_b];
    if !positive(interface_temperature_k) {
        return Err(invalid("contact interface temperature must be finite positive kelvin"));
    }
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
    // The old states may already have distinct temperatures. Compare every
    // retained coordinate in each role rather than assuming one old bulk point.
    for original in originals {
        if reference.axes().len()!=original.resolved().query_point().len()
            || !reference.axes().contains_key("T")
            || original.resolved().query_point().iter().any(|(axis,value)|axis!="T"
                && reference.axes().get(axis).is_none_or(|v|v.to_bits()!=value.to_bits()))
            || reference.axis_quantities()!=original.resolved().axis_quantities() {
            return Err(invalid("thermal refresh must retain every nonthermal coordinate and axis schema"));
        }
    }
    let update=|t|reference.clone().with("T",t).map_err(|e|ImpactError::Owner(e.to_string()));
    let points=[update(stores[0].state().temperature_k())?,update(stores[1].state().temperature_k())?];
    let interface_point=update(interface_temperature_k)?;
    let source=cards.resolve_at_temperatures([&points[0],&points[1]],&interface_point,
        contact.source().geometry(),contact.source().envelope())
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
