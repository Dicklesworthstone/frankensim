//! Keep selected interface and bulk claims through one composed physical solve.
//! Values remain live: no temperature, resistance or conductivity is cached.
use super::*;
use crate::circuit_enthalpy::material::PreparedMaterialNetwork;

/// One borrowed electrical/thermal contact network with fixed source choices.
/// The existing bulk preparation and interface evaluator own all material queries.
/// Source-domain failure cannot select another claim during Newton iteration.
/// Keep this across steps when source continuity is required; explicit preparation
/// at a new starting state permits new policy selections, not automatic fallback.
pub struct PreparedElectricalContactNetwork<'n,'a> {
    network:&'n ElectricalContactNetwork<'a>,
    bulk:PreparedMaterialNetwork<'n,'a>,
    pins:BTreeMap<u64,MaterialPropertySelection>,
}
impl PreparedElectricalContactNetwork<'_,'_> {
    /// Fresh material snapshots and contacts at the actual thermal trial state.
    /// All original identity, geometry, phase, quantity and source checks remain.
    pub fn resolve(&self,thermal:&DepositionSimulation<'_>,mut cancelled:impl FnMut()->bool)
        ->Result<ContactThermalCoefficients<MaterialNetworkSample>,ElectricalContactNetworkError> {
        if cancelled() {return Err(ElectricalContactNetworkError::Thermal(DepositionError::Cancelled));}
        let bulk=self.bulk.resolve(thermal,&mut cancelled).map_err(ElectricalContactNetworkError::Bulk)?;
        self.network.resolve_contacts(thermal,bulk,Some(&self.pins),cancelled)
    }
}
impl<'a> ElectricalContactNetwork<'a> {
    /// Select each original bulk and interface claim once at the supplied state.
    /// Later evaluations query those pins through the existing material owners,
    /// including when competing claims overlap. Unsupported pinned states refuse.
    /// Cancellation returns no partially prepared object or physical endpoint.
    pub fn prepare<'n>(&'n self,thermal:&DepositionSimulation<'_>,mut cancelled:impl FnMut()->bool)
        ->Result<PreparedElectricalContactNetwork<'n,'a>,ElectricalContactNetworkError> {
        use ElectricalContactNetworkError as Error;
        if cancelled() {return Err(Error::Thermal(DepositionError::Cancelled));}
        self.routes(thermal)?;
        let bulk=self.bulk.prepare_with_bases(thermal,self.electrical_bases,&mut cancelled).map_err(Error::Bulk)?;
        let sample=bulk.resolve(thermal,&mut cancelled).map_err(Error::Bulk)?;
        let initial=self.resolve_contacts(thermal,sample,None,&mut cancelled)?;
        let mut pins=BTreeMap::new();
        for contact in initial.contacts {
            if cancelled() {return Err(Error::Thermal(DepositionError::Cancelled));}
            let claim=contact.source().resolved().property(ELECTRICAL_CONTACT_RESISTANCE_PROPERTY)
                .ok_or(Error::Invalid("initial electrical contact property is absent"))?.answer().receipt.selected;
            pins.insert(contact.id(),MaterialPropertySelection::PinnedByProperty(vec![
                (ELECTRICAL_CONTACT_RESISTANCE_PROPERTY.to_owned(),claim)]));
        }
        if cancelled() {return Err(Error::Thermal(DepositionError::Cancelled));}
        Ok(PreparedElectricalContactNetwork {network:self,bulk,pins})
    }
}
