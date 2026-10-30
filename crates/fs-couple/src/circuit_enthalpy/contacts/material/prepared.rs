//! Keep selected interface and bulk claims through one composed physical solve.
//! Values remain live: no temperature, resistance or conductivity is cached.
use super::*;
use crate::circuit_enthalpy::material::{DirectionalConductionLink, PreparedMaterialNetwork};

mod extension;
mod field_storage;

/// One borrowed electrical/thermal contact network with fixed source choices.
/// The existing bulk preparation and interface evaluator own all material queries.
/// Source-domain failure cannot select another claim during Newton iteration.
/// Keep this across steps when source continuity is required; explicit preparation
/// at a new starting state permits new policy selections, not automatic fallback.
pub struct PreparedElectricalContactNetwork<'n,'a> {
    network:&'n ElectricalContactNetwork<'a>,
    bulk:PreparedMaterialNetwork<'n,'a>,
    pins:BTreeMap<u64,MaterialPropertySelection>,
    geometry:BTreeMap<u64,([CellId;2],f64)>,
}
impl PreparedElectricalContactNetwork<'_,'_> {
    /// Original principal source roles, material frames and link directions.
    pub fn directional_conduction(&self) -> impl Iterator<Item = &DirectionalConductionLink> {
        self.bulk.directional_conduction()
    }

    /// Original bulk heat destinations, without a second caller-authored map.
    pub fn resistor_destinations(&self) -> Vec<ResistorCell> {
        self.bulk.resistor_destinations()
    }

    /// Original two-sided interface partitions in stable contact-ID order.
    /// Thermal links remain the authority for side order and active cells.
    pub fn heat_destinations(&self, thermal: &DepositionSimulation<'_>)
        -> Result<Vec<ContactHeatCells>, ElectricalContactNetworkError> {
        self.network.heat_destinations(thermal)
    }

    /// Fresh material snapshots and contacts at the actual thermal trial state.
    /// All original identity, geometry, phase, quantity and source checks remain.
    pub fn resolve(&self,thermal:&DepositionSimulation<'_>,mut cancelled:impl FnMut()->bool)
        ->Result<ContactThermalCoefficients<MaterialNetworkSample>,ElectricalContactNetworkError> {
        if cancelled() {return Err(ElectricalContactNetworkError::Thermal(DepositionError::Cancelled));}
        self.admit_geometry(thermal,&mut cancelled)?;
        let bulk=self.bulk.resolve(thermal,&mut cancelled).map_err(ElectricalContactNetworkError::Bulk)?;
        self.network.resolve_contacts(thermal,bulk,Some(&self.pins),cancelled)
    }

    /// Resolve original electrical contacts and bulk sources while another
    /// physical owner supplies the remaining native thermal conductances.
    /// The bulk owner enforces a complete, disjoint link partition.
    pub fn resolve_with_link_conductances(
        &self, thermal: &DepositionSimulation<'_>,
        additional_links: &BTreeMap<u64, f64>, mut cancelled: impl FnMut() -> bool,
    ) -> Result<ContactThermalCoefficients<MaterialNetworkSample>, ElectricalContactNetworkError> {
        if cancelled() { return Err(ElectricalContactNetworkError::Thermal(DepositionError::Cancelled)); }
        self.admit_geometry(thermal, &mut cancelled)?;
        let bulk = self.bulk.resolve_with_link_conductances(thermal, additional_links, &mut cancelled)
            .map_err(ElectricalContactNetworkError::Bulk)?;
        self.network.resolve_contacts(thermal, bulk, Some(&self.pins), cancelled)
    }
}
impl<'a> ElectricalContactNetwork<'a> {
    /// Select each original bulk and interface claim once at the supplied state.
    /// Later evaluations query those pins through the existing material owners,
    /// including when competing claims overlap. Unsupported pinned states refuse.
    /// Cancellation returns no partially prepared object or physical endpoint.
    pub fn prepare<'n>(&'n self,thermal:&DepositionSimulation<'_>,cancelled:impl FnMut()->bool)
        ->Result<PreparedElectricalContactNetwork<'n,'a>,ElectricalContactNetworkError> {
        self.prepare_with_optional_links(thermal, &[], None, cancelled)
    }

    /// Retain the same original electrical and interface source selections for
    /// a heat assembly with additional, independently owned thermal links.
    /// Ordinary resolution still requires complete Fourier coverage.
    pub fn prepare_with_link_conductances<'n>(
        &'n self, thermal: &DepositionSimulation<'_>,
        additional_links: &BTreeMap<u64, f64>, cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedElectricalContactNetwork<'n, 'a>, ElectricalContactNetworkError> {
        self.prepare_with_optional_links(thermal, &[], Some(additional_links), cancelled)
    }

    /// Source principal conductivities in each declared link direction while
    /// retaining the original electrical contact laws and heat destinations.
    /// Direction descriptors are data; all values are re-queried at trial T.
    pub fn prepare_with_directional_conduction<'n>(
        &'n self, thermal: &DepositionSimulation<'_>,
        directions: &[DirectionalConductionLink], cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedElectricalContactNetwork<'n, 'a>, ElectricalContactNetworkError> {
        self.prepare_with_optional_links(thermal, directions, None, cancelled)
    }

    /// Combine oriented bulk links and disjoint links owned by other physical
    /// mechanisms, such as finite-fluid convection, on the same source timeline.
    pub fn prepare_with_directional_conduction_and_link_conductances<'n>(
        &'n self, thermal: &DepositionSimulation<'_>,
        directions: &[DirectionalConductionLink], additional_links: &BTreeMap<u64, f64>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedElectricalContactNetwork<'n, 'a>, ElectricalContactNetworkError> {
        self.prepare_with_optional_links(thermal, directions, Some(additional_links), cancelled)
    }

    fn prepare_with_optional_links<'n>(
        &'n self, thermal: &DepositionSimulation<'_>,
        directions: &[DirectionalConductionLink],
        additional_links: Option<&BTreeMap<u64, f64>>, mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedElectricalContactNetwork<'n, 'a>, ElectricalContactNetworkError> {
        use ElectricalContactNetworkError as Error;
        if cancelled() {return Err(Error::Thermal(DepositionError::Cancelled));}
        self.routes(thermal)?;
        let bulk = match additional_links {
            Some(links) => self.bulk.prepare_with_directional_conduction_and_link_conductances(
                thermal, self.electrical_bases, directions, links, &mut cancelled,
            ),
            None => self.bulk.prepare_with_directional_conduction(
                thermal, self.electrical_bases, directions, &mut cancelled),
        }.map_err(Error::Bulk)?;
        let sample = match additional_links {
            Some(links) => bulk.resolve_with_link_conductances(thermal, links, &mut cancelled),
            None => bulk.resolve(thermal, &mut cancelled),
        }.map_err(Error::Bulk)?;
        let initial=self.resolve_contacts(thermal,sample,None,&mut cancelled)?;
        let mut pins=BTreeMap::new();
        let mut geometry=BTreeMap::new();
        for contact in initial.contacts {
            if cancelled() {return Err(Error::Thermal(DepositionError::Cancelled));}
            let link=thermal.links().get(&contact.id())
                .ok_or(Error::Invalid("initial electrical contact link is absent"))?;
            geometry.insert(contact.id(),([link.a,link.b],link.area_m2));
            let claim=contact.source().resolved().property(ELECTRICAL_CONTACT_RESISTANCE_PROPERTY)
                .ok_or(Error::Invalid("initial electrical contact property is absent"))?.answer().receipt.selected;
            pins.insert(contact.id(),MaterialPropertySelection::PinnedByProperty(vec![
                (ELECTRICAL_CONTACT_RESISTANCE_PROPERTY.to_owned(),claim)]));
        }
        if cancelled() {return Err(Error::Thermal(DepositionError::Cancelled));}
        Ok(PreparedElectricalContactNetwork {network:self,bulk,pins,geometry})
    }
}
