//! Original interface sources in the native nonlinear electrical closure.
use super::*;
use crate::circuit_enthalpy::material::CellConductionResponse;
use crate::circuit_enthalpy::StepCoefficients;
use crate::material_circuit::conductivity::{ConductivitySample, ConductorResponse,
    ConductionCoordinate};
use crate::material_circuit::conductivity::contact::{ContactAxis, ContactCoordinate,
    ContactDirection, ContactResponse};
use fs_qty::QuantitySpec;

/// Operating coordinates of one original electrical interface. Its material
/// card, A/B cells, area, source selection and heat split remain network-owned.
#[derive(Clone, Copy, Debug)]
pub struct ContactConductionResponse<'a> {
    pub id: u64,
    /// Source A-to-B orientation relative to the containing graph branch.
    pub direction: ContactDirection,
    /// At most one voltage, current and current-density role, with explicit
    /// signed or magnitude semantics. Every contact needs one descriptor.
    pub axes: &'a [ContactAxis<'a>],
}

impl<'a> ElectricalContactNetwork<'a> {
    fn contact_response_sources<'s>(
        &'s self, responses: &'s [ContactConductionResponse<'_>],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<ContactResponse<'s>>, ElectricalContactNetworkError> {
        use ElectricalContactNetworkError as Error;
        if responses.len() != self.contacts.len() || responses.len() > MAX_SERIES_CONTACTS {
            return Err(Error::Invalid("one response per original electrical contact is required"));
        }
        let mut seen = std::collections::BTreeSet::new();
        let mut sources = Vec::with_capacity(responses.len());
        for response in responses {
            if cancelled() { return Err(Error::Thermal(DepositionError::Cancelled)); }
            let binding = self.contacts.iter().find(|contact| contact.link == response.id)
                .ok_or(Error::Invalid("contact response has no original interface binding"))?;
            if !seen.insert(response.id) || response.axes.len() > 3 {
                return Err(Error::Invalid("duplicate or unbounded electrical contact response"));
            }
            let context = binding.point.context;
            for (index, axis) in response.axes.iter().enumerate() {
                if !context.axes().contains_key(axis.name)
                    || context.axis_quantities().get(axis.name)
                        != Some(&QuantitySpec::dimensional(axis.coordinate.dimensions()))
                    || response.axes[..index].iter().any(|other| other.name == axis.name
                        || other.coordinate.dimensions() == axis.coordinate.dimensions()) {
                    return Err(Error::Invalid("contact response requires unique existing SI V/I/J coordinates"));
                }
            }
            for (name, quantity) in context.axis_quantities() {
                if [ContactCoordinate::Voltage.dimensions(), ContactCoordinate::Current.dimensions(),
                    ContactCoordinate::CurrentDensity.dimensions(), ConductionCoordinate::ElectricField.dimensions()]
                    .contains(&quantity.dims())
                    && !response.axes.iter().any(|axis| axis.name == name.as_str()) {
                    return Err(Error::Invalid("contact electrical coordinates cannot remain frozen or imply a field"));
                }
            }
            sources.push(ContactResponse { id: response.id, card: binding.point.interface,
                direction: response.direction, axes: response.axes });
        }
        Ok(sources)
    }

    /// Pin the same complete bulk/interface source owners for native E/J and
    /// interface V/I/J iteration. Directional Fourier links retain their original
    /// principal claims; an empty directions slice selects scalar conduction.
    pub fn prepare_conductivity<'n>(
        &'n self, thermal: &DepositionSimulation<'_>, responses: &[CellConductionResponse<'_>],
        contacts: &[ContactConductionResponse<'_>], directions: &[DirectionalConductionLink],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedElectricalContactNetwork<'n, 'a>, ElectricalContactNetworkError> {
        use ElectricalContactNetworkError as Error;
        if cancelled() { return Err(Error::Thermal(DepositionError::Cancelled)); }
        self.routes(thermal)?;
        let bases = self.electrical_bases.iter().copied().collect::<BTreeMap<_, _>>();
        if bases.len() != self.electrical_bases.len() || bases.len() != responses.len()
            || responses.iter().any(|response| bases.get(&response.branch) != Some(&response.basis)) {
            return Err(Error::Invalid("conductivity must retain original electrical property bases"));
        }
        let sources = self.contact_response_sources(contacts, &mut cancelled)?;
        let bulk = self.bulk.prepare_conductivity_with_directional_conduction(
            thermal, responses, directions, &mut cancelled).map_err(Error::Bulk)?;
        let sample = bulk.resolve_conductivity_seeds(thermal, responses, &mut cancelled).map_err(Error::Bulk)?;
        let initial = self.resolve_contacts_with_response(thermal, sample, None, &sources, &mut cancelled)?;
        self.retain_contact_sources(thermal, bulk, initial, cancelled)
    }
}

impl PreparedElectricalContactNetwork<'_, '_> {
    pub(crate) fn conductivity_sources<'s>(
        &'s self, responses: &'s [CellConductionResponse<'_>],
        contacts: &'s [ContactConductionResponse<'_>], mut cancelled: impl FnMut() -> bool,
    ) -> Result<(Vec<ConductorResponse<'s>>, Vec<ContactResponse<'s>>), ElectricalContactNetworkError> {
        let bulk = self.bulk.conductivity_sources(responses, &mut cancelled)
            .map_err(ElectricalContactNetworkError::Bulk)?;
        let contacts = self.network.contact_response_sources(contacts, &mut cancelled)?;
        Ok((bulk, contacts))
    }

    /// Re-query original pinned bulk/contact sources at the actual trial thermal
    /// state and authored electrical seeds. The native electrical iteration must
    /// close these seeds before any thermal candidate may be accepted.
    pub fn resolve_conductivity_seeds(
        &self, thermal: &DepositionSimulation<'_>, responses: &[CellConductionResponse<'_>],
        contacts: &[ContactConductionResponse<'_>], mut cancelled: impl FnMut() -> bool,
    ) -> Result<ContactThermalCoefficients<MaterialNetworkSample>, ElectricalContactNetworkError> {
        if cancelled() { return Err(ElectricalContactNetworkError::Thermal(DepositionError::Cancelled)); }
        self.admit_geometry(thermal, &mut cancelled)?;
        let sources = self.network.contact_response_sources(contacts, &mut cancelled)?;
        let bulk = self.bulk.resolve_conductivity_seeds(thermal, responses, &mut cancelled)
            .map_err(ElectricalContactNetworkError::Bulk)?;
        self.network.resolve_contacts_with_response(thermal, bulk, Some(&self.pins), &sources, cancelled)
    }

    pub(crate) fn refresh_conductivity_transport(
        &self, thermal: &DepositionSimulation<'_>, electrical: &[ConductivitySample],
        coefficients: StepCoefficients<MaterialNetworkSample>, cancelled: impl FnMut() -> bool,
    ) -> Result<StepCoefficients<MaterialNetworkSample>, ElectricalContactNetworkError> {
        self.bulk.refresh_conductivity_transport(thermal, electrical, coefficients, cancelled)
            .map_err(ElectricalContactNetworkError::Bulk)
    }
}
