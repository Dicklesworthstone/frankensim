//! Original interface sources on the native explicit finite-cell timeline.
//!
//! Wrap the existing BULK-ONLY coefficient owner. Its fresh conductance map
//! receives the contact resistance once, before the native heat and CFL owner
//! consumes it. All other hooks retain their original owner and error types.
use super::{
    ContactInterfaceState, ContactLinkSample, MAX_MATERIAL_NETWORK_CONTACTS, MaterialContactLink,
    apply_material_contacts, resolve_contact,
};
use crate::ConductionError;
use fs_matdb::ClaimSelection;
use fs_thermal::enthalpy::{
    CellGeometry, CellId, ContactPatch, DepositedCell, DepositionError, DepositionObserver,
    DepositionSimulation, Event, HeatTransfer, MaterialId,
};
use std::collections::BTreeMap;

/// Preserve the original failing physical owner, including checkpoint errors.
#[derive(Clone, Debug, PartialEq)]
pub enum MaterialContactObserverError<E> {
    Thermal(DepositionError),
    Bulk(E),
    Contact(ConductionError),
    Binding(&'static str),
}
impl<E> From<DepositionError> for MaterialContactObserverError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl<E: core::fmt::Display> core::fmt::Display for MaterialContactObserverError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Thermal(error) => write!(f, "material contact thermal state: {error}"),
            Self::Bulk(error) => write!(f, "material contact bulk owner: {error}"),
            Self::Contact(error) => write!(f, "material contact source: {error}"),
            Self::Binding(reason) => write!(f, "material contact binding: {reason}"),
        }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for MaterialContactObserverError<E> {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Thermal(error) => Some(error),
            Self::Bulk(error) => Some(error),
            Self::Contact(error) => Some(error),
            Self::Binding(_) => None,
        }
    }
}

#[derive(Clone, Copy, PartialEq)]
struct CellIdentity {
    material: MaterialId,
    body: u64,
    segment: u64,
    mass_kg: f64,
    geometry: CellGeometry,
}
impl From<&DepositedCell> for CellIdentity {
    fn from(cell: &DepositedCell) -> Self {
        Self {
            material: cell.material_id,
            body: cell.body_id,
            segment: cell.segment_id,
            mass_kg: cell.mass_kg,
            geometry: cell.geometry,
        }
    }
}

#[derive(Clone)]
struct FixedContact {
    patch: ContactPatch,
    cells: [CellIdentity; 2],
}

/// Source-conditioned interface resistance in series with an existing observer.
///
/// Each selected link must receive a fresh explicit BULK-ONLY conductance from
/// `inner`; do not wrap an already contact-corrected owner. Untargeted links,
/// sources, reservoirs, pair transfers and interval-dependent work are forwarded.
/// No energy, clock, solver or thermal conductance is persisted into the cells.
/// The native timeline owns heat exchange, conservation, timestep bounds and
/// rollback. The underlying observer must obey its deep-Clone contract.
///
/// Geometry, mass and source identities remain bound. Reprepare after a contact
/// area, coverage, topology or material change. Temperature and phase follow the
/// actual caloric cells; other source coordinates stay explicitly supplied.
/// This retains the existing explicit coefficient-freezing approximation and
/// does not infer mechanical contact pressure, interfacial storage or spreading.
#[derive(Clone)]
pub struct MaterialContactObserver<'a, O> {
    inner: O,
    contacts: Vec<MaterialContactLink<'a>>,
    fixed: Vec<FixedContact>,
    samples: Vec<ContactLinkSample>,
    endpoint: Vec<ContactInterfaceState>,
}

impl<'a, O: DepositionObserver> MaterialContactObserver<'a, O> {
    /// Admit the original ordered sources and pin each selected claim. A cloned
    /// bulk observer is probed in begin/sources/conductances order; the retained
    /// observer's physical history and counters do not advance during attachment.
    pub fn new(
        inner: O,
        thermal: &DepositionSimulation<'_>,
        mut contacts: Vec<MaterialContactLink<'a>>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, MaterialContactObserverError<O::Error>> {
        checkpoint()?;
        if contacts.is_empty() || contacts.len() > MAX_MATERIAL_NETWORK_CONTACTS {
            return Err(MaterialContactObserverError::Binding(
                "a bounded nonempty contact set is required",
            ));
        }
        contacts.sort_by_key(|contact| contact.link);
        let mut probe = inner.clone();
        probe.begin().map_err(MaterialContactObserverError::Bulk)?;
        let _ = probe
            .sources(thermal, &mut checkpoint)
            .map_err(MaterialContactObserverError::Bulk)?;
        let mut bulk = probe
            .conductances(thermal, &mut checkpoint)
            .map_err(MaterialContactObserverError::Bulk)?;
        let samples = compose(thermal, &contacts, &mut bulk, &mut checkpoint)?;
        let mut fixed = Vec::with_capacity(contacts.len());
        let mut endpoint = Vec::with_capacity(contacts.len());
        for (contact, sample) in contacts.iter_mut().zip(&samples) {
            checkpoint()?;
            contact.selection = ClaimSelection::Pinned(sample.interface.receipt().selected);
            let patch = thermal.contact_patch(contact.link).ok_or(
                MaterialContactObserverError::Binding("admitted contact disappeared"),
            )?;
            let cells = [patch.link.a, patch.link.b]
                .map(|cell| CellIdentity::from(&thermal.cells()[&cell]));
            fixed.push(FixedContact { patch, cells });
            endpoint.push(ContactInterfaceState {
                link: sample.link,
                cells: sample.cells,
                phases: sample.phases,
                interface: sample.interface.clone(),
                area_m2: sample.area_m2,
            });
        }
        checkpoint()?;
        Ok(Self {
            inner,
            contacts,
            fixed,
            samples,
            endpoint,
        })
    }

    pub fn inner(&self) -> &O {
        &self.inner
    }
    pub fn into_inner(self) -> O {
        self.inner
    }
    /// Canonical link order, with the original source selections now pinned.
    pub fn bindings(&self) -> &[MaterialContactLink<'a>] {
        &self.contacts
    }
    /// Attachment evaluation, then the most recent starting-state coefficients.
    /// These are not end-state coefficients or an integrated heat receipt.
    pub fn last_samples(&self) -> &[ContactLinkSample] {
        &self.samples
    }
    /// Actual final caloric states and contact sources, independent of frozen G.
    pub fn endpoint_sources(&self) -> &[ContactInterfaceState] {
        &self.endpoint
    }

    fn validate_bindings(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), MaterialContactObserverError<O::Error>> {
        for (contact, fixed) in self.contacts.iter().zip(&self.fixed) {
            checkpoint()?;
            let patch = thermal.contact_patch(contact.link).ok_or(
                MaterialContactObserverError::Binding("bound thermal contact is inactive"),
            )?;
            // Total authored G may change: this observer replaces it with the
            // fresh bulk/contact series law. Physical geometry may not change.
            if patch.link.id != fixed.patch.link.id
                || patch.link.a != fixed.patch.link.a
                || patch.link.b != fixed.patch.link.b
                || patch.link.area_m2 != fixed.patch.link.area_m2
                || patch.covered_area_m2 != fixed.patch.covered_area_m2
            {
                return Err(MaterialContactObserverError::Binding(
                    "bound contact geometry or ordered sides changed",
                ));
            }
            for (side, id) in [patch.link.a, patch.link.b].into_iter().enumerate() {
                let cell =
                    thermal
                        .cells()
                        .get(&id)
                        .ok_or(MaterialContactObserverError::Binding(
                            "bound contact cell is inactive",
                        ))?;
                if CellIdentity::from(cell) != fixed.cells[side] {
                    return Err(MaterialContactObserverError::Binding(
                        "bound contact cell inventory or geometry changed",
                    ));
                }
            }
        }
        checkpoint()?;
        Ok(())
    }

    fn refresh_endpoint(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), MaterialContactObserverError<O::Error>> {
        self.validate_bindings(thermal, checkpoint)?;
        let mut endpoint = Vec::with_capacity(self.contacts.len());
        for contact in &self.contacts {
            checkpoint()?;
            endpoint.push(
                resolve_contact(thermal, contact).map_err(MaterialContactObserverError::Contact)?,
            );
        }
        checkpoint()?;
        self.endpoint = endpoint;
        Ok(())
    }
}

fn compose<E>(
    thermal: &DepositionSimulation<'_>,
    contacts: &[MaterialContactLink<'_>],
    bulk: &mut BTreeMap<u64, f64>,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<Vec<ContactLinkSample>, MaterialContactObserverError<E>> {
    let mut stopped = None;
    let result = apply_material_contacts(thermal, contacts, bulk, || match checkpoint() {
        Ok(()) => false,
        Err(error) => {
            stopped = Some(error);
            true
        }
    });
    if let Some(error) = stopped {
        return Err(MaterialContactObserverError::Thermal(error));
    }
    result.map_err(MaterialContactObserverError::Contact)
}

impl<O: DepositionObserver> DepositionObserver for MaterialContactObserver<'_, O> {
    type Error = MaterialContactObserverError<O::Error>;

    fn begin(&mut self) -> Result<(), Self::Error> {
        self.inner
            .begin()
            .map_err(MaterialContactObserverError::Bulk)
    }
    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.inner
            .sources(thermal, checkpoint)
            .map_err(MaterialContactObserverError::Bulk)
    }
    fn conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<u64, f64>, Self::Error> {
        self.validate_bindings(thermal, checkpoint)?;
        let mut bulk = self
            .inner
            .conductances(thermal, checkpoint)
            .map_err(MaterialContactObserverError::Bulk)?;
        let samples = compose(thermal, &self.contacts, &mut bulk, checkpoint)?;
        checkpoint()?;
        self.samples = samples;
        Ok(bulk)
    }
    fn reservoir_conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.inner
            .reservoir_conductances(thermal, checkpoint)
            .map_err(MaterialContactObserverError::Bulk)
    }
    fn heat_transfers(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<HeatTransfer>, Self::Error> {
        self.inner
            .heat_transfers(thermal, checkpoint)
            .map_err(MaterialContactObserverError::Bulk)
    }
    fn limit_timestep_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        proposed_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, Self::Error> {
        self.inner
            .limit_timestep_s(thermal, proposed_s, checkpoint)
            .map_err(MaterialContactObserverError::Bulk)
    }
    fn requires_step_interval(&self) -> bool {
        self.inner.requires_step_interval()
    }
    fn step_sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        start_s: f64,
        dt_s: f64,
        end_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.inner
            .step_sources(thermal, start_s, dt_s, end_s, checkpoint)
            .map_err(MaterialContactObserverError::Bulk)
    }
    fn event(
        &mut self,
        event: &Event,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), Self::Error> {
        self.inner
            .event(event, thermal)
            .map_err(MaterialContactObserverError::Bulk)?;
        self.refresh_endpoint(thermal, &mut || Ok(()))
    }
    fn step(
        &mut self,
        start_s: f64,
        duration_s: f64,
        old: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.inner
            .step(start_s, duration_s, old, thermal, checkpoint)
            .map_err(MaterialContactObserverError::Bulk)?;
        self.refresh_endpoint(thermal, checkpoint)
    }
}
