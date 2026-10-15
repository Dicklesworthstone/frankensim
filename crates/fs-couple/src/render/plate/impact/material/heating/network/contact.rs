//! Contact area, thermal interface and ambient coverage share one physical link.
use super::*;
use fs_conduction::interface::{InterfaceResistance, SeriesThermalResistance};
use fs_matdb::{QueryPoint, SelectionPolicy};
use fs_thermal::enthalpy::{ContactUpdate, ThermalLink};

/// One source-bound mechanical contact owns one ordered thermal link.
/// The interface resistance is queried from that contact's ORIGINAL interface
/// card and interface-law state. It is not another bulk-material property.
#[derive(Clone, Debug, PartialEq)]
pub struct ContactNetworkLink {
    pub contact: String,
    pub link_id: u64,
    pub cells: [CellId; 2],
    pub selection: SelectionPolicy,
}

/// Exact source and geometry used by this sample's network heat transfer.
#[derive(Clone, Debug, PartialEq)]
pub struct ContactNetworkReport {
    pub contact: String,
    pub link_id: u64,
    pub area_m2: f64,
    pub interface: InterfaceResistance,
    /// None means physically separated, not a large invented resistance.
    pub resistance: Option<SeriesThermalResistance>,
}
#[derive(Debug)]
pub struct ContactNetworkFrame {
    pub network: NetworkImpactFrame,
    pub contacts: Vec<ContactNetworkReport>,
}

impl<'m> HeatedImpactSystem<'m> {
    /// Couple every selected mechanical footprint to the shared thermal network.
    ///
    /// Resolve R'' through the original `InterfaceResistance` owner at the exact
    /// retained interface point, and form R''/area through its existing series
    /// operator. Both thermal curves must identify the original ordered bulk
    /// cards. The existing physical circle/ellipse accessor supplies area; a
    /// pressure-equivalent finite-gap footprint is not admitted as transfer area.
    ///
    /// Replace all changed contacts together before advancing the network. An
    /// open contact exposes both covered areas for the original convection law;
    /// recontact covers those areas again without resetting any enthalpy. Existing
    /// unrelated links remain. Link IDs cannot be reused for different partners,
    /// asymmetric footprints cannot be silently converted, and an observer may
    /// not override the selected links' source-derived conductances.
    ///
    /// Areas and interface coefficients are held at the proposed mechanical
    /// endpoint for this thermal interval. This first-order composition is NOT
    /// event-resolved area integration, thermal expansion, spreading resistance,
    /// flash-temperature resolution or continuous loaded thermomechanical work.
    /// Material-source refresh remains unloaded-only. Every resulting state,
    /// source, observer and contact topology is accepted together or discarded.
    pub fn step_with_contact_network<'n, O: DepositionObserver>(
        &mut self, thermal: &mut DepositionSimulation<'n>, observer: &mut O,
        bindings: &[ReservoirCell], links: &[ContactNetworkLink], config: StepConfig,
        forces: &[f64], gate: &CancelGate,
        resolve: impl FnMut(&MaterialContact, &[CaloricReservoir<'m>])
            -> Result<Option<SourcePointContact>, ImpactError>,
    ) -> Result<ContactNetworkFrame, ImpactError>
    where O::Error: core::fmt::Display {
        if links.len()>32 {return Err(invalid("physical thermal contacts exceed the 32-contact ceiling"));}
        let (network, contacts)=self.step_with_network_preparation(thermal,observer,bindings,
            config,forces,gate,resolve,|trial,network|prepare(trial,network,bindings,links,gate))?;
        Ok(ContactNetworkFrame {network,contacts})
    }
}

fn prepare<'n>(trial:&ImpactThermalTrial<'_>,thermal:&DepositionSimulation<'n>,
    bindings:&[ReservoirCell],links:&[ContactNetworkLink],gate:&CancelGate)
    ->Result<(DepositionSimulation<'n>,Vec<ContactNetworkReport>,BTreeSet<u64>),ImpactError> {
    let mut names=BTreeSet::new();let mut ids=BTreeSet::new();
    let mut changes=Vec::with_capacity(links.len());let mut reports=Vec::with_capacity(links.len());
    // Stable link order also fixes area-subtraction roundoff for several patches.
    let mut ordered:Vec<_>=links.iter().collect();ordered.sort_by_key(|link|link.link_id);
    for link in ordered {
        if gate.is_requested() {return Err(ImpactError::Cancelled);}
        if !names.insert(link.contact.as_str()) || !ids.insert(link.link_id) || link.cells[0]==link.cells[1] {
            return Err(invalid("mechanical thermal contacts need unique keys/links and distinct partners"));
        }
        let source=trial.contact_source(&link.contact)?;
        let normal=source.model().normal_state();let elastic=normal.elastic();
        let bulk=[elastic.surface_a_state(),elastic.surface_b_state()];
        for (side,&cell) in link.cells.iter().enumerate() {
            let binding=bindings.iter().find(|b|b.cell==cell).ok_or_else(||invalid("unbound physical contact cell"))?;
            let input=trial.reservoirs().iter().find(|r|r.key==binding.reservoir)
                .ok_or_else(||invalid("unbound physical contact reservoir"))?;
            if input.thermal.curve.material_card_identity()!=bulk[side].resolved().card_identity() {
                return Err(invalid("thermal link partners do not match the ordered contact material cards"));
            }
        }
        let mut point=QueryPoint::new();
        for (axis,value) in normal.interface_state().query_point() {
            point=match normal.interface_state().axis_quantities().get(axis) {
                Some(&quantity)=>point.with_quantity(axis,quantity,*value),
                None=>point.with(axis,*value),
            }.map_err(|e|ImpactError::Owner(e.to_string()))?;
        }
        let interface=InterfaceResistance::from_card(&link.contact,elastic.interface().card(),&point,link.selection)
            .map_err(|e|ImpactError::Owner(e.to_string()))?;
        let area=trial.contact_patch_area_m2(&link.contact)?;
        let resistance=if area==0.0 {None} else {
            Some(SeriesThermalResistance::new(vec![interface.term_for_area(&link.contact,area)
                .map_err(|e|ImpactError::Owner(e.to_string()))?]).map_err(|e|ImpactError::Owner(e.to_string()))?)
        };
        let expected=thermal.links().get(&link.link_id).copied();
        if expected.is_some_and(|old|old.a!=link.cells[0] || old.b!=link.cells[1]) {
            return Err(invalid("physical contact cannot reuse a thermal link with different ordered partners"));
        }
        let replacement=if let Some(resistance)=&resistance {
            let conductance=resistance.budget().value_k_per_w.recip();
            if !positive(conductance) {return Err(invalid("physical contact conductance is unrepresentable"));}
            Some(ThermalLink {id:link.link_id,a:link.cells[0],b:link.cells[1],area_m2:area,conductance_w_k:conductance})
        }else{None};
        changes.push(ContactUpdate {id:link.link_id,expected,replacement});
        reports.push(ContactNetworkReport {contact:link.contact.clone(),link_id:link.link_id,
            area_m2:area,interface,resistance});
    }
    let next=thermal.replace_contacts(&changes,||if gate.is_requested(){Err(DepositionError::Cancelled)}else{Ok(())})
        .map_err(|e|if gate.is_requested(){ImpactError::Cancelled}else{ImpactError::Owner(e.to_string())})?;
    Ok((next,reports,ids))
}
