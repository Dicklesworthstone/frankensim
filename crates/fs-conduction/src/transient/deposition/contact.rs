//! State-dependent contacts in series with existing finite-cell bulk transport.
//! Contact resistance, series composition, caloric storage and the thermal
//! residual retain their original owners. This module supplies coefficients,
//! not another thermal or electrical time integrator.
use std::collections::BTreeMap;
use super::{CellId, DepositionSimulation, MaterialId, ThermalLaw};
use crate::{ConductionError, InterfaceResistance, SeriesThermalResistance,
    ThermalResistanceTerm, ResistanceUncertainty, AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS};
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState};
use fs_material::thermal::interface::ThermalInterfacePoint;
use fs_matdb::ClaimSelection;
use fs_qty::QuantitySpec;

mod observer;
pub use observer::{MaterialContactObserver, MaterialContactObserverError};
/// Prepared original interface sources for native implicit chain/graph transport.
pub mod implicit;

/// Maximum contacts queried during one coefficient evaluation.
pub const MAX_MATERIAL_NETWORK_CONTACTS: usize = 256;

/// Original caloric owner of one ordered side of a thermal contact.
#[derive(Clone, Copy)]
pub struct ContactCell<'a> {
    /// Existing active cell, not electrical branch order.
    pub cell: CellId,
    /// Original entry in the thermal law table.
    pub material_id: MaterialId,
    /// SAME borrowed caloric instance used by that entry.
    pub curve: &'a EquilibriumEnthalpyPhaseCurve,
}

/// One existing link with an explicit ordered interface law.
#[derive(Clone, Copy)]
pub struct MaterialContactLink<'a> {
    /// Existing link. Its A/B order and physical area remain authoritative.
    pub link: u64,
    /// Caloric owners, in the link AND interface card's A/B order.
    pub sides: [ContactCell<'a>; 2],
    /// Same bulk cards and full two-sided thermal/phase source conditions.
    pub point: ThermalInterfacePoint<'a>,
    /// Original selection policy or exact pinned claim.
    pub selection: ClaimSelection,
}

/// An original thermal interface's pressure supplied by its mechanical owner.
/// Link identity, ordered materials, source selection and thermal geometry stay
/// with the existing contact. Only the explicitly named pressure axis evolves.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalContactPressure<'a> {
    /// Original material thermal-link identity.
    pub link: u64,
    /// Existing static-pressure role, distinct from both thermal states.
    pub pressure_axis: &'a str,
    /// Finite nonnegative compression in pascals.
    pub pressure_pa: f64,
}

/// Current geometry of one original planar, two-half-cell thermal path.
/// A mechanical owner supplies both distances and the common interface area;
/// these values never change the reference caloric inventory or source context.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalLinkMetric {
    /// Original native thermal-link identity, retaining its A/B order.
    pub link: u64,
    /// Current positive interface area [m2].
    pub area_m2: f64,
    /// Current positive Fourier distances in original A/B order [m].
    pub distances_m: [f64; 2],
}

impl ThermalLinkMetric {
    /// Admit a finite, noncollapsed path. The supplying mechanical owner must
    /// separately establish correspondence and applicability of the reduction.
    pub fn validate(&self) -> Result<(), ConductionError> {
        if !self.area_m2.is_finite() || self.area_m2 <= 0.0
            || self.distances_m.iter().any(|d| !d.is_finite() || *d <= 0.0) {
            return Err(invalid("current thermal area and both half-paths must be positive finite values"));
        }
        Ok(())
    }
}

/// Complete contact result for one trial, retained with its bulk source payload.
#[derive(Clone, Debug, PartialEq)]
pub struct ContactLinkSample {
    /// Original link identity.
    pub link: u64,
    /// Actual ordered cells and endpoint phase states.
    pub cells: [CellId; 2],
    pub phases: [EquilibriumPhaseState; 2],
    /// Original contact-property result, including its full query receipt.
    pub interface: InterfaceResistance,
    /// Existing physical contact area [m2]. No electrical-area substitution.
    pub area_m2: f64,
    /// Bulk-only conductance supplied by the original material/transport owner.
    pub bulk_conductance_w_k: f64,
    /// Bulk and contact series result, or zero for explicit isolation [W/K].
    pub conductance_w_k: f64,
    /// Exact original resistance-owner report. None only at explicit isolation.
    /// The bulk term is declared; its original material receipts remain in the
    /// caller's bulk payload. This does not fabricate a new bulk-source receipt.
    pub series: Option<SeriesThermalResistance>,
}

/// Source state at an actual thermal endpoint, independent of any held bulk G.
#[derive(Clone, Debug, PartialEq)]
pub struct ContactInterfaceState {
    pub link: u64,
    pub cells: [CellId; 2],
    pub phases: [EquilibriumPhaseState; 2],
    pub interface: InterfaceResistance,
    pub area_m2: f64,
}

/// Add contact resistance to freshly evaluated BULK-ONLY link conductances.
///
/// Invoke inside the existing backward-Euler or circuit/enthalpy coefficient
/// callback so R'' follows the SAME endpoint enthalpies as k and electrical
/// resistivity. Every selected link must have an explicit bulk coefficient.
/// Unselected links are unchanged. Do not pass a previously contact-corrected
/// map: that would count the same physical contact twice.
///
/// The original series owner composes 1/G_bulk + R''/A. Neither bulk layers nor
/// the interface overwrite one another. Zero area or zero bulk conductance
/// gives exactly zero heat, but still validates both caloric identities and the
/// interface source. R'' itself remains finite and positive; an absent claim
/// never becomes perfect contact. All selected overrides are staged before any
/// map update, so refusal or cancellation leaves the caller's map unchanged.
///
/// Fixed geometry, uniform cells, negligible contact thermal storage. Supplied
/// pressure/history/gap conditions remain explicit, not mechanically evolved.
/// Curves and supported two-coordinate surfaces may respond to the declared
/// evolving thermal/phase axes and typed fixed external conditions. History,
/// field and geometry response axes require their own evolving owners.
/// This is not a contact-area, adhesion, spreading-resistance or phase-mixture
/// solver. Its nonlinear and temporal accuracy belong to the calling integrator.
pub fn apply_material_contacts(
    thermal: &DepositionSimulation<'_>, contacts: &[MaterialContactLink<'_>],
    bulk_conductances: &mut BTreeMap<u64, f64>, cancelled: impl FnMut() -> bool,
) -> Result<Vec<ContactLinkSample>, ConductionError> {
    apply_material_contacts_with_pressures(thermal, contacts, &[], bulk_conductances, cancelled)
}

/// The same original bulk/contact series composition at mechanically supplied
/// interface pressures. Unselected contacts retain their authored conditions.
/// Every pressure names a unique existing contact and a distinct original
/// static-pressure coordinate; other state/history coordinates remain intact.
///
/// Both thermal states, source pins, reference areas and bulk conductances keep
/// their original owners. This does not infer contact mechanics or transform
/// thermal geometry. All source/series evaluation finishes before publication.
pub fn apply_material_contacts_with_pressures(
    thermal: &DepositionSimulation<'_>, contacts: &[MaterialContactLink<'_>],
    pressures: &[ThermalContactPressure<'_>],
    bulk_conductances: &mut BTreeMap<u64, f64>, cancelled: impl FnMut() -> bool,
) -> Result<Vec<ContactLinkSample>, ConductionError> {
    apply_material_contacts_with_metrics(thermal, contacts, pressures, &[], bulk_conductances, cancelled)
}

/// The same series law on current mechanically owned thermal paths.
/// Bulk conductances MUST already use these same metrics and fresh original
/// conductivities. This owner applies R''/A once, with the identical current A.
/// Noncontact path metrics are admitted but receive no interface correction.
/// No area is inferred from electrical resistance, and no cell mass is changed.
pub fn apply_material_contacts_with_metrics(
    thermal: &DepositionSimulation<'_>, contacts: &[MaterialContactLink<'_>],
    pressures: &[ThermalContactPressure<'_>], metrics: &[ThermalLinkMetric],
    bulk_conductances: &mut BTreeMap<u64, f64>, mut cancelled: impl FnMut() -> bool,
) -> Result<Vec<ContactLinkSample>, ConductionError> {
    check(&mut cancelled)?;
    if contacts.len() > MAX_MATERIAL_NETWORK_CONTACTS || contacts.len() > thermal.links().len()
        || pressures.len() > contacts.len() || metrics.len() > thermal.links().len() {
        return Err(invalid("complete bounded active contact bindings required"));
    }
    let mut ordered: Vec<_> = contacts.iter().collect();
    ordered.sort_by_key(|c| c.link);
    let mut pressure_by_link = BTreeMap::new();
    for pressure in pressures {
        check(&mut cancelled)?;
        if ordered.binary_search_by_key(&pressure.link, |contact| contact.link).is_err()
            || pressure_by_link.insert(pressure.link, pressure).is_some() {
            return Err(invalid("mechanical pressure must name each selected original thermal contact once"));
        }
    }
    let mut metric_by_link = BTreeMap::new();
    for metric in metrics {
        check(&mut cancelled)?;
        metric.validate()?;
        if !thermal.links().contains_key(&metric.link)
            || metric_by_link.insert(metric.link, metric).is_some() {
            return Err(invalid("current thermal metrics must name unique original active links"));
        }
    }
    let mut samples = Vec::with_capacity(ordered.len());
    for (index, binding) in ordered.iter().enumerate() {
        check(&mut cancelled)?;
        if index > 0 && ordered[index-1].link == binding.link {
            return Err(invalid("a physical link cannot receive duplicate contact resistance"));
        }
        let mut source = resolve_contact_at_pressure(
            thermal, binding, pressure_by_link.get(&binding.link).copied())?;
        if let Some(metric) = metric_by_link.get(&binding.link) {
            source.area_m2 = metric.area_m2;
        }
        let bulk = *bulk_conductances.get(&binding.link).ok_or_else(|| invalid("contact lacks explicit bulk conductance"))?;
        if !bulk.is_finite() || bulk < 0.0 || !source.area_m2.is_finite() || source.area_m2 < 0.0 {
            return Err(invalid("finite nonnegative bulk conductance and contact area required"));
        }
        let (conductance, series) = if bulk == 0.0 || source.area_m2 == 0.0 { (0.0, None) } else {
            let bulk_term = ThermalResistanceTerm::declared("bulk", 1.0/bulk,
                ResistanceUncertainty::Unstated, "original bulk-only endpoint coefficient; source retained by caller")?;
            let series = SeriesThermalResistance::new(vec![bulk_term, source.interface.term_for_area("contact", source.area_m2)?])?;
            let g = 1.0/series.budget().value_k_per_w;
            if !g.is_finite() || g <= 0.0 { return Err(invalid("series contact conductance is unrepresentable")); }
            (g, Some(series))
        };
        samples.push(ContactLinkSample { link:binding.link, cells:source.cells, phases:source.phases, interface:source.interface,
            area_m2:source.area_m2, bulk_conductance_w_k:bulk, conductance_w_k:conductance, series });
    }
    check(&mut cancelled)?;
    for sample in &samples {
        // Every key was checked above. No fallible work follows publication.
        *bulk_conductances.get_mut(&sample.link).expect("validated coefficient") = sample.conductance_w_k;
    }
    Ok(samples)
}

// Keep initial, callback and endpoint source evaluation under one owner.
// Endpoint admission must not invent a bulk coefficient merely to query R''.
fn resolve_contact(
    thermal: &DepositionSimulation<'_>, binding: &MaterialContactLink<'_>,
) -> Result<ContactInterfaceState, ConductionError> {
    resolve_contact_at_pressure(thermal, binding, None)
}

fn resolve_contact_at_pressure(
    thermal: &DepositionSimulation<'_>, binding: &MaterialContactLink<'_>,
    pressure: Option<&ThermalContactPressure<'_>>,
) -> Result<ContactInterfaceState, ConductionError> {
    let link = thermal.links().get(&binding.link).ok_or_else(|| invalid("inactive thermal contact"))?;
    if [binding.sides[0].cell, binding.sides[1].cell] != [link.a, link.b] {
        return Err(invalid("caloric sides must retain the original ordered link endpoints"));
    }
    let phase = |side: ContactCell<'_>| {
        let cell = thermal.cells().get(&side.cell).ok_or_else(|| invalid("inactive contact cell"))?;
        if cell.material_id != side.material_id
            || !thermal.uses_material_law(side.material_id, side.curve as &dyn ThermalLaw) {
            return Err(invalid("contact must use each cell's actual borrowed caloric owner"));
        }
        side.curve.state_at_specific_enthalpy(cell.specific_h_j_kg)
            .map_err(|e| invalid(&e.to_string()))
    };
    let phases = [phase(binding.sides[0])?, phase(binding.sides[1])?];
    let point = if let Some(pressure) = pressure {
        binding.point.at_pressure(phases, pressure.pressure_axis, pressure.pressure_pa)
    } else {
        binding.point.at(phases)
    }.map_err(|e| invalid(&e.to_string()))?;
    let name = format!("thermal-link-{}", binding.link);
    let interface = match binding.selection {
        ClaimSelection::Policy(policy) => InterfaceResistance::from_card(&name, binding.point.interface, &point, policy)?,
        ClaimSelection::Pinned(id) => InterfaceResistance::from_card_pinned(&name, binding.point.interface, &point, id)?,
    };
    let claim = binding.point.interface.claims().claim(interface.receipt().selected)
        .ok_or_else(|| invalid("selected contact source is absent"))?;
    if claim.key.quantity() != QuantitySpec::dimensional(AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS) {
        return Err(invalid("static dimensional contact resistance is required"));
    }
    if !claim.value.axes().all(|(axis, dims)| binding.point.is_response_axis(axis, dims)) {
        return Err(invalid("contact response requires declared thermal/phase or typed fixed external coordinates"));
    }
    Ok(ContactInterfaceState { link: binding.link, cells: [link.a, link.b], phases,
        interface, area_m2: link.area_m2 })
}
fn invalid(what: &str) -> ConductionError {
    ConductionError::Interface { interface:"material-network-contact".into(), what:what.into(),
        fix:"supply original ordered states, supported contact data and fresh bulk-only conductances".into() }
}
fn check(cancelled: &mut impl FnMut() -> bool) -> Result<(), ConductionError> {
    if cancelled() { Err(ConductionError::Cancelled { stage:"material-network-contact", at:0 }) } else { Ok(()) }
}
