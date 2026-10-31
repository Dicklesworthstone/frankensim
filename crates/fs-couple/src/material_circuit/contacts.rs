//! Electrical contact drops and heat on the original material circuit branches.
//! The interface owner supplies specific resistance; physical area supplies its
//! geometric reduction. No fictitious bulk material or contact heat split.
use super::*;
use fs_material::electrical_interface::ResolvedElectricalInterface;
use std::collections::BTreeSet;

/// Maximum independently sourced contacts across the bounded circuit graph.
pub const MAX_SERIES_CONTACTS: usize = 128;

/// One zero-storage contact in series with an existing bulk resistor branch.
/// The electrical direction follows the branch; source A/B order describes its
/// two physical surfaces. Their assignment and temperature closure are explicit.
#[derive(Clone, Debug, PartialEq)]
pub struct SeriesElectricalContact {
    id: u64,
    branch: usize,
    source: ResolvedElectricalInterface,
    area_m2: f64,
    resistance_ohm: f64,
}
impl SeriesElectricalContact {
    /// Reduce the original nonnegative specific resistance by positive area.
    /// Explicit source-backed zero resistance is allowed; missing source data is
    /// not. Zero area is NOT a large resistance: circuit opening is unsupported.
    pub fn new(id: u64, branch: usize, source: ResolvedElectricalInterface, area_m2: f64)
        -> Result<Self, MaterialCircuitError>
    {
        if !positive(area_m2) {
            return Err(MaterialCircuitError::Invalid("electrical contact requires positive finite area"));
        }
        let specific = source.resistance_ohm_m2();
        let resistance_ohm = specific/area_m2;
        if !resistance_ohm.is_finite() || resistance_ohm < 0.0 || (specific != 0.0 && resistance_ohm == 0.0) {
            return Err(MaterialCircuitError::Invalid("specific electrical contact resistance/area is unrepresentable"));
        }
        Ok(Self { id, branch, source, area_m2, resistance_ohm })
    }
    /// Stable contact identity, independent of vector order.
    pub const fn id(&self) -> u64 { self.id }
    /// Original graph branch containing this series drop.
    pub const fn branch(&self) -> usize { self.branch }
    /// Exact source, ordered surface states and original query/receipts.
    pub const fn source(&self) -> &ResolvedElectricalInterface { &self.source }
    /// Physical conducting area [m2].
    pub const fn area_m2(&self) -> f64 { self.area_m2 }
    /// Total contact resistance [ohm], not bulk resistivity.
    pub const fn resistance_ohm(&self) -> f64 { self.resistance_ohm }

    pub(super) fn heat(&self, current_a: f64, duration_s: f64) -> Result<ContactHeat, MaterialCircuitError> {
        let voltage_v = current_a*self.resistance_ohm;
        let power_w = voltage_v*current_a;
        let heat_j = duration_s*power_w;
        if ![voltage_v,power_w,heat_j].iter().all(|v| v.is_finite()) || power_w < 0.0
            || (current_a != 0.0 && self.resistance_ohm != 0.0 && (voltage_v == 0.0 || power_w == 0.0))
            || (power_w != 0.0 && heat_j == 0.0) {
            return Err(MaterialCircuitError::Balance("electrical contact drop or Joule heat is unrepresentable"));
        }
        Ok(ContactHeat { contact:self.clone(), current_a, voltage_v, power_w, heat_j })
    }
}

/// Actual contact work, separate from heat generated inside either bulk body.
#[derive(Clone, Debug, PartialEq)]
pub struct ContactHeat {
    /// Exact source and geometry used in this physical interval.
    pub contact: SeriesElectricalContact,
    /// Signed current in the containing graph branch [A].
    pub current_a: f64,
    /// Signed contact voltage drop in that direction [V].
    pub voltage_v: f64,
    /// Nonnegative contact Joule power [W].
    pub power_w: f64,
    /// Whole-step contact heat [J]; a thermal consumer must deposit it ONCE.
    pub heat_j: f64,
}

pub(super) fn admit(contacts: &mut [SeriesElectricalContact]) -> Result<(), MaterialCircuitError> {
    if contacts.len() > MAX_SERIES_CONTACTS {
        return Err(MaterialCircuitError::Invalid("electrical contact work capacity exceeded"));
    }
    contacts.sort_by_key(SeriesElectricalContact::id);
    if contacts.windows(2).any(|p| p[0].id == p[1].id) {
        return Err(MaterialCircuitError::Invalid("duplicate electrical contact identity"));
    }
    Ok(())
}

// The caller has stamped the original bulk resistors first. Never feed a
// previous total resistance through this reduction a second time.
pub(super) fn lower(graph: &mut CircuitGraph, contacts: &[SeriesElectricalContact]) -> Result<(), MaterialCircuitError> {
    let mut seen = BTreeSet::new();
    if contacts.len() > MAX_SERIES_CONTACTS { return Err(MaterialCircuitError::Invalid("too many electrical contacts")); }
    for contact in contacts {
        if !seen.insert(contact.id) { return Err(MaterialCircuitError::Invalid("duplicate electrical contact")); }
        let branch = graph.branches.get_mut(contact.branch)
            .ok_or(MaterialCircuitError::Invalid("electrical contact branch is absent"))?;
        let Branch::Resistor { ohms } = &mut branch.2 else {
            return Err(MaterialCircuitError::Invalid("series contact must bind an existing bulk resistor branch"));
        };
        let total = *ohms+contact.resistance_ohm;
        if !positive(total) || !positive(1.0/total) || (contact.resistance_ohm != 0.0 && total == *ohms) {
            return Err(MaterialCircuitError::Invalid("contact resistance is lost or overflows the series branch"));
        }
        *ohms = total;
    }
    Ok(())
}

impl MaterialCircuit {
    /// Explicitly add/remove/rebind contacts without advancing charge or flux.
    /// Storage material records survive. The original descriptor must admit the
    /// same conserved coordinates with the new algebraic constraints. No contact
    /// capacitance, arc energy, ideal-switch impulse or mechanical work is inferred.
    /// Use this after `new_with_storage` to compose source-bound C/L and contacts.
    pub fn rebind_contacts_candidate(&self, contacts: Vec<SeriesElectricalContact>, inputs: &[f64])
        -> Result<Self, MaterialCircuitError>
    {
        let mut result = Self::new_with_contacts_and_mutual_inductance(self.graph.clone(),
            self.resistors.clone(), contacts, self.mutual.clone(),
            self.fluxes_wb.clone(), self.charges_c.clone(), inputs)?;
        result.storage = self.storage.clone();
        result.magnetic_couplings = self.magnetic_couplings.clone();
        Ok(result)
    }
}
