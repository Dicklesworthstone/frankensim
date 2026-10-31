//! Face-local thermal resistance in the existing conservative energy equation.

use std::collections::{BTreeMap, BTreeSet};
use fs_exec::CancelGate;
use super::{ChtError, Compact, FlowField, ThermalSetup, Voxel, VoxelDomain};

/// One zero-thickness thermal joint on a Cartesian cell face.
///
/// Additional to the bulk half cells: the two surface temperatures differ by
/// R times the heat flux. Constant bulks give G = area/(wa/(2*ka)+R+wb/(2*kb));
/// native k(T) curves use their integrated half-cell laws. A source owner may
/// retain ordered surfaces and resolve R at both local temperatures. Same-material
/// joints, solid/fluid walls and stagnant fluid pairs are allowed. No contact-area,
/// pressure, phase or interfacial storage law is inferred.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FaceContactResistance {
    /// Two distinct face-adjacent cells, not material indices.
    pub cells: [usize; 2],
    /// Additional area-specific resistance in square metres kelvin per watt.
    /// Must be finite and non-negative. Explicit zero denotes perfect contact.
    pub resistance_m2_k_w: f64,
}

pub(super) fn prepare(
    domain: &VoxelDomain,
    flow: &FlowField,
    setup: &ThermalSetup,
    compact: &Compact,
    gate: &CancelGate,
) -> Result<BTreeMap<(usize, usize), f64>, ChtError> {
    super::super::poll(gate)?;
    let refuse = |reason: String| ChtError::InvalidInput {
        field: "thermal.face_contacts", reason,
    };
    if setup.face_contacts.len() > domain.cell_count().saturating_mul(3) {
        return Err(refuse("contact count exceeds Cartesian face count".into()));
    }
    let pairs: BTreeSet<_> = setup.contacts.iter().map(|r| {
        let (a, b) = r.materials;
        (a.min(b), a.max(b))
    }).collect();
    let mut field = BTreeMap::new();
    for (index, contact) in setup.face_contacts.iter().enumerate() {
        if index.is_multiple_of(256) { super::super::poll(gate)?; }
        let [a, b] = contact.cells;
        if a >= domain.cell_count() || b >= domain.cell_count() || a == b {
            return Err(refuse(format!("contact {index}: needs two distinct domain cells")));
        }
        let face = (0..6).find(|&f| domain.neighbor(a, f) == Some(b))
            .ok_or_else(|| refuse(format!("contact {index}: cells are not face-adjacent")))?;
        let key = (a.min(b), a.max(b));
        let r = contact.resistance_m2_k_w;
        if !r.is_finite() || r < 0.0 {
            return Err(refuse(format!("contact {index}: resistance must be finite and non-negative")));
        }
        if field.insert(key, r).is_some() {
            return Err(refuse(format!("contact {index}: duplicate face, including reversed order")));
        }
        if let (Voxel::Solid(ma), Voxel::Solid(mb)) = (domain.voxel_at(a), domain.voxel_at(b))
            && pairs.contains(&(ma.min(mb), ma.max(mb)))
        {
            return Err(refuse(format!("contact {index}: face also has a material-pair contact law")));
        }
        if compact.owner(a).is_some() || compact.owner(b).is_some() {
            return Err(refuse(format!("contact {index}: compact links replace ordinary face transport")));
        }
        if flow.outward(domain, a, face) != 0.0 {
            return Err(refuse(format!("contact {index}: stationary thermal contacts cannot carry through-flow")));
        }
    }
    super::super::poll(gate)?;
    Ok(field)
}
