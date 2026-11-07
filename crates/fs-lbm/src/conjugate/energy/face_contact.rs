//! Face-local thermal resistance in the existing conservative energy equation.

use std::collections::{BTreeMap, BTreeSet};
use fs_exec::CancelGate;
use super::{ChtError, FlowField, ThermalSetup, Voxel, VoxelDomain};

/// One zero-thickness thermal joint on a Cartesian cell face.
///
/// Additional to the bulk half cells: the two surface temperatures differ by
/// R times the heat flux. Constant bulks give G = area/(wa/(2*ka)+R+wb/(2*kb));
/// native k(T) curves use their integrated half-cell laws. A source owner may
/// retain ordered surfaces and resolve R at both local temperatures. Same-material
/// joints, solid/fluid walls and stagnant fluid pairs are allowed. No contact-area,
/// pressure, phase or interfacial storage law is inferred.
/// A compact component's case/board contact is additional to its declared
/// junction resistor; the collapsed package contributes no fictitious bulk leg.
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
    compact_owner: &[usize],
    gate: &CancelGate,
) -> Result<BTreeMap<(usize, usize), f64>, ChtError> {
    prepare_inner(domain, Some(flow), setup, compact_owner, gate)
}

/// The conduction-only owner has no through-flow or collapsed package cells.
pub(crate) fn prepare_conduction(
    domain: &VoxelDomain,
    setup: &ThermalSetup,
    gate: &CancelGate,
) -> Result<BTreeMap<(usize, usize), f64>, ChtError> {
    prepare_inner(domain, None, setup, &[], gate)
}

fn prepare_inner(
    domain: &VoxelDomain,
    flow: Option<&FlowField>,
    setup: &ThermalSetup,
    compact_owner: &[usize],
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
        compact_face(domain, setup, compact_owner, contact.cells)?;
        if flow.is_some_and(|flow| flow.outward(domain, a, face) != 0.0) {
            return Err(refuse(format!("contact {index}: stationary thermal contacts cannot carry through-flow")));
        }
    }
    super::super::poll(gate)?;
    Ok(field)
}

/// The actual package branch behind one case/board surface. Only a surface
/// contact is admitted: collapsed internal faces and adiabatic package sides
/// remain closed, and two compact components cannot share a physical link.
#[derive(Clone, Copy)]
pub(super) struct CompactFace {
    pub alias: usize,
    pub outside: usize,
    pub axis: usize,
    pub package_first: bool,
    pub package_resistance_m2_k_w: f64,
}

pub(super) fn compact_face(
    domain: &VoxelDomain, setup: &ThermalSetup, owner: &[usize], cells: [usize; 2],
) -> Result<Option<CompactFace>, ChtError> {
    let invalid = |reason: &str| ChtError::InvalidInput {
        field: "thermal.face_contacts", reason: reason.into(),
    };
    let [a, b] = cells;
    if a >= domain.cell_count() || b >= domain.cell_count() || a == b {
        return Err(invalid("contact needs two distinct domain cells"));
    }
    let face = (0..6).find(|&f| domain.neighbor(a, f) == Some(b))
        .ok_or_else(|| invalid("contact cells are not face-adjacent"))?;
    let component = |c: usize| owner.get(c).copied().filter(|&i| i != usize::MAX);
    let (index, alias, outside, outward, package_first) = match (component(a), component(b)) {
        (None, None) => return Ok(None),
        (Some(i), None) => (i, a, b, face, true),
        (None, Some(i)) => (i, b, a, face ^ 1, false),
        (Some(_), Some(_)) => return Err(invalid(
            "compact contacts require one case/board surface and one ordinary outside cell",
        )),
    };
    let part = &setup.compact_components[index];
    let board = part.board_face as usize;
    if outward / 2 != board / 2 {
        return Err(invalid("compact sides are adiabatic; contacts require a case or board face"));
    }
    let resistance = if outward == board {
        part.junction_to_board_k_w
    } else {
        part.junction_to_case_k_w
    } * part.face_area(domain);
    if !resistance.is_finite() || resistance <= 0.0 {
        return Err(invalid("compact branch area resistance is unrepresentable"));
    }
    Ok(Some(CompactFace { alias, outside, axis: outward / 2,
        package_first, package_resistance_m2_k_w: resistance }))
}
