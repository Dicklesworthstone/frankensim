//! Matching imperfect-contact contribution to the same equilibrated majorant.
//!
//! The energy space has independent H1 traces on the two solids and
//! a(u,w) includes integral jump(u)*jump(w)/R. For an equilibrated flux
//! with opposite normal traces, integration by parts contributes
//! integral (q.n - jump(v)/R)*jump(e), e=u-v. Cauchy-Schwarz in the
//! product of the volume, Robin and contact energy spaces therefore adds
//! integral (R*q.n - jump(v))^2/R to the SAME squared energy majorant.
//! It does not assume that a discontinuous temperature is globally H1.
//!
//! Geometry binding below makes each contact one graph edge. The existing
//! conservative forest supplies a single face flux with opposite signs; the
//! affine-source lifting has zero normal trace and leaves this unchanged.
use super::*;

pub(super) fn jump(face: &Face, partner: [usize; 3], field: &[f64]) -> [Iv; 3] {
    std::array::from_fn(|i| Iv::point(field[face.vertices[i]]).sub(Iv::point(field[partner[i]])))
}

/// Link exact duplicated exterior triangles, never merge their vertex values.
/// Every pair is reciprocal with identical R; a mismatched/one-sided contact
/// cannot silently become insulation, perfect contact or a boundary anchor.
pub(super) fn link(
    vertices: &[[f64; 3]], cells: &mut [Cell], mut faces: Vec<Face>,
    keep_going: &mut impl FnMut() -> bool,
) -> Result<Vec<Face>, TetError> {
    let mut lookup = BTreeMap::new();
    for (i, face) in faces.iter().enumerate() {
        poll(keep_going)?;
        lookup.insert(face.vertices, i);
    }
    let mut pairs = Vec::new();
    for (i, face) in faces.iter().enumerate() {
        poll(keep_going)?;
        let Some(BoundaryCondition::Contact { partner, resistance }) = face.condition else { continue; };
        let &j = lookup.get(&partner).ok_or(TetError::Invalid("contact partner is not an extracted face"))?;
        if i == j || face.sides.len() != 1 || faces[j].sides.len() != 1 {
            return Err(TetError::Invalid("contact must pair two distinct exterior triangles"));
        }
        match faces[j].condition {
            Some(BoundaryCondition::Contact { partner: back, resistance: other })
                if back == face.vertices && other == resistance => {}
            _ => return Err(TetError::Invalid("contact declarations must be reciprocal with identical resistance")),
        }
        if dot(face.normal_area, faces[j].normal_area).hi >= 0.0 {
            return Err(TetError::Invalid("contact solids must lie on opposite sides of their common face"));
        }
        let mut aligned = [0; 3];
        for (slot, &vertex) in face.vertices.iter().enumerate() {
            aligned[slot] = *partner.iter().find(|&&v| vertices[v] == vertices[vertex])
                .ok_or(TetError::Unsupported("contact requires exact matching physical triangles"))?;
        }
        // Canonical geometric face order, independent of declaration order,
        // chooses ownership. Each valid pair is visited twice but linked once.
        if i < j { pairs.push((i, j, aligned, resistance)); }
    }
    let mut removed = vec![false; faces.len()];
    for (a, b, aligned, resistance) in pairs {
        poll(keep_going)?;
        let other = faces[b].sides[0];
        faces[a].sides.push(other);
        faces[a].condition = Some(BoundaryCondition::Contact { partner: aligned, resistance });
        cells[other.0].faces[other.1] = a;
        removed[b] = true;
    }
    let mut remap = vec![usize::MAX; faces.len()];
    let mut linked = Vec::with_capacity(faces.len());
    for (old, face) in faces.into_iter().enumerate() {
        poll(keep_going)?;
        if !removed[old] {
            remap[old] = linked.len();
            linked.push(face);
        }
    }
    for cell in cells {
        poll(keep_going)?;
        for f in &mut cell.faces {
            *f = remap[*f];
            if *f == usize::MAX { return Err(TetError::Invalid("unlinked contact face")); }
        }
    }
    poll(keep_going)?;
    Ok(linked)
}
