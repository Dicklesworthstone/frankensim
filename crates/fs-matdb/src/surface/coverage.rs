//! Continuous rectangle support in an already admitted source complex.

use super::{BTreeMap, MatDbError, Sign, contains, coordinate_admitted, orient2d, vertices};

// Separating axes for a segment and an OPEN axis-aligned rectangle are the
// two coordinate axes and the segment normal. Strict projection overlap on
// each coordinate, and rectangle corners on both sides of the supporting
// line, are necessary and sufficient. Corner/edge-only contact is excluded.
// No rounded clipping/intersection point or geometric tolerance is used.
fn enters_interior(a: [f64; 2], b: [f64; 2], lower: [f64; 2], upper: [f64; 2]) -> bool {
    if (0..2).any(|axis| {
        a[axis].max(b[axis]) <= lower[axis] || a[axis].min(b[axis]) >= upper[axis]
    }) {
        return false;
    }
    let mut positive = false;
    let mut negative = false;
    for corner in [lower, [upper[0], lower[1]], upper, [lower[0], upper[1]]] {
        match orient2d(a, b, corner) {
            Sign::Positive => positive = true,
            Sign::Negative => negative = true,
            Sign::Zero => {},
        }
    }
    positive && negative
}

fn interior_coordinate(lower: f64, upper: f64) -> Option<f64> {
    // Source/query bounds keep this subtraction and addition finite. A rounded
    // endpoint is not an interior witness, even for a one-ULP-wide request.
    let middle = lower + (upper - lower) * 0.5;
    if lower < middle && middle < upper && coordinate_admitted(middle) {
        Some(middle)
    } else if lower < 0.0 && 0.0 < upper {
        // Zero is admitted by the original exact-predicate owner; a tiny
        // nonzero rounded midpoint may lie outside its no-underflow domain.
        Some(0.0)
    } else {
        None
    }
}

/// Prove that declared source triangles cover a positive-area closed box.
/// Source admission has already established a conforming nonoverlapping
/// complex: every edge is shared by two opposite cells or is exposed.
/// Its finite triangle union is closed. If no exposed edge enters the
/// connected open rectangle, that rectangle cannot cross the union boundary.
/// One covered interior witness therefore proves its entire interior covered;
/// closure supplies the rectangle boundary. Interior holes and disconnected
/// islands have exposed edges and cannot be filled by this proof.
///
/// Isolated samples never provide the witness. All exact predicates retain
/// their original coordinate domain, and an unresolved interior witness
/// refuses. Source admission bounds work by 1024 cells and 3072 edge records.
/// Degenerate source-coordinate boxes use the existing point/line owners.
pub(crate) fn rectangle_support(
    samples: &[([f64; 2], f64)],
    triangles: &[[u32; 3]],
    lower: [f64; 2],
    upper: [f64; 2],
) -> Result<(), MatDbError> {
    let refusal = |reason| MatDbError::UnsupportedEvaluation { reason };
    if !(0..2).all(|axis| lower[axis] < upper[axis]
        && coordinate_admitted(lower[axis]) && coordinate_admitted(upper[axis]))
    {
        return Err(refusal("surface rectangle requires admitted positive coordinate spans"));
    }
    let [Some(x), Some(y)] = [
        interior_coordinate(lower[0], upper[0]),
        interior_coordinate(lower[1], upper[1]),
    ] else {
        return Err(refusal("surface rectangle interior witness is numerically unresolved"));
    };
    if !triangles.iter().any(|&triangle| contains(vertices(samples, triangle), [x, y])) {
        return Err(refusal("surface rectangle interior is not covered by declared source cells"));
    }

    let mut edges = BTreeMap::<[u32; 2], usize>::new();
    for triangle in triangles {
        for edge in 0..3 {
            let (a, b) = (triangle[edge], triangle[(edge + 1) % 3]);
            *edges.entry([a.min(b), a.max(b)]).or_default() += 1;
        }
    }
    for ([a, b], count) in edges {
        if count == 1 && enters_interior(samples[a as usize].0, samples[b as usize].0, lower, upper) {
            return Err(refusal("an exposed source edge enters the requested surface rectangle"));
        }
    }
    Ok(())
}
