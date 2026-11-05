//! Admission of the original, publicly mutable vector-space representation.

use super::super::{
    Chart, Family, MAX_R, VecSpace, cross3, dot3, nedelec_entity_dofs, rt_entity_dofs,
};
use super::{MAX_LOCAL, MAX_MONOMIAL};
use crate::weighted::{
    CellWeight, WeightedAssemblyLimits, WeightedError, poll, validate_affine_scalar,
    validate_tensor,
};

pub(super) struct Admitted {
    pub(super) dofs: Vec<usize>,
    pub(super) vol6: Vec<f64>,
    pub(super) local_count: usize,
}

fn invalid(cell: usize, reason: &'static str) -> WeightedError {
    WeightedError::InvalidGeometry { cell, reason }
}

pub(super) fn admit(
    space: &VecSpace<'_>,
    positions: &[[f64; 3]],
    weights: &[CellWeight],
    limits: WeightedAssemblyLimits,
    max_storage_bytes: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Admitted, WeightedError> {
    poll(cancelled)?;
    if !(1..=MAX_R).contains(&space.r) {
        return Err(invalid(
            0,
            "vector-family order is outside the supported range",
        ));
    }
    let complex = space.complex;
    let cells = complex.tets.len();
    // Preparation does not stage triplets. Assembly admits that extent itself.
    limits.check(cells, space.ndof, 0)?;
    if positions.len() != complex.vertex_count || space.elements.len() != cells {
        return Err(invalid(
            0,
            "position or element-basis extent disagrees with the complex",
        ));
    }
    if weights.len() != cells {
        return Err(WeightedError::CoefficientCount {
            expected: cells,
            actual: weights.len(),
        });
    }
    let (pe, pf, pc) = match space.family {
        Family::Nedelec => nedelec_entity_dofs(space.r),
        Family::Rt => {
            let (pf, pc) = rt_entity_dofs(space.r);
            (0, pf, pc)
        }
    };
    let face_off = complex
        .edges
        .len()
        .checked_mul(pe)
        .ok_or(WeightedError::SizeOverflow)?;
    let cell_off = complex
        .faces
        .len()
        .checked_mul(pf)
        .and_then(|count| face_off.checked_add(count))
        .ok_or(WeightedError::SizeOverflow)?;
    let ndof = cells
        .checked_mul(pc)
        .and_then(|count| cell_off.checked_add(count))
        .ok_or(WeightedError::SizeOverflow)?;
    if (
        space.per_edge,
        space.per_face,
        space.per_cell,
        space.face_off,
        space.cell_off,
        space.ndof,
    ) != (pe, pf, pc, face_off, cell_off, ndof)
    {
        return Err(invalid(
            0,
            "vector-family entity counts or offsets disagree with the original layout",
        ));
    }
    // r <= MAX_R makes these fixed local extents small before any allocation.
    let local_count = 6 * pe + 4 * pf + pc;
    let monomial_count = (space.r + 1) * (space.r + 2) * (space.r + 3) / 6;
    if local_count > MAX_LOCAL || monomial_count > MAX_MONOMIAL {
        return Err(invalid(
            0,
            "vector-family local extent exceeds the prepared kernel",
        ));
    }
    let flat_count = cells
        .checked_mul(local_count)
        .ok_or(WeightedError::SizeOverflow)?;
    let order = space.r + 3;
    let quadrature_count = order.checked_pow(3).ok_or(WeightedError::SizeOverflow)?;
    let mut bytes = 0usize;
    for (count, width) in [
        (flat_count, core::mem::size_of::<usize>()),
        (cells, core::mem::size_of::<f64>()),
        (ndof, core::mem::size_of::<f64>()),
        (quadrature_count, core::mem::size_of::<([f64; 4], f64)>()),
    ] {
        bytes = count
            .checked_mul(width)
            .and_then(|part| bytes.checked_add(part))
            .ok_or(WeightedError::SizeOverflow)?;
    }
    if bytes > max_storage_bytes {
        return Err(WeightedError::LimitExceeded {
            resource: "prepared vector material mass bytes",
            requested: bytes,
            limit: max_storage_bytes,
        });
    }
    for (i, edge) in complex.edges.iter().enumerate() {
        poll(cancelled)?;
        if edge[0] >= edge[1]
            || edge[1] as usize >= complex.vertex_count
            || (i > 0 && complex.edges[i - 1] >= *edge)
        {
            return Err(invalid(
                0,
                "edge table is not canonical or has an invalid vertex",
            ));
        }
    }
    for (i, face) in complex.faces.iter().enumerate() {
        poll(cancelled)?;
        if face[0] >= face[1]
            || face[1] >= face[2]
            || face[2] as usize >= complex.vertex_count
            || (i > 0 && complex.faces[i - 1] >= *face)
        {
            return Err(invalid(
                0,
                "face table is not canonical or has an invalid vertex",
            ));
        }
    }
    // Match monomials(r) without an infallible temporary allocation.
    let mut monomials = [[0usize; 3]; MAX_MONOMIAL];
    let mut m = 0;
    for total in 0..=space.r {
        for a in (0..=total).rev() {
            for b in (0..=(total - a)).rev() {
                monomials[m] = [a, b, total - a - b];
                m += 1;
            }
        }
    }
    poll(cancelled)?;
    let mut dofs = Vec::new();
    dofs.try_reserve_exact(flat_count)
        .map_err(|_| WeightedError::Allocation)?;
    let mut vol6 = Vec::new();
    vol6.try_reserve_exact(cells)
        .map_err(|_| WeightedError::Allocation)?;
    for (cell, tet) in complex.tets.iter().enumerate() {
        poll(cancelled)?;
        for (i, &vertex) in tet.iter().enumerate() {
            if vertex as usize >= positions.len() || tet[..i].contains(&vertex) {
                return Err(invalid(
                    cell,
                    "tetrahedron has an invalid or repeated vertex",
                ));
            }
        }
        let corners = tet.map(|vertex| positions[vertex as usize]);
        if corners.iter().flatten().any(|x| !x.is_finite()) {
            return Err(invalid(cell, "tetrahedron coordinates are not finite"));
        }
        let element = &space.elements[cell];
        if element.funcs.len() != local_count
            || element.monos.as_slice() != &monomials[..monomial_count]
        {
            return Err(invalid(
                cell,
                "basis functions or monomial table disagree with the original layout",
            ));
        }
        for function in &element.funcs {
            poll(cancelled)?;
            if function.deg != space.r
                || function.comp.iter().any(|component| {
                    component.len() != monomial_count
                        || component.iter().any(|coefficient| !coefficient.is_finite())
                })
            {
                return Err(invalid(
                    cell,
                    "basis polynomial degree, extent, or coefficient is invalid",
                ));
            }
        }
        let original_chart = Chart::of(&corners);
        if !element.chart.h.is_finite()
            || element.chart.h <= 0.0
            || element.chart.x0.iter().any(|x| !x.is_finite())
            || element.chart.h != original_chart.h
            || element.chart.x0 != original_chart.x0
        {
            return Err(invalid(
                cell,
                "basis chart is invalid or disagrees with the supplied corners",
            ));
        }
        // Chart agreement is a consistency check, not a proof of the caller's
        // original Vandermonde duality or basis-to-geometry correspondence.
        let e1 = core::array::from_fn(|axis| corners[1][axis] - corners[0][axis]);
        let e2 = core::array::from_fn(|axis| corners[2][axis] - corners[0][axis]);
        let e3 = core::array::from_fn(|axis| corners[3][axis] - corners[0][axis]);
        let volume6 = dot3(e1, cross3(e2, e3)).abs();
        if !volume6.is_finite() || volume6 <= 0.0 {
            return Err(invalid(
                cell,
                "tetrahedron volume is not finite and strictly positive",
            ));
        }
        match &weights[cell] {
            CellWeight::Scalar(value) if !value.is_finite() || *value <= 0.0 => {
                return Err(WeightedError::InvalidCoefficient {
                    cell,
                    reason: "scalar coefficient must be finite and strictly positive",
                });
            }
            CellWeight::Scalar(_) => {}
            CellWeight::AffineScalar(values) => validate_affine_scalar(values, cell)?,
            CellWeight::Tensor(tensor) => validate_tensor(tensor, cell)?,
        }
        // This is element_dofs' exact native traversal, directly into the
        // fallibly reserved flat buffer. Missing topology returns an error
        // instead of reaching that method's binary-search expect or Vec growth.
        for p in 0..4 {
            for q in (p + 1)..4 {
                let mut edge = [tet[p], tet[q]];
                edge.sort_unstable();
                let index = complex.edges.binary_search(&edge).map_err(|_| {
                    invalid(cell, "tetrahedron edge is missing from the canonical table")
                })?;
                for k in 0..pe {
                    dofs.push(index * pe + k);
                }
            }
        }
        for omit in 0..4 {
            let mut face = [0u32; 3];
            let mut next = 0;
            for (i, &vertex) in tet.iter().enumerate() {
                if i != omit {
                    face[next] = vertex;
                    next += 1;
                }
            }
            face.sort_unstable();
            let index = complex.faces.binary_search(&face).map_err(|_| {
                invalid(cell, "tetrahedron face is missing from the canonical table")
            })?;
            for k in 0..pf {
                dofs.push(face_off + index * pf + k);
            }
        }
        for k in 0..pc {
            dofs.push(cell_off + cell * pc + k);
        }
        vol6.push(volume6);
    }
    poll(cancelled)?;
    Ok(Admitted {
        dofs,
        vol6,
        local_count,
    })
}
