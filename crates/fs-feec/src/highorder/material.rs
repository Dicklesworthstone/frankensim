//! Piecewise-constant materials on the hierarchical tetrahedral H1 space.
//!
//! Mass pairs scalar values; diffusion pairs physical gradients with a scalar
//! or Cartesian SPD tensor. Both reuse the existing oriented basis and exact-
//! degree Duffy rule. No incidence, boundary conditions, or gauge is changed.
//! The matrix carries dual-functional units, not a new primal cochain degree.
//!
//! Orders 1..=6 bound each basis/quadrature tile. Construction admits geometry,
//! coefficients and global extents before numerical work. Legacy basis-list
//! and quadrature helpers make small, order-bounded allocations; large scratch
//! and result allocations are fallible. Checkpoints include quadrature rows,
//! sorting, duplicate reduction and publication. No worker tasks are spawned.
//! Input geometry is checked for finite metrics, not certified against an
//! embedding. Material positivity alone is not a PDE stability certificate.

use super::simplex::{LocalFn, SimplexSpace, duffy_quadrature, entity_dofs};
use crate::weighted::{
    CellWeight, WeightedAssemblyLimits, WeightedError, assemble_triplets, poll, validate_tensor,
};
use crate::whitney::ElementGeometry;
use fs_qty::Dims;
use fs_sparse::Csr;

/// Largest admitted polynomial order; at most 84 local basis functions.
pub const MAX_MATERIAL_ORDER: usize = 6;
const MAX_LOCAL: usize = 84;

/// The physical bulk pairing; diffusion is not a scalar mass matrix.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum H1Pairing {
    /// Integral of coefficient times test value times trial value.
    Mass,
    /// Integral of test gradient dotted with coefficient times trial gradient.
    Diffusion,
}

/// An assembled H1 operator, including dimensions of its dual output.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialMatrix {
    matrix: Csr,
    pairing: H1Pairing,
    coefficient_dims: Dims,
    entry_dims: Dims,
}

impl MaterialMatrix {
    /// Matrix in the existing vertex/edge/face/interior H1 ordering.
    #[must_use]
    pub const fn matrix(&self) -> &Csr {
        &self.matrix
    }
    /// Which physical bilinear form was assembled.
    #[must_use]
    pub const fn pairing(&self) -> H1Pairing {
        self.pairing
    }
    /// Dimensions of the cellwise material coefficient.
    #[must_use]
    pub const fn coefficient_dims(&self) -> Dims {
        self.coefficient_dims
    }
    /// Mass entries multiply the coefficient by m^3; diffusion entries by m.
    #[must_use]
    pub const fn entry_dims(&self) -> Dims {
        self.entry_dims
    }
    /// Dimensions of the dual test-functional result, given primal field units.
    pub fn output_dims(&self, input: Dims) -> Result<Dims, WeightedError> {
        self.entry_dims
            .checked_plus(input)
            .ok_or(WeightedError::DimensionOverflow)
    }
}

/// Admitted, reusable material pairing over an existing H1 space.
///
/// Geometry and coefficients are borrowed immutably. Updating a nonlinear
/// iterate requires a newly admitted operator; this type does not supply a
/// nonlinear law, its derivative, or any history-dependent constitutive state.
pub struct MaterialOperator<'a, 'c> {
    space: &'a SimplexSpace<'c>,
    geometry: &'a ElementGeometry,
    weights: &'a [CellWeight],
    pairing: H1Pairing,
    coefficient_dims: Dims,
    entry_dims: Dims,
    limits: WeightedAssemblyLimits,
    local_count: usize,
}

impl<'a, 'c> MaterialOperator<'a, 'c> {
    /// Admit bounded order, topology/layout, finite metrics and material laws.
    ///
    /// The caller retains ownership of its cancellation scope; for example,
    /// pass `&mut || cx.checkpoint().is_err()`. No numerical work is published
    /// on refusal. `max_triplets` applies to assembly, not to this borrowed plan.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn new(
        space: &'a SimplexSpace<'c>,
        geometry: &'a ElementGeometry,
        pairing: H1Pairing,
        weights: &'a [CellWeight],
        coefficient_dims: Dims,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, WeightedError> {
        poll(cancelled)?;
        let invalid = |cell, reason| WeightedError::InvalidGeometry { cell, reason };
        if !(1..=MAX_MATERIAL_ORDER).contains(&space.r) {
            return Err(invalid(0, "material H1 order must be 1..=6"));
        }
        let complex = space.complex;
        let cells = complex.tets.len();
        limits.check(cells, space.ndof, 0)?;
        let (pe, pf, pc) = entity_dofs(space.r);
        let face_off = complex
            .edges
            .len()
            .checked_mul(pe)
            .and_then(|n| n.checked_add(complex.vertex_count))
            .ok_or(WeightedError::SizeOverflow)?;
        let cell_off = complex
            .faces
            .len()
            .checked_mul(pf)
            .and_then(|n| n.checked_add(face_off))
            .ok_or(WeightedError::SizeOverflow)?;
        let ndof = cells
            .checked_mul(pc)
            .and_then(|n| n.checked_add(cell_off))
            .ok_or(WeightedError::SizeOverflow)?;
        if (
            space.per_edge,
            space.per_face,
            space.per_cell,
            space.edge_off,
            space.face_off,
            space.cell_off,
            space.ndof,
        ) != (pe, pf, pc, complex.vertex_count, face_off, cell_off, ndof)
        {
            return Err(invalid(0, "H1 layout does not match its order and complex"));
        }
        if weights.len() != cells {
            return Err(WeightedError::CoefficientCount {
                expected: cells,
                actual: weights.len(),
            });
        }
        if geometry.vol_signed.len() != cells || geometry.grads.len() != cells {
            return Err(invalid(0, "geometry arrays must match tetrahedra"));
        }
        let length_power = match pairing {
            H1Pairing::Mass => 3,
            H1Pairing::Diffusion => 1,
        };
        let entry_dims = coefficient_dims
            .checked_plus(Dims([length_power, 0, 0, 0, 0, 0]))
            .ok_or(WeightedError::DimensionOverflow)?;
        for (i, edge) in complex.edges.iter().enumerate() {
            poll(cancelled)?;
            if edge[0] >= edge[1]
                || edge[1] as usize >= complex.vertex_count
                || (i > 0 && complex.edges[i - 1] >= *edge)
            {
                return Err(invalid(0, "noncanonical edge table"));
            }
        }
        for (i, face) in complex.faces.iter().enumerate() {
            poll(cancelled)?;
            if face[0] >= face[1]
                || face[1] >= face[2]
                || face[2] as usize >= complex.vertex_count
                || (i > 0 && complex.faces[i - 1] >= *face)
            {
                return Err(invalid(0, "noncanonical face table"));
            }
        }
        for (cell, (&tet, weight)) in complex.tets.iter().zip(weights).enumerate() {
            poll(cancelled)?;
            match weight {
                CellWeight::Scalar(value) if value.is_finite() && *value > 0.0 => {}
                CellWeight::Scalar(_) => {
                    return Err(WeightedError::InvalidCoefficient {
                        cell,
                        reason: "scalar must be finite and strictly positive",
                    });
                }
                CellWeight::Tensor(tensor) if pairing == H1Pairing::Diffusion => {
                    validate_tensor(tensor, cell)?
                }
                CellWeight::Tensor(_) => {
                    return Err(WeightedError::InvalidCoefficient {
                        cell,
                        reason: "H1 mass requires a scalar coefficient",
                    });
                }
            }
            let volume = geometry.vol_signed[cell].abs();
            if volume == 0.0
                || !(6.0 * volume).is_finite()
                || geometry.grads[cell]
                    .iter()
                    .flatten()
                    .any(|v| !v.is_finite())
            {
                return Err(invalid(cell, "nonfinite or degenerate element metrics"));
            }
            for i in 0..4 {
                if tet[i] as usize >= complex.vertex_count || tet[..i].contains(&tet[i]) {
                    return Err(invalid(cell, "missing or repeated tetrahedron vertex"));
                }
                for j in i + 1..4 {
                    let key = [tet[i].min(tet[j]), tet[i].max(tet[j])];
                    if complex.edges.binary_search(&key).is_err() {
                        return Err(invalid(cell, "tetrahedron edge is absent"));
                    }
                }
                if pf > 0 {
                    let mut face = [0; 3];
                    let mut next = 0;
                    for (j, &vertex) in tet.iter().enumerate() {
                        if i != j {
                            face[next] = vertex;
                            next += 1;
                        }
                    }
                    face.sort_unstable();
                    if complex.faces.binary_search(&face).is_err() {
                        return Err(invalid(cell, "tetrahedron face is absent"));
                    }
                }
            }
        }
        poll(cancelled)?;
        Ok(Self {
            space,
            geometry,
            weights,
            pairing,
            coefficient_dims,
            entry_dims,
            limits,
            local_count: 4 + 6 * pe + 4 * pf + pc,
        })
    }

    /// Assemble the material-weighted matrix with deterministic cell reduction.
    ///
    /// One local upper triangle is evaluated and mirrored, making symmetry
    /// bit-exact. Coefficient discontinuities are integrated before shared
    /// degrees of freedom are summed. A diffusion matrix may have a constant
    /// nullspace; positive material does not make an unconstrained solve SPD.
    pub fn assemble(
        &self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MaterialMatrix, WeightedError> {
        poll(cancelled)?;
        let nl = self.local_count;
        let triplets = self
            .space
            .complex
            .tets
            .len()
            .checked_mul(nl * nl)
            .ok_or(WeightedError::SizeOverflow)?;
        self.limits
            .check(self.space.complex.tets.len(), self.space.ndof, triplets)?;
        let mut entries = Vec::new();
        entries
            .try_reserve_exact(triplets)
            .map_err(|_| WeightedError::Allocation)?;
        let mut local = zeros(nl * nl)?;
        let quad = duffy_quadrature(self.space.r + 2);
        let mut values = [[0.0; 3]; MAX_LOCAL];
        for (cell, &weight) in self.weights.iter().enumerate() {
            poll(cancelled)?;
            local.fill(0.0);
            let dofs = self.space.element_dofs(cell);
            let (scale, tensor) = normalized(weight);
            for &(lam, w) in &quad {
                self.sample(cell, lam, &dofs, &mut values, cancelled)?;
                let factor = w * (6.0 * self.geometry.vol_signed[cell].abs());
                for i in 0..nl {
                    poll(cancelled)?;
                    for j in i..nl {
                        let value =
                            factor.mul_add(pair(tensor, values[i], values[j]), local[i * nl + j]);
                        if !value.is_finite() {
                            return Err(WeightedError::NonFiniteEntry);
                        }
                        local[i * nl + j] = value;
                    }
                }
            }
            for (i, &(gi, _)) in dofs.iter().enumerate() {
                poll(cancelled)?;
                for (j, &(gj, _)) in dofs.iter().enumerate().skip(i) {
                    let value = scale * local[i * nl + j];
                    if !value.is_finite() {
                        return Err(WeightedError::NonFiniteEntry);
                    }
                    if i == j && value <= 0.0 {
                        return Err(WeightedError::InvalidGeometry {
                            cell,
                            reason: "positive material diagonal is not representable",
                        });
                    }
                    entries.push((gi, gj, value));
                    if i != j {
                        entries.push((gj, gi, value));
                    }
                }
            }
        }
        let matrix = assemble_triplets(self.space.ndof, entries, cancelled)?;
        Ok(MaterialMatrix {
            matrix,
            pairing: self.pairing,
            coefficient_dims: self.coefficient_dims,
            entry_dims: self.entry_dims,
        })
    }

    fn sample(
        &self,
        cell: usize,
        lam: [f64; 4],
        dofs: &[(usize, LocalFn)],
        values: &mut [[f64; 3]; MAX_LOCAL],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), WeightedError> {
        for (i, (_, basis)) in dofs.iter().enumerate() {
            poll(cancelled)?;
            values[i] = match self.pairing {
                H1Pairing::Mass => [basis.eval(lam, self.space.r), 0.0, 0.0],
                H1Pairing::Diffusion => {
                    let mut gradient = [0.0; 3];
                    for (a, derivative) in basis.d_lambda(lam, self.space.r).into_iter().enumerate()
                    {
                        for (c, value) in gradient.iter_mut().enumerate() {
                            *value = derivative.mul_add(self.geometry.grads[cell][a][c], *value);
                        }
                    }
                    gradient
                }
            };
            if values[i].iter().any(|v| !v.is_finite()) {
                return Err(WeightedError::NonFiniteEntry);
            }
        }
        Ok(())
    }
}

fn zeros(n: usize) -> Result<Vec<f64>, WeightedError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(n)
        .map_err(|_| WeightedError::Allocation)?;
    values.resize(n, 0.0);
    Ok(values)
}

fn normalized(weight: CellWeight) -> (f64, [[f64; 3]; 3]) {
    match weight {
        CellWeight::Scalar(value) => (value, [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        CellWeight::Tensor(tensor) => {
            let scale = tensor.iter().flatten().fold(0.0_f64, |s, v| s.max(v.abs()));
            (scale, tensor.map(|row| row.map(|v| v / scale)))
        }
    }
}

fn pair(tensor: [[f64; 3]; 3], left: [f64; 3], right: [f64; 3]) -> f64 {
    let mut result = 0.0;
    for i in 0..3 {
        let row = tensor[i][0].mul_add(
            right[0],
            tensor[i][1].mul_add(right[1], tensor[i][2] * right[2]),
        );
        result = left[i].mul_add(row, result);
    }
    result
}
