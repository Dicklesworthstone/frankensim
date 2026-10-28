//! Prepared edge-cochain projection onto the complete coexact space.
//!
//! For selected edge columns D of d1 and the original positive barycentric
//! edge metric M, solve P w = D x, P = D M^-1 D^T, then x_c=M^-1 D^T w.
//! P is POSITIVE SEMIDEFINITE, not globally SPD: its face nullspace is harmless
//! only because D x is compatible. The caller owns a checked solver, budgets,
//! and true projection residual acceptance. No face-potential uniqueness is
//! claimed. The recovered x_c is M-orthogonal to ALL ker(D), including harmonic
//! cochains; no guessed Betti number or incomplete random basis is required.

use super::{ElementGeometry, TetComplex, hodge_diagonal_barycentric};
use crate::assembly::stiffness;
use fs_sparse::{Csr, ops};

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CoexactProjectionError {
    Invalid(&'static str),
    Cancelled,
}
impl core::fmt::Display for CoexactProjectionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "coexact projection: {self:?}")
    }
}
impl std::error::Error for CoexactProjectionError {}
type Result<T> = core::result::Result<T, CoexactProjectionError>;
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(CoexactProjectionError::Cancelled)
    } else {
        Ok(())
    }
}

/// Original discrete-curl columns and metric, with zero face rows removed.
/// Geometry must be the original checked geometry of this tetrahedral complex.
/// Preparation uses the existing Hodge and sparse-composition owners.
pub struct PreparedCoexact1 {
    derivative: Csr,
    transpose: Csr,
    matrix: Csr,
    metric: Vec<f64>,
    active_faces: Vec<usize>,
}
impl PreparedCoexact1 {
    pub fn new(
        mesh: &TetComplex,
        positions: &[[f64; 3]],
        geometry: &ElementGeometry,
        free_edges: &[usize],
        maximum_product_terms: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        use CoexactProjectionError::Invalid;
        poll(cancelled)?;
        if positions.len() != mesh.vertex_count
            || geometry.vol_signed.len() != mesh.tets.len()
            || free_edges.windows(2).any(|p| p[0] >= p[1])
            || free_edges.iter().any(|&e| e >= mesh.edges.len())
        {
            return Err(Invalid(
                "original geometry and strictly ordered unique edge selection required",
            ));
        }
        let full_metric = hodge_diagonal_barycentric(mesh, positions, geometry, 1);
        let mut metric = Vec::with_capacity(free_edges.len());
        let mut index = vec![usize::MAX; mesh.edges.len()];
        for (i, &edge) in free_edges.iter().enumerate() {
            poll(cancelled)?;
            let value = full_metric[edge];
            if !value.is_finite() || value <= 0.0 || !value.recip().is_finite() {
                return Err(Invalid("positive representable coexact metric required"));
            }
            metric.push(value);
            index[edge] = i;
        }
        let incidence = mesh.d1();
        let mut row_ptr = vec![0];
        let mut columns = Vec::new();
        let mut values = Vec::new();
        let mut counts = vec![0usize; free_edges.len()];
        let mut active_faces = Vec::new();
        for (face, row) in incidence.rows.iter().enumerate() {
            poll(cancelled)?;
            let mut selected: Vec<_> = row
                .iter()
                .filter_map(|&(edge, sign)| {
                    (index[edge] != usize::MAX).then_some((index[edge], f64::from(sign)))
                })
                .collect();
            if selected.is_empty() {
                continue;
            }
            selected.sort_unstable_by_key(|&(edge, _)| edge);
            for (edge, value) in selected {
                columns.push(edge);
                values.push(value);
                counts[edge] += 1;
            }
            row_ptr.push(columns.len());
            active_faces.push(face);
        }
        let mut terms = 0usize;
        for count in counts {
            poll(cancelled)?;
            terms = count
                .checked_mul(count)
                .and_then(|n| terms.checked_add(n))
                .ok_or(Invalid("coexact product work overflow"))?;
            if terms > maximum_product_terms {
                return Err(Invalid("coexact sparse-product work budget exceeded"));
            }
        }
        let derivative = Csr::from_parts(
            active_faces.len(),
            free_edges.len(),
            row_ptr,
            columns,
            values,
        );
        let transpose = ops::transpose(&derivative);
        let inverse = Csr::from_parts(
            metric.len(),
            metric.len(),
            (0..=metric.len()).collect(),
            (0..metric.len()).collect(),
            metric.iter().map(|v| v.recip()).collect(),
        );
        poll(cancelled)?;
        let matrix = stiffness(&transpose, &inverse);
        for r in 0..matrix.nrows() {
            poll(cancelled)?;
            if matrix.row(r).1.iter().any(|v| !v.is_finite()) {
                return Err(Invalid("unrepresentable coexact operator"));
            }
        }
        Ok(Self {
            derivative,
            transpose,
            matrix,
            metric,
            active_faces,
        })
    }
    pub const fn matrix(&self) -> &Csr {
        &self.matrix
    }
    pub const fn derivative(&self) -> &Csr {
        &self.derivative
    }
    pub fn metric(&self) -> &[f64] {
        &self.metric
    }
    pub fn active_faces(&self) -> &[usize] {
        &self.active_faces
    }
    /// Compatible RHS D x; face coefficients remain an auxiliary gauge solve.
    pub fn rhs(&self, x: &[f64], cancelled: &mut impl FnMut() -> bool) -> Result<Vec<f64>> {
        apply(&self.derivative, x, cancelled)
    }
    /// Dual cochain D^T w, orthogonal to the entire curl kernel.
    pub fn recover_dual(
        &self,
        w: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<f64>> {
        apply(&self.transpose, w, cancelled)
    }
    /// Primal coexact cochain M^-1 D^T w; the caller checks projection residuals.
    pub fn recover(&self, w: &[f64], cancelled: &mut impl FnMut() -> bool) -> Result<Vec<f64>> {
        let mut x = self.recover_dual(w, cancelled)?;
        for (value, weight) in x.iter_mut().zip(&self.metric) {
            poll(cancelled)?;
            *value /= weight;
            if !value.is_finite() {
                return Err(CoexactProjectionError::Invalid(
                    "unrepresentable coexact recovery",
                ));
            }
        }
        Ok(x)
    }
}
fn apply(matrix: &Csr, x: &[f64], cancelled: &mut impl FnMut() -> bool) -> Result<Vec<f64>> {
    use CoexactProjectionError::Invalid;
    if x.len() != matrix.ncols() || x.iter().any(|v| !v.is_finite()) {
        return Err(Invalid("complete finite cochain required"));
    }
    let mut result = Vec::with_capacity(matrix.nrows());
    for r in 0..matrix.nrows() {
        poll(cancelled)?;
        let (cols, vals) = matrix.row(r);
        let value = cols.iter().zip(vals).map(|(&c, &v)| v * x[c]).sum::<f64>();
        if !value.is_finite() {
            return Err(Invalid("unrepresentable coexact application"));
        }
        result.push(value);
    }
    Ok(result)
}
