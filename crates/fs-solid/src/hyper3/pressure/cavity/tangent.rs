//! Sparse geometric pressure stiffness plus one source acoustic rank-one term.
//!
//! The force Jacobian is J = (p-pe) Hess(V) - w w^T, where
//! w = sqrt(rho*c^2/V) grad(V). The first term retains the original surface
//! incidence; the second couples every boundary coordinate without creating
//! their dense Cartesian product. Storage is O(surface nonzeros + solid DOFs).
//! Applying the full Jacobian has the same linear complexity and fixed order.

use fs_exec::Cx;
use fs_math::det;
use fs_sparse::Csr;

use super::{HyperTetError, checkpoint, finite, invalid};

const TILE: usize = 256;

/// Full sealed-cavity force Jacobian on the original solid displacement chart.
///
/// This is a sparse matrix plus an explicitly retained rank-one contribution.
/// The dense acoustic coupling is never assembled. Entry access remains useful
/// for small reference checks; iterative solvers should use apply_into and reuse
/// their caller-owned output and scratch vectors.
#[derive(Debug, Clone)]
pub struct TetCavityGasTangent {
    geometric: Csr,
    acoustic_vector: Vec<f64>,
}

impl TetCavityGasTangent {
    pub(super) fn new(
        geometric: Csr,
        mut gradient: Vec<f64>,
        pressure_volume_derivative: f64,
        cx: &Cx<'_>,
    ) -> Result<Self, HyperTetError> {
        checkpoint(cx)?;
        if geometric.nrows() != gradient.len() || geometric.ncols() != gradient.len()
            || !pressure_volume_derivative.is_finite() || pressure_volume_derivative >= 0.0
        {
            return Err(invalid("cavity tangent needs the original chart and negative source dp/dV"));
        }
        let scale = finite(det::sqrt(-pressure_volume_derivative))?;
        let mut maximum = 0.0_f64;
        for (i, value) in gradient.iter_mut().enumerate() {
            if i % TILE == 0 { checkpoint(cx)?; }
            *value = finite(scale * *value)?;
            maximum = maximum.max(value.abs());
        }
        let rank_bound = finite(maximum * maximum)?;
        // Establish finite entry access without an O(n^2) product inspection.
        // The bound is conservative only at arithmetic overflow.
        for row in 0..geometric.nrows() {
            if row % TILE == 0 { checkpoint(cx)?; }
            let (_, values) = geometric.row(row);
            for (i, &value) in values.iter().enumerate() {
                if i % TILE == 0 { checkpoint(cx)?; }
                finite(value.abs() + rank_bound)?;
            }
        }
        checkpoint(cx)?;
        Ok(Self { geometric, acoustic_vector: gradient })
    }

    #[must_use]
    pub fn nrows(&self) -> usize { self.geometric.nrows() }

    #[must_use]
    pub fn ncols(&self) -> usize { self.geometric.ncols() }

    /// Sparse current-surface pressure stiffness, excluding compressibility.
    #[must_use]
    pub fn geometric_part(&self) -> &Csr { &self.geometric }

    /// w in J = J_geometry - w w^T, in the original full solid DOF ordering.
    /// This factor can also be used by a caller's low-rank preconditioner.
    #[must_use]
    pub fn acoustic_rank_one_vector(&self) -> &[f64] { &self.acoustic_vector }

    /// Full Jacobian entry, including off-graph gas compressibility coupling.
    #[must_use]
    pub fn get(&self, row: usize, column: usize) -> f64 {
        self.geometric.get(row, column)
            - self.acoustic_vector[row] * self.acoustic_vector[column]
    }

    /// Allocate one result and apply the complete pressure Jacobian.
    pub fn apply(&self, direction: &[f64], cx: &Cx<'_>) -> Result<Vec<f64>, HyperTetError> {
        checkpoint(cx)?;
        let mut result = vec![0.0; self.nrows()];
        self.apply_candidate(direction, &mut result, cx)?;
        Ok(result)
    }

    /// Heap-free full Jacobian action using caller-owned scratch.
    ///
    /// Every slice has one entry per full solid DOF. Refusal or cancellation
    /// leaves output unchanged; scratch may contain an unfinished candidate.
    /// Neither source coefficients nor geometry can change during this action.
    pub fn apply_into(
        &self,
        direction: &[f64],
        output: &mut [f64],
        scratch: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), HyperTetError> {
        checkpoint(cx)?;
        if output.len() != self.nrows() {
            return Err(invalid("cavity tangent output has the wrong solid chart"));
        }
        self.apply_candidate(direction, scratch, cx)?;
        output.copy_from_slice(scratch);
        Ok(())
    }

    fn apply_candidate(
        &self, direction: &[f64], candidate: &mut [f64], cx: &Cx<'_>,
    ) -> Result<(), HyperTetError> {
        if direction.len() != self.ncols() || candidate.len() != self.nrows() {
            return Err(invalid("cavity tangent direction or scratch has the wrong solid chart"));
        }
        let (mut projection, mut correction) = (0.0, 0.0);
        for (i, (&weight, &value)) in self.acoustic_vector.iter().zip(direction).enumerate() {
            if i % TILE == 0 { checkpoint(cx)?; }
            finite(value)?;
            let adjusted = finite(weight * value)? - correction;
            let next = finite(projection + adjusted)?;
            correction = (next - projection) - adjusted;
            projection = next;
        }
        for (row, output) in candidate.iter_mut().enumerate() {
            if row % TILE == 0 { checkpoint(cx)?; }
            let (columns, values) = self.geometric.row(row);
            let mut sparse = 0.0_f64;
            // The same ascending-column fused accumulation as Csr::spmv, with
            // cancellation inside long rows rather than an unbounded kernel.
            for (i, (&column, &value)) in columns.iter().zip(values).enumerate() {
                if i % TILE == 0 { checkpoint(cx)?; }
                sparse = finite(value.mul_add(direction[column], sparse))?;
            }
            *output = finite(sparse - finite(self.acoustic_vector[row] * projection)?)?;
        }
        checkpoint(cx)?;
        Ok(())
    }
}
