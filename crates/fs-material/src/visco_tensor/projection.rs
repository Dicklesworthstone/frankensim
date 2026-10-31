//! Work-conjugate projection of the original three-dimensional material.
//! A column is a physical tensor strain per displacement coordinate. This is
//! fixed-kinematics quadrature, NOT plane-stress elimination or a loss fit.
use super::*;

mod transfer;
pub use transfer::TensorThermalState;
mod eigenstrain;
pub use eigenstrain::TensorThermalExpansion;

// Shared by observable-memory admission and the advertised capacity bound.
const THERMAL_EXPANSION_ENVELOPE_GUARD: f64 = 1.0 + 4096.0 * f64::EPSILON;

/// Cold projection limits; all columns and all six rows of every branch survive.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TensorProjectionBudget {
    /// Maximum displacement coordinates, in 1..=256.
    pub max_coordinates: usize,
    /// Bound on (branches+1+kv)*36*n + (1+kv)*6*n^2 scalar products,
    /// where kv=1 when the material supplies a Kelvin viscosity, otherwise 0.
    /// Thermal derivative projection separately bounds (branches+1)*42*n^2.
    pub max_scalar_terms: usize,
}

/// Derivatives of the original stored tensors through the SAME strain columns
/// and volume as a mechanical projection. Matrices are symmetric row-major;
/// signed and indefinite derivatives are valid and are never Gram-factored.
#[derive(Clone, Debug, PartialEq)]
pub struct ProjectedTensorTemperatureDerivative {
    equilibrium: Vec<f64>,
    branches: Vec<Vec<f64>>,
}
impl ProjectedTensorTemperatureDerivative {
    /// Complete equilibrium derivative, in N/(m K) for displacement columns.
    pub fn equilibrium(&self) -> &[f64] { &self.equilibrium }
    /// Complete Maxwell derivatives, in the original source branch order.
    pub fn branches(&self) -> &[Vec<f64>] { &self.branches }
}

/// One volume-scaled Maxwell factor A = sqrt(volume) B S in caller coordinates.
/// Its internal state z has energy |z|^2/2 [J], evolves as
/// z' = A displacement' - z/tau, and exerts resisting force A^T z.
#[derive(Clone, Debug, PartialEq)]
pub struct ProjectedPronyBranch {
    rows: Vec<Vec<f64>>,
    tau_s: f64,
}
impl ProjectedPronyBranch {
    /// Six unchanged factor rows, each with one entry per displacement coordinate.
    pub fn rows(&self) -> &[Vec<f64>] { &self.rows }
    /// Original material relaxation time [s].
    pub fn relaxation_time_s(&self) -> f64 { self.tau_s }
}

/// One material quadrature contribution bound to its original law and strain map.
/// Displacements are in metres when columns have units 1/m; the stiffness is
/// then N/m, factor rows sqrt(J)/m, and conjugate generalized forces newtons.
#[derive(Clone, Debug, PartialEq)]
pub struct ProjectedTensorProny {
    model: TensorProny,
    columns: Vec<Voigt>,
    volume_m3: f64,
    equilibrium: Vec<f64>,
    viscosity: Vec<f64>,
    viscosity_rows: Option<Vec<Vec<f64>>>,
    branches: Vec<ProjectedPronyBranch>,
    free_strain_mandel: Voigt,
}
impl TensorProny {
    /// Project a full 3D material through fixed tensor-shear strain columns.
    /// Positive volume includes the physical quadrature weight. All projected
    /// cross terms and rank-deficient branches are retained. No independent
    /// eigensolve, branchwise plane-stress approximation or damping ratio enters.
    /// A collection of these contributions can be assembled in one global basis.
    ///
    /// # Errors
    /// Invalid geometry/columns/budget, overflow, or cancellation. Partial
    /// factors are never published, and the material remains immutable.
    pub fn project_strain(&self, columns: &[Voigt], volume_m3: f64,
        budget: TensorProjectionBudget, mut cancelled: impl FnMut() -> bool)
        -> Result<ProjectedTensorProny>
    {
        poll(&mut cancelled)?;
        let n = columns.len();
        if !(1..=256).contains(&budget.max_coordinates) || n == 0 || n > budget.max_coordinates
            || !volume_m3.is_finite() || volume_m3 <= 0.0 || columns.iter().any(|v| !finite(v)) {
            return Err(invalid("projection needs bounded finite strain columns and positive volume"));
        }
        let tensors = 1 + usize::from(self.kv_factor.is_some());
        let work = (self.branches().len()+tensors).checked_mul(36).and_then(|v| v.checked_mul(n))
            .and_then(|v| n.checked_mul(n).and_then(|w| w.checked_mul(6*tensors)).and_then(|w| v.checked_add(w)))
            .ok_or_else(|| invalid("tensor projection work overflow"))?;
        if work > budget.max_scalar_terms { return Err(invalid("tensor projection exceeds scalar-work budget")); }
        let volume_root = det::sqrt(volume_m3);
        let mut project = |factor: &Tangent6| -> Result<Vec<Vec<f64>>> {
            let mut rows = vec![vec![0.0; n]; 6];
            for (j, column) in columns.iter().enumerate() {
                poll(&mut cancelled)?;
                let strain = mandel(column);
                for i in 0..6 {
                    rows[i][j] = volume_root*dot(&factor[i], &strain);
                    if !rows[i][j].is_finite() { return Err(invalid("projected material factor overflow")); }
                }
            }
            Ok(rows)
        };
        let equilibrium_rows = project(&self.equilibrium_gram)?;
        let mut branches = Vec::with_capacity(self.branches().len());
        for (branch, factor) in self.branches().iter().zip(&self.factors) {
            let gram_rows = branch.resolved_gram_rows(factor);
            branches.push(ProjectedPronyBranch {
                rows: project(&gram_rows)?,
                tau_s: branch.tau_s,
            });
        }
        let viscous_rows = self.kv_factor.as_ref().map(|l| {
            let rows = core::array::from_fn(|i|core::array::from_fn(|j|l[j][i]));
            project(&rows)
        }).transpose()?;
        let mut viscosity = vec![0.0; n*n];
        let mut equilibrium = vec![0.0; n*n];
        for i in 0..n {
            poll(&mut cancelled)?;
            for j in i..n {
                let value = equilibrium_rows.iter().fold(0.0, |s, row| s+row[i]*row[j]);
                if !value.is_finite() { return Err(invalid("projected equilibrium stiffness overflow")); }
                equilibrium[i*n+j] = value; equilibrium[j*n+i] = value;
                if let Some(rows) = &viscous_rows {
                    let value = rows.iter().fold(0.0,|s,row|s+row[i]*row[j]);
                    if !value.is_finite() {return Err(invalid("projected Kelvin viscosity overflow"));}
                    viscosity[i*n+j] = value; viscosity[j*n+i] = value;
                }
            }
        }
        poll(&mut cancelled)?;
        Ok(ProjectedTensorProny { model:self.clone(), columns:columns.to_vec(), volume_m3,
            equilibrium, viscosity, viscosity_rows: viscous_rows, branches, free_strain_mandel: [0.0; 6] })
    }
}
impl ProjectedTensorProny {
    /// Number of original displacement coordinates; no rank truncation occurs.
    pub fn coordinates(&self) -> usize { self.columns.len() }
    /// Original immutable material, including its small-strain validity limit.
    pub fn material(&self) -> &TensorProny { &self.model }
    /// Positive physical quadrature volume [m^3].
    pub fn volume_m3(&self) -> f64 { self.volume_m3 }
    /// Original tensor-shear strain columns, in unchanged coordinate order.
    pub fn strain_columns(&self) -> &[Voigt] { &self.columns }
    /// Full symmetric row-major equilibrium stiffness in the caller's coordinates.
    pub fn equilibrium_stiffness(&self) -> &[f64] { &self.equilibrium }
    /// Full symmetric Kelvin damping matrix [N s/m] in the same coordinates.
    /// Absent viscosity is exactly zero; signed cross terms are never dropped.
    pub fn viscosity_matrix(&self) -> &[f64] { &self.viscosity }
    /// Original projected Kelvin Gram rows, in sqrt(N s/m) for displacement
    /// coordinates. They retain zero/rank-deficient directions without a second
    /// factorization or eigenvalue repair by the mechanical consumer.
    pub fn viscosity_gram_rows(&self) -> Option<&[Vec<f64>]> { self.viscosity_rows.as_deref() }
    /// All volume-scaled memory factors in original material branch order.
    pub fn branches(&self) -> &[ProjectedPronyBranch] { &self.branches }

    /// Project dC/dT by V S^T (dC/dT) S using the original tensor-shear columns.
    /// The source owner establishes that the derivative belongs to this law
    /// and material frame. No strain recovery is needed: consumers retaining
    /// physical spring strains can use varying rank-deficient Maxwell tensors
    /// without losing their otherwise hidden strain history.
    ///
    /// This is a derivative of stored free energy only. Kelvin viscosity and
    /// relaxation-time derivatives contribute no mechanical entropy.
    pub fn thermal_stiffness_derivatives(&self,
        derivative: &source::thermal::TensorTemperatureDerivative,
        budget: TensorProjectionBudget, mut cancelled: impl FnMut() -> bool)
        -> Result<ProjectedTensorTemperatureDerivative>
    {
        poll(&mut cancelled)?;
        let n = self.coordinates();
        let count = 1 + self.branches.len();
        let work = n.checked_mul(n).and_then(|v| v.checked_mul(42))
            .and_then(|v| v.checked_mul(count))
            .ok_or_else(|| invalid("thermal projection work overflow"))?;
        if !(1..=256).contains(&budget.max_coordinates) || n > budget.max_coordinates
            || work > budget.max_scalar_terms
            || derivative.branches_pa_k.len() != self.branches.len() {
            return Err(invalid("thermal projection needs complete branches within its original work budget"));
        }
        let mut project = |matrix: &Tangent6| -> Result<Vec<f64>> {
            let mut result = vec![0.0; n*n];
            for i in 0..n {
                poll(&mut cancelled)?;
                let left = mandel(&self.columns[i]);
                for j in i..n {
                    let right = mandel(&self.columns[j]);
                    let value = self.volume_m3 * dot(&left, &apply(matrix, &right));
                    if !value.is_finite() { return Err(invalid("projected thermal derivative overflow")); }
                    result[i*n+j] = value; result[j*n+i] = value;
                }
            }
            Ok(result)
        };
        let equilibrium = project(&derivative.equilibrium_pa_k)?;
        let branches = derivative.branches_pa_k.iter().map(project)
            .collect::<Result<Vec<_>>>()?;
        poll(&mut cancelled)?;
        Ok(ProjectedTensorTemperatureDerivative { equilibrium, branches })
    }

    /// Reconstruct physical strain and enforce the original material domain.
    /// This checks a supplied state, not every point along a time interval.
    pub fn strain(&self, displacement: &[f64]) -> Result<Voigt> {
        if displacement.len() != self.coordinates() || displacement.iter().any(|v| !v.is_finite()) {
            return Err(invalid("displacements must match the finite material strain map"));
        }
        let mut out = [0.0; 6];
        for (column, &value) in self.columns.iter().zip(displacement) {
            for i in 0..6 { out[i] += column[i]*value; }
        }
        self.model.validate_mandel_strain(&mandel(&out))
            .map_err(|_| invalid("projected strain exceeds original small-strain domain"))?;
        if self.free_strain_mandel.iter().any(|x| *x != 0.0) {
            let elastic = core::array::from_fn(|i| WEIGHTS[i]*out[i] - self.free_strain_mandel[i]);
            self.model.validate_mandel_strain(&elastic)
                .map_err(|_| invalid("projected elastic strain exceeds original small-strain domain"))?;
        }
        Ok(out)
    }

    /// Transpose the SAME strain map: f_i = volume * stress : column_i.
    /// This is resisting/internal force, not a sign-flipped external load.
    pub fn conjugate_force(&self, stress: &Voigt) -> Result<Vec<f64>> {
        if !finite(stress) { return Err(invalid("nonfinite projected stress")); }
        let result: Vec<_> = self.columns.iter().map(|column|
            self.volume_m3*crate::tensor::contract(stress, column)).collect();
        if result.iter().any(|v| !v.is_finite()) { return Err(invalid("projected stress force overflow")); }
        Ok(result)
    }
}

#[cfg(test)]
mod tests;
