//! Directional diffusion without isotropic substitution. Tensor definiteness
//! and its inverse are enclosed using scale-normalized outward Sylvester minors.
//! The graph projection is only a proposal; the shared conservative forest and
//! exact simplex moments establish the same majorant as for scalar diffusion.
use super::*;

/// A symmetric positive-definite conductivity tensor in the mesh coordinates.
pub type ConductivityTensor = [[f64; 3]; 3];

/// Linear anisotropic diffusion on the same admitted tetrahedral domain.
/// Each tensor is element-constant. Boundary fluxes use the OUTWARD normal.
/// Entries denote their exact binary64 values, not uncertain material data.
#[derive(Debug, Clone, Copy)]
pub struct TensorTetProblem<'a> {
    pub vertices: &'a [[f64; 3]],
    pub tets: &'a [[usize; 4]],
    pub conductivity: &'a [ConductivityTensor],
    pub source: &'a [f64],
    pub boundary: &'a [BoundaryFace],
}

#[derive(Debug, Clone)]
pub(super) struct Tensor {
    matrix: [[Iv; 3]; 3],
    inverse: [[Iv; 3]; 3],
}

impl Tensor {
    fn new(k: ConductivityTensor) -> Result<Self, TetError> {
        if k.iter().flatten().any(|v| !v.is_finite())
            || (0..3).any(|i| k[i][i] <= 0.0 || (0..3).any(|j| k[i][j] != k[j][i])) {
            return Err(TetError::Invalid("finite exactly symmetric conductivity with positive diagonal required"));
        }
        // Scaling avoids underflow/overflow of determinants from physical units.
        // Divisions are themselves outward: rounded normalized entries must not
        // be mistaken for the exact matrix being certified.
        let scale = k[0][0].max(k[1][1]).max(k[2][2]);
        let m = k.map(|row| row.map(|v| Iv::point(v).div_pos(Iv::point(scale))));
        let [a, b, c] = m[0];
        let d = m[1][1];
        let e = m[1][2];
        let f = m[2][2];
        let c00 = d.mul(f).sub(e.sq());
        let c01 = c.mul(e).sub(b.mul(f));
        let c02 = b.mul(e).sub(c.mul(d));
        let c11 = a.mul(f).sub(c.sq());
        let c12 = b.mul(c).sub(a.mul(e));
        let c22 = a.mul(d).sub(b.sq());
        let det = a.mul(c00).add(b.mul(c01)).add(c.mul(c02));
        if [a, c22, det].iter().any(|v| v.is_unbounded() || v.lo <= 0.0) {
            return Err(TetError::Unsupported("positive definiteness not established by outward principal minors"));
        }
        let adj = [[c00, c01, c02], [c01, c11, c12], [c02, c12, c22]];
        let inverse = adj.map(|row| row.map(|v| v.div_pos(det).div_pos(Iv::point(scale))));
        if inverse.iter().flatten().any(Iv::is_unbounded) { return Err(TetError::Unbounded); }
        Ok(Self { matrix: k.map(|row| row.map(Iv::point)), inverse })
    }

    fn apply(&self, v: [Iv; 3]) -> [Iv; 3] {
        self.matrix.map(|row| dot(row, v))
    }

    fn inverse_quadratic(&self, v: [Iv; 3]) -> Iv {
        let mut q = Iv::zero();
        for i in 0..3 {
            q = q.add(self.inverse[i][i].mul(v[i].sq()));
            for j in (i+1)..3 {
                q = q.add(self.inverse[i][j].mul(v[i]).mul(v[j]).scale_pos(2.0));
            }
        }
        q
    }

    fn defect_integral(&self, nodal: &[[Iv; 3]; 4], volume: Iv) -> Iv {
        let sum = std::array::from_fn(|d| nodal.iter().fold(Iv::zero(), |s, v| s.add(v[d])));
        let diagonal = nodal.iter().fold(Iv::zero(), |s, &v| s.add(self.inverse_quadratic(v)));
        let mut value = volume.mul(diagonal.add(self.inverse_quadratic(sum)))
            .div_pos(Iv::point(20.0));
        // The exact quadratic form is nonnegative by the admitted SPD proof.
        // Do not clamp the upper endpoint or substitute a diagonal metric.
        if !value.is_unbounded() && value.hi >= 0.0 { value.lo = value.lo.max(0.0); }
        value
    }
}

#[derive(Clone, Copy)]
pub(super) enum Conductivity<'a> {
    Scalar(&'a [f64]),
    Tensor(&'a [Tensor]),
}
impl Conductivity<'_> {
    pub(super) fn len(self) -> usize {
        match self { Self::Scalar(k) => k.len(), Self::Tensor(k) => k.len() }
    }
    pub(super) fn validate(self, keep_going: &mut impl FnMut() -> bool) -> Result<(), TetError> {
        if let Self::Scalar(values) = self {
            for chunk in values.chunks(256) {
                poll(keep_going)?;
                if chunk.iter().any(|v| !v.is_finite() || *v <= 0.0) {
                    return Err(TetError::Invalid("finite positive conductivity required"));
                }
            }
        }
        Ok(())
    }
    pub(super) fn apply(self, e: usize, v: [Iv; 3]) -> [Iv; 3] {
        match self { Self::Scalar(k) => scale(v, Iv::point(k[e])), Self::Tensor(k) => k[e].apply(v) }
    }
    pub(super) fn flux_proposal(self, e: usize, gradient: [Iv; 3], normal_area: [Iv; 3]) -> Result<f64, TetError> {
        match self {
            Self::Scalar(k) => Ok(-k[e] * midpoint(dot(gradient, normal_area))?),
            Self::Tensor(k) => Ok(-midpoint(dot(k[e].apply(gradient), normal_area))?),
        }
    }
    pub(super) fn normal_coefficient(self, e: usize, normal_area: [Iv; 3], area: Iv) -> Result<f64, TetError> {
        match self {
            Self::Scalar(k) => Ok(k[e]),
            Self::Tensor(k) => {
                let value = midpoint(dot(normal_area, k[e].apply(normal_area)).div_pos(area.sq()))?;
                if value <= 0.0 { Err(TetError::Unbounded) } else { Ok(value) }
            }
        }
    }
    pub(super) fn defect_integral(self, e: usize, nodal: &[[Iv; 3]; 4], volume: Iv) -> Iv {
        match self {
            Self::Scalar(k) => {
                let squared: Vec<Iv> = (0..3).map(|d| {
                    integral_square(&nodal.map(|v| v[d]), volume, 20.0)
                }).collect();
                squared.into_iter().fold(Iv::zero(), Iv::add).div_pos(Iv::point(k[e]))
            }
            Self::Tensor(k) => k[e].defect_integral(nodal, volume),
        }
    }
    pub(super) fn bilinear_integral(self, e: usize, a: [Iv; 3], b: [Iv; 3], volume: Iv) -> Iv {
        match self {
            Self::Scalar(k) => volume.mul(Iv::point(k[e])).mul(dot(a,b)),
            Self::Tensor(k) => volume.mul(dot(k[e].apply(a), b)),
        }
    }
}

impl TensorTetProblem<'_> {
    fn prepare(&self, budget: FluxBudget, keep_going: &mut impl FnMut() -> bool) -> Result<Vec<Tensor>, TetError> {
        poll(keep_going)?;
        if self.tets.len() > budget.max_cells || self.vertices.len() > budget.max_cells.saturating_mul(4)
            || budget.max_iterations > 1_000_000 { return Err(TetError::Budget); }
        if self.conductivity.len() != self.tets.len() { return Err(TetError::Invalid("tensor conductivity length")); }
        let mut tensors = Vec::new();
        tensors.try_reserve_exact(self.tets.len()).map_err(|_| TetError::Budget)?;
        for &k in self.conductivity { poll(keep_going)?; tensors.push(Tensor::new(k)?); }
        Ok(tensors)
    }
    fn problem<'a>(&'a self, tensors: &'a [Tensor]) -> Problem<'a> {
        Problem { vertices: self.vertices, tets: self.tets, conductivity: Conductivity::Tensor(tensors),
            source: self.source, boundary: self.boundary }
    }
}

/// Bound error in the full tensor energy norm, including Robin trace energy.
/// Inconclusive positive-definiteness/inverse arithmetic refuses; tensors are
/// never symmetrized, clipped or replaced by their diagonal/isotropic mean.
pub fn tensor_energy_bound(
    problem: &TensorTetProblem<'_>, candidate: &[f64], budget: FluxBudget,
    mut keep_going: impl FnMut() -> bool,
) -> Result<EnergyBound, TetError> {
    let tensors = problem.prepare(budget, &mut keep_going)?;
    energy_bound_impl(&problem.problem(&tensors), candidate, budget, &mut keep_going)
}

/// Bound the exact declared weighted integral with the SAME tensor operator
/// for both primal and dual; retain the complete outward residual correction.
pub fn tensor_goal_bound(
    problem: &TensorTetProblem<'_>, candidate: &[f64], dual_candidate: &[f64],
    weights: &[f64], budget: FluxBudget, mut keep_going: impl FnMut() -> bool,
) -> Result<GoalBound, TetError> {
    let tensors = problem.prepare(budget, &mut keep_going)?;
    goal::goal_bound_impl(&problem.problem(&tensors), candidate, dual_candidate, weights, budget, &mut keep_going)
}

/// Whole-domain mean using the tensor unit-source dual and outward volume.
/// Does not establish a point maximum or material/model uncertainty bound.
pub fn tensor_mean_bound(
    problem: &TensorTetProblem<'_>, candidate: &[f64], dual_candidate: &[f64],
    budget: FluxBudget, mut keep_going: impl FnMut() -> bool,
) -> Result<MeanBound, TetError> {
    let tensors = problem.prepare(budget, &mut keep_going)?;
    goal::mean_bound_impl(&problem.problem(&tensors), candidate, dual_candidate, budget, &mut keep_going)
}

#[cfg(test)]
mod tests;
