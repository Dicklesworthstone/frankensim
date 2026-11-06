//! Local spinodal modes and finite Gibbs tangent-plane distances.
//!
//! The local test examines ALL composition directions on the present-species
//! support. Congruence with sqrt(x) turns the ideal Hessian into the identity:
//! K = I + diag(sqrt(x)) chi diag(sqrt(x)), restricted perpendicular to sqrt(x).
//! This preserves inertia without overflowing 1/x for dilute components. A
//! Householder basis and bounded cyclic Jacobi solve give every eigenvalue.
//! They are entropy-metric eigenvalues, not Euclidean Hessian eigenvalues.
//!
//! Local positive curvature does not exclude metastability. A negative finite
//! tangent-plane distance supplies a separate instability witness. Testing
//! finitely many trial compositions cannot establish global stability; see
//! Michelsen, Fluid Phase Equilibria 9 (1982), 1-19,
//! <https://doi.org/10.1016/0378-3812(82)85001-2>.

use super::{MulticomponentSolutionErrorV1, MulticomponentSolutionStateV1, sum};
use fs_math::det::sqrt;

/// Explicit bounded work, scratch, and numerical resolution for the local test.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MulticomponentStabilityControlV1 {
    pub max_sweeps: usize,
    pub max_workspace_bytes: usize,
    pub relative_eigensolve_tolerance: f64,
    /// Absolute tolerance in the dimensionless ideal-entropy metric.
    pub stability_tolerance: f64,
}
impl Default for MulticomponentStabilityControlV1 {
    fn default() -> Self {
        Self {
            max_sweeps: 32,
            max_workspace_bytes: 512 * 1024,
            relative_eigensolve_tolerance: 1e-12,
            stability_tolerance: 1e-10,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum MulticomponentStabilityErrorV1 {
    InvalidControl,
    WorkspaceBudget { required_bytes: usize },
    WorkBudget { completed_sweeps: usize },
    Arithmetic(&'static str),
    Cancelled,
}
impl core::fmt::Display for MulticomponentStabilityErrorV1 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "multicomponent local stability: {self:?}")
    }
}
impl std::error::Error for MulticomponentStabilityErrorV1 {}
type Result<T> = core::result::Result<T, MulticomponentStabilityErrorV1>;

/// Numerical local assessment, never a global phase-equilibrium certificate.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MulticomponentLocalStabilityKindV1 {
    NoCompositionDirections,
    LocallyStable,
    Spinodal,
    Unresolved,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentLocalStabilityV1 {
    kind: MulticomponentLocalStabilityKindV1,
    eigenvalues: Box<[f64]>,
    direction: Box<[f64]>,
    numerical_uncertainty: f64,
    negative_modes: usize,
    completed_sweeps: usize,
}
impl MulticomponentLocalStabilityV1 {
    pub const fn kind(&self) -> MulticomponentLocalStabilityKindV1 {
        self.kind
    }
    /// Sorted eigenvalues on the present-species simplex in the ideal metric.
    pub fn eigenvalues(&self) -> &[f64] {
        &self.eigenvalues
    }
    pub fn smallest_eigenvalue(&self) -> Option<f64> {
        self.eigenvalues.first().copied()
    }
    /// Smallest-mode tangent in canonical species order; sum(d_i)=0 and
    /// sum(d_i^2/x_i)=1 to roundoff. Absent components remain exactly zero.
    /// The all-zero vector is returned when only one species is present.
    pub fn direction(&self) -> &[f64] {
        &self.direction
    }
    /// Remaining off-diagonal norm plus a floating-point assembly allowance.
    /// This is a numerical uncertainty estimate, not an interval certificate.
    pub const fn numerical_uncertainty(&self) -> f64 {
        self.numerical_uncertainty
    }
    pub const fn negative_modes(&self) -> usize {
        self.negative_modes
    }
    pub const fn completed_sweeps(&self) -> usize {
        self.completed_sweeps
    }
}

impl MulticomponentSolutionStateV1 {
    /// Evaluate g(y)/(RT) - sum_i y_i mu_i(x)/(RT), using one original law.
    /// A negative value is a finite demixing witness, including a metastable
    /// reference with positive local curvature. The trial must stay on the
    /// reference's present-species support; no absent species is manufactured.
    /// Pure-species standard terms cancel at the same temperature and pressure.
    pub fn tangent_plane_distance(
        &self,
        trial_fractions: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> core::result::Result<f64, MulticomponentSolutionErrorV1> {
        let trial = self
            .model()
            .evaluate(self.temperature(), trial_fractions, &mut cancelled)?;
        let mut distance = super::Sum::default();
        for i in 0..trial.mole_fractions().len() {
            super::poll(&mut cancelled)?;
            let y = trial.mole_fractions()[i];
            if y == 0.0 {
                continue;
            }
            let reference_mu = self.reduced_chemical_potentials()[i].ok_or(
                MulticomponentSolutionErrorV1::Invalid(
                    "tangent-plane trial adds an absent reference species",
                ),
            )?;
            let trial_mu = trial.reduced_chemical_potentials()[i].ok_or(
                MulticomponentSolutionErrorV1::Unrepresentable("trial chemical potential"),
            )?;
            distance.add(y * (trial_mu - reference_mu));
        }
        super::poll(&mut cancelled)?;
        super::finite(distance.value(), "Gibbs tangent-plane distance")
    }

    /// Examine the complete present-species simplex, with O(N^2) scratch and
    /// bounded O(max_sweeps*N^3) arithmetic. No concentration floor or 1/x
    /// matrix is formed. A positive result establishes only LOCAL curvature.
    pub fn local_stability(
        &self,
        control: MulticomponentStabilityControlV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentLocalStabilityV1> {
        check_cancel(&mut cancelled)?;
        if control.max_sweeps == 0
            || control.max_sweeps > 256
            || !control.relative_eigensolve_tolerance.is_finite()
            || control.relative_eigensolve_tolerance <= 0.0
            || control.relative_eigensolve_tolerance >= 1.0
            || !control.stability_tolerance.is_finite()
            || control.stability_tolerance < 0.0
        {
            return Err(MulticomponentStabilityErrorV1::InvalidControl);
        }
        let fractions = self.mole_fractions();
        let active: Vec<_> = fractions
            .iter()
            .enumerate()
            .filter_map(|(i, &x)| (x > 0.0).then_some(i))
            .collect();
        let n = active.len();
        let m = n.saturating_sub(1);
        let required_bytes =
            (n * n + 2 * m * m + 10 * n + 2 * fractions.len()) * core::mem::size_of::<f64>();
        if required_bytes > control.max_workspace_bytes {
            return Err(MulticomponentStabilityErrorV1::WorkspaceBudget { required_bytes });
        }
        if m == 0 {
            return Ok(MulticomponentLocalStabilityV1 {
                kind: MulticomponentLocalStabilityKindV1::NoCompositionDirections,
                eigenvalues: Box::new([]),
                direction: vec![0.0; fractions.len()].into_boxed_slice(),
                numerical_uncertainty: 0.0,
                negative_modes: 0,
                completed_sweeps: 0,
            });
        }
        let roots: Vec<_> = active.iter().map(|&i| sqrt(fractions[i])).collect();
        let root_norm = sqrt(sum(roots.iter().map(|v| v * v)));
        let anchor = (0..n)
            .max_by(|&i, &j| roots[i].total_cmp(&roots[j]).then(j.cmp(&i)))
            .expect("positive composition");
        let mut house: Vec<_> = roots.iter().map(|v| v / root_norm).collect();
        house[anchor] += 1.0;
        let hnorm = sqrt(sum(house.iter().map(|v| v * v)));
        for v in &mut house {
            *v /= hnorm;
        }
        let indices: Vec<_> = (0..n).filter(|&i| i != anchor).collect();
        let mut b = vec![0.0; n * n];
        let mut scale = 1.0_f64;
        for i in 0..n {
            check_cancel(&mut cancelled)?;
            for j in i + 1..n {
                let value = (self.chi(active[i], active[j]) * roots[i]) * roots[j];
                if !value.is_finite() {
                    return Err(MulticomponentStabilityErrorV1::Arithmetic(
                        "entropy Hessian",
                    ));
                }
                b[i * n + j] = value;
                b[j * n + i] = value;
                scale = scale.max(value.abs());
            }
        }
        for value in &mut b {
            *value /= scale;
        }
        let bh: Vec<_> = (0..n)
            .map(|i| sum((0..n).map(|j| b[i * n + j] * house[j])))
            .collect();
        let hbh = sum(house.iter().zip(&bh).map(|(h, b)| h * b));
        let mut a = vec![0.0; m * m];
        let mut vectors = vec![0.0; m * m];
        for (row, &i) in indices.iter().enumerate() {
            check_cancel(&mut cancelled)?;
            vectors[row * m + row] = 1.0;
            for (column, &j) in indices.iter().enumerate() {
                let ideal = if row == column { 1.0 / scale } else { 0.0 };
                a[row * m + column] = ideal + b[i * n + j]
                    - 2.0 * (house[i] * bh[j] + bh[i] * house[j])
                    + 4.0 * hbh * house[i] * house[j];
            }
        }
        let assembly_norm = a.iter().map(|v| v.abs()).fold(1.0, f64::max);
        let mut completed = 0;
        let mut residual = off_diagonal_norm(&a, m);
        while residual > control.relative_eigensolve_tolerance * assembly_norm {
            if completed == control.max_sweeps {
                return Err(MulticomponentStabilityErrorV1::WorkBudget {
                    completed_sweeps: completed,
                });
            }
            for p in 0..m {
                check_cancel(&mut cancelled)?;
                for q in p + 1..m {
                    let apq = a[p * m + q];
                    if apq == 0.0 {
                        continue;
                    }
                    let tau = (a[q * m + q] - a[p * m + p]) / (2.0 * apq);
                    let t = if tau.is_finite() {
                        let at = tau.abs();
                        let hyp = if at > 1.0 {
                            at * sqrt(1.0 + (1.0 / at) * (1.0 / at))
                        } else {
                            sqrt(1.0 + at * at)
                        };
                        let sign = if tau < 0.0 { -1.0 } else { 1.0 };
                        if at > f64::MAX / 2.0 {
                            sign * (0.5 / at)
                        } else {
                            sign / (at + hyp)
                        }
                    } else {
                        0.0
                    };
                    if t == 0.0 {
                        continue;
                    }
                    let c = 1.0 / sqrt(1.0 + t * t);
                    let s = t * c;
                    a[p * m + p] -= t * apq;
                    a[q * m + q] += t * apq;
                    a[p * m + q] = 0.0;
                    a[q * m + p] = 0.0;
                    for k in 0..m {
                        if k != p && k != q {
                            let kp = a[k * m + p];
                            let kq = a[k * m + q];
                            let np = c * kp - s * kq;
                            let nq = s * kp + c * kq;
                            a[k * m + p] = np;
                            a[p * m + k] = np;
                            a[k * m + q] = nq;
                            a[q * m + k] = nq;
                        }
                        let vp = vectors[k * m + p];
                        let vq = vectors[k * m + q];
                        vectors[k * m + p] = c * vp - s * vq;
                        vectors[k * m + q] = s * vp + c * vq;
                    }
                }
            }
            completed += 1;
            residual = off_diagonal_norm(&a, m);
        }
        check_cancel(&mut cancelled)?;
        let mut order: Vec<_> = (0..m).collect();
        order.sort_by(|&i, &j| a[i * m + i].total_cmp(&a[j * m + j]).then(i.cmp(&j)));
        let eigenvalues: Vec<_> = order.iter().map(|&i| a[i * m + i] * scale).collect();
        let uncertainty = (residual + 256.0 * f64::EPSILON * n as f64 * assembly_norm) * scale;
        if eigenvalues.iter().any(|v| !v.is_finite()) || !uncertainty.is_finite() {
            return Err(MulticomponentStabilityErrorV1::Arithmetic(
                "stability spectrum",
            ));
        }
        let column = order[0];
        let projection = sum(indices
            .iter()
            .enumerate()
            .map(|(row, &i)| house[i] * vectors[row * m + column]));
        let mut direction = vec![0.0; fractions.len()];
        for (row, &i) in indices.iter().enumerate() {
            direction[active[i]] =
                roots[i] * (vectors[row * m + column] - 2.0 * house[i] * projection);
        }
        direction[active[anchor]] = -sum(active
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != anchor)
            .map(|(_, &i)| direction[i]));
        let norm = sqrt(sum(active.iter().enumerate().map(|(i, &j)| {
            let v = direction[j] / roots[i];
            v * v
        })));
        let first = direction.iter().copied().find(|v| *v != 0.0).unwrap_or(1.0);
        let orient = if first < 0.0 { -norm } else { norm };
        for value in &mut direction {
            *value /= orient;
        }
        if !norm.is_finite() || norm <= 0.0 || direction.iter().any(|v| !v.is_finite()) {
            return Err(MulticomponentStabilityErrorV1::Arithmetic(
                "stability direction",
            ));
        }
        let resolution = uncertainty + control.stability_tolerance;
        let negative_modes = eigenvalues.iter().filter(|&&v| v < -resolution).count();
        let kind = if negative_modes > 0 {
            MulticomponentLocalStabilityKindV1::Spinodal
        } else if eigenvalues[0] > resolution {
            MulticomponentLocalStabilityKindV1::LocallyStable
        } else {
            MulticomponentLocalStabilityKindV1::Unresolved
        };
        check_cancel(&mut cancelled)?;
        Ok(MulticomponentLocalStabilityV1 {
            kind,
            eigenvalues: eigenvalues.into_boxed_slice(),
            direction: direction.into_boxed_slice(),
            numerical_uncertainty: uncertainty,
            negative_modes,
            completed_sweeps: completed,
        })
    }
}

fn off_diagonal_norm(matrix: &[f64], n: usize) -> f64 {
    sqrt(sum((0..n).flat_map(|i| {
        (i + 1..n).map(move |j| 2.0 * matrix[i * n + j] * matrix[i * n + j])
    })))
}
fn check_cancel(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(MulticomponentStabilityErrorV1::Cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests;
