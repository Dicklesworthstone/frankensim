//! Geometry maps for the original Maxwell internal-strain storage. This file
//! supplies no time integrator: every coordinate is advanced by the existing
//! port-Hamiltonian owner, with its original work and resistance ledger.
use super::{PhsError, PortHamiltonian};
use fs_math::det;
use std::rc::Rc;

/// One scalar strain and its exact derivatives on the original state prefix.
/// The map is immutable geometry, never an evolving constitutive history.
/// Compatible units require `stiffness * strain^2` to be energy. Implementations
/// fill the entire derivative output and return nonfinite values on invalid
/// arguments. Source validity and admissible strain/history domains belong to
/// the composing physical owner, as for linear `RelaxationBranch`.
pub trait RelaxationStrainMap: core::fmt::Debug {
    fn state_dimension(&self) -> usize;
    fn strain(&self, state: &[f64]) -> f64;
    fn gradient(&self, state: &[f64], output: &mut [f64]);
    fn hessian_vector(&self, state: &[f64], direction: &[f64], output: &mut [f64]);
    /// Explicit declaration of an at-most-quadratic STRAIN map. This does not
    /// make its Maxwell energy quadratic in the complete state. It permits a
    /// cancellation-stable exact Gonzalez remainder using the same derivatives.
    fn is_quadratic_strain(&self) -> bool { false }
    /// Optional exact linear polynomial row of a zero-constant quadratic
    /// strain: epsilon(x) = b.x + quadratic(x), with epsilon(0)=0.
    /// This is original source algebra, never a sampled tangent. It lets the
    /// complete axial/history small-signal block propagate exactly while all
    /// geometric coupling stays implicit. None keeps the existing split.
    fn zero_constant_linear_strain(&self) -> Option<Vec<f64>> { None }
}

/// A Maxwell arm on an explicitly supplied, possibly nonlinear strain map.
/// `z=sqrt(stiffness)*viscous_strain` is appended to the original state. Its
/// energy and loss are the SAME scalar law as the linear projection API.
#[derive(Clone, Debug)]
pub struct MappedRelaxationBranch {
    pub strain: Rc<dyn RelaxationStrainMap>,
    pub stiffness: f64,
    pub relaxation_time_s: f64,
}

/// Stable midpoint work remainder for `(scale * epsilon(x) - z)^2 / 2`.
///
/// This uses the original declared quadratic strain and its exact derivatives;
/// it does not approximate the strain or advance a constitutive history.
/// `offsets` contains the actual scalar z coordinate at a, b and mid. A
/// physical equilibrium strain energy `k * epsilon^2 / 2` uses scale 1,
/// zero offsets, and multiplies the returned remainder by its original k.
/// Midpoint rounding is retained. No subtraction of total stored energies is
/// used to recover small cubic/quartic work. Finite complete coordinates,
/// positive finite scale and an explicit quadratic declaration are required.
pub fn quadratic_strain_midpoint_energy_remainder(
    strain: &dyn RelaxationStrainMap,
    scale: f64,
    a: &[f64],
    b: &[f64],
    mid: &[f64],
    offsets: [f64; 3],
) -> Option<f64> {
    let n = strain.state_dimension();
    if n == 0
        || a.len() != n
        || b.len() != n
        || mid.len() != n
        || !scale.is_finite()
        || scale <= 0.0
        || a.iter()
            .chain(b)
            .chain(mid)
            .chain(&offsets)
            .any(|v| !v.is_finite())
        || !strain.is_quadratic_strain()
    {
        return None;
    }
    let [za, zb, zm] = offsets;
    let delta: Vec<_> = b.iter().zip(a).map(|(b, a)| b - a).collect();
    let left: Vec<_> = a.iter().zip(mid).map(|(a, m)| a - m).collect();
    let right: Vec<_> = b.iter().zip(mid).map(|(b, m)| b - m).collect();
    let mut g = vec![0.0; mid.len()];
    let mut hd = vec![0.0; mid.len()];
    let mut ha = vec![0.0; mid.len()];
    let mut hb = vec![0.0; mid.len()];
    strain.gradient(mid, &mut g);
    strain.hessian_vector(mid, &delta, &mut hd);
    strain.hessian_vector(mid, &left, &mut ha);
    strain.hessian_vector(mid, &right, &mut hb);
    let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(a, b)| a * b).sum::<f64>();
    let rounding: Vec<_> = left.iter().zip(&right).map(|(a, b)| a + b).collect();
    // Retain midpoint rounding, not only the ideal symmetric-endpoint
    // formula. With exact midpoint, cross=0 and the first factor is
    // Delta x^T Hess(sqrt(k)*epsilon) Delta x / 8.
    let cross = 0.5 * scale * dot(&rounding, &hd);
    let change = scale * dot(&g, &delta) - (zb - za) + cross;
    let offset = scale * (dot(&g, &rounding) + 0.5 * (dot(&left, &ha) + dot(&right, &hb)))
        - ((za - zm) + (zb - zm));
    let effort = scale * strain.strain(mid) - zm;
    let result = (0.5 * offset) * change + effort * cross;
    result.is_finite().then_some(result)
}

pub(super) enum RelaxationCoupling {
    Linear(Vec<f64>),
    Mapped { strain: Rc<dyn RelaxationStrainMap>, scale: f64 },
}
impl RelaxationCoupling {
    pub(super) fn mapped_remainder(&self,a:&[f64],b:&[f64],mid:&[f64],
        za:f64,zb:f64,zm:f64)->Option<f64>
    {
        let Self::Mapped {strain,scale}=self else {return None;};
        quadratic_strain_midpoint_energy_remainder(strain.as_ref(), *scale, a, b, mid, [za,zb,zm])
    }
    pub(super) fn strain(&self, x: &[f64]) -> f64 {
        match self {
            Self::Linear(row) => row.iter().zip(x).map(|(b,x)|b*x).sum(),
            Self::Mapped { strain, scale } => scale*strain.strain(x),
        }
    }
    pub(super) fn add_gradient(&self, x: &[f64], effort: f64, out: &mut [f64]) {
        match self {
            Self::Linear(row) => {
                for (o,b) in out.iter_mut().zip(row) { *o += b*effort; }
            }
            Self::Mapped { strain, scale } => {
                let mut gradient=vec![0.0;x.len()];
                strain.gradient(x,&mut gradient);
                for (o,b) in out.iter_mut().zip(gradient) { *o += (scale*b)*effort; }
            }
        }
    }
}
// Differentiate the SAME full-state Maxwell energy used by the original
// gradient, including prestress curvature and every explicit memory row.
// Linear arms need no scratch; mapped arms reuse two vectors across all arms.
impl super::ContinuousStorageTangent for super::RelaxationStorage {
    fn apply(&self, state: &[f64], direction: &[f64], out: &mut [f64]) -> bool {
        let Some(dimension) = self.base_dim.checked_add(self.couplings.len()) else {
            return false;
        };
        if state.len() != dimension
            || direction.len() != dimension
            || out.len() != dimension
            || state.iter().chain(direction).any(|v| !v.is_finite())
        {
            return false;
        }
        let Some(base) = self.base.continuous_tangent() else {
            return false;
        };
        out.fill(f64::NAN);
        if !base.apply(
            &state[..self.base_dim],
            &direction[..self.base_dim],
            &mut out[..self.base_dim],
        ) || out[..self.base_dim].iter().any(|v| !v.is_finite())
        {
            return false;
        }
        let mut mapped_scratch: Option<(Vec<f64>, Vec<f64>)> = None;
        for (arm, coupling) in self.couplings.iter().enumerate() {
            let memory = self.base_dim + arm;
            let mut change = -direction[memory];
            match coupling {
                RelaxationCoupling::Linear(row) => {
                    for (&weight, &v) in row.iter().zip(&direction[..self.base_dim]) {
                        change += weight * v;
                    }
                    for (value, &weight) in out[..self.base_dim].iter_mut().zip(row) {
                        *value += weight * change;
                    }
                }
                RelaxationCoupling::Mapped { strain, scale } => {
                    let (gradient, curvature) = mapped_scratch.get_or_insert_with(|| {
                        (vec![f64::NAN; self.base_dim], vec![f64::NAN; self.base_dim])
                    });
                    gradient.fill(f64::NAN);
                    curvature.fill(f64::NAN);
                    let x = &state[..self.base_dim];
                    let v = &direction[..self.base_dim];
                    strain.gradient(x, gradient);
                    strain.hessian_vector(x, v, curvature);
                    let effort = scale * strain.strain(x) - state[memory];
                    if !effort.is_finite() {
                        return false;
                    }
                    for i in 0..self.base_dim {
                        gradient[i] *= scale;
                        curvature[i] *= scale;
                        if !gradient[i].is_finite() || !curvature[i].is_finite() {
                            return false;
                        }
                        change += gradient[i] * v[i];
                    }
                    for i in 0..self.base_dim {
                        out[i] += gradient[i] * change + effort * curvature[i];
                    }
                }
            }
            out[memory] = -change;
        }
        out.iter().all(|v| v.is_finite())
    }
}

impl PortHamiltonian {
    /// Append nonlinear Maxwell strain maps through the original H/R owner.
    ///
    /// Each arm stores `(sqrt(k)*epsilon(x)-z)^2/2`, gives the exact reaction
    /// `-sqrt(k)*(sqrt(k)*epsilon-z)*gradient(epsilon)` on canonical momenta,
    /// and has `R_zz=1/tau`. Existing ports, state order and original storage
    /// remain. This is a nonlinear pHS system. An explicitly quadratic strain
    /// map retains the original base quadratic block and each z²/2 storage for
    /// the existing exponential discrete-gradient owner. An explicit
    /// zero-constant linear strain row also retains its exact Maxwell dyad and
    /// cross-memory terms; all geometric strain/history coupling remains implicit. The complete energy is never marked quadratic.
    /// The ordinary discrete-gradient owner is unchanged: its accepted loss is
    /// `dt*e_z^2/tau`. The exponential owner instead integrates the SAME positive
    /// resistance along its accepted held-remainder trajectory. Only its linear
    /// propagation is exact; nonlinear accuracy still needs time refinement.
    /// Mapped evaluation currently uses bounded temporary derivative vectors;
    /// this reference nonlinear path does not claim allocation-free execution.
    /// If the original storage declares a continuous tangent, the complete
    /// relaxed storage retains it. Its exact action includes the geometric
    /// curvature under stored stress and every mechanical/memory cross term.
    /// Two derivative vectors are reused across mapped arms per action;
    /// purely linear relaxation arms require no derivative scratch.
    ///
    /// Geometry, material bands, clocks and history bounds must be admitted by
    /// the caller. Nonlinear discrete-gradient evolution does not imply that
    /// history remains in a convex endpoint-strain hull. No finite-strain or
    /// material authority follows just from implementing this map.
    pub fn with_mapped_relaxation_branches(self, branches: Vec<MappedRelaxationBranch>)
        -> Result<Self, PhsError>
    {
        if branches.is_empty() {return Ok(self);}
        let mut couplings=Vec::with_capacity(branches.len());
        let mut rates=Vec::with_capacity(branches.len());
        for branch in branches {
            if branch.strain.state_dimension()!=self.n {
                return Err(PhsError::Dimension {what:"relaxation strain map state dimensions"});
            }
            let scale=det::sqrt(branch.stiffness);
            let rate=1.0/branch.relaxation_time_s;
            if !branch.stiffness.is_finite() || branch.stiffness<=0.0
                || !branch.relaxation_time_s.is_finite() || branch.relaxation_time_s<=0.0
                || !scale.is_finite() || scale<=0.0 || !rate.is_finite() || rate<=0.0 {
                return Err(PhsError::RelaxationParameters);
            }
            couplings.push(RelaxationCoupling::Mapped {strain:branch.strain,scale});
            rates.push(rate);
        }
        self.with_relaxation_couplings(couplings,rates)
    }
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod exponential_tests;

#[cfg(test)]
mod square_work_tests;
