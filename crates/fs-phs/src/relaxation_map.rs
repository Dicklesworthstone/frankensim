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

pub(super) enum RelaxationCoupling {
    Linear(Vec<f64>),
    Mapped { strain: Rc<dyn RelaxationStrainMap>, scale: f64 },
}
impl RelaxationCoupling {
    pub(super) fn mapped_remainder(&self,a:&[f64],b:&[f64],mid:&[f64],
        za:f64,zb:f64,zm:f64)->Option<f64>
    {
        let Self::Mapped {strain,scale}=self else {return None;};
        if !strain.is_quadratic_strain() {return None;}
        let delta:Vec<_>=b.iter().zip(a).map(|(b,a)|b-a).collect();
        let left:Vec<_>=a.iter().zip(mid).map(|(a,m)|a-m).collect();
        let right:Vec<_>=b.iter().zip(mid).map(|(b,m)|b-m).collect();
        let mut g=vec![0.0;mid.len()];let mut hd=vec![0.0;mid.len()];
        let mut ha=vec![0.0;mid.len()];let mut hb=vec![0.0;mid.len()];
        strain.gradient(mid,&mut g);strain.hessian_vector(mid,&delta,&mut hd);
        strain.hessian_vector(mid,&left,&mut ha);strain.hessian_vector(mid,&right,&mut hb);
        let dot=|a:&[f64],b:&[f64]|a.iter().zip(b).map(|(a,b)|a*b).sum::<f64>();
        let rounding:Vec<_>=left.iter().zip(&right).map(|(a,b)|a+b).collect();
        // Retain midpoint rounding, not only the ideal symmetric-endpoint
        // formula. With exact midpoint, cross=0 and the first factor is
        // Delta x^T Hess(sqrt(k)*epsilon) Delta x / 8.
        let cross=0.5*scale*dot(&rounding,&hd);
        let change=scale*dot(&g,&delta)-(zb-za)+cross;
        let offset=scale*(dot(&g,&rounding)+0.5*(dot(&left,&ha)+dot(&right,&hb)))
            -((za-zm)+(zb-zm));
        let effort=scale*strain.strain(mid)-zm;
        Some((0.5*offset)*change+effort*cross)
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
impl PortHamiltonian {
    /// Append nonlinear Maxwell strain maps through the original H/R owner.
    ///
    /// Each arm stores `(sqrt(k)*epsilon(x)-z)^2/2`, gives the exact reaction
    /// `-sqrt(k)*(sqrt(k)*epsilon-z)*gradient(epsilon)` on canonical momenta,
    /// and has `R_zz=1/tau`. Existing ports, state order and original storage
    /// remain. This is a nonlinear pHS system: a quadratic/exact-flow claim is
    /// deliberately unavailable, and the existing discrete-gradient solver
    /// advances mechanical and memory coordinates together. Its accepted
    /// dissipated energy is `dt*e_z^2/tau`, using the full discrete effort.
    /// Mapped evaluation currently uses bounded temporary derivative vectors;
    /// this reference nonlinear path does not claim allocation-free execution.
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
