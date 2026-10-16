//! Passive anisotropic generalized-Maxwell material memory (IRP M0.2-W2a).
//!
//! Constitutive coefficients use Mandel order [xx, yy, zz, sqrt(2)xy,
//! sqrt(2)yz, sqrt(2)zx]. Strain/stress at the time-domain boundary use the
//! crate's tensor-shear Voigt convention. These are NOT engineering shears.
//! Each branch is C_j = B_j^T B_j, allowing positive-semidefinite, coupled and
//! rank-deficient tensors without eigenvalue clipping or simultaneous
//! diagonalization. q_j' = B_j epsilon' - q_j/tau_j, sigma_j = B_j^T q_j.
//! Free energy is (|B_eq epsilon|^2 + sum |q_j|^2)/2; heat rate is
//! sum |q_j|^2/tau_j. No empirical damping ratio is inferred.
//!
//! The ramp update integrates prescribed piecewise-linear strain exactly
//! (up to floating-point arithmetic), including work and heat separately.
//! It is not an unconditional stability claim for a coupled time integrator.
//! Parameters and the small-strain validity bound are caller declarations;
//! no measured-material, temperature-shift or fractional-lowering claim.
use crate::{Tangent6, Voigt};
use fs_math::det;
use std::sync::Arc;

mod reduction;
pub use reduction::{Complex3, HarmonicPlateSection, PlaneStressModulus};

/// Maximum retained relaxation branches per material.
pub const MAX_BRANCHES: usize = 256;
const WEIGHTS: Voigt = [1.0, 1.0, 1.0, core::f64::consts::SQRT_2,
    core::f64::consts::SQRT_2, core::f64::consts::SQRT_2];

/// Constitutive refusal; failed advances leave the complete history unchanged.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum TensorPronyError {
    /// Invalid data, an exceeded validity bound or unrepresentable arithmetic.
    Invalid(&'static str),
    /// Cooperative cancellation, distinct from material inadmissibility.
    Cancelled,
}
impl core::fmt::Display for TensorPronyError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(what) => write!(f, "FS-MAT-TENSOR-PRONY: {what}"),
            Self::Cancelled => write!(f, "FS-MAT-TENSOR-PRONY-CANCELLED"),
        }
    }
}
impl std::error::Error for TensorPronyError {}
type Result<T> = core::result::Result<T, TensorPronyError>;
fn invalid(what: &'static str) -> TensorPronyError { TensorPronyError::Invalid(what) }
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() { Err(TensorPronyError::Cancelled) } else { Ok(()) }
}
fn dot(a: &Voigt, b: &Voigt) -> f64 {
    a.iter().zip(b).fold(0.0, |sum, (a,b)| sum + a*b)
}
fn apply(a: &Tangent6, x: &Voigt) -> Voigt { core::array::from_fn(|i| dot(&a[i], x)) }
fn transpose_add(b: &Tangent6, q: &Voigt, out: &mut Voigt) {
    for i in 0..6 { for j in 0..6 { out[i] += b[j][i]*q[j]; } }
}
fn mandel(v: &Voigt) -> Voigt { core::array::from_fn(|i| WEIGHTS[i]*v[i]) }
fn tensor(v: &Voigt) -> Voigt { core::array::from_fn(|i| v[i]/WEIGHTS[i]) }
fn finite(v: &Voigt) -> bool { v.iter().all(|x| x.is_finite()) }
fn gram(b: &Tangent6) -> Result<Tangent6> {
    if b.iter().flatten().any(|v| !v.is_finite()) { return Err(invalid("nonfinite Gram factor")); }
    let mut c = [[0.0;6];6];
    for i in 0..6 { for j in i..6 {
        let value = (0..6).fold(0.0, |sum,k| sum+b[k][i]*b[k][j]);
        if !value.is_finite() { return Err(invalid("unrepresentable stiffness")); }
        c[i][j] = value; c[j][i] = value;
    } }
    Ok(c)
}
// Strict fixed-order Cholesky: no tolerance-based promotion of negative or
// zero pivots. Rank-deficient branch inputs use from_gram_rows instead.
fn positive_factor(c: &Tangent6) -> Result<Tangent6> {
    for i in 0..6 { for j in 0..6 {
        if !c[i][j].is_finite() || c[i][j] != c[j][i] {
            return Err(invalid("stiffness must be finite and symmetric in Mandel coordinates"));
        }
    } }
    let mut l = [[0.0;6];6];
    for i in 0..6 { for j in 0..=i {
        let mut value = c[i][j];
        for k in 0..j { value -= l[i][k]*l[j][k]; }
        if i == j {
            if !value.is_finite() || value <= 0.0 { return Err(invalid("stiffness must be positive definite")); }
            l[i][j] = det::sqrt(value);
        } else { l[i][j] = value/l[j][j]; }
        if !l[i][j].is_finite() { return Err(invalid("unrepresentable stiffness factor")); }
    } }
    Ok(core::array::from_fn(|i| core::array::from_fn(|j| l[j][i])))
}

/// One Maxwell tensor. Factor rows have units sqrt(Pa); tau is in seconds.
#[derive(Clone, Debug, PartialEq)]
pub struct TensorPronyBranch { factor: Tangent6, stiffness: Tangent6, tau_s: f64 }
impl TensorPronyBranch {
    /// Construct any PSD tensor, including zero/rank-deficient tensors, as B^T B.
    pub fn from_gram_rows(factor: Tangent6, tau_s: f64) -> Result<Self> {
        if !tau_s.is_finite() || tau_s <= 0.0 { return Err(invalid("relaxation time must be finite and positive")); }
        let stiffness = gram(&factor)?;
        Ok(Self { factor, stiffness, tau_s })
    }
    /// Convenience for a strictly SPD Mandel stiffness; no repair is performed.
    pub fn from_stiffness_mandel(stiffness: Tangent6, tau_s: f64) -> Result<Self> {
        Self::from_gram_rows(positive_factor(&stiffness)?, tau_s)
    }
    /// Actual factor-defined Mandel stiffness [Pa].
    pub fn stiffness_mandel(&self) -> &Tangent6 { &self.stiffness }
    /// Relaxation time [s].
    pub fn relaxation_time_s(&self) -> f64 { self.tau_s }
}

#[derive(Clone, Debug, PartialEq)]
struct Parameters { equilibrium: Tangent6, equilibrium_factor: Tangent6,
    branches: Vec<TensorPronyBranch>, strain_limit: f64 }

/// Immutable material shared by its independent, model-bound histories.
#[derive(Clone, Debug, PartialEq)]
pub struct TensorProny { parameters: Arc<Parameters> }
impl TensorProny {
    /// Admit an SPD equilibrium tensor and up to MAX_BRANCHES PSD branches.
    /// `strain_limit` bounds the Frobenius norm of the small-strain tensor.
    /// Factorization changes the represented equilibrium by rounding only.
    pub fn new_mandel(equilibrium: Tangent6, branches: Vec<TensorPronyBranch>,
        strain_limit: f64) -> Result<Self> {
        if branches.len() > MAX_BRANCHES || !strain_limit.is_finite() || strain_limit <= 0.0 {
            return Err(invalid("invalid branch count or positive finite strain bound"));
        }
        let equilibrium_factor = positive_factor(&equilibrium)?;
        let equilibrium = gram(&equilibrium_factor)?;
        // The instantaneous tensor must be representable too.
        let mut instantaneous = equilibrium;
        for branch in &branches { for i in 0..6 { for j in 0..6 {
            instantaneous[i][j] += branch.stiffness[i][j];
            if !instantaneous[i][j].is_finite() { return Err(invalid("instantaneous stiffness overflow")); }
        } } }
        Ok(Self { parameters: Arc::new(Parameters { equilibrium, equilibrium_factor, branches, strain_limit }) })
    }
    /// Convert an existing tensor-shear stress/strain tangent (e.g. elastic::
    /// IsotropicElastic::stiffness) without doubling shear stiffness.
    pub fn from_tensor_shear(equilibrium: Tangent6, branches: Vec<TensorPronyBranch>,
        strain_limit: f64) -> Result<Self> {
        let c = core::array::from_fn(|i| core::array::from_fn(|j|
            (WEIGHTS[i]/WEIGHTS[j])*equilibrium[i][j]));
        Self::new_mandel(c, branches, strain_limit)
    }
    /// Retained branches, in caller order; none are pruned or fitted.
    pub fn branches(&self) -> &[TensorPronyBranch] { &self.parameters.branches }
    /// Actual equilibrium Mandel stiffness [Pa].
    pub fn equilibrium_mandel(&self) -> &Tangent6 { &self.parameters.equilibrium }
    /// Virgin zero-strain, relaxed history. No state can be rebound to another law.
    pub fn history(&self) -> TensorPronyHistory {
        TensorPronyHistory { model: self.clone(), strain: [0.0;6], stress: [0.0;6],
            q: vec![[0.0;6];self.branches().len()], stored: 0.0, work: 0.0, heat: 0.0 }
    }
    /// Full storage/loss tensors for exp(+i omega t), matching visco::modulus.
    /// The loss matrix is PSD; its OFF-DIAGONAL entries need not be positive.
    pub fn modulus_mandel(&self, omega: f64) -> Result<(Tangent6, Tangent6)> {
        if !omega.is_finite() || omega < 0.0 { return Err(invalid("frequency must be finite and nonnegative")); }
        let mut storage = self.parameters.equilibrium;
        let mut loss = [[0.0;6];6];
        for branch in self.branches() {
            let z = omega*branch.tau_s;
            // Stable at both ends, including overflow of positive omega*tau.
            let (re, im) = if z <= 1.0 { (z*(z/(1.0+z*z)), z/(1.0+z*z)) }
                else { let u=1.0/z; (1.0/(1.0+u*u), u/(1.0+u*u)) };
            for i in 0..6 { for j in 0..6 {
                storage[i][j] += re*branch.stiffness[i][j];
                loss[i][j] += im*branch.stiffness[i][j];
            } }
        }
        if storage.iter().flatten().chain(loss.iter().flatten()).any(|v| !v.is_finite()) {
            return Err(invalid("modulus overflow"));
        }
        Ok((storage,loss))
    }
    /// Exact endpoint stress derivative for a prescribed linear strain ramp.
    /// Returned tensor uses the crate's tensor-shear Voigt convention.
    pub fn tangent(&self, dt_s: f64) -> Result<Tangent6> {
        if !dt_s.is_finite() || dt_s < 0.0 { return Err(invalid("time step must be finite and nonnegative")); }
        let mut c = self.parameters.equilibrium;
        for branch in self.branches() {
            let r = dt_s/branch.tau_s;
            let b = if r == 0.0 {1.0} else {-det::expm1(-r)/r};
            for i in 0..6 { for j in 0..6 { c[i][j] += b*branch.stiffness[i][j]; } }
        }
        for i in 0..6 { for j in 0..6 { c[i][j] *= WEIGHTS[j]/WEIGHTS[i]; } }
        if c.iter().flatten().any(|v| !v.is_finite()) { return Err(invalid("tangent overflow")); }
        Ok(c)
    }
}

/// Exact-ramp endpoint and independently integrated energy densities [J/m^3].
#[derive(Clone, Debug, PartialEq)]
pub struct TensorPronyStep {
    /// Endpoint Cauchy stress in tensor-shear Voigt order [Pa].
    pub stress: Voigt,
    /// Endpoint total recoverable energy density.
    pub stored_energy: f64,
    /// External stress : strain-increment work on this step (can be negative).
    pub work: f64,
    /// Dashpot heat on this step, integrated as sums of squares.
    pub dissipated: f64,
    /// work - change in stored energy - heat; no balancing term is inserted.
    pub balance_residual: f64,
}

/// Private branch coordinates prevent shape/model mismatches and manual damage.
#[derive(Clone, Debug, PartialEq)]
pub struct TensorPronyHistory {
    model: TensorProny, strain: Voigt, stress: Voigt, q: Vec<Voigt>,
    stored: f64, work: f64, heat: f64,
}
impl TensorPronyHistory {
    /// Current physical strain in tensor-shear Voigt order.
    pub fn strain(&self) -> Voigt { tensor(&self.strain) }
    /// Current physical stress in tensor-shear Voigt order [Pa].
    pub fn stress(&self) -> Voigt { tensor(&self.stress) }
    /// Current recoverable energy density [J/m^3].
    pub fn stored_energy(&self) -> f64 { self.stored }
    /// Total supplied mechanical work density [J/m^3].
    pub fn work(&self) -> f64 { self.work }
    /// Total generated heat density [J/m^3].
    pub fn dissipated(&self) -> f64 { self.heat }
    /// Advance a prescribed strain ramp. dt=0 means an elastic jump with zero
    /// heat. Cancellation/invalid input/overflow publishes no partial state.
    pub fn advance(&mut self, strain: Voigt, dt_s: f64,
        mut cancelled: impl FnMut() -> bool) -> Result<TensorPronyStep> {
        poll(&mut cancelled)?;
        let next_strain = mandel(&strain);
        let norm2 = dot(&next_strain,&next_strain);
        if !dt_s.is_finite() || dt_s < 0.0 || !finite(&next_strain) || !norm2.is_finite()
            || det::sqrt(norm2) > self.model.parameters.strain_limit {
            return Err(invalid("invalid time step or strain outside declared small-strain domain"));
        }
        let delta = core::array::from_fn(|i| next_strain[i]-self.strain[i]);
        let eq = &self.model.parameters.equilibrium_factor;
        let old_eq = apply(eq,&self.strain); let delta_eq = apply(eq,&delta);
        let next_eq = apply(eq,&next_strain);
        let mut stress = [0.0;6]; transpose_add(eq,&next_eq,&mut stress);
        let mut stored = 0.5*dot(&next_eq,&next_eq);
        let mut work = (0..6).fold(0.0, |sum,i| sum+(old_eq[i]+0.5*delta_eq[i])*delta_eq[i]);
        let mut heat = 0.0;
        let mut q = self.q.clone();
        for (q,branch) in q.iter_mut().zip(self.model.branches()) {
            poll(&mut cancelled)?;
            let d = apply(&branch.factor,&delta);
            for i in 0..6 {
                let (next,mean,loss) = ramp(q[i],d[i],dt_s/branch.tau_s);
                q[i]=next; work += mean*d[i]; heat += loss;
            }
            stored += 0.5*dot(q,q); transpose_add(&branch.factor,q,&mut stress);
        }
        let residual = work-(stored-self.stored)-heat;
        let total_work=self.work+work; let total_heat=self.heat+heat;
        if !finite(&stress) || ![stored,work,heat,residual,total_work,total_heat].iter().all(|v|v.is_finite())
            || q.iter().any(|v|!finite(v)) || heat < 0.0 {
            return Err(invalid("nonfinite constitutive response or energy"));
        }
        poll(&mut cancelled)?;
        let report=TensorPronyStep {stress:tensor(&stress),stored_energy:stored,
            work,dissipated:heat,balance_residual:residual};
        self.strain=next_strain; self.stress=stress; self.q=q;
        self.stored=stored; self.work=total_work; self.heat=total_heat;
        Ok(report)
    }
}

// Independent unit-stiffness Maxwell channel in factor coordinates. The same
// analytical ramp used by scalar visco::GeneralizedMaxwell is the test oracle.
// Completing the square keeps heat nonnegative without work-minus-storage
// cancellation: A(q0+h*delta)^2 + R*delta^2.
fn ramp(q0: f64, delta: f64, r: f64) -> (f64,f64,f64) {
    if r == 0.0 { return (q0+delta,q0+0.5*delta,0.0); }
    let a=det::exp(-r); let s=-det::expm1(-r); let b=s/r; let h=b/(1.0+a);
    let (c,remainder)=if r<0.1 {
        let c=0.5+r*(-1.0/6.0+r*(1.0/24.0+r*(-1.0/120.0+r*(1.0/720.0
            +r*(-1.0/5040.0+r*(1.0/40320.0+r*(-1.0/362880.0+r/3628800.0)))))));
        let r2=r*r;
        (c,r*(1.0/12.0+r2*(-1.0/120.0+r2*(17.0/20160.0
            +r2*(-31.0/362880.0+r2*691.0/79833600.0)))))
    } else { ((1.0-b)/r,(1.0-2.0*h)/r) };
    let shifted=q0+h*delta;
    let loss=(0.5*s*(1.0+a))*shifted*shifted+(remainder*delta)*delta;
    (a*q0+b*delta,b*q0+c*delta,loss)
}

#[cfg(test)]
mod tests;
