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

mod projection;
pub use projection::{ProjectedPronyBranch, ProjectedTensorProny, ProjectedTensorTemperatureDerivative,
    TensorProjectionBudget, TensorThermalExpansion, TensorThermalState};

mod reduction;
pub use reduction::{Complex3, HarmonicPlateSection, PlaneStressModulus};

/// Bind complete hereditary tensors to original material-card state points.
pub mod source;

/// Exact equal-time coordinates with original material-history reconstruction.
pub mod coordinates;

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

pub use crate::visco::tensor::{PronyBranch, PronyBranch as TensorPronyBranch, TensorProny};

impl PronyBranch {
    /// Construct any PSD tensor, including zero/rank-deficient tensors, as B^T B.
    pub fn from_gram_rows(factor: Tangent6, tau_s: f64) -> Result<Self> {
        if !tau_s.is_finite() || tau_s <= 0.0 {
            return Err(invalid("relaxation time must be finite and positive"));
        }
        let stiffness = gram(&factor)?;
        Ok(Self {
            c_pa: stiffness,
            tau_s,
            gram_rows: Some(factor),
        })
    }
    /// Convenience for a strictly SPD Mandel stiffness; no repair is performed.
    pub fn from_stiffness_mandel(stiffness: Tangent6, tau_s: f64) -> Result<Self> {
        Self::from_gram_rows(positive_factor(&stiffness)?, tau_s)
    }
    /// Actual factor-defined Mandel stiffness [Pa].
    pub fn stiffness_mandel(&self) -> &Tangent6 {
        &self.c_pa
    }
    /// Relaxation time [s].
    pub fn relaxation_time_s(&self) -> f64 {
        self.tau_s
    }
}

impl TensorProny {
    /// Admit an SPD equilibrium tensor and up to MAX_BRANCHES PSD branches.
    /// `strain_limit` bounds the Frobenius norm of the small-strain tensor.
    /// Factorization changes the represented equilibrium by rounding only.
    pub fn new_mandel(
        equilibrium: Tangent6,
        branches: Vec<TensorPronyBranch>,
        strain_limit: f64,
    ) -> Result<Self> {
        if branches.len() > MAX_BRANCHES || !strain_limit.is_finite() || strain_limit <= 0.0 {
            return Err(invalid(
                "invalid branch count or positive finite strain bound",
            ));
        }
        let equilibrium_gram = positive_factor(&equilibrium)?;
        let c_eq = gram(&equilibrium_gram)?;
        let c_eq_factor =
            core::array::from_fn(|i| core::array::from_fn(|j| equilibrium_gram[j][i]));
        // The instantaneous tensor must be representable too.
        let mut instantaneous = c_eq;
        let mut factors = Vec::with_capacity(branches.len());
        let mut resolved_branches = Vec::with_capacity(branches.len());
        for branch in branches {
            if !branch.tau_s.is_finite() || branch.tau_s <= 0.0 {
                return Err(invalid("relaxation time must be finite and positive"));
            }
            let gram_rows = match branch.gram_rows {
                Some(g) => g,
                None => {
                    let l = crate::visco::tensor::psd_factor(&branch.c_pa)
                        .map_err(|_| invalid("branch stiffness must be positive semidefinite"))?;
                    core::array::from_fn(|i| core::array::from_fn(|j| l[j][i]))
                }
            };
            let stiffness = gram(&gram_rows)?;
            let factor = core::array::from_fn(|i| core::array::from_fn(|j| gram_rows[j][i]));
            for i in 0..6 {
                for j in 0..6 {
                    instantaneous[i][j] += stiffness[i][j];
                    if !instantaneous[i][j].is_finite() {
                        return Err(invalid("instantaneous stiffness overflow"));
                    }
                }
            }
            factors.push(factor);
            resolved_branches.push(PronyBranch {
                c_pa: stiffness,
                tau_s: branch.tau_s,
                gram_rows: Some(gram_rows),
            });
        }
        let branch_tags = crate::visco::tensor::default_branch_tags(resolved_branches.len());
        Ok(Self {
            c_eq,
            c_eq_factor,
            equilibrium_gram,
            kv_pa_s: None,
            kv_factor: None,
            branches: resolved_branches,
            factors,
            band_rad_s: (0.0, f64::MAX),
            strain_limit,
            iso_k3_eq_pa: None,
            branch_tags,
        })
    }
    /// Convert an existing tensor-shear stress/strain tangent (e.g. elastic::
    /// IsotropicElastic::stiffness) without doubling shear stiffness.
    pub fn from_tensor_shear(
        equilibrium: Tangent6,
        branches: Vec<TensorPronyBranch>,
        strain_limit: f64,
    ) -> Result<Self> {
        let c = core::array::from_fn(|i| {
            core::array::from_fn(|j| (WEIGHTS[i] / WEIGHTS[j]) * equilibrium[i][j])
        });
        Self::new_mandel(c, branches, strain_limit)
    }
    /// Actual equilibrium Mandel stiffness [Pa].
    pub fn equilibrium_mandel(&self) -> &Tangent6 {
        &self.c_eq
    }
    /// Original Frobenius-norm ceiling on the small-strain tensor.
    /// Infinite when the frequency-domain constructor supplied no strain cap.
    pub const fn maximum_strain_norm(&self) -> f64 { self.strain_limit }

    /// Virgin zero-strain, relaxed history. No state can be rebound to another law.
    pub fn history(&self) -> TensorPronyHistory {
        TensorPronyHistory {
            model: self.clone(),
            strain: [0.0; 6],
            stress: [0.0; 6],
            q: vec![[0.0; 6]; self.branches().len()],
            stored: 0.0,
            work: 0.0,
            heat: 0.0,
        }
    }
    /// Full storage/loss tensors for exp(+i omega t), matching visco::modulus.
    /// The loss matrix is PSD; its OFF-DIAGONAL entries need not be positive.
    pub fn modulus_mandel(&self, omega: f64) -> Result<(Tangent6, Tangent6)> {
        if !omega.is_finite() || omega < 0.0 {
            return Err(invalid("frequency must be finite and nonnegative"));
        }
        let mut storage = self.c_eq;
        let mut loss = [[0.0; 6]; 6];
        if let Some(kv) = &self.kv_pa_s {
            for i in 0..6 {
                for j in 0..6 {
                    loss[i][j] += omega * kv[i][j];
                }
            }
        }
        for branch in self.branches() {
            let z = omega * branch.tau_s;
            // Stable at both ends, including overflow of positive omega*tau.
            let (re, im) = if z <= 1.0 {
                (z * (z / (1.0 + z * z)), z / (1.0 + z * z))
            } else {
                let u = 1.0 / z;
                (1.0 / (1.0 + u * u), u / (1.0 + u * u))
            };
            for i in 0..6 {
                for j in 0..6 {
                    storage[i][j] += re * branch.c_pa[i][j];
                    loss[i][j] += im * branch.c_pa[i][j];
                }
            }
        }
        if storage
            .iter()
            .flatten()
            .chain(loss.iter().flatten())
            .any(|v| !v.is_finite())
        {
            return Err(invalid("modulus overflow"));
        }
        Ok((storage, loss))
    }
    /// Exact endpoint stress derivative for a prescribed linear strain ramp.
    /// Returned tensor uses the crate's tensor-shear Voigt convention.
    /// A nonzero Kelvin viscosity requires positive dt; its tangent is C_kv/dt.
    pub fn tangent(&self, dt_s: f64) -> Result<Tangent6> {
        if !dt_s.is_finite() || dt_s < 0.0 {
            return Err(invalid("time step must be finite and nonnegative"));
        }
        let mut c = self.c_eq;
        for branch in self.branches() {
            let r = dt_s / branch.tau_s;
            let b = if r == 0.0 { 1.0 } else { -det::expm1(-r) / r };
            for i in 0..6 {
                for j in 0..6 {
                    c[i][j] += b * branch.c_pa[i][j];
                }
            }
        }
        if let Some(kv) = self.kv_pa_s().filter(|c| c.iter().flatten().any(|v| *v != 0.0)) {
            if dt_s == 0.0 { return Err(invalid("Kelvin viscosity has no finite instantaneous strain tangent")); }
            for i in 0..6 { for j in 0..6 { c[i][j] += kv[i][j]/dt_s; } }
        }
        for i in 0..6 {
            for j in 0..6 {
                c[i][j] *= WEIGHTS[j] / WEIGHTS[i];
            }
        }
        if c.iter().flatten().any(|v| !v.is_finite()) {
            return Err(invalid("tangent overflow"));
        }
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
    /// Advance a prescribed strain ramp. Without Kelvin viscosity, dt=0 is an
    /// elastic jump. Nonzero Kelvin viscosity refuses instantaneous strain
    /// changes; positive-duration ramps include its stress and dashpot heat.
    /// Cancellation/invalid input/overflow publishes no partial state.
    pub fn advance(&mut self, strain: Voigt, dt_s: f64,
        mut cancelled: impl FnMut() -> bool) -> Result<TensorPronyStep> {
        poll(&mut cancelled)?;
        let next_strain = mandel(&strain);
        if !dt_s.is_finite() || dt_s < 0.0 {
            return Err(invalid("invalid time step or strain outside declared small-strain domain"));
        }
        self.model.validate_mandel_strain(&next_strain)
            .map_err(|_| invalid("invalid time step or strain outside declared small-strain domain"))?;
        let delta = core::array::from_fn(|i| next_strain[i]-self.strain[i]);
        let eq = &self.model.equilibrium_gram;
        let old_eq = apply(eq,&self.strain); let delta_eq = apply(eq,&delta);
        let next_eq = apply(eq,&next_strain);
        let mut stress = [0.0;6]; transpose_add(eq,&next_eq,&mut stress);
        let mut stored = 0.5*dot(&next_eq,&next_eq);
        let mut work = (0..6).fold(0.0, |sum,i| sum+(old_eq[i]+0.5*delta_eq[i])*delta_eq[i]);
        let mut heat = 0.0;
        let mut q = self.q.clone();
        for ((q, branch), factor) in q
            .iter_mut()
            .zip(self.model.branches())
            .zip(&self.model.factors)
        {
            poll(&mut cancelled)?;
            let gram_rows = branch.resolved_gram_rows(factor);
            let d = apply(&gram_rows, &delta);
            let ramp = crate::visco::PronyRamp::new(dt_s / branch.tau_s);
            for i in 0..6 {
                let (next, mean, loss) = ramp.apply(q[i], d[i]);
                q[i]=next; work += mean*d[i]; heat += loss;
            }
            stored += 0.5*dot(q,q); transpose_add(&gram_rows,q,&mut stress);
        }
        if let Some(kv) = self.model.kv_pa_s().filter(|c| c.iter().flatten().any(|v| *v != 0.0)) {
            if dt_s == 0.0 {
                if delta.iter().any(|v| *v != 0.0) {
                    return Err(invalid("Kelvin viscosity requires positive time for a strain change"));
                }
            } else {
                let rate = delta.map(|v| v/dt_s);
                let viscous = apply(kv,&rate);
                let factor = self.model.kv_factor.as_ref()
                    .ok_or_else(||invalid("missing admitted Kelvin viscosity factor"))?;
                // Independent nonnegative integral dt * |L^T strain_rate|^2.
                // Kelvin viscosity stores no energy and has no memory state.
                let mut kelvin_heat = 0.0;
                for i in 0..6 {
                    stress[i] += viscous[i];
                    let channel = (0..6).fold(0.0,|s,k|s+factor[k][i]*rate[k]);
                    kelvin_heat += (dt_s*channel)*channel;
                }
                work += dot(&viscous,&delta);
                heat += kelvin_heat;
            }
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

#[cfg(test)]
fn ramp(q0: f64, delta: f64, r: f64) -> (f64, f64, f64) {
    crate::visco::prony_ramp(q0, delta, r)
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod kelvin_tests;
