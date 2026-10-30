//! Small-signal response about an admitted equilibrium of the ORIGINAL storage.
//!
//! Constant J/R/G give E dx_dot = (J-R) W0 dx + G du, with
//! W0=diag(H''(x0)_diff,I_alg). The exact analytic Hessian action is the same
//! callback contract used by StepWorkspace::step_into_analytic_controlled.
//! No finite differences, fitted modal damping, descriptor elimination or
//! replacement nonlinear law. All existing harmonic/pole/contour/residue owners
//! consume the resulting complete linear image.
use super::{HarmonicBudget, HarmonicError, HarmonicModel, MAX_HARMONIC_DIMENSION};
use crate::{DescriptorPortHamiltonian, PortHamiltonian, Storage};

/// A stationary point and its held input in the source's ORIGINAL coordinates.
#[derive(Debug, Clone, Copy)]
pub struct EquilibriumPoint<'a> {
    pub state: &'a [f64],
    pub input: &'a [f64],
    /// Positive reference magnitude for EACH original flow/constraint equation.
    /// The units are those of that equation, NOT state units or a common SI norm.
    pub balance_scales: &'a [f64],
    /// Require max_i |((J-R)e0+G u0)_i|/balance_scales_i <= this value.
    /// Scales are fixed caller data; they are never enlarged after failure.
    pub relative_tolerance: f64,
}

/// Complete-state numerical limits. Preparation includes its extra Hessian and
/// source-observation scratch in harmonic.max_dense_bytes. Callback internals
/// are caller-owned: the action cap is not a wall-time or callback-memory bound.
#[derive(Debug, Clone, Copy)]
pub struct LinearizationBudget {
    pub harmonic: HarmonicBudget,
    /// Exactly n full-state unit-direction actions are requested; no truncation.
    pub max_hessian_actions: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub enum LinearizationError {
    Invalid { what: &'static str },
    NotEquilibrium { row: usize, scaled_residual: f64, tolerance: f64 },
    HessianRefused { column: usize },
    Harmonic(HarmonicError),
    Cancelled,
}
impl core::fmt::Display for LinearizationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "equilibrium linearization refused: {self:?}")
    }
}
impl std::error::Error for LinearizationError {}
impl From<HarmonicError> for LinearizationError {
    fn from(error: HarmonicError) -> Self { Self::Harmonic(error) }
}
fn invalid(what: &'static str) -> LinearizationError { LinearizationError::Invalid { what } }
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), LinearizationError> {
    if cancelled() { Err(LinearizationError::Cancelled) } else { Ok(()) }
}
fn add(sum: &mut f64, correction: &mut f64, term: f64) {
    let y = term - *correction;
    let next = *sum + y;
    *correction = (next - *sum) - y;
    *sum = next;
}

/// Immutable LOCAL image, with operating-point observations kept separate from
/// perturbations. HarmonicResponse.state/output are dx/dy, not x0+dx/y0+dy.
/// The image's quadratic energy is the second variation of H, NOT total energy.
/// A symmetric indefinite Hessian is retained so saddles/unstable poles remain
/// visible; it is NOT admitted as positive storage or repaired to stability.
/// Validity for finite oscillation amplitude must be established by the source.
/// This is not a periodic-orbit/Floquet, contact-switching or limit-cycle model.
pub struct EquilibriumLinearization {
    model: HarmonicModel,
    state: Vec<f64>,
    input: Vec<f64>,
    effort: Vec<f64>,
    output: Vec<f64>,
    balance: Vec<f64>,
    maximum_scaled_residual: f64,
    stored_energy_j: f64,
}
impl EquilibriumLinearization {
    /// Linearize a complete nonlinear ODE at a caller-resolved stationary point.
    /// `hessian` must be the exact, read-only H''(state)*direction of THIS source,
    /// including every coupled/internal coordinate. False or incomplete output
    /// refuses. At a nonsmooth point the caller must refuse unless a genuine
    /// two-sided derivative exists; an arbitrary active-branch choice is not one.
    pub fn from_ode(
        source: &PortHamiltonian, point: EquilibriumPoint<'_>, budget: LinearizationBudget,
        hessian: &dyn Fn(&[f64], &[f64], &mut [f64]) -> bool,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, LinearizationError> {
        Self::prepare(source.n, source.n, source.m, &source.j, &source.r, &source.g,
            source.storage.as_ref(), point, budget, hessian, &mut cancelled)
    }

    /// Same linearization with the original algebraic efforts retained. The
    /// callback sees FULL state vectors and must return zero Hessian rows and
    /// columns for multipliers. No inverse of E and no constraint regularization.
    /// The source storage wrapper evaluates only its differential coordinates.
    pub fn from_descriptor(
        source: &DescriptorPortHamiltonian, point: EquilibriumPoint<'_>, budget: LinearizationBudget,
        hessian: &dyn Fn(&[f64], &[f64], &mut [f64]) -> bool,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, LinearizationError> {
        Self::prepare(source.n, source.n_diff, source.m, &source.j, &source.r, &source.g,
            source.storage.as_ref(), point, budget, hessian, &mut cancelled)
    }

    /// All existing AC, isolated-pole, contour and participation APIs, in LOCAL
    /// perturbation coordinates. No discretized time-step frequencies are used.
    pub fn model(&self) -> &HarmonicModel { &self.model }
    pub fn equilibrium_state(&self) -> &[f64] { &self.state }
    pub fn equilibrium_input(&self) -> &[f64] { &self.input }
    pub fn equilibrium_effort(&self) -> &[f64] { &self.effort }
    pub fn equilibrium_output(&self) -> &[f64] { &self.output }
    pub fn equilibrium_balance(&self) -> &[f64] { &self.balance }
    pub fn maximum_scaled_residual(&self) -> f64 { self.maximum_scaled_residual }
    pub fn equilibrium_stored_energy_j(&self) -> f64 { self.stored_energy_j }
    /// Full original-coordinate Hessian, including zero algebraic rows/columns.
    pub fn tangent_storage_matrix(&self) -> &[f64] { &self.model.q }

    #[allow(clippy::too_many_arguments)]
    fn prepare(
        n: usize, nd: usize, m: usize, j: &[f64], r: &[f64], g: &[f64], storage: &dyn Storage,
        point: EquilibriumPoint<'_>, budget: LinearizationBudget,
        hessian: &dyn Fn(&[f64], &[f64], &mut [f64]) -> bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, LinearizationError> {
        poll(cancelled)?;
        let h = budget.harmonic;
        // Use wide integer arithmetic only after the dense dimension hard cap.
        if n == 0 || nd > n || n > MAX_HARMONIC_DIMENSION || n > h.max_dimension
            || n > budget.max_hessian_actions {
            return Err(invalid("complete state exceeds the dimension or Hessian-action budget"));
        }
        let (wide_n, wide_m) = (n as u128, m as u128);
        let bytes = 112*wide_n*wide_n + 16*wide_n*wide_m + 512*wide_n + 128*wide_m;
        let work = wide_n*wide_n*(wide_n+wide_m+16);
        if bytes > h.max_dense_bytes as u128 || work > h.max_dense_work as u128 {
            return Err(invalid("linearization exceeds complete preparation payload/work budget"));
        }
        let ports = n.checked_mul(m).ok_or_else(|| invalid("port dimension overflow"))?;
        if j.len()!=n*n || r.len()!=n*n || g.len()!=ports
            || j.iter().chain(r).chain(g).any(|x|!x.is_finite())
            || point.state.len()!=n || point.input.len()!=m || point.balance_scales.len()!=n
            || point.state.iter().chain(point.input).any(|x|!x.is_finite())
            || point.balance_scales.iter().any(|x|!x.is_finite() || *x<=0.0)
            || !point.relative_tolerance.is_finite() || point.relative_tolerance<=0.0
            || point.relative_tolerance>=1.0 || !h.relative_residual_tolerance.is_finite()
            || h.relative_residual_tolerance<=0.0 || h.relative_residual_tolerance>=1.0 {
            return Err(invalid("finite complete source/point data and positive accuracy/scales required"));
        }
        let energy = storage.hamiltonian(point.state);
        let mut effort = vec![f64::NAN;n];
        storage.gradient(point.state,&mut effort);
        poll(cancelled)?;
        if !energy.is_finite() || effort.iter().any(|x|!x.is_finite()) {
            return Err(invalid("original storage refused the operating point"));
        }
        if effort[nd..].iter().any(|&x|x!=0.0) {
            return Err(invalid("algebraic coordinates cannot carry stored effort"));
        }
        effort[nd..].copy_from_slice(&point.state[nd..]);
        let mut balance = vec![0.0;n];
        let mut maximum = 0.0_f64;
        for i in 0..n {
            poll(cancelled)?;
            let mut correction = 0.0;
            for k in 0..n { add(&mut balance[i],&mut correction,(j[i*n+k]-r[i*n+k])*effort[k]); }
            for p in 0..m { add(&mut balance[i],&mut correction,g[i*m+p]*point.input[p]); }
            let scaled = balance[i].abs()/point.balance_scales[i];
            if !balance[i].is_finite() || !scaled.is_finite() {
                return Err(invalid("original equilibrium balance is not representable"));
            }
            if scaled>point.relative_tolerance {
                return Err(LinearizationError::NotEquilibrium {
                    row:i,scaled_residual:scaled,tolerance:point.relative_tolerance,
                });
            }
            maximum=maximum.max(scaled);
        }
        let mut matrix=vec![0.0;n*n];
        let mut direction=vec![0.0;n];
        let mut column=vec![f64::NAN;n];
        for k in 0..n {
            poll(cancelled)?;
            direction[k]=1.0;
            column.fill(f64::NAN);
            if !hessian(point.state,&direction,&mut column) || column.iter().any(|x|!x.is_finite()) {
                return Err(LinearizationError::HessianRefused { column:k });
            }
            direction[k]=0.0;
            for i in 0..n { matrix[i*n+k]=column[i]; }
        }
        poll(cancelled)?;
        crate::require_symmetric(&matrix,n,"operating-point storage Hessian")
            .map_err(|_|invalid("analytic storage Hessian must be symmetric; no repair is performed"))?;
        for i in 0..n { for k in 0..n {
            if (i>=nd || k>=nd) && matrix[i*n+k]!=0.0 {
                return Err(invalid("storage Hessian cannot couple algebraic coordinates"));
            }
        }}
        let mut output=vec![0.0;m];
        for p in 0..m {
            poll(cancelled)?;
            let mut correction=0.0;
            for i in 0..n { add(&mut output[p],&mut correction,g[i*m+p]*effort[i]); }
        }
        if output.iter().any(|x|!x.is_finite()) { return Err(invalid("equilibrium port output overflow")); }
        // This private transient adapter describes ONLY the local second variation.
        // It never replaces the source storage and is never exposed for stepping.
        let tangent=LocalQuadratic { n, matrix };
        let model=HarmonicModel::prepare(n,nd,m,j,r,g,&tangent,h,cancelled)?;
        poll(cancelled)?;
        Ok(Self { model,state:point.state.to_vec(),input:point.input.to_vec(),effort,output,balance,
            maximum_scaled_residual:maximum,stored_energy_j:energy })
    }
}

struct LocalQuadratic { n:usize, matrix:Vec<f64> }
impl Storage for LocalQuadratic {
    fn quadratic_part(&self)->Option<crate::linear_exact::QuadraticPart> {
        Some(crate::linear_exact::QuadraticPart { matrix:self.matrix.clone(),exact:true })
    }
    fn hamiltonian(&self,x:&[f64])->f64 {
        0.5*(0..self.n).map(|i|x[i]*(0..self.n)
            .map(|k|self.matrix[i*self.n+k]*x[k]).sum::<f64>()).sum::<f64>()
    }
    fn gradient(&self,x:&[f64],out:&mut[f64]) {
        for i in 0..self.n { out[i]=(0..self.n).map(|k|self.matrix[i*self.n+k]*x[k]).sum(); }
    }
}
