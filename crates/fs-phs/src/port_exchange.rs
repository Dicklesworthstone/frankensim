//! Prepared collective lossless exchange of two unit-energy coordinate banks.
//!
//! For `left' = -L^T right`, `right' = L left`, thin SVD supplies independent
//! commuting rotations. The ENTIRE coupling is evolved together, not a sequence
//! of physical-port rotations whose splitting error depends on an arbitrary
//! basis/order. Reuse fs-la's SVD and fs-math's deterministic trigonometry.
//!
//! This is the exact isolated coupling flow to the admitted decomposition and
//! floating-point error, NOT an exact propagator for a surrounding split system.
//! Physical modal frequencies, damping, displacement and time belong to callers.
use crate::{PhsError, PreparedStepError};
use fs_la::factor::svd_jacobi;
use fs_math::det;

/// Explicit cold-work and rate-resolution admission, not a physical cutoff.
#[derive(Clone, Copy, Debug)]
pub struct PortExchangeBudget {
    /// Maximum number of coordinates in the left bank.
    pub max_left: usize,
    /// Maximum number of coordinates in the right bank.
    pub max_right: usize,
    /// Bound on `60 * max(left,right) * min(left,right)^2` SVD sweep terms.
    pub max_setup_terms: usize,
    /// Bound on `dt * ||[0,-L^T;L,0]||_infinity` (positive and at most one).
    /// Accuracy of any enclosing operator splitting remains a separate matter.
    pub maximum_dt_coupling: f64,
}

/// Signed floating-point defect of a conservative exchange, NEVER dissipation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PortExchangeRecord {
    /// Change in the combined quadratic energy, in caller-supplied energy units.
    pub energy_defect: f64,
}

/// Cold-prepared full coupling with allocation-free, transactional application.
///
/// Coordinate storage is `H=(||left||^2+||right||^2)/2`. A physical caller must
/// mass/energy normalize its velocities first. All original coordinates remain;
/// neither a small singular value nor the uncoupled complement is truncated.
/// The SVD is checked for orthogonality and componentwise reconstruction BEFORE
/// publication. Work in the trusted cold SVD is bounded by `max_setup_terms`;
/// its internal sweeps are not interruptible. Application polls between rows.
#[derive(Debug)]
pub struct PreparedPortExchange {
    left: usize,
    right: usize,
    directions: usize,
    left_basis: Vec<f64>,
    right_basis: Vec<f64>,
    sine: Vec<f64>,
    cosine_minus_one: Vec<f64>,
    delta_left: Vec<f64>,
    delta_right: Vec<f64>,
    candidate_left: Vec<f64>,
    candidate_right: Vec<f64>,
}

fn bad(what: &'static str) -> PhsError { PhsError::Dimension { what } }
fn zeroed(n: usize) -> Result<Vec<f64>, PhsError> {
    let mut v=Vec::new();
    v.try_reserve_exact(n).map_err(|_|bad("port exchange capacity"))?;
    v.resize(n,0.0); Ok(v)
}
fn cancelled<F: FnMut()->bool>(poll: &mut F) -> Result<(), PreparedStepError> {
    if poll() { Err(PreparedStepError::Cancelled) } else { Ok(()) }
}

impl PreparedPortExchange {
    /// Prepare one duration of collective power exchange. `coupling` is the
    /// row-major RIGHT-by-LEFT matrix L, in inverse time units. Negative entries
    /// are physical signed ports. A zero matrix is an exact identity operation.
    ///
    /// # Errors
    /// Bad dimensions, budgets, finite/rate bounds, capacity, or unresolved SVD.
    /// Failed reconstruction never causes a rank cutoff or pairwise fallback.
    pub fn new(left: usize, right: usize, coupling: &[f64], dt: f64,
        budget: PortExchangeBudget) -> Result<Self, PhsError>
    {
        let rank=left.min(right);let tall=left.max(right);
        let entries=left.checked_mul(right).ok_or_else(||bad("port exchange extent"))?;
        let work=rank.checked_mul(rank).and_then(|x|x.checked_mul(tall))
            .and_then(|x|x.checked_mul(60)).ok_or_else(||bad("port exchange setup extent"))?;
        if left==0 || right==0 || left>budget.max_left || right>budget.max_right
            || coupling.len()!=entries || work>budget.max_setup_terms
            || !dt.is_finite() || dt<=0.0
            || !budget.maximum_dt_coupling.is_finite()
            || budget.maximum_dt_coupling<=0.0 || budget.maximum_dt_coupling>1.0
            || coupling.iter().any(|x|!x.is_finite()) {
            return Err(bad("port exchange dimensions, finite inputs or work budget"));
        }
        let mut columns=zeroed(left)?;let mut maximum=0.0_f64;let mut scale=0.0_f64;
        for row in coupling.chunks_exact(left) {
            let mut sum=0.0;
            for (j,&x) in row.iter().enumerate() {sum+=x.abs();columns[j]+=x.abs();scale=scale.max(x.abs());}
            maximum=maximum.max(sum);
        }
        maximum=maximum.max(columns.into_iter().fold(0.0_f64,f64::max));
        if !maximum.is_finite() || !((maximum*dt).is_finite())
            || maximum*dt>budget.maximum_dt_coupling {
            return Err(bad("port exchange exceeds the declared rate-resolution bound"));
        }
        let mut result=Self {left,right,directions:0,left_basis:Vec::new(),right_basis:Vec::new(),
            sine:Vec::new(),cosine_minus_one:Vec::new(),delta_left:Vec::new(),delta_right:Vec::new(),
            candidate_left:zeroed(left)?,candidate_right:zeroed(right)?};
        if scale==0.0 {return Ok(result);}
        // Normalize before the existing SVD so its norm squares cannot overflow.
        // Use the transpose for a wide matrix; no padded fictitious channels.
        let mut normalized=zeroed(entries)?;
        for i in 0..right {for j in 0..left {
            let value=coupling[i*left+j]/scale;
            if coupling[i*left+j]!=0.0 && value==0.0 {
                return Err(bad("port exchange coefficient underflows SVD scaling"));
            }
            normalized[if right>=left {i*left+j}else{j*right+i}]=value;
        }}
        let svd=svd_jacobi(&normalized,tall,rank);
        if svd.sigma.iter().any(|s|!s.is_finite() || *s<0.0)
            || svd.u.iter().chain(&svd.v).any(|x|!x.is_finite()) {
            return Err(bad("port exchange SVD has nonfinite factors"));
        }
        let (lb,rb)=if right>=left {(svd.v,svd.u)} else {(svd.u,svd.v)};
        let tolerance=256.0*f64::EPSILON*(tall+rank+1) as f64;
        for (basis,rows) in [(&lb,left),(&rb,right)] {
            for i in 0..rank {for j in 0..=i {
                // A zero singular column has no active rotation. Its U column
                // may be exactly zero; do not invent a vector for the nullspace.
                if svd.sigma[i]==0.0 || svd.sigma[j]==0.0 {continue;}
                let dot=(0..rows).map(|r|basis[r*rank+i]*basis[r*rank+j]).sum::<f64>();
                if !dot.is_finite() || (dot-if i==j {1.0}else{0.0}).abs()>tolerance {
                    return Err(bad("port exchange SVD directions are not orthonormal"));
                }
            }}
        }
        for i in 0..right {for j in 0..left {
            let wanted=coupling[i*left+j]/scale;let mut value=0.0;let mut magnitude=wanted.abs();
            for k in 0..rank {
                let term=rb[i*rank+k]*svd.sigma[k]*lb[j*rank+k];value+=term;magnitude+=term.abs();
            }
            if !value.is_finite() || !magnitude.is_finite() || (value-wanted).abs()>tolerance*magnitude {
                return Err(bad("port exchange SVD fails original componentwise reconstruction"));
            }
        }}
        let duration=scale*dt;
        result.sine=zeroed(rank)?;result.cosine_minus_one=zeroed(rank)?;
        for (k,&s) in svd.sigma.iter().enumerate() {
            let angle=s*duration;
            if !angle.is_finite() || s!=0.0 && angle==0.0 {
                return Err(bad("port exchange rotation is not representable"));
            }
            result.sine[k]=det::sin(angle);
            // Preserve small rotations without subtracting two numbers near 1.
            result.cosine_minus_one[k]=-2.0*det::sin(0.5*angle).powi(2);
        }
        result.left_basis=lb;result.right_basis=rb;result.directions=rank;
        result.delta_left=zeroed(rank)?;result.delta_right=zeroed(rank)?;
        Ok(result)
    }

    /// Apply the prepared isolated coupling flow in place. All state is published
    /// together, after finite/energy checks. Scratch is private; no allocations.
    ///
    /// # Errors
    /// Invalid state shape, nonfinite energy or failed conservative roundoff gate.
    pub fn apply(&mut self, left: &mut [f64], right: &mut [f64])
        -> Result<PortExchangeRecord, PreparedStepError>
    {
        self.apply_controlled(left,right,||false)
    }

    /// Same operation with bounded row-boundary cancellation and exact retry.
    /// No failed/cancelled call changes either caller-owned bank.
    ///
    /// # Errors
    /// The same refusals as [`Self::apply`], plus cancellation.
    pub fn apply_controlled<F: FnMut()->bool>(&mut self, left: &mut [f64], right: &mut [f64],
        mut poll: F) -> Result<PortExchangeRecord, PreparedStepError>
    {
        cancelled(&mut poll)?;
        if left.len()!=self.left || right.len()!=self.right
            || left.iter().chain(right.iter()).any(|x|!x.is_finite()) {
            return Err(bad("port exchange state dimensions or finiteness").into());
        }
        let before=left.iter().chain(right.iter()).map(|x|0.5*x*x).sum::<f64>();
        if !before.is_finite() {return Err(bad("port exchange energy overflow").into());}
        if self.directions==0 {
            cancelled(&mut poll)?;
            return Ok(PortExchangeRecord {energy_defect:0.0});
        }
        let rank=self.directions;
        self.delta_left.fill(0.0);self.delta_right.fill(0.0);
        for (r,&x) in left.iter().enumerate() {
            cancelled(&mut poll)?;
            for k in 0..rank {self.delta_left[k]+=self.left_basis[r*rank+k]*x;}
        }
        for (r,&x) in right.iter().enumerate() {
            cancelled(&mut poll)?;
            for k in 0..rank {self.delta_right[k]+=self.right_basis[r*rank+k]*x;}
        }
        for k in 0..rank {
            let a=self.delta_left[k];let b=self.delta_right[k];
            self.delta_left[k]=self.cosine_minus_one[k]*a-self.sine[k]*b;
            self.delta_right[k]=self.sine[k]*a+self.cosine_minus_one[k]*b;
        }
        // Update the ORIGINAL vectors. Reconstructing a tall bank from thin U
        // would discard its uncoupled complement, including acoustic history.
        for (r,&x) in left.iter().enumerate() {
            cancelled(&mut poll)?;
            let mut delta=0.0;for k in 0..rank {delta+=self.left_basis[r*rank+k]*self.delta_left[k];}
            self.candidate_left[r]=x+delta;
        }
        for (r,&x) in right.iter().enumerate() {
            cancelled(&mut poll)?;
            let mut delta=0.0;for k in 0..rank {delta+=self.right_basis[r*rank+k]*self.delta_right[k];}
            self.candidate_right[r]=x+delta;
        }
        let after=self.candidate_left.iter().chain(&self.candidate_right).map(|x|0.5*x*x).sum::<f64>();
        let defect=after-before;
        let tolerance=512.0*f64::EPSILON*(self.left+self.right+rank+1) as f64*(before+after);
        if !after.is_finite() || !tolerance.is_finite() || defect.abs()>tolerance {
            return Err(bad("collective port exchange fails conservative energy admission").into());
        }
        cancelled(&mut poll)?;
        left.copy_from_slice(&self.candidate_left);right.copy_from_slice(&self.candidate_right);
        Ok(PortExchangeRecord {energy_defect:defect})
    }
}

#[cfg(test)]
#[path="port_exchange_tests.rs"]
mod tests;
