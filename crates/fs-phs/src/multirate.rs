//! Physical subcycling with work-conjugate held interface inputs.
//!
//! Each complete subsystem advances on its own clock. The interface solves
//! against ACTUAL window-average outputs, making its two works equal/opposite.
//! This is a second-order held-interface discretization for smooth converged
//! systems, not the exact flow of their interconnection. Refine the macro window
//! as well as local clocks: subcycling alone cannot remove interface phase error.
use crate::linear_exact::{LinearSplitBudget, PreparedExactZoh, energetic_scale};
use crate::{PhsError, PortHamiltonian, PreparedStepError, StepWorkspace, storage_energy_change};
use fs_la::factor::{Lu, lu};

/// Maximum native physical substeps per common window.
pub const MAX_SUBSTEPS: usize = 16;
/// Maximum simultaneous power-port pairs in the dense interface solve.
pub const MAX_INTERFACE_PORTS: usize = 32;

/// Explicit choice of the original subsystem's time integrator.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum IntegrationMethod {
    /// Complete exponential dynamics and independent time-integrated heat.
    ExactQuadratic,
    /// Original Gonzalez method with a caller-declared bound on all local pole
    /// rates, including nonlinear tangent stiffness in the operating domain.
    DiscreteGradient { maximum_rate_rad_s: f64, maximum_rate_dt: f64 },
}

/// Native subsystem clock and an optional independent contact-resolution bound.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Clock {
    pub substeps: usize,
    pub method: IntegrationMethod,
    /// Optional T_contact/N [s]; also applies to exact linear local stepping.
    pub maximum_contact_step_s: Option<f64>,
}

/// Canonical skew connection: u_left = -average(y_right), u_right = average(y_left).
/// The caller supplies compatible conjugate units. Physical transformer/area/
/// gear ratios belong to the original owners, not to numerical input scaling.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PortPair {
    pub left: usize,
    pub right: usize,
    /// Positive numerical reference with units of left input/right output.
    pub input_scale: f64,
}

/// Explicit bounded interface solve and independent physical energy admission.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SolveConfig {
    /// At most 64 Newton updates; zero admits only an already solved interface.
    pub max_iterations: usize,
    pub relative_tolerance: f64,
    /// Absolute energy allowance [J], never raised silently.
    pub absolute_energy_tolerance_j: f64,
    /// Relative to a conservative coordinate-dependent numerical norm of
    /// x_i*grad(H)_i, actual energy changes, work and heat; invariant under
    /// adding a constant energy datum. That norm is not thermodynamic storage.
    pub relative_energy_tolerance: f64,
}
impl Default for SolveConfig {
    fn default() -> Self { Self { max_iterations: 16, relative_tolerance: 1e-13,
        absolute_energy_tolerance_j: 0.0, relative_energy_tolerance: 1e-12 } }
}

/// Every refusal leaves caller-owned subsystem states unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum MultirateError {
    Invalid(&'static str),
    Subsystem { side: &'static str, error: PhsError },
    NotConverged { iterations: usize, residual: f64 },
    Energy { side: &'static str, residual_j: f64, tolerance_j: f64 },
    Cancelled,
}
impl core::fmt::Display for MultirateError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "multirate coupling: {self:?}")
    }
}
impl std::error::Error for MultirateError {}

/// Complete original-state evolution over one physical window.
#[derive(Clone, Debug, PartialEq)]
pub struct SubsystemStep {
    pub state: Vec<f64>,
    pub average_output: Vec<f64>,
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    /// Compensated sum of original per-step constitutive energy changes.
    pub stored_change_j: f64,
    /// Independently evaluated physical heat, never a balance remainder.
    pub dissipated_j: f64,
    /// Complete native input work, including the interface contribution.
    pub supplied_j: f64,
    pub substeps: usize,
}

/// Accepted states, independent source work, interface exchange, and heat.
#[derive(Clone, Debug, PartialEq)]
pub struct CoupledStep {
    pub left: SubsystemStep,
    pub right: SubsystemStep,
    pub external_work_j: [f64; 2],
    /// Signed work entering each subsystem through its connection.
    pub exchange_work_j: [f64; 2],
    pub balance_residual_j: f64,
    pub interface_residual: f64,
    pub iterations: usize,
}

#[derive(Default, Clone, Copy)]
struct Sum { high: f64, low: f64 }
impl Sum {
    fn add(&mut self, value: f64) {
        let next = self.high + value;
        self.low += if self.high.abs() >= value.abs() {
            (self.high-next)+value
        } else { (value-next)+self.high };
        self.high = next;
    }
    fn value(self) -> f64 { self.high+self.low }
}
fn sum(values: impl IntoIterator<Item=f64>) -> f64 {
    let mut total = Sum::default();
    for value in values { total.add(value); }
    total.value()
}
fn finite(values: &[f64]) -> bool { values.iter().all(|x| x.is_finite()) }
fn positive(x: f64) -> bool { x.is_finite() && x > 0.0 }
fn nonnegative(x: f64) -> bool { x.is_finite() && x >= 0.0 }
fn check(cancelled: &mut impl FnMut() -> bool) -> Result<(), MultirateError> {
    if cancelled() { Err(MultirateError::Cancelled) } else { Ok(()) }
}
fn native(side: &'static str, error: PreparedStepError) -> MultirateError {
    match error { PreparedStepError::Cancelled => MultirateError::Cancelled,
        PreparedStepError::Solver(error) => MultirateError::Subsystem { side, error } }
}

enum Evolution<'a> { Exact(PreparedExactZoh<'a>), Gonzalez(StepWorkspace) }
struct LocalClock<'a> {
    system: &'a PortHamiltonian,
    clock: Clock,
    dt: f64,
    evolution: Evolution<'a>,
    side: &'static str,
}
impl<'a> LocalClock<'a> {
    fn new(system: &'a PortHamiltonian, clock: Clock, window: f64,
        side: &'static str) -> Result<Self, MultirateError> {
        if !(1..=MAX_SUBSTEPS).contains(&clock.substeps) {
            return Err(MultirateError::Invalid("physical substeps must lie in 1..=16"));
        }
        let dt = window/clock.substeps as f64;
        if !positive(dt) || clock.maximum_contact_step_s.is_some_and(|bound|
            !positive(bound) || dt > bound) {
            return Err(MultirateError::Invalid("invalid or unresolved physical contact clock"));
        }
        let evolution = match clock.method {
            IntegrationMethod::ExactQuadratic => Evolution::Exact(
                PreparedExactZoh::new(system, dt, LinearSplitBudget {
                    max_modes: system.n, max_omega_dt: 0.04,
                }).map_err(|error| MultirateError::Subsystem { side, error })?),
            IntegrationMethod::DiscreteGradient { maximum_rate_rad_s, maximum_rate_dt } => {
                if !nonnegative(maximum_rate_rad_s) || !positive(maximum_rate_dt)
                    || maximum_rate_dt >= core::f64::consts::PI
                    || maximum_rate_rad_s*dt > maximum_rate_dt {
                    return Err(MultirateError::Invalid("unresolved declared local pole-rate band"));
                }
                Evolution::Gonzalez(StepWorkspace::new(system)
                    .map_err(|error| MultirateError::Subsystem { side, error })?)
            }
        };
        Ok(Self { system, clock, dt, evolution, side })
    }

    fn evolve(&mut self, initial: &[f64], input: &[f64],
        cancelled: &mut impl FnMut() -> bool) -> Result<SubsystemStep, MultirateError> {
        check(cancelled)?;
        let sys = self.system;
        if initial.len()!=sys.n || input.len()!=sys.m || !finite(initial) || !finite(input) {
            return Err(MultirateError::Invalid("complete finite subsystem state and input required"));
        }
        let mut state = initial.to_vec();
        let mut output = vec![Sum::default(); sys.m];
        let (mut change, mut heat, mut work) = (Sum::default(),Sum::default(),Sum::default());
        for _ in 0..self.clock.substeps {
            check(cancelled)?;
            let (next, y, dh, loss, supplied) = match &mut self.evolution {
                Evolution::Exact(flow) => {
                    let r = flow.step(&state, input)
                        .map_err(|error| MultirateError::Subsystem { side: self.side, error })?;
                    (r.x,r.y,r.delta_h,r.dissipated,r.supplied)
                }
                Evolution::Gonzalez(workspace) => {
                    let (mut next, mut y) = (vec![0.0;sys.n],vec![0.0;sys.m]);
                    let r = workspace.step_into_controlled(sys,&state,input,self.dt,
                        &mut next,&mut y,&mut *cancelled).map_err(|e| native(self.side,e))?;
                    (next,y,r.delta_h,r.dissipated,r.supplied)
                }
            };
            for (total,value) in output.iter_mut().zip(y) { total.add(value/self.clock.substeps as f64); }
            change.add(dh); heat.add(loss); work.add(supplied); state = next;
        }
        let result = SubsystemStep { stored_before_j: sys.hamiltonian(initial),
            stored_after_j: sys.hamiltonian(&state), state,
            average_output: output.into_iter().map(Sum::value).collect(),
            stored_change_j: change.value(), dissipated_j: heat.value(), supplied_j: work.value(),
            substeps: self.clock.substeps };
        if !finite(&result.state) || !finite(&result.average_output) || !finite(&[
            result.stored_before_j,result.stored_after_j,result.stored_change_j,
            result.dissipated_j,result.supplied_j]) {
            return Err(MultirateError::Invalid("nonfinite original subsystem trajectory"));
        }
        check(cancelled)?;
        Ok(result)
    }
}

struct Trial { left: SubsystemStep, right: SubsystemStep, residual: Vec<f64>, norm: f64 }

/// Prepared two-clock coupling. It borrows immutable physical owners and retains
/// only numerical scratch/last accepted input guess. Caller state is published
/// only in a returned successful result. Quadratic pairs cache one affine LU;
/// nonlinear pairs reuse the original implicit solver with bounded interface Newton.
pub struct MultiratePair<'a> {
    left: LocalClock<'a>,
    right: LocalClock<'a>,
    pairs: Vec<PortPair>,
    window: f64,
    solve: SolveConfig,
    affine: Option<Lu>,
    guess: Vec<f64>,
}
impl<'a> MultiratePair<'a> {
    #[allow(clippy::too_many_arguments)]
    pub fn new(left: &'a PortHamiltonian, right: &'a PortHamiltonian,
        pairs: &[PortPair], duration_s: f64, left_clock: Clock, right_clock: Clock,
        solve: SolveConfig) -> Result<Self, MultirateError> {
        if !positive(duration_s) || pairs.is_empty() || pairs.len()>MAX_INTERFACE_PORTS
            || solve.max_iterations>64 || !positive(solve.relative_tolerance)
            || !nonnegative(solve.absolute_energy_tolerance_j)
            || !nonnegative(solve.relative_energy_tolerance) {
            return Err(MultirateError::Invalid("invalid window, interface count or solve budget"));
        }
        for (i,pair) in pairs.iter().enumerate() {
            if pair.left>=left.m || pair.right>=right.m || !positive(pair.input_scale)
                || pairs[..i].iter().any(|p| p.left==pair.left || p.right==pair.right) {
                return Err(MultirateError::Invalid("invalid, duplicate or unscaled interface port"));
            }
        }
        let mut result = Self { left:LocalClock::new(left,left_clock,duration_s,"left")?,
            right:LocalClock::new(right,right_clock,duration_s,"right")?, pairs:pairs.to_vec(),
            window:duration_s, solve, affine:None, guess:vec![0.0;pairs.len()] };
        if left_clock.method==IntegrationMethod::ExactQuadratic
            && right_clock.method==IntegrationMethod::ExactQuadratic {
            let k = pairs.len();
            let (xl,xr,ul,ur) = (vec![0.0;left.n],vec![0.0;right.n],vec![0.0;left.m],vec![0.0;right.m]);
            let mut jac = vec![0.0;k*k];
            for col in 0..k {
                let mut z = vec![0.0;k]; z[col]=1.0;
                let trial = result.evaluate(&xl,&xr,&ul,&ur,&z,&mut || false)?;
                for row in 0..k { jac[row*k+col]=trial.residual[row]; }
            }
            if !finite(&jac) { return Err(MultirateError::Invalid("nonfinite affine interface response")); }
            result.affine=Some(lu(&jac,k).map_err(|_| MultirateError::Invalid("singular affine interface"))?);
        }
        Ok(result)
    }

    fn evaluate(&mut self, xl:&[f64], xr:&[f64], external_l:&[f64], external_r:&[f64],
        z:&[f64], cancelled:&mut impl FnMut()->bool) -> Result<Trial,MultirateError> {
        check(cancelled)?;
        let mut ul=external_l.to_vec();
        for (pair,&value) in self.pairs.iter().zip(z) { ul[pair.left]+=pair.input_scale*value; }
        let left=self.left.evolve(xl,&ul,cancelled)?;
        let mut ur=external_r.to_vec();
        for pair in &self.pairs { ur[pair.right]+=left.average_output[pair.left]; }
        let right=self.right.evolve(xr,&ur,cancelled)?;
        let residual:Vec<f64>=self.pairs.iter().zip(z).map(|(pair,&value)|
            value+right.average_output[pair.right]/pair.input_scale).collect();
        let scale=self.pairs.iter().zip(z).fold(1.0_f64,|s,(pair,&value)|
            s.max(value.abs()).max((right.average_output[pair.right]/pair.input_scale).abs()));
        let norm=residual.iter().map(|v|v.abs()).fold(0.0_f64,f64::max)/scale;
        if !finite(&residual) || !norm.is_finite() || !scale.is_finite() {
            return Err(MultirateError::Invalid("nonfinite interface equation"));
        }
        Ok(Trial {left,right,residual,norm})
    }

    pub fn step(&mut self, left:&[f64], right:&[f64], left_external:&[f64],
        right_external:&[f64]) -> Result<CoupledStep,MultirateError> {
        self.step_controlled(left,right,left_external,right_external,|| false)
    }

    /// Cooperative cancellation is checked between all native steps and Newton
    /// probes and before publication. Individual dense solves/callbacks are
    /// bounded by their owners but not preemptible. No hard latency is promised.
    pub fn step_controlled(&mut self, left:&[f64], right:&[f64], left_external:&[f64],
        right_external:&[f64], mut cancelled:impl FnMut()->bool) -> Result<CoupledStep,MultirateError> {
        check(&mut cancelled)?;
        if left.len()!=self.left.system.n || right.len()!=self.right.system.n
            || left_external.len()!=self.left.system.m || right_external.len()!=self.right.system.m
            || !finite(left) || !finite(right) || !finite(left_external) || !finite(right_external) {
            return Err(MultirateError::Invalid("complete finite paired state and external input required"));
        }
        let mut z=self.guess.clone();
        let k=self.pairs.len();
        for iteration in 0..=self.solve.max_iterations {
            let trial=self.evaluate(left,right,left_external,right_external,&z,&mut cancelled)?;
            if trial.norm<=self.solve.relative_tolerance {
                let result=self.finish(left,right,left_external,right_external,&z,trial,iteration)?;
                check(&mut cancelled)?;
                self.guess=z;
                return Ok(result);
            }
            if iteration==self.solve.max_iterations {
                return Err(MultirateError::NotConverged {iterations:iteration,residual:trial.norm});
            }
            let mut correction=trial.residual;
            if let Some(affine)=&self.affine { affine.solve(&mut correction); }
            else {
                let mut jac=vec![0.0;k*k];
                for col in 0..k {
                    check(&mut cancelled)?;
                    let h=1e-6*(1.0+z[col].abs());
                    let (mut plus,mut minus)=(z.clone(),z.clone());
                    plus[col]+=h; minus[col]-=h;
                    if !h.is_finite() || plus[col]==z[col] || minus[col]==z[col] {
                        return Err(MultirateError::Invalid("unrepresentable interface Jacobian probe"));
                    }
                    let a=self.evaluate(left,right,left_external,right_external,&plus,&mut cancelled)?;
                    let b=self.evaluate(left,right,left_external,right_external,&minus,&mut cancelled)?;
                    for row in 0..k { jac[row*k+col]=(a.residual[row]-b.residual[row])/(2.0*h); }
                }
                if !finite(&jac) { return Err(MultirateError::Invalid("nonfinite interface Jacobian")); }
                lu(&jac,k).map_err(|_| MultirateError::Invalid("singular nonlinear interface"))?.solve(&mut correction);
            }
            for (value,delta) in z.iter_mut().zip(correction) { *value-=delta; }
            if !finite(&z) { return Err(MultirateError::Invalid("nonfinite interface Newton update")); }
        }
        unreachable!("bounded interface loop always returns")
    }

    fn gate(&self, side:&'static str, residual:f64, scale:f64) -> Result<(),MultirateError> {
        let tolerance=self.solve.absolute_energy_tolerance_j+self.solve.relative_energy_tolerance*scale;
        if !finite(&[residual,scale,tolerance]) || residual.abs()>tolerance {
            Err(MultirateError::Energy {side,residual_j:residual,tolerance_j:tolerance})
        } else { Ok(()) }
    }

    #[allow(clippy::too_many_arguments)]
    fn finish(&self, left_initial:&[f64], right_initial:&[f64], ul:&[f64], ur:&[f64],
        z:&[f64], trial:Trial, iterations:usize) -> Result<CoupledStep,MultirateError> {
        let mut scales=[0.0;2];
        for (index,(owner,initial,result)) in [(&self.left,left_initial,&trial.left),
            (&self.right,right_initial,&trial.right)].into_iter().enumerate() {
            let sys=owner.system;
            let change=storage_energy_change(sys.storage.as_ref(),initial,&result.state,
                &mut vec![0.0;sys.n],&mut vec![0.0;sys.n]);
            let scale=energetic_scale(sys,initial).max(energetic_scale(sys,&result.state))
                .max(result.stored_change_j.abs()).max(result.dissipated_j.abs()).max(result.supplied_j.abs());
            self.gate(owner.side,change-result.stored_change_j,scale)?;
            self.gate(owner.side,sum([result.stored_change_j,result.dissipated_j,-result.supplied_j]),scale)?;
            if result.dissipated_j<0.0 { self.gate(owner.side,result.dissipated_j,scale)?; }
            scales[index]=scale;
        }
        let external=[self.window*sum(ul.iter().zip(&trial.left.average_output).map(|(u,y)|u*y)),
            self.window*sum(ur.iter().zip(&trial.right.average_output).map(|(u,y)|u*y))];
        let exchange=[self.window*sum(self.pairs.iter().zip(z).map(|(pair,value)|
            value*pair.input_scale*trial.left.average_output[pair.left])),
            self.window*sum(self.pairs.iter().map(|pair|
                trial.left.average_output[pair.left]*trial.right.average_output[pair.right]))];
        for i in 0..2 { scales[i]=scales[i].max(external[i].abs()).max(exchange[i].abs()); }
        self.gate("left work split",sum([external[0],exchange[0],-trial.left.supplied_j]),scales[0])?;
        self.gate("right work split",sum([external[1],exchange[1],-trial.right.supplied_j]),scales[1])?;
        let scale=sum(scales);
        self.gate("interface exchange",sum(exchange),scale)?;
        let residual=sum([trial.left.stored_change_j,trial.right.stored_change_j,
            trial.left.dissipated_j,trial.right.dissipated_j,-external[0],-external[1]]);
        self.gate("combined energy",residual,scale)?;
        Ok(CoupledStep {left:trial.left,right:trial.right,external_work_j:external,
            exchange_work_j:exchange,balance_residual_j:residual,interface_residual:trial.norm,iterations})
    }
}
