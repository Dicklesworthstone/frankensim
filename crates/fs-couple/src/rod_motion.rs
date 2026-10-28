//! Native large-rotation rod dynamics by symmetric Hamiltonian splitting.
//!
//! fs-solid owns the source-bound strain energy and analytical forces; fs-time
//! owns free rigid rotation. This module composes half momentum kicks, exact
//! translational drift and native gyroscopic rotation, then the second kick.
//! On smooth conservative fixed-step paths this is a second-order symplectic
//! method, NOT an exactly energy-conserving discrete-gradient method. Explicit
//! elastic-wave stability and temporal refinement remain caller obligations.
//! The energy gate is a signed splitting/source-quadrature defect, never heat.
//!
//! Sources, section dimensions, reference inertia and supports stay fixed.
//! Centerline tensile links use declared endpoint mass lumping. Eccentric links,
//! moving supports, contact, material damping, damage, thermal evolution and
//! time-dependent source replacement require other mechanical models.

use std::sync::Arc;
use fs_exec::CancelGate;
use fs_solid::rod::material_field::SpatialRodError;
use fs_solid::rod::material_field::motion::{PreparedRodMotion, RodMotionPose, RodMotionResponse};
use fs_time::{LieStepError, So3, Vec3, rigid_body_momentum_step};

/// A nodal force and couple held in the WORLD frame for this interval.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct RodNodalWrench {
    pub force_n: [f64; 3],
    pub torque_nm: [f64; 3],
}

#[derive(Clone, Copy, Debug)]
pub struct RodMotionControls {
    /// Explicit upper step bound, not an inferred stability limit.
    pub maximum_step_s: f64,
    pub maximum_steps: usize,
    /// Bound on |delta kinetic + delta potential - source quadrature| [J].
    pub energy_defect_limit_j: f64,
    /// Upper bound on each free rotor's swept angle, in (0, pi/2].
    pub maximum_rotation_rad: f64,
    /// Fixed origin for orbital angular momentum and applied-moment accounting.
    pub balance_origin_m: [f64; 3],
}

#[derive(Clone, Debug)]
pub enum RodMotionError {
    Solid(SpatialRodError),
    Rotation(LieStepError),
    Invalid(&'static str),
    Energy { defect_j: f64, limit_j: f64 },
}
impl core::fmt::Display for RodMotionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Solid(e) => write!(f,"rod motion: {e}"),
            Self::Rotation(e) => write!(f,"rod motion: {e}"),
            Self::Invalid(s) => write!(f,"rod motion: {s}"),
            Self::Energy { defect_j,limit_j } => write!(f,
                "rod motion energy defect {defect_j:e} J exceeds {limit_j:e} J"),
        }
    }
}
impl std::error::Error for RodMotionError {}
impl From<SpatialRodError> for RodMotionError {
    fn from(e: SpatialRodError) -> Self { Self::Solid(e) }
}
impl From<LieStepError> for RodMotionError {
    fn from(e: LieStepError) -> Self { Self::Rotation(e) }
}

#[derive(Clone, Debug, PartialEq)]
pub struct RodMotionObservables {
    pub kinetic_energy_j: f64,
    /// Native elastic/link storage plus gravitational potential.
    pub potential_energy_j: f64,
    pub linear_momentum_kg_m_s: [f64; 3],
    /// Orbital momentum about the requested origin PLUS all section spins.
    pub angular_momentum_kg_m2_s: [f64; 3],
}

#[derive(Clone, Debug, PartialEq)]
pub struct RodSupportImpulse {
    /// Impulse ON the rod, world coordinates [N s].
    pub linear_n_s: [f64; 3],
    /// Pure support-couple impulse ON the rod, world coordinates [N m s].
    pub angular_nm_s: [f64; 3],
}

#[derive(Clone, Debug)]
pub struct RodMotionStep {
    pub before: RodMotionObservables,
    pub after: RodMotionObservables,
    /// Symmetric kick quadrature of APPLIED wrench power. It is exact for each
    /// quadratic-in-momentum kick, but not exact continuous-trajectory work.
    /// Gravity is in potential energy, never counted again in this subaccount.
    pub source_work_quadrature_j: f64,
    pub integration_defect_j: f64,
    /// Includes applied wrenches, gravity and stationary-support impulses.
    pub external_linear_impulse_n_s: [f64; 3],
    pub external_angular_impulse_nm_s: [f64; 3],
    pub linear_balance_residual_n_s: [f64; 3],
    pub angular_balance_residual_nm_s: [f64; 3],
    /// Node order; zero entries denote no support impulse, not missing nodes.
    pub support_impulses: Vec<RodSupportImpulse>,
}

/// Accepted trajectory; cloning shares sources and copies only motion/state.
#[derive(Clone, Debug)]
pub struct RodMotion {
    model: Arc<PreparedRodMotion>,
    pose: RodMotionPose,
    linear: Vec<[f64; 3]>,
    spin: Vec<[f64; 3]>,
    response: RodMotionResponse,
    time_s: f64,
    steps: usize,
}

impl RodMotion {
    /// Initialize at an explicitly supplied admitted pose and velocity field.
    /// Stationary supports require exactly zero velocity in constrained DOFs.
    pub fn new(model: Arc<PreparedRodMotion>, pose: RodMotionPose,
        velocity_m_s: &[[f64;3]], angular_velocity_body_rad_s: &[[f64;3]],
        time_s: f64, gate: &CancelGate) -> Result<Self,RodMotionError> {
        checkpoint(gate)?;
        let n=model.inertia().len();
        if !time_s.is_finite() || velocity_m_s.len()!=n || angular_velocity_body_rad_s.len()!=n
            || velocity_m_s.iter().flatten().chain(angular_velocity_body_rad_s.iter().flatten())
                .any(|v| !v.is_finite()) {
            return Err(RodMotionError::Invalid("matching finite velocities and initial time required"));
        }
        let response=model.evaluate(&pose,gate)?;
        let mut linear=vec![[0.0;3];n];
        let mut spin=vec![[0.0;3];n];
        for node in 0..n {
            checkpoint(gate)?;
            for c in 0..3 {
                if (model.fixed_positions()[node][c] && velocity_m_s[node][c]!=0.0)
                    || (model.fixed_frames()[node] && angular_velocity_body_rad_s[node][c]!=0.0) {
                    return Err(RodMotionError::Invalid("stationary supports require zero supported velocity"));
                }
                linear[node][c]=model.inertia()[node].mass_kg()*velocity_m_s[node][c];
                spin[node][c]=model.inertia()[node].principal_kg_m2()[c]*angular_velocity_body_rad_s[node][c];
            }
        }
        let state=Self { model,pose,linear,spin,response,time_s,steps:0 };
        state.observe([0.0;3])?;
        checkpoint(gate)?;
        Ok(state)
    }
    pub fn model(&self) -> &PreparedRodMotion { &self.model }
    pub fn pose(&self) -> &RodMotionPose { &self.pose }
    pub fn linear_momenta_kg_m_s(&self) -> &[[f64;3]] { &self.linear }
    pub fn body_spins_kg_m2_s(&self) -> &[[f64;3]] { &self.spin }
    pub const fn time_s(&self) -> f64 { self.time_s }
    pub const fn accepted_steps(&self) -> usize { self.steps }

    pub fn observe(&self, origin_m:[f64;3]) -> Result<RodMotionObservables,RodMotionError> {
        finite(&origin_m)?;
        let mut kinetic=0.0;
        let mut linear=[0.0;3];
        let mut angular=[0.0;3];
        for (node,inertia) in self.model.inertia().iter().enumerate() {
            let p=self.linear[node]; let spin=self.spin[node];
            finite(&p)?; finite(&spin)?;
            for c in 0..3 {
                kinetic+=0.5*p[c]*(p[c]/inertia.mass_kg())
                    +0.5*spin[c]*(spin[c]/inertia.principal_kg_m2()[c]);
                linear[c]+=p[c];
            }
            let orbital=cross(sub(self.pose.positions_m[node],origin_m),p);
            let world=rotate(self.pose.frames[node],spin)?;
            for c in 0..3 { angular[c]+=orbital[c]+world[c]; }
        }
        finite(&[kinetic,self.response.potential_energy_j])?;
        finite(&linear)?; finite(&angular)?;
        Ok(RodMotionObservables { kinetic_energy_j:kinetic,
            potential_energy_j:self.response.potential_energy_j,
            linear_momentum_kg_m_s:linear,angular_momentum_kg_m2_s:angular })
    }

    /// One bounded candidate; failure leaves self bit-for-bit reusable.
    /// Only one fresh native force evaluation is required after initialization.
    pub fn step_candidate(&self, dt:f64, loads:&[RodNodalWrench], controls:RodMotionControls,
        gate:&CancelGate) -> Result<(Self,RodMotionStep),RodMotionError> {
        checkpoint(gate)?;
        if !dt.is_finite() || dt<=0.0 || 0.5*dt==0.0
            || !controls.maximum_step_s.is_finite() || controls.maximum_step_s<=0.0 || dt>controls.maximum_step_s
            || !controls.energy_defect_limit_j.is_finite() || controls.energy_defect_limit_j<0.0
            || !controls.maximum_rotation_rad.is_finite() || controls.maximum_rotation_rad<=0.0
            || controls.maximum_rotation_rad>core::f64::consts::FRAC_PI_2
            || self.steps>=controls.maximum_steps || loads.len()!=self.model.inertia().len() {
            return Err(RodMotionError::Invalid("invalid rod interval, work gate or step/rotation budget"));
        }
        for load in loads { finite(&load.force_n)?; finite(&load.torque_nm)?; }
        let next_time=self.time_s+dt;
        if !next_time.is_finite() || next_time<=self.time_s {
            return Err(RodMotionError::Invalid("rod clock cannot represent the requested advance"));
        }
        let before=self.observe(controls.balance_origin_m)?;
        let mut candidate=self.clone();
        let mut accounts=KickAccounts::new(loads.len());
        candidate.kick(0.5*dt,loads,controls.balance_origin_m,&mut accounts,gate)?;
        for node in 0..loads.len() {
            checkpoint(gate)?;
            let inertia=candidate.model.inertia()[node];
            for c in 0..3 {
                if !candidate.model.fixed_positions()[node][c] {
                    let shift=dt*(candidate.linear[node][c]/inertia.mass_kg());
                    let old=candidate.pose.positions_m[node][c];
                    let next=old+shift;
                    if !shift.is_finite() || !next.is_finite() || (shift!=0.0 && next==old) {
                        return Err(RodMotionError::Invalid("rod translation step is unrepresentable at this world position"));
                    }
                    candidate.pose.positions_m[node][c]=next;
                }
            }
            if !candidate.model.fixed_frames()[node] {
                let pi=candidate.spin[node];
                let moments=inertia.principal_kg_m2();
                let swept=dt*(fs_math::det::hypot(fs_math::det::hypot(pi[0],pi[1]),pi[2])
                    /moments.into_iter().fold(f64::INFINITY,f64::min));
                if !swept.is_finite() || swept>controls.maximum_rotation_rad {
                    return Err(RodMotionError::Invalid("rotor swept-angle budget exceeded; reduce the interval"));
                }
                let (frame,spin)=rigid_body_momentum_step(candidate.pose.frames[node],v(pi),v(moments),dt)?;
                candidate.pose.frames[node]=frame;
                candidate.spin[node]=a(spin);
            }
        }
        candidate.response=candidate.model.evaluate(&candidate.pose,gate)?;
        candidate.kick(0.5*dt,loads,controls.balance_origin_m,&mut accounts,gate)?;
        let after=candidate.observe(controls.balance_origin_m)?;
        let defect=(after.kinetic_energy_j-before.kinetic_energy_j)
            +(after.potential_energy_j-before.potential_energy_j)-accounts.work;
        finite(&[accounts.work,defect])?;
        if defect.abs()>controls.energy_defect_limit_j {
            return Err(RodMotionError::Energy { defect_j:defect,limit_j:controls.energy_defect_limit_j });
        }
        finite(&accounts.linear)?; finite(&accounts.angular)?;
        let linear_residual=sub(sub(after.linear_momentum_kg_m_s,before.linear_momentum_kg_m_s),accounts.linear);
        let angular_residual=sub(sub(after.angular_momentum_kg_m2_s,before.angular_momentum_kg_m2_s),accounts.angular);
        finite(&linear_residual)?; finite(&angular_residual)?;
        candidate.time_s=next_time;
        candidate.steps=self.steps.checked_add(1).ok_or(RodMotionError::Invalid("rod step count overflow"))?;
        checkpoint(gate)?;
        Ok((candidate,RodMotionStep { before,after,source_work_quadrature_j:accounts.work,
            integration_defect_j:defect,external_linear_impulse_n_s:accounts.linear,
            external_angular_impulse_nm_s:accounts.angular,linear_balance_residual_n_s:linear_residual,
            angular_balance_residual_nm_s:angular_residual,support_impulses:accounts.supports }))
    }

    fn kick(&mut self,h:f64,loads:&[RodNodalWrench],origin:[f64;3],accounts:&mut KickAccounts,
        gate:&CancelGate) -> Result<(),RodMotionError> {
        for (node,load) in loads.iter().enumerate() {
            checkpoint(gate)?;
            let frame=self.pose.frames[node];
            let mass=self.model.inertia()[node].mass_kg();
            let inertia=self.model.inertia()[node].principal_kg_m2();
            let external_body=rotate(frame.inverse(),load.torque_nm)?;
            let mut support_f=[0.0;3]; let mut support_t=[0.0;3];
            for c in 0..3 {
                let total=self.response.force_n[node][c]+load.force_n[c];
                if self.model.fixed_positions()[node][c] { support_f[c]=-h*total; }
                else {
                    let old=self.linear[node][c];
                    self.linear[node][c]=h.mul_add(total,old);
                    accounts.work+=h*load.force_n[c]*(0.5*(old/mass)+0.5*(self.linear[node][c]/mass));
                }
                let torque=self.response.torque_body_nm[node][c]+external_body[c];
                if self.model.fixed_frames()[node] { support_t[c]=-h*torque; }
                else {
                    let old=self.spin[node][c];
                    self.spin[node][c]=h.mul_add(torque,old);
                    accounts.work+=h*external_body[c]*(0.5*(old/inertia[c])+0.5*(self.spin[node][c]/inertia[c]));
                }
            }
            let support_t=rotate(frame,support_t)?;
            let gravity=self.model.source().gravity_nodal_forces_n().map_or([0.0;3],|f|f[node]);
            let impulse=std::array::from_fn(|c|h*(load.force_n[c]+gravity[c])+support_f[c]);
            let orbital=cross(sub(self.pose.positions_m[node],origin),impulse);
            for c in 0..3 {
                accounts.linear[c]+=impulse[c];
                accounts.angular[c]+=orbital[c]+h*load.torque_nm[c]+support_t[c];
                accounts.supports[node].linear_n_s[c]+=support_f[c];
                accounts.supports[node].angular_nm_s[c]+=support_t[c];
            }
            finite(&self.linear[node])?; finite(&self.spin[node])?;
            finite(&accounts.supports[node].linear_n_s)?; finite(&accounts.supports[node].angular_nm_s)?;
        }
        Ok(())
    }
}

struct KickAccounts {
    work:f64,linear:[f64;3],angular:[f64;3],supports:Vec<RodSupportImpulse>,
}
impl KickAccounts {
    fn new(n:usize)->Self { Self { work:0.0,linear:[0.0;3],angular:[0.0;3],
        supports:vec![RodSupportImpulse { linear_n_s:[0.0;3],angular_nm_s:[0.0;3] };n] } }
}
fn checkpoint(gate:&CancelGate)->Result<(),RodMotionError> {
    if gate.is_requested() { Err(SpatialRodError::Cancelled.into()) } else { Ok(()) }
}
fn finite(values:&[f64])->Result<(),RodMotionError> {
    if values.iter().any(|v| !v.is_finite()) { Err(RodMotionError::Invalid("nonfinite rod motion arithmetic")) }
    else { Ok(()) }
}
fn v(x:[f64;3])->Vec3 { Vec3::new(x[0],x[1],x[2]) }
fn a(v:Vec3)->[f64;3] { [v.x,v.y,v.z] }
fn rotate(r:So3,x:[f64;3])->Result<[f64;3],RodMotionError> {
    Ok(a(r.rotate(v(x)).map_err(LieStepError::from)?))
}
fn sub(a:[f64;3],b:[f64;3])->[f64;3] { std::array::from_fn(|i|a[i]-b[i]) }
fn cross(a:[f64;3],b:[f64;3])->[f64;3] { [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]] }
