//! Full tensor-Maxwell memory coupled to structural inertia and force ports.
//!
//! Every material point supplies a fixed 3D strain map and positive quadrature
//! volume through fs-material. Equilibrium stiffness is assembled ONCE; each
//! original factor row becomes a retained internal coordinate z with energy
//! z^2/2. With displacement u, momentum p and admitted SPD inertia M:
//! u' = M^-1 p, p' = -K_eq u - C_kv M^-1 p - sum A_j^T z_j + force,
//! z_j' = A_j M^-1 p - z_j/tau_j.
//! These are the original constitutive equations, not scalar modal damping.
//!
//! The existing fs-phs exact-ZOH propagator advances the entire coupled state
//! and independently integrates port work and dashpot loss. This module adds
//! no time integrator, mode truncation, empirical loss or radiation observer.
//! Scope: fixed kinematics, diagonal or consistent physical inertia, small strain, fixed
//! material state. Plane-stress Maxwell points have an explicit alternate
//! constructor; neither path is a nonlinear shell or a radiation model.
use fs_exec::CancelGate;
use fs_material::visco_tensor::{ProjectedTensorProny, TensorPronyError};
use fs_phs::{PhsError, PortHamiltonian, QuadraticStorage, StepRecord};
use fs_phs::linear_exact::{LinearSplit, LinearSplitBudget, split_quadratic, step_render,
    MAX_LINEAR_EXPONENTIAL_BLOCK_DIM};

/// Mesh-backed constant-strain membranes with original material memory.
pub mod membrane;
/// Original tetrahedral quadrature and consistent inertia in the same runtime.
pub mod solid;

mod inertia;
use inertia::Inertia;

mod heat;
pub use heat::MaterialHeatStep;

mod plane_stress;
pub use plane_stress::{PlaneStressPoint, PlaneStressRelaxation};

/// Explicit capacity; the unchanged exact propagator imposes its own dense cap.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TensorMaterialBudget {
    /// Maximum displacement + momentum + material-memory coordinates.
    pub max_state_coordinates: usize,
    /// Maximum quadrature contributions (each retains its original material).
    pub max_material_points: usize,
}

/// Refusal without publishing partial mechanics or constitutive memory.
#[derive(Debug)]
pub enum TensorMaterialError {
    /// Invalid shapes, inertia, clock, arithmetic or capacity.
    Invalid(&'static str),
    /// Original material or strain-map admission failed.
    Material(TensorPronyError),
    /// The original dense SPD factorization refused consistent inertia.
    Mass(fs_la::factor::FactorError),
    /// Original exact propagation/admission failed; never retry another solver.
    Propagation(PhsError),
    /// Observed before publishing a complete state.
    Cancelled,
}
impl core::fmt::Display for TensorMaterialError {
    fn fmt(&self, f:&mut core::fmt::Formatter<'_>)->core::fmt::Result {
        match self {
            Self::Invalid(what)=>write!(f,"tensor material dynamics: {what}"),
            Self::Material(error)=>write!(f,"tensor material dynamics: {error}"),
            Self::Mass(error)=>write!(f,"tensor material inertia: {error}"),
            Self::Propagation(error)=>write!(f,"tensor material dynamics: {error}"),
            Self::Cancelled=>write!(f,"tensor material dynamics cancelled"),
        }
    }
}
impl std::error::Error for TensorMaterialError {}
type Result<T> = core::result::Result<T,TensorMaterialError>;
fn bad(what:&'static str)->TensorMaterialError {TensorMaterialError::Invalid(what)}
fn poll(gate:&CancelGate)->Result<()> {
    if gate.is_requested() {Err(TensorMaterialError::Cancelled)} else {Ok(())}
}
fn finite(values:&[f64])->bool {values.iter().all(|v|v.is_finite())}

/// Explicit constitutive prehistory for a supplied initial displacement.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum InitialMemory {
    /// Held long enough for all Maxwell branches to relax at that displacement.
    Relaxed,
    /// Instantaneous jump from a zero-strain relaxed material. Full-3D points
    /// initialize spring z_j=A_j*u; plane-stress points initialize viscous a_j=0.
    Instantaneous,
}

/// Immutable full-dimensional structural/material operator with force/velocity
/// ports in the original coordinate order. Dense capacity refuses, never truncates.
pub struct TensorMaterialBody {
    system: PortHamiltonian,
    split: LinearSplit,
    inertia: Inertia,
    points: Vec<ProjectedTensorProny>,
    plane_points: Vec<PlaneStressPoint>,
    plane_offsets: Vec<usize>,
}
impl TensorMaterialBody {
    /// Assemble from physical displacement coordinates [m], diagonal masses [kg],
    /// and strain columns [1/m]. No extra elastic model may be added for the SAME
    /// material: each point already supplies its complete equilibrium stiffness.
    ///
    /// Cancellation is polled between points and before/after existing dense
    /// preparation. This does not promise an interruption inside the eigensolve.
    pub fn new(masses_kg:Vec<f64>, points:Vec<ProjectedTensorProny>,
        budget:TensorMaterialBudget, gate:&CancelGate)->Result<Self>
    {
        poll(gate)?;
        let n=masses_kg.len();
        let states=full_3d_extent(n,&points,budget)?;
        let inertia=Inertia::diagonal(masses_kg)?;
        Self::from_inertia(inertia,points,states,gate)
    }

    /// Use the COMPLETE symmetric positive-definite mass matrix [kg], row-major.
    /// Momentum is M*v and kinetic energy is p^T M^-1 p/2; force ports remain
    /// conjugate to the original coordinate velocities. No diagonal mass is
    /// substituted. Capacity is checked before dense inversion or propagation.
    pub fn from_mass_matrix(coordinates:usize,mass_kg:Vec<f64>,points:Vec<ProjectedTensorProny>,
        budget:TensorMaterialBudget,gate:&CancelGate)->Result<Self>
    {
        poll(gate)?;
        let states=full_3d_extent(coordinates,&points,budget)?;
        let inertia=Inertia::consistent(coordinates,mass_kg,gate)?;
        Self::from_inertia(inertia,points,states,gate)
    }

    fn from_inertia(inertia:Inertia,points:Vec<ProjectedTensorProny>,states:usize,
        gate:&CancelGate)->Result<Self>
    {
        let n=inertia.dimension();
        poll(gate)?;
        let mut q=vec![0.0;states*states];
        let mut j=vec![0.0;states*states];
        let mut r=vec![0.0;states*states];
        let mut g=vec![0.0;states*n];
        for i in 0..n {
            for k in 0..n {q[(n+i)*states+n+k]=inertia.inverse_entry(i,k);}
            j[i*states+n+i]=1.0; j[(n+i)*states+i]=-1.0;
            g[(n+i)*n+i]=1.0;
        }
        let mut z=2*n;
        for point in &points {
            poll(gate)?;
            for i in 0..n {for k in 0..n {
                q[i*states+k]+=point.equilibrium_stiffness()[i*n+k];
                // grad H with respect to momentum is velocity. The complete
                // projected Kelvin matrix dissipates v^T C_kv v, once only.
                r[(n+i)*states+n+k]+=point.viscosity_matrix()[i*n+k];
            }}
            for branch in point.branches() {
                let rate=branch.relaxation_time_s().recip();
                if !rate.is_finite() || rate<=0.0 {return Err(bad("material relaxation rate is unrepresentable"));}
                for row in branch.rows() {
                    q[z*states+z]=1.0; r[z*states+z]=rate;
                    for i in 0..n {
                        j[z*states+n+i]=row[i]; j[(n+i)*states+z]=-row[i];
                    }
                    z+=1;
                }
            }
        }
        if !finite(&q) || !finite(&j) || !finite(&r) {return Err(bad("assembled tensor-memory operator overflow"));}
        let storage=QuadraticStorage::new(q,states).map_err(TensorMaterialError::Propagation)?;
        let system=PortHamiltonian::new(states,n,j,r,g,Box::new(storage))
            .map_err(TensorMaterialError::Propagation)?;
        poll(gate)?;
        let split=split_quadratic(&system,LinearSplitBudget {max_modes:states,
            max_omega_dt:core::f64::consts::PI}).map_err(TensorMaterialError::Propagation)?;
        poll(gate)?;
        Ok(Self {system,split,inertia,points,plane_points:Vec::new(),plane_offsets:Vec::new()})
    }
    /// Original physical coordinate count.
    pub fn coordinates(&self)->usize {self.inertia.dimension()}
    /// Full physical inertia in the unchanged coordinate order [kg].
    pub fn mass_matrix_kg(&self)->Vec<f64> {self.inertia.matrix()}
    /// Full retained state count, including every material factor row.
    pub fn state_dimension(&self)->usize {self.system.state_dim()}
    /// Original full-3D quadrature, in assembly order; empty for plane-stress bodies.
    pub fn material_points(&self)->&[ProjectedTensorProny] {&self.points}
    /// Original force-input / velocity-output port-Hamiltonian system.
    pub fn port_system(&self)->&PortHamiltonian {&self.system}
    /// Transfer the SAME operator to existing fs-phs graph composition. The
    /// consuming graph then owns state, clocks and material-domain admission.
    pub fn into_port_system(self)->PortHamiltonian {self.system}

    fn admit_state(&self,x:&[f64])->Result<()> {
        if x.len()!=self.state_dimension() || !finite(x) {return Err(bad("nonfinite or incomplete tensor-memory state"));}
        for point in &self.points {
            point.strain(&x[..self.coordinates()]).map_err(TensorMaterialError::Material)?;
        }
        for (point,&offset) in self.plane_points.iter().zip(&self.plane_offsets) {
            point.strain(&x[..self.coordinates()],&x[offset..offset+point.memory_coordinates()])?;
        }
        if !self.system.hamiltonian(x).is_finite() {return Err(bad("material-body energy overflow"));}
        Ok(())
    }

    /// Initialize without inventing a preload or resetting an existing history.
    /// Supplied displacement/velocity and prehistory are explicit; work/heat
    /// counters start here, while the nonzero initial energy is retained.
    pub fn history(&self,displacement_m:&[f64],velocity_m_s:&[f64],
        memory:InitialMemory,gate:&CancelGate)->Result<TensorMaterialHistory<'_>>
    {
        poll(gate)?;
        let n=self.coordinates();
        if displacement_m.len()!=n || velocity_m_s.len()!=n || !finite(displacement_m) || !finite(velocity_m_s) {
            return Err(bad("initial displacement/velocity must match physical coordinates"));
        }
        let mut state=vec![0.0;self.state_dimension()];state[..n].copy_from_slice(displacement_m);
        state[n..2*n].copy_from_slice(&self.inertia.momentum(velocity_m_s));
        if memory==InitialMemory::Instantaneous {
            let mut at=2*n;
            for point in &self.points {for branch in point.branches() {
                poll(gate)?;
                for row in branch.rows() {
                    state[at]=row.iter().zip(displacement_m).fold(0.0,|s,(a,u)|s+a*u);at+=1;
                }
            }}
        }
        for (point,&offset) in self.plane_points.iter().zip(&self.plane_offsets) {
            poll(gate)?;
            // Plane-stress memory stores viscous coordinates a, not spring z.
            // Instantaneous deformation has a=0; relaxed deformation has a=B e.
            if memory==InitialMemory::Relaxed {
                for (i,row) in point.relaxed_memory.iter().enumerate() {
                    state[offset+i]=row.iter().zip(displacement_m).fold(0.0,|s,(a,u)|s+a*u);
                }
            }
        }
        self.admit_state(&state)?;poll(gate)?;
        let initial_energy=self.system.hamiltonian(&state);
        Ok(TensorMaterialHistory {body:self,state,time:0.0,steps:0,work:0.0,heat:0.0,initial_energy})
    }
}

fn full_3d_extent(n:usize,points:&[ProjectedTensorProny],budget:TensorMaterialBudget)->Result<usize> {
    if n==0 || n>128 || budget.max_state_coordinates==0
        || budget.max_state_coordinates>MAX_LINEAR_EXPONENTIAL_BLOCK_DIM
        || budget.max_material_points==0 || budget.max_material_points>4096
        || points.is_empty() || points.len()>budget.max_material_points
        || points.iter().any(|p|p.coordinates()!=n) {
        return Err(bad("bounded matching material points and finite positive inertia required"));
    }
    let branches=points.iter().try_fold(0_usize,|sum,p|sum.checked_add(p.branches().len()))
        .ok_or_else(||bad("material branch count overflow"))?;
    let states=branches.checked_mul(6).and_then(|v|v.checked_add(2*n))
        .ok_or_else(||bad("material state extent overflow"))?;
    let block=states.checked_add(n+1).and_then(|v|v.checked_mul(2))
        .ok_or_else(||bad("material exponential extent overflow"))?;
    if states>budget.max_state_coordinates || block>MAX_LINEAR_EXPONENTIAL_BLOCK_DIM {
        return Err(bad("complete tensor-memory dynamics exceed dense capacity; no coordinates removed"));
    }
    Ok(states)
}

/// Accepted state borrows its immutable body; memory cannot be rebound to a
/// different material. Cloning explicitly forks the complete accepted history.
#[derive(Clone)]
pub struct TensorMaterialHistory<'a> {
    body:&'a TensorMaterialBody,
    state:Vec<f64>,
    time:f64,
    steps:u64,
    work:f64,
    heat:f64,
    initial_energy:f64,
}
impl TensorMaterialHistory<'_> {
    /// Read-only [displacement, momentum, material coordinates]. Full-3D
    /// points use spring z; plane-stress points use viscous a coordinates.
    pub fn state(&self)->&[f64] {&self.state}
    /// Accepted physical displacements [m].
    pub fn displacement_m(&self)->&[f64] {&self.state[..self.body.coordinates()]}
    /// Accepted physical velocities [m/s].
    pub fn velocity_m_s(&self)->Vec<f64> {
        let n=self.body.coordinates();self.body.inertia.velocity(&self.state[n..2*n])
    }
    /// Accepted physical time [s].
    pub fn time_s(&self)->f64 {self.time}
    /// Number of accepted positive-duration steps.
    pub fn steps(&self)->u64 {self.steps}
    /// Current kinetic + equilibrium + recoverable branch energy [J].
    pub fn stored_energy_j(&self)->f64 {self.body.system.hamiltonian(&self.state)}
    /// Integrated externally supplied work since initialization [J].
    pub fn supplied_work_j(&self)->f64 {self.work}
    /// Independently integrated material dashpot heat since initialization [J].
    pub fn dissipated_heat_j(&self)->f64 {self.heat}
    /// Complete history balance, not a correction applied to any channel [J].
    pub fn balance_residual_j(&self)->f64 {self.stored_energy_j()-self.initial_energy+self.heat-self.work}

    /// Advance the full coupled system under a constant force during dt.
    /// The existing exact-ZOH propagator owns state, work and heat integration.
    /// Endpoint strains are checked against every original material domain;
    /// excursions between endpoints are not certified. Cancellation is observed
    /// before/after dense propagation, with no wall-clock latency claim.
    /// Any refusal preserves state, memory, clock and accumulated work/heat.
    pub fn advance(&mut self,force_n:&[f64],dt_s:f64,gate:&CancelGate)->Result<StepRecord> {
        poll(gate)?;
        if !dt_s.is_finite() || dt_s<=0.0 || force_n.len()!=self.body.coordinates() || !finite(force_n) {
            return Err(bad("positive finite clock and complete finite force vector required"));
        }
        let time=self.time+dt_s;
        let steps=self.steps.checked_add(1).ok_or_else(||bad("material-body sample counter overflow"))?;
        if !time.is_finite() || time<=self.time {return Err(bad("material-body time is unrepresentable"));}
        let record=step_render(&self.body.system,&self.body.split,&self.state,force_n,dt_s)
            .map_err(TensorMaterialError::Propagation)?;
        poll(gate)?;self.body.admit_state(&record.x)?;
        let work=self.work+record.supplied;let heat=self.heat+record.dissipated;
        if !work.is_finite() || !heat.is_finite() {return Err(bad("material-body work/heat counter overflow"));}
        poll(gate)?;
        self.state=record.x.clone();self.time=time;self.steps=steps;self.work=work;self.heat=heat;
        Ok(record)
    }
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod kelvin_tests;
