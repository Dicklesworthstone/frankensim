//! Nonlinear shell/film impact and hysteretic support composition.
//!
//! fs-plate owns geometric reduction, fs-dcontact owns elastic contact,
//! fs-material owns wool-felt history, and fs-phs owns the ONLY time integrator.
//! This reference host combines those physical participants. Initial stick
//! velocity is kinetic energy, not a prescribed strike-force waveform. The
//! built-in discrete-gradient Newton currently allocates and uses a dense FD
//! Jacobian: this is NOT a hard-real-time or allocation-free audio callback.
//! Its purpose is a physics reference for a subsequently measured fast image.
//! Returned states are mechanics, NOT fabricated microphone pressure.
use std::cell::RefCell;
use std::rc::Rc;
use fs_dcontact::{ContactStorage,Obstacle};
use fs_exec::CancelGate;
use fs_material::fiber::{Uniaxial,WoolFeltState};
use fs_phs::{PortHamiltonian,Storage};
use fs_plate::shell::reduction::ShellReduction;
use crate::modal_acoustic_time::ModalAcousticState;

/// Solver-derived causal listener pressure and streamed audio.
pub mod audio;
pub mod felt;
pub mod striker;
use felt::FeltPad;

/// Bounded stack scratch for this reference host, not a cymbal adequacy claim.
pub const MAX_IMPACT_MODES:usize=64;
/// Mechanical potential in a caller-declared mass-normalized basis.
#[derive(Debug,Clone)]
pub enum BodyPotential {
    /// Actual nonlinear curved-shell membrane and DKT energy.
    Shell(ShellReduction),
    /// Linear reduced pencil (e.g. tensioned film). Zero frequency is an
    /// explicitly free inertial coordinate, not a low-frequency oscillator.
    Linear(Vec<f64>),
}
impl BodyPotential {
    fn count(&self)->usize {match self {Self::Shell(s)=>s.mode_count(),Self::Linear(w)=>w.len()}}
    fn omegas(&self)->&[f64] {match self {Self::Shell(s)=>s.omegas(),Self::Linear(w)=>w}}
    fn potential(&self,q:&[f64])->f64 {match self {
        Self::Shell(s)=>s.potential(q),Self::Linear(w)=>w.iter().zip(q).map(|(w,q)|0.5*(w*q).powi(2)).sum(),
    }}
    fn gradient(&self,q:&[f64],g:&mut[f64]) {match self {
        Self::Shell(s)=>s.gradient(q,g),Self::Linear(w)=>{for ((g,w),q) in g.iter_mut().zip(w).zip(q) {*g=w*w*q;}}
    }}
}
/// One body in the concatenated mechanical basis.
#[derive(Debug,Clone)]
pub struct ImpactBody {
    /// Potential, from an actual upstream geometric reduction where available.
    pub potential:BodyPotential,
    /// Initial mass-normalized displacement and velocity, one per coordinate.
    pub initial:Vec<ModalAcousticState>,
    /// Explicit viscous coefficients [1/s], one per coordinate (2*zeta*omega).
    pub damping_per_s:Vec<f64>,
}
impl ImpactBody {
    /// An explicit fixed-axis effective mass. Use the returned 1/sqrt(m) weight
    /// in contacts; the supplied launch velocity is not rescaled or clipped.
    pub fn free_mass(mass_kg:f64,position_m:f64,velocity_m_s:f64)->Result<(Self,f64),ImpactError> {
        if !mass_kg.is_finite() || mass_kg<=0.0 || !position_m.is_finite() || !velocity_m_s.is_finite() {
            return Err(invalid("free striker requires finite positive mass and finite physical motion"));
        }
        let root=mass_kg.sqrt();
        Ok((Self{potential:BodyPotential::Linear(vec![0.0]),initial:vec![ModalAcousticState{
            displacement_m_sqrt_kg:position_m*root,velocity_m_sqrt_kg_per_s:velocity_m_s*root}],
            damping_per_s:vec![0.0]},1.0/root))
    }
}
/// Small-signal sealed volume coupling. The signed areas integrate physical
/// boundary displacement; changing either head loads the other through this H.
/// No cavity wave modes, leakage, radiation or thermal loss are inferred.
#[derive(Debug,Clone)]
pub struct VolumeSpring {
    /// Fluid bulk modulus [Pa], e.g. rho*c^2 from the same declared gas state.
    pub bulk_modulus_pa:f64,
    /// Reference enclosed volume [m^3].
    pub volume_m3:f64,
    /// Oriented modal surface integrals [m^2/sqrt(kg)], all mechanical modes.
    pub areas:Vec<f64>,
}
/// Explicit solve/physical limits. Tolerances reject, never repair energy.
#[derive(Debug,Clone,Copy)]
pub struct ImpactConfig {
    /// Fixed mechanical step [s]; oversampling is caller controlled.
    pub dt_s:f64,
    /// Total accepted sample budget.
    pub max_steps:u64,
    /// Absolute energy ceiling [J].
    pub maximum_energy_j:f64,
    /// Energy equation tolerance [J].
    pub energy_absolute_tolerance_j:f64,
    /// Relative energy equation tolerance.
    pub energy_relative_tolerance:f64,
    /// Maximum absolute generalized external force [N/sqrt(kg)].
    pub maximum_generalized_force:f64,
}
/// One completely accepted mechanical sample.
#[derive(Debug,Clone,Copy,Default)]
pub struct ImpactFrame {
    /// One-based accepted sample number.
    pub sample:u64,
    /// End-of-step time [s].
    pub time_s:f64,
    /// Actual mechanical/contact/felt/creep/cavity energy [J].
    pub stored_energy_j:f64,
    /// Total viscous plus irreversible felt-crush energy [J].
    pub dissipated_energy_j:f64,
    /// Only the irreversible conditioning/crush part [J].
    pub felt_crush_loss_j:f64,
    /// External port work supplied during this step [J].
    pub supplied_work_j:f64,
    /// Uncorrected total balance [J].
    pub balance_residual_j:f64,
    /// Newton residual disclosed by fs-phs, not an error bound on sound.
    pub solver_residual:f64,
}
/// Typed reference-host refusal. No accepted state/history changes on error.
#[derive(Debug)]
pub enum ImpactError {
    /// Admission or physical-limit failure.
    Invalid(&'static str),
    /// Existing owner rejected the solve or constitutive construction.
    Owner(String),
    /// Gate requested cancellation before publication.
    Cancelled,
    /// Lifetime accepted-step budget is exhausted.
    Budget,
    /// Independent complete-window energy gate failed.
    Energy { residual_j:f64,tolerance_j:f64 },
}
impl core::fmt::Display for ImpactError {
    fn fmt(&self,f:&mut core::fmt::Formatter<'_>)->core::fmt::Result {write!(f,"percussion mechanics: {self:?}")}
}
impl std::error::Error for ImpactError {}
fn invalid(what:&'static str)->ImpactError {ImpactError::Invalid(what)}
#[derive(Clone)]
struct Pad {spec:FeltPad,creep_start:usize}
struct MechanicalStorage {
    bodies:Vec<BodyPotential>,modes:usize,pads:Vec<Pad>,volumes:Vec<VolumeSpring>,
    histories:Rc<RefCell<Vec<WoolFeltState>>>,
}
impl Pad {
    fn strain(&self,x:&[f64],modes:usize)->f64 {
        let mut compression=self.spec.precompression_m;
        for (i,b) in self.spec.weights.iter().enumerate() {compression+=b*x[2*i];}
        for (i,branch) in self.spec.creep.iter().enumerate() {
            compression-=x[2*modes+self.creep_start+i]/branch.stiffness_n_m.sqrt();
        }
        compression/self.spec.thickness_m
    }
}
impl Storage for MechanicalStorage {
    fn hamiltonian(&self,x:&[f64])->f64 {
        let mut q=[0.0;MAX_IMPACT_MODES];for i in 0..self.modes {q[i]=x[2*i];}
        let mut energy=(0..self.modes).map(|i|0.5*x[2*i+1]*x[2*i+1]).sum::<f64>();
        let mut offset=0;for body in &self.bodies {energy+=body.potential(&q[offset..offset+body.count()]);offset+=body.count();}
        for volume in &self.volumes {let v=volume.areas.iter().zip(&q).map(|(a,q)|a*q).sum::<f64>();
            energy+=0.5*(volume.bulk_modulus_pa/volume.volume_m3)*v*v;}
        let history=self.histories.borrow();
        for (pad,h) in self.pads.iter().zip(history.iter()) {energy+=pad.spec.path_energy(pad.strain(x,self.modes),h);}
        for &z in &x[2*self.modes..] {energy+=0.5*z*z;}
        energy
    }
    fn gradient(&self,x:&[f64],out:&mut[f64]) {
        out.fill(0.0);let mut q=[0.0;MAX_IMPACT_MODES];let mut g=[0.0;MAX_IMPACT_MODES];
        for i in 0..self.modes {q[i]=x[2*i];out[2*i+1]=x[2*i+1];}
        let mut offset=0;for body in &self.bodies {let n=body.count();body.gradient(&q[offset..offset+n],&mut g[offset..offset+n]);offset+=n;}
        for volume in &self.volumes {let v=volume.areas.iter().zip(&q).map(|(a,q)|a*q).sum::<f64>();
            for (i,a) in volume.areas.iter().enumerate() {g[i]+=(volume.bulk_modulus_pa/volume.volume_m3)*v*a;}}
        let history=self.histories.borrow();
        for (pad,h) in self.pads.iter().zip(history.iter()) {
            let force=pad.spec.force(pad.strain(x,self.modes),h);
            for (i,b) in pad.spec.weights.iter().enumerate() {g[i]+=force*b;}
            for (i,branch) in pad.spec.creep.iter().enumerate() {
                let index=2*self.modes+pad.creep_start+i;out[index]=x[index]-force/branch.stiffness_n_m.sqrt();
            }
        }
        for i in 0..self.modes {out[2*i]=g[i];}
    }
}
/// Persistent physical reference. Mode order is body order then body mode order;
/// all contact, pad and volume weights use this same concatenation.
/// Local RefCell history is immutable during each owner solve and commits only
/// after all gates. This host is single-thread-owned, not a shared audio service.
pub struct ImpactSystem {
    system:PortHamiltonian,x:Vec<f64>,pads:Vec<Pad>,histories:Rc<RefCell<Vec<WoolFeltState>>>,
    modes:usize,config:ImpactConfig,sample:u64,
}
impl ImpactSystem {
    /// Compose real body storage, elastic contacts, felt patches and fluid volume.
    /// Initial precompression/conditioning are declared initial energy. A loaded
    /// static equilibrium is NOT silently manufactured by this constructor.
    /// Contact-internal damping is explicitly refused here; Hunt-Crossley needs
    /// a simultaneous dissipative port, not insertion in the elastic potential.
    pub fn new(bodies:Vec<ImpactBody>,contacts:Vec<Obstacle>,pads:Vec<FeltPad>,volumes:Vec<VolumeSpring>,
        config:ImpactConfig)->Result<Self,ImpactError> {
        let modes=bodies.iter().try_fold(0usize,|n,b|n.checked_add(b.potential.count()))
            .ok_or_else(||invalid("mode count overflow"))?;
        if modes==0 || modes>MAX_IMPACT_MODES || contacts.len()>32 || pads.len()>16 || volumes.len()>8
            || config.max_steps==0 || config.max_steps>(1u64<<53)
            || ![config.dt_s,config.maximum_energy_j,config.energy_absolute_tolerance_j,
                config.energy_relative_tolerance,config.maximum_generalized_force].iter().all(|v|v.is_finite() && *v>0.0)
            || config.energy_relative_tolerance>=1.0 || !(config.dt_s*config.max_steps as f64).is_finite() {
            return Err(invalid("impact needs bounded nonempty modes and finite positive time/energy/work limits"));
        }
        let mut x=vec![0.0;2*modes];let mut damping=Vec::with_capacity(modes);let mut potentials=Vec::new();
        let mut offset=0;
        for body in bodies {
            let n=body.potential.count();
            if n==0 || body.initial.len()!=n || body.damping_per_s.len()!=n
                || body.potential.omegas().iter().any(|w|!w.is_finite() || *w<0.0 || *w*config.dt_s>=0.9*core::f64::consts::PI)
                || body.damping_per_s.iter().any(|d|!d.is_finite() || *d<0.0) {
                return Err(invalid("body mode/state/damping dimensions or linear Nyquist guard failed"));
            }
            for (i,s) in body.initial.iter().enumerate() {x[2*(offset+i)]=s.displacement_m_sqrt_kg;x[2*(offset+i)+1]=s.velocity_m_sqrt_kg_per_s;}
            damping.extend(body.damping_per_s);offset+=n;potentials.push(body.potential);
        }
        let mut admitted=Vec::with_capacity(contacts.len());
        for ob in contacts {
            if ob.internal_loss()!=0.0 || ob.provenance().trim().is_empty() || ob.n_points()>4096 {
                return Err(invalid("impact contacts require explicit elastic provenance; nonzero internal loss is not silently ignored"));
            }
            admitted.push(Obstacle::new(ob.collocation().to_vec(),ob.n_points(),modes,ob.gaps().to_vec(),
                ob.weights().to_vec(),ob.stiffness(),ob.alpha(),ob.provenance().to_string()).map_err(|e|ImpactError::Owner(e.to_string()))?);
        }
        for v in &volumes {if ![v.bulk_modulus_pa,v.volume_m3].iter().all(|x|x.is_finite() && *x>0.0)
            || !(v.bulk_modulus_pa/v.volume_m3).is_finite() || v.areas.len()!=modes || v.areas.iter().any(|a|!a.is_finite()) {
            return Err(invalid("volume spring requires positive finite fluid/volume and reciprocal modal areas"));}}
        let mut retained=Vec::with_capacity(pads.len());let mut histories=Vec::with_capacity(pads.len());
        for spec in pads {
            spec.validate(modes)?;let start=x.len()-2*modes;x.resize(x.len()+spec.creep.len(),0.0);
            let pad=Pad{spec,creep_start:start};let strain=pad.strain(&x,modes);
            if !strain.is_finite() || strain>pad.spec.law.eps_densify {return Err(invalid("initial felt exceeds densification validity"));}
            histories.push(pad.spec.history_at(strain.max(pad.spec.prior_maximum_strain)));retained.push(pad);
        }
        if x.iter().any(|v|!v.is_finite()) {return Err(invalid("initial impact state is nonfinite"));}
        let dim=x.len();let mut j=vec![0.0;dim*dim];let mut r=vec![0.0;dim*dim];let mut g=vec![0.0;dim*modes];
        for i in 0..modes {j[(2*i)*dim+2*i+1]=1.0;j[(2*i+1)*dim+2*i]=-1.0;
            r[(2*i+1)*dim+2*i+1]=damping[i];g[(2*i+1)*modes+i]=1.0;}
        for pad in &retained {for (i,b) in pad.spec.creep.iter().enumerate() {let index=2*modes+pad.creep_start+i;r[index*dim+index]=b.stiffness_n_m/b.viscosity_n_s_m;}}
        let histories=Rc::new(RefCell::new(histories));
        let storage=MechanicalStorage{bodies:potentials,modes,pads:retained.clone(),volumes,histories:Rc::clone(&histories)};
        let storage=ContactStorage::new(Box::new(storage),modes,admitted).map_err(|e|ImpactError::Owner(e.to_string()))?;
        let system=PortHamiltonian::new(dim,modes,j,r,g,Box::new(storage)).map_err(|e|ImpactError::Owner(e.to_string()))?;
        let energy=system.hamiltonian(&x);
        if !energy.is_finite() || energy<0.0 || energy>config.maximum_energy_j {return Err(invalid("initial impact energy exceeds admission"));}
        Ok(Self{system,x,pads:retained,histories,modes,config,sample:0})
    }
    /// Accepted mass-normalized q,p; Kelvin coordinates follow the 2*modes prefix.
    #[must_use]
    pub fn state(&self)->&[f64] {&self.x}
    /// Accepted sample count, unaffected by failed trials.
    #[must_use]
    pub const fn samples(&self)->u64 {self.sample}
    /// Actual current storage including conditioning and creep.
    #[must_use]
    pub fn stored_energy_j(&self)->f64 {self.system.hamiltonian(&self.x)}
    /// Copy one accepted material history, not a new fitted material.
    #[must_use]
    pub fn felt_history(&self,pad:usize)->Option<WoolFeltState> {self.histories.borrow().get(pad).cloned()}
    /// Current felt strain and compression force [N]. Endpoint constitutive
    /// values are not asserted to be the step's averaged contact reaction.
    #[must_use]
    pub fn felt_observation(&self,pad:usize)->Option<(f64,f64)> {
        let p=self.pads.get(pad)?;let strain=p.strain(&self.x,self.modes);
        Some((strain,p.spec.force(strain,&self.histories.borrow()[pad])))
    }
    /// Advance one held-force sample through the existing discrete-gradient
    /// solve. Cancellation is checked before and after it, not inside its Newton.
    /// All physical state and felt histories remain unchanged on refusal.
    pub fn step(&mut self,external:&[f64],gate:&CancelGate)->Result<ImpactFrame,ImpactError> {
        if gate.is_requested() {return Err(ImpactError::Cancelled);}
        if self.sample>=self.config.max_steps {return Err(ImpactError::Budget);}
        if external.len()!=self.modes || external.iter().any(|v|!v.is_finite() || v.abs()>self.config.maximum_generalized_force) {
            return Err(invalid("external generalized force shape or ceiling failed"));
        }
        let before=self.stored_energy_j();
        let record=fs_phs::step(&self.system,&self.x,external,self.config.dt_s).map_err(|e|ImpactError::Owner(e.to_string()))?;
        if record.x.iter().chain(&record.y).any(|v|!v.is_finite()) {return Err(invalid("impact solve left finite state"));}
        let frozen=self.system.hamiltonian(&record.x);let mut crush=0.0;
        let mut candidate=self.histories.borrow().clone();
        for (pad,h) in self.pads.iter().zip(&mut candidate) {
            let strain=pad.strain(&record.x,self.modes);
            if !strain.is_finite() || strain>pad.spec.law.eps_densify {return Err(invalid("felt trial exceeds densification validity"));}
            let new=pad.spec.law.update_state(strain,h);
            let loss=pad.spec.path_energy(strain,h)-pad.spec.recovered(strain,&new);
            if !loss.is_finite() || loss < -64.0*f64::EPSILON*frozen.abs() {return Err(invalid("felt history update creates energy"));}
            crush+=loss;*h=new;
        }
        let after=frozen-crush;let dissipated=record.dissipated+crush;
        let residual=after-before+dissipated-record.supplied;
        let tolerance=self.config.energy_absolute_tolerance_j+self.config.energy_relative_tolerance
            *(before.abs()+after.abs()+dissipated.abs()+record.supplied.abs());
        if ![after,dissipated,residual,record.supplied,record.solver_residual,tolerance].iter().all(|v|v.is_finite())
            || after<0.0 || after>self.config.maximum_energy_j || record.dissipated<0.0 {
            return Err(invalid("impact candidate exceeds finite energy limits"));
        }
        if residual.abs()>tolerance {return Err(ImpactError::Energy{residual_j:residual,tolerance_j:tolerance});}
        if gate.is_requested() {return Err(ImpactError::Cancelled);}
        self.x=record.x;*self.histories.borrow_mut()=candidate;self.sample+=1;
        Ok(ImpactFrame{sample:self.sample,time_s:self.sample as f64*self.config.dt_s,stored_energy_j:after,
            dissipated_energy_j:dissipated,felt_crush_loss_j:crush,supplied_work_j:record.supplied,
            balance_residual_j:residual,solver_residual:record.solver_residual})
    }
}
