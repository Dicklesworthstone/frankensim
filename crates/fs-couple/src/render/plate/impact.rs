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
use fs_dcontact::Obstacle;
use material::refresh::LiveContactStorage as ContactStorage;
use fs_exec::CancelGate;
use fs_material::fiber::{Uniaxial,WoolFeltState};
use fs_phs::{PortHamiltonian,Storage};
use fs_plate::shell::reduction::ShellReduction;
use crate::modal_acoustic_time::ModalAcousticState;

/// Solver-derived causal listener pressure and streamed audio.
pub mod audio;
/// Reciprocal interior acoustic modes in the existing nonlinear execution.
pub mod cavity;
/// Finite-mass, compression-only pads and opposed compliant jaws.
pub mod compliant;
/// Physical, spatially localized viscous attachments.
pub mod damping;
/// Source-bound point contact in the same nonlinear potential/port equation.
pub mod material;
/// Reciprocal configuration-dependent squeeze-film pressure ports.
pub mod squeeze;
/// Compressible thin-gap gas storage in the same coupled time equation.
pub mod gas_film;
/// Hereditary material loss in the existing shared mechanical time owner.
pub mod relaxation;
/// Passive multiport acoustic reaction in the same nonlinear time equation.
pub mod radiation;
/// Prepared modal/contact image for explicitly linear bodies.
pub mod linear;
/// Statically relaxed geometric stretching of prestressed films.
pub mod membrane;
/// Force-driven uniform edge strain, coupled to the original nonlinear head.
pub mod dilation;
/// Geometric tension modulation of fixed-end filaments, using fs-nlmodal.
pub mod string;
/// Finite-mass translating supports with reciprocal string inertia and motion.
pub mod supported;
mod tangent;
mod prepared;
pub use prepared::{PreparedImpactSystem, ImpactSubstepConfig, ImpactSubstepReport, SubsteppedImpactSystem};
pub mod felt;
pub mod striker;
use felt::FeltPad;

/// Bounded scratch for joint nonlinear films and distributed wire banks.
/// This is a work ceiling, not a physical bandwidth or real-time claim.
pub const MAX_IMPACT_MODES:usize=256;
/// Mechanical potential in a caller-declared mass-normalized basis.
#[derive(Debug,Clone)]
pub enum BodyPotential {
    /// Actual nonlinear curved-shell membrane and DKT energy.
    Shell(ShellReduction),
    /// Geometry-owned planar-film stretching, not an authored pitch bend.
    Membrane(membrane::MembranePotential),
    /// Original nonlinear film with an additional finite-inertia radial edge.
    Dilation(dilation::DilationPotential),
    /// Fixed-end Kirchhoff--Carrier filament with its original EI and modal basis.
    String(string::StringPotential),
    /// Wires carried by one moving support, with exact kinetic cross terms.
    SupportedStrings(supported::SupportedStrings),
    /// Linear reduced pencil (e.g. tensioned film). Zero frequency is an
    /// explicitly free inertial coordinate, not a low-frequency oscillator.
    Linear(Vec<f64>),
}
impl BodyPotential {
    fn count(&self)->usize {match self {Self::Shell(s)=>s.mode_count(),Self::Membrane(s)=>s.reduction().mode_count(),Self::Dilation(s)=>s.reduction().mode_count(),Self::String(s)=>s.mode_count(),Self::SupportedStrings(s)=>s.mode_count(),Self::Linear(w)=>w.len()}}
    fn omegas(&self)->&[f64] {match self {Self::Shell(s)=>s.omegas(),Self::Membrane(s)=>s.reduction().omegas(),Self::Dilation(s)=>s.clock_bounds(),Self::String(s)=>s.omegas(),Self::SupportedStrings(s)=>s.clock_bounds(),Self::Linear(w)=>w}}
    fn potential(&self,q:&[f64])->f64 {match self {
        Self::Shell(s)=>s.potential(q),Self::Membrane(s)=>s.reduction().potential(q),Self::Dilation(s)=>s.reduction().potential(q),Self::String(s)=>s.potential(q),Self::SupportedStrings(s)=>s.potential(q),Self::Linear(w)=>w.iter().zip(q).map(|(w,q)|0.5*(w*q).powi(2)).sum(),
    }}
    fn gradient(&self,q:&[f64],g:&mut[f64]) {match self {
        Self::Shell(s)=>s.gradient(q,g),Self::Membrane(s)=>s.reduction().gradient(q,g),Self::Dilation(s)=>s.reduction().gradient(q,g),Self::String(s)=>s.gradient(q,g),Self::SupportedStrings(s)=>s.gradient(q,g),Self::Linear(w)=>{for ((g,w),q) in g.iter_mut().zip(w).zip(q) {*g=w*w*q;}}
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
    /// The complete retained linear block exceeds the clock's one-cent band.
    TimeResolution(fs_phs::linear_exact::ResolvedBandError),
    /// Admission or physical-limit failure.
    Invalid(&'static str),
    /// Existing owner rejected the solve or constitutive construction.
    Owner(String),
    /// Prepared numerical owner refusal, without allocating an error string.
    PreparedSolve(fs_phs::PhsError),
    /// A source-bound normal contact left its declared physical operating box.
    MaterialContact { contact: usize, source: fs_contact::normal_patch::NormalPatchError },
    /// Gate requested cancellation before publication.
    Cancelled,
    /// An output tick exhausted its explicit internal nonlinear-solve budget.
    SubstepBudget { attempted_solves: usize, accepted_substeps: usize, deepest_level: u8 },
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
    mechanical:Rc<MechanicalStorage>,contact:Rc<ContactStorage>,contact_loss:bool,
    material_contacts:Vec<material::MaterialContact>,
    relaxation:Option<relaxation::Memory>,
    radiation:Option<radiation::Memory>,
    squeeze_films:Vec<squeeze::ResistiveFilm>,
    gas_film:Option<Rc<gas_film::Memory>>,
    // Immutable shared laws plus body/start addresses; no duplicate mesh data.
    membranes:Vec<(usize,usize,membrane::MembranePotential)>,
    strings:Vec<(usize,usize,string::StringPotential)>,
    supports:Vec<(usize,usize,supported::SupportedStrings)>,
}
impl ImpactSystem {
    /// Compose real body storage, elastic contacts, felt patches and fluid volume.
    /// Initial precompression/conditioning are declared initial energy. A loaded
    /// static equilibrium is NOT silently manufactured by this constructor.
    /// Hunt-Crossley loss enters the same implicit equation through the existing
    /// contact owner's passive port, never as elastic storage or lagged forcing.
    pub fn new(bodies:Vec<ImpactBody>,contacts:Vec<Obstacle>,pads:Vec<FeltPad>,volumes:Vec<VolumeSpring>,
        config:ImpactConfig)->Result<Self,ImpactError> {
        Self::new_with_dampers(bodies,contacts,pads,volumes,Vec::new(),config)
    }
    /// Compose the same storage/contact/felt laws with localized viscous ports.
    /// Their full reciprocal resistance enters the existing implicit solve and
    /// its dissipation ledger. No modal diagonalization or output fade is used.
    pub fn new_with_dampers(bodies:Vec<ImpactBody>,contacts:Vec<Obstacle>,pads:Vec<FeltPad>,
        volumes:Vec<VolumeSpring>,dampers:Vec<damping::ViscousDamper>,
        config:ImpactConfig)->Result<Self,ImpactError> {
        let modes=bodies.iter().try_fold(0usize,|n,b|n.checked_add(b.potential.count()))
            .ok_or_else(||invalid("mode count overflow"))?;
        // A distributed cavity supplies one spring per retained pressure mode, up to its own basis limit.
        if modes==0 || modes>MAX_IMPACT_MODES || contacts.len()>32 || pads.len()>16
            || volumes.len()>cavity::MAX_CAVITY_PRESSURE_MODES
            || config.max_steps==0 || config.max_steps>(1u64<<53)
            || ![config.dt_s,config.maximum_energy_j,config.energy_absolute_tolerance_j,
                config.energy_relative_tolerance,config.maximum_generalized_force].iter().all(|v|v.is_finite() && *v>0.0)
            || config.energy_relative_tolerance>=1.0 || !(config.dt_s*config.max_steps as f64).is_finite() {
            return Err(invalid("impact needs bounded nonempty modes and finite positive time/energy/work limits"));
        }
        damping::validate(&dampers,modes)?;
        let mut x=vec![0.0;2*modes];let mut damping=Vec::with_capacity(modes);let mut potentials=Vec::new();
        let mut offset=0;let mut membranes=Vec::new();let mut strings=Vec::new();let mut supports=Vec::new();
        let mut stiffness_bound=0.0_f64;
        for (body_index,body) in bodies.into_iter().enumerate() {
            let n=body.potential.count();
            if n==0 || body.initial.len()!=n || body.damping_per_s.len()!=n
                || body.potential.omegas().iter().any(|w|!w.is_finite() || *w<0.0 || *w*config.dt_s>=0.9*core::f64::consts::PI)
                || body.damping_per_s.iter().any(|d|!d.is_finite() || *d<0.0) {
                return Err(invalid("body mode/state/damping dimensions or linear Nyquist guard failed"));
            }
            for &omega in body.potential.omegas() {stiffness_bound=stiffness_bound.max(omega*omega);}
            for (i,s) in body.initial.iter().enumerate() {x[2*(offset+i)]=s.displacement_m_sqrt_kg;x[2*(offset+i)+1]=s.velocity_m_sqrt_kg_per_s;}
            if let BodyPotential::Membrane(film)=&body.potential {
                film.observe_interleaved(&x,offset)?;
                membranes.push((body_index,offset,film.clone()));
            }
            if let BodyPotential::String(string)=&body.potential {
                string.observe_interleaved(&x,offset)?;
                strings.push((body_index,offset,string.clone()));
            }
            if let BodyPotential::SupportedStrings(support)=&body.potential {
                if body.damping_per_s.iter().any(|d|*d!=0.0) {
                    return Err(invalid("supported wire damping must use its transformed internal resistance"));
                }
                support.observe_interleaved(&x,offset)?;
                supports.push((body_index,offset,support.clone()));
            }
            damping.extend(body.damping_per_s);offset+=n;potentials.push(body.potential);
        }
        let mut admitted=Vec::with_capacity(contacts.len());let mut contact_loss=false;
        for ob in contacts {
            if ob.provenance().trim().is_empty() || ob.n_points()>4096 {
                return Err(invalid("impact contacts require explicit provenance and bounded collocation"));
            }
            let checked=Obstacle::new(ob.collocation().to_vec(),ob.n_points(),modes,ob.gaps().to_vec(),
                ob.weights().to_vec(),ob.stiffness(),ob.alpha(),ob.provenance().to_string())
                .and_then(|checked|checked.with_internal_loss(ob.internal_loss()))
                .map_err(|e|ImpactError::Owner(e.to_string()))?;
            contact_loss|=checked.internal_loss()>0.0;admitted.push(checked);
        }
        for v in &volumes {if ![v.bulk_modulus_pa,v.volume_m3].iter().all(|x|x.is_finite() && *x>0.0)
            || !(v.bulk_modulus_pa/v.volume_m3).is_finite() || v.areas.len()!=modes || v.areas.iter().any(|a|!a.is_finite()) {
            return Err(invalid("volume spring requires positive finite fluid/volume and reciprocal modal areas"));}}
        // PSD dyads from every volume spring add to the original solid pencil.
        // Their sum B^T B (rows sqrt(K/V) a) shares its nonzero spectrum with the Gram
        // matrix B B^T. Gershgorin on that Gram matrix keeps the bound near the largest
        // spring when springs act on nearly orthogonal coordinates, as a distributed
        // cavity's pressure modes do; the trace bound grows with their count. Both are
        // upper bounds, so the clock admits by the smaller.
        let mut trace_bound=stiffness_bound;
        for volume in &volumes {
            trace_bound+=(volume.bulk_modulus_pa/volume.volume_m3)
                *volume.areas.iter().map(|a|a*a).sum::<f64>();
        }
        let roots:Vec<f64>=volumes.iter().map(|v|(v.bulk_modulus_pa/v.volume_m3).sqrt()).collect();
        let gram_bound=stiffness_bound+volumes.iter().zip(&roots).map(|(u,ru)|volumes.iter().zip(&roots)
            .map(|(v,rv)|(ru*rv*u.areas.iter().zip(&v.areas).map(|(a,b)|a*b).sum::<f64>()).abs()).sum::<f64>())
            .fold(0.0_f64,f64::max);
        stiffness_bound=trace_bound.min(gram_bound);
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
        damping::add_resistance(&dampers,modes,dim,&mut r)?;
        for (_,offset,support) in &supports {support.add_resistance(*offset,dim,&mut r)?;}
        let resistance_bound=r.chunks_exact(dim).map(|row|row.iter().map(|x|x.abs()).sum::<f64>())
            .fold(0.0_f64,f64::max);
        fs_phs::linear_exact::admit_midpoint_band(stiffness_bound.sqrt()+resistance_bound,config.dt_s)
            .map_err(ImpactError::TimeResolution)?;
        let histories=Rc::new(RefCell::new(histories));
        let mechanical=Rc::new(MechanicalStorage{bodies:potentials,modes,pads:retained.clone(),volumes,histories:Rc::clone(&histories)});
        mechanical.validate_rigid_motion(&x)?;
        let contact=Rc::new(ContactStorage::new(Box::new(tangent::SharedStorage(Rc::clone(&mechanical))),modes,admitted)
            .map_err(|e|ImpactError::Owner(e.to_string()))?);
        let system=PortHamiltonian::new(dim,modes,j,r,g,Box::new(tangent::SharedStorage(Rc::clone(&contact))))
            .map_err(|e|ImpactError::Owner(e.to_string()))?;
        let energy=system.hamiltonian(&x);
        if !energy.is_finite() || energy<0.0 || energy>config.maximum_energy_j {return Err(invalid("initial impact energy exceeds admission"));}
        Ok(Self{system,x,pads:retained,histories,modes,config,sample:0,membranes,strings,supports,mechanical,contact,contact_loss,material_contacts:Vec::new(),radiation:None,relaxation:None,squeeze_films:Vec::new(),gas_film:None})
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
    /// Small-angle rigid rotation [rad] of an explicitly enabled shell body.
    /// The same physical rotation is bounded before initial or endpoint state
    /// acceptance. Ordinary shell images and non-shell bodies return None.
    pub fn shell_rotation(&self,body:usize)->Option<[f64;3]> {
        let BodyPotential::Shell(shell)=self.mechanical.bodies.get(body)? else {return None;};
        let offset:usize=self.mechanical.bodies[..body].iter().map(BodyPotential::count).sum();
        let mut q=[0.0;MAX_IMPACT_MODES];
        for (i,value) in q[..shell.mode_count()].iter_mut().enumerate() {*value=self.x[2*(offset+i)];}
        shell.rigid_rotation(&q[..shell.mode_count()]).ok().flatten()
    }
    /// Physical stretching diagnostics for a body in original construction order.
    /// Returns None for non-membrane bodies; observing consumes no mechanical time.
    pub fn membrane_observation(&self,body:usize)->Option<membrane::MembraneObservation> {
        let (_,offset,film)=self.membranes.iter().find(|(index,_,_)|*index==body)?;
        film.observe_interleaved(&self.x,*offset).ok()
    }
    /// Physical tension/stretching diagnostics in original body order.
    /// None denotes a body without the explicit nonlinear-string potential.
    pub fn string_observation(&self,body:usize)->Option<string::StringObservation> {
        let (_,offset,string)=self.strings.iter().find(|(index,_,_)|*index==body)?;
        string.observe_interleaved(&self.x,*offset).ok()
    }
    /// Global coordinate and inverse-root-mass weight for a physical force
    /// on a translating support. No direct force on its relative wire modes.
    pub fn support_force_port(&self,body:usize)->Option<(usize,f64)> {
        let (_,offset,support)=self.supports.iter().find(|(index,_,_)|*index==body)?;
        Some((offset+support.support_coordinate(),support.force_weight()))
    }
    /// Physical support and carried-wire observations by original body index.
    pub fn support_observation(&self,body:usize)->Option<supported::SupportObservation> {
        let (_,offset,support)=self.supports.iter().find(|(index,_,_)|*index==body)?;
        support.observe_interleaved(&self.x,*offset).ok()
    }
    /// Physical radial force address and canonical weight for an actuated film.
    /// Head, contact and acoustic source indices retain their original prefix.
    pub fn dilation_force_port(&self,body:usize)->Option<(usize,f64)> {
        let BodyPotential::Dilation(film)=self.mechanical.bodies.get(body)? else {return None;};
        let offset=self.mechanical.bodies[..body].iter().map(BodyPotential::count).sum::<usize>();
        Some((offset+film.reduction().rim_coordinate(),film.reduction().force_weight()))
    }
    /// Uniform edge motion/tension of one original body, without advancing time.
    pub fn dilation_observation(&self,body:usize)->Option<dilation::DilationObservation> {
        let BodyPotential::Dilation(film)=self.mechanical.bodies.get(body)? else {return None;};
        let offset=self.mechanical.bodies[..body].iter().map(BodyPotential::count).sum::<usize>();
        let n=film.reduction().mode_count();let mut q=[0.0;MAX_IMPACT_MODES];
        for i in 0..n {q[i]=self.x[2*(offset+i)];}film.observe(&q[..n]).ok()
    }
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
    /// solve. Loss-free reference execution checks cancellation before and after
    /// Newton; the implicit contact-loss path also polls inside its workspace.
    /// All physical state and felt histories remain unchanged on refusal.
    pub fn step(
        &mut self,
        external: &[f64],
        gate: &CancelGate,
    ) -> Result<ImpactFrame, ImpactError> {
        if gate.is_requested() {return Err(ImpactError::Cancelled);}
        if self.sample>=self.config.max_steps {return Err(ImpactError::Budget);}
        if external.len()!=self.modes || external.iter().any(|v|!v.is_finite() || v.abs()>self.config.maximum_generalized_force) {
            return Err(invalid("external generalized force shape or ceiling failed"));
        }
        let before=self.stored_energy_j();
        if self.has_nonlinear_dissipation() {
            // Reference execution owns temporary scratch; prepare() retains it.
            // Both use the SAME equation, port law, and physical acceptance gate.
            let mut workspace=fs_phs::StepWorkspace::new(&self.system).map_err(ImpactError::PreparedSolve)?;
            let mut candidate=vec![0.0;self.x.len()];let mut output=vec![0.0;self.modes];
            let ledger=workspace.step_into_dissipative_controlled(&self.system,&self.x,external,
                self.config.dt_s,&mut candidate,&mut output,
                &|x,e,out|self.dissipative_flow_into(x,e,out),None,None,||gate.is_requested())
                .map_err(|e|match e {fs_phs::PreparedStepError::Cancelled=>ImpactError::Cancelled,
                    fs_phs::PreparedStepError::Solver(e)=>ImpactError::PreparedSolve(e)})?;
            let mut histories=self.histories.borrow().clone();
            return self.accept_step(&mut candidate,&mut histories,before,ledger,gate);
        }
        let mut record=fs_phs::step(&self.system,&self.x,external,self.config.dt_s).map_err(|e|ImpactError::Owner(e.to_string()))?;
        if record.x.iter().chain(&record.y).any(|v|!v.is_finite()) {return Err(invalid("impact solve left finite state"));}
        let ledger=fs_phs::PreparedStepRecord {delta_h:record.delta_h,dissipated:record.dissipated,
            supplied:record.supplied,newton_iters:record.newton_iters,solver_residual:record.solver_residual,
            solver_certificate:record.solver_certificate};
        let mut histories=self.histories.borrow().clone();
        self.accept_step(&mut record.x,&mut histories,before,ledger,gate)
    }
}

impl MechanicalStorage {
    // Share one allocation-free model-domain check across initial construction,
    // reference, prepared, analytic and substepped endpoint acceptance. It runs
    // before felt/material history or the output clock can be committed.
    fn validate_rigid_motion(&self,state:&[f64])->Result<(),ImpactError> {
        if state.len()<2*self.modes {return Err(invalid("shell motion state is incomplete"));}
        let mut q=[0.0;MAX_IMPACT_MODES];let mut offset=0;
        for body in &self.bodies {
            let n=body.count();
            if let BodyPotential::Shell(shell)=body {
                for (i,value) in q[..n].iter_mut().enumerate() {*value=state[2*(offset+i)];}
                shell.rigid_rotation(&q[..n]).map_err(|e|ImpactError::Owner(e.to_string()))?;
            }
            if let BodyPotential::Dilation(film)=body {
                for (i,value) in q[..n].iter_mut().enumerate() {*value=state[2*(offset+i)];}
                film.observe(&q[..n])?;
            }
            offset+=n;
        }
        Ok(())
    }
}

#[cfg(test)]
mod contact_loss_tests;

#[cfg(test)]
mod integrator_band_tests {
    use super::*;
    #[test]
    fn full_impact_band_includes_fluid_springs_and_spatial_dampers() {
        let cap=fs_phs::linear_exact::MIDPOINT_ONE_CENT_ANGULAR_STEP;
        let make=|dt| ImpactSystem::new_with_dampers(vec![ImpactBody {
            potential:BodyPotential::Linear(vec![1.0]),
            initial:vec![ModalAcousticState::default()],damping_per_s:vec![0.0],
        }],vec![],vec![],vec![VolumeSpring {bulk_modulus_pa:3.0,volume_m3:1.0,areas:vec![1.0]}],
        vec![damping::ViscousDamper {weights:vec![1.0],damping_n_s_m:1.0}],ImpactConfig {
            dt_s:dt,max_steps:1,maximum_energy_j:1e6,energy_absolute_tolerance_j:1e-12,
            energy_relative_tolerance:1e-10,maximum_generalized_force:1e6,
        });
        assert!(make(cap/3.0).is_ok());
        assert!(matches!(make(cap/3.0+1e-12),Err(ImpactError::TimeResolution(_))));
    }

    fn config(dt_s:f64)->ImpactConfig {
        ImpactConfig {dt_s,max_steps:1,maximum_energy_j:1e6,energy_absolute_tolerance_j:1e-12,
            energy_relative_tolerance:1e-10,maximum_generalized_force:1e6}
    }

    #[test]
    fn orthogonal_cavity_springs_admit_the_clock_of_their_own_top_frequency() {
        // Twelve springs, one per coordinate, as a cavity's pressure modes: every coordinate
        // rings at sqrt(1 + 3) = 2 rad/s. The trace bound would demand sqrt(1 + 12*3).
        let cap=fs_phs::linear_exact::MIDPOINT_ONE_CENT_ANGULAR_STEP;
        let make=|dt| ImpactSystem::new(vec![ImpactBody {
            potential:BodyPotential::Linear(vec![1.0;12]),
            initial:vec![ModalAcousticState::default();12],damping_per_s:vec![0.0;12],
        }],vec![],vec![],(0..12).map(|j|VolumeSpring {bulk_modulus_pa:3.0,volume_m3:1.0,
            areas:(0..12).map(|i|if i==j {1.0} else {0.0}).collect()}).collect(),config(dt));
        assert!(make(cap/2.0).is_ok());
        assert!(matches!(make(cap/2.0+1e-12),Err(ImpactError::TimeResolution(_))));
    }

    #[test]
    fn more_than_eight_volume_springs_store_exact_energy_gradient_and_tangent() {
        use fs_phs::Storage;
        let (modes,springs)=(12,20);
        let omega:Vec<f64>=(0..modes).map(|i|1.0+0.5*i as f64).collect();
        let volumes:Vec<_>=(0..springs).map(|j|VolumeSpring {bulk_modulus_pa:2.0+j as f64,volume_m3:1.0+0.25*j as f64,
            areas:(0..modes).map(|i|((3*i+5*j)%7) as f64-3.0).collect()}).collect();
        let initial:Vec<_>=(0..modes).map(|i|ModalAcousticState {
            displacement_m_sqrt_kg:((2*i)%5) as f64*1e-3-2e-3,velocity_m_sqrt_kg_per_s:((3*i)%4) as f64*1e-2}).collect();
        let system=ImpactSystem::new(vec![ImpactBody {potential:BodyPotential::Linear(omega.clone()),
            initial:initial.clone(),damping_per_s:vec![0.0;modes]}],vec![],vec![],volumes.clone(),config(1e-4))
            .unwrap_or_else(|e|panic!("{e}"));
        let q:Vec<f64>=initial.iter().map(|s|s.displacement_m_sqrt_kg).collect();
        let stiffness:Vec<f64>=volumes.iter().map(|v|v.bulk_modulus_pa/v.volume_m3).collect();
        let strain:Vec<f64>=volumes.iter().map(|v|v.areas.iter().zip(&q).map(|(a,q)|a*q).sum()).collect();
        let energy=initial.iter().zip(&omega).map(|(s,w)|0.5*(s.velocity_m_sqrt_kg_per_s.powi(2)
            +(w*s.displacement_m_sqrt_kg).powi(2))).sum::<f64>()
            +stiffness.iter().zip(&strain).map(|(k,s)|0.5*k*s*s).sum::<f64>();
        assert!((system.stored_energy_j()-energy).abs()<=1e-12*energy);
        let mut gradient=vec![0.0;2*modes];
        system.mechanical.gradient(system.state(),&mut gradient);
        for i in 0..modes {
            let force=omega[i]*omega[i]*q[i]+volumes.iter().zip(&stiffness).zip(&strain)
                .map(|((v,k),s)|k*s*v.areas[i]).sum::<f64>();
            assert!((gradient[2*i]-force).abs()<=1e-12*force.abs().max(1e-9),"mode {i}: {} against {force}",gradient[2*i]);
            assert_eq!(gradient[2*i+1].to_bits(),initial[i].velocity_m_sqrt_kg_per_s.to_bits());
        }
        // The energy is quadratic, so the tangent's action equals the gradient's change exactly.
        let direction:Vec<f64>=(0..2*modes).map(|k|((7*k)%9) as f64*1e-3-4e-3).collect();
        let mut action=vec![0.0;2*modes];
        assert!(system.hessian_vector(system.state(),&direction,&mut action));
        let shifted:Vec<f64>=system.state().iter().zip(&direction).map(|(x,d)|x+d).collect();
        let mut moved=vec![0.0;2*modes];
        system.mechanical.gradient(&shifted,&mut moved);
        let scale=action.iter().fold(0.0_f64,|m,v|m.max(v.abs()));
        for k in 0..2*modes {
            assert!((moved[k]-gradient[k]-action[k]).abs()<=1e-10*scale,"coordinate {k}");
        }
    }
}
