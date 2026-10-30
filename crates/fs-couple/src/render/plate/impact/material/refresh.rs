//! Replace stateless material contacts only where their stored energy is zero.
//! Every storage wrapper and tangent keeps the same live contact owner. No body,
//! acoustic state, material memory, canonical coordinate or solver is rebuilt.
use std::cell::RefCell;
use std::rc::Rc;
use fs_contact::normal_patch::source::SourcePointContact;
use fs_dcontact::{ContactStorage as Owner, DContactError, Obstacle};
use fs_exec::CancelGate;
use fs_phs::Storage;
use super::super::{ImpactSystem, ImpactError, invalid, tangent};

/// One stateless constitutive replacement at the existing physical contact.
/// Geometry, contact map, gap and body inertia are deliberately unchanged.
#[derive(Clone, Debug, PartialEq)]
pub struct ContactSourceUpdate {
    pub contact: String,
    pub source: SourcePointContact,
}

// A single indirection is necessary because acoustic/relaxation storage wraps
// the contact owner. Replacing only ImpactSystem.contact would leave those
// Hamiltonians using stale stiffness while tangents/forces used the new source.
// Read borrows cover one original owner call. Mutation is restricted to cold,
// exclusive transaction boundaries; no callback runs under a mutable borrow.
pub(crate) struct LiveContactStorage(RefCell<Owner>);
impl LiveContactStorage {
    pub(crate) fn new(inner:Box<dyn Storage>,modes:usize,obstacles:Vec<Obstacle>) -> Result<Self,DContactError> {
        Owner::new(inner,modes,obstacles).map(|owner|Self(RefCell::new(owner)))
    }
    pub(crate) fn hessian_vector_with<F>(&self,x:&[f64],d:&[f64],out:&mut[f64],f:F)->bool
    where F:FnOnce(&[f64],&[f64],&mut[f64])->bool {
        self.0.borrow().hessian_vector_with(x,d,out,f)
    }
    pub(crate) fn dissipative_flow_into(&self,x:&[f64],e:&[f64],out:&mut[f64])->bool {
        self.0.borrow().dissipative_flow_into(x,e,out)
    }
    pub(crate) fn dissipative_flow_tangent_into(&self,x:&[f64],e:&[f64],dx:&[f64],de:&[f64],out:&mut[f64])->bool {
        self.0.borrow().dissipative_flow_tangent_into(x,e,dx,de,out)
    }
    pub(crate) fn dissipated_power_by_obstacle_into(&self,x:&[f64],e:&[f64],out:&mut[f64])->bool {
        self.0.borrow().dissipated_power_by_obstacle_into(x,e,out)
    }
}
impl Storage for LiveContactStorage {
    fn hamiltonian(&self,x:&[f64])->f64 {self.0.borrow().hamiltonian(x)}
    fn gradient(&self,x:&[f64],out:&mut[f64]) {self.0.borrow().gradient(x,out);}
}

pub(crate) struct ContactRefresh {
    owner:Owner,
    sources:Vec<super::MaterialContact>,
    loss:bool,
}
impl ImpactSystem {
    /// Refresh selected source laws at separated/touching, uncompressed contacts.
    ///
    /// Every changed contact must have EXACTLY zero indentation in the accepted
    /// state. No tolerance makes a compressed contact unloaded. Otherwise changing
    /// stiffness would require constitutive work/thermal accounting not supplied
    /// by this operation. Curvature, map, gap, source order, body masses and all
    /// physical state/history remain unchanged. The new model checks the current
    /// rate and timestep; its full source envelope was admitted at construction.
    ///
    /// All updates are admitted before any publication. A bad key, duplicate,
    /// loaded contact, geometry change, source failure or cancellation changes
    /// nothing. This operation advances no time and adds/removes no energy.
    /// Cold allocation is bounded by the existing 32-contact system limit.
    pub fn refresh_unloaded_contacts(&mut self,updates:Vec<ContactSourceUpdate>,gate:&CancelGate)
        ->Result<usize,ImpactError> {
        let count=updates.len();
        if let Some(next)=self.stage_contact_refresh(updates,gate)? {self.swap_contact_refresh(next);}
        Ok(count)
    }

    pub(crate) fn stage_contact_refresh(&self,updates:Vec<ContactSourceUpdate>,gate:&CancelGate)
        ->Result<Option<ContactRefresh>,ImpactError> {
        if gate.is_requested() {return Err(ImpactError::Cancelled);}
        if updates.len()>32 {return Err(invalid("contact refresh exceeds the 32-contact limit"));}
        if updates.is_empty() {return Ok(None);}
        let mut sources=self.material_contacts.clone();
        let mut raw=self.contact.0.borrow().obstacles().to_vec();
        if sources.len()!=raw.len() {return Err(invalid("refresh requires complete source-bound contact ownership"));}
        let mut visited=vec![false;sources.len()];
        for update in updates {
            if gate.is_requested() {return Err(ImpactError::Cancelled);}
            let i=sources.iter().position(|c|c.key()==update.contact)
                .ok_or_else(||invalid("unknown source contact in refresh"))?;
            if visited[i] {return Err(invalid("a contact may be refreshed only once per transaction"));}
            visited[i]=true;
            let c=&mut sources[i];
            let (d,v)=c.kinematics(&self.x,&self.x)
                .map_err(|source|ImpactError::MaterialContact {contact:i,source})?;
            if d!=0.0 {return Err(invalid("loaded contact refresh requires an explicit thermomechanical work model"));}
            if update.source.geometry()!=c.source().geometry() {
                return Err(invalid("source refresh cannot change physical contact geometry"));
            }
            update.source.validate_kinematics(d,v,self.config.dt_s)
                .map_err(|source|ImpactError::MaterialContact {contact:i,source})?;
            c.source=update.source;
            raw[i]=c.lower(self.modes)?;
        }
        let loss=raw.iter().any(|o|o.internal_loss()>0.0);
        let owner=Owner::new(Box::new(tangent::SharedStorage(Rc::clone(&self.mechanical))),self.modes,raw)
            .map_err(|e|ImpactError::Owner(e.to_string()))?;
        // Only the original mechanical prefix belongs to the contact base;
        // acoustic/relaxation/gas states are retained by their existing wrappers.
        let base=2*self.modes+self.pads.iter().map(|p|p.spec.creep.len()).sum::<usize>();
        let before=self.contact.hamiltonian(&self.x[..base]);
        let after=owner.hamiltonian(&self.x[..base]);
        if !before.is_finite() || !after.is_finite() || before.to_bits()!=after.to_bits() {
            return Err(invalid("unloaded contact replacement changed stored energy"));
        }
        if gate.is_requested() {return Err(ImpactError::Cancelled);}
        Ok(Some(ContactRefresh {owner,sources,loss}))
    }

    // Infallible swap also supplies rollback for a later coupled step refusal.
    // The complete source descriptors and the numerical owner move together.
    pub(crate) fn swap_contact_refresh(&mut self,mut next:ContactRefresh)->ContactRefresh {
        std::mem::swap(&mut *self.contact.0.borrow_mut(),&mut next.owner);
        std::mem::swap(&mut self.material_contacts,&mut next.sources);
        std::mem::swap(&mut self.contact_loss,&mut next.loss);
        next
    }
}

// Continuous local analysis is isolated from the existing refresh transaction.
// The original contact owner supplies every coefficient and physical callback.
mod resonance {
    use super::{ImpactError, ImpactSystem, LiveContactStorage, invalid};
    use fs_phs::harmonic::HarmonicError;
    use fs_phs::harmonic::equilibrium::{EquilibriumLinearization, EquilibriumPoint,
        LinearizationBudget, LinearizationError, NonlinearDissipation};

    impl LiveContactStorage {
        // Two-sided derivatives only. This checks the original contact sites,
        // not the inner body's law or a finite-amplitude admissibility interval.
        fn has_two_sided_contact_tangent(&self,x:&[f64],effort:&[f64])->bool {
            let owner=self.0.borrow();let modes=owner.n_modes();
            if modes.checked_mul(2).is_none_or(|n|x.len()<n)
                || x.len()!=effort.len() || x.iter().chain(effort).any(|v|!v.is_finite()) {return false;}
            for ob in owner.obstacles() {
                if ob.n_points().checked_mul(modes)!=Some(ob.collocation().len())
                    || ob.gaps().len()!=ob.n_points() || ob.weights().len()!=ob.n_points()
                    || ob.collocation().iter().chain(ob.gaps()).chain(ob.weights()).any(|v|!v.is_finite())
                    || !ob.stiffness().is_finite() || ob.stiffness()<0.
                    || !ob.alpha().is_finite() || ob.alpha()<1.
                    || !ob.internal_loss().is_finite() || ob.internal_loss()<0.
                    || ob.weights().iter().any(|v|*v<0.) {return false;}
                if ob.stiffness()==0. {continue;}
                for i in 0..ob.n_points() {
                    if ob.weights()[i]==0. {continue;}
                    let row=&ob.collocation()[i*modes..(i+1)*modes];
                    if row.iter().all(|b|*b==0.) {continue;}
                    // Same ordered scalar displacement sum as penetration_at,
                    // retaining the sign to distinguish separation from touching.
                    let mut displacement=0.;
                    for (k,&b) in row.iter().enumerate() {displacement+=b*x[2*k];}
                    let penetration=displacement-ob.gaps()[i];
                    if !penetration.is_finite() {return false;}
                    if penetration==0. && ob.alpha()==1. {return false;}
                    if penetration>0. && ob.internal_loss()>0. {
                        let rate=row.iter().enumerate().map(|(k,b)|b*effort[2*k+1]).sum::<f64>();
                        let factor=ob.internal_loss()*rate;
                        if !factor.is_finite() || factor== -1. {return false;}
                    }
                }
            }
            true
        }
    }

    impl ImpactSystem {
        /// Analyze a resolved stationary state of this complete physical model.
        /// Input is the original generalized force in body-mode order. Source
        /// coordinates, acoustic states, Maxwell memory, full viscous attachments,
        /// contact stiffness and nonlinear film/contact damping remain present.
        ///
        /// This read-only operation does not advance the accepted clock, commit
        /// a configuration, or refresh material/history. The result is an owned
        /// local snapshot; subsequent motion or refresh needs another preparation.
        /// The original nonlinear time equation remains authoritative for finite
        /// amplitude, contact switching and constitutive-history evolution.
        ///
        /// Irreversibly conditioned felt is refused: its loading history is not
        /// an explicit smooth state of this image. Penalty/cap corners refuse.
        /// Source validity is checked at the point, not certified over a finite
        /// oscillation envelope. Additional loss-action capacity is explicit.
        pub fn linearize_equilibrium(
            &self,point:EquilibriumPoint<'_>,budget:LinearizationBudget,
            max_loss_tangent_actions:usize,mut cancelled:impl FnMut()->bool,
        )->Result<EquilibriumLinearization,ImpactError> {
            if cancelled() {return Err(ImpactError::Cancelled);}
            if point.state.len()!=self.x.len() || point.input.len()!=self.modes
                || point.state.iter().chain(point.input).any(|x|!x.is_finite())
                || point.input.iter().any(|f|f.abs()>self.config.maximum_generalized_force) {
                return Err(invalid("resonance point requires complete finite coordinates and admitted force"));
            }
            if !self.pads.is_empty() {
                return Err(invalid("irreversible felt histories need a smooth explicit-state linearization"));
            }
            if point.state.len()>budget.harmonic.max_dimension
                || point.state.len()>budget.max_hessian_actions {
                return Err(invalid("complete mechanical resonance state exceeds its declared capacity"));
            }
            let has_loss=self.has_nonlinear_dissipation();
            if has_loss && point.state.len()>max_loss_tangent_actions {
                return Err(invalid("complete nonlinear-loss tangent exceeds its action capacity"));
            }
            self.mechanical.validate_rigid_motion(point.state)?;
            for (_,offset,film) in &self.membranes {film.observe_interleaved(point.state,*offset)?;}
            for (_,offset,string) in &self.strings {string.observe_interleaved(point.state,*offset)?;}
            for (_,offset,support) in &self.supports {support.observe_interleaved(point.state,*offset)?;}
            self.validate_squeeze_configuration(point.state)?;
            let energy=self.system.hamiltonian(point.state);
            if !energy.is_finite() || energy<0. || energy>self.config.maximum_energy_j {
                return Err(invalid("resonance point exceeds original physical energy admission"));
            }
            let mut effort=vec![f64::NAN;point.state.len()];
            self.system.storage().gradient(point.state,&mut effort);
            if cancelled() {return Err(ImpactError::Cancelled);}
            if effort.iter().any(|x|!x.is_finite()) {
                return Err(invalid("original resonance-point effort is not representable"));
            }
            self.validate_material_contacts(point.state,&effort)?;
            if !self.contact.has_two_sided_contact_tangent(point.state,&effort) {
                return Err(invalid("contact point lacks a two-sided local storage/loss derivative"));
            }
            let hessian=|x:&[f64],d:&[f64],out:&mut[f64]|self.hessian_vector(x,d,out);
            let result=if has_loss {
                let flow=|x:&[f64],e:&[f64],out:&mut[f64]|self.dissipative_flow_into(x,e,out);
                let tangent=|x:&[f64],e:&[f64],dx:&[f64],de:&[f64],out:&mut[f64]|
                    self.dissipative_flow_tangent_into(x,e,dx,de,out);
                EquilibriumLinearization::from_dissipative_ode(&self.system,point,budget,&hessian,
                    NonlinearDissipation {flow:&flow,tangent:&tangent,max_tangent_actions:max_loss_tangent_actions},
                    &mut cancelled)
            } else {
                EquilibriumLinearization::from_ode(&self.system,point,budget,&hessian,&mut cancelled)
            };
            result.map_err(|error|match error {
                LinearizationError::Cancelled|LinearizationError::Harmonic(HarmonicError::Cancelled)=>ImpactError::Cancelled,
                other=>ImpactError::Owner(other.to_string()),
            })
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use crate::render::plate::impact::{BodyPotential,ImpactBody,ImpactConfig};
        use crate::render::plate::impact::relaxation::InitialMemory;
        use crate::render::plate::impact::squeeze::{FilmCell,FilmChannel,FilmLimits,GapPort,ResistiveFilm};
        use crate::modal_acoustic_time::ModalAcousticState;
        use super::super::LiveContactStorage as ContactStorage;
        use fs_dcontact::Obstacle;
        use fs_exec::CancelGate;
        use fs_math::c64::C64;
        use fs_phs::{QuadraticStorage,RelaxationBranch};
        use fs_phs::harmonic::{HarmonicBudget,PoleOptions};
        use std::cell::Cell;

        fn config()->ImpactConfig {ImpactConfig {dt_s:1e-5,max_steps:1000,maximum_energy_j:100.,
            energy_absolute_tolerance_j:1e-12,energy_relative_tolerance:1e-7,maximum_generalized_force:1000.}}
        fn budget()->LinearizationBudget {LinearizationBudget {harmonic:HarmonicBudget {
            max_dimension:32,max_dense_bytes:2_000_000,max_dense_work:10_000_000,
            relative_residual_tolerance:1e-11},max_hessian_actions:32}}
        fn body(w:Vec<f64>,q:&[f64],d:Vec<f64>)->ImpactBody {
            ImpactBody {potential:BodyPotential::Linear(w),damping_per_s:d,
                initial:q.iter().map(|&q|ModalAcousticState {
                    displacement_m_sqrt_kg:q,velocity_m_sqrt_kg_per_s:0.,
                }).collect()}
        }
        fn obstacle(row:Vec<f64>,alpha:f64,chi:f64)->Obstacle {
            let n=row.len();
            Obstacle::new(row,1,n,vec![0.],vec![1.],1000.,alpha,
                "manufactured resonance fixture, not a measured material".into()).unwrap()
                .with_internal_loss(chi).unwrap()
        }
        fn one(preload:f64)->ImpactSystem {
            ImpactSystem::new(vec![body(vec![10.],&[preload],vec![0.2])],
                vec![obstacle(vec![1.],1.5,0.3)],vec![],vec![],config()).unwrap()
        }
        fn image(system:&ImpactSystem,force:&[f64])->EquilibriumLinearization {
            let scales=vec![10.;system.state().len()];
            system.linearize_equilibrium(EquilibriumPoint {state:system.state(),input:force,
                balance_scales:&scales,relative_tolerance:1e-12},budget(),32,||false).unwrap()
        }
        fn near(a:C64,b:C64) {assert!((a-b).abs()<2e-10*(1.+b.abs()),"{a:?} != {b:?}");}

        #[test]
        fn g1_preload_changes_contact_frequency_and_decay_from_the_original_law() {
            let mut previous=None;
            for preload in [0.01_f64,0.02,0.04] {
                let system=one(preload);
                let force=[100.*preload+1000.*preload.powf(1.5)];
                let local=image(&system,&force);
                let k=100.+1500.*preload.sqrt();let d=0.2+300.*preload.powf(1.5);
                let exact=C64::new(-d/2.,-(k-d*d/4.).sqrt());
                let pole=local.model().resonance_pole(exact,&[C64::ONE,exact],&PoleOptions {
                    spectral_scale:30.,step_tolerance:1e-13,residual_tolerance:1e-13,
                    max_iterations:24,max_dimension:32,max_dense_bytes:2_000_000,max_work:10_000_000,
                },||false).unwrap();
                near(pole.laplace_pole_per_s,exact);
                assert!((local.tangent_storage_matrix()[0]-k).abs()<1e-10);
                let t=local.model().nonlinear_resisting_tangent().unwrap();
                assert!((t[3]-(d-0.2)).abs()<1e-12);
                assert_eq!(local.equilibrium_nonlinear_dissipated_power_w(),0.);
                assert!(pole.dissipated_power_w>0. && pole.decay_power_residual_w.abs()<1e-10);
                if let Some((frequency,decay))=previous {assert!(exact.im.abs()>frequency && -exact.re>decay);}
                previous=Some((exact.im.abs(),-exact.re));assert_eq!(system.samples(),0);
            }
        }

        #[test]
        fn g1_loaded_contact_and_maxwell_memory_match_independent_dynamic_stiffness() {
            let preload=0.02_f64;
            let system=one(preload).with_relaxation_branches(vec![RelaxationBranch {
                projection:vec![1.,0.],stiffness:30.,relaxation_time_s:0.2,
            }],InitialMemory::Relaxed,1).unwrap();
            let local=image(&system,&[100.*preload+1000.*preload.powf(1.5)]);
            assert_eq!(local.equilibrium_state().len(),3);
            let k=100.+1500.*preload.sqrt();let d=0.2+300.*preload.powf(1.5);
            let force=C64::new(0.7,-0.2);
            for w in [5.,10.,20.,80.] {
                let s=C64::new(0.,-w);let memory=C64::ONE+s.scale(0.2);
                let q=force/(s*s+s.scale(d)+C64::from_re(k)+s.scale(6.)/memory);
                let out=local.model().factor(w,||false).unwrap().response(&[force],||false).unwrap();
                near(out.state[0],q);near(out.state[1],s*q);
                near(out.state[2],q.scale(30.0_f64.sqrt())/memory);
                assert!(out.power_residual_w.abs()<1e-11);
            }
            assert_eq!(system.relaxation_branch_count(),1);
            assert_eq!(system.relaxation_observation().unwrap().stored_energy_j,0.);
        }

        #[test]
        fn g3_signed_two_body_contact_and_reynolds_film_share_complete_resonance_equations() {
            let q=[0.0001,-0.0001];let penetration=0.0002_f64;
            let mut system=ImpactSystem::new(vec![body(vec![10.,15.],&q,vec![0.2,0.4])],
                vec![obstacle(vec![1.,-1.],1.5,100.)],vec![],vec![],config()).unwrap();
            let gap=GapPort {reference_m:0.001,closure:vec![1.,-1.]};
            let film=ResistiveFilm::new(vec![FilmCell {area_m2:0.002,gap:gap.clone()}],
                vec![FilmChannel {from:0,to:None,width_m:0.05,length_m:0.01,gap}],2,1.8e-5,
                FilmLimits {minimum_cell_gap_m:1e-6,maximum_gap_m:0.01,maximum_pressure_pa:1000.}).unwrap();
            system=system.with_squeeze_film(film).unwrap();
            let contact=1000.*penetration.powf(1.5);
            let local=image(&system,&[100.*q[0]+contact,225.*q[1]-contact]);
            let kc=1500.*penetration.sqrt();
            let film_damping=0.002_f64.powi(2)*12.*1.8e-5*0.01/(0.05*0.0008_f64.powi(3));
            let damping=100.*contact+film_damping;
            let t=local.model().nonlinear_resisting_tangent().unwrap();
            for (i,j,sign) in [(1,1,1.),(1,3,-1.),(3,1,-1.),(3,3,1.)] {
                assert!((t[i*4+j]-sign*damping).abs()<1e-12);
            }
            let f=[C64::new(1e-4,2e-5),C64::new(-4e-5,7e-5)];
            for w in [5.,10.,15.,30.] {
                let s=C64::new(0.,-w);
                let a=s*s+s.scale(0.2+damping)+C64::from_re(100.+kc);
                let d=s*s+s.scale(0.4+damping)+C64::from_re(225.+kc);
                let b=-s.scale(damping)-C64::from_re(kc);let determinant=a*d-b*b;
                let x=[(d*f[0]-b*f[1])/determinant,(a*f[1]-b*f[0])/determinant];
                let out=local.model().factor(w,||false).unwrap().response(&f,||false).unwrap();
                for i in 0..2 {near(out.state[2*i],x[i]);near(out.output[i],s*x[i]);}
                let v=[s*x[0],s*x[1]];
                let heat=0.5*(0.2*v[0].norm_sq()+0.4*v[1].norm_sq()+damping*(v[0]-v[1]).norm_sq());
                assert!((out.dissipated_power_w-heat).abs()<1e-10*heat);
                assert!(out.power_residual_w.abs()<1e-11*heat);
            }
        }

        #[test]
        fn g0_only_genuine_two_sided_contact_tangents_enter_the_spectrum() {
            let owner=|alpha|ContactStorage::new(Box::new(QuadraticStorage::new(
                vec![100.,0.,0.,1.],2).unwrap()),1,vec![obstacle(vec![1.],alpha,0.5)]).unwrap();
            let kink=owner(1.);
            assert!(!kink.has_two_sided_contact_tangent(&[0.,0.],&[0.,0.]));
            assert!(kink.has_two_sided_contact_tangent(&[-0.001,0.],&[-0.1,0.]));
            assert!(!kink.has_two_sided_contact_tangent(&[0.01,-2.],&[1.,-2.]));
            assert!(kink.has_two_sided_contact_tangent(&[0.01,-1.99],&[1.,-1.99]));
            assert!(owner(1.5).has_two_sided_contact_tangent(&[0.,0.],&[0.,0.]));
            let system=ImpactSystem::new(vec![body(vec![10.],&[0.],vec![0.2])],
                vec![obstacle(vec![1.],1.,0.)],vec![],vec![],config()).unwrap();
            assert!(system.linearize_equilibrium(EquilibriumPoint {state:system.state(),input:&[0.],
                balance_scales:&[1.;2],relative_tolerance:1e-12},budget(),32,||false).is_err());
            assert_eq!(system.samples(),0);
        }

        #[test]
        fn g0_g4_preparation_and_refused_trials_leave_the_original_trajectory_unchanged() {
            let mut system=one(0.02);let mut untouched=one(0.02);
            let force=[2.+1000.*0.02_f64.powf(1.5)];let before=system.state().to_vec();
            let energy=system.stored_energy_j();let scales=[10.;2];let calls=Cell::new(0);
            let point=EquilibriumPoint {state:&before,input:&force,balance_scales:&scales,relative_tolerance:1e-12};
            system.linearize_equilibrium(point,budget(),32,||{calls.set(calls.get()+1);false}).unwrap();
            let total=calls.get();
            for cut in [1,total/2,total] {
                let count=Cell::new(0);
                assert!(system.linearize_equilibrium(point,budget(),32,||{count.set(count.get()+1);count.get()>=cut}).is_err());
            }
            assert!(system.linearize_equilibrium(point,budget(),1,||false).is_err());
            let wrong=[0.];
            assert!(system.linearize_equilibrium(EquilibriumPoint {input:&wrong,..point},budget(),32,||false).is_err());
            assert_eq!(system.state(),before);assert_eq!(system.samples(),0);
            assert_eq!(system.stored_energy_j().to_bits(),energy.to_bits());
            let gate=CancelGate::new_clock_free();let excitation=[force[0]+1e-4];
            for _ in 0..8 {system.step(&excitation,&gate).unwrap();untouched.step(&excitation,&gate).unwrap();}
            assert_eq!(system.state(),untouched.state());assert_eq!(system.samples(),untouched.samples());
        }
    }
}
