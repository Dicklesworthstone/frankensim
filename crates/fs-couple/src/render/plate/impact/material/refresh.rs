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
