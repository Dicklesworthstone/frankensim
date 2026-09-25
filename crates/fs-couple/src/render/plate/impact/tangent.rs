//! Exact storage tangents for the composed percussion model, no new physics.
use std::rc::Rc;
use fs_material::fiber::Uniaxial;
use fs_phs::Storage;
use super::{BodyPotential, MechanicalStorage, ImpactSystem, MAX_IMPACT_MODES};

// Storage and its derivative share the SAME immutable geometry/contact owner
// and RefCell history. No duplicate shell mesh, mode basis or material state.
pub(super) struct SharedStorage<T>(pub Rc<T>);
impl<T: Storage> Storage for SharedStorage<T> {
    fn hamiltonian(&self,x:&[f64])->f64 { self.0.hamiltonian(x) }
    fn gradient(&self,x:&[f64],out:&mut[f64]) { self.0.gradient(x,out); }
}
impl BodyPotential {
    fn hessian_vector(&self,q:&[f64],d:&[f64],out:&mut[f64]) {
        match self {
            Self::Shell(s)=>s.hessian_vector(q,d,out),
            Self::Membrane(s)=>s.reduction().hessian_vector(q,d,out),
            Self::String(s)=>s.hessian_vector(q,d,out),
            Self::SupportedStrings(s)=>s.hessian_vector(q,d,out),
            Self::Linear(w)=>{for ((o,w),d) in out.iter_mut().zip(w).zip(d) {*o=w*w*d;}}
        }
    }
}
impl MechanicalStorage {
    fn hessian_vector(&self,x:&[f64],direction:&[f64],out:&mut[f64])->bool {
        let expected=2*self.modes+self.pads.iter().map(|p|p.spec.creep.len()).sum::<usize>();
        if x.len()!=expected || direction.len()!=expected || out.len()!=expected
            || x.iter().chain(direction).any(|v|!v.is_finite()) {return false;}
        out.fill(0.0);
        let (mut q,mut d,mut hd)=([0.0;MAX_IMPACT_MODES],[0.0;MAX_IMPACT_MODES],[0.0;MAX_IMPACT_MODES]);
        for i in 0..self.modes {q[i]=x[2*i];d[i]=direction[2*i];out[2*i+1]=direction[2*i+1];}
        // Mechanical potentials do not depend on momentum. Skip whole bodies
        // with a zero directional displacement, not tiny nonzero couplings.
        let mut offset=0;
        for body in &self.bodies {
            let end=offset+body.count();
            if d[offset..end].iter().any(|d|*d!=0.0) {
                body.hessian_vector(&q[offset..end],&d[offset..end],&mut hd[offset..end]);
            }
            offset=end;
        }
        for v in &self.volumes {
            let rate=v.areas.iter().zip(&d).map(|(b,d)|b*d).sum::<f64>();
            let force=(v.bulk_modulus_pa/v.volume_m3)*rate;
            for (i,b) in v.areas.iter().enumerate() {hd[i]+=b*force;}
        }
        let history=self.histories.borrow();
        for (pad,h) in self.pads.iter().zip(history.iter()) {
            let mut dc=pad.spec.weights.iter().zip(&d).map(|(b,d)|b*d).sum::<f64>();
            for (i,branch) in pad.spec.creep.iter().enumerate() {
                dc-=direction[2*self.modes+pad.creep_start+i]/branch.stiffness_n_m.sqrt();
            }
            // fs-material owns both loading and conditioned unloading tangent.
            let df=pad.spec.area_m2/pad.spec.thickness_m
                *pad.spec.law.tangent(pad.strain(x,self.modes),h)*dc;
            for (i,b) in pad.spec.weights.iter().enumerate() {hd[i]+=df*b;}
            for (i,branch) in pad.spec.creep.iter().enumerate() {
                let index=2*self.modes+pad.creep_start+i;
                out[index]=direction[index]-df/branch.stiffness_n_m.sqrt();
            }
        }
        for i in 0..self.modes {out[2*i]=hd[i];}
        out.iter().all(|v|v.is_finite())
    }
}
impl ImpactSystem {
    pub(super) fn hessian_vector(&self,x:&[f64],d:&[f64],out:&mut[f64])->bool {
        let base=self.relaxation.as_ref().map_or(self.x.len(),|memory|memory.base_dim)
            .min(self.radiation.as_ref().map_or(self.x.len(),|air|air.base_dim))
            .min(self.gas_film.as_ref().map_or(self.x.len(),|gas|gas.base_dim));
        if x.len()!=self.x.len() || d.len()!=x.len() || out.len()!=x.len()
            || x.iter().chain(d).any(|v|!v.is_finite()) {return false;}
        if !self.contact.hessian_vector_with(&x[..base],&d[..base],&mut out[..base],
            |x,d,out|self.mechanical.hessian_vector(x,d,out)) {return false;}
        out[base..].fill(0.0);
        if let Some(memory)=&self.relaxation {memory.add_hessian(d,out);}
        if let Some(air)=&self.radiation {air.add_hessian(d,out);}
        if self.gas_film.as_ref().is_some_and(|gas|!gas.add_hessian(x,d,out)) {return false;}
        out.iter().all(|v|v.is_finite())
    }
}

#[cfg(test)]
mod tests;
