//! Per-body linearization reuse, bound to the very same immutable storage owner.
//! Only shell/film potential Hessians are retained. Felt, contact, relaxation,
//! gas-film and radiation terms remain in the original live composed action.
use std::rc::Rc;
use super::{BodyPotential,MechanicalStorage};

struct BodyMatrix {
    point:Vec<f64>,
    matrix:Vec<f64>,
    valid:bool,
}

pub(in crate::render::plate::impact) struct BodyTangents {
    owner:Rc<MechanicalStorage>,
    entries:Vec<Option<BodyMatrix>>,
    scratch:Vec<[f64;3]>,
    assemblies:usize,
    actions:usize,
}

impl BodyTangents {
    // The admitted host bounds the sum of body modes by MAX_IMPACT_MODES.
    // Therefore sum(n_i^2) <= MAX_IMPACT_MODES^2: no acoustic memory matrix is
    // duplicated. This constructor is cold, like the original solver prepare.
    pub(in crate::render::plate::impact) fn new(owner:Rc<MechanicalStorage>)->Self {
        let mut maximum=0usize;
        let entries=owner.bodies.iter().map(|body| {
            if matches!(body,BodyPotential::Shell(_)|BodyPotential::Membrane(_)) {
                let n=body.count();maximum=maximum.max(n);
                Some(BodyMatrix {point:vec![0.0;n],matrix:vec![0.0;n*n],valid:false})
            }else{None}
        }).collect();
        Self {owner,entries,scratch:vec![[0.0;3];maximum],assemblies:0,actions:0}
    }
    pub(in crate::render::plate::impact) fn belongs_to(&self,owner:&Rc<MechanicalStorage>)->bool {
        Rc::ptr_eq(&self.owner,owner)
    }
    pub(in crate::render::plate::impact) fn reset_counts(&mut self) {self.assemblies=0;self.actions=0;}
    pub(in crate::render::plate::impact) fn counts(&self)->(usize,usize) {(self.assemblies,self.actions)}

    pub(super) fn apply(&mut self,index:usize,q:&[f64],d:&[f64],out:&mut[f64])->bool {
        let Some(body)=self.owner.bodies.get(index) else {return false;};
        let n=body.count();
        if q.len()!=n || d.len()!=n || out.len()!=n || q.iter().chain(d).any(|v|!v.is_finite()) {
            return false;
        }
        let Some(entry)=self.entries[index].as_mut() else {
            body.hessian_vector(q,d,out);
            return out.iter().all(|v|v.is_finite());
        };
        if d.iter().all(|v|*v==0.0) {out.fill(0.0);return true;}
        // No approximate key, hashing, previous-step extrapolation or threshold.
        // A rejected trial may populate numerical scratch, never physical state.
        if !entry.valid || entry.point.iter().zip(q).any(|(a,b)|a.to_bits()!=b.to_bits()) {
            entry.valid=false;
            self.assemblies=self.assemblies.saturating_add(1);
            match body {
                BodyPotential::Shell(s)=>s.hessian_matrix(q,&mut entry.matrix,&mut self.scratch[..n]),
                BodyPotential::Membrane(s)=>s.reduction().hessian_matrix(q,&mut entry.matrix,&mut self.scratch[..n]),
                _=>return false, // Layout and owner are immutable after preparation.
            }
            if entry.matrix.iter().any(|v|!v.is_finite()) {return false;}
            entry.point.copy_from_slice(q);entry.valid=true;
        }
        self.actions=self.actions.saturating_add(1);
        out.fill(0.0);
        // Analytic Newton supplies basis columns. Do not turn a single-column
        // lookup back into a dense product by multiplying all its exact zeros.
        // General directions retain increasing-column accumulation order.
        for (col,&value) in d.iter().enumerate() {
            if value==0.0 {continue;}
            for row in 0..n {out[row]+=entry.matrix[row*n+col]*value;}
        }
        out.iter().all(|v|v.is_finite())
    }
}
