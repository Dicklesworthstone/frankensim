//! A finite-mass translating support carrying both ends of tensioned filaments.
//!
//! Relative string motion and support motion have kinetic cross terms. Retain
//! them by an exact triangular mass normalization, not by moving contact gaps:
//! w_i(s)=X+sum(phi_in(s)*q_in), z_in=q_in+c_in*X, z_C=sqrt(M_r)*X,
//! c_in=integral(mu_i*phi_in ds), M_r=M_C+sum(mu_i L_i)-sum(c_in^2).
//! Then T=|zdot|^2/2. Existing string energy acts on q=z-c*X; its transpose
//! supplies support reactions. The same transformation applies to contact and
//! relative modal damping. No new constitutive law or time integrator.
use std::{ops::Range, rc::Rc};
use fs_dcontact::Obstacle;
use fs_math::det;
use crate::modal_acoustic_time::ModalAcousticState;
use super::{BodyPotential, ImpactBody, ImpactError, MAX_IMPACT_MODES, invalid};
use super::linear::wire::{LineContact, WireSpan};
use super::string::{StringObservation, StringPotential, StringStretching};

/// One physical rail carrying both endpoints of every supplied wire. Positive
/// translation has the same direction as wire displacement; geometry owns it.
#[derive(Clone, Copy, Debug)]
pub struct TranslatingSupport {
    /// Rail/actuator effective mass only [kg]. Wire mass is added exactly once.
    pub mass_kg: f64,
    /// Optional ground stiffness [N/m], with rest at X=0.
    pub stiffness_n_m: f64,
    /// Optional ground resistance [N s/m].
    pub damping_n_s_m: f64,
    pub initial_position_m: f64,
    pub initial_velocity_m_s: f64,
    /// Symmetric |X| limit [m]; refusal, not a stop or clamp.
    pub maximum_travel_m: f64,
    /// Relative wire-slope bound in (0,0.3], including linear wires.
    pub maximum_slope: f64,
}
impl TranslatingSupport {
    pub fn validate(self) -> Result<(), ImpactError> {
        if [self.mass_kg,self.maximum_travel_m,self.maximum_slope].iter().any(|v| !v.is_finite() || *v<=0.0)
            || [self.stiffness_n_m,self.damping_n_s_m].iter().any(|v| !v.is_finite() || *v<0.0)
            || !self.initial_position_m.is_finite() || !self.initial_velocity_m_s.is_finite()
            || self.initial_position_m.abs()>self.maximum_travel_m || self.maximum_slope>0.3 {
            return Err(invalid("moving support requires explicit finite mass, motion, losses and travel/slope bounds"));
        }
        Ok(())
    }
}

/// Aggregate observations, with all energy already included in the host ledger.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SupportObservation {
    pub position_m: f64,
    pub velocity_m_s: f64,
    pub maximum_wire_slope: f64,
    pub maximum_wire_tension_n: f64,
    pub wire_stretching_energy_j: f64,
}

#[derive(Debug)]
struct Data {
    support: TranslatingSupport,
    wires: Vec<WireSpan>,
    potentials: Vec<StringPotential>,
    ranges: Vec<Range<usize>>,
    translation: Vec<f64>,
    damping: Vec<f64>,
    weight: f64,
    total_mass_kg: f64,
    // A Gershgorin upper bound for the coupled SMALL-SIGNAL pencil, not its
    // eigenfrequencies. Only the host's cold Nyquist admission consumes this.
    clock_bounds: Vec<f64>,
}

/// Immutable shared reduction. Relative string coordinates come first, followed
/// by exactly one support coordinate. Only the physical string laws are reused;
/// no additional string mass, tension spring, prescribed motion or reset.
#[derive(Clone, Debug)]
pub struct SupportedStrings(Rc<Data>);
impl SupportedStrings {
    pub fn new(wires: &[WireSpan], stretching: Option<StringStretching>, support: TranslatingSupport)
        -> Result<Self, ImpactError>
    {
        support.validate()?;
        if let Some(s)=stretching {s.validate()?;}
        let count=wires.iter().try_fold(0usize,|n,w|n.checked_add(w.damping_per_s.len()))
            .ok_or_else(||invalid("supported wire mode count overflow"))?;
        if wires.is_empty() || wires.len()>32 || count==0 || count>=MAX_IMPACT_MODES {
            return Err(invalid("supported wires exceed the bounded complete mechanical basis"));
        }
        let law=StringStretching {axial_rigidity_n:stretching.map_or(0.0,|s|s.axial_rigidity_n),
            maximum_slope:stretching.map_or(support.maximum_slope,|s|s.maximum_slope.min(support.maximum_slope))};
        let mut potentials=Vec::with_capacity(wires.len());let mut ranges=Vec::with_capacity(wires.len());
        let mut translation=Vec::with_capacity(count);let mut damping=Vec::with_capacity(count);
        let mut total_mass=support.mass_kg;let mut residual_mass=support.mass_kg;let mut frequencies=Vec::with_capacity(count);
        for wire in wires {
            let n=wire.damping_per_s.len();
            let body=wire.stretching_body(vec![ModalAcousticState::default();n],law)?;
            let BodyPotential::String(potential)=body.potential else {return Err(invalid("wire preparation changed its storage image"));};
            let mass=wire.linear_density_kg_m*wire.length_m();
            let start=translation.len();let mut represented=0.0;
            for k in 1..=n {
                let c=if k%2==0 {0.0}else{2.0*det::sqrt(2.0*mass)/(k as f64*core::f64::consts::PI)};
                represented+=c*c;translation.push(c);
            }
            let remainder=mass-represented;
            if !mass.is_finite() || mass<=0.0 || !remainder.is_finite() || remainder<=0.0 {
                return Err(invalid("supported wire translation mass is unrepresentable"));
            }
            // The unresolved sine tail still translates with the support.
            // Discarding it would change total rigid-translation kinetic energy.
            residual_mass+=remainder;total_mass+=mass;ranges.push(start..start+n);
            frequencies.extend_from_slice(potential.omegas());potentials.push(potential);
            damping.extend_from_slice(&wire.damping_per_s);
        }
        let weight=1.0/det::sqrt(residual_mass);
        if !total_mass.is_finite() || !weight.is_finite() || weight<=0.0 {
            return Err(invalid("moving support mass normalization overflow"));
        }
        let mut support_row=support.stiffness_n_m*weight*weight;let mut bound=0.0_f64;
        for (&omega,&c) in frequencies.iter().zip(&translation) {
            let a=c*weight;let k=omega*omega;
            bound=bound.max(k*(1.0+a.abs()));support_row+=k*(a*a+a.abs());
        }
        bound=bound.max(support_row);
        if !bound.is_finite() {return Err(invalid("moving support linear clock bound overflow"));}
        let result=Self(Rc::new(Data {support,wires:wires.to_vec(),potentials,ranges,translation,damping,weight,
            total_mass_kg:total_mass,clock_bounds:vec![det::sqrt(bound);count+1]}));
        // Check initial energy/motion representation before returning an adapter.
        result.observe_interleaved(&result.initial_state(),0)?;
        Ok(result)
    }
    pub fn mode_count(&self)->usize {self.0.translation.len()+1}
    pub fn wire_count(&self)->usize {self.0.wires.len()}
    pub fn support_coordinate(&self)->usize {self.mode_count()-1}
    /// Generalized force is F * this weight, applied ONLY at support_coordinate.
    pub fn force_weight(&self)->f64 {self.0.weight}
    pub fn total_mass_kg(&self)->f64 {self.0.total_mass_kg}
    pub fn clock_bounds(&self)->&[f64] {&self.0.clock_bounds}
    fn initial_state(&self)->Vec<f64> {
        let mut state=vec![0.0;2*self.mode_count()];let s=self.0.support;
        for (i,&c) in self.0.translation.iter().enumerate() {
            state[2*i]=c*s.initial_position_m;state[2*i+1]=c*s.initial_velocity_m_s;
        }
        let i=self.support_coordinate();state[2*i]=s.initial_position_m/self.0.weight;
        state[2*i+1]=s.initial_velocity_m_s/self.0.weight;state
    }
    /// Initial relative modes are at rest; the entire support and all wires have
    /// the declared rigid translation/velocity. Internal damping is transformed
    /// by the host; adding diagonal modal drag here would double-count it.
    pub fn body(&self)->ImpactBody {
        let x=self.initial_state();
        ImpactBody {potential:BodyPotential::SupportedStrings(self.clone()),
            initial:x.chunks_exact(2).map(|s|ModalAcousticState {
                displacement_m_sqrt_kg:s[0],velocity_m_sqrt_kg_per_s:s[1]}).collect(),
            damping_per_s:vec![0.0;self.mode_count()]}
    }
    fn relative(&self,z:&[f64],q:&mut[f64])->bool {
        if z.len()!=self.mode_count() || q.len()!=self.support_coordinate()
            || z.iter().any(|v|!v.is_finite()) {return false;}
        let x=z[self.support_coordinate()]*self.0.weight;
        for ((q,z),c) in q.iter_mut().zip(z).zip(&self.0.translation) {*q=*z-c*x;}
        q.iter().all(|v|v.is_finite())
    }
    pub(super) fn potential(&self,z:&[f64])->f64 {
        let mut q=[0.0;MAX_IMPACT_MODES];let n=self.support_coordinate();
        if !self.relative(z,&mut q[..n]) {return f64::NAN;}
        let x=z[n]*self.0.weight;
        let mut energy=0.5*self.0.support.stiffness_n_m*x*x;
        for (p,r) in self.0.potentials.iter().zip(&self.0.ranges) {energy+=p.potential(&q[r.clone()]);}
        energy
    }
    pub(super) fn gradient(&self,z:&[f64],out:&mut[f64]) {
        let mut q=[0.0;MAX_IMPACT_MODES];let n=self.support_coordinate();
        if out.len()!=self.mode_count() || !self.relative(z,&mut q[..n]) {out.fill(f64::NAN);return;}
        for (p,r) in self.0.potentials.iter().zip(&self.0.ranges) {p.gradient(&q[r.clone()],&mut out[r.clone()]);}
        out[n]=self.0.support.stiffness_n_m*z[n]*self.0.weight*self.0.weight
            -self.0.weight*self.0.translation.iter().zip(&out[..n]).map(|(c,g)|c*g).sum::<f64>();
    }
    pub(super) fn hessian_vector(&self,z:&[f64],d:&[f64],out:&mut[f64]) {
        let (mut q,mut v)=([0.0;MAX_IMPACT_MODES],[0.0;MAX_IMPACT_MODES]);let n=self.support_coordinate();
        if out.len()!=self.mode_count() || !self.relative(z,&mut q[..n]) || !self.relative(d,&mut v[..n]) {
            out.fill(f64::NAN);return;
        }
        for (p,r) in self.0.potentials.iter().zip(&self.0.ranges) {p.hessian_vector(&q[r.clone()],&v[r.clone()],&mut out[r.clone()]);}
        out[n]=self.0.support.stiffness_n_m*d[n]*self.0.weight*self.0.weight
            -self.0.weight*self.0.translation.iter().zip(&out[..n]).map(|(c,g)|c*g).sum::<f64>();
    }
    /// Preserve the existing line-contact law and receiver row. Its relative
    /// wire row is transformed to ABSOLUTE motion before the implicit solve.
    #[allow(clippy::too_many_arguments)]
    pub fn contact(&self,wire:usize,line:&LineContact,receiver_shapes:&[Vec<f64>],
        receiver:Range<usize>,first:usize,total:usize)->Result<Obstacle,ImpactError>
    {
        let r=self.0.ranges.get(wire).ok_or_else(||invalid("unknown supported wire"))?;
        let end=first.checked_add(self.mode_count()).ok_or_else(||invalid("supported wire address overflow"))?;
        if end>total || receiver.start<end && first<receiver.end {
            return Err(invalid("supported wire and receiver ranges must be disjoint and complete"));
        }
        let ob=self.0.wires[wire].contact(line,receiver_shapes,receiver,first+r.start..first+r.end,total)?;
        let mut rows=ob.collocation().to_vec();
        for row in rows.chunks_exact_mut(total) {
            let cross=row[first+r.start..first+r.end].iter().zip(&self.0.translation[r.clone()])
                .map(|(b,c)|b*c).sum::<f64>();
            row[first+self.support_coordinate()]=(-1.0-cross)*self.0.weight;
        }
        Obstacle::new(rows,ob.n_points(),total,ob.gaps().to_vec(),ob.weights().to_vec(),
            ob.stiffness(),ob.alpha(),ob.provenance().into())
            .and_then(|o|o.with_internal_loss(ob.internal_loss())).map_err(|e|ImpactError::Owner(e.to_string()))
    }
    // Internal relative modal loss is part of the body, not a collection of
    // user-added localized mufflers. Do not consume/increase that separate cap.
    pub(super) fn add_resistance(&self,first:usize,dim:usize,r:&mut[f64])->Result<(),ImpactError> {
        let n=self.support_coordinate();let p=2*(first+n)+1;
        if p>=dim || r.len()!=dim*dim {return Err(invalid("supported damping dimensions"));}
        for i in 0..n {
            let j=2*(first+i)+1;let a=self.0.translation[i]*self.0.weight;let drag=self.0.damping[i];
            r[j*dim+j]+=drag;r[j*dim+p]-=drag*a;r[p*dim+j]-=drag*a;r[p*dim+p]+=drag*a*a;
        }
        r[p*dim+p]+=self.0.support.damping_n_s_m*self.0.weight*self.0.weight;
        if r.iter().any(|v|!v.is_finite()) {return Err(invalid("supported damping overflow"));}Ok(())
    }
    pub fn string_observation(&self,wire:usize,z:&[f64])->Result<StringObservation,ImpactError> {
        let r=self.0.ranges.get(wire).ok_or_else(||invalid("unknown supported wire"))?;
        let mut q=[0.0;MAX_IMPACT_MODES];
        if !self.relative(z,&mut q[..self.support_coordinate()]) {return Err(invalid("invalid supported wire state"));}
        self.0.potentials[wire].observe(&q[r.clone()])
    }
    pub fn observe_interleaved(&self,state:&[f64],first:usize)->Result<SupportObservation,ImpactError> {
        let end=first.checked_add(self.mode_count()).and_then(|n|n.checked_mul(2))
            .ok_or_else(||invalid("supported observation address overflow"))?;
        if state.len()<end || state[2*first..end].iter().any(|v|!v.is_finite()) {
            return Err(invalid("supported observation requires complete finite motion"));
        }
        let n=self.support_coordinate();let x=state[2*(first+n)]*self.0.weight;
        if x.abs()>self.0.support.maximum_travel_m {return Err(invalid("support travel exceeds its declared validity"));}
        let mut q=[0.0;MAX_IMPACT_MODES];
        for i in 0..n {q[i]=state[2*(first+i)]-self.0.translation[i]*x;}
        let mut result=SupportObservation {position_m:x,velocity_m_s:state[2*(first+n)+1]*self.0.weight,
            ..SupportObservation::default()};
        for (p,r) in self.0.potentials.iter().zip(&self.0.ranges) {
            let observation=p.observe(&q[r.clone()])?;
            result.maximum_wire_slope=result.maximum_wire_slope.max(observation.slope_bound);
            result.maximum_wire_tension_n=result.maximum_wire_tension_n.max(observation.tension_n);
            result.wire_stretching_energy_j+=observation.stretching_energy_j;
        }
        if !result.wire_stretching_energy_j.is_finite() {return Err(invalid("supported wire energy overflow"));}
        Ok(result)
    }
}

#[cfg(test)]
#[path = "supported_tests.rs"]
mod tests;
