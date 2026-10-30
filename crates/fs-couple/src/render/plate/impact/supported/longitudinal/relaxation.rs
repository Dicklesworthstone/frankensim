//! Hereditary strain on the original moving support, with no second time owner.
//! Each physical layer/arm has independent mean and retained-cosine memories.
//! The material law, exact support inertia pullback and endpoint chord remain
//! the same whether excitation is a force, hammer contact or initial pluck.
use super::*;
use crate::render::plate::impact::relaxation::InitialMemory;
use fs_phs::{MappedRelaxationBranch, RelaxationStrainMap};

/// One positive relaxing axial rigidity from a physical layer and its time.
/// Geometry and equilibrium EA remain in the original supported string.
#[derive(Clone, Debug, PartialEq)]
pub struct SupportedRelaxationArm {
    pub string: usize,
    pub axial_rigidity_n: f64,
    pub relaxation_time_s: f64,
}

#[derive(Clone, Debug)]
pub(in crate::render::plate::impact) struct SupportedRelaxation {
    arms: Vec<SupportedRelaxationArm>,
    pub(in crate::render::plate::impact) initial: InitialMemory,
    states: usize,
    stiffness_bound: f64,
    // (complete source omega ceiling, exact source step). This is retained on
    // copied bodies so generic impact composition cannot discard authority.
    source_domain: Option<(f64,f64)>,
}

impl SupportedStrings {
    /// Attach independent full-strain Maxwell arms as immutable material data.
    /// Every ImpactSystem made from this body automatically retains their state;
    /// copying bodies into a contact assembly cannot silently drop the law.
    /// The original equilibrium stiffness and preload are never replaced.
    ///
    /// History coefficients obey the SAME explicit strain envelope as the
    /// retained strain: |mean|+sum|cosine| <= maximum_strain. The original
    /// discrete-gradient owner refuses an endpoint outside that domain. Its
    /// nonlinear update is not asserted to preserve a convex strain hull.
    /// The complete hereditary tension must remain positive over the span.
    /// A later radiation attachment rechecks the combined oscillator/coupling
    /// bound. Nonlinear squeeze-film loading needs an additional strain/history
    /// envelope and is currently refused by the original impact owner.
    pub fn with_longitudinal_relaxation(mut self, arms: Vec<SupportedRelaxationArm>,
        initial: InitialMemory, maximum_states: usize) -> Result<Self, ImpactError>
    {
        if self.0.axial_relaxation.is_some() || !(1..=256).contains(&maximum_states) {
            return Err(invalid("supported material memory needs one cold attachment and a 1..=256 state ceiling"));
        }
        let mut states=0usize;let mut bound=0.0;
        for arm in &arms {
            let Some(SupportedPotential::Longitudinal(p))=self.0.potentials.get(arm.string) else {
                return Err(invalid("supported memory requires the original longitudinal strain image"));
            };
            if !arm.axial_rigidity_n.is_finite() || arm.axial_rigidity_n<=0.0
                || !arm.relaxation_time_s.is_finite() || arm.relaxation_time_s<=0.0 {
                return Err(invalid("supported memory needs positive finite relaxing rigidity and time"));
            }
            states=states.checked_add(1+p.model.axial_modes().len())
                .filter(|n|*n<=maximum_states)
                .ok_or_else(||invalid("complete supported material histories exceed their declared ceiling"))?;
            for channel in 0..=p.model.axial_modes().len() {
                let map=SupportedStrainMap {support:self.clone(),string:arm.string,channel,
                    offset:0,base_dim:2*self.mode_count()};
                let c=arm.axial_rigidity_n*self.0.wires[arm.string].length_m()
                    *if channel==0 {1.0} else {0.5};
                if !c.is_finite() || c<=0.0 {return Err(invalid("supported memory energy normalization overflows"));}
                // Include the geometric stress Hessian as well as B^T C B.
                // Both epsilon and admitted viscous strain are bounded here.
                bound+=c*map.stiffness_bound(2.0*p.maximum_strain)?;
            }
        }
        if !bound.is_finite() {return Err(invalid("supported hereditary stiffness bound overflows"));}
        let data=Rc::make_mut(&mut self.0);
        if bound>0.0 {
            for omega in &mut data.clock_bounds {*omega=det::sqrt(*omega**omega+bound);}
            if data.clock_bounds.iter().any(|w|!w.is_finite()) {return Err(invalid("supported memory clock bound overflows"));}
        }
        data.axial_relaxation=Some(SupportedRelaxation {
            arms,initial,states,stiffness_bound:bound,source_domain:None,
        });
        Ok(self)
    }
    /// Added full-strain material tangent envelope [s^-2], including retained
    /// hereditary stress and the original canonical support-coordinate map.
    pub fn longitudinal_relaxation_stiffness_bound_per_s2(&self) -> f64 {
        self.0.axial_relaxation.as_ref().map_or(0.0,|m|m.stiffness_bound)
    }
    pub(in crate::render::plate::impact) fn memory_configuration(&self) -> Option<&SupportedRelaxation> {
        self.0.axial_relaxation.as_ref()
    }
    pub(crate) fn bind_longitudinal_relaxation_domain(mut self,omega_ceiling:f64,dt_s:f64)
        -> Result<Self,ImpactError>
    {
        if omega_ceiling.is_nan() || omega_ceiling<=0.0 || !dt_s.is_finite() || dt_s<=0.0 {
            return Err(invalid("supported source memory needs positive frequency support and a finite physical clock"));
        }
        let Some(memory)=Rc::make_mut(&mut self.0).axial_relaxation.as_mut() else {
            return Err(invalid("cannot bind source authority without original material memory"));
        };
        if memory.source_domain.is_some() {return Err(invalid("supported material authority cannot be replaced"));}
        memory.source_domain=Some((omega_ceiling,dt_s));Ok(self)
    }
    pub(in crate::render::plate::impact) fn memory_branches(&self, offset:usize,base_dim:usize)
        -> Result<Vec<MappedRelaxationBranch>,ImpactError>
    {
        let Some(config)=&self.0.axial_relaxation else {return Ok(vec![]);};
        if 2*(offset+self.mode_count())>base_dim {return Err(invalid("supported memory lost its original state prefix"));}
        let mut branches=Vec::with_capacity(config.states);
        for arm in &config.arms {
            let SupportedPotential::Longitudinal(p)=&self.0.potentials[arm.string] else {unreachable!()};
            for channel in 0..=p.model.axial_modes().len() {
                branches.push(MappedRelaxationBranch {strain:Rc::new(SupportedStrainMap {
                    support:self.clone(),string:arm.string,channel,offset,base_dim,
                }),stiffness:arm.axial_rigidity_n*self.0.wires[arm.string].length_m()
                    *if channel==0 {1.0} else {0.5},relaxation_time_s:arm.relaxation_time_s});
            }
        }
        Ok(branches)
    }
    pub(in crate::render::plate::impact) fn admit_memory(&self, offset:usize,state:&[f64],
        memory:&[f64]) -> Result<(),ImpactError>
    {
        let Some(config)=&self.0.axial_relaxation else {return Ok(());};
        if memory.len()!=config.states || memory.iter().any(|v|!v.is_finite()) {
            return Err(invalid("supported material history requires every finite original coordinate"));
        }
        let mut z=vec![0.0;self.mode_count()];
        for (i,z) in z.iter_mut().enumerate() {
            *z=*state.get(2*(offset+i)).ok_or_else(||invalid("supported history state prefix is incomplete"))?;
        }
        let mut q=vec![0.0;self.support_coordinate()];
        if !self.relative(&z,&mut q) {return Err(invalid("supported history has invalid original coordinates"));}
        let mut strains=Vec::with_capacity(self.wire_count());let mut tensions=Vec::new();
        for (i,p) in self.0.potentials.iter().enumerate() {
            let SupportedPotential::Longitudinal(p)=p else {strains.push(vec![]);tensions.push(vec![]);continue;};
            let mut local=[0.0;MAX_IMPACT_MODES];
            let count=self.local_coordinates(i,&q,z[self.support_coordinate()]*self.0.weight,&mut local);
            p.observe(&local[..count])?;
            let strain=p.strain(&local[..count])?;
            let mut epsilon=vec![strain.mean];epsilon.extend(strain.harmonics);
            let mut tension:Vec<_>=epsilon.iter().map(|v|p.ea*v).collect();tension[0]+=p.tension;
            strains.push(epsilon);tensions.push(tension);
        }
        let mut at=0;
        for arm in &config.arms {
            let SupportedPotential::Longitudinal(p)=&self.0.potentials[arm.string] else {unreachable!()};
            let epsilon=&strains[arm.string];let length=self.0.wires[arm.string].length_m();
            let mut history_bound=0.0;
            for (channel,&epsilon) in epsilon.iter().enumerate() {
                let c=arm.axial_rigidity_n*length*if channel==0 {1.0} else {0.5};
                let viscous=memory[at]/det::sqrt(c);at+=1;history_bound+=viscous.abs();
                tensions[arm.string][channel]+=arm.axial_rigidity_n*(epsilon-viscous);
            }
            if !history_bound.is_finite() || history_bound>p.maximum_strain {
                return Err(invalid("supported material memory exceeds its declared strain-history envelope"));
            }
        }
        if tensions.iter().filter(|t|!t.is_empty()).any(|t|t.iter().any(|v|!v.is_finite())
            || t[0]-t[1..].iter().map(|v|v.abs()).sum::<f64>()<=0.0) {
            return Err(invalid("complete hereditary support stress loses local string tension"));
        }
        Ok(())
    }
}
impl SupportedRelaxation {
    pub(in crate::render::plate::impact) fn has_source_domain(&self)->bool {self.source_domain.is_some()}
    pub(in crate::render::plate::impact) fn admit_domain(&self,omega:f64,dt:f64)->Result<(),ImpactError> {
        if self.source_domain.is_some_and(|(hi,step)|omega>hi || dt.to_bits()!=step.to_bits()) {
            return Err(invalid("composed impact object exceeds original supported material frequency/clock authority"));
        }
        Ok(())
    }
}

/// Geometry only. The actual z history is held once in ImpactSystem's state.
#[derive(Clone, Debug)]
struct SupportedStrainMap {
    support: SupportedStrings,
    string: usize,
    channel: usize,
    offset: usize,
    base_dim: usize,
}
impl SupportedStrainMap {
    fn potential(&self) -> &LongitudinalPotential {
        let SupportedPotential::Longitudinal(p)=&self.support.0.potentials[self.string] else {unreachable!()};p
    }
    fn local(&self,state:&[f64]) -> Option<Vec<f64>> {
        if state.len()!=self.base_dim {return None;}
        let s=&self.support;let mut z=vec![0.0;s.mode_count()];
        for (i,z) in z.iter_mut().enumerate() {*z=state[2*(self.offset+i)];}
        let mut q=vec![0.0;s.support_coordinate()];if !s.relative(&z,&mut q) {return None;}
        let mut local=[0.0;MAX_IMPACT_MODES];
        let count=s.local_coordinates(self.string,&q,z[s.support_coordinate()]*s.0.weight,&mut local);
        Some(local[..count].to_vec())
    }
    fn pull(&self,local:&[f64],out:&mut[f64]) {
        out.fill(0.0);let s=&self.support;let range=&s.0.storage_ranges[self.string];
        for (i,&g) in range.clone().zip(local) {out[2*(self.offset+i)]=g;}
        // Equal/opposite endpoint work includes the relative-mode inertia lift.
        let endpoint=local[range.len()]-range.clone().zip(local)
            .map(|(i,g)|s.0.translation[i]*g).sum::<f64>();
        out[2*(self.offset+s.support_coordinate())]=s.0.weight*endpoint;
    }
    fn local_gradient(&self,q:&[f64],hessian:bool) -> Option<Vec<f64>> {
        let p=self.potential();let n=p.scales.len();let mut weights=vec![0.0;p.model.axial_modes().len()];
        let mean=if self.channel==0 {1.0} else {weights[self.channel-1]=1.0;0.0};
        let slopes=p.slopes(q);let mut g=vec![0.0;slopes.len()];
        p.model.geometric_strain_gradient(&slopes,mean,&weights,&mut g).ok()?;
        let mut out=vec![0.0;p.count()];
        for plane in 0..2 {for mode in 0..n {
            out[plane*n+mode]=p.scales[mode]*g[plane*(n+1)+mode+1];
        }}
        let last=out.len()-1;out[last]=p.chord[0]*g[0]+p.chord[1]*g[n+1];
        if !hessian {
            out[last]+=p.chord[2]*mean;
            for (i,(&scale,&weight)) in p.model.axial_slope_scales().iter().zip(&weights).enumerate() {
                out[p.omegas.len()+i]=scale*weight;
            }
        }
        Some(out)
    }
    fn stiffness_bound(&self,residual_bound:f64) -> Result<f64,ImpactError> {
        let p=self.potential();let mut gradient_square=0.0;let mut rows=vec![0.0;self.base_dim];
        for mode in 0..self.support.mode_count() {
            let mut unit=vec![0.0;self.base_dim];unit[2*(self.offset+mode)]=1.0;
            let local=self.local(&unit).ok_or_else(||invalid("memory tangent map lost its original coordinates"))?;
            let slopes=p.slopes(&local);let n=p.scales.len()+1;
            let slope_bound=(0..n).map(|i|det::hypot(slopes[i],slopes[n+i])).sum::<f64>();
            let axial_bound=(p.chord[2]*local[local.len()-1]).abs()+p.model.axial_slope_scales().iter()
                .zip(&local[p.omegas.len()..local.len()-1]).map(|(a,q)|(a*q).abs()).sum::<f64>();
            let bound=(if self.channel==0 {1.0} else {2.0})*(axial_bound+p.maximum_slope*slope_bound);
            gradient_square+=bound*bound;
            let mut column=vec![0.0;self.base_dim];
            self.pull(&self.local_gradient(&local,true).ok_or_else(||invalid("memory Hessian map failed"))?,&mut column);
            for (r,c) in rows.iter_mut().zip(column) {*r+=c.abs();}
        }
        let result=gradient_square+residual_bound*rows.into_iter().fold(0.0_f64,f64::max);
        if !result.is_finite() {return Err(invalid("supported memory tangent envelope overflows"));}Ok(result)
    }
}
impl RelaxationStrainMap for SupportedStrainMap {
    fn is_quadratic_strain(&self)->bool {true}
    fn state_dimension(&self)->usize {self.base_dim}
    fn strain(&self,state:&[f64])->f64 {
        let Some(q)=self.local(state) else {return f64::NAN;};
        let Ok(s)=self.potential().strain(&q) else {return f64::NAN;};
        if self.channel==0 {s.mean} else {s.harmonics[self.channel-1]}
    }
    fn gradient(&self,state:&[f64],out:&mut[f64]) {
        if out.len()!=self.base_dim {out.fill(f64::NAN);return;}
        let Some(q)=self.local(state) else {out.fill(f64::NAN);return;};
        let Some(g)=self.local_gradient(&q,false) else {out.fill(f64::NAN);return;};self.pull(&g,out);
    }
    fn hessian_vector(&self,_state:&[f64],direction:&[f64],out:&mut[f64]) {
        if out.len()!=self.base_dim {out.fill(f64::NAN);return;}
        let Some(d)=self.local(direction) else {out.fill(f64::NAN);return;};
        let Some(h)=self.local_gradient(&d,true) else {out.fill(f64::NAN);return;};self.pull(&h,out);
    }
}

#[cfg(test)]
mod tests;
