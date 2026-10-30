//! Work-conjugate material strain-rate loss on the original moving support.
//! The retained strain is quadratic, so its midpoint derivative gives the exact
//! strain increment. fs-phs supplies that same accepted coordinate velocity and
//! owns the resulting positive loss, without a separate damping update.
use super::*;

impl SupportedStrings {
    /// Whether positive finite-strain resistance requires the nonlinear owner.
    /// An explicitly zero source keeps its receipt and the elastic solver path.
    pub fn has_finite_strain_axial_viscosity(&self) -> bool {
        self.0.finite_strain_axial_viscosity
            && self.0.axial_viscosity_n_s.as_ref().is_some_and(|v|v.iter().any(|r|*r>0.0))
    }

    pub(in crate::render::plate::impact) fn viscous_reaction(&self,z:&[f64],v:&[f64],out:&mut[f64]) -> bool {
        self.viscous_action(z,v,None,out)
    }

    pub(in crate::render::plate::impact) fn viscous_reaction_tangent(&self,z:&[f64],v:&[f64],
        dz:&[f64],dv:&[f64],out:&mut[f64]) -> bool
    {
        self.viscous_action(z,v,Some((dz,dv)),out)
    }

    fn viscous_action(&self,z:&[f64],v:&[f64],direction:Option<(&[f64],&[f64])>,out:&mut[f64]) -> bool {
        let n=self.support_coordinate();
        let (mut q,mut velocity)=([0.0;MAX_IMPACT_MODES],[0.0;MAX_IMPACT_MODES]);
        if out.len()!=self.mode_count() || !self.relative(z,&mut q[..n])
            || !self.relative(v,&mut velocity[..n]) {return false;}
        out.fill(0.0);
        if !self.has_finite_strain_axial_viscosity() {return true;}
        let (mut dq,mut dvelocity)=([0.0;MAX_IMPACT_MODES],[0.0;MAX_IMPACT_MODES]);
        if let Some((dz,dv))=direction {
            if !self.relative(dz,&mut dq[..n]) || !self.relative(dv,&mut dvelocity[..n]) {return false;}
        }
        let (mut local,mut speed,mut force)=([0.0;MAX_IMPACT_MODES],[0.0;MAX_IMPACT_MODES],[0.0;MAX_IMPACT_MODES]);
        let (mut dlocal,mut dspeed)=([0.0;MAX_IMPACT_MODES],[0.0;MAX_IMPACT_MODES]);
        let mut endpoint=0.0;
        let Some(rigidities)=&self.0.axial_viscosity_n_s else {return false;};
        for (i,(&rigidity,r)) in rigidities.iter().zip(&self.0.storage_ranges).enumerate() {
            if rigidity==0.0 {continue;}
            let SupportedPotential::Longitudinal(potential)=&self.0.potentials[i] else {return false;};
            let count=self.local_coordinates(i,&q,z[n]*self.0.weight,&mut local);
            self.local_coordinates(i,&velocity,v[n]*self.0.weight,&mut speed);
            let local_direction=if let Some((dz,dv))=direction {
                self.local_coordinates(i,&dq,dz[n]*self.0.weight,&mut dlocal);
                self.local_coordinates(i,&dvelocity,dv[n]*self.0.weight,&mut dspeed);
                Some((&dlocal[..count],&dspeed[..count]))
            } else {None};
            if !potential.viscous_action(&local[..count],&speed[..count],local_direction,
                rigidity,&mut force[..count]) {return false;}
            out[r.clone()].copy_from_slice(&force[..r.len()]);
            endpoint+=force[r.len()];
        }
        // Pull the local physical endpoint force back through the SAME kinetic
        // coordinate map as storage. Its -c_i F_i term cannot be discarded.
        out[n]=self.0.weight*(endpoint-dot(&self.0.translation,&out[..n]));
        out.iter().all(|v|v.is_finite())
    }
}

impl LongitudinalPotential {
    fn viscous_action(&self,q:&[f64],v:&[f64],direction:Option<(&[f64],&[f64])>,
        rigidity:f64,out:&mut[f64]) -> bool
    {
        let count=self.count();let bend=self.omegas.len();let last=count-1;
        if q.len()!=count || v.len()!=count || out.len()!=count
            || q.iter().chain(v).any(|x|!x.is_finite()) {return false;}
        let slopes=self.slopes(q);let speed=self.slopes(v);
        let mut rates=vec![0.0;self.model.axial_modes().len()];
        let Ok(mut mean)=self.model.strain_rate(&slopes,&speed,&v[bend..last],&mut rates) else {return false;};
        mean+=self.chord[2]*v[last];
        let c=2.0*self.model.strain_energy_coefficient()*(rigidity/self.ea);
        if !c.is_finite() || c<=0.0 {return false;}
        let mut weights:Vec<_>=rates.iter().map(|rate|0.5*c*rate).collect();
        let mut gradient=vec![0.0;slopes.len()];
        if let Some((dq,dv))=direction {
            if dq.len()!=count || dv.len()!=count || dq.iter().chain(dv).any(|x|!x.is_finite()) {return false;}
            let ds=self.slopes(dq);let dspeed=self.slopes(dv);
            let mut drates=vec![0.0;rates.len()];
            let Ok(mut dmean)=self.model.strain_rate(&slopes,&dspeed,&dv[bend..last],&mut drates) else {return false;};
            dmean+=self.chord[2]*dv[last];
            // Differentiate B(q)v. The second contribution is the bilinear
            // geometric strain derivative; the axial/endpoint map is linear.
            let zero=vec![0.0;rates.len()];let mut cross=vec![0.0;rates.len()];
            let Ok(cross_mean)=self.model.strain_rate(&ds,&speed,&zero,&mut cross) else {return false;};
            dmean+=cross_mean;
            for (rate,extra) in drates.iter_mut().zip(cross) {*rate+=extra;}
            // The other product-rule term is (dB[dq])^T C B(q)v. This includes
            // the support's transverse chord slopes in exactly the same owner.
            if self.model.geometric_strain_gradient(&ds,c*mean,&weights,&mut gradient).is_err() {return false;}
            mean=dmean;
            for (weight,rate) in weights.iter_mut().zip(drates) {*weight=0.5*c*rate;}
        }
        let mut primary=vec![0.0;slopes.len()];
        if self.model.geometric_strain_gradient(&slopes,c*mean,&weights,&mut primary).is_err() {return false;}
        for (gradient,primary) in gradient.iter_mut().zip(primary) {*gradient+=primary;}
        let n=self.scales.len();
        for plane in 0..2 {for mode in 0..n {
            out[plane*n+mode]=-self.scales[mode]*gradient[plane*(n+1)+mode+1];
        }}
        for (i,(&scale,&weight)) in self.model.axial_slope_scales().iter().zip(&weights).enumerate() {
            out[bend+i]=-scale*weight;
        }
        out[last]=-(self.chord[0]*gradient[0]+self.chord[1]*gradient[n+1]+self.chord[2]*c*mean);
        out.iter().all(|v|v.is_finite())
    }
}
