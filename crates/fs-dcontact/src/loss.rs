//! The existing Hunt--Crossley law as a passive implicit state-flow port.
//! No contact potential, quadrature, coefficient, or constitutive history changes.
use super::{ContactStorage, Obstacle, det};

impl Obstacle {
    // Shared by the legacy modal-force API and the new allocation-free port.
    pub(super) fn loss_force_at(&self, i:usize, penetration:f64, rate:f64)->f64 {
        let elastic=self.weights[i]*self.stiffness*det::pow(penetration,self.alpha);
        // Form the dimensionless factor before multiplying by elastic force.
        // Preserve NaN propagation and the existing non-attractive unloading cap.
        let factor=self.internal_loss*rate;
        elastic*if factor < -1.0 {-1.0}else{factor}
    }
}

impl ContactStorage {
    fn loss_dimensions(&self,x:&[f64],e:&[f64],out:&[f64])->bool {
        let Some(prefix)=self.n_modes.checked_mul(2) else {return false;};
        if x.len()<prefix || e.len()!=x.len() || out.len()!=x.len()
            || x.iter().chain(e).any(|v|!v.is_finite()) {return false;}
        // Also refuse malformed unsafe/raw obstacles rather than index them.
        self.obstacles.iter().all(|ob| ob.n_points.checked_mul(self.n_modes)==Some(ob.collocation.len())
            && ob.gaps.len()==ob.n_points && ob.weights.len()==ob.n_points
            && ob.internal_loss.is_finite() && ob.internal_loss>=0.0
            && ob.stiffness.is_finite() && ob.stiffness>=0.0
            && ob.alpha.is_finite() && ob.alpha>=1.0)
    }

    /// Fill the positive resisting flow D for `xdot=(J-R)e-D+Gu`.
    /// The caller supplies the SAME effort e used by its implicit step equation.
    /// Only mechanical momentum entries are loaded. Additional internal states
    /// (for example Kelvin memory) are accepted and receive exact zeros.
    /// Each point contributes `weight*K*p_+^alpha*max(chi*rate,-1)*b`;
    /// `rate=b dot e_momentum`, so `e dot D` is nonnegative pointwise.
    /// This is the existing loss increment, NOT the conservative reaction.
    ///
    /// No allocation or physical state/history mutation. False means the scratch
    /// output must be discarded; shape failure leaves it untouched. A time owner
    /// must evaluate this inside its implicit equation, not lag it as an input.
    #[must_use]
    pub fn dissipative_flow_into(&self,x:&[f64],e:&[f64],out:&mut[f64])->bool {
        if !self.loss_dimensions(x,e,out) {return false;}
        out.fill(0.0);
        for ob in &self.obstacles {
            if ob.internal_loss==0.0 {continue;}
            for i in 0..ob.n_points {
                let p=ob.penetration_at(self.n_modes,x,i);
                if !p.is_finite() {return false;}
                if p<=0.0 {continue;}
                let b=&ob.collocation[i*self.n_modes..(i+1)*self.n_modes];
                let rate=b.iter().enumerate().map(|(k,b)| b*e[2*k+1]).sum::<f64>();
                let f=ob.loss_force_at(i,p,rate);
                if !rate.is_finite() || !f.is_finite() {return false;}
                for (k,b) in b.iter().enumerate() {out[2*k+1]+=b*f;}
            }
        }
        out.iter().all(|v|v.is_finite())
    }

    /// Analytic action `D_x*dx + D_e*de` for the above port. State and effort
    /// directions are independent: e may contain a Gonzalez correction.
    /// At separation the inactive derivative is zero. At the unloading-cap
    /// corner the unsaturated one-sided derivative is selected. Away from these
    /// corners this is the exact derivative, including both penetration and rate.
    /// Same scratch/refusal and allocation contract as `dissipative_flow_into`.
    #[must_use]
    pub fn dissipative_flow_tangent_into(&self,x:&[f64],e:&[f64],dx:&[f64],de:&[f64],out:&mut[f64])->bool {
        if !self.loss_dimensions(x,e,out) || dx.len()!=x.len() || de.len()!=x.len()
            || dx.iter().chain(de).any(|v|!v.is_finite()) {return false;}
        out.fill(0.0);
        for ob in &self.obstacles {
            if ob.internal_loss==0.0 {continue;}
            for i in 0..ob.n_points {
                let p=ob.penetration_at(self.n_modes,x,i);
                if !p.is_finite() {return false;}
                if p<=0.0 {continue;}
                let b=&ob.collocation[i*self.n_modes..(i+1)*self.n_modes];
                let mut rate=0.0;let mut dp=0.0;let mut dv=0.0;
                for (k,b) in b.iter().enumerate() {rate+=b*e[2*k+1];dp+=b*dx[2*k];dv+=b*de[2*k+1];}
                let factor=ob.internal_loss*rate;
                let elastic=ob.weights[i]*ob.stiffness*det::pow(p,ob.alpha);
                let slope=ob.weights[i]*ob.stiffness*ob.alpha*det::pow(p,ob.alpha-1.0);
                let df=if factor < -1.0 {-slope*dp}
                    else {slope*dp*factor+elastic*(ob.internal_loss*dv)};
                if ![rate,dp,dv,df].iter().all(|v|v.is_finite()) {return false;}
                for (k,b) in b.iter().enumerate() {out[2*k+1]+=b*df;}
            }
        }
        out.iter().all(|v|v.is_finite())
    }
}
