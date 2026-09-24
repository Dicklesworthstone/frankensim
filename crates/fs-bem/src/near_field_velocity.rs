//! Particle velocity and ideal first-order microphone observation of the SAME
//! exterior Green field. No new boundary equation or finite-difference stencil.
use super::*;

// Derivatives with respect to the observation point x, not the source y.
// For u=(x-y)/r: grad_x S=G' u and
// grad_x D=-(G''-G'/r)(n_y.u)u-(G'/r)n_y.
#[derive(Clone,Copy,Default)]
pub(super) struct Derivative {
    s:[C64;3], d:[C64;3], abs_s:[f64;3], abs_d:[f64;3],
    error_s:[f64;3], error_d:[f64;3],
}
impl Derivative {
    pub(super) fn accumulate(&mut self,k:f64,g:C64,delta:Point,r:f64,normal:Point,jac:f64) -> Result<(),HelmholtzError> {
        let u=delta.map(|v|v/r);let a=dot(normal,u);let inv=1./r;
        let first=g*C64::new(-inv,k);
        let radial=g*C64::new(3.*inv*inv-k*k,-3.*k*inv);
        let tangent=first.scale(inv);
        for c in 0..3 {
            let ds=first.scale(u[c]);let dd=-(radial.scale(a*u[c])+tangent.scale(normal[c]));
            if ![ds.re,ds.im,dd.re,dd.im].iter().all(|x|x.is_finite()) {
                return Err(bad("nonfinite near-field derivative kernel"));
            }
            self.s[c]=self.s[c]+ds.scale(jac);self.d[c]=self.d[c]+dd.scale(jac);
            self.abs_s[c]+=jac*ds.abs();self.abs_d[c]+=jac*dd.abs();
        }
        Ok(())
    }
    pub(super) fn add(&mut self,b:Self) {
        for c in 0..3 {
            self.s[c]=self.s[c]+b.s[c];self.d[c]=self.d[c]+b.d[c];
            self.abs_s[c]+=b.abs_s[c];self.abs_d[c]+=b.abs_d[c];
            self.error_s[c]+=b.error_s[c];self.error_d[c]+=b.error_d[c];
        }
    }
    pub(super) fn accept(&mut self,coarse:&Self,tolerance:f64)->bool {
        let mut admitted=true;
        for c in 0..3 {
            self.error_s[c]=(self.s[c]-coarse.s[c]).abs();
            self.error_d[c]=(self.d[c]-coarse.d[c]).abs();
            admitted &= [self.s[c].re,self.s[c].im,self.d[c].re,self.d[c].im,
                self.abs_s[c],self.abs_d[c],self.error_s[c],self.error_d[c]].iter().all(|x|x.is_finite())
                && self.error_s[c]<=tolerance*self.abs_s[c] && self.error_d[c]<=tolerance*self.abs_d[c];
        }
        admitted
    }
}
impl<'a> Geometry<'a> {
    /// Integrate pressure AND its analytic spatial gradient on one adaptive
    /// tree. The six derivative kernels have their own discrepancy admission;
    /// convergence of pressure alone does not establish velocity accuracy.
    /// The original geometry, separation, phase, depth and total work bounds
    /// remain. One kernel evaluation includes all selected components.
    /// # Errors
    /// Same as `prepare`, including failure of any derivative component.
    pub fn prepare_velocity(&self,k:f64,medium:Medium,options:Options)->Result<Prepared<'a>,HelmholtzError> {
        self.prepare_fields(k,medium,options,true)
    }
}

/// Co-located pressure and fluid particle velocity, not surface velocity or
/// microphone-diaphragm motion. All phasors use exp(-i omega t).
pub struct VectorObservation {
    /// Scalar pressure and its existing quadrature discrepancy estimate [Pa].
    pub scalar:Observation,
    /// Cartesian particle velocity [m/s], in the geometry's reference frame.
    pub particle_velocity_m_s:Vec<[C64;3]>,
    /// Componentwise propagated quadrature discrepancies [m/s]. Not bounds on
    /// BEM discretization, source-model error or cancellation-conditioned error.
    pub quadrature_error_estimate_m_s:Vec<[f64;3]>,
}
impl Prepared<'_> {
    /// Euler momentum: v=grad(p)/(i omega rho) for exp(-i omega t).
    /// Reuses the prepared kernels across all boundary solutions at this k.
    /// # Errors
    /// Scalar identity/finite checks, missing derivative rows, or overflow.
    pub fn evaluate_velocity(&self,solution:&RadiationSolution)->Result<VectorObservation,HelmholtzError> {
        if !self.with_velocity {return Err(bad("particle velocity requires prepare_velocity"));}
        let scalar=self.evaluate(solution)?;
        let omega_rho=self.k*self.medium.sound_speed*self.medium.density;
        let inverse=C64::new(0.,-1./omega_rho);
        let mut velocities=Vec::with_capacity(self.rows.len());let mut errors=Vec::with_capacity(self.rows.len());
        for row in &self.rows {
            let(mut re,mut im,mut cr,mut ci,mut error)=([0.;3],[0.;3],[0.;3],[0.;3],[0.;3]);
            for ((cell,&p),&v) in row.iter().zip(&solution.pressure).zip(&solution.velocity) {
                let q=v*C64::new(0.,omega_rho);
                for c in 0..3 {
                    let term=cell.gradient.d[c]*p-cell.gradient.s[c]*q;
                    compensated(&mut re[c],&mut cr[c],term.re);compensated(&mut im[c],&mut ci[c],term.im);
                    error[c]+=p.abs()*cell.gradient.error_d[c]+q.abs()*cell.gradient.error_s[c];
                }
            }
            let value:[C64;3]=std::array::from_fn(|c|C64::new(re[c],im[c])*inverse);
            let error=error.map(|e|e/omega_rho);
            if value.iter().any(|v|!v.re.is_finite()||!v.im.is_finite()) || error.iter().any(|e|!e.is_finite()) {
                return Err(bad("nonfinite near-field particle velocity or discrepancy"));
            }
            velocities.push(value);errors.push(error);
        }
        Ok(VectorObservation {scalar,particle_velocity_m_s:velocities,quadrature_error_estimate_m_s:errors})
    }
}

/// Ideal coincident pressure/velocity receiver, normalized to unit on-axis
/// plane-wave pressure sensitivity. The axis points FROM microphone TO its
/// front source, opposite the incoming wave's propagation direction.
/// alpha=1: omni; 1/2: cardioid; 0: figure eight. No diaphragm, electronic
/// response, low-frequency rolloff, noise, saturation or probe backreaction.
#[derive(Clone,Copy,Debug,PartialEq)]
pub struct FirstOrder { fraction:f64, front:[f64;3] }
impl Default for FirstOrder {
    fn default()->Self {Self {fraction:1.,front:[0.,0.,1.]}}
}
impl FirstOrder {
    /// Require a finite fraction in [0,1] and an explicit unit front axis.
    /// # Errors
    /// Invalid controls; axes are never normalized or guessed silently.
    pub fn new(pressure_fraction:f64,front_axis:[f64;3])->Result<Self,HelmholtzError> {
        if !pressure_fraction.is_finite() || !(0.0..=1.0).contains(&pressure_fraction)
            || front_axis.iter().any(|x|!x.is_finite()) || (norm(front_axis)-1.).abs()>1e-12 {
            return Err(bad("first-order receiver needs alpha in [0,1] and a unit front axis"));
        }
        Ok(Self {fraction:pressure_fraction,front:front_axis})
    }
    /// Pressure fraction alpha.
    pub fn pressure_fraction(self)->f64 {self.fraction}
    /// Unit axis pointing toward the front source.
    pub fn front_axis(self)->[f64;3] {self.front}
    /// Pa-equivalent output alpha*p-(1-alpha)*rho*c*(front.v).
    /// Near-field reactive velocity produces proximity response; no assumed
    /// source distance or plane-wave approximation is used in this projection.
    /// # Errors
    /// Nonfinite fields, invalid medium, or overflowing output.
    pub fn observe(self,pressure:C64,velocity:[C64;3],medium:Medium)->Result<C64,HelmholtzError> {
        if !medium.density.is_finite() || medium.density<=0. || !medium.sound_speed.is_finite() || medium.sound_speed<=0.
            || !pressure.re.is_finite() || !pressure.im.is_finite()
            || velocity.iter().any(|v|!v.re.is_finite()||!v.im.is_finite()) {
            return Err(bad("invalid first-order receiver field or medium"));
        }
        if self.fraction==1. {return Ok(pressure);}
        let axial=(0..3).fold(C64::ZERO,|a,c|a+velocity[c].scale(self.front[c]));
        let out=pressure.scale(self.fraction)-axial.scale((1.-self.fraction)*medium.density*medium.sound_speed);
        if !out.re.is_finite() || !out.im.is_finite() {return Err(bad("first-order receiver overflow"));}
        Ok(out)
    }
}
