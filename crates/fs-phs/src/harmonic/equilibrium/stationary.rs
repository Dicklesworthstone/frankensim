//! Original stationary equations and matrix-free analytic actions for an outer
//! nonlinear solver. This module owns no Newton iteration or time integration.
//! Numerical coordinates z have x_i=state_scales_i*z_i; residual rows use fixed
//! balance_scales. Both are caller data, never adjusted to manufacture convergence.
use super::*;

/// Fixed physical input and coordinate/row scales for one stationary solve.
#[derive(Clone, Copy, Debug)]
pub struct StationaryCoordinates<'a> {
    pub input: &'a [f64],
    /// Positive physical magnitude of each original state/constraint multiplier.
    pub state_scales: &'a [f64],
    /// Positive magnitude of each original flow/constraint equation.
    pub balance_scales: &'a [f64],
}

/// Read-only original source, including every constraint and material state.
/// Callbacks must describe this exact source and remain pure for the complete
/// solve. A domain predicate admits configurations, not equilibrium residuals.
/// No claim that a stationary solution exists, is unique, or is stable is made.
pub struct StationaryProblem<'a> {
    n: usize,
    nd: usize,
    m: usize,
    j: &'a [f64],
    r: &'a [f64],
    g: &'a [f64],
    storage: &'a dyn Storage,
    coordinates: StationaryCoordinates<'a>,
    budget: LinearizationBudget,
    hessian: &'a dyn Fn(&[f64], &[f64], &mut [f64]) -> bool,
    loss: Option<NonlinearDissipation<'a>>,
    domain: &'a dyn Fn(&[f64]) -> bool,
}

/// Scratch and one exact-point cache; not an accepted physical trajectory.
/// Repeated Krylov actions at one point reuse its storage gradient and resisting
/// force. Cache comparison is bitwise, with no tolerance or quantized state.
pub struct StationaryWorkspace<'p,'s> {
    owner: &'p StationaryProblem<'s>,
    z: Vec<f64>,
    x: Vec<f64>,
    effort: Vec<f64>,
    residual: Vec<f64>,
    dx: Vec<f64>,
    de: Vec<f64>,
    loss: Vec<f64>,
    action: Vec<f64>,
    valid: bool,
}

impl<'a> StationaryProblem<'a> {
    pub fn from_ode(
        source: &'a PortHamiltonian, coordinates: StationaryCoordinates<'a>,
        budget: LinearizationBudget,
        hessian: &'a dyn Fn(&[f64], &[f64], &mut [f64]) -> bool,
        loss: Option<NonlinearDissipation<'a>>, domain: &'a dyn Fn(&[f64]) -> bool,
    ) -> Result<Self, LinearizationError> {
        Self::new(source.n, source.n, source.m, &source.j, &source.r, &source.g,
            source.storage.as_ref(), coordinates, budget, hessian, loss, domain)
    }

    pub fn from_descriptor(
        source: &'a DescriptorPortHamiltonian, coordinates: StationaryCoordinates<'a>,
        budget: LinearizationBudget,
        hessian: &'a dyn Fn(&[f64], &[f64], &mut [f64]) -> bool,
        loss: Option<NonlinearDissipation<'a>>, domain: &'a dyn Fn(&[f64]) -> bool,
    ) -> Result<Self, LinearizationError> {
        Self::new(source.n, source.n_diff, source.m, &source.j, &source.r, &source.g,
            source.storage.as_ref(), coordinates, budget, hessian, loss, domain)
    }

    #[allow(clippy::too_many_arguments)]
    fn new(n: usize, nd: usize, m: usize, j: &'a [f64], r: &'a [f64], g: &'a [f64],
        storage: &'a dyn Storage, coordinates: StationaryCoordinates<'a>,
        budget: LinearizationBudget,
        hessian: &'a dyn Fn(&[f64], &[f64], &mut [f64]) -> bool,
        loss: Option<NonlinearDissipation<'a>>, domain: &'a dyn Fn(&[f64]) -> bool,
    ) -> Result<Self, LinearizationError> {
        let h=budget.harmonic;
        if n==0 || nd>n || n>MAX_HARMONIC_DIMENSION || n>h.max_dimension
            || n>budget.max_hessian_actions || loss.is_some_and(|l|n>l.max_tangent_actions)
            || !h.relative_residual_tolerance.is_finite()
            || h.relative_residual_tolerance<=0.0 || h.relative_residual_tolerance>=1.0 {
            return Err(invalid("stationary source exceeds complete-state/action capacity"));
        }
        let (nn,mm)=(n as u128,m as u128);
        // Include one matrix-free workspace beside the final dense local image.
        let bytes=128*nn*nn+16*nn*mm+640*nn+128*mm;
        let work=nn*nn*(nn+mm+24);
        if bytes>h.max_dense_bytes as u128 || work>h.max_dense_work as u128 {
            return Err(invalid("stationary preparation exceeds physical payload/work capacity"));
        }
        if j.len()!=n*n || r.len()!=n*n || n.checked_mul(m)!=Some(g.len())
            || j.iter().chain(r).chain(g).any(|v|!v.is_finite())
            || coordinates.input.len()!=m || coordinates.input.iter().any(|v|!v.is_finite())
            || coordinates.state_scales.len()!=n || coordinates.balance_scales.len()!=n
            || coordinates.state_scales.iter().chain(coordinates.balance_scales)
                .any(|v|!v.is_finite() || *v<=0.0) {
            return Err(invalid("stationary coordinates require complete finite inputs and positive scales"));
        }
        Ok(Self{n,nd,m,j,r,g,storage,coordinates,budget,hessian,loss,domain})
    }

    pub fn dimension(&self)->usize {self.n}

    /// Fresh scratch borrows this immutable problem. Cross-wiring another source
    /// of the same dimension is refused before reading a cached point.
    pub fn workspace(&self)->StationaryWorkspace<'_,'a> {
        StationaryWorkspace {owner:self,z:vec![0.;self.n],x:vec![0.;self.n],effort:vec![0.;self.n],
            residual:vec![0.;self.n],dx:vec![0.;self.n],de:vec![0.;self.n],
            loss:vec![0.;self.n],action:vec![0.;self.n],valid:false}
    }

    pub fn scaled_state(&self,state:&[f64])->Result<Vec<f64>,LinearizationError> {
        if state.len()!=self.n || state.iter().any(|v|!v.is_finite()) {
            return Err(invalid("complete finite stationary seed required"));
        }
        let z:Vec<_>=state.iter().zip(self.coordinates.state_scales).map(|(x,s)|x/s).collect();
        if z.iter().zip(state).any(|(v,x)|!v.is_finite() || (*x!=0.0 && *v==0.0)) {
            return Err(invalid("stationary seed scaling is not representable"));
        }
        Ok(z)
    }

    fn map_state(&self,z:&[f64],out:&mut[f64])->Result<(),LinearizationError> {
        if z.len()!=self.n || out.len()!=self.n || z.iter().any(|v|!v.is_finite()) {
            return Err(invalid("complete finite numerical stationary state required"));
        }
        for ((o,v),s) in out.iter_mut().zip(z).zip(self.coordinates.state_scales) {
            *o=v*s;
            if !o.is_finite() || (*v!=0.0 && *o==0.0) {
                return Err(invalid("physical stationary state is not representable"));
            }
        }
        if !(self.domain)(out) {return Err(invalid("stationary trial outside original source domain"));}
        Ok(())
    }

    fn point(&self,z:&[f64],w:&mut StationaryWorkspace<'_,'_>,cancel:&mut impl FnMut()->bool)
        ->Result<(),LinearizationError> {
        poll(cancel)?;
        if !std::ptr::eq(self,w.owner) {return Err(invalid("stationary workspace belongs to another problem"));}
        if w.valid && z.len()==self.n && z.iter().zip(&w.z).all(|(a,b)|a.to_bits()==b.to_bits()) {
            return Ok(());
        }
        w.valid=false;
        self.map_state(z,&mut w.x)?;
        let energy=self.storage.hamiltonian(&w.x);
        w.effort.fill(f64::NAN);
        self.storage.gradient(&w.x,&mut w.effort);
        poll(cancel)?;
        if !energy.is_finite() || w.effort.iter().any(|v|!v.is_finite())
            || w.effort[self.nd..].iter().any(|v|*v!=0.0) {
            return Err(invalid("stationary storage effort is invalid or contains algebraic storage"));
        }
        w.effort[self.nd..].copy_from_slice(&w.x[self.nd..]);
        if let Some(law)=self.loss {
            // Reuse the time/equilibrium law's original conjugate-power admission.
            let (flow,_)=law.evaluate(&w.x,&w.effort,cancel)?;
            w.loss.copy_from_slice(&flow);
        } else {w.loss.fill(0.0);}
        for i in 0..self.n {
            poll(cancel)?;
            let (mut value,mut correction)=(0.0,0.0);
            for j in 0..self.n {add(&mut value,&mut correction,
                (self.j[i*self.n+j]-self.r[i*self.n+j])*w.effort[j]);}
            add(&mut value,&mut correction,-w.loss[i]);
            for p in 0..self.m {add(&mut value,&mut correction,self.g[i*self.m+p]*self.coordinates.input[p]);}
            w.residual[i]=value/self.coordinates.balance_scales[i];
            if !value.is_finite() || !w.residual[i].is_finite() || (value!=0.0 && w.residual[i]==0.0) {
                return Err(invalid("stationary balance scaling is not representable"));
            }
        }
        poll(cancel)?;
        w.z.copy_from_slice(z);w.valid=true;
        Ok(())
    }

    /// Complete original stationary balance, divided by the declared row scales.
    /// Caller output is untouched on failure; only private scratch may change.
    pub fn residual_into(&self,z:&[f64],w:&mut StationaryWorkspace<'_,'_>,out:&mut[f64],
        mut cancelled:impl FnMut()->bool)->Result<(),LinearizationError> {
        if out.len()!=self.n {return Err(invalid("stationary residual output shape"));}
        self.point(z,w,&mut cancelled)?;
        poll(&mut cancelled)?;out.copy_from_slice(&w.residual);Ok(())
    }

    /// Exact scaled Jacobian action: ((J-R) W dx - D_x dx - D_e W dx)/row_scale.
    /// One source Hessian action per call; no assembled Jacobian or inverse E/Q.
    pub fn jacobian_action_into(&self,z:&[f64],direction:&[f64],w:&mut StationaryWorkspace<'_,'_>,
        out:&mut[f64],mut cancelled:impl FnMut()->bool)->Result<(),LinearizationError> {
        if direction.len()!=self.n || out.len()!=self.n || direction.iter().any(|v|!v.is_finite()) {
            return Err(invalid("stationary Jacobian direction/output shape"));
        }
        self.point(z,w,&mut cancelled)?;
        for i in 0..self.n {w.dx[i]=direction[i]*self.coordinates.state_scales[i];}
        if w.dx.iter().zip(direction).any(|(v,d)|!v.is_finite() || (*d!=0.0 && *v==0.0)) {
            return Err(invalid("stationary direction overflow/underflow"));
        }
        w.de.fill(f64::NAN);
        if !(self.hessian)(&w.x,&w.dx,&mut w.de) || w.de.iter().any(|v|!v.is_finite())
            || w.de[self.nd..].iter().any(|v|*v!=0.0) {
            return Err(invalid("source refused the complete stationary Hessian action"));
        }
        poll(&mut cancelled)?;
        w.de[self.nd..].copy_from_slice(&w.dx[self.nd..]);
        if let Some(law)=self.loss {
            w.loss.fill(f64::NAN);
            if !(law.tangent)(&w.x,&w.effort,&w.dx,&w.de,&mut w.loss)
                || w.loss.iter().any(|v|!v.is_finite()) {
                return Err(invalid("source refused the complete stationary loss action"));
            }
        } else {w.loss.fill(0.0);}
        for i in 0..self.n {
            poll(&mut cancelled)?;
            let (mut value,mut correction)=(0.0,0.0);
            for j in 0..self.n {add(&mut value,&mut correction,
                (self.j[i*self.n+j]-self.r[i*self.n+j])*w.de[j]);}
            add(&mut value,&mut correction,-w.loss[i]);
            w.action[i]=value/self.coordinates.balance_scales[i];
            if !value.is_finite() || !w.action[i].is_finite() || (value!=0.0 && w.action[i]==0.0) {
                return Err(invalid("stationary Jacobian action overflow/underflow"));
            }
        }
        poll(&mut cancelled)?;out.copy_from_slice(&w.action);Ok(())
    }

    /// Independently re-evaluate the original equilibrium and construct the
    /// existing complete resonance image. A solver's convergence flag is not
    /// trusted. This also checks symmetry and all algebraic Hessian columns.
    pub fn linearize(&self,z:&[f64],relative_tolerance:f64,mut cancelled:impl FnMut()->bool)
        ->Result<EquilibriumLinearization,LinearizationError> {
        poll(&mut cancelled)?;
        let mut state=vec![0.;self.n];self.map_state(z,&mut state)?;
        EquilibriumLinearization::prepare(self.n,self.nd,self.m,self.j,self.r,self.g,self.storage,
            EquilibriumPoint{state:&state,input:self.coordinates.input,
                balance_scales:self.coordinates.balance_scales,relative_tolerance},
            self.budget,self.hessian,self.loss,&mut cancelled)
    }
}
