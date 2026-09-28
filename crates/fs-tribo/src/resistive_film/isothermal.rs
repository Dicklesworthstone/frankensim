//! Compressible isothermal storage and mass transport on the existing film graph.
//!
//! Unlike the pressure-eliminated image, a sealed pocket retains its gas mass.
//! Each state z is energy-scaled mass: m = sqrt(p0 V0) z / (R T).
//! Free energy relative to ambient is H = p0 V [u log u - u + 1], u=p/p0.
//! The q gradient gives reciprocal gauge-pressure force; the z gradient is
//! sqrt(p0 V0) log(p/p0). Positive logarithmic-mean mobility gives exactly the
//! isothermal Poiseuille mass flux G (p_i^2-p_j^2)/(2 R T) at that gradient.
//!
//! No time integrator lives here. The caller supplies its joint discrete-gradient
//! effort to the dissipative port, so mass, motion and the free-energy balance
//! advance together. Isothermal, ideal-gas, no-slip, laminar thin-gap assumptions
//! remain required; there is no fluid inertia, thermal evolution or noise source.
use super::{FilmError, ResistiveFilm, MAX_CELLS, MAX_PORTS, dot, fail};

/// Explicit, fixed thermal and pressure reservoir; no gas-name lookup.
#[derive(Clone, Copy, Debug)]
pub struct AmbientGas {
    pub pressure_pa: f64,
    pub temperature_k: f64,
    pub specific_gas_constant_j_kg_k: f64,
}
#[derive(Clone, Copy, Debug, Default)]
pub struct GasObservation {
    /// Relative isothermal free energy, already part of joint mechanical storage.
    pub free_energy_j: f64,
    pub mass_kg: f64,
    pub minimum_pressure_pa: f64,
    pub maximum_pressure_pa: f64,
}
#[derive(Clone, Debug)]
pub struct IsothermalFilm {
    geometry: ResistiveFilm,
    gas: AmbientGas,
    initial_gap: Vec<f64>,
    scales: Vec<f64>,
}
struct State {
    gaps: [f64; MAX_CELLS],
    pressure: [f64; MAX_CELLS],
    observation: GasObservation,
}
impl IsothermalFilm {
    /// Cold initialisation at ambient pressure on the actual supplied geometry.
    /// No pressure solve or drainage requirement: closed channels can trap mass.
    /// The geometry's gauge-pressure limit still applies, but is not restricted
    /// to the small-pressure approximation required by the incompressible image.
    pub fn at_ambient(geometry: ResistiveFilm, q: &[f64], gas: AmbientGas) -> Result<Self, FilmError> {
        geometry.validate_configuration(q)?;
        let rt=gas.temperature_k*gas.specific_gas_constant_j_kg_k;
        if !gas.pressure_pa.is_finite() || gas.pressure_pa<=0.0 || !gas.pressure_pa.powi(2).is_finite()
            || !gas.temperature_k.is_finite() || gas.temperature_k<=0.0
            || !gas.specific_gas_constant_j_kg_k.is_finite() || gas.specific_gas_constant_j_kg_k<=0.0
            || !rt.is_finite() || rt<=0.0 {
            return fail("isothermal film requires positive finite absolute pressure, temperature and gas constant");
        }
        let mut initial_gap=Vec::with_capacity(geometry.cells.len());let mut scales=Vec::new();
        for cell in &geometry.cells {
            let h=cell.gap.reference_m-dot(&cell.gap.closure,q);
            let energy=gas.pressure_pa*cell.area_m2*h;
            if !energy.is_finite() || energy<=0.0 || !(energy/rt).is_finite() || energy/rt<=0.0 {
                return fail("isothermal film reference mass or energy is unrepresentable");
            }
            initial_gap.push(h);scales.push(energy.sqrt());
        }
        Ok(Self{geometry,gas,initial_gap,scales})
    }
    pub fn port_count(&self)->usize {self.geometry.port_count()}
    pub fn cell_count(&self)->usize {self.geometry.cell_count()}
    /// Positive scaled masses at the declared initial ambient state.
    pub fn initial_state(&self)->&[f64] {&self.scales}
    pub fn observe(&self,q:&[f64],z:&[f64])->Result<GasObservation,FilmError> {
        Ok(self.state(q,z)?.observation)
    }
    fn state(&self,q:&[f64],z:&[f64])->Result<State,FilmError> {
        self.geometry.validate_configuration(q)?;
        if z.len()!=self.cell_count() || z.iter().any(|v|!v.is_finite()||*v<=0.0) {
            return fail("isothermal film requires positive finite cell masses");
        }
        let mut s=State{gaps:[0.0;MAX_CELLS],pressure:[0.0;MAX_CELLS],
            observation:GasObservation{minimum_pressure_pa:f64::INFINITY,..Default::default()}};
        let p0=self.gas.pressure_pa;let rt=self.gas.temperature_k*self.gas.specific_gas_constant_j_kg_k;
        for (i,c) in self.geometry.cells.iter().enumerate() {
            let h=c.gap.reference_m-dot(&c.gap.closure,q);s.gaps[i]=h;
            let u=(z[i]/self.scales[i])*(self.initial_gap[i]/h);let p=p0*u;
            if !p.is_finite() || p<=0.0 || !p.powi(2).is_finite()
                || (p-p0).abs()>self.geometry.limits.maximum_pressure_pa {
                return fail("isothermal film absolute pressure or gauge domain exceeded");
            }
            s.pressure[i]=p;
            s.observation.free_energy_j+=p0*c.area_m2*h*relative_energy(u);
            s.observation.mass_kg+=(self.scales[i]*z[i])/rt;
            s.observation.minimum_pressure_pa=s.observation.minimum_pressure_pa.min(p);
            s.observation.maximum_pressure_pa=s.observation.maximum_pressure_pa.max(p);
        }
        if !s.observation.free_energy_j.is_finite() || s.observation.free_energy_j<0.0
            || !s.observation.mass_kg.is_finite() || s.observation.mass_kg<=0.0 {
            return fail("isothermal film storage overflow");
        }
        Ok(s)
    }
    /// Analytic storage derivatives. Mechanical reaction is minus gq.
    pub fn gradient_into(&self,q:&[f64],z:&[f64],gq:&mut[f64],gz:&mut[f64])->Result<(),FilmError> {
        if gq.len()!=self.port_count() || gz.len()!=self.cell_count() {return fail("film gradient shape");}
        let s=self.state(q,z)?;let mut a=[0.0;MAX_PORTS];let mut b=[0.0;MAX_CELLS];
        for (i,c) in self.geometry.cells.iter().enumerate() {
            let force=c.area_m2*(s.pressure[i]-self.gas.pressure_pa);
            for (j,w) in c.gap.closure.iter().enumerate(){a[j]+=w*force;}
            b[i]=self.scales[i]*log_ratio(s.pressure[i],self.gas.pressure_pa);
        }
        if a[..gq.len()].iter().chain(&b[..gz.len()]).any(|x|!x.is_finite()) {return fail("film gradient overflow");}
        gq.copy_from_slice(&a[..self.port_count()]);gz.copy_from_slice(&b[..self.cell_count()]);Ok(())
    }
    pub fn hessian_into(&self,q:&[f64],z:&[f64],dq:&[f64],dz:&[f64],
        hq:&mut[f64],hz:&mut[f64])->Result<(),FilmError> {
        if dq.len()!=self.port_count() || dz.len()!=self.cell_count() || hq.len()!=dq.len() || hz.len()!=dz.len()
            || dq.iter().chain(dz).any(|x|!x.is_finite()) {return fail("film Hessian direction/shape");}
        let s=self.state(q,z)?;let mut a=[0.0;MAX_PORTS];let mut b=[0.0;MAX_CELLS];
        for (i,c) in self.geometry.cells.iter().enumerate() {
            let relative=dz[i]/z[i]+dot(&c.gap.closure,dq)/s.gaps[i];
            let dp=s.pressure[i]*relative;
            for (j,w) in c.gap.closure.iter().enumerate(){a[j]+=c.area_m2*w*dp;}
            b[i]=self.scales[i]*relative;
        }
        if a[..hq.len()].iter().chain(&b[..hz.len()]).any(|x|!x.is_finite()) {return fail("film Hessian overflow");}
        hq.copy_from_slice(&a[..self.port_count()]);hz.copy_from_slice(&b[..self.cell_count()]);Ok(())
    }
    /// Positive resisting mass-state flow at a supplied discrete-gradient effort.
    /// Closed passages have EXACTLY zero mobility. Interior exchange conserves
    /// mass even for arbitrary discrete efforts, not only endpoint gradients.
    pub fn flow_into(&self,q:&[f64],z:&[f64],effort:&[f64],out:&mut[f64])->Result<f64,FilmError> {
        if effort.len()!=self.cell_count() || out.len()!=effort.len() || effort.iter().any(|e|!e.is_finite()) {
            return fail("film mass-flow effort/shape");
        }
        let s=self.state(q,z)?;let mut candidate=[0.0;MAX_CELLS];let mut power=0.0;
        for (k,c) in self.geometry.channels.iter().enumerate() {
            let h=c.gap.reference_m-dot(&c.gap.closure,q);if h<=0.0 {continue;}
            let pi=s.pressure[c.from];let pj=c.to.map_or(self.gas.pressure_pa,|j|s.pressure[j]);
            let mean=pressure_mean(pi,pj,0.0,0.0).0;
            let mobility=self.geometry.factors[k]*h*h*h*mean;
            let delta=effort[c.from]/self.scales[c.from]-c.to.map_or(0.0,|j|effort[j]/self.scales[j]);
            let flux=mobility*delta;candidate[c.from]+=flux/self.scales[c.from];
            if let Some(j)=c.to {candidate[j]-=flux/self.scales[j];}
            power+=mobility*delta*delta;
        }
        if !power.is_finite() || power<0.0 || candidate[..out.len()].iter().any(|x|!x.is_finite()) {
            return fail("film mass-flow overflow");
        }
        out.copy_from_slice(&candidate[..self.cell_count()]);Ok(power)
    }
    /// Exact directional derivative of both gap conductance and pressure mobility.
    #[allow(clippy::too_many_arguments)]
    pub fn flow_tangent_into(&self,q:&[f64],z:&[f64],e:&[f64],dq:&[f64],dz:&[f64],de:&[f64],
        out:&mut[f64])->Result<(),FilmError> {
        let n=self.cell_count();
        if e.len()!=n || dz.len()!=n || de.len()!=n || out.len()!=n || dq.len()!=self.port_count()
            || e.iter().chain(dq).chain(dz).chain(de).any(|x|!x.is_finite()) {return fail("film flow tangent shape");}
        let s=self.state(q,z)?;let mut dp=[0.0;MAX_CELLS];let mut result=[0.0;MAX_CELLS];
        for (i,c) in self.geometry.cells.iter().enumerate() {
            dp[i]=s.pressure[i]*(dz[i]/z[i]+dot(&c.gap.closure,dq)/s.gaps[i]);
        }
        for (k,c) in self.geometry.channels.iter().enumerate() {
            let h=c.gap.reference_m-dot(&c.gap.closure,q);if h<=0.0 {continue;}
            let dh=-dot(&c.gap.closure,dq);
            let pi=s.pressure[c.from];let pj=c.to.map_or(self.gas.pressure_pa,|j|s.pressure[j]);
            let (mean,dm)=pressure_mean(pi,pj,dp[c.from],c.to.map_or(0.0,|j|dp[j]));
            let g=self.geometry.factors[k]*h*h*h;
            let dg=3.0*self.geometry.factors[k]*h*h*dh;
            let delta=e[c.from]/self.scales[c.from]-c.to.map_or(0.0,|j|e[j]/self.scales[j]);
            let dd=de[c.from]/self.scales[c.from]-c.to.map_or(0.0,|j|de[j]/self.scales[j]);
            let df=(dg*mean+g*dm)*delta+g*mean*dd;
            result[c.from]+=df/self.scales[c.from];if let Some(j)=c.to {result[j]-=df/self.scales[j];}
        }
        if result[..n].iter().any(|x|!x.is_finite()) {return fail("film flow tangent overflow");}
        out.copy_from_slice(&result[..n]);Ok(())
    }
}
fn log_ratio(a:f64,b:f64)->f64 {
    let difference=a-b;
    if difference.abs()<0.5*b {(difference/b).ln_1p()}else{a.ln()-b.ln()}
}
fn relative_energy(u:f64)->f64 {
    let x=u-1.0;
    if x.abs()<1e-4 {x*x*(0.5+x*(-1.0/6.0+x*(1.0/12.0+x*(-1.0/20.0+x/30.0))))}
    else {u*u.ln()-u+1.0}
}
// L(p_i^2,p_j^2), and its directional derivative. Series remove the removable
// singularity at equal pressure without perturbing the pressure or mobility.
fn pressure_mean(a:f64,b:f64,da:f64,db:f64)->(f64,f64) {
    let d=log_ratio(a,b);let d2=d*d;
    let (s,ds)=if d.abs()<1e-3 {
        (1.0+d2*(1.0/6.0+d2*(1.0/120.0+d2/5040.0)),
         d*(1.0/3.0+d2*(1.0/30.0+d2/840.0)))
    } else {(d.sinh()/d,(d*d.cosh()-d.sinh())/d2)};
    (a*b*s,(da*b+a*db)*s+a*b*ds*(da/a-db/b))
}
#[cfg(test)]
mod tests;
