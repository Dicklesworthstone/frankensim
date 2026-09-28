//! Compressible squeeze-film storage inside the original percussion time owner.
//!
//! fs-tribo supplies geometry, ideal-gas free energy and mass mobility. fs-phs
//! advances those masses WITH mechanical motion, contacts and felt history.
//! No pressure reset at seal/reopening, separate gas step, or extra damping law.
pub use fs_tribo::resistive_film::isothermal::{AmbientGas, GasObservation};
use fs_tribo::resistive_film::{ResistiveFilm, MAX_CELLS};
use fs_tribo::resistive_film::isothermal::IsothermalFilm;
use fs_phs::Storage;
use std::rc::Rc;
use super::{ImpactError, ImpactSystem, MAX_IMPACT_MODES, invalid};

pub(super) struct Memory {
    pub(super) base_dim: usize,
    modes: usize,
    law: IsothermalFilm,
}
impl Memory {
    fn end(&self) -> usize { self.base_dim + self.law.cell_count() }
    fn coordinates<'a>(&self, x: &'a [f64], q: &mut [f64; MAX_IMPACT_MODES]) -> Option<&'a [f64]> {
        if x.len() < self.end() { return None; }
        for i in 0..self.modes { q[i] = x[2*i]; }
        Some(&x[self.base_dim..self.end()])
    }
    pub(super) fn observe(&self, x: &[f64]) -> Result<GasObservation, ImpactError> {
        let mut q=[0.0;MAX_IMPACT_MODES];
        let z=self.coordinates(x,&mut q).ok_or_else(||invalid("short gas-film state"))?;
        self.law.observe(&q[..self.modes],z).map_err(|e|invalid(e.0))
    }
    fn add_gradient(&self, x: &[f64], out: &mut [f64]) -> bool {
        if out.len()!=x.len() {return false;}
        let mut q=[0.0;MAX_IMPACT_MODES];let mut gq=[0.0;MAX_IMPACT_MODES];let mut gz=[0.0;MAX_CELLS];
        let Some(z)=self.coordinates(x,&mut q) else {return false;};
        if self.law.gradient_into(&q[..self.modes],z,&mut gq[..self.modes],
            &mut gz[..self.law.cell_count()]).is_err() {return false;}
        for i in 0..self.modes {out[2*i]+=gq[i];}
        out[self.base_dim..self.end()].copy_from_slice(&gz[..self.law.cell_count()]);
        out.iter().all(|v|v.is_finite())
    }
    pub(super) fn add_hessian(&self, x:&[f64], d:&[f64], out:&mut[f64]) -> bool {
        if d.len()!=x.len() || out.len()!=x.len() {return false;}
        let mut q=[0.0;MAX_IMPACT_MODES];let mut dq=[0.0;MAX_IMPACT_MODES];
        let mut hq=[0.0;MAX_IMPACT_MODES];let mut hz=[0.0;MAX_CELLS];
        let Some(z)=self.coordinates(x,&mut q) else {return false;};
        for i in 0..self.modes {dq[i]=d[2*i];}
        if self.law.hessian_into(&q[..self.modes],z,&dq[..self.modes],&d[self.base_dim..self.end()],
            &mut hq[..self.modes],&mut hz[..self.law.cell_count()]).is_err() {return false;}
        for i in 0..self.modes {out[2*i]+=hq[i];}
        for (o,h) in out[self.base_dim..self.end()].iter_mut().zip(hz) {*o+=h;}
        out.iter().all(|v|v.is_finite())
    }
    pub(super) fn add_flow(&self, x:&[f64], effort:&[f64], out:&mut[f64]) -> bool {
        if effort.len()!=x.len() || out.len()!=x.len() {return false;}
        let mut q=[0.0;MAX_IMPACT_MODES];let mut flow=[0.0;MAX_CELLS];
        let Some(z)=self.coordinates(x,&mut q) else {return false;};
        if self.law.flow_into(&q[..self.modes],z,&effort[self.base_dim..self.end()],
            &mut flow[..self.law.cell_count()]).is_err() {return false;}
        for (o,f) in out[self.base_dim..self.end()].iter_mut().zip(flow) {*o+=f;}
        out.iter().all(|v|v.is_finite())
    }
    pub(super) fn add_flow_tangent(&self, x:&[f64], e:&[f64], dx:&[f64], de:&[f64], out:&mut[f64]) -> bool {
        if e.len()!=x.len() || dx.len()!=x.len() || de.len()!=x.len() || out.len()!=x.len() {return false;}
        let mut q=[0.0;MAX_IMPACT_MODES];let mut dq=[0.0;MAX_IMPACT_MODES];let mut flow=[0.0;MAX_CELLS];
        let Some(z)=self.coordinates(x,&mut q) else {return false;};
        for i in 0..self.modes {dq[i]=dx[2*i];}
        let range=self.base_dim..self.end();
        if self.law.flow_tangent_into(&q[..self.modes],z,&e[range.clone()],&dq[..self.modes],
            &dx[range.clone()],&de[range.clone()],&mut flow[..self.law.cell_count()]).is_err() {return false;}
        for (o,f) in out[range].iter_mut().zip(flow) {*o+=f;}
        out.iter().all(|v|v.is_finite())
    }
}
struct WithGas {
    base: Box<dyn Storage>,
    gas: Rc<Memory>,
}
impl Storage for WithGas {
    fn hamiltonian(&self,x:&[f64])->f64 {
        let Ok(gas)=self.gas.observe(x) else {return f64::NAN;};
        self.base.hamiltonian(&x[..self.gas.base_dim])+gas.free_energy_j
    }
    fn gradient(&self,x:&[f64],out:&mut[f64]) {
        if x.len()<self.gas.end() || out.len()!=x.len() {out.fill(f64::NAN);return;}
        out.fill(0.0);
        self.base.gradient(&x[..self.gas.base_dim],&mut out[..self.gas.base_dim]);
        if !self.gas.add_gradient(x,out) {out.fill(f64::NAN);}
    }
}
impl ImpactSystem {
    /// Cold-attach one compressible graph, initially at explicit ambient pressure.
    ///
    /// Adds one energy-scaled positive mass per fluid cell. All existing state,
    /// felt histories and mechanical force/source indices retain their prefix.
    /// Pressure reaction is a storage gradient, while mass exchange contributes
    /// nonnegative dissipation at the SAME discrete-gradient effort as contact.
    /// Ambient pressure/temperature define a fixed reservoir: the total ledger
    /// includes relative gas free energy, not isolated thermal internal energy.
    ///
    /// Gas may precede or follow radiation/material memory. At most 64 gas cells
    /// and 1024 total scalar states are admitted; no cell is dropped to fit.
    /// Simultaneous resistance-only and compressible images refuse rather than
    /// double-counting the same interface. No attachment after time has advanced.
    pub fn with_compressible_squeeze_film(mut self, geometry:ResistiveFilm, gas:AmbientGas)
        -> Result<Self,ImpactError>
    {
        if self.sample!=0 || self.gas_film.is_some() || !self.squeeze_films.is_empty()
            || geometry.port_count()!=self.modes {
            return Err(invalid("compressible film needs one cold exclusive full-basis gas graph"));
        }
        let base_dim=self.x.len();
        let dim=base_dim.checked_add(geometry.cell_count()).ok_or_else(||invalid("gas state overflow"))?;
        if dim>1024 {return Err(invalid("gas film exceeds joint 1024-state work ceiling"));}
        let mut q=[0.0;MAX_IMPACT_MODES];for i in 0..self.modes {q[i]=self.x[2*i];}
        let law=IsothermalFilm::at_ambient(geometry,&q[..self.modes],gas).map_err(|e|invalid(e.0))?;
        let memory=Rc::new(Memory{base_dim,modes:self.modes,law});
        let storage_memory=Rc::clone(&memory);
        self.system=self.system.with_appended_storage(memory.law.cell_count(),move |base,_|
            Box::new(WithGas{base,gas:storage_memory})).map_err(ImpactError::PreparedSolve)?;
        self.x.extend_from_slice(memory.law.initial_state());
        let energy=self.stored_energy_j();
        if !energy.is_finite() || energy<0.0 || energy>self.config.maximum_energy_j {
            return Err(invalid("initial gas storage exceeds mechanical energy ceiling"));
        }
        self.gas_film=Some(memory);Ok(self)
    }
    /// Absolute pressure, physical mass and relative free energy at accepted time.
    /// Observation neither solves pressure nor advances a separate gas history.
    pub fn gas_film_observation(&self)->Result<Option<GasObservation>,ImpactError> {
        self.gas_film.as_ref().map(|gas|gas.observe(&self.x)).transpose()
    }
}

#[cfg(test)]
mod tests;
