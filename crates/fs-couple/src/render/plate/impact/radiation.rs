//! Passive multiport radiation memory inside the existing joint impact solve.
//!
//! Each signed row l realizes Z(s) = l^T l s/(s^2+2*zeta*omega*s+omega^2).
//! This is a model class, not evidence of agreement with any acoustic geometry.
//! The same class also feeds the piano's split runtime; here it is composed into
//! the original Gonzalez equation without an extra time integrator or split.
use super::{ImpactError, ImpactSystem, invalid};
use fs_math::c64::C64;

/// Shared fixed-dictionary PSD-residue fit, also used by the piano runtime.
pub mod fit;

/// One unit-mass acoustic oscillator and its signed mechanical velocity ports.
#[derive(Clone, Debug)]
pub struct Pole {
    pub omega: f64,
    pub zeta: f64,
    /// Generalized mechanical velocities -> acoustic modal force [1/s].
    pub coupling: Vec<f64>,
}
/// Complete reciprocal positive-real impedance in the declared velocity basis.
#[derive(Clone, Debug)]
pub struct Model {
    pub ports: usize,
    pub poles: Vec<Pole>,
}
impl Model {
    pub fn validate(&self) -> Result<(), String> {
        if !(1..=32).contains(&self.ports) || self.poles.is_empty() || self.poles.len()>1024
            || self.poles.iter().any(|p| !p.omega.is_finite() || p.omega<=0.
                || !p.zeta.is_finite() || !(0.0..1.0).contains(&p.zeta)
                || p.coupling.len()!=self.ports || p.coupling.iter().any(|x|!x.is_finite())) {
            return Err("invalid complete passive acoustic pole/port model".into());
        }
        Ok(())
    }
    /// Continuous impedance in exp(-i omega t), not the finite-step transfer.
    /// A zero-damping pole exactly at the query frequency refuses.
    pub fn impedance(&self, omega:f64) -> Result<Vec<C64>, String> {
        self.validate()?;
        if !omega.is_finite() || omega<=0. {return Err("invalid acoustic frequency".into());}
        let mut z=vec![C64::ZERO;self.ports*self.ports];
        for p in &self.poles {
            let den=C64::new(p.omega*p.omega-omega*omega,-2.*p.zeta*p.omega*omega);
            if den.abs()==0. {return Err("unresolved lossless acoustic pole".into());}
            let h=C64::new(0.,-omega)/den;
            for i in 0..self.ports {for j in 0..self.ports {
                z[i*self.ports+j]=z[i*self.ports+j]+h.scale(p.coupling[i]*p.coupling[j]);
            }}
        }
        if z.iter().any(|v|!v.re.is_finite()||!v.im.is_finite()) {return Err("acoustic impedance overflow".into());}
        Ok(z)
    }
}

/// These endpoint quantities are ALREADY included in the system's total ledger.
#[derive(Clone, Copy, Debug, Default)]
pub struct RadiationObservation {
    pub poles: usize,
    pub stored_energy_j: f64,
    /// Instantaneous acoustic damping power, not interval dissipation.
    pub dissipated_power_w: f64,
}
pub(super) struct Memory { pub(super) base_dim:usize, omega:Vec<f64>, rates:Vec<f64> }
impl Memory {
    pub(super) fn add_hessian(&self,d:&[f64],out:&mut[f64]) {
        for (i,w) in self.omega.iter().enumerate() {
            let q=self.base_dim+2*i;out[q]=w*w*d[q];out[q+1]=d[q+1];
        }
    }
    fn observe(&self,x:&[f64])->RadiationObservation {
        let mut result=RadiationObservation{poles:self.omega.len(),..Default::default()};
        for (i,(&w,&rate)) in self.omega.iter().zip(&self.rates).enumerate() {
            let q=self.base_dim+2*i;let v=x[q+1];
            result.stored_energy_j+=0.5*((w*x[q]).powi(2)+v*v);
            result.dissipated_power_w+=rate*v*v;
        }
        result
    }
}
impl ImpactSystem {
    /// Cold-attach a complete passive load to selected original mechanical modes.
    ///
    /// The existing modal-bank owner builds acoustic storage, damping and ports;
    /// `with_parallel_load` connects them to the SAME implicit time equation.
    /// All mechanical/contact/felt/material states and external force addresses
    /// retain their prefix. Acoustic q/p pairs follow the current state and start
    /// at rest. No acoustic coordinate becomes a striker, cavity, or radiation
    /// observer input. All loss/work/energy gates and full-state rollback apply.
    ///
    /// Distinct source_modes bind model columns without reordering or truncating
    /// them. At most 256 poles and 1024 total scalar states are admitted here;
    /// maximum_poles is the caller's additional explicit cold-work ceiling.
    /// Poles obey the existing Nyquist guard; coupling infinity-norm * dt <= .25
    /// is a separate resolution guard, not a temporal error bound. A completely
    /// zero load preserves the original trajectory and state dimensions exactly.
    ///
    /// # Errors
    /// Invalid/duplicate source indices, nonfinite or unresolved model, exceeded
    /// work limits, or repeat/late attachment. No failed fit is replaced by zeros.
    pub fn with_radiation_load(mut self,model:&Model,source_modes:&[usize],maximum_poles:usize)
        -> Result<Self,ImpactError>
    {
        model.validate().map_err(ImpactError::Owner)?;
        if self.sample!=0 || self.radiation.is_some() || maximum_poles>256
            || model.poles.len()>maximum_poles || source_modes.len()!=model.ports
            || source_modes.iter().enumerate().any(|(i,&k)|k>=self.modes||source_modes[..i].contains(&k)) {
            return Err(invalid("radiation load needs one cold complete mapping within its 256-pole ceiling"));
        }
        let dim=self.x.len().checked_add(2*model.poles.len()).ok_or_else(||invalid("radiation state overflow"))?;
        if dim>1024 {return Err(invalid("radiation load exceeds the joint 1024-state work envelope"));}
        let mut col_sum=vec![0.;model.ports];let mut row_max=0.0_f64;
        for p in &model.poles {
            if !p.omega.powi(2).is_finite() || !(2.*p.zeta*p.omega).is_finite()
                || p.omega*self.config.dt_s>=0.9*core::f64::consts::PI {
                return Err(invalid("radiation pole exceeds finite mechanical-rate bounds"));
            }
            row_max=row_max.max(p.coupling.iter().map(|x|x.abs()).sum());
            for (sum,l) in col_sum.iter_mut().zip(&p.coupling) {*sum+=l.abs();}
        }
        let norm=row_max.max(col_sum.into_iter().fold(0.0_f64,f64::max));
        if !norm.is_finite() || norm*self.config.dt_s>0.25 {
            return Err(invalid("radiation coupling exceeds mechanical-rate resolution"));
        }
        if norm==0. {return Ok(self);}
        let omega:Vec<_>=model.poles.iter().map(|p|p.omega).collect();
        let zeta:Vec<_>=model.poles.iter().map(|p|p.zeta).collect();
        let drives:Vec<Vec<_>>=(0..model.ports).map(|j|model.poles.iter().map(|p|p.coupling[j]).collect()).collect();
        let refs:Vec<_>=drives.iter().map(Vec::as_slice).collect();
        let load=fs_phs::modal_bank_ports(&omega,&zeta,&refs).map_err(ImpactError::PreparedSolve)?;
        self.system=self.system.with_parallel_load(load,source_modes).map_err(ImpactError::PreparedSolve)?;
        let base_dim=self.x.len();self.x.resize(dim,0.);
        self.radiation=Some(Memory {base_dim,omega,rates:model.poles.iter().map(|p|2.*p.zeta*p.omega).collect()});
        Ok(self)
    }
    pub fn radiation_observation(&self)->Option<RadiationObservation> {
        self.radiation.as_ref().map(|air|air.observe(&self.x))
    }
}

#[cfg(test)]
#[path = "radiation_tests.rs"]
mod tests;
