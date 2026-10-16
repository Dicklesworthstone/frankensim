//! Exterior pressure reacting on compact necks in the SAME impact equation.
//!
//! A supplied neck can either contain an effective exterior end correction,
//! or leave exterior inertia/loss to an impedance load, never both. This module
//! records that physical choice; it does not guess or subtract an end correction.
//! The existing cavity and radiation owners retain all state and integration.
use super::{CavityCoupling,CavityNeck,ImpactError,invalid};
use super::super::super::{ImpactSystem,radiation::Model};
use fs_exec::CancelGate;

impl CavityCoupling {
    /// Append necks whose supplied length and resistance EXCLUDE exterior
    /// radiation. Length still includes the duct and any declared interior-end
    /// correction; resistance includes only retained internal losses. Inputs
    /// are caller declarations, not deductions from radius or a measured card.
    ///
    /// The ordinary `with_necks` keeps its original effective-length/reservoir
    /// semantics and arithmetic. Here the identical compact-neck law supplies
    /// the interior participant of a coupled load. Until a load is attached,
    /// it is the explicitly unloaded, zero-exterior-pressure comparison.
    /// All existing compactness, finite-value, cancellation and work gates apply.
    pub fn with_necks_for_exterior_load(self,necks:Vec<CavityNeck>,gate:&CancelGate)
        ->Result<Self,ImpactError> {
        let first=self.necks.len();
        let mut result=self.with_necks(necks,gate)?;
        for neck in &mut result.necks[first..] {neck.exterior_load=true;}
        Ok(result)
    }

    /// Whether the declared neck impedance leaves exterior radiation to a load.
    /// This is NOT evidence that a load has been attached or successfully fitted.
    pub fn neck_accepts_exterior_load(&self,index:usize)->Result<bool,ImpactError> {
        self.necks.get(index).map(|n|n.exterior_load)
            .ok_or_else(||invalid("unknown cavity neck termination"))
    }

    /// Admit the complete exterior source map before attaching pressure feedback.
    /// Every opening must use the internal-only impedance convention and appear
    /// exactly once. Interior pressure-basis inertias are not surface sources.
    /// Structural source selection remains with the actual boundary geometry.
    ///
    /// The columns are unit generalized velocities, including neck momentum:
    /// Q = volume_weight * p. A uniform mouth therefore needs the area-normalized
    /// surface row volume_weight / mouth_area, not unit physical volume flow.
    pub fn admit_exterior_sources(&self,source_modes:&[usize])->Result<(),ImpactError> {
        if source_modes.is_empty() || source_modes.iter().enumerate().any(|(i,&mode)|
            source_modes[..i].contains(&mode) || (mode>=self.structural
                && !self.necks.iter().any(|n|n.coordinate==mode))) {
            return Err(invalid("exterior sources must be distinct structural coordinates or actual neck flows"));
        }
        if self.necks.iter().any(|n|!n.exterior_load || !source_modes.contains(&n.coordinate)) {
            return Err(invalid("every radiating neck needs internal-only impedance and its complete exterior source"));
        }
        Ok(())
    }

    /// Attach one full multiport load to the existing solid/interior-air/neck
    /// system. The original force ports, coordinates, histories and clock remain
    /// in place. Signed head-mouth and mouth-mouth reactions share the SAME
    /// Gonzalez solve; no prescribed-flow feedback lag or additional integrator.
    ///
    /// This consumes a system built from this cavity's coordinate layout. Every
    /// neck must have been declared with `with_necks_for_exterior_load`; an old
    /// effective-length neck refuses rather than double-counting exterior mass
    /// or loss. Fitting to a particular BEM boundary remains the caller's job.
    pub fn attach_exterior_radiation(&self,system:ImpactSystem,model:&Model,
        source_modes:&[usize],maximum_poles:usize)->Result<ImpactSystem,ImpactError> {
        self.admit_exterior_sources(source_modes)?;
        if system.modes!=self.total {
            return Err(invalid("exterior cavity load lost its original mechanical coordinate layout"));
        }
        system.with_radiation_load(model,source_modes,maximum_poles)
    }

    /// Effective mouth-averaged exterior pressure [Pa] exerted by the attached
    /// impedance. Positive pressure opposes positive OUTWARD neck volume flow.
    /// This is the time-domain load pressure, not the pressure at a microphone.
    /// With no load it is zero, making the unloaded comparison explicit.
    ///
    /// F_neck = -volume_weight * p_exterior, hence F_neck * p_neck
    /// = -p_exterior * Q. Cross-coupled head radiation is retained in F_neck.
    /// It must not be added to the external-force or loss ledger a second time.
    pub fn neck_exterior_pressure(&self,system:&ImpactSystem,index:usize)->Result<f64,ImpactError> {
        if system.modes!=self.total {
            return Err(invalid("exterior neck observation lost its mechanical layout"));
        }
        let neck=self.necks.get(index).ok_or_else(||invalid("unknown exterior neck pressure port"))?;
        let pressure=-system.radiation_force(neck.coordinate)?/neck.volume_weight;
        if !pressure.is_finite() {return Err(invalid("exterior neck pressure overflow"));}
        Ok(pressure)
    }
}

#[cfg(test)]
#[path="exterior_tests.rs"]
mod tests;
