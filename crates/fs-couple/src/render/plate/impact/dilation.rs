//! Admit physical uniform edge motion to the existing joint impact equation.
//! fs-plate owns the full nonlinear potential and mixed tangent. This adapter
//! owns only the time-image limits and physical force/observation units.
use std::sync::Arc;
use fs_plate::{ModePair,PlateError,shell::head::{TensionedDisk,dilation::DilationReduction,
    nonlinear::MembraneReductionBudget}};
use super::{ImpactError,MAX_IMPACT_MODES,invalid};

#[derive(Clone,Copy,Debug)]
pub struct DilationConfig {
    /// TOTAL effective inertia seen by radial edge displacement [kg].
    pub effective_mass_kg:f64,
    /// Extra hardware stiffness, in parallel with computed film stiffness [N/m].
    pub hardware_stiffness_n_m:f64,
    /// Small fixed-reference radial strain cap; must lie in (0,.05).
    pub maximum_abs_strain:f64,
    /// Maximum transverse slope; must lie in (0,1).
    pub maximum_slope:f64,
}
#[derive(Clone,Debug)]
pub struct DilationPotential {
    reduction:Arc<DilationReduction>,
    clock_bounds:Vec<f64>,
    maximum_abs_strain:f64,
    maximum_slope:f64,
}
#[derive(Clone,Copy,Debug,PartialEq)]
pub struct DilationObservation {
    pub rim_displacement_m:f64,
    pub edge_strain:f64,
    /// Uniform baseline tension, not total local vibration-induced stress.
    pub baseline_tension_n_m:f64,
    pub maximum_slope:f64,
}
impl DilationPotential {
    /// Preserve the original film basis and add one physical edge coordinate.
    /// Static preload is balanced at zero; external effort is incremental work.
    pub fn from_disk(film:&TensionedDisk,modes:&[ModePair],config:DilationConfig,
        budget:MembraneReductionBudget)->Result<Self,ImpactError> {
        if modes.is_empty() || modes.len()>=MAX_IMPACT_MODES
            || !config.maximum_abs_strain.is_finite() || !(0.0<config.maximum_abs_strain && config.maximum_abs_strain<0.05)
            || !config.maximum_slope.is_finite() || !(0.0<config.maximum_slope && config.maximum_slope<1.0) {
            return Err(invalid("actuated membrane needs bounded modes and explicit small strain/slope caps"));
        }
        let reduction=DilationReduction::new(film,modes,config.effective_mass_kg,
            config.hardware_stiffness_n_m,budget,|n,k,rhs| {
                let factor=fs_la::factor::cholesky(k,n).map_err(|_|PlateError::BadSection {
                    what:"actuated film's in-plane stiffness is not positive definite"})?;
                let mut out=rhs.to_vec();for u in &mut out {factor.solve(u);}Ok(out)
            }).map_err(|e|ImpactError::Owner(e.to_string()))?;
        let clock_bounds=reduction.clock_bounds(config.maximum_abs_strain)
            .map_err(|e|ImpactError::Owner(e.to_string()))?;
        Ok(Self {reduction:Arc::new(reduction),clock_bounds,
            maximum_abs_strain:config.maximum_abs_strain,maximum_slope:config.maximum_slope})
    }
    pub fn reduction(&self)->&DilationReduction {&self.reduction}
    pub fn clock_bounds(&self)->&[f64] {&self.clock_bounds}
    pub fn observe(&self,q:&[f64])->Result<DilationObservation,ImpactError> {
        let n=self.reduction.rim_coordinate();
        if q.len()!=n+1 || q.iter().any(|v|!v.is_finite()) {return Err(invalid("incomplete actuated membrane state"));}
        let edge_strain=self.reduction.edge_strain(q);
        let maximum_slope=self.reduction.membrane().maximum_slope(&q[..n]);
        let tension=self.reduction.baseline_tension_n_m(q);
        let energy=self.reduction.potential(q);
        if ![edge_strain,maximum_slope,tension,energy].iter().all(|v|v.is_finite())
            || edge_strain.abs()>self.maximum_abs_strain || maximum_slope>self.maximum_slope
            || tension<=0.0 || energy<0.0 {
            return Err(invalid("actuated membrane exceeds strain/slope or loses positive baseline tension"));
        }
        Ok(DilationObservation {rim_displacement_m:self.reduction.rim_displacement_m(q),
            edge_strain,baseline_tension_n_m:tension,maximum_slope})
    }
}
#[cfg(test)]
mod tests;
