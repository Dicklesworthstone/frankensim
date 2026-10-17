//! Geometry sampling uses the shared passive multiport fitter, not a second fit.
use super::{engine::radiation::{PortBasis,fit_projected},exterior_geometry::{Boundary,Specification,Samples},exterior_loading};
#[cfg(test)]
use super::engine::radiation::validate_complete_fit;
use fs_math::c64::C64;
pub use fs_couple::render::plate::impact::radiation::fit::{Fit,fit};

/// Exterior loaded playback fits a complete load with the shared fixed dictionary,
/// whose capacity stays at 32 ports; larger boards use explicit acoustic ports.
const MAX_COMPLETE_PORTS: usize = 32;

/// Solve one modal BEM batch per frequency. Those SAME fields feed the load
/// and receiver fits. Pressure/acceleration = pressure/velocity * i/omega.
pub fn prepare(boundary:&Boundary,spec:&Specification,mechanics_rate:u32)->Result<(Fit,Samples),String> {
    let r=boundary.weights.len();let omega=spec.omega();
    if r>MAX_COMPLETE_PORTS || omega.len()<33 {return Err("loaded playback requires <=32 complete board modes and >=33 frequencies".into());}
    if mechanics_rate<8_000 || 8.*omega[omega.len()-1]>=0.9*std::f64::consts::PI*f64::from(mechanics_rate) {
        return Err("loaded playback's off-band storage poles exceed the declared mechanical-rate guard".into());
    }
    let prepared=exterior_loading::LoadingSweep::new(boundary,spec,&omega)?;
    let delays_s=prepared.delays_s().to_vec();
    let mut matrices=Vec::with_capacity(omega.len());
    let mut values=vec![vec![Vec::with_capacity(omega.len());r];spec.receivers.len()];
    let mut minimum_ppw=f64::INFINITY;let mut maximum_condition_lower_bound=0.0_f64;
    for &w in &omega {
        let field=prepared.sample(w)?;
        matrices.push(field.impedance);minimum_ppw=minimum_ppw.min(field.minimum_ppw);
        maximum_condition_lower_bound=maximum_condition_lower_bound.max(field.condition_lower_bound);
        for (channel,row) in field.receiver_transfer.iter().enumerate() {
            for (input,&h) in row.iter().enumerate() {values[channel][input].push(h*C64::new(0.,1./w));}
        }
    }
    let fitted=fit(&omega,&matrices,r)?;
    Ok((fitted,Samples {omega,values,delays_s,minimum_ppw,maximum_condition_lower_bound}))
}

/// The fitted acoustic coordinates are independent of the complete structural
/// basis. Error fields in `fit` refer to the lifted COMPLETE BEM load, including
/// projection error, not just the retained acoustic matrix.
pub struct ProjectedFit {
    pub fit:Fit,
    pub basis:PortBasis,
    /// Relative area-weighted normal-velocity projection error. This geometric
    /// diagnostic does not replace the full complex and power acceptance gates.
    pub surface_rms_error:f64,
}

/// Retain this exterior geometry owner's existing panel ceiling when routing
/// geometry-only selection through the shared acoustic coordinate owner.
fn geometry_basis(weights:&[Vec<f64>],areas:&[f64],maximum:usize)
    ->Result<(PortBasis,f64),String> {
    super::engine::radiation::geometry_basis(weights,areas,maximum,
        super::exterior_geometry::MAX_PANELS)
}

/// Keep every structural coordinate and receiver input, fitting only the
/// explicitly bounded geometry-derived radiating subspace. Each frequency is
/// solved once by the existing complete modal BEM batch. An insufficient rank
/// refuses after full-load validation; it never silently drops a mechanical
/// mode, relaxes the resistance gate, or substitutes an output-only load.
pub fn prepare_projected(boundary:&Boundary,spec:&Specification,mechanics_rate:u32,
    max_ports:usize)->Result<(ProjectedFit,Samples),String> {
    let r=boundary.weights.len();let omega=spec.omega();
    if !(33..=257).contains(&omega.len()) || omega.len()%2==0
        || omega.iter().enumerate().any(|(i,w)|!w.is_finite() || *w<=0.
            || (i>0 && *w<=omega[i-1])) {
        return Err("projected loaded playback requires an odd 33..257-frequency increasing positive grid".into());
    }
    if mechanics_rate<8_000 || 8.*omega[omega.len()-1]>=0.9*std::f64::consts::PI*f64::from(mechanics_rate) {
        return Err("loaded playback's off-band storage poles exceed the declared mechanical-rate guard".into());
    }
    let (basis,surface_rms_error)=geometry_basis(&boundary.weights,boundary.surface.areas(),max_ports)?;
    let prepared=exterior_loading::LoadingSweep::new(boundary,spec,&omega)?;
    let delays_s=prepared.delays_s().to_vec();
    let mut matrices=Vec::with_capacity(omega.len());
    let mut values=vec![vec![Vec::with_capacity(omega.len());r];spec.receivers.len()];
    let mut minimum_ppw=f64::INFINITY;let mut maximum_condition_lower_bound=0.0_f64;
    for &w in &omega {
        let field=prepared.sample(w)?;
        matrices.push(field.impedance);minimum_ppw=minimum_ppw.min(field.minimum_ppw);
        maximum_condition_lower_bound=maximum_condition_lower_bound.max(field.condition_lower_bound);
        for (channel,row) in field.receiver_transfer.iter().enumerate() {
            for (input,&h) in row.iter().enumerate() {values[channel][input].push(h*C64::new(0.,1./w));}
        }
    }
    let fit=fit_projected(&omega,&matrices,&basis)?;
    Ok((ProjectedFit {fit,basis,surface_rms_error},
        Samples {omega,values,delays_s,minimum_ppw,maximum_condition_lower_bound}))
}

#[cfg(test)]
#[path="radiation_port_fit_tests.rs"]
mod port_tests;
