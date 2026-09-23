//! Geometry sampling uses the shared passive multiport fitter, not a second fit.
use super::{exterior_geometry::{Boundary,Specification,Samples},exterior_loading};
use fs_math::c64::C64;
pub use fs_couple::render::plate::impact::radiation::fit::{Fit,fit};
const MAX_PORTS:usize=32;

/// Solve one modal BEM batch per frequency. Those SAME fields feed the load
/// and receiver fits. Pressure/acceleration = pressure/velocity * i/omega.
pub fn prepare(boundary:&Boundary,spec:&Specification,mechanics_rate:u32)->Result<(Fit,Samples),String> {
    let r=boundary.weights.len();let omega=spec.omega();
    if r>MAX_PORTS || omega.len()<33 {return Err("loaded playback requires <=32 complete board modes and >=33 frequencies".into());}
    if mechanics_rate<8_000 || 8.*omega[omega.len()-1]>=0.9*std::f64::consts::PI*f64::from(mechanics_rate) {
        return Err("loaded playback's off-band storage poles exceed the declared mechanical-rate guard".into());
    }
    let mut matrices=Vec::with_capacity(omega.len());
    let mut values=vec![vec![Vec::with_capacity(omega.len());r];spec.receivers.len()];
    let mut minimum_ppw=f64::INFINITY;let mut maximum_condition_lower_bound=0.0_f64;
    for &w in &omega {
        let field=exterior_loading::sample(boundary,spec,w)?;
        matrices.push(field.impedance);minimum_ppw=minimum_ppw.min(field.minimum_ppw);
        maximum_condition_lower_bound=maximum_condition_lower_bound.max(field.condition_lower_bound);
        for (channel,row) in field.receiver_transfer.iter().enumerate() {
            for (input,&h) in row.iter().enumerate() {values[channel][input].push(h*C64::new(0.,1./w));}
        }
    }
    let fitted=fit(&omega,&matrices,r)?;
    let delays_s=spec.receivers.iter().map(|p| {
        let radius=(0..3).fold(0.0_f64,|a,j|a.hypot(p[j]-boundary.center[j]));
        (radius-boundary.radius)/spec.medium.sound_speed
    }).collect();
    Ok((fitted,Samples {omega,values,delays_s,minimum_ppw,maximum_condition_lower_bound}))
}
