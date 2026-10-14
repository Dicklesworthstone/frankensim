//! The existing reciprocal image, with shared source work and near receivers.
//! No second CLI selection or alternate acoustic equation on failure.
use super::{Error, Medium, C64, OUTPUT_RATE, admit_power, phase, field};
use fs_bem::{near_field::Options, panel3d::SpherePanels};
use fs_exec::CancelGate;

pub(super) struct Samples {
    pub impedances:Vec<Vec<C64>>,
    pub pressure:Vec<Vec<C64>>,
    pub delay_s:f64,
    pub minimum_ppw:f64,
}

/// Cold geometry -> one complete source batch per frequency -> work AND pressure.
/// Separate from fitting so regressions can exercise the actual BEM owner.
#[allow(clippy::too_many_arguments)]
pub(super) fn sample(surface:&SpherePanels,weights:&[Vec<f64>],radius_m:f64,
    medium:Medium,omega:&[f64],receiver:[f64;3],maximum_work:u64,
    gate:&CancelGate)->Result<Samples,Error> {
    if gate.is_requested() {return Err("reciprocal air preparation cancelled".into());}
    let n=weights.len();
    if n==0 || n>32 || omega.is_empty() || omega.len()>257
        || omega.iter().enumerate().any(|(i,w)|!w.is_finite() || *w<=0.0 || i>0 && *w<=omega[i-1])
        || !medium.density.is_finite() || medium.density<=0.0
        || !medium.sound_speed.is_finite() || medium.sound_speed<=0.0
        || receiver.iter().any(|x|!x.is_finite())
        || !radius_m.is_finite() || radius_m<=0.0
        || weights.iter().any(|r|r.len()!=surface.areas().len() || r.iter().any(|v|!v.is_finite())) {
        return Err("reciprocal air needs complete finite rows, medium, receiver and increasing frequencies".into());
    }
    let k_max=omega[omega.len()-1]/medium.sound_speed;
    let source=field::prepare(surface,k_max,n,omega.len(),maximum_work,gate)?;
    let scene=source.receivers(&[receiver],64.0*f64::EPSILON*radius_m)?;
    // Actual triangle clearance and exterior winding replace the bounding-sphere
    // exclusion. The fixed positive clearance is geometric, not fitted latency.
    let travel=scene.clearances_m()[0]/medium.sound_speed;
    if !travel.is_finite() || travel>2.0 {return Err("reciprocal microphone travel exceeds two seconds".into());}
    // The delay owner requires >=2 samples. Keep the entire propagation in the
    // transfer for closer microphones, rather than inventing extra transit time.
    let delay_s=if travel.next_down()>=2.0/f64::from(OUTPUT_RATE) {travel.next_down()} else {0.0};
    let mut pressure=vec![vec![C64::ZERO;omega.len()];n];
    let mut impedances=vec![vec![C64::ZERO;n*n];omega.len()];
    let rows:Vec<_>=weights.iter().map(Vec::as_slice).collect();
    let (mut residual,mut flux,mut quadrature,mut ppw)=(0.0_f64,0.0_f64,0.0_f64,f64::INFINITY);
    for index in (0..omega.len()).rev() {
        if gate.is_requested() {return Err("reciprocal air sweep cancelled".into());}
        let w=omega[index];let k=w/medium.sound_speed;
        // Prepare each physical receiver ONCE per frequency, not once per mode.
        let observation=scene.prepare(k,medium,Options::default(),false,||gate.is_requested())?;
        let batch=source.solve_velocity_ports(k,medium,&rows,||gate.is_requested())?;
        admit_power(batch.impedance(),n)?;
        impedances[index].copy_from_slice(batch.impedance());
        for (j,solution) in batch.solutions().iter().enumerate() {
            if gate.is_requested() {return Err("reciprocal microphone evaluation cancelled".into());}
            let observed=observation.evaluate(solution)?;
            // Z is velocity-to-force. Only the receiver transfer is converted
            // from p/v to p/a, AFTER observing the original P1/P0 densities.
            pressure[j][index]=phase(observed.pressure[0]*C64::new(0.0,1.0/w),-w*delay_s);
            if !pressure[j][index].re.is_finite() || !pressure[j][index].im.is_finite() {
                return Err("reciprocal acceleration pressure overflow".into());
            }
            residual=residual.max(solution.backward_residual);
            flux=flux.max(solution.flux_balance);
            quadrature=quadrature.max(observed.quadrature_error_estimate_pa[0]/w);
            ppw=ppw.min(solution.panels_per_wavelength);
        }
    }
    if gate.is_requested() {return Err("reciprocal air cancelled before sample publication".into());}
    eprintln!("reciprocal kettle receiver: minimum_ppw={ppw}, backward_residual={residual}, flux_defect={flux}, maximum_acceleration_quadrature_estimate={quadrature}, clearance_m={}, explicit_delay_s={delay_s}; actual P1/P0 field, no extra range gain or calibrated/full-band claim",scene.clearances_m()[0]);
    Ok(Samples {impedances,pressure,delay_s,minimum_ppw:ppw})
}

#[cfg(test)]
#[path="reciprocal_air_tests.rs"]
mod tests;
