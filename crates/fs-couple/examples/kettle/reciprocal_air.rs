//! Optional variational exterior: the SAME P1/P0 fields do work and reach the mic.
//! Retains the collocation comparison; never retries a refused solve with another
//! formulation, repairs a pressure matrix, or feeds mean traces to a P0 observer.
use super::{Bake, Boundary, Error, Settings, Spec, Medium, C64, OUTPUT_RATE,
    admit_power, observer, phase, fit, radiation};
use fs_bem::{near_field::Options, panel3d::SpherePanels,
    reciprocal_radiation::{Budget, Prepared}};
use fs_exec::CancelGate;
use fs_plate::{ModePair, shell::head::TensionedDisk};
use core::f64::consts::TAU;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Controls {
    pub angular_refinement: u32,
    pub maximum_sweep_work: u64,
}
impl Controls {
    fn validate(self) -> Result<(), Error> {
        if self.angular_refinement>2 || self.maximum_sweep_work==0
            || self.maximum_sweep_work>10_000_000_000_000 {
            return Err("reciprocal air requires angular refinement 0..2 and 1..10000000000000 sweep work units".into());
        }
        Ok(())
    }
}
/// Consume just this option, only after validating both controls. Work units are
/// the existing solver's estimate, NOT milliseconds, flops or a callback budget.
pub fn option(args: &mut Vec<String>) -> Result<Option<Controls>, Error> {
    let mut positions=args.iter().enumerate().filter(|(_,a)|a.as_str()=="--reciprocal-air");
    let Some((i,_))=positions.next() else { return Ok(None); };
    if positions.next().is_some() || args.len()-i<3 {
        return Err("supply --reciprocal-air ANGULAR_REFINEMENT MAX_SWEEP_WORK once".into());
    }
    let controls=Controls {angular_refinement:args[i+1].parse()?,maximum_sweep_work:args[i+2].parse()?};
    controls.validate()?;
    args.drain(i..i+3);Ok(Some(controls))
}
struct Samples { impedances:Vec<Vec<C64>>, pressure:Vec<Vec<C64>>, delay_s:f64 }

// Kept at the sampling seam so focused tests can use actual coarse BEM geometry
// without constructing a synthetic substitute observer or a full drum eigensolve.
fn sample(boundary:Boundary,medium:Medium,omega:&[f64],receiver:[f64;3],controls:Controls,
    gate:&CancelGate)->Result<Samples,Error> {
    controls.validate()?;
    if gate.is_requested() {return Err("reciprocal air preparation cancelled".into());}
    let n=boundary.weights.len();
    if n==0 || n>32 || omega.is_empty() || omega.len()>257
        || omega.iter().enumerate().any(|(i,w)|!w.is_finite() || *w<=0.0 || i>0 && *w<=omega[i-1])
        || !medium.density.is_finite() || medium.density<=0.0
        || !medium.sound_speed.is_finite() || medium.sound_speed<=0.0
        || receiver.iter().any(|x|!x.is_finite())
        || !boundary.radius_m.is_finite() || boundary.radius_m<=0.0
        || boundary.weights.iter().any(|r|r.len()!=boundary.triangles.len() || r.iter().any(|v|!v.is_finite())) {
        return Err("reciprocal air needs complete finite source rows, medium, receiver and increasing frequencies".into());
    }
    let k_max=omega[omega.len()-1]/medium.sound_speed;
    let surface=SpherePanels::from_triangles(boundary.triangles)?;
    let per_frequency=controls.maximum_sweep_work/omega.len() as u64;
    let source=Prepared::new(&surface,k_max,Budget {maximum_unknowns:4096,
        maximum_cache_bytes:256*1024*1024,maximum_work:per_frequency,
        angular_refinement:controls.angular_refinement},||gate.is_requested())?;
    let work=source.velocity_port_work(n)?;
    if work>per_frequency {return Err("complete reciprocal source/projection sweep exceeds the supplied work budget".into());}
    let scene=source.receivers(&[receiver],64.0*f64::EPSILON*boundary.radius_m)?;
    // No enclosing-sphere exclusion: use actual triangle clearance. A very
    // close mic retains full propagation in the fitted transfer instead of
    // forcing an unsupported sub-two-sample explicit delay.
    let travel=scene.clearances_m()[0]/medium.sound_speed;
    if !travel.is_finite() || travel>2.0 {return Err("reciprocal microphone travel exceeds two seconds".into());}
    let delay_s=if travel.next_down()>=2.0/f64::from(OUTPUT_RATE) {travel.next_down()} else {0.0};
    let mut values=vec![vec![C64::ZERO;omega.len()];n];
    let mut impedances=vec![vec![C64::ZERO;n*n];omega.len()];
    let rows:Vec<_>=boundary.weights.iter().map(Vec::as_slice).collect();
    let (mut residual,mut flux,mut quadrature,mut ppw)=(0.0_f64,0.0_f64,0.0_f64,f64::INFINITY);
    for index in (0..omega.len()).rev() {
        if gate.is_requested() {return Err("reciprocal air sweep cancelled".into());}
        let w=omega[index];let k=w/medium.sound_speed;
        // Receiver rows are assembled once per frequency, shared across sources.
        let observation=scene.prepare(k,medium,Options::default(),false,||gate.is_requested())?;
        let batch=source.solve_velocity_ports(k,medium,&rows,||gate.is_requested())?;
        admit_power(batch.impedance(),n)?;
        impedances[index].copy_from_slice(batch.impedance());
        for (j,solution) in batch.solutions().iter().enumerate() {
            if gate.is_requested() {return Err("reciprocal microphone evaluation cancelled".into());}
            let observed=observation.evaluate(solution)?;
            // The library solves UNIT VELOCITY. Pressure observation preserves
            // P1/P0 densities, THEN converts p/v to p/a with 1/(-i*w)=i/w.
            values[j][index]=phase(observed.pressure[0]*C64::new(0.0,1.0/w),-w*delay_s);
            if !values[j][index].re.is_finite() || !values[j][index].im.is_finite() {
                return Err("reciprocal acceleration pressure overflow".into());
            }
            residual=residual.max(solution.backward_residual);
            flux=flux.max(solution.flux_balance);
            quadrature=quadrature.max(observed.quadrature_error_estimate_pa[0]/w);
            ppw=ppw.min(solution.panels_per_wavelength);
        }
    }
    if gate.is_requested() {return Err("reciprocal air cancelled before sample publication".into());}
    eprintln!("reciprocal kettle air: unknowns={}, cache_bytes={}, sweep_work={}, angular_refinement={}, min_panels_per_wavelength={ppw}, backward_residual={residual}, flux_defect={flux}, maximum_acceleration_receiver_quadrature_estimate={quadrature}, microphone_clearance_m={}, explicit_delay_s={delay_s}; fixed P1/P0 reference scene, not calibrated or continuous-band certified",
        source.unknowns(),source.cache_bytes(),work*omega.len() as u64,controls.angular_refinement,scene.clearances_m()[0]);
    Ok(Samples {impedances,pressure:values,delay_s})
}

pub fn bake(spec:&Spec,film:&TensionedDisk,modes:&[ModePair],settings:Settings,
    controls:Controls,gate:&CancelGate)->Result<(radiation::Model,Bake),Error> {
    settings.validate(spec)?;controls.validate()?;
    if modes.iter().any(|m|m.lambda.sqrt()/TAU<settings.band_hz[0] || m.lambda.sqrt()/TAU>settings.band_hz[1]) {
        return Err("reciprocal acoustic band must cover all retained head modes".into());
    }
    let boundary=Boundary::new(spec,film,modes)?;
    let count=2*settings.intervals+1;
    let omega:Vec<_>=(0..count).map(|i|TAU*(settings.band_hz[0]
        +(settings.band_hz[1]-settings.band_hz[0])*i as f64/(count-1) as f64)).collect();
    let samples=sample(boundary,Medium{density:spec.air[0],sound_speed:spec.air[1]},
        &omega,settings.receiver_m,controls,gate)?;
    if gate.is_requested() {return Err("reciprocal load fit cancelled".into());}
    let fitted=fit::fit(&omega,&samples.impedances,modes.len())?;
    let mut filters=Vec::with_capacity(modes.len());
    for row in &samples.pressure {
        if gate.is_requested() {return Err("reciprocal pressure fit cancelled".into());}
        filters.push(observer(&omega,row,settings.order)?);
    }
    if gate.is_requested() {return Err("reciprocal fit cancelled before publication".into());}
    eprintln!("reciprocal kettle fit: poles={}, complex_peak/RMS={}/{}, resistance_peak/RMS={}/{}; independent pressure audit and complete mechanical ports retained",
        fitted.model.poles.len(),fitted.peak_error,fitted.rms_error,fitted.resistance_peak_error,fitted.resistance_rms_error);
    Ok((fitted.model,Bake{filters,delay_s:samples.delay_s}))
}

#[cfg(test)]
#[path="reciprocal_air_tests.rs"]
mod tests;
