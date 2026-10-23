//! Feedback and observation of one stationary acoustic scene, never an output EQ.
use super::*;
use crate::Mechanics;
pub(super) use fs_couple::render::plate::impact::radiation::{Model,fit};
use fs_couple::render::plate::impact::radiation::RadiationObservation;

pub fn option(args:&mut Vec<String>)->Result<bool,Error> {
    let count=args.iter().filter(|s|s.as_str()=="--radiation-feedback").count();
    if count>1 {return Err("--radiation-feedback may be supplied only once".into());}
    args.retain(|s|s!="--radiation-feedback");Ok(count==1)
}
pub fn admit_command(selected:bool,command:&str,has_neck:bool)->Result<(),Error> {
    if selected && (!matches!(command,"splash-wav"|"splash-mic"|"drum-wav"|"drum-mic"|
        "drum-stretch-wav"|"drum-stretch-mic"|"snare-wav"|"snare-mic"|"snare-off-wav"|"snare-off-mic"|
        "hihat-wav"|"hihat-mic") || (has_neck && !matches!(command,
            "drum-wav"|"drum-mic"|"drum-stretch-wav"|"drum-stretch-mic"|
            "snare-wav"|"snare-mic"|"snare-off-wav"|"snare-off-mic"))) {
        return Err("radiation feedback needs a nonlinear-capable pressure command; necks require the coupled vent selection and internal-only impedance".into());
    }
    Ok(())
}

pub(super) fn admit_boundary(boundary:&Boundary,spec:radiation_spec::Spec)->Result<(),Error> {
    spec.validate()?;
    if boundary.weights.is_empty() || boundary.weights.len()>32 || !(8..=64).contains(&spec.training_intervals)
        || 8.*core::f64::consts::TAU*spec.band_hz[1]*MECHANICAL_DT>=0.9*core::f64::consts::PI {
        return Err("feedback requires the complete <=32-source boundary, 33..257 samples, and resolved off-band acoustic poles".into());
    }
    Ok(())
}

// Diagonal radiation-power checks do NOT establish multiport passivity. Admit
// the real representation of the full Hermitian power form through fs-phs's
// existing numerical PSD check. Never clip, symmetrize or replace the raw Z.
fn admit_power(z:&[C64],n:usize)->Result<(),Error> {
    if !(1..=32).contains(&n) || z.len()!=n*n || z.iter().any(|v|!v.re.is_finite()||!v.im.is_finite()) {
        return Err("nonfinite or incomplete radiation impedance".into());
    }
    let m=2*n;let mut h=vec![0.;m*m];
    for i in 0..n {for j in 0..n {
        let v=(z[i*n+j]+z[j*n+i].conj()).scale(0.5);
        h[i*m+j]=v.re;h[(i+n)*m+j+n]=v.re;
        h[i*m+j+n]=-v.im;h[(i+n)*m+j]=v.im;
    }}
    fs_phs::QuadraticStorage::new(h,m).map_err(|e|format!("BEM multiport power form refused; refine the acoustic boundary: {e}"))?;
    Ok(())
}

/// Acceleration-driven BEM p_j = i/omega * p_velocity_j. Integrate the
/// conjugate surface work and multiply by -i*omega to obtain resisting Z.
/// A real shape row is used once to drive velocity and once to transpose work.
pub(super) fn project(surface:&SpherePanels,rows:&[Vec<f64>],solutions:&[RadiationSolution],omega:f64)
    ->Result<Vec<C64>,Error> {
    let traces:Vec<&[C64]>=solutions.iter().map(|s|s.pressure.as_slice()).collect();
    project_pressures(surface,rows,&traces,omega)
}
// Both acoustic images return panel mean pressure for work against the same
// P0 velocity rows. Only the work projection is shared; their field observers
// must retain their own basis. The full Hermitian admission is unchanged.
pub(super) fn project_pressures(surface:&SpherePanels,rows:&[Vec<f64>],pressures:&[&[C64]],omega:f64)
    ->Result<Vec<C64>,Error> {
    let n=rows.len();let panels=surface.areas().len();
    if !(1..=32).contains(&n) || pressures.len()!=n || !omega.is_finite() || omega<=0.
        || rows.iter().any(|b|b.len()!=panels||b.iter().any(|x|!x.is_finite()))
        || pressures.iter().any(|p|p.len()!=panels) {
        return Err("radiation load projection lost its full physical boundary basis".into());
    }
    let mut z=vec![C64::ZERO;n*n];
    for (j,pressure) in pressures.iter().enumerate() {for (i,row) in rows.iter().enumerate() {
        let mut force=C64::ZERO;
        for ((p,area),b) in pressure.iter().zip(surface.areas()).zip(row) {force=force+p.scale(area*b);}
        z[i*n+j]=force*C64::new(0.,-omega);
    }}
    admit_power(&z,n)?;Ok(z)
}

// All fields are private: a loaded observer can only be obtained together with
// its admitted load, from one BEM scene, before numerical preparation/playing.
pub(crate) struct Prepared {
    baked:Vec<Bake>,
    sources:Vec<usize>,
}
fn admit_instrument(e:&Experiment)->Result<(),Error> {
    match &e.system {
        Mechanics::Reference(s) if s.samples()==0 && s.radiation_observation().is_none()=>{},
        _=>return Err("attach radiation before playing/numerical preparation; linear-only modal/snare images refuse instead of silently converting".into()),
    }
    if let Some(air)=&e.air {
        let boundary=e.acoustics.as_ref().ok_or("missing exterior feedback boundary")?;
        // This refuses effective-length necks, omitted/duplicate mouth columns,
        // and accidental radiation from interior standing-wave coordinates.
        air.coupling.admit_exterior_sources(&boundary.state_modes)?;
    }
    Ok(())
}
// Cold attachment uses the same map that produced both impedance and pressure.
// Kept separate from fitting so the mechanics seam can also be checked against
// an explicitly declared passive test load, without substituting one in playback.
fn attach_model(mut e:Experiment,model:&Model,sources:&[usize])->Result<Experiment,Error> {
    admit_instrument(&e)?;
    if e.acoustics.as_ref().is_none_or(|b|b.state_modes.as_slice()!=sources) {
        return Err("radiation load lost the source map used for boundary work".into());
    }
    let Mechanics::Reference(system)=e.system else {unreachable!("admitted before preparation")};
    let system=match &e.air {
        Some(air)=>air.coupling.attach_exterior_radiation(system,model,sources,256)?,
        None=>system.with_radiation_load(model,sources,256)?,
    };
    e.system=Mechanics::Reference(system);
    Ok(e)
}
pub(crate) fn prepare(e:Experiment,frames:usize,scale:f64,receivers:&[Receiver],
    spec:radiation_spec::Spec,gate:&CancelGate)->Result<(Experiment,Prepared),Error> {
    admit_instrument(&e)?;admit_render(&e,frames,scale)?;
    // Feedback must respect the SAME compact-neck chart as one-way observation,
    // including user-selected bands beyond the aperture builder's default band.
    spec.admit_necks(&e)?;
    let boundary=e.acoustics.as_ref().ok_or("missing feedback boundary")?;
    let sources=boundary.state_modes.clone();
    let reciprocal=e.air.as_ref().is_some_and(|a|a.coupling.neck_count()>0);
    let (baked,model)=bake_scene_image(boundary,receivers,spec,e.ambient.medium(),gate,true,reciprocal)?;
    if gate.is_requested() {return Err("radiation feedback cancelled before attachment".into());}
    let model=model.ok_or("missing admitted feedback model")?;
    let e=attach_model(e,&model,&sources)?;
    Ok((e,Prepared{baked,sources}))
}
fn observation(s:&Mechanics)->Option<RadiationObservation> {
    match s {Mechanics::Reference(s)=>s.radiation_observation(),Mechanics::Nonlinear(s)=>s.radiation_observation(),
        Mechanics::Substepped(s)=>s.radiation_observation(),Mechanics::Driven{inner,..}=>observation(inner),
        Mechanics::Prepared(_)=>None}
}
impl Prepared {
    pub(crate) fn render(self,e:&mut Experiment,frames:usize,scale:f64,gate:&CancelGate)->Result<Vec<u8>,Error> {
        if e.acoustics.as_ref().is_none_or(|b|b.state_modes!=self.sources) {
            return Err("feedback observer lost its mechanical source addresses".into());
        }
        let wav=render_baked_with_gate(e,frames,scale,&self.baked,gate)?;
        if let Some(air)=observation(&e.system) {
            eprintln!("radiation feedback: poles={}, acoustic_storage_j={}, endpoint_dissipation_W={}; ALREADY included in total mechanics; stationary linear acoustic fit, not full-band calibration",
                air.poles,air.stored_energy_j,air.dissipated_power_w);
        }
        Ok(wav)
    }
}

#[cfg(test)]
#[path="radiation_feedback_tests.rs"]
mod tests;

#[cfg(test)]
#[path="radiation_feedback_vent_tests.rs"]
mod vent_tests;
