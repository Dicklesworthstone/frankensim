//! One source, one mechanical/decimation clock, one or two physical receivers.
//! BEM boundary solves are shared; receiver response, fitting, delays and filter
//! histories remain independent. Mono uses this same path with one receiver.
use super::*;
use fs_couple::pcm_wav::encode_pcm16_wav_interleaved;

#[path = "observer_fit.rs"]
mod observer_fit;

#[path = "radiation_spec.rs"]
pub mod radiation_spec;

#[path = "radiation_feedback.rs"]
pub mod feedback;

/// Same comma-separated SI position syntax as grand_piano. Remove only this
/// option, transactionally: bad/duplicate input leaves all arguments unchanged.
pub fn option(args:&mut Vec<String>)->Result<Option<[f64;3]>,Error> {
    let mut matches=args.iter().enumerate().filter(|(_,a)|a.as_str()=="--microphone-right");
    let Some((index,_))=matches.next() else {return Ok(None);};
    if matches.next().is_some() {return Err("--microphone-right may be supplied only once".into());}
    let text=args.get(index+1).ok_or("--microphone-right needs x_m,y_m,z_m")?;
    let values=text.split(',').map(str::parse::<f64>).collect::<Result<Vec<_>,_>>()?;
    if values.len()!=3 || values.iter().any(|x|!x.is_finite()) {
        return Err("--microphone-right needs three finite comma-separated metre coordinates".into());
    }
    let position=[values[0],values[1],values[2]];
    args.drain(index..index+2);
    Ok(Some(position))
}
pub fn admit_command(right:Option<[f64;3]>,command:&str)->Result<(),Error> {
    if right.is_some() && !matches!(command,"splash-mic"|"drum-mic"|"drum-stretch-mic"|
        "drum-modal-mic"|"snare-mic"|"snare-off-mic") {
        return Err("--microphone-right requires a finite-point -mic command; no implicit CSV/far-field conversion".into());
    }
    Ok(())
}

fn bake_receivers(boundary:&Boundary,receivers:&[Receiver])->Result<Vec<Bake>,Error> {
    bake_receivers_with_spec(boundary,receivers,radiation_spec::Spec::default(),&CancelGate::new_clock_free())
}

fn bake_receivers_with_spec(boundary:&Boundary,receivers:&[Receiver],spec:radiation_spec::Spec,
    gate:&CancelGate)->Result<Vec<Bake>,Error> {
    Ok(bake_scene_with_spec(boundary,receivers,spec,gate,false)?.0)
}

// The same BEM fields feed both pressure observers and (when selected) the
// complete force/velocity matrix. One-way callers never construct a load.
fn bake_scene_with_spec(boundary:&Boundary,receivers:&[Receiver],spec:radiation_spec::Spec,
    gate:&CancelGate,load:bool)->Result<(Vec<Bake>,Option<feedback::Model>),Error> {
    if gate.is_requested() {return Err("radiation preparation cancelled".into());}
    if load {feedback::admit_boundary(boundary,spec)?;}
    let count=boundary.weights.len();let panels=boundary.triangles.len();
    let (prepared_panels,work)=spec.work(panels,count,receivers.len())?;
    if !(1..=2).contains(&receivers.len()) || count==0 || count>MAX_INPUTS
        || count!=boundary.state_modes.len() || panels==0 || panels>spec.max_panels
        || boundary.weights.iter().any(|b|b.len()!=panels || b.iter().any(|v|!v.is_finite())) {
        return Err("exterior radiation exceeds its receiver/shape/work contract".into());
    }
    let radius=boundary.triangles.iter().flatten().map(|p|p.iter().map(|x|x*x).sum::<f64>().sqrt())
        .fold(0.0_f64,f64::max);
    let medium=Medium::air();let dt=1.0/f64::from(OUTPUT_RATE);
    let refined=curved_aperture::uniform_refinement(boundary,spec.subdivisions,spec.max_panels,gate)?;
    let surface=SpherePanels::from_triangles(refined.triangles.clone())?;
    // Admit EVERY receiver against the actual boundary before a source solve.
    // Legacy finite/far fields retain their original formulas and timing.
    let observations=receivers::Scene::new(&surface,receivers,radius,medium,dt,gate)?;
    let mut baked=Vec::with_capacity(receivers.len());
    for (channel,&receiver) in receivers.iter().enumerate() {
        let range_m=receiver.position().iter().map(|x|x*x).sum::<f64>().sqrt();
        baked.push(Bake {filters:Vec::with_capacity(count),range_m,medium,
            propagation_delay_s:observations.delay(channel),pressure_gain:observations.gain(channel)});
    }
    let omega=spec.frequencies()?;
    let mut values=vec![vec![vec![C64::new(0.0,0.0);omega.len()];count];receivers.len()];
    let mut impedances=if load {vec![vec![C64::ZERO;count*count];omega.len()]}else{Vec::new()};
    eprintln!("radiation preparation: source_panels={panels}, prepared_panels={prepared_panels}, subdivisions={}, dense_work_units={work}, band_hz={:?}; same polyhedral source geometry, no mechanical refinement",spec.subdivisions,spec.band_hz);
    let mut ppw=f64::INFINITY;let mut condition=0.0_f64;
    // One formulation for the entire transfer; stitching different discrete
    // operators at kR=0.5 creates a numerical jump that a causal fit cannot fix.
    let formulation=if omega[omega.len()-1]*radius/medium.sound_speed<0.5 {
        Formulation::PlainCbie
    }else{Formulation::BurtonMiller};
    // Solve the highest frequency first: the BEM owner's wavelength guard can
    // reject an underresolved band before spending work on its lower samples.
    for index in (0..omega.len()).rev() {
        if gate.is_requested() {return Err("radiation preparation cancelled".into());}
        let w=omega[index];let k=w/medium.sound_speed;
        // Independent bounded receiver quadrature, prepared ONCE per frequency.
        let observation=observations.prepare(k,gate)?;
        let fields=acceleration_fields(&refined.weights,w);
        let field_refs:Vec<&[C64]>=fields.iter().map(Vec::as_slice).collect();
        // Geometry, factorization and all source-mode solves are independent of
        // the observation point. Do not repeat them per microphone.
        let solutions=solve_radiation_batch(&surface,k,medium,&field_refs,formulation)?;
        if solutions.len()!=count {return Err("BEM source batch lost a mechanical input".into());}
        for (input,solution) in solutions.iter().enumerate() {
            if !solution.radiated_power_roundoff_interval.1.is_finite()
                || !solution.panels_per_wavelength.is_finite()
                || !solution.condition_lower_bound.is_finite()
                || solution.radiated_power_roundoff_interval.1<0.0 {
                return Err("BEM reports invalid diagnostics or negative radiation power beyond roundoff; refine the acoustic solve".into());
            }
            ppw=ppw.min(solution.panels_per_wavelength);condition=condition.max(solution.condition_lower_bound);
            for channel in 0..receivers.len() {
                values[channel][input][index]=observation.response(channel,solution)?;
            }
        }
        if load {impedances[index]=feedback::project(&surface,&refined.weights,&solutions,w)?;}
    }
    let fitted=if load {
        if gate.is_requested() {return Err("radiation load fitting cancelled".into());}
        let f=feedback::fit::fit(&omega,&impedances,count)?;
        eprintln!("passive radiation fit: ports={}, poles={}, complex_peak/RMS={}/{}, resistance_peak/RMS={}/{}; complete signed matrices, independent held-out samples, no rank truncation",
            count,f.model.poles.len(),f.peak_error,f.rms_error,f.resistance_peak_error,f.resistance_rms_error);
        Some(f.model)
    }else{None};
    for (channel,bake) in baked.iter_mut().enumerate() {
        let mut maximum=0.0_f64;let mut rms=0.0_f64;
        for row in &values[channel] {
            if gate.is_requested() {return Err("radiation fitting cancelled".into());}
            let (filter,report)=observer_fit::fit(&omega,row,dt,spec.max_order)?;
            maximum=maximum.max(report.audit_maximum);rms=rms.max(report.audit_rms);
            eprintln!("observer fit: selected_order={}, attempts={}, selection_max={}, selection_rms={}, independent_audit_max={}, independent_audit_rms={}",
                report.order,report.attempts,report.selection_maximum,report.selection_rms,report.audit_maximum,report.audit_rms);
            bake.filters.push(filter);
        }
        eprintln!("receiver {channel} BEM bake: panels={prepared_panels}, inputs={count}, band_hz={:?}, training={}, order_selection={}, independent_audit={}, min_panels_per_wavelength={ppw}, condition_lower_bound_max={condition}, max_error={maximum}, worst_input_rms={rms}",
            spec.band_hz,spec.training_intervals+1,spec.training_intervals,2*spec.training_intervals);
        eprintln!("receiver={:?}; independent causal filters from shared source solves; linear undeformed acoustics, passive_feedback={load}; propagation_delay_s={}, pressure_gain={}",
            receivers[channel],bake.propagation_delay_s,bake.pressure_gain);
    }
    Ok((baked,fitted))
}

fn admit_render(experiment:&Experiment,frames:usize,full_scale_pa:f64)->Result<(),Error> {
    let boundary=experiment.acoustics.as_ref().ok_or("missing acoustic boundary")?;
    if frames==0 || frames>480000 || !full_scale_pa.is_finite() || full_scale_pa<=0.0
        || boundary.state_modes.is_empty() || boundary.state_modes.len()>MAX_INPUTS
        || boundary.state_modes.iter().any(|&k|k>=experiment.force.len()
            || k>=experiment.system.state().len()/2) {
        return Err("pressure export requires bounded frames, finite positive full-scale and valid state addresses".into());
    }
    Ok(())
}

pub(crate) fn render_receivers(experiment:&mut Experiment,frames:usize,full_scale_pa:f64,
    receivers:&[Receiver])->Result<Vec<u8>,Error> {
    render_receivers_with_spec(experiment,frames,full_scale_pa,receivers,
        radiation_spec::Spec::default(),&CancelGate::new_clock_free())
}

pub(crate) fn render_receivers_with_spec(experiment:&mut Experiment,frames:usize,full_scale_pa:f64,
    receivers:&[Receiver],spec:radiation_spec::Spec,gate:&CancelGate)->Result<Vec<u8>,Error> {
    admit_render(experiment,frames,full_scale_pa)?;
    spec.admit_necks(experiment)?;
    let baked=bake_receivers_with_spec(experiment.acoustics.as_ref().unwrap(),receivers,spec,gate)?;
    if gate.is_requested() {return Err("radiation preparation cancelled before mechanics".into());}
    render_baked_with_gate(experiment,frames,full_scale_pa,&baked,gate)
}

// Private so callers cannot replace geometry/BEM admission with an authored
// transfer. Tests can isolate frame plumbing with explicitly labeled fixtures.
fn render_baked(experiment:&mut Experiment,frames:usize,full_scale_pa:f64,baked:&[Bake])->Result<Vec<u8>,Error> {
    render_baked_with_gate(experiment,frames,full_scale_pa,baked,&CancelGate::new_clock_free())
}
fn render_baked_with_gate(experiment:&mut Experiment,frames:usize,full_scale_pa:f64,
    baked:&[Bake],gate:&CancelGate)->Result<Vec<u8>,Error> {
    admit_render(experiment,frames,full_scale_pa)?;
    let boundary=experiment.acoustics.as_ref().unwrap();
    let count=boundary.state_modes.len();let channels=baked.len();
    if !(1..=2).contains(&channels) || baked.iter().any(|b|b.filters.len()!=count
        || b.filters.iter().any(|f|f.t_s.to_bits()!=(1.0/f64::from(OUTPUT_RATE)).to_bits())) {
        return Err("receiver bank must retain the full source basis and output clock".into());
    }
    let mut observers=baked.iter().map(Bake::runtime).collect::<Result<Vec<_>,_>>()?;
    let mut decimator=Decimator::new(SUBSTEPS,count)?;
    let mut block=vec![0.0;SUBSTEPS*count];
    let x=experiment.system.state();
    let mut previous:Vec<f64>=boundary.state_modes.iter().map(|&k|x[2*k+1]).collect();
    let mut pressure=Vec::with_capacity(frames*channels);
    for _ in 0..frames {
        for frame in block.chunks_exact_mut(count) {
            experiment.system.step(&experiment.force,gate)?;
            let x=experiment.system.state();
            for ((sample,previous),&k) in frame.iter_mut().zip(&mut previous).zip(&boundary.state_modes) {
                let velocity=x[2*k+1];*sample=(velocity-*previous)/MECHANICAL_DT;*previous=velocity;
            }
        }
        // ONE complete mechanics block and decimator state for every receiver.
        // Neither observation can advance contact, a force program or gas time.
        let acceleration=decimator.preview(&block)?;
        let mut frame=[0.0;2];
        for (value,observer) in frame.iter_mut().zip(&mut observers) {*value=observer.step(acceleration)?;}
        decimator.commit();pressure.extend_from_slice(&frame[..channels]);
    }
    // This offline front door publishes nothing on any mechanics/observer/PCM
    // error. Accepted mechanics may have advanced: discard on failure, no retry.
    let (wav,clips)=encode_pcm16_wav_interleaved(&pressure,OUTPUT_RATE,channels as u16,full_scale_pa)?;
    let peak=pressure.iter().map(|p|p.abs()).fold(0.0_f64,f64::max);
    eprintln!("pressure WAV: frames={frames}, channels={channels}, rate_hz={OUTPUT_RATE}, full_scale_pa={full_scale_pa}, clips={clips}, peak_pa={peak}, decimator_delay_frames={}; no channel normalization; acceleration is step-average tagged at step end",decimator.delay_output_frames());
    for (channel,bake) in baked.iter().enumerate() {
        eprintln!("receiver {channel} origin_range_m={}, explicit_delay_s={}; other propagation remains in the fitted transfer",bake.range_m,bake.propagation_delay_s);
    }
    Ok(wav)
}

#[cfg(test)]
#[path = "stereo_tests.rs"]
mod tests;
