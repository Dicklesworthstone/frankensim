//! Profile-defined, single-headed sealed chamber using the existing physics owners.
//! cargo run --release -p fs-couple --example kettle -- instrument.fsk 4096
//! Default output is mechanical CSV; --wav selects coupled exterior pressure.
use std::io::Write;
use fs_exec::CancelGate;
type Error = Box<dyn std::error::Error>;
mod spec;
mod model;
mod exterior;
mod boundary;
mod audio;
mod tuning;
mod execution;
#[cfg(test)]
mod tests;

fn main() -> Result<(), Error> {
    let mut args:Vec<_>=std::env::args().skip(1).collect();
    if args.len()<2 {return Err("usage: kettle SPEC.fsk COUNT [--wav ...] [--edge-tensioner M K C STRAIN] [--edge-force SAMPLE N]".into());}
    let spec=spec::Spec::load(&args[0])?;
    let reciprocal=audio::reciprocal::option(&mut args)?;
    let execution=execution::option(&mut args)?;
    let tensioner=tuning::option(&mut args,spec.maximum_slope)?;
    if args.len()!=2 && args.len()!=11 {return Err("usage: kettle SPEC.fsk COUNT [--wav FULL_SCALE_PA MIN_HZ MAX_HZ INTERVALS ORDER X_M Y_M Z_M] plus explicit edge options".into());}
    if reciprocal.is_some() && args.len()!=11 {return Err("--reciprocal-air requires --wav; it does not alter CSV mechanics".into());}
    let steps: u64 = args[1].parse()?;
    if steps > 10_000_000 { return Err("at most ten million requested mechanical steps".into()); }
    let gate = CancelGate::new_clock_free();
    eprintln!("nonlinear execution: max_substep_depth={}, max_solve_attempts_per_tick={}, Newton_backtracking_halvings={}; unchanged output clock and physical limits, no real-time qualification",
        execution.substeps.max_depth,execution.substeps.max_attempts,execution.backtracking_halvings);
    if args.len()==11 {
        let settings=audio::Settings::read(&args[2..],&spec)?;
        if steps==0 || steps>480_000 { return Err("audio count must be 1..480000 output frames".into()); }
        let ticks=steps.checked_mul(u64::from(spec.rate_hz/audio::OUTPUT_RATE)).ok_or("audio clock overflow")?;
        let mut baked=None;
        let mut model=model::Model::with_execution(&spec,ticks,&gate,tensioner,execution,|film,modes| {
            let (load,observer)=match reciprocal {
                Some(controls)=>audio::reciprocal::bake(&spec,film,modes,settings,controls,&gate)?,
                None=>audio::bake(&spec,film,modes,settings,&gate)?,
            };
            baked=Some(observer); Ok(Some(load))
        })?;
        let wav=audio::render(&mut model,&baked.ok_or("missing acoustic bake")?,&spec,settings,steps as usize,&gate)?;
        std::io::stdout().lock().write_all(&wav)?;
        return Ok(());
    }
    let mut model=model::Model::with_execution(&spec,steps,&gate,tensioner,execution,|_,_|Ok(None))?;
    eprintln!("supplied single-head chamber: volume_m3={}, head_mass_kg={}, head_modes={}, air_modes={}, rate_hz={}; rigid wall, nonlinear stretching, independent felt histories; no calibration or real-time claim",
        model.air.volume_m3(), model.film.mass_kg, model.modes.len(), model.cavity.cavity_modes(),spec.rate_hz);
    eprintln!("uncoupled_head_hz={:?}; uncoupled_air_hz={:?}; these are NOT coupled instrument pitches",
        model.modes.iter().map(|m| m.lambda.sqrt()/std::f64::consts::TAU).collect::<Vec<_>>(),
        model.air.modes().iter().map(|m| m.omega_rad_s/std::f64::consts::TAU).collect::<Vec<_>>());
    let observation = model.air.values_at([0.0,0.0,
        (spec.profile[0].z_m+spec.profile.last().ok_or("missing rim")?.z_m)/2.0])?;
    let mut out = std::io::BufWriter::new(std::io::stdout().lock());
    write!(out,"time_s,mallet_down_m,mallet_velocity_m_s,head_down_m,uniform_internal_pa,point_internal_pa,energy_j,loss_j,felt_crush_j,balance_j")?;
    if model.mechanics.dilation_force_port(1).is_some() {write!(out,",edge_outward_m,edge_baseline_tension_n_m,supplied_work_j")?;}
    writeln!(out)?;
    // One slot per retained cavity mode: the meridian basis keeps up to 256 (model.rs), not the 8 it once did.
    let mut pressure = vec![0.0;model.cavity.cavity_modes()];
    for _ in 0..steps {
        let frame=model.step(&gate)?;
        model.cavity.pressures_into(model.mechanics.state(),&mut pressure)?;
        let local = model.cavity.pressure_at(model.mechanics.state(),&observation)?;
        write!(out,"{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}",
            frame.time_s, model.tip_weight*model.mechanics.state()[0], model.tip_weight*model.mechanics.state()[1],
            model.head_displacement(), pressure[0],local,frame.stored_energy_j,frame.dissipated_energy_j,
            frame.felt_crush_loss_j,frame.balance_residual_j)?;
        if let Some(edge)=model.mechanics.dilation_observation(1) {
            write!(out,",{:.12e},{:.12e},{:.12e}",edge.rim_displacement_m,edge.baseline_tension_n_m,frame.supplied_work_j)?;
        }
        writeln!(out)?;
    }
    out.flush()?;
    Ok(())
}
