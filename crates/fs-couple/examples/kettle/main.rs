//! Profile-defined, single-headed sealed chamber using the existing physics owners.
//! cargo run --release -p fs-couple --example kettle -- instrument.fsk 4096
//! Output is mechanics/internal pressure, NOT microphone audio. See README.md.
use std::io::Write;
use fs_exec::CancelGate;
type Error = Box<dyn std::error::Error>;
mod spec;
mod model;
#[cfg(test)]
mod tests;

fn main() -> Result<(), Error> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() != 2 { return Err("usage: kettle SPEC.fsk MECHANICAL_STEPS".into()); }
    let spec = spec::Spec::load(&args[0])?;
    let steps: u64 = args[1].parse()?;
    if steps > 10_000_000 { return Err("at most ten million requested mechanical steps".into()); }
    let gate = CancelGate::new_clock_free();
    let mut model = model::Model::new(&spec,steps,&gate)?;
    eprintln!("supplied single-head chamber: volume_m3={}, head_mass_kg={}, head_modes={}, air_modes={}, rate_hz={}; rigid wall, nonlinear stretching, independent felt histories; no calibration or real-time claim",
        model.air.volume_m3(), model.film.mass_kg, model.modes.len(), model.cavity.cavity_modes(),spec.rate_hz);
    eprintln!("uncoupled_head_hz={:?}; uncoupled_air_hz={:?}; these are NOT coupled instrument pitches",
        model.modes.iter().map(|m| m.lambda.sqrt()/std::f64::consts::TAU).collect::<Vec<_>>(),
        model.air.modes().iter().map(|m| m.omega_rad_s/std::f64::consts::TAU).collect::<Vec<_>>());
    let observation = model.air.values_at([0.0,0.0,
        (spec.profile[0].z_m+spec.profile.last().ok_or("missing rim")?.z_m)/2.0])?;
    let mut out = std::io::BufWriter::new(std::io::stdout().lock());
    writeln!(out,"time_s,mallet_down_m,mallet_velocity_m_s,head_down_m,uniform_internal_pa,point_internal_pa,energy_j,loss_j,felt_crush_j,balance_j")?;
    for _ in 0..steps {
        let frame = model.mechanics.step(&model.force,&gate)?;
        let mut pressure = [0.0;8];
        model.cavity.pressures_into(model.mechanics.state(),&mut pressure[..model.cavity.cavity_modes()])?;
        let local = model.cavity.pressure_at(model.mechanics.state(),&observation)?;
        writeln!(out,"{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}",
            frame.time_s, model.tip_weight*model.mechanics.state()[0], model.tip_weight*model.mechanics.state()[1],
            model.head_displacement(), pressure[0],local,frame.stored_energy_j,frame.dissipated_energy_j,
            frame.felt_crush_loss_j,frame.balance_residual_j)?;
    }
    out.flush()?;
    Ok(())
}
