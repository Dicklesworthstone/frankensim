//! Select the numerical owner without changing supplied head or wire laws.
use super::Error;

pub fn option(args:&mut Vec<String>)->Result<bool,Error> {
    let count=args.iter().filter(|s|s.as_str()=="--head-stretching").count();
    if count>1 {return Err("--head-stretching may be supplied only once".into());}
    args.retain(|s|s!="--head-stretching");Ok(count==1)
}
fn snare(command:&str)->bool {
    matches!(command,"snare"|"snare-wav"|"snare-mic"|"snare-off"|"snare-off-wav"|"snare-off-mic")
}
pub fn admit_command(stretching:bool,command:&str)->Result<(),Error> {
    if stretching && !snare(command) {
        return Err("--head-stretching applies to snare[-off][-wav|-mic]; use drum-stretch for a drum without wires".into());
    }
    Ok(())
}
pub fn admit_prepared_command(prepared:bool,nonlinear:bool,command:&str)->Result<(),Error> {
    admit_command(nonlinear,command)?;
    if prepared && !(matches!(command,"splash"|"splash-wav"|"splash-mic"|
        "drum"|"drum-wav"|"drum-mic"|"drum-stretch"|"drum-stretch-wav"|"drum-stretch-mic")
        || snare(command)) {
        return Err("nonlinear preparation requires splash, drum, drum-stretch, or snare; drum-modal keeps its explicitly selected modal owner".into());
    }
    Ok(())
}
pub fn admit_image(linear_prepared:bool,nonlinear:bool)->Result<(),Error> {
    if nonlinear && linear_prepared {
        return Err("nonlinear head/wire material, felt contact and moving supports require coupled nonlinear-capable mechanics".into());
    }
    // The joint contact owner also admits linear body potentials. In particular,
    // analytic/substepped integration and radiation memory must not force a
    // caller to invent head stretching, material memory or flexible hardware.
    Ok(())
}

#[cfg(test)]
mod tests;
