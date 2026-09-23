//! Select physical nonlinearity explicitly; never erase supplied head or wire laws.
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
pub fn admit_prepared_command(prepared:bool,stretching:bool,command:&str)->Result<(),Error> {
    admit_command(stretching,command)?;
    if prepared && !(matches!(command,"splash"|"splash-wav"|"splash-mic"|
        "drum"|"drum-wav"|"drum-mic"|"drum-stretch"|"drum-stretch-wav"|"drum-stretch-mic")
        || stretching && snare(command)) {
        return Err("nonlinear preparation requires splash, drum, drum-stretch, or snare with head/wire stretching or a moving carrier; no silent conversion of linear modal mechanics".into());
    }
    Ok(())
}
pub fn admit_image(linear_prepared:bool,wires:bool,stretching:bool)->Result<(),Error> {
    if stretching && linear_prepared {
        return Err("stretching heads/wires and moving supports require coupled nonlinear-capable mechanics, not the linear modal image".into());
    }
    if wires && !linear_prepared && !stretching {
        return Err("a fully linear snare keeps its original prepared modal image; select head/wire stretching or a carrier explicitly".into());
    }
    Ok(())
}

#[cfg(test)]
mod tests;
