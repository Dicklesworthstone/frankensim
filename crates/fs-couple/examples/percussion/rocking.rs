//! Optional bounded rigid heave/tilt on the existing elastic shell and clock.
use super::{Error, Mechanics};
use std::io::Write;

pub fn option(args:&mut Vec<String>)->Result<Option<f64>,Error> {
    let mut positions=args.iter().enumerate().filter(|(_,a)|a.as_str()=="--cymbal-rocking");
    let Some((at,_))=positions.next() else {return Ok(None);};
    if positions.next().is_some() {return Err("--cymbal-rocking may be supplied once".into());}
    let limit:f64=args.get(at+1).ok_or("--cymbal-rocking needs MAXIMUM_ANGLE_RAD")?.parse()?;
    if !limit.is_finite() || limit<=0.0 || limit>0.1 {
        return Err("cymbal rocking requires a maximum angle in (0,0.1] rad".into());
    }
    args.drain(at..at+2);Ok(Some(limit))
}

pub fn admit_command(selected:Option<f64>,command:&str)->Result<(),Error> {
    if selected.is_some() && !matches!(command,"splash"|"splash-wav"|"splash-mic"|
        "hihat"|"hihat-wav"|"hihat-mic") {
        return Err("--cymbal-rocking requires a single or paired cymbal command".into());
    }
    Ok(())
}

pub fn observation(system:&Mechanics,body:usize)->Option<[f64;3]> {
    match system {
        Mechanics::Reference(s)=>s.shell_rotation(body),
        Mechanics::Nonlinear(s)=>s.shell_rotation(body),
        Mechanics::Substepped(s)=>s.shell_rotation(body),
        Mechanics::Driven{inner,..}=>observation(inner,body),
        Mechanics::Prepared(_)=>None,
    }
}

pub fn header(out:&mut impl Write,label:&str)->Result<(),Error> {
    for axis in ["x","y","z"] {write!(out,",{label}_rigid_rotation_{axis}_rad")?;}
    Ok(())
}

pub fn row(system:&Mechanics,body:usize,lower_shell:bool,out:&mut impl Write)->Result<(),Error> {
    let mut rotation=observation(system,body).ok_or("missing admitted cymbal rigid rotation")?;
    // A lower hi-hat shell is placed by a proper pi rotation about world x.
    // Axial rotation vectors transform by that same matrix, not a reflection.
    if lower_shell {rotation[1]=-rotation[1];rotation[2]=-rotation[2];}
    for value in rotation {write!(out,",{value:.17e}")?;}
    Ok(())
}

#[cfg(test)]
#[path="rocking_tests.rs"]
mod tests;
