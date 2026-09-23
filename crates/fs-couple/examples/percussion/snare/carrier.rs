//! Authored force on the same translating support that carries the wire mass.
//! Input admission and force interpolation reuse the existing mechanical owners.
use super::{Error, SnareSet, TranslatingSupport};
use super::super::{Mechanics, mechanics::drive::{Input, Program}};
use fs_couple::render::plate::impact::supported::SupportObservation;
use std::io::Read;

const MAX_BYTES: u64=64*1024;

pub struct Specification {
    pub support: TranslatingSupport,
    pub program: Program,
}
impl Specification {
    pub fn parse(text:&str)->Result<Self,Error> {
        if text.len() as u64>MAX_BYTES {return Err("snare carrier file exceeds 64 KiB".into());}
        let mut rows=text.lines().enumerate().filter_map(|(i,line)| {
            let row=line.split('#').next().unwrap_or("").trim();
            (!row.is_empty()).then_some((i+1,row))
        });
        if rows.next().map(|(_,row)|row)!=Some("frankensim-snare-carrier-v1") {
            return Err("snare carrier requires frankensim-snare-carrier-v1 header".into());
        }
        let mut support=None;let mut initial=None;let mut forces=String::new();
        for (line,row) in rows {
            let fields:Vec<_>=row.split(',').map(str::trim).collect();
            let bad=||format!("snare carrier line {line}: unknown, duplicate, malformed or nonfinite record");
            let number=|i:usize|->Result<f64,Error> {
                let value=fields[i].parse::<f64>().map_err(|_|bad())?;
                if !value.is_finite() {return Err(bad().into());}Ok(value)
            };
            match (fields[0],fields.len()) {
                ("support",6) if support.is_none()=>support=Some([number(1)?,number(2)?,number(3)?,number(4)?,number(5)?]),
                ("initial",3) if initial.is_none()=>initial=Some([number(1)?,number(2)?]),
                ("force",3)=>{
                    forces.push_str(fields[1]);forces.push(',');forces.push_str(fields[2]);forces.push('\n');
                }
                _=>return Err(bad().into()),
            }
        }
        let [mass_kg,stiffness_n_m,damping_n_s_m,maximum_travel_m,maximum_slope]=support.ok_or("missing carrier support record")?;
        let [initial_position_m,initial_velocity_m_s]=initial.ok_or("missing carrier initial motion record")?;
        let support=TranslatingSupport {mass_kg,stiffness_n_m,damping_n_s_m,maximum_travel_m,
            maximum_slope,initial_position_m,initial_velocity_m_s};
        support.validate()?;
        // No parallel force law: preserve the existing finite, increasing-time,
        // zero-endpoint program and exact per-tick impulse integration.
        Ok(Self {support,program:Program::parse(&forces)?})
    }
    pub fn load(path:&str)->Result<Self,Error> {
        let mut text=String::new();std::fs::File::open(path)?.take(MAX_BYTES+1).read_to_string(&mut text)?;
        Self::parse(&text)
    }
    pub fn into_input(self,system:&Mechanics,body:usize)->Result<Input,Error> {
        let (coordinate,tip_weight)=force_port(system,body).ok_or("missing physical snare carrier force port")?;
        Ok(Input {program:self.program,coordinate,tip_weight})
    }
}

pub fn option(args:&mut Vec<String>)->Result<Option<String>,Error> {
    let mut path=None;let mut i=0;
    while i<args.len() {
        if args[i]!="--snare-carrier" {i+=1;continue;}
        if path.is_some() {return Err("--snare-carrier may be supplied only once".into());}
        let value=args.get(i+1).ok_or("--snare-carrier requires an input file")?;
        if value.starts_with("--") {return Err("--snare-carrier needs a file, not another option".into());}
        path=Some(value.clone());args.drain(i..i+2);
    }
    Ok(path)
}

/// A carrier is independent of the supplied wire card and head law. Non-snare
/// requests refuse before opening the file rather than acquiring an unused input.
pub fn select(path:Option<&str>,mut snare:Option<SnareSet>)->Result<(Option<SnareSet>,Option<Specification>),Error> {
    let Some(path)=path else {return Ok((snare,None));};
    let bank=snare.as_mut().ok_or("--snare-carrier applies only to snare[-off][-wav|-mic]")?;
    if bank.carrier.is_some() {return Err("snare carrier has already been supplied".into());}
    let spec=Specification::load(path)?;bank.carrier=Some(spec.support);bank.mode_count()?;
    Ok((snare,Some(spec)))
}

pub fn force_port(system:&Mechanics,body:usize)->Option<(usize,f64)> {
    match system {
        Mechanics::Reference(s)=>s.support_force_port(body),Mechanics::Nonlinear(s)=>s.support_force_port(body),
        Mechanics::Substepped(s)=>s.support_force_port(body),Mechanics::Prepared(_)=>None,
        Mechanics::Driven {inner,..}=>force_port(inner,body),
    }
}
pub fn observe(system:&Mechanics,body:usize)->Option<SupportObservation> {
    match system {
        Mechanics::Reference(s)=>s.support_observation(body),Mechanics::Nonlinear(s)=>s.support_observation(body),
        Mechanics::Substepped(s)=>s.support_observation(body),Mechanics::Prepared(_)=>None,
        Mechanics::Driven {inner,..}=>observe(inner,body),
    }
}

#[cfg(test)]
mod tests;
