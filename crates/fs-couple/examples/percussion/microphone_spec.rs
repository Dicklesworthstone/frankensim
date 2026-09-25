//! Strict, bounded physical receiver selection for the existing pressure path.
//! Omitted patterns mean ideal omni. No inferred axes, gain or capsule model.
use super::{Error, Microphone, Receiver};
use fs_bem::near_field::{FirstOrder, Options};
use std::{collections::BTreeMap, io::Read};

const MAX_BYTES: usize = 8192;

pub struct Spec { receivers: Vec<Receiver> }
impl Spec {
    pub fn parse(text: &str) -> Result<Self, Error> {
        if text.len()>MAX_BYTES {return Err("microphone file exceeds 8 KiB".into());}
        let mut header=false;let mut controls=None;let mut points=Vec::new();
        let mut patterns=BTreeMap::new();
        for (line, raw) in text.lines().enumerate() {
            let row=raw.split('#').next().unwrap_or("").trim();if row.is_empty(){continue;}
            let bad=||->Error{format!("microphone line {}: invalid or duplicate record",line+1).into()};
            if !header {
                if row!="frankensim-microphones-v1" {return Err(bad());}
                header=true;continue;
            }
            let fields:Vec<_>=row.split(',').map(str::trim).collect();
            let real=|i:usize|->Result<f64,Error>{
                let value:f64=fields[i].parse().map_err(|_|bad())?;
                if !value.is_finite(){return Err(bad());}Ok(value)
            };
            match fields[0] {
                "near_field" if fields.len()==5 && controls.is_none()=>{
                    let minimum_clearance_m=real(1)?;
                    let options=Options {relative_tolerance:real(2)?,
                        maximum_depth:fields[3].parse().map_err(|_|bad())?,
                        maximum_kernel_evaluations:fields[4].parse().map_err(|_|bad())?};
                    controls=Some((minimum_clearance_m,options));
                }
                "receiver" if fields.len()==4 && points.len()<2=>points.push([real(1)?,real(2)?,real(3)?]),
                "pattern" if fields.len()==6=>{
                    let index:usize=fields[1].parse().map_err(|_|bad())?;
                    if index>=2 || patterns.contains_key(&index){return Err(bad());}
                    let pattern=FirstOrder::new(real(2)?,[real(3)?,real(4)?,real(5)?])?;
                    patterns.insert(index,pattern);
                }
                _=>return Err(bad()),
            }
        }
        if !header || points.is_empty() || patterns.keys().any(|&index|index>=points.len()) {
            return Err("microphone file needs one or two receivers and valid pattern indices".into());
        }
        let (clearance,options)=controls.ok_or("microphone file needs explicit near_field controls")?;
        let receivers=points.into_iter().enumerate().map(|(i,position)| {
            Microphone::new(position,clearance,options,patterns.get(&i).copied().unwrap_or_default())
                .map(Receiver::NearField)
        }).collect::<Result<Vec<_>,_>>()?;
        Ok(Self{receivers})
    }
    pub fn load(path:&str)->Result<Self,Error> {
        let mut text=String::new();
        std::fs::File::open(path)?.take((MAX_BYTES+1) as u64).read_to_string(&mut text)?;
        Self::parse(&text)
    }
    pub fn admit_command(&self,command:&str,positional_receiver:bool,right_receiver:bool)->Result<(),Error> {
        if !matches!(command,"splash-mic"|"drum-mic"|"drum-stretch-mic"|"drum-modal-mic"|
            "snare-mic"|"snare-off-mic"|"hihat-mic") {
            return Err("--microphone-spec requires an existing finite-point -mic command".into());
        }
        if positional_receiver || right_receiver {
            return Err("--microphone-spec supplies every receiver; do not combine it with positional XYZ or --microphone-right".into());
        }
        Ok(())
    }
    pub fn into_receivers(self)->Vec<Receiver>{self.receivers}
}

/// Failed admission leaves arguments unchanged. Only successful cold input is
/// removed; geometry and source construction still have not advanced mechanics.
pub fn option(args:&mut Vec<String>)->Result<Option<Spec>,Error> {
    let indices:Vec<_>=args.iter().enumerate().filter_map(|(i,s)|(s=="--microphone-spec").then_some(i)).collect();
    if indices.len()>1 {return Err("--microphone-spec may be supplied only once".into());}
    let Some(&i)=indices.first() else{return Ok(None);};
    let path=args.get(i+1).filter(|s|!s.is_empty()&&!s.starts_with("--"))
        .ok_or("--microphone-spec needs an input path; see MICROPHONES.md")?;
    let spec=Spec::load(path)?;
    args.drain(i..i+2);Ok(Some(spec))
}

#[cfg(test)]
#[path="receivers/input_tests.rs"]
mod tests;
