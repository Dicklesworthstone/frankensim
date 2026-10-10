//! Explicit proportional flexural memory on the actual curved-shell operator.
//! The existing material and impact owners supply all storage and time evolution.
use super::{Error, ImpactSystem};
use fs_couple::render::plate::impact::relaxation::{InitialMemory,ShellBendingSpectrum};
use fs_material::visco::GeneralizedMaxwell;
use std::io::Read;

const MAX_BYTES:usize=65536;
#[derive(Clone,Debug)]
pub struct Spec {
    shells:[Option<GeneralizedMaxwell>;2],
    initial:InitialMemory,
    band_hz:[f64;2],
    paired:bool,
}
impl Spec {
    pub fn read(text:&str,paired:bool)->Result<Self,Error> {
        if text.len()>MAX_BYTES {return Err("shell relaxation exceeds 64 KiB".into());}
        let (mut header,mut replace,mut initial,mut band)=(false,false,None,None);
        let mut shells:[Option<GeneralizedMaxwell>;2]=[None,None];
        for (line,raw) in text.lines().enumerate() {
            let row=raw.split('#').next().unwrap_or("").trim();if row.is_empty(){continue;}
            let bad=||->Error{format!("shell relaxation line {}: invalid, duplicate or missing record",line+1).into()};
            if !header {
                if row!="frankensim-shell-relaxation-v1" {return Err(bad());}
                header=true;continue;
            }
            let f:Vec<_>=row.split(',').map(str::trim).collect();
            let number=|i:usize|->Result<f64,Error>{let x:f64=f[i].parse().map_err(|_|bad())?;
                if !x.is_finite(){return Err(bad());}Ok(x)};
            let side=|name:&str|->Result<usize,Error>{match (paired,name) {
                (false,"single")|(true,"upper")=>Ok(0),(true,"lower")=>Ok(1),_=>Err(bad()),
            }};
            match f[0] {
                "intrinsic_loss" if f.len()==2 && !replace && f[1]=="replace"=>replace=true,
                "initial" if f.len()==2 && initial.is_none()=>initial=Some(match f[1] {
                    "relaxed"=>InitialMemory::Relaxed,"unrelaxed"=>InitialMemory::Unrelaxed,_=>return Err(bad()),
                }),
                "band_hz" if f.len()==3 && band.is_none()=>{
                    let b=[number(1)?,number(2)?];if b[0]<0. || b[1]<=b[0] {return Err(bad());}band=Some(b);
                }
                "shell" if f.len()==2=>{
                    let i=side(f[1])?;if shells[i].is_some(){return Err(bad());}
                    // These are ratios to actual equilibrium bending stiffness,
                    // including heterogeneous sections, not a modulus in Pa.
                    shells[i]=Some(GeneralizedMaxwell::new(1.,vec![])?);
                }
                "branch" if f.len()==4=>{
                    let i=side(f[1])?;let law=shells[i].as_mut().ok_or_else(bad)?;
                    if law.terms.len()==8 {return Err(bad());}
                    let mut terms=law.terms.clone();terms.push((number(2)?,number(3)?));
                    *law=GeneralizedMaxwell::new(1.,terms)?;
                }
                _=>return Err(bad()),
            }
        }
        if !replace || shells.iter().all(Option::is_none) {
            return Err("shell relaxation needs explicit intrinsic_loss,replace and at least one selected shell".into());
        }
        for law in shells.iter().flatten() {
            if !(1.+law.terms.iter().map(|t|t.0).sum::<f64>()).is_finite() {
                return Err("instantaneous flexural stiffness ratio overflows".into());
            }
        }
        Ok(Self{shells,paired,initial:initial.ok_or("missing initial shell memory")?,
            band_hz:band.ok_or("missing shell material band")?})
    }
    pub fn load(path:&str,paired:bool)->Result<Self,Error> {
        let mut text=String::new();std::fs::File::open(path)?.take((MAX_BYTES+1) as u64).read_to_string(&mut text)?;
        Self::read(&text,paired)
    }
    pub fn selected(&self,side:usize)->bool {self.shells.get(side).is_some_and(Option::is_some)}
    pub fn admit_windows(&self,windows:&[[f64;2]])->Result<(),Error> {
        if windows.len()!=if self.paired {2}else{1} {
            return Err("shell material selection does not match the single/paired instrument".into());
        }
        for (i,window) in windows.iter().enumerate() {
            if self.selected(i) && (window.iter().any(|v|!v.is_finite())
                || window[0]<self.band_hz[0] || window[1]>self.band_hz[1]) {
                return Err("requested shell modes exceed the supplied flexural material band".into());
            }
        }
        Ok(())
    }
    pub fn attach(&self,system:ImpactSystem)->Result<ImpactSystem,Error> {
        let spectra:Vec<_>=self.shells.iter().enumerate().filter_map(|(i,law)|law.as_ref().map(|law|
            ShellBendingSpectrum{body:i+1,branches:law.terms.clone(),band_hz:self.band_hz})).collect();
        Ok(system.with_shell_bending_relaxation(&spectra,self.initial,256)?)
    }
}

pub fn option(args:&mut Vec<String>)->Result<Option<String>,Error> {
    let positions:Vec<_>=args.iter().enumerate().filter_map(|(i,s)|(s=="--shell-relaxation").then_some(i)).collect();
    if positions.len()>1 {return Err("--shell-relaxation may be supplied once".into());}
    let Some(&i)=positions.first() else{return Ok(None);};
    let path=args.get(i+1).filter(|s|!s.is_empty()&&!s.starts_with("--"))
        .ok_or("--shell-relaxation needs a material file")?.clone();
    args.drain(i..i+2);Ok(Some(path))
}
pub fn admit_command(selected:bool,command:&str)->Result<(),Error> {
    if selected && !matches!(command,"splash"|"splash-wav"|"splash-mic"|"hihat"|"hihat-wav"|"hihat-mic") {
        return Err("shell relaxation requires single or paired cymbals; use --head-relaxation for drumheads".into());
    }
    Ok(())
}

#[cfg(test)]
#[path="shell_relaxation_tests.rs"]
mod tests;
